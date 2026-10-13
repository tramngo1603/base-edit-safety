"""Select complete experiments and account for their declared resource instances.

Question satisfaction comes from the source-bound question evaluator. Selection
never changes the evidence records. Exact and bounded searches report different
completion guarantees, including when their declared limits interrupt a search.
"""
from __future__ import annotations

from collections import Counter, OrderedDict
from copy import deepcopy
from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
from itertools import combinations
import json
from pathlib import Path
import random
import re
import resource
import sys
import time

from .common import PopulationError
from .comparison_contracts import make_record, resource_compatible
from .comparison_identity import canonical_bytes, content_id


SELECTORS = ("complete_bundle_greedy", "donor_first", "adapted_coverage",
             "canonical_bundle", "random_bundle", "exact_tiny")
PLANNER_VERSION = "complete-experiments-v2"
HARD_EXACT_ACTION_CAP = 22
HARD_EXACT_RESOURCE_CAP = 128
HARD_EXACT_STATE_CAP = 2_000_000
_RUNTIME_FIELDS = {"runtime_seconds", "process_peak_rss_bytes", "memory_measurement"}
_AUDIT_HASH_EXCLUDED_FIELDS = _RUNTIME_FIELDS | {"plan_id", "plan_record"}


@dataclass(frozen=True)
class PlannerLimits:
    """All search bounds are positive and appear in the returned plan."""

    beam_width: int
    max_expansions: int
    max_alternatives_per_question: int
    exact_action_cap: int
    exact_resource_cap: int = 128
    exact_state_cap: int = 2_000_000

    def __post_init__(self):
        for name, value in asdict(self).items():
            if type(value) is not int or value < 1:
                raise PopulationError("planning." + name + " must be a positive integer")
        for name, ceiling in (("exact_action_cap", HARD_EXACT_ACTION_CAP),
                              ("exact_resource_cap", HARD_EXACT_RESOURCE_CAP),
                              ("exact_state_cap", HARD_EXACT_STATE_CAP)):
            if getattr(self, name) > ceiling:
                raise PopulationError(f"planning.{name} exceeds the implementation limit {ceiling}")

    @classmethod
    def from_policy(cls, policy):
        try:
            return cls(**{name: policy[name] for name in (
                "beam_width", "max_expansions", "max_alternatives_per_question", "exact_action_cap")},
                **{name: policy[name] for name in ("exact_resource_cap", "exact_state_cap") if name in policy})
        except KeyError as exc:
            raise PopulationError("Missing explicit planning limit: " + str(exc.args[0])) from exc


@dataclass(frozen=True)
class _Selection:
    actions: tuple = ()
    assignments: tuple = ()

    def assignment_rows(self):
        return [{"action_id": aid, "assignment_key": key, "resource_id": rid}
                for aid, key, rid in self.assignments]


def _questions_api():
    from . import comparison_questions
    return comparison_questions


def _reprepare(question_set):
    # A caller-supplied eligibility flag, index or cached objective is not proof.
    return _questions_api().prepare_question_set(deepcopy(question_set.study),
                                                question_set.records)


def _config(question_set):
    return question_set.study.get("config", question_set.study)


def _selection(action_ids, assignments=()):
    actions = tuple(sorted(action_ids))
    if len(set(actions)) != len(actions):
        raise PopulationError("A selected action may occur only once")
    rows = []
    for row in assignments:
        if set(row) != {"action_id", "assignment_key", "resource_id"} or not all(
                isinstance(row[key], str) and row[key] for key in row):
            raise PopulationError("Resource assignments require action_id, assignment_key and resource_id")
        rows.append((row["action_id"], row["assignment_key"], row["resource_id"]))
    if len({(aid, key) for aid, key, _ in rows}) != len(rows):
        raise PopulationError("Each action/assignment key must choose exactly one resource instance")
    return _Selection(actions, tuple(sorted(rows)))


def _eligible(question_set, action_id):
    return bool(question_set.eligibility.get(action_id, {}).get("eligible", False))


def _account(question_set, selected, budget, *, partial_assignments=False):
    """Recompute costs and summed capacity; fixed sharing names never merge IDs."""
    actions, resources = question_set.actions, question_set.resources
    assigned = {(aid, key): rid for aid, key, rid in selected.assignments}
    needed, used, demands, problems = set(), set(), Counter(), []
    action_resources = {}
    unknown = set(selected.actions) - actions.keys()
    if unknown:
        problems.append({"code": "unknown_action", "action_ids": sorted(unknown)})
    for aid in selected.actions:
        if aid not in actions:
            continue
        action = actions[aid]
        if not _eligible(question_set, aid):
            problems.append({"code": "ineligible_action", "action_id": aid,
                             "reason_codes": list(question_set.eligibility.get(aid, {}).get("reason_codes", []))})
        missing = set(action["dependency_action_ids"]) - set(selected.actions)
        if missing:
            problems.append({"code": "missing_prerequisites", "action_id": aid,
                             "action_ids": sorted(missing)})
        for demand in action["resource_demands"]:
            if "resource_id" in demand:
                rid = demand["resource_id"]
            else:
                key = aid, demand["assignment_key"]
                needed.add(key)
                rid = assigned.get(key)
                if rid is None:
                    if not partial_assignments:
                        problems.append({"code": "missing_resource_assignment", "action_id": aid,
                                         "assignment_key": key[1]})
                    continue
                if rid not in demand["eligible_resource_ids"]:
                    problems.append({"code": "resource_not_in_declared_alternatives", "action_id": aid,
                                     "assignment_key": key[1], "resource_id": rid})
                    continue
            if rid not in resources:
                problems.append({"code": "unknown_resource", "resource_id": rid})
                continue
            reasons = resource_compatible(action, resources[rid], question_set.index)
            if reasons:
                problems.append({"code": "incompatible_resource", "action_id": aid,
                                 "resource_id": rid, "reason_codes": reasons})
            used.add(rid)
            action_resources.setdefault(aid, set()).add(rid)
            demands[rid, demand["dimension"]] += demand["quantity"]
    for aid, key in sorted(assigned.keys() - needed):
        problems.append({"code": "extraneous_resource_assignment", "action_id": aid, "assignment_key": key})
    if not partial_assignments:
        for aid in set(selected.actions) & actions.keys():
            action = actions[aid]
            if action["operation_type"] not in {"measurement_panel", "locus_assay", "intended_assay", "functional_assay", "mock_control"}:
                continue
            closure, pending = set(), [aid]
            while pending:
                parent = pending.pop()
                if parent in closure or parent not in actions:
                    continue
                closure.add(parent)
                pending.extend(actions[parent]["dependency_action_ids"])
            actual_types = {resources[rid]["resource_type"] for parent in closure
                            for rid in action_resources.get(parent, ())}
            required_types = {"study_setup", "donor_activation", "specimen_preparation", "sequencing_batch"}
            required_types.add("mock_control" if action["operation_type"] == "mock_control" else "assay_library")
            if action["design_id"]["value"] is not None and action["operation_type"] != "mock_control":
                required_types.add("treatment")
            if action["replicate_position"] > 1:
                required_types.add("additional_replicate")
            if required_types - actual_types:
                problems.append({"code": "assigned_cost_packages_missing", "action_id": aid,
                                 "resource_types": sorted(required_types - actual_types)})
    active = used | {rid for rid, resource in resources.items()
                    if resource["activation_rule"] == "always" or
                    (resource["activation_rule"] == "nonempty_plan" and selected.actions)}
    units = {(row["cost_unit"], row["cost_scale"]) for row in resources.values()}
    if len(units) > 1:
        problems.append({"code": "mixed_cost_units_or_scales"})
    capacity_rows = []
    for (rid, dimension), quantity in sorted(demands.items()):
        cap = next((row["capacity"] for row in resources[rid]["capacities"]
                    if row["dimension"] == dimension), None)
        if cap is None:
            problems.append({"code": "undeclared_capacity_dimension", "resource_id": rid,
                             "dimension": dimension})
        elif cap["kind"] == "finite" and quantity > cap["quantity"]:
            problems.append({"code": "resource_capacity_exceeded", "resource_id": rid,
                             "dimension": dimension, "demand": quantity, "capacity": cap["quantity"]})
        capacity_rows.append({"resource_id": rid, "dimension": dimension,
                              "quantity": quantity, "capacity": deepcopy(cap)})
    cost = sum(resources[rid]["cost"] for rid in active)
    if cost > budget:
        problems.append({"code": "budget_exceeded", "total_cost": cost, "budget": budget})
    return {"feasible": not problems, "total_cost": cost, "budget": budget,
            "cost_unit": next(iter(units))[0] if units else _config(question_set)["policies"]["planning"]["cost_unit"],
            "cost_scale": next(iter(units))[1] if units else 1,
            "active_resource_ids": sorted(active), "resource_demands": capacity_rows,
            "cost_breakdown": [{"resource_id": rid, "resource_type": resources[rid]["resource_type"],
                                "sharing_key": resources[rid]["sharing_key"], "cost": resources[rid]["cost"],
                                "activation_rule": resources[rid]["activation_rule"]} for rid in sorted(active)],
            "problems": problems}


def evaluate_plan(question_set, selected_action_ids=(), resource_assignments=(), *, budget=None):
    """Validate a proposed selection afresh, ignoring any supplied solver result."""
    return _evaluate_prepared(_reprepare(question_set), selected_action_ids,
                              resource_assignments, budget=budget)


def _evaluate_prepared(qset, selected_action_ids=(), resource_assignments=(), *, budget=None):
    """Recompute a selection against inputs owned by the current operation."""
    if budget is None:
        budget = _config(qset)["policies"]["planning"]["total_budget"]
    _budget(budget)
    selected = _selection(selected_action_ids, resource_assignments)
    accounting = _account(qset, selected, budget)
    try:
        requirements = _questions_api().evaluate_requirements(qset, selected.actions)
    except PopulationError as exc:
        if accounting["feasible"]:
            raise
        requirements = _questions_api().evaluate_requirements(qset, ())
        accounting["problems"].append({"code": "selected_requirements_invalid", "reason": str(exc)})
    required = sorted(qid for qid, question in qset.questions.items() if question["hard_required"])
    unmet = sorted(set(required) - set(requirements["completed_question_ids"]))
    return {**accounting, "requirements": requirements, "unmet_hard_required_question_ids": unmet,
            "meets_request": accounting["feasible"] and not unmet}


def _budget(value):
    if type(value) is not int or value < 0:
        raise PopulationError("planning.total_budget must be a nonnegative integer")


class _Search:
    def __init__(self, question_set, budget, limits):
        self.qset, self.budget, self.limits = question_set, budget, limits
        self.evaluations = OrderedDict()
        self.accounts = OrderedDict()
        self.evaluation_cache_size = max(1, min(256, 100_000 // max(1, len(question_set.slots))))
        self.bundle_catalog = {}
        self.closures = {}
        self._coverage = None

    def coverage(self, selected):
        if self._coverage is None:
            self._coverage = _Coverage(self.qset)
        return self._coverage.units(selected.actions)

    def evaluate(self, selected):
        if selected.actions not in self.evaluations:
            self.evaluations[selected.actions] = _questions_api().evaluate_requirements(self.qset, selected.actions)
            if len(self.evaluations) > self.evaluation_cache_size:
                self.evaluations.popitem(last=False)
        self.evaluations.move_to_end(selected.actions)
        return self.evaluations[selected.actions]

    def account(self, selected):
        if selected not in self.accounts:
            self.accounts[selected] = _account(self.qset, selected, self.budget)
            if len(self.accounts) > 256:
                self.accounts.popitem(last=False)
        self.accounts.move_to_end(selected)
        return self.accounts[selected]

    def closure(self, action_ids):
        closed = set()
        for aid in action_ids:
            if aid not in self.closures:
                self.closures[aid] = tuple(_questions_api().dependency_closure(self.qset, (aid,)))
            closed.update(self.closures[aid])
        return tuple(sorted(closed))

    def donor_bundles(self, question_id):
        if question_id not in self.bundle_catalog:
            self.bundle_catalog[question_id] = _questions_api().eligible_donor_bundles(self.qset, question_id)
        return self.bundle_catalog[question_id]

    def assignments(self, action_ids, base, stats, limit):
        """Enumerate real alternative instances, pruning exact partial demand."""
        groups = {}
        fixed = {(aid, key): rid for aid, key, rid in base.assignments}
        for aid in action_ids:
            action = self.qset.actions[aid]
            for demand in action["resource_demands"]:
                if "resource_id" in demand:
                    continue
                key = aid, demand["assignment_key"]
                compatible = {rid for rid in demand["eligible_resource_ids"]
                              if not resource_compatible(action, self.qset.resources[rid], self.qset.index)}
                groups[key] = groups.get(key, compatible) & compatible
        if any(key not in groups or rid not in groups[key] for key, rid in fixed.items()):
            return
        pending = sorted(groups.keys() - fixed.keys())
        initial = _Selection(tuple(sorted(action_ids)), tuple(sorted((a, k, r) for (a, k), r in fixed.items())))
        if not _account(self.qset, initial, self.budget, partial_assignments=True)["feasible"]:
            stats["constraint_rejections"] += 1
            return
        def visit(position, values):
            if stats["expansions"] >= limit:
                stats["expansion_limit_reached"] = True
                return
            if position == len(pending):
                stats["expansions"] += 1
                selected = _Selection(tuple(sorted(action_ids)), tuple(sorted((a, k, r) for (a, k), r in values.items())))
                if self.account(selected)["feasible"]:
                    yield selected
                else:
                    stats["constraint_rejections"] += 1
                return
            key = pending[position]
            for rid in sorted(groups[key]):
                if stats["expansions"] >= limit:
                    stats["expansion_limit_reached"] = True
                    return
                stats["expansions"] += 1
                next_values = {**values, key: rid}
                partial = _Selection(tuple(sorted(action_ids)), tuple(sorted((a, k, r) for (a, k), r in next_values.items())))
                if _account(self.qset, partial, self.budget, partial_assignments=True)["feasible"]:
                    yield from visit(position + 1, next_values)
                else:
                    stats["constraint_rejections"] += 1
        yield from visit(0, fixed)


def _stats(method, limits):
    return {"method": method, "limits": asdict(limits), "expansions": 0,
            "constraint_rejections": 0, "biologically_ineligible_actions": 0,
            "beam_pruned": 0, "alternative_pruned": 0, "deduplicated": 0,
            "expansion_limit_reached": False, "exhaustive": False,
            "search_space_exhausted": False, "status": "searching"}


def _canonical(selected):
    return selected.actions, selected.assignments


def _construction_limited(search, question_ids):
    return sorted(qid for qid in question_ids if search.qset.questions[qid].get("construction_status") == "excluded_by_declared_limit")


def _completed_support(search, selected, base):
    """Keep complete experimental requirements and their prerequisites, not fill."""
    evaluation = search.evaluate(selected)
    needed = set(base.actions)
    for qid in evaluation["completed_question_ids"]:
        donors = set(evaluation["questions"][qid]["selected_donor_ids"])
        for slot_id in search.qset.questions[qid]["requirement_slot_ids"]:
            if search.qset.slots[slot_id]["donor_id"] not in donors:
                continue
            assignment = evaluation["slot_assignments"].get(slot_id, {})
            if "action_id" in assignment:
                needed.add(assignment["action_id"])
    actions = search.closure(needed)
    return _Selection(actions, tuple(row for row in selected.assignments if row[0] in actions))


def _exact(search, base, target_ids=(), *, optimize=False):
    """Enumerate action subsets and assignments under explicit hard bounds."""
    stats = _stats("exact_subsets", search.limits)
    excluded = _construction_limited(search, target_ids)
    if excluded:
        stats.update(status="input_construction_limited", excluded_question_ids=excluded,
                     reason_codes=["question_actions_excluded_by_declared_construction_limit"])
        return [], stats
    eligible = sorted(aid for aid in search.qset.actions if _eligible(search.qset, aid))
    stats["biologically_ineligible_actions"] = len(search.qset.actions) - len(eligible)
    if len(eligible) > search.limits.exact_action_cap or len(search.qset.resources) > search.limits.exact_resource_cap:
        stats["status"] = "exact_size_limit_exceeded"
        return [], stats
    remaining = sorted(set(eligible) - set(base.actions))
    results, found = [], 0
    keep = search.limits.max_alternatives_per_question
    maximum = min(search.limits.max_expansions, search.limits.exact_state_cap)
    for count in range(len(remaining) + 1):
        for extra in combinations(remaining, count):
            if stats["expansions"] >= maximum:
                stats["expansion_limit_reached"] = True
                break
            stats["expansions"] += 1
            chosen = tuple(sorted((*base.actions, *extra)))
            try:
                if search.closure(chosen) != chosen:
                    stats["constraint_rejections"] += 1
                    continue
            except PopulationError:
                stats["constraint_rejections"] += 1
                continue
            for selected in search.assignments(chosen, base, stats, maximum):
                evaluation = search.evaluate(selected)
                if not set(target_ids) <= set(evaluation["completed_question_ids"]):
                    continue
                selected = _completed_support(search, selected, base)
                evaluation = search.evaluate(selected)
                if any(row[1] == selected for row in results):
                    stats["deduplicated"] += 1
                    continue
                found += 1
                rank = ((-evaluation["weighted_completion"],) if optimize else ()) + (
                    search.account(selected)["total_cost"], _canonical(selected))
                results.append((rank, selected))
                results.sort(key=lambda row: row[0])
                del results[keep:]
        if stats["expansion_limit_reached"]:
            break
    stats["search_space_exhausted"] = not stats["expansion_limit_reached"]
    stats["exhaustive"] = stats["search_space_exhausted"]
    stats["status"] = "exhaustive" if stats["exhaustive"] else "expansion_limit_reached"
    if optimize and _construction_limited(search, search.qset.questions):
        stats.update(exhaustive=False, status="input_construction_limited",
                     excluded_question_ids=_construction_limited(search, search.qset.questions),
                     reason_codes=["objective_question_universe_has_unconstructed_actions"])
    stats["alternative_pruned"] = max(0, found - keep)
    stats["complete_alternatives_found"] = found
    return [row[1] for row in results], stats


@dataclass(frozen=True)
class _Branch:
    selection: _Selection
    # A branch remembers donor alternatives until all their slots are covered.
    donors: tuple = ()


def _donor_state(evaluation, qid, donor, arm):
    return next((row for row in evaluation["questions"][qid]["donors"]
                 if row["donor_id"] == donor and row["arm_id"] == arm), None)


def _bounded_step(search, branch, target_ids):
    evaluation = search.evaluate(branch.selection)
    complete = set(evaluation["completed_question_ids"])
    for qid in sorted(target_ids):
        if qid in complete:
            continue
        question = search.qset.questions[qid]
        picked = [row for row in branch.donors if row[0] == qid]
        options = [row for row in search.donor_bundles(qid) if row["eligible"]]
        unmet = []
        for _, arm, donor, _ in picked:
            donor_state = _donor_state(evaluation, qid, donor, arm)
            donor_options = next(row for row in options if row["donor_id"] == donor and row["arm_id"] == arm)
            missing = set(donor_state["unmet_slot_ids"]) if donor_state is not None else {
                row["slot_id"] for row in donor_options["slots"]}
            unmet.extend(row for row in donor_options["slots"] if row["slot_id"] in missing)
        if unmet:
            slot = min(unmet, key=lambda row: row["slot_id"])
            return "slot", qid, slot
        # Choose donor identities before filling their canonical slots. The
        # separate branches keep more expensive compatible alternatives alive.
        for arm in sorted(question["arms"], key=lambda row: row["arm_id"]):
            chosen = [row for row in picked if row[1] == arm["arm_id"]]
            if len(chosen) < arm["required_independent_donors"]:
                groups = {row[3] for row in picked}
                people = {row[2] for row in picked}
                candidates = [row for row in options if row["arm_id"] == arm["arm_id"]
                              and row["donor_id"] not in people and row["independence_group_id"] not in groups]
                return "donor", qid, sorted(candidates, key=lambda row: (row["donor_id"], row["arm_id"]))
        # If all chosen donor slots are covered but the joint evaluator rejects
        # reuse or covariate pairing, this branch is not a complete experiment.
        return "blocked", qid, []
    return "complete", None, []


def _unresolved(search, branch, target_ids):
    evaluation = search.evaluate(branch.selection)
    count = 0
    for qid in target_ids:
        if qid in evaluation["completed_question_ids"]:
            continue
        question = search.qset.questions[qid]
        options = search.donor_bundles(qid)
        picked = [row for row in branch.donors if row[0] == qid]
        for _, arm, donor, _ in picked:
            current = _donor_state(evaluation, qid, donor, arm)
            count += len(current["unmet_slot_ids"]) if current is not None else 1
        for arm in question["arms"]:
            missing_donors = arm["required_independent_donors"] - sum(row[1] == arm["arm_id"] for row in picked)
            lengths = [len(row["slots"]) for row in options if row["eligible"] and row["arm_id"] == arm["arm_id"]]
            count += max(0, missing_donors) * (min(lengths) if lengths else 1)
    return count


def _bounded(search, base, target_ids):
    stats = _stats("deterministic_beam", search.limits)
    stats["biologically_ineligible_actions"] = sum(not _eligible(search.qset, aid) for aid in search.qset.actions)
    beam, complete, seen = [_Branch(base)], {}, set()
    maximum = search.limits.max_expansions
    while beam and stats["expansions"] < maximum:
        next_beam = []
        for branch in beam:
            kind, qid, choices = _bounded_step(search, branch, target_ids)
            if kind == "complete":
                complete[branch.selection] = branch.selection
                continue
            if kind == "blocked":
                stats["constraint_rejections"] += 1
                continue
            branches = []
            if kind == "donor":
                for donor in choices:
                    if stats["expansions"] >= maximum:
                        stats["expansion_limit_reached"] = True
                        break
                    stats["expansions"] += 1
                    row = qid, donor["arm_id"], donor["donor_id"], donor["independence_group_id"]
                    branches.append(_Branch(branch.selection, tuple(sorted((*branch.donors, row)))))
            else:
                for aid in sorted(choices["action_ids"]):
                    if aid in branch.selection.actions:
                        continue
                    try:
                        actions = search.closure((*branch.selection.actions, aid))
                    except PopulationError:
                        stats["constraint_rejections"] += 1
                        continue
                    for selected in search.assignments(actions, branch.selection, stats, maximum):
                        branches.append(_Branch(selected, branch.donors))
            for candidate in branches:
                if candidate in seen:
                    stats["deduplicated"] += 1
                    continue
                seen.add(candidate)
                next_beam.append(candidate)
            if stats["expansions"] >= maximum:
                stats["expansion_limit_reached"] = True
                break
        next_beam.sort(key=lambda branch: (_unresolved(search, branch, target_ids),
                                          search.account(branch.selection)["total_cost"] - search.account(base)["total_cost"],
                                          _canonical(branch.selection), branch.donors))
        stats["beam_pruned"] += max(0, len(next_beam) - search.limits.beam_width)
        beam = next_beam[:search.limits.beam_width]
    # A final expansion can complete a bundle exactly at the expansion bound.
    for branch in beam:
        if set(target_ids) <= set(search.evaluate(branch.selection)["completed_question_ids"]):
            complete[branch.selection] = branch.selection
    values = sorted(complete, key=lambda selected: (search.account(selected)["total_cost"], _canonical(selected)))
    stats["complete_alternatives_found"] = len(values)
    stats["alternative_pruned"] = max(0, len(values) - search.limits.max_alternatives_per_question)
    stats["search_space_exhausted"] = not beam and not stats["expansion_limit_reached"]
    # Even an empty beam is not a certificate: only the explicit subset search
    # implements complete enumeration of all actions and resource assignments.
    stats["status"] = "expansion_limit_reached" if stats["expansion_limit_reached"] else (
        "beam_pruned" if stats["beam_pruned"] else "bounded_search_finished")
    return values[:search.limits.max_alternatives_per_question], stats


def _construct(search, base, target_ids):
    excluded = _construction_limited(search, target_ids)
    if excluded:
        stats = _stats("input_construction_limited", search.limits)
        stats.update(status="input_construction_limited", excluded_question_ids=excluded,
                     reason_codes=["question_actions_excluded_by_declared_construction_limit"])
        return [], stats
    eligible_count = sum(_eligible(search.qset, aid) for aid in search.qset.actions)
    if eligible_count <= search.limits.exact_action_cap and len(search.qset.resources) <= search.limits.exact_resource_cap:
        return _exact(search, base, target_ids)
    return _bounded(search, base, target_ids)


def _ratio(gain, cost, identity):
    return (0 if cost == 0 else 1, -Fraction(gain, cost) if cost else Fraction(0), -gain, cost, identity)


class _Coverage:
    """Coverage lookups owned by one revalidated planning snapshot."""

    def __init__(self, qset):
        self.actions = qset.actions
        carried, configurations = set(), {}
        for row in qset.records:
            if row["record_type"] == "comparison_donor_configuration" and row["carrier_status"] in {
                    "supported_carrier", "modeled_carrier"}:
                carried.add((row["donor_id"], row["design_configuration_id"]))
            elif row["record_type"] == "comparison_design_configuration" and row["role"] == "off_target":
                key = row["design_id"], row["sequence_configuration_id"]
                configurations.setdefault(key, set()).add(row["design_configuration_id"])
        self.carried = frozenset(carried)
        self.configurations = {key: frozenset(values) for key, values in configurations.items()}
        self.by_action = {}

    def units(self, action_ids):
        covered = set()
        for aid in action_ids:
            if aid not in self.by_action:
                action, units = self.actions[aid], set()
                if action["operation_type"] in {"locus_assay", "measurement_panel"}:
                    sequence_ids = {action["target_sequence_id"]["value"]}
                    sequence_ids.update(sid for row in action.get("target_sequences", [])
                                        for sid in row["sequence_configuration_ids"])
                    for sid in sequence_ids:
                        key = action["design_id"]["value"], sid
                        units.update(cid for cid in self.configurations.get(key, ())
                                     if (action["donor_id"], cid) in self.carried)
                self.by_action[aid] = frozenset(units)
            covered.update(self.by_action[aid])
        return covered


def _augmentations(search, base, trace):
    completed = set(search.evaluate(base)["completed_question_ids"])
    unique = {}
    for qid in sorted(search.qset.questions.keys() - completed):
        bundles, stats = _construct(search, base, (qid,))
        trace.append({"purpose": "optional_bundle", "question_ids": [qid], **stats})
        for selected in bundles:
            if search.evaluate(selected)["weighted_completion"] > search.evaluate(base)["weighted_completion"]:
                unique[selected] = selected
    return sorted(unique, key=_canonical)


def _one_addition(search, base, action_ids, trace, purpose):
    stats = _stats("resource_assignment_enumeration", search.limits)
    try:
        actions = search.closure((*base.actions, *action_ids))
    except PopulationError:
        stats["constraint_rejections"] += 1
        trace.append({"purpose": purpose, **stats})
        return []
    result = list(search.assignments(actions, base, stats, search.limits.max_expansions))
    stats["search_space_exhausted"] = not stats["expansion_limit_reached"]
    stats["status"] = "expansion_limit_reached" if stats["expansion_limit_reached"] else "finished"
    trace.append({"purpose": purpose, **stats})
    return result


def _select_optional(search, base, selector, seed, trace):
    selected, steps, rng = base, [], random.Random(seed)
    while True:
        before = search.evaluate(selected)
        cost_before = search.account(selected)["total_cost"]
        choices = []
        if selector in {"complete_bundle_greedy", "canonical_bundle", "random_bundle"}:
            candidates = _augmentations(search, selected, trace)
            for candidate in candidates:
                gain = search.evaluate(candidate)["weighted_completion"] - before["weighted_completion"]
                incremental = search.account(candidate)["total_cost"] - cost_before
                choices.append((_ratio(gain, incremental, _canonical(candidate)), candidate))
            if choices and selector == "canonical_bundle":
                choices.sort(key=lambda item: _canonical(item[1]))
            elif choices and selector == "random_bundle":
                choices.sort(key=lambda item: _canonical(item[1]))
                rng.shuffle(choices)
            else:
                choices.sort(key=lambda item: item[0])
        elif selector == "donor_first":
            filled = set(before["satisfied_slot_ids"])
            donor_actions = {}
            for qid in sorted(search.qset.questions.keys() - set(before["completed_question_ids"])):
                for donor in search.donor_bundles(qid):
                    if not donor["eligible"]:
                        continue
                    for slot in donor["slots"]:
                        if slot["slot_id"] not in filled:
                            donor_actions.setdefault(donor["donor_id"], set()).update(slot["action_ids"])
            for donor, aids in sorted(donor_actions.items()):
                remaining = aids - set(selected.actions)
                if not remaining:
                    continue
                for candidate in _one_addition(search, selected, remaining, trace, "donor_first"):
                    gain = len(set(search.evaluate(candidate)["satisfied_slot_ids"]) - filled)
                    incremental = search.account(candidate)["total_cost"] - cost_before
                    if gain > 0:
                        choices.append((_ratio(gain, incremental, (donor, _canonical(candidate))), candidate))
            choices.sort(key=lambda item: item[0])
        elif selector == "adapted_coverage":
            covered = search.coverage(selected)
            for aid, action in sorted(search.qset.actions.items()):
                if aid in selected.actions or not _eligible(search.qset, aid) or action["operation_type"] not in {"locus_assay", "measurement_panel"}:
                    continue
                for candidate in _one_addition(search, selected, (aid,), trace, "adapted_coverage"):
                    gain = len(search.coverage(candidate) - covered)
                    incremental = search.account(candidate)["total_cost"] - cost_before
                    if gain > 0:
                        choices.append((_ratio(gain, incremental, (aid, _canonical(candidate))), candidate))
            choices.sort(key=lambda item: item[0])
        if not choices:
            break
        previous = selected
        selected = choices[0][1]
        steps.append({"added_action_ids": sorted(set(selected.actions) - set(previous.actions)),
                      "candidate_order": [content_id("comparison_bundle", {
                          "action_ids": list(row.actions), "assignments": row.assignment_rows()}) for _, row in choices],
                      "resource_assignments": selected.assignment_rows(),
                      "incremental_cost": search.account(selected)["total_cost"] - cost_before,
                      "newly_completed_question_ids": sorted(set(search.evaluate(selected)["completed_question_ids"]) -
                                                             set(before["completed_question_ids"])),
                      "weighted_question_gain": search.evaluate(selected)["weighted_completion"] - before["weighted_completion"]})
    return selected, steps


def plan_experiments(question_set, *, budget=None, selector=None, limits=None,
                     random_seed=None, previous_plan_id=None):
    """Construct a reproducible plan; failed hard requirements remain explicit."""
    started = time.monotonic()
    return _plan_prepared(_reprepare(question_set), budget=budget, selector=selector,
                          limits=limits, random_seed=random_seed,
                          previous_plan_id=previous_plan_id, started=started)


def _planning_source_hash(qset):
    """Hash the complete scientific input without retaining serialized records."""
    source_config = {"mode": qset.mode, "coverage_basis": qset.coverage_basis,
                     "planning_policy": deepcopy(_config(qset)["policies"]["planning"])}
    include_predictions = _config(qset)["policies"]["comparison"]["include_predictions"]
    source_payload = {"identity_version": "question-records-and-planning-policy-v2",
                      "study": source_config, "sources": {
                        key: {field: value for field, value in row.items() if field not in {"path", "bytes"}}
                        for key, row in qset.study.get("sources", {}).items()
                        if row["role"] != "prediction_inputs" or include_predictions},
                      "records": sorted(qset.records, key=lambda row: row["record_id"])}
    return hashlib.sha256(canonical_bytes(source_payload)).hexdigest()


def _plan_prepared(qset, *, budget=None, selector=None, limits=None,
                   random_seed=None, previous_plan_id=None, source_sha256=None, started=None):
    """Plan from an owned snapshot; supported public entry points validate first."""
    started = time.monotonic() if started is None else started
    policy = _config(qset)["policies"]["planning"]
    budget = policy["total_budget"] if budget is None else budget
    selector = policy["algorithm"] if selector is None else selector
    random_seed = policy["seed"] if random_seed is None else random_seed
    limits = PlannerLimits.from_policy(policy) if limits is None else (
        PlannerLimits(**limits) if isinstance(limits, dict) else limits)
    _budget(budget)
    if not isinstance(limits, PlannerLimits):
        raise PopulationError("Planner limits must be explicit PlannerLimits or a matching object")
    if selector not in SELECTORS:
        raise PopulationError("Unsupported experiment selector: " + str(selector))
    if type(random_seed) is not int or random_seed < 0:
        raise PopulationError("planning.seed must be a nonnegative integer")
    if previous_plan_id is not None and (not isinstance(previous_plan_id, str) or not re.fullmatch(
            r"comparison_experiment_plan:v1:[0-9a-f]{64}", previous_plan_id)):
        raise PopulationError("previous_plan_id must name a canonical experiment-plan revision")
    if policy["fill_policy"] != "none":
        raise PopulationError("Only the explicit no-fill policy is implemented")
    search, empty, trace = _Search(qset, budget, limits), _Selection(), []
    baseline = search.evaluate(empty)
    required = sorted(qid for qid, question in qset.questions.items() if question["hard_required"])
    selected, seed_selection, alternatives, steps = empty, empty, [], []
    seed_report = {"required_question_ids": required, "status": "not_required", "total_cost": search.account(empty)["total_cost"],
                   "selected_action_ids": [], "resource_assignments": [], "alternatives": []}
    exhaustive, status = False, "heuristic_complete"
    if selector == "exact_tiny":
        alternatives, stats = _exact(search, empty, required, optimize=True)
        trace.append({"purpose": "exact_objective", "question_ids": required, **stats})
        exhaustive = stats["exhaustive"]
        if alternatives:
            selected = alternatives[0]
            status = "optimal" if exhaustive else "feasible_search_incomplete"
        else:
            status = "infeasible_proven" if exhaustive else "no_feasible_bundle_found"
        seed_report["status"] = "exact_joint_optimization"
        seed_report["total_cost"] = None
    else:
        if required and not set(required) <= set(baseline["completed_question_ids"]):
            alternatives, stats = _construct(search, empty, required)
            trace.append({"purpose": "mandatory_seed", "question_ids": required, **stats})
            seed_report["status"] = "feasible" if alternatives else (
                "infeasible_proven" if stats["exhaustive"] else "no_feasible_bundle_found")
            seed_report["alternatives"] = [{"selected_action_ids": list(row.actions),
                                             "resource_assignments": row.assignment_rows(),
                                             "total_cost": search.account(row)["total_cost"]} for row in alternatives]
            if alternatives:
                seed_selection = selected = alternatives[0]
            else:
                status = seed_report["status"]
        elif required:
            seed_report["status"] = "already_satisfied"
        seed_report.update({"total_cost": search.account(seed_selection)["total_cost"],
                            "selected_action_ids": list(seed_selection.actions),
                            "resource_assignments": seed_selection.assignment_rows()})
        if status == "heuristic_complete" and search.account(selected)["feasible"]:
            selected, steps = _select_optional(search, selected, selector, random_seed, trace)
    evaluation = search.evaluate(selected)
    accounting = search.account(selected)
    unmet_required = sorted(set(required) - set(evaluation["completed_question_ids"]))
    if not accounting["feasible"]:
        status = "infeasible_proven" if all(row["code"] == "budget_exceeded" for row in accounting["problems"]) else "invalid_plan"
    certificate = None
    if status == "infeasible_proven":
        certificate = ({"kind": "unavoidable_overhead", "minimum_cost": accounting["total_cost"],
                        "resource_ids": accounting["active_resource_ids"]} if not selected.actions and not accounting["feasible"]
                       else {"kind": "exhaustive_declared_action_and_instance_catalog",
                             "search_indexes": [i for i, row in enumerate(trace) if row["exhaustive"] and
                                                row["purpose"] in {"mandatory_seed", "exact_objective"}]})
    planned_states = {}
    for qid, question in qset.questions.items():
        state = deepcopy(evaluation["questions"][qid])
        if state["planning_status"] == "already_satisfied":
            state["planning_status"] = "already_satisfied"
        elif qid in evaluation["completed_question_ids"]:
            state["planning_status"] = "fully_planned"
        elif qid in unmet_required:
            state["planning_status"] = status if status in {"infeasible_proven", "no_feasible_bundle_found"} else "partially_planned"
        if question.get("construction_status") == "excluded_by_declared_limit":
            state["planning_status"] = "no_feasible_bundle_found"
            state["remaining_gaps"] = sorted(set(state.get("remaining_gaps", [])) | {"question_actions_excluded_by_declared_construction_limit"})
        state["evidence_status"] = baseline["questions"][qid]["evidence_status"]
        planned_states[qid] = state
    source_sha256 = _planning_source_hash(qset) if source_sha256 is None else source_sha256
    implementation_files = (Path(__file__), Path(_questions_api().__file__),
                            Path(__file__).with_name("comparison_contracts.py"))
    implementation_sha256 = hashlib.sha256(canonical_bytes({path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                                                          for path in implementation_files})).hexdigest()
    already_satisfied = sorted(qid for qid, state in evaluation["questions"].items()
                               if state["planning_status"] == "already_satisfied")
    existing_weight = sum(qset.questions[qid]["priority_weight"] for qid in already_satisfied)
    result = {"artifact_type": "comparison_experiment_plan", "schema_version": "1", "planner_version": PLANNER_VERSION,
              "source_sha256": source_sha256, "previous_plan_id": previous_plan_id,
              "implementation_sha256": implementation_sha256,
              "selector": selector, "random_seed": random_seed, "random_generator": "python.Random-MT19937",
              "limits": asdict(limits), "budget": budget, "fill_policy": "none",
              "status": status, "optimality_established": exhaustive and status == "optimal",
              "infeasibility_certificate": certificate,
              "meets_request": accounting["feasible"] and not unmet_required,
              "exploratory_partial": bool(unmet_required), "unmet_hard_required_question_ids": unmet_required,
              "construction_limited_question_ids": _construction_limited(search, qset.questions),
              "selected_action_ids": list(selected.actions), "resource_assignments": selected.assignment_rows(),
              "accounting": accounting, "mandatory_seed": seed_report, "selection_steps": steps,
              "searches": trace, "questions": planned_states, "requirements": evaluation,
              "objective_name": "weighted_question_completion",
              "weighted_completion": evaluation["weighted_completion"],
              "unweighted_completion": evaluation["unweighted_completion"],
              "already_satisfied_question_ids": already_satisfied,
              "existing_evidence_weight": existing_weight,
              "newly_planned_weight": evaluation["weighted_completion"] - existing_weight,
              "baseline_satisfied_question_ids": baseline["completed_question_ids"],
              "baseline_evidence_weight": baseline["weighted_completion"],
              "displaced_baseline_question_ids": sorted(set(baseline["completed_question_ids"]) - set(already_satisfied)),
              "coverage_definition": "one_unit_per_actual_carrier_off_target_design_configuration_v1",
              "covered_configuration_ids": sorted(search.coverage(selected))}
    result["plan_record"] = _plan_record(result)
    result["plan_id"] = result["plan_record"]["plan_id"]
    result["runtime_seconds"] = time.monotonic() - started
    result["process_peak_rss_bytes"] = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) * (1 if sys.platform == "darwin" else 1024)
    result["memory_measurement"] = "Cumulative peak for this Python process, including work before this plan"
    return result


def _plan_record(audit):
    payload = {key: value for key, value in audit.items() if key not in _AUDIT_HASH_EXCLUDED_FIELDS}
    names = ("planner_version", "source_sha256", "implementation_sha256", "selector", "random_seed", "limits",
             "budget", "fill_policy", "status", "optimality_established", "meets_request", "selected_action_ids",
             "resource_assignments", "weighted_completion", "unweighted_completion", "existing_evidence_weight", "newly_planned_weight")
    values = {key: deepcopy(audit[key]) for key in names}
    previous, seed_cost = audit["previous_plan_id"], audit["mandatory_seed"]["total_cost"]
    values.update({"previous_plan_id": {"value": previous, "reason": "First plan revision" if previous is None else None},
                   "mandatory_seed_cost": {"value": seed_cost, "reason": "Exact joint optimization has no separate seed" if seed_cost is None else None},
                   "question_states": [{"question_id": qid, "planning_status": state["planning_status"],
                                        "evidence_status": state["evidence_status"]} for qid, state in sorted(audit["questions"].items())],
                   "audit_sha256": hashlib.sha256(canonical_bytes(payload)).hexdigest()})
    values.update({key: deepcopy(audit["accounting"][key]) for key in (
        "total_cost", "cost_unit", "cost_scale", "active_resource_ids")})
    return make_record("experiment_plan", **values)


def _verify_plan_identity(plan):
    expected = _plan_record(plan)
    if plan.get("plan_record") != expected or plan.get("plan_id") != expected["plan_id"]:
        raise PopulationError("Plan revision identity or audit binding does not match its content")


def verify_plan_revision(question_set, plan):
    """Replay declared selection rules and recompute every result from inputs."""
    fresh = plan_experiments(question_set, budget=plan["budget"], selector=plan["selector"],
                             limits=plan["limits"], random_seed=plan["random_seed"],
                             previous_plan_id=plan["previous_plan_id"])
    return _compare_replayed_plan(fresh, plan)


def _verify_prepared_revision(qset, plan, *, source_sha256):
    fresh = _plan_prepared(qset, budget=plan["budget"], selector=plan["selector"],
                           limits=plan["limits"], random_seed=plan["random_seed"],
                           previous_plan_id=plan["previous_plan_id"], source_sha256=source_sha256)
    return _compare_replayed_plan(fresh, plan)


def _compare_replayed_plan(fresh, plan):
    expected = {key: value for key, value in fresh.items() if key not in _RUNTIME_FIELDS}
    supplied = {key: value for key, value in plan.items() if key not in _RUNTIME_FIELDS}
    if canonical_bytes(expected) != canonical_bytes(supplied):
        raise PopulationError("Plan differs from source-bound recomputation of selection, requirements or costs")
    return fresh


def write_plan_revision(plan, directory):
    """Write a content-addressed revision without overwriting a historical plan."""
    _verify_plan_identity(plan)
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    if plan.get("previous_plan_id") is not None:
        previous = directory / (plan["previous_plan_id"].rsplit(":", 1)[1] + ".json")
        if not previous.is_file():
            raise PopulationError("Previous plan revision must exist in the revision directory")
        prior = json.loads(previous.read_text())
        if prior.get("plan_id") != plan["previous_plan_id"]:
            raise PopulationError("Previous plan revision content has changed")
        _verify_plan_identity(prior)
    path = directory / (plan["plan_id"].rsplit(":", 1)[1] + ".json")
    data = canonical_bytes({key: value for key, value in plan.items() if key not in _RUNTIME_FIELDS}) + b"\n"
    try:
        with path.open("xb") as handle:
            handle.write(data)
    except FileExistsError:
        if path.read_bytes() != data:
            raise PopulationError("Historical plan revision has changed")
    return path
