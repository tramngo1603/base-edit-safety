"""Select explicit assays, charging activation once per biological donor."""
from __future__ import annotations

from fractions import Fraction
from itertools import combinations
import random

from .common import PopulationError
from .contracts import ALGORITHMS, validate_request, validate_schema, planner_schema


def eligible_actions(request):
    donors = {d["id"]: d for d in request["donors"]}
    result = []
    for action in request["actions"]:
        donor = donors[action["donor_id"]]
        allowed = donor["use_permitted"] and (
            donor["availability"] == "verified_available" if request["mode"] == "research_sample_planning"
            else donor["simulation_eligible"])
        if action["eligibility"] == "eligible" and allowed:
            result.append(action)
    return sorted(result, key=lambda a: a["id"])


def cost(request, selected):
    donors = {a["donor_id"] for a in selected}
    return sum(a["assay_cost"] for a in selected) + sum(
        d["activation_cost"] for d in request["donors"] if d["id"] in donors)


def feasible(request, selected):
    b = request["budget"]
    return cost(request, selected) <= b["max_cost"] and (
        b["max_donors"] is None or len({a["donor_id"] for a in selected}) <= b["max_donors"]) and (
        b["max_assays"] is None or len(selected) <= b["max_assays"])


def coverage(request, selected):
    existing = {u["id"]: set() for u in request["units"]}
    planned = {u["id"]: set() for u in request["units"]}
    for e in request["existing_usable_evidence"]:
        existing[e["unit_id"]].add(e["donor_id"])
    for action in selected:
        for uid in action["covered_units"]:
            planned[uid].add(action["donor_id"])
    return existing, planned


def objective(request, selected):
    existing, planned = coverage(request, selected)
    return sum((Fraction(u["weight"] * min(u["required_donors"],
                    len(existing[u["id"]] | planned[u["id"]])), u["required_donors"])
                for u in request["units"]), Fraction())


def exact_tiny(request, actions):
    """Find the exact optimum by exhaustive selection with an independent objective calculation."""
    if len(actions) > 12:
        raise PopulationError("exact_tiny accepts no more than 12 eligible actions")
    activations = {d["id"]: d["activation_cost"] for d in request["donors"]}
    best, best_key = [], None
    for size in range(len(actions) + 1):
        for subset in combinations(actions, size):
            people = {a["donor_id"] for a in subset}
            total = sum(a["assay_cost"] for a in subset) + sum(activations[d] for d in people)
            budget = request["budget"]
            if total > budget["max_cost"] or (budget["max_donors"] is not None and len(people) > budget["max_donors"]) \
                    or (budget["max_assays"] is not None and size > budget["max_assays"]):
                continue
            value = Fraction()
            for unit in request["units"]:
                witnesses = {e["donor_id"] for e in request["existing_usable_evidence"] if e["unit_id"] == unit["id"]}
                witnesses.update(a["donor_id"] for a in subset if unit["id"] in a["covered_units"])
                value += Fraction(unit["weight"] * min(unit["required_donors"], len(witnesses)), unit["required_donors"])
            key = (-value, total, len(people), size, tuple(a["id"] for a in subset))
            if best_key is None or key < best_key:
                best, best_key = list(subset), key
    return best


def select(request, algorithm="greedy_cost", seed=None):
    validate_request(request)
    if algorithm not in ALGORITHMS:
        raise PopulationError("Unsupported algorithm")
    if algorithm == "random_order" and type(seed) is not int:
        raise PopulationError("random_order requires an integer seed")
    actions = eligible_actions(request)
    order = [a["id"] for a in actions]
    selected = []
    if algorithm == "exact_tiny":
        selected = exact_tiny(request, actions)
    elif algorithm in {"canonical_order", "random_order"}:
        if algorithm == "random_order":
            random.Random(seed).shuffle(actions)
        order = [a["id"] for a in actions]
        for a in actions:
            if feasible(request, selected + [a]) and objective(request, selected + [a]) > objective(request, selected):
                selected.append(a)
    else:
        remaining = list(actions)
        while remaining:
            choices = []
            current = objective(request, selected)
            for a in remaining:
                if not feasible(request, selected + [a]):
                    continue
                gain = objective(request, selected + [a]) - current
                if gain <= 0:
                    continue
                incremental = cost(request, selected + [a]) - cost(request, selected)
                if algorithm == "greedy_coverage":
                    key = (-gain, incremental, a["id"])
                else:
                    key = (0 if incremental == 0 else 1,
                           -gain / incremental if incremental else Fraction(), -gain, incremental, a["id"])
                choices.append((key, a))
            if not choices:
                break
            chosen = min(choices, key=lambda item: item[0])[1]
            selected.append(chosen)
            remaining.remove(chosen)
    return selected, order


def plan(request, algorithm="greedy_cost", seed=None):
    selected, order = select(request, algorithm, seed)
    before, after = objective(request, []), objective(request, selected)
    existing, planned = coverage(request, selected)
    units = []
    for u in sorted(request["units"], key=lambda u: u["id"]):
        count = len(existing[u["id"]] | planned[u["id"]])
        units.append({"unit_id": u["id"], "required_donors": u["required_donors"],
                      "existing_donors": sorted(existing[u["id"]]), "planned_donors": sorted(planned[u["id"]]),
                      "combined_donor_count": count, "fully_satisfied": count >= u["required_donors"]})
    if selected:
        status = "planned"
    elif all(u["fully_satisfied"] for u in units):
        status = "no_remaining_demand"
    elif not any(feasible(request, [a]) for a in eligible_actions(request)):
        status = "no_feasible_action"
    else:
        status = "no_positive_gain"
    result = {"record_type": "donor_panel_result", "schema_version": "1",
              **{key: request[key] for key in ("mode", "coverage_basis", "protocol_id", "universe_id", "cost_unit")},
              "algorithm": algorithm, "plan_status": status,
              "selected_actions": sorted(a["id"] for a in selected),
              "selected_donors": sorted({a["donor_id"] for a in selected}),
              "total_cost": cost(request, selected), "unused_budget": request["budget"]["max_cost"] - cost(request, selected),
              "objective_numerator": after.numerator, "objective_denominator": after.denominator,
              "gain_numerator": (after - before).numerator, "gain_denominator": (after - before).denominator,
              "fully_covered_units": [u["unit_id"] for u in units if u["fully_satisfied"]],
              "uncovered_units": [u["unit_id"] for u in units if not u["fully_satisfied"]], "unit_coverage": units}
    validate_schema(result, planner_schema("result"))
    return result, order
