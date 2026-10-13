"""Independent exhaustive checks for small, source-bound synthetic plans.

This module deliberately imports no package code. It supports explicit synthetic
sequence truth and prospective actions, with configuration/any-resolved genotype
predicates and literal single-base alleles. Unsupported semantics fail rather than acquire
an exhaustive label. The action, instance and state bounds are hard limits.
"""
from collections import Counter
from itertools import combinations, product


ACTION_CAP = 16
INSTANCE_CAP = 64
STATE_CAP = 300_000


def _value(value):
    return value.get("value") if isinstance(value, dict) else value


class TinyOracle:
    def __init__(self, study, records):
        if study["mode"] != "synthetic_test" or study["coverage_basis"] != "synthetic_truth":
            raise ValueError("Independent tiny oracle is restricted to declared synthetic truth")
        self.study = study
        self.index = {row["record_id"]: row for row in records}
        self.records = list(records)
        self.actions = self.kind("action")
        self.resources = self.kind("resource")
        self.questions = self.kind("question")
        self.slots = self.kind("requirement_slot")
        self.donors = self.kind("donor")
        if len(self.actions) > ACTION_CAP or len(self.resources) > INSTANCE_CAP:
            raise ValueError("Independent oracle hard action/resource bound exceeded")
        if self.kind("observation_revision"):
            raise ValueError("This independent oracle does not model existing observations")
        for question in self.questions.values():
            predicate = question["genotype_predicate"]
            if predicate["kind"] not in {"configuration", "any_resolved", "allele"}:
                raise ValueError("Independent oracle predicate is outside its declared scope")
            if predicate["kind"] == "allele" and (len(predicate["variants"]) != 1 or
                    len(predicate["variants"][0]["reference"]) != 1 or
                    len(predicate["variants"][0]["alternates"][predicate["variants"][0]["selected_alt_index"] - 1]) != 1 or
                    predicate["phase_requirement"] != "not_required"):
                raise ValueError("Independent allele oracle requires one unphased single-base allele")
            if question["matching_covariates"] or question["template"] == "evidence_replication":
                raise ValueError("Independent oracle does not certify covariate or historical replication semantics")
            if not question["evidence_reuse"]["allow_across_questions"] and len(self.questions) > 1:
                raise ValueError("Independent oracle does not certify cross-question measurement exclusion")
        self.action_eligible = {aid: self.eligible(action) for aid, action in self.actions.items()}
        for aid in self.actions:
            self.dependencies(aid, set())

    def kind(self, kind):
        return {key: row for key, row in self.index.items() if row["record_type"] == "comparison_" + kind}

    def dependencies(self, aid, visiting):
        if aid in visiting:
            raise ValueError("Cyclic action prerequisites")
        if aid not in self.actions:
            raise ValueError("Missing action prerequisite")
        result = {aid}
        for child in self.actions[aid]["dependency_action_ids"]:
            result.update(self.dependencies(child, visiting | {aid}))
        return result

    def sequences(self, donor, locus):
        region = self.index[locus]
        copies = {}
        for occurrence in self.kind("sequence_occurrence").values():
            if occurrence["donor_id"] != donor:
                continue
            sequence = self.index[occurrence["sequence_configuration_id"]]
            if sequence["support_status"] != "supported" or sequence["contig"] != region["contig"]:
                continue
            positions = {point["pos0"] for point in sequence["path"] if point["kind"] == "reference"}
            if not set(range(region["start"], region["end"])) <= positions:
                continue
            copies.setdefault(occurrence["copy_index"], set()).add(sequence["record_id"])
        return copies

    def membership(self, question, donor, arm):
        predicate = question["genotype_predicate"]
        copies = self.sequences(donor, predicate["locus_id"])
        if set(copies) != {0, 1}:
            return False
        if arm["membership"] == "all_resolved":
            return True
        if predicate["kind"] == "allele":
            variant = predicate["variants"][0]
            alternate = variant["alternates"][variant["selected_alt_index"] - 1]
            dosage = 0
            for possibilities in copies.values():
                if len(possibilities) != 1:
                    return False
                sequence = self.index[next(iter(possibilities))]
                if sequence["strand"] != "+":
                    raise ValueError("Independent allele oracle requires genomic-forward paths")
                base = next((point["base"] for point in sequence["path"] if point["kind"] == "reference"
                             and point["pos0"] == variant["pos0"]), None)
                if base not in {variant["reference"], alternate}:
                    return False
                dosage += base == alternate
        else:
            target = _value(predicate["sequence_configuration_id"])
            if target is None:
                return False
            dosage = sum(target in value for value in copies.values())
        if arm["membership"] == "noncarrier":
            return dosage == 0
        wanted = _value(predicate["dosage"])
        return dosage > 0 if wanted is None else dosage == wanted

    def compatible(self, action, resource):
        compatibility = resource["compatibility"]
        for field, actual in (("donor_ids", action["donor_id"]), ("design_ids", _value(action["design_id"])),
                              ("context_ids", action["context_id"]), ("batch_ids", action["batch_id"]),
                              ("assay_capability_ids", action["assay_capability_id"])):
            if compatibility[field] and actual not in compatibility[field]:
                return False
        for field in ("endpoint_ids", "locus_ids"):
            if compatibility[field] and not set(action[field]) <= set(compatibility[field]):
                return False
        if compatibility["timepoints"]:
            timepoint = self.index[action["context_id"]]["measurement_timepoint"]
            number = _value(timepoint["value"])
            label = None if number is None else str(number) + " " + timepoint["unit"]
            if label not in compatibility["timepoints"]:
                return False
        return True

    def eligible(self, action):
        specimen = self.index[action["specimen_id"]]
        assay = self.index[action["assay_capability_id"]]
        if specimen["donor_id"] != action["donor_id"] or specimen["availability"] not in {"verified", "simulation_assumption"}:
            return False
        if _value(specimen["permitted_use"]) is None or action["context_id"] not in specimen["context_ids"]:
            return False
        if assay["feasibility"] not in {"verified", "simulation_assumption"} or action["context_id"] not in assay["context_ids"]:
            return False
        if not set(action["endpoint_ids"]) <= set(assay["endpoint_ids"]) or not set(action["locus_ids"]) <= set(assay["locus_ids"]):
            return False
        if len(action["locus_ids"]) > assay["max_targets"]:
            return False
        if action["operation_type"] in {"measurement_panel", "locus_assay", "intended_assay", "mock_control"}:
            for target in action.get("target_sequences", []):
                copies = self.sequences(action["donor_id"], target["locus_id"])
                actual = set().union(*copies.values()) if copies else set()
                if set(copies) != {0, 1} or set(target["sequence_configuration_ids"]) != actual:
                    return False
            if not set(assay["qc_requirements"]) <= set(action.get("performed_qc", [])):
                return False
        for design in [_value(action["design_id"])]:
            if design is not None:
                for locus in self.index[design]["intended_locus_ids"]:
                    if set(self.sequences(action["donor_id"], locus)) != {0, 1}:
                        return False
        return True

    def assignment_options(self, selected):
        groups = {}
        for aid in selected:
            for demand in self.actions[aid]["resource_demands"]:
                if "resource_id" in demand:
                    continue
                key = aid, demand["assignment_key"]
                allowed = set(demand["eligible_resource_ids"])
                groups[key] = groups.get(key, allowed) & allowed
        keys = sorted(groups)
        for values in product(*(sorted(groups[key]) for key in keys)):
            yield {key: value for key, value in zip(keys, values)}

    def cost(self, selected, assignments):
        if any(not self.action_eligible[aid] or not self.dependencies(aid, set()) <= set(selected) for aid in selected):
            return None
        for aid in selected:
            action = self.actions[aid]
            available = [self.actions[other] for other in self.dependencies(aid, set())]
            for requirement in action["control_requirements"]:
                valid = [other for other in available if other["donor_id"] == action["donor_id"]
                         and all(other[field] == action[field] for field in (
                             "specimen_id", "context_id", "batch_id", "experiment_id", "replicate_position"))
                         and other["assay_capability_id"] == requirement["assay_capability_id"]
                         and (requirement["role"] in other.get("performed_control_roles", [])
                              or (requirement["role"] == "mock" and other["operation_type"] == "mock_control"))
                         and set(requirement["endpoint_ids"]) <= set(other["endpoint_ids"])
                         and set(requirement["locus_ids"]) <= set(other["locus_ids"])
                         and set(requirement["required_qc"]) <= set(other.get("performed_qc", []))]
                if not valid:
                    return None
        consumption = Counter()
        for aid in selected:
            for demand in self.actions[aid]["resource_demands"]:
                rid = demand.get("resource_id")
                if rid is None:
                    rid = assignments.get((aid, demand["assignment_key"]))
                    if rid not in demand["eligible_resource_ids"]:
                        return None
                if rid not in self.resources or not self.compatible(self.actions[aid], self.resources[rid]):
                    return None
                consumption[rid, demand["dimension"]] += demand["quantity"]
        for (rid, dimension), amount in consumption.items():
            caps = {cap["dimension"]: cap["capacity"] for cap in self.resources[rid]["capacities"]}
            if dimension not in caps or (caps[dimension]["kind"] == "finite" and amount > caps[dimension]["quantity"]):
                return None
        active = {rid for rid, _ in consumption}
        active.update(rid for rid, resource in self.resources.items() if resource["activation_rule"] == "always"
                      or (resource["activation_rule"] == "nonempty_plan" and selected))
        if len({(row["cost_unit"], row["cost_scale"]) for row in self.resources.values()}) > 1:
            return None
        return sum(self.resources[rid]["cost"] for rid in active)

    def action_matches(self, action, slot):
        if action["action_id"] not in slot["allowed_action_ids"]:
            return False
        if (action["donor_id"] != slot["donor_id"] or action["context_id"] != slot["context_id"] or
                _value(action["design_id"]) != _value(slot["design_id"]) or
                action["replicate_position"] != slot["replicate_position"]):
            return False
        if slot["endpoint_id"] not in action["endpoint_ids"]:
            return False
        locus = _value(slot["locus_id"])
        if locus is not None and locus not in action["locus_ids"]:
            return False
        pairs = {(row["endpoint_id"], row["locus_id"]) for row in action.get("endpoint_locus_pairs", [])}
        if pairs and (slot["endpoint_id"], locus) not in pairs:
            return False
        role = slot["control_role"]
        if role == "measurement":
            if action["operation_type"] not in {"measurement_panel", "locus_assay", "intended_assay", "functional_assay"}:
                return False
        elif role not in set(action.get("performed_control_roles", [])) | ({"mock"} if action["operation_type"] == "mock_control" else set()):
            return False
        return True

    def complete_slot_matrix(self, question, slots):
        present = {(row["design_id"]["value"], row["endpoint_id"], row["locus_id"]["value"],
                    row["context_id"], row["replicate_position"], row["control_role"]) for row in slots}
        pairs = question.get("endpoint_locus_pairs") or [
            {"endpoint_id": endpoint, "locus_id": locus}
            for endpoint in question["endpoint_ids"] for locus in question["locus_ids"]]
        intended = {locus for did in question["design_ids"] for locus in self.index[did]["intended_locus_ids"]}
        for did in question["design_ids"]:
            own_loci = set(self.index[did]["intended_locus_ids"])
            if own_loci - {pair["locus_id"] for pair in pairs}:
                return False
            for pair in pairs:
                eid, lid = pair["endpoint_id"], pair["locus_id"]
                if lid in intended and lid not in own_loci:
                    continue
                for context in question["context_ids"]:
                    for repeat in range(1, question["technical_replicates"] + 1):
                        measure_role = "intended_activity" if lid in own_loci else "measurement"
                        if any((did, eid, lid, context, repeat, role) not in present for role in (measure_role, "assay_qc")):
                            return False
                        if not any((mock_design, eid, lid, context, repeat, "mock") in present for mock_design in (None, did)):
                            return False
        return True

    def donor_choices(self, question, selected):
        arms = sorted(question["arms"], key=lambda row: row["arm_id"])
        per_arm = []
        for arm in arms:
            eligible = []
            for donor, record in sorted(self.donors.items()):
                group = _value(record["independence_group_id"])
                if group is None or not self.membership(question, donor, arm):
                    continue
                slots = [row for row in self.slots.values() if row["question_id"] == question["record_id"]
                         and row["donor_id"] == donor and row["arm_id"] == arm["arm_id"]]
                if not slots or not self.complete_slot_matrix(question, slots):
                    continue
                choices = {slot["slot_id"]: [aid for aid in selected if self.action_matches(self.actions[aid], slot)] for slot in slots}
                if any(not aids for aids in choices.values()):
                    continue
                eligible.append((donor, group, choices))
            per_arm.append(list(combinations(eligible, arm["required_independent_donors"])))
        for groups in product(*per_arm):
            donors = [row for arm in groups for row in arm]
            if len({row[0] for row in donors}) != len(donors) or len({row[1] for row in donors}) != len(donors):
                continue
            yield donors

    def question_complete(self, question, selected):
        return next(self.donor_choices(question, selected), None) is not None

    def evaluate(self, selected, assignments):
        cost = self.cost(selected, assignments)
        completed = sorted(qid for qid, question in self.questions.items() if self.question_complete(question, selected)) if cost is not None else []
        return {"total_cost": cost, "completed_question_ids": completed,
                "weighted_completion": sum(self.questions[qid]["priority_weight"] for qid in completed),
                "unweighted_completion": len(completed)}

    def has_unrelated_actions(self, selected, assignments, completed):
        """Try deletions independently; do not use production slot assignments."""
        selected = set(selected)
        for aid in selected:
            drop = {aid}
            while True:
                expanded = drop | {other for other in selected
                                   if set(self.actions[other]["dependency_action_ids"]) & drop}
                if expanded == drop:
                    break
                drop = expanded
            remaining = selected - drop
            result = self.evaluate(remaining, {key: rid for key, rid in assignments.items() if key[0] in remaining})
            if result["total_cost"] is not None and result["completed_question_ids"] == completed:
                return True
        return False

    def solve(self, budget):
        actions = sorted(self.actions)
        required = {qid for qid, question in self.questions.items() if question["hard_required"]}
        examined, best, minimum_full_cost = 0, None, None
        for count in range(len(actions) + 1):
            for selected in combinations(actions, count):
                for assignments in self.assignment_options(selected):
                    examined += 1
                    if examined > STATE_CAP:
                        return {"exhaustive": False, "status": "state_limit_reached", "states": examined - 1,
                                "best": best, "action_cap": ACTION_CAP, "instance_cap": INSTANCE_CAP, "state_cap": STATE_CAP}
                    result = self.evaluate(selected, assignments)
                    if result["total_cost"] is None:
                        continue
                    if len(result["completed_question_ids"]) == len(self.questions):
                        minimum_full_cost = result["total_cost"] if minimum_full_cost is None else min(minimum_full_cost, result["total_cost"])
                    if result["total_cost"] > budget or not required <= set(result["completed_question_ids"]):
                        continue
                    if self.has_unrelated_actions(selected, assignments, result["completed_question_ids"]):
                        continue
                    rows = [{"action_id": key[0], "assignment_key": key[1], "resource_id": value} for key, value in sorted(assignments.items())]
                    candidate = {**result, "selected_action_ids": list(selected), "resource_assignments": rows}
                    key = (-candidate["weighted_completion"], candidate["total_cost"], tuple(selected), tuple(sorted(assignments.items())))
                    if best is None or key < best["rank"]:
                        best = {**candidate, "rank": key}
        if best is not None:
            best.pop("rank")
        return {"exhaustive": True, "status": "optimal" if best is not None else "infeasible_proven",
                "states": examined, "best": best, "minimum_full_cost": minimum_full_cost,
                "action_cap": ACTION_CAP, "instance_cap": INSTANCE_CAP, "state_cap": STATE_CAP}
