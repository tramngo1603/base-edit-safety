"""Build and evaluate complete prospective questions from exact donor evidence."""
from __future__ import annotations

from collections import Counter, defaultdict
from copy import deepcopy
from dataclasses import dataclass, field
from itertools import combinations

from .comparison_outcomes import outcome_reasons
from .common import PopulationError, ResourceLimit
from .comparison_contracts import active_observations, validate_catalog
from .comparison_identity import canonical_bytes, content_id
from .comparison_provenance import verify_sources
from .identity import reverse_complement


class _QuestionIndex(dict):
    def __init__(self, records):
        super().__init__((row["record_id"], row) for row in records)
        self.kinds = defaultdict(dict)
        self.occurrences = defaultdict(list)
        self.support = defaultdict(set)
        self.assessments = defaultdict(list)
        self.copy_cache = {}
        for row in records:
            self.kinds[row["record_type"]][row["record_id"]] = row
            if row["record_type"] == "comparison_sequence_occurrence":
                self.occurrences[(row["donor_id"], self[row["sequence_configuration_id"]]["locus_id"])].append(row)
            elif row["record_type"] == "comparison_donor_configuration":
                self.support[row["sequence_occurrence_id"]].add(row["carrier_status"].split("_", 1)[0])
            elif row["record_type"] == "comparison_intended_assessment":
                self.assessments[(row["donor_id"], row["locus_id"] )].append(row)


def _kind(index, name):
    if isinstance(index, _QuestionIndex):
        return index.kinds["comparison_" + name]
    return {key: row for key, row in index.items() if row["record_type"] == "comparison_" + name}


def _pairs(record):
    return record.get("endpoint_locus_pairs", [{"endpoint_id": endpoint, "locus_id": locus}
        for endpoint in record["endpoint_ids"] for locus in record["locus_ids"]])


def required_design_measurements(index, question, design_id):
    """Unscoped question pairs are required under every declared design.

    A locus can be intended for one design and a common candidate for another.
    Only explicit intended endpoint pairs at this design's own loci receive the
    intended-activity role; every other declared pair remains a measurement.
    """
    own_intended = set(index[design_id]["intended_locus_ids"])
    intended_pairs = {(endpoint, locus) for requirement in question["control_requirements"]
                      if requirement["role"] == "intended_activity"
                      for endpoint in requirement["endpoint_ids"] for locus in requirement["locus_ids"]}
    return [{**pair, "control_role": "intended_activity" if
             pair["locus_id"] in own_intended and (pair["endpoint_id"], pair["locus_id"]) in intended_pairs
             else "measurement"} for pair in _pairs(question)]


def _biological_slot(slot):
    return (slot["donor_id"], slot["design_id"]["value"], slot["context_id"], slot["endpoint_id"],
            slot["locus_id"]["value"], slot["control_role"])


def _actual_copies(index, donor_id, locus_id, coverage_basis):
    """Require source-established support for both chromosome copies at a locus."""
    supported = {"supported"} if coverage_basis != "retained_callset_model" else {"supported", "modeled"}
    if isinstance(index, _QuestionIndex) and (donor_id, locus_id, coverage_basis) in index.copy_cache:
        return index.copy_cache[(donor_id, locus_id, coverage_basis)]
    occurrences = index.occurrences[(donor_id, locus_id)] if isinstance(index, _QuestionIndex) else [row for row in index.values() if row["record_type"] == "comparison_sequence_occurrence"
                   and row["donor_id"] == donor_id and index[row["sequence_configuration_id"]]["locus_id"] == locus_id]
    assignments, reasons = {}, []
    for occurrence in occurrences:
        sequence_id = occurrence["sequence_configuration_id"]
        statuses = set(index.support[occurrence["record_id"]]) if isinstance(index, _QuestionIndex) else {row["carrier_status"].split("_", 1)[0] for row in index.values()
            if row["record_type"] == "comparison_donor_configuration" and row["sequence_occurrence_id"] == occurrence["record_id"]}
        assessment_rows = index.assessments[(donor_id, locus_id)] if isinstance(index, _QuestionIndex) else _kind(index, "intended_assessment").values()
        assessments = [row for row in assessment_rows if
                       row["donor_id"] == donor_id and row["locus_id"] == locus_id
                       and row["sequence_configuration_id"]["value"] == sequence_id]
        statuses.update(row["support_status"] for row in assessments)
        sequence = index[sequence_id]
        if not statuses & supported or statuses & {"unresolved", "not"} or occurrence["phase_status"] in {"putative", "unresolved"} or set(sequence["sequence"]) - set("ACGT"):
            continue
        previous = assignments.get(occurrence["copy_index"])
        if previous is not None and previous["sequence"]["record_id"] != sequence_id:
            reasons.append("conflicting_actual_copy_sequences")
        else:
            assignments[occurrence["copy_index"]] = {"sequence": sequence, "occurrence": occurrence}
    if set(assignments) != {0, 1}:
        reasons.append("both_actual_copies_are_not_resolved")
    result = assignments, sorted(set(reasons))
    if isinstance(index, _QuestionIndex):
        index.copy_cache[(donor_id, locus_id, coverage_basis)] = result
    return result


def _allele_on_path(sequence, variant):
    """Read a literal allele through its actual reference and insertion path."""
    contig, start, reference = variant["contig"], variant["pos0"], variant["reference"]
    if sequence["contig"] != contig:
        return None
    path = sequence["path"] if sequence["strand"] == "+" else [
        {**point, "base": reverse_complement(point["base"])} for point in reversed(sequence["path"])]
    end = start + len(reference)
    positions = [point["pos0"] for point in path if point["kind"] == "reference"]
    if not positions or min(positions) > start or max(positions) < end - 1:
        return None
    observed = "".join(point["base"] for point in path if
        point["kind"] == "reference" and start <= point["pos0"] < end or
        point["kind"] == "insertion" and start < point["boundary0"] <= end)
    return observed == variant["alternates"][variant["selected_alt_index"] - 1]


def donor_membership(index, question, donor_id, mode, coverage_basis):
    predicate = question["genotype_predicate"]
    reasons = []
    if predicate["coverage_basis"] != coverage_basis:
        reasons.append("genotype_coverage_basis_mismatch")
    copies, copy_reasons = _actual_copies(index, donor_id, predicate["locus_id"], coverage_basis)
    reasons.extend(copy_reasons)
    donor = index[donor_id]
    group = donor["independence_group_id"]["value"]
    assumptions = []
    if group is None:
        if mode != "research_sample_planning" and question["relatedness_policy"] == "declared_simulation_assumption":
            group = "declared-simulation-independent:" + donor_id
            assumptions.append("unknown_relatedness_assumed_independent_for_simulation")
        else:
            reasons.append("independence_group_unresolved")
    if mode == "research_sample_planning" and question["relatedness_policy"] != "strict":
        reasons.append("research_mode_requires_strict_relatedness")
    covariates = {item["name"]: item["value"]["value"] for item in donor.get("covariates", [])}
    if any(covariates.get(name) is None for name in question["matching_covariates"]):
        reasons.append("required_matching_covariate_unresolved")
    membership, dosage = None, None
    if not copy_reasons:
        if predicate["kind"] == "any_resolved":
            membership = "all_resolved"
        else:
            states = []
            for copy in (0, 1):
                sequence = copies[copy]["sequence"]
                if predicate["kind"] == "configuration":
                    target = predicate["sequence_configuration_id"]["value"]
                    states.append([sequence["record_id"] == target] if target is not None else [None])
                else:
                    states.append([_allele_on_path(sequence, variant) for variant in predicate["variants"]])
            if not states[0] or any(value is None for state in states for value in state):
                reasons.append("genotype_predicate_not_resolved_at_actual_path")
            else:
                phase = predicate["phase_requirement"]
                occurrences = [copies[copy]["occurrence"] for copy in (0, 1)]
                if phase != "not_required" and any(row["phase_status"] not in {"confirmed", "homozygous"} for row in occurrences):
                    reasons.append("required_haplotype_phase_unresolved")
                if phase == "confirmed_trans":
                    if len(states[0]) != 2:
                        reasons.append("trans_predicate_requires_two_alleles")
                    matched = len(states[0]) == 2 and (states[0] == [True, False] and states[1] == [False, True] or states[0] == [False, True] and states[1] == [True, False])
                    dosage = 1 if matched else 0
                else:
                    dosage = sum(all(state) for state in states)
                required = predicate["dosage"]["value"]
                membership = "carrier" if dosage > 0 and (required is None or dosage == required) else "noncarrier" if dosage == 0 else None
                if membership is None:
                    reasons.append("allele_dosage_outside_declared_carrier_group")
    roles = {requirement["role"] for requirement in question["control_requirements"]}
    for role in ("intended_activity", "mock", "assay_qc"):
        if role not in roles:
            reasons.append("mandatory_question_control_not_declared:" + role)
    allowed = {"supported", "modeled"} if coverage_basis == "retained_callset_model" else {"supported"}
    for design_id in question["design_ids"]:
        for locus in index[design_id]["intended_locus_ids"]:
            local = index.assessments[(donor_id, locus)] if isinstance(index, _QuestionIndex) else _kind(index, "intended_assessment").values()
            assessments = [row for row in local if (row["donor_id"], row["design_id"], row["locus_id"]) == (donor_id, design_id, locus)]
            if not assessments or any(row["support_status"] not in allowed or row["sequence_configuration_id"]["value"] is None for row in assessments):
                reasons.append("required_intended_locus_unresolved:" + locus)
            if _actual_copies(index, donor_id, locus, coverage_basis)[1]:
                reasons.append("required_intended_copies_unresolved:" + locus)
    return {"membership": membership, "dosage": dosage, "independence_group_id": group,
            "eligible": not reasons, "reason_codes": sorted(set(reasons)), "assumptions": assumptions,
            "covariates": {name: covariates.get(name) for name in question["matching_covariates"]}}


@dataclass
class QuestionSet:
    study: dict
    records: list
    index: dict
    questions: dict
    slots: dict
    actions: dict
    resources: dict
    eligibility: dict = field(default_factory=dict)
    donor_eligibility: dict = field(default_factory=dict)
    options: dict = field(default_factory=dict)
    observation_tokens: dict = field(default_factory=dict)
    evidence_rejections: dict = field(default_factory=dict)

    @property
    def mode(self):
        return self.study["config"]["mode"]

    @property
    def coverage_basis(self):
        return self.study["config"]["coverage_basis"]


def dependency_closure(question_set, action_ids):
    result, visiting = set(), set()
    def visit(identifier):
        if identifier in visiting:
            raise PopulationError("Cyclic experimental action dependencies")
        if identifier in result:
            return
        if identifier not in question_set.actions:
            raise PopulationError("Selected experimental action is not declared")
        visiting.add(identifier)
        for parent in question_set.actions[identifier]["dependency_action_ids"]:
            visit(parent)
        visiting.remove(identifier)
        result.add(identifier)
    for identifier in action_ids:
        visit(identifier)
    return frozenset(result)


def _performs_role(action, role):
    if role not in action.get("performed_control_roles", []):
        return False
    if role == "mock":
        return action["operation_type"] == "mock_control"
    if role == "intended_activity":
        return action["operation_type"] in {"measurement_panel", "intended_assay", "locus_assay"} and action["design_id"]["value"] is not None
    if role == "assay_qc":
        return action["operation_type"] in {"measurement_panel", "intended_assay", "locus_assay", "functional_assay", "mock_control", "assay_control"}
    return action["operation_type"] == "assay_control"


def _control_matches(index, action, control, requirement):
    if any(action[key] != control[key] for key in ("donor_id", "specimen_id", "context_id", "batch_id", "experiment_id", "replicate_position")):
        return False
    if control["design_id"]["value"] not in {None, action["design_id"]["value"]}:
        return False
    if control["assay_capability_id"] != requirement["assay_capability_id"]:
        return False
    if not _performs_role(control, requirement["role"]) or set(requirement["required_qc"]) - set(control.get("performed_qc", [])):
        return False
    requested = {(endpoint, locus) for endpoint in requirement["endpoint_ids"] for locus in requirement["locus_ids"]}
    return requested <= {(pair["endpoint_id"], pair["locus_id"]) for pair in _pairs(control)}


def _action_eligibility(question_set, action):
    index, reasons = question_set.index, []
    from .comparison_domain import target_reasons
    evidence_policy = question_set.study["config"]["policies"]["evidence"]
    for endpoint_id in action["endpoint_ids"] or [None]:
        reasons += target_reasons(evidence_policy, index, design_id=action["design_id"]["value"],
            context_id=action["context_id"], endpoint_id=endpoint_id)
    if action["eligibility_status"] != "eligible":
        reasons.extend(action["reason_codes"] or ["declared_action_not_eligible"])
    specimen, assay = index[action["specimen_id"]], index[action["assay_capability_id"]]
    if specimen["availability"] not in ({"verified"} if question_set.mode == "research_sample_planning" else {"verified", "simulation_assumption"}):
        reasons.append("specimen_availability_unresolved_or_unavailable")
    if specimen.get("use_authorized") is not True or specimen["permitted_use"]["value"] is None:
        reasons.append("specimen_use_authorization_unresolved")
    if action["context_id"] not in specimen["context_ids"]:
        reasons.append("specimen_context_incompatible")
    if assay["feasibility"] not in ({"verified"} if question_set.mode == "research_sample_planning" else {"verified", "simulation_assumption"}):
        reasons.append("assay_feasibility_not_established")
    if question_set.mode == "research_sample_planning" and question_set.coverage_basis != "callable_supported":
        reasons.append("research_mode_requires_callable_source_support")
    measuring = action["operation_type"] in {"measurement_panel", "locus_assay", "intended_assay", "functional_assay", "mock_control", "assay_control"}
    for endpoint_id in action["endpoint_ids"]:
        endpoint = index[endpoint_id]
        if endpoint["locus_resolution"] == "copy_specific" and assay["resolution"] != "copy_specific":
            reasons.append("assay_does_not_resolve_required_genomic_copy")
        if endpoint["allele_resolution"] == "required" and not assay["allele_resolved"]:
            reasons.append("assay_does_not_resolve_required_allele")
    if measuring:
        if action["operation_type"] != "mock_control":
            declared_roles = {requirement["role"] for requirement in action["control_requirements"]}
            if set(assay["required_control_roles"]) - declared_roles:
                reasons.append("mandatory_assay_control_requirement_omitted")
        supplied = {row["locus_id"]: set(row["sequence_configuration_ids"]) for row in action.get("target_sequences", [])}
        if action["target_sequence_id"]["value"] is not None:
            sequence = index[action["target_sequence_id"]["value"]]
            supplied.setdefault(sequence["locus_id"], set()).add(sequence["record_id"])
        for locus_id in action["locus_ids"]:
            copies, unresolved = _actual_copies(index, action["donor_id"], locus_id, question_set.coverage_basis)
            reasons.extend("target_" + reason + ":" + locus_id for reason in unresolved)
            expected = {copy["sequence"]["record_id"] for copy in copies.values()}
            if supplied.get(locus_id) != expected:
                reasons.append("action_target_is_not_actual_donor_sequence:" + locus_id)
        if set(assay["qc_requirements"]) - set(action.get("performed_qc", [])):
            reasons.append("mandatory_assay_qc_operation_unpriced_or_undeclared")
        closure = dependency_closure(question_set, [action["record_id"]])
        for requirement in action["control_requirements"]:
            if not any(_control_matches(index, action, index[identifier], requirement) for identifier in closure):
                reasons.append("compatible_control_prerequisite_missing:" + requirement["role"])
        if action["design_id"]["value"] is not None and action["operation_type"] != "mock_control":
            if not any(index[identifier]["operation_type"] == "treatment" and all(index[identifier][key] == action[key] for key in ("donor_id", "specimen_id", "design_id", "context_id", "batch_id", "experiment_id")) for identifier in closure):
                reasons.append("matching_treatment_preparation_missing")
        resource_types = {index[resource_id]["resource_type"] for identifier in closure
            for demand in index[identifier]["resource_demands"]
            for resource_id in ([demand["resource_id"]] if "resource_id" in demand else demand["eligible_resource_ids"])}
        required = {"study_setup", "donor_activation", "specimen_preparation", "sequencing_batch"}
        required.add("mock_control" if action["operation_type"] == "mock_control" else "assay_library")
        if action["design_id"]["value"] is not None and action["operation_type"] != "mock_control":
            required.add("treatment")
        if action["replicate_position"] > 1:
            required.add("additional_replicate")
        if required - resource_types:
            reasons.append("required_cost_packages_undeclared:" + ",".join(sorted(required - resource_types)))
    return {"eligible": not reasons, "reason_codes": sorted(set(reasons))}


def _replication_reasons(question_set, question, action=None, observation=None):
    if question["template"] != "evidence_replication":
        return []
    index, policy = question_set.index, question["evidence_reuse"]
    originals = [index[identifier] for identifier in question["motivating_observation_ids"]]
    revisions = [row for row in index.values() if row["record_type"] == "comparison_observation_revision" and row["source_record_id"] in question["motivating_observation_ids"]]
    reasons = []
    experiment = action["experiment_id"] if action else index[observation["source_record_id"]]["source_experiment"]
    if experiment in {row["source_experiment"] for row in originals} | set(policy["excluded_experiment_ids"]):
        reasons.append("replication_reuses_original_experiment")
    for dimension in policy["require_new"]:
        if dimension == "experiment":
            continue
        old_values = []
        for revision in revisions:
            if dimension in {"specimen", "donor"}:
                old_values.append(revision["subject"][dimension + "_id"]["value"])
            else:
                executed = revision["executed_action_id"]["value"]
                old_values.append(index[executed]["batch_id"] if executed else None)
        if action:
            actual = action[dimension + "_id"]
        elif dimension in {"specimen", "donor"}:
            actual = observation["subject"][dimension + "_id"]["value"]
        else:
            executed = observation["executed_action_id"]["value"]
            actual = index[executed]["batch_id"] if executed else None
        if not old_values or any(value is None for value in old_values) or actual is None:
            reasons.append("replication_" + dimension + "_independence_unresolved")
        elif actual in old_values:
            reasons.append("replication_reuses_original_" + dimension)
    return reasons


def _action_slot_compatible(question_set, action, slot):
    if not question_set.eligibility[action["record_id"]]["eligible"]:
        return False
    if (action["donor_id"], action["design_id"]["value"], action["context_id"], action["replicate_position"]) != (slot["donor_id"], slot["design_id"]["value"], slot["context_id"], slot["replicate_position"]):
        return False
    if {"endpoint_id": slot["endpoint_id"], "locus_id": slot["locus_id"]["value"]} not in _pairs(action):
        return False
    role = slot["control_role"]
    if role == "measurement":
        if action["operation_type"] not in {"measurement_panel", "locus_assay", "intended_assay", "functional_assay"}:
            return False
    elif not _performs_role(action, role):
        return False
    question = question_set.questions[slot["question_id"]]
    if role != "measurement" and not any(requirement["role"] == role and
            requirement["assay_capability_id"] == action["assay_capability_id"] and
            slot["endpoint_id"] in requirement["endpoint_ids"] and
            slot["locus_id"]["value"] in requirement["locus_ids"] and
            set(requirement["required_qc"]) <= set(action.get("performed_qc", []))
            for requirement in question["control_requirements"]):
        return False
    return not _replication_reasons(question_set, question, action=action)


def _observation_slot_reasons(question_set, observation, slot):
    question, index = question_set.questions[slot["question_id"]], question_set.index
    policy, reasons = question["evidence_reuse"], outcome_reasons(observation)
    from .comparison_domain import biological_reasons
    reasons += biological_reasons(question_set.study["config"]["policies"]["evidence"], index,
        design_id=slot["design_id"]["value"], context_id=slot["context_id"], endpoint_id=slot["endpoint_id"])
    if not policy["allow_existing"] or question_set.coverage_basis == "retained_callset_model":
        reasons.append("existing_empirical_evidence_not_allowed")
    excluded = set(policy["excluded_observation_ids"]) | set(question["motivating_observation_ids"] if question["template"] == "evidence_replication" else [])
    source = index[observation["source_record_id"]]
    if observation["source_record_id"] in excluded or source["source_experiment"] in policy["excluded_experiment_ids"]:
        reasons.append("excluded_original_measurement_or_experiment")
    token = question_set.observation_tokens[observation["record_id"]]
    if any(question_set.observation_tokens.get(other["record_id"]) == token for other in index.values()
           if other["record_type"] == "comparison_observation_revision" and other["source_record_id"] in excluded):
        reasons.append("equivalent_export_of_excluded_measurement")
    if not observation["applicability"]["empirical_donor_credit"] or observation["control_status"] not in {"passed", "not_applicable"} or observation["qc_status"] != "passed" or observation["normalized_value"]["value"] is None:
        reasons.append("measurement_quality_or_exact_donor_applicability_failed")
    if (observation["subject"]["donor_id"]["value"], observation["design_id"]["value"], observation["context_id"]["value"], observation["endpoint_id"]) != (slot["donor_id"], slot["design_id"]["value"], slot["context_id"], slot["endpoint_id"]):
        reasons.append("observation_donor_design_context_endpoint_mismatch")
    locus = slot["locus_id"]["value"]
    if locus is not None and locus not in observation["locus_ids"]:
        reasons.append("observation_locus_mismatch")
    if locus is not None:
        copies, unknowns = _actual_copies(index, slot["donor_id"], locus, question_set.coverage_basis)
        if unknowns or observation["sequence_configuration_id"]["value"] not in {copy["sequence"]["record_id"] for copy in copies.values()}:
            reasons.append("observation_actual_donor_sequence_mismatch")
    if slot["control_role"] != "measurement":
        executed = observation["executed_action_id"]["value"]
        if executed is None or not _performs_role(index[executed], slot["control_role"]):
            reasons.append("completed_control_operation_not_established")
        elif (index[executed]["donor_id"], index[executed]["design_id"]["value"], index[executed]["context_id"]) != (slot["donor_id"], slot["design_id"]["value"], slot["context_id"]) or {"endpoint_id":slot["endpoint_id"],"locus_id":locus} not in _pairs(index[executed]):
            reasons.append("completed_control_operation_is_not_for_this_requirement")
        elif index[executed]["experiment_id"] != source["source_experiment"]:
            reasons.append("completed_control_experiment_does_not_match_source")
        elif not any(requirement["role"] == slot["control_role"] and requirement["assay_capability_id"] == index[executed]["assay_capability_id"] and slot["endpoint_id"] in requirement["endpoint_ids"] and locus in requirement["locus_ids"] and set(requirement["required_qc"]) <= set(index[executed].get("performed_qc", [])) for requirement in question["control_requirements"]):
            reasons.append("completed_control_does_not_meet_declared_assay_qc_requirement")
    reasons.extend(_replication_reasons(question_set, question, observation=observation))
    return sorted(set(reasons))


def prepare_question_set(study, records):
    """Validate and independently recompute the choices offered to the planner."""
    verify_sources(study)
    merged = {}
    for row in [*study["config"]["catalogs"], *records]:
        old = merged.get(row["record_id"])
        if old is not None and old != row:
            raise PopulationError("Question inputs contain conflicting canonical records")
        merged[row["record_id"]] = deepcopy(row)
    index = _QuestionIndex(list(validate_catalog(list(merged.values())).values()))
    qset = QuestionSet(study, list(index.values()), index, *(_kind(index, name) for name in ("question", "requirement_slot", "action", "resource")))
    for identifier, action in qset.actions.items():
        qset.eligibility[identifier] = _action_eligibility(qset, action)
    for identifier in qset.actions:
        unavailable = [parent for parent in dependency_closure(qset, [identifier]) if not qset.eligibility[parent]["eligible"]]
        if unavailable:
            qset.eligibility[identifier] = {"eligible": False, "reason_codes": sorted(set(qset.eligibility[identifier]["reason_codes"] + ["ineligible_dependency:" + parent for parent in unavailable if parent != identifier]))}
    for question in qset.questions.values():
        if question.get("construction_status") == "excluded_by_declared_limit":
            continue
        for donor_id in _kind(index, "donor"):
            membership = donor_membership(index, question, donor_id, qset.mode, qset.coverage_basis)
            qset.donor_eligibility[(question["record_id"], donor_id)] = membership
    equivalents = {}
    for row in _kind(index, "observation_equivalence").values():
        for identifier in row["source_record_ids"]:
            equivalents[identifier] = row["record_id"]
    for observation in _kind(index, "observation_revision").values():
        source = index[observation["source_record_id"]]
        qset.observation_tokens[observation["record_id"]] = equivalents.get(source["record_id"]) or content_id("comparison_actual_measurement", {
            "study": source["source_study"], "experiment": source["source_experiment"], "replicate": source["source_replicate"],
            "donor": observation["subject"]["donor_id"]["value"], "specimen": observation["subject"]["specimen_id"]["value"],
            "design": observation["design_id"]["value"], "context": observation["context_id"]["value"],
            "endpoint": observation["endpoint_id"], "loci": sorted(observation["locus_ids"]), "sequence": observation["sequence_configuration_id"]["value"]})
    active = {row["record_id"] for row in active_observations(list(index.values()))}
    for identifier, slot in qset.slots.items():
        membership = qset.donor_eligibility[(slot["question_id"], slot["donor_id"])]
        arm = next(arm for arm in qset.questions[slot["question_id"]]["arms"] if arm["arm_id"] == slot["arm_id"])
        eligible = membership["eligible"] and (arm["membership"] == membership["membership"] or arm["membership"] == "all_resolved")
        actions = [action_id for action_id in slot["allowed_action_ids"] if eligible and _action_slot_compatible(qset, qset.actions[action_id], slot)]
        observations = []
        for observation_id in slot["allowed_observation_ids"]:
            reasons = _observation_slot_reasons(qset, index[observation_id], slot)
            if observation_id not in active:
                reasons.append("superseded_observation_revision")
            if not reasons and eligible:
                observations.append(observation_id)
            elif reasons:
                qset.evidence_rejections[(identifier, observation_id)] = sorted(set(reasons))
        qset.options[identifier] = {"action_ids": sorted(actions), "observation_ids": sorted(observations)}
    return qset


def slot_options(question_set, slot_id):
    return deepcopy(question_set.options[slot_id])


def eligible_donor_bundles(question_set, question_id):
    question = question_set.questions[question_id]
    groups = defaultdict(list)
    for identifier in question["requirement_slot_ids"]:
        slot = question_set.slots[identifier]
        groups[(slot["donor_id"], slot["arm_id"])].append(identifier)
    result = []
    for (donor_id, arm_id), slot_ids in sorted(groups.items()):
        membership = question_set.donor_eligibility[(question_id, donor_id)]
        arm = next(arm for arm in question["arms"] if arm["arm_id"] == arm_id)
        reasons = list(membership["reason_codes"])
        reasons.extend(_required_slot_gaps(question_set, question, slot_ids))
        if arm["membership"] != "all_resolved" and arm["membership"] != membership["membership"]:
            reasons.append("donor_does_not_meet_declared_arm_predicate")
        if any(not any(question_set.options[identifier].values()) for identifier in slot_ids):
            reasons.append("required_slot_has_no_eligible_measurement_or_action")
        result.append({"donor_id": donor_id, "arm_id": arm_id, "independence_group_id": membership["independence_group_id"],
                       "eligible": not reasons, "reason_codes": sorted(set(reasons)), "assumptions": membership["assumptions"],
                       "slots": [{"slot_id": identifier, **slot_options(question_set, identifier)} for identifier in sorted(slot_ids)]})
    return result


def _required_slot_gaps(question_set, question, slot_ids):
    """Rebuild the whole scientific panel independently of supplied slot lists."""
    actual = {(slot["design_id"]["value"], slot["endpoint_id"], slot["locus_id"]["value"],
               slot["context_id"], slot["replicate_position"], slot["control_role"])
              for slot in (question_set.slots[identifier] for identifier in slot_ids)}
    reasons = []
    for design in question["design_ids"]:
        measurements = required_design_measurements(question_set.index, question, design)
        own_intended = set(question_set.index[design]["intended_locus_ids"])
        if own_intended - {pair["locus_id"] for pair in measurements if pair["control_role"] == "intended_activity"}:
            reasons.append("question_omits_required_intended_locus")
        for pair in measurements:
            endpoint, locus = pair["endpoint_id"], pair["locus_id"]
            role = pair["control_role"]
            for context in question["context_ids"]:
                for repeat in range(1, question["technical_replicates"] + 1):
                    for required_role in (role, "assay_qc"):
                        if (design, endpoint, locus, context, repeat, required_role) not in actual:
                            reasons.append("required_paired_or_control_slot_omitted:" + required_role)
                    if not any((mock_design, endpoint, locus, context, repeat, "mock") in actual for mock_design in (None, design)):
                        reasons.append("required_paired_or_control_slot_omitted:mock")
    return sorted(set(reasons))


def _selected_donors(question_set, question, donors):
    arms = sorted(question["arms"], key=lambda arm: arm["arm_id"])
    choices = []
    for arm in arms:
        complete = [donor for donor in donors if donor["arm_id"] == arm["arm_id"] and donor["complete"]]
        required = arm["required_independent_donors"]
        if len(complete) < required:
            return None
        choices.append((sorted(complete, key=lambda donor: donor["donor_id"]), required))
    visited = 0
    def visit(position, selected, groups, blocks):
        nonlocal visited
        if position == len(choices):
            return sorted(donor["donor_id"] for donor in selected)
        for candidate in combinations(*choices[position]):
            visited += 1
            if visited > question_set.study["config"]["policies"]["planning"]["max_expansions"]:
                raise ResourceLimit("Independent donor matching exceeded declared max_expansions")
            local = [donor["independence_group_id"] for donor in candidate]
            if None in local or len(set(local)) != len(local) or groups & set(local):
                continue
            block = Counter(canonical_bytes(question_set.donor_eligibility[(question["record_id"], donor["donor_id"])]["covariates"]) for donor in candidate)
            if question["matching_covariates"] and blocks and block != blocks[0]:
                continue
            result = visit(position + 1, [*selected,*candidate], groups | set(local), [*blocks,block])
            if result is not None:
                return result
        return None
    return visit(0, [], set(), [])


def _requirement_state(question_set, assignments):
    questions, completed = {}, []
    for identifier, question in sorted(question_set.questions.items()):
        donors = []
        for bundle in eligible_donor_bundles(question_set, identifier):
            slot_ids = [slot["slot_id"] for slot in bundle["slots"]]
            unmet = sorted(set(slot_ids) - assignments.keys())
            donors.append({key: bundle[key] for key in ("donor_id", "arm_id", "independence_group_id", "eligible", "reason_codes")} |
                          {"complete": bundle["eligible"] and not unmet, "slot_ids": slot_ids, "unmet_slot_ids": unmet})
        selected = _selected_donors(question_set, question, donors)
        complete = selected is not None
        own_slots = set(question["requirement_slot_ids"])
        observed = {slot for slot in own_slots & assignments.keys() if "observation_id" in assignments[slot]}
        measured_donors = [{**donor, "complete": donor["eligible"] and set(donor["slot_ids"]) <= observed} for donor in donors]
        observed_complete = _selected_donors(question_set, question, measured_donors) is not None
        evidence = "evaluable" if observed_complete else "partial_measurements" if observed else "no_measurements"
        if not observed and any(slot in own_slots for slot, _ in question_set.evidence_rejections):
            evidence = "insufficient_quality"
        planning = "already_satisfied" if observed_complete else "fully_planned" if complete else "partially_planned" if own_slots & assignments.keys() else "blocked" if not any(donor["eligible"] for donor in donors) else "unplanned"
        gaps = [] if complete else list(question["remaining_gaps"]) if not question["requirement_slot_ids"] else ["unmet_independent_donor_or_paired_slot_requirements"]
        questions[identifier] = {"complete": complete, "evidence_status": evidence, "planning_status": planning,
            "donors": donors, "selected_donor_ids": selected or [], "remaining_gaps": gaps}
        if complete:
            completed.append(identifier)
    return {"questions": questions, "completed_question_ids": completed,
            "satisfied_slot_ids": sorted(assignments), "slot_assignments": deepcopy(assignments),
            "weighted_completion": sum(question_set.questions[identifier]["priority_weight"] for identifier in completed),
            "unweighted_completion": len(completed)}


def evaluate_requirements(question_set, selected_action_ids=(), *, observation_assignments=None):
    """Find a feasible joint assignment; repeats never reuse one actual measurement.

    Ambiguous assignments are checked exactly under the declared expansion bound.
    Exceeding that bound raises a resource error, never a false completeness claim.
    """
    selected = set(selected_action_ids)
    if dependency_closure(question_set, selected) != selected:
        raise PopulationError("Selected actions omit a required prerequisite")
    if any(not question_set.eligibility[identifier]["eligible"] for identifier in selected):
        raise PopulationError("Selected action fails source-bound experimental eligibility")
    options = {}
    for identifier, slot in question_set.slots.items():
        candidates = [{"observation_id": observation} for observation in question_set.options[identifier]["observation_ids"]]
        candidates += [{"action_id": action} for action in question_set.options[identifier]["action_ids"] if action in selected]
        if observation_assignments is not None and identifier in observation_assignments:
            requested = observation_assignments[identifier]
            requested = {"observation_id": requested} if isinstance(requested, str) else requested
            if requested not in candidates:
                raise PopulationError("Explicit slot assignment is not source-compatible")
            candidates = [requested]
        options[identifier] = candidates
    def token(slot_id, assignment):
        if "observation_id" in assignment:
            return "observation", question_set.observation_tokens[assignment["observation_id"]]
        return "action", assignment["action_id"], *_biological_slot(question_set.slots[slot_id])
    def share(left, right):
        a, b = question_set.slots[left], question_set.slots[right]
        if a["question_id"] == b["question_id"]:
            # A passed QC annotation belongs to its actual measured endpoint.
            # It does not create a second technical replicate or mock result.
            same_measurement = _biological_slot(a)[:-1] == _biological_slot(b)[:-1] and a["replicate_position"] == b["replicate_position"]
            return left == right or same_measurement and "assay_qc" in {a["control_role"], b["control_role"]}
        return (_biological_slot(a) == _biological_slot(b)
            and question_set.questions[a["question_id"]]["evidence_reuse"]["allow_across_questions"]
            and question_set.questions[b["question_id"]]["evidence_reuse"]["allow_across_questions"])
    users = defaultdict(set)
    for slot_id, candidates in options.items():
        for candidate in candidates:
            users[token(slot_id, candidate)].add(slot_id)
    assignments, contested = {}, []
    for slot_id, candidates in sorted(options.items()):
        if not candidates:
            continue
        conflicts = any(not share(slot_id, other) for candidate in candidates for other in users[token(slot_id, candidate)])
        if conflicts:
            contested.append(slot_id)
        else:
            assignments[slot_id] = candidates[0]
    if not contested:
        return _requirement_state(question_set, assignments)
    limit = question_set.study["config"]["policies"]["planning"]["max_expansions"]
    visited, best, best_score = 0, None, None
    occupied = defaultdict(list)
    for slot_id, assignment in assignments.items():
        occupied[token(slot_id, assignment)].append(slot_id)
    def visit(position):
        nonlocal visited, best, best_score
        visited += 1
        if visited > limit:
            raise ResourceLimit("Exact evidence/repeat assignment exceeded declared max_expansions")
        if position == len(contested):
            result = _requirement_state(question_set, assignments)
            hard = {identifier for identifier, question in question_set.questions.items() if question["hard_required"]}
            completed_hard = hard & set(result["completed_question_ids"])
            score = (completed_hard == hard, len(completed_hard),
                     sum(question_set.questions[identifier]["priority_weight"] for identifier in completed_hard),
                     result["weighted_completion"], result["unweighted_completion"], len(assignments))
            if best_score is None or score > best_score:
                best, best_score = result, score
            return
        slot_id = contested[position]
        for assignment in options[slot_id]:
            key = token(slot_id, assignment)
            if all(share(slot_id, other) for other in occupied[key]):
                assignments[slot_id] = assignment
                occupied[key].append(slot_id)
                visit(position + 1)
                occupied[key].pop()
                del assignments[slot_id]
        if observation_assignments is None or slot_id not in observation_assignments:
            visit(position + 1)
    visit(0)
    if best is None:
        raise PopulationError("Explicit observation assignments reuse an incompatible measurement")
    return best


def generate_questions(study, records):
    """Generate canonical questions and operations from one pinned package policy."""
    from .comparison_question_generation import generate_questions as generate
    return generate(study, records)
