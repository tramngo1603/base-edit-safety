"""Render verified comparative records without recomputing scientific decisions."""
from __future__ import annotations

from collections import Counter, defaultdict
from copy import deepcopy
import csv
import io
import json
import re

from .common import PopulationError
from .comparison_contracts import active_observations, validate_catalog, validate_record
from .comparison_evidence import _quality
from .comparison_domain import assess_domain, biological_reasons
from .comparison_questions import dependency_closure
from .comparison_provenance import verify_sources
from .experiment_planning import _verify_plan_identity


REPORT_FORMAT_VERSION = "2"
OBSERVATION_OUTCOME_CONTRACT = "comparison_observation_revision:v2"
ACCEPTED_OBSERVATION_HISTORY_VERSIONS = ["1", "2"]
HUMAN_RECORD_FIELDS = (
    ("reviewer_attestation", "reviewer attestation"),
    ("funding_disclosure", "funding disclosure"),
    ("customer_disclosure", "customer disclosure"),
)
NOT_RECORDED_REASON = "No {label} is present in the report inputs"


def _report_format(study):
    """Return the export contract, captured code identity and missing review-record statuses."""
    return {"format_version": REPORT_FORMAT_VERSION,
            "observation_outcome_contract": OBSERVATION_OUTCOME_CONTRACT,
            "accepted_observation_history_versions": list(ACCEPTED_OBSERVATION_HISTORY_VERSIONS),
            "producing_code": deepcopy(study["code"]),
            "human_records": {key: {"status": "not_recorded", "reason": NOT_RECORDED_REASON.format(label=label)}
                              for key, label in HUMAN_RECORD_FIELDS}}


def _records(artifacts):
    return [row for name, rows in artifacts.items() if name.endswith(".jsonl") for row in rows]


def _kind(records, name):
    return sorted((row for row in records if row["record_type"] == "comparison_" + name),
                  key=lambda row: row["record_id"])


def _value(value):
    return value.get("value") if isinstance(value, dict) and "value" in value else value


def _text(value):
    if isinstance(value, dict) and "value" in value:
        return "unknown (" + str(value["reason"]) + ")" if value["value"] is None else _text(value["value"])
    if value is None:
        return "unknown"
    if isinstance(value, bool):
        return "yes" if value else "no"
    return str(value)


def _cell(value):
    if value is None:
        return "null"
    if isinstance(value, (dict, list, bool)):
        return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False)
    return str(value)


def _tsv(columns, rows):
    stream = io.StringIO(newline="")
    writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
    writer.writerow(columns)
    for row in rows:
        writer.writerow([_cell(row.get(name)) for name in columns])
    return stream.getvalue()


def _table(headers, rows):
    def clean(value):
        return _text(value).replace("|", "\\|").replace("\r", " ").replace("\n", " ")
    return "\n".join(["| " + " | ".join(headers) + " |", "| " + " | ".join("---" for _ in headers) + " |",
                      *["| " + " | ".join(clean(value) for value in row) + " |" for row in rows]]) + "\n"


def _unique_summaries(records, kind, keys):
    result = {}
    for row in _kind(records, kind):
        key = tuple(row[field] for field in keys)
        if key in result:
            raise PopulationError("Reporting received multiple current summaries for one comparison unit")
        result[key] = row
    return result


class _Readable:
    def __init__(self, study, index):
        self.study, self.index = study, index
        self.policy = study["config"]["policies"]["report"]
        donors = _kind(index.values(), "donor")
        self.donors = {row["record_id"]: row["canonical_label"] if self.policy["include_donor_identifiers"]
                       else "Donor " + str(i + 1) for i, row in enumerate(donors)}
        self.questions = {row["record_id"]: "Q" + str(i + 1)
                          for i, row in enumerate(_kind(index.values(), "question"))}
        self.actions = {row["record_id"]: "Operation " + str(i + 1)
                        for i, row in enumerate(_kind(index.values(), "action"))}
        self.replacements = {}
        if not self.policy["include_donor_identifiers"]:
            for row in index.values():
                kind = row["record_type"]
                if kind == "comparison_donor":
                    for value in (row["record_id"], row["canonical_label"], *row["sample_ids"]):
                        self.replacements[value] = self.donors[row["record_id"]]
                    if _value(row["independence_group_id"]):
                        self.replacements[_value(row["independence_group_id"])] = "[relationship group withheld]"
                elif kind == "comparison_specimen":
                    for key in ("record_id", "source_specimen_id"):
                        self.replacements[row[key]] = "[specimen identifier withheld]"
                elif kind == "comparison_observation_revision":
                    for value in (row["subject"]["source_subject_id"], row["independent_unit_id"],
                                  row["relationship_group"], row["pairing_key"]):
                        if _value(value):
                            self.replacements[_value(value)] = "[source subject identifier withheld]"
        if not self.policy["include_local_paths"]:
            self.replacements[str(study["config_path"])] = "[local path withheld]"
            for source in study["sources"].values():
                self.replacements[source["path"]] = "[local path withheld]"
            for source in study["config"]["input_refs"]:
                self.replacements[source["path"]] = "[local path withheld]"
        self.pattern = re.compile(r"(?<!\w)(?:" + "|".join(re.escape(value) for value in
            sorted(self.replacements, key=lambda value: (-len(value), value))) + r")(?!\w)") if self.replacements else None

    def label(self, identifier):
        if identifier is None:
            return "design-neutral control"
        if identifier in self.donors:
            return self.donors[identifier]
        if identifier in self.questions:
            return self.questions[identifier]
        if identifier in self.actions:
            return self.actions[identifier]
        row = self.index.get(identifier, {})
        return self.redact(next((row[key] for key in ("display_name", "locus_name", "name") if key in row), identifier))

    def context(self, identifier):
        if identifier is None:
            return "context unknown"
        row = self.index[identifier]
        timepoint = row["measurement_timepoint"]
        return "; ".join([self.redact(row["cell_type"]), self.redact(row["delivery_method"]),
                          self.value(timepoint["value"]) + " " + self.redact(timepoint["unit"]), self.redact(row["condition"])])

    def value(self, value):
        """Preserve numbers even when a source happens to use a numeric donor ID."""
        if isinstance(value, dict) and "value" in value:
            return "unknown (" + self.redact(value["reason"]) + ")" if value["value"] is None else self.value(value["value"])
        if isinstance(value, str) and not re.fullmatch(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", value):
            return self.redact(value)
        return _text(value)

    def explain(self, text):
        text = self.redact(text)
        return re.sub(r"comparison_(?:design|locus|intended_locus|endpoint_definition):v1:[0-9a-f]{64}",
                      lambda match: self.label(match.group(0)), text)

    def redact(self, text):
        text = str(text)
        if not self.policy["include_local_paths"]:
            text = re.sub(r"file://[^\s|<>`]+", "[local path withheld]", text)
            text = re.sub(r"(?<![\w:/])/(?:[^\s|<>`]+)", "[local path withheld]", text)
            text = re.sub(r"[A-Za-z]:[\\/][^\s|<>`]+", "[local path withheld]", text)
        if self.pattern is not None:
            text = self.pattern.sub(lambda match: self.replacements[match.group(0)], text)
        if not self.policy["include_donor_identifiers"]:
            text = re.sub(r"comparison_donor:v1:[0-9a-f]{64}", "[donor identifier withheld]", text)
        return text


def _experimental_rows(question_set, plan, index):
    questions_by_action, roles = defaultdict(set), defaultdict(set)
    for slot_id, assignment in plan["requirements"]["slot_assignments"].items():
        if "action_id" not in assignment:
            continue
        slot = question_set.slots[slot_id]
        roles[assignment["action_id"]].add(slot["control_role"])
        for identifier in dependency_closure(question_set, [assignment["action_id"]]):
            questions_by_action[identifier].add(slot["question_id"])
    alternatives = {(row["action_id"], row["assignment_key"]): row["resource_id"]
                    for row in plan["resource_assignments"]}
    rows, resource_users = [], defaultdict(set)
    for identifier in plan["selected_action_ids"]:
        action = question_set.actions[identifier]
        demands = []
        for demand in action["resource_demands"]:
            rid = demand["resource_id"] if "resource_id" in demand else alternatives[(identifier, demand["assignment_key"])]
            demands.append({"resource_id": rid, "dimension": demand["dimension"], "quantity": demand["quantity"],
                            "assignment_key": demand.get("assignment_key")})
            resource_users[rid].add(identifier)
        targets = action.get("target_sequences", [])
        if not targets and _value(action["target_sequence_id"]) is not None:
            seq = index[_value(action["target_sequence_id"])]
            targets = [{"locus_id": seq["locus_id"], "sequence_configuration_ids": [seq["record_id"]]}]
        exact = [{"locus_id": target["locus_id"], "sequences": [
            {key: index[sid][key] for key in ("record_id", "contig", "strand", "sequence", "path", "variants", "support_status", "support_reason")}
            for sid in sorted(target["sequence_configuration_ids"])]} for target in targets]
        qids = sorted(questions_by_action[identifier])
        rows.append({"action_id": identifier, "operation_type": action["operation_type"], "donor_id": action["donor_id"],
            "specimen_id": action["specimen_id"], "design_id": _value(action["design_id"]),
            "context_id": action["context_id"], "batch_id": action["batch_id"], "experiment_id": action["experiment_id"],
            "replicate_position": action["replicate_position"], "assay_capability_id": action["assay_capability_id"],
            "endpoint_locus_pairs": action.get("endpoint_locus_pairs", []), "target_sequences": exact,
            "control_requirements": action["control_requirements"], "performed_control_roles": action.get("performed_control_roles", []),
            "performed_qc": action.get("performed_qc", []), "satisfied_roles": sorted(roles[identifier]),
            "dependency_action_ids": action["dependency_action_ids"], "resource_assignments": demands,
            "question_ids": qids, "rationales": [question_set.questions[q]["rationale"] for q in qids],
            "arms": [{"question_id": qid, "arm_id": donor["arm_id"]} for qid in qids
                     for donor in plan["questions"][qid]["donors"] if donor["donor_id"] == action["donor_id"]]})
    return rows, resource_users


def _unmet_rows(question_set, plan):
    rows = []
    for qid, state in sorted(plan["questions"].items()):
        for reason in state["remaining_gaps"]:
            rows.append({"kind": "question", "question_id": qid, "reason": reason,
                         "planning_status": state["planning_status"], "question_complete": state["complete"]})
        if qid in plan["unmet_hard_required_question_ids"]:
            rows.append({"kind": "mandatory_question", "question_id": qid, "reason": "mandatory_question_not_completed",
                         "planning_status": state["planning_status"], "question_complete": False,
                         "detail": {"solver_status": plan["status"], "budget": plan["budget"],
                                    "certificate": plan["infeasibility_certificate"]}})
        for donor in state["donors"]:
            for reason in donor["reason_codes"]:
                rows.append({"kind": "donor_option", "question_id": qid, "donor_id": donor["donor_id"],
                             "arm_id": donor["arm_id"], "reason": reason, "question_complete": state["complete"]})
            if not state["complete"]:
                for sid in donor["unmet_slot_ids"]:
                    rows.append({"kind": "candidate_slot", "question_id": qid, "donor_id": donor["donor_id"],
                                 "arm_id": donor["arm_id"], "slot_id": sid, "reason": "slot_in_donor_option_not_filled",
                                 "question_complete": False, "detail": question_set.slots[sid]})
    for aid, eligibility in sorted(question_set.eligibility.items()):
        for reason in eligibility["reason_codes"]:
            rows.append({"kind": "action_option", "action_id": aid, "reason": reason})
    rows.extend({"kind": "accounting", "reason": problem["code"], "detail": problem} for problem in plan["accounting"]["problems"])
    return rows


def render_report(study, *, search, observations, evidence, comparison, questions, plan, question_set):
    """Return flat report files after the caller verifies the plan revision.

    The caller must load verified stage ancestors and call verify_plan_revision.
    This renderer checks record/audit identity, then presents their recorded
    decisions, eligibility, assignments and costs without rerunning selection.
    The stage writer supplies run-manifest.json and protected output permissions.
    """
    verify_sources(study)
    if question_set.study["config"] != study["config"]:
        raise PopulationError("Report and question set use different study configurations")
    validate_record(plan["plan_record"])
    _verify_plan_identity(plan)
    if set(plan["questions"]) != set(question_set.questions):
        raise PopulationError("Report plan does not describe the supplied questions")
    merged = dict(question_set.index)
    for artifacts in (search, observations, evidence, comparison, questions):
        for row in _records(artifacts):
            if row["record_id"] in merged and merged[row["record_id"]] != row:
                raise PopulationError("Report ancestors contain conflicting canonical records")
            merged[row["record_id"]] = row
    index = validate_catalog(list(merged.values()))
    current = comparison.get("comparison-summaries.jsonl", [])
    designs = _unique_summaries(current, "design_summary", ("design_id",))
    donor_designs = _unique_summaries(current, "donor_design_summary", ("donor_id", "design_id"))
    decisions = _kind(current, "decision")
    reader = _Readable(study, index)
    applications = _kind(_records(evidence), "evidence_applicability")
    active = active_observations(list(index.values()))
    if study["config"]["coverage_basis"] == "retained_callset_model" and any(
            row["assessment"]["empirical_donor_credit"] for row in applications):
        raise PopulationError("Retained-callset report contains forbidden empirical donor credit")
    domain = index[study["config"]["policies"]["comparison"]["domain_id"]]
    profiles = {index[identifier]["design_id"]: index[identifier]
                for identifier in study["config"]["policies"]["search"]["profile_ids"]}
    design_rows = []
    count_fields = ("complete_locus_count", "complete_sequence_configuration_count", "complete_alignment_count")
    for did in sorted(study["config"]["policies"]["search"]["design_ids"]):
        summary, design = designs.get((did,), {}), index[did]
        profile = profiles[did]
        row = {**summary, "design_id": did, "design_name": design["display_name"],
               "editor": index[design["editor_definition_id"]]["variant_name"],
               "intended_locus_ids": design["intended_locus_ids"], "common_requested_domain_id": domain["record_id"],
               "search_profile_id": profile["record_id"], "pam_patterns": profile["pam_patterns"],
               "max_mismatches": profile["max_mismatches"], "max_dna_bulges": profile["max_dna_bulges"],
               "max_rna_bulges": profile["max_rna_bulges"],
               "scope_status": "not_available" if not summary else "complete" if summary["search_unresolved_donors"] == 0
               and _value(summary["complete_alignment_count"]) is not None else "partial_or_unresolved",
               "evidence_endpoint_ids": sorted({index[oid]["endpoint_id"] for oid in summary.get("evidence_observation_ids", [])}),
               "search_scopes": [index[sid] for sid in summary.get("search_scope_ids", [])]}
        for key in count_fields:
            value = summary.get(key, {"value": None, "reason": "comparison summary not available"})
            row[key], row[key + "_reason"] = value["value"], value["reason"]
        design_rows.append(row)
    donor_rows = []
    for summary in donor_designs.values():
        donor_rows.append({**summary, "complete_configuration_count": summary["complete_configuration_count"]["value"],
            "complete_configuration_count_reason": summary["complete_configuration_count"]["reason"],
            "intended_assessments": [index[key] for key in summary["intended_assessment_ids"]]})
    evidence_rows = []
    active_ids = {row["record_id"] for row in active}
    evidence_policy = study["config"]["policies"]["evidence"]
    for application in applications:
        assessment = application["assessment"]
        observation = index[application["observation_id"]]
        gate = _quality(observation, evidence_policy)
        current_gate_reasons = list(gate["reason_codes"])
        current_gate_reasons += biological_reasons(evidence_policy, index,
            design_id=application["design_id"], context_id=application["context_id"], endpoint_id=application["endpoint_id"])
        if study["config"]["mode"] != "research_sample_planning":
            current_gate_reasons.append("study_mode_disallows_empirical_donor_credit")
        if not evidence_policy["allow_empirical_donor_credit"]:
            current_gate_reasons.append("study_policy_disallows_empirical_donor_credit")
        current_credit = bool(assessment["empirical_donor_credit"]) and observation["record_id"] in active_ids \
            and gate["status"] == "applicable" and not current_gate_reasons
        evidence_rows.append({**application, "donor_id": _value(application["donor_id"]),
            "classification": assessment["classification"], "empirical_donor_credit": assessment["empirical_donor_credit"],
            "identity": assessment["identity"], "context": assessment["context"], "endpoint": assessment["endpoint"],
            "quality": assessment["quality"], "reason_codes": assessment["reason_codes"],
            "normalized_value": observation["normalized_value"], "normalized_unit": observation["normalized_unit"],
            "source_record_id": observation["source_record_id"],
            "adjudicated_outcome": observation.get("adjudicated_outcome") or "legacy_unadjudicated",
            "adjudication_reason": observation.get("adjudication_reason") or
                "Legacy version-1 revision requires explicit adjudication; conclusive credit withheld",
            "source_reported_outcome": observation["source_reported_outcome"],
            "detection_limit_value": observation["detection_limit"]["value"],
            "detection_limit_unit": observation["detection_limit"]["unit"],
            "normalized_value_reason": observation["normalized_value"]["reason"],
            "observation_active_revision": observation["record_id"] in active_ids,
            "current_empirical_donor_credit": current_credit,
            "current_gate_reason_codes": sorted(set(current_gate_reasons))})
    out_of_scope = [deepcopy(row) for row in evidence.get("published-alignment-replay.json", [])
                    if row["disposition"] != "recovered_exact_alignment"]
    for row in evidence.get("unresolved-observations.json", []):
        out_of_scope.append({**row, "disposition": "unresolved_observation_context", "reason_codes": [row["reason"]],
                             "observation": index[row["observation_id"]]})
    question_rows = []
    for qid, question in sorted(question_set.questions.items()):
        state = plan["questions"][qid]
        question_rows.append({**question, "planning_status": state["planning_status"],
            "evidence_status": state["evidence_status"], "remaining_gaps": state["remaining_gaps"],
            "selected_donor_ids": state["selected_donor_ids"], "donor_options": state["donors"],
            "construction_status": question.get("construction_status", "complete")})
    action_rows, resource_users = _experimental_rows(question_set, plan, index)
    resource_rows = [{**index[row["resource_id"]], **row, "cost_unit": plan["accounting"]["cost_unit"],
                      "cost_scale": plan["accounting"]["cost_scale"], "charged_once": True,
                      "used_by_action_ids": sorted(resource_users[row["resource_id"]]),
                      "actual_demands": [d for d in plan["accounting"]["resource_demands"] if d["resource_id"] == row["resource_id"]]}
                     for row in plan["accounting"]["cost_breakdown"]]
    if len({row["resource_id"] for row in resource_rows}) != len(resource_rows) or sum(row["cost"] for row in resource_rows) != plan["accounting"]["total_cost"]:
        raise PopulationError("Report resource charges disagree with the plan accounting")
    intended = _kind(index.values(), "intended_assessment")
    effects = _kind(index.values(), "variant_effect")
    predictions = _kind(index.values(), "prediction")
    gaps = _unmet_rows(question_set, plan)
    columns = {
        "design-comparison.tsv": ("design_id design_name editor intended_locus_ids common_requested_domain_id search_profile_id pam_patterns max_mismatches max_dna_bulges max_rna_bulges scope_status donor_count intended_evaluable_donors intended_unresolved_donors search_evaluable_donors search_unresolved_donors observed_locus_count observed_sequence_configuration_count observed_alignment_count complete_locus_count complete_locus_count_reason complete_sequence_configuration_count complete_sequence_configuration_count_reason complete_alignment_count complete_alignment_count_reason supported_carrier_donors modeled_carrier_donors evidence_observation_ids evidence_endpoint_ids variant_effect_ids reason_codes search_scopes", design_rows),
        "donor-design-matrix.tsv": ("donor_id design_id domain_id intended_evaluability search_status observed_configuration_count observed_alignment_count complete_configuration_count complete_configuration_count_reason evidence_availability evidence_observation_ids nominated_configuration_ids nominated_alignment_ids intended_assessments search_scope_ids reason_codes", donor_rows),
        "evidence-applicability.tsv": ("record_id observation_id source_record_id design_id context_id endpoint_id donor_id sequence_configuration_id classification empirical_donor_credit identity context endpoint quality reason_codes normalized_value normalized_unit adjudicated_outcome adjudication_reason source_reported_outcome detection_limit_value detection_limit_unit normalized_value_reason observation_active_revision current_empirical_donor_credit current_gate_reason_codes", evidence_rows),
        "out-of-scope-evidence.tsv": ("source_key observation_id source_artifact_id source_locator design_id source_positive source_reported_label disposition reason_codes alignment exact_alignment_ids related_alignment_ids mapping_sha256 empirical_target_donor_credit observation", out_of_scope),
        "questions.tsv": ("question_id template hypothesis decision design_ids endpoint_locus_pairs context_ids arms pairing technical_replicates genotype_predicate relatedness_policy matching_covariates control_requirements evidence_reuse priority_weight hard_required construction_status planning_status evidence_status selected_donor_ids donor_options requirement_slot_ids remaining_gaps rationale", question_rows),
        "experimental-plan.tsv": ("action_id operation_type donor_id specimen_id design_id context_id batch_id experiment_id replicate_position assay_capability_id arms endpoint_locus_pairs target_sequences control_requirements performed_control_roles performed_qc satisfied_roles dependency_action_ids resource_assignments question_ids rationales", action_rows),
        "resources.tsv": ("resource_id resource_type cost cost_unit cost_scale charged_once activation_rule sharing_key compatibility capacities actual_demands pool_id instance_index max_pool_instances used_by_action_ids", resource_rows),
        "unmet-requirements.tsv": ("kind question_id question_complete donor_id arm_id slot_id action_id planning_status reason detail", gaps),
        "intended-assessments.tsv": ("record_id donor_id design_id locus_id support_status pam_compatible sequence_configuration_id alignment_ids reference_differences reason_codes", intended),
        "variant-effects.tsv": ("record_id donor_id status donor_configuration_id reference_configuration_id donor_alignment_id reference_alignment_id source_variants reason", effects),
    }
    result = {name: _tsv(names.split(), rows) for name, (names, rows) in columns.items()}
    report_format = _report_format(study)
    result.update({"report-format.json": report_format,
        "observations.jsonl": sorted(deepcopy(_records(observations)), key=lambda row: (row["record_type"], row["record_id"])),
        "plan.json": deepcopy(plan), "plan-records.jsonl": [deepcopy(plan["plan_record"])],
        "search-audit.json": deepcopy(search.get("search-audit.json", {})),
        "source-evidence-summary.json": deepcopy(evidence.get("source-evidence-summary.json", [])),
        "paired-source-descriptions.json": deepcopy(comparison.get("paired-source-descriptions.json", [])),
        "sequence-annotations.json": deepcopy(comparison.get("sequence-annotations.json", [])),
        "annotation-manifest.json": deepcopy(comparison.get("annotation-manifest.json", [])),
        "predictions.jsonl": deepcopy(predictions)})
    result["comparison-report.md"] = _markdown(study, reader, index, decisions, design_rows, active,
        evidence_rows, out_of_scope, question_rows, action_rows, resource_rows, intended, effects,
        predictions, gaps, plan, comparison, evidence.get("source-evidence-summary.json", []), report_format)
    domain_assessment = assess_domain(study["config"], search)
    if domain_assessment is not None:
        result["search-audit.json"]["supported_domain"] = domain_assessment
        section = ["## Proposed supported domain", "", "Profile: `" + domain_assessment["profile_id"] + "`.", ""]
        for key, label in (("computational_support", "Computational support"),
                           ("coverage_completeness", "Coverage completeness"),
                           ("biological_applicability", "Biological applicability")):
            claim = domain_assessment[key]
            section += [label + ": **" + claim["status"] + "**. " + "; ".join(claim["reason_codes"]), ""]
        section += ["Declared-interval execution: " + domain_assessment["coverage_completeness"]["declared_scope_execution"] +
                    ". Interval counts cover only the reported reference and alternate scope.", "",
                    "Source measurements, adjudicated outcomes and paired descriptive values retain their own source details.", ""]
        title, rest = result["comparison-report.md"].split("\n", 1)
        result["comparison-report.md"] = title + "\n\n" + "\n".join(section) + rest
    verify_sources(study)
    return result


def _markdown(study, reader, index, decisions, designs, observations, evidence, out_of_scope,
              questions, actions, resources, intended, effects, predictions, gaps, plan, comparison, source_summaries,
              report_format):
    label, safe = reader.label, reader.explain
    lines = ["# Comparative research report", "", "## Research decision", ""]
    recommendations = {
        "observed_difference_in_matched_context": "Review the measured difference in its matched context when prioritizing further testing.",
        "no_difference_detected_within_scope": "No difference was detected in the recorded matched measurements and sequence scope; unresolved questions remain explicit.",
        "sequence_based_hypothesis": "Investigate the genotype-associated sequence concern with the declared experiments.",
        "tradeoff_requires_experiment": "Run a matched comparison to resolve the recorded sequence and evidence tradeoff.",
        "insufficient_comparable_evidence": "Run a matched comparison because current evidence cannot choose between these designs.",
    }
    for decision in decisions:
        lines += [safe(decision["question"]), "", ", ".join(label(d) for d in decision["design_ids"]) + ": **" +
                  decision["state"].replace("_", " ") + "**. " + recommendations[decision["state"]], "",
                  "Recorded reasons: " + safe("; ".join(decision["reason_codes"])) + ".", ""]
    if not decisions:
        lines += ["No current comparative decision record was supplied. The experimental questions below retain their declared purpose.", ""]
    lines += ["## Designs and scope", "", "Execution mode: **" + study["config"]["mode"] + "**; coverage basis: **" + study["config"]["coverage_basis"] + "**.", ""]
    domain = index[study["config"]["policies"]["comparison"]["domain_id"]]
    intervals = domain["intervals"]
    reference = index[domain["reference_bundle_id"]]
    lines += [safe(domain["description"]), "", "Declared reference assembly: **" + safe(reference["assembly"]) + "**.", "",
              f"The shared requested domain contains {len(intervals)} declared intervals and {sum(r['end'] - r['start'] for r in intervals)} interval bases. Per-design evaluated regions, skips and retention are retained in search-audit.json and design-comparison.tsv.", ""]
    if study["config"]["mode"] == "synthetic_test":
        lines += ["This mock example uses the specimens, capabilities and prices declared in its inputs.", ""]
    elif study["config"]["coverage_basis"] == "retained_callset_model":
        lines += ["This public-genome simulation uses hypothetical specimens. Empirical target-donor credit is zero. Modeled sequence support remains separate from measured observations.", ""]
    lines += [_table(["Design", "Editor", "PAM", "Mismatch / DNA gap / RNA gap caps", "Scope", "Donors evaluable / total", "Unresolved donors", "Complete loci", "Complete sequences", "Complete alignments"], [
        [label(r["design_id"]), safe(r["editor"]), ", ".join(r["pam_patterns"]), " / ".join(str(r[k]) for k in ("max_mismatches", "max_dna_bulges", "max_rna_bulges")), r["scope_status"],
         _text(r.get("search_evaluable_donors")) + " / " + _text(r.get("donor_count")), r.get("search_unresolved_donors"),
         *[reader.value({"value": r[k], "reason": r[k + "_reason"]}) for k in ("complete_locus_count", "complete_sequence_configuration_count", "complete_alignment_count")]] for r in designs]),
        "Counts describe the declared sequence search. A missing complete count is unknown, including when some results were retained before a partial search.", "",
        "## Observations, sequence findings and predictions", ""]
    if observations:
        def outcome_cell(o):
            outcome = o.get("adjudicated_outcome")
            if outcome is None:
                return "legacy revision (explicit adjudication required)"
            if outcome == "tested_no_signal":
                return "tested no signal (within stated detection limit)"
            return outcome.replace("_", " ")
        lines += [_table(["Source", "Design", "Context", "Endpoint", "Raw source value", "Normalized value", "Controls / QC",
                          "Adjudicated outcome", "Adjudication reason", "Source outcome (verbatim)", "Detection limit"], [
            [safe(index[o["source_record_id"]]["source_study"]), label(_value(o["design_id"])) if _value(o["design_id"]) is not None else reader.value(o["design_id"]), reader.context(_value(o["context_id"])), label(o["endpoint_id"]),
             reader.value(index[o["source_record_id"]]["raw_value"]) + " " + safe(index[o["source_record_id"]]["raw_unit"]),
             reader.value(o["normalized_value"]) + " " + safe(o["normalized_unit"]), o["control_status"] + " / " + o["qc_status"],
             safe(outcome_cell(o)), safe(o.get("adjudication_reason") or "Legacy version-1 revision requires explicit adjudication; conclusive credit withheld"),
             reader.value(o["source_reported_outcome"]),
             reader.value(o["detection_limit"]["value"]) + " " + safe(o["detection_limit"]["unit"])] for o in observations]),
            "Adjudications are explicit four-state outcomes: observed positive, tested no signal, inconclusive and untested. "
            "Inconclusive and untested records withhold conclusive empirical credit even with a numeric value and passed QC. "
            "No-signal results include their stated detection limit. "
            "Censored or missing values retain their typed outcome and reason but withhold conclusive credit.", ""]
    else:
        lines += ["No active measured observations were supplied. Missing efficacy evidence is unknown.", ""]
    if source_summaries:
        lines += ["Each observation row is a source measurement. Technical repeats do not establish independent donors; missing source-unit metadata remains unknown.", "",
            _table(["Source", "Design", "Endpoint", "Measurements", "Reported independent units", "Measurements with unknown unit"], [
                [safe(row["source_study"]), label(row["design_id"]) if row["design_id"] is not None else "unresolved design",
                 label(row["endpoint_id"]), row["measurements"], row["reported_independent_units"], row["measurements_with_unknown_unit"]]
                for row in source_summaries])]
    classes = Counter(row["classification"] for row in evidence)
    stored_credit = sum(bool(row["empirical_donor_credit"]) for row in evidence)
    current_credit = sum(bool(row["current_empirical_donor_credit"]) for row in evidence)
    superseded = sum(not row["observation_active_revision"] for row in evidence)
    lines += ["Applicability assessments: " + ("; ".join(f"{key.replace('_', ' ')}: {value}" for key, value in sorted(classes.items())) or "none supplied") + ". "
              "Each assessment checks the donor identity, context, endpoint and assay quality separately.", "",
              "Current empirical target-donor credits, counting only active observation revisions that pass the shared outcome and quality gate: " + str(current_credit) + ". "
              "Current study permission must also enable empirical credit in research specimen planning mode. "
              "Stored credit bits in all retained applicability records: " + str(stored_credit) + "; a historical or non-conforming stored bit cannot grant current credit."
              + (" " + str(superseded) + " retained records reference superseded observation revisions and earn no current credit." if superseded else ""), ""]
    for pair in comparison.get("paired-source-descriptions.json", []):
        lines += ["Matched source difference (" + " minus ".join(label(d) for d in pair["design_ids"]) + "): **" + pair["left_minus_right"] + " " + safe(pair["unit"]) + "**, " + label(pair["endpoint_id"]) + "; " + reader.context(pair["context_id"]) + ". " + safe(pair["interpretation"]) + ".", ""]
    changed = [r for r in intended if r["reference_differences"] or _value(r["pam_compatible"]) is not True or r["support_status"] in {"unresolved", "not_evaluated"}]
    lines += ["Intended loci were assessed independently of nomination. The protected intended-assessments.tsv retains every assessment, including lost PAMs, indels and unresolved sequences.", ""]
    if changed:
        lines += [_table(["Donor alias", "Design", "Intended locus", "Support", "PAM compatible", "Differences / reasons"], [
            [label(r["donor_id"]), label(r["design_id"]), label(r["locus_id"]), r["support_status"], reader.value(r["pam_compatible"]),
             str(len(r["reference_differences"])) + " recorded source variant(s); " + safe("; ".join(r["reason_codes"]))] for r in changed])]
    lines += ["Recorded donor/reference sequence effects: " + ("; ".join(f"{key}: {value}" for key, value in sorted(Counter(r["status"] for r in effects).items())) or "none supplied") + ". These are sequence findings.", ""]
    annotations = comparison.get("sequence-annotations.json", [])
    annotation_sources = comparison.get("annotation-manifest.json", [])
    if annotation_sources or any(row["features"] for row in annotations):
        lines += [str(sum(bool(row["features"]) for row in annotations)) + " sequence configurations have supplied gene or regulatory coordinate overlaps. These are coordinate annotations only. Exact overlaps and source versions remain in sequence-annotations.json and annotation-manifest.json.", "",
            _table(["Annotation source", "Release", "Assembly", "Binding"], [
                [safe(row["source"]), safe(row.get("release", "specified by the source record")), safe(row["assembly"]),
                 safe(row.get("sha256", row.get("binding", "unknown")))] for row in annotation_sources])]
    if predictions:
        lines += [_table(["Prediction provider", "Model/version", "Track", "Score", "Uncertainty"], [
            [safe(r["provider"]), safe(r["model"] + "/" + r["model_version"]), safe(r["track"]), reader.value(r["score"]) + " " + safe(r["unit"]), reader.value(r["uncertainty"])] for r in predictions]),
            "Predictions are separate computational records and do not satisfy measured-evidence requirements.", ""]
    else:
        lines += ["No prediction records were supplied.", ""]
    lines += ["Published evidence outside the admitted profile or with unresolved mapping remains in out-of-scope-evidence.tsv: " + str(len(out_of_scope)) + " retained entries.", "",
              "## Proposed experiments", "", f"Plan status: **{plan['status']}**. Completed experimental structures: **{plan['unweighted_completion']}**; frozen priority weight covered: **{plan['weighted_completion']}**. Existing evidence weight: {plan['existing_evidence_weight']}; newly planned weight: {plan['newly_planned_weight']}.", "",
              "Planning and measurement status are shown separately below.", ""]
    for q in questions:
        state = plan["questions"][q["record_id"]]
        arms = "; ".join(safe(a["arm_id"]) + ": " + str(a["required_independent_donors"]) + " independent donor(s)" for a in q["arms"])
        chosen = "; ".join(label(d["donor_id"]) + " in " + safe(d["arm_id"]) for d in state["donors"] if d["donor_id"] in state["selected_donor_ids"])
        lines += ["**" + label(q["record_id"]) + " — " + q["template"].replace("_", " ") + "**. " + safe(q["hypothesis"]), "",
                  safe(q["decision"] + " " + q["rationale"]), "",
                  "Requirements: " + arms + "; " + str(q["technical_replicates"]) + " technical measurement(s) per endpoint; " + q["pairing"].replace("_", " ") + ". Selected donors: " + (chosen or "none") + ".", "",
                  "Planning: **" + q["planning_status"] + "**; evidence: **" + q["evidence_status"] + "**; priority: " + str(q["priority_weight"]) + "; hard required: " + _text(q["hard_required"]) + ".", ""]
    if actions:
        lines += [_table(["Operation", "Donor alias / arm", "Design", "Type", "Readouts", "Actual target DNA", "QC / controls", "Prerequisites"], [
            [label(a["action_id"]), label(a["donor_id"]) + "; " + safe(", ".join(sorted({x["arm_id"] for x in a["arms"]}))), label(a["design_id"]), a["operation_type"].replace("_", " "),
             "; ".join(label(p["endpoint_id"]) + " at " + label(p["locus_id"]) for p in a["endpoint_locus_pairs"]),
             "; ".join(label(t["locus_id"]) + ": " + ", ".join(s["sequence"] for s in t["sequences"]) for t in a["target_sequences"]),
             safe(", ".join(a["performed_qc"] + a["performed_control_roles"])) or "preparation only", ", ".join(label(d) for d in a["dependency_action_ids"]) or "none"] for a in actions])]
    else:
        lines += ["No new operations were selected.", ""]
    subtotals = defaultdict(int)
    for row in resources:
        subtotals[row["resource_type"]] += row["cost"]
    lines += ["## Cost and unresolved requirements", "", _table(["Resource package", "Distinct instance subtotal"], sorted(subtotals.items())),
              "Total: **" + str(plan["accounting"]["total_cost"]) + " " + safe(plan["accounting"]["cost_unit"]) + "** against budget **" + str(plan["budget"]) + "** (declared cost scale " + str(plan["accounting"]["cost_scale"]) + "). Each active named resource is charged once; resources.tsv lists its actual consumers and capacity demands.", "",
              "Selector: " + plan["selector"] + "; global optimum established: " + _text(plan["optimality_established"]) + "; requested mandatory questions met: " + _text(plan["meets_request"]) + ". Search limits and any infeasibility certificate remain in plan.json.", ""]
    if plan["construction_limited_question_ids"]:
        lines += [str(len(plan["construction_limited_question_ids"])) + " questions were excluded from action construction by the declared limit and remain unplanned.", ""]
    if gaps:
        lines += ["Recorded gaps: " + safe("; ".join(sorted({row["reason"] for row in gaps}))) + ". Candidate-slot rows in unmet-requirements.tsv describe donor alternatives; they do not require every alternative donor to be selected.", ""]
    else:
        lines += ["The selected plan has no recorded requirement gaps. Returned measurements and their controls still require reconciliation.", ""]
    code = report_format["producing_code"]
    lines += ["## Report format and human review records", "",
              "Report format version: **" + report_format["format_version"] + "**; observation outcome contract: `" +
              report_format["observation_outcome_contract"] + "`; accepted observation history versions: " +
              ", ".join(report_format["accepted_observation_history_versions"]) + ". "
              "The machine-readable copy is [report-format.json](report-format.json).", "",
              "Producing code identity captured when the study was resolved: package **" + safe(code["package_version"]) +
              "**, source SHA-256 `" + safe(code["source_sha256"]) + "`, Python " + safe(code["python"]) +
              ". It is not recomputed at render time.", ""]
    for key, label in HUMAN_RECORD_FIELDS:
        record = report_format["human_records"][key]
        lines += [label.capitalize() + ": **not recorded**. " + safe(record["reason"]) + ".", ""]
    lines += ["Detailed donor, sequence and source records remain protected local research files. Donor aliases and aggregation do not guarantee anonymity.", "",
              "Exports: [design comparison](design-comparison.tsv), [questions](questions.tsv), [experimental operations](experimental-plan.tsv), [resources](resources.tsv), [unmet requirements](unmet-requirements.tsv), [observations](observations.jsonl), [evidence applicability](evidence-applicability.tsv), [out-of-scope evidence](out-of-scope-evidence.tsv), [plan](plan.json), [report format](report-format.json). The stage manifest binds source hashes, policies, code and verified ancestors.", ""]
    return "\n".join(lines)
