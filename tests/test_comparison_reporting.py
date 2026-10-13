"""Reports preserve measured values, incomplete scope and complete plan costs."""
from copy import deepcopy
import csv
import io
import json
from pathlib import Path

import pytest

from base_edit_safety.population.common import PopulationError
from base_edit_safety.population.comparison import local_search
from base_edit_safety.population.comparison_evidence import compare_with_evidence, reconcile_evidence
from base_edit_safety.population.comparison_provenance import resolve_study
from base_edit_safety.population.comparison_questions import prepare_question_set
from base_edit_safety.population.comparison_reporting import render_report
from base_edit_safety.population.experiment_planning import plan_experiments, verify_plan_revision
from base_edit_safety.population.observations import normalize_observations
from planning_cases import full_fixture, alternatives_case, remake, write_case
from test_comparison_questions import observation_records

ROOT = Path(__file__).resolve().parents[1]


def table(output, name):
    return list(csv.DictReader(io.StringIO(output[name]), delimiter="\t"))


def render(qset, *, plan=None, **artifacts):
    values = {name: {} for name in ("search", "observations", "evidence", "comparison", "questions")}
    values.update(artifacts)
    plan = plan_experiments(qset) if plan is None else plan
    verify_plan_revision(qset, plan)
    return render_report(qset.study, **values, plan=plan, question_set=qset)


@pytest.fixture
def fixture(tmp_path):
    return full_fixture(tmp_path / "section9")


@pytest.fixture(scope="module")
def sequence_case():
    study = resolve_study(ROOT / "examples/comparison/synthetic/study.json")
    search = local_search(study)
    observations = {"observations.jsonl": []}
    evidence = reconcile_evidence(study, search, observations)
    comparison = compare_with_evidence(study, search, evidence, observations)
    rows = [row for artifacts in (search, observations, evidence, comparison)
            for name, records in artifacts.items() if name.endswith(".jsonl") for row in records]
    qset = prepare_question_set(study, rows)
    return qset, {"search": search, "observations": observations, "evidence": evidence, "comparison": comparison}


def test_section9_exports_unique_costs_actual_noncarrier_dna_and_separate_states(fixture):
    qset = fixture["question_set"]
    before = deepcopy(qset.records)
    plan = plan_experiments(qset)
    original_plan = deepcopy(plan)
    output = render(qset, plan=plan)
    required = {"comparison-report.md", "design-comparison.tsv", "donor-design-matrix.tsv", "observations.jsonl",
        "evidence-applicability.tsv", "out-of-scope-evidence.tsv", "questions.tsv", "experimental-plan.tsv",
        "resources.tsv", "unmet-requirements.tsv", "plan.json", "plan-records.jsonl", "search-audit.json",
        "report-format.json"}
    assert required <= output.keys()
    assert "run-manifest.json" not in output
    resources = table(output, "resources.tsv")
    assert len({r["resource_id"] for r in resources}) == len(resources)
    assert sum(int(r["cost"]) for r in resources) == 18
    operations = table(output, "experimental-plan.tsv")
    assert len(operations) == 10
    assert sum(r["operation_type"] == "measurement_panel" for r in operations) == 4
    assert sum(r["operation_type"] == "mock_control" for r in operations) == 2
    for operation in operations:
        assert json.loads(operation["question_ids"])
        if operation["donor_id"] == fixture["labels"]["donors"]["N"]:
            target = next(t for t in json.loads(operation["target_sequences"]) if t["locus_id"] == fixture["labels"]["loci"]["L"])
            assert [s["sequence"] for s in target["sequences"]] == ["ACGTACGTACGTACGTACGTAGG"]
    question = table(output, "questions.tsv")[0]
    assert question["planning_status"] == "fully_planned" and question["evidence_status"] == "no_measurements"
    assert output["plan.json"] == original_plan and output["plan-records.jsonl"] == [plan["plan_record"]]
    assert qset.records == before and plan == original_plan


def test_report_format_declares_contract_captured_identity_and_missing_human_records(fixture):
    qset = fixture["question_set"]
    output = render(qset)
    record = output["report-format.json"]
    assert record["format_version"] == "2"
    assert record["observation_outcome_contract"] == "comparison_observation_revision:v2"
    assert record["accepted_observation_history_versions"] == ["1", "2"]
    assert record["producing_code"] == qset.study["code"]
    assert set(record["human_records"]) == {"reviewer_attestation", "funding_disclosure", "customer_disclosure"}
    for row in record["human_records"].values():
        assert row["status"] == "not_recorded"
        assert "is present in the report inputs" in row["reason"]
        assert "none" not in row["reason"].lower().replace("not ", "")
    text = output["comparison-report.md"]
    assert "Report format version: **2**" in text
    assert "`comparison_observation_revision:v2`" in text
    assert qset.study["code"]["source_sha256"] in text
    assert text.count("**not recorded**") == 3
    assert "[report-format.json](report-format.json)" in text


def test_report_format_keeps_a_deliberately_different_captured_study_identity(fixture, tmp_path):
    from base_edit_safety.population.comparison_questions import prepare_question_set
    from base_edit_safety.population.provenance import code_identity
    qset = write_case(fixture, tmp_path / "other-code", deepcopy(fixture["question_set"].records))
    other = {"package_version": "0.0.0-review-fixture", "source_sha256": "f" * 64, "python": "0.0.0"}
    assert other != code_identity()
    study = resolve_study(tmp_path / "other-code" / "study.json", code=other)
    output = render(prepare_question_set(study, qset.records))
    record = output["report-format.json"]
    assert record["producing_code"] == other
    assert record["producing_code"] != code_identity()
    text = output["comparison-report.md"]
    assert "0.0.0-review-fixture" in text and "f" * 64 in text
    assert qset.study["code"]["source_sha256"] not in text


def test_budget17_keeps_the_failure_and_does_not_publish_a_successful_contrast(fixture):
    qset = fixture["question_set"]
    plan = plan_experiments(qset, budget=17, selector="exact_tiny")
    output = render(qset, plan=plan)
    assert not output["plan.json"]["meets_request"]
    assert table(output, "experimental-plan.tsv") == []
    assert table(output, "questions.tsv")[0]["planning_status"] == "infeasible_proven"
    assert any(row["kind"] == "mandatory_question" for row in table(output, "unmet-requirements.tsv"))
    assert "Completed experimental structures: **0**" in output["comparison-report.md"]


def test_only_current_comparison_summaries_are_exported_and_unknown_is_not_zero(sequence_case):
    qset, artifacts = sequence_case
    output = render(qset, **artifacts)
    current = [r for r in artifacts["comparison"]["comparison-summaries.jsonl"] if r["record_type"] == "comparison_design_summary"]
    earlier = [r for r in artifacts["search"]["comparison-records.jsonl"] if r["record_type"] == "comparison_design_summary"]
    assert {r["record_id"] for r in current} != {r["record_id"] for r in earlier}
    exported = table(output, "design-comparison.tsv")
    assert len(exported) == 2
    assert {r["design_id"] for r in exported} == {r["design_id"] for r in current}
    assert any(r["complete_alignment_count"] == "null" and r["complete_alignment_count_reason"] != "null" for r in exported)
    assert "unknown (incomplete donor search scope)" in output["comparison-report.md"]
    assert len(table(output, "donor-design-matrix.tsv")) == 6


def test_intended_pam_loss_and_indel_survive_without_nomination(sequence_case):
    qset, artifacts = sequence_case
    output = render(qset, **artifacts)
    intended = table(output, "intended-assessments.tsv")
    assert any(json.loads(row["pam_compatible"])["value"] is False for row in intended)
    assert any(row["support_status"] == "unresolved" for row in intended)
    assert any(json.loads(row["alignment_ids"]) == [] and row["support_status"] == "supported" for row in intended)
    text = output["comparison-report.md"]
    assert "intended_alignment_outside_nomination_profile" in text
    assert "Intended loci were assessed independently of nomination" in text


def test_duplicate_current_summary_is_an_error_instead_of_a_double_count(sequence_case):
    qset, artifacts = sequence_case
    changed = deepcopy(artifacts)
    summary = next(r for r in changed["comparison"]["comparison-summaries.jsonl"] if r["record_type"] == "comparison_design_summary")
    changed["comparison"]["comparison-summaries.jsonl"].append(remake(summary, reason_codes=summary["reason_codes"] + ["another_current_revision"]))
    with pytest.raises(PopulationError, match="multiple current summaries"):
        render(qset, **changed)


def test_source_decimal_and_inapplicable_evidence_remain_visible():
    study = resolve_study(ROOT / "examples/comparison/observations/study.json")
    observations = normalize_observations(study)
    evidence = reconcile_evidence(study, {}, observations)
    qset = prepare_question_set(study, observations["observations.jsonl"] + evidence["evidence-applicability.jsonl"])
    evidence["published-alignment-replay.json"] = [{"source_key": "artificial-source-positive", "design_id": study["config"]["policies"]["search"]["design_ids"][0],
        "source_positive": True, "source_reported_label": "positive", "disposition": "outside_declared_search_policy",
        "reason_codes": ["source_bulges_exceed_profile"], "alignment": {"guide": "A-C", "target": "AGC"},
        "exact_alignment_ids": [], "related_alignment_ids": [], "empirical_target_donor_credit": False}]
    output = render(qset, observations=observations, evidence=evidence)
    assert output["observations.jsonl"] == sorted(observations["observations.jsonl"], key=lambda r:(r["record_type"],r["record_id"]))
    text = output["comparison-report.md"]
    assert "9.8000000000000007" in text and "0.098000000000000007" in text
    rows = table(output, "evidence-applicability.tsv")
    assert rows and all(row["empirical_donor_credit"] == "false" for row in rows)
    excluded = table(output, "out-of-scope-evidence.tsv")[0]
    assert excluded["source_positive"] == "true" and "source_bulges_exceed_profile" in excluded["reason_codes"]
    assert excluded["empirical_target_donor_credit"] == "false"


def test_default_markdown_redacts_identifiers_and_paths_but_protected_tables_keep_targets(fixture, tmp_path):
    rows = deepcopy(fixture["question_set"].records)
    identifier = fixture["labels"]["donors"]["C"]
    for row in rows:
        if row["record_type"] == "comparison_question":
            row["hypothesis"] += " synthetic-donor-C synthetic-sample-C " + identifier + " /private/protected-study/notes.txt C:\\protected\\notes.txt file:///private/other-study/notes.txt"
    qset = write_case(fixture, tmp_path / "private-case", rows)
    output = render(qset)
    text = output["comparison-report.md"]
    for secret in (identifier, "synthetic-donor-C", "synthetic-sample-C", "/private/protected-study", "C:\\protected", "file:///private/other-study"):
        assert secret not in text
    assert "Donor 1" in text and "local path withheld" in text
    assert identifier in output["experimental-plan.tsv"]
    assert "ACGTGCGTACGTACGTACGTAGG" in text


def test_actual_alternative_resource_instances_are_charged_once(fixture, tmp_path):
    qset = alternatives_case(fixture, tmp_path / "alternatives")
    output = render(qset)
    resources = table(output, "resources.tsv")
    sequencing = [row for row in resources if row["resource_type"] == "sequencing_batch"]
    assert sorted(int(row["cost"]) for row in sequencing) == [3, 5]
    assert sum(int(row["cost"]) for row in resources) == 26
    used = {r["resource_id"] for r in sequencing}
    actual = {d["resource_id"] for row in table(output, "experimental-plan.tsv")
              for d in json.loads(row["resource_assignments"]) if d["assignment_key"] == "sequencing"}
    assert actual == used


def test_reordered_record_inputs_leave_readable_scientific_exports_unchanged(sequence_case):
    qset, artifacts = sequence_case
    original = render(qset, **artifacts)
    reordered = deepcopy(artifacts)
    for value in reordered.values():
        for name, rows in value.items():
            if name.endswith(".jsonl"):
                rows.reverse()
    changed = render(qset, **reordered)
    for name in original:
        if name.endswith((".tsv", ".md")):
            assert changed[name] == original[name]


def test_changed_plan_audit_is_rejected_before_rendering(fixture):
    qset = fixture["question_set"]
    plan = plan_experiments(qset)
    plan["weighted_completion"] = 999
    with pytest.raises(PopulationError, match="identity or audit|completion weights"):
        render_report(qset.study, search={}, observations={}, evidence={}, comparison={}, questions={}, plan=plan, question_set=qset)


def test_source_bytes_are_rechecked_before_reporting(fixture):
    qset = fixture["question_set"]
    plan = plan_experiments(qset)
    source = qset.study["config_path"].parent / "synthetic-genomes.json"
    source.write_text(source.read_text() + " ")
    with pytest.raises(PopulationError, match="source changed"):
        render_report(qset.study, search={}, observations={}, evidence={}, comparison={}, questions={}, plan=plan, question_set=qset)


def test_numeric_subject_identifiers_do_not_change_scientific_values(fixture, tmp_path):
    rows = deepcopy(fixture["question_set"].records) + observation_records(fixture, credit=False)
    for row in rows:
        if row["record_type"] == "comparison_donor":
            row["sample_ids"].append("0.2")
    qset = write_case(fixture, tmp_path / "numeric-subject", rows)
    observations = {"observations.jsonl": [r for r in qset.records if r["record_type"] in {
        "comparison_source_observation", "comparison_observation_revision", "comparison_observation_state"}]}
    output = render(qset, observations=observations)
    assert "0.2 fraction" in output["comparison-report.md"]
    assert "**18 simulation_units**" in output["comparison-report.md"]


def test_public_mode_preserves_zero_empirical_credit():
    study = resolve_study(ROOT / "examples/comparison/workflow/mode-declarations/public-positive/study.json")
    qset = prepare_question_set(study, [])
    output = render(qset)
    assert "Empirical target-donor credit is zero" in output["comparison-report.md"]
    assert table(output, "questions.tsv")[0]["evidence_status"] == "no_measurements"


def test_coordinate_annotations_are_exported_without_becoming_predictions(sequence_case):
    qset, artifacts = sequence_case
    supplied = deepcopy(artifacts)
    sequence = next(r for r in qset.records if r["record_type"] == "comparison_sequence_configuration")
    annotations = [{"sequence_configuration_id": sequence["record_id"], "features": [
        {"feature_type": "regulatory", "name": "artificial region", "contig": sequence["contig"],
         "start": 0, "end": 23, "source_id": "Artificial annotation", "release": "fixture-v1"}],
        "interpretation": "Coordinate annotation only; no activity or outcome estimate"}]
    manifest = [{"source": "Artificial annotation", "release": "fixture-v1", "assembly": "artificial", "sha256": "a" * 64}]
    supplied["comparison"].update({"sequence-annotations.json": annotations, "annotation-manifest.json": manifest})
    output = render(qset, **supplied)
    assert output["sequence-annotations.json"] == annotations and output["annotation-manifest.json"] == manifest
    assert "coordinate annotations only" in output["comparison-report.md"]
    assert "Exact overlaps and source versions remain in sequence-annotations.json and annotation-manifest.json" in output["comparison-report.md"]
    assert output["predictions.jsonl"] == []


def _credit_case(fixture, directory, rows, extra=(), *, mode="research_sample_planning", allow=True):
    # Declared preconditions in invented software fixtures, not actual specimens or assays.
    rows = deepcopy(rows)
    if mode == "research_sample_planning":
        for row in rows:
            if row["record_type"] == "comparison_specimen":
                row["availability"] = "verified"
            elif row["record_type"] == "comparison_assay_capability":
                row["feasibility"] = "verified"
    qset = write_case(fixture, directory, rows)
    path = directory / "study.json"
    config = json.loads(path.read_text())
    config.update(mode=mode, coverage_basis="callable_supported" if mode == "research_sample_planning" else "synthetic_truth")
    config["policies"]["evidence"]["allow_empirical_donor_credit"] = allow
    path.write_text(json.dumps(config, indent=2) + "\n")
    from base_edit_safety.population.comparison_provenance import resolve_study
    from base_edit_safety.population.comparison_questions import prepare_question_set
    return prepare_question_set(resolve_study(path), list(extra))


def _outcome_bundle(fixture, outcome, *, value, credit, source_text=None, limit=0.01):
    from base_edit_safety.population.common import known, unknown
    source, revision, _state = observation_records(fixture, export=outcome, credit=credit)
    if value != 0.2:
        source = remake(source, raw_value=known(str(value)),
                        raw_row_text="Artificial measurement: " + str(value) + " fraction")
    revision = remake(revision, source_record_id=source["record_id"], adjudicated_outcome=outcome,
        adjudication_reason="Artificial acceptance adjudication: " + outcome,
        source_reported_outcome=known(outcome if source_text is None else source_text),
        normalized_value=known(value),
        detection_limit={"value": known(limit) if limit is not None else unknown("Not reported"), "unit": "fraction"})
    state = remake(_state, source_record_id=source["record_id"],
                   active_revision_id=revision["record_id"], revision_ids=[revision["record_id"]])
    return [source, revision, state]


def _targeted_evidence(qset, revisions):
    from base_edit_safety.population.comparison_evidence import assess_observation
    config = qset.study["config"]
    applications = [assess_observation(revision, qset.index,
        design_id=revision["design_id"]["value"], context_id=revision["context_id"]["value"],
        endpoint_id=revision["endpoint_id"], donor_id=revision["subject"]["donor_id"]["value"],
        sequence_configuration_id=revision["sequence_configuration_id"]["value"],
        mode=config["mode"], coverage_basis=config["coverage_basis"], policy=config["policies"]["evidence"])
        for revision in revisions]
    return {"evidence-applicability.jsonl": applications}


def _render_with_observations(fixture, directory, bundles):
    qset = _credit_case(fixture, directory, deepcopy(fixture["question_set"].records), extra=bundles)
    artifacts = {"observations.jsonl": bundles}
    revisions = [row for row in bundles if row["record_type"] == "comparison_observation_revision"]
    return qset, render(qset, observations=artifacts, evidence=_targeted_evidence(qset, revisions))


def test_report_exports_all_four_outcomes_with_verbatim_source_and_limits(fixture, tmp_path):
    bundles = []
    for outcome, credit in [("observed_positive", True), ("tested_no_signal", True),
                            ("inconclusive", False), ("untested", False)]:
        value = 0.0 if outcome == "tested_no_signal" else 0.2
        bundles += _outcome_bundle(fixture, outcome, value=value, credit=credit)
    qset = _credit_case(fixture, tmp_path / "four-states", deepcopy(fixture["question_set"].records), extra=bundles)
    artifacts = {"observations.jsonl": bundles}
    revisions = [row for row in bundles if row["record_type"] == "comparison_observation_revision"]
    applications = _targeted_evidence(qset, revisions)["evidence-applicability.jsonl"]
    assert sum(row["assessment"]["empirical_donor_credit"] for row in applications) == 2
    output = render(qset, observations=artifacts, evidence={"evidence-applicability.jsonl": applications})
    rows = table(output, "evidence-applicability.tsv")
    by_outcome = {row["adjudicated_outcome"]: row for row in rows}
    assert set(by_outcome) == {"observed_positive", "tested_no_signal", "inconclusive", "untested"}
    for outcome, row in by_outcome.items():
        assert json.loads(row["source_reported_outcome"])["value"] == outcome
        assert row["adjudication_reason"].endswith(outcome)
        assert json.loads(row["detection_limit_value"])["value"] == 0.01
        assert row["detection_limit_unit"] == "fraction"
        assert row["normalized_value_reason"] == "null"
        assert row["observation_active_revision"] == "true"
    assert by_outcome["observed_positive"]["current_empirical_donor_credit"] == "true"
    assert by_outcome["tested_no_signal"]["current_empirical_donor_credit"] == "true"
    for outcome in ("inconclusive", "untested"):
        assert by_outcome[outcome]["current_empirical_donor_credit"] == "false"
        assert json.loads(by_outcome[outcome]["normalized_value"])["value"] == 0.2
        assert "outcome_" + outcome in json.loads(by_outcome[outcome]["current_gate_reason_codes"])
    text = output["comparison-report.md"]
    assert "tested no signal (within stated detection limit)" in text
    assert "Inconclusive and untested records withhold conclusive empirical credit" in text
    assert "Current empirical target-donor credits, counting only active observation revisions that pass the shared outcome and quality gate: 2." in text
    assert "Stored credit bits in all retained applicability records: 2" in text


def test_incompatible_limits_and_contradictory_labels_stay_visible_without_credit(fixture, tmp_path):
    bundles = _outcome_bundle(fixture, "tested_no_signal", value=0.005, credit=False, limit=None)
    bundles += _outcome_bundle(fixture, "observed_positive", value=0.2, credit=False, source_text="untested")
    qset, output = _render_with_observations(fixture, tmp_path / "visible-failures", bundles)
    rows = table(output, "evidence-applicability.tsv")
    assert len(rows) == 2
    missing_limit = next(row for row in rows if row["adjudicated_outcome"] == "tested_no_signal")
    assert json.loads(missing_limit["detection_limit_value"])["value"] is None
    assert "no_signal_without_usable_detection_limit" in json.loads(missing_limit["current_gate_reason_codes"])
    assert missing_limit["current_empirical_donor_credit"] == "false"
    contradicted = next(row for row in rows if row["adjudicated_outcome"] == "observed_positive")
    assert json.loads(contradicted["source_reported_outcome"])["value"] == "untested"
    assert "adjudication_contradicts_source_outcome" in json.loads(contradicted["current_gate_reason_codes"])
    assert contradicted["current_empirical_donor_credit"] == "false"
    assert "quality gate: 0." in output["comparison-report.md"]


def test_superseded_and_legacy_stored_credit_cannot_enter_current_totals(fixture, tmp_path):
    from base_edit_safety.population.comparison_contracts import make_record
    from base_edit_safety.population.common import known
    source, revision, state = observation_records(fixture, export="legacy")
    legacy_fields = {key: deepcopy(value) for key, value in revision.items() if key not in {
        "record_type", "schema_version", "record_id", "observation_id", "adjudicated_outcome", "adjudication_reason"}}
    legacy_fields.update(normalization_schema_version="1", source_reported_outcome=known("untested"))
    legacy = make_record("observation_revision", version="1", **legacy_fields)
    newer = remake(revision, adjudicated_outcome="inconclusive",
        adjudication_reason="Latest artificial review replaces the earlier reading",
        applicability={**deepcopy(revision["applicability"]), "empirical_donor_credit": False})
    state = remake(state, active_revision_id=newer["record_id"], revision_ids=[legacy["record_id"], newer["record_id"]])
    bundles = [source, legacy, newer, state]
    qset = _credit_case(fixture, tmp_path / "superseded", deepcopy(fixture["question_set"].records), extra=bundles)
    artifacts = {"observations.jsonl": bundles}
    evidence = _targeted_evidence(qset, [newer])
    # Simulate a dishonest historical application: stored credit pointing at the superseded v1 revision.
    current_app = next(row for row in evidence["evidence-applicability.jsonl"])
    stale_assessment = {**deepcopy(current_app["assessment"]), "classification": "exact_donor_configuration",
                        "empirical_donor_credit": True, "quality": {"status": "applicable", "reason_codes": []}}
    stale_app = remake(current_app, observation_id=legacy["record_id"], assessment=stale_assessment)
    evidence["evidence-applicability.jsonl"].append(stale_app)
    output = render(qset, observations=artifacts, evidence=evidence)
    rows = table(output, "evidence-applicability.tsv")
    stale_row = next(row for row in rows if row["observation_id"] == legacy["record_id"])
    assert stale_row["adjudicated_outcome"] == "legacy_unadjudicated"
    assert stale_row["empirical_donor_credit"] == "true"
    assert stale_row["observation_active_revision"] == "false"
    assert stale_row["current_empirical_donor_credit"] == "false"
    active_row = next(row for row in rows if row["observation_id"] == newer["record_id"])
    assert active_row["current_empirical_donor_credit"] == "false"
    text = output["comparison-report.md"]
    assert "quality gate: 0." in text
    assert "Stored credit bits in all retained applicability records: 1" in text
    assert "1 retained records reference superseded observation revisions and earn no current credit." in text


@pytest.mark.parametrize("mode,allow,reason", [
    ("synthetic_test", True, "study_mode_disallows_empirical_donor_credit"),
    ("research_sample_planning", False, "study_policy_disallows_empirical_donor_credit"),
])
def test_current_study_permission_cannot_reuse_stored_research_credit(fixture, tmp_path, mode, allow, reason):
    bundles = _outcome_bundle(fixture, "observed_positive", value=0.2, credit=True)
    previous = _credit_case(fixture, tmp_path / "previous", fixture["question_set"].records, extra=bundles)
    revisions = [row for row in bundles if row["record_type"] == "comparison_observation_revision"]
    evidence = _targeted_evidence(previous, revisions)
    assert evidence["evidence-applicability.jsonl"][0]["assessment"]["empirical_donor_credit"] is True
    current = _credit_case(fixture, tmp_path / "current", fixture["question_set"].records,
                           extra=bundles, mode=mode, allow=allow)
    output = render(current, observations={"observations.jsonl": bundles}, evidence=evidence)
    row = table(output, "evidence-applicability.tsv")[0]
    assert row["empirical_donor_credit"] == "true"
    assert row["observation_active_revision"] == "true"
    assert row["current_empirical_donor_credit"] == "false"
    assert reason in json.loads(row["current_gate_reason_codes"])
    assert "quality gate: 0." in output["comparison-report.md"]
    assert "Stored credit bits in all retained applicability records: 1" in output["comparison-report.md"]
