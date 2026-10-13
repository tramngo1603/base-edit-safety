"""Check accepted and rejected mock inputs for the proposed domain."""
from copy import deepcopy
import csv
import io
import json
from pathlib import Path
import shutil

import pytest

from base_edit_safety.pipeline.comparison import main
from base_edit_safety.population.common import PopulationError, checksum
from base_edit_safety.population.comparison_contracts import validate_study
from base_edit_safety.population.comparison_domain import (
    PROFILE_ID, DomainAbstention, assess_domain, profile_definition,
)
from base_edit_safety.population.comparison_evidence import assess_observation
from base_edit_safety.population.comparison_provenance import effective_inputs, resolve_study
from test_comparison_evidence import POLICY, ROWS, credited_observation, one, remake
from test_comparison_workflow import _reidentify, _remap, invoke

ROOT = Path(__file__).resolve().parents[1]


def write(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")
    return path


@pytest.fixture
def profile_case(tmp_path):
    folder = tmp_path / "artificial-domain"
    shutil.copytree(ROOT / "examples/comparison/workflow", folder)
    config = json.loads((folder / "study.json").read_text())
    before = deepcopy(config)
    definition = profile_definition()
    context = next(r for r in config["catalogs"] if r["record_type"] == "comparison_experimental_context")
    for key, value in definition["context"].items():
        context[key] = {"value": value, "reason": None} if key == "cell_state" else value
    for row in config["catalogs"]:
        if row["record_type"] == "comparison_editor_definition":
            row.update(deepcopy(definition["editor"]))
            row["editing_window"] = deepcopy(definition["editing_window"])
        elif row["record_type"] == "comparison_search_profile":
            row.update(max_mismatches=6, max_total_edits=6)
    designs = [r for r in config["catalogs"] if r["record_type"] == "comparison_design"]
    # Distinct mock guides keep two designs after choosing one exact editor.
    designs[1]["guide_sequence"] = "T" + designs[1]["guide_sequence"][1:]
    endpoints = [r["record_id"] for r in config["catalogs"] if r["record_type"] == "comparison_endpoint_definition"]
    config["policies"]["evidence"]["supported_domain"] = {
        "profile_id": PROFILE_ID, "context_id": context["record_id"], "endpoint_ids": endpoints}
    config = _reidentify(config)
    mapping = {a["record_id"]: b["record_id"] for a, b in zip(before["catalogs"], config["catalogs"])}
    config["catalogs"] = list({r["record_id"]: r for r in config["catalogs"]}.values())
    for ref in config["input_refs"]:
        if ref["role"] in {"extraction_map", "experiment_policy"}:
            path = folder / ref["path"]
            write(path, _remap(json.loads(path.read_text()), mapping))
            ref["sha256"] = checksum(path)
    path = write(folder / "study.json", config)
    validate_study(config)
    return path, config


def row(config, kind):
    return next(r for r in config["catalogs"] if r["record_type"] == "comparison_" + kind)


def changed_case(config, change):
    config = deepcopy(config)
    editor, context = row(config, "editor_definition"), row(config, "experimental_context")
    profile, endpoint = row(config, "search_profile"), row(config, "endpoint_definition")
    if change == "editor":
        editor["variant_name"] = "ABE8e-SpCas9"
    elif change == "cas":
        editor["cas_component"] = "SpCas9-NG"
    elif change == "conversion":
        editor["conversion"].update(substrate_base="C", product_base="T")
    elif change == "context":
        context["cell_type"] = "HSPC"
    elif change == "material":
        context["cell_state"]["value"] = "cell-free"
    elif change == "delivery":
        context["delivery_method"] = "RNP electroporation"
    elif change == "endpoint":
        endpoint["name"] = "HbF_fraction"
    elif change == "mismatches":
        profile.update(max_mismatches=7, max_total_edits=7)
    elif change == "total":
        profile["max_total_edits"] = 7
    elif change == "dna_bulge":
        profile["max_dna_bulges"] = 1
    elif change == "rna_bulge":
        profile["max_rna_bulges"] = 1
    elif change == "pam":
        editor["pam_patterns"] = ["NAG"]
        for candidate in config["catalogs"]:
            if candidate["record_type"] == "comparison_search_profile":
                candidate["pam_patterns"] = ["NAG"]
    elif change == "strand":
        profile["strands"] = ["+"]
    elif change == "window":
        editor["editing_window"]["start_1based"] = 1
    elif change == "native":
        profile["backend_name"], profile["backend_version"] = "cas_offinder", "2.4.1"
    elif change in {"chrX", "chrY", "chrM"}:
        row(config, "domain")["intervals"][0]["contig"] = change
    elif change == "prediction":
        config["policies"]["predictions"]["enabled"] = True
    return _reidentify(config)


@pytest.mark.parametrize("change,reason", [
    ("editor", "domain_unsupported_editor"), ("cas", "domain_unsupported_editor"),
    ("conversion", "domain_unsupported_editor"), ("context", "domain_unsupported_context"),
    ("material", "domain_unsupported_context"), ("delivery", "domain_unsupported_context"),
    ("endpoint", "domain_unsupported_endpoint"), ("mismatches", "domain_unsupported_search_settings"),
    ("total", "domain_unsupported_search_settings"), ("dna_bulge", "domain_unsupported_search_settings"),
    ("rna_bulge", "domain_unsupported_search_settings"), ("pam", "domain_unsupported_search_settings"),
    ("strand", "domain_unsupported_search_settings"), ("window", "domain_unsupported_window_annotation"),
    ("native", "domain_native_backend_not_integrated"),
    ("chrX", "domain_sex_chromosome_reconstruction_unsupported"),
    ("chrY", "domain_sex_chromosome_reconstruction_unsupported"),
    ("chrM", "domain_nonprimary_contig_excluded"), ("prediction", "domain_predictions_excluded"),
])
def test_unsupported_requests_abstain_before_execution(profile_case, change, reason, capsys, tmp_path):
    path, config = profile_case
    config = changed_case(config, change)
    with pytest.raises(DomainAbstention) as exc:
        validate_study(config)
    assert reason in exc.value.assessment["request"]["reason_codes"]
    write(path, config)
    output = tmp_path / "must-not-exist"
    assert main(["search", "--config", str(path), "--output", str(output)]) == 2
    result = json.loads(capsys.readouterr().err)
    assert result["status"] == "abstained" and not output.exists()
    assert reason in result["supported_domain"]["request"]["reason_codes"]


def test_no_caller_approval_override_and_no_silent_profile_fallback(profile_case):
    _, config = profile_case
    for field, value in (("approved", True), ("profile_id", "unreviewed-v2")):
        changed = deepcopy(config)
        changed["policies"]["evidence"]["supported_domain"][field] = value
        with pytest.raises(PopulationError):
            validate_study(changed)


def test_scope_completion_does_not_promote_coverage_or_biology(profile_case):
    _, config = profile_case
    assessment = assess_domain(config)
    assert assessment["computational_support"]["status"] == "supported_for_artificial_checks"
    assert not assessment["release_complete"]
    assert "domain_callability_not_established" in assessment["coverage_completeness"]["reason_codes"]
    assert "domain_alternate_xy_ploidy_unsupported" in assessment["coverage_completeness"]["reason_codes"]
    assert "domain_exposure_conflict_unresolved" in assessment["biological_applicability"]["reason_codes"]
    assert "domain_dose_unresolved" in assessment["biological_applicability"]["reason_codes"]
    # Even a declared all-contig search and a caller's completion flag cannot close a genome claim.
    changed = deepcopy(config)
    row(changed, "domain")["intervals"] = [{"contig": c, "start": 0, "end": 100} for c in assessment["coverage_completeness"]["required_primary_contigs"]]
    fake = {"rows.jsonl": [{"record_type": "comparison_search_scope", "execution_status": "complete"}]}
    result = assess_domain(changed, fake)
    assert result["coverage_completeness"]["status"] == "incomplete"
    assert result["biological_applicability"]["status"] == "abstain"
    assert not result["release_complete"]


@pytest.mark.parametrize("outcome,value", [("observed_positive", "0.2"), ("tested_no_signal", "0"), ("inconclusive", "0.2"), ("untested", "0.2")])
def test_four_states_and_source_facts_survive_applicability_abstention(outcome, value):
    definition = profile_definition()
    editor = remake(one("editor_definition"), **definition["editor"], editing_window=definition["editing_window"])
    design = remake(one("design"), editor_definition_id=editor["record_id"])
    values = deepcopy(definition["context"])
    values["cell_state"] = {"value": values["cell_state"], "reason": None}
    context = remake(one("experimental_context"), **values)
    endpoint = remake(one("endpoint_definition"), required_context_ids=[context["record_id"]])
    observation = credited_observation()
    observation.update(design_id={"value": design["record_id"], "reason": None},
        context_id={"value": context["record_id"], "reason": None}, endpoint_id=endpoint["record_id"],
        normalized_value={"value": value, "reason": None}, adjudicated_outcome=outcome,
        adjudication_reason="Artificial explicit state", source_reported_outcome={"value": outcome, "reason": None})
    observation["detection_limit"] = {"value": {"value": 0.01, "reason": None}, "unit": "fraction"}
    before = deepcopy(observation)
    index = {r["record_id"]: r for r in [*ROWS, editor, design, context, endpoint]}
    args = dict(design_id=design["record_id"], context_id=context["record_id"], endpoint_id=endpoint["record_id"],
        donor_id=one("donor")["record_id"], sequence_configuration_id=one("sequence_configuration")["record_id"],
        mode="research_sample_planning", coverage_basis="callable_supported")
    policy = {**POLICY, "supported_domain": {"profile_id": PROFILE_ID, "context_id": context["record_id"], "endpoint_ids": [endpoint["record_id"]]}}
    result = assess_observation(observation, index, **args, policy=policy)["assessment"]
    assert not result["empirical_donor_credit"] and result["context"]["status"] == "unresolved"
    assert "domain_scientific_approval_pending" in result["reason_codes"]
    assert observation == before
    historical = assess_observation(observation, index, **args, policy=POLICY)["assessment"]
    assert historical["empirical_donor_credit"] == (outcome in {"observed_positive", "tested_no_signal"})
    assert result["quality"] == historical["quality"]


def test_profile_binding_is_a_stage_dependency(profile_case):
    path, config = profile_case
    study = resolve_study(path)
    original = effective_inputs(study, "search")
    historical = deepcopy(study)
    del historical["config"]["policies"]["evidence"]["supported_domain"]
    assert "supported_domain" not in effective_inputs(historical, "search")
    assert original != effective_inputs(historical, "search")
    changed = deepcopy(study)
    row(changed["config"], "experimental_context")["exposure"]["value"] = {"value": None, "reason": "artificial missing exposure"}
    changed["config"] = _reidentify(changed["config"])
    assert original != effective_inputs(changed, "search")


def test_artificial_profile_workflow_keeps_three_claims_separate(profile_case, tmp_path):
    path, config = profile_case
    result = invoke("validate", "--config", path)
    assert result["supported_domain"]["biological_applicability"]["status"] == "abstain"
    for command, parents in (("search", ()), ("observations", ()),
            ("evidence", ("search", "observations")), ("compare", ("search", "evidence")),
            ("questions", ("comparison",)), ("plan", ("questions", "evidence")),
            ("report", ("comparison", "questions", "plan"))):
        arguments = [command, "--config", path, "--output", tmp_path / command]
        for parent in parents:
            arguments += ["--" + parent, tmp_path / ("compare" if parent == "comparison" else parent)]
        invoke(*arguments)
    search = json.loads((tmp_path / "search/search-audit.json").read_text())
    report_audit = json.loads((tmp_path / "report/search-audit.json").read_text())
    domain = search["supported_domain"]
    assert domain == report_audit["supported_domain"]
    assert domain["coverage_completeness"]["declared_scope_execution"] == "complete"
    assert domain["coverage_completeness"]["status"] == "incomplete"
    assert domain["computational_support"]["status"] == "supported_for_artificial_checks"
    report = (tmp_path / "report/comparison-report.md").read_text()
    assert "Biological applicability: **abstain**" in report
    assert "Computational support: **supported_for_artificial_checks**" in report
    assert "Coverage completeness: **incomplete**" in report
    assert "domain_exposure_conflict_unresolved" in report
    rows = list(csv.DictReader(io.StringIO((tmp_path / "report/evidence-applicability.tsv").read_text()), delimiter="\t"))
    assert rows and all(r["current_empirical_donor_credit"] == "false" for r in rows)
    assert all("domain_scientific_approval_pending" in r["reason_codes"] for r in rows)
    assert all("domain_scientific_approval_pending" in r["current_gate_reason_codes"] for r in rows)
    retained = (tmp_path / "report/observations.jsonl").read_text()
    assert "9.8000000000000007" in retained and '"inconclusive"' in retained
    summaries = [json.loads(line) for line in (tmp_path / "compare/comparison-summaries.jsonl").read_text().splitlines()]
    decisions = [r for r in summaries if r["record_type"] == "comparison_decision"]
    assert decisions and all("domain_scientific_approval_pending" in r["reason_codes"] for r in decisions)
    assert all(r["state"] not in {"observed_difference_in_matched_context", "no_difference_detected_within_scope"} for r in decisions)


@pytest.mark.parametrize("field", ["exposure", "dose", "measurement_timepoint"])
def test_missing_exposure_quantities_cannot_be_promoted(profile_case, field):
    _, config = profile_case
    row(config, "experimental_context")[field]["value"] = {"value": None, "reason": "Artificial missing quantity"}
    config = _reidentify(config)
    validate_study(config)
    result = assess_domain(config)
    assert result["computational_support"]["status"] == "supported_for_artificial_checks"
    assert result["biological_applicability"]["status"] == "abstain"
    assert "domain_" + field + "_unresolved" in result["biological_applicability"]["reason_codes"]


def test_declared_callability_alone_cannot_close_missing_execution(profile_case):
    _, config = profile_case
    config["coverage_basis"] = "callable_supported"
    config["input_refs"].append({"role": "callability"})
    # Check the supplied declarations before execution.
    result = assess_domain(config)
    assert result["coverage_completeness"]["callability"] == "not_established"
    assert "domain_callability_not_established" in result["coverage_completeness"]["reason_codes"]
    assert not result["release_complete"]


def test_real_mode_cannot_activate_the_artificial_profile(profile_case):
    _, config = profile_case
    config["mode"], config["coverage_basis"] = "benchmark_simulation", "retained_callset_model"
    with pytest.raises(DomainAbstention) as exc:
        validate_study(config)
    assert "domain_real_processing_pending" in exc.value.assessment["request"]["reason_codes"]


def test_abstention_does_not_echo_rejected_private_contig_labels(profile_case, capsys):
    path, config = profile_case
    marker = "private-subject-contig-never-print"
    row(config, "domain")["intervals"][0]["contig"] = marker
    write(path, _reidentify(config))
    assert main(["validate", "--config", str(path)]) == 2
    captured = capsys.readouterr()
    assert marker not in captured.out + captured.err
    result = json.loads(captured.err)
    assert result["status"] == "abstained"
    assert "domain_nonprimary_contig_excluded" in result["supported_domain"]["request"]["reason_codes"]
    assert result["supported_domain"]["coverage_completeness"]["declared_contig_count"] == 2
