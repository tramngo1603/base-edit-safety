"""Check stage file hashes, protected output paths and installed commands."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import shutil

import pytest

from base_edit_safety.pipeline import population as pipeline
from base_edit_safety.pipeline.donor_panel import execute
from base_edit_safety.population.benchmarks import benchmark
from base_edit_safety.population.common import PopulationError, checksum, json_text, read_json, read_jsonl
from base_edit_safety.population.contracts import validate_config
from base_edit_safety.population.provenance import resolve_study, load_stage
from base_edit_safety.population.resolver import resolve_request

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def completed(tmp_path):
    source = tmp_path/"inputs"
    shutil.copytree(ROOT/"examples/population/synthetic", source)
    study = resolve_study(source/"study-config.json")
    root = tmp_path/"run"
    pipeline.build(study, root/"population")
    pipeline.search(study, root/"population", root/"search")
    pipeline.reconcile_stage(study, root/"search", root/"evidence")
    request = pipeline.actions_stage(study, root/"population", root/"search", root/"evidence", root/"panel-request.json")
    result = execute(root/"panel-request.json", root/"planning-manifest.json", "greedy_cost", root/"plan")
    return study, root, request, result


def test_synthetic_truth_hits_controls_plan_and_private_public_separation(completed):
    study, root, request, result = completed
    candidates = read_jsonl(root/"search/candidate-configurations.jsonl")
    assert len(candidates) == 4 and len(request["units"]) == 3
    assert sum(c["role"] == "intended_control" for c in candidates) == 1
    assert len(result["selected_actions"]) == 3 and result["total_cost"] == 4
    assert len(result["selected_donors"]) == 1 and result["objective_numerator"] == 3
    truth = read_json(ROOT/"tests/fixtures/population/expected.json")
    by_id = {c["configuration_id"]: c for c in candidates}
    hits = read_jsonl(root/"search/search-hits.jsonl")
    for donor, expected in truth["expected_donor_hits"].items():
        observed = {(by_id[h["configuration_id"]]["contig"], by_id[h["configuration_id"]]["strand"])
                    for h in hits if h["donor_id"] == donor and by_id[h["configuration_id"]]["role"] == "candidate_offtarget"}
        assert observed == {(r["contig"], r["strand"]) for r in expected}
    summary = benchmark(study["config"], request, study["run_id"], candidates=candidates,
                        variants=read_jsonl(root/"population/normalized-variants.jsonl"))
    pipeline.save_stage(study, "benchmark", root/"benchmark", {"benchmark-summary.json": summary})
    pipeline.report_stage(study, root, root/"population-report.md")
    text = (root/"population-report.md").read_text()
    assert "SYN_REF" not in text and "SYN_HET" not in text and str(root) not in text
    assert "proposed assays" in text and "Non-detection applies only to the matching assay" in text
    assert len(summary["rows"]) == 412
    assert {r["all_units"] for r in summary["rows"]} == {3}
    assert "Paired search scope" in text and "Retained variant strata" in text
    assert "Usable positive links: 0" in text and "Unmatched source measurements: 0" in text
    assert sum(summary["rare_variant_strata"]["retained_alt_alleles"].values()) == 3


def test_assay_review_table_contains_sequences_paths_and_evidence_requirements(completed):
    import csv
    study, root, request, result = completed
    with (root/"plan/assay-review.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == len(result["selected_actions"]) == 3
    assert {r["donor_id"] for r in rows} == set(result["selected_donors"])
    assert {r["contig"] for r in rows} == {"toy_gain", "toy_loss", "toy_reverse"}
    assert {r["strand"] for r in rows} == {"+", "-"}
    for row in rows:
        assert row["guide_spacer"] == "ACGTTGCAATGCTAGCTAAC"
        assert len(row["target_spacer"]) == 20 and row["target_pam"].endswith("GG")
        assert len(json.loads(row["target_path_json"])) == 23
        assert row["assay_id"] == study["config"]["assay_policy"]["assay_id"]
        assert row["prior_distinct_donors"] == "0" and row["remaining_donor_requirement_before_plan"] == "1"
        assert row["declared_feasibility"] == "declared_feasible"
        assert "distinct donors" in row["evidence_gap"] and "greedy_cost" in row["selection_rationale"]
        assert row["allele_resolution_required"] == row["locus_resolution_required"] == "True"
        for variant in json.loads(row["linked_variant_records_json"]):
            assert variant["selected_donor_call"]["call_status"] == "complete"
    assert (root/"plan/assay-review.tsv").stat().st_mode & 0o777 == 0o600


def test_variant_strata_use_called_subset_denominators_without_changing_plans(completed):
    from base_edit_safety.population.benchmarks import variant_strata
    study, _, request, _ = completed
    variants = [{"variant_id": "SYN_rare", "frequency": {"AC": [0, 1, 2, 2, 10], "AN": 400}},
                {"variant_id": "SYN_common", "frequency": {"AC": [10], "AN": 200}},
                {"variant_id": "SYN_unknown", "frequency": {"AC": [0], "AN": 0}}]
    candidates = [{"configuration_id": u["configuration_id"], "variant_links": links}
                  for u, links in zip(request["units"], [["SYN_rare", "SYN_common"], ["SYN_unknown"], []])]
    strata, units = variant_strata(variants, candidates, request["units"])
    assert strata["retained_alt_alleles"] == {"unknown_denominator": 1, "not_observed": 1, "singleton": 1,
        "rare_below_1_percent": 2, "frequency_1_to_5_percent": 1, "frequency_at_least_5_percent": 1}
    assert units[request["units"][2]["id"]] == {"no_linked_retained_variant"}
    config = deepcopy(study["config"])
    config["benchmark_policy"].update(selectors=["greedy_cost"], budget_grid=[config["planning_policy"]["max_cost"]])
    original = benchmark(config, request, "test:run")
    with_strata = benchmark(config, request, "test:run", variants=variants, candidates=candidates)
    assert with_strata["rows"][0]["selected_actions"] == original["rows"][0]["selected_actions"]
    assert all(r["all_units"] == 1 for r in with_strata["rows"][0]["variant_strata"])


def test_tampered_stages_and_changed_config_are_rejected(completed):
    study, root, _, _ = completed
    p = root/"search/donor-carriage.jsonl"
    p.write_text(p.read_text()+"\n")
    with pytest.raises(PopulationError, match="checksum"):
        load_stage(study, root/"search", "search")
    c = read_json(study["config_path"]); c["search_profile"]["max_mismatches"] = 1
    study["config_path"].write_text(json_text(c))
    with pytest.raises(PopulationError, match="Stale"):
        load_stage(resolve_study(study["config_path"]), root/"population", "build")


@pytest.mark.parametrize("mutation", ["eligibility", "invented_evidence", "possible_credit", "budget", "context"])
def test_cli_recomputes_content_even_with_updated_request_checksum(completed, mutation):
    study, root, request, _ = completed
    if mutation == "eligibility":
        next(a for a in request["actions"] if a["eligibility"] == "ineligible")["eligibility"] = "eligible"
    elif mutation == "invented_evidence":
        request["existing_usable_evidence"] = [{"unit_id": request["units"][0]["id"], "donor_id": "invented", "evidence_id": "stale"}]
    elif mutation == "possible_credit":
        request["actions"][0]["covered_units"] = [u["id"] for u in request["units"]]
    elif mutation == "budget":
        request["budget"]["max_cost"] += 1
    else:
        request["actions"][0]["resolution_contract_id"] = "other-locus"
    (root/"panel-request.json").write_text(json_text(request))
    manifest = read_json(root/"planning-manifest.json")
    manifest["request"]["sha256"] = checksum(root/"panel-request.json")
    (root/"planning-manifest.json").write_text(json_text(manifest))
    with pytest.raises(PopulationError):
        resolve_request(root/"panel-request.json", root/"planning-manifest.json")


def test_invalid_template_unknown_fields_and_alias_collision():
    c = read_json(ROOT/"examples/population/synthetic/study-config.json")
    for path, value in [("mode", "REPLACE_mode"), ("unexpected", True)]:
        x = deepcopy(c); x[path] = value
        with pytest.raises(PopulationError):
            validate_config(x)
    c["reference_policy"]["contig_alias_map"]["alias"] = "toy_gain"
    with pytest.raises(PopulationError, match="collision"):
        validate_config(c)


@pytest.mark.parametrize("selectors", [[], ["greedy_cost", "greedy_cost"]])
def test_empty_or_duplicate_selectors_fail_before_benchmark(selectors):
    config = read_json(ROOT/"examples/population/synthetic/study-config.json")
    config["benchmark_policy"]["selectors"] = selectors
    with pytest.raises(PopulationError):
        validate_config(config)
    with pytest.raises(PopulationError, match="nonempty, unique"):
        benchmark(config, read_json(ROOT/"examples/population/planner/request.json"), "test:run")


def test_explicit_overwrite_input_overlap_and_resource_cap_outputs(completed):
    study, root, request, _ = completed
    with pytest.raises(PopulationError, match="exists"):
        pipeline.build(study, root/"population")
    with pytest.raises(PopulationError, match="overlap"):
        pipeline.build(study, study["config_path"].parent, overwrite=True)
    c = read_json(study["config_path"]); c["resource_limits"]["max_configurations"] = 1
    study["config_path"].write_text(json_text(c))
    out = root/"capped"
    assert pipeline.main(["build", "--config", str(study["config_path"]), "--output-dir", str(out)]) == 4
    manifest = read_json(out/"run-manifest.json")
    assert manifest["execution_status"] == "partial"
    assert read_json(out/"execution-failure.json")["configuration_coverage"] == "unknown"
    assert not (out/"sequence-configurations.jsonl").exists()


def test_report_rejects_a_changed_plan(completed):
    study, root, request, result = completed
    pipeline.save_stage(study, "benchmark", root/"benchmark", {"benchmark-summary.json": benchmark(study["config"], request, study["run_id"])})
    result["total_cost"] = 0
    (root/"plan/panel-plan.json").write_text(json_text(result))
    with pytest.raises(PopulationError, match="provenance"):
        pipeline.report_stage(study, root, root/"summary.md")


def test_deterministic_semantic_outputs_from_second_execution(completed):
    study, root, request, _ = completed
    other = root.parent/"second"
    pipeline.build(study, other/"population")
    pipeline.search(study, other/"population", other/"search")
    for name in ("candidate-configurations.jsonl", "search-hits.jsonl", "donor-carriage.jsonl", "search-audit.json"):
        assert (root/"search"/name).read_bytes() == (other/"search"/name).read_bytes()


def test_record_validation_rejects_identity_and_unknown_field_changes(completed):
    from base_edit_safety.population.contracts import validate_record
    _, root, _, _ = completed
    row = read_jsonl(root/"search/candidate-configurations.jsonl")[0]
    broken = deepcopy(row); broken["patient_risk"] = 0
    with pytest.raises(PopulationError):
        validate_record(broken)
    row["configuration_id"] = "cfg:"+"0"*64
    with pytest.raises(PopulationError, match="identity"):
        validate_record(row)


def test_winner_only_external_output_never_becomes_all_configuration_truth():
    from base_edit_safety.population.external_contract import validate_external_contract
    c = {"record_type": "external_search_contract", "schema_version": "1", "backend_id": "synthetic_external",
         "executable_sha256": "0"*64, "source_revision": None, "native_schema": "synthetic", "argv": ["synthetic_backend"],
         "assembly": "synthetic", "reference_sha256": "1"*64, "phase_handling": "not_retained",
         "all_configuration_retention": False, "row_retention": "winner_only", "verified_fixture_ids": [],
         "stages": [{"name": "search", "status": "completed", "exit_code": 0}]}
    result = validate_external_contract(c)
    assert result["execution_status"] == "completed" and result["configuration_coverage"] == "unknown"
    c["stages"][0].update(status="failed", exit_code=1)
    assert validate_external_contract(c)["execution_status"] == "failed"


def test_all_units_denominator_retains_unavailable_configurations(completed):
    study, _, request, _ = completed
    for action in request["actions"]:
        action["eligibility"] = "requires_review"
    result = benchmark(study["config"], request, study["run_id"])
    assert all(r["all_units"] == 3 and r["actionable_units"] == 0 for r in result["rows"])
