"""Frozen selector requests preserve source, budget and control contracts."""
from copy import deepcopy
import os
from pathlib import Path
import stat

import pytest

from base_edit_safety.population.common import PopulationError, ResourceLimit, checksum, json_text, read_json
from base_edit_safety.population.comparison_benchmark import run_benchmark
from base_edit_safety.population.comparison_provenance import resolve_study, save_stage
from base_edit_safety.population.comparison_questions import prepare_question_set
from base_edit_safety.population.experiment_planning import SELECTORS, _plan_record
from planning_cases import full_fixture, optional_trap_case


def _write(path, value):
    path.write_text(json_text(value))


def _reference(path, base, *, directory=False):
    return {"path": os.path.relpath(path, base),
            "sha256": checksum(path / "run-manifest.json" if directory else path)}


def _request(root, *, budgets=(17, 18), selectors=SELECTORS, seeds=(7, 19), qset=None):
    root.mkdir(parents=True, exist_ok=True)
    qset = qset or full_fixture(root / "inputs")["question_set"]
    study = qset.study
    stages = {name: root / "stages" / name for name in (
        "search", "observation_normalization", "evidence_applicability", "comparison", "questions")}
    # The mock study's catalog already contains its sequences. Empty stages
    # supply the parent records needed to test planning.
    save_stage(study, "search", stages["search"], {"source-records.jsonl": []})
    save_stage(study, "observation_normalization", stages["observation_normalization"], {"observations.jsonl": []})
    save_stage(study, "evidence_applicability", stages["evidence_applicability"], {"evidence-records.jsonl": []},
               parents={name: stages[name] for name in ("search", "observation_normalization")})
    save_stage(study, "comparison", stages["comparison"], {"comparison-records.jsonl": []},
               parents={name: stages[name] for name in ("search", "evidence_applicability")})
    catalogs = {row["record_id"] for row in study["config"]["catalogs"]}
    save_stage(study, "questions", stages["questions"],
               {"question-records.jsonl": [row for row in qset.records if row["record_id"] not in catalogs]},
               parents={"comparison": stages["comparison"]})
    schedule = {"budgets": list(budgets), "selectors": list(selectors),
                "random_seeds": list(seeds) if "random_bundle" in selectors else []}
    frozen = {"schema_version": "1", "status": "frozen", "policy_version": "artificial-selector-fixture-v1",
              "code_snapshot": study["code"]["source_sha256"], "study_sha256": checksum(study["config_path"]),
              "selector_benchmark": schedule,
              "scientific_scope": "Artificial source-declared question fixture; no public-panel benchmark claim"}
    _write(root / "frozen.json", frozen)
    request = {"schema_version": "1", "policy_version": frozen["policy_version"],
               "study": _reference(study["config_path"], root),
               "questions": _reference(stages["questions"], root, directory=True),
               "evidence": _reference(stages["evidence_applicability"], root, directory=True),
               "frozen_protocol": _reference(root / "frozen.json", root), **deepcopy(schedule)}
    path = root / "request.json"
    _write(path, request)
    return path, qset, stages


def test_every_selector_pays18_and_no_selector_completes_at17(tmp_path):
    request, qset, _ = _request(tmp_path / "case")
    output = tmp_path / "results"
    result = run_benchmark(request, output)
    assert result["status"] == "completed" and result["result_count"] == 14
    rows = read_json(output / "benchmark-summary.json")
    for row in rows:
        complete = row["budget"] == 18
        assert row["unweighted_completion"] == row["weighted_completion"] == int(complete)
        assert row["total_cost"] == (18 if complete else 0)
        assert row["meets_request"] == complete
        assert row["unevaluable_question_fraction"] == 1
        assert row["evidence_state_counts"] == {"no_measurements": 1}
        assert row["exact_weighted_optimality_gap"] == (0 if complete else None)
        if complete:
            assert row["filled_catalog_slot_count"] == row["catalog_slot_count"] == 20
            assert row["engaged_intended_completeness"] == row["engaged_control_completeness"] == 1
    seeds = {row["random_seed"] for row in rows if row["selector"] == "random_bundle"}
    assert seeds == {7, 19}
    before = {path.name: path.read_bytes() for path in output.iterdir()}
    assert run_benchmark(request, output, reuse=True) == result
    assert before == {path.name: path.read_bytes() for path in output.iterdir()}
    assert stat.S_IMODE(output.stat().st_mode) == 0o700
    assert all(stat.S_IMODE(path.stat().st_mode) == 0o600 for path in output.iterdir())
    assert all(row["record_type"] != "comparison_observation_revision" for row in qset.records)


def test_mandatory_seed_remains_identical_against_a_high_weight_optional_trap(tmp_path):
    fixture = full_fixture(tmp_path / "source")
    qset = optional_trap_case(fixture, tmp_path / "optional-trap")
    request, _, _ = _request(tmp_path / "case", budgets=(18,), selectors=SELECTORS[:-1], qset=qset)
    output = tmp_path / "results"
    result = run_benchmark(request, output)
    plans = [read_json(output / item["file"]) for item in result["results"]]
    required = {qid for qid, row in qset.questions.items() if row["hard_required"]}
    assert len(required) == 1
    for plan in plans:
        assert plan["mandatory_seed"] == plans[0]["mandatory_seed"]
        assert plan["mandatory_seed"]["total_cost"] == 18
        assert plan["accounting"]["total_cost"] == 18
        assert set(plan["requirements"]["completed_question_ids"]) == required
        assert plan["weighted_completion"] == 1


@pytest.mark.parametrize("field,value,match", [
    ("budgets", [18], "frozen selector policy"),
    ("budgets", [17, True], "integer bounds"),
    ("selectors", ["complete_bundle_greedy", "complete_bundle_greedy"], "duplicate"),
    ("selectors", ["unregistered-selector"], "unsupported selector"),
    ("random_seeds", [19], "frozen selector policy"),
    ("policy_version", "post-outcome-policy", "frozen selector policy"),
    ("undeclared", True, "unexpected fields"),
])
def test_request_cannot_expand_or_replace_the_frozen_policy(tmp_path, field, value, match):
    path, _, _ = _request(tmp_path / "case")
    request = read_json(path)
    request[field] = value
    _write(path, request)
    with pytest.raises(PopulationError, match=match):
        run_benchmark(path, tmp_path / "results")
    assert not (tmp_path / "results").exists()


@pytest.mark.parametrize("field,value,match", [
    ("status", "draft", "status frozen"),
    ("code_snapshot", "0" * 64, "frozen code snapshot"),
    ("study_sha256", "0" * 64, "frozen study checksum"),
])
def test_frozen_protocol_must_bind_the_current_study_and_code(tmp_path, field, value, match):
    path, _, _ = _request(tmp_path / "case")
    frozen_path = path.parent / "frozen.json"
    frozen = read_json(frozen_path)
    frozen[field] = value
    _write(frozen_path, frozen)
    request = read_json(path)
    request["frozen_protocol"]["sha256"] = checksum(frozen_path)
    _write(path, request)
    with pytest.raises(PopulationError, match=match):
        run_benchmark(path, tmp_path / "results")


def test_later_evidence_cannot_be_substituted_for_the_question_ancestry(tmp_path):
    path, qset, stages = _request(tmp_path / "case", budgets=(18,), selectors=("canonical_bundle",))
    alternate = tmp_path / "later-evidence"
    save_stage(qset.study, "evidence_applicability", alternate,
               {"evidence-records.jsonl": [], "changed-assessment-note.txt": "A distinct retained assessment\n"},
               parents={name: stages[name] for name in ("search", "observation_normalization")})
    request = read_json(path)
    request["evidence"] = _reference(alternate, path.parent, directory=True)
    _write(path, request)
    with pytest.raises(PopulationError, match="evidence used to generate"):
        run_benchmark(path, tmp_path / "results")


def test_source_changes_and_rehashed_false_plan_are_rejected_on_reuse(tmp_path):
    path, qset, stages = _request(tmp_path / "case", budgets=(18,), selectors=("canonical_bundle",))
    output = tmp_path / "results"
    run_benchmark(path, output)
    source = Path(next(row["path"] for row in qset.study["sources"].values() if row["role"] == "variants"))
    original = source.read_bytes()
    source.write_bytes(original + b"\n")
    with pytest.raises(PopulationError, match="checksum changed"):
        run_benchmark(path, output, reuse=True)
    source.write_bytes(original)
    plan_path = output / "plan-0001.json"
    plan = read_json(plan_path)
    plan["accounting"]["total_cost"] = 0
    plan["plan_record"] = _plan_record(plan)
    plan["plan_id"] = plan["plan_record"]["plan_id"]
    _write(plan_path, plan)
    manifest_path = output / "run-manifest.json"
    manifest = read_json(manifest_path)
    manifest["outputs"][plan_path.name] = checksum(plan_path)
    manifest["results"][0]["plan_id"] = plan["plan_id"]
    _write(manifest_path, manifest)
    with pytest.raises(PopulationError, match="source-bound recomputation"):
        run_benchmark(path, output, reuse=True)


def test_changed_ancestor_or_result_bytes_never_reuse_as_a_complete_run(tmp_path):
    path, _, stages = _request(tmp_path / "case", budgets=(18,), selectors=("canonical_bundle",))
    output = tmp_path / "results"
    run_benchmark(path, output)
    summary = output / "benchmark-summary.tsv"
    summary.write_text(summary.read_text().replace("18", "19", 1))
    with pytest.raises(PopulationError, match="artifact changed"):
        run_benchmark(path, output, reuse=True)
    ancestor = stages["search"] / "source-records.jsonl"
    ancestor.write_text("\n")
    with pytest.raises(PopulationError, match="artifact changed"):
        run_benchmark(path, output, reuse=True)


def test_exact_size_limit_is_retained_without_an_optimality_claim(tmp_path):
    fixture = full_fixture(tmp_path / "source")
    study_path = fixture["study"]["config_path"]
    config = read_json(study_path)
    config["policies"]["planning"]["exact_action_cap"] = 1
    _write(study_path, config)
    qset = prepare_question_set(resolve_study(study_path), fixture["records"])
    path, _, _ = _request(tmp_path / "case", budgets=(18,), selectors=("exact_tiny",), qset=qset)
    output = tmp_path / "results"
    run_benchmark(path, output)
    row = read_json(output / "benchmark-summary.json")[0]
    assert not row["optimality_established"]
    assert row["exact_weighted_optimality_gap"] is None
    assert row["search_status_counts"] == {"exact_size_limit_exceeded": 1}
    assert row["status"] == "no_feasible_bundle_found"


def test_existing_and_overlapping_outputs_fail_before_selection(tmp_path, monkeypatch):
    path, _, stages = _request(tmp_path / "case", budgets=(18,), selectors=("canonical_bundle",))
    output = tmp_path / "exists"
    output.mkdir()
    def unexpected(*args, **kwargs):
        raise AssertionError("Selection must not start for an invalid output path")
    monkeypatch.setattr("base_edit_safety.population.comparison_benchmark._plan_prepared", unexpected)
    for target, match in ((output, "exists"), (stages["search"] / "child", "overlap"), (path.parent, "overlap")):
        with pytest.raises(PopulationError, match=match):
            run_benchmark(path, target)


def test_interruption_and_source_changes_during_execution_leave_no_completed_output(tmp_path, monkeypatch):
    import base_edit_safety.population.comparison_benchmark as benchmark
    path, _, stages = _request(tmp_path / "case", budgets=(17, 18), selectors=("canonical_bundle",))
    production = benchmark._plan_prepared
    calls = 0
    def interrupted(*args, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise ResourceLimit("Declared execution cap reached")
        return production(*args, **kwargs)
    monkeypatch.setattr(benchmark, "_plan_prepared", interrupted)
    output = tmp_path / "interrupted"
    with pytest.raises(ResourceLimit):
        run_benchmark(path, output)
    assert not output.exists()
    def changed(*args, **kwargs):
        plan = production(*args, **kwargs)
        (stages["search"] / "source-records.jsonl").write_text("\n")
        return plan
    monkeypatch.setattr(benchmark, "_plan_prepared", changed)
    output = tmp_path / "changed"
    with pytest.raises(PopulationError, match="input changed during execution"):
        run_benchmark(path, output)
    assert not output.exists()


def test_case_archive_keeps_completed_case_and_truthful_remaining_status(tmp_path, monkeypatch):
    import base_edit_safety.population.comparison_benchmark as benchmark
    path, _, _ = _request(tmp_path / "case", budgets=(0, 17, 18), selectors=("canonical_bundle",))
    original = benchmark._plan_prepared
    calls = 0

    def interrupted(*args, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise KeyboardInterrupt
        return original(*args, **kwargs)
    monkeypatch.setattr(benchmark, "_plan_prepared", interrupted)
    output, archive = tmp_path / "results", tmp_path / "cases"
    with pytest.raises(KeyboardInterrupt):
        run_benchmark(path, output, save_cases=archive)
    assert not output.exists()
    progress = read_json(archive / "progress.json")
    assert progress["status"] == "interrupted"
    assert [row["status"] for row in progress["cases"]] == ["completed", "incomplete", "not_started"]
    case = archive / "case-0001"
    manifest = read_json(case / "run-manifest.json")
    assert manifest["outputs"]["plan.json"] == checksum(case / "plan.json")
    assert manifest["budget"] == 0 and manifest["status"] == "completed"
    assert stat.S_IMODE(archive.stat().st_mode) == 0o700
    assert stat.S_IMODE((case / "plan.json").stat().st_mode) == 0o600
    before = {str(p.relative_to(archive)): p.read_bytes() for p in archive.rglob("*") if p.is_file()}
    with pytest.raises(PopulationError, match="Output exists"):
        run_benchmark(path, output, save_cases=archive)
    assert before == {str(p.relative_to(archive)): p.read_bytes() for p in archive.rglob("*") if p.is_file()}


def test_finished_case_archive_matches_the_complete_benchmark(tmp_path):
    path, _, _ = _request(tmp_path / "case", budgets=(17, 18), selectors=("canonical_bundle",))
    output, archive = tmp_path / "results", tmp_path / "cases"
    run_benchmark(path, output, save_cases=archive)
    progress = read_json(archive / "progress.json")
    assert progress["status"] == "schedule_completed"
    assert progress["completed_output_manifest_sha256"] == checksum(output / "run-manifest.json")
    for index in range(1, 3):
        assert (archive / f"case-{index:04d}" / "plan.json").read_bytes() == (output / f"plan-{index:04d}.json").read_bytes()


@pytest.mark.parametrize("target", ["output", "input", "existing"])
def test_case_archive_path_is_checked_before_any_selection(tmp_path, monkeypatch, target):
    import base_edit_safety.population.comparison_benchmark as benchmark
    path, _, stages = _request(tmp_path / "case", budgets=(18,), selectors=("canonical_bundle",))
    output = tmp_path / "results"
    existing = tmp_path / "existing"
    existing.mkdir()
    archive = {"output": output / "child", "input": stages["search"] / "child", "existing": existing}[target]

    def unexpected(*args, **kwargs):
        raise AssertionError("Selection must not start for an invalid archive path")
    monkeypatch.setattr(benchmark, "_plan_prepared", unexpected)
    with pytest.raises(PopulationError, match="overlap|exists"):
        run_benchmark(path, output, save_cases=archive)


@pytest.mark.parametrize("notes_available", [False, True])
def test_failed_progress_write_preserves_original_interruption(tmp_path, monkeypatch, notes_available):
    import base_edit_safety.population.comparison_benchmark as benchmark

    class Interrupted(KeyboardInterrupt):
        def __getattribute__(self, name):
            if name == "add_note" and not notes_available:
                raise AttributeError("Exception notes are unavailable")
            return super().__getattribute__(name)

    original = Interrupted("original interruption")
    path, _, _ = _request(tmp_path / "case", budgets=(18,), selectors=("canonical_bundle",))

    def interrupted(*args, **kwargs):
        raise original

    def failed_status(*args, **kwargs):
        raise OSError("progress disk unavailable")

    monkeypatch.setattr(benchmark, "_plan_prepared", interrupted)
    monkeypatch.setattr(benchmark._CaseOutput, "failed", failed_status)
    archive = tmp_path / "cases"
    with pytest.raises(Interrupted) as caught:
        run_benchmark(path, tmp_path / "results", save_cases=archive)
    assert caught.value is original
    assert original.case_archive_path == str(archive.resolve())
    assert (archive / "progress.json").is_file()
    if notes_available and hasattr(BaseException, "add_note"):
        assert original.__notes__ == ["Could not update case progress: progress disk unavailable"]
