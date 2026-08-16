"""Filesystem and safety-boundary tests for tuning and held-out entry points."""

import json
from pathlib import Path

import pytest

from base_edit_safety.evaluation import runner
from base_edit_safety.evaluation.runner import run_held_out, run_tuning
from base_edit_safety.evaluation.safety import SafetyViolation, verify_freeze_tag


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "tests" / "fixtures" / "evaluation"
REGISTRY = ROOT / "validation" / "dataset_registry.md"
PROTOCOL = ROOT / "validation" / "validation-protocol.md"
OVERLAP = ROOT / "validation" / "han_overlap_map.json"


def _clock():
    values = iter(("2026-08-16T10:00:00+00:00", "2026-08-16T10:00:01+00:00"))
    return lambda: next(values)


def _tuning_kwargs(tmp_path):
    return {
        "measurement_paths": [FIXTURES / "tuning_measurements.tsv"],
        "pipeline_output_paths": [FIXTURES / "tuning_pipeline.json"],
        "golden_set_path": FIXTURES / "tuning_golden.json",
        "output_dir": tmp_path / "tuning-output",
        "registry_path": REGISTRY,
        "protocol_path": PROTOCOL,
        "han_overlap_map_path": OVERLAP,
        "fixture_mode": True,
        "fixture_root": ROOT / "tests" / "fixtures",
        "repo_root": ROOT,
        "clock": _clock(),
    }


def _held_kwargs(tmp_path):
    return {
        "measurement_paths": [FIXTURES / "held_out_measurements.tsv"],
        "pipeline_output_paths": [FIXTURES / "held_out_pipeline.json"],
        "golden_set_path": FIXTURES / "held_out_golden.json",
        "output_dir": tmp_path / "held-output",
        "freeze_tag": "v0.0.0-validation-fixture",
        "registry_path": REGISTRY,
        "protocol_path": PROTOCOL,
        "han_overlap_map_path": OVERLAP,
        "fixture_mode": True,
        "fixture_root": ROOT / "tests" / "fixtures",
        "repo_root": ROOT,
        "clock": _clock(),
    }


def test_tuning_runner_writes_section_9_skeleton_numbers_and_manifest(tmp_path):
    artifacts = run_tuning(**_tuning_kwargs(tmp_path))
    assert artifacts.report_markdown.is_file()
    assert artifacts.metrics_json.is_file()
    assert artifacts.run_manifest_json.is_file()

    metrics = json.loads(artifacts.metrics_json.read_text())
    assert {row["metric_id"] for row in metrics["metrics"]} == {f"M{i}" for i in range(1, 10)}
    m1 = next(row for row in metrics["metrics"] if row["dataset_id"] == "D1" and row["metric_id"] == "M1")
    assert m1["value"] == 0.8 and m1["status"] == "pass"

    manifest = json.loads(artifacts.run_manifest_json.read_text())
    assert manifest["repository"]["commit_sha"]
    assert manifest["started_at_utc"] == "2026-08-16T10:00:00+00:00"
    assert manifest["ended_at_utc"] == "2026-08-16T10:00:01+00:00"
    assert manifest["tool_versions"]["fixture-predictor"] == "1.2.3"
    assert len(manifest["input_checksums"]) == 6
    assert all(len(item["sha256"]) == 64 for item in manifest["input_checksums"])

    report = artifacts.report_markdown.read_text()
    for heading in (
        "## 1. Summary", "## 2. Per-dataset", "## 3. Failure log",
        "## 4. Held-out results", "## 5. Limitations", "## 6. Appendix A",
        "## 7. Appendix B", "## 8. Appendix C", "## 9. Human sign-off",
    ):
        assert heading in report
    assert "tuning-only and contains no held-out concordance" in report


def test_tuning_entry_point_rejects_locked_held_out_role(tmp_path):
    kwargs = _tuning_kwargs(tmp_path)
    kwargs.update(
        measurement_paths=[FIXTURES / "held_out_measurements.tsv"],
        pipeline_output_paths=[FIXTURES / "held_out_pipeline.json"],
        golden_set_path=FIXTURES / "held_out_golden.json",
        clock=_clock(),
    )
    with pytest.raises(SafetyViolation, match="tuning entry point refuses"):
        run_tuning(**kwargs)


def test_real_path_is_blocked_while_protocol_and_registry_are_pre_freeze(tmp_path, monkeypatch):
    monkeypatch.delenv("CI", raising=False)
    monkeypatch.delenv("BASE_EDIT_SAFETY_FIXTURES_ONLY", raising=False)
    kwargs = _tuning_kwargs(tmp_path)
    kwargs["fixture_mode"] = False
    with pytest.raises(SafetyViolation, match="before protocol lock"):
        run_tuning(**kwargs)


def test_ci_refuses_non_fixture_evaluation_before_opening_inputs(tmp_path, monkeypatch):
    monkeypatch.setenv("CI", "true")
    kwargs = _tuning_kwargs(tmp_path)
    kwargs["fixture_mode"] = False
    with pytest.raises(SafetyViolation, match="CI may execute.*fixture mode"):
        run_tuning(**kwargs)


def test_fixture_mode_rejects_primary_input_outside_fixture_tree(tmp_path):
    outside = tmp_path / "not-a-fixture.tsv"
    outside.write_text("this file must never be opened as evaluation data\n")
    kwargs = _tuning_kwargs(tmp_path)
    kwargs["measurement_paths"] = [outside]
    with pytest.raises(SafetyViolation, match="non-fixture primary input"):
        run_tuning(**kwargs)


def test_held_out_is_one_shot_and_lock_carries_timestamp_and_input_hashes(tmp_path, monkeypatch):
    monkeypatch.setattr(runner, "verify_freeze_tag", lambda _root, _tag: "a" * 40)
    kwargs = _held_kwargs(tmp_path)
    artifacts = run_held_out(**kwargs)
    assert artifacts.held_out_lock is not None and artifacts.held_out_lock.is_file()
    lock_before = artifacts.held_out_lock.read_bytes()
    lock = json.loads(lock_before)
    assert lock["created_at_utc"] == "2026-08-16T10:00:00+00:00"
    assert lock["freeze_tag"] == "v0.0.0-validation-fixture"
    assert len(lock["input_checksums"]) == 6
    assert len(lock["combined_input_sha256"]) == 64

    report = artifacts.report_markdown.read_text()
    assert "single permitted held-out run" in report
    assert "Clean points" in report and "Partially entangled points" in report
    assert "| 2 | 2 | 6 |" in report

    kwargs["clock"] = _clock()
    with pytest.raises(SafetyViolation, match="one-shot.*already consumed"):
        run_held_out(**kwargs)
    assert artifacts.held_out_lock.read_bytes() == lock_before


def test_held_out_requires_an_existing_git_tag(tmp_path):
    with pytest.raises(SafetyViolation, match="missing/non-commit freeze tag"):
        verify_freeze_tag(ROOT, "definitely-no-such-validation-tag")
