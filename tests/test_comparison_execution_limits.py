"""Search ceilings fail explicitly without changing scientific admission rules."""
from copy import deepcopy
from dataclasses import asdict
import json
from pathlib import Path

import pytest

from base_edit_safety.population import comparison, external_search
from base_edit_safety.population.adapters import crisprme_plus as adapter
from base_edit_safety.population.common import PopulationError, ResourceLimit, checksum
from base_edit_safety.population.comparison import (
    HostBudget, LocalLimits, host_limits_from_study, local_limits_from_study,
    local_search, normalize_search,
)
from base_edit_safety.population.comparison_contracts import validate_study
from base_edit_safety.population.comparison_inputs import prepare_inputs
from base_edit_safety.population.comparison_provenance import resolve_study, stage_key
from base_edit_safety.population.external_search import (
    NativeLimits, native_job_binding, native_limits_from_study, run_native,
    search_prepared, verify_native_job,
)
from test_comparison_external import FIXTURES, synthetic_inputs
from test_comparison_sequences import EXAMPLE, case


def configuration():
    return json.loads((EXAMPLE / "study.json").read_text())


@pytest.mark.parametrize("limits", [
    {"unknown": {}}, {"host": {"extra": 1}}, {"local": {"max_resident_bytes": 0}},
    {"host": {"timeout_seconds": True}}, {"local": {"timeout_seconds": 1.5}},
    {"native": {"cpus": 5}}, {"native": {"memory_bytes": 7 * 1024 ** 3}},
    {"native": {"disk_bytes": 1024}}, {"native": {"pids": 513}},
    {"native": {"max_rows": 0}}, {"native": {"timeout_seconds": "10"}},
])
def test_execution_configuration_rejects_unknown_or_unsafe_values(limits):
    config = configuration()
    config["policies"]["search"]["execution_limits"] = limits
    with pytest.raises(PopulationError):
        validate_study(config)


def test_optional_execution_limits_keep_defaults_and_bind_search_identity():
    study = resolve_study(EXAMPLE / "study.json")
    assert local_limits_from_study(study) == host_limits_from_study(study) == LocalLimits()
    assert native_limits_from_study(study) == NativeLimits()
    changed = deepcopy(study)
    changed["config"]["policies"]["search"]["execution_limits"] = {
        "local": {"timeout_seconds": 20}, "host": {"max_resident_bytes": 1024 ** 3},
        "native": {"max_rows": 7, "cpus": 1}}
    validate_study(changed["config"])
    assert local_limits_from_study(changed) == LocalLimits(timeout_seconds=20)
    assert host_limits_from_study(changed) == LocalLimits(max_resident_bytes=1024 ** 3)
    assert native_limits_from_study(changed) == NativeLimits(max_rows=7, cpus=1)
    assert stage_key(study, "search") != stage_key(changed, "search")
    assert stage_key(study, "observation_normalization") == stage_key(changed, "observation_normalization")


@pytest.mark.parametrize("role,stage", [("plan_history", "planning"),
                                       ("observation_history", "observation_normalization")])
def test_history_inputs_use_their_declared_stage(role, stage):
    config = configuration()
    config["input_refs"].append({"input_id": "historical_input", "role": role,
        "path": "historical.jsonl", "sha256": "a" * 64, "format": "jsonl", "stages": [stage]})
    validate_study(config)
    config["input_refs"][-1]["stages"] = ["search"]
    with pytest.raises(PopulationError, match="stage"):
        validate_study(config)


@pytest.mark.parametrize("kind", ["local", "host"])
def test_configured_memory_cap_fails_before_source_preparation(tmp_path, monkeypatch, kind):
    study = case(tmp_path, change_config=lambda config: config["policies"]["search"].update(
        execution_limits={kind: {"max_resident_bytes": 1}}))
    monkeypatch.setattr(comparison, "prepare_inputs", lambda *_: pytest.fail("Preparation ran after the memory cap"))
    with pytest.raises(ResourceLimit, match="memory ceiling") as error:
        local_search(study)
    audit = error.value.execution_audit
    assert audit["execution_status"] == "partial" and audit["retention_status"] == "incomplete"
    assert audit["reason_codes"] == [kind + "_memory_cap"]
    assert audit["effective_limits"][kind]["max_resident_bytes"] == 1
    assert audit["resource_usage"]["peak_resident_bytes"] > 1


def test_shared_host_deadline_covers_native_source_reconstruction(tmp_path, monkeypatch):
    prepared = prepare_inputs(resolve_study(EXAMPLE / "study.json"))
    original = prepared.iter_windows
    host = HostBudget(LocalLimits(timeout_seconds=1))
    def delayed_windows():
        for window in original():
            host.started -= 2
            yield window
    monkeypatch.setattr(prepared, "iter_windows", delayed_windows)
    monkeypatch.setattr(external_search, "run_native", lambda *_args, **_kwargs: pytest.fail("Native search ran after the host cap"))
    with pytest.raises(ResourceLimit) as error:
        search_prepared(prepared, output_dir=tmp_path, host_budget=host)
    assert error.value.execution_audit["phase"] == "native_source_reconstruction"
    assert error.value.execution_audit["reason_codes"] == ["host_timeout"]
    assert not list(tmp_path.iterdir())


def test_normalization_checks_hit_materialization_and_keeps_native_bindings():
    prepared = prepare_inputs(resolve_study(EXAMPLE / "study.json"))
    host = HostBudget(LocalLimits(timeout_seconds=1))
    binding = {"path": "/synthetic/evidence", "stage_key": "1" * 64,
               "result_sha256": "2" * 64, "completion_sha256": None}
    def delayed_hits():
        host.started -= 2
        yield {"window_id": "never_retained"}
    with pytest.raises(ResourceLimit) as error:
        normalize_search(prepared, delayed_hits(), [{"native_manifest": [binding]}], host_budget=host)
    audit = error.value.execution_audit
    assert audit["phase"] == "normalization_hit_materialization"
    assert audit["reason_codes"] == ["host_timeout"]
    assert audit["native_manifest"] == [binding]


def test_normalization_deadline_also_covers_source_window_replay(monkeypatch):
    prepared = prepare_inputs(resolve_study(EXAMPLE / "study.json"))
    original = prepared.iter_windows
    host = HostBudget(LocalLimits(timeout_seconds=1))
    def delayed_windows():
        for window in original():
            host.started -= 2
            yield window
    monkeypatch.setattr(prepared, "iter_windows", delayed_windows)
    with pytest.raises(ResourceLimit) as error:
        normalize_search(prepared, [], [], host_budget=host)
    assert error.value.execution_audit["phase"] == "normalization_source_reconstruction"


@pytest.mark.parametrize("operation", ["parse", "derive", "project"])
def test_host_cap_interrupts_native_alignment_processing_without_editing_raw_rows(operation):
    windows, design, editor, profile = synthetic_inputs("gapped")
    aliases = {f"chrW{i:08d}": window for i, window in enumerate(sorted(windows, key=lambda row: row["window_id"]))}
    path = FIXTURES / "gapped-ngg.targets.tsv"
    original_hash = checksum(path)
    parsed = adapter.parse_targets([path], windows=aliases, design=design, editor=editor, profile=profile)
    host = HostBudget(LocalLimits(timeout_seconds=1))
    calls = 0
    def check():
        nonlocal calls
        calls += 1
        if calls == 3:
            host.started -= 2
        host.check("native_" + operation)
    with pytest.raises(ResourceLimit) as error:
        if operation == "parse":
            adapter.parse_targets([path], windows=aliases, design=design, editor=editor,
                                  profile=profile, check_limits=check)
        elif operation == "derive":
            adapter.extend_leading_dna(parsed, windows=aliases, profile=profile, check_limits=check)
        else:
            _, projection = external_search._sequence_aliases(windows, 2)
            external_search._project_alignments(parsed, projection, profile, check)
    assert calls == 3 and error.value.execution_audit["reason_codes"] == ["host_timeout"]
    assert checksum(path) == original_hash
    assert parsed["retention_status"] == "complete"


def test_native_caps_bind_arguments_and_cache_without_changing_profile(tmp_path, monkeypatch):
    windows, design, editor, profile = synthetic_inputs("gapped")
    windows = windows[:1]
    windows[0]["source_genotypes"] = [{"gt": "0/0"}]
    profile["max_window_variants"] = 0
    original_profile = deepcopy(profile)
    monkeypatch.setattr(external_search.subprocess, "Popen", lambda *_args, **_kwargs: pytest.fail("No evaluable windows should launch Docker"))
    limits = NativeLimits(cpus=1, max_rows=1, timeout_seconds=9)
    first = run_native(windows, design, editor, profile, output_dir=tmp_path, limits=limits,
                       source_equivalence={"status": "not_applicable_synthetic"})
    path = Path(first["artifact_path"])
    payload = json.loads((path / "input.json").read_text())
    assert profile == original_profile == payload["profile"]
    assert payload["effective_row_cap"] == first["effective_row_cap"] == 1
    assert payload["limits"] == asdict(limits)
    assert first["effective_limits"]["host"] == asdict(LocalLimits())
    assert {"--pull=never", "--network=none", "--cpus=1", "--memory=4294967296"} <= set(first["argv"])
    command = payload["commands"][-1]
    assert command[command.index("--mm") + 1] == str(profile["max_mismatches"])
    assert command[command.index("--bDNA") + 1] == str(profile["max_dna_bulges"])
    assert command[command.index("--bRNA") + 1] == str(profile["max_rna_bulges"])
    assert command[command.index("--max-total-edits") + 1] == str(profile["max_total_edits"])
    assert not (path / "COMPLETE.json").exists()
    assert verify_native_job(native_job_binding(first))["execution_status"] == "partial"
    second = run_native(windows, design, editor, profile, output_dir=tmp_path,
                        limits=NativeLimits(max_rows=profile["max_rows"] + 1),
                        source_equivalence={"status": "not_applicable_synthetic"})
    assert second["effective_row_cap"] == profile["max_rows"]
    assert second["stage_key"] != first["stage_key"]


def test_native_attempt_interrupted_by_host_cap_retains_verifiable_partial_evidence(tmp_path, monkeypatch):
    windows, design, editor, profile = synthetic_inputs()
    windows = windows[:1]
    windows[0]["source_genotypes"] = [{"gt": "0/0"}]
    profile["max_window_variants"] = 0
    class ExpireOnSerialization(HostBudget):
        def check(self, phase=None):
            if phase == "native_input_serialization":
                self.started -= 2
            super().check(phase)
    host = ExpireOnSerialization(LocalLimits(timeout_seconds=1))
    monkeypatch.setattr(external_search.subprocess, "Popen", lambda *_args, **_kwargs: pytest.fail("Expired input serialization should never launch Docker"))
    result = run_native(windows, design, editor, profile, output_dir=tmp_path, host_budget=host,
                        source_equivalence={"status": "not_applicable_synthetic"})
    assert result["execution_status"] == "partial" and result["retention_status"] == "truncated"
    assert "host_timeout" in result["reason_codes"] and result["hits"] == []
    binding = native_job_binding(result)
    assert binding["completion_sha256"] is None and host.native_manifest == [binding]
    retained = verify_native_job(binding)
    assert retained["execution_status"] == "partial"
    assert "raw/container-links.json" in retained["artifact_manifest"]
    assert (Path(result["artifact_path"]) / "source-windows.jsonl").is_file()


def test_multi_pam_aggregation_cannot_bypass_configured_retention_cap(tmp_path, monkeypatch):
    study = case(tmp_path, change_config=lambda config: config["policies"]["search"].update(
        execution_limits={"native": {"max_rows": 1}}))
    prepared = prepare_inputs(study)
    prepared.profiles = deepcopy(prepared.profiles[:1])
    prepared.profiles[0]["profile"]["pam_patterns"] = ["AGG", "TGG"]
    def retained_job(windows, design, editor, profile, *, pam, **_kwargs):
        # This isolates aggregation. Native parser and admission are tested with
        # retained executable output in test_comparison_external.py.
        return {"hits": [{"window_id": windows[0]["window_id"], "pam_sequence": pam}],
            "execution_status": "complete", "retention_status": "complete", "exit_code": 0,
            "evaluated_window_ids": [window["window_id"] for window in windows],
            "skipped_windows": [], "reason_codes": [], "pam": pam}
    monkeypatch.setattr(external_search, "run_native", retained_job)
    monkeypatch.setattr(external_search, "native_job_binding", lambda job: {"test_pam": job["pam"]})
    result = search_prepared(prepared, output_dir=tmp_path / "native")
    assert len(result["jobs"]) == 2 and len(result["hits"]) == 1
    execution = result["execution"][0]
    assert execution["effective_row_cap"] == 1
    assert execution["execution_status"] == "partial" and execution["retention_status"] == "truncated"
    assert execution["reason_codes"] == ["native_profile_row_cap"]


def test_effective_native_row_cap_survives_canonical_scope_mapping():
    study = resolve_study(EXAMPLE / "study.json")
    result = local_search(study)
    execution = deepcopy(result["search-audit.json"]["execution"])
    for row in execution:
        row.update(effective_row_cap=1, execution_status="partial", retention_status="truncated",
                   reason_codes=["native_profile_row_cap"])
    normalized = normalize_search(prepare_inputs(study), [], execution)
    scopes = [row for row in normalized["comparison-records.jsonl"] if row["record_type"] == "comparison_search_scope"]
    assert scopes and all(row["backend_limits"]["max_rows"] == 1 for row in scopes)
    assert all(row["execution_status"] == "partial" and row["retention_status"] == "truncated" for row in scopes)
