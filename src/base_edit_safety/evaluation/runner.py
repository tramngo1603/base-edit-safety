"""Run tuning and held-out evaluations after required pre-run checks."""

from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path
import re
from typing import Any, Callable, Iterable, Mapping

from .inputs import LoadedInputs, load_all
from .manifest import build_manifest, checksum_inputs, create_one_shot_lock, git_state
from .metrics import evaluate
from .report import ReportArtifacts, write_report_artifacts
from .safety import (
    SafetyViolation,
    assert_ci_fixture_mode,
    assert_dataset_roles,
    assert_fixture_inputs,
    assert_protocol_frozen,
    assert_registry_qc_admitted,
    parse_registry,
    verify_freeze_tag,
)


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_REGISTRY = REPO_ROOT / "validation" / "dataset_registry.md"
DEFAULT_PROTOCOL = REPO_ROOT / "validation" / "validation-protocol.md"
DEFAULT_HAN_OVERLAP_MAP = REPO_ROOT / "validation" / "han_overlap_map.json"
DEFAULT_FIXTURE_ROOT = REPO_ROOT / "tests" / "fixtures"


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _named_inputs(
    measurement_paths: Iterable[str | Path],
    pipeline_output_paths: Iterable[str | Path],
    golden_set_path: str | Path,
    registry_path: str | Path,
    protocol_path: str | Path,
    han_overlap_map_path: str | Path,
) -> list[tuple[str, str | Path]]:
    named = [("ingested_dataset_table", path) for path in measurement_paths]
    named.extend(("tool_pipeline_output", path) for path in pipeline_output_paths)
    named.extend(
        [
            ("golden_set", golden_set_path),
            ("locked_dataset_registry", registry_path),
            ("validation_protocol", protocol_path),
            ("han_overlap_map", han_overlap_map_path),
        ]
    )
    return named


def _registered_dataset_ids(inputs: LoadedInputs) -> set[str]:
    return {
        dataset_id
        for row in (*inputs.measurements, *inputs.predictions, *inputs.golden)
        if re.fullmatch(r"D\d+", dataset_id := str(row.get("dataset_id", "")).strip())
    }


def _prepare(
    *,
    fixture_mode: bool,
    fixture_root: str | Path,
    primary_paths: Iterable[str | Path],
    registry_path: str | Path,
    protocol_path: str | Path,
) -> dict[str, Any]:
    assert_ci_fixture_mode(fixture_mode=fixture_mode)
    if fixture_mode:
        assert_fixture_inputs(primary_paths, fixture_root)
    registry = parse_registry(registry_path)
    if not fixture_mode:
        # Complete protocol and registry checks before opening measurement or prediction tables.
        assert_protocol_frozen(protocol_path)
        assert_registry_qc_admitted(registry)
    return registry


def _assert_locked_control_paths(
    *,
    fixture_mode: bool,
    repo_root: str | Path,
    registry_path: str | Path,
    protocol_path: str | Path,
    han_overlap_map_path: str | Path,
) -> None:
    if fixture_mode:
        return
    root = Path(repo_root).resolve(strict=True)
    expected = {
        "registry": root / "validation" / "dataset_registry.md",
        "protocol": root / "validation" / "validation-protocol.md",
        "Han overlap map": root / "validation" / "han_overlap_map.json",
    }
    supplied = {
        "registry": Path(registry_path),
        "protocol": Path(protocol_path),
        "Han overlap map": Path(han_overlap_map_path),
    }
    wrong = {
        name: str(path)
        for name, path in supplied.items()
        if path.resolve(strict=True) != expected[name].resolve(strict=True)
    }
    if wrong:
        raise SafetyViolation(
            f"non-fixture evaluation refuses substituted locked control files: {wrong}"
        )


def _load_overlap(path: str | Path) -> dict[str, Any]:
    value = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError("Han overlap map must be a JSON object")
    return value


def run_tuning(
    *,
    measurement_paths: Iterable[str | Path],
    pipeline_output_paths: Iterable[str | Path],
    golden_set_path: str | Path,
    output_dir: str | Path,
    registry_path: str | Path = DEFAULT_REGISTRY,
    protocol_path: str | Path = DEFAULT_PROTOCOL,
    han_overlap_map_path: str | Path = DEFAULT_HAN_OVERLAP_MAP,
    fixture_mode: bool = False,
    fixture_root: str | Path = DEFAULT_FIXTURE_ROOT,
    repo_root: str | Path = REPO_ROOT,
    tool_versions: Mapping[str, Any] | None = None,
    clock: Callable[[], str] = _now,
) -> ReportArtifacts:
    """Run tuning metrics; every registered dataset on this path must have role ``tuning``."""
    measurements = tuple(measurement_paths)
    predictions = tuple(pipeline_output_paths)
    primary = (*measurements, *predictions, golden_set_path)
    _assert_locked_control_paths(
        fixture_mode=fixture_mode,
        repo_root=repo_root,
        registry_path=registry_path,
        protocol_path=protocol_path,
        han_overlap_map_path=han_overlap_map_path,
    )
    registry = _prepare(
        fixture_mode=fixture_mode,
        fixture_root=fixture_root,
        primary_paths=primary,
        registry_path=registry_path,
        protocol_path=protocol_path,
    )
    named = _named_inputs(
        measurements, predictions, golden_set_path, registry_path, protocol_path,
        han_overlap_map_path,
    )
    checksums = checksum_inputs(named)
    started = clock()
    inputs = load_all(measurements, predictions, golden_set_path)
    dataset_ids = _registered_dataset_ids(inputs)
    assert_dataset_roles(dataset_ids, registry, expected_role="tuning")
    if not fixture_mode:
        expected = {dataset_id for dataset_id, entry in registry.items() if entry.role == "tuning"}
        if dataset_ids != expected:
            raise SafetyViolation(
                f"non-fixture tuning must cover the complete locked tuning set {sorted(expected)}; "
                f"got {sorted(dataset_ids)}"
            )

    result = evaluate(
        inputs,
        han_overlap_map=_load_overlap(han_overlap_map_path),
        enforce_audit_sample_size=True,
    )
    ended = clock()
    versions = dict(inputs.tool_versions)
    versions.update(tool_versions or {})
    manifest = build_manifest(
        mode="tuning",
        started_at_utc=started,
        ended_at_utc=ended,
        repo_root=repo_root,
        input_checksums=checksums,
        tool_versions=versions,
        fixture_mode=fixture_mode,
    )
    return write_report_artifacts(
        output_dir,
        result,
        manifest,
        mode="tuning",
        protocol_text=Path(protocol_path).read_text(encoding="utf-8"),
        registry_text=Path(registry_path).read_text(encoding="utf-8"),
    )


def run_held_out(
    *,
    measurement_paths: Iterable[str | Path],
    pipeline_output_paths: Iterable[str | Path],
    golden_set_path: str | Path,
    output_dir: str | Path,
    freeze_tag: str,
    registry_path: str | Path = DEFAULT_REGISTRY,
    protocol_path: str | Path = DEFAULT_PROTOCOL,
    han_overlap_map_path: str | Path = DEFAULT_HAN_OVERLAP_MAP,
    lock_file: str | Path | None = None,
    fixture_mode: bool = False,
    fixture_root: str | Path = DEFAULT_FIXTURE_ROOT,
    repo_root: str | Path = REPO_ROOT,
    tool_versions: Mapping[str, Any] | None = None,
    clock: Callable[[], str] = _now,
) -> ReportArtifacts:
    """Consume and run the protocol-v1 held-out set exactly once."""
    measurements = tuple(measurement_paths)
    predictions = tuple(pipeline_output_paths)
    primary = (*measurements, *predictions, golden_set_path)
    _assert_locked_control_paths(
        fixture_mode=fixture_mode,
        repo_root=repo_root,
        registry_path=registry_path,
        protocol_path=protocol_path,
        han_overlap_map_path=han_overlap_map_path,
    )
    registry = _prepare(
        fixture_mode=fixture_mode,
        fixture_root=fixture_root,
        primary_paths=primary,
        registry_path=registry_path,
        protocol_path=protocol_path,
    )
    held_out_ids = {dataset_id for dataset_id, entry in registry.items() if entry.role == "held-out"}
    assert_dataset_roles(held_out_ids, registry, expected_role="held-out")
    freeze_commit = verify_freeze_tag(repo_root, freeze_tag)
    if not fixture_mode:
        repository = git_state(repo_root)
        if repository["commit_sha"] != freeze_commit:
            raise SafetyViolation(
                "held-out run must execute the exact freeze-tag commit: "
                f"HEAD={repository['commit_sha']}, tag={freeze_commit}"
            )
        if repository["working_tree_dirty"] is not False:
            raise SafetyViolation("held-out run refuses a dirty or unreadable working tree")
    lock_path = Path(lock_file) if lock_file is not None else Path(output_dir) / "held_out_one_shot.lock.json"
    if lock_path.exists():
        raise SafetyViolation(
            f"held-out evaluation is one-shot and was already consumed; lock exists: {lock_path}"
        )
    named = _named_inputs(
        measurements, predictions, golden_set_path, registry_path, protocol_path,
        han_overlap_map_path,
    )
    checksums = checksum_inputs(named)
    started = clock()

    # Create the non-overwritable run record before loading prediction-measurement comparison tables.
    create_one_shot_lock(
        lock_path,
        timestamp_utc=started,
        freeze_tag=freeze_tag,
        freeze_tag_commit=freeze_commit,
        input_checksums=checksums,
    )

    inputs = load_all(measurements, predictions, golden_set_path)
    dataset_ids = _registered_dataset_ids(inputs)
    assert_dataset_roles(dataset_ids, registry, expected_role="held-out")
    if dataset_ids != held_out_ids:
        raise SafetyViolation(
            f"held-out one-shot must cover exactly {sorted(held_out_ids)}; got {sorted(dataset_ids)}"
        )
    overlap = _load_overlap(han_overlap_map_path)
    result = evaluate(
        inputs,
        han_overlap_map=overlap,
        enforce_audit_sample_size=True,
    )
    # Report both protocol-§9 source-overlap groups, including a group with zero observations.
    if not result.leakage_strata:
        raise SafetyViolation("held-out result has no D5 Han leakage-stratum report")
    for values in result.leakage_strata.values():
        if not {"clean", "partially_entangled"} <= set(values):
            raise SafetyViolation("held-out result omitted a required Han leakage stratum")

    ended = clock()
    versions = dict(inputs.tool_versions)
    versions.update(tool_versions or {})
    manifest = build_manifest(
        mode="held_out",
        started_at_utc=started,
        ended_at_utc=ended,
        repo_root=repo_root,
        input_checksums=checksums,
        tool_versions=versions,
        fixture_mode=fixture_mode,
        freeze_tag=freeze_tag,
        freeze_tag_commit=freeze_commit,
    )
    return write_report_artifacts(
        output_dir,
        result,
        manifest,
        mode="held_out",
        protocol_text=Path(protocol_path).read_text(encoding="utf-8"),
        registry_text=Path(registry_path).read_text(encoding="utf-8"),
        held_out_lock=lock_path,
    )
