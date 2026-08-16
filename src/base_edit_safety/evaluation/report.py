"""Render the protocol-§9 report skeleton and machine-readable artifacts."""

from __future__ import annotations

from dataclasses import dataclass
import json
import os
from pathlib import Path
from typing import Any, Mapping

from .manifest import atomic_write_json
from .models import EvaluationResult, MetricResult


@dataclass(frozen=True)
class ReportArtifacts:
    report_markdown: Path
    metrics_json: Path
    run_manifest_json: Path
    held_out_lock: Path | None = None


def _cell(value: Any) -> str:
    return str(value).replace("|", "\\|").replace("\n", " ")


def _number(value: float | None) -> str:
    if value is None:
        return "N/A"
    return f"{value:.6g}"


def _threshold(metric: MetricResult) -> str:
    return json.dumps(metric.threshold, sort_keys=True, separators=(",", ":"))


def _section(text: str, number: int, next_number: int) -> str:
    marker = f"## {number}."
    next_marker = f"\n## {next_number}."
    if marker not in text:
        return "Section unavailable in protocol input."
    return text.split(marker, 1)[1].split(next_marker, 1)[0].strip()


def render_report(
    result: EvaluationResult,
    manifest: Mapping[str, Any],
    *,
    mode: str,
    protocol_text: str,
    registry_text: str,
) -> str:
    lines = [
        "# Validation report — HBG1/HBG2 Base-Edit Safety Triage",
        "",
        f"Run mode: `{mode}`. Generated: `{manifest['ended_at_utc']}`.",
        "All values below were inserted from `metrics.json`; no result field is hand-edited.",
        "",
        "## 1. Summary: metrics table with pass/fail",
        "",
        "| Dataset | Assay | Cell type | Metric | Editor | Leakage stratum | Value | n/N | Threshold | Status |",
        "| --- | --- | --- | --- | --- | --- | ---: | ---: | --- | --- |",
    ]
    for metric in result.metrics:
        ratio = (
            f"{metric.numerator}/{metric.denominator}"
            if metric.numerator is not None and metric.denominator is not None
            else "N/A"
        )
        lines.append(
            "| " + " | ".join(
                _cell(value)
                for value in (
                    metric.dataset_id,
                    metric.assay,
                    metric.cell_type,
                    metric.metric_id,
                    metric.editor_class or "—",
                    metric.leakage_stratum or "all",
                    _number(metric.value),
                    ratio,
                    _threshold(metric),
                    metric.status,
                )
            ) + " |"
        )

    lines.extend([
        "",
        "## 2. Per-dataset concordance: plots + tables",
        "",
        "The numeric concordance tables are below. Plot generation is deliberately dependency-free in v1; "
        "the paired-point counts and coefficients are the authoritative output.",
    ])
    contexts = sorted({(m.dataset_id, m.assay, m.cell_type) for m in result.metrics})
    for dataset_id, assay, cell_type in contexts:
        lines.extend([
            "",
            f"### {_cell(dataset_id)} — {_cell(assay)} — {_cell(cell_type)}",
            "",
            "| Metric | Editor | Leakage stratum | Value | Numerator | Denominator | Components | Status | Reason |",
            "| --- | --- | --- | ---: | ---: | ---: | --- | --- | --- |",
        ])
        for metric in result.metrics:
            if (metric.dataset_id, metric.assay, metric.cell_type) != (dataset_id, assay, cell_type):
                continue
            lines.append(
                "| " + " | ".join(
                    _cell(value)
                    for value in (
                        metric.metric_id,
                        metric.editor_class or "—",
                        metric.leakage_stratum or "all",
                        _number(metric.value),
                        metric.numerator if metric.numerator is not None else "N/A",
                        metric.denominator if metric.denominator is not None else "N/A",
                        json.dumps(metric.components, sort_keys=True, separators=(",", ":")),
                        metric.status,
                        metric.reason or "",
                    )
                ) + " |"
            )

    lines.extend([
        "",
        "## 3. Failure log: every failure, F-code, and disposition",
        "",
    ])
    if result.failures:
        lines.extend([
            "| Dataset | Assay | Cell type | Metric | Editor | Leakage stratum | Subject | F-code | Disposition | Human note |",
            "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
        ])
        for failure in result.failures:
            lines.append(
                "| " + " | ".join(
                    _cell(value)
                    for value in (
                        failure.dataset_id,
                        failure.assay,
                        failure.cell_type,
                        failure.metric_id,
                        failure.editor_class or "—",
                        failure.leakage_stratum or "all",
                        failure.subject,
                        failure.failure_code,
                        failure.disposition,
                        failure.human_note,
                    )
                ) + " |"
            )
    else:
        lines.append("No misses or failed aggregate metrics were emitted.")

    lines.extend([
        "",
        "## 4. Held-out results",
        "",
    ])
    if mode == "held_out":
        lines.extend([
            "This was the single permitted held-out run. The one-shot lock was written before "
            "concordance inputs were loaded.",
            "",
            "| Dataset / assay / cell stratum | Clean points | Partially entangled points | Han-map keys |",
            "| --- | ---: | ---: | ---: |",
        ])
        if result.leakage_strata:
            for context, counts in sorted(result.leakage_strata.items()):
                lines.append(
                    f"| {_cell(context)} | {counts['clean']} | {counts['partially_entangled']} | "
                    f"{counts['overlap_map_key_count']} |"
                )
        else:
            lines.append("| No D5 strata were computable | 0 | 0 | 0 |")
        lines.extend([
            "",
            "`clean` means the position/editor-context key is absent from "
            "`validation/han_overlap_map.json`; `partially_entangled` means it is present.",
        ])
    else:
        lines.append("Not run. This report is tuning-only and contains no held-out concordance.")

    lines.extend([
        "",
        "## 5. Limitations (from protocol §8)",
        "",
        _section(protocol_text, 8, 9),
        "",
        "## 6. Appendix A: dataset registry",
        "",
        registry_text.strip(),
        "",
        "## 7. Appendix B: provenance dump",
        "",
        "| Input role | Path | SHA-256 | Bytes |",
        "| --- | --- | --- | ---: |",
    ])
    for item in manifest["input_checksums"]:
        lines.append(
            f"| {_cell(item['role'])} | `{_cell(item['path'])}` | `{item['sha256']}` | {item['size_bytes']} |"
        )
    lines.extend([
        "",
        f"Combined input SHA-256: `{manifest['combined_input_sha256']}`.",
        "",
        "## 8. Appendix C: environment",
        "",
        f"Repository commit: `{manifest['repository']['commit_sha']}`; "
        f"working tree dirty: `{manifest['repository']['working_tree_dirty']}`.",
        "",
        "| Tool/runtime | Version |",
        "| --- | --- |",
    ])
    for name, version in sorted(manifest["tool_versions"].items()):
        lines.append(f"| {_cell(name)} | `{_cell(version)}` |")
    if manifest.get("freeze_tag"):
        lines.extend([
            "",
            f"Freeze tag: `{manifest['freeze_tag']}` at `{manifest['freeze_tag_commit']}`.",
        ])
    lines.extend([
        "",
        "## 9. Human sign-off",
        "",
        "- [ ] I hand-verified every golden-set expected value against primary sources.",
        "- [ ] I verified the tuning/held-out split was respected.",
        "- [ ] I reviewed every F1–F4/F7 fix and its re-run.",
        "- [ ] I reviewed the held-out run manifest (single run, correct commit).",
        "- [ ] The numbers in §4 are the numbers produced by the pipeline, unedited.",
        "",
        "Signed: ____________________ — __________",
        "",
    ])
    return "\n".join(lines)


def _atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        temporary.write_text(text, encoding="utf-8")
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def write_report_artifacts(
    output_dir: str | Path,
    result: EvaluationResult,
    manifest: Mapping[str, Any],
    *,
    mode: str,
    protocol_text: str,
    registry_text: str,
    held_out_lock: str | Path | None = None,
) -> ReportArtifacts:
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    metrics_path = destination / "metrics.json"
    manifest_path = destination / "run_manifest.json"
    report_path = destination / "validation_report.md"
    atomic_write_json(metrics_path, result.as_dict())
    atomic_write_json(manifest_path, manifest)
    _atomic_write_text(
        report_path,
        render_report(
            result, manifest, mode=mode, protocol_text=protocol_text, registry_text=registry_text
        ),
    )
    return ReportArtifacts(
        report_markdown=report_path,
        metrics_json=metrics_path,
        run_manifest_json=manifest_path,
        held_out_lock=Path(held_out_lock) if held_out_lock is not None else None,
    )
