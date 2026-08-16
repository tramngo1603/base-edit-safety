"""Load measurement tables, pipeline output, and hand-verified reference cases.

The evaluation runner accepts TSV, CSV, and JSON files. An expected row carries `dataset_id`,
`assay`, `cell_type`, `guide_id`,
and ``record_type``.  The loader also recognizes the normalized D1/D3/D4/D5
columns already emitted by this repository and current enumeration/design
record JSON envelopes.  It never downloads or executes a prediction tool.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
import json
from pathlib import Path
from typing import Any, Iterable, Mapping


@dataclass(frozen=True)
class LoadedInputs:
    measurements: tuple[dict[str, Any], ...]
    predictions: tuple[dict[str, Any], ...]
    golden: tuple[dict[str, Any], ...]
    dataset_metadata: dict[str, dict[str, Any]] = field(default_factory=dict)
    tool_versions: dict[str, str] = field(default_factory=dict)


def _clean_row(row: Mapping[str, Any]) -> dict[str, Any]:
    return {
        str(key): (value.strip() if isinstance(value, str) else value)
        for key, value in row.items()
        if key is not None
    }


def read_rows(path: str | Path) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Read a table/document and return rows plus its top-level metadata."""
    source = Path(path)
    suffix = source.suffix.lower()
    if suffix in {".tsv", ".csv"}:
        delimiter = "\t" if suffix == ".tsv" else ","
        with source.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter=delimiter)
            if reader.fieldnames is None:
                raise ValueError(f"input table has no header: {source}")
            return [_clean_row(row) for row in reader], {}
    if suffix != ".json":
        raise ValueError(f"unsupported evaluation input format: {source}")

    document = json.loads(source.read_text(encoding="utf-8"))
    if isinstance(document, list):
        return [_clean_row(row) for row in document], {}
    if not isinstance(document, dict):
        raise ValueError(f"JSON evaluation input must be an object or array: {source}")
    for key in ("records", "rows", "measurements", "golden_set"):
        if key in document:
            rows = document[key]
            if not isinstance(rows, list) or not all(isinstance(row, dict) for row in rows):
                raise ValueError(f"{source}: {key!r} must be an array of objects")
            metadata = {name: value for name, value in document.items() if name != key}
            return [_clean_row(row) for row in rows], metadata
    raise ValueError(f"JSON evaluation input has no records/rows array: {source}")


_DEFAULT_CONTEXT = {
    "D1": ("CHANGE-seq-BE", "not_reported"),
    "D2": ("Digenome-seq / targeted follow-up", "CD34+ HSPCs"),
    "D3": ("CHANGE-seq", "cell-free primary T-cell genomic DNA"),
    "D4": ("amplicon NGS", "not_reported"),
    "D5": ("targeted HBG amplicon sequencing", "primary CD34+ HSPCs"),
}


def _context(row: dict[str, Any]) -> dict[str, Any]:
    dataset = str(row.get("dataset_id", "")).strip()
    if not dataset:
        raise ValueError(f"evaluation row has no dataset_id: {row}")
    default_assay, default_cell = _DEFAULT_CONTEXT.get(dataset, ("not_reported", "not_reported"))
    return {
        **row,
        "dataset_id": dataset,
        "assay": row.get("assay") or default_assay,
        "cell_type": row.get("cell_type") or default_cell,
        "guide_id": row.get("guide_id") or row.get("target_id") or row.get("run_accession") or "",
    }


def _infer_measurement_type(row: dict[str, Any]) -> str:
    explicit = str(row.get("record_type", "")).strip()
    if explicit:
        return explicit
    if "detected_by_changeseq_be" in row or "changeseq_reads" in row:
        return "off_target_site"
    if "strand_normalized_c_to_t_frequency" in row:
        return "bystander"
    if "measured_editing_efficiency_percent" in row:
        return "on_target_outcome"
    raise ValueError(
        "measurement row has no record_type and does not match a supported ingested schema: "
        f"{row}"
    )


def _normalize_measurement(row: dict[str, Any]) -> dict[str, Any]:
    normalized = _context(row)
    normalized["record_type"] = _infer_measurement_type(normalized)
    if normalized["record_type"] == "off_target_site":
        if "measured_positive" not in normalized:
            if "detected_by_changeseq_be" in normalized:
                normalized["measured_positive"] = normalized["detected_by_changeseq_be"]
            elif "changeseq_reads" in normalized:
                normalized["measured_positive"] = normalized["changeseq_reads"]
            else:
                # Published detected-site-only tables contain positives by construction.
                normalized["measured_positive"] = True
    elif normalized["record_type"] == "bystander":
        normalized.setdefault(
            "measured_value", normalized.get("strand_normalized_c_to_t_frequency", "")
        )
    elif normalized["record_type"] == "on_target_outcome":
        normalized.setdefault(
            "measured_value", normalized.get("measured_editing_efficiency_percent", "")
        )
    return normalized


def load_measurements(paths: Iterable[str | Path]) -> tuple[list[dict[str, Any]], dict[str, dict[str, Any]]]:
    rows: list[dict[str, Any]] = []
    metadata: dict[str, dict[str, Any]] = {}
    for path in paths:
        source_rows, document_metadata = read_rows(path)
        rows.extend(_normalize_measurement(row) for row in source_rows)
        supplied = document_metadata.get("dataset_metadata", {})
        if isinstance(supplied, dict):
            for dataset_id, values in supplied.items():
                if isinstance(values, dict):
                    metadata.setdefault(str(dataset_id), {}).update(values)
    return rows, metadata


def _extract_tool_versions(metadata: Mapping[str, Any]) -> dict[str, str]:
    versions: dict[str, str] = {}
    direct = metadata.get("tool_versions", {})
    if isinstance(direct, dict):
        versions.update({str(name): str(value) for name, value in direct.items()})

    reproducibility = metadata.get("reproducibility")
    if not isinstance(reproducibility, dict):
        meta = metadata.get("meta")
        reproducibility = meta.get("reproducibility") if isinstance(meta, dict) else None
    if isinstance(reproducibility, dict):
        for tool in reproducibility.get("tools", []):
            if not isinstance(tool, dict) or not tool.get("name"):
                continue
            version = tool.get("version") or tool.get("commit") or "version_not_reported"
            versions[str(tool["name"])] = str(version)
    return versions


def _flag_from_design(record: Mapping[str, Any], needle: str) -> bool:
    flags = record.get("hbg_paralog_flags", {})
    return isinstance(flags, dict) and any(
        needle in str(value).lower() for value in flags.get("flags", [])
    )


def _design_record_rows(record: Mapping[str, Any]) -> list[dict[str, Any]]:
    candidate = record.get("candidate", {})
    if not isinstance(candidate, dict):
        raise ValueError("malformed design record: candidate must be an object")
    guide_id = candidate.get("guide_id") or record.get("provenance", {}).get("guide_id")
    editor_class = candidate.get("editor_class")
    base = {
        "dataset_id": record.get("dataset_id", ""),
        "assay": record.get("assay", ""),
        "cell_type": record.get("cell_type", ""),
        "guide_id": guide_id or "",
        "editor_class": editor_class or "",
    }
    ch = record.get("ch_annotation", {})
    ch_hits = ch.get("driver_offtargets", []) if isinstance(ch, dict) else []
    required = record.get("required_validation", [])
    rows = [
        {
            **base,
            "record_type": "guide_result",
            "status": "complete",
            "dual_engagement_flag": _flag_from_design(record, "also hits the paralog"),
            "intergenic_deletion_warning": any(
                isinstance(item, dict) and item.get("validates") == "hbg_paralog_rearrangement"
                for item in required
            ),
            "chip_driver_overlap_flag": bool(ch_hits),
            "evidence_tier": record.get("evidence_tier", ""),
        }
    ]
    spectrum = record.get("predicted_edit_spectrum", {})
    if isinstance(spectrum, dict):
        for bystander in spectrum.get("bystanders", []):
            if not isinstance(bystander, dict):
                continue
            allele = bystander.get("allele_key", {})
            rows.append(
                {
                    **base,
                    "record_type": "bystander_prediction",
                    "offset": allele.get("promoter_offset", bystander.get("promoter_offset", "")),
                    "ref": allele.get("promoter_sense_ref", ""),
                    "alt": allele.get("promoter_sense_alt", ""),
                    "predicted_value": bystander.get("predicted_rate"),
                }
            )
    return rows


def _enumeration_record_rows(records: list[Mapping[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    emitted_guides: set[str] = set()
    emitted_candidates: set[tuple[str, str, str, str]] = set()
    for record in records:
        provenance = record.get("provenance", {})
        provenance = provenance if isinstance(provenance, dict) else {}
        guide_id = str(provenance.get("guide_id", ""))
        base = {
            "dataset_id": record.get("dataset_id", ""),
            "assay": record.get("assay", ""),
            "cell_type": record.get("cell_type", ""),
            "guide_id": guide_id,
            "editor_class": record.get("editor_class", ""),
        }
        if record.get("on_target_or_bystander") == "bystander":
            rows.append(
                {
                    **base,
                    "record_type": "bystander_prediction",
                    "offset": -abs(int(record["promoter_offset"])),
                    "ref": record.get("promoter_sense_ref", ""),
                    "alt": record.get("promoter_sense_alt", ""),
                    "predicted_value": record.get("predicted_rate"),
                    "chrom": record.get("chrom", ""),
                    "pos_1based": record.get("pos_1based", ""),
                }
            )
        if guide_id and guide_id not in emitted_guides:
            emitted_guides.add(guide_id)
            ch_hits = provenance.get("ch_driver_offtargets", [])
            dual = bool(provenance.get("co_targets_paralog"))
            rows.append(
                {
                    **base,
                    "record_type": "guide_result",
                    "status": "complete",
                    "dual_engagement_flag": dual,
                    "intergenic_deletion_warning": dual,
                    "chip_driver_overlap_flag": bool(ch_hits),
                    "evidence_tier": "",
                }
            )
        for rank, site in enumerate(provenance.get("offtarget_sites_low_mismatch", []), start=1):
            if not isinstance(site, dict):
                continue
            unique = (guide_id, str(site.get("chrom")), str(site.get("start")), str(site.get("strand")))
            if unique in emitted_candidates:
                continue
            emitted_candidates.add(unique)
            rows.append(
                {
                    **base,
                    "record_type": "off_target_candidate",
                    "chrom": site.get("chrom", ""),
                    "start": site.get("start", ""),
                    "strand": site.get("strand", ""),
                    "rank": rank,
                    "candidate_list_complete": False,
                }
            )
    return rows


def _normalize_prediction_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    if rows and all(row.get("record_type") == "hbg_design_record" for row in rows):
        return [converted for row in rows for converted in _design_record_rows(row)]
    if rows and any("on_target_or_bystander" in row for row in rows):
        return _enumeration_record_rows(rows)
    normalized = []
    for row in rows:
        if not row.get("record_type"):
            raise ValueError(f"pipeline-output row has no record_type: {row}")
        # Empty dataset, assay, or cell fields receive context from the matching measurement or
        # reference-case row during the join.
        normalized.append(
            {
                **row,
                "dataset_id": row.get("dataset_id", ""),
                "assay": row.get("assay", ""),
                "cell_type": row.get("cell_type", ""),
                "guide_id": row.get("guide_id") or row.get("target_id") or "",
            }
        )
    return normalized


def load_pipeline_outputs(
    paths: Iterable[str | Path],
) -> tuple[list[dict[str, Any]], dict[str, dict[str, Any]], dict[str, str]]:
    rows: list[dict[str, Any]] = []
    dataset_metadata: dict[str, dict[str, Any]] = {}
    versions: dict[str, str] = {}
    for path in paths:
        source_rows, metadata = read_rows(path)
        rows.extend(_normalize_prediction_rows(source_rows))
        versions.update(_extract_tool_versions(metadata))
        supplied = metadata.get("dataset_metadata", {})
        if isinstance(supplied, dict):
            for dataset_id, values in supplied.items():
                if isinstance(values, dict):
                    dataset_metadata.setdefault(str(dataset_id), {}).update(values)
    return rows, dataset_metadata, versions


def load_golden_set(path: str | Path) -> list[dict[str, Any]]:
    rows, _metadata = read_rows(path)
    normalized = []
    for row in rows:
        item = _context(row)
        if not item.get("guide_id"):
            raise ValueError(f"golden-set row has no guide_id: {row}")
        normalized.append(item)
    return normalized


def load_all(
    measurement_paths: Iterable[str | Path],
    pipeline_output_paths: Iterable[str | Path],
    golden_set_path: str | Path,
) -> LoadedInputs:
    measurements, measurement_meta = load_measurements(measurement_paths)
    predictions, prediction_meta, versions = load_pipeline_outputs(pipeline_output_paths)
    for dataset_id, values in prediction_meta.items():
        measurement_meta.setdefault(dataset_id, {}).update(values)
    return LoadedInputs(
        measurements=tuple(measurements),
        predictions=tuple(predictions),
        golden=tuple(load_golden_set(golden_set_path)),
        dataset_metadata=measurement_meta,
        tool_versions=versions,
    )
