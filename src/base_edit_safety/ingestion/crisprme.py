"""Check and read supported CRISPRme+ integrated-result files.

CRISPRme+ runs separately. This reader accepts its recorded native TSV format and an older
publication workbook used for validation. Both use the same row checks and discard sample columns.
"""
from __future__ import annotations

from collections import Counter
from copy import deepcopy
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
import csv
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Iterable, Iterator

from .xlsx import iter_rows as iter_xlsx_rows
from .xlsx import sheet_names


CRISPRME_VERSION = "2.3.2"
CRISPRME_SCORING_MODEL = "CFD"
CRISPRME_TSV_FORMAT_SCHEMA = "crisprme-2.3.2-integrated-results-v1"
CRISPRME_XLSX_FORMAT_SCHEMA = "crisprme-publication-sg1617-integrated-results-v1"
SUPPORTED_FORMAT_SCHEMAS = {
    CRISPRME_TSV_FORMAT_SCHEMA,
    CRISPRME_XLSX_FORMAT_SCHEMA,
}

REQUIRED_COLUMNS = frozenset(
    {
        "Spacer+PAM",
        "Chromosome",
        "Start_coordinate_(highest_CFD)",
        "Strand_(highest_CFD)",
        "Aligned_spacer+PAM_(highest_CFD)",
        "Aligned_protospacer+PAM_REF_(highest_CFD)",
        "Aligned_protospacer+PAM_ALT_(highest_CFD)",
        "PAM_(highest_CFD)",
        "Mismatches_(highest_CFD)",
        "Bulges_(highest_CFD)",
        "Bulge_type_(highest_CFD)",
        "REF/ALT_origin_(highest_CFD)",
        "PAM_creation_(highest_CFD)",
        "CFD_score_(highest_CFD)",
        "CFD_score_REF_(highest_CFD)",
        "CFD_score_ALT_(highest_CFD)",
        "Variant_info_spacer+PAM_(highest_CFD)",
        "Variant_info_genome_(highest_CFD)",
        "Variant_MAF_(highest_CFD)",
        "Variant_rsID_(highest_CFD)",
        "Not_found_in_REF",
    }
)

SKIP_REASONS = frozenset(
    {
        "invalid_coordinate",
        "invalid_strand",
        "invalid_alignment",
        "invalid_numeric_value",
        "guide_mismatch",
        "unsupported_bulge_type",
        "inconsistent_span",
    }
)

_DECLARATION_FIELDS = {
    "record_type",
    "schema_version",
    "format_schema",
    "native_input",
    "selected_sheet",
    "guide_id",
    "guide_protospacer",
    "genome_build",
    "coordinate_system",
    "crisprme_version",
    "scoring_model",
    "variant_databases",
    "pam_patterns",
    "max_mismatches",
    "max_dna_bulge",
    "max_rna_bulge",
    "retention",
    "source_command",
}
_RETENTION_FIELDS = {
    "completeness",
    "top_n",
    "score_cutoff",
    "merge_window_bp",
    "all_candidate_coordinates_retained",
}
_IUPAC = frozenset("ACGTRYSWKMBDHVN")
_DNA_ALIGNMENT = frozenset("ACGTN-")
_PROHIBITED_COLUMN_PREFIX = "variant_samples"
_IUPAC_BASES = {
    "A": frozenset("A"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "T": frozenset("T"),
    "R": frozenset("AG"),
    "Y": frozenset("CT"),
    "S": frozenset("GC"),
    "W": frozenset("AT"),
    "K": frozenset("GT"),
    "M": frozenset("AC"),
    "B": frozenset("CGT"),
    "D": frozenset("AGT"),
    "H": frozenset("ACT"),
    "V": frozenset("ACG"),
    "N": frozenset("ACGT"),
}


class CrisprmeImportError(ValueError):
    """Raised when a CRISPRme input cannot be interpreted without guessing."""


class _RowFailure(ValueError):
    def __init__(self, reason: str):
        if reason not in SKIP_REASONS:  # pragma: no cover - internal invariant
            raise AssertionError(f"unknown CRISPRme row failure reason: {reason}")
        self.reason = reason
        super().__init__(reason)


@dataclass(frozen=True)
class CrisprmeImportResult:
    """Checked rows without sample fields, plus aggregate import counts."""

    rows: tuple[dict, ...]
    declaration: dict
    import_summary: dict
    sheets_ignored: dict
    privacy_filtering_applied: bool
    original_declaration: dict | None = None


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        handle = path.open("rb")
    except OSError as exc:
        raise CrisprmeImportError(f"cannot read CRISPRme input {path.name}") from exc
    with handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _required_text(value: object, field: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise CrisprmeImportError(f"CRISPRme declaration field {field} must be non-empty text")
    return value.strip()


def _nonnegative_int(value: object, field: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise CrisprmeImportError(
            f"CRISPRme declaration field {field} must be a non-negative integer"
        )
    return value


def _exact_object(value: object, fields: set[str], label: str) -> dict:
    if not isinstance(value, dict) or set(value) != fields:
        raise CrisprmeImportError(
            f"CRISPRme declaration {label} must contain exactly: {', '.join(sorted(fields))}"
        )
    return value


def _validate_score_cutoff(value: object) -> dict | None:
    if value is None:
        return None
    block = _exact_object(value, {"score_name", "minimum"}, "retention.score_cutoff")
    _required_text(block["score_name"], "retention.score_cutoff.score_name")
    minimum = block["minimum"]
    if isinstance(minimum, bool) or not isinstance(minimum, (int, float)):
        raise CrisprmeImportError(
            "CRISPRme declaration field retention.score_cutoff.minimum must be finite"
        )
    if not math.isfinite(float(minimum)):
        raise CrisprmeImportError(
            "CRISPRme declaration field retention.score_cutoff.minimum must be finite"
        )
    return block


def validate_import_declaration(
    declaration: object, native_path: Path | None = None
) -> dict:
    """Check the recorded settings and, when supplied, the native file checksum."""
    payload = deepcopy(_exact_object(declaration, _DECLARATION_FIELDS, "root"))
    if payload["record_type"] != "crisprme_import_declaration" or payload["schema_version"] != "1":
        raise CrisprmeImportError("unsupported CRISPRme import declaration record or schema")

    format_schema = payload["format_schema"]
    if format_schema not in SUPPORTED_FORMAT_SCHEMAS:
        raise CrisprmeImportError("unsupported CRISPRme format_schema")

    native = _exact_object(payload["native_input"], {"file", "sha256"}, "native_input")
    filename = _required_text(native["file"], "native_input.file")
    if Path(filename).name != filename:
        raise CrisprmeImportError(
            "CRISPRme declaration native_input.file must be a basename"
        )
    if native_path is not None and filename != native_path.name:
        raise CrisprmeImportError(
            "CRISPRme declaration native_input.file must be the supplied input basename"
        )
    checksum = native["sha256"]
    if not isinstance(checksum, str) or not re.fullmatch(r"[0-9a-fA-F]{64}", checksum):
        raise CrisprmeImportError("CRISPRme declaration native_input.sha256 is invalid")
    normalized_checksum = checksum.lower()
    if native_path is not None:
        actual_checksum = _sha256(native_path)
        if normalized_checksum != actual_checksum:
            raise CrisprmeImportError("CRISPRme input checksum does not match its declaration")
        normalized_checksum = actual_checksum
    native["sha256"] = normalized_checksum

    guide_id = _required_text(payload["guide_id"], "guide_id")
    protospacer = payload["guide_protospacer"]
    if (
        not isinstance(protospacer, str)
        or not re.fullmatch(r"[ACGT]{20}", protospacer)
    ):
        raise CrisprmeImportError(
            "CRISPRme declaration guide_protospacer must be exactly 20 uppercase DNA bases"
        )
    _required_text(payload["genome_build"], "genome_build")
    if payload["coordinate_system"] not in {
        "zero_based_half_open",
        "one_based_inclusive",
    }:
        raise CrisprmeImportError("unsupported CRISPRme coordinate_system")
    version = _required_text(payload["crisprme_version"], "crisprme_version")
    scoring_model = _required_text(payload["scoring_model"], "scoring_model")
    if scoring_model != CRISPRME_SCORING_MODEL:
        raise CrisprmeImportError(
            f"CRISPRme scoring_model must be {CRISPRME_SCORING_MODEL}"
        )

    selected_sheet = payload["selected_sheet"]
    if format_schema == CRISPRME_TSV_FORMAT_SCHEMA:
        if selected_sheet is not None:
            raise CrisprmeImportError("native CRISPRme TSV declarations require selected_sheet null")
        if version != CRISPRME_VERSION:
            raise CrisprmeImportError(
                f"native CRISPRme TSV contract requires version {CRISPRME_VERSION}"
            )
        if not filename.endswith("integrated_results.tsv"):
            raise CrisprmeImportError(
                "native CRISPRme input must be an integrated_results.tsv file"
            )
    else:
        expected_sheet = f"{guide_id}_integrated_results"
        if selected_sheet != expected_sheet:
            raise CrisprmeImportError(
                "publication XLSX declaration must select its guide integrated-results sheet"
            )

    databases = payload["variant_databases"]
    if not isinstance(databases, list) or not databases:
        raise CrisprmeImportError("CRISPRme declaration variant_databases must be a non-empty list")
    for index, database in enumerate(databases):
        database = _exact_object(database, {"name", "release"}, f"variant_databases[{index}]")
        _required_text(database["name"], f"variant_databases[{index}].name")
        _required_text(database["release"], f"variant_databases[{index}].release")

    patterns = payload["pam_patterns"]
    if not isinstance(patterns, list) or not patterns:
        raise CrisprmeImportError("CRISPRme declaration pam_patterns must be a non-empty list")
    if len(set(patterns)) != len(patterns):
        raise CrisprmeImportError("CRISPRme declaration pam_patterns contains duplicates")
    for pattern in patterns:
        if (
            not isinstance(pattern, str)
            or not pattern
            or pattern != pattern.upper()
            or set(pattern) - _IUPAC
        ):
            raise CrisprmeImportError("CRISPRme declaration contains an invalid PAM pattern")

    for field in ("max_mismatches", "max_dna_bulge", "max_rna_bulge"):
        _nonnegative_int(payload[field], field)

    retention = _exact_object(payload["retention"], _RETENTION_FIELDS, "retention")
    completeness = retention["completeness"]
    if completeness not in {"complete", "truncated_top_n", "filtered"}:
        raise CrisprmeImportError("unsupported CRISPRme retention.completeness")
    top_n = retention["top_n"]
    if completeness == "truncated_top_n":
        if isinstance(top_n, bool) or not isinstance(top_n, int) or top_n < 1:
            raise CrisprmeImportError(
                "truncated_top_n CRISPRme input requires a positive retention.top_n"
            )
    elif top_n is not None:
        raise CrisprmeImportError(
            "CRISPRme retention.top_n must be null unless completeness is truncated_top_n"
        )
    _validate_score_cutoff(retention["score_cutoff"])
    _nonnegative_int(retention["merge_window_bp"], "retention.merge_window_bp")
    if not isinstance(retention["all_candidate_coordinates_retained"], bool):
        raise CrisprmeImportError(
            "CRISPRme retention.all_candidate_coordinates_retained must be boolean"
        )
    if completeness == "complete" and (
        top_n is not None or retention["score_cutoff"] is not None
    ):
        raise CrisprmeImportError(
            "complete CRISPRme input cannot declare top-N or score truncation"
        )
    if completeness != "complete" and retention["all_candidate_coordinates_retained"]:
        raise CrisprmeImportError(
            "incomplete CRISPRme input cannot claim all candidate coordinates were retained"
        )

    _required_text(payload["source_command"], "source_command")
    return payload


def load_import_declaration(path: Path, native_path: Path) -> dict:
    """Read and check the required CRISPRme+ input and settings record."""
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise CrisprmeImportError(f"cannot read CRISPRme import declaration {path.name}") from exc
    return validate_import_declaration(payload, native_path)


def _validate_header(header: list[str], source_name: str) -> None:
    if not header or not any(header):
        raise CrisprmeImportError(f"{source_name}: missing CRISPRme header")
    if len(header) != len(set(header)):
        raise CrisprmeImportError(f"{source_name}: duplicate CRISPRme header")
    missing = sorted(REQUIRED_COLUMNS - set(header))
    if missing:
        raise CrisprmeImportError(
            f"{source_name}: missing required CRISPRme columns: {', '.join(missing)}"
        )


def _safe_row(header: list[str], values: list[str]) -> dict[str, str]:
    """Copy only the permitted fields; discard sample-level values."""
    return {
        field: (values[index] if index < len(values) else "").strip()
        for index, field in enumerate(header)
        if field in REQUIRED_COLUMNS
        and not field.casefold().startswith(_PROHIBITED_COLUMN_PREFIX)
    }


def _iter_tsv(path: Path) -> Iterator[tuple[int, dict[str, str], bool]]:
    try:
        handle = path.open(encoding="utf-8-sig", newline="")
    except OSError as exc:
        raise CrisprmeImportError(f"cannot read CRISPRme input {path.name}") from exc
    with handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as exc:
            raise CrisprmeImportError(f"{path.name}: missing CRISPRme header") from exc
        header = [field.strip() for field in header]
        _validate_header(header, path.name)
        for line_number, values in enumerate(reader, start=2):
            if not any(value.strip() for value in values):
                continue
            yield line_number, _safe_row(header, values), len(values) != len(header)


def _iter_xlsx(
    path: Path, selected_sheet: str
) -> Iterator[tuple[int, dict[str, str], bool]]:
    try:
        available = sheet_names(path)
    except Exception as exc:
        raise CrisprmeImportError(f"cannot read CRISPRme workbook {path.name}") from exc
    if selected_sheet not in available:
        # Available names are intentionally not rendered: private sheet names may contain an ID.
        raise CrisprmeImportError("declared CRISPRme integrated-results sheet is absent")
    try:
        rows = iter_xlsx_rows(path, selected_sheet)
        header = next(rows)
    except StopIteration as exc:
        raise CrisprmeImportError(f"{path.name}: missing CRISPRme header") from exc
    except Exception as exc:
        raise CrisprmeImportError(f"cannot read CRISPRme workbook {path.name}") from exc
    header = [field.strip() for field in header]
    _validate_header(header, path.name)
    try:
        for row_number, values in enumerate(rows, start=2):
            if not any(value.strip() for value in values):
                continue
            yield row_number, _safe_row(header, values), len(values) > len(header)
    except Exception as exc:
        raise CrisprmeImportError(f"cannot read CRISPRme workbook {path.name}") from exc


def _ignored_sheet_summary(path: Path, selected_sheet: str) -> dict:
    """Count skipped sheet types without retaining their names."""
    try:
        available = sheet_names(path)
    except Exception as exc:
        raise CrisprmeImportError(f"cannot read CRISPRme workbook {path.name}") from exc
    categories = set()
    ignored_count = 0
    for name in available:
        if name == selected_sheet:
            continue
        ignored_count += 1
        if name.endswith("_private_targets"):
            categories.add("private_targets")
        elif name.endswith("_altMerge"):
            categories.add("altMerge")
        else:
            categories.add("other")
    return {"count": ignored_count, "categories": sorted(categories)}


def _parse_integer(value: str, reason: str) -> int:
    try:
        parsed = Decimal(value)
    except (InvalidOperation, ValueError):
        raise _RowFailure(reason) from None
    if not parsed.is_finite() or parsed != parsed.to_integral_value():
        raise _RowFailure(reason)
    return int(parsed)


def _parse_float(value: str) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        raise _RowFailure("invalid_numeric_value") from None
    if not math.isfinite(parsed):
        raise _RowFailure("invalid_numeric_value")
    return parsed


def _null_text(value: str) -> str | None:
    stripped = value.strip()
    return None if not stripped or stripped.upper() == "NA" else stripped


def _alignment(value: str, *, required: bool) -> str | None:
    parsed = _null_text(value)
    if parsed is None:
        if required:
            raise _RowFailure("invalid_alignment")
        return None
    if set(parsed.upper()) - _DNA_ALIGNMENT:
        raise _RowFailure("invalid_alignment")
    return parsed


def _span(alignment: str) -> int:
    span = sum(base != "-" for base in alignment)
    if span < 1:
        raise _RowFailure("inconsistent_span")
    return span


def _observed_mismatches(guide_alignment: str, target_alignment: str) -> int:
    """Count known, non-gap base differences in the aligned guide region."""
    guide_position = 0
    mismatches = 0
    for guide_base, target_base in zip(
        guide_alignment.upper(), target_alignment.upper()
    ):
        if guide_base != "-":
            guide_position += 1
        if (
            1 <= guide_position <= 20
            and guide_base != "-"
            and target_base != "-"
            and guide_base != "N"
            and target_base != "N"
            and guide_base != target_base
        ):
            mismatches += 1
    return mismatches


def _pam_matches(pam: str, pattern: str) -> bool:
    return len(pam) == len(pattern) and all(
        base in _IUPAC_BASES[code] for base, code in zip(pam, pattern)
    )


def _internal_start(stored_start: int, coordinate_system: str) -> int:
    if coordinate_system == "zero_based_half_open":
        if stored_start < 0:
            raise _RowFailure("invalid_coordinate")
        return stored_start
    if stored_start < 1:
        raise _RowFailure("invalid_coordinate")
    return stored_start - 1


def _state(
    allele: str,
    alignment: str,
    cfd_text: str,
    chrom: str,
    start: int,
    strand: str,
    origin: str,
    rsid_raw: str | None,
    maf_raw: str | None,
    variant_spacer_raw: str | None,
    variant_genome_raw: str | None,
) -> dict:
    end = start + _span(alignment)
    return {
        "allele": allele,
        "alignment": alignment,
        "effective_interval": {
            "chrom": chrom,
            "start_0based": start,
            "end_0based_exclusive": end,
            "strand": strand,
        },
        "CFD": _parse_float(cfd_text),
        "variant_rsid_raw": rsid_raw,
        "variant_maf_raw": maf_raw,
        "variant_info_raw": {
            "spacer_pam": variant_spacer_raw,
            "genome": variant_genome_raw,
        },
        "origin_selected": allele == origin,
    }


def reference_mapped_states(states: list[dict]) -> list[dict]:
    """Keep ALT indel sequences without claiming a reference-coordinate mapping."""
    mapped = deepcopy(states)
    reference = next((state for state in mapped if state["allele"] == "ref"), None)
    for state in mapped:
        if state["allele"] == "ref":
            state["mapping_status"] = "reference_alignment"
        elif reference is None:
            state["effective_interval"] = None
            state["mapping_status"] = "unsupported_without_reference"
        elif (
            len(state["alignment"]) != len(reference["alignment"])
            or any(
                (a == "-") != (b == "-")
                for a, b in zip(state["alignment"], reference["alignment"])
            )
        ):
            state["effective_interval"] = None
            state["mapping_status"] = "unsupported_indel_alignment"
        else:
            state["effective_interval"] = deepcopy(reference["effective_interval"])
            state["mapping_status"] = "reference_alignment"
    for state in mapped:
        state["genomic_forward_changes"] = None
        if state["effective_interval"] is None or reference is None:
            continue
        ref_sequence = reference["alignment"].upper()
        alt_sequence = state["alignment"].upper()
        if "N" in ref_sequence or "N" in alt_sequence:
            continue
        interval = reference["effective_interval"]
        forward = interval["strand"] == "+"
        position = interval["start_0based"] if forward else interval["end_0based_exclusive"] - 1
        changes = []
        complement = str.maketrans("ACGT", "TGCA")
        for ref_base, alt_base in zip(ref_sequence, alt_sequence):
            if ref_base == "-":
                continue
            if ref_base != alt_base:
                changes.append({
                    "chrom": interval["chrom"], "position_0based": position,
                    "ref": ref_base if forward else ref_base.translate(complement),
                    "alt": alt_base if forward else alt_base.translate(complement),
                })
            position += 1 if forward else -1
        state["genomic_forward_changes"] = sorted(changes, key=lambda item: item["position_0based"])
    return mapped


def _record_id(row: dict) -> str:
    stable = {
        key: row[key]
        for key in (
            "guide_id",
            "genome_build",
            "chrom",
            "start_0based",
            "end_0based_exclusive",
            "strand",
            "guide_alignment",
            "selected_alignment",
            "canonical_alignment",
            "allele_states",
            "mismatch_count",
            "dna_bulge",
            "rna_bulge",
        )
    }
    digest = hashlib.sha256(
        json.dumps(stable, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()[:20]
    return f"crisprme:{row['guide_id']}:{digest}"


def _normalize_row(raw: dict[str, str], declaration: dict) -> dict:
    spacer_pam = raw["Spacer+PAM"].upper()
    if len(spacer_pam) < 20 or spacer_pam[:20] != declaration["guide_protospacer"]:
        raise _RowFailure("guide_mismatch")
    if set(spacer_pam) - _IUPAC:
        raise _RowFailure("guide_mismatch")

    guide_alignment = _alignment(raw["Aligned_spacer+PAM_(highest_CFD)"], required=True)
    assert guide_alignment is not None
    if guide_alignment.replace("-", "").upper()[:20] != declaration["guide_protospacer"]:
        raise _RowFailure("guide_mismatch")

    strand = raw["Strand_(highest_CFD)"].strip()
    if strand not in {"+", "-"}:
        raise _RowFailure("invalid_strand")
    chrom = raw["Chromosome"].strip()
    if not chrom or any(character.isspace() for character in chrom):
        raise _RowFailure("invalid_coordinate")
    stored_start = _parse_integer(
        raw["Start_coordinate_(highest_CFD)"], "invalid_coordinate"
    )
    start = _internal_start(stored_start, declaration["coordinate_system"])

    ref_alignment = _alignment(
        raw["Aligned_protospacer+PAM_REF_(highest_CFD)"], required=False
    )
    alt_alignment = _alignment(
        raw["Aligned_protospacer+PAM_ALT_(highest_CFD)"], required=False
    )
    if ref_alignment is None and alt_alignment is None:
        raise _RowFailure("invalid_alignment")
    for alignment in (ref_alignment, alt_alignment):
        if alignment is not None and len(alignment) != len(guide_alignment):
            raise _RowFailure("invalid_alignment")

    origin = raw["REF/ALT_origin_(highest_CFD)"].strip().lower()
    if origin not in {"ref", "alt"}:
        raise _RowFailure("invalid_alignment")
    selected_alignment = ref_alignment if origin == "ref" else alt_alignment
    if selected_alignment is None:
        raise _RowFailure("invalid_alignment")

    pam = raw["PAM_(highest_CFD)"].strip().upper()
    if (
        not pam
        or set(pam) - frozenset("ACGT")
        or len(spacer_pam) != 20 + len(pam)
        or len(guide_alignment.replace("-", "")) != 20 + len(pam)
        or "-" in guide_alignment[-len(pam):]
        or selected_alignment[-len(pam):].upper() != pam
        or not any(_pam_matches(pam, pattern) for pattern in declaration["pam_patterns"])
    ):
        raise _RowFailure("invalid_alignment")

    mismatch_count = _parse_integer(
        raw["Mismatches_(highest_CFD)"], "invalid_numeric_value"
    )
    bulges = _parse_integer(raw["Bulges_(highest_CFD)"], "invalid_numeric_value")
    if mismatch_count < 0 or bulges < 0:
        raise _RowFailure("invalid_numeric_value")
    bulge_type = raw["Bulge_type_(highest_CFD)"].strip().upper()
    if bulge_type not in {"X", "DNA", "RNA"}:
        raise _RowFailure("unsupported_bulge_type")
    if (bulge_type == "X" and bulges != 0) or (bulge_type != "X" and bulges < 1):
        raise _RowFailure("inconsistent_span")
    guide_gaps = guide_alignment.count("-")
    selected_gaps = selected_alignment.count("-")
    gap_topology = {
        "X": (0, 0),
        "DNA": (bulges, 0),
        "RNA": (0, bulges),
    }
    if (guide_gaps, selected_gaps) != gap_topology[bulge_type]:
        raise _RowFailure("inconsistent_span")
    if _observed_mismatches(guide_alignment, selected_alignment) != mismatch_count:
        raise _RowFailure("invalid_alignment")
    dna_bulge = bulges if bulge_type == "DNA" else 0
    rna_bulge = bulges if bulge_type == "RNA" else 0
    if (
        mismatch_count > declaration["max_mismatches"]
        or dna_bulge > declaration["max_dna_bulge"]
        or rna_bulge > declaration["max_rna_bulge"]
    ):
        raise _RowFailure("invalid_numeric_value")

    canonical_alignment = ref_alignment or alt_alignment
    assert canonical_alignment is not None
    # Without REF, retain the reported start as an anchor; the target span is unknown.
    end = start + (_span(ref_alignment) if ref_alignment is not None else 1)
    canonical_basis = "ref" if ref_alignment is not None else "alt_only"

    rsid_raw = _null_text(raw["Variant_rsID_(highest_CFD)"])
    maf_raw = _null_text(raw["Variant_MAF_(highest_CFD)"])
    variant_spacer_raw = _null_text(raw["Variant_info_spacer+PAM_(highest_CFD)"])
    variant_genome_raw = _null_text(raw["Variant_info_genome_(highest_CFD)"])
    states = []
    if ref_alignment is not None:
        states.append(
            _state(
                "ref",
                ref_alignment,
                raw["CFD_score_REF_(highest_CFD)"],
                chrom,
                start,
                strand,
                origin,
                rsid_raw,
                maf_raw,
                variant_spacer_raw,
                variant_genome_raw,
            )
        )
    if alt_alignment is not None:
        states.append(
            _state(
                "alt",
                alt_alignment,
                raw["CFD_score_ALT_(highest_CFD)"],
                chrom,
                start,
                strand,
                origin,
                rsid_raw,
                maf_raw,
                variant_spacer_raw,
                variant_genome_raw,
            )
        )
    states = reference_mapped_states(states)
    selected_cfd = _parse_float(raw["CFD_score_(highest_CFD)"])
    selected_state = next(state for state in states if state["origin_selected"])
    if not math.isclose(selected_cfd, selected_state["CFD"], rel_tol=1e-12, abs_tol=1e-15):
        raise _RowFailure("invalid_numeric_value")

    normalized = {
        "guide_id": declaration["guide_id"],
        "genome_build": declaration["genome_build"],
        "chrom": chrom,
        "start_0based": start,
        "end_0based_exclusive": end,
        "strand": strand,
        "guide_alignment": guide_alignment,
        "canonical_alignment": canonical_alignment,
        "selected_alignment": selected_alignment,
        "canonical_interval_basis": canonical_basis,
        "pam": pam,
        "mismatch_count": mismatch_count,
        "dna_bulge": dna_bulge,
        "rna_bulge": rna_bulge,
        "selected_cfd": selected_cfd,
        "selected_origin": origin,
        "allele_states": states,
        "variant_rsid_raw": rsid_raw,
        "variant_maf_raw": maf_raw,
        "pam_creation_raw": _null_text(raw["PAM_creation_(highest_CFD)"]),
        "not_found_in_ref_raw": _null_text(raw["Not_found_in_REF"]),
    }
    normalized["source_record_id"] = _record_id(normalized)
    return normalized


def import_crisprme(native_path: Path, declaration_path: Path) -> CrisprmeImportResult:
    """Import one declared CRISPRme result without retaining sample-level values."""
    declaration = load_import_declaration(declaration_path, native_path)
    if declaration["format_schema"] == CRISPRME_TSV_FORMAT_SCHEMA:
        native_rows: Iterable[tuple[int, dict[str, str], bool]] = _iter_tsv(native_path)
        sheets_ignored = {"count": 0, "categories": []}
    else:
        sheets_ignored = _ignored_sheet_summary(native_path, declaration["selected_sheet"])
        native_rows = _iter_xlsx(native_path, declaration["selected_sheet"])

    rows = []
    skipped = Counter()
    input_rows = 0
    seen_ids = set()
    for _row_number, raw, structural_error in native_rows:
        input_rows += 1
        if structural_error:
            skipped["invalid_alignment"] += 1
            continue
        try:
            row = _normalize_row(raw, declaration)
            if row["source_record_id"] in seen_ids:
                raise _RowFailure("invalid_alignment")
        except _RowFailure as exc:
            skipped[exc.reason] += 1
            continue
        seen_ids.add(row["source_record_id"])
        rows.append(row)

    if not rows:
        breakdown = ", ".join(f"{key}={skipped[key]}" for key in sorted(skipped))
        suffix = f" ({breakdown})" if breakdown else ""
        raise CrisprmeImportError(
            f"{native_path.name}: no valid CRISPRme rows were emitted{suffix}"
        )
    summary = {
        "input_data_rows": input_rows,
        "emitted_nomination_rows": len(rows),
        "skipped_rows": sum(skipped.values()),
        "skipped_by_reason": dict(sorted(skipped.items())),
    }
    original_declaration = deepcopy(declaration) if skipped else None
    if skipped:
        retention = declaration["retention"]
        retention["all_candidate_coordinates_retained"] = False
        if retention["completeness"] == "complete":
            retention["completeness"] = "filtered"
    return CrisprmeImportResult(
        tuple(rows),
        declaration,
        summary,
        sheets_ignored,
        True,
        original_declaration,
    )
