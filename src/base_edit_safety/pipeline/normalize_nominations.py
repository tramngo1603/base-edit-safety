"""Convert supported candidate files to the common candidate TSV format.

Read Cas-OFFinder, CRISPRme+, CHANGE-seq-BE identified tables or experimental BED6 lists.
Identify intended targets from their coordinate file and write companion source details with
input and reference checksums. This step does not run the search tools.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re

from .. import __version__
from ..ingestion.crisprme import (
    CRISPRME_XLSX_FORMAT_SCHEMA,
    CrisprmeImportError,
    CrisprmeImportResult,
    import_crisprme,
)
from .candidate_panel import (
    CandidatePanelError,
    SearchProfile,
    _pam_matches,
    _sha256,
    _site_key,
    load_guides,
    load_intended_loci,
    validate_crispor_source_declaration,
    validate_reference_sequences,
)
from .environment import run_stamp


NOMINATION_COLUMNS = [
    "guide_id",
    "nomination_kind",
    "source_method",
    "source_version_or_citation",
    "genome_build",
    "chrom",
    "start_0based",
    "end_0based_exclusive",
    "strand",
    "site_role",
    "site_sequence",
    "pam",
    "mismatch_count",
    "dna_bulge",
    "rna_bulge",
    "score_name",
    "score_value",
    "editable_base_in_window",
    "population_variant_id",
    "ch_driver_gene",
    "source_record_id",
    "search_profile",
    "source_metadata_json",
]

CHANGE_BE_REQUIRED = {
    "#Chromosome",
    "Start",
    "End",
    "Nuclease_Read_Count",
    "Strand",
    "Site_Sequence",
    "Site_Substitution_Number",
}


def _stable_record_id(prefix: str, guide_id: str, native_row: str) -> str:
    digest = hashlib.sha256(native_row.encode("utf-8")).hexdigest()[:16]
    return f"{prefix}:{guide_id}:{digest}"


def _intended_keys(loci: list[dict]) -> set[tuple]:
    return {_site_key(locus) for locus in loci}


def _role(
    guide_id: str,
    reference_id: str,
    chrom: str,
    start: int,
    end: int,
    intended: set[tuple],
) -> str:
    key = (guide_id, reference_id, chrom, start, end)
    # Use the target file's chromosome names here; the ledger checks aliases again.
    # Different source spellings must not change which sites count as intended.
    manifest_match = next(
        (
            locus_key
            for locus_key in intended
            if locus_key[0] == guide_id
            and locus_key[1] == reference_id
            and locus_key[3:] == (start, end)
            and _chrom_equivalent(locus_key[2], chrom)
        ),
        None,
    )
    return "intended" if manifest_match else "candidate"


def _chrom_equivalent(left: str, right: str) -> bool:
    def token(value: str) -> str:
        value = value.strip()
        return value[3:] if value.lower().startswith("chr") else value

    return token(left).upper().replace("MT", "M") == token(right).upper().replace("MT", "M")


def _editable(guide: dict, sequence: str) -> bool | None:
    concrete = sequence.replace("-", "").upper()
    spacer_length = len(guide["spacer_sequence"])
    if len(concrete) < spacer_length or set(concrete[:spacer_length]) - set("ACGT"):
        return None
    base = {"ABE": "A", "CBE": "C"}.get(guide["editor_class"].upper())
    if not base:
        return None
    start = guide["editing_window_start_1based"] - 1
    end = guide["editing_window_end_1based_inclusive"]
    return base in concrete[:spacer_length][start:end]


def _blank_row() -> dict[str, object]:
    return {column: "" for column in NOMINATION_COLUMNS}


def normalize_cas_offinder(
    native_path: Path,
    guides: list[dict],
    intended_loci: list[dict],
    profile: SearchProfile,
    source_version: str,
    guide_override: str | None = None,
) -> list[dict]:
    """Parse the official nine-column Cas-OFFinder output without filtering rows."""
    guide_by_id = {guide["guide_id"]: guide for guide in guides}
    intended = _intended_keys(intended_loci)
    rows = []
    seen_native = set()
    with native_path.open(encoding="utf-8", newline="") as handle:
        for line_number, line in enumerate(handle, start=1):
            native = line.rstrip("\r\n")
            if not native or native.startswith("#"):
                continue
            fields = native.split("\t")
            if len(fields) != 9:
                raise CandidatePanelError(
                    f"{native_path}:{line_number}: Cas-OFFinder output requires exactly 9 columns"
                )
            native_id, bulge_type, query, dna, chrom, start_text, strand, mm_text, bulge_text = fields
            guide_id = guide_override or native_id
            guide = guide_by_id.get(guide_id)
            if guide is None:
                raise CandidatePanelError(
                    f"{native_path}:{line_number}: Cas-OFFinder ID {guide_id!r} is absent from guides"
                )
            query_ungapped = query.replace("-", "").upper()
            if query_ungapped[: len(guide["spacer_sequence"])] != guide["spacer_sequence"]:
                raise CandidatePanelError(
                    f"{native_path}:{line_number}: Cas-OFFinder query does not match guide {guide_id!r}"
                )
            if bulge_type not in {"X", "DNA", "RNA"}:
                raise CandidatePanelError(
                    f"{native_path}:{line_number}: bulge type must be X, DNA, or RNA"
                )
            if strand not in {"+", "-"}:
                raise CandidatePanelError(
                    f"{native_path}:{line_number}: strand must be + or -"
                )
            try:
                start = int(start_text)
                mismatches = int(mm_text)
                bulge_size = int(bulge_text)
            except ValueError as exc:
                raise CandidatePanelError(
                    f"{native_path}:{line_number}: location, mismatches, and bulge size must be integers"
                ) from exc
            if min(start, mismatches, bulge_size) < 0:
                raise CandidatePanelError(
                    f"{native_path}:{line_number}: numeric fields must be non-negative"
                )
            if bulge_type == "X" and bulge_size != 0:
                raise CandidatePanelError(
                    f"{native_path}:{line_number}: X rows must have bulge size 0"
                )
            dna_sequence = dna.upper()
            if set(dna_sequence) - set("ACGTN-"):
                raise CandidatePanelError(
                    f"{native_path}:{line_number}: DNA sequence contains unsupported symbols"
                )
            interval_length = len(dna_sequence.replace("-", ""))
            if interval_length == 0:
                raise CandidatePanelError(f"{native_path}:{line_number}: empty DNA sequence")
            end = start + interval_length
            pam_length = len(guide["pam_class"])
            ungapped_dna = dna_sequence.replace("-", "")
            if pam_length < 1 or len(ungapped_dna) <= pam_length:
                raise CandidatePanelError(
                    f"{native_path}:{line_number}: cannot derive an observed PAM"
                )
            pam = ungapped_dna[-pam_length:]
            if not any(_pam_matches(pam, pattern) for pattern in profile.pam_patterns):
                raise CandidatePanelError(
                    f"{native_path}:{line_number}: observed PAM {pam!r} is outside the search profile"
                )
            record_id = _stable_record_id("casoff", guide_id, native)
            if record_id in seen_native:
                raise CandidatePanelError(
                    f"{native_path}:{line_number}: exact duplicate Cas-OFFinder row"
                )
            seen_native.add(record_id)
            row = _blank_row()
            row.update(
                {
                    "guide_id": guide_id,
                    "nomination_kind": "computational",
                    "source_method": "Cas-OFFinder",
                    "source_version_or_citation": source_version,
                    "genome_build": profile.reference_id,
                    "chrom": chrom,
                    "start_0based": start,
                    "end_0based_exclusive": end,
                    "strand": strand,
                    "site_role": _role(
                        guide_id, profile.reference_id, chrom, start, end, intended
                    ),
                    "site_sequence": dna_sequence,
                    "pam": pam,
                    "mismatch_count": mismatches,
                    "dna_bulge": bulge_size if bulge_type == "DNA" else 0,
                    "rna_bulge": bulge_size if bulge_type == "RNA" else 0,
                    "editable_base_in_window": (
                        _editable(guide, dna_sequence) if bulge_type == "X" else None
                    ),
                    "source_record_id": record_id,
                    "search_profile": profile.name,
                    "source_metadata_json": json.dumps(
                        {"cas_offinder_bulge_type": bulge_type, "native_query": query},
                        sort_keys=True,
                        separators=(",", ":"),
                    ),
                }
            )
            profile.validate_computational_nomination(
                {
                    "search_profile": profile.name,
                    "genome_build": profile.reference_id,
                    "mismatch_count": mismatches,
                    "dna_bulge": row["dna_bulge"],
                    "rna_bulge": row["rna_bulge"],
                    "pam": pam,
                    "population_variant_id": None,
                    "score_name": None,
                    "score_value": None,
                },
                f"{native_path}:{line_number}",
            )
            rows.append(row)
    if not rows:
        raise CandidatePanelError(f"{native_path}: no Cas-OFFinder rows found")
    return rows


def _read_change_be(path: Path) -> list[tuple[str, dict]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or [])
        missing = sorted(CHANGE_BE_REQUIRED - fields)
        if missing:
            raise CandidatePanelError(
                f"{path}: CHANGE-seq-BE table is missing columns: {', '.join(missing)}"
            )
        result = []
        for line_number, row in enumerate(reader, start=2):
            if None in row:
                raise CandidatePanelError(f"{path}:{line_number}: row has more values than columns")
            cleaned = {key: (value or "").strip() for key, value in row.items()}
            native = "\t".join(cleaned.get(name, "") for name in reader.fieldnames or [])
            result.append((native, cleaned))
        return result


def _read_bed6(path: Path) -> list[tuple[str, dict]]:
    result = []
    with path.open(encoding="utf-8", newline="") as handle:
        for line_number, line in enumerate(handle, start=1):
            native = line.rstrip("\r\n")
            if not native or native.startswith(("#", "track ", "browser ")):
                continue
            fields = native.split("\t")
            if len(fields) < 6:
                raise CandidatePanelError(f"{path}:{line_number}: BED input requires at least 6 columns")
            chrom, start, end, name, score, strand = fields[:6]
            result.append(
                (
                    native,
                    {
                        "#Chromosome": chrom,
                        "Start": start,
                        "End": end,
                        "Genomic Coordinate": name,
                        "Nuclease_Read_Count": score if score not in {"", "."} else "",
                        "Strand": strand,
                        "Site_Sequence": "",
                        "Site_Substitution_Number": "",
                    },
                )
            )
    return result


def normalize_experimental(
    native_path: Path,
    native_format: str,
    guide_id: str,
    guides: list[dict],
    intended_loci: list[dict],
    reference_id: str,
    source_method: str,
    source_version: str,
) -> list[dict]:
    """Normalize a CHANGE-seq-BE identified table or a BED6 nomination file."""
    guide_by_id = {guide["guide_id"]: guide for guide in guides}
    guide = guide_by_id.get(guide_id)
    if guide is None:
        raise CandidatePanelError(f"guide_id {guide_id!r} is absent from the guide manifest")
    intended = _intended_keys(intended_loci)
    native_rows = _read_change_be(native_path) if native_format == "change-seq-be" else _read_bed6(native_path)
    if not native_rows:
        raise CandidatePanelError(f"{native_path}: no experimental nomination rows found")

    rows = []
    seen_native = set()
    for index, (native, item) in enumerate(native_rows, start=1):
        try:
            start = int(item["Start"])
            end = int(item["End"])
        except ValueError as exc:
            raise CandidatePanelError(
                f"{native_path}:row {index}: Start and End must be integers"
            ) from exc
        if start < 0 or end <= start:
            raise CandidatePanelError(
                f"{native_path}:row {index}: coordinates must be zero-based, half-open"
            )
        strand = item["Strand"]
        if strand not in {"+", "-"}:
            raise CandidatePanelError(f"{native_path}:row {index}: strand must be + or -")
        mismatch_text = item.get("Site_Substitution_Number", "")
        if mismatch_text:
            try:
                mismatch = int(mismatch_text)
            except ValueError as exc:
                raise CandidatePanelError(
                    f"{native_path}:row {index}: Site_Substitution_Number must be an integer"
                ) from exc
            if mismatch < 0:
                raise CandidatePanelError(
                    f"{native_path}:row {index}: mismatch count must be non-negative"
                )
        else:
            mismatch = ""
        read_count_text = item.get("Nuclease_Read_Count", "")
        if read_count_text:
            try:
                read_count = float(read_count_text)
            except ValueError as exc:
                raise CandidatePanelError(
                    f"{native_path}:row {index}: read count/BED score must be numeric"
                ) from exc
            if read_count < 0:
                raise CandidatePanelError(
                    f"{native_path}:row {index}: read count/BED score must be non-negative"
                )
        else:
            read_count = None
        sequence = item.get("Site_Sequence", "").upper()
        if set(sequence) - set("ACGTN-"):
            raise CandidatePanelError(
                f"{native_path}:row {index}: Site_Sequence contains unsupported symbols"
            )
        record_id = _stable_record_id(
            "changebe" if native_format == "change-seq-be" else "bed6", guide_id, native
        )
        if record_id in seen_native:
            raise CandidatePanelError(f"{native_path}:row {index}: exact duplicate native row")
        seen_native.add(record_id)
        metadata = {
            "native_format": native_format,
            "native_name": item.get("Genomic Coordinate") or None,
        }
        for source_name, output_name in (
            ("Nuclease_Read_Count", "nuclease_read_count"),
            ("Control_Read_Count", "control_read_count"),
            ("Nuclease_overlap_bp_list", "nuclease_overlap_bp_list"),
            ("Control_overlap_bp_list", "control_overlap_bp_list"),
            ("File_Name", "replicate_file"),
            ("Cell", "biological_context_label"),
        ):
            if item.get(source_name):
                metadata[output_name] = (
                    Path(item[source_name]).name
                    if source_name == "File_Name"
                    else item[source_name]
                )
        row = _blank_row()
        row.update(
            {
                "guide_id": guide_id,
                "nomination_kind": "experimental",
                "source_method": source_method,
                "source_version_or_citation": source_version,
                "genome_build": reference_id,
                "chrom": item["#Chromosome"],
                "start_0based": start,
                "end_0based_exclusive": end,
                "strand": strand,
                "site_role": _role(
                    guide_id,
                    reference_id,
                    item["#Chromosome"],
                    start,
                    end,
                    intended,
                ),
                "site_sequence": sequence,
                "mismatch_count": mismatch,
                "score_name": (
                    "Nuclease_Read_Count" if native_format == "change-seq-be" else "BED_score"
                ) if read_count is not None else "",
                "score_value": read_count if read_count is not None else "",
                "editable_base_in_window": _editable(guide, sequence),
                "source_record_id": record_id,
                "source_metadata_json": json.dumps(
                    metadata, sort_keys=True, separators=(",", ":")
                ),
            }
        )
        rows.append(row)
    return rows


def _validate_sha256(value: str, label: str) -> str:
    lowered = value.lower()
    if not re.fullmatch(r"[0-9a-f]{64}", lowered):
        raise CandidatePanelError(f"{label} must be a 64-character SHA-256 digest")
    return lowered


def _crisprme_source_version(declaration: dict) -> str:
    name = (
        "CRISPRme"
        if declaration["format_schema"] == CRISPRME_XLSX_FORMAT_SCHEMA
        else "CRISPRme+"
    )
    return f"{name} {declaration['crisprme_version']}"


def normalize_crisprme(
    imported: CrisprmeImportResult,
    guides: list[dict],
    intended_loci: list[dict],
    profile: SearchProfile,
) -> list[dict]:
    """Convert checked CRISPRme+ rows to the common candidate TSV format."""
    declaration = imported.declaration
    guide_by_id = {guide["guide_id"]: guide for guide in guides}
    guide = guide_by_id.get(declaration["guide_id"])
    if guide is None:
        raise CandidatePanelError(
            "CRISPRme declaration guide_id is absent from the guide manifest"
        )
    if guide["spacer_sequence"] != declaration["guide_protospacer"]:
        raise CandidatePanelError(
            "CRISPRme declaration protospacer does not match the guide manifest"
        )
    if declaration["genome_build"] != profile.reference_id:
        raise CandidatePanelError(
            "CRISPRme declaration genome build does not match the search profile"
        )
    expected_limits = (
        profile.max_mismatches,
        profile.max_dna_bulge,
        profile.max_rna_bulge,
    )
    declared_limits = (
        declaration["max_mismatches"],
        declaration["max_dna_bulge"],
        declaration["max_rna_bulge"],
    )
    if declared_limits != expected_limits:
        raise CandidatePanelError(
            "CRISPRme mismatch/bulge limits do not match the search profile"
        )
    if set(declaration["pam_patterns"]) != set(profile.pam_patterns):
        raise CandidatePanelError("CRISPRme PAM patterns do not match the search profile")
    if not profile.variant_aware:
        raise CandidatePanelError("CRISPRme variant-aware import requires a variant-aware profile")

    intended = _intended_keys(intended_loci)
    rows = []
    for imported_row in imported.rows:
        metadata = {
            "crisprme": {
                "format_schema": declaration["format_schema"],
                "guide_alignment": imported_row["guide_alignment"],
                "selected_alignment": imported_row["selected_alignment"],
                "selected_origin": imported_row["selected_origin"],
                "canonical_interval_basis": imported_row["canonical_interval_basis"],
                "allele_states": imported_row["allele_states"],
                "variant_maf_raw": imported_row["variant_maf_raw"],
                "pam_creation_raw": imported_row["pam_creation_raw"],
                "not_found_in_ref_raw": imported_row["not_found_in_ref_raw"],
            }
        }
        row = _blank_row()
        row.update(
            {
                "guide_id": imported_row["guide_id"],
                "nomination_kind": "computational",
                "source_method": "CRISPRme",
                "source_version_or_citation": _crisprme_source_version(declaration),
                "genome_build": imported_row["genome_build"],
                "chrom": imported_row["chrom"],
                "start_0based": imported_row["start_0based"],
                "end_0based_exclusive": imported_row["end_0based_exclusive"],
                "strand": imported_row["strand"],
                "site_role": _role(
                    imported_row["guide_id"],
                    imported_row["genome_build"],
                    imported_row["chrom"],
                    imported_row["start_0based"],
                    imported_row["end_0based_exclusive"],
                    intended,
                ),
                # The canonical/ref alignment supports reference validation.  The selected
                # ref/alt state and its score remain explicit in source metadata.
                "site_sequence": (
                    imported_row["canonical_alignment"]
                    if imported_row["canonical_interval_basis"] == "ref" else ""
                ),
                "pam": imported_row["pam"],
                "mismatch_count": imported_row["mismatch_count"],
                "dna_bulge": imported_row["dna_bulge"],
                "rna_bulge": imported_row["rna_bulge"],
                "score_name": "CFD",
                "score_value": imported_row["selected_cfd"],
                "editable_base_in_window": "",
                "population_variant_id": imported_row["variant_rsid_raw"] or "",
                "source_record_id": imported_row["source_record_id"],
                "search_profile": profile.name,
                "source_metadata_json": json.dumps(
                    metadata, sort_keys=True, separators=(",", ":")
                ),
                # Used by reference validation before write_normalized drops non-schema keys.
                "source_metadata": metadata,
            }
        )
        profile.validate_computational_nomination(
            {
                "search_profile": profile.name,
                "genome_build": profile.reference_id,
                "mismatch_count": imported_row["mismatch_count"],
                "dna_bulge": imported_row["dna_bulge"],
                "rna_bulge": imported_row["rna_bulge"],
                "pam": imported_row["pam"],
                "population_variant_id": imported_row["variant_rsid_raw"],
                "score_name": "CFD",
                "score_value": imported_row["selected_cfd"],
            },
            f"CRISPRme source_record_id={imported_row['source_record_id']}",
        )
        rows.append(row)
    return rows


def write_normalized(rows: list[dict], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=NOMINATION_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in NOMINATION_COLUMNS})


def write_provenance(
    output_path: Path,
    normalized_path: Path,
    native_path: Path,
    mode: str,
    source_method: str,
    source_version: str,
    source_command: str,
    reference_id: str,
    reference_path: Path,
    guide_path: Path,
    intended_path: Path,
    profile_path: Path | None,
    row_count: int,
    source_declarations: dict | None = None,
    import_summary: dict | None = None,
    *,
    crispor_declaration: dict | None = None,
    sheets_ignored: dict | None = None,
    privacy_filtering_applied: bool | None = None,
    original_source_declaration: dict | None = None,
) -> None:
    index_path = Path(str(reference_path) + ".fai")
    reference_sha256 = _sha256(reference_path)
    if (source_declarations is None) != (import_summary is None):
        raise CandidatePanelError(
            "source declarations and import summary must be supplied together"
        )
    if source_declarations is not None and crispor_declaration is not None:
        raise CandidatePanelError(
            "CRISPRme and CRISPOR declarations require separate source provenance files"
        )
    if source_declarations is not None and (
        sheets_ignored is None or privacy_filtering_applied is not True
    ):
        raise CandidatePanelError(
            "CRISPRme provenance requires sheets_ignored and privacy_filtering_applied=true"
        )
    if crispor_declaration is not None:
        crispor_declaration = validate_crispor_source_declaration(crispor_declaration)
    payload = {
        "record_type": "nomination_source_provenance",
        "schema_version": "2" if source_declarations is not None else "1",
        "normalized_output": {
            "file": normalized_path.name,
            "sha256": _sha256(normalized_path),
            "row_count": row_count,
        },
        "native_input": {"file": native_path.name, "sha256": _sha256(native_path)},
        "source_run": {
            "method": source_method,
            "version_or_citation": source_version,
            "command": source_command,
            "reference_id": reference_id,
            "reference_sha256": _validate_sha256(reference_sha256, "reference_sha256"),
            "reference_checksum_status": "recorded",
        },
        "reference_fasta": {
            "file": reference_path.name,
            "sha256": reference_sha256,
            "index_file": index_path.name,
            "index_sha256": _sha256(index_path),
        },
        "adapter": {
            "mode": mode,
            "base_edit_safety_version": __version__,
        },
        "guide_manifest": {"file": guide_path.name, "sha256": _sha256(guide_path)},
        "intended_locus_manifest": {
            "file": intended_path.name,
            "sha256": _sha256(intended_path),
        },
        "search_profile": (
            {"file": profile_path.name, "sha256": _sha256(profile_path)}
            if profile_path is not None
            else None
        ),
        "generated_at_utc": run_stamp(),
    }
    if source_declarations is not None:
        payload["source_declarations"] = source_declarations
        payload["import_summary"] = import_summary
        payload["sheets_ignored"] = sheets_ignored
        payload["privacy_filtering_applied"] = privacy_filtering_applied
        if original_source_declaration is not None:
            payload["original_source_declaration"] = original_source_declaration
    elif crispor_declaration is not None:
        payload["source_declarations"] = {"crispor": crispor_declaration}
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _common(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--input", type=Path, required=True, help="native source file")
    parser.add_argument("--guides", type=Path, required=True, help="guide manifest TSV")
    parser.add_argument(
        "--intended-loci", type=Path, required=True, help="intended-locus manifest TSV"
    )
    parser.add_argument("--source-version", required=True, help="tool version or source citation")
    parser.add_argument(
        "--source-command",
        required=True,
        help="exact generation/acquisition command, or an explicit not-applicable reason",
    )
    parser.add_argument(
        "--reference-fasta",
        type=Path,
        required=True,
        help="reference FASTA with a samtools .fai index",
    )
    parser.add_argument("--output", type=Path, required=True, help="normalized nomination TSV")
    parser.add_argument(
        "--provenance-output", type=Path, required=True, help="source-run provenance JSON"
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="mode", required=True)

    cas = subparsers.add_parser("cas-offinder", help="normalize official nine-column output")
    _common(cas)
    cas.add_argument("--profile", type=Path, required=True, help="search-profile JSON")
    cas.add_argument(
        "--guide-id",
        help="single guide override when native output uses numeric IDs; otherwise native IDs are used",
    )

    experimental = subparsers.add_parser(
        "experimental", help="normalize CHANGE-seq-BE identified or BED6 nominations"
    )
    _common(experimental)
    experimental.add_argument(
        "--format", choices=("change-seq-be", "bed6"), required=True
    )
    experimental.add_argument("--guide-id", required=True)
    experimental.add_argument("--reference-id", required=True)
    experimental.add_argument(
        "--source-method", required=True, help="for example CHANGE-seq-BE or CIRCLE-seq"
    )

    crisprme = subparsers.add_parser(
        "crisprme", help="normalize a declared CRISPRme integrated-results file"
    )
    crisprme.add_argument("--input", type=Path, required=True, help="native CRISPRme TSV or validation XLSX")
    crisprme.add_argument(
        "--declaration", type=Path, required=True, help="mandatory CRISPRme import declaration JSON"
    )
    crisprme.add_argument("--guides", type=Path, required=True, help="guide manifest TSV")
    crisprme.add_argument(
        "--intended-loci", type=Path, required=True, help="intended-locus manifest TSV"
    )
    crisprme.add_argument("--profile", type=Path, required=True, help="search-profile JSON")
    crisprme.add_argument(
        "--reference-fasta",
        type=Path,
        required=True,
        help="reference FASTA with a samtools .fai index",
    )
    crisprme.add_argument("--output", type=Path, required=True, help="normalized nomination TSV")
    crisprme.add_argument(
        "--provenance-output", type=Path, required=True, help="source-run provenance JSON"
    )

    args = parser.parse_args(argv)
    try:
        guides = load_guides(args.guides)
        loci = load_intended_loci(args.intended_loci)
        source_declarations = None
        import_summary = None
        if args.mode == "cas-offinder":
            profile = SearchProfile.load(args.profile)
            rows = normalize_cas_offinder(
                args.input,
                guides,
                loci,
                profile,
                args.source_version,
                args.guide_id,
            )
            reference_id = profile.reference_id
            source_method = "Cas-OFFinder"
            profile_path = args.profile
            source_version = args.source_version
            source_command = args.source_command
            adapter_mode = "cas-offinder"
        elif args.mode == "experimental":
            rows = normalize_experimental(
                args.input,
                args.format,
                args.guide_id,
                guides,
                loci,
                args.reference_id,
                args.source_method,
                args.source_version,
            )
            reference_id = args.reference_id
            source_method = args.source_method
            profile_path = None
            source_version = args.source_version
            source_command = args.source_command
            adapter_mode = args.format
        else:
            profile = SearchProfile.load(args.profile)
            imported = import_crisprme(args.input, args.declaration)
            rows = normalize_crisprme(imported, guides, loci, profile)
            reference_id = profile.reference_id
            source_method = "CRISPRme"
            source_version = _crisprme_source_version(imported.declaration)
            source_command = imported.declaration["source_command"]
            profile_path = args.profile
            source_declarations = {"crisprme": imported.declaration}
            import_summary = imported.import_summary
            adapter_mode = "crisprme"
        validate_reference_sequences(rows, args.reference_fasta)
        write_normalized(rows, args.output)
        write_provenance(
            args.provenance_output,
            args.output,
            args.input,
            adapter_mode,
            source_method,
            source_version,
            source_command,
            reference_id,
            args.reference_fasta,
            args.guides,
            args.intended_loci,
            profile_path,
            len(rows),
            source_declarations,
            import_summary,
            original_source_declaration=(imported.original_declaration if args.mode == "crisprme" else None),
            sheets_ignored=(imported.sheets_ignored if args.mode == "crisprme" else None),
            privacy_filtering_applied=(
                imported.privacy_filtering_applied if args.mode == "crisprme" else None
            ),
        )
    except (CandidatePanelError, CrisprmeImportError, OSError) as exc:
        parser.error(str(exc))
    print(f"wrote {len(rows)} normalized nominations to {args.output}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
