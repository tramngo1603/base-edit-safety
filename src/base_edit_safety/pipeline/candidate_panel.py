"""Combine candidate lists and experiment records for review.

Match supplied measurements to reported DNA versions at exact sites. Return candidate
records with separate source results, unresolved fields and researcher panel decisions.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from dataclasses import dataclass
import hashlib
import json
import math
import platform
from pathlib import Path
from typing import Iterable

from .. import __version__, anchors
from ..measurement_units import editing_detection_limit
from ..ingestion.crisprme import (
    SKIP_REASONS,
    SUPPORTED_FORMAT_SCHEMAS,
    CrisprmeImportError,
    validate_import_declaration,
    reference_mapped_states,
)
from .environment import repo_state, run_stamp


LEDGER_SCHEMA_VERSION = "3"

CRISPOR_MAX_MISMATCHES = 4
CRISPOR_FOUR_MISMATCH_MIN_CFD = 0.10

CRISPOR_DECLARATION_FIELDS = frozenset(
    {
        "record_type",
        "schema_version",
        "guide_id",
        "guide_protospacer",
        "genome_build",
        "pam_patterns",
        "max_mismatches",
        "max_dna_bulge",
        "max_rna_bulge",
        "retention",
    }
)
CRISPOR_RETENTION_FIELDS = frozenset(
    {"completeness", "top_n", "score_cutoff", "all_candidate_coordinates_retained"}
)

DECISION_REASONS = frozenset(
    {
        "not_reviewed",
        "low_mm_high_cfd",
        "bulge_site",
        "variant_created_pam",
        "variant_reduced_mm",
        "ch_driver_overlap",
        "genic_context",
        "lab_nominated",
        "no_base_in_usual_window",
        "redundant_with",
        "other_note",
    }
)

GUIDE_FIELDS = {
    "guide_id",
    "spacer_sequence",
    "editor",
    "editor_class",
    "pam_class",
    "editing_window_start_1based",
    "editing_window_end_1based_inclusive",
    "intended_target_description",
    "locus_system",
}
GUIDE_REQUIRED = GUIDE_FIELDS

INTENDED_LOCUS_FIELDS = {
    "guide_id",
    "genome_build",
    "chrom",
    "start_0based",
    "end_0based_exclusive",
    "strand",
    "locus_name",
}
INTENDED_LOCUS_REQUIRED = INTENDED_LOCUS_FIELDS

NOMINATION_FIELDS = {
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
}
NOMINATION_REQUIRED = {
    "guide_id",
    "nomination_kind",
    "source_method",
    "source_version_or_citation",
    "genome_build",
    "chrom",
    "start_0based",
    "end_0based_exclusive",
    "source_record_id",
}

MEASUREMENT_CONTEXT_FIELDS = {
    "editor", "measurement_kind", "allele", "allele_sequence", "allele_evidence",
    "allele_resolution", "configuration_id", "haplotype_id", "phase_status", "qc_note",
}

MEASUREMENT_FIELDS = {
    *MEASUREMENT_CONTEXT_FIELDS,
    "guide_id",
    "genome_build",
    "chrom",
    "start_0based",
    "end_0based_exclusive",
    "strand",
    "assay",
    "source_version_or_citation",
    "status",
    "value",
    "unit",
    "limit_of_detection_value",
    "limit_of_detection_unit",
    "biological_context",
    "source_record_id",
    "replicate_id",
    "treated_group",
    "control_group",
    "read_depth",
    "qc_status",
    "statistical_test",
    "p_value",
    "adjusted_p_value",
}
MEASUREMENT_REQUIRED = {
    "guide_id",
    "genome_build",
    "chrom",
    "start_0based",
    "end_0based_exclusive",
    "assay",
    "source_version_or_citation",
    "status",
    "biological_context",
    "source_record_id",
}

DECISION_FIELDS = {
    "guide_id",
    "genome_build",
    "chrom",
    "start_0based",
    "end_0based_exclusive",
    "decision",
    "reason",
    "note",
    "related_site_id",
}
DECISION_REQUIRED = DECISION_FIELDS - {"note", "related_site_id"}

IUPAC = {
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


class CandidatePanelError(ValueError):
    """Raised when an input cannot be interpreted without guessing."""


def validate_crispor_source_declaration(value: object) -> dict:
    """Check the recorded CRISPOR search settings and rules for keeping results."""
    if not isinstance(value, dict) or set(value) != CRISPOR_DECLARATION_FIELDS:
        raise CandidatePanelError(
            "CRISPOR source declaration must contain exactly: "
            + ", ".join(sorted(CRISPOR_DECLARATION_FIELDS))
        )
    if value["record_type"] != "crispor_search_declaration" or value["schema_version"] != "1":
        raise CandidatePanelError("unsupported CRISPOR source declaration record or schema")
    guide_id = _required_text(value["guide_id"], "CRISPOR declaration guide_id")
    protospacer = value["guide_protospacer"]
    if not isinstance(protospacer, str) or len(protospacer) != 20 or set(protospacer) - set("ACGT"):
        raise CandidatePanelError(
            "CRISPOR declaration guide_protospacer must be exactly 20 uppercase DNA bases"
        )
    genome_build = _required_text(
        value["genome_build"], "CRISPOR declaration genome_build"
    )
    patterns = value["pam_patterns"]
    if not isinstance(patterns, list) or not patterns or len(patterns) != len(set(patterns)):
        raise CandidatePanelError(
            "CRISPOR declaration pam_patterns must be a non-empty list without duplicates"
        )
    normalized_patterns = [_validate_pam_pattern(pattern) for pattern in patterns]
    limits = {
        field: _nonnegative_int(value[field], f"CRISPOR declaration {field}")
        for field in ("max_mismatches", "max_dna_bulge", "max_rna_bulge")
    }
    if limits["max_dna_bulge"] != 0 or limits["max_rna_bulge"] != 0:
        raise CandidatePanelError("CRISPOR declaration bulge limits must both be zero")

    retention = value["retention"]
    if not isinstance(retention, dict) or set(retention) != CRISPOR_RETENTION_FIELDS:
        raise CandidatePanelError(
            "CRISPOR retention declaration must contain exactly: "
            + ", ".join(sorted(CRISPOR_RETENTION_FIELDS))
        )
    completeness = retention["completeness"]
    if completeness not in {"complete", "truncated_top_n", "filtered"}:
        raise CandidatePanelError("unsupported CRISPOR retention completeness")
    top_n = retention["top_n"]
    if completeness == "truncated_top_n":
        if isinstance(top_n, bool) or not isinstance(top_n, int) or top_n < 1:
            raise CandidatePanelError(
                "truncated_top_n CRISPOR provenance requires a positive retention.top_n"
            )
    elif top_n is not None:
        raise CandidatePanelError(
            "CRISPOR retention.top_n must be null unless completeness is truncated_top_n"
        )
    score_cutoff = retention["score_cutoff"]
    if score_cutoff is not None:
        if not isinstance(score_cutoff, dict) or set(score_cutoff) != {"score_name", "minimum"}:
            raise CandidatePanelError(
                "CRISPOR retention.score_cutoff must contain score_name and minimum"
            )
        _required_text(
            score_cutoff["score_name"], "CRISPOR retention.score_cutoff.score_name"
        )
        minimum = score_cutoff["minimum"]
        if isinstance(minimum, bool) or not isinstance(minimum, (int, float)):
            raise CandidatePanelError("CRISPOR retention score minimum must be finite")
        if not math.isfinite(float(minimum)):
            raise CandidatePanelError("CRISPOR retention score minimum must be finite")
    all_retained = retention["all_candidate_coordinates_retained"]
    if not isinstance(all_retained, bool):
        raise CandidatePanelError(
            "CRISPOR retention.all_candidate_coordinates_retained must be boolean"
        )
    if completeness == "complete" and (
        top_n is not None or score_cutoff is not None or all_retained is not True
    ):
        raise CandidatePanelError(
            "complete CRISPOR provenance cannot declare truncation or omitted coordinates"
        )
    if completeness != "complete" and all_retained:
        raise CandidatePanelError(
            "incomplete CRISPOR provenance cannot claim all candidate coordinates were retained"
        )
    return {
        **value,
        "guide_id": guide_id,
        "genome_build": genome_build,
        "pam_patterns": normalized_patterns,
        **limits,
        "retention": {
            **retention,
            "score_cutoff": (
                None
                if score_cutoff is None
                else {
                    "score_name": score_cutoff["score_name"].strip(),
                    "minimum": float(score_cutoff["minimum"]),
                }
            ),
        },
    }


@dataclass(frozen=True)
class SearchProfile:
    """Search settings that imported candidate rows must satisfy."""

    name: str
    max_mismatches: int
    max_dna_bulge: int
    max_rna_bulge: int
    pam_patterns: tuple[str, ...]
    variant_aware: bool
    intended_use: str
    reference_id: str
    reference_provenance: str
    variant_source: str | None
    retention_rule: str
    retention_score_name: str | None
    minimum_score: float | None

    @classmethod
    def load(cls, path: Path) -> "SearchProfile":
        try:
            raw = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise CandidatePanelError(f"cannot read search profile {path}: {exc}") from exc
        if not isinstance(raw, dict):
            raise CandidatePanelError("search profile must be a JSON object")

        required = {
            "name",
            "max_mismatches",
            "max_dna_bulge",
            "max_rna_bulge",
            "pam_patterns",
            "variant_aware",
            "intended_use",
            "reference_id",
            "reference_provenance",
            "variant_source",
            "retention_rule",
            "retention_score_name",
            "minimum_score",
        }
        missing = sorted(required - raw.keys())
        unknown = sorted(raw.keys() - required)
        if missing or unknown:
            details = []
            if missing:
                details.append(f"missing fields: {', '.join(missing)}")
            if unknown:
                details.append(f"unknown fields: {', '.join(unknown)}")
            raise CandidatePanelError("invalid search profile (" + "; ".join(details) + ")")

        name = _required_text(raw["name"], "search profile name")
        intended_use = _required_text(raw["intended_use"], "search profile intended_use")
        reference_id = _required_text(raw["reference_id"], "search profile reference_id")
        reference_provenance = _required_text(
            raw["reference_provenance"], "search profile reference_provenance"
        )
        retention_rule = _required_text(raw["retention_rule"], "search profile retention_rule")
        retention_score_name = raw["retention_score_name"]
        minimum_score = raw["minimum_score"]
        if (retention_score_name is None) != (minimum_score is None):
            raise CandidatePanelError(
                "retention_score_name and minimum_score must both be set or both be null"
            )
        if retention_score_name is not None:
            retention_score_name = _required_text(
                retention_score_name, "search profile retention_score_name"
            )
            if isinstance(minimum_score, bool):
                raise CandidatePanelError("minimum_score must be a finite number")
            try:
                minimum_score = float(minimum_score)
            except (TypeError, ValueError) as exc:
                raise CandidatePanelError("minimum_score must be a finite number") from exc
            if not math.isfinite(minimum_score):
                raise CandidatePanelError("minimum_score must be a finite number")
        max_mismatches = _nonnegative_int(raw["max_mismatches"], "max_mismatches")
        max_dna_bulge = _nonnegative_int(raw["max_dna_bulge"], "max_dna_bulge")
        max_rna_bulge = _nonnegative_int(raw["max_rna_bulge"], "max_rna_bulge")
        if not isinstance(raw["variant_aware"], bool):
            raise CandidatePanelError("variant_aware must be true or false")
        variant_source = raw["variant_source"]
        if raw["variant_aware"]:
            variant_source = _required_text(variant_source, "variant_source")
        elif variant_source is not None:
            raise CandidatePanelError("variant_source must be null when variant_aware is false")
        if not isinstance(raw["pam_patterns"], list) or not raw["pam_patterns"]:
            raise CandidatePanelError("pam_patterns must be a non-empty JSON list")
        patterns = tuple(_validate_pam_pattern(item) for item in raw["pam_patterns"])
        return cls(
            name=name,
            max_mismatches=max_mismatches,
            max_dna_bulge=max_dna_bulge,
            max_rna_bulge=max_rna_bulge,
            pam_patterns=patterns,
            variant_aware=raw["variant_aware"],
            intended_use=intended_use,
            reference_id=reference_id,
            reference_provenance=reference_provenance,
            variant_source=variant_source,
            retention_rule=retention_rule,
            retention_score_name=retention_score_name,
            minimum_score=minimum_score,
        )

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "max_mismatches": self.max_mismatches,
            "max_dna_bulge": self.max_dna_bulge,
            "max_rna_bulge": self.max_rna_bulge,
            "pam_patterns": list(self.pam_patterns),
            "variant_aware": self.variant_aware,
            "intended_use": self.intended_use,
            "reference_id": self.reference_id,
            "reference_provenance": self.reference_provenance,
            "variant_source": self.variant_source,
            "retention_rule": self.retention_rule,
            "retention_score_name": self.retention_score_name,
            "minimum_score": self.minimum_score,
        }

    def validate_computational_nomination(self, nomination: dict, location: str) -> None:
        """Check that an imported row falls within the declared computational search."""
        if nomination.get("search_profile") != self.name:
            raise CandidatePanelError(
                f"{location}: computational row search_profile must be {self.name!r}"
            )
        if nomination.get("genome_build") != self.reference_id:
            raise CandidatePanelError(
                f"{location}: genome_build must match profile reference_id {self.reference_id!r}"
            )
        mismatch_count = nomination.get("mismatch_count")
        if mismatch_count is None:
            raise CandidatePanelError(f"{location}: computational row requires mismatch_count")
        if mismatch_count > self.max_mismatches:
            raise CandidatePanelError(
                f"{location}: mismatch_count {mismatch_count} exceeds declared maximum "
                f"{self.max_mismatches}"
            )
        for field, maximum in (
            ("dna_bulge", self.max_dna_bulge),
            ("rna_bulge", self.max_rna_bulge),
        ):
            value = nomination.get(field) or 0
            if value > maximum:
                raise CandidatePanelError(
                    f"{location}: {field} {value} exceeds declared maximum {maximum}"
                )
        pam = nomination.get("pam")
        if not pam:
            raise CandidatePanelError(f"{location}: computational row requires pam")
        if not any(_pam_matches(pam, pattern) for pattern in self.pam_patterns):
            allowed = ", ".join(self.pam_patterns)
            raise CandidatePanelError(f"{location}: PAM {pam!r} does not match {allowed}")
        if nomination.get("population_variant_id") and not self.variant_aware:
            raise CandidatePanelError(
                f"{location}: population_variant_id is present but profile is not variant-aware"
            )
        if self.minimum_score is not None:
            if nomination.get("score_name") != self.retention_score_name:
                raise CandidatePanelError(
                    f"{location}: profile requires score_name {self.retention_score_name!r}"
                )
            score = nomination.get("score_value")
            if score is None or score < self.minimum_score:
                raise CandidatePanelError(
                    f"{location}: score_value must be at least {self.minimum_score}"
                )


def _required_text(value: object, field: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise CandidatePanelError(f"{field} must be a non-empty string")
    return value.strip()


def _nonnegative_int(value: object, field: str) -> int:
    if isinstance(value, bool):
        raise CandidatePanelError(f"{field} must be a non-negative integer")
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        raise CandidatePanelError(f"{field} must be a non-negative integer") from exc
    if parsed < 0 or str(parsed) != str(value).strip():
        raise CandidatePanelError(f"{field} must be a non-negative integer")
    return parsed


def _optional_int(value: str, field: str, location: str) -> int | None:
    if value == "":
        return None
    try:
        return _nonnegative_int(value, field)
    except CandidatePanelError as exc:
        raise CandidatePanelError(f"{location}: {exc}") from exc


def _optional_float(value: str, field: str, location: str) -> float | None:
    if value == "":
        return None
    try:
        parsed = float(value)
    except ValueError as exc:
        raise CandidatePanelError(f"{location}: {field} must be numeric") from exc
    if not math.isfinite(parsed):
        raise CandidatePanelError(f"{location}: {field} must be finite")
    return parsed


def _optional_bool(value: str, field: str, location: str) -> bool | None:
    if value == "":
        return None
    lowered = value.lower()
    if lowered == "true":
        return True
    if lowered == "false":
        return False
    raise CandidatePanelError(f"{location}: {field} must be true, false, or blank")


def _validate_pam_pattern(value: object) -> str:
    pattern = _required_text(value, "PAM pattern").upper()
    invalid = sorted(set(pattern) - IUPAC.keys())
    if invalid:
        raise CandidatePanelError(f"PAM pattern {pattern!r} has invalid symbols: {''.join(invalid)}")
    return pattern


def _pam_matches(pam: str, pattern: str) -> bool:
    observed = pam.upper()
    if len(observed) != len(pattern) or any(base not in "ACGT" for base in observed):
        return False
    return all(base in IUPAC[code] for base, code in zip(observed, pattern))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_tsv(path: Path, allowed: set[str], required: set[str]) -> Iterable[tuple[int, dict]]:
    try:
        handle = path.open(encoding="utf-8-sig", newline="")
    except OSError as exc:
        raise CandidatePanelError(f"cannot read {path}: {exc}") from exc
    with handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise CandidatePanelError(f"{path}: missing TSV header")
        if len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise CandidatePanelError(f"{path}: duplicate TSV header")
        fields = set(reader.fieldnames)
        missing = sorted(required - fields)
        unknown = sorted(fields - allowed)
        if missing or unknown:
            details = []
            if missing:
                details.append(f"missing columns: {', '.join(missing)}")
            if unknown:
                details.append(f"unknown columns: {', '.join(unknown)}")
            raise CandidatePanelError(f"{path}: " + "; ".join(details))
        for line_number, raw in enumerate(reader, start=2):
            if None in raw:
                raise CandidatePanelError(f"{path}:{line_number}: row has more values than columns")
            row = {field: (value or "").strip() for field, value in raw.items()}
            if not any(row.values()):
                continue
            empty_required = sorted(field for field in required if not row.get(field))
            if empty_required:
                raise CandidatePanelError(
                    f"{path}:{line_number}: empty required fields: {', '.join(empty_required)}"
                )
            yield line_number, row


def _canonical_chrom(value: str) -> str:
    """Normalize standard human chromosome aliases while leaving named contigs unchanged."""
    chrom = value.strip()
    if not chrom:
        raise CandidatePanelError("chrom must be a non-empty string")
    token = chrom[3:] if chrom.lower().startswith("chr") else chrom
    upper = token.upper()
    if upper in {"X", "Y"}:
        return f"chr{upper}"
    if upper in {"M", "MT"}:
        return "chrM"
    if token.isdigit() and 1 <= int(token) <= 22:
        return f"chr{int(token)}"
    return chrom


def _fasta_index(reference_path: Path) -> tuple[Path, dict[str, tuple[int, int, int, int]]]:
    """Load a samtools-compatible FASTA index, keyed by canonical chromosome name."""
    index_path = Path(str(reference_path) + ".fai")
    try:
        lines = index_path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise CandidatePanelError(
            f"cannot read {index_path}; create it with 'samtools faidx {reference_path.name}'"
        ) from exc
    entries = {}
    for line_number, line in enumerate(lines, start=1):
        fields = line.split("\t")
        if len(fields) < 5:
            raise CandidatePanelError(f"{index_path}:{line_number}: invalid FASTA index row")
        name = fields[0]
        try:
            length, offset, line_bases, line_width = map(int, fields[1:5])
        except ValueError as exc:
            raise CandidatePanelError(
                f"{index_path}:{line_number}: invalid numeric FASTA index fields"
            ) from exc
        if min(length, offset, line_bases, line_width) < 0 or line_bases < 1 or line_width < line_bases:
            raise CandidatePanelError(f"{index_path}:{line_number}: invalid FASTA index geometry")
        canonical = _canonical_chrom(name)
        if canonical in entries:
            raise CandidatePanelError(
                f"{index_path}: multiple FASTA records collapse to chromosome alias {canonical!r}"
            )
        entries[canonical] = (length, offset, line_bases, line_width)
    if not entries:
        raise CandidatePanelError(f"{index_path}: empty FASTA index")
    return index_path, entries


def _fetch_fasta_interval(
    handle,
    entry: tuple[int, int, int, int],
    start: int,
    end: int,
    location: str,
) -> str:
    length, offset, line_bases, line_width = entry
    if start < 0 or end <= start or end > length:
        raise CandidatePanelError(
            f"{location}: interval {start}-{end} lies outside the indexed reference sequence"
        )
    first_byte = offset + (start // line_bases) * line_width + (start % line_bases)
    final_base = end - 1
    final_byte = offset + (final_base // line_bases) * line_width + (final_base % line_bases)
    handle.seek(first_byte)
    raw = handle.read(final_byte - first_byte + 1)
    sequence = b"".join(raw.split()).decode("ascii").upper()
    expected = end - start
    if len(sequence) != expected or set(sequence) - set("ACGTN"):
        raise CandidatePanelError(f"{location}: reference FASTA interval could not be read exactly")
    return sequence


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def validate_reference_sequences(nominations: list[dict], reference_path: Path) -> dict:
    """Verify imported site sequences against an indexed reference FASTA in source orientation."""
    index_path, entries = _fasta_index(reference_path)
    reference_sha256 = _sha256(reference_path)
    index_sha256 = _sha256(index_path)
    try:
        handle = reference_path.open("rb")
    except OSError as exc:
        raise CandidatePanelError(f"cannot read reference FASTA {reference_path}: {exc}") from exc
    with handle:
        for nomination in nominations:
            location = (
                f"{nomination.get('source_file', reference_path.name)} "
                f"source_record_id={nomination.get('source_record_id')!r}"
            )
            metadata = nomination.get("source_metadata")
            if metadata is None and nomination.get("source_metadata_json"):
                try:
                    metadata = json.loads(nomination["source_metadata_json"])
                except json.JSONDecodeError:
                    metadata = None
            crisprme = metadata.get("crisprme") if isinstance(metadata, dict) else None
            if isinstance(crisprme, dict) and crisprme.get("canonical_interval_basis") == "alt_only":
                entry = entries.get(_canonical_chrom(nomination["chrom"]))
                if entry is None or nomination["end_0based_exclusive"] != nomination["start_0based"] + 1:
                    raise CandidatePanelError(f"{location}: ALT-only row requires a valid start anchor")
                _fetch_fasta_interval(
                    handle, entry, nomination["start_0based"],
                    nomination["end_0based_exclusive"], location,
                )
                nomination["reference_validation"] = {
                    "status": "reported_start_anchor_only",
                    "reference_fasta": reference_path.name,
                    "reference_sha256": reference_sha256,
                    "reference_index": index_path.name,
                    "reference_index_sha256": index_sha256,
                }
                continue
            sequence = nomination.get("site_sequence")
            strand = nomination.get("strand")
            if not sequence:
                nomination["reference_validation"] = {
                    "status": "not_evaluated_missing_site_sequence"
                }
                continue
            if not strand:
                nomination["reference_validation"] = {"status": "not_evaluated_missing_strand"}
                continue
            entry = entries.get(_canonical_chrom(nomination["chrom"]))
            if entry is None:
                raise CandidatePanelError(
                    f"{location}: chromosome {nomination['chrom']!r} is absent from the FASTA index"
                )
            observed = sequence.replace("-", "").upper()
            interval_length = nomination["end_0based_exclusive"] - nomination["start_0based"]
            if len(observed) != interval_length:
                raise CandidatePanelError(
                    f"{location}: ungapped site_sequence length does not equal the genomic interval"
                )
            reference = _fetch_fasta_interval(
                handle,
                entry,
                nomination["start_0based"],
                nomination["end_0based_exclusive"],
                location,
            )
            oriented = reference if strand == "+" else _reverse_complement(reference)
            if any(source != "N" and source != expected for source, expected in zip(observed, oriented)):
                raise CandidatePanelError(
                    f"{location}: site_sequence does not match the declared reference and strand"
                )
            nomination["reference_validation"] = {
                "status": "partial_unknown_sequence" if "N" in observed else "verified",
                "reference_fasta": reference_path.name,
                "reference_sha256": reference_sha256,
                "reference_index": index_path.name,
                "reference_index_sha256": index_sha256,
            }
    return {
        "file": reference_path.name,
        "sha256": reference_sha256,
        "index_file": index_path.name,
        "index_sha256": index_sha256,
    }


def _parse_site(row: dict, location: str) -> dict:
    start = _optional_int(row["start_0based"], "start_0based", location)
    end = _optional_int(row["end_0based_exclusive"], "end_0based_exclusive", location)
    assert start is not None and end is not None
    if end <= start:
        raise CandidatePanelError(f"{location}: end must be greater than start")
    strand = row.get("strand", "")
    if strand not in {"", "+", "-"}:
        raise CandidatePanelError(f"{location}: strand must be +, -, or blank")
    source_chrom = row["chrom"]
    site = {
        "guide_id": row["guide_id"],
        "genome_build": row["genome_build"],
        "chrom": _canonical_chrom(source_chrom),
        "start_0based": start,
        "end_0based_exclusive": end,
        "strand": strand or None,
    }
    if site["chrom"] != source_chrom:
        site["source_chrom"] = source_chrom
    return site


def _site_key(record: dict) -> tuple[str, str, str, int, int]:
    return (
        record["guide_id"],
        record["genome_build"],
        record["chrom"],
        record["start_0based"],
        record["end_0based_exclusive"],
    )


def _site_id(key: tuple[str, str, str, int, int]) -> str:
    encoded = "\x1f".join(str(part) for part in key).encode("utf-8")
    return "site-" + hashlib.sha256(encoded).hexdigest()[:16]


def _validate_crisprme_metadata(
    metadata: dict, site: dict, site_sequence_upper: str, location: str
) -> None:
    """Check the permitted CRISPRme+ fields without retaining sample identifiers."""
    if set(metadata) != {"crisprme"} or not isinstance(metadata["crisprme"], dict):
        raise CandidatePanelError(f"{location}: invalid CRISPRme source_metadata_json")
    block = metadata["crisprme"]
    expected = {
        "format_schema",
        "guide_alignment",
        "selected_alignment",
        "selected_origin",
        "canonical_interval_basis",
        "allele_states",
        "variant_maf_raw",
        "pam_creation_raw",
        "not_found_in_ref_raw",
    }
    if set(block) != expected or block["format_schema"] not in SUPPORTED_FORMAT_SCHEMAS:
        raise CandidatePanelError(f"{location}: invalid CRISPRme metadata fields")
    for field in ("guide_alignment", "selected_alignment"):
        value = block[field]
        if not isinstance(value, str) or not value or set(value.upper()) - set("ACGTN-"):
            raise CandidatePanelError(f"{location}: invalid CRISPRme {field}")
    if block["selected_origin"] not in {"ref", "alt"}:
        raise CandidatePanelError(f"{location}: invalid CRISPRme selected_origin")
    if block["canonical_interval_basis"] not in {"ref", "alt_only"}:
        raise CandidatePanelError(
            f"{location}: invalid CRISPRme canonical_interval_basis"
        )
    for field in ("variant_maf_raw", "pam_creation_raw", "not_found_in_ref_raw"):
        if block[field] is not None and not isinstance(block[field], str):
            raise CandidatePanelError(f"{location}: invalid CRISPRme {field}")
    states = block["allele_states"]
    if not isinstance(states, list) or not 1 <= len(states) <= 2:
        raise CandidatePanelError(f"{location}: invalid CRISPRme allele_states")
    alleles = set()
    selected = []
    canonical_alignment = None
    for state in states:
        state_fields = {
            "allele",
            "alignment",
            "effective_interval",
            "CFD",
            "variant_rsid_raw",
            "variant_maf_raw",
            "variant_info_raw",
            "origin_selected",
        }
        if (not isinstance(state, dict) or not state_fields <= set(state)
                or set(state) - state_fields - {"mapping_status", "genomic_forward_changes"}):
            raise CandidatePanelError(f"{location}: invalid CRISPRme allele state")
        allele = state["allele"]
        if allele not in {"ref", "alt"} or allele in alleles:
            raise CandidatePanelError(f"{location}: invalid CRISPRme allele label")
        alleles.add(allele)
        alignment = state["alignment"]
        if (
            not isinstance(alignment, str)
            or not alignment
            or set(alignment.upper()) - set("ACGTN-")
            or len(alignment) != len(block["guide_alignment"])
        ):
            raise CandidatePanelError(f"{location}: invalid CRISPRme allele alignment")
        interval = state["effective_interval"]
        if interval is None:
            if allele != "alt":
                raise CandidatePanelError(f"{location}: reference allele needs an interval")
        else:
            if not isinstance(interval, dict) or set(interval) != {
                "chrom",
                "start_0based",
                "end_0based_exclusive",
                "strand",
            }:
                raise CandidatePanelError(f"{location}: invalid CRISPRme allele interval")
            if (
                _canonical_chrom(str(interval["chrom"])) != site["chrom"]
                or interval["start_0based"] != site["start_0based"]
                or not isinstance(interval["end_0based_exclusive"], int)
                or interval["end_0based_exclusive"] <= interval["start_0based"]
                or interval["strand"] != site["strand"]
            ):
                raise CandidatePanelError(f"{location}: inconsistent CRISPRme allele interval")
        cfd = state["CFD"]
        if isinstance(cfd, bool) or not isinstance(cfd, (int, float)) or not math.isfinite(cfd):
            raise CandidatePanelError(f"{location}: invalid CRISPRme allele CFD")
        if state["variant_rsid_raw"] is not None and not isinstance(
            state["variant_rsid_raw"], str
        ):
            raise CandidatePanelError(f"{location}: invalid CRISPRme variant rsID")
        if state["variant_maf_raw"] is not None and not isinstance(
            state["variant_maf_raw"], str
        ):
            raise CandidatePanelError(f"{location}: invalid CRISPRme variant MAF")
        variant_info = state["variant_info_raw"]
        if not isinstance(variant_info, dict) or set(variant_info) != {
            "spacer_pam",
            "genome",
        }:
            raise CandidatePanelError(f"{location}: invalid CRISPRme variant info")
        if not isinstance(state["origin_selected"], bool):
            raise CandidatePanelError(f"{location}: invalid CRISPRme origin flag")
        if state["origin_selected"]:
            selected.append(state)
        if (
            block["canonical_interval_basis"] == "ref" and allele == "ref"
        ) or (
            block["canonical_interval_basis"] == "alt_only" and allele == "alt"
        ):
            canonical_alignment = alignment
            if (block["canonical_interval_basis"] == "ref"
                    and interval["end_0based_exclusive"] != site["end_0based_exclusive"]):
                raise CandidatePanelError(
                    f"{location}: CRISPRme canonical allele interval differs from site identity"
                )
    for state, mapped in zip(states, reference_mapped_states(states)):
        if "mapping_status" in state and (
            state["mapping_status"] != mapped["mapping_status"]
            or state["effective_interval"] != mapped["effective_interval"]
        ):
            raise CandidatePanelError(f"{location}: unsupported CRISPRme coordinate assertion")
        if (
            "genomic_forward_changes" in state
            and state["genomic_forward_changes"] != mapped["genomic_forward_changes"]
        ):
            raise CandidatePanelError(f"{location}: inconsistent genomic-forward allele changes")
    if (
        len(selected) != 1
        or selected[0]["allele"] != block["selected_origin"]
        or selected[0]["alignment"] != block["selected_alignment"]
        or canonical_alignment is None
        or (block["canonical_interval_basis"] == "ref"
            and canonical_alignment.upper() != site_sequence_upper)
        or (block["canonical_interval_basis"] == "alt_only"
            and (site_sequence_upper or site["end_0based_exclusive"] != site["start_0based"] + 1))
    ):
        raise CandidatePanelError(f"{location}: inconsistent CRISPRme selected/canonical state")


def load_guides(path: Path) -> list[dict]:
    """Load the guide and editor definitions referenced by every evidence row."""
    guides = []
    seen = set()
    for line_number, row in _read_tsv(path, GUIDE_FIELDS, GUIDE_REQUIRED):
        location = f"{path}:{line_number}"
        guide_id = row["guide_id"]
        if guide_id in seen:
            raise CandidatePanelError(f"{location}: duplicate guide_id {guide_id!r}")
        seen.add(guide_id)
        spacer = row["spacer_sequence"].upper()
        if not spacer or any(base not in "ACGT" for base in spacer):
            raise CandidatePanelError(f"{location}: spacer_sequence must contain only A/C/G/T")
        start = _optional_int(
            row["editing_window_start_1based"], "editing_window_start_1based", location
        )
        end = _optional_int(
            row["editing_window_end_1based_inclusive"],
            "editing_window_end_1based_inclusive",
            location,
        )
        assert start is not None and end is not None
        if start < 1 or end < start or end > len(spacer):
            raise CandidatePanelError(
                f"{location}: editing window must lie within the spacer and use 1-based "
                "inclusive positions"
            )
        locus_system = row["locus_system"]
        if locus_system not in {"HBG1_HBG2", "other"}:
            raise CandidatePanelError(
                f"{location}: locus_system must be HBG1_HBG2 or other"
            )
        pam_class = _validate_pam_pattern(row["pam_class"])
        guides.append(
            {
                "guide_id": guide_id,
                "spacer_sequence": spacer,
                "editor": row["editor"],
                "editor_class": row["editor_class"],
                "pam_class": pam_class,
                "editing_window_start_1based": start,
                "editing_window_end_1based_inclusive": end,
                "intended_target_description": row["intended_target_description"],
                "locus_system": locus_system,
            }
        )
    if not guides:
        raise CandidatePanelError("at least one guide definition is required")
    return sorted(guides, key=lambda item: item["guide_id"])


def load_intended_loci(path: Path) -> list[dict]:
    """Load the only coordinates permitted to carry the intended-target role."""
    loci = []
    seen_keys = set()
    seen_names = set()
    for line_number, row in _read_tsv(path, INTENDED_LOCUS_FIELDS, INTENDED_LOCUS_REQUIRED):
        location = f"{path}:{line_number}"
        locus = {**_parse_site(row, location), "locus_name": row["locus_name"]}
        if locus["strand"] is None:
            raise CandidatePanelError(f"{location}: intended locus requires strand")
        key = _site_key(locus)
        name_key = (locus["guide_id"], locus["locus_name"])
        if key in seen_keys:
            raise CandidatePanelError(f"{location}: duplicate intended-locus coordinate")
        if name_key in seen_names:
            raise CandidatePanelError(
                f"{location}: duplicate intended locus name {locus['locus_name']!r} "
                f"for guide {locus['guide_id']!r}"
            )
        seen_keys.add(key)
        seen_names.add(name_key)
        loci.append(locus)
    if not loci:
        raise CandidatePanelError("at least one intended locus is required")
    return sorted(loci, key=lambda item: (_site_key(item), item["locus_name"]))


def _validate_intended_loci(
    guides: list[dict], intended_loci: list[dict], profile: SearchProfile
) -> list[dict]:
    """Require intended coordinates and the correct spacing for paired HBG targets."""
    guide_ids = {guide["guide_id"] for guide in guides}
    locus_guide_ids = {locus["guide_id"] for locus in intended_loci}
    unknown = sorted(locus_guide_ids - guide_ids)
    missing = sorted(guide_ids - locus_guide_ids)
    if unknown:
        raise CandidatePanelError(
            "intended-locus manifest references unknown guide_id values: " + ", ".join(unknown)
        )
    if missing:
        raise CandidatePanelError(
            "every guide requires at least one intended locus; missing: " + ", ".join(missing)
        )
    wrong_reference = sorted(
        {locus["genome_build"] for locus in intended_loci if locus["genome_build"] != profile.reference_id}
    )
    if wrong_reference:
        raise CandidatePanelError(
            "intended loci must use profile reference_id "
            f"{profile.reference_id!r}; found: {', '.join(wrong_reference)}"
        )

    requirements = []
    by_guide: dict[str, list[dict]] = {}
    for locus in intended_loci:
        by_guide.setdefault(locus["guide_id"], []).append(locus)
    guide_by_id = {guide["guide_id"]: guide for guide in guides}
    for guide_id, loci in sorted(by_guide.items()):
        hbg = {locus["locus_name"]: locus for locus in loci if locus["locus_name"] in {"HBG1", "HBG2"}}
        declared_hbg = guide_by_id[guide_id]["locus_system"] == "HBG1_HBG2"
        if not hbg and not declared_hbg:
            continue
        if set(hbg) != {"HBG1", "HBG2"}:
            raise CandidatePanelError(
                f"guide {guide_id!r} declares an HBG intended locus but must declare both HBG1 and HBG2"
            )
        left, right = hbg["HBG1"], hbg["HBG2"]
        if left["chrom"] != right["chrom"]:
            raise CandidatePanelError(f"guide {guide_id!r}: HBG1/HBG2 intended loci must share a chromosome")
        if left["strand"] != right["strand"]:
            raise CandidatePanelError(f"guide {guide_id!r}: HBG1/HBG2 intended loci must share a strand")
        if left["end_0based_exclusive"] - left["start_0based"] != (
            right["end_0based_exclusive"] - right["start_0based"]
        ):
            raise CandidatePanelError(f"guide {guide_id!r}: HBG1/HBG2 intended intervals must have equal length")
        spacing = right["start_0based"] - left["start_0based"]
        if spacing != anchors.PARALOG_SPACING_BP:
            raise CandidatePanelError(
                f"guide {guide_id!r}: HBG1/HBG2 intended intervals are {spacing} bp apart; "
                f"expected {anchors.PARALOG_SPACING_BP}"
            )
        requirements.append(
            {
                "requirement_id": f"{guide_id}:hbg_duplicated_locus",
                "guide_id": guide_id,
                "locus_system": "HBG1/HBG2 duplicated promoter",
                "paralog_spacing_bp": anchors.PARALOG_SPACING_BP,
                "intended_locus_names": ["HBG1", "HBG2"],
                "concern": "co-editing can produce an approximately 4.9-kb intergenic deletion or rearrangement",
                "required_confirmation": [
                    {
                        "validates": "hbg_paralog_rearrangement",
                        "assay": "long-read or spanning-amplicon sequencing across HBG1-HBG2",
                        "status": "required",
                    },
                    {
                        "validates": "hbg_paralog_cotargeting",
                        "assay": "locus-resolving amplicon sequencing that distinguishes HBG1 and HBG2",
                        "status": "required",
                    },
                ],
            }
        )
    return requirements


def load_nominations(paths: list[Path], profile: SearchProfile) -> list[dict]:
    nominations = []
    seen_source_ids: dict[tuple[str, str, str, str], str] = {}
    seen_rows: dict[str, str] = {}
    for path in paths:
        for line_number, row in _read_tsv(path, NOMINATION_FIELDS, NOMINATION_REQUIRED):
            location = f"{path}:{line_number}"
            site = _parse_site(row, location)
            kind = row["nomination_kind"].lower()
            if kind not in {"computational", "experimental"}:
                raise CandidatePanelError(
                    f"{location}: nomination_kind must be computational or experimental"
                )
            role = row.get("site_role", "").lower()
            if role not in {"", "candidate", "intended"}:
                raise CandidatePanelError(
                    f"{location}: site_role must be candidate, intended, or blank"
                )

            mismatch = _optional_int(row.get("mismatch_count", ""), "mismatch_count", location)
            dna_bulge = _optional_int(row.get("dna_bulge", ""), "dna_bulge", location)
            rna_bulge = _optional_int(row.get("rna_bulge", ""), "rna_bulge", location)
            score = _optional_float(row.get("score_value", ""), "score_value", location)
            editable = _optional_bool(
                row.get("editable_base_in_window", ""), "editable_base_in_window", location
            )
            score_name = row.get("score_name", "")
            if (score is None) != (score_name == ""):
                raise CandidatePanelError(
                    f"{location}: score_name and score_value must be supplied together"
                )
            sequence = row.get("site_sequence", "").upper()
            if sequence and any(base not in "ACGTN-" for base in sequence):
                raise CandidatePanelError(
                    f"{location}: site_sequence must contain only A/C/G/T/N and alignment gaps (-)"
                )
            metadata_text = row.get("source_metadata_json", "")
            try:
                source_metadata = json.loads(metadata_text) if metadata_text else {}
            except json.JSONDecodeError as exc:
                raise CandidatePanelError(
                    f"{location}: source_metadata_json must be a JSON object"
                ) from exc
            if not isinstance(source_metadata, dict):
                raise CandidatePanelError(
                    f"{location}: source_metadata_json must be a JSON object"
                )
            if row["source_method"].casefold() == "crisprme":
                _validate_crisprme_metadata(source_metadata, site, sequence, location)

            nomination = {
                **site,
                "nomination_kind": kind,
                "source_method": row["source_method"],
                "source_version_or_citation": row["source_version_or_citation"],
                "source_file": path.name,
                "site_role": role or None,
                "site_sequence": sequence or None,
                "pam": row.get("pam", "").upper() or None,
                "mismatch_count": mismatch,
                "dna_bulge": 0 if kind == "computational" and dna_bulge is None else dna_bulge,
                "rna_bulge": 0 if kind == "computational" and rna_bulge is None else rna_bulge,
                "score_name": score_name or None,
                "score_value": score,
                "editable_base_in_window": editable,
                "population_variant_id": row.get("population_variant_id", "") or None,
                "ch_driver_gene": row.get("ch_driver_gene", "") or None,
                "source_record_id": row.get("source_record_id", "") or None,
                "search_profile": row.get("search_profile", "") or None,
                "source_metadata": source_metadata,
                "annotation_validation": {
                    "population_variant_id": (
                        "imported_not_verified"
                        if row.get("population_variant_id", "")
                        else "not_present"
                    ),
                    "ch_driver_gene": (
                        "imported_not_verified" if row.get("ch_driver_gene", "") else "not_present"
                    ),
                },
            }
            if kind == "computational":
                profile.validate_computational_nomination(nomination, location)
            source_key = (
                kind,
                nomination["source_method"],
                nomination["source_version_or_citation"],
                nomination["source_record_id"],
            )
            if source_key in seen_source_ids:
                raise CandidatePanelError(
                    f"{location}: duplicate nomination source_record_id; first used at "
                    f"{seen_source_ids[source_key]}"
                )
            fingerprint = json.dumps(nomination, sort_keys=True, separators=(",", ":"))
            if fingerprint in seen_rows:
                raise CandidatePanelError(
                    f"{location}: exact duplicate nomination row; first used at "
                    f"{seen_rows[fingerprint]}"
                )
            seen_source_ids[source_key] = location
            seen_rows[fingerprint] = location
            nominations.append(nomination)
    if not nominations:
        raise CandidatePanelError("at least one nomination row is required")
    return nominations


def load_measurements(path: Path | None) -> list[dict]:
    if path is None:
        return []
    measurements = []
    seen_source_ids: dict[tuple[str, str, str], str] = {}
    seen_rows: dict[str, str] = {}
    for line_number, row in _read_tsv(path, MEASUREMENT_FIELDS, MEASUREMENT_REQUIRED):
        location = f"{path}:{line_number}"
        site = _parse_site(row, location)
        status = row["status"].lower()
        if status not in {"edited", "not_detected", "inconclusive"}:
            raise CandidatePanelError(
                f"{location}: status must be edited, not_detected, or inconclusive"
            )
        value = _optional_float(row.get("value", ""), "value", location)
        unit = row.get("unit", "")
        if (value is None) != (unit == ""):
            raise CandidatePanelError(f"{location}: value and unit must be supplied together")
        lod_value = _optional_float(
            row.get("limit_of_detection_value", ""), "limit_of_detection_value", location
        )
        lod_unit = row.get("limit_of_detection_unit", "")
        if (lod_value is None) != (lod_unit == ""):
            raise CandidatePanelError(
                f"{location}: limit_of_detection_value and limit_of_detection_unit must be supplied together"
            )
        if lod_value is not None and lod_value <= 0:
            raise CandidatePanelError(f"{location}: limit_of_detection_value must be greater than zero")
        if status == "not_detected" and lod_value is None:
            raise CandidatePanelError(
                f"{location}: not_detected measurement requires numeric limit_of_detection_value and unit"
            )
        read_depth = _optional_int(row.get("read_depth", ""), "read_depth", location)
        p_value = _optional_float(row.get("p_value", ""), "p_value", location)
        adjusted_p_value = _optional_float(
            row.get("adjusted_p_value", ""), "adjusted_p_value", location
        )
        for field, parsed in (("p_value", p_value), ("adjusted_p_value", adjusted_p_value)):
            if parsed is not None and not 0 <= parsed <= 1:
                raise CandidatePanelError(f"{location}: {field} must lie inside [0, 1]")
        qc_status = row.get("qc_status", "").lower()
        if qc_status not in {"", "passed", "failed", "not_reported"}:
            raise CandidatePanelError(
                f"{location}: qc_status must be passed, failed, not_reported, or blank"
            )
        context = {field: row.get(field, "") or None for field in MEASUREMENT_CONTEXT_FIELDS}
        for field, allowed in (
            ("measurement_kind", {"editing", "consequence", "not_reported"}),
            ("allele", {"ref", "alt", "unknown"}),
            ("allele_resolution", {"resolved", "mixed", "not_reported"}),
            ("phase_status", {"phased", "unphased", "not_reported"}),
        ):
            if context[field] is not None:
                context[field] = context[field].lower()
                if context[field] not in allowed:
                    raise CandidatePanelError(f"{location}: invalid {field}: {context[field]!r}")
        if context["allele_sequence"] is not None:
            context["allele_sequence"] = context["allele_sequence"].upper()
            if set(context["allele_sequence"]) - set("ACGTN"):
                raise CandidatePanelError(f"{location}: allele_sequence must be genomic-forward A/C/G/T/N")
        measurement = {
            **site,
            **context,
            "assay": row["assay"],
            "source_version_or_citation": row["source_version_or_citation"],
            "status": status,
            "value": value,
            "unit": unit or None,
            "limit_of_detection": (
                {"value": lod_value, "unit": lod_unit} if lod_value is not None else None
            ),
            "biological_context": row["biological_context"],
            "source_record_id": row["source_record_id"],
            "replicate_id": row.get("replicate_id", "") or None,
            "treated_group": row.get("treated_group", "") or None,
            "control_group": row.get("control_group", "") or None,
            "read_depth": read_depth,
            "qc_status": qc_status or None,
            "statistical_test": row.get("statistical_test", "") or None,
            "p_value": p_value,
            "adjusted_p_value": adjusted_p_value,
            "source_file": path.name,
        }
        source_key = (
            measurement["assay"],
            measurement["source_version_or_citation"],
            measurement["source_record_id"],
        )
        if source_key in seen_source_ids:
            raise CandidatePanelError(
                f"{location}: duplicate measurement source_record_id; first used at "
                f"{seen_source_ids[source_key]}"
            )
        fingerprint = json.dumps(measurement, sort_keys=True, separators=(",", ":"))
        if fingerprint in seen_rows:
            raise CandidatePanelError(
                f"{location}: exact duplicate measurement row; first used at "
                f"{seen_rows[fingerprint]}"
            )
        seen_source_ids[source_key] = location
        seen_rows[fingerprint] = location
        measurements.append(measurement)
    return measurements


def load_decisions(path: Path | None) -> dict[tuple[str, str, str, int, int], dict]:
    if path is None:
        return {}
    decisions = {}
    for line_number, row in _read_tsv(path, DECISION_FIELDS, DECISION_REQUIRED):
        location = f"{path}:{line_number}"
        site = _parse_site(row, location)
        decision = row["decision"].lower()
        if decision not in {"include", "exclude", "undecided"}:
            raise CandidatePanelError(
                f"{location}: decision must be include, exclude, or undecided"
            )
        reason = row["reason"].lower()
        if reason not in DECISION_REASONS:
            raise CandidatePanelError(
                f"{location}: decision reason must be one of "
                + ", ".join(sorted(DECISION_REASONS))
            )
        note = row.get("note", "")
        if len(note) > 140 or "\t" in note or "\n" in note or "\r" in note:
            raise CandidatePanelError(
                f"{location}: decision note must be at most 140 characters with no tabs or newlines"
            )
        related_site_id = row.get("related_site_id", "") or None
        if reason == "redundant_with":
            if related_site_id is None:
                raise CandidatePanelError(
                    f"{location}: redundant_with requires related_site_id"
                )
        elif related_site_id is not None:
            raise CandidatePanelError(
                f"{location}: related_site_id is allowed only for redundant_with"
            )
        key = _site_key(site)
        if key in decisions:
            raise CandidatePanelError(f"{location}: duplicate panel decision for the same site")
        decisions[key] = {
            "status": decision,
            "reason": reason,
            "note": note or None,
            "related_site_id": related_site_id,
            "source_file": path.name,
        }
    return decisions


def load_source_provenance(
    paths: list[Path] | None,
    nomination_paths: list[Path],
    nominations: list[dict],
    profile: SearchProfile,
    profile_path: Path,
    guide_path: Path,
    intended_locus_path: Path,
    reference_record: dict | None,
) -> list[dict]:
    """Check the companion source-details JSON for every converted candidate file."""
    if not paths:
        raise CandidatePanelError(
            "at least one --source-provenance file is required for nomination inputs"
        )
    expected_by_name = {path.name: path for path in nomination_paths}
    if len(expected_by_name) != len(nomination_paths):
        raise CandidatePanelError("nomination input basenames must be unique")
    nomination_counts = Counter(item["source_file"] for item in nominations)
    by_name = {}
    required_top = {
        "record_type",
        "schema_version",
        "normalized_output",
        "native_input",
        "source_run",
        "adapter",
        "guide_manifest",
        "intended_locus_manifest",
        "search_profile",
        "generated_at_utc",
        "reference_fasta",
    }
    for path in paths:
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise CandidatePanelError(f"cannot read source provenance {path}: {exc}") from exc
        if not isinstance(payload, dict):
            raise CandidatePanelError(f"{path}: source provenance must be a JSON object")
        schema_version = payload.get("schema_version")
        if schema_version == "1":
            expected_fields = required_top | (
                {"source_declarations"} if "source_declarations" in payload else set()
            )
        elif schema_version == "2":
            expected_fields = required_top | {
                "source_declarations",
                "import_summary",
                "sheets_ignored",
                "privacy_filtering_applied",
            }
            if "original_source_declaration" in payload:
                expected_fields.add("original_source_declaration")
        else:
            expected_fields = set()
        if not expected_fields or set(payload) != expected_fields:
            raise CandidatePanelError(
                f"{path}: invalid source-provenance fields; expected exactly "
                + ", ".join(sorted(expected_fields or required_top))
            )
        if payload["record_type"] != "nomination_source_provenance":
            raise CandidatePanelError(f"{path}: unsupported source-provenance record or schema")
        normalized = payload["normalized_output"]
        native = payload["native_input"]
        source_run = payload["source_run"]
        if not isinstance(normalized, dict) or set(normalized) != {"file", "sha256", "row_count"}:
            raise CandidatePanelError(f"{path}: invalid normalized_output provenance")
        normalized_name = Path(_required_text(normalized["file"], "normalized output file")).name
        if normalized_name != normalized["file"]:
            raise CandidatePanelError(f"{path}: provenance files must use basenames, not paths")
        nomination_path = expected_by_name.get(normalized_name)
        if nomination_path is None:
            raise CandidatePanelError(
                f"{path}: normalized output {normalized_name!r} is not a --nominations input"
            )
        if normalized_name in by_name:
            raise CandidatePanelError(
                f"{path}: duplicate source provenance for {normalized_name!r}"
            )
        if normalized.get("sha256") != _sha256(nomination_path):
            raise CandidatePanelError(
                f"{path}: normalized output checksum does not match {normalized_name}"
            )
        if normalized.get("row_count") != nomination_counts[normalized_name]:
            raise CandidatePanelError(
                f"{path}: normalized output row_count does not match parsed nominations"
            )
        if schema_version == "1" and "source_declarations" in payload:
            declarations = payload["source_declarations"]
            if not isinstance(declarations, dict) or set(declarations) != {"crispor"}:
                raise CandidatePanelError(f"{path}: invalid source_declarations provenance")
            try:
                declarations["crispor"] = validate_crispor_source_declaration(
                    declarations["crispor"]
                )
            except CandidatePanelError as exc:
                raise CandidatePanelError(
                    f"{path}: invalid CRISPOR source declaration: {exc}"
                ) from exc
        if schema_version == "2":
            declarations = payload["source_declarations"]
            if not isinstance(declarations, dict) or set(declarations) != {"crisprme"}:
                raise CandidatePanelError(f"{path}: invalid source_declarations provenance")
            try:
                declaration = validate_import_declaration(declarations["crisprme"])
            except CrisprmeImportError as exc:
                raise CandidatePanelError(
                    f"{path}: invalid CRISPRme source declaration: {exc}"
                ) from exc
            if declaration["native_input"] != native:
                raise CandidatePanelError(f"{path}: invalid CRISPRme source declaration")
            declarations["crisprme"] = declaration
            summary = payload["import_summary"]
            expected_summary = {
                "input_data_rows",
                "emitted_nomination_rows",
                "skipped_rows",
                "skipped_by_reason",
            }
            if not isinstance(summary, dict) or set(summary) != expected_summary:
                raise CandidatePanelError(f"{path}: invalid CRISPRme import_summary")
            counts = [
                summary["input_data_rows"],
                summary["emitted_nomination_rows"],
                summary["skipped_rows"],
            ]
            if any(isinstance(value, bool) or not isinstance(value, int) or value < 0 for value in counts):
                raise CandidatePanelError(f"{path}: invalid CRISPRme import counts")
            reason_counts = summary["skipped_by_reason"]
            if (
                not isinstance(reason_counts, dict)
                or any(
                    not isinstance(reason, str)
                    or reason not in SKIP_REASONS
                    or isinstance(count, bool)
                    or not isinstance(count, int)
                    or count < 1
                    for reason, count in reason_counts.items()
                )
                or sum(reason_counts.values()) != summary["skipped_rows"]
                or summary["input_data_rows"]
                != summary["emitted_nomination_rows"] + summary["skipped_rows"]
                or summary["emitted_nomination_rows"] != normalized["row_count"]
            ):
                raise CandidatePanelError(f"{path}: inconsistent CRISPRme import_summary")
            effective = declarations["crisprme"]
            if summary["skipped_rows"] and (
                effective["retention"]["completeness"] == "complete"
                or effective["retention"]["all_candidate_coordinates_retained"]
            ):
                raise CandidatePanelError(f"{path}: incomplete import cannot claim completeness")
            if "original_source_declaration" in payload:
                original = validate_import_declaration(payload["original_source_declaration"])
                if {**original, "retention": effective["retention"]} != effective:
                    raise CandidatePanelError(f"{path}: original and effective declarations disagree")
            ignored = payload["sheets_ignored"]
            if (
                not isinstance(ignored, dict)
                or set(ignored) != {"count", "categories"}
                or isinstance(ignored["count"], bool)
                or not isinstance(ignored["count"], int)
                or ignored["count"] < 0
                or not isinstance(ignored["categories"], list)
                or ignored["categories"] != sorted(set(ignored["categories"]))
                or any(
                    category not in {"altMerge", "private_targets", "other"}
                    for category in ignored["categories"]
                )
            ):
                raise CandidatePanelError(f"{path}: invalid sheets_ignored provenance")
            if payload["privacy_filtering_applied"] is not True:
                raise CandidatePanelError(
                    f"{path}: CRISPRme provenance must record privacy_filtering_applied=true"
                )
        for block_name, block in (
            ("native_input", native),
            ("guide_manifest", payload["guide_manifest"]),
            ("intended_locus_manifest", payload["intended_locus_manifest"]),
        ):
            if not isinstance(block, dict) or set(block) != {"file", "sha256"}:
                raise CandidatePanelError(f"{path}: invalid {block_name} provenance")
            filename = _required_text(block["file"], f"{block_name} file")
            if Path(filename).name != filename:
                raise CandidatePanelError(f"{path}: provenance files must use basenames, not paths")
            if not isinstance(block["sha256"], str) or not re_full_sha256(block["sha256"]):
                raise CandidatePanelError(f"{path}: {block_name} sha256 is invalid")
        if payload["guide_manifest"] != {
            "file": guide_path.name,
            "sha256": _sha256(guide_path),
        }:
            raise CandidatePanelError(f"{path}: guide-manifest provenance does not match this run")
        if payload["intended_locus_manifest"] != {
            "file": intended_locus_path.name,
            "sha256": _sha256(intended_locus_path),
        }:
            raise CandidatePanelError(
                f"{path}: intended-locus provenance does not match this run"
            )
        profile_block = payload["search_profile"]
        if profile_block is not None:
            if profile_block != {"file": profile_path.name, "sha256": _sha256(profile_path)}:
                raise CandidatePanelError(f"{path}: search-profile provenance does not match this run")
        reference_block = payload["reference_fasta"]
        if reference_record is None:
            if reference_block is not None:
                raise CandidatePanelError(
                    f"{path}: synthetic run must not claim an unchecked reference FASTA"
                )
        else:
            expected_reference_block = {
                "file": reference_record["file"],
                "sha256": reference_record["sha256"],
                "index_file": reference_record["index_file"],
                "index_sha256": reference_record["index_sha256"],
            }
            if reference_block != expected_reference_block:
                raise CandidatePanelError(
                    f"{path}: reference FASTA provenance does not match this run"
                )
        if not isinstance(source_run, dict) or set(source_run) != {
            "method",
            "version_or_citation",
            "command",
            "reference_id",
            "reference_sha256",
            "reference_checksum_status",
        }:
            raise CandidatePanelError(f"{path}: invalid source_run provenance")
        for field in ("method", "version_or_citation", "command", "reference_id"):
            _required_text(source_run[field], f"source_run {field}")
        source_rows = [
            item for item in nominations if item["source_file"] == normalized_name
        ]
        if profile.reference_id != "SYNTHETIC":
            methods = {item["source_method"] for item in source_rows}
            versions = {item["source_version_or_citation"] for item in source_rows}
            if methods != {source_run["method"]} or versions != {
                source_run["version_or_citation"]
            }:
                raise CandidatePanelError(
                    f"{path}: source_run method/version does not match its normalized rows"
                )
        if source_run["reference_id"] != profile.reference_id:
            raise CandidatePanelError(f"{path}: source_run reference_id does not match the profile")
        checksum_status = source_run["reference_checksum_status"]
        if checksum_status == "recorded":
            if not isinstance(source_run["reference_sha256"], str) or not re_full_sha256(
                source_run["reference_sha256"]
            ):
                raise CandidatePanelError(f"{path}: source_run reference_sha256 is invalid")
        elif checksum_status == "not_applicable_synthetic":
            if source_run["reference_sha256"] is not None or profile.reference_id != "SYNTHETIC":
                raise CandidatePanelError(
                    f"{path}: a missing reference checksum is allowed only for the synthetic fixture"
                )
        else:
            raise CandidatePanelError(f"{path}: invalid reference_checksum_status")
        adapter = payload["adapter"]
        if not isinstance(adapter, dict) or set(adapter) != {
            "mode",
            "base_edit_safety_version",
        }:
            raise CandidatePanelError(f"{path}: invalid adapter provenance")
        if schema_version == "2" and adapter["mode"] != "crisprme":
            raise CandidatePanelError(f"{path}: source provenance v2 is reserved for CRISPRme")
        by_name[normalized_name] = payload
    missing = sorted(set(expected_by_name) - set(by_name))
    if missing:
        raise CandidatePanelError(
            "missing source provenance for nomination inputs: " + ", ".join(missing)
        )
    return [by_name[path.name] for path in nomination_paths]


def re_full_sha256(value: str) -> bool:
    """Return whether a value is exactly one lowercase- or uppercase-hex SHA-256 digest."""
    return len(value) == 64 and all(character in "0123456789abcdefABCDEF" for character in value)


def _validate_nomination_content(
    nominations: list[dict], guides: list[dict], intended_loci: list[dict]
) -> None:
    """Recheck target roles and sequence-derived fields in converted candidate rows."""
    guide_by_id = {guide["guide_id"]: guide for guide in guides}
    intended_by_key = {_site_key(locus): locus for locus in intended_loci}
    edit_base = {"ABE": "A", "CBE": "C"}
    for nomination in nominations:
        key = _site_key(nomination)
        intended = intended_by_key.get(key)
        expected_role = "intended" if intended else "candidate"
        if nomination["site_role"] is not None and nomination["site_role"] != expected_role:
            locus_text = f" ({intended['locus_name']})" if intended else ""
            raise CandidatePanelError(
                f"{nomination['source_file']} source_record_id={nomination['source_record_id']!r}: "
                f"site_role must be {expected_role!r}{locus_text} according to the intended-locus manifest"
            )
        nomination["site_role"] = expected_role
        if intended and nomination.get("strand") and nomination["strand"] != intended["strand"]:
            raise CandidatePanelError(
                f"{nomination['source_file']} source_record_id={nomination['source_record_id']!r}: "
                "strand conflicts with the intended-locus manifest"
            )

        if nomination["nomination_kind"] != "computational":
            nomination["sequence_assertion_status"] = (
                "not_evaluated_experimental_nomination"
            )
            continue
        guide = guide_by_id.get(nomination["guide_id"])
        if guide is None:
            continue  # the unknown-guide guard reports this with the complete set later
        if not nomination.get("strand"):
            raise CandidatePanelError(
                f"{nomination['source_file']} source_record_id={nomination['source_record_id']!r}: "
                "computational nomination requires strand for reference verification"
            )
        if nomination["source_method"].casefold() == "crisprme":
            crisprme = nomination["source_metadata"]["crisprme"]
            guide_alignment = crisprme["guide_alignment"].replace("-", "").upper()
            if guide_alignment[: len(guide["spacer_sequence"])] != guide["spacer_sequence"]:
                raise CandidatePanelError(
                    f"{nomination['source_file']} "
                    f"source_record_id={nomination['source_record_id']!r}: "
                    "CRISPRme guide alignment conflicts with the guide manifest"
                )
            if nomination.get("editable_base_in_window") is not None:
                raise CandidatePanelError(
                    f"{nomination['source_file']} "
                    f"source_record_id={nomination['source_record_id']!r}: "
                    "CRISPRme editability must be derived by the ledger"
                )
            nomination["sequence_assertion_status"] = (
                "verified_crisprme_import_contract"
            )
            continue
        if (nomination.get("dna_bulge") or 0) or (nomination.get("rna_bulge") or 0):
            if nomination.get("editable_base_in_window") is not None:
                raise CandidatePanelError(
                    f"{nomination['source_file']} source_record_id={nomination['source_record_id']!r}: "
                    "editable_base_in_window must be blank for a gapped alignment"
                )
            nomination["sequence_assertion_status"] = (
                "not_evaluated_gapped_alignment_requires_alignment_aware_adapter"
            )
            continue
        sequence = nomination.get("site_sequence")
        pam = nomination.get("pam") or ""
        if not sequence:
            raise CandidatePanelError(
                f"{nomination['source_file']} source_record_id={nomination['source_record_id']!r}: "
                "ungapped computational nomination requires site_sequence"
            )
        if set(sequence) - set("ACGT"):
            raise CandidatePanelError(
                f"{nomination['source_file']} source_record_id={nomination['source_record_id']!r}: "
                "ungapped computational site_sequence must contain concrete A/C/G/T bases"
            )
        expected_length = len(guide["spacer_sequence"]) + len(pam)
        if len(sequence) != expected_length or not sequence.endswith(pam):
            raise CandidatePanelError(
                f"{nomination['source_file']} source_record_id={nomination['source_record_id']!r}: "
                "site_sequence length/PAM does not match the guide and observed PAM"
            )
        spacer = sequence[: len(guide["spacer_sequence"])]
        observed_mismatches = sum(
            left != right for left, right in zip(guide["spacer_sequence"], spacer)
        )
        if nomination["mismatch_count"] != observed_mismatches:
            raise CandidatePanelError(
                f"{nomination['source_file']} source_record_id={nomination['source_record_id']!r}: "
                f"mismatch_count is {nomination['mismatch_count']}, but sequence comparison gives "
                f"{observed_mismatches}"
            )
        source_base = edit_base.get(guide["editor_class"].upper())
        if source_base:
            start = guide["editing_window_start_1based"] - 1
            end = guide["editing_window_end_1based_inclusive"]
            observed_editable = source_base in spacer[start:end]
            asserted = nomination.get("editable_base_in_window")
            if asserted is not None and asserted != observed_editable:
                raise CandidatePanelError(
                    f"{nomination['source_file']} source_record_id={nomination['source_record_id']!r}: "
                    "editable_base_in_window conflicts with the guide, editor, and site sequence"
                )
            nomination["editable_base_in_window"] = observed_editable
        nomination["sequence_assertion_status"] = "verified_ungapped_against_guide"


def _validate_measurement_intended_strands(
    measurements: list[dict], intended_loci: list[dict]
) -> None:
    intended_by_key = {_site_key(locus): locus for locus in intended_loci}
    for measurement in measurements:
        intended = intended_by_key.get(_site_key(measurement))
        if intended and measurement.get("strand") and measurement["strand"] != intended["strand"]:
            raise CandidatePanelError(
                f"{measurement['source_file']} source_record_id={measurement['source_record_id']!r}: "
                "measurement strand conflicts with the intended-locus manifest"
            )


def _new_site(record: dict) -> dict:
    key = _site_key(record)
    return {
        "site_id": _site_id(key),
        "guide_id": record["guide_id"],
        "genome_build": record["genome_build"],
        "chrom": record["chrom"],
        "start_0based": record["start_0based"],
        "end_0based_exclusive": record["end_0based_exclusive"],
        "strands": set(),
        "site_role": None,
        "intended_locus": None,
        "identity_warnings": [],
        "locus_requirement_ids": [],
        "required_confirmation": [],
        "nominations": [],
        "measurements": [],
    }


def _measurement_summary(measurements: list[dict]) -> str:
    if not measurements:
        return "not_measured"
    statuses = {item["status"] for item in measurements if item.get("qc_status") == "passed"}
    if not statuses:
        return "unusable"
    if len(statuses) > 1:
        return "mixed"
    return next(iter(statuses))


def _evidence_state(site: dict) -> str:
    if not site["nominations"]:
        if site["site_role"] == "intended":
            return (
                "intended_control_measurement_without_nomination"
                if site["measurements"]
                else "intended_control_not_nominated"
            )
        return "measurement_without_nomination"
    if site["measurements"]:
        summary = site["measurement_summary"]
        if summary in {"edited", "not_detected", "mixed"}:
            return f"measured_{summary}"
        return {"not_measured": "measurement_does_not_cover_scenario",
                "partial_coverage": "measurement_partial_coverage",
                "coverage_unresolved": "measurement_coverage_unresolved",
                "unusable": "measurement_unusable"}[summary]
    kinds = {item["nomination_kind"] for item in site["nominations"]}
    if "experimental" in kinds:
        return "experimentally_nominated_unmeasured"
    if site["site_role"] == "intended":
        return "intended_control_unmeasured"
    return "computational_candidate_unmeasured"


def _review_group(site: dict) -> tuple[str, list[str]]:
    reasons = []
    nominations = site["nominations"]
    if site["site_role"] == "intended":
        reasons.append("declared intended target/control")
        return "intended_control", reasons
    if site["measurements"]:
        if not nominations:
            reasons.append("measurement has no matching nomination and needs reconciliation")
            return "measurement_reconciliation", reasons
        if any(item.get("qc_status") == "failed" for item in site["measurements"]):
            reasons.append("at least one imported measurement has failed QC")
            return "measurement_qc_attention", reasons
        if site["measurement_summary"] == "mixed":
            reasons.append("recorded assays or cell contexts have different results")
            return "measurement_reconciliation", reasons
        if site["measurement_summary"] not in {"edited", "not_detected"}:
            reasons.append("measurement coverage is incomplete, unresolved or unusable")
            return "measurement_reconciliation", reasons
        reasons.append("published or imported measurement is available")
        return "measured_site", reasons
    if any(item["nomination_kind"] == "experimental" for item in nominations):
        reasons.append("experimentally nominated but not measured")
        return "experimental_nomination", reasons
    if any(item.get("ch_driver_gene") for item in nominations):
        reasons.append("imported clonal-hematopoiesis driver-overlap label requires review")
        return "annotation_flagged_candidate", reasons
    if any(item.get("population_variant_id") for item in nominations):
        reasons.append("population-variant-aware candidate")
    if any((item.get("dna_bulge") or 0) > 0 or (item.get("rna_bulge") or 0) > 0 for item in nominations):
        reasons.append("gapped or bulged candidate")
    if reasons:
        return "expanded_search_candidate", reasons
    mismatches = [item["mismatch_count"] for item in nominations if item["mismatch_count"] is not None]
    if mismatches:
        lowest = min(mismatches)
        reasons.append(f"lowest imported mismatch count: {lowest}")
        reasons.append("mismatch band is a review queue, not a biological-risk estimate")
        if lowest <= 2:
            return "sequence_match_0_to_2_mismatches", reasons
        if lowest <= 4:
            return "sequence_match_3_to_4_mismatches", reasons
        return "sequence_match_5_plus_mismatches", reasons
    else:
        reasons.append("computational nomination without a reported mismatch count")
    return "computational_nomination_unclassified", reasons


def _computational_source_bucket(site: dict) -> str | None:
    methods = {
        item["source_method"].casefold()
        for item in site["nominations"]
        if item["nomination_kind"] == "computational"
    }
    has_crispor = "crispor" in methods
    has_crisprme = "crisprme" in methods
    if has_crispor and has_crisprme:
        return "both"
    if has_crisprme:
        return "crisprme_only"
    if has_crispor:
        return "crispor_only"
    return None


def _crisprme_metadata_for_site(site: dict) -> dict | None:
    blocks = [
        item["source_metadata"]["crisprme"]
        for item in site["nominations"]
        if item["source_method"].casefold() == "crisprme"
    ]
    if not blocks:
        return None
    fingerprints = {
        json.dumps(block, sort_keys=True, separators=(",", ":")) for block in blocks
    }
    return blocks[0] if len(fingerprints) == 1 else None


def _mismatches_for_alignment(guide_alignment: str, target_alignment: str) -> int:
    guide_position = 0
    mismatches = 0
    for guide_base, target_base in zip(guide_alignment.upper(), target_alignment.upper()):
        if guide_base != "-":
            guide_position += 1
        if (
            1 <= guide_position <= 20
            and guide_base != "-"
            and target_base != "-"
            and guide_base != target_base
        ):
            mismatches += 1
    return mismatches


def _pam_for_alignment(guide_alignment: str, target_alignment: str) -> str:
    guide_position = 0
    bases = []
    for guide_base, target_base in zip(guide_alignment.upper(), target_alignment.upper()):
        if guide_base != "-":
            guide_position += 1
        if guide_position > 20 and guide_base != "-" and target_base != "-":
            bases.append(target_base)
    return "".join(bases)


def _allele_pam_change(block: dict, guide: dict | None) -> str:
    states = {state["allele"]: state for state in block["allele_states"]}
    if guide is None or not {"ref", "alt"} <= set(states):
        return "not_assessable"
    ref = states["ref"]["alignment"].upper()
    alt = states["alt"]["alignment"].upper()
    if len(ref) != len(alt) or any((a == "-") != (b == "-") for a, b in zip(ref, alt)):
        return "unsupported_indel_alignment"
    ref_pam = _pam_for_alignment(block["guide_alignment"], ref)
    alt_pam = _pam_for_alignment(block["guide_alignment"], alt)
    if not ref_pam or not alt_pam or "N" in ref_pam + alt_pam:
        return "unknown_sequence"
    ref_valid = _pam_matches(ref_pam, guide["pam_class"])
    alt_valid = _pam_matches(alt_pam, guide["pam_class"])
    if not ref_valid and alt_valid:
        return "created"
    if ref_valid and not alt_valid:
        return "destroyed"
    return "retained" if ref_valid else "neither_matches"


def _why_crispor_missed(site: dict, guide: dict | None) -> list[str]:
    if site["source_bucket"] != "crisprme_only":
        return []
    configurations = [
        item for item in site["nominations"]
        if item["source_method"].casefold() == "crisprme"
    ]
    if len(configurations) > 1:
        reasons = set()
        for item in configurations:
            reasons.update(_why_crispor_missed({**site, "nominations": [item]}, guide))
        return sorted(reasons) or ["unexplained"]
    reasons = set()
    nominations = [
        item
        for item in site["nominations"]
        if item["source_method"].casefold() == "crisprme"
    ]
    metadata = _crisprme_metadata_for_site(site)
    for nomination in nominations:
        if (nomination.get("dna_bulge") or 0) or (nomination.get("rna_bulge") or 0):
            reasons.add("bulge")
        mismatch_count = nomination.get("mismatch_count")
        if mismatch_count is not None and mismatch_count > CRISPOR_MAX_MISMATCHES:
            reasons.add("deep_mismatch")
        if guide is not None and nomination.get("pam") and not _pam_matches(
            nomination["pam"], guide["pam_class"]
        ):
            reasons.add("relaxed_pam")

    if metadata is not None:
        states = {state["allele"]: state for state in metadata["allele_states"]}
        if metadata["canonical_interval_basis"] == "alt_only" or "ref" not in states:
            reasons.add("variant_allele_only")
        if str(metadata.get("not_found_in_ref_raw")).strip().casefold() in {"1", "true", "yes", "y"}:
            reasons.add("variant_allele_only")
        ref = states.get("ref")
        alt = states.get("alt")
        if ref is not None and alt is not None:
            if _allele_pam_change(metadata, guide) == "created":
                reasons.add("variant_created_pam")
            ref_mismatches = _mismatches_for_alignment(
                metadata["guide_alignment"], ref["alignment"]
            )
            alt_mismatches = _mismatches_for_alignment(
                metadata["guide_alignment"], alt["alignment"]
            )
            if (
                "N" not in (ref["alignment"] + alt["alignment"]).upper()
                and alt_mismatches < ref_mismatches
            ):
                reasons.add("variant_reduced_mm")

    for nomination in nominations:
        mismatch_count = nomination.get("mismatch_count")
        in_crispor_search_space = (
            mismatch_count == CRISPOR_MAX_MISMATCHES
            and not (nomination.get("dna_bulge") or 0)
            and not (nomination.get("rna_bulge") or 0)
            and "variant_allele_only" not in reasons
            and "relaxed_pam" not in reasons
        )
        if (
            in_crispor_search_space
            and (nomination.get("score_name") or "").casefold() == "cfd"
            and nomination.get("score_value") is not None
            and nomination["score_value"] < CRISPOR_FOUR_MISMATCH_MIN_CFD
        ):
            reasons.add("below_crispor_retention")
    if not reasons:
        reasons.add("unexplained")
    return sorted(reasons)


def _pattern_is_subset(left: str, right: str) -> bool:
    return len(left) == len(right) and all(
        IUPAC[left_base] <= IUPAC[right_base]
        for left_base, right_base in zip(left, right)
    )


def _crisprme_declarations_by_guide(source_provenance: list[dict] | None) -> dict[str, dict]:
    by_guide = {}
    for record in source_provenance or []:
        if record.get("schema_version") != "2":
            continue
        declaration = record["source_declarations"]["crisprme"]
        guide_id = declaration["guide_id"]
        if guide_id in by_guide and by_guide[guide_id] != declaration:
            raise CandidatePanelError(
                f"multiple incompatible CRISPRme declarations for guide {guide_id!r}"
            )
        by_guide[guide_id] = declaration
    return by_guide


def _crispor_declarations_by_guide(source_provenance: list[dict] | None) -> dict[str, dict]:
    by_guide = {}
    for record in source_provenance or []:
        declarations = record.get("source_declarations")
        if not isinstance(declarations, dict) or "crispor" not in declarations:
            continue
        declaration = declarations["crispor"]
        guide_id = declaration["guide_id"]
        if guide_id in by_guide and by_guide[guide_id] != declaration:
            raise CandidatePanelError(
                f"multiple incompatible CRISPOR declarations for guide {guide_id!r}"
            )
        by_guide[guide_id] = declaration
    return by_guide


def _comparison_is_armed(
    site: dict,
    guide: dict | None,
    crisprme_declaration: dict | None,
    crispor_declaration: dict | None,
) -> bool:
    if guide is None or crisprme_declaration is None or crispor_declaration is None:
        return False
    crisprme_retention = crisprme_declaration["retention"]
    crispor_retention = crispor_declaration["retention"]
    pam_covered = all(
        any(
            _pattern_is_subset(pattern, candidate)
            for candidate in crisprme_declaration["pam_patterns"]
        )
        for pattern in crispor_declaration["pam_patterns"]
    )
    return all(
        (
            crisprme_declaration["guide_id"] == guide["guide_id"],
            crispor_declaration["guide_id"] == guide["guide_id"],
            crisprme_declaration["guide_protospacer"]
            == crispor_declaration["guide_protospacer"]
            == guide["spacer_sequence"],
            crisprme_declaration["genome_build"]
            == crispor_declaration["genome_build"]
            == site["genome_build"],
            pam_covered,
            crispor_declaration["max_mismatches"]
            <= crisprme_declaration["max_mismatches"],
            crispor_declaration["max_dna_bulge"]
            <= crisprme_declaration["max_dna_bulge"],
            crispor_declaration["max_rna_bulge"]
            <= crisprme_declaration["max_rna_bulge"],
            crisprme_retention["completeness"] == "complete",
            crisprme_retention["top_n"] is None,
            crisprme_retention["score_cutoff"] is None,
            crisprme_retention["merge_window_bp"] == 0,
            crisprme_retention["all_candidate_coordinates_retained"] is True,
            crispor_retention["completeness"] == "complete",
            crispor_retention["top_n"] is None,
            crispor_retention["score_cutoff"] is None,
            crispor_retention["all_candidate_coordinates_retained"] is True,
        )
    )


def _cannot_assess_editability(allele: str, reason: str) -> dict:
    return {
        "allele": allele,
        "bases": [],
        "summary": "cannot_assess",
        "reason": reason,
    }


def _assess_allele_window(
    allele: str,
    guide_alignment: str,
    target_alignment: str,
    interval: dict | None,
    editor_class: str,
    window_start: int | None = None,
    window_end: int | None = None,
) -> dict:
    edit_base = {"ABE": "A", "CBE": "C"}.get(editor_class.upper())
    if edit_base is None:
        return _cannot_assess_editability(allele, "unknown_editor_class")
    if (
        isinstance(window_start, bool) or not isinstance(window_start, int)
        or isinstance(window_end, bool) or not isinstance(window_end, int)
        or not 1 <= window_start <= window_end <= 20
    ):
        return _cannot_assess_editability(allele, "missing_or_invalid_declared_window")
    if not guide_alignment or not target_alignment:
        return _cannot_assess_editability(allele, "missing_alignment")
    if len(guide_alignment) != len(target_alignment):
        return _cannot_assess_editability(allele, "alignment_length_mismatch")
    if set(guide_alignment.upper()) - set("ACGTN-") or set(target_alignment.upper()) - set(
        "ACGTN-"
    ):
        return _cannot_assess_editability(allele, "unsupported_symbol")
    if "N" in target_alignment.upper() or "N" in guide_alignment.replace("-", "").upper()[:20]:
        return _cannot_assess_editability(allele, "unknown_sequence")
    if sum(base != "-" for base in guide_alignment) < 20:
        return _cannot_assess_editability(allele, "guide_axis_mismatch")

    if interval is None:
        return _cannot_assess_editability(allele, "unresolved_reference_mapping")
    strand = interval["strand"]
    coordinate = (
        interval["start_0based"]
        if strand == "+"
        else interval["end_0based_exclusive"] - 1
    )
    guide_position = 0
    bases = []
    for guide_base, target_base in zip(guide_alignment.upper(), target_alignment.upper()):
        if guide_base != "-":
            guide_position += 1
        target_coordinate = coordinate if target_base != "-" else None
        if (
            guide_base != "-"
            and target_base == edit_base
            and window_start <= guide_position <= window_end
            and target_coordinate is not None
        ):
            bases.append(
                {
                    "position": guide_position,
                    "genomic_position_0based": target_coordinate,
                    "base": target_base,
                    "genomic_forward_base": (
                        target_base if strand == "+" else _reverse_complement(target_base)
                    ),
                    "region": "declared_window",
                }
            )
        if target_base != "-":
            coordinate += 1 if strand == "+" else -1
    expected_final = (
        interval["end_0based_exclusive"]
        if strand == "+"
        else interval["start_0based"] - 1
    )
    if coordinate != expected_final:
        return _cannot_assess_editability(allele, "alignment_length_mismatch")
    summary = "base_in_declared_window" if bases else "no_base_in_declared_window"
    return {"allele": allele, "bases": bases, "summary": summary, "reason": None}


def _alignment_block_for_site(site: dict, guide: dict | None) -> dict | None:
    crisprme = _crisprme_metadata_for_site(site)
    if crisprme is not None:
        return crisprme
    if guide is None:
        return None
    nominations = [
        item
        for item in site["nominations"]
        if item["nomination_kind"] == "computational" and item.get("site_sequence")
    ]
    if not nominations:
        experimental = [n for n in site["nominations"] if n["nomination_kind"] == "experimental"]
        if not experimental or any(
            n.get("reference_validation", {}).get("status") != "verified"
            or n.get("strand") not in {"+", "-"}
            or not n.get("site_sequence") or set(n["site_sequence"]) - set("ACGT")
            or len(n["site_sequence"]) != 23
            or site["end_0based_exclusive"] - site["start_0based"] != 23
            or n.get("dna_bulge") not in {None, 0}
            or n.get("rna_bulge") not in {None, 0}
            or n.get("population_variant_id")
            or len(n.get("pam") or "") != 3
            or n["site_sequence"][-3:] != n["pam"]
            for n in experimental
        ):
            return None
        if len({(n["site_sequence"], n["strand"]) for n in experimental}) != 1:
            return None
        nominations = experimental
    nomination = nominations[0]
    target = nomination["site_sequence"]
    native_query = nomination.get("source_metadata", {}).get("native_query")
    guide_alignment = native_query or (
        guide["spacer_sequence"] + "N" * len(nomination.get("pam") or "")
    )
    if len(guide_alignment) != len(target):
        return None
    return {
        "guide_alignment": guide_alignment,
        "allele_states": [
            {
                "allele": "ref",
                "alignment": target,
                "effective_interval": {
                    "chrom": site["chrom"],
                    "start_0based": site["start_0based"],
                    "end_0based_exclusive": site["end_0based_exclusive"],
                    "strand": nomination["strand"],
                },
                "CFD": (
                    nomination["score_value"]
                    if (nomination.get("score_name") or "").casefold() == "cfd"
                    else None
                ),
                "variant_rsid_raw": nomination.get("population_variant_id"),
                "variant_maf_raw": None,
                "variant_info_raw": {"spacer_pam": None, "genome": None},
                "origin_selected": True,
            }
        ],
    }


def _alignment_blocks_for_site(site: dict, guide: dict | None) -> list[dict]:
    """Keep each reported alignment and DNA version with its source row."""
    blocks = []
    for nomination in site["nominations"]:
        if nomination["source_method"].casefold() != "crisprme":
            continue
        block = nomination["source_metadata"]["crisprme"]
        identifier = nomination.get("source_record_id") or hashlib.sha256(
            json.dumps(block, sort_keys=True).encode("utf-8")
        ).hexdigest()[:20]
        blocks.append({
            **block, "configuration_id": identifier,
            "allele_states": reference_mapped_states(block["allele_states"]),
            "phase_status": "not_provided_in_integrated_contract",
            "genotype_status": "not_imported",
            "configuration_coverage": "reported_highest_CFD_only",
            "interpretation": "reported_sequence_only",
        })
    if blocks:
        return blocks
    block = _alignment_block_for_site(site, guide)
    return [{**block, "configuration_id": "reference"}] if block else []


def _assess_site_windows(site: dict, guide: dict | None) -> dict:
    blocks = _alignment_blocks_for_site(site, guide)
    if not blocks or guide is None:
        result = [_cannot_assess_editability("unknown", "missing_alignment")]
    else:
        result = []
        alignment_axes = {
            (block["guide_alignment"], tuple(
                state["alignment"] for state in block["allele_states"] if state["allele"] == "ref"
            )) for block in blocks
        }
        for block in blocks:
            for state in block["allele_states"]:
                if len(alignment_axes) > 1:
                    assessment = _cannot_assess_editability(
                        state["allele"], "unresolved_alternative_alignment"
                    )
                elif (state["allele"] == "alt"
                        and len(state.get("genomic_forward_changes") or []) > 1):
                    assessment = _cannot_assess_editability(
                        state["allele"], "phase_not_provided_for_multiple_changes"
                    )
                else:
                    assessment = _assess_allele_window(
                        state["allele"], block["guide_alignment"], state["alignment"],
                        state["effective_interval"], guide["editor_class"],
                        guide.get("editing_window_start_1based"),
                        guide.get("editing_window_end_1based_inclusive"),
                    )
                result.append({**assessment, "configuration_id": block["configuration_id"]})
    if any(item["summary"] == "cannot_assess" for item in result):
        site_summary = "cannot_assess"
        label = "At least one allele cannot be assessed from the available sequence and mapping."
    elif any(item["bases"] for item in result):
        site_summary = "base_in_declared_window"
        label = "A substrate base occurs in the declared window; activity is not established."
    else:
        site_summary = "no_base_in_declared_window"
        label = "No substrate base in the declared window; this does not exclude editing elsewhere."
    return {"allele_states": result, "site_summary": site_summary, "label": label}


def _annotate_v3_candidates(
    candidates: list[dict],
    guides: list[dict] | None,
    source_provenance: list[dict] | None,
) -> None:
    guide_by_id = {guide["guide_id"]: guide for guide in (guides or [])}
    crisprme_declaration_by_guide = _crisprme_declarations_by_guide(source_provenance)
    crispor_declaration_by_guide = _crispor_declarations_by_guide(source_provenance)
    for site in candidates:
        guide = guide_by_id.get(site["guide_id"])
        site["source_bucket"] = _computational_source_bucket(site)
        blocks = _alignment_blocks_for_site(site, guide)
        for block in blocks:
            block["pam_change"] = _allele_pam_change(block, guide)
            if "variant_maf_raw" in block:
                declaration = crisprme_declaration_by_guide.get(site["guide_id"], {})
                block["frequency"] = {
                    "raw": block["variant_maf_raw"],
                    "semantics": "uninterpreted_upstream_MAF", "allele": None, "value": None,
                    "population": None, "AC": None, "AN": None, "callability": None,
                    "sources": declaration.get("variant_databases", []),
                }
        site["allele_configurations"] = blocks
        site["allele_states"] = [
            {**state, "configuration_id": block["configuration_id"],
             "guide_alignment": block["guide_alignment"]}
            for block in blocks for state in block["allele_states"]
        ]
        site["why_crispor_missed"] = _why_crispor_missed(site, guide)
        warnings = []
        if site["source_bucket"] == "both":
            reconciliation_status = "matched"
        elif site["source_bucket"] == "crisprme_only":
            reconciliation_status = (
                "single_source_unexplained"
                if "unexplained" in site["why_crispor_missed"]
                else "single_source_explained"
            )
        elif site["source_bucket"] == "crispor_only":
            armed = _comparison_is_armed(
                site,
                guide,
                crisprme_declaration_by_guide.get(site["guide_id"]),
                crispor_declaration_by_guide.get(site["guide_id"]),
            )
            reconciliation_status = (
                "qc_warning_missing_from_crisprme" if armed else "comparison_not_armed"
            )
            if armed:
                warnings.append(
                    "CRISPOR site is absent from a declared complete, compatible CRISPRme import."
                )
        else:
            reconciliation_status = "not_applicable"
        site["reconciliation"] = {
            "status": reconciliation_status,
            "warnings": warnings,
        }
        site["editability"] = _assess_site_windows(site, guide)


def _has_reported_text(value: object) -> bool:
    return isinstance(value, str) and value.strip().casefold() not in {
        "", "unknown", "not_reported", "not reported", "na", "n/a", ".",
        "none", "null", "-", "missing",
    }


def _match_measurements_to_site_alleles(site: dict, guide: dict | None) -> list[dict]:
    """Check every measurement against each reported DNA version at this exact site."""
    scenarios = []
    experimental_model_unsupported = (
        not site["allele_states"] and bool(site.get("nominations"))
        and all(n["nomination_kind"] == "experimental" for n in site["nominations"])
    )
    states = site["allele_states"] or [{"allele": "unknown", "configuration_id": "unknown"}]
    for state in states:
        interval = state.get("effective_interval")
        sequence = state.get("alignment", "").upper()
        mapping_reason = None
        if experimental_model_unsupported:
            mapping_reason = "experimental_candidate_model_unsupported"
        elif state["allele"] == "alt" and len(state.get("genomic_forward_changes") or []) > 1:
            mapping_reason = "candidate_phase_unresolved"
        elif state["allele"] == "alt" and state.get("genomic_forward_changes") is None:
            mapping_reason = "candidate_allele_identity_unresolved"
        elif not sequence or set(sequence) - set("ACGT"):
            mapping_reason = "candidate_sequence_or_alignment_unresolved"
        elif not interval or interval.get("strand") not in {"+", "-"}:
            mapping_reason = "candidate_mapping_unresolved"
        elif any(interval[field] != site[field] for field in ("chrom", "start_0based", "end_0based_exclusive")):
            mapping_reason = "candidate_interval_mismatch"
        elif interval["end_0based_exclusive"] - interval["start_0based"] != len(sequence):
            mapping_reason = "candidate_mapping_unresolved"
        if mapping_reason is None and interval["strand"] == "-":
            sequence = _reverse_complement(sequence)
        scenario = {
            "configuration_id": state["configuration_id"], "allele": state["allele"],
            "sequence_genomic_forward": sequence if mapping_reason is None else None,
            "mapping_reason": mapping_reason, "measurements": [],
        }
        usable = []
        for measurement in site["measurements"]:
            identity_reasons = [mapping_reason] if mapping_reason else []
            for field, expected in (
                ("editor", (guide or {}).get("editor")),
                ("allele", state["allele"]),
                ("allele_sequence", scenario["sequence_genomic_forward"]),
                ("strand", interval.get("strand") if interval else None),
            ):
                observed = measurement.get(field)
                if field == "strand":
                    known = observed in ("+", "-") and expected in ("+", "-")
                else:
                    known = _has_reported_text(observed) and _has_reported_text(expected)
                if not known:
                    identity_reasons.append(f"{field}_unknown")
                elif observed != expected:
                    identity_reasons.append(f"{field}_mismatch")
            if measurement.get("allele_resolution") != "resolved":
                identity_reasons.append("allele_resolution_unknown")
            if not _has_reported_text(measurement.get("allele_evidence")):
                identity_reasons.append("allele_evidence_missing")
            if (measurement.get("configuration_id")
                    and measurement["configuration_id"] != state["configuration_id"]):
                identity_reasons.append("configuration_mismatch")
            quality_reasons = []
            if measurement.get("qc_status") != "passed":
                quality_reasons.append("qc_failed" if measurement.get("qc_status") == "failed"
                                       else "qc_not_reported")
            if measurement.get("status") == "inconclusive":
                quality_reasons.append("inconclusive_result")
            for field in ("assay", "biological_context", "treated_group", "control_group"):
                if not _has_reported_text(measurement.get(field)):
                    quality_reasons.append(f"{field}_missing")
            if measurement.get("measurement_kind") != "editing":
                quality_reasons.append("editing_measurement_not_established")
            if measurement["status"] == "not_detected":
                limit = editing_detection_limit(measurement.get("limit_of_detection"))
                if limit["reason"]:
                    quality_reasons.append(limit["reason"])
            match_status = ("mismatch" if any(r.endswith("_mismatch") for r in identity_reasons)
                            else "unknown" if identity_reasons else "matched")
            is_usable = not identity_reasons and not quality_reasons
            scenario["measurements"].append({
                "source_record_id": measurement["source_record_id"],
                "source_version_or_citation": measurement["source_version_or_citation"],
                "assay": measurement["assay"], "match_status": match_status,
                "usable_for_editing": is_usable,
                "reasons": identity_reasons + quality_reasons,
            })
            if is_usable:
                usable.append(measurement)
        if usable:
            scenario["summary"] = _measurement_summary(usable)
        elif any(m["match_status"] == "matched" for m in scenario["measurements"]):
            scenario["summary"] = "unusable"
        elif any(m["match_status"] == "unknown" for m in scenario["measurements"]):
            scenario["summary"] = "coverage_unresolved"
        else:
            scenario["summary"] = "not_measured"
        labels = {
            "not_detected": "No editing detected within the recorded assay limits",
            "edited": "Editing reported in the matched assay",
            "mixed": "Different results; review each assay and cell context",
            "not_measured": "Untested for this DNA version",
            "unusable": "No usable editing result",
            "coverage_unresolved": "Cannot establish which DNA version was tested",
        }
        scenario["evidence_label"] = labels[scenario["summary"]]
        if experimental_model_unsupported:
            question = (
                "This workflow does not model an allele for this experimental nomination. "
                "Review the candidate representation using the nomination records. "
                "Separately check the tested DNA version and measurement quality; "
                "reference identity does not establish sample alleles."
            )
        elif mapping_reason:
            question = "Resolve the candidate sequence, genomic mapping or phase before matching measurements."
        elif scenario["summary"] == "not_measured":
            question = "Does editing occur for this DNA version with this guide/editor in the relevant cells?"
        elif scenario["summary"] == "coverage_unresolved":
            question = "Confirm the editor, tested sequence and allele-specific coverage from the measurement source."
        elif scenario["summary"] == "unusable":
            question = "Resolve missing controls, quality checks or inconclusive results; obtain usable editing evidence."
        else:
            question = "Review each assay's cell context and detection limit, and whether consequences were measured."
        scenario["follow_up_questions"] = [question]
        scenarios.append(scenario)
    return scenarios


def _coverage_summary(scenarios: list[dict]) -> str:
    statuses = {scenario["summary"] for scenario in scenarios}
    measured = statuses & {"edited", "not_detected", "mixed"}
    if measured and statuses - measured:
        return "partial_coverage"
    if len(statuses) == 1:
        return next(iter(statuses))
    return "mixed" if measured else "coverage_unresolved"




def _add_strand(site: dict, record: dict, evidence_label: str) -> None:
    strand = record.get("strand")
    if not strand:
        return
    if site["strands"] and strand not in site["strands"]:
        raise CandidatePanelError(
            f"conflicting strands at {site['guide_id']} {site['chrom']}:"
            f"{site['start_0based']}-{site['end_0based_exclusive']}; "
            f"existing={','.join(sorted(site['strands']))}, {evidence_label}={strand}"
        )
    site["strands"].add(strand)


def _flag_one_base_coordinate_conflicts(sites: dict[tuple, dict]) -> None:
    """Flag one-base coordinate disagreements without choosing which row is correct."""
    grouped: dict[tuple[str, str, str], list[dict]] = {}
    for site in sites.values():
        grouped.setdefault(
            (site["guide_id"], site["genome_build"], site["chrom"]), []
        ).append(site)
    for group in grouped.values():
        ordered = sorted(group, key=lambda site: (site["start_0based"], site["end_0based_exclusive"]))
        for index, left in enumerate(ordered):
            for right in ordered[index + 1 :]:
                if right["start_0based"] - left["start_0based"] > 1:
                    break
                if (
                    abs(right["start_0based"] - left["start_0based"]) <= 1
                    and abs(right["end_0based_exclusive"] - left["end_0based_exclusive"]) <= 1
                ):
                    left_text = f"{left['chrom']}:{left['start_0based']}-{left['end_0based_exclusive']}"
                    right_text = f"{right['chrom']}:{right['start_0based']}-{right['end_0based_exclusive']}"
                    warning = (
                        "probable one-base coordinate-convention disagreement with "
                        f"{right_text}; not merged"
                    )
                    reverse_warning = (
                        "probable one-base coordinate-convention disagreement with "
                        f"{left_text}; not merged"
                    )
                    left["identity_warnings"].append(warning)
                    right["identity_warnings"].append(reverse_warning)


def merge_evidence(
    nominations: list[dict],
    measurements: list[dict],
    decisions: dict,
    intended_loci: list[dict] | None = None,
    locus_requirements: list[dict] | None = None,
    guides: list[dict] | None = None,
    source_provenance: list[dict] | None = None,
) -> list[dict]:
    sites: dict[tuple[str, str, str, int, int], dict] = {}
    requirement_by_guide = {
        requirement["guide_id"]: requirement for requirement in (locus_requirements or [])
    }
    for locus in intended_loci or []:
        key = _site_key(locus)
        site = sites.setdefault(key, _new_site(locus))
        site["site_role"] = "intended"
        site["intended_locus"] = locus
        _add_strand(site, locus, "intended-locus manifest")
        requirement = requirement_by_guide.get(locus["guide_id"])
        if requirement:
            site["locus_requirement_ids"].append(requirement["requirement_id"])
            site["required_confirmation"] = requirement["required_confirmation"]
    for nomination in nominations:
        key = _site_key(nomination)
        site = sites.setdefault(key, _new_site(nomination))
        role = nomination["site_role"]
        if site["site_role"] is not None and site["site_role"] != role:
            raise CandidatePanelError(
                f"conflicting site_role declarations for {nomination['guide_id']} "
                f"{nomination['chrom']}:{nomination['start_0based']}-"
                f"{nomination['end_0based_exclusive']}"
            )
        site["site_role"] = role
        _add_strand(site, nomination, "nomination")
        site["nominations"].append(nomination)

    for measurement in measurements:
        key = _site_key(measurement)
        site = sites.setdefault(key, _new_site(measurement))
        _add_strand(site, measurement, "measurement")
        site["measurements"].append(measurement)

    _flag_one_base_coordinate_conflicts(sites)

    unknown_decisions = sorted(set(decisions) - set(sites))
    if unknown_decisions:
        key = unknown_decisions[0]
        raise CandidatePanelError(
            "panel decision does not match a nomination or measurement: "
            f"{key[0]} {key[2]}:{key[3]}-{key[4]}"
        )
    known_site_ids = {site["site_id"] for site in sites.values()}
    for key, decision in decisions.items():
        related_site_id = decision.get("related_site_id")
        if related_site_id is None:
            continue
        if related_site_id not in known_site_ids:
            raise CandidatePanelError(
                f"panel decision related_site_id {related_site_id!r} does not match a site"
            )
        if related_site_id == sites[key]["site_id"]:
            raise CandidatePanelError(
                "panel decision redundant_with cannot reference the same site"
            )

    result = []
    for key in sorted(sites):
        site = sites[key]
        if site["site_role"] is None:
            site["site_role"] = "candidate"
        site["strands"] = sorted(site["strands"])
        site["identity_warnings"] = sorted(set(site["identity_warnings"]))
        site["locus_requirement_ids"] = sorted(set(site["locus_requirement_ids"]))
        site["nominations"].sort(
            key=lambda item: (
                item["nomination_kind"], item["source_method"], item.get("source_record_id") or ""
            )
        )
        site["measurements"].sort(
            key=lambda item: (item["assay"], item["status"], item.get("source_record_id") or "")
        )
        site["panel_decision"] = decisions.get(
            key,
            {
                "status": "undecided",
                "reason": "not_reviewed",
                "note": None,
                "related_site_id": None,
                "source_file": None,
            },
        )
        result.append(site)
    _annotate_v3_candidates(result, guides, source_provenance)
    guide_by_id = {guide["guide_id"]: guide for guide in (guides or [])}
    for site in result:
        site["reported_measurement_statuses"] = sorted({m["status"] for m in site["measurements"]})
        site["measurement_coverage"] = _match_measurements_to_site_alleles(site, guide_by_id.get(site["guide_id"]))
        site["measurement_summary"] = _coverage_summary(site["measurement_coverage"])
        site["follow_up_questions"] = list(dict.fromkeys([
            *(question for scenario in site["measurement_coverage"] for question in scenario["follow_up_questions"]),
            *site["identity_warnings"],
            *(item["assay"] for item in site["required_confirmation"]),
        ]))
        site["evidence_state"] = _evidence_state(site)
        site["review_group"], site["review_reasons"] = _review_group(site)
    return result


def _input_record(kind: str, path: Path) -> dict:
    return {"kind": kind, "file": path.name, "sha256": _sha256(path)}


def compare_ledgers(previous: dict, current: dict) -> dict:
    """List exact-site differences between two v3 records without adding removed candidates."""
    for label, ledger in (("previous", previous), ("current", current)):
        if (
            not isinstance(ledger, dict)
            or ledger.get("record_type") != "candidate_panel_ledger"
            or ledger.get("schema_version") != "3"
            or not isinstance(ledger.get("candidates"), list)
        ):
            raise CandidatePanelError(f"{label} diff input must be a candidate ledger v3")
    previous_guides = {item["guide_id"] for item in previous.get("guides", [])}
    current_guides = {item["guide_id"] for item in current.get("guides", [])}
    if previous_guides != current_guides:
        raise CandidatePanelError("ledger diff requires identical guide sets")
    previous_builds = {item["genome_build"] for item in previous["candidates"]}
    current_builds = {item["genome_build"] for item in current["candidates"]}
    if previous_builds != current_builds:
        raise CandidatePanelError("ledger diff requires identical genome-build sets")

    def by_id(ledger: dict, label: str) -> dict[str, dict]:
        mapped = {item.get("site_id"): item for item in ledger["candidates"]}
        if None in mapped or len(mapped) != len(ledger["candidates"]):
            raise CandidatePanelError(f"{label} ledger contains missing or duplicate site IDs")
        return mapped

    old = by_id(previous, "previous")
    new = by_id(current, "current")
    old_ids = set(old)
    new_ids = set(new)
    rebucketed = [
        {
            "site_id": site_id,
            "old_source_bucket": old[site_id].get("source_bucket"),
            "new_source_bucket": new[site_id].get("source_bucket"),
        }
        for site_id in sorted(old_ids & new_ids)
        if old[site_id].get("source_bucket") != new[site_id].get("source_bucket")
    ]
    return {
        "added_site_ids": sorted(new_ids - old_ids),
        "removed_site_ids": sorted(old_ids - new_ids),
        "rebucketed": rebucketed,
    }


def build_ledger(
    profile_path: Path,
    guide_path: Path,
    intended_locus_path: Path,
    nomination_paths: list[Path],
    measurement_path: Path | None = None,
    decision_path: Path | None = None,
    source_provenance_paths: list[Path] | None = None,
    reference_fasta_path: Path | None = None,
    previous_ledger_path: Path | None = None,
    panel_notes: str | None = None,
) -> dict:
    """Read, check and combine candidate lists, measurements and panel decisions."""
    if not nomination_paths:
        raise CandidatePanelError("at least one --nominations file is required")
    profile = SearchProfile.load(profile_path)
    guides = load_guides(guide_path)
    intended_loci = load_intended_loci(intended_locus_path)
    locus_requirements = _validate_intended_loci(guides, intended_loci, profile)
    nominations = load_nominations(nomination_paths, profile)
    if profile.reference_id == "SYNTHETIC" and reference_fasta_path is None:
        reference_record = None
        for nomination in nominations:
            nomination["reference_validation"] = {
                "status": "not_applicable_synthetic_fixture"
            }
    else:
        if reference_fasta_path is None:
            raise CandidatePanelError(
                "--reference-fasta is required for non-synthetic candidate panels"
            )
        reference_record = validate_reference_sequences(nominations, reference_fasta_path)
    source_provenance = load_source_provenance(
        source_provenance_paths,
        nomination_paths,
        nominations,
        profile,
        profile_path,
        guide_path,
        intended_locus_path,
        reference_record,
    )
    measurements = load_measurements(measurement_path)
    decisions = load_decisions(decision_path)
    if panel_notes is not None:
        if (
            not isinstance(panel_notes, str)
            or len(panel_notes) > 500
            or "\t" in panel_notes
            or "\n" in panel_notes
            or "\r" in panel_notes
        ):
            raise CandidatePanelError(
                "panel_notes must be at most 500 characters with no tabs or newlines"
            )
        panel_notes = panel_notes.strip() or None
    all_paths = [profile_path, guide_path, intended_locus_path, *nomination_paths]
    if measurement_path is not None:
        all_paths.append(measurement_path)
    if decision_path is not None:
        all_paths.append(decision_path)
    all_paths.extend(source_provenance_paths or [])
    if previous_ledger_path is not None:
        all_paths.append(previous_ledger_path)
    if reference_fasta_path is not None:
        all_paths.extend([reference_fasta_path, Path(str(reference_fasta_path) + ".fai")])
    duplicate_basenames = sorted(
        name for name, count in Counter(path.name for path in all_paths).items() if count > 1
    )
    if duplicate_basenames:
        raise CandidatePanelError(
            "input basenames must be unique so checksums remain unambiguous: "
            + ", ".join(duplicate_basenames)
        )
    wrong_reference = sorted(
        {
            item["genome_build"]
            for item in [*nominations, *measurements]
            if item["genome_build"] != profile.reference_id
        }
        | {key[1] for key in decisions if key[1] != profile.reference_id}
    )
    if wrong_reference:
        raise CandidatePanelError(
            "all evidence must use profile reference_id "
            f"{profile.reference_id!r}; found: {', '.join(wrong_reference)}"
        )
    known_guides = {item["guide_id"] for item in guides}
    referenced_guides = {
        item["guide_id"] for item in [*nominations, *measurements]
    } | {key[0] for key in decisions} | {item["guide_id"] for item in intended_loci}
    unknown_guides = sorted(referenced_guides - known_guides)
    if unknown_guides:
        raise CandidatePanelError(
            "evidence references guide_id values absent from the guide manifest: "
            + ", ".join(unknown_guides)
        )
    _validate_nomination_content(nominations, guides, intended_loci)
    _validate_measurement_intended_strands(measurements, intended_loci)
    candidates = merge_evidence(
        nominations,
        measurements,
        decisions,
        intended_loci,
        locus_requirements,
        guides,
        source_provenance,
    )

    inputs = [
        _input_record("search_profile", profile_path),
        _input_record("guide_manifest", guide_path),
        _input_record("intended_locus_manifest", intended_locus_path),
    ]
    inputs.extend(_input_record("nominations", path) for path in nomination_paths)
    if measurement_path is not None:
        inputs.append(_input_record("measurements", measurement_path))
    if decision_path is not None:
        inputs.append(_input_record("panel_decisions", decision_path))
    inputs.extend(
        _input_record("nomination_source_provenance", path)
        for path in (source_provenance_paths or [])
    )
    if reference_fasta_path is not None:
        inputs.extend(
            [
                _input_record("reference_fasta", reference_fasta_path),
                _input_record("reference_fasta_index", Path(str(reference_fasta_path) + ".fai")),
            ]
        )
    if previous_ledger_path is not None:
        inputs.append(_input_record("previous_candidate_ledger", previous_ledger_path))

    import_summaries = [
        {
            "normalized_output": record["normalized_output"]["file"],
            **record["import_summary"],
        }
        for record in source_provenance
        if record["schema_version"] == "2"
    ]

    ledger = {
        "record_type": "candidate_panel_ledger",
        "schema_version": LEDGER_SCHEMA_VERSION,
        "search_profile": profile.to_dict(),
        "guides": guides,
        "intended_loci": intended_loci,
        "locus_requirements": locus_requirements,
        "nomination_source_provenance": source_provenance,
        "import_summaries": import_summaries,
        "export_metadata": {"panel_notes": panel_notes},
        "diff_summary": None,
        "reference_validation": (
            reference_record
            if reference_record is not None
            else {"status": "not_applicable_synthetic_fixture"}
        ),
        "inputs": inputs,
        "summary": {
            "candidate_sites": len(candidates),
            "nomination_rows": len(nominations),
            "measurement_rows": len(measurements),
            "by_measurement_coverage": dict(sorted(Counter(c["measurement_summary"] for c in candidates).items())),
            "identity_warning_sites": sum(bool(c["identity_warnings"]) for c in candidates),
            "by_source_bucket": dict(
                sorted(
                    Counter(
                        c["source_bucket"]
                        for c in candidates
                        if c["source_bucket"] is not None
                    ).items()
                )
            ),
            "by_reconciliation_status": dict(
                sorted(
                    Counter(c["reconciliation"]["status"] for c in candidates).items()
                )
            ),
            "by_evidence_state": dict(sorted(Counter(c["evidence_state"] for c in candidates).items())),
            "by_review_group": dict(
                sorted(Counter(c["review_group"] for c in candidates).items())
            ),
            "by_panel_decision": dict(
                sorted(Counter(c["panel_decision"]["status"] for c in candidates).items())
            ),
        },
        "candidates": candidates,
        "interpretation": {
            "nomination": "A reason to inspect or test a site; not proof that editing occurs.",
            "measurement": (
                "An assay result in its stated biological context; not_detected is bounded by the "
                "assay and its detection limit."
            ),
            "panel_decision": (
                "A researcher-supplied decision for targeted follow-up."
            ),
        },
        "limitations": [
            "This workflow reconciles imported candidates; it does not perform genome-wide nomination.",
            "Search completeness is limited by the external method and the predeclared profile.",
            "An unmeasured nomination is a candidate, not an off-target finding.",
            "Measurement coverage uses declared source metadata; genotype and assay quality are not independently audited.",
            "A usable negative is bounded by its tested allele, guide/editor, assay, cells, controls and detection limit.",
            "Consequence measurements are retained separately and do not establish editing measurement coverage.",
            "An excluded site is absent from this panel, not demonstrated safe.",
            "One-base coordinate warnings require source reconciliation and are never auto-merged.",
        ],
        "reproducibility": {
            "generated_at_utc": run_stamp(),
            "base_edit_safety_version": __version__,
            "python_version": platform.python_version(),
            "repository": repo_state(),
        },
    }
    if previous_ledger_path is not None:
        try:
            previous = json.loads(previous_ledger_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise CandidatePanelError(
                f"cannot read previous candidate ledger {previous_ledger_path}"
            ) from exc
        ledger["diff_summary"] = compare_ledgers(previous, ledger)
    return ledger


PANEL_COLUMNS = [
    "site_id",
    "guide_id",
    "genome_build",
    "chrom",
    "start_0based",
    "end_0based_exclusive",
    "strands",
    "site_role",
    "source_bucket",
    "reconciliation_status",
    "reconciliation_warnings",
    "why_crispor_missed",
    "evidence_state",
    "review_group",
    "review_reasons",
    "identity_warnings",
    "locus_requirements",
    "required_confirmation",
    "nomination_kinds",
    "source_methods",
    "source_versions_or_citations",
    "source_scores",
    "mismatch_counts",
    "editable_base_in_window",
    "population_variant_ids",
    "population_frequency_raw",
    "allele_states_json",
    "editability_summary",
    "editability_bases_json",
    "ch_driver_genes",
    "annotation_warnings",
    "measurement_summary",
    "reported_measurement_statuses",
    "measurement_coverage_json",
    "follow_up_questions",
    "assays",
    "replicate_ids",
    "read_depths",
    "measurement_qc_statuses",
    "detection_limits",
    "panel_decision",
    "decision_reason",
    "decision_note",
    "related_site_id",
    "interpretation_warning",
]


def _joined(values: Iterable[object]) -> str:
    return ";".join(sorted({str(value) for value in values if value is not None and value != ""}))


def panel_rows(candidates: list[dict]) -> list[dict]:
    rows = []
    for site in candidates:
        nominations = site["nominations"]
        measurements = site["measurements"]
        rows.append(
            {
                "site_id": site["site_id"],
                "guide_id": site["guide_id"],
                "genome_build": site["genome_build"],
                "chrom": site["chrom"],
                "start_0based": site["start_0based"],
                "end_0based_exclusive": site["end_0based_exclusive"],
                "strands": _joined(site["strands"]),
                "site_role": site["site_role"],
                "source_bucket": site["source_bucket"] or "",
                "reconciliation_status": site["reconciliation"]["status"],
                "reconciliation_warnings": _joined(site["reconciliation"]["warnings"]),
                "why_crispor_missed": _joined(site["why_crispor_missed"]),
                "evidence_state": site["evidence_state"],
                "review_group": site["review_group"],
                "review_reasons": _joined(site["review_reasons"]),
                "identity_warnings": _joined(site["identity_warnings"]),
                "locus_requirements": _joined(site["locus_requirement_ids"]),
                "required_confirmation": _joined(
                    item["assay"] for item in site["required_confirmation"]
                ),
                "nomination_kinds": _joined(item["nomination_kind"] for item in nominations),
                "source_methods": _joined(item["source_method"] for item in nominations),
                "source_versions_or_citations": _joined(
                    item["source_version_or_citation"] for item in [*nominations, *measurements]
                ),
                "source_scores": _joined(
                    f"{item['source_method']}:{item['score_name']}={item['score_value']}"
                    for item in nominations
                    if item.get("score_name") is not None
                    and item.get("score_value") is not None
                ),
                "mismatch_counts": _joined(item["mismatch_count"] for item in nominations),
                "editable_base_in_window": _joined(
                    str(item["editable_base_in_window"]).lower()
                    for item in nominations
                    if item["editable_base_in_window"] is not None
                ),
                "population_variant_ids": _joined(
                    item["population_variant_id"] for item in nominations
                ),
                "population_frequency_raw": _joined(
                    state.get("variant_maf_raw") for state in site["allele_states"]
                ),
                "allele_states_json": json.dumps(
                    site["allele_states"], sort_keys=True, separators=(",", ":")
                ),
                "editability_summary": site["editability"]["site_summary"],
                "editability_bases_json": json.dumps(
                    site["editability"]["allele_states"],
                    sort_keys=True,
                    separators=(",", ":"),
                ),
                "ch_driver_genes": _joined(item["ch_driver_gene"] for item in nominations),
                "annotation_warnings": _joined(
                    f"{field}: imported, not independently verified"
                    for item in nominations
                    for field, status in item.get("annotation_validation", {}).items()
                    if status == "imported_not_verified"
                ),
                "measurement_summary": site["measurement_summary"],
                "reported_measurement_statuses": _joined(site["reported_measurement_statuses"]),
                "measurement_coverage_json": json.dumps(site["measurement_coverage"], sort_keys=True),
                "follow_up_questions": _joined(site["follow_up_questions"]),
                "assays": _joined(item["assay"] for item in measurements),
                "replicate_ids": _joined(item.get("replicate_id") for item in measurements),
                "read_depths": _joined(item.get("read_depth") for item in measurements),
                "measurement_qc_statuses": _joined(
                    item.get("qc_status") for item in measurements
                ),
                "detection_limits": _joined(
                    f"{item['limit_of_detection']['value']} {item['limit_of_detection']['unit']}"
                    for item in measurements
                    if item["limit_of_detection"] is not None
                ),
                "panel_decision": site["panel_decision"]["status"],
                "decision_reason": site["panel_decision"]["reason"],
                "decision_note": site["panel_decision"].get("note") or "",
                "related_site_id": site["panel_decision"].get("related_site_id") or "",
                "interpretation_warning": (
                    "Panel exclusion is not evidence that a site is safe."
                ),
            }
        )
    return rows


def write_outputs(ledger: dict, ledger_path: Path, panel_path: Path) -> None:
    ledger_path.parent.mkdir(parents=True, exist_ok=True)
    panel_path.parent.mkdir(parents=True, exist_ok=True)
    ledger_path.write_text(json.dumps(ledger, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    with panel_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=PANEL_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(panel_rows(ledger["candidates"]))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", type=Path, required=True, help="predeclared search-profile JSON")
    parser.add_argument("--guides", type=Path, required=True, help="guide and editor manifest TSV")
    parser.add_argument(
        "--intended-loci",
        type=Path,
        required=True,
        help="canonical intended-locus manifest TSV",
    )
    parser.add_argument(
        "--nominations",
        type=Path,
        action="append",
        required=True,
        help="normalized nomination TSV; repeat for multiple sources",
    )
    parser.add_argument("--measurements", type=Path, help="optional normalized measurement TSV")
    parser.add_argument("--decisions", type=Path, help="optional explicit panel-decision TSV")
    parser.add_argument(
        "--reference-fasta",
        type=Path,
        help="indexed reference FASTA; required unless reference_id is SYNTHETIC",
    )
    parser.add_argument(
        "--source-provenance",
        type=Path,
        action="append",
        required=True,
        help="source-run JSON for a normalized nomination file; repeat to cover every input",
    )
    parser.add_argument(
        "--previous-ledger",
        type=Path,
        help="optional candidate-ledger v3 to compare by exact site identity",
    )
    parser.add_argument(
        "--panel-notes",
        help="optional export-level sensitivity or depth context (maximum 500 characters)",
    )
    parser.add_argument("--output", type=Path, required=True, help="candidate-ledger JSON output")
    parser.add_argument("--panel-output", type=Path, required=True, help="panel-review TSV output")
    args = parser.parse_args(argv)
    try:
        ledger = build_ledger(
            args.profile,
            args.guides,
            args.intended_loci,
            args.nominations,
            args.measurements,
            args.decisions,
            args.source_provenance,
            args.reference_fasta,
            args.previous_ledger,
            args.panel_notes,
        )
        write_outputs(ledger, args.output, args.panel_output)
    except CandidatePanelError as exc:
        parser.error(str(exc))
    print(
        f"wrote {ledger['summary']['candidate_sites']} candidate sites to {args.output} "
        f"and {args.panel_output}"
    )
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
