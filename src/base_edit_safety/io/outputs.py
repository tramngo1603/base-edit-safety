"""Read and write versioned pipeline artifacts under ``data/outputs``.

The schema version prevents an older output file from being read with newer field meanings. Bump
it whenever a shape change could make an existing file ambiguous.
"""
from __future__ import annotations

import dataclasses
import json
import os

from ..pipeline.datatypes import EditRecord

# 4: multi-edit enumeration records carry an allele-keyed joint_on_target_rate and its model,
#    method, confidence, and caveat provenance. BE-DICT-backed records also carry uniform-v1 score
#    normalization provenance and raw mass; run metadata records the policy and any exclusions.
# 3: enumeration records carry `offtarget_sites_low_mismatch`, the per-site coordinates for
#    <=3-mismatch off-targets. Without these a run cannot be compared against a published panel,
#    which is how eleven Ravi sites for 61forw sat uncomparable against our own twelve.
# 2: records carry `confidence` (model_predicted | not_evaluated) and a rate of None where the
#    predictor was never asked, the header key is `scope`, and meta carries a reproducibility block.
# 1: anything written before that, including every file in data/outputs as of 2026-07-26.
ARTIFACT_SCHEMA_VERSION = "4"

# Header keys that were removed at version 2. Their presence dates a file on its own.
_RETIRED_HEADER_KEYS = ("disclaimer", "claim_disclaim")


class StaleArtifact(RuntimeError):
    """Raised when an output file on disk was written by a different build of this code."""

# Goes at the top of every output file, so a file that gets passed around still says what it is.
SCOPE = {
    "what_this_is": (
        "Base-editing guide and editor records with predicted edits, candidate off-targets, "
        "evidence links, and follow-up assays."
    ),
    "what_this_is_not": [
        "A replacement for experimental confirmation or review of the underlying evidence.",
    ],
    "confirmation_assays": "amplicon NGS for on-target, CHANGE-seq-BE or GUIDE-seq for off-target.",
}

OUTPUT_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))),
    "data",
    "outputs",
)


def output_path(filename: str) -> str:
    return os.path.join(OUTPUT_DIR, filename)


def enumeration_path() -> str:
    """The enumeration file that later steps read back."""
    return output_path("enumeration_HBG1_HBG2.json")


def records_to_dicts(records: list[EditRecord]) -> list[dict]:
    return [dataclasses.asdict(r) for r in records]


def _write_json(filename: str, body: dict, meta: dict | None) -> str:
    """Write {artifact_schema_version, scope, meta, **body} as indented JSON under data/outputs."""
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    path = output_path(filename)
    payload = {"artifact_schema_version": ARTIFACT_SCHEMA_VERSION,
               "scope": SCOPE, "meta": meta or {}, **body}
    with open(path, "w") as handle:
        json.dump(payload, handle, indent=2)
    return path


def schema_complaints(doc: dict) -> list[str]:
    """Every way this document disagrees with what the current writer produces.

    Returns all of them rather than the first, because a file written by an old build usually
    disagrees in several ways at once and the full list is what dates it.
    """
    found = []
    version = doc.get("artifact_schema_version")
    if version is None:
        found.append(
            f"no artifact_schema_version (current writer emits {ARTIFACT_SCHEMA_VERSION!r}); "
            f"written before the field existed"
        )
    elif version != ARTIFACT_SCHEMA_VERSION:
        found.append(f"artifact_schema_version is {version!r}, current is {ARTIFACT_SCHEMA_VERSION!r}")

    for key in _RETIRED_HEADER_KEYS:
        if key in doc:
            found.append(f"carries retired header key {key!r}, replaced by 'scope' at version 2")
    if "scope" not in doc:
        found.append("no 'scope' header")

    meta = doc.get("meta") or {}
    step = str(meta.get("step", ""))
    if "spine" in step:
        found.append(f"meta.step is {step!r}; 'spine' is a retired step label")
    if not (meta.get("reproducibility") or {}):
        found.append("meta.reproducibility is empty, so the run cannot be traced to its inputs")

    if "records" not in doc:
        found.append("no 'records' array")
    return found


def read_enumeration() -> dict:
    """Load the enumeration file, refusing anything a different build wrote.

    Every later step reads this file. Letting a stale one through means every downstream number is
    computed from records whose meaning has since changed, which is exactly what happened here.
    """
    path = enumeration_path()
    if not os.path.isfile(path):
        raise StaleArtifact(
            f"no enumeration file at {path}. Run "
            f"`python -m base_edit_safety.pipeline.enumerate_guides` first; it needs CRISPOR and "
            f"BE-DICT configured (see docs/setup.md)."
        )
    with open(path) as handle:
        doc = json.load(handle)
    complaints = schema_complaints(doc)
    if complaints:
        raise StaleArtifact(
            f"{path} was not written by this build of the code:\n"
            + "\n".join(f"  - {c}" for c in complaints)
            + "\n\nRegenerate it with `python -m base_edit_safety.pipeline.enumerate_guides`. "
              "Do not validate results against it: the records predate the current meaning of "
              "predicted_rate and confidence."
        )
    return doc


def write_records(records: list[EditRecord], filename: str, meta: dict | None = None) -> str:
    """Write edit records. Returns the path."""
    return _write_json(filename, {"record_count": len(records), "records": records_to_dicts(records)}, meta)


def write_design_records(record_set: dict, filename: str, meta: dict | None = None) -> str:
    """Write the design records, the main output."""
    return _write_json(filename, record_set, meta)
