"""Write versioned pipeline records under `data/outputs`.

Each output carries `ARTIFACT_SCHEMA_VERSION`. `read_enumeration()` rejects a different version so
older fields cannot be interpreted under the current field definitions. Increase the version when a field or
its meaning changes in a way that makes an older file incompatible.
"""
from __future__ import annotations

import dataclasses
import json
import os

from ..pipeline.datatypes import EditRecord

# 4: multi-edit enumeration records carry an allele-keyed joint_on_target_rate and its model,
#    method, confidence, and source fields. BE-DICT-backed records also carry uniform-v1 score
#    score-scaling details and the unscaled total; run metadata stores totals and exclusions.
#    Design-record version 2 reads the joint value.
# 3: enumeration records carry `offtarget_sites_low_mismatch`, the per-site coordinates for
#    <=3-mismatch off-targets. Without these a run cannot be compared against a published panel,
#    which permits comparison of the eleven Ravi sites with the twelve CRISPOR sites for 61forw.
# 2: records carry `confidence` (model_predicted | not_evaluated) and a rate of None where the
#    predictor was never asked, the header key is `scope`, and meta carries a reproducibility block.
# 1: files written before the version-2 fields were introduced.
ARTIFACT_SCHEMA_VERSION = "4"

# Header keys that were removed at version 2. Their presence dates a file on its own.
_RETIRED_HEADER_KEYS = ("disclaimer", "claim_disclaim")


class StaleArtifact(RuntimeError):
    """Raised when an output file on disk was written by a different build of this code."""

# Include scope in every output so copied files retain their intended use and limits.
SCOPE = {
    "what_this_is": (
        "A shortlist of base-editing guide and editor designs for prioritizing laboratory tests. "
        "Ordering uses the provisional triage rubric described in the output."
    ),
    "what_this_is_not": [
        "Not a safety clearance: no design here is called safe.",
        "Not an efficacy ranking: efficiency is an eligibility threshold, not the ranking target.",
        "Not a per-person prediction: the clonal-hematopoiesis layer flags gene membership and "
        "uses it as a triage flag, and produces no patient-level estimate.",
        "Not measured: on-target rates are model predictions and off-targets are candidate sites.",
    ],
    "thresholds": (
        "The efficiency and weighting thresholds are starting values taken from allelic-editing "
        "literature. Re-estimate them with the final guide, editor, process, and a measured "
        "potency assay."
    ),
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
    """List each way a document differs from the current required output fields.

    Return all differences so the error identifies every incompatible field in one pass.
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
        found.append(f"meta.step is {step!r}; 'spine' was renamed to 'single_guide'")
    if not (meta.get("reproducibility") or {}):
        found.append("meta.reproducibility is empty, so the run cannot be traced to its inputs")

    if "records" not in doc:
        found.append("no 'records' array")
    return found


def read_enumeration() -> dict:
    """Load an enumeration file only when its version matches the current fields.

    Later steps share this file, so rejecting older incompatible schemas prevents calculations from using
    fields whose meaning has changed.
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


def write_credibility(payload_body: dict, filename: str, meta: dict | None = None) -> str:
    """Write the divergence and uncertainty report."""
    return _write_json(filename, payload_body, meta)


def write_ranking(ranked: list[dict], filename: str, meta: dict | None = None) -> str:
    """Write the ranked shortlist."""
    return _write_json(filename, {"guide_count": len(ranked), "ranking": ranked}, meta)


def write_design_records(record_set: dict, filename: str, meta: dict | None = None) -> str:
    """Write the design records, the main output."""
    return _write_json(filename, record_set, meta)
