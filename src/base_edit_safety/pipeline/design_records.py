"""Write the primary review output: one design record per candidate guide.

Each record leads with the evidence and the assays it needs, not with a score.
"""
from __future__ import annotations

from ..design_record import build_design_records
from ..io import outputs
from ..scoring.rubric import RubricWeights


def run_design_records(record_dicts: list[dict] | None = None,
                       weights: RubricWeights | None = None,
                       empirical: dict | None = None, write: bool = True) -> dict:
    """Build records from the supplied rows, or read the current enumeration file."""
    weights = weights or RubricWeights()
    if record_dicts is None:
        record_dicts = outputs.read_enumeration()["records"]

    record_set = build_design_records(record_dicts, weights, empirical)
    if write:
        outputs.write_design_records(
            record_set, filename="design_records_HBG1_HBG2.json",
            meta={"step": "design records", "count": record_set["count"]},
        )
    return record_set


if __name__ == "__main__":
    rs = run_design_records()
    print(f"design records: {rs['count']}\n")
    for r in rs["records"]:
        c = r["candidate"]
        edits = ",".join(f"-{e['promoter_offset']}" for e in c["target"]["intended_edits"])
        print(f"  {c['guide_id']:<22} {c['gene']} {edits:<10} "
              f"tier={r['evidence_tier']:<18} registry={r['registry_match']['relationship']:<12} "
              f"paralog={r['hbg_paralog_flags']['flagged']!s:<5} "
              f"assays={len(r['required_validation'])}")
