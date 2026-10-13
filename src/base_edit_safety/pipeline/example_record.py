"""Build a small no-data example of the HBG review record.

The example uses the checked-in GRCh38 anchor table, invokes no predictor, and contains no
paper-derived rows or measurements. It exists so a new reviewer can inspect the record shape
before installing CRISPOR or BE-DICT.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
import json
from pathlib import Path

from .. import anchors
from ..design_record import build_design_records
from ..io.outputs import records_to_dicts
from .datatypes import NOT_EVALUATED


EXAMPLE_OFFSET = 113
EXAMPLE_EDITOR = "ABE8e"


def build_example() -> dict:
    """Return two review records: the HBG1 and HBG2 copies of one intended edit."""
    selected = [
        record
        for record in anchors.load_in_scope_records()
        if record.promoter_offset == EXAMPLE_OFFSET and record.editor_class == "ABE"
    ]
    if {record.gene for record in selected} != {"HBG1", "HBG2"}:
        raise RuntimeError("the example requires both HBG1 and HBG2 anchor records")

    example_records = []
    for record in selected:
        other_gene = "HBG2" if record.gene == "HBG1" else "HBG1"
        provenance = {
            **record.provenance,
            "source": "checked-in GRCh38 anchor table; no external data",
            "guide_id": f"example_{record.gene}_minus_{EXAMPLE_OFFSET}_{EXAMPLE_EDITOR}",
            "editor": EXAMPLE_EDITOR,
            "protospacer": None,
            "pam": None,
            "pam_class": "NGG",
            "strand": "-",
            "predictor_runs": None,
            "rate_semantics": None,
            "co_targets_paralog": True,
            "co_target_gene": other_gene,
            "offtarget": {
                "source": "not_run_in_synthetic_example",
                "mm_histogram": {},
                "cfd_weighted_burden": None,
                "paralog_count": 1,
                "offtarget_count_excl_paralog": 0,
            },
            "offtarget_sites_low_mismatch": [],
            "reason": "synthetic record example; no outcome predictor or off-target search was run",
        }
        example_records.append(
            replace(
                record,
                predicted_rate=None,
                model_name=None,
                confidence=NOT_EVALUATED,
                consequence=NOT_EVALUATED,
                provenance=provenance,
            )
        )

    result = build_design_records(records_to_dicts(example_records))
    for design in result["records"]:
        registry_match = design["registry_match"]
        for edit_match in registry_match["by_edit"]:
            for matched_entry in edit_match["matched_entries"]:
                matched_entry.pop("citations", None)
        registry_match["note"] = (
            "Aggregate relationship labels are shown for record-shape inspection. Citation and "
            "source-audit details are intentionally omitted from this synthetic example."
        )
    result["example"] = {
        "synthetic": True,
        "external_files_used": False,
        "predictor_run": False,
        "checked_in_sources": ["GRCh38 anchor table", "evidence registry"],
        "registry_detail": "aggregate_only",
        "purpose": "record-shape inspection only",
    }
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        help="write the JSON example to this path; otherwise print it to standard output",
    )
    args = parser.parse_args(argv)
    result = build_example()
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
        print(f"wrote {result['count']} synthetic design records to {args.output}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
