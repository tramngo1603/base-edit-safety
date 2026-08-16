"""Build a self-contained JSON view of predicted bystander edits.

The module reads enumeration output, anchor definitions, and bystander notes, then writes a display
file. It does not modify the enumeration or scoring results.

A bystander is an unintended edit near the target and inside the editing window. An off-target is a
different genomic location. Values come from the enumeration output; a missing value remains `null`.
"""
from __future__ import annotations

import hashlib
import json

from . import outputs
from .. import anchors
from ..bystanders import allele_context, bystander_function, rule_for

_ENUMERATION_FILE = outputs.enumeration_path()
_DISPLAY_FILE = outputs.output_path("bystander_display.json")

EDITOR_NOTE = ("a bystander is an unintended edit near the target, inside the editing window. "
               "It is not an off-target, which would be elsewhere in the genome.")

# Positions with exact-allele HSPC evidence. See `docs/coordinates.md`.
PRIMARY_ANCHORS = {175, 113, 123, 124}

# How the positions are grouped for display; -123/-124 belong together.
TARGET_GROUPS = [
    ("-198", [198]),
    ("-175", [175]),
    ("-113", [113]),
    ("-123/-124", [123, 124]),
    ("-117", [117]),
    ("-114", [114]),
]

# Notes per position, drawn from the anchor table and the rubric.
TARGET_NOTES = {
    "-123/-124": ("A synthetic double mutant. KLF1 binding was shown for the -123/-124 pair "
                  "(Ravi 2022), not for either position alone, so both are intended edits here."),
    "-175": ("The main bystander here breaks the TAL1 E-box. Editor choice requires comparing "
             "functional effects: ABE8e can be preferred despite more bystanders (Mayuranathan 2023). "
             "Only ABE8e was run here, so compare editors in the lab."),
    "-113": ("The neighbouring -112 bystander breaks the GATA1 motif created by the intended edit; -116 "
             "helps. What a bystander does matters more than how many there are."),
    "-114": ("The neighbouring C edits at -114 and -115 are both part of breaking the distal "
             "BCL11A motif, which is the point."),
}


def _rate(value, model_name):
    """Wrap a predicted rate so it stays labelled as a prediction. Null if there is no value."""
    return {
        "value": value,
        "model_predicted": True,
        "measured": False,
        "source": model_name,
        "note": "BE-DICT prediction from cell-line and library data, not measured in HSCs",
    }


def _anchors_by_offset() -> dict:
    return {p.promoter_offset: p for p in anchors.IN_SCOPE}


def _functional_source(promoter_offset: int, ref: str, alt: str) -> str | None:
    """Return the publication identifiers stored on the matching bystander rule."""
    rule = rule_for(promoter_offset, ref, alt)
    if rule is None:
        return None
    return f"bystanders.RULES: PMID {', PMID '.join(rule.source_pmid)} ({rule.confidence})"


def build_display(records: list[dict]) -> dict:
    by_offset = _anchors_by_offset()

    # Index the HBG1 records by guide. HBG2 mirrors them, since rates are per protospacer.
    by_guide: dict[str, list[dict]] = {}
    for r in records:
        if r["gene"] != "HBG1":
            continue
        by_guide.setdefault(r["provenance"]["guide_id"], []).append(r)

    def representative_guide(offsets: list[int]) -> str | None:
        best, best_rate = None, None
        for gid, recs in by_guide.items():
            ot = [x for x in recs if x["on_target_or_bystander"] == "on_target"
                  and x["promoter_offset"] in offsets]
            if not ot:
                continue
            rate = max((x["predicted_rate"] or 0.0) for x in ot)
            if best is None or rate > best_rate:
                best, best_rate = gid, rate
        return best

    targets = []
    for label, offsets in TARGET_GROUPS:
        ap = by_offset[offsets[0]]
        gid = representative_guide(offsets)
        recs = by_guide.get(gid, []) if gid else []
        model_name = recs[0]["model_name"] if recs else None

        intended = []
        for off in offsets:
            a = by_offset[off]
            ot = [x for x in recs if x["on_target_or_bystander"] == "on_target"
                  and x["promoter_offset"] == off]
            rate_val = ot[0]["predicted_rate"] if ot else None
            intended.append({
                "position": f"-{off}",
                "genomic": {"chrom": anchors.CHROM,
                            "hbg1_pos": a.hbg1_pos, "hbg2_pos": a.hbg2_pos,
                            "ref_alt": f"{a.genomic_plus_ref}>{a.genomic_plus_alt}"},
                "sense": f"{a.promoter_sense_ref}>{a.promoter_sense_alt} (promoter/sense strand)",
                "predicted_rate": _rate(rate_val, model_name),
            })

        on_allele = allele_context(recs)
        bystanders = []
        for x in recs:
            if x["on_target_or_bystander"] != "bystander":
                continue
            off = x["promoter_offset"]
            sref, salt = x["promoter_sense_ref"], x["promoter_sense_alt"]
            cls, mech = bystander_function(x["gene"], off, sref, salt, on_allele)
            bystanders.append({
                "position": f"-{off}",
                "genomic": {"chrom": x["chrom"], "hbg1_pos": x["pos_1based"],
                            "ref_alt": f"{x['genomic_plus_ref']}>{x['genomic_plus_alt']}"},
                "sense": f"{x['promoter_sense_ref']}>{x['promoter_sense_alt']} (promoter/sense strand)",
                "predicted_rate": _rate(x["predicted_rate"], x["model_name"]),
                "functional_class": cls,
                "mechanism": mech,
                "functional_source": _functional_source(off, sref, salt),
            })
        bystanders.sort(key=lambda b: -(b["predicted_rate"]["value"] or 0.0))

        targets.append({
            "label": label,
            "gene_loci": ["HBG1", "HBG2"],
            "editor_class": ap.editor_class,
            "intended_mechanism": ap.mechanism,
            "is_primary_anchor": any(o in PRIMARY_ANCHORS for o in offsets),
            "variant_db_status": {"HBG1": by_offset[offsets[0]].hbg1_status,
                                  "HBG2": by_offset[offsets[0]].hbg2_status},
            "representative_guide": gid,
            "intended_edits": intended,
            "bystanders": bystanders,
            "note": TARGET_NOTES.get(label),
        })

    return {
        "meta": {
            "generated_from": {
                "enumeration_file": "data/outputs/enumeration_HBG1_HBG2.json",
                "run_ref_sha256": _enumeration_sha256(),
                "pipeline_step": "a view of the enumeration results; the pipeline is unchanged",
            },
            "scope": outputs.SCOPE,
            "editor_note": EDITOR_NOTE,
        },
        "targets": targets,
    }


def _enumeration_sha256() -> str:
    with open(_ENUMERATION_FILE, "rb") as handle:
        return hashlib.sha256(handle.read()).hexdigest()[:16]


def write_display() -> str:
    records = outputs.read_enumeration()["records"]
    display = build_display(records)
    with open(_DISPLAY_FILE, "w") as handle:
        json.dump(display, handle, indent=2)
    return _DISPLAY_FILE


if __name__ == "__main__":
    print(write_display())
