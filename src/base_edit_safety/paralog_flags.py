"""Flags for editing inside the HBG1/HBG2 duplication, and the assays those flags call for.

HBG1 and HBG2 came from a recent duplication. Their proximal promoters are nearly identical and the
genes sit about 4.9 kb apart on chr11. Editing there can produce things a generic off-target tool
cannot see at all: both copies edited at once, and large deletions or rearrangements between the two
genes (the ~4.9 kb inter-HBG deletion is the known one). This module flags candidates in that region
and lists the assays that would settle it.
"""
from __future__ import annotations

from . import anchors


def paralog_flags(record: dict) -> dict:
    """Flag one edit record for paralog trouble, and list the assays it needs."""
    prov = record.get("provenance", {})
    co_targets = bool(prov.get("co_targets_paralog"))
    in_promoter_window = record.get("gene") in ("HBG1", "HBG2")

    flags = []
    if in_promoter_window:
        flags.append("edits the duplicated HBG1/HBG2 proximal promoter (near-identical sequence)")
    if co_targets:
        flags.append(f"guide also hits the paralog ({prov.get('co_target_gene')})")

    flagged = bool(flags)
    assays = []
    if flagged:
        assays = [
            {"validates": "hbg_paralog_rearrangement",
             "assay": "long-read or amplicon sequencing spanning HBG1-HBG2, to catch the ~4.9 kb "
                      "deletion or a rearrangement", "status": "required"},
            {"validates": "hbg_paralog_cotargeting",
             "assay": "amplicon NGS that can tell HBG1 and HBG2 apart, measuring editing at each",
             "status": "required"},
        ]
    return {
        "flagged": flagged,
        "paralog_spacing_bp": anchors.PARALOG_SPACING_BP,
        "flags": flags,
        "required_confirmation": assays,
    }
