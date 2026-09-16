"""Flag the need to distinguish the two HBG copies and check larger changes.

HBG1 and HBG2 have similar promoters about 4.9 kb apart. A reference-sequence search alone
cannot establish which copy edits or whether a deletion or rearrangement occurs between them.
This module records those concerns and the assays needed to check them.
"""
from __future__ import annotations

from . import anchors


def paralog_flags(record: dict) -> dict:
    """Flag the HBG copy checks and assays needed for one edit record."""
    prov = record.get("provenance", {})
    value = prov.get("co_targets_paralog")
    if value is None or str(value).strip().casefold() in {"", "false", "0", "no", "n"}:
        co_targets = False
    elif str(value).strip().casefold() in {"true", "1", "yes", "y"}:
        co_targets = True
    else:
        raise ValueError("co_targets_paralog must be an explicit boolean or missing")
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
