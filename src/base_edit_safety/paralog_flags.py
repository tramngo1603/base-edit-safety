"""Flag HBG1/HBG2 duplication concerns and list assays that can resolve them.

HBG1 and HBG2 came from a duplication. Their promoter sequences are highly similar, and their
transcription start sites are 4,924 bp apart on chromosome 11. A sequence search does not establish
whether both copies were edited or whether a deletion or rearrangement formed between them. This
module flags both questions and lists laboratory tests that can answer them.
"""
from __future__ import annotations

from . import anchors, assay_protocols


def paralog_flags(record: dict) -> dict:
    """Return HBG gene-copy flags and the assays required to resolve them."""
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
        # Include a published primer set when available because the assay label alone does not
        # specify how to distinguish the near-identical HBG1 and HBG2 sequences.
        for item in assays:
            published = assay_protocols.protocols_for(item["validates"])
            if published:
                item["published_protocols"] = published
    return {
        "flagged": flagged,
        "paralog_spacing_bp": anchors.PARALOG_SPACING_BP,
        "flags": flags,
        "required_confirmation": assays,
    }
