"""Store published primer sets for required laboratory checks.

The primers distinguish HBG1 from HBG2, detect the deletion between them, or amplify a listed
off-target site. They come from Ravi 2022 Supplementary file 2. The deletion qPCR primers retain
that table's attribution to Chang Li et al., Blood 2018.

A primer set is a protocol, not a result. Its presence does not show that a candidate was tested and
does not change the candidate's evidence level. Adapter overhangs are stored separately from the
sequence that binds the target.
"""
from __future__ import annotations

from dataclasses import dataclass

# Illumina adapter overhangs used on the NGS and off-target primers in this table.
_ADAPTER_F = "tacacgacgctcttccgatct"
_ADAPTER_R = "agacgtgtgctcttccgatct"
_ADAPTER_SHORT = "gctcttccgatct"


@dataclass(frozen=True)
class Primer:
    name: str
    target_specific: str          # sequence that binds the target
    adapter: str = ""             # overhang as published, empty when there is none

    @property
    def as_published(self) -> str:
        return f"{self.adapter}{self.target_specific}"


@dataclass(frozen=True)
class Protocol:
    """A published primer set for one checklist item."""

    validates: str                # matches the `validates` key on a required_validation item
    method: str
    source_pmid: str
    source_table: str
    attribution: str              # original primer authors, when different from the cited table
    primers: tuple
    amplicon_size: str = ""
    note: str = ""


PARALOG_SEPARATION = Protocol(
    validates="hbg_paralog_cotargeting",
    method="gene-specific PCR: a forward primer unique to each paralog against a shared reverse",
    source_pmid="35147495", source_table="Supplementary file 2",
    attribution="Ravi et al. 2022",
    primers=(
        Primer("HBG1 F", "ccacagtacctgccaaagaa"),
        Primer("HBG2 F", "ccatagtatctggtaaagagca"),
        Primer("HBG1/2 R", "ggcgtctggactaggag"),
    ),
    amplicon_size="940 bp",
    note="The two forward primers distinguish the HBG genes; the reverse primer is shared.",
)

PARALOG_SEPARATION_NGS = Protocol(
    validates="hbg_paralog_cotargeting",
    method="amplicon NGS across the promoter, with a paralog-distinguishing forward primer",
    source_pmid="35147495", source_table="Supplementary file 2",
    attribution="Ravi et al. 2022",
    primers=(
        Primer("NGS 4 F", "cggctgacaaaagaagtcct", _ADAPTER_SHORT),
        Primer("NGS 2 R", "gtgaaatgacccatggcgtc", _ADAPTER_SHORT),
    ),
    note="Supplementary file 1 pairs NGS 4 F with NGS 2 R for HBG1/HBG2 separation. The primers "
         "appear in different rows of the source table.",
)

INTER_HBG_DELETION = Protocol(
    validates="hbg_paralog_rearrangement",
    method="qPCR across the deletion junction, with a control amplicon outside it",
    source_pmid="35147495", source_table="Supplementary file 2",
    attribution="Chang Li et al., Blood 2018, as cited by Ravi 2022",
    primers=(
        Primer("Del_qPCR_F", "aggggctcaacgaagaaaaagtgt"),
        Primer("Del_qPCR_R", "cacttcattgtagttaccgtggaaaga"),
        Primer("Del_ctrl_qPCR_F", "aaatgaatcagcagaggctcac"),
        Primer("Del_ctrl_qPCR_R", "atgcactaacatccaactatacaaaa"),
    ),
    note="Detects the ~4.9 kb inter-HBG deletion the paralog layer warns about. qPCR quantifies "
         "deletion alleles; it does not characterise a rearrangement, so the long-read arm of that "
         "checklist item still stands.",
)

# Per-site amplicons for the eleven off-target loci in Ravi Supplementary file 3. The primer names
# use the names "G11 OT1" through "G11 OT11". Those names and the sequence alignment both assign
# the panel to gRNA 11.
_OFFTARGET_PRIMERS = {
    "OT1":  ("ggcagaggggacacatcagt", "tgcgctatgatgttgggtatgtcc", "192bp"),
    "OT2":  ("acttaatgctcacagattggttcaatcaggc", "agccaatcaggttcaagatggcag", "190bp"),
    "OT3":  ("gtctactggaaagcccatttgcatagtaag", "gtgagagagagaggacttctgag", "200bp"),
    "OT4":  ("tcctcagttccaagccttggg", "gtcatccccaatccacaaactcac", "194bp"),
    "OT5":  ("cacacgtgtcttatctgtcacctc", "ctgaccactcttgcaactccatg", "193bp"),
    "OT6":  ("atcagggaaaaacctgccatgtgc", "gcaattcctctgttcaagacaggatgtag", "191bp"),
    "OT7":  ("gggcaaggggaaacttcatttgca", "gagaggcgacttctgagagg", "192bp"),
    "OT8":  ("ttacaggcatgaaccaccacatcc", "gtcagagaactgagcctaatggag", "194bp"),
    "OT9":  ("cacttgagatcaggatttggagaccag", "gcgcaatcttgggtcactgca", "192bp"),
    "OT10": ("ggcccacagattggttagatcag", "ctacttgggactaggcatgttcag", "194bp"),
    "OT11": ("gggggctactctaataaaaagactttaggaac",
             "gttccttatcagcttatgtagattttgggctgag", "200bp"),
}

OFFTARGET_AMPLICONS = tuple(
    Protocol(
        validates="spacer_homology_cas_style",
        method=f"amplicon NGS at off-target site {site}",
        source_pmid="35147495", source_table="Supplementary file 2",
        attribution="Ravi et al. 2022",
        primers=(Primer(f"G11 {site} F", fwd, _ADAPTER_F),
                 Primer(f"G11 {site} R", rev, _ADAPTER_R)),
        amplicon_size=size,
        note="Amplifies one of the eleven sites analysed for gRNA 11, named 61forw here.",
    )
    for site, (fwd, rev, size) in _OFFTARGET_PRIMERS.items()
)

ALL = (PARALOG_SEPARATION, PARALOG_SEPARATION_NGS, INTER_HBG_DELETION) + OFFTARGET_AMPLICONS


def _as_dict(p: Protocol) -> dict:
    return {
        "method": p.method,
        "source": f"PMID {p.source_pmid} {p.source_table}",
        "attribution": p.attribution,
        "primers": [{"name": x.name, "sequence": x.target_specific,
                     **({"adapter_overhang": x.adapter} if x.adapter else {})}
                    for x in p.primers],
        **({"amplicon_size": p.amplicon_size} if p.amplicon_size else {}),
        **({"note": p.note} if p.note else {}),
    }


def protocols_for(validates: str) -> list[dict]:
    """Return published primer sets for a checklist item, or an empty list if none are recorded."""
    return [_as_dict(p) for p in ALL if p.validates == validates]


def offtarget_primers_for(site_name: str) -> dict | None:
    """The amplicon primers for one named site in the gRNA 11 off-target panel."""
    for p in OFFTARGET_AMPLICONS:
        if p.method.endswith(site_name):
            return _as_dict(p)
    return None
