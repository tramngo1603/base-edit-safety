"""Compare CRISPOR candidate sites with published off-target panels for the same guide.

The workflow's candidate list comes from one CRISPOR run. A published panel for the same guide
provides a second search method, so differences identify sites that one method did not return.

These panels contain sites selected for testing by another method. A site list is not a measurement,
so importing one must not change `base_editor_context` from `not_assessed_here`. Detection results,
when available, belong to a separate results table. Agreement between two candidate searches does
not establish that editing occurred.

The comparison asks whether CRISPOR returns the sites in the published panel. A missing site marks a
coverage gap in the CRISPOR search. The comparison cannot determine whether an extra site reflects
a method difference or a false positive.
"""
from __future__ import annotations

from dataclasses import dataclass

from . import assay_protocols


@dataclass(frozen=True)
class PanelSite:
    name: str
    protospacer: str        # 20 nt as published, on the strand the study reported
    pam: str
    chrom: str
    pos_1based: int
    mismatch_count: int


@dataclass(frozen=True)
class PublishedPanel:
    """Sites a published study selected for off-target analysis of one guide."""

    source_pmid: str
    source_table: str
    guide_spacer: str               # guide sequence stored by this workflow
    guide_names: dict               # what each source calls the same molecule
    selection_basis: str
    is_measurement: bool            # False for candidate-site panels; see the module docstring
    sites: tuple

    @property
    def mismatch_counts(self) -> dict:
        counts: dict = {}
        for s in self.sites:
            counts[s.mismatch_count] = counts.get(s.mismatch_count, 0) + 1
        return counts


# Ravi 2022, Supplementary file 3, "The targets analyzed for DNA off-target".
#
# The table carries no guide column. All eleven sites were assigned here by aligning each published
# protospacer against every guide in Ravi Supplementary file 1: all eleven are 3-mismatch off-targets
# of gRNA 11 in reverse-complement orientation. No other stored guide has the same three-mismatch
# relationship. This workflow names gRNA 11 as 61forw,
# one of the four designs above the editing-activity threshold.
#
# The table has no editing-rate column and no confirmed/refuted column, so it records what was looked
# at and not what was found.
RAVI_2022_GRNA11 = PublishedPanel(
    source_pmid="35147495",
    source_table="Supplementary file 3",
    guide_spacer="GGCAAGGCTGGCCAACCCAT",
    guide_names={"ravi": "gRNA 11", "ours": "61forw"},
    selection_basis="sites selected by the study for DNA off-target analysis; all 3 mismatches",
    is_measurement=False,
    sites=(
        PanelSite("OT1",  "TTGGGTTGCCCAGACTTGCC", "AGG", "chr1",  25_039_865,  3),
        PanelSite("OT2",  "GTGGGGTGGCCAGCCTTTCC", "TGG", "chr1",  235_427_633, 3),
        PanelSite("OT3",  "GTGGGGTGGCCAGCCTTCCC", "TGG", "chr2",  64_808_323,  3),
        PanelSite("OT4",  "ATGGGGTGTGCAGCCTTGCC", "TGG", "chr12", 113_177_162, 3),
        PanelSite("OT5",  "AGGGGTTGGCCCGCCTTCCC", "AGG", "chr16", 674_938,     3),
        PanelSite("OT6",  "TTGGCTTGGCCAGACTTGCC", "CGG", "chr9",  110_575_343, 3),
        PanelSite("OT7",  "GTGGGGTGGCCAGCCTTCCC", "AGG", "chr9",  118_816_545, 3),
        PanelSite("OT8",  "ATGGGTTGTCCAGCGTTGCT", "TGG", "chrX",  150_946_920, 3),
        PanelSite("OT9",  "ATGGGTTGGCCAGGCTGGTC", "TGG", "chr6",  7_187_976,   3),
        PanelSite("OT10", "GTGGGGTGGCCAGCCTTCCC", "TGG", "chr10", 119_547_882, 3),
        PanelSite("OT11", "ATGGGTTGACCAACCTGGCC", "AGG", "chr13", 86_143_790,  3),
    ),
)

PANELS = (RAVI_2022_GRNA11,)
_BY_SPACER = {p.guide_spacer: p for p in PANELS}


def panel_for(protospacer: str) -> PublishedPanel | None:
    """Return the published panel recorded for this guide, if present."""
    return _BY_SPACER.get(protospacer.upper())


def compare_to_panel(our_sites: list[dict], panel: PublishedPanel,
                     position_tolerance: int = 25) -> dict:
    """Return which published panel sites the current search found or missed.

    `our_sites` are per-site dictionaries from the guide-design adapter and require `chrom` and
    `start`. Matching uses coordinates within a fixed tolerance because the two sources report site
    ends and PAM positions differently. Sequence matching would also require the reference sequence
    at each location.

    A miss identifies a published candidate site absent from the single-tool search and adds that
    location to follow-up review.
    """
    ours = []
    for s in our_sites:
        chrom, start = s.get("chrom"), s.get("start")
        if chrom is None or start is None:
            continue
        try:
            ours.append((chrom if str(chrom).startswith("chr") else f"chr{chrom}", int(start)))
        except (TypeError, ValueError):
            continue

    found, missed = [], []
    for site in panel.sites:
        hit = next((o for o in ours
                    if o[0] == site.chrom and abs(o[1] - site.pos_1based) <= position_tolerance),
                   None)
        (found if hit else missed).append(site.name)

    return {
        "panel": f"PMID {panel.source_pmid} {panel.source_table}",
        "guide": panel.guide_names,
        "panel_size": len(panel.sites),
        "found_by_our_search": found,
        "missed_by_our_search": missed,
        "our_site_count": len(ours),
        "is_measurement": panel.is_measurement,
        "note": ("Both lists contain candidates from different search methods. Agreement does not "
                 "show that editing occurred, and this comparison does not change the base_editor_context status, "
                 "which remains not_assessed_here until these loci are sequenced."),
        # Attach published amplicons for missed sites so follow-up can test those locations.
        "amplicons_for_missed_sites": {
            name: assay_protocols.offtarget_primers_for(name) for name in missed
        },
    }
