"""Flag an edit whose coordinate falls inside a listed CH driver gene.

CH means clonal hematopoiesis. This module performs a location check against the 12-gene boostDM-CH
panel. It reports the matched gene, coordinate, sequence mismatch count, and whether the site is in
an exon or intron. It does not estimate a person's risk, disease rate, or time to an event.
`tests/test_no_risk_field.py` enforces that limit.

The boostDM-CH gene panel is CC-BY-NC. Gene boundaries come from the fixed Ensembl release recorded
below. Scores for individual DNA changes are not included.
"""
from __future__ import annotations

from typing import Optional

from .pipeline.datatypes import EditRecord

# boostDM-CH driver-gene panel (12 genes; scores from the source project are CC-BY-NC).
CH_DRIVER_GENES = frozenset(
    {"ASXL1", "CHEK2", "DNMT3A", "GNAS", "IDH2", "MDM4", "PPM1D", "SF3B1", "SRSF2", "TET2",
     "TP53", "U2AF1"}
)

# GRCh38 gene spans, 1-based and inclusive, looked up in Ensembl by symbol.
#
# The release is fixed because gene boundaries can change between releases. Compared with release
# 115, four of the twelve spans differ, including SRSF2 by 524 bp; compared with release 110, nine
# differ by up to 586 bp. `tests/test_ch_gene_spans.py` pins these values so a release update requires
# an explicit code change and matching fixed test data.
#
# Each span covers the full gene interval, including all transcripts and untranslated regions. This
# allows one coordinate check to cover every listed transcript.
CH_GENE_SPANS_ENSEMBL_RELEASE = "116"
CH_GENE_SPANS_ASSEMBLY = "GRCh38.p14"
CH_GENE_SPANS = {
    "ASXL1": ("chr20", 32358266, 32439319),
    "CHEK2": ("chr22", 28687738, 28742422),
    "DNMT3A": ("chr2", 25227855, 25342590),
    "GNAS": ("chr20", 58839669, 58911192),
    "IDH2": ("chr15", 90083045, 90102511),
    "MDM4": ("chr1", 204516369, 204558120),
    "PPM1D": ("chr17", 60599735, 60666293),
    "SF3B1": ("chr2", 197388515, 197435079),
    "SRSF2": ("chr17", 76734115, 76737919),
    "TET2": ("chr4", 105145875, 105279816),
    "TP53": ("chr17", 7661779, 7687546),
    "U2AF1": ("chr21", 43092956, 43107605),
}
CH_DRIVER_SET_PROVENANCE = (
    f"boostDM-CH 12-gene panel (CC-BY-NC); gene spans from Ensembl release "
    f"{CH_GENE_SPANS_ENSEMBL_RELEASE} on {CH_GENE_SPANS_ASSEMBLY}"
)


def _norm_chrom(chrom: str) -> str:
    return chrom if chrom.startswith("chr") else f"chr{chrom}"


def gene_at(chrom: str, pos_1based: int) -> Optional[str]:
    """The CH-driver gene containing this coordinate, or None."""
    chrom = _norm_chrom(chrom)
    for gene, (gchrom, start, end) in CH_GENE_SPANS.items():
        if chrom == gchrom and start <= pos_1based <= end:
            return gene
    return None


def annotate_record(record: EditRecord) -> EditRecord:
    """Flag whether this edit's own position sits in a CH-driver gene."""
    gene = gene_at(record.chrom, record.pos_1based)
    record.ch_driver_flag = gene is not None
    record.ch_driver_gene = gene
    record.provenance = {
        **record.provenance,
        "ch_driver_set": CH_DRIVER_SET_PROVENANCE,
        "ch_annotation": "gene membership decided by coordinate",
    }
    return record


def annotate_records(records: list[EditRecord]) -> list[EditRecord]:
    """Flag every record."""
    return [annotate_record(r) for r in records]


def ch_driver_offtargets(offtarget_sites: list[dict]) -> list[dict]:
    """Keep the off-target sites whose coordinate falls in a CH-driver gene.

    Each input site needs `chrom` and `start`. The optional `locusDesc` field supplies the feature
    label. Return the matched gene with the search tool's mismatch and score fields.
    """
    hits: list[dict] = []
    for site in offtarget_sites:
        chrom = site.get("chrom")
        start = site.get("start")
        if chrom is None or start is None:
            continue
        try:
            pos = int(start)
        except (TypeError, ValueError):
            continue
        # Membership is decided by coordinate alone. A tool-provided location label such as
        # A label such as "intergenic:TP53-WRAP53" names flanking genes. Use the coordinate, not
        # that label, to set the overlap flag.
        gene = gene_at(chrom, pos)
        if gene is not None:
            locus_desc = site.get("locusDesc", "") or ""
            # CRISPOR prefixes the locus with the feature: exon, intron, or intergenic (UTRs come
            # through as exon).
            feature = locus_desc.split(":", 1)[0] if ":" in locus_desc else "unknown"
            hits.append({
                "ch_driver_gene": gene,
                "chrom": _norm_chrom(chrom),
                "start": pos,
                "strand": site.get("strand"),
                "mismatch_count": site.get("mismatchCount"),
                "cfd_offtarget_score": site.get("cfdOfftargetScore"),
                "locus_desc": locus_desc,
                "feature": feature,                # exon | intron | intergenic, from the coordinate
                "boostdm_ch_score": None,          # only filled in for coding sites, once scores exist
            })
    return hits


# Coding changes have different meanings across the listed genes. The table separates genes linked
# to loss of function from genes linked to changes at named residues. These are review categories,
# not numerical scores.
CH_DRIVER_MECHANISM = {
    "ASXL1": "lof_tumor_suppressor",
    "CHEK2": "lof_tumor_suppressor",
    "TET2": "lof_tumor_suppressor",
    "DNMT3A": "loss_or_hypomorphic_with_R882_hotspot",
    "TP53": "tumor_suppressor_mixed_lof_dominant_negative",
    "PPM1D": "truncating_activating",       # truncating the C-terminus stabilizes PPM1D
    "SF3B1": "spliceosome_hotspot",
    "SRSF2": "spliceosome_hotspot",
    "U2AF1": "spliceosome_hotspot",
    "IDH2": "gof_hotspot",                  # R140 / R172
    "GNAS": "gof_hotspot",                  # R201
    "MDM4": "copy_number_or_overexpression_driver",   # a single-base change does not represent the
                                                      # copy-number mechanism
}


def ch_consequence_class(gene: str, feature: str | None) -> dict:
    """Read the coordinate feature together with the gene's mechanism.

    Distinguish genes where coding disruption is relevant from genes associated with specific
    residues. The return value contains a category, not a numerical or clinical estimate.
    """
    mechanism = CH_DRIVER_MECHANISM.get(gene, "unknown")
    feat = (feature or "unknown").lower()
    coding = feat in ("exon", "utr")
    coding_relevant = {"lof_tumor_suppressor", "loss_or_hypomorphic_with_R882_hotspot",
                       "tumor_suppressor_mixed_lof_dominant_negative"}
    hotspot = {"gof_hotspot", "spliceosome_hotspot", "truncating_activating"}
    if mechanism in coding_relevant:
        relevance = "coding_disruption_relevant" if coding else "non_coding_lower_concern"
    elif mechanism in hotspot:
        relevance = "hotspot_specific" if coding else "non_coding_lower_concern"
    elif mechanism == "copy_number_or_overexpression_driver":
        relevance = "snv_consequence_poorly_represented"
    else:
        relevance = "unknown"
    return {
        "gene": gene,
        "feature": feat,
        "driver_mechanism": mechanism,
        "consequence_relevance": relevance,
    }
