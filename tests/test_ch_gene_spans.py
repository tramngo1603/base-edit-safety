"""The twelve CH-driver gene spans are frozen, so an upstream refresh cannot change them quietly.

These boundaries decide whether an off-target counts as landing in a driver gene. They are not stable
across Ensembl releases: four of the twelve differ against release 115, nine against release 110, by
up to 586 bp. Without a pinned release and a fixture, a refresh would silently move flags.

If one of these fails, that is not a bug to paper over. Work out which release the new value comes
from, update the pin and the fixture together, and check what it does to the flags.
"""
from __future__ import annotations

from base_edit_safety import annotation

# Verified against Ensembl release 116 / GRCh38.p14, by stable ID.
EXPECTED = {
    "ASXL1":  ("chr20", 32_358_266, 32_439_319, "ENSG00000171456"),
    "CHEK2":  ("chr22", 28_687_738, 28_742_422, "ENSG00000183765"),
    "DNMT3A": ("chr2",  25_227_855, 25_342_590, "ENSG00000119772"),
    "GNAS":   ("chr20", 58_839_669, 58_911_192, "ENSG00000087460"),
    "IDH2":   ("chr15", 90_083_045, 90_102_511, "ENSG00000182054"),
    "MDM4":   ("chr1", 204_516_369, 204_558_120, "ENSG00000198625"),
    "PPM1D":  ("chr17", 60_599_735, 60_666_293, "ENSG00000170836"),
    "SF3B1":  ("chr2", 197_388_515, 197_435_079, "ENSG00000115524"),
    "SRSF2":  ("chr17", 76_734_115, 76_737_919, "ENSG00000161547"),
    "TET2":   ("chr4", 105_145_875, 105_279_816, "ENSG00000168769"),
    "TP53":   ("chr17",  7_661_779, 7_687_546, "ENSG00000141510"),
    "U2AF1":  ("chr21", 43_092_956, 43_107_605, "ENSG00000160201"),
}


def test_the_release_is_pinned():
    assert annotation.CH_GENE_SPANS_ENSEMBL_RELEASE == "116"
    assert annotation.CH_GENE_SPANS_ASSEMBLY == "GRCh38.p14"
    assert "116" in annotation.CH_DRIVER_SET_PROVENANCE, "outputs must say which release they used"


def test_every_span_matches_the_pinned_release():
    for gene, (chrom, start, end, _ensg) in EXPECTED.items():
        assert annotation.CH_GENE_SPANS[gene] == (chrom, start, end), gene


def test_the_panel_is_exactly_these_twelve_genes():
    assert set(annotation.CH_GENE_SPANS) == set(EXPECTED)
    assert annotation.CH_DRIVER_GENES == set(EXPECTED)


def test_spans_are_well_formed():
    for gene, (_chrom, start, end) in annotation.CH_GENE_SPANS.items():
        assert start < end, gene


def test_the_one_observed_off_target_still_falls_in_dnmt3a():
    """The real run's only CH hit. If a span moves, this is where it shows up first."""
    assert annotation.gene_at("chr2", 25_264_110) == "DNMT3A"


def test_a_site_just_outside_a_span_is_not_flagged():
    chrom, start, end, _ = EXPECTED["TP53"]
    assert annotation.gene_at(chrom, start - 1) is None
    assert annotation.gene_at(chrom, end + 1) is None
    assert annotation.gene_at(chrom, start) == "TP53"
    assert annotation.gene_at(chrom, end) == "TP53"
