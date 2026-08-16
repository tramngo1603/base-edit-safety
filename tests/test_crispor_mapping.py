"""The CRISPOR adapter works out where a guide sits from the reference window.

No external tool and no network: the reference window is supplied by hand, and both a
forward-strand and a reverse-strand guide are checked.
"""
from base_edit_safety.adapters.crispor import CrisporGuideDesigner
from base_edit_safety.config import CrisporConfig
from base_edit_safety.io.reference import ReferenceWindow

# HBG1 proximal-promoter window, GRCh38 plus strand, chr11:5,249,937-5,250,067.
HBG1_WIN = ReferenceWindow(
    genome_build="GRCh38", chrom="chr11", start_1based=5_249_937, end_1based=5_250_067,
    plus_sequence=(
        "GACTATTGGTCAAGTTTGCCTTGTCAAGGCTATTGGTCAAGGCAAGGCTGGCCAACCCATGG"
        "GTGGAGTTTAGCCAGGGACCGTTTCAGACAGATATTTGCATTGAGATAGTGTGGGGAAGGGGCCCCCAA"
    ),
)

DESIGNER = CrisporGuideDesigner(config=CrisporConfig("x", "y", "hg38", "z"))


def _row(guide_id, target_seq):
    return {
        "guideId": guide_id, "targetSeq": target_seq,
        "mitSpecScore": "45", "cfdSpecScore": "90", "offtargetCount": "100",
        "targetGenomeGeneLocus": "intergenic:HBG1-HBG2",
    }


def test_reverse_strand_guide_maps_to_high_coord_5p_end():
    # The -113 ABE guide: protospacer CTTGACCAATAGCCTTGACA, PAM AGG, minus strand.
    cand = DESIGNER._row_to_candidate(_row("19rev", "CTTGACCAATAGCCTTGACAAGG"), HBG1_WIN, "HBG1")
    assert cand.strand == "-"
    assert cand.proto_5p_genomic == 5_249_977
    # protospacer position 8 -> genomic 5,249,970 = HBG1 -113
    guide = cand.for_editor("ABE", "ABE8e")
    assert guide.genomic_coord(8) == 5_249_970
    assert guide.promoter_offset(8) == 113


def test_forward_strand_guide_maps_to_low_coord_5p_end():
    plus = HBG1_WIN.plus_sequence
    target = plus[20:43]  # an arbitrary 23-mer present on the plus strand
    cand = DESIGNER._row_to_candidate(_row("21forw", target), HBG1_WIN, "HBG1")
    assert cand.strand == "+"
    assert cand.proto_5p_genomic == 5_249_937 + 20
    # base at protospacer position 1 equals the plus-strand reference base there
    assert cand.protospacer[0] == plus[20]


def test_offtarget_metadata_attached():
    cand = DESIGNER._row_to_candidate(_row("19rev", "CTTGACCAATAGCCTTGACAAGG"), HBG1_WIN, "HBG1")
    assert cand.offtarget["offtarget_count"] == "100"
    assert cand.offtarget["mit_specificity"] == "45"
    assert cand.offtarget["source"] == "CRISPOR"


def test_editable_base_detection():
    cand = DESIGNER._row_to_candidate(_row("19rev", "CTTGACCAATAGCCTTGACAAGG"), HBG1_WIN, "HBG1")
    assert cand.has_editable_base("ABE")  # has A in window
    assert cand.has_editable_base("CBE")  # has C in window
