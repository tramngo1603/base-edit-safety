"""The anchor table loads with both allele orientations, and the spacing check works."""
import pytest

from base_edit_safety import anchors
from base_edit_safety.anchors import load_in_scope_records


def test_spacing_holds():
    anchors.check_paralog_spacing()  # must not raise


def test_wrong_spacing_stops_the_run(monkeypatch):
    # A wrong TSS separation has to raise, not slip through.
    monkeypatch.setattr(anchors, "HBG2_TSS", anchors.HBG1_TSS + anchors.PARALOG_SPACING_BP + 1)
    with pytest.raises(AssertionError):
        anchors.check_paralog_spacing()


def test_loader_emits_both_paralogs_for_every_position():
    records = load_in_scope_records()
    assert len(records) == 2 * len(anchors.IN_SCOPE)
    genes = {r.gene for r in records}
    assert genes == {"HBG1", "HBG2"}


def test_dual_orientation_populated_everywhere():
    for r in load_in_scope_records():
        assert r.genomic_plus_ref and r.genomic_plus_alt
        assert r.promoter_sense_ref and r.promoter_sense_alt
        # sense base is the complement of the genomic-plus base for these minus-strand genes
        comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
        assert r.promoter_sense_ref == comp[r.genomic_plus_ref]
        assert r.promoter_sense_alt == comp[r.genomic_plus_alt]


def test_coordinates_match_the_doc():
    by_key = {(r.gene, r.promoter_offset): r for r in load_in_scope_records()}
    # HBG2 -114 = 5,254,895, genomic G>A / sense C>T
    r = by_key[("HBG2", 114)]
    assert r.pos_1based == 5_254_895
    assert (r.genomic_plus_ref, r.genomic_plus_alt) == ("G", "A")
    assert (r.promoter_sense_ref, r.promoter_sense_alt) == ("C", "T")
    # HBG1 -113 = 5,249,970, genomic T>C / sense A>G
    r = by_key[("HBG1", 113)]
    assert r.pos_1based == 5_249_970
    assert (r.genomic_plus_ref, r.genomic_plus_alt) == ("T", "C")
    assert (r.promoter_sense_ref, r.promoter_sense_alt) == ("A", "G")


def test_variant_db_status_is_per_gene():
    by_key = {(r.gene, r.promoter_offset): r for r in load_in_scope_records()}
    # HBG2 -175 has a record (rs63750654). HBG1 -175 does not, and must not borrow it.
    assert by_key[("HBG2", 175)].variant_db_status == anchors.DB
    assert by_key[("HBG1", 175)].variant_db_status == anchors.UNVERIFIED
    # An rsID belongs to one gene and is never copied to the other. HBG1 -175 has its own dbSNP id;
    # what it must never carry is HBG2's.
    assert by_key[("HBG2", 175)].provenance.get("dbsnp") == "rs63750654"
    assert by_key[("HBG1", 175)].provenance.get("dbsnp") == "rs34153902"
    assert by_key[("HBG1", 175)].provenance.get("dbsnp") != "rs63750654"
    # A dbSNP id on an unverified row is marked as such, so it cannot be read as confirmation.
    assert by_key[("HBG1", 175)].variant_db_status == anchors.UNVERIFIED
    assert "no ClinVar record" in by_key[("HBG1", 175)].provenance["dbsnp_note"]
    assert "dbsnp_note" not in by_key[("HBG2", 175)].provenance   # HBG2 -175 is db_confirmed

    # Neither gene has a -113 database record, so neither status can claim one.
    assert by_key[("HBG1", 113)].variant_db_status == anchors.UNVERIFIED
    assert by_key[("HBG2", 113)].variant_db_status == anchors.UNVERIFIED
    assert "dbsnp" not in by_key[("HBG2", 113)].provenance
    # The screen hits are marked as such.
    assert by_key[("HBG1", 124)].variant_db_status == anchors.SYNTHETIC
    assert by_key[("HBG1", 123)].variant_db_status == anchors.SYNTHETIC
