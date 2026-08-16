"""The CH layer decides gene membership by coordinate, and produces no outcome number.

It does emit coordinate, mismatch count, feature and mechanism alongside the gene. Those are
deliberate; see the field list in annotation.py. What it must never emit is a risk, an
incidence, a time-to-event, or a per-person number.
"""
from base_edit_safety.annotation import CH_DRIVER_GENES, annotate_record, ch_driver_offtargets
from base_edit_safety.pipeline.datatypes import EditRecord


def _record(chrom, pos, gene="HBG1"):
    return EditRecord(
        genome_build="GRCh38", chrom=chrom, pos_1based=pos,
        genomic_plus_ref="T", genomic_plus_alt="C",
        promoter_sense_ref="A", promoter_sense_alt="G",
        gene=gene, gene_strand="-", promoter_offset=113,
        editor_class="ABE", on_target_or_bystander="on_target",
    )


def test_driver_panel_is_the_twelve_genes():
    assert CH_DRIVER_GENES == {
        "ASXL1", "CHEK2", "DNMT3A", "GNAS", "IDH2", "MDM4", "PPM1D", "SF3B1", "SRSF2",
        "TET2", "TP53", "U2AF1",
    }


def test_promoter_edit_is_not_flagged():
    r = annotate_record(_record("chr11", 5_249_970))  # HBG1 -113
    assert r.ch_driver_flag is False
    assert r.ch_driver_gene is None
    assert "ch_driver_set" in r.provenance


def test_edit_inside_ch_gene_is_flagged():
    r = annotate_record(_record("chr2", 25_264_110, gene="DNMT3A"))  # within DNMT3A span
    assert r.ch_driver_flag is True
    assert r.ch_driver_gene == "DNMT3A"


def test_offtarget_in_ch_gene_detected():
    sites = [{"chrom": "chr2", "start": "25264110", "strand": "+", "mismatchCount": "4",
              "cfdOfftargetScore": "0.1", "locusDesc": "exon:DNMT3A"}]
    hits = ch_driver_offtargets(sites)
    assert len(hits) == 1
    assert hits[0]["ch_driver_gene"] == "DNMT3A"
    assert hits[0]["boostdm_ch_score"] is None  # per-mutation scores need the CC-BY-NC download


def test_site_near_a_ch_gene_is_not_flagged():
    # This site sits outside the TP53 span, and the tool names it "intergenic:TP53-WRAP53".
    # The name mentions TP53, but the coordinate decides, so it stays unflagged.
    sites = [{"chrom": "chr17", "start": "7700000", "strand": "+", "mismatchCount": "3",
              "locusDesc": "intergenic:TP53-WRAP53"}]
    assert ch_driver_offtargets(sites) == []


def test_ch_output_carries_no_risk_field():
    sites = [{"chrom": "chr2", "start": "25264110", "locusDesc": "exon:DNMT3A"}]
    hit = ch_driver_offtargets(sites)[0]
    forbidden = ("clonal", "expansion", "leukemia", "lifetime", "risk", "trajectory")
    for key in hit:
        assert not any(bad in key.lower() for bad in forbidden)
