"""Turning a predicted spectrum into records: rates, on-target vs bystander, and allele strands.

No external predictor and no network: the guide and the outcomes are built by hand.
"""
from base_edit_safety.pipeline.datatypes import (
    MODEL_PREDICTED, NOT_EVALUATED, EditSpectrum, Guide, OutcomeHaplotype, spectrum_to_records,
)

# Validated HBG1 -113 ABE protospacer (minus strand), 5' end at genomic 5,249,977.
PROTO = "CTTGACCAATAGCCTTGACA"
GUIDE = Guide(
    guide_id="HBG1_m113_ABE",
    protospacer=PROTO,
    pam="AGG",
    strand="-",
    genome_build="GRCh38",
    chrom="chr11",
    gene="HBG1",
    gene_strand="-",
    tss_1based=5_249_857,
    editor_class="ABE",
    editor_name="ABE8e",
    edit_from="A",
    edit_to="G",
    proto_5p_genomic=5_249_977,
)


def _spectrum():
    # pos5 edited (A->G) at 0.20; pos8 (target -113) at 0.10; both at 0.05; unedited 0.65.
    def mut(seq, pos):
        return seq[: pos - 1] + "G" + seq[pos:]
    haps = [
        OutcomeHaplotype(PROTO, 0.65),
        OutcomeHaplotype(mut(PROTO, 5), 0.20),
        OutcomeHaplotype(mut(PROTO, 8), 0.10),
        OutcomeHaplotype(mut(mut(PROTO, 5), 8), 0.05),
    ]
    return EditSpectrum(guide=GUIDE, model_name="test", haplotypes=haps,
                        provenance={"ensemble_runs": 5})


def test_genomic_mapping_minus_strand():
    # protospacer position 8 maps to genomic 5,249,970 (HBG1 -113)
    assert GUIDE.genomic_coord(8) == 5_249_970
    assert GUIDE.promoter_offset(8) == 113


def test_on_target_alleles_match_the_coordinate_table():
    a = GUIDE.alleles_for(8)
    assert (a["genomic_plus_ref"], a["genomic_plus_alt"]) == ("T", "C")
    assert (a["promoter_sense_ref"], a["promoter_sense_alt"]) == ("A", "G")


def test_marginals_sum_over_haplotypes():
    rates = _spectrum().per_position_rate()
    assert abs(rates[5] - 0.25) < 1e-9   # 0.20 + 0.05
    assert abs(rates[8] - 0.15) < 1e-9   # 0.10 + 0.05


def test_reduction_classifies_on_target_and_bystander():
    records = spectrum_to_records(_spectrum(), on_target_offset=113)
    roles = {r.promoter_offset: r.on_target_or_bystander for r in records}
    assert roles[113] == "on_target"
    # position 5 -> genomic 5,249,973 -> offset 116 -> bystander
    assert roles[116] == "bystander"
    # Inside the window, emit positions with a nonzero predicted rate plus the intended position.
    inside = [r for r in records if r.confidence == "model_predicted"]
    assert all(r.predicted_rate > 0 or r.promoter_offset == 113 for r in inside)


def test_positions_outside_the_window_are_reported_as_not_evaluated():
    """Report positions outside EDIT_WINDOW as not evaluated.

    These positions retain a record so a missing prediction is not confused with a predicted zero.
    """
    records = spectrum_to_records(_spectrum(), on_target_offset=113)
    by_pos = {r.provenance["protospacer_position"]: r for r in records}
    for pos in (11, 18, 20):   # editable A's on this protospacer, all past the window edge at 10
        assert pos in by_pos, f"protospacer position {pos} produced no record at all"
        r = by_pos[pos]
        assert r.confidence == NOT_EVALUATED
        assert r.predicted_rate is None
    # In-window positions with predicted edits retain their predicted status.
    for pos in (5, 8):
        assert by_pos[pos].confidence == MODEL_PREDICTED
        assert by_pos[pos].predicted_rate is not None
    # Position 9 is an editable base inside the window with predicted rate zero. Omit this bystander,
    # while retaining positions that were outside the prediction request.
    assert 9 not in by_pos


def test_real_predicted_zero_inside_the_window_is_model_predicted():
    """An evaluated zero is retained for an intended edit and is not a missing prediction."""
    records = spectrum_to_records(_spectrum(), on_target_offset=112)
    by_pos = {r.provenance["protospacer_position"]: r for r in records}

    zero = by_pos[9]
    assert zero.on_target_or_bystander == "on_target"
    assert zero.predicted_rate == 0.0
    assert zero.confidence == MODEL_PREDICTED

    outside = by_pos[11]
    assert outside.predicted_rate is None
    assert outside.confidence == NOT_EVALUATED


def test_records_carry_provenance_and_no_risk_field():
    records = spectrum_to_records(_spectrum(), on_target_offset=113)
    for r in records:
        assert r.model_name == "test"
        assert "guide_id" in r.provenance
        assert not hasattr(r, "lifetime_risk")
        assert not hasattr(r, "clonal_expansion")
