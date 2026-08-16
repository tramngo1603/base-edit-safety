"""Test functional bystander terms, efficiency tiers, off-target weights, and CH flags."""
from base_edit_safety.bystanders import UNKNOWN, allele_context, bystander_function
from base_edit_safety.scoring.rubric import RubricWeights, WEIGHTS_STATUS, score_guides


def _rec(guide_id, offset, role, rate, hist=None, ch=None, sense="A>G"):
    hist = hist or {}
    sref, salt = sense.split(">")
    return {
        "gene": "HBG1", "editor_class": "ABE", "promoter_offset": offset,
        "promoter_sense_ref": sref, "promoter_sense_alt": salt,
        "on_target_or_bystander": role, "predicted_rate": rate,
        "confidence": "not_evaluated" if rate is None else "model_predicted",
        "provenance": {
            "guide_id": guide_id,
            "offtarget": {
                "offtarget_count": str(sum(hist.values()) + 1),  # +1 paralog
                "mit_specificity": "40",
                "mm_histogram": hist,
                "cfd_weighted_burden": 0.0,
                "paralog_count": 1,
            },
            "ch_driver_offtargets": ch or [],
        },
    }


def _ch(gene, mm, feature):
    return {"ch_driver_gene": gene, "mismatch_count": str(mm), "feature": feature,
            "chrom": "chr2", "start": "25264110"}


def test_open_weights_are_listed():
    assert "Starting values" in WEIGHTS_STATUS and RubricWeights().provisional is True
    knobs = RubricWeights().as_dict()["weights_to_review"]
    assert any("off_w_mm" in k for k in knobs)
    assert any("ch_gate_mm_threshold" in k for k in knobs)


ON_113 = {(-113, "A", "G"): "present"}
ON_175 = {(-175, "A", "G"): "present"}


def test_lookup_needs_the_allele_not_just_the_offset():
    """-110 is an ABE bystander as A>G and an out-of-scope NF-Y transversion as A>C."""
    assert bystander_function("HBG1", 110, "A", "G", ON_113)[0] == "adverse"
    assert bystander_function("HBG1", 110, "A", "C", ON_113)[0] == UNKNOWN
    assert bystander_function("HBG1", 999, "A", "G", ON_113)[0] == UNKNOWN


def test_classes_are_conditional_on_the_rest_of_the_allele():
    """Require the -113 edit before applying conditional GATA-site bystander rules."""
    for offset, expected in ((112, "adverse"), (110, "adverse"), (116, "weakly_productive")):
        assert bystander_function("HBG1", offset, "A", "G", ON_113)[0] == expected
        # same edit, no -113 on the allele: no claim
        assert bystander_function("HBG1", offset, "A", "G", {})[0] == UNKNOWN


def test_181_is_neutral_only_while_173_is_unedited():
    """Mayuranathan: A3G5A11 and A3G5G11 are both productive, G3G5G11 is not."""
    assert bystander_function("HBG1", 181, "A", "G", ON_175)[0] == "neutral"
    both = {**ON_175, (-173, "A", "G"): "present"}
    assert bystander_function("HBG1", 181, "A", "G", both)[0] == UNKNOWN
    # and the harm is attributed to -173 itself
    assert bystander_function("HBG1", 173, "A", "G", ON_175)[0] == "adverse"


def test_an_unresolved_dependency_declines_rather_than_guessing():
    """Leave -112 unresolved when the required -113 edit is unresolved."""
    unresolved = {(-113, "A", "G"): "unresolved"}
    cls, note = bystander_function("HBG1", 112, "A", "G", unresolved)
    assert cls == UNKNOWN and "unresolved" in note


def test_allele_context_treats_a_missing_rate_as_unresolved():
    ctx = allele_context([
        _rec("g", 113, "on_target", 0.4),          # edited
        _rec("g", 112, "bystander", 0.0),          # predicted zero: really absent
        _rec("g", 110, "bystander", None),         # never asked: unresolved
    ])
    assert ctx[(-113, "A", "G")] == "present"
    assert ctx[(-112, "A", "G")] == "absent"
    assert ctx[(-110, "A", "G")] == "unresolved"


def test_function_not_count_drives_bystander_penalty():
    adv = score_guides([_rec("adv", 113, "on_target", 0.4, {}),
                        _rec("adv", 112, "bystander", 0.2, {})])[0]
    prod = score_guides([_rec("prod", 113, "on_target", 0.4, {}),
                         _rec("prod", 116, "bystander", 0.2, {})])[0]
    assert adv.bystander_penalty > 0 and prod.bystander_penalty == 0.0
    assert adv.bystander_breakdown.get("adverse") == 0.2
    assert prod.bystander_breakdown.get("weakly_productive") == 0.2


def test_offtarget_penalty_is_mismatch_weighted_not_flat():
    # 200 sites at 4 mismatches should cost far less than 2 sites at 2 mismatches
    many_4mm = score_guides([_rec("far", 113, "on_target", 0.4, {"4": 200})])[0]
    few_2mm = score_guides([_rec("near", 113, "on_target", 0.4, {"2": 2})])[0]
    assert few_2mm.terms["offtarget_penalty"] > many_4mm.terms["offtarget_penalty"]
    # the raw counts stay visible, with the paralog hit already dropped
    assert many_4mm.terms["offtarget_mm_histogram"] == {"4": 200}
    assert many_4mm.terms["offtarget_paralog_excluded"] is True


def test_offtarget_weights_configurable():
    recs = [_rec("g", 113, "on_target", 0.4, {"4": 100})]
    base = score_guides(recs)[0].terms["offtarget_penalty"]
    heavier = score_guides(recs, RubricWeights(off_w_mm4=0.1))[0].terms["offtarget_penalty"]
    assert heavier > base


def test_ch_gate_only_fires_on_a_likely_exonic_hit():
    recs = [_rec("clean", 113, "on_target", 0.4, {}),
            _rec("exonic", 113, "on_target", 0.49, {}, ch=[_ch("DNMT3A", 1, "exon")])]
    by = {s.guide_id: s for s in score_guides(recs)}
    assert by["exonic"].ch_gate is True
    assert by["exonic"].tier_label == "CH_DRIVER_FLAGGED"
    assert by["clean"].rank < by["exonic"].rank


def test_intronic_or_distant_ch_hit_becomes_a_note():
    # The recorded 61forw/DNMT3A case is intronic with four mismatches, so it is noted but not gated.
    s = score_guides([_rec("g", 123, "on_target", 0.19, {}, ch=[_ch("DNMT3A", 4, "intron")])])[0]
    assert s.ch_gate is False
    assert s.ch_flag_noted is True
    assert "not gating" in s.ch_flag_reason
    assert s.tier_label != "CH_DRIVER_FLAGGED"
    assert s.efficacy_band == "PRECLINICAL_ACTIVITY_ONLY"


def test_ch_gate_threshold_is_tunable():
    rec = [_rec("g", 113, "on_target", 0.4, {}, ch=[_ch("DNMT3A", 3, "exon")])]
    assert score_guides(rec)[0].ch_gate is True                              # default threshold is 3
    assert score_guides(rec, RubricWeights(ch_gate_mm_threshold=2))[0].ch_gate is False


def test_efficacy_tiers_and_reward():
    def band(e):
        return score_guides([_rec("g", 113, "on_target", e, {})])[0].efficacy_band
    assert band(0.7) == "STRONG_CLINICAL_RANGE"
    assert band(0.4) == "MEETS_CLINICAL_PRECEDENT"
    assert band(0.2) == "PRECLINICAL_ACTIVITY_ONLY"
    assert band(0.1) == "BELOW_EVIDENCED_ACTIVITY"
    hi = score_guides([_rec("hi", 113, "on_target", 0.55, {}),
                       _rec("lo", 175, "on_target", 0.32, {})])
    assert {s.guide_id: s.rank for s in hi}["hi"] == 1


def test_all_intended_edits_not_evaluated_return_a_score():
    score = score_guides([
        _rec("unknown", 113, "on_target", None, {}),
        _rec("unknown", 175, "on_target", None, {}),
    ])[0]

    assert score.on_target_efficiency is None
    assert score.efficacy_band == "NOT_EVALUATED"
    assert score.meets_clinical_floor is None
    assert score.tier == RubricWeights().tier_not_evaluated
    assert score.tier_label == "NOT_EVALUATED"
    assert score.terms["unevaluated_anchors"] == 2


def test_evaluated_zero_is_not_treated_as_not_evaluated():
    score = score_guides([_rec("zero", 113, "on_target", 0.0, {})])[0]

    assert score.on_target_efficiency == 0.0
    assert score.efficacy_band == "BELOW_EVIDENCED_ACTIVITY"
    assert score.meets_clinical_floor is False
    assert score.tier == RubricWeights().tier_below_activity
    assert score.terms["unevaluated_anchors"] == 0


def test_same_input_gives_the_same_scores_and_no_risk_field():
    recs = [_rec("a", 113, "on_target", 0.4, {"4": 10}),
            _rec("a", 112, "bystander", 0.1, {"4": 10})]
    assert [s.safety_score for s in score_guides(recs)] == [s.safety_score for s in score_guides(recs)]
    d = score_guides(recs)[0].as_dict()
    forbidden = ("clonal", "expansion", "leukemia", "lifetime", "risk", "trajectory")
    for key in d:
        assert not any(bad in key.lower() for bad in forbidden)


def test_guides_without_anchor_excluded():
    assert score_guides([_rec("byst_only", 116, "bystander", 0.2, {})]) == []


def _rec_cons(gid, off, role, rate, sev):
    r = _rec(gid, off, role, rate, {})
    r["consequence"] = {"severity": sev, "classification": "functionally-consequential",
                        "source": "AlphaGenome"}
    return r


def test_consequence_penalizes_nonproductive_bystanders_only():
    # A bystander that breaks the motif (-112) is penalized; one that helps (-116) is not.
    adv = score_guides([_rec_cons("a", 113, "on_target", 0.4, 0.7),
                        _rec_cons("a", 112, "bystander", 0.2, 0.7)])[0]
    prod = score_guides([_rec_cons("p", 113, "on_target", 0.4, 0.7),
                         _rec_cons("p", 116, "bystander", 0.2, 0.7)])[0]
    # -116 is weakly_productive, which is still the point of the edit, so no consequence penalty
    assert adv.terms["consequence_penalty"] > 0
    assert prod.terms["consequence_penalty"] == 0.0


def _ch_fc(mm, classification, feature="intron"):
    return {"ch_driver_gene": "DNMT3A", "mismatch_count": str(mm), "feature": feature,
            "chrom": "chr2", "start": "25264110",
            "ch_functional_consequence": {"classification": classification, "severity": 0.6}}


def test_ch_gate_uses_the_alphagenome_call():
    # An accepted functional class and one mismatch select the CH-flagged tier despite an intron label.
    g1 = score_guides([_rec("g1", 113, "on_target", 0.4, {}, ch=[_ch_fc(1, "functionally-consequential")])])[0]
    assert g1.ch_gate is True
    # One mismatch with a low-consequence class remains a note and does not change the tier.
    g2 = score_guides([_rec("g2", 113, "on_target", 0.4, {}, ch=[_ch_fc(1, "low-consequence")])])[0]
    assert g2.ch_gate is False and g2.ch_flag_noted is True
    # Four mismatches exceed the default limit of three, so the site remains a note.
    g3 = score_guides([_rec("g3", 113, "on_target", 0.4, {}, ch=[_ch_fc(4, "functionally-consequential")])])[0]
    assert g3.ch_gate is False and g3.ch_flag_noted is True
    assert g3.ch_functional_consequence["classification"] == "functionally-consequential"


def test_the_stand_in_call_does_not_gate_on_its_own():
    # An intronic hit whose only functional signal is the stand-in A>G call falls back to the
    # coordinate, so it does not change tier unless configuration permits stand-in scores.
    hit = {"ch_driver_gene": "DNMT3A", "mismatch_count": "1", "feature": "intron",
           "chrom": "chr2", "start": "25264110",
           "ch_functional_consequence": {"classification": "functionally-consequential",
                                         "severity": 0.7, "is_proxy": True}}
    rec = [_rec("g", 113, "on_target", 0.4, {}, ch=[hit])]
    assert score_guides(rec)[0].ch_gate is False
    assert score_guides(rec, RubricWeights(ch_gate_trust_proxy_consequence=True))[0].ch_gate is True
