"""The checks on the ranking: how it differs from a simple sort, and what is still unmeasured."""
from base_edit_safety.scoring.credibility import annotate_uncertainty, divergence_report
from base_edit_safety.scoring.rubric import score_guides
from base_edit_safety.io.outputs import SCOPE
from base_edit_safety.pipeline.datatypes import JOINT_ON_TARGET_CAVEAT


def _rec(guide_id, offset, role, rate, hist=None, ch=None):
    hist = hist or {}
    return {
        "gene": "HBG1", "editor_class": "ABE", "promoter_offset": offset,
        "promoter_sense_ref": "A", "promoter_sense_alt": "G",
        "on_target_or_bystander": role, "predicted_rate": rate,
        "provenance": {
            "guide_id": guide_id,
            "offtarget": {"offtarget_count": str(sum(hist.values()) + 1), "mit_specificity": "40",
                          "mm_histogram": hist, "cfd_weighted_burden": 0.0, "paralog_count": 1},
            "ch_driver_offtargets": ch or [],
        },
    }


def _ch(gene, mm, feature):
    return {"ch_driver_gene": gene, "mismatch_count": str(mm), "feature": feature,
            "chrom": "chr2", "start": "25264110"}


def test_a_fast_but_gated_guide_gets_moved_down():
    # Sorting by efficiency would put the CH-gated guide first. The ranking must not.
    recs = ([_rec("clean", 113, "on_target", 0.4, {"2": 1})]
            + [_rec("dirty_ch", 113, "on_target", 0.5, {}, ch=[_ch("DNMT3A", 1, "exon")])])
    rep = divergence_report(recs)
    row = {r["guide_id"]: r for r in rep["per_design"]}
    assert row["dirty_ch"]["naive_efficiency_rank"] == 1     # it is the most efficient
    assert row["dirty_ch"]["rubric_rank"] > row["clean"]["rubric_rank"]
    assert row["dirty_ch"]["delta_vs_naive_efficiency"] > 0  # positive means it moved down
    assert "CH-driver" in row["dirty_ch"]["rationale"]
    assert rep["n_changed_vs_naive_efficiency"] >= 1


def test_divergence_reports_count_and_movers():
    recs = [_rec("a", 113, "on_target", 0.4, {"4": 500}),
            _rec("b", 175, "on_target", 0.4, {"2": 3})]
    rep = divergence_report(recs)
    assert rep["n_designs"] == 2
    assert len(rep["biggest_movers_vs_naive_efficiency"]) >= 1


def test_uncertainty_names_the_assay_for_each_gap():
    s = score_guides([_rec("g", 113, "on_target", 0.4, {"2": 1, "4": 100})])[0].as_dict()
    u = annotate_uncertainty(s)
    assert set(u["source_caveats"]) == {
        "on_target_edit_rate", "off_target_candidates", "regulatory_consequence",
        "ch_feature_call", "efficacy_floor"}
    assert "Predicted by BE-DICT" in u["source_caveats"]["on_target_edit_rate"]
    assert any("amplicon NGS" in t for t in u["testable_predictions"])
    assert any("<=2 mismatches" in t for t in u["testable_predictions"])
    assert "not a measured rate" in u["confidence_basis"]


def test_uncertainty_mentions_a_noted_ch_hit():
    s = score_guides([_rec("g", 123, "on_target", 0.19, {}, ch=[_ch("DNMT3A", 4, "intron")])])[0].as_dict()
    u = annotate_uncertainty(s)
    assert any("CH-driver off-target noted" in t for t in u["testable_predictions"])


def test_uncertainty_does_not_mislabel_a_missing_joint_as_outside_the_window():
    records = [
        _rec("pair", 123, "on_target", 0.45),
        _rec("pair", 124, "on_target", 0.35),
    ]
    score = score_guides(records)[0].as_dict()
    testable = " ".join(annotate_uncertainty(score)["testable_predictions"])

    assert "Joint on-target editing" in testable
    assert "joint haplotype probability was unavailable or invalid" in testable
    assert "every intended edit" not in testable
    assert "outside the window" not in testable


def test_uncertainty_labels_a_multi_edit_value_as_joint_and_not_measured():
    records = [
        _rec("pair", 123, "on_target", 0.45),
        _rec("pair", 124, "on_target", 0.35),
    ]
    allele_keys = [
        {"promoter_offset": -123, "promoter_sense_ref": "A", "promoter_sense_alt": "G"},
        {"promoter_offset": -124, "promoter_sense_ref": "A", "promoter_sense_alt": "G"},
    ]
    for record in records:
        record["provenance"].update({
            "joint_on_target_rate": 0.20,
            "joint_on_target_allele_keys": allele_keys,
            "joint_on_target_confidence": "model_predicted",
            "joint_on_target_method": "sum_haplotypes_containing_all_intended_edits",
            "joint_on_target_model": "BE-DICT-bystander:fixture",
            "predictor_runs": 5,
            "joint_on_target_caveat": JOINT_ON_TARGET_CAVEAT,
        })
    score = score_guides(records)[0].as_dict()
    testable = " ".join(annotate_uncertainty(score)["testable_predictions"])

    assert "Model-predicted joint on-target haplotype frequency" in testable
    assert "this is not measured" in testable


def test_uncertainty_never_states_a_risk():
    s = score_guides([_rec("g", 113, "on_target", 0.4, {"2": 1})])[0].as_dict()
    blob = str(annotate_uncertainty(s)).lower().replace("clonal-hematopoiesis", "")
    for bad in ("clonal", "expansion", "leukemia", "lifetime", "risk", "trajectory"):
        assert bad not in blob


def test_scope_note_says_what_this_is_and_is_not():
    assert "shortlist" in SCOPE["what_this_is"]
    joined = " ".join(SCOPE["what_this_is_not"]).lower()
    assert "not a safety clearance" in joined
    assert "not an efficacy ranking" in joined
    assert "not a per-person prediction" in joined
    assert "potency assay" in SCOPE["thresholds"]
