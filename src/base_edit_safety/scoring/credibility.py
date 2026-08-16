"""Describe ranking differences and unresolved measurements.

`divergence_report` compares the rubric ordering with two defined single-variable orderings: editing
efficiency and raw off-target count. Differences show where other rubric terms change the order.

`annotate_uncertainty` lists unmeasured fields and a named assay for each design. This module does
not calculate clonal expansion or long-term outcomes.
"""
from __future__ import annotations

from .rubric import GuideScore, RubricWeights, score_guides


# --- What is not known, per design ----------------------------------------------------------------

# These hold for every design, because they are limits of the prediction sources themselves.
SOURCE_CAVEATS = {
    "on_target_edit_rate": (
        "Predicted by BE-DICT, which was trained on cell lines and libraries, not on hematopoietic "
        "stem cells. A multi-edit value is a model-predicted joint haplotype frequency, not a "
        "measurement. Amplicon NGS in CD34+ HSPCs would measure it."
    ),
    "off_target_candidates": (
        "Candidate sites from a CRISPOR search of the reference genome, up to 4 mismatches, with no "
        "bulges and no population variants. Most of the count is 4-mismatch sites, which receive a "
        "lower rubric weight. CHANGE-seq-BE or GUIDE-seq would provide measured off-target evidence."
    ),
    "regulatory_consequence": (
        "Predicted by AlphaGenome (non-commercial), and only as a magnitude of regulatory "
        "disturbance: no direction and no clinical outcome. A reporter or functional assay would "
        "measure the effect."
    ),
    "ch_feature_call": (
        "Whether a site is exonic, intronic, or intergenic comes from its coordinate. What that "
        "site does comes from AlphaGenome. A noted CH off-target needs targeted sequencing plus a "
        "functional check."
    ),
    "efficacy_floor": (
        "The tier boundaries come from allelic editing in mixed blood and marrow cells, and an "
        "allelic percentage does not convert into a fraction of edited HSCs. Re-estimate with the "
        "final guide, editor, and process, using a measured potency assay."
    ),
}


def annotate_uncertainty(score: dict) -> dict:
    """List what one ranked design still needs measured, and how."""
    testable = []

    eff = score["on_target_efficiency"]
    band = score["efficacy_band"]
    efficiency_basis = score.get("terms", {}).get("efficiency_basis")
    if eff is None:
        if efficiency_basis == "joint_on_target_rate_not_evaluated":
            testable.append(
                f"Joint on-target editing at {score['anchors']} was not evaluated: a complete, "
                f"allele-matched joint haplotype probability was unavailable or invalid, and the "
                f"marginal rates cannot substitute for it. Amplicon NGS would measure the joint "
                f"haplotype; predictor and source review would resolve the missing output."
            )
        else:
            testable.append(
                f"On-target editing at {score['anchors']} was not predicted: every intended edit "
                f"for this design sits outside the requested prediction window, so there is "
                f"no rate to confirm or refute. Amplicon NGS would measure it, but the gap is in the "
                f"tool."
            )
    elif efficiency_basis == "joint_on_target_rate":
        testable.append(
            f"Model-predicted joint on-target haplotype frequency ~{round(eff * 100, 1)}% at "
            f"{score['anchors']} (band: {band}); this is not measured, and amplicon NGS in the "
            f"target cells would measure it."
        )
    else:
        testable.append(
            f"Predicted on-target editing ~{round(eff * 100, 1)}% at {score['anchors']} "
            f"(band: {band}); amplicon NGS in the target cells would measure it."
        )

    hist = score.get("offtarget_mm_histogram", {}) or {}
    low_mm = {m: n for m, n in hist.items() if int(m) <= 2}
    if low_mm:
        testable.append(
            f"Most likely off-targets: {sum(low_mm.values())} site(s) at <=2 mismatches "
            f"({low_mm}); CHANGE-seq-BE or GUIDE-seq would test these sites."
        )
    testable.append(
        f"Off-target sites are weighted by how likely they are; {hist.get('4', 0)} of them have 4 "
        f"mismatches; an assay is required before treating any site as an observed edit."
    )

    cons_sev = score.get("terms", {}).get("consequence_severity_nonproductive_bystanders")
    if cons_sev:
        testable.append(
            f"Predicted regulatory effect of the bystanders that are not productive (AlphaGenome, "
            f"severity-weighted {round(cons_sev, 3)}); a reporter assay would measure the effect."
        )

    if score.get("ch_gate"):
        testable.append(
            f"A CH-driver off-target in {score['ch_genes']} drops this design to the bottom tier "
            f"(predicted, few mismatches, exonic); it needs off-target sequencing and a functional check."
        )
    elif score.get("ch_flag_noted"):
        testable.append(
            f"CH-driver off-target noted ({score['ch_genes']}) but not gating "
            f"({score.get('ch_flag_reason')}); confirm the coordinate, the mismatches, and the effect."
        )

    return {
        "source_caveats": SOURCE_CAVEATS,
        "testable_predictions": testable,
        "confidence_basis": (
            "Every term is a model prediction or a candidate site from a reference-genome search: "
            "not a measured rate, not a confirmed off-target, and not a per-person estimate."
        ),
    }


# --- How the ranking compares to a simple sort ----------------------------------------------------

def _ranked_ids(scores: list[GuideScore], key, reverse: bool) -> list[str]:
    ordered = sorted(scores, key=lambda s: (key(s), s.guide_id), reverse=reverse)
    # Ties break on guide_id, so the order is the same on every run.
    return [s.guide_id for s in ordered]


def _rationale(s: GuideScore, eff_delta: int) -> str:
    """One line on why this design landed where it did, compared to sorting by efficiency."""
    bits = []
    if s.ch_gate:
        bits.append(f"CH-driver off-target gates it (few mismatches, exonic: {','.join(s.ch_genes)})")
    elif s.ch_flag_noted:
        bits.append(f"CH-driver off-target noted ({','.join(s.ch_genes)}, not gating)")
    if s.tier_label == "BELOW_EVIDENCED_ACTIVITY":
        bits.append("below the clinical efficiency anchor")
    elif s.tier_label == "PRECLINICAL_ACTIVITY_ONLY":
        bits.append("preclinical activity only")
    hist = s.offtarget_mm_histogram or {}
    low_mm = sum(n for m, n in hist.items() if int(m) <= 2)
    if low_mm:
        bits.append(f"{low_mm} off-target(s) at <=2 mismatches")
    if s.bystander_penalty > 0:
        adv = s.bystander_breakdown.get("adverse", 0)
        if adv:
            bits.append("a bystander that breaks the motif")
    direction = "moved down for safety" if eff_delta > 0 else ("moved up" if eff_delta < 0 else "unchanged")
    return f"{direction}: " + ("; ".join(bits) if bits else "the safety terms reorder it")


def divergence_report(records: list[dict], weights: RubricWeights | None = None) -> dict:
    """Compare the ranking against sorting by efficiency and by raw off-target count."""
    scores = score_guides(records, weights or RubricWeights())
    rubric_rank = {s.guide_id: s.rank for s in scores}

    # A design with no efficiency cannot be ranked by efficiency; it sorts last rather than being
    # read as zero.
    eff_order = _ranked_ids(scores, reverse=True,
                            key=lambda s: (s.on_target_efficiency is not None,
                                           s.on_target_efficiency or 0.0))
    eff_rank = {gid: i + 1 for i, gid in enumerate(eff_order)}

    # Sorting by raw count places designs with fewer candidate sites earlier.
    cnt_order = _ranked_ids(scores, key=lambda s: s.offtarget_count, reverse=False)
    cnt_rank = {gid: i + 1 for i, gid in enumerate(cnt_order)}

    rows = []
    for s in scores:
        # A positive delta means the rubric ranks it worse than the simple sort would.
        eff_delta = rubric_rank[s.guide_id] - eff_rank[s.guide_id]
        cnt_delta = rubric_rank[s.guide_id] - cnt_rank[s.guide_id]
        rows.append({
            "guide_id": s.guide_id,
            "rubric_rank": rubric_rank[s.guide_id],
            "naive_efficiency_rank": eff_rank[s.guide_id],
            "naive_offtarget_count_rank": cnt_rank[s.guide_id],
            "delta_vs_naive_efficiency": eff_delta,
            "delta_vs_naive_count": cnt_delta,
            "tier": s.tier_label,
            "on_target_efficiency": s.on_target_efficiency,
            "offtarget_count": s.offtarget_count,
            "offtarget_penalty": s.terms["offtarget_penalty"],
            "rationale": _rationale(s, eff_delta),
        })

    changed_eff = [r for r in rows if r["delta_vs_naive_efficiency"] != 0]
    changed_cnt = [r for r in rows if r["delta_vs_naive_count"] != 0]
    movers = sorted(rows, key=lambda r: abs(r["delta_vs_naive_efficiency"]), reverse=True)[:6]

    return {
        "n_designs": len(rows),
        "n_changed_vs_naive_efficiency": len(changed_eff),
        "n_changed_vs_naive_offtarget_count": len(changed_cnt),
        "biggest_movers_vs_naive_efficiency": movers,
        "per_design": rows,
        "note": ("The simple sorts each look at one number. The rubric puts safety tiers first, "
                 "judges bystanders by what they do, weights off-targets by how likely they are, "
                 "and gates on CH drivers. Where the orders differ is where that matters."),
    }
