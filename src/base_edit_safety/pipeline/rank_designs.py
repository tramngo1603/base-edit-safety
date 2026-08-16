"""Rank designs, then compare the result with single-variable orderings.

Both steps read the same enumeration file and use the same `RubricWeights` values.
"""
from __future__ import annotations

from ..io import outputs
from ..scoring.credibility import annotate_uncertainty, divergence_report
from ..scoring.rubric import EFFICACY_NOTE, RubricWeights, WEIGHTS_STATUS, score_guides


def _load_enumeration() -> list[dict]:
    return outputs.read_enumeration()["records"]


def run_ranking(record_dicts: list[dict] | None = None,
                weights: RubricWeights | None = None, write: bool = True) -> list[dict]:
    """Return designs ordered by triage tier and the terms used for each score."""
    weights = weights or RubricWeights()
    if record_dicts is None:
        record_dicts = _load_enumeration()

    ranked = [s.as_dict() for s in score_guides(record_dicts, weights)]
    if write:
        outputs.write_ranking(
            ranked,
            filename="ranking_HBG1_HBG2.json",
            meta={
                "step": "5 (safety scoring and ranking)",
                "weights_status": WEIGHTS_STATUS,
                "weights": weights.as_dict(),
                "efficacy_note": EFFICACY_NOTE,
                "bystander_scoring": "by what the bystander does to the motif, not by how many there are",
                "offtarget_scoring": "weighted by mismatch count, so 4-mismatch sites count for little; "
                                     "the 0-mismatch paralog hit is dropped; no hard count cutoff, "
                                     "because that would gate on noise",
                "ch_gate": "a predicted CH-driver off-target drops a design to the bottom tier only if "
                           "it has few mismatches and sits in an exon or UTR. Otherwise it becomes a "
                           "visible note. Either way it means go and measure it.",
                "consequence_term": "penalizes the predicted regulatory effect of bystanders that are "
                                    "not productive, once annotate_consequence.py has filled it in. "
                                    "Enumeration alone leaves consequence 'deferred', so the term is "
                                    "zero until then.",
                "ranking_unit": "one design = one guide and editor that hits a targeted position",
                "score_interpretation": "Read the tier and the `terms` first. safety_score orders "
                                        "designs inside a tier and nothing more.",
            },
        )
    return ranked


def run_credibility(record_dicts: list[dict] | None = None,
                    weights: RubricWeights | None = None, write: bool = True) -> dict:
    """Compare the ranking to simple sorts, and note what is still unmeasured per design."""
    weights = weights or RubricWeights()
    if record_dicts is None:
        record_dicts = _load_enumeration()

    scores = [s.as_dict() for s in score_guides(record_dicts, weights)]
    divergence = divergence_report(record_dicts, weights)
    annotated = [{**s, "uncertainty": annotate_uncertainty(s)} for s in scores]

    body = {
        "divergence": divergence,
        "ranked_designs_with_uncertainty": annotated,
    }
    if write:
        outputs.write_credibility(
            body,
            filename="credibility_HBG1_HBG2.json",
            meta={
                "step": "6 (checks on the ranking)",
                "weights_status": WEIGHTS_STATUS,
                "simple_sorts_compared": ["on-target efficiency alone", "raw off-target count alone"],
            },
        )
    return body


if __name__ == "__main__":
    ranked = run_ranking()
    print(f"weights: {WEIGHTS_STATUS}")
    print(f"ranked designs: {len(ranked)}\n")
    hdr = "{:>3} {:<26} {:<18} {:<13} {:>6} {:>8} {:>9} {:>10} {:>8}".format(
        "#", "tier", "guide", "anchors", "eff", "bystPen", "offtPen", "CHflag", "score")
    print(hdr); print("-" * len(hdr))
    for s in ranked:
        chf = (",".join(s["ch_genes"]) + ("*GATE" if s["ch_gate"] else ("*noted" if s["ch_flag_noted"] else ""))) if s["ch_genes"] else "-"
        print("{:>3} {:<26} {:<18} {:<13} {:>6} {:>8} {:>9} {:>10} {:>8}".format(
            s["rank"], s["tier_label"], s["guide_id"], str(s["anchors"]),
            s["on_target_efficiency"], s["terms"]["bystander_penalty"],
            s["terms"]["offtarget_penalty"], chf, s["safety_score"]))

    body = run_credibility()
    d = body["divergence"]
    print(f"\n=== how the ranking differs from a simple sort: {d['n_designs']} designs ===")
    print(f"moved vs sorting by efficiency: {d['n_changed_vs_naive_efficiency']} | "
          f"moved vs sorting by off-target count: {d['n_changed_vs_naive_offtarget_count']}\n")
    print("biggest moves vs sorting by efficiency (rank here <- rank by efficiency):")
    for m in d["biggest_movers_vs_naive_efficiency"]:
        print(f"  {m['guide_id']:<20} rank#{m['rubric_rank']:<3} byEff#{m['naive_efficiency_rank']:<3} "
              f"byCount#{m['naive_offtarget_count_rank']:<3} | {m['rationale']}")
