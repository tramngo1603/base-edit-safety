"""Assign a review order without making a clinical safety claim.

The order uses four inputs: the known effect of each bystander edit, predicted editing efficiency,
sequence-matched off-target sites, and overlap with a listed CH driver gene. A CH-gene overlap sends
a design to a separate review tier until laboratory evidence resolves it. It does not predict
clonal expansion or a long-term outcome. Each result includes the values that produced its rank.

The numeric weights are starting values. They remain marked as provisional until the validation
protocol locks them.
"""
from __future__ import annotations

from dataclasses import dataclass, field, asdict

from ..bystanders import (
    ADVERSE, NEUTRAL, PRODUCTIVE, UNKNOWN, WEAKLY_PRODUCTIVE,
    allele_context, bystander_function,
)
from ..pipeline.datatypes import JOINT_ON_TARGET_CAVEAT, JOINT_ON_TARGET_METHOD

# The efficiency thresholds and other weights remain provisional until source review and calibration
# are complete. Their status is carried into generated records.
WEIGHTS_STATUS = (
    "Starting values informed by published data. They are not protocol-locked and remain listed "
    "under weights_to_review."
)

# The efficiency threshold supports comparison and is not a clinical criterion.
EFFICACY_NOTE = (
    "Efficiency values come from BE-DICT, which was trained on cell lines and libraries. The tier "
    "boundaries come from allelic editing measured in mixed blood and marrow cells, not sorted "
    "LT-HSCs, and an allelic percentage does not convert into a fraction of edited HSCs. Use the "
    "floor to compare designs against each other. Multi-edit efficiency is the model-predicted "
    "joint haplotype frequency, not a measurement."
)
NOT_EVALUATED_BAND = "NOT_EVALUATED"

# Bystander rules use the gene, position, exact DNA change, and any required companion edits. A
# position alone cannot distinguish two different DNA changes.


@dataclass(frozen=True)
class RubricWeights:
    """Store every configurable ranking parameter in one read-only data object."""

    # ===== Provisional efficiency boundaries. See EFFICACY_NOTE. ================================
    clinical_floor: float = 0.30      # comparison boundary; not a clinical minimum
    preclinical_floor: float = 0.16   # lower comparison boundary
    strong_band_low: float = 0.60     # upper reporting band; not a clinical claim
    strong_band_high: float = 0.80
    # ===== Bystander weights. These values require protocol review. =============================
    w_bystander: float = 1.0          # overall multiplier on the bystander penalty
    penalty_adverse: float = 1.0      # per unit of edit probability at a motif-breaking position
    penalty_unknown: float = 0.3      # weight for a bystander with no recorded functional class
    penalty_productive: float = 0.0   # productive bystanders are not penalized
    # `weakly_productive` is distinct from `productive`: -116 reduces BCL11A binding
    # 3.7x where -114 reduces it 55x. The 0.0 value remains provisional pending calibration review.
    penalty_weakly_productive: float = 0.0
    penalty_neutral: float = 0.0
    # ===== Off-target penalty, weighted by mismatch count rather than raw site count. =============
    # The 0-mismatch HBG gene-copy hit is removed earlier. A 4-mismatch site receives less weight
    # than a 1- or 2-mismatch site. Weight values remain subject to protocol review.
    w_offtarget: float = 1.0
    off_w_mm1: float = 8.0    # 1 mismatch receives the largest weight
    off_w_mm2: float = 1.0    # 2 mismatches, the reference unit
    off_w_mm3: float = 0.1
    off_w_mm4: float = 0.002  # 4 mismatches receive the smallest weight
    # No fixed cutoff is applied to candidate-site count or variant allele frequency because the
    # boundary must be set for the chosen assay and its measured background signal.
    # ===== Efficiency term above the threshold, preventing an inactive guide from ranking first. ==
    w_efficiency_reward: float = 0.5  # applied to max(0, efficiency - preclinical_floor)
    # ===== Regulatory-consequence penalty, filled in once AlphaGenome has run. =====
    # Penalizes predicted regulatory severity for bystanders that are not productive. Productive
    # bystanders are the point of the edit and are left alone.
    w_consequence: float = 2.0
    # ===== When a CH-driver off-target should drop a design to the bottom tier. =====
    # The default requires at most three mismatches and an exon or untranslated-region location.
    # Other CH-gene sites remain visible without changing the tier. These thresholds require review.
    ch_gate_mm_threshold: int = 3          # maximum mismatches for the CH-flagged tier
    ch_gate_requires_exon: bool = True     # require an exon or UTR for the CH-flagged tier
    ch_gate_trust_proxy_consequence: bool = False  # allow a stand-in A>G score to affect the tier

    provisional: bool = True          # set False only after documented calibration review

    # Lower tier numbers are reviewed first. A qualifying CH flag overrides other terms.
    tier_preferred: int = 0
    tier_preclinical: int = 1
    tier_below_activity: int = 2
    tier_ch_flagged: int = 3
    # Use this tier when none of the intended edits has a predictor result.
    tier_not_evaluated: int = 4
    tier_labels = {
        0: "PREFERRED",
        1: "PRECLINICAL_ACTIVITY_ONLY",
        2: "BELOW_EVIDENCED_ACTIVITY",
        3: "CH_DRIVER_FLAGGED",
        4: "NOT_EVALUATED",
    }

    def as_dict(self) -> dict:
        # Return parameters and items awaiting review. Numeric tier codes are used only for sorting.
        internal = {"provisional", "tier_preferred", "tier_preclinical", "tier_below_activity",
                    "tier_ch_flagged", "tier_not_evaluated"}
        d = {"status": WEIGHTS_STATUS if self.provisional else "calibrated"}
        d.update({k: v for k, v in asdict(self).items() if k not in internal})
        d["weights_to_review"] = [
            "clinical_floor (a comparison boundary, not a clinical minimum)",
            "penalty_adverse", "penalty_unknown",
            "off_w_mm1..4 (off-target mismatch weighting)",
            "w_efficiency_reward (ranking weight above the activity boundary)",
            "ch_gate_mm_threshold", "ch_gate_requires_exon",
            "ch_gate_trust_proxy_consequence (temporary A>G substitute)",
        ]
        return d


@dataclass
class GuideScore:
    guide_id: str
    gene: str
    editor_class: str
    anchors: list
    # Single-position rate or joint multi-edit rate; None when the required value was not evaluated.
    on_target_efficiency: float | None
    efficacy_band: str
    meets_clinical_floor: bool | None     # None when the efficiency is unknown
    bystander_penalty: float
    bystander_breakdown: dict              # edit rate added up per bystander class
    offtarget_count: int                   # raw count from the tool, paralog included
    offtarget_mm_histogram: dict           # mismatches (>=1) -> number of sites, paralog dropped
    offtarget_cfd_burden: float            # CFD summed over non-paralog sites, a second view
    mit_specificity: str
    ch_gate: bool                          # a CH-driver off-target drops this to the bottom tier
    ch_flag_noted: bool                    # a CH-driver candidate is noted without changing tier
    ch_flag_reason: object
    ch_genes: list
    ch_functional_consequence: object      # filled in by annotate_consequence.py
    consequence_severity: float
    tier: int
    tier_label: str
    safety_score: float
    terms: dict = field(default_factory=dict)
    rank: int = 0

    def as_dict(self) -> dict:
        d = asdict(self)
        # Lead with the tier and the broken-down terms. safety_score orders designs inside a tier;
        # it is not a clinical safety measure, so it comes after the supporting terms.
        primary = ["guide_id", "gene", "editor_class", "anchors", "tier", "tier_label",
                   "efficacy_band", "terms"]
        ordered = {k: d[k] for k in primary if k in d}
        ordered.update({k: v for k, v in d.items() if k not in ordered})
        return ordered


def _to_int(value, default=0) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def _efficacy_band(eff: float, w: RubricWeights) -> str:
    if eff >= w.strong_band_low:
        return "STRONG_CLINICAL_RANGE"
    if eff >= w.clinical_floor:
        return "MEETS_CLINICAL_PRECEDENT"
    if eff >= w.preclinical_floor:
        return "PRECLINICAL_ACTIVITY_ONLY"
    return "BELOW_EVIDENCED_ACTIVITY"


def _allele_key(record: dict) -> tuple[int, object, object]:
    """Return the stable sense-strand key linking a joint rate to intended edits."""
    ref = record.get("promoter_sense_ref")
    alt = record.get("promoter_sense_alt")
    if (not isinstance(ref, str) or ref not in "ACGT"
            or not isinstance(alt, str) or alt not in "ACGT" or ref == alt):
        raise ValueError(
            f"invalid intended allele key at {-abs(record['promoter_offset'])}: {ref!r}>{alt!r}"
        )
    return (
        -abs(record["promoter_offset"]),
        ref,
        alt,
    )


def _on_target_efficiency(on_target: list[dict]) -> tuple[float | None, str]:
    """Select the only efficiency quantity that represents the complete intended design."""
    if len(on_target) == 1:
        return on_target[0]["predicted_rate"], "single_edit_rate"

    # A multi-edit design is one haplotype-level unit. Missing joint data is not permission to use a
    # maximum or minimum of separate position rates: the required joint value remains unevaluated.
    joint_payloads = []
    for record in on_target:
        provenance = record.get("provenance") or {}
        if "joint_on_target_rate" not in provenance:
            continue
        payload = (
            provenance.get("joint_on_target_rate"),
            provenance.get("joint_on_target_confidence"),
            provenance.get("joint_on_target_method"),
            tuple(
                (key.get("promoter_offset"), key.get("promoter_sense_ref"),
                 key.get("promoter_sense_alt"))
                for key in (provenance.get("joint_on_target_allele_keys") or [])
            ),
            provenance.get("joint_on_target_model"),
            provenance.get("predictor_runs"),
            provenance.get("joint_on_target_caveat"),
        )
        joint_payloads.append(payload)

    if not joint_payloads or len(joint_payloads) != len(on_target):
        return None, "joint_on_target_rate_not_evaluated"
    if any(payload != joint_payloads[0] for payload in joint_payloads[1:]):
        raise ValueError("inconsistent joint_on_target_rate provenance across one design")

    joint_rate, confidence, method, allele_keys, model, predictor_runs, caveat = joint_payloads[0]
    expected_keys = {_allele_key(record) for record in on_target}
    if len(expected_keys) != len(on_target):
        raise ValueError("multi-edit design contains a duplicate intended allele key")
    if len(set(allele_keys)) != len(allele_keys) or set(allele_keys) != expected_keys:
        raise ValueError(
            "joint_on_target_rate allele keys do not match the design's intended edits"
        )
    if confidence == "not_evaluated":
        if joint_rate is not None:
            raise ValueError("not_evaluated joint_on_target_rate must be None")
        return None, "joint_on_target_rate_not_evaluated"
    if confidence != "model_predicted":
        raise ValueError(f"unsupported joint_on_target_rate confidence {confidence!r}")
    if method != JOINT_ON_TARGET_METHOD:
        raise ValueError(f"unsupported joint_on_target_rate method {method!r}")
    if not isinstance(model, str) or not model:
        raise ValueError("model-predicted joint_on_target_rate needs model provenance")
    if not isinstance(predictor_runs, int) or isinstance(predictor_runs, bool) or predictor_runs < 1:
        raise ValueError("model-predicted joint_on_target_rate needs a positive predictor run count")
    if caveat != JOINT_ON_TARGET_CAVEAT:
        raise ValueError("joint_on_target_rate is missing its model-predicted-not-measured caveat")
    if joint_rate is None:
        raise ValueError("model_predicted joint_on_target_rate cannot be None")
    try:
        joint_rate = float(joint_rate)
    except (TypeError, ValueError):
        raise ValueError(f"joint_on_target_rate is not numeric: {joint_rate!r}")
    if not 0.0 <= joint_rate <= 1.0:
        raise ValueError(f"joint_on_target_rate {joint_rate} outside [0, 1]")

    marginals = [record["predicted_rate"] for record in on_target]
    if any(rate is None for rate in marginals):
        raise ValueError("model-predicted joint rate cannot accompany an unevaluated intended edit")
    if joint_rate > min(float(rate) for rate in marginals) + 1e-6:
        raise ValueError("joint_on_target_rate exceeds an intended edit's marginal rate")
    return joint_rate, "joint_on_target_rate"


def score_guides(records: list[dict], weights: RubricWeights | None = None) -> list[GuideScore]:
    """Rank the designs that hit a position this tool targets, safest first."""
    w = weights or RubricWeights()
    pen_by_class = {
        ADVERSE: w.penalty_adverse, PRODUCTIVE: w.penalty_productive,
        WEAKLY_PRODUCTIVE: w.penalty_weakly_productive,
        NEUTRAL: w.penalty_neutral, UNKNOWN: w.penalty_unknown,
    }

    by_guide: dict[str, list[dict]] = {}
    for r in records:
        by_guide.setdefault(r["provenance"]["guide_id"], []).append(r)

    scores: list[GuideScore] = []
    for guide_id, recs in by_guide.items():
        on_target = [r for r in recs if r["on_target_or_bystander"] == "on_target"]
        if not on_target:
            continue
        bystanders = [r for r in recs if r["on_target_or_bystander"] == "bystander"]
        prov = recs[0]["provenance"]
        offtarget = prov.get("offtarget", {})
        ch_hits = prov.get("ch_driver_offtargets", [])

        # A single-edit design uses its position rate. A multi-edit design uses the predicted
        # frequency of complete sequences containing every intended edit. Do not substitute the
        # highest or lowest separate position rate when that combined value is missing.
        evaluated = [r["predicted_rate"] for r in on_target if r["predicted_rate"] is not None]
        on_target_efficiency, efficiency_basis = _on_target_efficiency(on_target)
        unevaluated_anchors = len(on_target) - len(evaluated)

        # Record companion edits on the same DNA molecule. Conditional bystander rules use this
        # state. A position outside the predictor request remains unresolved rather than absent.
        edits_on_allele = allele_context(recs)

        # Multiply each bystander edit rate by the weight for its functional class.
        breakdown = {ADVERSE: 0.0, PRODUCTIVE: 0.0, WEAKLY_PRODUCTIVE: 0.0,
                     NEUTRAL: 0.0, UNKNOWN: 0.0}
        bystander_penalty = 0.0
        for r in bystanders:
            rate = r["predicted_rate"] or 0.0
            if rate <= 0:
                continue
            cls, _ = bystander_function(r["gene"], r["promoter_offset"],
                                        r.get("promoter_sense_ref"), r.get("promoter_sense_alt"),
                                        edits_on_allele)
            breakdown[cls] += rate
            bystander_penalty += pen_by_class[cls] * rate
        bystander_penalty *= w.w_bystander

        offtarget_count = _to_int(offtarget.get("offtarget_count"))
        # Weight off-target candidates by mismatch count. The HBG gene-copy hit is already excluded.
        mm_weight = {1: w.off_w_mm1, 2: w.off_w_mm2, 3: w.off_w_mm3, 4: w.off_w_mm4}
        histogram = offtarget.get("mm_histogram", {}) or {}
        offt_tier_contrib = {}
        offt_weighted = 0.0
        for k, n in histogram.items():
            try:
                m, count = int(k), int(n)
            except (TypeError, ValueError):
                raise ValueError(f"malformed off-target histogram for {guide_id}: {k!r} -> {n!r} "
                                 f"(expected int mismatch-count -> int site-count)")
            if m < 0 or count < 0:
                raise ValueError(f"malformed off-target histogram for {guide_id}: negative value {m}:{count}")
            wt = mm_weight.get(m, w.off_w_mm4)  # the search caps at 4; use that weight above the cap
            offt_tier_contrib[str(m)] = round(wt * count, 6)
            offt_weighted += wt * count
        cfd_burden = float(offtarget.get("cfd_weighted_burden") or 0.0)

        # A CH-driver candidate changes the tier only when its mismatch count is at or below the
        # configured limit and it has either an accepted functional class or an exon/UTR coordinate.
        gating_hits, noted_hits = [], []
        for h in ch_hits:
            mm = _to_int(h.get("mismatch_count"), 99)
            cfc = h.get("ch_functional_consequence") or {}
            classification = (cfc.get("classification") or "").lower()
            # A stand-in A>G score does not affect the tier unless configuration permits stand-in scores.
            trust_call = bool(classification) and not (cfc.get("is_proxy") and not w.ch_gate_trust_proxy_consequence)
            if trust_call:
                functionally_adverse = classification == "functionally-consequential"
            else:  # without a trusted consequence class, use the coordinate feature
                functionally_adverse = ((h.get("feature") or "").lower() in ("exon", "utr")
                                        or not w.ch_gate_requires_exon)
            gates = (mm <= w.ch_gate_mm_threshold) and functionally_adverse
            (gating_hits if gates else noted_hits).append(h)
        ch_gate = len(gating_hits) > 0
        ch_flag_noted = (not ch_gate) and len(noted_hits) > 0
        ch_functional_consequence = ch_hits[0].get("ch_functional_consequence") if ch_hits else None
        ch_flag_reason = None
        if ch_flag_noted:
            h = noted_hits[0]
            cls = ((h.get("ch_functional_consequence") or {}).get("classification") or "n/a")
            ch_flag_reason = (f"CH-driver off-target in {h['ch_driver_gene']}, noted but not gating "
                              f"({h.get('mismatch_count')}-mm, {h.get('feature')}, consequence={cls}). "
                              f"Confirm it before use.")
        ch_genes = sorted({h["ch_driver_gene"] for h in ch_hits})

        # Regulatory severity of the bystanders that are not productive.
        consequence_severity = 0.0
        for r in bystanders:
            rate = r["predicted_rate"] or 0.0
            cons = r.get("consequence")
            sev = cons.get("severity") if isinstance(cons, dict) else None
            if sev is None or rate <= 0:
                continue
            cls, _ = bystander_function(r["gene"], r["promoter_offset"],
                                        r.get("promoter_sense_ref"), r.get("promoter_sense_alt"),
                                        edits_on_allele)
            if cls in (PRODUCTIVE, WEAKLY_PRODUCTIVE):
                continue  # the rubric does not penalize a bystander classified as productive
            consequence_severity += sev * rate

        if on_target_efficiency is None:
            band = NOT_EVALUATED_BAND
            meets_clinical_floor = None
            tier = w.tier_ch_flagged if ch_gate else w.tier_not_evaluated
        else:
            band = _efficacy_band(on_target_efficiency, w)
            meets_clinical_floor = on_target_efficiency >= w.clinical_floor
            if ch_gate:
                tier = w.tier_ch_flagged
            elif on_target_efficiency >= w.clinical_floor:
                tier = w.tier_preferred
            elif on_target_efficiency >= w.preclinical_floor:
                tier = w.tier_preclinical
            else:
                tier = w.tier_below_activity

        offt_pen = w.w_offtarget * offt_weighted
        cons_pen = w.w_consequence * consequence_severity
        eff_reward = (0.0 if on_target_efficiency is None
                      else w.w_efficiency_reward * max(0.0, on_target_efficiency - w.preclinical_floor))
        safety_score = -(bystander_penalty + offt_pen + cons_pen) + eff_reward

        scores.append(GuideScore(
            guide_id=guide_id,
            gene=recs[0]["gene"],
            editor_class=recs[0]["editor_class"],
            anchors=sorted({-abs(r["promoter_offset"]) for r in on_target}),
            on_target_efficiency=(None if on_target_efficiency is None
                                  else round(on_target_efficiency, 6)),
            efficacy_band=band,
            meets_clinical_floor=meets_clinical_floor,
            bystander_penalty=round(bystander_penalty, 6),
            bystander_breakdown={k: round(v, 6) for k, v in breakdown.items() if v > 0},
            offtarget_count=offtarget_count,
            offtarget_mm_histogram=histogram,
            offtarget_cfd_burden=round(cfd_burden, 4),
            mit_specificity=offtarget.get("mit_specificity"),
            ch_gate=ch_gate,
            ch_flag_noted=ch_flag_noted,
            ch_flag_reason=ch_flag_reason,
            ch_genes=ch_genes,
            ch_functional_consequence=ch_functional_consequence,
            consequence_severity=consequence_severity,
            tier=tier,
            tier_label=RubricWeights.tier_labels[tier],
            # Rounded so it does not look more precise than it is. The tier and `terms` are what to
            # read; this number only orders designs within a tier.
            safety_score=round(safety_score, 3),
            terms={
                "bystander_penalty": round(bystander_penalty, 6),
                "bystander_breakdown": {k: round(v, 6) for k, v in breakdown.items() if v > 0},
                "offtarget_penalty": round(offt_pen, 6),
                "offtarget_mm_histogram": histogram,           # raw counts, paralog dropped
                "offtarget_penalty_by_mm": offt_tier_contrib,  # what each mismatch tier contributed
                "offtarget_cfd_burden": round(cfd_burden, 4),  # the CFD-weighted view
                "offtarget_paralog_excluded": True,
                "consequence_penalty": round(cons_pen, 6),
                "consequence_severity_nonproductive_bystanders": round(consequence_severity, 6),
                "efficiency_reward": round(eff_reward, 6),
                "efficacy_band": band,
                "efficiency_basis": efficiency_basis,
                "unevaluated_anchors": unevaluated_anchors,
                "ch_gate": ch_gate,
                "ch_flag_noted": ch_flag_noted,
                "ch_flag_reason": ch_flag_reason,
            },
        ))

    # Sort by tier first, then safety score, then efficiency as a tie-break.
    # Within a tier, use the higher rubric score and then higher efficiency. Designs without an
    # efficiency value occupy a separate tier, so a numeric value is not compared with a missing one.
    scores.sort(key=lambda s: (s.tier, -s.safety_score,
                               0.0 if s.on_target_efficiency is None else -s.on_target_efficiency))
    for i, s in enumerate(scores, start=1):
        s.rank = i
    return scores
