"""Build a review record for each candidate guide.

Keep predictions, candidate sites, HBG copy checks, gene overlaps, evidence sources and required
assays together. These fields remain separate; no combined safety score is calculated.
"""
from __future__ import annotations

from . import annotation, anchors, evidence, hbg_registry, offtarget_kinds, paralog_flags
from .pipeline import environment
from .pipeline.datatypes import (
    JOINT_ON_TARGET_CAVEAT,
    JOINT_ON_TARGET_METHOD,
    MODEL_PREDICTED,
    MODEL_PREDICTED_FALLBACK,
    NOT_EVALUATED,
)

SCHEMA_VERSION = "7"
_MECHANISM_BY_OFFSET = {p.promoter_offset: p.mechanism for p in anchors.IN_SCOPE}
_REL_ORDER = (hbg_registry.NO_MATCH, hbg_registry.NEARBY,
              hbg_registry.SAME_MOTIF, hbg_registry.EXACT_EDIT)


def _edit_str(record: dict) -> dict:
    return {
        "genomic_plus": f"{record['genomic_plus_ref']}>{record['genomic_plus_alt']}",
        "promoter_sense": f"{record['promoter_sense_ref']}>{record['promoter_sense_alt']}",
    }


def _allele_key(record: dict) -> dict:
    """Identify an intended or bystander change by promoter-strand position, REF and ALT."""
    ref = record.get("promoter_sense_ref")
    alt = record.get("promoter_sense_alt")
    if (not isinstance(ref, str) or ref not in "ACGT"
            or not isinstance(alt, str) or alt not in "ACGT" or ref == alt):
        raise ValueError(
            f"invalid allele key at {-abs(record['promoter_offset'])}: {ref!r}>{alt!r}"
        )
    return {
        "promoter_offset": -abs(record["promoter_offset"]),
        "promoter_sense_ref": ref,
        "promoter_sense_alt": alt,
    }


def _promoter_location(record: dict) -> dict:
    """Public promoter location with one unambiguous signed convention.

    Negative offsets are upstream of the transcription start site. The positive magnitude is
    exposed under a different field name so it cannot be mistaken for the signed coordinate.
    """
    distance = abs(record["promoter_offset"])
    return {"promoter_offset": -distance, "upstream_distance_bp": distance}


def _joint_edit_summary(on_target: list[dict]) -> dict:
    """Return the validated joint rate for a multi-edit design, if the predictor supplied one."""
    if len(on_target) < 2:
        return {}

    expected_keys = {
        (key["promoter_offset"], key["promoter_sense_ref"], key["promoter_sense_alt"])
        for key in (_allele_key(record) for record in on_target)
    }
    fallback_provenance = on_target[0].get("provenance") or {}

    def not_evaluated() -> dict:
        return {
            "joint_on_target_rate": None,
            "joint_on_target_rate_provenance": {
                "allele_keys": [_allele_key(record) for record in on_target],
                "method": JOINT_ON_TARGET_METHOD,
                "model": on_target[0].get("model_name"),
                "predictor_runs": fallback_provenance.get("predictor_runs"),
                "confidence": NOT_EVALUATED,
                "caveat": JOINT_ON_TARGET_CAVEAT,
                "rate_semantics": fallback_provenance.get("rate_semantics"),
            },
        }

    payloads = []
    for record in on_target:
        provenance = record.get("provenance") or {}
        if "joint_on_target_rate" not in provenance:
            continue
        allele_keys = tuple(
            (key.get("promoter_offset"), key.get("promoter_sense_ref"),
             key.get("promoter_sense_alt"))
            for key in (provenance.get("joint_on_target_allele_keys") or [])
        )
        payloads.append({
            "rate": provenance.get("joint_on_target_rate"),
            "confidence": provenance.get("joint_on_target_confidence"),
            "method": provenance.get("joint_on_target_method"),
            "allele_keys": allele_keys,
            "model": provenance.get("joint_on_target_model"),
            "predictor_runs": provenance.get("predictor_runs"),
            "caveat": provenance.get("joint_on_target_caveat"),
            "rate_semantics": provenance.get("joint_on_target_rate_semantics"),
        })

    if len(payloads) != len(on_target):
        return not_evaluated()
    if any(payload != payloads[0] for payload in payloads[1:]):
        raise ValueError("inconsistent joint_on_target_rate provenance across one design")

    payload = payloads[0]
    allele_keys = payload["allele_keys"]
    if len(set(allele_keys)) != len(allele_keys) or set(allele_keys) != expected_keys:
        raise ValueError("joint_on_target_rate allele keys do not match the intended edits")
    if payload["method"] != JOINT_ON_TARGET_METHOD:
        raise ValueError(f"unsupported joint_on_target_rate method {payload['method']!r}")
    if payload["caveat"] != JOINT_ON_TARGET_CAVEAT:
        raise ValueError("joint_on_target_rate is missing its required provenance note")
    if payload["rate_semantics"] != "model_relative":
        raise ValueError("joint_on_target_rate must be labelled model_relative")

    confidence = payload["confidence"]
    rate = payload["rate"]
    if confidence == NOT_EVALUATED:
        if rate is not None:
            raise ValueError("not_evaluated joint_on_target_rate must be None")
        return not_evaluated()
    if confidence not in (MODEL_PREDICTED, MODEL_PREDICTED_FALLBACK):
        raise ValueError(f"unsupported joint_on_target_rate confidence {confidence!r}")
    if not isinstance(payload["model"], str) or not payload["model"]:
        raise ValueError("model-predicted joint_on_target_rate needs model provenance")
    runs = payload["predictor_runs"]
    if not isinstance(runs, int) or isinstance(runs, bool) or runs < 1:
        raise ValueError("model-predicted joint_on_target_rate needs a positive predictor run count")
    try:
        rate = float(rate)
    except (TypeError, ValueError):
        raise ValueError(f"joint_on_target_rate is not numeric: {rate!r}")
    if not 0.0 <= rate <= 1.0:
        raise ValueError(f"joint_on_target_rate {rate} outside [0, 1]")
    marginals = [record.get("predicted_rate") for record in on_target]
    if any(value is None for value in marginals):
        raise ValueError("model-predicted joint rate cannot accompany an unevaluated intended edit")
    if rate > min(float(value) for value in marginals) + 1e-6:
        raise ValueError("joint_on_target_rate exceeds an intended edit's marginal rate")

    return {
        "joint_on_target_rate": round(rate, 6),
        "joint_on_target_rate_provenance": {
            "allele_keys": [
                {"promoter_offset": offset, "promoter_sense_ref": ref,
                 "promoter_sense_alt": alt}
                for offset, ref, alt in allele_keys
            ],
            "method": payload["method"],
            "model": payload["model"],
            "predictor_runs": runs,
            "confidence": confidence,
            "caveat": payload["caveat"],
            "rate_semantics": payload["rate_semantics"],
        },
    }


def _strongest(values, order) -> str:
    present = [v for v in values if v in order]
    return max(present, key=order.index) if present else order[0]


def _validation_checklist(bystanders, off_triggers, flags, ch_annot) -> list[dict]:
    """List the experiments needed to assess this design."""
    items = [{"validates": "on_target_edit_rate",
              "assay": "amplicon NGS in the target cell type (CD34+ HSPCs)", "status": "required"}]
    if bystanders:
        items.append({"validates": "bystander_spectrum",
                      "assay": "amplicon NGS covering the whole editing window",
                      "status": "required"})
    items.extend(off_triggers)
    items.extend(flags.get("required_confirmation", []))
    if ch_annot:
        items.append({"validates": "ch_driver_offtarget",
                      "assay": "targeted off-target sequencing plus a functional check at the "
                               "CH-driver site", "status": "required"})
    return items


def _limitations(records: list[dict]) -> list[str]:
    limitations = [
        "BE-DICT and any executed CRISPRon-ABE fallback supply model-relative on-target and "
        "bystander rates, not calibrated probabilities; confirm them in the target cell type.",
        "The CRISPOR search covers mismatch-only sites in the reference genome. Bulges and "
        "population variants need separate searches, while editor-context and RNA effects remain "
        "separate fields.",
        "The CH annotation reports site overlap and gene context without estimating patient "
        "outcomes.",
        "The current target scope is the HBG1/HBG2 proximal promoter. Clinical-program matches "
        "compare motif and mechanism when exact positions are unavailable.",
    ]
    pam_classes = {
        (record.get("provenance") or {}).get("pam_class", "NGG") for record in records
    }
    if any(pam_class != "NGG" for pam_class in pam_classes):
        limitations[0] += (
            " For non-NGG designs, outcome and bystander rates are not evaluated because the "
            "available BE-DICT models were trained for NGG designs."
        )
    fallback_records = [
        record for record in records
        if (record.get("provenance") or {}).get("analysis_role") == "sensitivity_fallback"
    ]
    if any(
        (record.get("provenance") or {}).get("fallback_status") == "executed"
        for record in fallback_records
    ):
        limitations.append(
            "Structurally unstable BE-DICT ensembles retain their exclusion ledger. Executed "
            "CRISPRon-ABE fallback rates are independent model predictions; rejected BE-DICT "
            "scores are neither revived nor blended."
        )
    elif fallback_records:
        limitations.append(
            "Structurally unstable BE-DICT ensembles carry no BE-DICT rate; affected guides are "
            "marked for CRISPRon sensitivity/fallback analysis."
        )
    return limitations


def build_design_records(records: list[dict], empirical: dict | None = None) -> dict:
    """Build design records from annotated guide-search results.

    records contains intended and bystander EditRecords as dictionaries with source details.
    empirical optionally maps guide IDs to measurements. Only attached measurements can set
    empirically_measured. This design label does not establish per-site or per-allele coverage.
    """
    empirical = empirical or {}

    by_guide: dict[str, list[dict]] = {}
    for r in records:
        by_guide.setdefault(r["provenance"]["guide_id"], []).append(r)

    out = []
    for guide_id, recs in by_guide.items():
        on_target = sorted((r for r in recs if r["on_target_or_bystander"] == "on_target"),
                           key=lambda r: abs(r["promoter_offset"]))
        if not on_target:
            continue  # a guide only gets a record if it hits a position this tool targets
        bystanders = [r for r in recs if r["on_target_or_bystander"] == "bystander"]
        prov = recs[0]["provenance"]
        ch_hits = prov.get("ch_driver_offtargets", [])
        base = on_target[0]

        # Compare every intended edit, not just the first, so a pair like -123/-124 stays intact.
        by_edit = []
        for a in on_target:
            mech = _MECHANISM_BY_OFFSET.get(abs(a["promoter_offset"]))
            m = hbg_registry.relationship(a["gene"], a["promoter_offset"], a["editor_class"], mech)
            by_edit.append({"allele_key": _allele_key(a),
                            **_promoter_location(a), "relationship": m["relationship"],
                            "evidence_strength": m["evidence_strength"],
                            "matched_entries": m["matched_entries"]})
        best_rel = _strongest([m["relationship"] for m in by_edit], _REL_ORDER)
        best_strength = _strongest([m["evidence_strength"] for m in by_edit],
                                   list(hbg_registry.STRENGTH_ORDER))
        empirical_support = empirical.get(guide_id, [])
        tier = evidence.evidence_tier(best_rel, empirical_support)

        off_kinds = offtarget_kinds.classify_offtargets(prov.get("offtarget", {}))
        off_triggers = offtarget_kinds.assay_triggers(off_kinds)
        flags = paralog_flags.paralog_flags(base)
        ch_annot = [
            {**{k: h.get(k) for k in ("ch_driver_gene", "chrom", "start", "mismatch_count", "feature")},
             "consequence_class": annotation.ch_consequence_class(h.get("ch_driver_gene"), h.get("feature"))}
            for h in ch_hits
        ]

        predicted_edit_spectrum = {
            "on_target": [{"allele_key": _allele_key(a),
                           **_promoter_location(a),
                           "predicted_rate": a.get("predicted_rate"),
                           "model": a.get("model_name"),
                           "confidence": a.get("confidence"),
                           "rate_semantics": (a.get("provenance") or {}).get("rate_semantics")}
                          for a in on_target],
            "bystanders": [{"allele_key": _allele_key(b),
                            **_promoter_location(b),
                            "predicted_rate": b.get("predicted_rate"),
                            "confidence": b.get("confidence"),
                            "rate_semantics": (b.get("provenance") or {}).get("rate_semantics"),
                            "genomic_plus": _edit_str(b)["genomic_plus"]} for b in bystanders],
        }
        predicted_edit_spectrum.update(_joint_edit_summary(on_target))

        design_provenance = {
            "guide_id": guide_id,
            "model_name": base.get("model_name"),
            "predictor_runs": prov.get("predictor_runs"),
            "protospacer": prov.get("protospacer"),
            "pam": prov.get("pam"),
            "pam_class": prov.get("pam_class", "NGG"),
            "strand": prov.get("strand"),
            "normalization_applied": prov.get("normalization_applied"),
            "structural_stability": prov.get("structural_stability"),
            "upstream_mass_anomaly": prov.get("upstream_mass_anomaly"),
            "raw_mass": prov.get("raw_mass"),
            "raw_mass_per_member": prov.get("raw_mass_per_member"),
            "rate_semantics": prov.get("rate_semantics"),
        }
        for key in (
            "predictor",
            "predictor_role",
            "model_relative",
            "prediction_tier",
            "rate_transform",
            "target_context_30mer",
            "analysis_role",
            "fallback_predictor",
            "fallback_status",
            "fallback_executed_date",
            "fallback_model_version",
            "fallback_model_member_count",
            "fallback_model_weights",
            "fallback_predicted_editing_efficiency_percent",
            "fallback_predicted_editing_efficiency_rate",
            "fallback_outcome_spectrum",
            "fallback_context_deduplicated",
            "fallback_shared_input_guide_ids",
            "primary_predictor",
            "primary_predictor_status",
            "primary_predictor_exclusion_link",
            "primary_predictor_exclusion_reason",
            "rejected_bedict_scores_used",
            "fallback_blended_with_primary",
        ):
            if prov.get(key) is not None:
                design_provenance[key] = prov[key]
        if prov.get("reason"):
            design_provenance["reason"] = prov["reason"]

        design = {
            "record_type": "hbg_design_record",
            "schema_version": SCHEMA_VERSION,
            "candidate": {
                "guide_id": guide_id, "gene": base["gene"], "editor_class": base["editor_class"],
                "editor_name": prov.get("editor"),
                "target": {"locus": "HBG1/HBG2 proximal promoter", "chrom": base["chrom"],
                           "intended_edits": [{"allele_key": _allele_key(a),
                                               **_promoter_location(a),
                                               "pos_1based": a["pos_1based"],
                                               "mechanism": _MECHANISM_BY_OFFSET.get(abs(a["promoter_offset"])),
                                               "edit": _edit_str(a)} for a in on_target]},
            },
            "evidence_tier": tier,
            "registry_match": {
                "relationship": best_rel,
                "evidence_strength": best_strength,
                "by_edit": by_edit,
                "note": ("Edit relationship and evidence strength are reported separately. "
                         "The empirically_measured tier is reserved for measurements attached "
                         "to this record."),
            },
            "predicted_edit_spectrum": predicted_edit_spectrum,
            "off_target": {
                "kinds": off_kinds,
                "candidate_sites": prov.get("offtarget_sites_low_mismatch", []),
                "assay_triggers": off_triggers,
            },
            "hbg_paralog_flags": flags,
            "ch_annotation": {"driver_offtargets": ch_annot,
                              "note": "for review and monitoring; no estimate is produced"},
            "empirical_support": empirical_support,
            "required_validation": _validation_checklist(bystanders, off_triggers, flags, ch_annot),
            "provenance": design_provenance,
        }
        if prov.get("prediction_tier") is not None:
            design["prediction_tier"] = prov["prediction_tier"]
        out.append(design)

    return {
        "record_type": "hbg_design_record_set",
        "schema_version": SCHEMA_VERSION,
        "count": len(out),
        "records": out,
        "limitations": _limitations(records),
        "reproducibility": {
            "generated_at": environment.run_stamp(),
            "repo": environment.repo_state(),
            "genome_build": anchors.GENOME_BUILD,
            "models": sorted({r.get("model_name") for r in records if r.get("model_name")}),
        },
    }
