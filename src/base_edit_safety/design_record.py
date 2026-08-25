"""Build one review record per candidate guide.

Each record keeps the predicted edit spectrum, off-target candidate summary, HBG paralog checks,
CH-driver overlap, evidence links, and follow-up assays together. It does not collapse those fields
into a composite score.
"""
from __future__ import annotations

from . import annotation, anchors, evidence, hbg_registry, offtarget_kinds, paralog_flags
from .pipeline import environment
from .pipeline.datatypes import JOINT_ON_TARGET_CAVEAT, JOINT_ON_TARGET_METHOD, NOT_EVALUATED

SCHEMA_VERSION = "3"
_MECHANISM_BY_OFFSET = {p.promoter_offset: p.mechanism for p in anchors.IN_SCOPE}
_REL_ORDER = (hbg_registry.NO_MATCH, hbg_registry.NEARBY,
              hbg_registry.SAME_MOTIF, hbg_registry.EXACT_EDIT)


def _edit_str(record: dict) -> dict:
    return {
        "genomic_plus": f"{record['genomic_plus_ref']}>{record['genomic_plus_alt']}",
        "promoter_sense": f"{record['promoter_sense_ref']}>{record['promoter_sense_alt']}",
    }


def _allele_key(record: dict) -> dict:
    """Canonical sense-strand identity for one intended or bystander edit."""
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

    confidence = payload["confidence"]
    rate = payload["rate"]
    if confidence == NOT_EVALUATED:
        if rate is not None:
            raise ValueError("not_evaluated joint_on_target_rate must be None")
        return not_evaluated()
    if confidence != "model_predicted":
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
        },
    }


def _strongest(values, order) -> str:
    present = [v for v in values if v in order]
    return max(present, key=order.index) if present else order[0]


def _validation_checklist(bystanders, off_triggers, flags, ch_annot) -> list[dict]:
    """The lab tests this candidate needs before anyone acts on it."""
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


def _limitations() -> list[str]:
    return [
        "BE-DICT estimates the on-target and bystander rates; confirm them in the target cell type.",
        "The CRISPOR search covers mismatch-only sites in the reference genome. Bulges and "
        "population variants need separate searches, while editor-context and RNA effects remain "
        "separate fields.",
        "The CH annotation reports site overlap and gene context without estimating patient "
        "outcomes.",
        "The current target scope is the HBG1/HBG2 proximal promoter. Clinical-program matches "
        "compare motif and mechanism when exact positions are unavailable.",
    ]


def build_design_records(records: list[dict], empirical: dict | None = None) -> dict:
    """Build the design records from annotated enumeration records.

    records:   enumeration EditRecords as dicts, on-target and bystander, with provenance.
    empirical: optional {guide_id: [measurement, ...]}. Attaching real measurements moves a
               candidate to empirically_measured. Nothing else can: not published work, not a model.
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
                            "promoter_offset": a["promoter_offset"], "relationship": m["relationship"],
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
                           "promoter_offset": a["promoter_offset"],
                           "predicted_rate": a.get("predicted_rate"),
                           "model": a.get("model_name"),
                           "confidence": a.get("confidence")} for a in on_target],
            "bystanders": [{"allele_key": _allele_key(b),
                            "promoter_offset": b["promoter_offset"],
                            "predicted_rate": b.get("predicted_rate"),
                            "genomic_plus": _edit_str(b)["genomic_plus"]} for b in bystanders],
        }
        predicted_edit_spectrum.update(_joint_edit_summary(on_target))

        out.append({
            "record_type": "hbg_design_record",
            "schema_version": SCHEMA_VERSION,
            "candidate": {
                "guide_id": guide_id, "gene": base["gene"], "editor_class": base["editor_class"],
                "editor_name": prov.get("editor"),
                "target": {"locus": "HBG1/HBG2 proximal promoter", "chrom": base["chrom"],
                           "intended_edits": [{"allele_key": _allele_key(a),
                                               "promoter_offset": a["promoter_offset"],
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
            "off_target": {"kinds": off_kinds, "assay_triggers": off_triggers},
            "hbg_paralog_flags": flags,
            "ch_annotation": {"driver_offtargets": ch_annot,
                              "note": "for review and monitoring; no estimate is produced"},
            "empirical_support": empirical_support,
            "required_validation": _validation_checklist(bystanders, off_triggers, flags, ch_annot),
            "provenance": {"guide_id": guide_id, "model_name": base.get("model_name"),
                           "predictor_runs": prov.get("predictor_runs"),
                           "protospacer": prov.get("protospacer"), "pam": prov.get("pam"),
                           "strand": prov.get("strand")},
        })

    return {
        "record_type": "hbg_design_record_set",
        "schema_version": SCHEMA_VERSION,
        "count": len(out),
        "records": out,
        "limitations": _limitations(),
        "reproducibility": {
            "generated_at": environment.run_stamp(),
            "repo": environment.repo_state(),
            "genome_build": anchors.GENOME_BUILD,
            "models": sorted({r.get("model_name") for r in records if r.get("model_name")}),
        },
    }
