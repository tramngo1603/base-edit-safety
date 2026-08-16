"""Build the primary review record for each candidate guide.

A record contains the available candidate information: the predicted on-target and bystander edits, the
off-target results separated by category, HBG gene-copy flags, CH-driver overlap, comparison with the
HBG registry, evidence level, and required laboratory tests. The score is included only as a
tiebreaker after the evidence and required tests.

A design that installs more than one edit, such as the -123/-124 pair, keeps every intended edit,
and the registry comparison runs over all of them rather than only the highest-rate one.
"""
from __future__ import annotations

from . import (annotation, anchors, evidence, hbg_registry, offtarget_kinds,
               offtarget_reference, paralog_flags)
from .pipeline import environment
from .pipeline.datatypes import JOINT_ON_TARGET_CAVEAT, JOINT_ON_TARGET_METHOD, NOT_EVALUATED
from .scoring.rubric import RubricWeights, score_guides

SCHEMA_VERSION = "2"
_MECHANISM_BY_OFFSET = {p.promoter_offset: p.mechanism for p in anchors.IN_SCOPE}
_REL_ORDER = (hbg_registry.NO_MATCH, hbg_registry.NEARBY,
              hbg_registry.SAME_MOTIF, hbg_registry.EXACT_EDIT)


def _edit_str(record: dict) -> dict:
    return {
        "genomic_plus": f"{record['genomic_plus_ref']}>{record['genomic_plus_alt']}",
        "promoter_sense": f"{record['promoter_sense_ref']}>{record['promoter_sense_alt']}",
    }


def _allele_key(record: dict) -> dict:
    """Return the stable sense-strand identity for one intended or bystander edit."""
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


def _strongest(values, order) -> str:
    present = [v for v in values if v in order]
    return max(present, key=order.index) if present else order[0]


def _validation_checklist(bystanders, off_triggers, flags, ch_annot) -> list[dict]:
    """Return the laboratory tests required before experimental follow-up."""
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


def _limitations(multi_edit: bool = False) -> list[str]:
    limitations = [
        "On-target and bystander rates come from BE-DICT, not from measurements in HSCs.",
        "Off-target candidates come from a Cas9-style search of the reference genome. It finds "
        "mismatch-only sites: no bulged sites and no population-variant-created sites. Measured "
        "assays (GUIDE-seq / CHANGE-seq-BE) detect bulged sites; population-variant off-targets "
        "need a variant-aware module that does not exist in v1. Off-targets specific to the base "
        "editor, guide-independent deaminase activity, and RNA editing are not checked here.",
        "The CH layer records which driver gene a site falls in, plus coordinate and mechanism "
        "context. It produces no number about clonal expansion, leukemia, or long-term outcome.",
        "Matches to clinical programs are motif and mechanism analogs. Exact clinical positions are "
        "not always public. Published work is analog support, not a measurement from this run.",
        "Scope is the HBG1/HBG2 proximal promoter only.",
    ]
    if multi_edit:
        limitations.insert(1, JOINT_ON_TARGET_CAVEAT)
    return limitations


def build_design_records(records: list[dict], weights: RubricWeights | None = None,
                         empirical: dict | None = None) -> dict:
    """Build the design records from annotated enumeration records.

    `records` contains on-target and bystander `EditRecord` dictionaries with source details.
    `empirical` optionally maps a guide ID to measurements produced for this run. A non-empty entry
    assigns `empirically_measured`; published studies and predictions do not assign that level.
    """
    weights = weights or RubricWeights()
    empirical = empirical or {}
    scores = {s.guide_id: s for s in score_guides(records, weights)}

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
        score = scores.get(guide_id)
        base = on_target[0]

        # Compare every intended edit so a joint target such as -123/-124 remains one design.
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
        measured_in = "this_run" if tier == evidence.EMPIRICALLY_MEASURED else "literature"

        # If a published study evaluated the same guide, compare its candidate-site list with the
        # current CRISPOR list. Both lists come from search methods rather than measurements. A site
        # absent from the CRISPOR list identifies a search-coverage gap.
        published_panel = None
        panel = offtarget_reference.panel_for(prov.get("protospacer") or "")
        if panel is not None:
            published_panel = offtarget_reference.compare_to_panel(
                prov.get("offtarget_sites_low_mismatch", []), panel)

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
            "model_disagreement": {"status": "single_predictor",
                                   "note": "BE-DICT is the only configured outcome predictor; no cross-predictor comparison was run"},
        }
        if len(on_target) > 1:
            joint_evaluated = (
                score is not None
                and score.terms.get("efficiency_basis") == "joint_on_target_rate"
            )
            joint_prov = next(
                (record.get("provenance") or {} for record in on_target
                 if "joint_on_target_rate" in (record.get("provenance") or {})),
                prov,
            )
            joint_rate = score.on_target_efficiency if joint_evaluated else None
            predicted_edit_spectrum.update({
                "joint_on_target_rate": joint_rate,
                "joint_on_target_rate_provenance": {
                    "allele_keys": joint_prov.get("joint_on_target_allele_keys")
                                   or [_allele_key(a) for a in on_target],
                    "method": joint_prov.get("joint_on_target_method") or JOINT_ON_TARGET_METHOD,
                    "model": joint_prov.get("joint_on_target_model") or base.get("model_name"),
                    "predictor_runs": joint_prov.get("predictor_runs"),
                    "confidence": (joint_prov.get("joint_on_target_confidence")
                                   if joint_evaluated else NOT_EVALUATED),
                    "caveat": joint_prov.get("joint_on_target_caveat") or JOINT_ON_TARGET_CAVEAT,
                },
            })

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
                "measured_in": measured_in,
                "by_edit": by_edit,
                "note": ("Edit relationship and evidence level are reported separately. Published "
                         "work, including an exact HSPC edit, can assign analog_supported; "
                         "empirically_measured requires a measurement attached to this run."),
            },
            "predicted_edit_spectrum": predicted_edit_spectrum,
            "off_target": {"kinds": off_kinds, "assay_triggers": off_triggers,
                           "published_panel": published_panel},
            "hbg_paralog_flags": flags,
            "ch_annotation": {"driver_offtargets": ch_annot,
                              "note": "for review and monitoring; no estimate is produced"},
            "empirical_support": empirical_support,
            "required_validation": _validation_checklist(bystanders, off_triggers, flags, ch_annot),
            "scoring": {
                "note": "Tiebreaker only; review the evidence level and required assays first.",
                **(score.as_dict() if score else {}),
            },
            "provenance": {"guide_id": guide_id, "model_name": base.get("model_name"),
                           "predictor_runs": prov.get("predictor_runs"),
                           "protospacer": prov.get("protospacer"), "pam": prov.get("pam"),
                           "strand": prov.get("strand")},
            "limitations": _limitations(multi_edit=len(on_target) > 1),
        })

    return {
        "record_type": "hbg_design_record_set",
        "schema_version": SCHEMA_VERSION,
        "count": len(out),
        "records": out,
        "reproducibility": {
            "generated_at": environment.run_stamp(),
            "repo": environment.repo_state(),
            "genome_build": anchors.GENOME_BUILD,
            "models": sorted({r.get("model_name") for r in records if r.get("model_name")}),
            "config": {"rubric_weights": weights.as_dict()},
        },
    }
