"""Search the HBG promoter windows and request supported edit predictions.

Run the external guide search, select compatible editors and make per-position records from
supported predictions. Mark intended bases, retain candidate sites, check HBG copies and gene
overlap, and write results under data/outputs. Leave unsupported predictions unassessed.
"""
from __future__ import annotations

from collections import defaultdict
import sys

from .. import anchors
from ..adapters.bedict import (
    NORMALIZATION_APPLIED,
    RATE_SEMANTICS,
    RAW_MASS_MAX,
    RAW_MASS_MIN,
    SUPERSEDED_RAW_MASS_MAX,
    SUPERSEDED_RAW_MASS_MIN,
    BedictBystanderPredictor,
)
from ..adapters.crispor import CrisporGuideDesigner
from ..adapters.crispron_abe import (
    FALLBACK_EXECUTED_DATE,
    FALLBACK_STATUS_EXECUTED,
    PREDICTOR_ROLE as CRISPRON_PREDICTOR_ROLE,
    RATE_SEMANTICS as CRISPRON_RATE_SEMANTICS,
    CrispronAbeFallbackPredictor,
)
from ..annotation import annotate_records as annotate_ch, ch_driver_offtargets
from ..config import crispron_abe_configured
from ..io import outputs
from ..io import reference
from . import environment
from .datatypes import (
    NON_NGG_MODEL_SCOPE_REASON,
    NOT_EVALUATED,
    EditRecord,
    spectrum_to_records,
    structural_records_without_rates,
)

# Our benchmarks include measured positives at four mismatches, while GUIDE-seq and CHANGE-seq can
# report sites through roughly six mismatches. Keep every 1-3-mm site and the stronger 4-mm sites;
# 5-mm and gapped-site enumeration remain documented out of scope. The full mismatch histogram is
# independent of this coordinate-retention policy and remains unchanged.
PERSIST_ALL_SITES_UP_TO_MISMATCHES = 3
PERSIST_4MM_CFD_FLOOR = 0.10

# The promoter windows, -210 to -80 from the cap site, on GRCh38. These coordinates are the point of
# this tool, not a placeholder to generalize later. See docs/coordinates.md.
WINDOWS = {
    "HBG1": ("chr11", 5_249_937, 5_250_067),
    "HBG2": ("chr11", 5_254_861, 5_254_991),
}
EDITOR_NAME = {"ABE": "ABE8e", "CBE": "BE4max"}


def _anchor_coords_by_gene_editor() -> dict[tuple[str, str], set[int]]:
    """The genomic coordinates this tool targets, keyed by (gene, editor class)."""
    coords: dict[tuple[str, str], set[int]] = defaultdict(set)
    for pos in anchors.IN_SCOPE:
        coords[("HBG1", pos.editor_class)].add(pos.hbg1_pos)
        coords[("HBG2", pos.editor_class)].add(pos.hbg2_pos)
    return coords


def _low_mismatch_sites(sites: list[dict]) -> list[dict]:
    """Sites retained for review, trimmed to the fields a coordinate cross-check needs."""
    kept = []
    for s in sites:
        try:
            mm = int(s.get("mismatchCount"))
        except (TypeError, ValueError):
            continue
        raw_cfd = s.get("cfdOfftargetScore")
        try:
            cfd = float(raw_cfd)
        except (TypeError, ValueError):
            cfd = None
        if mm > PERSIST_ALL_SITES_UP_TO_MISMATCHES:
            if mm != 4 or cfd is None or cfd < PERSIST_4MM_CFD_FLOOR:
                continue
        kept.append({
            "chrom": s.get("chrom"),
            "start": s.get("start"),
            "strand": s.get("strand"),
            "mismatch_count": mm,
            "cfd_offtarget_score": raw_cfd,
            "locus_desc": s.get("locusDesc"),
            "identity_status": (
                "requires_intended_locus_coordinate_match" if mm == 0 else "candidate"
            ),
        })

    def sort_key(site: dict) -> tuple:
        try:
            cfd_desc = -float(site.get("cfd_offtarget_score"))
        except (TypeError, ValueError):
            cfd_desc = float("inf")
        return (site["mismatch_count"], cfd_desc, str(site["chrom"]), str(site["start"]))

    return sorted(kept, key=sort_key)


def run_enumeration(write: bool = True) -> list[EditRecord]:
    """Search both promoter windows and return annotated EditRecords."""
    anchors.check_paralog_spacing()  # stop here rather than downstream if a coordinate has shifted
    designer = CrisporGuideDesigner()
    predictor = BedictBystanderPredictor() if designer.pam == "NGG" else None
    anchor_coords = _anchor_coords_by_gene_editor()

    # Find candidates per gene, and note which protospacers show up in both promoters.
    candidates_by_gene = {}
    protospacers_by_gene: dict[str, set[str]] = {}
    ambiguous: list[str] = []   # guides that mapped to more than one spot: dropped, but reported
    included_guide_count = 0
    prediction_exclusions: list[dict] = []
    fallback_routes: list[dict] = []
    fallback_route_by_guide: dict[str, dict] = {}
    fallback_predictor: CrispronAbeFallbackPredictor | None = None
    fallback_execution_count = 0
    for gene, (chrom, start, end) in WINDOWS.items():
        cands = designer.design(chrom, start, end, gene)
        ambiguous.extend(getattr(designer, "last_ambiguous", []))
        candidates_by_gene[gene] = cands
        protospacers_by_gene[gene] = {c.protospacer for c in cands}

    all_records: list[EditRecord] = []
    for gene, cands in candidates_by_gene.items():
        other_gene = "HBG2" if gene == "HBG1" else "HBG1"
        # Pair every candidate with the editors that have something to edit in it.
        editor_guides = []
        candidate_by_guide_id = {}
        for cand in cands:
            for editor_class in ("ABE", "CBE"):
                if cand.has_editable_base(editor_class):
                    guide = cand.for_editor(editor_class, EDITOR_NAME[editor_class])
                    editor_guides.append(guide)
                    candidate_by_guide_id[guide.guide_id] = cand

        if predictor is not None:
            batch = predictor.predict_many(editor_guides)
            included_guide_count += len(batch.spectra)
            for exclusion in batch.exclusions:
                prediction_exclusions.append(exclusion)
                if exclusion.get("fallback"):
                    route = {
                        "guide_id": exclusion["guide_id"],
                        **exclusion["fallback"],
                        "source_exclusion_status": exclusion["status"],
                    }
                    fallback_routes.append(route)
                    fallback_route_by_guide[exclusion["guide_id"]] = route
                print(f"BE-DICT EXCLUDED: {exclusion['reason']}", file=sys.stderr)

            exclusions_by_guide = {
                exclusion["guide_id"]: exclusion for exclusion in batch.exclusions
            }
            for guide in editor_guides:
                if guide.guide_id not in batch.spectra:
                    exclusion = exclusions_by_guide.get(guide.guide_id)
                    if not exclusion or not exclusion.get("fallback"):
                        continue
                    cand = candidate_by_guide_id[guide.guide_id]
                    on_coords = anchor_coords.get((gene, guide.editor_class), set())
                    fallback = exclusion["fallback"]
                    execute_fallback = (
                        fallback.get("predictor") == "CRISPRon-ABE"
                        and crispron_abe_configured()
                    )
                    if execute_fallback:
                        if fallback_predictor is None:
                            fallback_predictor = CrispronAbeFallbackPredictor()
                        fallback_spectrum = fallback_predictor.predict_many([guide])[guide.guide_id]
                        records = spectrum_to_records(
                            fallback_spectrum, on_target_coords=on_coords
                        )
                        included_guide_count += 1
                        fallback_execution_count += 1
                        route = fallback_route_by_guide[guide.guide_id]
                        route.update({
                            "original_route_status": fallback["status"],
                            "status": FALLBACK_STATUS_EXECUTED,
                            "executed_date": FALLBACK_EXECUTED_DATE,
                            "predictor_role": CRISPRON_PREDICTOR_ROLE,
                            "model_name": fallback_spectrum.model_name,
                        })
                        fallback_provenance = fallback_spectrum.provenance
                        for record in records:
                            record.provenance = {
                                **record.provenance,
                                "predictor": fallback_provenance["predictor"],
                                "predictor_role": fallback_provenance["predictor_role"],
                                "model_relative": True,
                                "prediction_tier": "predicted_only",
                                "rate_semantics": CRISPRON_RATE_SEMANTICS,
                                "rate_transform": fallback_provenance["rate_transform"],
                                "target_context_30mer": fallback_provenance[
                                    "target_context_30mer"
                                ],
                                "fallback_model_version": fallback_provenance["model_version"],
                                "fallback_model_member_count": fallback_provenance[
                                    "model_member_count"
                                ],
                                "fallback_model_weights": fallback_provenance["model_weights"],
                                "fallback_predicted_editing_efficiency_percent": (
                                    fallback_provenance[
                                        "predicted_editing_efficiency_percent"
                                    ]
                                ),
                                "fallback_predicted_editing_efficiency_rate": (
                                    fallback_provenance[
                                        "predicted_editing_efficiency_rate"
                                    ]
                                ),
                                "fallback_outcome_spectrum": fallback_provenance[
                                    "outcome_spectrum"
                                ],
                                "fallback_context_deduplicated": fallback_provenance[
                                    "deduplicated_identical_context"
                                ],
                                "fallback_shared_input_guide_ids": fallback_provenance[
                                    "shared_input_guide_ids"
                                ],
                                "analysis_role": fallback["analysis_role"],
                                "fallback_predictor": fallback["predictor"],
                                "fallback_status": FALLBACK_STATUS_EXECUTED,
                                "fallback_executed_date": FALLBACK_EXECUTED_DATE,
                                "primary_predictor": predictor.name,
                                "primary_predictor_status": exclusion["status"],
                                "primary_predictor_exclusion_link": {
                                    "ledger": "meta.prediction_exclusions",
                                    "guide_id": exclusion["guide_id"],
                                    "status": exclusion["status"],
                                },
                                "primary_predictor_exclusion_reason": exclusion["reason"],
                                "structural_stability": exclusion.get("structural_stability"),
                                "upstream_mass_anomaly": exclusion.get(
                                    "upstream_mass_anomaly"
                                ),
                                "raw_mass": exclusion.get("raw_mass"),
                                "raw_mass_per_member": exclusion.get("raw_mass_per_member"),
                                "normalization_applied": None,
                                "rejected_bedict_scores_used": False,
                                "fallback_blended_with_primary": False,
                                "reason": exclusion["reason"],
                            }
                    else:
                        records = structural_records_without_rates(
                            guide,
                            on_target_coords=on_coords,
                            reason=exclusion["reason"],
                            extra_provenance={
                                "predictor_runs": predictor.runs,
                                "normalization_applied": None,
                                "structural_stability": exclusion.get("structural_stability"),
                                "upstream_mass_anomaly": exclusion.get("upstream_mass_anomaly"),
                                "raw_mass": exclusion.get("raw_mass"),
                                "raw_mass_per_member": exclusion.get("raw_mass_per_member"),
                                "rate_semantics": None,
                                "analysis_role": fallback["analysis_role"],
                                "fallback_predictor": fallback["predictor"],
                                "fallback_status": fallback["status"],
                                "primary_predictor": predictor.name,
                                "primary_predictor_status": exclusion["status"],
                                "primary_predictor_exclusion_link": {
                                    "ledger": "meta.prediction_exclusions",
                                    "guide_id": exclusion["guide_id"],
                                    "status": exclusion["status"],
                                },
                            },
                        )
                    co_targets = guide.protospacer in protospacers_by_gene[other_gene]
                    ch_offtargets = ch_driver_offtargets(cand.offtarget_sites)
                    low_mm_sites = _low_mismatch_sites(cand.offtarget_sites)
                    for record in records:
                        record.provenance = {
                            **record.provenance,
                            "offtarget": cand.offtarget,
                            "offtarget_sites_low_mismatch": low_mm_sites,
                            "co_targets_paralog": co_targets,
                            "co_target_gene": other_gene if co_targets else None,
                            "ch_driver_offtargets": ch_offtargets,
                        }
                    all_records.extend(records)
                    continue
                cand = candidate_by_guide_id[guide.guide_id]
                spectrum = batch.spectra[guide.guide_id]
                on_coords = anchor_coords.get((gene, guide.editor_class), set())
                records = spectrum_to_records(spectrum, on_target_coords=on_coords)
                co_targets = guide.protospacer in protospacers_by_gene[other_gene]
                ch_offtargets = ch_driver_offtargets(cand.offtarget_sites)
                low_mm_sites = _low_mismatch_sites(cand.offtarget_sites)
                for record in records:
                    record.provenance = {
                        **record.provenance,
                        "offtarget": cand.offtarget,
                        "offtarget_sites_low_mismatch": low_mm_sites,
                        "co_targets_paralog": co_targets,
                        "co_target_gene": other_gene if co_targets else None,
                        "ch_driver_offtargets": ch_offtargets,
                    }
                all_records.extend(records)
        else:
            included_guide_count += len(editor_guides)
            for guide in editor_guides:
                cand = candidate_by_guide_id[guide.guide_id]
                on_coords = anchor_coords.get((gene, guide.editor_class), set())
                records = structural_records_without_rates(
                    guide,
                    on_target_coords=on_coords,
                    reason=NON_NGG_MODEL_SCOPE_REASON,
                )
                co_targets = guide.protospacer in protospacers_by_gene[other_gene]
                ch_offtargets = ch_driver_offtargets(cand.offtarget_sites)
                low_mm_sites = _low_mismatch_sites(cand.offtarget_sites)
                for record in records:
                    record.provenance = {
                        **record.provenance,
                        "offtarget": cand.offtarget,
                        "offtarget_sites_low_mismatch": low_mm_sites,
                        "co_targets_paralog": co_targets,
                        "co_target_gene": other_gene if co_targets else None,
                        "ch_driver_offtargets": ch_offtargets,
                    }
                all_records.extend(records)

    for record in all_records:
        record.consequence = NOT_EVALUATED
        record.provenance = {**record.provenance, "consequence_source": NOT_EVALUATED}
    all_records = annotate_ch(all_records)

    if write:
        on_target_n = sum(1 for r in all_records if r.on_target_or_bystander == "on_target")
        ch_genes_hit = sorted({
            h["ch_driver_gene"]
            for r in all_records
            for h in r.provenance.get("ch_driver_offtargets", [])
        })
        tools = [designer.tool_metadata()]
        if predictor is not None:
            tools.append(predictor.tool_metadata())
        if fallback_predictor is not None:
            tools.append(fallback_predictor.tool_metadata())
        repro = environment.reproducibility(
            tools=tools,
            editors=EDITOR_NAME,
            ensemble_runs=predictor.runs if predictor is not None else None,
            reference_windows={g: reference.window_provenance(c, s, e)
                               for g, (c, s, e) in WINDOWS.items()},
        )
        outputs.write_records(
            all_records,
            filename="enumeration_HBG1_HBG2.json",
            meta={
                "step": "2 (guides and their edits) + 4 (CH-driver annotation)",
                "windows": WINDOWS,
                "guide_design_tool": designer.name,
                "predictor": predictor.name if predictor is not None else NOT_EVALUATED,
                "prediction_scope_reason": (
                    None if predictor is not None else NON_NGG_MODEL_SCOPE_REASON
                ),
                "pam_class": designer.pam,
                "editors": EDITOR_NAME,
                "record_count": len(all_records),
                "on_target_records": on_target_n,
                "bystander_records": len(all_records) - on_target_n,
                "ch_driver_set": "boostDM-CH 12-gene panel (CC-BY-NC)",
                "ch_driver_offtarget_genes_hit": ch_genes_hit,
                "ambiguous_candidate_count": len(ambiguous),   # this many guides were dropped
                "ambiguous_candidates": ambiguous,
                "normalization_applied": NORMALIZATION_APPLIED if predictor is not None else None,
                "rate_semantics": RATE_SEMANTICS if predictor is not None else None,
                "structural_stability_required": predictor is not None,
                "raw_mass_allowed_range": (
                    [RAW_MASS_MIN, RAW_MASS_MAX] if predictor is not None else None
                ),
                "upstream_mass_anomaly_reference_range": (
                    [SUPERSEDED_RAW_MASS_MIN, SUPERSEDED_RAW_MASS_MAX]
                    if predictor is not None else None
                ),
                "included_guide_count": included_guide_count,
                "prediction_exclusion_count": len(prediction_exclusions),
                "prediction_exclusions": prediction_exclusions,
                "malformed_guide_exclusion_count": sum(
                    exclusion["status"] == "excluded_raw_mass_outside_authorized_range"
                    for exclusion in prediction_exclusions
                ),
                "malformed_guide_exclusions": [
                    exclusion for exclusion in prediction_exclusions
                    if exclusion["status"] == "excluded_raw_mass_outside_authorized_range"
                ],
                "fallback_route_count": len(fallback_routes),
                "fallback_routes": fallback_routes,
                "fallback_execution_count": fallback_execution_count,
                "fallback_executed_date": (
                    FALLBACK_EXECUTED_DATE if fallback_execution_count else None
                ),
                "reproducibility": repro,
            },
        )
    return all_records


if __name__ == "__main__":
    recs = run_enumeration()
    on = [r for r in recs if r.on_target_or_bystander == "on_target"]
    ch_hits = {}
    for r in recs:
        for h in r.provenance.get("ch_driver_offtargets", []):
            gid = r.provenance.get("guide_id")
            ch_hits.setdefault((r.gene, gid, h["ch_driver_gene"]),
                               (h["chrom"], h["start"], h["mismatch_count"]))
    print(f"total records: {len(recs)}  |  on-target records: {len(on)}  "
          f"|  edits inside CH-driver genes: {sum(1 for r in recs if r.ch_driver_flag)}")
    print("off-targets landing in CH-driver genes (guide -> gene):")
    for (gene, gid, ch_gene), (chrom, start, mm) in sorted(ch_hits.items()):
        print(f"  {gid} ({gene}) -> off-target in {ch_gene} at {chrom}:{start} ({mm} mismatches)")
    print("\ntargeted positions that were hit (guide -> position):")
    seen = set()
    for r in sorted(on, key=lambda r: (r.gene, r.promoter_offset)):
        key = (r.gene, r.promoter_offset, r.provenance.get("guide_id"))
        if key in seen:
            continue
        seen.add(key)
        ot = r.provenance.get("offtarget", {})
        print(
            f"  {r.gene} {-abs(r.promoter_offset):>4} ({r.chrom}:{r.pos_1based}) "
            f"{r.editor_class} rate={r.predicted_rate} via {r.provenance.get('guide_id')} "
            f"| also hits paralog={r.provenance.get('co_targets_paralog')} "
            f"| offtargets={ot.get('offtarget_count')} mitSpec={ot.get('mit_specificity')}"
        )
