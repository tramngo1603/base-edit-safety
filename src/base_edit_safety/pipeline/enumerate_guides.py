"""Enumerate guides across both promoter windows and predict their edited sequences.

For the HBG1 and HBG2 promoter windows this runs the external guide-design tool, which also does a
genome-wide off-target search. For every guide the editor can act on, it predicts the edit outcomes,
turns those into per-position EditRecords, and marks which positions this tool was aiming for. It
then attaches the off-target summary, notes whether the guide also hits the paralog, flags
off-targets that land in CH-driver genes, leaves the consequence field for the optional annotation
step, and writes the records to `data/outputs`.
"""
from __future__ import annotations

from collections import defaultdict
import math
import sys

from .. import anchors
from ..adapters.bedict import (
    NORMALIZATION_APPLIED,
    RAW_MASS_MAX,
    RAW_MASS_MIN,
    BedictBystanderPredictor,
)
from ..adapters.consequence import mark_consequence_deferred
from ..adapters.crispor import CrisporGuideDesigner
from ..annotation import annotate_records as annotate_ch, ch_driver_offtargets
from ..io import outputs
from ..io import reference
from . import environment
from .datatypes import EditRecord, spectrum_to_records

# Store per-site coordinates for sites with up to three mismatches. Four-mismatch sites account for
# 142 of 155 sites in one guide and receive a weight of 0.002 in the rubric. The retained coordinates
# support comparisons with published panels and site-specific assays; summary counts still include
# the four-mismatch sites.
PERSIST_SITES_UP_TO_MISMATCHES = 3

# Fixed promoter windows from -210 to -80 relative to the transcription start site on GRCh38. See
# `docs/coordinates.md`.
WINDOWS = {
    "HBG1": ("chr11", 5_249_937, 5_250_067),
    "HBG2": ("chr11", 5_254_861, 5_254_991),
}
EDITOR_NAME = {"ABE": "ABE8e", "CBE": "BE4max"}


def _percentile(sorted_values: list[float], fraction: float) -> float:
    """Linearly interpolated percentile over an already sorted, non-empty sequence."""
    position = (len(sorted_values) - 1) * fraction
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return sorted_values[lower]
    weight = position - lower
    return sorted_values[lower] * (1.0 - weight) + sorted_values[upper] * weight


def _raw_mass_distribution(guides: list[dict]) -> dict:
    """Summarize the predictor's input score totals and scaled values for one run."""
    ordered_guides = sorted(guides, key=lambda row: row["guide_id"])
    values = sorted(float(row["raw_mass"]) for row in ordered_guides)
    if not values:
        return {
            "normalization_applied": NORMALIZATION_APPLIED,
            "allowed_raw_mass": {"minimum": RAW_MASS_MIN, "maximum": RAW_MASS_MAX},
            "guide_count": 0,
            "included_guide_count": 0,
            "excluded_guide_count": 0,
            "summary": None,
            "histogram": [],
            "guides": [],
        }

    mean = math.fsum(values) / len(values)
    population_stddev = math.sqrt(
        math.fsum((value - mean) ** 2 for value in values) / len(values)
    )
    histogram_specs = [
        ("<0.75", lambda value: value < 0.75),
        ("[0.75,0.90)", lambda value: 0.75 <= value < 0.90),
        ("[0.90,0.95)", lambda value: 0.90 <= value < 0.95),
        ("[0.95,1.00)", lambda value: 0.95 <= value < 1.00),
        ("[1.00,1.05)", lambda value: 1.00 <= value < 1.05),
        ("[1.05,1.10)", lambda value: 1.05 <= value < 1.10),
        ("[1.10,1.25]", lambda value: 1.10 <= value <= 1.25),
        (">1.25", lambda value: value > 1.25),
    ]
    return {
        "normalization_applied": NORMALIZATION_APPLIED,
        "allowed_raw_mass": {"minimum": RAW_MASS_MIN, "maximum": RAW_MASS_MAX},
        "guide_count": len(ordered_guides),
        "included_guide_count": sum(row["status"] == "included" for row in ordered_guides),
        "excluded_guide_count": sum(row["status"] != "included" for row in ordered_guides),
        "summary": {
            "minimum": values[0],
            "q1": _percentile(values, 0.25),
            "median": _percentile(values, 0.50),
            "q3": _percentile(values, 0.75),
            "maximum": values[-1],
            "mean": mean,
            "population_stddev": population_stddev,
        },
        "histogram": [
            {"raw_mass_range": label, "guide_count": sum(predicate(value) for value in values)}
            for label, predicate in histogram_specs
        ],
        "guides": ordered_guides,
    }


def _anchor_coords_by_gene_editor() -> dict[tuple[str, str], set[int]]:
    """The genomic coordinates this tool targets, keyed by (gene, editor class)."""
    coords: dict[tuple[str, str], set[int]] = defaultdict(set)
    for pos in anchors.IN_SCOPE:
        coords[("HBG1", pos.editor_class)].add(pos.hbg1_pos)
        coords[("HBG2", pos.editor_class)].add(pos.hbg2_pos)
    return coords


def _low_mismatch_sites(sites: list[dict]) -> list[dict]:
    """Return coordinates and comparison fields for sites with up to three mismatches."""
    kept = []
    for s in sites:
        try:
            mm = int(s.get("mismatchCount"))
        except (TypeError, ValueError):
            continue
        if mm == 0 or mm > PERSIST_SITES_UP_TO_MISMATCHES:
            continue          # 0-mm is the intended paralog, not an off-target
        kept.append({"chrom": s.get("chrom"), "start": s.get("start"), "strand": s.get("strand"),
                     "mismatch_count": mm, "cfd_offtarget_score": s.get("cfdOfftargetScore"),
                     "locus_desc": s.get("locusDesc")})
    return sorted(kept, key=lambda x: (x["mismatch_count"], str(x["chrom"]), str(x["start"])))


def run_enumeration(write: bool = True) -> list[EditRecord]:
    """Find every guide over both promoters and return all the annotated EditRecords."""
    anchors.check_paralog_spacing()  # stop before prediction if a coordinate has shifted
    designer = CrisporGuideDesigner()
    predictor = BedictBystanderPredictor()
    anchor_coords = _anchor_coords_by_gene_editor()

    # Find candidates per gene, and note which protospacers show up in both promoters.
    candidates_by_gene = {}
    protospacers_by_gene: dict[str, set[str]] = {}
    ambiguous: list[str] = []   # guides excluded after mapping to more than one location
    raw_mass_guides: list[dict] = []
    malformed_guide_exclusions: list[dict] = []
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

        batch = predictor.predict_many(editor_guides)
        for spectrum in batch.spectra.values():
            raw_mass_guides.append({
                "guide_id": spectrum.guide.guide_id,
                "gene": spectrum.guide.gene,
                "editor_class": spectrum.guide.editor_class,
                "editor": spectrum.guide.editor_name,
                "raw_mass": spectrum.provenance["raw_mass"],
                "status": "included",
            })
        for exclusion in batch.exclusions:
            malformed_guide_exclusions.append(exclusion)
            raw_mass_guides.append({
                "guide_id": exclusion["guide_id"],
                "gene": exclusion["gene"],
                "editor_class": exclusion["editor_class"],
                "editor": exclusion["editor"],
                "raw_mass": exclusion["raw_mass"],
                "status": exclusion["status"],
            })
            print(f"EXCLUDED: {exclusion['reason']}", file=sys.stderr)

        for guide in editor_guides:
            if guide.guide_id not in batch.spectra:
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

    all_records = mark_consequence_deferred(all_records)
    all_records = annotate_ch(all_records)

    if write:
        on_target_n = sum(1 for r in all_records if r.on_target_or_bystander == "on_target")
        ch_genes_hit = sorted({
            h["ch_driver_gene"]
            for r in all_records
            for h in r.provenance.get("ch_driver_offtargets", [])
        })
        repro = environment.reproducibility(
            tools=[designer.tool_metadata(), predictor.tool_metadata()],
            editors=EDITOR_NAME,
            ensemble_runs=predictor.runs,
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
                "predictor": predictor.name,
                "editors": EDITOR_NAME,
                "record_count": len(all_records),
                "on_target_records": on_target_n,
                "bystander_records": len(all_records) - on_target_n,
                "ch_driver_set": "boostDM-CH 12-gene panel (CC-BY-NC)",
                "ch_driver_offtarget_genes_hit": ch_genes_hit,
                "ambiguous_candidate_count": len(ambiguous),   # number of excluded ambiguous guides
                "ambiguous_candidates": ambiguous,
                "normalization_applied": NORMALIZATION_APPLIED,
                "raw_mass_allowed_range": [RAW_MASS_MIN, RAW_MASS_MAX],
                "raw_mass_distribution": _raw_mass_distribution(raw_mass_guides),
                "malformed_guide_exclusion_count": len(malformed_guide_exclusions),
                "malformed_guide_exclusions": malformed_guide_exclusions,
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
    print(f"\ntotal records: {len(recs)}  |  on-target records: {len(on)}")
    print("targeted positions that were hit (guide -> position):")
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
