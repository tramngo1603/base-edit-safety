"""Render aggregate summaries separately from protected donor-site tables."""
from __future__ import annotations

import csv
import io


def tsv(columns, rows):
    stream = io.StringIO()
    writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue()


def plan_tables(request, result, catalog=()):
    actions = {a["id"]: a for a in request["actions"]}
    specimens = {a["id"]: a["specimen_id"] for a in catalog}
    donors = {d["id"]: d for d in request["donors"]}
    rows, activated = [], set()
    for aid in result["selected_actions"]:
        a = actions[aid]
        activation = 0 if a["donor_id"] in activated else donors[a["donor_id"]]["activation_cost"]
        activated.add(a["donor_id"])
        rows.append({"action_id": aid, "donor_id": a["donor_id"], "specimen_id": specimens.get(aid), "assay_cost": a["assay_cost"],
                     "activation_cost_charged": activation, "covered_units": ",".join(a["covered_units"]),
                     "coverage_basis": request["coverage_basis"], "status": "planned"})
    from .planning import eligible_actions
    eligible = eligible_actions(request)
    uncovered = []
    for uid in result["uncovered_units"]:
        carriers = [a for a in eligible if uid in a["covered_units"]]
        reason = "budget_limits_or_distinct_donor_demand" if carriers else "no_eligible_action"
        uncovered.append({"unit_id": uid, "reason": reason})
    return {"donor-site-plan.tsv": tsv(["action_id", "donor_id", "specimen_id", "assay_cost", "activation_cost_charged", "covered_units", "coverage_basis", "status"], rows),
            "uncovered-configurations.tsv": tsv(["unit_id", "reason"], uncovered)}


def assay_review_table(request, result, catalog, candidates, guides, units, evidence, assay, variants):
    """Join protected source records into one row per proposed donor-site assay."""
    from .common import canonical
    actions = {a["id"]: a for a in catalog}
    candidates = {c["configuration_id"]: c for c in candidates}
    guides = {g["guide_id"]: g for g in guides}
    units = {u["id"]: u for u in units}
    variants = {v["variant_id"]: v for v in variants}
    coverage = {u["unit_id"]: u for u in result["unit_coverage"]}
    columns = ["action_id", "donor_id", "specimen_id", "unit_id", "coverage_basis", "status", "selector",
               "contig", "start_0based", "end_0based_exclusive", "strand", "guide_id", "guide_spacer", "editor",
               "target_spacer", "target_pam", "sequence_genomic_forward", "configuration_id", "target_path_json",
               "linked_variant_records_json", "assay_id", "biological_context_json", "declared_feasibility",
               "allele_resolution_required", "locus_resolution_required", "declared_allele_resolved", "declared_locus_resolved",
               "resolution_contract_id", "max_detection_limit_fraction", "required_distinct_donors", "prior_distinct_donors",
               "remaining_donor_requirement_before_plan", "excluded_evidence_reasons", "evidence_gap", "selection_rationale",
               "assay_cost", "interpretation"]
    rows = []
    for aid in result["selected_actions"]:
        action = actions[aid]
        # Current catalogs contain one configuration and one unit per assay.
        uid, = action["covered_units"]
        unit, candidate = units[uid], candidates[action["configuration_id"]]
        guide = guides[candidate["guide_id"]]
        prior = len(coverage[uid]["existing_donors"])
        missing = max(0, unit["required_donors"] - prior)
        linked_variants = [{"variant_id": vid, "contig": variants[vid]["contig"], "pos0": variants[vid]["pos0"],
                            "ref": variants[vid]["ref"], "original_alts": variants[vid]["original_alts"],
                            "selected_donor_call": variants[vid]["calls"].get(action["donor_id"])}
                           for vid in candidate["variant_links"]]
        reasons = sorted({reason for e in evidence if e["unit_id"] == uid and not e["usable"] for reason in e["reasons"]})
        gap = ("Conditional sequence-model demand; empirical measurements do not reduce it."
               if request["coverage_basis"] == "retained_callset_model" else
               "Matched usable measurements do not yet meet the required number of distinct donors.")
        rows.append({"action_id": aid, "donor_id": action["donor_id"], "specimen_id": action["specimen_id"],
            "unit_id": uid, "coverage_basis": request["coverage_basis"], "status": "planned", "selector": result["algorithm"],
            "contig": candidate["contig"], "start_0based": candidate["bounding_interval"][0],
            "end_0based_exclusive": candidate["bounding_interval"][1], "strand": candidate["strand"],
            "guide_id": candidate["guide_id"], "guide_spacer": guide["spacer_sequence"], "editor": guide["editor"],
            "target_spacer": candidate["protospacer_guide_orientation"], "target_pam": candidate["pam_guide_orientation"],
            "sequence_genomic_forward": candidate["sequence_genomic_forward"], "configuration_id": candidate["configuration_id"],
            "target_path_json": canonical(candidate["target_path"]).decode(),
            "linked_variant_records_json": canonical(linked_variants).decode(), "assay_id": action["assay_id"],
            "biological_context_json": canonical(assay["context"]).decode(), "declared_feasibility": assay["feasibility"],
            "allele_resolution_required": unit["allele_resolution_required"], "locus_resolution_required": unit["locus_resolution_required"],
            "declared_allele_resolved": assay["allele_resolved"], "declared_locus_resolved": assay["locus_resolved"],
            "resolution_contract_id": action["resolution_contract_id"],
            "max_detection_limit_fraction": assay.get("max_detection_limit_fraction"),
            "required_distinct_donors": unit["required_donors"], "prior_distinct_donors": prior,
            "remaining_donor_requirement_before_plan": missing, "excluded_evidence_reasons": ",".join(reasons),
            "evidence_gap": gap, "selection_rationale": f"Selected by {result['algorithm']} within the declared budget to address the remaining distinct-donor requirement; eligibility follows the declared carriage, specimen and assay contracts.",
            "assay_cost": action["assay_cost"], "interpretation": action["interpretation"]})
    return tsv(columns, rows)


def aggregate_report(config, search_audit, request, result, benchmark_summary, guides=(), intended=(), changes=(), actions=(), evidence=None):
    from collections import Counter
    exclusions = Counter()
    for row in search_audit["executions"]:
        exclusions.update({k: v for k, v in row["counts"].items() if k not in {"hits", "evaluated_windows"}})
    scope = ", ".join(f"{s['contig']}:[{s['start']},{s['end']})" for s in config["scope"]["requested_intervals"])
    p = config["planning_policy"]
    action_reasons = Counter(reason for a in actions for reason in a["eligibility_reasons"])
    guide_text = "; ".join(f"{g['spacer_sequence']} / {g['editor']} (positions {g['editing_window_start_1based']}–{g['editing_window_end_1based_inclusive']})" for g in guides)
    loci_text = "; ".join(f"{l['locus_name']}: {l['chrom']}:[{l['start_0based']},{l['end_0based_exclusive']}) {l['strand']}" for l in intended)
    frequencies = [r for r in search_audit["configuration_frequencies"] if r["basis"] == config["coverage_basis"]]
    paired = search_audit["paired_window_scope"]
    paired_counts = Counter()
    for row in paired["rows"]:
        paired_counts.update(row["unpaired_reasons"])
    lines = ["# Population search and proposed assays", "", f"Mode: {config['mode']}. Coverage basis: {config['coverage_basis']}.", "",
             f"Guide / editor and declared editing window: {guide_text}.", f"Intended loci: {loci_text}.",
             f"Variant source release: {config['input_sources']['variants']['release']}.", "",
             f"Search scope (zero-based, half-open): {scope}.",
             f"PAM patterns: {', '.join(config['search_profile']['pam_patterns'])}. Maximum spacer mismatches: {config['search_profile']['max_mismatches']}.",
             "Both strands; 20-base guides; three-base PAMs; no alignment bulges. Targets spanning scope boundaries are excluded.", "",
             f"Distinct nominated configurations: {search_audit['distinct_configurations']}. Off-target coverage units: {len(request['units'])}.",
             "Intended controls remain in search records and receive no off-target coverage credit.", "",
             f"Plan status: {result['plan_status']}. Proposed assays: {len(result['selected_actions'])}. Donors with proposed assays: {len(result['selected_donors'])}.",
             f"Activation and assay cost: {result['total_cost']} {request['cost_unit']}. Fixed study overhead: {p['fixed_overhead']}. Unused budget: {result['unused_budget']}.",
             f"Total declared cost including overhead: {result['total_cost']+p['fixed_overhead']}. Total study budget: {p['max_cost']}.",
             f"Fully satisfied configuration requirements after the plan: {len(result['fully_covered_units'])}/{len(request['units'])}. Unmet: {len(result['uncovered_units'])}.",
             f"Weighted planned gain: {result['gain_numerator']}/{result['gain_denominator']}. Existing distinct donor-unit credits: {sum(len(u['existing_donors']) for u in result['unit_coverage'])}.", "",
             "Planned coverage describes the proposed assays. Conditional model coverage assumes reference sequence at unlisted positions in the retained callset.", "",
             "Specimen availability and assay capability use the study's supplied declarations.",
             f"Evaluable donor–configuration pairs under the selected basis: {sum(r['evaluable_donors'] for r in frequencies)}. Unknown pairs: {sum(r['unknown_donors'] for r in frequencies)}.",
             f"Paired-evaluable nomination-comparison records: {sum(r['paired_evaluable'] for r in changes)}. Unpaired records: {sum(not r['paired_evaluable'] for r in changes)}. Counts are per donor and sequence configuration.", "",
             "## Paired search scope", "", paired["definition"], "",
             f"Coverage basis: {paired['coverage_basis']}. Union of oriented windows: {sum(r['union_oriented_windows'] for r in paired['rows'])}. "
             f"Paired and evaluable: {sum(r['paired_evaluable_oriented_windows'] for r in paired['rows'])}. "
             f"Unpaired: {sum(r['unpaired_oriented_windows'] for r in paired['rows'])}.", "",
             "| Unpaired reason | Oriented windows |", "| --- | ---: |"]
    lines.extend(f"| {key} | {value} |" for key, value in sorted(paired_counts.items()))
    lines += ["", "## Search exclusions", "", "| Reason | Oriented windows |", "| --- | ---: |"]
    lines.extend(f"| {key} | {value} |" for key, value in sorted(exclusions.items()))
    lines += ["", "## Assay eligibility exclusions", "", "Counts refer to donor–configuration actions; one action can have several reasons.", "",
              "| Reason | Actions |", "| --- | ---: |"]
    lines.extend(f"| {key} | {value} |" for key, value in sorted(action_reasons.items()))
    if evidence is not None:
        lines += ["", "## Measurement evidence", "",
                  f"Usable positive links: {evidence['usable_positive_links']}. Usable non-detection links: {evidence['usable_nondetection_links']}.",
                  f"Conflicting donor–unit pairs: {evidence['conflicting_donor_unit_pairs']}. Excluded links: {evidence['excluded_links']}.",
                  f"Unmatched source measurements: {evidence['unmatched_source_measurements']}.",
                  "Links are measurement-to-unit comparisons, not independent donors. Conflicting outcomes remain separate; no result is averaged.",
                  "Retained-callset model plans do not receive empirical measurement credit.", "",
                  "| Exclusion reason | Links |", "| --- | ---: |"]
        lines.extend(f"| {key} | {value} |" for key, value in evidence["excluded_link_reasons"].items())
        lines += ["", "| Unmatched reason | Source measurements |", "| --- | ---: |"]
        lines.extend(f"| {key} | {value} |" for key, value in evidence["unmatched_source_reasons"].items())
    lines += ["", "## Selector comparison", "", "All selectors use the same universe, carriage, action eligibility, costs and evidence. This is a transductive study.", "",
              "| Selector | Study budget | Runs | Mean used action budget | Mean fully satisfied units | All units | Actionable units |", "| --- | ---: | ---: | ---: | ---: | ---: | ---: |"]
    groups = {}
    for row in benchmark_summary["rows"]:
        groups.setdefault((row["algorithm"], row["study_budget"]), []).append(row)
    for (algorithm, budget), rows in sorted(groups.items()):
        lines.append(f"| {algorithm} | {budget} | {len(rows)} | {sum(r['activation_and_assay_cost'] for r in rows)/len(rows):.3f} | {sum(r['fully_satisfied_all_units'] for r in rows)/len(rows):.3f} | {rows[0]['all_units']} | {rows[0]['actionable_units']} |")
    strata = benchmark_summary["rare_variant_strata"]
    lines += ["", "## Retained variant strata", ""]
    if strata["status"] == "reported":
        lines += [strata["scope"], "", "| Stratum | Subset rule | Retained ALT alleles |", "| --- | --- | ---: |"]
        lines.extend(f"| {label} | {strata['definitions'][label]} | {count} |"
                     for label, count in strata["retained_alt_alleles"].items())
        lines += ["", "Per-selector results retain all-unit and actionable-unit coverage within each linked-variant stratum. "
                  "These descriptive strata do not change weights, the frozen universe or assay selection.", "",
                  "| Selector | Study budget | Linked-variant stratum | All units | Actionable units | Mean fully satisfied units |", "| --- | ---: | --- | ---: | ---: | ---: |"]
        for (algorithm, budget), rows in sorted(groups.items()):
            for label in strata["definitions"]:
                values = [s for r in rows for s in r["variant_strata"] if s["stratum"] == label]
                if values and values[0]["all_units"]:
                    lines.append(f"| {algorithm} | {budget} | {label} | {values[0]['all_units']} | {values[0]['actionable_units']} | "
                                 f"{sum(v['fully_satisfied_all_units'] for v in values)/len(values):.3f} |")
    else:
        lines.append("Variant strata were not computed because source files were not supplied.")
    lines += ["", "Non-detection applies only to the matching assay and its detection limit.",
              "This summary omits donor IDs, sample IDs, individual memberships and local paths. Aggregates still require source-term and sharing review.", ""]
    return "\n".join(lines)
