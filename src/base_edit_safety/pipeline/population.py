"""Bounded population search, evidence reconciliation and assay catalogs."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys
import time

from .candidate_panel import CandidatePanelError, load_guides, load_intended_loci, _has_reported_text
from ..population.actions import make_actions, make_request, make_units
from ..population.benchmarks import benchmark
from ..population.common import (PopulationError, ResourceLimit, check_output, checksum, json_text,
    jsonl_text, read_json, read_jsonl, record, write_file)
from ..population.evidence import reconcile, resolve_measurements, representation_replay, unmatched_measurements, evidence_summary
from ..population.haplotypes import build_materializations
from ..population.identity import pam_matches, reverse_complement
from ..population.provenance import artifact_file, load_stage, resolve_file, resolve_study, save_stage
from ..population.reporting import aggregate_report
from ..population.resolver import resolve_request
from ..population.search import nomination_changes, scan
from ..population.variants import callable_intervals, canonical_donors, donors_from_tsv, load_vcf


def definitions(study):
    config, paths = study["config"], study["paths"]
    guides = load_guides(paths["guides"])
    loci = load_intended_loci(paths["intended_loci"])
    for guide in guides:
        if not _has_reported_text(guide["editor"]):
            raise PopulationError("Guide requires a specific editor definition")
        if len(guide["spacer_sequence"]) != 20 or len(guide["pam_class"]) != 3:
            raise PopulationError("Local search requires 20-base guides and three-base PAMs")
        if guide["editor_class"] not in {"ABE", "CBE", "nuclease"}:
            raise PopulationError("Unknown editor class")
        group = [l for l in loci if l["guide_id"] == guide["guide_id"]]
        if not group:
            raise PopulationError("Guide has no intended locus declaration")
        if guide["locus_system"] == "HBG1_HBG2":
            by_name = {l["locus_name"]: l for l in group}
            if {l["locus_name"] for l in group} != {"HBG1", "HBG2"} or len(group) != 2 \
                    or group[0]["chrom"] != group[1]["chrom"] \
                    or group[0]["strand"] != group[1]["strand"] \
                    or by_name["HBG2"]["start_0based"]-by_name["HBG1"]["start_0based"] != 4924:
                raise PopulationError("HBG intended loci must preserve the two distinct coordinates and 4,924-base spacing")
        for locus in group:
            if locus["genome_build"] != config["reference_policy"]["genome_build"]:
                raise PopulationError("Intended-locus build mismatch")
            seq = study["reference"].fetch(locus["chrom"], locus["start_0based"], locus["end_0based_exclusive"])
            seq = seq if locus["strand"] == "+" else reverse_complement(seq)
            if len(seq) != 23 or seq[:20] != guide["spacer_sequence"] or not pam_matches(seq[20:], [guide["pam_class"]]):
                raise PopulationError("Intended guide sequence or PAM disagrees with pinned reference")
    if set(config["assay_policy"]["guide_ids"]) - {g["guide_id"] for g in guides}:
        raise PopulationError("Assay refers to an unknown guide")
    if any(l["guide_id"] not in {g["guide_id"] for g in guides} for l in loci):
        raise PopulationError("Intended-locus guide does not resolve")
    return guides, loci


def build(study, output, overwrite=False):
    started = time.perf_counter()
    c, p, ref, rid = study["config"], study["paths"], study["reference"], study["run_id"]
    guides, loci = definitions(study)
    donors = donors_from_tsv(p["donors"], c["coverage_basis"])
    if len(canonical_donors(donors)) > c["resource_limits"]["max_donors"]:
        raise ResourceLimit("max_donors exceeded before genotype ingestion")
    variants, audit, ingestion = load_vcf(p["variants"], c, ref, donors, rid)
    calls = callable_intervals(p["callable_intervals"], donors, ref, c["input_sources"].get("callable_intervals"))
    materials, counts = build_materializations(c, ref, variants, donors, rid)
    return save_stage(study, "build", output, {
        "donor-inventory.jsonl": [record("donor_inventory", rid, **d) for d in donors],
        "guides.json": guides, "intended-loci.json": loci,
        "assay-resolution.json": record("assay_resolution_contract", rid, **c["assay_policy"]),
        "normalized-variants.jsonl": variants, "genotype-audit.jsonl": audit,
        "reference-callability.jsonl": [record("reference_callability", rid, **r) for r in calls],
        "sequence-configurations.jsonl": materials,
        "population-audit.json": record("population_audit", rid, ingestion=ingestion, reconstruction=counts)},
        overwrite=overwrite, timing=time.perf_counter()-started)


def search(study, population_dir, output, overwrite=False):
    started = time.perf_counter()
    pop, _ = load_stage(study, population_dir, "build")
    c, rid = study["config"], study["run_id"]
    candidates, hits, carriage, audit = scan(c, study["reference"], pop["sequence-configurations.jsonl"],
        pop["guides.json"], pop["intended-loci.json"], pop["donor-inventory.jsonl"], pop["reference-callability.jsonl"], rid)
    changes = nomination_changes(c, candidates, pop["sequence-configurations.jsonl"], pop["donor-inventory.jsonl"], pop["guides.json"], carriage, rid,
                                 pop["reference-callability.jsonl"])
    return save_stage(study, "search", output, {"candidate-configurations.jsonl": candidates, "search-hits.jsonl": hits,
        "donor-carriage.jsonl": carriage, "search-audit.json": audit, "nomination-changes.jsonl": changes},
        parents=[population_dir], overwrite=overwrite, timing=time.perf_counter()-started)


def reconcile_stage(study, search_dir, output, overwrite=False):
    artifacts, _ = load_stage(study, search_dir, "search")
    c, rid, paths = study["config"], study["run_id"], study["paths"]
    candidates = artifacts["candidate-configurations.jsonl"]
    units = make_units(c, candidates, rid)
    measurements = resolve_measurements(read_json(paths["candidate_ledger"]) if paths["candidate_ledger"] else None,
        read_json(paths["measurement_linkage"]) if paths["measurement_linkage"] else None,
        study["sources"].get("candidate_ledger", {}).get("sha256"))
    guides, _ = definitions(study)
    links = reconcile(c, units, candidates, guides, measurements, rid)
    unmatched = unmatched_measurements(units, candidates, guides, measurements, links, rid)
    replay = representation_replay(read_json(paths["candidate_ledger"]), study["reference"], guides, rid) if paths["candidate_ledger"] else []
    return save_stage(study, "reconcile", output, {"coverage-units.jsonl": units, "evidence-links.jsonl": links,
        "representation-replay.jsonl": replay, "unmatched-measurements.json": unmatched,
        "evidence-audit.json": record("evidence_audit", rid, source_measurements=len(measurements),
            usable_donor_unit_links=sum(e["usable"] for e in links),
            out_of_universe_measurements=len(unmatched["rows"]), summary=evidence_summary(links, unmatched["rows"]),
            source_measurements_retained=measurements)}, parents=[search_dir], overwrite=overwrite)


def actions_stage(study, population_dir, search_dir, evidence_dir, output, overwrite=False):
    pop, _ = load_stage(study, population_dir, "build")
    found, _ = load_stage(study, search_dir, "search")
    evidence, _ = load_stage(study, evidence_dir, "reconcile")
    c, rid = study["config"], study["run_id"]
    candidates, units, donors = found["candidate-configurations.jsonl"], evidence["coverage-units.jsonl"], pop["donor-inventory.jsonl"]
    actions = make_actions(c, units, candidates, donors, found["donor-carriage.jsonl"], rid)
    request = make_request(c, units, donors, actions, evidence["evidence-links.jsonl"])
    output = Path(output).resolve()
    manifest_path = output.parent/"planning-manifest.json"
    catalog = output.parent/"assay-actions.jsonl"
    targets = [output, manifest_path, catalog]
    if len(set(targets)) != len(targets):
        raise PopulationError("Planner output paths collide")
    inputs = [study["config_path"], population_dir, search_dir, evidence_dir] + [s["path"] for s in study["sources"].values()]
    for target in targets:
        check_output(target, inputs, overwrite)
    write_file(output, json_text(request), inputs=inputs, overwrite=overwrite)
    write_file(catalog, jsonl_text(actions), inputs=inputs, overwrite=overwrite)
    refs = {"config": study["config_path"], "coverage_units": Path(evidence_dir)/"coverage-units.jsonl",
        "donor_carriage": Path(search_dir)/"donor-carriage.jsonl", "donor_inventory": Path(population_dir)/"donor-inventory.jsonl",
        "assay_resolution": Path(population_dir)/"assay-resolution.json", "action_catalog": catalog,
        "evidence_links": Path(evidence_dir)/"evidence-links.jsonl", "candidate_configurations": Path(search_dir)/"candidate-configurations.jsonl",
        "search_manifest": Path(search_dir)/"run-manifest.json", "evidence_manifest": Path(evidence_dir)/"run-manifest.json"}
    manifest = {"record_type": "donor_panel_manifest", "schema_version": "1", "run_id": rid,
        "mode": c["mode"], "coverage_basis": c["coverage_basis"], "request": artifact_file(output, output.parent),
        "synthetic_declaration": None, "source_population_manifest": artifact_file(Path(population_dir)/"run-manifest.json", output.parent),
        "artifacts": [artifact_file(path, output.parent, role) for role, path in sorted(refs.items())]}
    write_file(manifest_path, json_text(manifest), inputs=inputs, overwrite=overwrite)
    resolve_request(output, manifest_path)
    return request


def report_stage(study, run_dir, output, overwrite=False):
    root = Path(run_dir)
    population_data, _ = load_stage(study, root/"population", "build")
    search_data, _ = load_stage(study, root/"search", "search")
    evidence_data, _ = load_stage(study, root/"evidence", "reconcile")
    benchmark_data, _ = load_stage(study, root/"benchmark", "benchmark")
    request, manifest, _ = resolve_request(root/"panel-request.json", root/"planning-manifest.json")
    plan = read_json(root/"plan"/"panel-plan.json")
    provenance = read_json(root/"plan"/"plan-provenance.json")
    if provenance.get("run_id") != study["run_id"] or provenance.get("result_sha256") != checksum(root/"plan"/"panel-plan.json") \
            or provenance.get("input_hashes", {}).get("request") != checksum(root/"panel-request.json") \
            or provenance.get("input_hashes", {}).get("manifest") != checksum(root/"planning-manifest.json"):
        raise PopulationError("Plan provenance is stale or does not match the reported result")
    from ..population.planning import plan as calculate
    if calculate(request, plan["algorithm"], provenance["seed"])[0] != plan:
        raise PopulationError("Reported plan differs from recomputed result")
    text = aggregate_report(study["config"], search_data["search-audit.json"], request, plan, benchmark_data["benchmark-summary.json"],
                            population_data["guides.json"], population_data["intended-loci.json"],
                            search_data["nomination-changes.jsonl"], read_jsonl(root/"assay-actions.jsonl"),
                            evidence_data["evidence-audit.json"]["summary"])
    write_file(output, text, inputs=[root/"plan", root/"search", root/"benchmark", root/"population", root/"evidence",
                                   root/"panel-request.json", root/"planning-manifest.json", root/"assay-actions.jsonl",
                                   study["config_path"]] + [s["path"] for s in study["sources"].values()], overwrite=overwrite)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    for command in ("validate", "build", "search", "reconcile", "actions", "benchmark", "report"):
        p = sub.add_parser(command)
        p.add_argument("--config", type=Path, required=True)
        p.add_argument("--overwrite", action="store_true")
        if command in {"search", "actions"}:
            p.add_argument("--population-dir", type=Path, required=True)
        if command in {"reconcile", "actions"}:
            p.add_argument("--search-dir", type=Path, required=True)
        if command == "actions":
            p.add_argument("--evidence-dir", type=Path, required=True)
        if command in {"actions", "report"}:
            p.add_argument("--output", type=Path, required=True)
        elif command != "validate":
            p.add_argument("--output-dir", type=Path, required=True)
        if command == "benchmark":
            p.add_argument("--request", type=Path, required=True)
        if command == "report":
            p.add_argument("--run-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    study = None
    try:
        study = resolve_study(args.config)
        if args.command == "validate":
            definitions(study)
            donors = donors_from_tsv(study["paths"]["donors"], study["config"]["coverage_basis"])
            load_vcf(study["paths"]["variants"], study["config"], study["reference"], donors, study["run_id"])
        elif args.command == "build":
            build(study, args.output_dir, args.overwrite)
        elif args.command == "search":
            search(study, args.population_dir, args.output_dir, args.overwrite)
        elif args.command == "reconcile":
            reconcile_stage(study, args.search_dir, args.output_dir, args.overwrite)
        elif args.command == "actions":
            actions_stage(study, args.population_dir, args.search_dir, args.evidence_dir, args.output, args.overwrite)
        elif args.command == "benchmark":
            request, manifest, _ = resolve_request(args.request, args.request.parent/"planning-manifest.json")
            if manifest["run_id"] != study["run_id"]:
                raise PopulationError("Benchmark request belongs to a different study")
            pop_path = resolve_file(manifest["source_population_manifest"], args.request.parent)
            search_path = resolve_file(next(r for r in manifest["artifacts"] if r["role"] == "search_manifest"), args.request.parent)
            population_data, _ = load_stage(study, pop_path.parent, "build")
            search_data, _ = load_stage(study, search_path.parent, "search")
            result = benchmark(study["config"], request, study["run_id"],
                variants=population_data["normalized-variants.jsonl"], candidates=search_data["candidate-configurations.jsonl"])
            protected = [args.request, args.request.parent/"planning-manifest.json"]
            protected += [args.request.parent/ref["path"] for ref in manifest["artifacts"]]
            save_stage(study, "benchmark", args.output_dir, {"benchmark-summary.json": result},
                       input_files=protected, overwrite=args.overwrite)
        else:
            report_stage(study, args.run_dir, args.output, args.overwrite)
    except ResourceLimit as exc:
        if study is not None and hasattr(args, "output_dir"):
            try:
                save_stage(study, args.command, args.output_dir,
                    {"execution-failure.json": record("execution_failure", study["run_id"], code="resource_cap",
                        reason=str(exc), affected_scope=study["config"]["scope"], row_retention="none",
                        configuration_coverage="unknown", search_completeness="not_established")},
                    overwrite=args.overwrite, status="partial")
            except (PopulationError, OSError):
                pass
        print("Partial execution: " + str(exc), file=sys.stderr)
        return 4
    except (PopulationError, CandidatePanelError, OSError, ValueError, KeyError, TypeError) as exc:
        print("Population input rejected: " + (str(exc) if isinstance(exc, PopulationError) else "invalid or unavailable input; consult the declared schema"), file=sys.stderr)
        return 2
    print("Population " + args.command + " completed within declared scope.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
