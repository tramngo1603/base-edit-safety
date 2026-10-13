"""Command line entry point for source-bound donor-site planning."""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import sys

from ..population.common import PopulationError, identity, json_text, read_json, read_jsonl, record, write_bundle
from ..population.contracts import ALGORITHMS
from ..population.planning import plan
from ..population.provenance import resolve_file, runtime
from ..population.reporting import plan_tables, assay_review_table
from ..population.resolver import resolve_request


def execute(request_path, manifest_path, algorithm, output_dir, seed=None, overwrite=False):
    request, manifest, inputs = resolve_request(request_path, manifest_path)
    result, order = plan(request, algorithm, seed)
    result_text = json_text(result)
    provenance = record("donor_panel_provenance", manifest["run_id"],
        input_hashes=inputs, result_sha256=hashlib.sha256(result_text.encode()).hexdigest(),
        algorithm=algorithm, seed=seed, visit_order=order, runtime=runtime(),
        fingerprint=identity("plan", {"request": request, "algorithm": algorithm, "seed": seed}),
        validation_outcomes={"source_checksums": "passed", "identifier_resolution": "passed",
                             "coverage_and_eligibility_recomputed": "synthetic_declaration" if manifest["synthetic_declaration"] else "passed",
                             "budget_and_distinct_donor_accounting": "passed"})
    catalog = next((read_jsonl(resolve_file(ref, Path(manifest_path).parent)) for ref in manifest["artifacts"]
                    if ref["role"] == "action_catalog"), [])
    files = {"panel-plan.json": result_text, "plan-provenance.json": json_text(provenance), **plan_tables(request, result, catalog)}
    protected = [request_path, manifest_path]
    protected += [Path(manifest_path).parent / ref["path"] for ref in manifest["artifacts"]]
    if manifest["source_population_manifest"]:
        upstream_path = Path(manifest_path).parent/manifest["source_population_manifest"]["path"]
        protected.append(upstream_path)
        protected.extend(r["path"] for r in read_json(upstream_path)["sources"].values())
        refs = {ref["role"]: resolve_file(ref, Path(manifest_path).parent) for ref in manifest["artifacts"]}
        from ..population.provenance import load_stage, resolve_study
        population, _ = load_stage(resolve_study(refs["config"]), upstream_path.parent, "build")
        files["assay-review.tsv"] = assay_review_table(request, result, catalog,
            read_jsonl(refs["candidate_configurations"]), population["guides.json"], read_jsonl(refs["coverage_units"]),
            read_jsonl(refs["evidence_links"]), population["assay-resolution.json"], population["normalized-variants.jsonl"])
    if manifest["synthetic_declaration"]:
        protected.append(Path(manifest_path).parent/manifest["synthetic_declaration"]["path"])
    write_bundle(output_dir, files, inputs=protected, overwrite=overwrite)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--request", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--algorithm", choices=sorted(ALGORITHMS), default="greedy_cost")
    parser.add_argument("--seed", type=int)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(argv)
    try:
        result = execute(args.request, args.manifest, args.algorithm, args.output_dir, args.seed, args.overwrite)
    except (PopulationError, OSError, ValueError, KeyError, TypeError) as exc:
        print("Planning failed: " + (str(exc) if isinstance(exc, PopulationError) else "invalid or unavailable input"), file=sys.stderr)
        return 2
    print(f"Plan {result['plan_status']}: {len(result['selected_actions'])} proposed assays; cost {result['total_cost']}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
