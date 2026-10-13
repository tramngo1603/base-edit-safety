"""Run source-bound comparative studies and plan complete research experiments."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import sys

from ..population.common import (PopulationError, ResourceLimit, check_output,
    checksum, json_text, write_bundle)
from ..population.comparison_provenance import (
    load_stage, load_stage_records, resolve_study, save_stage,
)


class _Parser(argparse.ArgumentParser):
    def error(self, message):
        # Argument values may be protected paths or identifiers.
        raise PopulationError("Command arguments are invalid; use --help for required options")


def _stage_arguments(commands, name, help_text, parents=()):
    command = commands.add_parser(name, help=help_text)
    command.add_argument("--config", required=True)
    command.add_argument("--output", "--output-dir", dest="output", required=True)
    command.add_argument("--reuse", action="store_true", help="Verify and reuse an existing immutable stage")
    for parent in parents:
        command.add_argument("--" + parent, required=True)
    return command


def _parent(directory, manifest, stage):
    if stage not in manifest["parents"]:
        raise PopulationError("Required comparison stage prerequisite is absent")
    return (Path(directory) / manifest["parents"][stage]["path"]).resolve().parent


def _same_parent(directory, manifest, stage, supplied):
    expected = manifest["parents"].get(stage)
    if expected is None or checksum(Path(supplied) / "run-manifest.json") != expected["sha256"]:
        raise PopulationError("Supplied stage differs from the verified prerequisite")


def _question_set(study, directory, evidence=None):
    from ..population.comparison_questions import prepare_question_set
    records, artifacts, manifest = load_stage_records(study, directory, "questions")
    comparison_path = _parent(directory, manifest, "comparison")
    _, comparison_manifest = load_stage(study, comparison_path, "comparison")
    if evidence is not None:
        _same_parent(comparison_path, comparison_manifest, "evidence_applicability", evidence)
    return prepare_question_set(study, records), artifacts, manifest


def _previous_plan(args, study):
    if bool(args.previous_plan) != bool(args.previous_config):
        raise PopulationError("A previous plan requires both --previous-plan and --previous-config")
    if not args.previous_plan:
        return None
    paths = [Path(args.previous_config).resolve(),
             Path(args.previous_plan).resolve() / "run-manifest.json",
             Path(args.previous_plan).resolve() / "plan.json"]
    declared = {Path(source["path"]) for source in study["sources"].values()
                if source["role"] == "plan_history"}
    if not set(paths) <= declared:
        raise PopulationError("Previous configuration, manifest and plan must be pinned plan_history inputs")
    previous_study = resolve_study(args.previous_config)
    artifacts, manifest = load_stage(previous_study, args.previous_plan, "planning")
    qset, _, _ = _question_set(previous_study,
        _parent(args.previous_plan, manifest, "questions"),
        _parent(args.previous_plan, manifest, "evidence_applicability"))
    from ..population.experiment_planning import verify_plan_revision
    verify_plan_revision(qset, artifacts["plan.json"])
    return artifacts["plan.json"]["plan_id"]


def _search_backend(args, study):
    active = set(study["config"]["policies"]["search"]["profile_ids"])
    names = {row["backend_name"] for row in study["config"]["catalogs"]
             if row["record_id"] in active and row["record_type"] == "comparison_search_profile"}
    backend = args.backend or ("local" if names == {"comparative_local"} else
                               "crisprme_plus" if names == {"crisprme_plus"} else None)
    expected = "comparative_local" if backend == "local" else "crisprme_plus"
    if backend is None or names != {expected}:
        raise PopulationError("Selected backend does not match every active search profile")
    return backend


def _search(args, study):
    from ..population.comparison import HostBudget, host_limits_from_study, local_search, normalize_search
    from ..population.comparison_inputs import prepare_inputs
    backend = _search_backend(args, study)
    if backend == "local":
        if args.native_output:
            raise PopulationError("--native-output applies only to the configured native backend")
        return local_search(study)
    from ..population.external_search import native_limits_from_study, search_prepared
    if not args.native_output:
        raise PopulationError("Native search requires a separate --native-output audit directory")
    check_output(args.native_output, [study["config_path"], args.output,
        *[source["path"] for source in study["sources"].values()]], overwrite=True)
    budget = HostBudget(host_limits_from_study(study))
    budget.check("source_preparation")
    prepared = prepare_inputs(study)
    budget.check("source_preparation")
    result = search_prepared(prepared, output_dir=args.native_output,
                            limits=native_limits_from_study(study), host_budget=budget)
    hits, execution = result["hits"], result["execution"]
    del result
    return normalize_search(prepared, hits, execution, host_budget=budget)


def _execute_stage(args, study):
    stages = {"search": "search", "observations": "observation_normalization",
              "evidence": "evidence_applicability", "compare": "comparison",
              "questions": "questions", "plan": "planning", "report": "report"}
    stage = stages[args.command]
    if args.reuse:
        if stage == "search":
            _search_backend(args, study)
        artifacts, manifest = load_stage(study, args.output, stage)
        for parent in ("search", "observations", "evidence", "comparison", "questions", "plan"):
            supplied = getattr(args, parent, None)
            actual = {"observations": "observation_normalization", "evidence": "evidence_applicability",
                      "plan": "planning"}.get(parent, parent)
            if supplied is not None and actual in manifest["parents"]:
                _same_parent(args.output, manifest, actual, supplied)
        if stage in {"planning", "report"}:
            from ..population.experiment_planning import verify_plan_revision
            qset, _, _ = _question_set(study, _parent(args.output, manifest, "questions"))
            verify_plan_revision(qset, artifacts["plan.json"])
            if stage == "planning" and artifacts["plan.json"]["previous_plan_id"] != _previous_plan(args, study):
                raise PopulationError("Reused plan has a different previous revision")
        return manifest
    inputs = [study["config_path"], *[source["path"] for source in study["sources"].values()]]
    inputs.extend(value for name in ("search", "observations", "evidence", "comparison", "questions", "plan")
                  if (value := getattr(args, name, None)) is not None)
    if Path(args.output).exists():
        raise PopulationError("Output exists; use --reuse to verify it or choose a new directory")
    check_output(args.output, inputs)
    parents = {}
    if args.command == "search":
        artifacts = _search(args, study)
    elif args.command == "observations":
        from ..population.observations import normalize_observations
        if args.previous:
            path = Path(args.previous).resolve()
            if not any(source["role"] in {"raw_observations", "observation_history"}
                       and source["format"] == "jsonl" and Path(source["path"]) == path
                       for source in study["sources"].values()):
                raise PopulationError("Previous observations must be a declared, hashed canonical history JSONL input")
        artifacts = normalize_observations(study)
    elif args.command == "evidence":
        from ..population.comparison_evidence import reconcile_evidence
        search, _ = load_stage(study, args.search, "search")
        observations, _ = load_stage(study, args.observations, "observation_normalization")
        artifacts = reconcile_evidence(study, search, observations)
        parents = {"search": args.search, "observation_normalization": args.observations}
    elif args.command == "compare":
        from ..population.comparison_evidence import compare_with_evidence
        search, _ = load_stage(study, args.search, "search")
        evidence, manifest = load_stage(study, args.evidence, "evidence_applicability")
        _same_parent(args.evidence, manifest, "search", args.search)
        observations, _ = load_stage(study,
            _parent(args.evidence, manifest, "observation_normalization"), "observation_normalization")
        artifacts = compare_with_evidence(study, search, evidence, observations)
        parents = {"search": args.search, "evidence_applicability": args.evidence}
    elif args.command == "questions":
        from ..population.comparison_questions import generate_questions
        records, _, _ = load_stage_records(study, args.comparison, "comparison")
        artifacts = generate_questions(study, records)
        parents = {"comparison": args.comparison}
    elif args.command == "plan":
        from ..population.experiment_planning import plan_experiments
        qset, _, _ = _question_set(study, args.questions, args.evidence)
        plan = plan_experiments(qset, previous_plan_id=_previous_plan(args, study))
        artifacts = {"plan.json": plan, "plan-records.jsonl": [plan["plan_record"]]}
        parents = {"questions": args.questions, "evidence_applicability": args.evidence}
    else:
        from ..population.comparison_reporting import render_report
        from ..population.experiment_planning import verify_plan_revision
        plan_artifacts, plan_manifest = load_stage(study, args.plan, "planning")
        _same_parent(args.plan, plan_manifest, "questions", args.questions)
        evidence_path = _parent(args.plan, plan_manifest, "evidence_applicability")
        qset, questions, question_manifest = _question_set(study, args.questions, evidence_path)
        _same_parent(args.questions, question_manifest, "comparison", args.comparison)
        plan = plan_artifacts["plan.json"]
        verify_plan_revision(qset, plan)
        comparison, comparison_manifest = load_stage(study, args.comparison, "comparison")
        search, _ = load_stage(study, _parent(args.comparison, comparison_manifest, "search"), "search")
        evidence, evidence_manifest = load_stage(study, evidence_path, "evidence_applicability")
        observations, _ = load_stage(study,
            _parent(evidence_path, evidence_manifest, "observation_normalization"), "observation_normalization")
        artifacts = render_report(study, search=search, observations=observations,
            evidence=evidence, comparison=comparison, questions=questions, plan=plan, question_set=qset)
        parents = {"comparison": args.comparison, "questions": args.questions, "planning": args.plan}
    return save_stage(study, stage, args.output, artifacts, parents=parents)


def _safe_reason(exc):
    if not isinstance(exc, PopulationError):
        return "Required local input could not be read" if isinstance(exc, (OSError, UnicodeError)) else \
               "A local artifact has missing or invalid fields"
    reason = str(exc)
    reason = re.sub(r"; unexpected fields \[.+\]", "; unexpected fields supplied", reason)
    reason = re.sub(r"(Invalid native [^:]+:).+", r"\1 [protected value]", reason)
    reason = re.sub(r"(artifact changed:).+", r"\1 [protected artifact]", reason)
    reason = re.sub(r"(?:/Users/|/private/|/tmp/|[A-Za-z]:\\)[^\s,;]+", "[protected path]", reason)
    return reason


def _retain_partial(args, study, exc, *, interrupted=False):
    if getattr(args, "command", None) == "benchmark" and getattr(args, "save_cases", None):
        archive = Path(args.save_cases).resolve()
        return getattr(exc, "case_archive_path", None) == str(archive) and (archive / "progress.json").is_file()
    if study is None or not getattr(args, "output", None) or Path(args.output).exists():
        return False
    audit = {"status": "interrupted" if interrupted else "partial", "command": args.command,
             "snapshot_run_id": study["snapshot_run_id"], "code": study["code"],
             "reason_codes": ["interrupted" if interrupted else "resource_limit"],
             "execution_audit": getattr(exc, "execution_audit", None)}
    write_bundle(args.output, {"stage-status.json": json_text(audit)},
        inputs=[study["config_path"], *[source["path"] for source in study["sources"].values()]])
    return True


def main(argv=None):
    parser = _Parser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    validate = commands.add_parser("validate", help="Validate local input declarations and comparison contracts")
    validate.add_argument("--config", required=True)
    computational = commands.add_parser("computational", help="Compare reference and sample sequences within one bounded region")
    computational.add_argument("--config", required=True)
    computational.add_argument("--output", required=True)
    search = _stage_arguments(commands, "search", "Search the declared donor/domain inputs")
    search.add_argument("--backend", choices=("local", "crisprme_plus"), help="Must match the configured profiles")
    search.add_argument("--native-output", help="Retain native commands, inputs and raw results separately")
    observations = _stage_arguments(commands, "observations", "Import pinned measurements and revision history")
    observations.add_argument("--previous", help="Declared earlier observation JSONL whose revisions must be retained")
    _stage_arguments(commands, "evidence", "Reconcile measured observations against active designs", ("search", "observations"))
    _stage_arguments(commands, "compare", "Compare design scope and measured evidence", ("search", "evidence"))
    _stage_arguments(commands, "questions", "Generate complete experimental requirements", ("comparison",))
    plan = _stage_arguments(commands, "plan", "Choose complete experiments under the configured budget", ("questions", "evidence"))
    plan.add_argument("--previous-plan", help="Earlier verified planning stage declared as plan_history")
    plan.add_argument("--previous-config", help="Pinned original configuration for the earlier plan")
    _stage_arguments(commands, "report", "Write the decision report and protected detailed exports", ("comparison", "questions", "plan"))
    benchmark = commands.add_parser("benchmark", help="Compare selectors using a frozen protocol and verified stages")
    benchmark.add_argument("--protocol", required=True)
    benchmark.add_argument("--output", "--output-dir", dest="output", required=True)
    benchmark.add_argument("--reuse", action="store_true")
    benchmark.add_argument("--save-cases", help="Preserve each checked case in a new, separate directory")
    args = study = None
    try:
        args = parser.parse_args(argv)
        if args.command == "computational":
            from ..population.computational import run
            manifest = run(args.config, args.output)
            print(json.dumps({key: manifest[key] for key in ("status", "run_id", "biological_status")}, sort_keys=True))
            return 0
        if args.command == "benchmark":
            from ..population.comparison_benchmark import run_benchmark
            manifest = run_benchmark(args.protocol, args.output, reuse=args.reuse, save_cases=args.save_cases)
            print(json.dumps({key: manifest[key] for key in ("status", "benchmark_id", "result_count")}, sort_keys=True))
            return 0 if manifest["status"] == "completed" else 4
        study = resolve_study(args.config)
        result = {"status": "completed", "command": args.command,
                  "snapshot_run_id": study["snapshot_run_id"],
                  "catalog_records": len(study["config"]["catalogs"]), "local_sources": len(study["sources"])}
        from ..population.comparison_domain import assess_domain
        domain = assess_domain(study["config"])
        if domain is not None and args.command == "validate":
            result["supported_domain"] = domain
        if args.command != "validate":
            manifest = _execute_stage(args, study)
            result.update(status=manifest["status"], stage_key=manifest["stage_key"], artifacts=len(manifest["outputs"]))
        print(json.dumps(result, sort_keys=True, allow_nan=False))
        return 0 if result["status"] == "completed" else 4
    except (ResourceLimit, KeyboardInterrupt) as exc:
        interrupted = isinstance(exc, KeyboardInterrupt)
        try:
            audit_retained = _retain_partial(args, study, exc, interrupted=interrupted)
        except (PopulationError, OSError, UnicodeError):
            audit_retained = False
        print(json.dumps({"status": "interrupted" if interrupted else "partial",
            "reason_code": "interrupted" if interrupted else "resource_limit",
            "partial_audit_write_succeeded": audit_retained,
            "reason": "Execution interrupted" if interrupted else _safe_reason(exc)}), file=sys.stderr)
        return 4
    except (ValueError, KeyError, TypeError, OSError, UnicodeError) as exc:
        from ..population.comparison_domain import DomainAbstention
        if isinstance(exc, DomainAbstention):
            # Rejected contig labels may contain protected identifiers or paths.
            coverage = exc.assessment["coverage_completeness"]
            coverage["declared_contig_count"] = len(coverage.pop("declared_contigs"))
            print(json.dumps({"status": "abstained", "reason_code": "unsupported_domain",
                "supported_domain": exc.assessment}, sort_keys=True), file=sys.stderr)
            return 2
        print(json.dumps({"status": "invalid_input", "reason_code": "invalid_input",
                          "reason": _safe_reason(exc)}), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
