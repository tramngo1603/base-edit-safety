"""Compare frozen selectors on one verified question and evidence catalog."""
from __future__ import annotations

from collections import Counter
from dataclasses import asdict
from datetime import datetime, timezone
from io import StringIO
from pathlib import Path
import csv
import hashlib
import math
import re
import time

from .common import (PopulationError, check_output, checksum, fields, integer,
                     json_text, jsonl_text, nonempty, parse_json, read_json, write_bundle,
                     write_file)
from .comparison_identity import content_id
from .comparison_provenance import load_stage_records, resolve_study, verify_sources
from .comparison_questions import prepare_question_set
from .experiment_planning import (SELECTORS, PlannerLimits, _evaluate_prepared,
                                  _plan_prepared, _planning_source_hash,
                                  _verify_prepared_revision)
from .provenance import code_identity


MAX_BENCHMARK_RUNS = 256
_REQUEST_FIELDS = {"schema_version", "policy_version", "study", "questions", "evidence",
                   "frozen_protocol", "budgets", "selectors", "random_seeds"}
_SCHEDULE_FIELDS = {"budgets", "selectors", "random_seeds"}


def _sha(value, label):
    if not isinstance(value, str) or re.fullmatch("[0-9a-f]{64}", value) is None:
        raise PopulationError(label + " must be a lowercase SHA-256 checksum")
    return value


def _array(value, label):
    if not isinstance(value, list):
        raise PopulationError(label + " must be an array")
    if len({json_text(item) for item in value}) != len(value):
        raise PopulationError(label + " must not contain duplicate values")
    return value


def _schedule(value, label):
    budgets = _array(value["budgets"], label + ".budgets")
    selectors = _array(value["selectors"], label + ".selectors")
    seeds = _array(value["random_seeds"], label + ".random_seeds")
    if not budgets or not selectors:
        raise PopulationError(label + " requires at least one budget and selector")
    for budget in budgets:
        integer(budget, label + ".budgets[]")
    for selector in selectors:
        if selector not in SELECTORS:
            raise PopulationError(label + ".selectors contains an unsupported selector")
    for seed in seeds:
        integer(seed, label + ".random_seeds[]")
    if bool(seeds) != ("random_bundle" in selectors):
        raise PopulationError(label + ".random_seeds must be nonempty exactly when random_bundle is selected")
    count = len(budgets) * (len(selectors) - ("random_bundle" in selectors) + len(seeds))
    if count > MAX_BENCHMARK_RUNS:
        raise PopulationError(label + " exceeds the 256-run benchmark limit")


def _reference(value, base, label, *, directory=False):
    fields(value, {"path", "sha256"}, label=label)
    nonempty(value["path"], label + ".path")
    _sha(value["sha256"], label + ".sha256")
    path = (base / value["path"]).resolve()
    source = path / "run-manifest.json" if directory else path
    if not source.is_file() or checksum(source) != value["sha256"]:
        raise PopulationError(label + " is unavailable or its checksum changed")
    return path


def _ancestry(directory, manifest):
    """Retain every verified dependency, including shared ancestor directories."""
    result = {Path(directory).resolve() / "run-manifest.json": manifest}
    pending = list(result.items())
    while pending:
        path, current = pending.pop()
        for reference in current["parents"].values():
            parent = (path.parent / reference["path"]).resolve()
            if checksum(parent) != reference["sha256"]:
                raise PopulationError("Benchmark ancestor changed during loading")
            if parent not in result:
                result[parent] = read_json(parent)
                pending.append((parent, result[parent]))
    return result


def _load_request(path):
    started = time.monotonic()
    path = Path(path).resolve()
    raw_request = path.read_bytes()
    request = parse_json(raw_request.decode("utf-8"))
    fields(request, _REQUEST_FIELDS, label="benchmark request")
    if request["schema_version"] != "1":
        raise PopulationError("benchmark request.schema_version must be 1")
    nonempty(request["policy_version"], "benchmark request.policy_version")
    _schedule(request, "benchmark request")
    paths = {name: _reference(request[name], path.parent, "benchmark request." + name,
                             directory=name in {"questions", "evidence"})
             for name in ("study", "questions", "evidence", "frozen_protocol")}
    raw_frozen = paths["frozen_protocol"].read_bytes()
    if hashlib.sha256(raw_frozen).hexdigest() != request["frozen_protocol"]["sha256"]:
        raise PopulationError("Frozen protocol changed during loading")
    frozen = parse_json(raw_frozen.decode("utf-8"))
    required = {"schema_version", "status", "policy_version", "code_snapshot",
                "study_sha256", "selector_benchmark"}
    if not isinstance(frozen, dict) or required - frozen.keys():
        raise PopulationError("Frozen benchmark protocol is missing its selector contract")
    if frozen["schema_version"] != "1" or frozen["status"] != "frozen":
        raise PopulationError("Benchmark protocol must have schema_version 1 and status frozen")
    fields(frozen["selector_benchmark"], _SCHEDULE_FIELDS, label="frozen protocol.selector_benchmark")
    _schedule(frozen["selector_benchmark"], "frozen protocol.selector_benchmark")
    if frozen["selector_benchmark"] != {key: request[key] for key in _SCHEDULE_FIELDS} \
            or frozen["policy_version"] != request["policy_version"]:
        raise PopulationError("Benchmark request differs from the frozen selector policy")
    if _sha(frozen["study_sha256"], "frozen protocol.study_sha256") != request["study"]["sha256"]:
        raise PopulationError("Benchmark study differs from the frozen study checksum")
    study = resolve_study(paths["study"])
    if _sha(frozen["code_snapshot"], "frozen protocol.code_snapshot") != study["code"]["source_sha256"]:
        raise PopulationError("Benchmark implementation differs from the frozen code snapshot")
    qrows, _, questions = load_stage_records(study, paths["questions"], "questions")
    erows, _, evidence = load_stage_records(study, paths["evidence"], "evidence_applicability")
    ancestry = _ancestry(paths["questions"], questions)
    if not any(row["stage"] == "evidence_applicability" and row["stage_key"] == evidence["stage_key"]
               and row["outputs"] == evidence["outputs"] for row in ancestry.values()):
        raise PopulationError("Benchmark evidence must be the evidence used to generate the questions")
    ancestry.update(_ancestry(paths["evidence"], evidence))
    records = {}
    for row in [*qrows, *erows]:
        if row["record_id"] in records and records[row["record_id"]] != row:
            raise PopulationError("Benchmark ancestry contains conflicting canonical records")
        records[row["record_id"]] = row
    qset = prepare_question_set(study, list(records.values()))
    # Preparation owns its copies. Release the loaded rows before hashing them.
    del qrows, erows, records
    planning_source_sha256 = _planning_source_hash(qset)
    bindings = {"request_sha256": hashlib.sha256(raw_request).hexdigest(), "frozen_protocol_sha256": request["frozen_protocol"]["sha256"],
                "study_sha256": request["study"]["sha256"], "snapshot_run_id": study["snapshot_run_id"],
                "code": study["code"], "questions": {"sha256": request["questions"]["sha256"],
                    "stage_key": questions["stage_key"], "outputs": questions["outputs"]},
                "evidence": {"sha256": request["evidence"]["sha256"],
                    "stage_key": evidence["stage_key"], "outputs": evidence["outputs"]}}
    benchmark_id = content_id("comparison_benchmark", {"request": request, "bindings": bindings})
    inputs = [path, *paths.values(), *[Path(row["path"]) for row in study["sources"].values()],
              *[p.parent for p in ancestry]]
    hashes = {path: bindings["request_sha256"], paths["study"]: request["study"]["sha256"],
              paths["frozen_protocol"]: request["frozen_protocol"]["sha256"]}
    for manifest_path, manifest in ancestry.items():
        hashes[manifest_path] = checksum(manifest_path)
        hashes.update({manifest_path.parent / name: sha for name, sha in manifest["outputs"].items()})
    return {"request": request, "frozen": frozen, "study": study, "qset": qset,
            "bindings": bindings, "benchmark_id": benchmark_id, "inputs": inputs, "hashes": hashes,
            "planning_source_sha256": planning_source_sha256,
            "request_setup_seconds": time.monotonic() - started}


def _runs(context):
    request = context["request"]
    seed = context["study"]["config"]["policies"]["planning"]["seed"]
    return [(budget, selector, chosen) for budget in request["budgets"] for selector in request["selectors"]
            for chosen in (request["random_seeds"] if selector == "random_bundle" else [seed])]


def _verify_current(context):
    verify_sources(context["study"])
    if code_identity() != context["study"]["code"]:
        raise PopulationError("Benchmark implementation changed during execution")
    for path, expected in context["hashes"].items():
        if not path.is_file() or checksum(path) != expected:
            raise PopulationError("Benchmark input changed during execution")


def _evaluate(qset, plan):
    verify_sources(qset.study)
    checked = _evaluate_prepared(qset, plan["selected_action_ids"], plan["resource_assignments"], budget=plan["budget"])
    if checked["requirements"] != plan["requirements"] or checked["meets_request"] != plan["meets_request"] \
            or checked["unmet_hard_required_question_ids"] != plan["unmet_hard_required_question_ids"] \
            or {key: checked[key] for key in plan["accounting"]} != plan["accounting"]:
        raise PopulationError("Benchmark plan disagrees with the common requirement and cost evaluator")
    for key in ("runtime_seconds", "process_peak_rss_bytes"):
        value = plan[key]
        if type(value) not in (int, float) or not math.isfinite(value) or value < 0:
            raise PopulationError("Benchmark plan has an invalid runtime measurement")
    return checked


def _seed_fairness(plans):
    seeds = {}
    for plan in plans:
        if plan["selector"] == "exact_tiny":
            continue
        previous = seeds.setdefault(plan["budget"], plan["mandatory_seed"])
        if previous != plan["mandatory_seed"]:
            raise PopulationError("Benchmark selectors did not use the same mandatory-question seed")


def _fraction(numerator, denominator):
    return numerator / denominator if denominator else None


def _metrics(qset, plan, exact):
    satisfied = set(plan["requirements"]["satisfied_slot_ids"])
    engaged = {(qset.slots[sid]["question_id"], qset.slots[sid]["donor_id"]) for sid in satisfied}
    engaged_slots = [row for row in qset.slots.values() if (row["question_id"], row["donor_id"]) in engaged]
    intended = {row["slot_id"] for row in engaged_slots if row["control_role"] == "intended_activity"}
    controls = {row["slot_id"] for row in engaged_slots if row["control_role"] not in {"measurement", "intended_activity"}}
    questions = plan["questions"]
    variant = {qid for qid, question in qset.questions.items() if question["template"] in {"variant_association", "design_by_variant"}}
    evidence_states = Counter(state["evidence_status"] for state in questions.values())
    donor_statuses = list(qset.donor_eligibility.values())
    unresolved = sum(row["membership"] is None for row in donor_statuses)
    ineligible = sum(not row["eligible"] for row in donor_statuses)
    reasons = Counter(reason for row in donor_statuses for reason in row["reason_codes"])
    action_reasons = Counter(reason for row in qset.eligibility.values() for reason in row["reason_codes"])
    search = plan["searches"]
    optimum = exact.get(plan["budget"])
    gap = optimum["weighted_completion"] - plan["weighted_completion"] if optimum is not None else None
    if gap is not None and gap < 0:
        raise PopulationError("Benchmark plan exceeds a claimed exact optimum")
    return {"plan_id": plan["plan_id"], "budget": plan["budget"], "selector": plan["selector"],
        "random_seed": plan["random_seed"], "status": plan["status"], "meets_request": plan["meets_request"],
        "weighted_completion": plan["weighted_completion"], "unweighted_completion": plan["unweighted_completion"],
        "question_count": len(questions), "total_question_weight": sum(row["priority_weight"] for row in qset.questions.values()),
        "existing_evidence_weight": plan["existing_evidence_weight"], "newly_planned_weight": plan["newly_planned_weight"],
        "catalog_slot_count": len(qset.slots), "filled_catalog_slot_count": len(satisfied),
        "engaged_question_donor_bundles": len(engaged), "engaged_bundle_slot_count": len(engaged_slots),
        "engaged_intended_slots": len(intended), "filled_engaged_intended_slots": len(intended & satisfied),
        "engaged_control_slots": len(controls), "filled_engaged_control_slots": len(controls & satisfied),
        "engaged_intended_completeness": _fraction(len(intended & satisfied), len(intended)),
        "engaged_control_completeness": _fraction(len(controls & satisfied), len(controls)),
        "total_cost": plan["accounting"]["total_cost"], "cost_unit": plan["accounting"]["cost_unit"],
        "cost_scale": plan["accounting"]["cost_scale"], "active_resource_count": len(plan["accounting"]["active_resource_ids"]),
        "mandatory_seed_cost": plan["mandatory_seed"]["total_cost"],
        "variant_question_count": len(variant),
        "completed_variant_questions": sum(questions[qid]["complete"] for qid in variant),
        "unresolved_variant_questions": sum(not questions[qid]["complete"] for qid in variant),
        "covered_configuration_count": len(plan["covered_configuration_ids"]),
        "coverage_definition": plan["coverage_definition"], "evidence_state_counts": dict(sorted(evidence_states.items())),
        "unevaluable_question_fraction": _fraction(len(questions) - evidence_states["evaluable"], len(questions)),
        "donor_question_pair_count": len(donor_statuses), "unresolved_membership_pairs": unresolved,
        "unresolved_membership_fraction": _fraction(unresolved, len(donor_statuses)),
        "ineligible_donor_question_pairs": ineligible, "ineligible_donor_question_fraction": _fraction(ineligible, len(donor_statuses)),
        "donor_eligibility_reason_counts": dict(sorted(reasons.items())),
        "ineligible_action_count": sum(not row["eligible"] for row in qset.eligibility.values()),
        "action_eligibility_reason_counts": dict(sorted(action_reasons.items())),
        "construction_limited_question_count": len(plan["construction_limited_question_ids"]),
        "search_expansions": sum(row["expansions"] for row in search),
        "beam_pruned": sum(row["beam_pruned"] for row in search),
        "alternative_pruned": sum(row["alternative_pruned"] for row in search),
        "search_status_counts": dict(sorted(Counter(row["status"] for row in search).items())),
        "searches_exhaustive": all(row["exhaustive"] for row in search),
        "expansion_limit_reached": any(row["expansion_limit_reached"] for row in search),
        "limits": plan["limits"], "runtime_seconds": plan["runtime_seconds"],
        "process_peak_rss_bytes": plan["process_peak_rss_bytes"], "memory_measurement": plan["memory_measurement"],
        "optimality_established": plan["optimality_established"], "exact_weighted_optimality_gap": gap,
        "exact_gap_reason": None if gap is not None else "No optimal exact solve at this budget in the frozen schedule"}


def _arm_rows(qset, plan):
    rows = []
    for qid, question in sorted(qset.questions.items()):
        state = plan["requirements"]["questions"][qid]
        selected = set(state["selected_donor_ids"])
        for arm in question["arms"]:
            donors = [row for row in state["donors"] if row["arm_id"] == arm["arm_id"]]
            rows.append({"plan_id": plan["plan_id"], "question_id": qid, "arm_id": arm["arm_id"],
                "required_independent_donors": arm["required_independent_donors"],
                "eligible_donor_count": sum(row["eligible"] for row in donors),
                "complete_independence_group_count": len({row["independence_group_id"] for row in donors
                    if row["complete"] and row["independence_group_id"] is not None}),
                "selected_independent_donor_count": sum(row["donor_id"] in selected for row in donors),
                "question_complete": state["complete"], "evidence_status": state["evidence_status"]})
    return rows


def _tsv(rows, columns):
    output = StringIO()
    writer = csv.DictWriter(output, fieldnames=columns, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({key: json_text(row[key]).strip().replace("\n", "") if isinstance(row[key], (dict, list))
                         else "" if row[key] is None else row[key] for key in columns})
    return output.getvalue()


def _files(context, plans):
    qset = context["qset"]
    exact = {plan["budget"]: plan for plan in plans if plan["selector"] == "exact_tiny" and plan["optimality_established"]}
    summaries = [_metrics(qset, plan, exact) for plan in plans]
    arms = [row for plan in plans for row in _arm_rows(qset, plan)]
    slots = [{"plan_id": plan["plan_id"], "slot_id": sid, "question_id": row["question_id"],
              "donor_id": row["donor_id"], "arm_id": row["arm_id"], "control_role": row["control_role"],
              "assignment": plan["requirements"]["slot_assignments"].get(sid),
              "filled": sid in plan["requirements"]["slot_assignments"]}
             for plan in plans for sid, row in sorted(qset.slots.items())]
    files = {f"plan-{i:04d}.json": json_text(plan) for i, plan in enumerate(plans, 1)}
    files.update({"benchmark-request.json": json_text(context["request"]),
        "frozen-protocol.json": json_text(context["frozen"]),
        "plan-records.jsonl": jsonl_text([plan["plan_record"] for plan in plans]),
        "benchmark-summary.json": json_text(summaries),
        "benchmark-summary.tsv": _tsv(summaries, list(summaries[0])),
        "question-arms.tsv": _tsv(arms, ["plan_id", "question_id", "arm_id", "required_independent_donors",
            "eligible_donor_count", "complete_independence_group_count", "selected_independent_donor_count",
            "question_complete", "evidence_status"]),
        "requirement-slots.tsv": _tsv(slots, ["plan_id", "slot_id", "question_id", "donor_id", "arm_id",
            "control_role", "assignment", "filled"])})
    return files


def _manifest(context, plans, files):
    return {"artifact_type": "comparison_benchmark_manifest", "schema_version": "1", "status": "completed",
        "benchmark_id": context["benchmark_id"], "bindings": context["bindings"], "result_count": len(plans),
        "request_setup_seconds": context["request_setup_seconds"],
        "case_runtime_scope": "Planning on shared prepared inputs; request setup is reported separately",
        "mandatory_seed_comparison": "identical_within_budget_for_all_nonexact_selectors",
        "results": [{"budget": plan["budget"], "selector": plan["selector"], "random_seed": plan["random_seed"],
                     "plan_id": plan["plan_id"], "file": f"plan-{i:04d}.json"} for i, plan in enumerate(plans, 1)],
        "outputs": {name: hashlib.sha256(value.encode("utf-8")).hexdigest() for name, value in files.items()}}


def _check_case_policy(context, plan, expected):
    required = {"budget", "selector", "random_seed", "previous_plan_id", "limits",
                "runtime_seconds", "process_peak_rss_bytes", "memory_measurement"}
    if not isinstance(plan, dict) or required - plan.keys() or \
            (plan.get("budget"), plan.get("selector"), plan.get("random_seed")) != expected:
        raise PopulationError("Benchmark plan does not match the frozen schedule")
    limits = asdict(PlannerLimits.from_policy(context["study"]["config"]["policies"]["planning"]))
    if plan["limits"] != limits or plan["previous_plan_id"] is not None:
        raise PopulationError("Benchmark plan changed the frozen planning limits or revision policy")


def _reuse(context, directory):
    manifest = read_json(directory / "run-manifest.json")
    setup = manifest.get("request_setup_seconds")
    if type(setup) not in (int, float) or not math.isfinite(setup) or setup < 0:
        raise PopulationError("Benchmark has an invalid setup runtime measurement")
    context = {**context, "request_setup_seconds": setup}
    plans = []
    for i, (budget, selector, seed) in enumerate(_runs(context), 1):
        path = directory / f"plan-{i:04d}.json"
        if not path.is_file():
            raise PopulationError("Benchmark result is incomplete")
        plan = read_json(path)
        _check_case_policy(context, plan, (budget, selector, seed))
        verify_sources(context["study"])
        _verify_prepared_revision(context["qset"], plan, source_sha256=context["planning_source_sha256"])
        _evaluate(context["qset"], plan)
        plans.append(plan)
    _seed_fairness(plans)
    files = _files(context, plans)
    expected = _manifest(context, plans, files)
    if manifest != expected:
        raise PopulationError("Benchmark manifest does not match the verified request and results")
    if {path.name for path in directory.iterdir()} != set(files) | {"run-manifest.json"}:
        raise PopulationError("Benchmark output contains missing or unexpected files")
    for name, text in files.items():
        path = directory / name
        if not path.is_file() or path.read_bytes() != text.encode("utf-8"):
            raise PopulationError("Benchmark artifact changed or is unavailable")
    _verify_current(context)
    return manifest


class _CaseOutput:
    """Keep completed cases immutable while the full schedule is still running."""

    def __init__(self, directory, context, output):
        self.directory = Path(directory).resolve()
        self.context = context
        self.progress = {"artifact_type": "comparison_benchmark_progress", "schema_version": "1",
                         "benchmark_id": context["benchmark_id"], "status": "incomplete", "cases": [
            {"ordinal": index, "budget": budget, "selector": selector, "random_seed": seed,
             "status": "not_started"}
            for index, (budget, selector, seed) in enumerate(_runs(context), 1)]}
        write_bundle(directory, {
            "schedule.json": json_text({"artifact_type": "comparison_benchmark_case_archive", "schema_version": "1",
                "benchmark_id": context["benchmark_id"], "bindings": context["bindings"],
                "request_setup_seconds": context["request_setup_seconds"],
                "purpose": "Preserve cases; this archive is not a completed benchmark",
                "schedule": _runs(context)}),
            "progress.json": json_text(self.progress)}, inputs=[*context["inputs"], output])

    def _progress(self):
        self.progress["updated_at"] = datetime.now(timezone.utc).isoformat()
        write_file(self.directory / "progress.json", json_text(self.progress), overwrite=True)

    def start(self, index):
        self.progress["cases"][index - 1].update(status="running", started_at=datetime.now(timezone.utc).isoformat())
        self._progress()

    def save(self, index, plan):
        text = json_text(plan)
        manifest = {"artifact_type": "comparison_benchmark_case", "schema_version": "1", "status": "completed",
                    "benchmark_id": self.context["benchmark_id"], "bindings": self.context["bindings"],
                    "ordinal": index, "budget": plan["budget"], "selector": plan["selector"],
                    "random_seed": plan["random_seed"], "plan_id": plan["plan_id"],
                    "outputs": {"plan.json": hashlib.sha256(text.encode("utf-8")).hexdigest()}}
        directory = self.directory / f"case-{index:04d}"
        write_bundle(directory, {"plan.json": text, "run-manifest.json": json_text(manifest)},
                     inputs=self.context["inputs"])
        self.progress["cases"][index - 1].update(status="completed", plan_id=plan["plan_id"],
            completed_at=datetime.now(timezone.utc).isoformat(), manifest_sha256=checksum(directory / "run-manifest.json"))
        self._progress()

    def failed(self, exc):
        interrupted = isinstance(exc, (KeyboardInterrupt, SystemExit))
        self.progress["status"] = "interrupted" if interrupted else "failed"
        self.progress["reason_code"] = type(exc).__name__
        for case in self.progress["cases"]:
            if case["status"] == "running":
                case["status"] = "incomplete" if interrupted else "failed"
                case["reason_code"] = type(exc).__name__
        self._progress()

    def complete(self, output):
        self.progress["status"] = "schedule_completed"
        self.progress["completed_output_manifest_sha256"] = checksum(output / "run-manifest.json")
        self._progress()


def run_benchmark(protocol_path, output_dir, *, reuse=False, save_cases=None):
    """Execute a pinned selector request, or fully verify a completed result.

    The referenced scientific protocol must have been frozen before evaluation.
    A late request cannot establish when that earlier freeze occurred.
    """
    if reuse and save_cases is not None:
        raise PopulationError("Completed benchmark reuse does not create a case archive")
    context = _load_request(protocol_path)
    output = Path(output_dir).resolve()
    check_output(output_dir, context["inputs"], overwrite=True)
    if reuse and not output.exists():
        raise PopulationError("Benchmark reuse requires an existing completed output")
    if not reuse and output.exists():
        raise PopulationError("Benchmark output exists; use --reuse or choose a new directory")
    if output.exists():
        if not output.is_dir():
            raise PopulationError("Benchmark output must be a directory")
        return _reuse(context, output)
    cases = _CaseOutput(save_cases, context, output) if save_cases is not None else None
    plans = []
    try:
        for index, (budget, selector, seed) in enumerate(_runs(context), 1):
            if cases:
                cases.start(index)
            verify_sources(context["study"])
            plan = _plan_prepared(context["qset"], budget=budget, selector=selector, random_seed=seed,
                                  source_sha256=context["planning_source_sha256"])
            _check_case_policy(context, plan, (budget, selector, seed))
            _evaluate(context["qset"], plan)
            plans.append(plan)
            _seed_fairness(plans)
            if cases:
                _verify_current(context)
                cases.save(index, plan)
        files = _files(context, plans)
        manifest = _manifest(context, plans, files)
        files["run-manifest.json"] = json_text(manifest)
        _verify_current(context)
        write_bundle(output_dir, files, inputs=context["inputs"], overwrite=False)
        if cases:
            cases.complete(output)
        return manifest
    except BaseException as exc:
        if cases:
            exc.case_archive_path = str(cases.directory)
            try:
                cases.failed(exc)
            except (OSError, ValueError) as status_error:
                add_note = getattr(exc, "add_note", None)
                if callable(add_note):
                    add_note("Could not update case progress: " + str(status_error))
        raise
