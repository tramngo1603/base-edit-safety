"""One owned catalog preserves public validation and the frozen case policy."""
from copy import deepcopy
from dataclasses import asdict

import pytest

from base_edit_safety.population import comparison_benchmark as benchmark
from base_edit_safety.population import experiment_planning as planning
from base_edit_safety.population.common import PopulationError, json_text
from planning_cases import full_fixture, write_case
from planning_oracle import TinyOracle
from test_comparison_benchmark import _request


def deterministic(plan):
    return {key: value for key, value in plan.items()
            if key not in {"runtime_seconds", "process_peak_rss_bytes"}}


def mutable_ids(value):
    found = set()

    def visit(item):
        if isinstance(item, (dict, list, set)):
            if id(item) in found:
                return
            found.add(id(item))
        if isinstance(item, dict):
            for child in item.values():
                visit(child)
        elif isinstance(item, (list, tuple, set)):
            for child in item:
                visit(child)
    visit(value)
    return found


def scientific_state(qset):
    return (qset.study, qset.records, qset.index, qset.eligibility, qset.options,
            qset.donor_eligibility, qset.observation_tokens, qset.evidence_rejections)


@pytest.mark.parametrize("field", [*planning.PlannerLimits.__dataclass_fields__, "previous_plan_id"])
def test_rehashed_override_plan_cannot_change_frozen_benchmark_policy(tmp_path, field):
    path, qset, _ = _request(tmp_path / "case", budgets=(18,), selectors=("exact_tiny",))
    output = tmp_path / "results"
    benchmark.run_benchmark(path, output)
    context = benchmark._load_request(path)
    limits = asdict(planning.PlannerLimits.from_policy(qset.study["config"]["policies"]["planning"]))
    previous = None
    if field == "previous_plan_id":
        previous = "comparison_experiment_plan:v1:" + "0" * 64
    else:
        limits[field] = 1 if limits[field] != 1 else 2
    changed = planning.plan_experiments(qset, budget=18, selector="exact_tiny", limits=limits,
                                        previous_plan_id=previous)
    files = benchmark._files(context, [changed])
    manifest = benchmark._manifest(context, [changed], files)
    for name, text in {**files, "run-manifest.json": json_text(manifest)}.items():
        (output / name).write_text(text)
    with pytest.raises(PopulationError, match="frozen planning limits or revision policy"):
        benchmark.run_benchmark(path, output, reuse=True)


def test_benchmark_prepares_and_hashes_once_but_replays_every_case(tmp_path, monkeypatch):
    path, _, _ = _request(tmp_path / "case", budgets=(17, 18), selectors=planning.SELECTORS, seeds=(7, 19, 43))
    counts = {"prepare": 0, "hash": 0, "replay": 0, "sources": 0}
    for name, counter in (("prepare_question_set", "prepare"), ("_planning_source_hash", "hash"),
                          ("_verify_prepared_revision", "replay"), ("verify_sources", "sources")):
        original = getattr(benchmark, name)

        def counted(*args, _original=original, _counter=counter, **kwargs):
            counts[_counter] += 1
            return _original(*args, **kwargs)
        monkeypatch.setattr(benchmark, name, counted)
    output = tmp_path / "results"
    manifest = benchmark.run_benchmark(path, output)
    assert manifest["result_count"] == 16
    assert counts["prepare"] == counts["hash"] == 1
    assert counts["sources"] >= 2 * manifest["result_count"]
    assert benchmark.run_benchmark(path, output, reuse=True) == manifest
    assert counts["prepare"] == counts["hash"] == 2
    assert counts["replay"] == manifest["result_count"]


def test_prepared_cases_match_public_plans_and_do_not_share_result_state(tmp_path):
    qset = full_fixture(tmp_path / "source")["question_set"]
    before = deepcopy(scientific_state(qset))
    source_hash = planning._planning_source_hash(qset)
    oracle = TinyOracle(qset.study["config"], qset.records)
    first = None
    for budget, selector, seed in [(18, "complete_bundle_greedy", 7), (17, "complete_bundle_greedy", 7),
                                   *[(18, selector, seed) for selector in planning.SELECTORS
                                     for seed in ((7, 19, 43) if selector == "random_bundle" else (7,))],
                                   (18, "complete_bundle_greedy", 7)]:
        plan = planning._plan_prepared(qset, budget=budget, selector=selector, random_seed=seed,
                                      source_sha256=source_hash)
        public = planning.plan_experiments(qset, budget=budget, selector=selector, random_seed=seed)
        assert deterministic(plan) == deterministic(public)
        checked = planning._evaluate_prepared(qset, plan["selected_action_ids"], plan["resource_assignments"], budget=budget)
        assert checked["requirements"] == plan["requirements"]
        independent = oracle.evaluate(plan["selected_action_ids"], {
            (row["action_id"], row["assignment_key"]): row["resource_id"]
            for row in plan["resource_assignments"]})
        assert independent["weighted_completion"] == plan["weighted_completion"]
        assert independent["total_cost"] == plan["accounting"]["total_cost"]
        assert not mutable_ids(scientific_state(qset)) & mutable_ids((plan, checked))
        assert scientific_state(qset) == before
        if first is None:
            first = deterministic(plan)
    assert deterministic(plan) == first


def test_returned_gap_and_ineligibility_lists_cannot_change_the_owned_catalog(tmp_path):
    fixture = full_fixture(tmp_path / "source")
    rows = deepcopy(fixture["question_set"].records)
    question = next(row for row in rows if row["record_type"] == "comparison_question")
    question.update(construction_status="excluded_by_declared_limit", requirement_slot_ids=[],
                    planning_status="blocked", remaining_gaps=["requirement_construction_shortlist_exclusion"])
    qset = write_case(fixture, tmp_path / "excluded", [r for r in rows if r["record_type"] != "comparison_requirement_slot"])
    before = deepcopy(scientific_state(qset))
    plan = planning._plan_prepared(qset)
    next(iter(plan["requirements"]["questions"].values()))["remaining_gaps"].clear()
    assert scientific_state(qset) == before
    aid = next(iter(qset.actions))
    qset.eligibility[aid] = {"eligible": False, "reason_codes": ["source_unavailable"]}
    checked = planning._evaluate_prepared(qset, [aid])
    problem = next(row for row in checked["problems"] if row["code"] == "ineligible_action")
    problem["reason_codes"].clear()
    assert qset.eligibility[aid]["reason_codes"] == ["source_unavailable"]


def test_public_boundaries_rebuild_forged_indexes_and_reject_changed_records(tmp_path):
    qset = full_fixture(tmp_path / "source")["question_set"]
    original = planning.plan_experiments(qset, budget=18)
    forged = deepcopy(qset)
    forged.actions.clear()
    forged.eligibility.clear()
    forged.options.clear()
    assert deterministic(planning.plan_experiments(forged, budget=18)) == deterministic(original)
    assert planning.evaluate_plan(forged, original["selected_action_ids"], budget=18)["meets_request"]
    assert deterministic(planning.verify_plan_revision(forged, original)) == deterministic(original)
    forged.records[0]["record_id"] = "changed-without-a-canonical-identity"
    for call in (lambda: planning.plan_experiments(forged), lambda: planning.evaluate_plan(forged),
                 lambda: planning.verify_plan_revision(forged, original)):
        with pytest.raises(PopulationError):
            call()
