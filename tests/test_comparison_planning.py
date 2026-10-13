"""Complete experimental bundles, exact instance costs and bounded-search claims."""
from copy import deepcopy
import json
from types import SimpleNamespace

import pytest

from base_edit_safety.population.common import PopulationError
from base_edit_safety.population.comparison_contracts import validate_catalog, validate_record
from base_edit_safety.population.experiment_planning import (
    PlannerLimits, _Search, _Selection, _account, _bounded, _construct, _exact,
    _plan_record, _select_optional, _selection,
    SELECTORS, evaluate_plan, plan_experiments, verify_plan_revision, write_plan_revision,
)
from planning_cases import full_fixture, joint_seed_case, alternatives_case, optional_trap_case, split_design_questions, write_case
from planning_oracle import TinyOracle


def accounting_fixture():
    compatibility = {field: [] for field in (
        "donor_ids", "design_ids", "context_ids", "batch_ids", "assay_capability_ids",
        "timepoints", "endpoint_ids", "locus_ids")}
    resources = {}
    for identifier, cost, capacity in (("panel-1", 3, 2), ("panel-2", 3, 2), ("setup", 2, None)):
        resources[identifier] = {"resource_id": identifier, "resource_type": "assay_library",
            "cost": cost, "cost_unit": "synthetic-units", "cost_scale": 1, "sharing_key": "declared-panel",
            "capacities": [{"dimension": "libraries", "capacity": {"kind": "finite", "quantity": capacity}
                             if capacity is not None else {"kind": "unbounded", "justification": "Declared setup"}}],
            "activation_rule": "on_use", "compatibility": deepcopy(compatibility)}
    resources["setup"]["activation_rule"] = "nonempty_plan"
    actions = {}
    for aid in ("a", "b", "c"):
        actions[aid] = {"action_id": aid, "operation_type": "preparation", "donor_id": "donor", "design_id": {"value": "design"},
            "context_id": "context", "batch_id": "batch", "assay_capability_id": "assay",
            "endpoint_ids": ["endpoint"], "locus_ids": ["locus"], "dependency_action_ids": [],
            "resource_demands": [{"eligible_resource_ids": ["panel-1", "panel-2"], "dimension": "libraries",
                                  "quantity": 1, "assignment_key": "measurement"}]}
    return SimpleNamespace(actions=actions, resources=resources,
        eligibility={aid: {"eligible": True, "reason_codes": []} for aid in actions}, index={},
        study={"policies": {"planning": {"cost_unit": "synthetic-units"}}})


def choose(actions, *instances):
    return _selection(actions, [{"action_id": aid, "assignment_key": "measurement", "resource_id": rid}
                                 for aid, rid in zip(actions, instances)])


def test_shared_instance_is_charged_once_and_overflow_requires_another_paid_instance():
    qset = accounting_fixture()
    shared = _account(qset, choose("ab", "panel-1", "panel-1"), 5)
    assert shared["feasible"] and shared["total_cost"] == 5
    overflow = _account(qset, choose("abc", "panel-1", "panel-1", "panel-1"), 100)
    assert not overflow["feasible"]
    assert "resource_capacity_exceeded" in {row["code"] for row in overflow["problems"]}
    two = _account(qset, choose("abc", "panel-1", "panel-1", "panel-2"), 8)
    assert two["feasible"] and two["total_cost"] == 8
    assert not _account(qset, choose("abc", "panel-1", "panel-1", "panel-2"), 7)["feasible"]


def test_assigned_instance_pays_its_own_cost_and_never_borrows_other_capacity():
    qset = accounting_fixture()
    qset.resources["panel-2"]["cost"] = 11
    expensive = _account(qset, choose("a", "panel-2"), 13)
    assert expensive["feasible"] and expensive["total_cost"] == 13
    assert not _account(qset, choose("a", "panel-2"), 5)["feasible"]
    qset.resources["panel-1"]["capacities"][0]["capacity"]["quantity"] = 0
    assert not _account(qset, choose("a", "panel-1"), 100)["feasible"]


def test_empty_plan_overhead_and_positive_zero_cost_use_have_exact_activation_rules():
    qset = accounting_fixture()
    assert _account(qset, _selection(()), 0)["total_cost"] == 0
    qset.resources["setup"]["activation_rule"] = "always"
    assert _account(qset, _selection(()), 2)["total_cost"] == 2
    assert not _account(qset, _selection(()), 1)["feasible"]
    qset.resources["setup"]["cost"] = 0
    qset.resources["panel-1"]["cost"] = 0
    assert _account(qset, choose("a", "panel-1"), 0)["feasible"]


def test_resource_names_prices_and_membership_flags_cannot_make_incompatible_sharing_valid():
    qset = accounting_fixture()
    qset.resources["panel-1"]["compatibility"]["donor_ids"] = ["different-donor"]
    bad = _account(qset, choose("a", "panel-1"), 100)
    assert not bad["feasible"]
    assert "incompatible_resource" in {row["code"] for row in bad["problems"]}
    assert _account(qset, choose("a", "panel-2"), 100)["feasible"]


def test_alternative_assignment_is_required_once_and_only_from_its_declared_catalog():
    qset = accounting_fixture()
    assert not _account(qset, _selection("a"), 100)["feasible"]
    assert not _account(qset, choose("a", "unlisted-instance"), 100)["feasible"]
    with pytest.raises(PopulationError, match="exactly one"):
        _selection("a", [{"action_id": "a", "assignment_key": "measurement", "resource_id": rid}
                          for rid in ("panel-1", "panel-2")])
    assert not _account(qset, choose("", "panel-1"), 1)["total_cost"]
    extras = _selection((), [{"action_id": "a", "assignment_key": "measurement", "resource_id": "panel-1"}])
    assert not _account(qset, extras, 100)["feasible"]


def test_multiple_dimensions_of_one_assignment_cannot_use_different_instances():
    qset = accounting_fixture()
    for resource in (qset.resources["panel-1"], qset.resources["panel-2"]):
        resource["capacities"].append({"dimension": "targets", "capacity": {"kind": "finite", "quantity": 2}})
    for action in qset.actions.values():
        action["resource_demands"].append({"eligible_resource_ids": ["panel-1", "panel-2"],
            "dimension": "targets", "quantity": 2, "assignment_key": "measurement"})
    assert _account(qset, choose("a", "panel-1"), 5)["feasible"]
    assert not _account(qset, choose("ab", "panel-1", "panel-1"), 8)["feasible"]
    assert _account(qset, choose("ab", "panel-1", "panel-2"), 8)["feasible"]


def test_missing_prerequisites_and_ineligible_actions_do_not_receive_a_valid_account():
    qset = accounting_fixture()
    qset.actions["b"]["dependency_action_ids"] = ["a"]
    assert not _account(qset, choose("b", "panel-1"), 100)["feasible"]
    assert _account(qset, choose("ab", "panel-1", "panel-1"), 100)["feasible"]
    qset.eligibility["a"] = {"eligible": False, "reason_codes": ["unresolved_source_sequence"]}
    assert not _account(qset, choose("ab", "panel-1", "panel-1"), 100)["feasible"]


def test_cost_scales_cannot_be_summed_without_conversion():
    qset = accounting_fixture()
    qset.resources["panel-2"]["cost_scale"] = 100
    assert not _account(qset, choose("a", "panel-1"), 100)["feasible"]


@pytest.mark.parametrize("field,value", [("beam_width", 0), ("max_expansions", True),
    ("max_alternatives_per_question", -1), ("exact_action_cap", 23),
    ("exact_resource_cap", 129), ("exact_state_cap", 2_000_001)])
def test_search_bounds_are_mandatory_positive_and_cannot_expand_the_hard_tiny_limits(field, value):
    values = {"beam_width": 10, "max_expansions": 100, "max_alternatives_per_question": 3, "exact_action_cap": 16}
    values[field] = value
    with pytest.raises(PopulationError):
        PlannerLimits(**values)


def test_missing_beam_protocol_values_fail_at_the_field():
    with pytest.raises(PopulationError, match="beam_width"):
        PlannerLimits.from_policy({})


def mathematical_question_fixture(requirements, prices, *, weights=None, capacities=None):
    """Build a fixed slot graph for planner tests."""
    qset = accounting_fixture()
    qset.study = {"config": {"policies": {"planning": {"cost_unit": "synthetic-units", "max_expansions": 100000}}}}
    template = qset.actions["a"]
    resource = qset.resources["panel-1"]
    qset.actions, qset.resources, qset.slots, qset.questions = {}, {}, {}, {}
    qset.eligibility, qset.donor_eligibility, qset.options = {}, {}, {}
    qset.evidence_rejections, qset.observation_tokens, qset.records = {}, {}, []
    for aid, (rid, cost) in prices.items():
        qset.actions[aid] = {**deepcopy(template), "action_id": aid, "replicate_position": 1,
                            "resource_demands": [{"resource_id": rid, "dimension": "libraries", "quantity": 1}]}
        qset.eligibility[aid] = {"eligible": True, "reason_codes": []}
        qset.resources[rid] = {**deepcopy(resource), "resource_id": rid, "cost": cost,
            "capacities": [{"dimension": "libraries", "capacity": {"kind": "finite", "quantity": (capacities or {}).get(rid, 20)}}]}
    for qid, donors in requirements.items():
        qset.questions[qid] = {"record_id": qid, "priority_weight": (weights or {}).get(qid, 1),
            "hard_required": True, "requirement_slot_ids": [], "matching_covariates": [],
            "design_ids": [],
            "arms": [{"arm_id": "resolved", "membership": "all_resolved", "required_independent_donors": 1}],
            "evidence_reuse": {"allow_across_questions": True}}
        for donor, aids in donors.items():
            qset.donor_eligibility[qid, donor] = {"eligible": True, "reason_codes": [], "assumptions": [],
                "membership": "all_resolved", "independence_group_id": donor, "covariates": {}}
            for number, aid in enumerate(aids):
                sid = qid + "/" + donor + "/" + str(number)
                qset.actions[aid]["donor_id"] = donor
                qset.slots[sid] = {"slot_id": sid, "question_id": qid, "donor_id": donor,
                    "arm_id": "resolved", "design_id": {"value": aid}, "endpoint_id": "endpoint",
                    "context_id": "context", "locus_id": {"value": "locus"}, "control_role": "measurement"}
                qset.slots[sid]["replicate_position"] = 1
                qset.questions[qid]["requirement_slot_ids"].append(sid)
                qset.options[sid] = {"action_ids": [aid], "observation_ids": []}
    return qset


def test_joint_seed_keeps_the_shared_feasible_alternative_instead_of_two_conflicting_cheapest_bundles():
    qset = mathematical_question_fixture({"A": {"C": ["ca"], "N": ["na"]},
                                         "B": {"C": ["cb"], "N": ["nb"]}},
        {"ca": ("C-panel", 2), "cb": ("C-panel", 2), "na": ("N-panel", 5), "nb": ("N-panel", 5)},
        capacities={"C-panel": 1, "N-panel": 2})
    search = _Search(qset, 5, PlannerLimits(32, 100000, 8, 16))
    a, _ = _construct(search, _Selection(), ["A"])
    b, _ = _construct(search, _Selection(), ["B"])
    assert a[0].actions == ("ca",) and b[0].actions == ("cb",)
    assert not search.account(_selection(("ca", "cb")))["feasible"]
    joint, stats = _construct(search, _Selection(), ["A", "B"])
    assert stats["exhaustive"] and joint[0].actions == ("na", "nb")
    bounded, stats = _bounded(search, _Selection(), ["A", "B"])
    assert bounded[0].actions == ("na", "nb") and not stats["exhaustive"]


@pytest.mark.parametrize("selector", ["complete_bundle_greedy", "donor_first", "adapted_coverage", "canonical_bundle", "random_bundle"])
def test_all_nonexact_selectors_preserve_a_feasible_required_seed_against_high_ratio_optional_work(selector):
    qset = mathematical_question_fixture({"required": {"D": ["h1", "h2"]}, "optional": {"O": ["o"]}},
        {"h1": ("h1", 1), "h2": ("h2", 1), "o": ("o", 1)}, weights={"required": 1, "optional": 1000})
    qset.questions["optional"]["hard_required"] = False
    search = _Search(qset, 2, PlannerLimits(32, 100000, 8, 16))
    seeds, _ = _construct(search, _Selection(), ["required"])
    assert seeds[0].actions == ("h1", "h2")
    selected, _ = _select_optional(search, seeds[0], selector, 19, [])
    assert selected.actions == ("h1", "h2")
    assert search.evaluate(selected)["completed_question_ids"] == ["required"]


def test_complete_bundle_adds_complementary_actions_with_zero_isolated_gain_and_prefers_zero_cost():
    qset = mathematical_question_fixture({"free": {"D": ["a", "b"]}, "paid": {"P": ["p"]}},
        {"a": ("free", 0), "b": ("free", 0), "p": ("paid", 1)})
    search = _Search(qset, 1, PlannerLimits(32, 100000, 8, 16))
    assert search.evaluate(_selection(["a"]))["weighted_completion"] == 0
    assert search.evaluate(_selection(["b"]))["weighted_completion"] == 0
    selected, steps = _select_optional(search, _Selection(), "complete_bundle_greedy", 0, [])
    assert steps[0]["added_action_ids"] == ["a", "b"] and steps[0]["incremental_cost"] == 0
    assert selected.actions == ("a", "b", "p")


def test_bounded_miss_and_interrupted_exact_enumeration_never_become_exhaustive():
    qset = mathematical_question_fixture({"Q": {"D": ["a", "b"]}}, {"a": ("a", 1), "b": ("b", 1)})
    bounded = _Search(qset, 2, PlannerLimits(1, 1, 1, 1))
    candidates, stats = _construct(bounded, _Selection(), ["Q"])
    assert not candidates and not stats["exhaustive"] and stats["expansion_limit_reached"]
    exact = _Search(qset, 2, PlannerLimits(2, 1, 1, 3))
    candidates, stats = _exact(exact, _Selection(), ["Q"])
    assert not candidates and not stats["exhaustive"] and stats["expansion_limit_reached"]


@pytest.fixture
def source_case(tmp_path):
    return full_fixture(tmp_path / "literal-section9-source")


@pytest.mark.parametrize("selector", SELECTORS)
def test_source_bound_section9_all_selectors_complete_at18_and_preserve_no_measurements(source_case, selector):
    qset = source_case["question_set"]
    original = deepcopy(qset.records)
    plan = plan_experiments(qset, budget=18, selector=selector)
    assert plan["meets_request"] and plan["weighted_completion"] == plan["unweighted_completion"] == 1
    assert plan["accounting"]["total_cost"] == 18
    assert len(plan["selected_action_ids"]) == 10
    assert all(row["planning_status"] == "fully_planned" and row["evidence_status"] == "no_measurements" for row in plan["questions"].values())
    assert qset.records == original
    assert evaluate_plan(qset, plan["selected_action_ids"], plan["resource_assignments"], budget=18)["meets_request"]
    assert plan["mandatory_seed"]["total_cost"] == (None if selector == "exact_tiny" else 18)
    assert validate_record(plan["plan_record"]) == plan["plan_record"]
    validate_catalog([*qset.records, plan["plan_record"]])


def test_source_bound_section9_independent_exhaustive_oracle_proves18_minimum_and17_infeasible(source_case):
    qset = source_case["question_set"]
    oracle = TinyOracle(qset.study["config"], qset.records)
    expected = oracle.solve(18)
    assert expected["exhaustive"] and expected["states"] == 1024 and expected["minimum_full_cost"] == 18
    plan = plan_experiments(qset, budget=18, selector="exact_tiny")
    assert plan["optimality_established"]
    for key in ("selected_action_ids", "resource_assignments", "weighted_completion", "unweighted_completion"):
        assert plan[key] == expected["best"][key]
    assert plan["accounting"]["total_cost"] == expected["best"]["total_cost"]
    impossible = oracle.solve(17)
    assert impossible["exhaustive"] and impossible["status"] == "infeasible_proven" and impossible["best"] is None
    plan17 = plan_experiments(qset, budget=17)
    assert plan17["status"] == "infeasible_proven" and not plan17["meets_request"]
    assert plan17["infeasibility_certificate"]["kind"] == "exhaustive_declared_action_and_instance_catalog"
    for label, aid in source_case["labels"]["actions"].items():
        if label.endswith(":measurement_panel"):
            subset = oracle.dependencies(aid, set())
            assert oracle.evaluate(subset, {})["weighted_completion"] == 0


def test_source_bound_forced_beam_finds18_but_only_reports_bounded_miss_at17(source_case):
    qset = source_case["question_set"]
    limits = PlannerLimits(64, 10000, 8, 1)
    plan = plan_experiments(qset, budget=18, limits=limits)
    assert plan["meets_request"] and plan["accounting"]["total_cost"] == 18
    assert not plan["optimality_established"]
    missing = plan_experiments(qset, budget=17, limits=limits)
    assert missing["status"] == "no_feasible_bundle_found" and missing["infeasibility_certificate"] is None
    limited = plan_experiments(qset, budget=18, limits=PlannerLimits(1, 1, 1, 1))
    assert limited["status"] == "no_feasible_bundle_found" and not limited["meets_request"]
    assert any(row["expansion_limit_reached"] for row in limited["searches"])


def test_source_bound_joint_seed_uses_shared_alternative_and_independent_oracle_agrees(source_case, tmp_path):
    qset = joint_seed_case(source_case, tmp_path / "joint-seed")
    oracle = TinyOracle(qset.study["config"], qset.records).solve(5)
    plan = plan_experiments(qset, budget=5)
    assert plan["meets_request"] and plan["weighted_completion"] == 2 and plan["mandatory_seed"]["total_cost"] == 5
    assert oracle["exhaustive"] and oracle["best"]["selected_action_ids"] == plan["selected_action_ids"]
    assert {qset.actions[aid]["donor_id"] for aid in plan["selected_action_ids"]} == {source_case["labels"]["donors"]["N"]}
    search = _Search(qset, 5, PlannerLimits(64, 100000, 16, 16))
    separately = [_construct(search, _Selection(), [qid])[0][0] for qid in sorted(qset.questions)]
    assert [search.account(selected)["total_cost"] for selected in separately] == [2, 2]
    naive = _selection(set().union(*(set(selected.actions) for selected in separately)))
    assert not search.account(naive)["feasible"]
    bounded = plan_experiments(qset, budget=5, limits=PlannerLimits(64, 10000, 16, 1))
    assert bounded["meets_request"] and bounded["accounting"]["total_cost"] == 5


@pytest.mark.parametrize("selector", SELECTORS[:-1])
def test_source_bound_optional_high_weight_cannot_spend_a_required_seed(source_case, tmp_path, selector):
    qset = optional_trap_case(source_case, tmp_path / "optional-trap")
    separate = [aid for aid, action in qset.actions.items() if action["experiment_id"] == "optional-new-experiment"]
    alone = evaluate_plan(qset, separate, budget=18)
    assert alone["feasible"] and alone["total_cost"] == 6 and alone["requirements"]["weighted_completion"] == 1000
    plan = plan_experiments(qset, budget=18, selector=selector)
    assert plan["meets_request"] and plan["mandatory_seed"]["total_cost"] == 18
    assert plan["weighted_completion"] == 1 and plan["unweighted_completion"] == 1
    assert len(plan["selected_action_ids"]) == 10
    assert all(qset.actions[aid]["experiment_id"] != "optional-new-experiment" for aid in plan["selected_action_ids"])


def test_source_bound_alternative_instances_pay_and_consume_the_same_instance(source_case, tmp_path):
    qset = alternatives_case(source_case, tmp_path / "alternative-instances")
    plan = plan_experiments(qset, budget=26, selector="exact_tiny")
    assert plan["meets_request"] and plan["accounting"]["total_cost"] == 26
    assert len(plan["resource_assignments"]) == 6
    amounts = {}
    for row in plan["resource_assignments"]:
        amounts[row["resource_id"]] = amounts.get(row["resource_id"], 0) + 1
    assert sorted(amounts.values()) == [3, 3]
    expected = TinyOracle(qset.study["config"], qset.records).solve(26)
    assert expected["exhaustive"] and expected["minimum_full_cost"] == 26
    assert expected["best"]["resource_assignments"] == plan["resource_assignments"]
    bad = deepcopy(plan["resource_assignments"])
    for row in bad:
        row["resource_id"] = next(iter(amounts))
    assert not evaluate_plan(qset, plan["selected_action_ids"], bad, budget=26)["feasible"]
    assert plan_experiments(qset, budget=25)["status"] == "infeasible_proven"
    bounded = plan_experiments(qset, budget=26, limits=PlannerLimits(64, 10000, 16, 1))
    assert bounded["meets_request"] and bounded["accounting"]["total_cost"] == 26


def test_optional_trap_independent_oracle_retains_the_required_lower_weight_plan(source_case, tmp_path):
    qset = optional_trap_case(source_case, tmp_path / "oracle-optional-trap")
    expected = TinyOracle(qset.study["config"], qset.records).solve(18)
    assert expected["exhaustive"] and expected["best"]["weighted_completion"] == 1
    plan = plan_experiments(qset, budget=18, selector="exact_tiny")
    assert plan["optimality_established"] and plan["selected_action_ids"] == expected["best"]["selected_action_ids"]


def test_source_bound_question_construction_limit_cannot_prove_infeasibility_or_optimality(source_case, tmp_path):
    rows = deepcopy(source_case["question_set"].records)
    question = next(row for row in rows if row["record_type"] == "comparison_question")
    question.update(construction_status="excluded_by_declared_limit", requirement_slot_ids=[], planning_status="blocked",
                    remaining_gaps=["requirement_construction_shortlist_exclusion"])
    rows = [row for row in rows if row["record_type"] != "comparison_requirement_slot"]
    qset = write_case(source_case, tmp_path / "construction-excluded", rows)
    for selector in ("complete_bundle_greedy", "exact_tiny"):
        plan = plan_experiments(qset, selector=selector)
        assert plan["status"] == "no_feasible_bundle_found" and not plan["optimality_established"]
        assert plan["infeasibility_certificate"] is None and plan["construction_limited_question_ids"]
        assert all(row["planning_status"] == "no_feasible_bundle_found" for row in plan["questions"].values())
    question["hard_required"] = False
    qset_optional = write_case(source_case, tmp_path / "optional-construction-excluded", rows)
    exact = plan_experiments(qset_optional, selector="exact_tiny")
    assert exact["status"] == "feasible_search_incomplete" and not exact["optimality_established"]


def test_seeded_random_plan_replays_and_relocating_identical_source_bytes_preserves_identity(source_case, tmp_path):
    optional = write_case(source_case, tmp_path / "two-optional", split_design_questions(source_case, required=False), budget=10)
    first = plan_experiments(optional, selector="random_bundle", random_seed=73)
    second = plan_experiments(optional, selector="random_bundle", random_seed=73)
    assert first["plan_id"] == second["plan_id"] and first["selection_steps"] == second["selection_steps"]
    assert first["selection_steps"] and all(step["candidate_order"] for step in first["selection_steps"])
    relocated = full_fixture(tmp_path / "relocated-identical-sources")
    original = plan_experiments(source_case["question_set"], budget=18)
    replay = plan_experiments(relocated["question_set"], budget=18)
    assert original["source_sha256"] == replay["source_sha256"] and original["plan_id"] == replay["plan_id"]


def test_plan_revisions_are_immutable_and_resealed_solver_claims_fail_recomputation(source_case, tmp_path):
    qset = source_case["question_set"]
    plan = plan_experiments(qset, budget=18)
    path = write_plan_revision(plan, tmp_path / "revisions")
    assert write_plan_revision(plan, path.parent) == path
    assert verify_plan_revision(qset, plan)["plan_id"] == plan["plan_id"]
    altered = deepcopy(plan)
    altered["accounting"]["total_cost"] = 1
    altered["plan_record"] = _plan_record(altered)
    altered["plan_id"] = altered["plan_record"]["plan_id"]
    with pytest.raises(PopulationError, match="recomputation"):
        verify_plan_revision(qset, altered)
    altered["accounting"]["total_cost"] = 2
    with pytest.raises(PopulationError, match="identity"):
        write_plan_revision(altered, path.parent)
    revised = plan_experiments(qset, budget=17, previous_plan_id=plan["plan_id"])
    second = write_plan_revision(revised, path.parent)
    assert second != path and json.loads(path.read_text())["accounting"]["total_cost"] == 18
    assert json.loads(second.read_text())["previous_plan_id"] == plan["plan_id"]
