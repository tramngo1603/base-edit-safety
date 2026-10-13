"""Budget conservation, distinct donors and the deliberate greedy counterexample."""
from copy import deepcopy
from fractions import Fraction
from pathlib import Path
import random

import pytest

from base_edit_safety.population.common import PopulationError, parse_json, read_json
from base_edit_safety.population.contracts import validate_request
from base_edit_safety.population.planning import plan
from base_edit_safety.population.resolver import resolve_request

FIXTURE = Path(__file__).resolve().parents[1]/"examples/population/planner"


def request():
    return read_json(FIXTURE/"request.json")


@pytest.mark.parametrize("algorithm", ["greedy_cost", "exact_tiny"])
def test_activation_trap(algorithm):
    result, _ = plan(request(), algorithm)
    assert result == read_json(FIXTURE/f"expected-{algorithm}.json")


@pytest.mark.parametrize("algorithm", ["greedy_cost", "greedy_coverage", "canonical_order", "random_order", "exact_tiny"])
def test_permutation_invariance(algorithm):
    original = request()
    expected = plan(original, algorithm, 7)
    for seed in range(10):
        r = deepcopy(original)
        for field in ("units", "donors", "actions"):
            random.Random(seed).shuffle(r[field])
        assert plan(r, algorithm, 7) == expected


@pytest.mark.parametrize("field,value", [("max_cost", -1), ("max_cost", True), ("max_donors", -1), ("max_assays", "1")])
def test_budget_types(field, value):
    r = request(); r["budget"][field] = value
    with pytest.raises(PopulationError):
        plan(r)


def test_unknown_fields_possible_links_and_model_physical_planning_rejected():
    for mutate in [lambda r: r.update(hard_coverage=True),
                   lambda r: r["actions"][0]["possible_only_units"].append("u1"),
                   lambda r: r.update(mode="research_sample_planning", coverage_basis="retained_callset_model"),
                   lambda r: r["actions"].append(r["actions"][0]),
                   lambda r: r["donors"][0].update(activation_cost=True)]:
        r = request(); mutate(r)
        with pytest.raises(PopulationError):
            plan(r)


def test_prior_evidence_empty_plan_and_distinct_donor_fraction():
    r = request()
    r["existing_usable_evidence"] = [{"unit_id": u["id"], "donor_id": "historical", "evidence_id": "synthetic_prior"} for u in r["units"]]
    result, _ = plan(r)
    assert result["plan_status"] == "no_remaining_demand"
    assert result["total_cost"] == 0 and not result["selected_donors"]
    r["units"][0]["required_donors"] = 3
    r["existing_usable_evidence"] *= 2
    r["budget"]["max_cost"] = 0
    result, _ = plan(r)
    assert Fraction(result["objective_numerator"], result["objective_denominator"]) == Fraction(7, 3)


def test_empty_status_precedence_and_free_action():
    r = request(); r["budget"]["max_cost"] = 0
    assert plan(r)[0]["plan_status"] == "no_feasible_action"
    r["donors"][2]["activation_cost"] = 0
    r["actions"][3]["assay_cost"] = 0
    assert plan(r)[0]["plan_status"] == "no_positive_gain"
    r["actions"][4]["assay_cost"] = 0
    result, _ = plan(r)
    assert result["selected_actions"] == ["C3"] and result["total_cost"] == 0


def test_max_assays_and_donors_charge_activation():
    r = request(); r["budget"].update(max_cost=6, max_donors=1, max_assays=2)
    result, _ = plan(r, "exact_tiny")
    assert result["selected_actions"] == ["A1", "A2"] and result["total_cost"] == 4
    r["budget"]["max_assays"] = 0
    assert plan(r)[0]["plan_status"] == "no_feasible_action"


def test_duplicate_assays_and_chromosomes_do_not_add_donors():
    r = request(); r["units"] = [r["units"][0]]; r["units"][0]["required_donors"] = 2
    r["actions"] = [r["actions"][0], {**r["actions"][0], "id": "replicate"}]
    result, _ = plan(r)
    assert len(result["selected_actions"]) == 1
    assert result["unit_coverage"][0]["combined_donor_count"] == 1
    assert result["objective_denominator"] == 2


def test_physical_availability_and_permission_are_enforced():
    r = request(); r.update(mode="research_sample_planning", coverage_basis="callable_supported")
    for donor in r["donors"]:
        donor["availability"] = "unverified"
    assert plan(r)[0]["plan_status"] == "no_feasible_action"
    r["donors"][0]["availability"] = "verified_available"
    r["donors"][0]["use_permitted"] = False
    assert plan(r)[0]["plan_status"] == "no_feasible_action"


def test_exact_tiny_is_bounded():
    r = request(); r["actions"] = [{**r["actions"][0], "id": str(i)} for i in range(13)]
    with pytest.raises(PopulationError, match="12"):
        plan(r, "exact_tiny")


def test_schema_and_mandatory_manifest():
    with pytest.raises(PopulationError):
        parse_json('{"cost": 1, "cost": 2}')
    with pytest.raises(PopulationError):
        parse_json('{"cost": NaN}')
    with pytest.raises(PopulationError):
        parse_json('{"cost": 1e999}')
    assert resolve_request(FIXTURE/"request.json", FIXTURE/"planning-manifest.json")[0] == request()


def test_random_tiny_problems_never_beat_independent_oracle():
    for seed in range(40):
        rng = random.Random(seed)
        r = request()
        for donor in r["donors"]:
            donor["activation_cost"] = rng.randrange(4)
        for action in r["actions"]:
            action["assay_cost"] = rng.randrange(4)
            action["possible_only_units"] = []
            action["covered_units"] = [u["id"] for u in r["units"] if rng.choice([True, False])]
        for unit in r["units"]:
            unit["required_donors"] = rng.randrange(1, 4)
            unit["weight"] = rng.randrange(1, 5)
        oracle, _ = plan(r, "exact_tiny")
        optimum = Fraction(oracle["objective_numerator"], oracle["objective_denominator"])
        for algorithm in ("greedy_cost", "greedy_coverage", "random_order", "canonical_order"):
            result, _ = plan(r, algorithm, seed)
            assert Fraction(result["objective_numerator"], result["objective_denominator"]) <= optimum
            assert result["total_cost"] + result["unused_budget"] == r["budget"]["max_cost"]
