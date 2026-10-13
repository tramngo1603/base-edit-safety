"""Hand-calculated mock evaluations and checks for split leakage."""
from copy import deepcopy
from pathlib import Path

import pytest

from base_edit_safety.population.common import PopulationError, read_json
from base_edit_safety.population.comparison_evaluation import (
    dependency_groups, evaluate, interval, main, make_splits, metrics, validate_splits,
)

FIXTURE = Path(__file__).parent / "fixtures/comparison/evaluation"


def inputs():
    rows, labels, methods, protocol = [read_json(FIXTURE / (name + ".json"))
                                     for name in ("candidates", "labels", "methods", "protocol")]
    return rows, labels, methods, protocol, make_splits(rows, 17, "0", "0")


def report(values=None):
    return evaluate(*(values or inputs()))


def test_hand_calculated_enumeration_ranking_and_matched_sensitivity():
    result = report()["strata"][0]
    proposed, baseline = [result["methods"][name] for name in ("proposed", "baseline")]
    assert proposed["eligible_positives"] == baseline["eligible_positives"] == 4
    assert proposed["enumeration_sensitivity"] == 3 / 4
    assert proposed["missed_enumeration_positives"] == 1
    assert proposed["enumerated_unscored_positives"] == 1
    assert proposed["at_budget"] == {"assays_used": 4, "positive": 2, "tested_no_signal": 0,
        "inconclusive": 1, "untested": 1, "conclusive_count": 2, "recall": 0.5, "precision_conclusive": 1.0, "conclusive_outcome_coverage": 0.5}
    assert baseline["at_budget"]["recall"] == 0.5
    assert baseline["at_budget"]["precision_conclusive"] == 0.5
    assert baseline["enumeration_sensitivity"] == 1
    assert proposed["matched_sensitivity"]["0.5"]["assays"] == 4
    assert baseline["matched_sensitivity"]["0.5"]["assays"] == 4
    assert proposed["matched_sensitivity"]["0.75"]["status"] == "unreachable"
    assert baseline["matched_sensitivity"]["0.75"]["assays"] == 6
    assert baseline["matched_sensitivity"]["1"]["assays"] == 8
    assert proposed["precision_recall"][-1]["recall"] == 0.5
    assert len(proposed["precision_recall"]) == 7
    comparison = result["matched_sensitivity_comparisons"]["proposed"]
    assert comparison["0.5"]["assay_difference_method_minus_baseline"] == 0
    assert comparison["0.75"]["status"] == "not_comparable"
    assert comparison["0.75"]["assay_difference_method_minus_baseline"] is None
    assert proposed["uncertainty"]["recall_at_budget"]["independent_clusters"] == 3
    assert proposed["uncertainty"]["recall_at_budget"]["defined_replicates"] == 200
    assert proposed["uncertainty"]["precision_conclusive_at_budget"]["reason"] == "fewer_than_three_relevant_dependency_groups"


def test_paired_recall_interval_uses_the_same_hand_calculated_draws(monkeypatch):
    values = inputs()
    assignments = validate_splits(values[0], values[4])
    groups = [assignments[key][0] for key in ("c1", "c4", "c6")]
    draws = iter(groups * 67)
    monkeypatch.setattr("random.Random.choices", lambda self, population, k: [next(draws)] * k)
    result = report(values)["strata"][0]
    paired = result["recall_at_budget_comparisons"]["proposed"]
    # Repeating group A: (1 - 1)/3 = 0; B: (1 - 3)/3 = -2/3;
    # C: (0 - 1)/6 = -1/6. Repeated unresolved rows consume the budget too.
    assert paired == {
        "baseline_method": "baseline", "status": "comparable", "budget": 4,
        "eligible_positives": 4, "method_recall": 0.5, "baseline_recall": 0.5,
        "recall_difference_method_minus_baseline": 0,
        "uncertainty": {"independent_clusters": 3, "replicates": 200,
            "defined_replicates": 200, "undefined_replicates": 0,
            "lower": -2 / 3, "upper": 0, "reason": None},
    }
    current, baseline = [result["methods"][name]["uncertainty"]["recall_at_budget"]
                         for name in ("proposed", "baseline")]
    assert paired["uncertainty"]["lower"] != current["lower"] - baseline["upper"]


def test_hand_calculated_positive_recall_interval_is_a_software_result(monkeypatch):
    values = inputs()
    values[3]["budget_per_editor_assay"] = 2
    for method, order in zip(values[2], [
        ["c1", "c4", "c6", "c8", "c2", "c3", "c5", "c7"],
        ["c2", "c3", "c5", "c7", "c1", "c4", "c6", "c8"],
    ]):
        method["enumerated_ids"] = order
        method["scores"] = [{"candidate_id": key, "score": -i} for i, key in enumerate(order)]
    assignments = validate_splits(values[0], values[4])
    draws = iter([assignments["c1"][0], assignments["c6"][0]] * 100)
    monkeypatch.setattr("random.Random.choices", lambda self, population, k: [next(draws)] * k)
    result = report(values)
    paired = result["strata"][0]["recall_at_budget_comparisons"]["proposed"]
    # The method selects two positives; baseline selects two unresolved/no-signal
    # rows. Repeated A has three positives, repeated C six: differences 2/3, 1/3.
    assert paired["recall_difference_method_minus_baseline"] == 0.5
    assert (paired["uncertainty"]["lower"], paired["uncertainty"]["upper"]) == (1 / 3, 2 / 3)
    assert result["data_kind"] == "artificial"
    assert "no measured accuracy or assay savings claim" in result["scope"]


@pytest.mark.parametrize("baseline,method,expected", [
    ("baseline", "proposed", -1 / 4), ("proposed", "baseline", 1 / 4),
])
def test_paired_recall_sign_and_denominator_include_missed_positives(baseline, method, expected):
    values = inputs()
    values[3].update(baseline_method=baseline, budget_per_editor_assay=6)
    result = report(values)["strata"][0]
    paired = result["recall_at_budget_comparisons"][method]
    assert set(result["recall_at_budget_comparisons"]) == {method}
    assert paired["recall_difference_method_minus_baseline"] == expected
    assert paired["eligible_positives"] == 4
    assert result["methods"]["proposed"]["missed_enumeration_positives"] == 1
    assert result["methods"]["proposed"]["enumerated_unscored_positives"] == 1


@pytest.mark.parametrize("groups", [1, 2])
def test_paired_recall_withholds_interval_with_insufficient_positive_groups(groups):
    values = list(inputs())
    for row in values[0]:
        if groups == 1 or row["study_id"] != "artificial-study-3":
            row["family_id"] = "artificial-shared-family"
    values[4] = make_splits(values[0], 17, "0", "0")
    paired = report(values)["strata"][0]["recall_at_budget_comparisons"]["proposed"]
    assert paired["status"] == "comparable"
    assert paired["recall_difference_method_minus_baseline"] == 0
    assert paired["uncertainty"]["independent_clusters"] == groups
    assert paired["uncertainty"]["lower"] is paired["uncertainty"]["upper"] is None
    assert paired["uncertainty"]["reason"] == "fewer_than_three_relevant_dependency_groups"


def test_paired_recall_has_no_difference_when_original_denominator_is_undefined():
    values = inputs()
    for label in values[1]:
        label["adjudicated_outcome"] = "inconclusive"
    paired = report(values)["strata"][0]["recall_at_budget_comparisons"]["proposed"]
    assert paired["status"] == "not_comparable"
    assert paired["eligible_positives"] == 0
    assert paired["method_recall"] is paired["baseline_recall"] is None
    assert paired["recall_difference_method_minus_baseline"] is None
    assert paired["uncertainty"]["defined_replicates"] == 0
    assert paired["uncertainty"]["undefined_replicates"] == 200
    assert paired["uncertainty"]["reason"] == "original_metric_undefined_or_unreachable"
    assert paired["uncertainty"]["lower"] is paired["uncertainty"]["upper"] is None


def test_paired_recall_keeps_undefined_resamples_with_three_positive_groups(monkeypatch):
    values = list(inputs())
    values[0].append({**values[0][0], "candidate_id": "c9", "study_id": "artificial-study-4",
                      "guide_id": "artificial-guide-4", "family_id": "artificial-family-4"})
    values[1].append({**values[1][0], "candidate_id": "c9", "adjudicated_outcome": "untested"})
    values[4] = make_splits(values[0], 17, "0", "0")
    assignments = validate_splits(values[0], values[4])
    empty, positive = assignments["c9"][0], assignments["c1"][0]
    draws = iter([empty, positive] * 100)
    monkeypatch.setattr("random.Random.choices", lambda self, population, k: [next(draws)] * k)
    paired = report(values)["strata"][0]["recall_at_budget_comparisons"]["proposed"]
    assert paired["status"] == "comparable"
    assert paired["eligible_positives"] == 4
    assert paired["uncertainty"]["independent_clusters"] == 3
    assert paired["uncertainty"]["defined_replicates"] == 100
    assert paired["uncertainty"]["undefined_replicates"] == 100
    assert paired["uncertainty"]["reason"] == "undefined_replicates_present"
    assert paired["uncertainty"]["lower"] is paired["uncertainty"]["upper"] is None


@pytest.mark.parametrize("budget,expected", [(0, 0), (4, -0.5)])
def test_paired_recall_retains_denominator_for_empty_scores_and_zero_budget(budget, expected):
    values = inputs()
    values[3]["budget_per_editor_assay"] = budget
    values[2][0]["scores"] = []
    paired = report(values)["strata"][0]["recall_at_budget_comparisons"]["proposed"]
    assert paired["status"] == "comparable"
    assert paired["eligible_positives"] == 4
    assert paired["recall_difference_method_minus_baseline"] == expected


def test_paired_recall_preserves_exact_score_order_and_is_deterministic():
    values = inputs()
    values[3]["budget_per_editor_assay"] = 1
    values[2][0]["scores"] = [
        {"candidate_id": "c2", "score": "1.000000000000000000000000000000"},
        {"candidate_id": "c4", "score": "1.000000000000000000000000000001"},
    ]
    first = report(values)
    assert first == report(values)
    paired = first["strata"][0]["recall_at_budget_comparisons"]["proposed"]
    assert paired["method_recall"] == 1 / 4
    assert paired["baseline_recall"] == 0
    assert paired["recall_difference_method_minus_baseline"] == 1 / 4
    values[2].reverse()
    assert report(values)["strata"] == first["strata"]


def test_splits_are_label_free_deterministic_and_order_independent():
    rows, *_ = inputs()
    a = make_splits(rows, 45, "0.34", "0.34")
    assert a == make_splits(list(reversed(rows)), 45, "0.34", "0.34")
    assert {g["partition"] for g in a["groups"]} == {"training", "development", "evaluation"}
    assert len(validate_splits(rows, a)) == 8
    a["groups"][0]["partition"] = "evaluation"
    with pytest.raises(PopulationError, match="manifest"):
        validate_splits(rows, a)


def test_transitive_ineligible_bridge_links_distinct_studies_and_guides_before_filter():
    values = list(inputs())
    rows = values[0]
    bridge = {**rows[0], "candidate_id": "bridge", "study_id": "bridge-study",
              "guide_id": rows[0]["guide_id"], "family_id": rows[3]["family_id"], "eligible": False}
    rows.append(bridge)
    groups = list(dependency_groups(rows).values())
    assert len(groups) == 2
    assert any({"c1", "c4", "bridge"} <= set(group) for group in groups)
    # The old split is invalid even though the new row is ineligible.
    with pytest.raises(PopulationError, match="manifest"):
        report(values)
    values[1].append({**values[1][0], "candidate_id": "bridge"})
    values[4] = make_splits(rows, 17, "0", "0")
    result = report(values)
    assert result["evaluation_dependency_groups"] == 2
    assert result["evaluation_candidates"] == 8
    assert result["strata"][0]["methods"]["proposed"]["eligible_positives"] == 4
    assert result["strata"][0]["methods"]["proposed"]["uncertainty"]["recall_at_budget"]["lower"] is None


def test_training_and_ineligible_positives_do_not_enter_evaluation_denominator():
    values = list(inputs())
    values[0][0]["eligible"] = False
    values[4] = make_splits(values[0], 17, "0.34", "0")
    assignments = validate_splits(values[0], values[4])
    expected = sum(row["adjudicated_outcome"] == "observed_positive" and
        values[0][i]["eligible"] and assignments[row["candidate_id"]][1] == "evaluation"
        for i, row in enumerate(values[1]))
    result = report(values)["strata"][0]
    assert all(method["eligible_positives"] == expected for method in result["methods"].values())
    assert result["recall_at_budget_comparisons"]["proposed"]["eligible_positives"] == expected
    assert expected < 4


def test_ties_charge_each_assay_and_use_candidate_id_not_input_order():
    values = list(inputs())
    values[3]["budget_per_editor_assay"] = 1
    values[2][0]["scores"] = [{"candidate_id": "c3", "score": "1.00"}, {"candidate_id": "c1", "score": "1"}]
    result = report(values)["strata"][0]["methods"]["proposed"]
    assert result["at_budget"]["positive"] == result["at_budget"]["assays_used"] == 1
    values[2][0]["scores"].reverse()
    assert report(values)["strata"][0]["methods"]["proposed"] == result
    values[3]["budget_per_editor_assay"] = 4
    result = report(values)["strata"][0]["methods"]["proposed"]
    assert result["budget_unfilled"] == 2
    assert result["at_budget"]["assays_used"] == 2
    assert result["at_budget"]["conclusive_outcome_coverage"] == 0.5


def test_cluster_resampling_preserves_multiplicity_in_cost_and_denominator():
    rows, labels, methods, _, _ = inputs()
    label_index = {r["candidate_id"]: r for r in labels}
    result = metrics({"c1", "c3"}, label_index, methods[0], 3, ["1"], {"c1": 2, "c3": 2})
    # c3 (inconclusive) twice, then two copies of positive c1.
    assert result["eligible_positives"] == 2
    assert result["ranked_candidates"] == 4
    assert result["at_budget"]["inconclusive"] == 2
    assert result["at_budget"]["positive"] == 1
    assert result["at_budget"]["recall"] == 0.5
    assert result["at_budget"]["conclusive_outcome_coverage"] == 1 / 3
    assert result["matched_sensitivity"]["1"]["assays"] == 4


def test_paired_draws_are_shared_across_methods_and_strata():
    values = list(inputs())
    values[2][0] = {**deepcopy(values[2][1]), "method_id": "proposed"}
    # Add a second assay with the same dependency groups and outcomes.
    for row, label in zip(deepcopy(values[0]), deepcopy(values[1])):
        old = row["candidate_id"]
        row.update(candidate_id="z" + old, assay="artificial-assay-B")
        label["candidate_id"] = row["candidate_id"]
        values[0].append(row)
        values[1].append(label)
        for method in values[2]:
            method["enumerated_ids"].append(row["candidate_id"])
            score = next(score for score in method["scores"] if score["candidate_id"] == old)
            method["scores"].append({**score, "candidate_id": row["candidate_id"]})
    values[4] = make_splits(values[0], 17, "0", "0")
    result = report(values)
    assert result["evaluation_dependency_groups"] == 3
    first, second = result["strata"]
    assert first["methods"]["proposed"] == second["methods"]["proposed"]
    for stratum in result["strata"]:
        paired_recall = stratum["recall_at_budget_comparisons"]["proposed"]
        assert paired_recall["recall_difference_method_minus_baseline"] == 0
        assert paired_recall["uncertainty"]["lower"] == paired_recall["uncertainty"]["upper"] == 0
        for comparison in stratum["matched_sensitivity_comparisons"]["proposed"].values():
            assert comparison["uncertainty"]["lower"] == comparison["uncertainty"]["upper"] == 0


def test_undefined_denominators_and_unreachable_resamples_are_not_dropped():
    values = list(inputs())
    for label in values[1]:
        if label["candidate_id"] not in {"c1"}:
            label["adjudicated_outcome"] = "untested"
    result = report(values)["strata"][0]["methods"]["proposed"]
    uncertainty = result["uncertainty"]["recall_at_budget"]
    assert uncertainty["undefined_replicates"] > 0
    assert uncertainty["defined_replicates"] + uncertainty["undefined_replicates"] == 200
    assert uncertainty["lower"] is None
    assert interval([1, None, 0], 3)["reason"] == "undefined_replicates_present"
    for label in values[1]:
        label["adjudicated_outcome"] = "untested"
    result = report(values)["strata"][0]["methods"]["proposed"]
    assert result["enumeration_sensitivity"] is None
    assert result["at_budget"]["conclusive_outcome_coverage"] == 0
    assert result["at_budget"]["recall"] is None
    assert result["at_budget"]["precision_conclusive"] is None
    assert result["matched_sensitivity"]["0"]["assays"] == 0
    assert result["matched_sensitivity"]["0.5"]["status"] == "undefined_no_eligible_positives"
    assert result["uncertainty"]["recall_at_budget"]["undefined_replicates"] == 200


@pytest.mark.parametrize("mutation,message", [
    (lambda v: v[0][0].update(family_id="unknown"), "Missing identity"),
    (lambda v: v[1].pop(), "Every candidate"),
    (lambda v: v[1][0].update(adjudicated_outcome="negative"), "Invalid adjudicated"),
    (lambda v: v[1].append(deepcopy(v[1][0])), "Duplicate"),
    (lambda v: v[2][0]["scores"].append(deepcopy(v[2][0]["scores"][0])), "Duplicate"),
    (lambda v: v[2][0]["scores"][0].update(score="NaN"), "finite"),
    (lambda v: v[2][0]["scores"][0].update(score=True), "finite"),
    (lambda v: v[2][0]["scores"][0].update(candidate_id="c8"), "enumerated"),
    (lambda v: v[2][0]["enumerated_ids"].append("unlisted"), "unknown"),
    (lambda v: v[3].update(data_kind="real"), "artificial"),
    (lambda v: v[3].update(budget_per_editor_assay=-1), "bounds"),
    (lambda v: v[3].update(baseline_method="missing"), "Baseline"),
    (lambda v: v[3].update(targets=["1.01"]), "between"),
])
def test_invalid_or_incomplete_inputs_fail_closed(mutation, message):
    values = list(inputs())
    mutation(values)
    with pytest.raises(PopulationError, match=message):
        report(values)


def test_zero_budget_and_empty_score_lists_keep_positive_denominators():
    values = list(inputs())
    values[3]["budget_per_editor_assay"] = 0
    values[2][0]["scores"] = []
    result = report(values)["strata"][0]["methods"]["proposed"]
    assert result["eligible_positives"] == 4
    assert result["at_budget"]["recall"] == 0
    assert result["at_budget"]["precision_conclusive"] is None
    assert result["matched_sensitivity"]["1"]["status"] == "unreachable"


def test_saved_fixture_split_and_cli_receipts(tmp_path):
    rows, *_ = inputs()
    saved = read_json(FIXTURE / "splits.json")
    assert saved == make_splits(rows, 17, "0", "0")
    split_output = tmp_path / "split"
    assert main(["split", "--candidates", str(FIXTURE / "candidates.json"), "--training-fraction", "0",
                 "--development-fraction", "0", "--output", str(split_output)]) == 0
    assert read_json(split_output / "splits.json") == saved
    command = ["evaluate"]
    for role in ("candidates", "labels", "methods", "protocol", "splits"):
        command += ["--" + role, str(FIXTURE / (role + ".json"))]
    output = tmp_path / "evaluation"
    assert main(command + ["--output", str(output)]) == 0
    assert read_json(output / "evaluation.json") == report()
    receipt = read_json(output / "receipt.json")
    assert len(receipt["inputs"]) == 5
    assert receipt["runtime"]["source_sha256"]
    assert receipt["elapsed_seconds"] >= 0
    with pytest.raises(SystemExit) as exc:
        main(command + ["--output", str(output)])
    assert exc.value.code == 2


def test_decimal_thresholds_and_scores_do_not_round_into_false_ties():
    rows, labels, methods, _, _ = inputs()
    label_index = {row["candidate_id"]: row for row in labels}
    method = {**methods[0], "scores": [
        {"candidate_id": "c1", "score": "1.000000000000000000000000000000"},
        {"candidate_id": "c4", "score": "1.000000000000000000000000000001"},
        {"candidate_id": "c6", "score": "0"}]}
    label_index["c1"]["adjudicated_outcome"] = "tested_no_signal"
    result = metrics({"c1", "c4", "c6"}, label_index, method, 1, ["1"])
    assert result["at_budget"]["positive"] == 1
    label_index["c1"]["adjudicated_outcome"] = "observed_positive"
    target = "0.333333333333333333333333333334"
    result = metrics({"c1", "c4", "c6"}, label_index, method, 1, [target])
    assert result["matched_sensitivity"][target]["required_positives"] == 2
    assert result["matched_sensitivity"][target]["assays"] == 2


def test_split_fraction_uses_exact_floor():
    rows, *_ = inputs()
    split = make_splits(rows, 17, "0.999999999999999999999999999999", "0")
    assert sum(group["partition"] == "training" for group in split["groups"]) == 2
    assert sum(group["partition"] == "evaluation" for group in split["groups"]) == 1


def test_unreachable_original_cannot_gain_interval_from_reachable_resamples(monkeypatch):
    values = inputs()
    assignments = validate_splits(values[0], values[4])
    group = assignments["c1"][0]
    # A possible cluster draw contains only the first study; all its positives are scored.
    monkeypatch.setattr("random.Random.choices", lambda self, population, k: [group] * k)
    stratum = report(values)["strata"][0]
    method = stratum["methods"]["proposed"]
    uncertainty = method["uncertainty"]["assays_at_sensitivity"]["0.75"]
    assert uncertainty["defined_replicates"] == 200
    assert uncertainty["lower"] is None
    assert uncertainty["reason"] == "original_metric_undefined_or_unreachable"
    paired = stratum["matched_sensitivity_comparisons"]["proposed"]["0.75"]
    assert paired["status"] == "not_comparable"
    assert paired["uncertainty"]["defined_replicates"] == 200
    assert paired["uncertainty"]["lower"] is None
    assert paired["uncertainty"]["reason"] == "original_metric_undefined_or_unreachable"


def test_simulated_slots_precision_and_coverage_do_not_claim_completed_assays():
    endpoints = interval(list(range(200)), 3)
    assert (endpoints["lower"], endpoints["upper"]) == (4, 194)
    result = report()
    assert result["format_version"] == "3"
    assert result["assay_count_basis"] == "simulated_testing_slots_including_currently_untested"
    assert result["completed_assays"] is None
    assert result["completed_assays_reason"] == "not_recorded_in_evaluation_inputs"
    methods = result["strata"][0]["methods"]
    proposed = methods["proposed"]["at_budget"]
    baseline = methods["baseline"]["at_budget"]
    assert proposed["assays_used"] == baseline["assays_used"] == 4
    assert proposed["untested"] == proposed["inconclusive"] == 1
    assert (proposed["precision_conclusive"], proposed["conclusive_outcome_coverage"]) == (1, 0.5)
    assert (baseline["precision_conclusive"], baseline["conclusive_outcome_coverage"]) == (0.5, 1)
    curve = methods["proposed"]["precision_recall"]
    assert curve[0]["conclusive_outcome_coverage"] is None
    assert curve[1]["precision_conclusive"] is None
    assert curve[1]["conclusive_outcome_coverage"] == 0
    assert curve[2]["precision_conclusive"] == 1
    assert curve[2]["conclusive_outcome_coverage"] == 0.5
    for method in methods.values():
        for prefix in method["precision_recall"]:
            slots = prefix["assays_used"]
            assert slots == prefix["conclusive_count"] + prefix["inconclusive"] + prefix["untested"]
            assert prefix["conclusive_outcome_coverage"] == (prefix["conclusive_count"] / slots if slots else None)


@pytest.mark.parametrize("numeric_literal", ["0.500000000000000000000000000001", "5e-1", "0.5"])
def test_cli_rejects_inexact_json_targets_but_preserves_decimal_strings(tmp_path, capsys, numeric_literal):
    import json
    protocol = read_json(FIXTURE / "protocol.json")
    target = "0.500000000000000000000000000001"
    protocol["targets"] = [target]
    text = json.dumps(protocol)
    path = tmp_path / "protocol.json"
    path.write_text(text.replace('"' + target + '"', numeric_literal))
    command = ["evaluate"]
    for role in ("candidates", "labels", "methods", "splits"):
        command += ["--" + role, str(FIXTURE / (role + ".json"))]
    command += ["--protocol", str(path)]
    output = tmp_path / "numeric"
    with pytest.raises(SystemExit) as exc:
        main(command + ["--output", str(output)])
    assert exc.value.code == 2
    assert "decimal string" in capsys.readouterr().err
    assert not output.exists()
    path.write_text(text)
    output = tmp_path / "string"
    assert main(command + ["--output", str(output)]) == 0
    matched = read_json(output / "evaluation.json")["strata"][0]["methods"]["proposed"]["matched_sensitivity"][target]
    assert matched == {"required_positives": 3, "status": "unreachable", "assays": None}


def test_cli_rejects_fractional_json_scores_and_split_fractions(tmp_path, capsys):
    import json
    methods = read_json(FIXTURE / "methods.json")
    methods[0]["scores"][0]["score"] = "1.000000000000000000000000000001"
    path = tmp_path / "methods.json"
    path.write_text(json.dumps(methods).replace('"1.000000000000000000000000000001"', '1.000000000000000000000000000001'))
    command = ["evaluate"]
    for role in ("candidates", "labels", "protocol", "splits"):
        command += ["--" + role, str(FIXTURE / (role + ".json"))]
    with pytest.raises(SystemExit) as exc:
        main(command + ["--methods", str(path), "--output", str(tmp_path / "score-output")])
    assert exc.value.code == 2
    assert "score must use an integer or decimal string" in capsys.readouterr().err
    assert not (tmp_path / "score-output").exists()
    split = read_json(FIXTURE / "splits.json")
    split["training_fraction"] = 0.0
    with pytest.raises(PopulationError, match="training fraction must use an integer or decimal string"):
        validate_splits(read_json(FIXTURE / "candidates.json"), split)
    # Runtime is descriptive metadata; fractional JSON values remain allowed there.
    values = list(inputs())
    values[2][0]["runtime_seconds"] = 0.125
    assert report(values)["strata"][0]["methods"]["proposed"]["declared_method_runtime_seconds_total"] == 0.125
