"""Check the paired recall CLI with hand-calculated mock inputs."""
from pathlib import Path

from base_edit_safety.population import comparison_evaluation
from base_edit_safety.population.common import read_json

FIXTURE = Path(__file__).parent / "fixtures/comparison/recall-improvement"
HAND_DIFFERENCE = 1 / 3

def evaluate_cli(case, output):
    command = ["evaluate"]
    for role in ("candidates", "labels", "methods", "protocol", "splits"):
        command += ["--" + role, str(FIXTURE / case / (role + ".json"))]
    assert comparison_evaluation.main(command + ["--output", str(output)]) == 0
    return read_json(output / "evaluation.json")


def test_hand_calculated_signed_recall_difference_survives_the_cli(tmp_path):
    report = evaluate_cli("paired-difference", tmp_path / "evaluation")
    assert report["format_version"] == "3"
    (stratum,) = report["strata"]
    comparison = stratum["recall_at_budget_comparisons"]["proposed"]
    assert comparison["baseline_method"] == "baseline"
    assert comparison["status"] == "comparable"
    assert comparison["budget"] == 1
    assert comparison["eligible_positives"] == 3
    assert comparison["method_recall"] == HAND_DIFFERENCE
    assert comparison["baseline_recall"] == 0
    assert comparison["recall_difference_method_minus_baseline"] == HAND_DIFFERENCE
    uncertainty = comparison["uncertainty"]
    assert uncertainty["reason"] is None
    assert uncertainty["defined_replicates"] == uncertainty["replicates"]
    # The fixed seed includes both possible differences in its percentile bounds.
    assert (uncertainty["lower"], uncertainty["upper"]) == (0, HAND_DIFFERENCE)
    assert stratum["matched_sensitivity_comparisons"]["proposed"]["1"]["baseline_method"] == "baseline"


def test_zero_positive_denominator_reports_null_difference_not_zero_improvement(tmp_path):
    report = evaluate_cli("zero-denominator", tmp_path / "evaluation")
    (stratum,) = report["strata"]
    comparison = stratum["recall_at_budget_comparisons"]["proposed"]
    assert comparison["status"] == "not_comparable"
    assert comparison["eligible_positives"] == 0
    assert comparison["method_recall"] is None
    assert comparison["baseline_recall"] is None
    assert comparison["recall_difference_method_minus_baseline"] is None
    assert comparison["uncertainty"]["reason"] == "original_metric_undefined_or_unreachable"
    assert comparison["uncertainty"]["lower"] is None and comparison["uncertainty"]["upper"] is None
