"""Check protocol-v1 M1--M9 calculations against fixed expected results."""

import json
from pathlib import Path

import pytest

from base_edit_safety.evaluation import evaluate, spearman_rho
from base_edit_safety.evaluation.inputs import LoadedInputs, load_all


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "tests" / "fixtures" / "evaluation"
OVERLAP = ROOT / "validation" / "han_overlap_map.json"


def _tuning_inputs():
    return load_all(
        [FIXTURES / "tuning_measurements.tsv"],
        [FIXTURES / "tuning_pipeline.json"],
        FIXTURES / "tuning_golden.json",
    )


def _metric(result, metric_id, dataset_id, *, leakage=None):
    matches = [
        metric for metric in result.metrics
        if metric.metric_id == metric_id
        and metric.dataset_id == dataset_id
        and metric.leakage_stratum == leakage
        and (metric.editor_class is not None or metric_id not in {"M3", "M4"})
    ]
    assert len(matches) == 1, matches
    return matches[0]


def test_spearman_known_answers_and_ties():
    assert spearman_rho([1, 2, 4, 3], [10, 20, 30, 40]) == pytest.approx(0.8)
    assert spearman_rho([1, 2, 2, 4], [10, 20, 20, 40]) == pytest.approx(1.0)
    assert spearman_rho([1, 1, 1], [1, 2, 3]) is None
    assert spearman_rho([1], [1]) is None
    with pytest.raises(ValueError, match="equal length"):
        spearman_rho([1, 2], [1])


def test_five_measured_four_predicted_gives_point_eight_recall_and_rank_quality():
    result = evaluate(_tuning_inputs())
    m1 = _metric(result, "M1", "D1")
    assert (m1.numerator, m1.denominator, m1.value, m1.status) == (4, 5, 0.8, "pass")
    assert m1.threshold["value"] == 0.8  # CHANGE-seq-BE is the biochemical-assay class.

    m2 = _metric(result, "M2", "D1")
    assert (m2.numerator, m2.denominator, m2.value, m2.status) == (4, 5, 0.8, "pass")
    assert m2.components["candidate_precision"] is None
    assert m2.components["candidate_precision_status"] == "not_computable"

    missed = [failure for failure in result.failures if failure.subject == "s5"]
    assert {failure.metric_id for failure in missed} == {"M1", "M2"}
    assert {failure.failure_code for failure in missed} == {"F5"}


def test_bystander_outcome_and_golden_metrics_have_known_answers():
    result = evaluate(_tuning_inputs())
    assert _metric(result, "M3", "D4").value == pytest.approx(0.8)
    assert _metric(result, "M3", "D4").status == "pass"
    assert _metric(result, "M4", "D4").value == pytest.approx(1.0)
    for metric_id in ("M5", "M6", "M7", "M8", "M9"):
        metric = _metric(result, metric_id, "D4")
        assert metric.status == "pass", metric
        assert metric.value == pytest.approx(1.0)
    assert _metric(result, "M8", "D4").denominator == 20
    assert _metric(result, "M8", "D4").components["sample_size_satisfied"] is True
    m6 = _metric(result, "M6", "D4")
    assert m6.components["sensitivity"] == 1.0
    assert m6.components["specificity"] == 1.0


def test_every_dataset_stratum_reports_all_nine_metric_ids_including_na():
    result = evaluate(_tuning_inputs())
    for dataset_id in {metric.dataset_id for metric in result.metrics}:
        assert {metric.metric_id for metric in result.metrics if metric.dataset_id == dataset_id} == {
            f"M{number}" for number in range(1, 10)
        }


def test_held_out_han_map_produces_clean_and_partially_entangled_results():
    inputs = load_all(
        [FIXTURES / "held_out_measurements.tsv"],
        [FIXTURES / "held_out_pipeline.json"],
        FIXTURES / "held_out_golden.json",
    )
    result = evaluate(inputs, han_overlap_map=json.loads(OVERLAP.read_text()))
    aggregate = _metric(result, "M3", "D5")
    clean = _metric(result, "M3", "D5", leakage="clean")
    entangled = _metric(result, "M3", "D5", leakage="partially_entangled")
    assert aggregate.value == pytest.approx(0.8)
    assert clean.value == pytest.approx(-1.0)
    assert clean.status == "fail"
    assert entangled.value == pytest.approx(1.0)
    assert entangled.status == "pass"
    # M9 accepts either a complete record or a recorded error. A missing record fails the metric.
    assert _metric(result, "M9", "D5").value == pytest.approx(1.0)
    assert _metric(result, "M9", "D5").status == "pass"
    counts = next(iter(result.leakage_strata.values()))
    assert counts["clean"] == 2
    assert counts["partially_entangled"] == 2


def test_incomplete_candidate_artifact_is_never_misreported_as_low_recall():
    source = _tuning_inputs()
    predictions = tuple(
        {**row, "candidate_list_complete": False}
        if row["record_type"] == "off_target_candidate" else row
        for row in source.predictions
    )
    inputs = LoadedInputs(
        source.measurements, predictions, source.golden,
        dataset_metadata={}, tool_versions=source.tool_versions,
    )
    result = evaluate(inputs)
    assert _metric(result, "M1", "D1").status == "not_computable"
    assert _metric(result, "M2", "D1").status == "not_computable"
    assert "incomplete candidate list" in _metric(result, "M1", "D1").reason
