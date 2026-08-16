"""Calculate validation-protocol metrics M1--M9 without file or network access."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
import math
from typing import Any, Mapping, Sequence

from .inputs import LoadedInputs
from .models import EvaluationResult, FailureRecord, MetricResult
from .thresholds import LOCKED_THRESHOLDS, METRIC_DATASETS, LockedThresholds, threshold_for_m1


_FAILURE_CODES = frozenset({"F1", "F2", "F3", "F4", "F5", "F6", "F7"})


@dataclass(frozen=True, order=True)
class _Context:
    dataset_id: str
    assay: str
    cell_type: str


def _ctx(row: Mapping[str, Any]) -> _Context:
    return _Context(
        str(row.get("dataset_id", "")).strip(),
        str(row.get("assay", "not_reported") or "not_reported").strip(),
        str(row.get("cell_type", "not_reported") or "not_reported").strip(),
    )


def _matches_context(row: Mapping[str, Any], context: _Context) -> bool:
    return all(
        not str(row.get(field, "")).strip() or str(row.get(field)).strip() == expected
        for field, expected in (
            ("dataset_id", context.dataset_id),
            ("assay", context.assay),
            ("cell_type", context.cell_type),
        )
    )


def _as_bool(value: Any, *, label: str) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return value != 0
    normalized = str(value).strip().lower()
    if normalized in {"true", "yes", "y", "1", "pass", "positive", "detected"}:
        return True
    if normalized in {"false", "no", "n", "0", "0.0", "fail", "negative", "not_detected"}:
        return False
    try:
        return float(normalized) != 0
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} is not boolean-like: {value!r}") from exc


def _optional_bool(value: Any, *, label: str) -> bool | None:
    if value is None or str(value).strip() == "":
        return None
    return _as_bool(value, label=label)


def _number(value: Any, *, label: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} is not numeric: {value!r}") from exc
    if not math.isfinite(number):
        raise ValueError(f"{label} must be finite: {value!r}")
    return number


def _ranks(values: Sequence[float]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    ranks = [0.0] * len(values)
    cursor = 0
    while cursor < len(order):
        end = cursor + 1
        while end < len(order) and values[order[end]] == values[order[cursor]]:
            end += 1
        average_rank = (cursor + 1 + end) / 2.0
        for index in order[cursor:end]:
            ranks[index] = average_rank
        cursor = end
    return ranks


def spearman_rho(predicted: Sequence[float], measured: Sequence[float]) -> float | None:
    """Spearman's rho with average ranks for ties; undefined cases return ``None``."""
    if len(predicted) != len(measured):
        raise ValueError("Spearman vectors must have equal length")
    if len(predicted) < 2:
        return None
    x = _ranks([_number(value, label="predicted Spearman value") for value in predicted])
    y = _ranks([_number(value, label="measured Spearman value") for value in measured])
    mean_x = sum(x) / len(x)
    mean_y = sum(y) / len(y)
    dx = [value - mean_x for value in x]
    dy = [value - mean_y for value in y]
    denominator = math.sqrt(sum(value * value for value in dx) * sum(value * value for value in dy))
    if denominator == 0:
        return None
    return sum(left * right for left, right in zip(dx, dy)) / denominator


def _failure(
    row: Mapping[str, Any],
    metric_id: str,
    context: _Context,
    subject: str,
    *,
    default_code: str,
    editor_class: str | None = None,
    leakage_stratum: str | None = None,
) -> FailureRecord:
    code = str(row.get("failure_code") or default_code).strip().upper()
    if code not in _FAILURE_CODES:
        raise ValueError(f"invalid protocol §5 failure code {code!r} for {subject}")
    if code in {"F1", "F2", "F3", "F4", "F7"}:
        default_disposition = "bug: fix and re-run tuning only"
    else:
        default_disposition = "finding: report; do not engineer around after freeze"
    return FailureRecord(
        metric_id=metric_id,
        dataset_id=context.dataset_id,
        assay=context.assay,
        cell_type=context.cell_type,
        editor_class=editor_class,
        leakage_stratum=leakage_stratum,
        failure_code=code,
        subject=subject,
        disposition=str(row.get("failure_disposition") or default_disposition),
        human_note=str(row.get("failure_note") or "Human review of the automatic classification is required."),
    )


def _site_key(row: Mapping[str, Any]) -> tuple[Any, ...]:
    if str(row.get("site_id", "")).strip():
        return ("site_id", str(row["site_id"]).strip())
    chrom = str(row.get("chrom", "")).strip()
    if chrom:
        for field, adjustment in (
            ("start_0based", 0),
            ("start", 0),
            ("pos_1based", -1),
            ("genomic_position_1based", -1),
            ("start_1based", -1),
        ):
            if str(row.get(field, "")).strip():
                position = int(float(row[field])) + adjustment
                # M1 asks whether a genomic site appears anywhere in the candidate list. Strand
                # is useful diagnostic metadata but does not make the same locus a different site.
                return ("coordinate", chrom, position)
    sequence = row.get("offtarget_sequence") or row.get("sequence")
    if str(sequence or "").strip():
        return ("sequence", str(sequence).replace("-", "").upper())
    raise ValueError(f"off-target row has no stable site identity: {dict(row)}")


def _guide_site_key(row: Mapping[str, Any]) -> tuple[Any, ...]:
    guide_id = str(row.get("guide_id", "")).strip()
    if not guide_id:
        raise ValueError(f"off-target row has no guide_id: {dict(row)}")
    return (guide_id, *_site_key(row))


def _point_key(row: Mapping[str, Any]) -> tuple[Any, ...]:
    guide = str(row.get("guide_id", "")).strip()
    editor = str(row.get("editor_class", "")).strip().upper()
    if not guide or not editor:
        raise ValueError(f"concordance row needs guide_id and editor_class: {dict(row)}")
    point = row.get("point_id") or row.get("site_id")
    if str(point or "").strip():
        return (guide, editor, "point", str(point).strip())
    if all(str(row.get(field, "")).strip() for field in ("offset", "ref", "alt")):
        return (
            guide,
            editor,
            "allele",
            -abs(int(float(row["offset"]))),
            str(row["ref"]).upper(),
            str(row["alt"]).upper(),
        )
    return (guide, editor, *_site_key(row))


def _subject(row: Mapping[str, Any], fallback: str) -> str:
    for field in ("site_id", "point_id", "guide_id", "target_id"):
        if str(row.get(field, "")).strip():
            return str(row[field]).strip()
    return fallback


def _metadata_bool(inputs: LoadedInputs, dataset_id: str, field: str) -> bool | None:
    value = inputs.dataset_metadata.get(dataset_id, {}).get(field)
    return _optional_bool(value, label=f"{dataset_id} metadata {field}")


def _candidate_map(
    inputs: LoadedInputs, context: _Context
) -> tuple[dict[tuple[Any, ...], int], bool]:
    candidates = [
        row
        for row in inputs.predictions
        if row.get("record_type") == "off_target_candidate" and _matches_context(row, context)
    ]
    ranks: dict[tuple[Any, ...], int] = {}
    complete = _metadata_bool(inputs, context.dataset_id, "candidate_list_complete")
    for index, row in enumerate(candidates, start=1):
        rank = int(float(row.get("rank") or index))
        if rank < 1:
            raise ValueError(f"candidate rank must be >=1: {row}")
        key = _guide_site_key(row)
        ranks[key] = min(rank, ranks.get(key, rank))
        row_complete = _optional_bool(
            row.get("candidate_list_complete"), label="candidate_list_complete"
        )
        if row_complete is False:
            complete = False
    return ranks, complete is not False


def _metric(
    metric_id: str,
    context: _Context,
    *,
    status: str,
    value: float | None,
    threshold: dict[str, Any],
    numerator: int | None = None,
    denominator: int | None = None,
    editor_class: str | None = None,
    leakage_stratum: str | None = None,
    components: dict[str, Any] | None = None,
    reason: str | None = None,
) -> MetricResult:
    return MetricResult(
        metric_id=metric_id,
        dataset_id=context.dataset_id,
        assay=context.assay,
        cell_type=context.cell_type,
        editor_class=editor_class,
        leakage_stratum=leakage_stratum,
        status=status,  # type: ignore[arg-type]
        value=value,
        threshold=threshold,
        numerator=numerator,
        denominator=denominator,
        components=components or {},
        reason=reason,
    )


def _off_target_metrics(
    inputs: LoadedInputs,
    context: _Context,
    thresholds: LockedThresholds,
) -> tuple[list[MetricResult], list[FailureRecord]]:
    measured = [
        row
        for row in inputs.measurements
        if row.get("record_type") == "off_target_site" and _ctx(row) == context
    ]
    positives = [
        row
        for row in measured
        if _as_bool(row.get("measured_positive", True), label="measured_positive")
    ]
    candidates, complete = _candidate_map(inputs, context)
    threshold_spec = {"operator": ">=", "value": None}
    m1_threshold: float | None = None
    m1_threshold_reason: str | None = None
    try:
        m1_threshold, assay_class = threshold_for_m1(
            context.assay,
            next((str(row.get("assay_class")) for row in measured if row.get("assay_class")), None),
        )
        threshold_spec = {"operator": ">=", "value": m1_threshold, "assay_class": assay_class}
    except ValueError as exc:
        m1_threshold_reason = str(exc)

    if not positives:
        reason = "no measured positive off-target sites in this stratum"
        return (
            [
                _metric("M1", context, status="not_computable", value=None,
                        threshold=threshold_spec, denominator=0, reason=reason),
                _metric(
                    "M2", context, status="not_computable", value=None,
                    threshold={"operator": ">=", "value": thresholds.m2_rank_quality,
                               "top_k": thresholds.m2_top_k},
                    denominator=0, reason=reason,
                ),
            ],
            [],
        )
    if not complete:
        reason = "pipeline output declares an incomplete candidate list"
        return (
            [
                _metric("M1", context, status="not_computable", value=None,
                        threshold=threshold_spec, denominator=len(positives), reason=reason),
                _metric(
                    "M2", context, status="not_computable", value=None,
                    threshold={"operator": ">=", "value": thresholds.m2_rank_quality,
                               "top_k": thresholds.m2_top_k},
                    denominator=len(positives), reason=reason,
                ),
            ],
            [],
        )

    matched = [row for row in positives if _guide_site_key(row) in candidates]
    top = [
        row
        for row in positives
        if candidates.get(_guide_site_key(row), thresholds.m2_top_k + 1) <= thresholds.m2_top_k
    ]
    recall = len(matched) / len(positives)
    rank_quality = len(top) / len(positives)
    failures: list[FailureRecord] = []
    for row in positives:
        key = _guide_site_key(row)
        if key not in candidates:
            failures.append(
                _failure(row, "M1", context, _subject(row, repr(key)), default_code="F6")
            )
        if candidates.get(key, thresholds.m2_top_k + 1) > thresholds.m2_top_k:
            failures.append(
                _failure(row, "M2", context, _subject(row, repr(key)), default_code="F6")
            )

    true_negatives = _metadata_bool(inputs, context.dataset_id, "true_negatives_available")
    row_tn = {
        value
        for row in measured
        if (value := _optional_bool(row.get("true_negatives_available"), label="true_negatives_available"))
        is not None
    }
    if len(row_tn) > 1:
        raise ValueError(f"conflicting true_negatives_available values in {context}")
    if row_tn:
        true_negatives = row_tn.pop()
    top_candidates = {key for key, rank in candidates.items() if rank <= thresholds.m2_top_k}
    positive_keys = {_guide_site_key(row) for row in positives}
    precision = (
        len(top_candidates & positive_keys) / len(top_candidates)
        if true_negatives and top_candidates
        else None
    )
    precision_status = "computed" if precision is not None else "not_computable"
    precision_reason = None
    if precision is None:
        precision_reason = (
            "no candidates in the top 100" if true_negatives else "true negatives unavailable"
        )

    m1_status = (
        "not_computable"
        if m1_threshold is None
        else ("pass" if recall >= m1_threshold else "fail")
    )

    return (
        [
            _metric(
                "M1", context,
                status=m1_status,
                value=recall,
                threshold=threshold_spec,
                numerator=len(matched), denominator=len(positives),
                components={"missed_site_count": len(positives) - len(matched)},
                reason=m1_threshold_reason,
            ),
            _metric(
                "M2", context,
                status="pass" if rank_quality >= thresholds.m2_rank_quality else "fail",
                value=rank_quality,
                threshold={"operator": ">=", "value": thresholds.m2_rank_quality,
                           "top_k": thresholds.m2_top_k},
                numerator=len(top), denominator=len(positives),
                components={
                    "rank_quality_at_100": rank_quality,
                    "candidate_precision": precision,
                    "candidate_precision_status": precision_status,
                    "candidate_precision_reason": precision_reason,
                    "top_candidate_count": len(top_candidates),
                },
            ),
        ],
        failures,
    )


def _prediction_values(
    inputs: LoadedInputs,
    context: _Context,
    record_type: str,
    editor_class: str,
) -> dict[tuple[Any, ...], float]:
    values: dict[tuple[Any, ...], float] = {}
    for row in inputs.predictions:
        if row.get("record_type") != record_type or not _matches_context(row, context):
            continue
        if str(row.get("editor_class", "")).upper() != editor_class:
            continue
        key = _point_key(row)
        value = _number(row.get("predicted_value"), label=f"{record_type} predicted_value")
        if key in values and not math.isclose(values[key], value, rel_tol=0, abs_tol=1e-12):
            raise ValueError(f"conflicting predictions for {key}: {values[key]} vs {value}")
        values[key] = value
    return values


def _leakage_key(row: Mapping[str, Any]) -> tuple[int, str, str, str]:
    missing = [field for field in ("offset", "ref", "alt", "editor_class")
               if not str(row.get(field, "")).strip()]
    if missing:
        raise ValueError(f"held-out bystander row lacks Han-stratification fields {missing}: {row}")
    return (
        -abs(int(float(row["offset"]))),
        str(row["ref"]).upper(),
        str(row["alt"]).upper(),
        str(row["editor_class"]).upper(),
    )


def _correlation_metric(
    inputs: LoadedInputs,
    context: _Context,
    metric_id: str,
    measurement_type: str,
    prediction_type: str,
    threshold: float,
    *,
    overlap_keys: frozenset[tuple[int, str, str, str]] | None,
) -> tuple[list[MetricResult], list[FailureRecord], dict[str, int]]:
    measurements = [
        row for row in inputs.measurements
        if row.get("record_type") == measurement_type and _ctx(row) == context
    ]
    editor_classes = sorted({str(row.get("editor_class", "")).upper() for row in measurements if row.get("editor_class")})
    if not editor_classes:
        empty_results = [
            _metric(metric_id, context, status="not_computable", value=None,
                    threshold={"operator": ">=", "value": threshold},
                    reason=f"no {measurement_type} rows in this stratum")
        ]
        if overlap_keys is not None and metric_id == "M3":
            empty_results.extend(
                _metric(
                    metric_id, context, status="not_computable", value=None,
                    threshold={"operator": ">=", "value": threshold},
                    leakage_stratum=leakage,
                    reason=f"no points in the {leakage} leakage stratum",
                )
                for leakage in ("clean", "partially_entangled")
            )
        return (
            empty_results,
            [],
            {},
        )

    results: list[MetricResult] = []
    failures: list[FailureRecord] = []
    leakage_counts: defaultdict[str, int] = defaultdict(int)
    for editor_class in editor_classes:
        editor_rows = [
            row for row in measurements
            if str(row.get("editor_class", "")).upper() == editor_class
        ]
        predictions = _prediction_values(inputs, context, prediction_type, editor_class)
        partitions: list[tuple[str | None, list[dict[str, Any]]]] = [(None, editor_rows)]
        if overlap_keys is not None and metric_id == "M3":
            clean, entangled = [], []
            for row in editor_rows:
                leakage = (
                    "partially_entangled" if _leakage_key(row) in overlap_keys else "clean"
                )
                leakage_counts[leakage] += 1
                (entangled if leakage == "partially_entangled" else clean).append(row)
            partitions.extend([("clean", clean), ("partially_entangled", entangled)])

        for leakage_stratum, rows in partitions:
            predicted_values: list[float] = []
            measured_values: list[float] = []
            missing: list[dict[str, Any]] = []
            for row in rows:
                key = _point_key(row)
                if key not in predictions:
                    missing.append(row)
                    continue
                predicted_values.append(predictions[key])
                measured_values.append(_number(row.get("measured_value"), label="measured_value"))
            rho = spearman_rho(predicted_values, measured_values)
            if not rows:
                status = "not_computable"
                reason = "no points in this leakage stratum"
            elif missing:
                status = "fail"
                reason = f"{len(missing)} measured point(s) have no prediction"
            elif rho is None:
                status = "not_computable"
                reason = "Spearman rho requires at least two pairs and non-constant ranked vectors"
            else:
                status = "pass" if rho >= threshold else "fail"
                reason = None
            results.append(
                _metric(
                    metric_id, context, status=status, value=rho,
                    threshold={"operator": ">=", "value": threshold},
                    numerator=len(predicted_values), denominator=len(rows),
                    editor_class=editor_class, leakage_stratum=leakage_stratum,
                    components={"paired_points": len(predicted_values),
                                "missing_predictions": len(missing)},
                    reason=reason,
                )
            )
            for row in missing:
                failures.append(
                    _failure(
                        row, metric_id, context, _subject(row, repr(_point_key(row))),
                        default_code="F6", editor_class=editor_class,
                        leakage_stratum=leakage_stratum,
                    )
                )
            if status == "fail" and not missing:
                failures.append(
                    _failure(
                        {}, metric_id, context, f"aggregate rho={rho}", default_code="F6",
                        editor_class=editor_class, leakage_stratum=leakage_stratum,
                    )
                )
    return results, failures, dict(leakage_counts)


def _guide_result(inputs: LoadedInputs, golden: Mapping[str, Any], context: _Context) -> Mapping[str, Any] | None:
    guide_id = str(golden.get("guide_id", "")).strip()
    expected_editor = str(golden.get("editor_class", "")).strip().upper()
    candidates = [
        row for row in inputs.predictions
        if row.get("record_type") == "guide_result"
        and str(row.get("guide_id", "")).strip() == guide_id
        and _matches_context(row, context)
        and (
            not expected_editor
            or not str(row.get("editor_class", "")).strip()
            or str(row.get("editor_class", "")).strip().upper() == expected_editor
        )
    ]
    if not candidates:
        return None
    # Prefer the row with the most filled dataset, assay, and cell-type fields. Equally complete rows
    # must agree.
    candidates.sort(
        key=lambda row: sum(bool(str(row.get(field, "")).strip())
                            for field in ("dataset_id", "assay", "cell_type")),
        reverse=True,
    )
    best_score = sum(bool(str(candidates[0].get(field, "")).strip())
                     for field in ("dataset_id", "assay", "cell_type"))
    best = [row for row in candidates if sum(bool(str(row.get(field, "")).strip())
                                             for field in ("dataset_id", "assay", "cell_type")) == best_score]
    if len(best) > 1 and any(dict(row) != dict(best[0]) for row in best[1:]):
        raise ValueError(f"conflicting guide results for {guide_id} in {context}")
    return best[0]


def _complete_result(result: Mapping[str, Any] | None) -> bool:
    if result is None:
        return False
    status = str(result.get("status", "complete")).strip().lower()
    declared = _optional_bool(result.get("record_complete"), label="record_complete")
    return status in {"complete", "completed", "success"} and declared is not False


def _actual_flag(result: Mapping[str, Any] | None, names: Sequence[str]) -> bool:
    if not _complete_result(result):
        return False
    for name in names:
        if name in result and str(result[name]).strip() != "":
            return _as_bool(result[name], label=name)
    nested = result.get("flags")
    if isinstance(nested, Mapping):
        for name in names:
            if name in nested and str(nested[name]).strip() != "":
                return _as_bool(nested[name], label=f"flags.{name}")
    return False


def _expected_flag(row: Mapping[str, Any], names: Sequence[str], *, default: bool) -> bool:
    for name in names:
        if name in row and str(row[name]).strip() != "":
            return _as_bool(row[name], label=name)
    for container_name in ("expected_flags", "expected"):
        nested = row.get(container_name)
        if isinstance(nested, Mapping):
            for name in names:
                if name in nested and str(nested[name]).strip() != "":
                    return _as_bool(nested[name], label=f"{container_name}.{name}")
    return default


def _golden_metrics(
    inputs: LoadedInputs,
    context: _Context,
    thresholds: LockedThresholds,
    *,
    enforce_audit_sample_size: bool,
    total_audit_count: int,
) -> tuple[list[MetricResult], list[FailureRecord]]:
    golden = [row for row in inputs.golden if _ctx(row) == context]
    if not golden:
        return [], []
    joined = [(row, _guide_result(inputs, row, context)) for row in golden]
    metrics: list[MetricResult] = []
    failures: list[FailureRecord] = []

    dual_positive = [
        item for item in joined
        if _expected_flag(
            item[0], ("expected_dual_engagement", "dual_engagement", "co_targets_paralog"),
            default=False,
        )
    ]
    dual_hits = sum(
        _actual_flag(result, ("dual_engagement_flag", "dual_engagement", "co_targets_paralog"))
        for _row, result in dual_positive
    )
    if dual_positive:
        value = dual_hits / len(dual_positive)
        metrics.append(_metric(
            "M5", context, status="pass" if value >= thresholds.m5_sensitivity else "fail",
            value=value, threshold={"operator": ">=", "value": thresholds.m5_sensitivity,
                                    "requirement": "no false negatives"},
            numerator=dual_hits, denominator=len(dual_positive),
        ))
        for row, result in dual_positive:
            if not _actual_flag(result, ("dual_engagement_flag", "dual_engagement", "co_targets_paralog")):
                failures.append(_failure(row, "M5", context, str(row["guide_id"]), default_code="F6"))
    else:
        metrics.append(_metric(
            "M5", context, status="not_applicable", value=None,
            threshold={"operator": ">=", "value": thresholds.m5_sensitivity},
            reason="golden stratum has no dual-engaging positive guides",
        ))

    deletion_expectations = []
    for item in joined:
        row = item[0]
        dual_expected = item in dual_positive
        deletion_expected = _expected_flag(
            row,
            ("expected_intergenic_deletion_warning", "intergenic_deletion_warning", "deletion_warning"),
            default=dual_expected,
        )
        if dual_expected and not deletion_expected:
            raise ValueError(
                f"golden guide {row['guide_id']} is dual-engaging but explicitly expects no "
                "intergenic-deletion warning"
            )
        deletion_expectations.append((item, deletion_expected))
    deletion_positive = [item for item, expected in deletion_expectations if expected]
    deletion_negative = [item for item, expected in deletion_expectations if not expected]
    deletion_hits = sum(
        _actual_flag(result, ("intergenic_deletion_warning", "deletion_warning"))
        for _row, result in deletion_positive
    )
    deletion_true_negatives = sum(
        not _actual_flag(result, ("intergenic_deletion_warning", "deletion_warning"))
        for _row, result in deletion_negative
    )
    sensitivity = deletion_hits / len(deletion_positive) if deletion_positive else None
    specificity = (
        deletion_true_negatives / len(deletion_negative) if deletion_negative else None
    )
    m6_computable = sensitivity is not None and specificity is not None
    m6_pass = (
        m6_computable
        and sensitivity >= thresholds.m6_sensitivity
        and specificity >= thresholds.m6_specificity
    )
    metrics.append(_metric(
        "M6", context,
        status=("pass" if m6_pass else "fail") if m6_computable else "not_computable",
        value=sensitivity,
        threshold={
            "sensitivity": {"operator": ">=", "value": thresholds.m6_sensitivity},
            "specificity": {"operator": ">=", "value": thresholds.m6_specificity},
        },
        numerator=deletion_hits if deletion_positive else None,
        denominator=len(deletion_positive) if deletion_positive else None,
        components={
            "sensitivity": sensitivity,
            "sensitivity_numerator": deletion_hits,
            "sensitivity_denominator": len(deletion_positive),
            "specificity": specificity,
            "specificity_numerator": deletion_true_negatives,
            "specificity_denominator": len(deletion_negative),
        },
        reason=None if m6_computable else "M6 needs both dual-engaging positives and golden negatives",
    ))
    for row, result in deletion_positive:
        if not _actual_flag(result, ("intergenic_deletion_warning", "deletion_warning")):
            failures.append(_failure(row, "M6", context, str(row["guide_id"]), default_code="F6"))
    for row, result in deletion_negative:
        if _actual_flag(result, ("intergenic_deletion_warning", "deletion_warning")):
            failures.append(_failure(row, "M6", context, str(row["guide_id"]), default_code="F6"))

    chip_positive = [
        item for item in joined
        if _expected_flag(
            item[0],
            ("expected_chip_driver_overlap", "chip_driver_overlap", "ch_driver_overlap"),
            default=False,
        )
    ]
    chip_hits = sum(
        _actual_flag(result, ("chip_driver_overlap_flag", "ch_driver_overlap_flag", "chip_flag"))
        for _row, result in chip_positive
    )
    if chip_positive:
        value = chip_hits / len(chip_positive)
        metrics.append(_metric(
            "M7", context, status="pass" if value >= thresholds.m7_sensitivity else "fail",
            value=value, threshold={"operator": ">=", "value": thresholds.m7_sensitivity,
                                    "requirement": "zero false negatives"},
            numerator=chip_hits, denominator=len(chip_positive),
        ))
        for row, result in chip_positive:
            if not _actual_flag(result, ("chip_driver_overlap_flag", "ch_driver_overlap_flag", "chip_flag")):
                failures.append(_failure(row, "M7", context, str(row["guide_id"]), default_code="F6"))
    else:
        metrics.append(_metric(
            "M7", context, status="not_applicable", value=None,
            threshold={"operator": ">=", "value": thresholds.m7_sensitivity},
            reason="golden stratum has no CHIP-driver positive guides",
        ))

    audits = [item for item in joined if str(item[0].get("expected_evidence_tier", "")).strip()]
    correct = sum(
        _complete_result(result)
        and str(result.get("evidence_tier", "")).strip() == str(row["expected_evidence_tier"]).strip()
        for row, result in audits
    )
    if audits:
        value = correct / len(audits)
        # Protocol M8 samples 20 records across the full manual evidence review, then reports results
        # by dataset group. The 20-record minimum applies to the full sample, not to each group.
        observed_enough = total_audit_count >= thresholds.m8_audit_sample_size
        accepted_sample_size = observed_enough or not enforce_audit_sample_size
        status = "pass" if value >= thresholds.m8_accuracy and accepted_sample_size else "fail"
        metrics.append(_metric(
            "M8", context, status=status, value=value,
            threshold={"operator": ">=", "value": thresholds.m8_accuracy,
                       "sample_size": thresholds.m8_audit_sample_size},
            numerator=correct, denominator=len(audits),
            components={"required_sample_size": thresholds.m8_audit_sample_size,
                        "total_audit_sample_size": total_audit_count,
                        "sample_size_satisfied": observed_enough,
                        "sample_size_enforced": enforce_audit_sample_size},
            reason=(
                None
                if accepted_sample_size
                else f"manual audit has {total_audit_count} total records; protocol requires {thresholds.m8_audit_sample_size}"
            ),
        ))
        for row, result in audits:
            if not _complete_result(result) or str(result.get("evidence_tier", "")).strip() != str(row["expected_evidence_tier"]).strip():
                failures.append(_failure(row, "M8", context, str(row["guide_id"]), default_code="F7"))
        if not accepted_sample_size:
            failures.append(
                _failure(
                    {}, "M8", context,
                    f"manual audit sample size {total_audit_count}/{thresholds.m8_audit_sample_size}",
                    default_code="F7",
                )
            )
    else:
        metrics.append(_metric(
            "M8", context, status="not_computable", value=None,
            threshold={"operator": ">=", "value": thresholds.m8_accuracy,
                       "sample_size": thresholds.m8_audit_sample_size},
            reason="no manually audited evidence-tier expectations in this golden stratum",
        ))
    return metrics, failures


def _m9(
    inputs: LoadedInputs,
    context: _Context,
    thresholds: LockedThresholds,
) -> tuple[MetricResult, list[FailureRecord]]:
    expected: dict[str, Mapping[str, Any]] = {}
    for row in inputs.measurements:
        if _ctx(row) == context and str(row.get("guide_id", "")).strip():
            expected.setdefault(str(row["guide_id"]).strip(), row)
    for row in inputs.golden:
        if _ctx(row) == context:
            expected.setdefault(str(row["guide_id"]).strip(), row)
    configured = inputs.dataset_metadata.get(context.dataset_id, {}).get("registry_guides", [])
    if isinstance(configured, list):
        for guide in configured:
            expected.setdefault(str(guide), {"guide_id": str(guide)})

    if not expected:
        return (
            _metric(
                "M9", context, status="not_computable", value=None,
                threshold={"operator": ">=", "value": thresholds.m9_coverage}, denominator=0,
                reason="no registry or golden-set guides in this stratum",
            ),
            [],
        )
    accounted = 0
    failures = []
    for guide_id, source_row in expected.items():
        result = _guide_result(inputs, {**source_row, "guide_id": guide_id}, context)
        status = str(result.get("status", "complete") if result else "missing").strip().lower()
        error = "" if result is None else str(
            result.get("error") or result.get("error_message") or result.get("logged_error") or ""
        ).strip()
        valid = _complete_result(result) or (
            status in {"error", "failed", "failure"} and bool(error)
        )
        if valid:
            accounted += 1
        else:
            failures.append(_failure(source_row, "M9", context, guide_id, default_code="F7"))
    value = accounted / len(expected)
    return (
        _metric(
            "M9", context, status="pass" if value >= thresholds.m9_coverage else "fail",
            value=value, threshold={"operator": ">=", "value": thresholds.m9_coverage},
            numerator=accounted, denominator=len(expected),
            components={"complete_or_explicit_error": accounted,
                        "silent_failure_count": len(expected) - accounted},
        ),
        failures,
    )


def _placeholder(metric_id: str, context: _Context, thresholds: LockedThresholds) -> MetricResult:
    reasons = {
        "M1": "metric applies only to D1, D2, and D3",
        "M2": "metric applies only to D1, D2, and D3",
        "M3": "metric applies only to D4 and D5",
        "M4": "metric applies only to D4 and D5 and only if v1 emits the outcome",
        "M5": "no golden-set rows in this stratum",
        "M6": "no golden-set rows in this stratum",
        "M7": "no golden-set rows in this stratum",
        "M8": "no golden-set rows in this stratum",
    }
    threshold: dict[str, Any] = {"source": "validation-protocol.md §4"}
    return _metric(
        metric_id, context, status="not_applicable", value=None, threshold=threshold,
        reason=reasons.get(metric_id, "metric has no applicable records in this stratum"),
    )


def evaluate(
    inputs: LoadedInputs,
    *,
    thresholds: LockedThresholds = LOCKED_THRESHOLDS,
    han_overlap_map: Mapping[str, Any] | None = None,
    enforce_audit_sample_size: bool = True,
) -> EvaluationResult:
    """Compute M1--M9 without performing filesystem or external-tool operations."""
    contexts = sorted({_ctx(row) for row in (*inputs.measurements, *inputs.golden)})
    if not contexts:
        raise ValueError("evaluation inputs contain no measurement or golden-set contexts")

    overlap_keys: frozenset[tuple[int, str, str, str]] | None = None
    if han_overlap_map is not None:
        if han_overlap_map.get("schema") != "han_overlap_map/v1":
            raise ValueError("unsupported Han overlap-map schema")
        raw_keys = han_overlap_map.get("partially_entangled_stratum_keys")
        if not isinstance(raw_keys, list):
            raise ValueError("Han overlap map has no partially_entangled_stratum_keys array")
        overlap_keys = frozenset(
            (int(key[0]), str(key[1]).upper(), str(key[2]).upper(), str(key[3]).upper())
            for key in raw_keys
            if isinstance(key, list) and len(key) == 4
        )
        if len(overlap_keys) != len(raw_keys):
            raise ValueError("malformed or duplicate Han overlap-map stratum key")

    metrics: list[MetricResult] = []
    failures: list[FailureRecord] = []
    leakage_summary: dict[str, Any] = {}
    total_audit_count = sum(
        bool(str(row.get("expected_evidence_tier", "")).strip()) for row in inputs.golden
    )
    for context in contexts:
        if context.dataset_id in METRIC_DATASETS["M1"]:
            new_metrics, new_failures = _off_target_metrics(inputs, context, thresholds)
            metrics.extend(new_metrics)
            failures.extend(new_failures)

        if context.dataset_id in METRIC_DATASETS["M3"]:
            new_metrics, new_failures, counts = _correlation_metric(
                inputs, context, "M3", "bystander", "bystander_prediction",
                thresholds.m3_spearman_rho,
                overlap_keys=overlap_keys if context.dataset_id == "D5" else None,
            )
            metrics.extend(new_metrics)
            failures.extend(new_failures)
            if context.dataset_id == "D5" and overlap_keys is not None:
                key = f"{context.dataset_id}|{context.assay}|{context.cell_type}"
                leakage_summary[key] = {
                    "clean": counts.get("clean", 0),
                    "partially_entangled": counts.get("partially_entangled", 0),
                    "overlap_map_key_count": len(overlap_keys),
                }
            new_metrics, new_failures, _unused = _correlation_metric(
                inputs, context, "M4", "on_target_outcome", "on_target_prediction",
                thresholds.m4_spearman_rho, overlap_keys=None,
            )
            metrics.extend(new_metrics)
            failures.extend(new_failures)

        new_metrics, new_failures = _golden_metrics(
            inputs, context, thresholds,
            enforce_audit_sample_size=enforce_audit_sample_size,
            total_audit_count=total_audit_count,
        )
        metrics.extend(new_metrics)
        failures.extend(new_failures)
        m9, m9_failures = _m9(inputs, context, thresholds)
        metrics.append(m9)
        failures.extend(m9_failures)

        present = {metric.metric_id for metric in metrics if (
            metric.dataset_id, metric.assay, metric.cell_type
        ) == (context.dataset_id, context.assay, context.cell_type)}
        for metric_id in (f"M{number}" for number in range(1, 10)):
            if metric_id not in present:
                metrics.append(_placeholder(metric_id, context, thresholds))

    metrics.sort(key=lambda result: (
        result.dataset_id, result.assay, result.cell_type, result.metric_id,
        result.editor_class or "", result.leakage_stratum or "",
    ))
    failures.sort(key=lambda failure: (
        failure.dataset_id, failure.assay, failure.cell_type, failure.metric_id,
        failure.subject, failure.leakage_stratum or "",
    ))
    return EvaluationResult(tuple(metrics), tuple(failures), leakage_summary)
