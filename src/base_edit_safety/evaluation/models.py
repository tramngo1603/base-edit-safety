"""Data objects emitted by the validation metrics runner."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Literal


MetricStatus = Literal["pass", "fail", "not_computable", "not_applicable"]


@dataclass(frozen=True)
class MetricResult:
    """Store one metric for one dataset, assay, cell type, editor, and source-overlap group."""

    metric_id: str
    dataset_id: str
    assay: str
    cell_type: str
    status: MetricStatus
    value: float | None
    threshold: dict[str, Any]
    numerator: int | None = None
    denominator: int | None = None
    editor_class: str | None = None
    leakage_stratum: str | None = None
    components: dict[str, Any] = field(default_factory=dict)
    reason: str | None = None

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class FailureRecord:
    """Store one protocol-section-5 failure for a record or combined metric."""

    metric_id: str
    dataset_id: str
    assay: str
    cell_type: str
    failure_code: str
    subject: str
    disposition: str
    human_note: str
    editor_class: str | None = None
    leakage_stratum: str | None = None

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class EvaluationResult:
    """Store complete M1--M9 output in fields that can be written as JSON."""

    metrics: tuple[MetricResult, ...]
    failures: tuple[FailureRecord, ...]
    leakage_strata: dict[str, Any] = field(default_factory=dict)

    def as_dict(self) -> dict[str, Any]:
        return {
            "schema": "base_edit_safety/evaluation-results/v1",
            "metrics": [metric.as_dict() for metric in self.metrics],
            "failures": [failure.as_dict() for failure in self.failures],
            "leakage_strata": self.leakage_strata,
        }
