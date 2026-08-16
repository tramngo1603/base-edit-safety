"""Protocol-v1 validation metrics and controlled evaluation runners.

The public API has two entry points. `evaluate` calculates M1--M9 from values already in memory.
`runner` handles file access, checks dataset roles, writes reports, and enforces the single permitted
held-out run.
"""

from .metrics import evaluate, spearman_rho
from .models import EvaluationResult, MetricResult
from .thresholds import LOCKED_THRESHOLDS, LockedThresholds

__all__ = [
    "EvaluationResult",
    "LOCKED_THRESHOLDS",
    "LockedThresholds",
    "MetricResult",
    "evaluate",
    "spearman_rho",
]
