"""Fixed numeric defaults from protocol section 4.

The repository protocol is still a template and has not completed its lock procedure. These values
are therefore used for synthetic tests and become research-data thresholds only after the
checks in `base_edit_safety.evaluation.safety` pass.
No CLI accepts threshold overrides.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from types import MappingProxyType
from typing import Any


@dataclass(frozen=True)
class LockedThresholds:
    m1_guide_seq_class_recall: float = 0.90
    m1_biochemical_recall: float = 0.80
    m2_rank_quality: float = 0.80
    m2_top_k: int = 100
    m3_spearman_rho: float = 0.60
    m4_spearman_rho: float = 0.60
    m5_sensitivity: float = 1.0
    m6_sensitivity: float = 1.0
    m6_specificity: float = 0.95
    m7_sensitivity: float = 1.0
    m8_accuracy: float = 0.95
    m8_audit_sample_size: int = 20
    m9_coverage: float = 1.0

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


LOCKED_THRESHOLDS = LockedThresholds()

# Dataset applicability is part of §4 and is locked with the numeric values.
METRIC_DATASETS = MappingProxyType(
    {
        "M1": frozenset({"D1", "D2", "D3"}),
        "M2": frozenset({"D1", "D2", "D3"}),
        "M3": frozenset({"D4", "D5"}),
        "M4": frozenset({"D4", "D5"}),
    }
)


def threshold_for_m1(assay: str, assay_class: str | None = None) -> tuple[float, str]:
    """Return the assay-sensitive M1 threshold without guessing unknown assays."""
    normalized_class = (assay_class or "").strip().lower().replace("-", "_")
    if normalized_class in {"guide_seq", "guide_seq_class", "cellular"}:
        return LOCKED_THRESHOLDS.m1_guide_seq_class_recall, "GUIDE-seq-class"
    if normalized_class in {"biochemical", "cell_free"}:
        return LOCKED_THRESHOLDS.m1_biochemical_recall, "biochemical"

    normalized_assay = assay.lower().replace("_", "-")
    if "guide-seq" in normalized_assay and "change-seq" not in normalized_assay:
        return LOCKED_THRESHOLDS.m1_guide_seq_class_recall, "GUIDE-seq-class"
    if any(
        token in normalized_assay
        for token in ("change-seq", "changeseq", "digenome", "circle-seq", "site-seq")
    ):
        return LOCKED_THRESHOLDS.m1_biochemical_recall, "biochemical"
    raise ValueError(
        f"M1 assay class is unknown for {assay!r}; provide assay_class="
        "'guide_seq_class' or 'biochemical' in the ingested table"
    )
