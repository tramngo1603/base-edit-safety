"""Assign one of three evidence levels to a candidate.

Each candidate has exactly one state. It remains `predicted_only` unless it matches an HBG registry
entry (`analog_supported`) or the caller attaches a measurement from the current run
(`empirically_measured`). A prediction cannot substitute for a missing measurement.
"""
from __future__ import annotations

PREDICTED_ONLY = "predicted_only"
ANALOG_SUPPORTED = "analog_supported"
EMPIRICALLY_MEASURED = "empirically_measured"
EVIDENCE_TIERS = (PREDICTED_ONLY, ANALOG_SUPPORTED, EMPIRICALLY_MEASURED)

# Registry relationships that count as analog support (see hbg_registry.relationship).
_ANALOG_RELATIONSHIPS = frozenset({"exact_edit", "same_motif", "nearby"})


def evidence_tier(registry_relationship: str | None, empirical_measurements) -> str:
    """Pick the tier for one candidate.

    Prefer a measurement from this run, then a registry match, then prediction-only evidence.
    """
    if empirical_measurements:
        return EMPIRICALLY_MEASURED
    if registry_relationship in _ANALOG_RELATIONSHIPS:
        return ANALOG_SUPPORTED
    return PREDICTED_ONLY
