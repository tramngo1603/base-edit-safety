"""How much evidence stands behind a candidate.

Every candidate lands in exactly one of three states. A missing measurement never quietly becomes a
model number: a candidate stays `predicted_only` unless it matches the HBG registry
(`analog_supported`) or someone attached real measurements from this run (`empirically_measured`).
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

    Real measurements win. Failing that, a registry match. Failing that, predicted only.
    """
    if empirical_measurements:
        return EMPIRICALLY_MEASURED
    if registry_relationship in _ANALOG_RELATIONSHIPS:
        return ANALOG_SUPPORTED
    return PREDICTED_ONLY
