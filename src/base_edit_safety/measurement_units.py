"""Normalize reported limits for measurements of editing frequency."""
from __future__ import annotations

import math
from decimal import Decimal


def editing_detection_limit(reported, maximum_fraction=None):
    """Keep the reported value intact and assess its quantity, units and bounds."""
    result = {"quantity": "editing_fraction", "normalized_fraction": None,
              "maximum_fraction": maximum_fraction, "status": "unresolved", "reason": None}
    if not isinstance(reported, dict) or type(reported.get("value")) not in (int, float) \
            or not math.isfinite(reported["value"]) or reported["value"] <= 0 \
            or not isinstance(reported.get("unit"), str) or not reported["unit"].strip():
        result["reason"] = "numeric_detection_limit_missing"
        return result
    unit = reported["unit"].strip().casefold()
    scales = {"percent": 100, "percentage": 100, "%": 100, "fraction": 1, "proportion": 1}
    if unit not in scales:
        result["reason"] = "detection_limit_unit_unsupported"
        return result
    # Shift the recorded decimal point exactly, without rounding at a boundary.
    sign, digits, exponent = Decimal(str(reported["value"])).as_tuple()
    fraction = Decimal((sign, digits, exponent - (2 if scales[unit] == 100 else 0)))
    if not 0 < fraction <= 1:
        result.update(status="invalid", reason="detection_limit_out_of_bounds")
        return result
    result["normalized_fraction"] = float(fraction)
    if result["normalized_fraction"] == 0:
        result.update(normalized_fraction=None, reason="detection_limit_below_numeric_range")
        return result
    if maximum_fraction is not None and fraction > Decimal(str(maximum_fraction)):
        result.update(status="insufficient", reason="detection_limit_exceeds_study_maximum")
    else:
        result["status"] = "accepted"
    return result
