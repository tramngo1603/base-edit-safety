"""Conclusive editing evidence requires an explicit, consistent adjudication."""
from decimal import Decimal

from ..measurement_units import editing_detection_limit

OUTCOMES = ("observed_positive", "tested_no_signal", "inconclusive", "untested")


def outcome_reasons(observation):
    """Return reasons to withhold credit without modifying any source value.

    Version-1 records remain readable history, even if they store old credit.
    Source text is not a classifier: only exact four-state labels are checked
    for a contradiction with the separate adjudication.
    """
    outcome = observation.get("adjudicated_outcome")
    if observation.get("schema_version") == "1" or outcome not in OUTCOMES:
        return ["legacy_outcome_requires_adjudication"]
    reasons = []
    if not isinstance(observation.get("adjudication_reason"), str) or not observation["adjudication_reason"].strip():
        reasons.append("outcome_adjudication_reason_missing")
    if outcome in {"inconclusive", "untested"}:
        reasons.append("outcome_" + outcome)
    reported = observation["source_reported_outcome"]["value"]
    if isinstance(reported, str) and reported.strip().casefold() in OUTCOMES \
            and reported.strip().casefold() != outcome:
        reasons.append("adjudication_contradicts_source_outcome")
    value = observation["normalized_value"]["value"]
    number = None if value is None else Decimal(str(value))
    if number is not None and observation["normalized_unit"] == "fraction" and not 0 <= number <= 1:
        reasons.append("measured_fraction_outside_bounds")
    if outcome == "observed_positive" and number is not None and number <= 0:
        reasons.append("positive_outcome_without_positive_measurement")
    if outcome in {"observed_positive", "tested_no_signal"}:
        threshold = _detection_threshold(observation)
        if outcome == "tested_no_signal":
            if threshold is None:
                reasons.append("no_signal_without_usable_detection_limit")
            elif number is not None and not 0 <= number <= threshold:
                reasons.append("no_signal_measurement_outside_detection_limit")
        elif threshold is not None and number is not None and number <= threshold:
            reasons.append("positive_measurement_not_above_detection_limit")
    return reasons


def _detection_threshold(observation):
    limit = observation["detection_limit"]
    amount = limit["value"]["value"]
    if observation["normalized_unit"] == "fraction":
        checked = editing_detection_limit({"value": amount, "unit": limit["unit"]})
        if checked["status"] != "accepted":
            return None
        threshold = Decimal(str(amount))
        if limit["unit"].strip().casefold() in {"percent", "percentage", "%"}:
            threshold /= 100
        return threshold
    if type(amount) not in (int, float) or limit["unit"] != observation["normalized_unit"]:
        return None
    threshold = Decimal(str(amount))
    return threshold if threshold.is_finite() and threshold > 0 else None
