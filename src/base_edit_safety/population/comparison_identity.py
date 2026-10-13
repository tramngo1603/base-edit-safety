"""Content identities for the comparative records, separate from V3 identities."""
from __future__ import annotations

from decimal import Decimal
import hashlib
import json
import math
import unicodedata

from .common import PopulationError


def _normalized(value):
    if isinstance(value, str):
        return unicodedata.normalize("NFC", value)
    if type(value) is dict:
        result = {}
        for key, item in value.items():
            if not isinstance(key, str):
                raise PopulationError("Canonical JSON keys must be strings")
            key = unicodedata.normalize("NFC", key)
            if key in result:
                raise PopulationError("Object keys collide after Unicode normalization")
            result[key] = _normalized(item)
        return result
    if type(value) is list:
        return [_normalized(item) for item in value]
    if type(value) in (type(None), bool, int):
        return value
    if type(value) is float and math.isfinite(value):
        return value
    raise PopulationError("Canonical JSON requires finite JSON values")


def _encode(value):
    if value is None:
        return "null"
    if type(value) is bool:
        return "true" if value else "false"
    if type(value) is int:
        return str(value)
    if type(value) is float:
        # Equal finite numeric values have one spelling, including signed zero.
        text = format(Decimal(str(value)), "f")
        text = text.rstrip("0").rstrip(".") if "." in text else text
        return "0" if text in ("-0", "") else text
    if isinstance(value, str):
        return json.dumps(value, ensure_ascii=False)
    if isinstance(value, list):
        return "[" + ",".join(_encode(item) for item in value) + "]"
    return "{" + ",".join(_encode(key) + ":" + _encode(value[key]) for key in sorted(value)) + "}"


def _standard_encoding_eligible(value):
    # Keep decimal numbers and deep structures on the original encoding path.
    # Iterators use memory proportional to depth, even for very wide arrays.
    stack = [iter((value,))]
    while stack:
        try:
            item = next(stack[-1])
        except StopIteration:
            stack.pop()
            continue
        if type(item) is float:
            return False
        if type(item) in (dict, list):
            if len(stack) > 64:
                return False
            stack.append(iter(item.values() if type(item) is dict else item))
    return True


def canonical_bytes(value, *, set_fields=()):
    """Encode UTF-8 JSON; sort only named top-level mathematical set fields.

    Paths, alignment columns, ALT indexes and replicates keep their supplied order.
    Nested set fields are declared by the record schemas and handled before this
    function is called.
    """
    value = _normalized(value)
    for key in set_fields:
        if key in value:
            items = value[key]
            if not isinstance(items, list):
                raise PopulationError("A declared set field must be an array")
            encoded = [_encode(item) for item in items]
            if len(set(encoded)) != len(encoded):
                raise PopulationError("Duplicate value in a declared set")
            value[key] = [item for _, item in sorted(zip(encoded, items), key=lambda pair: pair[0])]
    if _standard_encoding_eligible(value):
        return json.dumps(value, ensure_ascii=False, sort_keys=True,
                          separators=(",", ":"), allow_nan=False).encode("utf-8")
    return _encode(value).encode("utf-8")


def content_id(record_type, payload, version="1", *, set_fields=()):
    """Hash both the namespace and payload; never reuse a V3 hash namespace."""
    if not record_type.startswith("comparison_") or (version != "1" and
            not (record_type == "comparison_observation_revision" and version == "2")):
        raise PopulationError("Unknown comparative identity namespace/version")
    namespace = f"{record_type}:v{version}"
    data = namespace.encode("ascii") + b"\0" + canonical_bytes(payload, set_fields=set_fields)
    return namespace + ":" + hashlib.sha256(data).hexdigest()
