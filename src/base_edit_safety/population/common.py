"""Strict file reading, content identities, and protected output writes."""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import tempfile


class PopulationError(ValueError):
    """An input does not meet the declared contract."""


class ResourceLimit(PopulationError):
    """A bounded execution cannot finish within its declared limits."""


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False, allow_nan=False).encode("utf-8")


def digest(value):
    return hashlib.sha256(canonical(value)).hexdigest()


def identity(prefix, value):
    return prefix + ":" + digest(value)


def checksum(path, *, check_limits=lambda: None):
    """Return a file's SHA-256, checking resource limits before each block."""
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while True:
            check_limits()
            block = handle.read(1048576)
            if not block:
                break
            h.update(block)
    return h.hexdigest()


def save_json(path, value):
    """Write indented JSON to a local file with owner-only permissions."""
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    Path(path).chmod(0o600)


def known(value):
    """Return a supplied value with no missing-value reason."""
    return {"value": value, "reason": None}


def unknown(reason):
    """Return a null value with its unresolved or missing-value reason."""
    return {"value": None, "reason": reason}


def _pairs(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise PopulationError("Duplicate JSON key")
        result[key] = value
    return result


def _constant(value):
    raise PopulationError("Non-finite JSON number")


def _float(value):
    number = float(value)
    if not math.isfinite(number):
        raise PopulationError("Non-finite JSON number")
    return number


def parse_json(text):
    try:
        return json.loads(text, object_pairs_hook=_pairs, parse_constant=_constant, parse_float=_float)
    except (ValueError, TypeError) as exc:
        raise PopulationError("Invalid JSON: " + str(exc)) from exc


def read_json(path):
    return parse_json(Path(path).read_text(encoding="utf-8"))


def read_jsonl(path):
    return [parse_json(line) for line in Path(path).read_text(encoding="utf-8").splitlines()
            if line.strip()]


def read_tsv(path, columns):
    with Path(path).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or set(reader.fieldnames) != set(columns) \
                or len(reader.fieldnames) != len(columns):
            raise PopulationError("TSV columns do not match the required contract")
        rows = list(reader)
    if any(None in row or any(v is None for v in row.values()) for row in rows):
        raise PopulationError("TSV row width does not match its header")
    return rows


def fields(value, required, optional=(), label="object"):
    if not isinstance(value, dict):
        raise PopulationError(label + " must be an object")
    missing = set(required) - value.keys()
    extra = value.keys() - set(required) - set(optional)
    if missing or extra:
        raise PopulationError(f"{label}: missing fields {sorted(missing)}; unexpected fields {sorted(extra)}")


def integer(value, label, minimum=0, maximum=None):
    if type(value) is not int or value < minimum or (maximum is not None and value > maximum):
        raise PopulationError(label + " is outside its integer bounds")
    return value


def nonempty(value, label):
    if not isinstance(value, str) or not value.strip() or any(c in value for c in "\t\n\r"):
        raise PopulationError(label + " must be a nonempty single-line string")
    return value


def unique(rows, key="id"):
    result = {}
    for row in rows:
        item = nonempty(row[key], "identifier")
        if item in result:
            raise PopulationError("Duplicate identifier")
        result[item] = row
    return result


def record(kind, run_id, **values):
    return {"record_type": kind, "schema_version": "1", "run_id": run_id,
            "provenance": [run_id], **values}


def json_text(value):
    return json.dumps(value, sort_keys=True, indent=2, ensure_ascii=False, allow_nan=False) + "\n"


def jsonl_text(rows):
    return "".join(canonical(row).decode("utf-8") + "\n" for row in rows)


def check_output(path, inputs=(), overwrite=False):
    target = Path(path).resolve()
    for source in inputs:
        source = Path(source).resolve()
        if source == target or target in source.parents or source in target.parents:
            raise PopulationError("Input and output paths overlap")
    if target.exists() and not overwrite:
        raise PopulationError("Output exists; use --overwrite to replace it")
    if Path(path).is_symlink():
        raise PopulationError("Output cannot be a symbolic link")


def write_bundle(directory, files, *, inputs=(), overwrite=False):
    """Commit a complete output directory after every file has been prepared."""
    target = Path(directory).resolve()
    check_output(directory, inputs, overwrite)
    if target.exists() and not target.is_dir():
        raise PopulationError("Output directory is an existing file")
    target.parent.mkdir(parents=True, exist_ok=True)
    temp = Path(tempfile.mkdtemp(prefix=".bes-", dir=target.parent))
    backup = None
    try:
        for name, content in files.items():
            if Path(name).name != name:
                raise PopulationError("Output artifact must be a simple filename")
            p = temp / name
            p.write_text(content, encoding="utf-8")
            p.chmod(0o600)
        if target.exists():
            backup = Path(tempfile.mkdtemp(prefix=".bes-backup-", dir=target.parent))
            backup.rmdir()
            os.replace(target, backup)
        try:
            os.replace(temp, target)
        except OSError:
            if backup is not None:
                os.replace(backup, target)
                backup = None
            raise
    finally:
        if temp.exists():
            shutil.rmtree(temp)
        if backup is not None:
            shutil.rmtree(backup)


def write_file(path, text, *, inputs=(), overwrite=False):
    check_output(path, inputs, overwrite)
    target = Path(path).resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=".bes-", dir=target.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            handle.write(text)
        os.replace(name, target)
    finally:
        if os.path.exists(name):
            os.unlink(name)
