"""Bind comparison stages to their local sources and direct dependencies."""
from __future__ import annotations

from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path

from .common import (PopulationError, checksum, fields, json_text, jsonl_text,
                     read_json, read_jsonl, write_bundle)
from .provenance import code_identity


STAGE_PARENTS = {
    "search": (),
    "observation_normalization": (),
    "evidence_applicability": ("search", "observation_normalization"),
    "predictions": ("search",),
    "comparison": ("search", "evidence_applicability"),
    "questions": ("comparison",),
    "planning": ("questions", "evidence_applicability"),
    "report": ("comparison", "questions", "planning"),
}

STAGE_POLICIES = {
    "search": ("search",),
    "observation_normalization": ("evidence",),
    "evidence_applicability": ("evidence",),
    "predictions": ("predictions",),
    "comparison": ("comparison",),
    "questions": ("questions",),
    "planning": ("planning",),
    "report": ("report",),
}

STAGE_CATALOGS = {
    "search": {"design", "editor_definition", "intended_locus", "locus", "reference_bundle",
               "domain", "search_profile", "donor"},
    "observation_normalization": {"design", "editor_definition", "intended_locus",
        "design_alias", "experimental_context", "endpoint_definition", "source_artifact",
        "source_observation", "observation_revision", "observation_state", "observation_equivalence"},
    "evidence_applicability": {"design", "editor_definition", "experimental_context",
        "endpoint_definition", "donor", "specimen", "assay_capability", "observation_state",
        "observation_equivalence", "source_artifact"},
    "predictions": {"design", "editor_definition", "experimental_context", "endpoint_definition"},
    "comparison": {"design", "editor_definition", "experimental_context", "endpoint_definition"},
    "questions": {"question", "requirement_slot", "action", "resource", "donor",
                  "experimental_context", "endpoint_definition", "design", "editor_definition",
                  "intended_locus", "locus", "assay_capability", "specimen", "source_artifact",
                  "sequence_configuration", "sequence_occurrence", "design_configuration",
                  "donor_configuration", "intended_assessment", "observation_revision",
                  "source_observation", "observation_state", "observation_equivalence",
                  "evidence_applicability", "prediction"},
    "planning": {"question", "requirement_slot", "action", "resource", "assay_capability",
                 "specimen", "donor", "experimental_context", "endpoint_definition"},
    "report": set(),
}


def _identity(kind, payload):
    from .comparison_identity import content_id
    return content_id("comparison_" + kind, payload)


def _validate_records(rows, references):
    from .comparison_contracts import validate_catalog, validate_record
    for row in rows:
        validate_record(row)
    merged = {r["record_id"]: r for r in references}
    for row in rows:
        old = merged.get(row["record_id"])
        if old is not None and old != row:
            raise PopulationError("Conflicting comparison records share an identity")
        merged[row["record_id"]] = row
    validate_catalog(list(merged.values()))


def _merge_records(*groups):
    rows = {}
    for group in groups:
        for row in group:
            identifier = row["record_id"]
            if identifier in rows and rows[identifier] != row:
                raise PopulationError("Conflicting comparison records share an identity")
            rows[identifier] = row
    return [rows[key] for key in sorted(rows)]


def resolve_study(path, *, code=None):
    """Resolve declared local files without downloading or invoking external tools."""
    from .comparison_contracts import validate_study
    config_path = Path(path).resolve()
    config = validate_study(read_json(config_path))
    sources = {}
    for ref in config["input_refs"]:
        source = (config_path.parent / ref["path"]).resolve()
        if not source.is_file():
            raise PopulationError("Required comparison input is unavailable: " + ref["role"])
        actual = checksum(source)
        if actual != ref["sha256"]:
            raise PopulationError("Comparison input checksum changed: " + ref["role"])
        sources[ref["input_id"]] = {**ref, "path": str(source), "sha256": actual,
                                    "bytes": source.stat().st_size}
    code = deepcopy(code if code is not None else code_identity())
    payload = {"config": config, "sources": {k: v["sha256"] for k, v in sources.items()},
               "code": code}
    return {"config_path": config_path, "config": config, "sources": sources,
            "code": code, "snapshot_run_id": _identity("snapshot_run", payload)}


def verify_sources(study):
    """Recheck bytes before using a study that was resolved earlier."""
    if read_json(study["config_path"]) != study["config"]:
        raise PopulationError("Comparison configuration changed after resolution")
    from .comparison_domain import require_supported_request
    require_supported_request(study["config"])
    for source in study["sources"].values():
        path = Path(source["path"])
        if not path.is_file() or checksum(path) != source["sha256"]:
            raise PopulationError("Comparison source changed after resolution: " + source["role"])


def effective_inputs(study, stage):
    if stage not in STAGE_PARENTS:
        raise PopulationError("Unknown comparison stage")
    config = study["config"]
    records = [r for r in config["catalogs"]
               if r["record_type"].removeprefix("comparison_") in STAGE_CATALOGS[stage]]
    records.sort(key=lambda r: (r["record_type"], r["record_id"]))
    sources = {key: {k: source[k] for k in ("role", "sha256", "format")}
               for key, source in study["sources"].items() if stage in source["stages"]}
    policies = {k: config["policies"][k] for k in STAGE_POLICIES[stage]}
    if stage == "questions":
        policies["requirement_evaluation"] = {key: config["policies"]["planning"][key]
                                             for key in ("max_expansions", "cost_unit")}
    result = {"mode": config["mode"], "coverage_basis": config["coverage_basis"],
              "catalogs": records, "sources": sources, "policies": policies}
    from .comparison_domain import assess_domain
    domain = assess_domain(config)
    if domain is not None:
        result["supported_domain"] = domain
    return result


def stage_key(study, stage, parent_keys=None, parent_outputs=None):
    parent_keys = parent_keys or {}
    parent_outputs = parent_outputs or {}
    required = set(STAGE_PARENTS.get(stage, ()))
    allowed = required | ({"predictions"} if stage == "comparison" else set())
    if not required <= set(parent_keys) or set(parent_keys) - allowed:
        raise PopulationError("Comparison stage dependencies do not match its contract")
    if set(parent_outputs) != set(parent_keys):
        raise PopulationError("Comparison stage requires the actual upstream artifact hashes")
    if stage == "comparison" and study["config"]["policies"]["comparison"]["include_predictions"] \
            != ("predictions" in parent_keys):
        raise PopulationError("Prediction dependency does not match comparison policy")
    return _identity("stage", {"stage": stage, "effective_inputs": effective_inputs(study, stage),
                               "parents": parent_keys, "parent_outputs": parent_outputs,
                               "code": study["code"]})


def _manifest(path):
    value = read_json(path)
    fields(value, {"record_type", "schema_version", "record_id", "stage", "stage_key",
        "snapshot_run_id", "status", "effective_inputs", "code", "parents", "outputs",
        "created_at"}, label="comparison stage manifest")
    if value["record_type"] != "comparison_stage_manifest" or value["schema_version"] != "1" \
            or value["status"] != "completed" or value["record_id"] != value["stage_key"]:
        raise PopulationError("Incomplete or incompatible comparison stage")
    if not isinstance(value["parents"], dict) or not isinstance(value["outputs"], dict):
        raise PopulationError("Invalid comparison manifest references")
    for key in ("record_id", "stage", "stage_key", "snapshot_run_id", "created_at"):
        if not isinstance(value[key], str) or not value[key]:
            raise PopulationError("Invalid comparison manifest " + key)
    try:
        if datetime.fromisoformat(value["created_at"]).tzinfo is None:
            raise ValueError("Missing timezone")
    except ValueError as exc:
        raise PopulationError("Invalid comparison manifest timestamp") from exc
    return value


def _json_record(value):
    """Apply the same registry to canonical records saved as JSON or JSONL."""
    if isinstance(value, dict) and "record_type" in value:
        from .comparison_contracts import validate_record
        validate_record(value)
        return [value]
    return []


def _verify_native_evidence(artifacts):
    """Native raw evidence remains a checked dependency of normalized stages."""
    audit = artifacts.get("search-audit.json")
    if audit is None:
        return
    from .external_search import verify_native_job
    for execution in audit["execution"]:
        bindings = execution.get("native_manifest")
        if bindings is not None:
            if not isinstance(bindings, list) or not bindings:
                raise PopulationError("Native search requires its retained job bindings")
            for binding in bindings:
                result = verify_native_job(binding)
                if execution["execution_status"] == "complete" and result["execution_status"] != "complete":
                    raise PopulationError("Incomplete native job cannot support a complete normalized search")


def _load_stage(study, directory, expected, visited, cache):
    """Verify each shared ancestor once during this one loading operation."""
    directory = Path(directory).resolve()
    if directory in visited:
        raise PopulationError("Comparison stage dependency cycle")
    key = (directory, expected)
    if key in cache:
        return cache[key]
    visited = visited | {directory}
    manifest = _manifest(directory / "run-manifest.json")
    if manifest["stage"] != expected:
        raise PopulationError("Comparison stage has the wrong role")
    parent_keys, parent_outputs, ancestor_groups = {}, {}, []
    for stage, ref in manifest["parents"].items():
        fields(ref, {"path", "sha256", "stage_key"}, label="upstream comparison manifest")
        path = (directory / ref["path"]).resolve()
        if path.name != "run-manifest.json" or not path.is_file() or checksum(path) != ref["sha256"]:
            raise PopulationError("Comparison ancestor changed or is unavailable")
        artifacts, upstream, ancestor_records = _load_stage(study, path.parent, stage, visited, cache)
        if upstream["stage_key"] != ref["stage_key"]:
            raise PopulationError("Comparison ancestor key mismatch")
        parent_keys[stage] = upstream["stage_key"]
        parent_outputs[stage] = upstream["outputs"]
        ancestor_groups.append(ancestor_records)
    expected_key = stage_key(study, expected, parent_keys, parent_outputs)
    if manifest["stage_key"] != expected_key or manifest["effective_inputs"] != effective_inputs(study, expected) \
            or manifest["code"] != study["code"]:
        raise PopulationError("Comparison stage is stale for the current inputs")
    artifacts, rows = {}, []
    for name, expected_hash in manifest["outputs"].items():
        if Path(name).name != name or name == "run-manifest.json":
            raise PopulationError("Invalid comparison artifact name")
        path = directory / name
        if not path.is_file() or checksum(path) != expected_hash:
            raise PopulationError("Comparison stage artifact changed or is unavailable")
        if name.endswith(".jsonl"):
            artifacts[name] = read_jsonl(path)
            rows.extend(artifacts[name])
        elif name.endswith(".json"):
            artifacts[name] = read_json(path)
            rows.extend(_json_record(artifacts[name]))
        else:
            artifacts[name] = path.read_text(encoding="utf-8")
    ancestors = _merge_records(*ancestor_groups)
    _validate_records(rows, [*study["config"]["catalogs"], *ancestors])
    _verify_native_evidence(artifacts)
    result = artifacts, manifest, _merge_records(ancestors, rows)
    cache[key] = result
    return result


def load_stage(study, directory, expected):
    """Verify sources, ancestors, bytes and record contracts before cache reuse."""
    verify_sources(study)
    artifacts, manifest, _ = _load_stage(study, directory, expected, set(), {})
    verify_sources(study)
    return artifacts, manifest


def load_stage_records(study, directory, expected):
    """Return verified canonical ancestry, current output files and their manifest.

    Catalog records remain in the study. Repeated ancestry is deduplicated by
    content identity; a second call performs fresh byte and source checks.
    """
    verify_sources(study)
    artifacts, manifest, records = _load_stage(study, directory, expected, set(), {})
    verify_sources(study)
    return records, artifacts, manifest


def save_stage(study, stage, directory, artifacts, *, parents=None):
    """Write an immutable, complete stage atomically after validating its records."""
    import hashlib
    import os

    verify_sources(study)
    parents = parents or {}
    references = list(study["config"]["catalogs"])
    parent_refs, parent_keys, parent_outputs, cache = {}, {}, {}, {}
    for name, path in parents.items():
        data, manifest, ancestor_records = _load_stage(study, path, name, set(), cache)
        manifest_path = Path(path).resolve() / "run-manifest.json"
        parent_keys[name] = manifest["stage_key"]
        parent_outputs[name] = manifest["outputs"]
        parent_refs[name] = {"path": os.path.relpath(manifest_path, Path(directory).resolve()),
                             "sha256": checksum(manifest_path), "stage_key": manifest["stage_key"]}
        references.extend(ancestor_records)
    key = stage_key(study, stage, parent_keys, parent_outputs)
    rows = [row for name, records in artifacts.items() if name.endswith(".jsonl") for row in records]
    rows.extend(row for name, value in artifacts.items() if name.endswith(".json")
                for row in _json_record(value))
    _validate_records(rows, references)
    _verify_native_evidence(artifacts)
    files = {}
    for name, value in artifacts.items():
        if Path(name).name != name or name == "run-manifest.json":
            raise PopulationError("Invalid comparison artifact name")
        files[name] = jsonl_text(value) if name.endswith(".jsonl") else \
                      json_text(value) if name.endswith(".json") else value
        if not isinstance(files[name], str):
            raise PopulationError("Comparison text artifact must be a string")
    execution = artifacts.get("search-audit.json", {}).get("execution", []) if stage == "search" else []
    complete = all(row["execution_status"] == "complete" and row["retention_status"] == "complete"
                   and row["exit_code"] == 0 for row in execution)
    manifest = {"record_type": "comparison_stage_manifest", "schema_version": "1",
        "record_id": key, "stage": stage, "stage_key": key,
        "snapshot_run_id": study["snapshot_run_id"], "status": "completed" if complete else "partial",
        "effective_inputs": effective_inputs(study, stage), "code": study["code"],
        "parents": parent_refs,
        "outputs": {name: hashlib.sha256(value.encode("utf-8")).hexdigest() for name, value in files.items()},
        "created_at": datetime.now(timezone.utc).isoformat()}
    files["run-manifest.json"] = json_text(manifest)
    inputs = [study["config_path"], *[s["path"] for s in study["sources"].values()], *parents.values()]
    write_bundle(directory, files, inputs=inputs, overwrite=False)
    return manifest
