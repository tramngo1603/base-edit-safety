"""Resolve local sources and reject stage reuse after any input change."""
from __future__ import annotations

from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path
import platform
import subprocess
import sys

from .common import (PopulationError, ResourceLimit, checksum, digest, identity, json_text, jsonl_text,
                     read_json, read_jsonl, record, write_bundle)
from .contracts import validate_config, validate_record
from .variants import Reference


def code_identity():
    root = Path(__file__).resolve().parents[1]
    content = {str(p.relative_to(root)): checksum(p) for p in sorted(root.rglob("*"))
               if p.suffix in {".py", ".json"}}
    return {"package_version": version("base-edit-safety"), "source_sha256": digest(content),
            "python": platform.python_version()}


def runtime():
    root = Path(__file__).resolve().parents[3]
    state = {"revision": None, "dirty": None}
    if (root / ".git").exists():
        try:
            revision = subprocess.run(["git", "rev-parse", "HEAD"], cwd=root, capture_output=True, text=True, timeout=5, check=True)
            status = subprocess.run(["git", "status", "--porcelain"], cwd=root, capture_output=True, text=True, timeout=5, check=True)
            state = {"revision": revision.stdout.strip(), "dirty": bool(status.stdout.strip())}
        except (OSError, subprocess.SubprocessError):
            pass
    return {**code_identity(), **state, "argv": sys.argv,
            "timestamp_utc": datetime.now(timezone.utc).isoformat()}


def resolve_study(path):
    path = Path(path).resolve()
    config = validate_config(read_json(path))
    if sum(i["end"]-i["start"] for i in config["scope"]["materialized_intervals"]) > config["resource_limits"]["max_reference_bases"]:
        raise ResourceLimit("max_reference_bases exceeded before reference access")
    paths = {k: (path.parent / v).resolve() if v is not None else None for k, v in config["inputs"].items()}
    sources = {}
    for role, p in paths.items():
        if p is None:
            continue
        if not p.is_file():
            raise PopulationError("Required local input is unavailable: " + role)
        actual = checksum(p)
        expected = config["input_sources"][role]["sha256"]
        if expected is not None and actual != expected:
            raise PopulationError("Source checksum mismatch: " + role)
        sources[role] = {"path": str(p), "sha256": actual, "bytes": p.stat().st_size,
                         "source": config["input_sources"][role]}
    reference = Reference(paths["reference"], config["reference_policy"], config["coverage_basis"] == "synthetic_truth")
    if reference.index_sha256:
        sources["reference_index"] = {"path": str(paths["reference"])+".fai", "sha256": reference.index_sha256,
                                       "bytes": Path(str(paths["reference"])+".fai").stat().st_size}
    for interval in config["scope"]["materialized_intervals"]:
        reference.fetch(interval["contig"], interval["start"], interval["end"])
    semantic = {"config": config, "sources": {k: v["sha256"] for k, v in sources.items()},
                "code": code_identity()}
    return {"config_path": path, "config": config, "paths": paths, "sources": sources,
            "reference": reference, "config_hash": digest(config), "run_id": identity("run", semantic)}


def artifact_file(path, relative_to=None, role=None):
    path = Path(path).resolve()
    import os
    result = {"path": os.path.relpath(path, relative_to) if relative_to else str(path), "sha256": checksum(path)}
    if role:
        result["role"] = role
    return result


def resolve_file(ref, base):
    if not isinstance(ref, dict) or set(ref) not in ({"path", "sha256"}, {"path", "sha256", "role"}):
        raise PopulationError("Invalid artifact reference")
    path = (Path(base) / ref["path"]).resolve()
    if not path.is_file() or checksum(path) != ref["sha256"]:
        raise PopulationError("Artifact is missing or its checksum changed")
    return path


def save_stage(study, stage, directory, artifacts, *, parents=(), input_files=(), overwrite=False, status="completed", timing=None):
    files = {}
    for name, value in artifacts.items():
        if name.endswith(".jsonl"):
            for row in value:
                validate_record(row)
        files[name] = jsonl_text(value) if name.endswith(".jsonl") else json_text(value)
    import hashlib
    hashes = {name: hashlib.sha256(text.encode("utf-8")).hexdigest() for name, text in files.items()}
    manifest = record("population_run_manifest", study["run_id"], stage=stage, execution_status=status,
        config_hash=study["config_hash"], config=study["config"], sources=study["sources"], outputs=hashes,
        upstream_manifests=[artifact_file(Path(p)/"run-manifest.json") for p in parents],
        runtime=runtime(), elapsed_seconds=timing,
        reference_bounds={c: study["reference"].bounds(c) for c in study["config"]["reference_policy"]["permitted_contigs"]})
    files["run-manifest.json"] = json_text(manifest)
    inputs = [study["config_path"]] + [r["path"] for r in study["sources"].values()] + list(parents) + list(input_files)
    write_bundle(directory, files, inputs=inputs, overwrite=overwrite)
    return manifest


def load_stage(study, directory, expected):
    directory = Path(directory).resolve()
    manifest = read_json(directory/"run-manifest.json")
    if manifest.get("record_type") != "population_run_manifest" or manifest.get("schema_version") != "1" \
            or manifest.get("run_id") != study["run_id"] or manifest.get("stage") != expected \
            or manifest.get("execution_status") != "completed" or manifest.get("config_hash") != study["config_hash"]:
        raise PopulationError("Stale, failed or incompatible upstream stage")
    for ref in manifest["upstream_manifests"]:
        resolve_file(ref, directory)
    artifacts = {}
    for name, sha in manifest["outputs"].items():
        if Path(name).name != name:
            raise PopulationError("Invalid stage artifact name")
        path = resolve_file({"path": name, "sha256": sha}, directory)
        artifacts[name] = read_jsonl(path) if name.endswith(".jsonl") else read_json(path)
        if name.endswith(".jsonl"):
            for row in artifacts[name]:
                validate_record(row)
                if row["run_id"] != study["run_id"]:
                    raise PopulationError("Population row belongs to a different run")
    return artifacts, manifest
