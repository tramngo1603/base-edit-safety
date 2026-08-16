"""Checksums, one-shot locking, and run-manifest construction."""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
from typing import Any, Iterable, Mapping

from .safety import SafetyViolation
from .thresholds import LOCKED_THRESHOLDS


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def checksum_inputs(named_paths: Iterable[tuple[str, str | Path]]) -> list[dict[str, Any]]:
    checksums = []
    for role, raw_path in named_paths:
        path = Path(raw_path)
        if not path.is_file():
            raise FileNotFoundError(path)
        checksums.append(
            {
                "role": role,
                "path": str(path),
                "sha256": sha256_file(path),
                "size_bytes": path.stat().st_size,
            }
        )
    checksums.sort(key=lambda item: (item["role"], item["path"]))
    return checksums


def combined_input_hash(checksums: Iterable[Mapping[str, Any]]) -> str:
    canonical = json.dumps(list(checksums), sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(canonical).hexdigest()


def git_state(repo_root: str | Path) -> dict[str, Any]:
    def run(*arguments: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            ["git", "-C", str(repo_root), *arguments], capture_output=True, text=True,
            timeout=10, check=False,
        )

    revision = run("rev-parse", "HEAD")
    status = run("status", "--porcelain")
    commit = revision.stdout.strip() if revision.returncode == 0 else None
    return {
        "commit_sha": commit,
        "commit_lookup": "ok" if commit else "unavailable",
        "working_tree_dirty": bool(status.stdout.strip()) if status.returncode == 0 else None,
    }


def runtime_versions(extra: Mapping[str, Any] | None = None) -> dict[str, str]:
    try:
        package_version = importlib.metadata.version("base-edit-safety")
    except importlib.metadata.PackageNotFoundError:
        package_version = "0.0.0+source-tree"
    versions = {
        "base-edit-safety": package_version,
        "python": platform.python_version(),
        "python_implementation": platform.python_implementation(),
        "platform": platform.platform(),
    }
    if extra:
        versions.update({str(name): str(value) for name, value in extra.items()})
    return dict(sorted(versions.items()))


def build_manifest(
    *,
    mode: str,
    started_at_utc: str,
    ended_at_utc: str,
    repo_root: str | Path,
    input_checksums: list[dict[str, Any]],
    tool_versions: Mapping[str, Any],
    fixture_mode: bool,
    freeze_tag: str | None = None,
    freeze_tag_commit: str | None = None,
) -> dict[str, Any]:
    return {
        "schema": "base_edit_safety/evaluation-run-manifest/v1",
        "mode": mode,
        "fixture_mode": fixture_mode,
        "started_at_utc": started_at_utc,
        "ended_at_utc": ended_at_utc,
        "repository": git_state(repo_root),
        "freeze_tag": freeze_tag,
        "freeze_tag_commit": freeze_tag_commit,
        "tool_versions": runtime_versions(tool_versions),
        "input_checksums": input_checksums,
        "combined_input_sha256": combined_input_hash(input_checksums),
        "locked_thresholds": LOCKED_THRESHOLDS.as_dict(),
        "argv": list(sys.argv),
    }


def atomic_write_json(path: str | Path, value: Mapping[str, Any]) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(f".{destination.name}.{os.getpid()}.tmp")
    try:
        temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)


def create_one_shot_lock(
    path: str | Path,
    *,
    timestamp_utc: str,
    freeze_tag: str,
    freeze_tag_commit: str,
    input_checksums: list[dict[str, Any]],
) -> dict[str, Any]:
    """Record use of the single held-out run before measurement values are loaded."""
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema": "base_edit_safety/held-out-one-shot-lock/v1",
        "created_at_utc": timestamp_utc,
        "freeze_tag": freeze_tag,
        "freeze_tag_commit": freeze_tag_commit,
        "input_checksums": input_checksums,
        "combined_input_sha256": combined_input_hash(input_checksums),
        "state": "consumed",
        "statement": "The protocol-v1 held-out evaluation has been opened; a second run is forbidden.",
    }
    encoded = (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode("utf-8")
    try:
        descriptor = os.open(destination, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    except FileExistsError as exc:
        raise SafetyViolation(
            f"held-out evaluation is one-shot and was already consumed; lock exists: {destination}"
        ) from exc
    try:
        view = memoryview(encoded)
        while view:
            written = os.write(descriptor, view)
            view = view[written:]
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    return payload
