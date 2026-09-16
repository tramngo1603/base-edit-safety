"""Records what a run depended on, so an output file can be traced back to what produced it.

That means tool paths and git commits where they can be found, the genome build and index, the
editors, how many predictor runs were averaged, and checksums of the reference windows. None of it
is essential to the run, so nothing here raises: if a lookup fails, that shows up in the output.
"""
from __future__ import annotations

import datetime
import os
import subprocess

from .. import anchors

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))


def git_commit(path: str | None) -> str | None:
    """Short git commit of a checkout, or None if it cannot be read."""
    if not path:
        return None
    try:
        done = subprocess.run(["git", "-C", path, "rev-parse", "--short", "HEAD"],
                              capture_output=True, text=True, timeout=5)
    except Exception:
        return None
    return (done.stdout.strip() or None) if done.returncode == 0 else None


def tool_environment(name: str, repo_root: str | None = None,
                     python_executable: str | None = None, genome: str | None = None,
                     **extra) -> dict:
    """What one external tool was: its paths, its genome index, and its git commit if there is one.

    A failed lookup never breaks the run. It is written into the output as unavailable rather than
    quietly left out.
    """
    env: dict = {"name": name}
    if repo_root is not None:
        env["repo_root"] = repo_root
        commit = git_commit(repo_root)
        env["commit"] = commit
        env["commit_lookup"] = "ok" if commit else "unavailable"
    if python_executable is not None:
        env["python_executable"] = python_executable
    if genome is not None:
        env["genome_index"] = genome
    env.update(extra)
    return env


def reproducibility(tools: list[dict], editors: dict, ensemble_runs: int | None,
                    reference_windows: dict) -> dict:
    """Collect the run settings and source versions for the output metadata."""
    return {
        "genome_build": anchors.GENOME_BUILD,
        "tools": tools,
        "editors": editors,
        "ensemble_runs": ensemble_runs,
        "reference_windows": reference_windows,
    }


def repo_state() -> dict:
    """This repository's git commit, and whether it had uncommitted changes at the time."""
    commit = git_commit(_REPO_ROOT)
    dirty = None
    try:
        done = subprocess.run(["git", "-C", _REPO_ROOT, "status", "--porcelain"],
                              capture_output=True, text=True, timeout=5)
        dirty = bool(done.stdout.strip()) if done.returncode == 0 else None
    except Exception:
        dirty = None
    return {"commit": commit, "dirty": dirty, "commit_lookup": "ok" if commit else "unavailable"}


def run_stamp() -> str:
    """UTC timestamp for an output file (ISO 8601)."""
    return datetime.datetime.now(datetime.timezone.utc).isoformat()
