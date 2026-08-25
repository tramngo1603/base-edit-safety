"""Runtime configuration for external prediction tools.

You install the external predictors yourself, and thin adapters call them. You give their
locations through environment variables, so no machine-specific path is committed.
See requirements-external.txt for the tools and their licenses.
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path


_REPO_ROOT = Path(__file__).resolve().parents[2]
_DOTENV_PATH = _REPO_ROOT / ".env"


def _setting(name: str) -> str | None:
    """Read an environment setting, with the ignored project ``.env`` as a local fallback."""
    value = os.environ.get(name)
    if value:
        return value
    try:
        lines = _DOTENV_PATH.read_text().splitlines()
    except OSError:
        return None
    for raw in lines:
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("export "):
            line = line[len("export "):].lstrip()
        key, separator, candidate = line.partition("=")
        if separator and key.strip() == name:
            candidate = candidate.strip()
            if len(candidate) >= 2 and candidate[0] == candidate[-1] and candidate[0] in "\"'":
                candidate = candidate[1:-1]
            return candidate or None
    return None


class ExternalToolNotConfigured(RuntimeError):
    """Raised when code asks for an external predictor but its location is not set."""


@dataclass(frozen=True)
class BedictConfig:
    """Where the installed base-editing outcome predictor lives.

    repo_root        directory of the installed predictor checkout (holds trained_models/).
    python_executable python interpreter of the separate environment the predictor was installed into.
    """

    repo_root: str
    python_executable: str


def load_bedict_config() -> BedictConfig:
    """Read the external predictor location from environment variables.

    Environment variables:
      BE_SAFETY_BEDICT_REPO    path to the installed predictor checkout
      BE_SAFETY_BEDICT_PYTHON  path to that environment's python interpreter
    """
    repo = _setting("BE_SAFETY_BEDICT_REPO")
    py = _setting("BE_SAFETY_BEDICT_PYTHON")
    if not repo or not py:
        raise ExternalToolNotConfigured(
            "Base-editing outcome predictor not configured. Set BE_SAFETY_BEDICT_REPO "
            "(installed checkout containing trained_models/) and BE_SAFETY_BEDICT_PYTHON "
            "(its environment's interpreter)."
        )
    if not os.path.isdir(repo):
        raise ExternalToolNotConfigured(f"BE_SAFETY_BEDICT_REPO is not a directory: {repo}")
    if not os.path.isfile(py):
        raise ExternalToolNotConfigured(f"BE_SAFETY_BEDICT_PYTHON is not a file: {py}")
    return BedictConfig(repo_root=repo, python_executable=py)


def bedict_configured() -> bool:
    """True when the external predictor environment can be found (used to gate integration tests)."""
    try:
        load_bedict_config()
        return True
    except ExternalToolNotConfigured:
        return False


@dataclass(frozen=True)
class CrisporConfig:
    """Where the installed guide-design tool and its genome directory live.

    repo_root         directory of the installed guide-design checkout (holds crispor.py).
    python_executable python interpreter of the environment it was installed into.
    genome            genome identifier present under genome_dir (e.g. hg38).
    genome_dir        directory holding the indexed genome(s).
    """

    repo_root: str
    python_executable: str
    genome: str
    genome_dir: str


def load_crispor_config() -> CrisporConfig:
    """Read the external guide-design tool location from environment variables.

    Environment variables:
      BE_SAFETY_CRISPOR_REPO    path to the installed guide-design checkout (contains crispor.py)
      BE_SAFETY_CRISPOR_PYTHON  path to that environment's python interpreter
      BE_SAFETY_CRISPOR_GENOME  genome identifier (default: hg38)
      BE_SAFETY_CRISPOR_GENOMEDIR  genome directory (default: <repo>/genomes)
    """
    repo = _setting("BE_SAFETY_CRISPOR_REPO")
    py = _setting("BE_SAFETY_CRISPOR_PYTHON")
    if not repo or not py:
        raise ExternalToolNotConfigured(
            "Guide-design tool not configured. Set BE_SAFETY_CRISPOR_REPO (checkout containing "
            "crispor.py) and BE_SAFETY_CRISPOR_PYTHON (its environment's interpreter)."
        )
    if not os.path.isfile(os.path.join(repo, "crispor.py")):
        raise ExternalToolNotConfigured(f"crispor.py not found under BE_SAFETY_CRISPOR_REPO: {repo}")
    if not os.path.isfile(py):
        raise ExternalToolNotConfigured(f"BE_SAFETY_CRISPOR_PYTHON is not a file: {py}")
    genome = _setting("BE_SAFETY_CRISPOR_GENOME") or "hg38"
    genome_dir = _setting("BE_SAFETY_CRISPOR_GENOMEDIR") or os.path.join(repo, "genomes")
    return CrisporConfig(repo_root=repo, python_executable=py, genome=genome, genome_dir=genome_dir)


def crispor_configured() -> bool:
    """True when the external guide-design tool can be found (used to gate integration tests)."""
    try:
        load_crispor_config()
        return True
    except ExternalToolNotConfigured:
        return False
