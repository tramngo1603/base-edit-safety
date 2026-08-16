"""Pre-run checks for the registry, protocol lock, CI, and held-out evaluation."""

from __future__ import annotations

from dataclasses import dataclass
import os
from pathlib import Path
import re
import subprocess
from typing import Iterable

from .thresholds import LOCKED_THRESHOLDS


class SafetyViolation(RuntimeError):
    """Raised before a disallowed evaluation can read prediction-measurement comparisons."""


@dataclass(frozen=True)
class RegistryEntry:
    dataset_id: str
    assay: str
    cell_type: str
    role: str
    qc_status: str
    raw_cells: tuple[str, ...]


def parse_registry(path: str | Path) -> dict[str, RegistryEntry]:
    entries: dict[str, RegistryEntry] = {}
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if not re.match(r"^\|\s*D\d+\b", line):
            continue
        cells = tuple(cell.strip() for cell in line.strip().strip("|").split("|"))
        if len(cells) < 11:
            raise SafetyViolation(f"malformed dataset-registry row: {line}")
        dataset_id = cells[0].split()[0]
        if dataset_id in entries:
            raise SafetyViolation(f"duplicate dataset-registry row for {dataset_id}")
        role = cells[9].strip().lower()
        if role not in {"tuning", "held-out"}:
            raise SafetyViolation(f"registry role for {dataset_id} is not locked: {cells[9]!r}")
        entries[dataset_id] = RegistryEntry(
            dataset_id=dataset_id,
            assay=cells[3],
            cell_type=cells[4],
            role=role,
            qc_status=cells[10].strip().lower(),
            raw_cells=cells,
        )
    if not entries:
        raise SafetyViolation(f"no dataset rows found in registry: {path}")
    return entries


def assert_dataset_roles(
    dataset_ids: Iterable[str],
    registry: dict[str, RegistryEntry],
    *,
    expected_role: str,
) -> None:
    selected = sorted(set(dataset_ids))
    if not selected:
        raise SafetyViolation("evaluation selected no registered dataset IDs")
    unknown = [dataset_id for dataset_id in selected if dataset_id not in registry]
    if unknown:
        raise SafetyViolation(f"datasets absent from the locked registry: {unknown}")
    wrong = {
        dataset_id: registry[dataset_id].role
        for dataset_id in selected
        if registry[dataset_id].role != expected_role
    }
    if wrong:
        raise SafetyViolation(
            f"{expected_role} entry point refuses dataset role(s): {wrong}; "
            "the locked registry split cannot be overridden"
        )


def assert_protocol_frozen(protocol_path: str | Path) -> None:
    """Reject research-data evaluation while required protocol fields remain unfilled."""
    text = Path(protocol_path).read_text(encoding="utf-8")
    header = text.split("## 1.", 1)[0]
    complaints = []
    if re.search(r"Status:\s*.*\bDRAFT\b", header, re.IGNORECASE):
        complaints.append("status is DRAFT")
    frozen = re.search(r"^- Date frozen:\s*(.+)$", header, re.MULTILINE)
    if frozen is None or "{{" in frozen.group(1) or not re.search(r"\d{4}-\d{2}-\d{2}", frozen.group(1)):
        complaints.append("Date frozen is not concrete")
    release = re.search(r"^- Code release under test:\s*(.+)$", header, re.MULTILINE)
    if release is None or "{{" in release.group(1):
        complaints.append("code release/tag is not concrete")
    section4 = text.split("## 4.", 1)[1].split("\n## 5.", 1)[0] if "## 4." in text else ""
    if not section4 or "{{" in section4:
        complaints.append("§4 still contains unlocked template placeholders")
    t = LOCKED_THRESHOLDS
    expected_fragments = {
        "M1": (f"{t.m1_guide_seq_class_recall:.2f}", f"{t.m1_biochemical_recall:.2f}"),
        "M2": (str(t.m2_top_k), f"{t.m2_rank_quality:.2f}"),
        "M3": (f"{t.m3_spearman_rho:.1f}",),
        "M4": (f"{t.m4_spearman_rho:.1f}",),
        "M5": (f"{t.m5_sensitivity:.0%}",),
        "M6": (f"{t.m6_sensitivity:.0%}", f"{t.m6_specificity:.2f}"),
        "M7": (f"{t.m7_sensitivity:.0%}",),
        "M8": (str(t.m8_audit_sample_size), f"{t.m8_accuracy:.0%}"),
        "M9": (f"{t.m9_coverage:.0%}",),
    }
    section_rows = {
        cells[0]: line
        for line in section4.splitlines()
        if line.startswith("| M")
        and len(cells := [cell.strip() for cell in line.strip().strip("|").split("|")]) >= 5
    }
    for metric_id, fragments in expected_fragments.items():
        row = section_rows.get(metric_id, "")
        if not row or any(fragment not in row for fragment in fragments):
            complaints.append(f"§4 {metric_id} does not match the locked code threshold")
    if complaints:
        raise SafetyViolation(
            "non-fixture evaluation is forbidden before protocol lock: " + "; ".join(complaints)
        )


def assert_registry_qc_admitted(registry: dict[str, RegistryEntry]) -> None:
    blocked = {
        dataset_id: entry.qc_status
        for dataset_id, entry in registry.items()
        if entry.qc_status not in {"pass", "warn"}
    }
    if blocked:
        raise SafetyViolation(
            "non-fixture evaluation is forbidden until every registry QC status is pass/warn: "
            f"{blocked}"
        )


def _truthy_environment(name: str) -> bool:
    return os.environ.get(name, "").strip().lower() in {"1", "true", "yes", "on"}


def assert_ci_fixture_mode(*, fixture_mode: bool) -> None:
    if (_truthy_environment("CI") or _truthy_environment("BASE_EDIT_SAFETY_FIXTURES_ONLY")) and not fixture_mode:
        raise SafetyViolation("CI may execute the evaluation runner only in fixture mode")


def assert_fixture_inputs(paths: Iterable[str | Path], fixture_root: str | Path) -> None:
    root = Path(fixture_root).resolve(strict=True)
    for raw_path in paths:
        path = Path(raw_path).resolve(strict=True)
        try:
            path.relative_to(root)
        except ValueError as exc:
            raise SafetyViolation(
                f"fixture mode refuses non-fixture primary input {path}; expected it below {root}"
            ) from exc


_TAG_RE = re.compile(r"[A-Za-z0-9][A-Za-z0-9._/+\-]*")


def verify_freeze_tag(repo_root: str | Path, tag: str) -> str:
    if not _TAG_RE.fullmatch(tag) or ".." in tag or tag.endswith("."):
        raise SafetyViolation(f"invalid freeze-tag name: {tag!r}")
    completed = subprocess.run(
        ["git", "-C", str(repo_root), "rev-parse", "--verify", f"refs/tags/{tag}^{{commit}}"],
        capture_output=True,
        text=True,
        timeout=10,
        check=False,
    )
    commit = completed.stdout.strip()
    if completed.returncode != 0 or not re.fullmatch(r"[0-9a-fA-F]{40,64}", commit):
        raise SafetyViolation(f"held-out run refuses missing/non-commit freeze tag {tag!r}")
    return commit.lower()
