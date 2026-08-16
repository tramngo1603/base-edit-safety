"""Require the tuning/held-out assignment and lock date to match the protocol version.

The registry `role` column and protocol-section-3 lock date are tied to the protocol version.
Changing a role or lock date requires a new protocol version and a matching entry below.
"""
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "validation" / "dataset_registry.md"
PROTOCOL = ROOT / "validation" / "validation-protocol.md"

# protocol version -> (locked role map, §3 lock date).
# Changing the split requires a new protocol version and a new entry; existing entries remain fixed.
LOCKED_SPLIT_BY_PROTOCOL_VERSION = {
    "v1": (
        {"D1": "tuning", "D2": "tuning", "D3": "tuning", "D4": "tuning", "D5": "held-out"},
        "2026-08-14",
    ),
}


def _registry_roles() -> dict[str, str]:
    roles = {}
    for line in REGISTRY.read_text().splitlines():
        if not line.startswith("| D"):
            continue
        cells = [cell.strip() for cell in line.strip().strip("|").split("|")]
        roles[cells[0].split()[0]] = cells[9]
    return roles


def _protocol_version(text: str) -> str:
    for line in text.splitlines():
        if line.strip().startswith("- Protocol version:"):
            return line.split(":", 1)[1].strip().split()[0]
    raise AssertionError("protocol version field is missing from validation-protocol.md")


def _section_3(text: str) -> str:
    return text.split("## 3. ", 1)[1].split("\n## 4.", 1)[0]


def test_registry_roles_are_fully_assigned_no_tbd():
    roles = _registry_roles()
    assert set(roles) == {"D1", "D2", "D3", "D4", "D5"}
    assert "TBD" not in roles.values(), roles


def test_locked_split_is_pinned_to_the_protocol_version():
    text = PROTOCOL.read_text()
    version = _protocol_version(text)
    assert version in LOCKED_SPLIT_BY_PROTOCOL_VERSION, (
        f"protocol version {version!r} has no frozen-split entry; add one when bumping the version"
    )
    expected_roles, expected_lock_date = LOCKED_SPLIT_BY_PROTOCOL_VERSION[version]

    # Registry role column must match the entry for the current protocol version exactly.
    assert _registry_roles() == expected_roles

    # §3 must carry the same lock date, keyed to the same version. A role or date change without a
    # version bump leaves this mismatched and fails.
    section3 = _section_3(text)
    assert f"Split locked: {expected_lock_date}" in section3
    assert f"protocol version {version}" in section3
    assert "Tuning set: D1, D2, D3, D4" in section3
    assert "Held-out set: D5" in section3
