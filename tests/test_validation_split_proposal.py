"""Check the public record of the locked tuning and held-out assignment."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PROPOSAL = ROOT / "validation" / "SPLIT_PROPOSAL.md"
REGISTRY = ROOT / "validation" / "dataset_registry.md"


def _read(path: Path) -> str:
    return path.read_text()


def _table_rows(text: str) -> dict[str, list[str]]:
    rows = {}
    for line in text.splitlines():
        if not line.startswith("| D"):
            continue
        cells = [cell.strip() for cell in line.strip().strip("|").split("|")]
        dataset_id = cells[0].split()[0]
        if dataset_id not in {"D1", "D2", "D3", "D4", "D5"}:
            continue
        rows[dataset_id] = cells
    return rows


def test_every_dataset_has_one_locked_role():
    rows = _table_rows(_read(PROPOSAL))
    assert set(rows) == {"D1", "D2", "D3", "D4", "D5"}
    assert {dataset_id: cells[1] for dataset_id, cells in rows.items()} == {
        "D1": "tuning",
        "D2": "tuning",
        "D3": "tuning",
        "D4": "tuning",
        "D5": "held-out",
    }
    assert all(len(cells) == 4 for cells in rows.values())
    assert all(cells[3] == "excluded pending QC" for cells in rows.values())


def test_lock_record_and_registry_roles_agree():
    proposal = _read(PROPOSAL)
    proposal_roles = {dataset_id: cells[1] for dataset_id, cells in _table_rows(proposal).items()}
    registry_roles = {dataset_id: cells[9] for dataset_id, cells in _table_rows(_read(REGISTRY)).items()}

    assert proposal_roles == registry_roles
    assert "**Status:** locked on 2026-08-14 by Tram Ngo, protocol owner." in proposal
    assert "Locked by: Tram Ngo | Date: 2026-08-14" in proposal


def test_held_out_scope_is_file_level_and_context_matched():
    d5_reason = _table_rows(_read(PROPOSAL))["D5"][2]
    assert "file-level subset" in d5_reason
    assert "primary CD34+ HSPCs" in d5_reason
    assert "HUDEP-2 and nuclease-comparator records are excluded" in d5_reason


def test_safeguards_and_unverified_items_are_explicit():
    text = _read(PROPOSAL)
    for required in (
        "Every measurement and table join must include the edited allele",
        "D5 may be run once only",
        "patient risk, incidence, or clinical-outcome estimates",
        "requires a new protocol version",
        "No D1-D5 row has passed all registry admission checks",
        "No validation metric is claimed",
    ):
        assert required in text
