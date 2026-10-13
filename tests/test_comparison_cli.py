"""The comparison entry point must validate local sources without disclosing them."""
import json
from pathlib import Path
import shutil

from base_edit_safety.pipeline.comparison import main


ROOT = Path(__file__).resolve().parents[1]


def test_validate_outputs_only_status_identity_and_counts(capsys):
    assert main(["validate", "--config", str(ROOT / "examples/comparison/contracts/study.json")]) == 0
    captured = capsys.readouterr()
    result = json.loads(captured.out)
    assert set(result) == {"status", "command", "snapshot_run_id", "catalog_records", "local_sources"}
    assert result["status"] == "completed"
    assert str(ROOT) not in captured.out + captured.err


def test_missing_protected_path_is_not_printed(tmp_path, capsys):
    marker = "private-specimen-name-never-print"
    assert main(["validate", "--config", str(tmp_path / marker / "study.json")]) == 2
    captured = capsys.readouterr()
    assert json.loads(captured.err)["status"] == "invalid_input"
    assert marker not in captured.out + captured.err
    assert str(tmp_path) not in captured.out + captured.err


def test_changed_input_is_rejected_without_disclosing_file_path(tmp_path, capsys):
    folder = tmp_path / "protected-donor-directory"
    shutil.copytree(ROOT / "examples/comparison/contracts", folder)
    config = json.loads((folder / "study.json").read_text())
    path = folder / config["input_refs"][0]["path"]
    path.write_text(path.read_text() + "\n")
    assert main(["validate", "--config", str(folder / "study.json")]) == 2
    captured = capsys.readouterr()
    assert "checksum" in json.loads(captured.err)["reason"]
    assert str(folder) not in captured.err


def test_observation_cli_runs_the_actual_importer_and_verifies_reuse(tmp_path, capsys):
    config = ROOT / "examples/comparison/observations/study.json"
    output = tmp_path / "private-observations"
    args = ["observations", "--config", str(config), "--output", str(output)]
    assert main(args) == 0
    result = json.loads(capsys.readouterr().out)
    assert result["artifacts"] == 2
    rows = [json.loads(line) for line in (output / "observations.jsonl").read_text().splitlines()]
    revision = next(r for r in rows if r["record_type"] == "comparison_observation_revision")
    assert revision["normalized_value"]["value"] == "0.098000000000000007"
    assert main(args + ["--reuse"]) == 0
    captured = capsys.readouterr()
    assert str(output) not in captured.out + captured.err
    assert "synthetic-donor" not in captured.out + captured.err
