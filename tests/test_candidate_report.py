"""Static reports compose design-record v7 and candidate-ledger v3 without embedding either."""
from copy import deepcopy
import json
from pathlib import Path

import pytest

from base_edit_safety.pipeline.candidate_panel import build_ledger
from base_edit_safety.pipeline.candidate_report import (
    CandidateReportError,
    main,
    render_report,
)


ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples" / "candidate-panel"


def ledger() -> dict:
    return build_ledger(
        EXAMPLE / "search-profile.json",
        EXAMPLE / "guides.tsv",
        EXAMPLE / "intended-loci.tsv",
        [EXAMPLE / "nominations.tsv"],
        EXAMPLE / "measurements.tsv",
        EXAMPLE / "panel-decisions.tsv",
        [EXAMPLE / "source-provenance.json"],
        panel_notes="Synthetic depth context.",
    )


def design_records() -> dict:
    return {
        "record_type": "hbg_design_record_set",
        "schema_version": "7",
        "count": 1,
        "records": [
            {
                "record_type": "hbg_design_record",
                "schema_version": "7",
                "candidate": {
                    "guide_id": "synthetic_HBG_guide",
                    "gene": "HBG1",
                    "editor_name": "ABE8e",
                    "target": {
                        "intended_edits": [{"promoter_offset": -113}],
                    },
                },
                "evidence_tier": "predicted_only",
                "hbg_paralog_flags": {"flagged": True},
            }
        ],
        "limitations": [],
        "reproducibility": {},
    }


def test_report_has_fixed_sections_and_leaves_inputs_unchanged():
    designs = design_records()
    candidates = ledger()
    before_designs = deepcopy(designs)
    before_candidates = deepcopy(candidates)

    report = render_report(designs, candidates)

    assert "## Design summary" in report
    assert "## Nomination summary" in report
    assert "## Candidate evidence table" in report
    assert "## Measurements and follow-up" in report
    assert "## Panel decision summary" in report
    assert "## Reproducibility and run diff" in report
    assert "Synthetic depth context." in report
    assert "zero-based and half-open" in report
    assert designs == before_designs
    assert candidates == before_candidates


def test_report_refuses_design_set_missing_a_ledger_guide():
    designs = design_records()
    designs["records"][0]["candidate"]["guide_id"] = "another-guide"

    with pytest.raises(CandidateReportError, match="missing ledger guides"):
        render_report(designs, ledger())


def test_report_cli_records_both_input_checksums(tmp_path, capsys):
    design_path = tmp_path / "design-records.json"
    ledger_path = tmp_path / "candidate-ledger.json"
    output_path = tmp_path / "candidate-report.md"
    design_path.write_text(json.dumps(design_records()), encoding="utf-8")
    ledger_path.write_text(json.dumps(ledger()), encoding="utf-8")

    status = main(
        [
            "--design-records",
            str(design_path),
            "--candidate-ledger",
            str(ledger_path),
            "--output",
            str(output_path),
        ]
    )

    report = output_path.read_text()
    assert status == 0
    assert "`design-records.json` (`sha256:" in report
    assert "`candidate-ledger.json` (`sha256:" in report
    assert "wrote candidate review report" in capsys.readouterr().out


def test_report_run_manifest_records_package_and_python_versions():
    candidates = ledger()

    report = render_report(design_records(), candidates)

    runtime = candidates["reproducibility"]
    assert (
        f"- Package: base-edit-safety {runtime['base_edit_safety_version']}" in report
    )
    assert f"- Python: {runtime['python_version']}" in report
