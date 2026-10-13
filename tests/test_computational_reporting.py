"""Check that reports match computed results for a hand-specified mock example."""
from copy import deepcopy
import csv
import io
import json
from pathlib import Path
import shutil

from base_edit_safety.population.computational import compute, run
from base_edit_safety.population.computational_reporting import COLUMNS, render
from base_edit_safety.population.common import checksum
from base_edit_safety.pipeline.comparison import main


EXAMPLE = Path(__file__).resolve().parents[1] / "examples/comparison/computational"


def fixture():
    return compute(EXAMPLE / "request.json")[0]


def test_hand_specified_example():
    expected = json.loads((EXAMPLE / "expected.json").read_text())
    result = fixture()
    for key, value in expected["counts"].items():
        assert result["summary"][key] == value
    case = expected["variant_created_candidate"]
    row = next(r for r in result["candidates"] if r["start0"] == case["start0"] and r["sample_id"] == case["sample_id"])
    for key in ("start0", "end0", "strand", "pam", "change", "local_copies"):
        assert row[key] == case[key]
    assert row["target_sequence"] == expected["guide"]
    assert row["reference_target_sequence"] == expected["guide"] + case["reference_pam"]
    for sample in expected["unknown_at_created_site"]:
        row = next(r for r in result["candidates"] if r["start0"] == case["start0"] and r["sample_id"] == sample)
        assert row["model_carriage"] == "unknown"
        assert row["change"] == "unpaired_unevaluable"
    assert all(row["evidence_supported_carriage"] == expected["evidence_supported_carriage"] for row in result["candidates"])
    assert result["domain"]["biological_status"] == expected["biological_status"]


def test_tsv_is_complete_and_ordered():
    result = fixture()
    rendered = render(result)
    assert set(rendered) == {"report.md", "candidates.tsv", "warnings.tsv"}
    reader = csv.DictReader(io.StringIO(rendered["candidates.tsv"]), delimiter="\t")
    assert reader.fieldnames == list(COLUMNS)
    rows = list(reader)
    assert len(rows) == len(result["candidates"])
    for source, row in zip(result["candidates"], rows):
        assert row["candidate_id"] == source["candidate_id"] and row["sample_id"] == source["sample_id"]
        for key in ("local_copies", "phase_status", "target_path", "variant_ids"):
            assert json.loads(row[key]) == source[key]
        assert row["reference_target_sequence"] == (source["reference_target_sequence"] or "unknown")
    warnings = list(csv.DictReader(io.StringIO(rendered["warnings.tsv"]), delimiter="\t"))
    assert len(warnings) == len(result["warnings"])
    for source, row in zip(result["warnings"], warnings):
        assert row["code"] == source["code"] and row["detail"] == source["detail"] and int(row["count"]) == source["count"]


def test_readable_report_keeps_abstention_and_all_warnings():
    result = fixture()
    text = render(result)["report.md"]
    assert "**abstain**" in text and "**artificial**" in text
    assert "genotypes at those positions remain unknown" in text
    assert "Candidates are sequence matches within the declared region." in text
    assert "unordered pair" in text and "run-manifest.json" in text
    assert "genotype" in text.lower() and "paired" in text
    for row in result["warnings"]:
        assert row["code"].replace("_", "\\_") in text


def test_renderer_does_not_mutate_or_replace_unknowns():
    result = fixture()
    before = deepcopy(result)
    first = render(result)
    assert result == before and first == render(result)
    result["summary"]["reference_candidates"] = None
    assert "reference\\_candidates | unknown" in render(result)["report.md"]


def test_markdown_escaping_and_tsv_quoting():
    result = fixture()
    result["warnings"].append({"sample_id": None, "code": "odd|`code", "count": 1,
                               "detail": 'One\tcell\n<script> & "quoted" | text'})
    rendered = render(result)
    assert "<script>" not in rendered["report.md"]
    assert "&lt;script&gt;" in rendered["report.md"]
    assert "odd\\|\\`code" in rendered["report.md"]
    row = list(csv.DictReader(io.StringIO(rendered["warnings.tsv"]), delimiter="\t"))[-1]
    assert row["detail"] == result["warnings"][-1]["detail"]


def test_one_command_bundle_and_repeated_scientific_hashes(tmp_path):
    a, b = tmp_path / "first", tmp_path / "second"
    first, second = run(EXAMPLE / "request.json", a), run(EXAMPLE / "request.json", b)
    assert first["run_id"] == second["run_id"] and first["outputs"] == second["outputs"]
    assert first["code"] == second["code"]
    assert first["biological_status"] == "abstain"
    for name, info in first["outputs"].items():
        assert checksum(a / name) == info["sha256"]
        assert (a / name).read_bytes() == (b / name).read_bytes()
        assert (a / name).stat().st_mode & 0o777 == 0o600
    assert a.stat().st_mode & 0o777 == 0o700
    assert first["resource_usage"]["peak_resident_bytes"] > 0
    assert "report.md" in first["outputs"] and "candidates.tsv" in first["outputs"]


def test_existing_output_is_preserved(tmp_path, capsys):
    out = tmp_path / "out"
    run(EXAMPLE / "request.json", out)
    before = checksum(out / "report.md")
    assert main(["computational", "--config", str(EXAMPLE / "request.json"), "--output", str(out)]) == 2
    assert "choose a new directory" in capsys.readouterr().err
    assert checksum(out / "report.md") == before


def test_input_mutation_during_render_prevents_commit(tmp_path, monkeypatch, capsys):
    from base_edit_safety.population import computational_reporting
    source = tmp_path / "inputs"
    shutil.copytree(EXAMPLE, source)
    original = computational_reporting.render

    def changed(result):
        with (source / "variants.vcf").open("a") as handle:
            handle.write("\n")
        return original(result)

    monkeypatch.setattr(computational_reporting, "render", changed)
    out = tmp_path / "out"
    assert main(["computational", "--config", str(source / "request.json"), "--output", str(out)]) == 2
    assert "Input changed" in capsys.readouterr().err and not out.exists()


def test_bundle_size_limit_prevents_commit(tmp_path, monkeypatch, capsys):
    from base_edit_safety.population import computational
    monkeypatch.setitem(computational.LIMITS, "max_output_bytes", 100)
    out = tmp_path / "out"
    assert main(["computational", "--config", str(EXAMPLE / "request.json"), "--output", str(out)]) == 4
    assert "max_output_bytes" in capsys.readouterr().err and not out.exists()
