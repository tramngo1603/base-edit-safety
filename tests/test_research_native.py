"""Mock checks for native completion, output writing and resource monitoring."""
import hashlib
import json
from pathlib import Path
import subprocess
import shlex
import sys

import pytest

from base_edit_safety.population.common import PopulationError, ResourceLimit
from base_edit_safety.population import research_native as native

GUIDES = [{"guide_id": "g1", "spacer_sequence": "A"*20, "pam_class": "NGG", "editor_class": "unassessed"}]
CONFIG = {"assembly": "artificial", "reference": {"sha256": "0"*64},
          "backend": {"device": "C"}, "max_mismatches": 2}
LIMITS = {"timeout_seconds": 10, "max_resident_bytes": 2*1024**3,
          "max_output_bytes": 1024**2, "max_native_rows": 100}


def row(dna=None, mm="0", **updates):
    fields = {"query": "A"*20+"NNN", "chrom": "artificial_1", "pos": "0",
              "dna": dna or "A"*20+"AGG", "strand": "+", "mm": mm, "guide": "g1"}
    fields.update(updates)
    return "\t".join(fields.values())+"\n"


def sealed(tmp_path, content):
    raw = tmp_path/"raw.tsv"
    raw.write_text(content)
    native.save_json(tmp_path/"execution.json", {"status": "complete", "exit_code": 0,
        "version": "2.4.1", "raw_bytes": raw.stat().st_size,
        "raw_sha256": hashlib.sha256(raw.read_bytes()).hexdigest()})
    return raw


def test_unknown_spacer_and_pam_follow_exact_writer_rules(tmp_path):
    raw = sealed(tmp_path, row("n"+"A"*19+"AGG", "1") + row("A"*20+"NGG") + row())
    rows, unknown, audit = native.parsed_rows(raw, GUIDES, CONFIG, LIMITS)
    assert len(rows) == 1 and len(unknown) == 2
    assert audit["raw_rows"] == 3
    assert json.loads(rows[0]["source_metadata_json"])["native_line"] == 3
    assert raw.read_text().startswith("A"*20+"NNN")


@pytest.mark.parametrize("updates", [
    {"query": "BOGUS"}, {"guide": "missing"}, {"strand": "?"}, {"pos": "-1"},
    {"mm": "-9"}, {"mm": "0"}, {"dna": "N"+"A"*19+"AGG"},
    {"dna": "n"+"A"*19+"ANg"}, {"dna": "N"}, {"chrom": "white space"},
])
def test_malformed_unknown_rows_are_rejected(tmp_path, updates):
    raw = sealed(tmp_path, row(**{"dna": "n"+"A"*19+"AGG", "mm": "1", **updates}))
    with pytest.raises(PopulationError):
        native.parsed_rows(raw, GUIDES, CONFIG, LIMITS)


def test_duplicate_unknown_rows_are_rejected(tmp_path):
    raw = sealed(tmp_path, row("A"*20+"NGG")*2)
    with pytest.raises(PopulationError, match="Duplicate"):
        native.parsed_rows(raw, GUIDES, CONFIG, LIMITS)


def test_zero_output_requires_unchanged_successful_receipt(tmp_path):
    raw = tmp_path/"raw.tsv"
    raw.write_text("")
    with pytest.raises(PopulationError, match="receipt"):
        native.parsed_rows(raw, GUIDES, CONFIG, LIMITS)
    sealed(tmp_path, "")
    assert native.parsed_rows(raw, GUIDES, CONFIG, LIMITS)[0] == []
    raw.write_text(row())
    with pytest.raises(PopulationError, match="unchanged"):
        native.parsed_rows(raw, GUIDES, CONFIG, LIMITS)


def test_parser_caps_and_callback_are_applied(tmp_path):
    raw = sealed(tmp_path, row()+row(pos="1"))
    with pytest.raises(ResourceLimit, match="row count"):
        native.parsed_rows(raw, GUIDES, CONFIG, {**LIMITS, "max_native_rows": 1})
    checks = []
    native.parsed_rows(raw, GUIDES, CONFIG, LIMITS, check_limits=lambda: checks.append(True))
    assert len(checks) >= 8
    raw = sealed(tmp_path, "X"*65537)
    with pytest.raises(PopulationError, match="line length"):
        native.parsed_rows(raw, GUIDES, CONFIG, LIMITS)


def fake_binary(tmp_path, body):
    binary = tmp_path/"artificial-native"
    binary.write_text("#!/bin/sh\n'''exec' " + shlex.quote(sys.executable) + " \"$0\" \"$@\"\n' '''\nimport sys,time\nfrom pathlib import Path\n"
                      "if len(sys.argv)==1:\n print('Cas-OFFinder v2.4.1 (artificial test)')\n sys.exit(0)\n"+body+"\n")
    binary.chmod(0o700)
    fasta = tmp_path/"artificial.fa"
    fasta.write_text(">artificial_1\n"+"A"*40+"\n")
    return binary, fasta


def test_runner_seals_successful_zero_output(tmp_path, monkeypatch):
    monkeypatch.setattr(native, "_process_group_rss", lambda group: 1024)
    binary, fasta = fake_binary(tmp_path, "Path(sys.argv[3]).write_text('')")
    raw, audit = native.run_native(fasta, tmp_path/"job", GUIDES, CONFIG, binary, LIMITS)
    assert audit["status"] == "complete" and audit["exit_code"] == 0
    assert audit["peak_resident_bytes"] > 0
    assert audit["version"] == "2.4.1"
    assert native.parsed_rows(raw, GUIDES, CONFIG, LIMITS)[0] == []


def test_failed_process_retains_execution_receipt(tmp_path, monkeypatch):
    monkeypatch.setattr(native, "_process_group_rss", lambda group: 0)
    binary, fasta = fake_binary(tmp_path, "sys.exit(7)")
    with pytest.raises(PopulationError, match="failed"):
        native.run_native(fasta, tmp_path/"job", GUIDES, CONFIG, binary, LIMITS)
    audit = json.loads((tmp_path/"job/execution.json").read_text())
    assert audit["status"] == "failed" and audit["exit_code"] == 7
    assert "error" in audit


def test_watchdog_retains_failure_receipt_and_stops_process(tmp_path, monkeypatch):
    monkeypatch.setattr(native, "_process_group_rss", lambda group: LIMITS["max_resident_bytes"]+1)
    binary, fasta = fake_binary(tmp_path, "time.sleep(30)")
    with pytest.raises(ResourceLimit, match="Running native"):
        native.run_native(fasta, tmp_path/"job", GUIDES, CONFIG, binary, LIMITS)
    audit = json.loads((tmp_path/"job/execution.json").read_text())
    assert audit["status"] == "failed"
    assert audit["sampled_process_group_peak_bytes"] > LIMITS["max_resident_bytes"]


def test_parent_callback_failure_is_retained_before_launch(tmp_path):
    binary, fasta = fake_binary(tmp_path, "time.sleep(30)")
    def stop():
        raise ResourceLimit("parent wall cap")
    with pytest.raises(ResourceLimit, match="parent wall"):
        native.run_native(fasta, tmp_path/"job", GUIDES, CONFIG, binary, LIMITS, check_limits=stop)
    audit = json.loads((tmp_path/"job/execution.json").read_text())
    assert audit["status"] == "failed" and audit["exit_code"] is None


def test_process_monitor_sums_only_native_process_group(monkeypatch):
    monkeypatch.setattr(native.subprocess, "run", lambda *a, **kw: subprocess.CompletedProcess(a, 0, "20 100\n20 300\n21 900\n", ""))
    assert native._process_group_rss(20) == 400*1024


def test_version_banner_survives_usage_device_listing_failure(tmp_path, monkeypatch):
    monkeypatch.setattr(native, "_process_group_rss", lambda group: 0)
    binary, fasta = fake_binary(tmp_path, "Path(sys.argv[3]).write_text('')")
    binary.write_text(binary.read_text().replace("sys.exit(0)", "sys.exit(1)"))
    raw, audit = native.run_native(fasta, tmp_path/"job", GUIDES, CONFIG, binary, LIMITS)
    assert audit["version_probe_exit_code"] == 1
    assert audit["exit_code"] == 0 and audit["status"] == "complete"


def test_native_timeout_and_output_caps_retain_receipts(tmp_path, monkeypatch):
    monkeypatch.setattr(native, "_process_group_rss", lambda group: 0)
    binary, fasta = fake_binary(tmp_path, "Path(sys.argv[3]).write_text('X'*20000);time.sleep(30)")
    with pytest.raises(ResourceLimit, match="artifacts"):
        native.run_native(fasta, tmp_path/"large", GUIDES, CONFIG, binary, {**LIMITS, "max_output_bytes": 10000})
    assert json.loads((tmp_path/"large/execution.json").read_text())["status"] == "failed"
    binary, fasta = fake_binary(tmp_path, "time.sleep(30)")
    with pytest.raises(ResourceLimit, match="timeout"):
        native.run_native(fasta, tmp_path/"slow", GUIDES, CONFIG, binary, {**LIMITS, "timeout_seconds": 1})
    assert json.loads((tmp_path/"slow/execution.json").read_text())["status"] == "failed"


def test_parser_counts_retained_artifacts_against_output_cap(tmp_path):
    raw = sealed(tmp_path, row())
    existing = sum(p.stat().st_size for p in tmp_path.iterdir())
    with pytest.raises(ResourceLimit, match="parsing artifacts"):
        native.parsed_rows(raw, GUIDES, CONFIG, {**LIMITS, "max_output_bytes": existing})


@pytest.mark.parametrize("target,reason", [
    ("cccccacGAcCcctaccgGtTGG", "mismatch count"),
    ("aagAcTA              ", "23 ACGTN"),
])
def test_observed_unmodified_gpu_discrepancies_fail_closed(tmp_path, target, reason):
    # These rows reproduce a failed search on invented multi-contig input.
    guide = {**GUIDES[0], "spacer_sequence": "TTAATTAGAACATCGGAAGG"}
    raw = sealed(tmp_path, row(query=guide["spacer_sequence"]+"NNN", dna=target, mm="0"))
    with pytest.raises(PopulationError, match=reason):
        native.parsed_rows(raw, [guide], {**CONFIG, "max_mismatches": 0}, LIMITS)
