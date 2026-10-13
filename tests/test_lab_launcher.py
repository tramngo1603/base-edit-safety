"""Checks for local input handling and truthful example completion."""
import json
import csv
from pathlib import Path
import sys
from types import SimpleNamespace

import pytest

from base_edit_safety import lab


@pytest.mark.parametrize("text", ["", "ACGT", "A" * 20 + "\n" + "A" * 20,
                                  "guide_id\tspacer_sequence\n=SUM()\t" + "A" * 20,
                                  "A" * 19 + "N"])
def test_invalid_guides_are_rejected_before_search(text):
    with pytest.raises(ValueError):
        lab.read_guides(text)


def test_plain_sequences_and_named_tsv():
    assert lab.read_guides("acgattcgatgctacgtagc\n")[0]["spacer_sequence"] == "ACGATTCGATGCTACGTAGC"
    assert lab.read_guides("guide_id\tspacer_sequence\nlab_1\tACGATTCGATGCTACGTAGC")[0]["guide_id"] == "lab_1"


def test_run_never_overwrites_existing_folder(tmp_path, monkeypatch):
    monkeypatch.setattr(lab, "preflight", lambda *a, **k: {})
    marker = tmp_path / "keep.txt"
    marker.write_text("preserve this result")
    with pytest.raises(ValueError, match="already exists"):
        lab.run_job("small", tmp_path, tmp_path / "backend")
    assert marker.read_text() == "preserve this result"


def test_corrupt_public_cache_is_not_reused(tmp_path, monkeypatch):
    monkeypatch.setattr(lab, "preflight", lambda *a, **k: {})
    (tmp_path / "reference.fa").write_text("changed")
    config = {"reference": {"path": "reference.fa", "sha256": "0" * 64},
              "guides": {}, "sample": {"variants": {}}}
    (tmp_path / "request.json").write_text(json.dumps(config))
    with pytest.raises(ValueError, match="Cached public input changed"):
        lab.prepare_public(tmp_path, tmp_path / "backend")


@pytest.mark.parametrize("changed_index", ["reference.fa.fai", "sample.vcf.gz.tbi"])
def test_changed_public_indexes_are_not_reused(tmp_path, monkeypatch, changed_index):
    monkeypatch.setattr(lab, "preflight", lambda *a, **k: {})
    for name in ("reference.fa", "guides.tsv", "sample.vcf.gz", "reference.fa.fai", "sample.vcf.gz.tbi"):
        (tmp_path / name).write_text("fixture")
    config = {"reference": lab.pin(tmp_path / "reference.fa"), "guides": lab.pin(tmp_path / "guides.tsv"),
              "sample": {"variants": lab.pin(tmp_path / "sample.vcf.gz")}}
    (tmp_path / "request.json").write_text(json.dumps(config))
    original = lab.checksum
    def digest(path):
        if Path(path).name == changed_index:
            return "0" * 64
        if str(path).endswith(".fai"):
            return "aeaf6f0a845900acea315871a62982711cebc6da5d982f4085b78ad1ba03a0ef"
        return original(path)
    monkeypatch.setattr(lab, "checksum", digest)
    with pytest.raises(ValueError, match="Cached public index changed"):
        lab.prepare_public(tmp_path, tmp_path / "backend")


def test_old_public_cache_replaces_only_legacy_star_masks(tmp_path, monkeypatch):
    monkeypatch.setattr(lab, 'preflight', lambda *a, **k: {})
    for name in ('reference.fa','reference.fa.fai','guides.tsv','sample.vcf.gz','sample.vcf.gz.tbi','backend'):
        (tmp_path/name).write_text('invented fixture')
    config = {'reference':lab.pin(tmp_path/'reference.fa'), 'guides':lab.pin(tmp_path/'guides.tsv'),
        'sample':{'variants':lab.pin(tmp_path/'sample.vcf.gz'), 'excluded_contigs':[
            {'contig':'chr1','reason':'unresolved_spanning_deletion_extent'},
            {'contig':'chr2','reason':'explicit separate input limitation'}]}}
    cached = json.dumps(config)
    (tmp_path/'request.json').write_text(cached)
    original = lab.checksum
    expected = {'.fai':'aeaf6f0a845900acea315871a62982711cebc6da5d982f4085b78ad1ba03a0ef',
                '.tbi':'2ca4fccd341f78714aaaf344b3c89867e8c6057f9ebbc79eca31c730a9702a29'}
    monkeypatch.setattr(lab, 'checksum', lambda p: expected.get(Path(p).suffix) or original(p))
    updated = lab.prepare_public(tmp_path,tmp_path/'backend',progress=lambda _:None)
    assert updated['sample']['excluded_contigs'] == [config['sample']['excluded_contigs'][1]]
    assert 'GIAB HG003' in updated['sample']['ploidy_source']
    assert (tmp_path/'request.json').read_text() == cached


def test_example_mismatch_cannot_report_success(tmp_path, monkeypatch):
    monkeypatch.setattr(lab, "scientific_summary", lambda p: {"reference_loci": 0})
    expected = tmp_path / "expected.json"
    expected.write_text('{"reference_loci": 1}')
    with pytest.raises(ValueError, match="differ from the accepted baseline"):
        lab.verify_example(tmp_path, expected)


def comparison_fixture(folder):
    folder.mkdir()
    state = dict(status="complete", genotype_status="complete_conditional_model", reference_loci=1,
                 sample_loci=2, candidate_configurations=2, changed_sequence_loci=1,
                 locus_changes={"retained": 1, "gained": 1}, outputs={})
    coverage = dict(guide_count=1, searched_contigs=["chr1"], sample_status="complete_conditional_model",
                    sample_scope_exclusions=[], sample_undeclared_ploidy_intervals=[], unsupported_scope=[],
                    reference_native_accounting={"execution_receipt_sha256": "old runtime", "parsed_rows": 1})
    lab.save_json(folder / "run-manifest.json", state)
    lab.save_json(folder / "coverage.json", coverage)
    lab.save_json(folder / "request.json", {"assembly": "GRCh38"})
    candidates = [dict(configuration_id="a", target_sequence="ACGT", sample_copies="[0,1]",
                       phase_status='["explicit_ps"]', coordinate_path='[{"pos0":1},{"pos0":2}]', evidence_links="private"),
                  dict(configuration_id="b", target_sequence="ATGT", sample_copies="[1]",
                       phase_status='["explicit_ps"]', coordinate_path='[{"pos0":1},{"pos0":2}]', evidence_links="private")]
    lab.research.write_tsv(folder / "candidates.tsv", candidates, list(candidates[0]))
    change = dict(sample_configurations='["a","b"]', reference_configurations='["a"]', status="retained")
    lab.research.write_tsv(folder / "locus-changes.tsv", [change], list(change))
    lab.research.write_tsv(folder / "exclusions.tsv", [dict(reason="unknown_phase", start=1)], ["reason", "start"])
    return folder


def rewrite_rows(path, transform):
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields, rows = reader.fieldnames, list(reader)
    lab.research.write_tsv(path, transform(rows), fields)


def test_scientific_comparison_ignores_only_display_order_and_run_receipt(tmp_path):
    output = comparison_fixture(tmp_path / "run")
    before = lab.scientific_summary(output)
    def reorder(rows):
        for row in rows:
            row["sample_copies"] = json.dumps(list(reversed(json.loads(row["sample_copies"]))))
            row["evidence_links"] = ""
        return list(reversed(rows))
    rewrite_rows(output / "candidates.tsv", reorder)
    rewrite_rows(output / "locus-changes.tsv", lambda rows: [{**rows[0], "sample_configurations": '["b","a"]'}])
    coverage = json.loads((output / "coverage.json").read_text())
    coverage["reference_native_accounting"]["execution_receipt_sha256"] = "new runtime"
    lab.save_json(output / "coverage.json", coverage)
    assert lab.scientific_summary(output) == before


@pytest.mark.parametrize("field,value", [
    ("target_sequence", "CCGT"), ("sample_copies", "[1]"), ("phase_status", '["unknown"]'),
    ("coordinate_path", '[{"pos0":2},{"pos0":1}]'),
])
def test_scientific_comparison_preserves_sequence_copy_phase_and_path(tmp_path, field, value):
    output = comparison_fixture(tmp_path / "run")
    before = lab.scientific_summary(output)
    rewrite_rows(output / "candidates.tsv", lambda rows: [{**rows[0], field: value}, *rows[1:]])
    assert lab.scientific_summary(output) != before


def test_scientific_comparison_retains_duplicate_rows_and_coverage_counts(tmp_path):
    output = comparison_fixture(tmp_path / "run")
    before = lab.scientific_summary(output)
    rewrite_rows(output / "exclusions.tsv", lambda rows: rows + rows)
    assert lab.scientific_summary(output)["scientific_rows_sha256"] != before["scientific_rows_sha256"]
    coverage = json.loads((output / "coverage.json").read_text())
    coverage["reference_native_accounting"]["parsed_rows"] = 2
    lab.save_json(output / "coverage.json", coverage)
    assert lab.scientific_summary(output)["coverage_sha256"] != before["coverage_sha256"]


def test_failed_verification_keeps_report_and_manifest_consistent(tmp_path):
    output = comparison_fixture(tmp_path / "run")
    expected = tmp_path / "expected.json"
    expected.write_text('{"reference_loci": 2}')
    state = lab.record_verification(output, expected)
    assert state["status"] == "complete"
    assert state["example_verification"]["status"] == "mismatched"
    assert "Example verification failed" in (output / "report.html").read_text()
    assert state["outputs"]["report.html"]["sha256"] == lab.checksum(output / "report.html")
    assert json.loads((output / "run-manifest.json").read_text()) == state


def test_cli_verification_failure_returns_nonzero_and_result_location(tmp_path, monkeypatch, capsys):
    from base_edit_safety.launcher import main
    monkeypatch.setattr(lab, "run_job", lambda *a, **kw: (tmp_path, {
        "status": "complete", "example_verification": {"status": "mismatched"}}))
    assert main(["--bundle", str(tmp_path), "--run", "small", "--output", str(tmp_path)]) == 2
    assert json.loads(capsys.readouterr().out)["output"] == str(tmp_path)


def test_resource_checks_stop_before_large_work(tmp_path, monkeypatch):
    import psutil
    monkeypatch.setattr(lab.platform, "system", lambda: "Darwin")
    monkeypatch.setattr(lab.platform, "machine", lambda: "arm64")
    monkeypatch.setattr(lab.platform, "mac_ver", lambda: ("14.4.1", "", ""))
    monkeypatch.setattr(lab.sys, "version_info", (3, 12, 12))
    monkeypatch.setattr(lab.os, "access", lambda *a: True)
    monkeypatch.setattr(Path, "exists", lambda p: True)
    monkeypatch.setattr(psutil, "virtual_memory", lambda: SimpleNamespace(total=18*lab.GIB, available=lab.GIB))
    with pytest.raises(ValueError, match="Close other applications"):
        lab.preflight(tmp_path, large=True)
    monkeypatch.setattr(psutil, "virtual_memory", lambda: SimpleNamespace(total=18*lab.GIB, available=8*lab.GIB))
    monkeypatch.setattr(lab.shutil, "disk_usage", lambda p: SimpleNamespace(free=3*lab.GIB))
    with pytest.raises(ValueError, match="need 25 GiB"):
        lab.preflight(tmp_path, large=True)


def test_report_escapes_input_and_keeps_sequence_change_separate(tmp_path):
    from base_edit_safety.lab_report import write_report
    state = dict(status="complete_reference_sample_failed", reference_loci=10, sample_loci=0,
                 candidate_configurations=10, changed_sequence_loci=0, locus_changes={"unresolved": 10})
    coverage = dict(guide_count=1, searched_contigs=["chr1"], sample_status="failed_sample_scope",
                    sample_scope_exclusions=[dict(contig="chr1", reason="<script>bad</script>")],
                    sample_undeclared_ploidy_intervals=[], unsupported_scope=["bulges"])
    write_report(tmp_path, state, coverage, {"assembly": "GRCh38"})
    text = (tmp_path / "report.html").read_text()
    assert "sample analysis incomplete" in text
    assert "<script>" not in text
    assert "&lt;script&gt;" in text
    assert "sequence-changes.tsv" in text and "locus-changes.tsv" in text


def test_own_vcf_requires_explicit_sample_ploidy(tmp_path, monkeypatch):
    import pysam
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n")
    Path(str(ref)+".fai").write_text("chr1\t4\t6\t4\t5\n")
    monkeypatch.setattr(lab, "GRCH38_LENGTHS", {"chr1": 4})
    monkeypatch.setattr(lab, "base_request", lambda *a: {})
    with pytest.raises(ValueError, match="Confirm two copies"):
        lab.prepare_own(tmp_path, ref, guide_text="ACGATTCGATGCTACGTAGC", reference=ref,
                        vcf=tmp_path / "private.vcf", sample="S", diploid=False)


def test_local_vcf_preparation_keeps_original_alt_indexes_and_phase(tmp_path, monkeypatch):
    import pysam
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\n" + "A"*100 + "\n")
    pysam.faidx(str(ref))
    vcf = tmp_path / "original.vcf"
    vcf.write_text('##fileformat=VCFv4.2\n##reference=GRCh38\n##contig=<ID=chr1,length=100>\n'
                   '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                   '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">\n'
                   '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n'
                   'chr1\t20\t.\tA\tC,G\t.\tPASS\t.\tGT:PS\t2|0:20\n')
    original = vcf.read_bytes()
    monkeypatch.setattr(lab, "GRCH38_LENGTHS", {"chr1":100})
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    request = lab.prepare_own(inputs, ref, guide_text="ACGATTCGATGCTACGTAGC", reference=ref,
                              vcf=vcf, sample="S", diploid=True)
    with pysam.VariantFile(request["sample"]["variants"]["path"]) as handle:
        row, = list(handle)
        assert row.alts == ("C", "G")
        assert row.samples["S"]["GT"] == (2,0)
        assert row.samples["S"]["PS"] == 20
        assert row.samples["S"].phased
    assert vcf.read_bytes() == original


@pytest.mark.skipif(sys.platform != "darwin", reason="Native Mac launcher")
def test_relocatable_tk_and_all_input_modes(tmp_path):
    from base_edit_safety.launcher import Launcher, configure_tk
    configure_tk()
    import tkinter as tk
    root = tk.Tk()
    root.withdraw()
    try:
        app = Launcher(root, tmp_path)
        for flow, kind in (("example", "small-reference"), ("example", "small"), ("example", "public"), ("own", None)):
            app.flow.set(flow)
            if kind:
                app.example_kind.set(kind)
            app.update_flow()
            root.update_idletasks()
            packed = app.own_frame if flow == "own" else app.example_frame
            assert packed.winfo_manager() == "pack"
            assert (app.public_note.winfo_manager() == "pack") == (kind == "public")
        app.lab_mode, app.own_kwargs, app.run_labels = "small", {}, {}
        app.show("run")
        assert app.start.cget("text") == "Start analysis"
        app.vcf.set("sample.vcf")
        app.diploid.set(True)
        app.clear_vcf()
        assert not app.vcf.get() and not app.diploid.get()
        import time
        app.busy, app.started = True, time.monotonic()
        app.events.put(("error", ("Example verification failed", "details")))
        app.poll()
        assert not app.busy
        assert app.failure_text.cget("text") == "Example verification failed"
        assert "Keep this window open" not in app.elapsed.cget("text")
        assert str(app.start.cget("state")) == "normal"
    finally:
        root.destroy()
