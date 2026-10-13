"""Native parsing checks and a separate exhaustive alignment oracle.

The oracle uses only Python strings and its own PAM/strand logic. It does not
call production reconstruction, alignment, admission, or identity functions.
"""
from copy import deepcopy
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import shutil

import pytest

from base_edit_safety.population.adapters import crisprme_plus as adapter
from base_edit_safety.population.common import PopulationError, checksum, digest
from base_edit_safety.population.external_search import (NativeLimits, _cached, _commands, run_native,
                                                        native_job_binding, verify_native_job)


FIXTURES = Path(__file__).parent / "fixtures/comparison/external"


def oracle_reverse(sequence):
    return "".join({"A": "T", "T": "A", "C": "G", "G": "C"}[base] for base in reversed(sequence))


def oracle_alignments(guide, target, max_mm, max_dna, max_rna, max_total):
    """Enumerate every column path, including adjacent and terminal gaps."""
    @lru_cache(None)
    def visit(i, j, mm, dna, rna):
        if mm > max_mm or dna > max_dna or rna > max_rna or mm + dna + rna > max_total:
            return ()
        if i == len(guide) and j == len(target):
            return (("", ""),)
        paths = []
        if i < len(guide) and j < len(target):
            paths.extend((guide[i] + a, target[j] + b) for a, b in
                         visit(i + 1, j + 1, mm + (guide[i] != target[j]), dna, rna))
        if j < len(target):
            paths.extend(("-" + a, target[j] + b) for a, b in visit(i, j + 1, mm, dna + 1, rna))
        if i < len(guide):
            paths.extend((guide[i] + a, "-" + b) for a, b in visit(i + 1, j, mm, dna, rna + 1))
        return tuple(paths)
    return set(visit(0, 0, 0, 0, 0))


def exhaustive_candidates(windows, guide, profile):
    codes = {"A": "A", "C": "C", "G": "G", "T": "T", "N": "ACGT", "R": "AG", "Y": "CT",
             "S": "CG", "W": "AT", "K": "GT", "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG"}
    candidates = set()
    pam_lengths = {len(pattern) for pattern in profile["pam_patterns"]}
    assert len(pam_lengths) == 1
    pam_length = next(iter(pam_lengths))
    for window in windows:
        sequence = window["sequence"]
        for strand in profile["strands"]:
            oriented = sequence if strand == "+" else oracle_reverse(sequence)
            for length in range(len(guide) - profile["max_rna_bulges"], len(guide) + profile["max_dna_bulges"] + 1):
                for left in range(len(oriented) - length - pam_length + 1):
                    target = oriented[left:left + length]
                    pam = oriented[left + length:left + length + pam_length]
                    if not any(all(base in codes[code] for base, code in zip(pam, pattern))
                               for pattern in profile["pam_patterns"]):
                        continue
                    for aligned_guide, aligned_target in oracle_alignments(guide, target, profile["max_mismatches"],
                            profile["max_dna_bulges"], profile["max_rna_bulges"], profile["max_total_edits"]):
                        start = left if strand == "+" else len(sequence) - left - length - pam_length
                        end = start + length + pam_length
                        candidates.add((window["window_id"], start, end, strand, aligned_guide, aligned_target, pam))
    return candidates


def synthetic_inputs(kind="ungapped", pam="NGG"):
    data = json.loads((FIXTURES / "domain.json").read_text())
    windows = []
    for case in data["cases"]:
        target = case["spacer"] + case["pam"]
        if case["strand"] == "-":
            target = oracle_reverse(target)
        sequence = data["padding"] + target + data["padding"]
        path, position = [], 1000
        for index, base in enumerate(sequence):
            if index == len(data["padding"]) + case.get("insertion_index", -1000):
                path.append({"kind": "insertion", "boundary0": position, "offset0": 0, "base": base})
                continue
            if index == len(data["padding"]) + case.get("deleted_index", -1000):
                position += 1
            path.append({"kind": "reference", "pos0": position, "base": base})
            position += 1
        variants, genotypes = [], []
        if "selected_alt_index" in case:
            variants = [{"contig": case["locus"], "pos0": 1040, "reference": "G", "alternates": ["C", "A"],
                         "selected_alt_index": 2, "source_artifact_id": "synthetic_genotypes", "source_row": 4}]
            genotypes = [{"original_alleles": ["G", "C", "A"], "selected_indexes": [2, 1], "gt": "2|1",
                          "phase_key": "synthetic_block", "phase_status": "confirmed", "call_status": "complete", "source_record_index": 4}]
        elif "insertion_index" in case:
            variants = [{"contig": case["locus"], "pos0": 1039, "reference": "C", "alternates": ["CA"],
                         "selected_alt_index": 1, "source_artifact_id": "synthetic_genotypes", "source_row": 5}]
            genotypes = [{"original_alleles": ["C", "CA"], "selected_indexes": [1, 0], "gt": "1|0",
                          "phase_key": "synthetic_block", "phase_status": "confirmed", "call_status": "complete", "source_record_index": 5}]
        elif "deleted_index" in case:
            variants = [{"contig": case["locus"], "pos0": 1036, "reference": "CC", "alternates": ["C"],
                         "selected_alt_index": 1, "source_artifact_id": "synthetic_genotypes", "source_row": 6}]
            genotypes = [{"original_alleles": ["CC", "C"], "selected_indexes": [1, 0], "gt": "1|0",
                          "phase_key": "synthetic_block", "phase_status": "confirmed", "call_status": "complete", "source_record_index": 6}]
        windows.append({"window_id": case["name"], "sequence": sequence, "path": path,
                        "source_material_id": "synthetic:" + case["name"], "reference_bundle_id": "synthetic",
                        "contig": case["locus"], "interval": {"contig": case["locus"], "start": 1000, "end": position},
                        "donor_id": "synthetic_donor", "sample_id": "synthetic_donor", "copy_index": 0,
                        "phase_status": "confirmed", "phase_block": "synthetic_block", "source_variants": variants,
                        "source_genotypes": genotypes,
                        "coverage_basis": "synthetic_truth", "support_status": "supported", "source_version": "synthetic-v1",
                        "source_refs": [{"source_artifact_id": "synthetic", "row": len(windows) + 1}]})
    profile = {**data[kind], "pam_patterns": [pam], "strands": ["+", "-"],
               "total_edit_accounting": "mismatch_plus_bulged_bases", "max_rows": 100000,
               "max_window_variants": 20, "retention": "all", "source_equivalence": "not_applicable_synthetic",
               "backend_version": adapter.VERSION, "backend_digest": adapter.IMAGE,
               "backend_name": "crisprme_plus",
               "search_profile_id": "synthetic_" + kind + "_" + pam, "source_panel_id": "synthetic"}
    design = {"guide_sequence": data["guide"], "design_id": "synthetic_guide"}
    editor = {"guide_geometry": {"spacer_length": 20, "pam_side": "3prime", "pam_length": len(pam)}}
    return windows, design, editor, profile


def hit_key(hit):
    return tuple(hit[k] for k in ("window_id", "start", "end", "strand", "guide_alignment", "target_alignment", "pam_sequence"))


def test_native_policy_is_pinned_and_counts_bases():
    _, design, editor, profile = synthetic_inputs()
    adapter.validate_policy(design, editor, profile)
    adapter.validate_policy(design, editor, {**profile, "backend_digest": adapter.IMAGE.rsplit(":", 1)[1]})
    for change in ({"backend_name": "local"}, {"backend_version": "v9"}, {"backend_digest": "sha256:unknown"},
                   {"total_edit_accounting": "mismatch_plus_bulge_events"}):
        with pytest.raises(PopulationError):
            adapter.validate_policy(design, editor, dict(profile, **change))


def test_native_exact_coordinates_and_reverse_path():
    windows, design, editor, profile = synthetic_inputs()
    window = windows[2]
    line = "X\t" + design["guide_sequence"] + "NNN\t" + design["guide_sequence"] + "TGG\tchrW0\t32\t32\t-\t0\t0\t0"
    hit = adapter.parse_row(line, windows={"chrW0": window}, guide=design["guide_sequence"],
                            geometry=editor["guide_geometry"], profile=profile, source_ref={"line": 1})
    assert hit["target_indexes"] == list(range(54, 34, -1))
    assert hit["pam_indexes"] == [34, 33, 32]
    with pytest.raises(PopulationError, match="sequence"):
        adapter.parse_row(line.replace("\t32\t32\t", "\t31\t32\t"), windows={"chrW0": window},
                          guide=design["guide_sequence"], geometry=editor["guide_geometry"], profile=profile, source_ref={})


def test_native_malformed_counts_and_own_pam_admission():
    windows, design, editor, profile = synthetic_inputs()
    line = "X\t" + design["guide_sequence"] + "NNN\t" + design["guide_sequence"] + "AGA\tchrW0\t32\t32\t+\t0\t0\t0"
    assert adapter.parse_row(line, windows={"chrW0": windows[-1]}, guide=design["guide_sequence"],
        geometry=editor["guide_geometry"], profile=profile, source_ref={}) is None
    with pytest.raises(PopulationError, match="counts"):
        adapter.parse_row(line[:-1] + "1", windows={"chrW0": windows[-1]}, guide=design["guide_sequence"],
            geometry=editor["guide_geometry"], profile=profile, source_ref={})


def test_native_retention_cap_is_explicit(tmp_path):
    windows, design, editor, profile = synthetic_inputs()
    profile["max_rows"] = 1
    row = "X\t" + design["guide_sequence"] + "NNN\t" + design["guide_sequence"] + "AGG\tchrW0\t32\t32\t+\t0\t0\t0\n"
    path = tmp_path / "targets.txt"
    path.write_text(row * 2)
    parsed = adapter.parse_targets([path], windows={"chrW0": windows[0]}, design=design, editor=editor, profile=profile)
    assert parsed["retention_status"] == "truncated"
    assert parsed["raw_rows"] == 2
    assert path.read_text() == row * 2


def test_cache_rejects_changed_raw_file(tmp_path):
    path = tmp_path / "raw.txt"
    path.write_text("retained")
    payload = {"fixture": "retained native input"}
    (tmp_path / "input.json").write_text(json.dumps(payload))
    key = digest(payload)
    manifest = {"stage_key": key, "execution_status": "complete", "retention_status": "complete",
                "artifact_manifest": {name: {"sha256": checksum(tmp_path / name)} for name in ("raw.txt", "input.json")}}
    (tmp_path / "result.json").write_text(json.dumps(manifest))
    (tmp_path / "COMPLETE.json").write_text(json.dumps({"stage_key": key, "result_sha256":
        hashlib.sha256((tmp_path / "result.json").read_bytes()).hexdigest()}))
    assert _cached(tmp_path, key)["reused"]
    path.write_text("changed")
    with pytest.raises(PopulationError, match="changed"):
        _cached(tmp_path, key)


def test_native_arguments_keep_explicit_gaps_and_retention():
    _, design, _, profile = synthetic_inputs("gapped")
    commands = _commands(design, profile, "NGG", 2)
    assert commands[0][1] == "build-index-only"
    assert "--ci-cd-test" in commands[-1]
    for flag, expected in (("--mm", "1"), ("--bDNA", "1"), ("--bRNA", "1"), ("--max-total-edits", "2")):
        assert commands[-1][commands[-1].index(flag) + 1] == expected
    with pytest.raises(PopulationError):
        NativeLimits(cpus=5).validate()


def test_independent_oracle_retains_alternative_and_adjacent_gaps():
    alignments = oracle_alignments("AA", "AA", 0, 1, 1, 2)
    assert ("AA", "AA") in alignments
    assert ("A-A", "AA-") in alignments
    assert ("-AA", "A-A") in alignments
    assert ("AA-", "A-A") in alignments
    assert len(alignments) > 4


def test_distinct_loci_and_genomic_insertion_paths_stay_separate():
    windows, design, _, profile = synthetic_inputs("gapped")
    candidates = exhaustive_candidates(windows, design["guide_sequence"], profile)
    exact = [row for row in candidates if row[4] == row[5] == design["guide_sequence"] and row[3] == "+"]
    assert {"exact_plus", "same_sequence_other_locus"} <= {row[0] for row in exact}
    insertion = next(w for w in windows if w["window_id"] == "genomic_insertion")
    assert len([p for p in insertion["path"] if p["kind"] == "insertion"]) == 1
    assert any(row[0] == insertion["window_id"] and "-" in row[4] for row in candidates)


@pytest.mark.parametrize("kind,pam,filename", [
    ("ungapped", "NGG", "ungapped-ngg.targets.tsv"),
    ("gapped", "NGG", "gapped-ngg.targets.tsv"),
    ("ungapped", "NG", "ungapped-ng.targets.tsv")])
def test_actual_retained_native_rows_against_independent_complete_sets(kind, pam, filename):
    windows, design, editor, profile = synthetic_inputs(kind, pam)
    aliases = {f"chrW{i:08d}": w for i, w in enumerate(sorted(windows, key=lambda row: row["window_id"]))}
    raw = adapter.parse_targets([FIXTURES / filename], windows=aliases, design=design, editor=editor, profile=profile)
    expected = exhaustive_candidates(windows, design["guide_sequence"], profile)
    if kind == "gapped":
        missing = expected - {hit_key(h) for h in raw["hits"]}
        assert len(missing) == 9 and all(row[4].startswith("-") for row in missing)
    result = adapter.extend_leading_dna(raw, windows=aliases, profile=profile)
    assert {hit_key(h) for h in result["hits"]} == expected
    for hit in result["hits"]:
        if hit["alignment_origin"] == "derived":
            assert hit["derivation"]["version"] == adapter.LEADING_DNA_TRANSFORM
            assert all("Derived alignment" in ref["note"] for ref in hit["native_row_refs"])


def test_high_density_is_partial_and_input_metadata_changes_job_key(tmp_path):
    windows, design, editor, profile = synthetic_inputs()
    windows = windows[:1]
    windows[0]["source_genotypes"] = [{"gt": "0/0"}]
    profile["max_window_variants"] = 0
    kwargs = {"output_dir": tmp_path, "source_equivalence": {"status": "not_applicable_synthetic"}}
    first = run_native(windows, design, editor, profile, **kwargs)
    assert first["execution_status"] == "partial"
    assert first["evaluated_window_ids"] == [] and first["hits"] == []
    assert "high_density_windows_skipped" in first["reason_codes"]
    assert first["skipped_windows"][0]["interval"] == windows[0]["interval"]
    assert verify_native_job(native_job_binding(first))["execution_status"] == "partial"
    repeated = run_native(windows, design, editor, profile, **kwargs)
    assert repeated["stage_key"] == first["stage_key"] and not repeated["reused"]
    assert repeated["artifact_path"] != first["artifact_path"]
    changed = deepcopy(windows)
    changed[0]["phase_block"] = "different_source_phase"
    second = run_native(changed, design, editor, profile, **kwargs)
    assert second["stage_key"] != first["stage_key"]
    third = run_native(windows, design, editor, {**profile, "max_total_edits": 2}, **kwargs)
    assert third["stage_key"] != first["stage_key"]


def test_unsupported_verified_donor_claim_fails_before_execution(tmp_path):
    windows, design, editor, profile = synthetic_inputs()
    with pytest.raises(PopulationError, match="source reconstruction audit"):
        run_native(windows, design, editor, profile, output_dir=tmp_path,
                   source_equivalence={"status": "verified"})


def test_guard_coordinates_are_rebased_and_generated_dna_is_excluded():
    windows, design, editor, profile = synthetic_inputs()
    window = deepcopy(windows[0])
    window["native_guard"] = {"prefix": design["guide_sequence"] + "AGG", "suffix": "AAAA"}
    row = "X\t" + design["guide_sequence"] + "NNN\t" + design["guide_sequence"] + "AGG\tchrW0\t55\t55\t+\t0\t0\t0"
    hit = adapter.parse_row(row, windows={"chrW0": window}, guide=design["guide_sequence"],
                            geometry=editor["guide_geometry"], profile=profile, source_ref={})
    assert hit["start"] == 32 and hit["end"] == 55
    assert hit["native_coordinate_transform"]["subtract_offset"] == 23
    guard_row = row.replace("\t55\t55\t", "\t0\t0\t")
    assert adapter.parse_row(guard_row, windows={"chrW0": window}, guide=design["guide_sequence"],
                            geometry=editor["guide_geometry"], profile=profile, source_ref={}) is None


def test_stage_binding_detects_simultaneously_resealed_native_result(tmp_path):
    payload = {"fixture": "native audit checksum linkage"}
    (tmp_path / "input.json").write_text(json.dumps(payload))
    result = {"stage_key": digest(payload), "execution_status": "complete", "retention_status": "complete",
              "artifact_manifest": {"input.json": {"sha256": checksum(tmp_path / "input.json")}}, "hits": []}
    result_path, seal = tmp_path / "result.json", tmp_path / "COMPLETE.json"
    result_path.write_text(json.dumps(result))
    seal.write_text(json.dumps({"stage_key": result["stage_key"], "result_sha256": checksum(result_path)}))
    binding = native_job_binding({**result, "artifact_path": str(tmp_path)})
    assert verify_native_job(binding)["hits"] == []
    result["hits"] = [{"changed": True}]
    result_path.write_text(json.dumps(result))
    seal.write_text(json.dumps({"stage_key": result["stage_key"], "result_sha256": checksum(result_path)}))
    with pytest.raises(PopulationError, match="enclosing-stage checksum"):
        verify_native_job(binding)


@pytest.mark.parametrize("mutation", ["raw", "result_and_seal"])
def test_native_tampering_invalidates_search_and_descendant_reuse(tmp_path, mutation):
    from base_edit_safety.population.comparison_provenance import resolve_study, save_stage, load_stage
    source = Path(__file__).resolve().parents[1] / "examples/comparison/contracts"
    shutil.copytree(source, tmp_path / "inputs")
    study = resolve_study(tmp_path / "inputs/study.json", code={"test": "native evidence chain"})
    native_dir = tmp_path / "native"
    native_dir.mkdir()
    payload = {"fixture": "interface-only native audit binding"}
    (native_dir / "input.json").write_text(json.dumps(payload))
    (native_dir / "raw.txt").write_text("retained native bytes")
    result = {"stage_key": digest(payload), "execution_status": "complete", "retention_status": "complete",
              "artifact_manifest": {name: {"sha256": checksum(native_dir / name)} for name in ("input.json", "raw.txt")}, "hits": []}
    result_path, seal = native_dir / "result.json", native_dir / "COMPLETE.json"
    result_path.write_text(json.dumps(result))
    seal.write_text(json.dumps({"stage_key": result["stage_key"], "result_sha256": checksum(result_path)}))
    binding = native_job_binding({**result, "artifact_path": str(native_dir)})
    execution = {"execution_status": "complete", "retention_status": "complete", "exit_code": 0,
                 "native_manifest": [binding]}
    search, observations, evidence = (tmp_path / p for p in ("search", "observations", "evidence"))
    save_stage(study, "search", search, {"search-audit.json": {"execution": [execution]}})
    save_stage(study, "observation_normalization", observations, {"check.txt": "measurement interface fixture"})
    save_stage(study, "evidence_applicability", evidence, {"check.txt": "derived interface fixture"},
               parents={"search": search, "observation_normalization": observations})
    load_stage(study, evidence, "evidence_applicability")
    if mutation == "raw":
        (native_dir / "raw.txt").write_text("changed")
    else:
        result["hits"] = [{"tampered": True}]
        result_path.write_text(json.dumps(result))
        seal.write_text(json.dumps({"stage_key": result["stage_key"], "result_sha256": checksum(result_path)}))
    for path, stage in ((search, "search"), (evidence, "evidence_applicability")):
        with pytest.raises(PopulationError, match="changed"):
            load_stage(study, path, stage)
