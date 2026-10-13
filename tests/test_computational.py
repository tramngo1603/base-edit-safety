"""Check computational limits with hand-specified mock inputs."""
from copy import deepcopy
import hashlib
import json

import pytest

from base_edit_safety.population.common import PopulationError, ResourceLimit, checksum, json_text
from base_edit_safety.population.computational import compute, validate_request
from base_edit_safety.pipeline.comparison import main


GUIDE = "ACGTTAGCTACGATCGTACA"


def request(tmp_path, *, sequence=None, variants=None, samples=None):
    """Write invented bases; assertions below specify truths without running a search."""
    samples = samples or ["S1", "S2", "S3"]
    sequence = sequence or "T" * 30 + GUIDE + "AGC" + "T" * 25 + GUIDE + "TGG" + "T" * 40
    region = {"contig": "artificial_1", "start": 1000, "end": 1000 + len(sequence)}
    reference = {"record_type": "reference_windows", "schema_version": "1", "genome_build": "artificial",
                 "windows": [{**region, "sequence": sequence, "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest()}]}
    ref = tmp_path / "reference.json"
    ref.write_text(json_text(reference))
    vcf = tmp_path / "variants.vcf"
    rows = variants if variants is not None else [
        "artificial_1\t1036\tspacer\tA\tC\t.\tPASS\t.\tGT:PS\t0|0:10\t0|0:10\t0/1:.",
        "artificial_1\t1053\tcreated\tC\tG\t.\tPASS\t.\tGT:PS\t1|0:10\t./.:.\t0/1:."]
    vcf.write_text("##fileformat=VCFv4.3\n##reference=artificial\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                   + "\t".join(samples) + "\n" + "\n".join(rows) + "\n")
    config = {"format_version": "1", "purpose": "computational_only", "input_kind": "artificial", "assembly": "artificial",
              "region": region, "guide": {"id": "test", "spacer": GUIDE, "pam": "NGG", "max_mismatches": 0},
              "samples": samples, "reference": {"path": "reference.json", "sha256": checksum(ref)},
              "variants": {"path": "variants.vcf", "sha256": checksum(vcf), "normalization": "invented literal alleles"},
              "genotype_policy": {"phase": "explicit_ps", "min_gq": None, "min_dp": None,
                                  "accepted_site_filters": ["PASS"], "accepted_sample_filters": ["PASS", "."]}}
    path = tmp_path / "request.json"
    path.write_text(json_text(config))
    return path


def rewrite(path, change):
    config = json.loads(path.read_text())
    change(config)
    path.write_text(json_text(config))


def rewrite_vcf(path, change):
    vcf = path.parent / "variants.vcf"
    vcf.write_text(change(vcf.read_text()))
    rewrite(path, lambda c: c["variants"].update(sha256=checksum(vcf)))


def test_created_candidate_reference_and_unknowns(tmp_path):
    result, artifacts, _, _ = compute(request(tmp_path))
    created = next(r for r in result["candidates"] if r["sample_id"] == "S1" and r["start0"] == 1030)
    assert (created["end0"], created["strand"], created["target_sequence"], created["pam"]) == (1053, "+", GUIDE, "AGG")
    assert created["change"] == "pam_gained"
    assert created["reference_target_sequence"] == GUIDE + "AGC"
    assert created["local_copies"] == [0]
    assert created["model_carriage"] == "present"
    assert created["evidence_supported_carriage"] == "unknown"
    assert created["target_path"] == [{"kind": "reference", "pos0": p} for p in range(1030, 1053)]
    assert result["summary"]["reference_candidates"] == 1
    assert result["summary"]["modeled_candidates"] == 2
    ref = next(r for r in result["candidates"] if r["sample_id"] == "S1" and r["start0"] == 1078)
    assert ref["change"] == "shared" and ref["pam"] == "TGG"
    for sample in ("S2", "S3"):
        row = next(r for r in result["candidates"] if r["sample_id"] == sample and r["start0"] == 1030)
        assert row["model_carriage"] == "unknown" and row["change"] == "unpaired_unevaluable"
    codes = {w["code"] for w in result["warnings"]}
    assert {"missing", "unresolved_unphased", "reference_callability_unavailable", "donor_linkage_unavailable"} <= codes
    assert result["domain"]["biological_status"] == "abstain"
    assert all(r["evidence_supported_carriage"] == "unknown" for r in result["candidates"])
    raw = [json.loads(line) for line in artifacts["candidate-configurations.jsonl"].splitlines()]
    assert all(r["editability"]["status"] == "not_assessed_computational_only" for r in raw)


def test_repeated_science_is_identical(tmp_path):
    path = request(tmp_path)
    a, first, _, _ = compute(path)
    b, second, _, _ = compute(path)
    assert a == b and first == second


def test_raw_evidence_preserves_existing_record_contracts(tmp_path):
    from base_edit_safety.population.contracts import validate_record
    _, artifacts, _, _ = compute(request(tmp_path))
    for name, content in artifacts.items():
        if name.endswith(".jsonl"):
            for line in content.splitlines():
                validate_record(json.loads(line))


@pytest.mark.parametrize("change,reason", [
    (lambda c: c.update(purpose="biological_prediction"), "computational_only"),
    (lambda c: c.update(input_kind="synthetic_test"), "input_kind"),
    (lambda c: c.update(assembly="GRCh38"), "Artificial inputs"),
    (lambda c: c["guide"].update(pam="NAG"), "NGG only"),
    (lambda c: c["guide"].update(spacer="A" * 19), "20-base"),
    (lambda c: c["guide"].update(max_mismatches=7), "integer bounds"),
    (lambda c: c["guide"].update(max_mismatches=True), "integer bounds"),
    (lambda c: c["region"].update(end=11001), "integer bounds"),
    (lambda c: c.update(samples=["S1"] * 9), "1–8"),
    (lambda c: c.update(samples=["S1", "S1"]), "Duplicate"),
    (lambda c: c["genotype_policy"].update(phase="invent_phase"), "phase policy"),
    (lambda c: c["genotype_policy"].update(min_dp=-1), "integer bounds"),
    (lambda c: c["genotype_policy"].update(accepted_site_filters=["all"]), "dot"),
    (lambda c: c.update(editor="ABE7.10"), "unexpected fields"),
])
def test_invalid_config_fails(tmp_path, change, reason):
    path = request(tmp_path)
    rewrite(path, change)
    with pytest.raises(PopulationError, match=reason):
        compute(path)


@pytest.mark.parametrize("change,reason", [
    (lambda s: s.replace("##reference=artificial", "##reference=GRCh37"), "reference declaration"),
    (lambda s: s.replace("1053\tcreated\tC", "1053\tcreated\tA"), "REF disagrees"),
    (lambda s: s.replace("1|0:10", "2|0:10"), "allele index"),
    (lambda s: s.replace("\tS3\n", "\tS4\n"), "sample IDs"),
])
def test_invalid_vcf_fails_without_report(tmp_path, capsys, change, reason):
    path = request(tmp_path)
    rewrite_vcf(path, change)
    out = tmp_path / "out"
    assert main(["computational", "--config", str(path), "--output", str(out)]) == 2
    assert reason in capsys.readouterr().err and not out.exists()


def test_hash_and_reference_window_bindings(tmp_path):
    path = request(tmp_path)
    rewrite(path, lambda c: c["reference"].update(sha256="0" * 64))
    with pytest.raises(PopulationError, match="checksum differs"):
        compute(path)
    path = request(tmp_path)
    rewrite(path, lambda c: c["region"].update(start=1001))
    with pytest.raises(PopulationError, match="must equal"):
        compute(path)


def test_input_size_and_row_limits(tmp_path):
    path = request(tmp_path)
    vcf = path.parent / "variants.vcf"
    vcf.write_text(" " * (1024 ** 2 + 1))
    with pytest.raises(ResourceLimit, match="max_input_bytes"):
        compute(path)
    path = request(tmp_path)
    rewrite_vcf(path, lambda s: s + s.splitlines()[-1] + "\n" + (s.splitlines()[-1] + "\n") * 999)
    with pytest.raises(ResourceLimit, match="max_variant_rows"):
        compute(path)


def test_real_domain_is_distinct_without_processing_real_data(tmp_path):
    artificial = json.loads(request(tmp_path).read_text())
    for kind in ("public_genome", "research_sample"):
        declaration = deepcopy(artificial)
        declaration.update(input_kind=kind, assembly="GRCh38")
        declaration["region"]["contig"] = "chr11"
        assert validate_request(declaration)["input_kind"] == kind
        declaration["region"]["contig"] = "chrX"
        with pytest.raises(PopulationError, match="primary autosome"):
            validate_request(declaration)
    # Check that declared input kinds match the supplied fixture.
    rewrite(tmp_path / "request.json", lambda c: c.update(input_kind="public_genome"))
    with pytest.raises(PopulationError, match="GRCh38"):
        compute(tmp_path / "request.json")


def test_multiallelic_index_and_minus_coordinates(tmp_path):
    # Reverse complement of GUIDE+AGC is GCT+reverse-complement(GUIDE).
    rc = "TGTACGATCGTAGCTAACGT"
    seq = "T" * 30 + "GCT" + rc + "T" * 30
    path = request(tmp_path, sequence=seq, samples=["S1"], variants=[
        "artificial_1\t1031\tmulti\tG\tA,C\t.\tPASS\t.\tGT:PS\t2|0:1"])
    result, artifacts, _, _ = compute(path)
    hit = next(r for r in result["candidates"] if r["start0"] == 1030)
    assert (hit["strand"], hit["end0"], hit["pam"], hit["change"]) == ("-", 1053, "AGG", "pam_gained")
    assert [p["pos0"] for p in hit["target_path"]] == list(range(1052, 1029, -1))
    variant = json.loads(artifacts["normalized-variants.jsonl"].splitlines()[0])
    assert variant["original_alts"] == ["A", "C"]
    assert variant["calls"]["S1"]["selected_indexes"] == [2, 0]


def test_inserted_path_stays_unpaired(tmp_path):
    # Insert G after the final A of GUIDE+AG, yielding an inserted PAM base.
    path = request(tmp_path, sequence="T" * 30 + GUIDE + "AG" + "T" * 30, samples=["S1"], variants=[
        "artificial_1\t1051\tins\tA\tAG\t.\tPASS\t.\tGT:PS\t1|0:1"])
    result, _, _, _ = compute(path)
    hit = next(r for r in result["candidates"] if r["start0"] == 1030)
    assert hit["target_sequence"] == GUIDE and hit["pam"] == "AGG"
    assert hit["change"] == "indel_path_requires_separate_comparison"
    assert hit["reference_target_sequence"] is None
    assert {p["kind"] for p in hit["target_path"]} == {"reference", "insertion"}
    assert hit["evidence_supported_carriage"] == "unknown"


def test_missing_and_unsupported_calls_visible_with_zero_hits(tmp_path):
    path = request(tmp_path, sequence="T" * 80, samples=["S1"], variants=[
        "artificial_1\t1031\tunknown\tT\t<DEL>\t.\tPASS\t.\tGT:PS\t1|0:1"])
    result, artifacts, _, _ = compute(path)
    assert result["candidates"] == []
    assert result["summary"]["excluded_oriented_windows"] > 0
    assert any(w["code"] == "unresolved_variant_footprint" for w in result["warnings"])
    material = json.loads(artifacts["sequence-configurations.jsonl"].splitlines()[1])
    assert material["excluded_intervals"][0]["reason"] == "unsupported_carried_allele"


def test_single_unphased_site_is_only_an_unordered_pair(tmp_path):
    path = request(tmp_path)
    rewrite_vcf(path, lambda s: s.replace("1|0:10", "1/0:."))
    result, _, _, _ = compute(path)
    row = next(r for r in result["candidates"] if r["sample_id"] == "S1" and r["start0"] == 1030)
    assert row["change"] == "pam_gained" and row["phase_status"] == ["single_heterozygous_site"]


def test_phase_blocks_and_quality_thresholds_never_silently_fill(tmp_path):
    path = request(tmp_path)
    rewrite_vcf(path, lambda s: s.replace("0|0:10", "0|1:20", 1))
    result, _, _, _ = compute(path)
    assert not any(r["start0"] == 1030 for r in result["candidates"])
    assert any(w["sample_id"] == "S1" and w["code"] == "unresolved_phase_blocks" for w in result["warnings"])
    path = request(tmp_path)
    rewrite(path, lambda c: c["genotype_policy"].update(min_dp=20))
    result, artifacts, _, _ = compute(path)
    assert result["summary"]["modeled_candidates"] == 1  # Far reference-fill candidate only.
    calls = [json.loads(line) for line in artifacts["genotype-audit.jsonl"].splitlines()]
    assert all("dp_missing_or_failed" in row["reasons"] for row in calls)
    assert all(r["evidence_supported_carriage"] == "unknown" for r in result["candidates"])


def test_reference_n_is_not_searchable_and_zero_result_has_warning(tmp_path):
    path = request(tmp_path, sequence="N" * 40, variants=[], samples=["S1"])
    result, _, _, _ = compute(path)
    assert not result["candidates"]
    assert result["summary"]["evaluated_oriented_windows"] == 0
    assert result["summary"]["excluded_oriented_windows"] == 3 * 18 * 2
    assert any(w["code"] == "sequence_unknown" for w in result["warnings"])


def test_reference_only_candidate_does_not_confirm_sample_genotype(tmp_path):
    path = request(tmp_path, sequence=GUIDE + "AGG", variants=[], samples=["S1"])
    result, _, _, _ = compute(path)
    row, = result["candidates"]
    assert row["change"] == "shared" and row["local_copies"] == [0, 1]
    assert row["model_carriage"] == "present" and row["evidence_supported_carriage"] == "unknown"
    assert row["phase_status"] == ["no_heterozygote_in_retained_model"]


def test_duplicate_json_fields_rejected(tmp_path):
    path = request(tmp_path)
    path.write_text(path.read_text().replace('"purpose": "computational_only",', '"purpose": "computational_only", "purpose": "artificial",'))
    with pytest.raises(PopulationError, match="Duplicate JSON key"):
        compute(path)


def test_resource_failure_never_writes_completed_report(tmp_path, monkeypatch, capsys):
    from base_edit_safety.population import computational
    path = request(tmp_path)
    monkeypatch.setitem(computational.LIMITS, "max_resident_bytes", 1)
    out = tmp_path / "out"
    assert main(["computational", "--config", str(path), "--output", str(out)]) == 4
    assert '"status": "partial"' in capsys.readouterr().err
    assert not out.exists()
