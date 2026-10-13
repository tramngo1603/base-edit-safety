"""Integrated sequence comparisons using hand-specified mock cases."""
import csv
import gzip
import hashlib
import json
from pathlib import Path

import pytest

from base_edit_safety.population import research
from base_edit_safety.population.common import PopulationError

GUIDE = "ACGATTCGATGCTACGTAGC"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def request(tmp_path, variants="", copies=2):
    seq = "T" * 40 + GUIDE + "TGG" + "T" * 40 + GUIDE + "AGC" + "T" * 40
    ref = tmp_path / "reference.fa"
    ref.write_text(">artificial_1\n" + seq + "\n")
    guides = tmp_path / "guides.tsv"
    guides.write_text("guide_id\tspacer_sequence\ng1\t" + GUIDE + "\n")
    vcf = tmp_path / "sample.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n##reference=artificial\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n" + variants)
    binary = tmp_path / "backend"
    binary.write_text("artificial backend identity")
    binary.chmod(0o700)
    config = {"format_version": "1", "input_kind": "artificial", "assembly": "artificial",
              "reference": {"path": str(ref), "sha256": sha(ref)},
              "guides": {"path": str(guides), "sha256": sha(guides)},
              "backend": {"path": str(binary), "sha256": sha(binary), "version": "2.4.1", "device": "G0"},
              "contigs": ["artificial_1"], "max_mismatches": 6,
              "sample": {"id": "S", "variants": {"path": str(vcf), "sha256": sha(vcf)},
                         "normalization": "hand-specified artificial SNVs",
                         "genotype_policy": {"phase": "explicit_ps", "min_gq": None, "min_dp": None,
                                             "accepted_site_filters": ["PASS", "."], "accepted_sample_filters": ["PASS", "."]},
                         "ploidy_intervals": [{"contig": "artificial_1", "start": 0, "end": len(seq), "copies": copies}]}}
    path = tmp_path / "request.json"
    path.write_text(json.dumps(config))
    return path, config, seq


def install_rows(monkeypatch, reference_rows, sample_rows=()):
    def native(fasta, directory, guides, config, binary, limits, **kwargs):
        directory.mkdir()
        path = directory / "raw.tsv"
        path.write_text("opaque completed fixture")
        return path, {"exit_code": 0, "wall_seconds": 1, "peak_resident_bytes": 1000}
    def parsed(raw, guides, config, limits, **kwargs):
        rows = reference_rows if raw.parent.name == "reference-native" else sample_rows
        if callable(rows):
            rows = rows(raw.parent.parent)
        return list(rows), [], {"rows_normalized": len(rows)}
    monkeypatch.setattr(research, "run_native", native)
    monkeypatch.setattr(research, "parsed_rows", parsed)


def variant_alias(root, copy=0):
    return (root / "variant-search.fa").read_text().splitlines()[copy*2][1:]


def row(contig, start, target, strand="+"):
    return {"chrom": contig, "start_0based": start, "strand": strand, "site_sequence": target, "guide_id": "g1"}


def read_rows(path):
    return list(csv.DictReader(path.open(), delimiter="\t"))


def test_variant_created_site_without_reference_hit(tmp_path, monkeypatch):
    path, config, seq = request(tmp_path, "artificial_1\t126\t.\tC\tG\t.\tPASS\t.\tGT\t1/1\n")
    def sample_rows(root):
        cluster = json.loads(gzip.open(root / "variant-windows.jsonl.gz", "rt").readline())
        offset = 103 - cluster["interval"]["start"]
        return [row(variant_alias(root, 0), offset, GUIDE + "AGG"), row(variant_alias(root, 1), offset, GUIDE + "AGG")]
    install_rows(monkeypatch, [], sample_rows)
    result = research.run(path, tmp_path / "out")
    assert result["reference_loci"] == 0
    assert result["locus_changes"] == {"gained": 1}
    assert result["candidate_configurations"] == 1


def test_sequence_change_is_not_locus_gain(tmp_path, monkeypatch):
    path, config, seq = request(tmp_path, "artificial_1\t41\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\n")
    def sample_rows(root):
        return [row(variant_alias(root, 0), 22, GUIDE + "TGG"), row(variant_alias(root, 1), 22, "C" + GUIDE[1:] + "TGG")]
    install_rows(monkeypatch, [row("artificial_1", 40, GUIDE + "TGG")], sample_rows)
    result = research.run(path, tmp_path / "out")
    assert result["locus_changes"] == {"retained": 1}
    assert result["changed_sequence_loci"] == 1
    assert result["candidate_configurations"] == 2


def test_missing_call_cannot_establish_loss(tmp_path, monkeypatch):
    path, _, _ = request(tmp_path, "artificial_1\t41\t.\tA\tC\t.\tPASS\t.\tGT\t./.\n")
    install_rows(monkeypatch, [row("artificial_1", 40, GUIDE + "TGG")])
    result = research.run(path, tmp_path / "out")
    assert result["locus_changes"] == {"unresolved": 1}
    assert read_rows(tmp_path / "out/exclusions.tsv")[0]["reason"] == "missing"


def test_haploid_scope_keeps_one_copy(tmp_path, monkeypatch):
    path, _, _ = request(tmp_path, copies=1)
    install_rows(monkeypatch, [row("artificial_1", 40, GUIDE + "TGG")])
    result = research.run(path, tmp_path / "out")
    assert result["locus_changes"] == {"retained": 1}
    assert json.loads(read_rows(tmp_path / "out/candidates.tsv")[0]["sample_copies"]) == [0]


def test_same_copy_adjacent_intervals_are_continuous(tmp_path, monkeypatch):
    path, config, seq = request(tmp_path)
    config["sample"]["ploidy_intervals"] = [{"contig": "artificial_1", "start": a, "end": b, "copies": 2} for a,b in [(0,50),(50,len(seq))]]
    path.write_text(json.dumps(config))
    install_rows(monkeypatch, [row("artificial_1", 40, GUIDE + "TGG")])
    assert research.run(path, tmp_path / "out")["locus_changes"] == {"retained": 1}


def test_no_sample_and_existing_output_protection(tmp_path, monkeypatch):
    path, config, _ = request(tmp_path)
    del config["sample"]
    path.write_text(json.dumps(config))
    install_rows(monkeypatch, [row("artificial_1", 40, GUIDE + "TGG")])
    assert research.run(path, tmp_path / "out")["locus_changes"] == {"reference_only_run": 1}
    before = (tmp_path / "out/run-manifest.json").read_bytes()
    with pytest.raises(PopulationError, match="Output must be new"):
        research.run(path, tmp_path / "out")
    assert (tmp_path / "out/run-manifest.json").read_bytes() == before


def test_unbounded_variant_excludes_its_chromosome(tmp_path, monkeypatch):
    path, _, _ = request(tmp_path, "artificial_1\t41\t.\tA\t<DEL>\t.\tPASS\t.\tGT\t0/1\n")
    install_rows(monkeypatch, [row("artificial_1", 40, GUIDE + "TGG")])
    result = research.run(path, tmp_path / "out")
    assert result["status"] == "complete"
    assert result["locus_changes"] == {"unresolved": 1}


def test_resource_gap_clears_candidate_crossing_its_boundary(tmp_path, monkeypatch):
    path, config, _ = request(tmp_path, "artificial_1\t81\t.\tT\tA\t.\tPASS\t.\tGT\t1/1\n")
    config["limits"] = {"max_cluster_bases": 40}
    path.write_text(json.dumps(config))
    install_rows(monkeypatch, [row("artificial_1", 40, GUIDE + "TGG")])
    result = research.run(path, tmp_path / "out")
    assert result["locus_changes"] == {"unresolved": 1}
    assert result["sample_loci"] == 0
    assert json.loads(read_rows(tmp_path / "out/candidates.tsv")[0]["sample_copies"]) == []


def test_later_unbounded_event_removes_sample_only_indel_paths(tmp_path, monkeypatch):
    path, _, _ = request(tmp_path,
        "artificial_1\t41\t.\tA\tAC\t.\tPASS\t.\tGT\t1/1\n"
        "artificial_1\t146\t.\tT\t<DEL>\t.\tPASS\t.\tGT\t1/1\n")
    install_rows(monkeypatch, [row("artificial_1", 40, GUIDE + "TGG")],
        lambda root: [row(variant_alias(root, 0), 23, "C" + GUIDE[1:] + "TGG")])
    result = research.run(path, tmp_path / "out")
    assert result["locus_changes"] == {"unresolved": 1}
    assert result["candidate_configurations"] == 1
    assert result["sample_loci"] == 0


@pytest.mark.parametrize('chrom,copies', [('chrX', 0),('chrX', 1),('chrX', 2),('chrY', 0),('chrY', 1),('chrY', 2)])
def test_named_xy_requires_source_and_uses_only_declared_copies(tmp_path, monkeypatch, chrom, copies):
    import pysam
    path, config, seq = request(tmp_path, copies=copies)
    reference = Path(config['reference']['path'])
    reference.write_text('>'+chrom+'\n'+seq+'\n')
    pysam.faidx(str(reference))
    config['reference']['sha256'] = sha(reference)
    config.update(input_kind='research_sample', assembly='GRCh38', contigs=[chrom])
    config['sample']['ploidy_intervals'][0]['contig'] = chrom
    path.write_text(json.dumps(config))
    with pytest.raises(PopulationError, match='source-supported ploidy_source'):
        research.prepare(path)
    config['sample']['ploidy_source'] = 'Invented test declaration; entire miniature contig is one copy-count region.'
    variants = Path(config['sample']['variants']['path'])
    variants.write_text(variants.read_text().replace('##reference=artificial', '##reference=GRCh38'))
    config['sample']['variants']['sha256'] = sha(variants)
    path.write_text(json.dumps(config))
    install_rows(monkeypatch, [row(chrom, 40, GUIDE+'TGG')])
    result = research.run(path, tmp_path / 'out')
    assert result['locus_changes'] == ({'retained':1} if copies else {'unresolved':1})
    assert json.loads(read_rows(tmp_path / 'out/candidates.tsv')[0]['sample_copies']) == list(range(copies))


def test_matched_allele_recovery_is_separate_from_locus_recovery(tmp_path, monkeypatch):
    path, config, seq = request(tmp_path, "artificial_1\t41\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\n")
    evidence = tmp_path / "evidence.tsv"
    columns = ["case_id", "guide_id", "contig", "start0", "end0", "strand", "target_sequence",
               "stratum", "outcome", "citation", "source_locator", "eligibility", "reason", "sample_id"]
    with evidence.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(columns)
        writer.writerow(["invented_case", "g1", "artificial_1", 40, 63, "+", "C" + GUIDE[1:] + "TGG",
                         "invented", "invented_allele", "invented", "invented", "supported", "", "S"])
    config["evidence"] = {"path": str(evidence), "sha256": sha(evidence)}
    path.write_text(json.dumps(config))
    install_rows(monkeypatch, [row("artificial_1", 40, GUIDE + "TGG")],
                 lambda root: [row(variant_alias(root, 0), 22, GUIDE + "TGG"), row(variant_alias(root, 1), 22, "C" + GUIDE[1:] + "TGG")])
    research.run(path, tmp_path / "out")
    result = read_rows(tmp_path / "out/evidence.tsv")[0]
    assert result["reference_locus_nominated"] == "True"
    assert result["sample_locus_nominated"] == "True"
    assert result["reference_recovered"] == "False"
    assert result["sample_sequence_nominated"] == "True"
    assert result["matched_sample_recovery"] == "recovered"
    coverage = json.loads((tmp_path / "out/coverage.json").read_text())
    assert coverage["variant_alt_class_counts"] == {"snv": 1}


def test_artificial_contig_cannot_escape_native_directory(tmp_path):
    path, config, _ = request(tmp_path)
    config["contigs"] = ["artificial_../../outside"]
    path.write_text(json.dumps(config))
    with pytest.raises(PopulationError, match="Unsupported or duplicate"):
        research.prepare(path)


@pytest.mark.parametrize("variants", ["", "artificial_1\t41\t.\tA\tC\t.\tPASS\t.\tGT\t0/0\n"])
def test_sample_contig_mask_prevents_reference_fill_even_without_variants(tmp_path, monkeypatch, variants):
    path, config, seq = request(tmp_path, variants)
    config["sample"]["excluded_contigs"] = [{"contig": "artificial_1", "reason": "unresolved spanning deletion"}]
    path.write_text(json.dumps(config))
    install_rows(monkeypatch, [row("artificial_1", 40, GUIDE + "TGG")])
    result = research.run(path, tmp_path / "out")
    assert result["locus_changes"] == {"unresolved": 1}
    assert result["sample_loci"] == 0
    candidate, = read_rows(tmp_path / "out/candidates.tsv")
    assert json.loads(candidate["sample_copies"]) == []
    coverage = json.loads((tmp_path / "out/coverage.json").read_text())
    assert coverage["sample_accounting_complete"] is True
    assert coverage["sample_ploidy_intervals"] == config["sample"]["ploidy_intervals"]
    assert coverage["sample_undeclared_ploidy_intervals"] == []
    assert coverage["sample_scope_exclusions"] == [{"contig": "artificial_1", "start": 0, "end": len(seq), "reason": "sample_contig_excluded: unresolved spanning deletion"}]


def test_masked_stars_and_gain_variants_are_counted_without_sample_claims(tmp_path, monkeypatch):
    variants = "artificial_1\t41\t.\tA\t*\t.\tPASS\t.\tGT\t0/1\nartificial_1\t126\t.\tC\tG\t.\tPASS\t.\tGT\t1/1\n"
    path, config, seq = request(tmp_path, variants)
    config["sample"]["excluded_contigs"] = [{"contig": "artificial_1", "reason": "unresolved spanning deletion"}]
    path.write_text(json.dumps(config))
    def forbidden_sample_search(root):
        raise AssertionError("Excluded sample contig reached native sample search")
    install_rows(monkeypatch, [row("artificial_1", 40, GUIDE + "TGG")], forbidden_sample_search)
    result = research.run(path, tmp_path / "out")
    assert result["genotype_status"] == "complete_conditional_model"
    assert result["locus_changes"] == {"unresolved": 1}
    coverage = json.loads((tmp_path / "out/coverage.json").read_text())
    assert coverage["sample_accounting"]["source_vcf_rows"] == 2
    assert coverage["sample_accounting"]["excluded_contig_source_rows"] == 2
    assert coverage["variant_alt_class_counts"] == {"spanning_deletion_marker": 1, "snv": 1}
    assert coverage["variant_call_status_counts"] == {"complete": 2}


def test_masked_contig_breakend_remains_explicitly_excluded(tmp_path, monkeypatch):
    path, config, _ = request(tmp_path, "artificial_1\t41\t.\tA\t*\t.\tPASS\t.\tGT\t0/1\nartificial_1\t126\t.\tC\tC]artificial_1:100]\t.\tPASS\t.\tGT\t1/1\n")
    config["sample"]["excluded_contigs"] = [{"contig": "artificial_1", "reason": "unresolved spanning deletion"}]
    path.write_text(json.dumps(config))
    install_rows(monkeypatch, [row("artificial_1", 40, GUIDE + "TGG")])
    result = research.run(path, tmp_path / "out")
    assert result["genotype_status"] == "complete_conditional_model"
    assert result["locus_changes"] == {"unresolved": 1}
    assert json.loads((tmp_path / "out/coverage.json").read_text())["sample_accounting_complete"] is True


@pytest.mark.parametrize("masks", [[{"contig": "artificial_1", "reason": ""}],
                                    [{"contig": "artificial_other", "reason": "outside"}],
                                    [{"contig": "artificial_1", "reason": "x"}] * 2])
def test_sample_contig_masks_require_unique_scope_and_reason(tmp_path, masks):
    path, config, _ = request(tmp_path)
    config["sample"]["excluded_contigs"] = masks
    path.write_text(json.dumps(config))
    with pytest.raises(PopulationError, match="excluded_contigs"):
        research.prepare(path)


def test_masked_matched_evidence_is_scope_unresolved_and_other_donor_unlinked(tmp_path):
    columns = ["case_id", "guide_id", "contig", "start0", "end0", "strand", "target_sequence",
               "stratum", "outcome", "citation", "source_locator", "eligibility", "reason", "sample_id"]
    evidence = tmp_path / "evidence.tsv"
    with evidence.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for sample in ["matched", "other"]:
            writer.writerow(dict(zip(columns, [sample, "g1", "artificial_1", 40, 63, "+", GUIDE + "TGG",
                                                "artificial", "unassessed", "artificial", "fixture", "supported", "", sample])))
    sample = {"id": "matched", "excluded_contigs": [{"contig": "artificial_1", "reason": "unresolved star extent"}]}
    rows = research._join_evidence(evidence, {}, sample, "complete_conditional_model")
    assert [r["matched_sample_recovery"] for r in rows] == ["unresolved_sample_scope", "unlinked"]
    assert all(not r["sample_sequence_nominated"] and not r["sample_locus_nominated"] for r in rows)
