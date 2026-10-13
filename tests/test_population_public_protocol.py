"""Check fixed policies and crop boundaries with mock acquisition records."""
from copy import deepcopy
import hashlib
from pathlib import Path
import sys
import types

import pytest

from base_edit_safety.population.common import PopulationError, checksum, json_text, read_json
from base_edit_safety.population.contracts import validate_config
from base_edit_safety.population.public_data import load_protocol, prepare, require_normalization_context, scientific_config

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT/"examples/population/public-protocol.json"
TEMPLATE = ROOT/"examples/population/synthetic/study-config.json"


@pytest.mark.parametrize("section,field,value", [
    ("planning_policy", "default_weight", 7), ("planning_policy", "default_required_donors", 2),
    ("planning_policy", "max_donors", 1), ("planning_policy", "max_assays", 1),
    ("benchmark_policy", "selectors", ["canonical_order"]), ("benchmark_policy", "seeds", [7]),
    ("genotype_policy", "accepted_site_filters", ["PASS"]), ("genotype_policy", "accepted_sample_filters", ["PASS"]),
    ("genotype_policy", "min_dp", 30), ("genotype_policy", "min_gq", 60),
    ("genotype_policy", "phase", "explicit_ps_only"), ("variant_policy", "max_allele_length", 1),
    ("assay_policy", "allele_resolved", False), ("assay_policy", "feasibility", "infeasible"),
    ("inputs", "candidate_ledger", "unrelated-ledger.json"),
])
def test_template_cannot_change_frozen_science(tmp_path, section, field, value):
    protocol = load_protocol(PROTOCOL)
    template = read_json(TEMPLATE)
    template[section][field] = value
    path = tmp_path/"template.json"; path.write_text(json_text(template))
    c = scientific_config(protocol, path, "bcftools test")
    assert c == scientific_config(protocol, TEMPLATE, "bcftools test")
    assert c["planning_policy"]["default_weight"] == c["planning_policy"]["default_required_donors"] == 1
    assert c["genotype_policy"]["min_dp"] is c["genotype_policy"]["min_gq"] is None
    assert c["genotype_policy"]["phase"] == "vcf_implicit_phase_set"
    assert c["benchmark_policy"]["selectors"] == ["greedy_cost", "greedy_coverage", "canonical_order", "random_order"]
    assert c["benchmark_policy"]["seeds"] == list(range(100))


def test_execution_caps_remain_adjustable_but_protocol_mutations_fail(tmp_path):
    c = read_json(TEMPLATE); c["resource_limits"]["max_hit_rows"] = 20
    path = tmp_path/"template.json"; path.write_text(json_text(c))
    assert scientific_config(load_protocol(PROTOCOL), path, "bcftools test")["resource_limits"]["max_hit_rows"] == 20
    protocol = read_json(PROTOCOL); protocol["weight"] = 7
    path.write_text(json_text(protocol))
    with pytest.raises(PopulationError, match="unchanged"):
        load_protocol(path)


def artificial_acquisition(output):
    output.mkdir()
    protocol = read_json(PROTOCOL); lo, hi = protocol["reference_acquisition_interval"]
    seq = "A"*(hi-lo)
    for start in protocol["intended_starts_0based"].values():
        offset = start-lo; seq = seq[:offset]+protocol["guide_sequence"]+"AGG"+seq[offset+23:]
    (output/"frozen-protocol.json").write_bytes(PROTOCOL.read_bytes())
    (output/"reference-acquired.json").write_text(json_text({"record_type": "reference_windows", "schema_version": "1", "genome_build": "GRCh38",
        "windows": [{"contig": "chr11", "start": lo, "end": hi, "sequence": seq, "sequence_sha256": hashlib.sha256(seq.encode()).hexdigest()}]}))
    samples = ["SYN_"+str(i) for i in range(8)]
    (output/"selected-samples.json").write_text(json_text(samples))
    (output/"bounded-original.vcf").write_text('##fileformat=VCFv4.3\n##contig=<ID=chr11,length=135086622>\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+"\t".join(samples)+"\n"
        +f'chr11\t{lo+251}\tSYN_deletion\tAA\tA\t.\tPASS\t.\tGT\t'+"\t".join(["0|1"]*8)+"\n")
    (output/"acquisition.json").write_text(json_text({"files": {p.name: checksum(p) for p in output.iterdir()}}))


def test_prepare_rejects_normalization_at_crop_edge(tmp_path, monkeypatch):
    output = tmp_path/"acquired"; artificial_acquisition(output)
    pysam = types.ModuleType("pysam"); pysam.bcftools = types.ModuleType("pysam.bcftools"); pysam.version = types.ModuleType("pysam.version")
    pysam.faidx = lambda *args: None; pysam.version.__bcftools_version__ = "synthetic_stub"
    def norm(*args, **kwargs):
        Path(args[args.index("-o")+1]).write_text("normalization_window\t1\tsource_record_1\tAA\tA\t.\tPASS\t.\tGT\t0|1\n")
    pysam.bcftools.norm = norm
    for name, module in [("pysam", pysam), ("pysam.bcftools", pysam.bcftools), ("pysam.version", pysam.version)]:
        monkeypatch.setitem(sys.modules, name, module)
    with pytest.raises(PopulationError, match="cropped reference boundary"):
        prepare(PROTOCOL, output, TEMPLATE)
    assert read_json(output/"normalization-failure.json")["status"] == "unresolved"
    assert not (output/"variants.vcf").exists() and not (output/"study-config.json").exists()


@pytest.mark.parametrize("ref,alt", [("AA", "A"), ("A", "AA"), ("A", "C,AA")])
def test_boundary_guard_requires_upstream_context_for_each_indel(ref, alt):
    row = ["normalization_window", "1", "SYN", ref, alt]
    with pytest.raises(PopulationError, match="upstream"):
        require_normalization_context(row, 1000)
    require_normalization_context(row, 0)
    require_normalization_context(["normalization_window", "2", "SYN", ref, alt], 1000)
    require_normalization_context(["normalization_window", "1", "SYN", "A", "C"], 1000)
