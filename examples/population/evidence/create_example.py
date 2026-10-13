"""Create a mock source ledger, donor links and proposed assay plan.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
from pathlib import Path
import shutil

from base_edit_safety.pipeline import candidate_panel as cp, normalize_nominations as nn, population
from base_edit_safety.pipeline.donor_panel import execute
from base_edit_safety.population.common import checksum, json_text, read_json, read_jsonl
from base_edit_safety.population.provenance import resolve_study
from base_edit_safety.population.reporting import tsv

EXAMPLES = Path(__file__).resolve().parents[2]


def write_source_ledger(source, candidate, guide, config):
    """Build a schema-3 ledger from mock BED6 and measurement inputs."""
    build = config["reference_policy"]["genome_build"]
    source_guide = deepcopy(guide)
    (source/"source-guides.tsv").write_text(tsv(sorted(cp.GUIDE_FIELDS), [source_guide]))
    locus = {"guide_id": guide["guide_id"], "genome_build": build, "chrom": "SYN_source_control",
             "start_0based": 0, "end_0based_exclusive": 23, "strand": "+", "locus_name": "source_control"}
    (source/"source-intended.tsv").write_text(tsv(sorted(cp.INTENDED_LOCUS_FIELDS), [locus]))
    profile = read_json(EXAMPLES/"candidate-panel/search-profile.json")
    profile.update(name="synthetic_source", max_mismatches=20, max_dna_bulge=0, max_rna_bulge=0,
                   variant_aware=False, reference_id=build, variant_source=None)
    (source/"source-profile.json").write_text(json_text(profile))
    measurement = {"editor": guide["editor"], "guide_id": guide["guide_id"], "genome_build": build,
        "chrom": candidate["contig"], "start_0based": 0, "end_0based_exclusive": 23, "strand": "+",
        "allele": "alt", "allele_sequence": candidate["sequence_genomic_forward"], "allele_resolution": "resolved",
        "allele_evidence": "synthetic resolved allele", "qc_status": "passed", "treated_group": "synthetic treated",
        "control_group": "synthetic control", "biological_context": "synthetic cells", "assay": "synthetic assay",
        "measurement_kind": "editing", "status": "not_detected", "limit_of_detection_value": 0.1,
        "limit_of_detection_unit": "percent", "source_record_id": "SYN_record", "source_version_or_citation": "synthetic"}
    (source/"source-measurements.tsv").write_text(tsv(sorted(cp.MEASUREMENT_FIELDS), [measurement]))
    fasta, index, offset = "", "", 0
    for contig, sequence in [(candidate["contig"], candidate["sequence_genomic_forward"]),
                             ("SYN_source_control", guide["spacer_sequence"]+"TGG")]:
        header = ">"+contig+"\n"
        fasta += header+sequence+"\n"
        index += f"{contig}\t23\t{offset+len(header)}\t23\t24\n"
        offset = len(fasta)
    (source/"source-reference.fa").write_text(fasta)
    (source/"source-reference.fa.fai").write_text(index)
    (source/"source.bed").write_text(candidate["contig"]+"\t0\t23\tSYN_nomination\t0\t+\n")
    status = nn.main(["experimental", "--input", str(source/"source.bed"), "--guides", str(source/"source-guides.tsv"),
        "--intended-loci", str(source/"source-intended.tsv"), "--source-version", "synthetic",
        "--source-command", "not applicable: synthetic fixture", "--reference-fasta", str(source/"source-reference.fa"),
        "--output", str(source/"source-nominations.tsv"), "--provenance-output", str(source/"source-provenance.json"),
        "--format", "bed6", "--guide-id", guide["guide_id"], "--reference-id", build,
        "--source-method", "synthetic experimental assay"])
    if status != 0:
        raise RuntimeError("Artificial nomination import failed")
    ledger = cp.build_ledger(source/"source-profile.json", source/"source-guides.tsv", source/"source-intended.tsv",
        [source/"source-nominations.tsv"], source/"source-measurements.tsv", reference_fasta_path=source/"source-reference.fa",
        source_provenance_paths=[source/"source-provenance.json"])
    (source/"ledger.json").write_text(json_text(ledger))
    site = next(s for s in ledger["candidates"] if s["measurements"])
    link = {k: measurement[k] for k in ("source_version_or_citation", "assay", "source_record_id")}
    link.update(site_id=site["site_id"], source_person_id="SYN_historical", configuration_id=candidate["configuration_id"],
        target_path=candidate["target_path"], locus_resolved=True, source_assertion="synthetic truth",
        biological_context_id=config["assay_policy"]["biological_context_id"], context=config["assay_policy"]["context"])
    (source/"linkage.json").write_text(json_text({"record_type": "measurement_linkage", "schema_version": "1",
        "ledger_sha256": checksum(source/"ledger.json"), "people": [{"source_person_id": "SYN_historical",
        "canonical_donor_id": "SYN_historical", "identity_status": "resolved", "source_assertion": "synthetic person"}],
        "links": [link]}))


def run(config_path, output):
    study = resolve_study(config_path)
    population.build(study, output/"population")
    population.search(study, output/"population", output/"search")
    population.reconcile_stage(study, output/"search", output/"evidence")
    population.actions_stage(study, output/"population", output/"search", output/"evidence", output/"panel-request.json")
    return execute(output/"panel-request.json", output/"planning-manifest.json", "greedy_cost", output/"plan")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True, help="New local directory for all artificial inputs and outputs")
    args = parser.parse_args()
    root = args.output_dir.resolve()
    root.mkdir(mode=0o700, parents=True, exist_ok=False)
    source = root/"inputs"
    shutil.copytree(EXAMPLES/"population/synthetic", source)
    config_path = source/"study-config.json"
    config = read_json(config_path)
    before = run(config_path, root/"without-evidence")
    candidate = next(c for c in read_jsonl(root/"without-evidence/search/candidate-configurations.jsonl")
                     if c["contig"] == "toy_gain")
    guide = cp.load_guides(source/"guides.tsv")[0]
    write_source_ledger(source, candidate, guide, config)
    config["inputs"].update(candidate_ledger="ledger.json", measurement_linkage="linkage.json")
    for key in ("candidate_ledger", "measurement_linkage"):
        config["input_sources"][key] = {**config["input_sources"]["guides"],
                                        "source_id": "synthetic_"+key, "sha256": checksum(source/config["inputs"][key])}
    evidence_config = source/"study-with-evidence.json"
    evidence_config.write_text(json_text(config))
    after = run(evidence_config, root/"with-evidence")
    print(f"Artificial example: {len(before['selected_actions'])} proposed assays before evidence, {len(after['selected_actions'])} after.")
    print("Inspect inputs/ledger.json, inputs/linkage.json, with-evidence/evidence/ and with-evidence/plan/assay-review.tsv.")


if __name__ == "__main__":
    main()
