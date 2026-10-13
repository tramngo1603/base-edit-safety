"""Explicit acquisition and normalization for the frozen, bounded public protocol."""
from __future__ import annotations

import argparse
import hashlib
from importlib.resources import files as package_files
import json
from pathlib import Path
import sys
import urllib.request

from .common import PopulationError, checksum, identity, json_text, parse_json, read_json, write_file

BASE = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/"
VCF = "1kGP_high_coverage_Illumina.chr11.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
PEDIGREE = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/1kGP.3202_samples.pedigree_info.txt"


def load_protocol(path):
    """This acquisition helper implements one versioned, frozen study."""
    protocol = read_json(path)
    frozen = parse_json(package_files(__package__).joinpath("protocols", "hbg-bounded-nygc-20220422-v1.json").read_text())
    if protocol != frozen:
        raise PopulationError("This helper requires the unchanged hbg-bounded-nygc-20220422-v1 protocol")
    return protocol


def scientific_config(protocol, example_config, normalizer):
    """Inherit execution limits only; all scientific choices belong to the frozen protocol."""
    context = {"cell_type": "primary_CD34_HSPC_scenario", "delivery": "electroporation_scenario",
               "condition": "declared_research_scenario_not_observed"}
    from .contracts import resolution_id
    assay = {"biological_context_id": identity("context", context), "specimen_context_id": "scenario_hspc",
        "measurement_class": "editing", "assay_id": "scenario_allele_locus_resolved_amplicon",
        "allele_resolved": True, "locus_resolved": True, "feasibility": "declared_feasible",
        "assay_cost": protocol["assay_cost"], "control_costs_included": False, "extra_context_bases": 0,
        "context": context, "guide_ids": ["HBG_m198_ABE8e"],
        "locus_scope": {"mode": "all_admitted_loci_simulation", "loci": []}}
    assay["resolution_contract_id"] = resolution_id(assay)
    interval = {"contig": protocol["contig"], "start": protocol["requested_interval"][0], "end": protocol["requested_interval"][1]}
    return {"record_type": "population_study_config", "schema_version": "1", "study_id": protocol["study_id"],
        "mode": protocol["mode"], "coverage_basis": protocol["coverage_basis"],
        "inputs": {"reference": "reference-windows.json", "variants": "variants.vcf", "donors": "donors.tsv",
                   "guides": "guides.tsv", "intended_loci": "intended-loci.tsv", "callable_intervals": None,
                   "candidate_ledger": None, "measurement_linkage": None}, "input_sources": {},
        "reference_policy": {"genome_build": protocol["assembly"], "contig_alias_map": {protocol["contig"]: protocol["contig"]},
                             "permitted_contigs": [protocol["contig"]]},
        "variant_policy": {"normalizer": normalizer, "max_allele_length": 50, "unsupported_footprint": "mark_affected_scope_unknown"},
        "genotype_policy": {"accepted_site_filters": ["PASS", "."], "accepted_sample_filters": [".", "PASS"],
            "reference_evidence": "reference_imputed_scenario", "min_dp": None, "min_gq": None,
            "phase": "vcf_implicit_phase_set", "ploidy": 2},
        "scope": {"requested_intervals": [interval], "materialized_intervals": [interval]},
        "search_profile": {**protocol["search"], "spacer_length": 20}, "assay_policy": assay,
        "planning_policy": {"cost_unit": protocol["cost_unit"], "max_cost": max(protocol["study_budget_grid"]),
            "fixed_overhead": protocol["fixed_control_overhead"], "default_weight": protocol["weight"],
            "default_required_donors": protocol["required_distinct_donors"], "max_donors": None, "max_assays": None,
            "algorithm": "greedy_cost"},
        "benchmark_policy": {"budget_grid": protocol["study_budget_grid"], "relatedness_policy": protocol["relatedness_policy"],
            "selectors": protocol["selectors"], "seeds": list(range(protocol["seed_range_inclusive"][0], protocol["seed_range_inclusive"][1]+1)),
            "endpoints": protocol["endpoints"], "split": "transductive"},
        "resource_limits": read_json(example_config)["resource_limits"]}


def require_normalization_context(row, window_start):
    """A left-aligned indel at a crop edge may still move with more upstream reference."""
    if window_start > 0 and int(row[1]) == 1 and any(len(alt) != len(row[3]) for alt in row[4].split(",")):
        raise PopulationError("Normalization reached the cropped reference boundary; acquire more genuine upstream reference context before preparing these inputs")


def retrieve(url, path, limit=10000000):
    with urllib.request.urlopen(url, timeout=30) as response:
        data = response.read(limit+1)
    if len(data) > limit:
        raise PopulationError("Public metadata response exceeded acquisition limit")
    path.write_bytes(data)
    path.chmod(0o600)
    return data


def fetch_window(contig, lo, hi, output):
    url = f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom={contig};start={lo};end={hi}"
    payload = json.loads(retrieve(url, output))
    sequence = payload["dna"].upper()
    if len(sequence) != hi-lo or set(sequence)-set("ACGTN"):
        raise PopulationError("Public reference response has invalid coordinates or sequence")
    return {"contig": contig, "start": lo, "end": hi, "sequence": sequence,
            "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest()}, url


def acquire(protocol_path, output):
    try:
        import pysam
    except ImportError as exc:
        raise PopulationError("Install the genomics extra to retrieve indexed public variants") from exc
    protocol = load_protocol(protocol_path)
    if output.exists():
        raise PopulationError("Acquisition output already exists; use a new protected directory")
    output.mkdir(parents=True, mode=0o700)
    (output/"frozen-protocol.json").write_bytes(Path(protocol_path).read_bytes())
    for name in ("20220804_manifest.txt", "README_1kGP_phased_panel_110722.pdf", VCF+".tbi"):
        retrieve(BASE+name, output/name)
    retrieve(PEDIGREE, output/"pedigree.txt")
    lo, hi = protocol["reference_acquisition_interval"]
    window, url = fetch_window("chr11", lo, hi, output/"reference-api.json")
    for start in protocol["intended_starts_0based"].values():
        seq = window["sequence"][start-lo:start-lo+23]
        if seq[:20] != protocol["guide_sequence"] or seq[21:] != "GG":
            raise PopulationError("Frozen guide or PAM disagrees with retrieved reference")
    with pysam.VariantFile(BASE+VCF, index_filename=str(output/(VCF+".tbi"))) as source:
        samples = sorted(source.header.samples)[:8]
        if len(samples) != 8:
            raise PopulationError("Source has fewer than eight samples")
        source.subset_samples(samples)
        with pysam.VariantFile(str(output/"bounded-original.vcf"), "w", header=source.header) as dest:
            for n, row in enumerate(source.fetch("chr11", lo, hi), 1):
                if n > 10000:
                    raise PopulationError("Public subset exceeds the predeclared row cap")
                dest.write(row)
    (output/"selected-samples.json").write_text(json_text(samples))
    (output/"reference-acquired.json").write_text(json_text({"record_type": "reference_windows", "schema_version": "1", "genome_build": "GRCh38", "windows": [window]}))
    manifest = {"record_type": "public_acquisition", "schema_version": "1", "protocol_sha256": checksum(protocol_path),
                "source_vcf_url": BASE+VCF, "source_reference_url": url, "pysam_version": pysam.__version__,
                "reference_scope": "interval_sequence_only", "full_reference_checksum": None,
                "source_vcf_checksum_scope": "upstream manifest supplies full-file MD5; local SHA-256 binds the bounded subset",
                "files": {p.name: checksum(p) for p in output.iterdir() if p.is_file()}}
    (output/"acquisition.json").write_text(json_text(manifest))
    for path in output.iterdir():
        path.chmod(0o600)


def prepare(protocol_path, output, example_config):
    """Use bcftools on an explicitly translated reference window; retain original GT indexes."""
    import pysam
    import pysam.bcftools
    import pysam.version
    protocol = load_protocol(protocol_path)
    if checksum(protocol_path) != checksum(output/"frozen-protocol.json"):
        raise PopulationError("Protocol changed after acquisition")
    acquisition = read_json(output/"acquisition.json")
    for name, expected in acquisition["files"].items():
        if Path(name).name != name or checksum(output/name) != expected:
            raise PopulationError("Acquired source bytes changed before normalization")
    raw = output/"bounded-original.vcf"
    lines = raw.read_text().splitlines()
    records = [line.split("\t") for line in lines if not line.startswith("#")]
    header = next(line for line in lines if line.startswith("#CHROM"))
    reference_path = output/"reference-windows.json"
    reference = read_json(output/"reference-acquired.json")
    window = reference["windows"][0]
    lo, hi = window["start"], window["end"]
    original, literal, unsupported = {}, [], []
    for index, row in enumerate(records):
        pos, ref, alts = int(row[1])-1, row[3], row[4].split(",")
        oid = "source_record_" + str(index+1)
        original[oid] = {"original_record_index": index+1, "fields": row[:9], "original_alt_indexes": list(range(1, len(alts)+1))}
        if not any(w["start"] <= pos and pos+len(ref) <= w["end"] for w in reference["windows"]):
            extra, url = fetch_window("chr11", pos, pos+len(ref), output/(oid+"-reference-api.json"))
            if extra["sequence"] != ref:
                raise PopulationError("Boundary-overlapping source REF disagrees with reference")
            reference["windows"].append(extra)
        if any(set(a)-set("ACGT") for a in alts) or max([len(ref)]+list(map(len, alts))) > 50:
            unsupported.append(row)
            continue
        if not (lo < pos and pos+len(ref) < hi):
            raise PopulationError("Normalizer requires additional real reference context")
        changed = list(row)
        changed[0], changed[1], changed[2], changed[7] = "normalization_window", str(pos-lo+1), oid, "."
        literal.append(changed)
    reference_path.write_text(json_text(reference))
    fasta = output/"normalization-window.fa"
    fasta.write_text(">normalization_window\n"+window["sequence"]+"\n")
    pysam.faidx(str(fasta))
    meta = [line for line in lines if line.startswith(("##FORMAT=", "##FILTER="))]
    local = output/"normalization-local.vcf"
    local.write_text("##fileformat=VCFv4.3\n##contig=<ID=normalization_window,length="+str(hi-lo)+">\n"+"\n".join(meta)+"\n"+header+"\n"+"".join("\t".join(row)+"\n" for row in literal))
    normalized = output/"normalized-local.vcf"
    argv = ["-f", str(fasta), "-c", "e", "-Ov", "-o", str(normalized), str(local)]
    pysam.bcftools.norm(*argv, catch_stdout=False)
    normalized_rows, links = [], []
    for line in normalized.read_text().splitlines():
        if line.startswith("#"):
            continue
        row = line.split("\t"); oid = row[2]
        try:
            require_normalization_context(row, lo)
        except PopulationError as exc:
            write_file(output/"normalization-failure.json", json_text({"record_type": "normalization_failure",
                "schema_version": "1", "status": "unresolved", "reason": str(exc), "source_record": oid,
                "reference_start": lo, "source_sha256": checksum(raw)}), overwrite=True)
            raise
        source = original[oid]
        if len(row[4].split(",")) != len(source["original_alt_indexes"]):
            raise PopulationError("Normalizer altered original ALT-index mapping")
        row[0], row[1] = "chr11", str(int(row[1])+lo)
        row[2] = source["fields"][2]
        links.append({"source": source, "normalized_fields": row[:9], "alt_index_mapping": source["original_alt_indexes"]})
        normalized_rows.append(row)
    normalized_rows += unsupported
    normalized_rows.sort(key=lambda r: (int(r[1]), r[3], r[4]))
    norm_version = pysam.version.__bcftools_version__
    # Keep source INFO and FORMAT declarations for untouched unsupported rows.
    meta = [line for line in lines if line.startswith(("##INFO=", "##FORMAT=", "##FILTER=", "##contig="))]
    vcf = output/"variants.vcf"
    vcf.write_text("##fileformat=VCFv4.3\n##reference=GRCh38\n##bcftools_normVersion="+norm_version+"\n"+"\n".join(meta)+"\n"+header+"\n"+"".join("\t".join(r)+"\n" for r in normalized_rows))
    (output/"normalization-provenance.json").write_text(json_text({"normalizer": "bcftools "+norm_version,
        "argv": argv, "source_sha256": checksum(raw), "output_sha256": checksum(vcf),
        "reference_window": {"contig": "chr11", "start": lo, "end": hi, "sequence_sha256": window["sequence_sha256"]},
        "coordinate_transform": "local POS = genomic POS - window start; genomic POS restored after normalization; no ALT splitting",
        "literal_records": len(literal), "unsupported_records_retained": len(unsupported), "record_mapping": links}))
    config = scientific_config(protocol, example_config, "bcftools "+norm_version)
    samples = read_json(output/"selected-samples.json")
    pedigree_path = output/"pedigree.txt"
    families = {}
    if pedigree_path.exists():
        lines = pedigree_path.read_text().splitlines()
        if lines[0].split() != ["sampleID", "fatherID", "motherID", "sex"]:
            raise PopulationError("Public pedigree columns changed")
        parents = {}
        def find(person):
            parents.setdefault(person, person)
            while parents[person] != person:
                person = parents[person]
            return person
        for line in lines[1:]:
            sample, father, mother, sex = line.split()
            find(sample)
            for relative in (father, mother):
                if relative != "0":
                    parents[find(relative)] = find(sample)
        groups = {}
        for person in parents:
            groups.setdefault(find(person), []).append(person)
        families = {person: identity("family", sorted(group)) for group in groups.values() for person in group}
    from .variants import DONOR_COLUMNS
    from .reporting import tsv
    inventory = [{"donor_id": "public_"+str(i+1), "genotype_sample_id": sample, "source_kind": "public_genome", "family_id": families.get(sample, ""),
        "specimen_id": "scenario_specimen_"+str(i+1), "specimen_context_id": "scenario_hspc", "availability": "unverified",
        "use_permitted": "true", "activation_cost": protocol["activation_cost"], "simulation_eligible": "true", "sampling_stratum": ""} for i, sample in enumerate(samples)]
    (output/"donors.tsv").write_text(tsv(DONOR_COLUMNS, inventory))
    (output/"guides.tsv").write_text("guide_id\tspacer_sequence\teditor\teditor_class\tpam_class\tediting_window_start_1based\tediting_window_end_1based_inclusive\tintended_target_description\tlocus_system\nHBG_m198_ABE8e\t"+protocol["guide_sequence"]+"\tABE8e\tABE\tNGG\t4\t8\tExisting reference HBG -198 design; simulation only.\tHBG1_HBG2\n")
    (output/"intended-loci.tsv").write_text("guide_id\tgenome_build\tchrom\tstart_0based\tend_0based_exclusive\tstrand\tlocus_name\n"+"".join(f"HBG_m198_ABE8e\tGRCh38\tchr11\t{start}\t{start+23}\t+\t{gene}\n" for gene,start in protocol["intended_starts_0based"].items()))
    urls = {"reference": "https://api.genome.ucsc.edu/", "variants": BASE+VCF, "donors": BASE, "guides": "https://github.com/tramngo1603/base-edit-safety", "intended_loci": "https://github.com/tramngo1603/base-edit-safety"}
    config["input_sources"] = {k: {"source_id": "public_bounded_"+k, "source_class": "public_genome" if k in {"variants", "reference", "donors"} else "study_definition",
        "release": protocol["release"] if k in {"variants", "donors"} else "UCSC hg38 interval" if k == "reference" else "frozen_hbg_m198_protocol_v1",
        "url": urls[k], "use_notes": "Public genomic research computation; specimen availability and assay capability are simulation assumptions. Individual linkage is protected locally.",
        "sha256": checksum(output/v)} for k,v in config["inputs"].items() if v is not None}
    from .contracts import validate_config
    validate_config(config)
    (output/"study-config.json").write_text(json_text(config))
    for p in output.iterdir():
        if p.is_file():
            p.chmod(0o600)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--protocol", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--example-config", type=Path, required=True)
    p.add_argument("--prepare-existing", action="store_true")
    args = p.parse_args(argv)
    try:
        if not args.prepare_existing:
            acquire(args.protocol, args.output_dir)
        prepare(args.protocol, args.output_dir, args.example_config)
    except Exception as exc:
        print("Required public-data acquisition or normalization failed; no benchmark completion claimed.", file=sys.stderr)
        print(str(exc) if isinstance(exc, PopulationError) else type(exc).__name__, file=sys.stderr)
        return 3
    print("Bounded public inputs prepared. Run the population workflow using the protected study-config.json.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
