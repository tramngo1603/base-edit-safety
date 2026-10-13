"""A bounded sequence comparison using the existing reconstruction and local search."""
from __future__ import annotations

from collections import Counter
from datetime import datetime, timezone
import hashlib
from pathlib import Path
import re

from .common import (PopulationError, ResourceLimit, check_output, checksum, fields,
                     identity, integer, json_text, jsonl_text, nonempty, read_json, write_bundle)
from .comparison import HostBudget, LocalLimits
from .haplotypes import build_materializations
from .provenance import code_identity
from .search import scan, nomination_changes
from .variants import Reference, load_vcf


LIMITS = {"max_donors": 8, "max_reference_bases": 10000, "max_variant_rows": 1000,
          "max_allele_length": 50, "max_reconstructed_bases": 250000,
          "max_configurations": 17, "max_hit_rows": 10000,
          "max_input_bytes": 1024 ** 2, "max_output_bytes": 64 * 1024 ** 2,
          "timeout_seconds": 60, "max_resident_bytes": 512 * 1024 ** 2}


def _label(value, label):
    nonempty(value, label)
    if len(value) > 100 or not re.fullmatch(r"[A-Za-z0-9_.-]+", value):
        raise PopulationError(label + " must be 1–100 letters, digits, underscores, dots or hyphens")
    return value


def validate_request(config):
    """Validate the computational request's fields, input kind, geometry and genotype policy."""
    fields(config, {"format_version", "purpose", "input_kind", "assembly", "region", "guide",
                    "samples", "reference", "variants", "genotype_policy"}, label="computational request")
    if config["format_version"] != "1" or config["purpose"] != "computational_only":
        raise PopulationError("Require format_version 1 and purpose computational_only")
    if config["input_kind"] not in ("artificial", "public_genome", "research_sample"):
        raise PopulationError("Unknown input_kind; real inputs must retain their real source kind")
    region = config["region"]
    fields(region, {"contig", "start", "end"}, label="region")
    contig = _label(region["contig"], "region contig")
    integer(region["start"], "region start")
    integer(region["end"], "region end", region["start"] + 23, region["start"] + LIMITS["max_reference_bases"])
    if config["input_kind"] == "artificial":
        if config["assembly"] != "artificial" or not contig.startswith("artificial_"):
            raise PopulationError("Artificial inputs require assembly artificial and an artificial_ contig")
    elif config["assembly"] != "GRCh38" or not re.fullmatch(r"(?:chr)?(?:[1-9]|1[0-9]|2[0-2])", contig):
        raise PopulationError("Real computational inputs require GRCh38 and one primary autosome")
    guide = config["guide"]
    fields(guide, {"id", "spacer", "pam", "max_mismatches"}, label="guide")
    _label(guide["id"], "guide id")
    if not isinstance(guide["spacer"], str) or not re.fullmatch(r"[ACGT]{20}", guide["spacer"]) or guide["pam"] != "NGG":
        raise PopulationError("Computational route supports one 20-base ACGT spacer and NGG only")
    integer(guide["max_mismatches"], "max_mismatches", 0, 6)
    samples = config["samples"]
    if not isinstance(samples, list) or not 1 <= len(samples) <= LIMITS["max_donors"]:
        raise PopulationError("Require 1–8 explicit VCF sample IDs")
    for sample in samples:
        _label(sample, "sample id")
    if len(set(samples)) != len(samples):
        raise PopulationError("Duplicate sample ID")
    for role in ("reference", "variants"):
        source = config[role]
        fields(source, {"path", "sha256"} | ({"normalization"} if role == "variants" else set()), label=role)
        nonempty(source["path"], role + " path")
        if not isinstance(source["sha256"], str) or not re.fullmatch(r"[0-9a-f]{64}", source["sha256"]):
            raise PopulationError(role + " requires a lowercase SHA-256")
    nonempty(config["variants"]["normalization"], "variant normalization provenance")
    policy = config["genotype_policy"]
    fields(policy, {"phase", "min_gq", "min_dp", "accepted_site_filters", "accepted_sample_filters"}, label="genotype_policy")
    if policy["phase"] not in ("explicit_ps", "vcf_implicit_phase_set"):
        raise PopulationError("Require explicit_ps or vcf_implicit_phase_set phase policy")
    for key in ("min_gq", "min_dp"):
        if policy[key] is not None:
            integer(policy[key], key, 0)
    for key in ("accepted_site_filters", "accepted_sample_filters"):
        values = policy[key]
        if not isinstance(values, list) or not values or any(v not in ("PASS", ".") for v in values) or len(set(values)) != len(values):
            raise PopulationError(key + " must explicitly list PASS and/or dot; dot is not PASS")
    return config


def _small_file(path):
    path = Path(path)
    if not path.is_file():
        raise PopulationError("Required local input must be a regular file")
    if path.stat().st_size > LIMITS["max_input_bytes"]:
        raise ResourceLimit("max_input_bytes exceeded before ingestion")
    return path


def _vcf_envelope(path, config):
    """Bound all samples and rows, including those outside the requested interval."""
    headers, builds, row_count = [], [], 0
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#CHROM\t"):
                headers.append(line.rstrip("\r\n").split("\t")[9:])
            elif line.startswith("##reference="):
                builds.append(line.strip().split("=", 1)[1])
            elif line.strip() and not line.startswith("#"):
                row_count += 1
                if row_count > LIMITS["max_variant_rows"]:
                    raise ResourceLimit("max_variant_rows exceeded before ingestion")
    if len(headers) != 1 or len(headers[0]) != len(config["samples"]) or set(headers[0]) != set(config["samples"]):
        raise PopulationError("VCF header must contain exactly the declared sample IDs once")
    if builds != [config["assembly"]]:
        raise PopulationError("VCF requires exactly one reference declaration matching the configured assembly")


def prepare(config_path):
    path = _small_file(Path(config_path).resolve())
    config = validate_request(read_json(path))
    paths, inputs = {"config": path}, [{"role": "config", "sha256": checksum(path), "size_bytes": path.stat().st_size}]
    for role, suffix in (("reference", ".json"), ("variants", ".vcf")):
        source = config[role]
        p = _small_file((path.parent / source["path"]).resolve())
        if p.suffix != suffix:
            raise PopulationError("Require one reference-window JSON and one uncompressed .vcf")
        sha = checksum(p)
        if sha != source["sha256"]:
            raise PopulationError(role + " checksum differs from the declared input")
        paths[role] = p
        inputs.append({"role": role, "sha256": sha, "size_bytes": p.stat().st_size})
    data = read_json(paths["reference"])
    if not isinstance(data, dict) or not isinstance(data.get("windows"), list) or len(data["windows"]) != 1:
        raise PopulationError("Require exactly one reference window")
    window = data["windows"][0]
    if not isinstance(window, dict) or {key: window.get(key) for key in ("contig", "start", "end")} != config["region"]:
        raise PopulationError("Reference window must equal the declared region")
    policy = {"genome_build": config["assembly"], "contig_alias_map": {config["region"]["contig"]: config["region"]["contig"]},
              "permitted_contigs": [config["region"]["contig"]]}
    reference = Reference(paths["reference"], policy, synthetic=config["input_kind"] == "artificial")
    _vcf_envelope(paths["variants"], config)
    # Apply the same sequence model to every source kind, including mock inputs.
    # No absent VCF row or reference base supplies evidence of a person's genotype.
    legacy = {"scope": {"requested_intervals": [config["region"]], "materialized_intervals": [config["region"]]},
              "coverage_basis": "retained_callset_model", "resource_limits": dict(LIMITS),
              "reference_policy": policy,
              "input_sources": {"variants": {"source_id": "vcf:" + config["variants"]["sha256"]}},
              "genotype_policy": {**config["genotype_policy"], "reference_evidence": "reference_imputed_scenario"},
              "variant_policy": {"max_allele_length": LIMITS["max_allele_length"], "normalizer": config["variants"]["normalization"]},
              "search_profile": {"max_mismatches": config["guide"]["max_mismatches"], "pam_patterns": ["NGG"],
                                 "strands": ["+", "-"], "max_dna_bulge": 0, "max_rna_bulge": 0}}
    donors = [{"donor_id": sample, "genotype_sample_id": sample, "source_kind": config["input_kind"]}
              for sample in sorted(config["samples"])]
    guide = {"guide_id": config["guide"]["id"], "spacer_sequence": config["guide"]["spacer"],
             "pam_class": "NGG", "purpose": "computational_only"}
    # Source paths and output location do not change scientific identities.
    semantic_config = {key: value for key, value in config.items() if key not in {"reference", "variants"}}
    semantic_config["samples"] = sorted(config["samples"])
    semantic_config["reference_sha256"] = config["reference"]["sha256"]
    semantic_config["variants_sha256"] = config["variants"]["sha256"]
    semantic_config["normalization"] = config["variants"]["normalization"]
    run_id = identity("computational", semantic_config)
    return config, paths, inputs, reference, legacy, donors, guide, run_id


def _warnings(config, genotype, audit, ingest):
    rows = []

    def add(code, detail, count=1, sample=None):
        rows.append({"sample_id": sample, "code": code, "detail": detail, "count": count})

    add("biological_abstention", "No assay measurements or detection limits were supplied.")
    add("reference_callability_unavailable", "No callability or confident-reference evidence was supplied. All sample sequences are conditional retained-callset models.")
    add("donor_linkage_unavailable", "VCF sample IDs identify sequence inputs only. Physical donor, specimen and assay linkage are not supplied.")
    add("bounded_scope", "Only full 23-base windows within the declared region were searched, on both strands; no bulges or outside-region coverage.")
    if config["genotype_policy"]["min_dp"] is None or config["genotype_policy"]["min_gq"] is None:
        add("quality_threshold_unset", "At least one DP/GQ threshold is unset; this does not establish high-quality genotype calls.")
    if "." in config["genotype_policy"]["accepted_site_filters"]:
        add("site_filter_dot_accepted", "Dot is accepted as declared and retained separately; it does not mean PASS.")
    if config["genotype_policy"]["phase"] == "vcf_implicit_phase_set":
        add("implicit_phase_assumption", "Phased calls without PS share the declared VCF implicit phase set; cross-block phase was not measured here.")
    counts = Counter()
    for row in genotype:
        if row["call_status"] != "complete":
            counts[(row["donor_id"], row["call_status"])] += 1
        if row["phase_status"] == "unresolved_unphased":
            counts[(row["donor_id"], "unresolved_unphased")] += 1
    for (sample, code), count in sorted(counts.items()):
        add(code, "See genotype-audit.jsonl. One heterozygote resolves an unordered pair; nearby heterozygotes require shared phase.", count, sample)
    excluded = Counter()
    for execution in audit["executions"]:
        for reason, count in execution["counts"].items():
            if reason not in {"hits", "evaluated_windows"} and count:
                excluded[(execution["donor_id"], reason)] += count
    for (sample, code), count in sorted(excluded.items(), key=lambda item: (item[0][0] or "", item[0][1])):
        add(code, "Oriented windows excluded from search; see search-audit.json and sequence-configurations.jsonl.", count, sample)
    unpaired = Counter()
    for row in audit["paired_window_scope"]["rows"]:
        for reason, count in row["unpaired_reasons"].items():
            unpaired[(row["donor_id"], reason)] += count
    for (sample, reason), count in sorted(unpaired.items()):
        add("comparison_" + reason, "Oriented paths lack an evaluable exact reference/copy comparison; this is not evidence of absence.", count, sample)
    if ingest["outside_scope_rows"]:
        add("outside_scope_variant_rows", "VCF rows outside this region were not reconstructed or REF-validated.", ingest["outside_scope_rows"])
    return sorted(rows, key=lambda row: (row["sample_id"] or "", row["code"], row["detail"]))


def compute(config_path, *, budget=None):
    budget = budget or HostBudget(LocalLimits(LIMITS["timeout_seconds"], LIMITS["max_resident_bytes"]))
    budget.check("input_validation")
    config, paths, inputs, reference, legacy, donors, guide, run_id = prepare(config_path)
    variants, genotype, ingest = load_vcf(paths["variants"], legacy, reference, donors, run_id)
    if any(row["call_status"] == "invalid" for row in genotype):
        raise PopulationError("Invalid VCF genotype syntax or allele index")
    budget.check("reconstruction")
    materials, material_audit = build_materializations(legacy, reference, variants, donors, run_id)
    budget.check("local_search")
    candidates, hits, carriage, audit = scan(legacy, reference, materials, [guide], [], donors, [], run_id, check_limits=budget.check)
    budget.check("reference_comparison")
    changes = nomination_changes(legacy, candidates, materials, donors, [guide], carriage, run_id, check_limits=budget.check)
    by_candidate = {row["configuration_id"]: row for row in candidates}
    by_carriage = {(row["configuration_id"], row["donor_id"]): row for row in carriage}
    hit_lookup = {}
    for hit in hits:
        hit_lookup.setdefault((hit["configuration_id"], hit["donor_id"]), []).append(hit)
    rows = []
    for change in changes:
        cid, sample = change["configuration_id"], change["donor_id"]
        c, carry = by_candidate[cid], by_carriage[(cid, sample)]
        sample_hits = hit_lookup.get((cid, sample), [])
        rows.append({"candidate_id": cid, "sample_id": sample, "contig": c["contig"],
                     "start0": c["bounding_interval"][0], "end0": c["bounding_interval"][1], "strand": c["strand"],
                     "target_sequence": c["protospacer_guide_orientation"], "pam": c["pam_guide_orientation"],
                     "mismatches": sum(a != b for a, b in zip(guide["spacer_sequence"], c["protospacer_guide_orientation"])),
                     "reference_target_sequence": change["reference_counterpart"], "change": change["change"],
                     "comparison_basis": "retained_callset_model", "model_carriage": carry["model_carriage"],
                     "evidence_supported_carriage": carry["evidence_supported_carriage"],
                     "local_copies": sorted({hit["local_copy"] for hit in sample_hits}),
                     "phase_status": ["no_heterozygote_in_retained_model" if value == "homozygous" else value
                                      for value in carry["phase_status"]], "pairing_reason": change["pairing_reason"],
                     "target_path": c["target_path"], "variant_ids": c["variant_links"]})
    rows.sort(key=lambda r: (r["contig"], r["start0"], r["end0"], r["strand"], r["candidate_id"], r["sample_id"]))
    result = {"format_version": "1", "purpose": "computational_only", "input_kind": config["input_kind"], "run_id": run_id,
              "scope": {"assembly": config["assembly"], "region": config["region"], "guide": config["guide"],
                        "sample_ids": sorted(config["samples"]), "genotype_policy": config["genotype_policy"],
                        "coordinate_convention": "0-based half-open; VCF POS 1-based"},
              "domain": {"computational_status": "bounded_complete", "biological_status": "abstain",
                         "reasons": ["no_biological_validation", "no_assay_evidence_or_detection_limits", "no_donor_or_specimen_linkage"]},
              "limits": dict(LIMITS), "inputs": inputs,
              "summary": {"reference_candidates": len({h["configuration_id"] for h in hits if h["donor_id"] is None}),
                          "modeled_candidates": len({h["configuration_id"] for h in hits if h["donor_id"] is not None}),
                          "comparison_rows": len(rows), "variant_rows": len(variants),
                          "unresolved_call_rows": sum(r["call_status"] != "complete" for r in genotype),
                          "evaluated_oriented_windows": sum(e["evaluated_windows"] for e in audit["executions"]),
                          "excluded_oriented_windows": sum(n for e in audit["executions"] for key, n in e["counts"].items() if key not in {"hits", "evaluated_windows"}),
                          "comparison_basis": "retained_callset_model"},
              "candidates": rows, "warnings": _warnings(config, genotype, audit, ingest)}
    artifacts = {"result.json": json_text(result), "search-audit.json": json_text(audit),
                 "population-audit.json": json_text({"ingestion": ingest, "materialization": material_audit})}
    for name, records in (("normalized-variants", variants), ("genotype-audit", genotype),
                          ("sequence-configurations", materials), ("candidate-configurations", candidates),
                          ("search-hits", hits), ("sample-carriage", carriage), ("nomination-changes", changes)):
        artifacts[name + ".jsonl"] = jsonl_text(records)
    budget.check("report_preparation")
    return result, artifacts, paths, budget


def run(config_path, output):
    if Path(output).exists():
        raise PopulationError("Output exists; choose a new directory")
    check_output(output, [config_path])
    result, artifacts, paths, budget = compute(config_path)
    check_output(output, paths.values())
    from .computational_reporting import render
    artifacts.update(render(result))
    # Verify sources again before committing a complete, immutable output bundle.
    for source in result["inputs"]:
        if checksum(paths[source["role"]]) != source["sha256"]:
            raise PopulationError("Input changed during computation; no complete result retained")
    manifest = {"format_version": "1", "status": "completed", "run_id": result["run_id"],
                "purpose": "computational_only", "biological_status": "abstain", "code": code_identity(),
                "completed_at_utc": datetime.now(timezone.utc).isoformat(), "resource_usage": budget.metrics(),
                "inputs": result["inputs"],
                "outputs": {name: {"sha256": hashlib.sha256(content.encode("utf-8")).hexdigest(),
                                   "size_bytes": len(content.encode("utf-8"))} for name, content in sorted(artifacts.items())}}
    artifacts["run-manifest.json"] = json_text(manifest)
    if sum(len(content.encode("utf-8")) for content in artifacts.values()) > LIMITS["max_output_bytes"]:
        raise ResourceLimit("max_output_bytes exceeded; no complete result retained")
    budget.check("output_commit")
    write_bundle(output, artifacts, inputs=paths.values())
    return manifest
