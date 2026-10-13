"""Primary-reference and sample sequence searches with pinned Cas-OFFinder."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import gzip
from datetime import datetime, timezone
from functools import lru_cache
import json
import os
from pathlib import Path
import re
import resource
import sys
import time

from .common import PopulationError, ResourceLimit, checksum, identity, save_json
from .haplotypes import local_status
from .identity import bounding_interval, canonical_path, candidate_identity, reverse_complement
from .search import evaluation_span, window_key
from .variants import Reference
from .provenance import code_identity
from .research_genotypes import iter_variant_clusters, materialize_cluster
from .research_native import run_native, parsed_rows
from .research_coverage import Coverage, load_confidence, merge
from ..pipeline.candidate_panel import CandidatePanelError

PRIMARY = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
DEFAULT_LIMITS = {"max_guides": 32, "max_variant_rows": 8000000,
                  "max_cluster_bases": 100000, "max_native_rows": 2000000,
                  "timeout_seconds": 21600, "max_resident_bytes": 10 * 1024**3,
                  "max_output_bytes": 20 * 1024**3}
CANDIDATE_COLUMNS = ["configuration_id", "locus_id", "guide_id", "contig", "start0", "end0", "strand",
                     "target_sequence", "mismatches", "reference_present", "sample_copies",
                     "phase_status", "coordinate_path", "evidence_links"]
CHANGE_COLUMNS = ["locus_id", "guide_id", "contig", "start0", "end0", "strand", "status",
                  "sequence_changed", "reference_configurations", "sample_configurations", "reason"]


def write_tsv(path, rows, columns):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(v, sort_keys=True, separators=(",", ":"))
                             if isinstance(v, (list, dict)) else v for k, v in row.items()})
    Path(path).chmod(0o600)


def _pinned_file(base, spec, role):
    if not isinstance(spec, dict) or set(spec) != {"path", "sha256"}:
        raise PopulationError(role + " requires path and sha256")
    path = (base / spec["path"]).resolve()
    if not re.fullmatch(r"[a-f0-9]{64}", spec["sha256"]):
        raise PopulationError(role + " requires a full SHA-256")
    if not path.is_file() or checksum(path) != spec["sha256"]:
        raise PopulationError(role + " checksum mismatch or missing file")
    return path


def prepare(config_path):
    """Validate a request and return its config, paths, hashes, reference, guides, contigs and limits."""
    config_path = Path(config_path).resolve()
    config = json.loads(config_path.read_text())
    required = {"format_version", "input_kind", "assembly", "reference", "guides", "backend"}
    allowed = required | {"contigs", "max_mismatches", "sample", "evidence", "limits"}
    if required - set(config) or set(config) - allowed or config["format_version"] != "1":
        raise PopulationError("Invalid research request fields or format_version")
    artificial = config["input_kind"] == "artificial"
    if config["input_kind"] not in {"artificial", "public_genome", "research_sample"}:
        raise PopulationError("Declare the actual input_kind")
    if config["assembly"] != ("artificial" if artificial else "GRCh38"):
        raise PopulationError("Real inputs require GRCh38; artificial inputs require artificial")
    contigs = config.get("contigs", PRIMARY)
    if not contigs or len(set(contigs)) != len(contigs) or any(
            not isinstance(c, str) or (not re.fullmatch(r"artificial_[A-Za-z0-9_.-]+", c) if artificial else c not in PRIMARY)
            for c in contigs):
        raise PopulationError("Unsupported or duplicate search contigs")
    mm = config.get("max_mismatches", 6)
    if type(mm) is not int or not 0 <= mm <= 6:
        raise PopulationError("max_mismatches must be an integer from 0 to 6")
    limits = dict(DEFAULT_LIMITS)
    for key, value in config.get("limits", {}).items():
        if key not in limits or type(value) is not int or value < 1:
            raise PopulationError("Invalid resource limit")
        limits[key] = value
    paths = {"config": config_path}
    for role in ("reference", "guides"):
        paths[role] = _pinned_file(config_path.parent, config[role], role)
    reference = Reference(paths["reference"], {"genome_build": config["assembly"],
                          "contig_alias_map": {c: c for c in contigs}, "permitted_contigs": contigs},
                          synthetic=artificial)
    if not reference.index and not artificial:
        raise PopulationError("Real reference must be an uncompressed indexed FASTA")
    with paths["guides"].open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["guide_id", "spacer_sequence"]:
            raise PopulationError("Guide TSV columns must be guide_id and spacer_sequence")
        guides = list(reader)
    if not 1 <= len(guides) <= limits["max_guides"]:
        raise ResourceLimit("Guide count exceeds the configured panel limit")
    if len({g["guide_id"] for g in guides}) != len(guides) or len({g["spacer_sequence"] for g in guides}) != len(guides):
        raise PopulationError("Guide IDs and spacers must be unique; retain study aliases in the evidence table")
    for guide in guides:
        if not re.fullmatch(r"[A-Za-z0-9_.-]{1,100}", guide["guide_id"]) or not re.fullmatch(r"[ACGT]{20}", guide["spacer_sequence"]):
            raise PopulationError("Require safe guide IDs and 20-base uppercase ACGT spacers")
        guide.update(pam_class="NGG", purpose="computational_only", editor_class="unassessed")
    backend = config["backend"]
    if set(backend) != {"path", "sha256", "version", "device"} or backend["version"] != "2.4.1":
        raise PopulationError("Backend must pin Cas-OFFinder 2.4.1, executable hash and device")
    paths["backend"] = _pinned_file(config_path.parent, {k: backend[k] for k in ("path", "sha256")}, "backend")
    if not os.access(paths["backend"], os.X_OK) or not re.fullmatch(r"[CG](?:[0-9]+(?:,[0-9]+)*)?", backend["device"]):
        raise PopulationError("Invalid executable or Cas-OFFinder device")
    sample = config.get("sample")
    if sample is not None:
        sample_required = {"id", "variants", "normalization", "genotype_policy", "ploidy_intervals"}
        if not sample_required <= set(sample) or set(sample) - sample_required - {"excluded_contigs", "ploidy_source", "confidence_bed", "confidence_source"}:
            raise PopulationError("Invalid sample fields")
        if not re.fullmatch(r"[A-Za-z0-9_.-]{1,100}", sample["id"]) or not sample["normalization"]:
            raise PopulationError("Sample ID and actual normalization provenance are required")
        paths["variants"] = _pinned_file(config_path.parent, sample["variants"], "variants")
        policy = sample["genotype_policy"]
        if set(policy) != {"phase", "min_gq", "min_dp", "accepted_site_filters", "accepted_sample_filters"}:
            raise PopulationError("Invalid genotype policy fields")
        if policy["phase"] not in {"explicit_ps", "vcf_implicit_phase_set"}:
            raise PopulationError("Declare the actual phase policy")
        for key in ("min_gq", "min_dp"):
            if policy[key] is not None and (type(policy[key]) is not int or policy[key] < 0):
                raise PopulationError("Invalid genotype quality threshold")
        for key in ("accepted_site_filters", "accepted_sample_filters"):
            if not policy[key] or len(set(policy[key])) != len(policy[key]) or set(policy[key]) - {"PASS", "."}:
                raise PopulationError("Declare accepted PASS and/or dot filters")
        exclusions = sample.get("excluded_contigs", [])
        if not isinstance(exclusions, list) or any(not isinstance(row, dict) or set(row) != {"contig", "reason"}
                or row["contig"] not in contigs or not isinstance(row["reason"], str) or not row["reason"].strip()
                for row in exclusions) or len({row["contig"] for row in exclusions}) != len(exclusions):
            raise PopulationError("Sample excluded_contigs requires unique searched contigs and explicit reasons")
        if "ploidy_source" in sample and (not isinstance(sample["ploidy_source"], str) or not sample["ploidy_source"].strip()):
            raise PopulationError("Ploidy source must describe the supporting sample record")
        if "confidence_bed" in sample:
            if not isinstance(sample.get("confidence_source"), str) or not sample["confidence_source"].strip():
                raise PopulationError("Confidence BED requires its source and interpretation")
            paths["confidence_bed"] = _pinned_file(config_path.parent, sample["confidence_bed"], "confidence_bed")
        intervals = sample["ploidy_intervals"]
        if not isinstance(intervals, list):
            raise PopulationError("ploidy_intervals must be an explicit list")
        for i, row in enumerate(intervals):
            if set(row) != {"contig", "start", "end", "copies"} or row["contig"] not in contigs:
                raise PopulationError("Invalid ploidy interval")
            if any(type(row[k]) is not int for k in ("start", "end", "copies")) or row["copies"] not in (0, 1, 2):
                raise PopulationError("Only explicitly declared 0, 1 or 2 copies are supported")
            if not artificial and row["contig"] in {"chrX", "chrY"} and not sample.get("ploidy_source"):
                raise PopulationError("X/Y ploidy requires an explicit source-supported ploidy_source, including PAR treatment")
            lo, hi = reference.bounds(row["contig"])[0]
            if not lo <= row["start"] < row["end"] <= hi:
                raise PopulationError("Ploidy interval is outside the reference")
            if any(other["contig"] == row["contig"] and other["start"] < row["end"] and row["start"] < other["end"] for other in intervals[:i]):
                raise PopulationError("Overlapping ploidy intervals")
    if "evidence" in config:
        paths["evidence"] = _pinned_file(config_path.parent, config["evidence"], "evidence")
    pins = {k: checksum(v) for k, v in paths.items()}
    if reference.index_sha256:
        paths["reference_index"] = Path(str(paths["reference"]) + ".fai")
        pins["reference_index"] = reference.index_sha256
    return config, paths, pins, reference, guides, contigs, limits


def _rss():
    scale = 1 if sys.platform == "darwin" else 1024
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * scale


def _reference_fasta(reference, contigs, path, check_limits=lambda: None):
    coverage = []
    path.mkdir(mode=0o700)
    for contig in contigs:
        with (path / (contig + ".fa")).open("w") as handle:
            bounds = reference.bounds(contig)
            if len(bounds) != 1:
                raise PopulationError("Require one complete interval per search contig")
            lo, hi = bounds[0]
            if lo != 0:
                raise PopulationError("Primary-reference searches require complete contigs starting at zero")
            handle.write(">" + contig + "\n")
            run, valid, unknown_bases = 0, 0, 0
            acgt_intervals = []
            for start in range(lo, hi, 1000000):
                check_limits()
                seq = reference.fetch(contig, start, min(start + 1000000, hi))
                if set(seq) - set("ACGTN"):
                    raise PopulationError("Reference contains bases outside ACGTN")
                handle.write(seq + "\n")
                for part in re.finditer(r"[ACGT]+|N+", seq):
                    value = part[0]
                    if value[0] == "N":
                        valid += max(0, run - 22)
                        run = 0
                        unknown_bases += len(value)
                    else:
                        run += len(value)
                        a, b = start+part.start(), start+part.end()
                        if acgt_intervals and acgt_intervals[-1][1] == a:
                            acgt_intervals[-1] = (acgt_intervals[-1][0], b)
                        else:
                            acgt_intervals.append((a,b))
            valid += max(0, run - 22)
            total = max(0, hi - 22)
            coverage.append({"contig": contig, "bases": hi, "unknown_bases": unknown_bases, "acgt_intervals": acgt_intervals,
                             "oriented_windows_per_guide": total * 2,
                             "evaluable_oriented_windows_per_guide": valid * 2,
                             "unknown_oriented_windows_per_guide": (total - valid) * 2})
    return coverage


def _candidate(reference, guide, contig, strand, target, path):
    cid, _ = candidate_identity(reference.sha256, contig, strand, guide, target, path)
    lo, hi = bounding_interval(path)
    forward = path if strand == "+" else path[::-1]
    key = window_key(forward)
    locus = identity("locus", [guide["spacer_sequence"], contig, strand, key])
    return {"configuration_id": cid, "locus_id": locus, "guide_id": guide["guide_id"],
            "contig": contig, "start0": lo, "end0": hi, "strand": strand,
            "target_sequence": target, "mismatches": sum(a != b for a, b in zip(target[:20], guide["spacer_sequence"])),
            "coordinate_path": canonical_path(path), "reference_present": False,
            "sample_copies": [], "phase_status": [], "evidence_links": [], "indel_path": not isinstance(key, int)}


def _sample_ploidy(sample, contig, lo, hi):
    if sample and any(row["contig"] == contig for row in sample.get("excluded_contigs", [])):
        return None
    if not sample:
        return None
    cursor, copies = lo, None
    for row in sorted(sample["ploidy_intervals"], key=lambda r: (r["contig"], r["start"])):
        if row["contig"] != contig or row["end"] <= cursor:
            continue
        if row["start"] > cursor or copies is not None and copies != row["copies"]:
            return None
        copies = row["copies"]
        cursor = row["end"]
        if cursor >= hi:
            return copies
    return None


def _reference_candidates(rows, reference, by_guide):
    configs = {}
    for row in rows:
        contig, lo, strand = row["chrom"], row["start_0based"], row["strand"]
        forward = reference.fetch(contig, lo, lo + 23)
        target = forward if strand == "+" else reverse_complement(forward)
        if target != row["site_sequence"]:
            raise PopulationError("Native candidate does not match pinned reference coordinates")
        path = [{"kind": "reference", "pos0": i} for i in range(lo, lo + 23)]
        if strand == "-":
            path.reverse()
        c = _candidate(reference, by_guide[row["guide_id"]], contig, strand, target, path)
        c["reference_present"] = True
        if c["configuration_id"] in configs:
            raise PopulationError("Duplicate reference candidate identity")
        configs[c["configuration_id"]] = c
    return configs


def _exact_copy_status(material, lo, hi):
    ok, reason = local_status(material, lo, hi)
    if not ok:
        return False, reason, None
    positions = {p["pos0"]: i for i, p in enumerate(material["sequence_path"]) if p["kind"] == "reference"}
    offset = positions.get(lo)
    if offset is None:
        return False, "indel_path_not_shared", None
    path = material["sequence_path"][offset:offset + 23]
    if len(path) != 23 or window_key(path) != lo:
        return False, "indel_path_not_shared", None
    seq = material["sequence_genomic_forward"][offset:offset + 23]
    if set(seq) - set("ACGT"):
        return False, "sequence_unknown", None
    return True, reason, seq


def _join_evidence(path, configs, sample=None, genotype_status="not_requested"):
    """Read outcome coordinates only after both candidate searches have finished."""
    if path is None:
        return []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"case_id", "guide_id", "contig", "start0", "end0", "strand", "target_sequence",
                    "stratum", "outcome", "citation", "source_locator", "eligibility", "reason"}
        if not required <= set(reader.fieldnames or []):
            raise PopulationError("Evidence table is missing required columns")
        rows = list(reader)
    lookup, locus_lookup = defaultdict(list), defaultdict(list)
    for c in configs.values():
        if c["indel_path"]:
            continue
        locus_lookup[(c["guide_id"], c["contig"], str(c["start0"]), str(c["end0"]), c["strand"])].append(c)
        lookup[(c["guide_id"], c["contig"], str(c["start0"]), str(c["end0"]), c["strand"], c["target_sequence"])].append(c)
    seen = set()
    for row in rows:
        if row["case_id"] in seen:
            raise PopulationError("Duplicate evidence case_id")
        seen.add(row["case_id"])
        key = tuple(row[k] for k in ("guide_id", "contig", "start0", "end0", "strand", "target_sequence"))
        hits = lookup.get(key, []) if row["eligibility"] == "supported" else []
        locus_hits = locus_lookup.get(key[:5], []) if row["eligibility"] == "supported" else []
        row["reference_locus_nominated"] = any(c["reference_present"] for c in locus_hits)
        row["sample_locus_nominated"] = any(c["sample_copies"] for c in locus_hits)
        row["reference_recovered"] = any(c["reference_present"] for c in hits)
        row["sample_sequence_nominated"] = any(c["sample_copies"] for c in hits)
        row["matched_sample_recovery"] = (
            "unlinked" if not sample or row.get("sample_id") != sample["id"]
            else "unresolved_sample_run" if genotype_status != "complete_conditional_model"
            else "unresolved_sample_scope" if any(e["contig"] == row["contig"] for e in sample.get("excluded_contigs", []))
            else "recovered" if row["sample_sequence_nominated"]
            else "not_nominated_or_unresolved")
        row["configuration_ids"] = [c["configuration_id"] for c in hits]
        for c in hits:
            c["evidence_links"].append({"case_id": row["case_id"], "citation": row["citation"], "source_locator": row["source_locator"]})
    return rows


def run(config_path, output, *, progress=lambda message: None):
    """Search a pinned request, write its output files to a new directory and return the run manifest."""
    started = time.monotonic()
    progress("Checking input files and checksums")
    config, paths, pins, reference, guides, contigs, limits = prepare(config_path)
    output = Path(output).resolve()
    if output.exists() or any(output == p or output in p.parents for p in paths.values()):
        raise PopulationError("Output must be new and must not contain an input")
    output.mkdir(mode=0o700, parents=True)
    os.chmod(output, 0o700)
    def check_limits():
        if time.monotonic() - started > limits["timeout_seconds"]:
            raise ResourceLimit("Workflow exceeded timeout_seconds")
        if _rss() > limits["max_resident_bytes"]:
            raise ResourceLimit("Workflow exceeded max_resident_bytes")
    run_id = identity("research", {"config": config, "inputs": pins})
    save_json(output / "request.json", config)
    state = {"status": "running", "run_id": run_id, "started_utc": datetime.now(timezone.utc).isoformat(),
             "inputs": pins, "source": code_identity(), "limits": limits, "purpose": "computational_only", "biological_benefit": "unestablished"}
    save_json(output / "run-manifest.json", state)
    try:
        by_guide = {g["guide_id"]: g for g in guides}
        ref_fasta = output / "reference-search"
        progress("Preparing reference chromosomes")
        coverage = _reference_fasta(reference, contigs, ref_fasta, check_limits)
        reference_search_bytes = sum(p.stat().st_size for p in ref_fasta.iterdir())
        progress("Searching reference chromosomes with Cas-OFFinder")
        raw, ref_execution = run_native(ref_fasta, output / "reference-native", guides, config, paths["backend"], limits, check_limits=check_limits)
        rows, unknown_rows, ref_accounting = parsed_rows(raw, guides, config, limits, check_limits=check_limits)
        configs = _reference_candidates(rows, reference, by_guide)
        exclusions = []
        for row in unknown_rows:
            forward = reference.fetch(row[1], int(row[2]), int(row[2]) + 23)
            if (forward if row[4] == "+" else reverse_complement(forward)) != row[3].upper():
                raise PopulationError("Unknown native target disagrees with reference")
        progress("Reconstructing sample variant neighborhoods")
        sample = config.get("sample")
        confidence = load_confidence(paths["confidence_bed"], reference) if "confidence_bed" in paths else None
        coverage_counter = Coverage(coverage, sample, confidence)
        dynamic_exclusions = {}
        statuses = {}
        refs_by_contig = defaultdict(list)
        for c in configs.values():
            refs_by_contig[c["contig"]].append(c)
            copies = _sample_ploidy(sample, c["contig"], c["start0"], c["end0"])
            statuses[c["locus_id"]] = {"evaluable": copies in (1, 2), "reason": "conditional_reference_fill" if copies in (1, 2) else "sample_scope_or_ploidy_unknown"}
            if copies in (1, 2):
                c["sample_copies"] = list(range(copies))
                c["phase_status"] = ["no_heterozygote_in_retained_model"]
        for values in refs_by_contig.values():
            values.sort(key=lambda c: c["start0"])
        sample_execution = None
        sample_accounting = Counter()
        variant_classes, call_statuses, phase_statuses = Counter(), Counter(), Counter()
        sample_by_contig = defaultdict(Counter)
        genotype_status = "not_requested"
        if sample:
            genotype_status = "complete_conditional_model"
            map_path, var_fasta = output / "variant-windows.jsonl.gz", output / "variant-search.fa"
            pointers = defaultdict(int)
            try:
                with map_path.open("wb") as mapping, var_fasta.open("w") as fasta:
                    for cluster in iter_variant_clusters(paths["variants"], reference, sample["id"], sample["genotype_policy"],
                            sample["ploidy_intervals"], max_rows=limits["max_variant_rows"], max_cluster_bases=limits["max_cluster_bases"],
                            check_limits=check_limits, excluded_contigs=[r["contig"] for r in sample.get("excluded_contigs", [])]):
                        check_limits()
                        interval = cluster["interval"]
                        sample_accounting["variant_clusters"] += 1
                        accounted = set(cluster.get("accounted_variant_ids", [v["variant_id"] for v in cluster["variants"]]))
                        accounted_variants = [v for v in cluster["variants"] if v["variant_id"] in accounted]
                        sample_accounting["variant_rows"] += len(accounted_variants)
                        sample_accounting["source_vcf_rows"] += len(cluster["source_rows"])
                        sample_by_contig[interval["contig"]]["variant_clusters"] += 1
                        sample_by_contig[interval["contig"]]["variant_rows"] += len(accounted_variants)
                        sample_by_contig[interval["contig"]]["source_vcf_rows"] += len(cluster["source_rows"])
                        sample_by_contig[interval["contig"]].update(cluster.get("variant_accounting", {}))
                        for gap in cluster.get("unbounded_exclusions", []):
                            dynamic_exclusions.setdefault(gap["contig"], []).append(gap)
                        for variant in accounted_variants:
                            call = variant["calls"][sample["id"]]
                            call_statuses[call["call_status"]] += 1
                            phase_statuses[call["phase_status"]] += 1
                            for alt, edit in zip(variant["original_alts"], variant["edits"][1:]):
                                kind = ("spanning_deletion_marker" if alt == "*" else "symbolic_or_nonliteral" if set(alt) - set("ACGT") else
                                        "literal_over_resource_limit" if max(len(variant["ref"]), len(alt)) > limits["max_cluster_bases"] else
                                        "complex_length_change" if edit is None else
                                        "snv" if len(alt) == len(variant["ref"]) == 1 else edit["kind"])
                                variant_classes[kind] += 1
                        if cluster.get("sample_contig_excluded"):
                            sample_accounting["excluded_contig_source_rows"] += len(cluster["source_rows"])
                            sample_by_contig[interval["contig"]]["excluded_source_rows"] += len(cluster["source_rows"])
                            continue
                        materials = materialize_cluster(reference, cluster, sample["id"], run_id)
                        coverage_counter.cluster(cluster, materials)
                        if materials:
                            sites = materials[0]["heterozygous_sites"]
                            keys = [h["phase_key"] for h in sites]
                            if len(sites) > 1 and (any(k is None for k in keys) or any(k != keys[0] for k in keys[1:])):
                                sample_accounting["neighborhoods_with_uncertain_cross_variant_phase"] += 1
                                exclusions.append({"contig": interval["contig"], "start": min(h["start"] for h in sites),
                                    "end": max(h["end"] for h in sites),
                                    "reason": "phase_scope: only windows spanning multiple unphased sites or different phase sets are excluded"})
                        if cluster.get("exclusions"):
                            exclusions.append({**interval, "reason": json.dumps(cluster["exclusions"], sort_keys=True)})
                        if mapping.tell() + fasta.tell() + reference_search_bytes + raw.stat().st_size > limits["max_output_bytes"]:
                            raise ResourceLimit("Reconstruction output exceeded max_output_bytes")
                        offset = mapping.tell()
                        block = gzip.compress((json.dumps(cluster, separators=(",", ":"), sort_keys=True) + "\n").encode(), compresslevel=1, mtime=0)
                        mapping.write(block)
                        for material in materials:
                            copy = material["local_copy"]
                            sequence = material["sequence_genomic_forward"]
                            fasta.write(f">v{offset}z{len(block)}c{copy}\n{sequence}\n")
                            sample_accounting["reconstructed_bases"] += len(sequence)
                            sample_by_contig[interval["contig"]]["reconstructed_bases"] += len(sequence)
                            sample_accounting["materializations"] += 1
                            sample_accounting["potential_oriented_variant_windows_per_guide"] += 2 * max(0, len(sequence) - 22)
                            for gap in material["excluded_intervals"]:
                                exclusions.append({"contig": interval["contig"], **gap, "local_copy": copy})
                        references = refs_by_contig[interval["contig"]]
                        pointer = pointers[interval["contig"]]
                        while pointer < len(references) and references[pointer]["end0"] <= interval["start"]:
                            pointer += 1
                        pointers[interval["contig"]] = pointer
                        for candidate_index in range(pointer, len(references)):
                            c = references[candidate_index]
                            if c["start0"] >= interval["end"]:
                                break
                            if not cluster.get("resource_excluded") and (not interval["start"] <= c["start0"] or c["end0"] > interval["end"]):
                                continue
                            c["sample_copies"], c["phase_status"] = [], []
                            checks = [_exact_copy_status(m, c["start0"], c["end0"]) for m in materials]
                            statuses[c["locus_id"]] = {"evaluable": bool(checks) and all(x[0] for x in checks),
                                "reason": "local_reconstruction_resource_limit" if cluster.get("resource_excluded") else
                                ";".join(sorted({x[1] for x in checks if not x[0]})) or ("paired_exact_path" if checks else "sample_scope_or_ploidy_unknown")}
                if var_fasta.stat().st_size:
                    var_raw, sample_execution = run_native(var_fasta, output / "sample-native", guides, config, paths["backend"], limits, check_limits=check_limits)
                    var_rows, var_unknown, accounting = parsed_rows(var_raw, guides, config, limits, check_limits=check_limits)
                    sample_accounting.update({"native_" + k: v for k, v in accounting.items() if type(v) is int})
                    with map_path.open("rb") as mapping:
                        @lru_cache(maxsize=16)
                        def material_for(offset, length):
                            mapping.seek(offset)
                            cluster = json.loads(gzip.decompress(mapping.read(length)))
                            return materialize_cluster(reference, cluster, sample["id"], run_id)
                        for row in var_rows:
                            check_limits()
                            match = re.fullmatch(r"v([0-9]+)z([0-9]+)c([01])", row["chrom"])
                            if not match:
                                raise PopulationError("Unknown native variant-window alias")
                            mats = material_for(int(match[1]), int(match[2]))
                            material = next((m for m in mats if m["local_copy"] == int(match[3])), None)
                            if material is None:
                                raise PopulationError("Native sample copy is absent from declared ploidy")
                            offset, strand = row["start_0based"], row["strand"]
                            path = material["sequence_path"][offset:offset + 23]
                            seq = material["sequence_genomic_forward"][offset:offset + 23]
                            target = seq if strand == "+" else reverse_complement(seq)
                            if len(path) != 23 or target != row["site_sequence"]:
                                raise PopulationError("Native variant sequence or coordinate disagreement")
                            lo, hi = evaluation_span(path, material)
                            ok, phase = local_status(material, lo, hi)
                            if not ok:
                                exclusions.append({"contig": material["contig"], "start": lo, "end": hi, "reason": phase})
                                sample_accounting["phase_or_footprint_excluded_native_rows"] += 1
                                continue
                            if strand == "-":
                                path.reverse()
                            candidate = _candidate(reference, by_guide[row["guide_id"]], material["contig"], strand, target, path)
                            c = configs.setdefault(candidate["configuration_id"], candidate)
                            c["sample_copies"] = sorted(set(c["sample_copies"]) | {material["local_copy"]})
                            c["phase_status"] = sorted(set(c["phase_status"]) | {"no_heterozygote_in_retained_model" if phase == "homozygous" else phase})
                        for row in var_unknown:
                            match = re.fullmatch(r"v([0-9]+)z([0-9]+)c([01])", row[1])
                            if not match:
                                raise PopulationError("Unknown variant-window alias in excluded native row")
                            mat = next(m for m in material_for(int(match[1]), int(match[2])) if m["local_copy"] == int(match[3]))
                            seq = mat["sequence_genomic_forward"][int(row[2]):int(row[2]) + 23]
                            if (seq if row[4] == "+" else reverse_complement(seq)) != row[3].upper():
                                raise PopulationError("Unknown variant target disagrees with materialized sequence")
            except PopulationError as exc:
                genotype_status = "failed_sample_scope"
                exclusions.append({"contig": "all_sample_scope", "start": "", "end": "", "reason": str(exc)})
                configs = {k: c for k, c in configs.items() if c["reference_present"]}
                for c in configs.values():
                    c["sample_copies"], c["phase_status"] = [], []
                    statuses[c["locus_id"]] = {"evaluable": False, "reason": "sample_run_failed"}
        # An unbounded event invalidates earlier conditional fills on its entire contig.
        for contig, gaps in dynamic_exclusions.items():
            exclusions.extend(gaps)
            reasons = sorted({g["reason"] for g in gaps})
            reason = "; ".join(reasons)
            coverage_counter.exclude_contig(contig, reason)
            exclusions.append({"contig": contig, "start": 0, "end": reference.bounds(contig)[0][1], "reason": reason})
            for c in configs.values():
                if c["contig"] == contig:
                    c["sample_copies"], c["phase_status"] = [], []
                    statuses[c["locus_id"]] = {"evaluable": False, "reason": reason}
        configs = {k:c for k,c in configs.items() if c["reference_present"] or c["sample_copies"]}
        coverage_rows, coverage_regions, coverage_windows = coverage_counter.finish(sample_by_contig, genotype_status == "failed_sample_scope")
        comparable = {c: [] for c in contigs}
        for row in coverage_windows:
            if row["status"] == "evaluated":
                comparable[row["contig"]].append((row["start"], row["end"]))
        from bisect import bisect_right
        comparable_ends = {c: [b for a,b in spans] for c,spans in comparable.items()}
        def comparable_start(contig, start):
            spans = comparable[contig]
            i = bisect_right(comparable_ends[contig], start)
            return i < len(spans) and spans[i][0] <= start < spans[i][1]
        progress("Comparing loci and target sequences")
        loci = defaultdict(list)
        for c in configs.values():
            loci[c["locus_id"]].append(c)
        changes = []
        for locus, members in sorted(loci.items()):
            check_limits()
            c = members[0]
            refs = [m for m in members if m["reference_present"]]
            samples = [m for m in members if m["sample_copies"]]
            reason = statuses.get(locus, {}).get("reason", "reference_not_nominated")
            if not sample:
                status = "reference_only_run"
            elif c["indel_path"]:
                status, reason = "unpaired_indel_path", "exact_reference_path_unavailable"
            elif not comparable_start(c["contig"], c["start0"]):
                status, reason = "unresolved", "outside_comparable_reference_windows"
            elif refs and samples:
                status = "retained"
            elif samples:
                ref_sequence = reference.fetch(c["contig"], c["start0"], c["end0"])
                if set(ref_sequence) - set("ACGT"):
                    status, reason = "unresolved", "reference_sequence_unknown"
                else:
                    status = "gained"
            elif statuses[locus]["evaluable"]:
                status = "lost"
            else:
                status = "unresolved"
            changes.append({"locus_id": locus, **{k: c[k] for k in ("guide_id", "contig", "start0", "end0", "strand")},
                "status": status, "sequence_changed": bool(status == "retained" and refs and samples and {m["target_sequence"] for m in refs} != {m["target_sequence"] for m in samples}),
                "reference_configurations": [m["configuration_id"] for m in refs], "sample_configurations": [m["configuration_id"] for m in samples], "reason": reason})
        check_limits()
        evidence = _join_evidence(paths.get("evidence"), configs, sample, genotype_status)
        candidates = sorted(configs.values(), key=lambda c: (c["guide_id"], c["contig"], c["start0"], c["strand"], c["target_sequence"]))
        write_tsv(output / "candidates.tsv", candidates, CANDIDATE_COLUMNS)
        write_tsv(output / "locus-changes.tsv", changes, CHANGE_COLUMNS)
        write_tsv(output / "sequence-changes.tsv",
                  [c for c in changes if c["sequence_changed"]], CHANGE_COLUMNS)
        unknown_ploidy = []
        sample_masks = [{"contig": r["contig"], "start": 0, "end": reference.bounds(r["contig"])[0][1],
                         "reason": "sample_contig_excluded: " + r["reason"]}
                        for r in sample.get("excluded_contigs", [])] if sample else []
        sample_masks.extend({"contig": c, "start": 0, "end": reference.bounds(c)[0][1],
                             "reason": "; ".join(sorted({g["reason"] for g in gaps}))}
                            for c,gaps in dynamic_exclusions.items())
        exclusions.extend(sample_masks)
        if sample:
            for contig in contigs:
                cursor = 0
                for interval in sorted((r for r in sample["ploidy_intervals"] if r["contig"] == contig), key=lambda r: r["start"]):
                    if cursor < interval["start"]:
                        unknown_ploidy.append({"contig": contig, "start": cursor, "end": interval["start"], "reason": "sample_ploidy_undeclared"})
                    if interval["copies"] == 0:
                        exclusions.append({"contig": contig, "start": interval["start"], "end": interval["end"], "reason": "zero_declared_copies"})
                    cursor = interval["end"]
                end = reference.bounds(contig)[0][1]
                if cursor < end:
                    unknown_ploidy.append({"contig": contig, "start": cursor, "end": end, "reason": "sample_ploidy_undeclared"})
            exclusions.extend(unknown_ploidy)
        exclusions = list({json.dumps(r, sort_keys=True): r for r in exclusions}.values())
        write_tsv(output / "exclusions.tsv", exclusions, ["contig", "start", "end", "reason", "variant_id", "local_copy", "boundary_evidence", "missing_information", "source_contig"])
        coverage_columns = list(coverage_rows[0]) if coverage_rows else ["contig"]
        for row in coverage_rows:
            for key in row:
                if key not in coverage_columns:
                    coverage_columns.append(key)
        write_tsv(output / "coverage.tsv", coverage_rows, coverage_columns)
        write_tsv(output / "coverage-regions.tsv", coverage_regions, ["contig", "start", "end", "status", "reasons"])
        write_tsv(output / "coverage-windows.tsv", coverage_windows, ["contig", "start", "end", "status", "reasons"])
        if evidence:
            write_tsv(output / "evidence.tsv", evidence, list(evidence[0]))
        counts = dict(Counter(c["status"] for c in changes))
        coverage_result = {"reference": [{k:v for k,v in r.items() if k != "acgt_intervals"} for r in coverage],
            "chromosomes": coverage_rows,
            "coverage_definitions": {
                "bases": "Reference-coordinate bases with a modeled state on every declared copy; full chromosome lengths remain the denominator. A resolved deletion counts as known absence, not an unknown base.",
                "comparable_windows": "All possible 23-base reference windows, both strands, per guide, before PAM/mismatch filtering; reference and every declared sample copy must have an exact evaluable coordinate path.",
                "sample_windows": "Copy-specific reconstructed windows plus conditional reference fill. Unknown alleles retain masked reference placeholders; their true sequence lengths are unknown. Nominal reference-copy counts remain explicit.",
                "regions": "coverage-regions.tsv partitions every reference base. coverage-windows.tsv partitions every reference window START; [start,end) denotes starts of 23-base windows, not their full footprints.",
                "exclusion_reasons": "Reason counts can overlap; total excluded windows use their union. Indel windows may be searchable but are outside exact reference-path comparisons.",
                "confidence": "Benchmark BED membership is reported separately; it is not evidence of a genotype, phase, or donor copy at every covered base. Confidently genotyped bases are not established by this variants-only input."},
            "sample_ploidy_source": sample.get("ploidy_source", "not_provided") if sample else None,
            "confidence_source": sample.get("confidence_source") if sample else None,
            "unbounded_sample_events": dynamic_exclusions, "guide_count": len(guides), "searched_contigs": contigs,
            "unsearched_primary_contigs": [c for c in PRIMARY if c not in contigs] if config["assembly"] == "GRCh38" else [],
            "unsupported_scope": ["alternate contigs", "mitochondrial genome", "bulged alignments", "structural rearrangements", "biological editing prediction"],
            "sample_status": genotype_status,
            "sample_accounting_complete": genotype_status == "complete_conditional_model",
            "variant_disposition_note": "Counts are records, not allele copies. A multiallelic/diploid record can be applied on one copy and unsupported or reference on another; disposition columns can overlap. Whole-contig excluded records are not counted as applied.",
            "sample_scope_exclusions": sample_masks,
            "sample_scope_basis": "Local masks use supported event bounds or accepted flanking called bases; explicit and unbounded contig masks suppress sample claims independently of declared ploidy. All input rows remain accounted for.",
            "reference_native_input_layout": "one complete contig per FASTA in one native directory search",
            "sample_ploidy_intervals": sample["ploidy_intervals"] if sample else [],
            "sample_accounting": dict(sample_accounting), "reference_native_accounting": ref_accounting,
            "variant_alt_class_counts": dict(variant_classes), "variant_call_status_counts": dict(call_statuses),
            "variant_phase_status_counts": dict(phase_statuses),
            "sample_by_contig": {c: dict(v) for c, v in sample_by_contig.items()},
            "sample_undeclared_ploidy_intervals": unknown_ploidy,
            "exclusion_summary": dict(Counter(r["reason"] for r in exclusions)),
            "supported_variant_classes": "literal substitutions and normalized simple indels within declared reconstruction limits; compatible spanning-deletion markers; original multiallelic indexes retained",
            "phase_scope_note": "phase_scope rows describe neighborhoods; only target windows that span incompatible phase are excluded, even if the native search nominates none",
            "sample_basis": "conditional retained-callset sequence model; VCF absence is not a measured reference genotype",
            "reference_callability": "not_assessed", "evidence_supported_carriage": "unknown"}
        save_json(output / "coverage.json", coverage_result)
        state.update(status="complete" if genotype_status != "failed_sample_scope" else "complete_reference_sample_failed",
            finished_utc=datetime.now(timezone.utc).isoformat(), wall_seconds=time.monotonic() - started,
            parent_peak_resident_bytes=_rss(), reference_execution=ref_execution, sample_execution=sample_execution,
            reference_loci=sum(bool(c["reference_configurations"]) for c in changes), sample_loci=sum(bool(c["sample_configurations"]) for c in changes),
            candidate_configurations=len(candidates), locus_changes=counts, changed_sequence_loci=sum(c["sequence_changed"] for c in changes),
            coverage_basis="retained_callset_model", genotype_status=genotype_status)
        check_limits()
        if code_identity() != state["source"]:
            raise PopulationError("Package source changed during execution")
        for role, path in paths.items():
            if checksum(path) != pins[role]:
                raise PopulationError("Input changed during execution: " + role)
        report = f"""# Reference and sample sequence search

Status: {state['status']}. Assembly: {config['assembly']}. Guides: {len(guides)}. Search: 20 bases, NGG, both strands, at most {config.get('max_mismatches', 6)} mismatches, no bulges.

Reference nominated loci: {state['reference_loci']}. Sample nominated paths: {state['sample_loci']}. Locus accounting: {json.dumps(counts, sort_keys=True)}. Retained loci with changed nominated sequences: {state['changed_sequence_loci']}. Indel paths without an exact reference counterpart remain unpaired.

Open `candidates.tsv` for exact sequences, coordinates, copy labels and evidence links; `locus-changes.tsv` separates locus changes from sequence changes. `coverage.tsv` reports per-chromosome denominators, evaluated windows and remaining gaps; `coverage-regions.tsv` and `coverage-windows.tsv` partition bases and reference-window starts. `coverage.json` defines the units, and `exclusions.tsv` retains event reasons. When supplied, `evidence.tsv` retains every source case and its exact sequence recovery. Published coordinates were read only after search completion and did not select search windows.

Sample status: {genotype_status}. Reference fill at absent VCF records is conditional on the retained callset. Ploidy applies only in declared intervals; missing or incompatible calls and phase blocks remain unresolved. Copy labels describe local configurations and do not establish long-range phase. Reference callability and evidence-supported donor carriage are unknown.

Linked evidence retains separate biochemical and cellular outcomes. Unmeasured candidates keep their untested status.

Full workflow wall time, parent peak resident memory, and each native job’s runtime and memory are recorded in `run-manifest.json`. Retain this directory and the checksum-pinned inputs for replay.
"""
        progress("Writing reports and checking completed output files")
        report += "\n## Computational coverage\n\nAll window counts below include both strands, per guide, before PAM and mismatch filtering. Comparable windows require exact paths on every declared copy. Variant dispositions can overlap between copies.\n\n| Chromosome | Modeled / total reference bases | Comparable / total reference windows | Evaluated sample copy windows | VCF rows | Records applied / stars resolved / unsupported |\n| --- | ---: | ---: | ---: | ---: | ---: |\n"
        for row in coverage_rows:
            dispositions = " / ".join(f"{row.get(k,0):,}" for k in ("applied_variant_records", "resolved_spanning_deletion_variant_records", "unsupported_variant_records"))
            report += f"| {row['contig']} | {row['modeled_reference_bases']:,} / {row['reference_bases']:,} | {row['comparable_oriented_windows_per_guide']:,} / {row['reference_oriented_windows_per_guide']:,} | {row['sample_evaluated_copy_oriented_windows_per_guide']:,} | {row.get('source_vcf_rows',0):,} | {dispositions} |\n"
        report += "\nProcessing a chromosome does not resolve every base or window. Benchmark-region membership, when supplied, is separate from computational coverage; confident genotypes outside explicit calls are not established. See coverage.tsv for exclusions and variant dispositions.\n"
        (output / "report.md").write_text(report)
        from ..lab_report import write_report
        write_report(output, state, coverage_result, config)
        state["outputs"] = {str(p.relative_to(output)): {"sha256": checksum(p), "bytes": p.stat().st_size}
                            for p in sorted(output.rglob("*")) if p.is_file() and p.name != "run-manifest.json"}
        state["wall_seconds"] = time.monotonic() - started
        state["parent_peak_resident_bytes"] = _rss()
        if sum(v["bytes"] for v in state["outputs"].values()) > limits["max_output_bytes"]:
            raise ResourceLimit("Output bundle exceeded max_output_bytes")
        check_limits()
        save_json(output / "run-manifest.json", state)
        for path in output.rglob("*"):
            path.chmod(0o700 if path.is_dir() else 0o600)
        return state
    except BaseException as exc:
        state.update(status="failed", error=str(exc), wall_seconds=time.monotonic() - started)
        if (output / "report.md").exists():
            (output / "report.md").rename(output / "report.unfinished.md")
        (output / "report.md").write_text("# Incomplete run\n\nThe workflow failed. See run-manifest.json and the retained native logs. These outputs do not establish search completion.\n")
        (output / "report.html").unlink(missing_ok=True)
        save_json(output / "run-manifest.json", state)
        raise


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="Checksum-bound research request JSON")
    parser.add_argument("--output", required=True, help="New private output directory")
    args = parser.parse_args(argv)
    try:
        result = run(args.config, args.output)
    except (PopulationError, CandidatePanelError, ValueError, OSError) as exc:
        print(json.dumps({"status": "failed", "reason": str(exc)}), file=sys.stderr)
        return 4 if isinstance(exc, ResourceLimit) else 2
    print(json.dumps({"status": result["status"], "reference_loci": result["reference_loci"],
                      "sample_loci": result["sample_loci"], "output": str(Path(args.output).resolve())}))
    return 0 if result["status"] == "complete" else 3


if __name__ == "__main__":
    raise SystemExit(main())
