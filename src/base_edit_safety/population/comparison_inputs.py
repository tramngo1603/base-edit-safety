"""Resolve pinned genomic inputs and stream exact, locally supported donor windows."""
from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field
from pathlib import Path

from ..pipeline.candidate_panel import load_guides, load_intended_loci
from .common import PopulationError, ResourceLimit, checksum, read_json, unknown
from .comparison_contracts import make_record, validate_catalog
from .comparison_identity import canonical_bytes, content_id
from .comparison_provenance import verify_sources
from .contracts import validate_config as validate_population_config
from .haplotypes import local_status, reconstruct
from .identity import reverse_complement
from .search import supported
from .variants import Reference, allele_edit, callable_intervals, canonical_donors, donors_from_tsv, load_vcf


def effective_profiles(study):
    """Return the exact active design/editor/profile triples in stable order."""
    config = study["config"]
    catalog = {row["record_id"]: row for row in config["catalogs"]}
    return [{"design": catalog[profile["design_id"]],
             "editor": catalog[catalog[profile["design_id"]]["editor_definition_id"]],
             "profile": profile}
            for profile in sorted((catalog[identifier] for identifier in config["policies"]["search"]["profile_ids"]),
                                  key=lambda row: row["design_id"])]


def path_point(point):
    fields = ("kind", "pos0", "base") if point["kind"] == "reference" else ("kind", "boundary0", "offset0", "base")
    return {name: point[name] for name in fields}


def oriented_path(path, strand):
    if strand == "+":
        return [path_point(point) for point in path]
    return [{**path_point(point), "base": reverse_complement(point["base"])} for point in reversed(path)]


def path_bounds(path):
    reference = [point["pos0"] for point in path if point["kind"] == "reference"]
    boundaries = [point["boundary0"] for point in path if point["kind"] == "insertion"]
    if not reference:
        if not boundaries:
            raise PopulationError("An empty sequence has no coordinate path")
        # Insertion-only sequences retain their reference anchor in the locus span.
        return max(0, min(boundaries) - 1), max(boundaries) + 1
    return min(reference + boundaries), max(reference + [boundary - 1 for boundary in boundaries]) + 1


def support_bounds(path, variants):
    lo, hi = path_bounds(path)
    boundaries = {point["boundary0"] for point in path if point["kind"] == "insertion"}
    for variant in variants if boundaries else []:
        edit = allele_edit(variant["pos0"], variant["reference"], variant["alternates"][variant["selected_alt_index"] - 1], 50)
        if edit is not None and edit["kind"] == "insertion" and edit["start"] in boundaries:
            lo, hi = min(lo, variant["pos0"]), max(hi, variant["pos0"] + len(variant["reference"]))
    return lo, hi


def overlaps(a, b):
    return a["contig"] == b["contig"] and a["start"] < b["end"] and b["start"] < a["end"]


def contains(outer, inner):
    return outer["contig"] == inner["contig"] and outer["start"] <= inner["start"] and inner["end"] <= outer["end"]


def _subtract(interval, allowed):
    cursor, missing = interval["start"], []
    for start, end in sorted(allowed):
        if cursor < start:
            missing.append({"start": cursor, "end": min(start, interval["end"]), "reason": "reference_callability_unavailable"})
        cursor = max(cursor, end)
    if cursor < interval["end"]:
        missing.append({"start": cursor, "end": interval["end"], "reason": "reference_callability_unavailable"})
    return [row for row in missing if row["start"] < row["end"]]


def _merge_ranges(ranges):
    merged = []
    for start, end in sorted(ranges):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def _phase_segments(material):
    """Find maximal genomic stretches that never join incompatible phase blocks."""
    interval = material["interval"]
    sites = sorted(material["heterozygous_sites"], key=lambda row: (row["start"], row["end"], row["variant_id"]))
    if not sites:
        return [(interval["start"], interval["end"])], []
    keys = [canonical_bytes(site["phase_key"]) if site["phase_key"] is not None else site["variant_id"].encode() for site in sites]
    segments, conflicts = set(), []
    for index, site in enumerate(sites):
        left = max((other["end"] for j, other in enumerate(sites[:index]) if keys[j] != keys[index]), default=interval["start"])
        right = min((other["start"] for j, other in enumerate(sites[index + 1:], index + 1) if keys[j] != keys[index]), default=interval["end"])
        if left < right:
            segments.add((left, right))
        for j in range(index + 1, len(sites)):
            if keys[index] != keys[j]:
                conflicts.append({"start": site["start"], "end": sites[j]["end"], "reason": "unresolved_phase_between_heterozygous_sites"})
    return sorted(segments), conflicts


def _split_known_path(material):
    path = material["sequence_path"]
    result, start = [], None
    for index in range(len(path) + 1):
        known_base = index < len(path) and path[index]["base"] in "ACGT"
        if known_base and start is None:
            start = index
        if not known_base and start is not None:
            result.append(path[start:index])
            start = None
    return result


@dataclass
class PreparedInputs:
    study: dict
    legacy_config: dict
    reference: Reference
    variants: list
    donors: list
    donor_records: dict
    source_records: list
    variant_source_id: str
    callable_rows: list
    legacy_guides: list
    legacy_intended: list
    reference_bundle_id: str
    material_audit: list = field(default_factory=list)
    source_equivalence_audit: dict = field(default_factory=dict)
    _variant_cache: dict = field(default_factory=dict)

    def __post_init__(self):
        self.catalog = validate_catalog(self.study["config"]["catalogs"])
        self.profiles = effective_profiles(self.study)
        self.domain = self.catalog[self.study["config"]["policies"]["search"]["domain_id"]]
        self._by_contig = {}
        self._unbounded = []
        for variant in self.variants:
            self._by_contig.setdefault(variant["contig"], []).append(variant)
            if variant["unbounded_footprint"]:
                self._unbounded.append(variant)
        self._legacy_by_canonical = {row["record_id"]: label for label, row in self.donor_records.items()}
        self._samples_by_label = {row["donor_id"]: row["genotype_sample_id"] for row in self.donors}

    def sample_id(self, donor_id):
        return self._samples_by_label[self._legacy_by_canonical[donor_id]]

    def variants_for(self, interval):
        key = interval["contig"], interval["start"], interval["end"]
        if key not in self._variant_cache:
            rows = [v for v in self._by_contig.get(interval["contig"], [])
                    if v["pos0"] < interval["end"] and interval["start"] < v["end"]]
            rows.extend(v for v in self._unbounded if v not in rows)
            self._variant_cache[key] = sorted(rows, key=lambda v: (v["contig"], v["pos0"], v["variant_id"]))
        return self._variant_cache[key]

    def reconstruct(self, interval, donor_id=None, copy_index=0):
        label = self._legacy_by_canonical[donor_id] if donor_id is not None else None
        return reconstruct(self.reference, interval, self.variants_for(interval), label, copy_index,
                           self.study["snapshot_run_id"])

    def sequence_support(self, material, interval):
        valid, reason = local_status(material, interval["start"], interval["end"])
        unknown_base = any(point["base"] == "N" and interval["start"] <= point.get("pos0", point.get("boundary0")) < interval["end"] for point in material["sequence_path"])
        if not valid or unknown_base:
            return "unresolved", reason if not valid else "sequence_unknown"
        if material["donor_id"] is None:
            return "supported", reason
        basis = self.legacy_config["coverage_basis"]
        if basis == "retained_callset_model":
            return "modeled", "retained_callset_reference_fill"
        if not supported(self.legacy_config, material["donor_id"], interval["contig"], interval["start"], interval["end"], self.callable_rows):
            return "unresolved", "reference_callability_unavailable"
        return "supported", reason

    def variant_metadata(self, interval, donor_id, copy_index):
        if donor_id is None:
            return [], [], []
        label = self._legacy_by_canonical[donor_id]
        sample = self._samples_by_label[label]
        variants, genotypes, refs = [], [], []
        for variant in self.variants_for(interval):
            call = variant["calls"][label]
            genotypes.append({"variant_id": variant["variant_id"], "contig": variant["contig"],
                "pos0": variant["pos0"], "reference": variant["ref"],
                "original_alleles": [variant["ref"], *variant["original_alts"]],
                "selected_indexes": call["selected_indexes"], "gt": call["gt"],
                "phase_key": call["phase_key"], "phase_status": call["phase_status"],
                "call_status": call["call_status"], "source_record_index": variant["source_record_index"]})
            refs.append({"source_artifact_id": self.variant_source_id,
                         "locator": {"kind": "table", "sheet": "VCF", "table": "variant_rows_excluding_headers",
                                     "row": variant["source_record_index"], "column": sample,
                                     "cell": f"variant_{variant['source_record_index']}:{sample}"},
                         "note": "Exact retained VCF row and canonical donor mapping"})
            selected = call["selected_indexes"]
            if call["call_status"] == "complete" and selected is not None and selected[copy_index] > 0:
                index = selected[copy_index]
                # Canonical literal variants are kept separate from unsupported footprints.
                if index <= len(variant["original_alts"]) and all(set(alt) <= set("ACGTN") for alt in variant["original_alts"]):
                    variants.append({"contig": variant["contig"], "pos0": variant["pos0"], "reference": variant["ref"],
                        "alternates": variant["original_alts"], "selected_alt_index": index,
                        "source_artifact_id": self.variant_source_id, "source_row": variant["source_record_index"]})
        return variants, genotypes, refs

    def _window(self, material, path, donor_id, copy_index, interval):
        sequence = "".join(point["base"] for point in path)
        lo, hi = path_bounds(path)
        span = {"contig": material["contig"], "start": lo, "end": hi}
        source_variants, _, _ = self.variant_metadata(material["interval"], donor_id, copy_index)
        lo, hi = support_bounds(path, source_variants)
        span = {"contig": material["contig"], "start": lo, "end": hi}
        variants, genotypes, refs = self.variant_metadata(span, donor_id, copy_index)
        sites = [site for site in material["heterozygous_sites"] if site["start"] < hi and lo < site["end"]]
        phase_keys = [site["phase_key"] for site in sites]
        phase = "homozygous" if not sites else "confirmed" if all(key is not None for key in phase_keys) else "unphased"
        phase_block = None if not sites or any(key is None for key in phase_keys) else canonical_bytes(phase_keys[0]).decode("utf-8")
        status, _ = self.sequence_support(material, span)
        if status == "unresolved":
            return None
        payload = {"source_material_id": material["sequence_configuration_id"],
            "reference_bundle_id": self.reference_bundle_id, "contig": material["contig"],
            "interval": span, "sequence": sequence, "path": [path_point(point) for point in path],
            "donor_id": donor_id, "sample_id": self.sample_id(donor_id) if donor_id else None,
            "copy_index": copy_index, "phase_status": phase, "phase_block": phase_block,
            "coverage_basis": self.legacy_config["coverage_basis"], "support_status": status,
            "source_version": next(source["sha256"] for source in self.study["sources"].values() if source["role"] == "variants"),
            "source_variants": variants, "source_genotypes": genotypes, "source_refs": refs}
        payload["window_id"] = "w" + content_id("comparison_prepared_window", payload).rsplit(":", 1)[1]
        validate_prepared_window(payload)
        return payload

    def iter_windows(self):
        """Stream maximal exact windows; missing sequence and phase never become reference."""
        self.material_audit = []
        limits = self.legacy_config["resource_limits"]
        requested = self.domain["intervals"]
        intervals = [interval for interval in self.legacy_config["scope"]["materialized_intervals"]
                     if any(overlaps(interval, domain) for domain in requested)]
        maximum_query = max(item["editor"]["guide_geometry"]["spacer_length"] + item["editor"]["guide_geometry"]["pam_length"] + item["profile"]["max_dna_bulges"] for item in self.profiles)
        minimum_query = min(max(1, item["editor"]["guide_geometry"]["spacer_length"] - item["profile"]["max_rna_bulges"]) + item["editor"]["guide_geometry"]["pam_length"] for item in self.profiles)
        bases, count = 0, 0
        for interval in intervals:
            for donor_id in [None, *sorted(self._legacy_by_canonical)]:
                for copy_index in range(1 if donor_id is None else 2):
                    material = self.reconstruct(interval, donor_id, copy_index)
                    bases += len(material["sequence_genomic_forward"])
                    count += 1
                    if bases > limits["max_reconstructed_bases"] or count > limits["max_configurations"]:
                        raise ResourceLimit("Comparative reconstruction exceeded its declared source limits")
                    segments, conflicts = _phase_segments(material)
                    skipped = [{"start": gap["start"], "end": gap["end"], "reason": gap["reason"]} for gap in material["excluded_intervals"]]
                    for conflict in conflicts:
                        path_length = sum(1 for point in material["sequence_path"]
                            if conflict["start"] <= point.get("pos0", point.get("boundary0")) < conflict["end"])
                        if path_length <= maximum_query:
                            skipped.append(conflict)
                    if donor_id is not None and self.legacy_config["coverage_basis"] == "callable_supported":
                        label = self._legacy_by_canonical[donor_id]
                        allowed = _merge_ranges((max(interval["start"], row["start_0based"]), min(interval["end"], row["end_0based_exclusive"])) for row in self.callable_rows if row["donor_id"] == label and row["contig"] == interval["contig"] and row["start_0based"] < interval["end"] and interval["start"] < row["end_0based_exclusive"])
                        skipped.extend(_subtract(interval, allowed))
                        segments = sorted({(max(a, lo), min(b, hi)) for a, b in segments for lo, hi in allowed if max(a, lo) < min(b, hi)})
                    emitted, evaluated, window_intervals = set(), [], {}
                    for lo, hi in segments:
                        segment_interval = {"contig": interval["contig"], "start": lo, "end": hi}
                        segment = material if segment_interval == interval else self.reconstruct(segment_interval, donor_id, copy_index)
                        for path in _split_known_path(segment):
                            if len(path) < minimum_query:
                                continue
                            window = self._window(segment, path, donor_id, copy_index, segment_interval)
                            if window is None:
                                skipped.append({"start": lo, "end": hi, "reason": "unresolved_source_anchor_or_phase"})
                                continue
                            if window["window_id"] not in emitted:
                                emitted.add(window["window_id"])
                                evaluated.append(window["interval"])
                                window_intervals[window["window_id"]] = window["interval"]
                                yield window
                    if "N" in material["sequence_genomic_forward"] and not material["excluded_intervals"]:
                        skipped.append({"start": interval["start"], "end": interval["end"], "reason": "reference_sequence_unknown"})
                    self.material_audit.append({"donor_id": donor_id, "copy_index": copy_index,
                        "interval": interval, "evaluated_intervals": evaluated,
                        "skipped_intervals": sorted(skipped, key=lambda row: (row["start"], row["end"], row["reason"])),
                        "window_count": len(emitted), "source_material_id": material["sequence_configuration_id"]})
                    self.material_audit[-1]["window_intervals"] = window_intervals
        self.source_equivalence_audit = {
            "status": "verified_reconstructed_windows", "verified": True,
            "source_panel_ids": sorted({row["source_panel_id"] for row in self.donor_records.values()}),
            "donor_sample_membership": [{"donor_id": row["record_id"], "sample_ids": row["sample_ids"]} for row in sorted(self.donor_records.values(), key=lambda row: row["record_id"])],
            "input_sha256": {key: source["sha256"] for key, source in sorted(self.study["sources"].items())},
            "genotype_policy": deepcopy(self.legacy_config["genotype_policy"]),
            "variant_policy": deepcopy(self.legacy_config["variant_policy"]),
            "ordered_alt_rule": "original VCF ALT indexes and GT copy order retained",
            "phase_rule": "maximal windows never join incompatible phase blocks",
            "missingness_rule": "missing or rejected calls remain excluded; never reference-filled",
            "coverage_basis": self.legacy_config["coverage_basis"],
            "reconstructed_bases": bases, "materializations": count,
            "material_audit": deepcopy(self.material_audit)}
        self.source_equivalence_audit["auxiliary_reference_checks"] = getattr(self, "reference_validation_audit", [])


def validate_prepared_window(window):
    required = {"window_id", "source_material_id", "reference_bundle_id", "contig", "interval", "sequence", "path", "donor_id", "sample_id", "copy_index", "phase_status", "phase_block", "coverage_basis", "support_status", "source_version", "source_variants", "source_genotypes", "source_refs"}
    if not isinstance(window, dict) or set(window) != required:
        raise PopulationError("Prepared window fields do not match the native sequence contract")
    sequence = window["sequence"]
    if not isinstance(sequence, str) or not sequence or set(sequence) - set("ACGT"):
        raise PopulationError("Prepared native windows require fully known ACGT sequence")
    if "".join(point["base"] for point in window["path"]) != sequence:
        raise PopulationError("Prepared window sequence/path mismatch")
    if window["support_status"] not in {"supported", "modeled"} or window["phase_status"] not in {"confirmed", "unphased", "homozygous"}:
        raise PopulationError("Native input cannot claim unresolved phase/sequence is complete")
    payload = {key: value for key, value in window.items() if key != "window_id"}
    expected = "w" + content_id("comparison_prepared_window", payload).rsplit(":", 1)[1]
    if window["window_id"] != expected:
        raise PopulationError("Prepared window no longer matches its immutable content identity")
    return window


class _ValidationReference:
    """Permit separately pinned REF checks without adding searchable sequence."""

    def __init__(self, reference, auxiliary):
        self.reference, self.auxiliary = reference, auxiliary
        self.build, self.sha256 = reference.build, reference.sha256
        self.checked_intervals = []

    def fetch(self, contig, start, end):
        try:
            return self.reference.fetch(contig, start, end)
        except PopulationError as original:
            if start < 0 or end <= start:
                raise
            pieces = []
            primary_bounds = self.reference.bounds(contig) if contig in self.reference.windows or contig in self.reference.index else []
            for lo, hi in primary_bounds:
                lo, hi = max(start, lo), min(end, hi)
                if lo < hi:
                    pieces.append((lo, hi, self.reference.fetch(contig, lo, hi), None))
            for source, reference in self.auxiliary:
                if contig not in reference.windows and contig not in reference.index:
                    continue
                for lo, hi in reference.bounds(contig):
                    lo, hi = max(start, lo), min(end, hi)
                    if lo < hi:
                        pieces.append((lo, hi, reference.fetch(contig, lo, hi), source["sha256"]))
            cursor, sequence, used = start, "", []
            for lo, hi, bases, source_sha256 in sorted(pieces, key=lambda part: (part[0], part[1])):
                if lo > cursor:
                    raise original
                overlap = min(cursor, hi) - lo
                if overlap and sequence[lo - start:lo - start + overlap] != bases[:overlap]:
                    raise PopulationError("Pinned REF validation sources disagree on overlapping bases")
                if hi > cursor:
                    if source_sha256 is not None:
                        used.append({"contig": contig, "start": cursor, "end": hi,
                            "source_sha256": source_sha256, "purpose": "VCF REF validation only"})
                    sequence += bases[cursor - lo:]
                    cursor = hi
            if cursor != end:
                raise original
            self.checked_intervals.extend(used)
            return sequence


def _validate_regional_bundle(study, path, bundle):
    actual = {(row["contig"], row["start"], row["end"]): row for row in read_json(path)["windows"]}
    declared = {(row["interval"]["contig"], row["interval"]["start"], row["interval"]["end"]): row for row in bundle["windows"]}
    if bundle["identity_kind"] != "indexed_windows" or actual.keys() != declared.keys():
        raise PopulationError("Regional reference intervals differ from their canonical bundle")
    sources = {row["sha256"]: row for row in study["sources"].values()}
    for coordinates, row in declared.items():
        source = sources.get(row["source_sha256"])
        if source is None or source["role"] != "reference":
            raise PopulationError("Every canonical reference-window source must be explicitly pinned")
        payload = read_json(source["path"]) if source["format"] == "json" else None
        if isinstance(payload, dict) and payload.get("record_type") == "reference_windows":
            matches = [window for window in payload["windows"]
                       if (window["contig"], window["start"], window["end"]) == coordinates]
            if payload["genome_build"] != bundle["assembly"] or len(matches) != 1:
                raise PopulationError("Reference-window source does not resolve the exact assembly and interval")
            sequence = matches[0]["sequence"]
        elif isinstance(payload, dict) and {"genome", "chrom", "start", "end", "dna"} <= payload.keys():
            assembly = "GRCh38" if payload["genome"] == "hg38" else payload["genome"]
            if assembly != bundle["assembly"] or (payload["chrom"], payload["start"], payload["end"]) != coordinates:
                raise PopulationError("UCSC reference source has a different assembly or interval")
            sequence = payload["dna"].upper()
        else:
            raise PopulationError("Unsupported indexed reference-window source format")
        if sequence != actual[coordinates]["sequence"]:
            raise PopulationError("Reconstructed reference window differs from its pinned source bytes")


def prepare_inputs(study):
    """Bridge a pinned V3 genomic source configuration into comparison inputs.

    Every nested path, including a FASTA index read implicitly by V3, must already
    be bound by an explicit comparison input reference. No acquisition occurs.
    """
    verify_sources(study)
    transformations = []
    for source in study["sources"].values():
        if source["role"] == "source_transformation" and source["format"] == "json":
            candidate = read_json(source["path"])
            if isinstance(candidate, dict) and candidate.get("record_type") == "population_study_config":
                transformations.append((source, candidate))
    if len(transformations) != 1:
        raise PopulationError("Comparison search requires exactly one pinned population-study genomic source configuration")
    transformation, source_config = transformations[0]
    source_config = validate_population_config(source_config)
    config = study["config"]
    if (source_config["mode"], source_config["coverage_basis"]) != (config["mode"], config["coverage_basis"]):
        raise PopulationError("Genomic source and comparison mode/coverage basis differ")
    directory = Path(transformation["path"]).parent
    pinned_paths = {Path(source["path"]).resolve(): source for source in study["sources"].values()}
    paths = {}
    expected_roles = {"reference": "reference", "variants": "variants", "donors": "donor_map", "callable_intervals": "callability"}
    for key, name in source_config["inputs"].items():
        if name is None:
            paths[key] = None
            continue
        path = (directory / name).resolve()
        declared = pinned_paths.get(path)
        if declared is None or checksum(path) != declared["sha256"]:
            raise PopulationError("Nested genomic input is not explicitly checksum-bound: " + key)
        if key in expected_roles and declared["role"] != expected_roles[key]:
            raise PopulationError("Genomic input has the wrong comparison dependency role: " + key)
        expected = source_config["input_sources"][key]["sha256"]
        if expected is not None and expected != declared["sha256"]:
            raise PopulationError("Nested genomic source checksum disagrees with comparison input")
        paths[key] = path
    reference_index = Path(str(paths["reference"]) + ".fai")
    if reference_index.exists() and reference_index.resolve() not in pinned_paths:
        raise PopulationError("Auto-read FASTA index requires its own pinned comparison input")
    if sum(interval["end"] - interval["start"] for interval in source_config["scope"]["materialized_intervals"]) > source_config["resource_limits"]["max_reference_bases"]:
        raise ResourceLimit("Genomic source exceeds its pinned reference-base ceiling")
    reference = Reference(paths["reference"], source_config["reference_policy"], config["coverage_basis"] == "synthetic_truth")
    catalog = validate_catalog(config["catalogs"])
    domain = catalog[config["policies"]["search"]["domain_id"]]
    bundle = catalog[domain["reference_bundle_id"]]
    if bundle["assembly"] != reference.build:
        raise PopulationError("Comparison reference bundle assembly differs from the actual genomic source")
    if paths["reference"].suffix == ".json":
        _validate_regional_bundle(study, paths["reference"], bundle)
    elif bundle["identity_kind"] != "complete_reference_files" or bundle["complete_reference_sha256"]["value"] != reference.sha256:
        raise PopulationError("Complete reference bytes differ from their canonical bundle")
    for interval in domain["intervals"]:
        if not any(contains(materialized, interval) for materialized in source_config["scope"]["materialized_intervals"]):
            raise PopulationError("Shared comparison domain exceeds the pinned materialized scope")
        reference.fetch(interval["contig"], interval["start"], interval["end"])
    donors = canonical_donors(donors_from_tsv(paths["donors"], config["coverage_basis"]))
    if len(donors) > source_config["resource_limits"]["max_donors"]:
        raise ResourceLimit("Donor count exceeds the pinned genomic-source ceiling")
    by_sample = {}
    for donor in (row for row in catalog.values() if row["record_type"] == "comparison_donor"):
        for sample in donor["sample_ids"]:
            if sample in by_sample:
                raise PopulationError("One genotype sample maps to multiple canonical donors")
            by_sample[sample] = donor
    if set(by_sample) != {row["genotype_sample_id"] for row in donors}:
        raise PopulationError("Canonical donor catalog and pinned genomic panel have different sample membership")
    donor_records = {row["donor_id"]: by_sample[row["genotype_sample_id"]] for row in donors}
    source_records = []
    for role in ("variants", "donors"):
        source = pinned_paths[paths[role]]
        metadata = source_config["input_sources"][role]
        source_records.append(make_record("source_artifact", sha256=source["sha256"],
            source_id=metadata["source_id"], release=metadata["release"],
            uri=metadata["url"] or "urn:bes:local-genomic-source:" + role,
            media_type="text/vcf" if role == "variants" else "text/tab-separated-values", scope="complete_file",
            parent_artifact_id=unknown("original pinned genomic input"),
            locator={"kind": "file", "member": paths[role].name},
            access="synthetic" if config["mode"] == "synthetic_test" else "restricted" if config["mode"] == "research_sample_planning" else "public",
            use_terms=metadata["use_notes"]))
    auxiliary = []
    for source in study["sources"].values():
        if source["role"] != "reference_validation":
            continue
        if source["format"] != "json":
            raise PopulationError("Auxiliary REF validation requires explicit reference-window JSON")
        data = read_json(source["path"])
        contigs = sorted({row["contig"] for row in data.get("windows", [])})
        policy = {"genome_build": bundle["assembly"], "contig_alias_map": {contig: contig for contig in contigs},
                  "permitted_contigs": contigs}
        checked = Reference(source["path"], policy)
        for row in data["windows"]:
            if any(overlaps(row, interval) for interval in domain["intervals"]):
                raise PopulationError("Auxiliary REF validation sequence must remain outside the searched domain")
        auxiliary.append((source, checked))
    validation_reference = _ValidationReference(reference, auxiliary)
    variants, _, _ = load_vcf(paths["variants"], source_config, validation_reference, donors, study["snapshot_run_id"])
    callability = callable_intervals(paths["callable_intervals"], donors, reference,
                                    source_config["input_sources"].get("callable_intervals"))
    legacy_guides = load_guides(paths["guides"])
    legacy_intended = load_intended_loci(paths["intended_loci"])
    prepared = PreparedInputs(study, source_config, reference, variants, donors, donor_records,
        source_records, source_records[0]["record_id"], callability, legacy_guides, legacy_intended, bundle["record_id"])
    prepared.reference_validation_audit = sorted({canonical_bytes(row): row for row in validation_reference.checked_intervals}.values(), key=canonical_bytes)
    validate_intended_loci(prepared)
    return prepared


def validate_intended_loci(prepared):
    """Check source coordinates and stated edits without requiring nomination admission."""
    for item in prepared.profiles:
        design, editor = item["design"], item["editor"]
        geometry = editor["guide_geometry"]
        if set(design["guide_sequence"]) - set("ACGT"):
            raise PopulationError("Active guide sequence requires known ACGT bases")
        for identifier in design["intended_locus_ids"]:
            locus = prepared.catalog[identifier]
            if locus["reference_bundle_id"] != prepared.reference_bundle_id:
                raise PopulationError("Intended locus refers to a different reference bundle")
            if not any(contains(interval, locus) for interval in prepared.legacy_config["scope"]["materialized_intervals"]):
                raise PopulationError("Intended locus is outside the pinned genomic source scope")
            sequence = prepared.reference.fetch(locus["contig"], locus["start"], locus["end"])
            path = [{"kind": "reference", "pos0": locus["start"] + i, "base": base} for i, base in enumerate(sequence)]
            path = oriented_path(path, locus["strand"])
            if locus["path"] and locus["path"] != path:
                raise PopulationError("Intended-locus source path disagrees with the pinned reference")
            if len(sequence) != geometry["spacer_length"] + geometry["pam_length"]:
                raise PopulationError("Intended locus source geometry differs from its exact editor definition")
        for edit in design["desired_edit_definitions"]:
            locus = prepared.catalog[edit["locus_id"]]
            base = prepared.reference.fetch(locus["contig"], edit["pos0"], edit["pos0"] + 1)
            if edit["strand"] == "-":
                base = reverse_complement(base)
            if edit["reference_base"] != base:
                raise PopulationError("Desired-edit source base differs from the pinned reference")
    return prepared
