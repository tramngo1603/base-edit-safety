"""Exhaustive ungapped scanning of bounded reconstructed sequences on both strands."""
from __future__ import annotations

from collections import Counter, defaultdict

from .common import ResourceLimit, identity, record
from .haplotypes import local_status
from .identity import bounding_interval, candidate_identity, definition_hashes, pam_matches, reverse_complement
from .variants import canonical_donors


def evaluation_span(path, material):
    """Include target bases and the REF spans of variants establishing inserted bases."""
    lo, hi = bounding_interval(path)
    boundaries = {p["boundary0"] for p in path if p["kind"] == "insertion"}
    for variant in material["variant_intervals"] if boundaries else []:
        if boundaries.intersection(variant["insertion_boundaries"]):
            lo, hi = min(lo, variant["start"]), max(hi, variant["end"])
    return lo, hi


def sequence_evaluable(material, lo, hi):
    ok, phase = local_status(material, lo, hi)
    if ok and any(p["kind"] == "reference" and lo <= p["pos0"] < hi and p["base"] == "N"
                  for p in material["sequence_path"]):
        return False, "reference_sequence_unknown"
    return ok, phase


def supported(config, donor_id, contig, lo, hi, callable_rows):
    if config["coverage_basis"] == "synthetic_truth":
        return True
    cursor = lo
    for row in sorted(callable_rows, key=lambda r: r["start_0based"]):
        if row["donor_id"] != donor_id or row["contig"] != contig:
            continue
        if row["start_0based"] <= cursor:
            cursor = max(cursor, row["end_0based_exclusive"])
        if cursor >= hi:
            return True
    return False


def role_for(contig, bounds, strand, guide, intended):
    for locus in intended:
        if (locus["guide_id"] == guide["guide_id"] and locus["chrom"] == contig
                and locus["strand"] == strand and bounds == [locus["start_0based"], locus["end_0based_exclusive"]]):
            return {"HBG1": "intended_hbg1", "HBG2": "intended_hbg2"}.get(locus["locus_name"], "intended_control")
    return "candidate_offtarget"


def window_key(path):
    """Identify a forward path independently of its bases or nomination status."""
    first = path[0].get("pos0")
    if first is not None and all(p["kind"] == "reference" and p["pos0"] == first+i for i, p in enumerate(path)):
        return first
    return tuple(("r", p["pos0"]) if p["kind"] == "reference" else
                 ("i", p["boundary0"], p["offset0"]) for p in path)


def pairing_status(config, donor_id, contig, ref_status, statuses, callable_rows):
    """Use the same exact-path eligibility for scope counts and nomination comparisons."""
    if ref_status is None:
        return False, "variant_only_indel_path"
    if not ref_status[0]:
        return False, "reference_" + ref_status[1]
    if len(statuses) != 2 or any(s is None for s in statuses):
        return False, "donor_indel_path_not_shared"
    if any(not s[0] for s in statuses):
        return False, next(s[1] for s in statuses if not s[0])
    if config["coverage_basis"] == "callable_supported" and not supported(
            config, donor_id, contig, ref_status[2], ref_status[3], callable_rows):
        return False, "reference_callability_unavailable"
    return True, "paired_exact_path"


def window_assessments(config, material):
    """Account for all requested paths before PAM or mismatch admission."""
    windows = {}
    path, sequence = material["sequence_path"], material["sequence_genomic_forward"]
    for offset in range(max(0, len(path)-22)):
        target = path[offset:offset+23]
        lo, hi = evaluation_span(target, material)
        if not any(s["contig"] == material["contig"] and s["start"] <= lo and hi <= s["end"]
                   for s in config["scope"]["requested_intervals"]):
            continue
        key = window_key(target)
        valid, reason = local_status(material, lo, hi)
        if valid and "N" in sequence[offset:offset+23]:
            valid, reason = False, "sequence_unknown"
        windows[key] = (valid, reason, lo, hi)
    return windows


def paired_window_scope(config, materials, guides, donors, callable_rows):
    """Compare exact paths on the reference and both copies, including windows with no hit."""
    partitions = defaultdict(list)
    for m in materials:
        interval = m["interval"]
        partitions[(m["contig"], interval["start"], interval["end"])].append(m)
    rows = []
    for (contig, start, end), group in sorted(partitions.items()):
        ref = next(m for m in group if m["donor_id"] is None)
        reference_windows = window_assessments(config, ref)
        for donor in donors:
            copies = [window_assessments(config, m) for m in group if m["donor_id"] == donor["donor_id"]]
            donor_paths = set().union(*(set(w) for w in copies))
            universe = set(reference_windows) | donor_paths
            excluded, paired_starts = Counter(), []
            for key in universe:
                ref_status = reference_windows.get(key)
                statuses = [w.get(key) for w in copies]
                paired, reason = pairing_status(config, donor["donor_id"], contig, ref_status, statuses, callable_rows)
                if paired:
                    paired_starts.append(key)
                    continue
                excluded[reason] += 2
            ranges = []
            for pos in sorted(paired_starts):
                if ranges and ranges[-1]["end"] == pos:
                    ranges[-1]["end"] += 1
                else:
                    ranges.append({"start": pos, "end": pos+1})
            for guide in sorted(guides, key=lambda g: g["guide_id"]):
                rows.append({"donor_id": donor["donor_id"], "guide_id": guide["guide_id"],
                    "partition": {"contig": contig, "start": start, "end": end},
                    "reference_oriented_windows": 2*len(reference_windows),
                    "variant_oriented_windows": 2*len(donor_paths), "union_oriented_windows": 2*len(universe),
                    "paired_evaluable_oriented_windows": 2*len(paired_starts),
                    "unpaired_oriented_windows": sum(excluded.values()),
                    "unpaired_reasons": dict(sorted(excluded.items())), "paired_reference_start_intervals": ranges})
    return {"coverage_basis": config["coverage_basis"], "rows": rows,
        "definition": "Per donor, guide and partition: union of exact reference/donor paths, counted on both strands before hit filtering. "
                      "Paired paths require an evaluable reference and both donor copies under the selected basis. "
                      "Start intervals are half-open ranges of reference window starts, not target spans; indel paths without an exact counterpart remain unpaired."}


def scan(config, reference, materials, guides, intended, donors, callable_rows, run_id, *, check_limits=lambda: None):
    donors = canonical_donors(donors)
    profile = config["search_profile"]
    profile_id = identity("profile", profile)
    candidates, hits, execution = {}, [], []
    hit_copies = defaultdict(set)
    candidate_sources = defaultdict(set)
    material_by_donor = defaultdict(list)
    for material in materials:
        check_limits()
        material_by_donor[material["donor_id"]].append(material)
        seq, path = material["sequence_genomic_forward"], material["sequence_path"]
        for guide in guides:
            counts = Counter()
            guide_sha, editor_sha = definition_hashes(guide)
            for offset in range(max(0, len(seq)-22)):
                if offset % 128 == 0:
                    check_limits()
                forward = seq[offset:offset+23]
                forward_path = path[offset:offset+23]
                lo, hi = evaluation_span(forward_path, material)
                if not any(s["contig"] == material["contig"] and s["start"] <= lo and hi <= s["end"]
                           for s in config["scope"]["requested_intervals"]):
                    counts["boundary_or_padding_windows"] += 2
                    continue
                valid, phase = local_status(material, lo, hi)
                if not valid or "N" in forward:
                    counts[phase if not valid else "sequence_unknown"] += 2
                    continue
                counts["evaluated_windows"] += 2
                for strand, target, target_path in (("+", forward, forward_path),
                                                    ("-", reverse_complement(forward), forward_path[::-1])):
                    mismatch_positions = [i+1 for i, (a, b) in enumerate(zip(guide["spacer_sequence"], target[:20])) if a != b]
                    if len(mismatch_positions) > profile["max_mismatches"] or not pam_matches(target[20:], profile["pam_patterns"]):
                        continue
                    cid, payload = candidate_identity(reference.sha256, material["contig"], strand, guide, target, target_path)
                    if cid not in candidates:
                        bounds = bounding_interval(target_path)
                        computational = guide.get("purpose") == "computational_only"
                        substrate = None if computational else {"ABE": "A", "CBE": "C"}.get(guide["editor_class"])
                        window = [] if computational else [guide["editing_window_start_1based"], guide["editing_window_end_1based_inclusive"]]
                        positions = [i for i in range(window[0], window[1]+1) if target[i-1] == substrate] if substrate else []
                        candidates[cid] = record("candidate_configuration", run_id, configuration_id=cid,
                            reference_id=reference.sha256, identity_payload=payload, contig=material["contig"], strand=strand,
                            target_path=payload["target_path"], bounding_interval=bounds, evaluation_interval=[lo, hi],
                            sequence_genomic_forward=forward, protospacer_guide_orientation=target[:20], pam_guide_orientation=target[20:],
                            guide_id=guide["guide_id"], guide_definition_id=guide_sha, editor_definition_id=editor_sha,
                            role="sequence_candidate" if computational else role_for(material["contig"], bounds, strand, guide, intended),
                            editability={"window": window, "substrate": substrate, "positions": positions,
                                         "status": "not_assessed_computational_only" if computational else "assessed" if substrate else "not_applicable"},
                            variant_links=[], legacy_site_links=[])
                    candidate_sources[cid].add(material["sequence_configuration_id"])
                    local_variants = {v["variant_id"] for v in material["variant_intervals"] if lo < v["end"] and v["start"] < hi}
                    candidates[cid]["variant_links"] = sorted(set(candidates[cid]["variant_links"]) | local_variants)
                    hit_copies[(material["donor_id"], cid)].add((material["interval"]["start"], material["local_copy"]))
                    hits.append(record("search_hit", run_id, configuration_id=cid, profile_id=profile_id,
                        donor_id=material["donor_id"], local_copy=material["local_copy"],
                        sequence_configuration_id=material["sequence_configuration_id"], partition=material["interval"],
                        mismatch_count=len(mismatch_positions), mismatch_positions_1based=mismatch_positions,
                        guide_alignment=guide["spacer_sequence"]+"NNN", target_alignment=target,
                        pam_assessment="matches_declared_profile", backend_score=None, score_scope=None,
                        phase_status=phase, execution_status="completed"))
                    counts["hits"] += 1
                    if len(hits) > config["resource_limits"]["max_hit_rows"]:
                        raise ResourceLimit("max_hit_rows exceeded; complete search not retained")
            execution.append({"donor_id": material["donor_id"], "local_copy": material["local_copy"],
                              "guide_id": guide["guide_id"], "interval": material["interval"],
                              "execution_status": "completed", "hit_count": counts["hits"],
                              "evaluated_windows": counts["evaluated_windows"], "counts": dict(sorted(counts.items()))})
    for cid, candidate in candidates.items():
        candidate["source_sequence_ids"] = sorted(candidate_sources[cid])
    carriage = []
    for cid, candidate in sorted(candidates.items()):
        check_limits()
        lo, hi = candidate["evaluation_interval"]
        for donor in donors:
            did = donor["donor_id"]
            applicable = [m for m in material_by_donor[did] if m["contig"] == candidate["contig"]
                          and m["interval"]["start"] <= lo and hi <= m["interval"]["end"]]
            statuses = [sequence_evaluable(m, lo, hi) for m in applicable]
            evaluated = len(applicable) == 2 and all(s[0] for s in statuses)
            copies = len(hit_copies[(did, cid)])
            call_support = evaluated and supported(config, did, candidate["contig"], lo, hi, callable_rows)
            row = record("donor_carriage", run_id, donor_id=did, configuration_id=cid,
                         evidence_supported_carriage=("confirmed" if copies else "absent") if call_support else "unknown",
                         model_carriage=("present" if copies else "absent") if evaluated else "unknown",
                         supported_copy_count=copies if call_support else None,
                         reference_evidence="synthetic_truth" if config["coverage_basis"] == "synthetic_truth" else
                         "callable_contract" if call_support else "reference_fill_conditional_on_retained_calls",
                         phase_status=sorted({s[1] for s in statuses}),
                         source_call_ids=sorted({c for m in applicable for c in m["source_call_ids"]}),
                         sequence_configuration_ids=sorted({m["sequence_configuration_id"] for m in applicable}),
                         model_id=identity("model", config["genotype_policy"]) if evaluated else None,
                         reasons=sorted({s[1] for s in statuses if not s[0]}) + ([] if call_support else ["complete_sequence_support_unavailable"]),
                         evaluated_scope={"contig": candidate["contig"], "start": lo, "end": hi})
            carriage.append(row)
    frequencies = []
    for cid in sorted(candidates):
        rows = [r for r in carriage if r["configuration_id"] == cid]
        for basis, key, yes, no in (("synthetic_truth" if config["coverage_basis"] == "synthetic_truth" else "callable_supported", "evidence_supported_carriage", "confirmed", "absent"),
                                     ("retained_callset_model", "model_carriage", "present", "absent")):
            n = sum(r[key] in {yes, no} for r in rows)
            k = sum(r[key] == yes for r in rows)
            frequencies.append({"configuration_id": cid, "basis": basis, "carrier_donors": k,
                                "evaluable_donors": n, "unknown_donors": len(rows)-n, "frequency": k/n if n else None})
    audit = record("search_audit", run_id, profile_id=profile_id, execution_status="completed", row_retention="all",
        configuration_coverage="supported_local_windows_only", search_completeness="bounded_supported_sequence_enumeration",
        requested_scope=config["scope"]["requested_intervals"], materialized_scope=config["scope"]["materialized_intervals"],
        omitted_scope={"outside_requested_intervals": "not_searched", "boundary_spanning_targets": "not_searched",
                       "unknown_or_unsupported_windows": "see execution counts"}, executions=execution,
        configuration_frequencies=frequencies, total_hit_rows=len(hits), distinct_configurations=len(candidates))
    audit["paired_window_scope"] = paired_window_scope(config, materials, guides, donors, callable_rows)
    return list(candidates.values()), hits, carriage, audit


def nomination_changes(config, candidates, materials, donors, guides, carriage, run_id, callable_rows=(), *, check_limits=lambda: None):
    """Compare exact short paths; keep unpaired and indel mappings explicitly unresolved."""
    result = []
    donors = canonical_donors(donors)
    lookup = {(r["donor_id"], r["configuration_id"]): r for r in carriage}
    guides = {g["guide_id"]: g for g in guides}
    refs = [m for m in materials if m["donor_id"] is None]
    windows = {}
    profile = config["search_profile"]
    for c in candidates:
        check_limits()
        for donor in donors:
            did = donor["donor_id"]
            lo, hi = c["evaluation_interval"]
            current = lookup[(did, c["configuration_id"])]
            mismatch_change, substrate_change = None, None
            basis_key = "model_carriage" if config["coverage_basis"] == "retained_callset_model" else "evidence_supported_carriage"
            ref = next((m for m in refs if m["contig"] == c["contig"] and m["interval"]["start"] <= lo and hi <= m["interval"]["end"]), None)
            pairing_reason = "reference_unavailable"
            paired = False
            if ref is not None:
                partition = ref["interval"]
                cache_key = (did, c["contig"], partition["start"], partition["end"])
                if cache_key not in windows:
                    copies = [m for m in materials if m["donor_id"] == did and m["contig"] == c["contig"] and m["interval"] == partition]
                    windows[cache_key] = (window_assessments(config, ref), [window_assessments(config, m) for m in copies])
                reference_windows, copies = windows[cache_key]
                path = c["target_path"] if c["strand"] == "+" else c["target_path"][::-1]
                key = window_key(path)
                paired, pairing_reason = pairing_status(config, did, c["contig"], reference_windows.get(key),
                    [copy.get(key) for copy in copies], callable_rows)
            if pairing_reason in {"variant_only_indel_path", "donor_indel_path_not_shared"}:
                change, counterpart = "indel_path_requires_separate_comparison", None
            elif not paired:
                change, counterpart = "unpaired_unevaluable", None
            else:
                counterpart = ref["sequence_genomic_forward"][lo-ref["interval"]["start"]:hi-ref["interval"]["start"]]
                if c["strand"] == "-":
                    counterpart = reverse_complement(counterpart)
                guide = guides[c["guide_id"]]
                ref_hit = pam_matches(counterpart[20:], profile["pam_patterns"]) and sum(a != b for a, b in zip(guide["spacer_sequence"], counterpart[:20])) <= profile["max_mismatches"]
                if current[basis_key] in {"present", "confirmed"}:
                    target = c["protospacer_guide_orientation"] + c["pam_guide_orientation"]
                    mismatch_change = sum(a != b for a, b in zip(guide["spacer_sequence"], target[:20])) - sum(a != b for a, b in zip(guide["spacer_sequence"], counterpart[:20]))
                    substrate = None if guide.get("purpose") == "computational_only" else {"ABE": "A", "CBE": "C"}.get(guide["editor_class"])
                    if substrate:
                        start, end = guide["editing_window_start_1based"]-1, guide["editing_window_end_1based_inclusive"]
                        substrate_change = target[start:end].count(substrate) - counterpart[start:end].count(substrate)
                    change = ("shared" if target == counterpart else "sequence_changed") if ref_hit else "pam_gained" if not pam_matches(counterpart[20:], profile["pam_patterns"]) else "mismatch_threshold_gained"
                elif ref_hit:
                    donor_sequences = []
                    for m in materials:
                        if m["donor_id"] == did and m["contig"] == c["contig"]:
                            points = [p for p in m["sequence_path"] if p["kind"] == "reference" and lo <= p["pos0"] < hi]
                            inserted = any(p["kind"] == "insertion" and lo < p["boundary0"] < hi for p in m["sequence_path"])
                            if len(points) == 23 and not inserted:
                                seq = "".join(p["base"] for p in points)
                                donor_sequences.append(seq if c["strand"] == "+" else reverse_complement(seq))
                    change = "pam_lost" if len(donor_sequences) == 2 and all(not pam_matches(s[20:], profile["pam_patterns"]) for s in donor_sequences) else "configuration_changed"
                else:
                    change = "absent_in_both_comparisons"
            result.append(record("nomination_change", run_id, donor_id=did, configuration_id=c["configuration_id"],
                                 change=change, paired_evaluable=paired, reference_counterpart=counterpart,
                                 pam_gained=change == "pam_gained", pam_lost=change == "pam_lost",
                                 mismatch_change=mismatch_change, substrate_count_change=substrate_change,
                                 pairing_reason=pairing_reason,
                                 comparison_scope="same_backend_guide_profile_and_intervals"))
    return result
