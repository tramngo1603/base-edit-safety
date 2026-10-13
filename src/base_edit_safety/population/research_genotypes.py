"""Stream local variant models without inferring reference calls or sex ploidy."""
from __future__ import annotations

import gzip
import re
from collections import deque

from .common import PopulationError, ResourceLimit, checksum, identity, record
from .haplotypes import reconstruct, research_variant_plan
from .variants import allele_edit, parse_call


def _ploidy_rows(rows, reference):
    result = {}
    for row in rows:
        contig, start, end, copies = (row[k] for k in ("contig", "start", "end", "copies"))
        if type(start) is not int or type(end) is not int or type(copies) is not int \
                or start < 0 or end <= start or copies not in (0, 1, 2):
            raise PopulationError("Ploidy intervals require nonnegative half-open bounds and 0, 1 or 2 copies")
        bounds = reference.bounds(contig)
        if not any(lo <= start < end <= hi for lo, hi in bounds):
            raise PopulationError("Ploidy interval lies outside the reference")
        result.setdefault(contig, []).append(dict(row))
    for contig, intervals in result.items():
        intervals.sort(key=lambda r: r["start"])
        merged = []
        for row in intervals:
            if merged and row["start"] < merged[-1]["end"]:
                raise PopulationError("Overlapping ploidy intervals")
            if merged and row["start"] == merged[-1]["end"] and row["copies"] == merged[-1]["copies"]:
                merged[-1]["end"] = row["end"]
            else:
                merged.append(dict(row))
        result[contig] = merged
    return result


def _interval_ploidy(rows, contig, start, end):
    for row in rows.get(contig, []):
        if row["start"] <= start and end <= row["end"]:
            return row["copies"]
    return None


def _validate_normalization(reference, contig, pos, ref, alts, edits):
    """Reject non-normalized literal alleles, checking multiallelic indels jointly."""
    if "*" in alts:
        # The marker is tied to this source position. Moving one literal ALT
        # independently would change its relationship to that marker.
        return
    symbolic = any(set(a) - set("ACGT") for a in alts)
    if len(alts) > 1 and not symbolic:
        alleles = [ref] + alts
        has_indel = any(e and e["kind"] in {"insertion", "deletion"} for e in edits[1:])
        suffix = len({a[-1] for a in alleles}) == 1
        prefix = len({a[0] for a in alleles}) == 1
        trim = min(map(len, alleles)) > 1
        if has_indel and ((suffix and (pos > 0 or trim)) or (prefix and trim)):
            raise PopulationError("Indel record is not jointly normalized against the reference")
    else:
        for edit in edits[1:]:
            if edit and edit["kind"] in {"insertion", "deletion"}:
                changed = edit["sequence"] if edit["kind"] == "insertion" else reference.fetch(contig, edit["start"], edit["end"])
                previous = reference.fetch(contig, edit["start"]-1, edit["start"]) if edit["start"] > pos else None
                if changed and previous == changed[-1]:
                    raise PopulationError("Indel is not left normalized against the reference")


def _rows(path, reference, sample, policy, ploidy, max_rows, max_allele_length, check_limits, excluded_contigs):
    source = "vcf:" + checksum(path, check_limits=check_limits)
    header, build, seen_contigs, previous, row_count = None, False, set(), None, 0
    aliases_seen = {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            check_limits()
            if line.startswith("##reference="):
                if build or line.rstrip().split("=", 1)[1] != reference.build:
                    raise PopulationError("VCF requires exactly one matching reference declaration")
                build = True
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM\t"):
                fields = line.rstrip("\r\n").split("\t")
                if header is not None or fields[:9] != ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]:
                    raise PopulationError("Invalid or duplicate VCF sample header")
                header = fields[9:]
                if len(header) != len(set(header)) or sample not in header:
                    raise PopulationError("Requested sample is missing or ambiguous in VCF")
                sample_column = 9 + header.index(sample)
                continue
            if not line.strip():
                continue
            if header is None or not build:
                raise PopulationError("VCF requires reference and sample headers before records")
            row_count += 1
            if row_count > max_rows:
                raise ResourceLimit("max_variant_rows exceeded during streaming ingestion")
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9 + len(header):
                raise PopulationError("VCF row width does not match header")
            chrom, pos_text, name, ref, alt_text, qual, filt, info_text, fmt = fields[:9]
            contig = reference.aliases.get(chrom)
            if contig is None:
                raise PopulationError("VCF contig is absent from the explicit alias map")
            if aliases_seen.setdefault(contig, chrom) != chrom:
                raise PopulationError("VCF contig alias collision")
            if not re.fullmatch(r"[0-9]+", pos_text) or int(pos_text) < 1:
                raise PopulationError("VCF POS must be a positive integer")
            pos = int(pos_text)-1
            if previous is not None:
                if contig == previous[0] and pos < previous[1] or contig != previous[0] and contig in seen_contigs:
                    raise PopulationError("VCF must be sorted with each contig in one block")
            seen_contigs.add(contig)
            previous = contig, pos
            if not ref or set(ref) - set("ACGTN"):
                raise PopulationError("VCF REF must be a literal reference allele")
            if reference.fetch(contig, pos, pos+len(ref)) != ref:
                raise PopulationError("VCF REF disagrees with pinned reference")
            alts = alt_text.split(",")
            if any(not a for a in alts) or len(alts) != len(set(alts)):
                raise PopulationError("VCF ALT indexes are empty or duplicated")
            formats = fmt.split(":")
            if any(not key for key in formats) or len(formats) != len(set(formats)):
                raise PopulationError("Empty or duplicate VCF FORMAT key")
            info = {}
            for item in info_text.split(";"):
                if item == ".":
                    continue
                key, value = item.split("=", 1) if "=" in item else (item, True)
                if not key or key in info:
                    raise PopulationError("Empty or duplicate VCF INFO key")
                info[key] = value
            if "END" in info and (not isinstance(info["END"], str) or not re.fullmatch(r"[0-9]+", info["END"]) or int(info["END"]) < pos+len(ref)):
                raise PopulationError("VCF END must bound the complete reference allele")
            end = max(pos+len(ref), int(info.get("END", pos+len(ref))))
            if not any(lo <= pos < end <= hi for lo, hi in reference.bounds(contig)):
                raise PopulationError("VCF footprint exceeds reference bounds")
            copies = _interval_ploidy(ploidy, contig, pos, pos+len(ref))
            sample_fields = fields[sample_column].split(":")
            gt_column = formats.index("GT") if "GT" in formats else len(sample_fields)
            raw_gt = sample_fields[gt_column] if gt_column < len(sample_fields) else "."
            raw_indexes = re.split(r"[/|]", raw_gt)
            if any(value != "." and (not re.fullmatch(r"[0-9]+", value) or int(value) > len(alts))
                   for value in raw_indexes) or "/" in raw_gt and "|" in raw_gt:
                raise PopulationError("Invalid VCF genotype syntax or allele index")
            call = parse_call(fields[sample_column], formats, filt, policy, source, sample, contig,
                              expected_ploidy=copies if copies in (1, 2) else 2)
            indexes = call["selected_indexes"]
            if call["call_status"] == "invalid" or indexes is not None and any(i > len(alts) for i in indexes):
                raise PopulationError("Invalid VCF genotype syntax or allele index")
            edits = [None] + [allele_edit(pos, ref, alt, max_allele_length) for alt in alts]
            _validate_normalization(reference, contig, pos, ref, alts, edits)
            unbounded_indexes = [i for i, alt in enumerate(alts, 1)
                                 if alt == "*" or alt != "." and set(alt)-set("ACGTN") and (
                                     "END" not in info or any(mark in alt for mark in "[].")
                                     or alt.startswith(("<BND", "<TRA")))]
            unbounded = bool(unbounded_indexes)
            # GT can identify an uncarried allele without establishing the
            # number of chromosome copies. Rejected or missing calls cannot.
            unknown_ploidy_indexes = (unbounded_indexes if call["call_status"] == "filtered" or "." in raw_indexes
                                      else sorted(set(map(int, raw_indexes)).intersection(unbounded_indexes)))
            for index, edit in enumerate(edits[1:], 1):
                if edit is None or edit["kind"] != "substitution":
                    continue
                left, right = 0, len(ref)
                while left < right and ref[left] == alts[index-1][left]:
                    left += 1
                while right > left and ref[right-1] == alts[index-1][right-1]:
                    right -= 1
                edits[index] = {"start": pos+left, "end": pos+right,
                                "sequence": alts[index-1][left:right], "kind": "substitution"}
            if copies not in (1, 2):
                call.update(call_status="unknown_ploidy" if copies is None else "zero_ploidy", selected_indexes=None,
                            phase_status="unsupported", phase_key=None,
                            reasons=["unknown_ploidy" if copies is None else "zero_ploidy"])
            vid = identity("var", [reference.sha256, contig, pos, ref, alts])
            call["call_id"] = identity("call", [source, vid, sample, call])
            counted = call["selected_indexes"] if call["call_status"] == "complete" else []
            counts = [counted.count(i) for i in range(1, len(alts)+1)]
            variant = record("normalized_variant", source, variant_id=vid, source_id=source,
                             source_record_index=row_count, original_id=name, contig=contig, pos0=pos,
                             end=end, unbounded_footprint=unbounded, ref=ref, original_alts=alts,
                             unbounded_allele_indexes=unbounded_indexes,
                             **({"unknown_ploidy_unbounded_indexes": unknown_ploidy_indexes} if copies is None else {}),
                             edits=edits, calls={sample: call}, normalization="validated_literal_alleles",
                             frequency={"AC": counts, "AN": len(counted),
                                        "AF": [n/len(counted) if counted else None for n in counts],
                                        "evaluable_donors": int(bool(counted)), "subset": "declared_sample"})
            yield variant, row_count
        if header is None or not build:
            raise PopulationError("VCF requires reference and sample headers")


def _with_star_bounds(rows, reference, sample, ploidy, max_pending_rows):
    """Bound unresolved deletions by bases called present on every copy.

    A missing VCF row is never an anchor. Retaining each interval between
    accepted anchors keeps bounded gaps and their neighboring variants in one
    reconstruction, even when the input records are far apart.
    """
    queued, pending, previous_anchor = deque(), [], None
    previous_region, anchor_variant = None, None
    for variant, row in rows:
        contig, pos = variant["contig"], variant["pos0"]
        region = next(((contig, r["start"], r["end"], r["copies"])
                       for r in ploidy.get(contig, [])
                       if r["start"] <= pos and pos+len(variant["ref"]) <= r["end"]), None)
        if region != previous_region or contig != (previous_region[0] if previous_region else None):
            pending.clear()
            while queued:
                yield queued.popleft()
            previous_anchor, anchor_variant = None, None
            previous_region = region
        call = variant["calls"][sample]
        indexes = call["selected_indexes"]
        alleles = [variant["ref"]] + variant["original_alts"]
        anchor = call["call_status"] == "complete" and region is not None and region[3] in (1, 2) \
            and all(len(alleles[i]) == len(variant["ref"]) and not set(alleles[i])-set("ACGT") for i in indexes)
        if anchor:
            remaining = []
            for star, left, left_variant in pending:
                if pos >= star["pos0"] + len(star["ref"]):
                    star["spanning_deletion_bounds"] = {
                        "start": left["pos0"]+1, "end": pos,
                        "basis": "bounded_by_flanking_called_bases",
                        "left_anchor": left,
                        "right_anchor": {"pos0": pos, "variant_id": variant["variant_id"],
                                         "source_record_index": row},
                    }
                    left_variant["_cluster_through_end"] = max(pos, left_variant.get("_cluster_through_end", pos))
                else:
                    remaining.append((star, left, left_variant))
            pending = remaining
        star_indexes = {i+1 for i, alt in enumerate(variant["original_alts"]) if alt == "*"}
        needs_bounds = star_indexes and (call["call_status"] != "complete" or star_indexes.intersection(indexes))
        if needs_bounds and call["call_status"] == "complete":
            prior = [v for v, _ in queued] + [variant]
            scope = {"contig": contig, "start": min(v["pos0"] for v in prior),
                     "end": max(v["end"] for v in prior)}
            proof = research_variant_plan(reference, scope, prior, sample, len(indexes))
            if "resolved_spanning_deletion" in proof["variant_dispositions"].get(variant["variant_id"], []):
                needs_bounds = False
        if needs_bounds \
                and previous_anchor is not None and previous_anchor["pos0"] < pos:
            pending.append((variant, previous_anchor, anchor_variant))
        if anchor and not pending:
            while queued:
                yield queued.popleft()
        if anchor:
            previous_anchor = {"pos0": pos+len(variant["ref"])-1, "variant_id": variant["variant_id"],
                               "source_record_index": row}
            anchor_variant = variant
        queued.append((variant, row))
        if len(queued) > max_pending_rows:
            raise ResourceLimit("Pending spanning-deletion records exceed the cluster record limit")
    while queued:
        yield queued.popleft()


def iter_variant_clusters(vcf_path, reference, sample, policy, ploidy_intervals, *,
                          max_rows=8000000, max_cluster_bases=100000, max_allele_length=None,
                          check_limits=lambda: None, excluded_contigs=()):
    """Yield bounded variant neighborhoods independently of reference nominations.

    Each cluster retains input row numbers, original ALT indexes and sample calls.
    Its reference padding includes 22 bases beyond the complete outer REF
    footprints, including deleted bases. A zero copy count marks an excluded cluster; exclusions distinguish
    absent copies from unknown or changing ploidy. Results remain conditional on
    the retained callset. Exhaust the iterator before declaring sample completion.
    """
    if max_allele_length is None:
        max_allele_length = max_cluster_bases
    if any(type(value) is not int or value < 1 for value in (max_rows, max_cluster_bases, max_allele_length)):
        raise PopulationError("Streaming limits must be positive integers")
    check_limits()
    ploidy = _ploidy_rows(ploidy_intervals, reference)
    if any(not isinstance(c, str) or c not in reference.aliases.values() for c in excluded_contigs):
        raise PopulationError("Excluded sample contig is outside the reference")
    excluded_contigs = frozenset(excluded_contigs)
    cluster = None

    def finish(current):
        interval = current["interval"]
        boundaries = {interval["start"], interval["end"]}
        for row in ploidy.get(interval["contig"], []):
            boundaries.update(x for x in (row["start"], row["end"]) if interval["start"] < x < interval["end"])
        points = sorted(boundaries)
        for start, end in zip(points, points[1:]):
            scope = {"contig": interval["contig"], "start": start, "end": end}
            copies = _interval_ploidy(ploidy, interval["contig"], start, end)
            included = [v for v in current["variants"]
                        if min(v["pos0"], v.get("spanning_deletion_bounds", {}).get("start", v["pos0"])) < end
                        and start < max(v["end"], v.get("spanning_deletion_bounds", {}).get("end", v["end"]))]
            result = {"interval": scope, "variants": included, "copies": copies if copies is not None else 0,
                      "exclusions": [], "source_rows": [n for n, pos in current["row_positions"] if start <= pos < end],
                      "accounted_variant_ids": [v["variant_id"] for v in included if start <= v["pos0"] < end],
                      "sample_contig_excluded": interval["contig"] in excluded_contigs}
            if copies is None:
                result["exclusions"].append({"reason": "unknown_or_changing_ploidy", **scope})
            elif copies == 0:
                result["exclusions"].append({"reason": "zero_ploidy", **scope})
            plan = research_variant_plan(reference, scope, included, sample, result["copies"])
            result["unbounded_exclusions"] = plan["unbounded_exclusions"]
            result["resource_excluded"] = end-start > max_cluster_bases
            if result["resource_excluded"]:
                result["exclusions"].append({**scope, "reason": "local_reconstruction_resource_limit",
                                              "limit": "max_cluster_bases", "maximum": max_cluster_bases})
            result["variant_dispositions"] = {vid: states for vid, states in plan["variant_dispositions"].items()
                                                if vid in result["accounted_variant_ids"]}
            if result["resource_excluded"]:
                result["variant_dispositions"] = {vid: ["unsupported"] for vid in result["accounted_variant_ids"]}
            result["variant_accounting"] = {"variant_records": len(result["accounted_variant_ids"])}
            for state in ("applied", "reference", "resolved_spanning_deletion", "unsupported"):
                result["variant_accounting"][state + "_variant_records"] = sum(
                    state in states for states in result["variant_dispositions"].values())
            yield result

    rows = _rows(vcf_path, reference, sample, policy, ploidy, max_rows, max_allele_length, check_limits, excluded_contigs)
    for variant, row_number in _with_star_bounds(rows, reference, sample, ploidy, max_cluster_bases):
        contig, pos, end = variant["contig"], variant["pos0"], variant["end"]
        bounds = next((b for b in reference.bounds(contig) if b[0] <= pos < end <= b[1]), None)
        if bounds is None:
            raise PopulationError("Variant footprint is outside reference bounds")
        footprint = variant.get("spanning_deletion_bounds", {})
        lo = max(bounds[0], min(pos, footprint.get("start", pos))-22)
        hi = min(bounds[1], max(end, footprint.get("end", end), variant.pop("_cluster_through_end", end))+22)
        if cluster is not None and (cluster["interval"]["contig"] != contig or lo > cluster["interval"]["end"]):
            yield from finish(cluster)
            cluster = None
        if cluster is None:
            cluster = {"interval": {"contig": contig, "start": lo, "end": hi}, "variants": [],
                       "copies": 0, "exclusions": [], "source_rows": [], "row_positions": [], "by_id": {}}
        existing = cluster["by_id"].get(variant["variant_id"])
        if existing is not None:
            if existing["calls"] != variant["calls"]:
                raise PopulationError("Discordant duplicate VCF record")
        else:
            cluster["variants"].append(variant)
            cluster["by_id"][variant["variant_id"]] = variant
        cluster["source_rows"].append(row_number)
        cluster["row_positions"].append((row_number, pos))
        cluster["interval"]["start"] = min(cluster["interval"]["start"], lo)
        cluster["interval"]["end"] = max(cluster["interval"]["end"], hi)
        if len(cluster["source_rows"]) > max_cluster_bases:
            raise ResourceLimit("Cluster record limit exceeded; sample variant search is incomplete")
    check_limits()
    if cluster is not None:
        yield from finish(cluster)


def materialize_cluster(reference, cluster, sample, run_id):
    """Build only the declared copies, retaining local phase and excluded calls."""
    if cluster["copies"] not in (0, 1, 2):
        raise PopulationError("Cluster copies must be explicit")
    if cluster.get("sample_contig_excluded") or cluster.get("resource_excluded"):
        return []
    plan = research_variant_plan(reference, cluster["interval"], cluster["variants"], sample, cluster["copies"])
    return [reconstruct(reference, cluster["interval"], cluster["variants"], sample, copy, run_id, plan=copy_plan)
            for copy, copy_plan in enumerate(plan["copies"])]
