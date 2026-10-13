"""Reconstruct local chromosome copies while retaining every unresolved gap."""
from __future__ import annotations

import re

from .common import ResourceLimit, identity, record
from .variants import allele_edit, canonical_donors


def intersects(start, end, lo, hi):
    return start < hi and lo < end


def _edit_overlap(a, b):
    """Compare changed bases, keeping shared VCF anchors out of the conflict."""
    if a == b:
        return False
    if a["start"] == a["end"] and b["start"] == b["end"]:
        return a["start"] == b["start"]
    if a["start"] == a["end"]:
        return b["start"] < a["start"] < b["end"]
    if b["start"] == b["end"]:
        return a["start"] < b["start"] < a["end"]
    return intersects(a["start"], a["end"], b["start"], b["end"])


def _unbounded_event_exclusions(reference, variant, indexes):
    """Retain every admitted chromosome that an unresolved adjacency can affect."""
    contig = variant["contig"]
    targets, missing = {contig}, set()
    for index in indexes:
        alt = variant["original_alts"][index-1]
        if any(mark in alt for mark in "[].") or alt.startswith(("<BND", "<TRA")):
            mate = re.search(r"[\[\]]([^\[\]]+):[0-9]+[\[\]]", alt)
            if mate is not None:
                # A named partner outside the admitted reference is known,
                # though its sequence is outside this search's scope.
                if mate[1] in reference.aliases:
                    targets.add(reference.aliases[mate[1]])
                missing.add("complete rearranged allele sequence and breakpoint definitions")
            else:
                targets.update(reference.aliases.values())
                missing.add("supported breakend partner, breakpoint and rearranged allele definitions")
        elif alt == "*":
            missing.add("flanking accepted called bases or a compatible covering deletion")
        else:
            missing.add("a defensible affected interval and allele sequence")
    return [{"contig": target, "start": lo, "end": hi, "source_contig": contig,
             "reason": "unbounded_unsupported_allele", "variant_id": variant["variant_id"],
             "missing_information": "; ".join(sorted(missing))}
            for target in sorted(targets) for lo, hi in reference.bounds(target)]


def research_variant_plan(reference, interval, variants, donor_id, copies):
    """Resolve local allele relationships without changing the recorded genotypes.

    A star identifies an absent base on a chromosome copy, not an edit. A
    unique covering deletion therefore constrains the two records' relative
    orientation. Other phase blocks remain separate in window evaluation.
    """
    contig, start, end = (interval[k] for k in ("contig", "start", "end"))
    entries = {}
    parents, parity = {}, {}

    def group(v, call):
        indexes = call["selected_indexes"]
        if call["call_status"] != "complete" or len(indexes) != 2 or indexes[0] == indexes[1]:
            return None
        key = tuple(call["phase_key"]) if call["phase_key"] is not None else (v["variant_id"],)
        parents.setdefault(key, key)
        parity.setdefault(key, 0)
        return key

    def root(key):
        if parents[key] != key:
            ancestor, flip = root(parents[key])
            parity[key] ^= flip
            parents[key] = ancestor
        return parents[key], parity[key]

    def constrain(a, b, flip):
        ar, ap = root(a)
        br, bp = root(b)
        if ar == br:
            return ap ^ bp == flip
        parents[br], parity[br] = ar, ap ^ bp ^ flip
        return True

    for v in variants:
        call = v["calls"][donor_id]
        entries[v["variant_id"]] = {"variant": v, "call": call, "group": group(v, call),
                                      "stars": {}, "reason": None, "bounds": None, "orphan": False}
    constraints = []
    for vid, entry in entries.items():
        v, call = entry["variant"], entry["call"]
        indexes = call["selected_indexes"]
        star_indexes = {i + 1 for i, alt in enumerate(v["original_alts"]) if alt == "*"}
        if not star_indexes or call["call_status"] == "complete" and not star_indexes.intersection(indexes):
            continue
        carried = []
        for other_id, other in entries.items():
            if other_id == vid:
                continue
            ov, oc = other["variant"], other["call"]
            if ov["contig"] != contig or not ov["pos0"] < v["pos0"] < ov["pos0"] + len(ov["ref"]):
                continue
            covered = []
            for i, alt in enumerate(ov["original_alts"], 1):
                if not set(alt) - set("ACGT") and len(alt) < len(ov["ref"]):
                    edit = allele_edit(ov["pos0"], ov["ref"], alt, max(len(alt), len(ov["ref"])))
                    if edit is not None and edit["kind"] == "deletion" and edit["start"] <= v["pos0"] < edit["end"]:
                        covered.append(i)
            if oc["call_status"] == "complete" and any(i in covered for i in oc["selected_indexes"]):
                carried.append((other_id, other, [i in covered for i in oc["selected_indexes"]]))
        # One ancestor with matching dosage gives a unique local pairing. Keep
        # ambiguous overlapping ancestors bounded, without choosing one.
        resolved = False
        if call["call_status"] == "complete" and len(carried) == 1:
            ancestor_id, ancestor, deleted = carried[0]
            selected_stars = [i in star_indexes for i in indexes]
            if len(deleted) == len(selected_stars) and sum(deleted) == sum(selected_stars):
                if all(deleted):
                    resolved = True
                elif len(deleted) == 2:
                    flip = int(deleted[0] != selected_stars[0])
                    resolved = constrain(ancestor["group"], entry["group"], flip)
                    if resolved:
                        constraints.append({"variant_ids": [ancestor_id, vid],
                                            "relation": "spanning_deletion_copy_constraint"})
                if resolved:
                    entry["stars"] = {i: ancestor_id for i in star_indexes}
        if not resolved:
            entry["orphan"] = not carried
            entry["reason"] = ("unresolved_spanning_deletion" if call["call_status"] == "complete"
                               else call["call_status"])
            bound = v.get("spanning_deletion_bounds")
            if bound is not None:
                entry["bounds"] = (min(bound["start"], v["pos0"]),
                                   max(bound["end"], v["pos0"] + len(v["ref"])))

    def unbounded_indexes(entry):
        v, call = entry["variant"], entry["call"]
        possible = v.get("unbounded_allele_indexes", [])
        if call["call_status"] == "unknown_ploidy":
            possible = v.get("unknown_ploidy_unbounded_indexes", possible)
        elif call["call_status"] == "zero_ploidy":
            possible = []
        elif call["call_status"] == "complete":
            possible = sorted(set(possible).intersection(call["selected_indexes"]))
        return [i for i in possible if v["original_alts"][i-1] != "*"
                or entry["bounds"] is None and i not in entry["stars"]]

    unbounded, plans, dispositions = [], [], {vid: set() for vid in entries}
    for copy in range(copies):
        plan = {"gaps": [], "heterozygous_sites": [], "accepted": [], "calls": [],
                "changes": [], "phase_constraints": constraints}
        selected, present, reference_calls = [], [], set()
        for vid, entry in entries.items():
            v, call = entry["variant"], entry["call"]
            indexes = call["selected_indexes"]
            reason = entry["reason"]
            flip = root(entry["group"])[1] if entry["group"] is not None else 0
            index = indexes[copy ^ flip] if call["call_status"] == "complete" and len(indexes) == copies else None
            if entry["orphan"] and entry["bounds"] is not None and index is not None \
                    and (index == 0 or v["original_alts"][index-1] != "*"):
                reason = None
            if call["call_status"] != "complete":
                reason = call["call_status"]
            elif len(indexes) != copies:
                reason = "unsupported_ploidy"
            star = index in entry["stars"] if index is not None else False
            edit = v["edits"][index] if index else None
            if reason is None and index and edit is None and not star:
                reason = "unsupported_carried_allele"
            if reason is None and index and (v["pos0"] < start or v["pos0"] + len(v["ref"]) > end):
                reason = "variant_crosses_ploidy_boundary"
            if reason is None and star:
                ancestor = entries[entry["stars"][index]]
                ac = ancestor["call"]
                aflip = root(ancestor["group"])[1] if ancestor["group"] is not None else 0
                ai = ac["selected_indexes"][copy ^ aflip]
                if ancestor["variant"]["edits"][ai] is None:
                    reason = "spanning_deletion_ancestor_unsupported"
                    entry["bounds"] = (ancestor["variant"]["pos0"], ancestor["variant"]["end"])
            lo, hi = entry["bounds"] or (v["pos0"], v["end"] if reason else v["pos0"] + len(v["ref"]))
            possible_unbounded = unbounded_indexes(entry)
            is_unbounded = reason is not None and bool(possible_unbounded)
            if is_unbounded:
                for exclusion in _unbounded_event_exclusions(reference, v, possible_unbounded):
                    if exclusion not in unbounded:
                        unbounded.append(exclusion)
                lo, hi = start, end
            if v["contig"] != contig or not intersects(lo, hi, start, end):
                continue
            plan["calls"].append(call["call_id"])
            if call["call_status"] == "complete" and len(indexes) == 2 and indexes[0] != indexes[1]:
                phase_lo, phase_hi = entry["bounds"] if entry["orphan"] and entry["bounds"] is not None \
                    else (v["pos0"], v["pos0"] + len(v["ref"]))
                plan["heterozygous_sites"].append({"start": phase_lo, "end": phase_hi,
                                                  "phase_key": call["phase_key"], "variant_id": vid})
            if index == 0 or edit is not None:
                spans = [(v["pos0"], v["pos0"] + len(v["ref"]))]
                if edit is not None and edit["kind"] == "deletion":
                    spans = [(v["pos0"], edit["start"]), (edit["end"], v["pos0"] + len(v["ref"]))]
                present.append((vid, [(a, b) for a, b in spans if a < b]))
            if reason:
                gap = {"start": max(start, lo), "end": min(end, hi), "reason": reason, "variant_id": vid}
                if entry["bounds"] is not None and v.get("spanning_deletion_bounds") is not None:
                    gap["boundary_evidence"] = v["spanning_deletion_bounds"]
                plan["gaps"].append(gap)
                dispositions[vid].add("unsupported")
            elif star:
                dispositions[vid].add("resolved_spanning_deletion")
                plan["accepted"].append(v)
            elif index:
                selected.append((edit, vid))
                plan["accepted"].append(v)
            else:
                reference_calls.add(vid)
        conflicts = set()
        for i, (a, aid) in enumerate(selected):
            for b, bid in selected[i+1:]:
                if _edit_overlap(a, b):
                    conflicts.update((aid, bid))
            if a["kind"] == "deletion":
                for bid, spans in present:
                    if bid == aid or not any(intersects(a["start"], a["end"], lo, hi) for lo, hi in spans):
                        continue
                    ag, bg = entries[aid]["group"], entries[bid]["group"]
                    if ag is None or bg is None or root(ag)[0] == root(bg)[0]:
                        conflicts.update((aid, bid))
        seen_edits = set()
        for edit, vid in selected:
            if vid in conflicts:
                v = entries[vid]["variant"]
                plan["gaps"].append({"start": max(start, v["pos0"]), "end": min(end, v["pos0"] + len(v["ref"])),
                                     "reason": "overlapping_carried_edits", "variant_id": vid})
                dispositions[vid].add("unsupported")
                continue
            edit_key = (edit["start"], edit["end"], edit["kind"], edit["sequence"])
            if edit_key not in seen_edits:
                plan["changes"].append((edit, vid))
                seen_edits.add(edit_key)
            dispositions[vid].add("applied")
        selected_ids = {vid for _, vid in selected}
        for vid in conflicts - selected_ids:
            v = entries[vid]["variant"]
            plan["gaps"].append({"start": max(start, v["pos0"]), "end": min(end, v["pos0"] + len(v["ref"])),
                                 "reason": "overlapping_carried_edits", "variant_id": vid})
            dispositions[vid].add("unsupported")
        for vid in reference_calls - conflicts:
            dispositions[vid].add("reference")
        plan["changes"].sort(key=lambda item: (item[0]["start"], item[0]["end"], item[1]))
        plans.append(plan)
    if not copies:
        dispositions = {vid: {"unsupported"} for vid in entries}
        for entry in entries.values():
            possible_unbounded = unbounded_indexes(entry)
            if possible_unbounded:
                unbounded.extend(_unbounded_event_exclusions(reference, entry["variant"], possible_unbounded))
    accounting = {"variant_records": len(entries)}
    for state in ("applied", "reference", "resolved_spanning_deletion", "unsupported"):
        accounting[state + "_variant_records"] = sum(state in states for states in dispositions.values())
    return {"copies": plans, "unbounded_exclusions": unbounded, "variant_accounting": accounting,
            "variant_dispositions": {vid: sorted(states) for vid, states in dispositions.items()}}


def reconstruct(reference, interval, variants, donor_id, copy_index, run_id, *, plan=None):
    contig, start, end = interval["contig"], interval["start"], interval["end"]
    if plan is not None:
        gaps = plan["gaps"]
        heterozygotes = plan["heterozygous_sites"]
        accepted = plan["accepted"]
        calls = plan["calls"]
        changes = plan["changes"]
    else:
        gaps, heterozygotes, accepted, calls = [], [], [], []
        for v in variants if donor_id is not None else []:
            call = v["calls"][donor_id]
            reason = None
            indexes = call["selected_indexes"]
            if call["call_status"] != "complete":
                reason = call["call_status"]
            elif any(i != 0 and v["edits"][i] is None for i in indexes):
                reason = "unsupported_carried_allele"
            # Only an unknown or carried unsupported allele can have a global footprint.
            # A supported call uses its own contig and literal REF span.
            global_gap = reason is not None and v["unbounded_footprint"]
            footprint_end = v["end"] if reason else v["pos0"] + len(v["ref"])
            if not global_gap and (v["contig"] != contig or not intersects(v["pos0"], footprint_end, start, end)):
                continue
            calls.append(call["call_id"])
            lo, hi = (start, end) if global_gap else (max(start, v["pos0"]), min(end, footprint_end))
            if reason is None and any(indexes) and (v["pos0"] < start or v["pos0"] + len(v["ref"]) > end):
                reason = "variant_crosses_materialized_boundary"
            if reason:
                gaps.append({"start": lo, "end": hi, "reason": reason, "variant_id": v["variant_id"]})
                continue
            if len(indexes) == 2 and indexes[0] != indexes[1]:
                heterozygotes.append({"start": v["pos0"], "end": v["pos0"] + len(v["ref"]),
                                      "phase_key": call["phase_key"], "variant_id": v["variant_id"]})
            if any(indexes):
                accepted.append(v)
        conflicts = set()
        for i, a in enumerate(accepted):
            for b in accepted[i+1:]:
                if intersects(a["pos0"], a["pos0"]+len(a["ref"]), b["pos0"], b["pos0"]+len(b["ref"])):
                    conflicts.update([a["variant_id"], b["variant_id"]])
        for v in accepted:
            if v["variant_id"] in conflicts:
                gaps.append({"start": v["pos0"], "end": v["pos0"]+len(v["ref"]),
                             "reason": "overlapping_variants", "variant_id": v["variant_id"]})
        changes = []
        for v in accepted:
            if v["variant_id"] in conflicts:
                continue
            idx = v["calls"][donor_id]["selected_indexes"][copy_index]
            if idx:
                changes.append((v["edits"][idx], v["variant_id"]))
        changes.sort(key=lambda e: (e[0]["start"], e[0]["end"], e[1]))
    sequence, path, deletions = [], [], []

    def append_reference(lo, hi, bases=None):
        bases = reference.fetch(contig, lo, hi) if bases is None and hi > lo else bases or ""
        for pos, base in zip(range(lo, hi), bases):
            if any(g["start"] <= pos < g["end"] for g in gaps):
                base = "N"
            sequence.append(base)
            path.append({"kind": "reference", "pos0": pos, "base": base})

    cursor = start
    for edit, vid in changes:
        append_reference(cursor, edit["start"])
        if edit["kind"] == "substitution":
            append_reference(edit["start"], edit["end"], edit["sequence"])
        elif edit["kind"] == "insertion":
            insertion_id = identity("ins", [reference.sha256, contig, edit["start"], edit["sequence"]])
            for offset, base in enumerate(edit["sequence"]):
                path.append({"kind": "insertion", "boundary0": edit["start"], "offset0": offset,
                             "insertion_sequence_id": insertion_id, "base": base})
                sequence.append(base)
        else:
            deletions.append({"start": edit["start"], "end": edit["end"], "variant_id": vid})
        cursor = edit["end"]
    append_reference(cursor, end)
    sequence = "".join(sequence)
    payload = {"reference_id": reference.sha256, "contig": contig, "interval": interval,
               "sequence_genomic_forward": sequence, "sequence_path": path, "deletion_events": deletions,
               "excluded_intervals": sorted(gaps, key=lambda g: (g["start"], g["end"], g["variant_id"])),
               "heterozygous_sites": heterozygotes,
               **({"phase_constraints": plan["phase_constraints"]} if plan is not None else {})}
    return record("sequence_configuration", run_id, sequence_configuration_id=identity("seq", payload),
                  donor_id=donor_id, local_copy=copy_index, **payload,
                  sequence_status="partial" if gaps else "local_phase_requires_window_check" if len(heterozygotes) > 1 else "resolved",
                  phase_status="reference" if donor_id is None else "local_window_evaluation",
                  reference_evidence="pinned_reference" if donor_id is None else "evaluated_per_target",
                  normalized_variant_keys=sorted(v["variant_id"] for v in accepted), source_call_ids=sorted(calls),
                  variant_intervals=[{"variant_id": v["variant_id"], "start": v["pos0"],
                                      "end": v["pos0"]+len(v["ref"]),
                                      "insertion_boundaries": sorted({v["edits"][i]["start"]
                                          for i in v["calls"][donor_id]["selected_indexes"]
                                          if i and v["edits"][i] and v["edits"][i]["kind"] == "insertion"})}
                                     for v in accepted])


def local_status(material, start, end):
    """A single heterozygote resolves an unordered pair; larger sets require common phase."""
    if start == end:
        start, end = start-1, end+1
    interval = material["interval"]
    if start < interval["start"] or end > interval["end"]:
        return False, "outside_materialized_scope"
    if any(intersects(start, end, g["start"], g["end"]) for g in material["excluded_intervals"]):
        return False, "unresolved_variant_footprint"
    sites = [h for h in material["heterozygous_sites"] if intersects(start, end, h["start"], h["end"])]
    if not sites:
        return True, "homozygous"
    if len(sites) == 1:
        return True, "single_heterozygous_site"
    keys = [h["phase_key"] for h in sites]
    if material.get("phase_constraints"):
        # Only constraints whose two records affect this window participate.
        # They do not erase a phase boundary in a different nearby window.
        ids = {h["variant_id"] for h in sites}
        neighbors = {vid: set() for vid in ids}
        for i, a in enumerate(sites):
            for b in sites[i+1:]:
                if a["phase_key"] is not None and a["phase_key"] == b["phase_key"]:
                    neighbors[a["variant_id"]].add(b["variant_id"])
                    neighbors[b["variant_id"]].add(a["variant_id"])
        used_constraint = False
        for constraint in material["phase_constraints"]:
            a, b = constraint["variant_ids"]
            if a in ids and b in ids:
                neighbors[a].add(b)
                neighbors[b].add(a)
                used_constraint = True
        seen, pending = set(), [next(iter(ids))]
        while pending:
            vid = pending.pop()
            if vid not in seen:
                seen.add(vid)
                pending.extend(neighbors[vid] - seen)
        if used_constraint and seen == ids:
            return True, "resolved_spanning_deletion_phase_constraint"
    if any(k is None for k in keys):
        return False, "unresolved_unphased"
    if any(k != keys[0] for k in keys[1:]):
        return False, "unresolved_phase_blocks"
    return True, "phased_explicit_ps" if keys[0][3] == "explicit" else "phased_implicit_ps"


def build_materializations(config, reference, variants, donors, run_id):
    donors = canonical_donors(donors)
    limits = config["resource_limits"]
    intervals = config["scope"]["materialized_intervals"]
    if len(donors) > limits["max_donors"]:
        raise ResourceLimit("max_donors exceeded")
    if sum(i["end"] - i["start"] for i in intervals) > limits["max_reference_bases"]:
        raise ResourceLimit("max_reference_bases exceeded")
    if len(intervals) * (1 + 2*len(donors)) > limits["max_configurations"]:
        raise ResourceLimit("max_configurations exceeded before reconstruction")
    materials, bases = [], 0
    for interval in intervals:
        for donor in [None] + [d["donor_id"] for d in donors]:
            for copy in range(1 if donor is None else 2):
                m = reconstruct(reference, interval, variants, donor, copy, run_id)
                bases += len(m["sequence_genomic_forward"])
                if bases > limits["max_reconstructed_bases"]:
                    raise ResourceLimit("max_reconstructed_bases exceeded")
                materials.append(m)
    return materials, {"materializations": len(materials), "reconstructed_bases": bases}
