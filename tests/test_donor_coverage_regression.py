"""Independent donor-coverage oracle and fixtures (second regression corpus).

Three-way agreement pattern for every case:
  hand-authored EXPECTED literals (this file) == independent oracle recomputation == pipeline.

The oracle restates the accepted contract semantics from first principles
(interval arithmetic over GT positions); it imports no pipeline code. Contract
rules restated independently:

R1  A complete call's carried supported literal alleles apply per copy.
R2  An incomplete call (missing, partial, filtered, ploidy mismatch) gaps the
    record footprint on every declared copy: the GT itself is not trusted.
R3  A complete call carrying an unsupported literal allele (complex or beyond
    the reconstruction limit) gaps the record footprint on carrying copies only.
R4  A carried spanning-deletion star is resolved when exactly one other record
    has a complete call carrying a pure deletion whose trimmed deleted span
    covers the star record's first base and the per-GT-position deleted/starred
    dosage matches. Resolution pairs the two records' GT positions (with a
    parity flip when the orders differ) and resolves window phase across them.
    Otherwise the star is unresolved; when flanking complete all-present calls
    bracket it within one ploidy region, the possible-deletion envelope between
    the anchors gaps the star-carrying GT position's copy only. Without bounds,
    the event excludes the whole contig on every copy and names its missing
    information (the deletion record or flanking anchors).
R5  Per-copy consistency: an applied deletion may not cover a record where the
    same copy asserts a literal base, and applied edits may not overlap on one
    copy. Violations drop the edits and gap the union footprint on that copy.
    Records in separate explicit phase sets cannot contradict one copy.
R6  Ploidy comes only from explicit intervals: unknown or changing ploidy over
    a cluster means zero materialized copies with an explicit exclusion; zero
    declared copies excludes; GT arity must match the declared ploidy.
R7  Window phase: a 23-base window is evaluable on a copy when it touches no
    gap and its heterozygous sites are zero, one, share one explicit phase key,
    or are connected by star-resolution constraints.

Cases needing the bounded-star machinery are RECOVERED_CASES; PRESERVED_CASES
must hold at every revision.
"""
import hashlib
import json
import re
from pathlib import Path

import pytest

FIXTURES = Path(__file__).resolve().parent / "fixtures" / "donor-coverage-regression"
PAD = 22
WINDOW = 23
POLICY = {"phase": "explicit_ps", "min_gq": None, "min_dp": None,
          "accepted_site_filters": ["PASS", "."], "accepted_sample_filters": ["."]}


def generate(length, seed):
    state = seed
    out = []
    for _ in range(length):
        state = (state * 1103515245 + 12345) % 2**31
        out.append("ACGT"[(state >> 16) % 4])
    return "".join(out)


class FakeReference:
    def __init__(self, *contigs):
        self.sequences = {c: (FIXTURES / "references" / (c + ".fa")).read_text().splitlines()[1].strip()
                          for c in contigs}
        self.build = "artificial"
        self.aliases = {c: c for c in contigs}
        self.sha256 = hashlib.sha256(
            "".join(self.sequences[c] for c in sorted(self.sequences)).encode()).hexdigest()

    def fetch(self, contig, start, end):
        assert 0 <= start < end <= len(self.sequences[contig])
        return self.sequences[contig][start:end]

    def bounds(self, contig):
        return [[0, len(self.sequences[contig])]]


# ---------------------------------------------------------------------------
# Independent oracle
# ---------------------------------------------------------------------------

def parse_vcf(text):
    """Minimal strict VCF parse; returns row dicts in file order."""
    rows = []
    for line in text.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split("\t")
        info = {}
        for item in f[7].split(";"):
            if item and item != ".":
                k, _, v = item.partition("=")
                info[k] = v
        fmt = f[8].split(":")
        sample = dict(zip(fmt, f[9].split(":")))
        gt = sample.get("GT", ".")
        rows.append({"contig": f[0], "pos0": int(f[1]) - 1, "ref": f[3],
                     "alts": f[4].split(","), "filt": f[6], "info": info,
                     "gt": gt, "gt_elems": re.split(r"[/|]", gt),
                     "phased": "|" in gt, "ps": sample.get("PS", "."),
                     "end": int(info["END"]) if "END" in info else int(f[1]) - 1 + len(f[3]),
                     "rid": len(rows)})
    return rows


def trimmed_deletion(pos0, ref, alt):
    """Deleted reference interval [start, end) for a pure deletion, else None."""
    if not alt or set(alt) - set("ACGT") or len(alt) >= len(ref):
        return None
    p = 0
    while p < min(len(ref), len(alt)) and ref[p] == alt[p]:
        p += 1
    r, a = ref[p:], alt[p:]
    s = 0
    while s < len(r) and s < len(a) and r[-1 - s] == a[-1 - s]:
        s += 1
    if a[:len(a) - s or None]:
        return None
    r2 = r[:len(r) - s or None]
    return (pos0 + p, pos0 + p + len(r2))


def allele_edit(pos0, ref, alt):
    """Supported literal edit for an allele, else None (independent restatement)."""
    if not alt or set(alt) - set("ACGT"):
        return None
    if len(ref) == len(alt):
        return {"kind": "substitution", "start": pos0, "end": pos0 + len(ref), "sequence": alt}
    span = trimmed_deletion(pos0, ref, alt)
    if span is not None:
        return {"kind": "deletion", "start": span[0], "end": span[1], "sequence": ""}
    p = 0
    while p < min(len(ref), len(alt)) and ref[p] == alt[p]:
        p += 1
    r, a = ref[p:], alt[p:]
    s = 0
    while s < len(r) and s < len(a) and r[-1 - s] == a[-1 - s]:
        s += 1
    if r[:len(r) - s or None]:
        return None
    return {"kind": "insertion", "start": pos0 + p, "end": pos0 + p,
            "sequence": a[:len(a) - s or None]}


def call_status(row, copies):
    """Independent call model: (status, indexes per GT position)."""
    elems = row["gt_elems"]
    if row["gt"] == "." or all(e == "." for e in elems):
        return "missing", None
    if any(e == "." for e in elems):
        return "partial_missing", None
    if any(not e.isdigit() for e in elems) or ("/" in row["gt"] and "|" in row["gt"]):
        return "invalid", None
    if copies in (1, 2) and len(elems) != copies:
        return "unsupported_ploidy", None
    if row["filt"] not in POLICY["accepted_site_filters"]:
        return "filtered", None
    return "complete", [int(e) for e in elems]


def interval_copies(ploidy, contig, start, end):
    for row in ploidy:
        if row["contig"] == contig and row["start"] <= start and end <= row["end"]:
            return row["copies"]
    return None


def _phase_key(row):
    if row["phased"] and row["ps"] != ".":
        return ("explicit", row["ps"])
    return None


def resolve_stars(rows):
    """R4: per-row resolution, parity flips, envelopes and unbounded markers.

    Returns {rid: {"resolved": (ancestor_rid, flip) | None,
                   "envelope": (start, end) | None, "unbounded": bool}}.
    Anchors are complete calls in a declared ploidy region whose carried
    alleles are all literal and the same length as REF (present on every copy).
    """
    out = {}
    n = len(rows)
    for i, row in enumerate(rows):
        star_indexes = {k + 1 for k, a in enumerate(row["alts"]) if a == "*"}
        if not star_indexes:
            continue
        status, idx = call_status(row, row.get("_region_copies") or 2)
        resolution = None
        if status == "complete":
            carried = []
            for j, other in enumerate(rows):
                if j == i or other["pos0"] >= row["pos0"] or \
                        row["pos0"] >= other["pos0"] + len(other["ref"]):
                    continue
                ostatus, oidx = call_status(other, other.get("_region_copies") or 2)
                if ostatus != "complete":
                    continue
                covering = [k for k in oidx
                            if k and trimmed_deletion(other["pos0"], other["ref"], other["alts"][k - 1])
                            and trimmed_deletion(other["pos0"], other["ref"], other["alts"][k - 1])[0]
                            <= row["pos0"] < trimmed_deletion(other["pos0"], other["ref"], other["alts"][k - 1])[1]]
                if covering:
                    carried.append((j, [k in covering for k in oidx]))
            if len(carried) == 1:
                j, deleted = carried[0]
                starred = [k in star_indexes for k in idx]
                if len(deleted) == len(starred) and sum(deleted) == sum(starred):
                    flip = 0 if all(deleted) else int(deleted[0] != starred[0])
                    resolution = (rows[j]["rid"], flip)
        if resolution is not None:
            out[row["rid"]] = {"resolved": resolution, "envelope": None, "unbounded": False}
            continue
        # unresolved: bound by flanking anchors within the same ploidy region
        left, right = None, None
        for j in range(i - 1, -1, -1):
            if is_anchor(rows[j]) and rows[j]["pos0"] < row["pos0"] \
                    and same_region(rows[j], row):
                left = rows[j]
                break
        for j in range(i + 1, n):
            if is_anchor(rows[j]) and rows[j]["pos0"] >= row["pos0"] + len(row["ref"]) \
                    and same_region(rows[j], row):
                right = rows[j]
                break
        envelope = None
        cluster_span = None
        if left is not None and right is not None and row.get("_region") is not None:
            envelope = (left["pos0"] + len(left["ref"]), right["pos0"])
            cluster_span = (left["rid"], right["pos0"] + len(right["ref"]))
        out[row["rid"]] = {"resolved": None, "envelope": envelope,
                           "unbounded": envelope is None, "cluster_span": cluster_span}
    return out


def same_region(a, b):
    ra, rb = a.get("_region"), b.get("_region")
    if ra is None or rb is None:
        return False
    return (ra["contig"], ra["start"], ra["end"], ra["copies"]) == \
        (rb["contig"], rb["start"], rb["end"], rb["copies"])


def is_anchor(row):
    status, idx = call_status(row, row.get("_region_copies") or 2)
    if status != "complete" or row.get("_region") is None:
        return False
    return all(len((row["alts"][i - 1] if i else row["ref"])) == len(row["ref"])
               and not set((row["alts"][i - 1] if i else row["ref"])) - set("ACGT")
               for i in idx)


def oracle_model(seq, contig, rows, ploidy):
    """Return clusters with per-copy materials, gaps, het sites and exclusions."""
    length = len(seq)
    rows = sorted(rows, key=lambda r: (r["pos0"], r["ref"], tuple(r["alts"])))
    for row in rows:
        region = next((r for r in ploidy
                       if r["contig"] == contig and r["start"] <= row["pos0"]
                       and row["pos0"] + len(row["ref"]) <= r["end"]), None)
        row["_region"] = region
        row["_region_copies"] = region["copies"] if region else None
    stars = resolve_stars(rows)
    # Cluster with envelope-extended spans: a bounded star's span covers its
    # envelope, and its left anchor's span extends through the right anchor, so
    # the neighborhood stays in one reconstruction.
    clusters = []
    anchor_ext = {}
    for row in rows:
        cspan = stars.get(row["rid"], {}).get("cluster_span")
        if cspan:
            anchor_ext[cspan[0]] = max(anchor_ext.get(cspan[0], 0), cspan[1])
    for row in rows:
        info = stars.get(row["rid"], {})
        env = info.get("envelope")
        lo = max(0, min(row["pos0"], env[0] if env else row["pos0"]) - PAD)
        hi = min(length, max(row["end"], env[1] if env else row["end"],
                             anchor_ext.get(row["rid"], 0) or row["end"]) + PAD)
        if clusters and lo <= clusters[-1]["end"]:
            clusters[-1]["end"] = max(clusters[-1]["end"], hi)
            clusters[-1]["rows"].append(row)
        else:
            clusters.append({"contig": contig, "start": lo, "end": hi, "rows": [row]})
    results = []
    for cluster in clusters:
        copies = interval_copies(ploidy, contig, cluster["start"], cluster["end"])
        out = {"interval": {"contig": contig, "start": cluster["start"], "end": cluster["end"]},
               "copies": copies if copies is not None else 0, "materials": {}, "exclusions": [],
               "unbounded_exclusions": [], "masked_copies": set()}
        if copies is None:
            out["exclusions"].append({"reason": "unknown_or_changing_ploidy", **out["interval"]})
            results.append(out)
            continue
        if copies == 0:
            out["exclusions"].append({"reason": "zero_ploidy", **out["interval"]})
            results.append(out)
            continue
        calls = {row["rid"]: call_status(row, copies) for row in cluster["rows"]}
        parity = {}
        constraints = []
        for row in cluster["rows"]:
            info = stars.get(row["rid"])
            if info and info["resolved"] is not None:
                anc, flip = info["resolved"]
                parity[row["rid"]] = flip
                constraints.append((anc, row["rid"]))
        # Unbounded events: whole-cluster gaps on every copy plus contig exclusion.
        unbounded_rows = [row for row in cluster["rows"]
                          if stars.get(row["rid"], {}).get("unbounded")]
        for copy in range(copies):
            edits, gaps, het, present = [], [], [], []
            for row in cluster["rows"]:
                status, idx = calls[row["rid"]]
                if status != "complete":
                    gaps.append((row["pos0"], row["end"], "call_" + status))  # R2
                    continue
                i = idx[copy ^ parity.get(row["rid"], 0)]
                info = stars.get(row["rid"])
                if len(idx) == 2 and idx[0] != idx[1]:
                    env = info["envelope"] if info and info["envelope"] is not None else None
                    het.append((*(env or (row["pos0"], row["pos0"] + len(row["ref"]))),
                                _phase_key(row), row["rid"]))
                if i and row["alts"][i - 1] == "*":
                    if info and info["resolved"] is not None:
                        continue  # resolved star: no edit, no gap
                    if info and info["envelope"] is not None:
                        gaps.append((*info["envelope"], "unresolved_spanning_deletion"))
                    continue
                if i == 0:
                    present.append((row["pos0"], row["pos0"] + len(row["ref"]), row["rid"]))
                else:
                    edit = allele_edit(row["pos0"], row["ref"], row["alts"][i - 1])
                    if edit is None:
                        gaps.append((row["pos0"], row["pos0"] + len(row["ref"]),
                                     "unsupported_carried_allele"))  # R3
                        continue
                    edit["ref_span"] = (row["pos0"], row["pos0"] + len(row["ref"]))
                    edit["rid"] = row["rid"]
                    edits.append(edit)
                    if edit["kind"] == "deletion":
                        present.extend([(row["pos0"], edit["start"], row["rid"]),
                                        (edit["end"], row["pos0"] + len(row["ref"]), row["rid"])])
                    else:
                        present.append((edit["start"], edit["end"], row["rid"]))
            # R5: per-copy conflicts (records in separate explicit phase sets cannot
            # contradict one copy).
            dropped, conflicted = set(), set()
            for a_i, a in enumerate(edits):
                for b in edits[a_i + 1:]:
                    if a["start"] < b["ref_span"][1] and b["start"] < a["ref_span"][1] and \
                            intersects(a["start"], a["end"], b["start"], b["end"]):
                        ka = _phase_key(cluster_row(cluster, a["rid"]))
                        kb = _phase_key(cluster_row(cluster, b["rid"]))
                        if ka is None or kb is None or ka == kb:
                            lo = min(a["ref_span"][0], b["ref_span"][0])
                            hi = max(a["ref_span"][1], b["ref_span"][1])
                            gaps.append((lo, hi, "overlapping_variants"))
                            dropped.update([a["rid"], b["rid"]])
                            conflicted.update([a["rid"], b["rid"]])
                if a["kind"] == "deletion":
                    for lo, hi, rid in present:
                        if rid == a["rid"] or not intersects(a["start"], a["end"], lo, hi):
                            continue
                        ka = _phase_key(cluster_row(cluster, a["rid"]))
                        kb = _phase_key(cluster_row(cluster, rid))
                        if ka is None or kb is None or ka == kb:
                            covered = cluster_row(cluster, rid)
                            gaps.append((min(a["ref_span"][0], covered["pos0"]),
                                         max(a["ref_span"][1], covered["pos0"] + len(covered["ref"])),
                                         "overlapping_variants"))
                            dropped.add(a["rid"])
                            conflicted.update([a["rid"], rid])
            edits = [e for e in edits if e["rid"] not in dropped]
            if unbounded_rows:
                gaps.append((cluster["start"], cluster["end"], "unresolved_spanning_deletion"))
            sequence, path = [], []
            cursor = cluster["start"]
            for e in sorted(edits, key=lambda x: (x["start"], x["end"], x["rid"])):
                for pos in range(cursor, e["start"]):
                    sequence.append(seq[pos])
                    path.append(pos)
                if e["kind"] == "substitution":
                    for k, base in enumerate(e["sequence"]):
                        sequence.append(base)
                        path.append(e["start"] + k)
                elif e["kind"] == "insertion":
                    sequence.extend(e["sequence"])
                    path.extend([None] * len(e["sequence"]))
                cursor = e["end"]
            for pos in range(cursor, cluster["end"]):
                sequence.append(seq[pos])
                path.append(pos)
            text = "".join(sequence)
            for gs, ge, _reason in gaps:
                clipped = (max(gs, cluster["start"]), min(ge, cluster["end"]))
                for pos in range(*clipped):
                    if pos in path:
                        text = text[:path.index(pos)] + "N" + text[path.index(pos) + 1:]
            out["materials"][copy] = {"sequence": text, "gaps": sorted(set(gaps)),
                                      "het": het, "path": path, "constraints": constraints,
                                      "rows": cluster["rows"]}
        if unbounded_rows:
            out["unbounded_exclusions"].append({
                "contig": contig, "start": 0, "end": length, "reason": "unbounded_unsupported_allele",
                "missing_information": "flanking accepted called bases or a compatible covering deletion"})
        results.append(out)
    return results


def cluster_row(cluster, rid):
    return next(r for r in cluster["rows"] if r["rid"] == rid)


def intersects(a, b, c, d):
    return a < d and c < b


def oracle_window_status(material, start, end):
    """R7: evaluable flag and reason for one 23-base window on one copy."""
    if start < material["interval"][0] or end > material["interval"][1]:
        return False, "outside_materialized_scope"
    for gs, ge, _ in material["gaps"]:
        if start < ge and gs < end:
            return False, "gap"
    sites = [h for h in material["het"] if start < h[1] and h[0] < end]
    if not sites:
        return True, "homozygous"
    if len(sites) == 1:
        return True, "single_heterozygous_site"
    keys = [h[2] for h in sites]
    if all(k is not None and k == keys[0] for k in keys):
        return True, "phased_explicit_ps"
    # Star-resolution constraints connect the two involved records' sites.
    if material["constraints"]:
        neighbors = {h[3]: set() for h in sites}
        for i, a in enumerate(sites):
            for b in sites[i + 1:]:
                if a[2] is not None and a[2] == b[2]:
                    neighbors[a[3]].add(b[3])
                    neighbors[b[3]].add(a[3])
        used = False
        for a_rid, b_rid in material["constraints"]:
            if a_rid in neighbors and b_rid in neighbors:
                neighbors[a_rid].add(b_rid)
                neighbors[b_rid].add(a_rid)
                used = True
        seen, pending = set(), [sites[0][3]]
        while pending:
            cur = pending.pop()
            if cur not in seen:
                seen.add(cur)
                pending.extend(neighbors[cur] - seen)
        if used and len(seen) == len(sites):
            return True, "resolved_by_spanning_deletion_constraint"
    return False, "unresolved_phase"


# ---------------------------------------------------------------------------
# Hand-authored expectations
# ---------------------------------------------------------------------------
# Derived by hand from the fixture references; the oracle must reproduce them
# independently. "clusters" = (start, end, copies); "sequences" keyed by copy
# (single-cluster cases) or (cluster, copy), "len:N" asserts length only;
# "gaps" per copy as (start, end, reason); "het" = (start, end, key|None)
# shared by every copy of the cluster; "probes" = (window_start, evaluable)
# per copy or shared; "probe_bases" = {(copy, pos): letter | "N" | None(ref)};
# "probe_absent" = positions absent from a copy's coordinate path;
# "masked_positions" = (start, end) spans that must be N-masked.

EXPECTED = {
    "c01_star_explained_same_position": {
        "behavior": "recovered",
        "clusters": [(78, 133, 2)],
        "sequences": {
            0: "TTCCGTGTGGGACGCGAAGGTTCAGGAAACTGCCATTTTCTCAGT",
            1: "TTCCGTGTGGGACGCGAAGGTTCTTCATATACGAGGAAACTGCCATTTTCTCAGT"},
        "gaps": {0: [], 1: []},
        "het": [(100, 111, ("explicit", "1441")), (105, 106, None)],
        "probes": [(79, True), (88, True), (95, True), (110, True)],
        "probe_absent": {0: [101, 105, 109], 1: []},
        "probe_bases": {(1, 105): "T"},
    },
    "c02_star_explained_large_deletion": {
        "behavior": "recovered",
        "clusters": [(178, 343, 2)],
        "sequences": {0: "len:45", 1: "len:165"},
        "gaps": {0: [], 1: []},
        "het": [(200, 321, ("explicit", "2442")), (250, 251, None)],
        "probes": {0: [(178, True), (240, True), (320, True)],
                   1: [(178, True), (240, True), (320, True)]},
        "probe_bases": {(0, 199): None, (1, 250): "A", (1, 320): None},
        "probe_absent": {0: [201, 250, 320]},
    },
    "c03_star_orphan_single": {
        "behavior": "recovered",
        "clusters": [(328, 473, 2)],
        "sequences": {0: "len:145", 1: "len:145"},
        "gaps": {0: [], 1: [(351, 450, "unresolved_spanning_deletion")]},
        "het": [(351, 450, None)],
        "probes": {0: [(328, True), (400, True), (449, True), (450, True)],
                   1: [(328, True), (329, False), (449, False), (450, True)]},
        "probe_bases": {(0, 400): "G", (1, 350): None, (1, 449): "N", (1, 450): None},
        "masked_positions": {1: [(351, 450)]},
        "exclusion_reasons": [],
    },
    "c04_star_orphan_cluster": {
        "behavior": "recovered",
        "clusters": [(428, 583, 2)],
        "sequences": {0: "len:155", 1: "len:155"},
        "gaps": {0: [(451, 560, "unresolved_spanning_deletion")], 1: []},
        "het": [(451, 560, None), (451, 560, None), (451, 560, None)],
        "probes": {0: [(428, True), (429, False), (559, False), (560, True)],
                   1: [(428, True), (429, False), (559, False), (560, True)]},
        "probe_bases": {(0, 449): None, (0, 500): "N", (0, 560): None,
                        (1, 500): "T", (1, 508): "G", (1, 517): "G"},
        "masked_positions": {0: [(451, 560)]},
        "exclusion_reasons": [],
    },
    "c05_star_multiallelic_same_record": {
        "behavior": "recovered",
        "clusters": [(528, 683, 2)],
        "sequences": {0: "len:148", 1: "len:155"},
        "gaps": {0: [], 1: [(551, 660, "unresolved_spanning_deletion")]},
        "het": [(551, 660, None)],
        "probes": {0: [(528, True), (551, True), (659, True)], 1: [(528, True), (551, False), (660, True)]},
        "probe_bases": {(0, 600): None, (1, 600): "N"},
        "probe_absent": {0: [601, 607]},
        "masked_positions": {1: [(551, 660)]},
        "exclusion_reasons": [],
    },
    "c06_star_filtered_support": {
        "behavior": "recovered",
        "clusters": [(628, 773, 2)],
        "sequences": {0: "len:145", 1: "len:145"},
        "gaps": {0: [(651, 750, "unresolved_spanning_deletion"), (700, 711, "call_filtered")],
                 1: [(700, 711, "call_filtered")]},
        "het": [(651, 750, None)],
        "probes": {0: [(628, True), (629, False), (749, False), (750, True)],
                   1: [(628, True), (678, False), (711, True), (750, True)]},
        "probe_bases": {(0, 650): None, (0, 660): "N", (1, 699): None, (1, 705): "N", (1, 711): None},
        "masked_positions": {0: [(651, 750)], 1: [(700, 711)]},
        "exclusion_reasons": [],
    },
    "c07_overlapping_same_copies_conflict": {
        "behavior": "preserved",
        "clusters": [(779, 843, 2)],
        "sequences": {0: "len:64", 1: "len:64"},
        "gaps": {0: [(801, 821, "overlapping_variants")], 1: [(801, 821, "overlapping_variants")]},
        "het": [],
        "masked_positions": {0: [(801, 821)], 1: [(801, 821)]},
        "probe_bases": {(0, 800): None, (0, 810): "N", (1, 810): "N"},
    },
    "c08a_overlapping_split_copies_starred": {
        "behavior": "recovered",
        "clusters": [(879, 943, 2)],
        "sequences": {0: "len:45", 1: "len:64"},
        "gaps": {0: [], 1: []},
        "het": [(901, 921, ("explicit", "88")), (910, 911, None)],
        "probes": [(888, True), (895, True)],
        "probe_bases": {(1, 910): "C", (0, 921): None},
        "probe_absent": {0: [902, 910, 920]},
    },
    "c08b_overlapping_split_copies_contradictory": {
        "behavior": "recovered",
        "clusters": [(879, 943, 2)],
        "sequences": {0: "len:64", 1: "len:64"},
        "gaps": {0: [(901, 921, "overlapping_variants")], 1: []},
        "het": [(901, 921, ("explicit", "88")), (910, 911, ("explicit", "88"))],
        "probes": {0: [(895, False)], 1: [(895, True)]},
        "probe_bases": {(0, 910): "N", (1, 910): "C"},
        "masked_positions": {0: [(901, 921)]},
    },
    "c09_missing_call": {
        "behavior": "preserved",
        "clusters": [(978, 1023, 2)],
        "sequences": {0: "len:45", 1: "len:45"},
        "gaps": {0: [(1000, 1001, "call_missing")], 1: [(1000, 1001, "call_missing")]},
        "het": [],
        "probes": [(978, False), (999, False)],
        "probe_bases": {(0, 999): None, (0, 1000): "N", (1, 1000): "N"},
    },
    "c10_partial_missing_call": {
        "behavior": "preserved",
        "clusters": [(1028, 1073, 2)],
        "sequences": {0: "len:45", 1: "len:45"},
        "gaps": {0: [(1050, 1051, "call_partial_missing")],
                 1: [(1050, 1051, "call_partial_missing")]},
        "het": [],
        "probes": [(1028, False)],
    },
    "c11_phase_block_boundary": {
        "behavior": "preserved",
        "clusters": [(1078, 1139, 2)],
        "sequences": {
            0: "TGCGATAACGAGGCCATATGAGTTTTACAGGTGTCCGGAATGACAATGCCGGCGGGAGTAA",
            1: "TGCGATAACGAGGCCATATGAGATTTACAGCTGTCCGGCATGACAATGCCGGCGGGAGTAA"},
        "gaps": {0: [], 1: []},
        "het": [(1100, 1101, ("explicit", "31")), (1108, 1109, ("explicit", "31")),
                (1116, 1117, ("explicit", "47"))],
        "probes": [(1085, True), (1093, True), (1094, False), (1110, True)],
    },
    "c12_supported_indel_paths": {
        "behavior": "preserved",
        "clusters": [(128, 181, 2)],
        "sequences": {
            0: "TCAGTCTCTCGCGCTCGGGTCCAGGAGAGCGCAGAAGTAAGGACTAACCACCCCT",
            1: "TCAGTCTCTCGCGCTCGGGTCCAAGAGCGCCGAAGTAAGGACTAACCACCCCT"},
        "gaps": {0: [], 1: []},
        "het": [(150, 151, ("explicit", "52")), (158, 159, ("explicit", "52"))],
        "probes": [(135, True), (136, True)],
        "insertion_bases": {0: ("GG", 23)},
    },
    "c13_complex_allele_bounded": {
        "behavior": "recovered",
        "clusters": [(328, 376, 2)],
        "sequences": {0: "len:48", 1: "len:48"},
        "gaps": {0: [], 1: [(350, 354, "unsupported_carried_allele")]},
        "het": [(350, 354, None)],
        "probes": {0: [(328, True), (353, True)], 1: [(328, False), (353, False)]},
        "probe_bases": {(0, 351): None, (1, 351): "N"},
        "masked_positions": {1: [(350, 354)]},
    },
    "c14_star_orphan_unbounded": {
        "behavior": "recovered",
        "clusters": [(678, 723, 2)],
        "sequences": {0: "len:45", 1: "len:45"},
        "gaps": {0: [(678, 723, "unresolved_spanning_deletion")],
                 1: [(678, 723, "unresolved_spanning_deletion")]},
        "het": [(700, 701, None)],
        "probes": [(678, False), (700, False)],
        "probe_bases": {(0, 700): "N", (1, 700): "N"},
        "masked_positions": {0: [(678, 723)], 1: [(678, 723)]},
        "exclusion_reasons": ["unbounded_unsupported_allele"],
    },
    "x01_xy_region_specific_ploidy": {
        "behavior": "preserved",
        "contig": "chrX",
        "clusters": [(178, 223, 2), (578, 623, 1)],
        "sequences": {(0, 0): "ACTCGTTTTGCGAATTACGCTGCGTTTTATTGCAGTAACAGGACG",
                      (0, 1): "ACTCGTTTTGCGAATTACGCTGGGTTTTATTGCAGTAACAGGACG",
                      (1, 0): "GAAACAAATGTTATTGCGAGGTGACCTCTGTCCTTAGATTACCAG"},
        "gaps": {(0, 0): [], (0, 1): [], (1, 0): []},
        "het": {(0, 0): [(200, 201, None)], (0, 1): [(200, 201, None)], (1, 0): []},
    },
    "x02_xy_diploid_call_in_haploid_region": {
        "behavior": "preserved",
        "contig": "chrX",
        "clusters": [(578, 623, 1)],
        "sequences": {(0, 0): "len:45"},
        "gaps": {(0, 0): [(600, 601, "call_unsupported_ploidy")]},
        "probe_bases": {((0, 0), 600): "N"},
    },
    "x03_xy_undeclared_ploidy": {
        "behavior": "preserved",
        "contig": "chrY",
        "clusters": [(78, 123, 0)],
        "sequences": {},
        "exclusion_reasons": ["unknown_or_changing_ploidy"],
    },
    "x04_xy_zero_copies": {
        "behavior": "preserved",
        "contig": "chrY",
        "clusters": [(78, 123, 0)],
        "sequences": {},
        "exclusion_reasons": ["zero_ploidy"],
    },
    "x05_xy_haploid_star_explained": {
        "behavior": "recovered",
        "contig": "chrX",
        "clusters": [(278, 331, 1)],
        "sequences": {0: "GGAGTTACATCCCATACTTCGTACGAGAGACGTTCACAAACCAAC"},
        "gaps": {0: []},
        "het": [],
        "probe_absent": {0: [301, 305, 308]},
    },
    "x06_xy_par_star_explained": {
        "behavior": "recovered",
        "contig": "chrX",
        "clusters": [(188, 241, 2)],
        "sequences": {0: "CGAATTACGCTGCGTTTTATTGCGACGCCTGCGCCTGATGTTCAT",
                      1: "CGAATTACGCTGCGTTTTATTGCAGTAACAGGACGCCTGCGCCTGATGTTCAT"},
        "gaps": {0: [], 1: []},
        "het": [(210, 219, ("explicit", "99")), (214, 215, ("explicit", "99"))],
        "probes": [(196, True), (200, True)],
        "probe_absent": {0: [211, 214, 218]},
    },
}

RECOVERED_CASES = sorted(k for k, v in EXPECTED.items() if v["behavior"] == "recovered")
PRESERVED_CASES = sorted(k for k, v in EXPECTED.items() if v["behavior"] == "preserved")


def load_case(name):
    manifest = json.loads((FIXTURES / "manifest.json").read_text())
    contig = EXPECTED[name].get("contig", "artificial_1")
    text = (FIXTURES / "cases" / (name + ".vcf")).read_text()
    ploidy = manifest["cases"][name]["ploidy"]
    return contig, text, ploidy


def _per_copy(mapping, n_clusters):
    out = {}
    for key, value in mapping.items():
        if isinstance(key, tuple):
            out[key] = value
        else:
            for ci in range(n_clusters):
                out[(ci, key)] = value
    return out


def _gaps_by_reason(gaps):
    """Merged masked positions per reason: representation-independent gap content."""
    merged = {}
    for gs, ge, reason in gaps:
        spans = merged.setdefault(reason, [])
        placed = False
        for span in spans:
            if span[0] <= ge and gs <= span[1]:
                span[0], span[1] = min(span[0], gs), max(span[1], ge)
                placed = True
                break
        if not placed:
            spans.append([gs, ge])
        changed = True
        while changed:
            changed = False
            for a in list(spans):
                for b in list(spans):
                    if a is not b and a[0] <= b[1] and b[0] <= a[1]:
                        a[0], a[1] = min(a[0], b[0]), max(a[1], b[1])
                        spans.remove(b)
                        changed = True
                        break
                if changed:
                    break
    return {r: sorted(tuple(s) for s in spans) for r, spans in merged.items()}


def _check_model(name, models):
    """Shared assertions: hand-authored EXPECTED vs a case model (oracle or pipeline)."""
    exp = EXPECTED[name]
    assert [(m["interval"]["start"], m["interval"]["end"], m["copies"]) for m in models] == \
        [tuple(c) for c in exp["clusters"]], name
    sequences = _per_copy(exp.get("sequences", {}), len(models))
    gaps = _per_copy(exp.get("gaps", {}), len(models))
    hets = exp.get("het", {})
    hets = _per_copy(hets, len(models)) if isinstance(hets, dict) else {
        (ci, h): hets for ci in range(len(models)) for h in models[ci]["materials"]}
    masked_positions = _per_copy(exp.get("masked_positions", {}), len(models))
    for ci, model in enumerate(models):
        assert set(model["materials"]) == set(range(model["copies"])) - set(model.get("masked_copies", ())), name
        for h, material in model["materials"].items():
            want = sequences.get((ci, h))
            if isinstance(want, str) and want.startswith("len:"):
                assert len(material["sequence"]) == int(want[4:]), (name, ci, h)
            elif want is not None:
                assert material["sequence"] == want, (name, ci, h)
            if (ci, h) in gaps:
                assert _gaps_by_reason(material["gaps"]) == _gaps_by_reason(
                    [tuple(g) for g in gaps[(ci, h)]]), (name, ci, h)
            if (ci, h) in hets:
                assert sorted(tuple(x[:3]) for x in material["het"]) == \
                    sorted(tuple(x) for x in hets[(ci, h)]), (name, ci, h)
            for gs, ge in masked_positions.get((ci, h), []):
                for pos in range(gs, ge):
                    assert material["base_at"](pos) == "N", (name, ci, h, pos)
        for key, letter in exp.get("probe_bases", {}).items():
            (ci, h), pos = key if isinstance(key[0], tuple) else ((0, key[0]), key[1])
            material = models[ci]["materials"].get(h)
            if material is None:
                raise AssertionError((name, "probe on missing copy", key))
            want = material["ref_base"](pos) if letter is None else letter
            assert material["base_at"](pos) == want, (name, key, want)
        for key, positions in exp.get("probe_absent", {}).items():
            ci, h = key if isinstance(key, tuple) else (0, key)
            for pos in positions:
                assert not models[ci]["materials"][h]["has_position"](pos), (name, key, pos)
        probes = exp.get("probes", [])
        probe_map = probes if isinstance(probes, dict) else {h: probes for h in model["materials"]}
        for h, plist in probe_map.items():
            if h not in model["materials"]:
                continue
            for start, ok in plist:
                assert models[ci]["materials"][h]["window_ok"](start, start + WINDOW)[0] == ok, \
                    (name, ci, h, start)
        if "insertion_bases" in exp:
            for h, (bases, offset) in exp["insertion_bases"].items():
                assert model["materials"][h]["sequence"][offset:offset + len(bases)] == bases, name
        assert sorted(e["reason"] for e in model.get("exclusions", [])) == \
            sorted(exp.get("exclusion_reasons", [])), name


# ---------------------------------------------------------------------------
# Adapters: one uniform model interface over oracle output and pipeline output
# ---------------------------------------------------------------------------

def _oracle_models(name):
    contig, text, ploidy = load_case(name)
    ref = FakeReference(contig)
    seq = ref.sequences[contig]
    models = []
    for cluster in oracle_model(seq, contig, parse_vcf(text), ploidy):
        materials = {}
        for h, m in cluster["materials"].items():
            path, gaps, het, constraints = m["path"], m["gaps"], m["het"], m["constraints"]
            interval = (cluster["interval"]["start"], cluster["interval"]["end"])
            payload = {"gaps": gaps, "het": het, "constraints": constraints,
                       "rows": m["rows"], "interval": interval}
            materials[h] = {
                "sequence": m["sequence"], "gaps": gaps, "het": het,
                "base_at": lambda pos, p=path, t=m["sequence"]: t[p.index(pos)] if pos in p else None,
                "ref_base": lambda pos, s=seq: s[pos],
                "has_position": lambda pos, p=path: pos in p,
                "window_ok": lambda s, e, pl=payload: oracle_window_status(pl, s, e),
            }
        models.append({"interval": cluster["interval"], "copies": cluster["copies"],
                       "materials": materials,
                       "exclusions": cluster["exclusions"] + cluster["unbounded_exclusions"],
                       "masked_copies": set()})
    return models


PIPELINE_REASON_MAP = {"missing": "call_missing", "partial_missing": "call_partial_missing",
                       "filtered": "call_filtered", "unsupported_ploidy": "call_unsupported_ploidy",
                       "overlapping_carried_edits": "overlapping_variants"}


def _pipeline_models(name, tmp_path):
    from base_edit_safety.population.haplotypes import local_status
    from base_edit_safety.population.research_genotypes import iter_variant_clusters, materialize_cluster
    contig, text, ploidy = load_case(name)
    ref = FakeReference(contig)
    vcf = tmp_path / (name + ".vcf")
    vcf.write_text(text)
    clusters = iter_variant_clusters(vcf, ref, "sample", POLICY, ploidy)
    models = []
    for cluster in clusters:
        interval = cluster["interval"]
        materials = {}
        for m in materialize_cluster(ref, cluster, "sample", "run"):
            path = [p.get("pos0") for p in m["sequence_path"]]
            seq = m["sequence_genomic_forward"]
            het = []
            for site in m["heterozygous_sites"]:
                key = site["phase_key"]
                het.append((site["start"], site["end"],
                            ("explicit", key[4]) if key and key[3] == "explicit" else None))
            gaps = [(g["start"], g["end"], PIPELINE_REASON_MAP.get(g["reason"], g["reason"]))
                    for g in m["excluded_intervals"]]
            materials[m["local_copy"]] = {
                "sequence": seq, "gaps": sorted(gaps), "het": sorted(het),
                "base_at": lambda pos, p=path, t=seq: t[p.index(pos)] if pos in p else None,
                "ref_base": lambda pos, c=contig, r=ref: r.fetch(c, pos, pos + 1),
                "has_position": lambda pos, p=path: pos in p,
                "window_ok": lambda s, e, mm=m: local_status(mm, s, e),
            }
        models.append({"interval": interval, "copies": cluster["copies"],
                       "materials": materials,
                       "exclusions": cluster["exclusions"] + cluster.get("unbounded_exclusions", []),
                       "masked_copies": set(range(cluster["copies"])) - set(materials)})
    return models


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_fixture_regeneration_and_uniqueness():
    for name, seed, length in (("artificial_1", 424242, 1200), ("chrX", 777, 800), ("chrY", 888, 400)):
        seq = generate(length, seed)
        on_disk = (FIXTURES / "references" / (name + ".fa")).read_text().splitlines()[1].strip()
        assert seq == on_disk, name
        forward = {seq[i:i + WINDOW] for i in range(len(seq) - WINDOW + 1)}
        assert len(forward) == len(seq) - WINDOW + 1, name + " repeated forward 23-mers"
        comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
        rc = lambda s: "".join(comp[b] for b in reversed(s))
        assert not {rc(seq[i:i + WINDOW]) for i in range(len(seq) - WINDOW + 1)} & forward, name


@pytest.mark.parametrize("name", sorted(EXPECTED))
def test_fixture_ref_alleles_match_reference(name):
    contig, text, _ploidy = load_case(name)
    ref = FakeReference(contig)
    for row in parse_vcf(text):
        assert ref.fetch(contig, row["pos0"], row["pos0"] + len(row["ref"])) == row["ref"], name


@pytest.mark.parametrize("name", sorted(EXPECTED))
def test_oracle_matches_hand_authored_expectations(name):
    _check_model(name, _oracle_models(name))


@pytest.mark.parametrize("name", PRESERVED_CASES)
def test_pipeline_matches_oracle_preserved(name, tmp_path):
    _check_model(name, _pipeline_models(name, tmp_path))


@pytest.mark.parametrize("name", RECOVERED_CASES)
def test_pipeline_matches_oracle_recovered(name, tmp_path):
    """Recovered behavior; runs against the integrated revision."""
    _check_model(name, _pipeline_models(name, tmp_path))
