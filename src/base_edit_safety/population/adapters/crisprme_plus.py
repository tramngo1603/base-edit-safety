"""Read the retained CRISPRitz rows emitted by pinned CRISPRme+ 2.5.5.

Native contigs contain already reconstructed donor copies. Their names have no
biological meaning; the retained input map supplies the original paths.
"""
from __future__ import annotations

from pathlib import Path
import re
from copy import deepcopy

from ..common import PopulationError, checksum, integer
from ..identity import IUPAC, matches_pam, reverse_complement


VERSION = "v2.5.5"
SOURCE_REVISION = "c005098e05ebbf83f049d9de5fd2598ceed65685"
IMAGE = "pinellolab/crisprme@sha256:3d5fc92e7c96cddef37e2d550e96966ac3d378966d93f2a1e9ea602258db15a6"
EXECUTABLE_HASHES = {
    "/opt/conda/bin/crisprme.py": "53ed8096b2f6af510641bf41fa7dc3e5b9662e8cd5f2208c1004456e1ff19fad",
    "/opt/conda/bin/crispritz.py": "14101f079a6e271f88530ee9f71496902fa24388df58da9552f5579dabb3f9ef",
    "/opt/conda/opt/crispritz/searchTST": "9939777e2c853d915731c13d242e5fb99c85f277c919fa2645638f6dea8a6374",
    "/opt/conda/opt/crispritz/buildTST": "138b0bdf84e658c9d5005ec72b4dbab71feb26ad4f9fdbf2f82e8a3dd4fe0a95",
    "/opt/conda/opt/crispritz/searchBruteForce": "72d660f04effd821361086a269ed0c4e3d5ba3cfa913a1782384b181fcc323c5",
}
HEADER = ("#Bulge type", "crRNA", "DNA", "Chromosome", "Position",
          "Cluster Position", "Direction", "Mismatches", "Bulge Size", "Total")
LEADING_DNA_TRANSFORM = "exact-leading-dna-extension-v1"
WINDOW_GUARD_TRANSFORM = "excluded-native-window-guards-v1"


def native_sequence(window):
    guard = window.get("native_guard", {})
    return guard.get("prefix", "") + window["sequence"] + guard.get("suffix", "")


def validate_policy(design, editor, profile):
    """Reject policies whose native counting or geometry has not been checked."""
    guide = design["guide_sequence"]
    geometry = editor["guide_geometry"]
    if not re.fullmatch("[ACGT]+", guide) or len(guide) != geometry["spacer_length"]:
        raise PopulationError("Native search requires an exact ACGT guide of the declared length")
    if geometry["spacer_length"] != 20 or geometry["pam_side"] != "3prime" \
            or geometry["pam_length"] not in (2, 3):
        raise PopulationError("Pinned native adapter currently supports 20-base guides with 2/3-base 3-prime PAMs")
    patterns = profile["pam_patterns"]
    if not patterns or any(len(p) != geometry["pam_length"] or set(p) - IUPAC.keys() for p in patterns):
        raise PopulationError("Native profile PAMs do not match the declared geometry")
    if profile["total_edit_accounting"] != "mismatch_plus_bulged_bases":
        raise PopulationError("CRISPRme+ counts bulged bases, not bulge events")
    for key in ("max_mismatches", "max_dna_bulges", "max_rna_bulges", "max_total_edits"):
        integer(profile[key], key)
    integer(profile["max_rows"], "max_rows", 1)
    integer(profile["max_window_variants"], "max_window_variants")
    if profile["retention"] != "all" or not profile["strands"] \
            or set(profile["strands"]) - {"+", "-"}:
        raise PopulationError("Native search requires all-row retention and explicit strands")
    if profile.get("backend_name") != "crisprme_plus" \
            or profile.get("backend_version") not in (VERSION, "2.5.5") \
            or profile.get("backend_digest") not in (IMAGE, IMAGE.split("@", 1)[1], IMAGE.rsplit(":", 1)[1]):
        raise PopulationError("Incompatible backend name, version or image digest; use the pinned crisprme_plus v2.5.5 profile")


def _number(value, name):
    if not re.fullmatch(r"0|[1-9][0-9]*", value):
        raise PopulationError(f"Invalid native {name}: {value!r}")
    return int(value)


def parse_row(line, *, windows, guide, geometry, profile, source_ref):
    """Validate one native alignment against the exact submitted sequence."""
    values = line.split("\t")
    if len(values) != len(HEADER):
        raise PopulationError("Native target row must have exactly ten columns (CRISPRme+ v2.5.5)")
    kind, aligned_guide, aligned_target, contig, position, cluster, strand, mm, bulges, total = values
    if kind not in {"X", "DNA", "RNA", "RNA,DNA"} or strand not in {"+", "-"}:
        raise PopulationError("Unrecognized native alignment type or strand")
    if contig not in windows:
        raise PopulationError("Native target contig is absent from the immutable input map")
    window = windows[contig]
    aligned_guide, aligned_target = aligned_guide.upper(), aligned_target.upper()
    if len(aligned_guide) != len(aligned_target) or not aligned_guide \
            or set(aligned_guide) - set("ACGTN-") or set(aligned_target) - set("ACGT-"):
        raise PopulationError("Native guide and target alignment columns are invalid")
    pam_length = geometry["pam_length"]
    if geometry["pam_side"] != "3prime":
        raise PopulationError("Unvalidated native PAM orientation")
    ga, ta = aligned_guide[:-pam_length], aligned_target[:-pam_length]
    pam = aligned_target[-pam_length:]
    if aligned_guide[-pam_length:] != "N" * pam_length or "-" in pam \
            or ga.replace("-", "") != guide or "N" in ga:
        raise PopulationError("Native row does not contain the submitted guide and PAM geometry")
    mismatches, dna, rna = [], [], []
    for index, (g, t) in enumerate(zip(ga, ta)):
        if g == t == "-":
            raise PopulationError("Native alignment contains a double gap")
        if g == "-":
            dna.append(index)
        elif t == "-":
            rna.append(index)
        elif g != t:
            mismatches.append(index)
    actual_kind = "RNA,DNA" if dna and rna else "DNA" if dna else "RNA" if rna else "X"
    counts = (len(mismatches), len(dna) + len(rna), len(mismatches) + len(dna) + len(rna))
    if actual_kind != kind or counts != (_number(mm, "mismatches"), _number(bulges, "bulge size"),
                                         _number(total, "total edits")):
        raise PopulationError("Native edit counts disagree with the retained alignment")
    start = _number(position, "position")
    cluster_position = _number(cluster, "cluster position")
    sequence = aligned_target.replace("-", "")
    end = start + len(sequence)
    submitted = native_sequence(window)
    if end > len(submitted):
        raise PopulationError("Native target extends beyond its submitted window")
    expected = submitted[start:end]
    if strand == "-":
        expected = reverse_complement(expected)
    if sequence != expected:
        raise PopulationError("Native target sequence does not match its retained window coordinates")
    offset = len(window.get("native_guard", {}).get("prefix", ""))
    if start < offset or end > offset + len(window["sequence"]):
        # Guard DNA only supplies index context. It is never donor sequence,
        # an admitted target, or evidence about skipped genomic scope.
        return None
    native_start, native_end = start, end
    start, end = start - offset, end - offset
    admitted = (strand in profile["strands"] and matches_pam(pam, profile["pam_patterns"])
                and len(mismatches) <= profile["max_mismatches"]
                and len(dna) <= profile["max_dna_bulges"] and len(rna) <= profile["max_rna_bulges"]
                and sum((len(mismatches), len(dna), len(rna))) <= profile["max_total_edits"])
    if not admitted:
        return None
    indexes = list(range(start, end))
    if strand == "-":
        indexes.reverse()
    return {"window_id": window["window_id"], "start": start, "end": end, "strand": strand,
            "guide_alignment": ga, "target_alignment": ta, "target_indexes": indexes[:-pam_length],
            "pam_sequence": pam, "pam_indexes": indexes[-pam_length:],
            "mismatch_positions": mismatches, "dna_bulge_positions": dna, "rna_bulge_positions": rna,
            "native_row_refs": [source_ref], "native_cluster_position": cluster_position,
            "alignment_origin": "native", "native_coordinate_transform": {
                "version": WINDOW_GUARD_TRANSFORM if offset else "identity",
                "native_start": native_start, "native_end": native_end, "subtract_offset": offset}}


def alignment_key(hit):
    return tuple(hit[key] for key in ("window_id", "start", "end", "strand", "guide_alignment",
                                      "target_alignment", "pam_sequence"))


def _no_limit_check():
    pass


def extend_leading_dna(parsed, *, windows, profile, check_limits=_no_limit_check):
    """Derive missing leading DNA-gap paths from retained native anchors.

    Removing all leading guide gaps from an admitted alignment leaves the full
    guide, the same PAM, and a shorter admitted target alignment. Conversely,
    prepending exact adjacent target bases to that anchor is valid when the DNA
    and total-edit budgets allow it. No local sequence search supplies anchors.
    Raw anchors and derived paths remain distinguishable in the result.
    """
    check_limits()
    result = deepcopy(parsed)
    check_limits()
    result["native_alignment_count"] = len(parsed["hits"])
    result["derived_alignment_count"] = 0
    result["alignment_transformation"] = LEADING_DNA_TRANSFORM
    if parsed["retention_status"] != "complete":
        return result
    by_window = {window["window_id"]: window for window in windows.values()}
    seen = {alignment_key(hit) for hit in result["hits"]}
    for anchor in parsed["hits"]:
        check_limits()
        window = by_window[anchor["window_id"]]
        spent = sum(len(anchor[field]) for field in ("mismatch_positions", "dna_bulge_positions", "rna_bulge_positions"))
        available = min(profile["max_dna_bulges"] - len(anchor["dna_bulge_positions"]),
                        profile["max_total_edits"] - spent)
        for count in range(1, available + 1):
            check_limits()
            if anchor["strand"] == "+":
                start, end = anchor["start"] - count, anchor["end"]
                indexes = list(range(start, anchor["start"]))
                if start < 0:
                    break
                prefix = window["sequence"][start:anchor["start"]]
            else:
                start, end = anchor["start"], anchor["end"] + count
                indexes = list(range(end - 1, anchor["end"] - 1, -1))
                if end > len(window["sequence"]):
                    break
                prefix = reverse_complement(window["sequence"][anchor["end"]:end])
            hit = deepcopy(anchor)
            hit.update(start=start, end=end, guide_alignment="-" * count + anchor["guide_alignment"],
                       target_alignment=prefix + anchor["target_alignment"],
                       target_indexes=indexes + anchor["target_indexes"], alignment_origin="derived")
            for field in ("mismatch_positions", "dna_bulge_positions", "rna_bulge_positions"):
                hit[field] = [pos + count for pos in anchor[field]]
            hit["dna_bulge_positions"] = list(range(count)) + hit["dna_bulge_positions"]
            key = alignment_key(hit)
            if key in seen:
                continue
            if len(result["hits"]) >= profile["max_rows"]:
                result["retention_status"] = "truncated"
                result["reason_codes"].append("derived_alignment_cap")
                return result
            derivation = {"version": LEADING_DNA_TRANSFORM, "anchor": {
                field: anchor[field] for field in ("window_id", "start", "end", "strand",
                    "guide_alignment", "target_alignment", "pam_sequence")},
                "added_target_indexes": indexes, "added_target_bases": prefix,
                "native_anchor_row_refs": deepcopy(anchor["native_row_refs"])}
            hit["derivation"] = derivation
            for ref in hit["native_row_refs"]:
                ref["note"] = ("Derived alignment; native row is its shorter anchor. " + LEADING_DNA_TRANSFORM
                               + ": prepended " + str(count) + " exact target base(s) as DNA gaps.")
                ref["transformation"] = LEADING_DNA_TRANSFORM
            result["hits"].append(hit)
            result["derived_alignment_count"] += 1
            seen.add(key)
    return result


def parse_targets(paths, *, windows, design, editor, profile, check_limits=_no_limit_check):
    """Keep every distinct admitted alignment and every raw row reference.

    The pinned pipeline removes the header while filtering PAMs. Both its
    headerless result and the engine's exact ten-column header are accepted.
    """
    hits, seen, raw_rows, rejected = [], {}, 0, 0
    for path in sorted(map(Path, paths)):
        check_limits()
        sha = checksum(path)
        check_limits()
        with path.open(encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, 1):
                check_limits()
                line = line.rstrip("\r\n")
                if line_number == 1 and tuple(line.split("\t")) == HEADER:
                    continue
                if not line:
                    raise PopulationError("Blank row in native targets file")
                raw_rows += 1
                if raw_rows > profile["max_rows"]:
                    return {"hits": hits, "raw_rows": raw_rows, "rejected_rows": rejected,
                            "retention_status": "truncated", "reason_codes": ["native_row_cap"]}
                ref = {"sha256": sha, "filename": path.name, "source_id": "sha256:" + sha,
                       "line": line_number}
                hit = parse_row(line, windows=windows, guide=design["guide_sequence"],
                                geometry=editor["guide_geometry"], profile=profile, source_ref=ref)
                if hit is None:
                    rejected += 1
                    continue
                key = alignment_key(hit)
                if key in seen:
                    seen[key]["native_row_refs"].extend(hit["native_row_refs"])
                else:
                    hits.append(hit)
                    seen[key] = hit
    return {"hits": hits, "raw_rows": raw_rows, "rejected_rows": rejected,
            "retention_status": "complete", "reason_codes": []}
