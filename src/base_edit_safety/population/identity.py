"""Candidate identity is independent of donors, scores and padded intervals."""
from __future__ import annotations

from .common import PopulationError, digest, identity

IUPAC = dict(zip("ACGTRYSWKMBDHVN", map(set, ["A", "C", "G", "T", "AG", "CT", "CG", "AT", "GT", "AC", "CGT", "AGT", "ACT", "ACG", "ACGT"])))


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def pam_matches(pam, patterns):
    return len(pam) == 3 and any(all(b in IUPAC[p] for b, p in zip(pam, pattern)) for pattern in patterns)


def matches_pam(sequence, patterns):
    """Return whether a sequence matches any same-length IUPAC PAM pattern."""
    return any(len(sequence) == len(pattern) and all(base in IUPAC[code]
               for base, code in zip(sequence, pattern)) for pattern in patterns)


def canonical_path(path):
    result = []
    for point in path:
        if point["kind"] == "reference":
            result.append({"kind": "reference", "pos0": point["pos0"]})
        elif point["kind"] == "insertion":
            result.append({"kind": "insertion", "boundary0": point["boundary0"], "offset0": point["offset0"]})
        else:
            raise PopulationError("Unknown sequence path entry")
    return result


def bounding_interval(path):
    positions = [p["pos0"] for p in path if p["kind"] == "reference"]
    if positions:
        return [min(positions), max(positions)+1]
    boundaries = {p["boundary0"] for p in path}
    if len(boundaries) != 1:
        raise PopulationError("Insertion-only target must have one normalized boundary")
    return [next(iter(boundaries))] * 2


def definition_hashes(guide):
    if guide.get("purpose") == "computational_only":
        return digest({"spacer_sequence": guide["spacer_sequence"]}), digest({
            "purpose": "computational_only", "pam_class": guide["pam_class"]})
    return digest({"spacer_sequence": guide["spacer_sequence"]}), digest({
        k: guide[k] for k in ("editor", "editor_class", "pam_class", "editing_window_start_1based", "editing_window_end_1based_inclusive")})


def candidate_identity(reference_sha256, contig, strand, guide, target, path):
    guide_sha, editor_sha = definition_hashes(guide)
    payload = {"identity_schema": "candidate_identity_1", "reference_sha256": reference_sha256,
               "contig": contig, "strand": strand, "guide_definition_sha256": guide_sha,
               "editor_definition_sha256": editor_sha, "protospacer_guide_orientation": target[:20],
               "pam_guide_orientation": target[20:], "target_path": canonical_path(path)}
    return identity("cfg", payload), payload
