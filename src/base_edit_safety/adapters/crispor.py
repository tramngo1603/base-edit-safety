"""Calls CRISPOR to find guides and search for their off-targets.

You install CRISPOR yourself, with its own indexed genome, under its own license (the core tool has
its own terms; the aligners it bundles are GPL). This adapter runs it as a separate process and
never copies its source, binaries, or genome here. Tell it where the tool lives with the environment
variables listed in config.py.

Each guide comes back as a CandidateGuide with no editor attached yet. For each one, the adapter
locates the protospacer and PAM in the reference window to work out its position and strand, then
attaches what CRISPOR said about its off-targets.
"""
from __future__ import annotations

import csv
import os
import subprocess
import tempfile

from ..config import CrisporConfig, load_crispor_config
from ..io.reference import fetch_window
from ..pipeline.datatypes import reverse_complement
from ..pipeline.datatypes import CandidateGuide

PROTO_LEN = 20
# Keep reported NAG/NGA matches even below CRISPOR's default MIT score cutoff of 1.0.
# Retain their scores for review; this does not establish an exhaustive search.
MIN_ALT_PAM_SCORE = 0.0
CRISPRON_CONTEXT_PREFIX = 4
CRISPRON_CONTEXT_SUFFIX = 3


class AmbiguousTargetMapping(ValueError):
    """Raised when a CRISPOR target maps to more than one site in the reference window."""


def _target_context_30mer(window, protospacer: str, pam: str, strand: str,
                          proto_5p_genomic: int) -> str | None:
    """Return the 4+20+3+3 target context used by CRISPRon-BE when it fits the window."""
    if len(protospacer) != PROTO_LEN or len(pam) != 3:
        return None
    if strand == "+":
        low = proto_5p_genomic - CRISPRON_CONTEXT_PREFIX
        high = proto_5p_genomic + PROTO_LEN + len(pam) + CRISPRON_CONTEXT_SUFFIX - 1
    else:
        low = proto_5p_genomic - (
            PROTO_LEN + len(pam) + CRISPRON_CONTEXT_SUFFIX - 1
        )
        high = proto_5p_genomic + CRISPRON_CONTEXT_PREFIX
    if low < window.start_1based or high > window.end_1based:
        return None
    start = low - window.start_1based
    context = window.plus_sequence[start:start + 30]
    if strand == "-":
        context = reverse_complement(context)
    expected_target = protospacer + pam
    if len(context) != 30 or context[4:27] != expected_target:
        raise ValueError(
            f"reference context does not contain mapped target {expected_target} at "
            f"{window.chrom}:{proto_5p_genomic} ({strand})"
        )
    return context


def _find_all(haystack: str, needle: str) -> list[int]:
    """Return every matching start index, including overlapping matches."""
    hits, idx = [], haystack.find(needle)
    while idx != -1:
        hits.append(idx)
        idx = haystack.find(needle, idx + 1)
    return hits


def _summarize_sequence_matches(sites: list[dict]) -> dict:
    """Count reported sequence matches and sum their CFD scores.

    Keep perfect matches until a separate coordinate check identifies the intended targets.
    These counts and score sums are not probabilities of binding or editing.
    """
    histogram: dict[str, int] = {}
    cfd_burden = 0.0
    perfect_match_candidate_count = 0
    for s in sites:
        try:
            mm = int(s.get("mismatchCount"))
        except (TypeError, ValueError):
            continue
        if mm == 0:
            perfect_match_candidate_count += 1
        histogram[str(mm)] = histogram.get(str(mm), 0) + 1
        try:
            cfd_burden += float(s.get("cfdOfftargetScore"))
        except (TypeError, ValueError):
            pass
    return {
        "mm_histogram": histogram,
        "cfd_weighted_burden": round(cfd_burden, 4),
        "perfect_match_candidate_count": perfect_match_candidate_count,
        "perfect_match_classification": "requires intended-locus coordinate match",
        # Kept for record compatibility. No paralog is inferred from mismatch count alone.
        "paralog_count": 0,
        "offtarget_count_excl_paralog": sum(histogram.values()),
    }


class CrisporGuideDesigner:
    """Find guides and ask CRISPOR for matching reference-genome sequences."""

    name = "CRISPOR"

    def __init__(self, config: CrisporConfig | None = None, pam: str | None = None,
                 tss_for_gene=None, min_alt_pam_score: float = MIN_ALT_PAM_SCORE):
        self.config = config or load_crispor_config()
        self.pam = (pam or self.config.pam).strip().upper()
        if not self.pam:
            raise ValueError("CRISPOR PAM pattern cannot be empty")
        self.min_alt_pam_score = min_alt_pam_score
        # Map a gene symbol to its TSS; defaults to the values in anchors.py.
        if tss_for_gene is None:
            from .. import anchors
            tss_for_gene = {"HBG1": anchors.HBG1_TSS, "HBG2": anchors.HBG2_TSS}
        self.tss_for_gene = tss_for_gene

    def design(self, chrom: str, start_1based: int, end_1based: int, gene: str) -> list[CandidateGuide]:
        window = fetch_window(chrom, start_1based, end_1based)
        rows, offtarget_rows = self._run_crispor(
            window.plus_sequence, f"{gene}_{chrom}_{start_1based}_{end_1based}"
        )
        sites_by_guide: dict[str, list[dict]] = {}
        for site in offtarget_rows:
            sites_by_guide.setdefault(site.get("guideId", ""), []).append(site)
        guides: list[CandidateGuide] = []
        self.last_ambiguous: list[str] = []   # flagged, not silently mapped to the first hit
        for row in rows:
            try:
                guide = self._row_to_candidate(row, window, gene, sites_by_guide)
            except AmbiguousTargetMapping as exc:
                self.last_ambiguous.append(str(exc))
                continue
            if guide is not None:
                guides.append(guide)
        return guides

    def tool_metadata(self) -> dict:
        """What this tool was for this run: paths, genome index, and git commit if there is one."""
        from ..pipeline import environment
        return environment.tool_environment(
            self.name, repo_root=self.config.repo_root,
            python_executable=self.config.python_executable, genome=self.config.genome, pam=self.pam,
            min_alt_pam_score=self.min_alt_pam_score,
        )

    def _run_crispor(self, window_seq: str, seq_id: str) -> tuple[list[dict], list[dict]]:
        with tempfile.TemporaryDirectory() as tmp:
            in_fa = os.path.join(tmp, "in.fa")
            out_tsv = os.path.join(tmp, "out.tsv")
            off_tsv = os.path.join(tmp, "off.tsv")
            with open(in_fa, "w") as handle:
                handle.write(f">{seq_id}\n{window_seq}\n")
            cmd = [
                self.config.python_executable, "crispor.py", self.config.genome,
                in_fa, out_tsv, "-p", self.pam, "--noEffScores",
                "--minAltPamScore", str(self.min_alt_pam_score),
                "--genomeDir", self.config.genome_dir, "-o", off_tsv,
            ]
            completed = subprocess.run(
                cmd, cwd=self.config.repo_root, capture_output=True, text=True
            )
            if completed.returncode != 0 or not os.path.isfile(out_tsv):
                raise RuntimeError(
                    "External guide-design tool failed.\n"
                    f"stderr (tail):\n{completed.stderr[-2000:]}"
                )
            with open(out_tsv) as handle:
                guide_rows = [dict(r) for r in csv.DictReader(handle, delimiter="\t")]
            offtarget_rows = []
            if os.path.isfile(off_tsv):
                with open(off_tsv) as handle:
                    offtarget_rows = [dict(r) for r in csv.DictReader(handle, delimiter="\t")]
            return guide_rows, offtarget_rows

    def _map_target(self, window, target_seq: str):
        """Find target_seq in the window.

        Returns (strand, genomic coordinate of the 5' end), or (None, None) if it is not there. If
        it turns up more than once, this raises rather than picking the first hit.
        """
        plus = window.plus_sequence
        plus_hits = _find_all(plus, target_seq)
        minus_hits = _find_all(reverse_complement(plus), target_seq)
        total = len(plus_hits) + len(minus_hits)
        if total == 0:
            return None, None
        if total > 1:
            raise AmbiguousTargetMapping(
                f"target {target_seq} appears at {total} places in the reference window "
                f"(plus={plus_hits}, minus={minus_hits}), so there is no way to say which one it is."
            )
        if plus_hits:
            return "+", window.start_1based + plus_hits[0]   # 5' end at low coord on plus strand
        return "-", window.end_1based - minus_hits[0]         # 5' end at high coord on minus strand

    def _row_to_candidate(self, row: dict, window, gene: str,
                          sites_by_guide: dict | None = None) -> CandidateGuide | None:
        target_seq = row.get("targetSeq")
        if not target_seq or len(target_seq) < PROTO_LEN + 1:
            return None
        protospacer = target_seq[:PROTO_LEN]
        pam = target_seq[PROTO_LEN:]

        strand, proto_5p_genomic = self._map_target(window, target_seq)
        if strand is None:
            return None  # not present in this window (e.g. the guide spans the edge)

        tool_guide_id = row.get("guideId", "")
        target_context_30mer = _target_context_30mer(
            window, protospacer, pam, strand, proto_5p_genomic
        )
        sites = (sites_by_guide or {}).get(tool_guide_id, [])
        offtarget = {
            "source": self.name,
            "mit_specificity": row.get("mitSpecScore"),
            "cfd_specificity": row.get("cfdSpecScore"),
            "offtarget_count": row.get("offtargetCount"),  # raw count, includes the paralog site
            "genome_locus": row.get("targetGenomeGeneLocus"),
            "tool_guide_id": tool_guide_id,
            **_summarize_sequence_matches(sites),
        }
        return CandidateGuide(
            guide_id=f"{gene}_{tool_guide_id or 'NA'}",
            protospacer=protospacer,
            pam=pam,
            strand=strand,
            genome_build=window.genome_build,
            chrom=window.chrom,
            gene=gene,
            gene_strand="-",
            tss_1based=self.tss_for_gene[gene],
            proto_5p_genomic=proto_5p_genomic,
            offtarget=offtarget,
            offtarget_sites=sites,
            pam_class=self.pam,
            target_context_30mer=target_context_30mer,
        )
