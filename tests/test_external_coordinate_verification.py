"""Verify HBG1/HBG2 promoter coordinates against independent fixed test data.

The fixed values were checked against the UCSC hg38 API in August 2026. `anchors.py` must reproduce
them. A failure indicates a change in `anchors.py`, the transcription-start-site constants, or the
reference assembly and requires source review before the expected values change.

The fixed tests require no network. One optional test fetches two guide-matched sequences
from UCSC and is skipped when the service is unavailable.
"""
from __future__ import annotations

import os
import socket

import pytest

from base_edit_safety import anchors


# offset -> (HBG1 plus-strand pos, HBG2 plus-strand pos, plus-strand ref base). Independently
# checked against api.genome.ucsc.edu (hg38) in August 2026 and copied here as fixed numbers. These
# are not recomputed from the TSS, so changing HBG1_TSS or HBG2_TSS causes the test to fail.
EXTERNALLY_VERIFIED = {
    113: (5_249_970, 5_254_894, "T"),
    114: (5_249_971, 5_254_895, "G"),
    117: (5_249_974, 5_254_898, "C"),
    122: (5_249_979, 5_254_903, "C"),
    123: (5_249_980, 5_254_904, "A"),
    124: (5_249_981, 5_254_905, "A"),
    175: (5_250_032, 5_254_956, "A"),
    198: (5_250_055, 5_254_979, "A"),
}


def _anchor_by_offset() -> dict[int, anchors.AnchorPosition]:
    by_offset: dict[int, anchors.AnchorPosition] = {}
    for pos in anchors.IN_SCOPE:
        # Each checked offset is unique within IN_SCOPE. Reject a duplicate instead of replacing the
        # earlier row in the dictionary.
        assert pos.promoter_offset not in by_offset, \
            f"duplicate IN_SCOPE offset {pos.promoter_offset}"
        by_offset[pos.promoter_offset] = pos
    return by_offset


def test_every_verified_offset_is_still_in_scope():
    by_offset = _anchor_by_offset()
    missing = sorted(set(EXTERNALLY_VERIFIED) - set(by_offset))
    assert not missing, f"offsets dropped from IN_SCOPE since external verification: {missing}"


def test_anchor_positions_match_externally_verified_coordinates():
    by_offset = _anchor_by_offset()
    for offset, (hbg1_pos, hbg2_pos, plus_ref) in EXTERNALLY_VERIFIED.items():
        pos = by_offset[offset]
        assert pos.hbg1_pos == hbg1_pos, f"offset {offset}: HBG1 {pos.hbg1_pos} != {hbg1_pos}"
        assert pos.hbg2_pos == hbg2_pos, f"offset {offset}: HBG2 {pos.hbg2_pos} != {hbg2_pos}"
        assert pos.genomic_plus_ref == plus_ref, \
            f"offset {offset}: plus-strand ref {pos.genomic_plus_ref!r} != {plus_ref!r}"


def test_tss_separation_is_the_verified_4924():
    assert anchors.HBG2_TSS - anchors.HBG1_TSS == 4924


# --- Optional network check for the Mayuranathan -198 guide-matched sequence ----------------------
#
# This 20-mer is the GRCh38 reference-derived -198 protospacer (Mayuranathan 2023, also reported by
# Ravi). It appears once per HBG gene copy on the plus strand at the coordinates below. This test
# fetches those bases from UCSC and skips when network access is unavailable.
MAYURANATHAN_198_PROTOSPACER = "GTGGGGAAGGGGCCCCCAAG"
MAYURANATHAN_198_STARTS = {"HBG1": 5_250_049, "HBG2": 5_254_973}


def _network_available(host: str = "api.genome.ucsc.edu", port: int = 443,
                       timeout: float = 2.0) -> bool:
    """Return whether UCSC is reachable and BES_OFFLINE is unset."""
    if os.environ.get("BES_OFFLINE"):
        return False
    try:
        socket.create_connection((host, port), timeout=timeout).close()
        return True
    except OSError:
        return False


@pytest.mark.skipif(
    not _network_available(),
    reason="offline: UCSC hg38 API not reachable (set BES_OFFLINE to force-skip)",
)
def test_mayuranathan_198_protospacer_is_at_the_verified_ucsc_coordinates():
    from urllib.error import URLError

    from base_edit_safety.io import reference

    n = len(MAYURANATHAN_198_PROTOSPACER)
    for gene, start in MAYURANATHAN_198_STARTS.items():
        try:
            # Disable the cache so the test checks the current UCSC response.
            window = reference.fetch_window("chr11", start, start + n - 1, use_cache=False)
        except (URLError, OSError) as exc:  # skip when the external service is unavailable
            pytest.skip(f"UCSC fetch failed for {gene}: {exc}")
        assert window.plus_sequence == MAYURANATHAN_198_PROTOSPACER, (
            f"{gene}: UCSC plus sequence at chr11:{start} is {window.plus_sequence!r}, "
            f"expected {MAYURANATHAN_198_PROTOSPACER!r}"
        )
