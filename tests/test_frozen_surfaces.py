"""Check scientific constants that can change only through documented review.

This file keeps the BE-DICT edit window, promoter windows, and Ensembl release in one fixed test.
The promoter coordinates are documented in `docs/coordinates.md`; full gene spans are covered in
`tests/test_ch_gene_spans.py`. A change must update both the implementation and expected values.
"""
from __future__ import annotations

from base_edit_safety import annotation
from base_edit_safety.pipeline import datatypes
from base_edit_safety.pipeline import enumerate_guides


def test_bedict_edit_window_is_frozen():
    # Guide-matched positions 3-10 inclusive are represented as range(3, 11).
    assert datatypes.EDIT_WINDOW == range(3, 11)
    assert list(datatypes.EDIT_WINDOW) == [3, 4, 5, 6, 7, 8, 9, 10]


def test_promoter_windows_are_frozen():
    # -210..-80 from each transcription start site on GRCh38 chromosome 11; see docs/coordinates.md.
    assert enumerate_guides.WINDOWS == {
        "HBG1": ("chr11", 5_249_937, 5_250_067),
        "HBG2": ("chr11", 5_254_861, 5_254_991),
    }


def test_ch_gene_spans_ensembl_release_is_frozen():
    # Per-gene expected values live in tests/test_ch_gene_spans.py. Fix the release because a new
    # release can change CH-driver gene boundaries.
    assert annotation.CH_GENE_SPANS_ENSEMBL_RELEASE == "116"
