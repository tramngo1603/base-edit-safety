"""Check the accepted edit window, promoter coordinates and annotation release.

The definitions are documented in docs/architecture.md and docs/coordinates.md.
The per-gene coordinate checks are in tests/test_ch_gene_spans.py.
"""
from __future__ import annotations

from base_edit_safety import annotation
from base_edit_safety.pipeline import datatypes
from base_edit_safety.pipeline import enumerate_guides


def test_bedict_edit_window_is_frozen():
    # Protospacer positions 3-10 inclusive == range(3, 11). See docs/architecture.md.
    assert datatypes.EDIT_WINDOW == range(3, 11)
    assert list(datatypes.EDIT_WINDOW) == [3, 4, 5, 6, 7, 8, 9, 10]


def test_promoter_windows_are_frozen():
    # -210..-80 from each cap site on GRCh38, chr11. See docs/coordinates.md.
    assert enumerate_guides.WINDOWS == {
        "HBG1": ("chr11", 5_249_937, 5_250_067),
        "HBG2": ("chr11", 5_254_861, 5_254_991),
    }


def test_ch_gene_spans_ensembl_release_is_frozen():
    # A different annotation release can change the gene-overlap boundaries.
    assert annotation.CH_GENE_SPANS_ENSEMBL_RELEASE == "116"
