"""Tripwires for the surfaces the 2026-08 external review froze.

None of these values is meant to drift on its own. Each is set deliberately and justified elsewhere:
the BE-DICT edit window in docs/architecture.md, the promoter windows in
docs/coordinates.md, the Ensembl pin (with its full per-gene fixture) in
tests/test_ch_gene_spans.py. This file gathers the bare constants in one place so that ANY future
diff touching them must change this test too, forcing the change to be deliberate rather than a
silent side effect.
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
    # The full per-gene span fixture lives in tests/test_ch_gene_spans.py; this pins the version
    # string alongside the other frozen surfaces, because a release bump silently moves the
    # CH-driver flag boundaries.
    assert annotation.CH_GENE_SPANS_ENSEMBL_RELEASE == "116"
