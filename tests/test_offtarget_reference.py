"""Compare published off-target candidate panels without treating them as measurements.

The panel that motivated this: Ravi 2022 Supplementary file 3 lists eleven sites for gRNA 11, which
is the workflow's 61forw. A prior run reported twelve sites at 3 mismatches for the same guide, but
comparison was unavailable because only CH-driver hits retained coordinates.
"""
from __future__ import annotations

from base_edit_safety import offtarget_reference as ref
from base_edit_safety.pipeline.enumerate_guides import (
    PERSIST_SITES_UP_TO_MISMATCHES, _low_mismatch_sites,
)

PANEL = ref.RAVI_2022_GRNA11


def test_the_panel_is_attached_to_the_guide_it_belongs_to():
    """The source table has no guide column; the assignment was derived by alignment."""
    assert PANEL.guide_names == {"ravi": "gRNA 11", "ours": "61forw"}
    assert ref.panel_for("GGCAAGGCTGGCCAACCCAT") is PANEL
    assert ref.panel_for("CTTGACCAATAGCCTTGACA") is None      # the -113 guide has no panel


def test_the_panel_is_recorded_as_a_nomination_not_a_measurement():
    """Do not treat a tested-site list as evidence that editing was detected.

    The source table lists examined sites but provides no editing rate or detection result.
    """
    assert PANEL.is_measurement is False
    assert "analysis" in PANEL.selection_basis or "analyzed" in PANEL.selection_basis


def test_importing_a_panel_does_not_upgrade_base_editor_context():
    """Two search methods returning a site do not establish a sequencing result."""
    result = ref.compare_to_panel([], PANEL)
    assert result["is_measurement"] is False
    assert "not_assessed_here" in result["note"]


def test_all_eleven_sites_are_three_mismatch():
    assert PANEL.mismatch_counts == {3: 11}
    assert len(PANEL.sites) == 11


def test_comparison_reports_what_our_search_missed():
    ours = [{"chrom": "chr1", "start": 25_039_865},        # OT1 exactly
            {"chrom": "chr9", "start": 118_816_550},        # OT7, 5 bp off, still a match
            {"chrom": "chr3", "start": 1_000_000}]          # not in the panel
    out = ref.compare_to_panel(ours, PANEL)
    assert "OT1" in out["found_by_our_search"]
    assert "OT7" in out["found_by_our_search"]
    assert len(out["missed_by_our_search"]) == 9
    assert out["our_site_count"] == 3


def test_a_far_away_coordinate_is_not_counted_as_a_match():
    ours = [{"chrom": "chr1", "start": 25_039_865 + 5_000}]
    assert ref.compare_to_panel(ours, PANEL)["found_by_our_search"] == []


def test_low_mismatch_sites_are_kept_and_the_noise_tail_is_not():
    """Retain the coordinates required for the published-panel comparison."""
    sites = [{"mismatchCount": "0", "chrom": "chr11", "start": "5254905"},   # the paralog
             {"mismatchCount": "2", "chrom": "chr1", "start": "100"},
             {"mismatchCount": "3", "chrom": "chr9", "start": "118816545"},
             {"mismatchCount": "4", "chrom": "chr7", "start": "9"}]
    kept = _low_mismatch_sites(sites)
    assert [k["mismatch_count"] for k in kept] == [2, 3]
    assert PERSIST_SITES_UP_TO_MISMATCHES == 3


def test_the_paralog_is_never_kept_as_an_off_target():
    kept = _low_mismatch_sites([{"mismatchCount": "0", "chrom": "chr11", "start": "5254905"}])
    assert kept == []
