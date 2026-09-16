"""The included no-data example preserves the central review distinctions."""
from base_edit_safety import anchors
from base_edit_safety.pipeline.example_record import build_example


def test_example_has_both_intended_hbg_copies_without_measurements_or_rates():
    result = build_example()

    assert result["example"] == {
        "synthetic": True,
        "external_files_used": False,
        "predictor_run": False,
        "checked_in_sources": ["GRCh38 anchor table", "evidence registry"],
        "registry_detail": "aggregate_only",
        "purpose": "record-shape inspection only",
    }
    assert result["count"] == 2
    assert {record["candidate"]["gene"] for record in result["records"]} == {
        "HBG1",
        "HBG2",
    }
    for record in result["records"]:
        intended = record["predicted_edit_spectrum"]["on_target"]
        assert len(intended) == 1
        assert intended[0]["predicted_rate"] is None
        assert intended[0]["confidence"] == "not_evaluated"
        assert record["empirical_support"] == []
        assert record["hbg_paralog_flags"]["flagged"] is True
        assert record["hbg_paralog_flags"]["paralog_spacing_bp"] == anchors.PARALOG_SPACING_BP
        assert record["off_target"]["candidate_sites"] == []
        assert all(
            "citations" not in matched_entry
            for edit_match in record["registry_match"]["by_edit"]
            for matched_entry in edit_match["matched_entries"]
        )
