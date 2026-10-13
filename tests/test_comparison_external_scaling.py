"""Identical search strings retain distinct donor and genomic source identities."""
from copy import deepcopy
import json

import pytest

from base_edit_safety.population.common import PopulationError, checksum, digest
from base_edit_safety.population.external_search import (
    SEQUENCE_PROJECTION, _project_alignments, _sequence_aliases,
    _verify_window_source, _window_source_digest, _write_window_sources,
)


def windows():
    return [
        {"window_id": "a", "sequence": "ACGT", "donor_id": "donor-a", "path": [{"pos0": 10}]},
        {"window_id": "b", "sequence": "ACGT", "donor_id": "donor-b", "path": [{"boundary0": 90}]},
        {"window_id": "c", "sequence": "AGGT", "donor_id": "donor-c", "path": [{"pos0": 10}]},
    ]


def test_deduplication_preserves_all_members_and_distinct_strings():
    original = windows()
    before = deepcopy(original)
    aliases, projection = _sequence_aliases(original, 2)
    assert len(aliases) == 2
    assert [row["window_ids"] for row in projection] == [["a", "b"], ["c"]]
    assert [row["sequence"] for row in aliases.values()] == ["ACGT", "AGGT"]
    assert original == before


def test_projected_rows_keep_the_native_anchor_and_each_source_identity():
    _, projection = _sequence_aliases(windows(), 2)
    anchor = {"window_id": "a", "native_row_refs": [{"line": 3}],
              "derivation": {"anchor": {"window_id": "a"}, "added_target_indexes": [1]}}
    parsed = {"hits": [anchor], "reason_codes": [], "retention_status": "complete"}
    result = _project_alignments(parsed, projection, {"max_rows": 10})
    assert [row["window_id"] for row in result["hits"]] == ["a", "b"]
    assert result["hits"][1]["sequence_projection"]["representative_window_id"] == "a"
    assert result["hits"][1]["sequence_projection"]["target_window_id"] == "b"
    assert result["hits"][1]["derivation"]["anchor"]["window_id"] == "a"
    assert parsed["hits"][0]["native_row_refs"] == [{"line": 3}]
    assert result["projected_alignment_count"] == 2


def test_projection_capacity_never_claims_complete_retention():
    _, projection = _sequence_aliases(windows(), 2)
    parsed = {"hits": [{"window_id": "a", "native_row_refs": []}], "reason_codes": [], "retention_status": "complete"}
    result = _project_alignments(parsed, projection, {"max_rows": 1})
    assert len(result["hits"]) == 1
    assert result["retention_status"] == "truncated"
    assert result["reason_codes"] == ["projected_alignment_cap"]


def test_streamed_windows_and_projection_are_bound_to_effective_input(tmp_path):
    source = windows()
    path = tmp_path / "source-windows.jsonl"
    expected = _window_source_digest(source)
    _write_window_sources(path, source, expected)
    assert checksum(path) == expected
    assert [json.loads(line) for line in path.read_text().splitlines()] == source
    _, projection = _sequence_aliases(source, 2)
    mapping = tmp_path / "sequence-projection.json"
    mapping.write_text(json.dumps(projection))
    payload = {"adapter_contract": "materialized-donor-windows-v2",
               "window_source": {"filename": path.name, "sha256": expected, "count": 3},
               "sequence_projection": {"version": SEQUENCE_PROJECTION, "sha256": digest(projection)}}
    _verify_window_source(tmp_path, payload)
    mapping.write_text("[]")
    with pytest.raises(PopulationError, match="projection differs"):
        _verify_window_source(tmp_path, payload)
    mapping.write_text(json.dumps(projection))
    path.write_text(path.read_text().replace("donor-b", "donor-wrong"))
    with pytest.raises(PopulationError, match="window source differs"):
        _verify_window_source(tmp_path, payload)
