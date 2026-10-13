"""Reference checks may span a search boundary without enlarging the search."""
import hashlib
import json

import pytest

from base_edit_safety.population.common import PopulationError
from base_edit_safety.population.comparison_inputs import _ValidationReference
from base_edit_safety.population.variants import Reference


def reference(tmp_path, name, intervals):
    rows = [{"contig": "chrSYN", "start": start, "end": start + len(sequence),
             "sequence": sequence, "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest()}
            for start, sequence in intervals]
    path = tmp_path / name
    path.write_text(json.dumps({"record_type": "reference_windows", "schema_version": "1",
                                "genome_build": "GRCh38", "windows": rows}))
    policy = {"genome_build": "GRCh38", "contig_alias_map": {"chrSYN": "chrSYN"},
              "permitted_contigs": ["chrSYN"]}
    result = Reference(path, policy)
    return {"sha256": result.sha256}, result


def test_adjacent_validation_spans_do_not_expand_search_reference(tmp_path):
    _, primary = reference(tmp_path, "search.json", [(10, "ACG"), (16, "TCA")])
    source, auxiliary = reference(tmp_path, "support.json", [(8, "TT"), (13, "GAT")])
    checked = _ValidationReference(primary, [(source, auxiliary)])
    assert checked.fetch("chrSYN", 8, 19) == "TTACGGATTCA"
    assert checked.checked_intervals == [
        {"contig": "chrSYN", "start": 8, "end": 10, "source_sha256": source["sha256"], "purpose": "VCF REF validation only"},
        {"contig": "chrSYN", "start": 13, "end": 16, "source_sha256": source["sha256"], "purpose": "VCF REF validation only"},
    ]
    assert primary.bounds("chrSYN") == [[10, 13], [16, 19]]
    with pytest.raises(PopulationError, match="lacks required global coordinates"):
        primary.fetch("chrSYN", 8, 19)


def test_primary_only_ref_check_needs_no_auxiliary_credit(tmp_path):
    _, primary = reference(tmp_path, "search.json", [(10, "ACG")])
    checked = _ValidationReference(primary, [])
    assert checked.fetch("chrSYN", 10, 13) == "ACG"
    assert checked.checked_intervals == []


def test_gap_in_acquired_reference_cannot_be_filled(tmp_path):
    _, primary = reference(tmp_path, "search.json", [(10, "ACG"), (16, "TCA")])
    source, auxiliary = reference(tmp_path, "support.json", [(8, "TT"), (14, "AT")])
    checked = _ValidationReference(primary, [(source, auxiliary)])
    with pytest.raises(PopulationError, match="lacks required global coordinates"):
        checked.fetch("chrSYN", 8, 19)
    assert checked.checked_intervals == []


def test_conflicting_validation_sources_are_rejected(tmp_path):
    _, primary = reference(tmp_path, "search.json", [(10, "ACG")])
    first = reference(tmp_path, "support-a.json", [(13, "GAT")])
    second = reference(tmp_path, "support-b.json", [(15, "C")])
    checked = _ValidationReference(primary, [first, second])
    with pytest.raises(PopulationError, match="sources disagree"):
        checked.fetch("chrSYN", 10, 16)
    assert checked.checked_intervals == []


@pytest.mark.parametrize("start,end", [(10, 10), (12, 11), (-1, 11)])
def test_invalid_reference_interval_stays_invalid(tmp_path, start, end):
    _, primary = reference(tmp_path, "search.json", [(10, "ACG")])
    checked = _ValidationReference(primary, [])
    with pytest.raises(PopulationError, match="Invalid reference interval"):
        checked.fetch("chrSYN", start, end)


def test_unknown_contig_preserves_reference_error(tmp_path):
    _, primary = reference(tmp_path, "search.json", [(10, "ACG")])
    checked = _ValidationReference(primary, [])
    with pytest.raises(PopulationError, match="Reference contig does not resolve"):
        checked.fetch("chrMISSING", 10, 13)
