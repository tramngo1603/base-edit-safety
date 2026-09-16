"""An output file written by a different build must not load.

A stale enumeration file can outlive the code that created it. These tests pin the guard that
rejects that file before its records enter a later analysis.
"""
from __future__ import annotations

import json

import pytest

from base_edit_safety.io import outputs


def _current_doc(**overrides) -> dict:
    """A document shaped the way the current writer shapes one."""
    doc = {
        "artifact_schema_version": outputs.ARTIFACT_SCHEMA_VERSION,
        "scope": outputs.SCOPE,
        "meta": {"step": "2 (guides and their edits) + 4 (CH-driver annotation)",
                 "reproducibility": {"genome_build": "GRCh38"}},
        "records": [],
    }
    doc.update(overrides)
    return doc


def test_a_current_document_has_no_complaints():
    assert outputs.schema_complaints(_current_doc()) == []


def test_missing_version_is_caught():
    doc = _current_doc()
    del doc["artifact_schema_version"]
    complaints = outputs.schema_complaints(doc)
    assert any("no artifact_schema_version" in c for c in complaints)


def test_wrong_version_is_caught():
    complaints = outputs.schema_complaints(_current_doc(artifact_schema_version="1"))
    assert any("artifact_schema_version is '1'" in c for c in complaints)


@pytest.mark.parametrize("retired", ["disclaimer", "claim_disclaim"])
def test_retired_header_keys_are_caught(retired):
    complaints = outputs.schema_complaints(_current_doc(**{retired: "anything"}))
    assert any(retired in c and "retired" in c for c in complaints)


def test_spine_in_the_step_label_is_caught():
    doc = _current_doc()
    doc["meta"]["step"] = "2 (locus->guides->spine) + 4 (CH-driver annotation)"
    complaints = outputs.schema_complaints(doc)
    assert any("spine" in c for c in complaints)


def test_empty_reproducibility_is_caught():
    doc = _current_doc()
    doc["meta"]["reproducibility"] = {}
    assert any("reproducibility" in c for c in outputs.schema_complaints(doc))


def test_legacy_document_reports_all_schema_errors():
    """Report each incompatible field in a constructed legacy document."""
    stale = {
        "disclaimer": "PREDICTED SHORTLIST FOR WET-LAB CONFIRMATION ...",
        "claim_disclaim": {"what_this_is": "..."},
        "meta": {"step": "2 (locus->guides->spine) + 4 (CH-driver annotation)"},
        "record_count": 158,
        "records": [],
    }
    complaints = outputs.schema_complaints(stale)
    assert len(complaints) >= 4, complaints
    joined = " ".join(complaints)
    for tell in ("no artifact_schema_version", "disclaimer", "claim_disclaim", "spine",
                 "reproducibility", "scope"):
        assert tell in joined, f"the guard misses the {tell!r} tell"


def test_read_enumeration_raises_and_names_the_mismatch(tmp_path, monkeypatch):
    path = tmp_path / "enumeration_HBG1_HBG2.json"
    path.write_text(json.dumps({"disclaimer": "x", "meta": {"step": "spine"}, "records": []}))
    monkeypatch.setattr(outputs, "enumeration_path", lambda: str(path))
    with pytest.raises(outputs.StaleArtifact) as exc:
        outputs.read_enumeration()
    message = str(exc.value)
    assert "was not written by this build" in message
    assert "spine" in message
    assert "enumerate_guides" in message, "the error should say how to regenerate the file"


def test_read_enumeration_says_so_when_the_file_is_absent(tmp_path, monkeypatch):
    monkeypatch.setattr(outputs, "enumeration_path", lambda: str(tmp_path / "nope.json"))
    with pytest.raises(outputs.StaleArtifact, match="no enumeration file"):
        outputs.read_enumeration()


def test_writer_stamps_the_version(tmp_path, monkeypatch):
    monkeypatch.setattr(outputs, "OUTPUT_DIR", str(tmp_path))
    path = outputs.write_records([], filename="x.json", meta={"step": "t"})
    written = json.loads(open(path).read())
    assert written["artifact_schema_version"] == outputs.ARTIFACT_SCHEMA_VERSION
    assert "scope" in written
    assert not any(k in written for k in ("disclaimer", "claim_disclaim"))
