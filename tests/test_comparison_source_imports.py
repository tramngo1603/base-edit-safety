"""Pinned publication imports preserve source scope and reject mapping drift."""
from copy import deepcopy
from decimal import Decimal
import os
from pathlib import Path
from xml.etree import ElementTree as ET
from zipfile import ZipFile

import pytest

from base_edit_safety.population.common import PopulationError, read_json
from base_edit_safety.population.observations import merge_history
from base_edit_safety.population.source_imports import (
    SOURCES, build_source_map, import_published_source,
)


@pytest.mark.parametrize("source_name", sorted(SOURCES))
def test_unrelated_workbook_with_expected_value_count_is_rejected(tmp_path, source_name):
    source = SOURCES[source_name]
    path = tmp_path / "unrelated.xlsx"
    ns = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
    rows = "".join(f'<row r="{i}"><c r="A{i}"><v>1</v></c></row>'
                   for i in range(1, source["count"] + 1))
    with ZipFile(path, "w") as book:
        book.writestr(source["member"], f'<worksheet xmlns="{ns}"><sheetData>{rows}</sheetData></worksheet>')
    output = tmp_path / "map.json"
    with pytest.raises(PopulationError, match="different workbook"):
        import_published_source(source_name, path, output, [], [{}] * source["count"], [])
    assert not output.exists()


@pytest.mark.parametrize("source_name", sorted(SOURCES))
@pytest.mark.parametrize("delta", [-1, 1])
def test_published_source_count_is_a_fixed_contract(source_name, delta):
    with pytest.raises(PopulationError, match="fixed contract"):
        build_source_map(source_name, [], [{}] * (SOURCES[source_name]["count"] + delta), [])


def test_unregistered_published_source_is_not_guessed():
    with pytest.raises(PopulationError, match="Unregistered"):
        build_source_map("another_120_cell_workbook", [], [], [])


@pytest.fixture(scope="module")
def real_curation():
    value = os.environ.get("BES_COMPARISON_CURATION_RUNTIME")
    if not value:
        pytest.skip("Exact published artifacts remain in the private runtime; set BES_COMPARISON_CURATION_RUNTIME")
    runtime = Path(value).resolve()
    stage = runtime / "stages/03"
    paths = read_json(stage / "canonical-v2/artifact-paths.json")
    workbooks = {entry["role"]: Path(entry["path"]) for entry in paths.values()}
    return {
        "catalog": read_json(stage / "canonical-v2/catalog.json"),
        "bindings": read_json(stage / "canonical-v2/source-extraction-map.json"),
        "cells": {
            "liang_endov_hbg": read_json(stage / "preparation-v2/endov-120-source-cells.json"),
            "mayuranathan_figure4": read_json(stage / "preparation-v2/may-figure4-source-cells.json"),
        },
        "workbooks": workbooks,
    }


@pytest.mark.parametrize("source_name", sorted(SOURCES))
def test_actual_workbook_cells_import_exactly_without_donor_credit(tmp_path, real_curation, source_name):
    source = SOURCES[source_name]
    workbook = real_curation["workbooks"][source["role"]]
    cells = real_curation["cells"][source_name]
    output = import_published_source(source_name, workbook, tmp_path / "map.json",
        real_curation["catalog"], cells, real_curation["bindings"])
    raw = {row["record_id"]: row for row in output if row["record_type"] == "comparison_source_observation"}
    normalized = [row for row in output if row["record_type"] == "comparison_observation_revision"]
    assert len(raw) == len(normalized) == source["count"]
    # Read numeric cells independently from the original XML, without either extractor.
    namespace = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
    with ZipFile(workbook) as book:
        tree = ET.fromstring(book.read(source["member"]))
    xml_values = {cell.attrib["r"]: cell.findtext(namespace + "v")
                  for cell in tree.iter(namespace + "c")}
    expected = {cell["source_record_key"]: cell for cell in cells}
    for row in normalized:
        measured = raw[row["source_record_id"]]
        cell = expected[measured["measurement_key"]]
        assert measured["raw_value"]["value"] == xml_values[measured["locator"]["cell"]]
        assert Decimal(str(row["normalized_value"]["value"])) == Decimal(measured["raw_value"]["value"]) / 100
        assert measured["raw_unit"] == "percent" and row["normalized_unit"] == "fraction"
        assert measured["source_reported_label"]["value"] == cell["source_arm"]
        assert measured["source_replicate"] == cell["source_replicate"]
        assert row["subject"]["donor_id"]["value"] is None
        assert not row["applicability"]["empirical_donor_credit"]
        if cell["source_arm"] in {"GFP", "UT"}:
            assert row["design_id"]["value"] is None
            assert row["source_reported_outcome"]["value"] is None
    independent = {row["independent_unit_id"]["value"] for row in normalized
                   if row["independent_unit_id"]["value"] is not None}
    assert len(independent) == (0 if source_name == "liang_endov_hbg" else 3)
    merged = merge_history(output, output)
    assert merge_history(output, merged) == merged
    assert sum(row["record_type"] == "comparison_source_observation" for row in merged) == source["count"]


@pytest.mark.parametrize("source_name", sorted(SOURCES))
def test_reordered_source_rows_have_the_same_extraction_map(real_curation, source_name):
    expected = build_source_map(source_name, real_curation["catalog"],
        real_curation["cells"][source_name], real_curation["bindings"])
    actual = build_source_map(source_name, list(reversed(real_curation["catalog"])),
        list(reversed(real_curation["cells"][source_name])), list(reversed(real_curation["bindings"])))
    assert actual == expected


@pytest.mark.parametrize("mutation", ["duplicate_cell", "borrowed_replicate", "control_as_treated"])
def test_endov_source_mapping_drift_fails(real_curation, mutation):
    cells = deepcopy(real_curation["cells"]["liang_endov_hbg"])
    bindings = deepcopy(real_curation["bindings"])
    if mutation == "duplicate_cell":
        cells[1] = deepcopy(cells[0])
    elif mutation == "borrowed_replicate":
        cells[0]["replicate_header_cell"] = "D5"
        cells[0]["source_replicate"] = "Data 2"
    else:
        treated = next(row for row in bindings if row["source_record_key"] == cells[0]["source_record_key"])
        control = next(cell for cell in cells if cell["source_arm"] == "GFP")
        match = next(row for row in bindings if row["source_record_key"] == control["source_record_key"])
        match["design_id"] = treated["design_id"]
    with pytest.raises(PopulationError):
        build_source_map("liang_endov_hbg", real_curation["catalog"], cells, bindings)


@pytest.mark.parametrize("mutation", ["editor_swap", "dose_context_collapse", "endpoint_swap", "borrowed_donor_cell"])
def test_functional_mapping_preserves_regimen_endpoint_and_source_donor(real_curation, mutation):
    cells = deepcopy(real_curation["cells"]["mayuranathan_figure4"])
    bindings = deepcopy(real_curation["bindings"])
    by_key = {row["source_record_key"]: row for row in bindings}
    treated = next(cell for cell in cells if cell["source_arm"] == "ABE8e" and cell["source_panel"] == "Figure 4a")
    other = next(cell for cell in cells if cell["source_arm"] == "ABE7.10" and cell["source_panel"] == "Figure 4a")
    if mutation in {"editor_swap", "dose_context_collapse"}:
        field = "design_id" if mutation == "editor_swap" else "context_id"
        by_key[treated["source_record_key"]][field] = by_key[other["source_record_key"]][field]
    elif mutation == "endpoint_swap":
        hbf = next(cell for cell in cells if cell["endpoint"] == "HbF_fraction")
        by_key[hbf["source_record_key"]]["endpoint_id"] = by_key[treated["source_record_key"]]["endpoint_id"]
    else:
        cells[0]["subject_cell"] = "A7"
    with pytest.raises(PopulationError):
        build_source_map("mayuranathan_figure4", real_curation["catalog"], cells, bindings)


INDEX = {
    "design-1": {"record_type": "comparison_design", "intended_locus_ids": ["locus-hbg1", "locus-hbg2"]},
    "locus-hbg1": {"gene": "HBG1"},
    "locus-hbg2": {"gene": "HBG2"},
}


def _binding(design_id):
    return {"design_id": {"value": design_id, "reason": None if design_id else "Control arm"},
            "context_id": "comparison_experimental_context:v1:" + "0" * 64,
            "endpoint_id": "comparison_endpoint_definition:v1:" + "1" * 64}


def _liang_row(arm, site):
    return {"source_arm": arm, "source_replicate": "Data 1", "source_site": site,
            "source_subject": "HEK293T", "source_reported_validation_label": "Yes"}


def _may_row(arm):
    return {"source_arm": arm, "source_panel": "Figure 4a", "source_replicate": "rep#1",
            "source_subject": "artificial-donor"}


@pytest.mark.parametrize("name,row,design_id,outcome", [
    ("liang_endov_hbg", _liang_row("ABE7.10", "HBG2-TA"), "design-1", "inconclusive"),
    ("liang_endov_hbg", _liang_row("GFP", "HBG-OT1"), None, "untested"),
    ("mayuranathan_figure4", _may_row("ABE8e"), "design-1", "inconclusive"),
    ("mayuranathan_figure4", _may_row("UT"), None, "untested"),
])
def test_producer_adjudication_is_explicit_and_never_numeric(name, row, design_id, outcome):
    from base_edit_safety.population.source_imports import _normalization
    normalization = _normalization(name, row, _binding(design_id), INDEX)
    assert normalization["adjudicated_outcome"] == outcome
    assert isinstance(normalization["adjudication_reason"], str) and normalization["adjudication_reason"].strip()
    assert set(normalization) >= {"adjudicated_outcome", "adjudication_reason"}
    assert normalization["detection_limit"]["value"]["value"] is None
    assert normalization["normalized_unit"] == "fraction"


def test_treated_rows_keep_the_verbatim_source_label_without_credit():
    from base_edit_safety.population.common import known
    from base_edit_safety.population.source_imports import _normalization
    normalization = _normalization("liang_endov_hbg", _liang_row("ABE7.10", "HBG2-TA"), _binding("design-1"), INDEX)
    assert normalization["source_reported_outcome"] == known("Yes")
    assert normalization["locus_ids"] == ["locus-hbg2"]
    assert "withheld" in normalization["adjudication_reason"]
    may = _normalization("mayuranathan_figure4", _may_row("ABE8e"), _binding("design-1"), INDEX)
    assert may["source_reported_outcome"]["value"] is None
    assert "No categorical outcome label" in may["source_reported_outcome"]["reason"]
    assert may["locus_ids"] == ["locus-hbg1", "locus-hbg2"]


def test_control_arms_are_untested_not_negative():
    from base_edit_safety.population.source_imports import _normalization
    for name, row in [("liang_endov_hbg", _liang_row("GFP", "HBG-OT1")), ("mayuranathan_figure4", _may_row("UT"))]:
        normalization = _normalization(name, row, _binding(None), INDEX)
        assert normalization["adjudicated_outcome"] == "untested"
        assert "no treatment design" in normalization["adjudication_reason"]
        assert normalization["source_reported_outcome"]["value"] is None
        assert normalization["unresolved_source_design"]["value"].endswith("control; no treatment design")
