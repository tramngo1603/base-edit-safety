"""Observation imports preserve source cells, units, scope and revision history."""
from copy import deepcopy
import csv
import json
from pathlib import Path
from zipfile import ZipFile

import pytest

from base_edit_safety.population.common import PopulationError, checksum, known, read_json, unknown
from base_edit_safety.population.comparison_contracts import (
    active_observations, make_record, record_registry, validate_catalog,
)
from base_edit_safety.population.observations import GENERIC_COLUMNS, import_artifact, merge_history

ROOT = Path(__file__).resolve().parents[1]
CATALOG = read_json(ROOT / "tests/fixtures/comparison/canonical-records.json")


def record(kind):
    return deepcopy(next(row for row in CATALOG if row["record_type"] == "comparison_" + kind))


def remake(row, **changes):
    kind = row.pop("record_type")
    for key in ("schema_version", "record_id", record_registry()[kind]):
        row.pop(key)
    row.update(changes)
    return make_record(kind, **row)


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


def inputs(tmp_path, format_name="json", raw_value="20", source_label="treated"):
    raw = {"measurement_key": "imported:1", "source_study": "artificial experiment",
        "source_experiment": "experiment 1", "source_replicate": "repeat 1",
        "raw_value": raw_value, "raw_unit": "percent", "source_reported_label": source_label}
    path = tmp_path / ("source." + format_name)
    if format_name == "tsv":
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=GENERIC_COLUMNS, delimiter="\t")
            writer.writeheader()
            writer.writerow(raw)
    else:
        write_json(path, [raw])
    source = remake(record("source_artifact"), sha256=checksum(path),
        source_id="artificial source import", uri="synthetic:observation-import",
        locator={"kind": "file", "member": path.name})
    normalization = record("observation_revision")
    for key in ("record_type", "schema_version", "record_id", "observation_id", "source_record_id",
                "extraction_version", "mapping_version", "normalization_schema_version",
                "extraction_map_sha256", "normalized_value", "applicability"):
        normalization.pop(key)
    normalization["conversion"] = known({"method": "Percent divided by 100", "factor": 0.01, "offset": 0})
    entry = {key: raw[key] for key in ("measurement_key", "source_study", "source_experiment", "source_replicate", "raw_unit")}
    entry.update(locator={"kind": "file", "member": "rows/1"}, expected_raw_value=known(raw_value),
        source_reported_label=known(source_label), normalization=normalization,
        normalization_unavailable_reason=None, assertions=[])
    mapping = {"record_type": "comparison_extraction_map", "schema_version": "1",
        "source_artifact_id": source["record_id"], "artifact_format": format_name,
        "extraction_version": "test-extractor-v1", "mapping_version": "test-map-v1",
        "expected_count": 1, "rows": [entry]}
    map_path = tmp_path / "map.json"
    write_json(map_path, mapping)
    return path, map_path, CATALOG + [source], mapping


@pytest.mark.parametrize("format_name", ["json", "tsv"])
def test_generic_import_reads_raw_bytes_and_explicit_conversion(tmp_path, format_name):
    path, map_path, catalog, _ = inputs(tmp_path, format_name)
    output = import_artifact(path, map_path, catalog)
    source, revision, state = output
    assert source["raw_value"] == known("20")
    assert source["source_reported_label"] == known("treated")
    assert revision["normalized_value"] == known(0.2)
    assert revision["extraction_map_sha256"] == checksum(map_path)
    assert state["active_revision_id"] == revision["observation_id"]
    validate_catalog(catalog + output)


def test_source_map_and_catalog_reject_sequence_from_a_different_locus(tmp_path):
    path, map_path, catalog, mapping = inputs(tmp_path)
    original = record("intended_locus")
    other = remake(original, gene="HBG2", locus_name="Artificial second HBG copy",
        start=original["start"] + 1000, end=original["end"] + 1000,
        path=[{**point, "pos0": point["pos0"] + 1000} for point in original["path"]])
    catalog = [*catalog, other]
    source, revision, state = import_artifact(path, map_path, catalog)
    broken = remake(deepcopy(revision), locus_ids=[other["record_id"]])
    broken_state = remake(deepcopy(state), active_revision_id=broken["record_id"], revision_ids=[broken["record_id"]])
    with pytest.raises(PopulationError, match="sequence and reported locus"):
        validate_catalog(catalog + [source, broken, broken_state])
    mapping["rows"][0]["normalization"]["locus_ids"] = [other["record_id"]]
    write_json(map_path, mapping)
    with pytest.raises(PopulationError, match="sequence and reported locus"):
        import_artifact(path, map_path, catalog)
    mapping["rows"][0]["normalization"]["sequence_configuration_id"] = unknown("The source did not resolve the two HBG copies")
    mapping["rows"][0]["normalization"]["locus_ids"] = [record("intended_locus")["record_id"], other["record_id"]]
    write_json(map_path, mapping)
    _, unresolved, _ = import_artifact(path, map_path, catalog)
    assert not unresolved["applicability"]["empirical_donor_credit"]
    assert unresolved["sequence_configuration_id"]["value"] is None


@pytest.mark.parametrize("change", ["raw_value", "artifact", "unknown_field", "duplicate_key", "boolean_count"])
def test_source_mutations_and_malformed_maps_fail(tmp_path, change):
    path, map_path, catalog, mapping = inputs(tmp_path)
    if change == "raw_value":
        mapping["rows"][0]["expected_raw_value"] = known("21")
    elif change == "artifact":
        path.write_text(path.read_text() + " ")
    elif change == "unknown_field":
        mapping["rows"][0]["invented"] = 3
    elif change == "boolean_count":
        mapping["expected_count"] = True
    if change == "duplicate_key":
        map_path.write_text('{"schema_version":"1","schema_version":"1"}')
    else:
        write_json(map_path, mapping)
    with pytest.raises(PopulationError):
        import_artifact(path, map_path, catalog)


def test_zero_and_control_label_do_not_invent_detection_or_passed_quality(tmp_path):
    path, map_path, catalog, mapping = inputs(tmp_path, raw_value="0", source_label="GFP control")
    normalization = mapping["rows"][0]["normalization"]
    normalization.update(qc_status="unknown", control_status="unknown",
        detection_limit={"value": unknown("Not reported"), "unit": "fraction"})
    write_json(map_path, mapping)
    source, observation, _ = import_artifact(path, map_path, catalog)
    assert source["source_reported_label"] == known("GFP control")
    assert observation["normalized_value"] == known(0)
    assert not observation["applicability"]["empirical_donor_credit"]
    assert "zero_without_usable_detection_limit" in observation["applicability"]["reason_codes"]


def test_revision_update_preserves_one_active_measurement(tmp_path):
    path, map_path, catalog, mapping = inputs(tmp_path)
    first = import_artifact(path, map_path, catalog)
    mapping["extraction_version"] = "test-extractor-v2"
    write_json(map_path, mapping)
    second = import_artifact(path, map_path, catalog)
    merged = merge_history(second, first)
    assert sum(r["record_type"] == "comparison_source_observation" for r in merged) == 1
    assert sum(r["record_type"] == "comparison_observation_revision" for r in merged) == 2
    current = [r for r in active_observations(catalog + merged) if r["source_record_id"] == first[0]["record_id"]]
    assert current == [second[1]]
    assert merge_history(second, merged) == merged


def test_censored_value_remains_unavailable_with_explicit_reason(tmp_path):
    path, map_path, catalog, mapping = inputs(tmp_path, raw_value="<0.1")
    mapping["rows"][0]["normalization_unavailable_reason"] = "Censored; no point estimate reported"
    mapping["rows"][0]["normalization"]["measurement_kind"] = "censored"
    write_json(map_path, mapping)
    source, observation, _ = import_artifact(path, map_path, catalog)
    assert source["raw_value"]["value"] == "<0.1"
    assert observation["normalized_value"]["value"] is None
    assert not observation["applicability"]["empirical_donor_credit"]


def test_source_decimal_conversion_does_not_round_an_excel_value(tmp_path):
    path, map_path, catalog, _ = inputs(tmp_path, raw_value="9.8000000000000007")
    source, observation, _ = import_artifact(path, map_path, catalog)
    assert source["raw_value"] == known("9.8000000000000007")
    assert observation["normalized_value"] == known("0.098000000000000007")


def test_real_schema_three_ledger_bridge_preserves_a_measurement(tmp_path):
    from base_edit_safety.pipeline.candidate_panel import build_ledger
    folder = ROOT / "examples/candidate-panel"
    ledger = build_ledger(folder / "search-profile.json", folder / "guides.tsv", folder / "intended-loci.tsv",
        [folder / "nominations.tsv"], folder / "measurements.tsv",
        source_provenance_paths=[folder / "source-provenance.json"])
    _, map_path, catalog, mapping = inputs(tmp_path)
    path = tmp_path / "ledger.json"
    write_json(path, ledger)
    source = remake(deepcopy(catalog[-1]), sha256=checksum(path))
    catalog[-1] = source
    site = next(s for s in ledger["candidates"] if s["measurements"])
    measurement = site["measurements"][0]
    mapping.update(source_artifact_id=source["record_id"], artifact_format="candidate_ledger")
    entry = mapping["rows"][0]
    entry.update(locator={"kind": "file", "member": f"candidates/{site['site_id']}/measurements/0"},
        expected_raw_value=known(str(measurement["value"])), raw_unit=measurement["unit"],
        source_study=measurement["source_version_or_citation"], source_experiment=measurement["assay"],
        source_replicate=measurement["replicate_id"], source_reported_label=known(measurement["status"]))
    entry["normalization"].update(design_id=unknown("Legacy design has no verified canonical mapping"),
        unresolved_source_design=known(measurement["guide_id"]), context_id=unknown("Source context mapping is unavailable"),
        sequence_configuration_id=unknown("Source sequence has no canonical mapping"), locus_ids=[])
    write_json(map_path, mapping)
    source_row, observation, _ = import_artifact(path, map_path, catalog)
    assert json.loads(source_row["raw_row_text"]) == measurement
    assert source_row["raw_value"]["value"] == str(measurement["value"])
    assert not observation["applicability"]["empirical_donor_credit"]
    entry["source_replicate"] = "Different replicate absent from the source ledger"
    write_json(map_path, mapping)
    with pytest.raises(PopulationError, match="replicate identity"):
        import_artifact(path, map_path, catalog)
    entry["source_replicate"] = measurement["replicate_id"]
    entry["raw_unit"] = "fraction"
    write_json(map_path, mapping)
    with pytest.raises(PopulationError, match="unit or reported outcome"):
        import_artifact(path, map_path, catalog)


def test_workbook_uses_physical_rows_and_checks_metadata_cells(tmp_path):
    _, map_path, catalog, mapping = inputs(tmp_path)
    path = tmp_path / "sparse.xlsx"
    ns = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
    with ZipFile(path, "w") as book:
        book.writestr("xl/workbook.xml", f'<workbook xmlns="{ns}" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets><sheet name="Results" sheetId="1" r:id="rId1"/></sheets></workbook>')
        book.writestr("xl/_rels/workbook.xml.rels", '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Target="worksheets/sheet1.xml"/></Relationships>')
        book.writestr("xl/worksheets/sheet1.xml", f'<worksheet xmlns="{ns}"><sheetData><row r="1"><c r="A1" t="inlineStr"><is><t>Treated</t></is></c></row><row r="6"><c r="C6"><v>20</v></c></row></sheetData></worksheet>')
    artifact = remake(deepcopy(catalog[-1]), sha256=checksum(path))
    catalog[-1] = artifact
    mapping.update(source_artifact_id=artifact["record_id"], artifact_format="xlsx")
    locator = {"kind": "table", "sheet": "Results", "table": "xl/worksheets/sheet1.xml", "row": 6, "column": "C", "cell": "C6"}
    mapping["rows"][0]["locator"] = locator
    assertion = {"locator": {**locator, "row": 1, "column": "A", "cell": "A1"}, "expected_text": "Treated"}
    mapping["rows"][0]["assertions"] = [assertion]
    write_json(map_path, mapping)
    output = import_artifact(path, map_path, catalog)
    assert output[0]["locator"]["row"] == 6
    assert output[1]["normalized_value"] == known(0.2)
    assertion["expected_text"] = "Control"
    write_json(map_path, mapping)
    with pytest.raises(PopulationError, match="metadata cell"):
        import_artifact(path, map_path, catalog)
