"""Source-specific maps for the verified EndoV and HBG functional workbooks."""
from __future__ import annotations

from decimal import Decimal, InvalidOperation
import json
from pathlib import Path
import re

from .common import PopulationError, checksum, known, read_json, unknown
from .comparison_contracts import validate_catalog
from .observations import import_artifact


SOURCES = {
    "liang_endov_hbg": {
        "sha256": "1951af1f7a3dba33b195dcef62e5208d323f6308464a433d82ff2f4218f89177",
        "role": "liang_workbook", "sheet": "Figure 2d", "member": "xl/worksheets/sheet2.xml",
        "count": 120, "doi": "10.1038/s41467-018-07988-z",
    },
    "mayuranathan_figure4": {
        "sha256": "d9580e9d4ff141bd074b9988b37d98ce59ccd543b33b0dd8f441ef2e69ba1ed9",
        "role": "may_figures", "sheet": "Figure 4", "member": "xl/worksheets/sheet4.xml",
        "count": 135, "doi": "10.1038/s41588-023-01434-7",
    },
}


def _source(name):
    try:
        return SOURCES[name]
    except KeyError as exc:
        raise PopulationError("Unregistered published-source extraction") from exc


def _locator(source, address):
    match = re.fullmatch(r"([A-Z]+)([1-9][0-9]*)", address)
    if match is None:
        raise PopulationError("Source cell must be a physical A1 address")
    return {"kind": "table", "sheet": source["sheet"], "table": source["member"],
            "row": int(match[2]), "column": match[1], "cell": address}


def _layout(name):
    if name == "liang_endov_hbg":
        return {f"{column}{row}": {"arm": "ABE7.10" if column in "CDE" else "GFP",
                    "panel": "Figure 2d", "day": None}
                for row in range(6, 26) for column in "CDEFGH"}
    result = {}
    for panel, first, last, columns, day in [("Figure 4a", 4, 12, "CDEFGHIJK", 6),
            ("Figure 4b", 17, 25, "CDE", 14), ("Figure 4c", 30, 38, "CDE", 21)]:
        for row in range(first, last + 1):
            for i, column in enumerate(columns):
                arm = ("UT", "ABE7.10", "ABE8e")[i // 3 if panel.endswith("a") else i]
                result[f"{column}{row}"] = {"arm": arm, "panel": panel, "day": day}
    return result


def _check_binding(name, row, binding, expected, index):
    design_id = binding["design_id"]["value"]
    control = expected["arm"] in {"GFP", "UT"}
    if control != (design_id is None) or control != (row["draft_design_key"] is None):
        raise PopulationError("Source control/treatment and canonical design binding disagree")
    if design_id is not None:
        design = index.get(design_id)
        if design is None or design["record_type"] != "comparison_design":
            raise PopulationError("Published observation design is unresolved")
        editor = index[design["editor_definition_id"]]
        guide = "GTGGGGAAGGGGCCCCCAAG" if name == "liang_endov_hbg" else "AGATATTTGCATTGAGATAG"
        variant = "ABE7.10-SpCas9" if name == "liang_endov_hbg" else (
            "ABE8e-SpCas9-NG" if expected["arm"] == "ABE8e" else "ABE7.10-NG")
        if design["guide_sequence"] != guide or editor["variant_name"] != variant:
            raise PopulationError("Guide/editor binding differs from the exact published combination")
    context = index.get(binding["context_id"])
    endpoint = index.get(binding["endpoint_id"])
    if context is None or context["record_type"] != "comparison_experimental_context" or endpoint is None or endpoint["record_type"] != "comparison_endpoint_definition":
        raise PopulationError("Published source context or endpoint is unresolved")
    if binding["context_id"] not in endpoint["required_context_ids"]:
        raise PopulationError("Published source context is outside the endpoint contract")
    if name == "liang_endov_hbg":
        endpoint_name = "intended_edit_fraction" if row["source_site"].endswith("-TA") else "unintended_edit_fraction"
        timepoint, unit = 48, "hour"
    else:
        endpoint_name = row["endpoint"]
        timepoint, unit = expected["day"], "day"
        expected_endpoint = {"Figure 4b": "F_cell_fraction", "Figure 4c": "HbF_fraction"}.get(expected["panel"])
        if expected_endpoint is None:
            expected_endpoint = "intended_edit_fraction" if row["guide_position_1based"] == 5 else "bystander_edit_fraction"
        if endpoint_name != expected_endpoint or row["measurement_day"] != timepoint:
            raise PopulationError("Source panel, editing position or timepoint disagrees")
        expected_dose = {"UT": None, "ABE7.10": 5, "ABE8e": 8}[expected["arm"]]
        if context["dose"]["value"]["value"] != expected_dose:
            raise PopulationError("Source optimized dose contexts cannot be collapsed")
    if endpoint["name"] != endpoint_name or endpoint["unit"] != "fraction":
        raise PopulationError("Source endpoint or unit differs from its published readout")
    if context["measurement_timepoint"] != {"value": known(timepoint), "unit": unit}:
        raise PopulationError("Source measurement timepoint differs from the canonical context")


def _check_metadata_locations(name, cell):
    """A value cannot borrow valid metadata from another row or arm."""
    address = cell["locator"]["cell"]
    column = re.sub(r"\d", "", address)
    row = int(re.search(r"\d+$", address)[0])
    if name == "liang_endov_hbg":
        treated = column in "CDE"
        site = "HBG2-TA" if row == 6 else "HBG1-TA" if row == 7 else "HBG-OT" + str(row - 7)
        expected = {"arm_header_cell": "C4" if treated else "F4", "row_label_cell": "B" + str(row),
            "replicate_header_cell": column + "5", "source_replicate": "Data " + str(("CDE" if treated else "FGH").index(column) + 1),
            "source_site": site, "source_row_label": "HBG2 A8" if row == 6 else "HBG1 A8" if row == 7 else "OT" + str(row - 7),
            "source_subject": "HEK293T"}
    else:
        first, panel, header = (4, "Figure 4a", 2) if row <= 12 else (17, "Figure 4b", 16) if row <= 25 else (30, "Figure 4c", 29)
        group = "CDEFGHIJK".index(column) // 3 if panel.endswith("a") else "CDE".index(column)
        expected = {"source_panel": panel, "arm_header_cell": ("CFI"[group] if panel.endswith("a") else column) + str(header),
            "subject_cell": "A" + str(first + 3 * ((row - first) // 3)), "replicate_cell": "B" + str(row),
            "source_replicate": "rep#" + str((row - first) % 3 + 1),
            "source_position_cell": column + "3" if panel.endswith("a") else None,
            "guide_position_1based": (11, 5, 3)["CDEFGHIJK".index(column) % 3] if panel.endswith("a") else None}
    if any(cell[key] != value for key, value in expected.items()):
        raise PopulationError("Source metadata locations or replicate labels disagree with the measured cell")


def _normalization(name, row, binding, index):
    liang = name == "liang_endov_hbg"
    design_id = binding["design_id"]["value"]
    source_subject = "doi:" + SOURCES[name]["doi"] + ":" + row["source_subject"]
    independent = unknown("The cell-line replicate labels do not identify independent human donors.") if liang else known(source_subject)
    pairing = unknown("Cross-site/cross-arm cell-line replicate pairing is not independently established.") if liang else known(
        source_subject + ":" + row["source_panel"] + ":" + row["source_replicate"])
    loci = []
    if design_id is not None:
        intended = index[design_id]["intended_locus_ids"]
        if not liang:
            loci = intended
        elif row["source_site"].endswith("-TA"):
            gene = row["source_site"].removesuffix("-TA")
            loci = [identifier for identifier in intended if index[identifier]["gene"] == gene]
    outcome = known(row["source_reported_validation_label"]) if liang and design_id is not None else unknown(
        "The source treated-site validation label does not label the GFP control." if liang else "No categorical outcome label is reported for this numeric cell.")
    if design_id is None:
        adjudicated, adjudication_reason = "untested", ("The source control arm has no treatment design; "
            "the design editing outcome is untested in this source row.")
    elif liang:
        adjudicated, adjudication_reason = "inconclusive", ("The verbatim source validation label has no reviewed "
            "four-state adjudication; conclusive credit is withheld pending scientific review.")
    else:
        adjudicated, adjudication_reason = "inconclusive", ("No categorical outcome label is reported for this numeric "
            "cell; a measured number alone is not an adjudication, so conclusive credit is withheld.")
    return {
        "adjudicated_outcome": adjudicated, "adjudication_reason": adjudication_reason,
        "subject": {"level": "cell_line" if liang else "donor", "source_subject_id": known(source_subject),
            "donor_id": unknown("Source subject has no validated canonical genomic-panel donor linkage."),
            "specimen_id": unknown("Source specimen has no validated canonical research-inventory linkage."),
            "resolution_status": "unresolved", "reason": "Source-local identity is retained; no individual genomic or specimen identity is inferred."},
        "design_id": binding["design_id"],
        "unresolved_source_design": unknown("Source guide/editor combination is resolved.") if design_id else known(row["source_arm"] + " control; no treatment design"),
        "context_id": known(binding["context_id"]),
        "sequence_configuration_id": unknown("Source genotype and exact measured allele/path are not established; a reference target is not a measured donor sequence."),
        "locus_ids": loci, "endpoint_id": binding["endpoint_id"],
        "normalized_unit": "fraction", "conversion": known({"method": "Source percent divided by 100 exactly", "factor": 0.01, "offset": 0}),
        "measurement_kind": "individual", "independent_unit_id": independent,
        "relationship_group": unknown("Relationships among source participants are not supplied."), "pairing_key": pairing,
        "read_count": unknown("Numerator read/cell counts are not supplied for this source cell."),
        "denominator": unknown("Numeric denominator is not supplied for this source cell."),
        "sample_size": unknown("The cell is one source experiment measurement; no additional per-cell sample size is inferred."),
        "uncertainty": unknown("No uncertainty interval accompanies this individual source cell."),
        "assay_id": known("doi:" + SOURCES[name]["doi"] + ":" + ("Figure 2d deep sequencing" if liang else row["source_panel"])),
        "control_status": "unknown", "qc_status": "unknown",
        "detection_limit": {"value": unknown("Numeric detection limit is not reported for this source cell."), "unit": "fraction"},
        "quantification_limit": {"value": unknown("Numeric quantification limit is not reported for this source cell."), "unit": "fraction"},
        "source_reported_outcome": outcome,
        "executed_action_id": unknown("Historical publication measurement was not returned from a planned action in this study."),
    }


def build_source_map(source_name, catalogs, curated_cells, canonical_bindings):
    """Make an exact workbook map from reviewed cells and canonical identities."""
    source = _source(source_name)
    if not isinstance(curated_cells, list) or len(curated_cells) != source["count"]:
        raise PopulationError("Published source extraction count differs from its fixed contract")
    index = validate_catalog(catalogs)
    if not isinstance(canonical_bindings, list):
        raise PopulationError("Canonical source bindings must be an array")
    bindings = {}
    for binding in canonical_bindings:
        key = binding["source_record_key"]
        if key in bindings:
            raise PopulationError("Canonical source bindings repeat a measurement key")
        bindings[key] = binding
    expected = _layout(source_name)
    seen, rows, source_ids = set(), [], set()
    for cell in curated_cells:
        try:
            loc = cell["locator"]
            address = loc["cell"]
            if address in seen or address not in expected:
                raise PopulationError("Published source cell set is duplicate or outside its exact contract")
            seen.add(address)
            if (loc["artifact_role"], loc["artifact_sha256"], loc["sheet"], loc["worksheet_xml"]) != (
                    source["role"], source["sha256"], source["sheet"], source["member"]):
                raise PopulationError("Source fixture does not identify the exact published workbook/sheet")
            if cell["source_record_key"] != source["role"] + "|" + source["sheet"] + "|" + address:
                raise PopulationError("Published source cell and measurement key disagree")
            binding = bindings[cell["source_record_key"]]
            artifact_id = binding["source_location"]["source_artifact_id"]
            artifact = index.get(artifact_id)
            if artifact is None or artifact["record_type"] != "comparison_source_artifact" or artifact["sha256"] != source["sha256"]:
                raise PopulationError("Canonical artifact is not the pinned published workbook")
            source_ids.add(artifact_id)
            if not isinstance(cell["raw_value"], str) or cell["raw_unit"] != "percent" or binding["raw_value"] != cell["raw_value"] or binding["raw_unit"] != "percent":
                raise PopulationError("Source raw values or units disagree across reviewed mappings")
            value = Decimal(cell["raw_value"])
            if not value.is_finite() or not 0 <= value <= 100 or Decimal(cell["normalized_decimal"]) != value / 100:
                raise PopulationError("Source percent-to-fraction fixture is invalid")
            if cell["source_arm"].replace(" ", "") != expected[address]["arm"]:
                raise PopulationError("Treatment/control label disagrees with its source column")
            _check_metadata_locations(source_name, cell)
            _check_binding(source_name, cell, binding, expected[address], index)
            assertions = []
            def assert_cell(at, value):
                assertions.append({"locator": _locator(source, at), "expected_text": value})
            assert_cell(cell["arm_header_cell"], cell["source_arm"])
            if source_name == "liang_endov_hbg":
                assert_cell("A3", "Figure 2d")
                assert_cell(cell["row_label_cell"], cell["source_row_label"])
                assert_cell(cell["replicate_header_cell"], cell["source_replicate"])
                experiment = "Figure 2d:" + cell["source_site"] + ":" + cell["source_arm"]
            else:
                assert_cell({"Figure 4a": "A1", "Figure 4b": "A15", "Figure 4c": "A28"}[cell["source_panel"]], cell["source_panel"])
                assert_cell(cell["subject_cell"], cell["source_subject"])
                assert_cell(cell["replicate_cell"], cell["source_replicate"])
                if cell["source_position_cell"] is not None:
                    assert_cell(cell["source_position_cell"], "A" + str(cell["guide_position_1based"]))
                experiment = cell["source_panel"] + ":" + cell["source_subject"] + ":" + cell["source_replicate"]
            rows.append({"measurement_key": cell["source_record_key"], "locator": _locator(source, address),
                "expected_raw_value": known(cell["raw_value"]), "raw_unit": "percent",
                "source_study": "doi:" + source["doi"], "source_experiment": experiment,
                "source_replicate": cell["source_replicate"], "source_reported_label": known(cell["source_arm"]),
                "normalization": _normalization(source_name, cell, binding, index),
                "normalization_unavailable_reason": None, "assertions": assertions})
        except (KeyError, TypeError, InvalidOperation) as exc:
            raise PopulationError("Published source fixture or canonical binding is incomplete or invalid") from exc
    if set(expected) != seen or len(source_ids) != 1:
        raise PopulationError("Published source extraction is incomplete or crosses artifact identities")
    rows.sort(key=lambda row: (row["locator"]["row"], row["locator"]["column"]))
    return {"record_type": "comparison_extraction_map", "schema_version": "1", "source_artifact_id": next(iter(source_ids)),
        "artifact_format": "xlsx", "extraction_version": source_name + "-physical-cells-v1",
        "mapping_version": source_name + "-reviewed-outcomes-v2", "expected_count": source["count"], "rows": rows}


def import_published_source(source_name, workbook, map_path, catalogs, curated_cells, canonical_bindings):
    """Pin source bytes, preserve the map, and run the production observation importer."""
    source = _source(source_name)
    if checksum(workbook) != source["sha256"]:
        raise PopulationError("Published-source wrapper refuses a different workbook, regardless of row count")
    mapping = build_source_map(source_name, catalogs, curated_cells, canonical_bindings)
    path = Path(map_path)
    if path.exists():
        if read_json(path) != mapping:
            raise PopulationError("Preserve prior extraction maps; use a new revision path")
    else:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(mapping, ensure_ascii=False, indent=2, allow_nan=False) + "\n")
    return import_artifact(workbook, path, catalogs)
