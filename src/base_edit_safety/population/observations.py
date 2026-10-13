"""Import measured values from pinned files without inventing assay metadata."""
from __future__ import annotations

from copy import deepcopy
from decimal import Decimal, InvalidOperation, localcontext
import json
import re
from xml.etree import ElementTree as ET
from zipfile import ZipFile

from ..ingestion.xlsx import _NS, _shared_strings, _sheet_paths, _text
from ..measurement_units import editing_detection_limit
from .common import PopulationError, checksum, fields, integer, known, read_json, read_jsonl, read_tsv, unknown
from .comparison_contracts import make_record, record_schema, validate_catalog
from .comparison_outcomes import outcome_reasons


GENERIC_COLUMNS = ("measurement_key", "source_study", "source_experiment",
                   "source_replicate", "raw_value", "raw_unit", "source_reported_label")
_REVISION_COMPUTED = {"record_type", "schema_version", "record_id", "observation_id",
    "source_record_id", "extraction_version", "mapping_version", "normalization_schema_version",
    "extraction_map_sha256", "normalized_value", "applicability"}


def _decimal(value):
    if isinstance(value, bool):
        raise PopulationError("A boolean is not a measured number")
    try:
        number = Decimal(str(value))
    except InvalidOperation as exc:
        raise PopulationError("A numeric normalization requires a numeric source value") from exc
    if not number.is_finite():
        raise PopulationError("A measured number must be finite")
    return number


def _number(value):
    """Keep decimal text when a JSON float would change the measured value."""
    result = int(value) if value == value.to_integral() else float(value)
    if Decimal(str(result)) != value:
        return format(value, "f")
    return result


def _workbook(path):
    """Read physical cells, preserving gaps in row numbers and raw XML rows."""
    cells = {}
    with ZipFile(path) as book:
        shared = _shared_strings(book)
        for sheet, member in _sheet_paths(book).items():
            root = ET.fromstring(book.read(member))
            for row in root.findall("m:sheetData/m:row", _NS):
                raw_row = ET.tostring(row, encoding="unicode")
                for cell in row.findall("m:c", _NS):
                    address = cell.attrib["r"]
                    if (sheet, address) in cells:
                        raise PopulationError("Workbook repeats a physical cell address")
                    node = cell.find("m:v", _NS)
                    raw = "" if node is None or node.text is None else node.text
                    kind = cell.attrib.get("t")
                    if kind == "s":
                        value = shared[int(raw)]
                    elif kind == "inlineStr":
                        node = cell.find("m:is", _NS)
                        value = "" if node is None else _text(node)
                    elif kind == "b":
                        value = "true" if raw == "1" else "false"
                    else:
                        value = raw
                    cells[sheet, address] = value, raw_row, member
    return cells


def _generic_rows(path, format_name):
    if format_name == "tsv":
        rows = read_tsv(path, GENERIC_COLUMNS)
    else:
        rows = read_json(path)
        if not isinstance(rows, list):
            raise PopulationError("Generic observation JSON must be an array")
    result = {}
    for position, row in enumerate(rows, 1):
        fields(row, GENERIC_COLUMNS, label="generic source measurement")
        if any(not isinstance(row[key], str) or not row[key] for key in GENERIC_COLUMNS):
            raise PopulationError("Generic source fields must be nonempty strings, including numeric text")
        if row["measurement_key"] in result:
            raise PopulationError("Generic source repeats a measurement key")
        result[row["measurement_key"]] = row, position
    return result


def _locator_check(locator):
    from .contracts import validate_schema
    validate_schema(locator, record_schema("comparison_source_observation")["properties"]["locator"])


def _lookup_xlsx(cells, locator):
    _locator_check(locator)
    if locator["kind"] != "table":
        raise PopulationError("Workbook measurements require physical table/cell locators")
    match = re.fullmatch(r"([A-Z]+)([1-9][0-9]*)", locator["cell"])
    if not match or match[1] != locator["column"] or int(match[2]) != locator["row"]:
        raise PopulationError("Workbook row, column and cell locators disagree")
    try:
        value, row, member = cells[locator["sheet"], locator["cell"]]
    except KeyError as exc:
        raise PopulationError("Mapped physical workbook cell is missing") from exc
    if locator["table"] != member:
        raise PopulationError("Workbook sheet member differs from its extraction map")
    return value, row


def _ledger_rows(path):
    from ..pipeline.candidate_panel import compare_ledgers
    ledger = read_json(path)
    if ledger.get("record_type") != "candidate_panel_ledger" or ledger.get("schema_version") != "3":
        raise PopulationError("Observation bridge requires a candidate ledger schema 3")
    compare_ledgers(ledger, ledger)
    result = {}
    for site in ledger["candidates"]:
        for position, measurement in enumerate(site["measurements"]):
            member = f"candidates/{site['site_id']}/measurements/{position}"
            result[member] = measurement
    return result


def _resolved_mock(revision, index, source):
    """An explicit completed mock has a control identity without an editor design."""
    action = index.get(revision["executed_action_id"]["value"], {})
    if action.get("record_type") != "comparison_action" or action["operation_type"] != "mock_control" \
            or action["design_id"]["value"] is not None or "mock" not in action.get("performed_control_roles", []):
        return False
    subject = revision["subject"]
    if (subject["donor_id"]["value"], subject["specimen_id"]["value"], revision["context_id"]["value"]) != \
            (action["donor_id"], action["specimen_id"], action["context_id"]):
        return False
    source = source or index.get(revision.get("source_record_id"), {})
    if source.get("source_experiment") != action["experiment_id"]:
        return False
    pairs = action.get("endpoint_locus_pairs", [])
    sequence = revision["sequence_configuration_id"]["value"]
    return len(revision["locus_ids"]) == 1 and all(
        {"endpoint_id": revision["endpoint_id"], "locus_id": locus} in pairs and
        any(target["locus_id"] == locus and sequence in target["sequence_configuration_ids"]
            for target in action.get("target_sequences", [])) for locus in revision["locus_ids"])


def intrinsic_assessment(revision, index, *, source_observation=None):
    """Assess the source's own identity and quality, before any target-study reuse."""
    from .comparison_contracts import validate_observation_locus
    validate_observation_locus(revision, index)
    endpoint = index.get(revision["endpoint_id"])
    if endpoint is None or endpoint["record_type"] != "comparison_endpoint_definition":
        raise PopulationError("Observation endpoint is unresolved")
    identity = []
    if revision["design_id"]["value"] is None and not _resolved_mock(revision, index, source_observation):
        identity.append("source_design_unresolved")
    if endpoint["requires_sequence"] and revision["sequence_configuration_id"]["value"] is None:
        identity.append("source_sequence_unresolved")
    if endpoint["locus_resolution"] in {"copy_specific", "locus_specific"} and len(revision["locus_ids"]) != 1:
        identity.append("source_locus_unresolved")
    context = [] if revision["context_id"]["value"] is not None else ["source_context_unresolved"]
    quality = outcome_reasons(revision)
    if revision["qc_status"] != "passed":
        quality.append("source_qc_not_passed")
    if revision["control_status"] != "passed":
        quality.append("source_controls_not_passed")
    if revision["normalized_value"]["value"] is None:
        quality.append("measurement_unavailable_or_censored")
    value = revision["normalized_value"]["value"]
    if value is not None and Decimal(str(value)) == 0:
        limit = revision["detection_limit"]
        valid = editing_detection_limit({"value": limit["value"]["value"], "unit": limit["unit"]})["status"] == "accepted" \
            if revision["normalized_unit"] == "fraction" else limit["unit"] == revision["normalized_unit"] and (limit["value"]["value"] or 0) > 0
        if not valid:
            quality.append("zero_without_usable_detection_limit")
    if value is not None and revision["normalized_unit"] == "fraction" and not 0 <= Decimal(str(value)) <= 1 \
            and "measured_fraction_outside_bounds" not in quality:
        quality.append("measured_fraction_outside_bounds")
    subject = revision["subject"]
    donor_known = subject["level"] in {"donor", "specimen"} and subject["resolution_status"] == "resolved" \
        and subject["donor_id"]["value"] is not None
    exact = not identity and not context
    classification = ("exact_donor_configuration" if donor_known else "exact_sequence_context") \
        if exact else "unresolved"
    def dimension(reasons):
        return {"status": "unresolved" if reasons else "applicable", "reason_codes": reasons}
    reasons = identity + context + quality
    if not donor_known:
        reasons.append("no_resolved_individual_donor_linkage")
    return {"identity": dimension(identity), "context": dimension(context),
            "endpoint": dimension([]), "quality": dimension(quality),
            "classification": classification,
            "empirical_donor_credit": donor_known and exact and not quality,
            "reason_codes": sorted(set(reasons))}


def import_artifact(path, extraction_map_path, catalogs):
    """Verify source bytes and locators, then create one revision per measured cell."""
    mapping = read_json(extraction_map_path)
    fields(mapping, {"record_type", "schema_version", "source_artifact_id", "artifact_format",
        "extraction_version", "mapping_version", "expected_count", "rows"}, label="extraction map")
    if mapping["record_type"] != "comparison_extraction_map" or mapping["schema_version"] != "1":
        raise PopulationError("Unsupported observation extraction map")
    if mapping["artifact_format"] not in {"xlsx", "tsv", "json", "candidate_ledger"}:
        raise PopulationError("Unsupported observation artifact format")
    integer(mapping["expected_count"], "expected measurement count", minimum=1)
    if not isinstance(mapping["rows"], list) or len(mapping["rows"]) != mapping["expected_count"]:
        raise PopulationError("Extraction count differs from its declared contract")
    index = validate_catalog(catalogs)
    source = index.get(mapping["source_artifact_id"])
    if source is None or source["record_type"] != "comparison_source_artifact" or checksum(path) != source["sha256"]:
        raise PopulationError("Observation artifact does not match its pinned source")
    format_name = mapping["artifact_format"]
    raw = _workbook(path) if format_name == "xlsx" else _ledger_rows(path) \
        if format_name == "candidate_ledger" else _generic_rows(path, format_name)
    map_hash, output, seen, locations = checksum(extraction_map_path), [], set(), set()
    normalization_fields = set(record_schema("comparison_observation_revision")["required"]) - _REVISION_COMPUTED
    for entry in mapping["rows"]:
        fields(entry, {"measurement_key", "locator", "expected_raw_value", "raw_unit", "source_study",
            "source_experiment", "source_replicate", "source_reported_label", "normalization",
            "normalization_unavailable_reason", "assertions"}, label="extraction row")
        _locator_check(entry["locator"])
        from .contracts import validate_schema
        for key in ("expected_raw_value", "source_reported_label"):
            schema_key = "raw_value" if key == "expected_raw_value" else key
            validate_schema(entry[key], record_schema("comparison_source_observation")["properties"][schema_key])
        if entry["measurement_key"] in seen:
            raise PopulationError("Extraction map repeats a measurement key")
        seen.add(entry["measurement_key"])
        location = json.dumps(entry["locator"], sort_keys=True)
        if location in locations:
            raise PopulationError("Extraction map counts the same source cell twice")
        locations.add(location)
        if not isinstance(entry["assertions"], list):
            raise PopulationError("Source assertions must be an array")
        for assertion in entry["assertions"]:
            fields(assertion, {"locator", "expected_text"}, label="source cell assertion")
            if format_name != "xlsx" or _lookup_xlsx(raw, assertion["locator"])[0] != assertion["expected_text"]:
                raise PopulationError("Source metadata cell differs from its extraction map")
        if format_name == "xlsx":
            value, raw_row = _lookup_xlsx(raw, entry["locator"])
        elif format_name == "candidate_ledger":
            locator = entry["locator"]
            if locator["kind"] != "file" or locator["member"] not in raw:
                raise PopulationError("Ledger measurement locator is unresolved")
            measurement = raw[locator["member"]]
            raw_unit = measurement.get("unit")
            if raw_unit is None and measurement.get("value") is None:
                raw_unit = "unreported in source ledger"
            if entry["raw_unit"] != raw_unit or entry["source_reported_label"]["value"] != measurement.get("status"):
                raise PopulationError("Ledger unit or reported outcome differs from its extraction map")
            if entry["source_experiment"] != measurement.get("assay") or entry["source_study"] != measurement.get("source_version_or_citation"):
                raise PopulationError("Ledger source identity differs from its extraction map")
            if entry["source_replicate"] != measurement.get("replicate_id"):
                raise PopulationError("Ledger replicate identity differs from its extraction map")
            value = measurement.get("value")
            if isinstance(value, dict):
                value = value.get("value")
            value = None if value is None else str(value)
            raw_row = json.dumps(measurement, sort_keys=True, ensure_ascii=False, allow_nan=False)
        else:
            found = raw.get(entry["measurement_key"])
            if found is None:
                raise PopulationError("Mapped source measurement is missing")
            row, position = found
            locator = entry["locator"]
            if locator["kind"] != "file" or locator["member"] != "rows/" + str(position):
                raise PopulationError("Generic source row locator is inconsistent")
            for key in ("source_study", "source_experiment", "source_replicate", "raw_unit"):
                if entry[key] != row[key]:
                    raise PopulationError("Generic source metadata differs from its extraction map")
            if entry["source_reported_label"]["value"] != row["source_reported_label"]:
                raise PopulationError("Generic source label differs from its extraction map")
            value = row["raw_value"]
            raw_row = json.dumps(row, sort_keys=True, ensure_ascii=False, allow_nan=False)
        if entry["expected_raw_value"]["value"] != value:
            raise PopulationError("Raw source value differs from its extraction map")
        source_row = make_record("source_observation", source_artifact_id=source["record_id"],
            locator=entry["locator"], measurement_key=entry["measurement_key"],
            source_study=entry["source_study"], source_experiment=entry["source_experiment"],
            source_replicate=entry["source_replicate"], raw_row_text=raw_row,
            raw_value=entry["expected_raw_value"], raw_unit=entry["raw_unit"],
            source_reported_label=entry["source_reported_label"])
        normalization = deepcopy(entry["normalization"])
        explicit_outcome = "adjudicated_outcome" in normalization or "adjudication_reason" in normalization
        required = normalization_fields | {"adjudicated_outcome", "adjudication_reason"} if explicit_outcome else normalization_fields
        fields(normalization, required, label="observation normalization")
        if not explicit_outcome:
            normalization.update(adjudicated_outcome="inconclusive",
                adjudication_reason="Legacy extraction map has no explicit outcome adjudication; conclusive credit withheld")
        schema = record_schema("comparison_observation_revision", "2")["properties"]
        for key, item in normalization.items():
            validate_schema(item, schema[key], "normalization." + key)
        reason = entry["normalization_unavailable_reason"]
        if reason is not None:
            if not isinstance(reason, str) or not reason.strip():
                raise PopulationError("An unavailable normalization requires a reason")
            normalized = unknown(reason)
        else:
            number = _decimal(value)
            conversion = normalization["conversion"]["value"]
            if conversion is not None:
                factor, offset = _decimal(conversion["factor"]), _decimal(conversion["offset"])
                with localcontext() as context:
                    context.prec = max(50, sum(len(n.as_tuple().digits) + abs(n.as_tuple().exponent)
                                               for n in (number, factor, offset)) + 4)
                    number = number * factor + offset
            normalized = known(_number(number))
        normalization["normalized_value"] = normalized
        normalization["applicability"] = intrinsic_assessment(normalization, index, source_observation=source_row)
        revision = make_record("observation_revision", **normalization,
            source_record_id=source_row["record_id"], extraction_version=mapping["extraction_version"],
            mapping_version=mapping["mapping_version"], normalization_schema_version="2",
            extraction_map_sha256=map_hash)
        state = make_record("observation_state", source_record_id=source_row["record_id"],
            active_revision_id=revision["record_id"], revision_ids=[revision["record_id"]],
            supersession_reason="Initial normalization from a verified source file")
        output.extend((source_row, revision, state))
    if format_name in {"json", "tsv"} and set(raw) != seen:
        raise PopulationError("Generic source contains unmapped measurements")
    validate_catalog(catalogs + output)
    return output


def merge_history(current, previous):
    """Retain earlier revisions while replacing each source's active state."""
    states, rows = {}, {}
    for collection in (previous, current):
        for row in collection:
            if row["record_type"] == "comparison_observation_state":
                old = states.get(row["source_record_id"])
                revisions = sorted(set(row["revision_ids"] + (old["revision_ids"] if old else [])))
                states[row["source_record_id"]] = make_record("observation_state",
                    source_record_id=row["source_record_id"], active_revision_id=row["active_revision_id"],
                    revision_ids=revisions, supersession_reason="Latest declared normalization; earlier revisions retained")
            else:
                old = rows.get(row["record_id"])
                if old is not None and old != row:
                    raise PopulationError("A source or normalization identity has conflicting content")
                rows[row["record_id"]] = row
    return sorted([*rows.values(), *states.values()], key=lambda row: (row["record_type"], row["record_id"]))


def normalize_observations(study, previous=()):
    """Import every declared extraction map, keeping source and revision history."""
    from .comparison_provenance import verify_sources
    verify_sources(study)
    history = list(previous)
    for source in sorted(study["sources"].values(), key=lambda row: row["input_id"]):
        if source["role"] == "observation_history" or (
                source["role"] == "raw_observations" and source["format"] == "jsonl"):
            if source["format"] != "jsonl":
                raise PopulationError("Observation history requires canonical JSONL records")
            history.extend(read_jsonl(source["path"]))
    catalog = study["config"]["catalogs"]
    source_files = {source["sha256"]: source["path"] for source in study["sources"].values()
                    if source["role"] == "raw_observations"}
    index = {row["record_id"]: row for row in catalog}
    current = []
    for source in sorted(study["sources"].values(), key=lambda row: row["input_id"]):
        if source["role"] != "extraction_map":
            continue
        mapping = read_json(source["path"])
        artifact = index.get(mapping.get("source_artifact_id"), {})
        path = source_files.get(artifact.get("sha256"))
        if path is None:
            raise PopulationError("Extraction map has no declared raw observation input")
        current.extend(import_artifact(path, source["path"], catalog))
    rows = merge_history(current, history)
    validate_catalog(catalog + rows)
    return {"observations.jsonl": rows,
            "normalization-summary.json": {"source_measurements": sum(row["record_type"] == "comparison_source_observation" for row in rows),
                "active_revisions": sum(row["record_type"] == "comparison_observation_state" for row in rows),
                "retained_revisions": sum(row["record_type"] == "comparison_observation_revision" for row in rows)}}
