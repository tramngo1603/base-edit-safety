"""Strict records and references for comparative studies.

The registry is separate from the legacy population, ledger and design schemas.
It validates stored records, including their content IDs, before scientific use.
"""
from __future__ import annotations

from copy import deepcopy
from decimal import Decimal, InvalidOperation
from functools import lru_cache
from importlib.resources import files
import hashlib
import re

from .common import PopulationError, parse_json
from .contracts import validate_basis, validate_schema
from .comparison_identity import canonical_bytes, content_id
from .comparison_outcomes import outcome_reasons

INPUT_ROLE_STAGES = {
    "reference": ("search",),
    "reference_validation": ("search",),
    "variants": ("search",),
    "donor_map": ("search",), "callability": ("search",),
    "source_transformation": ("search",),
    "raw_observations": ("observation_normalization",),
    "extraction_map": ("observation_normalization",),
    "annotations": ("comparison",), "prediction_inputs": ("predictions",),
    "specimen_inventory": ("planning",),
    "published_alignments": ("evidence_applicability",),
    "experiment_policy": ("questions",),
    "question_inventory": ("questions", "planning"),
    "plan_history": ("planning",),
    "observation_history": ("observation_normalization",),
}


@lru_cache(maxsize=3)
def _schema(name):
    return parse_json(files(__package__).joinpath("schemas", f"comparison-{name}.schema.json").read_text())


def _expand(schema, definitions):
    if "$ref" in schema:
        return _expand(definitions[schema["$ref"].removeprefix("#/$defs/")], definitions)
    return {key: (_expand(value, definitions) if isinstance(value, dict) else
                  [_expand(item, definitions) if isinstance(item, dict) else item for item in value]
                  if isinstance(value, list) else value)
            for key, value in schema.items() if key != "$defs"}


@lru_cache(maxsize=None)
def record_schema(record_type, version="1"):
    if record_type == "comparison_observation_revision" and version == "2":
        return _schema("observation-v2")
    if version != "1" or record_type not in record_registry():
        raise PopulationError("Unknown comparison record type/version")
    source = _schema("record")
    return _expand(source["$defs"][record_type], source["$defs"])


@lru_cache(maxsize=1)
def record_registry():
    """Map version-1 record types to their primary ID field."""
    source = _schema("record")
    return {item["$ref"].removeprefix("#/$defs/"):
            source["$defs"][item["$ref"].removeprefix("#/$defs/")]["x-id-field"]
            for item in source["oneOf"]}


@lru_cache(maxsize=1)
def study_schema():
    source = deepcopy(_schema("study"))
    # Catalog records dispatch by type, avoiding a scan of every record schema.
    source["properties"]["catalogs"]["items"] = {"type": "object"}
    return _expand(source, {})


def _option(value, schema):
    if "oneOf" not in schema:
        return schema
    valid = []
    for option in schema["oneOf"]:
        try:
            validate_schema(value, option)
            valid.append(option)
        except PopulationError:
            continue
    if len(valid) != 1:
        raise PopulationError("Value does not match exactly one comparative contract")
    return valid[0]


def _walk(value, schema, visit, location="record"):
    schema = _option(value, schema)
    visit(value, schema, location)
    if isinstance(value, dict):
        for key, item in value.items():
            if key in schema.get("properties", {}):
                _walk(item, schema["properties"][key], visit, location + "." + key)
    elif isinstance(value, list):
        for index, item in enumerate(value):
            _walk(item, schema["items"], visit, f"{location}[{index}]")


def _check_value(value, schema, location):
    if isinstance(value, str) and "enum" not in schema and "const" not in schema:
        if any(part in location for part in (".raw_row_text", ".raw_value", ".source_reported_label", ".source_reported_outcome", ".legacy_payload_json")):
            return
        if not value.strip() or re.search(r"\b(?:REPLACE(?:_[A-Z0-9]+)*|TODO|TBD|PLACEHOLDER)\b", value, re.I):
            raise PopulationError(location + ": unresolved placeholder or empty text")
        if value.strip().lower() in {"unknown", "n/a", "na", "none", "null"}:
            raise PopulationError(location + ": use null with a reason for an unknown value")
    if isinstance(value, dict):
        if {"contig", "start", "end"} <= value.keys() and value["start"] >= value["end"]:
            raise PopulationError(location + ": interval must be nonempty and half open")
        if {"alternates", "selected_alt_index"} <= value.keys():
            if value["selected_alt_index"] > len(value["alternates"]):
                raise PopulationError(location + ": ALT index exceeds the ordered source alleles")
    if isinstance(value, list) and schema.get("x-set"):
        if len({canonical_bytes(item) for item in value}) != len(value):
            raise PopulationError(location + ": duplicate value after canonical normalization")


def _canonical_structure(value, schema):
    schema = _option(value, schema)
    if isinstance(value, dict):
        return {key: _canonical_structure(item, schema["properties"][key])
                for key, item in value.items()}
    if isinstance(value, list):
        values = [_canonical_structure(item, schema["items"]) for item in value]
        return sorted(values, key=canonical_bytes) if schema.get("x-set") else values
    return value


def identity_payload(record):
    schema = record_schema(record["record_type"], record["schema_version"])
    payload = {key: record[key] for key in schema["x-identity-fields"]
               if key in record or key in schema["required"]}
    return _canonical_structure(payload, schema)


def record_identity(record):
    return content_id(record["record_type"], identity_payload(record), record["schema_version"])


def make_record(kind, *, version=None, **values):
    """Build and validate one canonical record; validate_catalog checks its links."""
    kind = kind if kind.startswith("comparison_") else "comparison_" + kind
    version = version or ("2" if kind == "comparison_observation_revision" and "adjudicated_outcome" in values else "1")
    schema = record_schema(kind, version)
    reserved = {"record_type", "schema_version", "record_id", schema["x-id-field"]}
    if reserved & values.keys():
        raise PopulationError("Canonical ID/envelope fields are assigned by make_record")
    record = {"record_type": kind, "schema_version": version, **deepcopy(values)}
    try:
        identifier = record_identity(record)
    except (KeyError, TypeError) as exc:
        raise PopulationError("Missing or invalid identity field: " + str(exc)) from exc
    record["record_id"] = record[schema["x-id-field"]] = identifier
    if kind == "comparison_design":
        if "definition_hash" in values and values["definition_hash"] != identifier.rsplit(":", 1)[1]:
            raise PopulationError("Design definition hash does not match its biology")
        record["definition_hash"] = identifier.rsplit(":", 1)[1]
    return validate_record(record)


def _path(sequence, path, location):
    if "".join(point["base"] for point in path) != sequence:
        raise PopulationError(location + ": sequence differs from its coordinate path")
    coords = [{key: item for key, item in point.items() if key != "base"} for point in path]
    if len({canonical_bytes(point) for point in coords}) != len(coords):
        raise PopulationError(location + ": path repeats a coordinate")


def _assessment(value):
    for field in ("identity", "context", "endpoint", "quality"):
        part = value[field]
        if part["status"] != "applicable" and not part["reason_codes"]:
            raise PopulationError("A nonapplicable evidence dimension requires a reason")
    if value["empirical_donor_credit"] and (value["classification"] != "exact_donor_configuration"
            or any(value[key]["status"] != "applicable" for key in ("identity", "context", "endpoint", "quality"))):
        raise PopulationError("Empirical credit requires exact, applicable donor evidence")


def _record_semantics(row):
    kind = row["record_type"].removeprefix("comparison_")
    if kind == "reference_bundle":
        complete = row["identity_kind"] == "complete_reference_files"
        if complete != (row["complete_reference_sha256"]["value"] is not None) or complete == bool(row["windows"]):
            raise PopulationError("Whole-reference and indexed-window identities must remain distinct")
    elif kind == "editor_definition":
        geometry, window = row["guide_geometry"], row["editing_window"]
        if "start_1based" in window and not 1 <= window["start_1based"] <= window["end_1based"] <= geometry["spacer_length"]:
            raise PopulationError("Editor editing window is outside the guide geometry")
        if any(len(pam) != geometry["pam_length"] for pam in row["pam_patterns"]):
            raise PopulationError("Editor PAM pattern length differs from its geometry")
        if row["conversion"]["substrate_base"] == row["conversion"]["product_base"]:
            raise PopulationError("An editor must declare a base conversion")
    elif kind == "design":
        if row["definition_hash"] != row["record_id"].rsplit(":", 1)[1]:
            raise PopulationError("Design definition hash differs from its identity")
        if any(edit["locus_id"] not in row["intended_locus_ids"] for edit in row["desired_edit_definitions"]):
            raise PopulationError("Desired edit is outside the declared intended loci")
        if row["source_status"] == "published_verified" and not row["source_refs"]:
            raise PopulationError("A verified published design requires primary source references")
    elif kind == "sequence_configuration":
        _path(row["sequence"], row["path"], kind)
        if row["support_status"] in {"modeled", "putative", "unresolved"} and row["support_reason"] is None:
            raise PopulationError("Uncertain sequence support requires a reason")
    elif kind == "alignment":
        _path(row["pam_sequence"], row["pam_path"], kind)
        guide_position, mismatches, dna, rna = 0, [], [], []
        target_indexes = []
        for index, col in enumerate(row["columns"]):
            guide, target = col["guide"], col["target"]
            if guide is None and target is None:
                raise PopulationError("An alignment column cannot contain two gaps")
            if (guide is None) != (col["guide_position_1based"] is None) or (target is None) != (col["target_path_index"] is None):
                raise PopulationError("Alignment gaps and positions disagree")
            if guide is not None:
                guide_position += 1
                if col["guide_position_1based"] != guide_position:
                    raise PopulationError("Guide alignment positions must preserve sequence order")
            if target is not None:
                target_indexes.append(col["target_path_index"])
            if guide is None:
                dna.append(index)
            elif target is None:
                rna.append(guide_position)
            elif guide != target:
                mismatches.append(guide_position)
        if len(set(target_indexes)) != len(target_indexes) or target_indexes != sorted(target_indexes):
            raise PopulationError("Alignment target path indexes must be distinct and ordered")
        if set(mismatches) != set(row["mismatch_positions"]) or set(dna) != set(row["dna_bulge_positions"]) or set(rna) != set(row["rna_bulge_positions"]):
            raise PopulationError("Alignment mismatch/bulge annotations disagree with its exact columns")
    elif kind == "donor":
        if row["relatedness_policy"] == "declared_simulation_assumption" and row["independence_group_id"]["value"] is None:
            raise PopulationError("A simulation independence assumption must name its group")
    elif kind == "donor_configuration":
        status = row["carrier_status"]
        if status.startswith("supported_") and row["genotype_support"] not in {"called", "callable_reference"}:
            raise PopulationError("Missing or modeled genotypes cannot establish supported carriage")
        if status.startswith("modeled_") and row["genotype_support"] != "retained_callset_model":
            raise PopulationError("Modeled carriage requires its declared retained-callset basis")
        if status in {"unresolved", "not_evaluated"} and not row["reason_codes"]:
            raise PopulationError("Unresolved carriage requires a reason")
    elif kind == "search_scope":
        if row["execution_status"] == "complete" and (row["skipped_intervals"] or row["retention_status"] != "complete" or row["exit_code"]["value"] != 0):
            raise PopulationError("Skipped, truncated or failed search cannot be complete")
        if row["execution_status"] != "complete" and not row["reason_codes"]:
            raise PopulationError("Incomplete search requires a reason")
    elif kind == "observation_revision":
        _assessment(row["applicability"])
        if row["design_id"]["value"] is None and row["unresolved_source_design"]["value"] is None:
            raise PopulationError("Unresolved design must retain its source label")
        subject = row["subject"]
        if subject["level"] in {"pooled_group", "study_summary", "cell_line"} and subject["donor_id"]["value"] is not None:
            raise PopulationError("A pooled, summary or cell-line observation cannot acquire an individual donor")
        count, denominator = row["read_count"]["value"], row["denominator"]["value"]
        if count is not None and denominator is not None and count > denominator:
            raise PopulationError("Read count exceeds its measured denominator")
        if row["applicability"]["empirical_donor_credit"]:
            if subject["donor_id"]["value"] is None or row["control_status"] != "passed" or row["qc_status"] != "passed":
                raise PopulationError("Empirical donor credit requires an identified donor and passed controls/QC")
            if row["normalized_value"]["value"] is None:
                raise PopulationError("Missing measurement cannot receive empirical credit")
            if Decimal(str(row["normalized_value"]["value"])) == 0 and not ((row["detection_limit"]["value"]["value"] or 0) > 0):
                raise PopulationError("Measured zero without a usable detection limit earns no reassuring credit")
        if row["schema_version"] == "2" and row["applicability"]["empirical_donor_credit"] and outcome_reasons(row):
            raise PopulationError("Conclusive empirical credit requires a consistent adjudicated outcome: " + ", ".join(outcome_reasons(row)))
    elif kind == "evidence_applicability":
        _assessment(row["assessment"])
        if row["coverage_basis"] == "retained_callset_model" and row["assessment"]["empirical_donor_credit"]:
            raise PopulationError("Retained-callset simulations receive zero empirical donor credit")
    elif kind == "prediction":
        _path(row["window"]["sequence"], row["window"]["path"], "prediction window")
    elif kind == "question":
        if not row["requirement_slot_ids"]:
            reason = "requirement_construction_shortlist_exclusion" if row.get("construction_status") == "excluded_by_declared_limit" else "no_eligible_declared_donor"
            if row["planning_status"] != "blocked" or reason not in row["remaining_gaps"]:
                raise PopulationError("A question without slots must preserve its explicit construction or donor eligibility gap")
        if row.get("construction_status") == "excluded_by_declared_limit" and row["requirement_slot_ids"]:
            raise PopulationError("An excluded question cannot contain constructed requirement slots")
        if row["template"] in {"candidate_validation", "variant_association", "evidence_replication"} and len(row["design_ids"]) != 1:
            raise PopulationError("This question template requires exactly one design")
        if row["template"] in {"paired_design_comparison", "design_by_variant"}:
            if len(row["design_ids"]) != 2 or row["pairing"] != "within_donor":
                raise PopulationError("Design comparisons require two designs paired within the same donor")
        if row["template"] in {"variant_association", "design_by_variant"}:
            if {arm["membership"] for arm in row["arms"]} != {"carrier", "noncarrier"}:
                raise PopulationError("Genotype contrasts require explicit carrier and noncarrier arms")
            if row["genotype_predicate"]["kind"] == "any_resolved":
                raise PopulationError("A genotype contrast requires a literal allele or configuration predicate")
        if len({arm["arm_id"] for arm in row["arms"]}) != len(row["arms"]):
            raise PopulationError("Question arm IDs must be distinct")
        if row["template"] == "evidence_replication":
            policy = row["evidence_reuse"]
            if not policy["require_new"] or not set(row["motivating_observation_ids"]) <= set(policy["excluded_observation_ids"]):
                raise PopulationError("Replication requires an explicit new dimension and excludes motivating measurements")
    elif kind == "resource":
        if row["instance_index"] >= row["max_pool_instances"]:
            raise PopulationError("Resource instance index exceeds its declared finite pool")
        dimensions = [value["dimension"] for value in row["capacities"]]
        if len(dimensions) != len(set(dimensions)):
            raise PopulationError("Resource capacity dimensions must be distinct")
    elif kind == "action":
        if row["eligibility_status"] != "eligible" and not row["reason_codes"]:
            raise PopulationError("Ineligible or unresolved action requires a reason")
        if row["operation_type"] in {"locus_assay", "intended_assay"} and row["eligibility_status"] == "eligible" and row["target_sequence_id"]["value"] is None:
            raise PopulationError("A locus assay requires the donor's actual target sequence")
        if row["operation_type"] == "measurement_panel" and any(not row.get(field) for field in ("endpoint_locus_pairs", "target_sequences", "performed_qc")):
            raise PopulationError("A measurement panel requires exact endpoint/locus coverage, donor sequences and priced QC operations")
    elif kind == "experiment_plan":
        if row["weighted_completion"] != row["existing_evidence_weight"] + row["newly_planned_weight"]:
            raise PopulationError("Plan completion weights do not preserve existing and newly planned questions")
        if row["optimality_established"] != (row["status"] == "optimal"):
            raise PopulationError("Only an exact optimal plan can establish optimality")
        if row["meets_request"] and (row["total_cost"] > row["budget"] or row["status"] in {"infeasible_proven", "no_feasible_bundle_found", "invalid_plan"}):
            raise PopulationError("An infeasible or over-budget plan cannot meet the request")
        if len({state["question_id"] for state in row["question_states"]}) != len(row["question_states"]):
            raise PopulationError("Plan has duplicate question states")
        if len({(item["action_id"], item["assignment_key"]) for item in row["resource_assignments"]}) != len(row["resource_assignments"]):
            raise PopulationError("Plan assigns more than one instance to one action assignment key")
    elif kind == "identity_bridge":
        versions = {"population": "1", "candidate_ledger": "3", "design": "7"}
        if versions[row["legacy_family"]] != row["legacy_schema_version"]:
            raise PopulationError("Unsupported legacy family/version bridge")
        legacy = parse_json(row["legacy_payload_json"])
        if hashlib.sha256(canonical_bytes(legacy)).hexdigest() != row["legacy_payload_sha256"]:
            raise PopulationError("Legacy bridge payload no longer matches its recorded checksum")
    elif kind == "design_summary":
        if row["intended_evaluable_donors"] + row["intended_unresolved_donors"] != row["donor_count"] or row["search_evaluable_donors"] + row["search_unresolved_donors"] != row["donor_count"]:
            raise PopulationError("Summary evaluable and unknown donor counts must preserve the denominator")
        for field in ("supported_carrier_donors", "modeled_carrier_donors"):
            if row[field] > row["donor_count"]:
                raise PopulationError("Summary carrier count exceeds its donor denominator")
        for field in ("locus_count", "sequence_configuration_count", "alignment_count"):
            complete = row["complete_" + field]["value"]
            if complete is not None and complete != row["observed_" + field]:
                raise PopulationError("Complete and retained observed counts disagree")
            if row["search_unresolved_donors"] and complete is not None:
                raise PopulationError("Partial search totals must remain unknown, including apparent zero")
    elif kind == "donor_design_summary":
        if row["observed_configuration_count"] != len(row["nominated_configuration_ids"]) or row["observed_alignment_count"] != len(row["nominated_alignment_ids"]):
            raise PopulationError("Donor summary counts disagree with their exact retained records")
        complete = row["complete_configuration_count"]["value"]
        if (row["search_status"] == "complete") != (complete is not None) or complete is not None and complete != row["observed_configuration_count"]:
            raise PopulationError("Donor summary must distinguish complete zero from incomplete search")
        if bool(row["evidence_observation_ids"]) != (row["evidence_availability"] in {"contextual_only", "exact_donor_evidence"}):
            raise PopulationError("Donor evidence availability disagrees with its empirical observation references")
    elif kind == "decision":
        if row["state"] == "observed_difference_in_matched_context" and not row["observation_ids"]:
            raise PopulationError("An observed difference requires empirical observation references")
    if kind in {"question", "action", "assay_capability"} and "endpoint_locus_pairs" in row:
        pairs = row["endpoint_locus_pairs"]
        if {pair["endpoint_id"] for pair in pairs} != set(row["endpoint_ids"]) or {pair["locus_id"] for pair in pairs} != set(row["locus_ids"]):
            raise PopulationError("Exact endpoint/locus pairs disagree with declared endpoint/locus sets")


def validate_record(value):
    if not isinstance(value, dict):
        raise PopulationError("Comparison record must be an object")
    if not isinstance(value.get("record_type"), str) or not isinstance(value.get("schema_version"), str):
        raise PopulationError("Comparison records require string record_type and schema_version")
    schema = record_schema(value.get("record_type"), value.get("schema_version"))
    validate_schema(value, schema, value["record_type"])
    _walk(value, schema, _check_value)
    if value["record_id"] != value[schema["x-id-field"]] or value["record_id"] != record_identity(value):
        raise PopulationError("Record content and namespaced identity disagree")
    _record_semantics(value)
    return value


def validate_observation_locus(observation, index):
    """A resolved measured sequence cannot belong to a different reported locus."""
    identifier = observation["sequence_configuration_id"]["value"]
    if identifier is None:
        return
    sequence = index.get(identifier)
    if sequence is None or sequence["record_type"] != "comparison_sequence_configuration":
        raise PopulationError("Observation sequence reference is unresolved")
    if observation["locus_ids"] and sequence["locus_id"] not in observation["locus_ids"]:
        raise PopulationError("Observation sequence and reported locus disagree")


def _foreign_keys(row, schema, index):
    def check(value, subschema, location):
        if "x-record-types" in subschema:
            target = index.get(value)
            if target is None or target["record_type"] not in subschema["x-record-types"]:
                raise PopulationError(location + ": broken or wrong-type foreign key")
    _walk(row, schema, check)


def _acyclic(index, kind, field):
    visiting, visited = set(), set()
    def visit(identifier):
        if identifier in visiting:
            raise PopulationError("Cyclic " + kind + " dependencies")
        if identifier in visited:
            return
        visiting.add(identifier)
        values = index[identifier][field]
        if isinstance(values, dict):
            values = [] if values["value"] is None else [values["value"]]
        for target in values:
            visit(target)
        visiting.remove(identifier)
        visited.add(identifier)
    for identifier, row in index.items():
        if row["record_type"] == "comparison_" + kind:
            visit(identifier)


def resource_compatible(action, resource, index):
    """Return exact sharing/scope conflicts without assigning resource capacity."""
    reasons = []
    for field, actual in (("donor_ids", action["donor_id"]), ("design_ids", action["design_id"]["value"]),
                          ("context_ids", action["context_id"]), ("batch_ids", action["batch_id"]),
                          ("assay_capability_ids", action["assay_capability_id"])):
        allowed = resource["compatibility"][field]
        if allowed and actual not in allowed:
            reasons.append("Action violates its resource sharing compatibility: " + field)
    for field in ("endpoint_ids", "locus_ids"):
        allowed = resource["compatibility"][field]
        if allowed and set(action[field]) - set(allowed):
            reasons.append("Action exceeds its resource's declared endpoint/locus scope: " + field)
    allowed_times = resource["compatibility"]["timepoints"]
    if allowed_times:
        timepoint = index[action["context_id"]]["measurement_timepoint"]
        value = timepoint["value"]["value"]
        label = None if value is None else canonical_bytes(value).decode("utf-8") + " " + timepoint["unit"]
        if label not in allowed_times:
            reasons.append("Action timepoint is incompatible with its resource sharing rule")
    return reasons


def _catalog_semantics(index):
    get = index.__getitem__
    states, revisions, equivalences, resource_pools = {}, {}, {}, {}
    for row in index.values():
        kind = row["record_type"].removeprefix("comparison_")
        if kind == "reference_bundle":
            if set(row["source_checksums"]) != {get(identifier)["sha256"] for identifier in row["source_artifact_ids"]}:
                raise PopulationError("Reference content checksums differ from their source artifacts")
            if row["identity_kind"] == "indexed_windows" and set(row["source_checksums"]) != {window["source_sha256"] for window in row["windows"]}:
                raise PopulationError("Each indexed window must bind its actual regional source checksum")
            for identifier in row["source_artifact_ids"]:
                source = get(identifier)
                expected = "complete_file" if row["identity_kind"] == "complete_reference_files" else "indexed_window"
                if source["scope"] != expected:
                    raise PopulationError("Reference bundle uses a different source checksum scope")
            if row["identity_kind"] == "complete_reference_files" and len(row["source_artifact_ids"]) == 1 and row["complete_reference_sha256"]["value"] != get(row["source_artifact_ids"][0])["sha256"]:
                raise PopulationError("Complete reference checksum differs from its actual source file")
        elif kind == "design":
            editor = get(row["editor_definition_id"])
            if len(row["guide_sequence"]) != editor["guide_geometry"]["spacer_length"]:
                raise PopulationError("Design guide differs from editor geometry")
            for edit in row["desired_edit_definitions"]:
                locus = get(edit["locus_id"])
                if not locus["start"] <= edit["pos0"] < locus["end"]:
                    raise PopulationError("Desired edit coordinate is outside its intended locus")
        elif kind == "search_profile":
            editor = get(get(row["design_id"])["editor_definition_id"])
            if set(row["pam_patterns"]) - set(editor["pam_patterns"]):
                raise PopulationError("Search PAM policy is not declared by this design's editor")
        elif kind == "sequence_configuration":
            locus = get(row["locus_id"])
            if (row["contig"], row["reference_bundle_id"]) != (locus["contig"], locus["reference_bundle_id"]):
                raise PopulationError("Sequence content does not use its declared locus/reference")
        elif kind == "sequence_occurrence":
            if row["sample_id"] not in get(row["donor_id"])["sample_ids"]:
                raise PopulationError("Sequence sample is not mapped to its canonical donor")
        elif kind == "alignment":
            configuration = get(row["design_configuration_id"])
            design = get(configuration["design_id"])
            sequence = get(configuration["sequence_configuration_id"])
            if "".join(col["guide"] or "" for col in row["columns"]) != design["guide_sequence"]:
                raise PopulationError("Alignment uses a different guide from its design")
            for col in row["columns"]:
                if col["target"] is not None:
                    pos = col["target_path_index"]
                    if pos >= len(sequence["sequence"]) or col["target"] != sequence["sequence"][pos]:
                        raise PopulationError("Alignment target does not match its sequence path")
            if any(point not in sequence["path"] for point in row["pam_path"]):
                raise PopulationError("Alignment PAM path is not part of the target sequence")
        elif kind == "donor_configuration":
            occurrence = get(row["sequence_occurrence_id"])
            configuration = get(row["design_configuration_id"])
            if (row["donor_id"], row["sample_id"], row["copy_index"]) != (occurrence["donor_id"], occurrence["sample_id"], occurrence["copy_index"]):
                raise PopulationError("Donor/configuration linkage changes the source copy")
            if configuration["sequence_configuration_id"] != occurrence["sequence_configuration_id"]:
                raise PopulationError("Donor configuration does not use that donor's sequence")
            same = occurrence["sequence_configuration_id"] == row["carrier_of_sequence_id"]
            if row["carrier_status"].endswith("noncarrier") and same or row["carrier_status"].endswith("_carrier") and not same:
                raise PopulationError("Carrier membership disagrees with its defined sequence")
            if row["carrier_status"].startswith("supported_") and get(occurrence["sequence_configuration_id"])["support_status"] != "supported":
                raise PopulationError("Unresolved or modeled DNA cannot establish supported carriage")
        elif kind == "search_scope":
            profile = get(row["search_profile_id"])
            if row["design_id"] != profile["design_id"]:
                raise PopulationError("Search scope uses another design's profile")
        elif kind == "observation_revision":
            revisions.setdefault(row["source_record_id"], set()).add(row["record_id"])
            validate_observation_locus(row, index)
            endpoint = get(row["endpoint_id"])
            source = get(row["source_record_id"])
            if row["normalized_unit"] != endpoint["unit"]:
                raise PopulationError("Normalized observation unit differs from its endpoint")
            normalized = row["normalized_value"]["value"]
            if normalized is not None:
                try:
                    expected = Decimal(str(source["raw_value"]["value"]))
                except InvalidOperation as exc:
                    raise PopulationError("Numeric normalization requires a numeric raw source value") from exc
                conversion = row["conversion"]["value"]
                if conversion is not None:
                    expected = expected * Decimal(str(conversion["factor"])) + Decimal(str(conversion["offset"]))
                elif source["raw_unit"] != row["normalized_unit"]:
                    raise PopulationError("Changed observation units require an explicit conversion")
                if not expected.is_finite() or expected != Decimal(str(normalized)):
                    raise PopulationError("Normalized value differs from the preserved source/conversion")
            if row["applicability"]["empirical_donor_credit"] and endpoint["requires_sequence"] and row["sequence_configuration_id"]["value"] is None:
                raise PopulationError("Sequence-specific empirical credit requires its exact sequence")
        elif kind == "observation_state":
            source = row["source_record_id"]
            if source in states:
                raise PopulationError("Multiple active normalization states for one source measurement")
            states[source] = row
            if row["active_revision_id"] not in row["revision_ids"] or any(get(identifier)["source_record_id"] != source for identifier in row["revision_ids"]):
                raise PopulationError("Normalization state mixes source measurements or omits its active revision")
        elif kind == "observation_equivalence":
            if row["canonical_source_record_id"] not in row["source_record_ids"]:
                raise PopulationError("Equivalence must include its canonical source measurement")
            for identifier in row["source_record_ids"]:
                if identifier in equivalences:
                    raise PopulationError("Overlapping equivalence groups must be reconciled explicitly")
                equivalences[identifier] = row["record_id"]
        elif kind == "evidence_applicability":
            observation = get(row["observation_id"])
            if row["assessment"]["empirical_donor_credit"]:
                if observation["schema_version"] == "2" and outcome_reasons(observation):
                    raise PopulationError("Evidence credit contradicts the adjudicated observation outcome")
                if any(observation[key]["value"] != row[key] for key in ("design_id", "context_id")) or observation["endpoint_id"] != row["endpoint_id"] or observation["subject"]["donor_id"] != row["donor_id"]:
                    raise PopulationError("Empirical evidence cannot change design, context, endpoint or donor")
                if observation["sequence_configuration_id"] != row["sequence_configuration_id"]:
                    raise PopulationError("Empirical evidence cannot change the measured sequence")
                if not observation["applicability"]["empirical_donor_credit"]:
                    raise PopulationError("An applicability record cannot upgrade unusable empirical evidence")
        elif kind == "question":
            for identifier in row["requirement_slot_ids"]:
                if get(identifier)["question_id"] != row["record_id"]:
                    raise PopulationError("Question includes a slot from a different question")
        elif kind == "requirement_slot":
            question = get(row["question_id"])
            if row["record_id"] not in question["requirement_slot_ids"] or row["arm_id"] not in {arm["arm_id"] for arm in question["arms"]}:
                raise PopulationError("Requirement slot is not declared by its question/arm")
            if row["replicate_position"] > question["technical_replicates"]:
                raise PopulationError("Slot replicate exceeds the declared technical-repeat policy")
            for identifier in row["allowed_action_ids"]:
                action = get(identifier)
                if action["donor_id"] != row["donor_id"] or action["context_id"] != row["context_id"] or action["design_id"]["value"] != row["design_id"]["value"] or row["endpoint_id"] not in action["endpoint_ids"] or action["replicate_position"] != row["replicate_position"]:
                    raise PopulationError("Allowed action does not perform the declared donor/slot operation")
            for identifier in row["allowed_observation_ids"]:
                observation = get(identifier)
                if not question["evidence_reuse"]["allow_existing"] or observation["source_record_id"] in question["evidence_reuse"]["excluded_observation_ids"]:
                    raise PopulationError("Observation violates the question's evidence-reuse restrictions")
                if observation["subject"]["donor_id"]["value"] != row["donor_id"]:
                    raise PopulationError("Observation cannot stand in for another independent donor")
                if not observation["applicability"]["empirical_donor_credit"] or observation["design_id"]["value"] != row["design_id"]["value"] or observation["context_id"]["value"] != row["context_id"] or observation["endpoint_id"] != row["endpoint_id"]:
                    raise PopulationError("Observation is not applicable to the slot's empirical requirement")
                source = get(observation["source_record_id"])
                if source["source_experiment"] in question["evidence_reuse"]["excluded_experiment_ids"]:
                    raise PopulationError("A replication cannot reuse the excluded original experiment")
        elif kind == "resource":
            pool_id = row["pool_id"]["value"]
            if pool_id is None:
                if row["instance_index"] != 0 or row["max_pool_instances"] != 1:
                    raise PopulationError("A standalone resource must declare exactly one instance")
            else:
                pool = resource_pools.setdefault(pool_id, {})
                if row["instance_index"] in pool:
                    raise PopulationError("Duplicate resource instance index in a finite pool")
                if pool:
                    previous = next(iter(pool.values()))
                    if row["max_pool_instances"] != previous["max_pool_instances"]:
                        raise PopulationError("Resource pool has inconsistent maximum instance counts")
                    if any(row[field] != previous[field] for field in ("resource_type", "cost", "cost_unit", "cost_scale")):
                        raise PopulationError("Resource pool instances require one declared type and per-instance cost")
                pool[row["instance_index"]] = row
        elif kind == "action":
            if get(row["specimen_id"])["donor_id"] != row["donor_id"]:
                raise PopulationError("Action specimen belongs to a different donor")
            assay = get(row["assay_capability_id"])
            if row["context_id"] not in assay["context_ids"] or set(row["endpoint_ids"]) - set(assay["endpoint_ids"]):
                raise PopulationError("Action exceeds its assay's context/endpoint capability")
            if set(row["locus_ids"]) - set(assay["locus_ids"]) or len(row["locus_ids"]) > assay["max_targets"]:
                raise PopulationError("Action exceeds its assay's target scope/capacity")
            if "endpoint_locus_pairs" in row and "endpoint_locus_pairs" in assay:
                if {canonical_bytes(pair) for pair in row["endpoint_locus_pairs"]} - {canonical_bytes(pair) for pair in assay["endpoint_locus_pairs"]}:
                    raise PopulationError("Action endpoint/locus pair exceeds its exact assay capability")
            for target in row.get("target_sequences", []):
                if target["locus_id"] not in row["locus_ids"]:
                    raise PopulationError("Action target sequence is outside its declared loci")
                for identifier in target["sequence_configuration_ids"]:
                    sequence = get(identifier)
                    if sequence["locus_id"] != target["locus_id"]:
                        raise PopulationError("Action target sequence resolves another locus")
                    if row["eligibility_status"] == "eligible" and not any(candidate["record_type"] == "comparison_sequence_occurrence" and candidate["donor_id"] == row["donor_id"] and candidate["sequence_configuration_id"] == identifier for candidate in index.values()):
                        raise PopulationError("Panel sequence is not established for its actual donor")
            if row["target_sequence_id"]["value"] is not None and row["eligibility_status"] == "eligible":
                if not any(candidate["record_type"] == "comparison_sequence_occurrence" and candidate["donor_id"] == row["donor_id"] and candidate["sequence_configuration_id"] == row["target_sequence_id"]["value"] for candidate in index.values()):
                    raise PopulationError("Action sequence is not established for its actual donor")
            demands, alternatives = {}, {}
            for demand in row["resource_demands"]:
                key = ("fixed", demand["resource_id"]) if "resource_id" in demand else ("alternative", demand["assignment_key"])
                options = [demand["resource_id"]] if "resource_id" in demand else demand["eligible_resource_ids"]
                if len({get(identifier)["resource_type"] for identifier in options}) != 1:
                    raise PopulationError("Alternative resource instances must perform the same resource type")
                if key in alternatives and set(alternatives[key]) != set(options):
                    raise PopulationError("One action assignment key has contradictory eligible instances")
                alternatives[key] = options
                demands[(key, demand["dimension"])] = demands.get((key, demand["dimension"]), 0) + demand["quantity"]
            for (key, dimension), quantity in demands.items():
                for identifier in alternatives[key]:
                    resource = get(identifier)
                    capacities = {cap["dimension"]: cap["capacity"] for cap in resource["capacities"]}
                    if dimension not in capacities:
                        raise PopulationError("Action demands an undeclared resource capacity dimension")
                    capacity = capacities[dimension]
                    if capacity["kind"] == "finite" and quantity > capacity["quantity"]:
                        raise PopulationError("A single action exceeds its resource instance capacity")
                    reasons = resource_compatible(row, resource, index)
                    if reasons:
                        raise PopulationError(reasons[0])
        elif kind == "identity_bridge" and row["new_record_id"] not in index:
            raise PopulationError("Legacy bridge does not resolve its new record")
        elif kind == "experiment_plan":
            selected = set(row["selected_action_ids"])
            assignments = {(item["action_id"], item["assignment_key"]):item["resource_id"] for item in row["resource_assignments"]}
            required_assignments, active = set(), set()
            for identifier in selected:
                action = get(identifier)
                if set(action["dependency_action_ids"]) - selected:
                    raise PopulationError("Plan omits an action dependency")
                for demand in action["resource_demands"]:
                    if "resource_id" in demand:
                        active.add(demand["resource_id"])
                    else:
                        key = (identifier, demand["assignment_key"])
                        required_assignments.add(key)
                        if assignments.get(key) not in demand["eligible_resource_ids"]:
                            raise PopulationError("Plan has no eligible resource instance for an action assignment")
                        active.add(assignments[key])
            if required_assignments != assignments.keys():
                raise PopulationError("Plan resource assignments contain unused or missing choices")
            active.update(identifier for identifier, resource in index.items() if resource["record_type"] == "comparison_resource" and (resource["activation_rule"] == "always" or selected and resource["activation_rule"] == "nonempty_plan"))
            if active != set(row["active_resource_ids"]) or row["total_cost"] != sum(get(identifier)["cost"] for identifier in active):
                raise PopulationError("Plan total differs from its actually activated named resources")
            if any((get(identifier)["cost_unit"],get(identifier)["cost_scale"]) != (row["cost_unit"],row["cost_scale"]) for identifier in active):
                raise PopulationError("Plan mixes incompatible resource cost units or scales")
            complete = [state for state in row["question_states"] if state["planning_status"] in {"already_satisfied", "fully_planned"}]
            if row["unweighted_completion"] != len(complete) or row["weighted_completion"] != sum(get(state["question_id"])["priority_weight"] for state in complete):
                raise PopulationError("Plan completion counts disagree with its declared question states")
            if row["existing_evidence_weight"] != sum(get(state["question_id"])["priority_weight"] for state in complete if state["planning_status"] == "already_satisfied"):
                raise PopulationError("Plan existing evidence weight disagrees with completed observed requirements")
            if row["meets_request"] and any(get(state["question_id"])["hard_required"] and state not in complete for state in row["question_states"]):
                raise PopulationError("Plan claims to meet the request with an incomplete hard question")
        elif kind in {"design_summary", "donor_design_summary"}:
            scopes = [get(identifier) for identifier in row["search_scope_ids"]]
            for scope_id in row["search_scope_ids"]:
                scope = get(scope_id)
                if scope["design_id"] != row["design_id"] or get(scope["search_profile_id"])["domain_id"] != row["domain_id"]:
                    raise PopulationError("Summary search scope belongs to another design/domain")
                if kind == "donor_design_summary" and scope["donor_id"]["value"] != row["donor_id"]:
                    raise PopulationError("Donor summary includes another donor's search scope")
            if kind == "donor_design_summary":
                assessments = [get(identifier) for identifier in row["intended_assessment_ids"]]
                for identifier in row["intended_assessment_ids"]:
                    assessment = get(identifier)
                    if (assessment["donor_id"], assessment["design_id"]) != (row["donor_id"], row["design_id"]):
                        raise PopulationError("Donor summary includes another donor/design assessment")
                if any(get(identifier)["design_id"] != row["design_id"] for identifier in row["nominated_configuration_ids"]):
                    raise PopulationError("Donor summary includes another design's nominated configuration")
                if {assessment["locus_id"] for assessment in assessments} != set(get(row["design_id"])["intended_locus_ids"]):
                    raise PopulationError("Donor summary must assess every intended locus")
                intended = "unresolved" if any(assessment["support_status"] in {"unresolved", "not_evaluated"} for assessment in assessments) else "modeled" if any(assessment["support_status"] == "modeled" for assessment in assessments) else "supported"
                if row["intended_evaluability"] != intended:
                    raise PopulationError("Donor intended evaluability contradicts its assessments")
                if (row["search_status"] == "complete") != all(scope["execution_status"] == "complete" for scope in scopes):
                    raise PopulationError("Donor summary completeness contradicts its searched scope")
                if any(get(identifier)["design_configuration_id"] not in row["nominated_configuration_ids"] for identifier in row["nominated_alignment_ids"]):
                    raise PopulationError("Donor summary alignment does not resolve a listed nominated configuration")
            else:
                donor_ids = {scope["donor_id"]["value"] for scope in scopes if scope["donor_id"]["value"] is not None}
                complete_ids = {identifier for identifier in donor_ids if all(scope["execution_status"] == "complete" for scope in scopes if scope["donor_id"]["value"] == identifier)}
                if row["donor_count"] != len(donor_ids) or row["search_evaluable_donors"] != len(complete_ids):
                    raise PopulationError("Design summary denominator/completeness contradicts its searched donors")
    for source, identifiers in revisions.items():
        if source not in states or identifiers != set(states[source]["revision_ids"]):
            raise PopulationError("Every normalized revision requires an explicit active/history state")
    for row in index.values():
        if row["record_type"] == "comparison_requirement_slot":
            question = get(row["question_id"])
            excluded = question["evidence_reuse"]["excluded_observation_ids"]
            excluded_groups = {equivalences[item] for item in excluded if item in equivalences}
            for identifier in row["allowed_observation_ids"]:
                source = get(identifier)["source_record_id"]
                if states[source]["active_revision_id"] != identifier or equivalences.get(source) in excluded_groups:
                    raise PopulationError("Superseded or equivalent motivating observations cannot fill a slot")
    _acyclic(index, "action", "dependency_action_ids")
    _acyclic(index, "source_artifact", "parent_artifact_id")


def validate_catalog(records):
    """Validate complete cross-record references and return an ID index."""
    if not isinstance(records, list):
        raise PopulationError("A comparison catalog must be a record array")
    index = {}
    for row in records:
        validate_record(row)
        if row["record_id"] in index:
            raise PopulationError("Duplicate canonical ID; store aliases as alias records")
        index[row["record_id"]] = row
    for row in records:
        _foreign_keys(row, record_schema(row["record_type"], row["schema_version"]), index)
    _catalog_semantics(index)
    return index


def validate_study(value):
    """Validate a comparison configuration without acquiring or reading its inputs."""
    schema = study_schema()
    validate_schema(value, schema, "comparison study")
    # Catalogs have their own type-dispatched walk below.
    metadata = {key: item for key, item in value.items() if key != "catalogs"}
    _walk(metadata, schema, _check_value, "comparison study")
    validate_basis(value["mode"], value["coverage_basis"])
    index = validate_catalog(value["catalogs"])
    _foreign_keys(value["policies"], schema["properties"]["policies"], index)
    sources = set()
    for source in value["input_refs"]:
        if source["input_id"] in sources:
            raise PopulationError("Duplicate comparison input ID")
        sources.add(source["input_id"])
        if set(source["stages"]) != set(INPUT_ROLE_STAGES[source["role"]]):
            raise PopulationError("Input stages do not match the required dependencies for " + source["role"])
    policies = value["policies"]
    if policies["search"]["domain_id"] != policies["comparison"]["domain_id"]:
        raise PopulationError("Design comparisons require a shared declared search domain")
    design_ids = set(policies["search"]["design_ids"])
    profiles = [index[identifier] for identifier in policies["search"]["profile_ids"]]
    if {row["design_id"] for row in profiles} != design_ids or len(profiles) != len(design_ids):
        raise PopulationError("Every active comparison design requires exactly one effective search profile")
    if any(profile["domain_id"] != policies["search"]["domain_id"] for profile in profiles):
        raise PopulationError("A search profile falls outside the shared comparison domain")
    if value["coverage_basis"] == "retained_callset_model" and policies["evidence"]["allow_empirical_donor_credit"]:
        raise PopulationError("Retained-callset simulations cannot use empirical donor credits")
    for row in index.values():
        if row["record_type"] == "comparison_resource" and row["cost_unit"] != policies["planning"]["cost_unit"]:
            raise PopulationError("Resource costs must use the declared study cost unit")
    if policies["comparison"]["include_predictions"] and not policies["predictions"]["enabled"]:
        raise PopulationError("Prediction comparison requires an enabled prediction policy")
    for identifier in policies["predictions"]["shortlist_ids"]:
        if identifier not in index:
            raise PopulationError("Prediction shortlist contains an unresolved record")
    if value["mode"] == "research_sample_planning":
        if policies["questions"]["relatedness_policy"] != "strict":
            raise PopulationError("Research specimen planning cannot assume unrelatedness")
        for row in index.values():
            if row["record_type"] == "comparison_specimen" and row["availability"] == "simulation_assumption":
                raise PopulationError("Research specimen planning requires actual specimen availability")
    from .comparison_domain import require_supported_request
    require_supported_request(value)
    return value


def active_observations(records):
    """Return one active normalization per documented source measurement."""
    index = validate_catalog(records)
    aliases = {identifier: row["canonical_source_record_id"]
               for row in records if row["record_type"] == "comparison_observation_equivalence"
               for identifier in row["source_record_ids"]}
    states = {row["source_record_id"]: row for row in records
              if row["record_type"] == "comparison_observation_state"}
    result = []
    for source, state in sorted(states.items()):
        if aliases.get(source, source) == source:
            result.append(index[state["active_revision_id"]])
    return result


def bridge_legacy_record(legacy, new_record, *, family, legacy_record_id, reason,
                         relationship="incompatible_identity_semantics", mapping_version="1"):
    """Preserve a legacy payload and an explicit mapping; do not reinterpret it."""
    versions = {"population": "1", "candidate_ledger": "3", "design": "7"}
    if family not in versions or str(legacy.get("schema_version")) != versions[family]:
        raise PopulationError("Legacy bridge requires its declared original schema version")
    if family == "population":
        from .contracts import validate_record as validate_population_record
        validate_population_record(legacy)
    elif family == "candidate_ledger":
        from ..pipeline.candidate_panel import compare_ledgers
        compare_ledgers(legacy, legacy)
    elif legacy.get("record_type") not in {"hbg_design_record", "hbg_design_record_set"}:
        raise PopulationError("Expected a legacy design record or panel envelope")
    validate_record(new_record)
    payload = canonical_bytes(legacy)
    return make_record("identity_bridge", legacy_family=family,
                       legacy_schema_version=versions[family], legacy_record_id=legacy_record_id,
                       legacy_payload_json=payload.decode("utf-8"),
                       legacy_payload_sha256=hashlib.sha256(payload).hexdigest(),
                       new_record_id=new_record["record_id"], relationship=relationship,
                       mapping_version=mapping_version, reason=reason)


def restore_legacy_record(bridge):
    validate_record(bridge)
    if bridge["record_type"] != "comparison_identity_bridge":
        raise PopulationError("Expected a comparison identity bridge")
    return parse_json(bridge["legacy_payload_json"])
