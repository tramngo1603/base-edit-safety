"""Canonical identities, strict record families and scientific reference boundaries."""
from copy import deepcopy
from pathlib import Path
import json

import pytest

from base_edit_safety.population.common import PopulationError, parse_json, read_json
from base_edit_safety.population.comparison_contracts import (
    active_observations, bridge_legacy_record, make_record, record_identity,
    record_registry, record_schema, restore_legacy_record, validate_catalog,
    validate_record, validate_study,
)
from base_edit_safety.population.comparison_identity import canonical_bytes, content_id

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests/fixtures/comparison/canonical-records.json"
RECORDS = read_json(FIXTURE)
KINDS = tuple(record_registry())


def record(kind, rows=None):
    return next(row for row in (RECORDS if rows is None else rows)
                if row["record_type"] == "comparison_" + kind)


def remake(row, **changes):
    values = deepcopy(row)
    for key in ("record_type", "schema_version", "record_id", record_registry()[row["record_type"]]):
        values.pop(key, None)
    values.pop("definition_hash", None)
    values.update(changes)
    return make_record(row["record_type"], **values)


def rekey(rows):
    """Rebind test mutations through the identity DAG without production resolvers."""
    rows = deepcopy(rows)
    def replace(value, mapping):
        if isinstance(value, str):
            return mapping.get(value, value)
        if isinstance(value, list):
            return [replace(item, mapping) for item in value]
        if isinstance(value, dict):
            return {key: replace(item, mapping) for key, item in value.items()}
        return value
    for _ in range(32):
        changes = {row["record_id"]: record_identity(row) for row in rows
                   if row["record_id"] != record_identity(row)}
        if not changes:
            return rows
        rows = replace(rows, changes)
        for row in rows:
            if row["record_type"] == "comparison_design":
                row["definition_hash"] = row["record_id"].rsplit(":", 1)[1]
    raise AssertionError("Test identity graph did not converge")


@pytest.mark.parametrize("kind", KINDS)
def test_each_record_family_has_a_roundtrip_and_rejects_unknown_keys(kind):
    row = next(row for row in RECORDS if row["record_type"] == kind)
    assert validate_record(parse_json(json.dumps(row))) == row
    bad = deepcopy(row)
    bad["unregistered_scientific_value"] = 0
    with pytest.raises(PopulationError, match="unexpected fields"):
        validate_record(bad)


@pytest.mark.parametrize("kind", KINDS)
def test_each_record_family_rejects_missing_and_wrong_type_fields(kind):
    row = deepcopy(next(row for row in RECORDS if row["record_type"] == kind))
    schema = record_schema(kind)
    field = next(key for key in schema["required"]
                 if key not in {"record_type", "schema_version", "record_id", schema["x-id-field"]})
    missing = deepcopy(row)
    missing.pop(field)
    with pytest.raises(PopulationError):
        validate_record(missing)
    row[field] = {"invented": True}
    with pytest.raises(PopulationError):
        validate_record(row)


@pytest.mark.parametrize("kind", KINDS)
def test_each_record_family_rejects_foreign_namespaces_and_version_changes(kind):
    row = deepcopy(next(row for row in RECORDS if row["record_type"] == kind))
    row["record_id"] = "cfg:" + "a" * 64
    with pytest.raises(PopulationError):
        validate_record(row)
    row["schema_version"] = "99"
    with pytest.raises(PopulationError, match="type/version"):
        validate_record(row)


def test_complete_positive_catalog_and_study():
    assert set(row["record_type"] for row in RECORDS) == set(KINDS)
    assert len(validate_catalog(RECORDS)) == len(RECORDS)
    config = read_json(ROOT / "examples/comparison/contracts/study.json")
    assert validate_study(config) == config
    assert validate_catalog(list(reversed(RECORDS))) == validate_catalog(RECORDS)


def test_unreported_editor_window_stays_unknown():
    window = {"value": None, "reason": "The source does not establish an editor-wide activity window"}
    editor = remake(record("editor_definition"), editing_window=window)
    assert validate_record(editor)["editing_window"] == window
    assert editor["editor_definition_id"] != record("editor_definition")["editor_definition_id"]
    assert record_identity(record("editor_definition")) == record("editor_definition")["record_id"]


@pytest.mark.parametrize("window", [
    {"value": None},
    {"value": None, "reason": ""},
    {"value": None, "reason": "Unreported", "start_1based": 1, "end_1based": 20},
    {"start_1based": 0, "end_1based": 20, "basis": "Reported interval"},
])
def test_unknown_editor_window_requires_a_reason_and_cannot_mix_bounds(window):
    with pytest.raises(PopulationError):
        remake(record("editor_definition"), editing_window=window)


@pytest.mark.parametrize("value", ['{"a":1,"a":2}', '{"x":NaN}', '{"x":Infinity}', '{"x":1e999}'])
def test_strict_file_boundary_rejects_duplicate_keys_and_nonfinite_values(value):
    with pytest.raises(PopulationError):
        parse_json(value)


@pytest.mark.parametrize("value", [True, float("nan"), float("inf"), -1])
def test_costs_reject_booleans_nonfinite_values_and_negative_numbers(value):
    row = deepcopy(record("resource"))
    row["cost"] = value
    with pytest.raises(PopulationError):
        validate_record(row)


def test_nested_objects_are_closed_and_unknown_science_has_a_reason():
    row = deepcopy(record("experimental_context"))
    row["dose"]["value"] = {"value": None, "reason": None}
    with pytest.raises(PopulationError):
        validate_record(row)
    row = deepcopy(record("experimental_context"))
    row["dose"]["units"] = "fraction"
    with pytest.raises(PopulationError):
        validate_record(row)
    row["dose"].pop("units")
    row["cell_type"] = "REPLACE_CELL_TYPE"
    with pytest.raises(PopulationError, match="placeholder"):
        validate_record(row)


def test_canonical_numbers_unicode_and_only_declared_set_order():
    assert canonical_bytes({"x": 1.0, "y": -0.0}) == b'{"x":1,"y":0}'
    assert canonical_bytes({"name": "e\u0301"}) == canonical_bytes({"name": "\u00e9"})
    with pytest.raises(PopulationError, match="collide"):
        canonical_bytes({"e\u0301": 1, "\u00e9": 2})
    assert canonical_bytes({"set": [2, 1]}, set_fields=["set"]) == canonical_bytes({"set": [1, 2]}, set_fields=["set"])
    assert canonical_bytes({"ALT": ["A", "G"]}) != canonical_bytes({"ALT": ["G", "A"]})
    assert canonical_bytes({"replicates": [1, 2]}) != canonical_bytes({"replicates": [2, 1]})
    assert content_id("comparison_design", {}) != content_id("comparison_alignment", {})


def test_canonical_nested_text_and_scalar_spellings():
    value = {"text": "e\u0301\n\"\\", "values": [True, False, None, 9007199254740993]}
    assert canonical_bytes(value) == '{"text":"é\\n\\\"\\\\","values":[true,false,null,9007199254740993]}'.encode("utf-8")
    assert canonical_bytes({"values": [1e-7, -0.0, 1.0, 1e20]}) == b'{"values":[0.0000001,0,1,100000000000000000000]}'
    with pytest.raises(PopulationError, match="Duplicate"):
        canonical_bytes({"values": [1, 1.0]}, set_fields=["values"])
    with pytest.raises(PopulationError, match="finite JSON"):
        canonical_bytes({"nested": [{"value": float("nan")} ]})


def test_source_alias_deduplication_keeps_complete_design_biology():
    design = record("design")
    renamed = remake(design, display_name="another source alias", source_refs=[])
    assert renamed["design_id"] == design["design_id"]
    editor = remake(record("editor_definition"), variant_name="a different editor")
    changed = remake(design, editor_definition_id=editor["record_id"])
    assert changed["design_id"] != design["design_id"]
    edits = deepcopy(design["desired_edit_definitions"])
    edits[0]["alternate_base"] = "T"
    assert remake(design, desired_edit_definitions=edits)["design_id"] != design["design_id"]


def test_identical_hbg_copy_dna_retains_distinct_locus_path_identity():
    old = record("intended_locus")
    path = [{**point, "pos0": point["pos0"] + 1000} for point in old["path"]]
    other = remake(old, gene="HBG2", locus_name="synthetic-HBG2", start=1100, end=1123, path=path)
    sequence = record("sequence_configuration")
    new_sequence = remake(sequence, locus_id=other["record_id"], path=path)
    assert new_sequence["sequence"] == sequence["sequence"]
    assert new_sequence["record_id"] != sequence["record_id"]


def test_sequence_content_excludes_donor_labels_and_search_scores():
    sequence = record("sequence_configuration")
    occurrence = record("sequence_occurrence")
    other = remake(occurrence, sample_id="another-sample", source_version="another-source")
    assert other["sequence_configuration_id"] == sequence["record_id"]
    assert other["record_id"] != occurrence["record_id"]
    changed_support = remake(sequence, support_status="modeled", support_reason="Retained-callset simulation")
    assert changed_support["record_id"] == sequence["record_id"]
    with pytest.raises(PopulationError):
        remake(sequence, donor_id="another-donor")
    with pytest.raises(PopulationError):
        remake(sequence, mismatch_score=4)


def test_alignment_identity_preserves_columns_but_excludes_native_row_aliases():
    alignment = record("alignment")
    assert remake(alignment, native_row_refs=[])["record_id"] == alignment["record_id"]
    bad = deepcopy(alignment)
    bad["mismatch_positions"] = [1]
    with pytest.raises(PopulationError, match="annotations"):
        validate_record(bad)
    bad = deepcopy(alignment)
    bad["columns"][0]["guide"] = None
    bad["columns"][0]["target"] = None
    bad["record_id"] = bad["alignment_id"] = record_identity(bad)
    with pytest.raises(PopulationError, match="two gaps"):
        validate_record(bad)


def test_missing_or_filtered_genotypes_cannot_claim_supported_carriage():
    for status in ("missing", "filtered", "retained_callset_model"):
        with pytest.raises(PopulationError, match="supported carriage"):
            remake(record("donor_configuration"), genotype_support=status)


def test_partial_search_never_acquires_a_complete_claim():
    scope = record("search_scope")
    with pytest.raises(PopulationError, match="cannot be complete"):
        remake(scope, exit_code={"value": 1, "reason": None})
    with pytest.raises(PopulationError, match="cannot be complete"):
        remake(scope, retention_status="truncated")
    with pytest.raises(PopulationError, match="cannot be complete"):
        remake(scope, skipped_intervals=[{"interval": scope["requested_intervals"][0], "reason": "density cap"}])


def test_regional_checksum_cannot_be_passed_as_whole_reference():
    row = record("reference_bundle")
    with pytest.raises(PopulationError, match="indexed-window"):
        remake(row, identity_kind="indexed_windows", windows=[{"interval": {"contig": "chr11", "start": 100, "end": 123}, "source_sha256": row["source_checksums"][0]}])


def test_distinct_source_artifacts_do_not_change_reference_content_identity():
    source = record("source_artifact")
    other = remake(source, source_id="another-study-same-bytes", uri="urn:bes:synthetic:another-origin")
    assert other["sha256"] == source["sha256"]
    assert other["record_id"] != source["record_id"]
    bundle = record("reference_bundle")
    alias_bundle = remake(bundle, source_artifact_ids=[other["record_id"]])
    assert alias_bundle["record_id"] == bundle["record_id"]


def test_regional_window_identity_binds_each_window_to_its_own_bytes():
    bundle = record("reference_bundle")
    windows = [{"interval": {"contig": "chr11", "start": 100, "end": 123}, "source_sha256": "a" * 64},
               {"interval": {"contig": "chr11", "start": 200, "end": 223}, "source_sha256": "b" * 64}]
    first = remake(bundle, identity_kind="indexed_windows", complete_reference_sha256={"value": None, "reason": "regional indexed retrieval"}, windows=windows, source_checksums=["a" * 64, "b" * 64])
    windows = deepcopy(windows)
    windows[0]["source_sha256"], windows[1]["source_sha256"] = windows[1]["source_sha256"], windows[0]["source_sha256"]
    second = remake(first, windows=windows)
    assert first["record_id"] != second["record_id"]


def test_foreign_keys_reject_missing_records_and_predictions_as_observations():
    rows = deepcopy(RECORDS)
    alias = record("design_alias", rows)
    alias["design_id"] = "comparison_design:v1:" + "f" * 64
    with pytest.raises(PopulationError, match="foreign key"):
        validate_catalog(rekey(rows))
    rows = deepcopy(RECORDS)
    record("requirement_slot", rows)["allowed_observation_ids"] = [record("prediction", rows)["record_id"]]
    with pytest.raises(PopulationError, match="wrong-type foreign key"):
        validate_catalog(rows)
    rows = deepcopy(RECORDS)
    record("evidence_applicability", rows)["observation_id"] = record("prediction", rows)["record_id"]
    with pytest.raises(PopulationError, match="wrong-type foreign key"):
        validate_catalog(rekey(rows))


def test_parser_revision_supersedes_one_measurement_and_preserves_history():
    rows = deepcopy(RECORDS)
    old = record("observation_revision", rows)
    new = remake(old, extraction_version="fixture-extractor-v2")
    rows.append(new)
    state = record("observation_state", rows)
    state["revision_ids"].append(new["record_id"])
    state["active_revision_id"] = new["record_id"]
    state["supersession_reason"] = "Parser update with unchanged source measurement"
    rows = rekey(rows)
    assert new["source_record_id"] == old["source_record_id"]
    assert new["observation_id"] != old["observation_id"]
    assert [item["observation_id"] for item in active_observations(rows)] == [new["observation_id"]]
    broken = [item for item in rows if item["record_type"] != "comparison_observation_state"]
    with pytest.raises(PopulationError, match="active/history"):
        validate_catalog(broken)


def test_equal_values_do_not_establish_measurement_equivalence():
    rows = deepcopy(RECORDS)
    source = remake(record("source_observation"), source_replicate="new-independent-replicate")
    revision = remake(record("observation_revision"), source_record_id=source["record_id"])
    state = make_record("observation_state", source_record_id=source["record_id"], active_revision_id=revision["record_id"], revision_ids=[revision["record_id"]], supersession_reason="Initial normalization of a different actual replicate")
    rows.extend([source, revision, state])
    assert len(active_observations(rows)) == 2


def test_observation_normalization_preserves_the_source_value():
    rows = deepcopy(RECORDS)
    record("observation_revision", rows)["normalized_value"]["value"] = 0
    with pytest.raises(PopulationError, match="preserved source"):
        validate_catalog(rekey(rows))


def test_decimal_source_normalization_remains_exact_without_float_rounding():
    rows = deepcopy(RECORDS)
    source = record("source_observation", rows)
    source["raw_value"]["value"] = "9.8000000000000007"
    source["raw_unit"] = "percent"
    normalized = record("observation_revision", rows)
    normalized["normalized_value"]["value"] = "0.098000000000000007"
    normalized["conversion"] = {"value": {"method": "percent_to_fraction", "factor": 0.01, "offset": 0}, "reason": None}
    assert validate_catalog(rekey(rows))
    normalized["normalized_value"]["value"] = 0.098
    with pytest.raises(PopulationError, match="preserved source"):
        validate_catalog(rekey(rows))


@pytest.mark.parametrize("text", ["NaN", "Infinity", "-Infinity", "1.2.3", "", "01.0"])
def test_normalized_decimal_text_rejects_nonfinite_or_invalid_values(text):
    row = deepcopy(record("observation_revision"))
    row["normalized_value"]["value"] = text
    with pytest.raises(PopulationError):
        validate_record(rekey([row])[0])


def test_decimal_text_zero_still_requires_a_detection_limit_for_credit():
    row = deepcopy(record("observation_revision"))
    row["normalized_value"]["value"] = "0.000"
    row["applicability"]["empirical_donor_credit"] = True
    row["detection_limit"]["value"] = {"value": None, "reason": "not reported"}
    with pytest.raises(PopulationError, match="detection limit"):
        validate_record(rekey([row])[0])


def test_zero_without_detection_limit_and_pooled_donor_claims_are_rejected():
    row = deepcopy(record("observation_revision"))
    row["normalized_value"]["value"] = 0
    row["applicability"]["empirical_donor_credit"] = True
    row["detection_limit"]["value"] = {"value": None, "reason": "not reported"}
    with pytest.raises(PopulationError, match="detection limit"):
        validate_record(rekey([row])[0])
    row = deepcopy(record("observation_revision"))
    row["subject"]["level"] = "pooled_group"
    with pytest.raises(PopulationError, match="individual donor"):
        validate_record(rekey([row])[0])


def test_retained_callset_mode_cannot_claim_empirical_donor_credit():
    row = deepcopy(record("evidence_applicability"))
    row["coverage_basis"] = "retained_callset_model"
    row["assessment"]["empirical_donor_credit"] = True
    row["record_id"] = row["evidence_applicability_id"] = record_identity(row)
    with pytest.raises(PopulationError, match="zero empirical"):
        validate_record(row)


def test_observation_matching_cannot_change_editor_or_cell_context():
    rows = deepcopy(RECORDS)
    observation = record("observation_revision", rows)
    observation["applicability"]["empirical_donor_credit"] = True
    changed_context = remake(record("experimental_context", rows), cell_type="synthetic-HEK293T-context")
    rows.append(changed_context)
    assessment = record("evidence_applicability", rows)
    assessment["assessment"]["empirical_donor_credit"] = True
    assessment["context_id"] = changed_context["record_id"]
    with pytest.raises(PopulationError, match="cannot change design, context"):
        validate_catalog(rekey(rows))


def test_pairing_replication_and_explicit_finite_capacities():
    with pytest.raises(PopulationError, match="paired within"):
        remake(record("question"), template="paired_design_comparison")
    reuse = deepcopy(record("question")["evidence_reuse"])
    reuse["require_new"] = []
    with pytest.raises(PopulationError, match="Replication requires"):
        remake(record("question"), template="evidence_replication", evidence_reuse=reuse)
    with pytest.raises(PopulationError, match="finite pool"):
        remake(record("resource"), instance_index=1)
    rows = deepcopy(RECORDS)
    record("action", rows)["resource_demands"][0]["quantity"] = 2
    with pytest.raises(PopulationError, match="instance capacity"):
        validate_catalog(rows)


def test_shared_price_does_not_authorize_incompatible_resource_sharing():
    rows = deepcopy(RECORDS)
    record("action", rows)["batch_id"] = "another-batch"
    with pytest.raises(PopulationError, match="sharing compatibility"):
        validate_catalog(rekey(rows))


def test_resource_timepoint_requires_an_exact_known_context_timepoint():
    rows = deepcopy(RECORDS)
    record("resource", rows)["compatibility"]["timepoints"] = ["24 hour"]
    with pytest.raises(PopulationError, match="timepoint is incompatible"):
        validate_catalog(rekey(rows))
    rows = deepcopy(RECORDS)
    record("experimental_context", rows)["measurement_timepoint"]["value"] = {"value": None, "reason": "not reported"}
    with pytest.raises(PopulationError, match="timepoint is incompatible"):
        validate_catalog(rekey(rows))
    assert validate_catalog(RECORDS)


def test_one_action_sums_repeated_resource_dimension_demands():
    rows = deepcopy(RECORDS)
    resource = record("resource", rows)
    resource["capacities"][0]["capacity"]["quantity"] = 2
    action = record("action", rows)
    action["resource_demands"].append({**action["resource_demands"][0], "quantity": 2})
    with pytest.raises(PopulationError, match="instance capacity"):
        validate_catalog(rekey(rows))
    resource["capacities"][0]["capacity"]["quantity"] = 3
    assert validate_catalog(rekey(rows))


def test_finite_resource_pool_positions_and_declarations_are_consistent():
    rows = deepcopy(RECORDS)
    resource = record("resource", rows)
    duplicate = remake(resource, sharing_key="another-library-same-pool-position")
    with pytest.raises(PopulationError, match="Duplicate resource instance"):
        validate_catalog(rows + [duplicate])
    resource["max_pool_instances"] = 2
    second = remake(resource, sharing_key="second-library", instance_index=1)
    assert validate_catalog(rekey(rows + [second]))
    second["max_pool_instances"] = 3
    with pytest.raises(PopulationError, match="inconsistent maximum"):
        validate_catalog(rekey(rows + [second]))


def test_action_dependencies_are_acyclic():
    rows = deepcopy(RECORDS)
    action = record("action", rows)
    action["dependency_action_ids"] = [action["record_id"]]
    with pytest.raises(PopulationError, match="Cyclic action"):
        validate_catalog(rows)


def test_config_cannot_omit_actual_input_dependencies_or_invent_unknown_policies():
    config = read_json(ROOT / "examples/comparison/contracts/study.json")
    config["input_refs"][0]["stages"] = ["report"]
    with pytest.raises(PopulationError, match="required dependencies"):
        validate_study(config)


def test_source_only_design_retains_identity_without_entering_the_active_search():
    config = read_json(ROOT / "examples/comparison/contracts/study.json")
    editor = remake(record("editor_definition", config["catalogs"]), variant_name="synthetic-published-editor")
    source_design = remake(record("design", config["catalogs"]), editor_definition_id=editor["record_id"], display_name="source-only evidence design")
    config["catalogs"].extend([editor, source_design])
    validate_study(config)
    assert source_design["record_id"] not in config["policies"]["search"]["design_ids"]
    config["policies"]["search"]["design_ids"].append(source_design["record_id"])
    with pytest.raises(PopulationError, match="active comparison design"):
        validate_study(config)
    config = read_json(ROOT / "examples/comparison/contracts/study.json")
    config["policies"]["planning"]["unpriced_controls"] = True
    with pytest.raises(PopulationError, match="unexpected fields"):
        validate_study(config)


def test_explicit_v3_bridge_roundtrips_original_bytes_semantics_and_source_links():
    bridge = record("identity_bridge")
    original = read_json(ROOT / "tests/fixtures/comparison/legacy-roundtrip.json")
    restored = restore_legacy_record(bridge)
    assert restored == original
    from base_edit_safety.population.contracts import validate_record as validate_v3
    assert validate_v3(restored) == original
    target = validate_catalog(RECORDS)[bridge["new_record_id"]]
    assert target["sequence"] == original["sequence_genomic_forward"]
    assert target["path"] == original["sequence_path"]
    assert target["record_id"] != original["sequence_configuration_id"]
    assert original["reference_id"] in {row["sha256"] for row in RECORDS if row["record_type"] == "comparison_source_artifact"}
    bad = deepcopy(bridge)
    bad["legacy_payload_json"] = '{}'
    bad["record_id"] = bad["bridge_id"] = record_identity(bad)
    with pytest.raises(PopulationError, match="checksum"):
        restore_legacy_record(bad)


def test_design_schema7_and_ledger_schema3_are_preserved_as_explicit_archives():
    # Their source readers keep their existing semantics. No measured zero is added.
    design = {"record_type": "hbg_design_record_set", "schema_version": "7", "count": 0, "records": []}
    ledger = {"record_type": "candidate_panel_ledger", "schema_version": "3", "guides": [], "candidates": []}
    for family, original in (("design", design), ("candidate_ledger", ledger)):
        bridge = bridge_legacy_record(original, record("design"), family=family, legacy_record_id="synthetic-archive", reason="Archive-only fixture preserves the original schema")
        assert restore_legacy_record(bridge) == original


def test_donor_summary_cannot_turn_unassessed_evidence_into_empirical_credit():
    with pytest.raises(PopulationError, match="evidence availability"):
        remake(record("donor_design_summary"), evidence_availability="exact_donor_evidence")


def test_summary_support_must_agree_with_actual_referenced_assessments():
    rows = deepcopy(RECORDS)
    record("donor_design_summary", rows)["intended_evaluability"] = "unresolved"
    with pytest.raises(PopulationError, match="contradicts its assessments"):
        validate_catalog(rekey(rows))


def test_design_summary_denominator_must_agree_with_actual_donors():
    rows = deepcopy(RECORDS)
    summary = record("design_summary", rows)
    summary["donor_count"] += 1
    summary["intended_unresolved_donors"] += 1
    summary["search_unresolved_donors"] += 1
    for key in ("complete_locus_count", "complete_sequence_configuration_count", "complete_alignment_count"):
        summary[key] = {"value": None, "reason": "Incomplete hypothetical donor"}
    with pytest.raises(PopulationError, match="searched donors"):
        validate_catalog(rekey(rows))
