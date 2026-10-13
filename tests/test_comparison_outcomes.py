"""Outcome adjudication controls evidence, history, planning and stage reuse."""
from copy import deepcopy
import json
from pathlib import Path

import pytest

from base_edit_safety.population.common import (PopulationError, checksum, jsonl_text,
                                                known, read_json, unknown)
from base_edit_safety.population.comparison_contracts import (
    active_observations, make_record, record_identity, record_registry, validate_catalog, validate_record,
)
from base_edit_safety.population.comparison_evidence import paired_descriptions
from base_edit_safety.population.comparison_migration import main as migrate_main, migrate_observation_catalog
from base_edit_safety.population.comparison_outcomes import outcome_reasons
from base_edit_safety.population.comparison_provenance import load_stage, resolve_study, save_stage
from base_edit_safety.population.observations import import_artifact, intrinsic_assessment, merge_history
from test_comparison_evidence import ROWS, POLICY, assess, one, remake
from test_comparison_observations import inputs, write_json
from test_comparison_provenance import build_chain

ROOT = Path(__file__).resolve().parents[1]
LEGACY = read_json(ROOT / "tests/fixtures/comparison/legacy-outcomes-v1.json")
INDEX = {row["record_id"]: row for row in ROWS}


def observation(outcome, *, value=0.2, limit=0.01, unit="fraction", source=None):
    row = one("observation_revision")
    row.update(adjudicated_outcome=outcome, adjudication_reason="Explicit artificial test adjudication",
               normalized_value=known(value) if value is not None else unknown("No usable measurement"),
               source_reported_outcome=known(outcome if source is None else source),
               detection_limit={"value": known(limit) if limit is not None else unknown("Not reported"), "unit": unit})
    row["applicability"] = intrinsic_assessment(row, INDEX)
    return remake(row)


@pytest.mark.parametrize("outcome,value,credit", [
    ("observed_positive", 0.2, True), ("tested_no_signal", 0, True),
    ("inconclusive", 0.2, False), ("untested", 0.2, False),
])
def test_four_outcomes_control_intrinsic_and_reconciled_credit(outcome, value, credit):
    row = observation(outcome, value=value)
    assert row["schema_version"] == row["normalization_schema_version"] == "2"
    assert row["applicability"]["empirical_donor_credit"] is credit
    assert assess(row)["empirical_donor_credit"] is credit
    assert row["source_reported_outcome"] == known(outcome)
    assert row["normalized_value"] == known(value)
    assert validate_record(row) == row


@pytest.mark.parametrize("outcome", ["inconclusive", "untested"])
def test_numeric_nonconclusive_records_cannot_launder_stored_credit(outcome):
    row = observation(outcome)
    row["applicability"]["empirical_donor_credit"] = True
    assert not assess(row)["empirical_donor_credit"]
    row["applicability"]["quality"] = {"status": "applicable", "reason_codes": []}
    row["record_id"] = row["observation_id"] = record_identity(row)
    with pytest.raises(PopulationError, match="adjudicated outcome"):
        validate_record(row)


@pytest.mark.parametrize("value,limit,unit,credit,reason", [
    (0, None, "fraction", False, "no_signal_without_usable_detection_limit"),
    (0.005, None, "fraction", False, "no_signal_without_usable_detection_limit"),
    (0.005, 0, "fraction", False, "no_signal_without_usable_detection_limit"),
    (0.005, 1.1, "fraction", False, "no_signal_without_usable_detection_limit"),
    (0.005, 0.01, "reads", False, "no_signal_without_usable_detection_limit"),
    (0.02, 0.01, "fraction", False, "no_signal_measurement_outside_detection_limit"),
    (0.005, 0.01, "fraction", True, None), (0.005, 1, "percent", True, None),
    (0.01, 1, "percent", True, None),
])
def test_no_signal_requires_a_compatible_limit_even_for_nonzero_values(value, limit, unit, credit, reason):
    row = observation("tested_no_signal", value=value, limit=limit, unit=unit)
    assessment = assess(row, policy={**POLICY, "require_detection_limit_for_zero": False})
    assert assessment["empirical_donor_credit"] is credit
    if reason:
        assert reason in assessment["reason_codes"]
    assert row["normalized_value"]["value"] == value


@pytest.mark.parametrize("reported", ["inconclusive", " Untested ", "tested_no_signal"])
def test_conflicting_source_outcome_is_preserved_and_withholds_credit(reported):
    row = observation("observed_positive", source=reported)
    assert row["source_reported_outcome"]["value"] == reported
    assert "adjudication_contradicts_source_outcome" in assess(row)["reason_codes"]
    assert not assess(row)["empirical_donor_credit"]


@pytest.mark.parametrize("value", [None, 0, -0.1])
def test_positive_adjudication_cannot_replace_a_missing_or_nonpositive_measurement(value):
    row = observation("observed_positive", value=value)
    assert not row["applicability"]["empirical_donor_credit"]
    assert not assess(row)["empirical_donor_credit"]


def test_legacy_stored_credit_remains_archived_but_is_not_runtime_evidence():
    row = deepcopy(next(r for r in LEGACY if r["record_type"] == "comparison_observation_revision"))
    row["applicability"]["empirical_donor_credit"] = True
    assert validate_record(row) == row
    before = deepcopy(row)
    assert intrinsic_assessment(row, INDEX)["reason_codes"] == ["legacy_outcome_requires_adjudication"]
    assert not assess(row)["empirical_donor_credit"]
    assert row == before


def test_outcome_and_reason_are_bound_to_v2_content_identity():
    first = observation("inconclusive")
    second = remake(first, adjudicated_outcome="untested")
    third = remake(first, adjudication_reason="A separate documented adjudication")
    assert len({first["record_id"], second["record_id"], third["record_id"]}) == 3
    tampered = deepcopy(first)
    tampered["adjudicated_outcome"] = "untested"
    with pytest.raises(PopulationError, match="identity disagree"):
        validate_record(tampered)
    with pytest.raises(PopulationError):
        remake(first, adjudicated_outcome="negative")


def test_nonconclusive_values_do_not_create_paired_measured_differences():
    left = observation("inconclusive")
    left.update(pairing_key=known("artificial-pair"))
    editor = remake(one("editor_definition"), variant_name="Second test editor")
    design = remake(one("design"), editor_definition_id=editor["record_id"])
    right = deepcopy(left)
    right.update(record_id="second-test-observation", design_id=known(design["record_id"]), normalized_value=known(0.8))
    assert paired_descriptions([left, right], {**INDEX, design["record_id"]: design, editor["record_id"]: editor},
        {left["design_id"]["value"], design["record_id"]}, POLICY) == []


def test_old_extraction_map_emits_explicit_unresolved_v2_without_source_rewrite(tmp_path):
    path, map_path, catalog, mapping = inputs(tmp_path)
    normalization = mapping["rows"][0]["normalization"]
    normalization.pop("adjudicated_outcome")
    normalization.pop("adjudication_reason")
    write_json(map_path, mapping)
    source_bytes, map_bytes = path.read_bytes(), map_path.read_bytes()
    imported = import_artifact(path, map_path, catalog)
    revision = next(row for row in imported if row["record_type"] == "comparison_observation_revision")
    assert revision["adjudicated_outcome"] == "inconclusive"
    assert "Legacy extraction map" in revision["adjudication_reason"]
    assert not revision["applicability"]["empirical_donor_credit"]
    assert path.read_bytes() == source_bytes and map_path.read_bytes() == map_bytes


def test_migration_preserves_old_sources_revisions_and_invalidates_stale_products():
    before = deepcopy(LEGACY)
    result = migrate_observation_catalog(LEGACY)
    index = validate_catalog(result["records"])
    old = next(row for row in LEGACY if row["record_type"] == "comparison_observation_revision")
    new = index[result["receipt"]["revision_id_map"][old["record_id"]]]
    assert LEGACY == before
    for row in LEGACY:
        if row["record_type"] in {"comparison_source_observation", "comparison_source_artifact", "comparison_observation_revision"}:
            assert index[row["record_id"]] == row
    assert new["adjudicated_outcome"] == "inconclusive"
    assert new["source_reported_outcome"] == old["source_reported_outcome"]
    assert new["normalized_value"] == old["normalized_value"]
    assert active_observations(result["records"]) == [new]
    state = next(row for row in result["records"] if row["record_type"] == "comparison_observation_state")
    assert set(state["revision_ids"]) == {old["record_id"], new["record_id"]}
    old_application = next(row for row in LEGACY if row["record_type"] == "comparison_evidence_applicability")
    assert old_application["record_id"] in result["receipt"]["invalidated_record_ids"]
    assert old_application["record_id"] not in index
    assert result["superseded_states"] == [next(row for row in LEGACY if row["record_type"] == "comparison_observation_state")]
    again = migrate_observation_catalog(result["records"])
    assert again["records"] == result["records"]
    assert again["receipt"]["revision_id_map"] == {}


def test_explicit_migration_and_conflicting_revisions_keep_distinct_ids():
    old = next(row for row in LEGACY if row["record_type"] == "comparison_observation_revision")
    choice = {old["record_id"]: {"adjudicated_outcome": "observed_positive", "adjudication_reason": "Reviewed constructed positive fixture"}}
    resolved = migrate_observation_catalog(LEGACY, choice)
    positive = active_observations(resolved["records"])[0]
    assert positive["applicability"]["empirical_donor_credit"]
    unresolved = migrate_observation_catalog(LEGACY)
    assert positive["record_id"] != active_observations(unresolved["records"])[0]["record_id"]
    merged = merge_history(resolved["records"], unresolved["records"])
    assert active_observations(merged)[0] == positive
    latest_unresolved = merge_history(unresolved["records"], resolved["records"])
    assert active_observations(latest_unresolved)[0]["adjudicated_outcome"] == "inconclusive"
    assert not active_observations(latest_unresolved)[0]["applicability"]["empirical_donor_credit"]
    with pytest.raises(PopulationError, match="active version-1"):
        migrate_observation_catalog(resolved["records"], choice)
    with pytest.raises(PopulationError, match="active version-1"):
        migrate_observation_catalog(LEGACY, {"not-a-revision": choice[old["record_id"]]})
    conflicting = deepcopy(positive)
    conflicting["adjudication_reason"] = "Changed content under the same ID"
    with pytest.raises(PopulationError, match="conflicting content"):
        merge_history([positive], [conflicting])


def test_migration_command_keeps_original_bytes_and_refuses_overwrite(tmp_path):
    source = tmp_path / "old.json"
    source.write_text(json.dumps(LEGACY))
    original = source.read_bytes()
    target = tmp_path / "migrated"
    assert migrate_main(["--catalog", str(source), "--output", str(target)]) == 0
    assert source.read_bytes() == original
    assert validate_catalog(read_json(target / "catalog.json"))
    assert target.stat().st_mode & 0o777 == 0o700
    assert all(p.stat().st_mode & 0o777 == 0o600 for p in target.iterdir())
    with pytest.raises(PopulationError):
        migrate_main(["--catalog", str(source), "--output", str(target)])


def test_migrated_history_rejects_actual_stage_and_parent_reuse(tmp_path):
    source = tmp_path / "reference.fa"
    source.write_bytes((ROOT / "examples/comparison/contracts/reference.fa").read_bytes())
    history = tmp_path / "history.jsonl"
    history.write_text(jsonl_text([row for row in LEGACY if row["record_type"] in {
        "comparison_source_observation", "comparison_observation_revision", "comparison_observation_state", "comparison_observation_equivalence"}]))
    config = read_json(ROOT / "examples/comparison/contracts/study.json")
    config["input_refs"].append({"input_id": "outcome-history", "role": "observation_history", "path": str(history),
        "sha256": checksum(history), "format": "jsonl", "stages": ["observation_normalization"]})
    config_path = tmp_path / "study.json"
    write_json(config_path, config)
    code = {"package_version": "test", "source_sha256": "a" * 64, "python": "test"}
    old = resolve_study(config_path, code=code)
    outputs = build_chain(old, tmp_path / "old-stages")
    old_manifest = (outputs["planning"] / "run-manifest.json").read_bytes()
    migrated = tmp_path / "migrated-history.jsonl"
    migrated.write_text(jsonl_text([row for row in migrate_observation_catalog(LEGACY)["records"] if row["record_type"] in {
        "comparison_source_observation", "comparison_observation_revision", "comparison_observation_state", "comparison_observation_equivalence"}]))
    config["input_refs"][-1].update(path=str(migrated), sha256=checksum(migrated))
    write_json(config_path, config)
    current = resolve_study(config_path, code=code)
    load_stage(current, outputs["search"], "search")
    for stage in ["observation_normalization", "evidence_applicability", "comparison", "questions", "planning", "report"]:
        with pytest.raises(PopulationError, match="stale"):
            load_stage(current, outputs[stage], stage)
    with pytest.raises(PopulationError, match="stale"):
        save_stage(current, "evidence_applicability", tmp_path / "new-evidence", {"check.txt": "new"},
            parents={"search": outputs["search"], "observation_normalization": outputs["observation_normalization"]})
    assert (outputs["planning"] / "run-manifest.json").read_bytes() == old_manifest


def test_no_signal_values_below_the_limit_are_not_exact_paired_differences():
    left = observation("tested_no_signal", value=0.001, limit=0.01)
    left.update(pairing_key=known("artificial-pair"))
    editor = remake(one("editor_definition"), variant_name="Second test editor")
    design = remake(one("design"), editor_definition_id=editor["record_id"])
    right = deepcopy(left)
    right.update(record_id="second-test-observation", design_id=known(design["record_id"]), normalized_value=known(0.005))
    assert not outcome_reasons(left) and not outcome_reasons(right)
    assert paired_descriptions([left, right], {**INDEX, design["record_id"]: design, editor["record_id"]: editor},
        {left["design_id"]["value"], design["record_id"]}, POLICY) == []


def test_archived_positive_application_cannot_grant_current_summary_credit(tmp_path):
    from base_edit_safety.population.comparison_evidence import compare_with_evidence
    from test_comparison_contracts import rekey
    rows = deepcopy(LEGACY)
    next(row for row in rows if row["record_type"] == "comparison_observation_revision")["applicability"]["empirical_donor_credit"] = True
    next(row for row in rows if row["record_type"] == "comparison_evidence_applicability")["assessment"]["empirical_donor_credit"] = True
    config = read_json(ROOT / "examples/comparison/contracts/study.json")
    config["catalogs"] = rekey(rows)
    for ref in config["input_refs"]:
        ref["path"] = str((ROOT / "examples/comparison/contracts" / ref["path"]).resolve())
    config_path = tmp_path / "study.json"
    write_json(config_path, config)
    result = compare_with_evidence(resolve_study(config_path), {}, {}, {})
    summaries = [row for row in result["comparison-summaries.jsonl"] if row["record_type"] == "comparison_donor_design_summary"]
    assert summaries
    assert all(row["evidence_availability"] != "exact_donor_evidence" for row in summaries)


@pytest.mark.parametrize("case", ["observed_positive", "legacy_stored_credit", "latest_inconclusive", "untested"])
def test_direct_planning_does_not_borrow_stale_or_nonconclusive_credit(tmp_path, case):
    from test_comparison_questions import build_fixture, mutated, observation_records, policy_changed
    from base_edit_safety.population.comparison_questions import evaluate_requirements
    fixture = build_fixture(tmp_path / "fixture")
    qset, output = policy_changed(fixture, lambda policy: policy["rules"][0].update(allow_existing=True))
    bundle = {"study": qset.study, "records": output["question-records.jsonl"], "labels": fixture["labels"]}
    observations = observation_records(bundle)
    positive = next(row for row in observations if row["record_type"] == "comparison_observation_revision")
    state = next(row for row in observations if row["record_type"] == "comparison_observation_state")
    current = positive
    if case == "legacy_stored_credit":
        fields = {key: deepcopy(value) for key, value in positive.items() if key not in {
            "record_type", "schema_version", "record_id", "observation_id", "adjudicated_outcome", "adjudication_reason"}}
        fields.update(normalization_schema_version="1", source_reported_outcome=known("untested"))
        current = make_record("observation_revision", version="1", **fields)
        observations.remove(positive)
        observations.append(current)
    elif case != "observed_positive":
        fields = deepcopy(positive)
        fields.update(adjudicated_outcome="inconclusive" if case == "latest_inconclusive" else "untested")
        fields["applicability"]["empirical_donor_credit"] = False
        current = remake(fields)
        if case == "untested":
            observations.remove(positive)
        observations.append(current)
    observations.remove(state)
    observations.append(make_record("observation_state", source_record_id=current["source_record_id"],
        active_revision_id=current["record_id"],
        revision_ids=[row["record_id"] for row in observations if row["record_type"] == "comparison_observation_revision"],
        supersession_reason="Latest artificial outcome review"))
    def change(rows):
        rows.extend(observations)
        for row in rows:
            if row["record_type"] == "comparison_requirement_slot" and row["donor_id"] == fixture["labels"]["donors"]["C"] \
                    and row["design_id"]["value"] == fixture["labels"]["designs"]["A"] and row["control_role"] == "measurement":
                row["allowed_observation_ids"] = [current["record_id"]] if case in {"observed_positive", "legacy_stored_credit"} else []
    prepared, _ = mutated(bundle, change)
    assignments = evaluate_requirements(prepared)["slot_assignments"]
    count = sum("observation_id" in value for value in assignments.values())
    assert count == (1 if case == "observed_positive" else 0)
    if case != "observed_positive":
        expected = "legacy_outcome_requires_adjudication" if case == "legacy_stored_credit" else "outcome_" + current["adjudicated_outcome"]
        from base_edit_safety.population.comparison_questions import _observation_slot_reasons
        slot = next(row for row in prepared.slots.values() if row["donor_id"] == fixture["labels"]["donors"]["C"]
            and row["design_id"]["value"] == fixture["labels"]["designs"]["A"] and row["control_role"] == "measurement")
        assert expected in _observation_slot_reasons(prepared, current, slot)


@pytest.mark.parametrize("outcome", ["inconclusive", "untested"])
def test_nonconclusive_number_does_not_suppress_missing_evidence_question(tmp_path, outcome):
    from test_comparison_questions import build_fixture, observation_records
    from base_edit_safety.population.comparison_question_generation import _events
    from test_comparison_contracts import rekey
    fixture = build_fixture(tmp_path / "fixture")
    observations = observation_records(fixture, credit=False)
    row = next(row for row in observations if row["record_type"] == "comparison_observation_revision")
    row["adjudicated_outcome"] = outcome
    records = rekey([*fixture["study"]["config"]["catalogs"], *fixture["records"], *observations])
    index = validate_catalog(records)
    rule = {"design_ids": list(fixture["labels"]["designs"].values()),
            "candidate_endpoint_id": fixture["labels"]["endpoints"]["L"], "trigger": "missing_functional_evidence"}
    assert _events(index, rule, synthetic=True)


@pytest.mark.parametrize("missing", ["adjudicated_outcome", "adjudication_reason"])
def test_partially_adjudicated_extraction_map_is_an_error(tmp_path, missing):
    path, map_path, catalog, mapping = inputs(tmp_path)
    mapping["rows"][0]["normalization"].pop(missing)
    write_json(map_path, mapping)
    with pytest.raises(PopulationError, match="observation normalization"):
        import_artifact(path, map_path, catalog)


def test_censored_no_signal_without_numeric_measurement_is_preserved_but_not_credited():
    row = observation("tested_no_signal", value=None, limit=0.01, source="below stated limit")
    row["measurement_kind"] = "censored"
    assert not intrinsic_assessment(row, INDEX)["empirical_donor_credit"]
    result = assess(row)
    assert not result["empirical_donor_credit"]
    assert "measurement_unavailable_or_censored" in result["reason_codes"]
    assert row["adjudicated_outcome"] == "tested_no_signal"
    assert row["normalized_value"]["value"] is None
    assert row["source_reported_outcome"]["value"] == "below stated limit"


def test_positive_value_at_or_below_stated_limit_is_flagged_without_changing_it():
    row = observation("observed_positive", value=0.005, limit=0.01)
    assert not assess(row)["empirical_donor_credit"]
    assert "positive_measurement_not_above_detection_limit" in assess(row)["reason_codes"]
    assert row["normalized_value"] == known(0.005)


def test_out_of_range_fraction_cannot_reenter_through_stored_applicability():
    row = observation("observed_positive", value=2)
    assert "measured_fraction_outside_bounds" in outcome_reasons(row)
    assert not row["applicability"]["empirical_donor_credit"]
    row["applicability"] = deepcopy(observation("observed_positive")["applicability"])
    row["record_id"] = row["observation_id"] = record_identity(row)
    with pytest.raises(PopulationError, match="measured_fraction_outside_bounds"):
        validate_record(row)
    assert not assess(row)["empirical_donor_credit"]
