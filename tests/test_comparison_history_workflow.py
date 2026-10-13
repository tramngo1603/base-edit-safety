"""Import mock assay results and preserve plans across later QC review."""
from contextlib import redirect_stderr, redirect_stdout
from copy import deepcopy
import io
import json
from pathlib import Path
import shutil

import pytest

from base_edit_safety.pipeline.comparison import main
from base_edit_safety.population.common import checksum, json_text, known, read_json, read_jsonl, unknown
from base_edit_safety.population.comparison_contracts import (
    active_observations, make_record, record_identity, validate_catalog,
)
from base_edit_safety.population.comparison_provenance import load_stage_records, resolve_study
from base_edit_safety.population.provenance import code_identity


ROOT = Path(__file__).resolve().parents[1]


def _write(path, value):
    path.write_text(json_text(value))
    return path


def _invoke(arguments, log, *, expected=0):
    stdout, stderr = io.StringIO(), io.StringIO()
    with redirect_stdout(stdout), redirect_stderr(stderr):
        status = main([str(item) for item in arguments])
    log.append({"arguments": [str(item) for item in arguments], "exit_code": status,
                "stdout": stdout.getvalue(), "stderr": stderr.getvalue()})
    assert status == expected, stderr.getvalue()
    assert "/Users/" not in stdout.getvalue() + stderr.getvalue()
    assert "synthetic-donor-" not in stdout.getvalue() + stderr.getvalue()
    return json.loads(stdout.getvalue() if status == 0 else stderr.getvalue())


def _pipeline(config, directory, log, *, previous=None):
    directory.mkdir(parents=True, exist_ok=False)
    def invoke(arguments):
        try:
            return _invoke(arguments, log)
        finally:
            _write(directory.parent / "command-log.json", log)
    invoke(["validate", "--config", config])
    for command, parents in (
        ("search", ()), ("observations", ()), ("evidence", ("search", "observations")),
        ("compare", ("search", "evidence")), ("questions", ("comparison",)),
        ("plan", ("questions", "evidence")), ("report", ("comparison", "questions", "plan")),
    ):
        arguments = [command, "--config", config, "--output-dir", directory / command]
        for parent in parents:
            arguments += ["--" + parent, directory / ("compare" if parent == "comparison" else parent)]
        if command == "plan" and previous is not None:
            arguments += ["--previous-plan", previous["directory"] / "plan",
                          "--previous-config", previous["config"]]
        invoke(arguments)
    return {"directory": directory, "config": config,
            "plan": read_json(directory / "plan/plan.json")}


def _input_ref(identifier, role, path, *, format_name="json"):
    return {"input_id": identifier, "role": role, "path": str(path.resolve()),
            "sha256": checksum(path), "format": format_name,
            "stages": ["planning" if role == "plan_history" else "observation_normalization"]}


def _pin_history(config, previous):
    config["input_refs"] = [row for row in config["input_refs"]
                            if row["role"] not in {"plan_history", "observation_history"}]
    for label, path in (("configuration", previous["config"]),
                        ("manifest", previous["directory"] / "plan/run-manifest.json"),
                        ("audit", previous["directory"] / "plan/plan.json")):
        config["input_refs"].append(_input_ref("history-plan-" + label, "plan_history", path))
    config["input_refs"].append(_input_ref("history-observations", "observation_history",
        previous["directory"] / "observations/observations.jsonl", format_name="jsonl"))


def _completed_declarations(baseline):
    study = resolve_study(baseline["config"])
    search, _, _ = load_stage_records(study, baseline["directory"] / "search", "search")
    generated = read_jsonl(baseline["directory"] / "questions/question-records.jsonl")
    catalog = {row["record_id"]: deepcopy(row) for row in [*study["config"]["catalogs"], *search]}
    for row in generated:
        if row["record_type"] not in {"comparison_question", "comparison_requirement_slot", "comparison_action"}:
            catalog[row["record_id"]] = deepcopy(row)
    performed, mapping = [], {}
    for row in generated:
        if row["record_type"] != "comparison_action":
            continue
        executed = deepcopy(row)
        executed["experiment_id"] = "artificial-completed-execution-v1"
        executed["memberships"] = []
        identifier = record_identity(executed)
        mapping[row["record_id"]] = identifier
        executed.update(record_id=identifier, action_id=identifier)
        performed.append(executed)
    for row in performed:
        row["dependency_action_ids"] = [mapping[identifier] for identifier in row["dependency_action_ids"]]
        catalog[row["record_id"]] = row
    validate_catalog(list(catalog.values()))
    return catalog, performed, mapping


def _measurement_sources(directory, catalog, performed):
    template = read_json(ROOT / "examples/comparison/workflow/extraction-map.json")["rows"][0]["normalization"]
    raw, entries, declared = [], [], []
    for action in sorted(performed, key=lambda row: row["record_id"]):
        if action["operation_type"] not in {"measurement_panel", "mock_control"}:
            continue
        for pair in sorted(action["endpoint_locus_pairs"], key=lambda row: (row["locus_id"], row["endpoint_id"])):
            target = next(row for row in action["target_sequences"] if row["locus_id"] == pair["locus_id"])
            sequences = [catalog[identifier] for identifier in target["sequence_configuration_ids"]]
            sequence = min(sequences, key=lambda row: (not bool(row["variants"]), row["record_id"]))
            donor = catalog[action["donor_id"]]
            label = "Artificial software test: " + action["operation_type"]
            source = {"measurement_key": "artificial-completed-measurement-" + str(len(raw) + 1),
                "source_study": "Artificial software test; no real specimens or biological results",
                "source_experiment": action["experiment_id"], "source_replicate": str(action["replicate_position"]),
                "raw_value": "20", "raw_unit": "percent", "source_reported_label": label}
            normalization = deepcopy(template)
            normalization.update(design_id=deepcopy(action["design_id"]),
                unresolved_source_design=known("Explicit artificial design-neutral mock") if action["design_id"]["value"] is None
                    else unknown("The artificial treated design is explicitly resolved"),
                context_id=known(action["context_id"]), endpoint_id=pair["endpoint_id"],
                sequence_configuration_id=known(sequence["record_id"]), locus_ids=[pair["locus_id"]],
                executed_action_id=known(action["record_id"]), assay_id=known(action["assay_capability_id"]),
                control_status="passed", qc_status="passed", independent_unit_id=deepcopy(donor["independence_group_id"]),
                relationship_group=deepcopy(donor["independence_group_id"]), pairing_key=known(donor["canonical_label"]),
                source_reported_outcome=known("Artificial numeric test measurement"),
                adjudicated_outcome="observed_positive", adjudication_reason="Constructed positive software fixture",
                subject={"level": "specimen", "resolution_status": "resolved", "reason": None,
                    "source_subject_id": known("artificial-completed:" + donor["canonical_label"]),
                    "donor_id": known(action["donor_id"]), "specimen_id": known(action["specimen_id"])})
            raw.append(source)
            entries.append({key: source[key] for key in ("measurement_key", "source_study", "source_experiment", "source_replicate", "raw_unit")} |
                {"locator": {"kind": "file", "member": "rows/" + str(len(raw))}, "expected_raw_value": known("20"),
                 "source_reported_label": known(label), "normalization": normalization,
                 "normalization_unavailable_reason": None, "assertions": []})
            declared.append({"measurement_key": source["measurement_key"], "action_id": action["record_id"],
                             "operation_type": action["operation_type"], **pair,
                             "sequence_configuration_id": sequence["record_id"]})
    source_path = _write(directory / "artificial-measurements.json", raw)
    artifact = make_record("source_artifact", sha256=checksum(source_path),
        source_id="artificial-completed-panel-and-controls", release="software-test-v1",
        uri="urn:bes:artificial:completed-panel-and-controls", media_type="application/json", scope="complete_file",
        parent_artifact_id=unknown("Original artificial measurements for this software test"),
        locator={"kind": "file", "member": source_path.name}, access="synthetic",
        use_terms="Artificial software test values and specimens; no human data or biological inference")
    mapping = {"record_type": "comparison_extraction_map", "schema_version": "1",
        "source_artifact_id": artifact["record_id"], "artifact_format": "json",
        "extraction_version": "generic-json-artificial-v1", "mapping_version": "artificial-completed-panel-v1",
        "expected_count": len(entries), "rows": entries}
    map_path = _write(directory / "completed-extraction-map.json", mapping)
    return source_path, map_path, artifact, declared


def build_history_workflow(directory, *, stop_after_passed=False):
    """Run CLI stages with mock inputs."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=False)
    log, source = [], code_identity()
    inputs = directory / "initial-inputs"
    shutil.copytree(ROOT / "examples/comparison/workflow", inputs)
    policy = read_json(inputs / "experiment-policy.json")
    for rule in policy["rules"]:
        rule["allow_existing"] = True
    _write(inputs / "experiment-policy.json", policy)
    config = read_json(inputs / "study.json")
    next(row for row in config["input_refs"] if row["role"] == "experiment_policy")["sha256"] = checksum(inputs / "experiment-policy.json")
    _write(inputs / "study.json", config)
    baseline = _pipeline(inputs / "study.json", directory / "baseline", log)
    preserved = {path: checksum(path) for path in [baseline["config"], baseline["directory"] / "plan/plan.json",
        baseline["directory"] / "plan/run-manifest.json", baseline["directory"] / "observations/observations.jsonl"]}
    incoming = directory / "incoming-inputs"
    incoming.mkdir()
    catalog, performed, action_mapping = _completed_declarations(baseline)
    raw, map_path, artifact, declared = _measurement_sources(incoming, catalog, performed)
    catalog[artifact["record_id"]] = artifact
    current = deepcopy(resolve_study(baseline["config"])["config"])
    for row in current["input_refs"]:
        row["path"] = str((baseline["config"].parent / row["path"]).resolve())
    current["catalogs"] = list(catalog.values())
    current["input_refs"] += [_input_ref("incoming-artificial-raw", "raw_observations", raw),
                               _input_ref("incoming-artificial-map", "extraction_map", map_path)]
    _pin_history(current, baseline)
    current_path = _write(incoming / "study-passed.json", current)
    _write(incoming / "execution-declarations.json", {"scope": "Artificial software test only",
        "planned_to_performed_action_ids": action_mapping, "measurements": declared})
    passed = _pipeline(current_path, directory / "passed", log, previous=baseline)
    _write(directory / "command-log.json", log)
    if stop_after_passed:
        return {"baseline": baseline, "passed": passed, "source": source, "commands": log,
                "preserved": preserved, "measurements": declared}
    failed_mapping = read_json(map_path)
    failed_mapping["mapping_version"] = "artificial-completed-panel-v2-failed-control-review"
    failed_entry = next(row for row in failed_mapping["rows"] if row["normalization"]["design_id"]["value"] is None)
    failed_entry["normalization"]["qc_status"] = "failed"
    failed_map = _write(incoming / "failed-control-extraction-map.json", failed_mapping)
    corrected = deepcopy(current)
    ref = next(row for row in corrected["input_refs"] if row["input_id"] == "incoming-artificial-map")
    ref.update(path=str(failed_map.resolve()), sha256=checksum(failed_map))
    _pin_history(corrected, passed)
    corrected_path = _write(incoming / "study-failed-control.json", corrected)
    preserved.update({path: checksum(path) for path in [passed["config"], passed["directory"] / "plan/plan.json",
        passed["directory"] / "plan/run-manifest.json", passed["directory"] / "observations/observations.jsonl"]})
    failed = _pipeline(corrected_path, directory / "failed-control", log, previous=passed)
    _write(directory / "command-log.json", log)
    assert all(checksum(path) == expected for path, expected in preserved.items())
    assert code_identity() == source
    return {"baseline": baseline, "passed": passed, "failed": failed, "source": source, "commands": log,
            "preserved": preserved, "measurements": declared, "failed_measurement_key": failed_entry["measurement_key"]}


@pytest.fixture(scope="module")
def history_workflow(tmp_path_factory):
    return build_history_workflow(tmp_path_factory.mktemp("incoming-workflow") / "case")


def test_actual_imported_complete_panel_satisfies_the_contrast(history_workflow):
    baseline, passed = history_workflow["baseline"], history_workflow["passed"]
    assert baseline["plan"]["accounting"]["total_cost"] == 18
    assert {row["evidence_status"] for row in baseline["plan"]["questions"].values()} == {"no_measurements"}
    assert len(history_workflow["measurements"]) == 12
    assert passed["plan"]["unweighted_completion"] == 1
    assert passed["plan"]["accounting"]["total_cost"] == 0
    assert passed["plan"]["selected_action_ids"] == []
    assert {row["evidence_status"] for row in passed["plan"]["questions"].values()} == {"evaluable"}
    assert {row["planning_status"] for row in passed["plan"]["questions"].values()} == {"already_satisfied"}
    assert all("observation_id" in row for row in passed["plan"]["requirements"]["slot_assignments"].values())
    observations = read_jsonl(passed["directory"] / "observations/observations.jsonl")
    incoming = [row for row in observations if row["record_type"] == "comparison_source_observation"
                and row["source_study"].startswith("Artificial software test;")]
    assert len(incoming) == 12 and {row["raw_value"]["value"] for row in incoming} == {"20"}
    assert all("Artificial software test" in row["raw_row_text"] for row in incoming)
    source_ids = {row["record_id"] for row in incoming}
    assert {row["normalized_value"]["value"] for row in observations
            if row["record_type"] == "comparison_observation_revision" and row["source_record_id"] in source_ids} == {0.2}


def test_completed_neutral_mocks_do_not_become_treated_design_evidence(history_workflow):
    passed = history_workflow["passed"]
    rows = read_jsonl(passed["directory"] / "observations/observations.jsonl")
    mocks = {row["record_id"]: row for row in rows if row["record_type"] == "comparison_observation_revision"
             and row["design_id"]["value"] is None}
    assert len(mocks) == 4
    assert all(row["applicability"]["empirical_donor_credit"] for row in mocks.values())
    applications = [row for row in read_jsonl(passed["directory"] / "evidence/evidence-applicability.jsonl")
                    if row["observation_id"] in mocks]
    assert applications
    assert not any(row["assessment"]["empirical_donor_credit"] for row in applications)
    assert all(row["assessment"]["identity"]["status"] != "applicable" for row in applications)


def test_later_failed_control_does_not_borrow_a_passed_historical_revision(history_workflow):
    failed = history_workflow["failed"]
    rows = read_jsonl(failed["directory"] / "observations/observations.jsonl")
    source = next(row for row in rows if row["record_type"] == "comparison_source_observation"
                  and row["measurement_key"] == history_workflow["failed_measurement_key"])
    revisions = [row for row in rows if row["record_type"] == "comparison_observation_revision"
                 and row["source_record_id"] == source["record_id"]]
    assert {row["qc_status"] for row in revisions} == {"passed", "failed"}
    catalog = read_json(failed["config"])["catalogs"]
    active = next(row for row in active_observations([*catalog, *rows])
                  if row["source_record_id"] == source["record_id"])
    assert active["qc_status"] == "failed" and not active["applicability"]["empirical_donor_credit"]
    assert {row["evidence_status"] for row in failed["plan"]["questions"].values()} == {"partial_measurements"}
    assert failed["plan"]["meets_request"]
    assert failed["plan"]["accounting"]["total_cost"] > 0
    report = read_json(failed["directory"] / "report/plan.json")
    assert report == failed["plan"]
    assert all(row["evidence_status"] != "evaluable" for row in report["questions"].values())


def test_incoming_assessments_form_new_plan_revisions_without_rewriting_history(history_workflow):
    baseline, passed, failed = (history_workflow[name] for name in ("baseline", "passed", "failed"))
    assert passed["plan"]["previous_plan_id"] == baseline["plan"]["plan_id"]
    assert failed["plan"]["previous_plan_id"] == passed["plan"]["plan_id"]
    assert len({row["plan"]["plan_id"] for row in (baseline, passed, failed)}) == 3
    for stage in ("observations", "evidence", "compare", "questions", "plan", "report"):
        assert len({read_json(row["directory"] / stage / "run-manifest.json")["stage_key"]
                    for row in (baseline, passed, failed)}) == 3
    assert all(checksum(path) == expected for path, expected in history_workflow["preserved"].items())
    assert len(history_workflow["commands"]) == 24


MOCK_ADVERSARIES = ("missing_operation", "treatment_operation", "donor", "specimen", "context",
                    "source_experiment", "endpoint", "locus", "sequence", "failed_control")


def run_mock_adversary(history, directory, mutation):
    """Pass altered raw/map declarations through the actual observation command."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=False)
    passed = history["passed"]
    config = deepcopy(read_json(passed["config"]))
    catalog = {row["record_id"]: row for row in config["catalogs"]}
    source_ref = next(row for row in config["input_refs"] if row["input_id"] == "incoming-artificial-raw")
    map_ref = next(row for row in config["input_refs"] if row["input_id"] == "incoming-artificial-map")
    raw, mapping = read_json(source_ref["path"]), read_json(map_ref["path"])
    # Use the noncarrier's candidate-locus mock so another donor's alternate
    # sequence is available as an explicitly wrong same-locus sequence.
    for entry in mapping["rows"]:
        normal = entry["normalization"]
        if normal["design_id"]["value"] is not None:
            continue
        action = catalog[normal["executed_action_id"]["value"]]
        locus = normal["locus_ids"][0]
        targets = next(row["sequence_configuration_ids"] for row in action["target_sequences"] if row["locus_id"] == locus)
        alternatives = [row for row in catalog.values() if row["record_type"] == "comparison_sequence_configuration"
                        and row["locus_id"] == locus and row["record_id"] not in targets]
        if alternatives:
            break
    else:
        raise AssertionError("The artificial noncarrier mock requires a different donor's candidate sequence")
    if mutation == "missing_operation":
        normal["executed_action_id"] = unknown("No artificial performed operation was provided")
    elif mutation == "treatment_operation":
        normal["executed_action_id"] = known(next(row["record_id"] for row in catalog.values()
            if row["record_type"] == "comparison_action" and row["operation_type"] == "treatment"
            and row["donor_id"] == action["donor_id"]))
    elif mutation == "donor":
        normal["subject"]["donor_id"] = known(next(row["record_id"] for row in catalog.values()
            if row["record_type"] == "comparison_donor" and row["record_id"] != action["donor_id"]))
    elif mutation == "specimen":
        normal["subject"]["specimen_id"] = known(next(row["record_id"] for row in catalog.values()
            if row["record_type"] == "comparison_specimen" and row["record_id"] != action["specimen_id"]))
    elif mutation == "context":
        normal["context_id"] = known(next(row["record_id"] for row in catalog.values()
            if row["record_type"] == "comparison_experimental_context" and row["record_id"] != action["context_id"]))
    elif mutation == "source_experiment":
        entry["source_experiment"] = "different-artificial-executed-experiment"
        next(row for row in raw if row["measurement_key"] == entry["measurement_key"])["source_experiment"] = entry["source_experiment"]
    elif mutation == "endpoint":
        normal["endpoint_id"] = next(identifier for identifier in action["endpoint_ids"] if identifier != normal["endpoint_id"])
    elif mutation == "locus":
        normal["locus_ids"] = [next(identifier for identifier in action["locus_ids"] if identifier != locus)]
    elif mutation == "sequence":
        normal["sequence_configuration_id"] = known(alternatives[0]["record_id"])
    elif mutation == "failed_control":
        normal["control_status"] = "failed"
    else:
        raise AssertionError("Unknown artificial mock adversary")
    raw_path = _write(directory / "artificial-measurements.json", raw)
    original = catalog[mapping["source_artifact_id"]]
    artifact = make_record("source_artifact", **{key: value for key, value in original.items()
        if key not in {"record_type", "schema_version", "record_id", "source_artifact_id", "sha256", "source_id"}},
        sha256=checksum(raw_path), source_id="artificial-mock-boundary-" + mutation)
    config["catalogs"].append(artifact)
    mapping["source_artifact_id"] = artifact["record_id"]
    mapping["mapping_version"] = "artificial-mock-boundary-" + mutation
    map_path = _write(directory / "extraction-map.json", mapping)
    source_ref.update(path=str(raw_path.resolve()), sha256=checksum(raw_path))
    map_ref.update(path=str(map_path.resolve()), sha256=checksum(map_path))
    config["input_refs"] = [row for row in config["input_refs"] if row["role"] not in {"observation_history", "plan_history"}]
    path = _write(directory / "study.json", config)
    log = []
    expected = 2 if mutation == "locus" else 0
    result = _invoke(["observations", "--config", path, "--output", directory / "observations"], log, expected=expected)
    _write(directory / "command-log.json", log)
    if expected == 2:
        assert "sequence and reported locus" in result["reason"]
        assert not (directory / "observations").exists()
        checked = {"mutation": mutation, "status": "rejected_input", "reason": result["reason"]}
    else:
        rows = read_jsonl(directory / "observations/observations.jsonl")
        source = next(row for row in rows if row["record_type"] == "comparison_source_observation"
            and row["source_artifact_id"] == artifact["record_id"] and row["measurement_key"] == entry["measurement_key"])
        observation = next(row for row in rows if row["record_type"] == "comparison_observation_revision"
            and row["source_record_id"] == source["record_id"])
        assert not observation["applicability"]["empirical_donor_credit"]
        reason = "source_controls_not_passed" if mutation == "failed_control" else "source_design_unresolved"
        assert reason in observation["applicability"]["reason_codes"]
        checked = {"mutation": mutation, "status": "zero_empirical_credit",
                   "reason_codes": observation["applicability"]["reason_codes"]}
    _write(directory / "check-result.json", checked)
    return checked


@pytest.mark.parametrize("mutation", MOCK_ADVERSARIES)
def test_actual_importer_checks_each_completed_mock_binding(history_workflow, tmp_path, mutation):
    run_mock_adversary(history_workflow, tmp_path / mutation, mutation)
