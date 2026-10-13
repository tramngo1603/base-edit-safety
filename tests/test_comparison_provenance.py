"""Cache reuse must depend on scientific inputs rather than the whole budgeted run."""
from copy import deepcopy
import json
from pathlib import Path
import shutil

import pytest

from base_edit_safety.population.common import PopulationError, checksum
from base_edit_safety.population.comparison_provenance import (
    load_stage, resolve_study, save_stage, stage_key,
)


ROOT = Path(__file__).resolve().parents[1]
CODE = {"package_version": "test", "source_sha256": "a" * 64, "python": "test"}


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


@pytest.fixture
def study_path(tmp_path):
    source = ROOT / "examples/comparison/contracts"
    shutil.copytree(source, tmp_path / "inputs")
    return tmp_path / "inputs/study.json"


def resolve(path):
    return resolve_study(path, code=CODE)


def build_chain(study, root):
    outputs = {}
    dependencies = {
        "search": [], "observation_normalization": [],
        "evidence_applicability": ["search", "observation_normalization"],
        "comparison": ["search", "evidence_applicability"],
        "questions": ["comparison"], "planning": ["questions", "evidence_applicability"],
        "report": ["comparison", "questions", "planning"],
    }
    for stage, parents in dependencies.items():
        outputs[stage] = root / stage
        save_stage(study, stage, outputs[stage], {"check.txt": stage},
                   parents={p: outputs[p] for p in parents})
    return outputs


def test_budget_changes_reuse_search_and_reject_old_plan(study_path, tmp_path):
    old = resolve(study_path)
    outputs = build_chain(old, tmp_path / "cache")
    config = deepcopy(old["config"])
    config["policies"]["planning"]["total_budget"] += 1
    write_json(study_path, config)
    new = resolve(study_path)
    assert old["snapshot_run_id"] != new["snapshot_run_id"]
    assert stage_key(old, "search") == stage_key(new, "search")
    assert load_stage(new, outputs["search"], "search")[0]["check.txt"] == "search"
    load_stage(new, outputs["questions"], "questions")
    with pytest.raises(PopulationError, match="stale"):
        load_stage(new, outputs["planning"], "planning")
    with pytest.raises(PopulationError, match="stale"):
        load_stage(new, outputs["report"], "report")


def test_control_policy_invalidates_evidence_and_descendants(study_path, tmp_path):
    old = resolve(study_path)
    outputs = build_chain(old, tmp_path / "cache")
    config = deepcopy(old["config"])
    config["policies"]["evidence"]["policy_version"] = "changed-control-contract"
    write_json(study_path, config)
    new = resolve(study_path)
    load_stage(new, outputs["search"], "search")
    for stage in ["evidence_applicability", "comparison", "questions", "planning", "report"]:
        with pytest.raises(PopulationError, match="stale"):
            load_stage(new, outputs[stage], stage)


def test_catalog_order_does_not_change_search_key(study_path):
    old = resolve(study_path)
    config = deepcopy(old["config"])
    config["catalogs"].reverse()
    write_json(study_path, config)
    assert stage_key(old, "search") == stage_key(resolve(study_path), "search")


def test_source_change_invalidates_only_its_stage_dependencies(study_path, tmp_path):
    config = json.loads(study_path.read_text())
    raw = study_path.parent / "raw-observations.txt"
    raw.write_text("first immutable observation source\n")
    config["input_refs"].append({"input_id": "test-observations", "role": "raw_observations",
        "path": raw.name, "sha256": checksum(raw), "format": "text",
        "stages": ["observation_normalization"]})
    write_json(study_path, config)
    old = resolve(study_path)
    outputs = build_chain(old, tmp_path / "cache")
    raw.write_text("changed source values\n")
    with pytest.raises(PopulationError, match="checksum"):
        resolve(study_path)
    config["input_refs"][-1]["sha256"] = checksum(raw)
    write_json(study_path, config)
    new = resolve(study_path)
    load_stage(new, outputs["search"], "search")
    with pytest.raises(PopulationError, match="stale"):
        load_stage(new, outputs["observation_normalization"], "observation_normalization")


def test_cache_rejects_corrupted_bytes_and_unexpected_overwrite(study_path, tmp_path):
    study = resolve(study_path)
    output = tmp_path / "search"
    save_stage(study, "search", output, {"check.txt": "original"})
    with pytest.raises(PopulationError, match="exists"):
        save_stage(study, "search", output, {"check.txt": "replacement"})
    (output / "check.txt").write_text("changed")
    with pytest.raises(PopulationError, match="changed"):
        load_stage(study, output, "search")


def test_stage_rejects_bad_jsonl_and_input_overlap(study_path, tmp_path):
    study = resolve(study_path)
    with pytest.raises(PopulationError):
        save_stage(study, "search", tmp_path / "bad", {"records.jsonl": [
            {"record_type": "unregistered", "schema_version": "1", "record_id": "bad"}]})
    assert not (tmp_path / "bad").exists()
    with pytest.raises(PopulationError, match="overlap"):
        save_stage(study, "search", study_path.parent, {"check.txt": "bad"})


def test_missing_dependency_and_interrupted_stage_are_not_completed(study_path, tmp_path):
    study = resolve(study_path)
    with pytest.raises(PopulationError, match="dependencies"):
        save_stage(study, "planning", tmp_path / "plan", {"check.txt": "bad"})
    output = tmp_path / "search"
    save_stage(study, "search", output, {"check.txt": "original"})
    path = output / "run-manifest.json"
    manifest = json.loads(path.read_text())
    manifest["status"] = "interrupted"
    write_json(path, manifest)
    with pytest.raises(PopulationError, match="Incomplete"):
        load_stage(study, output, "search")


def test_source_change_after_resolution_is_detected_before_write(study_path, tmp_path):
    study = resolve(study_path)
    config = deepcopy(study["config"])
    config["policies"]["planning"]["total_budget"] += 1
    write_json(study_path, config)
    with pytest.raises(PopulationError, match="changed after resolution"):
        save_stage(study, "search", tmp_path / "search", {"check.txt": "bad"})


def test_unknown_source_stage_omission_cannot_make_cache_reusable(study_path):
    config = json.loads(study_path.read_text())
    raw = study_path.parent / "source.txt"
    raw.write_text("source\n")
    config["input_refs"].append({"input_id": "bad-stage", "role": "raw_observations",
        "path": raw.name, "sha256": checksum(raw), "format": "text", "stages": []})
    write_json(study_path, config)
    with pytest.raises(PopulationError):
        resolve(study_path)


def test_persisted_jsonl_cannot_bypass_registry_after_rehash(study_path, tmp_path):
    study = resolve(study_path)
    output = tmp_path / "search"
    row = deepcopy(study["config"]["catalogs"][0])
    save_stage(study, "search", output, {"records.jsonl": [row]})
    path = output / "records.jsonl"
    row["unregistered_scientific_claim"] = True
    path.write_text(json.dumps(row) + "\n")
    manifest_path = output / "run-manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["outputs"][path.name] = checksum(path)
    write_json(manifest_path, manifest)
    with pytest.raises(PopulationError):
        load_stage(study, output, "search")


def test_json_record_uses_the_same_registry(study_path, tmp_path):
    study = resolve(study_path)
    with pytest.raises(PopulationError):
        save_stage(study, "search", tmp_path / "bad-json", {
            "result.json": {"record_type": "fabricated_prediction_as_evidence",
                            "schema_version": "1", "record_id": "bad"}})


def test_different_parent_artifact_bytes_change_descendant_identity(study_path, tmp_path):
    study = resolve(study_path)
    first, second, third, observations = (tmp_path / name for name in
                                          ["search-a", "search-b", "search-c", "observations"])
    one = save_stage(study, "search", first, {"check.txt": "first actual result"})
    two = save_stage(study, "search", second, {"check.txt": "different actual result"})
    save_stage(study, "search", third, {"check.txt": "first actual result"})
    assert one["stage_key"] == two["stage_key"]
    save_stage(study, "observation_normalization", observations, {"check.txt": "source measurements"})
    children = [save_stage(study, "evidence_applicability", tmp_path / ("evidence-" + str(i)),
        {"check.txt": "derived fixture"},
        parents={"search": parent, "observation_normalization": observations})
        for i, parent in enumerate([first, second, third])]
    assert children[0]["stage_key"] != children[1]["stage_key"]
    assert children[0]["stage_key"] == children[2]["stage_key"]
