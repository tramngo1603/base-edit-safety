"""Evidence remains conditional on its measured endpoint, subject and context."""
from copy import deepcopy
import json
from pathlib import Path

import pytest

from base_edit_safety.population.common import PopulationError, checksum, known, read_json, unknown
from base_edit_safety.population.comparison_contracts import make_record, record_registry, validate_catalog
from base_edit_safety.population.comparison_evidence import (
    annotate_sequences, assess_observation, paired_descriptions,
    replay_published_alignments, source_summaries,
)
from base_edit_safety.population.comparison_provenance import resolve_study

ROOT = Path(__file__).resolve().parents[1]
ROWS = read_json(ROOT / "tests/fixtures/comparison/canonical-records.json")
POLICY = {"policy_version": "test", "require_controls": True, "require_qc": True,
    "require_detection_limit_for_zero": True, "allow_related_context": True, "allow_empirical_donor_credit": True}


def one(kind):
    return deepcopy(next(r for r in ROWS if r["record_type"] == "comparison_" + kind))


def remake(row, **changes):
    row = deepcopy(row)
    kind = row.pop("record_type")
    for key in ("schema_version", "record_id", record_registry()[kind], "definition_hash"):
        row.pop(key, None)
    row.update(changes)
    return make_record(kind, **row)


def assess(observation=None, extra=(), **changes):
    row = one("observation_revision") if observation is None else observation
    index = {r["record_id"]: r for r in [*ROWS, *extra]}
    args = {"design_id": one("design")["record_id"], "context_id": one("experimental_context")["record_id"],
        "endpoint_id": one("endpoint_definition")["record_id"], "donor_id": one("donor")["record_id"],
        "sequence_configuration_id": one("sequence_configuration")["record_id"],
        "mode": "research_sample_planning", "coverage_basis": "callable_supported", "policy": POLICY}
    args.update(changes)
    return assess_observation(row, index, **args)["assessment"]


def credited_observation():
    row = one("observation_revision")
    row["applicability"]["empirical_donor_credit"] = True
    return row


def test_exact_research_evidence_requires_full_target_and_supported_mode():
    row = credited_observation()
    assert assess(row)["empirical_donor_credit"]
    assert not assess(row, sequence_configuration_id=None)["empirical_donor_credit"]
    assert not assess(row, mode="benchmark_simulation", coverage_basis="retained_callset_model")["empirical_donor_credit"]
    assert not assess(row, donor_id=None)["empirical_donor_credit"]


def test_assessment_rejects_a_different_reported_hbg_copy():
    row = credited_observation()
    original = one("intended_locus")
    other = remake(original, gene="HBG2", locus_name="Artificial second HBG copy",
        start=original["start"] + 1000, end=original["end"] + 1000,
        path=[{**point, "pos0": point["pos0"] + 1000} for point in original["path"]])
    row["locus_ids"] = [other["record_id"]]
    with pytest.raises(PopulationError, match="sequence and reported locus"):
        assess(row, [other])


def test_abe710_cell_line_cannot_supply_abe8e_hspc_credit():
    source = credited_observation()
    editor = remake(one("editor_definition"), variant_name="ABE8e test variant")
    design = remake(one("design"), editor_definition_id=editor["record_id"])
    context = remake(one("experimental_context"), cell_type="HSPC", delivery_method="RNP electroporation")
    result = assess(source, [editor, design, context], design_id=design["record_id"], context_id=context["record_id"])
    assert not result["empirical_donor_credit"]
    assert result["identity"]["status"] == "related"
    assert result["context"]["status"] == "related"
    assert "different_complete_design" in result["reason_codes"]


def test_source_context_endpoint_remains_contextual_when_target_context_differs():
    source = credited_observation()
    endpoint = remake(one("endpoint_definition"), required_context_ids=[one("experimental_context")["record_id"]])
    source["endpoint_id"] = endpoint["record_id"]
    context = remake(one("experimental_context"), cell_type="HSPC", delivery_method="RNP electroporation")
    result = assess(source, [endpoint, context], endpoint_id=endpoint["record_id"], context_id=context["record_id"])
    assert result["classification"] == "related_design_context"
    assert result["endpoint"]["status"] == "related"
    assert not result["empirical_donor_credit"]
    strict = {**POLICY, "allow_related_context": False}
    result = assess(source, [endpoint, context], endpoint_id=endpoint["record_id"], context_id=context["record_id"], policy=strict)
    assert result["classification"] == "incompatible"


@pytest.mark.parametrize("limit", [unknown("Not reported"), known(0), known(1.1)])
def test_zero_without_usable_sensitivity_earns_no_credit(limit):
    row = credited_observation()
    row["normalized_value"] = known("0.000")
    row["detection_limit"] = {"value": limit, "unit": "fraction"}
    result = assess(row)
    assert not result["empirical_donor_credit"]
    assert "zero_without_usable_detection_limit" in result["reason_codes"]


def test_pooled_functional_value_cannot_resolve_a_donor_or_copy_specific_edit():
    row = credited_observation()
    endpoint = remake(one("endpoint_definition"), name="HbF_fraction", biological_meaning="HbF fraction in a pooled assay",
        denominator="Total hemoglobin", measurement_level="pooled_group", locus_resolution="not_applicable",
        allele_resolution="not_required", requires_sequence=False, required_context_ids=[])
    row.update(endpoint_id=endpoint["record_id"], sequence_configuration_id=unknown("Pooled functional assay"), locus_ids=[])
    row["subject"].update(level="pooled_group", donor_id=unknown("Pooled donors"), resolution_status="pooled", reason="Pooled donors")
    functional = assess(row, [endpoint], endpoint_id=endpoint["record_id"], sequence_configuration_id=None)
    assert not functional["empirical_donor_credit"]
    editing = assess(row, [endpoint])
    assert editing["endpoint"]["status"] == "incompatible"


def test_prediction_never_enters_reconciliation():
    with pytest.raises(PopulationError, match="measured observation"):
        assess(one("prediction"))


def test_repeated_technical_values_keep_one_biological_unit_and_conflicts():
    index = validate_catalog(ROWS)
    first = one("observation_revision")
    second = deepcopy(first)
    second["record_id"] = "test-second-measurement"
    second["measurement_kind"] = "technical_replicate"
    second["normalized_value"] = known(0.8)
    summary = source_summaries([first, second], index)[0]
    assert summary["measurements"] == 2
    assert summary["reported_independent_units"] == 1
    assert summary["source_donor_labels"] == 1
    assert summary["distinct_values_within_unit"] == [["0.2", "0.8"]]


def test_source_alignment_requires_actual_admission_and_exact_columns(tmp_path):
    source_path = tmp_path / "primary.txt"
    source_path.write_text("Artificial measured positive\n")
    artifact = remake(one("source_artifact"), sha256=checksum(source_path), source_id="Artificial alignment source")
    alignment = one("alignment")
    sequence = one("sequence_configuration")
    dc = one("design_configuration")
    guide = "".join(column["guide"] or "-" for column in alignment["columns"])
    target = "".join(column["target"] or "-" for column in alignment["columns"])
    positions = [point["pos0"] for point in sequence["path"]]
    expected = {"contig": sequence["contig"], "strand": sequence["strand"],
        "start": min(positions), "end": max(positions) + 1, "guide": guide, "target": target,
        "pam": alignment["pam_sequence"], "path": sequence["path"]}
    row = {"source_key": "artificial-positive", "design_id": one("design")["record_id"],
        "source_artifact_id": artifact["record_id"], "source_locator": {"kind": "file", "member": "line1"},
        "source_positive": True, "source_reported_label": "positive", "alignment": expected, "unresolved_reason": None}
    path = tmp_path / "published.json"
    mapping = {"schema_version": "1", "mapping_version": "test", "rows": [row]}
    path.write_text(json.dumps(mapping))
    config = read_json(ROOT / "examples/comparison/contracts/study.json")
    config["catalogs"].append(artifact)
    study = {"config": config, "sources": {"primary": {"role": "raw_observations", "path": str(source_path), "sha256": checksum(source_path)},
        "alignments": {"role": "published_alignments", "path": str(path), "sha256": checksum(path)}}}
    artifacts = {"records.jsonl": [alignment, sequence, dc], "search-audit.json": {"hits": []}}
    result = replay_published_alignments(study, artifacts)[0]
    assert result["source_positive"] and not result["exact_alignment_ids"]
    assert result["disposition"] != "recovered_exact_alignment"
    artifacts["search-audit.json"]["hits"] = [{"alignment_id": alignment["record_id"]}]
    result = replay_published_alignments(study, artifacts)[0]
    assert result["disposition"] == "recovered_exact_alignment"
    assert result["exact_alignment_ids"] == [alignment["record_id"]]
    assert not result["empirical_target_donor_credit"]
    # A reported RNA gap changes the source sequence/path and cannot be recovered
    # by an overlapping ungapped nomination.
    expected["target"] = "-" + target[1:]
    expected["end"] -= 1
    dna = expected["target"].replace("-", "") + expected["pam"]
    expected["path"] = [{"kind": "reference", "pos0": pos, "base": base}
                         for pos, base in zip(range(expected["start"], expected["end"]), dna)]
    path.write_text(json.dumps(mapping))
    study["sources"]["alignments"]["sha256"] = checksum(path)
    result = replay_published_alignments(study, artifacts)[0]
    assert result["disposition"] == "outside_declared_search_policy"
    assert "source_bulges_exceed_profile" in result["reason_codes"]
    assert not result["exact_alignment_ids"]
    assert result["related_alignment_ids"] == [alignment["record_id"]]
    row["design_id"] = None
    row["unresolved_reason"] = "The source target row does not uniquely identify its editor"
    path.write_text(json.dumps(mapping))
    study["sources"]["alignments"]["sha256"] = checksum(path)
    result = replay_published_alignments(study, artifacts)[0]
    assert result["source_positive"] and result["alignment"] == expected
    assert result["disposition"] == "unresolved_source_mapping"
    assert "source_design_unresolved" in result["reason_codes"]
    assert not result["exact_alignment_ids"] and not result["empirical_target_donor_credit"]
    row["unresolved_reason"] = None
    path.write_text(json.dumps(mapping))
    study["sources"]["alignments"]["sha256"] = checksum(path)
    with pytest.raises(PopulationError, match="source design requires"):
        replay_published_alignments(study, artifacts)


def test_paired_descriptions_require_same_regimen_and_distinct_unambiguous_arms():
    first = one("observation_revision")
    first.update(pairing_key=known("same-donor-experiment"),
                 normalized_value=known("0.123456789012345678901234567891"))
    second_editor = remake(one("editor_definition"), variant_name="Second artificial editor")
    second_design = remake(one("design"), editor_definition_id=second_editor["record_id"])
    second = deepcopy(first)
    second.update(record_id="second-artificial-measurement", design_id=known(second_design["record_id"]),
                  normalized_value=known("0.123456789012345678901234567890"))
    context = remake(one("experimental_context"), condition="Different treatment dose")
    index = {row["record_id"]: row for row in [*ROWS, second_editor, second_design, context]}
    designs = {first["design_id"]["value"], second_design["record_id"]}
    result = paired_descriptions([first, second], index, designs, POLICY)
    assert len(result) == 1
    assert result[0]["left_minus_right"] in {"0.000000000000000000000000000001", "-0.000000000000000000000000000001"}
    second["context_id"] = known(context["record_id"])
    assert paired_descriptions([first, second], index, designs, POLICY) == []
    second["context_id"] = first["context_id"]
    repeat = deepcopy(second)
    repeat["record_id"] = "third-artificial-measurement"
    assert paired_descriptions([first, second, repeat], index, designs, POLICY) == []


def test_local_regulatory_features_use_exact_path_coordinates_and_version(tmp_path):
    config = read_json(ROOT / "examples/comparison/contracts/study.json")
    config["input_refs"] = []
    sequence = one("sequence_configuration")
    position = sequence["path"][0]["pos0"]
    annotation = {"schema_version": "1", "assembly": one("reference_bundle")["assembly"],
        "source": "Artificial regulatory fixture", "release": "test-1", "features": [
            {"feature_type": "regulatory", "name": "included", "contig": sequence["contig"], "start": position, "end": position + 1},
            {"feature_type": "regulatory", "name": "excluded", "contig": sequence["contig"], "start": position - 1, "end": position}]}
    feature_path = tmp_path / "features.json"
    feature_path.write_text(json.dumps(annotation))
    config["input_refs"].append({"input_id": "local-features", "role": "annotations", "path": feature_path.name,
        "sha256": checksum(feature_path), "format": "json", "stages": ["comparison"]})
    config_path = tmp_path / "study.json"
    config_path.write_text(json.dumps(config))
    result = annotate_sequences(resolve_study(config_path), {"records.jsonl": [sequence]})
    assert [row["name"] for row in result["sequence-annotations.json"][0]["features"]] == ["included"]
    assert result["annotation-manifest.json"] == [{"source": annotation["source"], "release": "test-1",
        "sha256": checksum(feature_path), "assembly": annotation["assembly"]}]
    annotation["release"] = ""
    feature_path.write_text(json.dumps(annotation))
    config["input_refs"][0]["sha256"] = checksum(feature_path)
    config_path.write_text(json.dumps(config))
    with pytest.raises(PopulationError, match="release"):
        annotate_sequences(resolve_study(config_path), {"records.jsonl": [sequence]})
