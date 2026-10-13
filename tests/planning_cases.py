"""Planner edge cases built from the mock full panel."""
from copy import deepcopy
import importlib.util
import json
from pathlib import Path

from base_edit_safety.population.comparison_contracts import make_record, record_identity, record_registry, validate_catalog
from base_edit_safety.population.comparison_provenance import resolve_study
from base_edit_safety.population.comparison_questions import prepare_question_set


ROOT = Path(__file__).resolve().parents[1]


def full_fixture(directory, **kwargs):
    spec = importlib.util.spec_from_file_location("section9_source_fixture", ROOT / "examples/comparison/questions/build_fixture.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    fixture = module.build_fixture(directory, **kwargs)
    fixture["question_set"] = prepare_question_set(fixture["study"], fixture["records"])
    return fixture


def remake(row, **changes):
    values = deepcopy(row)
    for key in ("record_type", "schema_version", "record_id", record_registry()[row["record_type"]], "definition_hash"):
        values.pop(key, None)
    values.update(changes)
    return make_record(row["record_type"], **values)


def replace(value, mapping):
    if isinstance(value, str):
        return mapping.get(value, value)
    if isinstance(value, list):
        return [replace(row, mapping) for row in value]
    if isinstance(value, dict):
        return {key: replace(row, mapping) for key, row in value.items()}
    return value


def rekey(records):
    records = deepcopy(records)
    for _ in range(40):
        mapping = {row["record_id"]: record_identity(row) for row in records if row["record_id"] != record_identity(row)}
        if not mapping:
            return list({row["record_id"]: row for row in records}.values())
        records = replace(records, mapping)
        for row in records:
            if row["record_type"] == "comparison_design":
                row["definition_hash"] = row["record_id"].rsplit(":", 1)[1]
    raise AssertionError("Synthetic fixture identity graph did not converge")


def write_case(fixture, directory, records, *, budget=18):
    """Freeze an explicitly authored variant with the same literal genome truth."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    rows = rekey(records)
    for row in rows:
        if row["record_type"] == "comparison_action":
            row["memberships"] = [{"question_id": slot["question_id"], "slot_id": slot["slot_id"]}
                                  for slot in rows if slot["record_type"] == "comparison_requirement_slot"
                                  and row["action_id"] in slot["allowed_action_ids"]]
    validate_catalog(rows)
    config = deepcopy(fixture["study"]["config"])
    config["catalogs"] = rows
    config["study_id"] = "declared-planner-case-" + directory.name
    # The cases below supply their own questions, resources and costs, so remove
    # the base example's generator policy.
    config["input_refs"] = [row for row in config["input_refs"] if row["role"] != "experiment_policy"]
    for ref in config["input_refs"]:
        ref["path"] = fixture["study"]["sources"][ref["input_id"]]["path"]
    config["policies"]["questions"]["policy_version"] = "hand-authored-adversarial-v1"
    config["policies"]["questions"]["templates"] = sorted({row["template"] for row in rows if row["record_type"] == "comparison_question"})
    config["policies"]["planning"]["total_budget"] = budget
    (directory / "study.json").write_text(json.dumps(config, indent=2) + "\n")
    (directory / "records.json").write_text(json.dumps(rows, indent=2) + "\n")
    study = resolve_study(directory / "study.json")
    return prepare_question_set(study, [])


def split_design_questions(fixture, *, required=True):
    qset = fixture["question_set"]
    original = next(iter(qset.questions.values()))
    rows = [deepcopy(row) for row in qset.records if row["record_type"] not in {"comparison_question", "comparison_requirement_slot"}]
    for design in sorted(original["design_ids"]):
        question = remake(original, template="candidate_validation", hypothesis="Separate complete validation for " + design,
            design_ids=[design], arms=[{"arm_id": "resolved", "membership": "all_resolved", "required_independent_donors": 1}],
            pairing="none", hard_required=required,
            genotype_predicate={**deepcopy(original["genotype_predicate"]), "kind": "any_resolved", "variants": [],
                                "sequence_configuration_id": {"value": None, "reason": "Any resolved actual donor sequence"}})
        slots = [remake(slot, question_id=question["record_id"], arm_id="resolved") for slot in qset.slots.values()
                 if slot["design_id"]["value"] in {None, design}]
        question["requirement_slot_ids"] = [slot["slot_id"] for slot in slots]
        rows.extend([question, *slots])
    return rows


def joint_seed_case(fixture, directory):
    """Each cheap C panel fills a capacity-one instance; both N panels share5."""
    rows = split_design_questions(fixture)
    resources = {row["record_id"]: row for row in rows if row["record_type"] == "comparison_resource"}
    mapping, replacements = {}, []
    donors = fixture["labels"]["donors"]
    shared = {}
    for rid, resource in resources.items():
        if resource["resource_type"] != "assay_library":
            new = remake(resource, cost=0)
        else:
            donor = resource["compatibility"]["donor_ids"][0]
            if donor not in shared:
                carrier = donor == donors["C"]
                shared[donor] = remake(resource, cost=2 if carrier else 5,
                    sharing_key="finite-shared-panel:" + donor,
                    capacities=[{"dimension": "uses", "capacity": {"kind": "finite", "quantity": 1 if carrier else 2}}],
                    compatibility={**resource["compatibility"], "design_ids": []})
            new = shared[donor]
        mapping[rid] = new["record_id"]
        replacements.append(new)
    rows = replace([row for row in rows if row["record_type"] != "comparison_resource"], mapping)
    rows.extend({row["record_id"]: row for row in replacements}.values())
    return write_case(fixture, directory, rows, budget=5)


def alternatives_case(fixture, directory):
    """All six measured panels need two real sequencing instances of capacity3."""
    rows = deepcopy(fixture["question_set"].records)
    original = next(row for row in rows if row["record_type"] == "comparison_resource" and row["resource_type"] == "sequencing_batch")
    first = remake(original, cost=3, sharing_key="finite-sequencing-instance-one",
                   pool_id={"value": None, "reason": "Separately priced finite standalone instance"}, instance_index=0, max_pool_instances=1,
                   capacities=[{"dimension": "uses", "capacity": {"kind": "finite", "quantity": 3}}])
    second = remake(original, cost=5, sharing_key="finite-sequencing-instance-two",
                    pool_id={"value": None, "reason": "Separately priced finite standalone instance"}, instance_index=0, max_pool_instances=1,
                    capacities=[{"dimension": "uses", "capacity": {"kind": "finite", "quantity": 3}}])
    rows = [row for row in rows if row["record_id"] != original["record_id"]] + [first, second]
    for row in rows:
        if row["record_type"] != "comparison_action":
            continue
        row["resource_demands"] = [({"eligible_resource_ids": [first["record_id"], second["record_id"]],
            "dimension": demand["dimension"], "quantity": demand["quantity"], "assignment_key": "sequencing"}
            if demand.get("resource_id") == original["record_id"] else demand) for demand in row["resource_demands"]]
    return write_case(fixture, directory, rows, budget=26)


def optional_trap_case(fixture, directory):
    """A high-priority new experiment is affordable alone but not after seed18."""
    qset = fixture["question_set"]
    rows = deepcopy(qset.records)
    labels = fixture["labels"]
    donor, design = labels["donors"]["C"], labels["designs"]["A"]
    existing_treatment = qset.actions[labels["actions"]["C:A:treatment"]]
    existing_panel = qset.actions[labels["actions"]["C:A:measurement_panel"]]
    existing_mock = qset.actions[labels["actions"]["C:mock:mock_control"]]
    new_resources, resource_map = [], {}
    for old in (existing_treatment, existing_panel, existing_mock):
        for demand in old["resource_demands"]:
            if "resource_id" not in demand:
                continue
            resource = qset.resources[demand["resource_id"]]
            if resource["resource_type"] not in {"treatment", "assay_library", "mock_control"}:
                continue
            fresh = remake(resource, sharing_key=resource["sharing_key"] + ":optional-new-experiment", cost=1)
            resource_map[resource["record_id"]] = fresh["record_id"]
            new_resources.append(fresh)
    treatment = remake(existing_treatment, experiment_id="optional-new-experiment", memberships=[],
                       resource_demands=replace(existing_treatment["resource_demands"], resource_map))
    mock = remake(existing_mock, experiment_id="optional-new-experiment", memberships=[],
                  resource_demands=replace(existing_mock["resource_demands"], resource_map))
    panel = remake(existing_panel, experiment_id="optional-new-experiment", memberships=[],
        resource_demands=replace(existing_panel["resource_demands"], resource_map),
        dependency_action_ids=replace(existing_panel["dependency_action_ids"], {
            existing_treatment["record_id"]: treatment["record_id"], existing_mock["record_id"]: mock["record_id"]}))
    original = next(iter(qset.questions.values()))
    question = remake(original, template="candidate_validation", hypothesis="High-priority separate declared experiment",
        design_ids=[design], arms=[{"arm_id": "resolved", "membership": "all_resolved", "required_independent_donors": 1}],
        pairing="none", priority_weight=1000, hard_required=False,
        genotype_predicate={**deepcopy(original["genotype_predicate"]), "kind": "any_resolved", "variants": [],
                            "sequence_configuration_id": {"value": None, "reason": "Any resolved donor sequence"}})
    slots = [remake(slot, question_id=question["record_id"], arm_id="resolved",
        allowed_action_ids=replace(slot["allowed_action_ids"], {
            existing_panel["record_id"]: panel["record_id"], existing_mock["record_id"]: mock["record_id"]}))
        for slot in qset.slots.values() if slot["donor_id"] == donor and slot["design_id"]["value"] in {None, design}]
    question["requirement_slot_ids"] = [row["record_id"] for row in slots]
    # The original six-panel sequencing capacity is explicit; an extra panel is
    # also explicitly admitted in this separate cost experiment.
    for row in rows:
        if row["record_type"] == "comparison_resource" and row["resource_type"] == "sequencing_batch":
            row["capacities"][0]["capacity"]["quantity"] = 8
    rows.extend([*new_resources, treatment, mock, panel, question, *slots])
    return write_case(fixture, directory, rows, budget=18)
