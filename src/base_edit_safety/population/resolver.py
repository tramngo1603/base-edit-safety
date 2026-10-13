"""Resolve planner source records and recompute every action and prior-evidence link."""
from __future__ import annotations

from pathlib import Path

from .actions import make_actions, make_request, make_units
from .common import PopulationError, checksum, fields, read_json, record
from .contracts import validate_request
from .evidence import reconcile, resolve_measurements
from .provenance import load_stage, resolve_file, resolve_study


def resolve_request(request_path, manifest_path):
    request_path, manifest_path = Path(request_path).resolve(), Path(manifest_path).resolve()
    request, manifest = read_json(request_path), read_json(manifest_path)
    validate_request(request)
    fields(manifest, {"record_type", "schema_version", "run_id", "mode", "coverage_basis", "request",
                      "synthetic_declaration", "source_population_manifest", "artifacts"}, label="planning manifest")
    if manifest["record_type"] != "donor_panel_manifest" or manifest["schema_version"] != "1":
        raise PopulationError("Unsupported planning manifest")
    if request_path != resolve_file(manifest["request"], manifest_path.parent):
        raise PopulationError("Manifest refers to a different request")
    for key in ("mode", "coverage_basis"):
        if request[key] != manifest[key]:
            raise PopulationError("Planning manifest mode or coverage basis mismatch")
    resolved = {"request": checksum(request_path), "manifest": checksum(manifest_path)}
    if manifest["synthetic_declaration"] is not None:
        if request["mode"] != "synthetic_test" or request["coverage_basis"] != "synthetic_truth" \
                or manifest["source_population_manifest"] is not None or manifest["artifacts"]:
            raise PopulationError("Synthetic declaration cannot replace real study provenance")
        path = resolve_file(manifest["synthetic_declaration"], manifest_path.parent)
        declaration = read_json(path)
        fields(declaration, {"record_type", "schema_version", "real_specimens", "real_measurements", "coverage_unit_ids",
                            "configuration_ids", "resolution_contract_ids", "specimen_context_id", "assumptions"},
               label="synthetic declaration")
        if declaration["record_type"] != "synthetic_planning_declaration" or declaration["schema_version"] != "1" \
                or declaration["real_specimens"] is not False or declaration["real_measurements"] is not False:
            raise PopulationError("Synthetic fixtures must explicitly exclude real specimens and measurements")
        for field, expected in (("coverage_unit_ids", {u["id"] for u in request["units"]}),
                                ("configuration_ids", {u["configuration_id"] for u in request["units"]}),
                                ("resolution_contract_ids", {a["resolution_contract_id"] for a in request["actions"]})):
            if not isinstance(declaration[field], list) or set(declaration[field]) != expected or len(declaration[field]) != len(expected):
                raise PopulationError("Synthetic declaration does not resolve every request identifier")
        if not declaration["assumptions"] or any(d["specimen_context_id"] != declaration["specimen_context_id"] for d in request["donors"]):
            raise PopulationError("Synthetic context or assumptions missing")
        if request["existing_usable_evidence"]:
            raise PopulationError("This synthetic declaration has no source evidence records")
        resolved["synthetic_declaration"] = checksum(path)
        return request, manifest, resolved
    refs = {}
    for ref in manifest["artifacts"]:
        if ref.get("role") in refs:
            raise PopulationError("Duplicate manifest artifact role")
        refs[ref.get("role")] = resolve_file(ref, manifest_path.parent)
    expected_roles = {"config", "coverage_units", "donor_carriage", "donor_inventory", "assay_resolution",
                      "action_catalog", "evidence_links", "candidate_configurations", "search_manifest", "evidence_manifest"}
    if set(refs) != expected_roles or manifest["source_population_manifest"] is None:
        raise PopulationError("Planner needs all canonical upstream artifacts")
    pop_path = resolve_file(manifest["source_population_manifest"], manifest_path.parent)
    study = resolve_study(refs["config"])
    if study["run_id"] != manifest["run_id"]:
        raise PopulationError("Planning manifest has a stale study identity")
    population, _ = load_stage(study, pop_path.parent, "build")
    search, _ = load_stage(study, refs["search_manifest"].parent, "search")
    evidence, _ = load_stage(study, refs["evidence_manifest"].parent, "reconcile")
    for role, directory, filename in (
        ("coverage_units", refs["evidence_manifest"].parent, "coverage-units.jsonl"),
        ("evidence_links", refs["evidence_manifest"].parent, "evidence-links.jsonl"),
        ("donor_carriage", refs["search_manifest"].parent, "donor-carriage.jsonl"),
        ("candidate_configurations", refs["search_manifest"].parent, "candidate-configurations.jsonl"),
        ("donor_inventory", pop_path.parent, "donor-inventory.jsonl"),
        ("assay_resolution", pop_path.parent, "assay-resolution.json")):
        if refs[role] != directory/filename:
            raise PopulationError("Planner role is detached from its upstream stage")
    config = study["config"]
    if population["assay-resolution.json"] != record("assay_resolution_contract", study["run_id"], **config["assay_policy"]):
        raise PopulationError("Assay capability artifact differs from the study declaration")
    candidates, donors = search["candidate-configurations.jsonl"], population["donor-inventory.jsonl"]
    units = make_units(config, candidates, study["run_id"])
    if units != evidence["coverage-units.jsonl"]:
        raise PopulationError("Coverage universe or context was modified")
    paths = study["paths"]
    measurements = resolve_measurements(read_json(paths["candidate_ledger"]) if paths["candidate_ledger"] else None,
        read_json(paths["measurement_linkage"]) if paths["measurement_linkage"] else None,
        study["sources"].get("candidate_ledger", {}).get("sha256"))
    links = reconcile(config, units, candidates, population["guides.json"], measurements, study["run_id"])
    if links != evidence["evidence-links.jsonl"]:
        raise PopulationError("Evidence applicability differs from its source measurements")
    actions = make_actions(config, units, candidates, donors, search["donor-carriage.jsonl"], study["run_id"])
    from .common import read_jsonl
    if actions != read_jsonl(refs["action_catalog"]):
        raise PopulationError("Action eligibility or carriage coverage was modified")
    expected = make_request(config, units, donors, actions, links)
    if request != expected:
        raise PopulationError("Planner request differs from resolved costs, context, coverage or evidence")
    resolved.update({role: checksum(path) for role, path in refs.items()})
    resolved["source_population_manifest"] = checksum(pop_path)
    return request, manifest, resolved
