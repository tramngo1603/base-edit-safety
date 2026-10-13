"""Build a frozen coverage universe and explicit donor-site assay actions."""
from __future__ import annotations

from .common import PopulationError, identity, record, unique
from .contracts import validate_request
from .variants import canonical_donors


def make_units(config, candidates, run_id):
    assay = config["assay_policy"]
    rows = []
    for c in candidates:
        if c["role"] != "candidate_offtarget":
            continue
        payload = {"configuration_id": c["configuration_id"], "guide_definition_id": c["guide_definition_id"],
                   "editor_definition_id": c["editor_definition_id"], "biological_context_id": assay["biological_context_id"],
                   "biological_context": assay["context"],
                   "measurement_class": assay["measurement_class"], "allele_resolution_required": True,
                   "locus_resolution_required": True, "extra_context_bases": assay["extra_context_bases"]}
        if assay.get("max_detection_limit_fraction") is not None:
            payload["max_detection_limit_fraction"] = assay["max_detection_limit_fraction"]
        rows.append(record("coverage_unit", run_id, id=identity("unit", payload), **payload,
                           weight=config["planning_policy"]["default_weight"],
                           required_donors=config["planning_policy"]["default_required_donors"]))
    return sorted(rows, key=lambda r: r["id"])


def universe_id(config, units):
    return identity("universe", {"coverage_basis": config["coverage_basis"],
        "units": [{k: v for k, v in u.items() if k not in {"record_type", "schema_version", "run_id", "provenance"}}
                  for u in sorted(units, key=lambda r: r["id"])],
        "policy": "all_admitted_offtarget_configurations_v1"})


def make_actions(config, units, candidates, donors, carriage, run_id):
    candidates = unique(candidates, "configuration_id")
    audit = {(r["donor_id"], r["configuration_id"]): r for r in carriage}
    if len(audit) != len(carriage):
        raise PopulationError("Duplicate donor-carriage row")
    actions = []
    assay = config["assay_policy"]
    for donor in sorted(donors, key=lambda d: (d["donor_id"], d["specimen_id"])):
        for unit in units:
            c = candidates[unit["configuration_id"]]
            row = audit.get((donor["donor_id"], c["configuration_id"]))
            reasons = []
            eligibility = "eligible"
            if row is None:
                raise PopulationError("Unit donor-carriage link does not resolve")
            status = row["model_carriage"] if config["coverage_basis"] == "retained_callset_model" else row["evidence_supported_carriage"]
            positive = status in {"present", "confirmed"}
            if not positive:
                reasons.append("carriage_" + status)
                eligibility = "ineligible"
            if not donor["use_permitted"]:
                reasons.append("use_not_permitted")
                eligibility = "ineligible"
            if config["mode"] == "research_sample_planning":
                if donor["availability"] != "verified_available":
                    reasons.append("specimen_" + donor["availability"])
                    eligibility = "ineligible"
            elif not donor["simulation_eligible"]:
                reasons.append("scenario_excludes_donor")
                eligibility = "ineligible"
            if donor["specimen_context_id"] != assay["specimen_context_id"]:
                reasons.append("specimen_context_mismatch")
                eligibility = "ineligible"
            capability_ok = (assay["allele_resolved"] and assay["locus_resolved"]
                             and assay["feasibility"] == "declared_feasible")
            if c.get("guide_id") not in assay["guide_ids"]:
                reasons.append("assay_guide_not_declared")
                capability_ok = False
            locus = {"contig": c["contig"], "start": c["bounding_interval"][0], "end": c["bounding_interval"][1], "strand": c["strand"]}
            if assay["locus_scope"]["mode"] == "explicit_loci" and locus not in assay["locus_scope"]["loci"]:
                reasons.append("assay_locus_not_declared")
                capability_ok = False
            if not capability_ok:
                reasons.append("assay_requires_review" if assay["feasibility"] != "infeasible" else "assay_infeasible")
                if eligibility == "eligible":
                    eligibility = "requires_review" if assay["feasibility"] != "infeasible" else "ineligible"
            payload = {"donor_id": donor["donor_id"], "specimen_id": donor["specimen_id"],
                       "configuration_id": c["configuration_id"], "assay_id": assay["assay_id"],
                       "resolution_contract_id": assay["resolution_contract_id"],
                       "biological_context_id": assay["biological_context_id"]}
            actions.append(record("assay_action", run_id, id=identity("action", payload), **payload,
                locus={"contig": c["contig"], "interval": c["bounding_interval"], "strand": c["strand"]},
                guide_definition_id=c["guide_definition_id"], editor_definition_id=c["editor_definition_id"],
                coverage_basis=config["coverage_basis"], eligibility=eligibility, eligibility_reasons=reasons,
                assay_cost=assay["assay_cost"], covered_units=[unit["id"]] if positive and capability_ok else [],
                possible_only_units=[unit["id"]] if status == "possible" else [],
                interpretation="proposed_assay" if config["mode"] == "research_sample_planning" else "simulation_assumption"))
    return sorted(actions, key=lambda a: a["id"])


def make_request(config, units, donors, actions, evidence):
    p = config["planning_policy"]
    people = []
    for donor in canonical_donors(donors):
        specimens = [d for d in donors if d["donor_id"] == donor["donor_id"]
                     and d["specimen_context_id"] == config["assay_policy"]["specimen_context_id"]]
        availability = {d["availability"] for d in specimens if d["use_permitted"]}
        people.append({"id": donor["donor_id"], "activation_cost": donor["activation_cost"],
            "specimen_context_id": config["assay_policy"]["specimen_context_id"],
            "availability": "verified_available" if "verified_available" in availability else
                            "unverified" if "unverified" in availability else "unavailable",
            "use_permitted": any(d["use_permitted"] for d in specimens),
            "simulation_eligible": any(d["use_permitted"] and d["simulation_eligible"] for d in specimens)})
    request = {"record_type": "donor_panel_request", "schema_version": "1", "mode": config["mode"],
               "coverage_basis": config["coverage_basis"], "protocol_id": identity("protocol", config),
               "universe_id": universe_id(config, units), "cost_unit": p["cost_unit"],
               "budget": {"max_cost": p["max_cost"]-p["fixed_overhead"], "max_donors": p["max_donors"], "max_assays": p["max_assays"]},
               "units": [{k: u[k] for k in ("id", "configuration_id", "weight", "required_donors")} for u in units],
               "donors": people,
               "actions": [{k: a[k] for k in ("id", "donor_id", "assay_cost", "covered_units", "possible_only_units", "eligibility", "resolution_contract_id")} for a in actions],
               "existing_usable_evidence": [{"unit_id": e["unit_id"], "donor_id": e["donor_id"], "evidence_id": e["evidence_id"]}
                   for e in evidence if e["usable"] and config["coverage_basis"] != "retained_callset_model"]}
    return validate_request(request)
