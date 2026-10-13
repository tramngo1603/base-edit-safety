"""Enforce the proposed ABE7.10 domain without granting biological approval."""
from __future__ import annotations

from copy import deepcopy
from functools import lru_cache
import hashlib
from importlib.resources import files
import json

from .common import PopulationError

PROFILE_ID = "abe710-spcas9-hek293t-plasmid-v1"


@lru_cache(maxsize=1)
def _profile_bytes():
    return files(__package__).joinpath("protocols", PROFILE_ID + ".json").read_bytes()


def profile_definition():
    """Return a separate copy of the pinned proposal, never a caller approval."""
    return json.loads(_profile_bytes())


def binding(policy):
    value = policy.get("supported_domain")
    if value is not None and value.get("profile_id") != PROFILE_ID:
        raise PopulationError("Unknown supported-domain profile")
    return value


class DomainAbstention(PopulationError):
    """A structurally valid request falls outside the opt-in proposal."""

    def __init__(self, assessment):
        self.assessment = assessment
        super().__init__("Supported-domain request abstained: " + "; ".join(assessment["request"]["reason_codes"]))


def _editor_reasons(editor, profile):
    return ["domain_unsupported_editor"] if any(editor.get(k) != v for k, v in profile["editor"].items()) else []


def _context_reasons(context, profile):
    expected = profile["context"]
    return ["domain_unsupported_context"] if any(
        (context.get(k, {}).get("value") if k == "cell_state" else context.get(k)) != v
        for k, v in expected.items()) else []


def _endpoint_reasons(endpoint, context_id, profile):
    if (endpoint.get("name") not in profile["endpoint_names"] or endpoint.get("unit") != "fraction"
            or endpoint.get("measurement_level") != "reads" or not endpoint.get("requires_sequence")
            or endpoint.get("allele_resolution") != "required"
            or context_id not in endpoint.get("required_context_ids", [])):
        return ["domain_unsupported_endpoint"]
    return []


def target_reasons(policy, index, *, design_id=None, context_id=None, endpoint_id=None):
    """Reject targets outside the bound editor, material and readout."""
    declared = binding(policy)
    if declared is None:
        return []
    profile, reasons = profile_definition(), []
    if design_id is not None:
        design = index.get(design_id, {})
        reasons += _editor_reasons(index.get(design.get("editor_definition_id"), {}), profile)
    reasons += _context_reasons(index.get(context_id, {}), profile)
    if context_id != declared["context_id"]:
        reasons.append("domain_context_outside_bound_target")
    if endpoint_id is not None:
        reasons += _endpoint_reasons(index.get(endpoint_id, {}), context_id, profile)
        if endpoint_id not in declared["endpoint_ids"]:
            reasons.append("domain_endpoint_outside_bound_target")
    return sorted(set(reasons))


def biological_reasons(policy, index, *, design_id=None, context_id=None, endpoint_id=None):
    """Keep target applicability separate from source outcome and assay quality."""
    if binding(policy) is None:
        return []
    reasons = list(profile_definition()["pending_biological_reasons"])
    reasons += target_reasons(policy, index, design_id=design_id, context_id=context_id, endpoint_id=endpoint_id)
    context = index.get(context_id, {})
    for field in ("exposure", "dose", "measurement_timepoint"):
        if context.get(field, {}).get("value", {}).get("value") is None:
            reasons.append("domain_" + field + "_unresolved")
    return sorted(set(reasons))


def assess_domain(config, search=None):
    """Describe request support, actual interval execution and pending claims."""
    policy = config["policies"]["evidence"]
    declared = binding(policy)
    if declared is None:
        return None
    profile = profile_definition()
    index = {row["record_id"]: row for row in config["catalogs"]}
    search_policy = config["policies"]["search"]
    domain = index[search_policy["domain_id"]]
    computational, target, intended_contigs = [], [], set()
    for identifier in search_policy["profile_ids"]:
        search_profile = index[identifier]
        design = index[search_profile["design_id"]]
        editor = index[design["editor_definition_id"]]
        intended_contigs.update(index[key]["contig"] for key in design["intended_locus_ids"])
        computational += _editor_reasons(editor, profile)
        if len(design["guide_sequence"]) != 20:
            computational.append("domain_unsupported_guide_geometry")
        window = editor["editing_window"]
        if any(window.get(k) != profile["editing_window"][k] for k in ("start_1based", "end_1based")):
            computational.append("domain_unsupported_window_annotation")
        if (search_profile["pam_patterns"] != ["NGG"] or set(search_profile["strands"]) != {"+", "-"}
                or not 0 <= search_profile["max_mismatches"] <= 6
                or search_profile["max_total_edits"] != search_profile["max_mismatches"]
                or search_profile["max_dna_bulges"] != 0 or search_profile["max_rna_bulges"] != 0
                or search_profile["total_edit_accounting"] != "mismatch_plus_bulged_bases"
                or search_profile["retention"] != "all"):
            computational.append("domain_unsupported_search_settings")
        if (search_profile["backend_name"], search_profile["backend_version"]) != ("comparative_local", "v1"):
            computational.append("domain_native_backend_not_integrated")
    if config["mode"] != "synthetic_test":
        computational.append("domain_real_processing_pending")
    if config["policies"]["predictions"]["enabled"] or config["policies"]["comparison"]["include_predictions"]:
        computational.append("domain_predictions_excluded")
    contigs = sorted({row["contig"] for row in domain["intervals"]})
    # Mock autosomal inputs use the fixture naming convention.
    for contig in sorted(set(contigs) | intended_contigs):
        primary = contig.removeprefix("synthetic_") if config["mode"] == "synthetic_test" else contig
        if primary in {"chrX", "chrY"}:
            computational.append("domain_sex_chromosome_reconstruction_unsupported")
        elif primary not in profile["reference"]["primary_contigs"]:
            computational.append("domain_nonprimary_contig_excluded")
    context = index[declared["context_id"]]
    target += _context_reasons(context, profile)
    for endpoint_id in declared["endpoint_ids"]:
        target += _endpoint_reasons(index[endpoint_id], declared["context_id"], profile)
    biological = biological_reasons(policy, index, context_id=declared["context_id"])
    coverage = list(profile["pending_coverage_reasons"])
    scopes = [row for name, rows in (search or {}).items() if name.endswith(".jsonl")
              for row in rows if row["record_type"] == "comparison_search_scope"]
    donors = {row["record_id"] for row in config["catalogs"] if row["record_type"] == "comparison_donor"}
    expected_pairs = {(design, donor) for design in search_policy["design_ids"] for donor in donors}
    actual_pairs = {(row.get("design_id"), row.get("donor_id", {}).get("value")) for row in scopes}
    scope_status = "not_evaluated" if not scopes else "complete" if (
        expected_pairs and actual_pairs == expected_pairs and len(scopes) == len(expected_pairs)
        and all(row["execution_status"] == "complete" for row in scopes)) else "partial"
    callable_scope = config["coverage_basis"] == "callable_supported" and scope_status == "complete" and any(
        r["role"] == "callability" for r in config["input_refs"])
    if not callable_scope:
        coverage.append("domain_callability_not_established")
    if scope_status != "complete":
        coverage.append("domain_declared_scope_" + scope_status)
    request = sorted(set(computational + target))
    return {
        "profile_id": PROFILE_ID, "profile_sha256": hashlib.sha256(_profile_bytes()).hexdigest(),
        "binding": deepcopy(declared), "domain_id": domain["record_id"],
        "request": {"status": "abstain" if request else "supported_for_artificial_checks", "reason_codes": request},
        "computational_support": {"status": "abstain" if computational else "supported_for_artificial_checks",
                                  "reason_codes": sorted(set(computational))},
        "coverage_completeness": {"status": "incomplete", "reason_codes": sorted(set(coverage)),
            "declared_scope_execution": scope_status, "declared_contigs": contigs,
            "callability": "supported_in_declared_intervals" if callable_scope else "not_established",
            "scope_reason_codes": sorted({reason for row in scopes for reason in row.get("reason_codes", [])}),
            "required_primary_contigs": profile["reference"]["primary_contigs"],
            "reference_bundle_id": domain["reference_bundle_id"]},
        "biological_applicability": {"status": "abstain", "reason_codes": sorted(set(biological + target))},
        "release_complete": False,
    }


def require_supported_request(config):
    assessment = assess_domain(config)
    if assessment is not None and assessment["request"]["status"] == "abstain":
        raise DomainAbstention(assessment)
    return assessment
