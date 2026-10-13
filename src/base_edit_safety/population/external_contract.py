"""Validate external-search source details, retained rows and execution status."""
from __future__ import annotations

from .common import PopulationError, fields, nonempty


def validate_external_contract(contract):
    """Validate a version-1 contract and return execution and configuration-coverage status."""
    fields(contract, {"record_type", "schema_version", "backend_id", "executable_sha256", "source_revision",
                      "native_schema", "argv", "assembly", "reference_sha256", "phase_handling",
                      "all_configuration_retention", "row_retention", "verified_fixture_ids", "stages"},
           label="external search contract")
    if contract["record_type"] != "external_search_contract" or contract["schema_version"] != "1":
        raise PopulationError("Unsupported external search contract")
    for key in ("backend_id", "executable_sha256", "native_schema", "assembly", "reference_sha256", "phase_handling"):
        nonempty(contract[key], key)
    import re
    for key in ("executable_sha256", "reference_sha256"):
        if not re.fullmatch("[0-9a-f]{64}", contract[key]):
            raise PopulationError("External executable and reference identities require full SHA-256")
    if not isinstance(contract["argv"], list) or not contract["argv"] or any(not isinstance(a, str) or "\x00" in a for a in contract["argv"]):
        raise PopulationError("External commands require explicit argv arrays")
    if contract["row_retention"] not in {"all", "winner_only", "partial", "unknown"} \
            or type(contract["all_configuration_retention"]) is not bool:
        raise PopulationError("Invalid external retention declaration")
    if not isinstance(contract["verified_fixture_ids"], list) or not isinstance(contract["stages"], list):
        raise PopulationError("External verification and execution records are required")
    for stage in contract["stages"]:
        fields(stage, {"name", "status", "exit_code"}, label="external stage")
        if stage["status"] not in {"completed", "failed", "not_run"}:
            raise PopulationError("Unknown external stage status")
        if stage["status"] == "completed" and stage["exit_code"] != 0:
            raise PopulationError("Successful external stage requires exit code zero")
    complete = (contract["all_configuration_retention"] and contract["row_retention"] == "all"
                and bool(contract["verified_fixture_ids"]) and bool(contract["stages"])
                and all(s["status"] == "completed" for s in contract["stages"]))
    return {"execution_status": "failed" if any(s["status"] == "failed" for s in contract["stages"]) else
            "completed" if contract["stages"] and all(s["status"] == "completed" for s in contract["stages"]) else "not_established",
            "row_retention": contract["row_retention"],
            "configuration_coverage": "declared_all_fixture_verified" if complete else "unknown",
            "search_completeness": "requires_backend_specific_scope_validation"}
