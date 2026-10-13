"""Explicitly migrate legacy observation history without rewriting its sources."""
from copy import deepcopy

from .common import PopulationError, digest, fields
from .comparison_contracts import (
    _walk, make_record, record_schema, validate_catalog,
)
from .observations import intrinsic_assessment

MIGRATION_VERSION = "observation-outcomes-v1-to-v2"


def _references(row):
    references = set()
    def collect(value, schema, _location):
        if "x-record-types" in schema:
            references.add(value)
    _walk(row, record_schema(row["record_type"], row["schema_version"]), collect)
    return references


def migrate_observation_catalog(records, adjudications=None):
    """Return a new catalog, retained states and an explicit invalidation receipt.

    Adjudications are keyed by the exact active v1 revision ID. Each value has
    adjudicated_outcome and adjudication_reason. Missing adjudication withholds
    conclusive credit. Dependent scientific products must be recomputed.
    """
    index = validate_catalog(records)
    adjudications = {} if adjudications is None else adjudications
    if not isinstance(adjudications, dict):
        raise PopulationError("Migration adjudications must be keyed by legacy revision ID")
    states = [row for row in records if row["record_type"] == "comparison_observation_state"
              and index[row["active_revision_id"]]["schema_version"] == "1"]
    active = {row["active_revision_id"] for row in states}
    if set(adjudications) - active:
        raise PopulationError("Adjudication must name an active version-1 observation revision")
    additions, superseded, mapping = [], [], {}
    for state in sorted(states, key=lambda row: row["source_record_id"]):
        old = index[state["active_revision_id"]]
        choice = adjudications.get(old["record_id"], {
            "adjudicated_outcome": "inconclusive",
            "adjudication_reason": "Legacy revision has no explicit reviewed adjudication; conclusive credit withheld",
        })
        fields(choice, {"adjudicated_outcome", "adjudication_reason"}, label="migration adjudication")
        values = deepcopy(old)
        for key in ("record_type", "schema_version", "record_id", "observation_id", "applicability"):
            values.pop(key)
        values.update(choice, normalization_schema_version="2")
        values["applicability"] = intrinsic_assessment(values, index)
        new = make_record("observation_revision", version="2", **values)
        mapping[old["record_id"]] = new["record_id"]
        additions.append(new)
        additions.append(make_record("observation_state", source_record_id=state["source_record_id"],
            active_revision_id=new["record_id"], revision_ids=sorted(set(state["revision_ids"] + [new["record_id"]])),
            supersession_reason="Explicit outcome migration; original revisions retained without runtime credit"))
        superseded.append(deepcopy(state))
    old_states = {row["record_id"] for row in superseded}
    changed = set(mapping) | old_states
    invalidated = set()
    retained_history = {"comparison_source_observation", "comparison_observation_revision", "comparison_observation_state"}
    while True:
        found = {row["record_id"] for row in records if row["record_type"] not in retained_history
                 and _references(row) & changed} - invalidated
        if not found:
            break
        invalidated.update(found)
        changed.update(found)
    result = [deepcopy(row) for row in records if row["record_id"] not in old_states | invalidated] + additions
    # Refuse a migration that cannot retain a historical revision's dependencies.
    validate_catalog(result)
    result.sort(key=lambda row: (row["record_type"], row["record_id"]))
    return {"records": result, "superseded_states": superseded,
            "receipt": {"migration_version": MIGRATION_VERSION,
                "input_catalog_sha256": digest(records), "output_catalog_sha256": digest(result),
                "adjudications_sha256": digest(adjudications), "revision_id_map": mapping,
                "invalidated_record_ids": sorted(invalidated),
                "rebuild_stages": ["observation_normalization", "evidence_applicability", "comparison", "questions", "planning", "report"] if mapping else []}}


def main(argv=None):
    """Write a fresh protected migration bundle from an explicit JSON catalog."""
    import argparse
    from .common import checksum, json_text, jsonl_text, read_json, write_bundle
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", required=True, help="Complete versioned record array in JSON")
    parser.add_argument("--adjudications", help="Reviewed outcomes keyed by active legacy revision ID")
    parser.add_argument("--output", required=True, help="New protected output directory")
    args = parser.parse_args(argv)
    inputs = [args.catalog] + ([args.adjudications] if args.adjudications else [])
    before = {path: checksum(path) for path in inputs}
    result = migrate_observation_catalog(read_json(args.catalog),
        read_json(args.adjudications) if args.adjudications else None)
    if before != {path: checksum(path) for path in inputs}:
        raise PopulationError("Migration inputs changed during processing")
    result["receipt"]["input_file_sha256"] = before
    observations = [row for row in result["records"] if row["record_type"] in {
        "comparison_source_observation", "comparison_observation_revision", "comparison_observation_state",
        "comparison_observation_equivalence"}]
    write_bundle(args.output, {"catalog.json": json_text(result["records"]),
        "observations.jsonl": jsonl_text(observations),
        "superseded-states.json": json_text(result["superseded_states"]),
        "migration-receipt.json": json_text(result["receipt"])}, inputs=inputs)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
