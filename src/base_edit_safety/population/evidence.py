"""Keep sequence representation, experimental outcomes, and donor credit separate."""
from __future__ import annotations

from .common import PopulationError, fields, identity, nonempty, record, unique
from .identity import definition_hashes
from ..measurement_units import editing_detection_limit
from ..pipeline.candidate_panel import _has_reported_text


def representation_replay(ledger, reference, guides, run_id):
    """Describe guarded representation changes without rewriting the source ledger."""
    from copy import deepcopy
    from ..pipeline.candidate_panel import _alignment_block_for_site
    from .identity import reverse_complement
    guides = unique(ledger.get("guides", guides), "guide_id")
    rows = []
    for original in ledger["candidates"]:
        if original.get("allele_states") or not original.get("nominations") or any(
                n["nomination_kind"] != "experimental" for n in original["nominations"]):
            continue
        site = deepcopy(original)
        verified = True
        for nomination in site["nominations"]:
            try:
                seq = reference.fetch(site["chrom"], site["start_0based"], site["end_0based_exclusive"])
            except PopulationError:
                verified = False
                break
            oriented = seq if nomination.get("strand") == "+" else reverse_complement(seq)
            if nomination.get("site_sequence") != oriented or "N" in seq:
                verified = False
                break
            nomination["reference_validation"] = {"status": "verified", "reference_sha256": reference.sha256}
        block = _alignment_block_for_site(site, guides.get(site["guide_id"])) if verified else None
        rows.append(record("representation_replay", run_id, site_id=site["site_id"],
            before="experimental_candidate_model_unsupported", after="verified_reference_model" if block else "unresolved",
            reference_model=block, source_measurements_unchanged=True, new_donor_measurement_credits=0,
            interpretation="Reference representation only; no sample allele, QC, control or detection limit inferred."))
    return rows


def resolve_measurements(ledger, linkage, ledger_sha256):
    if ledger is None and linkage is None:
        return []
    if ledger.get("record_type") != "candidate_panel_ledger" or ledger.get("schema_version") != "3":
        raise PopulationError("Measurement integration requires a candidate ledger schema 3")
    source_guides = unique(ledger.get("guides", []), "guide_id")
    fields(linkage, {"record_type", "schema_version", "ledger_sha256", "people", "links"}, label="measurement linkage")
    if linkage["record_type"] != "measurement_linkage" or linkage["schema_version"] != "1" or linkage["ledger_sha256"] != ledger_sha256:
        raise PopulationError("Measurement linkage does not bind this ledger")
    people = {}
    for row in linkage["people"]:
        fields(row, {"source_person_id", "canonical_donor_id", "identity_status", "source_assertion"}, label="person linkage")
        nonempty(row["source_person_id"], "source person")
        nonempty(row["source_assertion"], "person linkage source assertion")
        if row["identity_status"] not in {"resolved", "unknown", "ambiguous"}:
            raise PopulationError("Invalid person identity status")
        if row["identity_status"] == "resolved":
            nonempty(row["canonical_donor_id"], "canonical biological donor")
        elif row["canonical_donor_id"] is not None:
            raise PopulationError("Unresolved person identity cannot name a canonical donor")
        if row["source_person_id"] in people:
            raise PopulationError("Duplicate source person linkage")
        people[row["source_person_id"]] = row
    raw = {}
    for site in ledger["candidates"]:
        for m in site["measurements"]:
            key = (site["site_id"], m["source_version_or_citation"], m["assay"], m["source_record_id"])
            if key in raw:
                raise PopulationError("Ambiguous ledger measurement identity")
            raw[key] = (site, m)
    linked = {}
    for link in linkage["links"]:
        fields(link, {"site_id", "source_version_or_citation", "assay", "source_record_id", "source_person_id",
                      "biological_context_id", "context", "configuration_id", "target_path", "locus_resolved", "source_assertion"}, label="measurement link")
        key = tuple(link[k] for k in ("site_id", "source_version_or_citation", "assay", "source_record_id"))
        if key not in raw or key in linked:
            raise PopulationError("Missing or duplicate linked measurement")
        person = people.get(link["source_person_id"])
        if link["source_person_id"] is not None and person is None:
            raise PopulationError("Measurement person does not resolve")
        if type(link["locus_resolved"]) is not bool:
            raise PopulationError("Locus resolution declaration must be boolean")
        nonempty(link["source_assertion"], "measurement linkage source assertion")
        linked[key] = {**link, "donor_id": person["canonical_donor_id"] if person else None,
                       "person_identity_status": person["identity_status"] if person else "unknown"}
    return [{"site": site, "measurement": m, "link": linked.get(key),
             "source_guide": source_guides.get(site["guide_id"])} for key, (site, m) in sorted(raw.items())]


def evidence_for_unit(config, unit, candidate, measurement, link, guide, run_id, *, source_guide=None):
    """Re-evaluate evidence from source fields; usability flags are never inputs."""
    reasons, context_reasons = [], []
    try:
        source_definitions = definition_hashes(source_guide) if source_guide else None
    except KeyError:
        source_definitions = None
    if source_definitions is None:
        reasons.append("source_guide_definition_unknown")
    else:
        for label, source, current in zip(("guide", "editor"), source_definitions, definition_hashes(guide)):
            if source != current:
                reasons.append("source_" + label + "_definition_mismatch")
    for key, expected in (("editor", guide["editor"]), ("guide_id", guide["guide_id"]),
                          ("chrom", candidate["contig"]), ("strand", candidate["strand"]),
                          ("start_0based", candidate["bounding_interval"][0]),
                          ("end_0based_exclusive", candidate["bounding_interval"][1]),
                          ("genome_build", config["reference_policy"]["genome_build"]),
                          ("allele_sequence", candidate["sequence_genomic_forward"])):
        value = measurement.get(key)
        if value is None:
            reasons.append(key + "_unknown")
        elif value != expected:
            reasons.append(key + "_mismatch")
    if measurement.get("allele_resolution") != "resolved":
        reasons.append("allele_resolution_unresolved")
    if not _has_reported_text(measurement.get("allele_evidence")):
        reasons.append("allele_evidence_missing")
    if link is None:
        reasons.append("measurement_linkage_missing")
        context_reasons.append("context_unresolved")
    else:
        if link["configuration_id"] != candidate["configuration_id"] or link["target_path"] != candidate["target_path"]:
            reasons.append("configuration_mismatch")
        if not link["locus_resolved"]:
            reasons.append("locus_resolution_unresolved")
        if link["person_identity_status"] != "resolved":
            reasons.append("biological_donor_identity_unresolved")
        if (link["biological_context_id"] != config["assay_policy"]["biological_context_id"]
                or link["context"] != config["assay_policy"]["context"]):
            context_reasons.append("biological_context_mismatch")
    quality = []
    if measurement.get("qc_status") != "passed":
        quality.append("qc_not_passed")
    for key in ("treated_group", "control_group", "biological_context", "assay"):
        if not _has_reported_text(measurement.get(key)):
            quality.append(key + "_missing")
    if measurement.get("measurement_kind") != unit["measurement_class"]:
        quality.append("measurement_class_mismatch")
    if measurement.get("status") not in {"edited", "not_detected"}:
        quality.append("outcome_not_usable")
    detection_limit = None
    if measurement.get("status") == "not_detected":
        detection_limit = editing_detection_limit(measurement.get("limit_of_detection"),
            config["assay_policy"].get("max_detection_limit_fraction"))
        if detection_limit["reason"]:
            quality.append(detection_limit["reason"])
    all_reasons = reasons + context_reasons + quality
    payload = {"unit_id": unit["id"], "source_measurement": measurement, "link": link, "source_guide": source_guide}
    return record("evidence_link", run_id, evidence_id=identity("evidence", payload), unit_id=unit["id"],
                  donor_id=link["donor_id"] if link else None, source_measurement=measurement,
                  source_linkage=link, source_guide=source_guide, outcome=measurement.get("status"),
                  identity_match="mismatch" if any(r.endswith("mismatch") for r in reasons) else "unknown" if reasons else "matched",
                  context_match="matched" if not context_reasons else "unknown" if link is None else "mismatch",
                  usable=not all_reasons, reasons=all_reasons, detection_limit_assessment=detection_limit)


def reconcile(config, units, candidates, guides, measurements, run_id):
    candidates = unique(candidates, "configuration_id")
    guides = unique(guides, "guide_id")
    result = []
    for unit in units:
        c = candidates[unit["configuration_id"]]
        for row in measurements:
            site = row["site"]
            if site["guide_id"] != c["guide_id"] or site["chrom"] != c["contig"] \
                    or [site["start_0based"], site["end_0based_exclusive"]] != c["bounding_interval"]:
                continue
            result.append(evidence_for_unit(config, unit, c, row["measurement"], row["link"], guides[c["guide_id"]], run_id,
                                            source_guide=row["source_guide"]))
    return sorted(result, key=lambda e: e["evidence_id"])


def unmatched_measurements(units, candidates, guides, measurements, links, run_id):
    """Explain retained source rows that have no current coverage-unit comparison."""
    matched = {identity("measurement", e["source_measurement"]) for e in links}
    current_guides = unique(guides, "guide_id")
    unit_configurations = {u["configuration_id"] for u in units}
    rows = []
    for row in measurements:
        if identity("measurement", row["measurement"]) in matched:
            continue
        site, source_guide = row["site"], row["source_guide"]
        reasons = []
        if site["guide_id"] not in current_guides:
            reasons.append("source_guide_id_not_in_study")
        if not any(c["configuration_id"] in unit_configurations and c["guide_id"] == site["guide_id"]
                   and c["contig"] == site["chrom"]
                   and c["bounding_interval"] == [site["start_0based"], site["end_0based_exclusive"]] for c in candidates):
            reasons.append("no_coverage_unit_for_source_guide_and_locus")
        same_definitions = sorted(gid for gid, guide in current_guides.items()
                                  if source_guide and definition_hashes(source_guide) == definition_hashes(guide))
        rows.append({"site_id": site["site_id"], "source_measurement": row["measurement"],
                     "source_guide": source_guide, "source_linkage": row["link"], "reasons": reasons,
                     "current_guide_ids_with_same_definitions": same_definitions})
    return record("unmatched_measurements", run_id, rows=rows,
        interpretation="Source labels are preserved. Matching definitions do not establish a guide alias or supply measurement credit.")


def evidence_summary(links, unmatched=()):
    """Count observations and conflicting donor-unit pairs without averaging outcomes."""
    from collections import Counter
    outcomes, reasons, pairs = Counter(), Counter(), {}
    for row in links:
        if row["usable"]:
            outcomes[row["outcome"]] += 1
            pairs.setdefault((row["donor_id"], row["unit_id"]), set()).add(row["outcome"])
        else:
            reasons.update(row["reasons"])
    unmatched_reasons = Counter(reason for row in unmatched for reason in row["reasons"])
    return {"usable_positive_links": outcomes["edited"], "usable_nondetection_links": outcomes["not_detected"],
            "conflicting_donor_unit_pairs": sum({"edited", "not_detected"} <= values for values in pairs.values()),
            "excluded_links": sum(not e["usable"] for e in links), "excluded_link_reasons": dict(sorted(reasons.items())),
            "unmatched_source_measurements": len(unmatched), "unmatched_source_reasons": dict(sorted(unmatched_reasons.items()))}
