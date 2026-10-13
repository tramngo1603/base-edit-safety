"""Compare measured evidence with explicit design, context and endpoint targets."""
from __future__ import annotations

from collections import defaultdict
from copy import deepcopy
from decimal import Decimal, localcontext
from itertools import combinations

from .. import annotation
from ..measurement_units import editing_detection_limit
from .common import PopulationError, fields, integer, known, read_json, unknown
from .comparison_contracts import (
    active_observations, make_record, record_registry, validate_catalog, validate_observation_locus,
)
from .comparison_outcomes import outcome_reasons
from .comparison_domain import assess_domain, biological_reasons


def _dimension(status="applicable", *reasons):
    return {"status": status, "reason_codes": sorted(set(reasons))}


def _zero(observation):
    value = observation["normalized_value"]["value"]
    return value is not None and Decimal(str(value)) == 0


def _quality(observation, policy):
    reasons = outcome_reasons(observation)
    if policy["require_controls"] and observation["control_status"] != "passed":
        reasons.append("controls_not_passed")
    if policy["require_qc"] and observation["qc_status"] != "passed":
        reasons.append("qc_not_passed")
    if observation["normalized_value"]["value"] is None:
        reasons.append("measurement_unavailable_or_censored")
    elif observation["normalized_unit"] == "fraction" and not 0 <= Decimal(str(observation["normalized_value"]["value"])) <= 1:
        reasons.append("measured_fraction_outside_bounds")
    if _zero(observation) and policy["require_detection_limit_for_zero"]:
        limit = observation["detection_limit"]
        value = limit["value"]["value"]
        if observation["normalized_unit"] == "fraction":
            result = editing_detection_limit({"value": value, "unit": limit["unit"]})
            if result["status"] != "accepted":
                reasons.append("zero_without_usable_detection_limit")
        elif limit["unit"] != observation["normalized_unit"] or value is None or value <= 0:
            reasons.append("zero_without_usable_detection_limit")
    return _dimension("unresolved", *reasons) if reasons else _dimension()


def _related_design(left, right, index):
    if left["guide_sequence"] == right["guide_sequence"]:
        return True
    left_genes = {index[key]["gene"] for key in left["intended_locus_ids"]}
    right_genes = {index[key]["gene"] for key in right["intended_locus_ids"]}
    return bool(left_genes & right_genes)


def assess_observation(observation, index, *, design_id, context_id, endpoint_id,
                       donor_id=None, sequence_configuration_id=None,
                       mode, coverage_basis, policy):
    """Assess one measured observation; predictions cannot enter this pathway."""
    if observation.get("record_type") != "comparison_observation_revision":
        raise PopulationError("Evidence reconciliation accepts measured observation revisions only")
    validate_observation_locus(observation, index)
    for identifier, kind in ((design_id, "design"), (context_id, "experimental_context"),
                              (endpoint_id, "endpoint_definition")):
        if index.get(identifier, {}).get("record_type") != "comparison_" + kind:
            raise PopulationError("Evidence target has an unresolved typed reference")
    identity = _dimension()
    source_design_id = observation["design_id"]["value"]
    if source_design_id is None:
        identity = _dimension("unresolved", "source_design_unresolved")
    elif source_design_id != design_id:
        related = _related_design(index[source_design_id], index[design_id], index)
        identity = _dimension("related" if related else "incompatible", "different_complete_design")
    endpoint = index[endpoint_id]
    if endpoint["requires_sequence"]:
        source_sequence = observation["sequence_configuration_id"]["value"]
        if source_sequence is None:
            identity["reason_codes"].append("source_sequence_unresolved")
            if identity["status"] == "applicable":
                identity["status"] = "unresolved"
        elif sequence_configuration_id is not None and source_sequence != sequence_configuration_id:
            identity = _dimension("incompatible", *identity["reason_codes"], "different_sequence_configuration")
        if endpoint["locus_resolution"] in {"copy_specific", "locus_specific"} and len(observation["locus_ids"]) != 1:
            identity["reason_codes"].append("required_locus_resolution_unavailable")
            if identity["status"] == "applicable":
                identity["status"] = "unresolved"
    subject = observation["subject"]
    resolved_donor = subject["resolution_status"] == "resolved" and subject["level"] in {"donor", "specimen"}
    if donor_id is not None and (not resolved_donor or subject["donor_id"]["value"] != donor_id):
        identity["reason_codes"].append("target_donor_linkage_unavailable")
        if identity["status"] == "applicable":
            identity["status"] = "unresolved"
    source_context = observation["context_id"]["value"]
    context = _dimension()
    if source_context is None:
        context = _dimension("unresolved", "source_context_unresolved")
    elif source_context != context_id:
        source = index[source_context]
        target = index[context_id]
        compatible_species = source["organism"] == target["organism"]
        context = _dimension("related" if compatible_species and policy["allow_related_context"] else "incompatible",
                             "different_experimental_context")
    endpoint_match = _dimension() if observation["endpoint_id"] == endpoint_id else \
        _dimension("incompatible", "different_endpoint_definition")
    if observation["endpoint_id"] == endpoint_id and endpoint["required_context_ids"] \
            and context_id not in endpoint["required_context_ids"]:
        contextual_use = policy["allow_related_context"] and source_context is not None \
            and index[source_context]["organism"] == index[context_id]["organism"]
        endpoint_match = _dimension("related" if contextual_use else "incompatible", "endpoint_context_not_supported")
    domain_reasons = biological_reasons(policy, index, design_id=design_id,
        context_id=context_id, endpoint_id=endpoint_id)
    if domain_reasons:
        domain_reasons += biological_reasons(policy, index, design_id=source_design_id,
            context_id=source_context, endpoint_id=observation["endpoint_id"])
        context["reason_codes"] = sorted(set(context["reason_codes"] + domain_reasons))
        if context["status"] == "applicable":
            context["status"] = "unresolved"
    quality = _quality(observation, policy)
    dimensions = (identity, context, endpoint_match, quality)
    statuses = {item["status"] for item in dimensions}
    source_has_donor = resolved_donor and subject["donor_id"]["value"] is not None
    if "incompatible" in statuses:
        classification = "incompatible"
    elif "related" in statuses:
        classification = "related_design_context"
    elif "unresolved" in {identity["status"], context["status"], endpoint_match["status"]}:
        classification = "unresolved"
    else:
        classification = "exact_donor_configuration" if source_has_donor and donor_id is not None else "exact_sequence_context"
    reasons = {reason for part in dimensions for reason in part["reason_codes"]}
    credit = classification == "exact_donor_configuration" and statuses == {"applicable"} \
        and policy["allow_empirical_donor_credit"] and mode == "research_sample_planning" \
        and coverage_basis != "retained_callset_model" and observation["applicability"]["empirical_donor_credit"]
    if endpoint["requires_sequence"] and sequence_configuration_id is None:
        credit = False
        reasons.add("exact_target_sequence_not_requested")
    if coverage_basis == "retained_callset_model":
        reasons.add("retained_callset_simulation_no_empirical_donor_credit")
    if not source_has_donor:
        reasons.add("no_resolved_individual_donor_linkage")
    if not credit:
        reasons.add("no_target_donor_assay_credit")
    assessment = {"identity": identity, "context": context, "endpoint": endpoint_match,
        "quality": quality, "classification": classification, "empirical_donor_credit": credit,
        "reason_codes": sorted(reasons)}
    return make_record("evidence_applicability", observation_id=observation["record_id"],
        design_id=design_id, context_id=context_id, endpoint_id=endpoint_id,
        donor_id=known(donor_id) if donor_id is not None else unknown("No individual target donor requested"),
        sequence_configuration_id=known(sequence_configuration_id) if sequence_configuration_id is not None else
            observation["sequence_configuration_id"], assessment=assessment,
        policy_version=policy["policy_version"], coverage_basis=coverage_basis)


def source_summaries(observations, index):
    """Count biological units within each source regimen; retain conflicting values."""
    groups = defaultdict(list)
    for observation in observations:
        source = index[observation["source_record_id"]]
        key = (source["source_study"], observation["design_id"]["value"],
               observation["context_id"]["value"], observation["endpoint_id"],
               source["source_reported_label"]["value"])
        groups[key].append(observation)
    result = []
    for key, rows in sorted(groups.items(), key=lambda item: str(item[0])):
        units = {row["independent_unit_id"]["value"] for row in rows if row["independent_unit_id"]["value"] is not None}
        donors = {row["subject"]["source_subject_id"]["value"] for row in rows
                  if row["subject"]["level"] in {"donor", "specimen"}
                  and row["subject"]["source_subject_id"]["value"] is not None}
        values = defaultdict(set)
        for row in rows:
            group = (row["independent_unit_id"]["value"], row["pairing_key"]["value"])
            if group[0] is not None and row["normalized_value"]["value"] is not None:
                values[group].add(str(row["normalized_value"]["value"]))
        result.append({"source_study": key[0], "design_id": key[1], "context_id": key[2],
            "endpoint_id": key[3], "source_label": key[4], "measurements": len(rows),
            "reported_independent_units": len(units), "source_donor_labels": len(donors),
            "measurements_with_unknown_unit": sum(r["independent_unit_id"]["value"] is None for r in rows),
            "distinct_values_within_unit": [sorted(v) for _, v in sorted(values.items(), key=lambda item: str(item[0])) if len(v) > 1],
            "observation_ids": sorted(r["record_id"] for r in rows),
            "interpretation": "Separate source measurements; no pooled effect estimate or inferred donor relationship"})
    return result


def paired_descriptions(observations, index, design_ids, policy):
    """Describe only unambiguous, explicitly paired measurements in one context."""
    grouped = defaultdict(lambda: defaultdict(list))
    for observation in observations:
        design = observation["design_id"]["value"]
        if observation.get("adjudicated_outcome") != "observed_positive" \
                or design not in design_ids or observation["pairing_key"]["value"] is None \
                or observation["independent_unit_id"]["value"] is None or _quality(observation, policy)["status"] != "applicable":
            continue
        endpoint = index[observation["endpoint_id"]]
        if endpoint["requires_sequence"] and observation["sequence_configuration_id"]["value"] is None:
            continue
        source = index[observation["source_record_id"]]
        key = (source["source_study"], source["source_experiment"], source["source_replicate"],
               observation["context_id"]["value"], observation["endpoint_id"],
               observation["independent_unit_id"]["value"], observation["pairing_key"]["value"],
               tuple(sorted(observation["locus_ids"])) if endpoint["requires_sequence"] else ())
        if key[3] is not None:
            grouped[key][design].append(observation)
    result = []
    for key, arms in sorted(grouped.items(), key=lambda item: str(item[0])):
        for left, right in combinations(sorted(arms), 2):
            if len(arms[left]) != 1 or len(arms[right]) != 1:
                continue
            a, b = arms[left][0], arms[right][0]
            left_value = Decimal(str(a["normalized_value"]["value"]))
            right_value = Decimal(str(b["normalized_value"]["value"]))
            with localcontext() as arithmetic:
                arithmetic.prec = max(28, len(left_value.as_tuple().digits), len(right_value.as_tuple().digits)) \
                    + abs(left_value.as_tuple().exponent - right_value.as_tuple().exponent) + 1
                difference = left_value - right_value
            result.append({"design_ids": [left, right], "context_id": key[3], "endpoint_id": key[4],
                "source_study": key[0], "source_experiment": key[1], "source_replicate": key[2],
                "independent_unit_id": key[5], "pairing_key": key[6],
                "observation_ids": [a["record_id"], b["record_id"]], "left_minus_right": format(difference, "f"),
                "unit": a["normalized_unit"],
                "interpretation": "Descriptive paired source values; no inferential test, pooled estimate or preferred design"})
    return result


def compare_with_evidence(study, search_artifacts, evidence_artifacts, observation_artifacts):
    """Return sequence summaries and design comparisons with matching measured evidence."""
    records = {r["record_id"]: r for r in study["config"]["catalogs"]}
    for artifacts in (search_artifacts, evidence_artifacts, observation_artifacts):
        for name, rows in artifacts.items():
            if name.endswith(".jsonl"):
                for row in rows:
                    old = records.get(row["record_id"])
                    if old is not None and old != row:
                        raise PopulationError("Comparison ancestors contain conflicting record content")
                    records[row["record_id"]] = row
    index = validate_catalog(list(records.values()))
    observations = active_observations(list(index.values()))
    config = study["config"]
    design_ids = set(config["policies"]["search"]["design_ids"])
    active_ids = {row["record_id"] for row in observations}
    applications = [r for r in index.values() if r["record_type"] == "comparison_evidence_applicability"
                    and r["observation_id"] in active_ids]
    by_design = {design: {r["observation_id"] for r in applications if r["design_id"] == design
                 and r["assessment"]["classification"] != "incompatible"} for design in design_ids}
    output, sequence_sets, search_complete = [], defaultdict(set), {}
    for hit in search_artifacts.get("search-audit.json", {}).get("hits", []):
        configuration = index[hit["design_configuration_id"]]
        if configuration["role"] != "intended":
            sequence_sets[configuration["design_id"]].add(configuration["sequence_configuration_id"])
    def revise(row, **changes):
        payload = deepcopy(row)
        kind = payload.pop("record_type")
        for key in ("schema_version", "record_id", record_registry()[kind]):
            payload.pop(key)
        payload.update(changes)
        return make_record(kind, **payload)
    for row in index.values():
        if row["record_type"] == "comparison_design_summary":
            search_complete[row["design_id"]] = row["complete_alignment_count"]["value"] is not None \
                and row["intended_unresolved_donors"] == 0
            output.append(revise(row, evidence_observation_ids=sorted(by_design[row["design_id"]]),
                reason_codes=sorted(set(row["reason_codes"]) - {"empirical_effectiveness_not_assessed"}
                    | {"evidence_context_and_quality_retained_separately", "nomination_counts_do_not_establish_a_preferred_design"})))
        elif row["record_type"] == "comparison_donor_design_summary":
            relevant = by_design[row["design_id"]]
            exact = {r["observation_id"] for r in applications if r["design_id"] == row["design_id"]
                     and r["donor_id"]["value"] == row["donor_id"] and r["assessment"]["empirical_donor_credit"]
                     and _quality(index[r["observation_id"]], config["policies"]["evidence"])["status"] == "applicable"
                     and not biological_reasons(config["policies"]["evidence"], index,
                         design_id=r["design_id"], context_id=r["context_id"], endpoint_id=r["endpoint_id"])}
            output.append(revise(row, evidence_observation_ids=sorted(relevant),
                evidence_availability="exact_donor_evidence" if exact else "contextual_only" if relevant else "none_in_current_sources"))
    pairs = paired_descriptions(observations, index, design_ids, config["policies"]["evidence"])
    domain = assess_domain(config, search_artifacts)
    domain_reasons = domain["biological_applicability"]["reason_codes"] if domain else []
    if domain:
        for pair in pairs:
            pair["biological_applicability"] = deepcopy(domain["biological_applicability"])
    for left, right in combinations(sorted(design_ids), 2):
        matching = [p for p in pairs if p["design_ids"] == [left, right]]
        measured = {identifier for pair in matching for identifier in pair["observation_ids"]}
        related = by_design[left] | by_design[right]
        sequence_difference = sequence_sets[left] != sequence_sets[right]
        if domain_reasons:
            state = "sequence_based_hypothesis" if sequence_difference else "insufficient_comparable_evidence"
            reason = "domain_biological_interpretation_withheld"
        elif matching and any(Decimal(pair["left_minus_right"]) != 0 for pair in matching):
            state, reason = "observed_difference_in_matched_context", "descriptive_difference_in_explicitly_paired_source_values"
        elif matching and not sequence_difference and search_complete.get(left) and search_complete.get(right):
            state, reason = "no_difference_detected_within_scope", "equal_paired_values_and_candidate_sets_within_declared_scope"
        elif sequence_difference and (matching or related):
            state, reason = "tradeoff_requires_experiment", "sequence_and_source_context_require_matched_experimental_comparison"
        elif sequence_difference:
            state, reason = "sequence_based_hypothesis", "candidate_sets_differ_under_the_exact_design_profiles"
        else:
            state, reason = "insufficient_comparable_evidence", "no_comparable_paired_source_measurements"
        output.append(make_record("decision", design_ids=[left, right], domain_id=config["policies"]["comparison"]["domain_id"],
            state=state, observation_ids=sorted(measured or related), prediction_ids=[],
            reason_codes=sorted(set([reason, "no_universal_safety_or_efficacy_ranking", *domain_reasons])),
            question="Which design should be compared next in controlled, context-matched donor experiments?"))
    return {"comparison-summaries.jsonl": sorted(output, key=lambda r: (r["record_type"], r["record_id"])),
            "paired-source-descriptions.json": pairs, **annotate_sequences(study, search_artifacts)}


def reconcile_evidence(study, search_artifacts, observation_artifacts, *, targets=None):
    """Retain every observation and assess it without filtering on search hits."""
    from .comparison_provenance import verify_sources
    verify_sources(study)
    records = list(study["config"]["catalogs"])
    for artifacts in (search_artifacts, observation_artifacts):
        records.extend(row for name, rows in artifacts.items() if name.endswith(".jsonl") for row in rows)
    index = validate_catalog(records)
    observations = active_observations(records)
    config = study["config"]
    applications, unresolved = {}, []
    for observation in observations:
        if targets is None:
            source_context = observation["context_id"]["value"]
            if source_context is None:
                unresolved.append({"observation_id": observation["record_id"], "reason": "source_context_unresolved"})
                continue
            selected = [{"design_id": design_id, "context_id": source_context,
                         "endpoint_id": observation["endpoint_id"]}
                        for design_id in config["policies"]["search"]["design_ids"]]
        else:
            selected = targets
        for target in selected:
            row = assess_observation(observation, index, **target, mode=config["mode"],
                coverage_basis=config["coverage_basis"], policy=config["policies"]["evidence"])
            applications[row["record_id"]] = row
    return {"evidence-applicability.jsonl": sorted(applications.values(), key=lambda r: r["record_id"]),
            "source-evidence-summary.json": source_summaries(observations, index),
            "unresolved-observations.json": unresolved,
            "published-alignment-replay.json": replay_published_alignments(study, search_artifacts),
            "evidence-scope.json": {"active_measurements": len(observations),
                "join_depends_on_nomination_admission": False,
                "default_context": "Each source's actual reported context; target experiments require their own assessment",
                "empirical_target_donor_credits": sum(row["assessment"]["empirical_donor_credit"] for row in applications.values())}}


def replay_published_alignments(study, search_artifacts):
    """Account for source alignments even when the active profile excludes them."""
    from .identity import matches_pam
    from .contracts import validate_schema
    from .comparison_contracts import record_schema
    catalog = {row["record_id"]: row for row in study["config"]["catalogs"]}
    records = {row["record_id"]: row for name, rows in search_artifacts.items() if name.endswith(".jsonl") for row in rows}
    index = {**catalog, **records}
    profiles = {catalog[key]["design_id"]: catalog[key] for key in study["config"]["policies"]["search"]["profile_ids"]}
    domain = catalog[study["config"]["policies"]["search"]["domain_id"]]
    candidates = defaultdict(list)
    admitted = {hit["alignment_id"] for hit in search_artifacts.get("search-audit.json", {}).get("hits", [])}
    for record in records.values():
        if record["record_type"] != "comparison_alignment" or record["record_id"] not in admitted:
            continue
        configuration = index[record["design_configuration_id"]]
        sequence = index[configuration["sequence_configuration_id"]]
        guide = "".join(column["guide"] or "-" for column in record["columns"])
        target = "".join(column["target"] or "-" for column in record["columns"])
        candidates[configuration["design_id"]].append((record, sequence, guide, target))
    source_hashes = {row["sha256"] for row in study["sources"].values()}
    output, keys = [], set()
    for source in study["sources"].values():
        if source["role"] != "published_alignments":
            continue
        mapping = read_json(source["path"])
        fields(mapping, {"schema_version", "mapping_version", "rows"}, label="published alignment mapping")
        if mapping["schema_version"] != "1" or not isinstance(mapping["mapping_version"], str) or not mapping["mapping_version"].strip() or not isinstance(mapping["rows"], list):
            raise PopulationError("Invalid published alignment mapping")
        for row in mapping["rows"]:
            fields(row, {"source_key", "design_id", "source_artifact_id", "source_locator", "source_positive",
                "source_reported_label", "alignment", "unresolved_reason"}, label="published alignment row")
            if not isinstance(row["source_key"], str) or not row["source_key"].strip() or row["source_key"] in keys:
                raise PopulationError("Published alignment source key is empty or repeated")
            keys.add(row["source_key"])
            if type(row["source_positive"]) is not bool:
                raise PopulationError("Published positivity must preserve an explicit source assertion")
            if not isinstance(row["source_reported_label"], str) or not row["source_reported_label"].strip():
                raise PopulationError("Published outcome must retain a nonempty source label")
            source_artifact = index.get(row["source_artifact_id"], {})
            if source_artifact.get("record_type") != "comparison_source_artifact" or source_artifact.get("sha256") not in source_hashes:
                raise PopulationError("Published alignment lacks its pinned source bytes")
            validate_schema(row["source_locator"], record_schema("comparison_source_observation")["properties"]["locator"])
            design_id = row["design_id"]
            if design_id is not None and (not isinstance(design_id, str) or index.get(design_id, {}).get("record_type") != "comparison_design"):
                raise PopulationError("Published alignment design is unresolved")
            if design_id is None and (not isinstance(row["unresolved_reason"], str) or not row["unresolved_reason"].strip()):
                raise PopulationError("Unresolved source design requires its mapping reason")
            expected = row["alignment"]
            result = {**row, "mapping_sha256": source["sha256"], "exact_alignment_ids": [],
                      "related_alignment_ids": [], "reason_codes": [], "empirical_target_donor_credit": False}
            if expected is None:
                if not isinstance(row["unresolved_reason"], str) or not row["unresolved_reason"].strip():
                    raise PopulationError("Unresolved source alignment requires a reason")
                result.update(disposition="unresolved_source_mapping", reason_codes=[row["unresolved_reason"]])
                output.append(result)
                continue
            if design_id is not None and row["unresolved_reason"] is not None:
                raise PopulationError("Resolved source alignment cannot retain an unresolved marker")
            fields(expected, {"contig", "strand", "start", "end", "guide", "target", "pam", "path"}, label="published alignment")
            integer(expected["start"], "published alignment start")
            integer(expected["end"], "published alignment end", minimum=expected["start"] + 1)
            if expected["strand"] not in {"+", "-"} or any(not isinstance(expected[k], str) or not expected[k] for k in ("contig", "guide", "target", "pam")):
                raise PopulationError("Invalid source alignment geometry")
            guide, target, pam = expected["guide"], expected["target"], expected["pam"]
            if len(guide) != len(target) or set(guide + target) - set("ACGT-") or set(pam) - set("ACGT") or any(a == b == "-" for a, b in zip(guide, target)):
                result.update(disposition="unsupported_representation", reason_codes=["unsupported_source_alignment_columns"])
                output.append(result)
                continue
            if design_id is not None and guide.replace("-", "") != index[design_id]["guide_sequence"]:
                raise PopulationError("Published alignment guide differs from its complete design")
            sequence = target.replace("-", "") + pam
            path = expected["path"]
            expected_positions = list(range(expected["start"], expected["end"]))
            if expected["strand"] == "-":
                expected_positions.reverse()
            if not isinstance(path, list) or len(path) != len(sequence) or len(path) != len(expected_positions):
                raise PopulationError("Published sequence/path length is inconsistent")
            for point, base, position in zip(path, sequence, expected_positions):
                fields(point, {"kind", "pos0", "base"}, label="published reference path")
                if point != {"kind": "reference", "pos0": position, "base": base}:
                    raise PopulationError("Published reference path differs from its mapped sequence")
            if design_id is None:
                result.update(disposition="unresolved_source_mapping", reason_codes=["source_design_unresolved", row["unresolved_reason"]])
                output.append(result)
                continue
            for alignment, candidate, actual_guide, actual_target in candidates[row["design_id"]]:
                if candidate["contig"] != expected["contig"] or candidate["strand"] != expected["strand"]:
                    continue
                positions = [p["pos0"] for p in candidate["path"] if p["kind"] == "reference"]
                same = candidate["path"] == path and actual_guide == guide and actual_target == target and alignment["pam_sequence"] == pam
                if same:
                    result["exact_alignment_ids"].append(alignment["record_id"])
                elif positions and min(positions) < expected["end"] and max(positions) >= expected["start"]:
                    result["related_alignment_ids"].append(alignment["record_id"])
            profile = profiles.get(row["design_id"])
            exclusions = []
            if profile is None:
                exclusions.append("source_design_outside_active_comparison")
            else:
                mismatches = sum(a != b and a != "-" and b != "-" for a, b in zip(guide, target))
                dna, rna = guide.count("-"), target.count("-")
                if mismatches > profile["max_mismatches"]:
                    exclusions.append("source_mismatches_exceed_profile")
                if dna > profile["max_dna_bulges"] or rna > profile["max_rna_bulges"]:
                    exclusions.append("source_bulges_exceed_profile")
                if mismatches + dna + rna > profile["max_total_edits"]:
                    exclusions.append("source_total_edits_exceed_profile")
                if not matches_pam(pam, profile["pam_patterns"]):
                    exclusions.append("source_pam_outside_profile")
                if expected["strand"] not in profile["strands"]:
                    exclusions.append("source_strand_outside_profile")
            if not any(interval["contig"] == expected["contig"] and interval["start"] <= expected["start"] and expected["end"] <= interval["end"] for interval in domain["intervals"]):
                exclusions.append("source_alignment_outside_domain")
            if result["exact_alignment_ids"]:
                result["disposition"] = "recovered_exact_alignment"
            elif exclusions:
                result.update(disposition="outside_declared_search_policy", reason_codes=exclusions)
            elif result["related_alignment_ids"]:
                result.update(disposition="related_alternative_alignment", reason_codes=["reported_alignment_not_recovered_exactly"])
            else:
                result.update(disposition="unresolved_source_mapping", reason_codes=["no_exact_or_related_result_in_declared_scope"])
            result["exact_alignment_ids"].sort()
            result["related_alignment_ids"].sort()
            output.append(result)
    return sorted(output, key=lambda row: row["source_key"])


def annotate_sequences(study, search_artifacts):
    """Reuse the pinned local gene panel and explicitly versioned local features."""
    from .comparison_provenance import verify_sources
    verify_sources(study)
    features = [{"feature_type": "gene", "name": gene, "contig": contig, "start": start - 1, "end": end,
        "source_id": "Ensembl", "release": annotation.CH_GENE_SPANS_ENSEMBL_RELEASE}
        for gene, (contig, start, end) in annotation.CH_GENE_SPANS.items()]
    manifests = [{"source": annotation.CH_DRIVER_SET_PROVENANCE,
                  "assembly": annotation.CH_GENE_SPANS_ASSEMBLY, "binding": "package source identity"}]
    domain_id = study["config"]["policies"]["comparison"]["domain_id"]
    catalog = {r["record_id"]: r for r in study["config"]["catalogs"]}
    reference = catalog[catalog[domain_id]["reference_bundle_id"]]
    for source in study["sources"].values():
        if source["role"] != "annotations":
            continue
        data = read_json(source["path"])
        fields(data, {"schema_version", "assembly", "source", "release", "features"}, label="local annotation")
        if data["schema_version"] != "1" or data["assembly"] != reference["assembly"]:
            raise PopulationError("Local annotation version/assembly differs from its reference")
        if any(not isinstance(data[key], str) or not data[key].strip() for key in ("source", "release")) \
                or not isinstance(data["features"], list):
            raise PopulationError("Local annotation requires a named source, release and feature list")
        for feature in data["features"]:
            fields(feature, {"feature_type", "name", "contig", "start", "end"}, label="annotation feature")
            integer(feature["start"], "annotation start")
            integer(feature["end"], "annotation end", minimum=feature["start"] + 1)
            if feature["feature_type"] not in {"gene", "regulatory"} or any(not isinstance(feature[k], str) or not feature[k].strip() for k in ("name", "contig")):
                raise PopulationError("Invalid local gene/regulatory feature")
            features.append({**feature, "source_id": data["source"], "release": data["release"]})
        manifests.append({"source": data["source"], "release": data["release"], "sha256": source["sha256"], "assembly": data["assembly"]})
    use_builtin = reference["assembly"] in {"GRCh38", "GRCh38.p14"}
    if not use_builtin:
        features = [row for row in features if row["source_id"] != "Ensembl"]
        manifests = manifests[1:]
    sequences = [row for name, rows in search_artifacts.items() if name.endswith(".jsonl") for row in rows
                 if row["record_type"] == "comparison_sequence_configuration"]
    rows = []
    for sequence in sequences:
        positions = [p["pos0"] for p in sequence["path"] if p["kind"] == "reference"]
        matches = [feature for feature in features if feature["contig"] == sequence["contig"]
                   and any(feature["start"] <= p < feature["end"] for p in positions)]
        rows.append({"sequence_configuration_id": sequence["record_id"], "features": matches,
            "interpretation": "Coordinate annotation only; no activity or outcome estimate"})
    return {"sequence-annotations.json": sorted(rows, key=lambda r: r["sequence_configuration_id"]),
            "annotation-manifest.json": manifests}
