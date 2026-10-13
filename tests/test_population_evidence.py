"""Exact source identity and assay context determine evidence applicability."""
from copy import deepcopy
import math
from pathlib import Path

import pytest

from base_edit_safety.pipeline.candidate_panel import _alignment_block_for_site
from base_edit_safety.population.actions import make_actions, make_request, make_units
from base_edit_safety.population.common import PopulationError, read_json
from base_edit_safety.population.evidence import evidence_for_unit, resolve_measurements, evidence_summary

ROOT = Path(__file__).resolve().parents[1]


def inputs():
    config = read_json(ROOT/"examples/population/synthetic/study-config.json")
    candidate = {"configuration_id": "cfg:synthetic", "guide_definition_id": "guide:synthetic", "editor_definition_id": "editor:synthetic",
                 "guide_id": "SYN_guide", "contig": "chrSYN", "strand": "+", "bounding_interval": [10, 33], "role": "candidate_offtarget",
                 "target_path": [{"kind": "reference", "pos0": i} for i in range(10, 33)],
                 "sequence_genomic_forward": "A"*20+"TGG"}
    guide = {"guide_id": "SYN_guide", "spacer_sequence": "A"*20, "editor": "SYN_ABE",
             "editor_class": "ABE", "pam_class": "NGG", "editing_window_start_1based": 4,
             "editing_window_end_1based_inclusive": 8}
    unit = make_units(config, [candidate], "test:run")[0]
    measurement = {"editor": "SYN_ABE", "guide_id": "SYN_guide", "chrom": "chrSYN", "strand": "+", "start_0based": 10,
                   "end_0based_exclusive": 33, "genome_build": "synthetic_v3_search_1", "allele_sequence": candidate["sequence_genomic_forward"],
                   "allele_resolution": "resolved", "allele_evidence": "Synthetic allele-resolved truth", "qc_status": "passed",
                   "treated_group": "synthetic treated", "control_group": "synthetic control", "biological_context": "synthetic cells",
                   "assay": "synthetic assay", "measurement_kind": "editing", "status": "not_detected",
                   "limit_of_detection": {"value": 0.1, "unit": "percent"}, "source_record_id": "SYN_record", "source_version_or_citation": "synthetic"}
    link = {"configuration_id": candidate["configuration_id"], "target_path": candidate["target_path"], "locus_resolved": True,
            "donor_id": "SYN_person", "person_identity_status": "resolved", "biological_context_id": config["assay_policy"]["biological_context_id"],
            "context": deepcopy(config["assay_policy"]["context"])}
    return config, unit, candidate, measurement, link, guide


def test_reference_negative_never_covers_alt():
    c, u, target, m, link, guide = inputs()
    assert evidence_for_unit(c, u, target, m, link, guide, "test:run", source_guide=guide)["usable"]
    target["sequence_genomic_forward"] = "C" + target["sequence_genomic_forward"][1:]
    e = evidence_for_unit(c, u, target, m, link, guide, "test:run", source_guide=guide)
    assert not e["usable"] and "allele_sequence_mismatch" in e["reasons"]


@pytest.mark.parametrize("field,value,reason", [
    ("spacer_sequence", "C"*20, "source_guide_definition_mismatch"),
    ("editor_class", "CBE", "source_editor_definition_mismatch"),
    ("pam_class", "NAG", "source_editor_definition_mismatch"),
    ("editing_window_end_1based_inclusive", 9, "source_editor_definition_mismatch"),
])
def test_source_guide_and_editor_definitions_must_match(field, value, reason):
    c, u, target, m, link, guide = inputs()
    source = {**guide, field: value}
    e = evidence_for_unit(c, u, target, m, link, guide, "test:run", source_guide=source)
    assert not e["usable"] and reason in e["reasons"]
    assert e["source_guide"] == source and e["source_measurement"] == m


def test_absent_source_guide_never_supplies_measurement_credit():
    c, u, target, m, link, guide = inputs()
    e = evidence_for_unit(c, u, target, m, link, guide, "test:run")
    assert not e["usable"] and "source_guide_definition_unknown" in e["reasons"]


@pytest.mark.parametrize("field,value,reason", [("editor", "another editor", "editor_mismatch"),
    ("allele_resolution", "mixed", "allele_resolution_unresolved"), ("qc_status", None, "qc_not_passed"),
    ("control_group", None, "control_group_missing"), ("limit_of_detection", None, "numeric_detection_limit_missing"),
    ("limit_of_detection", {"value": True, "unit": "percent"}, "numeric_detection_limit_missing")])
def test_unusable_measurement_is_retained(field, value, reason):
    c, u, target, m, link, guide = inputs(); m[field] = value
    e = evidence_for_unit(c, u, target, m, link, guide, "test:run", source_guide=guide)
    assert not e["usable"] and reason in e["reasons"]
    assert e["outcome"] == "not_detected"


def test_nonfinite_measurements_fail_serialization_instead_of_becoming_evidence():
    c, u, target, m, link, guide = inputs()
    m["limit_of_detection"]["value"] = float("nan")
    with pytest.raises(ValueError):
        evidence_for_unit(c, u, target, m, link, guide, "test:run", source_guide=guide)


@pytest.mark.parametrize("value,unit,fraction,reason", [
    (0.1, "percent", 0.001, None), (0.1, " % ", 0.001, None),
    (100, "percentage", 1, None), (0.001, "fraction", 0.001, None),
    (0.001, "Proportion", 0.001, None), (200, "percent", None, "detection_limit_out_of_bounds"),
    (1.1, "fraction", None, "detection_limit_out_of_bounds"),
    (0.1, "arbitrary dimension", None, "detection_limit_unit_unsupported"),
    (10, "copies", None, "detection_limit_unit_unsupported"),
])
def test_detection_limits_are_editing_fractions_with_valid_bounds(value, unit, fraction, reason):
    c, u, target, m, link, guide = inputs()
    m["limit_of_detection"] = {"value": value, "unit": unit}
    e = evidence_for_unit(c, u, target, m, link, guide, "test:run", source_guide=guide)
    assert e["usable"] is (reason is None)
    assert e["detection_limit_assessment"]["normalized_fraction"] == fraction
    assert e["detection_limit_assessment"]["reason"] == reason
    assert e["source_measurement"]["limit_of_detection"] == {"value": value, "unit": unit}


def test_study_sensitivity_requirement_changes_evidence_and_coverage_identity():
    from base_edit_safety.population.contracts import resolution_id, validate_config
    c, u, target, m, link, guide = inputs()
    c["assay_policy"]["max_detection_limit_fraction"] = 0.0005
    with pytest.raises(PopulationError, match="resolution identity"):
        validate_config(c)
    c["assay_policy"]["resolution_contract_id"] = resolution_id(c["assay_policy"])
    validate_config(c)
    sensitive_unit = make_units(c, [target], "test:run")[0]
    assert sensitive_unit["id"] != u["id"]
    e = evidence_for_unit(c, sensitive_unit, target, m, link, guide, "test:run", source_guide=guide)
    assert not e["usable"] and e["reasons"] == ["detection_limit_exceeds_study_maximum"]
    m["limit_of_detection"] = {"value": 0.05, "unit": "percent"}
    assert evidence_for_unit(c, sensitive_unit, target, m, link, guide, "test:run", source_guide=guide)["usable"]


@pytest.mark.parametrize("percent,fraction,accepted", [
    (0.07, 0.0007, True),
    (math.nextafter(0.07, -math.inf), math.nextafter(0.0007, -math.inf), True),
    (math.nextafter(0.07, math.inf), math.nextafter(0.0007, math.inf), False),
])
def test_inclusive_detection_limit_uses_decimal_values_without_tolerance(percent, fraction, accepted):
    from base_edit_safety.measurement_units import editing_detection_limit
    limits = [editing_detection_limit({"value": value, "unit": unit}, 0.0007)
              for value, unit in [(percent, "percent"), (fraction, "fraction")]]
    assert all((limit["status"] == "accepted") is accepted for limit in limits)
    assert all(limit["reason"] == (None if accepted else "detection_limit_exceeds_study_maximum") for limit in limits)
    if percent == 0.07:
        assert limits[0] == limits[1]


def test_detection_limit_conversion_does_not_depend_on_global_decimal_precision():
    from decimal import localcontext
    from base_edit_safety.measurement_units import editing_detection_limit
    with localcontext() as context:
        context.prec = 2
        limit = editing_detection_limit({"value": 0.07002, "unit": "percent"}, 0.0007)
    assert limit["status"] == "insufficient" and limit["normalized_fraction"] == 0.0007002


def test_detection_limit_underflow_stays_unresolved_instead_of_becoming_zero():
    from base_edit_safety.measurement_units import editing_detection_limit
    limit = editing_detection_limit({"value": math.ulp(0.0), "unit": "percent"}, 0.0007)
    assert limit["status"] == "unresolved" and limit["normalized_fraction"] is None
    assert limit["reason"] == "detection_limit_below_numeric_range"


@pytest.mark.parametrize("maximum", [0, -0.1, 1.1, True, "0.001"])
def test_invalid_study_detection_limit_is_rejected(maximum):
    from base_edit_safety.population.contracts import validate_config
    c, *_ = inputs()
    c["assay_policy"]["max_detection_limit_fraction"] = maximum
    with pytest.raises(PopulationError):
        validate_config(c)


@pytest.mark.parametrize("field", ["cell_type", "delivery", "condition"])
def test_context_mismatch_even_when_context_label_matches(field):
    c, u, target, m, link, guide = inputs(); link["context"][field] = "different"
    e = evidence_for_unit(c, u, target, m, link, guide, "test:run", source_guide=guide)
    assert not e["usable"] and e["context_match"] == "mismatch"


def test_unknown_person_cannot_add_distinct_donor_credit():
    c, u, target, m, link, guide = inputs(); link.update(donor_id=None, person_identity_status="ambiguous")
    e = evidence_for_unit(c, u, target, m, link, guide, "test:run", source_guide=guide)
    assert not e["usable"] and e["donor_id"] is None


def test_conflicting_outcomes_are_not_averaged():
    c, u, target, m, link, guide = inputs()
    negative = evidence_for_unit(c, u, target, m, link, guide, "test:run", source_guide=guide)
    m = {**m, "source_record_id": "SYN_second", "status": "edited"}
    positive = evidence_for_unit(c, u, target, m, link, guide, "test:run", source_guide=guide)
    assert negative["usable"] and positive["usable"]
    assert negative["evidence_id"] != positive["evidence_id"]
    assert {negative["outcome"], positive["outcome"]} == {"edited", "not_detected"}
    excluded = {**negative, "usable": False, "reasons": ["qc_not_passed"]}
    summary = evidence_summary([negative, positive, excluded], [{"reasons": ["source_guide_id_not_in_study"]}])
    assert summary == {"usable_positive_links": 1, "usable_nondetection_links": 1,
        "conflicting_donor_unit_pairs": 1, "excluded_links": 1, "excluded_link_reasons": {"qc_not_passed": 1},
        "unmatched_source_measurements": 1, "unmatched_source_reasons": {"source_guide_id_not_in_study": 1}}


def test_assay_resolution_and_context_do_not_follow_from_carriage():
    c, unit, candidate, *_ = inputs()
    donor = {"donor_id": "SYN", "specimen_id": "SYN_specimen", "use_permitted": True, "simulation_eligible": True,
             "specimen_context_id": c["assay_policy"]["specimen_context_id"], "availability": "verified_available", "activation_cost": 1}
    carriage = [{"donor_id": "SYN", "configuration_id": candidate["configuration_id"], "evidence_supported_carriage": "confirmed", "model_carriage": "present"}]
    c["assay_policy"]["locus_resolved"] = False
    action = make_actions(c, [unit], [candidate], [donor], carriage, "test:run")[0]
    assert action["eligibility"] == "requires_review" and not action["covered_units"]
    c["assay_policy"]["locus_resolved"] = True
    donor["specimen_context_id"] = "other_context"
    assert make_actions(c, [unit], [candidate], [donor], carriage, "test:run")[0]["eligibility"] == "ineligible"


@pytest.mark.parametrize("strand", ["+", "-"])
def test_guarded_experimental_reference_representation(strand):
    nomination = {"nomination_kind": "experimental", "site_sequence": "A"*20+"TGG", "pam": "TGG", "strand": strand,
                  "source_method": "synthetic assay", "reference_validation": {"status": "verified"}}
    site = {"chrom": "chrSYN", "start_0based": 10, "end_0based_exclusive": 33, "nominations": [nomination]}
    guide = {"spacer_sequence": "A"*20}
    block = _alignment_block_for_site(site, guide)
    assert block is not None and block["allele_states"][0]["allele"] == "ref"
    assert "donor_id" not in block
    for mutation in [{"reference_validation": {"status": "not_verified"}}, {"site_sequence": "N"*20+"TGG"},
                     {"strand": None}, {"pam": None}, {"population_variant_id": "synthetic_alt"}, {"dna_bulge": 1}]:
        site["nominations"] = [{**nomination, **mutation}]
        assert _alignment_block_for_site(site, guide) is None
    site["nominations"] = [nomination, {**nomination, "site_sequence": "C"*20+"TGG"}]
    assert _alignment_block_for_site(site, guide) is None


def test_hbg2_only_assay_cannot_cover_hbg1():
    c, unit, candidate, *_ = inputs()
    c["assay_policy"]["locus_scope"] = {"mode": "explicit_loci", "loci": [{"contig": "chrSYN", "start": 4934, "end": 4957, "strand": "+"}]}
    donor = {"donor_id": "SYN", "specimen_id": "SYN_specimen", "use_permitted": True, "simulation_eligible": True,
             "specimen_context_id": c["assay_policy"]["specimen_context_id"], "availability": "verified_available", "activation_cost": 1}
    rows = [{"donor_id": "SYN", "configuration_id": candidate["configuration_id"], "evidence_supported_carriage": "confirmed", "model_carriage": "present"}]
    action = make_actions(c, [unit], [candidate], [donor], rows, "test:run")[0]
    assert not action["covered_units"] and "assay_locus_not_declared" in action["eligibility_reasons"]


def test_placeholder_control_is_not_evidence():
    c, u, target, m, link, guide = inputs(); m["control_group"] = "not_reported"
    assert not evidence_for_unit(c, u, target, m, link, guide, "test:run", source_guide=guide)["usable"]


def test_representation_replay_keeps_original_measurements_and_zero_new_donor_credit():
    from base_edit_safety.population.evidence import representation_replay
    from base_edit_safety.population.identity import reverse_complement
    class Reference:
        sha256 = "1"*64
        def fetch(self, contig, start, end):
            seq = "A"*20+"TGG"
            return seq if start == 10 else reverse_complement(seq)
    sites = []
    for start, strand in [(10, "+"), (50, "-")]:
        sites.append({"site_id": "synthetic_"+strand, "guide_id": "SYN", "chrom": "chrSYN", "start_0based": start,
            "end_0based_exclusive": start+23, "allele_states": [], "measurements": [{"status": "not_detected", "qc_status": None}],
            "nominations": [{"nomination_kind": "experimental", "source_method": "synthetic", "site_sequence": "A"*20+"TGG", "pam": "TGG", "strand": strand}]})
    ledger = {"candidates": sites}; before = deepcopy(ledger)
    rows = representation_replay(ledger, Reference(), [{"guide_id": "SYN", "spacer_sequence": "A"*20}], "test:run")
    assert len(rows) == 2 and all(r["after"] == "verified_reference_model" for r in rows)
    assert all(r["new_donor_measurement_credits"] == 0 for r in rows)
    assert ledger == before
