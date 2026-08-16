"""Use complete joint edited sequences for multi-edit efficacy, not a position-level maximum."""
from __future__ import annotations

import copy

import pytest

from base_edit_safety.design_record import SCHEMA_VERSION, build_design_records
from base_edit_safety.io import outputs
from base_edit_safety.io.outputs import records_to_dicts
from base_edit_safety.pipeline.datatypes import (
    JOINT_ON_TARGET_CAVEAT,
    JOINT_ON_TARGET_METHOD,
    EditSpectrum,
    Guide,
    OutcomeHaplotype,
    spectrum_to_records,
)
from base_edit_safety.scoring.rubric import score_guides


PROTO = "GGCAAGGCTGGCCAACCCAT"
PAIR_GUIDE = Guide(
    guide_id="HBG1_pair_ABE8e",
    protospacer=PROTO,
    pam="GGG",
    strand="+",
    genome_build="GRCh38",
    chrom="chr11",
    gene="HBG1",
    gene_strand="-",
    tss_1based=5_249_857,
    editor_class="ABE",
    editor_name="ABE8e",
    edit_from="A",
    edit_to="G",
    proto_5p_genomic=5_249_977,
)
PAIR_POSITIONS = {4, 5}  # promoter alleles -123 T>C and -124 T>C
PAIR_COORDS = {PAIR_GUIDE.genomic_coord(pos) for pos in PAIR_POSITIONS}


def _mutate(seq: str, *positions: int) -> str:
    chars = list(seq)
    for pos in positions:
        chars[pos - 1] = "G"
    return "".join(chars)


def _spectrum() -> EditSpectrum:
    # Complete probability mass: marginals are 0.45 and 0.35; their joint is 0.20.
    return EditSpectrum(
        guide=PAIR_GUIDE,
        model_name="BE-DICT-bystander:fixture",
        haplotypes=[
            OutcomeHaplotype(PROTO, 0.40),
            OutcomeHaplotype(_mutate(PROTO, 4), 0.25),
            OutcomeHaplotype(_mutate(PROTO, 5), 0.15),
            OutcomeHaplotype(_mutate(PROTO, 4, 5), 0.20),
        ],
        provenance={"ensemble_runs": 5},
    )


def _pair_records() -> list[dict]:
    return records_to_dicts(spectrum_to_records(_spectrum(), on_target_coords=PAIR_COORDS))


def _on_target(records: list[dict]) -> list[dict]:
    return [record for record in records if record["on_target_or_bystander"] == "on_target"]


def test_joint_rate_is_no_greater_than_the_smaller_marginal():
    records = _pair_records()
    marginals = [record["predicted_rate"] for record in _on_target(records)]
    joint = records[0]["provenance"]["joint_on_target_rate"]

    assert joint == 0.20
    assert joint <= min(marginals)


def test_joint_rate_is_no_greater_than_the_larger_marginal():
    records = _pair_records()
    marginals = [record["predicted_rate"] for record in _on_target(records)]
    joint = records[0]["provenance"]["joint_on_target_rate"]

    assert joint <= max(marginals)


def test_joint_rate_requires_the_intended_allele_not_merely_any_change_at_each_position():
    wrong = list(PROTO)
    wrong[3] = "C"
    wrong[4] = "T"
    spectrum = EditSpectrum(
        guide=PAIR_GUIDE,
        model_name="fixture",
        haplotypes=[OutcomeHaplotype("".join(wrong), 1.0)],
    )

    assert set(spectrum.edited_positions(spectrum.haplotypes[0])) >= PAIR_POSITIONS
    assert spectrum.joint_edit_rate(PAIR_POSITIONS) == 0.0


def test_multi_edit_band_consumes_joint_rate_not_max_or_min():
    score = score_guides(_pair_records())[0]

    assert score.on_target_efficiency == 0.20
    assert score.on_target_efficiency != 0.45
    assert score.on_target_efficiency != 0.35
    assert score.efficacy_band == "PRECLINICAL_ACTIVITY_ONLY"
    assert score.terms["efficiency_basis"] == "joint_on_target_rate"


def test_single_edit_fallback_is_unchanged_and_adds_no_joint_field():
    one_coord = {PAIR_GUIDE.genomic_coord(4)}
    records = records_to_dicts(spectrum_to_records(_spectrum(), on_target_coords=one_coord))
    score = score_guides(records)[0]

    assert score.on_target_efficiency == 0.45
    assert score.efficacy_band == "MEETS_CLINICAL_PRECEDENT"
    assert score.terms["efficiency_basis"] == "single_edit_rate"
    assert all("joint_on_target_rate" not in record["provenance"] for record in records)
    design = build_design_records(records)["records"][0]
    assert "joint_on_target_rate" not in design["predicted_edit_spectrum"]
    assert JOINT_ON_TARGET_CAVEAT not in design["limitations"]


def test_multi_edit_design_record_keeps_marginals_and_allele_keyed_joint_provenance():
    design = build_design_records(_pair_records())["records"][0]
    spectrum = design["predicted_edit_spectrum"]
    joint_provenance = spectrum["joint_on_target_rate_provenance"]

    assert design["schema_version"] == SCHEMA_VERSION == "2"
    assert [row["predicted_rate"] for row in spectrum["on_target"]] == [0.45, 0.35]
    assert spectrum["joint_on_target_rate"] == 0.20
    assert joint_provenance["allele_keys"] == [
        {"promoter_offset": -123, "promoter_sense_ref": "T", "promoter_sense_alt": "C"},
        {"promoter_offset": -124, "promoter_sense_ref": "T", "promoter_sense_alt": "C"},
    ]
    assert joint_provenance["method"] == JOINT_ON_TARGET_METHOD
    assert joint_provenance["model"] == "BE-DICT-bystander:fixture"
    assert joint_provenance["predictor_runs"] == 5
    assert joint_provenance["confidence"] == "model_predicted"
    assert joint_provenance["caveat"] == JOINT_ON_TARGET_CAVEAT
    assert JOINT_ON_TARGET_CAVEAT in design["limitations"]
    assert all("allele_key" in row for row in design["registry_match"]["by_edit"])


def test_multi_edit_without_joint_data_is_not_evaluated_instead_of_using_a_marginal():
    records = _pair_records()
    for record in records:
        for key in list(record["provenance"]):
            if key.startswith("joint_on_target_"):
                del record["provenance"][key]

    score = score_guides(records)[0]

    assert score.on_target_efficiency is None
    assert score.efficacy_band == "NOT_EVALUATED"
    assert score.terms["efficiency_basis"] == "joint_on_target_rate_not_evaluated"


def test_partial_joint_provenance_is_not_enough_to_evaluate_a_multi_edit_design():
    records = _pair_records()
    on_target = _on_target(records)
    for key in list(on_target[1]["provenance"]):
        if key.startswith("joint_on_target_"):
            del on_target[1]["provenance"][key]

    score = score_guides(records)[0]
    design = build_design_records(records)["records"][0]

    assert score.on_target_efficiency is None
    assert score.terms["efficiency_basis"] == "joint_on_target_rate_not_evaluated"
    assert design["predicted_edit_spectrum"]["joint_on_target_rate"] is None
    assert design["predicted_edit_spectrum"]["joint_on_target_rate_provenance"]["confidence"] == \
        "not_evaluated"


def test_joint_provenance_must_name_each_intended_allele_once():
    records = _pair_records()
    duplicated = [
        {"promoter_offset": -123, "promoter_sense_ref": "T", "promoter_sense_alt": "C"},
        {"promoter_offset": -123, "promoter_sense_ref": "T", "promoter_sense_alt": "C"},
    ]
    for record in records:
        record["provenance"]["joint_on_target_allele_keys"] = duplicated

    with pytest.raises(ValueError, match="allele keys do not match"):
        score_guides(records)


@pytest.mark.parametrize(
    ("missing_key", "message"),
    [
        ("joint_on_target_model", "model provenance"),
        ("predictor_runs", "positive predictor run count"),
        ("joint_on_target_caveat", "model-predicted-not-measured caveat"),
    ],
)
def test_joint_rate_requires_complete_model_provenance(missing_key: str, message: str):
    records = _pair_records()
    for record in records:
        del record["provenance"][missing_key]

    with pytest.raises(ValueError, match=message):
        score_guides(records)


def test_joint_rate_above_a_marginal_is_rejected():
    records = copy.deepcopy(_pair_records())
    for record in records:
        record["provenance"]["joint_on_target_rate"] = 0.36

    with pytest.raises(ValueError, match="exceeds an intended edit's marginal"):
        score_guides(records)


def test_enumeration_schema_bumped_for_joint_provenance_contract():
    assert outputs.ARTIFACT_SCHEMA_VERSION == "4"
