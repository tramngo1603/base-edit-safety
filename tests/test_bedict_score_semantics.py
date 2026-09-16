"""Member-wise BE-DICT normalization, stability gating, and explicit fallback routing."""
from __future__ import annotations

from dataclasses import replace
import math

import pytest

from base_edit_safety.adapters.bedict import (
    NORMALIZATION_APPLIED,
    RATE_SEMANTICS,
    RAW_MASS_MAX,
    RAW_MASS_MIN,
    BedictBystanderPredictor,
    MalformedHaplotypeMass,
    _validated_haplotypes,
)
from base_edit_safety.config import BedictConfig
from base_edit_safety.pipeline.datatypes import Guide, spectrum_to_records


PROTO = "GCTAACGAACGTACGATCAA"
GUIDE = Guide(
    guide_id="HBG1_m113_ABE8e",
    protospacer=PROTO,
    pam="AGG",
    strand="-",
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


def _mutate(seq: str, *positions: int) -> str:
    changed = list(seq)
    for position in positions:
        changed[position - 1] = "G"
    return "".join(changed)


OUTCOME_SEQUENCES = [
    PROTO,
    _mutate(PROTO, 8),
    _mutate(PROTO, 5),
    _mutate(PROTO, 5, 8),
]


def _raw_scores(*scores: float) -> list[dict]:
    return [
        {"output_seq": sequence, "pred_score": score}
        for sequence, score in zip(OUTCOME_SEQUENCES, scores)
    ]


def _raw(first: float, second: float) -> list[dict]:
    return _raw_scores(first, second)


def _response(guide_id: str, member_scores: list[tuple[float, ...]]) -> dict:
    member_rows = [_raw_scores(*scores) for scores in member_scores]
    aggregate = [
        {
            "output_seq": sequence,
            "pred_score": math.fsum(rows[index]["pred_score"] for rows in member_rows)
            / len(member_rows),
        }
        for index, sequence in enumerate(OUTCOME_SEQUENCES[: len(member_scores[0])])
    ]
    return {
        "ensemble_runs": len(member_rows),
        "sequences": {guide_id: aggregate},
        "member_predictions": {
            guide_id: [
                {"member": f"run_{run_num}", "haplotypes": rows}
                for run_num, rows in enumerate(member_rows)
            ]
        },
    }


def _predictor() -> BedictBystanderPredictor:
    return BedictBystanderPredictor(config=BedictConfig("unused", "unused"))


@pytest.mark.parametrize(
    ("scores", "expected"),
    [
        ((0.8, 0.3), (0.8 / 1.1, 0.3 / 1.1)),
        ((0.6, 0.4), (0.6, 0.4)),
        ((1.0, 0.0), (1.0, 0.0)),
    ],
)
def test_one_member_normalization_is_applied_to_every_guide(scores, expected):
    haplotypes, raw_mass = _validated_haplotypes(GUIDE, _raw(*scores))

    assert raw_mass == pytest.approx(sum(scores))
    assert [hap.pred_score for hap in haplotypes] == pytest.approx(expected)
    assert math.fsum(hap.pred_score for hap in haplotypes) == 1.0


def test_stable_mass_anomaly_is_readmitted_with_member_normalization_and_flagged_rates(
    monkeypatch,
):
    member_scores = [
        (0.9, 0.5),
        (0.8, 0.4),
        (0.7, 0.4),
        (0.9, 0.6),
        (0.8, 0.3),
    ]
    response = _response(GUIDE.guide_id, member_scores)

    predictor = _predictor()
    monkeypatch.setattr(predictor, "_run", lambda request: response)
    batch = predictor.predict_many([GUIDE])
    assert batch.exclusions == []
    assert set(batch.spectra) == {GUIDE.guide_id}
    spectrum = batch.spectra[GUIDE.guide_id]
    expected_reference = math.fsum(
        reference / (reference + edited) for reference, edited in member_scores
    ) / len(member_scores)

    assert math.fsum(hap.pred_score for hap in spectrum.haplotypes) == 1.0
    assert spectrum.haplotypes[0].pred_score == pytest.approx(expected_reference)
    assert spectrum.provenance["normalization_applied"] == NORMALIZATION_APPLIED == \
        "per_member_then_ensemble_v1"
    assert spectrum.provenance["structural_stability"] == "stable"
    assert spectrum.provenance["upstream_mass_anomaly"] is True
    assert spectrum.provenance["raw_mass"] == pytest.approx(1.26)
    assert spectrum.provenance["raw_mass_per_member"] == pytest.approx({
        "run_0": 1.4,
        "run_1": 1.2,
        "run_2": 1.1,
        "run_3": 1.5,
        "run_4": 1.1,
    })
    assert spectrum.provenance["rate_semantics"] == RATE_SEMANTICS == "model_relative"

    records = spectrum_to_records(spectrum, on_target_offset=113)
    assert records
    for record in records:
        assert record.provenance["upstream_mass_anomaly"] is True
        assert record.provenance["raw_mass_per_member"] == \
            spectrum.provenance["raw_mass_per_member"]
        assert record.provenance["rate_semantics"] == "model_relative"


def test_mass_only_case_inside_new_envelope_no_longer_auto_excludes():
    haplotypes, raw_mass = _validated_haplotypes(GUIDE, _raw(0.8, 0.6))

    assert raw_mass == pytest.approx(1.4)
    assert math.fsum(hap.pred_score for hap in haplotypes) == 1.0


@pytest.mark.parametrize(
    "raw",
    [
        _raw(0.3, 0.2),
        _raw(1.0, 1.0),
    ],
)
def test_new_raw_mass_bounds_are_inclusive(raw):
    haplotypes, raw_mass = _validated_haplotypes(GUIDE, raw)

    assert raw_mass in (RAW_MASS_MIN, RAW_MASS_MAX)
    assert math.fsum(hap.pred_score for hap in haplotypes) == 1.0


@pytest.mark.parametrize(
    "raw",
    [
        _raw(0.3, 0.199999),
        _raw_scores(0.8, 0.7, 0.500001),
    ],
)
def test_raw_mass_outside_new_bounds_is_excluded_after_stability_check(raw):
    with pytest.raises(MalformedHaplotypeMass, match="after stability check"):
        _validated_haplotypes(GUIDE, raw)


@pytest.mark.parametrize("score", [-0.01, 1.01, float("nan"), float("inf")])
def test_non_finite_or_out_of_range_raw_score_is_rejected_before_normalization(score):
    with pytest.raises(ValueError, match=r"inside \[0, 1\]"):
        _validated_haplotypes(GUIDE, _raw(score, 0.5))


def test_unstable_rank_order_is_excluded_and_routed_to_crispron_abe(monkeypatch):
    response = _response(
        GUIDE.guide_id,
        [
            (0.8, 0.2),
            (0.4, 0.6),
            (0.7, 0.3),
            (0.75, 0.25),
            (0.45, 0.55),
        ],
    )
    predictor = _predictor()
    monkeypatch.setattr(predictor, "_run", lambda request: response)

    batch = predictor.predict_many([GUIDE])

    assert batch.spectra == {}
    assert len(batch.exclusions) == 1
    exclusion = batch.exclusions[0]
    assert exclusion["status"] == "excluded_structurally_unstable_ensemble"
    assert exclusion["structural_stability"] == "unstable"
    assert exclusion["normalization_applied"] is None
    assert exclusion["fallback"] == {
        "predictor": "CRISPRon-ABE",
        "analysis_role": "sensitivity_fallback",
        "status": "routed_not_executed",
    }
    assert exclusion["member_rank_order"]["run_0"] != \
        exclusion["member_rank_order"]["run_1"]


def test_stable_ensemble_outside_new_mass_gate_is_excluded(monkeypatch):
    malformed = replace(GUIDE, guide_id="HBG1_malformed_ABE8e")
    response = _response(malformed.guide_id, [(0.2, 0.2)] * 5)
    predictor = _predictor()
    monkeypatch.setattr(predictor, "_run", lambda request: response)

    batch = predictor.predict_many([malformed])

    assert batch.spectra == {}
    assert len(batch.exclusions) == 1
    exclusion = batch.exclusions[0]
    assert exclusion["guide_id"] == malformed.guide_id
    assert exclusion["raw_mass"] == pytest.approx(0.4)
    assert exclusion["raw_mass_per_member"] == pytest.approx({
        f"run_{run_num}": 0.4 for run_num in range(5)
    })
    assert exclusion["allowed_raw_mass"] == {"minimum": 0.5, "maximum": 2.0}
    assert exclusion["structural_stability"] == "stable"
    assert exclusion["normalization_applied"] is None
    assert exclusion["status"] == "excluded_raw_mass_outside_authorized_range"
