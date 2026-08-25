"""Uniform BE-DICT score normalization with explicit malformed-mass exclusions."""
from __future__ import annotations

from dataclasses import replace

import pytest

from base_edit_safety.adapters.bedict import (
    NORMALIZATION_APPLIED,
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


def _mutate(seq: str, position: int) -> str:
    return seq[: position - 1] + "G" + seq[position:]


def _raw(first: float, second: float) -> list[dict]:
    return [
        {"output_seq": PROTO, "pred_score": first},
        {"output_seq": _mutate(PROTO, 8), "pred_score": second},
    ]


def _predictor() -> BedictBystanderPredictor:
    return BedictBystanderPredictor(config=BedictConfig("unused", "unused"))


@pytest.mark.parametrize(
    ("scores", "expected"),
    [
        ((0.8, 0.3), (0.8 / 1.1, 0.3 / 1.1)),
        ((0.6, 0.4), (0.6, 0.4)),
        ((1.1, 0.0), (1.0, 0.0)),
    ],
)
def test_uniform_normalization_is_applied_to_every_guide(scores, expected):
    haplotypes, raw_mass = _validated_haplotypes(GUIDE, _raw(*scores))

    assert raw_mass == pytest.approx(sum(scores))
    assert [hap.pred_score for hap in haplotypes] == pytest.approx(expected)
    assert sum(hap.pred_score for hap in haplotypes) == pytest.approx(1.0)


def test_raw_mass_and_uniform_policy_are_preserved_on_every_derived_record():
    response = {
        "ensemble_runs": 5,
        "sequences": {GUIDE.guide_id: _raw(0.66, 0.44)},
    }
    spectrum = _predictor().spectrum_from_response(GUIDE, "ABE8e", response)
    records = spectrum_to_records(spectrum, on_target_offset=113)

    assert spectrum.provenance["normalization_applied"] == NORMALIZATION_APPLIED == "uniform_v1"
    assert spectrum.provenance["raw_mass"] == pytest.approx(1.1)
    assert records
    assert all(record.provenance["normalization_applied"] == "uniform_v1" for record in records)
    assert all(record.provenance["raw_mass"] == pytest.approx(1.1) for record in records)
    on_target = next(record for record in records if record.on_target_or_bystander == "on_target")
    assert on_target.predicted_rate == 0.4


@pytest.mark.parametrize("raw_mass", [0.75, 1.25])
def test_raw_mass_bounds_are_inclusive(raw_mass):
    first = min(raw_mass, 0.7)
    haplotypes, preserved = _validated_haplotypes(GUIDE, _raw(first, raw_mass - first))

    assert preserved == pytest.approx(raw_mass)
    assert sum(hap.pred_score for hap in haplotypes) == pytest.approx(1.0)


@pytest.mark.parametrize("raw_mass", [0.749999, 1.250001])
def test_raw_mass_outside_bounds_is_rejected_as_malformed(raw_mass):
    first = min(raw_mass, 0.7)
    with pytest.raises(MalformedHaplotypeMass, match="guide excluded as malformed") as exc:
        _validated_haplotypes(GUIDE, _raw(first, raw_mass - first))

    assert exc.value.raw_mass == pytest.approx(raw_mass)


@pytest.mark.parametrize("score", [-0.01, float("nan"), float("inf")])
def test_non_finite_or_negative_raw_score_is_rejected_before_normalization(score):
    with pytest.raises(ValueError, match="must be finite and non-negative"):
        _validated_haplotypes(GUIDE, _raw(score, 1.0))


def test_batch_exclusion_is_returned_with_raw_mass_instead_of_silently_dropped(monkeypatch):
    malformed = replace(GUIDE, guide_id="HBG1_malformed_ABE8e")
    response = {
        "ensemble_runs": 5,
        "sequences": {
            GUIDE.guide_id: _raw(0.66, 0.44),
            malformed.guide_id: _raw(0.4, 0.3),
        },
    }
    predictor = _predictor()
    monkeypatch.setattr(predictor, "_run", lambda request: response)

    batch = predictor.predict_many([GUIDE, malformed])

    assert set(batch.spectra) == {GUIDE.guide_id}
    assert len(batch.exclusions) == 1
    exclusion = batch.exclusions[0]
    assert exclusion["guide_id"] == malformed.guide_id
    assert exclusion["gene"] == "HBG1"
    assert exclusion["editor_class"] == "ABE"
    assert exclusion["editor"] == "ABE8e"
    assert exclusion["protospacer"] == PROTO
    assert exclusion["raw_mass"] == pytest.approx(0.7)
    assert exclusion["allowed_raw_mass"] == {"minimum": 0.75, "maximum": 1.25}
    assert exclusion["normalization_applied"] is None
    assert exclusion["status"] == "excluded_malformed_raw_mass"
    assert "guide excluded as malformed" in exclusion["reason"]
