"""Parsing and validation against hand-written synthetic fixture files.

Nothing here installs a tool, starts a subprocess, or touches the network.
"""
import csv
import json
import pathlib

import pytest

from base_edit_safety.adapters.bedict import BedictBystanderPredictor, _validated_haplotypes
from base_edit_safety.adapters.crispor import AmbiguousTargetMapping, CrisporGuideDesigner
from base_edit_safety.config import BedictConfig, CrisporConfig
from base_edit_safety.io.reference import ReferenceWindow
from base_edit_safety.pipeline.datatypes import Guide, reverse_complement, spectrum_to_records

FIXTURES = pathlib.Path(__file__).parent / "fixtures"

SYNTHETIC_PROTO = "GCTAACGAACGTACGATCAA"
SYNTHETIC_TARGET = SYNTHETIC_PROTO + "AGG"

# Synthetic sequence with the same coordinate geometry used by the mapping assertions.
HBG1_WIN = ReferenceWindow(
    genome_build="GRCh38", chrom="chr11", start_1based=5_249_937, end_1based=5_250_067,
    plus_sequence=reverse_complement("A" * 90 + SYNTHETIC_TARGET + "C" * 18),
)

GUIDE_M113 = Guide(
    guide_id="HBG1_m113_ABE", protospacer=SYNTHETIC_PROTO, pam="AGG", strand="-",
    genome_build="GRCh38", chrom="chr11", gene="HBG1", gene_strand="-", tss_1based=5_249_857,
    editor_class="ABE", editor_name="ABE8e", edit_from="A", edit_to="G", proto_5p_genomic=5_249_977,
)


def _predictor():
    # a config is passed in so __init__ does not go looking for a real predictor
    return BedictBystanderPredictor(config=BedictConfig("x", "y"))


def _designer():
    return CrisporGuideDesigner(config=CrisporConfig("x", "y", "hg38", "z"))


def _rows(name):
    with open(FIXTURES / name) as handle:
        return [dict(r) for r in csv.DictReader(handle, delimiter="\t")]


# --- BE-DICT response parsing + validation--------------------------------------

def test_bedict_response_parses_to_expected_records():
    response = json.loads((FIXTURES / "bedict_response_HBG1_m113_ABE.json").read_text())
    spectrum = _predictor().spectrum_from_response(GUIDE_M113, "ABE8e", response)
    records = spectrum_to_records(spectrum, on_target_offset=113)
    scored = [r for r in records if r.confidence == "model_predicted"]
    by_role = {r.on_target_or_bystander: r for r in scored}
    assert by_role["on_target"].promoter_offset == 113
    assert abs(by_role["on_target"].predicted_rate - 0.25) < 1e-9   # position 8: 0.20 + 0.05
    assert abs(by_role["bystander"].predicted_rate - 0.20) < 1e-9   # position 5: 0.15 + 0.05
    # Every scored record sits inside the window; the rest say the predictor was never asked.
    assert all(3 <= r.provenance["protospacer_position"] <= 10 for r in scored)
    unscored = [r for r in records if r.confidence == "not_evaluated"]
    assert unscored and all(r.predicted_rate is None for r in unscored)


def test_bedict_rejects_wrong_length_haplotype():
    with pytest.raises(ValueError, match="length"):
        _validated_haplotypes(GUIDE_M113, [{"output_seq": "ACGT", "pred_score": 0.5}])


def test_bedict_rejects_probability_out_of_range():
    with pytest.raises(ValueError, match="inside"):
        _validated_haplotypes(GUIDE_M113, [{"output_seq": "A" * 20, "pred_score": 1.5}])


def test_bedict_mass_above_superseded_envelope_is_normalized_inside_new_gate():
    raw = [{"output_seq": "A" * 20, "pred_score": 0.7}, {"output_seq": "C" * 20, "pred_score": 0.7}]
    haplotypes, raw_mass = _validated_haplotypes(GUIDE_M113, raw)

    assert raw_mass == pytest.approx(1.4)
    assert sum(haplotype.pred_score for haplotype in haplotypes) == pytest.approx(1.0)


# --- CRISPOR row parsing + off-target summary-----------------------------------

def test_crispor_guide_row_maps_and_summarizes_offtargets():
    guide_rows = _rows("crispor_guides.tsv")
    sites_by_guide = {}
    for s in _rows("crispor_offtargets.tsv"):
        sites_by_guide.setdefault(s["guideId"], []).append(s)
    row = next(r for r in guide_rows if r["guideId"] == "19rev")
    cand = _designer()._row_to_candidate(row, HBG1_WIN, "HBG1", sites_by_guide)
    assert cand.strand == "-" and cand.proto_5p_genomic == 5_249_977
    assert cand.offtarget["mm_histogram"] == {"0": 1, "2": 1, "4": 2}
    assert cand.offtarget["perfect_match_candidate_count"] == 1
    assert cand.offtarget["paralog_count"] == 0
    assert cand.offtarget["offtarget_count_excl_paralog"] == 4


def test_crispor_ambiguous_target_is_rejected():
    # A target that shows up twice in the window must not quietly map to the first hit.
    target = "ACGTACGTACGTACGTACGTAGG"  # 20 nt protospacer + AGG PAM
    win = ReferenceWindow("GRCh38", "chr11", 1, 200,
                          plus_sequence="TTTTT" + target + "TTTTT" + target + "TTTTT")
    row = {"guideId": "dup", "targetSeq": target, "mitSpecScore": "1", "cfdSpecScore": "1",
           "offtargetCount": "0", "targetGenomeGeneLocus": "x"}
    with pytest.raises(AmbiguousTargetMapping):
        _designer()._row_to_candidate(row, win, "HBG1")
