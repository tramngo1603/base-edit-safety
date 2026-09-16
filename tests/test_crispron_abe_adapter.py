"""CRISPRon-ABE fallback parsing and process isolation."""
from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from base_edit_safety.adapters.crispron_abe import (
    ABE8E_WEIGHTS,
    MODEL_NAME,
    PINNED_REVISION,
    CrispronAbeFallbackPredictor,
)
from base_edit_safety.config import CrispronAbeConfig
from base_edit_safety.pipeline.datatypes import (
    MODEL_PREDICTED_FALLBACK,
    Guide,
    spectrum_to_records,
)

CONTEXT = "GTCAAGGCAAGGCTGGCCAACCCATGGGTG"
PROTO = "AGGCAAGGCTGGCCAACCCA"
OUTCOMES = (
    "seq_id\ttarget\toutcome\tpred_eff\tpred_freq\n"
    "context_0_p_5\tGTCAAGGCAAGGCTGGCCAACCCATGGGTG\t"
    "GTCAAGGCGGGGCTGGCCAACCCATGGGTG\t59.67570877075195\t58.80533745550009\n"
    "context_0_p_5\tGTCAAGGCAAGGCTGGCCAACCCATGGGTG\t"
    "GTCAAGGCGAGGCTGGCCAACCCATGGGTG\t59.67570877075195\t0.5438450868339061\n"
    "context_0_p_5\tGTCAAGGCAAGGCTGGCCAACCCATGGGTG\t"
    "GTCAAGGCAGGGCTGGCCAACCCATGGGTG\t59.67570877075195\t0.3265262284179564\n"
)


def _guide(gene: str) -> Guide:
    hbg1 = gene == "HBG1"
    return Guide(
        guide_id=f"{gene}_60forw_ABE8e",
        protospacer=PROTO,
        pam="TGG",
        strand="+",
        genome_build="GRCh38",
        chrom="chr11",
        gene=gene,
        gene_strand="-",
        tss_1based=5_249_857 if hbg1 else 5_254_781,
        editor_class="ABE",
        editor_name="ABE8e",
        edit_from="A",
        edit_to="G",
        proto_5p_genomic=5_249_976 if hbg1 else 5_254_900,
        target_context_30mer=CONTEXT,
    )


def _config(tmp_path: Path) -> CrispronAbeConfig:
    repo = tmp_path / "crispron-BE"
    (repo / "bin").mkdir(parents=True)
    (repo / "bin" / "CRISPRonBE.sh").write_text("#!/bin/bash\n")
    python = tmp_path / "env" / "bin" / "python3"
    python.parent.mkdir(parents=True)
    python.write_text("")
    return CrispronAbeConfig(str(repo), str(python))


def test_duplicate_locus_context_is_executed_once_and_mapped_to_both_guides(tmp_path):
    calls = []

    def fake_run(command, **kwargs):
        calls.append((command, kwargs))
        fasta = Path(command[3]).read_text()
        assert fasta == f">context_0\n{CONTEXT}\n"
        output = Path(command[4])
        output.mkdir()
        (output / "crispronABE_prediction.tsv").write_text(OUTCOMES)
        return subprocess.CompletedProcess(command, 0, stdout="#All done\n", stderr="")

    predictor = CrispronAbeFallbackPredictor(
        config=_config(tmp_path),
        run_command=fake_run,
        revision_resolver=lambda _: PINNED_REVISION,
    )
    spectra = predictor.predict_many([_guide("HBG1"), _guide("HBG2")])

    assert len(calls) == 1
    command, kwargs = calls[0]
    assert command[0] == "/bin/bash"
    assert command[2] == "ABE"
    assert command[5] == ABE8E_WEIGHTS
    assert kwargs["cwd"].startswith("/")
    assert kwargs["env"]["PATH"].startswith(str(tmp_path / "env" / "bin"))
    assert set(spectra) == {"HBG1_60forw_ABE8e", "HBG2_60forw_ABE8e"}
    for spectrum in spectra.values():
        assert spectrum.model_name == MODEL_NAME
        assert sum(haplotype.pred_score for haplotype in spectrum.haplotypes) == pytest.approx(1.0)
        assert spectrum.per_position_rate()[5] == pytest.approx(0.5934918254233399)
        assert spectrum.per_position_rate()[6] == pytest.approx(0.5913186368391804)
        assert spectrum.provenance["predicted_editing_efficiency_rate"] == \
            pytest.approx(0.5967570877075195)
        assert spectrum.provenance["model_relative"] is True
        assert spectrum.provenance["deduplicated_identical_context"] is True
        assert spectrum.provenance["shared_input_guide_ids"] == [
            "HBG1_60forw_ABE8e", "HBG2_60forw_ABE8e"
        ]

    hbg1 = spectra["HBG1_60forw_ABE8e"]
    intended = {hbg1.guide.genomic_coord(5), hbg1.guide.genomic_coord(6)}
    records = spectrum_to_records(hbg1, on_target_coords=intended)
    on_target = [record for record in records if record.on_target_or_bystander == "on_target"]
    assert [record.predicted_rate for record in on_target] == [0.593492, 0.591319]
    assert all(record.confidence == MODEL_PREDICTED_FALLBACK for record in on_target)
    assert records[0].provenance["joint_on_target_rate"] == 0.588053
    assert records[0].provenance["joint_on_target_confidence"] == MODEL_PREDICTED_FALLBACK


def test_shell_success_without_prediction_file_is_rejected(tmp_path):
    def fake_run(command, **kwargs):
        del kwargs
        return subprocess.CompletedProcess(command, 0, stdout="failed", stderr="traceback")

    predictor = CrispronAbeFallbackPredictor(
        config=_config(tmp_path),
        run_command=fake_run,
        revision_resolver=lambda _: PINNED_REVISION,
    )

    with pytest.raises(RuntimeError, match="no accepted fallback output"):
        predictor.predict_many([_guide("HBG1")])


def test_unpinned_checkout_is_rejected_before_execution(tmp_path):
    with pytest.raises(RuntimeError, match="must be pinned"):
        CrispronAbeFallbackPredictor(
            config=_config(tmp_path),
            revision_resolver=lambda _: "not-the-pinned-revision",
        )


def test_context_must_match_guide_spacer_and_pam(tmp_path):
    guide = _guide("HBG1")
    guide = Guide(**{**guide.__dict__, "target_context_30mer": "A" * 30})
    predictor = CrispronAbeFallbackPredictor(
        config=_config(tmp_path),
        revision_resolver=lambda _: PINNED_REVISION,
    )

    with pytest.raises(ValueError, match="does not match spacer and PAM"):
        predictor.predict_many([guide])
