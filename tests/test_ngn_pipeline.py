"""Non-NGG designs keep structural analysis while rate prediction stays out of scope."""
from __future__ import annotations

from base_edit_safety import anchors
from base_edit_safety.adapters.bedict import BedictPredictionBatch
from base_edit_safety.design_record import build_design_records
from base_edit_safety.io.outputs import records_to_dicts
from base_edit_safety.pipeline import enumerate_guides
from base_edit_safety.pipeline.datatypes import (
    EDIT_WINDOW,
    MODEL_PREDICTED,
    MODEL_PREDICTED_FALLBACK,
    NON_NGG_MODEL_SCOPE_REASON,
    NOT_EVALUATED,
    CandidateGuide,
    EditSpectrum,
    OutcomeHaplotype,
)

PROTO = "GCTAACGAACGTACGATCAA"


class _Designer:
    name = "CRISPOR"

    def __init__(self, pam_class: str, observed_pam: str):
        self.pam = pam_class
        self.observed_pam = observed_pam
        self.last_ambiguous = []

    def design(self, chrom: str, start: int, end: int, gene: str):
        del start, end
        tss = anchors.HBG1_TSS if gene == "HBG1" else anchors.HBG2_TSS
        self.last_ambiguous = []
        return [
            CandidateGuide(
                guide_id=f"{gene}_fixture",
                protospacer=PROTO,
                pam=self.observed_pam,
                strand="-",
                genome_build="GRCh38",
                chrom=chrom,
                gene=gene,
                gene_strand="-",
                tss_1based=tss,
                proto_5p_genomic=tss + 120,
                offtarget={
                    "source": "CRISPOR",
                    "mit_specificity": "45",
                    "cfd_specificity": "80",
                    "offtarget_count": "2",
                    "mm_histogram": {"2": 1},
                    "cfd_weighted_burden": 0.2,
                    "paralog_count": 1,
                    "offtarget_count_excl_paralog": 1,
                },
                offtarget_sites=[
                    {
                        "chrom": "chr2",
                        "start": "25264110",
                        "strand": "+",
                        "mismatchCount": "2",
                        "cfdOfftargetScore": "0.2",
                        "locusDesc": "intron:DNMT3A",
                    }
                ],
                pam_class=self.pam,
                target_context_30mer="TTTT" + PROTO + self.observed_pam + "AAA",
            )
        ]

    def tool_metadata(self):
        return {"name": self.name, "pam": self.pam}


class _Predictor:
    name = "BE-DICT-bystander"
    runs = 5

    def __init__(self):
        self.calls = 0

    def tool_metadata(self):
        return {"name": self.name, "runs": self.runs}

    def predict_many(self, guides):
        self.calls += 1
        spectra = {}
        for guide in guides:
            editable = [
                position
                for position in EDIT_WINDOW
                if guide.protospacer[position - 1] == guide.edit_from
            ]
            position = 8 if guide.editor_class == "ABE" else editable[0]
            changed = list(guide.protospacer)
            changed[position - 1] = guide.edit_to
            spectra[guide.guide_id] = EditSpectrum(
                guide=guide,
                model_name=f"fixture:{guide.editor_name}",
                haplotypes=[
                    OutcomeHaplotype(guide.protospacer, 0.8),
                    OutcomeHaplotype("".join(changed), 0.2),
                ],
                provenance={
                    "ensemble_runs": self.runs,
                    "normalization_applied": "per_member_then_ensemble_v1",
                    "structural_stability": "stable",
                    "upstream_mass_anomaly": False,
                    "raw_mass": 1.0,
                    "raw_mass_per_member": {
                        f"run_{run_num}": 1.0 for run_num in range(self.runs)
                    },
                    "rate_semantics": "model_relative",
                },
            )
        return BedictPredictionBatch(spectra=spectra, exclusions=[])


class _FallbackPredictor(_Predictor):
    def predict_many(self, guides):
        batch = super().predict_many(guides)
        exclusions = []
        for guide in guides:
            if guide.editor_class != "ABE":
                continue
            del batch.spectra[guide.guide_id]
            exclusions.append({
                "guide_id": guide.guide_id,
                "reason": f"{guide.guide_id}: unstable member rank order; BE-DICT rates excluded",
                "status": "excluded_structurally_unstable_ensemble",
                "structural_stability": "unstable",
                "upstream_mass_anomaly": True,
                "raw_mass": 1.3,
                "raw_mass_per_member": {
                    f"run_{run_num}": 1.3 for run_num in range(self.runs)
                },
                "fallback": {
                    "predictor": "CRISPRon-ABE",
                    "analysis_role": "sensitivity_fallback",
                    "status": "routed_not_executed",
                },
            })
        batch.exclusions.extend(exclusions)
        return batch


class _ExecutedCrispronFallback:
    name = "CRISPRon-ABE"

    def __init__(self):
        self.calls = 0

    def tool_metadata(self):
        return {"name": self.name, "model_version": "fixture"}

    def predict_many(self, guides):
        self.calls += 1
        spectra = {}
        for guide in guides:
            changed = list(guide.protospacer)
            changed[7] = "G"
            spectra[guide.guide_id] = EditSpectrum(
                guide=guide,
                model_name="CRISPRon-ABE:fixture:ABE8e",
                haplotypes=[
                    OutcomeHaplotype(guide.protospacer, 0.6),
                    OutcomeHaplotype("".join(changed), 0.4),
                ],
                provenance={
                    "ensemble_runs": 1,
                    "prediction_confidence": MODEL_PREDICTED_FALLBACK,
                    "predictor": "CRISPRon-ABE",
                    "predictor_role": "fallback_for_bedict_instability",
                    "model_version": "CRISPRonBE_models-1.0",
                    "model_member_count": 25,
                    "model_weights": "0-0-0-0-1",
                    "model_relative": True,
                    "rate_semantics": "model_relative",
                    "rate_transform": "official_percent_to_fraction_v1",
                    "target_context_30mer": guide.target_context_30mer,
                    "predicted_editing_efficiency_percent": 40.0,
                    "predicted_editing_efficiency_rate": 0.4,
                    "outcome_spectrum": [{
                        "outcome_protospacer": "".join(changed),
                        "edited_protospacer_positions": [8],
                        "predicted_frequency_percent": 40.0,
                        "model_relative_rate": 0.4,
                    }],
                    "deduplicated_identical_context": False,
                    "shared_input_guide_ids": [guide.guide_id],
                },
            )
        return spectra


def _designs(records, editor_name):
    record_set = build_design_records(records_to_dicts(records))
    matches = [
        record
        for record in record_set["records"]
        if record["candidate"]["editor_name"] == editor_name
    ]
    return record_set, matches


def test_ngn_produces_complete_unscored_records_without_constructing_bedict(monkeypatch):
    designer = _Designer("NGN", "TGT")
    monkeypatch.setattr(enumerate_guides, "CrisporGuideDesigner", lambda: designer)

    def forbidden_predictor():
        raise AssertionError("BE-DICT must not be constructed for a non-NGG design")

    monkeypatch.setattr(enumerate_guides, "BedictBystanderPredictor", forbidden_predictor)

    records = enumerate_guides.run_enumeration(write=False)
    abe_records = [record for record in records if record.editor_class == "ABE"]

    assert abe_records
    assert {record.gene for record in abe_records} == {"HBG1", "HBG2"}
    assert any(record.on_target_or_bystander == "on_target" for record in abe_records)
    assert any(record.on_target_or_bystander == "bystander" for record in abe_records)
    assert all(record.predicted_rate is None for record in abe_records)
    assert all(record.confidence == NOT_EVALUATED for record in abe_records)
    assert all(record.provenance["pam_class"] == "NGN" for record in abe_records)
    assert all(record.provenance["reason"] == NON_NGG_MODEL_SCOPE_REASON
               for record in abe_records)
    assert all(record.provenance["co_targets_paralog"] for record in abe_records)
    assert all(record.provenance["offtarget"]["mm_histogram"] == {"2": 1}
               for record in abe_records)
    assert all(record.provenance["ch_driver_offtargets"][0]["ch_driver_gene"] == "DNMT3A"
               for record in abe_records)

    record_set, designs = _designs(records, "ABE8e")

    assert len(designs) == 2
    for design in designs:
        assert design["provenance"]["pam_class"] == "NGN"
        assert design["provenance"]["reason"] == NON_NGG_MODEL_SCOPE_REASON
        spectrum = design["predicted_edit_spectrum"]
        assert all(edit["predicted_rate"] is None for edit in spectrum["on_target"])
        assert all(edit["confidence"] == NOT_EVALUATED for edit in spectrum["on_target"])
        assert all(edit["predicted_rate"] is None for edit in spectrum["bystanders"])
        assert all(edit["confidence"] == NOT_EVALUATED for edit in spectrum["bystanders"])
        assert design["registry_match"]["relationship"] == "exact_edit"
        assert design["hbg_paralog_flags"]["flagged"] is True
        assert design["ch_annotation"]["driver_offtargets"][0]["ch_driver_gene"] == "DNMT3A"
    assert "non-NGG designs" in " ".join(record_set["limitations"])


def test_ngg_path_still_calls_bedict_and_keeps_model_rates(monkeypatch):
    designer = _Designer("NGG", "AGG")
    predictor = _Predictor()
    monkeypatch.setattr(enumerate_guides, "CrisporGuideDesigner", lambda: designer)
    monkeypatch.setattr(enumerate_guides, "BedictBystanderPredictor", lambda: predictor)

    records = enumerate_guides.run_enumeration(write=False)
    on_target = [
        record
        for record in records
        if record.editor_class == "ABE" and record.on_target_or_bystander == "on_target"
    ]

    assert predictor.calls == 2
    assert len(on_target) == 2
    assert all(record.predicted_rate == 0.2 for record in on_target)
    assert all(record.confidence == MODEL_PREDICTED for record in on_target)
    assert all(record.provenance["pam_class"] == "NGG" for record in on_target)
    assert all("reason" not in record.provenance for record in on_target)

    record_set, designs = _designs(records, "ABE8e")

    assert len(designs) == 2
    assert all(design["provenance"]["pam_class"] == "NGG" for design in designs)
    assert all("reason" not in design["provenance"] for design in designs)
    assert all(design["predicted_edit_spectrum"]["on_target"][0]["predicted_rate"] == 0.2
               for design in designs)
    assert "non-NGG designs" not in " ".join(record_set["limitations"])


def test_unstable_ngg_abe_is_retained_as_marked_crispron_fallback(monkeypatch):
    designer = _Designer("NGG", "AGG")
    predictor = _FallbackPredictor()
    monkeypatch.setattr(enumerate_guides, "CrisporGuideDesigner", lambda: designer)
    monkeypatch.setattr(enumerate_guides, "BedictBystanderPredictor", lambda: predictor)
    monkeypatch.setattr(enumerate_guides, "crispron_abe_configured", lambda: False)

    records = enumerate_guides.run_enumeration(write=False)
    fallback_records = [record for record in records if record.editor_class == "ABE"]

    assert fallback_records
    assert all(record.predicted_rate is None for record in fallback_records)
    assert all(record.confidence == NOT_EVALUATED for record in fallback_records)
    assert all(record.model_name is None for record in fallback_records)
    assert all(record.provenance["structural_stability"] == "unstable"
               for record in fallback_records)
    assert all(record.provenance["analysis_role"] == "sensitivity_fallback"
               for record in fallback_records)
    assert all(record.provenance["fallback_predictor"] == "CRISPRon-ABE"
               for record in fallback_records)
    assert all(record.provenance["fallback_status"] == "routed_not_executed"
               for record in fallback_records)
    assert all(record.provenance["raw_mass_per_member"] == {
        f"run_{run_num}": 1.3 for run_num in range(predictor.runs)
    } for record in fallback_records)

    record_set, designs = _designs(records, "ABE8e")
    assert len(designs) == 2
    assert all(design["provenance"]["analysis_role"] == "sensitivity_fallback"
               for design in designs)
    assert all(design["predicted_edit_spectrum"]["on_target"][0]["predicted_rate"] is None
               for design in designs)
    assert "CRISPRon sensitivity/fallback" in " ".join(record_set["limitations"])


def test_executed_crispron_fallback_repopulates_rates_and_keeps_bedict_exclusion(monkeypatch):
    designer = _Designer("NGG", "AGG")
    predictor = _FallbackPredictor()
    fallback = _ExecutedCrispronFallback()
    monkeypatch.setattr(enumerate_guides, "CrisporGuideDesigner", lambda: designer)
    monkeypatch.setattr(enumerate_guides, "BedictBystanderPredictor", lambda: predictor)
    monkeypatch.setattr(enumerate_guides, "crispron_abe_configured", lambda: True)
    monkeypatch.setattr(
        enumerate_guides, "CrispronAbeFallbackPredictor", lambda: fallback
    )
    written = {}

    def capture_write(records, filename, meta):
        written.update({"records": records, "filename": filename, "meta": meta})

    monkeypatch.setattr(enumerate_guides.outputs, "write_records", capture_write)

    records = enumerate_guides.run_enumeration(write=True)
    on_target = [
        record for record in records
        if record.editor_class == "ABE" and record.on_target_or_bystander == "on_target"
    ]

    assert fallback.calls == 2
    assert len(on_target) == 2
    assert all(record.predicted_rate == 0.4 for record in on_target)
    assert all(record.confidence == MODEL_PREDICTED_FALLBACK for record in on_target)
    assert all(record.model_name == "CRISPRon-ABE:fixture:ABE8e" for record in on_target)
    assert all(record.provenance["predictor"] == "CRISPRon-ABE" for record in on_target)
    assert all(record.provenance["predictor_role"] ==
               "fallback_for_bedict_instability" for record in on_target)
    assert all(record.provenance["model_relative"] is True for record in on_target)
    assert all(record.provenance["fallback_status"] == "executed" for record in on_target)
    assert all(record.provenance["fallback_executed_date"] == "2026-08-26"
               for record in on_target)
    assert all(record.provenance["structural_stability"] == "unstable"
               for record in on_target)
    assert all(record.provenance["raw_mass"] == 1.3 for record in on_target)
    assert all(record.provenance["primary_predictor_status"] ==
               "excluded_structurally_unstable_ensemble" for record in on_target)
    assert all(record.provenance["primary_predictor_exclusion_link"] == {
        "ledger": "meta.prediction_exclusions",
        "guide_id": record.provenance["guide_id"],
        "status": "excluded_structurally_unstable_ensemble",
    } for record in on_target)
    assert all(record.provenance["rejected_bedict_scores_used"] is False
               for record in on_target)
    assert all(record.provenance["fallback_blended_with_primary"] is False
               for record in on_target)

    record_set, designs = _designs(records, "ABE8e")
    assert len(designs) == 2
    assert all(design["prediction_tier"] == "predicted_only" for design in designs)
    assert all(design["evidence_tier"] != "empirically_measured" for design in designs)
    assert all(design["provenance"]["fallback_status"] == "executed"
               for design in designs)
    assert all(design["predicted_edit_spectrum"]["on_target"][0]["confidence"] ==
               MODEL_PREDICTED_FALLBACK for design in designs)
    assert "neither revived nor blended" in " ".join(record_set["limitations"])

    exclusions = written["meta"]["prediction_exclusions"]
    assert len(exclusions) == 2
    assert all(exclusion["fallback"]["status"] == "routed_not_executed"
               for exclusion in exclusions)
    routes = written["meta"]["fallback_routes"]
    assert len(routes) == 2
    assert all(route["status"] == "executed" for route in routes)
    assert all(route["original_route_status"] == "routed_not_executed" for route in routes)
    assert written["meta"]["fallback_execution_count"] == 2
