"""Run a separate BE-DICT installation to predict edited sequences and bystanders.

BE-DICT uses its own Python environment. This adapter starts a separate process and copies
neither source nor model weights into this repository. Set the paths listed in config.py.
"""
from __future__ import annotations

from dataclasses import dataclass
import json
import math
import os
import subprocess
import tempfile

from ..config import BedictConfig, load_bedict_config
from ..pipeline.datatypes import EDIT_WINDOW, EditSpectrum, Guide, OutcomeHaplotype

_RUNNER = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_bedict_runner.py")

# Which editor to assume when a guide does not name one.
DEFAULT_EDITOR = {"ABE": "ABE8e", "CBE": "BE4max"}

# Normalize each model's scores before averaging, then normalize the average for rounding error.
# Do this only after the models pass the agreement checks. The scores are relative model outputs,
# not calibrated probabilities of editing in cells.
NORMALIZATION_APPLIED = "per_member_then_ensemble_v1"
RATE_SEMANTICS = "model_relative"

# The old range marks unusual raw score sums. Exclusion uses the wider range below;
# a sum inside that range does not by itself exclude the guide.
SUPERSEDED_RAW_MASS_MIN = 0.75
SUPERSEDED_RAW_MASS_MAX = 1.25
RAW_MASS_MIN = 0.5
RAW_MASS_MAX = 2.0

STRUCTURAL_STABILITY_STABLE = "stable"
STRUCTURAL_STABILITY_UNSTABLE = "unstable"
FALLBACK_ANALYSIS_ROLE = "sensitivity_fallback"
FALLBACK_STATUS = "routed_not_executed"


def _fallback_route(guide: Guide) -> dict:
    editor = "ABE" if guide.editor_class == "ABE" else "CBE"
    return {
        "predictor": f"CRISPRon-{editor}",
        "analysis_role": FALLBACK_ANALYSIS_ROLE,
        "status": FALLBACK_STATUS,
    }


class MalformedHaplotypeMass(ValueError):
    """Raised when raw BE-DICT scores sum outside the accepted range."""

    def __init__(
        self,
        guide: Guide,
        raw_mass: float,
        raw_mass_per_member: dict[str, float] | None = None,
        detail: str | None = None,
    ):
        self.guide = guide
        self.raw_mass = raw_mass
        self.raw_mass_per_member = raw_mass_per_member or {}
        self.detail = detail or (
            f"raw haplotype score mass {raw_mass!r} outside "
            f"[{RAW_MASS_MIN}, {RAW_MASS_MAX}] after stability check; guide excluded"
        )
        super().__init__(f"{guide.guide_id}: {self.detail}")

    def exclusion_record(self, editor: str) -> dict:
        """Return the excluded guide, its raw scores and the reason."""
        return {
            "guide_id": self.guide.guide_id,
            "gene": self.guide.gene,
            "editor_class": self.guide.editor_class,
            "editor": editor,
            "protospacer": self.guide.protospacer,
            "raw_mass": self.raw_mass,
            "raw_mass_per_member": self.raw_mass_per_member,
            "allowed_raw_mass": {"minimum": RAW_MASS_MIN, "maximum": RAW_MASS_MAX},
            "structural_stability": STRUCTURAL_STABILITY_STABLE,
            "upstream_mass_anomaly": True,
            "rate_semantics": None,
            "normalization_applied": None,
            "status": "excluded_raw_mass_outside_authorized_range",
            "reason": str(self),
        }


class StructurallyUnstableEnsemble(ValueError):
    """Raised when model members are missing, inconsistent, or disagree on outcome order."""

    def __init__(
        self,
        guide: Guide,
        raw_mass: float,
        raw_mass_per_member: dict[str, float],
        member_rank_order: dict[str, list[str]],
        detail: str = "ensemble members do not have the same score-ranked haplotype order",
    ):
        self.guide = guide
        self.raw_mass = raw_mass
        self.raw_mass_per_member = raw_mass_per_member
        self.member_rank_order = member_rank_order
        self.detail = detail
        super().__init__(f"{guide.guide_id}: {detail}; BE-DICT rates excluded")

    def exclusion_record(self, editor: str) -> dict:
        """Return the exclusion and the recorded fallback choice."""
        anomaly = _has_upstream_mass_anomaly(
            [self.raw_mass, *self.raw_mass_per_member.values()]
        )
        return {
            "guide_id": self.guide.guide_id,
            "gene": self.guide.gene,
            "editor_class": self.guide.editor_class,
            "editor": editor,
            "protospacer": self.guide.protospacer,
            "raw_mass": self.raw_mass,
            "raw_mass_per_member": self.raw_mass_per_member,
            "allowed_raw_mass": {"minimum": RAW_MASS_MIN, "maximum": RAW_MASS_MAX},
            "structural_stability": STRUCTURAL_STABILITY_UNSTABLE,
            "member_rank_order": self.member_rank_order,
            "upstream_mass_anomaly": anomaly,
            "rate_semantics": None,
            "normalization_applied": None,
            "status": "excluded_structurally_unstable_ensemble",
            "fallback": _fallback_route(self.guide),
            "reason": str(self),
        }


@dataclass(frozen=True)
class _ValidatedEnsemble:
    haplotypes: list[OutcomeHaplotype]
    raw_mass: float
    raw_mass_per_member: dict[str, float]
    upstream_mass_anomaly: bool


def _has_upstream_mass_anomaly(masses: list[float]) -> bool:
    return any(
        mass < SUPERSEDED_RAW_MASS_MIN or mass > SUPERSEDED_RAW_MASS_MAX
        for mass in masses
    )


def _parse_raw_haplotypes(guide: Guide, raw: list[dict]) -> list[OutcomeHaplotype]:
    """Validate one raw member or aggregate without changing its scores."""
    if not isinstance(raw, list):
        raise ValueError(f"{guide.guide_id}: haplotypes must be a list")
    n = len(guide.protospacer)
    parsed: list[OutcomeHaplotype] = []
    seen_sequences: set[str] = set()
    for outcome in raw:
        seq = outcome.get("output_seq")
        try:
            score = float(outcome["pred_score"])
        except (KeyError, TypeError, ValueError):
            raise ValueError(
                f"{guide.guide_id}: non-numeric haplotype score {outcome.get('pred_score')!r}"
            )
        if not isinstance(seq, str) or len(seq) != n:
            raise ValueError(
                f"{guide.guide_id}: haplotype length "
                f"{len(seq) if isinstance(seq, str) else '?'} != protospacer length {n}"
            )
        if seq in seen_sequences:
            raise ValueError(f"{guide.guide_id}: duplicate haplotype sequence {seq!r}")
        if not math.isfinite(score) or score < 0.0 or score > 1.0:
            raise ValueError(
                f"{guide.guide_id}: raw haplotype score {score!r} must be finite and inside [0, 1]"
            )
        seen_sequences.add(seq)
        parsed.append(OutcomeHaplotype(output_seq=seq, pred_score=score))
    return parsed


@dataclass(frozen=True)
class BedictPredictionBatch:
    """Predictions plus the reasons for every guide excluded from the batch."""

    spectra: dict[str, EditSpectrum]
    exclusions: list[dict]


def _validated_haplotypes(
    guide: Guide, raw: list[dict]
) -> tuple[list[OutcomeHaplotype], float]:
    """Check the raw score sum, then normalize the outcome scores.

    Return the raw sum too, so later records show what was changed. For several model members,
    use _validated_ensemble_haplotypes to check their agreement before normalization.
    """
    raw_haplotypes = _parse_raw_haplotypes(guide, raw)
    raw_mass = math.fsum(hap.pred_score for hap in raw_haplotypes)
    if raw_mass < RAW_MASS_MIN or raw_mass > RAW_MASS_MAX:
        raise MalformedHaplotypeMass(guide, raw_mass)

    haplotypes = [
        OutcomeHaplotype(output_seq=hap.output_seq, pred_score=hap.pred_score / raw_mass)
        for hap in raw_haplotypes
    ]
    if haplotypes:
        final_mass = math.fsum(hap.pred_score for hap in haplotypes)
        last = haplotypes[-1]
        haplotypes[-1] = OutcomeHaplotype(
            output_seq=last.output_seq,
            pred_score=last.pred_score + (1.0 - final_mass),
        )
    return haplotypes, raw_mass


def _rank_order(haplotypes: list[OutcomeHaplotype]) -> list[str]:
    """Sequence order under descending score, with a deterministic tie break."""
    return [
        hap.output_seq
        for hap in sorted(haplotypes, key=lambda hap: (-hap.pred_score, hap.output_seq))
    ]


def _validated_ensemble_haplotypes(
    guide: Guide,
    raw_aggregate: list[dict],
    member_payloads: list[dict] | None,
    ensemble_runs: int,
) -> _ValidatedEnsemble:
    """Check member structure, normalize members, average them, and normalize the final scores."""
    aggregate = _parse_raw_haplotypes(guide, raw_aggregate)
    raw_mass = math.fsum(hap.pred_score for hap in aggregate)

    if member_payloads is None and ensemble_runs == 1:
        member_payloads = [{"member": "run_0", "haplotypes": raw_aggregate}]
    if not isinstance(member_payloads, list) or len(member_payloads) != ensemble_runs:
        raise StructurallyUnstableEnsemble(
            guide,
            raw_mass=raw_mass,
            raw_mass_per_member={},
            member_rank_order={},
            detail=(
                f"expected {ensemble_runs} member prediction sets, received "
                f"{len(member_payloads) if isinstance(member_payloads, list) else 0}"
            ),
        )

    scores_per_member: dict[str, dict[str, float]] = {}
    raw_mass_per_member: dict[str, float] = {}
    member_rank_order: dict[str, list[str]] = {}
    for index, payload in enumerate(member_payloads):
        if not isinstance(payload, dict):
            raise ValueError(f"{guide.guide_id}: ensemble member {index} must be an object")
        member = payload.get("member")
        if not isinstance(member, str) or not member or member in scores_per_member:
            raise StructurallyUnstableEnsemble(
                guide,
                raw_mass=raw_mass,
                raw_mass_per_member=raw_mass_per_member,
                member_rank_order=member_rank_order,
                detail="ensemble member labels are missing or duplicated",
            )
        parsed = _parse_raw_haplotypes(guide, payload.get("haplotypes"))
        scores_per_member[member] = {hap.output_seq: hap.pred_score for hap in parsed}
        raw_mass_per_member[member] = math.fsum(hap.pred_score for hap in parsed)
        member_rank_order[member] = _rank_order(parsed)

    expected_members = {f"run_{run_num}" for run_num in range(ensemble_runs)}
    if set(scores_per_member) != expected_members:
        raise StructurallyUnstableEnsemble(
            guide,
            raw_mass=raw_mass,
            raw_mass_per_member=raw_mass_per_member,
            member_rank_order=member_rank_order,
            detail="ensemble member labels do not match the requested production runs",
        )

    first_member = f"run_0"
    expected_sequences = set(scores_per_member[first_member])
    expected_order = member_rank_order[first_member]
    structurally_stable = all(
        set(scores) == expected_sequences and member_rank_order[member] == expected_order
        for member, scores in scores_per_member.items()
    )
    if set(hap.output_seq for hap in aggregate) != expected_sequences:
        structurally_stable = False
    if not structurally_stable:
        raise StructurallyUnstableEnsemble(
            guide,
            raw_mass=raw_mass,
            raw_mass_per_member=raw_mass_per_member,
            member_rank_order=member_rank_order,
        )

    outside_members = {
        member: mass
        for member, mass in raw_mass_per_member.items()
        if mass < RAW_MASS_MIN or mass > RAW_MASS_MAX
    }
    if raw_mass < RAW_MASS_MIN or raw_mass > RAW_MASS_MAX or outside_members:
        detail = (
            f"aggregate raw mass {raw_mass!r} or member masses {outside_members!r} outside "
            f"[{RAW_MASS_MIN}, {RAW_MASS_MAX}] after stability check; guide excluded"
        )
        raise MalformedHaplotypeMass(
            guide,
            raw_mass,
            raw_mass_per_member=raw_mass_per_member,
            detail=detail,
        )

    averaged_scores = []
    for sequence in expected_order:
        normalized_member_scores = [
            scores_per_member[f"run_{run_num}"][sequence]
            / raw_mass_per_member[f"run_{run_num}"]
            for run_num in range(ensemble_runs)
        ]
        averaged_scores.append(math.fsum(normalized_member_scores) / ensemble_runs)
    final_mass = math.fsum(averaged_scores)
    final_scores = [score / final_mass for score in averaged_scores]
    if final_scores:
        final_scores[-1] += 1.0 - math.fsum(final_scores)

    haplotypes = [
        OutcomeHaplotype(output_seq=sequence, pred_score=score)
        for sequence, score in zip(expected_order, final_scores)
    ]
    return _ValidatedEnsemble(
        haplotypes=haplotypes,
        raw_mass=raw_mass,
        raw_mass_per_member=raw_mass_per_member,
        upstream_mass_anomaly=_has_upstream_mass_anomaly(
            [raw_mass, *raw_mass_per_member.values()]
        ),
    )


class BedictBystanderPredictor:
    """Predicts a guide's edit outcomes by averaging over BE-DICT's ensemble of models."""

    name = "BE-DICT-bystander"

    def __init__(self, config: BedictConfig | None = None, runs: int = 5):
        self.config = config or load_bedict_config()
        self.runs = runs

    def predict(self, guide: Guide) -> EditSpectrum:
        """Predict outcomes for one guide."""
        batch = self.predict_many([guide])
        if batch.exclusions:
            raise ValueError(batch.exclusions[0]["reason"])
        return batch.spectra[guide.guide_id]

    def predict_many(self, guides: list[Guide]) -> BedictPredictionBatch:
        """Predict guides in one external call per editor.

        Each editor has its own models. Return excluded guides with their reasons alongside
        successful predictions, so a missing prediction cannot disappear from the result.
        """
        by_editor: dict[str, list[Guide]] = {}
        for guide in guides:
            editor = guide.editor_name or DEFAULT_EDITOR[guide.editor_class]
            by_editor.setdefault(editor, []).append(guide)

        spectra: dict[str, EditSpectrum] = {}
        exclusions: list[dict] = []
        for editor, group in by_editor.items():
            # Every guide for a given editor changes the same base, so take it from the first one.
            request = {
                "editor": editor,
                "edit_from": group[0].edit_from,
                "edit_to": group[0].edit_to,
                "edit_window": [min(EDIT_WINDOW), max(EDIT_WINDOW)],
                "runs": self.runs,
                "sequences": [{"seq_id": g.guide_id, "inp_seq": g.protospacer} for g in group],
            }
            response = self._run(request)
            for guide in group:
                try:
                    spectra[guide.guide_id] = self.spectrum_from_response(guide, editor, response)
                except (MalformedHaplotypeMass, StructurallyUnstableEnsemble) as exc:
                    exclusions.append(exc.exclusion_record(editor))
        return BedictPredictionBatch(spectra=spectra, exclusions=exclusions)

    def spectrum_from_response(self, guide: Guide, editor: str, response: dict) -> EditSpectrum:
        """Turn one predictor response into a checked EditSpectrum. No subprocess, so tests can use it."""
        raw = response.get("sequences", {}).get(guide.guide_id, [])
        try:
            ensemble_runs = int(response.get("ensemble_runs", self.runs))
        except (TypeError, ValueError):
            raise ValueError(f"{guide.guide_id}: ensemble_runs must be a positive integer")
        if ensemble_runs < 1:
            raise ValueError(f"{guide.guide_id}: ensemble_runs must be a positive integer")
        member_payloads = response.get("member_predictions", {}).get(guide.guide_id)
        validated = _validated_ensemble_haplotypes(
            guide,
            raw_aggregate=raw,
            member_payloads=member_payloads,
            ensemble_runs=ensemble_runs,
        )
        return EditSpectrum(
            guide=guide,
            model_name=f"{self.name}:{editor}",
            haplotypes=validated.haplotypes,
            provenance={
                "predictor": self.name,
                "editor": editor,
                "ensemble_runs": ensemble_runs,
                "normalization_applied": NORMALIZATION_APPLIED,
                "structural_stability": STRUCTURAL_STABILITY_STABLE,
                "upstream_mass_anomaly": validated.upstream_mass_anomaly,
                "raw_mass": validated.raw_mass,
                "raw_mass_per_member": validated.raw_mass_per_member,
                "rate_semantics": RATE_SEMANTICS,
                "license": "MIT",
                "execution": "separate process, in its own environment",
            },
        )

    def tool_metadata(self) -> dict:
        """What this predictor was for this run: paths, run count, and git commit if there is one."""
        from ..pipeline import environment
        return environment.tool_environment(
            self.name, repo_root=self.config.repo_root,
            python_executable=self.config.python_executable, ensemble_runs=self.runs,
        )

    def _run(self, request: dict) -> dict:
        with tempfile.TemporaryDirectory() as tmp:
            req_path = os.path.join(tmp, "request.json")
            resp_path = os.path.join(tmp, "response.json")
            with open(req_path, "w") as handle:
                json.dump(request, handle)
            completed = subprocess.run(
                [self.config.python_executable, _RUNNER, self.config.repo_root, req_path, resp_path],
                capture_output=True,
                text=True,
            )
            if completed.returncode != 0 or not os.path.isfile(resp_path):
                raise RuntimeError(
                    "External bystander predictor failed.\n"
                    f"stderr (tail):\n{completed.stderr[-2000:]}"
                )
            with open(resp_path) as handle:
                return json.load(handle)
