"""Calls BE-DICT to predict what a base editor does to a guide, including its bystander edits.

BE-DICT is MIT-licensed and must be installed in a separate Python environment. Its dependency
versions conflict with this package's versions, so the adapter calls it as a separate process. This
repository does not copy BE-DICT source code or weights. Configure its location with the environment
variables listed in `config.py`.
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

# Normalize BE-DICT scores for every included guide. A total outside these inclusive bounds marks a
# malformed response that is excluded and recorded. These bounds check score format; they are not
# efficacy or ranking thresholds.
NORMALIZATION_APPLIED = "uniform_v1"
RAW_MASS_MIN = 0.75
RAW_MASS_MAX = 1.25


class MalformedHaplotypeMass(ValueError):
    """Raised when a guide's total raw BE-DICT score is outside the allowed range."""

    def __init__(self, guide: Guide, raw_mass: float):
        self.guide = guide
        self.raw_mass = raw_mass
        super().__init__(
            f"{guide.guide_id}: raw haplotype score mass {raw_mass!r} outside "
            f"[{RAW_MASS_MIN}, {RAW_MASS_MAX}]; guide excluded as malformed"
        )

    def exclusion_record(self, editor: str) -> dict:
        """Return a machine-readable record for a guide excluded from a batch."""
        return {
            "guide_id": self.guide.guide_id,
            "gene": self.guide.gene,
            "editor_class": self.guide.editor_class,
            "editor": editor,
            "protospacer": self.guide.protospacer,
            "raw_mass": self.raw_mass,
            "allowed_raw_mass": {"minimum": RAW_MASS_MIN, "maximum": RAW_MASS_MAX},
            "normalization_applied": None,
            "status": "excluded_malformed_raw_mass",
            "reason": str(self),
        }


@dataclass(frozen=True)
class BedictPredictionBatch:
    """Store valid edit spectra and malformed-guide exclusion records from one call."""

    spectra: dict[str, EditSpectrum]
    exclusions: list[dict]


def _validated_haplotypes(
    guide: Guide, raw: list[dict]
) -> tuple[list[OutcomeHaplotype], float]:
    """Validate scores and divide each guide's scores by their total.

    The return value includes the input total and the scaled complete-sequence scores. A total
    outside the allowed range raises a specific error, which the batch caller records as an
    exclusion.
    """
    n = len(guide.protospacer)
    raw_haplotypes: list[OutcomeHaplotype] = []
    for o in raw:
        seq = o.get("output_seq")
        try:
            score = float(o["pred_score"])
        except (KeyError, TypeError, ValueError):
            raise ValueError(f"{guide.guide_id}: non-numeric haplotype probability {o.get('pred_score')!r}")
        if not isinstance(seq, str) or len(seq) != n:
            raise ValueError(f"{guide.guide_id}: haplotype length {len(seq) if isinstance(seq, str) else '?'} "
                             f"!= protospacer length {n}")
        if not math.isfinite(score) or score < 0.0:
            raise ValueError(
                f"{guide.guide_id}: raw haplotype score {score!r} must be finite and non-negative"
            )
        raw_haplotypes.append(OutcomeHaplotype(output_seq=seq, pred_score=score))

    raw_mass = math.fsum(hap.pred_score for hap in raw_haplotypes)
    if raw_mass < RAW_MASS_MIN or raw_mass > RAW_MASS_MAX:
        raise MalformedHaplotypeMass(guide, raw_mass)

    # Apply this division when `raw_mass` equals 1.0 as well, so every guide follows the same rule.
    haplotypes = [
        OutcomeHaplotype(output_seq=hap.output_seq, pred_score=hap.pred_score / raw_mass)
        for hap in raw_haplotypes
    ]
    return haplotypes, raw_mass


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
        """Predict outcomes for a batch of guides, using one predictor call per editor.

        Guides are grouped by editor because each editor uses a separate predictor. The return value
        always includes excluded guides, so callers can distinguish an exclusion from a missing
        result.
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
                except MalformedHaplotypeMass as exc:
                    exclusions.append(exc.exclusion_record(editor))
        return BedictPredictionBatch(spectra=spectra, exclusions=exclusions)

    def spectrum_from_response(self, guide: Guide, editor: str, response: dict) -> EditSpectrum:
        """Turn one predictor response into a checked EditSpectrum. No subprocess, so tests can use it."""
        raw = response.get("sequences", {}).get(guide.guide_id, [])
        haplotypes, raw_mass = _validated_haplotypes(guide, raw)
        return EditSpectrum(
            guide=guide,
            model_name=f"{self.name}:{editor}",
            haplotypes=haplotypes,
            provenance={
                "predictor": self.name,
                "editor": editor,
                "ensemble_runs": response.get("ensemble_runs", self.runs),
                "normalization_applied": NORMALIZATION_APPLIED,
                "raw_mass": raw_mass,
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
