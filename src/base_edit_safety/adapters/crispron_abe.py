"""Run the separately installed CRISPRon-ABE fallback predictor.

CRISPRon-BE is BSL-1.1 software and remains in its own checkout and Conda environment. This
adapter invokes the official shell entry point as a separate process, validates its TSV output,
and maps the reported 30-nt outcomes back to the repository's 20-nt guide spectra.
"""
from __future__ import annotations

import csv
import math
import os
import subprocess
import tempfile
from collections.abc import Callable

from ..config import CrispronAbeConfig, load_crispron_abe_config
from ..pipeline import environment
from ..pipeline.datatypes import (
    EDIT_WINDOW,
    MODEL_PREDICTED_FALLBACK,
    EditSpectrum,
    Guide,
    OutcomeHaplotype,
)

NAME = "CRISPRon-ABE"
PREDICTOR_ROLE = "fallback_for_bedict_instability"
FALLBACK_STATUS_EXECUTED = "executed"
FALLBACK_EXECUTED_DATE = "2026-08-26"
RATE_SEMANTICS = "model_relative"
RATE_TRANSFORM = "official_percent_to_fraction_v1"

CODE_VERSION = "v1.0"
PINNED_REVISION = "c5b28d4362b84c7e933e3715346ffbc1181fd6f2"
MODEL_VERSION = "CRISPRonBE_models-1.0"
MODEL_MEMBER_COUNT = 25
CRISPRON_VERSION = "1.0"
CRISPRON_V0_MODEL_VERSION = "CRISPRon_models_V0-1.0"
CRISPROFF_VERSION = "1.1.2"
ABE8E_WEIGHTS = "0-0-0-0-1"
MODEL_NAME = f"{NAME}:{MODEL_VERSION}:ABE8e"


def _full_git_revision(repo_root: str) -> str | None:
    try:
        completed = subprocess.run(
            ["git", "-C", repo_root, "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            timeout=5,
        )
    except Exception:
        return None
    return completed.stdout.strip() if completed.returncode == 0 else None


def _finite_percent(value: str, field: str, seq_id: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        raise ValueError(f"{seq_id}: CRISPRon-ABE {field} is not numeric: {value!r}")
    if not math.isfinite(number) or not 0.0 <= number <= 100.0:
        raise ValueError(f"{seq_id}: CRISPRon-ABE {field} {number!r} outside [0, 100]")
    return number


def _outcome_positions(target: str, outcome: str, seq_id: str) -> list[int]:
    if len(target) != 30 or len(outcome) != 30:
        raise ValueError(f"{seq_id}: CRISPRon-ABE target and outcome must both be 30 nt")
    positions = []
    for index, (ref, alt) in enumerate(zip(target, outcome)):
        if ref == alt:
            continue
        protospacer_position = index - 4 + 1
        if ref != "A" or alt != "G" or protospacer_position not in EDIT_WINDOW:
            raise ValueError(
                f"{seq_id}: unsupported CRISPRon-ABE change {ref}>{alt} at 30-mer "
                f"position {index + 1}"
            )
        positions.append(protospacer_position)
    if not positions:
        raise ValueError(f"{seq_id}: CRISPRon-ABE returned an unedited row as an edited outcome")
    return positions


def _parse_rows(guide: Guide, context: str, rows: list[dict]) -> EditSpectrum:
    """Check the output table and convert its outcomes to sequences with relative weights."""
    if not rows:
        raise ValueError(f"{guide.guide_id}: CRISPRon-ABE returned no edited outcomes")
    expected_seq_id = rows[0].get("seq_id", "")
    efficiencies = {
        _finite_percent(row.get("pred_eff"), "pred_eff", expected_seq_id) for row in rows
    }
    if len(efficiencies) != 1:
        raise ValueError(f"{guide.guide_id}: inconsistent CRISPRon-ABE pred_eff values")
    efficiency_percent = efficiencies.pop()

    seen_outcomes: set[str] = set()
    edited_haplotypes: list[OutcomeHaplotype] = []
    outcome_spectrum: list[dict] = []
    for row in rows:
        seq_id = row.get("seq_id", "")
        target = str(row.get("target", "")).upper()
        outcome = str(row.get("outcome", "")).upper()
        if target != context:
            raise ValueError(
                f"{guide.guide_id}: CRISPRon-ABE target {target!r} does not match "
                f"the submitted context"
            )
        if outcome in seen_outcomes:
            raise ValueError(f"{guide.guide_id}: duplicate CRISPRon-ABE outcome {outcome}")
        positions = _outcome_positions(target, outcome, seq_id)
        frequency_percent = _finite_percent(row.get("pred_freq"), "pred_freq", seq_id)
        seen_outcomes.add(outcome)
        outcome_protospacer = outcome[4:24]
        edited_haplotypes.append(
            OutcomeHaplotype(outcome_protospacer, frequency_percent / 100.0)
        )
        outcome_spectrum.append({
            "outcome_30mer": outcome,
            "outcome_protospacer": outcome_protospacer,
            "edited_protospacer_positions": positions,
            "predicted_frequency_percent": frequency_percent,
            "model_relative_rate": frequency_percent / 100.0,
        })

    frequency_sum = math.fsum(
        outcome["predicted_frequency_percent"] for outcome in outcome_spectrum
    )
    if not math.isclose(frequency_sum, efficiency_percent, rel_tol=0.0, abs_tol=1e-6):
        raise ValueError(
            f"{guide.guide_id}: edited outcome frequencies sum to {frequency_sum}, "
            f"not pred_eff {efficiency_percent}"
        )
    unedited_rate = 1.0 - efficiency_percent / 100.0
    haplotypes = [OutcomeHaplotype(guide.protospacer, unedited_rate), *edited_haplotypes]
    if not math.isclose(
        math.fsum(haplotype.pred_score for haplotype in haplotypes),
        1.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    ):
        raise ValueError(f"{guide.guide_id}: CRISPRon-ABE fractional spectrum does not sum to 1")

    return EditSpectrum(
        guide=guide,
        model_name=MODEL_NAME,
        haplotypes=haplotypes,
        provenance={
            "ensemble_runs": 1,
            "prediction_confidence": MODEL_PREDICTED_FALLBACK,
            "predictor": NAME,
            "predictor_role": PREDICTOR_ROLE,
            "model_version": MODEL_VERSION,
            "model_member_count": MODEL_MEMBER_COUNT,
            "model_weights": ABE8E_WEIGHTS,
            "model_relative": True,
            "rate_semantics": RATE_SEMANTICS,
            "rate_transform": RATE_TRANSFORM,
            "target_context_30mer": context,
            "predicted_editing_efficiency_percent": efficiency_percent,
            "predicted_editing_efficiency_rate": efficiency_percent / 100.0,
            "unedited_model_relative_rate": unedited_rate,
            "outcome_spectrum": outcome_spectrum,
        },
    )


class CrispronAbeFallbackPredictor:
    """Run CRISPRon-ABE v1.0 for ABE8e guides excluded by BE-DICT model-agreement checks."""

    name = NAME

    def __init__(
        self,
        config: CrispronAbeConfig | None = None,
        run_command: Callable | None = None,
        revision_resolver: Callable[[str], str | None] | None = None,
    ):
        self.config = config or load_crispron_abe_config()
        self._run_command = run_command or subprocess.run
        resolver = revision_resolver or _full_git_revision
        self.revision = resolver(self.config.repo_root)
        if self.revision != PINNED_REVISION:
            raise RuntimeError(
                f"CRISPRon-BE checkout must be pinned to {PINNED_REVISION}; "
                f"found {self.revision or 'unavailable'}"
            )

    def tool_metadata(self) -> dict:
        return environment.tool_environment(
            self.name,
            repo_root=self.config.repo_root,
            python_executable=self.config.python_executable,
            pinned_revision=PINNED_REVISION,
            code_version=CODE_VERSION,
            model_version=MODEL_VERSION,
            model_member_count=MODEL_MEMBER_COUNT,
            model_weights=ABE8E_WEIGHTS,
            crispron_version=CRISPRON_VERSION,
            crispron_v0_model_version=CRISPRON_V0_MODEL_VERSION,
            crisproff_version=CRISPROFF_VERSION,
            rate_semantics=RATE_SEMANTICS,
        )

    @staticmethod
    def _context_for(guide: Guide) -> str:
        context = (guide.target_context_30mer or "").upper()
        if len(context) != 30:
            raise ValueError(
                f"{guide.guide_id}: CRISPRon-ABE needs a 30-nt target context"
            )
        if context[4:24] != guide.protospacer or context[24:27] != guide.pam:
            raise ValueError(
                f"{guide.guide_id}: CRISPRon-ABE context does not match spacer and PAM"
            )
        return context

    def _run_one_context(self, context: str, context_number: int, guide: Guide) -> EditSpectrum:
        context_id = f"context_{context_number}"
        with tempfile.TemporaryDirectory(prefix="base-edit-safety-crispron-abe-") as tmp:
            fasta_path = os.path.join(tmp, "input.fasta")
            output_dir = os.path.join(tmp, "output")
            with open(fasta_path, "w") as handle:
                handle.write(f">{context_id}\n{context}\n")
            entrypoint = os.path.join(self.config.repo_root, "bin", "CRISPRonBE.sh")
            command = [
                "/bin/bash",
                entrypoint,
                "ABE",
                fasta_path,
                output_dir,
                ABE8E_WEIGHTS,
            ]
            process_env = os.environ.copy()
            process_env["PATH"] = (
                os.path.dirname(self.config.python_executable)
                + os.pathsep
                + process_env.get("PATH", "")
            )
            process_env["PYTHONNOUSERSITE"] = "1"
            completed = self._run_command(
                command,
                cwd=tmp,
                env=process_env,
                capture_output=True,
                text=True,
                timeout=900,
            )
            result_path = os.path.join(output_dir, "crispronABE_prediction.tsv")
            if completed.returncode != 0 or not os.path.isfile(result_path):
                raise RuntimeError(
                    f"CRISPRon-ABE failed for {guide.guide_id}; no accepted fallback output.\n"
                    f"stdout (tail):\n{completed.stdout[-2000:]}\n"
                    f"stderr (tail):\n{completed.stderr[-2000:]}"
                )
            with open(result_path) as handle:
                all_rows = list(csv.DictReader(handle, delimiter="\t"))
        expected_seq_id = f"{context_id}_p_5"
        rows = [row for row in all_rows if row.get("seq_id") == expected_seq_id]
        if not rows:
            raise RuntimeError(
                f"CRISPRon-ABE output for {guide.guide_id} did not contain {expected_seq_id}"
            )
        return _parse_rows(guide, context, rows)

    def predict_many(self, guides: list[Guide]) -> dict[str, EditSpectrum]:
        """Predict each unique 30-mer once and return a guide-specific spectrum for every guide."""
        by_context: dict[str, list[Guide]] = {}
        for guide in guides:
            if guide.editor_class != "ABE" or guide.editor_name != "ABE8e":
                raise ValueError(
                    f"{guide.guide_id}: CRISPRon-ABE fallback is restricted to ABE8e guides"
                )
            by_context.setdefault(self._context_for(guide), []).append(guide)

        spectra: dict[str, EditSpectrum] = {}
        for context_number, (context, context_guides) in enumerate(by_context.items()):
            template = self._run_one_context(context, context_number, context_guides[0])
            shared_ids = [guide.guide_id for guide in context_guides]
            for guide in context_guides:
                spectra[guide.guide_id] = EditSpectrum(
                    guide=guide,
                    model_name=template.model_name,
                    haplotypes=list(template.haplotypes),
                    provenance={
                        **template.provenance,
                        "deduplicated_identical_context": len(context_guides) > 1,
                        "shared_input_guide_ids": shared_ids,
                    },
                )
        return spectra
