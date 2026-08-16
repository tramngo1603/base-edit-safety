"""Record or predict an edit's effect on gene regulation.

`mark_consequence_deferred` records that regulatory effects were not scored in the current run.
`AlphaGenomeConsequenceAnnotator` calls the optional external service and stores its predicted
scores. `pipeline/annotate_consequence.py` runs that optional step.
"""
from __future__ import annotations

import json
import os

from ..pipeline.datatypes import EditRecord

DEFERRED_VALUE = "deferred"
_DEFERRED_SOURCE = "deferred-consequence-passthrough"


def mark_consequence_deferred(records: list[EditRecord]) -> list[EditRecord]:
    """Mark regulatory effects as not scored for this run; update records in place."""
    for record in records:
        record.consequence = DEFERRED_VALUE
        record.provenance = {
            **record.provenance,
            "consequence_source": _DEFERRED_SOURCE,
            "consequence_note": "consequence not scored in this run; annotate_consequence.py fills it in",
        }
    return records


# --- Optional AlphaGenome annotator ---------------------------------------------------------------
# AlphaGenome is non-commercial under its own terms, and those terms say its outputs must not be used
# to train a model. The external package must be installed separately; this adapter calls its client
# without copying source code or weights. The API key is read from ALPHAGENOME_API_KEY and is not
# logged, printed, cached, or included in output.
#
# Severity is the mean absolute `quantile_score` across the service's recommended output tracks. It
# describes the predicted size of a regulatory change, not its direction or clinical effect.

DEFAULT_FUNCTIONAL_THRESHOLD = 0.5
_SEQUENCE_LENGTH = "SEQUENCE_LENGTH_16KB"
_CACHE_PATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))),
    "data", "inputs", "alphagenome_cache.json",
)


class AlphaGenomeConsequenceAnnotator:
    """Scores how much one variant disturbs gene regulation, through the AlphaGenome client."""

    name = "AlphaGenome"

    def __init__(self, functional_threshold: float = DEFAULT_FUNCTIONAL_THRESHOLD, use_cache: bool = True):
        self.functional_threshold = functional_threshold  # provisional; lock before validation
        self.use_cache = use_cache
        self._model = None
        self._cache = self._load_cache() if use_cache else {}

    # --- External client and cache ----------------------------------------------
    def _ensure_model(self):
        if self._model is None:
            from alphagenome.models import dna_client  # delayed import for an optional dependency
            key = os.environ.get("ALPHAGENOME_API_KEY")
            if not key:
                raise RuntimeError("ALPHAGENOME_API_KEY is not set, so the model cannot be reached.")
            self._model = dna_client.create(key)  # the key is passed to the client and not persisted here
        return self._model

    def _load_cache(self) -> dict:
        try:
            with open(_CACHE_PATH) as handle:
                return json.load(handle)
        except (FileNotFoundError, json.JSONDecodeError):
            return {}

    def _save_cache(self):
        if not self.use_cache:
            return
        os.makedirs(os.path.dirname(_CACHE_PATH), exist_ok=True)
        with open(_CACHE_PATH, "w") as handle:
            json.dump(self._cache, handle)  # cache contains scores but no credential

    # --- scoring ----------------------------------------------------------------
    def score(self, chrom: str, pos_1based: int, ref: str, alt: str) -> dict:
        """Return the predicted regulatory-effect score and class for one variant."""
        vkey = f"{chrom}:{pos_1based}:{ref}>{alt}"
        if self.use_cache and vkey in self._cache:
            return self._cache[vkey]

        import numpy as np
        from alphagenome.models import dna_client, variant_scorers
        from alphagenome.data import genome

        model = self._ensure_model()
        variant = genome.Variant(chromosome=chrom, position=pos_1based,
                                 reference_bases=ref, alternate_bases=alt)
        interval = variant.reference_interval.resize(getattr(dna_client, _SEQUENCE_LENGTH))
        scores = model.score_variant(interval=interval, variant=variant, variant_scorers=None)
        df = variant_scorers.tidy_scores(scores)
        q = np.abs(df["quantile_score"].astype(float).to_numpy())
        q = q[~np.isnan(q)]
        severity = float(q.mean()) if q.size else None
        top = df.iloc[int(np.nanargmax(np.abs(df["quantile_score"].astype(float))))]
        result = {
            "source": self.name,
            "license": "non-commercial; outputs must not be used to train models",
            "model_predicted": True,
            "measured": False,
            "severity": round(severity, 4) if severity is not None else None,
            "severity_metric": "mean |quantile_score| over the recommended tracks",
            "classification": self._classify(severity),
            "top_output_type": str(top.get("output_type")),
            "top_biosample": str(top.get("biosample_name")),
            "n_tracks": int(len(df)),
        }
        if self.use_cache:
            self._cache[vkey] = result
            self._save_cache()
        return result

    def _classify(self, severity) -> str:
        if severity is None:
            return "unknown"
        return "functionally-consequential" if severity >= self.functional_threshold else "low-consequence"

    def annotate(self, record: EditRecord) -> EditRecord:
        """Fill in one record's consequence, scoring the plus-strand allele."""
        result = self.score(record.chrom, record.pos_1based,
                            record.genomic_plus_ref, record.genomic_plus_alt)
        record.consequence = result
        record.provenance = {
            **record.provenance,
            "consequence_source": self.name,
            "consequence_note": "regulatory consequence predicted by AlphaGenome",
        }
        return record
