"""Optionally fill the regulatory-effect field using AlphaGenome.

The step reads the enumeration file, scores each distinct edit, and writes scores in fields that were
previously `deferred`. When configured for them, it also scores CH-driver off-target locations for
review. Scores are cached and duplicate requests are removed. The API key is read from the
environment and is not saved.

Coding off-target sites require a coding-variant method rather than this regulatory method. The
current promoter workflow contains no such site, so that path is outside this module.
"""
from __future__ import annotations

import json

from ..adapters.consequence import AlphaGenomeConsequenceAnnotator
from ..io import outputs

_ENUMERATION_FILE = outputs.enumeration_path()


def run_consequence(annotator: AlphaGenomeConsequenceAnnotator | None = None,
                    write: bool = True, score_ch_offtargets: bool = True) -> dict:
    """Score every record's consequence, and each CH off-target site, in the enumeration file."""
    annotator = annotator or AlphaGenomeConsequenceAnnotator()
    doc = outputs.read_enumeration()   # refuses a file written by a different build
    records = doc["records"]

    # 1) Regulatory consequence per edit, using the plus-strand allele.
    n_scored = 0
    for r in records:
        r["consequence"] = annotator.score(
            r["chrom"], r["pos_1based"], r["genomic_plus_ref"], r["genomic_plus_alt"])
        r["provenance"] = {
            **r["provenance"],
            "consequence_source": annotator.name,
            "consequence_note": "regulatory consequence predicted by AlphaGenome",
        }
        n_scored += 1

    # 2) Score CH-driver off-target locations with a stand-in edit. A sequence match does not specify
    # which base an editor would change, so this step scores a fixed A>G change. The result describes
    # regulatory sensitivity at the location, not the effect of a specific off-target edit. The
    # `is_proxy` field marks this limit. Set `score_ch_offtargets=False` to skip the step.
    n_ch = 0
    if score_ch_offtargets:
        for r in records:
            for hit in r["provenance"].get("ch_driver_offtargets", []):
                if hit.get("ch_functional_consequence"):
                    continue
                cons = annotator.score(hit["chrom"], int(hit["start"]), "A", "G")
                hit["ch_functional_consequence"] = {
                    **cons,
                    "is_proxy": True,
                    "proxy_edit": "A>G",
                    "note": "a stand-in A>G edit at the off-target site because the specific "
                            "off-target base change is unknown. Regulatory severity only.",
                }
                n_ch += 1

    doc["meta"] = {**doc.get("meta", {}),
                   "consequence_layer": "filled in by AlphaGenome (non-commercial; outputs must not "
                                        "be used to train models)",
                   "consequence_scored_edits": n_scored,
                   "ch_offtarget_consequence_scored": n_ch}
    if write:
        with open(_ENUMERATION_FILE, "w") as handle:
            json.dump(doc, handle, indent=2)
    return {"scored_edits": n_scored, "ch_offtargets_scored": n_ch}


if __name__ == "__main__":
    print(run_consequence())
