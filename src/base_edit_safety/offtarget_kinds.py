"""The kinds of off-target a base editor can produce, and which ones this repo actually checks.

"Off-target" covers several different things, and a Cas9-style sequence search answers only one of
them. This module keeps them apart. For the kinds nothing here checks, it says so and names the
assay that would, so a sequence search is never mistaken for a measured result.

  spacer_homology_cas_style   sites that look like the protospacer, found by a Cas9-style tool
                              (CRISPOR) on the reference genome. Candidates, not measurements.
  base_editor_context         off-targets specific to the editor and its editing window. Measured
                              with CHANGE-seq-BE, GUIDE-seq2, or similar. Not computed here.
  guide_independent_deaminase deaminase activity that happens without the guide. Not computed here.
  rna_off_target              RNA editing (C>U / A>I) across the transcriptome. Not computed here.
"""
from __future__ import annotations

SPACER_HOMOLOGY = "spacer_homology_cas_style"
BASE_EDITOR_CONTEXT = "base_editor_context"
GUIDE_INDEPENDENT_DEAMINASE = "guide_independent_deaminase"
RNA_OFF_TARGET = "rna_off_target"

NOT_ASSESSED = "not_assessed_here"


def classify_offtargets(offtarget_summary: dict | None) -> dict:
    """Sort a guide's off-target findings into the four kinds.

    Only the sequence-match kind gets filled in, from the Cas9-style tool. The other three are
    marked as not checked, each with the assay that would check it.
    """
    summary = offtarget_summary or {}
    return {
        SPACER_HOMOLOGY: {
            "source": summary.get("source"),
            "status": "reference_genome_sequence_match",
            "mm_histogram": summary.get("mm_histogram", {}),
            "cfd_weighted_burden": summary.get("cfd_weighted_burden"),
            "mit_specificity": summary.get("mit_specificity"),
            "caveat": "sites that look like the protospacer on the reference genome, with no bulges "
                      "and no population variants. Candidates only.",
        },
        BASE_EDITOR_CONTEXT: {"status": NOT_ASSESSED,
                              "assay": "CHANGE-seq-BE or GUIDE-seq2 (base-editor context)"},
        GUIDE_INDEPENDENT_DEAMINASE: {"status": NOT_ASSESSED,
                                      "assay": "a genomic deaminase assay such as Detect-seq, plus a "
                                               "no-guide control"},
        RNA_OFF_TARGET: {"status": NOT_ASSESSED,
                         "assay": "RNA-seq editome (C>U / A>I across the transcriptome)"},
    }


def assay_triggers(kinds: dict) -> list[dict]:
    """List the assays this guide still needs. Every unchecked kind adds one."""
    triggers = []
    for kind in (BASE_EDITOR_CONTEXT, GUIDE_INDEPENDENT_DEAMINASE, RNA_OFF_TARGET):
        entry = kinds.get(kind, {})
        if entry.get("status") == NOT_ASSESSED:
            triggers.append({"validates": kind, "assay": entry["assay"], "status": "required"})
    homology = kinds.get(SPACER_HOMOLOGY, {})
    if homology.get("mm_histogram"):
        triggers.append({"validates": SPACER_HOMOLOGY,
                         "assay": "off-target sequencing (CHANGE-seq-BE / GUIDE-seq2) at the "
                                  "candidate sites", "status": "required"})
    return triggers
