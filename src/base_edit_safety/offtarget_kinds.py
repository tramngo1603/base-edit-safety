"""Keep different possible editing effects separate.

spacer_homology_cas_style lists CRISPOR reference-sequence candidates, not measurements.
base_editor_context covers editor-specific effects that need experimental testing.
guide_independent_deaminase covers DNA changes without guide targeting.
rna_off_target covers RNA editing.

This module fills only the sequence-search category. It marks the others as unassessed and
lists follow-up assays. Candidate-nomination assays do not establish cellular editing.
"""
from __future__ import annotations

SPACER_HOMOLOGY = "spacer_homology_cas_style"
BASE_EDITOR_CONTEXT = "base_editor_context"
GUIDE_INDEPENDENT_DEAMINASE = "guide_independent_deaminase"
RNA_OFF_TARGET = "rna_off_target"

NOT_ASSESSED = "not_assessed_here"


def classify_offtargets(offtarget_summary: dict | None) -> dict:
    """Group the guide's candidate information by possible effect.

    Fill the sequence-match category from the search. Mark the other categories as unassessed
    and list the required assays.
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
