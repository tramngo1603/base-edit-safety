"""Keep four kinds of off-target activity separate.

CRISPOR checks only DNA sites whose sequence resembles the guide. Those sites are candidates, not
measured edits. The workflow marks three other kinds as `not_assessed_here`: DNA activity specific
to the base editor, DNA activity by the editing enzyme without its guide, and unintended RNA edits.
Each unchecked kind names a laboratory test used to measure it.
"""
from __future__ import annotations

SPACER_HOMOLOGY = "spacer_homology_cas_style"
BASE_EDITOR_CONTEXT = "base_editor_context"
GUIDE_INDEPENDENT_DEAMINASE = "guide_independent_deaminase"
RNA_OFF_TARGET = "rna_off_target"

NOT_ASSESSED = "not_assessed_here"


def classify_offtargets(offtarget_summary: dict | None) -> dict:
    """Sort results by off-target kind and mark the three unchecked kinds."""
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
    """List one required assay for each unchecked off-target category."""
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
