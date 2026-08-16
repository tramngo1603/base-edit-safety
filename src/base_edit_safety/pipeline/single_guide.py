"""Run the prediction path for one guide as an integration check and example.

It builds one NGG-PAM guide that puts a targeted position inside the editing window, predicts what
the editor does to it, turns that into per-position EditRecords, leaves the consequence field for
later, and writes the result to data/outputs. It needs the external predictor (see config.py).

Use `enumerate_guides.py` to evaluate every guide across both promoter windows.
"""
from __future__ import annotations

from .. import anchors
from ..adapters.bedict import BedictBystanderPredictor
from ..adapters.consequence import mark_consequence_deferred
from ..io import outputs
from ..io.reference import fetch_window
from . import environment
from .datatypes import (
    EDIT_SOURCE_BASE, EDIT_TARGET_BASE, EDIT_WINDOW, EditRecord, Guide,
    reverse_complement, spectrum_to_records,
)

PROTO_LEN = 20


def _gene_tss(gene: str) -> int:
    return anchors.HBG1_TSS if gene == "HBG1" else anchors.HBG2_TSS


def build_anchor_guide(gene: str, promoter_offset: int, editor_class: str,
                       editor_name: str | None = None) -> Guide:
    """Find a guide whose protospacer puts `promoter_offset` inside the editing window."""
    tss = _gene_tss(gene)
    target = tss + promoter_offset
    edit_from, edit_to = EDIT_SOURCE_BASE[editor_class], EDIT_TARGET_BASE[editor_class]
    # Grab enough sequence to fit a 20 nt protospacer plus a 3 nt PAM on either side of the target.
    win = fetch_window("chr11", target - 30, target + 30)
    plus = win.plus_sequence
    lo = win.start_1based

    def plus_sub(a: int, b: int) -> str:  # inclusive genomic [a, b]
        return plus[a - lo:b - lo + 1]

    candidates = []
    # Plus strand: the protospacer starts at the low coordinate, and the PAM follows at higher ones.
    for g in range(lo, win.end_1based - PROTO_LEN - 2):
        proto = plus_sub(g, g + PROTO_LEN - 1)
        pam = plus_sub(g + PROTO_LEN, g + PROTO_LEN + 2)
        if len(proto) != PROTO_LEN or pam[1:] != "GG":
            continue
        tpos = target - g + 1
        if tpos in EDIT_WINDOW and proto[tpos - 1] == edit_from:
            candidates.append(("+", proto, pam, g, tpos))
    # Minus strand: the protospacer starts at the high coordinate, and the PAM follows at lower ones.
    for g in range(lo + PROTO_LEN + 2, win.end_1based + 1):
        proto = reverse_complement(plus_sub(g - PROTO_LEN + 1, g))
        pam = reverse_complement(plus_sub(g - PROTO_LEN - 2, g - PROTO_LEN))
        if len(proto) != PROTO_LEN or pam[1:] != "GG":
            continue
        tpos = g - target + 1
        if tpos in EDIT_WINDOW and proto[tpos - 1] == edit_from:
            candidates.append(("-", proto, pam, g, tpos))

    if not candidates:
        raise ValueError(
            f"No NGG-PAM guide puts {gene} offset {promoter_offset} inside the editing window."
        )
    # Prefer the guide that puts the target nearest the middle of the window (position 6).
    strand, proto, pam, proto_5p, tpos = min(candidates, key=lambda c: abs(c[4] - 6))
    guide_id = f"{gene}_{'m' if promoter_offset >= 0 else 'p'}{abs(promoter_offset)}_{editor_class}"
    return Guide(
        guide_id=guide_id,
        protospacer=proto,
        pam=pam,
        strand=strand,
        genome_build=anchors.GENOME_BUILD,
        chrom=anchors.CHROM,
        gene=gene,
        gene_strand=anchors.GENE_STRAND,
        tss_1based=tss,
        editor_class=editor_class,
        editor_name=editor_name or "",
        edit_from=edit_from,
        edit_to=edit_to,
        proto_5p_genomic=proto_5p,
    )


def run_single_guide(
    gene: str = "HBG1",
    promoter_offset: int = 113,
    editor_class: str = "ABE",
    editor_name: str = "ABE8e",
    write: bool = True,
) -> list[EditRecord]:
    """Run one guide end to end and return its EditRecords."""
    guide = build_anchor_guide(gene, promoter_offset, editor_class, editor_name)
    predictor = BedictBystanderPredictor()
    spectrum = predictor.predict(guide)
    records = spectrum_to_records(spectrum, on_target_offset=promoter_offset)
    records = mark_consequence_deferred(records)
    if write:
        repro = environment.reproducibility(
            tools=[predictor.tool_metadata()],
            editors={editor_class: editor_name},
            ensemble_runs=predictor.runs,
            reference_windows={},
        )
        outputs.write_records(
            records,
            filename=f"single_guide_{guide.guide_id}.json",
            meta={
                "step": "1 (one guide, end to end)",
                "gene": gene,
                "anchor_offset": -abs(promoter_offset),
                "guide_id": guide.guide_id,
                "protospacer": guide.protospacer,
                "pam": guide.pam,
                "strand": guide.strand,
                "editor": editor_name,
                "predictor": spectrum.model_name,
                "reproducibility": repro,
            },
        )
    return records


if __name__ == "__main__":
    recs = run_single_guide()
    for r in recs:
        role = r.on_target_or_bystander
        print(
            f"{r.gene} offset {r.promoter_offset:>4} ({r.chrom}:{r.pos_1based}) "
            f"genomic {r.genomic_plus_ref}>{r.genomic_plus_alt} / sense "
            f"{r.promoter_sense_ref}>{r.promoter_sense_alt}  rate={r.predicted_rate}  [{role}]"
        )
