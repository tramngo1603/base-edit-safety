"""Store the HBG1/HBG2 coordinates and base changes used by this tool.

The values match docs/coordinates.md; this module does not look them up again. Each change
stores both promoter-strand and genomic-plus alleles because HBG1 and HBG2 are minus-strand
genes. Loading the table checks that matching positions remain 4,924 bp apart.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Optional

from .pipeline.datatypes import EditRecord

# Agamma means HBG1; Ggamma means HBG2. Use this mapping when connecting paper names to coordinates.
#
# authority = HGNC nomenclature: HBG1 = "hemoglobin subunit gamma 1" (Aγ);
# HBG2 = "hemoglobin subunit gamma 2" (Gγ).
# Sources: https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:4831 and
# https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:4832
# access_level=database_record; audit_status=NOT_REVIEWED; accessed=2026-08-12.
# The runtime status remains conservative until these source records complete review.
GAMMA_GENE_NAMING = {
    "Agamma": "HBG1",
    "Ggamma": "HBG2",
    "status": "project convention, no authority cited",
    "consequence_if_wrong": "every per-gene rsID, ClinVar record and citation attribution inverts",
}

GENOME_BUILD = "GRCh38"
CHROM = "chr11"
GENE_STRAND = "-"  # HBG1 and HBG2 are on the minus strand
HBG1_TSS = 5_249_857
HBG2_TSS = 5_254_781
PARALOG_SPACING_BP = 4924  # HBG2_TSS - HBG1_TSS, on both GRCh38 and GRCh37

# Default for lost_motif: nobody has checked, which is not the same as "there isn't one".
LOST_MOTIF_NOT_REVIEWED = "not_reviewed"


@dataclass(frozen=True)
class AnchorPosition:
    """One promoter position this tool targets, shared by both paralogs.

    HBG1 and HBG2 use the same alleles because their proximal promoters are near-identical. What
    differs per gene is whether a variant database lists the position, and under which rsID, so
    those two are stored separately for each gene.
    """

    promoter_offset: int                 # distance upstream of the cap site, as a positive number
    editor_class: Literal["ABE", "CBE"]
    genomic_plus_ref: str
    genomic_plus_alt: str
    promoter_sense_ref: str
    promoter_sense_alt: str
    mechanism: str                       # the motif the edit CREATES
    hbg1_status: str
    hbg2_status: str
    # Record any motif lost by the edit separately from a motif it creates.
    # The default means unreviewed, not that no motif is lost.
    lost_motif: str = LOST_MOTIF_NOT_REVIEWED
    hbg1_rsid: Optional[str] = None
    hbg2_rsid: Optional[str] = None

    @property
    def hbg1_pos(self) -> int:
        return HBG1_TSS + self.promoter_offset

    @property
    def hbg2_pos(self) -> int:
        return HBG2_TSS + self.promoter_offset


# Variant database status for each gene, from the recorded source review.
# db_confirmed: a ClinVar accession and an rsID for this gene.
# db_unverified: no matching ClinVar record in that review; a dbSNP rsID may still exist.
# synthetic: a constructed screen edit rather than a catalogued natural HPFH allele.
# An rsID is never copied from one paralog to the other.
#
# The position-based check on 2026-07-26 found rs1261339811 for HBG2 -198 and rs34153902
# for HBG1 -175. It found no matching dbSNP record for HBG1 -113, HBG2 -113 or HBG2 -117.
# These are dated search results; database records can change.
DB = "db_confirmed"
UNVERIFIED = "db_unverified"
SYNTHETIC = "synthetic"

IN_SCOPE: list[AnchorPosition] = [
    AnchorPosition(198, "ABE", "A", "G", "T", "C", "de novo KLF1 motif",
                   hbg1_status=DB, hbg2_status=UNVERIFIED,
                   hbg1_rsid="rs35710727", hbg2_rsid="rs1261339811"),
    AnchorPosition(175, "ABE", "A", "G", "T", "C", "TAL1 E-box (looping)",
                   # Confirmed second consequence, from Wienert's supplement: the same substitution
                   # completes the E-box (sense -175 to -170) and destroys an overlapping GATA1 core
                   # (sense -175 to -172) at the same base.
                   lost_motif="overlapping GATA1 core, sense -175 to -172",
                   hbg1_status=UNVERIFIED, hbg2_status=DB,
                   hbg1_rsid="rs34153902", hbg2_rsid="rs63750654"),
    AnchorPosition(124, "ABE", "A", "G", "T", "C", "synthetic screen hit; KLF1 (double-mutant probe)",
                   hbg1_status=SYNTHETIC, hbg2_status=SYNTHETIC),
    AnchorPosition(123, "ABE", "A", "G", "T", "C", "synthetic screen hit; KLF1 (double-mutant probe)",
                   hbg1_status=SYNTHETIC, hbg2_status=SYNTHETIC),
    AnchorPosition(117, "CBE", "C", "T", "G", "A", "distal CCAAT / BCL11A motif (tiled/motif-level)",
                   hbg1_status=DB, hbg2_status=UNVERIFIED, hbg1_rsid="rs35378915"),
    AnchorPosition(114, "CBE", "G", "A", "C", "T", "distal BCL11A motif (tiled/motif-level)",
                   hbg1_status=DB, hbg2_status=DB,
                   hbg1_rsid="rs281860601", hbg2_rsid="rs34809449"),
    # Neither gene has a -113 record of its own, so both stay db_unverified.
    AnchorPosition(113, "ABE", "T", "C", "A", "G", "de novo GATA1 motif",
                   # Confirmed second consequence, also from a supplement.
                   lost_motif="disrupts NF-Y",
                   hbg1_status=UNVERIFIED, hbg2_status=UNVERIFIED),
    # Ravi 2022 states in both Results and Discussion that a CBE converts -122 G>A as the ON-TARGET
    # site for gRNA-10/-11, inside the editing window. It was in no table at all: not an anchor, not
    # a bystander. Strand check passes, -122 is a G on the sense strand, so the CBE acts on the
    # antisense strand and reads out as sense G>A.
    AnchorPosition(122, "CBE", "C", "T", "G", "A",
                   "synthetic screen hit (Ravi 2022); CBE on-target for gRNA-10/-11",
                   hbg1_status=SYNTHETIC, hbg2_status=SYNTHETIC),
]


@dataclass(frozen=True)
class ExcludedPosition:
    """A position we know about and deliberately do not model, with the reason recorded."""

    promoter_offset: int
    genomic_plus_ref: str
    genomic_plus_alt: str
    promoter_sense_ref: str
    promoter_sense_alt: str
    reason: str
    note: str = ""


# Keyed by (offset, ref, alt) like everything else, because -110 exists twice: A>C here as an NF-Y
# transversion, and A>G in the bystander table as an ABE edit. Offset alone would collide.
OUT_OF_SCOPE: list[ExcludedPosition] = [
    ExcludedPosition(110, "T", "G", "A", "C", "transversion",
                     "NF-Y motif, but base editors do not make transversions. Not the same edit as "
                     "the -110 A>G ABE bystander, which is a different allele at the same offset."),
    ExcludedPosition(114, "G", "T", "C", "A", "transversion",
                     "Doerfler's cellular allele at this position. Not the repo's C>T."),
    ExcludedPosition(114, "G", "C", "C", "G", "transversion",
                     "Doerfler's EMSA oligo allele. Not the repo's C>T."),
    # Ravi lists these among that study's novel substitutions. Recorded rather than omitted so the
    # omission is a decision instead of an oversight; widening the window is a separate call.
    ExcludedPosition(69, "G", "A", "C", "T", "outside_modelling_window",
                     "Ravi 2022 novel substitution. Falls outside the -210 to -80 window."),
    ExcludedPosition(70, "G", "A", "C", "T", "outside_modelling_window",
                     "Ravi 2022 novel substitution. Falls outside the -210 to -80 window."),
]


@dataclass(frozen=True)
class OutsideWindowConversion:
    """A reported conversion outside this tool's supported prediction window."""

    promoter_offset: int
    promoter_sense_ref: str
    promoter_sense_alt: str
    editor_class: str
    source_pmid: str
    note: str


# No BE-DICT rates are reported outside the supported window. Keep reported outside-window
# conversions separate from predicted bystanders so this limit remains visible.
OBSERVED_OUTSIDE_WINDOW: list[OutsideWindowConversion] = [
    OutsideWindowConversion(
        117, "G", "A", "CBE", "35147495",
        "Ravi 2022 reports a real -117 G>A CBE conversion occurring outside the protospacer AND "
        "outside the editing window. Distinct from the -117 CBE anchor above, which is the intended "
        "within-window edit for a different guide. No prediction will ever contain this."),
]


def check_paralog_spacing() -> None:
    """Stop the run if matching HBG1/HBG2 positions are not exactly 4,924 bp apart.

    This catches a shifted row (it caught a -114/-113 off-by-one once). It does not tell you which
    gene a position belongs to, so never work a gene out by adding or subtracting 4,924.
    """
    if HBG2_TSS - HBG1_TSS != PARALOG_SPACING_BP:
        raise AssertionError(
            f"HBG1/HBG2 TSS are {HBG2_TSS - HBG1_TSS} bp apart, expected {PARALOG_SPACING_BP}."
        )
    for pos in IN_SCOPE:
        if pos.hbg2_pos - pos.hbg1_pos != PARALOG_SPACING_BP:
            raise AssertionError(
                f"offset {pos.promoter_offset}: HBG1/HBG2 positions are "
                f"{pos.hbg2_pos - pos.hbg1_pos} bp apart, expected {PARALOG_SPACING_BP}."
            )


def _record_for_gene(pos: AnchorPosition, gene: str) -> EditRecord:
    if gene == "HBG1":
        coord, status, rsid = pos.hbg1_pos, pos.hbg1_status, pos.hbg1_rsid
    elif gene == "HBG2":
        coord, status, rsid = pos.hbg2_pos, pos.hbg2_status, pos.hbg2_rsid
    else:
        raise ValueError(gene)
    provenance = {
        "source": "docs/coordinates.md",
        "promoter_offset_numbering": "counted back from the cap site: +1 is the first transcribed "
                                     "base, and there is no zero",
        "mechanism": pos.mechanism,
    }
    if rsid:
        provenance["dbsnp"] = rsid
        if status == UNVERIFIED:
            # An rsID here means somebody submitted the variant, not that anyone curated it. Saying
            # so on the record stops a reader treating the id as confirmation.
            provenance["dbsnp_note"] = "dbSNP submission only; no ClinVar record for this gene"
    return EditRecord(
        genome_build=GENOME_BUILD,
        chrom=CHROM,
        pos_1based=coord,
        genomic_plus_ref=pos.genomic_plus_ref,
        genomic_plus_alt=pos.genomic_plus_alt,
        promoter_sense_ref=pos.promoter_sense_ref,
        promoter_sense_alt=pos.promoter_sense_alt,
        gene=gene,
        gene_strand=GENE_STRAND,
        promoter_offset=pos.promoter_offset,
        editor_class=pos.editor_class,
        on_target_or_bystander="on_target",
        variant_db_status=status,
        provenance=provenance,
    )


def load_in_scope_records() -> list[EditRecord]:
    """Return the table as EditRecords, one per gene per position. Checks the spacing first."""
    check_paralog_spacing()
    records: list[EditRecord] = []
    for pos in IN_SCOPE:
        records.append(_record_for_gene(pos, "HBG1"))
        records.append(_record_for_gene(pos, "HBG2"))
    return records
