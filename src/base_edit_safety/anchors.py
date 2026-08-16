"""The HBG1/HBG2 promoter coordinates and the edit positions this tool covers.

The table below is copied from `docs/coordinates.md`. This module does not calculate coordinates; it
stores the reviewed table and converts its rows into `EditRecord` objects.

Every position stores two versions of the same allele. Papers write HBG promoter alleles on the
gene's own strand ("-198 T>C"); genome browsers and VCFs write them on the plus strand ("-198 A>G").
HBG1 and HBG2 sit on the minus strand, so the two forms differ. Storing both prevents a strand
orientation from being inferred later.

Loading the table first checks that matching HBG1 and HBG2 positions stay 4,924 bp apart. If they
don't, a row has shifted and the run stops.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Optional

from .pipeline.datatypes import EditRecord

# A-gamma denotes HBG1 and G-gamma denotes HBG2. Papers use the Greek names while this module
# uses HBG gene symbols, so citation matching depends on this explicit mapping.
#
# authority = HGNC nomenclature: HBG1 = "hemoglobin subunit gamma 1" (Aγ);
# HBG2 = "hemoglobin subunit gamma 2" (Gγ).
# Sources: https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:4831 and
# https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:4832
# access_level=database_record; audit_status=NOT_REVIEWED; accessed=2026-08-12. The runtime status
# below remains conservative because the cited records have not completed source review.
GAMMA_GENE_NAMING = {
    "Agamma": "HBG1",
    "Ggamma": "HBG2",
    "status": "HGNC database records listed in source comments; source review pending",
    "consequence_if_wrong": "HBG1 and HBG2 database identifiers and citations would be swapped",
}

GENOME_BUILD = "GRCh38"
CHROM = "chr11"
GENE_STRAND = "-"  # HBG1 and HBG2 are on the minus strand
HBG1_TSS = 5_249_857
HBG2_TSS = 5_254_781
PARALOG_SPACING_BP = 4924  # HBG2_TSS - HBG1_TSS, on both GRCh38 and GRCh37

# The default records that motif loss has not been reviewed; it does not claim that no motif is lost.
LOST_MOTIF_NOT_REVIEWED = "not_reviewed"


@dataclass(frozen=True)
class AnchorPosition:
    """Store one targeted promoter position for both HBG genes.

    The targeted HBG1 and HBG2 promoter sequences use the same reference and alternate bases. Each
    gene keeps its own variant-database status and rsID.
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
    # The motif the edit destroys. Reviewed examples can create one motif while disrupting another.
    # The default therefore means "not reviewed" rather than "no disrupted motif."
    lost_motif: str = LOST_MOTIF_NOT_REVIEWED
    hbg1_rsid: Optional[str] = None
    hbg2_rsid: Optional[str] = None

    @property
    def hbg1_pos(self) -> int:
        return HBG1_TSS + self.promoter_offset

    @property
    def hbg2_pos(self) -> int:
        return HBG2_TSS + self.promoter_offset


# How well a variant database supports each position, per `docs/coordinates.md`.
# db_confirmed   a ClinVar accession AND an rsID, both for this gene.
# db_unverified  no ClinVar record for this gene. Two positions have a dbSNP rsID; a dbSNP
#                submission is not a reviewed ClinVar record and does not promote
#                the status.
# synthetic      a screen hit, not a natural HPFH allele.
#
# rs1261339811 and rs34153902 are dbSNP submissions for the stated gene and coordinate. An rsID for
# one HBG gene is not assigned to the other HBG gene.
#
# All five `db_unverified` rows were searched by position on 2026-07-26. SPDI uses zero-based
# coordinates and HGVS uses one-based coordinates, so the position search checked for a one-base
# conversion error. It returned the same two rsIDs. HBG1 -113, HBG2 -113, and HBG2 -117 returned no
# dbSNP record in that search. Repeat the search only when updating the recorded access date.
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
    """Store a documented position excluded from prediction and the exclusion reason."""

    promoter_offset: int
    genomic_plus_ref: str
    genomic_plus_alt: str
    promoter_sense_ref: str
    promoter_sense_alt: str
    reason: str
    note: str = ""


# Keyed by (offset, ref, alt), because -110 exists twice: A>C here as an NF-Y
# transversion, and A>G in the bystander table as an ABE edit. Offset alone would collide.
OUT_OF_SCOPE: list[ExcludedPosition] = [
    ExcludedPosition(110, "T", "G", "A", "C", "transversion",
                     "NF-Y motif, but base editors do not make transversions. Not the same edit as "
                     "the -110 A>G ABE bystander, which is a different allele at the same offset."),
    ExcludedPosition(114, "G", "T", "C", "A", "transversion",
                     "Doerfler's cellular allele at this position. Not the repo's C>T."),
    ExcludedPosition(114, "G", "C", "C", "G", "transversion",
                     "Doerfler's EMSA oligo allele. Not the repo's C>T."),
    # Ravi lists these as newly reported substitutions. Recording the exclusion makes the window
    # boundary explicit; changing that boundary requires a separate review.
    ExcludedPosition(69, "G", "A", "C", "T", "outside_modelling_window",
                     "Substitution first reported by Ravi 2022. Falls outside the -210 to -80 window."),
    ExcludedPosition(70, "G", "A", "C", "T", "outside_modelling_window",
                     "Substitution first reported by Ravi 2022. Falls outside the -210 to -80 window."),
]


@dataclass(frozen=True)
class OutsideWindowConversion:
    """A conversion reported in the literature that this tool structurally cannot predict."""

    promoter_offset: int
    promoter_sense_ref: str
    promoter_sense_alt: str
    editor_class: str
    source_pmid: str
    note: str


# BE-DICT receives only positions inside the editing window, so it cannot predict a conversion
# outside that window. These observations remain separate from predicted within-window bystanders.
OBSERVED_OUTSIDE_WINDOW: list[OutsideWindowConversion] = [
    OutsideWindowConversion(
        117, "G", "A", "CBE", "35147495",
        "Ravi 2022 reports an observed -117 G>A CBE conversion outside the protospacer and "
        "outside the editing window. Distinct from the -117 CBE anchor above, which is the intended "
        "within-window edit for a different guide. No prediction will ever contain this."),
]


def check_paralog_spacing() -> None:
    """Stop the run if matching HBG1/HBG2 positions are not exactly 4,924 bp apart.

    This catches a shifted row, including an earlier -114/-113 off-by-one error. It does not identify
    the gene copy, so HBG1 or HBG2 must not be assigned by adding or subtracting 4,924.
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
            # An rsID records a dbSNP submission, not a reviewed ClinVar assertion. Preserve that
            # distinction so the identifier is not treated as confirmation.
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
