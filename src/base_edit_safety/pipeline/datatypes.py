"""The data types that move through the pipeline.

The flow is: a CandidateGuide (no editor picked yet) becomes a Guide (editor picked), a predictor
turns that into an EditSpectrum, and the spectrum reduces to one EditRecord per edited position.

Every record stores both allele orientations. One reference base is insufficient because papers write HBG
promoter alleles on the gene's own strand, genome files write them on the plus strand, and HBG1 and
HBG2 are minus-strand genes. See `docs/coordinates.md`.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Literal

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def complement(base: str) -> str:
    """Complement of a base, keeping the original case."""
    return base.translate(_COMPLEMENT)


def reverse_complement(seq: str) -> str:
    """Reverse complement of a sequence."""
    return seq.translate(_COMPLEMENT)[::-1]


# The supported BE-DICT edit window, counted from the end of the protospacer furthest from the PAM.
# The adapter passes these bounds to BE-DICT's haplotype processor; positions outside it have no
# prediction at all.
EDIT_WINDOW = range(3, 11)
EDIT_SOURCE_BASE = {"ABE": "A", "CBE": "C"}
EDIT_TARGET_BASE = {"ABE": "G", "CBE": "T"}

# How much a rate can be trusted. A position inside EDIT_WINDOW gets a number from the predictor; one
# outside it gets no number, and says so rather than reporting zero.
MODEL_PREDICTED = "model_predicted"
NOT_EVALUATED = "not_evaluated"

# A design-level joint rate is derived from the predictor's complete haplotypes, never from the
# individual position rates. These strings travel with the value so a later record does not turn
# a model output into an unlabeled measurement.
JOINT_ON_TARGET_METHOD = "sum_haplotypes_containing_all_intended_edits"
JOINT_ON_TARGET_CAVEAT = "model-predicted joint frequency, not measured."


@dataclass
class EditRecord:
    """Represent one predicted edit at one genomic position."""

    genome_build: str
    chrom: str
    pos_1based: int
    genomic_plus_ref: str
    genomic_plus_alt: str
    promoter_sense_ref: str
    promoter_sense_alt: str
    gene: str                      # "HBG1" | "HBG2" | "HBG1/HBG2_AMBIGUOUS"
    gene_strand: str               # "-"
    promoter_offset: int
    editor_class: Literal["ABE", "CBE"]
    on_target_or_bystander: Literal["on_target", "bystander"]
    predicted_rate: Optional[float] = None
    model_name: Optional[str] = None
    confidence: Optional[str] = None   # MODEL_PREDICTED | NOT_EVALUATED
    consequence: Optional[str] = None
    ch_driver_flag: bool = False
    ch_driver_gene: Optional[str] = None
    variant_db_status: str = "unset"   # see anchors.py: db_confirmed | db_unverified | synthetic
    provenance: dict = field(default_factory=dict)
    # There is deliberately no clonal-expansion or long-term-risk field. See annotation.py.


@dataclass(frozen=True)
class Guide:
    """A protospacer placed on the genome, with the editor already chosen.

    Positions inside the protospacer are numbered 1-20 from the end furthest from the PAM, which is
    what edit-outcome models expect. Turning one of those positions into a genomic coordinate needs
    the strand the protospacer sits on, so that is stored here too.
    """

    guide_id: str
    protospacer: str                       # 20 nt, 5'->3' on `strand`
    pam: str                               # immediately 3' of the protospacer on `strand`
    strand: Literal["+", "-"]
    genome_build: str
    chrom: str
    gene: str                              # explicit value; not inferred from HBG1/HBG2 spacing
    gene_strand: Literal["+", "-"]
    tss_1based: int                        # transcription start site of `gene`
    editor_class: Literal["ABE", "CBE"]
    editor_name: str                       # e.g. ABE8e, BE4max
    edit_from: str                         # base the editor changes (A for ABE, C for CBE)
    edit_to: str
    proto_5p_genomic: int = field(default=0)  # genomic coordinate of protospacer position 1

    def genomic_coord(self, proto_pos_1based: int) -> int:
        """Genomic coordinate of a protospacer position."""
        if self.strand == "+":
            return self.proto_5p_genomic + (proto_pos_1based - 1)
        return self.proto_5p_genomic - (proto_pos_1based - 1)

    def promoter_offset(self, proto_pos_1based: int) -> int:
        """How far a protospacer position sits from the gene's cap site."""
        return self.genomic_coord(proto_pos_1based) - self.tss_1based

    def editable_positions(self) -> list[int]:
        """Protospacer positions holding the base this editor changes."""
        return [i + 1 for i, b in enumerate(self.protospacer) if b == self.edit_from]

    def alleles_for(self, proto_pos_1based: int) -> dict[str, str]:
        """Both allele orientations for an edit at a protospacer position.

        The editor changes `edit_from` to `edit_to` on the protospacer strand, which gives the
        plus-strand alleles. Flipping again for a minus-strand gene gives the sense-strand ones.
        """
        if self.strand == "+":
            g_ref, g_alt = self.edit_from, self.edit_to
        else:
            g_ref, g_alt = complement(self.edit_from), complement(self.edit_to)
        if self.gene_strand == "+":
            s_ref, s_alt = g_ref, g_alt
        else:
            s_ref, s_alt = complement(g_ref), complement(g_alt)
        return {
            "genomic_plus_ref": g_ref,
            "genomic_plus_alt": g_alt,
            "promoter_sense_ref": s_ref,
            "promoter_sense_alt": s_alt,
        }


@dataclass(frozen=True)
class CandidateGuide:
    """A guide found by the guide-design tool, before an editor is chosen.

    It stores its genomic location and the tool's off-target results. Call `for_editor()` before
    predicting edited sequences.
    """

    guide_id: str
    protospacer: str
    pam: str
    strand: str
    genome_build: str
    chrom: str
    gene: str
    gene_strand: str
    tss_1based: int
    proto_5p_genomic: int
    offtarget: dict = field(default_factory=dict)
    offtarget_sites: list = field(default_factory=list)  # kept in memory, not all written out

    def has_editable_base(self, editor_class: str) -> bool:
        """Return whether the guide has an editable base inside the supported window."""
        source = EDIT_SOURCE_BASE[editor_class]
        return any(self.protospacer[p - 1] == source for p in EDIT_WINDOW if p - 1 < len(self.protospacer))

    def for_editor(self, editor_class: str, editor_name: str) -> Guide:
        """Pair this guide with an editor."""
        return Guide(
            guide_id=f"{self.guide_id}_{editor_name}",
            protospacer=self.protospacer,
            pam=self.pam,
            strand=self.strand,
            genome_build=self.genome_build,
            chrom=self.chrom,
            gene=self.gene,
            gene_strand=self.gene_strand,
            tss_1based=self.tss_1based,
            editor_class=editor_class,
            editor_name=editor_name,
            edit_from=EDIT_SOURCE_BASE[editor_class],
            edit_to=EDIT_TARGET_BASE[editor_class],
            proto_5p_genomic=self.proto_5p_genomic,
        )


@dataclass(frozen=True)
class OutcomeHaplotype:
    """Store one predicted output sequence and its predicted frequency."""

    output_seq: str
    pred_score: float


@dataclass
class EditSpectrum:
    """Store a predictor's possible output sequences and frequencies for one guide."""

    guide: Guide
    model_name: str
    haplotypes: list[OutcomeHaplotype]
    provenance: dict = field(default_factory=dict)

    def edited_positions(self, hap: OutcomeHaplotype) -> list[int]:
        """Protospacer positions this outcome changed."""
        ref = self.guide.protospacer
        out = hap.output_seq
        return [i + 1 for i in range(min(len(ref), len(out))) if ref[i] != out[i]]

    def per_position_rate(self) -> dict[int, float]:
        """Return each position's edit frequency summed across predicted output sequences."""
        rates: dict[int, float] = {p: 0.0 for p in self.guide.editable_positions()}
        for hap in self.haplotypes:
            for pos in self.edited_positions(hap):
                if pos in rates:
                    rates[pos] += hap.pred_score
        return rates

    def joint_edit_rate(self, proto_positions: set[int]) -> float:
        """Model-predicted frequency of haplotypes containing every requested edit.

        This sums complete predicted output sequences containing all requested edits. It is not
        reconstructed from separate position-level rates, and it is not a measurement.
        """
        if not proto_positions:
            raise ValueError("joint edit rate needs at least one protospacer position")
        invalid = proto_positions - set(self.guide.editable_positions())
        if invalid:
            raise ValueError(
                f"joint edit rate requested non-editable protospacer positions {sorted(invalid)} "
                f"for {self.guide.guide_id}"
            )
        return sum(
            hap.pred_score
            for hap in self.haplotypes
            if all(hap.output_seq[pos - 1] == self.guide.edit_to for pos in proto_positions)
        )


def spectrum_to_records(
    spectrum: EditSpectrum,
    on_target_offset: int | None = None,
    on_target_coords: set[int] | None = None,
    allow_overunity: bool = False,
) -> list[EditRecord]:
    """Turn a spectrum into one EditRecord per editable position.

    A position counts as on-target when it lands on a position this tool targets, matched either by
    promoter offset (on_target_offset, for a single guide) or by genomic coordinate
    (`on_target_coords`, when enumerating against the coordinate table). Other editable positions
    are bystanders.

    A per-position rate above 1 means the predictor's probabilities do not add up, which raises
    unless allow_overunity is set.

    Positions outside `EDIT_WINDOW` are not sent to the predictor. They use `predicted_rate=None` and
    `confidence=NOT_EVALUATED`. A numeric zero is reserved for a position the predictor evaluated and
    assigned zero. Keeping the unevaluated row makes this distinction visible.
    """
    guide = spectrum.guide
    rates = spectrum.per_position_rate()

    def is_on_target(proto_pos: int) -> bool:
        offset = guide.promoter_offset(proto_pos)
        coord = guide.genomic_coord(proto_pos)
        if on_target_coords is not None:
            return coord in on_target_coords
        return on_target_offset is not None and abs(offset) == abs(on_target_offset)

    intended_positions = {pos for pos in rates if is_on_target(pos)}
    joint_provenance: dict = {}
    if len(intended_positions) > 1:
        ordered_positions = sorted(
            intended_positions, key=lambda pos: abs(guide.promoter_offset(pos))
        )
        allele_keys = []
        for pos in ordered_positions:
            alleles = guide.alleles_for(pos)
            allele_keys.append({
                "promoter_offset": -abs(guide.promoter_offset(pos)),
                "promoter_sense_ref": alleles["promoter_sense_ref"],
                "promoter_sense_alt": alleles["promoter_sense_alt"],
            })
        evaluated_joint = all(pos in EDIT_WINDOW for pos in intended_positions)
        joint_rate = spectrum.joint_edit_rate(intended_positions) if evaluated_joint else None
        joint_provenance = {
            "joint_on_target_rate": round(joint_rate, 6) if joint_rate is not None else None,
            "joint_on_target_allele_keys": allele_keys,
            "joint_on_target_confidence": MODEL_PREDICTED if evaluated_joint else NOT_EVALUATED,
            "joint_on_target_method": JOINT_ON_TARGET_METHOD,
            "joint_on_target_model": spectrum.model_name,
            "joint_on_target_caveat": JOINT_ON_TARGET_CAVEAT,
        }

    records: list[EditRecord] = []
    for proto_pos, rate in sorted(rates.items()):
        evaluated = proto_pos in EDIT_WINDOW
        if evaluated and rate > 1.0 + 1e-6 and not allow_overunity:
            raise ValueError(
                f"per-position edit rate {rate:.4f} > 1 for {guide.guide_id} at protospacer position "
                f"{proto_pos}: the predictor's probabilities do not add up. Fix the predictor output, "
                f"or pass allow_overunity=True to keep the value and flag it."
            )
        offset = guide.promoter_offset(proto_pos)
        coord = guide.genomic_coord(proto_pos)
        position_is_on_target = is_on_target(proto_pos)
        # Omit an evaluated bystander with predicted rate zero. Keep an unevaluated position so the
        # output distinguishes "not evaluated" from a predicted zero.
        if evaluated and rate <= 0 and not position_is_on_target:
            continue
        confidence = MODEL_PREDICTED if evaluated else NOT_EVALUATED
        predicted_rate = round(rate, 6) if evaluated else None
        alleles = guide.alleles_for(proto_pos)
        role = "on_target" if position_is_on_target else "bystander"
        provenance = {
            "guide_id": guide.guide_id,
            "protospacer": guide.protospacer,
            "pam": guide.pam,
            "strand": guide.strand,
            "protospacer_position": proto_pos,
            "editor": guide.editor_name,
            "predictor_runs": spectrum.provenance.get("ensemble_runs"),
            "normalization_applied": spectrum.provenance.get("normalization_applied"),
            "raw_mass": spectrum.provenance.get("raw_mass"),
            **joint_provenance,
        }
        records.append(
            EditRecord(
                genome_build=guide.genome_build,
                chrom=guide.chrom,
                pos_1based=coord,
                genomic_plus_ref=alleles["genomic_plus_ref"],
                genomic_plus_alt=alleles["genomic_plus_alt"],
                promoter_sense_ref=alleles["promoter_sense_ref"],
                promoter_sense_alt=alleles["promoter_sense_alt"],
                gene=guide.gene,
                gene_strand=guide.gene_strand,
                promoter_offset=offset,
                editor_class=guide.editor_class,
                on_target_or_bystander=role,
                predicted_rate=predicted_rate,
                model_name=spectrum.model_name,
                confidence=confidence,
                variant_db_status="unset",
                provenance=provenance,
            )
        )
    return records
