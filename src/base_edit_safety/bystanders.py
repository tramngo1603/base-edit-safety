"""Classify additional edits made on the same predicted DNA molecule.

A rule needs the exact DNA change, not only its promoter position. For example, -110 A>G is an ABE
bystander, while -110 A>C is a different change. A rule can also depend on companion edits on the
same molecule. The -181 change is neutral only when -173 remains unedited.

Rules therefore use `(gene, offset, reference base, alternate base)` plus any required companion
edits. The guide and editor stay in the source record but do not decide which rule matches. ABE
changes A to G, while CBE changes C to T.

Offsets follow the gene-reading strand and count backward from the transcription start site. The
literature writes upstream positions as negative values; lookup code stores their absolute value.
"""
from __future__ import annotations

from dataclasses import dataclass, field

ADVERSE = "adverse"
PRODUCTIVE = "productive"
WEAKLY_PRODUCTIVE = "weakly_productive"
NEUTRAL = "neutral"
UNKNOWN = "unknown"

BOTH_GENES = ("HBG1", "HBG2")

# State of another edit on the same predicted sequence. `unresolved` is distinct from absence: a
# position outside the predictor's window has no rate, so a dependent rule cannot be evaluated.
PRESENT, ABSENT, UNRESOLVED = "present", "absent", "unresolved"


@dataclass(frozen=True)
class BystanderRule:
    """One bystander classification with its conditions and supporting source details."""

    genes: tuple
    offset: int                      # negative, sense strand, cap-site numbering
    ref: str                         # sense strand
    alt: str                         # sense strand
    functional_class: str
    note: str
    confidence: str
    # Companion edits on the same DNA molecule, each stored as `(offset, ref, alt)`.
    requires_present: tuple = ()
    requires_absent: tuple = ()
    # Supporting source fields are recorded but excluded from rule matching.
    observed_in_guide: tuple = ()
    observed_with_editor: tuple = ()
    observed_rate: str = ""
    source_pmid: tuple = ()
    evidence_class: str = ""

    @property
    def key(self) -> tuple:
        return (self.offset, self.ref, self.alt)


def _edit(offset: int, ref: str, alt: str) -> tuple:
    return (-abs(offset), ref, alt)


# The -113 ABE guide creates a GATA site spanning -113 to -110. -111 is a T, so -112 and -110 are the
# only ABE-editable positions inside it. Two papers report editing at both positions. Neither paper
# tests their functional effect together with -113A>G, so that effect has medium confidence.
_GATA_113 = _edit(113, "A", "G")
# The -175 ABE guide. Its spacer is published in Mayuranathan Supplementary Table 12 and was located
# in the GRCh38 reference on 2026-07-28: plus strand, 5' end at chr11:5,250,028. Every A in it maps to
# directly to an offset: A3 = -173, A5 = -175, and A11 = -181. These coordinates come from the
# sequence rather than an estimate based on the quoted "-175 (A5)" anchor.
_EBOX_175 = _edit(175, "A", "G")

# Spacers from Mayuranathan Supplementary Table 12, each located in the reference window.
PUBLISHED_SPACERS = {
    "sgRNA -175":     {"spacer": "AGATATTTGCATTGAGATAG", "pam": "NG", "strand": "+",
                       "five_prime_hbg1": 5_250_028,
                       "offsets": {3: -173, 5: -175, 11: -181, 15: -185, 17: -187, 19: -189}},
    "sgRNA -113":     {"spacer": "CTTGACCAATAGCCTTGACA", "pam": "NGG", "strand": "-",
                       "five_prime_hbg1": 5_249_977,
                       "offsets": {5: -116, 8: -113, 9: -112, 11: -110}},
    "sgRNA -187/189": {"spacer": "ATTGAGATAGTGTGGGGAAG", "pam": "NGG", "strand": "+",
                       "five_prime_hbg1": 5_250_038,
                       "offsets": {1: -181, 5: -185, 7: -187, 9: -189}},
    "sgRNA -198":     {"spacer": "GTGGGGAAGGGGCCCCCAA", "pam": "NGG", "strand": "+",
                       "five_prime_hbg1": 5_250_049, "offsets": {7: -198},
                       "note": "Mayuranathan Supplementary Table 12 reports "
                               "GTGGGGAAGGGGCCCCCAA (19 nt) with PAM class NGG. GRCh38 contains "
                               "the reference-derived 20 nt protospacer GTGGGGAAGGGGCCCCCAAG "
                               "followed by AGG; Ravi Supplementary file 1 independently reports "
                               "that 20-mer as gRNA 4. Accessible primary sources do not establish "
                               "whether the 19-mer was intentional or omitted a terminal G. "
                               "Concordant Ensembl and UCSC results close GRCh38 sequence and "
                               "coordinate validation for v0.1; this provenance discrepancy is "
                               "nonblocking."},
}

# Match published guides to the CRISPOR enumeration. Nineteen of the 22 distinct guide-matched
# sequences were screened by Ravi, so the papers evaluate the same sequences rather than only the
# same promoter positions. Established 2026-07-28 from Ravi Supplementary file 1 and Mayuranathan
# Supplementary Table 12.
GUIDE_IDENTITIES = {
    "CTTGACCAATAGCCTTGACA": {"ours": "19rev",   "ravi": "gRNA 2",  "mayuranathan": "sgRNA -113"},
    "AGGCAAGGCTGGCCAACCCA": {"ours": "60forw",  "ravi": "gRNA 10", "mayuranathan": None},
    "GGCAAGGCTGGCCAACCCAT": {"ours": "61forw",  "ravi": "gRNA 11", "mayuranathan": None},
    "ATTGAGATAGTGTGGGGAAG": {"ours": "122forw", "ravi": "gRNA 25", "mayuranathan": "sgRNA -187/189"},
    "GCATTGAGATAGTGTGGGGA": {"ours": "120forw", "ravi": "gRNA 24", "mayuranathan": None},
    "ATATTTGCATTGAGATAGTG": {"ours": "114forw", "ravi": "gRNA 3",  "mayuranathan": None},
    "TATTTGCATTGAGATAGTGT": {"ours": "115forw", "ravi": "gRNA 18", "mayuranathan": None},
    "CTTGTCAAGGCTATTGGTCA": {"ours": "40forw",  "ravi": "gRNA 7",  "mayuranathan": "sgRNA HBG"},
}


def catastrophic_exposure(guide, edit_window) -> list[tuple]:
    """Return high-impact positions that this guide can edit inside its supported window.

    In the cited study, disruption at -187/-189 prevented the -175 edit from increasing HbF and
    removed binding by GATA1, TAL1, LMO2, and LDB1. The result contains `(guide position, promoter
    offset)` pairs and is empty when neither position is reachable.
    """
    reachable = {c.offset for c in CATASTROPHIC}
    found = []
    for pos in guide.editable_positions():
        if pos not in edit_window:
            continue
        offset = -abs(guide.promoter_offset(pos))
        if offset in reachable:
            found.append((pos, offset))
    return found

RULES: tuple[BystanderRule, ...] = (
    BystanderRule(
        genes=BOTH_GENES, offset=-116, ref="A", alt="G",
        functional_class=WEAKLY_PRODUCTIVE,
        requires_present=(_GATA_113,),
        note="Downgraded from productive. Han 2023 measures BCL11A binding affinity reduced only "
             "3.7x by -116A>G, against 31x for -115C>T and 55x for -114C>T, and the ABE8e strategy "
             "installing it gave the lowest gamma-globin across the tested conditions despite >90% editing.",
        confidence="medium-high",
        observed_in_guide=("sgHBG-1", "gRNA-2"), observed_with_editor=("ABE8e",),
        observed_rate="92.2% (Han 2023, sgHBG-1); 66% (Ravi 2022, gRNA-2, CD34+ HSPCs)",
        source_pmid=("37989316", "35147495"), evidence_class="measured_rate_and_binding_affinity",
    ),
    BystanderRule(
        genes=BOTH_GENES, offset=-112, ref="A", alt="G",
        functional_class=ADVERSE,
        requires_present=(_GATA_113,),
        note="Inside the GATA site the -113 edit creates. Mayuranathan 2023 says the A9 bystander "
             "'should inhibit HbF induction by destroying the GATA motif created by the on-target "
             "edit'. That is a stated expectation, not a measurement: no paper has combined "
             "-112A>G with -113A>G. The adverse call rests on motif structure alone.",
        confidence="high for position and rate; medium for the adverse class",
        observed_in_guide=("sgHBG-1", "gRNA-2"), observed_with_editor=("ABE8e",),
        observed_rate="40.4% (Han 2023, sgHBG-1); 37% (Ravi 2022, gRNA-2, CD34+ HSPCs)",
        source_pmid=("37989316", "35147495", "37400614"), evidence_class="motif_structure_only",
    ),
    BystanderRule(
        genes=BOTH_GENES, offset=-110, ref="A", alt="G",
        functional_class=ADVERSE,
        requires_present=(_GATA_113,),
        note="Fourth base of the GATA site created by the -113 edit. Doerfler 2021 introduced "
             "-110A>G alone as a control "
             "and measured no HbF change by IE-HPLC and no accessibility change by ATAC-seq. That "
             "is the exact allele, tested alone in cells, with no effect. It does not test the "
             "bystander case, because it was never combined with -113A>G. The adverse call rests "
             "entirely on motif structure.",
        confidence="high for position and rate; medium for the adverse class because the combined "
                   "-110A>G and -113A>G change was not tested",
        observed_in_guide=("sgHBG-1", "gRNA-2"), observed_with_editor=("ABE8e",),
        observed_rate="80.0% (Han 2023, sgHBG-1); 31% (Ravi 2022, gRNA-2, CD34+ HSPCs)",
        source_pmid=("37989316", "35147495", "34341563"),
        evidence_class="motif_structure_only_with_contrary_isolated_null",
    ),
    # -114 and -115 are absent from the section-2 change table. Retain them only for supported CBE/tBE
    # BCL11A-motif sequence contexts rather than applying them as global offset rules.
    BystanderRule(
        genes=BOTH_GENES, offset=-114, ref="C", alt="T",
        functional_class=PRODUCTIVE,
        note="The intended distal-BCL11A-motif edit. Productive in the reported CBE and "
             "transformer-base-editor guide contexts, and natural -114 alleles support biological "
             "relevance. The evidence does not justify a global editor-independent offset rule, so "
             "it is keyed to this allele. Note the literature also carries -114 C>A (Doerfler, "
             "genomic) and C>G (Doerfler, EMSA oligo); neither is this edit.",
        confidence="high for this DNA change in the cited CBE and tBE contexts; not a general "
                   "position-only rule",
        observed_with_editor=("BE4max", "transformer base editor"),
        source_pmid=("37989316", "35147495"), evidence_class="supported_haplotype_context",
    ),
    BystanderRule(
        genes=BOTH_GENES, offset=-115, ref="C", alt="T",
        functional_class=PRODUCTIVE,
        requires_present=(_edit(114, "C", "T"),),
        note="Participates in productive BCL11A-motif-disrupting CBE/tBE edit mixtures. The cited "
             "work does not isolate -115 as a universally productive single edit, so it is "
             "conditioned on -114 being edited too rather than standing alone.",
        confidence="medium",
        observed_with_editor=("BE4max", "transformer base editor"),
        source_pmid=("37989316",), evidence_class="supported_haplotype_context",
    ),
    BystanderRule(
        genes=BOTH_GENES, offset=-173, ref="A", alt="G",
        functional_class=ADVERSE,
        requires_present=(_EBOX_175,),
        note="Mayuranathan 2023, designed clonal comparison: G3G5G11 significantly reduced %HbF "
             "against both A3G5A11 and A3G5G11, and the authors state the A3 edit likely "
             "interferes with HbF induction by -175A>G. -173 falls inside the created E-box "
             "(sense -175 to -170), which corroborates the offset.",
        confidence="high for the result; coordinate derived directly from the published guide "
                   "sequence",
        observed_in_guide=("sgRNA -175 = AGATATTTGCATTGAGATAG (Mayuranathan Supp Table 12); "
                           "plus strand, 5' at chr11:5,250,028, so A3 = -173",),
        observed_with_editor=("ABE8e",),
        observed_rate="not reported per position; read out as clonal haplotypes",
        source_pmid=("37400614",), evidence_class="designed_clonal_comparison",
    ),
    BystanderRule(
        genes=BOTH_GENES, offset=-181, ref="A", alt="G",
        functional_class=NEUTRAL,
        requires_present=(_EBOX_175,),
        requires_absent=(_edit(173, "A", "G"),),
        note="A3G5A11 and A3G5G11 are both productive edit combinations in Mayuranathan 2023, so "
             "the A11 edit does not reduce induction when -173 is unedited. When -173 and -181 are "
             "both edited, the adverse effect is assigned to -173. A second guide also reaches "
             "-181, but that experiment does not isolate its effect.",
        confidence="high for the A11 result; coordinate derived directly from the published guide "
                   "sequence",
        observed_in_guide=("sgRNA -175 = AGATATTTGCATTGAGATAG, plus strand 5' at chr11:5,250,028, "
                           "so A11 = -181",
                           "sgRNA -187/189 = ATTGAGATAGTGTGGGGAAG, plus strand 5' at "
                           "chr11:5,250,038, so A1 = -181 by a second, different guide"),
        observed_with_editor=("ABE8e", "ABE7.10"),
        observed_rate="not reported per position; read out as clonal haplotypes",
        source_pmid=("37400614", "34341563"), evidence_class="designed_clonal_comparison",
    ),
)

_BY_KEY: dict[tuple, BystanderRule] = {r.key: r for r in RULES}

# Positions where disruption caused a severe loss of function in the cited study. They are not
# bystanders of a guide currently in scope: Mayuranathan used a separate guide, and these positions
# fall outside the supported ABE window for the -175 guide. Keep them in the rule table so a future
# candidate that includes either position receives the recorded high-impact label.
@dataclass(frozen=True)
class CatastrophicPosition:
    genes: tuple
    offset: int
    element: str
    note: str
    source_pmid: tuple


CATASTROPHIC = (
    CatastrophicPosition(
        BOTH_GENES, -187, "upstream GATA1 motif (WGATAR), sense -190 to -185",
        "ABE7.10 disruption at this position produced clones in which -175G did not increase HbF. "
        "CUT&RUN detected no GATA1, TAL1, LMO2, or LDB1 occupancy. Reached by sgRNA -187/189 = "
        "ATTGAGATAGTGTGGGGAAG (A7), and also inside sgRNA -175 at A17.", ("37400614",)),
    CatastrophicPosition(
        BOTH_GENES, -189, "upstream GATA1 motif (WGATAR), sense -190 to -185",
        "The same experiment and result as -187. This position is A9 in sgRNA -187/189 and A19 in "
        "sgRNA -175. The latter guide spans promoter positions -175 through -189.", ("37400614",)),
)

# Motif geometry at -175, measured from Wienert 2015 Fig 1a. The same substitution completes one
# motif and destroys another at the same base, which is why a single-consequence label is incomplete.
MOTIF_GEOMETRY_175 = {
    "created_tal1_ebox": (-175, -170),
    "destroyed_gata1_core": (-175, -172),
    "redundant_gata1_site": (-190, -185),   # genomic GATA1 ChIP was unchanged at this second site
    "inside_both_ebox_and_destroyed_gata": (-174, -172),
}


def allele_context(records: list[dict]) -> dict:
    """Which edits sit on one guide's allele, as (offset, ref, alt) -> present | absent | unresolved.

    A rate above zero means present. A predicted rate of zero means absent. `None` means the position
    was outside the request sent to the predictor, so presence and dependent rules remain unresolved.
    """
    context: dict = {}
    for r in records:
        key = (-abs(r["promoter_offset"]), r.get("promoter_sense_ref"), r.get("promoter_sense_alt"))
        rate = r.get("predicted_rate")
        state = UNRESOLVED if rate is None else (PRESENT if rate > 0 else ABSENT)
        if context.get(key) != PRESENT:      # keep `present` if duplicate records disagree
            context[key] = state
    return context


def _presence(edit: tuple, edits_on_allele: dict) -> str:
    """Is this edit on the allele? Unknown is distinct from absent."""
    return edits_on_allele.get(edit, ABSENT)


def bystander_function(gene: str, promoter_offset: int, ref: str, alt: str,
                       edits_on_allele: dict | None = None) -> tuple[str, str]:
    """Return the class and note for one bystander, given its companion edits.

    `edits_on_allele` maps (negative offset, ref, alt) to PRESENT, ABSENT or UNRESOLVED. A rule whose
    conditions cannot be evaluated returns unknown. This prevents an unsupported class from
    changing the design's rank.
    """
    edits_on_allele = edits_on_allele or {}
    rule = _BY_KEY.get((-abs(promoter_offset), ref, alt))
    if rule is None:
        return UNKNOWN, "no rule for this edit; apply the configured unknown-edit penalty"
    if gene not in rule.genes:
        return UNKNOWN, f"rule at {rule.offset} {ref}>{alt} does not cover {gene}"

    for required in rule.requires_present:
        state = _presence(required, edits_on_allele)
        if state is UNRESOLVED or state == UNRESOLVED:
            return UNKNOWN, (f"depends on {required[0]} {required[1]}>{required[2]}, whose presence "
                             f"is unresolved; the rule cannot be evaluated")
        if state != PRESENT:
            return UNKNOWN, (f"only classified when {required[0]} {required[1]}>{required[2]} is "
                             f"also edited, and it is not")
    for forbidden in rule.requires_absent:
        state = _presence(forbidden, edits_on_allele)
        if state == UNRESOLVED:
            return UNKNOWN, (f"depends on {forbidden[0]} {forbidden[1]}>{forbidden[2]} being "
                             f"unedited, and its presence is unresolved")
        if state == PRESENT:
            return UNKNOWN, (f"{rule.functional_class} only while "
                             f"{forbidden[0]} {forbidden[1]}>{forbidden[2]} is unedited, and it is "
                             f"edited here")
    return rule.functional_class, rule.note


def rule_for(promoter_offset: int, ref: str, alt: str) -> BystanderRule | None:
    """Return the matching rule with its conditions and supporting source fields."""
    return _BY_KEY.get((-abs(promoter_offset), ref, alt))
