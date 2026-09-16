"""Compare a candidate with recorded HBG variant and editing evidence.

The list contains natural promoter variants, published HSPC editing work and clinical program
references. Match type (exact_edit, same_motif, nearby or no_match) stays separate from evidence
strength. Unknown program targets stay unknown. Source-review status limits each claim.
"""
from __future__ import annotations

from dataclasses import dataclass

from . import anchors

# How a candidate lines up with the list, structurally.
EXACT_EDIT = "exact_edit"
SAME_MOTIF = "same_motif"
NEARBY = "nearby"
NO_MATCH = "no_match"

# Evidence strength, weakest to strongest. Kept separate from the evidence tier in evidence.py, so
# something published elsewhere never turns into `empirically_measured`, which means measured here.
STRENGTH_ORDER = ("none", "analog", "motif", "published_hspc", "clinical")

# Mechanism support is separate from evidence that a variant exists.
# Use NOT_REVIEWED until the mechanism's source has been checked.
#   in_vitro_binding_only  gel-shift or competition work, no cellular confirmation that discriminates
#   binding_plus_cellular  in vitro binding plus a cellular assay that rules alternatives out
#   established            the mechanism is not seriously in question
#   not_reviewed           the paper behind the label has not completed source review
IN_VITRO_BINDING_ONLY = "in_vitro_binding_only"
BINDING_PLUS_CELLULAR = "binding_plus_cellular"
ESTABLISHED = "established"
NOT_REVIEWED = "not_reviewed"

# Whether anyone has read the cited paper.
FULL_TEXT_READ = "full_text_read"
NEVER_REVIEWED = "never_reviewed"


def _cite(pmid: str, role: str, audit: str = NOT_REVIEWED, note: str | None = None) -> dict:
    """Record a citation, the claim it supports and its source-review status."""
    entry = {"pmid": pmid, "role": role, "audit_status": audit}
    if note:
        entry["note"] = note
    return entry


@dataclass(frozen=True)
class RegistryEntry:
    entry_id: str
    label: str
    genes: tuple
    modality: str                 # base_editing | nuclease
    editor_class: str | None      # ABE | CBE | None (nucleases, or an unconfirmed editor class)
    offsets: tuple                # promoter offsets as positive numbers, or () if not position-specific
    motif: str                    # BCL11A | GATA1 | KLF1 | TAL1 | ...
    strength: str                 # one of STRENGTH_ORDER
    citations: tuple              # dicts from _cite()
    note: str
    mechanism_evidence: str = NOT_REVIEWED
    program_status: str | None = None   # clinical entries only: active | discontinued_<date>


# --- What system a claim was made in -------------------------------------------------------------
#
# Architecture belongs to the CLAIM, not the PMID. Wienert 2015 alone spans two systems, and half its
# claims come from each. Mayuranathan and Ravi each span an engineered line and native primary cells.
#
# Check `intact_paralogs`: a mechanism shown in a system where HBG2 was deleted
# cannot be projected onto HBG2. None of the five mechanism papers used an intact native locus.

@dataclass(frozen=True)
class LocusArchitecture:
    system: str
    description: str
    intact_paralogs: tuple          # which of HBG1/HBG2 were present and unmodified
    native: bool = False


WIENERT_MEL_BAC = LocusArchitecture(
    "MEL (mouse)", "human beta-globin BAC; Agamma coding replaced by ECFP; PGK-NeoR retained",
    intact_paralogs=())
WIENERT_K562 = LocusArchitecture(
    "K562", "Ggamma gene body deleted; tdTomato knock-in; NeoR retained",
    intact_paralogs=("HBG1",))
HUDEP2_DELTA_GGAMMA = LocusArchitecture(
    "HUDEP-2", "Ggamma deleted", intact_paralogs=("HBG1",))
DOERFLER_FUSED = LocusArchitecture(
    "HUDEP-2", "one chromosome's beta-cluster deleted; the gamma genes fused into a single "
               "HBG2-HBG1 gene", intact_paralogs=())
HUDEP2_INTACT = LocusArchitecture(
    "HUDEP-2", "unmodified locus in an immortalised erythroid line",
    intact_paralogs=("HBG1", "HBG2"))
PRIMARY_CD34 = LocusArchitecture(
    "primary CD34+ HSPCs", "native locus", intact_paralogs=("HBG1", "HBG2"), native=True)
PRIMARY_CD34_TRANSPLANT = LocusArchitecture(
    "primary CD34+ HSPCs with xenotransplantation", "native locus; engraftment demonstrated",
    intact_paralogs=("HBG1", "HBG2"), native=True)

# How far anyone got into a source. A claim may not assert more than its access allows.
ACCESS_ORDER = ("not_accessed", "abstract_only", "database_record", "full_text")
FULL_TEXT, ABSTRACT_ONLY, NOT_ACCESSED = "full_text", "abstract_only", "not_accessed"

# Record transplantation/engraftment evidence separately from an HSPC evidence label.
# An HSPC experiment alone does not show durability after transplantation.
NO_ENGRAFTMENT_DATA = "no_engraftment_data"
XENOTRANSPLANT = "xenotransplant_engraftment_shown"


def _claim(pmid: str, role: str, access: str = NOT_ACCESSED, audit: str = NOT_REVIEWED,
           architecture: LocusArchitecture | None = None, engraftment: bool = False,
           note: str | None = None) -> dict:
    """Record a claim, the experimental system and the source material reviewed."""
    entry = {"pmid": pmid, "role": role, "access_level": access, "audit_status": audit,
             "demonstrates_engraftment": engraftment}
    if architecture is not None:
        entry["locus_architecture"] = {
            "system": architecture.system,
            "description": architecture.description,
            "intact_paralogs": list(architecture.intact_paralogs),
            "native": architecture.native,
        }
    if note:
        entry["note"] = note
    return entry


def _registry_ref(identifier: str, role: str, audit: str = NOT_REVIEWED,
                  note: str | None = None) -> dict:
    """A registry ID or unresolved non-PMID source label, kept distinct from paper PMIDs."""
    entry = {"registry_id": identifier, "role": role, "audit_status": audit,
             "access_level": NOT_ACCESSED}
    if note:
        entry["note"] = note
    return entry


# --- Claims per allele ---------------------------------------------------------------------------
#
# Keyed by (offset, sense ref, sense alt), not by offset. -114 carries C>T here, C>A as Doerfler's
# cellular allele and C>G as its EMSA oligo; filing all three under "-114" attaches a paper about one
# allele to a claim about another.

# Doerfler splits two and two rather than sitting on four positions equally.
_DOERFLER_TESTED = lambda branch: _claim(
    "34341563", f"qualified_variant_specific_{branch}", FULL_TEXT, FULL_TEXT_READ, DOERFLER_FUSED,
    note="Retained only for the dependency branch actually tested, with the fusion-promoter caveat. "
         "Not proof of KLF1 binding, KLF1-ZBTB7A competition, or direct GATA1-versus-BCL11A "
         "displacement.")
_DOERFLER_FRAMEWORK = _claim(
    "34341563", "competition_framework_only", FULL_TEXT, FULL_TEXT_READ, DOERFLER_FUSED,
    note="Framework only at this position, recorded once at the disrupted-BCL11A-motif class.")

_CITATIONS = {
    (198, "T", "C"): (
        _claim("2430647", "human_haplotype_report", FULL_TEXT, FULL_TEXT_READ,
               note="n=1 haplotype also carrying -158 and +210. The authors explicitly decline to "
                    "attribute the phenotype to -198."),
        _claim("28659276", "mechanism_klf1", NOT_ACCESSED, NOT_REVIEWED, HUDEP2_DELTA_GGAMMA,
               note="KLF1 was nominated on erythroid abundance, not binding exclusivity; twelve "
                    "other KLFs bound the mutant probe."),
        _claim("37400614", "hspc_base_editing", NOT_ACCESSED, NOT_REVIEWED,
               PRIMARY_CD34_TRANSPLANT, engraftment=True),
        _DOERFLER_TESTED("upstream_gata1_dependent"),
    ),
    (175, "T", "C"): (
        _claim("2449926", "human_variant_report_hbg2", FULL_TEXT, FULL_TEXT_READ,
               note="n=1, in cis with -158 and 3'-flanking changes. Proposes an OCTAMER, not TAL1. "
                    "Must never be cited for a TAL1 claim; TAL1 belongs to 25971621."),
        _claim("25971621", "mechanism_tal1_ebox", FULL_TEXT, FULL_TEXT_READ, WIENERT_MEL_BAC,
               note="Measures NO gamma-globin and NO HbF: the gamma coding sequence is replaced by "
                    "a reporter, so every value called 'gamma-globin expression' is fluorescence. "
                    "The MEL ChIP amplicon almost certainly does not distinguish the paralogs while "
                    "MEL retains both promoters, so it averages a mutant and a WT promoter. Error "
                    "direction is conservative, so the conclusion holds but the attribution does "
                    "not: occupancy rises at the gamma promoter REGION, not demonstrably at the "
                    "Agamma -175T>C promoter."),
        _claim("25971621", "mechanism_tal1_ebox_second_system", FULL_TEXT, FULL_TEXT_READ,
               WIENERT_K562,
               note="Second architecture in the same paper. Same reporter caveat."),
        _claim("37400614", "hspc_base_editing", NOT_ACCESSED, NOT_REVIEWED,
               PRIMARY_CD34_TRANSPLANT, engraftment=True),
        _claim("7679879", "human_variant_report_hbg1", NOT_ACCESSED, NEVER_REVIEWED,
               note="Coleman 1993, HBG1 -175 in cis with HBG2 -158. 2 of 2 known Agamma carriers, "
                    "but phase is asserted by a compilation and all six underlying primaries are "
                    "unread."),
    ),
    (124, "T", "C"): (
        _claim("35147495", "synthetic_screen_hit", FULL_TEXT, FULL_TEXT_READ, HUDEP2_INTACT,
               note="Discovery screen was HUDEP-2, not HSPCs."),
        _claim("35147495", "allele_installed_in_hspcs", FULL_TEXT, FULL_TEXT_READ, PRIMARY_CD34,
               note="Exact allele installed and genotyped in primary CD34+ HSPCs at day 6, before "
                    "differentiation. NO xenotransplantation, so this does not show durable "
                    "editing in repopulating HSCs."),
    ),
    (123, "T", "C"): (
        _claim("35147495", "synthetic_screen_hit", FULL_TEXT, FULL_TEXT_READ, HUDEP2_INTACT),
        _claim("35147495", "allele_installed_in_hspcs", FULL_TEXT, FULL_TEXT_READ, PRIMARY_CD34,
               note="See -124: same paper, same absence of transplant data. KLF1 binding was shown "
                    "for the -123/-124 DOUBLE substitution only, and the ChIP negative control "
                    "moved similarly while the authors name alternative factors."),
    ),
    (122, "G", "A"): (
        _claim("35147495", "cbe_on_target_conversion", FULL_TEXT, FULL_TEXT_READ, HUDEP2_INTACT,
               note="Stated in both Results and Discussion as the on-target site for gRNA-10/-11. "
                    "CBE work in this paper was HUDEP-2, not HSPCs."),
    ),
    (117, "G", "A"): (
        _claim("2417646", "phenotype_concordance", FULL_TEXT, FULL_TEXT_READ,
               note="NOT a variant description. Its own introduction credits the -117 description "
                    "to Gelinas 1985 and Collins 1985, neither cited anywhere in this repo. There "
                    "is NO existence citation for -117 here."),
        _claim("3181130", "mechanism_distal_ccaat", NOT_ACCESSED, NOT_REVIEWED),
        _claim("37989316", "motif_base_editing", NOT_ACCESSED, NOT_REVIEWED, PRIMARY_CD34),
        _DOERFLER_FRAMEWORK,
    ),
    (114, "C", "T"): (
        _claim("1698280", "human_variant_report_hbg2", FULL_TEXT, FULL_TEXT_READ,
               note="HBG2 only, a mother-son pair. Agamma was NEVER SEQUENCED in any subject and no "
                    "Agamma sequence entered any assay, so Oner 1991's 'one can assume' at HBG1 "
                    "-114 has no HBG1 claim to assume from."),
        _claim("1704803", "correspondence_n2_siblings", FULL_TEXT, FULL_TEXT_READ,
               note="Two brothers, published as correspondence. Its paralog claim is an explicit "
                    "assumption in the paper's own words, not a measurement."),
        _claim("35147495", "regional_base_editing", FULL_TEXT, FULL_TEXT_READ, HUDEP2_INTACT,
               note="Ravi's -114 editing was HUDEP-2 only. -114C>T needs a CBE and no CBE was "
                    "applied to primary HSPCs in that paper, so 'HSPC base editing' is wrong here."),
        _claim("37989316", "motif_base_editing", NOT_ACCESSED, NOT_REVIEWED, PRIMARY_CD34),
        _claim("34341563", "competition_framework_only", FULL_TEXT, FULL_TEXT_READ, DOERFLER_FUSED,
               note="Framework only. Doerfler's cellular allele here is C>A and its EMSA oligo is "
                    "C>G; NEITHER is this repo's C>T. The CP1 mechanism rests on EMSA competition "
                    "with no antibody and no supershift, and DNaseI footprinting failed to resolve "
                    "the distal motif in either probe."),
    ),
    (113, "A", "G"): (
        _claim("23621512", "review_reporting_single_case", FULL_TEXT, FULL_TEXT_READ,
               note="A review, not a primary report. Its single -113 subject also carries the "
                    "beta-gene variant IVS 2-74 G>T, so the phenotype is not attributable to -113 "
                    "alone."),
        _claim("30617196", "mechanism_de_novo_gata1", NOT_ACCESSED, NOT_REVIEWED,
               HUDEP2_DELTA_GGAMMA),
        _claim("36006707", "exact_allele_base_editing", NOT_ACCESSED, NOT_REVIEWED, PRIMARY_CD34),
        _claim("37400614", "hspc_base_editing", NOT_ACCESSED, NOT_REVIEWED,
               PRIMARY_CD34_TRANSPLANT, engraftment=True),
        _DOERFLER_TESTED("proximal_nfy_independent"),
    ),
}


def durability_for(offset: int, ref: str, alt: str) -> str:
    """Return the recorded engraftment evidence separately from other evidence labels."""
    claims = _CITATIONS.get((abs(offset), ref, alt), ())
    return XENOTRANSPLANT if any(c.get("demonstrates_engraftment") for c in claims) \
        else NO_ENGRAFTMENT_DATA


def mechanism_projectable_onto(gene: str, offset: int, ref: str, alt: str) -> bool:
    """Can a mechanism claim here be projected onto this paralog?

    Only if some mechanism claim was made in a system where that paralog's locus was intact. A
    mechanism shown in a HUDEP-2 line with Ggamma deleted says nothing about HBG2.
    """
    for c in _CITATIONS.get((abs(offset), ref, alt), ()):
        if not c["role"].startswith("mechanism"):
            continue
        arch = c.get("locus_architecture")
        if arch and gene in arch["intact_paralogs"]:
            return True
    return False


# How well the mechanism named for each position is supported. Three have been read; all three turned
# out to rest on in vitro binding alone.
#   -198 KLF1     KLF1 was nominated on erythroid abundance, not on binding exclusivity, and twelve
#                 other KLFs bound the mutant probe.
#   -123/-124     the EMSA is good, but the ChIP negative control moved similarly and the authors
#                 name alternative factors themselves.
#   -114 CP1      EMSA competition only, with no antibody and no supershift, and DNaseI footprinting
#                 failed to resolve the distal motif in either probe.
# The rest are unread, so they stay NOT_REVIEWED rather than being assumed sound.
_MECHANISM_EVIDENCE = {
    198: IN_VITRO_BINDING_ONLY,
    124: IN_VITRO_BINDING_ONLY,
    123: IN_VITRO_BINDING_ONLY,
    114: IN_VITRO_BINDING_ONLY,
    175: NOT_REVIEWED,
    117: NOT_REVIEWED,
    113: NOT_REVIEWED,
}

# Evidence strength per HPFH position (see docs/coordinates.md).
_STRENGTH = {
    198: "published_hspc",   # HSPC ABE (Mayuranathan 2023)
    175: "published_hspc",   # strongest exact-allele HSPC ABE (Mayuranathan 2023)
    124: "published_hspc",   # the synthetic -123/-124 pair in HSPC-derived erythroid cells (Ravi 2022)
    123: "published_hspc",
    117: "motif",            # covered by a tiling screen, not the exact allele (Han 2023)
    114: "motif",            # covered by a tiling screen (Ravi 2022 / Han 2023)
    113: "published_hspc",   # exact -113 A>G ABE8e in patient CD34+ HSCs (Li 2022)
}


def _motif_for(mechanism: str) -> str:
    m = mechanism.upper()
    for motif in ("BCL11A", "GATA1", "KLF1", "TAL1", "CCAAT"):
        if motif in m:
            return "BCL11A" if motif == "CCAAT" else motif   # the distal CCAAT box is the BCL11A site
    return "unknown"


def _entry_id(pos) -> str:
    """Use the hpfh_ prefix only for an allele recorded as natural in the anchor table.

    The constructed -123/-124 screen changes receive a different prefix.
    """
    is_synthetic = {pos.hbg1_status == anchors.SYNTHETIC, pos.hbg2_status == anchors.SYNTHETIC}
    if len(is_synthetic) > 1:
        # Check both genes before assigning a shared natural/synthetic prefix.
        # Their database statuses can differ; a mixed natural/synthetic pair needs separate entries.
        raise ValueError(
            f"offset {pos.promoter_offset}: HBG1 is {pos.hbg1_status} and HBG2 is "
            f"{pos.hbg2_status}, so one entry_id cannot describe both paralogs. Make the id "
            f"gene-scoped before adding this row."
        )
    prefix = "synthetic" if is_synthetic.pop() else "hpfh"
    return f"{prefix}_-{pos.promoter_offset}_{pos.editor_class}"


def _hpfh_entries() -> list[RegistryEntry]:
    return [
        RegistryEntry(
            entry_id=_entry_id(pos),
            label=f"HBG -{pos.promoter_offset} {pos.editor_class} ({pos.mechanism})",
            genes=("HBG1", "HBG2"), modality="base_editing", editor_class=pos.editor_class,
            offsets=(pos.promoter_offset,), motif=_motif_for(pos.mechanism),
            strength=_STRENGTH.get(pos.promoter_offset, "motif"),
            citations=_CITATIONS.get(
                (pos.promoter_offset, pos.promoter_sense_ref, pos.promoter_sense_alt), ()),
            note=("A screen hit, not a natural allele. The source paper says HPFH-LIKE throughout, "
                  "including in its title."
                  if pos.hbg1_status == anchors.SYNTHETIC
                  else "A natural HPFH variant, published HSPC base-editing work, or both."),
            mechanism_evidence=_MECHANISM_EVIDENCE.get(pos.promoter_offset, NOT_REVIEWED),
        )
        for pos in anchors.IN_SCOPE
    ]


# Clinical HBG1/HBG2 promoter programs, recorded no more specifically than the available source
# supports because exact positions and motifs are not always public.
_CLINICAL = [
    RegistryEntry(
        "beam_101", "risto-cel, formerly BEAM-101 (Beam Therapeutics, BEACON)", ("HBG1", "HBG2"),
        "base_editing", None, (), "BCL11A", "clinical",
        (_registry_ref("NCT05456880", "trial_registration", note="BEACON"),
         _registry_ref("Beam Therapeutics BEAM-101 program materials (public)", "sponsor_material",
                       note="Not resolvable as written; needs a dated, specific source.")),
        "Base-edited CD34+ HSPCs in the HBG1/HBG2 promoter, aimed at stopping BCL11A from binding, "
        "the way HPFH variants do. Same approach and same target as this tool. Renamed to "
        "ristoglogene autogetemcel (risto-cel); the entry_id keeps the old name so nothing "
        "downstream breaks. Public reporting says the editor is an ABE. The supporting sponsor "
        "source has not completed review, so editor_class stays unasserted.",
        program_status="active",
    ),
    RegistryEntry(
        "reni_cel_edit301", "reni-cel / EDIT-301 (Editas, RUBY / EdiThal)", ("HBG1", "HBG2"),
        "nuclease", None, (), "BCL11A", "clinical",
        (_registry_ref("NCT04853576", "trial_registration", note="RUBY"),
         _registry_ref("NCT05444894", "trial_registration", note="EdiThal"),
         _registry_ref("Editas reni-cel program materials (8-K summary)", "sponsor_material",
                       note="Not resolvable as written; needs the SEC accession and a date.")),
        "AsCas12a cutting at the HBG1/HBG2 promoter BCL11A site. Same place and same idea, but a "
        "nuclease rather than a base editor, so it makes indels instead of point edits. "
        "DISCONTINUED: Editas ended development in December 2024 after failing to find a commercial "
        "partner, and both trials were stopped early for that reason, not for safety or efficacy. "
        "The human evidence still stands, which is why strength remains 'clinical'.",
        program_status="discontinued_2024_12",
    ),
    # These entries complete the known promoter-program list. Their cited sources have not
    # completed review, so every identifier they carry remains marked unread.
    RegistryEntry(
        "cs_101", "CS-101 / CS-206 (CorrectSequence Therapeutics)", ("HBG1", "HBG2"),
        "base_editing", None, (), "BCL11A", "clinical",
        (_registry_ref("NCT06024876", "trial_registration", NEVER_REVIEWED),
         _registry_ref("NCT06291961", "trial_registration", NEVER_REVIEWED),
         _registry_ref("NCT06565026", "trial_registration", NEVER_REVIEWED)),
        "A clinical base editor at the HBG1/HBG2 promoter BCL11A site, which makes it the closest "
        "analog anywhere to what this tool models. It was missing from this list entirely. "
        "Transformer base editor; the subclass is not public, so editor_class is None.",
        program_status="active",
    ),
    RegistryEntry(
        "otq923", "OTQ923 (Novartis)", ("HBG1", "HBG2"),
        "nuclease", None, (), "BCL11A", "clinical",
        (_registry_ref("NCT04443907", "trial_registration", NEVER_REVIEWED),
         _cite("37646679", "clinical_trial_report", NEVER_REVIEWED,
               "Sharma et al., NEJM 2023, Cas9 disruption of the HBG1 and HBG2 promoters. Cited by "
               "PMID because this repo ships no DOIs.")),
        "Cas9-edited CD34+ HSPCs with a targeted disruption of both gamma-globin promoters. Same "
        "locus as this tool, nuclease rather than base editor. Read out in 2023.",
        program_status="discontinued",
    ),
    RegistryEntry(
        "rm_001", "RM-001 (Reforgene Therapeutics / 瑞风生物, Guangzhou)",
        ("HBG1", "HBG2"), "nuclease", None, (), "unknown", "clinical",
        (_registry_ref(
            "sponsor materials (unopened)", "sponsor_material", NEVER_REVIEWED,
            "RM-001 trial registry identifiers were not verified. This placeholder needs a dated, "
            "specific source."
        ),),
        "Cas9 editing reported directly at the HBG1/HBG2 promoter region, distinct from "
        "BCL11A-enhancer programs; the exact promoter motif is unverified. The Phase I study is "
        "reported complete, and company "
        "announcements from late 2025 report a completed confirmatory trial now in follow-up. "
        "Those sponsor materials were not opened here, and no trial registry identifier is "
        "asserted.",
        program_status="active",
    ),
    RegistryEntry(
        "yolt_204", "YOLT-204 (YolTech Therapeutics)",
        ("HBG1", "HBG2"), "base_editing", None, (), "unknown", "clinical",
        (_registry_ref(
            "NCT06678165", "trial_registration", NEVER_REVIEWED,
            "Registry identifier awaiting source review."
        ),
         _registry_ref(
             "ChiCTR2500112517", "trial_registration", NEVER_REVIEWED,
             "Registry identifier awaiting source review."
         ),
         _registry_ref(
             "YOLT-204 public program materials", "sponsor_material", NEVER_REVIEWED,
             "A dated source has not been reviewed for this entry."
         ),
         _registry_ref(
             "YOLT-204 preclinical publication",
             "preclinical_paper", NEVER_REVIEWED,
             "No PMID or DOI has been reviewed for this entry."
         )),
        "Reported in vivo base-editing program at the HBG regulatory region. The sources have not "
        "been reviewed, so the editor class, exact position, motif, participant details, doses, "
        "and outcomes remain unasserted.",
        program_status="active",
    ),
]

# Programs targeting the BCL11A erythroid enhancer, rather than an HBG promoter, fail the registry's
# target-first criterion: same therapeutic goal, different locus. Recorded here rather than deleted,
# because "we considered it and it does not qualify" is information.
MOVED_OUT_OF_REGISTRY = (
    {"label": "exa-cel / CTX001 (Casgevy)",
     "reason": "BCL11A erythroid enhancer, not the HBG1/HBG2 promoter",
     "note": "Useful as a goal-level comparator. It must not count toward how complete the "
             "promoter-program list is, which is what it was doing before."},
    {"label": "BRL-101 (Bioraylab / 邦耀生物)",
     "reason": "BCL11A erythroid enhancer, not the HBG1/HBG2 promoter",
     "citations": (_registry_ref(
         "Bioraylab BRL-101 sponsor/trial materials (unopened)", "sponsor_or_trial_material",
         NEVER_REVIEWED, "Unverified; needs source review and a dated, specific source."
     ),),
     "note": "Kept outside the promoter registry pending source review."},
    {"label": "ET-01 (EdiGene / 博雅辑因)",
     "reason": "BCL11A erythroid enhancer, not the HBG1/HBG2 promoter",
     "citations": (_registry_ref(
         "EdiGene ET-01 sponsor/trial materials (unopened)", "sponsor_or_trial_material",
         NEVER_REVIEWED, "Unverified; needs source review and a dated, specific source."
     ),),
     "note": "Kept outside the promoter registry pending source review."},
)

REGISTRY = _hpfh_entries() + _CLINICAL
_BY_ID = {e.entry_id: e for e in REGISTRY}


def entry(entry_id: str) -> RegistryEntry | None:
    return _BY_ID.get(entry_id)


def relationship(gene: str, promoter_offset: int, editor_class: str, mechanism: str | None = None) -> dict:
    """Compare a candidate against the list.

    Returns the closest structural match, the strongest evidence across everything it matched, and
    the individual matches.

    exact_edit: a base-editing entry at the same promoter offset with the same editor class.
    same_motif: same motif (BCL11A, say) at an HBG-promoter entry, but not the same edit.
    nearby:     within 10 bp of an entry's position.
    """
    offset = abs(promoter_offset)
    motif = _motif_for(mechanism) if mechanism else None
    matches: list[tuple[str, RegistryEntry]] = []
    for e in REGISTRY:
        if offset in e.offsets and e.modality == "base_editing" and e.editor_class == editor_class:
            matches.append((EXACT_EDIT, e))
        elif motif and motif != "unknown" and motif == e.motif and "HBG1" in e.genes:
            matches.append((SAME_MOTIF, e))
        elif e.offsets and any(abs(o - offset) <= 10 for o in e.offsets):
            matches.append((NEARBY, e))

    levels = {lvl for lvl, _ in matches}
    if EXACT_EDIT in levels:
        rel = EXACT_EDIT
    elif SAME_MOTIF in levels:
        rel = SAME_MOTIF
    elif matches:
        rel = NEARBY
    else:
        rel = NO_MATCH

    strengths = [e.strength for _, e in matches] or ["none"]
    strength = max(strengths, key=STRENGTH_ORDER.index)

    return {
        "relationship": rel,
        "evidence_strength": strength,
        "matched_entries": [{"entry_id": e.entry_id, "label": e.label, "relationship": lvl,
                             "strength": e.strength, "mechanism_evidence": e.mechanism_evidence,
                             "citations": list(e.citations)}
                            for lvl, e in matches],
        "note": "Clinical programs do not always publish their exact position, so a motif match is "
                "an analog, not the same edit.",
    }
