"""Compare candidates with documented HBG1/HBG2 edits and programs.

The registry contains natural hereditary-persistence-of-fetal-hemoglobin (HPFH) promoter variants,
published HBG-promoter editing in hematopoietic stem and progenitor cells (HSPCs), and clinical
HBG1/HBG2 programs.

A candidate receives two separate labels. The relationship label is `exact_edit`, `same_motif`,
`nearby`, or `no_match`. The evidence label states the support available for the matched entry.
Keeping them separate prevents an exact published HSPC edit from being reported as a nearby match.

When an accessible source does not give an exact target position or motif, those fields remain
unknown. Each citation records whether its source has been reviewed.
"""
from __future__ import annotations

from dataclasses import dataclass

from . import anchors

# Ways a candidate can match an entry in this registry.
EXACT_EDIT = "exact_edit"
SAME_MOTIF = "same_motif"
NEARBY = "nearby"
NO_MATCH = "no_match"

# Evidence levels, ordered from least to most support. In `evidence.py`, the separate
# `empirically_measured` label requires measurements from the current run.
STRENGTH_ORDER = ("none", "analog", "motif", "published_hspc", "clinical")

# Evidence for a biological mechanism is separate from evidence that a variant exists. A mechanism
# remains `NOT_REVIEWED` until its supporting paper has been checked.
#   in_vitro_binding_only  binding tests outside cells, without a cell-based confirmation
#   binding_plus_cellular  binding tests plus a cell-based test that excludes other explanations
#   established            supported by the reviewed evidence defined for this registry
#   not_reviewed           the supporting paper has not been checked for this claim
IN_VITRO_BINDING_ONLY = "in_vitro_binding_only"
BINDING_PLUS_CELLULAR = "binding_plus_cellular"
ESTABLISHED = "established"
NOT_REVIEWED = "not_reviewed"

# Whether the cited paper has been reviewed for this registry.
FULL_TEXT_READ = "full_text_read"
NEVER_REVIEWED = "never_reviewed"


def _cite(pmid: str, role: str, audit: str = NOT_REVIEWED, note: str | None = None) -> dict:
    """Return a citation with its claim role and review status.

    `role` states the supported claim type so clinical context is not used as a substitute for a
    variant description, phenotype comparison, or source review.
    """
    entry = {"pmid": pmid, "role": role, "audit_status": audit}
    if note:
        entry["note"] = note
    return entry


@dataclass(frozen=True)
class RegistryEntry:
    entry_id: str
    label: str
    genes: tuple
    modality: str                 # `base_editing` or `nuclease`
    editor_class: str | None      # `ABE`, `CBE`, or `None` when the class is not confirmed
    offsets: tuple                # positive promoter offsets; empty when no exact position is known
    motif: str                    # named DNA-binding site, such as BCL11A or GATA1
    strength: str                 # one of the values in `STRENGTH_ORDER`
    citations: tuple              # source records created by `_cite`
    note: str
    mechanism_evidence: str = NOT_REVIEWED
    program_status: str | None = None   # reported program state for clinical entries


# --- Experimental system for each claim ----------------------------------------------------------
#
# Store the experimental system on each claim, not once per paper. Wienert 2015 reports two systems.
# Mayuranathan and Ravi each report an engineered cell line and unmodified primary cells.
#
# `intact_paralogs` records which HBG gene copies were present and unchanged. A mechanism studied
# after HBG2 deletion cannot be assigned to HBG2. None of the five mechanism papers used an
# unchanged, naturally occurring HBG region.

@dataclass(frozen=True)
class LocusArchitecture:
    system: str
    description: str
    intact_paralogs: tuple          # HBG genes that were present and unchanged
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

# Source access level. A claim cannot exceed what the accessed material supports.
ACCESS_ORDER = ("not_accessed", "abstract_only", "database_record", "full_text")
FULL_TEXT, ABSTRACT_ONLY, NOT_ACCESSED = "full_text", "abstract_only", "not_accessed"

# Engraftment evidence is separate from general evidence strength. Ravi and Mayuranathan both have
# `published_hspc` strength, but only one includes transplantation.
NO_ENGRAFTMENT_DATA = "no_engraftment_data"
XENOTRANSPLANT = "xenotransplant_engraftment_shown"


def _claim(pmid: str, role: str, access: str = NOT_ACCESSED, audit: str = NOT_REVIEWED,
           architecture: LocusArchitecture | None = None, engraftment: bool = False,
           note: str | None = None) -> dict:
    """Return one assertion with its study system, source access, and review status."""
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


# --- Claims for each DNA change -----------------------------------------------------------------
#
# Use `(offset, reference base, alternate base)` as the key. At -114, different sources study C>T,
# C>A, and C>G. An offset-only key would attach evidence for one DNA change to another.

# Doerfler supports two positions on each dependency branch, rather than all four positions equally.
_DOERFLER_TESTED = lambda branch: _claim(
    "34341563", f"qualified_variant_specific_{branch}", FULL_TEXT, FULL_TEXT_READ, DOERFLER_FUSED,
    note="Retained only for the tested dependency branch, with the fusion-promoter limitation. "
         "Not proof of KLF1 binding, KLF1-ZBTB7A competition, or direct GATA1-versus-BCL11A "
         "displacement.")
_DOERFLER_FRAMEWORK = _claim(
    "34341563", "competition_framework_only", FULL_TEXT, FULL_TEXT_READ, DOERFLER_FUSED,
    note="Framework only at this position, recorded once at the disrupted-BCL11A-motif class.")

_CITATIONS = {
    (198, "T", "C"): (
        _claim("2430647", "human_haplotype_report", FULL_TEXT, FULL_TEXT_READ,
               note="One reported chromosome also carries -158 and +210. The source does not "
                    "attribute the reported phenotype to -198 alone."),
        _claim("28659276", "mechanism_klf1", NOT_ACCESSED, NOT_REVIEWED, HUDEP2_DELTA_GGAMMA,
               note="KLF1 was nominated on erythroid abundance, not binding exclusivity; twelve "
                    "other KLFs bound the mutant probe."),
        _claim("37400614", "hspc_base_editing", NOT_ACCESSED, NOT_REVIEWED,
               PRIMARY_CD34_TRANSPLANT, engraftment=True),
        _DOERFLER_TESTED("upstream_gata1_dependent"),
    ),
    (175, "T", "C"): (
        _claim("2449926", "human_variant_report_hbg2", FULL_TEXT, FULL_TEXT_READ,
               note="One reported chromosome also carries -158 and changes after the gene. This "
                    "source proposes an octamer mechanism, not TAL1; PMID 25971621 supports TAL1."),
        _claim("25971621", "mechanism_tal1_ebox", FULL_TEXT, FULL_TEXT_READ, WIENERT_MEL_BAC,
               note="The gamma coding sequence is replaced by a reporter, so the reported "
                    "gamma-globin expression value is fluorescence rather than gamma-globin or "
                    "HbF. The MEL ChIP amplicon does not establish which HBG promoter contributed "
                    "the signal. It therefore supports increased promoter occupancy but not a "
                    "specific assignment to HBG1 -175T>C."),
        _claim("25971621", "mechanism_tal1_ebox_second_system", FULL_TEXT, FULL_TEXT_READ,
               WIENERT_K562,
               note="Second architecture in the same paper. Same reporter caveat."),
        _claim("37400614", "hspc_base_editing", NOT_ACCESSED, NOT_REVIEWED,
               PRIMARY_CD34_TRANSPLANT, engraftment=True),
        _claim("7679879", "human_variant_report_hbg1", NOT_ACCESSED, NEVER_REVIEWED,
               note="Coleman 1993 reports HBG1 -175 on the same chromosome as HBG2 -158. Its phase "
                    "claim comes from a compilation; the six primary sources remain unreviewed."),
    ),
    (124, "T", "C"): (
        _claim("35147495", "synthetic_screen_hit", FULL_TEXT, FULL_TEXT_READ, HUDEP2_INTACT,
               note="Discovery screen was HUDEP-2, not HSPCs."),
        _claim("35147495", "allele_installed_in_hspcs", FULL_TEXT, FULL_TEXT_READ, PRIMARY_CD34,
               note="The exact allele was installed and genotyped in primary CD34+ HSPCs at day "
                    "6, before differentiation. The study did not include transplantation, so it "
                    "does not establish durable editing in repopulating HSCs."),
    ),
    (123, "T", "C"): (
        _claim("35147495", "synthetic_screen_hit", FULL_TEXT, FULL_TEXT_READ, HUDEP2_INTACT),
        _claim("35147495", "allele_installed_in_hspcs", FULL_TEXT, FULL_TEXT_READ, PRIMARY_CD34,
               note="The same paper reports no transplantation data for this edit. Its KLF1 "
                    "binding test uses the combined -123/-124 change. The ChIP control changed by "
                    "a similar amount, and the source names other possible factors."),
    ),
    (122, "G", "A"): (
        _claim("35147495", "cbe_on_target_conversion", FULL_TEXT, FULL_TEXT_READ, HUDEP2_INTACT,
               note="Stated in both Results and Discussion as the on-target site for gRNA-10/-11. "
                    "CBE work in this paper was HUDEP-2, not HSPCs."),
    ),
    (117, "G", "A"): (
        _claim("2417646", "phenotype_concordance", FULL_TEXT, FULL_TEXT_READ,
               note="This source compares phenotypes but is not the original variant description. "
                    "It credits Gelinas 1985 and Collins 1985, which are not included in this "
                    "registry. This entry therefore has no reviewed existence source for -117."),
        _claim("3181130", "mechanism_distal_ccaat", NOT_ACCESSED, NOT_REVIEWED),
        _claim("37989316", "motif_base_editing", NOT_ACCESSED, NOT_REVIEWED, PRIMARY_CD34),
        _DOERFLER_FRAMEWORK,
    ),
    (114, "C", "T"): (
        _claim("1698280", "human_variant_report_hbg2", FULL_TEXT, FULL_TEXT_READ,
               note="This report covers HBG2 in two related cases. HBG1 was not sequenced or used "
                    "in an assay, so the source does not support an HBG1 -114 claim."),
        _claim("1704803", "correspondence_n2_siblings", FULL_TEXT, FULL_TEXT_READ,
               note="Two brothers, published as correspondence. Its paralog claim is an explicit "
                    "assumption in the paper's own words, not a measurement."),
        _claim("35147495", "regional_base_editing", FULL_TEXT, FULL_TEXT_READ, HUDEP2_INTACT,
               note="The -114 editing experiment used HUDEP-2 cells. The paper did not apply a CBE "
                    "to primary HSPCs, so it does not support an HSPC base-editing label for this "
                    "DNA change."),
        _claim("37989316", "motif_base_editing", NOT_ACCESSED, NOT_REVIEWED, PRIMARY_CD34),
        _claim("34341563", "competition_framework_only", FULL_TEXT, FULL_TEXT_READ, DOERFLER_FUSED,
               note="Framework only. Doerfler's cellular change is C>A and its EMSA probe is C>G; "
                    "neither is the C>T change stored here. The CP1 mechanism rests on EMSA competition "
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
    """Return whether a cited experiment reported engraftment, separately from evidence level."""
    claims = _CITATIONS.get((abs(offset), ref, alt), ())
    return XENOTRANSPLANT if any(c.get("demonstrates_engraftment") for c in claims) \
        else NO_ENGRAFTMENT_DATA


def mechanism_projectable_onto(gene: str, offset: int, ref: str, alt: str) -> bool:
    """Return whether a mechanism was studied with the named HBG gene present and unchanged.

    A HUDEP-2 experiment performed after deletion of G-gamma, which is HBG2, cannot support an HBG2
    mechanism claim.
    """
    for c in _CITATIONS.get((abs(offset), ref, alt), ()):
        if not c["role"].startswith("mechanism"):
            continue
        arch = c.get("locus_architecture")
        if arch and gene in arch["intact_paralogs"]:
            return True
    return False


# Support for the mechanism assigned to each position. Three supporting papers have been reviewed;
# each provides binding tests outside cells but no cell-based confirmation.
#   -198 KLF1     KLF1 was nominated on erythroid abundance, not on binding exclusivity, and twelve
#                 other KLFs bound the mutant probe.
#   -123/-124     EMSA supports binding, but the ChIP negative control changed similarly and the authors
#                 name alternative factors themselves.
#   -114 CP1      EMSA competition only, with no antibody and no supershift, and DNaseI footprinting
#                 failed to resolve the distal motif in either probe.
# The other papers have not been reviewed, so their mechanism status remains `NOT_REVIEWED`.
_MECHANISM_EVIDENCE = {
    198: IN_VITRO_BINDING_ONLY,
    124: IN_VITRO_BINDING_ONLY,
    123: IN_VITRO_BINDING_ONLY,
    114: IN_VITRO_BINDING_ONLY,
    175: NOT_REVIEWED,
    117: NOT_REVIEWED,
    113: NOT_REVIEWED,
}

# Evidence strength per HPFH position (see `docs/coordinates.md`).
_STRENGTH = {
    198: "published_hspc",   # HSPC ABE (Mayuranathan 2023)
    175: "published_hspc",   # highest-support exact-allele HSPC ABE entry (Mayuranathan 2023)
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
    """Use `hpfh_` for natural HPFH alleles and `synthetic_` for screen-derived changes."""
    is_synthetic = {pos.hbg1_status == anchors.SYNTHETIC, pos.hbg2_status == anchors.SYNTHETIC}
    if len(is_synthetic) > 1:
        # One entry cannot cover both genes if their natural-versus-synthetic status differs. Use a
        # separate, gene-specific identifier before adding such a row.
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


# Reported clinical HBG1/HBG2 promoter programs. Unknown fields remain unset until source review.
_CLINICAL = [
    RegistryEntry(
        "beam_101", "risto-cel, formerly BEAM-101 (Beam Therapeutics, BEACON)", ("HBG1", "HBG2"),
        "base_editing", None, (), "BCL11A", "clinical",
        (_registry_ref("NCT05456880", "trial_registration", note="BEACON"),
         _registry_ref("Beam Therapeutics BEAM-101 program materials (public)", "sponsor_material",
                       note="Not resolvable as written; needs a dated, specific source.")),
        "Reported CD34+ HSPC base-editing program at the HBG1/HBG2 promoter BCL11A site. The "
        "stable entry identifier retains the former program name. The editor class remains unset "
        "until a specific sponsor source has been reviewed.",
        program_status="active",
    ),
    RegistryEntry(
        "reni_cel_edit301", "reni-cel / EDIT-301 (Editas, RUBY / EdiThal)", ("HBG1", "HBG2"),
        "nuclease", None, (), "BCL11A", "clinical",
        (_registry_ref("NCT04853576", "trial_registration", note="RUBY"),
         _registry_ref("NCT05444894", "trial_registration", note="EdiThal"),
         _registry_ref("Editas reni-cel program materials (8-K summary)", "sponsor_material",
                       note="Not resolvable as written; needs the SEC accession and a date.")),
        "Reported AsCas12a nuclease program at the HBG1/HBG2 promoter BCL11A site. A nuclease can "
        "create insertions or deletions, so this entry is a same-location comparison rather than a "
        "base-editing comparison. The recorded program state does not state a safety or efficacy "
        "result.",
        program_status="discontinued_2024_12",
    ),
    # The sources for these entries have not been reviewed. Their identifiers therefore retain the
    # `NEVER_REVIEWED` status and cannot support more detailed fields.
    RegistryEntry(
        "cs_101", "CS-101 / CS-206 (CorrectSequence Therapeutics)", ("HBG1", "HBG2"),
        "base_editing", None, (), "BCL11A", "clinical",
        (_registry_ref("NCT06024876", "trial_registration", NEVER_REVIEWED),
         _registry_ref("NCT06291961", "trial_registration", NEVER_REVIEWED),
         _registry_ref("NCT06565026", "trial_registration", NEVER_REVIEWED)),
        "Reported clinical base-editing program at the HBG1/HBG2 promoter BCL11A site. The exact "
        "editor class remains unset because the listed registry sources have not been reviewed.",
        program_status="active",
    ),
    RegistryEntry(
        "otq923", "OTQ923 (Novartis)", ("HBG1", "HBG2"),
        "nuclease", None, (), "BCL11A", "clinical",
        (_registry_ref("NCT04443907", "trial_registration", NEVER_REVIEWED),
         _cite("37646679", "clinical_trial_report", NEVER_REVIEWED,
               "Sharma et al., NEJM 2023, Cas9 disruption of the HBG1 and HBG2 promoters. Cited by "
               "PMID because this registry records no DOI for the source.")),
        "Reported Cas9 nuclease program that disrupts both HBG promoters. It is a same-location "
        "comparison, not a base-editing comparison.",
        program_status="discontinued",
    ),
    RegistryEntry(
        "rm_001", "RM-001 (Reforgene Therapeutics)",
        ("HBG1", "HBG2"), "nuclease", None, (), "unknown", "clinical",
        (_registry_ref(
            "RM-001 public program materials", "sponsor_material", NEVER_REVIEWED,
            "No dated source or trial registry identifier has been reviewed for this entry."
        ),),
        "Reported HBG1/HBG2 promoter nuclease program. No reviewed source supports an exact "
        "position or motif, so this entry cannot create a position or motif match.",
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

# Programs that target the BCL11A erythroid enhancer are outside this HBG-promoter registry. Each
# excluded entry records that location difference.
MOVED_OUT_OF_REGISTRY = (
    {"label": "exa-cel / CTX001 (Casgevy)",
     "reason": "BCL11A erythroid enhancer, not the HBG1/HBG2 promoter",
     "note": "Retained as a comparison with a different target location."},
    {"label": "BRL-101 (Bioraylab / 邦耀生物)",
     "reason": "BCL11A erythroid enhancer, not the HBG1/HBG2 promoter",
     "citations": (_registry_ref(
         "BRL-101 public program materials", "sponsor_or_trial_material",
         NEVER_REVIEWED, "No dated source has been reviewed for this entry."
     ),),
     "note": "Kept outside the promoter registry pending source review."},
    {"label": "ET-01 (EdiGene / 博雅辑因)",
     "reason": "BCL11A erythroid enhancer, not the HBG1/HBG2 promoter",
     "citations": (_registry_ref(
         "ET-01 public program materials", "sponsor_or_trial_material",
         NEVER_REVIEWED, "No dated source has been reviewed for this entry."
     ),),
     "note": "Kept outside the promoter registry pending source review."},
)

REGISTRY = _hpfh_entries() + _CLINICAL
_BY_ID = {e.entry_id: e for e in REGISTRY}


def entry(entry_id: str) -> RegistryEntry | None:
    return _BY_ID.get(entry_id)


def relationship(gene: str, promoter_offset: int, editor_class: str, mechanism: str | None = None) -> dict:
    """Compare a candidate against the list.

    Return the closest relationship, the highest evidence level among matched entries, and each
    individual match.

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
