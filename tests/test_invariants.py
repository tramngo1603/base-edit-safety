"""Cross-table consistency checks for allele keys, sources, and coordinates.

The allele-coverage check verifies that distinct alleles cannot overwrite one another. The source
access limit is enforced by
`test_a_claim_may_not_assert_more_than_its_access_allows`, which compares each claim's `audit_status`
with the recorded source access level.
"""
from __future__ import annotations

import dataclasses

import pytest

from base_edit_safety import anchors, bystanders, hbg_registry
from base_edit_safety.bystanders import ABSENT, PRESENT, UNKNOWN, UNRESOLVED, bystander_function


# --- Keep more than one DNA change at the same offset separate ----------------------------------

def test_minus_110_resolves_to_exactly_one_entry_per_allele():
    """-110 A>G is an ABE bystander; -110 A>C is an out-of-scope NF-Y transversion.

    An offset-only dictionary would overwrite one allele. The allele key must keep both distinct.
    """
    assert bystanders.rule_for(110, "A", "G") is not None
    assert bystanders.rule_for(110, "A", "C") is None

    excluded = [e for e in anchors.OUT_OF_SCOPE if e.promoter_offset == 110]
    assert len(excluded) == 1
    assert (excluded[0].promoter_sense_ref, excluded[0].promoter_sense_alt) == ("A", "C")

    # The two records are in different tables and describe different edits.
    assert not any(e.promoter_sense_alt == "G" for e in excluded)


def test_no_allele_appears_in_both_in_scope_and_out_of_scope():
    in_scope = {(p.promoter_offset, p.promoter_sense_ref, p.promoter_sense_alt)
                for p in anchors.IN_SCOPE}
    excluded = {(e.promoter_offset, e.promoter_sense_ref, e.promoter_sense_alt)
                for e in anchors.OUT_OF_SCOPE}
    assert in_scope & excluded == set()


def test_every_table_is_keyed_by_allele_not_offset():
    """Require every table to preserve a second allele at the same offset."""
    keys = [r.key for r in bystanders.RULES]
    assert len(keys) == len(set(keys)), "duplicate (offset, ref, alt) in the bystander table"
    # -114 carries C>T here and C>A / C>G in the excluded table: three alleles, one offset
    offsets = [k[0] for k in keys] + [e.promoter_offset for e in anchors.OUT_OF_SCOPE]
    assert offsets.count(-114) + offsets.count(114) >= 3


# --- Apply a conditional rule only when its required edit is present -----------------------------

def test_a_conditional_rule_never_fires_without_its_condition():
    for rule in bystanders.RULES:
        if not rule.requires_present:
            continue
        cls, note = bystander_function("HBG1", rule.offset, rule.ref, rule.alt, {})
        assert cls == UNKNOWN, (
            f"{rule.offset} {rule.ref}>{rule.alt} returned {cls} with none of its required edits "
            f"present; it must decline instead"
        )


def test_a_conditional_rule_declines_when_the_condition_is_unresolved():
    for rule in bystanders.RULES:
        for required in rule.requires_present:
            ctx = {required: UNRESOLVED}
            cls, note = bystander_function("HBG1", rule.offset, rule.ref, rule.alt, ctx)
            assert cls == UNKNOWN and "unresolved" in note


def test_requires_absent_is_honoured():
    """-181 is neutral only while -173 is unedited. With both edited it must not claim neutral."""
    rule = bystanders.rule_for(181, "A", "G")
    assert rule.requires_absent, "the -181 rule lost its requires_absent condition"
    ctx = {e: PRESENT for e in rule.requires_present}
    assert bystander_function("HBG1", 181, "A", "G", ctx)[0] == bystanders.NEUTRAL
    ctx.update({e: PRESENT for e in rule.requires_absent})
    assert bystander_function("HBG1", 181, "A", "G", ctx)[0] == UNKNOWN


def test_every_rule_states_a_confidence_and_a_source():
    for rule in bystanders.RULES:
        assert rule.confidence, f"{rule.offset} {rule.ref}>{rule.alt} has no confidence"
        assert rule.note, f"{rule.offset} {rule.ref}>{rule.alt} has no note"
        assert rule.source_pmid, f"{rule.offset} {rule.ref}>{rule.alt} has no source PMID"


def test_the_medium_confidence_calls_stay_medium():
    """Keep motif-only adverse classifications at their recorded medium confidence."""
    for offset in (112, 110):
        rule = bystanders.rule_for(offset, "A", "G")
        assert "medium" in rule.confidence, (
            f"-{offset} adverse uses motif structure and has not been tested with -113A>G; "
            f"the functional class must remain medium confidence"
        )
    # Supplementary Table 12 supplies the guide used to derive the -173 and -181 coordinates. Keep
    # the coordinate method and published sequence in each rule.
    for offset in (173, 181):
        rule = bystanders.rule_for(offset, "A", "G")
        assert "derived directly" in rule.confidence, f"-{offset} lacks its coordinate method"
        assert any("AGATATTTGCATTGAGATAG" in g for g in rule.observed_in_guide), (
            f"-{offset} must cite the published spacer it was derived from"
        )


def test_110_records_the_no_change_result_for_isolated_110():
    """Retain the reported no-change result for -110A>G tested alone."""
    rule = bystanders.rule_for(110, "A", "G")
    assert "measured no HbF change" in rule.note
    assert "34341563" in rule.source_pmid


# --- Derive a shared entry ID only when both HBG genes have the same status ----------------------

def test_entry_id_prefix_is_only_derived_where_the_two_statuses_agree():
    for pos in anchors.IN_SCOPE:
        both_synthetic = (pos.hbg1_status == anchors.SYNTHETIC,
                          pos.hbg2_status == anchors.SYNTHETIC)
        assert both_synthetic[0] == both_synthetic[1], (
            f"offset {pos.promoter_offset} is synthetic on one paralog and not the other, so one "
            f"entry_id cannot name both"
        )


def test_a_divergent_row_raises_rather_than_picking_one_paralog():
    """Reject one shared identifier when the two HBG genes have different statuses."""
    divergent = dataclasses.replace(anchors.IN_SCOPE[0],
                                    hbg1_status=anchors.SYNTHETIC,
                                    hbg2_status=anchors.DB)
    with pytest.raises(ValueError, match="cannot describe both paralogs"):
        hbg_registry._entry_id(divergent)


def test_synthetic_positions_are_named_synthetic():
    ids = {e.entry_id for e in hbg_registry.REGISTRY}
    for offset in (123, 124, 122):
        assert f"synthetic_-{offset}_" in " ".join(ids), offset


# --- Check the gene-name mapping and its source-review state -------------------------------------

def test_the_gamma_naming_convention_cites_hgnc_and_marks_review_pending():
    conv = anchors.GAMMA_GENE_NAMING
    assert conv["Agamma"] == "HBG1" and conv["Ggamma"] == "HBG2"
    assert "HGNC" in conv["status"]
    assert "review pending" in conv["status"]


# --- Citation-to-allele consistency ---------------------------------------------------------------

def test_no_citation_references_an_allele_absent_from_the_tables():
    known = {(p.promoter_offset, p.promoter_sense_ref, p.promoter_sense_alt)
             for p in anchors.IN_SCOPE}
    known |= {(-r.offset, r.ref, r.alt) for r in bystanders.RULES}
    for key in hbg_registry._CITATIONS:
        assert isinstance(key, tuple) and len(key) == 3, (
            f"_CITATIONS is still keyed by {key!r}; it must be keyed by allele so that evidence "
            f"about -114 C>A does not attach to a row about -114 C>T"
        )
        assert key in known


# --- Study-system structure, engraftment evidence, and disrupted motifs --------------------------

def test_architecture_is_recorded_on_the_claim_not_the_paper():
    """Store the experimental system on each claim when one paper covers two systems.

    Wienert 2015 is the case that forces this: MEL with a human BAC, and K562 with the Ggamma body
    deleted. Half its claims come from each.
    """
    wienert = [c for claims in hbg_registry._CITATIONS.values() for c in claims
               if c["pmid"] == "25971621"]
    systems = {c["locus_architecture"]["system"] for c in wienert}
    assert len(systems) >= 2, f"25971621 collapsed to one architecture: {systems}"


def test_no_mechanism_paper_used_an_intact_native_locus():
    """Mechanism and editing evidence come from structurally different systems."""
    for claims in hbg_registry._CITATIONS.values():
        for c in claims:
            if not c["role"].startswith("mechanism"):
                continue
            arch = c.get("locus_architecture")
            if arch is None:
                continue
            assert not arch["native"], (
                f"PMID {c['pmid']} is recorded as a mechanism claim in a native locus; every "
                f"mechanism paper checked so far used an engineered one"
            )


def test_a_mechanism_cannot_be_projected_onto_a_paralog_that_was_not_intact():
    # Wienert's K562 kept HBG1 and deleted HBG2; the MEL BAC replaced Agamma coding.
    assert hbg_registry.mechanism_projectable_onto("HBG1", 175, "T", "C") is True
    assert hbg_registry.mechanism_projectable_onto("HBG2", 175, "T", "C") is False
    # Martyn's HUDEP-2 was Ggamma-deleted, so the GATA1 mechanism has no HBG2 evidence at all.
    assert hbg_registry.mechanism_projectable_onto("HBG2", 113, "A", "G") is False


def test_the_transplant_axis_separates_two_rows_that_share_a_tier():
    """Separate transplantation evidence for two entries with the same general evidence level."""
    assert hbg_registry.durability_for(124, "T", "C") == hbg_registry.NO_ENGRAFTMENT_DATA
    assert hbg_registry.durability_for(113, "A", "G") == hbg_registry.XENOTRANSPLANT
    by_id = {e.entry_id: e for e in hbg_registry.REGISTRY}
    assert by_id["synthetic_-124_ABE"].strength == by_id["hpfh_-113_ABE"].strength == "published_hspc"


def test_lost_motif_defaults_to_not_reviewed_never_to_absent():
    for pos in anchors.IN_SCOPE:
        assert pos.lost_motif, f"-{pos.promoter_offset} has an empty lost_motif"
    unchecked = [p.promoter_offset for p in anchors.IN_SCOPE
                 if p.lost_motif == anchors.LOST_MOTIF_NOT_REVIEWED]
    assert unchecked, "every position claims a reviewed lost_motif despite incomplete source review"
    # the two confirmed cases, both found in a supplement rather than a main text
    by_offset = {p.promoter_offset: p for p in anchors.IN_SCOPE}
    assert "GATA1" in by_offset[175].lost_motif
    assert "NF-Y" in by_offset[113].lost_motif


def test_a_claim_may_not_assert_more_than_its_access_allows():
    for claims in hbg_registry._CITATIONS.values():
        for c in claims:
            if c["audit_status"] == hbg_registry.NEVER_REVIEWED:
                assert c["access_level"] == hbg_registry.NOT_ACCESSED, c
            if c["audit_status"] == hbg_registry.FULL_TEXT_READ:
                assert c["access_level"] == hbg_registry.FULL_TEXT, c
            assert c["access_level"] in hbg_registry.ACCESS_ORDER, c


def test_exa_cel_is_out_of_the_registry_and_the_reason_is_recorded():
    assert not any("bcl11a_enhancer" == e.entry_id for e in hbg_registry.REGISTRY)
    assert any("exa-cel" in m["label"] for m in hbg_registry.MOVED_OUT_OF_REGISTRY)
    # six correct promoter programs, with enhancer programs kept outside the registry
    assert len(hbg_registry._CLINICAL) == 6


# --- the published spacers, located in the reference ----------------------------------------------

def _window():
    import json, pathlib
    cache = pathlib.Path("data/inputs/reference_cache/hg38_chr11_5249937_5250067.json")
    if not cache.is_file():
        pytest.skip("requires the ignored GRCh38 reference-window cache")
    w = json.load(open(cache))
    return w["plus_sequence"], w["start_1based"]


def _revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def test_every_published_spacer_is_where_the_table_says_it_is():
    """Derive recorded offsets from each Mayuranathan Supplementary Table 12 guide sequence.

    This source resolved the earlier interpolation at -173 and -181. A coordinate change must still
    reproduce every recorded offset from the sequence.
    """
    plus, lo = _window()
    for name, g in bystanders.PUBLISHED_SPACERS.items():
        spacer = g["spacer"]
        idx = plus.find(spacer) if g["strand"] == "+" else plus.find(_revcomp(spacer))
        assert idx >= 0, f"{name} not found in the HBG1 window"
        if g["strand"] == "+":
            five_prime = lo + idx
            coord = lambda p: five_prime + (p - 1)
        else:
            five_prime = lo + idx + len(spacer) - 1
            coord = lambda p: five_prime - (p - 1)
        assert five_prime == g["five_prime_hbg1"], f"{name}: 5' end moved"
        for pos, expected in g["offsets"].items():
            assert spacer[pos - 1] == "A", f"{name} position {pos} is not an A"
            assert coord(pos) - anchors.HBG1_TSS == -expected, (
                f"{name} A{pos} maps to {coord(pos) - anchors.HBG1_TSS}, table says {-expected}"
            )


def test_the_repos_113_protospacer_is_mayuranathans_published_one():
    """Match the stored -113 guide sequence and coordinate to the published guide."""
    published = bystanders.PUBLISHED_SPACERS["sgRNA -113"]
    assert published["spacer"] == "CTTGACCAATAGCCTTGACA"
    assert published["five_prime_hbg1"] == 5_249_977
    # the four positions Han and Ravi both report are exactly this spacer's editable A's
    assert set(published["offsets"].values()) == {-116, -113, -112, -110}


def test_the_175_spacer_spans_the_catastrophic_positions_too():
    """Require the -175 guide to report its overlap with -187 and -189."""
    offsets = bystanders.PUBLISHED_SPACERS["sgRNA -175"]["offsets"]
    assert offsets[17] == -187 and offsets[19] == -189
    catastrophic = {c.offset for c in bystanders.CATASTROPHIC}
    assert {-187, -189} <= catastrophic


# --- guide identity, established from the two papers' supplements ---------------------------------

def test_our_guides_are_the_guides_the_papers_measured():
    """Match 19 stored guide sequences to Ravi's screened set, including both preferred designs.

    A measured rate can be linked to a design only while its guide sequence matches the publication.
    """
    ident = bystanders.GUIDE_IDENTITIES
    assert ident["AGGCAAGGCTGGCCAACCCA"]["ours"] == "60forw"      # PREFERRED tier
    assert ident["AGGCAAGGCTGGCCAACCCA"]["ravi"] == "gRNA 10"
    assert ident["GGCAAGGCTGGCCAACCCAT"]["ours"] == "61forw"
    # the -113 guide is the same molecule in three places
    m113 = ident["CTTGACCAATAGCCTTGACA"]
    assert m113["ours"] == "19rev" and m113["ravi"] == "gRNA 2"
    assert m113["mayuranathan"] == "sgRNA -113"
    assert bystanders.PUBLISHED_SPACERS["sgRNA -113"]["spacer"] == "CTTGACCAATAGCCTTGACA"


def test_the_catastrophic_check_fires_on_the_guide_that_demonstrated_it():
    """Match 122forw to Mayuranathan's sgRNA -187/189 and report its high-impact positions.

    The guide remains in the enumerated set, so the overlap check must flag it.
    """
    from base_edit_safety.pipeline.datatypes import Guide, EDIT_WINDOW
    guide = Guide(guide_id="HBG1_122forw_ABE8e", protospacer="ATTGAGATAGTGTGGGGAAG", pam="NGG",
                  strand="+", genome_build="GRCh38", chrom="chr11", gene="HBG1", gene_strand="-",
                  tss_1based=anchors.HBG1_TSS, editor_class="ABE", editor_name="ABE8e",
                  edit_from="A", edit_to="G", proto_5p_genomic=5_250_038)
    exposure = bystanders.catastrophic_exposure(guide, EDIT_WINDOW)
    assert exposure == [(7, -187), (9, -189)], exposure


def test_a_clear_guide_reports_no_exposure():
    from base_edit_safety.pipeline.datatypes import Guide, EDIT_WINDOW
    guide = Guide(guide_id="HBG1_19rev_ABE8e", protospacer="CTTGACCAATAGCCTTGACA", pam="AGG",
                  strand="-", genome_build="GRCh38", chrom="chr11", gene="HBG1", gene_strand="-",
                  tss_1based=anchors.HBG1_TSS, editor_class="ABE", editor_name="ABE8e",
                  edit_from="A", edit_to="G", proto_5p_genomic=5_249_977)
    assert bystanders.catastrophic_exposure(guide, EDIT_WINDOW) == []


def test_the_198_spacer_preserves_unresolved_nonblocking_provenance():
    published = bystanders.PUBLISHED_SPACERS["sgRNA -198"]
    published_spacer = "GTGGGGAAGGGGCCCCCAA"
    reference_protospacer = "GTGGGGAAGGGGCCCCCAAG"
    note = published.get("note", "")

    assert published["spacer"] == published_spacer
    assert len(published["spacer"]) == 19
    assert published["pam"] == "NGG"
    assert published["strand"] == "+"
    assert published["five_prime_hbg1"] == 5_250_049
    assert published["offsets"] == {7: -198}

    assert reference_protospacer == published_spacer + "G"
    assert reference_protospacer in note
    assert "followed by AGG" in note
    assert "gRNA 4" in note

    assert "do not establish whether" in note
    assert "almost certainly" not in note
    assert "dropped it in transcription" not in note

    assert "Ensembl and UCSC" in note
    assert "GRCh38" in note and "v0.1" in note
    assert "nonblocking" in note
