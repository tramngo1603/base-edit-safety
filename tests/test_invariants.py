"""Cross-table invariants. Each one exists because it was violated at least once.

All five invariants from the implementation brief are covered here. The two that needed the citation
tables re-keyed by allele are live now that section 4 has landed: the allele-coverage check sits with
the other keying tests, and the verification ceiling is enforced by
`test_a_claim_may_not_assert_more_than_its_access_allows`, which compares each claim's `audit_status`
against how far into the source anyone actually got.
"""
from __future__ import annotations

import dataclasses

import pytest

from base_edit_safety import anchors, bystanders, hbg_registry
from base_edit_safety.bystanders import ABSENT, PRESENT, UNKNOWN, UNRESOLVED, bystander_function


# --- invariant: one offset can hold several alleles, and they must never merge ------------------

def test_minus_110_resolves_to_exactly_one_entry_per_allele():
    """-110 A>G is an ABE bystander; -110 A>C is an out-of-scope NF-Y transversion.

    Keyed by offset alone these collide and one silently wins. This is the case that forced the
    re-key, so it gets its own test.
    """
    assert bystanders.rule_for(110, "A", "G") is not None
    assert bystanders.rule_for(110, "A", "C") is None

    excluded = [e for e in anchors.OUT_OF_SCOPE if e.promoter_offset == 110]
    assert len(excluded) == 1
    assert (excluded[0].promoter_sense_ref, excluded[0].promoter_sense_alt) == ("A", "C")

    # the two live in different tables and describe different edits
    assert not any(e.promoter_sense_alt == "G" for e in excluded)


def test_no_allele_appears_in_both_in_scope_and_out_of_scope():
    in_scope = {(p.promoter_offset, p.promoter_sense_ref, p.promoter_sense_alt)
                for p in anchors.IN_SCOPE}
    excluded = {(e.promoter_offset, e.promoter_sense_ref, e.promoter_sense_alt)
                for e in anchors.OUT_OF_SCOPE}
    assert in_scope & excluded == set()


def test_every_table_is_keyed_by_allele_not_offset():
    """If any table still collapses to offset, a second allele there would be lost silently."""
    keys = [r.key for r in bystanders.RULES]
    assert len(keys) == len(set(keys)), "duplicate (offset, ref, alt) in the bystander table"
    # -114 carries C>T here and C>A / C>G in the excluded table: three alleles, one offset
    offsets = [k[0] for k in keys] + [e.promoter_offset for e in anchors.OUT_OF_SCOPE]
    assert offsets.count(-114) + offsets.count(114) >= 3


# --- invariant: a conditional rule may not be applied to a record lacking the condition ---------

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
        assert rule.source_pmid, f"{rule.offset} {rule.ref}>{rule.alt} cites nothing"


def test_the_medium_confidence_calls_stay_medium():
    """These have needed re-flagging twice across sessions. Pin them so they cannot flatten."""
    for offset in (112, 110):
        rule = bystanders.rule_for(offset, "A", "G")
        assert "MEDIUM" in rule.confidence, (
            f"-{offset} adverse rests on motif structure alone and has never been tested with "
            f"-113A>G; its class must not be recorded as high confidence"
        )
    # -173 and -181 were interpolated until Supp Table 12 supplied the spacer. They are now read off
    # the sequence, so they must NOT slide back to a hedge, and must cite the spacer.
    for offset in (173, 181):
        rule = bystanders.rule_for(offset, "A", "G")
        assert "CONFIRMED" in rule.confidence, f"-{offset} lost its confirmation"
        assert "interpolated" not in rule.confidence.replace("rather than interpolated", ""), offset
        assert any("AGATATTTGCATTGAGATAG" in g for g in rule.observed_in_guide), (
            f"-{offset} must cite the published spacer it was derived from"
        )


def test_110_records_its_contrary_isolated_null():
    """Doerfler tested -110A>G alone and found nothing. That must stay visible on the row."""
    rule = bystanders.rule_for(110, "A", "G")
    assert "null" in rule.note.lower() or "NULL" in rule.note
    assert "34341563" in rule.source_pmid


# --- invariant: an entry id may only be derived where both paralogs agree ------------------------

def test_entry_id_prefix_is_only_derived_where_the_two_statuses_agree():
    for pos in anchors.IN_SCOPE:
        both_synthetic = (pos.hbg1_status == anchors.SYNTHETIC,
                          pos.hbg2_status == anchors.SYNTHETIC)
        assert both_synthetic[0] == both_synthetic[1], (
            f"offset {pos.promoter_offset} is synthetic on one paralog and not the other, so one "
            f"entry_id cannot name both"
        )


def test_a_divergent_row_raises_rather_than_picking_one_paralog():
    """The old code read HBG1's status alone. A divergent row must now be loud."""
    divergent = dataclasses.replace(anchors.IN_SCOPE[0],
                                    hbg1_status=anchors.SYNTHETIC,
                                    hbg2_status=anchors.DB)
    with pytest.raises(ValueError, match="cannot describe both paralogs"):
        hbg_registry._entry_id(divergent)


def test_synthetic_positions_are_named_synthetic():
    ids = {e.entry_id for e in hbg_registry.REGISTRY}
    for offset in (123, 124, 122):
        assert f"synthetic_-{offset}_" in " ".join(ids), offset


# --- invariant: the naming convention is recorded, not implicit ----------------------------------

def test_the_gamma_naming_convention_is_stated_and_marked_uncited():
    conv = anchors.GAMMA_GENE_NAMING
    assert conv["Agamma"] == "HBG1" and conv["Ggamma"] == "HBG2"
    assert "no authority cited" in conv["status"], (
        "the Agamma/Ggamma mapping underpins every per-gene claim in the repo and has no source; "
        "that must stay visible until one is cited"
    )


# --- the two invariants that need section 4 first -------------------------------------------------

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


# --- section 4: architecture, durability, lost motif ---------------------------------------------

def test_architecture_is_recorded_on_the_claim_not_the_paper():
    """One PMID can span several systems, so architecture cannot live on the citation.

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
    """Ravi and Mayuranathan would both read published_hspc on materially different evidence."""
    assert hbg_registry.durability_for(124, "T", "C") == hbg_registry.NO_ENGRAFTMENT_DATA
    assert hbg_registry.durability_for(113, "A", "G") == hbg_registry.XENOTRANSPLANT
    by_id = {e.entry_id: e for e in hbg_registry.REGISTRY}
    assert by_id["synthetic_-124_ABE"].strength == by_id["hpfh_-113_ABE"].strength == "published_hspc"


def test_lost_motif_defaults_to_not_reviewed_never_to_absent():
    for pos in anchors.IN_SCOPE:
        assert pos.lost_motif, f"-{pos.promoter_offset} has an empty lost_motif"
    unchecked = [p.promoter_offset for p in anchors.IN_SCOPE
                 if p.lost_motif == anchors.LOST_MOTIF_NOT_REVIEWED]
    assert unchecked, "every position claims a checked lost_motif, which nobody has done"
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
    """Each spacer from Mayuranathan Supp Table 12, located in GRCh38, gives its recorded offsets.

    This is what moved -173 and -181 off "interpolated". If a coordinate ever drifts, the offsets
    stop deriving and these rows go back to being inferences.
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


def test_locally_configured_spacer_is_not_embedded_in_the_package():
    assert "sgRNA -113" not in bystanders.PUBLISHED_SPACERS
    assert all(row.get("ours") != "19rev" for row in bystanders.GUIDE_IDENTITIES.values())


def test_the_175_spacer_spans_the_catastrophic_positions_too():
    """sgRNA -175 reaches -187 and -189 at A17 and A19, which is worth knowing before using it."""
    offsets = bystanders.PUBLISHED_SPACERS["sgRNA -175"]["offsets"]
    assert offsets[17] == -187 and offsets[19] == -189
    catastrophic = {c.offset for c in bystanders.CATASTROPHIC}
    assert {-187, -189} <= catastrophic


# --- guide identity, established from the two papers' supplements ---------------------------------

def test_our_guides_are_the_guides_the_papers_measured():
    """Selected built-in identities retain their published crosswalks."""
    ident = bystanders.GUIDE_IDENTITIES
    assert ident["AGGCAAGGCTGGCCAACCCA"]["ours"] == "60forw"
    assert ident["AGGCAAGGCTGGCCAACCCA"]["ravi"] == "gRNA 10"
    assert ident["GGCAAGGCTGGCCAACCCAT"]["ours"] == "61forw"


def test_the_catastrophic_check_fires_on_the_guide_that_demonstrated_it():
    """122forw IS Mayuranathan's sgRNA -187/189, the guide that abolished HbF induction.

    It sits in our enumerated set as an ordinary candidate. The check must catch it.
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
    guide = Guide(guide_id="synthetic_clear", protospacer="GCTAACGAACGTACGATCAA", pam="AGG",
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
