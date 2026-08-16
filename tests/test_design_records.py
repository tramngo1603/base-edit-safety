"""The design records: evidence tiers, shape, off-target kinds, paralog flags, and the checklist."""
from base_edit_safety import annotation, evidence, hbg_registry
from base_edit_safety.design_record import build_design_records

_FORBIDDEN = ("clonal", "expansion", "leukemia", "leukaemia", "mds", "lifetime", "risk",
              "trajectory", "incidence")


def _rec(guide_id, gene, offset, role, rate, *, ref="T", alt="C", sref="A", salt="G",
         pos=5_249_970, editor="ABE", offtarget=None, ch=None, co_target=None):
    prov = {"guide_id": guide_id, "editor": "ABE8e", "protospacer": "CTTGACCAATAGCCTTGACA",
            "pam": "AGG", "strand": "-", "predictor_runs": 5}
    if offtarget is not None:
        prov["offtarget"] = offtarget
    if ch is not None:
        prov["ch_driver_offtargets"] = ch
    if co_target is not None:
        prov["co_targets_paralog"] = True
        prov["co_target_gene"] = co_target
    return {"genome_build": "GRCh38", "chrom": "chr11", "pos_1based": pos,
            "genomic_plus_ref": ref, "genomic_plus_alt": alt,
            "promoter_sense_ref": sref, "promoter_sense_alt": salt,
            "gene": gene, "gene_strand": "-", "promoter_offset": offset,
            "editor_class": editor, "on_target_or_bystander": role, "predicted_rate": rate,
            "model_name": "BE-DICT-bystander:ABE8e", "consequence": "deferred", "provenance": prov}


def _records():
    offt = {"source": "CRISPOR", "offtarget_count": "10", "mit_specificity": "40",
            "mm_histogram": {"2": 1, "4": 5}, "cfd_weighted_burden": 0.3, "paralog_count": 1}
    ch = [{"ch_driver_gene": "TET2", "chrom": "chr4", "start": "105275000",
           "mismatch_count": "1", "feature": "exon"}]
    return [
        _rec("HBG1_m113_ABE", "HBG1", 113, "on_target", 0.45, offtarget=offt, ch=ch, co_target="HBG2"),
        _rec("HBG1_m113_ABE", "HBG1", 116, "bystander", 0.10, pos=5_249_973, offtarget=offt,
             ch=ch, co_target="HBG2"),
    ]


def test_evidence_tier_enforcement():
    rec = build_design_records(_records())["records"][0]
    rm = rec["registry_match"]
    assert rm["relationship"] == hbg_registry.EXACT_EDIT       # -113 ABE is a known HPFH edit
    assert rm["evidence_strength"] == "published_hspc"         # kept, not blurred into a weaker match
    assert rm["measured_in"] == "literature"                   # published, not measured here
    assert rec["evidence_tier"] == evidence.ANALOG_SUPPORTED   # a registry match, with no measurement
    rs2 = build_design_records(_records(),
                               empirical={"HBG1_m113_ABE": [{"assay": "amplicon", "edit_rate": 0.4}]})
    assert rs2["records"][0]["evidence_tier"] == evidence.EMPIRICALLY_MEASURED
    assert rs2["records"][0]["registry_match"]["measured_in"] == "this_run"


def test_no_measurement_means_no_measured_tier():
    rec = build_design_records(_records(), empirical={})["records"][0]
    assert rec["evidence_tier"] != evidence.EMPIRICALLY_MEASURED


def test_predicted_only_when_no_registry_and_no_empirical():
    rec = build_design_records([_rec("g", "HBG1", 999, "on_target", 0.4)])["records"][0]
    assert rec["evidence_tier"] == evidence.PREDICTED_ONLY


def test_design_record_schema():
    rec = build_design_records(_records())["records"][0]
    required = {"record_type", "schema_version", "candidate", "evidence_tier", "registry_match",
                "predicted_edit_spectrum", "off_target", "hbg_paralog_flags", "ch_annotation",
                "empirical_support", "required_validation", "scoring", "provenance",
                "limitations"}
    assert required <= set(rec)
    assert [e["promoter_offset"] for e in rec["candidate"]["target"]["intended_edits"]] == [113]


def test_offtarget_blind_spots_are_record_level_limitations():
    limitations = build_design_records(_records())["records"][0]["limitations"]
    required = (
        "Off-target candidates come from a Cas9-style search of the reference genome. It finds "
        "mismatch-only sites: no bulged sites and no population-variant-created sites. Measured "
        "assays (GUIDE-seq / CHANGE-seq-BE) detect bulged sites; population-variant off-targets "
        "need a variant-aware module that does not exist in v1."
    )
    assert any(limitation.startswith(required) for limitation in limitations)


def test_offtarget_kinds_stay_separate():
    kinds = build_design_records(_records())["records"][0]["off_target"]["kinds"]
    assert set(kinds) == {"spacer_homology_cas_style", "base_editor_context",
                          "guide_independent_deaminase", "rna_off_target"}
    assert kinds["spacer_homology_cas_style"]["mm_histogram"] == {"2": 1, "4": 5}
    assert kinds["base_editor_context"]["status"] == "not_assessed_here"
    assert kinds["rna_off_target"]["status"] == "not_assessed_here"


def test_paralog_flags_set():
    pf = build_design_records(_records())["records"][0]["hbg_paralog_flags"]
    assert pf["flagged"] is True
    assert any("also hits the paralog" in f for f in pf["flags"])
    assert pf["required_confirmation"]


def test_required_validation_checklist():
    validates = {i["validates"] for i in
                 build_design_records(_records())["records"][0]["required_validation"]}
    assert {"on_target_edit_rate", "bystander_spectrum", "base_editor_context", "rna_off_target",
            "hbg_paralog_rearrangement", "ch_driver_offtarget"} <= validates


def test_ch_consequence_reads_gene_and_feature_together():
    ch = build_design_records(_records())["records"][0]["ch_annotation"]["driver_offtargets"]
    cc = ch[0]["consequence_class"]
    assert cc["driver_mechanism"] == "lof_tumor_suppressor"             # TET2
    assert cc["consequence_relevance"] == "coding_disruption_relevant"  # and it is in an exon


def test_no_field_anywhere_claims_a_risk():
    def keys(o):
        if isinstance(o, dict):
            for k, v in o.items():
                yield k
                yield from keys(v)
        elif isinstance(o, list):
            for x in o:
                yield from keys(x)
    for key in keys(build_design_records(_records())):
        assert not any(bad in key.lower() for bad in _FORBIDDEN), key


def test_reproducibility_block_present():
    repro = build_design_records(_records())["reproducibility"]
    assert "generated_at" in repro
    assert "commit" in repro["repo"] and "dirty" in repro["repo"]
    assert "rubric_weights" in repro["config"]


def test_score_sits_under_scoring_not_at_the_top():
    rec = build_design_records(_records())["records"][0]
    assert "safety_score" not in rec            # the score is nested under its explanation
    scoring = rec["scoring"]
    assert "tiebreaker" in scoring["note"].lower() and "safety_score" in scoring


def test_a_design_with_two_edits_keeps_both():
    # The -123/-124 pair: both edits stay, and both get compared to the registry.
    recs = [_rec("HBG_pair", "HBG1", 123, "on_target", 0.30, pos=5_249_980),
            _rec("HBG_pair", "HBG1", 124, "on_target", 0.28, pos=5_249_981)]
    rec = build_design_records(recs)["records"][0]
    assert {e["promoter_offset"] for e in rec["candidate"]["target"]["intended_edits"]} == {123, 124}
    assert len(rec["registry_match"]["by_edit"]) == 2
    assert {a["relationship"] for a in rec["registry_match"]["by_edit"]} == {hbg_registry.EXACT_EDIT}
    assert len(rec["predicted_edit_spectrum"]["on_target"]) == 2


def test_ch_mechanisms_are_gene_specific():
    assert annotation.ch_consequence_class("TP53", "exon")["driver_mechanism"] == \
        "tumor_suppressor_mixed_lof_dominant_negative"
    assert annotation.ch_consequence_class("DNMT3A", "exon")["driver_mechanism"] == \
        "loss_or_hypomorphic_with_R882_hotspot"
    mdm4 = annotation.ch_consequence_class("MDM4", "exon")
    assert mdm4["driver_mechanism"] == "copy_number_or_overexpression_driver"
    assert mdm4["consequence_relevance"] == "snv_consequence_poorly_represented"


def _citation_text() -> str:
    return " ".join(str(v) for e in hbg_registry.REGISTRY for c in e.citations for v in c.values())


def test_registry_carries_no_unchecked_dois():
    assert "10.1056" not in _citation_text()   # an unchecked DOI has no business in an evidence record


def test_every_citation_states_its_role_and_whether_it_was_read():
    for entry in hbg_registry.REGISTRY:
        for c in entry.citations:
            assert {"role", "audit_status"} <= set(c), c
            # exactly one kind of identifier: a paper or a registry/company reference, never both
            assert ("pmid" in c) ^ ("registry_id" in c), c
            assert c["audit_status"] in (hbg_registry.FULL_TEXT_READ,
                                         hbg_registry.NEVER_REVIEWED,
                                         hbg_registry.NOT_REVIEWED), c


def test_coleman_is_still_unread_and_doerfler_is_not():
    """Preserve both citations and the review status recorded for each one.

    The full-text review of Doerfler 2021 supports two positions on each dependency branch. Coleman
    1993 and the six primary sources behind its phase claim remain unreviewed.
    """
    status = {c["pmid"]: c["audit_status"]
              for e in hbg_registry.REGISTRY for c in e.citations if "pmid" in c}
    assert status["7679879"] == hbg_registry.NEVER_REVIEWED
    assert status["34341563"] == hbg_registry.FULL_TEXT_READ


def test_synthetic_positions_do_not_claim_to_be_hpfh_alleles():
    ids = {e.entry_id for e in hbg_registry.REGISTRY}
    assert "synthetic_-123_ABE" in ids and "synthetic_-124_ABE" in ids
    assert not any(i.startswith("hpfh_-123") or i.startswith("hpfh_-124") for i in ids)


def test_mechanism_evidence_is_recorded_and_never_assumed():
    by_id = {e.entry_id: e for e in hbg_registry.REGISTRY if e.offsets}
    # Each of the three reviewed papers provides binding tests outside cells only.
    for entry_id in ("hpfh_-198_ABE", "synthetic_-123_ABE", "synthetic_-124_ABE", "hpfh_-114_CBE"):
        assert by_id[entry_id].mechanism_evidence == hbg_registry.IN_VITRO_BINDING_ONLY
    # Entries without a reviewed paper must retain the `NOT_REVIEWED` status.
    for entry_id in ("hpfh_-175_ABE", "hpfh_-117_CBE", "hpfh_-113_ABE"):
        assert by_id[entry_id].mechanism_evidence == hbg_registry.NOT_REVIEWED


def test_the_octamer_paper_carries_no_tal1_claim():
    """Keep the octamer and TAL1 claims attached to their respective sources."""
    for c in hbg_registry._CITATIONS[(175, "T", "C")]:
        if c["pmid"] == "2449926":
            assert "tal1" not in c["role"].lower()
            assert "octamer" in c.get("note", "").lower()


def test_promoter_program_list_has_expected_entries_and_excludes_enhancer_programs():
    """Require the recorded HBG1/HBG2 promoter programs and exclude enhancer programs.

    exa-cel is deliberately excluded from the count: it targets the BCL11A enhancer, not the HBG
    promoter, so it is a comparator rather than a member of this list.
    """
    by_id = {e.entry_id: e for e in hbg_registry.REGISTRY}
    for entry_id in ("beam_101", "reni_cel_edit301", "cs_101", "otq923"):
        assert entry_id in by_id, f"{entry_id} missing from the promoter-program list"
        assert "HBG1" in by_id[entry_id].genes and "HBG2" in by_id[entry_id].genes


def test_discontinued_programs_are_not_presented_as_running():
    by_id = {e.entry_id: e for e in hbg_registry.REGISTRY}
    assert by_id["reni_cel_edit301"].program_status == "discontinued_2024_12"
    assert "does not state a safety or efficacy result" in by_id["reni_cel_edit301"].note
    # Program state and evidence category are separate fields.
    assert by_id["reni_cel_edit301"].strength == "clinical"


def test_beam_101_still_does_not_assert_an_editor_class():
    """Leave editor class unset while the reported source remains unreviewed."""
    by_id = {e.entry_id: e for e in hbg_registry.REGISTRY}
    assert by_id["beam_101"].editor_class is None
    assert "risto-cel" in by_id["beam_101"].label


def test_unreviewed_programs_keep_every_citation_marked_unreviewed():
    by_id = {e.entry_id: e for e in hbg_registry.REGISTRY}
    for entry_id in ("cs_101", "otq923"):
        for c in by_id[entry_id].citations:
            assert c["audit_status"] == hbg_registry.NEVER_REVIEWED, c
