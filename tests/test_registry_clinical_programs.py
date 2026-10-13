"""Check the completeness and sources of HBG-promoter program records."""

from base_edit_safety import hbg_registry


def test_rm_001_and_yolt_204_are_cited_clinical_programs():
    by_id = {entry.entry_id: entry for entry in hbg_registry.REGISTRY}

    for entry_id in ("rm_001", "yolt_204"):
        entry = by_id[entry_id]
        assert entry.strength == "clinical"
        assert entry.offsets == ()
        assert entry.program_status == "active"
        assert entry.citations
        assert all(citation["access_level"] == hbg_registry.NOT_ACCESSED
                   for citation in entry.citations)
        assert all(citation["audit_status"] == hbg_registry.NEVER_REVIEWED
                   for citation in entry.citations)

    assert by_id["rm_001"].modality == "nuclease"
    assert by_id["yolt_204"].modality == "base_editing"
    assert by_id["yolt_204"].editor_class is None
    assert {citation["registry_id"] for citation in by_id["yolt_204"].citations} >= {
        "NCT06678165", "ChiCTR2500112517",
    }


def test_unopened_rm_001_and_yolt_204_sources_cannot_create_same_motif_matches():
    """An unverified BCL11A assignment must not become operational matching data."""
    by_id = {entry.entry_id: entry for entry in hbg_registry.REGISTRY}
    guarded_ids = {"rm_001", "yolt_204"}

    for entry_id in guarded_ids:
        entry = by_id[entry_id]
        assert entry.motif == "unknown"
        assert all(citation["access_level"] == hbg_registry.NOT_ACCESSED
                   for citation in entry.citations)
        assert all(citation["audit_status"] == hbg_registry.NEVER_REVIEWED
                   for citation in entry.citations)

    match = hbg_registry.relationship("HBG1", 300, "ABE", "BCL11A motif")
    matched_ids = {entry["entry_id"] for entry in match["matched_entries"]}
    assert guarded_ids.isdisjoint(matched_ids)


def test_every_clinical_program_has_an_explicit_status():
    for entry in hbg_registry._CLINICAL:
        assert entry.program_status == "active" or entry.program_status.startswith("discontinued")


def test_every_moved_out_program_has_a_reason():
    for entry in hbg_registry.MOVED_OUT_OF_REGISTRY:
        assert entry.get("reason", "").strip()

    by_label = {entry["label"]: entry for entry in hbg_registry.MOVED_OUT_OF_REGISTRY}
    for label in ("BRL-101 (Bioraylab / 邦耀生物)", "ET-01 (EdiGene / 博雅辑因)"):
        assert by_label[label]["reason"] == \
            "BCL11A erythroid enhancer, not the HBG1/HBG2 promoter"
        assert by_label[label]["citations"]
        assert all(citation["access_level"] == hbg_registry.NOT_ACCESSED
                   for citation in by_label[label]["citations"])
        assert all(citation["audit_status"] == hbg_registry.NEVER_REVIEWED
                   for citation in by_label[label]["citations"])


def test_editor_class_requires_a_full_text_read_citation():
    for entry in hbg_registry.REGISTRY:
        if entry.editor_class is None:
            continue
        assert any(
            citation["audit_status"] == hbg_registry.FULL_TEXT_READ
            for citation in entry.citations
        ), f"{entry.entry_id} asserts editor_class={entry.editor_class!r} without full-text support"


def test_every_clinical_entry_has_a_trial_registration_or_explicit_unread_sponsor_ref():
    """A 'clinical' strength claim must be backed by at least one trial-registration reference, or
    else by an explicitly unread sponsor placeholder (audit_status NEVER_REVIEWED) that names it
    needs a dated, specific source. RM-001 is the sole entry that qualifies only via the sponsor
    branch: no registry ID was verified for it, so it carries an unopened sponsor_material ref
    instead of a trial_registration. This stops a 'clinical' label from ever standing on nothing.
    """
    clinical = [e for e in hbg_registry.REGISTRY if e.strength == "clinical"]
    assert clinical, "no clinical entries found; the filter or STRENGTH_ORDER changed"
    for entry in clinical:
        has_trial_registration = any(
            citation.get("role") == "trial_registration" for citation in entry.citations
        )
        has_unread_sponsor = any(
            citation.get("role") == "sponsor_material"
            and citation.get("audit_status") == hbg_registry.NEVER_REVIEWED
            for citation in entry.citations
        )
        assert has_trial_registration or has_unread_sponsor, (
            f"{entry.entry_id}: strength 'clinical' needs a trial_registration ref or an explicit "
            f"unread (NEVER_REVIEWED) sponsor ref; it has neither"
        )
