"""Enumeration over both promoters, with predicted rates or explicit fallback for every guide.

Skipped unless both CRISPOR and BE-DICT are installed and configured.
"""
import pytest

from base_edit_safety.config import (
    bedict_configured,
    crispron_abe_configured,
    crispor_configured,
    load_crispor_config,
)
from base_edit_safety.pipeline.datatypes import NON_NGG_MODEL_SCOPE_REASON

pytestmark = pytest.mark.integration

if not crispor_configured():
    pytest.skip("external guide-design tool not configured", allow_module_level=True)

PAM_CLASS = load_crispor_config().pam
if PAM_CLASS == "NGG" and not bedict_configured():
    pytest.skip("NGG rate predictor not configured", allow_module_level=True)

from base_edit_safety import anchors  # noqa: E402
from base_edit_safety.pipeline.enumerate_guides import run_enumeration  # noqa: E402


@pytest.fixture(scope="module")
def records():
    return run_enumeration(write=False)


def test_guides_for_both_promoters(records):
    genes = {r.gene for r in records}
    assert genes == {"HBG1", "HBG2"}


def test_every_record_has_a_rate_or_explicit_fallback(records):
    assert records
    if PAM_CLASS == "NGG":
        fallback = [
            record for record in records
            if record.provenance.get("analysis_role") == "sensitivity_fallback"
        ]
        predicted = [record for record in records if record not in fallback]
        assert predicted and all(record.predicted_rate is not None for record in predicted)
        assert all(record.model_name.startswith("BE-DICT-bystander") for record in predicted)
        assert all(record.provenance.get("fallback_predictor") ==
                   f"CRISPRon-{record.editor_class}"
                   for record in fallback)
        if fallback and crispron_abe_configured():
            assert all(record.predicted_rate is not None for record in fallback
                       if record.confidence == "model_predicted_fallback")
            assert all(record.provenance.get("fallback_status") == "executed"
                       for record in fallback)
            assert all(record.provenance.get("predictor_role") ==
                       "fallback_for_bedict_instability" for record in fallback)
        else:
            assert all(record.predicted_rate is None for record in fallback)
            assert all(record.confidence == "not_evaluated" for record in fallback)
    else:
        assert all(r.predicted_rate is None for r in records)
        assert all(r.confidence == "not_evaluated" for r in records)
        assert all(r.provenance.get("reason") == NON_NGG_MODEL_SCOPE_REASON for r in records)


def test_anchor_on_targets_found_in_both_genes(records):
    on = [r for r in records if r.on_target_or_bystander == "on_target"]
    offsets_hbg1 = {r.promoter_offset for r in on if r.gene == "HBG1"}
    offsets_hbg2 = {r.promoter_offset for r in on if r.gene == "HBG2"}
    # every targeted position should turn up in both paralogs
    for off in {p.promoter_offset for p in anchors.IN_SCOPE}:
        assert off in offsets_hbg1
        assert off in offsets_hbg2


def test_cotargeting_and_offtargets_represented(records):
    assert any(r.provenance.get("co_targets_paralog") for r in records)
    assert any(r.provenance.get("offtarget", {}).get("offtarget_count") for r in records)


def test_consequence_is_explicitly_not_evaluated(records):
    assert all(r.consequence == "not_evaluated" for r in records)
    assert all(r.provenance.get("consequence_source") == "not_evaluated" for r in records)


def test_ch_flag_populated_and_promoter_edits_not_drivers(records):
    # Every record gets a flag, and no promoter edit is inside a CH-driver gene.
    assert all(isinstance(r.ch_driver_flag, bool) for r in records)
    assert all(r.ch_driver_flag is False for r in records)  # HBG1/HBG2 are not CH-driver genes


def test_ch_driver_offtarget_surfaced(records):
    # The -123/-124 guide has a DNMT3A off-target, which has to show up.
    genes = {
        h["ch_driver_gene"]
        for r in records
        for h in r.provenance.get("ch_driver_offtargets", [])
    }
    assert "DNMT3A" in genes
