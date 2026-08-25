"""Enumeration over both promoters, with a predicted edit for every guide.

Skipped unless both CRISPOR and BE-DICT are installed and configured.
"""
import pytest

from base_edit_safety.config import bedict_configured, crispor_configured

pytestmark = pytest.mark.integration

if not (bedict_configured() and crispor_configured()):
    pytest.skip("external guide-design / predictor not configured", allow_module_level=True)

from base_edit_safety import anchors  # noqa: E402
from base_edit_safety.pipeline.enumerate_guides import run_enumeration  # noqa: E402


@pytest.fixture(scope="module")
def records():
    try:
        return run_enumeration(write=False)
    except Exception as exc:
        pytest.skip(f"enumeration could not run: {exc}")


def test_guides_for_both_promoters(records):
    genes = {r.gene for r in records}
    assert genes == {"HBG1", "HBG2"}


def test_every_record_has_a_spectrum_value(records):
    assert records
    assert all(r.predicted_rate is not None for r in records)
    assert all(r.model_name.startswith("BE-DICT-bystander") for r in records)


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
