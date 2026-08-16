"""Run one guide end to end against a configured BE-DICT installation.

Skipped unless BE-DICT is configured (BE_SAFETY_BEDICT_REPO and BE_SAFETY_BEDICT_PYTHON) and the
reference sequence can be reached.
"""
import pytest

from base_edit_safety.config import bedict_configured

pytestmark = pytest.mark.integration

if not bedict_configured():
    pytest.skip("external base-editing predictor not configured", allow_module_level=True)

from base_edit_safety.pipeline.single_guide import run_single_guide  # noqa: E402


@pytest.fixture(scope="module")
def records():
    try:
        return run_single_guide(gene="HBG1", promoter_offset=113,
                                editor_class="ABE", editor_name="ABE8e", write=False)
    except Exception as exc:  # no network, or the predictor is not reachable
        pytest.skip(f"could not run: {exc}")


def test_spectrum_non_empty(records):
    assert len(records) >= 1


def test_on_target_present_and_correct(records):
    on = [r for r in records if r.on_target_or_bystander == "on_target"]
    assert len(on) == 1
    r = on[0]
    assert r.gene == "HBG1" and r.promoter_offset == 113
    assert r.pos_1based == 5_249_970
    assert (r.genomic_plus_ref, r.genomic_plus_alt) == ("T", "C")
    assert (r.promoter_sense_ref, r.promoter_sense_alt) == ("A", "G")


def test_bystander_burden_represented(records):
    bystanders = [r for r in records if r.on_target_or_bystander == "bystander"]
    assert len(bystanders) >= 1  # the -113 guide is predicted to make bystander edits


def test_provenance_and_consequence(records):
    for r in records:
        assert r.model_name.startswith("BE-DICT-bystander")
        assert "guide_id" in r.provenance
        assert r.consequence == "deferred"  # annotate_consequence.py fills this in later
        # A rate is present only where the position was sent to the predictor. Elsewhere it is None and
        # the record says so, rather than carrying a zero that looks like a prediction.
        if r.confidence == "not_evaluated":
            assert r.predicted_rate is None
        else:
            assert r.predicted_rate is not None


def test_gata_site_positions_all_reported(records):
    """The -113 A>G edit creates a GATA site spanning -113 to -110.

    -111 is a T, so -112 and -110 are the only ABE-editable positions inside it. Two papers measure
    -110 as edited in primary CD34+ cells, once at twice the rate of -112. Both must be reported.

    -110 sits at protospacer position 11, outside BE-DICT's supported 3-10 conversion window. It
    must remain present as an explicit not-evaluated record.
    """
    offsets = {r.promoter_offset for r in records}
    for expected in (116, 113, 112, 110):
        assert expected in offsets, (
            f"no record at -{expected}; emitted offsets were {sorted(offsets)}. "
            f"-110 is protospacer position 11 for this guide, one past the window edge at 10."
        )


def test_gata_site_rates_respect_the_supported_predictor_window(records):
    """An in-window rate is predicted; position 11 remains visibly unsupported, never zero."""
    by_offset = {r.promoter_offset: r for r in records}
    assert by_offset[112].confidence == "model_predicted"
    assert by_offset[112].predicted_rate is not None
    assert by_offset[110].confidence == "not_evaluated"
    assert by_offset[110].predicted_rate is None
