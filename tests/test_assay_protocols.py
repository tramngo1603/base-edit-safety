"""Published primer sets attached to the checklist items that had only been naming an assay.

Ravi 2022 Supplementary file 2 supplies protocols for BOTH paralog checklist items and for all
eleven off-target loci. A protocol documents how an assay can be run; it does not claim a result.
"""
from __future__ import annotations

import re

from base_edit_safety import assay_protocols as ap
from base_edit_safety import offtarget_reference as ref
from base_edit_safety import paralog_flags

DNA = re.compile(r"^[acgtACGT]+$")


def _flagged_record():
    return {"gene": "HBG1", "provenance": {"co_targets_paralog": True, "co_target_gene": "HBG2"}}


def test_both_paralog_checklist_items_now_carry_a_protocol():
    """Attach executable primer protocols to both HBG gene-copy checklist items."""
    items = {i["validates"]: i for i in paralog_flags.paralog_flags(_flagged_record())["required_confirmation"]}
    for key in ("hbg_paralog_cotargeting", "hbg_paralog_rearrangement"):
        assert items[key].get("published_protocols"), f"{key} still has no runnable protocol"


def test_the_paralog_separating_primers_are_the_gene_specific_ones():
    protocols = ap.protocols_for("hbg_paralog_cotargeting")
    names = {x["name"] for p in protocols for x in p["primers"]}
    assert {"HBG1 F", "HBG2 F", "HBG1/2 R"} <= names
    pcr = next(p for p in protocols if p.get("amplicon_size") == "940 bp")
    # the forwards differ per gene, the reverse is shared: that is what separates the paralogs
    fwd = [x["sequence"] for x in pcr["primers"] if x["name"].endswith("F")]
    assert len(fwd) == 2 and fwd[0] != fwd[1]


def test_the_deletion_qpcr_keeps_its_real_attribution():
    """Ravi published it, but it is Chang Li's assay. Absorbing that would be a citation error."""
    p = ap.protocols_for("hbg_paralog_rearrangement")[0]
    assert "Chang Li" in p["attribution"]
    assert "35147495" in p["source"]


def test_adapter_overhangs_are_separated_from_what_anneals():
    """Separate the adapter overhang from the target-binding sequence used in an order."""
    ngs = [p for p in ap.protocols_for("hbg_paralog_cotargeting") if "NGS" in str(p["primers"])]
    assert ngs
    for x in ngs[0]["primers"]:
        assert "adapter_overhang" in x
        assert not x["sequence"].startswith(x["adapter_overhang"])


def test_every_primer_is_dna():
    for p in ap.ALL:
        for primer in p.primers:
            assert DNA.match(primer.target_specific), (p.method, primer.name)
            if primer.adapter:
                assert DNA.match(primer.adapter), primer.name


def test_all_eleven_off_target_loci_have_amplicons():
    assert len(ap.OFFTARGET_AMPLICONS) == 11
    for site in ref.RAVI_2022_GRNA11.sites:
        assert ap.offtarget_primers_for(site.name) is not None, site.name


def test_the_primer_names_independently_confirm_the_guide_assignment():
    """The panel table had no guide column; the assignment was derived by alignment.

    The primer names read "G11 OT1".."G11 OT11", so the study names the guide itself. A derived
    claim and a second stated identifier should remain consistent.
    """
    ot1 = ap.offtarget_primers_for("OT1")
    assert all(x["name"].startswith("G11 ") for x in ot1["primers"])
    assert ref.RAVI_2022_GRNA11.guide_names["ravi"] == "gRNA 11"


def test_missed_sites_come_with_the_primers_to_check_them():
    """Attach a published amplicon protocol to each missed candidate site."""
    out = ref.compare_to_panel([], ref.RAVI_2022_GRNA11)
    assert len(out["missed_by_our_search"]) == 11
    amplicons = out["amplicons_for_missed_sites"]
    assert set(amplicons) == set(out["missed_by_our_search"])
    assert all(v is not None for v in amplicons.values())


def test_a_protocol_is_never_mistaken_for_a_result():
    """Do not represent a published protocol as a result from this workflow."""
    for p in ap.ALL:
        blob = (p.method + p.note + p.attribution).lower()
        for claimed in ("measured", "confirmed", "validated", "detected"):
            assert claimed not in blob, (p.method, claimed)
