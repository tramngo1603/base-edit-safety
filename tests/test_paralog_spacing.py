"""HBG1 and HBG2 positions must stay 4,924 bp apart. If they do not, a row has shifted.

The table below is copied independently from `docs/coordinates.md` so it can detect a typo in
`anchors.py`.
"""
from base_edit_safety import anchors

# (promoter_offset, hbg1_pos, hbg2_pos), copied from docs/coordinates.md
FROM_THE_DOC = [
    (198, 5_250_055, 5_254_979),
    (175, 5_250_032, 5_254_956),
    (124, 5_249_981, 5_254_905),
    (123, 5_249_980, 5_254_904),
    (122, 5_249_979, 5_254_903),
    (117, 5_249_974, 5_254_898),
    (114, 5_249_971, 5_254_895),
    (113, 5_249_970, 5_254_894),
]


def test_tss_separation():
    assert anchors.HBG2_TSS - anchors.HBG1_TSS == anchors.PARALOG_SPACING_BP


def test_every_position_in_the_code_matches_the_doc():
    in_code = {p.promoter_offset: (p.hbg1_pos, p.hbg2_pos) for p in anchors.IN_SCOPE}
    assert set(in_code) == {off for off, _, _ in FROM_THE_DOC}
    for off, hbg1, hbg2 in FROM_THE_DOC:
        assert in_code[off] == (hbg1, hbg2), f"offset {off}: code has {in_code[off]}"


def test_paired_positions_differ_by_the_spacing():
    for p in anchors.IN_SCOPE:
        assert p.hbg2_pos - p.hbg1_pos == anchors.PARALOG_SPACING_BP, p.promoter_offset


def test_positions_are_the_tss_plus_the_offset():
    for p in anchors.IN_SCOPE:
        assert p.hbg1_pos == anchors.HBG1_TSS + p.promoter_offset
        assert p.hbg2_pos == anchors.HBG2_TSS + p.promoter_offset
