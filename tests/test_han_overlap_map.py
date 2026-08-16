"""Check coverage of the Han 2023 source-overlap group used by protocol section 9.

`validation/han_overlap_map.json` lists every bystander rule and registry claim that cites Han 2023
(PMID 37989316), which is also the source of the D5 held-out data. These tests rebuild the set from
the source modules and require the committed map to match it. See `validation/HAN_OVERLAP_MAP.md`
and validation-protocol.md section 9.
"""
import json
from pathlib import Path

from base_edit_safety import bystanders, hbg_registry

ROOT = Path(__file__).resolve().parents[1]
MAP_JSON = ROOT / "validation" / "han_overlap_map.json"
MAP_MD = ROOT / "validation" / "HAN_OVERLAP_MAP.md"
HAN_PMID = "37989316"

# Editor class follows the allele change: an ABE gives A>G and a CBE gives C>T on the sense strand.
_EDITOR_BY_ALLELE = {("A", "G"): "ABE", ("T", "C"): "ABE", ("C", "T"): "CBE", ("G", "A"): "CBE"}


def _doc() -> dict:
    return json.loads(MAP_JSON.read_text())


def _source_bystander_rules() -> list:
    return [r for r in bystanders.RULES if HAN_PMID in r.source_pmid]


def _source_registry_entries() -> list:
    return [e for e in hbg_registry.REGISTRY
            if any(c.get("pmid") == HAN_PMID for c in e.citations)]


def _han_citation_keys() -> set:
    return {key for key, claims in hbg_registry._CITATIONS.items()
            if any(c.get("pmid") == HAN_PMID for c in claims)}


def test_map_files_exist():
    assert MAP_JSON.is_file()
    assert MAP_MD.is_file()


def test_json_covers_every_han_citing_bystander_rule():
    mapped = {(b["offset"], b["ref"], b["alt"]): b for b in _doc()["bystander_rules"]}
    source = _source_bystander_rules()
    # Exactly the Han-citing rules, no more, no fewer.
    assert {(r.offset, r.ref, r.alt) for r in source} == set(mapped)
    for r in source:
        b = mapped[(r.offset, r.ref, r.alt)]
        assert b["functional_class"] == r.functional_class, r.key
        assert b["evidence_class"] == r.evidence_class, r.key
        assert b["editor_class"] == _EDITOR_BY_ALLELE[(r.ref, r.alt)], r.key
        assert HAN_PMID in b["source_pmid"], r.key


def test_json_covers_every_han_citing_registry_entry():
    mapped = {e["entry_id"]: e for e in _doc()["registry_entries"]}
    source = _source_registry_entries()
    assert {e.entry_id for e in source} == set(mapped)
    for e in source:
        m = mapped[e.entry_id]
        assert m["offset"] == -abs(e.offsets[0]), e.entry_id
        assert m["editor_class"] == e.editor_class, e.entry_id
        assert m["strength"] == e.strength, e.entry_id
        assert m["motif"] == e.motif, e.entry_id
        assert m["mechanism_evidence"] == e.mechanism_evidence, e.entry_id
        assert m["han_citation"]["pmid"] == HAN_PMID, e.entry_id


def test_partially_entangled_keys_are_the_union_and_double_entanglement_is_flagged():
    doc = _doc()
    keys = {tuple(k) for k in doc["partially_entangled_stratum_keys"]}

    bystander_keys = {(r.offset, r.ref, r.alt, _EDITOR_BY_ALLELE[(r.ref, r.alt)])
                      for r in _source_bystander_rules()}
    han_cit_keys = _han_citation_keys()
    registry_keys = set()
    for e in _source_registry_entries():
        off = e.offsets[0]
        ref, alt = next((rf, al) for (o, rf, al) in han_cit_keys if o == off)
        registry_keys.add((-abs(off), ref, alt, e.editor_class))

    assert keys == bystander_keys | registry_keys
    # A key in both sources is doubly entangled and must be called out as such.
    assert {tuple(k) for k in doc["double_entangled_keys"]} == bystander_keys & registry_keys


def test_counts_and_header_are_consistent():
    doc = _doc()
    assert doc["pmid"] == HAN_PMID
    assert doc["counts"]["bystander_rules"] == len(doc["bystander_rules"]) == len(_source_bystander_rules())
    assert doc["counts"]["registry_entries"] == len(doc["registry_entries"]) == len(_source_registry_entries())
    assert doc["counts"]["distinct_keys"] == len(doc["partially_entangled_stratum_keys"])
    assert doc["counts"]["double_entangled_keys"] == len(doc["double_entangled_keys"])
    assert doc["key_fields"] == ["offset", "ref", "alt", "editor_class"]


def test_markdown_lists_every_key_and_flags_the_clean_stratum():
    text = MAP_MD.read_text()
    for offset, ref, alt, editor in _doc()["partially_entangled_stratum_keys"]:
        assert str(offset) in text, offset
        assert f"{ref}>{alt}" in text
        assert editor in text
    assert "hpfh_-114_CBE" in text          # the doubly entangled registry entry
    assert "No-shared-source group" in text
