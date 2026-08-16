# Han 2023 source-overlap map

**Status: authored 2026-08-14. This map does not admit a dataset, change the locked split, or run the
held-out evaluation.**

## What this is

The D5 held-out measurements (Han 2023, PMID 37989316, Figure S3C) come from a paper that also
informed several bystander rules and registry claims. Results at those positions are not fully
independent of rule development. Protocol section 9 therefore reports two groups separately:

- the **source-overlap group**, stored as `partially_entangled`, contains keys listed in this map;
- the **no-shared-source group**, stored as `clean`, contains all other D5 keys.

- Source of truth: [`bystanders.py`](../src/base_edit_safety/bystanders.py) `RULES` and
  [`hbg_registry.py`](../src/base_edit_safety/hbg_registry.py) `REGISTRY`, filtered to entries whose
  `source_pmid` or `citations` include PMID 37989316.
- Machine-readable copy: [`han_overlap_map.json`](han_overlap_map.json).
- `tests/test_han_overlap_map.py` requires the map to contain every Han-citing entry.

## Key convention

The key is `(offset, ref, alt, editor_class)`. Offsets use the promoter-sense strand and are written
as negative positions relative to the transcription start site. Registry offsets are converted to
that form. For bystander rules, the allele change determines the editor class: `A>G` or `T>C` maps
to ABE, and `C>T` or `G>A` maps to CBE. Registry entries use their stored `editor_class`.

## Bystander rules citing Han 2023 (5)

Keyed and classified from `bystanders.RULES`. `functional_class` and `evidence_class` are the rule's
own values, unchanged.

| offset | ref>alt | editor_class | functional_class | evidence_class | source PMIDs |
| ---: | :--- | :--- | :--- | :--- | :--- |
| -116 | A>G | ABE | weakly_productive | measured_rate_and_binding_affinity | 37989316, 35147495 |
| -115 | C>T | CBE | productive | supported_haplotype_context | 37989316 |
| -114 | C>T | CBE | productive | supported_haplotype_context | 37989316, 35147495 |
| -112 | A>G | ABE | adverse | motif_structure_only | 37989316, 35147495, 37400614 |
| -110 | A>G | ABE | adverse | motif_structure_only_with_contrary_isolated_null | 37989316, 35147495, 34341563 |

## Registry claims citing Han 2023 (2)

Keyed from `hbg_registry.REGISTRY` entries whose `citations` include a PMID 37989316 claim. Both are
HPFH promoter entries; no clinical program entry cites Han. The Han citation is `motif_base_editing`
in primary CD34+ HSPCs, recorded as `not_accessed` / `not_reviewed`, so the overlap is at the
motif-tiling level, not an audited per-allele measurement.

| offset | ref>alt | editor_class | entry_id | motif | strength | mechanism_evidence | Han citation (role / access / audit) |
| ---: | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| -117 | G>A | CBE | hpfh_-117_CBE | BCL11A | motif | not_reviewed | motif_base_editing / not_accessed / not_reviewed |
| -114 | C>T | CBE | hpfh_-114_CBE | BCL11A | motif | in_vitro_binding_only | motif_base_editing / not_accessed / not_reviewed |

## Source-overlap group (6 distinct keys)

The group contains every distinct key in the two tables above:

| offset | ref>alt | editor_class | present in |
| ---: | :--- | :--- | :--- |
| -117 | G>A | CBE | registry (hpfh_-117_CBE) |
| -116 | A>G | ABE | bystander rule |
| -115 | C>T | CBE | bystander rule |
| -114 | C>T | CBE | **bystander rule and registry (hpfh_-114_CBE)** |
| -112 | A>G | ABE | bystander rule |
| -110 | A>G | ABE | bystander rule |

`-114 C>T CBE` appears in both a bystander rule and a registry claim that cite Han 2023. It therefore
has two forms of source overlap and must be identified in the report.

## No-shared-source group

Any D5 `(offset, ref, alt, editor_class)` key not listed above belongs to the no-shared-source group.
Its rule did not use Han 2023. Protocol section 9 reports this group separately from keys whose rules
and measurements share that source.
