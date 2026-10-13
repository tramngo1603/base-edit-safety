# Importing the published HBG sources

`base_edit_safety.population.source_imports` creates extraction maps for two exact primary source files:

| Source key | Primary source file | Extraction |
| --- | --- | --- |
| `liang_endov_hbg` | Liang et al., DOI `10.1038/s41467-018-07988-z`, source-data workbook | `Figure 2d`, physical cells C6:H25: 120 values, including ABE7.10 and GFP controls |
| `mayuranathan_figure4` | Mayuranathan et al., DOI `10.1038/s41588-023-01434-7`, figure source-data workbook | `Figure 4a–c`: 81 editing, 27 F-cell and 27 HbF values |

The wrapper pins each workbook's SHA-256, exact sheet/XML member, cell set and count. An unrelated workbook with the same number of values is rejected. Guide/editor bindings must match the actual source combinations. A control cannot acquire a treatment design, and the two optimized editor-dose contexts cannot be collapsed.

Provide a reviewed cell fixture and its canonical source/design/context/endpoint bindings. Both are local runtime inputs. The wrapper builds strict maps and then calls the production `observations.import_artifact` reader, which reads the actual workbook bytes and verifies every expected value and metadata cell.

```python
from base_edit_safety.population.source_imports import import_published_source

records = import_published_source(
    "liang_endov_hbg",
    workbook_path,
    extraction_map_path,
    canonical_catalog,
    reviewed_source_cells,
    canonical_source_bindings,
)
```

`build_source_map(source_name, catalogs, curated_cells, canonical_bindings)` can prepare the strict map separately. `import_published_source` preserves an existing identical map and rejects overwriting a different map revision. Store an updated mapping under a new path and merge observation history through `observations.merge_history`.

The source fixture identifies each physical cell with `source_record_key`, the file checksum, sheet, XML member, A1 address, original numeric text, percent unit, source arm and replicate labels. Functional fixtures also identify source-local participants, panels, timepoints and the corresponding donor-label cells. The canonical bindings resolve the reviewed source file, design, context and endpoint. Metadata assertions are tied to the measured cell's actual row and column; a value cannot borrow a valid donor or replicate label from a neighboring row.

Every source value remains as its original decimal text. Percent-to-fraction conversion is explicit. When a JSON binary float would change the exact converted value, normalization retains decimal text. The treatment/control header remains the source label. Liang's separate treated-site Yes/No validation outcome is retained on treated observations; it is not assigned to GFP controls.

Each emitted normalization also carries the schema v2 outcome pair `adjudicated_outcome` and `adjudication_reason`. The mapping never derives an outcome from a measured number, a treatment name or a missing row. Treated rows are `inconclusive`: their verbatim source label (Liang's Yes/No validation call; May's numeric cells report no categorical label) has no reviewed four-state adjudication yet, so conclusive credit is withheld pending scientific review. GFP and UT control arms carry no treatment design and are `untested`; that is not a negative experiment. Neither workbook reports a numeric detection limit, so no source cell can be adjudicated `tested_no_signal`. If a reviewer later supplies a reviewed adjudication, it lands as a new mapping revision under a new path, not an edit to a preserved map. The current mapping revision is `<source>-reviewed-outcomes-v2`.

The individual values do not establish technical versus biological repeat classifications beyond their source labels. May Figure 4 retains three source donor units across its repeated experiments and endpoints. Those local participants are not joined to public genomic-panel donors. Liang's HEK293T replicates do not become independent human donors. Missing genotype, numeric assay limits, per-cell QC and specimen identity remain explicit unknowns. These source imports receive no empirical public-panel donor credits.

HbF fraction, F-cell fraction and editing fraction remain separate endpoints. The functional HbF result measures a combined response rather than an individual HBG copy. Mayuranathan's ABE7.10-NG and ABE8e-NG regimens used different optimized doses, so the imported observations do not establish an equal-dose intrinsic-editor effect.

The public tests contain no source measurements or source-participant identifiers. To verify the complete pipeline against the approved curation files, set `BES_COMPARISON_CURATION_RUNTIME` to the runtime directory and run:

```sh
PYTHONPATH=src python -m pytest tests/test_comparison_source_imports.py
```

Without that variable, the immutable-source boundary tests still run; the real-file checks are explicitly skipped. Keep original workbooks, source-cell fixtures, extraction maps and normalized individual records in approved storage outside Git.
