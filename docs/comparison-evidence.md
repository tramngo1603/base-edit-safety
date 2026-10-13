# Importing and using measured evidence

The comparison workflow reads pinned source files, verifies their extraction maps,
and creates immutable source measurements and versioned normalizations. Source
measurements remain available even when a search admits no matching candidate.

Run the mock import example after installing the package:

```sh
bes-compare observations --config examples/comparison/observations/study.json --output /tmp/bes-observation-example
bes-compare observations --config examples/comparison/observations/study.json --output /tmp/bes-observation-example --reuse
```

Use a new output directory for changed inputs. `--reuse` verifies the source bytes,
code identity, stored records and output hashes.

Each extraction map names one source file ID, its file format, extraction and
mapping versions, expected measurement count and exact rows. Each row identifies
the source cell or record, expected raw value, raw unit, experiment, replicate,
source label and declared normalization metadata. Unknown identity, assay quality,
controls and detection limits require reasons. The importer never supplies these
facts from a measured fraction or a treatment name.

The generic JSON format is an array of objects with exactly these string fields:
`measurement_key`, `source_study`, `source_experiment`, `source_replicate`,
`raw_value`, `raw_unit`, and `source_reported_label`. The TSV format uses those
same columns. Numeric source values are strings so their reported decimal text
survives. Generic row locators use `{"kind":"file","member":"rows/1"}` for the
first data row. Every generic row needs a mapping.

Workbook maps use the physical sheet, worksheet XML member, A1 cell, row and
column. Missing physical row numbers do not move subsequent measurements.
Optional metadata assertions verify treatment, donor or replicate label cells.
The original XML row is retained. A repeated source cell is rejected. The legacy
candidate-ledger bridge requires schema 3, preserves each measurement object and
checks its source study, assay, replicate, unit and reported outcome against the
map. It rejects a claimed replicate that the ledger does not report. A resolved
measured sequence must agree with the observation's reported locus.

Percent-to-fraction conversion uses an explicit factor of 0.01. A normalized
value uses a JSON number when its decimal value survives serialization; otherwise
it retains finite decimal text. For example, `9.8000000000000007` percent becomes
`"0.098000000000000007"` fraction. The raw value remains unchanged. Censored or
unavailable values retain their reason and do not become zero.

A new extraction or mapping version creates a new normalization of the same
source measurement. The active-state record selects one revision and retains the
complete history. Declare canonical JSONL with the hashed `observation_history`
input role. Every declared history is included automatically; `observations
--previous PATH` checks that the named file is among those inputs. Earlier
`raw_observations` inputs with format `jsonl` have the same history meaning.
Changing the flag alone cannot change the computed assessment or its cache key. A documented
observation-equivalence record can link duplicate exports of one experiment.
Equal numbers alone do not establish equivalence.

An explicitly completed design-neutral mock can satisfy its matched control
requirement without an editing-design ID. Import checks its executed mock
operation, specimen, donor, context, source experiment, endpoint, locus and
actual target sequence. It still requires passed controls and quality checks.
An unknown design label alone cannot establish a mock, and a mock does not
become a treated-design editing measurement during evidence reconciliation.

After search and normalization stages exist, reconcile evidence with:

```sh
bes-compare evidence --config STUDY.json --search SEARCH_DIRECTORY --observations OBSERVATION_DIRECTORY --output EVIDENCE_DIRECTORY
```

Applicability has separate identity, context, endpoint and quality assessments.
The default evidence view evaluates active designs in each source's actual
reported context. An intended experiment needs its own explicit context and
endpoint assessment. The library's `assess_observation()` also accepts an exact
target donor and sequence for that purpose.

Different editor definitions remain different designs. Different cell types,
delivery methods, doses and timepoints remain different contexts. Related data
may motivate a question without satisfying its donor-assay requirement. A pooled
HbF measurement cannot resolve a donor, an HBG copy or an editing-read fraction.
Measured zero needs a usable reported detection limit before it can provide
reassuring evidence. Public retained-callset simulations always receive zero
empirical donor-assay credits.

Source summaries count measurements, reported independent units and source donor
labels separately. Technical repeats do not become additional donors. Distinct
values within one unit remain visible; no pooled effect estimate is computed.
Predictions are rejected by the measured-observation entry point.

The `published_alignments` input role pins source-reported alignment mappings.
Replay compares exact columns, PAMs, sequences and genomic paths with admitted
alignments. Every mapped source result retains its positive/negative source label
and receives a disposition: exact recovery, a related alternative alignment,
outside policy, unsupported representation or unresolved mapping. A nearby
ungapped hit does not recover a reported bulged alignment. Replay does not assign
empirical donor credit.

A source row may leave `design_id` null with an `unresolved_reason` when the
published result does not identify a unique complete design. Its known alignment,
literal outcome and cell locator remain in the replay as an unresolved mapping.
Shared guide sequence alone cannot assign that result to a particular editor.

Local annotation reuses the package's versioned Ensembl CH-gene panel on its
supported GRCh38 assembly. Additional `annotations` inputs contain a schema version,
assembly, source, release and closed gene/regulatory feature intervals. Intervals
use zero-based, half-open coordinates. Output records bind the input checksum and
remain coordinate annotations, without estimating editing activity or outcomes.

## Explicit outcome adjudication and legacy history

New imports write observation revision schema `2`. Every revision has an
`adjudicated_outcome` and a nonempty `adjudication_reason`. The source's
`source_reported_outcome`, original decimal text and raw record stay unchanged.

| Adjudication | Conclusive evidence rule |
| --- | --- |
| `observed_positive` | Requires a usable positive measurement and the existing identity, context, endpoint, controls and QC checks. A value at or below a usable stated detection limit is flagged. |
| `tested_no_signal` | Requires a usable positive detection limit in compatible units and a numeric value from zero through that limit. |
| `inconclusive` | Withholds conclusive empirical credit, even with a numeric value and passed QC. |
| `untested` | No completed test is recorded; conclusive evidence remains unavailable. |

An exact four-state source label that conflicts with the adjudication is flagged.
Other source prose is preserved without guessing its meaning. A contradiction
never changes the measured value. No-signal values below the limit are not used
as exact numbers for paired measured-difference claims. This implementation also
withholds credit for a censored or missing numeric measurement, even when its
adjudication is no-signal and a limit is present; it retains the outcome, source
text and missing/censoring reason. Supporting such censored evidence requires a
separate reviewed extension.

An extraction map with neither new field receives an explicit inconclusive
revision and a legacy-adjudication reason. Supplying just one field is an error.
Supplying both uses the declared adjudication after validation. Importers must
not infer either positive or no-signal status from a number alone.

Version-1 revisions remain readable as history. Their old stored credit bits do
not grant current intrinsic, reconciliation, comparison or planner credit.
To migrate an existing complete JSON catalog, first retain its original file:

```sh
python -m base_edit_safety.population.comparison_migration \
  --catalog catalog-v1.json --adjudications reviewed-outcomes.json \
  --output migrated-outcomes
```

The optional adjudications file is an object keyed by exact active v1 revision
IDs. Each value contains only `adjudicated_outcome` and `adjudication_reason`.
Omit `--adjudications` to withhold conclusive credit for every unresolved legacy
revision. Unknown or already-migrated revision IDs are errors.

The new protected directory contains `catalog.json`, `observations.jsonl`,
`superseded-states.json` and `migration-receipt.json`. Original source records,
old revisions, source IDs and source bytes remain unchanged. New v2 IDs bind the
complete interpretation. Replacement states select the new revisions and retain
old IDs in their history. Superseded states are archived separately because one
catalog cannot declare two current states for the same source.

The receipt maps old active revision IDs to new IDs and lists dependent records
removed from the new catalog. Dependent assessments, summaries and plans must be
recomputed; replacing IDs inside an old conclusion is not migration. No text in
source records is searched and replaced. If invalidating a dependent record would
break a retained historical reference, migration rejects that catalog rather
than dropping history. Keep the original bundle for that case and separate its
historical products from the supported current observation inputs before review.

Pin the generated `observations.jsonl` as a checksummed `observation_history`
input with format `jsonl` in a new study. Supply the supporting records from the
new catalog, removing its observation/history records from the inline catalog to
avoid declaring the same source state twice. Update or remove old extraction maps
for migrated sources: rerunning an unadjudicated map intentionally creates a new
inconclusive normalization and must not supersede a reviewed migration silently.
Run a fresh observation-normalization stage and regenerate evidence, comparison,
questions, planning and reporting. Old stage directories stay immutable; changed
code, inputs and parent bindings reject their reuse. A migrated JSONL file is not
an old stage manifest. Retain prior plans as history, not as current evidence.
