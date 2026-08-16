# Data boundary

The public repository intentionally contains no patient-level records, protected health
information, raw sequencing reads, proprietary pharmaceutical datasets, downloaded article PDFs,
or held-out validation outcomes.

`data/inputs/` and `data/outputs/` are local runtime directories and are ignored by Git:

- `inputs/` may hold reference caches or source data obtained by the user under the source's own
  access and reuse terms.
- `outputs/` holds generated candidate records, run records, rankings, and validation reports.

Synthetic and deliberately minimal examples used by the test suite live under `tests/fixtures/` and
are labelled as synthetic there. Source accessions and checksums belong in the validation registry;
large or licensed source files do not belong in this repository.

Before sharing any output file, review it for consent or reuse restrictions, sample or donor
identifiers, absolute local paths, credentials, and unpublished results.
