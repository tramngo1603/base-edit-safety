# Local data

The `data/` directory is the local workspace for inputs and generated results. Its contents are
ignored by Git except for this file.

- `data/inputs/` may hold reference caches or source files obtained under their original access
  and reuse terms.
- `data/outputs/` holds candidate records, run manifests, and validation artifacts.

The repository includes only small, hand-written examples under `tests/fixtures/`. Before sharing
a generated file, check its source terms and remove local paths, credentials, donor identifiers,
and unpublished results.
