# Protocol-v1 metrics runner

The evaluation package has two deliberately separate entry points:

```sh
base-edit-safety-evaluate-tuning \
  --measurement NORMALIZED.tsv \
  --pipeline-output PIPELINE.json \
  --golden-set GOLDEN.json \
  --output-dir data/outputs/validation/tuning

base-edit-safety-evaluate-held-out \
  --measurement NORMALIZED_D5.tsv \
  --pipeline-output PIPELINE_D5.json \
  --golden-set GOLDEN_D5.json \
  --freeze-tag vX.Y.Z-validation \
  --output-dir data/outputs/validation/held-out
```

Both commands read the locked registry and protocol plus
`validation/han_overlap_map.json`. A non-fixture evaluation is rejected while the protocol contains
freeze placeholders or a registry row lacks a `pass` or `warn` status. The tuning command accepts
only datasets assigned to `tuning`. The held-out command requires an existing Git tag, accepts only
the registry's held-out set, and requires a clean checkout at that tag. Before reading values used
to compare predictions with measurements, it creates `held_out_one_shot.lock.json` in a single
filesystem operation. An existing lock causes any later attempt to fail.

Each normalized measurement row is TSV, CSV, or JSON and carries:

- `dataset_id`, `assay`, `cell_type`, `guide_id`, and `record_type`;
- `off_target_site`: `site_id` (or a coordinate), `measured_positive`, and optional
  `true_negatives_available`;
- `bystander`: `measured_value`, `editor_class`, and either `point_id` or
  `offset`/`ref`/`alt`;
- `on_target_outcome`: `measured_value`, `editor_class`, and a `point_id`.

The expected pipeline JSON uses a `records` array containing `off_target_candidate`,
`bystander_prediction`, `on_target_prediction`, and `guide_result` rows. A top-level
`dataset_metadata` object records whether candidate lists and true negatives are complete; a
top-level `tool_versions` object is copied to the file that records run inputs and versions. Current enumeration and design
record envelopes are also readable, but truncated low-mismatch candidate lists are explicitly
marked incomplete, so M1/M2 become not computable rather than misleadingly low.

The output directory contains `validation_report.md`, `metrics.json`, and `run_manifest.json`.
The report follows protocol section 9. It includes threshold status; results grouped by assay, cell
type, and editor; section 5 failure codes; limitations; source, run, and environment details; and
reviewer sign-off. Held-out M3 also separates points with no Han-source rule from points whose rule
and held-out measurement share the Han source.

For synthetic tests only, pass `--fixture-mode --fixture-root tests/fixtures`. CI is barred from
using the non-fixture path.
