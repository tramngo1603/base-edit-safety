# Try the example and read its output

I included this example so you can inspect the records without installing prediction tools or
downloading a genome or paper workbook. It uses the repository's HBG coordinate definitions.
It runs no model, includes no experimental measurement and leaves unrun checks clearly marked.

## Run the no-data example

After installing the package:

```bash
python -m base_edit_safety.pipeline.example_record \
  --output /tmp/base-edit-safety-example.json
python -m json.tool /tmp/base-edit-safety-example.json >/dev/null
```

The first command prints:

```text
wrote 2 synthetic design records to /tmp/base-edit-safety-example.json
```

The two records describe the same intended promoter change at HBG1 and HBG2. Both copies are
intended targets. The records keep their 4,924-bp separation and the need for assays that tell
the copies apart. `predicted_rate` is null and `confidence` is `not_evaluated` because no model ran.

## Main fields

| Field | Meaning |
| --- | --- |
| `candidate` | Guide, editor, location and intended base changes |
| `predicted_edit_spectrum` | Intended and bystander predictions; joint predictions keep the method and limits for changes on the same molecule |
| `off_target.kinds` | Separate sequence-match, editor-specific, guide-independent DNA and RNA concerns |
| `off_target.candidate_sites` | Search matches to review, not experimental results |
| `hbg_paralog_flags` | The two HBG copies, their spacing, possible deletions or rearrangements, and required assays |
| `ch_annotation` | Overlap with the recorded blood-cell clonal expansion genes, without a clinical-risk estimate |
| `registry_match` | Related published evidence; this is separate from a measurement of this design |
| `evidence_tier` | `predicted_only`, `analog_supported` or `empirically_measured`; the last requires attached measurements |
| `empirical_support` | Measurements attached to this design, kept separate from predictions |
| `required_validation` | Experiments needed to check intended editing, bystanders, other sites and the HBG copies |
| `provenance` | Sources, tools, model limits and reasons for missing values |
| `limitations` and `reproducibility` | What the record covers, the reference build, software revision and run time |

`promoter_offset` follows the literature's signed convention, such as `-113`.
`upstream_distance_bp` stores the positive distance used in calculations.

## Keep these distinctions when reading results

A nomination proposes a site. A prediction estimates an outcome. A measurement reports what an
experiment found under stated conditions. An experimentally nominated site without follow-up is
still an unmeasured candidate. A high CFD score is not a measurement, and an empty prediction is
not zero editing.

The candidate report checks measurements separately for each DNA version. A design's evidence
tier does not establish that all its candidate sites or alleles were tested.

## Run searches and predictions

After configuring the external tools in [setup](setup.md):

```bash
python -m base_edit_safety.pipeline.enumerate_guides
python -m base_edit_safety.pipeline.design_records
```

The first command searches the fixed HBG windows with CRISPOR and obtains supported predictions.
It writes `data/outputs/enumeration_HBG1_HBG2.json`. The second writes
`data/outputs/design_records_HBG1_HBG2.json`.

Generated files stay in the Git-ignored local workspace. They may contain local paths or source
data that are not approved for sharing. The [data policy](data-policy.md) excludes paper rows
and individual experimental values from public outputs.
