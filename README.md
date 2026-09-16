# Base Edit Safety

I built Base Edit Safety to help researchers check the evidence behind a base-editing guide:
**which DNA sites need attention, what was tested, and what still needs an experiment.**

My first use case is the **HBG1/HBG2 fetal-globin promoters**, which are studied as targets for
beta-thalassemia and sickle cell disease. This tool brings candidate sites, predictions and
experimental results into one report. It keeps them separate, so a predicted site cannot look
like a measured edit, and a missing result cannot look like a clean result.

## System design

![Base Edit Safety: check candidate identity, match results to the tested DNA, and review what remains unknown](docs/system-chart.png)

This tool checks the guide, editor, genome coordinates and DNA sequence before matching an
experimental result to a candidate. HBG1 and HBG2 stay separate even where their sequences match.
A negative result needs a stated detection limit, controls and passed quality checks. A test of
reference DNA does not establish what happens to a different DNA version carried by a person.
These checks use the supplied experiment details; they do not independently verify the experiment.

The chart's bulges are gaps in reported guide–DNA alignments. This version accepts supported
DNA-only or RNA-only gaps; it does not accept mixed DNA/RNA gaps or establish their physical
pairing. Window annotations for gapped alignments cover bases paired to guide positions;
unpaired DNA bases can be omitted. A "no base" label does not rule out an editable base at a bulge.
See the [candidate guide](docs/candidate-panel.md) and [import limits](docs/crisprme-import.md).

## Five-minute review

Install the package and run the included example:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[test]"
python -m pytest -q
python -m base_edit_safety.pipeline.example_record \
  --output /tmp/base-edit-safety-example.json
```

The example writes two design records, one for each HBG copy. It runs no prediction model and
uses no paper workbook or patient data. Prediction values stay empty because no model was run.
The [field guide](docs/example-run.md) explains the output.

## Candidate-panel workflow

Use this part of the tool when you already have candidate lists from searches or experiments.
It checks the input files, combines sites by identity, and records which DNA versions have
matching experimental evidence. Your choices about which sites to test stay separate from that
evidence.

Try the constructed example:

```bash
bes-candidate-panel \
  --profile examples/candidate-panel/search-profile.json \
  --guides examples/candidate-panel/guides.tsv \
  --intended-loci examples/candidate-panel/intended-loci.tsv \
  --nominations examples/candidate-panel/nominations.tsv \
  --source-provenance examples/candidate-panel/source-provenance.json \
  --measurements examples/candidate-panel/measurements.tsv \
  --decisions examples/candidate-panel/panel-decisions.tsv \
  --output /tmp/base-edit-safety-candidate-ledger.json \
  --panel-output /tmp/base-edit-safety-panel-review.tsv
```

This writes an eight-site JSON record and a TSV table. All example data are
constructed for software checks. Five example results lack the identity of the tested DNA, so
none can establish editing for a particular DNA version. Three other sites have no results.

The [candidate guide](docs/candidate-panel.md) explains the input files. It covers Cas-OFFinder,
CRISPRme+, CHANGE-seq-BE and experimental BED6 lists. Searches run separately; this import step
does not run those tools. To make a Markdown report, use `bes-candidate-report` with your
matching design records as described in that guide.

## What I use it to check

| Question | What the tool records |
| --- | --- |
| Which sites were proposed for testing? | Each source list, its search settings and any rejected or skipped rows |
| Which DNA version does a result apply to? | Guide, editor, genome location, strand, sequence and supplied evidence for the tested allele |
| What might be edited near the intended base? | Intended and bystander predictions, within the model's stated limits |
| Could both HBG copies be affected? | Separate HBG1/HBG2 locations and the need for assays that distinguish them and check larger deletions or rearrangements |
| Does a candidate overlap a blood-cell clonal expansion gene? | Gene overlap for review, without estimating anyone's clinical risk |
| What should be checked next? | Missing experiment details, untested DNA versions and suggested follow-up assays |

Both HBG promoters are intended targets in this use case. A site counts as intended only if its
coordinates match the supplied target list. An extra perfect match elsewhere stays in the
candidate list. This tool does not combine these checks into a single safety score.

## What I have checked so far

- In the CHANGE-seq-BE/Digenome-seq primary panel of 4,291 rows, 95 off-target rows passed the
  imported significance and quality flags. CHANGE-seq-BE had nominated 88/95; the combined
  methods had nominated 95/95 because that is how the panel was selected. A four-mismatch NGG
  sequence rule retained 20/95. Applying CFD and MIT average-rank cutoffs of 100 per guide to
  that measured panel retained 35/95 and 27/95. These are retrospective panel comparisons,
  not independent estimates of search sensitivity.
- Five HBG paper reruns accounted for the rows selected for each check. They tested whether
  this tool could recover source identities and records, not whether it could reproduce the
  experiments.
- Comparing the Mayuranathan CIRCLE-seq list with its targeted sequencing panel found 66
  four-mismatch, NGN candidates without that follow-up; 62 contain an adenine in the configured
  editing window. These were already in the authors' candidate list. They are **untested
  candidates, not new off-target findings**.
- A small real-data CRISPRme+ 2.3.2 run produced the two intended HBG reference sites. Its
  search and merge finished, but plotting failed; a separate unchanged export step produced
  the imported results. It did not demonstrate a real variant-created off-target.

The [validation summary](docs/validation.md) gives the denominators, source citations and limits.
I share the results of these checks, not copies of the papers' individual data rows or detailed
extracted tables. See the [data policy](docs/data-policy.md).

## Run guide searches and predictions

This separate path needs CRISPOR and, for supported NGG predictions, BE-DICT installed first.
Follow [setup](docs/setup.md), then run:

```bash
python -m base_edit_safety.pipeline.enumerate_guides
python -m base_edit_safety.pipeline.design_records
```

The commands write JSON under the local, Git-ignored `data/outputs/` directory. Non-NGG designs
can still be recorded, but their BE-DICT rates stay empty with `confidence=not_evaluated`.

## Limits to keep in mind

- The ordinary CRISPOR path uses the reference genome. It keeps sites with up to three
  mismatches, plus four-mismatch sites with CFD at least 0.10. It does not search population
  variants or gapped sites. Broader searches can be imported under their recorded settings.
- Variant-aware imports retain reported DNA versions. They do not prove that every possible
  version was found or determine a donor's genotype. Unknown sequence, alignment and phase
  remain unresolved. Imported frequency labels are not recalculated as population frequencies.
- Prediction scores are not calibrated probabilities of editing in cells. Missing predictions
  stay missing; they are never replaced with zero.
- Guide-independent DNA editing, RNA editing, delivery, dose and long-term outcomes need
  separate evidence. No output from this tool establishes biological or clinical safety.

## Documentation

- [Setup and external tools](docs/setup.md)
- [Example and output fields](docs/example-run.md)
- [Candidate files and evidence checks](docs/candidate-panel.md)
- [CRISPRme+ import](docs/crisprme-import.md)
- [Validation results](docs/validation.md)
- [How the parts fit together](docs/architecture.md)
- [GRCh38 coordinate definitions](docs/coordinates.md)
- [Terms used in the records](docs/glossary.md)
- [Public data policy](docs/data-policy.md)

## License

My project code is under the [Apache License 2.0](LICENSE). External tools, models, publications
and datasets keep their own licenses and terms. This repository does not redistribute paper
workbooks, individual experimental data points or tables extracted from them.
