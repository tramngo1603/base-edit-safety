# Install and run the tool

The core package uses the Python standard library. CRISPOR, BE-DICT, and the conditional
CRISPRon-ABE fallback are configured as separate tools and keep their own installation and license
terms.

## Install the package

The package requires Python 3.10 or newer. My recorded installation checks used Python 3.14.4
on macOS ARM64; those checks do not establish support on every listed version or platform.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[test]"
python -m pytest -q
```

The unit tests run without external predictors. Tests marked `integration` skip until their
corresponding tools are configured.

## Run the included example

The example uses the checked-in coordinate table and no external data or predictor:

```bash
python -m base_edit_safety.pipeline.example_record \
  --output /tmp/base-edit-safety-example.json
```

It writes one synthetic review record for each intended HBG promoter copy. See the
[field guide](example-run.md) for the expected states and field meanings.

The separate candidate-panel example also uses only the standard library and constructed test data:

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

See the [candidate-panel guide](candidate-panel.md). Real broader searches still require a
separately installed search engine. Included adapters normalize Cas-OFFinder, CHANGE-seq-BE
identified tables, experimental BED6 lists, and declared CRISPRme+ 2.3.2 integrated-results TSVs
while recording the source files, versions and settings. CRISPRme+ setup and import are documented separately in the
[CRISPRme+ import guide](crisprme-import.md).

## Connect external prediction tools

Install each external tool from its own distribution and review its current license and usage
terms. `requirements-external.txt` lists the integrations; it is not installed by pip.

| Tool | Role | Required for |
| --- | --- | --- |
| BE-DICT | Predicts edited sequences and their relative scores | NGG intended-edit and bystander predictions |
| CRISPOR | Finds guides and matching reference-genome sites | Guide search |
| CRISPRon-BE | Runs a separate ABE model when needed | ABE8e guides whose BE-DICT models disagree on outcome order |

Copy the template and replace paths locally:

```bash
cp .env.example .env
```

The local `.env` is ignored by Git. Equivalent process environment variables take precedence.

```bash
export BE_SAFETY_BEDICT_REPO=/path/to/bedict
export BE_SAFETY_BEDICT_PYTHON=/path/to/bedict/.venv/bin/python
export BE_SAFETY_CRISPOR_REPO=/path/to/crispor
export BE_SAFETY_CRISPOR_PYTHON=/path/to/crispor/.venv/bin/python
export BE_SAFETY_CRISPOR_GENOME=hg38
export BE_SAFETY_CRISPOR_GENOMEDIR=/path/to/crispor/genomes
export BE_SAFETY_CRISPOR_PAM=NGG
export BE_SAFETY_CRISPRON_ABE_REPO=/path/to/crispron-BE
export BE_SAFETY_CRISPRON_ABE_PYTHON=/path/to/crispronbe/bin/python
```

`BE_SAFETY_CRISPOR_PAM` defaults to `NGG` and may use another PAM pattern supported by the local
CRISPOR installation, such as `NG`, `NGA`, `NGAG`, or `NGN`. For a non-NGG design, CRISPOR still
enumerates guides and off-target candidates, and the remaining record checks still run. BE-DICT is
not called because its shipped models are NGG-trained; outcome and bystander rates are recorded as
`not_evaluated` instead.

CRISPRon-BE is conditional. Install it from its own distribution only after reviewing its BSL 1.1
license and the licenses of its CRISPRon/CRISPRoff dependencies. Keep its source, model archives,
and environment outside this repository. The adapter accepts the pinned v1.0 checkout, calls
`bin/CRISPRonBE.sh` in a separate process, and uses the authors' ABE8e weight vector
`0-0-0-0-1`. When configured, it runs for an ABE8e guide whose BE-DICT models disagree on the order of
predicted outcomes. Its result uses `confidence=model_predicted_fallback`; the original BE-DICT
exclusion stays in the output. These relative model scores are not measured editing rates. Without this configuration, the same guide remains
`routed_not_executed` with null rates.

## Run the workflow

```bash
python -m base_edit_safety.pipeline.enumerate_guides
python -m base_edit_safety.pipeline.design_records
```

Results are written under the ignored `data/outputs/` directory. The first reference fetch uses
network access and caches reference windows under `data/inputs/`.

Run the full suite after configuring BE-DICT and CRISPOR:

```bash
python -m pytest -q
```

See the [example walkthrough](example-run.md) for the output flow and the
[workflow overview](architecture.md) for the main checks.
