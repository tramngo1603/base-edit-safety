# Setup and running

The core package uses the Python standard library. CRISPOR and BE-DICT are configured as separate
tools and keep their own installation and license terms.

## Install the package

Python 3.10 or newer is required.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[test]"
python -m pytest -q
```

The unit tests run without external predictors. Tests marked `integration` skip until their
corresponding tools are configured.

## Configure predictor adapters

Install each external tool from its own distribution and review its current license and usage
terms. `requirements-external.txt` lists the integrations; it is not installed by pip.

| Tool | Role | Required for |
| --- | --- | --- |
| BE-DICT | Base-editing allele-spectrum prediction | enumeration |
| CRISPOR | Guide discovery and sequence-match off-target enumeration | enumeration |

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
```

## Run the workflow

```bash
python -m base_edit_safety.pipeline.enumerate_guides
python -m base_edit_safety.pipeline.design_records
```

Artifacts are written under the ignored `data/outputs/` directory. The first reference fetch uses
network access and caches reference windows under `data/inputs/`.

Run the full suite after configuring BE-DICT and CRISPOR:

```bash
python -m pytest -q
```

See the [example walkthrough](example-run.md) for the output flow and the
[workflow overview](architecture.md) for the main checks.
