# Setup and running

The core package uses the Python standard library. Small adapter modules call separately installed
prediction tools. This repository does not include those tools' source code, model weights,
reference indexes, or licenses.

## Install the package

Python 3.10 or newer is required.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[test]"
python -m pytest -q
```

Unit and synthetic-fixture tests run without external predictors. Tests marked `integration` skip
until their corresponding tools are configured.

## Configure external tools

Install each external tool from its official source and review its current license and usage
terms. `requirements-external.txt` is an inventory, not a pip requirements file.

| Tool | Role | Required for |
| --- | --- | --- |
| BE-DICT | Prediction of possible edited sequences and their frequencies | enumeration and single-guide runs |
| CRISPOR | Guide discovery and sequence-match off-target enumeration | enumeration and single-guide runs |
| AlphaGenome | Optional regulatory-consequence prediction | consequence annotation only |

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

For optional AlphaGenome annotation, supply `ALPHAGENOME_API_KEY` through the process environment.
The adapter reads the credential from the running process and does not write it to the local cache.

## Run the workflow

```bash
python -m base_edit_safety.pipeline.single_guide
python -m base_edit_safety.pipeline.enumerate_guides
python -m base_edit_safety.pipeline.annotate_consequence  # optional
python -m base_edit_safety.pipeline.rank_designs
python -m base_edit_safety.pipeline.design_records
```

Output files are written under `data/outputs/`, which is intentionally ignored. The first reference
fetch uses network access and caches reference windows under ignored `data/inputs/`.

Run the full suite after configuring BE-DICT and CRISPOR:

```bash
python -m pytest -q
```

See the [example walkthrough](example-run.md) for the output flow and
[architecture](architecture.md) for the safety boundaries.
