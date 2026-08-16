# Base Edit Safety

Research triage for CRISPR base-editing designs at the duplicated **HBG1/HBG2**
fetal-globin promoters, with an initial focus on sickle cell disease and beta-thalassemia research.

Base Edit Safety combines CRISPOR guide and off-target searches, BE-DICT edit-outcome predictions,
and an optional AlphaGenome regulatory-effect adapter. It turns those inputs into a reviewable
record for each guide and editor pair.

## What the record makes reviewable

| Review dimension | Recorded output |
| --- | --- |
| Intended and bystander edits | Predicted edited sequences, position-level edit chances, and explicit `not_evaluated` states |
| Off-target scope | Sequence-matched candidate sites kept separate from editor-specific, guide-independent, and RNA concerns |
| HBG1/HBG2 duplication | Uncertain gene-copy assignment, possible targeting of both copies, and deletion concerns that require assays able to distinguish HBG1 from HBG2 |
| Clonal hematopoiesis | Overlap flags for a fixed driver-gene list, without patient-risk or incidence claims |
| Evidence strength | `predicted_only`, `analog_supported`, or `empirically_measured` labels linked to their supporting source |
| Run history | Reference coordinates, tool and code versions, fields that could not be verified, and rules that limit claims to the available evidence |

## Intended users

- **Researchers** comparing which HBG1/HBG2 guide/editor designs should continue experimental follow-up.
- **Translational sponsors** reviewing evidence gaps, assay requirements, and reproducibility before
  committing validation resources.
- **Regulatory and safety reviewers** inspecting how predictions, published analogs, measurements,
  and unresolved claims are kept distinct.

One caveat I have to add: The tool supports triage and review; it does not approve a design or determine efficacy, toxicity,
or clinical readiness.

For those new to base editing, I have the [plain-language glossary](docs/glossary.md) that defines the
genetics and validation terms used in the records and documentation.

## Install

Python 3.10 or newer is required.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[test]"
python -m pytest -q
```

BE-DICT, CRISPOR, and AlphaGenome are optional external installations governed by their own terms. See
[the setup guide](docs/setup.md) for configuration.

## Example command

After configuring CRISPOR and BE-DICT, enumerate candidates across the fixed GRCh38 promoter
windows:

```bash
python -m base_edit_safety.pipeline.enumerate_guides
```

The command above writes a JSON file under `data/outputs/` with the tool, code, reference, and source
details needed to trace the run. Generated files are ignored by Git. The
[example walkthrough](docs/example-run.md) describes the fields without publishing current rankings
or held-out measurements.

## Validation status

Validation is underway. See
[validation status and safeguards](docs/validation.md).

## Checks that protect scientific claims

- `anchors.py` contains the fixed coordinate table; tests check GRCh38 alleles, strand conventions,
  reference bases, and the required 4,924 bp HBG1/HBG2 spacing.
- Predictor, reference, and code details travel with generated files.
- The held-out runner can execute once per locked code release and refuses datasets that have not passed
  the documented admission checks.

## Repository map

```text
src/base_edit_safety/   workflow, external-tool adapters, scoring, and evaluation checks
tests/                  unit, fixed-value, source-record, and synthetic-input tests
docs/                   setup, architecture, glossary, coordinates, example, and validation notes
data/README.md          data boundary; raw inputs and generated outputs are excluded
validation/             versioned protocol and locked split metadata
```

The [data policy](docs/data-policy.md) states what the public repository intentionally excludes.

## License

The project is licensed under [Apache License 2.0](LICENSE). External tools, models, publications,
and datasets retain their own licenses and terms; this repository does not relicense or redistribute
them.
