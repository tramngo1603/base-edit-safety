# Base Edit Safety

Base Edit Safety (BES) searches for potential off-target sites and helps researchers compare candidate DNA sequences with experimental evidence. For research use.

Version **0.5.0** includes a local researcher app, reference and sample-aware searches,
coverage reports, and candidate exports. The app searches 20-base guides with NGG PAMs,
both strands, up to six mismatches and no gaps. Missing genotypes, uncertain phasing,
unsupported ploidy and excluded regions remain visible in the results.

## Install

From this source checkout, use Python 3.12.12:

```sh
python3.12 -m venv .venv
.venv/bin/python -m pip install '.[lab]'
```

The core command-line tools support Python 3.10 or later. The `lab` extra pins pysam
0.24.1 and psutil 7.1.0. The researcher app also needs Tcl/Tk 9 and the reviewed
Cas-OFFinder 2.4.1 backend with the included [compatibility patch](tools/cas-offinder-2.4.1-compat.patch).
A wheel alone does not supply that backend.

## Launch the app

Set `BES_BUNDLE` to your local reviewed handoff folder, which contains
`backend/cas-offinder`, then run with the handoff's Python 3.12.12 runtime:

```sh
.venv/bin/bes-start --bundle "$BES_BUNDLE"
```

A complete offline handoff can instead be installed and opened with **Start.command**.
See the [app quickstart](handoff/QUICKSTART.html) for setup and the
[build instructions](handoff/ADVANCED.md) for dependencies and backend checks.
The tested native platform is macOS 14.4.1, Apple silicon, Apple M3 Pro with 18 GiB RAM.
Full-reference jobs require 16 GiB installed RAM, 4 GiB available RAM and 25 GiB free disk.

In **My analysis**, paste guides, select a reference and a results folder, then choose
**Review & run**. **Add sample variants** is optional. Select the VCF sample and confirm
its autosomal ploidy. **Try an example** includes reference-only, reference-plus-sample
and public HG003 runs. **Open results…** reopens a completed analysis.

## Small example

This mock command-line example runs without a native backend or network access.
Use a new output folder:

```sh
.venv/bin/bes-compare computational \
  --config examples/comparison/computational/request.json \
  --output ./example-results
```

Read `example-results/report.md` and `candidates.tsv`. The
[example guide](docs/computational-audit.md) explains its bounded search and expected results.
For the app's native example, choose **Try an example → Reference + sample**.

## Inputs and outputs

| Workflow | Inputs | Outputs |
| --- | --- | --- |
| Researcher app / `bes-search` | Guide text or TSV; GRCh38 FASTA; optional sample VCF and declared ploidy | `report.html`, `candidates.tsv`, `locus-changes.tsv`, exclusion and coverage tables, request and file hashes |
| Candidate review | Candidate lists, source details, measurements and panel decisions | Candidate JSON records and review tables |
| Design comparison | Study settings, donor sequences, experimental evidence and assay costs | Comparison reports and proposed assays |

In **Results**, inspect candidates, sample differences and coverage, then use **Export…**.
Export all candidates, every row matching the applied filters, or selected rows. All and
matching exports include rows beyond the 5,000-row display limit. Exports retain source
hashes, selection criteria and field descriptions in a new folder. Adding a sample starts
a new analysis and preserves the original results. Genomic inputs and results stay local.

## Documentation

- [App quickstart](handoff/QUICKSTART.html) and [advanced commands](handoff/ADVANCED.md)
- [Public example preparation](docs/research-quickstart.md), [input contract](docs/research-advanced.md) and [coverage definitions](docs/donor-coverage.md)
- [Candidate review](docs/candidate-panel.md), [donor-aware analysis](docs/population.md) and [design comparison](docs/comparison-workflow.md)
- [Architecture](docs/architecture.md), [coordinates](docs/coordinates.md) and [output fields](docs/example-run.md)
- [Validation](docs/validation.md), [historical benchmark](docs/research-benchmark.md) and [data policy](docs/data-policy.md)

## License

Project code uses the [Apache License 2.0](LICENSE). External tools, models and data
retain their own terms; see the [component and permission inventory](docs/permissions.md).
Genomic data, publisher workbooks and private review records are excluded from this repository.
