# External donor-window search

The `crisprme_plus` adapter runs CRISPRme+ v2.5.5 on exact, reconstructed donor-copy windows. The source VCF, donor membership, ordered ALT indexes, genotype filters, phase blocks and missingness rules remain authoritative. Native population summaries and sites-only frequencies do not establish donor carriage.

Install the backend explicitly before an offline run. The adapter uses the Linux ARM64 image below and verifies the native executable hashes before each job. It never pulls an image or downloads an index during analysis.

```sh
docker pull --platform linux/arm64 \
  pinellolab/crisprme@sha256:3d5fc92e7c96cddef37e2d550e96966ac3d378966d93f2a1e9ea602258db15a6
```

This image contains release source `c005098e05ebbf83f049d9de5fd2598ceed65685`. Its [release license](https://github.com/pinellolab/crisprme-plus/blob/c005098e05ebbf83f049d9de5fd2598ceed65685/LICENSE) describes academic use under AGPL-3.0 and separate commercial terms. Keep the image, scoring models and third-party data outside this package; they are not bundled. The CRISPRitz build comes from [v2.8.2 source](https://github.com/pinellolab/CRISPRitz/tree/bb33e75126b22eaec06c1bd5cd3c0b4c3cd91bd5); its wrapper still prints `2.8.1`, so the adapter pins binary hashes as well as version strings.

An external search profile must name `crisprme_plus`, version `v2.5.5`, and the exact image digest. Supported input geometry is a 20-base spacer with a two- or three-base 3-prime PAM. Mismatch, DNA-gap, RNA-gap and total-edit limits are explicit. Total edits mean mismatches plus bulged **bases**. Event-based counting, unknown versions and other geometry fail compatibility checks. Every design is searched with its own PAM patterns, then each retained row is checked against that design's exact policy.

```sh
python -m base_edit_safety.pipeline.comparison search \
  --config /absolute/path/to/native-study.json \
  --backend crisprme_plus \
  --native-output /absolute/path/to/native-audit \
  --output /absolute/path/to/comparison-search
```

The comparison study must pin the genomic inputs and source transformation used by `prepare_inputs`. A missing or unverified source-equivalence audit produces candidate-only, partial scope. A profile for the local engine cannot be relabeled as external by a CLI override.

For explicit execution limits or cancellation, use the Python interface:

```python
from base_edit_safety.population.comparison_inputs import prepare_inputs
from base_edit_safety.population.comparison_provenance import resolve_study
from base_edit_safety.population.external_search import NativeLimits, search_prepared

prepared = prepare_inputs(resolve_study("/absolute/path/to/native-study.json"))
result = search_prepared(
    prepared,
    output_dir="/absolute/path/to/native-audit",
    limits=NativeLimits(cpus=2, memory_bytes=4 * 1024**3,
                        disk_bytes=1024**3, timeout_seconds=600),
    cancelled=lambda: False,
)
```

Each job uses argument arrays, a unique disposable container, `--pull=never`, `--network=none`, CPU/RAM/process limits and a timeout/disk watchdog. The container mounts only its generated input directory. On macOS this is an isolated `/private/tmp` directory, avoiding Docker Desktop privacy restrictions on Desktop paths. Cancellation and timeout stop the named container, retain its available output and prevent successful cache reuse. Container memory/CPU measurements, adapter peak memory, wall time, disk use, exit state and exact arguments are recorded separately. Disk checks are periodic, so a brief overshoot is possible and is reported as failure or partial scope.

Every donor window remains in `source-windows.jsonl`, including its reference identity, locus, donor, sample, copy, genotype rows, phase, sequence and per-base path. The file's hash is part of the native job identity. Writing one window at a time avoids building a second large JSON string in memory.

The adapter submits each distinct complete DNA string as one native contig. `sequence-projection.json` records all original windows represented by that string. Each retained alignment is then associated with every member's own path using the unchanged string indexes. Equal HBG1/HBG2 DNA therefore retains separate locus identities, and equal donor DNA retains each occurrence. Raw native rows, distinct submitted strings and projected alignments have separate counts. The projected row ceiling can produce a partial result even when native retrieval completed.

Missing/filtered calls, callability gaps, incompatible phase blocks and high-density skips remain in the source audit. High-density exclusions occur before strings are grouped, so an excluded window cannot receive another donor's search coverage. `retained_callset_model` remains modeled coverage even when the searched material matches its retained source exactly. String grouping changes neither empirical evidence nor genotype support.

Two explicit transformations address behavior measured in the pinned native release:

1. `excluded-native-window-guards-v1` adds generated flanking DNA solely as native index context. This avoids native coordinate/strand errors for very short contigs. A hit is admitted only if every target and PAM base lies inside the original source window. Coordinates subtract the recorded prefix length; generated bases have no donor path and never enter evaluated genomic scope.
2. `exact-leading-dna-extension-v1` derives leading DNA-gap alignments from retained shorter native anchors. The native engine omits some of these alternatives. Removing leading guide gaps preserves the full guide and PAM and lowers the DNA/total budget; the adapter reverses that operation using only exact adjacent bases in the original window. It checks strand, bounds and all budgets. Derived rows carry `alignment_origin="derived"`, the complete transformation, added path indexes and explicit native-anchor references. They are not represented as verbatim native rows.

The unchanged exhaustive test oracle enumerates full domains on both strands, including adjacent/terminal gaps and alternative alignments. It imports no production reconstruction, alignment, admission or identity helper. Actual native conformance includes both PAM lengths, genomic insertions/deletions, distinct loci, two-base gaps, mixed gaps, tight total budgets and window boundaries. Raw disagreements remain available alongside the transformed comparisons. Check each new native version against the same cases.

The adapter passes `--ci-cd-test` because this release otherwise deletes `crispritz_targets` at normal completion. It retains all raw results before parsing the ten-column CRISPRitz target format. The pipeline's headerless form and the engine's exact header are supported. Aligned sequence, zero-based forward coordinates, strand, gap types and reported counts are checked directly against the submitted DNA. Cluster positions are retained as source details, not substituted for alignment coordinates. A genomic insertion is part of the input path; an alignment bulge is a gap between guide and target columns.

Final integrated tables contain selected representatives and optional columns; they are not the complete alignment source. Their CFD/CRISTA values and native reports remain unchanged in the raw directory. CFD and CRISTA remain external nuclease scores, without base-editor calibration.

Every completed job includes exact inputs, raw outputs, executable/index hashes, parsed native anchors, transformed alignments, execution status and a file manifest. `COMPLETE.json` seals the result checksum. Reuse requires the same complete input/source transformation, design/editor/profile, effective arguments, environment, limits and adapter code. Missing or altered bytes, partial attempts, interrupted execution and changed source/code/settings cannot be reused as complete jobs. An empty result requires the native release's explicit successful zero-hit signal and retained empty raw output.

The enclosing comparison stage also binds each native directory's stage key, result checksum and completion-seal checksum. Offline `verify_native_job` checks that external binding and every retained file. Editing both a native result and its seal therefore invalidates comparison-stage and descendant reuse. Incomplete native evidence has no completion seal and remains inspectable as incomplete evidence.

Run the local regression checks without Docker:

```sh
python -m pytest -q tests/test_comparison_external.py tests/test_comparison_external_scaling.py
```

Mock native rows and their image and source details are under `tests/fixtures/comparison/external`. Full native and real-source conformance reports are kept with the validation records, where raw source genotypes and actual coverage limits can be audited.
