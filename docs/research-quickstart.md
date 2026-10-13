# Donor coverage quickstart — 0.5.0

Use the pinned Python 3.12.12 runtime and patched Cas-OFFinder from the accepted 0.4.1 handoff. Install the new wheel in a separate environment:

```sh
/path/to/handoff/.runtime/python/bin/python3 -m venv .venv
.venv/bin/python -m pip install '/path/to/base_edit_safety-0.5.0-py3-none-any.whl[lab]'
```

Reuse the existing checksum-pinned genomic files on approved storage. From the source checkout:

```sh
.venv/bin/python tools/prepare_public_example.py \
  --workdir /protected-storage/coverage-preparation \
  --backend /path/to/patched/cas-offinder --device G0 \
  --guides examples/research/guides.tsv \
  --reuse-prepared /protected-storage/public-preparation/out \
  --max-output-bytes 8589934592
.venv/bin/bes-search \
  --config /protected-storage/coverage-preparation/request.json \
  --output /protected-storage/coverage-results
```

The preparation command verifies all four prepared reference/VCF files and retains the complete 23-guide panel. Choose new preparation and output directories. The 8 GiB output budget requires that much free space; larger panels or outputs may need a larger explicit budget. Use the installed copies under `.venv/share/base-edit-safety/` if working without the source checkout.

Open `report.html` for the summary and `coverage.tsv` for chromosome totals. `coverage-regions.tsv` and `coverage-windows.tsv` list evaluated and excluded reference bases and window starts. Processing all autosomes does not resolve every region. Exact deletion breakpoints, missing calls and phase gaps remain explicit. X/Y requires a compatible callset plus source-supported PAR/non-PAR ploidy declarations; those inputs are absent from this HG003 example.

The existing `bes-start` launcher remains available. Keep the 0.4.1 bundle as its historical release; use the new environment's commands for this release. See [coverage definitions](donor-coverage.md) and the [advanced input contract](research-advanced.md).
