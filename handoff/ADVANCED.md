# Advanced use and rebuilding

The existing 0.4.1 bundle uses `Start.command`. For the 0.5.0 wheel and donor coverage
workflow, follow the [updated quickstart](../docs/research-quickstart.md) in a separate environment.
The bundle installs CPython 3.12.12
(python-build-standalone build 20260203) and hash-pinned wheels from the bundle,
then runs the mock reference/sample workflow on the included native backend.
No system Python or shell environment is used. This path was tested in a new
local environment on macOS 14.4.1 arm64, Apple M3 Pro, 18 GiB RAM.

Run from the extracted bundle. These commands use the same preparation and search
functions as the launcher. Replace the example output and cache paths with new
local folders of your choice; paths with spaces are supported.

```sh
./Start.command --setup-only
.venv/bin/python -I -m base_edit_safety.launcher --bundle "$PWD" \
  --run small --output "$PWD/small-example"
.venv/bin/python -I -m base_edit_safety.launcher --bundle "$PWD" \
  --run public --cache "$PWD/public-cache" --output "$PWD/hg003-example"
```

For an already generated request, `bes-search --config request.json --output NEW_FOLDER`
remains available in `.venv/bin`. The launcher calls the same `research.run` function.
Keep replay requests with their input files; generated absolute paths name the local
files actually used, rather than a development checkout. The public cache's input
paths are relative to that cache and can be moved before another launcher run.

## Public data and expected results

`source.tar.gz` contains `tools/prepare_public_example.py` and its two helpers.
Installed copies are in `.venv/share/base-edit-safety/tools`. The sources, exact
hashes and preparation steps are in `PUBLIC_INPUTS.md`. The UI and preparation
script check resources before downloading. Preparation retains `source-manifest.json`
and `preparation-receipt.json`; searches retain `run-manifest.json`. These are the
existing release record formats.

`expected-hg003.json` records the accepted public workflow for its package version.
The original 0.4.1 bundle retains the 0.4.0 results recorded by revision
`8982108ac9622285e52a9a056b503252c1818ace`; the 0.5.0 package uses the completed donor-coverage run.
It compares all candidate rows, locus
changes, exclusions and scientific coverage. Candidate evidence links are excluded
because the private curation ledger is not distributed. JSON keys and membership
lists (configuration IDs and copy labels) are sorted. Each complete row is hashed;
sorted row hashes are then hashed together, retaining duplicate rows. This ignores
display order while preserving sequences, ordered coordinate paths, copy labels,
phase status, exclusion reasons and all other scientific fields. Coverage excludes
only `reference_native_accounting.execution_receipt_sha256`, which hashes a record
containing runtime and command paths. That receipt remains in the run's coverage
and file manifest. Runtime, timestamps, run identity, package identity and backend
command paths are recorded separately. The bundled example reproduces the
complete public sequence workflow.

Historical 0.4.0/0.4.1 results retained all 4,000,097 source rows. The 306 spanning-deletion markers triggered
whole-sample-chromosome exclusions on 15 autosomes. Modeled chromosomes were chr5,
8, 11, 12, 13, 14 and 21. Reference-only and sample counts have different genomic
scope. Reference loci on those seven autosomes number 191,771; 188,767 are retained,
2,681 gained, 1,126 lost, 3,701 unpaired indel paths, and 1,878 unresolved reference
loci. Another 514,817 reference loci are outside modeled sample scope. The combined
output has 715,143 sequence configurations and 2,164 retained loci with changed
sequences. These historical totals do not describe the new all-autosome processing.
The [donor coverage report](../docs/donor-coverage.md) explains current local exclusions,
unchanged denominators and reproduction. Missingness, phase and input exclusions remain explicit.

For other ploidy, phase policies or quality thresholds, consult the installed
`share/base-edit-safety/docs/research-advanced.md` input contract. The simple launcher
supports diploid autosomes only. It never infers X/Y copies from a name or sex label.

## Rebuild Cas-OFFinder

The included executable is permitted by the retained BSD license in
`licenses/Cas-OFFinder-LICENSE.txt`. It links only Apple's OpenCL framework,
libc++ and libSystem. Apple's framework remains an OS dependency and is not bundled.
The source archive pins commit `9816b94c20c4cba2e79b039e1e2a6dee684b7b66`.
`cas-offinder-2.4.1-compat.patch` is byte-identical to the cleanup release patch,
label `local-workgroup-id-tail-guard-aligned-upload-v3`. Its fixes cover workgroup
indexing, terminal access, aligned uploads and device-buffer access flags.

Install Apple Command Line Tools if needed using `xcode-select --install`, complete
the system dialog, then run:

```sh
./build-backend.sh "$PWD/native-rebuild"
```

The script uses included CMake 4.0.3, Release mode, Apple clang/clang++, deployment
target 14.4, and disables optional OpenMP discovery. The tested compiler is Apple
clang 15.0.0 (clang-1500.3.9.4). Exact commands and backend hash are in `RELEASE.json`.
Toolchain versions can change binary hashes; the native smoke example must pass
before using any rebuilt executable. No compiler is needed for ordinary installation.

## Source build and installed-package acceptance

All runtime, build and test wheels have versions and hashes in the lock files.
`dependencies.json` pins the source/runtime downloads. License and notice files
are retained under `licenses`, within the wheels and in the Python archive.

```sh
.venv/bin/python -I -m pip install --no-index --require-hashes \
  --find-links vendor/wheels -r build-requirements.lock
tar -xzf source.tar.gz
.venv/bin/python -I -m pip wheel --no-index --no-deps --no-build-isolation \
  --wheel-dir rebuilt-wheels ./source
```

To test the installed package, without importing from `source/src`:

```sh
BASE_EDIT_SAFETY_FIXTURES_ONLY=1 BES_OFFLINE=1 .venv/bin/python -I -m pytest \
  -q -o pythonpath= -m 'not integration' source/tests
```

Optional tests that require private curation, extra external tools or offline network
fixtures remain skipped. `RELEASE.json` records the final counts and the exact
environment actually tested. Earlier releases and their evidence are preserved.
