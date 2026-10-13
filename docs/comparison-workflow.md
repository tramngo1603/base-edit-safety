# Compare designs and plan experiments

Use `bes-compare` to ask which editing designs warrant further testing, which
genotype-associated sequence findings need investigation, and which complete
experiments fit a declared budget. The workflow keeps sequence findings,
measured observations, optional predictions and prospective plans separate.

Start with the [complete offline example](../examples/comparison/workflow/README.md).
It runs every core command without external tools and shows why a complete
18-unit contrast cannot be completed at 17 units.

## Declare the study

A study supplies canonical designs, editor geometry, intended loci, experimental
contexts and endpoints. Each active design has its own search profile over the
same declared domain. A design's exact sequence, editor and intended edit are
part of its identity; source aliases remain traceable. Keep HBG1 and HBG2 loci
separate even when their local sequences match.

Input references contain local paths, formats, roles and SHA-256 checksums.
Acquire references, public genotypes, source workbooks and extraction maps
separately. Pin their source releases and transformations before analysis.
`validate` checks the declarations, record links and supplied file bytes. It
does not independently verify a physical specimen, perform a genomic search or
establish that scientific metadata are true.

| Mode | What must be declared |
| --- | --- |
| `synthetic_test` / `synthetic_truth` | Mock inputs with known expected results |
| `benchmark_simulation` / `retained_callset_model` | Exact public callset and assumptions; missing calls stay unknown; no empirical donor credit |
| `benchmark_simulation` / `callable_supported` | Explicit source support, with hypothetical specimen assumptions labeled |
| `research_sample_planning` / `callable_supported` | Actual specimen availability, permitted use, verified assay feasibility, supported copies and strict independence declarations |

Unknown is not reference sequence or noncarriage. A resolved noncarrier requires
both actual copies. Public DNA availability does not establish that an assayable
specimen exists. The mock mode examples demonstrate these declaration
checks without claiming any real specimen exists.

See the [record contracts](comparison-contracts.md) and
[sequence input guide](comparison-sequences.md) for formats. The full genomic
bridge needs a pinned population source configuration in addition to the new
comparison catalog; a catalog-only declaration example is not a genomic input.

## Run and interpret the stages

| Command | Required stage arguments | Result |
| --- | --- | --- |
| `validate` | None | Contract and source-byte checks |
| `search` | None | Donor sequences, intended assessments, nominations and scope |
| `observations` | None | Source measurements, normalization revisions and active states |
| `evidence` | `--search`, `--observations` | Evidence applicability under current designs and contexts |
| `compare` | `--search`, `--evidence` | Current design/donor summaries and conditional decisions |
| `questions` | `--comparison` | Declared hypotheses, all required arms/controls and eligible operations |
| `plan` | `--questions`, `--evidence` | Complete prospective contrasts and resource charges |
| `report` | `--comparison`, `--questions`, `--plan` | Readable decision report and detailed exports |

Each command takes `--config`. Every stage command also takes `--output` (or
`--output-dir`). Dependencies are explicit so a new budget can refer to an
unchanged search without copying or regenerating it. Plan and report reuse
recomputes the scientific plan against the verified input records.

Search always assesses intended loci independently of nomination. A missing
nomination cannot erase an intended PAM loss or unresolved donor copy. Counts
refer to the declared genomes, domain and search profile; they are not patient
prevalence, editing probabilities or genome-wide completeness estimates.

The local backend runs offline. Native search requires profiles for the pinned
CRISPRme+ adapter and a separate `--native-output` directory. Backend selection
is inferred from active profiles; an explicit `--backend` must agree. Native
work uses an already installed pinned image with network access disabled and
image pulling disabled. See [native execution](comparison-external.md).

Declare optional `policies.search.execution_limits` before execution. The
`local` and `host` objects accept positive integer `timeout_seconds` and
`max_resident_bytes`. Defaults are 600 seconds and 2 GiB. Host elapsed time
covers source preparation, native waiting and normalization; peak memory is
for the whole Python process. Native limits include `timeout_seconds`,
`memory_bytes`, `cpus`, `disk_bytes`, `pids` and `max_rows`. They are bounded
execution limits, separate from each profile's mismatch, gap and row policies.
Recorded caps and partial statuses must be inspected before interpreting scope.

## Connect measured evidence

Use a pinned raw source file and an explicit extraction map. Preserve raw values,
units, source replicate IDs, reported labels and row/cell locators. Supply the
exact guide/editor, cell context, locus, tested sequence, donor/specimen linkage,
controls and quality information when known. Missing fields retain a reason.
Published results from another editor or cell type can motivate a question
without satisfying its empirical donor requirements.

Importing later measurements creates a new assessment. Declare earlier
canonical JSONL using the `observation_history` input role; all declared history
is included automatically. `observations --previous PATH` verifies that PATH is
a declared history input. Earlier `raw_observations` files marked `jsonl` also
mean canonical history. Keep the original source file declarations so old
revisions remain verifiable. The active state chooses the current normalization;
old source values and revisions remain available.

A completed assay needs observed passing controls and quality checks. Paying
for a future control cannot repair failed measurement quality. Fully planned
and empirically evaluable are separate state fields. Read the
[evidence guide](comparison-evidence.md) for exact, related, incompatible and
unresolved matches, including retained published positives outside search scope.

## Plan the next complete comparison

Declare biological sample counts, technical replicates, within-donor pairing,
independence groups, exact endpoint/locus pairs, compatible controls and evidence
reuse rules. Prices are integer units with an explicit scale. Distinct resource
instances have finite capacity; shared setup, donor, specimen and control costs
are charged once only when their compatibility rules permit sharing.

The planner first funds a jointly feasible mandatory-question seed. Optional
work cannot consume the budget needed by that seed. Five bounded selectors and
a tiny exact solver use the same eligibility and cost rules. All excluded
questions, missing slots, ineligible donors and enumeration limits remain
visible. A heuristic failing to find a bundle does not prove no experiment
exists. See [question requirements](comparison-questions.md) and
[planning semantics](comparison-planning.md).

Budget and algorithm choices belong in the study configuration, not unrecorded
command flags. To link a revised plan, add hashed `plan_history` inputs for the
original configuration, original plan directory's `run-manifest.json` and
`plan.json`, each with stage `planning`. Then pass `--previous-config` and
`--previous-plan` to `plan`. The command verifies that earlier study, its complete
ancestry and its plan before recording `previous_plan_id`. Keep both output
directories. Old plans are never rewritten by new observations or budgets.

Use [the benchmark command](comparison-benchmark-command.md) to compare selectors
against one frozen schedule. Freeze a real study's complete scientific protocol
before examining its donor/design outcomes. A later selector request binds
derived question/evidence files; it cannot establish an earlier freeze by itself.

## Replay, outputs and errors

The [input and file-hash guide](comparison-provenance.md) describes stage-specific cache
keys and the separate whole-run identity. A budget change invalidates planning
and reporting. A changed genotype, design or search profile invalidates search
and its descendants. Raw observations affect normalization and its descendants;
context/control policy and selected annotations affect their consuming stages.
A report/privacy change affects the report. Changed package source conservatively
invalidates comparative cached stages. V3 retains its original run-ID checks.

Completed output directories are immutable. Use `--reuse` for verified reuse
or choose a new directory. Retain partial-status records and command logs from
interrupted or capped work. The selector benchmark publishes its complete directory
after all scheduled cases finish. Add `--save-cases` with a new directory to retain
each checked plan and case progress before the whole schedule finishes. The case
archive supports inspection and audit; it does not certify or resume a full run.
Check the actual process before restarting: losing a monitor does not establish
that its computation stopped. Preserve the original run and its deadline, and
record any monitoring gap. Completed native jobs may be reused after their full
retained inputs, outputs and limits have been checked. See the
[benchmark guide](comparison-benchmark-command.md) for its recovery limits.

Commands return JSON status without donor IDs or local paths. Exit codes are
0 for a completed calculation, 2 for invalid/unavailable local input, 3 reserved
for acquisition failure, and 4 for capped or interrupted execution. Ordinary
analysis performs no acquisition. A completed plan with unmet requirements can
exit 0; inspect `meets_request`, planning status and unmet requirements.

The Markdown report begins with the decision and the evidence available for it.
[Report exports](comparison-reports.md) provide observations, applicability,
intended changes, exact assay sequences, controls, resource totals, unanswered
questions and search scope. Default prose omits direct donor identifiers and
local paths. Detailed result files stay local and are not part of the public package;
aggregation alone does not guarantee anonymity.
