# Compare experiment selectors

`bes-compare benchmark` compares declared selectors on one local question and
evidence catalog. It uses the same eligibility, controls, costs, resource
capacities and planning limits for every selector. It does not acquire data or
run genomic searches.

Freeze the scientific study and selector schedule before evaluating outcomes.
After the question and evidence stages finish, create a selector request that
pins their manifests and refers to that earlier frozen protocol. The command
checks the bindings. A request written later cannot prove when the scientific
protocol was frozen; retain the original freeze record and its supporting
source manifest separately.

Run the request with:

```sh
bes-compare benchmark --protocol selector-request.json --output-dir benchmark
bes-compare benchmark --protocol selector-request.json --output-dir benchmark --reuse
bes-compare benchmark --protocol selector-request.json --output-dir benchmark-new --save-cases benchmark-cases
```

The Python entry point is
`run_benchmark(protocol_path, output_dir, *, reuse=False, save_cases=None)` in
`base_edit_safety.population.comparison_benchmark`. It returns a manifest with
`benchmark_id`, `status` and `result_count`.

## Request and frozen protocol

The request is a closed JSON object. Unknown fields, duplicate JSON keys,
duplicate schedule values and missing required fields are errors. This example
is a template: replace every checksum with the SHA-256 of the stated file.

```json
{
  "schema_version": "1",
  "policy_version": "my-study-selector-policy-v1",
  "study": {"path": "study.json", "sha256": "STUDY_SHA256"},
  "questions": {"path": "run/questions", "sha256": "QUESTION_MANIFEST_SHA256"},
  "evidence": {"path": "run/evidence", "sha256": "EVIDENCE_MANIFEST_SHA256"},
  "frozen_protocol": {"path": "frozen-protocol.json", "sha256": "PROTOCOL_SHA256"},
  "budgets": [17, 18],
  "selectors": [
    "complete_bundle_greedy",
    "donor_first",
    "adapted_coverage",
    "canonical_bundle",
    "random_bundle",
    "exact_tiny"
  ],
  "random_seeds": [7, 19]
}
```

Paths resolve relative to the request file. `study` and `frozen_protocol` name
files. `questions` and `evidence` name stage directories; their checksums refer
to each directory's `run-manifest.json`. The evidence must be the assessment
used by the question stage. New observations require new dependent stages and
a new request.

The earlier scientific protocol must contain these fields:

```json
{
  "schema_version": "1",
  "status": "frozen",
  "policy_version": "my-study-selector-policy-v1",
  "code_snapshot": "PACKAGE_SOURCE_SHA256",
  "study_sha256": "STUDY_SHA256",
  "selector_benchmark": {
    "budgets": [17, 18],
    "selectors": [
      "complete_bundle_greedy",
      "donor_first",
      "adapted_coverage",
      "canonical_bundle",
      "random_bundle",
      "exact_tiny"
    ],
    "random_seeds": [7, 19]
  }
}
```

The scientific protocol may also contain its source-release, sample-selection,
domain, design, search, cost and validation declarations. Its
`selector_benchmark` object is closed. Its policy version and schedule must
exactly match the request, including order. `code_snapshot` is the
`source_sha256` from the comparative study's package code identity, not a Git
commit or a checksum of a separate review folder.

Budgets and seeds are nonnegative integers. Every selector runs at every
budget. `random_bundle` runs once per declared random seed; other selectors
use the study's pinned `policies.planning.seed`. Use an empty `random_seeds`
array when `random_bundle` is absent. One request may contain at most 256
selector/budget/seed runs.

## Fairness and retained results

All runs use one freshly verified canonical question set owned by the command.
Each case gets fresh search state and receives its budget, selector and seed
explicitly. Public planning and evaluation APIs still revalidate caller-supplied
inputs on every call. Each plan identity binds the actual case parameters.

For each budget, every non-exact selector must return the same mandatory
question seed, including the selected actions, resource assignments, cost and
search alternatives. A disagreement stops the benchmark. `exact_tiny` performs
joint optimization under the study's action, resource and enumeration caps;
it has no separately selected mandatory seed.

Every plan's actual actions and resource assignments are checked by the common
requirement and cost evaluator. The output directory contains:

| File | Contents |
| --- | --- |
| `plan-0001.json`, … | Complete plan audits, including charges, assignments, question states, search limits and pruning |
| `plan-records.jsonl` | Registered compact plan records bound to their full audits |
| `benchmark-summary.json` and `.tsv` | Completion, costs, controls, coverage, eligibility, runtime and memory metrics for every run |
| `question-arms.tsv` | Required and selected independent donor counts for each question arm |
| `requirement-slots.tsv` | Every catalog slot and its planned or observed assignment |
| `benchmark-request.json` | The parsed request |
| `frozen-protocol.json` | The parsed earlier scientific protocol |
| `run-manifest.json` | Input bindings, code identity, run schedule, plan identities and output checksums |

Catalog slot counts include unused donor alternatives. An engaged question/donor
bundle has at least one assigned slot. The intended-readout and control
completeness fractions concern these engaged bundles; an empty denominator is
`null`. Per-arm selected donor counts come from the complete common question
evaluation. Complete independence-group counts describe available complete
donor bundles and do not replace the question's joint matching decision.

Unknown membership and ineligibility fractions use question/donor pairs as the
denominator. Their reason counts remain separate from construction limits and
beam pruning. Evidence states remain unchanged by planning: fully planned
questions can still have no measurements. Configuration coverage is a secondary
metric with its definition retained in each audit.

The absolute weighted optimality gap is present only when a scheduled exact
solve established an optimum at the same budget. A capped or non-exhaustive
attempt does not establish an optimum or prove infeasibility. A completed
benchmark means every requested attempt returned its audit; it does not mean
every plan met its requirements or every search was exhaustive. Ties and
unfavorable results are retained.

The manifest records `request_setup_seconds` for shared input loading,
validation and source hashing. Each case's `runtime_seconds` measures planning
after that setup. Standalone public planner timing includes its own fresh input
preparation. Each plan also builds its own coverage lookup and retains per-action
coverage within that search. It does not share this cache with another case or
change donor sequence matching, costs, requirements or selector order.
Record total command duration separately, including result checks
and exports. Peak memory is
the cumulative high-water mark for that Python process, including earlier work.
It is not an isolated per-selector memory comparison. Scientific replay checks
the plan's calculations; it cannot reproduce the original elapsed time or
memory observation.

## Replay and protected output

The command checks source bytes, stage ancestry, implementation identity and
the early protocol before selection, and rechecks them before writing. It
writes the finished directory atomically, with directory mode `0700` and file
mode `0600`. Detailed outputs contain local research identifiers. Ordinary
command status reports should expose only the benchmark identity and counts.

An existing output requires `--reuse`. Reuse verifies all source and output
bytes, reruns every plan under its declared selector, and rebuilds all summary
values. A rehashed fabricated plan is rejected. Extra files, interrupted outputs
and stale results cannot be accepted as a completed benchmark. Reuse leaves the
original files unchanged. Reuse also checks all six planning limits against the
frozen study policy and requires the fresh benchmark's null prior-plan link.
Standalone plan revisions retain their explicit limits and revision history.

The schedule runs serially and publishes the complete output atomically.
With `--save-cases`, it also creates a separate protected archive. Each
`case-0001/`, and so on, contains an immutable plan and manifest after the
case passes policy, requirement, cost, fairness and input checks.
`progress.json` distinguishes completed, running, failed, incomplete and
not-started cases. A process killed without cleanup can leave a stale running
label; use actual process and termination records to resolve it. A completed
case can still be a plan that did not meet its requested questions.

The archive preserves work for inspection and independent audits. It is not a
completed benchmark and does not resume computation. An existing archive is
never overwritten. `--reuse` applies only to a completed benchmark and cannot be
combined with `--save-cases`. When case archiving is omitted, interrupted cases
remain unavailable and their internal progress cannot be inferred from absent
files. Neither mode extends an external execution budget.

A stopped monitor can leave the benchmark process running. Verify the original
process identity, command and output directory before starting another command.
Reconnect to that process where possible, keeping its original dated cutoff and
resource policy. A replacement monitor's samples cannot establish resource peaks
or uninterrupted execution during an unobserved interval. Preserve the gap and
distinguish its timing from the original measurement.

After confirmed termination, a new authorized attempt must retain the failed
attempt, bind the same request or an explicit amendment, and use a new output
directory. Such an attempt repeats the schedule. Existing completed results
remain immutable. Input/output overlap is rejected, including overlap with any
ancestor stage directory.
