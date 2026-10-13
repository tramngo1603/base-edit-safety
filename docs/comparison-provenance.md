# Comparison inputs and file hashes

The comparative workflow has its own contracts and cache. The V3 population workflow keeps
its existing whole-configuration run checks. A comparison cache cannot be substituted for a
V3 population stage.

`bes-compare validate --config study.json` reads local inputs, checks their declared SHA-256
hashes, and validates the comparison catalog. It does not acquire data. Its JSON result gives
the snapshot identity and record counts without listing donor identifiers or local paths.

A snapshot identity describes the complete request, input hashes and package code. A stage
identity describes only that stage's effective policies, relevant catalog records, input
hashes, upstream stage identities, actual upstream output hashes and code. Changing a budget creates a new snapshot and
invalidates planning and reporting. It leaves an otherwise identical search reusable.
Changing a control or evidence policy invalidates evidence and its descendants. A package
code change conservatively invalidates cached stages.

| Stage | Required predecessors |
| --- | --- |
| Search and intended assessment | Local source and design inputs |
| Observation normalization | Local source observations and extraction maps |
| Evidence applicability | Search, observation normalization |
| Predictions | Search and an explicit prediction policy |
| Comparison | Search, evidence applicability; predictions when selected |
| Questions | Comparison |
| Planning | Questions, evidence applicability |
| Report | Comparison, questions, planning |

The configured input role determines which stage consumes a file. A configuration cannot
omit that dependency to reuse a stale result. File hashes describe the actual local bytes;
they do not establish that an indexed window is a complete genome file.

The stage writer validates new JSONL record types through the comparison record registry,
including references to input and ancestor records. It prepares the full directory and then
moves it into place. It rejects overlapping input/output paths and existing output
directories. Use a new directory for a changed run; cached stages are immutable.

The stage reader checks current source bytes, the expected role and policy, every ancestor
manifest and output hash, and the record contracts. An interrupted, altered or incompatible
stage cannot be reused. Successful checksum checks establish that sources are unchanged; the
scientific evaluators must still reconstruct eligibility, evidence use and experimental
costs from the bound records.

`load_stage_records` returns current and ancestor canonical records after the
same verification. Shared ancestors are read once within a loading operation;
records are merged by content identity. A later call checks all bytes again.
This keeps complete CLI workflows from repeatedly loading the same ancestry
without allowing an earlier in-memory check to certify a later changed file.

Historical observation JSONL uses `observation_history` (normalization input).
The original configuration, manifest and audit of a linked prior plan use
`plan_history` (planning input). Their exact bytes participate in the consuming
stage's identity. See the [workflow guide](comparison-workflow.md) for replay
commands and immutable revision links.
