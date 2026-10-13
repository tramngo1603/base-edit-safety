The comparative planner selects complete experiments under an explicit budget.
Its objective is **weighted question completion**: the sum of the frozen weights
of questions whose entire declared experiment is already supported by applicable
measurements or covered by selected actions. A fully planned question still has
its original evidence state.

Prepare a question set with `comparison_questions.prepare_question_set(study,
records)`, then call `experiment_planning.plan_experiments(question_set)`. The
study's planning policy supplies the budget, cost unit, algorithm, random seed,
beam width, expansion bound, retained alternatives and exact action bound.
Those values appear in the returned plan. Explicit keyword arguments can select
another budget or algorithm without changing the source records.

`evaluate_plan` independently recomputes a proposed selection's resource use and
requirements from the supplied source records. It does not accept a solver's
reported cost, gain or completion flags as proof. `verify_plan_revision` also
replays the declared selection algorithm and checks its complete result.

Each action has fixed resource demands or a declared list of alternative
instances. Every alternative requires an output assignment containing
`action_id`, `assignment_key` and `resource_id`. Demands that share an assignment
key on one action use the same instance. The planner sums each instance's demand
by dimension, checks its compatibility, and charges that instance once. It never
combines capacities under an informal sharing name, charges a cheaper instance
while consuming another, or invents additional batches. An explicitly unbounded
dimension can be shared without a quantity limit. Other capacities are finite.

Resources declare whether their cost activates on use, on any nonempty plan, or
even for an empty plan. The same rule applies to every selector. The output lists
all active instances, their costs and summed demand. Cost units and integer
scales must agree. Dependencies must form an acyclic graph and selected actions
must include their prerequisite closure.

Before optional selection, the five heuristic selectors use the same joint
mandatory-question search. This searches combinations of donor alternatives,
controls, actions and resource assignments across all hard-required questions.
It chooses the lowest-cost returned joint seed, with canonical identities as the
tie break, and reports its cost and retained alternatives. An optional experiment
cannot consume the money needed by a feasible mandatory seed. A failed seed
leaves the unmet required questions explicit and returns a labeled partial
result. It does not spend the remaining budget on optional work.

The selectors are:

| Selector | Optional selection after the common mandatory seed |
| --- | --- |
| `complete_bundle_greedy` | Complete an outstanding question, including its paired donors, controls and prerequisites. Rank by newly completed question weight divided by exact incremental cost, then greater weight, lower cost and canonical action IDs. Positive-gain zero-cost bundles come first. |
| `donor_first` | Group each donor's remaining eligible slot actions with their prerequisites. Rank by newly filled required slots per incremental cost, then more slots, lower cost and canonical donor ID. Incomplete questions receive no full-question credit. |
| `adapted_coverage` | Add an eligible off-target locus assay or measurement panel with the greatest new configuration coverage per incremental cost, including all prerequisites. Ties use greater gain, lower cost and canonical action ID, as in V3's cost-aware selector. Actual donor carriage is required for configuration credit. |
| `canonical_bundle` | Use the same regenerated complete-bundle constructor and select its first affordable positive-gain augmentation in canonical order. |
| `random_bundle` | Shuffle that same canonical bundle list with the recorded Python MT19937 seed and choose its first augmentation. Record the order at every step. |
| `exact_tiny` | Enumerate all eligible action subsets and their resource assignments within the declared small bounds; maximize completed weight, then minimize cost and canonical IDs. This selector jointly optimizes hard requirements and optional work. |

Complete-bundle candidates are regenerated after every selected augmentation;
shared controls can change the next incremental cost. Already satisfied
questions are reported separately and do not manufacture new selection gain.
The coverage baseline uses one equal-weight unit per canonical off-target
design/sequence configuration carried by a measured donor. Activity and clinical weights
are not inputs to that baseline.
No selector silently fills unused budget with unrelated work. A tie or worse
result remains a valid benchmark result.

Small constructors use exhaustive subset enumeration. Larger constructors
retain explicit donor alternatives, expand their first unresolved canonical
slot, close prerequisites, and reject invalid identity, capacity and budget
branches. The beam ranks fewer unresolved slots, lower incremental cost, then
canonical action and assignment identities. Its width, expansion bound and
retained complete alternatives are mandatory policy values. Existing selected
instance assignments remain fixed during an augmentation; the initial joint
seed explores their alternatives together.

The exact implementation additionally caps the action universe at 22, declared
resource instances at 128, and state visits at 2,000,000. A policy may lower these
bounds. The effective state bound is the smaller of `max_expansions` and
`exact_state_cap`. Every result records which search exhausted, which branches
were pruned, and whether the expansion limit was reached. Resource-assignment
search visits count toward the same bound.

| Status | Meaning |
| --- | --- |
| `optimal` | Exact subset and instance enumeration exhausted and established the returned objective. |
| `heuristic_complete` | The chosen heuristic finished; this is not an optimality claim. |
| `feasible_search_incomplete` | An exact request found a feasible result but hit its search bound. |
| `infeasible_proven` | Exhaustive search found no jointly feasible mandatory plan in the declared action/instance catalog, or declared unavoidable overhead alone exceeds the budget. |
| `no_feasible_bundle_found` | The bounded or interrupted search found no satisfying mandatory plan. It has not proved infeasibility. |

A declared question-construction limit is also a search limitation. Questions
excluded by that limit remain visible, and a required excluded question cannot
receive a proven-infeasible label. Exact optimization over the retained catalog
does not establish a full-universe optimum when optional questions were excluded.

Per-question planning states and evidence states remain separate. Costs and
controls belong to the historical plan. Returned measurements require a new
evidence reconciliation and a new plan revision, including failed controls.
`write_plan_revision` writes content-addressed files exclusively, checks an
existing file rather than overwriting it, and verifies any linked preceding
revision. Runtime and the Python process's cumulative peak memory are reported
separately from the deterministic plan identity.

The test-only `tests/planning_oracle.py` implements a separate exhaustive
evaluator without importing production code. Its scope is prospective synthetic
truth, explicit configuration, resolved-locus or literal single-base allele predicates, at most
16 actions, 64 resource instances and 300,000 subset/assignment states. It checks
actual donor sequences, controls, paired slots, dependence groups, prerequisites,
instance assignments, capacity and costs. Unsupported observation, historical
replication or covariate semantics are rejected rather than certified. Production
tests separately exercise those broader question contracts. The mock 18-unit
example checks the expected operations and costs.

The returned dictionary is a full audit output (`artifact_type`), with the
registered, closed `comparison_experiment_plan` record in `plan_record`. Export
that record in `plan-records.jsonl` and the full audit in `plan.json`. The record
binds the source, implementation, assignments, costs, states and an audit hash.
That hash excludes exactly `plan_record`, `plan_id`, `runtime_seconds`,
`process_peak_rss_bytes` and `memory_measurement`, preventing a circular identity.
Stage manifests additionally bind the complete saved file bytes. Storage
checks identity and historical revisions; importing or reusing a plan also
requires `verify_plan_revision` to recompute its scientific and accounting claims.
The version-2 source identity binds the execution mode, coverage basis, planning
policy, canonical scientific records and consumed source declarations/hashes.
Presentation settings and unused prediction policy do not change it. Upstream
stage verification separately binds the generation and evidence policies and
their actual output bytes. Moving identical input bytes to another directory
preserves the plan identity; the stage manifest retains their actual paths.
