# Comparative research reports

The report presents the recorded comparison, source observations and proposed
experiments together. It starts with the current research decision and its
limitations. A fully planned experiment still has no returned measurements until
new observations satisfy its declared controls and readouts.

## Using the renderer

`base_edit_safety.population.comparison_reporting.render_report` accepts a
resolved study, the `search`, `observations`, `evidence`, `comparison` and
`questions` output dictionaries, a full `plan` audit and its reconstructed
`question_set`. It returns a flat mapping from filenames to text or JSON values.
The comparative stage writer supplies atomic writes, protected permissions and
`run-manifest.json`.

Load the ancestors through the verified stage loader and call
`verify_plan_revision(question_set, plan)` before rendering. That call checks the
selected experiment using the existing requirement and resource evaluators.
The renderer checks canonical records, the compact plan identity, the bound audit,
configuration agreement and current source bytes. It does not select experiments
or recalculate scientific decisions. A matching audit hash alone is not a
substitute for the required plan verification.

## Exports

| File | Contents |
| --- | --- |
| `comparison-report.md` | Current decision, measured values, sequence findings, planned contrasts, costs and unresolved requirements, plus the report format and human review record section below. |
| `report-format.json` | Report export contract: `format_version: "2"`, the observation outcome contract `comparison_observation_revision:v2`, accepted observation history versions `["1", "2"]`, the producing code identity captured when the study was resolved (not recomputed at render time), and the three human records below as explicitly `not_recorded`. |
| `design-comparison.tsv` | One current summary per active design, its own PAM and mismatch/gap rules, the shared requested domain, donor denominators and full recorded search scopes. |
| `donor-design-matrix.tsv` | Current donor/design summaries, evidence availability and complete intended assessments. |
| `observations.jsonl` | Original source records, normalized revisions and state history; exact decimal text is retained. Version-2 revisions carry the explicit `adjudicated_outcome` and `adjudication_reason`; legacy revisions stay readable. |
| `evidence-applicability.tsv` | Separate identity, context, endpoint and quality assessments, with exact or related classifications and donor-credit flags. Per observation: adjudicated outcome and its reason, verbatim source outcome, detection limit and unit, missing/censoring value reason, active-revision flag, the gate-rechecked current credit and the current gate reason codes. |
| `out-of-scope-evidence.tsv` | Published alignments not recovered exactly, including positives outside the search profile, plus unresolved observation mappings. |
| `questions.tsv` | Hypotheses, arms, endpoints, donor alternatives, controls, required repeats, priorities and separate planning/evidence states. |
| `experimental-plan.tsv` | Every selected operation, specimen, design, context, arm, actual target DNA/path, endpoint, control/QC role, prerequisite and assigned resource. |
| `resources.tsv` | Each activated resource instance once, its cost, capacities, consumers and actual demands. |
| `unmet-requirements.tsv` | Recorded question, donor, action, slot and accounting gaps; incomplete mandatory questions retain the solver result and certificate. |
| `intended-assessments.tsv` | Every intended assessment, independently of whether a search nomination passed. |
| `variant-effects.tsv` | Recorded donor/reference sequence effects and source variants. |
| `sequence-annotations.json`, `annotation-manifest.json` | Supplied gene/regulatory coordinate overlaps and their source binding. These are not activity or outcome estimates. |
| `predictions.jsonl` | Separate computational prediction records, when supplied. |
| `source-evidence-summary.json`, `paired-source-descriptions.json` | Existing source summaries and matched descriptive differences. |
| `plan.json` | Full immutable plan audit, including search limits, evidence states and any infeasibility certificate. |
| `plan-records.jsonl` | Registered compact plan record binding the non-runtime audit hash. |
| `search-audit.json` | Search execution, evaluated and skipped scope, retention and source details. |
| `run-manifest.json` | Added by the stage writer: source, configuration, code, output and verified-parent binding. |

TSV cells use literal `null` for missing values. Known zero remains `0`. Structured
cells contain JSON; complete counts have separate reason columns. The renderer
uses only the updated comparison summaries, so earlier search summaries are not
counted again. Retained results from a partial search do not become complete
counts. Lost intended PAMs, indels and unresolved intended sequences remain
visible even without a nominated alignment.

Each selected resource instance is charged once, including explicitly assigned
alternatives. The report preserves the verified plan total and does not add costs
per action. Candidate-slot gaps describe donor alternatives; they do not imply
that every alternative donor must be selected. A declared limit on constructing
actions remains distinct from biological ineligibility or proven infeasibility.

## Evidence and privacy

Donor measurements require matching identity, context, endpoint and quality.
Technical repeats are counted within their source donor. Unknown design linkage stays unknown unless
an explicit executed source operation establishes a design-neutral control.
Coordinate annotations and prediction records cannot satisfy measured-evidence
requirements. Retained-callset public simulations keep empirical target-donor
credit at zero and label prospective specimens as hypothetical.

The observations table shows each active revision's adjudicated outcome with its
reason, the verbatim source outcome and the detection limit. The four states are
`observed_positive`, `tested_no_signal`, `inconclusive` and `untested`.
Inconclusive and untested records withhold conclusive empirical credit even with
a numeric value and passed QC. A no-signal result is shown as qualified by its
stated detection limit, never as proof that editing is impossible, and no-signal
values do not produce exact paired numerical difference claims. Censored or
missing values keep their typed outcome and reason but withhold conclusive
credit. Reported credit totals count only current active observation revisions
that pass the shared outcome and quality gate; a stored credit bit on a
historical, superseded or legacy version-1 record cannot add to them, and any
such retained record stays visible with its status. Current credit additionally
requires the current study's permission: the study must run in
`research_sample_planning` mode with an evidence policy that allows empirical
donor credit. Otherwise the current gate reason codes carry
`study_mode_disallows_empirical_donor_credit` or
`study_policy_disallows_empirical_donor_credit` and the record earns no current
credit.

## Report format and human review records

`report-format.json` and the matching readable-report section declare the
export contract (`format_version: "2"`, observation outcome contract
`comparison_observation_revision:v2`, accepted observation history versions
`["1", "2"]`) and the producing code identity captured when the study was
resolved; the renderer never recomputes that identity. Reviewer attestation,
funding disclosure and customer disclosure appear as explicitly
`not_recorded` because they are absent from the report inputs. The format has
no intake schema for these records. Missing review or disclosure information
remains unresolved.

By default, Markdown uses donor aliases and removes direct donor, sample,
specimen and source-subject identifiers and local paths from free text. Scientific
numbers and target DNA retain their literal values. Detailed TSV/JSON files retain
the source-bound identifiers and sequences needed for local review and replay;
the stage writer protects them. Aliases and aggregation do not guarantee
anonymity. Explicit report policies control whether identifiers or local paths
may appear in readable text.
