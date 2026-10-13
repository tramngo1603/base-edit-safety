# Comparison records and identities

The comparative workflow uses version 1 for its study and most records. Observation
revisions also support version 2 for explicit outcome adjudication. This does not
change population schema 1, candidate-ledger schema 3 or HBG design schema 7.
Record types begin with `comparison_`; IDs include their record version, such as
`comparison_observation_revision:v2:<sha256>`.

Every canonical row has `record_type`, `schema_version`, `record_id` and its named
primary ID, such as `design_id` or `observation_id`. The two ID fields must agree.
`comparison-record.schema.json` defines the closed version-1 record registry;
`comparison-observation-v2.schema.json` defines the current observation revision.
Version-1 observation records remain readable history but cannot supply current
conclusive credit. See [outcomes and migration](comparison-evidence.md#explicit-outcome-adjudication-and-legacy-history).
`comparison-study.schema.json` defines the configuration. The Python annotations
in `comparison_types.py` describe the same field structures. File readers must
validate JSON before treating it as these types.

The contract example can be checked from the repository root after installation:

```sh
bes-compare validate --config examples/comparison/contracts/study.json
python -m pytest -q tests/test_comparison_contracts.py tests/test_comparison_provenance.py
```

This mock example checks input fields and links to source files.

The configuration contains an inline `catalogs` array, local `input_refs`, the
execution mode and coverage basis, and separate policies for search, evidence,
comparison, questions, planning, reporting and predictions. Each input declares
its exact SHA-256, format and role. Input stages must match the dependency mapping
for that role; a caller cannot label a reference input as report-only to bypass
search invalidation. The study resolver checks the actual local bytes.

`policies.search.design_ids` is the explicit active comparison set. Each active
design requires exactly one effective search profile in the shared domain.
Additional, fully defined designs can remain in the catalog solely to resolve
published observations. They do not enter the active search by being present.

All object fields are declared. Unknown keys, duplicate JSON keys, invalid
versions, booleans in numeric fields, non-finite values and unresolved placeholders
are errors. An unknown scientific value uses a reason-bearing object:

```json
{"value": null, "reason": "not reported in the source table"}
```

A known value uses `{"value": 0, "reason": null}`. Thus a measured zero cannot be
created by omitting a measurement. Structural nulls such as an alignment gap have
their own exact schema and do not need an unknown-value reason. Raw source text
can retain literal labels such as `NA`; these labels do not become normalized
numeric values.

An editor activity window can be a reported interval with its basis, or
`{"value": null, "reason": "The source does not establish an editor-wide activity window"}`.
An assay's summary positions do not establish a universal editor window. An
unknown window stays unknown during target assessment and does not exclude a
search hit or an observed edit.

The registry covers the following families:

| Purpose | Record types after the `comparison_` prefix |
| --- | --- |
| Sources and genome scope | `source_artifact`, `reference_bundle`, `domain`, `locus`, `intended_locus` |
| Designs and conditions | `editor_definition`, `design`, `design_alias`, `experimental_context`, `endpoint_definition`, `search_profile` |
| Exact DNA and search | `sequence_configuration`, `sequence_occurrence`, `design_configuration`, `alignment`, `donor`, `donor_configuration`, `search_scope`, `variant_effect`, `intended_assessment` |
| Measured evidence | `source_observation`, `observation_revision`, `observation_state`, `observation_equivalence`, `evidence_applicability` |
| Computational outputs | `prediction` |
| Comparative views | `design_summary`, `donor_design_summary`, `decision` |
| Prospective experiments | `question`, `requirement_slot`, `action`, `resource`, `assay_capability`, `specimen` |
| Legacy mappings | `identity_bridge` |

`locus` represents a general genomic location, including an off-target location.
`intended_locus` adds the intended endpoint and required copy resolution. A
sequence can link to either. The HBG1 and HBG2 locus records and coordinate paths
remain distinct even if their DNA strings match.

Canonical encoding uses UTF-8 JSON, NFC Unicode normalization, sorted object
keys, no optional whitespace and finite decimal numbers. Equal numeric values
such as `1` and `1.0` have the same encoding; signed zero is encoded as `0`.
Object keys that collide after Unicode normalization are rejected. The namespace
is included in the hash input as well as in the visible ID.

Only arrays marked `x-set` in the schema are sorted. These include declared ID
sets, PAM pattern sets and resource compatibility sets. DNA paths, alignment
columns, source ALT arrays and replicate order retain their given order. Callers
must not sort arbitrary lists before hashing.

| Identity | Included scientific content | Kept separately |
| --- | --- | --- |
| Editor definition | Exact variant, Cas component, PAMs, guide geometry, conversion and declared window | Display name and source citations |
| Design | Guide, complete editor-definition ID, intended loci and desired changes | Display and source aliases |
| Reference bundle | Assembly, full component checksums, window-to-checksum mapping and transformation version | Source labels and file locations |
| Sequence configuration | Reference bundle, locus, contig, strand, complete sequence and path | Donor labels, support assertions, source variants and search scores |
| Sequence occurrence | Sequence ID, source sample, copy, phase block/status and source version | Human donor labels and source row aliases |
| Design configuration | Design and sequence-content IDs | Intended/off-target reporting role |
| Alignment | Design configuration, exact columns and PAM sequence/path | Native row aliases; mismatch/bulge lists are checked against the columns |
| Source observation | Immutable source file ID, exact locator, measurement key and replicate key | Parser output |
| Normalized revision | Source observation plus extraction, mapping and schema versions and extraction-map checksum | Other revisions of the measurement |
| Question | Hypothesis, requirements, scientific policy, independent counts, controls, reuse and priority | Generated slot IDs, plan state and evidence state |
| Action | Actual specimen/design/context/endpoint/target/replicate/batch operation | Question membership, dependencies and resource assignment |
| Requirement slot | Question, actual donor, arm, endpoint, design, context, control role and repeat position | Allowed assignments |

Excluding generated references from question, slot and action identity avoids
circular hashes. Changes to the full stored records still invalidate stage
outputs through the stage manifest. A record's biological ID alone is
not a sufficient cache key. In particular, a resource reassignment or changed
control dependency must be bound by the planning stage even when the action's
operation ID remains stable.

Source files include their source identity, release and URI, as well as
their immutable byte hash and locator. Equal bytes from two different source
origins do not prove that their experiments are the same. Reference bundles
separately bind actual component checksums, so a citation alias does not change
the underlying DNA identity. Regional indexed windows contain their own
interval-to-checksum mapping and a null whole-reference checksum with a reason.
They cannot be presented as a complete-reference-file checksum.

An auxiliary `reference_validation` input is a separately pinned reference-window
JSON used only to check VCF REF alleles outside the searched domain. It has a
search-stage dependency but cannot supply searched or callable sequence. The
reader checks the exact search bundle separately and rejects auxiliary windows
that overlap the declared search domain. `published_alignments` inputs have an
evidence-applicability dependency and bind source alignment replay mappings.

Summary records retain observed nomination counts separately from complete
counts. A partial or unknown search has null complete counts with reasons,
including when the retained count is zero. Donor summaries reference every
intended locus assessment, exact search scopes, nominated configurations and
alignments. Their stated support and completeness must agree with those records.
They also distinguish evidence that has not been assessed, no evidence in current
sources, contextual evidence and exact donor evidence. Evidence references are
normalized observations, never predictions. The immutable stage manifest binds
each summary and decision to its current source and evidence state.

A parser update creates another `observation_revision` for the same
`source_observation`. `observation_state` records the complete revision history
and exactly one active revision. `active_observations()` returns only active
measurements. An observation equivalence record can document that different
source files contain the same actual measurement. Equal values alone never
create equivalence. Superseded normalizations and equivalent motivating
measurements cannot fill an allowed repeat slot.

The raw value, unit and exact source locator remain immutable source records.
Numeric normalization must preserve that value or declare its conversion factor
and offset. The contract checks these calculations. Source importers still have
to preserve exact source decimals: a normalized value can use strict finite decimal
text when the exact result cannot be represented as a JSON floating-point number.
Decimal text zero follows the same detection-limit rules as numeric zero. Importers still have
to verify the actual file locator and extraction against the original file.
The contract does not infer a detection limit, donor identity, independent unit,
assay success or copy resolution from a reported fraction.

Predictions have their own record type, model and request/response identities,
window construction basis and terms metadata. Foreign keys for empirical
observations accept only `observation_revision`; a prediction cannot satisfy
one of those references. Exact empirical credit also requires matching donor,
design, context, endpoint and sequence, with passed controls and QC. A measured
zero without a positive reported detection limit earns no reassuring donor
credit. Retained-callset simulations retain zero empirical donor credit.

Donors carry explicit independence-group and relatedness policies. Specimens do
not create additional independent donors. Resource records declare integer cost,
unit and scale, compatibility, activation, pool bounds and capacity dimensions.
Unbounded capacity needs an explicit justification. A matching price or informal
resource name is not a sharing rule. The record layer rejects an individually
oversized demand, including repeated entries for the same resource/dimension,
and incompatible resource assignment. Resource timepoint compatibility uses the
exact known context value and unit, written as `48 hour`, for example; unknown
timepoints cannot meet a specified timepoint rule. Pool instance indexes are
unique, with consistent maximum counts and per-instance costs. Summed selection capacity,
mandatory controls, pairing, independent counts and bundle completion are also
checked by the planning evaluator. Planning checks describe the proposed
experiment; observed control and QC results are assessed separately.

The Python entry points are:

```python
from base_edit_safety.population.comparison_contracts import (
    make_record, record_identity, validate_record, validate_catalog, validate_study,
    active_observations, bridge_legacy_record, restore_legacy_record,
)
```

`make_record(kind, **fields)` fills the envelope, primary IDs and design definition
hash, then validates the row. `validate_record()` checks one row's structure,
identity and local scientific constraints. `validate_catalog()` additionally
checks foreign-key types, cross-record consistency and dependency cycles, and
returns an ID-to-record index. `validate_study()` adds configuration, active
design, input-role and policy checks. It does not read or acquire input files.

To construct a question and its slots, compute `record_identity()` from its
scientific fields first, use that ID in the slots, then pass their IDs to
`make_record("question", ...)`. Generated slot IDs are excluded from the question
hash. A catalog must contain both sides of the final references before it is
valid for storage or use.

Legacy records continue through their existing readers. An identity bridge
preserves the complete original JSON, checksum, family/schema version, original
ID, new ID, mapping version and explicit reason. Restoring a bridge returns the
old record without changing nulls, source links or IDs.
The fixture `tests/fixtures/comparison/legacy-roundtrip.json` is an actual V3
synthetic sequence record. Its mapped comparative sequence retains the exact
DNA, path and reference-file checksum but receives a new namespaced content ID.
Generalized gapped alignments must likewise receive new comparative identities.

## Moving an existing study to the comparative workflow

Keep the original configuration and output directories. Existing population,
candidate-ledger and design commands still read their original schemas. An
identity bridge preserves a record for traceability; it does not supply the
additional information required by a comparison study.

Create a separate comparison configuration with explicit guide/editor designs,
per-design search profiles, intended loci and a common genomic domain. Declare
the coverage basis, exact input hashes and unresolved calls. Add the measured
evidence's actual context and sequence links, then declare the new question's
independent donor counts, pairing, controls, specimens, resource capacities and
costs. Do not translate a legacy configuration-coverage count into completed
experimental questions or reuse its nominal noncarriers without checking DNA.

Run the new stages into new directories. Their namespaced identities and stage
manifests preserve the new scientific meaning without changing historical
records. A budget-only change can reuse matching completed comparative search
and evidence stages; a changed package source requires newly bound stages under
the conservative code-identity policy. The [workflow guide](comparison-workflow.md)
describes revision and replay, and the offline example supplies a complete
artificial configuration.

The checked-in corpus exercises every registry family and includes adversarial
tests for unknown fields, types, identities and scientific reference boundaries.
These checks establish the contract layer. Live external search, real source
extraction, complete experiment selection and the public comparative benchmark
have separate acceptance gates.

## Proposed ABE7.10 supported-domain profile

Opt in through `policies.evidence.supported_domain` in an existing version-1 study:

```json
{
  "profile_id": "abe710-spcas9-hek293t-plasmid-v1",
  "context_id": "<existing comparison_experimental_context ID>",
  "endpoint_ids": ["<existing comparison_endpoint_definition ID>"]
}
```

Use actual content IDs, with the existing catalog and source bindings. This optional
field does not change historical record identities or activate itself for older studies.
The packaged `protocols/abe710-spcas9-hek293t-plasmid-v1.json` is versioned with the
source. Its hash and bound context/endpoint IDs enter every affected stage dependency.
Unknown profile versions and caller-supplied approval flags fail validation.

The exact editor fields are `variant_name: ABE7.10-SpCas9`, `editor_class: ABE`,
`cas_component: SpCas9`, guide-oriented A-to-G conversion, a 20-base spacer and a
three-base 3-prime NGG PAM. Context fields are `organism: Homo sapiens`,
`cell_type: HEK293T`, `cell_state.value: transformed cell line` and
`delivery_method: plasmid transfection`. These are explicit declarations, not inferred
aliases or verification of a construct. Keep dose, exposure and collection time in
separate existing quantity fields; unknowns remain reason-bearing values.

Only source-bound intended and guide-dependent unintended editing fractions at the
read level are admitted as target endpoints. Bind their required context and allele
resolution. Other endpoints, editors or materials abstain. Search uses both strands,
zero alignment bulges, all admitted matches and mismatch/total-edit caps equal and
at most six. A lower cap is a narrower declared search, not proof of six-mismatch
coverage. Positions 4–8 are a declared annotation heuristic; no other position is
removed as safe. Simple genomic indels retain the existing exact sequence path and
50-base allele ceiling; they do not enable alignment bulges.

This profile currently executes only the existing local `comparative_local` v1 path
in `synthetic_test` mode on mock autosomal inputs. Cas-OFFinder 2.4.1 is the
proposed subsequent native route; parser compatibility alone does not enable it.
Unsupported editor, bound context, endpoint, PAM, strand, mismatch/bulge, prediction,
backend or sex-chromosome requests return a structured CLI `abstained` result (exit 2)
before execution. Structurally invalid configurations remain `invalid_input`.
Real input modes cannot activate this profile yet. Do not relabel real data as artificial
or a non-diploid genome as diploid to pass a declaration check.

The existing `search-audit.json` and Markdown report keep three claims separate:

- Computational support means the declared request can be checked with mock inputs.
- Coverage completeness remains incomplete. A complete search of the declared intervals
  is reported separately. Full primary-reference coverage and alternate X/Y ploidy support
  are unverified; declared callability files alone do not establish completed coverage.
- Biological applicability abstains while scientific approval and the 48-hour/72-hour
  exposure conflict remain unresolved. Missing dose/exposure/timepoint fields add reasons.

Existing applicability assessments, decisions and current report credit carry these
reasons. Stored positive credit cannot bypass the current profile in comparison,
planning or reports. Prospective mock plans can still be generated within the
bound target. Unsupported action contexts/endpoints are ineligible. Raw measurements,
source text, all four adjudicated outcome states, source assay-quality checks and paired
source descriptions remain intact.
No profile run can claim `release_complete: true` in this version.

Scientific approval requires the qualified scientist's accepted construct/material,
exposure and endpoint contract plus reviewed source adjudication and sensitivity limits.
The data, custodian and rights decisions are recorded with the release permissions
note. They cannot be supplied by a boolean in a study. A later reviewed profile/version
must implement any accepted change without rewriting old observations or protocols.

Before native integration, pin the Cas-OFFinder executable, build/runtime and notices;
compare exact/mismatch/NGG/both-strand boundary cases against the tiny local oracle;
check query binding and zero-based positions, duplicate loci, indel-derived windows,
source-copy projection and row accounting; test empty success versus failed/truncated
execution. Then review reference/alternate contig, callability and ploidy accounting.
No installation, native run, real processing or V8 replay is authorized by this profile.
The paired interval for fixed-budget recall improvement remains a subsequent evaluator
dependency, with mismatch ranking mapped as `score = -mismatches`.
