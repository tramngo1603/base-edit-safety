# Donor-aware search and proposed assays

The population workflow searches reconstructed DNA in explicit intervals, links compatible
measurements, and selects donor–site assays under a fixed budget. The standard-library core
requires no external search tool. Specimen availability is supplied by the caller.
Population schemas use version `1`; design and candidate-ledger schemas use `7` and `3`.

## Reproduce the examples

The [mock workflow](../examples/population/evidence/README.md) includes a reference,
variants, donors, assay capabilities and costs. Outputs containing donor
linkage are created with restricted local permissions. Never commit real genotype inputs,
donor matrices or protected output tables.

The [measurement integration example](../examples/population/evidence/README.md) creates a
complete mock source ledger, person/configuration linkage and plans before and after
evidence. It runs without optional genomic tools or network access.

Run the deliberate greedy counterexample independently:

```bash
bes-donor-panel --request examples/population/planner/request.json \
  --manifest examples/population/planner/planning-manifest.json \
  --algorithm greedy_cost --output-dir /tmp/bes-v3-greedy
bes-donor-panel --request examples/population/planner/request.json \
  --manifest examples/population/planner/planning-manifest.json \
  --algorithm exact_tiny --output-dir /tmp/bes-v3-exact
```

Greedy chooses `B1,C3`, spends 4 of 6 cost units and covers 2 requirements. The independent exact
solver chooses `A1,A2,C3`, spends 6 and covers 3. Activating donor A once makes two assays jointly
attractive; the one-action greedy rule misses that combination. The exact solver is limited to
12 eligible actions and checks the results of the other selectors.

## Declare the study

Use [study-config.json](../examples/population/synthetic/study-config.json) as a complete small
example. Paths are relative to that configuration. The
[study schema](../src/base_edit_safety/population/schemas/population-study.schema.json) describes
serialized types; the validator also checks scientific relationships. Unexpected keys, duplicate
JSON keys, non-finite numbers, unresolved `REPLACE` values and boolean costs are rejected.

| Mode and coverage basis | Meaning |
| --- | --- |
| `synthetic_test` / `synthetic_truth` | Mock genomes and study settings for software checks. |
| `benchmark_simulation` / `retained_callset_model` | Assume reference at unlisted positions, conditional on the retained variant callset. Never fill over retained missing or filtered calls. |
| `benchmark_simulation` / `callable_supported` | Use complete sequence support; availability and assay costs may be a declared scenario. |
| `research_sample_planning` / `callable_supported` | Require sequence support, available appropriate specimens, permitted use, and explicitly feasible assays. |

A variants-only VCF cannot establish complete reference-call evidence. The conditional model
does not bound all possible human variation. The NYGC integrated resource omits singleton and
other filtered variants, contains related individuals and uses inferred phase. See the
[official release directory](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/)
and [source study](https://pmc.ncbi.nlm.nih.gov/articles/PMC9439720/).

The protected donor TSV has these columns:

```text
donor_id genotype_sample_id source_kind family_id specimen_id specimen_context_id availability use_permitted activation_cost simulation_eligible sampling_stratum
```

Use tabs, `true`/`false` booleans, explicit integer costs and blank optional family/stratum fields.
One row represents one specimen. Several rows may share a canonical donor ID and its genotype
sample ID, with distinct specimen IDs. Genotype mapping, source kind, family, activation cost
and sampling stratum must agree for that person. A genotype sample cannot map to two people.
Availability, permitted use and context are checked separately for each specimen. Genotypes,
frequency denominators and donor credit count each person once; activation is charged once
even when several specimens are selected. Family membership does not establish independent
validation units.

For supported reference calls, supply the following tab-delimited columns:

```text
donor_id contig start_0based end_0based_exclusive evidence_type source_id
```

Evidence types are `explicit_reference_calls` and `callable_mask_with_source_contract`. The
source entry must describe quality and nonreference-allele representation. Overlapping rows
must be reconciled upstream; a depth mask alone is insufficient.

The core reads bounded plain or gzip VCF. It preserves original multiallelic ALT order and GT
indexes, site FILTER, sample FT, GQ/DP policy and PS. Phased calls without PS use the explicitly
selected VCF convention. PS never phases a slash genotype. Multiple local heterozygotes need a
common source/sample/contig phase block. Windows outside an unresolved footprint can still be
evaluated. Unbounded unsupported alleles conservatively invalidate all materialized scope for
affected donors. A complete call carrying only supported alleles uses the record's own contig
and literal REF span, even when an uncarried ALT in that record is symbolic. Accepted complete
`0/0` calls do not mask other donors.

Public variants require a recorded bcftools normalization version. The acquisition helper
preserves original rows and ALT-index mappings without splitting multiallelic GTs. General
structural variants, gVCF conversion, sex chromosomes, mixed ploidy and de novo phasing are not
implemented. Allele counts use accepted complete calls at each variant; configuration frequencies
use donors evaluable for that complete target. Neither is a patient-prevalence estimate.

Literal multiallelic indels are checked as complete records. An individual ALT may be able to
move within a repeat even when the unsplit record is jointly normalized. That alone does not
invalidate the record or justify changing its original ALT indexes.

If left alignment reaches the start of a cropped reference, preparation fails with an explicit
`normalization-failure.json`. The indel needs more genuine upstream context; a crop-limited
placement is not accepted as globally normalized. A chromosome's actual start is not a crop.

Reference input is an indexed FASTA or checksum-bound windows with explicit global coordinates.
Small synthetic FASTA files also work without an index. VCF `##reference` must match the declared
assembly. A source preparation step must document any conversion of an upstream reference name.
Reference alleles are checked before search; no REF/ALT swap or build conversion is inferred.

## Sequence search and identity

The local backend scans all supported 23-base windows on both strands: a 20-base guide plus a
three-base IUPAC PAM, with a declared mismatch limit and no alignment bulges. Genomic indels can
create an ungapped guide match and are distinct from alignment bulges. No candidate coordinate
list, score threshold or editable-base filter limits primary nominations.

Paths retain reference positions, insertion boundaries and offsets. Deletion junctions retain
their skipped reference interval. An entire target inside an insertion has a zero-width reference
locator and an explicit insertion path. Reverse-strand paths follow guide orientation. These
paths must not be passed through legacy interval-only site identifiers.
Insertion checks use target bases and the original REF span of the variant establishing the
insertion. An unrelated right-flanking base is not required for phase or callability. A target
entirely inside an insertion can be evaluated at the end of a contig.

Candidate hashes bind the reference input, guide/editor definitions, oriented sequence and
path. Donor IDs, scores, search thresholds and padding are excluded. Identical HBG1 and HBG2
sequences remain distinct coordinates; only declared intended loci receive a control role.

Targets spanning requested boundaries are excluded. Provide enough materialized context for
the study you intend to run. Search covers the supplied intervals; provide padding and
chunks explicitly. Default guards are 32 donors, 100,000 reference bases,
5,000,000 reconstructed bases, 100,000 materializations, 100,000 variant rows and 250,000 hit rows.
A cap failure produces an explicit partial execution without a complete result. Unknown calls
can coexist with a completed scan of the remaining supported windows.

`bes-population search` uses the built-in bounded search. The separate external-contract
validator checks declared execution and row-retention status; backend-specific scope checks
remain separate. The CRISPRme highest-CFD importer is sample-free and retains no all-donor
configuration coverage.

The search audit also counts paired and unpaired windows before PAM or mismatch filtering.
For each donor, guide and partition, it takes the union of exact paths on the reference and
both donor copies. A paired window requires all three sequences to be evaluable under the
selected coverage basis. Indel paths without an exact reference counterpart remain unpaired;
their nominations remain in the search results. The protected audit records compressed ranges
of paired reference window starts, oriented-window denominators and exclusion reasons. A
start range describes window starts, not a merged target interval. Aggregate reports omit
individual linkage.

Nomination-change records use the same exact-path pairing rule as the window audit. An indel
inside either donor copy prevents a paired comparison for the reference path, including a
heterozygous insertion whose other copy still matches reference. Those comparisons retain an
explicit indel reason and do not count as paired PAM losses.

## Match evidence and declare assays

An optional schema-3 candidate ledger and a checksum-bound measurement-linkage JSON connect
existing observations. Linkage supplies canonical people, exact configurations and paths,
source assertions, locus resolution, and structured cell/delivery/condition context. Unknown or
ambiguous person identity cannot add donor multiplicity. Different source identities for the
same person must resolve to one canonical donor ID.

The linkage document has `record_type: measurement_linkage`, `schema_version: "1"`,
`ledger_sha256`, `people` and `links`. Each person contains `source_person_id`,
`canonical_donor_id`, `identity_status` (`resolved`, `unknown`, `ambiguous`) and `source_assertion`.
Unresolved people require a null canonical donor. Each link identifies `site_id`,
`source_version_or_citation`, `assay`, `source_record_id`, `source_person_id`,
`biological_context_id`, `context`, `configuration_id`, `target_path`, `locus_resolved` and
`source_assertion`. Missing links preserve observations as unresolved.

A usable non-detection needs exact compatible DNA, guide/editor, locus, context, controls,
passed QC and a positive detection limit for editing frequency in supported units. `percent`,
`percentage` and `%` accept values up to 100; `fraction` and `proportion` accept values up to 1.
Units ignore case and surrounding spaces. The evidence assessment normalizes these to a
fraction while preserving the source value and unit. Unknown units remain unresolved, and
out-of-range values cannot earn non-detection credit. No universal sensitivity cutoff is
assumed. A prespecified `assay_policy.max_detection_limit_fraction` can require a lower limit;
this changes the coverage and resolution identities. Omission or null means no additional
sensitivity threshold. The assay declaration and laboratory review must support the proposed
sensitivity requirement.
Unit conversion and the inclusive threshold comparison use the decimal spelling of recorded
numeric values. For example, 0.07 percent and 0.0007 fraction both meet a maximum of 0.0007.
No rounding tolerance relaxes that maximum. A limit too small to represent as a positive
normalized JSON number remains unresolved.

Source-ledger guide definitions
must match the current spacer and editor, class, PAM and editing-window definition. Matching
display IDs or a ledger checksum alone cannot establish compatibility; missing definitions
remain unresolved. The evidence record retains the source definition. Pooled or unknown-allele assays do
not resolve both copies. REF results cannot cover ALT DNA. Conflicting outcomes remain separate
records. A verified ordinary experimental reference sequence can now be represented without
inventing a donor allele or measurement quality; the separate representation replay preserves
the historical source ledger and reports zero new donor-measurement credits from that repair.

Guide IDs must also agree across prepared inputs. `unmatched-measurements.json` retains source
observations without a current coverage-unit comparison, explains why, and lists any current
guide IDs with the same complete guide/editor definitions. These suggestions supply no credit
and do not establish aliases. Prepare a documented source-name mapping before import and
preserve the original export. Matching names with different definitions still fail evidence
checks. See the complete mock example above.

The assay policy declares its guide IDs and locus scope. Physical planning requires explicit
exact assay loci. An HBG2-only capability cannot cover HBG1. Simulations can explicitly assume
coverage of all admitted loci. Capability and context IDs hash their actual definitions;
changed cell, delivery, condition or resolution declarations change the relevant identity.
Unknown feasibility goes to `requires_review` and earns no coverage. Each generated action
assesses one donor, specimen and configuration; there is no implicit free multiplexing.
The protected plan TSV identifies the selected specimen for each resolved assay action.
The additional `assay-review.tsv` joins donor/specimen, coordinates, strand, guide and target
spacers/PAM, sequence path, linked variants and donor calls, assay/context, resolution and
feasibility declarations, prior evidence gap, and selection rationale. Variant links describe
source records associated with the configuration; inspect each selected donor's call. These
are proposed experiments, with requirements to review, not completed assays.

## Plan and audit costs

For each configuration/context unit, the objective counts distinct donors with usable existing
evidence or proposed eligible assays, up to the required donor count. Priority weights are
integers, not risk probabilities. Two chromosomes, repeated source rows and technical replicates
in one person cannot increase donor multiplicity.

Cost equals selected donor activation costs plus selected assay costs. Declared control costs
must be in assay packages or a fixed study overhead, deducted before optimization. The report
restores that overhead. This is a declared cost model; procurement, shared prerequisites and
laboratory scheduling are outside scope.

`greedy_cost` recomputes marginal gain and activation cost after every selection. Ties use higher
gain, lower added cost, then canonical action ID. `greedy_coverage`, `canonical_order` and
`random_order` use the same eligibility and cost rules. Random input order is recorded along
with the seed and Python version. Stratified sampling and family-aware external evaluation
require a separate protocol and are not silently inferred here.

The [request](../src/base_edit_safety/population/schemas/donor-panel-request.schema.json) and
[result](../src/base_edit_safety/population/schemas/donor-panel-result.schema.json) schemas are
strict. The CLI also requires `planning-manifest.json`. It checks upstream bytes and identities,
recomputes units, action eligibility, costs and evidence applicability, and rejects fabricated
coverage even if someone updates the request checksum. The mock fixture
route accepts a synthetic declaration in place of biological source files.

Every result records exact objective fractions, existing and planned donors separately, unused
budget and unmet requirements. Empty results distinguish no remaining demand, no affordable
eligible action and no positive gain. No action is converted into a measurement record.

## Public-data protocol

[public-protocol.json](../examples/population/public-protocol.json) freezes the NYGC release,
eight-sample selection rule, bounded chr11 interval, existing HBG −198 guide, NGG/four-mismatch
profile, cost assumptions, seeds 0–99 and budget grid. It is a transductive retained-callset
simulation. Its full universe remains in every selector comparison, including units with no
eligible assays. Zero new nominations is a valid result.
The helper validates this exact packaged protocol and constructs all scientific settings from
its implementation. Only resource limits come from `--example-config`; editing template
weights, filters, phase rules, assay settings, selectors or seeds cannot alter the study.
Other study designs need a separately declared protocol. Benchmark selectors must be a
nonempty, unique list; an empty planning universe still produces the requested comparison rows.

Rare-variant reporting uses fixed descriptive bins from accepted complete diploid calls in
the declared donor subset: unknown AN, AC zero, singleton AC, and nonsingleton ALT frequencies
below 1%, from 1% to below 5%, or at least 5%. It reports retained ALT counts and each selector's
all-unit and actionable-unit coverage for linked variants in those bins. A unit can link to
several bins; those overlapping counts must not be added as distinct configurations. Unlinked
configurations have their own category. These bins do not alter weights or selection and do
not establish source-wide frequency, singleton completeness or patient prevalence.
With eight diploid donors, nonsingleton frequencies below 5% cannot be resolved.
A broader study needs predeclared regions and a justified donor selection rule. Preserve null
results and compare selectors under identical inputs rather than tuning a study for an advantage.

Acquisition is an explicit separate command. It requires network access and the optional
genomics dependency; ordinary pipeline commands access only local files:

```bash
python -m pip install ".[genomics]"
python -m base_edit_safety.population.public_data \
  --protocol examples/population/public-protocol.json \
  --example-config examples/population/synthetic/study-config.json \
  --output-dir data/inputs/population-v3-public
```

Then run the README stages with the generated local `study-config.json` and a separate output
root. The helper records actual checksums, versions and normalization arguments. It translates
VCF positions into a real reference window for bcftools and restores global coordinates after
normalization. It retains unsupported records and checks any boundary-overlapping REF bases
against separately fetched real sequence. A window checksum is never called a whole-genome
checksum. Original source rows, mapping and pedigree metadata stay in protected local storage.

Acquisition errors fail the required gate. They are not skipped or relabeled as a successful
public benchmark. No controlled-access data or laboratory work is requested.

The separate BCL11A/SpCas9 rs114518452 test uses the published locus and a retrieved UCSC
reference window, with explicitly synthetic genotypes. It checks the reverse-strand TGC→TGG
PAM change and three spacer mismatches. It is a known-case nuclease regression, not HBG
base-editor validation. See [Cancellieri et al.](https://doi.org/10.1038/s41588-022-01257-y).

Reproduce the separate public carrier count with the packaged, fixed chr2 protocol:

```bash
python -m base_edit_safety.population.known_case \
  --output-dir data/inputs/population-v3-carrier-recount
```

This explicitly acquires the exact allele with all source sample columns and needs the
genomics extra and network access. Alternatively, recount the previously saved bounded,
plain VCF entirely locally:

```bash
python -m base_edit_safety.population.known_case \
  --source-vcf /path/to/known-case-source.vcf \
  --expected-sha256 2698d96021709d5281bc5f0ab707c3aea4b507ecb4cbd8dddd99464bfab2ef1f \
  --output-dir data/outputs/population-v3-carrier-recount
```

That checksum identifies the historical local source subset, not the full upstream chromosome.
Fresh extraction may serialize its header differently; each run records its actual checksum.
The protocol requires GRCh38. The reader checks `##reference`, target-contig assembly metadata
and any declared chr2 length; incompatible or unrecognized reference declarations fail. Plain
GRCh38/hg38 names and reference filenames containing those names are recognized. A matching
local declaration establishes a declared assembly, not the file's data collection or release.

The historical subset has no reference-build header. Its exact checksum is bound to a recorded
public acquisition in the packaged protocol, supported by the [IGSR source collection
documentation](https://www.internationalgenome.org/announcements/3202-samples-at-high-coverage-from-NYGC/).
Only those recorded bytes or a direct acquisition by this helper can use that source contract.
Other files without assembly metadata fail; an arbitrary `--expected-sha256` proves byte
identity alone. Reacquire from the documented source instead of inserting a guessed build.

`summary.json` uses the generic `known_case_genotype_check` record type and separates
`expected_provenance` from `source_verification`. The latter reports direct acquisition,
a match to recorded acquisition bytes, or an unverified local source. Unverified local files
receive null source release and URL fields. Assembly declarations and acquisition evidence
are distinguished by `assembly_basis`. These checks are input-history and metadata checks, not a
whole-reference sequence verification or a cryptographic proof of the upstream data's accuracy.
The summary contains aggregate counts and the declared quality policy. Source sample columns
and `protected-carriers.json` stay protected locally. Dot FILTER is retained as unfiltered,
not upgraded to PASS. Complete diploid calls count by original ALT index; missing, partial and
sample-filtered calls are excluded. No depth, quality, callability or editing result is inferred.
Duplicate ALT entries and sample values exceeding the FORMAT fields are malformed inputs and
fail the run, including when that sample would otherwise be excluded by FT.

## Output and limitations

The [population row schema](../src/base_edit_safety/population/schemas/population-record.schema.json)
describes genotype audit, sequence paths, candidates, hits, carriage, evidence, units and actions.
Each stage binds the input hashes, code identity and upstream stage hashes in `run-manifest.json`.
Changing an input or policy requires rebuilding downstream outputs. Writes use temporary files
or directories, reject input/output overlap, and require `--overwrite` to replace existing output.

The default Markdown report omits individual donor/sample IDs, memberships and local paths.
The report separates usable positive/non-detection links, conflicting donor–unit pairs,
excluded links and unmatched source measurements, with reason counts. Multiple links can refer
to one donor; these counts must not be read as independent people.
`donor-site-plan.tsv`, `assay-review.tsv`, `uncovered-configurations.tsv`, unmatched source rows,
JSONL matrices and run manifests stay protected
locally. Small aggregates can still identify people; the report is not automatically
anonymous or approved for sharing.

Successful exit codes are `0`; invalid inputs return `2`, required public acquisition failures
return `3`, and resource-cap executions return `4`.
