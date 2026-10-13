# Local comparative sequence analysis

The local search reconstructs the same declared donor domain for every active
design. It then applies each design's exact spacer, PAM geometry, mismatch limits
and DNA/RNA bulge limits. An intended site is assessed independently of whether
the search admits a hit. A lost PAM, changed sequence length or unresolved phase
therefore remains visible in the intended assessment.

Intended assessments retain every fixed-locus alignment allowed by the declared
gap limits, plus an ungapped diagnostic when the lengths permit one. Donor and
reference admission consider gap alternatives as well. A variant-effect record
uses the canonical admitted alignment when one exists; the other alternatives
remain separate canonical records. An excluded ungapped diagnostic cannot hide
a nomination created through an admitted gapped alignment.

Run the mock example after installing the repository:

```sh
bes-compare search --config examples/comparison/synthetic/study.json --output /tmp/bes-comparison-synthetic
python -m pytest -q tests/test_comparison_sequences.py
```

The example has two designs with the same spacer and different PAM policies,
two mock HBG copies, forward and reverse targets, a PAM gain/loss, an
intended genomic insertion, unphased heterozygotes, missing calls and a jointly
normalized multiallelic deletion. `expected.json` contains independently derived
expected nominations and exact copy sequences.

The search writes three output files. `comparison-records.jsonl` contains canonical
sequence, occurrence, design-configuration, alignment, intended-assessment,
carriage, variant-effect, scope, summary and decision records. `search-audit.json`
retains actual execution metadata, exact hit/window links, exclusions and the
source-equivalence audit. `editability-annotations.json` records substrate and
potential bystander positions for each alignment. The pipeline binds all three
to its source/code manifest.

The donor summary's observed counts describe retained nominations. Complete
counts are null if any relevant source or backend scope is unresolved. A
fully evaluated zero remains a known zero. No nomination count selects a winner.
The initial decision records a sequence hypothesis or insufficient comparable
evidence; empirical measurements are assessed by the separate evidence stage.
Unreported editor activity windows remain null with their reasons. Window
membership is an annotation and never a filter on observed evidence or search
admission.

## Genomic source binding

`comparison_inputs.prepare_inputs(study)` uses a checksum-pinned V3 population
source configuration to load the reference, donor map, genotypes and callability
policy. This preserves the existing V3 genotype, phase, unsupported-footprint
and ordered-ALT rules. It does not run the old single-guide search or change its
schemas.

The comparative config declares the V3 JSON as a `source_transformation` input.
Every nested file read by that JSON must also have an explicit comparative input
reference with the correct role and actual SHA-256. An automatically read FASTA
index must be pinned as well. An unbound nested path fails before search. The
canonical donor/sample membership must exactly match the source TSV.

Complete FASTA references bind their complete file hash. Regional bundles must
contain exactly the same indexed intervals as the pinned search reference. Each
canonical window source is separately pinned and its actual sequence is checked
against the aggregate window. Supported source formats are explicit
`reference_windows` JSON and the retained UCSC sequence-response JSON. Citation
hashes alone cannot establish sequence equality.

Some retained VCF records have REF spans outside the search domain.
The optional `reference_validation` role supplies explicit same-assembly reference
windows for those checks. These inputs cannot overlap the search domain or add
searchable sequence. A REF span crossing a domain edge may use adjacent acquired
pieces from both sources. Validation requires continuous, agreeing sequence;
a gap or conflicting base fails. This joining is used only for REF checks and
does not alter the reference used for donor reconstruction. The exact auxiliary
pieces used are recorded with their source hashes.
Unsupported carried alleles and missing or filtered genotypes retain their
unknown footprints. A known reference genotype at an unsupported alternate does
not silently become a carried structural variant.

## Prepared windows and mapping

`PreparedInputs.iter_windows()` streams one bounded source interval and donor
copy at a time. It separates missing sequence, incompatible phase blocks and
uncallable regions. Adjacent callable intervals are merged before constructing
windows. Single unphased heterozygotes retain an unordered pair; multiple
heterozygotes require a common phase block to produce a complete local sequence.
Inserted bases require support for the REF span establishing the insertion.

Each emitted window carries an immutable content ID, exact ACGT sequence,
genomic-forward per-base path, reference identity, donor/sample/copy, source
genotypes in original ALT order, phase state, source locators and support basis.
Insertion points retain their true boundary and offset. Reference windows have
no donor ID. The reconstruction audit preserves every source and skipped region
and proves which donor samples and policies produced the submitted windows.

Sequence content is separate from donor occurrence, design configuration and
alignment identity. Equal DNA at HBG1 and HBG2 retains different locus IDs.
Different donors can share a sequence while retaining every copy and source
occurrence. Different valid gapped alignments remain distinct. A genomic indel
does not itself imply a guide-target bulge.

`comparison.normalize_search(prepared, hits, execution)` maps local or native
hits through those exact paths. It reconstructs and checks submitted window IDs,
then verifies every target/PAM index and base and reapplies the exact per-design
profile and shared domain. Native derivations retain both their source row and
an explicit transformation note. Native evaluated-window IDs and skipped-window
reasons control actual evaluated scope. Truncation, a nonzero exit, incomplete
source equivalence or a skipped relevant window prevents a complete claim.

## Bounds and checks

The local backend contract is `comparative_local`, version `v1`, with its declared
digest. Calling it with another backend's profile is an error. The source
configuration's donor, variant, reference, materialization and reconstructed-base
ceilings remain active. Per-design row and source-window variant ceilings also
apply. A high-density window is explicitly skipped; intended assessment remains
independent of that skip.

`local_search(study, limits=LocalLimits(...))` defaults to 600 seconds and 2 GiB
maximum process resident memory. It checks elapsed time and resident memory
during reconstruction and alignment enumeration and before returning complete
outputs. It retains actual elapsed time, CPU time, peak process resident bytes
and the enforced limits. Alignment-state and row overruns raise a resource error
without returning a complete result.

The focused suite checks exact expected admissions, independent intended losses,
copy paths, multiallelic order, unknowns, arbitrary PAM geometry, alternative gap
placements, skipped scope, deterministic biological results under source-row
reordering, source binding and resource failure. Native execution, source-panel
equivalence and the larger public benchmark retain their own independent gates.
