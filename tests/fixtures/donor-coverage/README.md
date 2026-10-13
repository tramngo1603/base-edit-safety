# Independent donor-coverage fixtures

These eight mock cases check donor reconstruction against independently written
sequences. `fixtures.json` stores the cases; `oracle.py` calculates the expected
results without using BES.

The cases cover supported and contradictory spanning deletions, bounded orphan
stars, multiallelic and overlapping records, missing calls, local phase, haploid
X/Y-style inputs, declared PAR/non-PAR copy changes, zero-copy regions, both
strands and indels. Every copy count comes from the fixture's explicit ploidy
intervals and source statement. Fixture copy labels are 1-based; the regression
adapter explicitly converts them to BES's 0-based labels.

The adapter compares every evaluated physical 23-base window, including indel
windows, with the declared sequence strings. It joins reconstructed neighborhoods
with reference context, applies gaps and local phase checks, and compares
unphased copies without assigning a long-range phase. Whole-contig exclusions
yield no evaluated windows; excluded regions never count as losses. Windows do
not cross declared ploidy-region boundaries.

Case 4 uses the explicitly selected `vcf_implicit_phase_set` policy because its
`|` genotypes have no PS field and its truth asserts a joint copy ordering. This
follows the implicit phase-set convention in the
[VCF specification](https://samtools.github.io/hts-specs/VCFv4.3.pdf). Other cases
use `explicit_ps`; the public HG003 workflow also retains `explicit_ps`.

From the repository root, run the adapter with:

```sh
python -m pytest -q tests/test_donor_coverage_oracle.py
```

To run the independent sequence-only enumerator directly:

```sh
cd tests/fixtures/donor-coverage
python oracle.py
```

The enumerator reads only the hand-declared haplotypes and masks. It does not
parse VCF or reconstruct genotypes. Expected totals are 18 NGG candidates and
18 hits under the fixed 20-nt, both-strand, ungapped, at-most-six-mismatch scope.
