# Mock computational example

This example uses mock DNA and three samples to show how BES handles missing
and unphased genotypes when a variant creates a new matching site.

Run from the repository root after installation:

```sh
bes-compare computational --config examples/comparison/computational/request.json --output /tmp/bes-computational-demo
```

The reference has one matching site. In `PHASED`, a variant changes the PAM from
`AGC` to `AGG`, creating a second matching site.

`MISSING` has no genotype call for the PAM change. `UNPHASED` has two variants,
but the data do not show whether they are on the same chromosome copy. BES
therefore cannot tell whether either sample has the new site.

The model shows the original site for all three samples, using the reference
sequence. This does not confirm their genotypes.

`report.md` and `candidates.tsv` show two sites and six rows: one row for each
site and sample. [expected.json](expected.json) lists the coordinates and
hand-calculated results used to check the software.

Use a new output directory for each run. Scientific results should match;
manifest timestamps and resource measurements may differ.
