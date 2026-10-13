# Mock Cas-OFFinder 2.4.1 inputs

These mock files test the Cas-OFFinder parser. `reference.fa` contains 400 bases;
all coordinates refer to that sequence.

The row layout follows snugel/cas-offinder
commit `9816b94c20c4cba2e79b039e1e2a6dee684b7b66` (tag 2.4.1) writes
tab-separated query, chromosome, zero-based start, matched DNA (mismatch
bases lowercase, reverse-complemented for the reverse strand), strand and
mismatch count, optionally followed by a query ID.

- `valid-seven-field.tsv`: three ID rows on both strands (intended perfect
  alpha site, beta reverse-strand single mismatch, alpha two mismatches).
- `valid-six-field.tsv` / `valid-six-field-alpha.tsv` /
  `valid-six-field-beta.tsv`: rows without IDs, linked to guides by the query maps.
- `query-map*.json`: query-to-guide mappings, including unknown
  guides, mismatched spacers and guides that share a spacer.
- `shared-guides.tsv`: two guides with identical spacers; a query map
  identifies which guide each row belongs to.
- `mixed-width.tsv`, `malformed-width.tsv`, `empty-field.tsv`,
  `bad-strand.tsv`, `bad-coordinate.tsv`, `bad-mismatch-count.tsv`,
  `bad-mismatch-text.tsv`, `bad-geometry.tsv`, `bad-pam.tsv`,
  `unknown-id-seven-field.tsv`, `duplicate-rows.tsv`, `nine-column-row.tsv`:
  rejection cases for widths, fields, strand, coordinates, mismatch counts,
  geometry, PAM, binding, duplicates and wrong declared layout.
- `uppercase_mismatch_claimed_zero.tsv`, `matching_base_marked_mismatch.tsv`,
  `extra_base_geometry.tsv`, `wrong_source_version.tsv`,
  `duplicate-query-map.json`: checks for an
  unmarked real mismatch, a matching base marked as a mismatch, a 24-base row
  for a 20+3 guide, a conflicting producer version, and a duplicated
  JSON key.
- `pam-mismatch.tsv`: a PAM mismatch (`NAG` mask, one reported mismatch)
  that the parser rejects because the mismatch count covers the spacer.
- `ambiguous-target.tsv`: an unresolved target base the parser cannot validate.
- `valid-seven-field-ngg-mask.tsv`: accepted rows with a compatible concrete
  `NGG` query PAM mask.
- `empty.tsv`, `comment-only.tsv`: rejected files with no result rows.
