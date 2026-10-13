# Independent regression corpus

These 21 mock cases check deletion handling, missing calls, phase and ploidy.
The original generator and expected results are retained from the independent
review. The generator and sequence checks run independently of BES; the adapter
compares BES results with the same expected values.

The preserved module's introductory rules have these limits:

- R4 describes the cases exercised here. Its oracle does not reject an
  opposite-copy star within the same explicit phase set. Production resolution
  must honor that constraint. Recovering only the known copy applies to a true
  orphan, not a contradictory carried deletion. The core tests and the separate
  eight-case independent corpus cover this distinction.
- R6's whole-cluster wording does not describe current boundary recovery.
  Production reconstruction splits neighborhoods at declared ploidy boundaries;
  these fixtures place their tested variants away from those boundaries.
- Window probes test gap and phase status over reference-coordinate spans.
  They do not establish that every probed span exists as a physical 23-base
  sequence on a deleted copy. The separate `donor-coverage` corpus enumerates
  physical windows and checks coverage totals.

Some original wording is historical: c02's 50-base cap is no longer the research
limit; its 120-base deletion is supported within the current declared limit.
The c04 generator comment's "ALT index 0" means the first ALT, whose GT index is
1. The generator's expected-file reference means
`test_donor_coverage_regression.py`.

The 64 checks passed against the installed package at runtime revision
`04c560fb8a81d059f2f56c8c590016126592d78b`. The original corpus was preserved
rather than rewriting its oracle or expectations during integration.
