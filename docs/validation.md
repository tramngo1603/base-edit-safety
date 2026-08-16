# Validation status and safeguards

Validation is in progress. The public repository contains the evaluation code and synthetic
fixtures with fixed expected results, but no held-out comparison between predictions and
measurements.

The v1 plan separates public measured datasets into a tuning allocation (D1-D4) and a single
held-out allocation (D5). That assignment was locked on 2026-08-14 before held-out outcome review.
Changing the assignment requires a new protocol version; tests pin the role map and lock date.

Additional gates remain intentionally closed until their evidence is complete:

- every dataset must pass checks covering its source, reuse terms, checksum, sample mapping,
  guide/editor mapping, reference-genome build, and overlap with development inputs;
- metric thresholds and the code release under test must be frozen before prediction-measurement
  comparisons are viewed;
- a non-fixture held-out run requires a clean checkout at the frozen validation tag;
- the single-run record is written before outcome values are loaded and cannot be overwritten; and
- failures are reported under the pre-specified categories; held-out data cannot be used to adjust
  the method after the run.

Accordingly, the current release makes no measured performance, clinical-safety, or regulatory
clearance claim. The versioned source document is
[`validation/validation-protocol.md`](../validation/validation-protocol.md); the locked split and
dataset registry record its change history. See [evaluation runner](evaluation.md) for the input
contract, fixture-only CI mode, and generated report structure.
