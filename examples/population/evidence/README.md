# Mock measurement example

Run from the repository root after installing the package:

```bash
python examples/population/evidence/create_example.py --output-dir /tmp/bes-v3-evidence-demo
```

Use a new output directory. The script creates all inputs and both complete plans using the
normal nomination importer, candidate-ledger builder and population pipeline. It needs no
external tools or network access. All inputs are mock data.

The first plan proposes three assays. The second includes a mock non-detection with a
detection limit of 0.1 percent for the `toy_gain` configuration in another mock donor.
It proposes two assays because that measurement supplies one matched donor–configuration
credit.

Inspect these generated files in order:

1. `inputs/source-measurements.tsv`: source guide, DNA version, assay, controls, QC and reported
   detection limit. `source-guides.tsv` preserves the full source guide/editor definitions.
2. `inputs/ledger.json`: schema-3 source ledger built from the mock BED6 nomination,
   reference and measurement TSV. Raw observations are retained.
3. `inputs/linkage.json`: checksum binding to that ledger, asserted mapping from source person
   to canonical donor, biological context, exact configuration/path and locus resolution.
   Real data need documented support for each assertion; equal labels do not establish identity.
4. `inputs/study-with-evidence.json`: complete runnable configuration with both optional
   evidence inputs and their checksums. `inputs/study-config.json` remains the baseline.
5. `with-evidence/evidence/evidence-links.jsonl`: source fields, normalized detection limit,
   identity/context decisions and any exclusion reasons. `unmatched-measurements.json` lists
   source observations that had no comparison with a current coverage unit.
6. `with-evidence/plan/assay-review.tsv`: readable targets and requirements for the two
   remaining proposed assays. Machine-readable plans and source details are beside it.

Guide IDs must agree across prepared inputs. A source name absent from the study appears in
the unmatched diagnostics, including current IDs with identical full definitions. This is a
review aid, not an automatic alias or measurement credit. Reconcile names upstream with an
explicit, documented mapping while preserving the original source export. Never equate guides
from their names alone; reused IDs with different spacers or editor definitions cannot transfer
evidence.
