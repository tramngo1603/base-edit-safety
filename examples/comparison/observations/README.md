# Mock observation import

This example imports a mock measurement and checks exact decimal conversion.
The extraction map has no outcome classification, so the result is `inconclusive`.
Assay quality and controls are also unknown. The original label and decimal text
are preserved.

From the repository root after installation:

```sh
bes-compare observations --config examples/comparison/observations/study.json --output /tmp/bes-observation-example
```

The source value, extraction map and normalized revision remain separate records.
See [the evidence guide](../../../docs/comparison-evidence.md) for formats, revision
history and applicability rules.
