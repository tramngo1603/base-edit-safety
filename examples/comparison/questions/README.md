# Mock comparison question

This mock example builds a carrier/noncarrier comparison with an 18-unit assay plan.

`study.json` and its pinned `experiment-policy.json` declare the source inputs, scientific counts, costs, controls and finite capacities. `records.json` contains the generated canonical capabilities, specimens, question, slots, operations and resources. `labels.json` maps canonical IDs to C/N, A/B, T/L and the simple expected 18-versus-17-unit structure. `synthetic-genomes.json` records the two exact chromosome-copy sequences. Coordinates use the included mock reference.

Rebuild an independent copy with `python examples/comparison/questions/build_fixture.py --output /tmp/bes-section9`. The callable builder returns `{study, records, labels}` and accepts `alternate_carrier=True` for a distinct C2 donor. The source-only fixture supplies explicit canonical sequence records; the full genomic CLI workflow is a separate example.

See [the question contract documentation](../../../docs/comparison-questions.md) for exact endpoint pairs, relatedness, evidence reuse, control matching and resource semantics. A future fully planned comparison remains `no_measurements` until applicable completed observations are imported.
