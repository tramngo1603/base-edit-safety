# Test fixtures

These files are written by hand. They are not captured tool output. They exist so the CRISPOR and
BE-DICT parsers, the validation, and the strand-aware mapping can be tested without installing
either tool. The shape is realistic (same columns, same JSON structure, plausible ranges) but the
numbers are made up.

They are not captured-output regression fixtures. Captured output from a CRISPOR run over the
HBG1/HBG2 windows or a BE-DICT ensemble response may be added only when its filename identifies it
as captured output and a separate regression test records how it was acquired. No captured output
is included now.

`evaluation/` contains protocol-v1 known-answer tables and tool-output JSON used by the M1--M9
runner tests. They are synthetic: D1 has five measured sites and four nominations (recall 0.8), D4
has fixed bystander/outcome vectors, reference-case flags, and exactly 20 evidence-level review answers for
M8. D5 has two points with no shared source and two source-overlap points. The evaluation input check
permits fixture-mode inputs only below this tree; none of these files is an ingested validation
dataset.
