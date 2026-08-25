# Example workflow

This walkthrough shows how a guide/editor analysis moves through the project.

## 1. Configure the external tools

Copy `.env.example` to `.env` and set the local BE-DICT and CRISPOR paths. The `.env` file stays
local.

## 2. Enumerate candidate edits

```bash
python -m base_edit_safety.pipeline.enumerate_guides
```

The command reads the fixed GRCh38 HBG1/HBG2 windows, obtains guide and sequence-match candidates
from CRISPOR, and obtains edit predictions from BE-DICT. It writes
`data/outputs/enumeration_HBG1_HBG2.json`.

Each record includes the guide, editor, locus, strand, intended and bystander alleles, off-target
categories, HBG1/HBG2 duplication flags, and tool details. Missing or unsupported predictions are
marked as `not_evaluated`.

## 3. Build review records

```bash
python -m base_edit_safety.pipeline.design_records
```

This writes `data/outputs/design_records_HBG1_HBG2.json`. Each design record groups its edit
predictions, evidence labels, locus checks, gene-overlap flags, and suggested follow-up assays.

Generated files stay under the ignored `data/` workspace. The small examples in `tests/fixtures/`
are available for code review and automated tests.
