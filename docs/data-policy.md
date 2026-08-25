# Public data and disclosure policy

This repository publishes the core analysis code, software tests, coordinate definitions, and a
short validation summary. It does not publish the underlying research workspace.

## Included

- source code and fixed configuration needed to inspect the workflow;
- small test inputs that cannot be interpreted as biological measurements;
- `docs/validation.md` as the only public validation document; and
- documentation of limitations, unknowns, and assay requirements.

## Excluded

- patient-, participant-, subject-, or donor-level records and all PHI;
- raw FASTQ/BAM/CRAM files and derived per-read or per-sample tables;
- the private `validation/` directory, detailed evaluation runner, dataset registry, protocol,
  source-overlap maps, and dataset-specific validation tests;
- paper-derived site lists, source workbooks, accession-level sample maps, and per-paper review
  notes;
- proprietary sponsor or pharmaceutical data;
- credentials, API keys, local `.env` files, and machine-specific paths;
- downloaded publication PDFs and supplementary files whose redistribution terms are external;
- draft manuscript figures, internal review packets, generated review sites, and working notes; and
- development or reserved prediction-measurement results.

Generated output remains local under `data/outputs/` unless a separate, deliberate release review
approves it. Public evidence labels describe the kind of support available; they are not safety,
efficacy, or clinical-readiness claims.
