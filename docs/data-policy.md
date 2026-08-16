# Public data and disclosure policy

This repository publishes workflow code, tests, synthetic inputs, coordinate definitions, and
validation rules. It does not publish the underlying research workspace.

## Included

- source code and fixed configuration needed to inspect the workflow;
- synthetic fixtures that cannot be interpreted as biological measurements;
- public accession identifiers and source records;
- the pre-specified validation design and locked held-out assignment; and
- documentation of limitations, unknowns, and assay requirements.

## Excluded

- patient-, participant-, subject-, or donor-level records and all PHI;
- raw FASTQ/BAM/CRAM files and derived per-read or per-sample tables;
- proprietary sponsor or pharmaceutical data;
- credentials, API keys, local `.env` files, and machine-specific paths;
- downloaded publication PDFs and supplementary files whose redistribution terms are external;
- draft manuscript figures, internal review packets, generated review sites, and working notes; and
- tuning or held-out prediction-measurement results that have not been intentionally released.

Generated output remains local under `data/outputs/` unless a separate, deliberate release review
approves it. Public evidence labels describe the kind of support available; they are not safety,
efficacy, or clinical-readiness claims.
