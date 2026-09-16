# What I share

I share this tool's code, constructed examples and the results of my checks, with citations to
the original studies. I do not include the papers' individual experimental data points, copied
tables or detailed extracts in the public repository. A paper being available online is not
permission to redistribute its files.

## Files that can be considered for release

- Source code and the reference coordinate definitions needed to run it.
- Small constructed examples and tests, clearly labelled as software checks rather than experiments.
- Instructions, limitations and source citations.
- Aggregate results of my checks in `docs/validation.md`, the public validation document.

An aggregate result must still say what was counted and what the check does not establish.
It must not present the authors' measurements as new experiments performed for this tool.

## Files that stay out of the public repository

- Paper workbooks, PDFs, supplements, copied figures and screenshots.
- Paper-derived site lists, individual measured values, extracted rows and detailed reports
  that reproduce them, even if this tool wrote the report.
- Patient or donor records, genotypes, sample maps, raw reads and per-read or per-sample tables.
- Private validation files, source checks, review notes and prediction–measurement comparisons.
- Proprietary data, credentials, local paths and unpublished material.

Downloaded inputs and detailed generated outputs stay in approved local storage outside Git.
The private development history is not a public release. Public updates use a separately
reviewed file list; generated outputs are excluded by default.

The project code's Apache-2.0 license does not cover external software, datasets or publications.
Record local-use permission, permission to export results, and redistribution permission
separately. Constructed inputs do not remove an external tool's license requirements.
