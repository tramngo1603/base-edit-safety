# Glossary

This page defines specialized terms used in the records and documentation. The definitions describe
how this repository uses each term; they are not clinical definitions.

- **ABE (adenine base editor):** an editor that can convert an A-T base pair to a G-C base pair
  within its supported editing window.
- **Allele:** one DNA sequence observed or predicted at a genomic position. An edited allele differs
  from the reference sequence at one or more positions.
- **Base editor:** a CRISPR-derived system designed to change a DNA base without making the standard
  double-strand cut used by many CRISPR nucleases.
- **Bystander edit:** an edit inside the intended target region that is not one of the requested
  changes.
- **CBE (cytosine base editor):** an editor that can convert a C-G base pair to a T-A base pair
  within its supported editing window.
- **CH (clonal hematopoiesis) driver gene:** a gene associated with expansion of a blood-cell clone
  when particular acquired variants occur. This tool reports gene overlap only; it does not estimate
  a person's risk.
- **Concordance:** agreement between a prediction and a measurement. The validation documents use
  this term for a defined metric, not as a general claim that the tool is correct.
- **Editing window:** the positions within a guide-targeted sequence for which an editor or predictor
  supports an edit.
- **Evidence level:** one of `predicted_only`, `analog_supported`, or `empirically_measured`. The
  label states the type of support attached to a record.
- **Guide:** the sequence that directs a CRISPR system to a matching DNA region.
- **Haplotype:** the complete set of edits present together on one predicted or measured DNA
  sequence. Position-by-position edit rates do not show which edits occur together.
- **Held-out set:** data reserved for a final evaluation and not used to choose thresholds, weights,
  or rules. The tuning set is the separate data used during development.
- **Hand-verified reference set (golden set):** cases for which a reviewer records the expected result
  directly from a primary source. The validation code uses the field name `golden_set`.
- **Locus:** a named region of the genome. In this project, HBG1 and HBG2 are distinct loci.
- **Off-target site:** a genomic location other than the intended target that a guide or editor may
  affect. A predicted site is a candidate for testing, not a measured edit.
- **PAM (protospacer-adjacent motif):** a short DNA sequence next to the guide-matched region that the
  CRISPR protein requires.
- **Paralog:** one of two related genes produced by gene duplication. HBG1 and HBG2 are paralogs with
  similar promoter sequences, so some assays and guides may not distinguish them.
- **Protospacer:** the genomic sequence matched by the guide, next to the PAM.
- **Provenance:** the source, tool version, code version, reference build, and other details needed to
  trace where a value came from.
- **Protocol lock (freeze):** the point after which the dataset assignment, thresholds, and code
  version cannot change for that protocol version.
- **QC (quality control):** technical checks on source identity, file integrity, sample mapping, and
  required metadata. Passing QC does not establish biological performance.
- **Reference build:** a named version of the reference genome. This project uses the GRCh38 primary
  assembly for its HBG coordinates.
- **Stratum:** a group reported separately during evaluation, such as results for one assay, cell
  type, or editor.
- **Source-overlap group:** held-out points whose source paper also informed a rule being evaluated.
  These points are reported separately because they are not fully independent of rule development.
- **Synthetic fixture:** test data created for software checks. It is not a biological measurement.
