# Terms used in the records

These definitions explain how these records use the terms in this tool.

- **ABE (adenine base editor):** an editor designed to change an A–T base pair to G–C.
- **Allele:** a version of the DNA sequence at a location. Here, `ref` means the reference
  sequence and `alt` means a reported alternative. A reported sequence is not proof of a donor's genotype.
- **Base editor:** a CRISPR system designed to change a base without the usual double-strand
  cut made by a Cas9 nuclease.
- **Bulge:** an unpaired base represented by a gap in a guide–DNA alignment. An alignment is
  a proposed pairing, not proof that the molecules pair that way.
- **Bystander edit:** a change near the intended base that was not requested.
- **Candidate ledger:** the JSON record of candidate sites, sources, results and review decisions.
- **CBE (cytosine base editor):** an editor designed to change a C–G base pair to T–A.
- **CH (clonal hematopoiesis) driver gene:** a gene in which particular acquired changes can
  help a blood-cell clone expand. BES flags overlaps with these genes.
- **Concordance:** agreement between a prediction and a measurement.
- **Editing window:** the guide positions assessed for an editor or model. Bases outside that
  window are outside that assessment, not proven impossible to edit.
- **Evidence level:** `predicted_only`, `analog_supported` or `empirically_measured` in a design
  record. Measurement coverage is recorded separately for each candidate DNA version.
- **Guide:** the sequence that directs a CRISPR system toward matching DNA.
- **Haplotype:** changes present together on one DNA molecule. Separate frequencies do not show
  whether changes occur together; that relationship is called phase.
- **Locus:** a location or region in the genome. HBG1 and HBG2 are separate loci.
- **Measurement coverage:** which reported DNA versions have matching, usable experiment records.
- **Nomination:** a site proposed by a search or assay for follow-up.
- **Off-target candidate:** a site outside the intended targets that may be affected and needs
  testing.
- **PAM:** the short DNA motif next to the guide-matched sequence required by the CRISPR protein.
- **Paralog:** a related gene copy. HBG1 and HBG2 have similar promoters, which some assays
  and guides may not distinguish.
- **Protospacer:** the DNA sequence matched by the guide, next to the PAM.
- **Provenance:** where a value came from, including the source file, settings and software versions.
- **QC (quality control):** checks of experiment or input-data quality, recorded with their
  pass/fail status and reasons.
- **Reference build:** a named version of the genome. The HBG definitions here use GRCh38.
- **Mock fixture:** a small input with known expected results, used to test the software.
- **Unmeasured candidate:** a proposed site without the relevant experimental follow-up. Its
  absence from a tested panel does not show whether it edits.
