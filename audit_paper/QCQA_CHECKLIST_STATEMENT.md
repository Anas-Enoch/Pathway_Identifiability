# Quality-control / quality-assurance checklist statement

Manuscript: Identifiability of independent sampling units in publicly deposited metabolomics data: a repository-level audit and documentary assessment (Enoch, submitted to *Metabolomics*).

*Metabolomics* asks authors to supply a QC/QA checklist with each submission. The checklist is written for studies that generate metabolomics measurements (sample handling, instrument QC, pooled QC samples, blanks, internal standards, batch correction, data processing). This manuscript generated no measurements:

- No biological samples were collected, stored, extracted or analysed.
- No mass-spectrometry or NMR data were acquired or processed; no raw, processed or derived measurement files were downloaded or opened (the protocols forbid it; Supplementary Methods S4).
- The study units are public repository accessions (metadata files, sample and assay sheets, investigation descriptions) and the published documents linked to them.

The checklist items that do apply, and where each is addressed:

| Checklist theme | Applies? | Where addressed |
|---|---|---|
| Study design and sample size justification | Yes (design-based) | Methods 2.2, 2.4 (census universe; random tranche of 40 of 371 under a fixed NIST beacon pulse; hypergeometric bound) |
| Pre-specification of analysis | Yes | Methods 2.1–2.6; Supplementary Methods S1, S3; protocols and code frozen by hash and signed before execution |
| Data provenance and versioning | Yes | Data Availability (commit identifiers); Supplementary Methods S4 (sealed evidence inventories, SHA-256 of every document) |
| Blinding | Yes (as applicable) | Methods 2.4, Supplementary Methods S3 (value-blind reading views: Results sections excluded from extraction) |
| Quality control of the extraction | Yes | Methods 2.7, Results 3.4, Supplementary Methods S5 (machine verification of every quotation; author review components A and B; disclosed deviations) |
| Sample collection, storage, preparation | Not applicable | No samples handled |
| Instrument QC, pooled QC, blanks, internal standards | Not applicable | No measurements acquired |
| Data processing, normalisation, batch correction | Not applicable | No measurement data processed |
| Metabolite identification / annotation levels | Not applicable | No metabolites reported |
| Data deposition | Yes | All protocols, code, locks, evidence records and grades public in the cited repository |

The author will complete the journal's checklist form marking the measurement items "not applicable" with the reason above, and pointing the design, pre-specification, provenance and extraction-QC items to the sections listed.
