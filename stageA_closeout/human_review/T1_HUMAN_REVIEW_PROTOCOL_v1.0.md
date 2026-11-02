# T1 human review of AI-assisted extraction: protocol v1.0 (frozen before selection)

Date: 2026-10-31. Approved design: owner decision of 2026-10-31 (10 of 17, plus a separate review of 12).
Scope: tranche T1 of documentary recovery. This review **does not change any frozen T1 record, grade or result**. Any disagreement is reported as a separately identified sensitivity analysis.

## Component A: independent re-extraction (blind to AI answers)

**Sampling frame (frozen, n = 17).** These are the T1 studies with at least one usable publication reading view in the sealed universe (`d6803c26…dbd2`) that were **not** among the 15 records the owner adjudicated:

MTBLS524, MTBLS718, MTBLS916, MTBLS1584, MTBLS1991, MTBLS2774, MTBLS3233, MTBLS6038, MTBLS7260, MTBLS8390, MTBLS9586, MTBLS10304, MTBLS11406, MTBLS12655, MTBLS13608, MTBLS13787, MTBLS13960.

**Selection rule (frozen before computation).**
- key(acc) = SHA-256 of the UTF-8 string `docrec-T1-human-review-v1|` + acc, as lowercase hex.
- Sort ascending; select the first 10.
- The rule depends only on the accession string. It was written and hashed before any key was computed.
- In the manuscript it is described as "a prespecified hash-ranked selection of 10 of 17 eligible studies", not as a probability sample.

**Blinding statement (for the manuscript, verbatim).** "The human reviewer was blinded to the AI-generated answers and supporting quotations but had prior knowledge of the overall study classifications."

**Materials.** For each selected study, only the sealed files in `stageA_closeout/docrec_out/store/<ACC>/`:
- `S1_i_Investigation.txt`;
- every `*.view.txt` listed on the answer sheet.

Until the sheet is complete, the reviewer does **not** open:
- `docrec_T1_evidence.jsonl`, `docrec_out/extraction/`, `docrec_T1_results.csv`;
- `DOCREC_T1_EXTRACTION_REVIEW_v1.0.md` or `DOCREC_T1_RESULT_v1.0.md`;
- the candidate file (it holds row counts);
- earlier chat messages about these studies.

**Questions.** These are the categories of extraction brief v1 (`b53bc3df…9a22`), §2, unchanged. The reviewer fills `T1_HUMAN_REVIEW_A_sheet.csv`.

- **Q3:**
  - `scope` ∈ {`selected_assay`, `study_enrollment`, `none`};
  - `count` (integer, when scope ≠ none);
  - `evidence_quote` (verbatim) and `source_file`.
  - `selected_assay` requires that the text ties the number of distinct participants to the samples analysed for the study's chosen matrix and platform (given on the sheet).
- **Q4:**
  - `answer` ∈ {`explicit_one`, `explicit_repeated`, `not_stated`};
  - `evidence_quote` and `source_file` (required unless `not_stated`).
  - Contradictory evidence is recorded as the answer the text supports most explicitly, with `contradiction=yes` and an explanation in `notes`. No new category is created.

**Agreement rules (frozen).**
- **Q3 scope:** agreement if the categories are identical.
- **Q3 count:** assessed only where both reviewer and AI chose `selected_assay`; agreement if the integers are equal.
- **Q4:** agreement if the categories are identical.
- **Reporting:**
  - k/10 for Q3 scope and for Q4, reported separately, never combined;
  - Q3 count as k/m, where m is the number of studies in which both chose `selected_assay`;
  - every disagreement listed, with both answers and both quotes.
- No kappa (n = 10 is too small for it to be stable). No extrapolation to all 40 T1 studies.
- No AI record is edited because of this review.

## Component B: non-blind review of the 12 final DOC_CONTRADICTED records

**Cases.**
- Repeated sampling (9): MTBLS89, 92, 200, 2291, 3800, 10634, 10972, 12326, 12372.
- Non-participant samples (3): MTBLS419, 2205, 12861.
- Three of the twelve have no publication text (MTBLS10972, 12326, 12861); for these the review uses the investigation file.

**Procedure.**
- The reviewer reads the source texts in `store/<ACC>/` and records a verdict on whether explicit evidence supports the classification:
  - repeated sampling of the same participants within the chosen assay;
  - multiple specimens from one individual;
  - a non-participant (cultured or non-human) sample set.
- The verdict is recorded before re-reading the AI evidence quotes.
- **Verdict** ∈ {`confirmed`, `disputed`, `unresolved`}, with a verbatim quote and the reason for any `disputed` or `unresolved` verdict.

**Report.** The number confirmed, disputed and unresolved, and the exact reasons. This is labelled **non-blind**: the reviewer adjudicated these records on 2026-10-31.

**Outside both components.** The three adjudication-rejected records (MTBLS24, 1122, 15433) are already DOC_UNRESOLVED, and their non-blind adjudication is on record.

## Sequence

1. This file is hashed.
2. The selection is computed and recorded.
3. The answer sheets are issued.
4. The owner completes A, then B.
5. The assembler scores A against the AI records and tabulates B.
6. A results note is written, labelled as a sensitivity analysis.
