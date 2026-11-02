# T1 human review, Component A: result v1.0

Date: 2026-11-01. Protocol: `T1_HUMAN_REVIEW_PROTOCOL_v1.0.md` (799865e8…f509). Sensitivity analysis only; no frozen record was changed.

## 1. Provenance and attestation

- Reviewer: the author (Anas Enoch). Attestation, verbatim, 2026-11-01 23:02 local: *"I read the 23 documents and wrote these Q3/Q4 answers and notes by myself."*
- Files: `T1_HUMAN_REVIEW_A_completed.numbers` (39f1145a0fe6147ed2cac5f5d99ccb6359b4040d53d6af4c2c238dda92e59466; the reviewer's original) and its CSV export `T1_HUMAN_REVIEW_A_completed.csv` (93465965e7345643474e53a5391db3819beafe2746934cb5c44ba9906cbb6f8b).
- Checks by the assembler before scoring: all ten rows use the frozen categories; every quotation occurs verbatim in a permitted sealed view; every cited source is on the row's permitted list. The reviewer's quotations and notes differ from the second-model sheet (4 of 20 quotation cells coincide, on sentences any reader would select).
- **Disclosures required in the manuscript:**
  1. The reviewer is the author, not an independent reviewer.
  2. Before this review the author had seen a second AI model's answers for the same ten studies (`T1_SECOND_MODEL_A_sheet.csv`), and a relabelled version of them, but not the original extractor's answers or quotations.
  3. The author had prior knowledge of the overall study classifications (as the protocol anticipated).
  4. The sheet was assembled by a script after the review; its `reviewed_utc` values are script timestamps (all within 0.2 s), so per-study completion times are not available.
- The assembler did not correct, suggest or discuss any answer before scoring.

## 2. Agreement (frozen rules, protocol §A)

| Question | Agree | Disagree | Denominator |
|---|---|---|---|
| Q3 scope | 6 | 4 | 10 |
| Q3 count (both `selected_assay`) | 4 | 0 | 4 |
| Q4 sampling structure | 10 | 0 | 10 |

Q4: both answered `not_stated` for all ten (a shared null; consistent with the frozen finding that no one-sample statement was identified in these views).

## 3. The four Q3 scope disagreements

| Study | Frozen | Author | What differs | Frozen grade |
|---|---|---|---|---|
| MTBLS718 | `selected_assay`: 325 (from S1 alone; Q2 not established) | `study_enrollment` | Both found that the documents do not tie a cohort to the chosen assay. The frozen extractor still recorded the S1 total (245 + 80) as assay-scoped; the author declined. The frozen record's Q2 = not established already produced UNRESOLVED. | UNRESOLVED (Q2 not established) |
| MTBLS11406 | `selected_assay`: 50 (20 + 20 + 10) | `none` | The author reads "50 samples" as a sample count, not a participant count; the frozen extractor accepted "20 patients…, 20 patients…, 10 healthy individuals". Substantive difference in reading. The T1 review record had already flagged this study's Q2 as doubtful. | COUNT (sensitivity only) |
| MTBLS2774 | `selected_assay`: 431 (330 + 101) | `study_enrollment` | Both read 330 + 101 individuals; the author would not assign the total to the chosen assay because the deposit does not say which cohorts it holds. The frozen record assumed both and was graded COUNT_EXCEEDS_ROWS (rows = 330). The author's reading is the more conservative one. | UNRESOLVED |
| MTBLS13787 | `selected_assay`: 134 | `none` | Same documents, same conflict (S1: 133; paper: 134 of 154). The frozen procedure records the published figure and lets the grader detect the mismatch; the author records the conflict as no established count. | UNRESOLVED |

Pattern: in all four, the author was more conservative than the frozen extractor about assigning a count to the selected assay. In none did the author find a count or a sampling statement that the frozen extractor missed. Three of the four frozen records were already UNRESOLVED; the fourth is a sensitivity-only COUNT grade.

## 4. Bearing on the primary findings

- No disagreement concerns Q4, so none bears on the STRICT criterion or on x = 0 / DR-7.
- Under the author's readings, the grades of the four studies would be UNRESOLVED (718, 2774, 13787 unchanged in effect; 11406 would move from COUNT to UNRESOLVED). The STRICT + COUNT sensitivity figure would fall from 5 to 4; the primary figure is unaffected.
- No reading by the author would raise any grade.

## 5. Limitations

- Ten studies from a prespecified hash-ranked selection of 17; agreement is not extrapolated to the tranche.
- The reviewer is the author and was partially exposed (§1). This is weaker than the protocol's intended independent reviewer and is reported as such.
- Component B (non-blind review of the 12 CONTRADICTED records) has not yet been performed.

## 6. Comparison with the second-model sheet (exploratory)

Where the second model's labels could be mapped, its answers and the author's coincide on all ten Q3 scopes (after the author's own mapping of the four unmapped labels) and all ten Q4 answers. This is reported for completeness; it does not add independent evidence, since the author had seen the second-model sheet.
