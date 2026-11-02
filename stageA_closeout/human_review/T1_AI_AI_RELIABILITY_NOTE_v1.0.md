# T1 reliability note v1.0: second-model re-extraction compared with the frozen records

Date: 2026-11-01. Status: **exploratory AI–AI reproducibility comparison. Not Component A. No human reference standard.**
Inputs: `T1_SECOND_MODEL_A_sheet.csv` (9d408b4c…850f8d; provenance in `PROVENANCE_T1_SECOND_MODEL_A.md`); frozen `docrec_T1_evidence.jsonl` (4ff8e617…2b08). Nothing frozen was changed.

## 1. What this measures and does not measure

- Agreement between two AI-assisted extraction procedures on 10 studies (a prespecified hash-ranked selection from the 17 tranche studies with publication text that the owner had not adjudicated).
- It measures reproducibility between procedures, not accuracy against a human reading. Agreement does not show either extraction is correct; disagreement does not show the frozen record is wrong.
- The second model did not receive the frozen protocol or the extraction brief, so its Q3 vocabulary differs. Its labels were mapped as follows, before scoring: `assay_scoped` → `selected_assay`; `not_stated` → `none`; `not_assay_scoped` → **unmapped** (it could mean `study_enrollment` or `none`); `conflicting_assay_counts` → **unmapped** (no frozen equivalent). Unmapped cases are reported as such, not forced.
- The ten are not a random sample of the 40; nothing here is extrapolated to the tranche.

## 2. Results

| Question | Agree | Disagree | Unmapped | Denominator |
|---|---|---|---|---|
| Q3 scope | 6 | 0 | 4 | 10 |
| Q3 count (both `selected_assay`) | 5 | 0 | — | 5 |
| Q4 sampling structure | 10 | 0 | 0 | 10 |

Q4: both procedures answered `not_stated` for all ten. This is consistent with the frozen finding that no one-sample-per-participant statement was identified in these views, and it is also the least informative kind of agreement (a shared null).

## 3. The four unmapped Q3 cases, read against the frozen records

| Study | Frozen | Second model | Reading |
|---|---|---|---|
| MTBLS7260 | `study_enrollment`: 1036 | `not_assay_scoped`, cites 173 + 863 | Same substance: both decline to scope the count to the assay. Frozen graded NO_ASSAY_SCOPED_COUNT. Probable agreement; unmapped by vocabulary only. |
| MTBLS718 | `selected_assay`: 325 (from S1 only; Q2 not established) | `not_assay_scoped` | Same outcome (frozen grade UNRESOLVED, Q2_NOT_ESTABLISHED). The frozen extractor recorded an S1 count while marking assay membership unestablished; the second model declined the count. Different handling of the same gap, no effect on grade. |
| MTBLS11406 | `selected_assay`: 50 (20 + 20 + 10, S3 paper) | `not_assay_scoped` ("sample counts, not participants") | **Substantive difference.** The second model treats "50 samples" as not a participant count; the frozen extractor accepted "20 patients with OBI, 20 patients…, 10 healthy individuals" as participants. The frozen grade is COUNT (sensitivity only) and was already flagged doubtful on Q2 in the T1 review record. Affects no primary figure. |
| MTBLS13787 | `selected_assay`: 134 (S3 paper), graded COUNT_EXCEEDS_ROWS (133 rows) | `conflicting_assay_counts` (133 vs 134 vs 154 vs 198) | Same documents, same conflict. The frozen procedure records a count and lets the grader detect the mismatch; the second model records the conflict directly. Same outcome: unresolved. |

## 4. Bearing on the primary findings

- No disagreement touches x = 0 or DR-7: both procedures found no one-sample statement in any of the ten (Q4), and STRICT needs one.
- The only substantive difference (MTBLS11406) concerns a COUNT grade, which is sensitivity-only and never entered the futility calculation.
- The mapped agreement is 6/6 on Q3 scope and 5/5 on Q3 count; the four unmapped cases resolve, on reading, to three outcome-equivalent handlings and one genuine difference.

## 5. Limitations

- Ten studies; a selection, not a sample of the tranche.
- No human reference standard. The planned independent human review (Component A) was not performed; the author adjudicated 15 CONTRADICTED/STRICT records on 2026-10-31, which is human but not independent of the AI evidence.
- The second model worked without the frozen definitions; the mapping in §1 is the assembler's, made after seeing both sets of answers, and is itself a judgement.
- Component B (non-blind review of the 12 contradicted records) has not been performed.
