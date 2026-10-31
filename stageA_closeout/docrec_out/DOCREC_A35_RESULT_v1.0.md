# Documentary recovery, stratum A (A35): V1–V4 pre-verification result v1.0

Date: 2026-10-31. Protocol `DOC_RECOVERY_PROTOCOL_v0.2.md` §7 (frozen), `docrec.py` v0.1.1.
Status: **final** (the code has no adjudication step for A35). Owner acceptance is pending (§8).

## 1. Result

| Status | n | Studies |
|---|---|---|
| VERIFIED (V1–V4 all pass) | **3** | MTBLS14865 (V3 = 368), MTBLS12684 (144), MTBLS15004 (85) |
| FAILED (any check fails) | 8 | MTBLS148, 9021, 6260, 13382, 12824, 15165, 11048, 13319 |
| UNVERIFIABLE | 8 | MTBLS13634, 14426, 20, 14358, 13827, 14421, 13136, 14376 |

- **N_A35,verified = 3 of 19.**
- **N_A35,documentary recovered = NOT ASSESSED.**
  - Protocol §7 sends verified studies to Q1–Q7 "in their own later tranche", but `docrec.py` v0.1.1 implements no such tranche.
  - Two of the three verified studies (MTBLS14865, MTBLS15004) have no publication in the evidence set, so rule 8 (S1-only evidence) would bar their recovery.
  - Only MTBLS12684 has a named publication. Its view states that each participant "provided a fasting plasma sample" (controls) or "provided a preoperative plasma sample" (patients). Whether that meets Q4 `explicit_one` is **UNKNOWN** until graded under Q1–Q7. It is reported here as a sensitivity observation only.
- **No effect on P4V or DR-7.** A35 can contribute at most 19 studies, and DR-7 concerns stratum M only.

**Why studies failed** (provisional reading of the FAILED rows):
- **Non-participant samples:**
  - MTBLS11048: cell lines.
  - MTBLS9021: mouse serum rows inside the assays.
  - MTBLS148: pooled urine aliquots.
  - MTBLS13319: one pooled plasma, aliquoted for method development.
- **Repeated sampling:**
  - MTBLS6260: 8 cyclists, 3 conditions × 3 time points.
  - MTBLS13382: dried blood spot and serum from the same 4 people.
  - MTBLS12824: muscle biopsies before and after exercise.
  - MTBLS15165: stored platelets sampled over days 0–7.

**Why studies were UNVERIFIABLE:** mostly there is no stated sample count (V3). MTBLS14376 also has an unresolved organism/protocol conflict.

## 2. Chain of custody

| Object | SHA-256 |
|---|---|
| Evidence universe A35 (seal signed by the owner before extraction) | a76e8c45e30df638e1312763daeed74f37ea37e5752cee5435125311856c704d |
| `docrec_authorisation.txt` after the A35 seal line | 27f5208f31c1465a4b7524ddc5ebf96d61b23d83658a3f91078bda3d3af18a0e |
| Extraction brief A35 v1 (fixed before any record) | c0990e22d957c9e67bebed258a792c82455c0789c3f0a450d5c69392e4fdcb38 |
| `docrec_A35_evidence.jsonl` (19 records) | b586d80605d8c41f8e4f18d23863d9859a4b5443d31feb1100978bb230885b79 |
| `docrec_A35_results.csv` | c4e4fb9b89118d32ae0928f502f2ec80144f2c4d6f76a8b5f70136583a117678 |
| `docrec_A35_report.json` | 78b575a685a30d13e139447885ceb729e3178d45303ffac6313b3518f640f9d5 |

- **Retrieval.** One fetch run (2026-10-31T10:50:34Z) made 124 requests, all HTTP 200.
- **Sealing.** No owner-supplied texts were registered (DR-9 optional papers for MTBLS11048, MTBLS20 and MTBLS6260 were not obtained). That choice was made before the seal and before any reading.
- **Copy check.** All 98 usable files were copied after the seal and matched the sealed hashes.
- **Grading agreement.** The owner's `grade --tranche A35` output (3 / 8 / 8 / 0 invalid) equals the assembler's independent grading.

## 3. How the records were made

- **Extractors.** Four Claude sub-agents of the assembler, 4–5 studies each, worked under one hashed brief. They read only the sealed S1 files (investigation, sample and assay files) and the two S2 views.
- **Blind to row counts.** The extractors were blind to the census row counts. Their validator used a universe copy without `n_rows`.
- **Privacy rule.** No sample names, source names or individual-level rows were quoted. Computed tallies appear in notes as numbers only.
- **No reviewer edits.** No record was edited after grading.

**Reviewer checks:**
- **The three VERIFIED records.** Their V3 counts equal the deposits' sample rows, which the extractors did not know. MTBLS15004 is the weakest: its count equates 85 patients with 85 samples via the word "baseline".
- **Duplicate suspicions from equal row counts were false.** MTBLS13382 (833 rows) is a dried-blood-spot method study, not HuMet (MTBLS89). MTBLS14358 (325 rows) is a colorectal-cancer stool study, not MTBLS718.

## 4. Deviation from the protocol text (fixed before reading)

- **The unit.** The census selected no assay for A35 studies: no biological rows survived rule M2, and the candidate file has no matrix or platform for them. The brief therefore applied V1–V4 to the biological samples in each study's sample file, i.e. the set its assays measure. Where assays measured different subsets, the largest shared set was used (MTBLS13382, MTBLS20, MTBLS15165).
- **The V3 counting rule.** V3 accepts only counts stated in text. Counting rows in a sample file does not satisfy it, which is how the frozen `check_count` works.

## 5. Findings (from the sealed A35 files)

1. **Template default mislabels participant samples as blanks** [CERTAIN, recomputed].
   - In 18 of 19 deposits, every sample row carries `Characteristics[Sample type] = experimental blank`. The 19th, MTBLS12824, has 94 of 95 rows so labelled.
   - This explains why the census QC-removal rule (M2) left no biological rows.
   - The label contradicts other fields in the same rows: disease-group factors and participant sample-collection protocols.
2. **Sample-to-measurement linkage defect in MTBLS14421** [CERTAIN, recomputed; counts only].
   - In both assay files, 50 of 53 rows carry a *different* sample's label in `Extract Name`, `MS Assay Name` and `Normalization Name`.
   - In 30 rows that label belongs to the other group: 15 control rows carry patient labels and 15 patient rows carry control labels.
   - From the deposit alone, a re-user cannot tell which measurement belongs to which participant group.
3. **One study, two deposits.** MTBLS14358 and MTBLS14376 share the title "Decoding the Human Gut Bacterial Plasmids in Colorectal Cancer" and one described cohort.
4. **Little publication text.** Only 2 of 19 A35 studies have a publication text in the evidence set.

## 6. Evidence status

- **CERTAIN:**
  - 3/19 verified, 8 failed, 8 unverifiable under the frozen code;
  - findings 5.1–5.4.
- **SUPPORTED:** the individual failure reasons (extractor reading, reviewed for the verified records only).
- **UNKNOWN:** whether any verified study is documentation-recoverable (Q1–Q7 not implemented for A35).
- **Not claimed:** anything about the true independence of these cohorts.

## 7. What this does not change

- Stage A L1 (= 0).
- The T1 result and DR-7.
- P4V v0.12.
- No external outcome was opened, and no measured value was read.

## 8. Owner acceptance

Accept by appending to `docrec_authorisation.txt`:
`DOCREC A35 RESULT ACCEPTED <date> record=<SHA-256 of this file> owner=<name>`
