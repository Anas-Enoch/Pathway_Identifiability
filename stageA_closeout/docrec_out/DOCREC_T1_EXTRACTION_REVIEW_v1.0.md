# Documentary recovery T1: extraction and review record, v1.0

Date: 2026-10-31. Tranche T1, 40 studies (simple random sample of 40 of the 371 stratum-M candidates; salt NISTv2:2-2005958).
Status: extraction complete; provisional grades computed; **owner adjudication of 15 CONTRADICTED records pending**. Futility is formally computed by `docrec.py grade` only after adjudication. As section 4 shows, no adjudication outcome can change x.

## 1. Chain of custody

| Item | SHA-256 |
|---|---|
| Lock `docrec_lock.json` (signed 2026-10-28T22:43:50Z) | 678933f642aafc205130b7f0d1243a86cd4191c7143ebd2f2b9179a927b859d2 |
| `docrec.py` v0.1.1 (frozen; also the validator's engine) | 419776ec8c378d5c275e1fcc7d046b179224d8d084362695541607079d520eb2 |
| Evidence universe (seal signed by the owner 2026-10-31, before any reading) | d6803c2666ebdc54f793659c4aec77b60c6c0c4f0a09c71f6b02363df074dbd2 |
| `docrec_authorisation.txt` after the seal signature | 76b4d0bc61a4d9b885b9e49a6019c8d92fb558f0da017ddbfd1ca36e3a25fcb2 |
| Extraction brief v1 (fixed before any record was written) | b53bc3df4e4887e21b23cf2f35af353a83e7562aa94444c37e4bfe1fc5a19a22 |
| Validator `validate.py` | 76c1ad4bce5b337de564db74cd4dfddc9371b3d8bbaef1f6772bf7cf725cf8c8 |
| `docrec_T1_evidence.jsonl` (40 records, rank order) | 4ff8e6177d3707fcf0728ede4f56ec31e031be2d2f2c8773d52ac11e686c2b08 |

- The 83 usable files (40 investigation files and 43 reading views) were copied from the owner's Mac after the seal. All 83 matched the sealed hashes.
- The owner's first `grade` run (before the evidence file existed) passed the lock, seal-signature, order and store-equals-universe checks. It stopped only at the missing evidence file.

## 2. How the records were extracted

- **Extractors.** Five Claude sub-agents of the assembler, 8 studies each, all working under the same brief. Each read only its studies' sealed files. No full XML/PDF, manifests, other studies or web were opened. Extractor names are recorded in each record.
- **Blind to row counts.** The extractors were not given the census `n_biological_rows`, and their validator used a copy of the universe with `n_rows` removed. Q5 (count = rows) was therefore decided only by the grader. The protocol does not require this; it changes no grading rule.
- **Validation.** Every record passed the frozen code's schema, verbatim-quote and count-form checks before review.
- **No reviewer edits.** The assembler made **no change to any extracted record** after seeing grades. All reviewer findings are recorded below, and in the adjudication sheet for CONTRADICTED records.

## 3. Provisional grades (frozen `grade_record`, real row counts)

| Grade (reason) | n | Studies |
|---|---|---|
| DOC_RECOVERED_STRICT | **0** | — |
| DOC_RECOVERED_COUNT (sensitivity only) | 5 | MTBLS916, 3233, 9586, 11406, 12655 |
| DOC_CONTRADICTED (REPEATED_SAMPLING) | 11 | MTBLS10972, 200, 3800, 1122, 2291, 10634, 24, 12372, 92, 12326, 89 |
| DOC_CONTRADICTED (NOT_PARTICIPANT_SAMPLES) | 4 | MTBLS15433, 419, 12861, 2205 |
| DOC_UNRESOLVED (Q2_NOT_ESTABLISHED) | 6 | MTBLS12722, 718, 10304, 12058, 8390, 1991 |
| DOC_UNRESOLVED (NO_ASSAY_SCOPED_COUNT) | 5 | MTBLS8920, 12944, 6038, 7260, 12692 |
| DOC_UNRESOLVED (ROWS_EXCEED_COUNT) | 3 | MTBLS524 (51 vs 58), 13960 (114 vs 121), 13608 (60 vs 64) |
| DOC_UNRESOLVED (COUNT_EXCEEDS_ROWS) | 3 | MTBLS1584 (100 vs 96), 2774 (431 vs 330), 13787 (134 vs 133) |
| DOC_UNRESOLVED (S1_ONLY_EVIDENCE) | 3 | MTBLS12968, 12636, 12590 |

## 4. Why x = 0 is not a matter of judgement

- **The binding requirement.** STRICT needs one-sample-per-participant evidence from a publication: either Q4 = explicit_one with a publication quote, or a Q7 sample-to-participant table.
- **What was searched.** Independently of the extractors, the reviewer searched all 43 publication views for one-sample phrasing:
  - "each/every participant|patient|subject|…";
  - "per participant|patient|…";
  - "one/single … sample|specimen|swab";
  - "a … sample … from each/all";
  - "single/one time point"; "cross-sectional"; "once".
- **What it found.** Not one sentence states that each participant contributed one sample to the analysed set.
  - The closest is MTBLS13960: "Stool samples were collected from each participant…". It does not state one sample, and that study's count (114) differs from its rows (121) in any case.
  - The views contain 16 tables. 15 are participant-characteristics tables; none maps samples to participants.
- **Why adjudication cannot change it.** Adjudication only confirms or rejects STRICT and CONTRADICTED grades. A rejection gives DOC_UNRESOLVED, never STRICT. **[Certain] x = 0 under the frozen protocol.**

## 5. DR-7 consequence (computed with the frozen code; formal output follows adjudication)

| Measure | x | One-sided 95% upper bound on the stratum-M pool (hypergeometric, N = 371) | DR-7 threshold | Futility reached |
|---|---|---|---|---|
| Primary (STRICT) | 0 | **25** | 94 | yes |
| Sensitivity (STRICT + COUNT) | 5 | 88 | 94 (not used by DR-7) | would be reached |

- With x = 0, the design-based bound says that, with 95% confidence, at most 25 of the 371 stratum-M candidates are documentation-recoverable under the STRICT standard. Even the weaker COUNT standard gives at most 88.
- These are confidence bounds, not logical bounds. Per DR-7, the result ends documentation-only recovery in stratum M as a P4V supply route.
- It is not proof that confirmatory validation is impossible. It says nothing about A35 verification or a future author-contact phase, which estimate different quantities.
- **The pool overstates distinct studies.** Q6 found that MTBLS1122 shares one sample set with MTBLS1124, 1129 and 1130, all four in the 390-candidate pool (ionisation-mode deposits of one study). Sibling deposits are also named: MTBLS6038/6039 (targeted and untargeted assays on different cohorts, overlap not stated) and MTBLS2291/2349 (plasma and swab). Any future count of "independent studies" from this pool must collapse such groups.

## 6. Reviewer findings

1. **View filter leak (docrec v0.1.1 limitation).**
   - `MTBLS9586/S2_PMC11338663.view.txt` contains "Table 2 DEMs identified between patients with HF and HC", a table of differential-metabolite statistics.
   - It entered because the caption matches "patients" and uses the abbreviation "DEMs", which escapes the `metabol|differential|signific` exclusion.
   - It is the only such table among the 16 table captions in the 43 views.
   - The extractor reported it and cited nothing from it (verified). These are published summary statistics, not deposited values. The exposure should still be declared if MTBLS9586 is ever considered under P4V's provenance gate.
   - The code is frozen and is not changed for T1. A later version would add caption-abbreviation handling (e.g. `\bDE[MPG]s?\b`).
2. **PDF view heuristic misses STAR+METHODS** (known before T1). Four graded records carry `VIEW_HEURISTIC_PDF`. None is affected by the gap, since none depends on a STAR+METHODS block.
3. **Brief limitation: the census assay file name was not supplied to the extractors.**
   - Only the matrix and platform class were given. This bears on two studies where the deposit has several assay files with different sample sets:
     - MTBLS718 (Q2 not established);
     - MTBLS2774 (the deposit description lists both cohorts, 330 + 101; the chosen assay's 330 rows equal the discovery cohort, but no document says which assay file holds which cohort).
   - Neither study has one-sample evidence, so neither could reach STRICT. The reviewer did not use the row count to rewrite Q2, because that would be circular (Q5 is computed from Q3).
4. **Rule-8 count gap (docrec v0.1.1 limitation).**
   - `check_count` accepts the count integer from any quote in Q3, including an S1 quote, provided Q3 also has a publication item. In MTBLS12655 the integer 650 appears only in S1.
   - The publication states the four disjoint components (182 + 364 + 52 + 52 = 650), so the COUNT grade is supported in substance.
   - A later version should require the integer or its derivation from publication quotes.
5. **MTBLS11406 COUNT is doubtful on Q2.**
   - The deposit labels its assays LC-MS but carries proteomics protocol text. The S3 paper (same title, accession quoted) says the metabolome was measured by NMR and the MS data are on iProX.
   - The participant set (20 + 20 + 10 = 50 = rows) is the same either way, but the platform element of Q2 conflicts.
   - COUNT is sensitivity only and is not adjudicated under §10. The reviewer would not rely on it.
6. **MTBLS24 is the one study whose documentation conflicts in a way that matters** (details in the adjudication sheet). It cannot reach STRICT: no publication states one sample per patient, and the deposit description states three time points.

## 7. Adjudication (owner)

- `docrec_T1_adjudication.csv` lists the 15 CONTRADICTED records with their exact `record_sha256`. The `confirmed` column is blank, for the owner to fill with `yes` or `no`.
- The columns `reviewer_recommendation` and `reviewer_reason` are the assembler's advice. The grader ignores them.
- Recommended: confirm 12 (MTBLS10972 marked "weaker"); reject 3 (MTBLS24, 1122, 15433).
- Rejected CONTRADICTED records become DOC_UNRESOLVED. If the owner follows all 15 recommendations, the final tally is STRICT 0, COUNT 5, CONTRADICTED 12, UNRESOLVED 23.

## 8. What T1 establishes, by evidence status

- **CERTAIN:**
  - Under the frozen protocol, 0 of 40 randomly sampled stratum-M studies are documentation-recoverable, and the DR-7 futility criterion is met.
  - 11 of 40 have no publication text in the evidence set.
  - No publication view states one sample per participant.
- **SUPPORTED** (pending adjudication): at least 12 of 40 (30%) document repeated sampling or non-participant samples in the chosen assay, so independence is not merely undocumented for those studies but contradicted.
- **HYPOTHESIS, not claimed:** that the repository contains few truly independent cohorts. T1 measures what documentation establishes (E), not true independence (I). The Stage A distinction applies unchanged.
