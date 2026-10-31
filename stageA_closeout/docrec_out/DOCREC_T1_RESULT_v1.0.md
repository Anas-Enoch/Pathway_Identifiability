# Documentary recovery, tranche T1: result record v1.0

Date: 2026-10-31. Protocol `DOC_RECOVERY_PROTOCOL_v0.2.md` (frozen), `docrec.py` v0.1.1, salt option (b).
Status: **final grades computed after owner adjudication. DR-7 futility criterion reached.** Owner acceptance pending (§9).

## 1. Result

- **Sample.** T1 is a simple random sample of 40 of the 371 stratum-M candidates (NIST Beacon v2 pulse 2-2005958, 2026-10-30T15:00:00Z).
- **x = 0.** None of the 40 is `DOC_RECOVERED_STRICT`.
- **Bound.** The one-sided 95% design-based (hypergeometric) upper bound on the number of STRICT-recoverable studies in stratum M is **25 of 371**. The Clopper–Pearson figure, reported for comparison, is 26.8.
- **DR-7.** The threshold is 94, so the criterion is reached. Documentation-only recovery in stratum M stops as a P4V supply route.
- **What DR-7 does not say.** It is an operational criterion. It is not proof that confirmatory validation is impossible, and it does not bound the true number of independent cohorts in the repository. It leaves A35 verification and any author-contact phase untouched.

| Final grade (reason) | n |
|---|---|
| DOC_RECOVERED_STRICT | 0 |
| DOC_RECOVERED_COUNT (sensitivity only) | 5 |
| DOC_CONTRADICTED (adjudicated): repeated sampling 9, non-participant samples 3 | 12 |
| DOC_UNRESOLVED (adjudication rejected) | 3 |
| DOC_UNRESOLVED (Q2 not established) | 6 |
| DOC_UNRESOLVED (no assay-scoped count) | 5 |
| DOC_UNRESOLVED (rows exceed count) | 3 |
| DOC_UNRESOLVED (count exceeds rows) | 3 |
| DOC_UNRESOLVED (S1-only evidence) | 3 |
| **Total** | **40** |

The sensitivity figure (STRICT + COUNT, overlap-collapsed) is 5, giving an upper bound of 88 (< 94). DR-7 does not use it.

## 2. Robustness of the futility conclusion

- **Requirement.** Futility fails only if x ≥ 6 (upper bound 99).
- **Evidence.** No reading view in the sealed set states that each participant contributed one sample. The search is recorded in `DOCREC_T1_EXTRACTION_REVIEW_v1.0.md` §4.
- **Adjudication.** Owner adjudication can only confirm or reject STRICT/CONTRADICTED grades; a rejection gives UNRESOLVED.
- **Unread publications.** Three publications could not be obtained or were not admissible:
  - MTBLS10972 (publisher access);
  - the named S2 paper of MTBLS2291;
  - MTBLS12944 (named preprint not retrieved; a later journal version was not named by the deposit and was excluded).
- Even if all three supplied STRICT evidence, x = 3 and the upper bound would be 65.
- **[Certain] The DR-7 outcome does not depend on any unread document, any adjudication decision, or any extraction judgement.**

## 3. Chain of custody

| Object | SHA-256 |
|---|---|
| Candidates `stageA_documentary_candidates_v1.0.csv` | c1b5a4893356c517b64ed51c8777c372c52dfbc4210bc86998e7b63be18ab0d7 |
| Protocol `DOC_RECOVERY_PROTOCOL_v0.2.md` | 664d222f7f0a90b1a79ecf471e1f3ea647a10f70ace563ca8b639d92be03ba4c |
| `docrec.py` v0.1.1 | 419776ec8c378d5c275e1fcc7d046b179224d8d084362695541607079d520eb2 |
| Lock `docrec_lock.json` (signed 2026-10-28T22:43:50Z; beacon_not_before 2026-10-30T15:00:00Z) | 678933f642aafc205130b7f0d1243a86cd4191c7143ebd2f2b9179a927b859d2 |
| Order `docrec_T1_order.csv` (recomputed byte-identically from the salt) | 081f94bf337a5b4cbd9f5f5956b0908f21a20f6689c2b55b20c5104092464485 |
| Evidence universe (sealed; seal signed before extraction) | d6803c2666ebdc54f793659c4aec77b60c6c0c4f0a09c71f6b02363df074dbd2 |
| `docrec_authorisation.txt` (fetch lines + seal line) | 76b4d0bc61a4d9b885b9e49a6019c8d92fb558f0da017ddbfd1ca36e3a25fcb2 |
| Extraction brief v1 | b53bc3df4e4887e21b23cf2f35af353a83e7562aa94444c37e4bfe1fc5a19a22 |
| `docrec_T1_evidence.jsonl` | 4ff8e6177d3707fcf0728ede4f56ec31e031be2d2f2c8773d52ac11e686c2b08 |
| `docrec_T1_adjudication.csv` (owner decisions) | f2e1a15d109389788368af84f6e19cd81ca065e7db29945605c7e3fc252ab727 |
| `docrec_T1_results.csv` | 43da6ea15e99fc8b0220d7b869a737a9cf090299e57a53a7d7ca3ea3d9a76255 |
| `docrec_T1_report.json` | 60514173d63057d4d76e5c699eecf0813978bc5c001f4269c9476b50402f6cb5 |

The provisional grades and record hashes in `docrec_T1_results.csv` were checked against an independent grading in the assembler's environment. All 40 were identical.

## 4. Execution history and deviations (all disclosed; none affects x)

1. **Superseded lock.**
   - A first lock (`0025f4bf…1e4b`, signed 22:40:19Z, beacon 2026-10-29T00:00Z) was overwritten when `plan` was re-run. `plan` does not refuse an existing signed lock.
   - The 29 Oct pulse (chain 2, pulse 2003618) was retrieved but **never used**.
   - The governing lock is `678933f6…`, signed before its own beacon time; its pulse index is consistent (2003618 + 2340 = 2005958).
   - The lock file is now read-only. A future `docrec.py` version should refuse to overwrite a signed lock.
2. **Retrieval.**
   - One T1 run plus retries. Three Europe PMC full-text requests returned HTTP 500 and were finalised after 3 attempts (MTBLS718, 7260, 89).
   - The owner supplied PDFs of the named S2 publications for these and two others (MTBLS92, 24) under DR-9. All were registered and hashed before the seal.
3. **No publication text** in the evidence set for 11 studies: MTBLS10972, 12058, 12326, 12590, 12636, 12692, 12722, 12861, 12944, 12968, 8920. These could only be UNRESOLVED or CONTRADICTED (§5.4).
4. **The PDF view heuristic omits "STAR+METHODS" headings.** No T1 record depends on such a block.
5. **Reading-view leak.**
   - The MTBLS9586 view contains a differential-metabolite table ("Table 2 DEMs …"), whose caption evades the exclusion pattern.
   - Nothing from it was cited. Its exposure is to be declared if MTBLS9586 is ever considered under a provenance gate.
6. **Assay file name not given to extractors.** The brief supplied only the census matrix and platform class. This bears on MTBLS718 and MTBLS2774, neither of which has one-sample evidence.
7. **Rule-8 count gap.** `check_count` accepts the count integer from an S1 quote when Q3 also has a publication item (MTBLS12655, a COUNT grade). The publication gives the disjoint components.
8. **Extraction and adjudication.**
   - Extraction was done by five Claude sub-agents under one hashed brief, blind to row counts, with no reviewer edit after grading.
   - Owner adjudication: 15 records; 12 confirmed, 3 rejected (MTBLS24, 1122, 15433).
   - All 15 decisions coincide with the assembler's written recommendations. The adjudication file carries the owner's note "owner reviewed evidence 2026-10-31".

## 5. Findings beyond the yield (from the T1 records)

1. **The binding gap is a reporting convention, not missing papers.**
   - 29 of 40 studies had at least one publication view.
   - None states that each participant contributed one sample to the analysed set, and none publishes a sample-to-participant table.
2. **Dependence is common where documentation exists.**
   - 12 of 40 (30%) have adjudicated explicit evidence of repeated sampling or non-participant samples in the chosen assay. Repeated sampling (9): time series, pre/post treatment, follow-up visits, several cysts per patient. Non-participant samples (3): cultured cell lines.
   - All entered the candidate pool through unique or degenerate source labels. **[Certain]** Unique source labels did not indicate independence for these studies.
3. **Counts often fail to reconcile.**
   - 13 studies not finally graded CONTRADICTED have an assay-scoped participant count with a publication item. 7 equal the deposited biological rows; 6 differ (3 lower, 3 higher).
   - None of the 7 equal-count studies has one-sample evidence:
     - 5 are COUNT;
     - MTBLS24 and MTBLS15433 are UNRESOLVED after adjudication.
4. **Duplicate deposits inside the pool.**
   - MTBLS1122, 1124, 1129 and 1130 are deposits of one sample set, all in the 390-candidate pool.
   - Sibling deposits are also named: MTBLS6038/6039 and MTBLS2291/2349.
   - The pool count overstates distinct studies.

## 6. Evidence status

- **CERTAIN:**
  - x = 0 and DR-7 reached under the frozen protocol;
  - robustness as in §2;
  - 11/40 without publication text;
  - no one-sample statement in any view.
- **SUPPORTED:** 12/40 contradicted by explicit documentation (owner-adjudicated).
- **HYPOTHESIS (not claimed):** that MetaboLights holds few truly independent cohorts. T1 measures what documentation establishes (E), not true independence (I).
- **UNKNOWN:**
  - the yield of an author-contact phase;
  - the A35 verification count (fetched 2026-10-31, not yet registered, sealed or extracted).

## 7. What this result does not change

- No change to any census rule or L1 count (Stage A: L1 = 0).
- No change to P4V v0.12.
- No external outcome was opened by any step of T1.
- No measured metabolite value was read.

## 8. Files

`docrec_out/`:
- `docrec_T1_order.csv`
- `docrec_fetch_log.tsv`
- `store/` (not published: full texts)
- `docrec_T1_evidence_universe.json`
- `docrec_T1_evidence.jsonl`
- `docrec_T1_adjudication.csv`
- `docrec_T1_results.csv`
- `docrec_T1_report.json`
- `DOCREC_T1_EXTRACTION_REVIEW_v1.0.md`
- `extraction/` (brief, validator, blind universe, batch definitions, raw batch outputs)

## 9. Owner acceptance

Accept this record by appending to `docrec_authorisation.txt`:
`DOCREC T1 RESULT ACCEPTED <date> record=<SHA-256 of this file> owner=<name>`
