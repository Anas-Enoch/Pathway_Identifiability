# Documentary recovery (§3.7): protocol draft v0.1

**Status: DRAFT for owner audit. Not frozen, not hashed, not executable.** 2026-10-28.

- This draft implements §3.7 of `EXTERNAL_DATA_FEASIBILITY_PROTOCOL_v0.7.8.md`: A31-D, extended by A35.
- Nothing in it may be applied to any study until three things exist: the owner's decisions DR-1 to DR-6, a hashed frozen version, and an authorisation line.
- It changes no census rule and no L1 count. Recovered studies are reported only as `L1_documentary_recovered`.
- Whether they count as B3 route-2 evidence for P4V is a separate P4V v0.13 decision (A36).

## 1. Population and strata

The input is `stageA_documentary_candidates_v1.0.csv` (SHA-256 `c1b5a4893356c517b64ed51c8777c372c52dfbc4210bc86998e7b63be18ab0d7`), with 390 rows in two strata:
- **Stratum M (371):** origin `M5_SOURCE_UNIQUE` (312) or `M5_SOURCE_DEGENERATE` (59). The tranche is drawn from this stratum.
- **Stratum A (19):** origin `A35_ALLQC_DIAGNOSTIC`. These studies first go through pre-verification V1–V4 (§5). Only those that pass are eligible for recovery. They are never pooled into stratum M's yield estimate.

## 2. Selection order (frozen before any evidence is examined)

- rank(acc) = SHA-256(`"stageA-docrec-v1|" + acc`), ascending, computed separately within each stratum.
- The rank is a function of the accession alone: not of study size, title, platform or likely documentation quality. This keeps the tranche yield an unbiased estimate for stratum M.
- The ranked list is written and hashed before any evidence is examined.
- **Tranche T1** is the first 40 of stratum M in rank order (DR-1: size within the owner's 30–50).
- No study is skipped or replaced. A study whose sources cannot be reached is graded `DOC_UNRESOLVED`, never swapped for the next one.

## 3. Admissible sources and access boundary

| Code | Source | T1 |
|---|---|---|
| S1 | The study's MetaboLights `i_`/`s_`/`a_` metadata already recorded by the census (checkpoint fields), with no new retrieval | yes |
| S2 | Publication(s) named in the investigation file (Study Publication DOI / PubMed ID): main text, methods, figure legends, and *descriptive* supplementary tables (participant or sample characteristics) | yes |
| S3 | Publications found by a frozen Europe PMC full-text query for the exact accession string | DR-2 |
| S4 | Author correspondence using a frozen question template; every reply filed verbatim | DR-3 (recommended: no in T1) |

The following are **forbidden**:
- any measured metabolite value;
- MAF / `m_` files, `FILES/`, raw or derived data;
- any supplementary file containing metabolite measurements. A mixed supplement is not opened.

Summary clinical statistics (for example the mean age in a Table 1) are not metabolite values and are admissible.

## 4. Questions (frozen; answered per study for the census-chosen assay `assay_file`)

- **Q1 Participants.** Are the biological samples of the chosen assay from living or deceased human participants? They must not be cell lines, organoids, pooled material, reference material or animals.
- **Q2 Assay membership.** Which cohort or sub-study described in the source corresponds to the chosen assay (matrix, platform, sample set)?
- **Q3 Participant count.** What is the exact integer number of distinct participants contributing to that assay?
- **Q4 Samples per participant.** Does each participant contribute exactly one biological sample to that assay? Any of the following answers no:
  - more than one visit, time point or pre/post sample;
  - repeated collection, or a stored unit sampled more than once;
  - a technical replicate entered as a separate biological row;
  - more than one tissue or matrix.
- **Q5 Correspondence.** Does Q3 equal, as an exact integer, the census's `n_biological_rows` for the chosen assay (stratum M), or the V3-verified count (stratum A)?
- **Q6 Overlap.** Do the participants overlap with those of another candidate? Indicators are the same cohort name, the same trial registration, the same title, or the same author group and recruitment period. This is a flag only. CONFIRMED overlap needs an explicit statement.
- **Q7 Mapping.** Is there a sample-to-participant mapping, either a published table or author-confirmed?

## 5. Stratum A pre-verification (sources S1–S2 only; before Q1–Q7)

| Check | Pass condition |
|---|---|
| V1 organism | every biological row human; mixed-organism values resolved by S1 assay/sample structure or S2 text, otherwise fail |
| V2 assay membership | the chosen assay's rows are identified as study samples, not blanks/QC, by S2 or by S1 fields other than the `blank` label |
| V3 biological count | an exact integer number of biological samples established by S1/S2 |
| V4 sampling structure | no evidence of repeated sampling (same criteria as Q4) |

A failure is recorded with its reason. Passing V1–V4 makes a study eligible for Q1–Q7; it is not a recovery.

## 6. Evidence rules and grades (the §3.7 grades, made mechanical)

- Every answer cites a verbatim sentence or table cell with its locator (source code, document, section / page / table, date accessed). Paraphrase is not evidence.
- A count must be an integer stated in the text, or a sum of stated integers for groups the text says are disjoint. Counts are never derived from percentages.

| Grade | Condition |
|---|---|
| `DOC_RECOVERED_STRICT` | Q1 yes; Q2 established; Q3 exact; Q5 equal; and Q4 = one sample per participant stated explicitly **or** shown by a Q7 mapping (§3.7 evidence a or b) |
| `DOC_RECOVERED_COUNT` | Q1, Q2, Q3, Q5 as above; Q4 not stated, with no contradicting evidence (§3.7 evidence c). Sensitivity only. |
| `DOC_CONTRADICTED` | evidence of more than one sample per participant in the chosen assay; or Q3 < biological rows; or Q1 no |
| `DOC_UNRESOLVED` | anything else, including no reachable source, Q2 not establishable, or Q3 > rows (subset deposited, membership unknown) |

- Q6 flags are reported beside the grade.
- Two studies with CONFIRMED overlap count once in any yield or pool figure.
- Grades are assigned by a frozen grading function from the structured answers, not by judgement. Its test suite must exist before freezing (§3.7 requirement).
- The owner adjudicates every STRICT and CONTRADICTED grade against the filed evidence before it is reported.

## 7. Yield report and futility rule (decided now; it does not alter §4–§6 afterwards)

Let x be the number of `DOC_RECOVERED_STRICT` among the 40 in T1. Report:
- x;
- the one-sided 95% Clopper–Pearson bounds;
- the projected stratum-M pool 371·L and 371·U.

These are *upper bounds* on P4V-eligible studies, because untouched provenance, RefMet naming, E1–E8 and informativeness still deflate the pool.

| x / 40 | U (one-sided 95%) | 371·U | 371·L |
|---|---|---|---|
| 5 | 0.245 | 90.9 | 19 |
| 6 | 0.275 | 101.9 | 25 |
| 16 | 0.542 | — | 99.9 |
| 19 | 0.615 | — | 125.3 |
| 25 | 0.753 | — | 179.1 |

- **Futility:** if 371·U < 94, documentary recovery stops as a P4V supply route, which happens at **x ≤ 5**. 94 is P4V v0.12 §5c's 0.95-probability screening pool at the most favourable tabulated scenario (w = 0.80, f = 0.05, q = 0.40). Below it, no tabulated scenario is reachable from this stratum even at the upper bound.
- **Not futile:** the owner decides on a further tranche. For reference, 371·L reaches 94 / 121 / 176 (the w = 0.80 / 0.75 / 0.70 pools) at x ≥ 16 / 19 / 25.
- The yield never changes the questions, standards or order.

## 8. Recording and outputs

- One JSON record per study: accession, stratum, rank, sources consulted (with access dates), Q1–Q7 answers with verbatim evidence and locators, grade, Q6 flags.
- `docrec_T1_results.csv`, `docrec_T1_evidence.jsonl`, SHA256SUMS.
- A result record with the §7 table.
- Identifiers kept are public accessions, DOIs and PMIDs. No participant-level information is recorded beyond counts.

## 9. Owner decisions before freezing

| Decision | Question | Recommended default |
|---|---|---|
| DR-1 | T1 size | 40 (stratum M only) |
| DR-2 | S3 (citing-literature search) in T1 | yes, query string frozen: `"MTBLSnnnn"` in Europe PMC full text, first 5 hits by date |
| DR-3 | S4 (author contact) in T1 | no; T1 measures what documents alone recover, and contact is a later tranche with a frozen template |
| DR-4 | Stratum A | V1–V4 on all 19 in parallel with T1; passers form their own tranche |
| DR-5 | Who executes | the assembler extracts evidence and fills the structured answers; the grading function grades; the owner adjudicates STRICT / CONTRADICTED |
| DR-6 | Coupling to P4V | none in this protocol; whether a STRICT grade satisfies B3 route 2 is decided in P4V v0.13 (A36) |
