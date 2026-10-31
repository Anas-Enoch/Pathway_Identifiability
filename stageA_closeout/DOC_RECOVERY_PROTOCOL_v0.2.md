# Documentary recovery (§3.7): protocol v0.2 (candidate for freezing)

**Status: DRAFT v0.2, not frozen.** 2026-10-28.

- v0.2 incorporates the owner's A36 decision memo: A36-3 approved; DR-1 to DR-6 approved; DR-7 approved with its statistical qualification; two corrections to the grading rules.
- It becomes executable only when all of the following hold:
  - DR-1b, DR-8 and DR-9 (§0) are settled;
  - this file, `docrec.py` and the candidate file are hashed into a signed lock;
  - an authorisation line exists for the tranche.
- It changes no census rule and no L1 count, and authorises no change to P4V v0.12.
- Implementation: `docrec.py` v0.1.1 (stdlib only; `python3 -I docrec.py selftest`, 146 checks, passing on Python 3.9 and 3.13).
- An independent pre-freeze review ran in two rounds and found 13 defects, then 7 introduced by the fixes. All 20 are fixed, and each has a regression check.

## 0. Decisions incorporated, and three still open

| Item | Content | Where |
|---|---|---|
| DR-1 | T1 = 40 studies from stratum M, never replaced | §2 |
| DR-2 | S1–S3; S3 query, provider, order, cap and de-duplication frozen; an S3 hit never establishes Q2 by itself | §3 |
| DR-3 | no author contact in T1; T1 estimates recoverability *from documentation*; a later contact phase estimates a different quantity and is reported separately | §11 |
| DR-4 | V1–V4 for all 19 A35 studies, separately; report N_A35,verified and N_A35,recovered separately | §7 |
| DR-5 | extractor → deterministic grading → owner adjudication; STRICT is provisional until the evidence and exact count are checked against the chosen assay | §9–10 |
| DR-6 | complete separation from P4V; recovery establishes participant independence only | §13 |
| DR-7 | futility x ≤ 5 as an operational criterion for documentation-only recovery in stratum M, not proof that confirmatory validation is impossible | §11 |
| C1 | the participant count must be documented for the selected assay and equal its biological rows; an enrollment count never suffices | §6, §9 |
| C2 | DOC_CONTRADICTED only on explicit evidence of dependence; absence of a statement is not contradiction | §9 |
| **DR-1b (open)** | how the T1 sample is drawn: **(a)** fixed salt → a *systematic hash-distributed feasibility sample*; bounds are model-based (exchangeability). **(b)** salt = a NIST Randomness Beacon v2 pulse at a time fixed in the signed lock, i.e. after freezing → a *simple random sample* of 40 of 371; exact hypergeometric bounds are design-based. **Recommended: (b).** It costs one beacon lookup, no one can choose the salt, and the futility threshold is x ≤ 5 under both (§11). | §2, §11 |
| **DR-8 (open)** | is the abstract part of the reading view? **Recommended: yes.** It is public and aggregate, often the only place the count appears in short papers, and holds no held-out data. | §4 |
| **DR-9 (open)** | may the owner supply full texts that are not open access (institutional access), registered and hashed before the document freeze? **Recommended: yes, S2 publications only.** Otherwise T1 measures open-access availability as much as documentation. | §3, §5 |

## 1. Population and strata

The input is `stageA_documentary_candidates_v1.0.csv` (390 rows, SHA-256 `c1b5a4893356c517b64ed51c8777c372c52dfbc4210bc86998e7b63be18ab0d7`). `docrec.py` refuses any other digest.

- **Stratum M (371):** `M5_SOURCE_UNIQUE` (312) + `M5_SOURCE_DEGENERATE` (59). `n_rows` is the census `n_biological_rows`.
- **Stratum A (19):** `A35_ALLQC_DIAGNOSTIC`. `n_rows` is all sample rows. The biological count is unknown until V3.

## 2. Selection order and sampling statement (DR-1, DR-1b)

- key(acc) = SHA-256(`salt + "|" + acc`) as lowercase hex; ascending sort within each stratum.
- T1 = the first 40 of stratum M.
- Salt:
  - option (a): the literal `stageA-docrec-v1`;
  - option (b): `NISTv2:<chainIndex>-<pulseIndex>:<outputValue>`, taken from the pulse returned by the URL that `plan` writes into the signed lock. That URL is `https://beacon.nist.gov/beacon/2.0/pulse/time/<beacon_not_before in Unix ms>`, which returns the pulse at that time or the next one. `beacon_not_before` is chosen later than the signing time. The owner opens the URL after that time, and anyone can re-open it to check the value.
- `docrec.py fetch` writes `docrec_<tranche>_order.csv` (stratum, rank, accession, origin, n_rows, key, in_T1), one file per tranche, and records its SHA-256 in the fetch-log header before the first request. `docrec.py order --salt …` is a network-free dry run.
- No study is skipped or replaced. Missing evidence gives `DOC_UNRESOLVED`.
- Wording of the result:
  - under (a): "a systematic, hash-distributed feasibility sample", with no claim of design-based unbiasedness;
  - under (b): "a simple random sample of 40 of the 371 stratum-M candidates".

## 3. Sources and the frozen retrieval procedure (DR-2)

- `docrec.py fetch` is run by the owner on his Mac. Neither the assembler's environment nor the linked sandbox can reach Europe PMC.
- Every request is checked against an allow-list and journaled (UTC, URL, HTTP status, bytes, `Content-Length`, elapsed time) with append + flush.
- Redirects are never followed.
- A body whose length differs from `Content-Length` is refused and retried once (the D2 lesson). The response size limit is 20 MB.
- **Outcome classes:**
  - *definite* (a property of the resource): 200; a redirect (never followed); 400; 404; 410; oversize;
  - *indeterminate*: everything else, including transport failure, 5xx, 429, 401/403/407/408 (which could be a firewall, proxy or timeout), a short body, and JSON or XML that does not parse (a possible truncation without `Content-Length`).
- A study is `complete` only if every request had a definite outcome. Otherwise the next `fetch` repeats it.
- After 3 attempts a study is finalised as it stands, with the note `FINALISED_AFTER_3_ATTEMPTS`. Its missing sources count as unavailable.
- 5 consecutive indeterminate outcomes on one host stop the run (exit 2). A run that ends this way does not count as an attempt for the studies it touched, so a long outage cannot finalise studies without their sources.
- `fetch` exits 3 while any member is incomplete.
- A refused URL halts the run (exit 1).

| Code | What | Exact request |
|---|---|---|
| S1 | ISA metadata | the study's directory index, then every `i_*.txt` it lists (normally one); stratum A also the `s_*.txt` and `a_*.txt` named in the investigation file. Only under `https://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/MTBLSn/`. `m_*`, `FILES/` and anything else are refused. |
| S2 | publications named in S1 | each `Study PubMed ID` or `Investigation PubMed ID` (digits) → `search?query=EXT_ID:<pmid> AND SRC:MED`; each `Study Publication DOI` or `Investigation Publication DOI` (lower-cased; `https://doi.org/` stripped; must be `10.<4–9 digits>/…`) → `search?query=DOI:"<doi>"`. Then for each record with a PMCID: `PMC<n>/fullTextXML`. |
| S3 | other publications citing the accession | `search?query="MTBLSn"` with `synonym=FALSE`, `resultType=core`, `format=json`, `pageSize=1000`, `cursorMark=*` (pages followed to at most 1000 hits). Sorted by (`firstPublicationDate`, source, id) ascending. A record that shares any identifier (DOI, PMID, PMCID or source:id) with an S2 record or with an earlier record is dropped. **The first 5** remaining are kept, regardless of text availability. Full text as for S2. An S3 article is usable only if its text contains the exact accession with non-alphanumeric boundaries. **An S3 article never establishes Q2 on its own.** It must describe the deposited assay itself. |

- The base URL is `https://www.ebi.ac.uk/europepmc/webservices/rest/`.
- Search responses and full texts are stored with their SHA-256.
- Owner-supplied texts (DR-9) are placed in `<store>/MTBLSn/manual/` as PDF. A text rendering (`pdftotext`) is stored beside each PDF, and both are registered with `docrec.py register` before the document freeze (§5). `register` runs only once every tranche member is complete, because a later fetch would rewrite the manifest.

**Forbidden:**
- any measured metabolite value;
- MAF / `m_` files, `FILES/`, raw or derived data;
- supplementary data files;
- reading any part of a document outside its reading view (§4).

## 4. Value-blind reading view

The extractor reads only `view.txt`. `docrec.py` builds it deterministically.

- **From JATS XML:**
  - the article title;
  - the abstract (if DR-8 = yes). In a structured abstract, the parts titled or labelled `Results`, `Findings`, `Conclusion(s)`, `Interpretation` or `Discussion` are dropped. This covers `sec` parts, paragraphs that open with a bold or italic label, and inline "Results:" labels; the text is cut at the label;
  - every **top-level** `sec`/`notes` (directly under `body` or `back`, looking through untitled, untyped wrappers) whose title contains the word `method(s)`/`methodology` (e.g. "Research Design and Methods", "STAR★Methods") or *begins* with a methods-type heading. The heading list is `materials`, `methods`, `methodology`, `experimental`, `research design`, `patients`, `participants`, `subjects`, `volunteers`, `study design/population/participants/subjects/cohort/setting`, `cohort`, `population`, `recruitment`, `enrolment`, `sample collection/preparation/acquisition`, `specimens`, `ethic…`, `data availability`, `availability of data`, optionally numbered or prefixed `online`. A top-level section also enters if its `sec-type` matches `method|subject|data-availability`. All descendants of an included section are included.
  - A top-level section whose title or type matches `result|discussion|conclusion|supplementar|finding|outcome` is skipped with all its descendants.
  - Inside an included section, only `result|discussion|conclusion|supplementar` exclude a subsection, so methods subsections such as "Outcome measures" stay.
  - Only paragraph-level text is taken. Figures, tables, supplementary material, boxed text, media and display formulas nested *inside* a paragraph are removed;
  - every `table-wrap` whose label or caption matches `characteristic|demographic|baseline|participant|patient|subject|clinical` and does not match `metabol|lipid|feature|compound|fold|differential|signific|biomarker|pathway|concentration|abundance|level|biochem|laborator`.
- **From an owner-supplied PDF text:**
  - the lines before the first `Introduction`/`Background` heading (title and abstract);
  - plus the lines from the first methods heading (`Materials and Methods` and its variants, `Research Design and Methods`, `(Online) Methods`, `Experimental …`, `Study design`, `Study population`, `Study participants`; a bare `Participants` line is not a heading, because tables use it) up to the next end heading: `Results`, `Discussion`, `Conclusion(s)`, `References`, `Acknowledg(e)ments`, `Funding`, `Author contributions`, `Supplementary …`, `Figure legends`, `Figures`, `Tables`.
  - A heading is a whole line of at most 60 characters, optionally numbered, in title or upper case, that *is* one of those phrases.
  - This is a heuristic; graded records citing such a view carry the flag `VIEW_HEURISTIC_PDF`.
- **S1 metadata** is read in full, since it holds no measured values.

The full documents are stored for audit only.

## 5. Document freeze

1. After `fetch`, plus `register` if DR-9 = yes, `docrec.py seal` writes `docrec_<tranche>_evidence_universe.json`. It lists every document and view, with SHA-256.
   - Seal refuses unless `docrec_<tranche>_order.csv` matches the order hash in the last fetch run for the tranche.
   - The universe must contain **exactly** the tranche members (40 for T1), every one `complete`.
   - A study whose index definitely failed (for example a 404) is sealed with its notes and graded on what exists.
2. The owner signs the seal.
3. Extraction starts only after that. Once a seal line is signed for a tranche, `fetch`, `register` and `seal` refuse that tranche, so no document can be added afterwards.
4. A study with no usable S2/S3 view is still extracted from S1 alone. It can only be `DOC_UNRESOLVED` (`Q2_NOT_ESTABLISHED` or `S1_ONLY_EVIDENCE`) or contradicted.

## 6. Questions (for the census-chosen assay of each study)

- **Q1 Participants.** The biological samples come from human participants, not cell lines, organoids, pooled or reference material, or animals: yes / no / not_stated.
- **Q2 Assay membership.** The source identifies the sample set analysed in the deposited assay (matrix, platform, cohort or sub-study): established / not_established.
- **Q3 Participant count (C1).** An integer `count` and its `scope`:
  - `selected_assay`: the source states the number of distinct participants whose samples make up that analysed set;
  - `study_enrollment`: an enrolled or recruited number not tied to the analysed set;
  - `none`.
  - The count must be stated as an integer, or be a sum of stated integers for groups the source says are disjoint.
- **Q4 Samples per participant.** `explicit_one` (the source states that each participant contributed one sample to the analysed set) / `explicit_repeated` (the source states more than one visit, time point, pre/post sample, replicate or matrix per participant in that set) / `not_stated`.
- **Q5 Correspondence** is computed, never extracted: Q3.count = `n_rows` (stratum M) or the V3 count (stratum A).
- **Q6 Overlap.** Flags (same title, cohort name, trial registration, or author group and period), plus `confirmed_overlap` only on an explicit statement.
- **Q7 Mapping.** `one_to_one` / `shows_repeated` / `none` (a published or deposited sample-to-participant table).

Every non-null answer needs at least one evidence item `{source, locator, quote}`.
- A quote has at least 3 tokens and contains no view marker (`[SECTION:` and so on).
- It must occur verbatim in a usable view, or in an S1 file, of the same study (whitespace normalised; case-sensitive). The full XML/PDF files cannot be cited.
- The Q3 count (and the V3 count) must appear as an integer in one of its quotes, not followed by a unit ("100 µL" does not count). Otherwise a `count_derivation` list of `{text, value}` items must be given:
  - each `text` is verbatim in a quote, and is either purely numeric (equal to its value) or purely English number words (parsed, and equal to its value; flagged `NUMBER_WORDS`);
  - a text mixing digits and words is rejected;
  - the values must sum to the count.
- `docrec.py` checks all of this. A record with any failure is `INVALID`, is returned to the extractor and is not graded.

## 7. Stratum A: V1–V4 (DR-4)

| Check | Pass condition |
|---|---|
| V1 organism | every row of the chosen assay is from a human participant |
| V2 assay membership | the chosen assay's rows are study samples, not blanks or QC, shown by S1 fields other than the `blank` label or by S2 |
| V3 biological count | an exact integer number of biological samples |
| V4 sampling structure | no evidence of repeated sampling; explicit repetition fails |

- Each check is `pass` / `fail` / `not_established`, with evidence.
- Status: `VERIFIED` if all four pass; `FAILED` if any fails; otherwise `UNVERIFIABLE`.
- Sources are S1 (including the `s_`/`a_` files) and S2. S3 is not searched for stratum A.
- Only VERIFIED studies go on to Q1–Q7, in their own later tranche.
- The report gives **N_A35,verified** and, later, **N_A35,documentary recovered** separately. No A35 study is treated as an independent cohort because of its origin.

## 8. Evidence record (one JSON object per line)

The record holds:
- `accession`, `tranche` (must equal the tranche), `extractor`, `extracted_utc` (`YYYY-MM-DDTHH:MM:SSZ`), `sources_consulted` (list of file names);
- `Q1` to `Q4`, `Q6` and `Q7` as above, each an object with `answer` (Q3: `count`, `scope`, optional `count_derivation`) and an `evidence` list. Q6 has `flags` (list), `confirmed_overlap` (a list of *other* accessions) and `evidence`. Stratum A records hold `V1` to `V4` instead (V3 with `count`);
- `notes`.

The schema is enforced by `docrec.py grade`.

## 9. Grading function (deterministic; rules applied in this order)

| Order | Condition | Grade / reason |
|---|---|---|
| 0 | any failed quote, citation of an unusable source, or schema error | `INVALID` (not graded) |
| 1 | Q2 ≠ established | `DOC_UNRESOLVED (Q2_NOT_ESTABLISHED)`: nothing about "the selected assay" is decided without Q2 |
| 2 | Q2 established and (Q4 = explicit_repeated or Q7 = shows_repeated) | `DOC_CONTRADICTED (REPEATED_SAMPLING)` |
| 3 | Q2 established and Q1 = no | `DOC_CONTRADICTED (NOT_PARTICIPANT_SAMPLES)` |
| 4 | Q1 ≠ yes | `DOC_UNRESOLVED (Q1_NOT_STATED)` |
| 5 | Q3.scope ≠ selected_assay | `DOC_UNRESOLVED (NO_ASSAY_SCOPED_COUNT)`: covers the owner's "120 enrolled, 100 samples" case |
| 6 | Q3.count > n_rows | `DOC_UNRESOLVED (COUNT_EXCEEDS_ROWS)` |
| 7 | Q3.count < n_rows | `DOC_UNRESOLVED (ROWS_EXCEED_COUNT)`: not contradicted, because unlabelled QC rows could explain it (C2) |
| 8 | Q2 or Q3 has no evidence item from a publication (S2, S3 or owner-supplied) | `DOC_UNRESOLVED (S1_ONLY_EVIDENCE)` |
| 8b | all publication items of Q2 are from S3 and none of their quotes contains the exact accession | `DOC_UNRESOLVED (Q2_S3_WITHOUT_ACCESSION)` |
| 9 | (Q4 = explicit_one with a publication item **in Q4**) or (Q7 = one_to_one with a publication item **in Q7**) | `DOC_RECOVERED_STRICT (PROVISIONAL)` |
| 10 | otherwise (Q4 = not_stated and Q7 = none, or one-sample evidence only from S1) | `DOC_RECOVERED_COUNT` (sensitivity only) |

**Independence of evidence (rule 8; §3.7 "obtain independent evidence; never weaken D1"; A29 rejected D1 relaxation).**
- S1 is the deposit that D1 already judged. It may identify the assay in Q2 and may supply contradiction evidence.
- A recovery grade needs the participant count (Q3) and its assay link (Q2) from a publication. A STRICT grade also needs its one-sample evidence (Q4 or Q7) from a publication.
- Otherwise a sentence in the study description, which R2 already tested, could recover the study under weaker rules.
- **S3 restriction (DR-2, rule 8b):** if every publication item of Q2 comes from S3 sources, at least one of those quotes must contain the exact accession. A paper that only mentions the accession elsewhere cannot tie its participants to the deposited assay. Rule 8b comes after the contradiction rules, so explicit repetition is still reported as contradiction.

Q6 flags are reported beside the grade. Studies with confirmed overlap count once in every yield figure.

## 10. Adjudication (DR-5)

- The owner checks every STRICT (PROVISIONAL) and every CONTRADICTED record against the filed evidence, the chosen assay and the exact count.
- Decisions are written to `docrec_<tranche>_adjudication.csv` (UTF-8, with or without BOM) with the columns `accession, provisional_grade, record_sha256, confirmed, note`. The provisional grade and record hash are copied from `docrec_<tranche>_results.csv`.
- A row applies only while both still match. A re-extracted record needs a new decision, and stale rows are listed in the report.
- Only confirmed STRICT records count as `DOC_RECOVERED_STRICT`.
- A rejected STRICT becomes `DOC_UNRESOLVED (ADJUDICATION_REJECTED)`, and a rejected CONTRADICTED becomes `DOC_UNRESOLVED`.

## 11. Yield, bounds and futility (DR-3, DR-7)

x is the number of adjudicated `DOC_RECOVERED_STRICT` in T1, counting confirmed overlaps once.

| x / 40 | CP upper × 371 (one-sided 95%) | hypergeometric upper K (N = 371) | CP lower × 371 | hypergeometric lower K |
|---|---|---|---|---|
| 5 | 90.9 | 88 | 18.8 | 20 |
| 6 | 101.9 | 99 | 25.0 | 27 |
| 15 | 191.9 | 189 | 91.7 | 94 |
| 16 | 201.2 | 198 | 99.9 | 103 |
| 19 | 228.2 | 225 | 125.3 | 128 |
| 25 | 279.3 | 277 | 179.1 | 182 |

- **DR-7 rule (prespecified, operational):** if the upper bound on the stratum-M pool is below 94, documentation-only recovery in stratum M stops scaling. Under either bound this happens exactly at **x ≤ 5**.
  - 94 is the 0.95-probability screening pool of P4V v0.12 §5c at w = 0.80, f = 0.05, q = 0.40, and so inherits those assumptions.
  - A confidence bound is not a mathematical bound on the true number.
  - Reaching the rule does not show that confirmatory validation is impossible.
  - It does not affect A35 verification or any future author-contact programme, which are different procedures with different yields.
- If x ≥ 6, the owner decides whether to run a further tranche.
- The reported bounds are those of DR-1b: hypergeometric under (b); CP, labelled model-based, under (a).
- The yield is "recoverable from documentation". It is never pooled with a later author-response yield.
- The report also gives STRICT + COUNT, overlap-collapsed, as a sensitivity figure. It never enters DR-7.
- Futility is not computed while any record is missing, invalid or awaiting adjudication.

## 12. Outputs and authorisation sequence

1. `plan` writes `docrec_lock.json`: the SHA-256 of the protocol, `docrec.py` and the candidates file, the salt option, DR-8, DR-9, and (under (b)) `beacon_not_before`, which must be in the future, and `beacon_url`. The owner signs it.
2. Under (b), after `beacon_not_before`, the owner opens `beacon_url` from the lock and copies `chainIndex`, `pulseIndex` and `outputValue` into the salt. The authorisation line is `DOCREC FETCH AUTHORISED <date> tranche=T1 candidates=<sha> protocol=<file> lock=<SHA-256 of the signed lock> salt=<salt> owner=<name>`. The last line for a tranche governs.
3. `fetch` → (`register`) → `seal`, which the owner signs with `DOCREC SEAL SIGNED <date> tranche=<T> universe=<sha> owner=<name>`.
4. Extraction (assembler) → `grade` (writes results and report; futility is computed only when no record is missing, invalid or awaiting adjudication) → owner adjudication → `grade` again.
5. Outputs: `docrec_<tranche>_order.csv`, `docrec_fetch_log.tsv`, `store/`, `docrec_<tranche>_evidence_universe.json`, `docrec_<tranche>_evidence.jsonl`, `docrec_<tranche>_results.csv`, `docrec_<tranche>_report.json`, and `SHA256SUMS` made by the owner with `shasum` at closure.
6. Only public identifiers are recorded (accessions, DOIs, PMIDs, PMCIDs), and no participant-level information beyond counts.

## 13. What recovery does not establish (DR-6)

A `DOC_RECOVERED_STRICT` study has documented participant independence in its chosen assay. It does **not** have any of the following:
- untouched evaluation provenance;
- enough RefMet-named features;
- compatible metabolite identities;
- informativeness for the acquisition task;
- study-level independence from other recovered studies beyond Q6.

P4V v0.13 must define those gates before any recovered study is admitted to a confirmatory analysis.
