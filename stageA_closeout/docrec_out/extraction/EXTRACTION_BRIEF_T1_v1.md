# Documentary recovery T1: extraction brief v1 (fixed before any record was written)

This brief governs the extraction of evidence records for tranche T1 under `DOC_RECOVERY_PROTOCOL_v0.2.md` (frozen; SHA-256 664d222f…ba4c). The protocol text is at `/home/claude/nova/census_in/deliver/DOC_RECOVERY_PROTOCOL_v0.2.md`. Read its §3 (the S3 note), §4, §6, §8 and §9 before starting. Where this brief and the protocol differ, the protocol wins.

You are the **extractor**. You read documents and fill structured answers with verbatim evidence. You do not grade. A deterministic grader and then the owner decide grades. The evidence set was sealed and signed before anyone read it; nothing may be added to it.

## 1. Reading boundary (hard rules)

- Read **only** the files listed under `usable_files` for your studies, at `/mnt/user-data/uploads/github.pathway/stageA_closeout/docrec_out/store/<ACCESSION>/<file>`.
  - `S1_i_Investigation.txt`: the MetaboLights investigation file (deposit metadata).
  - `S2_*.view.txt`: reading view of a publication named by the deposit.
  - `S2M_*.view.txt`: reading view of an owner-supplied PDF of a named publication (heuristic view).
  - `S3_*.view.txt`: reading view of a later paper that mentions the accession.
- Do **not** open anything else: no other studies' folders, no `.xml`, `.json`, `.pdf`, `index.html`, manifest, order, candidate, universe or result files, and no web access. Do not try to find out how many sample rows the deposit has; you are deliberately not told.
- No measured metabolite values exist in these files; if you ever meet one, do not use it.

## 2. What to extract (protocol §6), per study, for the census-chosen assay

The chosen assay is identified by `chosen_assay_matrix` and `chosen_assay_platform_class` (MS or NMR) in your batch file, plus the deposit title. If the i_ file lists several assays with that matrix and platform class (for example positive and negative ionisation, or HILIC and reversed phase), treat them as one sample set unless a source says the sample sets differ.

**Q1 Participants** — `yes` / `no` / `not_stated`.
- `yes`: a source says the samples of this matrix come from human participants (patients, volunteers, subjects, donors of biopsies or surgical specimens, children, women…).
- `no`: a source says the selected assay's samples are cell lines, organoids, in-vitro cultures (including primary cells cultured or treated in vitro), animals, or pooled/reference material.
- Evidence is required for `yes` and `no`.

**Q2 Assay membership** — `established` / `not_established`.
- `established` only if a source passage identifies the sample set analysed in the deposited assay: its matrix, its analytical platform (any MS or NMR variant) and the population (cohort, groups or sub-study).
- Use the S1 i_ file (study description, design, protocols, assay list) to decide what was deposited, then find the passage(s) that describe that set.
- If a paper analyses several sample sets (discovery/validation, several matrices, several cohorts), Q2 is established only if you can tell which set(s) the deposit holds. If the deposit is clearly all of them, Q2 covers all of them. If you cannot tell, answer `not_established` and explain in `notes`.
- **S3 papers** (later papers mentioning the accession) establish Q2 only if they describe the deposited assay itself: the paper's own data deposited under the accession, or an explicit re-analysis of that deposit with its samples described. If every Q2 publication item would come from S3, at least one Q2 quote must contain the exact accession (for example "deposited to MetaboLights under accession MTBLS1234"), which ties the described samples to the deposit. A paper that only lists the accession among datasets it compares with does not establish Q2.
- S1 quotes may be added to Q2 as well, but a recovery needs a publication item in Q2.

**Q3 Participant count** — object with `scope`, `count`, optional `count_derivation`, `evidence`.
- `scope = "selected_assay"`: the source states the number of distinct participants whose samples make up the analysed set for this matrix and platform. Examples: "Serum samples from 58 patients were analysed by LC–MS"; "30 SLE patients and 28 healthy controls … faecal samples were collected" (a sum of disjoint groups).
- An enrolled, recruited or screened number is `selected_assay` only when the text explicitly ties that number to the analysed samples ("samples from all 60 enrolled patients were analysed"). If the text gives an enrolled number and the analysed set may differ (exclusions, subsets, unstated), use `scope = "study_enrollment"`.
- A number of *samples* is not a number of *participants*. Use it only if the same passage also states the participants (for example "58 serum samples from 58 patients").
- `scope = "none"` with no count if nothing usable is stated.
- The count must appear as an integer in one of the Q3 quotes, not followed by a unit ("100 µL", "58 years" do not count). Otherwise give `count_derivation`: a list of `{"text": <verbatim number text from a quote>, "value": <int>}` for **disjoint** groups that the source states, summing to `count`. Each `text` is purely digits (e.g. "30") or purely English number words (e.g. "twenty-eight"), never mixed. Never subtract, never use percentages, never multiply.
- Report the count for the set that the deposit holds (from Q2). If the deposit holds several disjoint described sets, sum them.

**Q4 Samples per participant (in the analysed set)** — `explicit_one` / `explicit_repeated` / `not_stated`.
- `explicit_one`: the source states that each participant contributed one sample to the analysed set ("a single fasting blood sample was collected from each participant", "one stool sample per subject", "each patient provided one sample"). A cross-sectional design or "blood samples were collected at admission" is **not** enough without a one/single/each statement.
- `explicit_repeated`: the source states more than one sample per participant **in the analysed set**: several time points or visits, pre/post-intervention samples, longitudinal sampling, paired tumour and normal tissue from the same patient analysed together, several biopsies, duplicate biological samples, or two matrices that are both inside the selected assay.
- Samples of other matrices that belong to other assays do not make the selected assay repeated.
- Analytical replicate injections of the same sample ("each sample was injected in triplicate") are not additional samples. Do not use them for `explicit_repeated`; write `ANALYTICAL_REPLICATES_STATED` in `notes`.
- Pooled QC samples are not participant samples; ignore them.
- Evidence is required for `explicit_one` and `explicit_repeated`. A STRICT grade needs the one-sample evidence from a publication, so cite the publication sentence when there is one.

**Q6 Overlap** — `flags` (list of strings), `confirmed_overlap` (list of other accessions), `evidence`.
- Flags: `cohort:<name>`, `trial:<registration id>`, `authors_period:<lead group or first author; recruitment years>`, `biobank:<name>` — whatever the sources state. Keep them short.
- `confirmed_overlap` only when a source explicitly says the participants or samples are the same as those of another MetaboLights accession (name it); evidence required.

**Q7 Mapping** — `one_to_one` / `shows_repeated` / `none`.
- Only a table or list in the sources that maps each sample to a participant counts. Participant-characteristics tables are not mappings. Expect `none` almost always.

## 3. Evidence rules (enforced by the validator)

- Each evidence item is `{"source": <file name>, "locator": <where>, "quote": <verbatim text>}`.
- `source` is the bare file name (e.g. `S2_PMC6509220.view.txt` or `S1_i_Investigation.txt`).
- `quote` must occur verbatim in that file (whitespace differences are ignored; case and punctuation are not). Copy-paste; do not retype, abridge or join non-adjacent text with "…".
- At least 3 tokens; must not contain view markers: `[TITLE]`, `[ABSTRACT]`, `[SECTION:`, `[TABLE:`, `[FRONT]`, `[METHODS_PDF]`, `[VIEW_HEURISTIC_PDF]`.
- `locator`: e.g. `S2_PMC6509220.view.txt › SECTION: Materials and methods › Subjects`, or `ABSTRACT`, or `S1 › Study Description`.
- For `not_stated` / `not_established` / `none` answers leave `evidence` empty and explain in `notes`.

## 4. Record format (one JSON object per line, UTF-8)

```json
{"accession": "MTBLS0000", "tranche": "T1",
 "extractor": "Claude assembler, sub-agent batch K, brief v1",
 "extracted_utc": "2026-10-31T17:30:00Z",
 "sources_consulted": ["S1_i_Investigation.txt", "S2_PMC0000000.view.txt"],
 "Q1": {"answer": "yes", "evidence": [{"source": "...", "locator": "...", "quote": "..."}]},
 "Q2": {"answer": "established", "evidence": [ ... ]},
 "Q3": {"scope": "selected_assay", "count": 58,
        "count_derivation": [{"text": "30", "value": 30}, {"text": "28", "value": 28}],
        "evidence": [ ... ]},
 "Q4": {"answer": "not_stated", "evidence": []},
 "Q6": {"flags": [], "confirmed_overlap": [], "evidence": []},
 "Q7": {"answer": "none", "evidence": []},
 "notes": "What the deposit holds per S1; which described set matches; ambiguities; ANALYTICAL_REPLICATES_STATED if applicable."}
```

- Omit `count_derivation` when the count appears directly in a quote. With `scope = "none"`, omit `count`.
- `extracted_utc` is the real UTC time (`date -u +%Y-%m-%dT%H:%M:%SZ`).
- `sources_consulted` lists every file you actually read for that study.

## 5. Conduct

- Decide each answer from the text alone. When the text falls short of a definition, choose `not_stated` / `not_established` / `none` / `study_enrollment` and say why in `notes`. But never omit explicit evidence of repeated sampling or of non-participant samples: that must be recorded.
- Do not guess, infer from typical practice, or use outside knowledge of the study.
- Studies with only `S1_i_Investigation.txt` are still extracted from S1 alone.
- Run the validator until every record is VALID:
  `python3 -I /tmp/claude-0/-home-claude/4597b89f-b4a3-551d-b8bb-cf1d52b09264/scratchpad/extract/validate.py <your jsonl> <your batch json>`
  It checks schema, verbatim quotes and count forms. It does not grade, and VALID says nothing about the grade.
