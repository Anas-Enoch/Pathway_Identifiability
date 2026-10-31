# Documentary recovery, stratum A (A35): V1–V4 extraction brief v1 (fixed before any record was written)

This brief governs pre-verification of the 19 A35 studies under `DOC_RECOVERY_PROTOCOL_v0.2.md` §7 (frozen). The protocol is at `/home/claude/nova/census_in/deliver/DOC_RECOVERY_PROTOCOL_v0.2.md`. Read §3, §4, §6 (evidence rules) and §7 before starting. Where this brief and the protocol differ, the protocol wins.

You are the **extractor**. You fill structured answers with verbatim evidence. You do not grade. The evidence set was sealed and signed by the owner before anyone read it.

## 1. Why these studies are here

- The Stage A census found these 19 human deposits had no biological rows left after its QC-removal rule (rule M2, "no biological rows after QC removal").
- The census therefore did not select an assay, and its candidate file has no matrix or platform for them.
- V1–V4 ask whether the deposit's own metadata, or its named publication, shows that the deposit holds human participant samples, how many, and whether each participant appears once.

## 2. The unit (fixed here because the census gives no assay)

- V1–V4 refer to **the biological samples in the study's sample file** (`S1_s_*.txt`) that the study's assays measure.
- If different assays measure different subsets of samples, or different matrices, say so in `notes`, then answer for the set the assays share. If there is no shared set, answer for the largest set.
- QC, pooled-QC, blank, standard and reference rows are not biological samples.

## 3. Reading boundary (hard rules)

- Read **only** the files listed under `usable_files` for your studies, at `/mnt/user-data/uploads/github.pathway/stageA_closeout/docrec_out/store/<ACCESSION>/<file>`.
  - `S1_i_Investigation.txt` — the investigation file;
  - `S1_s_*.txt` — the sample file;
  - `S1_a_*.txt` — the assay files;
  - `S2_*.view.txt` — the reading view of a named publication, where one exists.
- Sample and assay files can be large. Tabulate them with short Python or awk scripts (column headers, distinct values per column, counts) instead of printing them whole.
- Do **not** open any other file or folder, and no web.
- Do not try to learn the census row counts; you are deliberately not told them.
- The assay files name raw data files but hold no measured values. Never open data files; none are present anyway.

## 4. Questions (protocol §7)

Each answer is `pass` / `fail` / `not_established`, with evidence. `pass` and `fail` need at least one evidence item.

**V1 organism (human participants).**
- `pass`: the sources state that the biological samples come from human participants (patients, volunteers, subjects, donors, children, women…), and no biological sample is shown to be non-human or cultured. An organism field reading "Homo sapiens" is not enough on its own, because cell lines are also Homo sapiens; it may be cited alongside a participant statement.
- `fail`: the biological samples are cell lines, organoids, in-vitro cultures, animals, or pooled or reference material.
- `not_established`: otherwise.

**V2 assay membership (study samples, not QC).**
- `pass`: the sources identify the rows as participants' biological samples, using S1 fields other than a column or label called `blank`, or using S2. Examples: a sample-collection protocol describing participant specimens; characteristics or factor values (disease group, sex, age group) on the sample rows; a statement in the description.
- `fail`: the rows are shown to be QC, blanks, standards or reference material only.
- `not_established`: otherwise.

**V3 biological count.**
- `pass` with `count`: the number of biological samples in the set is stated as an integer in a quote, or as a sum of stated integers for disjoint groups (`count_derivation`, exactly as in protocol §6).
  - It counts *samples*, not participants, unless the source equates the two.
  - **Counting rows in the sample file does not count:** the integer must appear in the quoted text.
  - Never subtract, never use percentages.
- `fail`: only if the sources state contradictory counts for the same set. Explain in `notes`.
- `not_established`: no stated count.

**V4 sampling structure.**
- `pass`: a source describes how samples were collected, and nothing indicates more than one sample per participant in the set. Example: "fasting plasma collected at enrolment". Cite that description.
- `fail`: the sources show more than one sample per participant in the set:
  - several time points or visits;
  - pre/post samples;
  - longitudinal sampling;
  - several specimens or matrices per participant within the set;
  - a time, visit or treatment factor in the sample file that varies within the same participant.
- **Repeated sample-name or source-name patterns are not evidence by themselves.** They may be bookkeeping, or a single value for all rows. A sample-file factor is evidence when its values (for example `baseline` / `week 12`) show that the same participants were sampled more than once, or the text says so.
- `not_established`: no sampling description.

## 5. Evidence rules (enforced by the validator)

- Evidence item: `{"source": <file name>, "locator": <where>, "quote": <verbatim text>}`.
  - The quote occurs verbatim in that file; whitespace differences are ignored, tabs count as spaces.
  - At least 3 tokens; no view markers (`[SECTION:`, `[TABLE:`, `[ABSTRACT]`, `[TITLE]`).
- **Privacy:** never quote sample names, source names, subject codes or any row that describes an individual (age, sex or clinical values for one person). Quote:
  - investigation-file text;
  - sample- or assay-file **column headers**;
  - **category values** that identify no one (for example `Factor Value[Time point]`, `baseline`, `Homo sapiens`, `blood plasma`);
  - publication text.
- Counts you compute from the files (distinct values, row tallies) go in `notes` as numbers only. They support your reasoning, but they are not evidence items.

## 6. Record format (one JSON object per line)

```json
{"accession": "MTBLS0000", "tranche": "A35",
 "extractor": "Claude assembler, sub-agent A35 batch K, brief v1",
 "extracted_utc": "2026-10-31T18:30:00Z",
 "sources_consulted": ["S1_i_Investigation.txt", "S1_s_MTBLS0000.txt", "S1_a_....txt"],
 "V1": {"answer": "pass", "evidence": [{"source": "...", "locator": "...", "quote": "..."}]},
 "V2": {"answer": "not_established", "evidence": []},
 "V3": {"answer": "pass", "count": 120, "evidence": [ ... ]},
 "V4": {"answer": "fail", "evidence": [ ... ]},
 "notes": "Which samples/assays form the set; computed tallies (numbers only); conflicts; reasons for not_established."}
```

- `V3.count` is present only when V3 = `pass`. Add `count_derivation` only when the integer is not in a quote.
- `extracted_utc` is the real UTC time (`date -u +%Y-%m-%dT%H:%M:%SZ`).

## 7. Conduct

- Decide from the text and the file structure alone. No outside knowledge of the study.
- When the evidence falls short, answer `not_established` and say why. But always record explicit evidence of non-human or cultured samples (V1 `fail`) and of repeated sampling (V4 `fail`).
- Run the validator until every record is VALID:
  `python3 -I /tmp/claude-0/-home-claude/4597b89f-b4a3-551d-b8bb-cf1d52b09264/scratchpad/extractA/validate.py <your jsonl> <your batch json>`
  It checks schema, verbatim quotes and count forms only.
