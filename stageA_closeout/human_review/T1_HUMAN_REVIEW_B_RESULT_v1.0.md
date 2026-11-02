# T1 human review, Component B: result v1.0

Date: 2026-11-02. Protocol: `T1_HUMAN_REVIEW_PROTOCOL_v1.0.md` (799865e8…f509), Component B. **Non-blind by design.** Sensitivity analysis only; no frozen record was changed.

## 1. Provenance and attestation

- Reviewer: the author. Attestation, verbatim, 2026-11-02 00:38 local: *"i read the source files or the 12 studies and wrote the verdicts, quotes and reasons myself"*.
- Files: `T1_HUMAN_REVIEW_B_completed.numbers` (d399d4bd90f3f4db8c85e5eb456e342608a9fb14ad2101419069f5d558905ca0; the reviewer's original, delivered under the filename `T1_AI_REVIEW_B_completed.numbers`) and its CSV export `T1_HUMAN_REVIEW_B_completed.csv` (14fdade9fb66ceba3f628bac78301af044add0fa414819a2cc73af7ab7cda939).
- Checks by the assembler: every quotation is verbatim in the cited sealed file; every cited file is on the row's permitted list; verdict vocabulary is the protocol's.
- **Disclosures:**
  1. The reviewer is the author, who adjudicated these same 12 records on 2026-10-31 (the protocol anticipated this: Component B is non-blind).
  2. The `reviewed_utc` values were not typed at completion: eleven of twelve end in `:29Z` and they are not chronological, so they are assembly artefacts and per-study timing is unavailable.
  3. The delivered filename contained "AI"; the author attests the content is his own.

## 2. Result

| Verdict | n | Studies |
|---|---|---|
| confirmed | 11 | MTBLS89, 92, 200, 2291, 3800, 10634, 12326, 12372 (repeated sampling); MTBLS419, 2205, 12861 (non-participant samples) |
| disputed | 0 | — |
| unresolved | 1 | MTBLS10972 |

In the 11 confirmed cases the author's quotation overlaps the frozen record's cited evidence, i.e. the same source sentence supports the classification.

## 3. The unresolved case

**MTBLS10972** (frozen: CONTRADICTED, REPEATED_SAMPLING; S1 only, no publication text). The author's reason: the investigation text describes a post-meal metabolomic response and changes at different times, but does not explicitly state that multiple plasma specimens from the same individuals were included in the deposited LC-MS assay. The frozen extractor's evidence was the time-course wording ("transiently decreased at 240 and 300 min", "elevated amino acids at 6 h", a `timepoint` study factor), and the extractor itself marked the answer "weaker"; the assembler's adjudication recommendation was "confirm (weaker)", and the author confirmed it on 2026-10-31.

On re-reading, the author now finds the evidence short of an explicit per-person repeated-sampling statement. Under the frozen protocol this would make the record DOC_UNRESOLVED rather than DOC_CONTRADICTED. The frozen grade stands; the change is reported here as a sensitivity result.

## 4. Bearing on the findings

- Repeated-sampling count: 9 of 40 (frozen) versus 8 of 40 under the author's re-reading (22.5% → 20%). Non-participant samples: 3 of 40, unchanged.
- No effect on x = 0, DR-7, or the STRICT count. The CONTRADICTED total would be 11 rather than 12.
- The manuscript should report both figures: the frozen 9/3 as the primary, and 8/3 as the author's re-review sensitivity.

## 5. Limitations

Non-blind; same reviewer as the original adjudication; no independent reviewer; timing artefacts as in §1.
