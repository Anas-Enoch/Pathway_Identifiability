# Stage A — MetaboLights census result record v1.2 (CLOSING RECORD)

**Stage A is closed (owner decision, 2026-10-28).** The automated MetaboLights census is complete. No further census
version will be built. Supersedes `STAGEA_CENSUS_RESULT_v1.1.md` (`bd5d2505c4ab97df8bf203c9c70e0ab2e502e2fdadc86a0b9de5ab5eae9d8753`), which superseded v1.0 (`bc86b6f925717143f842ba92e6df86b48073f541a887dc6c6ceb7d63eee2cb6f`).
The census outputs, locks and logs listed in v1.0 are unchanged, and so are the classification counts.

Changes from v1.0, made at the owner's direction:
1. The person-key statement is limited to inspected records.
2. The 183 all-blank finding is described as metadata ambiguity, not as evidence that the rows are biological samples.
3. D2 uses the owner's wording, and the observed rate is no longer treated as an estimate.
4. A35 is recorded: 19 studies are added to the documentary-recovery pool under their own origin flag.
5. The candidate file is recorded.

v1.2 changes one character from v1.1. The abbreviated hash of protocol v0.7.8 was misprinted as `b7e5aa95…f8f4` in v1.0 and v1.1. The correct value is `b7e5aa95392a6c1b3184a7fdf1db82db2eb547d02f1c2643deea316af9f3e8f4`, as in `STAGEA_FREEZE_RECORD_v1.5.md` and the signed `stageA_census_lock.json`. The archive check of 2026-10-29 found the misprint. Every other abbreviated hash in this record was re-verified against its file.

**Result.** Of 3582 MetaboLights accessions in the A17 universe (exemplars excluded), 3569 reached a terminal status.
**No study met the prespecified affirmative independence-evidence criteria: strict L1 = 0, sensitivity L1 = 0.** A35 does
not change either count. L2-confirmed and L3-confirmed: NOT_ASSESSED. D3 metadata units: 0.

## Provenance (unchanged from v1.0)

| Invocation | Census | Protocol | Sign-off | Lock file (SHA-256) |
|---|---|---|---|---|
| run 1 | v0.1.6 `7b36164b…5453` | v0.7.6 `dc4e313c…b117` | `Anas 2026-10-27T03:07:32Z` | `stageA_census_lock_run1_v0.1.6.json` `913602d3…d39bb` |
| retry pass | v0.1.7 `726ed984…477d` | v0.7.7 `f37e94c7…9922` | `Anas 2026-10-27T17:26:59Z` | `stageA_census_lock_retry_v0.1.7.json` `84c14ae1…f2aa2` |
| retry + A34 rescreen | v0.1.8 `2b9880a6…9744` | v0.7.8 `b7e5aa95…e8f4` | `Anas 2026-10-28T13:17:29Z` | `stageA_census_lock.json` `016cb35f…a1f2` |

Outputs in `github.pathway/stageA_census_full/`:
- cohorts `9c720f83…c7f`
- exclusions `51c72062…946`
- units `f6e7d338…880`
- overlap flags `071aefaa…488`
- summary `72c1f82c…52b`
- files `b90aafba…e2c`
- fetch log `242bb01e…e36`: 3 run headers, 21 050 requests, no access violation
- checkpoint `9062b713…112`: not for publication, because it holds pseudonymous hashes

## Outcome (unchanged)

| Status | Rule | n |
|---|---|---|
| EXCLUDED | M1 non-human / M4 > 1 matrix / M5 multiplicity / 3.4 no MS-NMR assay | 2321 / 67 / 60 / 1 |
| UNRESOLVED | M5 `SOURCE_UNIQUE` (≥ 50 sources: 312) | 616 |
| UNRESOLVED | M5 `SOURCE_DEGENERATE` (≥ 50 samples: 59) | 121 |
| UNRESOLVED | M2: every sample-sheet row QC-classified (D1) / ambiguous QC label (A22) / unidentified rows | 220 / 66 / 7 |
| UNRESOLVED | M3 duplicate `Sample Name` / M1 organism missing or conflicting / M4 organism part missing | 13 / 3 / 1 |
| FORMAT_REFUSED / UNRESOLVED_STRUCTURE / INDEX_UNRECOGNISED | 3.3 / A23 / 3.3 | 61 / 10 / 2 |
| TRANSPORT_FAILED / FETCH_HTTP_403 | non-terminal | 11 / 2 |
| **ELIGIBLE** | | **0** |

## The independence test (M5): 797 human, single-matrix MS (714) / NMR (83) studies

- **R1: no qualifying person-identifier columns among the inspected records.**
  - Inspected:
    - the 658 M5 studies whose source-node columns were tested against the frozen R1 pattern (616 `SOURCE_UNIQUE`, plus 42 degenerate or multiplicity records written by v0.1.7/v0.1.8);
    - the 220 all-QC studies rescreened under v0.1.8, which recorded `person_key_columns` before the M2 return.
  - Not inspected for person keys:
    - 139 v0.1.6 degenerate or multiplicity records;
    - every study that stopped before M5;
    - every column outside the R1 pattern.
  - No claim is made that person identifiers are absent from the studies or columns that were not inspected.
- **R2 was never satisfied.** Layer 1 passed in 18 of 797 studies. None of the 18 contained a T1–T4 sentence.

## Diagnosis D1 (A34): the 220 "all rows QC" studies

Recorded triggers (designated column : vocabulary token):

| Trigger | Studies |
|---|---|
| `Characteristics[Sample type]` : `blank` | 183 |
| `Characteristics[Organism]` : `reference compound` | 11 |
| `Characteristics[Sample type]` : `reference` / `pooled` / `quality control` / `solvent` / other | 26 |

In 183 studies every sample row is labelled `blank` in a sample-type field. **This suggests systematic metadata ambiguity, but does not prove that every affected row represents a biological sample.** The QC rule read the label as written.

Organisms among the 220: mouse 44, human 40 (19 with ≥ 50 rows), reference compound 21, rat 9, others.

Strict L1 = 0 does not depend on this ambiguity:
- R1: no person-key column was found in any of the 220.
- R2: no human study among them carries a CROSS_SECTIONAL design term, so R2 Layer 1 would fail.

## A35 (owner decision, 2026-10-28): APPROVED, a separately flagged documentary-recovery expansion

The documentary-recovery pool is **390** candidates:

| Origin flag | Census basis | n |
|---|---|---|
| `M5_SOURCE_UNIQUE` | UNRESOLVED, `SOURCE_UNIQUE`, ≥ 50 sources | 312 |
| `M5_SOURCE_DEGENERATE` | UNRESOLVED, `SOURCE_DEGENERATE`, ≥ 50 samples | 59 |
| `A35_ALLQC_DIAGNOSTIC` | latest record UNRESOLVED, rule M2, reason "no biological rows after QC removal", `organism_values_raw` contains `homo sapiens`, `n_sample_rows` ≥ 50 | 19 |

The candidate file is `github.pathway/stageA_closeout/stageA_documentary_candidates_v1.0.csv`:
- 390 rows, SHA-256 `c1b5a4893356c517b64ed51c8777c372c52dfbc4210bc86998e7b63be18ab0d7`;
- generated read-only from the latest checkpoint record per accession, with accessions in numeric order;
- columns: `accession, origin, census_grade, n_sources, n_rows, platform_class, matrix, census_version, title`;
- it holds public accession numbers and titles only, with no D3 hashes.

The A35 studies are not an L1 result. Their selection came from a post-census diagnostic, which is why they keep a separate origin flag. Before any documentary recovery, each one must have the following verified, from metadata and documents only:
- V1: organism;
- V2: assay membership;
- V3: the actual number of biological samples;
- V4: the sampling structure.

| Accession | Rows | Organism values | Trigger | Triage note (title or organism field only; **not evidence**) |
|---|---|---|---|---|
| MTBLS20 | 148 | blank, homo sapiens, reference compound | blank ×148 | non-human values present (V1) |
| MTBLS148 | 280 | homo sapiens | blank ×280 | collection-conditions study: repeated samples per donor possible (V4) |
| MTBLS6260 | 61 | homo sapiens | blank ×61 | ingestion intervention: time-course possible (V4) |
| MTBLS9021 | 287 | homo sapiens, mus musculus | blank ×287 | mixed organism (V1) |
| MTBLS11048 | 67 | homo sapiens, mus musculus, solvent blank | blank ×67 | mixed organism (V1) |
| MTBLS12684 | 144 | homo sapiens | blank ×144 | — |
| MTBLS12824 | 95 | homo sapiens | blank ×94, quality control ×1 | — |
| MTBLS13136 | 100 | homo sapiens | blank ×100 | — |
| MTBLS13319 | 61 | homo sapiens | blank ×61 | method study by title (V3, V4) |
| MTBLS13382 | 833 | homo sapiens | blank ×833 | — |
| MTBLS13634 | 100 | homo sapiens | blank ×100 | — |
| MTBLS13827 | 170 | homo sapiens | blank ×170 | — |
| MTBLS14358 | 325 | homo sapiens | blank ×325 | same title as MTBLS14376: possible participant overlap |
| MTBLS14376 | 126 | homo sapiens | blank ×126 | same title as MTBLS14358: possible participant overlap |
| MTBLS14421 | 53 | homo sapiens | blank ×53 | — |
| MTBLS14426 | 64 | homo sapiens | blank ×64 | — |
| MTBLS14865 | 368 | homo sapiens | blank ×368 | — |
| MTBLS15004 | 85 | homo sapiens | blank ×85 | — |
| MTBLS15165 | 70 | homo sapiens | blank ×70 | stored platelet units: repeated samples per unit possible (V4) |

All 19 are MS except MTBLS6260 (NMR).

## Diagnosis D2: transfer inconsistencies (limitation)

**Two observed transfer inconsistencies were identified, both associated with subsequently refused metadata; the
frequency of undetected truncation is unknown.**

The two cases are `s_MTBLS1873.txt` (70 965 B → 55 727 B) and `s_MTBLS8109.txt` (155 036 B → 143 726 B):
- both were shorter on the slower transfer;
- both were refused for a row-width mismatch.

They come from the 1010 files that happened to be fetched more than once. Those files are not a representative sample, and truncation exactly at a row boundary would pass the structural checks. **Any future retrieval pipeline should verify the expected transfer length where available** (for example `Content-Length`). This is an engineering lesson for later data acquisition, not a reason to repeat the census.

## Non-terminal studies (13, 0.36 %): limitation

There are 11 `TRANSPORT_FAILED` studies (MTBLS10740, MTBLS11413 …) and 2 persistent `FETCH_HTTP_403` studies (MTBLS15156, MTBLS15262). They are left as recorded and were not inspected for person keys.

## Admissible statements

- *Among 3569 MetaboLights accessions screened to a terminal status, no study met the prespecified affirmative
  independence-evidence criteria (strict and sensitivity L1 = 0). 797 human single-matrix MS/NMR studies reached the
  independence test. Among the inspected records, none carried a qualifying person-identifier column, and none documented
  a one-sample-per-participant design in the prespecified form.*
- *The census measures whether deposited metadata establish independence (E), not whether studies are independent (I).
  It found no E = 1 study and does not estimate the prevalence of I = 1.*
- Not admissible:
  - "MetaboLights contains no independent human cohort";
  - "no study has person identifiers";
  - any estimate of the true cohort count;
  - any truncation rate.

## Disposition

Stage A is CLOSED, with these limitations recorded:
- the 183 all-blank sample-type studies;
- the 13 non-terminal accessions;
- D2;
- the inspected-records scope of R1.

It may later be reported as a feasibility section or a supplementary analysis. Other decisions in force:
- A35 is approved (pool 390, with origin flags).
- Documentary recovery (§3.7) is approved for development only. Its questions, evidence standards and selection order are to be frozen before any study is examined, followed by an initial bounded tranche.
- The participant-identified controlled-access route is the owner's primary route. Its first assessment is in `EXTERNAL_ROUTE_ASSESSMENT_LIPID_ADVANCE_v1.0.md`.
