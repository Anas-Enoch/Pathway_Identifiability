# Public-cohort screen log (PUBLIC_COHORT_SCREEN_v0.1.md, frozen at commit f53d8b6, sha256 da7f80e4…)

Metadata only. No quantitative file was downloaded or opened for any entry. Each criterion carries PASS / FAIL /
UNKNOWN, a source URL, and the passage or field the verdict rests on. Entries are never deleted or reordered.
Retrieval by Claude (Fable 5.1) through web fetch of the pages the owner supplied; the owner's own reading of the
same pages on 4 Nov 2026 is recorded where it adds a fact.

## 0. Prior-contact register for E6 (programme history, to be checked against every candidate)

| Source | Records | Nature of contact |
|---|---|---|
| Stage A census (28 Oct 2026) | all MetaboLights accessions public on the census date (3,582) | metadata only; T1 documentary tranche of 40; A35 19 |
| B0 census / P4V B3 (17–25 Oct 2026) | Metabolomics Workbench universe on 17 Oct 2026 | metadata only |
| 186-record archive (Paper III/IV) | 188 Workbench analyses | values opened; method development |
| Paper I cohorts | ST000356, ST001849, ST001865, ST002829, ST003390, ST003506 | values opened |
| P4U / Stage 10K | ST004571, ST002733, ST003368, ST000568, ST004168 | values opened or metadata |
| Stage 16L (11 Oct 2026) | NHANES 2017–18 PHTHTE_J | values read for status audit |
| **Stage 16L real screen (11 Oct 2026)** — `stage16l_real_screen/SCREENING.md`, directly inspected in Terminal on 4 Nov | EDRN prostate-cancer metabolomics; Emory Cardiovascular Biobank (Dryad 10.5061/dryad.866t1g1mt); BIDMC MASLD; ADNI1 metabolomics | Metadata/documentation-only dataset-discovery screen; no outcome values accessed, model trained, or predictive comparison performed (source excerpts below). |
| NHANES-FA (4 Nov 2026) | NHANES 2011–12 FAS_G/DEMO_G | values opened (development); 2013–14 never opened |

**Source excerpts for Stage 16L real-screen prior contact** (verbatim from the Terminal display of `stage16l_real_screen/SCREENING.md`, dated 2026-10-11):

> **Verdict: no real dataset qualified. No outcome values accessed, no candidate model trained, no predictive comparison performed.** The repositories/articles were reviewed at the metadata/documentation level.
>
> ### EDRN prostate cancer metabolomics (NCI / Weill Cornell; Metabolomics Workbench DOI 10.21228/M86H7K) — highest priority for a *data-integrity investigation*, NOT Stage 16L PASS
>
> ### Emory Cardiovascular Biobank metabolomics (Dryad 10.5061/dryad.866t1g1mt) — two-cohort raw-data fallback
>
> ### BIDMC metabolomic MASLD cohort (PhysioNet version 1.0.0) — smaller public fallback
>
> ### ADNI1 targeted metabolomics (Scientific Data 2017) — attractive LOD documentation but controlled access

The source additionally states, regarding Emory: "**Strict Stage 16L FAIL/PENDING**: explicit per-cell detection states and threshold provenance not verified, feature/target comparability to frozen Paper IV unverified, independence from development list unverified, local hashes not recorded." The `stage16l_real_screen/screen_status.json` records `DRYAD_EMORY_866t1g1mt` with DOI `10.5061/dryad.866t1g1mt` and status `PENDING_RAW_AND_LOD_DOCUMENTATION`. These quotes establish prior metadata-only contact for the four named candidates; they do **not** establish contact with either OMIX009009 or Dryad `h44j0zpz0`. The Stage 16L screen and this quotation are distinct from the NHANES urinary-phthalate audit.

Repository grep on 4 Nov 2026 for `OMIX009009|OMIX0090|ngdc|dryad` (md/py/json/csv/txt, excluding .git and
virtual environments) matched: this protocol; `stage16l_real_screen/screen_status.json` and `SCREENING.md` (the
Emory Dryad record above); two MetaboLights investigation files under `stage10E_metadata_prescreen_raw`
(MTBLS1720, MTBLS15778) whose text references other Dryad deposits. `git log --all` matched nothing for
`omix|dryad`. Neither candidate below appears in any prior record.

## Entry 01 — NGDC OMIX009009 (preidentified before the freeze; slot 1 of 20)

Record: https://ngdc.cncb.ac.cn/omix/release/OMIX009009 (fetched 4 Nov 2026). Title: "Untargeted plasma metabolomics
for lung cancer in two prospective cohort studies". Organism Homo sapiens; sample number 704; data type "Metabolome
Data by Mass Spectrometry (MS)"; submitted 2025-02-13, released 2026-05-13; open access; BioProject PRJCA036104; no
BioSample; no publication listed; submitter Meng Zhu, Nanjing Medical University.

Description (verbatim segments): "This study was based on two prospective community cohorts of chronic
non-communicable diseases in Jiangsu Province, China." "Using a 1:1 matched nested case-control study design,"
"untargeted metabolomics analysis was performed on plasma samples collected at baseline from 704 participants."
"(analyzed using LC-MS) and lung cancer risk."

File list (verbatim): OMIX009009-01 `NEGa`, 704 samples, rar, 7.13 GB; -02 `NEGb`, 704, rar, 6.67 GB;
-03 `Phenotype`, 1, xlsx, 185.65 KB; -04 `POSa`, 704, rar, 4.01 GB; -05 `POSb`, 704, rar, 4.29 GB.

| Criterion | Verdict | Evidence |
|---|---|---|
| E1 participant structure | UNKNOWN | "plasma samples collected at baseline from 704 participants" (one baseline sample per participant is stated at enrolment level); the `Phenotype` xlsx is not described, so a sample-to-participant mapping for the deposited runs is not established; no BioSample records. |
| E2 ≥ 60 jointly measured metabolites | UNKNOWN (not assessable) | no processed matrix is deposited; a feature count would be the product of our own processing of the raw runs. |
| E3 participant-disjoint evaluation | UNKNOWN | "two prospective community cohorts" named, but the deposit does not state which samples belong to which cohort; the "1:1 matched nested case-control" design introduces case–control pair dependence that a partition must respect. |
| **E4 accessible quantified matrix** | **FAIL** | the only data files are four RAR archives of raw LC-MS acquisitions (two ionisation modes × two batches, 22 GB in all); no processed or annotated metabolite table is listed. Producing one requires peak picking, alignment and annotation by us, which E4 excludes as ad hoc raw-spectra reconstruction. |
| E5 feasibility | UNFAVOURABLE (recorded, not decisive) | 704 participants in total across both cohorts and both case/control arms; any cross-cohort split yields ≈ 350 per population, below the 1,500 development minimum used in NHANES-FA (the threshold for a new cohort would be fixed separately, but the order of magnitude is short). |
| E6 no prior contact | PASS | repository grep and git log, §0. |
| E7 shared platform | PASS (within deposit) | single LC-MS acquisition series for all 704 samples. |
| E8 no outcome-driven feature selection | PASS for features; note | untargeted acquisition, no selected panel; participant *selection* is outcome-driven (nested case–control for lung cancer), which affects E3/E5, not E8. |

**First failure: E4. Disposition: FAIL — not eligible under v0.1.** The record is a raw-data deposit without a
quantified matrix; it could become eligible only through a future release of a processed matrix with a sample
mapping, which would be a new record for a new screening cycle.

## Entry 02 — Dryad 10.5061/dryad.h44j0zpz0 (preidentified before the freeze; slot 2 of 20)

Record: https://datadryad.org/dataset/doi:10.5061/dryad.h44j0zpz0 (fetched 4 Nov 2026). Title: "Metabolomic aging and
mortality clocks based on ensemble learning identify systemic aging risk across the human lifespan in 18,476
individuals"; Pu J, Qian K, Xu W, Cao L, Chen Y, Jiang H, Huang H, Xu Y, Zhang W; published on Dryad 1 Sep 2026;
no related publication listed; licence not shown beyond Dryad's reuse statement.

Files (verbatim, version of 1 Sep 2026, 13.29 MB): `nom_diseased_val_7768_56.csv` 5.59 MB;
`nom_healthy_train_5252_56.csv` 3.78 MB; `nom_healthy_val_5456_56.csv` 3.92 MB; `README.md` 2.45 KB.

| Criterion | Verdict | Evidence |
|---|---|---|
| E1 participant structure | UNKNOWN | README: data "de-identified"; abstract: 18,476 participants aged 1–96 at six centres; whether each row is one individual and whether any individual appears in more than one matrix is not stated on the page. |
| **E2 ≥ 60 jointly measured metabolites** | **FAIL** | abstract: the study "identified 56 representative age-associated metabolites"; file names carry `_56`; the released matrices hold 56 features. |
| E3 participant-disjoint evaluation | PASS (provisional) | three released matrices: healthy training (5,252), healthy validation (5,456), diseased validation (7,768); disjointness between them is implied by the design, not stated. |
| E4 accessible quantified matrix | PASS | processed CSVs; README: spectra "smoothed, baseline-corrected, peak-detected, aligned, and filtered, then batch-corrected with pyCombat. Values are log-transformed and z-score normalized." |
| E5 feasibility | PASS on size | 18,476 individuals. |
| E6 no prior contact | PASS | §0 (the Emory Dryad record 866t1g1mt is a different deposit). |
| E7 shared platform | PASS, with a discrepancy noted | README: Bruker Autoflex III MALDI-TOF/TOF, positive-ion reflector, m/z 100–1000; abstract: "nano-assisted laser desorption/ionization mass spectrometry". One platform family (LDI-MS), two descriptions. |
| **E8 no outcome-driven feature selection** | **FAIL** | README: the 56 features "were selected by consensus of Random Forest, Lasso regression, and mutual information methods" for age association. The released panel is outcome-selected. |

**First failure: E2 (E8 also fails). Disposition: EXCLUDE.** A project repository reportedly mentions a 232-feature
matrix; it is not part of this Dryad release and does not alter this entry. If such a matrix is ever deposited, it
is a new record for a new cycle.

## Register

| Entry | Record | Preidentified | First failure | Disposition |
|---|---|---|---|---|
| 01 | NGDC OMIX009009 | yes | E4 | FAIL |
| 02 | Dryad 10.5061/dryad.h44j0zpz0 | yes | E2 (and E8) | EXCLUDE |

Slots consumed: 2 of 20. Discovery slots (18, three per family) await the hashed enumeration script. No candidate
has been authorised for any value access.
