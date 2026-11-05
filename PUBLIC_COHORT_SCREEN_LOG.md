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

## Entry 03 — Dryad 10.5061/dryad.gqnk98srt (Dryad discovery slot 1 of 3)

Record: https://datadryad.org/dataset/doi:10.5061/dryad.gqnk98srt (fetched 5 Nov 2026). "The multiomic landscape of
epidemiological factors contributing to preterm birth in low- and middle-income countries"; Espinosa C, Khan W,
Khanam R et al.; published on Dryad 19 May 2023; Stanford University; no related publication listed; licence not
named beyond Dryad's reuse statement. Selected by the frozen hash order of the Dryad index (sort key in
`screen_enum_v010_run1/screen_enum_dryad.json`).

Abstract (verbatim fragment): "Plasma samples from 231 participants were analyzed to generate proteomic,
metabolomic, and lipidomic datasets." Maternal covariates were collected from 13,841 pregnant women across five
sites; the omics subset is the 231.

Files (verbatim, 18.11 MB): `clinical_targets.csv` 6.78 KB; `epidemiological_model.R` 452 B;
`feature_correlations.csv` 2.94 MB; `feature_table.csv` 14.58 MB; `maternal_covariate_dictionary.xlsx` 28.90 KB;
`metabolome_key.csv` 263.36 KB; `multiomics_feature_key.csv` 287.92 KB; `README.txt` 3.84 KB.

| Criterion | Verdict | Evidence |
|---|---|---|
| E1 participant structure | UNKNOWN | "Plasma samples from 231 participants" is an enrolment-level statement; the record page shows empty Methods and Usage-notes sections and does not display the README, so a sample-to-participant mapping in `feature_table.csv` is not established. |
| E2 ≥ 60 metabolites | UNKNOWN (favourable) | a `metabolome_key.csv` of 263 KB and a 14.6 MB feature table are consistent with hundreds or thousands of features; no count is stated on the page and none was inferred by opening a file. |
| E3 participant-disjoint evaluation | UNKNOWN | five sites in one study; no second population is described; any split would be a within-cohort holdout, not cross-cohort transport. |
| E4 accessible quantified matrix | PASS | `feature_table.csv` with `metabolome_key.csv` and `multiomics_feature_key.csv` is a processed matrix with feature keys. |
| **E5 feasibility** | **FAIL** | the omics population is 231 participants in total. The acquisition comparison needs a development population and an untouched evaluation population; 231 split between them is far below the participant counts at which this programme has previously been able to observe any selector contrast (NHANES-FA development used 1,725 and observed none). No numeric MDE is computed here, as none can be computed from an unobserved contrast (protocol E5). |
| E6 no prior contact | PASS | §0. |
| E7 shared platform | UNKNOWN | the platform is not stated on the record page. |
| E8 no outcome-driven feature selection | PASS (as far as shown) | untargeted multiomic profiling; the page describes no selected panel in the deposited table. |

**Identifiable failures: E5. Disposition: FAIL.**

## Entry 04 — Dryad 10.5061/dryad.p2ngf1w1z (Dryad discovery slot 2 of 3)

Record: https://datadryad.org/dataset/doi:10.5061/dryad.p2ngf1w1z (fetched 5 Nov 2026). "Multiomics characterization
of acute child illness and mortality in Africa and South Asia"; Espinosa Bernal C, Njunge J, Tickell K et al.;
published 6 Jan 2026, updated 10 Feb 2026; no related publication listed. Abstract: 3,101 acutely ill children at
nine sites; a nested case-cohort of 1,008; 270 community children; an independent validation cohort of 100; samples
at admission and at discharge; plasma proteomics, serum metabolomics and lipidomics, stool metagenomics.

Files (verbatim, both versions, 97.3 KB in total): `packages.csv` 3.33 KB; `README.md` 8.89 KB;
`scripts.zip` 85.08 KB.

| Criterion | Verdict | Evidence |
|---|---|---|
| E1 participant structure | UNKNOWN | described at the study level ("admission and discharge" implies repeated sampling); no deposited sample table. |
| E2 ≥ 60 metabolites | UNKNOWN (not assessable) | no measurement file is deposited. |
| E3 participant-disjoint evaluation | FAVOURABLE if data existed | a nested case-cohort plus "an independent cohort of 100 children" is the structure the programme needs; it is described, not released. |
| **E4 accessible quantified matrix** | **FAIL** | the deposit contains only analysis scripts, a package list and a README (97 KB in all). No metabolite matrix, no sample table, no participant key. |
| E5 feasibility | not reached | — |
| E6 no prior contact | PASS | §0. |
| E7 shared platform | UNKNOWN | not stated on the record page. |
| E8 no outcome-driven feature selection | UNKNOWN | not stated. |

**Identifiable failures: E4. Disposition: FAIL.** The underlying study's design would satisfy E3; the deposit does
not release the data. A future deposit of its matrices would be a new record for a new cycle.

*Dryad discovery slot 3: no candidate.* The frozen index query returned 2 records in total
(`screen_enum_dryad.json`, total 2, both topical-eligible, both screened above). The slot is recorded as
**NO FURTHER ELIGIBLE RECORD IN INDEX**, not as unsearchable, and is not refilled from another family.

## Entry 05 — figshare 10.6084/m9.figshare.33080369.v1 (figshare discovery slot 1 of 3)

Record: https://api.figshare.com/v2/articles/33080369 (owner-fetched 5 Nov 2026; `figshare_33080369.json`, hashed in
`screen_enum_v011/meta/SHA256SUMS.txt`). "Dataset on the changes in the gut microbiome and serum metabolome in type 2
diabetic patients with cognitive decline"; posted 25 Jul 2026; licence CC BY 4.0; no related article DOI recorded.
Description: "we analyzed fecal microbiota and serum metabolomic profiles in healthy controls, patients with type 2
diabetes and DACD patients by using 16S rRNA sequencing and LC-MS-based metabolomics analyses".

Files (verbatim): `General Information.docx` 15.3 KB; `Protocol for LC-MS metabolomics analysis.docx` 15.6 KB;
`Sample list_Metabolomics raw data.xls` 23.6 KB; `CON_neg.rar` 1.41 GB; `T2D_neg.rar` 1.41 GB; `DACD_neg.rar`
1.40 GB; `CON_pos.rar` 1.78 GB; `T2D_pos.rar` 1.73 GB; `DACD_pos.rar` 1.70 GB.

| Criterion | Verdict | Evidence |
|---|---|---|
| E1 participant structure | UNKNOWN | three groups are named; no participant counts or sample-to-participant mapping appear in the record metadata; the `Sample list` xls is a manifest of raw files, not opened. |
| E2 ≥ 60 metabolites | UNKNOWN (not assessable) | no processed matrix; a feature count would come from our own processing. |
| E3 participant-disjoint evaluation | UNKNOWN | one study, three diagnostic groups; no second population. |
| **E4 accessible quantified matrix** | **FAIL** | the six data files are RAR archives of raw LC-MS acquisitions (≈ 9.4 GB), with two documentation files and a sample list. No processed metabolite table is deposited. Same disposition as entry 01. |
| E5 feasibility | not reached | — |
| E6 no prior contact | PASS | §0. |
| E7 shared platform | PASS (within deposit) | one LC-MS protocol, positive and negative modes. |
| E8 no outcome-driven feature selection | PASS | untargeted acquisition. |

**Identifiable failures: E4. Disposition: FAIL.**

## Entry 06 — figshare 10.3389/fmicb.2026.1801806.s003 (figshare discovery slot 2 of 3)

Record: https://api.figshare.com/v2/articles/32357331 (owner-fetched 5 Nov 2026). "Table 1_A comparative study of the
gut microbiome and fecal metabolome in hypertensive patients from middle-temperate and tropical cities of China:
Daqing and Haikou.xlsx"; posted 22 May 2026; CC BY 4.0; supplement of Frontiers in Microbiology article
10.3389/fmicb.2026.1801806. Description: "Fecal samples were collected from 28 hypertensive patients in Daqing (DQ
group) and 32 in Haikou (HK group), and analyzed using shotgun metagenomic sequencing and untargeted metabolomics."

Files (verbatim): one file, `Table 1_… .xlsx`, 12.1 KB.

| Criterion | Verdict | Evidence |
|---|---|---|
| E1 participant structure | UNKNOWN | 28 + 32 patients stated; the deposit is a single article table, with no sample-to-participant mapping described. |
| E2 ≥ 60 metabolites | UNKNOWN (not assessable) | a 12 KB table is not a metabolite matrix; its contents are not described beyond "Table 1". |
| E3 participant-disjoint evaluation | UNKNOWN | two cities are described, which would be two populations, but no deposited data distinguish them. |
| **E4 accessible quantified matrix** | **FAIL** | the record is one supplementary table of an article (12.1 KB). The study's metabolomic data are not deposited here. |
| **E5 feasibility** | **FAIL** | 60 participants in total across both cities. |
| E6 no prior contact | PASS | §0. |
| E7 shared platform | UNKNOWN | one untargeted platform implied; not specified in the record. |
| E8 no outcome-driven feature selection | UNKNOWN | not stated. |

**Identifiable failures: E4, E5. Disposition: FAIL.**

## Entry 07 — figshare 10.3389/fendo.2026.1878728.s002 (figshare discovery slot 3 of 3)

Record: https://api.figshare.com/v2/articles/33194799 (owner-fetched 5 Nov 2026). "Data Sheet 1_A lifestyle-driven
inflammatory load index and its circulating metabolomic signature identify susceptibility to hospital-recorded acute
respiratory infections: a population-based cohort study of 295,679 adults.docx"; posted 10 Aug 2026; CC BY 4.0;
supplement of Frontiers in Endocrinology article 10.3389/fendo.2026.1878728. Description: "This prospective cohort
study included 295,679 UK Biobank participants after prespecified exclusions… a nuclear magnetic resonance
(NMR)–derived metabolomic risk signature (MRS)."

Files (verbatim): one file, `Data Sheet 1_… .docx`, 5.97 MB.

| Criterion | Verdict | Evidence |
|---|---|---|
| E1 participant structure | UNKNOWN | the analysed cohort is UK Biobank; this deposit contains no participant or sample table. |
| E2 ≥ 60 metabolites | UNKNOWN (not assessable) | the deposit is a Word document. |
| E3 participant-disjoint evaluation | UNKNOWN | — |
| **E4 accessible quantified matrix** | **FAIL** | one supplementary `.docx`; the underlying measurements are UK Biobank data held under that resource's own access process and are not in this record. |
| E5 feasibility | not reached | — |
| E6 no prior contact | PASS | §0. |
| E7 shared platform | UNKNOWN | — |
| **E8 outcome-driven feature selection** | **FAIL (as far as described)** | the described quantity is a "metabolomic risk signature" derived for an outcome; a signature is an outcome-selected construct, not an acquisition panel. |

**Identifiable failures: E4, E8. Disposition: FAIL.**

*Observation, not a candidate.* This record points at UK Biobank NMR metabolomics, a resource with participant-level
identifiers and a large panel. It is **not** added to this cycle: the screened record is the supplement, the
underlying resource was not reached through the frozen enumeration, and its access route is an application, not a
public deposit. Any use of it would require its own prospective amendment and access approval, like LIPID/COMETS.

## Entry 08 — Metabolomics Workbench ST005211 (Workbench discovery slot 1 of 3)

Metadata: `rest/study/study_id/ST005211/{summary,analysis,factors,number_of_metabolites}` (owner-fetched 5 Nov 2026;
hashed in `screen_enum_v011/meta/SHA256SUMS.txt`). "Untargeted Metabolomic Profiling Reveals a Distinct Intestinal
Metabolic Signature in Necrotizing Enterocolitis"; University of Illinois Chicago; *Homo sapiens*; LC-MS;
`number_of_samples` 21; released 2026-10-21; CC BY-NC-ND.

| Criterion | Verdict | Evidence |
|---|---|---|
| E1 participant structure | PARTIAL / UNKNOWN | the `factors` table carries patient labels inside a treatment field (`Treatment:CONTROL PT 1`, `CONTROL PT 2`, `CONTROL PT 3`, …), so a participant key is partially recoverable; it is a free-text factor value, not a documented participant identifier, and the deposit does not state how many patients the 21 samples come from. |
| E2 ≥ 60 metabolites | PASS | `num_metabolites` 291 (AN008903) and 202 (AN008904). |
| **E3 participant-disjoint evaluation** | **FAIL** | one study, 21 samples; no second population. |
| E4 accessible quantified matrix | PASS (presumed) | Workbench deposits carry named-metabolite tables; not opened. |
| **E5 feasibility** | **FAIL** | 21 samples in total, drawn from fewer patients than that; a development population and an untouched evaluation population cannot both be formed. No numeric MDE is computed. |
| E6 no prior contact | PASS | §0; ST005211 is absent from the B0 universe (released after 17 Oct 2026). |
| E7 shared platform | PASS within the deposit | one LC-MS series. |
| E8 no outcome-driven feature selection | PASS | untargeted. |

Noted for the record: `sample_source` spans `jejunum`, `ileum`, `colon`, `cecum`, `descending colon` and the pair
`jejunum, ileum`, i.e. several intestinal sites, in places more than one per patient. The spelling `jejujum` appears
in one row. **Identifiable failures: E3, E5. Disposition: FAIL.**

## Entry 09 — Metabolomics Workbench ST005210 (Workbench discovery slot 2 of 3)

Metadata as above for ST005210. "Targeted and untargeted metabolomics and lipidomics of paired epicardial adipose
tissue and left ventricular myocardium from 47 human donor hearts"; University at Buffalo; *Homo sapiens*; LC-MS;
`number_of_samples` 94; released 2026-10-21.

| Criterion | Verdict | Evidence |
|---|---|---|
| E1 participant structure | UNKNOWN | no participant identifier field. `local_sample_id` values (`D50A`, `D12A`, `D25A`, …) suggest a donor number with a tissue letter, which would make the donor key recoverable by convention; the deposit documents no such convention, so the mapping is inferred, not established. |
| E2 ≥ 60 metabolites | PASS | `num_metabolites` 275 (AN008898), 206 (AN008899). |
| **E3 participant-disjoint evaluation** | **FAIL** | the 94 samples are two tissues from the same 47 donors (title: "paired epicardial adipose tissue and left ventricular myocardium from 47 human donor hearts"); the two groups are matrices, not independent populations. |
| E4 accessible quantified matrix | PASS (presumed) | not opened. |
| **E5 feasibility** | **FAIL** | 47 donors; 94 samples is not 94 independent units. |
| E6 no prior contact | PASS | §0; released after the B0 date. |
| **E7 shared platform** | fails as a two-population design | one platform, but the two sample groups are different tissues; cross-tissue transport is not the frozen task. |
| E8 no outcome-driven feature selection | PASS | targeted + untargeted panels, not outcome-selected. |

Noted for the record, and relevant to the documentation-audit manuscript: the `sample_source` field reads `_QC_` on
**all 94 rows**, while the actual matrix appears in the factor string (`Sample source:epicardial adipose tissue |
Tissue:…`). This is the same class of template-default sample-type annotation the audit paper documents for
MetaboLights `blank` labels, in a different repository. **Identifiable failures: E3, E5. Disposition: FAIL.**

## Entry 10 — Metabolomics Workbench ST005223 (Workbench discovery slot 3 of 3)

Metadata as above for ST005223. "Steady state metabolomics of fresh and 12-day stored red blood cells from mice
expression human canonical (human WT) G6PD or human G6PD variants…"; University of Colorado Denver;
**species: *Mus musculus***; `number_of_samples` 29; released 2026-10-21.

| Criterion | Verdict | Evidence |
|---|---|---|
| **E1 participant structure** | **FAIL** | the study organism is *Mus musculus*; there are no human participants. The 29 rows comprise 18 `red blood cells`, 7 `_BLANK_` and 4 `_QC_`. |
| E2 ≥ 60 metabolites | PASS numerically | 76 (AN008927), 88 (AN008928). |
| E3 / E5 | FAIL | 18 biological samples, one experiment. |
| E6 no prior contact | PASS | §0. |
| E7 / E8 | not reached | — |

**Identifiable failures: E1 (organism), E3, E5. Disposition: FAIL.**

**Recorded weakness of the v0.1 topical filter.** This record entered its slot because the frozen filter requires a
metabolomics term and a human/cohort term in title+summary, and the title contains "human canonical (human WT)
G6PD" describing a transgene, not the study organism. The filter does not read the `species` field at enumeration
time. The slot is **not** refilled and the entry stands with its FAIL; a species check at enumeration would be a
change to the selection rule and therefore a matter for a future prospective amendment, not this cycle.

## Register

| Entry | Record | Family / slot | Preidentified | Identifiable failures | Disposition |
|---|---|---|---|---|---|
| 01 | NGDC OMIX009009 | preidentified | yes | E4 | FAIL |
| 02 | Dryad 10.5061/dryad.h44j0zpz0 | preidentified | yes | E2, E8 | EXCLUDE |
| 03 | Dryad 10.5061/dryad.gqnk98srt | dryad 1/3 | no | E5 | FAIL |
| 04 | Dryad 10.5061/dryad.p2ngf1w1z | dryad 2/3 | no | E4 | FAIL |
| — | (no further record in the Dryad index) | dryad 3/3 | — | — | NO FURTHER ELIGIBLE RECORD |
| 05 | figshare 10.6084/m9.figshare.33080369.v1 | figshare 1/3 | no | E4 | FAIL |
| 06 | figshare 10.3389/fmicb.2026.1801806.s003 | figshare 2/3 | no | E4, E5 | FAIL |
| 07 | figshare 10.3389/fendo.2026.1878728.s002 | figshare 3/3 | no | E4, E8 | FAIL |
| 08 | Workbench ST005211 | workbench 1/3 | no | E3, E5 | FAIL |
| 09 | Workbench ST005210 | workbench 2/3 | no | E3, E5 | FAIL |
| 10 | Workbench ST005223 | workbench 3/3 | no | E1 (organism), E3, E5 | FAIL |

**Slot accounting (not a final cycle result).**

| Family | Slots | State |
|---|---|---|
| preidentified | 2 | screened: entries 01–02, both ineligible |
| Dryad | 3 | 2 screened (03, 04); 1 with no further record in the frozen index |
| figshare | 3 | screened: entries 05–07 |
| Metabolomics Workbench | 3 | screened: entries 08–10 |
| Zenodo | 3 | **pending**: enumeration failed in runs 1 and 2 (HTTP 400; cause identified 5 Nov as the unauthenticated page-size limit of 25, requested 100); script v0.1.2 corrects the request mechanics only |
| NGDC OMIX | 3 | **UNSEARCHABLE**: list endpoint returned HTTP 404; no public search interface reached |
| MetaboLights | 3 | **UNSEARCHABLE**: 34 accessions absent from the Stage A universe were identified from the LIST endpoint, but the per-accession `i_Investigation.txt` reads needed to apply the release-date rule were refused (run 1: 28 of 34 failed, `Connection refused`; run 2 with paced retries: `RemoteDisconnected`) |

Records screened to date: **10**, of which **0 qualify**. Three Zenodo slots are pending and six slots are
unsearchable. **No cycle result is recorded yet**; a cycle result requires the Zenodo slots to be run or recorded
as unsearchable, and it will state the unsearchable coverage explicitly, because a bounded search that could not
enumerate three of six families is not evidence about those families.

No quantitative measurement file has been downloaded or opened for any entry in this log.
