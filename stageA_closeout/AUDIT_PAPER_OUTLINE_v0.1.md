# Audit paper: outline v0.1 (target: *Metabolomics*, Springer)

Status: outline for owner approval. No LaTeX until it is approved. Date: 2026-10-31.
Every number below is taken from a hashed record and cited by file. Items marked **[OPEN]** must be settled before drafting.

## Working title

*Do public metabolomics repositories document participant independence? A repository-wide census and a randomised documentary audit*

## The one claim

Neither the deposited metadata nor the linked publications establish that each participant contributed one sample to the analysed set, which is the property that study-level external validation and meta-analysis require. Where sampling structure is documented, dependence is common. The paper measures **documented** independence (E), not true independence (I), and says so in the abstract.

## What is already known, and what is new

| Already published (cite) | New here |
|---|---|
| Poor minimum-information (MSI) compliance in metabolomics repositories (Spicer, Salek & Steinbeck 2017, *Sci Data*); metadata quality scoring (Ferreira et al. 2017); duplicate MetaboLights deposits (Ramos et al., MTBLS107–111) | A target property defined operationally (documented one sample per participant) and applied to a whole repository (3,582 accessions) |
| Sample-ID mismatches between sections of Workbench mwTab files (Powell & Moseley 2021, *Metabolites*) | A randomised documentary audit with design-based bounds, reading deposit metadata **and** linked methods sections |
| Hidden duplicates inflate cross-study validation in transcriptomics (Waldron et al. 2016, *JNCI*); sample misannotation in GEO (Toker et al. 2016) | Rates of explicit dependence (time series, pre/post, multiple specimens, cell lines) among human deposits that pass naive filters |
| dbGaP requires a subject–sample mapping file | A template-default label ("experimental blank") that makes participant samples look like QC blanks; a sample-to-measurement linkage defect |
| Reuse of MetaboLights data is common (Karaman et al. 2026, preprint; Mendez et al. 2019 as an example of row-as-participant reuse) | The consequence: study-level external validation from repositories is structurally blocked under an affirmative-evidence standard |

## Structure

### 1. Introduction (about 700 words)
- **Why it matters.** External validation needs independent study units, and each unit needs independent participants. Repository rows are not participants.
- **What is known.** Metadata quality is known to be poor, and duplicates are known to inflate validation in other omics.
- **The gap.** Nobody has measured whether metabolomics repositories document the structure that validation needs.
- **Aims.** (i) Census MetaboLights. (ii) Audit Metabolomics Workbench subject metadata. (iii) Run a randomised documentary audit with publications. (iv) Diagnose deposits that automated filters misclassify.

### 2. Methods
- **2.1 Definitions.**
  - Unit: the deposited assay.
  - Affirmative-evidence standard: what counts as evidence, and why unique sample labels do not.
  - E vs I.
  - Grades: STRICT / COUNT / CONTRADICTED / UNRESOLVED.
- **2.2 MetaboLights census.**
  - Universe and exclusions; rules M1–M5, R1 (person-key columns) and R2 (design sentence templates).
  - Value-blind access: no measured values or data files.
  - The frozen protocol and code, with hashes.
  - Source: `STAGEA_CENSUS_RESULT_v1.2.md`.
- **2.3 Metabolomics Workbench subject audit.**
  - The 104 candidate human studies and their factor tables; subject-key patterns.
  - Precise scope: REST `factors` endpoint only. mwTab files were not opened, because they embed values. The optional mwTab Subject field therefore was not assessed.
  - Source: `P4V_B3_EVIDENCE_AUDIT_RESULT_v1.0.md`.
- **2.4 Randomised documentary audit (T1).**
  - Population: 371 candidates.
  - Sampling: a NIST Randomness Beacon pulse fixed after the lock was signed, giving a simple random sample of 40.
  - Sources: deposit metadata, named publications, citing papers, and owner-supplied texts of named publications.
  - Value-blind reading views (methods, abstract with Results removed, participant tables).
  - Document freeze, then seal signature, then extraction.
  - Extraction blind to deposit row counts; deterministic grading with verbatim-quote verification; owner adjudication of every CONTRADICTED or STRICT grade.
  - A futility rule prespecified before sampling.
  - Source: `DOC_RECOVERY_PROTOCOL_v0.2.md`.
- **2.5 Pre-verification of misclassified deposits (A35).**
  - 19 human deposits whose every row an automated QC rule removed; checks V1–V4.
- **2.6 Statistics.** Exact hypergeometric bounds for a finite population (design-based); Clopper–Pearson for comparison.
- **2.7 Use of AI tools.**
  - Evidence extraction was performed by Claude (Anthropic) sub-agents under a fixed written brief.
  - Every quote was machine-verified against the sealed texts, and every adjudicable grade was checked by the author.
  - The brief, validator and raw outputs are public.
  - Required by Springer policy; reviewers will scrutinise it (see §6).
- **2.8 Integrity.**
  - Everything was hashed and signed before execution, but not registered externally. Say so.
  - No measured metabolite value was read by any step.
  - Full texts are not redistributed.

### 3. Results
- **3.1 Census.**
  - 3,569 of 3,582 accessions reached a terminal status. 797 human single-matrix MS/NMR studies reached the independence test.
  - No qualifying person-identifier column among inspected records (658 + 220). R2: 18 passed layer 1, and none had a design sentence.
  - L1 = 0. 183 studies label every sample row `blank`.
- **3.2 Workbench.** 0 of 104 studies have a subject/participant factor key; the most common keys are design factors. Disclose the detector limitation, which undercounts repeated measures.
- **3.3 Documentary audit.**
  - 0/40 STRICT; 95% upper bound 25/371.
  - Explicit dependence 12/40: 9 repeated sampling, 3 cell cultures.
  - 11/40 have no publication text; 29/40 have at least one view, and none states one sample per participant.
  - Assay-scoped counts: 7 of 13 reconcile with the deposited rows.
  - Duplicate deposits within the pool: MTBLS1122/1124/1129/1130.
  - The futility rule was met; this result does not depend on any judgement call (robustness analysis).
  - Source: `DOCREC_T1_RESULT_v1.0.md`.
- **3.4 A35.**
  - In 17 of 19 deposits, other fields identify the rows as participant samples (V2 pass), so the `experimental blank` label is wrong. It sits on every row in 18 of 19.
  - 3/19 verified; 8 failed (non-participant material or repeated sampling).
  - Linkage defect in MTBLS14421: in 30 of 53 rows, the measurement labels cross patient/control groups.
  - One study deposited twice (MTBLS14358/14376).
  - Source: `DOCREC_A35_RESULT_v1.0.md`.

### 4. Discussion
- **What follows for external validation.** Repository-based study-level validation is infeasible under an affirmative standard. Consortia with custodian-held identifiers (COMETS; Yu et al. 2019) and controlled access are the realistic routes.
- **What follows for re-users.** Treating rows as participants risks pseudo-replication; this links to Waldron et al.
- **Recommendations:**
  1. A mandatory subject identifier or subject–sample map at deposition, as dbGaP requires.
  2. Fix vocabulary defaults for sample type.
  3. Reporting guidelines should require a samples-per-participant statement. Check STROBE-MetEpi first **[OPEN]**.
  4. Repositories should flag related deposits.
- **Limitations:**
  - one repository for the randomised audit;
  - n = 40, though the design-based bound is the result;
  - E is not I;
  - AI-assisted extraction (§2.7);
  - frozen-code limitations disclosed in the records: the "DEMs" table leak, the STAR+METHODS gap, and the rule-8 count gap;
  - open-access dependence (11/40 without text).

### 5. Figures and tables
- **Fig. 1** Flow diagram: 3,582 → 797 → L1 = 0; 390 candidates → T1 40 (grades) and A35 19 (V-status).
- **Fig. 2** T1 grades by reason (stacked bar).
- **Fig. 3** Hypergeometric upper bound on the recoverable pool as a function of x, with x = 0 → 25 marked.
- **Table 1** Prior work against this study.
- **Table 2** Types of documented dependence, with counts and example designs (no participant data).
- **Table 3** Metadata defects found: blank-label default, linkage mismatch, duplicate deposits, count mismatches.

## 6. Required before drafting or submission

1. **[OPEN] STROBE-MetEpi checklist** (Kelly et al., Research Square 2026, revision requested by *Metabolomics*). Does it require samples-per-participant? This decides recommendation 3 and part of the framing.
2. **Independent human check of the AI extraction (proposed; not yet done).**
   - The author re-extracts Q4 (samples per participant) and Q3 scope for a random 10 of the 29 T1 studies with publication text. He works blind to the AI answers, with the subset drawn by a recorded seed.
   - Report agreement.
   - Reviewers are likely to ask for this. It cannot change x = 0, which rests on an exhaustive phrase search, but it tests the other categories.
3. **Citation verification.**
   - The Karaman et al. preprint DOI prefix (`10.64898/…`) looks unusual. Verify it or drop the citation.
   - Confirm El Abiead et al. 2025 and Gonçalves & Musen 2019.
4. **Wording of the Workbench finding:** the REST `factors` endpoint only; not "Workbench has no subject IDs".
5. **Data availability.** Publish the records, evidence JSONL, briefs, validators and code. Do not publish `store/` (copyrighted full texts) or correspondence.
6. **Authorship and AI statement** per Springer policy: the author is listed as Anas Enoch (UM6SS), and AI is not an author.

## Not in this paper

- The selector programme.
- P4V.
- Any claim about true cohort counts.
