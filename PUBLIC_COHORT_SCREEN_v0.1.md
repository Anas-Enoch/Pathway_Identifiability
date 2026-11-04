# PUBLIC_COHORT_SCREEN_v0.1

**Prospective amendment to EXTERNAL_VALIDATION_PREREG_STUB_v0.1.md §6**  
Date: 4 November 2026. Programme: frozen NOVA-S versus G1 acquisition, budget 1.  
**Status: proposed; NOT frozen until SHA-256 recorded and Git commit successfully pushed.**

## 1. Provenance and rationale

The original stub prohibits additional repository-search cycles. This document seeks one bounded exception to the *data-discovery restriction only*. It does not change the 186-record archive results, frozen selectors, hyperparameters, success criterion, or the preexisting LIPID/COMETS inquiries (15 December 2026 decision deadline). The amendment is motivated by lack of suitable immediately accessible external data. Before this amendment was written, two potential records were named: **NGDC OMIX009009** and a **Dryad aging dataset reportedly comprising 56 selected metabolite features**. Both are counted in the cap, disclosed as preidentified, and require a recorded disposition. No candidate-specific measurement values or NOVA/G1 results may be examined during discovery.

## 2. Candidate universe, reproducible order, and cap

**Cap: 20 distinct candidate records in one cycle**, inclusive of both preidentified records. Six source families: NGDC OMIX, Dryad, Zenodo, figshare, MetaboLights first released after **28 October 2026**, and Metabolomics Workbench first released after **17 October 2026** (B0 date). A source whose first-release date cannot be verified is not assumed new.

Two slots are fixed for the preidentified records. The remaining maximum 18 are **three discovery slots for each of the six families**, with no replacement between families. For each family, record the UTC enumeration timestamp, API/index endpoint, exact search expression, total hits and archived index response. Restrict indexing to records whose title/abstract/keywords mention metabolomics or lipidomics and public human participant cohorts. Enumerate candidate accessions without opening quantitative data, deduplicate by persistent ID/DOI and sort eligible index records by the lexicographic SHA-256 digest of UTF-8 `PUBLIC-COHORT-v0.1|<family>|<canonical-id>`. Screen the first three according to the fixed order. If a source cannot be enumerated reproducibly through public metadata, mark its three slots **unsearchable**, not manually substituted. Do not stop early upon finding promising records. Previously contacted candidates still occupy their slots.

## 3. Eight sequential metadata-only eligibility criteria

For each criterion record **PASS / FAIL / UNKNOWN**, a supporting URL/document reference, an exact relevant passage or field description, and whether the evidence applies to the *deposited assay matrix*. UNKNOWN does not qualify. Record every identifiable failure, not only the first.

**E1. Documented participant structure (binding).** An explicit person-to-sample mapping, or documentary evidence meeting the affirmative Stage A/T1 independence rules for the exact analysed assay, must establish sampling structure and independent units. A cohort enrollment count or unique sample labels alone are insufficient. Repeated specimens can qualify if a recoverable participant key permits participant-level clustering.

**E2. Compatible feature count.** At least 60 simultaneously measured and uniquely identifiable quantitative metabolites remain after the documented preprocessing exclusions, supporting 20 observed / 20 targets / 20 candidates. Do not infer this from assay names alone.

**E3. Participant-disjoint evaluation.** Document two nonoverlapping cohorts or a cohort that permits a frozen participant-disjoint evaluation partition. Label whether eventual evidence is cross-cohort transport or within-cohort holdout; these are not synonymous.

**E4. Accessible quantified matrix.** Obtainable processed metabolite matrix and sample mapping under lawful access terms, without requiring ad hoc raw-spectra reconstruction or unverifiable value-dependent processing.

**E5. Feasibility for the frozen task.** Metadata support the expected eligible participant counts, panel overlap and plausible precision. No numeric MDE may be inferred from an unobserved NOVA/G1 contrast; design-based power assessment and numerical GO/NO-GO thresholds are frozen separately before value access.

**E6. No prior programme use/contact.** Check accession, DOI, dataset title and aliases against existing Git history, Stage A/B0/T1, NHANES/Stage16 and NOVA development records. Previously value-opened or used-for-method-selection records are excluded. Document grep/search coverage and aliases.

**E7. Shared measurement platform.** Both populations use the same identifiable assay platform/family and demonstrably compatible analyte definitions and measurement units; a cross-platform harmonization requires a separate preregistered study.

**E8. No endpoint-driven feature selection that invalidates transfer.** The assay panel was not constructed by selecting features to predict a phenotype/outcome that contaminates the intended acquisition task. A 56-feature outcome-selected Dryad panel independently fails E2 and E8 if its initial description is verified.

## 4. Development feasibility and subsequent evaluation

The search itself cannot open quantitative measurements. A record meeting E1–E8 needs a separately hashed executable protocol committed *before* value access. Reuse frozen `p6_n2.py` SHA-256 `681ce852609e6ab4fcaf0f0ff603477c64e2351e10bd119811d45f008d4bbe8e`, its pinned dependencies, G1 and NOVA-S, τ=2, N2 nuisance grids, identical preprocessing and prediction loss. Do not tune a new selector.

For a qualifying panel of at least 60 metabolites, the prospective default is **five fixed salted 20/20/20 blocks** with participant-level pooled paired loss and sampling-cluster-respecting uncertainty; the five blocks and their folds are dependent, not independent study replications. Block salts, allocation rule and feasibility thresholds must be fully specified in the executable protocol before values are accessed. The original P4D panel specification has priority; this multi-block extension is declared rather than silently relabeled as identical to P4D. Do not transfer NHANES's restricted 10/10/9 assignment.

The development GO/NO-GO rule must include nondegeneracy: **zero NOVA-S departures, zero/nonfinite paired-loss SE, or failure of prespecified sample-size/power requirements ⇒ NO-GO**. The sign of the development acquisition contrast must not affect GO. An MDE ceiling (a candidate is 0.30 percentage points) must be justified for the cohort and fixed before development quantitative data are opened; 0.30 pp is not automatically transferable from NHANES. External evaluation data may be opened once only following a documented GO and executable freeze.

The original supportive external success rule remains: positive one-sided 97.5% lower bound of `Δ_ext = 1 − L_NOVA-S/L_G1` AND `Δ_ext ≥ 0.10 × Δ_G1→∅` under the same cohort, tasks and loss normalization. Record the G1 base-panel benefit, selector agreement, departure rate and uncertainty regardless of outcome. One independent cohort supports, but does not establish, generalizable superiority.

## 5. Full screen log and stopping rule

For every selected candidate log sequential ID, discovery timestamp, family, persistent ID/DOI, title, source URLs, preidentified yes/no, source-index snapshot, all E1–E8 findings and citations, first failure, final PASS/FAIL/UNKNOWN, aliases checked and any deviations. No deleting unfavourable records, refilling slots, reranking based on outcome or opening values during screening. If no candidate qualifies among 20, the cycle ends **NO ELIGIBLE COHORT**. Another search requires a new documented prospective amendment, not a quiet extension.

## 6. Preidentified records — entries to evaluate after freeze

| Entry | Candidate | Status at freeze | Prior exposure | Provisional issue |
|---|---|---|---|---|
| 01 | NGDC **OMIX009009** | Not screened | Mentioned before amendment | A reported 704 participants in two cohorts does **not** establish assay-specific participant linkage. E1, E2–E8 require documentary checking. |
| 02 | Dryad aging/metabolomics dataset (persistent DOI to be resolved from the already mentioned URL) | Not screened | Mentioned before amendment | Reportedly 56 outcome-selected features: provisional **E2 fail**, and evaluate E8; confirm against official metadata. |

**Scientific label:** This is a disclosed, prospective *data-access strategy amendment*. It does not rescue or retroactively validate prior negative results. Do not represent a subsequently identified cohort as wholly preregistered under the original, unamended stub.
