# N2 failure-regime diagnostic: specification v0.1 (DRAFT, not frozen)

Status: draft for the owner's review. It becomes frozen when the owner records the SHA-256 of this file and of the
script written to it, before the script is run. Classification: **exploratory, descriptive**, on the frozen N2
replay outputs of the 186-record development archive. It trains no selector, optimises no threshold, reads no
evaluation row, and supports no superiority claim. Its outputs are hypotheses for an external test, nothing more.

## 1. Question

Can information available *before* a measurement is chosen identify the records in which the nonlinear value
estimator's departures from G1 (NOVA-V pick ≠ G1 pick) realise an acquisition advantage? The question is about
realised held-out gain conditional on the paired statistic t, not about t itself: a covariate that predicts t but
not gain-given-t adds nothing a selector could use.

## 2. Inputs (all frozen; hashes checked by the script before any computation)

| Input | Content | SHA-256 |
|---|---|---|
| `n2_blockfolds.csv` | 4,760 block-folds: record, block, fold, j_V, j_G1, t, held-out loss sums and cell counts for both arms | b7934523cfc124510fdb72c0469089eecd8b8f4f3676b68552eb28618957d3bf |
| N0 per-record outputs | O-only predictive edge of kernel ridge over log-linear ridge (wide grid), per record | (owner to supply path and hash) |
| P4D design-model manifests | block membership, feature names and classes for the 952 blocks | (owner to supply; same files Gate 3A's loader read) |
| design rows of the 186 records | for the skewness covariate only (no fitting) | archive as verified by `f1_setup.py` |

No file produced by Paper III/IV's evaluation half, Gate 2, or the Gate D pool is read.

## 3. Units and dependence

- Analysis unit for all inference: the **record** (186). Block-folds within a record are dependent (same subjects
  across folds, same panel across blocks) and are never treated as independent observations.
- Every interval: record-clustered bootstrap, 10,000 resamples, seed 20261102, percentile intervals.
- Subgroups are reported with their record count; no statement is made about a subgroup with fewer than 30 records.

## 4. Outcome

For block-folds with j_V ≠ j_G1 (3,062):
- **Primary outcome** D = mean per-cell paired loss difference L_G1 − L_jV on the held-out fold (the difference
  estimand; linear in aggregation, as the estimand audit of 20 Oct records). Aggregated to record level as the
  cell-weighted mean over that record's differing block-folds.
- **Sensitivity**: the primary's own normalisation, 1 − Σ L_jV / Σ L_G1 within record.
- The per-block-fold ratio E[1 − L_jV/L_G1] is **not** used (estimand audit, finding 4.3).

## 5. Covariates (fixed; no additions after the first run)

| Covariate | Definition | Source | Why |
|---|---|---|---|
| C1 `log_n` | log of design-row count of the record | manifest | GDOA-0: nonlinear excess shrank with n |
| C2 `panel_class` | lipid-dominated (≥ 50 % of the 60 block features carry a lipid-class name under the frozen class rule of the manifests) vs broad | manifest | Paper 2 records note coherence confounded with study type |
| C3 `skew_O` | mean over the 20 base features of the sample skewness of the log-rule-transformed design values (log rule as in Gate 3A §2) | design rows | NBRA pre-check: untransformed intensity skew as an alternative explanation |
| C4 `edge_O` | N0's per-record O-only relative bounded-loss improvement of kernel ridge over log-linear ridge at the wide grid | N0 outputs | the one nonlinear signal the archive supports; measurable externally before acquisition |

Covariates are computed per record (C2 per block, then the record's majority class). All are computable in an
external cohort from its design partition alone.

## 6. Analyses (all pre-specified; all reported)

A1. **Descriptive table**: for each covariate, D and the sensitivity outcome by covariate tertile (C2: by class),
    with record counts and clustered intervals; separately for all differing block-folds and for departures (t ≥ 2).

A2. **Partial association (primary)**: record-level weighted least squares of D on t̄ (record mean t over differing
    block-folds) and C1–C4 jointly, weights = number of differing cells; coefficients with record-clustered
    bootstrap 95 % intervals; **Bonferroni over the four covariate coefficients** (97.5 % two-sided → 1.25 % each
    tail). The coefficient on t̄ is reported for context, not tested.

A3. **Does the covariate predict t or gain?** Same regression with t̄ as the outcome. A covariate that is
    associated with t̄ (A3) but not with D given t̄ (A2) is classified "predicts departure, not advantage".

A4. **Pre-stated subgroup check** for any covariate that passes A2: split at the covariate's design-row median
    (C2: by class); report the subgroup's mean D, its primary-normalisation equivalent, record count, and whether the
    subgroup's departures (t ≥ 2) give an equal-record mean relative improvement ≥ 0.29 % with clustered lower bound
    > 0. The median split is fixed here; no cut point is chosen from the data.

A5. **t ≥ 4 block-folds**: listed with their 8 records and the four covariates, as a table. No inference.

## 7. Pre-stated readings (chosen before the run)

- **"Regime identified (exploratory)"**: a covariate passes A2 at the Bonferroni level **and** its A4 subgroup has
  ≥ 30 records **and** meets the 0.29 % / LB > 0 condition. Consequence: the covariate defines a falsifiable
  hypothesis for an external cohort (the subgroup is predicted to show departures with positive realised gain; the
  complement is predicted not to). This is a hypothesis, not a finding.
- **"Predicts departure only"**: passes A3, fails A2. Consequence: the statistic t already carries the information;
  nothing for a selector to add.
- **"No regime identified"**: no covariate passes A2. Consequence: the heterogeneity in the replay is not explained
  by any pre-acquisition covariate in this list; the NOVA archive paper's negative framing stands and no new
  selector is proposed from this archive.
- A small-n-only association for C1 in the direction of larger apparent gain at smaller n is read as **estimation
  noise**, consistent with GDOA-0, not as a regime.

## 8. Integrity rules

- Script hash and this spec's hash recorded in `n2_regime_freeze_hashes.txt` before the run; the script refuses to
  run if `n2_blockfolds.csv` does not match b7934523….
- One run. A crashed run may be repeated with identical settings. Outputs: `n2_regime_summary.json` (aggregates and
  intervals), `n2_regime_records.csv` (per-record covariates and outcomes; kept with the development data).
- No covariate, cut point, outcome or test is added after the first run. Anything learned here that motivates a
  further analysis of the archive is written as a separate spec and labelled as a further use of the archive.
- The result record states the archive use count (this would be use 11) and that every quantity is exploratory.

## 9. What this diagnostic cannot do

It cannot show that any selector beats G1; it cannot confirm a regime; it cannot be used to choose a threshold or
an α for a hybrid. Its only legitimate output is a covariate-defined hypothesis, with its predicted sign, that an
external cohort with established participant structure could falsify.
