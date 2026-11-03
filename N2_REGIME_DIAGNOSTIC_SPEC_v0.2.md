# N2 failure-regime diagnostic: specification v0.2 (for freezing)

Status: v0.2 of 2026-11-03, issued on the owner's approval of v0.1 with six requirements (recorded in §10). It is
frozen when the owner records the SHA-256 of this file and of `n2_regime.py` in `n2_regime_freeze_hashes.txt` before
the script is run. Classification: **exploratory, descriptive, one run**, on frozen N2 and Gate 3A outputs of the
186-record development archive (archive use 11). It trains no selector, refits no model, optimises no threshold,
reads no evaluation row, adds no covariate after the first run, and supports no superiority claim. Its only
legitimate output is a covariate-defined hypothesis, with a predicted sign, for an external cohort to falsify.

## 1. Question

The programme has one development result in which departures from G1 realised a positive held-out gain (NOVA-S,
+0.112 %, N2), and two candidate explanations for it:

- **H-nonlinear**: the kernel estimator extracts conditional information the covariance criterion cannot
  represent, and the gain concentrates where that information is strong.
- **H-stability**: G1's own estimated ranking is unstable at these sample sizes (Gate 3A: the G1 pick re-estimated
  on 4/5 of the design rows equals the full-design pick in 3,376 of 4,760 block-folds, 70.9 %), and the gain
  concentrates where G1's pick is a near-tie that a second estimator happens to resolve better.

The diagnostic asks which, if either, is supported by pre-acquisition covariates **predicting realised held-out
gain conditional on the paired statistic t**. A covariate that predicts t but not gain-given-t adds nothing a
selector could use. The owner's distinction is adopted as a design constraint: G1 *selection* instability does not
imply *predictive* instability, since differently selected candidates may have nearly identical value; therefore
every test in §6 has realised gain as its outcome, and C5 is never tested against pick disagreement alone.

## 2. Inputs (frozen; hashes checked by the script before any computation)

| Input | Content | SHA-256 |
|---|---|---|
| `n2_blockfolds.csv` | 4,760 block-folds: record, block, outer fold, j_V, j_G1, t, held-out loss sums and cell counts for both arms (N2 replay) | b7934523cfc124510fdb72c0469089eecd8b8f4f3676b68552eb28618957d3bf |
| `gate3a_blocks.jsonl` | Gate 3A per block-fold output, including the fold-refit G1 pick and the full-design G1 pick | 6decb9fb0b9f568571a5b24f12c9dad88ef9f31be9653681e0877b2bf36e3747 |
| N0 per-record outputs | O-only relative bounded-loss improvement of kernel ridge over log-linear ridge at the wide grid, per record | (owner supplies path and hash before freezing) |
| P4D design-model manifests | block membership, feature names and classes, design-row counts, for the 952 blocks | (owner supplies path and hash before freezing) |
| design rows of the 186 records | for C3 only (no fitting; read by the frozen block loader of `uapa_f1.py`, SHA-256 4fce6a88…) | archive as verified by `f1_setup.py` |

Nothing produced by Paper III/IV's evaluation half, by Gate 2, or by the Gate D pool is read. The script refuses to
run if any listed hash does not match.

## 3. Units and dependence

- Analysis unit: the **record** (186). Block-folds within a record are dependent and are never treated as
  independent observations.
- **Clustering unit for every bootstrap: the underlying study.** Records are analyses (`ST…_AN…`); where two or more
  records share a study accession, they form one cluster. The script reports the number of clusters (≤ 186) before
  any result. 10,000 cluster resamples, seed 20261103, percentile intervals.
- A subgroup is reported only with its record count and cluster count; nothing is stated about a subgroup with fewer
  than 30 records.

## 4. Outcome

For block-folds with j_V ≠ j_G1 (3,062):
- **Primary outcome D**: mean per-cell paired loss difference L_G1 − L_jV on the held-out fold (the difference
  estimand, linear in aggregation; estimand audit of 20 Oct, finding 4.1). Aggregated to the record as the
  cell-weighted mean over the record's differing block-folds.
- **Sensitivity**: the primary's own normalisation, 1 − Σ L_jV / Σ L_G1 within record.
- The per-block-fold ratio E[1 − L_jV/L_G1] is not used (estimand audit, finding 4.3).

## 5. Covariates (exactly five; none added after the first run)

| Covariate | Definition | Source | Prespecified direction (sign of the A2 coefficient under the hypothesis it serves) |
|---|---|---|---|
| C1 `log_n` | log of the record's design-row count | manifests | H-nonlinear predicts **positive** (estimable nonlinear value grows with n); a **negative** coefficient is read as small-n estimation noise, consistent with GDOA-0 |
| C2 `panel_class` | 1 if lipid-dominated (≥ 50 % of the record's 60 block features carry a lipid-class name under the manifests' frozen class rule, majority over blocks), else 0 | manifests | no direction prespecified; two-sided |
| C3 `skew_O` | mean over the 20 base features of the sample skewness of the log-rule-transformed design values (log rule of Gate 3A §2, fitted on all design rows of the record) | design rows | no direction prespecified; two-sided |
| C4 `edge_O` | N0's per-record O-only relative bounded-loss improvement of kernel ridge over log-linear ridge (wide grid) | N0 outputs | H-nonlinear predicts **positive** |
| C5 `g1_stability` | share of the record's block-folds in Gate 3A whose fold-refit G1 pick equals the full-design G1 pick (0–1) | `gate3a_blocks.jsonl` | H-stability predicts **negative**: lower stability, larger realised gain conditional on t |

All five are computable in an external cohort from its design partition alone, before any acquisition.

## 6. Analyses (all prespecified; all reported whatever they show)

A1. **Descriptive table**: D and the sensitivity outcome by covariate tertile (C2: by class), record and cluster
    counts, cluster-bootstrap intervals; for all differing block-folds and for departures (t ≥ 2) separately.

A2. **Partial association (primary)**: record-level weighted least squares of D on t̄ (the record's mean t over
    differing block-folds) and C1–C5 jointly, weights = number of differing cells, covariates standardised to unit
    SD. Coefficients with **two-sided 99 % cluster-bootstrap intervals (Bonferroni over five covariates)**. The
    coefficient on t̄ is reported for context and is not a test. Collinearity check: the condition number of the
    standardised design and all pairwise covariate correlations are reported before the coefficients.

A3. **Covariate → t or covariate → gain?** The same regression with t̄ as the outcome. A covariate associated with
    t̄ (A3) but not with D given t̄ (A2) is classified "predicts departure, not advantage".

A4. **Predictive, not ranking, instability (the owner's distinction).** For C5 only, in addition to A2: among the
    record's block-folds in which the *Gate 3A* fold-refit G1 pick differed from the full-design pick, the realised
    held-out loss of the two G1 picks is compared from `gate3a_blocks.jsonl` (both are G1 picks under the same P1
    family, so the contrast measures what G1's instability costs in loss, not in rank). Reported as the equal-record
    mean relative difference with a cluster-bootstrap interval. If this quantity is near zero, selection instability
    is not predictive instability, and a C5 effect in A2 would have to be read as "near-ties resolved", not as "G1
    errors corrected".

A5. **Prespecified subgroup check** for any covariate that passes A2: split at the covariate's median over the 186
    records (C2: by class). For the subgroup on the hypothesised side (C1 high, C4 high, C5 low; two-sided
    covariates: the side with the larger D), report record and cluster counts, mean D, the primary-normalisation
    equivalent, and whether the subgroup's departures (t ≥ 2) give an equal-record mean relative improvement
    ≥ 0.29 % with a cluster-bootstrap lower 99 % bound > 0. The median split is fixed here; no cut point is chosen
    from the data.

A6. **The t ≥ 4 block-folds** (13, from 8 records): listed with record, study, and C1–C5. No inference.

## 7. Prestated readings (fixed before the run)

- **"Regime identified (exploratory): nonlinear"** — C4 (or C1 with positive sign) passes A2 at the Bonferroni
  level, its A5 subgroup has ≥ 30 records, and the subgroup meets 0.29 % with LB > 0. Hypothesis for an external
  test: a nonlinear departure policy gated on C4 (or n).
- **"Regime identified (exploratory): stability"** — C5 passes A2 with negative sign, its A5 subgroup has ≥ 30
  records and meets 0.29 % with LB > 0, **and** A4 shows a non-trivial loss cost of G1's own instability (interval
  excluding 0). Hypothesis for an external test: a stability-gated policy (abstain from, or average over, G1's pick
  where its fold-stability is low), which is a different mechanism from nonlinear information extraction.
- **"Stability without cost"** — C5 passes A2 but A4's interval includes 0: the gain arises where G1's pick is a
  near-tie among candidates of nearly equal value. Consequence: no policy is implied; the finding is that NOVA-S's
  gain is partly an artefact of tie-breaking, and the external-validation hypothesis is neither nonlinear nor
  stability-based.
- **"Predicts departure only"** — a covariate passes A3 but not A2. Consequence: t already carries the information.
- **"No regime identified"** — no covariate passes A2. Consequence: the NOVA archive paper's negative framing
  stands; no new selector is proposed from this archive; the next external-validation hypothesis is neither.
- Two covariates passing A2 is reported as both, with their partial coefficients; no ranking between hypotheses is
  made from this archive.

## 8. Integrity rules

- `n2_regime.py` is written to this specification; its SHA-256 and this file's are recorded in
  `n2_regime_freeze_hashes.txt` before the run. The script verifies every input hash in §2 and refuses to run on
  mismatch, reads no file outside §2, and writes `n2_regime_summary.json` (aggregates, intervals, counts),
  `n2_regime_records.csv` (per-record covariates and outcomes; stays with the development data) and a run manifest.
- Self-tests before the real run, all required to pass: (i) a synthetic planted association between a covariate and
  D is recovered with interval excluding 0; (ii) a null with no association yields at most 1 of 5 coefficients
  excluding 0 across 20 replicates at the 99 % level (expected 1 in 100 tests); (iii) the cluster bootstrap
  reproduces the clustering (two records of one study always resampled together).
- One run. A crashed run may be repeated with identical settings. No covariate, cut point, outcome, model or
  test is added after the first run; a further analysis of the archive is a separate spec and a counted use.
- The result record states: archive use 11; exploratory; cluster count; every quantity in §6 whatever its value.

## 9. What this diagnostic cannot do

It cannot show that any selector beats G1; it cannot confirm a regime; it cannot be used to choose a threshold, an
α for a hybrid, or a stability cut-off for deployment. Its output is one of the readings in §7, each of which names
the hypothesis, if any, that an external cohort with established participant structure would test.

## 10. Changes from v0.1 (owner's requirements of 3 Nov 2026)

1. Five covariates, C5 `g1_stability` added from the frozen Gate 3A block output. 2. C5 direction prespecified
(negative). 3. Multiplicity: Bonferroni over five, two-sided 99 % intervals (v0.1 had four and 97.5 %). 4. Bootstrap
clustering at the study level where records share a study (v0.1 clustered at the record). 5. Regime criterion
restated with the 99 % lower bound. 6. Inputs restricted to frozen N2 and Gate 3A outputs plus the fixed metadata
needed for C1–C4; no refitting. 7. A4 added to separate selection instability from predictive instability, and the
"stability without cost" reading added, on the owner's distinction. 8. Seed changed to 20261103.
