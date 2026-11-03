# N2 failure-regime diagnostic: specification v0.3 (for freezing)

Status: v0.3 of 2026-11-03, superseding v0.2 after inspection of the frozen file formats (changes in §11). Script: `n2_regime.py` v0.1.0. It is
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
| `gate3a_blocks.jsonl` | Gate 3A per-block output: per family the held-out losses, and per block `_g1_full_design_agreement` (folds whose refit G1 pick equals the full-design pick) and `_g1_folds` | 6decb9fb0b9f568571a5b24f12c9dad88ef9f31be9653681e0877b2bf36e3747 |
| `n0_records.csv` | N0 regularization parity, per record: bounded-loss sums for P1 and P3 at the frozen and wide grids on O alone | 6d3b0341b5fed7f81e6b0935ced5b0a0c7c8d0398b2e450a52258c11968b75ff |
| `paper4_design_models/<record>/design_model_manifest.json` (186 of the 188 directories) | `n_design`; per block the O, T and P feature names and the design-only per-feature median, mean and SD after design-median imputation (P4D, 29 Sep 2026) | recorded per record in `n2_regime_records.csv` at run time |

**No canonical matrix, design row or evaluation row is read.** Nothing produced by Paper III/IV's evaluation half,
by Gate 2, or by the Gate D pool is read. The script verifies the three pinned hashes and refuses to run on
mismatch; it also refuses to run unless it reproduces the frozen N2 primary (+0.1118 %) from the block-fold file.

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
| C2 `panel_class` | 1 if lipid-dominated: a block is lipid if ≥ 50 % of its 60 feature names (O ∪ T ∪ P) match the lexical lipid rule below; a record is lipid-dominated if ≥ 50 % of its blocks are | manifests (feature names) | no direction prespecified; two-sided |
| C3 `skew_O` | raw-scale skewness proxy: mean over the record's blocks of the mean over the 20 base features of 3·(mean − median)/SD, from the manifest's design-only per-feature statistics (Pearson's second skewness coefficient; no row is read) | manifests (design statistics) | no direction prespecified; two-sided |
| C4 `edge_O` | N0's per-record O-only relative bounded-loss improvement of kernel ridge over log-linear ridge (wide grid) | N0 outputs | H-nonlinear predicts **positive** |
| C5 `g1_stability` | Σ_blocks `_g1_full_design_agreement` / Σ_blocks `_g1_folds` for the record (0–1) | `gate3a_blocks.jsonl` | H-stability predicts **negative**: lower stability, larger realised gain conditional on t |

**Lexical lipid rule for C2** (fixed here; no claim of RefMet equivalence; disclosed as a new rule): a feature name
is lipid if it begins with one of `PC, PE, PI, PS, PG, PA, LPC, LPE, LPI, LPS, TG, TAG, DG, DAG, MG, MAG, CE, Cer,
HexCer, Hex2Cer, LacCer, GlcCer, SM, FA, CAR, CL, PIP, PIP2, PIP3` followed by a space, bracket, digit or hyphen
(case-insensitive), or contains any of `phosphatidyl, ceramide, sphingomyelin, sphingosine, sphinganine, triglycer,
triacylglycer, diglycer, diacylglycer, monoacylglycer, cholesteryl, cholesterol ester, lyso, acylcarnitine,
cardiolipin, ganglioside`, a sphingoid-base pattern `(d16:`/`(d18:`/`(t18:`, or a fatty-acyl tail `/NN:N)`. The
script's self-test checks the rule on 8 positive and 10 negative names (e.g. `Carnitine`, `Fatty acid binding
protein` are negative).

All five are computable in an external cohort from its design partition alone, before any acquisition; C3 and C5
need only summary statistics and a fold-refit of G1, never outcomes.

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

A4′. **Predictive, not ranking, instability (the owner's distinction).** The loss of the full-design G1 pick on the
    same fold is not stored in any frozen output, so the loss cost of G1's instability cannot be measured without
    refitting, which this spec forbids. Instead, for each block the **G1-plausible set** is the set of candidates
    G1 picked on any of the block's five folds (`jG1` in `n2_blockfolds.csv`). Each differing block-fold is
    classified *inside* (jV ∈ G1-plausible set: NOVA chose a candidate G1 itself chooses on a resample, i.e. a
    near-tie resolved) or *outside* (jV ∉ set: a candidate the covariance criterion never selects). For each side:
    D over all differing block-folds and over departures (t ≥ 2) with study-cluster 99 % intervals, and the
    contribution of that side's departures to the primary (the equal-record mean Δ_acq that the frozen NOVA-S
    pipeline yields when departures are allowed only on that side; the two contributions and the full primary are
    reported together). `Vhat_jV`/`Vhat_jG1` are not used: they are NOVA's values, not G1's.

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
  records and meets 0.29 % with LB > 0, **and** A4′ shows that inside-set departures realise positive gain (99 %
  interval of D for inside departures excluding 0 on the positive side). Hypothesis for an external test: a stability-gated policy (abstain from, or average over, G1's pick
  where its fold-stability is low), which is a different mechanism from nonlinear information extraction.
- **"Stability without cost"** — C5 passes A2 but the inside-set departures' interval in A4′ includes 0: the gain
  arises where G1's pick is a near-tie among candidates of nearly equal value. Consequence: no policy is implied; the finding is that NOVA-S's
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
- Self-tests (`n2_regime.py selftest`, synthetic data only), all required to pass before the real run: (i) planted
  C5 (negative) and C4 (positive) associations recovered with 99 % intervals excluding 0; (ii) a null with no
  association yields at most 3 of 100 covariate tests excluding 0 at the 99 % level across 20 replicates; (iii) the
  cluster bootstrap never splits a study's records; (iv) the NOVA-S primary reconstruction is consistent between
  the two code paths; (v) the lexical lipid rule classifies 8 positive and 10 negative names correctly.
- Reproduction gate: the run aborts unless the equal-record mean Δ_acq reconstructed from `n2_blockfolds.csv`
  equals the frozen N2 primary (+0.1118 %) to 1e-6.
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

## 11. Changes from v0.2 (after inspection of the frozen file formats, 3 Nov 2026)

1. A4 replaced by A4′: the full-design G1 pick's held-out loss is not stored in `gate3a_blocks.jsonl` (per-block
   output with `hat`, `g1`, `leak`, `O` losses and agreement counts only), so the loss cost of G1 instability is not
   computable without refitting; the inside/outside G1-plausible-set split from `n2_blockfolds.csv` (`jG1` per
   fold) is used instead and the "stability" and "stability without cost" readings are restated on it.
2. C2: a lexical lipid rule is specified, because the manifests carry feature names and no class, and B0's class
   rule needs live RefMet lookups.
3. C3: computed from the manifests' design-only per-feature median/mean/SD (raw-scale Pearson second skewness)
   instead of from design rows; the diagnostic now reads no row of any matrix.
4. C5: definition restated on the actual fields (`_g1_full_design_agreement`, `_g1_folds`).
5. Inputs table completed with the N0 file and its hash; reproduction gate and self-test list made concrete.
