# From prediction gain to measurement choice in budget-one metabolite reconstruction

Manuscript outline v0.3 · 10 November 2026

**Paper question.** When does a better predictive model or a different covariance
estimator improve the choice of one additional metabolite to measure for a
declared reconstruction task?

**Central finding.** At the operative shrinkage settings, frozen G1 changed
the Gaussian population winner in **24 of 36 NOVA-D conditions with a 2%
top-two gap** and **3 of 36 with a 20% gap**, all three at n = 128. At
G1-R's designated 2%-gap condition D1, decoupling regularisation raised
pick accuracy from **17.8% to 56.8%**, about threefold, while reducing mean
regret by **0.124 percentage points of no-acquisition risk**. The two
rules retained 98.3% and 99.1% of the oracle gain, respectively. G1-R's
designated no-harm criterion remained unestablished. The contrast between
pick accuracy and risk is the result the paper needs to explain.

## Abstract plan

1. Lead with ranking and risk on their respective scales. Shrinkage changed
   the population Gaussian winner in 24/36 conditions with a 2% gap and
   3/36 with a 20% gap, with the latter all at n = 128. At G1-R's
   designated 2%-gap condition, pick accuracy rose from 17.8% to 56.8%
   and mean regret fell by 0.124 percentage points of no-acquisition risk.
   The designated no-harm criterion remained unestablished.
2. Motivate the budget-one squared-loss reconstruction task. Give the signed
   relationship between affine and Bayes acquisition value and the matched
   Schur representation of joint covariance shrinkage as explanations of
   what each comparison measures.
3. State that NOVA-D missed both designated recovery criteria (0/500 NOVA-S
   departures in each). In both conditions the shrunk comparator evaluated
   at the population covariance already selected the nonlinear optimum,
   diminishing the opportunity against that comparator.
4. Report NHANES as development evidence. Each 2011–12 block had
   c_b = 0 in the observed refits, giving pooled C = 0. At the same size,
   the simulated probability
   of any disagreement was 13.6% under the Gaussian fit and 37.3% under
   empirical resampling. The largest model difference was in the EPA–GLA
   pair. The development sharpness screen stopped the proposed external
   evaluation before the 2013–14 files were opened.
5. State the limits of the paper. Evidence covers this reconstruction
   objective, the frozen simulation families, and the reported development
   checks. External calibration and external superiority remain unevaluated.

## 1. Introduction and estimands

- Start with the decision. Given observed panel O, target vector T and candidate
  analytes X_j, choose one candidate before measuring it. Keep loss, target
  scaling, budget and comparator explicit.
- Define three objects separately. Population acquisition value under a declared
  loss; oracle regret of a fitted selector; and held-out loss of the fitted
  prediction procedure. Refit disagreement is a fourth quantity that can be
  measured without evaluating prediction loss. It describes whether actions
  change across training sets.
- Explain the evidence sequence. NOVA supplies the estimand distinction, NOVA-D
  tests recovery, exact population checks diagnose a comparator ranking effect,
  G1-R tests finite-sample regret, and NHANES supplies development evidence.
  Cite Paper IV for objective dependence and the documentation audit for
  verification of sampling units. Leave their results and proofs in their own
  manuscripts.

## 2. Mathematical results

### 2.1 Affine versus Bayes acquisition value

Define A(S) = R_affine(S) − R_Bayes(S) under square-integrable targets and squared
loss. For candidate j the difference in population values is
V_Bayes,j − V_affine,j = A(O ∪ {X_j}) − A(O). Give the short projection argument.
The right-hand side is signed even though each approximation gap is
nonnegative. Include the recorded nonlinear proxy example where the added
analyte has affine predictive value and zero Bayes acquisition value.
Under joint Gaussianity and positive residual candidate variance, both values
reduce to the population Schur-complement score. Use V_affine for the
unregularised population quantity and G1/A0 for the frozen implementation
throughout the manuscript.

Source NOVA framework and derivation (`NOVA_PAPER_v0.5.md`, line 84).

### 2.2 What joint shrinkage changes

For unit-diagonal covariance, 0 ≤ α < 1 and observed-block ridge ε, give the
exact matched Schur representation. It uses conditioning ridge
(α + ε)/(1 − α), candidate variance floor α/(1 − α), a common score factor
1 − α, and the corresponding denominator clip rescaling. Explain the
difference between ranking distortion at the population covariance and the
finite-sample regret of the fitted decision. Reserve the operational ranking
counts for §3 so the algebra and the experiment stay distinct.

Source Equivalence derivation (`g1_shrinkage_equivalence.py`).

### 2.3 A finite-sample limit for a specified comparison

Sketch the residual-Wishart representation under iid Gaussian sampling,
unshrunk OLS conditioning, full-rank design and a fixed target metric.
Conditioning on the two-candidate scatter gives the known-parameter pair law.
Keep the derivation and numerical implementation details in the supplement.
In the main text, give the variance-perturbation fixed-pair example with free
candidate residual variances. At ν = 389 and level 0.05 the power ceiling is
approximately 7.42% for the stated D1-like alternative. This result concerns a
uniformly valid fixed-pair test in the specified model. It leaves selected-pair
and policy-level comparisons open. If retained, 2/√ν belongs in the supplement
as a planning scale for independent residual log-variance ratios.

Sources Residual-Wishart and pair law (`EXACT_WISHART_PAIR_NOTE_v0.1.md`) and Fixed-pair power ceiling (`CERTIFICATE_POWER_SCOPE_v0.1.md`, line 53).

## 3. Frozen synthetic studies

### 3.1 NOVA-D. Can the tested selector recover a planted advantage?

- Introduce the 72-condition design, with 100 replicates per condition and
  five folds per replicate. Distinguish the nonlinear population optimum, the
  Gaussian population comparator used in calibration, and the fitted frozen G1.
- Lead with the prespecified primary and positive control. NOVA-S departed
  from G1 in 0/500 fold decisions in each condition and recovered κ = 0 in
  both. Fold decisions within a replicate share training data. The null
  validity controls passed, and no near-tie threshold yielded a qualifying
  improvement.
- Separate the diagnostics from the decisions. The gate statistic stayed below
  its threshold, and ungated NOVA-V picked the nonlinear optimum less often than
  G1 in both designated conditions. Exact population checks later found a
  changed Gaussian argmax in 27/72 conditions at operative shrinkage and ridge.
  Report the stratified denominators of 24/36 at a 2% gap and 3/36 at a
  20% gap, with the latter all at n = 128. In both designated conditions the shrunk population comparator picked the
  nonlinear optimum. This qualifies the opportunity against the actual
  comparator. It does not identify the cause of individual finite-sample errors.

Sources Corrected NOVA-D record (`NOVA_D_RESULT_v1.2.md`, line 51) and Exact population ordering (`NOVA_D_EXACT_ORDERING_NOTE_v1.2.md`, line 39).

### 3.2 G1-R. Does decoupling regularisation improve decisions under both criteria?

- Introduce the 360-condition Gaussian simulation, 1,000 replicates per
  condition, with exact oracle regret for each selected candidate. Explain why
  this endpoint isolates decision quality from downstream predictor error.
- Show both designated conditions together. At D1 (2% gap), pick
  accuracy rose from 0.178 to 0.568. Mean normalised regret fell from
  0.002574 to 0.001332, a gain of 0.001242 or 0.1242 percentage points of
  no-acquisition risk. At oracle amplitude A★ = 0.15, the two rules kept
  98.3% and 99.1% of oracle gain. D1 met the prespecified improvement
  criterion. D2 did not establish the 10% no-harm criterion. Its Q2 lower
  bound was −0.000408 and its two-sided interval crossed zero
  [−0.000456, +0.000102]. The joint criterion was unmet. Harm beyond the
  margin was not established.
- Put the gain in context using the charter's analogous matched-loss bar
  as a descriptive scale, separate from G1-R's decision rule. With
  A_G1 = 0.15 − 0.002574, the bar 0.10 A_G1(1 − A_G1) is about 0.01257
  of no-acquisition risk. On the recorded mean-regret scale, even oracle
  selection would recover only 20.5% of that bar, and A2 recovered 9.9%.
  This headroom calculation holds the target, candidates, budget and
  Bayes loss fixed. It does not set an external relevance threshold for D1.
- Attribute the D1 population reordering to joint shrinkage at the operative
  setting. The 0.05 floor was inactive there. Report grid summaries and sparse
  contrasts descriptively. Put the selected grid maximum, if included at all,
  in a labelled supplement.

Source G1-R decision and designated conditions (`G1R_RESULT_v1.1.md`, line 22).

## 4. Development evidence and the decision to stop

### 4.1 Archive context and observed NHANES picks

- Use the versioned Paper III REML correction for the Gaussian-versus-variance
  benchmark. The estimate was +0.00867 [0.00659, 0.01075] across 186 records,
  below the 0.014 practical RMSE margin. The separate NOVA-S development
  contrast was +0.112% [+0.029, +0.221] against its 0.29% qualification scale.
  Keep loss units, reuse of records, and sampling-unit limitations explicit.
- In Gate 3A, the 29% is disagreement in a particular refit check (3,376 of
  4,760 block-folds agreed). It is distinct from the within-block all-pairs
  NHANES statistic C/50. Present them separately with their own denominators.
- In NHANES 2011–12, 1,725 eligible participants contributed to five blocks
  sharing five participant-disjoint held-out folds. The training sets
  overlap. G1 chose the same candidate within each block (observed C = 0).
  NOVA-S matched G1 in 25/25 block-folds. NOVA-V and NOVA-U each matched
  G1 in 24/25. The stability endpoint was chosen after observing C = 0.

### 4.2 What the development-fitted models predict

Describe the Gaussian and participant-resampling models, fitted to the same
2011–12 data, and the whole-cycle simulations that preserve the shared folds
and analytes. Report the 1,201 exact sizes from 1,200 to 2,400, with 2,000 cycles
per size per model. Treat C = 0 as an observation from this cycle. At
n = 1,725, 13.6% of Gaussian cycles and 37.3% of resampled cycles showed
within-cycle disagreement. All 25 picks matched the development picks
in 85.0% and 58.1% of cycles, respectively.

**Main-text Table 2.** At n = 1,725, reversal means a negative signed
difference for the pair fixed on the full development fit. Relative gap is
(leader − runner-up)/leader. Report Monte Carlo proportions per simulated
fold. Report the cycle-level attribution separately from these marginal
fixed-pair reversal proportions.

| Block | Development leader / runner-up | Full-fit gap | Gaussian reversal | Empirical resampling reversal |
|---:|---|---:|---:|---:|
| 1 | LBXTSA / LBXAR1 | 23.54% | 0% | 0% |
| 2 | LBXEPA / LBXGLA | 11.28% | 6.30% | 20.37% |
| 3 | LBXDTA / LBXMR1 | 20.80% | 0.04% | 3.57% |
| 4 | LBXEPA / LBXARA | 50.54% | 0% | 0% |
| 5 | LBXDTA / LBXARA | 25.46% | 0% | 0.05% |

Lead with EPA–GLA. This pair had the smallest full-fit relative gap and the
largest difference in reversal proportions. The supplied reproduction at
n = 1,725 matched the saved histograms and reversal counts. Among 272
Gaussian cycles with C > 0, 269 had disagreement in block 2 only and 3
in block 3 only. Among 746 empirical cycles with C > 0, the single-block
counts were 545 for block 2, 137 for block 3, and one for block 5, with 63
cycles showing disagreement in multiple blocks. In the local v1.1
reproduction, 28 Gaussian cycles (1.4%) and 92 empirical cycles (4.6%)
had C = 0 with a different unanimous pick somewhere. Use these counts to
distinguish stability within a cycle from recovery of the development picks.
The fitted-model differences do not identify the panel's true sampling
behaviour.

### 4.3 Why the external test did not run

- The cap U ≤ 5 was declared in an uncommitted draft before simulation. It
  required 1,950/2,000 Gaussian cycles at C ≤ 4. The development-size count
  was 1,904. The largest saved count, 1,946 at n = 2,273, was selected from
  1,201 Monte Carlo counts. A post-hoc quadratic-logit binomial fit gives
  expected shortfalls of about 67 cycles at n = 1,725, 27 at n = 2,273 and
  20 at n = 2,400. Under the fitted trend the approximate chance of any
  Monte Carlo crossing across the size grid was 41%. This diagnoses a
  weakness in screening many noisy estimates. The saved Gaussian upper
  limit was 6 at every size, so zero sizes passed the unchanged cap.
- Place the flat upper-limit plot and the two n = 1,725 histograms in the
  supplement. Give the independent integer check and the five empirical
  quantile corrections there, naming the runner defect and preserving the
  original files. The Gaussian screen result stands.
- Label the histogram comparison as post hoc. At n = 1,725, only 3.15% of
  empirical-model cycles exceed the actual Gaussian set [0, 6]. Across
  sizes the proportion is 0.80%–6.70%. A hypothetical upper endpoint of 4
  raises the development-size proportion to 17.55%. This is a descriptive
  model comparison using the saved fits. It supplies no external power
  estimate and leaves the cap unchanged.
- The stability protocol did not advance to a freeze or external calibration.
  The earlier NOVA-S superiority NO-GO is separate. NHANES 2013–14 remains
  unopened. Runner v0.1.0 never implemented the declared 100-cycle NOVA-S
  prediction. This is a separate implementation defect that would have
  needed resolution for a protocol-compliant freeze.

Sources `PAPER3_REML_CORRECTION.md` (Results),
`SELECTOR_CLAIM_EVIDENCE_MATRIX_v0.4.md` (§§1.5 and 3),
`NHANES_FA_RESULT_v1.1.md` (versioned wording correction),
`NHANES_FA_STABILITY_PREDICTION_PROTOCOL_v0.4.md` (§§0–3),
and the two v0.4 closeout records and independent check listed below.

## 5. Discussion

- Lead with the two results that could have gone the other way. The
  reordering was concentrated in the 2%-gap conditions. At D1 a roughly
  threefold increase in pick accuracy moved risk by only 0.124 percentage
  points of no-acquisition risk. The designated no-harm criterion remained
  unresolved at D2.
- Explain what the mathematics adds. A nonlinear proxy can matter to an
  affine predictor without increasing Bayes value. The equivalence explains
  how shrinkage changes conditional scores. The fixed-pair power calculation
  applies to its stated model and does not supply a universal limit on
  selectors, policies, or real panels.
- Interpret the development evidence. Agreement between selectors and
  stability across refits are different measurements. Observed unanimity
  in one cycle coexists with disagreement under both fitted models.
  The large EPA–GLA discrepancy and low post-hoc model discrimination
  explain why a compatibility verdict would have offered little separation
  between these fits. Neither model was evaluated on 2013–14.
- Close the budget-one reconstruction question at the evidence it produced.
  A panel, platform, or reflex-assay decision would define a different study.
  Keep the documentation audit and Paper IV independent of this manuscript.

## Planned displays and supplement

- **Table 1.** NOVA-D and G1-R designated results, with the comparator,
  endpoint, rule, and denominator for each study. Stratify the post-hoc
  27/72 population reordering count by 2% and 20% gaps. Report D1 pick
  accuracy and risk side by side.
- **Table 2.** The five-block fixed-pair gap and reversal table in §4.2.
- **Supplement.** The flat Gaussian upper-limit plot with the unchanged cap,
  the n = 1,725 model histograms, proof details and assumptions, frozen
  study identifiers, descriptive grid summaries, and the five empirical
  quantile corrections. Label frozen analyses, development screens, and
  post-hoc summaries separately in captions.

## Source and release notes

Repository references in this outline use repository-relative filenames.
The development-only stability sources are
`NHANES_FA_STABILITY_CLOSEOUT_v0.4.md`,
`NHANES_FA_STABILITY_CLOSEOUT_ADDENDUM_v0.4_r1.md`, and
`nhanes_fa_stability_independent_check_result.json`. The independent check is
`nhanes_fa_stability_independent_check.py`. The closeout gives the aggregate
prediction and screen hashes. The addendum gives the fit hash, execution
provenance, and the identities of the attached and local margin checks.
The v1.1 local reproduction for both models is
`nhanes_fa_stability_margin_check_result_local.json`.

The Paper III REML correction and NHANES v1.1 wording correction are
archived in commit `52e3e9b`. Any public package of the stability
materials should be limited to the protocol, code, aggregate records and
errata. Exclude
`dev_rows.npy`, `dev_ids.json`, and participant-level source files. The audit paper's v0.8 commit and tag
have been published on GitHub. Journal submission is a separate step.
