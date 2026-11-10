# NHANES fatty-acid G1 stability v0.4 — closeout addendum

10 November 2026. Revision 1 incorporates the post-hoc margin and
block-attribution checks. It supersedes the wording of the original
v0.4 addendum and preserves the development-only decision. It supplements
`NHANES_FA_STABILITY_CLOSEOUT_v0.4.md`. It uses the saved aggregate fit,
prediction and screen files. The attached independent checker reproduced
all 2,402 prediction sets from the integer histograms and reproduced the
four screen decisions. Its output matched
`nhanes_fa_stability_independent_check_result.json` byte for byte.

## Sharpness margin and disposition

The cap declared in the uncommitted v0.4 draft was **U ≤ 5**. The pooled
count cannot equal 5, so passing required at least 1,950 of 2,000 Gaussian
cycles with C ≤ 4. At the development size n = 1,725 there were **1,904**
such cycles (95.20%), 46 short. The highest count across the 1,201 sizes
was **1,946 at n = 2,273** (97.30%). That count is the maximum selected
from 1,201 independently seeded Monte Carlo counts.

A post-hoc binomial GLM (logit link, quadratic in sample size) fitted to
all 1,201 saved counts gives a larger shortfall at the selected size.
Pearson dispersion was 1.025, consistent with binomial Monte Carlo scatter
about the fitted trend. The following are descriptive estimates under that
trend. The fresh-run probabilities use a continuity-corrected normal
approximation to the binomial.

| Size n | Saved count | Fitted count | Fitted shortfall from 1,950 | Approximate chance of reaching 1,950 in a fresh run |
|---:|---:|---:|---:|---:|
| 1,725 | 1,904 | 1,883.1 | 66.9 (6.38 simulation SE) | 1.2 × 10⁻¹⁰ |
| 2,273 | 1,946 | 1,923.0 | 27.0 (3.14 simulation SE) | 0.10% |
| 2,400 | 1,940 | 1,929.9 | 20.1 (2.44 simulation SE) | 0.86% |

Under this fitted trend and the independent-seed calculation, the
approximate chance of at least one Monte Carlo crossing across all 1,201
sizes was **41%**. A rule that accepts any passing size is sensitive to
Monte Carlo fluctuation across the grid. In the observed run no crossing
occurred. The decision uses the saved counts and the original cap.

All 1,201 Gaussian upper limits were 6, and **no size passed**. The cap
stands. Milestone M closes at this development decision. The prediction
protocol never advanced to a freeze or to external calibration.

## Where the fitted models differ

At n = 1,725 the full-fit relative gap is (leader score − runner-up score) /
leader score. A reversal is a negative signed gap for the same
development-fixed pair in a simulated fold. The proportions are over
10,000 simulated folds per block per model (2,000 cycles × five folds);
folds within a cycle share training participants.

| Block | Development leader / runner-up | Relative full-fit gap | Gaussian reversal | Empirical resampling reversal |
|---:|---|---:|---:|---:|
| 1 | LBXTSA / LBXAR1 | 23.54% | 0% | 0% |
| 2 | LBXEPA / LBXGLA | 11.28% | 6.30% | 20.37% |
| 3 | LBXDTA / LBXMR1 | 20.80% | 0.04% | 3.57% |
| 4 | LBXEPA / LBXARA | 50.54% | 0% | 0% |
| 5 | LBXDTA / LBXARA | 25.46% | 0% | 0.05% |

Block 2 (EPA versus GLA) had the largest reversal rate under both models
and the largest difference between them. The observed 2011–12 refits had
c_b = 0 in each block, giving pooled C = 0. The stability endpoint was
chosen after that unanimity was observed.

The supplied n = 1,725 reproduction matched the saved histogram and
reversal counts for each model. The local empirical check used the
participant rows on the owner's machine. Its aggregate cycle patterns
give the following attribution.

| Model | Cycles with C > 0 | Block 2 only | Block 3 only | Block 5 only | Two or more blocks |
|---|---:|---:|---:|---:|---:|
| Gaussian | 272/2,000 (13.60%) | 269 | 3 | 0 | 0 |
| Empirical resampling | 746/2,000 (37.30%) | 545 | 137 | 1 | 63 |

The local v1.1 reproduction distinguishes unanimity from agreement
with the development picks. Each model has 2,000 cycles at n = 1,725.

| Cycle outcome | Gaussian | Empirical resampling |
|---|---:|---:|
| C = 0 | 1,728 (86.4%) | 1,254 (62.7%) |
| All 25 picks match the development picks | 1,700 (85.0%) | 1,162 (58.1%) |
| C = 0 with a different unanimous pick in at least one block | 28 (1.4%) | 92 (4.6%) |
| Block 2 unanimous on GLA | 28 | 98 |
| Block 3 unanimous on MR1 | 0 | 6 |

At least one pick differed from development in 15.0% of Gaussian cycles
and 41.9% of empirical cycles. Within-cycle disagreement occurred in 13.6%
and 37.3%, respectively. The differences of 1.4 and 4.6 percentage points
are cycles with within-block unanimity that still changed a development
pick. The per-block counts of 98 and 6 are not a partition of the 92
cycles with C = 0. They can include cycles in which another block
disagreed.

## Post-hoc discrimination check

This descriptive check uses the two fitted models' saved histograms and
does not enter the screen. At n = 1,725, 63 of 2,000 empirical-model
cycles (3.15%) fell above the Gaussian set [0, 6]. The proportion was
0.80%–6.70% across the exact sizes. At the development size, the saved
Gaussian histogram put 100% of its draws inside that set and the empirical
histogram put 96.85% inside. A hypothetical upper endpoint of 4 would put
351 of 2,000 empirical cycles (17.55%) above it at n = 1,725. This
hypothetical set was not the Gaussian prediction set and is not an
estimate of external power. The comparison illustrates how rarely the
actual set would distinguish the two fitted models at the development
size, even though their simulated disagreement rates differed.

## Provenance and deviations

- Protocol v0.4 and runner v0.1.0 were **uncommitted and unfrozen**.
  The cap was declared in an uncommitted draft before simulation. The
  private run copy records 2026-11-10T00:44:34Z file modification times
  for both files. The protocol and runner SHA-256 prefixes are
  `3bbce925` and `24b5daf3`, respectively. File timestamps and the
  execution transcript document ordering. There is no pre-run commit
  or preregistration timestamp for these files.
- The execution transcript records a passing development self-test,
  reproduction of all 25 recorded G1 picks with C = 0, and a separate
  scale-equivalence check on the actual training sets. It also records
  `bench` before the prediction run. At n = 1,725, the reported rates
  were 11 ms per Gaussian cycle (7.6 estimated core-hours) and 9 ms per
  empirical cycle (6.2 estimated core-hours). Separate self-test and
  benchmark logs were not retained in the run folder.
- The empty `predict-v0.4.log` dates to the detached launch at
  00:49:25Z, which exited before producing any parts. The foreground
  run started at 00:49:58Z and completed prediction, assembly and
  screening at 04:50:59Z. Progress and completion are recorded in
  the foreground execution transcript and the saved files. No simulation
  parts were reused from the failed detached launch.
- The five empirical upper-limit corrections in the original closeout are
  a **runner v0.1.0 defect** in `quantile_set`. Summing floating-point
  proportions placed the CDF just below 0.975 at the exact 1,950-count
  boundary, and the original test fixture did not exercise that boundary.
  Integer arithmetic recovers the five corrected upper limits and
  reproduces every Gaussian prediction set. Saved outputs are preserved.
- The missing 100-cycle descriptive NOVA-S prediction is a second
  runner v0.1.0 defect. The runner declares B_S = 100 and its prediction
  command never calls the agreement simulation or writes its output.
  A protocol-compliant freeze would require this deliverable. The
  runner's freeze command does not enforce its presence. The failed
  Gaussian sharpness screen already prevented a freeze.
- The fitted margin and attribution are post hoc. The attached
  `nhanes_fa_stability_margin_check_result.json` (v1.1) has SHA-256
  `474bf476d5e16f4b2f2ce090cbd41e44d3952555e06e08012d82bbb846484942`.
  The owner's local `nhanes_fa_stability_margin_check_result_local.json`
  (v1.1, both models) has SHA-256
  `cadbe90773532842ae2602d5ad0e3764046c1e3cbf7b76ece6bad71a9cff21fc`.
  The shared results agree. At n = 1,725, both models reproduced their
  saved histograms and fixed-pair reversal counts. This revision uses
  those supplied checks.

| Saved aggregate record | SHA-256 |
|---|---|
| `nhanes_fa_prediction_v0.4.json` | `dc8f83c845fb1b8be4454f7c83e0723220d0d9bb404d3750bd378e5bd14501fa` |
| `nhanes_fa_stab_out/screen_v0.4.json` | `cbf4a14d29f179aec03a301d5f2521eaa80054c411389b28d0183253dc8bcf39` |
| `nhanes_fa_stab_out/dev_fit.json` | `e9d894aa2e0b549ce4c02712cbbbcb2fd679ffddbac85ce9539420125dc46669` |

Participant rows and identifiers remain local. NHANES 2013–14 remains
unopened. This addendum creates no freeze, authorisation or commit.
