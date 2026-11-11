# NHANES 2011–2012 value-transport feasibility v0.2

**Decision.** The frozen reconstruction procedure supports planning a test of retained value for fixed measurement choices. The matched affine diagnostic was less precise and could not close that route. This record is development-only. It neither freezes a protocol nor authorises use of the reserved cycle.

## E1. Frozen procedure reproduced

The 2011–2012 files, frozen wrapper, blocks, saved fit and result were verified by SHA-256. The rerun used the original participant folds, G1 choices, target scaling, RBF predictor, inner tuning and survey-weighted pooling. All 25 G1 choices matched the saved record. Both pooled gains and their Taylor standard errors matched the recorded squared and bounded results at relative tolerance 10⁻⁷. The rerun computed all nine candidate predictions per block for the new E2 and E3 diagnostics. No participant-level output was written.

| Prediction instrument and loss | Pooled gain (%) | Design SE (pp) | Pointwise 95% interval (%) |
|---|---:|---:|---:|
| Frozen RBF, squared | 15.77 | 1.25 | [13.13, 18.42] |
| Frozen RBF, bounded | 16.74 | 1.74 | [13.07, 20.40] |
| Matched affine, squared (v0.1) | 11.46 | 3.27 | [4.57, 18.35] |

In the affine diagnostic, EPA's standard errors were 6.20 and 8.29 pp in blocks 2 and 4. In the frozen squared-loss rerun, those standard errors were 3.59 and 4.70 pp. This compares two prediction instruments. It leaves the cause of the difference unresolved.

**Proposed primary rule for a separate prospective protocol.** Fix the five development choices, their target scales, the RBF reconstruction recipe, the survey-weighted pooled squared-loss gain and the floor at **11.83%**. This is 75% of the frozen 2011–2012 point estimate. Use participant-disjoint refits of the frozen predictor recipe within the evaluation cycle while holding the measurement choices and 2011–2012 target scales fixed. This tests retention of measurement value with the same reconstruction recipe. Transporting model coefficients trained in 2011–2012 would require separate development validation.

Call the fixed benchmark exceeded when the one-sided 97.5% lower confidence bound is greater than 11.83%. If the corresponding upper bound is below the floor, the benchmark is missed. If the interval crosses the floor, retention remains unresolved. A comparison with the unknown true development gain would require propagating uncertainty from both cycles. This proposal uses the recorded point estimate as a fixed benchmark.

The noncentral t planning calculation assumes 17 design degrees of freedom and a true evaluation gain of 15.77%. Power is 84.2% at the development SE, 66.0% at 1.25 times that SE and 50.7% at 1.5 times. The approximate SE ceiling for 80% power under these assumptions is 1.33 pp. External precision, target distribution and fitted performance may differ.

## E2. Rank agreement

Fixed score ranks come from the full 2011–2012 covariance fit. Realised ranks use survey-weighted held-out gain for all nine candidates in each block. Each of the 5,000 bootstrap draws resampled PSUs within strata, using the same replicate weights for every candidate and block. Each resampled PSU carries all its participant rows together. Percentile ranges are conditional on the predictions, score ranks and folds. They describe survey resampling of the held-out loss contrasts. With nine candidates per block and possible near ties, percentile ranges can be asymmetric or even exclude the point estimate. Their nominal coverage for population ranks has not been validated.

| Instrument | Pooled Spearman | Conditional 95% percentile range | Mean within-block Spearman | Conditional 95% percentile range |
|---|---:|---:|---:|---:|
| Frozen RBF | 0.851 | [0.756, 0.879] | 0.810 | [0.713, 0.837] |
| Matched affine | 0.940 | [0.662, 0.938] | 0.890 | [0.603, 0.910] |

The pooled affine slope was 1.058. The pooled RBF slope was 1.210. The equal-block mean of the within-block Spearman correlations is the proposed secondary summary. It compares ranks within the actual choice sets and avoids relying on how the five blocks are scaled against one another. Report its conditional survey-resampling range and the five block correlations. E2 has no separate pass threshold.

## E3. Swap contrasts

Keep this analysis descriptive. The saved aggregate output retains the signed loss contrasts and pointwise intervals. Against the development plug-in reference, the counts of half-widths smaller than the bar, intervals entirely above it, and intervals entirely within the symmetric bar are **0, 4, 0** of 40 adjacent contrasts and **0, 29, 0** of 40 winner-versus-alternative contrasts. The comparisons share the same participants and choices. The counts do not carry simultaneous coverage.

## Status and files

The frozen NOVA-S feasibility outcome and the later stability-screen outcome remain as recorded. Submission is paused. No 2013–2014 data file was opened, no study was frozen, and no file was pushed. The bundle retains both development diagnostics and includes the scripts and aggregate outputs under Development_only.
