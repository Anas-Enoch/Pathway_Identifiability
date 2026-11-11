# NHANES 2011–2012 value-transport feasibility

**Disposition.** Keep the reserve closed. The development diagnostic supports a pooled positive-gain question at the observed scale. Precision is insufficient for a tight value-transport criterion or a two-sided classification of swaps at the programme reference scale. This is an exploratory 2011–2012 diagnostic. It supplies neither a frozen GO/NO-GO decision nor an estimate of transport across cycles.

## Design

The script reads only the hash-pinned 2011–2012 FAS_G and DEMO_G data and the frozen block file, wrapper and G1 module, plus the saved aggregate development fit. It reproduces all 45 full-development scores and verifies that the 25 training-fold picks match the frozen G1 implementation and the full-development leaders. No participant-level files are written.

The five existing participant-disjoint folds were used to fit a matched affine reconstruction procedure. Within each training fold, centring, standardisation, covariance shrinkage and the observed-block ridge were learned from training rows. Base-panel predictions used the regularised conditional mean. Each augmentation added the corresponding candidate residual times its conditional regression coefficient. The held-out outcomes were squared target errors on the training-fold target scale. This predictor differs from the prior RBF procedure behind the published 15.8% development estimate.

The five blocks were pooled at the participant level before taking a ratio. Ratios divide improvement by base-panel loss and use fatty-acid subsample weights. Taylor linearisation uses the full survey design with 1,725 eligible participants, 14 strata, 31 PSUs and 17 degrees of freedom. The intervals are pointwise, conditional on the cross-fitted predictions. They omit training-set variability, development selection of the pairs, cross-cycle shift and any adjustment for multiple comparisons.

For a diagnostic comparison, the block reference was set to 0.10 A(1 − A), with A the block's estimated held-out gain for its full-development G1 leader. This is a development plug-in bar. Class counts use the same positive threshold for both directions.

## Results

- **E1, pooled value.** The gain for the five development choices was **11.46%**, with a pointwise 95% interval of **4.57% to 18.35%**. A positive-gain question is distinguishable in development. The interval half-width of 6.89 percentage points is too wide to treat this as evidence that a transported gain could be verified within a tolerance near the programme reference scale.
- **E2, calibration.** The pooled fold-matched score prediction was **12.06%**. Realised-minus-predicted gain was **−0.60 percentage points** with a pointwise 95% interval of **−7.56 to +6.36 percentage points**. Across the 45 block-candidates the mean absolute point calibration error was **1.26 percentage points**. This is descriptive and uses the matched covariance model and its own base-risk scale. Near-zero average error can mask opposite errors across candidates.
- **E3, consequences of a swap.** Of 40 adjacent contrasts, **3** had a 95% interval half-width smaller than their block's bar, **7** had an interval entirely above the bar in the predicted direction and **0** had an interval entirely inside the symmetric bar. For 40 winner-versus-alternative contrasts, these counts were **0**, **20** and **0**. No pair had an interval entirely below the negative bar. The paired estimates share participants, analytes and a development ranking. These counts carry no simultaneous or confirmatory interpretation.

| Block | Leader | Block reference (pp) | Leader minus runner-up (pp) | Pointwise 95% interval (pp) |
|---|---|---:|---:|---:|
| 1 | LBXTSA | 1.49 | 4.29 | [1.92, 6.66] |
| 2 | LBXEPA | 0.45 | −6.59 | [−20.07, 6.89] |
| 3 | LBXDTA | 1.06 | 1.29 | [−1.24, 3.83] |
| 4 | LBXEPA | 0.79 | −0.85 | [−18.91, 17.21] |
| 5 | LBXDTA | 0.98 | 2.76 | [0.87, 4.66] |

Only block 1's leader-to-runner-up interval lies wholly above its reference. Blocks 2 and 4 have negative point contrasts and intervals crossing zero. The covariance ordering alone does not establish how these pairs rank on held-out loss.

## Reproduction

The bundle's Development_only directory contains the hash-pinned script and aggregate output. Run the script with NumPy, pandas and SciPy and pass the project root as its only argument. The original 2011–2012 XPT files remain in the sibling nhanes_data directory. Two invocations produced byte-identical aggregate JSON. The prospective estimand, tolerance and decision rules remain to be specified before any authorisation to use 2013–2014 data files.
