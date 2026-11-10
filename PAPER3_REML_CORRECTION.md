# Paper III: correcting the numerical failure in the frozen REML

Date: 2026-10-16. This is a numerical correction of the frozen estimator's implementation. It is
not a new analysis: the estimand, data, model and decision grid are unchanged.

**Inputs.**
- `paper3_inference/p3k_study_bootstrap.csv`, sha256 `c2e5562f…f5ea5`. This matches the hash in
  `p3k_inference_manifest.json`.
- The frozen script `p3k_bootstrap_reml_primary.py` (`e0c6568e…9576`); its `reml_fit` is identical
  to v1's.
- `p3_blocks_per_record.csv`, from `g0_blocks.jsonl`.

**Check script.** `p3k_reml_converged_check.py` v1.0.0, sha256
`7001dea669059842a75a61062d89c06efebce1aaca26b7b64d5ea2964df9918b`. Run here, it writes
`p3k_reml_converged_check.json` (`b277abcd…d711d`).

**Reproduced on the owner's Mac on 2026-10-16 at 20:23.**
- The script hash matched (`7001dea6…`).
- All seven rows were identical.
- The output JSON sha256 is `b277abcd9c29205af01bf166856b23ff3a816ac2de899d441d5c6ce8228d711d`, byte-identical
  to the run here.

**Owner's decision (2026-10-16).** The converged REML is reported as Paper III's primary result,
because it is the preregistered estimator computed to convergence. The output as executed and this
defect are disclosed as an implementation correction.

## Defect

The frozen `reml_fit` searches τ² on a linear scale:
- the bracket is [0, max(scale·1e6, 1)];
- the tolerance is xatol = scale·1e-12, where scale = max(var(D), median v).

ST005094 has D = 337,247. That gives var(D) ≈ 6.1e8, so xatol ≈ 6.1e-4. The tolerance is wider
than the true τ².

The exact replica of the frozen search reproduces the reported τ² bit for bit: 0.00040221477618907316.
On synthetic data with the same outlier, it returns this same value whatever the true τ² is.

## Restricted likelihood

The likelihood has an interior mode at τ² = 1.030e-4 (nll −537.4).

Four studies have v = 0 and identical D = 0 (ST001906, ST001908, ST002211, ST003031). Because of
them, the likelihood is formally unbounded as τ² → 0+. It drops below the interior mode only for
τ² < ~1e-150, which has no meaning in practice. The interior mode is reported. The fit that
excludes those four studies has no degeneracy and gives the same classification.

## Results (frozen grid; the CI is 95% Wald, as frozen)

| Analysis | k | τ̂ | μ̂ | 95% CI | 95% prediction interval | Grid |
|---|---|---|---|---|---|---|
| Frozen search, as run | 186 | 0.0201 | 0.01019 | [0.00660, 0.01377] | [−0.0293, 0.0497] | statistically positive, practically subthreshold |
| **Converged REML (interior mode)** | 186 | **0.0102** | **0.00867** | **[0.00659, 0.01075]** | **[−0.0113, 0.0287]** | statistically positive, practically subthreshold |
| Zero-variance studies excluded | 182 | 0.0103 | 0.00910 | [0.00695, 0.01125] | [−0.0112, 0.0294] | same |
| Records with RMSE > 10 excluded (n = 20) | 166 | 0.0105 | 0.00917 | [0.00698, 0.01136] | [−0.0115, 0.0298] | same |
| Both exclusions | 162 | 0.0106 | 0.00963 | [0.00737, 0.01190] | [−0.0113, 0.0306] | same |
| Cohort-clustered (13 runs collapsed) | 162 | 0.0110 | 0.00898 | [0.00655, 0.01141] | [−0.0127, 0.0307] | same |
| Cohort-clustered, RMSE > 10 excluded | 146 | 0.0113 | 0.00950 | [0.00696, 0.01204] | [−0.0128, 0.0318] | same |

The Knapp–Hartung intervals were computed for every row and do not change any classification.

## What changes in the reporting

- **Unchanged:** the classification, in every analysis.
- **The primary estimate is smaller than reported.** μ̂ = 0.0087, which is 62% of δ_min. The run
  as executed reported 0.0102 (73%).
- **The upper CI bound now sits clearly below δ_min:** 0.0108 against 0.014, where the run as
  executed reported 0.0138.
- **Between-study heterogeneity is half the reported value:** τ̂ = 0.010 instead of 0.020. The
  prediction interval narrows to [−0.011, 0.029].
  - The share of studies whose true effect is negative falls from about 31% to about 20% under the
    normal model.
  - That figure is still an upper bound, because the evaluation-only bootstrap leaves design-split
    noise inside τ².
- **The paper should report both:** the frozen output as executed, and the converged estimate,
  with this defect explained.
- **The prespecified equal-study mean (§13), 1813.3, is uninformative.** It is driven by
  ST005094's scale pathology and should be reported as such.

## Paper IV

`p4f_reml_primary.json` gives τ² = 1.35e-13 and μ = 4.85e-11, with 93 zero-variance studies.
With 93 studies at v = 0 and identical D = 0, the restricted likelihood falls by about 92 log units
per decade as τ² → 0. The degeneracy is therefore probably real, and Paper IV's switch to the
equal-study bootstrap was justified. This still needs confirming against `p4f_study_bootstrap.csv`.

## Within-study variance in the B0 setting (descriptive; 162 studies with v > 0 and RMSE ≤ 10)

The fit is: log v = −6.38 − 0.61 log n_eval − 0.13 log blocks.
- Standard errors: 0.25 for log n_eval and 0.22 for log blocks.
- R² = 0.04. Residual SD = 2.5 on the log scale, a factor of about 12 on v.
- Design size explains almost none of v; it is driven by the study itself.
- The block count barely matters. The 1/b scaling assumed in B0 v0.2 is not supported.

The B0 gate is therefore recomputed from the empirical distribution of v and the converged τ²
(B0 v0.3).
