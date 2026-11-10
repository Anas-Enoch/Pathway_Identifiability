# NHANES fatty-acid G1 stability prediction, v0.4 closeout

10 November 2026. Local record for the development-only stability prediction
under protocol v0.4. The screen completed at 2026-11-10T04:50:59Z.

## Decision

**Not sharp enough to run.** The Gaussian upper prediction limit was U = 6 at
every integer sample size from 1,200 through 2,400. The prespecified cap stays
at U ≤ 5. Zero of 1,201 sizes are assessable, so the screen did not permit a
freeze or an external run. This is a development sharpness decision. No external
calibration verdict was issued.

The six-worker run completed all 2,402 model-size parts, with 2,000 simulated
cycles per part (4,804,000 total). At the development size of 1,725, the observed
disagreement count C = 0 was inside the Gaussian prediction interval [0, 6].
The predicted modal picks matched all five development picks, and all five
full-fit signed gaps were positive. The empirical resampling interval at that
size was [0, 8]. The Monte Carlo proportions for C = 0 were 0.864 (Gaussian)
and 0.627 (empirical).

## Empirical upper-limit corrections

The table below is an erratum to the saved empirical-model output. At B = 2,000
the declared 97.5th-percentile rule selects the smallest count C whose
cumulative histogram count reaches 1,950. For the five entries below, the
cumulative count at the correct limit was exactly 1,950. The floating-point
CDF fell just below 0.975 at that count and the implementation selected a
higher limit. The lower limits are unchanged.

| Size n | Saved U | Correct U | Cumulative count below correct U | Cumulative count at correct U |
|---:|---:|---:|---:|---:|
| 1,437 | 10 | 8 | 1,919 | 1,950 |
| 1,476 | 10 | 8 | 1,919 | 1,950 |
| 2,005 | 8 | 6 | 1,643 | 1,950 |
| 2,013 | 8 | 6 | 1,635 | 1,950 |
| 2,169 | 8 | 6 | 1,631 | 1,950 |

The integer-count check found no Gaussian quantile discrepancies. These five
empirical corrections leave the Gaussian sharpness decision unchanged. They
were computed from saved integer histograms without rerunning the simulations.

## Saved outputs

- [Prediction table](nhanes_fa_prediction_v0.4.json), SHA-256 `dc8f83c845fb1b8be4454f7c83e0723220d0d9bb404d3750bd378e5bd14501fa`
- [Screen report](nhanes_fa_stab_out/screen_v0.4.json), SHA-256 `cbf4a14d29f179aec03a301d5f2521eaa80054c411389b28d0183253dc8bcf39`

The saved table, screen, and simulation parts retain their original values.
Participant-level outputs remain local. NHANES 2013–14 remains unopened. This
record creates no freeze, authorisation, or commit.
