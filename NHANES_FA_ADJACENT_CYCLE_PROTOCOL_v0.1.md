# Replication of fixed fatty-acid measurement value in an adjacent NHANES cycle

Protocol v0.1

## Scope and fixed choices

This protocol consolidates the 2011–12 feasibility record v0.2 and E1 planning addendum v0.1. The independent evaluation is in NHANES 2013–14, an adjacent cycle of the same survey. It asks how much held-out reconstruction loss the five fixed block choices reduce when the predictor recipe is refitted within that cycle. It does not test selector superiority or the transfer of fitted coefficients. The earlier NOVA-S NO-GO and stability sharpness closeout remain separate.

| Block | Fixed choice | Analyte |
|---:|---|---|
| 1 | LBXTSA | Tricosanoic acid |
| 2 | LBXEPA | Eicosapentaenoic acid (EPA) |
| 3 | LBXDTA | Docosatetraenoic acid (DTA) |
| 4 | LBXEPA | Eicosapentaenoic acid (EPA) |
| 5 | LBXDTA | Docosatetraenoic acid (DTA) |

The blocks, their roles and candidate order are the frozen v0.4 allocation. Use a common complete-case domain across their 29 analytes, age at least 12 and positive `WTFAS2YR`. Join fatty-acid and demographic records one to one by SEQN. The order of the fatty-acid file is preserved through the join and evaluation. Use the development target divisor per block and target, namely the maximum of the saved target standard deviation and saved floor. Candidate score ranks are the full 2011–12 G1 ranks from the aggregate fit. Rank ties use average ranks for Spearman and frozen candidate order for ordered comparisons.

Fit base-panel and base-plus-one-candidate RBF models under the frozen N2 recipe. Five outer folds use `folds_for` with cycle seeds 20261103 for the development dry run and 20261104 for evaluation. All blocks share the participant folds. The frozen inner tuning runs on each outer training set, with `20261019 + k` for zero-based outer fold `k`, and uses squared loss. No evaluation fold informs its fitted predictor or tuning. Each eligible participant contributes one held-out loss per block and arm. G1 refit picks are computed for the descriptive disagreement count. Only the fixed choices enter E1.

## E1 and the survey variance rule

Let `b_i` and `a_i` be participant `i`'s sums over the five blocks of squared target errors for base and fixed-choice reconstruction, respectively. E1 is `sum(w_i (b_i - a_i)) / sum(w_i b_i)` with fatty-acid subsample weights. Report the estimate and two-sided 95% Taylor design interval first. The variance conditions on the out-of-fold predictions and measures sampling variability of their loss contrast. It does not include additional refit variability or uncertainty in the development choices and ranks.

Use the full positive-weight fatty-acid file design, retaining children outside the age domain as zero-contribution rows. The analysis describes the complete-case domain and does not adjust for analyte item nonresponse. The Taylor influence of a row outside the eligible domain is zero. Sum influences within nested masked PSUs. Within each stratum use the usual with-replacement factor `m_h/(m_h−1)` on the sum of squared deviations of the full-design PSU totals from their mean. A stratum represented by just one PSU **in the eligible domain** still uses the full design and its zero-contribution PSU totals. This stratum adds zero to the degrees of freedom. A stratum with eligible participants but only one PSU **in the full positive-weight design** makes the run unassessable. A stratum with no eligible participants contributes zero and does not enter the degrees of freedom. There is no stratum merging or variance imputation. The t reference uses the number of nested PSUs represented in the domain minus the number of strata represented there. Report both counts and full and domain PSU support by stratum. The nominal 2013–14 planning value is 15. The actual domain counts determine the realised value. [NCHS analytic guidelines, §§3.2.3.1–3.2.3.2](https://wwwn.cdc.gov/Nchs/data/nhanes/analyticguidelines/11-16-analytic-guidelines.pdf) · [FAS_H analytic notes](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2013/DataFiles/FAS_H.htm) · [DEMO_H codebook](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2013/DataFiles/DEMO_H.htm)

Keep the full-precision benchmark `0.11828095881172637`, displayed as **11.83%**. It is 75% of the recorded development point estimate, `0.15770794508230182`. The 75% fraction came from a power illustration and has no independent assay-relevance justification. The fixed-benchmark comparison does not propagate uncertainty from the earlier cycle. Classify on the unrounded interval. A lower endpoint strictly above the benchmark gives **exceeded**. An upper endpoint strictly below it gives **missed**. Otherwise, including equality, the reading is **unresolved**. An integrity or design failure gives **unassessable** without an E1 class.

- **Exceeded.** In NHANES 2013–14, the five fixed block choices (three distinct analytes) reduced survey-weighted pooled squared reconstruction loss by [gain]% (95% design-based CI [lower] to [upper]%) under within-cycle refitting, exceeding the pre-declared benchmark of 11.83%.
- **Missed.** In NHANES 2013–14, the estimated reduction in survey-weighted pooled squared reconstruction loss for the five fixed block choices (three distinct analytes) was [gain]% (95% design-based CI [lower] to [upper]%) under within-cycle refitting, with the interval wholly below the pre-declared benchmark of 11.83%.
- **Unresolved.** In NHANES 2013–14, the estimated reduction in survey-weighted pooled squared reconstruction loss for the five fixed block choices (three distinct analytes) was [gain]% (95% design-based CI [lower] to [upper]%) under within-cycle refitting, and the interval included the pre-declared benchmark of 11.83%.

## Closed list of outputs

The evaluation writes one aggregate JSON. Its provenance fields contain the status, cycle, file hashes, runner and frozen-input hashes, software versions and aggregate design counts. Its scientific fields are limited to the following.

| Output | Contents |
|---|---|
| E1 | Pooled squared-loss gain, design SE, 95% interval, benchmark and class |
| Block gains | Five fixed-choice squared-loss gains, design SEs and pointwise 95% intervals |
| Bounded-loss sensitivity | Pooled gain, design SE and 95% interval after capping each squared target error at 25 before summing. No class |
| E2 | All 45 candidate gains, five within-block Spearman correlations and their equal-block mean. Report its conditional 2.5th–97.5th percentile range and the number of defined replicates. No threshold |
| E3 | Signed held-out squared-loss contrasts for the 40 adjacent development ranks and 40 leader-versus-alternative pairs, each divided by its block's weighted base loss. Positive values favour the higher-ranked candidate. Report design SEs and pointwise 95% intervals descriptively. No threshold or simultaneous coverage |
| Observed refit disagreement | Number of unequal G1 picks among the ten within-block pairs of outer-training refits, by block and total out of 50. Descriptive only |

For E2, use 5,000 rescaled PSU bootstrap replicates with NumPy seed 2026111002, sorted masked strata and PSUs and common replicate weights for every block and candidate. In a stratum with eligible rows, draw `m_h−1` full-design PSUs with replacement and rescale multiplicities by `m_h/(m_h−1)`. A sampled PSU carries all of its eligible participants. Retain predictions, folds and development ranks as fixed. Use linear interpolation for the percentile range. A block correlation is undefined if either rank vector is constant. A replicate is undefined if any block correlation is undefined or any bootstrap base denominator is zero. Report the point mean only if all five block correlations are defined and the range only if the point mean and all replicates are defined. The range is a resampling summary and its population-rank coverage has not been validated. E2 and E3 share participants and candidates. All other analyses after opening are exploratory.

## Integrity, dry runs and access sequence

A mismatch in a pinned input, missing file or required field, wrong release cycle, duplicate or malformed SEQN, missing demographic match or design variables, overlap with an eligible development SEQN, invalid age, weight or design code, fewer than ten eligible rows, invalid fold partition, unusable full-design stratum, nonpositive domain degrees of freedom, missing or nonfinite predictions or contrasts, a nonpositive weighted base denominator, or a zero primary design SE makes the run unassessable. Missing analyte values exclude rows under the common domain rule and do not by themselves make the run unassessable. An unassessable run records only the aggregate failure status and reason. No participant rows, predictions, SEQNs or fold allocations are written. An interrupted evaluation leaves the opening attempt in place and does not authorise a second run.

The runner is `nhanes_fa_adjacent_cycle.py`. Dry-run `--cycle G` against the pinned 2011–12 files and `--cycle synthetic` against independent XPT fixtures using the 2013–14 layout and fold seed. The only scientific input that changes is the cycle. The development path must reproduce all 25 frozen G1 refit picks, the recorded squared and bounded gain and SE at relative tolerance `10⁻⁷` and absolute tolerance `10⁻¹⁰`, and 1,725 eligible participants with 17 design degrees of freedom. The synthetic run exercises the same fits, pooling, variance and closed output list. Only a separately authorised real evaluation receives an E1 class.

Before downloading either evaluation file, publicly pin the protocol, runner and dry-run record at a Git commit. Obtain and verify an independent timestamp over a manifest binding that commit and its file digests. Keep the 2013–14 files sealed until a separate authorisation. Record source URLs, UTC download times and SHA-256 checksums upon download, before parsing. The opening authorisation JSON supplies `cycle` (`2013–14`), `protocol_sha256`, `runner_sha256`, `public_commit`, `timestamp_receipt_sha256`, `timestamp_utc`, `authorised_utc` and a `downloads` object keyed by `FAS_H.xpt` and `DEMO_H.xpt`, each with `download_started_utc` and `sha256`. The operator verifies the independent receipt against the manifest and public commit before preparing this record. The runner checks the recorded chronology and hashes and creates a local exclusive attempt record before touching either file. Run only once in the designated checkout, with `--cycle H` and the authorisation. Retain the attempt and result or failure together. No selector development or alternate rule follows an unresolved or unassessable result.

## Planning interpretation

At the development SE of 1.2535 percentage points and nominal 15 degrees of freedom, the conditional power is 83.6% if the true evaluation gain equals the recorded development estimate. It is 65.3% at 1.25 times that SE. The corresponding 80% SE ceiling is 1.315 percentage points (1.326 at 17 degrees of freedom).

In the planning sensitivity model, the true evaluation gain is normal around the development point estimate with assumed SD `τ`. With assumed evaluation SE `s`, an independent chi-square model for its estimator and `σ = √(s² + τ²)`, the outcome probabilities use a noncentral t with noncentrality `(development estimate − benchmark)/σ` and cutoffs `±t₀.₉₇₅,ν s/σ`.

| df | τ / development SE | Evaluation SE / development SE | Exceeded | Unresolved | Missed |
|---:|---:|---:|---:|---:|---:|
| 17 | 1.0 | 1.00 | 76.76% | 23.22% | 0.017% |
| 17 | 1.0 | 1.25 | 62.90% | 37.07% | 0.027% |
| 15 | 1.0 | 1.00 | 76.30% | 23.69% | 0.017% |
| 15 | 1.0 | 1.25 | 62.35% | 37.62% | 0.027% |
| 15 | 1.5 | 1.25 | 60.31% | 39.48% | 0.205% |

These are scenario probabilities. The development SE does not estimate cross-cycle shift. If the observed evaluation SE equals the development SE at 15 degrees of freedom, an observed gain must exceed about 14.50% to exceed the benchmark or fall below about 9.16% to miss it. Unresolved is an ordinary reading under these assumptions.
