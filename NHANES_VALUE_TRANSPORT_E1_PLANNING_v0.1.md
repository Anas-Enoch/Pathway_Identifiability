# NHANES fatty-acid measurement value in an adjacent survey cycle

Planning addendum v0.1 · 10 November 2026

This addendum refines the reporting rule in the 2011–12 feasibility record v0.2. It uses the recorded development estimate, its design standard error and public survey documentation. It is a draft planning record. No protocol has been frozen and no 2013–14 participant data file was opened.

## Primary outcome and fixed benchmark

Report the 2013–14 survey-weighted pooled reduction in squared reconstruction loss from adding the five choices fixed on 2011–12. The base-panel and augmented RBF models are refitted on the complements of participant-disjoint evaluation folds within 2013–14. The choices, 2011–12 target scales, predictor recipe and pooling rule are held fixed. The result addresses replication of measurement value in an adjacent cycle of the same survey, conditional on the reconstruction procedure as run.

Keep the benchmark at 75% of the frozen development point estimate of 0.15770794508230182. The full-precision benchmark is 0.11828095881172637, displayed as **11.83%**. The 75% multiplier came from a planning illustration. It has no independent assay-relevance justification. The benchmark remains a fixed number in the evaluation. An inference about retaining 75% of the true 2011–12 gain would require a two-cycle contrast with uncertainty from both estimates.

The primary report is the gain and its two-sided 95% design-based interval. The class is a predeclared reading of that interval. Compute the class before rounding. A lower endpoint strictly above the benchmark gives **exceeded**. An upper endpoint strictly below it gives **missed**. Otherwise, including endpoint equality, the result is **unresolved**. Data or design integrity failures would make the analysis unassessable and are separate from these three readings.

## Result sentences

- **Exceeded.** In NHANES 2013–14, the five fixed measurements reduced survey-weighted pooled squared reconstruction loss by [gain]% (95% design-based CI [lower] to [upper]%) under within-cycle refitting, exceeding the pre-declared benchmark of 11.83%.
- **Missed.** In NHANES 2013–14, the estimated reduction in survey-weighted pooled squared reconstruction loss was [gain]% (95% design-based CI [lower] to [upper]%) under within-cycle refitting, with the interval wholly below the pre-declared benchmark of 11.83%.
- **Unresolved.** In NHANES 2013–14, the estimated reduction in survey-weighted pooled squared reconstruction loss was [gain]% (95% design-based CI [lower] to [upper]%) under within-cycle refitting, and the interval included the pre-declared benchmark of 11.83%.

## Planning probabilities

The noncentral t calculation conditional on an evaluation gain equal to the development estimate gives 84.2% power with 17 degrees of freedom and 83.6% with 15, at the development standard error of 1.2535 percentage points. At 15 degrees of freedom and 1.25 times the development standard error, conditional power is 65.3%. The conditional standard error ceiling for 80% power is 1.315 percentage points at 15 degrees of freedom (1.326 at 17).

To examine uncertainty about the evaluation gain, take it to vary normally around the recorded development point estimate with an assumed standard deviation τ. Let s be the assumed evaluation standard error, independent of that gain variation, with the usual chi-square model for its estimator. Under this scenario model the evaluation estimate has marginal standard deviation √(s² + τ²). The exact outcome probabilities under these assumptions follow from a noncentral t with noncentrality equal to the mean-minus-benchmark difference divided by √(s² + τ²), and signed cutoffs ±t₀.₉₇₅,ν s/√(s² + τ²). This checks the accompanying simulation without Monte Carlo error. The assumed τ is a sensitivity parameter. The development standard error does not estimate cross-cycle variation.

| Degrees of freedom | Assumed τ / development SE | Evaluation SE / development SE | Exceeded | Unresolved | Missed |
|---:|---:|---:|---:|---:|---:|
| 17 | 1.0 | 1.00 | 76.76% | 23.22% | 0.017% |
| 17 | 1.0 | 1.25 | 62.90% | 37.07% | 0.027% |
| 15 | 1.0 | 1.00 | 76.30% | 23.69% | 0.017% |
| 15 | 1.0 | 1.25 | 62.35% | 37.62% | 0.027% |
| 15 | 1.5 | 1.25 | 60.31% | 39.48% | 0.205% |

The scenario probabilities are conditional on the stated assumptions. They are not estimated chances of success in 2013–14. With 15 degrees of freedom and an observed standard error equal to the development value, an observed gain would need to exceed about 14.50% for an exceeded reading or fall below about 9.16% for a missed reading. Unresolved is therefore an ordinary result under these planning scenarios and will be reported as such.

## Survey design documentation

Use **15 nominal degrees of freedom for planning**. The public DEMO_H codebook lists masked stratum codes 104–118 and PSU codes 1–2. NCHS documents 30 sampled PSUs for the cycle. The FAS_H documentation specifies the fatty-acid subsample weight WTFAS2YR and links to the masked design variables in DEMO_H. [DEMO_H codebook](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2013/DataFiles/DEMO_H.htm#SDMVSTRA) · [NCHS sample design report](https://www.cdc.gov/nchs/data/series/sr_02/sr02-190.pdf) · [FAS_H analytic notes](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2013/DataFiles/FAS_H.htm)

The codebooks do not give the joint counts for the eligible fatty-acid analysis domain. In the eventual authorised run, retain the full relevant subsample design for Taylor domain variance and tabulate the strata and nested PSUs represented in the eligible domain. Use the domain PSU count minus the domain stratum count for the primary t reference and report those counts and the realised degrees of freedom. Finalise integrity and unsupported-design handling before freezing. [NCHS analytic guidelines, §§3.2.3.1–3.2.3.2](https://wwwn.cdc.gov/Nchs/data/nhanes/analyticguidelines/11-16-analytic-guidelines.pdf)

## Secondary and descriptive outcomes

- **E2.** Report the equal-block mean of the five within-block Spearman correlations, the individual correlations, and the conditional survey-bootstrap percentile range. Use common replicate weights across all participants, candidates and blocks and resample PSUs within strata. Call the percentile range a resampling summary. It has no validated confidence-coverage guarantee for population ranks. There is no pass threshold. In 2011–12, block 2 had the lowest correlation, 0.58.
- **E3.** Report signed paired swap contrasts descriptively. The comparisons share participants and candidates. They do not determine the E1 class.

This addendum leaves the earlier manuscript, the frozen studies and the submission pause in place. It supplies reporting language and planning assumptions only.
