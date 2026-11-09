# Acquisition external-target charter v0.1 — translation audit only

9 November 2026. Retrospective decision support from versioned summaries. This
draft does not declare a target cohort, Π, Θ_cap, ε or τ, set eligibility rules,
modify the external stub, or authorize a study. The archived grid and raw runs
were not consulted; the G1-R numbers below are from its published-style result
record's summary table. No new data were accessed.

The [scope addendum](ACQUISITION_SUCCESSOR_SCOPE_v1.1.md) records the role of
this audit in the programme closeout and links to the programme synthesis.

## 1. Put the relevance bar and the headroom in the same units

Let L₀ be the no-additional-measurement held-out loss, L_G1 the held-out loss
after G1's acquisition, and L_h the alternative's loss, pooled on the same
evaluation units with the same weights and loss. The existing external stub is
for NOVA-S versus G1. Its relevance factor is a useful reference scale here,
not a successor threshold adopted by this audit.

    A_G1 = (L₀ − L_G1) / L₀
    Δ_ext = (L_G1 − L_h) / L_G1
    I₀ = (L_G1 − L_h) / L₀ = (1 − A_G1) Δ_ext

Thus Δ_ext ≥ 0.10 A_G1 is equivalent to I₀ ≥ 0.10 A_G1(1 − A_G1), assuming
positive denominators. This is a change of units only. The original lower-bound
requirements on baseline gain and superiority remain in the original criterion.

| Actual G1 gain A_G1 | Relevance bar in G1-loss units | Same bar in no-acquisition-loss units |
|---:|---:|---:|
| .05 | .005 | .00475 |
| .15 | .015 | .01275 |
| .20 | .020 | .01600 |
| .158, recorded NHANES development example | .0158 | .0133036 |

In particular, 1.58% of G1's remaining loss and 1.33036% of the
no-acquisition loss are the **same bar** at a baseline gain of 15.8%.

In G1-R the calibration axis A is **oracle** gain, A★ = G★/R₀. Under the
study's common Bayes loss, expected G1 gain is A_G1^Bayes = A★ − r₀, where r₀
is A0's expected normalized selection regret. The expected selection gain of a
replacement h is r₀ − r_h. Therefore, for a matched synthetic population and
this loss, a necessary oracle headroom condition for meeting the analogous bar is

    r₀ ≥ 0.10 (A★ − r₀) [1 − (A★ − r₀)].

This bound assumes positive baseline gain and the same target, candidate set,
budget and Bayes loss. It says nothing about differences in fitted-predictor
error, target scaling, sampling weights, assay noise or transport. R₀ is the
remaining risk with the observed panel; a fraction of R₀ is not a fraction of
explained or predictable variance. An external bridge must specify when the
Bayes selection contrast approximates the deployed predictive contrast.

### What the reported synthetic averages show

These are differences of the displayed **all-condition** means at each training
size, using the same 72 conditions for both arms. They are not maxima, estimates
stratified by A★, or guarantees at a given training size.

| Training n | Mean A0 regret | Mean A2 regret | Mean recovery, A0 − A2 |
|---:|---:|---:|---:|
| 102 | .00916 | .00589 | .00327 |
| 410 | .00376 | .00223 | .00153 |
| 820 | .00244 | .00139 | .00105 |
| 1640 | .00146 | .00090 | .00056 |
| 3280 | .00145 | .00059 | .00086 |

As a scale illustration, .00327 is about 69% of .00475, the bar if **actual
G1 gain** is .05. This juxtaposes a mixed-grid mean with a hypothetical bar;
it does not estimate recovery in the A★=.05 stratum. Dividing the all-condition
mean by the mean oracle amplitude instead gives .0245 at n=102 and .0115 at
n=410 in recovery per unit oracle amplitude. These are ratios of averages.
Amplitude-specific conclusions require matched estimates of recovery and the
G1 gain in each stratum.

The reported maximum of .028376 at n=410 and A★=.15 **exceeds the .01275
reference bar by a factor of 2.23**. This bar uses .15 as actual G1 gain for
the comparison. In the matched Bayes setting with positive G1 gain, .01275 is
an upper bound on the bar because A_G1^Bayes ≤ A★=.15. The reported maximum
was selected for being largest after the grid outcomes were seen. Its crossing
of the bar is precisely the outcome that needs independent confirmation before
using that condition to support a target or a candidate. Neither the means nor
the selected maximum define a target distribution or the benefit and harm budgets.

### External evidence and the target decision

| Proposed setting | What the versioned records establish | What remains open |
|---|---|---|
| NHANES serum fatty acids | The **2011–12 development cycle** had 1,725 eligible participants, 0/25 block-fold departures and 15.8% baseline gain. The reported feasibility decision was NO-GO for the NOVA-S comparison. The 2013–14 evaluation cycle was not downloaded in that study. | Neither the gain nor agreement is measured for 2013–14 in this record. There is no external headroom estimate or justified harm-only assignment to Θ_cap. No new data access is proposed. |
| Small targeted clinical panels | A development stability check recorded changes in G1 choices in 29% of block-fold comparisons. The documentary audits found no eligible route under their stated strict rules in the screened samples. | Instability is not improvement opportunity, and the record does not establish that n≈50–200 or low gain identifies a useful regime. A study with documented sampling and a concrete panel decision could be informative; unknown participant linkage does not warrant treating rows as independent even for supportive inference. |
| LIPID / COMETS | Prior project records considered these as distinct access routes and weighed governance, panel compatibility and independent study supply. | The historical route assessment is not a current eligibility or access check. Cohort size is not necessarily training n. A synthetic mean cannot rule out meaningful gains for these panels or establish that their use is warranted. No contact or allocation is proposed. |

The NOVA-S development gain was +0.112% in its own equal-record relative-loss
endpoint, below the 0.29% development detectability/qualification scale. That
scale is not a clinical usefulness threshold. Neither those two numbers nor
the NHANES agreement can be transported to the synthetic comparison by pooling
them or by assuming a common denominator across studies.

**Operational gate.** An external driver must name the measurement action,
the loss used to judge it, and the minimum change in that loss that would change
the action. Its owner must identify the population and panel, acceptable harm,
and an independent evaluation path with verifiable sampling units. Because the
bar is relative, a small measured G1 gain lowers its numerical value but does
not qualify a population for the gate. Eligibility comes from the external use
case. The 5% interpretation floor is **not an empirically validated biological
threshold**.

No such driver is established in the records used here. That supports pausing
successor design and preparing the programme synthesis. It does not establish
that no relevant use case exists, that only small cohorts are eligible, or that
regularization cannot help any external population. If an independent driver is
supplied, the prospective declarations must be made separately, without using
this retrospective screen to set Π, Θ_cap, ε or τ.

### Source records

- [External stub §§2–4](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/EXTERNAL_VALIDATION_PREREG_STUB_v0.1.md), [relevance confirmation](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/EXTERNAL_RELEVANCE_FACTOR_CONFIRMATION_v1.0.md) and [NHANES protocol, baseline denominator](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/NHANES_FA_PROTOCOL_v0.4.md:51).
- [G1-R summary and reported exploratory maximum](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/G1R_RESULT_v1.1.md:95); [oracle calibration in the implementation](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/g1r_run.py:143).
- [NHANES development record](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/NHANES_FA_RESULT_v1.0.md), [selector evidence matrix](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/SELECTOR_CLAIM_EVIDENCE_MATRIX_v0.4.md) and [historical external-route memo](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/EXTERNAL_VALIDATION_DECISION_MEMO_v1.0.md).
