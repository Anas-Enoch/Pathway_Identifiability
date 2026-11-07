# NOVA-D: exact population ordering of the frozen G1 rule — note v1.0 (7 Nov 2026)

**Status.** Read-only diagnostic artefact accompanying `NOVA_D_RESULT_v1.2.md` (owner's record, `60056edc…`),
whose §7.4 states that exact limiting rankings had not been checked. This note supplies that check. It is not a
new experiment, changes no frozen output, protocol or preregistration, and reruns no simulation. Prepared by
Claude (Fable 5.1); author of record Anas Enoch.

**Method.** `nova_d_exact_ordering.py` (sha256 `5143c8fc418c83ce657201f44e1074bbeb801840837c54627dc2bbe84563bea6`)
imports the frozen `nova_d_run.structure()` for each of the 72 conditions, builds the **exact** population
covariance of Z = [O, T, P] in closed form (Σ_OP = Σ_O aᵀ; Σ_TP = b Σ_O aᵀ + s²λ; Σ_TT = b Σ_O bᵀ + s²λλᵀ + ννᵀ +
σ²_ε I; Σ_PP = a Σ_O aᵀ + s² I; the Hermite term has E[u h(u/s)] = 0 and unit variance), and applies the frozen
`g1_pick` algebra to it — standardise, shrink S ← (1−α)S + α(tr S/d)I, condition on O **with the frozen 1e-6 ridge
on S_OO**, score_j = Σ_k C²_kj / v_j — at α = 0.05 (the rule's floor and therefore its large-n limit) and at α = 0.
No sampling is involved. The covariance builder is checked against the structure's equal-marginal-variance
construction: max |diag Σ_TT / c − 1| = 1.8 × 10⁻¹⁵ across all 72 conditions.
Output: `nova_d_out/exact_ordering.json`, computed on the owner's machine: sha256
`2b2fb3b14502dd2e623714d041d69c3aefb4d9d5f550a391a2a721d3c1c026a9`. A copy computed independently in the
drafting container (`3943ecb7…`) has identical integer, boolean and list fields in every row and differs only in
floating-point trailing digits (max |Δ| = 9.1 × 10⁻¹⁶, different BLAS); the owner-machine file is the record.

**Results (exact, all 72 conditions).**

| Quantity | Result |
|---|---|
| Frozen rule's population argmax ≠ Gaussian argmax | **15 of 72**: M3 — 11, 12, 13, 14; M2 — 21, 22, 24, 30, 32, 33, 35, 37, 39, 43, 44 |
| By Gaussian top-two gap | all 15 at Γ_rel = 0.02; **0 of 36** at Γ_rel = 0.20 (no M1 or M4 condition flips) |
| Where the eleven M2 flips land | on candidate 1 = j\*, the nonlinear optimum, in every case |
| Positive control (43) | **flips**: the frozen rule's limit pick is j\* |
| Primary (31) | does not flip; the frozen rule's top-two relative gap is 8.2 × 10⁻⁵ (the 2 % Gaussian gap compressed ~250-fold) |
| Same rule at α = 0 | population argmax equals the Gaussian argmax in **72 of 72** |
| Condition 3 (null, A 0.05), frozen top-4 | **[0, 3, 5, 1]** — the 200,000-row diagnostic in result v1.0/v1.1 reported [0, 3, 5, 10]; its fourth entry was finite-sample error. Candidate 10's 96 picks at n = 512 (α = 0.128 there) are finite-sample, not a population reordering |

**What this certifies.** Objective A only: with its α ≥ 0.05 floor, the frozen G1 rule is not a consistent
estimator of the population Gaussian acquisition ranking in this family — its limiting argmax differs from the
Gaussian argmax in 15 of 72 conditions, exclusively where the Gaussian top two are within 2 %, and the discrepancy
vanishes at α = 0 in every condition. In the positive control the limiting pick coincides with the nonlinear
optimum, which is why the finite-sample G1 reached j\* in 67 % of folds and why the realised opportunity for NOVA-S
against the frozen comparator was smaller than the designed ρ_pop = 0.40 (result §7.4).

**What this does not certify.** Anything about finite-sample decision loss (Objective B); any causal attribution
of the observed excess regret of the estimated G1 at small n to the shrinkage rule; anything about the archive's
71 % self-agreement. Those remain for the prospective G1 study, whose design decisions are recorded in result
v1.2 §8.3.
