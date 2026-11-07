# G1-R — Covariance regularisation and finite-sample acquisition decisions: protocol **draft** v0.1 (7 Nov 2026)

**Status: DRAFT for owner review. Not frozen, not hashed; no run may occur under a draft.**
Author of record: Anas Enoch. Drafted by Claude (Fable 5.1) under the owner's decisions of 7 Nov 2026
(`NOVA_D_RESULT_v1.2.md` §8.3; review of `NOVA_D_RESULT_v1.1`; review of `g1_shrinkage_equivalence.py`).
Open items for owner decision are in §13.

## 1. Question and scope

> When does covariance shrinkage improve finite-sample acquisition decisions, and when does its distortion of the
> conditional covariance induce persistent ranking errors?

Three evidential levels, kept apart throughout:
- **Objective A — population ranking consistency.** Does the rule's population argmax equal the Gaussian argmax?
  Decided by exact calculation (§8), not simulation.
- **Objective B — finite-sample decision regret.** Expected regret of the selected candidate at finite n_train.
  The simulation's only question.
- **Robustness to ill-conditioning.** Regret as a function of the observed-block condition number.

**Scope: Gaussian only.** Z = [O, T, P] jointly Gaussian; V_j = G_j; the oracle optimum is j\* = argmax_j G_j and
is the same for every arm. Not a test of nonlinear acquisition, not a NOVA-S rescue, and nothing about biological
data. The 186-record archive is not used; no NOVA-D output is rerun or reused except the frozen `structure`
machinery of `nova_d_calib.py` for calibration.

**Established before this draft (not re-tested here):** the frozen rule's population argmax differs from the
Gaussian argmax in 27 of the 72 NOVA-D conditions at their operative α (`NOVA_D_EXACT_ORDERING_NOTE_v1.2.md`),
and the frozen joint shrinkage is exactly a ridge-conditioned estimator (§4, A1) — `g1_shrinkage_equivalence.py`,
committed `05ec454`.

## 2. Hypotheses

**H1 (primary, Objective B where A0 is inconsistent).** At the designated condition D1, where the frozen rule's
population argmax at its operative α is not the Gaussian argmax, the decoupled arm A2 has lower expected regret
than the frozen rule A0.

**H2 (no-harm, Objective B where A0 is consistent).** At the designated condition D2, where the frozen rule's
population argmax equals the Gaussian argmax and the Gaussian gap is wide, A2's expected regret is not greater
than A0's by more than the prespecified margin (§10). The frozen rule's heavier shrinkage may *win* here; that is
the bias–variance trade-off the study must be able to report, not a failure of the study.

Superiority of decoupling is a hypothesis to test. No claim that conditioning first dominates at every n.

## 3. Primary outcome

For a selector that picks ĵ from a training sample, the **acquisition decision loss** is the exact population
regret of that pick,

    regret(ĵ) = (G_{j*} − G_ĵ) / R_∅ ,

computed from the condition's population covariance (Pythagoras: R_S = R_∅ − V_S, and V_j = G_j in this family).
No evaluation sample, no fitted predictor, no cross-fitting. Per condition and arm: mean regret over replicates;
paired differences A_k − A0 with a bootstrap over replicates; P(ĵ = j\*); the distribution of picks.
Covariance-estimation error is reported only as a secondary (Frobenius distance of each arm's effective matrix to
Σ), to exhibit — not assume — that lower estimation error need not mean lower regret.

## 4. Arms

All arms receive the same standardised training sample Z (columns scaled to unit sample variance, as the frozen
rule does) and the frozen scoring form score_j = Σ_k C²_kj / max(v_j, 10⁻¹²). ε = 10⁻⁶ is the frozen ridge on
the observed block. d = 60 is the dimension of the standardised matrix.

| Arm | Definition | Role |
|---|---|---|
| **A0 — frozen G1** | `p5_gate3a.g1_pick` (sha256 `4600e2f3…`), hash-verified, called unchanged. α_n = min(1, max(0.05, d/(d + n_train))). | comparator; retained on the record |
| **A1 — equivalence control** | ridge-conditioned: λ = (α_n + ε)/(1 − α_n) on the O-regression; floor α_n/(1 − α_n) on each residual variance; scores × (1 − α_n). Derived identity; **gate**: must equal A0's scores to 10⁻¹² on every fold-specific α, every condition's exact Σ, and every sample S (§11). | certifies the decomposition; produces no separate outcome |
| **A2 — decoupled conditional regularisation (primary)** | two independently regularised components: (i) conditioning ridge λ_c on S_OO chosen by **GCV** for the multivariate ridge regression of [T, P] on O (closed form); (ii) residual-variance stabilisation v_β,j = (1 − β)v_j + β·mean_j v_j with **β = 0 designated**; cross-covariances C are **not** shrunk. | primary |
| **A2′ — decoupled, λ_c = 0** | OLS conditioning (n_train > d_O = 20 always), β = 0 | secondary: the consistent, unregularised estimator |
| **A3 — vanishing floor** | frozen rule with α_n = min(1, d/(d + n_train)), no 0.05 floor | control; differs from A0 only when n_train ≥ 1140 |
| **A4 — Ledoit–Wolf** | joint shrinkage toward (tr S/d)·I at the Ledoit–Wolf (2004) analytical intensity, then the frozen conditioning | control: Frobenius-optimal shrinkage versus decision loss |

**β, specified before any invariance is claimed.** β acts on the residual variances only; its target is their
mean; admissible range [0, 1]. At β = 1 every denominator equals mean_j v_j and the score reduces to Σ_k C²_kj,
which does not tie generically. No ranking invariance in β is claimed; β ∈ {0, 0.5, 1} is a tested secondary
knob, with β = 0 the designated primary. Shrinking the full residual covariance toward identity is **not** done,
because at complete shrinkage its cross-covariances vanish and all scores tie.

**α = 1.** The frozen rule's min(1, ·) permits α = 1 only when d/(d + n_train) ≥ 1, impossible for n_train ≥ 1;
the implementation asserts α_n < 1 in every condition, and A1's λ is therefore always finite.

**Denominator clipping.** The frozen max(v_j, 10⁻¹²) is retained exactly in every arm. Under A1's (1 − α_n)
scaling, the clip acts at 10⁻¹²/(1 − α_n) in unscaled units; the self-test asserts the clip is inactive
(min v > 10⁻¹²) in every evaluated case, so the equivalence is tested on unclipped quantities and the clip's
scale transformation is recorded.

## 5. Family and calibration

Gaussian; the NOVA-D construction with H ≡ 0: O ~ N(0, Σ_O); X_j = a_jᵀO + u_j, u_j ~ N(0, s²_j); T_k = b_kᵀO +
Σ_j λ_kj u_j + ε_k. G_j = s²_j Σ_k λ²_kj; the G profile G_1 = 1, G_2 = 1 − Γ_rel, G_j = (1 − Γ_rel) r^{j−2}; equal
marginal target variances; irreducible-noise share ψ = 0.25; r and σ²_ε solved so that A hits its target exactly
(`nova_d_calib.build` with ρ_pop = 0). Structure fixed per condition from the study's own salt.

**Conditional variances (new axis, §6).** *Equal*: s²_j = 0.5 for all j (as NOVA-D), under which the floor cannot
reorder at λ = 0 but can once λ > 0. *Heterogeneous*: s²_j drawn once per condition from log-uniform[0.25, 0.75]
with mean 0.5, under which the floor reorders even at λ = 0. Both levels are run; the heterogeneous level is where
the variance-floor interaction (review of 7 Nov) is expected to be consequential.

## 6. Grid — frozen in full; every condition reported

| Axis | Values | Definition |
|---|---|---|
| **Dimension ratio** d / n_train | n_train ∈ {102, 410, 820, 1640, 3280} | **numerator d = 60**, the dimension of the standardised matrix being shrunk (d_O + d_T + p); **denominator n_train**, the number of rows the selector sees — a direct axis, not derived from folds. α_n(frozen) = 0.3704, 0.1277, 0.0682, **0.05 (floor binds, n_train ≥ 1140)**, 0.05 |
| **Observed-block condition number** κ(Σ_OO) | φ ∈ {0.0, 0.5, 0.8, 0.95} in AR(1) Σ_OO | **population** κ(Σ_OO) = λ_max/λ_min of the 20 × 20 AR(1) matrix, computed exactly and reported per φ; the **estimated** κ(S_OO) per replicate is recorded as a diagnostic and never used for selection |
| Gaussian top-two gap Γ_rel | 0.02, 0.05, 0.20 | near-tie to clear |
| Assay utility A | 0.05, 0.15, 0.20 | as NOVA-D |
| Conditional variances | equal, heterogeneous | §5 |
| **Conditions** | 5 × 4 × 3 × 3 × 2 = **360** | |
| Replicates | 1,000 per condition (360,000 training samples; ~1 ms each) | |

Feasibility: no kernel ridge, no evaluation rows; the whole grid is minutes, not days.

## 7. Designated conditions (prospective)

- **D1 (H1):** Γ_rel = 0.02, A = 0.15, φ = 0.5, equal variances, n_train = 410 — the operative-α analogue of
  NOVA-D's primary, where the frozen rule's population argmax is known to differ from the Gaussian argmax. The
  self-test (§11) must confirm the flip for D1's own structure before the freeze; if D1's structure does not flip,
  the designation moves to the nearest flipping condition by a prespecified order (§13).
- **D2 (H2):** Γ_rel = 0.20, A = 0.15, φ = 0.5, equal variances, n_train = 410 — no population flip; the
  frozen rule's extra shrinkage may lower variance.

## 8. Objective A — exact, before any run

For every condition and every arm with a data-independent regularisation (A0, A1, A3; A2 at λ_c ∈ {0, matched};
A2′), the population argmax is computed from the exact Σ at the arm's exact α_n or λ, for every n_train in the
grid. Reported as a table; no simulation involved. A4's and A2's data-adaptive settings have no fixed population
argmax; their limits as n_train → ∞ (intensity → 0) are the Gaussian argmax and are stated, not simulated.

## 9. Procedure per replicate

Draw n_train rows from the condition's DGP (seed = SHA-256 of salt | condition | replicate); standardise; run
A0 (frozen call), A1, A2, A2′, A3, A4; record each pick, each regret (exact), the A1–A0 score discrepancy, the
estimated κ(S_OO), and each arm's chosen regularisation (λ_c^GCV, LW intensity).

## 10. Decision rules, fixed before any run

Let Δ_k = mean regret(A0) − mean regret(A_k) at a condition, with a paired percentile-bootstrap 95 % interval over
replicates (B = 2,000, seeded).

| Claim | Condition | Required evidence |
|---|---|---|
| **H1 — decoupling lowers regret where A0 is inconsistent** | D1 | Δ_2 > 0 with interval excluding 0, **and** Δ_2 ≥ 0.10 × mean regret(A0) (relevance margin; §13) |
| **H2 — no harm where A0 is consistent** | D2 | lower bound of Δ_2 ≥ −0.10 × mean regret(A0) |
| Floor contribution | n_train ≥ 1640 | Δ_3 reported with interval; it is zero by construction for n_train < 1140 |
| Frobenius versus decision | all | Spearman(mean Frobenius error, mean regret) across arms within condition — descriptive |

Both H1 and H2 are required for "a genuinely effective alternative"; H1 alone establishes a regime, not a method.
The full grid is descriptive; the grid maximum of Δ_2 is labelled upward-biased. Significance is not a criterion
(replicates are unlimited): the margin carries the claim.

## 11. Self-test gates (all must pass before `freeze`; the run refuses otherwise)

1. `p5_gate3a.py` hash `4600e2f3…` verified before import; A0's picks are obtained by calling it, not a copy.
2. **A1 ≡ A0 to 10⁻¹²** (relative, on C, v and scores) at **every fold-specific α_n of every n_train in the grid**,
   on every condition's exact Σ **and** on 1,000 random sample covariances of the grid's dimensions; zero argmax
   mismatches. The clip is asserted inactive in every case; α_n < 1 asserted everywhere.
3. Objective A table (§8) computed for every arm at its exact setting and every n_train; D1 confirmed to flip
   and D2 confirmed not to, for their own structures.
4. Calibration: A verified to 10⁻⁹ against its target for all 360 conditions; the heterogeneous-variance level's
   s²_j recorded.
5. Sampler check (owner's NOVA-D safeguard): on six designated conditions, 400,000 rows, G_j, R_∅ and A
   compared with closed forms by batch-means z-tests (|z| < 4) and absolute caps; assert-only.
6. Exact-regret identity: for 10 random picks per designated condition, the closed-form regret equals
   (R_pick − R_{j*})/R_∅ computed from the analytic conditional-variance formulas to 10⁻¹².
7. A4's intensity reproduced against an independent implementation on 100 random matrices to 10⁻¹⁰.

## 12. Integrity

Everything in §4–§11 is hashed before the first run; no condition, arm, rule or margin is changed after any result;
every condition is reported; `run` refuses unless frozen and hashes match; one run per freeze; the archive stays
closed; A0 is retained as the comparator whatever the outcome. The study tests whether decoupling improves
finite-sample decision loss; it does not retune NOVA-S and produces no claim about nonlinear acquisition.

## 13. Open items for owner decision

1. The relevance margin 0.10 in §10 (a fraction of A0's own regret) — value and whether it applies to H2 as drafted.
2. The GCV rule for λ_c as the designated primary setting, versus λ_c = 0 (A2′) as primary with GCV secondary.
3. β target = mean of residual variances; admissible range [0, 1]; designated β = 0 — confirm or amend.
4. Heterogeneous conditional variances as log-uniform[0.25, 0.75] — confirm or specify otherwise.
5. n_train grid {102, 410, 820, 1640, 3280} and φ grid {0, 0.5, 0.8, 0.95}; 1,000 replicates.
6. The fallback order for D1 if its structure does not flip (proposal: same cell at Γ_rel = 0.02, next φ).
7. Whether a realised-loss secondary (fitted linear predictor on an independent sample) is wanted at D1 and D2
   as a check that exact regret tracks realised loss, or whether the exact outcome suffices.
