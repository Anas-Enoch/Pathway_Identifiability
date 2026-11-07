# G1-R — Covariance regularisation and finite-sample acquisition decisions: protocol **draft** v0.2 (7 Nov 2026)

**Status: DRAFT for owner review. Not frozen, not hashed; no run may occur under a draft.**
Author of record: Anas Enoch. Drafted by Claude (Fable 5.1). v0.1 (`00ccf557…`, committed `c761fe4`) retained.

**Change log (v0.1 → v0.2), from the owner's decisions of 7 Nov 2026.** §13 decisions 1–7 applied: complete
bootstrapped contrasts with the margin inside the statistic, and a prespecified sparse-regret rule (§10); GCV fully
specified — objective, weighting, intercept, penalty scaling, range, endpoints, tie-breaking (§4.2); β confirmed
(§4.3); heterogeneous variances given one complete construction (§5); grid and 1,000 replicates fixed, runtime claim
removed, no post-hoc extension (§6); D1/D2 designation by frozen seed and finite ordered fallback lists with exact
criteria only (§7); no realised-loss secondary (§3). Additional corrections: A3 differs from A0 only for
n_train > 1140 (§4.5); population tuning criteria for the data-adaptive arms are specified and evaluated, not
asserted (§8); Schur-complement versus residual-covariance forms distinguished and assigned (§4.1); Ledoit–Wolf
centering, normalisation and intensity = 1 handled (§4.6); Frobenius comparison removed; no outer folds in this
study — fold-size coverage retained only inside the inherited equivalence gate (§11).

## 1. Question and scope

> When does covariance shrinkage improve finite-sample acquisition decisions, and when does its distortion of the
> conditional covariance induce persistent ranking errors?

Three evidential levels, kept apart: **Objective A**, population ranking consistency, decided by exact calculation
(§8); **Objective B**, finite-sample decision regret, the simulation's only question; **robustness to
ill-conditioning**, regret as a function of the population condition number of the observed block.

**Scope: Gaussian only.** Z = [O, T, P] jointly Gaussian, V_j = G_j, and the oracle optimum j\* = argmax_j G_j is
the same for every arm. Not a test of nonlinear acquisition, not a NOVA-S rescue, nothing about biological data.
The 186-record archive is not used; no NOVA-D output is reused except the calibration machinery of
`nova_d_calib.py`. **The outcome isolates acquisition selection from predictor fitting; it does not measure the
risk of any fitted downstream predictor, and no such claim will be made.**

Established before this draft and not re-tested: the frozen rule's population argmax differs from the Gaussian
argmax in 27 of 72 NOVA-D conditions at operative α (`NOVA_D_EXACT_ORDERING_NOTE_v1.2.md`); frozen joint
shrinkage is exactly a ridge-conditioned Schur estimator (`g1_shrinkage_equivalence.py`, `05ec454`), for
0 ≤ α < 1 with μ = 1 and inactive clipping.

## 2. Hypotheses

**H1 (primary; Objective B where A0 is inconsistent).** At D1, where the frozen rule's population argmax at its
operative α is not the Gaussian argmax, A2 has lower expected regret than A0 by at least the margin of §10.

**H2 (no-harm; Objective B where A0 is consistent).** At D2, where the frozen rule's population argmax equals the
Gaussian argmax and the Gaussian gap is wide, A2's expected regret does not exceed A0's by more than the margin.
The frozen rule's heavier shrinkage may win at D2; that outcome is reportable as such.

Both are required for "a genuinely effective alternative". Superiority of decoupling is a hypothesis to test.

## 3. Primary outcome

For a pick ĵ from a training sample, the **acquisition decision loss** is the exact population regret

    r = (G_{j*} − G_ĵ) / R_∅ ,

from the condition's population covariance (R_S = R_∅ − V_S; V_j = G_j here). No evaluation sample, no fitted
predictor, no cross-fitting, no realised-loss secondary. Per condition and arm: mean regret; P(ĵ = j\*); the pick
distribution; the arm's chosen regularisation (λ_c^GCV, δ_LW); degenerate-tie count (§4.6).

## 4. Arms

All arms receive the same training sample Z (n_train × 60), standardised as the frozen rule standardises
(column mean removed, column sample SD with ddof = 1 divided out), and use the frozen scoring form
score_j = Σ_k C²_kj / max(v_j, 10⁻¹²), argmax with ties to the lowest index (NumPy's `argmax`). ε = 10⁻⁶ is the
frozen ridge on the observed block; d = 60. Covariances are sample covariances with ddof = 1, S = ZᵀZ/(n_train − 1)
on centred Z, as `np.cov` computes them.

### 4.1 Two conditional forms, distinguished
For a ridge λ ≥ 0 on the observed block:
- **Schur form:** C^S(λ) = S_TP − S_TO(S_OO + λI)⁻¹S_OP; v^S_j(λ) = [S_PP − S_PO(S_OO + λI)⁻¹S_OP]_jj.
- **Residual form:** B_Y = (S_OO + λI)⁻¹S_OY for Y ∈ {T, P}; residuals e_Y = Y − O B_Y; C^R(λ) = cov(e_T, e_P),
  v^R_j(λ) = var(e_P,j), sample covariances with ddof = 1.
They coincide at λ = 0 and differ for λ > 0 (the residual form carries B_TᵀS_OO B_P − B_TᵀS_OP − S_TO B_P + S_TP
rather than S_TP − S_TO B_P). Joint shrinkage induces the **Schur** form (A1). "Condition first" means the
**residual** form (A2). The two are separated by arm A2-S.

### 4.2 Arm definitions

| Arm | Definition | Form | Role |
|---|---|---|---|
| **A0 — frozen G1** | `p5_gate3a.g1_pick` (sha256 `4600e2f3…`), hash-verified, called unchanged; α_n = min(1, max(0.05, d/(d + n_train))) | Schur | comparator |
| **A1 — equivalence control** | λ = (α_n + ε)/(1 − α_n); floor α_n/(1 − α_n) added to v^S_j; scores × (1 − α_n). **Gate** (§11): equals A0 to 10⁻¹². | Schur | certifies the decomposition; no separate outcome |
| **A2 — decoupled (primary)** | λ_c by GCV (§4.2a) on the residual form; β = 0 (§4.3) | residual | primary |
| **A2-S** | λ_c by the same GCV; Schur form; β = 0 | Schur | secondary: separates form from tuning |
| **A2′** | λ_c = 0 (OLS conditioning; n_train > 20 always); β = 0 | both coincide | secondary: the consistent, unregularised estimator |
| **A3 — vanishing floor** | frozen rule with α_n = min(1, d/(d + n_train)); differs from A0 **only for n_train > 1140** (equal at 1140, where d/(d + 1140) = 0.05) | Schur | control |
| **A4 — Ledoit–Wolf** | §4.6; then the frozen Schur conditioning with ridge ε | Schur | control: Frobenius-optimal joint shrinkage versus decision loss |

**4.2a GCV, fully specified.** Multivariate ridge regression of Y = [T, P] (n_train × 40, centred columns, no
intercept) on O (n_train × 20, centred, no intercept), **one shared penalty** for all 40 responses (a conditional
covariance requires one conditioning of all of T and P on O; per-response penalties are inadmissible). Responses
are equally weighted (all columns have unit sample variance). Hat matrix H_γ = O(OᵀO + γI)⁻¹Oᵀ; objective

    GCV(γ) = (1/n_train) ‖Y − H_γ Y‖²_F / (1 − tr(H_γ)/n_train)² ,

the standard generalised cross-validation (leave-one-out approximation with the hat-matrix trace). **Penalty
scaling:** the covariance-form ridge λ_c relates to the regression penalty by γ = (n_train − 1)·λ_c, so that
(S_OO + λ_c I)⁻¹S_OY = (OᵀO + γI)⁻¹OᵀY exactly. **Search set:** λ_c ∈ {0} ∪ {10^{−6 + 8k/80}: k = 0,…,80}
(81 log-spaced points from 10⁻⁶ to 10², plus zero); GCV(0) is evaluated with H_0 the OLS hat matrix. **Minimiser:**
the λ_c with the smallest GCV; **ties** (|ΔGCV| < 10⁻¹² relative) resolved to the **smallest** λ_c. **Endpoints:**
if the minimiser is 10² the replicate is flagged `endpoint_hit` and the value is used as is; the search set is
never extended. Implemented in closed form from one eigendecomposition of OᵀO per replicate.

### 4.3 β, confirmed
β acts on the residual variances only: v_β,j = (1 − β)v_j + β·mean_j v_j; target the mean of the 20 residual
variances; admissible range [0, 1]; **β = 0 is the designated primary**; β ∈ {0.5, 1} are descriptive secondary
settings. Cross-covariances C are never shrunk. At β = 1 the score is Σ_k C²_kj, which does not tie generically.
No ranking invariance in β is claimed.

### 4.4 α = 1
The frozen min(1, ·) gives α = 1 only when d/(d + n_train) ≥ 1, impossible for n_train ≥ 1; asserted in every
condition; A1's λ is therefore finite everywhere.

### 4.5 Denominator clipping
max(v_j, 10⁻¹²) retained exactly in every arm. Under A1's (1 − α_n) scaling the clip acts at 10⁻¹²/(1 − α_n) in
unscaled units; the gate asserts min v > 10⁻¹² in every evaluated case so that equivalence is tested on
unclipped quantities, and the scale transformation is recorded.

### 4.6 Ledoit–Wolf, fully specified
On the standardised, centred Z: S_n = ZᵀZ/n_train (the 1/n convention of Ledoit & Wolf 2004, used only inside
the intensity formula; the factor n/(n − 1) cancels in C²/v); target μI with μ = tr(S_n)/d; intensity
δ = min(1, max(0, π̂/(n_train·γ̂))) with π̂ = Σ_{ij} (1/n_train)Σ_t (z_ti z_tj − S_n,ij)² and γ̂ = ‖S_n − μI‖²_F —
the Ledoit–Wolf (2004) estimator, as implemented in scikit-learn's `ledoit_wolf` with `assume_centered=True`,
against which the implementation is checked (§11.7). S_LW = (1 − δ)S_n + δμI, then the frozen Schur conditioning
with ridge ε and scoring. **δ = 1:** S_LW = μI, C = 0, all scores 0 after clipping; the pick resolves to index 0
by the tie rule, the replicate is flagged `degenerate_tie`, its regret is computed for that pick, and the
degenerate-tie rate is reported per arm and condition.

## 5. Family and calibration

Gaussian, H ≡ 0: O ~ N(0, Σ_O) with (Σ_O)_ab = φ^{|a−b|}; X_j = a_jᵀO + u_j, u_j ~ N(0, s²_j), Var(a_jᵀO) = 1 − s²_j;
T_k = b_kᵀO + Σ_j λ_kj u_j + ε_k. G_j = s²_j Σ_k λ²_kj with λ_{·j} = √(G_j/s²_j)·w_j, w_j a random unit vector;
G profile G_1 = 1, G_2 = 1 − Γ_rel, G_j = (1 − Γ_rel)r^{j−2}; equal marginal target variances via b_k; ψ = 0.25;
r and σ²_ε solved so that A hits its target exactly (`nova_d_calib.build`, ρ_pop = 0). Structure fixed per
condition from the study salt `G1R-V02|`.

**Conditional variances — one complete construction.** *Equal:* s²_j = 0.5 for all j. *Heterogeneous:* the fixed
geometric sequence s²_(k) = 0.25·3^{k/19}, k = 0,…,19 (0.2500, 0.2649, …, 0.7079, 0.7500; exact mean 0.4575; all
values inside [0.25, 0.75] by construction; no rescaling), assigned to candidates by a permutation drawn once per
condition from the structure seed, so that s²_j is not aligned with Gaussian rank. Candidate index equals Gaussian
rank in both levels.

## 6. Grid — frozen in full; every condition reported

| Axis | Values | Definition |
|---|---|---|
| **Dimension ratio** d/n_train | n_train ∈ {102, 410, 820, 1640, 3280} | numerator **d = 60** (dimension of the standardised matrix being shrunk); denominator **n_train** (rows the selector sees) — a direct axis; **there are no outer folds in this study**. Frozen α_n = 0.3704, 0.1277, 0.0682, 0.0500 (floor), 0.0500 (floor); A3's α = 0.0353 and 0.0180 at the two largest |
| **Condition number** κ(Σ_OO) | φ ∈ {0, 0.5, 0.8, 0.95} | **population** κ = λ_max/λ_min of the 20 × 20 AR(1) matrix: **1.00, 8.63, 64.7, 569.2**; the **estimated** κ(S_OO) per replicate is recorded as a diagnostic only |
| Gaussian top-two gap Γ_rel | 0.02, 0.05, 0.20 | |
| Assay utility A | 0.05, 0.15, 0.20 | |
| Conditional variances | equal, heterogeneous | §5 |
| **Conditions** | 5 × 4 × 3 × 3 × 2 = **360** | |
| **Replicates** | **1,000 per condition, fixed**; not extended after any result is seen | |

## 7. Designated conditions — by frozen seed and exact criteria only

Structural seeds are fixed by the salt before any evaluation. Each designation is the **first** entry of a finite
ordered list whose structure satisfies the exact population criterion at the listed n_train; if none qualifies,
the study **stops** and reports that, with no discretionary substitution.

**D1 (H1) — criterion:** A0's population argmax at its operative α_n ≠ Gaussian argmax (exact, §8).
Ordered list (Γ_rel = 0.02, A = 0.15, equal variances, n_train = 410, then φ ∈ 0.5, 0.8, 0.0, 0.95; then the same
four at A = 0.20; then at A = 0.05). Twelve candidates.

**D2 (H2) — criterion:** A0's population argmax at its operative α_n = Gaussian argmax **and** Γ_rel = 0.20.
Ordered list (Γ_rel = 0.20, A = 0.15, equal variances, n_train = 410, then φ ∈ 0.5, 0.8, 0.0, 0.95; then A = 0.20;
then A = 0.05). Twelve candidates.

D1 is **not** assumed to flip because its axes resemble NOVA-D's primary: this study's signal construction and
structural seeds differ, and the criterion is evaluated for D1's own structure.

## 8. Objective A — exact, before any run

For every condition, every n_train, and every arm with a data-independent setting (A0, A1, A2′, A3; A2 and A2-S at
λ_c = 0 and at the A1-matched λ), the population argmax from the exact Σ.

**Data-adaptive arms — population tuning criteria, specified and evaluated.**
- *A2/A2-S (GCV).* The population analogue of GCV is the population prediction error of the ridge regression,
  PE(λ) = tr[(B_λ − B_0)ᵀ Σ_OO (B_λ − B_0)] with B_λ = (Σ_OO + λI)⁻¹Σ_OY. Evaluated on each condition's exact Σ
  over the search set; the minimiser is reported. (It is λ = 0 — PE is nondecreasing in λ, e.g. 0, 1.3 × 10⁻⁵,
  0.129, 9.19, 175 at λ = 0, 10⁻⁴, 10⁻², 10⁻¹, 1 on a test structure — because at the population level ridge adds
  bias and removes no variance; this is *evaluated* per condition, not asserted.) The population argmax of A2 at
  that minimiser is reported.
- *A4 (Ledoit–Wolf).* The population intensity at n_train is δ_pop(n) = min(1, π_pop/(n·γ_pop)) with
  π_pop = Σ_{ij} Var(Z_iZ_j) = Σ_{ij}(Σ_iiΣ_jj + Σ²_ij) (Gaussian fourth moments, exact) and γ_pop = ‖Σ_std − μI‖²_F;
  the population argmax of A4 at δ_pop(n) is reported for each n_train in the grid.

## 9. Procedure per replicate

Seed = SHA-256(salt | condition index | replicate). Draw n_train rows; standardise; run A0 (frozen call), A1, A2,
A2-S, A2′, A3, A4 and the β ∈ {0.5, 1} variants of A2; record each pick, its exact regret, the A1–A0 relative
discrepancy, λ_c^GCV and its `endpoint_hit` flag, δ_LW and its `degenerate_tie` flag, and κ(S_OO).

## 10. Decision rules, fixed before any run

Let r_0 and r_2 be the per-replicate regrets of A0 and A2 at a condition. Both contrasts are bootstrapped **as
complete statistics**, margin included, so that baseline uncertainty is inside the interval; no confidence bound is
compared with a separately estimated margin:

    H1 at D1:  Q1 = mean(0.9·r_0 − r_2) ;  supported iff the percentile-bootstrap 95 % lower bound of Q1 > 0.
    H2 at D2:  Q2 = mean(1.1·r_0 − r_2) ;  supported iff the percentile-bootstrap 95 % lower bound of Q2 ≥ 0.

B = 2,000 paired resamples of replicates, seeded. The factors 0.9 and 1.1 are the provisional 10 % relative
margins (§13.1).

**Sparse-regret rule (prespecified).** A condition–arm pair is *sparse* if fewer than 10 of the 1,000 replicates
have r > 0. If either arm of a contrast is sparse, the bootstrap is not used for that contrast: the exact
Clopper–Pearson 95 % upper bound on P(ĵ ≠ j\*) is reported for each arm, the hypothesis is recorded as **not
assessable by bootstrap** (neither supported nor refuted), and the reason is stated. A degenerate bootstrap cannot
establish the absence of rare harmful picks, and this record does not pretend otherwise.

**Other reports.** A3 − A0 regret with interval at n_train ∈ {1640, 3280} (identically zero below). A2-S − A2
(form effect) and A2′ − A2 (tuning effect), descriptive. The full grid is descriptive; the grid maximum of
mean(r_0 − r_2) is labelled upward-biased. Significance is never a criterion.

## 11. Self-test gates — all must pass before `freeze`; `run` refuses otherwise

1. `p5_gate3a.py` hash `4600e2f3…` verified before import; A0's picks come from calling it.
2. **A1 ≡ A0 to 10⁻¹²** (relative, on C, v, scores; the implementation gate is 10⁻¹², tighter than the 10⁻¹⁰ in
   the committed derivation script) at every α_n of every grid n_train **and**, for inherited coverage, at
   α = d/(d + m) for m ∈ {102, 103, 409, 410, 819, 820} — the NOVA-D fold sizes, included only so that the
   equivalence record covers them, **not** because this study cross-fits — on every condition's exact Σ and on
   1,000 random sample covariances of the grid's dimensions; zero argmax mismatches; clipping inactive and
   α_n < 1 asserted in every case.
3. Objective A tables (§8) for every arm and n_train; D1 and D2 selected from their lists by exact criteria and
   recorded with the list position used; stop if none qualifies.
4. Calibration: A to 10⁻⁹ for all 360 conditions; each heterogeneous condition's permutation of s²_j recorded.
5. Sampler: on six designated conditions, 400,000 rows, G_j, R_∅, A compared with closed forms by batch-means
   z-tests (|z| < 4) and absolute caps; assert-only.
6. Exact-regret identity: for 10 random picks per designated condition, (G_{j*} − G_ĵ)/R_∅ equals
   (R_pick − R_{j*})/R_∅ from the analytic conditional-variance formulas to 10⁻¹².
7. GCV: the closed-form implementation equals direct leave-one-out-free evaluation of the objective on 100 random
   problems to 10⁻¹⁰, and the matched-penalty identity (S_OO + λ_cI)⁻¹S_OY = (OᵀO + γI)⁻¹OᵀY holds to 10⁻¹².
   Ledoit–Wolf: δ equals scikit-learn's `ledoit_wolf(assume_centered=True)` on 100 random matrices to 10⁻¹⁰.
8. Schur and residual forms coincide at λ = 0 to 10⁻¹² and differ at λ > 0 on a test structure (recorded).

## 12. Integrity

Everything in §4–§11 hashed before the first run; no condition, arm, rule, margin or replicate count changed after
any result; every condition reported; `run` refuses unless frozen and hashes match; one run per freeze; the archive
stays closed; A0 retained as comparator whatever the outcome. The study produces no claim about nonlinear
acquisition and does not retune NOVA-S.

## 13. Open items for owner decision

1. Margins 0.9 / 1.1 kept provisionally per your decision; confirm before freeze.
2. The GCV specification of §4.2a (shared penalty; search set; smallest-λ tie rule) — confirm or amend.
3. The sparse-regret threshold (10 of 1,000) and the Clopper–Pearson reporting in its place.
4. The D1/D2 ordered lists of §7 (twelve candidates each, order as written).
5. Whether A2-S (Schur-form decoupled) is wanted as a secondary arm, as drafted, to separate form from tuning.
6. The study salt string `G1R-V02|` (to be re-salted at freeze if the version changes).
