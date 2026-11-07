# G1-R — Covariance regularisation and finite-sample acquisition decisions: protocol **draft** v0.4 (7 Nov 2026)

**Change log (v0.3 → v0.4), from the owner's code review of 7 Nov 2026.** All gates are now enforced, not recorded:
(1) the A1 gate compares scaled C, v and scores under a defined norm (max-abs relative to the max-abs of the frozen
quantity) at the correctly scaled clipping threshold; (2) replicate-level gates — mirror argmax = frozen pick, A1
within 10⁻¹² on C, v and s, clip inactive, α < 1, picks in range — are **assertions that stop the run**; (3) `freeze`
binds the passing self-test report to the current file hashes and version, records the report's checksum, and
refuses to overwrite an existing freeze; `run` and `analyze` re-verify the report's checksum; (4) complete shrinkage
(δ = 1) is distinguished from exact score ties; (5) the secondary contrasts A2-S vs A2 and A2′ vs A2 are implemented
with sparse handling; (6) the GCV objective is checked direct-versus-SVD across the whole search set with the selected
penalty verified; (7) the sampler moment check is implemented on six fixed conditions with fixed seeds; (8) a
replicate-level integrity checker (`g1r_integrity.py`) is required before `analyze`; (9) §11 states what has been
verified rather than "all passed". Implementation: `g1r_run.py` **v0.2.0**.

**Status: DRAFT. Not frozen; no study run; no freeze record exists.** The implementation's self-test passes every
gate in §11; the freeze is a separate, later act on the owner's instruction.
Author of record: Anas Enoch. Drafted by Claude (Fable 5.1). v0.1 (`c761fe4`), v0.2 (`bf7d6f6`), v0.3 (`897b85c`) retained.

**Change log (v0.2 → v0.3).** §13 of v0.2: margins 0.9/1.1 confirmed,
both hypotheses required, margins inside the bootstrapped statistics (§10); GCV rule accepted (§4.2a); sparse-regret
rule accepted, Clopper–Pearson bound specified as one-sided 95 % on misselection probability and applied to every
contrast (§10); ordered designation lists accepted (§7); A2-S retained (§4.2); salt `G1R-V02|` fixed and not
re-salted on document versioning. Three substantive corrections: the GCV population calculation is relabelled an
**infinite-data population benchmark**, evaluated in population-standardised coordinates, and is not used to
certify the finite-n data-adaptive arm's ranking (§8); the Ledoit–Wolf population intensity is relabelled a
**deterministic diagnostic proxy**, not a finite-n expectation (§8); the Ledoit–Wolf covariance-scaling choice is
made explicit — intensity from the 1/n formula, applied to the ddof = 1 covariance, then the frozen ridge ε (§4.6).
Reporting: bootstrap lower bound defined (§10); "significance is never a criterion" replaced (§10); heterogeneous
mean stated as ≈ 0.457484 (§5); structural-seed serialisation and condition enumeration made explicit (§6, §9).
The owner's residual identity C^R = C^S − λB_TᵀB_P, v^R_j = v^S_j − λ‖B_P,j‖² is a self-test gate (§11.8).

## 1. Question and scope

> When does covariance shrinkage improve finite-sample acquisition decisions, and when does its distortion of the
> conditional covariance induce persistent ranking errors?

**Objective A** — population ranking consistency, by exact calculation (§8). **Objective B** — finite-sample
decision regret, the simulation's only question. **Ill-conditioning** — regret against the population condition
number of the observed block. Kept apart throughout.

**Scope: Gaussian only.** Z = [O, T, P] jointly Gaussian; V_j = G_j; the oracle optimum j\* = argmax_j G_j is the
same for every arm. Not a test of nonlinear acquisition, not a NOVA-S rescue, nothing about biological data. The
186-record archive is not used. **The outcome isolates acquisition selection from predictor fitting; it does not
measure the risk of any fitted downstream predictor, and no such claim will be made.**

Established before this draft, not re-tested: frozen joint shrinkage is exactly a ridge-conditioned Schur estimator
for 0 ≤ α < 1, μ = 1, inactive clipping (`g1_shrinkage_equivalence.py`, `05ec454`); the frozen rule's population
argmax differs from the Gaussian argmax in 27 of 72 NOVA-D conditions at operative α (`NOVA_D_EXACT_ORDERING_NOTE_v1.2.md`).

## 2. Hypotheses

**H1 (primary; Objective B where A0 is inconsistent).** At D1, where the frozen rule's population argmax at its
operative α is not the Gaussian argmax, A2 has lower expected regret than A0 by at least the margin (§10).

**H2 (no-harm; Objective B where A0 is consistent).** At D2, where the frozen rule's population argmax equals the
Gaussian argmax and the Gaussian gap is wide, A2's expected regret does not exceed A0's by more than the margin.
The frozen rule's heavier shrinkage may win at D2; that outcome is reportable as such.

Both are required. Superiority of decoupling is a hypothesis to test.

## 3. Primary outcome

Exact population regret of each pick, r = (G_{j\*} − G_ĵ)/R_∅, from the condition's population covariance
(R_S = R_∅ − V_S; V_j = G_j). No evaluation sample, no fitted predictor, no cross-fitting, no outer folds, no
realised-loss secondary. Per condition and arm: mean regret; P(ĵ = j\*); the one-sided 95 % Clopper–Pearson upper
bound on P(ĵ ≠ j\*); λ_c^GCV with endpoint flag; δ_LW with degenerate-tie flag; estimated κ(S_OO) (diagnostic).

## 4. Arms

All arms receive the same training sample (n_train × 60), centred and scaled by column sample SD (ddof = 1), as the
frozen rule standardises; S = ZᵀZ/(n_train − 1) on the centred, scaled Z (`np.cov`, ddof = 1). Scoring form
score_j = Σ_k C²_kj / max(v_j, 10⁻¹²); argmax with ties to the lowest index. ε = 10⁻⁶; d = 60.

### 4.1 Two conditional forms
For ridge λ ≥ 0 on the observed block, with B_Y = (S_OO + λI)⁻¹S_OY:
- **Schur:** C^S = S_TP − S_TO B_P; v^S_j = [S_PP − S_PO B_P]_jj.
- **Residual:** C^R = cov(T − OB_T, P − OB_P) = S_TP − S_TO B_P − B_TᵀS_OP + B_TᵀS_OO B_P; v^R_j likewise.
- **Exact relation (owner, 7 Nov):** C^R = C^S − λB_TᵀB_P; v^R_j = v^S_j − λ‖B_P,j‖². Coincide at λ = 0; differ for
  λ > 0. Joint shrinkage induces the Schur form (A1). "Condition first" is the residual form (A2).

### 4.2 Arm definitions

| Arm | Definition | Form | Role |
|---|---|---|---|
| **A0** | `p5_gate3a.g1_pick` (sha256 `4600e2f3…`), hash-verified, called unchanged; α_n = min(1, max(0.05, d/(d + n_train))) | Schur | comparator |
| **A1** | λ = (α_n + ε)/(1 − α_n); floor α_n/(1 − α_n) on v^S; C, v and scores × (1 − α_n). **Gate:** scaled C, scaled v and scores each within 10⁻¹² of the frozen quantities (max-abs relative to the frozen max-abs), argmax equal, clip inactive at the scaled threshold 10⁻¹²/(1 − α_n) | Schur | equivalence control; no separate outcome |
| **A2 (primary)** | λ_c by GCV (§4.2a); residual form; β = 0 | residual | primary |
| **A2-S** | same λ_c; Schur form; β = 0 | Schur | separates form from tuning |
| **A2′** | λ_c = 0; β = 0 (forms coincide) | both | consistent, unregularised |
| **A3** | frozen algebra at α = min(1, d/(d + n_train)), no floor; differs from A0 **only for n_train > 1140** | Schur | floor control |
| **A4** | Ledoit–Wolf (§4.6), then frozen Schur conditioning with ridge ε | Schur | Frobenius-optimal joint shrinkage |
| A2 β = 0.5, β = 1 | A2 with residual-variance stabilisation (§4.3) | residual | descriptive secondary |

A0's scores are not returned by `g1_pick`; the gate compares A1 with a verified mirror of the frozen algebra. **In every
replicate the following are assertions — a failure raises and stops the run, and nothing is recorded for that
replicate:** mirror argmax = `g1_pick`; A1 within 10⁻¹² on scaled C, scaled v and scores; clip inactive; α_n < 1 and
α_A3 < 1; every pick in [0, 20). A3 necessarily uses the mirror (its α is not the frozen one).

**4.2a GCV.** Multivariate ridge of Y = [T, P] (n × 40, centred, no intercept) on O (n × 20, centred, no intercept),
**one shared penalty** γ for all 40 responses (a conditional covariance requires one conditioning of T and P on O),
equal response weights (unit sample variances). GCV(γ) = (RSS(γ)/n)/(1 − tr H_γ/n)², H_γ = O(OᵀO + γI)⁻¹Oᵀ, computed
in closed form from the thin SVD of O: RSS(γ) = ‖Y − WWᵀY‖² + Σ_i (γ/(e_i + γ))²‖(WᵀY)_i‖². **Penalty scaling**
γ = (n_train − 1)·λ_c, so that (S_OO + λ_cI)⁻¹S_OY = (OᵀO + γI)⁻¹OᵀY exactly (gate §11.7). **Search set**
Λ = {0} ∪ {10^{−6+8k/80}: k = 0..80}. **Minimiser:** smallest GCV; a candidate replaces the incumbent only if it
improves by more than 10⁻¹² relative, so **ties resolve to the smallest λ_c**. **Endpoint:** λ_c = 10² flags
`endpoint_hit`; Λ is never extended.

### 4.3 β
v_β,j = (1 − β)v_j + β·mean_j v_j; target the mean of the 20 residual variances; range [0, 1]; **β = 0 designated**;
β ∈ {0.5, 1} descriptive. C never shrunk. No ranking invariance in β is claimed.

### 4.4 α = 1
Impossible for n_train ≥ 1; asserted in every replicate; A1's λ finite everywhere.

### 4.5 Clipping
max(v_j, 10⁻¹²) retained in every arm; under A1's scaling the clip acts at 10⁻¹²/(1 − α_n) unscaled; the gate asserts
the clip inactive in every evaluated case (it was, in all 9,520 evaluations at drafting).

### 4.6 Ledoit–Wolf — explicit scaling choice
Intensity δ from the Ledoit–Wolf (2004) formula on the centred, scaled Z with S_n = ZᵀZ/n and target
μ_n I, μ_n = tr S_n/d — as scikit-learn's `ledoit_wolf(assume_centered=True)`, against which it is checked to
10⁻¹⁰ (§11.7). **Scaling:** δ is then applied to the **ddof = 1** covariance, S_LW = (1 − δ)S + δ(tr S/d)I, so that
the frozen ridge ε acts on A0's scale (replacing S by cS would change the effective ridge to ε/c; this choice fixes
c = 1). **δ = 1 (complete shrinkage):** S_LW = μI, C = 0, all scores 0 after clipping; the pick resolves to index 0 by the tie
rule, the replicate is flagged `complete_shrinkage`, and its regret is computed for that pick. Separately, an **exact**
top tie (more than one index attaining the maximum score, exact equality, not `allclose`) is flagged `A4_exact_top_tie`.
Both rates are reported per condition.

## 5. Family and calibration

Gaussian, H ≡ 0: O ~ N(0, Σ_O), (Σ_O)_ab = φ^{|a−b|}; X_j = a_jᵀO + u_j, u_j ~ N(0, s²_j), Var(a_jᵀO) = 1 − s²_j;
T_k = b_kᵀO + Σ_j λ_kj u_j + ε_k. λ_{·j} = √(G_j/s²_j)·w_j (w_j a random unit vector), so G_j = s²_jΣ_kλ²_kj exactly.
G profile G_1 = 1, G_2 = 1 − Γ_rel, G_j = (1 − Γ_rel)r^{j−2}; candidate index = Gaussian rank. Equal marginal target
variances via b_k (c = 2·max_k Var(T_k|O)); ψ = 0.25; r and σ²_ε from `nova_d_calib.build` with ρ_pop = 0 so that
A is exact (gate §11.4: 2.8 × 10⁻¹⁷ at drafting).

**Conditional variances.** *Equal:* s²_j = 0.5. *Heterogeneous:* the fixed sequence s²_(k) = 0.25·3^{k/19},
k = 0..19 — 0.2500, 0.2649, …, 0.7079, 0.7500; mean ≈ 0.457484; all values in [0.25, 0.75] by construction; no
rescaling — assigned to candidates by a permutation drawn once per condition from the structure seed.

## 6. Grid, enumeration and seeds — frozen in full; every condition reported

| Axis | Values | Definition |
|---|---|---|
| Dimension ratio d/n_train | n_train ∈ {102, 410, 820, 1640, 3280} | d = 60 (dimension of the matrix shrunk); n_train = rows the selector sees; **no outer folds in this study**. Frozen α_n = 0.3704, 0.1277, 0.0682, 0.0500, 0.0500; A3's α = 0.0353, 0.0180 at the two largest |
| Condition number κ(Σ_OO) | φ ∈ {0, 0.5, 0.8, 0.95} | population κ = **1.00, 8.63, 64.7, 569.2** (exact); κ(S_OO) per replicate is a diagnostic |
| Γ_rel | 0.02, 0.05, 0.20 | |
| A | 0.05, 0.15, 0.20 | |
| Conditional variances | equal (0), heterogeneous (1) | §5 |
| **Conditions** | 5 × 4 × 3 × 3 × 2 = **360** | |
| **Replicates** | **1,000 per condition, fixed; never extended after any result** | |

**Enumeration (explicit).** Nested loops outer→inner: n_train, φ, Γ_rel, A, hetero, in the orders listed; idx = 0..359.
**Condition string:** `n_train=<int>|phi=<float>|gamma=<float>|A=<float>|hetero=<0|1>` with Python's default float
repr (e.g. `n_train=410|phi=0.8|gamma=0.02|A=0.15|hetero=0`).
**Seeds:** 60-bit integers from the first 15 hex digits of SHA-256(UTF-8 of `G1R-V02|` + parts joined by `|`):
structure seed = SHA-256(`G1R-V02|structure|<condition string>`); replicate seed =
SHA-256(`G1R-V02|rep|<condition string>|<replicate int>`). Implemented in `g1r_run.seed_for`.

## 7. Designated conditions — frozen lists, exact criteria, no discretion

Ordered lists (A outer, φ inner): for A ∈ (0.15, 0.20, 0.05) and φ ∈ (0.5, 0.8, 0.0, 0.95), with equal variances
and n_train = 410 — D1 list at Γ_rel = 0.02, D2 list at Γ_rel = 0.20; twelve candidates each.
**D1 criterion:** A0's population argmax at its operative α_n ≠ Gaussian argmax (exact). **D2 criterion:** equality.
The designation is the first qualifying entry; if none qualifies the study stops.

**Result of the frozen procedure at drafting (population calculation, not a simulation outcome):**
D1 = **condition 110** (Γ_rel 0.02, A 0.15, φ 0.8, equal, n_train 410) — list position **1**; position 0 (φ = 0.5,
the axes of NOVA-D's primary) does **not** flip for this study's structure, as the owner anticipated. Of the twelve
D1 candidates, exactly four flip — 110, 128, 108, 126 — all with φ ∈ {0.8, 0.95}: at n_train = 410 the frozen
rule's population inconsistency appears only in the ill-conditioned observed blocks. D2 = **condition 104**
(Γ_rel 0.20, A 0.15, φ 0.5) — position 0; all twelve D2 candidates satisfy the criterion.

## 8. Objective A — exact, before any run

For every condition and n_train: the population argmax of A0, A3, A2′, and of A2 and A2-S at the A1-matched λ, from
the exact standardised Σ. Recorded in `selftest.json["objective_A"]`.

**Data-adaptive arms — benchmarks, labelled as such.**
- *GCV — infinite-data population benchmark.* PE(λ) = tr[(B_λ − B_0)ᵀΣ_OO(B_λ − B_0)] in population-standardised
  coordinates, over Λ; its minimiser is reported per condition (zero everywhere at drafting — ridge adds bias and
  removes no variance at infinite n). **This is not the finite-n expected GCV objective, which includes estimation
  variability, and it is not used to certify A2's finite-n ranking.**
- *Ledoit–Wolf — deterministic diagnostic proxy.* δ_pop(n) = min(1, π_pop/(n·γ_pop)) with π_pop = Σ_ij(Σ_iiΣ_jj + Σ²_ij)
  (Gaussian fourth moments) and γ_pop = ‖Σ_std − μI‖²_F; A4's argmax at δ_pop(n) is reported. **This is a proxy, not
  the expectation of the estimated, clipped intensity after sample centring and standardisation**; it gives no
  exact finite-n population ranking of A4.

## 9. Procedure per replicate

Seed as §6. Draw n_train rows; centre and scale; S = `np.cov(ddof=1)`; A0 by the frozen call; mirror and A1 gates as
**assertions** (§4.2); λ_c by GCV; A2, A2-S, A2′, A2 at β = 0.5 and 1, A3, A4; record every pick, every exact regret,
the three A1 discrepancies, λ_c, `endpoint_hit`, δ_LW, `complete_shrinkage`, `A4_exact_top_tie`, κ(S_OO). A `GateError`
propagates out of `run`; the run halts and the partial log is preserved for diagnosis.

## 10. Decision rules, fixed before any run

Let r_0, r_k be per-replicate regrets of A0 and arm k at a condition. Contrasts are **complete bootstrapped
statistics with the margin inside**:

    Q1 = mean(0.9·r_0 − r_2),   Q2 = mean(1.1·r_0 − r_2),   plain = mean(r_0 − r_k).

**Bootstrap:** B = 2,000 paired resamples of the 1,000 replicates (with replacement, seeded); the **lower bound is the
5th percentile** of the resampled means (one-sided 95 %); the two-sided 95 % interval is the 2.5th–97.5th percentiles.
**H1 at D1:** supported iff LB₅(Q1) > 0. **H2 at D2:** supported iff LB₅(Q2) ≥ 0. The decisions require **both** a
relevance margin (inside the statistic) and an uncertainty bound (the quantile); neither alone suffices.

**Sparse-regret rule — applied to every contrast, primary and secondary.** An arm is *sparse* at a condition if fewer
than 10 of 1,000 replicates have r > 0. If A0 or arm k is sparse, that contrast is **not assessable by bootstrap**
(neither supported nor refuted); the one-sided 95 % Clopper–Pearson upper bound on P(ĵ ≠ j\*) — the Beta(k+1, n−k)
0.95-quantile with k misselections of n — is reported for each arm. It bounds misselection probability, not expected
regret, and the record says so.

**Secondary contrasts (implemented, with the same sparse handling):** A2-S versus A2 — the *form effect*, mean(r_A2 −
r_A2-S) with bootstrap interval, not assessable if either arm is sparse; A2′ versus A2 — the *tuning effect*,
mean(r_A2 − r_A2′), likewise. A3 − A0 at n_train ∈ {1640, 3280} (identically zero below). The grid maximum of plain
A2 gain is labelled upward-biased.

## 11. Self-test gates — what has been verified (`g1r_out/selftest.json`, v0.2.0)

Every gate below is implemented and was run at drafting; `ok: true` means all of them passed. It is a pass of the
implemented gate set and becomes freeze-readiness only on the owner's instruction to freeze.

1. `p5_gate3a.py` hash verified before import; A0 by call. Mirror argmax = `g1_pick` on 200 random samples: 0 mismatches.
2. **A1 ≡ A0** on 1,360 matrices (all 360 exact Σ + 1,000 random sample covariances) × 7 α (five grid α_n and, for
   inherited coverage only, d/(d+m), m ∈ {102, 103, 409, 410, 819, 820} — this study does not cross-fit): max relative
   discrepancy **C 3.4 × 10⁻¹⁵, v 1.2 × 10⁻¹⁵, scores 2.5 × 10⁻¹⁵**; 0 argmax mismatches; clip active in 0 of 9,520
   (scaled thresholds 1.05–1.59 × 10⁻¹² unscaled).
3. Objective A tables for all 360 conditions; designation by the frozen lists: D1 = 110 (position 1), D2 = 104
   (position 0).
4. Calibration: max |A_oracle − A_target| = 2.8 × 10⁻¹⁷ over 360.
5. **Sampler moment check** (owner's safeguard) on six fixed conditions with fixed seeds — D1 (110), D2 (104), and
   four chosen to cover both variance constructions and both conditioning extremes: 1 (n 102, φ 0, Γ 0.02, A 0.05,
   heterogeneous), 359 (n 3280, φ 0.95, Γ 0.20, A 0.20, heterogeneous), 204 (n 820, φ 0.95, Γ 0.05, A 0.05, equal),
   225 (n 1640, φ 0, Γ 0.05, A 0.15, heterogeneous); 400,000 rows each; sampled G_j, R_∅ and A against the analytic
   values by batch-means z (20 batches, |z| < 4) and absolute caps (0.02): all six pass — max |z| 1.43, max G error
   0.0071 (scale G_1 = 1), max A error 0.0005, max marginal-variance spread 0.0055; sampled argmax G = j\* in all six.
6. Exact-regret identity via the covariance route: 3.9 × 10⁻¹⁶ on 20 random picks at D1 and D2.
7. **GCV objective**: direct hat-matrix versus SVD closed form over the whole 82-point search set on 20 random
   problems: max relative difference 2.4 × 10⁻¹⁵; selected penalty identical in 20 of 20. Penalty-scaling identity
   (S_OO + λ_cI)⁻¹S_OY = (OᵀO + γI)⁻¹OᵀY: 9.7 × 10⁻¹⁷.
8. Ledoit–Wolf δ versus an independent direct-definition formula: 0 (exact); versus scikit-learn where available:
   0. Residual identity C^R = C^S − λB_TᵀB_P, v^R = v^S − λ‖B_P‖²: 2.2 × 10⁻¹⁶; forms coincide at λ = 0 (2.2 × 10⁻¹⁶)
   and differ at λ = 0.3. Clopper–Pearson at k = 0 versus 1 − 0.05^{1/1000}: < 10⁻¹².
9. Replicate gates exercised as assertions: 5 replicates at D1 and 5 at D2 passed.
10. Timing (measured, drafting container): 1.5 ms per replicate at n_train = 102; 13 ms at 3,280.

**Integrity checker** `g1r_integrity.py` (read-only): rejects unparsable or malformed records, duplicate or missing
(condition, replicate) keys, picks outside [0, 20), non-finite values, regrets inconsistent with the recorded pick
(recomputed exactly from the frozen structure), and records whose A1 discrepancies exceed 10⁻¹²; re-verifies the
freeze hashes and the self-test checksum. Exercised on fabricated structure-only records for nine failure modes
(complete grid PASS; one missing, one duplicate, invalid pick, inconsistent regret, A1 failure, non-finite value,
malformed record, corrupt line each FAIL with the correct diagnosis). `analyze` refuses without INTEGRITY PASS.

## 12. Integrity

Everything in §4–§11 hashed at freeze; no condition, arm, rule, margin, list or replicate count changed after any
result; every condition reported; the archive closed; A0 retained as comparator whatever the outcome; no claim about
nonlinear acquisition; NOVA-S not retuned. **Freeze binding:** `freeze` refuses if a freeze record already exists, if
any output exists, if the self-test report is absent or not a pass, or if the report's version or recorded file hashes
differ from the current `g1r_run.py`, `nova_d_calib.py` and `p5_gate3a.py`; it records the protocol hash, the three
file hashes, the self-test report's checksum, the designations and the moment-check conditions. `run` and `analyze`
re-verify all of these, including the report checksum. One run per freeze. `analyze` additionally refuses without an
INTEGRITY PASS from `g1r_integrity.py`. Implementation: `g1r_run.py` v0.2.0 (selftest | freeze | run [--cond a:b] |
analyze) with `g1r_integrity.py`, `nova_d_calib.py` and the frozen `p5_gate3a.py`.

## 13. Open items before freeze

None from the owner's review of 7 Nov 2026 remain unimplemented. The freeze awaits the owner's explicit instruction;
on it, the sequence is: `python3 -I g1r_run.py selftest` (regenerates the report against the current files) →
`python3 -I g1r_run.py freeze` → commit `g1r_out/FREEZE.json` and `g1r_out/selftest.json` → `run` → `g1r_integrity.py`
→ `analyze`.
