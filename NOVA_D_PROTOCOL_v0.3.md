# NOVA-D — Detectability of Nonlinear Acquisition Advantage: synthetic protocol v0.3 (5 Nov 2026)

**Status: ready to freeze on one owner confirmation (§13, item 1). No simulation may run before the freeze.**
Author of record: Anas Enoch. Drafted by Claude (Fable 5.1) under the owner's decisions of 5 Nov 2026 (two
memoranda). Implementation: `nova_d_run.py` v0.1.0 with `nova_d_calib.py`; both hashed at freeze.

## 1. Scope and research question

> How much of a known nonlinear acquisition advantage can NOVA-S recover at finite sample sizes, and does
> near-tie gating improve that recovery or obstruct it?

Three questions, answered separately and never sharing a success label:
1. **Existence** — under what joint distributions does a meaningful nonlinear acquisition advantage exist?
2. **Recoverability** — can the frozen NOVA-S estimator recover a positive fraction of it?
3. **Selection** — does near-tie screening increase useful departures, or merely estimation-driven disagreement?

A negative answer to (2) or (3) is a substantive result.

**What this is not.** Not an external validation. Synthetic success justifies renewed external testing; it does not
validate NOVA-S on biological data. A successful synthetic regime yields a precise biological prediction, not a
validation claim. NOVA-D can show that a mechanism generates useful nonlinear acquisition opportunities; it cannot
show that the mechanism is common in real metabolomics cohorts.

**Two inference levels, never conflated.** *Monte Carlo uncertainty* is variation across independently generated
training datasets under a fixed synthetic DGP; the bootstrap over replicates estimates it. *External predictive
uncertainty* is variation across independent participants in a real cohort; nothing here estimates it, and no
NOVA-D interval substitutes for participant-cluster inference in a future biological protocol.

**Frozen code, unchanged.** `p6_n2.py` v0.1.0 sha256 `681ce852609e6ab4fcaf0f0ff603477c64e2351e10bd119811d45f008d4bbe8e`
and `p5_gate3a.py` `4600e2f3e3bb652affadcc6f6935380d8b9a77a4237e7cb306f056f849478130`, copied by hash, verified
before import, and called exactly as the NHANES wrapper called them (self-test §11). The near-tie rule is a policy
**outside** the selector (§8). No archive data and no biological data of any kind are read.

## 2. Hypothesis and the obstruction identity

**H_D:** A design-only near-tie rule can increase expected useful nonlinear acquisition gain without an offsetting
increase in estimation-induced selection regret.

With V_j = G_j + H_j, j_G = argmax_j G_j, j\* = argmax_j V_j and Γ = G_{j_G} − G_{j\*} ≥ 0,

    V_{j*} − V_{j_G} = (H_{j*} − H_{j_G}) − Γ .

A beneficial departure requires the nonlinear-excess differential to exceed the Gaussian ranking obstruction. The
near-tie rule selects decisions where Γ is small; it does not create nonlinear structure and may amplify noise.

**Aim limitation.** The rule observes only Γ̂_obs = (Ĝ_(1) − Ĝ_(2))/Ĝ_(1), the gap between the two highest
*estimated* Gaussian scores; these need not be the candidates involved in the optimal nonlinear departure. M4 (§5)
measures the fraction of beneficial opportunities the screen excludes, separating two failures that must not be
conflated: the estimator cannot recover a genuine advantage, versus the screen excludes the decisions where the
advantage exists.

## 3. Endpoints

    A = (L_∅ − L_G1)/L_∅ ,    Δ = (L_G1 − L_policy)/L_G1 ,    ρ = Δ/A   (A > 0) ,    κ = ρ/ρ_pop   (ρ_pop > 0)

All are **pooled-ratio** estimates: loss sums over evaluation rows, folds and replicates are formed first, then the
ratios — never an unweighted average of replicate-specific ratios. Δ and A are always reported separately.

**Exact identity (Pythagoras, R_S = R_∅ − V_S):** for the family of §4, ρ_pop = (V_{j\*}/V_{j_G} − 1)/(1 − A)
exactly. Hence E[ρ_policy] ≤ ρ_pop for every policy, with equality only if ĵ = j\* almost surely. This is why the
boundary condition cannot be scored against 0.10 (§13, v0.2 correction) and why κ is the recovery measure:

| κ | meaning |
|---|---|
| 1 | full recovery of the population advantage |
| (0, 1) | partial recovery |
| 0 | no improvement over G1 |
| < 0 | harmful acquisition decisions |

κ is **undefined and not computed** for M1 and M3 (ρ_pop = 0). The assay-utility condition A ≥ 0.05 holds **by
construction** in every condition because A is a controlled grid dimension (§6); it is not a post-hoc filter.

**Secondary, all reported:** absolute Δ; A; p_dep = P(ĵ_NOVA-S ≠ ĵ_G1); p_ben|dep; the opportunity rate
P(j\* ≠ ĵ_G1); the oracle decomposition (§7); Spearman(t, realised gain) as a diagnostic only, never compared with
the archive's 0.16.

## 4. The data-generating family, fully specified

Architecture preserved from P4D: O ∈ R²⁰, T ∈ R²⁰, P = 20 candidates, generated as distinct objects (no block
allocation needed; the P4D salt and allocation conventions are retained as reference only).

- **Base panel.** O ~ N(0, Σ_O), (Σ_O)_{ab} = φ^{|a−b|}, φ = 0.5.
- **Candidates.** X_j = a_jᵀO + u_j, u_j ~ N(0, s²) independent of O and across j, s² = 0.5, a_j scaled so
  Var(a_jᵀO) = 1 − s². Jointly Gaussian; r_j = X_j − E[X_j|O] = u_j.
- **Targets.** T_k = b_kᵀO + Σ_j λ_{kj} u_j + Σ_j ν_{kj} h(u_j/s) + ε_k, ε_k ~ N(0, σ²_ε) iid, h the normalised
  second Hermite polynomial h(z) = (z² − 1)/√2: E h = 0, Var h = 1, **E[z h(z)] = 0**. Then exactly
  G_j = s²Σ_k λ²_{kj}, H_j = Σ_k ν²_{kj}, V_j = G_j + H_j.
- **Structure, fixed per condition** from its structure seed and held across replicates (owner decision):
  a_j and b_k as scaled random directions; λ_{·j} = √(G_j/s²)·w_j with w_j a random unit vector in R²⁰ (dense);
  ν nonzero on exactly one designated candidate, ν_{·j_nl} = √H·w_ν with w_ν a random unit vector.
- **Equal marginal target variances.** Var(T_k) = c for every k, with c = 2·max_k Var(T_k|O) and
  Var(b_kᵀO) = c − Var(T_k|O). Reason: the frozen pipeline standardises each target by its training SD; with
  equal marginal variances that standardisation is, in population, a common scalar that cancels in A, Δ and ρ,
  so the oracle quantities computed in raw units are the quantities the pipeline faces. Verified in §11.
- **Log rule.** Every feature has negative values, so the frozen log rule takes its identity branch on every
  feature. The implementation asserts this on every fold and refuses otherwise.
- **Recorded limitation.** H_j ≥ 0 for every j by construction. The general signed decomposition permits H_j < 0
  (the claim H_j ≥ 0 was withdrawn on the X0 linear-proxy counterexample), so this is a favourable subfamily.
  Owner's qualification, recorded verbatim in substance: failure in this favourable family is strong evidence
  against the tested near-tie strategy **within this family and grid**, not a universal refutation of nonlinear
  acquisition; confirmation here establishes nothing about performance when the signed correction is negative,
  and the manuscript must not extrapolate to the entire class of nonlinear acquisition problems.

## 5. Mechanisms

G_j follows a prespecified profile G_1 = 1, G_2 = 1 − Γ_rel, G_j = (1 − Γ_rel)·r^{j−2} (j ≥ 3), with r solved
per condition (§6). Candidate index equals Gaussian rank.

| Class | Construction | Purpose | Expected under H_D |
|---|---|---|---|
| **M1** Gaussian null | H ≡ 0; Γ_rel = 0.20 | G1 is population-optimal; every departure is noise | ρ̄ ≤ 0 |
| **M2** near-tie, advantage at the runner-up | Γ_rel = 0.02; H on candidate 2 | the regime the screen is designed for | ρ̄ > 0 |
| **M3** near-tie, no advantage | Γ_rel = 0.02; H ≡ 0 | instability control | ρ̄ ≤ 0 |
| **M4** advantage away from the G frontier | Γ_rel = 0.20; H on candidate **7** (Gaussian rank 7 of 20) | aim control | high screen-exclusion of opportunities |

M1 and M3 have H ≡ 0 by definition and exist only at ρ_pop = 0; they are not crossed with the nonlinear
dimension. M4 places the nonlinear optimum outside the top two Gaussian candidates by construction; the
implementation asserts that the designated candidate is the argmax of V and sits at rank 7 in G.

## 6. Grid — frozen in full; every condition reported

**Calibration by construction.** For each condition, r, H and σ²_ε are solved in closed form so that the
population oracle attains the targets exactly, under a prespecified irreducible-noise share ψ = 20σ²_ε/R_∅ = 0.25:
H = Γ_rel + ρ_pop(1 − A) (M2; the M4 analogue uses G_7), r from Σ_jG_j = (1 − ψ)/A − H, σ²_ε from the residual.
The solver (`nova_d_calib.py`) verifies A and ρ_pop to 10⁻⁹ and checks the M4 rank placement; it uses only
population quantities of a known distribution — no selector, no losses, no outcome.

| Dimension | Frozen values |
|---|---|
| Mechanism × advantage | M1 (ρ_pop = 0); M3 (ρ_pop = 0); M2 × ρ_pop ∈ {0.05, 0.10, 0.40}; M4 × ρ_pop ∈ {0.05, 0.10, 0.40} — 8 combinations |
| Training sample size n | **128, 512, 1024** (§13, item 1) |
| Assay utility A | **0.05, 0.15, 0.20** (0.40 infeasible: common ceiling 0.242 at ψ = 0.25, set by M4) |
| Conditions | 8 × 3 × 3 = **72**, all feasible and verified |
| Monte Carlo replicates | 100 per condition (7,200) |
| Evaluation sample | 5,000 newly generated observations per replicate |
| Near-tie threshold γ | 0.02, 0.05, 0.10, 0.25, ∞ — a policy evaluated within every condition (§8), not a cell dimension |

ρ_pop = 0.05 lies below the relevance threshold and 0.10 on it, so the grid straddles the criterion.

**Seeds.** Master salt `NOVA-D-SYNTH-V03|`; child seed = SHA-256(salt | mechanism | ρ_pop | n | A | replicate |
role), role ∈ {structure, design, evaluation, folds}. No dependence on any archive salt.

## 7. Per-replicate procedure

1. Draw n training observations and 5,000 evaluation observations from the condition's DGP (distinct seeds).
2. Standardise targets by the **training** SD with the frozen 1 % median floor; apply the same divisor to the
   evaluation rows. No evaluation information enters preprocessing, selection or tuning.
3. Run the frozen procedure on the training sample exactly as `p6_n2.run_block('rbf')` does: K_OUT = 5 outer
   folds; in each, G1's pick, the cross-fitted scores, NOVA-V/U/S picks, the paired statistic t and the claim at
   τ = 2, all from that fold's training rows. n is the size handed to the frozen pipeline, which cross-fits
   internally.
4. In the same fold, compute Γ̂_obs from the G1 score vector on the same training rows (§8).
5. **External evaluation:** for each distinct pick and for the base panel, predict the 5,000 evaluation rows with
   the fold's training rows at the hyperparameters the frozen inner CV chose; record the pooled squared loss.
   Every arm and every γ-policy is scored on the same evaluation rows — paired by construction.
6. **Oracle quantities** (the DGP is known): G_j, H_j, V_j, j\*, j_G, Γ, the population approximation gap
   (V_{j\*} − V_{j_G})/R_∅, and the per-pick selection regret (V_{j\*} − V_ĵ)/R_∅ for G1, NOVA-S and NOVA-V.

Within-sample held-out losses are also computed (as the frozen pipeline reports them) and are used only for the
equivalence self-test, never as an endpoint.

## 8. Screening as a policy over a common decision population

For each γ: **use NOVA-S's pick where Γ̂_obs ≤ γ, G1's pick otherwise.** Both arms are evaluated on the same
decision population in every condition; the screened-versus-unscreened comparison is paired and exact. γ = ∞ is
unscreened NOVA-S. If Ĝ_(1) = 0, Γ̂_obs is undefined and the policy uses G1. The policy sits outside the selector.

## 9. Analysis

Per condition and arm (NOVA-S, NOVA-V, NOVA-U, and each γ-policy): pooled Δ, A, ρ; κ where ρ_pop > 0; a
**paired bootstrap over replicates** (B = 2,000, seeded; each replicate's five folds resampled as a unit; numerator
and denominator resampled together) giving the one-sided 95 % lower bound and the two-sided 95 % interval of ρ,
of Δ, and of κ; the paired interval of ρ_screened − ρ_unscreened at each γ; p_dep, p_ben|dep, the opportunity
rate, the screen-exclusion fraction among opportunities, and the mean oracle regret per arm.

## 10. Claims and decision rules, fixed before any run

| Claim | Condition | Required evidence |
|---|---|---|
| **Recoverability (primary)** | M2, ρ_pop = 0.10, n = 512, A = 0.15 | LCB₉₅(ρ_NOVA-S) > 0; κ reported with interval |
| **Regime viability at the relevance threshold** | M2, ρ_pop = 0.40, n = 1024, A = 0.15 (**positive control**) | LCB₉₅(ρ_NOVA-S) > 0.10, equivalently κ > 0.25 |
| **Near-tie rule effectiveness** | any condition, any γ < ∞ | paired interval of ρ_screened − ρ_unscreened excludes 0 above |

A positive primary result establishes recoverability, **not** a practically relevant advantage of 0.10. A regime
in which ρ ≥ 0.10 without screening helping is evidence about the mechanism, not about the rule.

**Positive-control failure** is reported and diagnosed with the oracle quantities of §7 — implementation defect,
estimator calibration, finite-sample limitation and nuisance-estimation variance are all candidate explanations
— and is never attributed automatically to a coding defect. The remaining grid is reported regardless.

**Validity controls.** Any M1 or M3 condition with ρ̄ ≥ 0.10 means the readout is measuring instability rather
than information; reported as a diagnostic failure that voids the M2 results of the same run.

**Aim control.** M4's screen-exclusion fraction is reported as a limitation of the rule, not of nonlinear
acquisition.

**Descriptive.** All other conditions. The grid maximum of ρ may be reported as the most favourable predefined
condition, labelled **upward-biased**; any inferential claim on it must account for selecting the maximum. For
the full grid, intervals are descriptive unless a multiplicity-adjusted procedure has been prespecified (none is).

**Falsification.** If the primary and the positive control both fail, and no condition clears recoverability, the
frozen estimator does not recover a known nonlinear advantage in this family at these sample sizes; recorded and
the grid is not extended.

## 11. Integrity controls and self-test

The failure mode of synthetic work is DGP selection, not data reuse. Accordingly: everything in §4–§10 is hashed
before the first run; no condition is added, removed or re-parameterised after any result is seen (an extension is
a new version declaring what had been observed); every condition is reported; **significance is never a
criterion** (unlimited replicates make any nonzero effect significant) — relevance is the effect size against
0.10 and recoverability is κ; exploratory pre-selection statistics come only after §10 is decided, under a
separate name; `run` refuses unless frozen and the hashes match; one run per freeze.

**Self-test (`nova_d_run.py selftest`, must pass before freeze; produces no NOVA-D outcome):**
1. SHA-256 of both frozen modules verified before import.
2. The G1 score-vector mirror (verbatim body of the frozen `g1_pick`, used only for Γ̂_obs) has argmax equal to
   the frozen `g1_pick` **exactly** on 200 random blocks.
3. The fold loop reproduces the frozen `run_block('rbf')` on the frozen module's own synthetic blocks: identical
   picks, within-sample losses equal to 10⁻⁹ or better.
4. All 72 conditions feasible; oracle A and ρ_pop verified to 10⁻⁹; M4 rank placement verified.
5. **Sampler instantiates the calibrated distribution** (owner's safeguard): on six designated conditions,
   400,000 rows, independently estimated G_j, H_j, R_∅, A, ρ_pop, Hermite orthogonality E[z h(z)] and the
   equal-marginal-variance property are compared with the closed forms by batch-means z-tests (|z| < 4) and
   absolute caps (0.02 at the experiment's scale). Assert-only: no parameter is changed in response.
6. The log rule takes its identity branch on every feature (negatives present) — asserted here and on every fold.
7. Timing of the frozen solve at each n, for the cost projection.

**Result at drafting (5 Nov 2026):** PASS — 0/200 mirror mismatches; equivalence 1.8 × 10⁻¹⁶, 0 pick mismatches;
72/72 feasible; all six moment checks within tolerance (max |z| 3.4 over 20 orthogonality tests, A and ρ_pop
within 0.009 of target). Record: `nova_d_out/selftest.json`.

## 12. Computational feasibility (measured, owner's shell, 4 cores)

| n | inner solve | per replicate | 2,400 replicates |
|---|---|---|---|
| 128 | 41 ms | 26 s | 17 h |
| 512 | 36 ms | 33 s | 22 h |
| 1024 | 157 ms | 109 s | 73 h |
| 2048 | 762 ms | 492 s | 328 h |

The cost is cubic in n: the frozen procedure re-forms and eigendecomposes the kernel for each of 21 models in
each of 25 inner folds per replicate. With n = 2048 the grid is ~370 core-hours; with n = 1024 it is ~112
core-hours, about 30 h wall at four-way parallelism. The run is partitionable by condition (`--cond a:b`) and
resumable, because every (condition, replicate) is deterministic and independent. Cost is reduced **prospectively**
and never by stopping selectively.

## 13. Changes from v0.2 and the single open item

1. **n = 2048 → 1024** on measured cost (§12). Positive control moves to n = 1024. **Requires owner confirmation;
   the freeze waits on it.** No real-data route examined by this programme has supplied 2,048 one-per-participant
   rows (NHANES: 1,725 eligible), so {128, 512, 1024} also brackets the realistic range.
2. A = 0.40 → 0.20 (owner decision; infeasible, common ceiling 0.242).
3. Primary criterion LCB₉₅(ρ) > 0.10 → LCB₉₅(ρ) > 0 with κ as the recovery measure (owner decision); the 0.10
   criterion retained at the positive control (κ > 0.25).
4. Pooled-ratio estimation; paired bootstrap of numerator and denominator; Δ and A reported separately.
5. Sampler-instantiation self-test (owner's safeguard).
6. Equal-marginal-variance construction and the b_k scaling (§4), added because the frozen standardisation would
   otherwise reweight targets and detach the pipeline's losses from the oracle quantities.
7. Two inference levels stated (§1). Positive-control failure handling rewritten (§10). Owner's qualification on
   the H ≥ 0 limitation recorded (§4).

Freeze sequence on confirmation of item 1: `nova_d_run.py freeze` (records the protocol, implementation and
calibration hashes; refuses if any output exists or the self-test has not passed) → commit the freeze record →
`run`.
