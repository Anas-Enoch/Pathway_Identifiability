# NOVA-D — Detectability of Nonlinear Acquisition Advantage: synthetic protocol **draft** v0.2 (5 Nov 2026)

**Status: DRAFT for owner review. Not frozen, not hashed. No simulation may be run under a draft.**
Author of record: Anas Enoch. Drafted by Claude (Opus 5) under the owner's decisions of 5 Nov 2026.
Changes from v0.1 are listed in §13; three of them alter the owner's stated grid and are flagged for confirmation.

## 1. Scope

NOVA-D is a synthetic-data programme answering three questions:

1. **Existence** — under what joint distributions does a meaningful nonlinear acquisition advantage exist?
2. **Recoverability** — can the frozen NOVA-S estimator recover it at realistic sample sizes?
3. **Selection** — does near-tie screening increase useful departures, or merely estimation-driven disagreement?

A negative answer to (2) or (3) is a substantive result.

**What this is not.** Not an external validation. Owner's recorded boundary: synthetic success would justify renewed
external testing; it would not validate NOVA-S on biological data. A successful synthetic regime yields a precise
biological prediction, not a validation claim. NOVA-D can show that a mechanism generates useful nonlinear
acquisition opportunities; it cannot show that the mechanism is common in real metabolomics cohorts.

**Frozen code, unchanged.** `p6_n2.py` v0.1.0 sha256 `681ce852609e6ab4fcaf0f0ff603477c64e2351e10bd119811d45f008d4bbe8e`
and `p5_gate3a.py` `4600e2f3e3bb652affadcc6f6935380d8b9a77a4237e7cb306f056f849478130`, copied by hash and verified at
run time. Algorithms byte-identical. Preprocessing conventions, acquisition budget B = 1 and the evaluation
procedure are preserved. The near-tie rule operates **outside** the selector (§6). No archive data is read. No
biological data of any kind is read.

## 2. Hypothesis and the obstruction identity

**H_D:** A design-only near-tie rule can increase expected useful nonlinear acquisition gain without an offsetting
increase in estimation-induced selection regret.

With V_j = G_j + H_j, j_G = argmax_j G_j, j\* = argmax_j V_j and Γ = G_{j_G} − G_{j\*} ≥ 0,

    V_{j*} − V_{j_G} = (H_{j*} − H_{j_G}) − Γ.

A beneficial departure requires the nonlinear-excess differential to exceed the Gaussian ranking obstruction. The
near-tie rule selects decisions where Γ is small; it does not create nonlinear structure and may amplify noise.

**Aim limitation.** The rule observes only Γ̂_obs = (Ĝ_(1) − Ĝ_(2))/Ĝ_(1), the gap between the two highest
*estimated* Gaussian scores. These need not be the candidates involved in the optimal nonlinear departure. M4 (§5)
measures this directly, separating two distinct failures: the estimator cannot recover a genuine advantage, versus
the screen systematically excludes the decisions where the advantage exists.

## 3. Primary endpoint and its interpretation

    A = (L_∅ − L_G1)/L_∅ ,    Δ = (L_G1 − L_NOVA)/L_G1 ,    ρ = Δ/A   (defined for A > 0)

tested against **0.10**, the factor confirmed in `EXTERNAL_RELEVANCE_FACTOR_CONFIRMATION_v1.0.md`. Within a cell, Δ
and A are formed as pooled ratios over evaluation rows and replicates, then divided. Both are reported separately;
a favourable ρ never substitutes for the confidence-bound requirement on Δ.

**Assay-utility condition.** A ≥ 0.05, a prospective design choice, not an empirically established threshold. No
upper bound on A. In this protocol the condition is satisfied **by construction** in every cell, because A is a
controlled grid dimension (§6) rather than a post-hoc filter.

**Interpretation.** For the family of §4,

    ρ ≈ (1/(1 − A)) · (V_{j*}/V_{j_G} − 1) ,

so ρ ≥ 0.10 requires the true value-maximiser to exceed the Gaussian pick by roughly ten per cent **in value**;
combined with §2, it requires (H_{j\*} − H_{j_G} − Γ)/V_{j_G} ≥ 0.10(1 − A). This identity is the reason the
nonlinear dimension of the grid is parameterised in ρ_pop rather than in an arbitrary signal multiplier.

**Secondary, all reported:** absolute Δ; A; p_dep = P(ĵ_NOVA ≠ ĵ_G1); p_ben = P(R(ĵ_NOVA) < R(ĵ_G1));
E[R(ĵ_G1) − R(ĵ_NOVA)]; the oracle decomposition of §7. Spearman(t, realised gain) is a diagnostic only, never an
endpoint, and is not compared with the archive's 0.16.

## 4. The data-generating family, fully specified

Architecture preserved from P4D: O ∈ R²⁰ (base panel), T ∈ R²⁰ (targets), P = 20 candidates. The three roles are
generated as distinct objects, so no block-allocation rule is required; the P4D salt and allocation conventions are
retained as documented reference rules only and are **not** reused as the synthetic randomisation mechanism.

**Base panel.** O ~ N(0, Σ_O), (Σ_O)_{ab} = φ^{|a−b|}, φ = 0.5.

**Candidates.** X_j = a_jᵀO + u_j, with u_j ~ N(0, s²) independent across j and of O, s² = 0.5, and a_j scaled so
Var(a_jᵀO) = 1 − s². The a_j are drawn once per cell from the cell's seed and then held fixed across replicates.
Hence (O, X) is jointly Gaussian and r_j = X_j − E[X_j|O] = u_j.

**Targets.** For k = 1,…,20,

    T_k = b_kᵀO + Σ_j λ_{kj} u_j + η Σ_j ν_{kj} h(u_j/s) + ε_k ,   ε_k ~ N(0, σ_ε²) iid,

with h the normalised second Hermite polynomial h(z) = (z² − 1)/√2, so E h = 0, Var h = 1 and **E[z h(z)] = 0**.
The orthogonality is what makes G and H separately controllable.

**Consequences (exact, since u_j ⫫ u_{j'} and u_j ⫫ O):**

    E[T|O,X_j] − E[T|O] = λ_{·j} u_j + η ν_{·j} h(u_j/s)
    G_j = Σ_k λ_{kj}² s²        V_j = Σ_k (λ_{kj}² s² + η² ν_{kj}²)        H_j = η² Σ_k ν_{kj}² .

**Recorded limitation.** This family yields H_j ≥ 0 for every j by construction. The general signed decomposition
permits H_j < 0 — the claim H_j ≥ 0 was withdrawn on a linear-proxy counterexample (X0) — so this is a *favourable*
subfamily. Consequence for inference: a refutation of H_D here is strong, because the family is favourable to it;
a confirmation here is weak, and §1 governs what may be claimed from it.

## 5. Mechanisms

G_j is set by a prespecified profile with relative top gap Γ_rel = (G_(1) − G_(2))/G_(1); ν is nonzero on exactly
one designated candidate; η and σ_ε are solved for (§6).

| Class | Construction | Purpose | Expected under H_D |
|---|---|---|---|
| **M1** Gaussian null | ν ≡ 0, H ≡ 0; Γ_rel = 0.20 | G1 is population-optimal; every departure is noise | ρ̄ ≤ 0 |
| **M2** near-tie, advantage at the runner-up | Γ_rel = 0.02; ν nonzero on the G-runner-up | the regime the screen is designed for | ρ̄ > 0, rising with ρ_pop |
| **M3** near-tie, no advantage | Γ_rel = 0.02, ν ≡ 0, H ≡ 0 | instability control: is the readout picking up estimation instability? | ρ̄ ≤ 0 |
| **M4** advantage away from the G frontier | Γ_rel = 0.20; ν nonzero on the candidate at **G-rank r = 7 of 20** | aim control: does the screen exclude the decisions where the advantage lies? | high screen-exclusion of beneficial opportunities |

M4 places the nonlinear-optimal candidate outside the top two Gaussian candidates by construction, as the owner
required; without that the obstruction cannot be tested. M1 and M3 have H ≡ 0 by definition and therefore exist
only at ρ_pop = 0 — they are **not** crossed with the nonlinear dimension (§13, correction 1).

## 6. Grid — frozen in full; every cell reported

Both A and ρ_pop are **targeted by construction**: given the G profile, Γ_rel and the designated ν candidate, the
magnitude of ν and σ_ε are solved in closed form from the §4 identities so that the cell attains its target
(A, ρ_pop) at the population level. This calibration uses only population oracle quantities of a known
distribution — no selector, no losses, no outcome — and is itself hashed with the protocol.

| Dimension | Frozen values |
|---|---|
| Mechanism × advantage | M1 (ρ_pop = 0); M3 (ρ_pop = 0); M2 × ρ_pop ∈ {0.05, 0.10, 0.40}; M4 × ρ_pop ∈ {0.05, 0.10, 0.40} — **8 combinations** |
| Training sample size n | 128, 512, 2048 |
| Assay utility A (population oracle) | 0.05, 0.15, 0.40 |
| **Conditions** | 8 × 3 × 3 = **72** |
| Monte Carlo replicates | **100** per condition (7,200 experiments) |
| Evaluation population | **5,000** newly generated observations per replicate |
| Near-tie threshold γ | 0.02, 0.05, 0.10, 0.25, ∞ — evaluated **within** every condition as a policy (§8), not a cell dimension |

ρ_pop = 0.05 is below the relevance threshold and ρ_pop = 0.10 sits on it, so the grid straddles the criterion
rather than only exceeding it. Computational feasibility is to be checked before freezing; if the grid cannot be
run, cost is reduced **prospectively**, never by stopping selectively after favourable outcomes.

**Seeds.** New synthetic master salt `NOVA-D-SYNTH-V02`, with no dependence on any archive salt. Child seeds are
SHA-256 of `salt | mechanism | condition-index | replicate | role`, role ∈ {design, evaluation, structure}, so
every draw is independently reproducible and the synthetic experiment is not coupled to historical archive
construction.

## 7. Per-replicate procedure

1. Draw n training observations and 5,000 independent evaluation observations from the condition's DGP.
2. Run the frozen pipeline: G1 and NOVA-S picks from training rows only; losses on evaluation rows. Both arms
   receive identical training rows, identical blocks, identical preprocessing and an identical tuning budget.
3. Record L_∅, L_G1, L_NOVA, the picks, the departure indicator, t, and Γ̂_obs.
4. **Oracle quantities** (available because the DGP is known): G_j, V_j, H_j for all j; j\*, j_G; Γ; and
   - population approximation gap  R(j_G) − R(j\*)
   - finite-sample selection regret  R(ĵ_NOVA) − R(j\*)  and  R(ĵ_G1) − R(j\*).

## 8. Screening as a policy, not a filter

For each γ, the screened policy is: **use NOVA-S's pick where Γ̂_obs ≤ γ, G1's pick otherwise.** Both arms are
therefore evaluated over the *same* prespecified decision population in every condition, and the screened-versus-
unscreened comparison is exact and paired by construction. Filtering blocks out of the experiment would change the
estimand's population, so that an apparent rise in mean ρ could reflect easier blocks rather than better
acquisition; this formulation removes that artefact. γ = ∞ recovers unscreened NOVA-S. The policy sits outside the
selector; NOVA-S is unchanged.

## 9. Analysis

Per condition and γ: ρ̄ with a **paired bootstrap over replicates of (Δ, A) jointly** — percentile interval, the
denominator resampled with the numerator, never separately — plus Δ and A with their own intervals; p_dep, p_ben,
E[loss difference]; the oracle decomposition; and for M4 the screen's exclusion rate among beneficial
opportunities.

## 10. Claims and decision rules, fixed before any run

Two claims are reported separately. Neither is required for the other.

| Claim | Required evidence |
|---|---|
| **Regime viability** | frozen NOVA-S attains ρ̄ ≥ 0.10 with its paired bootstrap 95 % lower bound above 0.10 in a prespecified condition |
| **Near-tie rule effectiveness** | ρ̄_screened > ρ̄_unscreened at some prespecified γ, paired bootstrap interval of the difference excluding 0, over the common decision population of §8 |

A regime where both selectors reach ρ ≥ 0.10 without screening helping is evidence about the data-generating
mechanism; it is **not** evidence that the near-tie rule adds value.

**Three tiers, designated now:**
- **Positive control:** M2, ρ_pop = 0.40, n = 2048, A = 0.15. If regime viability fails here, the implementation
  or the estimator — not the mechanism — is the obstruction, and no other cell may be interpreted until that is
  diagnosed.
- **Primary scientific endpoint:** M2, ρ_pop = 0.10, n = 512, A = 0.15 — threshold-level advantage at the one
  empirically anchored utility value. This single prespecified condition carries the primary claim, so the primary
  estimate is not a maximum and needs no selection adjustment.
- **Descriptive:** all remaining conditions. The grid maximum may be reported as the most favourable predefined
  condition, labelled **upward-biased**; any significance claim on it must account for having selected the maximum.

**Validity controls.** If any M1 or M3 condition attains ρ̄ ≥ 0.10, the readout is measuring instability rather
than information. That is a diagnostic failure of the design, reported as such, and it voids the M2 results of the
same run rather than supporting them.

**Aim control.** M4 reports the fraction of beneficial nonlinear opportunities the screen excludes. A high fraction
is a limitation of the rule, not of nonlinear acquisition, and must be reported as the former.

**Falsification.** If no M2 condition clears regime viability across the entire frozen grid, H_D is refuted on its
own prespecified grid — and, per §4, in a subfamily favourable to it. The programme records that the near-tie route
does not reach the confirmed relevance threshold even under known nonlinearity, and stops. The grid is not
extended.

## 11. Integrity controls

The failure mode of synthetic work is DGP selection, not data reuse.

1. Mechanisms, DGP family, calibration code, grid, seeds, readout and decision rules are hashed before the first run.
2. No condition may be added, removed or re-parameterised after any result is seen. An extension is a new protocol
   version carrying an explicit declaration of what had already been observed.
3. Every condition is reported, including failures. None omitted, pooled away or post-selected.
4. **Significance is not a criterion.** With unlimited replicates any nonzero effect can be made significant, so
   conclusions are effect sizes against 0.10.
5. Exploratory pre-selection statistics beyond t may be examined only after §10 has been decided and recorded,
   under a separate name, never by re-running the frozen grid with a substituted statistic.
6. No simulation runs until this protocol and the executable implementation are both frozen and hashed.

## 12. Sequencing (owner's decision, recorded)

1. The 0.10 factor and the denominator interpretation, including the 5 % assay-utility condition — **done**
   (`EXTERNAL_RELEVANCE_FACTOR_CONFIRMATION_v1.0.md`, confirmed before any synthetic outcome existed).
2. Freeze NOVA-D v0.2: mechanisms, grid, paired comparisons, seeds, reporting rules.
3. Implement and hash.
4. Run the full grid with no outcome-dependent modification.
5. Interpret; formulate future questions.

Any eligibility rule for a future biological protocol that is revised in response to these results is **newly
developed**, not an inherited confirmatory criterion, and must be labelled so.

## 13. Changes from v0.1, and items needing owner confirmation

1. **M1 and M3 are no longer crossed with the nonlinear dimension.** They are defined by H ≡ 0, so cells such as
   "M1 at multiplier 1.5" are self-contradictory. Mechanism and advantage are nested: 8 combinations.
2. **A is a controlled grid dimension** at {0.05, 0.15, 0.40}, replacing the noise-SD dimension; σ_ε is solved for.
   This implements the owner's 5 % floor by construction in every cell instead of as a post-hoc filter, and
   prevents A from varying uncontrolled under a ratio endpoint. 0.15 is the one empirically anchored value
   (NHANES 15.8 %).
3. **The nonlinear dimension is parameterised as ρ_pop ∈ {0.05, 0.10, 0.40}**, the criterion's own currency, rather
   than as a signal multiplier {0, 0.5, 1.5}. Under §3's identity the owner's multipliers correspond to population
   advantages far above 0.10, so they would have probed recoverability only; ρ_pop straddles the threshold.
4. Screening is a **policy over a common decision population** (§8) rather than a block filter, which implements
   the owner's safeguard exactly and removes γ as a cell dimension.
5. Three-tier claim structure (§10) with a prespecified primary condition, so the primary estimate is not a grid
   maximum.
6. M4's nonlinear-optimal candidate is fixed at G-rank 7 of 20 (owner required "outside the top two").

Items 1–3 alter the owner's stated grid and require confirmation. The cell count (72) and the replicate and
evaluation budgets (100, 5,000) are unchanged from the owner's specification.

**Open for owner:** (i) confirm 1–3; (ii) Γ_rel values 0.02 / 0.20 and M4's rank 7; (iii) the primary condition
designated in §10; (iv) φ = 0.5, s² = 0.5 and the fixed-structure-per-cell convention in §4.
