# NOVA-D — Detectability of Nonlinear Acquisition Advantage: synthetic protocol **draft** v0.1 (5 Nov 2026)

**Status: DRAFT for owner review. Not frozen, not hashed, not executable yet.** No run may occur under a draft.
Author of record: Anas Enoch. Drafted by Claude (Opus 5). Open items for owner decision are listed in §11.

## 1. Scope, and what this is not

NOVA-D is a synthetic-data programme. It tests whether a **design-only near-tie rule** can raise the expected useful
nonlinear acquisition gain without an offsetting increase in estimation-induced selection regret.

It is not an external validation, and it is not evidence about biological data. Owner's recorded boundary, carried
here verbatim in substance: *synthetic success would justify renewed external testing; it would not, by itself,
validate NOVA-S on biological data.*

NOVA-S is unchanged. `p6_n2.py` v0.1.0 sha256 `681ce852609e6ab4fcaf0f0ff603477c64e2351e10bd119811d45f008d4bbe8e`
and `p5_gate3a.py` `4600e2f3e3bb652affadcc6f6935380d8b9a77a4237e7cb306f056f849478130` are copied by hash and
verified at run time, as in the NHANES wrapper. **The near-tie rule governs which blocks enter the experiment; it
never operates inside the selector.** Any pre-selection statistic explored under §9.4 is a separately named
experimental method, not a modification of the frozen NOVA-S validation claim.

The 186-record archive is not used. No biological data of any kind is read.

## 2. Hypothesis

**H_D:** A design-only near-tie rule can increase expected useful nonlinear acquisition gain without an offsetting
increase in estimation-induced selection regret.

The mechanism is an obstruction, not an opportunity. With V_j = G1_j + H_j and Γ = G1_{j_G1} − G1_{j*} ≥ 0 (where
j_G1 = argmax G1, j* = argmax V),

    V_{j*} − V_{j_G1} = (H_{j*} − H_{j_G1}) − Γ,

so a beneficial departure requires **H_{j*} − H_{j_G1} > Γ**. A large G-gap is a floor the nonlinear excess must
clear. Near-tie design selects decisions where that floor is low. It does not create nonlinear structure, and it
may amplify estimation noise. The sign of the net effect is the object of this protocol and is not assumed.

**Known aim limitation.** The rule can only screen on the *observable* gap between the top two estimated G1 scores,
Γ̂_obs = (Ĝ_(1) − Ĝ_(2))/Ĝ_(1). This equals the relevant Γ only when j* is the G1 runner-up. Where the nonlinear
advantage sits on a candidate far down the G1 ranking, the screen is aimed at the wrong gap. Mechanism M4 (§5) tests
this directly, so that a negative result in M2 is not mistaken for "nonlinear acquisition never helps" when it means
"the screen looks in the wrong place."

## 3. Primary readout

The primary quantity is the external criterion's own currency, so that a synthetic result is commensurable with the
confirmed external relevance threshold and across mechanisms of different scale:

    ρ = Δ / Δ_G1→∅ ,    Δ = 1 − L_NOVA-S / L_G1 ,    Δ_G1→∅ = 1 − L_G1 / L_O

tested against **0.10** (`EXTERNAL_RELEVANCE_FACTOR_CONFIRMATION_v1.0.md`). Within a cell, Δ and Δ_G1→∅ are formed
as pooled ratios over evaluation rows and replicates first, then divided — matching how the stub computes them.
Cells with Δ_G1→∅ ≤ 0 are reported separately and excluded from ρ: a candidate panel where measurement is worthless
is not an acquisition problem.

Absolute Δ, the departure rate p_dep, the beneficial-departure rate p_ben = P(R(ĵ_NOVA) < R(ĵ_G1)), and
E[R(ĵ_G1) − R(ĵ_NOVA)] are all reported as secondary. Spearman(t, realised gain) is reported as a diagnostic only;
it is not an endpoint and is not compared to the archive's 0.16.

## 4. Units and procedure per replicate

For cell c and replicate r (seed derived deterministically from a master salt fixed at freeze):

1. Draw a design sample of n_design participants and an independent evaluation sample of n_eval from the same DGP.
   n_eval is set large (§6) so that evaluation noise is negligible against selection noise: the object of study is
   selection regret, not evaluation error.
2. Construct blocks by the frozen P4D rule (20 O / 20 T / 20 P).
3. Apply the near-tie screen at threshold γ using **design rows only**: block admitted iff Γ̂_obs ≤ γ. γ = ∞ is the
   unscreened comparator and is run in every cell.
4. Run the frozen pipeline on admitted blocks: G1 and NOVA-S picks from design rows, losses on evaluation rows.
   Both arms receive identical design rows, blocks, preprocessing and tuning budget.
5. Record L_O, L_G1, L_NOVA-S, the picks, the departure indicator, and t.
6. **Oracle quantities**, available because the DGP is known: V_j and G1_j for every candidate, j*, j_G1^*, and the
   decomposition
   - population approximation gap  R(j_G1^*) − R(j*)
   - finite-sample selection regret  R(ĵ_NOVA) − R(j*)   and   R(ĵ_G1) − R(j*).

## 5. Mechanism classes

Base construction, chosen so G1 and H are separately controllable. O ~ N(0, Σ_O); candidates X_j = a_jᵀO + u_j with
u_j ~ N(0, s_j²) independent, so (O, X) is jointly Gaussian; target

    T = m(O) + Σ_j λ_j u_j + η·h(u_{j_nl}) + ε,   ε ~ N(0, σ²),

with h normalised to mean 0 and unit variance and **orthogonal to the identity**, E[u h(u)] = 0 (second Hermite
polynomial in u/s). Then for a purely linear candidate G1_j = λ_j² s_j² and H_j = 0; for a purely nonlinear
candidate G1_j = 0 and H_j = η². A mixed candidate carries both.

| Class | Construction | Purpose | Expected under H_D |
|---|---|---|---|
| **M1** Gaussian null | all λ free, η = 0 ⇒ H_j ≡ 0 | G1 is population-optimal; every departure is noise | ρ̄ ≤ 0 |
| **M2** Near-tie with advantage at the runner-up | j_a, j_b with G1_a ≈ G1_b, H_b = η² > 0, H_a = 0 ⇒ j* = j_b, Γ small | the regime the screen is designed for | ρ̄ > 0, rising with γ tightening |
| **M3** Near-tie without advantage | G1_a ≈ G1_b, H_a = H_b = 0 | instability control: does the procedure mistake estimation instability for information? | ρ̄ ≤ 0 |
| **M4** Advantage away from the G1 frontier | H large on a candidate with low G1 (pure nonlinear) | aim control: does the screen *exclude* the blocks where departure would help? | screen exclusion rate high ⇒ rule mis-aimed |

M1–M3 are the owner's three classes. M4 is added in this draft and is an open item (§11).

## 6. Parameter grid (frozen in full at freeze; every cell reported)

| Parameter | Values |
|---|---|
| n_design | 100, 250, 500, 1000, 2000 |
| d_O (base panel) | 10, 20 |
| p (candidates per block) | 20 |
| ζ = η² / max_j G1_j (nonlinear strength) | 0.8, 1.0, 1.25, 2.0 — ζ ≤ 1 means **no** population advantage exists |
| Γ_rel in the generating design | 0.01, 0.05, 0.20 |
| γ (near-tie screen) | 0.02, 0.05, 0.10, 0.25, ∞ |
| σ² | set so Δ_G1→∅ spans a plausible range; the achieved Δ_G1→∅ is reported per cell |
| Monte Carlo replicates R | 500 per cell |
| n_eval | 10,000 |

## 7. Analysis

Per cell: ρ̄ with a percentile bootstrap interval over replicates; p_dep, p_ben, E[loss difference]; the oracle
decomposition; the screen's admission and exclusion rates.

## 8. Decision rules, fixed before any run

**Primary (H_D supported) requires both:**
1. Some M2 cell with γ < ∞ attains ρ̄ ≥ 0.10 with its Monte Carlo 95 % lower bound above 0.10; **and**
2. in that same cell, the screened ρ̄ exceeds the unscreened (γ = ∞) ρ̄ with the paired bootstrap interval of the
   difference excluding 0. Condition 2 is required because the hypothesis concerns the *rule*, not the regime.

**Validity controls.** If any M1 or M3 cell attains ρ̄ ≥ 0.10, the readout is measuring instability rather than
information. That is a diagnostic failure of the design, reported as such, and it voids any M2 result in the same
run rather than supporting it.

**Aim control.** M4 reports the fraction of blocks with a genuine beneficial departure that the screen excludes. A
high fraction is recorded as a limitation of the rule, not of nonlinear acquisition.

**Falsification.** If no M2 cell clears both primary conditions over the entire frozen grid, H_D is refuted on its
own prespecified grid. The programme records that the near-tie route does not reach the external relevance
threshold even under a favourable DGP with known nonlinearity, and stops. The grid is not extended.

## 9. Integrity controls specific to synthetic work

The failure mode here is not data reuse — it is DGP selection. Accordingly:

1. Mechanisms, grid, seeds, readout and decision rules are hashed before the first run.
2. No cell may be added, removed or re-parameterised after any result is seen. A grid extension is a new protocol
   version carrying an explicit declaration of what had already been observed.
3. Every cell is reported, including losses. No cell is omitted, pooled away, or post-selected.
4. **Significance is not a criterion.** With unlimited replicates any nonzero effect can be made significant, so
   the conclusion is the effect size against 0.10, never a p-value.
5. §9.4 — exploratory pre-selection statistics beyond t may be examined, but only after §8 has been decided and
   recorded, under a separate name, and never by re-running the frozen grid with a new statistic substituted.

## 10. Relation to the suspended data search

If H_D is supported, the supporting M2 cells specify the DGP signature a qualifying external cohort must have —
near-tied top G1 scores, nonlinear excess exceeding the G-gap, sample size at or above the supporting cell's
n_design. That specification, not a repository sweep, is what a resumed Phase II would screen for. If H_D is
refuted, no resumed search is warranted on this hypothesis.

## 11. Open items for owner decision

1. Adopt M4 (aim control) or run the owner's three classes only.
2. The grid values in §6, R = 500, and n_eval = 10,000.
3. Whether condition 2 of §8 (screened must beat unscreened) is required, as drafted, or reported only.
4. §4 step 2: whether the P4D block rule and its master salt apply unchanged to synthetic candidates, or a
   synthetic block rule is specified instead.
5. Whether the Δ_G1→∅ eligibility window left open in
   `EXTERNAL_RELEVANCE_FACTOR_CONFIRMATION_v1.0.md` §4 should be fixed before or after NOVA-D reports, given that
   §10 makes NOVA-D's output an input to that choice.
