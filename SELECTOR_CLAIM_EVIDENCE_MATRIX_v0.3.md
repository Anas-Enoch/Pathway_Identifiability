# Acquisition-selector programme: claim–evidence matrix v0.3

Date: 2026-11-02; v0.2 of 2026-11-03 adds the Gate 3A per-family results (§1.5) and the G1 self-agreement finding (§3.3); v0.3 (3 Nov, later) adds the N2 failure-regime diagnostic (§1.6, use 11) and closes the archive for selector development (§3). Author of record: Anas Enoch; compiled by Claude (Fable 5.1) from the frozen result records named
below, opened in this session. Nothing here is a new computation; every number is copied from a record and labelled
with the record it came from. Evidence grades: **CERTAIN** = read from the named frozen record (hash recorded there);
**SUPPORTED** = stated in a record but the underlying output file was not opened here; **UNKNOWN** = the record was
not opened in this session.

Records opened: `GATE2_DECISION_RECORD.md`, `GATE3A_PROTOCOL_v1.0.md.md`, `N1_RESULT_RECORD.md`, `N2_RESULT_RECORD.md`,
`N2_REPLAY_RECORD.md`, `ESTIMAND_AUDIT_RECORD.md`, `X01_RESULT_RECORD.md`, `L0_HIERARCHY_RECORD.md`,
`HELDOUT_VALUE_OUTCOME_v1.md`, `EXTERNAL_VALIDATION_STATUS.md`, `GDOA0_SPEC.md`, `UAPA_F1_SPEC.md`,
`NOVA_PAPER_v0.4.md` (§1, §4, §6), project records `uapa_f1_result.md` (F1, B23, NBRA pre-check, GDOA-0 result,
theory-package reviews), `PAPER4_CORRECTIONS_P4U_STAGE10K.md` (H19–H20), `P4V_B3_EVIDENCE_AUDIT_RESULT_v1.0.md`,
`TRACKB_P4V_STATUS.md`. Not opened: Paper III/IV manuscripts and their result tables,
the DEP/DIAF/DEP-U/MS/GM theory freezes, the 9 Oct predictor-benchmark outputs.

## 0. The one dataset

Every selector comparison in the programme after Paper I was run on the **186-record Metabolomics Workbench
archive** (188 canonical matrices, 952 frozen blocks of 20 base / 20 target / 20 candidate metabolites; design rows
plus an archived evaluation half). Recorded uses of that archive, in order [CERTAIN from the records that list them]:

| # | Use | Rows touched | Date |
|---|---|---|---|
| 1 | Paper III held-out evaluation (gaussian = G1, variance, D-optimal, E-optimal; B = 1, 2, 5, 10) | design + evaluation | 29 Sep 02:36Z |
| 2 | Paper IV held-out primary (G1 vs RIG, B = 1) — both arms already present in #1 (185/186 studies exactly equal) | evaluation | 29 Sep 20:05Z (prereg 18:33Z) |
| 3 | 9 Oct retrospective predictor benchmark (quadratic learner selected in 24/746 blocks on the evaluation half) | design + evaluation | 9 Oct |
| 4 | UAPA F1 / B23 (joint A-optimal, nested CV selection; B = 1–3) | design | 9 Oct |
| 5 | GDOA-0 signal-ceiling gate | design | 9 Oct |
| 6 | Gate 2 (statistical feasibility; hindsight headroom, MDE) | evaluation cells of #1–#2 | 19 Oct |
| 7 | Gate 3A (predictor families P0–P3; acquisition regret vs G1) | design (+ evaluation half as secondary) | 19 Oct |
| 8 | N0 (regularization parity) | design | 20 Oct |
| 9 | NOVA N2 (NOVA-V/U/S, two nuisances) | design | 20 Oct |
| 10 | N2 replay (instrumented twin; threshold envelope) | design | 20 Oct |
| 11 | N2 failure-regime diagnostic (five pre-acquisition covariates vs realised gain given t; `N2_REGIME_DIAGNOSTIC_RESULT_v1.0.md`) | frozen outputs of #7, #8, #10 and P4D manifests only; no rows | 3 Nov |

The evaluation half is therefore not held out for anything decided after 29 Sep; the design rows have been used for
selector choice at least seven times. The reserved **Gate D pool** (512 clusters, OSF `osf.io/egnqy`, registered
19 Oct 17:27Z) remains sealed [CERTAIN, Gate 3A §9]. The internal 93-study paired assessment was never authorised
[CERTAIN, Track B status].

## 1. Matrix

Columns: definition; implementation and tuning; evaluation history; result; strongest admissible claim.
"Δ_acq" = equal-record mean of 1 − R(O ∪ ĵ)/R(O ∪ j_G1) with the same predictor in both arms; "floor" = the 0.29 %
MDE of Gate 2 at 512 clusters.

### 1.1 Paper I prioritisers (pathway-ambiguity programme, pre-G1)

| Selector | Definition | Implementation / tuning | Evaluation | Result | Admissible claim |
|---|---|---|---|---|---|
| `gnc_commutator` | operator-commutator surrogate of an oracle re-solve (ambiguity reduction) | fixed; no CV | 5 public cohorts, held-out excess gain in log-loss per pick, frozen protocol v1 + A1; confirmatory cohorts ST001849, ST002829 | **NOT_SUPPORTED**; paired difference vs difference-in-means negative in 5/5 cohorts, no interval excludes 0 [CERTAIN] | None for acquisition value. Surrogate fidelity to its own oracle only. |
| `uk_spectral`, `uk_summary` (U_k) | spectral / summary identifiability score | fixed | same | **NOT_SUPPORTED**; below random in 3–4 of 5 cohorts; ST000356 intervals below 0 [CERTAIN] | Negative result; oracle-tied results withdrawn. |
| `diffabundance`, `surrogate` | univariate label statistic and a prioritiser built on it | fixed | same | only prioritisers with CI > 0 (ST002829: +0.0029, +0.0033 nats) [CERTAIN] | Where there is value, a univariate label statistic captures as much as any structural prioritiser. |

Hindsight headroom in these cohorts: 0.014–0.054 nats. The whole question "which one metabolite to add within a
pathway" has little held-out value in these data.

### 1.2 The Gaussian baseline

| Selector | Definition | Implementation / tuning | Evaluation | Result | Admissible claim |
|---|---|---|---|---|---|
| **G1** | argmax_j c_jᵀc_j / v_j, conditional target–candidate covariance and candidate variance given O | z-scoring on training rows; S* = (1−α)S + αμI, α = min(1, max(0.05, d/(d+n))), d = 60; Schur complements, ridge 1e-6. **No tuned hyperparameter**; α is a fixed formula | #1, #2, #4, #5, #6, #7, #9, #10 above; re-estimated per training fold from Gate 3A on | best or tied with hindsight-best archived policy in 91/186 records; own gain over no measurement mean 3.1 %, median 1.9 % (bounded) [CERTAIN, Gate 2] | Population one-step optimum under Gaussian conditional model and isotropic squared loss (standard). On this archive no challenger has beaten it by the floor. |

Note on "properly tuned G1": the shrinkage weight α is a formula, not a CV-tuned quantity. No record opened here
tunes α or compares alternative covariance estimators (Ledoit–Wolf, graphical lasso, factor models) for G1 itself.
That is a gap in the baseline, not in the challengers [CERTAIN that no opened record does it; UNKNOWN whether an
unopened one does].

### 1.3 Paper III / IV criteria

| Selector | Definition | Implementation / tuning | Evaluation | Result | Admissible claim |
|---|---|---|---|---|---|
| `variance`, `E_optimal` | candidate conditional variance; smallest-eigenvalue criterion | fixed | #1, #6 | bounded mean relative to G1: −3.5 %, −0.82 % [CERTAIN, Gate 2] | Worse than G1. |
| **RIG** (= one-step target D-optimal; Paper III `D_optimal`) | log-det information gain on the target block; same rank-one update as G1 under a different metric | fixed; `RIG_reshrink` adds residual-correlation shrinkage λ = min(1, max(0.05, p/(p+n))) | #1, #2 (prespecified primary, both arms pre-existing), #6 | Paper IV primary: RIG superiority **not supported** (primary contrast 0.019297, reproduced in Gate 2); bounded mean +0.13 % (RIG), +0.14 % (reshrink) — half the floor [CERTAIN] | Theoretical: G1 and RIG score the same update; strict disagreement certifies Loewner incomparability (Paper IV, SUPPORTED). Empirical: no superiority. Held-out primary must be labelled "prespecified primary on previously evaluated outcomes", not preregistered [CERTAIN, H20]. |
| Paper IV directional-geometry prediction | sign of the G1–RIG disagreement predicted from a geometric witness | — | P4U (ST004571) and Stage 10K: supportive external tests with documented eligibility deviations | **supportive only**, not confirmatory [SUPPORTED, memory + corrections record]; Track B (P4V) was to test it on independent studies and is at EXECUTION NO-GO [CERTAIN] | The only external test attempted in the programme; it concerns directional geometry, not average superiority. |

### 1.4 Design-only selector gates (9 Oct)

| Selector | Definition | Implementation / tuning | Evaluation | Result | Admissible claim |
|---|---|---|---|---|---|
| AOPT (joint A-optimal, B = 2, 3) | exhaustive joint trace reduction | same covariance as G1 | #4, design rows, 4 outer folds | +0.03 % / −0.04 % bounded; gate (≥ 1 %, LB > 0) **failed**; greedy captures ≥ 96 % of the joint objective in every block (B23) [CERTAIN] | Joint optimisation adds nothing over greedy G1 at B ≤ 3 on this archive. |
| NESTED (inner-CV choice among candidate sets) | pick set with lowest inner-CV loss | linear ridge, λ ∈ {0.01, 0.1, 1} | #4 | +0.04 % [−0.51, +0.72] bounded at B = 1; **failed** at all budgets [CERTAIN] | Data-driven set selection does not convert apparent advantages (APPARENT median d = 0.14) into out-of-fold gain. |
| Quadratic learner (NBRA pre-check) | [O, m, m², O·m] vs linear | λ ∈ {0.01, 0.1, 1} | #3 (descriptive) | quadratic beats linear in 12 % of candidates, median −8.2 %; held-out RMSE diff vs tuned G1 −0.0066 [−0.0259, +0.0049], favouring G1 [CERTAIN] | Restricted nonlinear acquisition signal not estimable; does not grow with n. |
| **GDOA-0** (gate for NCPG / synergy-aware selection) | nonlinear acquisition excess v_nl − v_lin with all 20 candidates (bounds every budget by monotonicity) | kernel ridge K_lin + g_a K_add + g_f K_rbf, 7-point grid × λ ∈ {0.01, 0.1, 1}, inner CV | #5 | **FAILED** both populations: all-records bounded +0.018 (LB +0.004) but primary LB −0.095; n ≥ 200 bounded −0.005 [CERTAIN]. Excess concentrated in n < 50 records and shrinks with n. | No estimable nonlinear or synergy acquisition value on this archive. NCPG and interaction-aware selection were never built because their premise failed at the gate. |
| GDOA-2 certified selector (Hoeffding / paired-z) | certify an alternative over G1 when D̂(S) exceeds u | arithmetic only | none | u = 8–30× baseline MSE under Hoeffding; paired-z needs n_v ≈ 200–600 per study; 2/186 records have n ≥ 600 [CERTAIN] | Would return G1 in essentially every study. Not run. |

Descriptive lead from GDOA-0 (reported, not gated): the nonlinear learner beats the linear learner **on O alone**
(median +5.5 % bounded, positive in ~90 % of records). This is a predictor finding, not an acquisition finding.

### 1.5 Gate 2 and Gate 3A (19 Oct)

| Item | Result | Grade |
|---|---|---|
| Gate 2 hindsight headroom over G1 (best of 6 archived policies, evaluation cells) | 1.39 % [0.89, 1.99]; median 0.056 %; biased upward by hindsight selection | CERTAIN |
| Detectable effect on the Gate D pool | MDE = 0.29 % (K = 512), 0.33 % (K = 392), 0.47 % at K = 200 | CERTAIN |
| Gate 2 verdict | PASS (necessary, not sufficient) | CERTAIN |
| Gate 3A predictor diagnostic Δ_pred vs P1 (log-linear ridge) | P0 (raw ridge) −8.7 % [upper −5.8 %]; **P2 (quadratic expansion) −18.0 % [upper −14.9 %]**; P3 (kernel ridge) **+9.14 %** (recorded 95 % bound +10.7 %; interval excludes 0) | CERTAIN (`gate3a_summary.json`, sha256 2166d914…, located 3 Nov in `../decision_identifiability/gate3a_out/`, copied into the repo) |
| Gate 3A acquisition qualification Δ_acq (≥ 0.29 % and Bonferroni 98.33 % LB > 0) | P0 −0.545 % [LB −2.14 %]; **P1 +0.003 % [LB −0.36 %]**; **P2 −0.433 % [LB −1.69 %]**; **P3 −0.108 % [LB −0.58 %]**; `qualifying_families: []` | CERTAIN |
| Gate 3A diagnostics | A (hindsight headroom, upper bound): 4.3–6.5 % per family; B (split-half Spearman of V̂, median over blocks): P0 0.31, P1 0.41, P3 0.42, **P2 0.10**; C (nonlinear-only realised value vs P1): P3 +0.07 % [LB −0.49 %], P2 −7.3 %; leaky-selection control positive for every family (+3.2 to +13.9 %), so the test can see leakage; `evaluation_rows_read: false` | CERTAIN |
| **G1 self-agreement across folds** | the G1 pick re-estimated on 4/5 of the design rows equals the full-design G1 pick in **3,376 of 4,760 block-folds (71 %)** | CERTAIN |
| N0 regularization parity | the ridge had chosen λ = 10 (grid edge) in 92 % of folds; at a wide grid the kernel edge is +4.25 % [2.30, 5.93] squared; linear-kernel control = ridge to 1e-15 | CERTAIN |

Reading: half of the "nonlinear prediction advantage" was a regularization ceiling in the baseline. The half that
remains (+4.25 %) is a predictor result on the measured panel; it did not translate into acquisition (P3 tied G1).

### 1.6 NOVA family (19–20 Oct)

| Selector | Definition | Implementation / tuning | Evaluation | Result | Admissible claim |
|---|---|---|---|---|---|
| **NOVA-V** | argmax_j V̂_j, V̂_j = cross-fitted estimate of E‖E[T∣O,X_j] − E[T∣O]‖² via orthogonal score ψ̂ = 2d̂ᵀ(T − m̂₀) − ‖d̂‖² ≡ L_g − L_j | RBF kernel ridge on log-rule inputs, bandwidth {0.5,1,2}× median, ridge grid to 1e4, inner CV; 5 outer folds seed 20261019; additive and residualised-kernel nuisances as ablations | N1 synthetic (320 blocks); N2 archive #9 | N2: **−0.282 % [−0.724, +0.066]** (RBF); −0.920 % (additive). Agreement with G1 35.6 % [CERTAIN] | Reduces to G1 under joint Gaussianity (V_j = G1_j, population; N1/X0). **Loses to G1 on the archive.** Algebraically the orthogonalised score *is* the cross-fitted loss difference (Lemma 1, ESTIMAND audit F1) and Gate 3A's P3 *is* NOVA-V under the bounded loss and narrow grid — one estimator, two names. |
| NOVA-U | argmax_j (V̂_j − se_j) | same | N1, N2 | −0.209 % [−0.846, +0.247] [CERTAIN] | Uncertainty penalty alone does not rescue the estimator. |
| **NOVA-S** | NOVA-V's pick iff the paired per-row score difference against G1's pick exceeds 2 SE; else G1 | same; τ = 2 fixed before N2 | N1 (calibration: 1 spurious claim / 320), N2 (**pre-registered primary**), replay #10 | **Δ_acq = +0.112 % [+0.029, +0.221]**; departs in 10.5 % of block-folds (null rate 5.0 %); capped-loss sensitivity +0.135 % [+0.056, +0.243]; **below the 0.29 % floor; does not qualify** [CERTAIN] | First selector in the programme whose development interval excludes 0; the gain is produced by the gate removing harmful NOVA-V decisions, not by finding good ones (replay reading 4). Power on the Gate D pool for a true 0.11 %: 0.29 (one-sided 5 %). Development finding on a 9-times-used archive; protocol's own prior expectation was a tie. |
| NOVA-S threshold envelope (hindsight) | NOVA-S at τ ∈ {0, …, 4} from the frozen block-fold losses | none (replay) | #10 | max Δ_acq = **+0.158 % at τ = 1.5**; frozen τ = 2 not the optimum; no τ reaches the floor [CERTAIN] | Bounds what *any* departure rule based on the paired statistic t can achieve on this archive: ≤ 0.16 %, hindsight. |
| Replay heterogeneity (primary's estimand) | realised gain of NOVA-V over G1 by t bin | — | #10 + estimand audit | [2,3): +0.89 % [+0.34, +1.60]; [3,4): +0.87 %; ≥ 4: +5.69 % (13 block-folds); Spearman(t, gain) 0.16 [CERTAIN] | Pre-selection evidence predicts realised advantage (descriptive). The heterogeneity is real but thin; exploiting it needs a better pre-selection statistic, not a better threshold. |
| **NOVA-X** (certified) | depart from G1 only if a simultaneous one-sided lower bound on D_{j,g} is > 0 (multiplier bootstrap, q ≈ 2.2) | B = 2,000 | X0/X1 synthetic only; X2 on the archive **declined** | never more accurate than NOVA-S, sometimes worse (planted square regret 0.106 vs 0.035); simultaneous coverage 0.80–0.85 on the interaction DGP vs Bayes value; at or near nominal vs the fitted-procedure advantage (X1b) [CERTAIN] | Same mechanism as NOVA-S at a stricter threshold. The certificate covers the deployed decision (procedure advantage), not the population Bayes advantage; the nuisance bias is the obstruction. |
| Claim "H_j = V_j − G1_j ≥ 0" | — | — | X0 | **false** (linear-proxy counterexample: V₀ = 0, G1₀ = 5.71) [CERTAIN] | Withdrawn. Replaced by the exact signed decomposition V_j − G1_j = A_j − A_0 (X0b; Pythagoras; theorem-grade but standard). |
| N2 failure-regime diagnostic (archive use 11) | WLS of record-level realised gain D on t̄ and five covariates (log n, panel class, raw-scale skew proxy, N0 O-only kernel edge, Gate 3A G1 fold-stability); study-cluster bootstrap, 99 % Bonferroni | spec v0.3 `0d36d518…`, script v0.1.1 `f3e11962…` (v0.1.0 crashed before output, disclosed); reproduction gate +0.11179 % passed | one run, 179 records | **No regime identified.** `edge_O` −0.0005 [−0.0021, +0.0018]; `g1_stability` −0.0023 [−0.0058, +0.0007]; `log_n` +0.0018 [+0.0001, +0.0045] but high-n subgroup Δ_acq +0.083 % (LB −0.019 %); low-n tertile D −0.0028 [−0.0048, −0.0005]. A4′ inside-G1set D +0.0068 [+0.0045, +0.0096] **confounded** (reference set saw scored-fold targets) [CERTAIN] | Departures' harm is concentrated in small records and vanishes in larger ones; no record size shows a qualifying benefit; no pre-acquisition covariate locates a regime. The tie-breaker idea is a hypothesis with a named confound, not a finding. |
| L0 representation hierarchy | V^L ≤ … via nested classes L ⊆ K ⊆ H | kernel level must be cumulative (1 + zᵀz′ + rbf(O) + rbf(O,x_j)) or inclusions fail | synthetic | K-level has a capacity component (V^K = 0.80 for pure-noise candidates on the linear-proxy design) [CERTAIN] | The kernel level is a diagnostic, not a score. |

### 1.7 Theory packages reviewed 11 Oct (project record `uapa_f1_result.md`)

Assay-aware information / learnability / transport; Gaussian Assay Experiment Theory v0.2; Assay-Channel
Identifiability v0.3; Finite-Sample Acquisition Theory; Paired Acquisition Theory. Reviewer's verdicts [CERTAIN as
recorded]: proofs correct; **all standard** (Le Cam deficiency, L² projection, inverse-Wishart OLS excess risk,
Davidson–Szarek, Reiersøl/Fuller errors-in-variables identification, Pitman–Morgan, Hoeffding, Bonferroni);
bounds vacuous at archive sizes (G3 needs n > 68–80; I7 needs n ≥ 1,344 even at κ = 2); the proposed "next
advances" already exist (White 2000 Reality Check, Hansen 2005 SPA, Model Confidence Set, Bradley–Blackwood,
Diebold–Mariano). Under the Gaussian equal-dimension model E[deployed Δ] = γ(n−1)/(n−d−1) has the sign of γ, so
F1's ≈ 0 means γ ≈ 0, not an advantage hidden by learning error. DEP, DEP-U, DEP-MS, DEP-GM, DIAF freezes: not
opened here [UNKNOWN]; the 24 Oct decision records them as Track A theory with "novelty not claimed from freezing".

## 2. Novelty classification (Part II of the directive)

| Component | Class | Basis |
|---|---|---|
| G1 = conditional covariance variance reduction; population optimum under Gaussianity | established (standard Gaussian conditioning / A-optimal one-step) | — |
| V_j = E‖E[T∣O,X_j] − E[T∣O]‖² = Bayes-risk reduction under squared loss | established identity (law of total variance) | directive itself; X0 |
| V_j = G1_j under joint Gaussianity | specialisation (standard) | N1, X0 |
| V_j − G1_j = A_j − A_0 (signed projection decomposition) | correct, Pythagorean; not new as mathematics; useful as a diagnostic | X0b, L0 |
| Orthogonalised cross-fitted score for V_j | specialisation of Neyman-orthogonal / DML variable-importance estimation (Williamson et al.; Zhang & Janson "Floodgate") | GDOA0 spec names the lineage; ESTIMAND audit F1 shows ψ̂ ≡ L_g − L_j |
| G1 ≡ RIG as the same rank-one update; Loewner incomparability certificate | SUPPORTED as a theoretical observation of Paper IV; prior-art status not checked here | UNKNOWN |
| NOVA-S / NOVA-X gating (depart from a default only on a paired test) | computational procedure; the inferential pieces are Pitman–Morgan / Diebold–Mariano / SPA-type | 11 Oct review |
| Certification for the procedure advantage vs Bayes advantage (X1b) | a clarification of *which* estimand a bootstrap certificate covers; potentially a small methodological point | X01 addendum |
| Kernel ridge, cross-fitting, interaction modelling, alternative (bounded) losses | **all already tested** under the names P3 / GDOA-0 kernel / quadratic learner / bounded loss | §1.4–1.6 |
| Distributional robustness (covariate shift) | tested only synthetically (N1 covariate-shift DGP: G1 regret 1.56, NOVA-S 0.93–1.02); not on real data | N1 |
| Empirical claim: any selector beats G1 on held-out prediction | **unsupported** on the archive; **untested** externally | all of §1 |

## 3. What the evidence says, in five sentences

1. **No selector in the programme has beaten G1 at the floor on any data, and none has been tested on external data.**
   The best development result is NOVA-S at +0.11 % (LB +0.03 %), one third of the smallest effect the sealed pool
   could detect; the hindsight ceiling for its whole mechanism is +0.16 %.
2. **The nonlinear premise was tested four ways and failed each time on acquisition** (quadratic pre-check, GDOA-0,
   Gate 3A P2/P3, NOVA-V), while succeeding on prediction of the measured panel (+4–5 % after regularization parity).
   The archive supports "better modelling of O" and does not support "better choice of X_j".
3. **The failure regime of G1 that the records document is not "nonlinear relationships" but three other things:**
   (i) **G1's own selection instability**: refit on 80 % of a record's design rows, G1 changes its pick in 29 % of
   block-folds (Gate 3A, 3,376/4,760 agreement), so a large part of the between-arm noise in every comparison is G1
   disagreeing with itself, not a challenger finding something; (ii)
   a linear-proxy pathology (G1 assigns high value to a candidate carrying information O already determines
   nonlinearly; X0 counterexample, L0) and (iii) covariate shift (N1); neither (ii) nor (iii) has been shown to occur in
   the real archive at a detectable rate. Gate 3A also shows the quadratic expansion (P2) is *harmful* for prediction
   (−18 %) and acquisition (−0.43 %), which is the clearest finite-sample evidence in the programme that added model
   complexity worsens generalisation at these n. GDOA-0's bounded excess shrinking with n argues the apparent small-n nonlinear
   excess is estimation noise.
4. **The archive is exhausted as a selector-development resource.** It has been used ten times; every later use
   must be labelled development; the one untouched resource (Gate D, 512 clusters) can detect 0.29 % and nothing
   proposed reaches it.
5. **Archive closed for selector development (owner, 3 Nov 2026)** after use 11. The remaining route is external: a
   frozen prediction (`EXTERNAL_VALIDATION_PREREG_STUB_v0.1.md`) and the LIPID/COMETS access decision on 15 Dec 2026.
6. **The memo's hybrid (1−α)Ĝ1 + αV̂_nl is a soft relaxation of NOVA-S's hard gate.** Its reachable gain on this
   archive is bounded by the replay envelope unless it uses information other than the paired statistic t; a run of
   it on the archive would be the eleventh use and could not qualify under the Gate 2/3 rules the owner froze.

## 4. The exact empirical claim still requiring falsification

> On independently recruited human cohorts with participant-level sampling structure established, with G1 and the
> challenger given identical design observations and comparable tuning budgets, the challenger's equal-study mean
> relative reduction in held-out squared prediction loss at budget 1 is ≥ 0.29 % (or the MDE recomputed for the
> actual number of study units) with a lower confidence bound above 0.

No existing selector is a credible candidate for that test at present. The only development evidence that points
anywhere is the replay heterogeneity (pre-selection evidence predicts realised advantage), and the only external
test ever run (P4U / Stage 10K) concerned a different, directional claim and was supportive at best.

## 5. Open items that would change this matrix

- ~~Open the Gate 3A summary JSON~~ Done 3 Nov 2026 (§1.5). The summary, run manifest and secondary now live in `gate3a_out/` in the repository; `gate3a_records.csv` and `gate3a_blocks.jsonl` stay local per protocol §6.
- Confirm from the Paper IV freeze whether the "G1 ≡ RIG rank-one update" statement has prior art (Loewner-order
  comparisons of information criteria are classical in optimal design; Pukelsheim 1993 is the obvious reference).
- Record whether any frozen note tunes G1's shrinkage (α) or covariance estimator; if none does, "properly tuned
  G1" has not yet been established and is the cheapest missing control.
