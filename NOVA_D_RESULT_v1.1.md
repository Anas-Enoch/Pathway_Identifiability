# NOVA-D — Detectability of Nonlinear Acquisition Advantage: result record v1.1 (7 Nov 2026)

**Correction log (v1.0 → v1.1; frozen outputs untouched; v1.0 retained on disk).** Owner review of 7 Nov 2026:
(1) v1.0's outcome box stated that every positive signal occurs at n = 128, while its own §4 and §7.3 list
condition 59 at n = 512 — corrected throughout, and the attribution to comparator degradation is now established
cell by cell (§7.3, new table) rather than inferred. (2) "500 block-folds" is now stated as 100 replicates × 5
outer folds; folds within a replicate share a training sample and are not independent experiments. (3) The G1
shrinkage-floor finding is reframed from "defect" to the owner's distinction between population ranking
consistency (Objective A) and finite-sample decision loss (Objective B) (§7.4, §8.3). No number from the frozen
analysis changes.

**Outcome: both prespecified claims fail. H_D (the near-tie rule) is refuted on its own grid. At the primary
condition and at the positive control the frozen NOVA-S made zero departures from G1 in 500 block-folds each
(κ = 0 exactly; 100 replicates × 5 outer folds); the ungated NOVA-V recovered κ = 0.008 and κ = −0.007
respectively. No near-tie threshold improved on unscreened NOVA-S in any of the 72 conditions, and tight thresholds
hurt in 18 condition–γ pairs. The validity controls are clean for NOVA-S. Eight of the nine nonlinear conditions
in which NOVA-S has a positive lower bound are at n = 128 (condition 59, n = 512, is the exception); in every one
of the twelve positive-signal cells, checked individually, the G1 comparator's oracle regret exceeds NOVA-V's, so
the apparent advantages are comparator regret, not nonlinear recovery. A post-hoc diagnostic identifies a
mechanism in the frozen G1's shrinkage rule (§7.4).**

Author of record: Anas Enoch. Compiled by Claude (Fable 5.1). Protocol `NOVA_D_PROTOCOL_v0.3.md` frozen
2026-11-05T21:11:39Z (commit `a98890d`); run 5–7 Nov 2026 on the owner's machine, four partitions, ~36 h wall.

## 1. Provenance

| Item | SHA-256 |
|---|---|
| `NOVA_D_PROTOCOL_v0.3.md` | `df13d06496ce59a5319a6c614a815d5a1d316c1e6538fd699eda5d2dbfaf5d20` |
| `nova_d_run.py` v0.1.0 | `f2359f6b43000fadb68a8728b979731110257c7403b97ce8b3ef0eaee156ec52` |
| `nova_d_calib.py` | `4ab49dbc585e5cef4da6d0f2a31ca0d44725ffbaf0364f058fd8b63615b1ce29` |
| `p6_n2.py` v0.1.0 (NOVA-V/U/S, unchanged) | `681ce852609e6ab4fcaf0f0ff603477c64e2351e10bd119811d45f008d4bbe8e` |
| `p5_gate3a.py` (G1, unchanged) | `4600e2f3e3bb652affadcc6f6935380d8b9a77a4237e7cb306f056f849478130` |
| `nova_d_out/FREEZE.json` | `1cfa42dd61b6e8a19d319d833eb05884390b7f553567da874cbf49c38762adf8` |
| `nova_d_out/replicates.jsonl` (7,200 records, 10,694,157 bytes) | `e67719b59d7312de3a4275c5c5f11a25a6d5dcbf396a956b1d09e32e675e4f35` |
| `nova_d_out/integrity.json` — PASS: 7,200 parsed, 7,200 distinct keys, 0 duplicates, 0 missing, 0 inconsistent | `ae1b3179e4b766ba9e3f06cc0854fac6e52bf3e0abc38ad3dd743cf88cbdcd41` |
| `nova_d_out/nova_d_result.json` (frozen `analyze`) | `29b25189d129240b20fae06cdbe320f46fa8a03d28bf6a3dbe5f56ce8e46de53` |

Grid as frozen: 4 mechanisms, 8 mechanism–advantage combinations, n ∈ {128, 512, 1024}, A ∈ {0.05, 0.15, 0.20},
72 conditions × 100 replicates, 5,000 evaluation rows per replicate, γ ∈ {0.02, 0.05, 0.10, 0.25, ∞} as a
policy. Nothing was changed after the freeze. The `--partial` path of `analyze` was never invoked.

## 2. Prespecified decisions (verbatim from `nova_d_result.json["decisions"]`)

```
positive_control_regime_viability : false
primary_recoverability            : false
validity_M1_M3_breach             : []
near_tie_effective_conditions     : []
grid_max_rho_NOVA_S_upward_biased : [0.0478, condition 29]   (M2, rho_pop 0.10, n 128, A 0.20)
```

## 3. Primary and positive control

| | Primary — M2, ρ_pop 0.10, n 512, A 0.15 (idx 31) | Positive control — M2, ρ_pop 0.40, n 1024, A 0.15 (idx 43) |
|---|---|---|
| Criterion | LCB₉₅(ρ_NOVA-S) > 0 | LCB₉₅(ρ_NOVA-S) > 0.10 (κ > 0.25) |
| NOVA-S ρ, LCB₉₅, κ | 0.0000, 0.0000, **κ = 0.000** [0, 0] | 0.0000, 0.0000, **κ = 0.000** [0, 0] |
| NOVA-S departures (100 replicates × 5 folds) | **0 / 500** | **0 / 500** |
| NOVA-V ρ, LCB₉₅, κ | 0.0008, −0.0010, 0.008 | −0.0027, −0.0053, **−0.007** |
| NOVA-U κ | 0.015 | −0.006 |
| Realised A (evaluation rows) | 0.144 | 0.149 |
| Opportunity rate P(ĵ_G1 ≠ j\*) | 0.404 | 0.326 |
| Oracle regret /R_∅: G1, NOVA-S, NOVA-V | 0.0052, 0.0052, 0.0070 | 0.0166, 0.0166, 0.0245 |
| Population approximation gap (vs Gaussian j_G) | 0.0127 | 0.0510 |
| Screened policies, all γ | identical to NOVA-S (no departures to gate) | identical to NOVA-S |
| M4-style exclusion of opportunities at γ = 0.02 | 0.92 | 0.74 |

**Why there were no departures (from the run log, read-only).** In the positive control NOVA-V and G1 disagreed in
83 of 500 folds; across those 83 the paired statistic had **max t = 1.03** (median 0.21, 99th percentile 0.99)
against the gate τ = 2. In the primary, 75 disagreements, max t = 1.19. The gate was never within reach.

**Which selector found the optimum.** Positive control: G1 picked j\* in **67.4 %** of folds, NOVA-V in **52.0 %**.
Primary: G1 59.6 %, NOVA-V 45.4 %. The cross-fitted kernel scores identified the nonlinear-optimal candidate less
often than the linear criterion did; NOVA-V's oracle regret exceeds G1's in both conditions.

## 4. The whole grid

Fifty-four conditions carry a nonlinear advantage (M2, M4). Across them:

- NOVA-S made **no departure at all in 37 of 54**; its maximum departure rate was 0.09 (condition 29, n = 128).
- NOVA-S attained LCB₉₅(ρ) > 0 in nine conditions: 20, 27, 29, 36, 45, 54, 56, 64 (all **n = 128**) and 59
  (n = 512, κ = 0.030). Its maximum κ was 0.478 (condition 29, n = 128).
- NOVA-V: κ ranged from −0.94 to 2.09, median **0.007**; ρ < 0 (worse than G1) in **21 of 54**. Four conditions
  show κ > 1 — 27 (1.00), 29 (2.09), 54 (1.07), 56 (1.47) — **all at n = 128**. κ > 1 is impossible for genuine
  recovery against the population comparator (E[ρ_policy] ≤ ρ_pop); see §7.3.
- At n ≥ 512 in M2 and M4: NOVA-S κ ≤ 0.030 everywhere; NOVA-V κ has median ≈ 0.01 and is negative in half the
  conditions.

Eighteen null conditions (M1, M3; H ≡ 0):

- **NOVA-S: clean.** Maximum ρ = 0.0071 (condition 11, n = 128); no condition approaches 0.10. The validity
  control passes.
- **NOVA-V in nulls** reaches ρ = 0.118 with LCB₉₅ > 0.10 in condition 3 (M1, n = 512, A = 0.05) and 0.058 in
  conditions 11 and 12. Where no nonlinear structure exists, the ungated nonlinear estimator "beats" G1 by the
  same magnitudes it shows in M2/M4. Its positive conditions in M2/M4 are therefore not evidence of nonlinear
  recovery. (The validity rule in protocol §10 was written for NOVA-S, the frozen primary; this reading of
  NOVA-V is the control doing the job it was designed for.)

Near-tie rule (H_D):

- No (condition, γ) pair has a paired bootstrap interval for ρ_screened − ρ_unscreened with lower bound > 0.
- Eighteen pairs have **upper bound < 0** — screening hurt — all in conditions where unscreened NOVA-S had a
  positive ρ (n = 128): tightening γ removed the departures that were helping.
- M4 aim control: the screen at γ = 0.02 excluded **74–100 %** of beneficial opportunities in every M4 condition,
  and 57–94 % in M2. The rule is mis-aimed where the nonlinear optimum is away from the Gaussian frontier, as
  the protocol anticipated, and ineffective where it is at the frontier.

## 5. Answers to the three protocol questions

1. **Existence.** By construction, every M2/M4 condition contained a population nonlinear advantage of the stated
   ρ_pop; the sampler was verified against the closed forms (self-test). Existence is not in question.
2. **Recoverability.** The frozen NOVA-S estimator recovered **none** of it at n ≥ 512. The frozen NOVA-V
   recovered a median 1 % of it and was worse than G1 in half the conditions. At the positive control — a 40 %
   population advantage, 34 % value excess on one candidate, n = 1024 — the paired statistic never exceeded 1.03.
3. **Selection.** Near-tie gating did not increase useful departures anywhere and reduced them where any existed.

## 6. What this does and does not establish

- **Established, within this family and grid:** the frozen NOVA-S does not recover a known nonlinear acquisition
  advantage at n ≤ 1024 in a favourable (H ≥ 0) family; the near-tie policy does not help and can hurt. Per
  protocol §4, failure in a favourable family is strong evidence against the tested strategy *within this family
  and grid*; it is not a universal refutation of nonlinear acquisition (owner's qualification, recorded).
- **Not established:** behaviour at n = 2048 or beyond (owner's recorded caveat: the positive control's failure at
  n = 1024 does not establish failure at n = 2048); behaviour when the signed correction H is negative;
  behaviour of estimators other than the frozen ones; anything about biological cohorts.
- **Consistent with every prior programme result:** archive NOVA-S +0.112 % below the floor (use #9); Gate 3A no
  family qualified; N2 "no regime identified"; NHANES-FA 0/25 departures. NOVA-D adds the mechanism: the gate's
  statistic has signal-to-noise below 1 even where a large advantage exists.

## 7. Post-hoc diagnostics (not prespecified; labelled as such)

### 7.1 The gate statistic
See §3. At n ≥ 512, t < 1.2 in every disagreeing fold of the primary and positive control. A τ anywhere above
≈ 1 would also have produced zero departures; this is not a threshold-tuning problem.

### 7.2 Identification of the optimum
G1 identified j\* more often than NOVA-V at n ≥ 512 in both designated conditions (§3). The nonlinear estimator's
disadvantage is in locating the candidate, not only in certifying it.

### 7.3 The positive-signal cells, checked individually: comparator regret, not recovery
ρ is measured against the *estimated* G1's realised loss, so when the estimated G1 collapses ρ can exceed ρ_pop;
the bound E[ρ_policy] ≤ ρ_pop holds against the population comparator only. Every cell with a positive signal was
checked against the run log (read-only): the nine nonlinear conditions with LCB₉₅(ρ_NOVA-S) > 0, and the three
nulls with NOVA-V ρ > 0 (3, 11, 12). "#G1" = number of distinct candidates G1 picked across the 500 fold-decisions.

| idx | mech | ρ_pop | n | A | κ_S | κ_V | G1→j\* | V→j\* | #G1 | top G1 picks | regret G1 | regret V | ratio |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 3 | M1 | 0 | 512 | 0.05 | — | — | 0.17 | 0.54 | 17 | 3:103, 10:96, 0:86 | 0.0099 | 0.0054 | 1.85 |
| 11 | M3 | 0 | 128 | 0.20 | — | — | 0.23 | 0.51 | 3 | 1:323, 0:117 | 0.0114 | 0.0042 | 2.73 |
| 12 | M3 | 0 | 512 | 0.05 | — | — | 0.07 | 0.29 | 12 | 3:202, 5:67, 4:59 | 0.0051 | 0.0027 | 1.89 |
| 20 | M2 | 0.05 | 128 | 0.20 | 0.038 | 0.22 | 0.29 | 0.47 | 3 | 0:342, 1:143 | 0.0080 | 0.0051 | 1.56 |
| 27 | M2 | 0.10 | 128 | 0.05 | 0.187 | 1.00 | 0.02 | 0.10 | 17 | 7:108, 6:97, 14:56 | 0.0153 | 0.0116 | 1.32 |
| 29 | M2 | 0.10 | 128 | 0.20 | 0.478 | 2.09 | 0.15 | 0.48 | 3 | 0:232, **2:194**, 1:74 | 0.0433 | 0.0098 | 4.41 |
| 36 | M2 | 0.40 | 128 | 0.05 | 0.013 | 0.22 | 0.10 | 0.16 | 14 | 7:114, 9:82, 5:65 | 0.0256 | 0.0219 | 1.17 |
| 45 | M4 | 0.05 | 128 | 0.05 | 0.194 | 0.60 | 0.03 | 0.08 | 17 | 0:166, 2:101, 16:57 | 0.0112 | 0.0097 | 1.16 |
| 54 | M4 | 0.10 | 128 | 0.05 | 0.182 | 1.07 | 0.00 | 0.02 | 18 | 14:125, 4:74, 0:70 | 0.0165 | 0.0139 | 1.18 |
| 56 | M4 | 0.10 | 128 | 0.20 | 0.329 | 1.47 | 0.00 | 0.00 | 4 | 1:340, 0:154 | 0.0449 | 0.0254 | 1.77 |
| **59** | M4 | 0.10 | **512** | 0.20 | 0.030 | 0.86 | 0.00 | 0.00 | 2 | 0:301, 1:199 | 0.0319 | 0.0178 | 1.80 |
| 64 | M4 | 0.40 | 128 | 0.15 | 0.119 | 0.44 | 0.00 | 0.00 | 6 | 1:306, 0:150 | 0.0765 | 0.0594 | 1.29 |
| *ref* 31 | M2 | 0.10 | 512 | 0.15 | 0 | 0.01 | 0.60 | 0.45 | 2 | — | 0.0052 | 0.0070 | 0.74 |
| *ref* 43 | M2 | 0.40 | 1024 | 0.15 | 0 | −0.01 | 0.67 | 0.52 | 2 | — | 0.0166 | 0.0245 | 0.68 |

In all twelve positive cells G1's regret exceeds NOVA-V's; NOVA-V identifies the true optimum in ≤ 54 % of folds
and in 0–10 % in eight of them. In the weak-signal n = 128 cells G1 spreads its picks over 14–18 candidates. The
four κ > 1 values (27, 29, 54, 56) are all n = 128. Condition 59 (n = 512) is the same mechanism in a cleaner
form: M4 has a 20 % Gaussian gap, G1 chose the Gaussian runner-up in 40 % of folds, neither selector ever chose the
rank-7 nonlinear optimum, and NOVA-S's κ = 0.030 is NOVA correcting G1's rank-1/rank-2 confusion. In the two
reference rows, where G1 is stable, the regret ordering reverses. The attribution "apparent advantage = comparator
regret" is therefore established in every positive cell, not inferred from a subset.

### 7.4 Mechanism in the frozen G1: the shrinkage floor changes its population ordering
The frozen rule uses α = min(1, max(0.05, d/(d+n))) and shrinks the whole standardised covariance toward
(tr S/d)·I. Shrinking Σ_OO attenuates the conditioning on O, so base-panel-mediated *unconditional* covariance leaks
into the conditional score. Verified against 200,000-row samples (α at its 0.05 floor), comparing the frozen rule
with the same rule at α = 0:

| Condition | Population G ordering (top 4) | Frozen G1 at n = 200k, α = 0.05 | Same rule, α = 0 |
|---|---|---|---|
| 3 — null, A 0.05 | 0, 1, 2, 3 | **0, 3, 5, 10** (the four candidates G1 picked at n = 512) | 0, 1, 2, 4 |
| 43 — positive control | 0, 1, 2, 3 (gap 0.02) | **1, 0**, 2, 3 — the floor **reverses** the top pair; G1's limit pick is j\* | 0, 1, 2, 3 |
| 31 — primary | 0, 1, 2, 3 (gap 0.02) | 0, 1 with scores 0.7802 vs 0.7789 — gap compressed from 2 % to 0.17 % | 0, 1 (2 % gap kept) |

Consequences. (i) The frozen G1 is not a consistent estimator of the population Gaussian ranking: with the 0.05
floor, Σ̂_α → 0.95 Σ + 0.05 μI rather than Σ, the conditional covariance Σ_TP|O is a nonlinear function of the
blocks, and argmax_j G_j(Σ_0.05) can differ from argmax_j G_j(Σ) — the 200k-row rows above show that it does in
these distributions. **Whether that is a defect depends on the objective** (owner's distinction, adopted here):
for *Objective A*, estimating the population Gaussian ranking, the floor is an asymptotic bias and can yield an
inconsistent decision; for *Objective B*, minimising finite-sample acquisition loss, some asymptotic bias may be
worth its variance reduction, and a regularised selector can beat a consistent but noisy one at small n. NOVA-D
was not designed to test Objective B and this record does not decide it. What NOVA-D's own cells say about B is
suggestive only: in the weak-signal and small-n cells (§7.3) the heavily shrunk G1 (α up to 0.42 at n = 128) had
higher decision regret than NOVA-V, whose kernel-ridge nuisance uses no covariance shrinkage — a differently
regularised estimator did better there — but no cell compares shrinkage levels within the G1 rule itself. (ii) In the positive control the floor pulled
G1's own limit onto the nonlinear optimum, so the realised opportunity for NOVA-S **against the frozen comparator**
was smaller than the designed ρ_pop = 0.40: the finite-sample G1 reached j\* in 67 % of folds and NOVA-S never
corrected the remaining 33 %. The positive control is a valid failure of NOVA-S in that condition, but a weaker
test of "recovery against G1" than the design intended. The primary is unaffected in this respect (G1's limit
pick there is the Gaussian pick). (iii) This mechanism accounts for the n = 128 scatter in §7.3 and is consistent
with the archive finding that G1 re-estimated on 4/5 of the design rows agrees with its own full-design pick in only
71 % of block-folds (matrix §3.3): the comparator's instability is in part a property of its shrinkage rule. The
defensible research question is not "is the floor wrong" but **when covariance shrinkage improves finite-sample
acquisition decisions, and when its distortion of the conditional covariance induces persistent ranking errors.**

**Design oversight, disclosed (Claude).** The calibration solved the DGP against the population Gaussian ordering
G, and the self-test verified the sampler against the DGP's closed forms — not the frozen G1's population limit
against that ordering. The protocol's identity V_{j\*} − V_{j_G} = (H_{j\*} − H_{j_G}) − Γ uses j_G = argmax G;
the implemented comparator's limit can be a different candidate. This was correctable before the freeze and was
not caught. It does not alter the primary decision, the null controls, or the near-tie conclusion; it alters how
much the positive control tested.

## 8. Consequences for the programme

1. **H_D is refuted and the near-tie route is closed.** No further γ, no second grid.
2. **The external-validation hypothesis in `EXTERNAL_VALIDATION_PREREG_STUB_v0.1.md` now has a documented mechanism
   of failure, and lowering τ would not remove it** (NOVA-V's ranking is the deeper problem, §3): the NOVA-S gate statistic does not approach its threshold where a large nonlinear advantage is
   known to exist, at sample sizes at or above anything the programme's data routes can supply. The stub and the
   15 December LIPID/COMETS decision are unchanged by this record; what has changed is the prior that a cohort
   would be spent on a hypothesis with a known failure mechanism. Owner's position (7 Nov 2026): the stub and the
   15 December deadline remain unchanged; choosing not to proceed is a programme decision, not a modification of
   the preregistration; the existing evidence no longer offers a strong scientific justification for prioritising
   that validation, which does not make external failure certain.
3. **The most consequential finding is about the baseline, not the challenger — and it is a question, not yet a
   result.** The frozen G1's shrinkage floor changes its selection in exactly the near-tie and weak-signal regimes
   where acquisition decisions matter, at the population level (Objective A). Whether a different regularisation
   gives lower finite-sample decision loss (Objective B) is open. Owner's authorisation (7 Nov 2026): one new,
   strictly prospective programme centred on G1 and independent of NOVA-S, with **acquisition decision loss** as
   the primary outcome (not covariance-estimation error), distinguishing population ranking consistency,
   finite-sample decision regret and robustness to ill-conditioning, with three prespecified arms — a vanishing or
   zero floor; data-adaptive shrinkage including Ledoit–Wolf; conditional-first estimation (condition on O before
   regularising) — and the original G1 retained as a comparator. The archive stays closed; a new frozen
   simulation protocol or a genuinely independent dataset is required. The general principle that minimising
   covariance-estimation error need not minimise decision regret is known in decision theory; the contribution
   would be a specific, reproducible failure regime and, if one exists, an effective alternative.
4. **Matrix row (v0.5 §1.6):** "NOVA-D (synthetic, 72 conditions × 100 replicates, frozen code): primary κ = 0,
   positive control κ = 0, 0/500 departures each; max t 1.03 at n = 1024 with a 40 % population advantage; near-tie
   rule never helps, hurts in 18 pairs; NOVA-S clean in all 18 nulls; all positive signals at n = 128 are G1
   degradation; G1 shrinkage floor reorders candidates at the population level (diagnostic)."
5. **NOVA manuscript v0.5:** the limitations paragraph can now state the mechanism rather than the absence of
   evidence; the G1 floor finding belongs in its own subsection, labelled as a post-hoc diagnostic of NOVA-D.

## 9. Files to commit

`NOVA_D_RESULT_v1.1.md` (v1.0 retained); `nova_d_out/nova_d_result.json`, `nova_d_out/integrity.json`, `nova_d_out/selftest.json`
(the owner's machine's run, which the freeze checked), the four `run_*.log` files; and `nova_d_out/replicates.jsonl`
(10.7 MB, the primary output — commit it or, if size is a concern, commit a gzip and record both hashes).
