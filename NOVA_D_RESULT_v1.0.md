# NOVA-D — Detectability of Nonlinear Acquisition Advantage: result record v1.0 (7 Nov 2026)

**Outcome: both prespecified claims fail. H_D (the near-tie rule) is refuted on its own grid. At the primary
condition and at the positive control the frozen NOVA-S made zero departures from G1 in 500 block-folds each
(κ = 0 exactly); the ungated NOVA-V recovered κ = 0.008 and κ = −0.007 respectively. No near-tie threshold
improved on unscreened NOVA-S in any of the 72 conditions, and tight thresholds hurt in 18 condition–γ pairs. The
validity controls are clean for NOVA-S. Every positive signal in the grid occurs at n = 128 and is attributable to
degradation of the G1 comparator, not to nonlinear recovery. A post-hoc diagnostic identifies the mechanism of that
degradation in the frozen G1 itself (§7.4).**

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
| NOVA-S departures | **0 / 500** | **0 / 500** |
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

### 7.3 κ > 1 and the n = 128 "successes" are G1 degradation
In condition 29 (M2, ρ_pop 0.10, n 128, A 0.20), G1 picked candidate 0 in 46 % of folds, the nonlinear optimum
(candidate 1) in 15 %, and the Gaussian-rank-3 candidate in **39 %**; its oracle regret (0.043/R_∅) was four times
NOVA-V's (0.010). ρ is measured against the *estimated* G1's loss, so when the estimated G1 collapses, ρ can exceed
ρ_pop; the bound E[ρ_policy] ≤ ρ_pop holds against the population comparator only. All four κ > 1 conditions and
eight of the nine NOVA-S LCB > 0 conditions are n = 128 conditions with scattered G1 pick histograms. In the null
condition 3 (n = 512, A = 0.05) G1 picked the true optimum in 17 % of folds and spread the rest over candidates 3,
5 and 10 — the last with G ≈ 0.047, under 5 % of the top Gaussian value — which is why NOVA-V shows ρ = 0.118 in a
Gaussian world: it corrects G1, it does not find nonlinearity.

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

Consequences. (i) The frozen G1 is not the population Gaussian selector; with the 0.05 floor its large-n limit can
differ from G when the Gaussian gap is small or the signal weak. (ii) In the positive control the floor pulled
G1's own limit onto the nonlinear optimum, so the realised opportunity for NOVA-S **against the frozen comparator**
was smaller than the designed ρ_pop = 0.40: the finite-sample G1 reached j\* in 67 % of folds and NOVA-S never
corrected the remaining 33 %. The positive control is a valid failure of NOVA-S in that condition, but a weaker
test of "recovery against G1" than the design intended. The primary is unaffected in this respect (G1's limit
pick there is the Gaussian pick). (iii) This mechanism accounts for the n = 128 scatter in §7.3 and is consistent
with the archive finding that G1 re-estimated on 4/5 of the design rows agrees with its own full-design pick in only
71 % of block-folds (matrix §3.3): the comparator's instability is in part a property of its shrinkage rule.

**Design oversight, disclosed (Claude).** The calibration solved the DGP against the population Gaussian ordering
G, and the self-test verified the sampler against the DGP's closed forms — not the frozen G1's population limit
against that ordering. The protocol's identity V_{j\*} − V_{j_G} = (H_{j\*} − H_{j_G}) − Γ uses j_G = argmax G;
the implemented comparator's limit can be a different candidate. This was correctable before the freeze and was
not caught. It does not alter the primary decision, the null controls, or the near-tie conclusion; it alters how
much the positive control tested.

## 8. Consequences for the programme

1. **H_D is refuted and the near-tie route is closed.** No further γ, no second grid.
2. **The external-validation hypothesis in `EXTERNAL_VALIDATION_PREREG_STUB_v0.1.md` now has a documented mechanism
   of failure:** the NOVA-S gate statistic does not approach its threshold where a large nonlinear advantage is
   known to exist, at sample sizes at or above anything the programme's data routes can supply. The stub and the
   15 December LIPID/COMETS decision are unchanged by this record; what has changed is the prior that a cohort
   would be spent on a hypothesis with a known failure mechanism. Owner's decision.
3. **The most consequential finding is about the baseline, not the challenger.** The frozen G1's shrinkage floor
   changes its selection in exactly the near-tie and weak-signal regimes where acquisition decisions matter, and
   it does so at the population level. That is a reproducible property of the Paper IV rule, testable without any
   biological data, and it bears on every G1 comparison in the programme. Any test of alternatives (α without a
   floor; Ledoit–Wolf; conditioning before shrinkage) is a **new prospective protocol** — it must not be run on the
   closed archive, and it is not a NOVA-S rescue.
4. **Matrix row (v0.5 §1.6):** "NOVA-D (synthetic, 72 conditions × 100 replicates, frozen code): primary κ = 0,
   positive control κ = 0, 0/500 departures each; max t 1.03 at n = 1024 with a 40 % population advantage; near-tie
   rule never helps, hurts in 18 pairs; NOVA-S clean in all 18 nulls; all positive signals at n = 128 are G1
   degradation; G1 shrinkage floor reorders candidates at the population level (diagnostic)."
5. **NOVA manuscript v0.5:** the limitations paragraph can now state the mechanism rather than the absence of
   evidence; the G1 floor finding belongs in its own subsection, labelled as a post-hoc diagnostic of NOVA-D.

## 9. Files to commit

`NOVA_D_RESULT_v1.0.md`; `nova_d_out/nova_d_result.json`, `nova_d_out/integrity.json`, `nova_d_out/selftest.json`
(the owner's machine's run, which the freeze checked), the four `run_*.log` files; and `nova_d_out/replicates.jsonl`
(10.7 MB, the primary output — commit it or, if size is a concern, commit a gzip and record both hashes).
