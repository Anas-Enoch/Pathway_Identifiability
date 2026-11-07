# NOVA-D — Detectability of Nonlinear Acquisition Advantage: result record v1.2 (7 Nov 2026)

**Correction log (v1.1 → v1.2; documentation only).** The exact v1.0 and v1.1 records and
the v1.1 review note are retained. v1.2 corrects the matrix-row proposal, defines the diagnostic subset explicitly,
adds omitted null condition 13, includes the tied pick count in condition 3, and corrects the optimum-hit summary.
It distinguishes measured gain over estimated G1 from recovery against the population comparator, qualifies
causal and population-limit interpretations, and corrects the threshold and matrix references. The agreed G1
study decisions are recorded in §8.3. All changes use existing records read-only; no experiment, calibration,
bootstrap, or frozen analysis was rerun. Frozen outputs, code, and preregistrations are unchanged.

**Earlier correction history.** v1.1 corrected the outcome box's n = 128 claim to include condition 59 at n = 512,
clarified 500 decisions as 100 replicates × 5 outer folds, and introduced the Objective A/B distinction. Its claim
that the correction was complete throughout and established the mechanism in every positive cell was too strong;
the remaining issues documented in `NOVA_D_RESULT_v1.1_REVIEW.md` are addressed here. Within-replicate training
sets overlap; the 500 fold-decisions are not 500 independent experiments.

**Outcome: both prespecified claims fail. H_D (the near-tie rule) is refuted on its own grid. At the primary
condition and at the positive control the frozen NOVA-S made zero departures from G1 in 500 block-folds each
(κ = 0 exactly; 100 replicates × 5 outer folds); the ungated NOVA-V recovered κ = 0.008 and κ = −0.007
respectively. No near-tie threshold had a qualifying improvement on unscreened NOVA-S in any of the 72 conditions, and tight thresholds
hurt in 18 condition–γ pairs. The validity controls are clean for NOVA-S. Eight of the nine nonlinear conditions
in which NOVA-S has a positive lower bound are at n = 128 (condition 59, n = 512, is the exception). In the
thirteen-cell diagnostic subset defined in §7.3, G1's oracle regret exceeds NOVA-V's. That comparison does not
by itself identify a shrinkage mechanism or exclude a nonlinear contribution. Null examples show that an
advantage over estimated G1 need not represent nonlinear recovery. Large-sample post-hoc diagnostics suggest
shrinkage-induced ranking distortion (§7.4); finite-sample superiority of an alternative remains untested.**

Author of record: Anas Enoch. Compiled by Claude (Fable 5.1). Protocol `NOVA_D_PROTOCOL_v0.3.md` frozen
2026-11-05T21:11:39Z (commit `a98890d`); run 5–7 Nov 2026 on the owner's machine, four partitions, ~36 h wall.
v1.2 editorial correction: Codex, 7 Nov 2026, based on the preserved records and read-only checks.

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
policy. The frozen experiment was not changed after the freeze; these versioned records correct its reporting.
The `--partial` path of `analyze` was never invoked.

## 2. Prespecified decisions (transcribed from `nova_d_result.json["decisions"]`; maximum rounded)

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
  show κ > 1 — 27 (1.00137), 29 (2.09), 54 (1.07), 56 (1.47) — **all at n = 128**. These are measured
  against estimated G1. The population expectation bound against the population Gaussian comparator is
  E[ρ_policy] ≤ ρ_pop; it does not bound every finite-sample estimate or use estimated G1 as its comparator (§7.3).
- At n ≥ 512 in M2 and M4: NOVA-S κ ≤ 0.030 everywhere; NOVA-V κ has median ≈ 0.01 and is negative in half the
  conditions.

Eighteen null conditions (M1, M3; H ≡ 0):

- **NOVA-S: clean.** Maximum ρ = 0.0071 (condition 11, n = 128); no condition approaches 0.10. The validity
  control passes.
- **NOVA-V in nulls** reaches ρ = 0.118 with LCB₉₅ > 0.10 in condition 3 (M1, n = 512, A = 0.05) and 0.058 in
  conditions 11 and 12. Condition 13 also has a positive lower bound (ρ = 0.00421189,
  LCB₉₅ = 0.00104082). Twelve nulls have positive NOVA-V point estimates, but only 3, 11, 12, 13 have positive
  lower bounds. An advantage over estimated G1 can therefore occur without nonlinear structure; a positive
  contrast alone does not establish nonlinear recovery in M2/M4. (The validity rule in protocol §10 was written for NOVA-S, the frozen primary; this reading of
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
2. **Recoverability.** Both prespecified NOVA-S tests fail with κ = 0. Across nonlinear conditions at
   n ≥ 512, measured κ_S against estimated G1 is at most 0.030, not identically zero: condition 59 has
   κ_S = 0.0298954. Its picks support correction of G1's Gaussian rank-1/rank-2 confusion (§7.3), rather than
   selection of the nonlinear optimum. NOVA-V's median measured κ at n ≥ 512 is approximately 0.01 and it
   is worse than estimated G1 in half those conditions. These ratios are not isolated estimates of recovery
   against the population Gaussian comparator. At the positive control — a 40 % designed population advantage,
   34 % value excess on one candidate, n = 1024 — the maximum paired statistic was 1.03 at reported precision.
3. **Selection.** Near-tie gating did not increase useful departures anywhere and reduced them where any existed.

## 6. What this does and does not establish

- **Established, within this family and grid:** the frozen NOVA-S fails both prespecified recovery tests in
  a favourable (H ≥ 0) family; this does not assert zero measured gain in every cell. The near-tie policy has
  no qualifying improvement and can hurt. Per protocol §4, failure in a favourable family is strong evidence against the tested strategy *within this family
  and grid*; it is not a universal refutation of nonlinear acquisition (owner's qualification, recorded).
- **Not established:** behaviour at n = 2048 or beyond (owner's recorded caveat: the positive control's failure at
  n = 1024 does not establish failure at n = 2048); behaviour when the signed correction H is negative;
  behaviour of estimators other than the frozen ones; anything about biological cohorts.
- **Consistent with every prior programme result:** archive NOVA-S +0.112 % below the floor (use #9); Gate 3A no
  family qualified; N2 "no regime identified"; NHANES-FA 0/25 departures. NOVA-D adds the mechanism: the gate's
  statistic stays below τ = 2 in the designated primary and positive-control conditions despite their designed
  population advantages (maxima 1.19 and 1.03 at reported precision).

## 7. Post-hoc diagnostics (not prespecified; labelled as such)

### 7.1 The gate statistic
See §3. In every disagreeing fold of the primary and positive control, t < 1.2. Read-only inspection of
the archived fold statistics gives maxima 1.1923283718 (primary) and 1.0331181743 (positive control).
Any τ above the larger archived maximum, for example 1.2, would give zero departures in both cells.
The lower hypothetical thresholds were not evaluated here. Their possible effect would not resolve
NOVA-V's poorer optimum identification and higher oracle regret in these two conditions (§3).

### 7.2 Identification of the optimum
G1 identified j\* more often than NOVA-V at n ≥ 512 in both designated conditions (§3). The nonlinear estimator's
disadvantage is in locating the candidate, not only in certifying it.

### 7.3 Defined diagnostic subset: decision regret relative to estimated G1
ρ uses estimated G1's realised loss as its comparator. Thus its measured κ can exceed 1 without violating
the population expectation bound against the population Gaussian comparator. The table contains the nine
nonlinear conditions with LCB₉₅(ρ_NOVA-S) > 0 and all four null conditions with LCB₉₅(ρ_NOVA-V) > 0,
plus the primary and positive control as references. It is not an inventory of every positive point estimate
or every positive NOVA-V lower bound in nonlinear conditions. The latter occur in 21 cells:
22, 24, 25, 27, 29, 32, 36, 39, 40, 47, 50, 54, 56, 57, 58, 59, 60, 62, 64, 68, 69.
The twelve nulls with positive NOVA-V point estimates are 0, 2, 3, 5, 6, 9, 11, 12, 13, 14, 15, 17.

These entries are read-only summaries of `nova_d_result.json` and `replicates.jsonl`; no analysis or
experiment was rerun. “#G1” counts distinct G1 picks across 500 dependent fold-decisions; “top G1 picks”
lists the three leading candidates and any ties at third place (or all candidates when fewer than three).
Regrets are normalised by R_∅. Hit rates and κ values are rounded to three decimals.

| idx | mech | ρ_pop | n | A | κ_S | κ_V | G1→j\* | V→j\* | #G1 | top G1 picks | regret G1 | regret V | ratio |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 3 | M1 | 0.00 | 512 | 0.05 | — | — | 0.172 | 0.540 | 17 | 3:103, 10:96, 0:86, 5:86 | 0.0099 | 0.0054 | 1.85 |
| 11 | M3 | 0.00 | 128 | 0.20 | — | — | 0.234 | 0.512 | 3 | 1:323, 0:117, 2:60 | 0.0114 | 0.0042 | 2.73 |
| 12 | M3 | 0.00 | 512 | 0.05 | — | — | 0.072 | 0.288 | 12 | 3:202, 5:67, 4:59 | 0.0051 | 0.0027 | 1.89 |
| 13 | M3 | 0.00 | 512 | 0.15 | — | — | 0.324 | 0.548 | 2 | 1:338, 0:162 | 0.0020 | 0.0014 | 1.42 |
| 20 | M2 | 0.05 | 128 | 0.20 | 0.038 | 0.221 | 0.286 | 0.474 | 3 | 0:342, 1:143, 2:15 | 0.0080 | 0.0051 | 1.56 |
| 27 | M2 | 0.10 | 128 | 0.05 | 0.187 | 1.001 | 0.018 | 0.104 | 17 | 7:108, 6:97, 14:56 | 0.0153 | 0.0116 | 1.32 |
| 29 | M2 | 0.10 | 128 | 0.20 | 0.478 | 2.090 | 0.148 | 0.482 | 3 | 0:232, 2:194, 1:74 | 0.0433 | 0.0098 | 4.41 |
| 36 | M2 | 0.40 | 128 | 0.05 | 0.013 | 0.217 | 0.100 | 0.162 | 14 | 7:114, 9:82, 5:65 | 0.0256 | 0.0219 | 1.17 |
| 45 | M4 | 0.05 | 128 | 0.05 | 0.194 | 0.603 | 0.028 | 0.076 | 17 | 0:166, 2:101, 16:57 | 0.0112 | 0.0097 | 1.16 |
| 54 | M4 | 0.10 | 128 | 0.05 | 0.182 | 1.069 | 0.000 | 0.016 | 18 | 14:125, 4:74, 0:70 | 0.0165 | 0.0139 | 1.18 |
| 56 | M4 | 0.10 | 128 | 0.20 | 0.329 | 1.471 | 0.000 | 0.000 | 4 | 1:340, 0:154, 3:4 | 0.0449 | 0.0254 | 1.77 |
| 59 | M4 | 0.10 | 512 | 0.20 | 0.030 | 0.863 | 0.000 | 0.000 | 2 | 0:301, 1:199 | 0.0319 | 0.0178 | 1.80 |
| 64 | M4 | 0.40 | 128 | 0.15 | 0.119 | 0.439 | 0.002 | 0.002 | 6 | 1:306, 0:150, 2:17 | 0.0765 | 0.0594 | 1.29 |
| *ref* 31 | M2 | 0.10 | 512 | 0.15 | 0.000 | 0.008 | 0.596 | 0.454 | 2 | 1:298, 0:202 | 0.0052 | 0.0070 | 0.74 |
| *ref* 43 | M2 | 0.40 | 1024 | 0.15 | 0.000 | -0.007 | 0.674 | 0.520 | 2 | 1:337, 0:163 | 0.0166 | 0.0245 | 0.68 |

In all thirteen diagnostic cells G1's oracle regret exceeds NOVA-V's. NOVA-V identifies the optimum in at
most 54.8 % of folds (condition 13), and in at most 10 % in five cells: 45, 54, 56, 59, 64.
Condition 27 is 10.4 %, so it does not belong in the latter count. Among the displayed n = 128 cells,
27, 36, 45, 54 have 14–18 distinct G1 picks; the others have 3–6. The four κ_V > 1 conditions
(27, 29, 54, 56) are all at n = 128.

Condition 59 (n = 512) supports a narrower interpretation: M4 has a 20 % Gaussian gap, G1 chose the
Gaussian runner-up in 199/500 folds, and neither G1 nor NOVA-V chose the rank-7 nonlinear optimum.
NOVA-S's κ = 0.0298954 accompanies a reduction of oracle regret from 0.03192 to 0.03144, consistent
with correcting rank-1/rank-2 confusion. In the two reference rows, the regret ordering reverses.

Lower oracle regret documents a decision advantage relative to estimated G1 in this subset. It does not
alone isolate covariance shrinkage as the cause or exclude nonlinear contributions in the nonlinear cells.
The null examples show that an advantage over estimated G1 need not represent nonlinear recovery;
a matched prospective comparison is required to attribute finite-sample differences to regularization.

### 7.4 Shrinkage mechanism: structural concern and reported large-sample diagnostic
The frozen rule uses α = min(1, max(0.05, d/(d+n))) and shrinks the whole standardised covariance toward
(tr S/d)·I. Joint shrinkage changes both the cross-covariances and the observed-block inverse used in the
conditional score. This can alter candidate rankings. The reported v1.0/v1.1 post-hoc comparison used
200,000-row samples at α = 0.05 versus the same rule with α = 0. The values below are retained from
those records; this editorial correction does not reproduce that diagnostic or certify exact population rankings.

| Condition | Population G ordering (top 4) | Reported G1 at n = 200k, α = 0.05 | Same diagnostic, α = 0 |
|---|---|---|---|
| 3 — null, A 0.05 | 0, 1, 2, 3 | **0, 3, 5, 10** (leading four candidates seen at n = 512) | 0, 1, 2, 4 |
| 43 — positive control | 0, 1, 2, 3 (gap 0.02) | **1, 0**, 2, 3 — reversed top pair in this diagnostic | 0, 1, 2, 3 |
| 31 — primary | 0, 1, 2, 3 (gap 0.02) | 0, 1 with scores 0.7802 vs 0.7789 — gap compressed to 0.17 % | 0, 1 (2 % gap kept) |

With a persistent 0.05 floor, the standardised covariance estimator converges to 0.95 Σ + 0.05 μI
rather than Σ. Gaussian conditional covariance is a nonlinear function of its blocks, so its induced
ranking need not equal the population Gaussian ranking. The implementation also includes a fixed
1e-6 ridge in the observed-block solve, which an exact equivalence or population-limit calculation must retain.
The reported 200k-row diagnostic suggests ranking distortion in these distributions; it is finite-sample
evidence, not an exact proof of their limiting ordering. The α = 0 null row itself differs at rank four
from the stated population ordering, so complete ordering recovery should not be claimed from this table.

**Whether that is a defect depends on the objective.** For *Objective A*, estimating the population
Gaussian ranking, persistent regularization bias can yield an inconsistent decision. For *Objective B*,
minimising finite-sample acquisition loss, bias may be worthwhile for variance reduction. NOVA-D has no
matched comparison of regularization levels within G1 and does not establish a better rule for Objective B.

In the positive control, the reported large-sample shrunk rule selects the nonlinear optimum, suggesting
that the opportunity against the frozen comparator is smaller than the designed ρ_pop = 0.40 against
population Gaussian j_G. Finite-sample G1 selected j* in 67.4 % of folds; NOVA-S made no departures.
This remains a failure of the frozen NOVA-S test. The primary's reported large-sample leading candidate
agrees with population j_G, but exact limiting rankings were not checked here.

The archive's 71 % G1 self-agreement (3,376/4,760 folds; matrix v0.4 §1.5 and §3 item 3) demonstrates
instability. Its attribution to shrinkage, or attribution of all n = 128 scatter to this mechanism, remains
a hypothesis without a matched ablation. The prospective question is **when regularization improves
finite-sample acquisition decisions and when conditional-covariance distortion produces ranking errors.**

**Design oversight, disclosed (Claude).** The calibration solved the DGP against the population Gaussian ordering
G, and the self-test verified the sampler against the DGP's closed forms — not the frozen G1's population limit
against that ordering. The protocol's identity V_{j\*} − V_{j_G} = (H_{j\*} − H_{j_G}) − Γ uses j_G = argmax G;
the implemented comparator's limit can be a different candidate. This was correctable before the freeze and was
not caught. It does not alter the primary decision, the null controls, or the near-tie conclusion; it alters how
much the positive control tested.

## 8. Consequences for the programme

1. **H_D is refuted and the near-tie route is closed.** No further γ, no second grid.
2. **The external-validation hypothesis has adverse synthetic evidence, not a demonstrated external failure.**
   The designated synthetic tests show no NOVA-S departures, and NOVA-V ranks the optimum less often than G1.
   Merely lowering τ is therefore not supported as a remedy. The existing evidence no longer offers a strong
   scientific justification for prioritising that validation; it does not make external failure certain.
   `EXTERNAL_VALIDATION_PREREG_STUB_v0.1.md` and the 15 December LIPID/COMETS deadline remain unchanged.
   Choosing not to proceed is a programme decision, not a modification of the preregistration.
3. **The G1 question requires a new prospective study.** Population ranking consistency (Objective A) and
   finite-sample acquisition decision loss (Objective B) must remain distinct. The agreed primary approach is
   **decoupled conditional regularization**, with acquisition decision loss as the primary outcome. Superiority
   is a hypothesis to test, not a claim that conditioning first dominates joint shrinkage at every n.
   The scope is **Gaussian-only**. Retain original G1 as a comparator; vanishing/zero-floor and data-adaptive
   shrinkage (including Ledoit–Wolf) provide control or secondary arms. Grid axes include **dimension ratio**
   and **observed-block covariance condition number**; define the dimension numerator, effective training
   denominator, and population versus estimated condition number explicitly before freezing.

   **Before freezing**, derive and numerically test exact joint-versus-matched-conditional shrinkage equivalence,
   retaining scale factors and ridge terms. A conditional arm algebraically equivalent to joint shrinkage is an
   equivalence control. The primary decoupled arm must specify which components are regularized independently;
   matching only a scalar shrinkage coefficient is insufficient. No equivalence experiment is performed here.
   These are study decisions, not an amended or newly frozen protocol. The archive stays closed for selector
   development; use a new frozen simulation protocol or a genuinely independent dataset. This is independent
   of NOVA-S. The general distinction between covariance error and decision regret is known; any contribution
   must establish a specific failure regime and, if one exists, an effective alternative.
4. **Corrected matrix-row proposal (for a future v0.5 §1.6; v0.4 is unchanged):**
   “NOVA-D (synthetic, 72 conditions × 100 replicates, frozen code): primary and positive-control κ_S = 0,
   0/500 departures each (100 replicates × 5 dependent outer folds); positive-control max t = 1.03 at
   n = 1024 with designed population advantage 40 %. Near-tie screening has no qualifying improvement and
   hurts in 18 condition–threshold pairs; NOVA-S passes all 18 null controls. Nine nonlinear conditions have
   a positive NOVA-S lower bound: eight at n = 128 and condition 59 at n = 512. Four κ_V > 1 estimates
   occur at n = 128 against estimated G1; the population bound uses a different comparator. In the defined
   thirteen-cell diagnostic subset (nine NOVA-S-positive nonlinear cells plus four NOVA-V-positive nulls,
   positivity defined by lower bound > 0), G1 has greater oracle regret than NOVA-V. This does not identify
   the cause or exclude nonlinear contributions. Reported large-sample post-hoc diagnostics suggest
   shrinkage-induced ranking distortion; exact regularized population ordering and alternative-rule superiority remain
   unestablished by this record.”
5. **NOVA manuscript v0.5:** any future update should state the observed ranking and gate failures with
   their scope limits; the reported G1 shrinkage finding belongs in a separate subsection labelled post-hoc,
   with the distinction between large-sample diagnostic evidence and exact population results explicit.

## 9. Versioned correction and preserved files

This correction adds only `NOVA_D_RESULT_v1.2.md` and its SHA-256 sidecar
`NOVA_D_RESULT_v1.2.sha256`. `NOVA_D_RESULT_v1.0.md`, `NOVA_D_RESULT_v1.1.md`, and
`NOVA_D_RESULT_v1.1_REVIEW.md` remain unchanged. The review note describes v1.1 historically;
v1.2 supersedes its proposed wording, including the diagnostic subset's count and maximum hit rate.

The frozen output JSON, integrity and self-test records, run logs, replicates, code, matrix versions,
and preregistrations are already preserved and are not rewritten, compressed, or recommitted as changes.
No experiment, calibration, bootstrap, or frozen analysis was rerun to prepare this correction.
