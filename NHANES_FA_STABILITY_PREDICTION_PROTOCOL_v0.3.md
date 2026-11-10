# NHANES serum fatty acids: prospective plug-in prediction of G1 refit stability on the 2013–2014 cycle — protocol v0.3 (draft, not frozen)

Programme: Milestone M (consolidated acquisition paper). Date: 10 November 2026. Owner: Anas Enoch. Drafted by Claude
(Fable 5.1) from the frozen NHANES-FA records. Revises v0.2 on the owner's review (changes in §9).
Status: **draft for the owner's review. Nothing is hashed, no 2013–2014 file is opened, no prediction is computed.**

## 0. Relation to the closed superiority study, and the authorisation this protocol requires

- Protocol v0.4 (`fcc106d8…`, commit `2bb4990`) tested NOVA-S against G1 and returned **NO-GO** at its development
  gate (`NHANES_FA_RESULT_v1.0.md`, wording corrected in v1.1: NOVA-S 25/25 agreement with G1, NOVA-V and NOVA-U
  24/25; FAS_H/DEMO_H never downloaded). That NO-GO stands and governs the superiority question. **This protocol does
  not reopen it**, computes no held-out loss, and makes no claim about NOVA-S.
- It defines a **different endpoint**: whether a plug-in model fitted to NHANES 2011–2012 predicts the refit stability
  of the frozen, shrunk G1 implementation on the untouched 2013–2014 cycle. It requires its own authorisation: **the
  owner authorises one opening of `FAS_H.xpt` and `DEMO_H.xpt` for this endpoint, after this protocol, its
  implementation and its prediction file are hashed and committed.** The frozen v0.4 external command is not used.
- Design provenance: blocks, panel, eligibility, fold seeds and the G1 implementation are inherited unchanged from
  v0.4, itself amended after the 2011–2012 development result (v0.4 §13). The stability endpoint was chosen after that
  result showed unanimous within-block G1 picks. The prediction is prospective with respect to 2013–2014, which no
  programme file has read.

## 1. Question, estimand and primary outcome

**Question.** Does a plug-in model fitted to NHANES 2011–2012 predict, within a declared interval, how often the
frozen G1 rule changes its chosen candidate when refit across participant-disjoint training folds of NHANES
2013–2014?

**Primary outcome.** For block b ∈ {1,…,5} and outer fold k ∈ {0,…,4}, ĵ_{b,k} is the pick of the frozen
`p5_gate3a.g1_pick` trained on all eligible 2013–2014 participants not in fold k. Per block the **integer
disagreement count** over the ten fold pairs, and the pooled statistic:

    c_b = Σ_{k<l} 1[ ĵ_{b,k} ≠ ĵ_{b,l} ] ∈ {0, 4, 6, 7, 8, 9, 10},      C = Σ_b c_b ∈ {0,…,50},      D̄ = C/50.

The per-block values are the only ones attainable with five picks (partitions of 5: (5)→0, (4,1)→4, (3,2)→6,
(3,1,1)→7, (2,2,1)→8, (2,1,1,1)→9, (1,…,1)→10); C is a sum of five such values, so some multiples of 1/50 are
unattainable. Counting is done on integers; D̄ is a display scale. The five training sets of a block overlap in 60 %
of rows and the blocks share participants and analytes; C is one statistic of one cycle. Development value
(2011–2012): C = 0.

**What G1 is here.** `g1_pick` standardises its columns on the training rows, forms the sample covariance, applies the
frozen joint shrinkage α = min(1, max(0.05, d/(d+n))), d = 29 (the floor 0.05 binds at the relevant n), the ridge
1e-6, and returns the lowest index attaining the maximal conditional score. The pick depends on the training rows
only; the v0.4 target scaling cannot change it (tested, §6). NOVA arms do not enter the primary.

## 2. Prediction models, fitted to 2011–2012 only

Fitting data: the 1,725 eligible 2011–2012 participants (frozen `load_cycle` rule) and the 29 retained analytes, raw
released units. No 2013–2014 value enters any fit.

- **M_G (primary): multivariate Gaussian** with the 2011–2012 sample mean and covariance (ddof 1). It predicts the
  behaviour of the frozen, shrunk G1 implementation under Gaussian sampling of the development covariance.
- **M_E (secondary): empirical**, participants resampled with replacement from the 1,725 rows; reported beside M_G,
  never in its place; a discrepancy between the two is itself reported.

**Plug-in prediction by whole-cycle simulation at a size n.** Draw n rows from the model; identifiers 1…n; outer
folds by the frozen `folds_for(ids, 20261104)`; for each of the five frozen blocks and each fold compute `g1_pick`
on the training rows; record c_b, C, the picks and the signed differences of §2.3. B = 2,000 cycles per model and
size, with the per-size generator seed `sha256("NHANES-FA-STAB-v0.3|" + model + "|" + n)[:15]` so every prediction
is reproducible from the frozen file and code alone. The resulting interval is a **plug-in prediction interval**:
it carries the sampling variability of a cycle under the fitted model, including shared rows and analytes, not the
uncertainty of the fit.

**Exact-size predictions, frozen before opening.** Predictions are computed and frozen **at every integer size in
the verdict range [1,200, 2,400]** for both models (1,201 sizes × 2 models × 2,000 cycles; each cycle is 25
`g1_pick` calls on ≈ n × 29 arrays, so the full table is hours of compute on the owner's machine, run once before
the freeze). The applicable prediction is the one frozen at the **exact realised n_eligible**. No interpolation,
no bracketing, no union. A realised n_eligible outside the verdict range is declared **not assessable**; no verdict
is issued and nothing further is computed for the primary. The coarser 50-participant grid is retained only as the
development screen's summary axis (§3) and for plots.

**Discrete quantile convention (frozen).** For the simulated distribution of C at the applicable n, the 95 %
prediction set is [L, U] with L = min{c : P(C ≤ c) ≥ 0.025} and U = min{c : P(C ≤ c) ≥ 0.975}; by construction
P(L ≤ C ≤ U) ≥ 0.95 under the model. The verdict rule is: realised C ∈ [L, U] → *compatible*; C > U → *above*;
C < L → *below*. The exact simulated tail probabilities P(C ≥ realised) and P(C ≤ realised) are reported as
descriptive quantities beside the verdict.

**Predicted quantities (per model, per size), written to `nhanes_fa_prediction_v0.3.json` and hashed:**
1. the simulated distribution of C (counts over 0…50); [L, U]; P(C = 0);
2. per block: the distribution of G1's pick over the nine candidates, mode and modal probability;
3. per block, the **signed score difference on a pair fixed on the full 2011–2012 fit**: j★_b = G1's pick and j°_b =
   its runner-up by score on all 1,725 eligible development rows. In every simulated fold record
   Δs_{b,k} = s_b(j★_b) − s_b(j°_b), the order of the pair being kept fixed in every refit; report its distribution
   and the per-fold reversal probability P(Δs_{b,k} < 0). Scores are computed by a separate function `g1_scores`
   that reproduces the internals of `g1_pick` and returns the full score vector (§6 tests its argmax against the
   frozen function).
4. **G1–NOVA-S agreement, M_G only, at the single size n = 1,725**, from B_S = 100 simulated cycles through a
   **picks-only path** that reproduces `run_block`'s selection logic (inner cross-fitted scores, V̂, ŝe, the τ = 2
   claim rule) without any held-out loss evaluation; its equivalence to the frozen `run_block` picks is tested on
   synthetic blocks (§6). Because it is simulated at one size only, its comparison with 2013–2014 is **descriptive at
   every realised size**; no verdict attaches to it.

## 3. In-sample screens on the development cycle (before freezing; may stop the protocol)

Computed from 2011–2012 alone; recorded in the prediction file; in-sample screens of the fitted model, not external
evidence.
- (a) C = 0 must lie in M_G's 95 % prediction set at n = 1,725. If not, M_G fails on its own fitting data; the
  protocol is not frozen and that is recorded.
- (b) **Sharpness requirement.** U(M_G, n) ≤ 5 (at most five of fifty refit pairs) at every integer n in the verdict
  range. If it fails on part of the range, the verdict range is narrowed to the sizes where it holds and the
  narrowing is recorded before the freeze; if it holds nowhere, the test is recorded as **not sharp enough to run**
  and FAS_H is not opened. A wider interval could still be contradicted; the cap is a requirement that the prediction
  be worth testing.
- (c) the predicted modal pick in each block at n = 1,725 must equal the 2011–2012 pick (LBXTSA, LBXEPA, LBXDTA,
  LBXEPA, LBXDTA); a mismatch is a model failure on development data; the protocol is not frozen.
- (d) the fixed-pair signed difference on the **full development fit** must be strictly positive in every block
  (this is what fixes the pair; fold-level reversals are an outcome, not a screen).
Thresholds (U ≤ 5; modal agreement) are declared here, before any simulation is run.

## 4. Secondary outcomes (declared; none decides the test)

- S1. Realised per-block picks and c_b, individually, beside the predicted pick distributions.
- S2. Realised Δs_{b,k} on the development-fixed pairs and the number of folds with a sign reversal, beside the
  predicted distributions.
- S3. Realised G1–NOVA-S agreement over the 25 block-folds, computed by the picks-only path under this protocol's
  authorisation; **descriptive at every size**, compared with the n = 1,725 prediction as context only.
- No held-out loss and no loss contrast is computed under this protocol.

## 5. Interpretation, frozen with the prediction

- C **compatible**: the development-fitted plug-in model is compatible with the frozen G1's refit behaviour on an
  independent cycle of the same survey and panel. It does not identify why (dominant candidate versus an unresolved
  near-tie falling the same way — S2 bears on it descriptively) and does not transfer to other panels, sizes or
  selectors.
- C **above**: G1 is less stable on 2013–2014 than M_G predicts. Hypotheses recorded, none adjudicated: tails or
  dependence absent from M_G (compare M_E first); a population difference between cycles; a near-tie resolved
  differently.
- C **below** (possible only if L > 0): G1 is more stable than predicted; same status.
- **Not assessable** (n outside the verdict range) or **not sharp enough** (§3b): reported as such; no verdict.
- In every case: no statement about NOVA-S's value, about the information-limitation hypothesis in general, or about
  panels other than serum total fatty acids by GC-ECNI-MS in the NHANES fasting subsample aged 12+.

## 6. Implementation and freeze

- **Runner `nhanes_fa_stability.py`** (written after approval of this draft). Imports by hash the frozen
  `nhanes_fa_run.load_cycle` and `arrays`, `p5_gate3a.g1_pick` and `folds_for`, and `p6_n2.inner_scores`, `CLAIM_SE`,
  `LAMBDA_U` for the picks-only path; defines `g1_scores` (score vector; same standardisation, shrinkage, ridge and
  clip as `g1_pick`) and `novas_picks_only` (selection logic of `run_block` without loss evaluation). Adds no
  selection logic of its own. Commands: `fit`, `predict` (exact-size table over the verdict range, both models),
  `screen` (§3), `freeze` (hashes), `run-ext` (one opening under §0; refuses if a result exists; asserts SEQN
  disjointness from the development SEQN hash list; records file hashes; computes picks, c_b, C, Δs, S3), `compare`.
- Self-tests before hashing: (i) `g1_pick` on the five development blocks reproduces the 25 recorded picks exactly,
  with and without the v0.4 target scaling; (ii) `argmax(g1_scores) == g1_pick` on the 25 development training sets
  and on synthetic cases including exact ties (lowest index); (iii) `folds_for` on synthetic identifiers gives the
  frozen fold-size pattern; (iv) per-block counts c_b ∈ {0, 4, 6, 7, 8, 9, 10} on exhaustive synthetic pick patterns,
  and C computed on integers; (v) a synthetic Gaussian block with a planted near-tie gives c_b > 0 at the frequency
  expected from the planted gap; (vi) M_E with n = 1,725 and the identity resample reproduces C = 0 and the 25
  picks; (vii) `novas_picks_only` reproduces the frozen `run_block` picks (G1, NOVA-V, NOVA-U, NOVA-S) on the
  synthetic blocks of the v0.4 equivalence test; (viii) §3(d) on the full development fit.
- Freeze: SHA-256 of this protocol, `nhanes_fa_stability.py`, `nhanes_fa_prediction_v0.3.json`, and the unchanged
  inherited files (`nhanes_fa_run.py` `f85d5235…`, `nhanes_fa_blocks_v0.4.json` `e3e1d8ac…`,
  `nhanes_fa_analytes.json`, `nhanes_fa/frozen/p6_n2.py` `681ce852…`, `nhanes_fa/frozen/p5_gate3a.py`
  `4600e2f3…`) in `NHANES_FA_STABILITY_FREEZE_HASHES.txt`, committed; **then** the §0 authorisation is recorded;
  **then** FAS_H/DEMO_H are downloaded once with their hashes.
- Run: `run-ext` once (seed 20261104, five blocks); `compare` reads the frozen prediction at the exact realised n. A
  crashed run may be repeated with identical settings; nothing is re-seeded, re-allocated or re-fitted after FAS_H
  is read. Per-participant tables and NHANES files stay outside the repository.
- Result record `NHANES_FA_STABILITY_RESULT_v1.0.md`: realised n_eligible, realised c_b and C, [L, U] at that n, the
  tail probabilities, the §5 sentence that applies, S1–S3, file hashes, deviations.

## 7. Corrections to adjacent records (owner's finding, to be made before freezing)

Block 3, fold 1 of the 2011–2012 run: NOVA-V and NOVA-U chose `LBXMR1`, G1 and NOVA-S `LBXDTA`; agreement with G1
24/25 for NOVA-V and NOVA-U, 25/25 for NOVA-S; G1's within-block disagreement 0; development contrast and NO-GO
unchanged. `NHANES_FA_RESULT_v1.1.md` (owner's draft) carries the correction with v1.0 preserved; matching
corrections at `NOVA_PAPER_v0.5.md:377` and `SELECTOR_CLAIM_EVIDENCE_MATRIX_v0.4.md:119`.

## 8. Inherited and new

Inherited unchanged: eligibility, 29-analyte panel, five block allocations, fold seeds (20261103 dev, 20261104 ext),
`g1_pick`, `folds_for`, `inner_scores` and the claim constants, `load_cycle`, `arrays`. New: endpoint on integer
counts (§1), plug-in models with exact-size predictions and the discrete quantile convention (§2), in-sample
screens with the sharpness requirement (§3), frozen interpretation (§5), the separate runner with `g1_scores` and
the picks-only NOVA path (§6), the authorisation (§0).

## 9. Changes from v0.2 (owner's review of 9–10 November 2026)

1. Size rule: exact-size predictions frozen at every integer n in the verdict range; the applicable prediction is the
   one at the realised n; the bracketing union is removed (it was a sensitivity device, not a prediction interval at
   the realised size); out-of-range cycles not assessable; the 50-grid kept only as a screen/plot axis.
2. Discrete quantile convention and the verdict rule declared; tail probabilities reported descriptively.
3. Primary counted on integers (c_b ∈ {0,4,6,7,8,9,10}, C ∈ 0…50); D̄ = C/50 a display scale; self-test (iv) fixed
   accordingly.
4. Signed gap: pair fixed on the full 2011–2012 fit, order kept in every refit; positivity screen applies to the
   full fit (§3d); fold-level reversals retained as an outcome; `g1_scores` defined with an argmax-equivalence test
   including ties.
5. S3 made descriptive at every size and routed through a picks-only path with an equivalence test; no held-out loss
   computed anywhere under this protocol.
