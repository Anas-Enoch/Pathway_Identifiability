# NHANES serum fatty acids: prospective plug-in prediction of G1 refit stability on the 2013–2014 cycle — protocol v0.2 (draft, not frozen)

Programme: Milestone M (consolidated acquisition paper). Date: 9 November 2026. Owner: Anas Enoch. Drafted by Claude
(Fable 5.1) from the frozen NHANES-FA records. Revises v0.1 on the owner's review (changes in §9).
Status: **draft for the owner's review. Nothing is hashed, no 2013–2014 file is opened, no prediction is computed.**

## 0. Relation to the closed superiority study, and the authorisation this protocol requires

- Protocol v0.4 (`fcc106d8…`, commit `2bb4990`) tested NOVA-S against G1 and returned **NO-GO** at its development
  gate (`NHANES_FA_RESULT_v1.0.md`, corrected in wording by v1.1: 0/25 NOVA-S departures; NOVA-V and NOVA-U 24/25;
  FAS_H/DEMO_H never downloaded). That NO-GO stands and governs the superiority question. **This protocol does not
  reopen it**, computes no loss contrast, and makes no claim about NOVA-S.
- It defines a **different endpoint**: whether a plug-in model fitted to NHANES 2011–2012 predicts the refit
  stability of the frozen, shrunk G1 implementation on the untouched 2013–2014 cycle. It requires its own
  authorisation: **the owner authorises one opening of `FAS_H.xpt` and `DEMO_H.xpt` for this endpoint, after this
  protocol, its implementation and its prediction file are hashed and committed.** The frozen v0.4 external command
  is not used: it refuses to run without a recorded GO and would apply the superseded superiority decision if run.
- Design provenance: blocks, analyte panel, eligibility, fold seeds and the G1 implementation are inherited
  unchanged from v0.4, itself amended after the 2011–2012 development result (v0.4 §13). The stability endpoint was
  chosen after that result showed unanimous within-block G1 picks. The prediction is prospective with respect to
  2013–2014, which no programme file has read.

## 1. Question, estimand and primary outcome

**Question.** Does a plug-in model fitted to NHANES 2011–2012 predict, within a declared interval, how often the
frozen G1 rule changes its chosen candidate when refit across participant-disjoint training folds of NHANES
2013–2014?

**Primary outcome (realised, 2013–2014).** For block b ∈ {1,…,5} and outer fold k ∈ {0,…,4}, let ĵ_{b,k} be the pick of
the frozen `p5_gate3a.g1_pick` trained on all eligible 2013–2014 participants not in fold k. Within-block pairwise
disagreement and its equal-block average:

    D_b = (1/10) Σ_{k<l} 1[ ĵ_{b,k} ≠ ĵ_{b,l} ],      D̄ = (1/5) Σ_b D_b ∈ {0, 1/50, …, 1}.

The five training sets of a block overlap in 60 % of rows and the five blocks share participants and analytes; D̄
is one statistic of one cycle. Development value (2011–2012, `FAS_G_development_v04_result.json`): D̄ = 0.

**What G1 is here.** `g1_pick` standardises its columns on the training rows, forms the sample covariance, applies
the frozen joint shrinkage α = min(1, max(0.05, d/(d+n))) with d = 29 and the ridge 1e-6, and maximises the
conditional score. With n_train in the range relevant here the floor α = 0.05 binds. The pick depends on the training
rows only; the v0.4 target scaling cannot change it (verified by test, §6). NOVA arms do not enter the primary.

## 2. Prediction models, fitted to 2011–2012 only

Fitting data: the 1,725 eligible 2011–2012 participants (v0.4 §2 rule, applied by the frozen `load_cycle`) and the 29
retained analytes (`nhanes_fa_analytes.json`), raw released units. No 2013–2014 value enters any fit.

- **M_G (primary): multivariate Gaussian** with the 2011–2012 sample mean and covariance (ddof 1). It predicts the
  behaviour of the frozen, shrunk G1 implementation under Gaussian sampling of the development covariance; Bayes
  optimality of the Gaussian score belongs to the population oracle under its loss and is not what is tested.
- **M_E (secondary): empirical**, participants resampled with replacement from the 1,725 rows; reported beside M_G,
  never in its place. A discrepancy between the two predictions is itself reported.

**Plug-in prediction by whole-cycle simulation.** For a size n: draw n rows from the model; assign identifiers 1…n;
assign outer folds with the frozen `folds_for(ids, 20261104)` (the 2013–2014 seed; the assignment depends on
identifiers only through a permutation); for each of the five frozen blocks and each fold compute `g1_pick` on the
training rows; compute D̄ and the per-block quantities of §2.1–2.3. B = 2,000 cycles per model and size. The
resulting interval is a **plug-in prediction interval**: it carries the sampling variability of a cycle under the
fitted model, including the shared rows and analytes, and not the uncertainty of the fit itself.

**Sizes.** The grid N = {1,200, 1,250, …, 2,400} ∪ {1,725} is a computational device. The **verdict range** is
[1,200, 2,400]: a realised n_eligible outside it is declared **not assessable** and no verdict is issued. For a
realised n inside the range, the **applicable interval is the union of the 95 % intervals at the two grid points
bracketing n** (a single point if n is on the grid); the union is used so that rounding within a bin cannot
manufacture or remove a verdict. The development check uses the actual development size, n = 1,725.

**Predicted quantities (per model, per grid size), written to `nhanes_fa_prediction_v0.2.json` and hashed:**
1. the simulated distribution of D̄ over its 51 attainable values; the two-sided 95 % interval [q_{.025}, q_{.975}];
   P(D̄ = 0);
2. per block: the distribution of G1's pick over the nine candidates, with mode and modal probability;
3. per block, the **signed score difference on a pair fixed on development**: let j★_b be G1's pick and j°_b its
   runner-up by score on the full 2011–2012 eligible data; record in every simulated fold
   Δs_{b,k} = s_b(j★_b) − s_b(j°_b) (scores as computed inside `g1_pick`); report its simulated distribution and
   P(Δs < 0), so a reversal of the development ordering remains visible as a sign change;
4. for M_G only, at the grid size 1,725: the predicted G1–NOVA-S agreement rate from B_S = 100 full runs of the
   frozen per-block routine (`p6_n2.run_block`, all arms, frozen τ = 2). Because it is computed at one size only, its
   comparison with 2013–2014 is **descriptive unless the realised n_eligible lies within the bin bracketing 1,725**
   (1,700–1,750); otherwise it is reported without a verdict.

## 3. In-sample screens on the development cycle (before freezing; may stop the protocol)

Computed from 2011–2012 alone, at n = 1,725 for (a) and (c) and across the verdict range for (b); recorded in the
prediction file. They are in-sample screens of the fitted model, not external evidence.
- (a) D̄ = 0 must lie inside M_G's 95 % interval at n = 1,725. If not, M_G fails on its own fitting data; the protocol
  is not frozen and that is recorded.
- (b) **Sharpness requirement.** For every adjacent pair of grid points in the verdict range, the union interval's
  upper limit must satisfy q_{.975} ≤ 0.10 (at most five of fifty refit pairs disagreeing). A wider interval could
  still be contradicted; the cap is a requirement that the prediction be sharp enough to be worth testing. If it
  fails anywhere in the range, the range is narrowed to where it holds and the narrowing is recorded; if it holds
  nowhere, the test is recorded as **not sharp enough to run** and FAS_H is not opened.
- (c) the predicted modal pick in each block must equal the 2011–2012 pick (LBXTSA, LBXEPA, LBXDTA, LBXEPA, LBXDTA).
  A mismatch is a model failure on development data; the protocol is not frozen.
Thresholds (0.10; modal agreement) are declared here, before any simulation is run.

## 4. Secondary outcomes (declared; none decides the test)

- S1. Realised per-block picks and D_b, reported individually beside the predicted pick distributions.
- S2. Realised signed score differences Δs_{b,k} on the development-fixed pairs, beside their predicted distributions;
  the count of folds with a sign reversal.
- S3. Realised G1–NOVA-S agreement over the 25 block-folds, from the frozen per-block routine run by the new runner
  under this protocol's authorisation, compared with §2.4 only under the size condition stated there; otherwise
  descriptive.
- No loss contrast (Δ_ext, Δ_G1→∅) is computed or reported under this protocol.

## 5. Interpretation, frozen with the prediction

- Realised D̄ **inside** the applicable interval: the development-fitted plug-in model is compatible with the frozen
  G1's refit behaviour on an independent cycle of the same survey and panel. It does not identify why (a dominant
  candidate or an unresolved near-tie falling the same way — S2 bears on this descriptively) and does not transfer
  to other panels, sizes or selectors.
- Realised D̄ **above** the interval: G1 is less stable on 2013–2014 than M_G predicts. Hypotheses recorded, none
  adjudicated here: tails or dependence absent from M_G (compare M_E first); a population difference between cycles;
  a near-tie resolved differently.
- Realised D̄ **below** the interval (possible only if M_G predicts disagreement): G1 is more stable than predicted;
  same status.
- Not assessable (n outside the verdict range) or not sharp enough (§3b): reported as such; no verdict.
- In every case: no statement about NOVA-S's value, about the information-limitation hypothesis in general, or
  about panels other than serum total fatty acids by GC-ECNI-MS in the NHANES fasting subsample aged 12+.

## 6. Implementation and freeze

- **New runner `nhanes_fa_stability.py`** (written after approval of this draft). It imports, by hash, the frozen
  `nhanes_fa_run.load_cycle` and `arrays` (eligibility and panel), `p5_gate3a.g1_pick` and `folds_for`, and
  `p6_n2.run_block` (for S3 only); it adds no selection logic. Commands: `fit` (M_G and M_E from the local
  2011–2012 table), `predict` (writes the prediction file), `screen` (§3), `freeze` (hashes), `run-ext` (one
  opening of FAS_H/DEMO_H under this protocol's authorisation; refuses if a result file exists, asserts SEQN
  disjointness from the development SEQN hash list, records file hashes), `compare`.
- Self-tests before hashing: (i) `g1_pick` on the five 2011–2012 blocks reproduces the 25 recorded picks exactly,
  with and without the v0.4 target scaling; (ii) `folds_for` on synthetic identifiers gives the frozen fold-size
  pattern; (iii) D̄ ∈ {0, 1/50, …, 1}; (iv) a synthetic Gaussian block with a planted near-tie gives D̄ > 0 at the
  frequency expected from its planted gap (positive control for the statistic); (v) M_E with n = 1,725 and the
  identity resample reproduces D̄ = 0 and the 25 picks; (vi) the signed difference Δs on the development data is
  positive in every block for the fixed pairs.
- Freeze: SHA-256 of this protocol, `nhanes_fa_stability.py`, `nhanes_fa_prediction_v0.2.json`, and the unchanged
  inherited files (`nhanes_fa_run.py` `f85d5235…`, `nhanes_fa_blocks_v0.4.json` `e3e1d8ac…`,
  `nhanes_fa_analytes.json`, `nhanes_fa/frozen/p6_n2.py` `681ce852…`, `nhanes_fa/frozen/p5_gate3a.py`
  `4600e2f3…`) in `NHANES_FA_STABILITY_FREEZE_HASHES.txt`, committed; **then** the §0 authorisation is recorded;
  **then** FAS_H/DEMO_H are downloaded once with their hashes.
- Run: `run-ext` once (seed 20261104, five blocks); `compare`. A crashed run may be repeated with identical settings;
  nothing is re-seeded, re-allocated or re-fitted after FAS_H is read. Per-participant tables and NHANES files stay
  outside the repository.
- Result record `NHANES_FA_STABILITY_RESULT_v1.0.md`: realised n_eligible and bracketing grid sizes, realised D̄ and
  D_b, the applicable interval, the §5 sentence that applies, S1–S3, file hashes, deviations.

## 7. Corrections to adjacent records (owner's finding, to be made before freezing)

In block 3, fold 1 of the 2011–2012 run, NOVA-V and NOVA-U chose `LBXMR1` while G1 and NOVA-S chose `LBXDTA`;
NOVA-V and NOVA-U agree with G1 in 24/25 block-folds, NOVA-S in 25/25; G1's within-block disagreement is 0; the
development contrast and NO-GO are unchanged. `NHANES_FA_RESULT_v1.1.md` (owner's draft) carries the wording
correction with v1.0 and the saved output preserved; the same overstatement is to be corrected at
`NOVA_PAPER_v0.5.md:377` and `SELECTOR_CLAIM_EVIDENCE_MATRIX_v0.4.md:119`.

## 8. Inherited and new

Inherited unchanged: eligibility, 29-analyte panel, five block allocations, fold seeds (20261103 dev, 20261104 ext),
`g1_pick`, `folds_for`, `run_block`, `load_cycle`, `arrays`. New: endpoint (§1), plug-in models and size rule (§2),
in-sample screens with the sharpness requirement (§3), frozen interpretation (§5), the separate runner and its
authorisation (§0, §6).

## 9. Changes from v0.1 (owner's review of 9 November 2026)

1. Separate stability runner; the v0.4 external command is not used; E1 (loss contrasts) dropped entirely.
2. Grid kept as a computational device; verdict range declared; out-of-range cycles not assessable; applicable
   interval is the union over the bracketing grid points; sharpness checked across the range; development check at
   n = 1,725.
3. Interval named a plug-in prediction interval; model statement corrected (M_G predicts the frozen shrunk
   implementation; Bayes optimality is the oracle's); development checks named in-sample screens.
4. Gap diagnostic replaced by the signed score difference on a candidate pair fixed on development.
5. NOVA-S agreement prediction descriptive unless the realised size falls in the bin bracketing 1,725.
6. §7 updated with the owner's verified candidate names and the two adjacent records to correct.
