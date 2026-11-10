# NHANES serum fatty acids: prospective prediction of G1 refit stability on the 2013–2014 cycle — protocol v0.1 (draft, not frozen)

Programme: Milestone M (consolidated acquisition paper). Date: 9 November 2026. Owner: Anas Enoch. Drafted by Claude
(Fable 5.1) from the frozen NHANES-FA records; every constant below is taken from those records or declared here.
Status: **draft for the owner's review. Nothing is hashed, no NHANES file is opened, no prediction is computed.**

## 0. Relation to the closed NHANES-FA superiority study, and the authorisation this protocol requires

- Protocol v0.4 (`fcc106d8…`, commit `2bb4990`) tested NOVA-S against G1 and returned **NO-GO** at its development gate
  (`NHANES_FA_RESULT_v1.0.md`: 0/25 NOVA-S departures; FAS_H/DEMO_H never downloaded). That NO-GO stands. **This
  protocol does not reopen the superiority question**, applies no success criterion to Δ_ext, and makes no claim about
  NOVA-S's advantage.
- It defines a **different endpoint** — whether a model fitted to the 2011–2012 development cycle predicts the
  refit stability of the frozen G1 rule on the untouched 2013–2014 cycle — and therefore requires its own explicit
  authorisation: **the owner authorises a single opening of `FAS_H.xpt` and `DEMO_H.xpt` for this endpoint, once
  this protocol, its implementation and its prediction file are hashed.** Until that signature exists the files are
  not downloaded.
- Design provenance, stated so a reader cannot mistake this for an untouched design: the five blocks, analyte
  panel, eligibility, folds and seeds are inherited unchanged from v0.4, which was itself amended after the
  2011–2012 development result (v0.4 §13). The stability endpoint was chosen after that result showed unanimous
  within-block G1 picks in 2011–2012. The prediction is nevertheless prospective with respect to 2013–2014, which
  no programme file has read.

## 1. Question, estimand and primary outcome

**Question.** Does a model fitted to NHANES 2011–2012 predict, within a declared interval, how often the frozen G1
rule changes its chosen candidate when refit across participant-disjoint training folds of NHANES 2013–2014?

**Primary outcome (realised, 2013–2014).** For block b ∈ {1,…,5} and outer fold k ∈ {0,…,4}, let ĵ_{b,k} be G1's
pick (`p5_gate3a.g1_pick`, frozen) trained on all eligible 2013–2014 participants not in fold k. Within-block
pairwise disagreement

    D_b = (1/10) Σ_{k<l} 1[ ĵ_{b,k} ≠ ĵ_{b,l} ],      D̄ = (1/5) Σ_b D_b.

D̄ takes values in {0, 1/50, …, 1}. Blocks enter with equal weight; the five training sets of a block overlap in
60 % of rows and the five blocks share participants and analytes, so D̄ is one statistic of one cycle, not an
average of independent replicates. Development value (2011–2012, from `FAS_G_development_v04_result.json`): D̄ = 0
(every block: one candidate in all five folds).

**Why G1 alone decides the primary.** `g1_pick` standardises its own columns and uses only the training rows'
sample covariance with the frozen shrinkage α = min(1, max(0.05, d/(d+n))), d = 29; with n_train ≈ 1,100–1,800 the
floor α = 0.05 binds. The target scaling of v0.4 §4 does not change G1's pick (column standardisation inside
`g1_pick` removes it); the implementation test in §6 verifies this on the development data. The NOVA arms do not
enter the primary outcome.

## 2. Prediction models, fitted to 2011–2012 only

Fitting data: the 1,725 eligible 2011–2012 participants (v0.4 §2 eligibility) and the 29 retained analytes
(`nhanes_fa_analytes.json`), raw released units (µmol/L, CDC fill for below-LLOD as released). No 2013–2014 value
enters any fit.

- **M_G (primary): multivariate Gaussian** with the 2011–2012 sample mean and sample covariance (ddof 1) of the 29
  analytes. This is the model in which G1 is Bayes-optimal and in which the programme's exact results (Wishart
  residual law, pair law) hold; it is the model whose calibration the test examines.
- **M_E (secondary): empirical.** Participants resampled with replacement from the 1,725 rows. It carries the real
  marginal tails and dependence of the development cycle and is reported beside M_G, never in place of it. A
  discrepancy between M_G and M_E predictions is itself reported.

**Simulation of a cycle** of size n under a model: draw n rows; assign synthetic identifiers 1…n; assign outer folds
with the frozen `folds_for(ids, 20261104)` (the 2013–2014 seed; the assignment depends on the identifiers only
through a permutation, so synthetic identifiers are distributionally equivalent to SEQN); for each of the five
frozen blocks and each fold compute `g1_pick` on the training rows; compute D̄. Repeat B = 2,000 times per model and
per n.

**Size grid.** n_eligible for 2013–2014 is unknown until the files are opened and is not an outcome. Predictions are
computed for n ∈ N = {1,200, 1,300, …, 2,400} and stored; the applicable prediction is the grid value nearest to the
realised n_eligible (ties to the lower value). If the realised n_eligible lies outside [1,200, 2,400] the nearest
endpoint applies and the fact is reported.

**Predicted quantities (per model, per n), written to `nhanes_fa_prediction_v0.1.json` and hashed:**
1. the full simulated distribution of D̄ (counts over its 51 attainable values);
2. the two-sided 95 % prediction interval [q_{.025}, q_{.975}] and P(D̄ = 0);
3. per block: the distribution of G1's pick over the nine candidates (mode and its probability) and the
   standardised top-two score gap, (s_(1) − s_(2)) / SD_sim(s_(1) − s_(2)), as a descriptive resolution diagnostic;
4. for M_G only, at the single grid value nearest 1,725: the predicted G1–NOVA-S agreement rate from B_S = 100 full
   runs of the frozen pipeline (`nhanes_fa_run.py` v0.2.0 `mirror_block`, all arms), declared secondary because of
   its smaller B; its 95 % interval is reported with that B stated.

## 3. Informativeness check on the development cycle (before freezing; may stop the protocol)

Computed from 2011–2012 alone and recorded in the prediction file:
- (a) the development value D̄ = 0 must lie inside M_G's 95 % interval at n = 1,700 (the grid value nearest 1,725).
  If it does not, M_G is miscalibrated on its own fitting data; the protocol is not frozen and the finding is
  recorded as such.
- (b) the interval must be able to exclude a plausible external result: the protocol is frozen only if
  q_{.975}(M_G, n = 1,700) ≤ 0.10, i.e. the model must predict that at most five of the fifty refit pairs disagree.
  If the model's own interval is wider than that, a realised 2013–2014 value could not contradict it, and the test
  is recorded as **uninformative by design** rather than run.
- (c) the per-block predicted modes must coincide with the 2011–2012 picks (LBXTSA, LBXEPA, LBXDTA, LBXEPA, LBXDTA);
  a mismatch is recorded as a model failure on development data and the protocol is not frozen.
Both thresholds (0.10; the mode check) are declared here, before any simulation is run.

## 4. Secondary and exploratory outcomes (all declared; none decides the test)

- S1. G1–NOVA-S agreement rate on 2013–2014 (25 block-folds), compared with the M_G prediction of §2.4. The NOVA-S
  gate is the frozen one (τ = 2 claim rule); nothing is retuned.
- S2. Per-block realised picks and the realised standardised top-two gap of G1 scores, beside the predicted ones.
- S3. The realised distribution of D_b across the five blocks, reported individually.
- E1. Δ_ext and Δ_G1→∅ with their design-based intervals, **as computed by the frozen v0.4 pipeline and reported
  descriptively under the standing NO-GO**; no success criterion is applied and no superiority or inferiority
  statement is made. They are reported because the pipeline produces them and withholding computed numbers would
  be worse than labelling them.

## 5. Interpretation, frozen with the prediction

- Realised D̄ **inside** M_G's interval at the applicable n: the development-fitted Gaussian model is compatible with
  G1's refit behaviour on an independent cycle of the same survey and panel. It does not identify why G1 is stable
  (a dominant candidate, or an unresolved near-tie that happened to fall the same way) — S2 bears on that
  descriptively — and it does not transfer to other panels or sample sizes.
- Realised D̄ **above** the interval: G1 is less stable on 2013–2014 than the Gaussian fit predicts. Candidate
  explanations are recorded as hypotheses only (heavier tails or dependence not captured by M_G — M_E's prediction
  is the first thing to compare; a population difference between cycles; a near-tie resolved differently). None is
  adjudicated by this protocol.
- Realised D̄ **below** the interval (possible only if M_G predicts some disagreement): G1 is more stable than the
  model predicts; same status.
- In every case: no statement about NOVA-S's value, about the information-limitation hypothesis in general, or about
  panels other than serum total fatty acids by GC-ECNI-MS in the NHANES fasting subsample aged 12+. The broader
  claim that realistic panels are typically information-limited remains a hypothesis outside this test.

## 6. Implementation and freeze

- New file `nhanes_fa_predict.py` (to be written after the owner approves this draft): fits M_G and M_E from the
  local 2011–2012 participant table, simulates as in §2, writes `nhanes_fa_prediction_v0.1.json`, and after the
  2013–2014 run computes the realised D̄ and the comparison. It imports the frozen `p5_gate3a.g1_pick` and
  `folds_for` by hash and adds no selection logic.
- Self-tests before hashing: (i) `g1_pick` on the 2011–2012 blocks reproduces the 25 recorded picks exactly, with
  and without the v0.4 target scaling (scale invariance); (ii) fold assignment with synthetic identifiers has the
  frozen fold-size pattern; (iii) simulated D̄ is in {0, 1/50, …, 1}; (iv) a synthetic Gaussian block with a
  planted near-tie yields D̄ > 0 with the expected frequency (positive control for the statistic); (v) M_E with
  n = 1,725 and no resampling (identity permutation) reproduces D̄ = 0.
- Freeze: SHA-256 of this protocol, `nhanes_fa_predict.py`, `nhanes_fa_prediction_v0.1.json`, and the unchanged
  inherited files (`nhanes_fa_run.py` `f85d5235…`, `nhanes_fa_blocks_v0.4.json` `e3e1d8ac…`,
  `nhanes_fa_analytes.json`, `nhanes_fa/frozen/p6_n2.py` `681ce852…`, `nhanes_fa/frozen/p5_gate3a.py`
  `4600e2f3…`) recorded in `NHANES_FA_PREDICTION_FREEZE_HASHES.txt`, committed, **then** the owner's authorisation
  of §0 is recorded, **then** FAS_H/DEMO_H are downloaded with their hashes.
- Run: the frozen v0.4 pipeline once on 2013–2014 (seed 20261104, all five blocks, wrapper v0.2.0 unchanged); then
  `nhanes_fa_predict.py compare`. One run; a crashed run may be repeated with identical settings; nothing is
  re-seeded, re-allocated or re-fitted after FAS_H is read. Per-participant tables stay local; the NHANES files stay
  outside the repository.
- Result record `NHANES_FA_STABILITY_RESULT_v1.0.md`: realised n_eligible, applicable grid n, realised D̄ and D_b,
  the declared interval, the §5 sentence that applies, S1–S3, E1 with the NO-GO note, hashes of the downloaded
  files, deviations.

## 7. Correction owed to the development record (independent of this protocol)

`NHANES_FA_RESULT_v1.0.md` §2 heads its table "Candidate chosen by G1, NOVA-V, NOVA-U and NOVA-S in all 5 folds" and
§3.1 says the unrestricted NOVA-V "made identical acquisition decisions in 25 of 25 block-folds". The frozen result
file shows that in **block 3, fold 1, NOVA-V and NOVA-U chose candidate index 0 while G1 and NOVA-S chose index 6**
(`FAS_G_development_v04_result.json`, `inference.picks[10]`); NOVA-S followed G1 because the claim rule was not met.
NOVA-S agreement is therefore 25/25 as stated, NOVA-V and NOVA-U agreement is 24/25. The gate decision, Δ_dev = 0
and the NO-GO are unaffected. A v1.1 wording correction of the result record, and matching corrections wherever the
"all four selectors in 25/25" phrasing was carried (NOVA manuscript v0.5 limitations; claim matrix v0.4 row), should be
made before this protocol is frozen so that the two records agree.

## 8. What is inherited and what is new

Inherited unchanged from v0.4: eligibility, 29-analyte panel, five block allocations, fold seeds (20261103 dev,
20261104 ext), wrapper v0.2.0, frozen `p6_n2.py` and `p5_gate3a.py`, target-scaling convention, design-variance code
(used only for E1). New: the endpoint (§1), the prediction models and size grid (§2), the informativeness gate (§3),
the frozen interpretation (§5), `nhanes_fa_predict.py` and the prediction file (§6), and the authorisation clause (§0).
