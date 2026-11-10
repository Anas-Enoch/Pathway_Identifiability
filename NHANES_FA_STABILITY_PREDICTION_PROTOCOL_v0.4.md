# NHANES serum fatty acids: prospective plug-in prediction of G1 refit stability on the 2013–2014 cycle — protocol v0.4 (draft, not frozen)

Programme: Milestone M (consolidated acquisition paper). Date: 10 November 2026. Owner: Anas Enoch. Drafted by Claude
(Fable 5.1) from the frozen NHANES-FA records. Revises v0.3 on the owner's approval of 10 November 2026 (changes in
§9). Status: **draft; approved for the 2011–2012 development implementation (`fit`, `predict`, `screen`) only.
Nothing is hashed, no 2013–2014 file is opened, `run-ext` is not authorised.**

## 0. Relation to the closed superiority study, and the authorisation this protocol requires

- Protocol v0.4 (`fcc106d8…`, commit `2bb4990`) tested NOVA-S against G1 and returned **NO-GO** at its development
  gate (`NHANES_FA_RESULT_v1.0.md`, wording corrected in v1.1: NOVA-S 25/25 agreement with G1, NOVA-V and NOVA-U
  24/25; FAS_H/DEMO_H never downloaded). That NO-GO stands and governs the superiority question. **This protocol does
  not reopen it**, computes no held-out loss, and makes no claim about NOVA-S.
- It defines a **different endpoint**: whether a plug-in model fitted to NHANES 2011–2012 predicts the refit stability
  of the frozen, shrunk G1 implementation on the untouched 2013–2014 cycle. It requires its own authorisation, in
  the form fixed by §6: **a file `NHANES_FA_STABILITY_AUTHORISATION.txt`, written by the owner after the freeze,
  containing the word AUTHORISED and the frozen SHA-256 of this protocol, of the runner and of the prediction file.**
  The runner refuses to read either external file without it. The frozen v0.4 external command is not used.
- Design provenance: blocks, panel, eligibility, fold seeds and the G1 implementation are inherited unchanged from
  v0.4, itself amended after the 2011–2012 development result (v0.4 §13). The stability endpoint was chosen after that
  result showed unanimous within-block G1 picks. The prediction is prospective with respect to 2013–2014, which no
  programme file has read.

## 1. Question, estimand and primary outcome

**Question.** Does a plug-in model fitted to NHANES 2011–2012 predict, within a declared set, how often the frozen G1
rule changes its chosen candidate when refit across the five outer-fold training sets of NHANES 2013–2014?

**Primary outcome.** For block b ∈ {1,…,5} and outer fold k ∈ {0,…,4}, ĵ_{b,k} is the pick of the frozen
`p5_gate3a.g1_pick` trained on all eligible 2013–2014 participants not in fold k. Per block the **integer
disagreement count** over the ten fold pairs, and the pooled statistic:

    c_b = Σ_{k<l} 1[ ĵ_{b,k} ≠ ĵ_{b,l} ] ∈ {0, 4, 6, 7, 8, 9, 10},      C = Σ_b c_b ∈ {0,…,50},      D̄ = C/50.

The per-block values are the only ones attainable with five picks (partitions of 5: (5)→0, (4,1)→4, (3,2)→6,
(3,1,1)→7, (2,2,1)→8, (2,1,1,1)→9, (1,…,1)→10; verified exhaustively, §6 test iv); C is a sum of five such values, so
some multiples of 1/50 are unattainable. Counting is done on integers; D̄ is a display scale. The held-out folds are
participant-disjoint; the five training sets of a block are not — any two share about 60 % of their rows — and the
blocks share participants and analytes; C is one statistic of one cycle. Development value (2011–2012): C = 0.

**What G1 is here.** `g1_pick` standardises its columns on the training rows, forms the sample covariance, applies the
frozen joint shrinkage α = min(1, max(0.05, d/(d+n))), d = 29 (the floor 0.05 binds at the relevant n), the ridge
1e-6, and returns the lowest index attaining the maximal conditional score. The pick depends on the training rows
only; the v0.4 target scaling cannot change it (tested, §6 i–ii). NOVA arms do not enter the primary.

## 2. Prediction models, fitted to 2011–2012 only

Fitting data: the 1,725 eligible 2011–2012 participants (frozen `load_cycle` rule) and the 29 retained analytes, raw
released units. No 2013–2014 value enters any fit.

- **M_G (primary): multivariate Gaussian** with the 2011–2012 sample mean and covariance (ddof 1). It predicts the
  behaviour of the frozen, shrunk G1 implementation under Gaussian sampling of the development covariance.
- **M_E (secondary): empirical**, participants resampled with replacement from the 1,725 rows; reported beside M_G,
  never in its place; a discrepancy between the two is itself reported.

**Plug-in prediction by whole-cycle simulation at a size n.** Draw n rows from the model; identifiers 0…n−1 as
strings; outer folds by the frozen `folds_for(ids, 20261104)`; for each of the five frozen blocks and each fold
compute the G1 pick on the training rows (through `g1_scores`, whose lowest-index argmax is tested equal to
`g1_pick`, §6 ii); record c_b, C, the picks and the signed differences of §2.3. B = 2,000 cycles per model and size,
with the per-size generator seed `sha256("NHANES-FA-STAB-v0.4|" + model + "|" + n)[:15]`, so every entry of the
table is reproducible from the frozen fit file and code alone.

**What the interval is (calibration stated precisely).** For the simulated distribution of C at a size n, the
reported set is [L, U] with L = min{c : P̂(C ≤ c) ≥ 0.025} and U = min{c : P̂(C ≤ c) ≥ 0.975}, P̂ being the Monte
Carlo proportion over the B cycles. [L, U] is an **estimate** of the model's 95 % plug-in prediction set: under the
fitted model, with the fitted parameters taken as true, the set it estimates contains C with probability ≥ 0.95.
The Monte Carlo error of the endpoints at B = 2,000 is small but not zero (a tail proportion of 0.025 has standard
error ≈ 0.0035), and the set carries **no allowance for the uncertainty of the fit**: it describes the sampling
variability of a cycle under M_G (including the shared rows and analytes), not the probability that M_G is right.
The tail quantities P̂(C ≥ realised) and P̂(C ≤ realised) reported beside the verdict are likewise Monte Carlo
proportions under the fitted model, not p-values of any hypothesis test.

**Exact-size predictions, frozen before opening.** Predictions are computed at every integer size in the verdict
range [1,200, 2,400] for both models (1,201 sizes × 2 models × 2,000 cycles; 25 G1 refits per cycle). The applicable
prediction is the one frozen at the **exact realised n_eligible**. No interpolation, no bracketing, no union. A
realised n_eligible outside the verdict range is declared **not assessable**; no verdict is issued and nothing
further is computed for the primary.

**Workload, benchmarked before launch.** `bench` times ten cycles per model at n = 1,725 and extrapolates to the
full table; it is run and its output recorded before `predict` is launched. Reference measurement (10 November 2026,
Linux container, one core, numpy 2.5.3, synthetic 29-column data of the same shape): 18 ms per M_G cycle, 24 ms per
M_E cycle, i.e. ≈ 0.6 and 0.8 core-minutes per (model, size) and **≈ 28 core-hours for the complete table**; with k
worker processes the wall time divides by about k. `predict` is resumable: each (model, size) is written as its own
part file, completed parts are never recomputed, and `--sizes a:b` restricts a launch to a sub-range. **A partial
table has no inferential status**: `assemble` marks it incomplete, `screen` refuses it, nothing is frozen, and no
external opening can follow from it. If the complete table proves impractical on the owner's hardware, the
protocol is amended prospectively (a narrower verdict range or smaller B, declared and justified before any part
of the realised 2013–2014 size is known) or the test is recorded as not run; no partial table is reinterpreted.

**Verdict rule (frozen).** Realised C ∈ [L, U] at the realised n → *compatible*; C > U → *above*; C < L → *below*.

**Predicted quantities (per model, per size), written to `nhanes_fa_prediction_v0.4.json` and hashed:**
1. the simulated distribution of C (counts over 0…50); [L, U]; P̂(C = 0);
2. per block: the distribution of G1's pick over the nine candidates, mode and modal probability;
3. per block, the **signed score difference on a pair fixed on the full 2011–2012 fit**: j★_b = G1's pick and j°_b =
   its runner-up by score on all 1,725 eligible development rows. In every simulated fold record
   Δs_{b,k} = s_b(j★_b) − s_b(j°_b), the order of the pair being kept fixed in every refit; report its mean, standard
   deviation and the per-fold reversal probability P̂(Δs_{b,k} < 0); at n = 1,725 the individual values are kept
   and their hash recorded.
4. **G1–NOVA-S agreement, M_G only, at the single size n = 1,725**, from B_S = 100 simulated cycles through the
   **picks-only path** `novas_picks_only`, which reproduces `run_block`'s selection logic (inner cross-fitted
   scores, V̂, ŝe, the τ = 2 claim rule) without any held-out loss evaluation; its equivalence to the frozen
   `run_block` picks is tested on synthetic blocks (§6 vii). Because it is simulated at one size only, its
   comparison with 2013–2014 is **descriptive at every realised size**; no verdict attaches to it.

## 3. In-sample screens on the development cycle (before freezing; may stop the protocol)

Computed from 2011–2012 alone after the complete table exists; recorded in `screen_v0.4.json`; in-sample screens
of the fitted model, not external evidence.
- (a) C = 0 must lie in M_G's set [L, U] at n = 1,725. If not, M_G fails on its own fitting data; the protocol is not
  frozen and that is recorded.
- (b) **Sharpness requirement and the assessable-size set.** At a size n the prediction is *sharp* if U(M_G, n) ≤ 5
  (at most five of fifty refit pairs; with the attainable values this means U ∈ {0, 4}). `screen` evaluates this at
  every integer n of the verdict range and **stores the set of passing integers explicitly** as
  `assessable_sizes` in `screen_v0.4.json`, which is hashed in the freeze. A realised n_eligible inside the verdict
  range but **not in the stored set is declared not assessable**: no verdict is issued. If the set is empty the
  test is recorded as **not sharp enough to run** and FAS_H is not opened. A wider interval could still be
  contradicted; the cap is a requirement that the prediction be worth testing.
- (c) the predicted modal pick in each block at n = 1,725 must equal the 2011–2012 full-fit pick (LBXTSA, LBXEPA,
  LBXDTA, LBXEPA, LBXDTA); a mismatch is a model failure on development data; the protocol is not frozen.
- (d) the fixed-pair signed difference on the **full development fit** must be strictly positive in every block
  (this is what fixes the pair; fold-level reversals are an outcome, not a screen).
Thresholds (U ≤ 5; modal agreement) are declared here, before any simulation is run. `freeze` refuses unless
`screen` recorded (a), (c), (d) passed and the assessable set is non-empty.

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
- **Not assessable** (n outside the verdict range, or inside it but outside the stored assessable set) or **not
  sharp enough** (§3b, empty set): reported as such; no verdict.
- In every case: no statement about NOVA-S's value, about the information-limitation hypothesis in general, or about
  panels other than serum total fatty acids by GC-ECNI-MS in the NHANES fasting subsample aged 12+.

## 6. Implementation, freeze and the external guard

**Runner `nhanes_fa_stability.py` v0.1.0** (written 10 November 2026; synthetic self-tests pass in the container
against the clone at `340ba48`). It verifies by SHA-256 and imports the frozen `nhanes_fa_run.py` (`f85d5235…`),
which loads `nhanes_fa/frozen/p6_n2.py` (`681ce852…`) and `p5_gate3a.py` (`4600e2f3…`), and verifies
`nhanes_fa_blocks_v0.4.json` (`e3e1d8ac…`). From them it uses `load_cycle`, `arrays`, `SEED_DEV`, `g1_pick`,
`folds_for`, `inner_scores`, `K_OUT`, `SEED`, `CLAIM_SE`, `LAMBDA_U`, `FLOOR`, `synth`, `run_block`. It defines
`g1_scores` (score vector; same standardisation, shrinkage, ridge and clip as `g1_pick`), `novas_picks_only`
(selection logic of `run_block` without loss evaluation) and `frozen_cycle` (picks on a loaded NHANES frame exactly
as the frozen wrapper makes them: `arrays` scaling, `folds_for` on SEQN identifiers, `g1_pick(O, Y, X)` per fold, with
the `g1_scores` argmax required equal on every training set). It adds no selection logic of its own.

Commands: `selftest [--dev]`; `fit` (2011–2012: eligible rows, mean, covariance, full-fit picks and runner-up pairs,
reproduction of the 25 recorded development picks through `frozen_cycle` — refuses if they differ); `bench`;
`predict [--sizes a:b] [--workers k]` (resumable parts); `assemble` (parts → prediction file; refuses a part with
B ≠ 2,000; marks an incomplete table); `screen` (§3; stores the assessable set); `freeze` (hash record; never
overwritten); `run-ext` (§6 guard below); `compare`. `fit` writes `nhanes_fa_stab_out/dev_fit.json` (aggregate:
mean, covariance, picks, pairs, scores — committable) and `dev_rows.npy`, `dev_ids.json` (participant level — stay
local, never committed).

**Self-tests.** Synthetic, run in the container and to be re-run on the owner's machine before hashing:
(ii) `argmax(g1_scores) == g1_pick` on 30 random cases, on an exact duplicate-candidate tie (lowest index by both)
and under target rescaling; (iii) `folds_for` on 1,725 synthetic identifiers gives five folds of 345; (iv) c_b over
all 9⁵ pick patterns ∈ {0, 4, 6, 7, 8, 9, 10}; (v) `simulate_size` runs with a histogram summing to B, and a planted
population tie between two candidates gives C > 0 in 97 % of cycles; (vii) `novas_picks_only` reproduces the frozen
`run_block` picks on the v0.4 synthetic blocks (null, quadratic); quantile convention on a test histogram.
Development-data tests (`selftest --dev`, after `fit`, owner's machine): (i) the simulation path (raw rows,
standardised inside `g1_scores`) and the frozen path (v0.4-scaled targets) pick identically on the 25 development
training sets, and `fit` reproduced the recorded 25 picks with C = 0; (vi) the simulation path on the identity
resample with SEQN identifiers and the development seed reproduces the 25 picks and C = 0; (viii) §3(d).

**Freeze.** `NHANES_FA_STABILITY_FREEZE_HASHES.txt`: SHA-256 of this protocol, `nhanes_fa_stability.py`,
`nhanes_fa_prediction_v0.4.json`, `nhanes_fa_stab_out/dev_fit.json`, `nhanes_fa_stab_out/screen_v0.4.json` and the
unchanged inherited files (`nhanes_fa_run.py`, `nhanes_fa_blocks_v0.4.json`, `nhanes_fa_analytes.json`,
`nhanes_fa/frozen/p6_n2.py`, `nhanes_fa/frozen/p5_gate3a.py`); committed. Every later command re-verifies every
listed hash and refuses on any change.

**External guard (applies before the first read).** `run-ext`, in this order: (1) verify all frozen hashes; (2)
require `NHANES_FA_STABILITY_AUTHORISATION.txt` containing AUTHORISED and the frozen hashes of the protocol, the
runner and the prediction file; (3) refuse if any attempt is recorded as completed or a result file exists; (4)
refuse if a previous attempt exists with different inputs or settings (file paths, seed 20261104, protocol and
runner hashes); (5) **append the attempt to `nhanes_fa_stab_out/FAS_H_opening_attempt.json` with status `started`,
time and settings — before either external file is opened or hashed**; (6) only then hash FAS_H/DEMO_H, record the
input hashes on the attempt, refuse if they differ from a previous attempt's, load the cycle, assert SEQN
disjointness from the development SEQN hash list, compute picks, c_b, C, Δs and S3; (7) write the result and mark the
attempt `completed`. Any failure after step 5 marks the attempt `failed` with the error text; **a retry is permitted
only with identical inputs and settings**, is itself recorded, and nothing is re-seeded, re-allocated or re-fitted.
Guard behaviour (steps 1–5, failed-attempt recording, identical-retry rule, tamper refusal) was exercised on a
synthetic fit with absent external files. Per-participant tables and NHANES files stay outside the repository.

**Order of events.** `fit` → `bench` (record) → `predict` (complete table) → `assemble` → `screen` → owner review
→ commit protocol, runner, prediction, fit aggregate, screen → `freeze` → commit → owner writes the authorisation
file → FAS_H/DEMO_H downloaded once → `run-ext` → `compare` → `NHANES_FA_STABILITY_RESULT_v1.0.md` (realised
n_eligible, c_b, C, [L, U] at that n, tail proportions, the §5 sentence that applies, S1–S3, file hashes, attempt
record, deviations).

## 7. Corrections to adjacent records (owner's finding, to be made before freezing)

Block 3, fold 1 of the 2011–2012 run: NOVA-V and NOVA-U chose `LBXMR1`, G1 and NOVA-S `LBXDTA`; agreement with G1
24/25 for NOVA-V and NOVA-U, 25/25 for NOVA-S; G1's within-block disagreement 0; development contrast and NO-GO
unchanged. `NHANES_FA_RESULT_v1.1.md` (owner's draft) carries the correction with v1.0 preserved; matching
corrections at `NOVA_PAPER_v0.5.md:377` and `SELECTOR_CLAIM_EVIDENCE_MATRIX_v0.4.md:119`.

## 8. Inherited and new

Inherited unchanged: eligibility, 29-analyte panel, five block allocations, fold seeds (20261103 dev, 20261104 ext),
`g1_pick`, `folds_for`, `inner_scores` and the claim constants, `load_cycle`, `arrays`. New: endpoint on integer
counts (§1), plug-in models with exact-size predictions and the discrete quantile convention (§2), in-sample
screens with the sharpness requirement and the stored assessable set (§3), frozen interpretation (§5), the separate
runner with `g1_scores`, `frozen_cycle`, the picks-only NOVA path and the pre-read external guard (§6), the
authorisation file (§0, §6).

## 9. Changes from v0.3 (owner's approval of 10 November 2026)

1. Calibration of the interval stated precisely: [L, U] is a Monte Carlo estimate of the model's 95 % plug-in
   prediction set; tail quantities are Monte Carlo proportions under the fitted model, not p-values; no allowance
   for fit uncertainty (§2).
2. The assessable-size set is computed by `screen`, stored explicitly as integers and hashed in the freeze; a
   realised size inside the verdict range but outside the set is not assessable (§3b, §5).
3. Workload benchmarked (`bench`) before launching the table; reference timing recorded; `predict` resumable in
   parts; a partial table is descriptive only, is refused by `screen`, and cannot trigger a freeze or an opening;
   any change to range or B is a prospective amendment made before the realised size is known (§2).
4. External guard moved before the first read: authorisation tied to the frozen hashes; the opening attempt is
   recorded before either file is opened or hashed; failed attempts recorded; retries only with identical inputs and
   settings (§6).
5. §1 wording: held-out folds are participant-disjoint; the training sets share about 60 % of their rows.
6. Implementation details recorded as built: `frozen_cycle` for the development reproduction and the external
   picks; `fit` refuses if the recorded 25 picks are not reproduced; participant-level files named and kept local;
   test list updated to what the runner actually runs.
