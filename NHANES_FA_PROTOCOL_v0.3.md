# NHANES serum fatty-acid acquisition test: protocol v0.3 (for freezing)

Programme: decision-optimal metabolite acquisition (NOVA-S vs G1). Date: 3 November 2026. Owner: Anas Enoch.
Status: v0.3 of 4 Nov 2026 (changes in §12); frozen when hashed with the wrapper and the block in `NHANES_FA_FREEZE_HASHES.txt`. It is frozen when `n2_regime`-style hashes are recorded for this file and for the
wrapper `nhanes_fa_run.py` in `NHANES_FA_FREEZE_HASHES.txt`, **before any NHANES data file is downloaded**.
Classification of the eventual result: **restricted-panel, single-unit, supportive** test on a previously unused
public survey cycle. It is not the P4D 20/20/20 benchmark, not a P4V confirmatory unit, and not the test described
in `EXTERNAL_VALIDATION_PREREG_STUB_v0.1.md`, whose block rule it cannot satisfy. LIPID/COMETS (deadline 15 Dec
2026) remain the primary external route.

## 0. Provenance boundary

- Prior programme contact with NHANES: Stage 16L–16O (11 Oct 2026) read NHANES 2017–2018 `PHTHTE_J` (urinary
  phthalate metabolites) for an analytical-status audit; no prediction or selector was run; Stage 16M blocked any
  transport of the frozen Paper IV study-specific models. No programme file has read FAS_G or FAS_H [verified by
  repository grep, 3 Nov 2026].
- The 186-record development archive is closed (11 uses). Nothing in this protocol reads it.
- Documentation read before drafting (metadata only): CDC pages for FAS_G (2011–12) and FAS_H (2013–14), and the
  FAS_G laboratory method PDF. Both cycles release the same 30 analyte variables with the same labels, in µmol/L,
  by the same GC-ECNI-MS method; FAS_H states no method change. The FAS_G method text mentions 35 measured fatty
  acids; **the 30 released variables are the universe**; the discrepancy is recorded and not resolved here.

## 1. Hypothesis and estimand

H1: at acquisition budget B = 1, the frozen NOVA-S selector yields lower held-out predictive loss than the frozen G1
selector on NHANES 2013–2014, under identical blocks, folds, preprocessing and predictor.

**What is tested.** Both selectors are refit within each cycle on its training folds; nothing fitted on 2011–12 is
transported. The experiment is therefore a **procedural replication** of the frozen acquisition procedure on a new
survey cycle, not a transport of fitted models, and the manuscript describes it so. The five outer folds are
correlated evaluations within one cycle, not five studies.

Per participant i in the evaluation cycle, over the evaluation cells of the folds in which i is held out:
a_i = Σ_cells (L_G1 − L_S), c_i = Σ_cells L_G1 (floored squared loss in design-standardised target units, the N2
loss). **Primary estimand** Δ_ext = Σ_i w_i a_i / Σ_i w_i c_i with w_i = WTFAS2YR (a weighted relative loss
reduction; positive favours NOVA-S). **Baseline improvement** Δ_G1→∅ = 1 − Σ_i w_i Σ_cells L_G1 / Σ_i w_i Σ_cells L_O, where L_O is the loss of the same
frozen predictor fitted on the base panel O alone (the "O" arm of `p6_n2.py`), same folds and cells; its design-based
one-sided 97.5 % lower bound is computed by the same linearisation. Secondary: the unweighted ratio; the bounded
loss min(e², 25); the equal-fold mean; the agreement rate of the two selectors across folds and the NOVA-S departure
rate.

## 2. Data, eligibility, identity (fixed before download)

| Item | Rule |
|---|---|
| Files | 2011–12: `FAS_G.xpt`, `DEMO_G.xpt`; 2013–14: `FAS_H.xpt`, `DEMO_H.xpt` (CDC public). SHA-256 of each file recorded at download; `FAS_H`/`DEMO_H` are not downloaded until the GO of §6. |
| Variables read | FAS: SEQN, WTFAS2YR, the 30 `LBX*` analytes and their `LBD*LC` flags. DEMO: SEQN, RIDAGEYR, SDMVSTRA, SDMVPSU. No other variable. |
| Eligibility | RIDAGEYR ≥ 12 and WTFAS2YR > 0 and all retained analytes non-missing. Children 3–11 are excluded (different sampling rule and no fasting subsample). |
| Identity | one row per SEQN; a duplicate SEQN in either file aborts the run. SEQN is not shared across cycles (NHANES samples new participants each cycle), so no participant overlap is expected; the run asserts an empty intersection of SEQN sets. |
| Below-LLOD values | kept as released (CDC fill LLOD/√2) for retained analytes; the N2 pipeline's log rule applies to inputs as frozen. |
| Survey design | strata SDMVSTRA and PSU SDMVPSU of each cycle are carried into the variance estimator (§5). |

## 3. Analyte panel and block allocation (fixed before download; computed from codebook counts only)

- **Exclusion rule**: an analyte is excluded if its codebook count of below-LLOD results exceeds 5 % of results
  with values in **either** cycle (FAS_G LBXCAP: 1,862/3,873 = 48 % → excluded; FAS_H LBXLAR 116/3,931 = 3.0 % →
  retained unless FAS_G's count exceeds 5 %; the owner reads the remaining FAS_G counts from the codebook and lists
  every retained analyte in §3.1 before hashing). Owner's codebook check (3 Nov): LBXCAP excluded (48 %); LBXLAR retained (3.98 % FAS_G, 2.97 % FAS_H);
  LBXSD1 retained (0.75 % FAS_G). **Expected panel: 29 analytes** (SFA 10, MUFA 6, PUFA 13), subject to completing the
  check for all 30 in both cycles and listing them in §3.1.
- **Allocation** (one block, deterministic, exact role sizes). With m retained analytes the role sizes are
  |O| = ⌈m/3⌉, |T| = ⌈(m − |O|)/2⌉, |P| = m − |O| − |T| (m = 29 → 10 / 10 / 9). Algorithm, written to
  `nhanes_fa_block.json` and hashed before download:
  1. Within each CDC class (SFA, MUFA, PUFA) sort analytes by SHA-256 of `NHANES-FA-v0.2|<LBX variable name>`,
     ties (none expected) by variable name.
  2. Concatenate the three hash-ordered class lists in the fixed order SFA, MUFA, PUFA into one sequence.
  3. Assign position k (0-based) of that sequence to role k mod 3 with role order (O, T, P). Because consecutive
     positions cycle through the roles, each class is split across O, T, P as evenly as possible (per-class role
     counts differ by at most 1) and the global sizes are exact (m = 29 → 10 / 10 / 9).
  Tested on synthetic class counts 10/6/13, 7/7/7, 1/1/27, 10/6/12, 11/6/13: exact sizes and within-class spread
  ≤ 1 in every case (3 Nov). A first version that *interleaved* classes was rejected by the same test: with three
  classes and three roles it mapped each class to one role. The allocator uses no measured value and is not
  optimised for either selector.
- 3.1 Retained analytes: the 29 listed in `nhanes_fa_analytes.json` (hashed with this protocol). Owner's codebook
  verification (4 Nov 2026), below-LLOD counts FAS_G / FAS_H: capric 1,862/3,873 and 1,531/3,931 → **excluded**;
  lauric 169/4,249 and 116/3,904; myristoleic 1/3,911 and 1/4,311; stearidonic 29/3,866 and 11/3,970;
  eicosatrienoic 6/3,782 and 0/4,261; EPA 1/4,421 and 0/4,339; all other analytes 0 in both cycles → retained.
  Roles are those produced by `nhanes_fa_run.py alloc` from that file (`nhanes_fa_block.json`, hashed).

## 4. Frozen procedure

- Selectors and predictor: `p6_n2.py` v0.1.0, SHA-256 `681ce852609e6ab4fcaf0f0ff603477c64e2351e10bd119811d45f008d4bbe8e`,
  under `N2_PROTOCOL_v1.0.md` `6bbafb3c…`: G1 (Paper IV plug-in rule, re-estimated per training fold), NOVA-V/U/S
  with the RBF nuisance, N2 grids, log rule, floored loss, **τ = 2**. The decision rule, grids and threshold are not
  changed for NHANES. The wrapper `nhanes_fa_run.py` imports the frozen module, builds the single block from
  `nhanes_fa_block.json`, forms participant-disjoint folds, and calls the frozen per-block routine; it adds no
  selection logic. It is written against the frozen code and tested on synthetic data **before** any NHANES file
  exists on the machine, and its hash is recorded with this protocol.
- **Target standardisation (frozen N2 convention).** Each target is divided by its standard deviation over all
  eligible observations of the survey cycle, subject to the original 1 % median-standard-deviation floor
  (s_j = max(SD(Y_j), 0.01 · median_k SD(Y_k))). Scaling is performed before outer-fold construction and is
  identical across the G1, NOVA-S and observed-only arms. This preserves the preprocessing convention of the frozen
  N2 implementation. Because the scale estimates include observations subsequently assigned to held-out folds, the
  resulting evaluation is not strictly fold-isolated with respect to target normalisation; the convention is
  symmetric across arms and selection-independent, but the held-out targets' scales do set the relative weighting
  of targets in the loss. The estimated scaling constants are reported in the result record. Fold-specific scaling
  may be run only as a separately identified sensitivity analysis with its own code path; it never replaces the
  primary comparison.
- Folds: 5 outer folds over eligible SEQN (seed 20261103 for 2011–12; seed 20261104 for 2013–14), fixed by hash
  of SEQN so they are reproducible; inner tuning as in `p6_n2.py`. Both selectors see identical training rows.
- One run per cycle. A crashed run may be repeated with identical settings. The 2013–14 run happens once.

## 5. Variance and inference (fixed)

Design-based Taylor linearisation for the ratio Δ: with u_i = w_i (a_i − Δ c_i) / Σ_j w_j c_j, the variance is the
stratified with-replacement PSU estimator over (SDMVSTRA, SDMVPSU) of the cycle, degrees of freedom = PSUs − strata.
One-sided 97.5 % lower bound on Δ from the t distribution with those df.
- **Domain estimation.** Eligibility (§2) defines an analysis domain, not a new sample: every record of the cycle's
  FAS file keeps its stratum and PSU, ineligible records enter the linearisation with u_i = 0, and the PSU and
  stratum counts are those of the full file. Strata or PSUs are never dropped by filtering.
- **Interpretation.** The interval is conditional on the frozen cross-fitting procedure (fixed folds, fixed grids,
  fixed threshold): it measures survey sampling uncertainty of the weighted loss contrast for this procedure as run,
  not variability from alternative fold assignments or retraining. Procedure instability is reported separately and
  descriptively: the agreement and departure rates by fold, and the per-fold Δ values. No second seed is run.
- Sensitivity: unweighted participant-cluster bootstrap (10,000 resamples), reported alongside, never substituted.

## 6. Development cycle and the feasibility decision (2011–2012)

1. Download `FAS_G`, `DEMO_G`; record hashes; apply §2–§3; record n_eligible.
2. Run the frozen procedure once. Record Δ_dev, its design-based SE_dev and df, Δ_G1→∅,dev, agreement and departure
   rates. 2011–12 is **development data** from this point: it may be inspected freely, but no change to the
   selector, threshold, grids, loss, panel rule or allocation may follow from it; any defect found in the wrapper
   is fixed, documented, and the fix hashed before §7.
3. **MDE** (defined now): MDE_dev = (t_{0.975,df} + t_{0.80,df}) × SE_dev, i.e. the effect detectable with
   one-sided α = 0.025 and power 0.80 under the assumption that SE in 2013–14 equals SE_dev (similar n and design;
   the assumption is stated in the result record and checked after the fact).
4. **GO/NO-GO for opening FAS_H**: GO iff n_eligible,dev ≥ 1,500 and MDE_dev ≤ 0.30 percentage points. NO-GO
   otherwise; FAS_H is then not downloaded and the record states "NHANES-FA assessed; infeasible under v0.1".
   The ceiling 0.30 is not revised after MDE_dev is seen. Δ_dev's sign and size play no part in the decision.

## 7. Evaluation cycle (2013–2014), once

Download `FAS_H`, `DEMO_H`; record hashes; apply the identical §2–§3; assert SEQN disjointness from 2011–12; run the
frozen procedure once; compute §1 and §5.

**Success criterion (supportive), fixed now, three conditions in order:** (i) the baseline improvement
Δ_G1→∅,ext has a design-based one-sided 97.5 % lower bound > 0 (the measurement is worth something in this panel;
if not, the comparison is reported but no success is declared, whatever Δ_ext is); (ii) the design-based one-sided
97.5 % lower bound of Δ_ext > 0; (iii) Δ_ext ≥ 0.10 × Δ_G1→∅,ext. Failure: any condition unmet. Both are reported with the archive's 0.29 % floor as context only.
The result is reported whatever it is, with Δ_dev beside it; no subgroup, re-threshold, re-allocation or re-run is
permitted after FAS_H is read; anything further is labelled exploratory.

## 8. What a result means

- Positive: the programme's first supportive external acquisition result, on a restricted, strongly correlated
  fatty-acid panel in one survey population; it would need independent replication and a properly tuned G1
  (shrinkage α and covariance estimator) before any superiority claim.
- Negative: independent evidence that nonlinear predictive information does not translate into a better measurement
  choice in this panel; consistent with the archive. Not evidence against the population identities.
- Either way: no claim beyond serum total fatty acids by GC-MS in the NHANES fasting subsample aged 12+.

## 9. Integrity

No NHANES value is read before the hashes of this protocol, `nhanes_fa_run.py` and `nhanes_fa_block.json` are
recorded. Files read are limited to §2. FAS_H is not downloaded before the §6 GO is recorded with Δ_dev, SE_dev and
MDE_dev. The result record lists every deviation. All outputs (feasibility record, result record, run manifests
with input hashes, summary JSONs) are committed; per-participant loss tables stay local.

## 10. Open items before hashing

(a) §3.1: complete the 5 % check for all 30 analytes in both cycles and list the 29 with roles (owner);
(b) inspect `p6_n2.py` (copied into `nhanes_fa/frozen/`, hash 681ce852… verified) for a callable per-block routine;
if none exists, the adapter must reproduce the frozen procedure without changing selection rule, numerical
conventions or tuning, and a **synthetic equivalence test against the original implementation is mandatory**
(same synthetic block through both paths, identical picks and losses to 1e-10) before hashing;
(c) downloads: FAS_G/DEMO_G authorised by the owner (3 Nov) **after** this protocol, `nhanes_fa_block.json`, the
wrapper and its tests are hashed; FAS_H/DEMO_H authorised only on the §6 GO;
(d) constants confirmed by the owner (3 Nov): MDE ceiling 0.30 pp; success factor 0.10 with the baseline condition.

## 11. Changes from v0.1

1. Analyte check completed for capric, lauric, stearidonic; expected panel 29 (SFA 10 / MUFA 6 / PUFA 13).
2. Allocator corrected: per-class restart gave 11/9/9; class interleaving aliased classes to roles; the concatenated
   class-ordered sequence dealt (O, T, P) gives exact 10/10/9 with per-class spread ≤ 1; tie rule and tests specified.
3. Variance: domain estimation preserves strata/PSUs; interpretation conditional on the frozen procedure;
   procedure instability reported descriptively.
4. Success criterion gains the baseline condition Δ_G1→∅ lower bound > 0, and Δ_G1→∅ is defined with its denominator.
5. Procedural-replication framing stated in §1.
6. Wrapper equivalence test made mandatory; download authorisation sequence recorded.

## 12. Changes from v0.2

1. §4: target-standardisation paragraph adopted in the owner's approved wording (frozen N2 convention as primary;
   non-fold-isolation stated; constants reported; fold-specific scaling only as a separate sensitivity).
2. §3.1: analyte verification counts recorded; 29 retained; the machine-readable table is `nhanes_fa_analytes.json`.
3. Implementation reference: `nhanes_fa_run.py` v0.1.1 (v0.1.0 629d2c67… differed only in the default protocol path), 15 self-tests including exact
   equivalence with the frozen `run_block`; frozen copies `nhanes_fa/frozen/p6_n2.py` (681ce852…) and
   `p5_gate3a.py` (4600e2f3…).
