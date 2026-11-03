# N2 failure-regime diagnostic: result record v1.0 (3 Nov 2026)

Classification: **exploratory, descriptive; archive use 11**; no selector trained; no matrix row read; one completed
run. Verdict under the prestated readings of spec §7: **"No regime identified."** No covariate defines a subgroup in
which departures from G1 reach the 0.29 % floor with a positive lower bound. The next external-validation hypothesis
is therefore **neither** nonlinear structure nor record-level G1 stability as operationalised here. One block-level
pattern (A4′) is suggestive but is confounded by construction and is recorded as a hypothesis only.

## 1. Provenance

| Item | Value |
|---|---|
| Specification | `N2_REGIME_DIAGNOSTIC_SPEC_v0.3.md` SHA-256 `0d36d518917149e27891780de8b4ed0e6c6b5623b63d0fb530e2ac82609a7726` |
| Script | `n2_regime.py` v0.1.1 SHA-256 `f3e11962aeb764a9b336432cf3ae7bbf1dbfb772d08ff55261e1140ab0e3915d` (`n2_regime_freeze_hashes.txt`) |
| Inputs (hash-checked by the script) | `n2_blockfolds.csv` `b7934523…`; `gate3a_blocks.jsonl` `6decb9fb…`; `n0_records.csv` `6d3b0341…`; 186 `design_model_manifest.json` (hashes in the records CSV) |
| Output | `n2_regime_out/n2_regime_summary.json` SHA-256 `5b9c42b509b8ba6f8790b53208017ff2ded860ad4d5583eaa7d31504f2b268da`; `n2_regime_records.csv` (kept with the development data) |
| Reproduction gate | N2 primary reconstructed from the block-fold file: **+0.11179 %** (frozen: +0.1118 %) — passed |
| Self-tests | 6/6 passed on the owner's machine before the freeze |
| Run | 2026-11-03T20:23:29Z; seed 20261103; 10,000 cluster resamples; 99 % two-sided intervals (Bonferroni over five) |

**Disclosure.** `n2_regime.py` v0.1.0 (SHA-256 `9719dd52…`) was frozen first and crashed in A2 on the real run
(`LinAlgError: SVD did not converge`) before writing any output; the A1 descriptives had been computed in memory and
discarded. v0.1.1 changed only the solver path (normal equations with a least-squares fallback), the finite-value
mask, and the order in which the records CSV is written. The estimator is unchanged. Both hashes are in the freeze
file. The crash cause was not a non-finite value in the data (v0.1.1 reports 0 non-finite entries in every covariate
and in D) and is attributed to the LAPACK backend on the run machine; this is not verified.

## 2. Units

179 of 186 records had at least one block-fold with j_V ≠ j_G1 and enter the analysis; 7 records had none. No two
records share a study accession, so the study-cluster bootstrap reduces to a record bootstrap (179 clusters). Two
records have a non-finite mean t (block-folds with t = ∞, a zero paired standard error) and are excluded from A2/A3
(177 records used); they remain in A1, A4′ and A6.

## 3. Results

**A2, partial association with realised gain D given t̄** (per-SD coefficients; 99 % cluster-bootstrap intervals;
condition number 5.5; largest covariate correlation 0.28, between `log_n` and `g1_stability`):

| Term | Coefficient | 99 % interval | Excludes 0 | Prespecified direction |
|---|---|---|---|---|
| t̄ (context) | +0.00965 | [+0.00511, +0.01479] | yes | — |
| C1 `log_n` | +0.00180 | [+0.00010, +0.00453] | **yes** | + (H-nonlinear) — sign matches |
| C2 `panel_class` | +0.00001 | [−0.00178, +0.00152] | no | two-sided |
| C3 `skew_O` | −0.00023 | [−0.00157, +0.00170] | no | two-sided |
| C4 `edge_O` | −0.00050 | [−0.00212, +0.00182] | no | + (H-nonlinear) — not supported |
| C5 `g1_stability` | −0.00232 | [−0.00582, +0.00072] | no | − (H-stability) — sign matches, interval includes 0 |

**A3, association with t̄**: no covariate predicts the paired statistic (all five intervals include 0). None is
"predicts departure only".

**A5, prespecified median split for C1** (high-n side, 87 records): Δ_acq of departures = **+0.083 %**, 99 % lower
bound **−0.019 %**. Fails the floor (0.29 %) and the positive-lower-bound condition. **A1** shows what the C1 effect
is: the low-n tertile (63 records) has D = −0.0028 [−0.0048, −0.0005] (NOVA-V's departures are harmful at small n),
the high-n tertile +0.0029 [−0.0008, +0.0066]. This is the GDOA-0 pattern (small-n estimation noise) seen from the
other side; it is not a regime in which departures help.

**A4′, inside/outside the G1-plausible set** (jV ∈ set of G1's picks over the block's five folds):

| Side | Block-folds (departures) | Records | D, all differing [99 %] | D, departures [99 %] | Contribution to primary |
|---|---|---|---|---|---|
| inside | 726 (134) | 134 | **+0.0068 [+0.0045, +0.0096]** | +0.0101 [+0.0058, +0.0157] | +0.050 % |
| outside | 2,339 (368) | 173 | −0.0004 [−0.0043, +0.0037] | +0.0088 [−0.0026, +0.0233] | +0.062 % |

The two contributions sum to the primary (+0.112 %). The inside-set block-folds show positive realised gain even
without the t-gate; the outside set shows none without the gate and an uncertain gain with it.

**Confound, stated before any interpretation.** The G1-plausible set of a block is built from G1's picks on the
other four folds, each of which was fitted on training rows that *include* the scored fold. A candidate whose
correlation with the target is high in the scored fold is therefore more likely to be in the set *and* more likely
to have low loss on that fold. The inside stratum is thus selected with partial sight of the scored outcomes, and
its positive D is biased upward by an amount this run cannot measure. A4′ is **not** evidence that nonlinear value
resolves G1 near-ties; it is a pattern that a leak-free design (G1 sets built from nested folds that exclude the
scored fold entirely) would have to test. The spec's author (Claude) introduced this confound in v0.3 and did not
catch it before the freeze; it is recorded here as a design error of the diagnostic, not of the frozen N2 result.

**A6, t ≥ 4**: 16 block-folds from 9 records (the replay's 13 excluded the three with t = ∞); 12 outside, 4 inside;
realised gains from −28.5 % to +17.2 %; six of the sixteen come from one record (ST002081_AN003790, n_design
large, `edge_O` 0.094). No inference.

## 4. Readings (spec §7, mechanical)

- C1 `log_n`: passes A2 with the H-nonlinear sign, fails A5 → **not a regime**; the A1 tertiles show the association
  is driven by harm at small n, consistent with GDOA-0.
- C2, C3, C4, C5: **no association** with gain given t̄ and none with t̄.
- Stability reading not reached (C5 did not pass A2). A4′'s inside-set result would have been the qualifier; it is
  confounded and cannot serve.
- **Overall: "No regime identified."**

## 5. What follows

1. The NOVA archive paper's framing stands: a small development gain from selective departures, below the floor,
   with no pre-acquisition covariate that predicts where departures help. H-nonlinear is not supported by C4
   (the one covariate that measures nonlinear structure in O directly) and H-stability is not supported by C5 at
   the record level.
2. The only hypothesis this diagnostic generates — *nonlinear value as a tie-breaker among the candidates the
   covariance criterion itself finds plausible* — rests on a confounded stratification and may be stated only as a
   hypothesis with its confound named. Testing it on the archive would require a nested-fold design (archive use 12)
   and would still be development evidence; testing it externally requires a cohort, which the programme does not
   have. The owner's decision rule from 3 Nov ("nonlinear structure, selection stability, or neither") resolves to
   **neither**, on this evidence.
3. No new selector is proposed from this archive. No frozen number changes.

## 6. Files

`N2_REGIME_DIAGNOSTIC_SPEC_v0.3.md`, `n2_regime.py`, `n2_regime_freeze_hashes.txt`, `n2_regime_run.log`,
`n2_regime_out/n2_regime_summary.json` (commit); `n2_regime_out/n2_regime_records.csv` (local, with the development
data; not committed).
