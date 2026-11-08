# G1-R — Covariance regularisation and finite-sample acquisition decisions: result record v1.0 (8 Nov 2026)

**Status: frozen study executed once; archive closed.** Protocol `G1_REGULARISATION_PROTOCOL_v0.5.md` (SHA-256
`cf1b08ad721b5821403583c4d3a42e76faf2a7f0f0022b6cef0736b4b8214b27`); implementation `g1r_run.py` v0.3.0
(`93026388861c758799ed572b690789d0a5813b7bc9a20c77daf32fc288af333f`), `g1r_integrity.py` v0.2.0
(`c59c92a41b33864db41b9fe7163652156d7217acc8f1b4d1a4e370a1e3e7bad2`), `nova_d_calib.py`
(`4ab49dbc585e5cef4da6d0f2a31ca0d44725ffbaf0364f058fd8b63615b1ce29`), frozen comparator `nhanes_fa/frozen/p5_gate3a.py`
(`4600e2f3e3bb652affadcc6f6935380d8b9a77a4237e7cb306f056f849478130`). Implementation and protocol committed at
`4f0f1fba9275846881bacb7880170ff794f71fdc`; freeze record and bound self-test committed at
`833c07292e3803698be2150296c4e81cc1aacf23`. Nothing in §4–§12 of the protocol was changed after any result.
Author of record: Anas Enoch. Drafted by Claude (Fable 5.1) from the frozen outputs; corrections of 8 Nov 2026 by the
owner incorporated (§1, §3, §5).

---

## 1. Decision — the prespecified joint criterion is not met

Protocol §2 and §10 require **both** H1 at D1 and H2 at D2. The frozen `analyze` output (`g1r_out/g1r_result.json`,
SHA-256 `2fc9ac74ab33626947cbc051563dc9a70327bd08ca799a1154524eac00fe37f7`) records:

| Hypothesis | Designated condition | Statistic (margin inside) | Rule | Realised | Outcome |
|---|---|---|---|---|---|
| H1 | D1 = condition 110 | Q1 = mean(0.9·r₀ − r₂) | LB₅(Q1) > 0 | mean 0.000984, LB₅ **0.000888**, 95 % CI [0.000871, 0.001125] | **supported** |
| H2 | D2 = condition 104 | Q2 = mean(1.1·r₀ − r₂) | LB₅(Q2) ≥ 0 | mean −0.000174, LB₅ **−0.000408**, 95 % CI [−0.000456, +0.000102] | **not established** |

**H2 is not established; it is not refuted.** The one-sided lower bound is negative, so the no-harm criterion
(decoupled regularisation no worse than 1.1 × A0's regret at the designated non-flip condition) fails. The two-sided
95 % interval for Q2 crosses zero, so harm beyond the 10 % margin is not established either. The joint claim of the
protocol — that the decoupled conditional regularisation (A2) improves on the frozen G1 rule (A0) where A0's
shrinkage floor reorders the population argmax **and** does not harm it where it does not — is therefore
**not supported** by this study. A0 is retained as the comparator; the frozen G1 rule is not changed by this result.

`decisions` block, verbatim: `H1_D1: true`, `H2_D2: false`, `n_conditions_with_assessable_A2_contrast: 288`,
`grid_max_plain_A2_gain_upward_biased: [0.028375510489736855, 123]`.

## 2. Execution record

| Step | Realised | Identity |
|---|---|---|
| Self-test regenerated (owner's machine) | PASS, v0.3.0; D1 = [110, 1], D2 = [104, 0]; A1 ≡ A0 max relative discrepancy C 3.2e-15, v 8.5e-16, s 2.0e-15, 0 argmax mismatches; six moment checks OK; 1.25 ms / 6.88 ms per replicate at n 102 / 3280 | `g1r_out/selftest.json` `eeb12459b8f095c188ad49a15af7545656bedbd8d449d5e46ec21907420ec729` |
| Freeze | `frozen_utc 2026-11-07T22:21:00Z`; binds protocol, four file hashes, self-test checksum, salt `G1R-V02|`, 360 conditions × 1000 replicates, Λ of 82 points, margins 0.9/1.1, B = 2000, sparse minimum 10, gate tolerance 1e-12, moment conditions 110/104/1/359/204/225 | `g1r_out/FREEZE.json` `885cd87d9bcaeda6a610d1101ff7cf83109f72691734b5c473801fe65c9f8013` |
| Run (owner's machine, bridge shell; 14 resumable `--cond` chunks, 2026-11-08 00:32–00:53 UTC) | exit 0 on every chunk; no `GateError`; after each chunk the line count equalled the expected total and the last line parsed; a final full-range resume call wrote 0 records | `g1r_out/replicates.jsonl`, 360,000 records, **canonical uncompressed identity** `ea33587ecd641b15b43b01052db6174ee9cd15b680da8845d172ab0f90e33d19` (198,074,286 bytes); archive `g1r_out/replicates.jsonl.gz` `45b0c66582bb281e1cc95b823992ae40b91736db6210e0cf2fd00fd225b96be0` (22,231,428 bytes; `gzip -9 -n`; decompresses to the canonical hash, verified by owner and drafter) |
| Integrity (`g1r_integrity.py` v0.2.0) | **PASS**: 360,000 lines, 360,000 valid, 0 unparsable, 0 invalid, 0 duplicate, 0 missing, 0 extra, 0 incomplete conditions, 0 regret/pick mismatches; bound to log `ea33587e…`, freeze `885cd87d…`, checker `c59c92a4…` | `g1r_out/integrity.json` `d688b2d8a1c28a5b7c19b16d56ab90efa38337f87b4ba1098ea8fada6c42de56` |
| Analyze | exit 0 after verifying freeze, self-test checksum and the three integrity bindings | `g1r_out/g1r_result.json` `2fc9ac74ab33626947cbc051563dc9a70327bd08ca799a1154524eac00fe37f7` |
| Independent recomputation (drafting container, from the archive; round-trip hash `ea33587e…` verified) | script takes the log as data and the result file as the object under test, importing nothing from `g1r_run.py` for statistics (Clopper–Pearson via the Beta quantile; bootstrap re-implemented from the protocol's seed and consumption order on the paired differences): 14,400 arm statistics, all 6,723 assessable bootstrap contrasts (mean, LB₅, CI endpoints to 1e-12 relative), 360 Q1/Q2 margin identities, 720 secondary-contrast slots, decisions block — **0 mismatches** | `g1r_analyze_independent_check.py` `42862725831d0bf445f8f7ebb04431371c41ce25a6bbcc3dd32c34d8c27f1e4b` |

Chunk wall times (conditions → minutes): 0:72 1.5; 72:108 1.0; 108:144 1.0; 144:174 1.3; 174:216 1.7; 216:236 1.3;
236:262 1.7; 262:288 1.7; 288:300 1.4; 300:314 1.6; 314:328 1.6; 328:342 1.6; 342:356 1.6; 356:360 0.5.

## 3. Designated conditions in full

**D1 = condition 110** — n_train 410, φ 0.8 (κ(Σ_OO) = 64.7), Γ 0.02, A 0.15, equal variances; j\* = 0. First entry of
the frozen D1 list at which A0's population argmax flips (Objective A).

| Arm | mean regret | P(ĵ = j\*) | replicates with r > 0 | CP₉₅ upper on misselection |
|---|---|---|---|---|
| A0 (frozen G1) | 0.002574 | 0.178 | 822 | 0.842 |
| A2 (decoupled, GCV) | 0.001332 | 0.568 | 432 | 0.458 |
| A2-S (Schur form) | 0.001512 | 0.508 | 492 | 0.519 |
| A2′ (population-benchmark λ) | 0.001320 | 0.572 | 428 | 0.454 |
| A3 (floor removed) | 0.002574 | 0.178 | 822 | 0.842 |
| A4 (Ledoit–Wolf) | 0.001644 | 0.488 | 512 | 0.538 |

A2 vs A0: plain mean 0.001242, LB₅ 0.001128, CI [0.001113, 0.001391]; **Q1** mean 0.000984, LB₅ 0.000888,
CI [0.000871, 0.001125]; Q2 mean 0.001499, LB₅ 0.001382, CI [0.001363, 0.001660]. Secondary: form effect
mean(r_A2 − r_A2-S) = −0.000180 (LB₅ −0.000216; A2 better than A2-S here); tuning effect mean(r_A2 − r_A2′) = 0.000012
(LB₅ 0.000003).

**D2 = condition 104** — n_train 410, φ 0.5 (κ = 8.63), Γ 0.20, A 0.15, equal variances; j\* = 0. First entry of the
frozen D2 list at which A0's population argmax does not flip.

| Arm | mean regret | P(ĵ = j\*) | replicates with r > 0 | CP₉₅ upper on misselection |
|---|---|---|---|---|
| A0 | 0.000960 | 0.968 | 32 | 0.043 |
| A2 | 0.001230 | 0.959 | 41 | 0.053 |
| A2-S | 0.001080 | 0.964 | 36 | 0.047 |
| A2′ | 0.001230 | 0.959 | 41 | 0.053 |
| A3 | 0.000960 | 0.968 | 32 | 0.043 |
| A4 | 0.001050 | 0.965 | 35 | 0.046 |

A2 vs A0: plain mean −0.000270, LB₅ −0.000510, CI [−0.000540, 0.000000]; Q1 mean −0.000366, LB₅ −0.000585,
CI [−0.000639, −0.000111]; **Q2** mean −0.000174, LB₅ −0.000408, CI [−0.000456, +0.000102]. Secondary: form effect
+0.000150 (LB₅ 0.000060; A2-S better than A2 here); tuning effect 0.000000.

Both designated conditions sit at n_train = 410 and equal variances; the decision therefore concerns that regime and
is not a statement about the grid.

## 4. Grid summaries — descriptive only

Every condition is reported in `g1r_result.json`. The summaries below carry no inferential weight: no multiplicity
control was prespecified for them, the grid maximum is selection-biased by construction, and the decision rule is §1.

**Mean exact population regret, averaged over the 72 conditions at each n_train**

| n_train | A0 | A2 | A2-S | A2′ | A3 | A4 | A2 β=0.5 | A2 β=1 |
|---|---|---|---|---|---|---|---|---|
| 102 | 0.00916 | 0.00589 | 0.00619 | 0.00588 | 0.00916 | 0.00716 | 0.00788 | 0.01136 |
| 410 | 0.00376 | 0.00223 | 0.00230 | 0.00224 | 0.00376 | 0.00273 | 0.00396 | 0.00670 |
| 820 | 0.00244 | 0.00139 | 0.00147 | 0.00139 | 0.00244 | 0.00172 | 0.00335 | 0.00599 |
| 1640 | 0.00146 | 0.00090 | 0.00091 | 0.00091 | 0.00125 | 0.00101 | 0.00382 | 0.00783 |
| 3280 | 0.00145 | 0.00059 | 0.00059 | 0.00058 | 0.00086 | 0.00067 | 0.00341 | 0.00662 |

**P(ĵ = j\*), averaged over conditions**

| n_train | A0 | A2 | A2-S | A2′ | A3 | A4 |
|---|---|---|---|---|---|---|
| 102 | 0.410 | 0.483 | 0.475 | 0.484 | 0.410 | 0.446 |
| 410 | 0.598 | 0.643 | 0.645 | 0.643 | 0.598 | 0.637 |
| 820 | 0.629 | 0.716 | 0.706 | 0.717 | 0.629 | 0.683 |
| 1640 | 0.698 | 0.780 | 0.776 | 0.780 | 0.719 | 0.757 |
| 3280 | 0.706 | 0.834 | 0.833 | 0.834 | 0.786 | 0.821 |

**A2 vs A0 across assessable conditions** (288 of 360)

| n_train | assessable | plain mean > 0 | LB₅(plain) > 0 | UB₉₇.₅(plain) < 0 | LB₅(Q1) > 0 | LB₅(Q2) ≥ 0 | mean of plain means |
|---|---|---|---|---|---|---|---|
| 102 | 72 | 51 | 48 | 16 | 45 | 51 | 0.00327 |
| 410 | 68 | 37 | 36 | 26 | 34 | 39 | 0.00167 |
| 820 | 54 | 36 | 34 | 16 | 32 | 38 | 0.00098 |
| 1640 | 51 | 31 | 30 | 17 | 27 | 34 | 0.00064 |
| 3280 | 43 | 31 | 30 | 12 | 26 | 31 | 0.00126 |

Over all 288 assessable conditions: LB₅(Q1) > 0 in 164; LB₅(Q2) ≥ 0 in 193; both in 164; UB₉₇.₅(Q2) < 0 (harm beyond
the 10 % margin, descriptively) in 78. A2 improves on A0 on average at every n_train and is worse in roughly a quarter
to a third of assessable conditions at every n_train. The pattern at D2 — A0 already near-optimal, regularisation
costing a little — recurs across the grid and is why H2 was the binding hypothesis.

**Other arms vs A0** (assessable conditions; LB₅ > 0 means the arm is better than A0, UB₉₇.₅ < 0 worse): A2-S 189 / 81
of 290; A2′ 176 / 91 of 288; A4 186 / 75 of 292; A2 β=0.5 113 / 152 of 293; A2 β=1 68 / 195 of 290. The β arms are
harmful: A2 β=1 has higher mean regret than A0 at every n_train, A2 β=0.5 at n_train ≥ 410. A3 equals A0 at
n_train ≤ 820 (floor binding, as the protocol predicted) and differs at 1640 (55 assessable: 33 better, 13 worse; mean
+0.00027) and 3280 (44 assessable: 32 better, 11 worse; mean +0.00071).

**Grid maximum of plain A2 gain** (upward-biased by selection): condition 123 — n_train 410, φ 0.8, Γ 0.20, A 0.15,
heterogeneous — plain mean 0.028376; A0 P(ĵ = j\*) 0.339, A2 0.978.

## 5. Sparse contrasts — 72 of 360, in both directions

An arm is sparse at a condition when fewer than 10 of 1,000 replicates have r > 0 (protocol §10); a contrast with a
sparse arm is not assessable by bootstrap, and the Clopper–Pearson bound is reported instead. The 72 non-assessable
A2 contrasts decompose as:

| | count | by n_train | note |
|---|---|---|---|
| A0 and A2 both sparse | 41 | 820: 10, 1640: 13, 3280: 18 | both arms near-perfect; nothing to separate by regret |
| A2 sparse, A0 not | 17 | 820: 6, 1640: 5, 3280: 6 | A0 has 10–625 replicates with r > 0 while A2 has < 10: A2 better, not bootstrappable |
| A0 sparse, A2 not | 14 | 410: 4, 820: 2, 1640: 3, 3280: 5 | A2 has 17–444 replicates with r > 0 while A0 has < 10: A2 worse, not bootstrappable |

So 55 conditions have sparse A0 regret and 58 sparse A2 regret, with 41 in common. **Sparsity does not mean there is
nothing to improve**: 17 non-assessable conditions are ones where A2 removed essentially all of A0's misselection, and
14 are ones where it introduced misselection that A0 did not have. Sparse counts per arm over the grid: A0 55, A2 58,
A2-S 57, A2′ 58, A3 57, A4 56, A2 β=0.5 53, A2 β=1 48.

## 6. Secondary contrasts (descriptive)

Form effect mean(r_A2 − r_A2-S): 301 assessable; A2-S better (LB₅ > 0) in 122, A2 better (UB₉₇.₅ < 0) in 128; mean
−0.000111. The residual and Schur conditional forms differ at the condition level in both directions and are
indistinguishable on average. Tuning effect mean(r_A2 − r_A2′): 302 assessable; 48 vs 46; mean +0.0000027 — the
finite-n GCV penalty and the infinite-data population-benchmark penalty select essentially the same candidates.

## 7. Diagnostics

Endpoint-hit rate 0, complete-shrinkage rate 0, A4 exact-top-tie rate 0 at every n_train. GCV λ_c medians decrease
with n (condition medians ≈ 0.04–0.13 at n 102; ≈ 0.0013–0.004 at 3280). Ledoit–Wolf δ medians 0.171, 0.051, 0.026,
0.013, 0.007; κ(S_OO) medians 307, 197, 181, 172, 167 across n = 102 … 3280. A1 ≡ A0 held on every replicate: maximum
relative discrepancy over the grid C 3.7e-15, v 1.4e-15, s 3.8e-15.

## 8. What this result does and does not say

- Primary outcome is **exact population regret of the selection** under the frozen Gaussian family (protocol §3). It
  does not measure a fitted downstream predictor's risk, and no claim about predictive performance follows.
- Gaussian-only, synthetic, 60-dimensional, five sample sizes, one calibration family: no claim about nonlinear
  acquisition, about NOVA-S (not retuned), or about any external cohort.
- The decision in §1 is about the two designated conditions chosen by frozen population criteria. The descriptive
  grid pattern — average improvement with a substantial minority of harms, concentrated where A0 is already
  near-optimal — is a hypothesis for a future prospective design, not a finding of this one.
- The frozen G1 rule (A0) stays the programme's comparator. A2 is not adopted. Nothing here reopens NOVA-D or the
  closed archive.

## 9. Software validation (separate from scientific evidence)

Before the freeze, the pipeline was exercised on **fabricated** replicate logs in throwaway copies: a rejection
battery (15 type/domain + 2 non-finite + 3 structural cases; stale-PASS replacement on malformed input and on two
crash paths; seven analyze refusals after log, freeze, report or checker drift) on the owner's machine and in the
drafting container; a complete degenerate log exercising the all-sparse branch; and a complete non-degenerate log on
which `analyze`'s positive path ran and was reproduced independently with 0 mismatches. Those fixtures and their
"positive decisions" are software validation only and carry no scientific content. The independent recomputation of
§2 was then applied to the real log.

## 10. Archive

Committed: `g1r_out/replicates.jsonl.gz`, `g1r_out/integrity.json`, `g1r_out/g1r_result.json`, this record, and the
checksum manifest `g1r_out/SHA256SUMS_G1R_v1.0.txt`. The raw `g1r_out/replicates.jsonl` is kept locally (198 MB,
over the hosting file-size limit); its SHA-256 `ea33587e…` is the canonical identity bound by the integrity report, and
the archive is verified to decompress to it. One run per freeze; `freeze` refuses a second record; the archive is
closed.
