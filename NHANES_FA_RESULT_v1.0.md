# NHANES serum fatty-acid acquisition test: result record v1.0 (4 Nov 2026)

**Outcome: NO-GO at the development-cycle feasibility gate of protocol v0.4. The 2013–2014 evaluation cycle (FAS_H,
DEMO_H) was never downloaded. NHANES-FA is recorded as assessed and infeasible for discriminating NOVA-S from G1;
no further amendment is made (protocol v0.4 §6.4).**

## 1. Provenance

| Item | Value |
|---|---|
| Protocol | v0.3 (sha256 e7a342f4…, commit e40f2e4) for the first development run; v0.4 (fcc106d8…, commit 2bb4990) for the second, with the owner's declaration of the development-adapted amendment |
| Wrapper | v0.1.1 (10645154…) and v0.2.0 (f85d5235…); selection and prediction by the frozen `p6_n2.py` 681ce852… and `p5_gate3a.py` 4600e2f3…; exact equivalence with the frozen `run_block` tested on synthetic blocks |
| Blocks | v0.3 block 93c26236… (block 1); v0.4 five-block file e3e1d8ac… (salts fixed in the protocol; block 1 unchanged) |
| Data | NHANES 2011–2012 `FAS_G.xpt` 5a8a907c… and `DEMO_G.xpt` eaf0525d… (CDC public); files kept outside the repository; no other NHANES file read |
| Results | v0.3 run: `nhanes_fa_out/FAS_G_development_result.json` 43f5fe89…; v0.4 run: `FAS_G_development_v04_result.json` f7f9f5ca…; logs `nhanes_fa_dev.log`, `nhanes_fa_dev_v04.log` |
| Prior programme contact with NHANES | Stage 16L–16O (11 Oct 2026): 2017–18 urinary phthalates, status audit only; no fatty-acid file had been opened before 4 Nov 2026 |

## 2. Development cycle (NHANES 2011–2012)

- File rows 5,643; eligible (age ≥ 12, WTFAS2YR > 0, all 29 analytes present) **1,725**; 14 strata, 31 PSUs,
  design df 17; one row per SEQN confirmed; five participant-disjoint folds (seed 20261103), identical in every block.
- Target standardisation: frozen N2 convention (cycle-wide SD, 1 % median floor); constants in the result files.

| Block | Candidate chosen by G1, NOVA-V, NOVA-U and NOVA-S in all 5 folds | NOVA-S departures |
|---|---|---|
| 1 (v0.3) | LBXTSA, tricosanoic acid (23:0) | 0/5 |
| 2 | LBXEPA, eicosapentaenoic acid (20:5n-3) | 0/5 |
| 3 | LBXDTA, docosatetraenoic acid (22:4n-6) | 0/5 |
| 4 | LBXEPA | 0/5 |
| 5 | LBXDTA | 0/5 |

Pooled over the five blocks: departures **0/25**; Δ_dev = 0 exactly; SE = 0; MDE undefined (0). Baseline improvement
Δ_G1→∅ = **15.8 %** (design-based one-sided 97.5 % lower bound 13.1 %; v0.3 block alone 18.0 %, lower bound 12.1 %).

Gate v0.4 §6.4: (i) n ≥ 1,500 — met; (ii) ≥ 1 departure — **not met**; (iii) SE finite and > 0 — **not met**;
(iv) MDE ≤ 0.30 pp — not assessable. **Decision: NO-GO.** Δ_dev played no part.

The v0.3 gate had returned GO on the same zero contrast because it treated SE = 0 as arbitrarily high precision;
that defect is the reason for the v0.4 amendment and is disclosed in the protocol's declaration.

## 3. What this establishes

1. **A replication fact, not a selector comparison.** On an independent survey population with established
   participant structure, the frozen Gaussian criterion and the frozen nonlinear value estimator — including its
   unrestricted form NOVA-V, which disagreed with G1 in 64 % of archive block-folds — made identical acquisition
   decisions in 25 of 25 block-folds. Agreement is a legitimate outcome of the frozen rules; it is not evidence
   that the two policies have equal population risk, and it is not evidence of NOVA-S superiority or inferiority.
2. **Why there is no contrast.** In each block one candidate carries acquisition value far above the rest, and that
   candidate is a chain-length or series neighbour of targets in the block (23:0 for 22:0/24:0 targets; EPA and DTA
   for the n-3/n-6 series). Serum total fatty acids form a highly structured, near-compositional family in which
   the covariance criterion already identifies the dominant measurement; nonlinear information has no second
   candidate to promote. The single measurement is worth 13–18 % of the base-panel loss, so the problem is not
   that acquisition has no value here but that its value is concentrated.
3. **Scope.** Serum total fatty acids by GC-ECNI-MS, NHANES fasting subsample aged 12+, 29 analytes, one-step
   budget. No claim about broader metabolomic panels follows in either direction.
4. **Design lesson, recorded against the diagnostic's author (Claude).** A single block gave roughly a 57 % chance
   of zero departures under archive rates before any fatty-acid specificity; the v0.3 gate was degenerate on a
   zero contrast. Both were correctable before any evaluation data were touched, and both are disclosed.

## 4. Consequences

- FAS_H/DEMO_H are not downloaded; the NHANES route is closed for this question under v0.4.
- The external-validation stub (`EXTERNAL_VALIDATION_PREREG_STUB_v0.1.md`) and the 15 December 2026 LIPID/COMETS
  decision are unaffected; NHANES-FA was a restricted-panel supportive test and was never a P4V unit.
- The NOVA archive paper may cite this record as an independent replication of selector agreement in a structured
  panel; it must not be cited as a test of NOVA-S's advantage.
- Matrix v0.3 §1.6 receives one row: "NHANES-FA (public survey, 1,725 participants, 5 blocks): 0/25 departures;
  NO-GO; selectors coincide where one candidate dominates."

## 5. Files to commit

`NHANES_FA_RESULT_v1.0.md`, `nhanes_fa_out/FAS_G_development_v04_result.json`, `nhanes_fa_dev_v04.log`.
Per-participant tables (`*_participants.csv`) stay local. The NHANES data files stay outside the repository.
