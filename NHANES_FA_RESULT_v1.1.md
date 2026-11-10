# NHANES serum fatty-acid acquisition test: result record v1.1 (9 Nov 2026)

This draft wording correction follows v1.0. The saved development output is the source for the correction.
Its SHA-256 is `f7f9f5cad1c25edcd9e395ae4198b5997880d80c58552fee4aa77b91121c2101`. All numerical endpoints and the feasibility decision are unchanged.

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

| Block | Candidate chosen by G1 and NOVA-S in all five folds | NOVA-S departures |
|---|---|---|
| 1 (v0.3) | LBXTSA, tricosanoic acid (23:0) | 0/5 |
| 2 | LBXEPA, eicosapentaenoic acid (20:5n-3) | 0/5 |
| 3 | LBXDTA, docosatetraenoic acid (22:4n-6) | 0/5 |
| 4 | LBXEPA | 0/5 |
| 5 | LBXDTA | 0/5 |

In block 3, the first outer fold (zero-based fold index 0, `inference.picks[10]`), NOVA-V and NOVA-U
chose LBXMR1 (candidate index 0). G1 and NOVA-S chose LBXDTA (candidate index 6). Across the five blocks,
NOVA-S matched G1 in 25 of 25 block-folds. NOVA-V and NOVA-U each matched G1 in 24 of 25.

Pooled over the five blocks, NOVA-S departures were **0/25**. Δ_dev = 0 exactly; SE = 0; MDE undefined (0). Baseline improvement
Δ_G1→∅ = **15.8 %** (design-based one-sided 97.5 % lower bound 13.1 %; v0.3 block alone 18.0 %, lower bound 12.1 %).

Gate v0.4 §6.4: (i) n ≥ 1,500 — met; (ii) ≥ 1 departure — **not met**; (iii) SE finite and > 0 — **not met**;
(iv) MDE ≤ 0.30 pp — not assessable. **Decision: NO-GO.** Δ_dev played no part.

The v0.3 gate had returned GO on the same zero contrast because it treated SE = 0 as arbitrarily high precision;
that defect is the reason for the v0.4 amendment and is disclosed in the protocol's declaration.

## 3. What this establishes

1. **Selector agreement in the development cycle.** G1 and NOVA-S chose the same candidate in all 25
   block-folds. NOVA-V and NOVA-U each differed from G1 once. This is a record of the fitted choices,
   and it leaves the population risk contrast unresolved.
2. **Concentration of acquisition value.** The stable G1 picks were chain-length or series neighbours of
   the targets in their blocks (23:0 for 22:0/24:0 targets, and EPA and DTA for the n-3/n-6 series).
   The recorded 13–18 % gain over the base-panel loss shows that the selected measurement had value
   in development. The stability is consistent with a clear leader in these blocks.
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
  NO-GO; G1 and NOVA-S agree in 25 of 25 block-folds."

## 5. Record disposition

Preserve v1.0, the frozen development JSON and the run log. This revision changes the narrative record only.
