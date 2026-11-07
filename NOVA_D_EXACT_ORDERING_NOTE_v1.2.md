# NOVA-D: exact population ordering of the frozen G1 rule — note v1.2 (7 Nov 2026)

**Correction log (v1.1 → v1.2), from the owner's review of 7 Nov 2026; four precision items, no result changes.**
(a) Outer-fold training sizes are not a single int(0.8 n): with five folds they are 102/103 (n = 128), 409/410
(512), 819/820 (1024), giving α ∈ {0.3704, 0.3681}, {0.1279, 0.1277}, {0.0683, 0.0682}. Verified exactly: the
population argmax and the top-5 are identical at both training sizes in all 72 conditions, so every reported
ordering stands; the tables below quote the larger-α value. (b) The floor binds when d/(d + n_train) ≤ 0.05,
i.e. **n_train ≥ 1140**, not total n > 1140. (c) The primary's ordering at the floor is a **near-tie** with a
nonzero gap (0.7821196 vs 0.7820558, 8.16 × 10⁻⁵), not a tie. (d) Condition 3's G1 made 17 distinct picks at
n = 512; the four candidates named are its four **most frequent** picks.

**Supersedes v1.1 (`ac37d063…`, committed at `87a513a`) and v1.0 (`c1b31969…`, `3d13ebd`); both retained.**

Read-only diagnostic artefact accompanying
`NOVA_D_RESULT_v1.2.md` (owner's record, `60056edc…`). Not an experiment; changes no frozen output, protocol or
preregistration; reruns no simulation. Prepared by Claude (Fable 5.1); author of record Anas Enoch.

**Correction log (v1.0 → v1.1), from the owner's review of 7 Nov 2026.**
1. v1.0 evaluated the frozen rule only at its shrinkage floor α = 0.05 and called that "its large-n limit". The
   floor is the n → ∞ asymptote; the rule actually sets α_n = min(1, max(0.05, d/(d + n_train))) with d = 60 and
   n_train the outer-fold training size — 102/103, 409/410, 819/820 in the NOVA-D grid, giving **α ≈ 0.370
   (n = 128), 0.128 (512), 0.068 (1024)** (both fold sizes yield identical orderings; §(a) above) — and the floor
   binds only for n_train ≥ 1140. The population object that finite-sample picks at a given n scatter around is
   the rule applied to the exact covariance at α_n. v1.1 reports that ordering for every condition.
2. v1.0 said condition 3's candidate 10 "picks at n = 512 are finite-sample noise, not a population reordering".
   **Wrong.** At α = 0.128 the exact population top-5 of condition 3 is [0, 3, 5, **10**, 7] with scores 0.1991,
   0.1953, 0.1948, 0.1948, 0.1848: candidate 10 is fourth, within 2.2 % of the leader. G1's four most frequent
   picks at n = 512 (3:103, 10:96, 0:86, 5:86; 17 distinct candidates in all) are the population top four at the α
   in force, nearly tied.
3. v1.0 said the positive control's limiting pick being j\* "is why" the finite-sample G1 reached j\* in 67 % of
   folds. Too causal. The correct statement: at α = 0.068 the rule's population argmax is j\* with a 4.1 %
   top-two gap, and a 67 % hit rate is **consistent with** sampling variation around that ordering.
4. v1.0 said the primary (31) "does not flip". True at the floor (near-tie, gap 8.2 × 10⁻⁵); **false at its operative
   α = 0.128**, where the population argmax is j\* with a 3.2 % gap. Both designated conditions therefore had the
   frozen comparator's population pick equal to the nonlinear optimum at the α in force.
5. v1.0 said no condition with Γ_rel = 0.20 flips. True at the floor; false at α = 0.370 (three M4 conditions at
   n = 128). v1.0 said every M2 flip lands on j\*. True at the floor; false at α_n (four M2 flips land on a low-G
   candidate).

**Method.** `nova_d_exact_ordering_alpha_n.py` (companion to `nova_d_exact_ordering.py`, which it imports for the
exact covariance builder and the frozen-rule algebra with the 1e-6 ridge). For each of the 72 conditions it
evaluates the frozen rule on the exact population covariance at α_n and at the floor, with no sampling.

**Results (exact, all 72 conditions).**

| Quantity | At the operative α_n | At the floor α = 0.05 (v1.0) |
|---|---|---|
| Frozen rule's population argmax ≠ Gaussian argmax | **27 of 72** | 15 of 72 |
| — by Gaussian top-two gap | 24 at Γ_rel = 0.02; **3 at Γ_rel = 0.20** (M4 54, 56, 64; all n = 128, α = 0.370) | 15 at 0.02; 0 at 0.20 |
| — by n | 11 (n = 128), 9 (512), 7 (1024) | — |
| M1 (Γ_rel = 0.20, H ≡ 0) | 0 flips at any n | 0 |
| M3 (Γ_rel = 0.02, H ≡ 0) | 9–14 flip at n ≤ 512; none at n = 1024 | 11–14 |
| M2 flips | 18; **14 land on j\***, 4 do not (18, 27, 36 → candidate 7; 39 → candidate 3; all A = 0.05) | 11; all on j\* |
| Primary (31), α = 0.128 | population argmax **= j\***, gap 3.2 % | near-tie (gap 8.2 × 10⁻⁵, nonzero) |
| Positive control (43), α = 0.068 | population argmax **= j\***, gap 4.1 % | = j\* |
| Condition 3, α = 0.128 | top-5 [0, 3, 5, 10, 7], scores 0.1991 / 0.1953 / 0.1948 / 0.1948 / 0.1848 | top-4 [0, 3, 5, 1] |
| Same rule at α = 0 | population argmax = Gaussian argmax in 72 of 72 | 72 of 72 |

The nine nonlinear conditions with a positive NOVA-S lower bound, at their α_n: the frozen rule's population argmax
is the Gaussian pick in 20, 29, 45, 59; a low-G candidate (7 or 14) in 27, 36, 54; the Gaussian runner-up in 56, 64.
In none of the M4 cells is it the nonlinear optimum (candidate 6).

**What this certifies (Objective A only).** The frozen G1's population ordering is a function of α, hence of n. At
every n in the NOVA-D grid the rule's population argmax differs from the Gaussian argmax in a substantial fraction
of the near-tie conditions (Γ_rel = 0.02), and at n = 128 the heavy shrinkage (α = 0.370) overturns even a 20 %
Gaussian gap in three M4 conditions. The discrepancy vanishes at α = 0 in every condition. In **both** designated
conditions the frozen comparator's population pick at the operative α was the nonlinear optimum; the finite-sample
hit rates (60 % and 67 %) are consistent with sampling around those orderings. The realised opportunity for NOVA-S
*against the frozen comparator* was therefore diminished in the primary as well as in the positive control — a
further qualification of how much each tested, which does not alter the frozen decisions (κ_S = 0, zero
departures), the null controls, or the finding that NOVA-V identified j\* less often than G1.

**What this does not certify.** Finite-sample decision loss (Objective B); any causal attribution of the
estimated G1's excess regret to the shrinkage rule; anything about the archive's 71 % self-agreement. Those remain
for the prospective G1 study (result v1.2 §8.3). One design consequence for that study follows directly: its
self-test must compute each arm's population argmax **at every actual training size and regularisation setting
in the grid** — fold by fold where fold sizes differ — not at an asymptote.

Output: `nova_d_out/exact_ordering_alpha_n.json`, computed on the owner's machine
(`e649836972bb76e739584db34b3f09775cae157217e64b040071862cd42fa8b5`, committed at `87a513a`). It evaluates the
larger-α fold size; the smaller-α fold size gives identical argmax and top-5 in all 72 conditions (verified, §(a)).
