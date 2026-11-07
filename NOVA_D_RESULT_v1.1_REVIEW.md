# NOVA-D v1.1 verification and remaining errata

Review date: 7 November 2026. This note accompanies the exact supplied v1.1;
it does not replace or modify either result record, frozen outputs, or preregistrations.
It records a read-only review, not a new experiment or a new frozen protocol.

## Identity and scope

| File | Verified SHA-256 |
|---|---|
| `NOVA_D_RESULT_v1.0.md` | `992e6433eba4298cb4509693c4b80277e2d79b3a5cd5fd3aba4a00b30dc1bedb` |
| `NOVA_D_RESULT_v1.1.md` | `87d18bcb30c877e5c6f2ae322393bdfdcff10049c3f74e0f3eda5d83d66ede87` |

v1.0 is retained at commit `1f4890265fa1266c2f1dcf972229df7760713740` and on disk.
The comparison confirms that v1.1 changes the interpretation and adds diagnostics;
it does not change the provenance table or frozen prespecified decision values.
The frozen JSON decision flags remain false for primary recoverability and positive-control
viability, with no validity breaches or effective near-tie conditions.

The nine nonlinear conditions with NOVA-S lower bound above zero are
20, 27, 29, 36, 45, 54, 56, 59, 64. Eight have n = 128; condition 59 has n = 512.
The 100 replicates x 5 outer-fold clarification is present. Within-replicate training
sets overlap; the 500 decisions are not 500 independent experiments. The Objective A
(population ranking consistency) / Objective B (finite-sample acquisition loss)
distinction and unchanged external-validation stub and deadline are present.

## Remaining inconsistencies and limits

1. **The proposed matrix row in section 8 item 4 is unchanged from v1.0.** Its
   attribution of all positive n = 128 signals to G1 degradation remains too broad,
   omits condition 59, and does not specify which selector or positivity criterion
   defines a signal. The claimed correction "throughout" is therefore incomplete.
   Only matrix versions through v0.4 are present in the repository; the v0.5 row is
   a proposal, not an already updated matrix.
2. **The twelve-cell list is not exhaustive.** Section 7.3 says three nulls have
   NOVA-V rho > 0, but twelve null conditions have a positive point estimate:
   0, 2, 3, 5, 6, 9, 11, 12, 13, 14, 15, 17. Even requiring a positive lower bound
   leaves four nulls: 3, 11, 12, 13. Condition 13 (M3, n = 512, A = 0.15) has
   rho_V = 0.0042118905 and lower bound 0.0010408192; oracle regrets are
   0.002028 for G1 and 0.0014276738 for NOVA-V. There are also 21 nonlinear
   conditions with a positive NOVA-V lower bound, so "every positive-signal cell"
   cannot describe only the nine NOVA-S cells and three selected nulls.
3. **Condition 3's pick label is wrong.** Its displayed `0:86` should not be
   presented as the third-most-frequent pick: candidates 0 and 5 each have 86 picks.
   The full leading frequencies are 3:103, 10:96, 0:86, 5:86. The other displayed
   pick frequencies, distinct-pick counts, hit rates, kappas and regret values agree
   with the archived records at their displayed precision.
4. **"0-10% in eight" is incorrect.** In the twelve displayed cells, exactly
   five NOVA-V optimum-hit rates are at most 10%: 45, 54, 56, 59, 64.
   Condition 27 is 10.4%, which rounds to 0.10 but exceeds 10%. Even counting
   that rounded entry gives six, not eight. The maximum is correctly 54%.
5. **The regret comparison does not alone establish the claimed mechanism.**
   Lower NOVA-V oracle regret establishes a decision advantage relative to the
   estimated G1 in those cells; it does not isolate covariance shrinkage as its
   cause or exclude a nonlinear contribution. The null examples establish that
   an advantage over estimated G1 need not be nonlinear recovery. Condition 59
   supports the narrower rank-1/rank-2 correction account: G1 chooses 0 or 1,
   neither selector chooses the nonlinear optimum, and NOVA-S regret decreases
   from 0.03192 to 0.03144. A matched prospective comparison is needed for a
   causal claim about regularization. The 200,000-row diagnostic is finite-sample
   evidence; an exact population calculation would be needed to certify its
   stated population ordering independently.
6. **Other retained wording needs qualification.** Section 5's "none ... at
   n >= 512" is not literally zero measured relative gain: condition 59 has
   kappa_S = 0.0298954. It should distinguish measured gain over estimated G1
   from recovery against the population comparator. Section 7.1's assertion that
   any threshold above approximately 1 yields no departures is too loose given
   the reported maxima 1.03 and 1.19; zero departures in both designated cells
   is guaranteed only above their maximum (1.19 at reported precision).
   The corresponding "signal-to-noise below 1" summary also overstates the
   reported maximum-statistic evidence. Section 7.4 attributes the archive's
   71% self-agreement partly to shrinkage without a matched ablation; this remains
   a hypothesis. Its reference to matrix section 3.3 should be section 3 item 3
   (the count is also in section 1.5 of matrix v0.4). Section 9 lists outputs
   already committed with v1.0; this correction commit adds documentation only.

## Corrected matrix-row proposal

NOVA-D (synthetic, 72 conditions x 100 replicates, frozen code): primary and
positive-control kappa_S = 0, with 0/500 departures each (100 replicates x 5
dependent outer folds); positive-control max t = 1.03 at n = 1024 and designed
population advantage 40%. Near-tie screening has no qualifying improvement and
hurts in 18 condition-threshold pairs; NOVA-S passes all 18 null controls. Nine
nonlinear conditions have a positive NOVA-S lower bound: eight at n = 128 and
condition 59 at n = 512. The four kappa_V > 1 values are all at n = 128 and
cannot be interpreted as recovery against the population comparator. Estimated
G1 has greater oracle regret than NOVA-V in the twelve displayed diagnostic
cells; that list is not exhaustive. Post-hoc large-sample diagnostics indicate
shrinkage-induced ranking distortion; finite-sample superiority of an alternative
regularization rule remains untested.

## G1 prospective study decisions confirmed by the owner

- Primary approach: **decoupled conditional regularization**, with acquisition
  decision loss as the primary outcome. Superiority is a hypothesis to test,
  not an unconditional claim that conditioning first dominates at every n.
- Gaussian-only scope. Keep the original G1 as a comparator; vanishing-floor
  and data-adaptive shrinkage arms are controls/secondary comparisons.
- Grid axes include dimension ratio and the observed-block covariance condition
  number. Define the dimension-ratio numerator, effective training denominator,
  and population versus estimated condition number explicitly before freezing.
- Before freezing, derive and numerically test exact joint-versus-matched-conditional
  shrinkage equivalence, including all scale factors and ridge terms. A conditional
  arm algebraically equivalent to joint shrinkage is an equivalence control;
  the primary decoupled arm must specify the independently regularized components
  that distinguish it. Matching only a scalar shrinkage coefficient is insufficient.
- No reopening the closed archive for selector development, no NOVA-D reruns,
  no change to frozen outputs or existing preregistrations. These decisions do
  not constitute a new frozen protocol or an executed equivalence experiment.
