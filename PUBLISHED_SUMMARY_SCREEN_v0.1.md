# Published-summary certificate screen: retrospective planning only

8 November 2026. This addendum explicitly uses the versioned result record's
published summaries, with the user's authorization. Earlier population-only
notes remain unchanged. No raw replicate log, archive result JSON or new
baseline samples are used. No original record, implementation, preregistration
or frozen output is changed.

## 1. Provenance and decision on reuse

Read-only use of these summaries is appropriate for a clearly labelled post-hoc
planning screen. It does not provide an independent validation result or turn
sample averages into exact population moments.

Source: [G1R_RESULT_v1.1.md, section 3](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/G1R_RESULT_v1.1.md),
SHA-256 `3855d57ad8c8038c67d0c5160784ca13d07acdff23a00a1a2171e1d8aa8562d2`.
At D1 the record gives A0 mean regret .002574, optimal-pick frequency .178,
and 822 nonzero regrets out of 1,000. Its reported A2 mean improvement is
.001242. The bound below uses only the two baseline summaries and the declared
population gain profile; the A2 number is descriptive context.

No fresh baseline characterization was run. The proposed conclusion that a fresh
finite sample's answer "will not move" does not follow from the regret identity.

## 2. Correct the tail decomposition

Let q=1-p_0, a=f_01=.003, and
b=f_02=.0388368820321, the smallest gap among candidates g>=2. Write

    m=sum_(g>=2) p_g,
    T=sum_(g>=2) p_g f_0g,
    R=a(q-m)+T.

Then the excess above the all-runner-up reference is

    E=R-aq=sum_(g>=2) p_g(f_0g-a)=T-a m.

Therefore

    m<=E/(b-a),
    T=E+a m<=[b/(b-a)]E.

The denominator in the proposed tail-probability bound must be b-a, not b.
The excess E is not itself the tail's entire regret contribution.

Using the displayed numbers as plug-in moments gives

    qhat=.822,
    Ehat=.000108,
    m upper=.003013655,
    T upper=.000117041,
    runner-up probability lower=.818986345.

Thus the basic picture of overwhelmingly runner-up errors is supported by the
record, with this correction to the decomposition.

## 3. A simple aggregate bound now covers all beneficial destinations

Assume a uniformly valid positive-gain certificate, such as departure only when
LB(f_jg)>0, in the fixed-W Gaussian model allowing unknown candidate variances.
On a runner-up baseline pick, the only beneficial destination is arm 0. Its
probability is bounded by beta_01. On any g>=2 baseline pick, even an oracle
can recover no more than f_0g. Harmful departures can be discarded when taking
an upper bound. Consequently every such policy satisfies

    expected improvement <= a beta_01+T
      <= a beta_01+[b/(b-a)] [R-a(1-p_0)].

This argument does not assume that tail departures go to arm 0 and does not
multiply dependent selection/certification probabilities.

If R and p_0 were exact population moments, this would be a population bound.
Plugging in their recorded estimates instead gives a *plug-in value of the
bound*. With beta_01=.0741555 (nu=389), it is .000339508. With the full-data
oracle beta_01=.0749024 (n=410 under the stated null extension), it is .000341748.

An independent all-pair linear program, maximizing over reference distributions
with those two moments, reproduces the formula to numerical precision. The
relaxation is sharp at those moments: put all tail mass on candidate 2, whose
power cap is much larger than that mass. This is sharpness of the *relaxation*,
not construction of an attainable uniformly valid policy.

These numerical ceilings retain their experiment/model qualifications. A policy
using only residual-scatter fixed-pair tests has the nu-based cap; a general
full-data policy needs its full-data comparison. Known conditional variances,
positive permitted harm in the certificate, or an expected-cap-only policy
would change the relevant bound.

## 4. Rounding permits a stronger empirical reconstruction

Under usual nearest rounding of the published mean to six decimals,

    .0025735 <= empirical mean regret < .0025745.

Every g>=2 draw contributes at least (b-a)/1000 excess above a runner-up draw,
so there can be at most three such draws. Enumerating their possible gains
under the declared profile leaves a unique compatible empirical histogram:

| Candidate, zero-based | Inferred count |
|---|---:|
| 0 | 178 |
| 1 | 819 |
| 2 | 3 |
| 3 through 19 | 0 |

This implies an empirical mean .0025735106461, which rounds to the displayed
.002574. This histogram is **inferred from the summaries and the rounding
assumption**, not independently checked against raw outcomes. If the displayed
mean were merely an approximation with a different error tolerance, the unique
histogram inference would no longer follow.

Evaluating the earlier min-bound/all-pair relaxation at these inferred empirical
frequencies gives:

| Information convention | Plug-in value | Fraction of inferred empirical A0 mean regret |
|---|---:|---:|
| Residual-scatter certificate, nu=389 | .000338977 | 13.1718% |
| Full-data oracle extension, n=410 | .000341218 | 13.2588% |

Thus **about .00034, or 13.2%-13.3%**, is the appropriate retrospective planning
summary, rather than a strict .00033 or <=13% population assertion.

Keeping only the mean-rounding interval, without the discrete histogram
reconstruction, gives plug-in ranges [.000338966,.000340049] for nu=389 and
[.000341206,.000342290] for n=410. The practical rounding effect is small.

## 5. Sampling uncertainty is different from rounding uncertainty

The identities hold for true population probabilities and also for empirical
frequencies, but substituting the latter into a population probability bound
does not preserve its guarantee. The NP caps are population bounds; they are
not caps on arbitrary finite-sample empirical rejection counts.

The inferred three tail picks out of 1,000 are particularly uncertain. As one
simple pointwise illustration, their one-sided 95% Clopper-Pearson upper bound
is p_tail<=.00773524. A conservative uncertainty-aware bound is then

    expected improvement <= a beta_01+f_0,19*p_tail_upper,

about .001375 for the scatter-based convention. This is deliberately loose,
not a sharp bound, not a new programme-level confidence statement and not a
formal comparison with A2. It demonstrates why the .00034 plug-in value cannot
be called a guaranteed population ceiling from these summaries alone.

No uncertainty calibration has been optimized after seeing this result. A fresh
baseline sample would also fluctuate, especially in these rare tail picks.

## 6. What can be recorded as programme rationale

The corrected retrospective screen is near .00034, versus a reported A2
improvement .001242 (about 48% of the reported baseline mean). That is useful
descriptive planning evidence for deprioritizing a standalone positive-gain-
certificate development line. It is not a formal population dominance result:
the prospective harm constraint, other conditions and sampling uncertainty
remain part of the comparison, and no candidate certificate is implemented.

The plug-in screen is still about 1.32 times the 10% planning margin. It does
not exclude that target; it shows limited headroom in this optimistic screen.

Record the decision as a programme priority judgment, **not an impossibility
theorem or a newly established 13% population ceiling**. This use of the existing
result does not reopen its prespecified claims, tune a successor on the raw
archive or establish an independently tested challenger.

The common umbrella stays unchanged: average regret under a declared Pi with an
absolute expected-harm cap over a declared Theta_cap; normalization, epsilon
and tau are externally justified before future outcomes. A0 remains the
comparator; no challenger is adopted. H2 was not established for A2 and the
frozen joint criterion remains unmet.

Reproducibility: `published_summary_screen.py` reads the hash-verified Markdown
record and earlier population-arithmetic JSON, and writes
`published_summary_screen.json`. It performs no sampling and reads no raw run
or archived result JSON. The new note is separate from the frozen repository.
