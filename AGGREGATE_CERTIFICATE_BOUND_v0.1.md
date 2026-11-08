# Aggregating certificate power with a data-selected baseline

8 November 2026. Mathematical audit and deterministic population arithmetic.
Earlier notes remain unchanged. No archived outcomes are read for this analysis;
no A0 sampling, experiment rerun, Monte Carlo, coverage study or new selector.

## 1. The proposed product is not an upper bound

Let P be the fixed true population, g(D) the baseline's data-dependent pick,
h(D) the new pick, and

    B_g={g(D)=g}, p_g=P(B_g),
    E_0g={g(D)=g and a certificate authorizes departure to true-best arm 0},
    Delta_0g=(G_0-G_g)/R0.

For a uniformly valid simultaneous positive-gain certificate, the event E_0g
has probability at most alpha at the constructed boundary null for pair 0,g.
Therefore its probability at P is at most the pair's NP ceiling beta_0g.
It is also at most p_g, because E_0g is a subset of B_g. Thus

    P(E_0g)<=min(p_g,beta_0g),
    expected direct-to-best improvement
        <= sum_(g!=0) Delta_0g min(p_g,beta_0g).

In contrast, the proposed sum of p_g Delta_0g beta_0g assumes

    P(certificate | B_g)<=beta_0g.

The NP ceiling obtained earlier is unconditional. It supplies no such
conditional statement. From it one gets at most

    P(certificate | B_g)<=min(1,beta_0g/p_g), p_g>0.

This is an actual selection-dependence issue, not an error allocation convention.

### A concrete check of the failure

In the known simple-pair testing experiment take E to be the most-powerful
size-alpha rejection event and define a hypothetical baseline to pick g on an
event B containing E, with P(B)=.10, and arm 0 otherwise. At the D1-like
residual-scatter input P(E)=beta=.0741555. Then

    P(E and B)=.0741555,
    P(B)*beta=.00741555.

The top-two normalized gap is .003. The actual joint-event contribution can be
.000222467, while the proposed product is .0000222467. For this hypothetical
reference, mean regret is .0003. The invalid product would even fall below its
10% margin .00003 and falsely suggest exclusion. No assertion is made that A0
has this hypothetical selection event. It proves why a product bound requires
additional conditional analysis.

For the actual frozen A0, the dependence must be accounted for, not assumed away.

## 2. A full policy needs all beneficial destinations, not only arm 0

The preceding sum caps departures specifically to the true best arm. A general
policy may depart from g to some other arm j with G_j>G_g. Replacing the attained
gain by G_0-G_g does not make the event a test of G_0>G_g: under the boundary
null equalizing arms 0,g, another arm can legitimately have gain above both.

For the general positive-certificate policy let

    E_jg={baseline g, new pick j, positive-gain certificate},
    x_jg=P(E_jg), Delta_jg=(G_j-G_g)/R0.

For every pair with G_j>G_g,

    0<=x_jg<=beta_jg,
    sum_(j:G_j>G_g) x_jg<=p_g.

Ignoring all harmful departures gives the valid optimistic bound

    expected improvement <= U(p),
    U(p)=max_x sum_(g,j:G_j>G_g) Delta_jg*x_jg
         subject to those constraints.

This is a small linear-program relaxation. Each baseline row can be solved by
filling the highest-gain destinations first, subject to their individual power
caps. The individual ceilings need not be jointly attainable; simultaneous
coverage can impose additional constraints. Consequently U(p) above the desired
margin is inconclusive, not evidence that a concrete candidate achieves it.

A sharper selected-event bound is

    beta_jg^B(alpha)=sup_{E subset B_g, P_null,jg(E)<=alpha} P(E).

Its optimizer ranks likelihood ratios *inside B_g*. It requires the joint law
of A0's selection and the relevant likelihood statistic. It is not p_g times
the unselected NP power. Conditional coverage, if desired, would need the
selection-conditioned null and its own size statement.

These bounds assume certification of strictly positive gain, such as LB(f)>0.
A certificate allowing harm tau, or a policy required only to meet an expected-
harm cap, has a different null/constraint and is not automatically covered.

## 3. Population arithmetic completed without archive outcomes

For the D1-like profile Gamma=.02,A=.15, the calibration fixes

    r=.7562116868563323,
    G_0=1, G_1=.98, G_k=.98*r^(k-1) for k>=2,
    R0=1/A.

No loading seed is needed for the pairwise variance-null ceilings: each depends
only on G_g/G_j, alpha and the observation count. The script evaluates all 190
beneficial ordered-pair caps, using the repeated geometric ratios efficiently.

For departures specifically to arm 0, with alpha=.05:

| Baseline candidate, zero-based | Normalized gap | Ceiling with nu=389 | Ceiling with full n=410 |
|---:|---:|---:|---:|
| 1, runner-up | .003000 | .0741555 | .0749024 |
| 2 | .0388369 | .905437 | .918119 |
| 3 | .0659372 | .9999846 | .9999923 |

The all-pair calculation is in `aggregate_certificate_screen.json`. Independent
normal-mixture and conditional-chi-square-CDF integration representations agree
within 1.33e-13. The numerical decimals are approximations, not formally
interval-certified bounds. Hypothetical-reference greedy solutions were checked
against an independent linear-program solver.

The computation deliberately leaves

    actual_A0_pick_distribution = null,
    actual_A0_aggregate_bound = null.

Maximizing the relaxation over *all possible* reference distributions gives a
very loose bound about .14904. That value describes a relaxation allowing a
reference concentrated on a very poor arm; it is not A0's gain, a selector's
achievable performance, or a scientific result. This reference-free screen
does not exclude the direction.

### Why the existing exact-ordering tools do not supply p_g

[The population score routine](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/nova_d_exact_ordering.py:40)
applies the rule to an exact Sigma. Its argmax tables are population orderings,
not finite-sample pick probabilities. The new Wishart pair-law tool concerns
the unshrunk OLS residual rule; it also does not supply A0's joint 20-arm picks.

For Gaussian full rows, n-1 times the centered sample covariance has a
Wishart_(n-1)(Sigma) law. A0 standardizes that random matrix and then applies
joint shrinkage and the frozen ridge/score rule. In principle,

    p_g=integral 1{A0(Corr(V))=g} d Wishart_(n-1)(Sigma)(V).

This is indeed a population quantity. Evaluating this different, high-dimensional
integral is not already implemented by the exact-ordering tools. No frozen-run
pick frequencies have been inserted as exact probabilities, and no new sampling
of the frozen condition was performed under the name of a population calculation.

If the certificate's fixed-pair rejection statistics use only the residual
scatter, their nu-based ceilings remain valid after intersecting with A0's
selection event. For a general full-data certificate, use an appropriate
full-data experiment and null extension; the n=410 oracle ceiling is the
conservative comparison supplied here. Selection alone does not license a
change of degrees of freedom.

## 4. Replace the hard variance 'floor' with a power statement

The single-variance relative standard error sqrt(2/nu) is correct for Q_jj/nu.
For two independent residual candidates, however,

    Var(log(Q_jj/nu)-log(Q_gg/nu))=2*trigamma(nu/2),
    log-ratio SE is approximately 2/sqrt(nu).

It is about .10153 at nu=389 and .03504 at nu=3259, rather than the
single-variance values .07170 and .02477. Comparing log(1/(1-Gamma)) with
one log-ratio standard error gives a crossing near nu=1521 for Gamma=.05,
not nu=800. This is a diagnostic scale, not a resolution theorem.

The precise statement already supported by the NP arithmetic is more useful:
for Gamma=.02, alpha=.05, fixed-pair power in the stated unknown-variance model
cannot exceed roughly 6.01%, 7.42%, 8.70%, 10.78%, 14.27% at the five residual
degrees of freedom. Thus 95% power at 5% size is ruled out throughout this grid.
Zero power, complete lack of ranking information, and a causal explanation of
A2's gain do not follow. The bound does not apply to unconditional rules as a
type-I-error restriction they never claimed to satisfy.

Likewise a power ceiling does not identify a particular estimator's standard
error or its D2 false-departure probability. Those remain method-dependent.

## 5. Current decision

The proposed product sum cannot be recorded as a certificate exclusion bound.
For direct-to-best departures the correct elementary bound is the min sum;
for an arbitrary policy use the all-pair relaxation or a stronger joint analysis.
The necessary finite-sample reference probabilities are not supplied by the
current tools. The completed population arithmetic therefore does not establish
exclusion at D1 and does not establish an attainable certified candidate.

No standalone certificate programme is opened. The common prospective umbrella
remains average regret under a declared Pi with an absolute expected-harm cap
over a declared Theta_cap. Its external normalization, epsilon and tau remain
prerequisites. A0 stays the comparator; no challenger is established; H2 was
not established for A2 and the frozen joint criterion remains unmet.

Reproducibility: `aggregate_certificate_screen.py` writes
`aggregate_certificate_screen.json`, using population profile equations only.
It imports no programme code and reads no archived experiment outcomes.
