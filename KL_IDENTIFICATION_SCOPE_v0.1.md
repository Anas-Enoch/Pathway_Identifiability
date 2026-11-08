# KL identification, admissible populations and the successor objective

8 November 2026. Population mathematics and deterministic numerical arithmetic.
No experiment, sampling, coverage simulation, candidate tuning, new selector,
freeze, or alteration of closed outputs. Earlier notes are retained unchanged.

## 1. Record the KL identity, with its experiment and parameter class

For residual candidates Z~N(0,I), target model T=F Z+epsilon with
epsilon~N(0,sigma^2 I), define F_alt by exchanging the magnitudes of the first
two columns while preserving their unit directions u_1,u_2. The preceding
audit establishes

    tr(K_alt^-1 K)-k=||F-F_alt||_F^2/sigma^2,
    det K_alt=det K,
    kappa := KL(N(0,K) || N(0,K_alt))
          =||F-F_alt||_F^2/(2sigma^2)
          =(1-sqrt(1-Gamma))^2/sigma^2
          =80 A(1-sqrt(1-Gamma))^2.

KL is symmetric for this particular pair. Its conditional regression
interpretation gives the same formula directly. The two populations have equal
R0 and opposite normalized top-two gaps, +/-Delta, Delta=Gamma*A.

For s independent residual Gaussian rows, total KL is s*kappa. For the
residual scatter Q~Wishart_nu(K), total KL between the two scatter laws is
nu*kappa. Full scatter is sufficient for this simple zero-mean Gaussian
covariance testing problem.

An n*kappa statement for the *original full data* needs a specified joint
extension: keep the law of O and residual-candidate construction fixed, and
keep the common O coefficient in T fixed while changing F. This extension
indeed has n*kappa KL. It need not preserve every additional marginal-variance
constraint in the literal frozen generator. Conversely nu*kappa is the right
experiment for a rule using only the OLS residual scatter; it must not be
substituted for n*kappa in an unrestricted full-data lower bound.

### The labelled-family qualification is essential

The literal G1-R generator assigns

    G_j=CAL.G_at(j+1,r,gamma), j=0,...,19,

and thus always has jstar=0. See [g1r_run.py:131](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/g1r_run.py:131)
and [its returned optimum](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/g1r_run.py:144).
On this literal labelled class, the constant rule "choose candidate 0" has zero
identification error and regret. Therefore a positive lower bound for *every
rule on the literal frozen class* is false.

The constructive alternative was admissible for the unrestricted SPD confidence
region in the preceding audit. That alone does not make it a member of the
literal frozen labelled family. For an information lower bound, explicitly
define a model containing both populations: for example, a conditional Gaussian
family allowing these loading directions and permutations of the gain ranks.
Any symmetry restriction imposed on the rules is another assumption requiring
its own argument. Neither model enlargement nor a restriction to symmetric
rules may be implicit.

This qualification does not change the frozen comparison of the generic
algorithms. It limits the scope of the new identification theorem. A prospective
generic-selector study should declare candidate-label randomization or the
relevant symmetry explicitly.

## 2. The Bretagnolle-Huber table is a pairwise minimax/average bound

For a binary rule distinguishing two admissible populations P_0,P_1, let e_0
and e_1 be its two error probabilities. Bretagnolle-Huber gives

    (e_0+e_1)/2 >= (1/4) exp(-s*kappa),
    max(e_0,e_1) >= (1/4) exp(-s*kappa).

It does not lower-bound either individual error in isolation. A constant rule
can have e_0=0 and e_1=1. Thus this table is not a pointwise error floor for D1.

For completeness, let affinity a=integral sqrt(p_0 p_1). Jensen gives
a>=exp(-KL/2), and total variation satisfies TV<=sqrt(1-a^2). Therefore optimal
equal-prior error (1-TV)/2 is at least
(1-sqrt(1-exp(-KL)))/2 >= exp(-KL)/4. This proves the factors used here.

Using s=n for the common-nuisance full-data extension, the requested arithmetic is:

| Gamma | A | KL per row | First integer n with BH floor <=.05 | Floor at n=102 / 410 / 820 / 1640 / 3280 |
|---|---|---:|---:|---|
| .02 | .05 | .000404051 | 3984 | .23991 / .21183 / .17949 / .12887 / .06643 |
| .02 | .15 | .001212152 | 1328 | .22092 / .15209 / .09253 / .03424 / .00469 |
| .02 | .20 | .001616203 | 996 | .21200 / .12887 / .06643 / .01765 / .00125 |
| .05 | .15 | .007693572 | 210 | .11406 / .01067 / 4.55e-4 / 8.28e-7 / 2.75e-12 |
| .20 | .15 | .133747416 | 13 | 2.97e-7 / 3.83e-25 / 5.86e-49 / 1.37e-96 / 7.53e-192 |

The integer thresholds round upwards from the real-valued crossings. A lower
bound falling below .05 does not establish that .05 error is achievable. For
example, at Gamma=.05,A=.15,n=410 the BH floor is .01067, while the optimal
known-pair error below is about .1048.

The normalized minimax regret over the known pair is bounded below by Delta
times the corresponding testing bound: every nonoptimal arm loses at least
Delta in either population. This is an absolute decision-risk bound, not a
baseline-relative no-harm result.

## 3. The known-pair Bayes/minimax error has an exact low-dimensional law

Let rho=u_1'u_2, with |rho|<=1. Conditional on the observed residual candidates,
the log likelihood ratio ell=log(dP_0/dP_1) has laws

    ell | Z ~ N(D/2,D) under P_0,
    ell | Z ~ N(-D/2,D) under P_1,
    D=sum_i ||(F-F_alt)Z_i||^2/sigma^2.

For s rows,

    D = kappa[(1+rho)U+(1-rho)V],
    U,V independently chi-square_s.

The equal-prior likelihood-ratio rule has equal errors at the two populations,
so its Bayes error is also the minimax binary error:

    R_star=E Phi(-sqrt(D)/2).

Its value depends on rho, not only on total KL. A one-dimensional integral is

    R_star=(1/pi) integral_0^(pi/2)
       [(1+kappa(1+rho)/(4sin(theta)^2))
        (1+kappa(1-rho)/(4sin(theta)^2))]^(-s/2) dtheta.

No numerical integration is needed for tight angle-uniform bounds. The function
h(x)=Phi(-sqrt(x)/2) is convex for x>0. Its expectation under the weighted sum
of iid U,V is an even convex function of rho, hence increases with |rho|.
At rho=0 and |rho|=1 the expectation is a Student-t tail. Consequently

    Pr(t_(2s)>sqrt(s*kappa/2)) <= R_star
        <= Pr(t_s>sqrt(s*kappa/2)).

These are exact mathematical bounds, evaluated numerically by standard Student-t
CDF functions. They give:

| Population profiles | Experiment | Bounds on optimal known-pair error |
|---|---|---|
| D1-like: Gamma=.02,A=.15 | n=410 full rows | [.3091366,.3092034] |
| D1-like | nu=389 residual scatter | [.3137097,.3137781] |
| D2-like: Gamma=.20,A=.15 | n=410 full rows | [1.0421e-7,1.3121e-7] |
| D2-like | nu=389 residual scatter | [2.1309e-7,2.6536e-7] |

Thus the claimed contrast in *this testing problem* is supported. The
known-pair minimax normalized regret equals Delta*R_star: restricting to the
two leading arms loses nothing, and the likelihood-ratio rule equalizes their
risks. These are oracle known-hypothesis benchmarks, not risks of A0 or A2 and
not exact error probabilities for selecting among 20 unknown gain vectors.

## 4. Identification error does not imply zero-power certification

A test can abstain/defer to a baseline and control its false departures while
having positive power. The sum-of-errors lower bound does not prohibit this.

Here is a mathematical counterexample inside the *known two-point experiment*.
Take candidate 2 as a fixed reference. It is worse under P_0 and optimal under
P_1. For alpha>0, depart to candidate 1 when

    ell > -D/2+sqrt(D) z_(1-alpha).

This has conditional and unconditional false-departure probability exactly
alpha under P_1. Under P_0 its departure probability is

    E Phi(sqrt(D)-z_(1-alpha)) > 0.

At the orthogonal and collinear direction endpoints it is the noncentral-t CDF
at sqrt(2s*kappa), with noncentrality z_(1-alpha) and degrees of freedom 2s or s,
respectively. For D1-like n=410:

| alpha | Power, rho=0 | Power, |rho|=1 |
|---:|---:|---:|
| .05 | .25850 | .25846 |
| .01 | .09188 | .09190 |
| .001 | .01818 | .01819 |
| .0001 | .003249 | .003254 |

These endpoint values are not asserted to bracket every intermediate rho.
Under P_0 the expected regret improvement over this reference is Delta times
power; under P_1 expected harm is Delta*alpha. This supplies a valid one-sided
certificate over the two known hypotheses: its lower bound is +Delta on the
departure event and -Delta otherwise. It covers under P_0 always and under P_1
with probability 1-alpha.

It is **not a proposed programme selector**: it knows F,F_alt and sigma^2, and
its fixed reference is not A0. It refutes the logical inference that the
identification floor alone makes every certified departure impossible. At
very small failure budgets power can indeed be small; that is a quantitative,
budget-dependent question, not a zero-power theorem. At alpha=0, mutual absolute
continuity of the two Gaussian experiments does rule out a positive-power
zero-false-departure event.

Nor does a lower bound pairing D1 with its near-tie swap connect D1 to the
different wide-gap D2 population. Establishing an unavoidable H1/H2 trade-off
would require a constrained risk or power bound involving the actual risk
requirements and relevant populations. The displayed KL table does not provide
one. It therefore does not prove that H2 was impossible, or the wrong scientific
formalization, or that no rule can meet it while improving at D1. The completed
procedure's unmet joint criterion remains exactly that result.

Population rank reversal and the testing floor are compatible with the proposed
bias/variance interpretation. They do not themselves identify the cause of
A2's finite-sample gain, establish baseline optimality in resolvable regimes, or
show that all helpful departures occur in information-poor regimes.

## 5. A rigorous random-Q extension replaces the suggested substitution

Replacing K by Ktilde in the central witness is not automatic: a generic sample
candidate block is correlated and its target conditional noise is not sigma^2 I.
The pivot's standard deviation alone does not establish inclusion of an
adaptively constructed alternative.

For the fixed true pair, equal determinants give exactly

    T(K_alt;Q)-T(K;Q)=tr[(K_alt^-1-K^-1)Q]=2ell.

For s residual degrees of freedom, its mean and variance under K are

    E[2ell]=2s*kappa,
    Var[2ell]=8s*kappa+4s*kappa^2(1+rho^2).

The roughly 41 standard deviation of the full covariance pivot is a different
random quantity. The trace-functional variance identity independently verifies
the displayed likelihood-ratio variance.

There is also a useful distributional statement requiring no full-pivot cutoff.
If a region C has coverage >=1-delta at each population, total variation yields

    P_0(K_alt in C) >= P_1(K_alt in C)-TV(P_0,P_1)
                       >= 2R_star-delta,
    P_0(K and K_alt in C) >= 2R_star-2delta.

The same statement holds for intervals covering the opposite gain values. For
D1-like residual scatter nu=389 and delta=.05, any such valid region includes
the alternative under P_0 with probability at least .57742, and both populations
with probability at least .52742. These bounds substantiate substantial random
data ambiguity on the enlarged two-point class. They are not probability-one
claims and do not establish the optimal departure power of an A0-based rule.

## 6. A split-sample direction does not preserve the full target gain automatically

For W=I and a fixed unit direction a, the scalar projected gain is
(a'c_j)^2/v_j <= ||c_j||^2/v_j=G_j. A scalar regression pivot can give inference
for that projection. Equality requires alignment with c_j. Fixing the direction
on an independent split solves its selection problem, not the loss in the
unmeasured orthogonal components.

A lower bound on projected challenger gain remains a lower bound on full gain;
projected baseline gain is not an upper bound on full baseline gain. Thus an
automatic sqrt(20) reduction for a certificate of the *same full-target gain
difference* is not established. It needs an additional residual/alignment bound,
a structural assumption, or an explicitly changed target loss.

## 7. The proposed average-risk objective is coherent, with a genuinely new scope

Let R_h(theta)=E_theta r(h(D),theta), including sampling and any rule
randomization. A precise version is to seek

    E_(theta~Pi)[R_A0(theta)-R_h(theta)] >= epsilon,
    sup_(theta in Theta_cap)[R_h(theta)-R_A0(theta)] <= tau.

This is a meaningful prospective scientific question if Pi reflects a declared
target population of conditions, Theta_cap identifies the populations to be
protected, and epsilon/tau are tied to external relevance before new outcomes.
It is not equivalent to a per-dataset harm limit or to a high-probability
certificate; the distinction must be explicit. Absolute normalized regret also
requires a declared normalization/bridge to the external prediction scale.

The G1-R descriptive grid can motivate this question. It does not supply a
defensible condition distribution, population-level worst-case guarantee, cap
or tuned successor. A maximum over a finite tested grid is a maximum over that
grid, not a supremum over a continuous class. Finite-simulation uncertainty must
be controlled for the declared cap as well as for the average benefit.

Keeping the previously studied generic GCV arm and floor-free arm as candidates
is a reasonable prospective option; it is not adoption or retrospective
optimization. New conditions/seeds alone are insufficient if the criterion,
weights or cap are selected to preserve the old favorable pattern. Candidate
labels and population-class symmetry must also be declared prospectively.

**Programme decision:** certificate-selector development stops here as requested.
The pair law and regression pivots remain analysis tools. This is a scope and
priority decision supported by the failed width screens and information limits;
it is not recorded as a theorem that all certificates are impossible. No new
study is frozen or started. The proposed average-risk-with-absolute-cap question
is scientifically coherent but awaits its external target distribution and cap
before a protocol can be made executable. A0 remains the comparator, the archive
is closed, and no challenger is established.

Reproducibility: `check_two_point_identification.py` generates
`two_point_identification_arithmetic.json` without importing programme code or
reading closed results. It checks likelihood-ratio moments against independent
dense matrix trace identities (maximum absolute discrepancies below 3e-11).
