# Boundary calibration, realistic power bounds and symmetry

8 November 2026. Mathematical addendum to the earlier identification notes.
Earlier versions remain unchanged. Only known-population probability calculations
and algebra are performed: no Monte Carlo, coverage study, profile-likelihood
selector, experiment, freeze, archive analysis or frozen-file change.

## 1. Main finding: a closer boundary null gives a stronger power ceiling

The 26% power in the previous note belongs to a particular conditional-level
test of two *simple, opposite-gap populations*. The null there has f=-Delta,
not f=0. It cannot be transferred to an interval rule calibrated over f<=0.
It was also not proved to be an unconditional Neyman-Pearson power ceiling:
conditional error allocation is an additional restriction.

There is a cheap and sharper bound for the broad fixed-W Gaussian
unknown-covariance model used in the earlier Wishart discussions. It avoids
implementing a profiled-likelihood rule merely to screen its possible power.

### Construct an admissible boundary population

Use population-standardized residual candidate coordinates Z, with true
covariance I, and T=F Z+epsilon. Keep F and the independent target noise law
fixed. The leading gains under W=I are G_j=1 and G_g=1-Gamma. Change only the
two candidate variances to

    d_j=(G_j+G_g)/(2G_j)=(2-Gamma)/2,
    d_g=(G_j+G_g)/(2G_g)=(2-Gamma)/(2(1-Gamma)).

The conditional law T|Z is identical at the two populations. Under the new
population the gains are d_j G_j=d_g G_g=(G_j+G_g)/2, so f_jg=0. Their sum
and hence R0 are preserved. Other candidate variances and gains are unchanged.
Thus this is a boundary null for the same normalized gain comparison as well.

The entire KL comes from the two candidate variances:

    KL(P_positive || P_boundary)
       =1/2 [1/d_j+1/d_g-2+log(d_j d_g)]
       =1/2 log[(2-Gamma)^2/(4(1-Gamma))].

For Gamma=.02, d_j=.99, d_g=1.0102040816 and KL per row is
**0.0000510178053**, versus 0.00121215213 for the D1-like opposite-magnitude
swap. The new null is substantially closer. This identity is independent of
loading directions, other columns, target dimension and noise variance.

**Scope of admissibility:** candidate residual variances must be unknown and
allowed to vary. This null belongs to the unrestricted SPD covariance model.
It is not admissible if the model externally fixes these variances, or imposes
known exact equality across candidates. Such a restricted model needs a
different power analysis and an independently justified restriction. The literal
frozen labelled generator is not silently substituted for this inference model.

### Exact Neyman-Pearson reduction

For s residual Gaussian coordinates let S_j=sum Z_ij^2 and S_g=sum Z_ig^2.
The log likelihood ratio of the positive population against the boundary is

    ell=(s/2)log(d_j d_g)+(q/2)(S_j-S_g),
    q=Gamma/(2-Gamma)>0.

Consequently the most powerful size-alpha test rejects for S_j-S_g>c_alpha,
where

    P_boundary(d_j U-d_g V>c_alpha)=alpha,
    U,V independently chi-square_s.

Its power is

    beta_NP=P_positive(U-V>c_alpha).

Every uniformly level-alpha fixed-pair test over f<=0 has size at most alpha
at this particular boundary population. Neyman-Pearson therefore gives the
necessary ceiling power<=beta_NP at the positive population. No bootstrap,
studentization, profile optimizer or estimator can evade this bound while
retaining that size guarantee and model class. It is not an assertion that a
practical method attains the ceiling.

This supplies a concrete boundary calibration case for any proposed GLR. A
larger fixed-pair power claim would require a different size condition/model,
or would conflict with this simple-null bound.

The calculation is one-dimensional. For positive a,b,

    aU-bV has the same distribution as
       (a-b)X+2sqrt(abX)N,
    X~chi-square_s, N~N(0,1), independently.

The equality follows from the common moment generating function
[(1-2at)(1+2bt)]^(-s/2). Thus the required tail is

    P(aU-bV>c)
       =E Phi(((a-b)X-c)/(2sqrt(abX))).

`check_boundary_power.py` evaluates this expectation and solves for c_alpha.
The mathematical ceiling is exact; the displayed numerical evaluations are
approximations, not formally interval-certified decimal enclosures.

| s | beta_NP, alpha=.05 | beta_NP, alpha=.01 |
|---:|---:|---:|
| 81 | .0601412 | .0127283 |
| 389 | **.0741555** | **.0167215** |
| 799 | .0870283 | .0206458 |
| 1619 | .1077822 | .0274425 |
| 3259 | .1427113 | .0400899 |

For the specified common-nuisance full-data extension with n=410, the ceilings
are .0749024 and .0169427. The residual-scatter experiment uses s=nu=389.

Numerical checks: adaptive integration over chi-square quantile limits with
omitted mass <=2e-13 was cross-checked with normalized Gauss-Laguerre rules
of orders 48 and 80. Maximum discrepancy was 1.72e-12. An independent F-CDF
identity at c=0 agreed within 4.79e-13; the covariance-matrix KL agreed within
1.32e-15. Quadrature estimates and convergence are not themselves rigorous
error certificates; they support the rounded values used for this screen.

## 2. Why the previous oracle calculation does not transfer

Even holding candidate variances fixed, the known equal-gain midpoint is closer
than the opposite-gap population. If F_+,F_- are the two prior coefficient
matrices, let F_0=(F_++F_-)/2. Its two gains are equal. The coefficient
difference from F_+ is half as large, so KL and D in the conditional normal
likelihood-ratio law are one quarter of their opposite-endpoint values.

For D1-like n=410 and alpha=.05, the *conditional-level* test powers become:

| Known-regression comparison | Orthogonal loading directions | Collinear loading directions |
|---|---:|---:|
| Opposite-gap simple null | .258498 | .258463 |
| Equal-gain midpoint null | .125809 | .125795 |

These are feasible conditional-level test powers, not unconditional NP ceilings.
The midpoint test covers only that specified regression line, not the full
unknown-covariance composite null. Section 1's variance perturbation gives the
stronger and relevant ceiling when candidate variances are unknown.

A level-alpha interval/test need only have false-rejection probability <=alpha
over its null. Equality at f=0, or negligible probability farther inside the
null, does not follow without a specified method and nuisance model. Random
reference selection also requires the previously agreed simultaneous treatment.
The numerical gap -.03 alone does not establish a tail probability without its
sampling scale and calibration.

There is another distinct transfer issue: the earlier illustrative reference
was always the second candidate. A0 is data-dependent and can select other
candidates. A single known-pair power cannot be multiplied by a single gap to
obtain the full A0-relative expected improvement.

For a valid simultaneous procedure, Section 1's bound *does* apply to the
unconditional event that it selects this particular pair and certifies a
departure: at the boundary that event has probability at most alpha. At D1,
the top-two normalized gap is .003. Thus its optimistic top-two contribution is
at most about

    .003 * .0741555 = .000222467  (residual-scatter experiment),
    .003 * .0749024 = .000224707  (full-data oracle extension).

Both are below the quoted .000257 planning margin, rather than three times
larger. That margin comes from the previously reported baseline mean and is
not a newly known exact population risk. Other baseline picks and other
departures could contribute additional improvement. Therefore this arithmetic
neither reverses a frozen decision nor proves that a full A0-relative certified
policy cannot meet H1. It invalidates the supplied numerical demonstration that
it can. A positive harm allowance instead of LB(f)>0 would also define a
different null and require different bounds.

## 3. The symmetry argument has a valid repair, but is not immediate

For a permutation-equivariant decision rule and permutation-invariant loss,
R_h(theta)=R_h(pi theta). Hence its risk at a parameter equals its risk averaged
over pure candidate permutations of that same parameter. This part is correct,
when equivariance holds for the complete rule.

But the low-KL alternative is a *magnitude swap*:

    F=[a u,b v,...],    F_alt=[b u,a v,...],

whereas a pure candidate permutation gives

    F pi=[b v,a u,...].

They are generally different populations. Already with equal candidate variances,
the pure-permutation KL per row is

    (a^2+b^2-2ab rho)/sigma^2,

whereas the low-KL magnitude-swap value is (a-b)^2/sigma^2. At roughly
orthogonal directions the first can be large even for a tiny gain difference.
Candidate equivariance alone therefore does not equate the risk at a particular
D1 population with risk at its low-KL magnitude swap.

There is a correct minimax version within a **labelled structural class allowing
all loading directions**. In the conditional model define

    theta_bar: F_bar=[a v,b u,...].

Then pi theta_bar=theta_alt. If theta_bar is in the labelled class and h is
equivariant, e_h(theta_alt)=e_h(theta_bar). The two-point bound consequently
implies max(e_h(theta),e_h(theta_bar))>=R_star. An exchangeable prior on these
directions likewise yields a prior-averaged lower bound. This is a nuisance-
direction minimax/average statement, not a pointwise bound at the single frozen
theta. The literal finite seeded condition list need not contain theta_bar.
For full-data arms, additional observed-block and marginal-variance constraints
must be respected by the joint extension too.

### Engineering qualifications for the frozen arms

The score formulas commute with candidate permutations in real arithmetic;
unique-score maximizers then give equivariant picks. Exact ties require care.
[pick()](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/g1r_run.py:220)
uses lowest-index `np.argmax`, so it is not equivariant on tied scores.
At complete Ledoit-Wolf shrinkage all candidate scores tie. The clipping
intensity can equal one on an open set around isotropic sample covariances,
so an assertion that this is always a measure-zero exception is unavailable
over the broad Gaussian class. Zero observed ties in the completed grid do not
establish a universal symmetry theorem.

Also [the LW intensity](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/g1r_run.py:255)
uses fourth moments of the standardized data; A4 is not a function of the sample
covariance alone. It can still have permutation-invariant tuning, but the
data-level symmetry statement is the appropriate one.

Consequently a covariance-only Wishart bootstrap does not by itself reproduce
A4's full data-dependent tuning. This does not affect a new test specified to
use residual scatter alone, but it matters if covariance sufficiency is claimed
for every archived arm.

No frozen arm or tie rule is changed. A future symmetry theorem should specify
unique-maximizer scope or a permutation-compatible randomized tie convention.

## 4. Updated programme decision

The revised umbrella is accepted: there is no standalone certificate-selector
programme. Any candidate, certified or unconditional, is assessed against the
same prospective average-risk objective and absolute expected-harm cap:

    E_(theta~Pi)[R_A0(theta)-R_h(theta)] >= epsilon,
    sup_(theta in Theta_cap)[R_h(theta)-R_A0(theta)] <= tau.

The condition distribution, protected class, normalization, epsilon and tau
remain to be declared from the external scientific target before new outcomes.
A Bayes rule need not have uniform frequentist interval coverage; its declared-
prior appeal does not establish the worst-case expected-harm cap. That cap
requires its own stated validation or proof.

The new fixed-pair ceiling is a useful inexpensive screen for broad unknown-
covariance certificates. It does not eliminate all certified policies or define
their actual power. It also does not justify assuming they have near-zero harm
at D2. Their full reference-selection and comparison rule is still unspecified.
No profiled GLR has been implemented, since the 26% premise can already be
corrected by this boundary calculation and the statistical model must be fixed
before profiling it.

A0 remains the comparator; no challenger is established; the archive remains
closed. H2 was not established for A2, and the frozen joint criterion remains
unmet. This addendum introduces no new scientific result from that archive.

Reproducibility: `check_boundary_power.py` writes
`boundary_power_arithmetic.json`. It imports no programme code, reads no study
outputs, generates no random observations, and changes no frozen files.
