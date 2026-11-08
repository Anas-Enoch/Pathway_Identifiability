# Functional inference: exact width witnesses and bootstrap obligations

8 November 2026. Mathematical review and deterministic arithmetic only. This note
preserves `EXACT_WISHART_PAIR_NOTE_v0.1.md`. It implements no interval procedure,
bootstrap, quadrature replacement, coverage simulation, challenger or protocol.
No closed programme output is used as data or modified.

## 1. The full region is too wide in the near-tie regime: an exact witness

The proposed dimension warning is justified for this regime, but a Wilks-scale
dimension count alone is not a proof about a particular functional. A constructive
witness avoids both that inference and a full pivot-quantile calculation.

For Q~Wishart_nu(K), Ktilde=Q/nu, and k=40, use

    T(K;Q)=nu[tr(K^-1 Ktilde)-log det(K^-1 Ktilde)-k].

The previous Bartlett representation gives

    T=Y+sum_i [X_i-nu log X_i]+nu*k*(log nu-1),
    Y~chi-square_h, h=k(k-1)/2=780.

Because x-nu log x >= nu-nu log nu for every x>0, T>=Y in
this coupling. Thus, without evaluating T's quantile,

    q_T(.95) >= q_chi-square_780(.95) = 846.083516073059.

Smaller failure probabilities only increase the cutoff. The .05 used here is an
illustrative coverage level, not an adopted expected-harm budget.

Exact moments, with a_i=(nu-i+1)/2, are

    E T = nu[k log(nu/2)-sum_i digamma(a_i)],
    Var T = nu^2 sum_i trigamma(a_i)-2nu*k.

| n_train | nu | Exact E T | Exact sd T |
|---|---:|---:|---:|
| 102 | 81 | 1006.2984 | 50.3287 |
| 410 | 389 | 850.0191 | 41.9934 |
| 820 | 799 | 834.2082 | 41.2017 |
| 1640 | 1619 | 826.9201 | 40.8394 |
| 3280 | 3259 | 823.4159 | 40.6658 |

The 820-parameter Wilks heuristic is reasonable at the largest nu but misses
substantial finite-n displacement at nu=81. The proposed roughly 33% and 11%
eigenvalue slacks describe equal allocation across 40 directions. They are not
individual-direction bounds. Using just the cutoff lower bound, one direction
alone can have an eigenvalue of K^-1 Ktilde between about .04366 and 4.72866
at nu=389, or .44108 and 1.90310 at nu=3259. These are boundary solutions of
x-1-log x=q_lower/nu, not full quantiles or widths of gain intervals.

### A gain-specific reversal inside the region

Standardize the residual candidates by their population marginal scales. This is
an invertible candidate-coordinate transformation and does not change any gain.
The Gaussian conditional family can then be written

    Z~N(0,I_20), T=F Z+epsilon,
    epsilon~N(0,sigma^2 I_20), independently,
    K_F=[[FF'+sigma^2 I,F],[F',I]].

Take W=I; a fixed common positive target scaling cancels in normalized gaps.
Write the columns F_j=sqrt(G_j) u_j with ||u_j||=1. Set G_1=1,
G_2=1-Gamma. Form F_alt by changing only their magnitudes:

    F_alt,1=sqrt(G_2) u_1,
    F_alt,2=sqrt(G_1) u_2.

All other columns and sigma^2 are unchanged. This exchanges the two gains and
preserves R0=tr(K_TT), since the sum of column norms squared is unchanged. It
works for every orientation of the u_j and both candidate-variance constructions.

The determinants of K_F and K_Falt agree. Block Gaussian regression gives the
exact trace identity

    tr(K_Falt^-1 K_F)-k = ||F-F_alt||_F^2/sigma^2
                          =2(1-sqrt(1-Gamma))^2/sigma^2.

The alternative is admissible in the unrestricted SPD region under discussion.
This argument does not automatically apply to a different confidence region
that imposes additional mechanistic constraints and separately justifies them.

At the *population-center scatter* Q=nu*K_F,

    T(K_Falt;Q)=2nu(1-sqrt(1-Gamma))^2/sigma^2.

In the family calibration sigma^2=.25/(20*A), R0=1/A, hence

    T_alt=160nu*A*(1-sqrt(1-Gamma))^2,
    f_12/R0=Gamma*A, f_alt,12/R0=-Gamma*A.

For A=.20, which maximizes this cost over the stated amplitudes:

| nu | Gamma=.02 | Gamma=.05 | Gamma=.20 |
|---:|---:|---:|---:|
| 81 | .2618 | 1.6618 | 28.8894 |
| 389 | 1.2574 | 7.9808 | 138.7407 |
| 799 | 2.5827 | 16.3924 | 284.9712 |
| 1619 | 5.2333 | 33.2157 | 577.4322 |
| 3259 | 10.5344 | 66.8623 | 1162.3542 |

Every Gamma=.02 or .05 entry is below 846.0835, including the largest sample
size. The full confidence region therefore contains an opposite ordering at its
population-center input for all these near-tie profiles. No accurate optimization
or larger full-pivot cutoff can remove this particular ambiguity.

For Gamma=.02, A=.15, the normalized gap .003 becomes -.003 under the alternative;
the largest-n pivot cost is only 7.9008. A lower bound projected from this full
region cannot certify this departure with harm allowance tau<.003 at that input.

The deterministic arithmetic checks 18 profile/variance combinations at five
degrees of freedom: 86/90 have a reversal certified inside the region by the
cutoff lower bound alone. The other four are **inconclusive for this witness**,
not positive certificates. This is not 90 simulated experiments. The identity
applies analytically to every loading orientation; dense deterministic matrices
only check the block algebra (maximum trace-cost discrepancy 9.1e-15 per degree
of freedom).

**Scope:** this establishes an uninformative central-data diagnostic for the
specified comparisons. It is not a proof that every random observed Q produces
an uninformative interval, or that every possible exact inference method fails.
The argument does not evaluate A0's random pick or any selector's performance.

## 2. Exact simultaneous finite-sample validity is available without an 820-parameter region

The absence of a nuisance-free distribution for the raw paired gap does not imply
that all nuisance-free pivots relevant to gains are unavailable. A classical
Gaussian regression construction gives a counterexample to the proposed
impossibility statement.

For candidate j, let m=20, v_j=K_jj, b_j=K_Tj/v_j and
Psi_j=K_TT-v_j b_j b_j'. From the marginal Wishart block define

    w_j=Q_jj,
    bhat_j=Q_Tj/w_j,
    S_j=Q_TT-Q_Tj Q_jT/w_j.

Conditional on w_j, sqrt(w_j)(bhat_j-b_j)~N_m(0,Psi_j), independently of
S_j~Wishart_(nu-1)(Psi_j). Also w_j/v_j~chi-square_nu. Whitening therefore gives
the exact regression pivot

    [w_j*(nu-m)/m] (bhat_j-b_j)' S_j^-1 (bhat_j-b_j)
        ~ F_(m,nu-m),               nu>m.

For completeness, if z~N_m(0,I) and S~Wishart_d(I) independently, rotational
invariance gives z'S^-1 z distributed as chi-square_m/chi-square_(d-m+1),
with the two scalar variables independent. Set d=nu-1 to obtain the displayed
factor. This factor includes the divisor of the *scatter* S_j, not an estimated
covariance with an implicit degrees-of-freedom divisor.

For any eta>0 define

    c_eta = m/(nu-m) F_(m,nu-m)^(-1)(1-eta),
    C_b,j = {b: w_j(b-bhat_j)'S_j^-1(b-bhat_j)<=c_eta},
    v_L,j = w_j / chi-square_nu^(-1)(1-eta/2),
    v_U,j = w_j / chi-square_nu^(-1)(eta/2).

Each slope set and each variance interval has exactly its stated marginal
coverage. For P=20 choose eta=delta/(2P). A union bound makes all 40 events
hold with probability at least 1-delta, without any independence assumption
between candidates. Since G_j=v_j b_j'W b_j and W>=0, define

    L_j=v_L,j * min_(b in C_b,j) b'Wb,
    U_j=v_U,j * max_(b in C_b,j) b'Wb.

Then all 20 gains are simultaneously in [L_j,U_j] with probability at least
1-delta. Consequently [L_j-U_g,U_j-L_g] covers every difference, including
one indexed by A0's data-selected reference. This is a derived construction,
**not an implemented or adopted interval method**.

The quadratic projections can be bounded conservatively by norm inequalities
or optimized over the ellipsoid. For rank-deficient W one may first transform
to its nonzero target coordinates and use that dimension instead of m.

This does **not** claim a scalar pivotal law for the gap, equality coverage of
the final simultaneous set, optimal width, or a useful challenger. The relevant
distinction is exact finite-sample *validity* (coverage >=1-delta) versus exact
nominal equality. A safety certificate needs the former. Thus “no exact
finite-sample inference exists” and “bootstrap is the only available construction”
are not warranted.

An unknown normalizer R0 must also be covered. Fixed W does not make
tr(W K_TT) known. One option is to diagonalize W and simultaneously bound the
variances of its target coordinates using their chi-square pivots, allocating
additional error probability. Those bounds give an interval for R0, from which
ratio bounds follow. A population R0 used to express diagnostic widths is not
a data-available normalizer.

### Pricing this conservative construction before implementing it

Validity does not make these particular intervals attractive. At Q=nu*K_F,
S_j=nu*Psi_j and Psi_j>=sigma^2 I. Therefore C_b,j contains the Euclidean ball
of radius sqrt(c_eta*sigma^2/v_j) centered at the true b_j. Even *exactly*
projecting this ellipsoid cannot eliminate the gains attainable in that ball.
For W=I, write r=sqrt(c_eta*sigma^2) and

    a_L=nu/chi-square_nu^(-1)(1-eta/2),
    a_U=nu/chi-square_nu^(-1)(eta/2).

The product-set projection on G_j has lower endpoint no greater than
a_L*(sqrt(G_j)-r)_+^2 and upper endpoint no smaller than
a_U*(sqrt(G_j)+r)^2. These yield an orientation-independent lower bound on
the width of this specific gain-set construction.

At Gamma=.02, A=.15, delta=.05, the true normalized gap is .003:

| nu | Lower bound on projected pair-interval width, expressed in units of the true R0 |
|---:|---:|
| 81 | .75150 |
| 389 | .26860 |
| 799 | .18211 |
| 1619 | .12621 |
| 3259 | .08838 |

These widths concern the product of the individual slope/variance sets, without
extra joint constraints. They are not lower bounds on *all* possible exact
functional intervals. Converting by the population R0 here is diagnostic only.
Nevertheless they are sufficient to reject this elementary construction as a
promising near-tie certificate at its population-center input. No implementation
or coverage simulation is needed to learn that.

## 3. Bootstrap is an admissible candidate, not automatically second-order accurate

I do not require an equality-exact pivot before considering a successor.
Approximate simultaneous inference can be investigated with explicitly stated
coverage evidence and limitations. But Gaussian inputs and a smooth algebraic
formula do not themselves establish second-order bootstrap accuracy.

At c_j=0 the gain G_j=c_j'Wc_j/v_j has zero first derivative and a nonzero
second derivative. The ordinary centered bootstrap has a familiar failure at
such points. In scalar form, if sqrt(nu)*chat converges to Z~N(0,1) at c=0,
then nu*chat^2 converges to Z^2. A parametric bootstrap with
sqrt(nu)(chat*-chat) approximating an independent Z* instead gives

    nu[(chat*)^2-chat^2] -> 2 Z Z*+(Z*)^2  conditional on Z.

This conditional law depends on Z and is not the target chi-square law. Ordinary
delta studentization has a zero population denominator there. For example, the
plug-in standard error 2|chat|/sqrt(nu) gives a nonnegative limiting original
studentized statistic, whereas the usual bootstrap version can be negative.
This is a warning about standard constructions, not a theorem that every
modified or null-restricted bootstrap must fail.

The general first-order-degeneracy issue and possible modified bootstrap remedies
are studied by [Chen and Fang, *Inference on Functionals under First Order
Degeneracy*](https://arxiv.org/abs/1901.04861), especially the squared-mean example
and the standard-bootstrap failure analysis. The simple Gaussian calculation
above is also a direct derivation in the present setting.

Thus a proposed simultaneous bootstrap needs its statistic, centering,
studentizer, zero/near-zero handling and max-statistic calibration specified
before a coverage claim. Even pointwise second-order theory in a regular region
would not automatically be uniform near nulls, ties or nearly degenerate
studentizers. These are exactly appropriate prospective failure families.

There is a useful known bias correction: for the fixed-W Wishart scores,

    E[(Z_j-Z_g)/(nu-1)]=G_j-G_g.

This is action-neutral but may matter for interval centering. Action neutrality
does not imply inferential irrelevance. A bootstrap must also declare whether
its covariance plug-in is Q/nu or a differently scaled Schur estimate; the
Q/(n-1) versus Q/nu distinction cannot be ignored for gain magnitudes.

## 4. Simulation measures coverage; it does not answer it exactly

At a fixed declared K, finitely many coverage trials give a binomial estimate
and a confidence bound. They do not establish the exact coverage number or
uniform coverage over a continuum of K. Independent reproduction and a
prospective specification remain good requirements; neither changes this limit.

For example, with zero coverage failures among N trials, a one-sided 95% upper
bound on failure probability is 1-.05^(1/N). To bound it by 1e-4 requires
N>=29,956 at one point. Protecting 360 declared points simultaneously at 95%
confidence by Bonferroni requires 88,814 zero-failure trials *per point*;
for a 1e-5 budget it requires 888,180. These are illustrative arithmetic, not
an approved simulation size or a claim that zero failures will occur.

For a mean-harm target, normalized regret bounded by one gives the elementary
bound expected excess regret <= tau+delta when the simultaneous event guarantees
harm <=tau. The confidence failure budget must be chosen together with the
external absolute harm budget. A default 95% interval does not by itself certify
a 1e-4 mean-harm bound. Better bounds may be possible if harmful failure events
have structure, but that must be proved or independently bounded.

## 5. Quadrature recommendation: plausible improvement, not a precision certificate

Using density-adapted quadrature is reasonable. For X~chi-square_nu write
X/2~Gamma(a,1), a=nu/2; generalized Gauss-Laguerre has weight
z^(a-1) exp(-z). Use a=(nu-1)/2 for the second Bartlett chi-square and a
standard-normal Gauss-Hermite rule for Y. Such rules integrate weighted
polynomials through their prescribed degree exactly, not arbitrary conditional
ordering probabilities exactly.

Inverse-CDF maps have endpoint derivative singularities, but a composition with
the conditional probability can partly cancel them. They are a plausible source
of slow convergence, not an established exclusive cause of the observed drift.

A normalized Jacobi-matrix construction avoids forming Gamma(a), whose value
overflows ordinary floating point for the largest a here. For an order-q rule,
the Gamma(a,1) Jacobi matrix has diagonal a+2i, i=0,...,q-1, and off-diagonal
sqrt((i+1)(a+i)), i=0,...,q-2. Its eigenvalues are the nodes; the squared first
components of normalized eigenvectors are probability weights. The corresponding
standard-normal Jacobi matrix has zero diagonal and off-diagonal sqrt(i+1).
This puts the nodes near the mass automatically. Further centering/scaling may
improve conditioning; it is not a prerequisite simply to put nodes there.

Geometric convergence is not guaranteed just by real smoothness. An analytic
continuation/growth argument or a suitable verified error bound is still needed.
Nor can 60,000 Monte Carlo draws certify deterministic error of 1e-8: near
p=.7 their standard error is about .00187. A normal-approximation calculation
puts 99% Monte Carlo half-width 1e-8 at roughly 1.4e16 trials, not 60,000.
Simulation is a coarse independent check, not the high-precision reference.

No quadrature code was changed in this review. Improving known-K pair integration
can remain a numerical project, separately from proving interval coverage.

## 6. What the review changes

The full 820-parameter region should not be developed as the near-tie certificate:
the explicit central-data reversal witness already answers that practical screen.
The elementary exact regression-pivot product construction is valid but also
too wide in the same screen. Neither conclusion proves that all exact or
approximate functional inference is impossible.

The next candidate can be a simultaneous parametric bootstrap if its nonregular
boundary handling and harm/coverage budgets are specified first. Its first
mathematical test should include zero and local-to-zero cross-covariances, not
only strong regular alternatives. Only after that specification should a coverage
design be prospectively fixed and implemented. Width and coverage are separate
outputs; favorable width does not compensate for undercoverage.

For a matched fixed conditional K, the OLS Wishart law does not depend on the
observed covariance condition number. That axis can affect A0 and departure
events, but is not a separate nuisance for pure residual-Wishart interval
coverage. Changing K's random seed alongside that axis would confound this
invariance unless deliberately declared as an additional family variation.

Sample splitting is valid when the reference is fixed by an independent part
and inference uses the other part. It changes the reference rule's data, and
therefore the comparison, as the user noted. All-pairs/full-gain simultaneous
coverage can instead protect the original full-data A0 index.

A0 remains the comparator; no challenger is established; the G1-R archive stays
closed. A rigorous negative or positive result may be scientifically useful,
but publishability and novelty are not established by this calculation.

Reproducibility: `check_functional_inference_arithmetic.py` writes
`functional_inference_arithmetic.json`. It reads no programme files and runs no
random sampling. The full Stein-pivot quantile was not computed.
