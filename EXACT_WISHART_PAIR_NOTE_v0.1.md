# Exact Wishart pair law and the remaining inference problem

8 November 2026. Mathematics/prototype note. No programme experiment, new challenger, freeze, or alteration of closed outputs. The artificial Wishart checks below are numerical validation of the derivation only.

## 1. Accepted correction: the OLS Schur remainder has an exact distributional treatment

Use population covariance blocks A=Cov(O), B=Cov(U,O), D=Cov(U), U=(T,P), and C=B A^-1. Let sample covariance have divisor n-1, observed block Ahat, estimated regression Chat=Bhat Ahat^-1, Delta=Chat-C, and R=[-C,I]. Direct expansion yields

    Khat = R Sigmahat R' - Delta Ahat Delta'.

For centered data, let E denote the true Gaussian regression errors and P_O the rank-p projection onto the centered observed design. Then

    R Sigmahat R' = E' M_1 E/(n-1),
    Delta Ahat Delta' = E' P_O E/(n-1),
    Khat = E' M_(1,O) E/(n-1),
    E[Delta Ahat Delta'] = p K/(n-1).

Here M_1 centers and M_(1,O) removes the intercept and p observed regressors. The remainder is negative semidefinite in the expression for Khat. Conditional on O it is a projected Gaussian scatter; it is not an unknown Taylor remainder requiring a generic resolvent bound. Thus Q=(n-1)Khat has exactly Wishart_nu(K) law, nu=n-p-1.

**Correction to the earlier review:** for this OLS estimator, use the exact Wishart representation rather than a generic bound for the Schur estimation remainder. The nonlinear ratio still needs analysis, and this simplification does not apply unchanged to ridge conditioning or data-dependent target scaling.

One precision point: regression errors E have iid rows conditional on O. The n rows of M_(1,O) E do not: their conditional covariance is M_(1,O) tensor K. An orthogonal basis of its nu-dimensional residual space gives nu iid Gaussian coordinates. That is the basis of the Wishart proof. Full-rank Q requires nu>=dim(U); the distribution exists below that threshold too, as a singular Wishart distribution.

## 2. Exact pair law conditional on a two-candidate scatter block

For one pair j,g, order the relevant population covariance as (T,P_j,P_g). Define

    H = K_{(j,g),(j,g)}, B2=K_{T,(j,g)} H^-1,
    Psi=K_TT - B2 H B2',
    A=Q_{(j,g),(j,g)}.

For SPD K and nu>=2, A is nonsingular almost surely. Conditional on A,

    Q_{T,(j,g)} | A = B2 A + Psi^(1/2) N A^(1/2),

in distribution, where N has independent standard Gaussian entries and the factors are chosen to give row covariance Psi and column covariance A. No estimate of the parameters has been substituted here.

Set a=A_jj, b=A_jg, d=A_gg and r=b/sqrt(ad). Let

    u=Q_Tj/sqrt(a), v=Q_Tg/sqrt(d),
    mu_u=(B2 A)_j/sqrt(a), mu_v=(B2 A)_g/sqrt(d).

Their conditional covariances are Cov(u)=Cov(v)=Psi and Cov(u,v)=r Psi. For fixed W>=0,

    Z_j-Z_g = u' W u-v' W v.

Define transformed target vectors

    L=(u+v)/sqrt[2(1+r)], Rv=(u-v)/sqrt[2(1-r)].

They are conditionally independent Gaussians with covariance Psi and means l0 and r0 obtained by the same transformations. The difference is

    Z_j-Z_g = 2 sqrt(1-r^2) L' W Rv.

Let Psi=L0 L0' be a Cholesky factorization, and diagonalize L0' W L0=V diag(lambda_i) V'. Define alpha=V' L0^-1 l0 and beta=V' L0^-1 r0. Then, up to a positive common factor, the strict ordering event is

    sum_i lambda_i (N_i+alpha_i)(M_i+beta_i) > 0,

with all N_i and M_i independent standard normals. Equivalently it is a difference of independent weighted noncentral chi-square terms using (N_i+M_i)/sqrt(2) and (N_i-M_i)/sqrt(2).

The conditional characteristic function of this product-normal sum is

    phi(t) = product_i (1+lambda_i^2 t^2)^(-1/2)
        * exp{[2 i lambda_i t alpha_i beta_i
                 -lambda_i^2 t^2(alpha_i^2+beta_i^2)]
                /[2(1+lambda_i^2 t^2)]}.

For a nonzero metric, the law is continuous in this SPD setting. The conditional probability is obtained by Fourier inversion:

    P(Z_j>Z_g | A)=1/2+(1/pi) integral_0^infinity Im phi(t)/t dt.

For effective metric rank one there is a simpler expression:

    P = Phi(alpha)Phi(beta)+Phi(-alpha)Phi(-beta).

For W=0 every score ties, so the strict ordering probability is zero; the continuity-based formula must not be used in that degenerate case.

## 3. Marginalization is an additional integration, not automatically one-dimensional

The candidate block A itself has Wishart_2(nu,H) law. A Bartlett representation is

    A=L_H T T' L_H',
    T=[[sqrt(X),0],[Y,sqrt(Z)]],
    X~chi-square_nu, Y~N(0,1), Z~chi-square_(nu-1), independently.

Thus the unconditional ordering probability is a three-variable expectation of the conditional expression above. For multiple target directions the implementation additionally performs one-dimensional Fourier inversion at each outer point. This establishes a finite-dimensional deterministic calculation, **not a general reduction to one numerical integral**. A special-case reduction may be possible but has not been derived here.

`wishart_pair_law.py` implements this representation. The outer integral uses Gauss-Legendre rules on inverse-CDF uniforms; orders are explicit. Inner quadrature errors do not bound the outer integration error. Tiny metric eigenvalues are numerically truncated, and this too would require error control for a certified implementation. This is a known-parameter numerical prototype, not a fully precision-certified evaluator.

## 4. Artificial numerical checks and their limits

Three artificial covariances were specified directly, with nu=12. No frozen-family structure, archive observation, or programme run was used. Each deterministic probability was compared with 60,000 independent draws from SciPy's Wishart distribution. The probabilities at outer orders 6/10/16 were:

| Case | order 6 | order 10 | order 16 | Direct Wishart frequency | 99% binomial interval |
|---|---|---|---|---|---|
| scalar target | .70272722 | .70244141 | .70234305 | .69838333 | [.69353306,.70320208] |
| two targets | .55189961 | .55161323 | .55150930 | .55236667 | [.54712537,.55759965] |
| zero-cross-covariance symmetric case | .5 | .5 | .5 | .49751667 | [.49225087,.50278285] |

All finest-order probabilities fall within those Monte Carlo intervals. The conditional swap identity was checked at identical permuted scatter blocks. An initial assertion that swapping candidates must also complement a **finite-order unconditional quadrature** to 1e-5 failed: Bartlett coordinates change under a candidate swap, so discretized outer integration is not exactly invariant. At order 10 its discrepancy was .00022226 and .00032795 for the two non-null cases. That failure is retained here as evidence that outer numerical accuracy is not yet certified; the inappropriate assertion was replaced with a conditional same-scatter identity check, not a claim that unconditional integration was exact.

The non-null conditional Fourier error estimates were around 4.4e-10 after averaging outer weights, while order-10 to order-16 changes were around 1e-4. **It would be false to report the unconditional probabilities as accurate to 1e-10 or 1e-12.** Simulation compatibility is a check, not a proof of that precision or of statistical coverage. `wishart_pair_validation.json` records the settings, convergence, simulation intervals and scope.

## 5. A known-parameter distribution is not yet a confidence interval

Positive common rescaling of W was checked at factors 1e-20 and 1e20, and W=0 was checked separately. The implementation normalizes the metric before numerical rank detection; a tiny overall loss scale must not be mistaken for a zero metric. Relative truncation of very small eigenvalues remains an uncertified approximation.

The pair law is indexed by H, B2 and Psi, not just by the desired population gap f_jg=G_j-G_g. Conditioning on the *observed* A does not remove B2 or Psi. Plugging their estimates into this law gives a plug-in procedure whose finite-n coverage is unproven.

A simple illustration: with T independent of both candidates, f_jg=0. Changing the target variance rescales Z_j-Z_g despite leaving the gap zero. More generally the candidate correlation and conditional target covariance affect the law. This does not rule out constructing a pivot or exact nuisance-adjusted inference, but the raw statistic's known-K law alone does not provide it.

To claim finite-sample coverage one must do one of the following, for example:

- derive an appropriate nuisance-free pivot under a stated restricted model;
- invert tests with correctly controlled nuisance parameters over the composite null;
- construct a valid joint confidence region for K and project it onto all relevant gain differences.

An exact-region route exists in principle. For Khat=Q/nu and Q SPD, the invariant statistic

    T(K;Q)=nu[tr(K^-1 Khat)-log det(K^-1 Khat)-k]

has a distribution depending only on nu and k, since whitening reduces Q to Wishart_nu(I). Its exact quantile, rather than an asymptotic Wilks quantile, defines a joint region. Projection onto gain differences then requires verified optimization bounds. This is a possible inference construction, not a claim that the resulting region is informative or cheaply optimized.

The pivotal distribution itself has a direct scalar characteristic-function representation. With V~Wishart_nu(I), Bartlett decomposition gives independent X_i~chi-square_(nu-i+1), i=1,...,k, and Y~chi-square_h, h=k(k-1)/2, such that

    T = Y + sum_i [X_i-nu log(X_i)] + c,
    c = nu*k*(log(nu)-1).

(For k=1 take Y=0.) Writing a_i=(nu-i+1)/2, integration against each gamma density yields

    phi_T(t) = exp(i*t*c) (1-2*i*t)^(-h/2)
        * product_i {2^(-i*t*nu) Gamma(a_i-i*t*nu)
                    / [Gamma(a_i) (1-2*i*t)^(a_i-i*t*nu)]}.

Powers use the branch continuous from t=0. This supports one-dimensional Fourier evaluation of the pivotal CDF and quantile, with numerical tail/error control still to be supplied. Unlike the raw pair law, nuisance covariance parameters have genuinely disappeared. It is a concrete route to an exact joint confidence region, although **no quantile computation, projected interval or safety rule is implemented in this prototype**.

## 6. Adaptive reference and whole-selector performance remain distinct obligations

If g is fixed independently of the inference sample, simultaneous one-sided comparisons against the other 19 candidates can use a 19-way error allocation. If g(D) is selected from the same data, fixed-pair coverage plus a union bound over the *selected* 19 pairs is not automatically valid. Options include all 190 unordered pairs with two-sided intervals (or 380 directed one-sided comparisons), valid selection-conditional inference, an independently fixed reference, or projection from one valid joint covariance region. The latter gives simultaneous protection without a separate per-pair Bonferroni allocation.

The frozen A0 can remain the random reference under a valid all-pairs/joint-region event. Its own shrunk score need not have a Wishart law: the event concerns the population gains of every index, including whichever index A0 chooses. The inferential statistic and fixed target loss must still be specified correctly. Conditioning naively on A0's selected index changes the sampling law and does not solve the selection issue.

Similarly, known pairwise misordering probabilities are not automatically exact probability of identifying the best of 20 or exact selector regret. With unique best candidate j*,

    P(wrong pick)=P(union over j!=j* of {Z_j>Z_j*}),
    max_j p_j <= P(wrong pick) <= min(1,sum_j p_j).

The events share the targets and candidates and are dependent. Exact mean selection regret needs P(j_hat=j) for each j, not merely P(Z_j>Z_j*). Thus pair calculations can supply useful bounds, but should not be reported as a full exact regret calculation without the additional joint analysis.

Nor can departure probability of a data-dependent certified rule be inferred from interval widths at the population K alone. Its acceptance event, and A0's pick, depend jointly on the reference data. That joint probability needs its own analysis or a new, prospectively authorized evaluation.

Finally, exact error probabilities for one specified OLS rule answer a different question from a Le Cam lower bound over all rules. They sharpen the evaluation of that estimator; they do not replace a general identification lower bound.

## 7. Two further corrections to proposed wording

1. Candidate-independent **strictly increasing** score transforms preserve ordering. Merely nondecreasing transforms can introduce ties and change a lowest-index tie-broken pick; decreasing transforms reverse ordering. A common strictly increasing transform is the correct action-neutral class.
2. A Loewner sandwich does not generally bound the cross-covariance numerator multiplicatively. For K=I_2 and Khat=[[1,e],[e,1]], 0<e<1, we have (1-e)K <= Khat <= (1+e)K, yet true c^2/v=0 and estimated c^2/v=e^2>0. A universal “b/a compounded” gain-width formula would therefore be wrong without extra conditions. The earlier broad covariance table motivates caution, not a universal impossibility result for every gain gap.

The proposed warning about model constraints is accepted: if a true population is excluded by a mistaken mechanistic constraint, a certificate built assuming that constraint can lose coverage. Constraint validity and acquisition benefit are separate questions.

## 8. Revised next calculation

1. The known-K conditional pair law is now derived and prototyped. Improve and certify outer integration error if high precision is needed; no claim of one-dimensional unconditional integration yet.
2. Choose and prove the **inference construction**: identify nuisance treatment, simultaneous/random-reference protection, and the exact safety target. Only then do interval width and departure calculations have a defined meaning.
3. Evaluate known-parameter pair risks, bounds on the full selector, and the chosen certificate on new mathematical structures under a clearly declared prospective scope. Do not call the resulting rule adopted or promote it into external validation unless it meets the new design's criteria.

This accepts the user's simplification of the OLS remainder while declining the unsupported inference that an exact pair sampling law already supplies a finite-sample certificate with constants. The remaining challenge is smaller and better defined, but it is not merely numerical interval-width evaluation.

Coverage and expected harm must also use compatible budgets: when normalized regret is bounded by one, a valid simultaneous event giving harm at most tau with probability 1-delta yields expected excess regret at most tau+delta. A conventional 95% confidence level alone does not certify a much smaller expected-harm budget.
