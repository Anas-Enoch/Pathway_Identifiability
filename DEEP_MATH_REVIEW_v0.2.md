# Deep mathematical review of Desktop/books, v0.2

8 November 2026. Scope: mathematical analysis and source review. No new selector experiment, no new freeze, no archive reuse for tuning. This extends the earlier bounded review; it does not replace any programme record.

## 1. What was actually available and read

The directory `/Users/nour/Desktop/books` contains **19 PDFs / 6,832 PDF pages**. Seventeen have extractable text. Woodhouse's geometric-quantization file (164 pages) and the topological-manifolds file (406 pages) have effectively no extractable text and were not substantively reviewed. No OCR engine was available. Van der Vaart's *Asymptotic Statistics* is a 47-page excerpt; the Cambridge file is a 56-page preview of *Systems Biology: Constraint-based Reconstruction and Analysis*, not a complete book. The Olive statistics volume is not Vershynin's book. The Ortiz manifold text is a dissertation.

All 19 files were hashed and inventoried in `inventory.json`. Extraction and keyword search are not reading. The source ledger at the end distinguishes passages read from contents screened and unread scans. This review goes deeper into the relevant mathematics; it does not claim cover-to-cover reading of 6,832 pages.

The central conclusion remains: a useful successor needs an identifiable decision objective and informative uncertainty for that objective. Changing geometry, regularization or score bias alone is insufficient. Several apparently promising directions can be excluded algebraically before any new experimental grid.

## 2. Set the mathematical problem in fixed coordinates

Let O have dimension p, T dimension m, and P contain q candidates. In a joint Gaussian model with nonsingular observed covariance, write

    U = (T,P),
    K = Cov(U | O) = D - B A^-1 B',
    A = Cov(O), B = Cov(U,O), D = Cov(U).

For a fixed positive-semidefinite target metric W,

    c_j = K_Tj, v_j = K_jj > 0,
    G_j = c_j' W c_j / v_j,
    R0 = tr(W K_TT),
    r(j) = [max_k G_k-G_j]/R0.

For nondegenerate targets, 0 <= r(j) <= 1. This is population selection regret, not fitted-predictor risk or pathway-identification error. All derivations below assume the coordinates/metric are fixed unless explicitly stated otherwise. Random sample target standardization changes these calculations and cannot be silently substituted.

## 3. Exact Gaussian residual law: conditioning can remove a nuisance, with a degrees-of-freedom cost

**Derived result, classical ingredients; no novelty claim.** With n iid Gaussian rows and unknown mean, let M_O be the orthogonal residual-maker for a fitted intercept and all p observed variables. Suppose their design has rank p+1. Then rank(M_O) is

    nu = n-p-1.

Conditional on O, the multivariate Gaussian regression representation is U = 1 mu' + O B_reg + E, with E independent of O and iid rows of covariance K. Hence

    Q = U' M_O U = E' M_O E ~ Wishart_nu(K).

Proof: write M_O = V diag(I_nu,0) V' using an orthogonal V. Gaussian rotational invariance makes V'E a matrix with independent N(0,K) rows. Its first nu rows generate Q. The conditional law does not depend on O, so the same law holds unconditionally.

This is useful because it isolates conditional-covariance sampling from the unknown regression coefficients. It does not eliminate uncertainty: estimating the intercept and O regression consumes p+1 degrees of freedom. The full joint sample-covariance Schur complement is Q/(n-1), not Q/nu. Confusing the divisors changes scores and ridge scales.

The law holds for the **OLS residual projection**, an idempotent matrix. A ridge residual-maker is generally not idempotent: the residual scatter is a weighted sum of Gaussian outer products rather than this simple Wishart law. It also does not survive arbitrary data-dependent target scaling as written.

Source bridge: Drton et al. explains Gaussian Schur complements and Wishart geometry; Shapiro et al., PDF pp.114-115, explicitly identifies sample scatter's Wishart law and degrees of freedom. The regression-projection argument here supplies the particular conditional-residual version rather than assuming it follows from the word “Gaussian.”

## 4. An exact bias correction that cannot improve the acquisition choice

Use the residual scatter Q above, and define the raw candidate statistic

    Z_j = Q_Tj' W Q_Tj / Q_jj.

Let t = tr(W K_TT). Conditional on Q_jj = w,

    Q_Tj | w ~ Normal[(c_j/v_j)w, w(K_TT-c_j c_j'/v_j)].

Therefore

    E[Z_j | w] = (G_j/v_j)w + t-G_j,
    E Z_j = (nu-1)G_j+t.

Also E tr(W Q_TT) = nu t. For nu>1,

    U_j = [Z_j - tr(W Q_TT)/nu]/(nu-1)

is exactly unbiased for G_j.

**But U_j is the same positive affine transformation of Z_j for every j. Consequently argmax_j U_j = argmax_j Z_j on every data set.** This remains true with heterogeneous candidate variances. Unbiasedness does not solve the selection problem. A successor that merely implements this correction would make identical acquisition decisions.

This statement is limited to this fixed-metric, zero-ridge, OLS-residual setting. It is not an unbiasedness claim for the frozen, sample-standardized and jointly shrunk G1 rule. Nor does candidate-wise unbiasedness imply the maximum score is unbiased.

This is a concrete lesson from the books and the failed programme: prove whether an estimator modification can change an action before testing it. A common shift/scale can change a plotted score while leaving every pick unchanged.

## 5. Decision-sensitive geometry: a low-rank derivative inside a large covariance

For a symmetric perturbation E of the full covariance Sigma, differentiation of K gives

    dK[E] = E_UU - E_UO A^-1 B' - B A^-1 E_OU
                     + B A^-1 E_OO A^-1 B'
           = R E R',  R = [-B A^-1, I].

For one candidate,

    dG_j = 2 c_j' W dc_j/v_j - G_j dv_j/v_j.

Let H_j be the symmetric gradient in K coordinates: its T,j block is W c_j/v_j, its j,T block is the transpose, its j,j entry is -G_j/v_j, and its remaining entries are zero. Then

    dG_j = tr(H_j dK),
    rank(H_j) <= 2.

For f_jg = G_j-G_g,

    d f_jg = tr(H_full E),
    H_full = R' (H_j-H_g) R,
    rank(H_full) <= 4.

Thus the first-order comparison depends on at most four matrix directions, even though Sigma has many parameters. For fixed symmetric H_full and the unbiased covariance S from Gaussian rows,

    Var[tr(H_full(S-Sigma))]
       = 2 tr(H_full Sigma H_full Sigma)/(n-1).

Proof: diagonalize Sigma^(1/2) H_full Sigma^(1/2) and use independent Gaussian square moments in the Wishart representation. This is exact for the **linear functional**. Applying it to nonlinear f_jg is a delta-method approximation with a remainder.

Define V_jg = 2 tr(H_full Sigma H_full Sigma). At a regular point, a local standard error for the gap is approximately sqrt(V_jg/n). Amari's Fisher metric provides the same dual quadratic form: the leading KL distance to the linearized zero-gap boundary is f_jg^2/(2 V_jg), when V_jg>0 and the boundary is within the approximation's valid neighborhood.

This is more directly relevant than a generic covariance distance. It is also not a new safe selector yet:

- H_full and Sigma are unknown and estimated.
- The selected baseline g(D) and competitor are data-dependent, requiring simultaneous/selective control.
- The Schur remainder can grow near singular observed covariance or small candidate residual variance.
- At zero derivative, the first-order approximation fails; a higher-order expansion is needed.
- Pointwise asymptotic normality does not give a uniform finite-n guarantee near ranking ties.

Van der Vaart's excerpt, PDF pp.43-44, Theorem 3.1, explicitly requires differentiability. It justifies the smooth gap expansion at regular covariance points, not treating argmax itself as differentiable. Olive, PDF pp.44-45, likewise distinguishes an asymptotic law from an approximation whose finite-n usefulness must be assessed.

**Research implication:** a nonvacuous bound should exploit paired decision directions and control the remainder. This is a mathematical proof obligation, not permission to scan the closed grid for a promising gradient threshold.

## 6. Why a generic valid covariance region may still be too wide

Apply the Gaussian singular-value inequality in Vershynin, Corollary 5.35, to the OLS residual representation Q ~ Wishart_nu(K), with residual dimension k=m+q. Put Khat=Q/nu, t=sqrt(2 log(2/delta)), and

    h = (sqrt(k)+t)/sqrt(nu),
    a = max(0,1-h)^2, b=(1+h)^2.

With probability at least 1-delta,

    a K <= Khat <= b K

in Loewner order. When a>0, this yields Khat/b <= K <= Khat/a. If a=0, this construction supplies no finite upper covariance bound. These statements are coordinate-congruence invariant for raw, fixed-scale conditional covariance. They are not direct confidence bounds for a random sample correlation matrix.

For illustration, use p=20, k=40 and delta=0.05. These are analytic constants, not estimates or results from a new experiment:

| n | nu=n-21 | a | b | multiplicative width b/a |
|---|---|---|---|---|
| 102 | 81 | 0 | 4.0181 | no finite upper bound |
| 410 | 389 | 0.29335 | 2.12689 | 7.2504 |
| 820 | 799 | 0.46262 | 1.74198 | 3.7655 |
| 1640 | 1619 | 0.60111 | 1.49986 | 2.4952 |
| 3280 | 3259 | 0.70835 | 1.34181 | 1.8943 |

This is not an impossibility theorem. The sufficient spectral bound is conservative. Exact Wishart quantiles, lower-dimensional projections, structured covariance assumptions, or direct paired-gap concentration may improve it. However, invoking “Gaussian confidence region” does not by itself show that a 2% decision gap can be certified at n=102 or 410. The effect on G ratios must still be bounded; Loewner covariance order is not scalar score order.

This answers part of the economics question raised by NOVA-D: more Monte Carlo replicates make the evaluation of a rule more precise, but cannot increase the information in one reference sample from which the rule must choose.

External mathematical source, separate from the 19 local books: Roman Vershynin, *Introduction to the non-asymptotic analysis of random matrices*, arXiv:1011.3027, PDF p.21, Corollary 5.35: https://arxiv.org/pdf/1011.3027 . `concentration_scale.json` records the arithmetic above.

## 7. Re-encoding data cannot supply missing information, but can improve a restricted procedure

Amari's information monotonicity and sufficient-statistic discussion (PDF pp.64-66) and Watrous's relative-entropy channel monotonicity (PDF pp.288-289, Theorem 5.35) suggest a useful boundary on method proposals.

For any fixed statistic or representation h(D), every rule based on h(D) is also a rule based on D by composition. Therefore the best achievable risk using D cannot be worse than the best achievable risk using h(D), under the same loss and parameter class. This elementary decision argument is enough; a quantum theorem is unnecessary for classical metabolomics data.

This does **not** mean that features, kernels, topology or geometric encodings cannot improve a particular restricted algorithm at finite n. They can change approximation, regularization or computation. It means they do not add observations or justify superiority without specifying the restrictions they exploit. A physical quantum observation experiment would be a different statistical experiment; these data contain no such experiment.

### A concrete limitation of support topology

For O empty and scalar T, consider the Gaussian covariances for (T,X1,X2):

    Sigma_A = [[1,.6,.2],[.6,1,0],[.2,0,1]],
    Sigma_B = [[1,.2,.6],[.2,1,0],[.6,0,1]].

Both are SPD and have support R^3 with identical support homotopy type. The useful candidates differ: G_A=(.36,.04), G_B=(.04,.36). Any rule using only support topology cannot distinguish the choices. A randomized summary-only rule has worst-case regret at least .16 in this two-model example (R0=1).

This is a counterexample to support-topology sufficiency, not a refutation of every metric or weighted persistent-homology statistic. A weighted statistic could distinguish the distributions, but then its relevant information and statistical estimation still need examination. Hatcher's homotopy invariants answer topological equivalence questions; they do not encode arbitrary predictive utility.

## 8. Symmetry is valuable when it preserves the decision loss

Hall, PDF p.19, characterizes orthogonal transformations by inner-product preservation. Applied to the acquisition criterion with W=I, an orthogonal rotation Q of target coordinates sends c_j to Q c_j and leaves ||c_j||^2/v_j unchanged. An invertible transformation of O also leaves the population conditional covariance unchanged.

These provide engineering invariance tests. An algorithm should not manufacture a new scientific advantage by changing coordinates while claiming to solve the same invariant problem. Conversely, scalar ridge lambda I is not invariant under an arbitrary nonorthogonal change of observed coordinates unless the penalty is transformed with the coordinates. Whitening followed by an unchanged isotropic penalty is a different estimator, not a mere numerically stable rewrite.

For a nonorthogonal target transformation Q, preserving the original loss requires transforming W to Q^-T W Q^-1. Resetting W to I changes what counts as a useful acquisition. Lie or manifold optimization can help compute under declared constraints, but that does not supply a more favorable statistical experiment.

## 9. Identifiability of a parameter is different from identifiability of the best action

Pachter and Sturmfels' volume gives examples where different latent parameters produce the same observable distribution (PDF pp.30 and 171). Let observational equivalence be theta~theta' when P_theta=P_theta'.

A target functional phi(theta) is identifiable only if it is constant on each observational fiber. Likewise, a uniquely best acquisition is identifiable from the current experiment only if its index is constant on every such fiber. Full parameter identifiability is unnecessary if all equivalent parameters recommend the same action; conversely, parameter-estimation machinery cannot choose reliably between observationally identical states requiring different actions.

This can turn pathway ambiguity into a decision target, but only after specifying the map from physical states to observations and the loss of choosing the wrong action. Graph alignment accuracy or a smaller equivalence set is not automatically predictive-risk reduction. A genuine acquisition changes the observation experiment; re-encoding the same measured variables does not.

## 10. Biochemical constraints can add information only when they constrain the relevant variables

The Cambridge preview is Bernhard Palsson's *Systems Biology: Constraint-based Reconstruction and Analysis*. PDF p.32 explicitly states Sv=0 for the reaction-flux vector at steady state. It does not state Sx=0 for metabolite concentrations. PDF p.33 distinguishes stoichiometric null spaces, reaction pathways and conservation quantities. Cross-sectional metabolomics needs its own justified link to these dynamical constraints.

Suppose a genuinely justified linear constraint on a latent state z gives a subspace parameterization z=Nw, and let the current posterior/state covariance be C. If the prediction target is Lz and an assay observes h_j' z plus independent noise of variance sigma_j^2, Gaussian conditioning gives

    G_j = h_j' C L' W L C h_j / (h_j' C h_j + sigma_j^2).

This is the same Gaussian variance-reduction principle, with **additional scientifically justified structure and assay noise**, not a new surrogate criterion. Its possible advantage comes from valid extra information. A wrong constraint or noise model can create persistent bias and false confidence. Singular state covariance should be handled in the valid subspace; it must not be forced into an SPD formula that assumes a nonexistent inverse.

For compositional data the covariance singularity can instead be structural: centered log-ratios lie in a sum-zero subspace. The compositional-data volume, PDF p.47, contrasts clr's singular covariance with coordinate representations. Work in an appropriate fixed basis and declare zero/reference policies. If adding a measured component changes a sample-dependent normalization, the base panel and target may also change; that would confound acquisition benefit with a new observation/representation rule.

## 11. How these results change the next-step decision

The most useful outcomes of this pass are restrictions, not a claimed new winning algorithm:

1. **Do not test the common unbiased-score correction as a new selector.** It changes no picks in its valid setting.
2. **Do not promote a generic spectral confidence ball as a no-harm solution without checking its width.** The explicit example is already broad at moderate n.
3. **Study paired decision directions if pursuing a certificate.** The rank-four gradient gives a precise object. The next obligation is a finite-sample bound on its estimation and Schur remainder, not a new threshold grid.
4. **Use structure only when it supplies independent information about the acquisition target.** Flux conservation, pathway annotations, graph topology and assay normalization have different meanings.
5. **Keep the observation experiment and loss fixed when using invariances.** Coordinate-dependent shrinkage is an estimator change; target reweighting changes the objective.

A proposed successor is ready for a prospective protocol only after it has: a stated model and observation experiment, a baseline-relative risk target, a mathematically valid error bound or explicit limitation, useful constants at feasible n, an independent success/failure construction, and an implementation that preserves those assumptions. No such complete successor is established here.

The existing G1-R result remains: H1 supported at D1; H2 not established at D2; joint criterion unmet. Neither that outcome nor this mathematical review adopts A2 or reopens the closed archive.

## 12. Reading ledger for the 19 supplied PDFs

PDF page numbers refer to files in `/Users/nour/Desktop/books`. “Focused” means the listed passages were inspected; it does not mean the entire volume was read.

| ID | Source | Depth and relevant location | Role / limit |
|---|---|---|---|
| 00 | Hatcher, Algebraic Topology | Contents and PDF p.20, graph collapse/homotopy | Distinguish support topology from decision sufficiency |
| 01 | van der Vaart, Asymptotic Statistics | 47-page excerpt; focused PDF pp.43-44, 46, delta method | Smooth gap vs discontinuous selection; full book unavailable in this file |
| 02 | Bengtsson–Życzkowski, Geometry of Quantum States | Contents and PDF p.55, classical distinguishability | Divergence is not acquisition utility |
| 03 | Hall, Lie Groups, Lie Algebras, and Representations | Contents and PDF p.19, orthogonal group | Loss-preserving transformations and invariance tests |
| 04 | Palsson, Systems Biology preview | Focused PDF pp.27, 32-33 | Flux versus concentration constraints; only 56 pages supplied |
| 05 | Ortiz, Constrained Optimization on Manifolds | Introduction/outline, PDF p.16 | Numerical retractions optimize declared objectives, not statistical risk automatically |
| 06 | Woodhouse, Geometric Quantization | Image-only, not substantively read | No supporting theorem claimed |
| 07 | Golub–Van Loan, Matrix Computations, 4th ed. | Focused PDF pp.111, 143-144 in previous pass; 285-286 here | Conditional covariance, stable solves, sensitivity |
| 08 | Olive, High Dimensional Statistics | Focused PDF pp.44-45 | Finite-n usefulness of asymptotic approximations; not Vershynin |
| 09 | Amari, Information Geometry | Focused pp.37-39, 44, 174-175 previously; 64-66 here | Fisher geometry, monotonicity, sufficiency; not an acquisition theorem |
| 10 | Tsybakov, Introduction to Nonparametric Estimation | Focused PDF pp.90-99 in previous pass | Testing reductions and information limitations |
| 11 | Watrous, The Theory of Quantum Information | Focused PDF pp.288-289 | Channel monotonicity; do not infer quantum advantage for classical observations |
| 12 | Drton–Sturmfels–Sullivant, Lectures on Algebraic Statistics | Focused PDF pp.73-75, 148 previously; 149 here | Gaussian constraints, singularities and Wishart geometry |
| 13 | Shapiro–Dentcheva–Ruszczyński, Lectures on Stochastic Programming | Focused PDF pp.175, 214-215, 218 previously; 114-115 here | Validation/optimization separation and Wishart scatter |
| 14 | Pawlowsky-Glahn et al., Modeling and Analysis of Compositional Data | Focused PDF pp.23-24, 36-37 previously; 47 here | Fixed coordinates, singular clr covariance, declared relative-information model |
| 15 | Pachter–Sturmfels, Algebraic Statistics for Computational Biology | Focused PDF pp.30, 171 | Observational fibers, latent nonidentifiability and decision identifiability |
| 16 | Zhang, Complex Geometry | Contents and PDF p.40, Hodge theorem | Elliptic/Hodge hypotheses need a real observation-to-loss bridge |
| 17 | Introduction to Topological Manifolds | Image-only, not substantively read | No supporting theorem claimed |
| 18 | Nielsen–Chuang, Quantum Computation and Quantum Information | Contents/front matter indexed; data-processing sections located but not substantively read in this pass | Not used to support an acquisition claim |

No chapter is promoted as relevant merely because it is sophisticated. The results above are classical consequences or new-to-this-review derivations; their novelty relative to the literature has not been established. The external Gaussian concentration reference supplements the collection and is separately identified. Repository outputs and preregistrations were not edited or recomputed.
