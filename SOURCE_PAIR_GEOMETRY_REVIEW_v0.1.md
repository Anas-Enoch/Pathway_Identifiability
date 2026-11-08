# Source verification and decision geometry: Amari / Drton–Sturmfels–Sullivant

8 November 2026. Mathematical source review; no experiment, protocol change, or adopted method.

## Verified source copies

The two exact filenames supplied by the user are accessible locally. These files, not the earlier historical inventory, are the source copies for this supplement:

| Source | Local path | SHA-256 | PDF pages |
|---|---|---|---|
| Shun-ichi Amari, Information Geometry and Its Applications | `/Users/nour/Desktop/Information_Geometry_and_Its_Applications.pdf` | `c165e7e609499223841093bcc6b3211e14186d9dc3d5c1240392aadf60084f69` | 378 |
| Mathias Drton, Bernd Sturmfels, Seth Sullivant, Lectures on Algebraic Statistics | `/Users/nour/Desktop/Lectures on Algebraic Statistics.pdf` | `5a5decab942b69f49cefb6d88ee304fa10f7b225cda5b7caaed5a2ed2a55bcac` | 179 |

The first is the local corrected Amari volume. The second is a text dated September 27, 2008. Page numbers below identify these PDFs. The local-copy match establishes access to these exact named books; it does not establish access to every ChatGPT project attachment.

Passages directly inspected for this supplement: Amari PDF pp.37–39, 44, 174–175; Drton et al. PDF pp.73–75 and 148. Earlier inspection included Amari pp.15–16, 139. Other pages were searched/indexed, not claimed read cover to cover.

## 1. Amari: the appropriate geometry depends on the object being compared

Amari's generalized Pythagorean relations and projection results on PDF pp.37–39 concern divergence, dual coordinates, orthogonality, and model projections. PDF p.44 connects exponential-family geometry, KL divergence, and Fisher information. PDF pp.174–175 discusses estimation efficiency and coordinate dependence.

These are relevant to covariance uncertainty, but do not show that a divergence-optimal estimate minimizes acquisition regret. Projection claims need their model and regularity assumptions; global uniqueness or global optimization cannot be inferred merely from a local orthogonality condition. Likewise, a Cramér–Rao statement concerns the stipulated estimation setting and regularity, not the error probability of a discontinuous argmax decision.

Programme connection: G1-R has already shown why a more accurate or less biased population ranking rule need not satisfy the finite-sample no-harm criterion. Replacing Euclidean covariance distance by KL distance alone does not close that gap.

## 2. Drton et al.: covariance constraints must be interpreted inside the positive-definite domain

PDF p.73 derives Gaussian conditional independence from a zero Schur-complement block. PDF pp.74–75 describes the associated polynomial constraints and graphical models. A particularly useful warning on p.74 is that algebraic components can exist entirely outside the positive-definite covariance cone; allowing singular matrices changes some conditional-independence implications. PDF p.148 gives an example of nonstandard likelihood asymptotics at a model singularity.

Programme connection: a proposed algebraic reduction must preserve positive definiteness, nonzero denominators and the actual conditioning set. A graph encoding biochemical reactions is not automatically a valid conditional-independence graph for measured concentrations. Using its constraints to reduce variance needs a scientific observation model or independently justified restrictions.

An acquisition tie is not automatically a singularity of a statistical covariance model. The ambient Gaussian model can be regular at an SPD matrix even though its best-candidate map changes discontinuously there. Those are distinct phenomena; the source's singular-model likelihood results cannot be imported just because a ranking has a tie.

## 3. Derived connection: acquisition ranking boundaries are explicit algebraic objects

This section is a derivation for the programme, not a theorem quoted from either book and not a novelty claim.

In conditional-covariance coordinates K = Cov((T,P) | O), let c_j = K_Tj, v_j = K_jj > 0, and fix target metric W >= 0. Then

    G_j = c_j' W c_j / v_j.

The sign of G_j-G_g equals the sign of

    F_jg(K) = v_g c_j' W c_j - v_j c_g' W c_g.

Thus a pairwise ranking boundary is F_jg = 0 inside the SPD domain. For fixed W, F_jg is a cubic polynomial in these conditional entries. In full-covariance coordinates the Schur map is rational; determinant factors can be cleared only while preserving their known positive signs. A decision region is described by several such inequalities, not merely by a covariance-estimation ball.

The potentially useful geometric object is distance to a relevant decision boundary, rather than distance to a generic covariance target. For a fixed reference g and proposed j, a population diagnostic could be

    rho_jg(Sigma) = inf KL(N(0,Sigma) || N(0,Sigma'))
                   over SPD Sigma' with G_j(Sigma') <= G_g(Sigma').

This formula defines a diagnostic; it does not supply a global optimizer, confidence coverage, or a safe learning algorithm. It should not be evaluated on the closed archive to tune a successor.

For a small symmetric perturbation E of an SPD covariance,

    KL(N(0,Sigma) || N(0,Sigma+E))
       = (1/4) tr[(Sigma^-1 E)^2] + O(||E||^3),

with a locally controlled remainder away from singular covariance. The Fisher metric is therefore (1/2)tr(Sigma^-1 E Sigma^-1 E). The differential of f = G_j-G_g, projected onto a specified regular model's tangent space, gives a local decision sensitivity. This is an asymptotic/local approximation; near singular covariance, vanishing gradients, or poorly controlled remainders, it is not a finite-sample certificate.

## 4. What would make this different from failed NOVA certification?

A different metric or a new name would not suffice. A successor would need all of the following:

1. The comparison is against the actual data-dependent frozen baseline pick, not a substituted population oracle.
2. The covariance region is simultaneous and statistically calibrated under the stated observation model, including sample standardization and estimation of O-conditioning.
3. The computed bound is valid in the conservative direction even if the boundary-distance problem is solved approximately.
4. It bounds the relevant expected excess regret, or explicitly states a different safety target; a 95% event alone does not meet a tiny relative expected-regret margin.
5. It retains a useful set of improvement opportunities at feasible n rather than returning G1 almost everywhere.

These are unresolved obligations, not findings that the method works. They distinguish useful mathematics from repeating a gated-estimator experiment with a new geometry.

## 5. Current decision

The two books support a sharper mathematical question: can the geometry and algebra of acquisition decision boundaries produce a calibrated, informative bound on baseline-relative regret? They do not establish a new improved selector.

The strongest deliverable before another study would be either (a) such a bound with useful constants and explicit assumptions, or (b) an information/uncertainty calculation showing that the desired protection and improvement cannot both be certified at realistic sample sizes by the proposed construction.

G1-R's joint criterion remains unmet. H2 remains not established, not refuted. A0 remains the comparator; A2 is not adopted; the archive stays closed. No frozen files or study outputs were modified in preparing this supplement.
