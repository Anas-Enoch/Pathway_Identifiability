# Nonlinear predictive value and the limits of certified metabolite acquisition under partial observation

Draft v0.5, 4 Nov 2026, archive version (v0.5 replaces the §7 Limitations paragraph with the archive-closure, N2 regime-diagnostic and NHANES-FA findings of 3–4 Nov 2026; no other change from v0.4. v0.4 of 20 Oct 2026, revision after the estimand audit: §2 gains an estimand map and the definition of the frozen normalized endpoint Θ; §3 gains Proposition 3, an exact weighting identity and a delta-method representation connecting the paired-loss difference to Θ; §5.3 adds the representation hierarchy L0; §6.3 adds the N2 replay with its diagnostic recomputed on the primary's normalization; abstract and §7 adjusted so that no interval for Θ is attributed to the unnormalized score theorem). Every number is taken from a frozen record (hash table in S1); the replay recomputations are from the frozen block-fold file through the audit record. Status: for the owner's review; not frozen.

## Abstract

Choosing which metabolite to measure next is usually decided by a covariance criterion. For the one-step Gaussian
criterion G1, we ask whether nonlinear conditional dependence that covariance cannot see can be turned into better
acquisition. We separate three quantities that are often conflated: the population linear acquisition value G1_j,
the population Bayes acquisition value V_j (the reduction in squared-error Bayes risk from measuring X_j), and the
fitted-procedure advantage D^proc_{j,g} (the held-out risk difference between the predictor we would deploy with
candidate j and with G1's choice g). V_j = G1_j under joint Gaussianity, and V_j − G1_j = A_j − A_0, the change
in the linear-versus-Bayes approximation gap when X_j is added; a candidate can improve linear prediction greatly
while adding nothing to the Bayes predictor, so nonlinear prediction gains need not become acquisition gains. We
then study selectors that depart from G1 only when a cross-fitted paired risk-difference test supports the
departure. In synthetic data with exact population values, such selectors never departed wrongly (0 of 560
selector-blocks), recovered G1's choice under Gaussianity, and detected planted nonlinear value only when it
dominated. Simultaneous multiplier-bootstrap bounds exhibited substantially better coverage for fitted-procedure
acquisition advantages than for population Bayes advantages (0.90 against 0.70 on the hardest regime, nominal
0.95); the difference is consistent with a substantial nuisance-estimation contribution, while the residual
procedure-level shortfall shows that targeting fitted-procedure risk does not by itself establish valid
selective inference. On a 186-study public metabolomics
archive, after a regularization-parity audit removed half of a kernel predictor's apparent edge, the selective
rule gave a development acquisition gain over G1 of +0.11% of loss (equal-study mean of within-study relative
loss reduction; 97.5% study-bootstrap bounds +0.03% to +0.22%), below the prespecified 0.29% qualification
threshold, while unrestricted nonlinear selection was worse than G1. The normalized endpoint is connected to the
paired-loss difference of the theory by an exact weighting identity and a delta-method representation; its
interval rests on the study bootstrap, not on the score theorem. A replay of the frozen selector shows that the
realized gain of a departure increases with the pre-selection evidence, and that no threshold on that evidence
would have reached the qualification floor (hindsight maximum +0.16%). No external test was run; the reserved
evaluation pool remains sealed.

## 1. Introduction

Paper IV established that G1 and the information criterion RIG score the same rank-one covariance update under
different metrics, that strict disagreement between them certifies Loewner incomparability across task
weightings, and that on 186 public metabolomics studies the prespecified RIG-superiority hypothesis was not
supported. Under the Gaussian conditional model and isotropic squared loss, G1 is the population one-step optimum.
Any selector that beats G1 on G1's own loss must therefore exploit non-Gaussian conditional dependence or a
finite-sample estimation advantage.

A kernel predictor on the archive predicts the targets about 9% better than a log-linear ridge under a bounded
loss (Gate 3A), which invites the inference that nonlinear structure is present and exploitable for acquisition.
This paper shows why that inference fails in general and quantifies what remains. The paper is a decision-theoretic analysis of why nonlinear predictive improvements fail to translate reliably into measurement-acquisition improvements, supported by a reproducible computational investigation; it does not propose a general theory of valid certified acquisition. Three questions are kept apart:
which measurement most improves population linear prediction (G1_j); which most reduces Bayes risk (V_j); and which
improves the predictive procedure we have actually fitted (D^proc). They have different answers, different
estimators and different certification problems.

## 2. Framework

Targets T ∈ ℝ^q, base panel O, candidates X_1, …, X_M. For a set S of inputs let R_L(S) = inf over affine f of
E‖T − f(S)‖² and R_H(S) = inf over f ∈ L²(σ(S)) of E‖T − f(S)‖². Define the linear acquisition value
G1_j = R_L(O) − R_L(O ∪ {X_j}) and the Bayes acquisition value V_j = R_H(O) − R_H(O ∪ {X_j}). For a
fitted procedure producing predictors m̂_g, m̂_j on independent training data, the fitted-procedure advantage is
D^proc_{j,g} = E_new[ℓ(T, m̂_g(O, X_g)) − ℓ(T, m̂_j(O, X_j)) | D_train], the expectation over fresh data
conditional on the training sample. Under K-fold cross-fitting with training sets D_{train,k} and held-out
fractions w_k, the quantity the cross-fitted score estimates is the fold-averaged advantage
D^CF_{j,g} = Σ_k w_k E_new[ℓ(T, m̂_{g,k}) − ℓ(T, m̂_{j,k}) | D_{train,k}], which is the advantage of the
procedure as applied to samples of size (K − 1)n/K, not of a single predictor refitted on all n rows
(Proposition 2, §3).

Paper IV's G1_j = c_jᵀc_j / v_j, with c_j and v_j the conditional target–candidate covariance and candidate
variance given O computed by Schur complements, equals R_L(O) − R_L(O ∪ {X_j}) (verified to four decimals on the
counterexample of §3.3), so the results below apply to G1 as implemented.

*Estimand map.* Four summaries of the same paired held-out losses appear in this paper and must not be read as one.
With L_{j,kb} and L_{g,kb} the held-out loss sums of the candidate arm and the G1 arm in block-fold b of study k,
μ_{j,k} = Σ_b L_{j,kb} / n_k the mean per-cell loss in study k, δ_k = μ_{g,k} − μ_{j,k} and b_k = μ_{g,k} > 0:

| Object | Definition | Role |
|---|---|---|
| Paired-loss difference | D = E[L_g − L_j] per cell; D^proc, D^CF of this section | the object of Lemma 1 and Theorems 1–3 |
| Pooled relative improvement | Θ_pool = 1 − Σ_k μ_{j,k} / Σ_k μ_{g,k} = Σ_k δ_k / Σ_k b_k | reported for comparison only |
| Equal-study mean of within-study relative improvement | Θ = (1/K) Σ_k (1 − μ_{j,k}/μ_{g,k}) = (1/K) Σ_k δ_k / b_k, K = 186 | the frozen primary Δ_acq of N2, Gate 3A and N0 |
| Mean of block-fold ratios | E_b[1 − L_{j,kb}/L_{g,kb}] | an audit diagnostic only (§6.3); dominated by small-denominator folds |

The equal-study weighting of Θ is part of the scientific endpoint: each study counts once regardless of its loss
scale. Θ is the frozen primary; its inference is the study-level bootstrap stated in §4, and its relation to D is
Proposition 3. Two targets must also be kept apart for Θ: the observed collection of 186 studies, for which Θ is a
finite-population quantity and the bootstrap describes resampling variability of the collection, and a
hypothetical population of future studies from which the collection is treated as a sample. The archive is a
heterogeneous convenience collection, so the second reading is an assumption, not a fact of the design.

## 3. Theoretical results

Throughout, T ∈ L², and for a set S of inputs R_L(S) and R_H(S) are as in §2; P_{L_S} and P_{H_S} denote the
L² projections of T onto the affine span of S and onto L²(σ(S)). Since the affine span is a closed subspace of
L²(σ(S)), A(S) = R_L(S) − R_H(S) = E‖P_{H_S}T − P_{L_S}T‖² ≥ 0 by Pythagoras.

**Theorem 1 (conditional risk-reduction identity).** Let m_0 = E[T | O] and m_j = E[T | O, X_j]. Then
V_j = R_H(O) − R_H(O ∪ {X_j}) = E‖m_j − m_0‖².

*Proof.* R_H(O) = E‖T − m_0‖² and R_H(O ∪ {X_j}) = E‖T − m_j‖². Write T − m_0 = (T − m_j) + (m_j − m_0). Then
E‖T − m_0‖² = E‖T − m_j‖² + E‖m_j − m_0‖² + 2E[(T − m_j)ᵀ(m_j − m_0)]. The cross term is
2E[E[T − m_j | O, X_j]ᵀ(m_j − m_0)] = 0 because m_j − m_0 is σ(O, X_j)-measurable and E[T − m_j | O, X_j] = 0.
Subtracting gives the claim. □

**Proposition 1 (Gaussian specialization).** If (T, O, X_j) is jointly Gaussian with A = Cov(T | O),
c_j = Cov(T, X_j | O) and v_j = Var(X_j | O) > 0, then m_j − m_0 = (c_j / v_j)(X_j − E[X_j | O]) and
V_j = c_jᵀc_j / v_j = G1_j.

*Proof.* Under joint Gaussianity, conditional expectations are affine, so E[T | O, X_j] − E[T | O] is the
regression of T − E[T | O] on X_j − E[X_j | O], with coefficient c_j / v_j. Then
E‖(c_j / v_j)(X_j − E[X_j | O])‖² = (c_jᵀc_j / v_j²) · v_j. □ Numerically, max_j |V_j − G1_j| = 0.003 on the
Gaussian synthetic design (X0).

**Theorem 2 (signed approximation-gap decomposition).** V_j − G1_j = A(O ∪ {X_j}) − A(O).

*Proof.* V_j − G1_j = [R_H(O) − R_H(O ∪ {X_j})] − [R_L(O) − R_L(O ∪ {X_j})]
= [R_L(O ∪ {X_j}) − R_H(O ∪ {X_j})] − [R_L(O) − R_H(O)] = A(O ∪ {X_j}) − A(O). □

The difference is signed because A(O) need not vanish: a candidate can close the linear model's approximation gap
without adding information to the Bayes predictor. In the construction T = O_1² + ε, X_1 = O_1² + noise (X0b),
m_0 = O_1² already, so V_1 = 0, while A(O) = 6.07, A(O ∪ {X_1}) = 0.26 and G1_1 = 5.81; the Schur-complement
G1 of Paper IV (5.8101) equals R_L(O) − R_L(O ∪ {X_1}) (5.8101), so the decomposition applies to G1 as
implemented. Two consequences. (i) Nonlinear prediction gain, which is A(O ∪ {X_j}) − A(O) summed with the
linear gain, is not acquisition gain. (ii) The signed decomposition invalidates screening arguments that assume
nonlinear excess acquisition value is nonnegative. Valid screening remains possible: with g = argmax_j G1_j and
M_{j,g} = G1_g − G1_j ≥ 0, V_j − V_g = (A(O ∪ {X_j}) − A(O ∪ {X_g})) − M_{j,g}, so on the coverage event of
simultaneous bounds U(A_j − A_g) < L(M_{j,g}) candidate j cannot outperform g, with no sign assumption on
A_j − A_g; constructing sufficiently informative bounds is an unresolved statistical problem.

**Theorem 3 (conditional safety of certified departure).** Let D_{j,g} denote either the Bayes advantage
V_j − V_g or the fitted-procedure advantage D^CF_{j,g}, and let L_{j,g} be lower bounds.
(a) *Fixed reference.* If g is prespecified and Pr(D_{j,g} ≥ L_{j,g} for all j ≠ g) ≥ 1 − α, then the rule
π that departs to argmax_j L_{j,g} iff max_j L_{j,g} > τ satisfies Pr(π ≠ g and D_{π,g} ≤ τ) ≤ α.
(b) *Data-selected reference.* If ĝ = ĝ(data) and the bounds satisfy Pr(D_{j,ĝ} ≥ L_{j,ĝ} for all j ≠ ĝ)
≥ 1 − α, where the probability is over the joint law of the data including the selection of ĝ, the same
conclusion holds with ĝ in place of g. Coverage established for a fixed g does not imply the premise of (b).

*Proof.* On the coverage event, π ≠ g implies L_{π,g} > τ and hence D_{π,g} ≥ L_{π,g} > τ; a departure with
D_{π,g} ≤ τ therefore lies in the complement of the coverage event, whose probability is at most α. The
argument is identical in (b) once the coverage event is defined for the selected pair. □

The premise of (b) is the substance of the result. In the experiments of §5 the reference is the plug-in G1 pick
of the same block, so the measured coverage is of type (b) empirically, but no construction in this paper is
shown to satisfy (b).

**Lemma 1 (orthogonal-score identity and nuisance remainder).** With d̂ = m̂_j − m̂_0 and fixed functions
m̂_0, m̂_j of the appropriate arguments, ψ̂ = 2d̂ᵀ(T − m̂_0) − ‖d̂‖² = ‖T − m̂_0‖² − ‖T − m̂_j‖² pointwise,
and E[ψ̂] = V_j + ‖m̂_0 − m_0‖²_{L²} − ‖m̂_j − m_j‖²_{L²}.

*Proof.* The pointwise identity is ‖T − m̂_0‖² − ‖T − m̂_0 − d̂‖². For the expectation, write e_0 = m̂_0 − m_0,
e_j = m̂_j − m_j and d = m_j − m_0, so d̂ = d + e_j − e_0 and T − m̂_0 = (T − m_j) + d − e_0. Expanding
E[ψ̂] and using E[T − m_j | O, X_j] = 0 leaves V_j − 2E[dᵀe_0] + ‖e_0‖² − ‖e_j‖²; the term E[dᵀe_0] vanishes
because e_0 is σ(O)-measurable and E[d | O] = E[m_j | O] − m_0 = 0. □

The cross-fitted loss difference used for selection is therefore this score; its bias is the second-order
remainder B_{j,g} = ‖m̂_g − m_g‖² − ‖m̂_j − m_j‖² for the paired comparison, which cross-fitting does not
remove and a bootstrap centred on the estimated score cannot see.

**Proposition 2 (fold-averaged fitted-procedure estimand).** Under K-fold cross-fitting with i.i.d. rows, let
D̂_k be the mean paired score over the held-out rows of fold k and D_k^proc = E_new[ℓ(T, m̂_{g,k}) − ℓ(T, m̂_{j,k})
| D_{−k}], with D_{−k} the training sample excluding fold k. Then E[D̂_k | D_{−k}] = D_k^proc for each k, and
E[Σ_k w_k D̂_k] = E[D^CF_{j,g}] with D^CF as in §2. The cross-fitted score therefore estimates the advantage of the
procedure as applied to samples of size (K − 1)n/K, not of predictors refitted on all n rows.

*Proof.* The held-out rows of fold k are independent of D_{−k}, so each paired score has conditional expectation
D_k^proc given D_{−k}, and so does their mean. Taking expectations and averaging with the weights w_k gives the
second claim. □

The statement is fold-specific by necessity. A held-out row of fold k belongs to the training sets of the other
folds, so conditioning on all K training sets conditions on the whole sample and leaves nothing random; the
conditional unbiasedness holds for each fold given its own complement, and only the unconditional identity holds
for the average. The variance of the average involves the dependence between folds induced by overlapping
training sets, which is one of the candidate sources of the coverage shortfall in §5.

**Proposition 3 (from the paired-loss difference to the normalized endpoint).** With θ_k = δ_k / b_k as in §2:
(i) *Exact weighting identity.* Θ_pool − Θ = Cov_K(b_k, θ_k) / b̄, where Cov_K uses the 1/K convention and
b̄ = (1/K) Σ_k b_k.
(ii) *Delta-method representation.* Let (δ̂_k, b̂_k) estimate (δ_k, b_k) with b_k bounded away from zero, and
suppose the joint estimation error is asymptotically linear with influence functions ψ_{δ,k}, ψ_{b,k}. Then
θ̂_k − θ_k = (δ̂_k − δ_k)/b_k − (δ_k / b_k²)(b̂_k − b_k) + r_k with r_k second order, so the normalized score has
influence function ψ_{θ,k} = ψ_{δ,k}/b_k − (δ_k / b_k²) ψ_{b,k}, and Θ̂ − Θ = (1/K) Σ_k (θ̂_k − θ_k) inherits the
average.

*Proof.* (i) Θ_pool = Σ_k b_k θ_k / Σ_k b_k = (1/K) Σ_k b_k θ_k / b̄, and (1/K) Σ_k b_k θ_k = b̄ Θ + Cov_K(b, θ).
(ii) is the first-order Taylor expansion of h(δ, b) = δ/b, whose gradient is (1/b, −δ/b²), applied to the
asymptotically linear pair. □

Part (i) is exact and explains the discrepancy between the two normalizations on the archive: Θ = +0.1118% and
Θ_pool = +0.1405% for the primary selector, and the difference, 0.0287 percentage points, equals Cov_K(b, θ)/b̄
to the last digit on the frozen block-fold file. The pooled figure weights studies by their G1 loss; the
difference is a weighting choice, not noise. Part (ii) is concrete in the simplest case. For fixed trained
predictors m̂_g, m̂_j and an independent evaluation observation W = (O, X, T), write Y_g(W) = ℓ(T, m̂_g(O, X_g))
and Y_j(W) = ℓ(T, m̂_j(O, X_j)), so that b = E[Y_g] and δ = E[Y_g − Y_j] are means of fixed functions with
influence functions ψ_b(W) = Y_g(W) − b and ψ_δ(W) = Y_g(W) − Y_j(W) − δ, the latter being the orthogonal score
of Lemma 1 centred. Then

ψ_θ(W) = (Y_g(W) − Y_j(W) − δ)/b − (δ/b²)(Y_g(W) − b) = (Y_g(W) − Y_j(W) − θ Y_g(W)) / b,

the influence function of the normalized fitted-procedure risk improvement under independent evaluation and
fixed trained predictors. The expression does not transfer automatically to the full procedure, in which the
candidate is data-selected, the training sets of different folds overlap and the normalization is estimated on
the same held-out rows; that extension needs arguments not given here. The denominator condition must also be
uniform for uniform inference: b > 0 for the realized study is not a population statement, and the condition
required is inf_{P ∈ 𝒫} b(P) ≥ b_min > 0 over the model class, with the study bootstrap exposed to whether the
observed within-study G1 denominators are stable under resampling of small studies. Neither is verified
empirically here; both are assumptions. Part (ii) therefore does not by itself validate the study bootstrap used
for the interval of Θ: that requires independence across studies (or a stated dependence structure), consistent
within-study estimation of (δ_k, b_k), the uniform denominator condition, and an asymptotic regime in which the
second-order remainders are negligible against the study-level variability. In this paper the interval for
Θ is a study-bootstrap interval under those assumptions; it is not a consequence of Lemma 1 or Theorem 3, which
are stated for D.

## 4. Methods

*Selectors.* G1: Paper IV's plug-in rule (linear shrinkage, Schur complements), re-estimated on each training
split. NOVA-V: argmax of the cross-fitted score V̂_j. NOVA-U: argmax of V̂_j − se_j. NOVA-S: NOVA-V's pick if the
paired per-row difference against G1's pick exceeds 2 standard errors, else G1's pick. NOVA-X: pairwise
D̂_{j,g} for j ≠ g with a one-sided simultaneous critical value from a Gaussian multiplier bootstrap (B = 2,000)
of max_j √n(D̂_j − D_j)/s_j; depart to argmax_j L_{j,g} iff max_j L_{j,g} > 0, α = 0.05.

*Nuisances.* RBF kernel ridge on log-rule-standardized inputs (bandwidth {0.5, 1, 2} × median distance, ridge
grid up to 10⁴), hyperparameters by inner cross-validation per model; two alternatives were tried and retained as
ablations, an additive candidate-effect model and a residualized kernel Q_O K_j Q_O.

*Cross-fitting.* Five outer folds over subjects (seed 20261019); selection on the training folds by inner
cross-fitted scores; scoring on the held-out fold with the same predictor family in both arms, so that
Δ_acq = 1 − R(O ∪ {ĵ}) / R(O ∪ {j_G1}) isolates the acquisition decision from predictor quality.

*Loss.* Squared error in design-standardized units with a scale floor at 1% of the block's median target SD;
the bounded loss min(e², 25) of Gate 3A as sensitivity.

*Regularization parity (N0).* Prediction on the base panel only, six families, including the log-linear ridge at
the Gate 3A grid (λ ≤ 10) and at the wide grid, a linear-kernel control, and the kernel ridge at both grids.

## 5. Synthetic experiments

Designs with known conditional means (Gaussian, planted square at two strengths, interaction, redundant
candidates, n = 30, covariate shift, compositional closure), 5 base features, 3 targets, 8 candidates, n = 150;
population V_j by Monte Carlo on the true functions (compositional: an oracle kernel estimate, approximate).
Targets standardized to their population SD so that V_j and G1_j share the isotropic objective.

*Calibration (N1, 320 blocks).* V_j = G1_j to three decimals where the model is Gaussian; NOVA-S's regret equals
G1's on Gaussian data with no departures; one spurious departure in 320 blocks (small-sample regime). Where a
planted square or interaction dominates the candidate set, G1 misses it (regret 0.64, 0.06–0.13) and NOVA-S finds
it (0.00). With the weak square (V_1 = 0.17 against a linear best of 0.19) no nuisance detects it at n = 150,
NOVA-V is worse than G1, and NOVA-S defers.

*Certification (X1, 480 selector-blocks).* No false departure for NOVA-S or NOVA-X. NOVA-X was never more
accurate than NOVA-S and was worse on the planted square (regret 0.106 vs 0.035) and under covariate shift (0.541
vs 0.386): the simultaneous critical value (q ≈ 2.2) costs power against a false-departure rate that was already
zero. Simultaneous coverage of the Bayes advantage: 0.80–1.00, below nominal on the interaction regime.

*Which estimand is covered (X1b).* Evaluating the fold-fitted models on fresh samples gives the fitted-procedure
advantage. Simultaneous multiplier-bootstrap bounds exhibited substantially better coverage for fitted-procedure
acquisition advantages than for population Bayes advantages: at 100 blocks (S2), 0.95 against 0.94 on the planted
square and 0.90 against 0.70 on the interaction regime (nominal 0.95). Nevertheless, coverage remained below
nominal in the interaction regime (0.90 versus 0.95), demonstrating that targeting fitted-procedure risk does
not, by itself, establish valid selective inference. [Certain] the fitted-procedure interval undercovers there;
[Likely] the procedure-versus-Bayes difference is substantially attributable to nuisance-estimation bias (Lemma
1); [Unresolved] the remaining procedure-level shortfall may arise from selection of the G1 reference,
heavy-tailed paired scores, dependence induced by overlapping cross-fitting folds, or their interaction. Zero
false departures in 200 further blocks is evidence of conservative observed behaviour, not proof of finite-sample
false-departure control. The residualized
kernel matched the RBF nuisance in regret and had slightly lower coverage; the additive model was biased downward
and worse everywhere.

*Representation hierarchy (L0).* Theorem 2 compares two levels, affine and Bayes. An intermediate level makes the
decomposition multilevel: for nested function classes F ⊆ G on each input set, V^G_j − V^F_j = A_{F,G}(O ∪ {X_j})
− A_{F,G}(O) with A_{F,G}(S) = R_F(S) − R_G(S) ≥ 0, and the gaps telescope, A_{L,H} = A_{L,K} + A_{K,H}. Two
qualifications are needed before a kernel class can serve as the middle level, and both were found by
construction rather than assumed. First, a bare Gaussian-RBF RKHS is nested in neither sense: affine functions
are not in it (L_S ⊄ K_S) and a function of O alone is not in the RKHS on (O, X_j) (K_O ⊄ K_{O∪j}); with that
class the K-level value V^K_j was negative for most candidates in every design while the subtraction identity
still held to 10⁻¹⁷, so the sign structure, not the identity, is what nesting buys. The cumulative kernel
k_S = 1 + zᵀz′ + rbf(O) + [rbf(O, x_j) if j ∈ S] restores both inclusions, and with it both gaps are nonnegative
by Pythagoras. Second, the computed K-level is a Tikhonov-regularized risk at fixed λ, not the projection risk
onto the RKHS; the theoretical hierarchy (projections, exact) and the computational hierarchy (regularized fits,
approximate) must be kept distinct, and only the former satisfies the identities above exactly.

On the synthetic designs (population values; nested class; λ = 0.01): the three levels coincide under Gaussianity
(V^L = V^K = V^H to two decimals, Proposition 1 at every level); the nonlinear candidate's value rises
monotonically through the levels on the planted square (0.08 → 0.53 → 0.92) and the interaction (0.24 → 0.42 →
0.77); on the linear-proxy design the uninformative candidate is valued 5.73 at L, 2.39 at K and 0 at H, and every
level's ordering buys it. The K-level is not a score: on the linear-proxy design the seven pure-noise candidates
have V^L = V^H = 0 but V^K = 0.80 each, because adding the summand rbf(O, x_j) enlarges the class's capacity to
represent E[T | O] = O_1² whether or not x_j carries information. V^K_j is a class-dependent accessibility value
with a capacity component; only V^H measures information and only V^L is the covariance criterion. The K-level
locates where between linear and Bayes a candidate's apparent value arises, and that is its only use here. The
fitted cross-fitted value at n = 150 with a tuned penalty tracks V^H rather than any fixed-λ V^K (square 0.90
against V^H 0.92 and V^K 0.53), as a tuned procedure is a far less restricted class than a fixed-λ population
class.

## 6. Development archive

The 186-study Metabolomics Workbench archive of Paper IV, design rows only, frozen blocks (20 base, 20 targets,
20 candidates), Gate 3A's folds. The archive had been used by Paper III, Paper IV, a retrospective predictor
benchmark, Gate 2 and Gate 3A before these analyses; every result here is a development result.

*Regularization parity (N0).* The kernel predictor's +9.14% bounded-loss edge over the log-linear ridge (Gate
3A) reproduces exactly; the ridge had chosen λ = 10, the grid edge, in 92% of folds. Against the ridge at the wide
grid the edge is +4.25% squared [2.30, 5.93] and +4.50% bounded; the linear-kernel control equals the ridge to
10⁻¹⁵. Half of the apparent nonlinear prediction gain was the regularization ceiling.

*Acquisition (N2; pre-registered primary NOVA-S, RBF, floored squared loss).* Δ_acq = +0.112%, 97.5% bounds
[+0.029%, +0.221%], 186 records; qualification required ≥ 0.29% and was not met. NOVA-V −0.282% [−0.724, +0.066];
NOVA-U −0.209% [−0.846, +0.247]; additive nuisance −0.9% for both and −0.079% for NOVA-S. NOVA-S departed from
G1 in 10.5% of block-folds (null self-test rate 5.0%). The bounded-loss sensitivity agrees (+0.135% [+0.056,
+0.243]). Gate 3A's RBF family, which is NOVA-V under the bounded loss and the narrow grid, had tied G1 (−0.108%).

The narrow conclusion: selective departures from G1 produced a small positive held-out acquisition contrast on a
repeatedly used development archive; unrestricted nonlinear acquisition did not improve on G1. With the reserved
pool's detectable effect of 0.29%, a true effect of 0.11% would be detected with power about 0.29; no external
test was run and the pool remains sealed.

*Replay of the frozen selector (exploratory).* An instrumented twin of the frozen N2 run, gated on exact
reproduction of every block (952 of 952; the summary, claim rate and per-study RBF columns to the last digit),
logged for each of the 4,760 block-folds the inner paired statistic t of the NOVA-V pick against the G1 pick and
the held-out loss of both arms. In 3,062 block-folds the two picks differed; NOVA-S departs when t > 2. The table
gives, by t bin fixed before the run, the pooled within-bin relative loss improvement 1 − Σ L_jV / Σ L_G1 of the
NOVA-V pick over the G1 pick (study-clustered 95% intervals, 2,000 resamples) and the mean per-cell paired-loss
difference. The pooled within-bin ratio shares the primary's normalization form but not its equal-study
weighting; it answers a related question, not the primary's, and is labelled accordingly. The equal-study
within-bin version is given for completeness.

| t bin | block-folds (studies) | pooled relative improvement | mean per-cell L_G1 − L_jV | equal-study within-bin | share > 0 |
|---|---|---|---|---|---|
| [0, 1) | 1,482 (174) | −0.27% [−0.59, +0.12] | −0.0037 [−0.0056, −0.0017] | −0.89% | 0.45 |
| [1, 2) | 1,078 (150) | +0.37% [−0.04, +0.81] | +0.0015 [−0.0002, +0.0030] | −0.01% | 0.52 |
| [2, 3) departures | 394 (109) | +0.89% [+0.34, +1.60] | +0.0030 [+0.0003, +0.0057] | +0.60% | 0.62 |
| [3, 4) departures | 95 (46) | +0.87% [−0.08, +1.64] | +0.0022 [−0.0044, +0.0072] | +0.96% | 0.72 |
| ≥ 4 departures | 13 (8) | +5.69% [−0.48, +6.54] | +0.0259 [−0.0024, +0.0387] | +2.46% | 0.92 |

Spearman(t, per-cell difference) = 0.160; Spearman(t, block-fold ratio) = 0.155, study-clustered 95% [0.106,
0.203]. The first run of this diagnostic summarized the mean of block-fold ratios, which is dominated by folds with
small G1 loss and read the [2,3) bin as zero; the audit replaced it (S1, audit record) and the original table is
retained there as the record of what was run. The contributions of the departure bins to the primary, computed on
the primary's own estimand, are +0.086% for 2 ≤ t < 3, +0.020% for 3 ≤ t < 4 and +0.006% for t ≥ 4.

Hindsight envelope: the primary Θ that the frozen pipeline would have produced had NOVA-S departed at t > τ, from
the same block-fold losses: τ = 0, −0.280%; 0.5, −0.007%; 1, +0.123%; 1.5, +0.158%; 2 (frozen), +0.112%; 2.5,
+0.052%; 3, +0.025%; 4, +0.006%. The maximum over thresholds, +0.158%, is a hindsight quantity on data used eight
times and lies below the 0.29% qualification floor.

What the replay supports: the realized improvement of a departure increases broadly with the pre-selection
evidence, the [2,3) bin is positive under the pooled ratio, the mean difference and the equal-study summary, the
[0,1) bin is negative in mean with a pooled-ratio interval that includes zero, the top bin rests on 13 block-folds,
and no threshold on t reaches the qualification requirement. The gate removes decisions that would have been
harmful on average (bins below t = 2 under the equal-study summary) more than it finds large gains, which is why
NOVA-V alone is −0.28% and NOVA-S is +0.11%.

## 7. Discussion

Three points follow from keeping the estimands apart. First, prediction gain is not acquisition gain: the
decomposition shows that a candidate's effect on the linear model's approximation gap can dominate its Bayes
information, and the archive shows a 4–9% prediction edge coexisting with a tie in acquisition. Second, the
certification target matters: the same bootstrap bounds cover the fitted-procedure advantage far better than the
Bayes advantage (0.90 against 0.70 on the hardest regime), consistent with the nuisance remainder entering one
and not the other, and still fall short of nominal for the fitted-procedure advantage, so neither estimand is
certified by the present construction. Third, the gate, not
the estimator, produced the archive gain: the selective rules won by abstaining in 90% of decisions, and the
replay shows the abstentions fell mostly on decisions that were harmful on average. Fourth, the endpoint must be
named: the theory is stated for a paired-loss difference, the qualification rule for an equal-study mean of
within-study ratios, and the two are related by Proposition 3, under which the pooled and equal-study
normalizations differ by an exact covariance term. A manuscript that reports one interval for one endpoint
without saying which object it is for would have conflated them; this one says which inferential claim belongs
to which object.

**Limitations.** Monte Carlo resolution for coverage (±0.03 at 100 blocks, ±0.07 at 20); selection of g treated as
fixed in the bounds while the empirical coverage is of the selected pair; the cross-fitted estimand is an average
over fold-specific procedures. The 186-record archive has been used for method development eleven times, including
by every analysis reported here; nothing it can still produce is confirmatory, and the small size of every archive
effect relative to the detectable effect of any available pool (0.29 % at 512 clusters) stands. A frozen exploratory
diagnostic on the N2 replay (archive use 11; `N2_REGIME_DIAGNOSTIC_RESULT_v1.0.md`) found no pre-acquisition
covariate that locates a regime in which departures from G1 help: the O-only kernel edge and G1's own fold-stability
had no association with realised gain given the paired statistic, and the one association found, with design size,
reflects harm concentrated in small records that vanishes in larger ones; no record size in this archive shows a
qualifying benefit, and whether a benefit exists at sample sizes the archive does not contain is unresolved. That
diagnostic also stratified departures by whether the nonlinear pick lay in the set of candidates G1 chose on other
folds of the same block; the inside stratum showed positive realised gain, but the reference set was built from
folds whose training rows include the scored fold, so the stratification is outcome-dependent and is excluded from
any evidence on selector superiority. Finally, a restricted-panel test on an independent public survey (NHANES
2011–2012 serum fatty acids, 1,725 participants with documented identifiers, five prespecified blocks;
`NHANES_FA_RESULT_v1.0.md`, commit fa84662) produced unanimous agreement between G1 and all NOVA variants,
including the unrestricted NOVA-V, in 25 of 25 block-fold decisions, while the selected measurement reduced
base-panel loss by 16 %; the prespecified feasibility gate therefore returned NO-GO and the 2013–2014 evaluation
cycle was never opened. That is a replication of selector agreement where acquisition value is concentrated in one
chemically related candidate, not a test of NOVA-S's advantage; the 25 decisions share participants and analytes
and are not independent experiments.

## 8. Open problems

Simultaneous lower bounds D_{j,g} ≥ D̂_{j,g} − q_{1−α} ŝ_{j,g} − b̂_{j,g} − r_{n,j,g} with an explicit
nuisance-bias term b̂ and selection term r, valid over a stated model class at n ≈ 150 and informative enough that
the rule still departs when it should: via honest nuisance-error envelopes, higher-order bias correction with
independent splits, or finite-sample bounds on fitted-procedure risk differences (which changes the estimand and
must say so). And a confirmatory design whose detectable effect is commensurate with effects of order 0.1%.

A third problem is raised by the replay. Let Z be pre-selection information (here the statistic t) and define the
conditional policy value 𝒱(Z) = E[max{0, max_j E[D_{j,g} | Z]}] for the additive loss-difference objective; it is
the value of the best departure rule that sees Z. Its plug-in, the weighted sum of positive parts of the bin means
in §6.3, is about 0.00076 loss units per cell, roughly 0.12% of the G1 per-cell loss, the same magnitude as the
hindsight envelope. That plug-in is upward biased because the maximum and the positive part are convex, and it is
computed on the archive the selector was developed on; it is a descriptive oracle estimate, not a bound on
attainable policy value. A treatment of 𝒱(Z) with the optimism controlled, by independent evaluation or by
complexity-controlled bounds, is left to a separate note and to data not yet used.

## S1. Provenance (hashes)

| Item | SHA-256 |
|---|---|
| NOVA charter v0.2 | ed7720503cc996d202a01dabffae5736a99f50a310203cdedf25acb61130e049 |
| p6_n1_synthetic.py / n1_summary.json | e756ea4309a99fa173b54a2d76d18ea8eb98225a391cc933ea971877a75258ac / aef926df9c96d54f87f36b5ec89dfecb28b86a09ba026c71bd77a7aa8f4cc178 |
| N1_RESULT_RECORD.md | 254881a73545f07a… (package nova_n1.zip 7a8a91ace152f4823aac12b5a6c7cc9446742359543e95f70e251528ca09bced) |
| p6_n0_regparity.py / n0_summary.json | d9aaefb70b5b243824323b7ef24c56b22921abb4a1ee030ecd9fdc118b6aaef3 / 9e3a78f6e1e356e2892bfb90464062a1fc177fd305e8ca1d9f0d610d1cdc0b9f |
| N2_PROTOCOL_v1.0.md / p6_n2.py / n2_summary.json | 6bbafb3c3e1867d439675ff10211af7dce02602d32af2471724837ac8cd80e57 / 681ce852609e6ab4fcaf0f0ff603477c64e2351e10bd119811d45f008d4bbe8e / 4662e9c9b52c605bb8f79a20f5095acb9669043ed0619add39e37c8ead6b56d8 |
| N2_RESULT_RECORD.md | 569933e124fd43232cb0484167a1ada765956d793cd9c7a9a2eec4596a50ffce |
| p7_x1_synthetic.py (X0/X1 run) / x1_summary.json | a78de1bf5a27a1dfeeaa4defe621068d3fcaa0ecadcd313b48bd42e927ce1651 / 25ae6d3a5f6987e1ed947c1e6955a91d68d6fc154e86b3159b9c29f24b66c603 |
| p7_x1_synthetic.py (X0b/X1b) / x1b_summary.json | 20f8ab35a8e0db66f71c95c5768230502b6788b2edd71fdc4d4fdaf65baadbaa / f05031b507ceda96ef1f8aa53107c497e83af13e383d264e68e86573637dffa5 |
| X01_RESULT_RECORD.md (with addendum) | cd0c3f6bc3b21715e4ce19e0f2b3268eb36792c0b2f6dc0313bb9e4bd471c1fb |
| Gate 3A protocol / summary | ed56e78c70b764f44ce14e6788c78243d8676c9b20bce90eaf090013d61b7024 / 2166d914f7275701e85ad22db487f548495195d39ae69ca4b1bff3a52f4e7204 |
| Gate 2 decision record | a1c1acab0b5b6e04c246ea7ee5de8467a1c4fa68decc69f3990bf1e00f3a4079 |
| Gate D deposit (OSF egnqy) | 16bdb343d898dc28ab40113a6e68a21975430018a559afb60742387884570b7a |
| p8_hierarchy.py / hier_out/hierarchy.json (nested class) | 3e378e20e0e2420f1eeb29a8a16d22dbd4471b10c9577b948b78e343bfe0c175 / 47a695bfe1af4b227e64e1043b3ef407a4ee7a30ff5435ecd59777a33f09b4ff |
| L0_HIERARCHY_RECORD.md | a1413bd864306c36bdf07693d03afe1d422306c99cc6959b15bc09d0ade4c41e |
| p6_n2_replay.py / n2_replay_summary.json / n2_blockfolds.csv | e3a6b7381202984d5040df9fb69be74bf3d1d85b9cf58a0b20bbf1b578d7bbfe / d72f35d86c564d3186a73841445a3300bf71c420931c7c655c6a32be5aabe4ed / b7934523cfc124510fdb72c0469089eecd8b8f4f3676b68552eb28618957d3bf |
| N2_REPLAY_RECORD.md (with addendum) | 2380993915b8ff3f0ff0b406d0342ac4921333943281dbad55ba790e5c74fe4c |
| ESTIMAND_AUDIT_RECORD.md | c54e734e2e16690517a0ecabf9ba198e8a9394f6d2dd61aeab3eab7f98715e18 |

## S2. Coverage replication at 100 blocks

`x1b_100/x1b_summary.json` `a65a68ccfc5d7b81040ab9289b6728552c0618969620db88a85476ca09ef68b7` (same script and seed
offsets as X1b; RBF nuisance; α = 0.05).

| regime | coverage vs Bayes D | coverage vs procedure D | departures | false departures | mean max |D^proc − D^Bayes| |
|---|---|---|---|---|---|
| interaction | 0.70 | 0.90 | 0.10 | 0 | 0.175 |
| quadratic | 0.94 | 0.95 | 0.77 | 0 | 0.106 |

Binomial standard error at 100 blocks ≈ 0.03 (at 0.9) and 0.02 (at 0.95).
