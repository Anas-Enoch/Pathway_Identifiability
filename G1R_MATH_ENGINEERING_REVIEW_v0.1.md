# Acquisition after NOVA-D and G1-R: mathematical engineering review v0.1

Date: 8 November 2026. Status: source-grounded mathematical review and prospective decision memo; not a protocol, freeze, experiment, or adopted challenger.

## 1. Conclusion and limits of this review

The failed joint criterion in G1-R identifies a decision problem that a better covariance estimator alone does not solve. The next potentially defensible contribution is a finite-sample bound on acquisition regret or on harm relative to the actual estimated baseline, together with an explicit demonstration that the bound is informative at the intended sample sizes. A new gate without that demonstration risks repeating NOVA-S/X: safe-looking abstention, weak identification, or a certificate for the wrong estimand.

No new challenger is established by this review. The external selector-improvement slot remains empty. The closed archive cannot validate any construction proposed after reading its results.

This is a bounded review, not a claim to have read every uploaded book. The user identifies the intended sources as files uploaded to a ChatGPT project. Project metadata is visible, but this session has no tool that exposes those project attachments. The previous local inventory lists 23 PDFs and an EPUB; that inventory is evidence of an earlier inventory, not access to those books or proof of their contents. Nine locally accessible PDFs were extracted and indexed. Selected complete pages and relevant passages were examined, especially in Tsybakov, Golub and Van Loan, Drton et al., Shapiro et al., and the compositional-data volume. Two additional primary-source literatures were checked online. The uploaded corpus still needs direct access before an exhaustive source-by-source conclusion is justified.

`source_inventory.json` records local paths, hashes, page counts, and extraction status. `passage_index.json` is a search index, not a reading log. An apparently substantial book file can be a preview: the local *Asymptotic Statistics* has only 47 PDF pages. The extracted *High-Dimensional Statistics* identifies David J. Olive as author; it must not be attributed to another author from its generic filename. The local *Constrained Optimization on Manifolds* is a dissertation by Julián Ortiz, not a generic monograph with that title.

## 2. The objective must precede the method

Let O be already observed variables, T the target, and X_j a candidate. Under a fixed target metric W >= 0,

    V_j = E || E[T | O, X_j] - E[T | O] ||_W^2.

This is the reduction in population Bayes squared prediction risk. For Gaussian variables, write

    K = Cov((T,P) | O), c_j = K_Tj, v_j = K_jj,
    G_j = c_j' W c_j / v_j,
    r(j) = (max_k G_k - G_j) / R_empty.

These identities require fixed units and a defined W. A sample-standardized scoring algorithm may be evaluated against a population-standardized target, but those are different objects and must be declared. For a fitted predictor, its estimation error adds another term; exact selection regret is not its deployed prediction risk.

Three other objectives are distinct: covariance Frobenius error, information gain, and structural pathway identifiability. For example, Gaussian target information gain is

    (1/2) log[det(K_TT) / det(K_TT - c_j c_j'/v_j)],

whereas G_j is trace reduction for W = I. No theorem converts improvement in one into improvement in the other without further conditions. Topological stability, geometric fidelity, and scientific interpretability likewise need a bridge to the stated loss.

Drton, Sturmfels and Sullivant, PDF p.73, gives the Gaussian conditional covariance/Schur-complement characterization. Boyd and Vandenberghe, §§6.4 and 7.5, provide robust approximation and distinguish experiment-design criteria; those tools do not establish acquisition superiority in this programme.

## 3. Failure ledger: what a successor must explain

| Existing evidence | Consequence for a successor | What it does not establish |
|---|---|---|
| GDOA-0 / Gate 3A / N2: flexible prediction did not establish useful acquisition superiority | Separate improvement on O from improvement in choosing X_j; compare arms under the same predictor and tuning budget | Nonlinearity is universally useless |
| N2 replay: existing statistic's threshold envelope remained below the practical floor | Changing a threshold is insufficient motivation without new information or a materially better error bound | All possible data-dependent departure rules obey that envelope |
| NOVA-X: nuisance bias impaired Bayes-value coverage; simultaneous correction lost power | Specify what is certified, address bias as well as variance, and require both coverage and useful departures | A larger bootstrap alone will repair biased nuisance estimates |
| NOVA-D v1.2: both designated recovery tests failed; near-tie screening had no qualifying improvement | Prove that a proposed screening rule retains beneficial opportunities and that the estimator can identify candidates before certification | A universal refutation beyond the tested family/grid |
| Exact-ordering notes: operative shrinkage changes the baseline population ranking | Apply each arm's complete algebra at its actual training-size settings; retain all ridges and scales | Attribution of every finite-sample mistake to shrinkage |
| G1-R v1.1: H1 supported, H2 not established | Improvement where A0 is distorted is insufficient; near-optimal A0 must be protected | H2 refuted or harm beyond the margin established at D2 |
| G1-R sparse contrasts: 41 both sparse, 17 A2-only sparse, 14 A0-only sparse | Report both directions and distinguish low misselection counts from absence of decision-relevant differences | Every sparse contrast has no opportunity or no harm |
| Earlier source-based predictor report, 9 Oct | Kernel features, validation separation, and bias-variance arguments are already part of this programme | Renaming them creates an untested method or novelty |
| DEP-C synthetic certification record | Testing lower bounds and oracle comparisons have already been investigated | Its discrete-bin model supplies a quantitative bound for Gaussian acquisition without a new reduction |

The strongest new engineering requirement is a complete success case and failure case before a new freeze: not merely a population flip or a more stable score, but a bound with useful scale and a counterexample that the proposed rule handles for a stated reason.

## 4. What shrinkage actually changes

For a standardized covariance S with mean diagonal 1, 0 <= alpha < 1, and observed-block ridge epsilon, joint shrinkage gives the exact matched parameters

    lambda = (alpha + epsilon)/(1-alpha),
    f = alpha/(1-alpha),
    C_alpha = (1-alpha) C_Schur(lambda),
    v_alpha = (1-alpha) [v_Schur(lambda) + f],
    score_alpha = (1-alpha) ||C_Schur(lambda)||^2/[v_Schur(lambda)+f].

Clipping must transform with scale; alpha = 1 needs separate handling. The common positive factor does not affect the argmax. The regression residual form is different:

    C_residual = C_Schur - lambda B_T' B_P,
    v_residual,j = v_Schur,j - lambda ||B_P,j||^2.

Ridge changes the conditional cross-covariance and creates unequal ridge-conditioned candidate variances even when true conditional variances are equal. The floor can then interact with that inequality. Equal flip totals are insufficient to show equal selections. The matched Schur arm is an equivalence control, not a novel estimator; changing only its representation cannot improve its decisions.

These are established by the programme's derivation and exact-ordering checks. They motivate independent tuning, but the failed G1-R joint criterion prevents treating independent tuning as an adopted improvement.

## 5. Derivation: ranking stability requires a gap, not merely stable arithmetic

The following elementary bound is derived here to expose the bottleneck; no novelty is claimed.

Partition a population covariance as [[A,B'],[B,D]], with A the observed block, and K = D - B A^-1 B'. Let estimated blocks satisfy spectral errors at most eta, lambda_min(A) >= m > eta, and ||B|| <= M. Resolvent and product bounds give

    ||K_hat-K|| <= e_K
    e_K = eta + (2 M eta + eta^2)/(m-eta)
                + M^2 eta/[m(m-eta)].

Proof: expand (B+E_B)(A+E_A)^-1(B+E_B)' - B A^-1 B'. Bound the two linear cross terms, the quadratic E_B term, and the inverse difference using

    ||(A+E_A)^-1|| <= 1/(m-eta),
    ||(A+E_A)^-1-A^-1|| <= eta/[m(m-eta)].

For W = I, suppose ||c_j|| <= L, v_j >= b > e_K, and both cross-covariance and variance errors are bounded by e_K. Then

    |G_hat,j-G_j| <= e_G
    e_G = (2 L e_K + e_K^2)/(b-e_K)
          + L^2 e_K/[b(b-e_K)].

Consequently a unique population leading gap greater than 2 e_G preserves the leading candidate. This is a sufficient condition, not a necessary one. Unknown m, M, L and b need observable confidence bounds before this becomes an implementable statistical certificate.

Two lessons follow. First, near-singularity amplifies uncertainty in conditional quantities. Second, even numerically perfect linear algebra cannot recover a ranking whose population gap is smaller than statistical uncertainty. Shrinkage may reduce variance while changing the estimand; the bound does not endorse either zero shrinkage or joint shrinkage universally.

Golub and Van Loan, §2.6, PDF pp.111 onward (printed pp.87 onward), explains sensitivity of square systems; §3.2, PDF pp.143-144, gives the Schur complement. The expression above is an application of those tools to this decision problem, not a quoted theorem from the book. A numerical roundoff bound must never be substituted for a confidence bound on eta.

## 6. Derivation: a baseline-relative uncertainty rule and its real limitation

Let g(D) be the actual frozen baseline pick from the design data, not the population oracle. Let U(D) be a simultaneous confidence region for population covariance, with coverage at least 1-delta under an explicitly stated model class. For each candidate define

    L_j(D) = inf_{Sigma in U(D)} [G_j(Sigma)-G_g(D)(Sigma)] / R_empty(Sigma).

Accept only candidates with L_j >= -tau; otherwise return g(D). On the coverage event,

    r(new)-r(g(D)) <= tau.

Proof: substitute the true covariance into the infimum and use r(j)-r(g) = (G_g-G_j)/R_empty. The statement remains valid when g and the new candidate are chosen from the same data because the region covers the whole covariance and the comparison is simultaneous. If only marginal candidate intervals were constructed, this argument would not apply automatically.

For exact normalized Gaussian acquisition regret, 0 <= r(j) <= 1. Thus a worst-case conversion to an expectation bound is

    E[r(new)-r(g)] <= tau + delta.

This does NOT establish the frozen G1-R relative requirement

    E r(new) <= 1.1 E r(g).

To imply that requirement through this particular bound would require tau+delta <= 0.1 E r(g), and E r(g) is unknown. At the scale of D2's reported baseline mean, 0.000960, the permitted excess is only 0.000096. A generic delta = 0.05 certificate is therefore far too coarse through this expectation conversion. This arithmetic uses D2 only to illustrate scale; it must not tune a successor to the closed result.

This is a limitation of the bound, not an impossibility theorem for all safe methods. A tighter expected-tail bound, a justified restricted parameter class, or an absolute harm budget could change the situation. Changing the relative tolerance to an absolute one changes the scientific claim and cannot retrospectively rescue H2.

The construction is prospective and unimplemented. Its coverage, informativeness, numerical optimization and ability to depart remain unestablished. Minimized differences of covariance ratios with Schur complements are not automatically convex; do not claim an SDP or SOCP formulation without proving the reformulation and its direction of validity. Conservative relaxations need a certified bound in the safe direction.

Baseline-protected improvement is existing research, including SPIBB in batch reinforcement learning. Its model, coverage and sampling assumptions are not automatically satisfied by metabolite acquisition. The potentially useful work here is a calibrated, nonvacuous bound for this Gaussian covariance decision, not the generic instruction to retain the baseline under uncertainty.

## 7. Derivation: test information limits before a new simulation grid

For two Gaussian population models with different best candidates and the same normalized two-candidate gap Delta, any decision rule induces a binary test. If P0^n and P1^n are the design-data distributions,

    max_i E_i r(j_hat) >= (Delta/2) [1-TV(P0^n,P1^n)].

The sum of binary test error probabilities is at least 1-TV, and the larger regret is at least half their gap-weighted sum. For iid data Pinsker gives

    TV(P0^n,P1^n) <= sqrt(n KL(P0,P1)/2).

For equal-mean nonsingular Gaussians of dimension d,

    KL(P0,P1) = (1/2)[tr(Sigma1^-1 Sigma0)-d
                           + log(det Sigma1/det Sigma0)].

Use a lower bound truncated at zero. This calculation can determine whether particular beneficial and harmful regimes can be distinguished at the proposed n. It cannot prove that G1 is admissible, rule out every improvement on a particular baseline, or supply exact limits from a poorly chosen pair of alternatives.

Tsybakov, Chapter 2, §§2.2-2.4, PDF pp.90-99 (printed pp.79-88), supplies reduction to testing, total variation, KL product additivity and Pinsker's inequality. Best-arm-identification results reinforce the relevance of gaps and information, but standard bandit models have different observation structures: here a reference row may contain all candidate variables, and shared covariance induces dependencies. Importing their sample-complexity constants directly would be invalid.

The engineering task before a new grid is to instantiate a legally specified pair/family and calculate a useful bound. A generic theorem alone is not the contribution, and DEP-C already contains a related lower-bound programme in a different model.

## 8. What the source families contribute, and what they do not

| Source family | Directly useful transfer | Required bridge / reason not to promote it yet |
|---|---|---|
| Matrix computations | Stable solves, conditioning, Schur perturbation, independent algebra checks | Floating-point correctness is not sampling accuracy |
| Statistical decision theory / nonparametric lower bounds | Separate loss, parameter class and observation experiment; testing reductions | Need an acquisition-specific family and constants |
| Robust convex optimization | Design against declared uncertainty; certify relaxations | The acquisition-ratio problem may be nonconvex; uncertainty sets need statistical coverage |
| Stochastic programming | Sample-optimum optimism and independent validation of chosen solutions | Simulation repetitions do not replace independent biological cohorts |
| Gaussian algebraic statistics | Exact covariance constraints and conditional-independence identities | Polynomial identities provide no lower-regret selector by themselves |
| Compositional data analysis | Define whether signal is relative; choose coordinates, reference and zeros policy | Do not impose log-ratios on absolute concentrations without a scientific model; measuring a new component can change the representation |
| Information geometry | Distinguish parameter geometry, divergence and loss; local sensitivity | A nicer metric is not an acquisition-loss theorem |
| Kernels / deep learning | Candidate function classes and regularization mechanisms | Already represented in prior tests; no new information source established |
| Topology / quantum / Lie theory | Potentially relevant only with a measurable structure and an explicit decision bridge | The uploaded originals are not directly available here; no blanket irrelevance claim, but no demonstrated acquisition bridge in the records read |

Shapiro, Dentcheva and Ruszczyński, PDF p.175 (printed p.163), Proposition 5.6, explains sample-optimum downward bias. PDF pp.214-215 and 218 (printed pp.202-203 and 206), §5.6, separates constructing a solution from validating its objective/optimality gap, including constrained and minimax problems. This directly supports retaining the selection-bias label on the grid maximum.

Pawlowsky-Glahn, Egozcue and Tolosana-Delgado, PDF pp.23-24 (printed pp.8-9), defines compositions through relative information and proportional equivalence. It is a reason to identify the measurement model, not permission to change the target loss after observing errors.

## 9. Candidate triage after this review

1. **Best justified mathematics-first question:** can a simultaneous, model-valid confidence region yield a useful baseline-relative expected-regret bound at the intended n, including near-optimal baselines? Establish its scale before constructing a new selector grid. This is open, not an endorsed challenger.
2. **Supporting analytical question:** derive a local information lower bound for beneficial-vs-harmful regimes. Stop if the demanded discrimination is unattainable at feasible n, or if any valid certificate is forced to return G1 almost everywhere. A no-go bound can be a useful final result.
3. **Possible separate representation question:** are the intended assays absolute, relative, censored, noisy, or transported across platforms? A scientifically supported change of observation model might supply information missing from the current criterion. It would be a different prospective question and needs separate data rights and validation.
4. **Not presently justified:** another nonlinear nuisance grid, another threshold scan, retuning beta on G1-R, a geometry-derived score without a loss bridge, or carrying A2 into external validation as an established improvement.

## 10. Minimum engineering specification before any successor freeze

- Define the observation experiment, target metric/units, budget, and the difference between selection regret and fitted-predictor risk.
- State the protected baseline as an executable data-dependent procedure. Distinguish per-condition guarantees from an average over a declared population of conditions.
- Specify absolute versus relative harm and confidence-event versus expected-risk guarantees. Include rare-error cases and zero/near-zero baseline regret.
- Prove any equivalence, reduction, calibration and optimization-bound direction. Retain every ridge, divisor, standardization and clip.
- Evaluate population behavior at every actual training size and regularization setting; include ridge-induced candidate-variance heterogeneity and a separate form/tuning control where relevant.
- Derive a realistic information or sample-size calculation before selecting a large grid. Specify a minimum nontrivial departure/opportunity-retention requirement; abstention alone cannot qualify as useful improvement.
- Use new structures/seeds and independently specified success and failure families. The existing outputs may explain past failures but may not select successor thresholds or redefine H2.
- Keep the proven fail-closed implementation boundary: current hashes, immutable freeze, atomic bound integrity reports, rejection tests, positive-path tests, independent arithmetic, and exact source provenance.
- Fix inference before execution. Joint requirements must both pass; failure to establish noninferiority is not evidence of inferiority. Sparse/rare events need appropriate bounded-risk or exact treatment rather than a degenerate bootstrap.
- Carry numerical and scientific null results forward. Existing books or established methods are not evidence of novelty; complete prior-art review is separate from mathematical correctness.

No experiments or closed-archive selector computations were performed for this memo. The auxiliary deterministic algebra checks use artificial matrices only, not programme samples or outputs. No repository file or frozen record was edited.

## 11. References and source access

Direct local passages used:

- A. B. Tsybakov, *Introduction to Nonparametric Estimation*, Springer, 2009. Local `/Users/nour/Desktop/Introduction to Nonparametric Estimation.pdf`; PDF pp.88-99, especially 90-99.
- G. H. Golub and C. F. Van Loan, *Matrix Computations*. Local `/Users/nour/Desktop/Golub_Van_Loan_Matrix_Computations.pdf`; §2.6, PDF p.111; §3.2, PDF pp.143-144. Supplied-file pagination used rather than assuming an edition.
- M. Drton, B. Sturmfels and S. Sullivant, *Lectures on Algebraic Statistics*, local dated 27 September 2008; `/Users/nour/Desktop/Lectures on Algebraic Statistics.pdf`, PDF p.73, Proposition 3.1.13.
- A. Shapiro, D. Dentcheva and A. Ruszczyński, *Lectures on Stochastic Programming: Modeling and Theory*, local 2009 file; `/Users/nour/Desktop/Lectures on Stochastic.pdf`, PDF pp.175, 214-215, 218.
- V. Pawlowsky-Glahn, J. J. Egozcue and R. Tolosana-Delgado, *Modeling and Analysis of Compositional Data*, 2015; `/Users/nour/Desktop/Modeling and Analysis of Compositional Data.pdf`, PDF pp.22-24, 36-37.
- S. Amari, *Information Geometry and Its Applications*, local corrected 2020 printing, `/Users/nour/Desktop/Information_Geometry_and_Its_Applications.pdf`, PDF pp.15-16 and 139. Limited scope: divergence/metric distinction and local geometry; no acquisition theorem inferred.

Primary online references inspected:

- Boyd and Vandenberghe, *Convex Optimization*, author-hosted PDF, §§6.4 and 7.5: https://stanford.edu/~boyd/cvxbook/bv_cvxbook.pdf . Robust-approximation passage PDF pp.335-337; source does not establish that our covariance-ratio comparison is convex.
- Laroche, Trichelair and Tachet des Combes, *Safe Policy Improvement with Baseline Bootstrapping*, ICML 2019: https://proceedings.mlr.press/v97/laroche19a.html . Existing baseline-protected improvement in batch RL, not an acquisition guarantee.
- Kaufmann, Cappé and Garivier, *On the Complexity of Best-Arm Identification in Multi-Armed Bandit Models*, JMLR 2016: https://www.jmlr.org/beta/papers/v17/kaufman16a.html . Context for information/gap lower bounds; no direct transfer of bandit sampling assumptions claimed.

Programme records read: `G1R_RESULT_v1.1.md`; `NOVA_D_RESULT_v1.2.md`; `SELECTOR_CLAIM_EVIDENCE_MATRIX_v0.4.md` (historical snapshot, not current exhaustive truth); `X01_RESULT_RECORD.md`; `DEPC_RESULT_RECORD.md`; `GDOA0_SPEC.md`; `paper4_predictor_benchmark/SOURCE_READING_AND_PROTOTYPE_REPORT.md`; `paper4_predictor_benchmark/source_inventory.json`, all under `/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/`.

To finish the requested exhaustive uploaded-source audit, provide the source-project PDFs/EPUB as accessible attachments or an exported local folder. The next version should add a per-source reading ledger, theorem/assumption mapping, unresolved proof obligations, and a source-completeness statement. This version is useful without that access, but must not be represented as the completed exhaustive audit.
