## Figure 1. Structural ambiguity in metabolic pathways under partial metabolomics
Incomplete metabolomics panels induce structural ambiguity at the pathway level. Observed metabolites (solid blue), panel-missing metabolites (hollow grey), and limit-of-detection (LOD)-censored metabolites (dashed red) are explicitly represented within a single pathway graph. Multiple, functionally distinct pathway interpretations—such as divergent flux routing at branch points—remain equally compatible with the same observed data. This ambiguity motivates an identifiability-centric approach that prioritizes resolving uncertainty rather than enumerating possible completions.

## Figure 3. Distributional node features encode uncertainty from missingness and censoring
Each metabolite node is represented by a distributional state parameterized by a mean and variance, reflecting assay noise, panel-missing uncertainty, or LOD censoring. Observed metabolites exhibit tightly concentrated distributions, while unmeasured or censored metabolites carry broader or truncated uncertainty. This representation enables uncertainty-aware pathway comparison and prevents missing data from being collapsed into point estimates or heuristic imputations.

## Figure 4. Johnson–Lindenstrauss–stabilized fused Gromov–Wasserstein alignment
Pathway comparison across conditions is performed using fused Gromov–Wasserstein (FGW) alignment over topology and node features. Direct alignment in high-dimensional feature spaces is unstable under sparse and heterogeneous measurements. Random projection via Johnson–Lindenstrauss embeddings stabilizes the geometry of node features, yielding consistent transport plans and reproducible pathway correspondences across runs and conditions.

## Figure 6. Measurement impact estimation without enumerating pathway completions
For underdetermined pathways, candidate next measurements are ranked by their estimated impact on reducing pathway ambiguity. Measurement impact is approximated as the product of local sensitivity of the underdetermination functional and the uncertainty associated with the metabolite. This estimator enables principled experimental prioritization without enumerating all possible metabolite completions or pathway instantiations.

## Figure 8. Regret-based evaluation of measurement recommendation quality
Measurement recommendations are evaluated under a synthetic masking protocol using regret as the primary metric. Regret quantifies the gap between the ambiguity reduction achieved by the recommended measurement and the best possible measurement revealed in hindsight. The full framework consistently achieves lower regret than random, centrality-based, or uncertainty-only baselines, demonstrating effective prioritization of informative measurements.

## Figure 9. Generalization across pathway topologies and missingness regimes
The framework generalizes across diverse pathway topologies, sizes, and missingness patterns. For each pathway, observed and latent metabolites are explicitly represented, and a calibrated underdetermination score is computed. Consistent behavior across heterogeneous graphs demonstrates that the proposed identifiability and measurement prioritization principles are not pathway-specific but apply broadly across metabolic network structures.
