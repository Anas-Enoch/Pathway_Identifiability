# Pathway Identifiability under Partial Metabolomics

**Regime-Aware Operator Selection and Commutator-Based Measurement Prioritization**

*Anas Enoch, MD — Mohammed VI University of Health Sciences (UM6SS), Casablanca, Morocco*

---

## What this repository is

This repository implements a **dataset-conditioned pathway underdetermination framework** for partial metabolomics.

The central question is:

> Under incomplete metabolite coverage, which additional measurement is expected to reduce pathway-level structural ambiguity most effectively?

The framework is built around four linked components:

1. **Pathway-level feature construction** — condition-specific node feature matrices from metabolomics cohorts
2. **Cross-condition alignment** — entropic fused Gromov–Wasserstein (FGW) transport between pathway states
3. **Regime-aware operator selection** — stability benchmarking under transition-constrained perturbations
4. **Metabolite prioritization** — commutator-based surrogate ranked against an exact reveal oracle

> **Note on terminology.** The term *pathway underdetermination* is used throughout rather than *pathway identifiability* to avoid confusion with the formal dynamical-systems concept (structural/practical identifiability in the sense of Raue et al. 2009). The framework quantifies operator-level alignment ambiguity under partial observability; the relationship to formal identifiability is an open question.

---

## What this repository is not

- A generic pathway enrichment toolkit
- A causal inference or flux-balance pipeline
- A claim that JL projection always improves alignment
- A system with biologically validated measurement recommendations

The primary benchmark measures **surrogate fidelity**: how efficiently the commutator approximates the expensive oracle re-solve. Biological validation of specific recommendations requires prospective experimental follow-up, which is outside the scope of this computational study.

---

## Repository layout

```
.
├── results/
│   ├── data/                   # processed metabolomics tables + pathway mappings
│   ├── scripts/
│   │   ├── active/             # main pipeline scripts
│   │   └── *.py                # additional benchmark scripts
│   └── results/                # benchmark output CSVs
├── figures/                    # manuscript figures (Fig1–Fig9)
├── codes/                      # figure-generation assets
├── ST001849/                   # external validation cohort txt files (place here)
├── ST002829/                   # external validation cohort txt files (place here)
├── README.md
├── references.bib
└── requirements.txt
```

---

## Cohorts

### Primary benchmark cohorts (three datasets)

| Dataset | Disease context | Source |
|---------|----------------|--------|
| ST000356 | Breast cancer serum vs control | Metabolomics Workbench PR000284 |
| ST003390 | Type 2 diabetes mellitus | Metabolomics Workbench PR002101 |
| ST003506 | Breast cancer-related lymphedema vs control | Metabolomics Workbench PR002152 |

### External validation cohorts (two datasets)

These cohorts are evaluated as **computational stress tests of generalisability**, not as biological validation. Both use PCA-derived severity labels and correlation-derived pathway structure without curated edges.

| Dataset | Disease context | Platform | n | Pathways | Source |
|---------|----------------|----------|---|----------|--------|
| ST001849 | COVID-19 severity (mild vs severe) | Untargeted LC/MS | 322 | 19 | Sindelar et al. 2021, PR001166 |
| ST002829 | COVID-19 severity (LEOCC cohort) | Metabolon 1108-metabolite | 609 | 21 | Mathew et al. 2021, PR001818 |

---

## Installation

```bash
# Create environment
python3 -m venv .env
source .env/bin/activate

# Install dependencies
pip install -r requirements.txt
# or explicitly:
pip install numpy pandas scipy scikit-learn POT tqdm matplotlib reportlab
```

---

## Pipeline execution

### Step 1 — Build pathway-level features

```bash
python3 results/scripts/active/compute_pathway_features.py
```

Loads each cohort, parses condition labels, handles the ST000356 transposed format, restricts metabolites to each pathway, and writes condition-specific feature and structure files.

### Step 2 — Compute pathway-level U_k scores

```bash
python3 results/scripts/active/compute_Uk_real.py
```

Computes the underdetermination functional from the FGW transport plan (transport entropy + spectral-gap structural term).

### Step 3 — Compute FGW alignments

```bash
python3 results/scripts/active/compute_fgw_alignment.py
```

Aligns case vs control states for each pathway under each preprocessing operator.

### Step 4 — Run the stability benchmark

```bash
python3 results/scripts/active/run_jl_stability_benchmark.py
```

Evaluates operator stability under three perturbation families: `noise`, `dropout`, `bootstrap`.

Operators compared: `none`, `l2`, `jl`, `randproj`, `pca_fixed`, `pca_var95`

Inspect results:

```python
import pandas as pd
df = pd.read_csv("results/results/jl_stability_benchmark.csv")
print(
    df.groupby("method")[["cv_fgw","transport_drift","delta_U","top3_jaccard","rank_tau"]]
      .mean().round(4).sort_values("cv_fgw").to_string()
)
```

> **Operator selection rule:** choose the method with lowest `cv_fgw` and highest `rank_tau`. In the low-node pathway regime tested here, non-projection operators (`none`, `l2`) are stable; `jl` and `randproj` are unstable.

### Step 5 — Run the primary multi-pathway benchmark

```bash
python3 results/scripts/active/run_real_multi_pathway_benchmark.py
```

Runs the metabolite-prioritization benchmark across ST000356, ST003390, ST003506.

Inspect hard-subset results (n_hidden ≥ 2):

```python
import pandas as pd
df = pd.read_csv("results/results/real_multi_pathway_results.csv")
hard = df[df["n_hidden"] >= 2]
print(
    hard.groupby("predictor_method")[["regret","nregret","top1","top3","rank_tau"]]
        .mean().round(4).sort_values(["regret","top1"], ascending=[True,False]).to_string()
)
```

### Step 6 — Run external validation benchmarks

**ST001849:**

```bash
# Place ST001849_AN*.txt files in ST001849/
python3 results/scripts/parse_ST001849.py          # generates benchmark_ready CSV
python3 results/scripts/active/run_ST001849_benchmark.py
```

**ST002829:**

```bash
# Place ST002829_AN004619–AN004622.txt files in ST002829/
python3 results/scripts/parse_ST002829.py          # generates benchmark_ready CSV
python3 results/scripts/build_ST002829_pathway_mapping.py
python3 results/scripts/active/run_ST002829_benchmark.py
```

Results land in `results/results/ST001849_benchmark_results.csv` and `results/results/ST002829_benchmark_results.csv`.

---

## Benchmark design

### What the benchmark measures

The primary benchmark evaluates **surrogate fidelity**: how accurately each candidate prioritization rule recovers the oracle ranking

```
m* = argmax_m ΔU_k(m)
```

where `ΔU_k(m)` is computed by the same underdetermination functional. The oracle and the surrogate are both derived from `U_k`. High fidelity means the commutator efficiently approximates the expensive oracle re-solve — it does not mean the recommendations are biologically validated.

### Predictors compared

| Predictor | Description |
|-----------|-------------|
| `gnc_commutator` | Masked-state operator-commutator score (**best overall**) |
| `surrogate` | Gradient-based underdetermination sensitivity |
| `variance` | Sample-level metabolite variance |
| `diffabundance` | Differential abundance between conditions |
| `degree` | Graph degree centrality in masked pathway |
| `random` | Uniform random selection (null baseline) |
| `mb2d_transport` | Transport-inspired exploratory predictor |

> **Note on baselines:** No existing method is purpose-built for the measurement-prioritization problem as formulated here. The baselines are adapted proxies (variance, differential abundance, centrality), and comparisons measure improvement over reasonable heuristics rather than over an optimal purpose-built competitor.

### Evaluation metrics

| Metric | Meaning |
|--------|---------|
| `regret` | U_k gap between oracle and predictor choice |
| `nregret` | Normalized regret |
| `top1` | Predictor matches oracle top-1 |
| `top3` | Predictor overlaps oracle top-3 |
| `rank_tau` | Kendall τ between predictor and oracle ranking |

> **Recommended subset:** `n_hidden >= 2`. Single-hidden-metabolite cases are near-trivial and inflate apparent performance.

---

## Key results

### Primary benchmark (3 cohorts, 29,750 hard-subset rows)

| Cohort | Comm top-1 | Var top-1 | Comm τ | Var τ |
|--------|-----------|-----------|--------|-------|
| ST000356 (breast cancer) | 0.969 | 0.831 | 0.815 | 0.086 |
| ST003390 (type 2 diabetes) | 0.954 | 0.898 | 0.710 | 0.365 |
| ST003506 (lymphedema) | 0.828 | 0.732 | 0.411 | 0.024 |

### External validation (2 COVID cohorts)

| Cohort | Comm top-1 | Var top-1 | Comm τ | Var τ |
|--------|-----------|-----------|--------|-------|
| ST001849 (LC/MS, n=322) | 0.382 | 0.292 | 0.243 | 0.116 |
| ST002829 (Metabolon, n=609) | 0.348 | 0.255 | 0.172 | 0.072 |

**Operating regime boundary:** commutator advantage is preserved in pathways with |M_k| ≤ 40 nodes and degrades to δτ ≈ −0.006 in large lipid pathways (|M_k| > 40). This threshold is reproduced independently in both external cohorts across different platforms, and is predicted to be resolvable with curated edge topology from LIPID MAPS or Reactome.

---

## Output files

| File | Contents |
|------|----------|
| `results/results/jl_stability_benchmark.csv` | Operator stability metrics by method and perturbation |
| `results/results/real_multi_pathway_results.csv` | Trial-level oracle benchmark (3 primary cohorts) |
| `results/results/ST001849_benchmark_results.csv` | External validation benchmark — ST001849 |
| `results/results/ST002829_benchmark_results.csv` | External validation benchmark — ST002829 |

---

## Figures

| File | Description |
|------|-------------|
| `Fig1_revised_pipeline.png` | Four-layer pipeline overview |
| `Fig2_Structural_Ambiguity.pdf` | Structural ambiguity under partial observability |
| `Fig3_operator_stability_heatmap.png` | Operator stability under perturbation regimes |
| `Fig4_jl_vs_randproj_validation.png` | Fixed-scale normalization validation |
| `Fig5_global_hard_subset_benchmark.png` | Global oracle-recovery performance |
| `Fig6_per_dataset_benchmark.png` | Per-dataset commutator vs variance |
| `Fig7_pathway_heterogeneity.png` | Pathway-level heterogeneity of performance |
| `Fig8_external_validation.pdf` | COVID-19 external validation (3-panel) |
| `Fig9_commutator_mechanism.png` | Commutator mechanism schematic |

---

## Important implementation notes

**ST000356 parsing.** This dataset is stored in a non-standard transposed format with metadata rows. Handled explicitly in `compute_pathway_features.py` and `run_real_multi_pathway_benchmark.py`.

**Condition labels.** ST001849 and ST002829 use `severe`/`mild` labels. An explicit override is needed:
```python
CASE_CTRL_OVERRIDE = {
    "ST001849": ("severe", "mild"),
    "ST002829": ("severe", "mild"),
}
```

**FGW regularization.** All results use `ε = 0.5` (Sinkhorn entropic regularization). Results are stable for ε ∈ [0.25, 1.0]; `ε = 0.1` causes convergence failures in ~18% of pathway solves. Do not normalize the cost matrix by its maximum — use fixed-scale normalization (`M /= sqrt(d)`) to preserve differences between preprocessing methods.

**Structure matrices.** When explicit pathway edges are unavailable, structure falls back to a correlation-distance matrix. This reduces commutator discriminability in large lipid pathways.

**Small pathways.** Pathways with |M_k| ≤ 3 are excluded from benchmarking. Pathways with 4–6 nodes should be interpreted cautiously for projection-based and rank-based metrics.

**Path resolution.** Benchmark scripts use walk-up path resolution (searching up to 6 parent directories for `results/data/`) and run correctly from any nesting depth including `results/scripts/active/`.

---

## Reproducibility checklist

- [ ] Environment created and dependencies installed
- [ ] Raw cohort CSVs placed in `results/data/`
- [ ] `core_pathway_mapping.csv` in `results/data/`
- [ ] `ST001849_pathway_mapping.csv` in `results/data/` (for external validation)
- [ ] `ST002829_pathway_mapping.csv` in `results/data/` (for external validation)
- [ ] `CASE_CTRL_OVERRIDE` set correctly for each dataset
- [ ] `ε = 0.5` in all FGW calls
- [ ] Operator selected by lowest `cv_fgw` from stability benchmark
- [ ] Hard subset (`n_hidden >= 2`) used for final evaluation

---

## Citation

If you use this repository, please cite:

> Anas Enoch. *Pathway Identifiability under Partial Metabolomics: Regime-Aware Operator Selection and Commutator-Based Measurement Prioritization.* Manuscript in preparation, 2025.

External cohort citations:

> Sindelar et al. *Longitudinal metabolomics of human plasma reveals prognostic markers of COVID-19 disease severity.* Cell Reports Medicine, 2(8):100369, 2021.

> Mathew et al. *Nucleotide, phospholipid, and kynurenine metabolites are robustly associated with COVID-19 severity.* Metabolomics Workbench ST002829, 2021–2022.

---

## Contact

**Anas Enoch, MD**
Mohammed VI University of Health Sciences (UM6SS), Casablanca, Morocco
`anas_nour@um5.ac.ma`

This repository is under active development. If scripts, figures, and CSVs disagree, the CSV outputs in `results/results/` and the scripts in `results/scripts/active/` are the source of truth.
