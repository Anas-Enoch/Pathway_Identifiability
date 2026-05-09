# Commutator-Based Prioritization of Metabolite Measurements under Partial Observability

**Anas Enoch, MD — Mohammed VI University of Health Sciences (UM6SS), Casablanca, Morocco**

---

## What this repository is

This repository is a **reproducibility package** for a computational methods paper.
It implements a pathway underdetermination framework for partial metabolomics:
given an incomplete metabolite panel and two biological conditions, the method ranks
unmeasured metabolites by their expected reduction of pathway-level structural
ambiguity if added to the panel.

> **Terminology.** The term *pathway underdetermination* is used throughout to
> avoid confusion with the formal dynamical-systems concept of structural
> identifiability (Raue et al. 2009). The framework quantifies alignment ambiguity
> under partial observability; the relationship to formal identifiability is an open
> question.

**This is not:**

- A generic pathway enrichment toolkit
- A causal inference or flux-balance pipeline
- A system with prospectively validated measurement recommendations

The primary benchmark measures **surrogate fidelity**: how accurately the
commutator approximates the oracle metabolite ranking. Biological validation
of specific recommendations requires prospective experimental follow-up outside
the scope of this study.

---

## Repository layout

```
.
├── results/
│   ├── data/                              # processed datasets + pathway mappings
│   │   ├── processed_metabolite_matrix_ST000356.csv
│   │   ├── ST003390_processed.csv
│   │   ├── ST003506_serum_processed.csv
│   │   ├── ST001865_standard.csv          # hypoxia perturbation dataset
│   │   ├── ST001865_metabolite_mapping.csv
│   │   └── core_pathway_mapping.csv
│   ├── scripts/
│   │   └── active/                        # all benchmark scripts
│   │       ├── compute_pathway_features.py
│   │       ├── compute_fgw_alignment.py
│   │       ├── compute_Uk_real.py
│   │       ├── run_jl_stability_benchmark.py
│   │       ├── run_real_multi_pathway_benchmark.py  # primary benchmark
│   │       ├── run_msea_mummichog_benchmark.py      # metabolomics-native proxies
│   │       ├── run_extended_baselines_benchmark.py  # MI / Bayesian OED / active acq
│   │       ├── run_acquisition_simulation.py        # sequential acquisition sim
│   │       ├── run_ST001865_perturbation_seeded.py  # hypoxia perturbation (seeded)
│   │       └── run_ST001865_perturbation.py         # hypoxia perturbation (unseeded)
│   └── results/                           # all benchmark output CSVs
│       ├── real_multi_pathway_results.csv
│       ├── jl_stability_benchmark.csv
│       ├── msea_mummichog_results.csv
│       ├── extended_baselines_results.csv
│       ├── simulation_delta_uk.csv
│       ├── simulation_pathway_tau.csv
│       ├── simulation_auc.csv
│       └── ST001865_perturbation_results.csv
├── figures/                               # manuscript figures (Fig1–Fig10)
├── ST001849/                              # place raw external cohort files here
├── ST002829/                              # place raw external cohort files here
├── README.md
├── references.bib
└── requirements.txt
```

---

## Setup

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip setuptools wheel
python -m pip install numpy pandas scipy scikit-learn POT matplotlib
```

All scripts resolve paths by walking up from their own location to find
`results/data/` and can be run from any directory depth.

---

## Cohorts

### Primary benchmark (three datasets)

| Dataset | Disease context | Source |
|---------|----------------|--------|
| ST000356 | Breast cancer serum vs control | Metabolomics Workbench PR000284 |
| ST003390 | Type 2 diabetes mellitus | Metabolomics Workbench PR002101 |
| ST003506 | Breast cancer-related lymphedema vs control | Metabolomics Workbench PR002152 |

### External stress-test cohorts (two datasets)

These cohorts are evaluated as **computational stress tests of generalisability**,
not as biological validation.

| Dataset | Context | Platform | n | Source |
|---------|---------|----------|---|--------|
| ST001849 | COVID-19 severity | Untargeted LC/MS | 322 | Sindelar et al. 2021 |
| ST002829 | COVID-19 severity (LEOCC) | Metabolon | 609 | Mathew et al. 2021 |

### Perturbation-oriented validation (one dataset)

| Dataset | Context | n | Source |
|---------|---------|---|--------|
| ST001865 | Hypoxia vs normoxia | 16 | Metabolomics Workbench |

---

## Reproducing the final manuscript benchmarks

Run the following scripts in order. Each is self-contained and writes its output
to `results/results/`. Scripts 1–2 must run before 5 (they build the data
structures the later scripts depend on).

### 1. Primary multi-cohort benchmark

```bash
python results/scripts/active/run_real_multi_pathway_benchmark.py
```

**Output:** `results/results/real_multi_pathway_results.csv`
29,750 masking-trial rows across ST000356, ST003390, ST003506. Contains
oracle delta, regret, top-1/top-3 recovery, and Kendall τ for seven predictors.

### 2. Operator stability benchmark

```bash
python results/scripts/active/run_jl_stability_benchmark.py
```

**Output:** `results/results/jl_stability_benchmark.csv`
Stability metrics (CV-FGW, transport drift, Jaccard top-3) for six preprocessing
operators under three perturbation families (noise, dropout, bootstrap).

### 3. Metabolomics-native proxy baselines

```bash
python results/scripts/active/run_msea_mummichog_benchmark.py
```

**Output:** `results/results/msea_mummichog_results.csv`
Evaluates two metabolomics-native proxies — an MSEA/ORA-adapted enrichment
heuristic and a mummichog-inspired activity-weighted score — under the same
masking protocol as the primary benchmark.
Runtime: approximately 5–10 minutes. Pathways with more than 25 metabolites
are skipped (quadratic FGW cost).

### 4. Information-theoretic and Bayesian baselines

```bash
python results/scripts/active/run_extended_baselines_benchmark.py
```

**Output:** `results/results/extended_baselines_results.csv`
Evaluates three stronger baselines:
- `mutual_info`: KSG mutual information between each candidate metabolite
  and the binary condition label.
- `bayes_oed`: Bayesian OED heuristic — Schur-complement variance reduction,
  approximating expected information gain to the pathway-state posterior.
- `active_acq`: Greedy active feature acquisition — logistic-regression
  cross-entropy reduction when each candidate is added.

### 5. Sequential acquisition simulation

```bash
python results/scripts/active/run_acquisition_simulation.py
```

**Outputs:**
- `results/results/simulation_delta_uk.csv` — ΔU_k per reveal step per pathway
- `results/results/simulation_pathway_tau.csv` — pathway ranking stability
- `results/results/simulation_auc.csv` — downstream classification AUC

Simulates sequential metabolite acquisition on ST000356 at 40% masking.
Compares commutator, variance, and random strategies over 50 trials.
Parameters: `INITIAL_MASK_RATE=0.40`, `MIN_MET=3`, `MAX_NODES=50`, `N_STEPS=5`.

### 6. ST001865 hypoxia perturbation validation

```bash
python results/scripts/active/run_ST001865_perturbation_seeded.py
```

**Output:** `results/results/ST001865_perturbation_results.csv`

This is the **seeded, reproducible version** (`GLOBAL_SEED=42`,
`np.random.seed(42)`). Always use this script for results reported in the
manuscript. The unseeded variant (`run_ST001865_perturbation.py`) produces
slightly different τ values across runs due to FGW internal initialization.

---

## ST001865 perturbation-oriented validation

### Dataset preparation

The raw Metabolomics Workbench file `MSdata_ST001865_2.txt` arrives in
transposed format with condition labels embedded in a metadata row.
A dedicated conversion script (`convert_ST001865_standard.py`) was used to:

- parse the metadata row and extract Hypoxia/Normoxia labels,
- transpose the matrix to sample-row format,
- normalise metabolite identifiers and resolve duplicates,
- export `ST001865_standard.csv` (16 samples × 110 features).

The processed file is included in `results/data/`. The condition override
`CASE_CTRL_OVERRIDE = {"ST001865": ("Hypoxia", "Normoxia")}` is applied
in the benchmark script.

### Perturbation benchmark design

The benchmark applies the same masking-based oracle protocol as the primary
benchmark, with masking rates ρ ∈ {0.40, 0.50, 0.60} and 20 trials per rate
(60 total). Seven competing strategies are evaluated on each retained pathway.

### Pathway coverage and non-trivial pathways

Intersection of the 110 measured metabolites with the core pathway map
retained four pathways. Three (Alanine–Aspartate–Glutamate Metabolism,
Glycerophospholipid Metabolism, Valine–Leucine–Isoleucine Degradation)
each contained three metabolites, producing exactly one hidden metabolite
under any masking rate — a forced-selection regime in which all strategies
are trivially equivalent.

**Only Arginine and Proline Metabolism** (6 metabolites, 2–3 hidden nodes)
constitutes a non-trivial perturbation benchmark.

> **Honest limitation.** Because three of four retained pathways are
> informationally degenerate for method comparison, the perturbation result
> is concentrated in a single pathway. The Citrulline finding (below)
> is a concentrated rather than broad signal.

### Key result

On Arginine and Proline Metabolism under hypoxia (60 trials, seeded run):

| Strategy | Top-1 | Kendall τ | nRegret |
|----------|-------|-----------|---------|
| `gnc_commutator` | **0.833** | **0.709** | **0.061** |
| `random` | 0.617 | 0.233 | 0.142 |
| `diffabundance` | 0.467 | 0.067 | 0.207 |
| `variance` | 0.283 | −0.178 | 0.247 |

The commutator significantly outperforms variance (z=6.07, p<0.001) and
differential abundance (z=4.21, p<0.001).

**Citrulline recovery.** The reveal-defined oracle selected Citrulline —
a structurally central node at the arginine/urea cycle interface — as the
top-priority metabolite in 45% of trials (27/60). When the oracle's choice
was Citrulline, each method recovered that choice as follows:

| Method | Citrulline recovery | Fisher exact vs commutator |
|--------|--------------------|-----------------------------|
| `gnc_commutator` | 92.6% (25/27) | — |
| `random` | 66.7% (18/27) | p=0.039 |
| `variance` | 3.7% (1/27) | p<0.001, OR=325 |
| `diffabundance` | 3.7% (1/27) | p<0.001, OR=325 |
| `surrogate` | 3.7% (1/27) | p<0.001, OR=325 |

Citrulline's marginal abundance shift under hypoxia is moderate; its
commutator score is high because its structural position causes large
perturbations to the pathway alignment geometry when revealed. This
demonstrates the distinction between the two acquisition philosophies:
abundance-based heuristics rank by marginal signal magnitude; the operator
framework ranks by structural effect on pathway observability.

---

## Key results summary

### Primary benchmark (3 cohorts, 29,750 hard-subset instances)

| Predictor | Regret | Top-1 | Top-3 | Kendall τ |
|-----------|--------|-------|-------|-----------|
| `gnc_commutator` | **0.018** | **0.774** | **0.821** | **0.599** |
| `surrogate` | 0.024 | 0.631 | 0.728 | 0.401 |
| `bayes_oed` | 0.165 | 0.568 | 0.697 | 0.100 |
| `mummichog_proxy` | 0.147 | 0.574 | 0.703 | 0.204 |
| `active_acq` | 0.159 | 0.516 | 0.702 | 0.017 |
| `variance` | 0.048 | 0.511 | 0.663 | 0.086 |
| `diffabundance` | 0.051 | 0.489 | 0.641 | 0.074 |
| `msea_proxy` | 0.216 | 0.465 | 0.697 | −0.030 |
| `mutual_info` | 0.257 | 0.372 | 0.696 | −0.517 |
| `degree` | 0.083 | 0.371 | 0.541 | 0.063 |
| `random` | 0.091 | 0.312 | 0.514 | 0.000 |

**Three-tier τ hierarchy:** operator-geometric (0.40–0.60) > connectivity-weighted
differential (0.10–0.20) > label-marginal (−0.52 to 0.09).

### External stress-test cohorts

| Cohort | Comm top-1 | Var top-1 | Comm τ | Var τ |
|--------|-----------|-----------|--------|-------|
| ST001849 (LC/MS, n=322) | 0.382 | 0.292 | 0.243 | 0.116 |
| ST002829 (Metabolon, n=609) | 0.348 | 0.255 | 0.172 | 0.072 |

**Operating-regime boundary.** The commutator advantage is preserved for
pathways with |M_k| ≤ 40 nodes and degrades (δτ ≈ −0.006) for large lipid
pathways (|M_k| > 40). This threshold is reproduced independently across
both external cohorts on different platforms.

---

## Important implementation notes

**ST000356 format.** This dataset is stored in a non-standard transposed
format with metadata rows. Handled explicitly via `fix_st000356()` in all
scripts.

**Condition label overrides.** Required for datasets with non-standard
labels:
```python
CASE_CTRL_OVERRIDE = {
    "ST001849": ("severe", "mild"),
    "ST002829": ("severe", "mild"),
    "ST001865": ("Hypoxia", "Normoxia"),
}
```

**FGW regularization.** All results use ε=0.5. Use fixed-scale cost
normalization (`M /= sqrt(d)`, not `M /= M.max()`) to preserve differences
between preprocessing operators. Do not use `SINK_MAX_ITER` above 300 for
the baseline scripts — 5000 causes multi-hour runtimes with no accuracy gain.

**node_features convention.** All benchmark scripts use
`node_features(X) → (n_nodes, 2)` returning [mean, std] per node.
This ensures `cdist(X_s, X_t)` always receives matching feature dimensions
regardless of case/control sample-count imbalance.

**Seeded runs.** Use `run_ST001865_perturbation_seeded.py` (not the
unseeded variant) for all results reported in the manuscript. The seeded
script sets `np.random.seed(42)` globally before any computation to fix
FGW internal initialization.

**Hard subset.** Always filter to `n_hidden >= 2` for final evaluation.
Single-hidden-metabolite trials are near-trivial and inflate apparent
performance. The primary benchmark table reports hard-subset results only.

---

## Output files reference

| File | Contents |
|------|----------|
| `real_multi_pathway_results.csv` | Primary benchmark, 3 cohorts, all predictors |
| `jl_stability_benchmark.csv` | Operator stability under perturbation regimes |
| `msea_mummichog_results.csv` | MSEA and mummichog proxy results |
| `extended_baselines_results.csv` | MI, Bayesian OED, active acquisition results |
| `simulation_delta_uk.csv` | ΔU_k per step per pathway per strategy |
| `simulation_pathway_tau.csv` | Pathway ranking stability (Kendall τ) |
| `simulation_auc.csv` | Classification AUC per step per strategy |
| `ST001865_perturbation_results.csv` | Hypoxia perturbation validation, seeded |

---

## Figures

| File | Description |
|------|-------------|
| `Fig1_revised_pipeline.png` | End-to-end measurement-prioritization pipeline |
| `Fig2_Structural_Ambiguity.pdf` | Structural ambiguity under partial observability |
| `Fig3_operator_stability_heatmap.png` | Operator stability under perturbation regimes |
| `Fig4_jl_vs_randproj_validation.png` | Fixed-scale normalization validation |
| `Fig5_global_hard_subset_benchmark.png` | Global oracle-recovery performance |
| `Fig6_per_dataset_benchmark.png` | Per-dataset commutator vs variance |
| `Fig7_pathway_heterogeneity.png` | Pathway-level performance heterogeneity |
| `Fig8_external_validation.png` | COVID-19 external stress-test (3-panel) |
| `Fig9_commutator_mechanism.png` | Commutator mechanism schematic |
| `Fig_Citrulline_validation.png` | Hypoxia perturbation Citrulline recovery (3-panel) |

---

## Reproducibility checklist

- [ ] Environment created and dependencies installed
- [ ] Raw cohort CSVs placed in `results/data/`
- [ ] `core_pathway_mapping.csv` in `results/data/`
- [ ] `ST001865_standard.csv` in `results/data/`
- [ ] `CASE_CTRL_OVERRIDE` set correctly per dataset
- [ ] `SINK_MAX_ITER = 300` in all baseline scripts
- [ ] `node_features` returns `(n_nodes, 2)` in all scripts
- [ ] `np.random.seed(42)` set in ST001865 seeded script
- [ ] Hard subset (`n_hidden >= 2`) used for final evaluation

---

## GitHub reproducibility note

This repository provides scripts and processed benchmark outputs for
reproducing the reported computational results. It is a reproducibility
package for a computational methods paper, not a polished standalone
software application. Scripts may require minor path adjustments if the
directory layout differs from the structure above.

---

## Citation

If you use this repository, please cite:

> Anas Enoch. *Commutator-Based Prioritization of Metabolite Measurements
> under Partial Observability.* Manuscript in preparation, 2025.

External cohort citations:

> Sindelar et al. *Longitudinal metabolomics of human plasma reveals
> prognostic markers of COVID-19 disease severity.* Cell Reports Medicine,
> 2(8):100369, 2021.

> Mathew et al. *Nucleotide, phospholipid, and kynurenine metabolites are
> robustly associated with COVID-19 severity.* Metabolomics Workbench
> ST002829, 2021–2022.

---

## Contact

**Anas Enoch, MD**
Mohammed VI University of Health Sciences (UM6SS), Casablanca, Morocco

This repository is under active development. If scripts and CSV outputs
disagree, the scripts in `results/scripts/active/` are the source of truth.
