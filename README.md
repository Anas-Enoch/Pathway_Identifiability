# Commutator-Based Prioritization of Metabolite Measurements under Partial Observability

**Anas Enoch, MD — Mohammed VI University of Health Sciences (UM6SS), Casablanca, Morocco**

---

## What this repository is

This repository is a **reproducibility package** for a computational methods paper addressing
a practical problem in metabolomics:

> When a metabolomics panel covers only a subset of pathway metabolites,
> which additional measurement would most reduce pathway-level structural ambiguity?

The framework ranks unmeasured metabolites by their expected contribution to restoring
pathway observability — evaluated against a reveal-defined oracle across five public cohorts.

**This repository is not** a production-ready metabolomics software platform, an automated
clinical decision system, or a tool for causal biomarker discovery. It provides scripts and
processed benchmark outputs sufficient to reproduce the reported computational results.

> **Terminology.** *Pathway underdetermination* is used throughout to avoid confusion with
> the formal dynamical-systems concept of structural identifiability (Raue et al. 2009).
> The framework quantifies alignment ambiguity under partial observability.

---

## Benchmark architecture

```
PRIMARY LAYER — Multi-cohort identifiability benchmark (core reproducibility)
─────────────────────────────────────────────────────────────────────────────
  ST000356   Breast cancer serum vs control          }
  ST003390   Type 2 diabetes mellitus                }  Principal benchmark
  ST003506   Breast cancer-related lymphedema        }

  ST001849   COVID-19 severity — untargeted LC/MS    }  External generalization
  ST002829   COVID-19 severity — Metabolon           }  stress tests

SECONDARY LAYER — Extended baseline comparison
──────────────────────────────────────────────
  mummichog proxy, MSEA proxy,
  mutual information, Bayesian OED, active feature acquisition,
  variance, differential abundance, degree, random

TERTIARY LAYER — Perturbation-oriented validation
──────────────────────────────────────────────────
  ST001865   Hypoxia vs normoxia (mechanistically anchored stress test)
             — additional biological extension, not the primary benchmark
```

The three-cohort primary benchmark is the core of the paper.
The two external cohorts test cross-platform generalisability.
ST001865 is a separate perturbation-oriented module.

---

## Repository layout

```
.
├── results/
│   ├── data/                              # processed datasets + pathway mappings
│   │   ├── processed_metabolite_matrix_ST000356.csv
│   │   ├── processed_metabolite_matrix_ST003390.csv
│   │   ├── ST003506_serum_processed.csv
│   │   ├── ST001849_benchmark_ready.csv
│   │   ├── ST002829_benchmark_ready.csv
│   │   ├── ST001865_standard.csv
│   │   └── core_pathway_mapping.csv
│   ├── scripts/
│   │   ├── active/                        # all runnable benchmark scripts
│   │   │   ├── run_real_multi_pathway_benchmark.py  ← PRIMARY
│   │   │   ├── run_jl_stability_benchmark.py
│   │   │   ├── run_ST001849_benchmark.py
│   │   │   ├── run_ST002829_benchmark.py
│   │   │   ├── run_msea_mummichog_benchmark.py
│   │   │   ├── run_extended_baselines_benchmark.py
│   │   │   ├── run_acquisition_simulation.py
│   │   │   └── run_ST001865_perturbation_seeded.py
│   │   └── prep/                          # data preparation scripts
│   │       ├── parse_ST002829.py
│   │       ├── parse_st003390.py
│   │       └── build_ST002829_pathway_mapping.py
│   └── results/                           # all benchmark output CSVs
├── figures/                               # manuscript figures (Fig1–Fig10)
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

Run all commands from the **repository root** unless otherwise noted.

---

## Minimal quickstart

Run the principal multi-cohort benchmark:

```bash
python results/scripts/active/run_real_multi_pathway_benchmark.py
```

This reproduces the primary manuscript benchmark across three cohorts
(ST000356, ST003390, ST003506) using the commutator surrogate and all
competing baselines.

**Output:** `results/results/real_multi_pathway_results.csv`

Inspect results:

```python
import pandas as pd
df   = pd.read_csv("results/results/real_multi_pathway_results.csv")
hard = df[df["n_hidden"] >= 2]   # hard subset: ≥2 hidden metabolites per trial
print(
    hard.groupby("predictor_method")
        [["regret", "nregret", "top1", "top3", "rank_tau"]]
        .mean().round(4)
        .sort_values("rank_tau", ascending=False)
        .to_string()
)
```

For external generalization results, see [External Cohort Benchmark](#external-cohort-benchmark) below.

---

## Benchmark cohort structure

### Primary benchmark cohorts

These three cohorts constitute the core benchmark and produce the 29,750
hard-subset evaluation instances reported in Table 1.

| Dataset | Disease context | n | Source |
|---------|----------------|---|--------|
| ST000356 | Breast cancer serum vs control | 134 | Metabolomics Workbench PR000284 |
| ST003390 | Type 2 diabetes mellitus | 300 | Metabolomics Workbench PR002101 |
| ST003506 | Breast cancer-related lymphedema vs control | 34 | Metabolomics Workbench PR002152 |

### External generalization cohorts

These two cohorts are evaluated as **computational stress tests of cross-platform
generalisability**. They are not included in the primary benchmark count.

| Dataset | Context | Platform | n | Source |
|---------|---------|----------|---|--------|
| ST001849 | COVID-19 severity | Untargeted LC/MS | 322 | Sindelar et al. 2021, PR001166 |
| ST002829 | COVID-19 severity (LEOCC) | Metabolon | 609 | Mathew et al. 2021, PR001818 |

### Perturbation validation dataset

| Dataset | Context | n | Source |
|---------|---------|---|--------|
| ST001865 | Hypoxia vs normoxia | 16 | Metabolomics Workbench ST001865 |

---

## PRIMARY multi-cohort benchmark

### The oracle and what the benchmark measures

The benchmark evaluates **surrogate fidelity**: how accurately each candidate
method approximates the reveal-defined oracle ranking.

For each pathway k and each hidden metabolite m, the oracle acquisition score is:

```
ΔU_k(m) = U_k(masked panel) − U_k(panel with m revealed)
```

A higher ΔU_k(m) means revealing m produces a larger reduction in pathway-level
structural ambiguity. The oracle-optimal metabolite is:

```
m* = argmax_m ΔU_k(m)
```

The benchmark evaluates how well each surrogate method recovers this oracle
ranking without recomputing the full FGW re-solve for every candidate.

### Run the primary benchmark

```bash
# Reproduces primary multi-cohort benchmark — commutator vs all baselines
python results/scripts/active/run_real_multi_pathway_benchmark.py
```

**Output:** `results/results/real_multi_pathway_results.csv`

### Operator stability benchmark

Before the primary benchmark, select the stable preprocessing operator:

```bash
python results/scripts/active/run_jl_stability_benchmark.py
```

**Output:** `results/results/jl_stability_benchmark.csv`

Evaluates six preprocessing operators (none, l2, jl, randproj, pca\_fixed,
pca\_var95) under three perturbation families. In the tested pathway regime,
non-projection operators (`none`, `l2`) are stable; `jl` and `randproj` are not.
`METHOD_DEFAULT="none"` is used throughout.

### Expected outputs

The primary benchmark produces the following per-trial metrics:

| Metric | Description |
|--------|-------------|
| `regret` | U_k gap between oracle and predictor choice |
| `nregret` | Normalised regret |
| `top1` | Predictor matches oracle top-1 metabolite |
| `top3` | Predictor overlaps oracle top-3 metabolites |
| `rank_tau` | Kendall τ between predictor and oracle full ranking |

Summary outputs are in `results/results/`. Combined manuscript benchmark summary (primary cohorts,
hard subset n_hidden ≥ 2, all predictors):

| Predictor | Regret | Top-1 | Top-3 | τ |
|-----------|--------|-------|-------|---|
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

Three-tier τ hierarchy: operator-geometric (0.40–0.60) > connectivity-weighted
differential (0.10–0.20) > label-marginal (−0.52 to 0.09).

---

## External cohort benchmark

### Preprocessing (run once)

The external cohorts require conversion into benchmark-ready sample × metabolite
matrices before running benchmark scripts. The processed CSVs are already included
in `results/data/`; re-run preprocessing only if you modify the raw data.

```bash
# ST001849/ and ST002829/ are local-only raw download folders.
# They are not committed to the repo. Skip this block if the processed
# benchmark_ready CSVs are already present in results/data/.
python results/scripts/prep/parse_ST002829.py
python results/scripts/prep/build_ST002829_pathway_mapping.py
python results/scripts/prep/parse_st003390.py   # if re-processing ST003390
```

Additional cohort-specific preprocessing scripts may be required depending
on raw dataset format. ST001849 benchmark-ready CSVs are included in
`results/data/` and do not require separate preprocessing.

Preprocessing steps: normalise metabolite naming, harmonise condition labels,
generate benchmark-ready CSV matrices, and construct pathway-overlap-compatible
metabolite panels. The processed matrices are then consumed by the benchmark scripts below.

### Run external benchmarks

```bash
python results/scripts/active/run_ST001849_benchmark.py
python results/scripts/active/run_ST002829_benchmark.py
```

**Outputs:** `ST001849_benchmark_results.csv`, `ST002829_benchmark_results.csv`

| Cohort | Comm τ | Var τ | Platform |
|--------|--------|-------|----------|
| ST001849 (n=322) | 0.243 | 0.116 | LC/MS |
| ST002829 (n=609) | 0.172 | 0.072 | Metabolon |

**Operating-regime boundary:** the advantage of structural prioritization is
preserved for pathways with |M_k| ≤ 40 nodes and degrades (δτ ≈ −0.006) for
large lipid pathways. This boundary is reproduced independently across both
cohorts on different platforms.

---

## SECONDARY extended baseline comparison

### Metabolomics-native proxy baselines

```bash
# Reproduces metabolomics-native enrichment proxy baselines (Table 1)
python results/scripts/active/run_msea_mummichog_benchmark.py
```

**Output:** `results/results/msea_mummichog_results.csv`

Evaluates two proxies adapted from pathway enrichment methods:

| Predictor | Top-1 | τ |
|-----------|-------|---|
| `mummichog_proxy` | 0.574 | 0.204 |
| `msea_proxy` | 0.465 | −0.030 |

### Information-theoretic and Bayesian baselines

```bash
# Reproduces MI / Bayesian OED / active-acquisition baselines (Table 1)
python results/scripts/active/run_extended_baselines_benchmark.py
```

**Output:** `results/results/extended_baselines_results.csv`

| Predictor | Top-1 | τ | Design rationale |
|-----------|-------|---|-----------------|
| `bayes_oed` | 0.568 | 0.100 | Schur-complement variance reduction |
| `active_acq` | 0.516 | 0.017 | Greedy logistic-regression loss reduction |
| `mutual_info` | 0.372 | −0.517 | KSG mutual information vs condition label |

`mutual_info` τ = −0.517 is the strongest negative result: information against
the condition label actively inverts the oracle ranking. This occurs because the
oracle is defined by pathway-alignment geometry, not classification signal.

### Sequential acquisition simulation

```bash
python results/scripts/active/run_acquisition_simulation.py
```

**Outputs:** `simulation_delta_uk.csv`, `simulation_pathway_tau.csv`, `simulation_auc.csv`

Simulates sequential metabolite acquisition on ST000356 at 40% masking (50 trials).
At each reveal step, oracle top-1 recovery measures whether each strategy selects
the metabolite that most restores pathway structural observability.

| Strategy | Oracle top-1 (step 1) |
|----------|----------------------|
| commutator | 71.5% |
| variance | 67.5% |
| random | 51.0% |

Commutator vs random: z = 4.30, p < 0.001. Commutator vs variance: directionally
consistent but not significant at n = 4 pathways (p = 0.38). Pathway ranking
stability at 40% masking: τ = 0.707.

---

## TERTIARY perturbation-oriented validation (ST001865)

**This is an additional perturbation-oriented validation module, not the primary benchmark.**

ST001865 provides a mechanistically anchored stress test: it asks whether the
framework recovers structurally informative metabolites under a controlled biological
perturbation, as opposed to a cross-sectional disease association.

```bash
python results/scripts/active/run_ST001865_perturbation_seeded.py
```

**Output:** `results/results/ST001865_perturbation_results.csv`

Always use the seeded variant (`GLOBAL_SEED=42`, `np.random.seed(42)`) for
reproducible results. The unseeded script produces different τ values across runs
due to FGW internal initialization.

### Dataset

ST001865 contrasts hypoxia (n=8) with normoxia (n=8) in a matched design.
Hypoxia induces coordinated rewiring of arginine and nitrogen metabolism —
a biologically interpretable perturbation context.

### Coverage and scope

Four pathways met the coverage threshold after intersection with the pathway map.
Three (three metabolites each) had exactly one hidden metabolite under any masking
rate — forced selection, all strategies trivially equivalent. Only **Arginine and
Proline Metabolism** (six metabolites, two to three hidden) produces a non-trivial
benchmark.

> **Scope:** The perturbation result is concentrated in a single pathway.
> This is a limitation of metabolite coverage, not of the method.

### Results (Arginine and Proline Metabolism, 60 trials, seeded)

| Strategy | Top-1 | τ | vs commutator |
|----------|-------|---|---------------|
| `gnc_commutator` | **0.833** | **0.709** | — |
| `random` | 0.617 | 0.233 | p = 0.008 |
| `diffabundance` | 0.467 | 0.067 | p < 0.001 |
| `variance` | 0.283 | −0.178 | p < 0.001 |

Performance is stable at 40–50% masking (top-1 = 0.95, τ ≈ 0.88) and degrades
at 60% masking (top-1 = 0.60), consistent with the operating-regime boundary
identified in the primary benchmark.

### Citrulline recovery under perturbation masking

The oracle selected Citrulline — a structurally central metabolite at the
arginine/urea cycle interface — as the top-priority acquisition target in 45%
of trials (27/60). When the oracle's top choice was Citrulline, each method
recovered that choice as follows:

| Method | Recovery | Fisher exact vs commutator |
|--------|----------|---------------------------|
| `gnc_commutator` | **92.6%** (25/27) | — |
| `random` | 66.7% (18/27) | p = 0.039 |
| `variance` | 3.7% (1/27) | p < 0.001, OR = 325 |
| `diffabundance` | 3.7% (1/27) | p < 0.001, OR = 325 |

The 25× gap between structural prioritization (92.6%) and abundance-based
methods (3.7%) illustrates the central distinction: abundance-based methods rank
by marginal signal magnitude; the operator-based framework ranks by structural
effect on pathway observability. Citrulline's marginal abundance shift under
hypoxia is moderate, but its structural position causes large perturbations to
the pathway alignment geometry when revealed.

---

## Output files reference

| File | Layer | Contents |
|------|-------|----------|
| `real_multi_pathway_results.csv` | Primary | Main benchmark, 3 cohorts, all predictors |
| `jl_stability_benchmark.csv` | Primary | Operator stability |
| `ST001849_benchmark_results.csv` | Primary | External generalization — ST001849 |
| `ST002829_benchmark_results.csv` | Primary | External generalization — ST002829 |
| `msea_mummichog_results.csv` | Secondary | Enrichment proxy baselines |
| `extended_baselines_results.csv` | Secondary | MI / OED / active acquisition |
| `simulation_delta_uk.csv` | Secondary | ΔU_k per reveal step |
| `simulation_pathway_tau.csv` | Secondary | Pathway ranking stability |
| `simulation_auc.csv` | Secondary | Classification AUC per step |
| `ST001865_perturbation_results.csv` | Tertiary | Perturbation validation |

---

## Implementation notes

**ST000356 format.** Non-standard transposed format with metadata rows. Handled
via `fix_st000356()` in all scripts.

**Condition label overrides.** Required for non-standard datasets:
```python
CASE_CTRL_OVERRIDE = {
    "ST001849": ("severe", "mild"),
    "ST002829": ("severe", "mild"),
    "ST001865": ("Hypoxia", "Normoxia"),
}
```

**FGW regularization.** All results use ε = 0.5 and fixed-scale normalization
(`M /= sqrt(d)`, not `M /= M.max()`). Keep `SINK_MAX_ITER` at 300; higher values
give no accuracy gain and cause multi-hour runtimes.

**node_features convention.** Returns `(n_nodes, 2)` — mean and standard deviation
per node — so that `cdist(X_s, X_t)` receives matching dimensions regardless of
case/control sample-count imbalance.

**Hard subset.** Always filter to `n_hidden ≥ 2` for final evaluation.
Single-hidden-metabolite trials are near-trivial and inflate apparent performance.

**Seeded run.** Use `run_ST001865_perturbation_seeded.py` for all results
reported in the manuscript (`GLOBAL_SEED = 42`).

---

## Reproducibility checklist

- [ ] Environment created and dependencies installed
- [ ] Raw cohort CSVs in `results/data/`
- [ ] `core_pathway_mapping.csv` in `results/data/`
- [ ] `CASE_CTRL_OVERRIDE` set correctly per dataset
- [ ] `SINK_MAX_ITER = 300` in all baseline scripts
- [ ] Hard subset (`n_hidden ≥ 2`) used for final evaluation
- [ ] `np.random.seed(42)` confirmed in ST001865 seeded script

---

## Citation

> Anas Enoch. *Commutator-Based Prioritization of Metabolite Measurements
> under Partial Observability.* Manuscript in preparation, 2026.

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

If scripts and CSV outputs disagree, the scripts in `results/scripts/active/`
are the source of truth.
