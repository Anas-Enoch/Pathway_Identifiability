# Commutator-Based Prioritization of Metabolite Measurements under Partial Observability

**Reproducibility package for the manuscript:**

> Anas Enoch. *Commutator-Based Prioritization of Metabolite Measurements under
> Partial Observability.* Manuscript in preparation, 2026.

**Author:** Anas Enoch, MD — Mohammed VI University of Health Sciences (UM6SS), Casablanca, Morocco

-----

## Scientific context

Metabolomics experiments rarely measure all metabolites in a given pathway.
Annotation rates of 10–30% are typical in untargeted studies, and targeted panels cover
only a fraction of most KEGG or Reactome pathways. This incomplete coverage creates a
structural problem: multiple distinct biochemical configurations can remain consistent
with the same observed measurements, leaving pathway activity ambiguous.

Existing computational tools — enrichment analysis (MSEA, ORA), pathway activity scoring
(mummichog), and imputation methods — do not address the practical question that follows
from this ambiguity:

> **Given my current metabolite panel, which single additional metabolite should I measure
> to most reduce pathway-level structural ambiguity?**

This repository implements a computational framework that answers this question directly.
It ranks unmeasured metabolites by their expected contribution to resolving pathway
underdetermination, using an efficient operator-commutator surrogate evaluated against
a reveal-defined oracle across five public metabolomics cohorts.

-----

## Methodological scope and limitations

This is a **computational methods paper**. The following limitations apply to all results:

- **No prospective acquisition experiment was performed.** The framework is evaluated
  by masking metabolites from otherwise complete datasets and measuring how accurately
  surrogate rankings recover the oracle. No new experimental measurements were made.
- **No causal biological intervention claims are made.** Prioritization recommendations
  are computationally motivated hypotheses about ambiguity reduction, not experimentally
  confirmed biological findings.
- **The framework complements rather than replaces enrichment analysis.** MSEA and ORA
  answer *which pathways are altered*; this framework answers *which metabolite to measure
  next* under incomplete coverage.
- **Surrogate fidelity is not biological validation.** High agreement with the
  reveal-defined oracle means the commutator efficiently approximates the oracle
  re-solve; it does not mean the recommended metabolite is a biologically validated measurement target.

-----

## Repository structure

```
.
├── results/
│   ├── data/                   # Processed sample × metabolite matrices and pathway maps
│   │   ├── ST000356_transposed.csv            # Breast cancer (134 samples, 157 metabolites)
│   │   ├── ST003390_processed.csv             # Type 2 diabetes (300 samples)
│   │   ├── ST003506_serum_processed.csv       # Lymphedema (34 samples)
│   │   ├── ST001849_benchmark_ready.csv       # COVID-19 LC/MS — processed (322 samples)
│   │   ├── ST001849_pathway_mapping.csv
│   │   ├── ST001865_standard.csv              # Hypoxia vs normoxia (16 samples)
│   │   ├── ST001865_metabolite_mapping.csv
│   │   ├── ST002829_benchmark_ready.csv       # COVID-19 Metabolon — processed (609 samples)
│   │   ├── ST002829_pathway_mapping.csv
│   │   └── core_pathway_mapping.csv           # KEGG-derived pathway membership map
│   │
│   ├── scripts/
│   │   ├── active/             # All runnable benchmark scripts
│   │   │   ├── run_real_multi_pathway_benchmark.py    # PRIMARY benchmark
│   │   │   ├── run_jl_stability_benchmark.py          # Operator selection
│   │   │   ├── run_ST001849_benchmark.py              # External cohort ST001849
│   │   │   ├── run_ST002829_benchmark.py              # External cohort ST002829
│   │   │   ├── run_msea_mummichog_benchmark.py        # Enrichment proxy baselines
│   │   │   ├── run_extended_baselines_benchmark.py    # MI / OED / active-acq baselines
│   │   │   ├── run_acquisition_simulation.py          # Sequential acquisition simulation
│   │   │   ├── run_ST001865_perturbation_seeded.py    # Perturbation validation (seeded)
│   │   │   └── run_ST001865_perturbation.py           # Perturbation validation (unseeded)
│   │   └── prep/               # One-time data preparation
│   │       ├── parse_ST002829.py
│   │       ├── parse_st003390.py
│   │       └── build_ST002829_pathway_mapping.py
│   │
│   └── results/                # Precomputed benchmark outputs (CSVs)
│
├── figures/                    # Manuscript figures (Fig1–Fig10)
├── references.bib
├── requirements.txt
└── README.md
```

-----

## Benchmark architecture

The evidence in this manuscript is organized in five layers:

```
1. PRIMARY BENCHMARK
   Three-cohort oracle-recovery benchmark — core manuscript claims
   ST000356 · ST003390 · ST003506  →  run_real_multi_pathway_benchmark.py

2. EXTERNAL VALIDATION
   Cross-platform generalisability stress tests
   ST001849 · ST002829  →  run_ST001849_benchmark.py / run_ST002829_benchmark.py

3. PERTURBATION-ORIENTED VALIDATION
   Mechanistically anchored biological extension
   ST001865  →  run_ST001865_perturbation_seeded.py

4. ACQUISITION SIMULATION
   Sequential measurement-selection corroboration
   ST000356  →  run_acquisition_simulation.py

5. EXTENDED BASELINE BENCHMARKING
   Stronger competing methods
   ST000356 · ST003390 · ST003506  →  run_extended_baselines_benchmark.py
                                       run_msea_mummichog_benchmark.py
```

-----

## Cohorts

### Primary benchmark cohorts (three datasets)

These three cohorts constitute the core benchmark. Results across them are
aggregated into the 29,750 hard-subset evaluation instances reported in Table 1.

|Dataset |Condition                 |n  |Metabolites|Source                         |
|--------|--------------------------|---|-----------|-------------------------------|
|ST000356|Breast cancer vs control  |134|157        |Metabolomics Workbench PR000284|
|ST003390|Type 2 diabetes vs control|300|~500       |Metabolomics Workbench PR002101|
|ST003506|Lymphedema vs control     |34 |~200       |Metabolomics Workbench PR002152|

### External validation cohorts (two datasets)

These cohorts are evaluated as **distribution-shift stress tests**, not as primary
benchmark data. They test whether relative method ordering is preserved when the
platform, disease context, and metabolite vocabulary all differ from the primary cohorts.

|Dataset |Condition        |n  |Platform        |Source              |
|--------|-----------------|---|----------------|--------------------|
|ST001849|COVID-19 severity|322|Untargeted LC/MS|Sindelar et al. 2021|
|ST002829|COVID-19 severity|609|Metabolon       |Mathew et al. 2021  |

Both use PCA-derived severity labels and correlation-derived pathway structure
(no curated edges). Absolute performance is lower than primary cohorts; the key
finding is that relative method ranking and the operating-regime boundary are
reproduced across platforms.

### Perturbation-oriented validation (one dataset)

|Dataset |Condition          |n |Source                |
|--------|-------------------|--|----------------------|
|ST001865|Hypoxia vs normoxia|16|Metabolomics Workbench|

ST001865 is an **external public dataset**, not a prospective experiment generated
by the authors. See [Perturbation-Oriented Validation](#perturbation-oriented-validation) below.

-----

## Setup

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip setuptools wheel
python -m pip install numpy pandas scipy scikit-learn POT matplotlib
```

Run all commands from the **repository root**.

-----

## End-to-end reproduction

The following commands reproduce all reported computational results in execution order.

### Step 0 — Data preparation (run once; processed CSVs already in `results/data/`)

```bash
# Only required if re-processing raw downloads from Metabolomics Workbench.
# Processed benchmark-ready CSVs are already present in results/data/.
python results/scripts/prep/parse_ST002829.py
python results/scripts/prep/parse_st003390.py
python results/scripts/prep/build_ST002829_pathway_mapping.py
```

Additional cohort-specific preprocessing may be needed depending on raw dataset
format. ST001849 benchmark-ready CSVs are pre-generated and included in `results/data/`.

### Step 1 — Operator stability benchmark (operator selection)

```bash
python results/scripts/active/run_jl_stability_benchmark.py
```

Evaluates six preprocessing operators under three perturbation families.
Selects `METHOD_DEFAULT="none"` as the stable operator for the primary benchmark.
Output: `results/results/jl_stability_benchmark.csv`

### Step 2 — Primary multi-cohort benchmark

```bash
# Reproduces Table 1 and Figures 5–7
python results/scripts/active/run_real_multi_pathway_benchmark.py
```

Core benchmark across ST000356, ST003390, ST003506. For each masking trial,
computes the oracle metabolite `m* = argmax_m ΔU_k(m)` and evaluates all
competing predictors against it.
Output: `results/results/real_multi_pathway_results.csv`

Inspect results:

```python
import pandas as pd
df   = pd.read_csv("results/results/real_multi_pathway_results.csv")
hard = df[df["n_hidden"] >= 2]          # hard subset: ≥2 hidden metabolites per trial
print(
    hard.groupby("predictor_method")
        [["regret","nregret","top1","top3","rank_tau"]]
        .mean().round(4)
        .sort_values("rank_tau", ascending=False)
        .to_string()
)
```

### Step 3 — External cohort benchmarks

```bash
# Reproduces Table 3 and Figure 8
python results/scripts/active/run_ST001849_benchmark.py
python results/scripts/active/run_ST002829_benchmark.py
```

Outputs: `ST001849_benchmark_results.csv`, `ST002829_benchmark_results.csv`

### Step 4 — Extended baseline benchmarks

```bash
# Reproduces MI / Bayesian OED / active-acquisition rows of Table 1
python results/scripts/active/run_extended_baselines_benchmark.py

# Reproduces MSEA proxy / mummichog proxy rows of Table 1
python results/scripts/active/run_msea_mummichog_benchmark.py
```

Outputs: `extended_baselines_results.csv`, `msea_mummichog_results.csv`

### Step 5 — Acquisition simulation

```bash
# Reproduces Section 3.4
python results/scripts/active/run_acquisition_simulation.py
```

Sequential metabolite-acquisition simulation on ST000356 at 40% masking.
Compares commutator, variance, and random selection over 50 trials.
Outputs: `simulation_delta_uk.csv`, `simulation_pathway_tau.csv`, `simulation_auc.csv`

### Step 6 — Perturbation-oriented validation

```bash
# Reproduces Section 3.6 — use the seeded script for manuscript-exact results
python results/scripts/active/run_ST001865_perturbation_seeded.py
```

Output: `results/results/ST001865_perturbation_results.csv`

Always use `run_ST001865_perturbation_seeded.py` (sets `GLOBAL_SEED=42`,
`np.random.seed(42)`) for results matching the manuscript. The unseeded variant
`run_ST001865_perturbation.py` produces variable τ values across runs due to
FGW internal initialization.

### Step 7 — Figure generation

Figures were generated from the CSV outputs above using custom matplotlib scripts
in `results/scripts/active/`. The pregenerated figures are in `figures/` and are
included in the submission ZIP. To regenerate, run the compute scripts in Steps 1–6
and pass the resulting CSVs to the figure scripts.

-----

## Expected outputs

### Combined manuscript benchmark summary

Primary cohorts combined, hard subset (n_hidden ≥ 2), 29,750 evaluation instances:

|Group      |Predictor        |Regret   |Top-1    |Top-3    |Kendall τ|
|-----------|-----------------|---------|---------|---------|---------|
|Ours       |`gnc_commutator` |**0.018**|**0.774**|**0.821**|**0.599**|
|Ours       |`surrogate`      |0.024    |0.631    |0.728    |0.401    |
|Info./Bayes|`bayes_oed`      |0.165    |0.568    |0.697    |0.100    |
|Info./Bayes|`active_acq`     |0.159    |0.516    |0.702    |0.017    |
|Info./Bayes|`mutual_info`    |0.257    |0.372    |0.696    |−0.517   |
|Heuristic  |`mummichog_proxy`|0.147    |0.574    |0.703    |0.204    |
|Heuristic  |`variance`       |0.048    |0.511    |0.663    |0.086    |
|Heuristic  |`diffabundance`  |0.051    |0.489    |0.641    |0.074    |
|Heuristic  |`msea_proxy`     |0.216    |0.465    |0.697    |−0.030   |
|Heuristic  |`degree`         |0.083    |0.371    |0.541    |0.063    |
|Heuristic  |`random`         |0.091    |0.312    |0.514    |0.000    |

Three-tier τ hierarchy: operator-geometric (0.40–0.60) > connectivity-weighted
differential (0.10–0.20) > label-marginal (−0.52 to 0.09).

### External validation

|Cohort                     |Comm top-1|Var top-1|Comm τ|Var τ|
|---------------------------|----------|---------|------|-----|
|ST001849 (LC/MS, n=322)    |0.382     |0.292    |0.243 |0.116|
|ST002829 (Metabolon, n=609)|0.348     |0.255    |0.172 |0.072|

-----

## Perturbation-oriented validation

ST001865 is a public untargeted metabolomics dataset from Metabolomics Workbench
contrasting hypoxia (n=8) with normoxia (n=8) in a matched cell-line design.
**It is not a prospective experiment generated by the authors.** It is used as an
external perturbation dataset to test whether the acquisition framework prioritizes
structurally informative metabolites under biologically controlled metabolic rewiring.

Hypoxia induces coordinated rewiring of arginine and nitrogen metabolism. After
intersection with the pathway map, four pathways met the minimum metabolite coverage
threshold. Three had only one hidden metabolite under any masking rate (forced
selection — all strategies trivially equivalent). Only **Arginine and Proline
Metabolism** (6 metabolites, 2–3 hidden) produced a non-trivial benchmark.

**Main result (60 trials, seeded, GLOBAL_SEED=42):**

|Strategy        |Top-1    |Kendall τ|
|----------------|---------|---------|
|`gnc_commutator`|**0.833**|**0.709**|
|`random`        |0.617    |0.233    |
|`diffabundance` |0.467    |0.067    |
|`variance`      |0.283    |−0.178   |

**Citrulline recovery:** the reveal-defined oracle selected Citrulline — a structurally
central metabolite at the arginine/urea cycle interface — in 45% of trials (27/60).
When oracle = Citrulline, the commutator recovered that choice in 92.6% of cases (25/27),
vs 3.7% for variance and differential abundance (Fisher exact OR = 325, p < 0.001).

This result illustrates the conceptual distinction between abundance-based and
structural prioritization: the commutator ranks by structural effect on pathway
observability; abundance methods rank by marginal signal magnitude. Citrulline
separates these two philosophies sharply in this dataset.

> **Scope:** the perturbation result is concentrated in a single pathway due to
> sparse metabolite coverage in ST001865. This is a coverage limitation, not a
> method limitation.

-----

## Known operating boundary

The commutator advantage is pathway-size dependent:

- **Preserved** in pathways with approximately |M_k| ≤ 40 nodes.
- **Degrades** in large lipid-dominated pathways (|M_k| > 40), where high
  inter-metabolite correlation flattens the correlation-distance structure matrix.
  When all metabolites appear equidistant, the commutator’s operator-geometry signal
  saturates and variance-based baselines become competitive.

This boundary is reproduced quantitatively (δτ ≈ −0.006 at |M_k| > 40) across
both external cohorts on different platforms (LC/MS and Metabolon), confirming it
as a regime property rather than a dataset artefact.

**Predicted fix:** providing curated biochemical edge topology from LIPID MAPS or
Reactome — rather than correlation-derived structure matrices — is expected to restore
commutator discriminability in the large-pathway regime. This is a directly testable
prediction requiring no new experimental data.

-----

## Implementation notes

**Hard subset.** Always filter to `n_hidden ≥ 2` for evaluation. Single-hidden-metabolite
trials are near-trivial (any method achieves perfect recovery by chance) and inflate
apparent performance.

**FGW parameters.** All results use entropic regularization ε = 0.5 and fixed-scale
cost normalization (`M /= sqrt(d)`). Do not use `M /= M.max()` — this destroys scale
information needed for operator comparison. Keep `SINK_MAX_ITER = 300`; higher values
provide no accuracy gain and cause multi-hour runtimes.

**node_features convention.** Returns `(n_nodes, 2)` [mean, std] per node, ensuring
`cdist(X_s, X_t)` receives matching feature dimensions regardless of case/control
sample-count imbalance.

**ST000356 format.** Stored in a non-standard transposed format with metadata rows.
Handled automatically via `fix_st000356()` in all benchmark scripts.

**Condition label overrides.** Required for non-standard label strings:

```python
CASE_CTRL_OVERRIDE = {
    "ST001849": ("severe", "mild"),
    "ST002829": ("severe", "mild"),
    "ST001865": ("Hypoxia", "Normoxia"),
}
```

**Seeded perturbation script.** Always use `run_ST001865_perturbation_seeded.py`
for manuscript results. `run_ST001865_perturbation.py` is available for exploratory
use but produces variable τ values between runs.

-----

## Output files reference

|File                               |Step|Contents                                     |
|-----------------------------------|----|---------------------------------------------|
|`jl_stability_benchmark.csv`       |1   |Operator stability — CV, drift, τ by method  |
|`real_multi_pathway_results.csv`   |2   |Primary benchmark — all predictors, 3 cohorts|
|`ST001849_benchmark_results.csv`   |3   |External validation — ST001849               |
|`ST002829_benchmark_results.csv`   |3   |External validation — ST002829               |
|`extended_baselines_results.csv`   |4   |MI / OED / active-acquisition results        |
|`msea_mummichog_results.csv`       |4   |Enrichment proxy baselines                   |
|`simulation_delta_uk.csv`          |5   |ΔU_k per reveal step per pathway             |
|`simulation_pathway_tau.csv`       |5   |Pathway ranking stability (Kendall τ)        |
|`simulation_auc.csv`               |5   |Classification AUC per reveal step           |
|`ST001865_perturbation_results.csv`|6   |Perturbation validation — seeded run         |

-----

## Reproducibility checklist

- [ ] Dependencies installed in `.venv`
- [ ] All datasets present in `results/data/`
- [ ] `core_pathway_mapping.csv` present in `results/data/`
- [ ] `SINK_MAX_ITER = 300` in all benchmark scripts
- [ ] `n_hidden >= 2` filter applied before evaluation
- [ ] `run_ST001865_perturbation_seeded.py` used for Section 3.6 results

-----

## Figures

|File                                   |Figure|Description                                    |
|---------------------------------------|------|-----------------------------------------------|
|`Fig1_revised_pipeline.png`            |1     |End-to-end measurement-prioritization pipeline |
|`Fig2_Structural_Ambiguity.pdf`        |2     |Structural ambiguity under partial metabolomics|
|`Fig3_operator_stability_heatmap.png`  |3     |Operator stability under perturbation regimes  |
|`Fig4_jl_vs_randproj_validation.png`   |4     |Fixed-scale normalization validation           |
|`Fig5_global_hard_subset_benchmark.png`|5     |Global oracle-recovery performance             |
|`Fig6_per_dataset_benchmark.png`       |6     |Per-dataset commutator vs variance             |
|`Fig7_pathway_heterogeneity.png`       |7     |Pathway-level performance heterogeneity        |
|`Fig8_commutator_mechanism.png`        |8     |External validation across five cohorts        |
|`Fig9_Citrulline_validation.png`       |9–10  |Commutator mechanism and Citrulline recovery   |

-----

## Citation

```
Anas Enoch. Commutator-Based Prioritization of Metabolite Measurements under
Partial Observability. Manuscript in preparation, 2026.
```

External cohort citations:

> Sindelar M et al. *Longitudinal metabolomics of human plasma reveals prognostic
> markers of COVID-19 disease severity.* Cell Reports Medicine 2(8):100369, 2021.
> Metabolomics Workbench ST001849.

> Mathew D et al. *Nucleotide, phospholipid, and kynurenine metabolites are robustly
> associated with COVID-19 severity.* Metabolomics Workbench ST002829, 2021.

-----

## Contact

**Anas Enoch, MD**
Mohammed VI University of Health Sciences (UM6SS), Casablanca, Morocco

Scripts in `results/scripts/active/` are the source of truth if outputs disagree.