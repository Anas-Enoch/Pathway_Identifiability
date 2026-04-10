# Pathway Identifiability under Partial Metabolomics

This repository contains code, processed data tables, pathway mappings, and benchmark outputs for a framework that studies **pathway identifiability under incomplete metabolite observation**.

The central question is:

> under partial metabolomics, which additional metabolite measurements are expected to reduce pathway-level ambiguity most effectively?

The repository therefore focuses on **measurement prioritization**, **retained-pathway analysis**, and **dataset-conditioned identifiability benchmarking**, rather than on generic pathway enrichment alone.

---

## Current benchmark scope

The current benchmark includes three public metabolomics cohorts:

- **ST000356**
- **ST003390**
- **ST003506**

A key result already visible in the current outputs is that pathway observability is **not fixed across datasets**. Some pathways are reproducibly retained across cohorts, whereas others are strongly cohort-dependent.

---

## Repository structure

```text
results/
├── data/
│   ├── core_pathway_mapping.csv
│   ├── processed_metabolite_matrix_ST000356.csv
│   ├── ST003390_processed.csv
│   ├── ST003506_serum_processed.csv
│   └── raw/
│
├── results/
│   ├── core_benchmark_results.csv
│   ├── pathway_scores_summary.csv
│   ├── pathway_coverage_summary.csv
│   ├── pathway_coverage_summary_ST003390.csv
│   ├── retained_pathways_by_dataset.csv
│   ├── baseline_trial_results.csv
│   ├── baseline_regret_summary.csv
│   ├── baseline_topk_summary.csv
│   ├── fig_pathway_coverage_clean.png
│   └── fig_retained_pathways_heatmap.png
│
└── scripts/
    ├── run_core_benchmark.py
    ├── run_baseline_benchmark.py
    ├── plot_pathway_coverage.py
    ├── screen_st003390_overlap.py
    ├── build_retained_pathways_by_dataset.py
    ├── plot_retained_pathways_heatmap.py
    ├── parse_st003390.py
    └── clean_st003390.py

Main data files

results/data/core_pathway_mapping.csv

Expanded pathway mapping used for cross-dataset coverage and identifiability screening.

results/data/processed_metabolite_matrix_ST000356.csv

Processed metabolite matrix for cohort ST000356.

results/data/ST003390_processed.csv

Processed metabolite matrix for cohort ST003390.

results/data/ST003506_serum_processed.csv

Processed metabolite matrix for cohort ST003506.

⸻

Main benchmark summary files

results/results/pathway_coverage_summary.csv

Pathway coverage summary for the original benchmark cohorts.

results/results/pathway_coverage_summary_ST003390.csv

Coverage summary for ST003390 against the current pathway mapping.

results/results/retained_pathways_by_dataset.csv

Cross-dataset summary of retained pathways, detected metabolite counts, and benchmark inclusion flags.

results/results/core_benchmark_results.csv

Core benchmark outputs for ambiguity-reduction / regret-style evaluation.

results/results/pathway_scores_summary.csv

Pathway-level summary outputs from the current benchmark pipeline.

results/results/baseline_trial_results.csv

Trial-level benchmark output for baseline comparison across the three cohorts.

results/results/baseline_regret_summary.csv

Aggregated regret summary across pathways, datasets, masking rates, and methods.

results/results/baseline_topk_summary.csv

Aggregated top-k success summary across pathways, datasets, masking rates, and methods.

⸻

Main figures

results/results/fig_pathway_coverage_clean.png

Pathway coverage plot for the benchmark cohorts.

results/results/fig_retained_pathways_heatmap.png

Cross-dataset heatmap showing pathway coverage heterogeneity across ST000356, ST003390, and ST003506.

⸻

Current biological interpretation

The current benchmark shows that:
	•	Arginine and Proline Metabolism is reproducibly retained across all three cohorts
	•	Glycine Serine Threonine Metabolism also remains recurrent across datasets
	•	Alanine Aspartate Glutamate Metabolism and Valine Leucine Isoleucine Degradation are usable in multiple cohorts
	•	Pyrimidine Metabolism is strong in some datasets and absent in others
	•	Glycolysis / Gluconeogenesis and TCA Cycle remain highly dataset-dependent

This is why the framework is formulated as a dataset-conditioned identifiability benchmark, not a one-size-fits-all pathway scoring method.

⸻

Baseline benchmark

The repository includes a reproducible cross-dataset benchmark comparing the proposed prioritization rule against several baseline strategies.

Methods currently compared
	•	proposed
	•	imputation
	•	degree
	•	random

Benchmark design
	•	three public metabolomics cohorts
	•	retained pathways only
	•	repeated masking experiments
	•	regret-based evaluation
	•	top-k recovery evaluation

Current benchmark interpretation

In the current three-cohort benchmark, the proposed structural prioritization rule clearly outperforms random and degree-based baselines and improves upon the imputation-based comparator. The proposed and imputation methods still agree in many trials, but the structural rule diverges in a non-trivial subset of cases and yields lower regret overall.


## Current Benchmark and Ablation Status

This repository contains both real-data instability analyses and synthetic benchmarking utilities.

### Real-data analyses
The following outputs are derived from real metabolomics cohorts and are used for pathway-ranking instability analysis under controlled masking:

- `ranking_stability_by_method.csv`
- `pathway_instability_vs_Uk.csv`
- processed cohort matrices in `data/` and `results/Data/`

These files support the manuscript’s real-data instability analyses.

### Synthetic benchmarking utilities
Several scripts in `results/scripts/` generate synthetic or simulated benchmarking outputs for stress testing, regret visualization, and pipeline debugging. In particular:

- `generate_multi_pathway_results.py`
- `generate_large_synthetic_dataset.py`

These scripts are intended for controlled benchmarking and figure generation, not as independent biological ground-truth validation.

### Important note on component ablations
The manuscript defines the composite underdetermination score as:

`U_k = α H(T) + β Var(A) + γ SGI`

where:
- `H(T)` = transport entropy
- `Var(A)` = alignment instability
- `SGI` = structural growth index

The current repository includes an ablation pipeline (`run_component_ablation.py`) and summary outputs:

- `ablation_summary.csv`
- `ablation_scores_by_pathway.csv`

However, unless explicitly regenerated from the full JL-FGW pathway pipeline, the current `Uk_components.csv` file should be interpreted as a provisional component decomposition used to validate the ablation workflow, not as the final biologically grounded export of `H(T)`, `Var(A)`, and `SGI`.

Accordingly, current ablation outputs should be interpreted as diagnostic support for component-correlation analysis rather than definitive evidence of independent biological contribution of each component. A fully grounded ablation requires direct export of the three component terms from the real `compute_Uk(...)` implementation.

⸻
## How to run the current workflow

### 1. Run the core benchmark
python3 results/scripts/run_core_benchmark.py

### 2. Generate the pathway coverage figure
python3 results/scripts/plot_pathway_coverage.py

### 3. Screen overlap for ST003390
python3 results/scripts/screen_st003390_overlap.py

### 4. Build the retained-pathway summary
python3 results/scripts/build_retained_pathways_by_dataset.py

### 5. Generate the cross-dataset heatmap
python3 results/scripts/plot_retained_pathways_heatmap.py

### 6. Run the baseline benchmark
python3 results/scripts/run_baseline_benchmark.py

### 9. Build pathway-level feature tables for real JL-FGW benchmarking
python3 results/scripts/compute_pathway_features.py

### 10. Compute real FGW alignments under preprocessing variants
python3 results/scripts/compute_fgw_alignment.py

### 11. Export real Uk components
python3 results/scripts/compute_Uk_real.py

### 12. Run the real multi-pathway benchmark
python3 results/scripts/run_real_multi_pathway_benchmark.py

### 13. Run the JL stability benchmark
python3 results/scripts/run_jl_stability_benchmark.py

This script generates:

- `results/results/baseline_trial_results.csv`
- `results/results/baseline_regret_summary.csv`
- `results/results/baseline_topk_summary.csv`

### 7. Compute pathway ranking instability under controlled masking

This evaluates how pathway rankings change when metabolites are randomly masked.

python3 results/scripts/generate_pathway_instability_csv.py

Outputs:

- `ranking_stability_by_method.csv`
- `pathway_instability_vs_Uk.csv`

### 8. Run component ablation analysis

This evaluates the relative contribution of the three components of the composite identifiability score:

U_k = α H(T) + β Var(A) + γ SGI

python3 results/scripts/run_component_ablation.py

Outputs:

- `results/Results/ablation_summary.csv`
- `results/Results/ablation_scores_by_pathway.csv`

## Current biological interpretation

Across the evaluated metabolomics cohorts, several pathways appear recurrently identifiable under the current benchmark configuration:

- **Arginine and Proline Metabolism** is consistently retained across all three cohorts.
- **Glycine, Serine, and Threonine Metabolism** also appears recurrent across datasets.
- **Alanine, Aspartate, and Glutamate Metabolism** and **Valine, Leucine, and Isoleucine Degradation** are usable in multiple cohorts.
- **Pyrimidine Metabolism** shows strong signal in some datasets but not others.
- **Glycolysis/Gluconeogenesis** and the **TCA cycle** remain strongly dataset-dependent.

This variability is expected: metabolomics coverage differs substantially between cohorts. For this reason the framework is formulated as a **dataset-conditioned identifiability benchmark**, rather than a universal pathway scoring system.
⸻

Manuscript-facing interpretation

This repository supports a manuscript that makes the following narrower but stronger claim:

the method is designed for measurement prioritization and pathway disambiguation under partial metabolomics, rather than as a generic replacement for all enrichment-based pathway analysis methods.

The benchmark therefore emphasizes:
	•	pathway coverage under incomplete observation
	•	retained-pathway analysis
	•	ambiguity reduction
	•	dataset-conditioned identifiability
	•	benchmark comparison against explicit baseline strategies

⸻

Public data provenance

This repository uses secondary analysis of public metabolomics datasets from Metabolomics Workbench.

The current benchmark includes:
	•	ST000356
	•	ST003390
	•	ST003506

Processed benchmark tables in this repository are derived from those public datasets.

⸻

Status

Current status of the repository:
	•	pathway mapping expanded
	•	three-cohort benchmark assembled
	•	cross-dataset retained-pathway summary generated
	•	baseline benchmark pipeline added
	•	benchmark outputs committed
	•	manuscript figures updated
	•	preprocessing and plotting scripts tracked in the repository

Still planned:
	•	broader comparator benchmark against stronger pathway-analysis families
	•	addition of at least one more cohort with stronger central-carbon coverage
	•	expanded cross-dataset evaluation

⸻
Requirements

Install dependencies with:
pip install -r requirements.txt

If needed, create a virtual environment first:
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt


# Author

Anas Enoch, MD\
Mohammed VI University of Health Sciences (UM6SS)\
Casablanca, Morocco\
anas_enoch@um5.ac.ma
