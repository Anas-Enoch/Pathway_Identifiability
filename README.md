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

⸻
How to run the current workflow

1. Run the core benchmark
python3 results/scripts/run_core_benchmark.py

2. Generate the pathway coverage figure
python3 results/scripts/plot_pathway_coverage.py

3. Screen overlap for ST003390
python3 results/scripts/screen_st003390_overlap.py

4. Build the retained-pathway summary
python3 results/scripts/build_retained_pathways_by_dataset.py

5. Generate the cross-dataset heatmap
python3 results/scripts/plot_retained_pathways_heatmap.py

6. Run the baseline benchmark
python3 results/scripts/run_baseline_benchmark.py
This script generates:
	•	results/results/baseline_trial_results.csv
	•	results/results/baseline_regret_summary.csv
	•	results/results/baseline_topk_summary.csv

Current biological interpretation

The current benchmark already shows:
	•	Arginine and Proline Metabolism is reproducibly retained across all three cohorts
	•	Glycine Serine Threonine Metabolism also remains recurrent across datasets
	•	Alanine Aspartate Glutamate Metabolism and Valine Leucine Isoleucine Degradation are usable in multiple cohorts
	•	Pyrimidine Metabolism is strong in some datasets and absent in others
	•	Glycolysis / Gluconeogenesis and TCA Cycle remain highly dataset-dependent

This is exactly why the framework is formulated as a dataset-conditioned identifiability benchmark, not a one-size-fits-all pathway scoring method.

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
