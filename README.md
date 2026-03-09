# Pathway Identifiability under Partial Metabolomics

This repository contains code, processed data tables, pathway mappings, and benchmark outputs for a framework that studies **pathway identifiability under incomplete metabolite observation**.

The core idea is simple:

> when metabolomics coverage is partial, pathway interpretation can become structurally ambiguous; this repository provides a benchmark framework to quantify that ambiguity and prioritize additional metabolite measurements that reduce it.

---

## What this repository currently contains

The repository includes:

- processed metabolomics matrices for public cohorts
- an expanded pathway mapping file
- cross-dataset pathway coverage summaries
- retained-pathway benchmark tables
- figures used in the manuscript
- scripts for preprocessing, overlap screening, benchmark generation, and plotting

---

## Current benchmark scope

The current benchmark is dataset-conditioned and includes three public metabolomics cohorts:

- **ST000356**
- **ST003390**
- **ST003506**

The benchmark currently focuses on pathway coverage, retained-pathway analysis, and measurement-oriented identifiability benchmarking under partial observation.

A key result already visible in the current outputs is that pathway observability is **not fixed across datasets**. Some pathways are reproducibly covered across cohorts, whereas others are strongly cohort-dependent.

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
│   ├── fig_pathway_coverage_clean.png
│   └── fig_retained_pathways_heatmap.png
│
└── scripts/
    ├── run_core_benchmark.py
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

Main result files

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

⸻

Main figures

results/results/fig_pathway_coverage_clean.png

Pathway coverage plot for the benchmark cohorts.

results/results/fig_retained_pathways_heatmap.png

Cross-dataset heatmap showing pathway coverage heterogeneity across ST000356, ST003390, and ST003506.

⸻

Scripts

results/scripts/run_core_benchmark.py

Runs the core pathway identifiability benchmark.

results/scripts/plot_pathway_coverage.py

Generates the pathway coverage figure.

results/scripts/screen_st003390_overlap.py

Screens metabolite overlap between ST003390 and the current pathway mapping.

results/scripts/build_retained_pathways_by_dataset.py

Builds the cross-dataset retained-pathway summary table.

results/scripts/plot_retained_pathways_heatmap.py

Generates the cross-dataset pathway coverage heatmap.

results/scripts/parse_st003390.py

Extracts the ST003390 metabolite block from the raw source file.

results/scripts/clean_st003390.py

Converts ST003390 into the benchmark-ready processed matrix.

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
	•	manuscript figures updated
	•	preprocessing and plotting scripts tracked in the repository

Still planned:
	•	broader comparator benchmark against stronger baseline families
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
