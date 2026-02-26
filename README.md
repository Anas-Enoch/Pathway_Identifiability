# Pathway Identifiability under Partial Metabolomics

## Quantifying Structural Ambiguity for Measurement Prioritization

This repository contains the code and processed data supporting the
manuscript:

**"Quantifying Pathway Identifiability under Partial Metabolomics for
Measurement Prioritization."**

The study introduces a structural instability framework for pathway
ranking under partial metabolite observability and proposes the **Uk
score** as a structural ambiguity metric for measurement prioritization.

------------------------------------------------------------------------

# Biological Dataset

Metabolomics data analyzed in this study are publicly available at:

Metabolomics Workbench (NIH Common Fund NMDR)\
Project ID: PR002152\
Study ID: ST003506\
DOI: https://doi.org/10.21228/M8FR6S

Title:\
*Comparison of serum and interstitial fluid from patients with breast
cancer-related lymphedema and healthy control subjects with NMR-based
metabolomics.*

This repository contains processed serum metabolomics matrices derived
from the above dataset.

------------------------------------------------------------------------

# Repository Structure

    github.pathway/
    │
    ├── codes/                  # Figure generation scripts
    ├── figures/                # Manuscript and supplementary figures (PDF)
    ├── results/
    │   ├── data/               # Processed biological data and derived metrics
    │   ├── results/            # Synthetic benchmarking outputs
    │   └── scripts/            # Instability analysis and benchmarking scripts
    └── requirements.txt

------------------------------------------------------------------------

# Data Description

### processed_metabolite_matrix.csv

Cleaned quantitative metabolite matrix derived from NMR serum data.

### pathway_mapping.csv

Mapping between measured metabolites and curated metabolic pathway sets.

### Uk_scores.csv

Structural ambiguity score (Uk) computed for each pathway.

### ranking_stability_by_method.csv

Replicate-level Kendall τ ranking stability values.

### pathway_instability_vs_Uk.csv

Per-pathway instability index and corresponding Uk score values.

### multi_pathway_results.csv

Synthetic multi-pathway benchmarking results evaluating regret under
masking regimes.

------------------------------------------------------------------------

# Reproducing Main Figures

Install dependencies:

    pip install -r requirements.txt

Example figure generation:

    python codes/make_fig_pathway_ranking_instability.py

------------------------------------------------------------------------

# Author

Anas Enoch, MD\
Mohammed VI University of Health Sciences (UM6SS)\
Casablanca, Morocco\
anas_enoch@um5.ac.ma
