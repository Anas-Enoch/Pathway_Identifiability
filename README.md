# Pathway Identifiability under Partial Metabolomics

This repository implements a **dataset-conditioned metabolomic pathway identifiability framework**.

The central question is:

> under partial metabolomics, which additional metabolite measurements are expected to reduce pathway-level ambiguity most effectively?

The framework is built around four linked components:

1. **pathway-level feature construction**
2. **cross-condition pathway alignment with fused Gromov–Wasserstein (FGW)**
3. **stability benchmarking under structured perturbations**
4. **non-trivial metabolite prioritization against an exact reveal oracle**

The repository is therefore **not** a generic pathway enrichment pipeline. It is a **stability-gated identifiability benchmark**.

---

## Current repository layout

```text
.
├── codes/
├── figures/
├── results/
│   ├── data/
│   ├── results/
│   └── scripts/
├── README.md
├── references.bib
└── requirements.txt

Main directories
	•	results/data/
	•	processed cohort-level metabolomics tables
	•	pathway mapping file used for pathway restriction
	•	results/scripts/
	•	executable Python scripts for each pipeline stage
	•	results/results/
	•	benchmark outputs and intermediate result tables
	•	figures/
	•	manuscript figures
	•	codes/
	•	figure-generation assets or exported figure-related files

⸻
Cohorts currently used

The current benchmark uses three public metabolomics cohorts:
	•	ST000356
	•	ST003390
	•	ST003506

These datasets are processed into condition-aware metabolite matrices and then restricted to pathway-specific metabolite subsets.

⸻

Conceptual pipeline

The implemented pipeline is:
Raw metabolomics
→ condition parsing
→ pathway restriction
→ node feature construction
→ optional preprocessing operator
→ FGW alignment
→ pathway identifiability functional U_k
→ perturbation-based stability benchmark
→ metabolite prioritization benchmark against exact oracle

What the pipeline does

For each pathway:
	1.	split the dataset into case and control
	2.	restrict to metabolites observed in that pathway
	3.	build condition-specific node feature matrices
	4.	construct a pathway structure matrix
	5.	apply a selected operator / preprocessing method
	6.	align case vs control with FGW
	7.	compute a pathway identifiability score (U_k)
	8.	evaluate stability under perturbations
	9.	run a masking benchmark to determine which hidden metabolite most reduces ambiguity when revealed


⸻

Main scripts

1. Build pathway-level features
python3 results/scripts/active/compute_pathway_features.py

This script:
	•	loads each cohort
	•	parses the condition labels
	•	handles the special ST000356 transposed format
	•	restricts metabolites to each pathway
	•	writes pathway-specific condition feature files and structure files

2. Compute FGW alignments
python3 results/scripts/active/compute_fgw_alignment.py

This script:
	•	reads the pathway feature files
	•	aligns case and control states for each pathway
	•	applies one of several preprocessing operators
	•	writes alignment outputs

3. Compute pathway-level identifiability scores

```bash
python3 results/scripts/active/compute_Uk_real.py

This script computes the pathway-level identifiability functional (U_k) from the aligned pathway state.

In the current implementation, (U_k) is derived from:
	•	the FGW transport plan
	•	transport entropy
	•	a structural spectral-gap term computed from the pathway structure matrix

This script provides the scoring layer used by downstream benchmarking and analysis.

4. Run the stability benchmark
python3 results/scripts/active/run_jl_stability_benchmark.py

This script evaluates operator stability under perturbation regimes.

Implemented perturbation families include:
	•	noise
	•	dropout
	•	bootstrap

Implemented operator / preprocessing methods include:
	•	none
	•	l2
	•	jl
	•	randproj
	•	pca_fixed
	•	pca_var95

The benchmark outputs pathway-level metrics such as:
	•	cv_fgw
	•	transport_drift
	•	delta_U
	•	top3_jaccard
	•	rank_tau

This stage is used as an operator-selection layer: unstable operators should not be carried forward into the downstream identifiability benchmark.

5. Run the real multi-pathway benchmark
python3 results/scripts/run_real_multi_pathway_benchmark.py
This script runs the non-trivial metabolite-prioritization benchmark.

For each pathway and masking replicate, it:
	1.	masks a subset of metabolites
	2.	computes the exact oracle by full reveal:
[
\Delta U_k(m)=U_k(\text{masked})-U_k(\text{masked}+m)
]
	3.	ranks hidden metabolites by exact reveal to define the oracle
	4.	evaluates non-oracle predictors against that oracle

Oracle benchmark design

The benchmark is only meaningful when the predictor is not identical to the oracle.

The current benchmark includes multiple predictors:
	•	gnc_commutator
	•	surrogate
	•	variance
	•	diffabundance
	•	degree
	•	random
	•	mb2d_transport (exploratory / transport-inspired)

Core evaluation metrics

For each predictor:
	•	regret
	•	nregret
	•	top1
	•	top3
	•	rank_tau

Recommended evaluation subset

The benchmark is most informative on the non-trivial masking subset: n_hidden >= 2
Single-hidden-metabolite cases are much easier and can inflate apparent performance.


⸻

Recommended workflow

Run the pipeline from the repository root.

Step 1 — create environment
python3 -m venv .env
source .env/bin/activate

Step 2 — install dependencies
pip install -r requirements.txt
If needed, install explicitly: pip install numpy pandas scipy scikit-learn POT tqdm

Step 3 — build features
python3 results/scripts/active/compute_pathway_features.py

Step 4 — compute pathway-level U_k
python3 results/scripts/active/compute_Uk_real.py

Step 5 — compute alignments
python3 results/scripts/active/compute_fgw_alignment.py

Step 6 — run stability benchmark
python3 results/scripts/active/run_jl_stability_benchmark.py


Step 7 — inspect operator stability
python3 - <<'PY'
import pandas as pd
df = pd.read_csv("results/results/jl_stability_benchmark.csv")
print(
    df.groupby("method")[["cv_fgw","transport_drift","delta_U","top3_jaccard","rank_tau"]]
      .mean()
      .round(4)
      .sort_values("cv_fgw")
      .to_string()
)
PY

Step 8 — run the real multi-pathway benchmark
python3 results/scripts/active/run_real_multi_pathway_benchmark.py

Step 9 — summarize the hard subset
python3 - <<'PY'
import pandas as pd
df = pd.read_csv("results/results/real_multi_pathway_results.csv")
hard = df[df["n_hidden"] >= 2].copy()

print("Rows:", len(hard))
print(
    hard.groupby("predictor_method")[["regret","nregret","top1","top3","rank_tau"]]
        .mean()
        .round(4)
        .sort_values(["regret","top1"], ascending=[True, False])
        .to_string()
)
PY

Step 10 — inspect per-dataset behavior
python3 - <<'PY'
import pandas as pd
df = pd.read_csv("results/results/real_multi_pathway_results.csv")
hard = df[df["n_hidden"] >= 2].copy()

print(
    hard.groupby(["dataset","predictor_method"])[["regret","nregret","top1","top3","rank_tau"]]
        .mean()
        .round(4)
        .to_string()
)
PY

Step 11 — inspect pathway heterogeneity
python3 - <<'PY'
import pandas as pd
df = pd.read_csv("results/results/real_multi_pathway_results.csv")
hard = df[df["n_hidden"] >= 2].copy()

print(
    hard.groupby(["dataset","pathway","predictor_method"])[["regret","nregret","top1","top3","rank_tau"]]
        .mean()
        .round(4)
        .to_string()
)
PY



⸻

Current interpretation of the pipeline

1. Operator stability comes first

The stability benchmark is not decorative. It is a selection step.

If an operator is unstable under perturbation, it should not be trusted in the downstream benchmark.

In the current low-node pathway regime, non-projection operators are the stable ones, whereas aggressive projection can be unstable.

2. The masking benchmark is now non-trivial

The benchmark is no longer an internal oracle self-check. Predictor scores are evaluated against an exact reveal oracle.

This makes the benchmark a real metabolite-prioritization task.

3. The commutator-based masked-state surrogate is the strongest current predictor

The current best-performing non-oracle prioritization rule is: gnc_commutator

It is built from a masked-state operator-commutator score and is intended to detect which hidden metabolite would most disrupt the current pathway operator when revealed.

4. Performance remains pathway-dependent

Some pathways are structurally easy, some are hard, and predictor performance is not uniform across cohorts. This is a feature of the problem, not a bug.

⸻

Important implementation notes

ST000356 parsing

ST000356 is stored in a non-standard transposed / metadata-row format.
This is handled explicitly in compute_pathway_features.py and run_real_multi_pathway_benchmark.py.

FGW numerical stability

FGW runs can be numerically sensitive. The scripts use entropic regularization and defensive handling of invalid solves.

Structure matrices

When explicit pathway edges are unavailable, the structure matrix falls back to a correlation-distance construction.

Small pathways

Pathways with very small node counts are harder for:
	•	projection-based preprocessing
	•	graded dropout severity
	•	rank-based evaluation

This should be kept in mind when interpreting results.

⸻

Main output files

Stability benchmark
results/results/jl_stability_benchmark.csv


Contains pathway-level stability summaries for each operator and perturbation regime.

Real multi-pathway benchmark
results/results/real_multi_pathway_results.csv

Contains trial-level oracle benchmark results across datasets, pathways, mask rates, and predictor methods.



⸻

Figures

The figures/ directory contains manuscript figures summarizing:
	•	revised pipeline overview
	•	operator stability heatmap
	•	JL vs random projection validation
	•	global hard-subset benchmark
	•	per-dataset benchmark comparison
	•	pathway heterogeneity
	•	commutator mechanism schematic

These figures reflect the current benchmark state and should be interpreted together with the CSV outputs in results/results/.

⸻

What this repository is and is not

This repository is
	•	a metabolomic pathway identifiability benchmark
	•	a stability-gated operator-selection framework
	•	a non-trivial oracle benchmark for metabolite prioritization

This repository is not
	•	a generic enrichment toolkit
	•	a causal pathway inference engine
	•	a transport phenotyping pipeline
	•	a universal claim that JL always improves alignment

⸻

Citation

If you use this repository, cite the associated manuscript or preprint once finalized.

A placeholder citation can be added here once the manuscript title and venue are fixed.

⸻

Contact / maintenance

This repository is under active methodological revision.
If figures, scripts, and outputs disagree, the CSV outputs in results/results/ and the scripts in results/scripts/ should be treated as the source of truth.

Two important notes before you paste this in:

The public repo page I inspected still shows the older README and older repository-structure description, so this replacement is intentionally bringing the README up to the current pipeline rather than preserving the outdated text.  [oai_citation:1‡GitHub](https://github.com/Anas-Enoch/Pathway_Identifiability)

Also, the repo root on GitHub shows `codes/`, `figures/`, and `results/` as the main directories, so the README above is aligned to that top-level structure.  [oai_citation:2‡GitHub](https://github.com/Anas-Enoch/Pathway_Identifiability)

If you want, I can turn this into a **cleaner journal-style README** with a shorter front page and a separate “Reproducibility” section.
