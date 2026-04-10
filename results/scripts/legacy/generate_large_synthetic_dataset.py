import numpy as np
import pandas as pd

np.random.seed(42)

# ---- PARAMETERS ----
n_samples = 40
n_metabolites = 60
n_pathways = 30
metabs_per_pathway = 6

# ---- Generate metabolite names ----
metabolites = [f"M{i}" for i in range(n_metabolites)]

# ---- Generate sample matrix ----
conditions = ["tumor"] * 20 + ["normal"] * 20
X = pd.DataFrame(
    np.random.normal(0, 1, (n_samples, n_metabolites)),
    columns=metabolites
)
X.insert(0, "condition", conditions)
X.insert(0, "sample_id", [f"S{i}" for i in range(n_samples)])

X.to_csv("x.csv", index=False)

# ---- Generate pathways ----
rows = []
for k in range(n_pathways):
    chosen = np.random.choice(metabolites, metabs_per_pathway, replace=False)
    for m in chosen:
        rows.append({"pathway": f"Pathway_{k}", "metabolite": m})

pd.DataFrame(rows).to_csv("pathways.csv", index=False)

# ---- Generate Uk scores ----
Uk_df = pd.DataFrame({
    "pathway": [f"Pathway_{k}" for k in range(n_pathways)],
    "Uk": np.random.uniform(0, 1, n_pathways)
})
Uk_df.to_csv("Uk_scores.csv", index=False)

print("Synthetic dataset generated.")
