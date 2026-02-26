import numpy as np
import pandas as pd
from scipy.stats import ttest_ind, kendalltau

# -------------------------
# CONFIG
# -------------------------

MASK_RATE = 0.3
R = 50
MIN_MET_PER_PATHWAY = 3
EPS = 1e-8
np.random.seed(0)

# -------------------------
# LOAD DATA
# -------------------------

X = pd.read_csv("X.csv")
pathways = pd.read_csv("pathways.csv")
Uk_df = pd.read_csv("Uk_scores.csv")

met_cols = [c for c in X.columns if c not in ["sample_id", "condition"]]

# Baseline ranking (no masking)
def compute_mean_activity(df):
    scores = {}
    for k, group in pathways.groupby("pathway"):
        mets = list(set(group["metabolite"]) & set(df.columns))
        if len(mets) < MIN_MET_PER_PATHWAY:
            continue

        tumor = df[df.condition=="tumor"][mets].mean().mean()
        normal = df[df.condition=="normal"][mets].mean().mean()
        scores[k] = tumor - normal
    return scores

def compute_enrichment(df):
    scores = {}
    for k, group in pathways.groupby("pathway"):
        mets = list(set(group["metabolite"]) & set(df.columns))
        if len(mets) < MIN_MET_PER_PATHWAY:
            continue

        d_stats = []
        for m in mets:
            tumor = df[df.condition=="tumor"][m].dropna()
            normal = df[df.condition=="normal"][m].dropna()
            if len(tumor)>2 and len(normal)>2:
                t,_ = ttest_ind(tumor, normal, equal_var=False)
                d_stats.append(t)
        if len(d_stats)>0:
            scores[k] = np.mean(d_stats)
    return scores

baseline_mean = compute_mean_activity(X)
baseline_enrich = compute_enrichment(X)

baseline_rank_mean = pd.Series(baseline_mean).abs().sort_values(ascending=False)
baseline_rank_enrich = pd.Series(baseline_enrich).abs().sort_values(ascending=False)

# -------------------------
# STORAGE
# -------------------------

stability_records = []
pathway_rank_drift = {k: [] for k in baseline_rank_mean.index}

# -------------------------
# MASKING LOOP
# -------------------------

for r in range(R):

    X_mask = X.copy()

    for i in range(len(X_mask)):
        mask = np.random.rand(len(met_cols)) < MASK_RATE
        X_mask.loc[i, np.array(met_cols)[mask]] = np.nan

    mean_scores = compute_mean_activity(X_mask)
    enrich_scores = compute_enrichment(X_mask)

    rank_mean = pd.Series(mean_scores).abs().sort_values(ascending=False)
    rank_enrich = pd.Series(enrich_scores).abs().sort_values(ascending=False)

    # Kendall stability
    tau_mean,_ = kendalltau(
        baseline_rank_mean.index[:len(rank_mean)],
        rank_mean.index
    )

    tau_enrich,_ = kendalltau(
        baseline_rank_enrich.index[:len(rank_enrich)],
        rank_enrich.index
    )

    stability_records.append({"method":"Mean", "tau":tau_mean})
    stability_records.append({"method":"Enrichment", "tau":tau_enrich})

    # Per-pathway rank drift
    for k in pathway_rank_drift.keys():
        if k in rank_mean.index:
            drift = abs(
                baseline_rank_mean.index.get_loc(k) -
                rank_mean.index.get_loc(k)
            )
            pathway_rank_drift[k].append(drift)

# -------------------------
# SAVE STABILITY CSV
# -------------------------

stability_df = pd.DataFrame(stability_records)
stability_df.to_csv("ranking_stability_by_method.csv", index=False)

# -------------------------
# Compute pathway instability
# -------------------------

instab_list = []
for k, drifts in pathway_rank_drift.items():
    if len(drifts)>0:
        instab_list.append({
            "pathway": k,
            "Instab": np.mean(drifts)
        })

instab_df = pd.DataFrame(instab_list)

# Merge Uk
final_df = instab_df.merge(Uk_df, on="pathway", how="inner")
final_df.to_csv("pathway_instability_vs_Uk.csv", index=False)

print("CSV files generated successfully.")
