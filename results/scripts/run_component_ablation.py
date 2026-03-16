import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr, kendalltau

# =========================================================
# CONFIG
# =========================================================

COMPONENTS_FILE = "results/Data/Uk_components.csv"
INSTABILITY_FILE = "results/Data/pathway_instability_vs_Uk.csv"

OUT_SUMMARY = "results/Results/ablation_summary.csv"
OUT_SCORES = "results/Results/ablation_scores_by_pathway.csv"

# Default composite weights used in manuscript
ALPHA = 1.0   # H_transport
BETA = 1.0    # alignment_instability
GAMMA = 1.0   # SGI

EPS = 1e-12


# =========================================================
# HELPERS
# =========================================================

def zscore(x: pd.Series) -> pd.Series:
    std = x.std(ddof=0)
    if std < EPS:
        return pd.Series(np.zeros(len(x)), index=x.index)
    return (x - x.mean()) / std


def safe_corr(x: pd.Series, y: pd.Series, method: str = "spearman"):
    x = pd.to_numeric(x, errors="coerce")
    y = pd.to_numeric(y, errors="coerce")
    mask = x.notna() & y.notna()
    x = x[mask]
    y = y[mask]

    if len(x) < 3:
        return np.nan, np.nan

    if method == "spearman":
        r, p = spearmanr(x, y)
    elif method == "pearson":
        r, p = pearsonr(x, y)
    elif method == "kendall":
        r, p = kendalltau(x, y)
    else:
        raise ValueError(f"Unknown method: {method}")

    return r, p


def safe_topk_overlap(score_a: pd.Series, score_b: pd.Series, k: int = 5) -> float:
    a = set(score_a.sort_values(ascending=False).head(k).index)
    b = set(score_b.sort_values(ascending=False).head(k).index)
    denom = len(a | b)
    if denom == 0:
        return np.nan
    return len(a & b) / denom


def normalize_weights(a, b, c):
    s = a + b + c
    if s < EPS:
        return 0.0, 0.0, 0.0
    return a / s, b / s, c / s


# =========================================================
# LOAD DATA
# =========================================================

components = pd.read_csv(COMPONENTS_FILE)
instab = pd.read_csv(INSTABILITY_FILE)

required_component_cols = {"pathway", "H_transport", "SGI", "alignment_instability"}
missing = required_component_cols - set(components.columns)
if missing:
    raise ValueError(
        f"{COMPONENTS_FILE} is missing required columns: {missing}\n"
        f"Expected at least: {sorted(required_component_cols)}"
    )

if "Instab" not in instab.columns:
    raise ValueError(
        f"{INSTABILITY_FILE} must contain an 'Instab' column.\n"
        f"Found columns: {list(instab.columns)}"
    )

# Keep only needed columns
components = components[["pathway", "H_transport", "SGI", "alignment_instability"]].copy()
instab = instab[["pathway", "Instab"]].copy()

# Merge
df = components.merge(instab, on="pathway", how="inner")

if df.empty:
    raise ValueError("Merged dataframe is empty. Check that pathway names match across files.")

# Numeric cast
for col in ["H_transport", "SGI", "alignment_instability", "Instab"]:
    df[col] = pd.to_numeric(df[col], errors="coerce")

df = df.dropna(subset=["H_transport", "SGI", "alignment_instability", "Instab"]).copy()

if df.empty:
    raise ValueError("No valid numeric rows remain after cleaning.")

# =========================================================
# STANDARDIZE COMPONENTS
# =========================================================
# Important:
# If the components are on different scales, raw weighted sums are misleading.
# Standardizing is the safest default for ablation comparison.

df["H_z"] = zscore(df["H_transport"])
df["SGI_z"] = zscore(df["SGI"])
df["INST_z"] = zscore(df["alignment_instability"])

a, b, c = normalize_weights(ALPHA, BETA, GAMMA)

# =========================================================
# BUILD SCORES
# =========================================================

# Single-component
df["score_entropy_only"] = df["H_z"]
df["score_sgi_only"] = df["SGI_z"]
df["score_instability_only"] = df["INST_z"]

# Full composite
df["score_full"] = a * df["H_z"] + b * df["INST_z"] + c * df["SGI_z"]

# Leave-one-out
ab, bb, cb = normalize_weights(0.0, BETA, GAMMA)
df["score_without_entropy"] = bb * df["INST_z"] + cb * df["SGI_z"]

ab, bb, cb = normalize_weights(ALPHA, BETA, 0.0)
df["score_without_sgi"] = ab * df["H_z"] + bb * df["INST_z"]

ab, bb, cb = normalize_weights(ALPHA, 0.0, GAMMA)
df["score_without_instability"] = ab * df["H_z"] + cb * df["SGI_z"]

# Optional raw composite if you want to inspect it
df["score_full_raw"] = (
    ALPHA * df["H_transport"] +
    BETA * df["alignment_instability"] +
    GAMMA * df["SGI"]
)

# =========================================================
# SUMMARY METRICS
# =========================================================

score_cols = [
    "score_entropy_only",
    "score_sgi_only",
    "score_instability_only",
    "score_full",
    "score_without_entropy",
    "score_without_sgi",
    "score_without_instability",
]

summary_rows = []

full_score = df.set_index("pathway")["score_full"]

for col in score_cols:
    score = df.set_index("pathway")[col]
    inst = df.set_index("pathway")["Instab"]

    spearman_r, spearman_p = safe_corr(score, inst, method="spearman")
    pearson_r, pearson_p = safe_corr(score, inst, method="pearson")
    kendall_r, kendall_p = safe_corr(score, inst, method="kendall")

    # Agreement with full composite ranking
    if col != "score_full":
        rank_tau, rank_tau_p = safe_corr(score, full_score, method="kendall")
        top5_overlap = safe_topk_overlap(score, full_score, k=5)
        top10_overlap = safe_topk_overlap(score, full_score, k=10)
    else:
        rank_tau, rank_tau_p = 1.0, 0.0
        top5_overlap = 1.0
        top10_overlap = 1.0

    summary_rows.append({
        "score": col,
        "n_pathways": len(df),
        "spearman_r_with_instability": spearman_r,
        "spearman_p": spearman_p,
        "pearson_r_with_instability": pearson_r,
        "pearson_p": pearson_p,
        "kendall_r_with_instability": kendall_r,
        "kendall_p": kendall_p,
        "kendall_tau_vs_full_composite": rank_tau,
        "kendall_tau_vs_full_composite_p": rank_tau_p,
        "top5_jaccard_vs_full": top5_overlap,
        "top10_jaccard_vs_full": top10_overlap,
        "mean_score": float(score.mean()),
        "std_score": float(score.std(ddof=0)),
    })

summary_df = pd.DataFrame(summary_rows)

# Sort by strongest Spearman magnitude
summary_df["abs_spearman"] = summary_df["spearman_r_with_instability"].abs()
summary_df = summary_df.sort_values(
    by=["abs_spearman", "score"],
    ascending=[False, True]
).drop(columns=["abs_spearman"])

# =========================================================
# SAVE OUTPUTS
# =========================================================

# Save pathway-level score table
score_export_cols = [
    "pathway",
    "H_transport",
    "SGI",
    "alignment_instability",
    "Instab",
    "score_entropy_only",
    "score_sgi_only",
    "score_instability_only",
    "score_full",
    "score_without_entropy",
    "score_without_sgi",
    "score_without_instability",
]

df[score_export_cols].to_csv(OUT_SCORES, index=False)
summary_df.to_csv(OUT_SUMMARY, index=False)

# =========================================================
# PRINT HUMAN-READABLE SUMMARY
# =========================================================

print("\n=== Component Ablation Summary ===\n")
print(summary_df.to_string(index=False))

best_row = summary_df.iloc[0]
print("\nBest correlation with pathway instability:")
print(f"  Score:     {best_row['score']}")
print(f"  Spearman:  {best_row['spearman_r_with_instability']:.4f}")
print(f"  p-value:   {best_row['spearman_p']:.3e}")

print(f"\nSaved:")
print(f"  - {OUT_SUMMARY}")
print(f"  - {OUT_SCORES}")
