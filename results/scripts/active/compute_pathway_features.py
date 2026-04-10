#!/usr/bin/env python3
"""
compute_pathway_features.py
────────────────────────────
Build per-pathway node-feature matrices from processed cohort CSVs
and core_pathway_mapping.csv.

Feature strategy (critical for JL to be non-trivial):
    PRIMARY:  per-sample abundance profile → (n_nodes, n_samples_per_condition)
              This gives d = n_samples ≈ 20–100, the natural JL regime.
    FALLBACK: rich distributional embedding → (n_nodes, 8)
              [μ, σ², std, q25, q50, q75, skewness, CV]
              Used when n_samples < 6  (too few for stable per-sample features).

Condition split:
    Each dataset is split into condition groups (case vs control).
    Features are computed *per condition* so that downstream FGW aligns
    case-pathway vs control-pathway — not self-alignment.

Output
    results/results/jl_features/{dataset}__{pathway}__{condition}.csv
        rows  = metabolites,  cols = feature dimensions
    results/results/jl_features/{dataset}__{pathway}__structure.npy
    results/results/jl_features/manifest.csv
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse.csgraph import shortest_path

# ──────────────────────────────────────────────
#  Config
# ──────────────────────────────────────────────
REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR  = REPO_ROOT / "results" / "data"
OUT_DIR   = REPO_ROOT / "results" / "results" / "jl_features"
OUT_DIR.mkdir(parents=True, exist_ok=True)

DATASETS = {
    "ST000356": DATA_DIR / "processed_metabolite_matrix_ST000356.csv",
    "ST003390": DATA_DIR / "ST003390_processed.csv",
    "ST003506": DATA_DIR / "ST003506_serum_processed.csv",
}
PATHWAY_MAP = DATA_DIR / "core_pathway_mapping.csv"
MIN_MET_PER_PATHWAY = 3
MIN_SAMPLES_FOR_PROFILE = 6     # below this, use distributional fallback

# Columns that are metadata, not metabolite abundances
META_COLS = {
    "sample_id", "sampleid", "sample",
    "condition", "group", "label", "cohort",
    "class", "subject", "subject_id", "sex", "age",
}

# Column names to try for condition/group detection (priority order)
CONDITION_CANDIDATES = [
    "condition", "group", "label",
    "phenotype", "class", "diagnosis", "status", "cohort",
]


# ──────────────────────────────────────────────
#  Helpers
# ──────────────────────────────────────────────
def detect_metabolite_columns(df: pd.DataFrame) -> list[str]:
    return [c for c in df.columns if c.lower().strip() not in META_COLS]


def detect_cond_col(df: pd.DataFrame) -> str | None:
    """Find the condition/group column (case-insensitive)."""
    lower_map = {c.lower().strip(): c for c in df.columns}
    for cand in CONDITION_CANDIDATES:
        if cand in lower_map:
            return lower_map[cand]
    return None


def safe_name(s: str) -> str:
    return s.replace("/", "_").replace(" ", "_").replace(":", "_")


# ──────────────────────────────────────────────
#  Pathway mapping
# ──────────────────────────────────────────────
def load_pathway_mapping(path):
    df = pd.read_csv(path)
    df.columns = [c.strip().lower() for c in df.columns]

    # Accept both old and current schema
    if "metabolite_name" in df.columns and "pathway_name" in df.columns:
        df = df.rename(columns={
            "metabolite_name": "metabolite",
            "pathway_name": "pathway",
        })

    assert {"metabolite", "pathway"}.issubset(df.columns), (
        f"Expected 'metabolite' and 'pathway'; got {df.columns.tolist()}"
    )

    pw_members = (
        df.groupby("pathway")["metabolite"]
        .apply(lambda x: sorted(set(map(str, x))))
        .to_dict()
    )

    # If no explicit edge information exists, return empty edge dict
    pw_edges = {pw: [] for pw in pw_members}

    return pw_members, pw_edges


# ──────────────────────────────────────────────
#  Structure matrix
# ──────────────────────────────────────────────
def build_structure_matrix(nodes, edges):
    idx = {n: i for i, n in enumerate(nodes)}
    n = len(nodes)
    adj = np.zeros((n, n))
    for u, v in edges:
        if u in idx and v in idx:
            adj[idx[u], idx[v]] = 1.0
            adj[idx[v], idx[u]] = 1.0
    dist = shortest_path(adj, method="D", unweighted=True)
    dist[np.isinf(dist)] = float(n)
    return dist


def correlation_distance_matrix(X):
    """d(i,j) = sqrt(2(1 − ρ_ij)).  X is (n_nodes, n_features)."""
    if X.shape[0] < 2 or X.shape[1] < 2:
        return np.zeros((X.shape[0], X.shape[0]))
    R = np.corrcoef(X)
    R = np.clip(R, -1.0, 1.0)
    return np.sqrt(np.maximum(2.0 * (1.0 - R), 0.0))


# ──────────────────────────────────────────────
#  Node features
# ──────────────────────────────────────────────
def build_node_features_profile(sub):
    """Per-sample profile, z-scored per metabolite.

    Input:  (n_samples, n_nodes)
    Output: (n_nodes, n_samples)
    """
    eps = 1e-12
    X = sub.T.copy()                               # (n_nodes, n_samples)
    mean = np.nanmean(X, axis=1, keepdims=True)
    std  = np.nanstd(X, axis=1, keepdims=True)
    X = (X - mean) / (std + eps)
    X = np.nan_to_num(X, nan=0.0)
    return X


def build_node_features_distributional(sub: np.ndarray) -> np.ndarray:
    """Rich distributional embedding (fallback for small n_samples).

    Input:  (n_samples, n_nodes)
    Output: (n_nodes, 8)
    """
    eps = 1e-12
    mu  = np.nanmean(sub, axis=0)
    var = np.nanvar(sub, axis=0)
    std = np.sqrt(var + eps)
    q25 = np.nanpercentile(sub, 25, axis=0)
    q50 = np.nanpercentile(sub, 50, axis=0)
    q75 = np.nanpercentile(sub, 75, axis=0)
    skew = np.nanmean(((sub - mu) / (std + eps)) ** 3, axis=0)
    cv   = std / (np.abs(mu) + eps)
    result = np.column_stack([mu, var, std, q25, q50, q75, skew, cv])
    result = np.nan_to_num(result, nan=0.0)
    return result


def build_node_features(sub: np.ndarray) -> np.ndarray:
    """Dispatch: use per-sample profile if enough samples, else distributional."""
    n_samples = sub.shape[0]
    if n_samples >= MIN_SAMPLES_FOR_PROFILE:
        return build_node_features_profile(sub)
    else:
        return build_node_features_distributional(sub)


# ──────────────────────────────────────────────
#  ST000356 special handling
# ──────────────────────────────────────────────
def fix_st000356(df: pd.DataFrame) -> pd.DataFrame:
    """Handle ST000356's transposed format with embedded condition row.

    The CSV has:
        - Column 0: some ID column (or index)
        - Column 1: Metabolite_name
        - Columns 2+: sample columns
        - Row 0 (data): condition labels embedded in the sample columns

    We need to:
        1. Extract condition labels from row 0
        2. Drop that metadata row
        3. Transpose so samples become rows
        4. Attach condition column
    """
    if "Metabolite_name" not in df.columns:
        return df

    # Identify sample columns (everything after Metabolite_name)
    met_col_idx = df.columns.get_loc("Metabolite_name")
    sample_cols = df.columns[met_col_idx + 1:]

    # Row 0 contains condition labels for each sample column
    cond_row = df.iloc[0][sample_cols]

    def extract_label(x):
        x = str(x).strip().lower()

        if any(k in x for k in ["control", "healthy", "normal"]):
            return "control"

        if any(k in x for k in ["cancer", "tumor", "tumour", "breast cancer"]):
            return "cancer"

        if any(k in x for k in ["qc", "blank", "pooled", "unknown", "na", "nan"]):
            return "unknown"

        return "unknown"

    conditions = cond_row.apply(extract_label)
    bad = cond_row[conditions == "unknown"]
    if len(bad) > 0:
        print("\n[ST000356] Unknown condition labels detected:")
        print(bad.to_string())
        
    # Drop the condition-label row, keep only metabolite data rows
    df_data = df.iloc[1:].copy()

    # Set metabolite names as index, keep only sample columns
    df_data = df_data.set_index("Metabolite_name")[sample_cols]

    # Force numeric (some cells may still be strings)
    df_data = df_data.apply(pd.to_numeric, errors="coerce")

    # Transpose: now rows = samples, columns = metabolites
    df_t = df_data.T.reset_index(drop=True)

    # Attach condition — lengths now guaranteed to match because
    # both conditions and df_t have len(sample_cols) rows
    df_t["condition"] = conditions.values
    df_t = df_t[df_t["condition"] != "unknown"].reset_index(drop=True)

    return df_t


# ──────────────────────────────────────────────
#  Main
# ──────────────────────────────────────────────
def main():
    pw_members, pw_edges = load_pathway_mapping(PATHWAY_MAP)
    manifest: list[dict] = []

    for ds_name, ds_path in DATASETS.items():
        if not ds_path.exists():
            print(f"[WARN] {ds_path} not found — skipping {ds_name}")
            continue

        df = pd.read_csv(ds_path)

        # ST000356 has a transposed format with embedded condition row
        if ds_name == "ST000356":
            df = fix_st000356(df)

        met_cols = detect_metabolite_columns(df)
        met_set  = set(met_cols)
        cond_col = detect_cond_col(df)          # ← was detect_condition_column (bug)

        if cond_col is None:
            print(f"[WARN] {ds_name}: no condition column found — "
                  f"treating entire cohort as one group")
            conditions = {"all": df}
        else:
            groups = df[cond_col].unique()
            conditions = {str(g): df[df[cond_col] == g] for g in groups}
            print(f"[{ds_name}] condition column='{cond_col}', "
                  f"groups={list(conditions.keys())}")

        print(f"[{ds_name}] {len(met_cols)} metabolite columns, "
              f"{len(df)} total samples")

        for pw_name, members in pw_members.items():
            available = [m for m in members if m in met_set]
            if len(available) < MIN_MET_PER_PATHWAY:
                continue

            pw_safe = safe_name(pw_name)

            # Structure matrix (shared across conditions)
            sub_all = df[available].values.astype(np.float64)
            if pw_name in pw_edges and pw_edges[pw_name]:
                C = build_structure_matrix(available, pw_edges[pw_name])
            else:
                C = correlation_distance_matrix(sub_all.T)

            struct_file = f"{ds_name}__{pw_safe}__structure.npy"
            np.save(OUT_DIR / struct_file, C)

            # Features per condition
            for cond_name, cond_df in conditions.items():
                sub = cond_df[available].values.astype(np.float64)
                if sub.shape[0] < 2:
                    continue

                features = build_node_features(sub)
                n_feat = features.shape[1]

                cond_safe = safe_name(cond_name)
                feat_file = f"{ds_name}__{pw_safe}__{cond_safe}.csv"

                feat_df = pd.DataFrame(
                    features,
                    index=available,
                    columns=[f"f{i}" for i in range(n_feat)],
                )
                feat_df.index.name = "metabolite"
                feat_df.to_csv(OUT_DIR / feat_file)

                manifest.append(dict(
                    dataset=ds_name,
                    pathway=pw_name,
                    condition=cond_name,
                    n_nodes=len(available),
                    n_features=n_feat,
                    n_samples=sub.shape[0],
                    feature_file=feat_file,
                    structure_file=struct_file,
                    metabolites=";".join(available),
                ))

        n_entries = sum(1 for r in manifest if r["dataset"] == ds_name)
        print(f"  → {n_entries} (pathway × condition) entries")

    pd.DataFrame(manifest).to_csv(OUT_DIR / "manifest.csv", index=False)
    print(f"\n[done] {len(manifest)} entries → {OUT_DIR}/")


if __name__ == "__main__":
    main()
