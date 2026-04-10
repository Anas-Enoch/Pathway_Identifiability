#!/usr/bin/env python3
"""
compute_fgw_alignment.py
────────────────────────
Read per-condition pathway feature CSVs from jl_features/, apply
preprocessing, and compute CROSS-CONDITION FGW alignment.

Critical fix: aligns case vs control (not self-alignment).
This is what FGW is for — soft correspondence between pathway
states under different biological conditions.

Output
    results/results/jl_alignments/{dataset}__{pathway}__{method}.npz
    Keys: distance, transport_plan, method, dataset, pathway,
          condition_source, condition_target

Methods:
    none       – raw features (d = n_samples or d = 8)
    jl         – JL-scaled Gaussian projection  (1/√d scaling)
    randproj   – unscaled Gaussian projection   (no 1/√d)
    pca_fixed  – PCA to same d as JL target
    pca_var95  – PCA retaining 95% variance
    l2         – L2 row-normalisation
"""

from __future__ import annotations

from pathlib import Path
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from sklearn.decomposition import PCA
from sklearn.preprocessing import normalize

import ot
import ot.gromov

# ──────────────────────────────────────────────
#  Config
# ──────────────────────────────────────────────
REPO_ROOT    = Path(__file__).resolve().parents[2]
FEATURE_DIR  = REPO_ROOT / "results" / "results" / "jl_features"
OUT_DIR      = REPO_ROOT / "results" / "results" / "jl_alignments"
OUT_DIR.mkdir(parents=True, exist_ok=True)

METHODS      = ["none", "jl", "randproj", "pca_fixed", "pca_var95", "l2"]
JL_DIM       = 8
RANDOM_STATE = 42
ALPHA        = 0.5
SINKHORN_REG = 0.05
EPS          = 1e-12


# ──────────────────────────────────────────────
#  Preprocessing
# ──────────────────────────────────────────────
def preprocess(
    X_s: np.ndarray,
    X_t: np.ndarray,
    method: str,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    """Apply shared preprocessing to source and target features.

    IMPORTANT: projection is fitted on the concatenation of both sides
    so that source and target live in the same projected space.
    """
    if method == "none":
        return X_s.copy(), X_t.copy()

    d_orig = X_s.shape[1]
    n_s    = X_s.shape[0]
    X_all  = np.vstack([X_s, X_t])

    if method == "jl":
        d = min(JL_DIM, d_orig)
        R = rng.standard_normal((d_orig, d)) / np.sqrt(d)
        Z = X_all @ R
        return Z[:n_s], Z[n_s:]

    if method == "randproj":
        d = min(JL_DIM, d_orig)
        R = rng.standard_normal((d_orig, d))           # no 1/√d
        Z = X_all @ R
        return Z[:n_s], Z[n_s:]

    if method == "pca_fixed":
        d = min(JL_DIM, d_orig, X_all.shape[0])
        pca = PCA(n_components=d, random_state=RANDOM_STATE)
        Z = pca.fit_transform(X_all)
        return Z[:n_s], Z[n_s:]

    if method == "pca_var95":
        max_d = min(d_orig, X_all.shape[0])
        pca = PCA(n_components=0.95, svd_solver="full")
        Z = pca.fit_transform(X_all)
        return Z[:n_s], Z[n_s:]

    if method == "l2":
        return (normalize(X_s, norm="l2", axis=1),
                normalize(X_t, norm="l2", axis=1))

    raise ValueError(f"Unknown method: {method}")


# ──────────────────────────────────────────────
#  FGW
# ──────────────────────────────────────────────
def compute_fgw(
    X_s: np.ndarray, X_t: np.ndarray,
    C_s: np.ndarray, C_t: np.ndarray,
    alpha: float = ALPHA, reg: float = SINKHORN_REG,
) -> tuple[float, np.ndarray]:
    n_s, n_t = X_s.shape[0], X_t.shape[0]
    p = np.ones(n_s) / n_s
    q = np.ones(n_t) / n_t

    M = cdist(X_s, X_t, metric="sqeuclidean")
    M_max = M.max()
    if M_max > 0:
        M /= M_max

    def _norm(C):
        mx = C.max()
        return C / mx if mx > 0 else C

    T, log = ot.gromov.entropic_fused_gromov_wasserstein(
        M, _norm(C_s), _norm(C_t), p, q,
        loss_fun="square_loss", epsilon=reg,
        alpha=alpha, log=True,
    )
    dist = log.get("fgw_dist", float(np.sum(T * M)))
    return float(dist), T


# ──────────────────────────────────────────────
#  Main
# ──────────────────────────────────────────────
def main():
    manifest_path = FEATURE_DIR / "manifest.csv"
    if not manifest_path.exists():
        print(f"[ERROR] {manifest_path} not found — "
              f"run compute_pathway_features.py first")
        return

    manifest = pd.read_csv(manifest_path)
    rng = np.random.default_rng(RANDOM_STATE)
    n_total = 0

    # Group by (dataset, pathway) → list of conditions
    grouped = manifest.groupby(["dataset", "pathway"])

    for (ds, pw), grp in grouped:
        conditions = grp["condition"].tolist()
        if len(conditions) < 2:
            print(f"  [{ds}] {pw}: only 1 condition ({conditions}) — "
                  f"skipping (need case vs control)")
            continue

        # Load structure (shared)
        struct_file = grp.iloc[0]["structure_file"]
        C = np.load(FEATURE_DIR / struct_file).astype(np.float64)

        # Load features per condition
        feat_by_cond: dict[str, np.ndarray] = {}
        for _, row in grp.iterrows():
            feat_df = pd.read_csv(FEATURE_DIR / row["feature_file"], index_col=0)
            feat_by_cond[row["condition"]] = feat_df.values.astype(np.float64)

        # Align each pair of conditions
        for cond_s, cond_t in combinations(conditions, 2):
            X_s_raw = feat_by_cond[cond_s]
            X_t_raw = feat_by_cond[cond_t]

            # Ensure same dimensionality (pad shorter with zeros if needed —
            # this handles unequal sample sizes between conditions)
            d_s, d_t = X_s_raw.shape[1], X_t_raw.shape[1]
            if d_s != d_t:
                d_max = max(d_s, d_t)
                if d_s < d_max:
                    X_s_raw = np.hstack([X_s_raw, np.zeros((X_s_raw.shape[0], d_max - d_s))])
                if d_t < d_max:
                    X_t_raw = np.hstack([X_t_raw, np.zeros((X_t_raw.shape[0], d_max - d_t))])

            pw_safe = pw.replace("/", "_").replace(" ", "_").replace(":", "_")

            for method in METHODS:
                X_s_p, X_t_p = preprocess(X_s_raw, X_t_raw, method, rng)
                dist, T = compute_fgw(X_s_p, X_t_p, C, C)

                out_name = f"{ds}__{pw_safe}__{method}.npz"
                np.savez(
                    OUT_DIR / out_name,
                    distance=np.array(dist),
                    transport_plan=T,
                    method=method,
                    dataset=ds,
                    pathway=pw,
                    condition_source=cond_s,
                    condition_target=cond_t,
                )
                n_total += 1

            print(f"  [{ds}] {pw}: {cond_s} vs {cond_t} → {len(METHODS)} methods")

    print(f"\n[done] {n_total} alignment files → {OUT_DIR}/")


if __name__ == "__main__":
    main()
