#!/usr/bin/env python3
"""
compute_Uk_real.py
──────────────────
Consume cross-condition alignment .npz files and compute:

    U_k = α · H(T) + β · solver_instability + γ · SGI

All three components are computed in the SAME feature space:
    - H(T) comes from the transport plan already produced by
      compute_fgw_alignment.py (preprocessed features).
    - solver_instability re-solves FGW N_INIT times on the SAME
      preprocessed features with tiny cost-matrix jitter, measuring
      CV of the FGW distance.  This quantifies how sensitive the
      OT solution is to numerical perturbation in the projected
      geometry — NOT data-level robustness (which is tested separately
      in run_jl_stability_benchmark.py).
    - SGI is topology-only (independent of features).

Output
    results/data/Uk_components.csv
    Columns: dataset, pathway, method,
             H_transport, solver_instability, SGI, Uk
"""

from __future__ import annotations

from pathlib import Path

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
REPO_ROOT  = Path(__file__).resolve().parents[2]
ALIGN_DIR  = REPO_ROOT / "results" / "results" / "jl_alignments"
FEAT_DIR   = REPO_ROOT / "results" / "results" / "jl_features"
OUT_DIR    = REPO_ROOT / "results" / "data"
OUT_FILE   = OUT_DIR / "Uk_components.csv"

ALPHA = 1.0
BETA  = 1.0
GAMMA = 1.0
EPS   = 1e-12

N_INIT       = 10
SINKHORN_REG = 0.05
ALPHA_FGW    = 0.5
JL_DIM       = 8
RANDOM_STATE = 42


# ──────────────────────────────────────────────
#  Preprocessing  (must match compute_fgw_alignment.py exactly)
# ──────────────────────────────────────────────
def _pad(X_s, X_t):
    d_s, d_t = X_s.shape[1], X_t.shape[1]
    if d_s == d_t:
        return X_s, X_t
    d = max(d_s, d_t)
    if d_s < d:
        X_s = np.hstack([X_s, np.zeros((X_s.shape[0], d - d_s))])
    if d_t < d:
        X_t = np.hstack([X_t, np.zeros((X_t.shape[0], d - d_t))])
    return X_s, X_t


def preprocess(X_s, X_t, method, rng):
    """Apply the SAME preprocessing used by compute_fgw_alignment.py."""
    if method == "none":
        return X_s.copy(), X_t.copy()
    X_s, X_t = _pad(X_s, X_t)
    d_orig = X_s.shape[1]
    n_s = X_s.shape[0]
    X_all = np.vstack([X_s, X_t])

    if method == "jl":
        d = min(JL_DIM, d_orig)
        R = rng.standard_normal((d_orig, d)) / np.sqrt(d)
        Z = X_all @ R
        return Z[:n_s], Z[n_s:]
    if method == "randproj":
        d = min(JL_DIM, d_orig)
        R = rng.standard_normal((d_orig, d))
        Z = X_all @ R
        return Z[:n_s], Z[n_s:]
    if method == "pca_fixed":
        d = min(JL_DIM, d_orig, X_all.shape[0])
        pca = PCA(n_components=d, random_state=RANDOM_STATE)
        Z = pca.fit_transform(X_all)
        return Z[:n_s], Z[n_s:]
    if method == "pca_var95":
        max_d = min(d_orig, X_all.shape[0])
        pca = PCA(n_components=min(0.95, max_d), svd_solver="full")
        Z = pca.fit_transform(X_all)
        return Z[:n_s], Z[n_s:]
    if method == "l2":
        return (normalize(X_s, norm="l2", axis=1),
                normalize(X_t, norm="l2", axis=1))
    raise ValueError(f"Unknown method: {method}")


# ──────────────────────────────────────────────
#  Component metrics
# ──────────────────────────────────────────────
def transport_entropy(T: np.ndarray) -> float:
    """H(T) = −Σ T_ij log(T_ij)  for T_ij > 0."""
    t = T.ravel()
    t = t[t > EPS]
    return float(-np.sum(t * np.log(t)))


def solver_instability(
    X_s: np.ndarray, X_t: np.ndarray,
    C_s: np.ndarray, C_t: np.ndarray,
    n_init: int = N_INIT,
) -> float:
    """CV of FGW distance over n_init re-solves with cost-matrix jitter.

    Measures numerical sensitivity of the entropic FGW solver in the
    PREPROCESSED feature space.  This is solver instability, not
    data-level robustness.
    """
    n_s, n_t = X_s.shape[0], X_t.shape[0]
    p = np.ones(n_s) / n_s
    q = np.ones(n_t) / n_t

    M = cdist(X_s, X_t, "sqeuclidean")
    mx = M.max()
    if mx > 0:
        M = M / mx

    def _norm(C):
        mx = C.max()
        return C / mx if mx > 0 else C

    C_sn, C_tn = _norm(C_s), _norm(C_t)
    dists = []

    for seed in range(n_init):
        rng = np.random.default_rng(seed)
        noise = rng.normal(0, 1e-6, size=M.shape)
        try:
            _, log = ot.gromov.entropic_fused_gromov_wasserstein(
                M + noise, C_sn, C_tn, p, q,
                loss_fun="square_loss", epsilon=SINKHORN_REG,
                alpha=ALPHA_FGW, log=True,
            )
            d = log.get("fgw_dist", float(np.sum(
                log.get("T", np.zeros_like(M)) * (M + noise))))
            dists.append(float(d))
        except Exception:
            continue

    arr = np.array(dists)
    if len(arr) < 2:
        return 0.0
    return float(arr.std() / (arr.mean() + EPS))


def spectral_gap_index(C: np.ndarray) -> float:
    """SGI = λ₂ / λ_max of the graph Laplacian.  Topology-only."""
    pos = C[C > 0]
    if len(pos) == 0:
        return 0.0
    sigma = float(np.median(pos))
    A = np.exp(-C ** 2 / (2 * sigma ** 2 + EPS))
    np.fill_diagonal(A, 0.0)
    L = np.diag(A.sum(1)) - A
    ev = np.sort(np.linalg.eigvalsh(L))
    return float(ev[1] / (ev[-1] + EPS)) if len(ev) > 1 else 0.0


# ──────────────────────────────────────────────
#  Main
# ──────────────────────────────────────────────
def main():
    if not ALIGN_DIR.exists():
        print(f"[ERROR] {ALIGN_DIR} not found — "
              f"run compute_fgw_alignment.py first")
        return

    manifest = pd.read_csv(FEAT_DIR / "manifest.csv")

    # Build lookups
    feat_lookup: dict[tuple, Path] = {}
    struct_lookup: dict[tuple, Path] = {}
    for _, row in manifest.iterrows():
        key = (row["dataset"], row["pathway"], row["condition"])
        feat_lookup[key] = FEAT_DIR / row["feature_file"]
        struct_lookup[(row["dataset"], row["pathway"])] = (
            FEAT_DIR / row["structure_file"]
        )

    rows: list[dict] = []

    for npz_path in sorted(ALIGN_DIR.glob("*.npz")):
        data   = np.load(npz_path, allow_pickle=True)
        ds     = str(data["dataset"])
        pw     = str(data["pathway"])
        method = str(data["method"])
        T      = data["transport_plan"]
        cond_s = str(data["condition_source"])
        cond_t = str(data["condition_target"])

        key_s      = (ds, pw, cond_s)
        key_t      = (ds, pw, cond_t)
        struct_key = (ds, pw)

        if key_s not in feat_lookup or key_t not in feat_lookup:
            continue
        if struct_key not in struct_lookup:
            continue

        # Load RAW features
        X_s_raw = pd.read_csv(
            feat_lookup[key_s], index_col=0,
        ).values.astype(np.float64)
        X_t_raw = pd.read_csv(
            feat_lookup[key_t], index_col=0,
        ).values.astype(np.float64)
        C = np.load(struct_lookup[struct_key]).astype(np.float64)

        # ── CRITICAL: preprocess with the SAME method that produced T ──
        rng = np.random.default_rng(RANDOM_STATE)
        X_s_p, X_t_p = preprocess(X_s_raw, X_t_raw, method, rng)

        # ── Component 1: H(T) from the existing transport plan ──
        H = transport_entropy(T)

        # ── Component 2: solver instability in PREPROCESSED space ──
        instab = solver_instability(X_s_p, X_t_p, C, C)

        # ── Component 3: SGI (topology-only) ──
        sgi = spectral_gap_index(C)

        # ── Composite ──
        Uk = ALPHA * H + BETA * instab + GAMMA * sgi

        rows.append(dict(
            dataset=ds, pathway=pw, method=method,
            H_transport=round(H, 6),
            solver_instability=round(instab, 6),
            SGI=round(sgi, 6),
            Uk=round(Uk, 6),
        ))

    out = pd.DataFrame(rows)
    out.to_csv(OUT_FILE, index=False)
    print(f"[done] {len(out)} rows → {OUT_FILE}")


if __name__ == "__main__":
    main()
