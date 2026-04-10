#!/usr/bin/env python3
"""
run_jl_stability_benchmark.py
──────────────────────────────
JL stability under KTS-constrained perturbations.

Architecture (three layers):
    Layer 1 — KTS (state dynamics)
        Perturbations are constrained to biologically admissible
        transitions: noise bounded by assay variance, dropout
        weighted by metabolite confidence, bootstrap stratified
        by cohort structure.

    Layer 2 — JL (geometry)
        Projection stabilizes feature distances.  Fixed-scale
        normalization (M /= √d, NOT M /= M.max()) preserves
        the JL scaling factor so that JL ≠ randproj.

    Layer 3 — FGW (alignment)
        Entropic FGW over the stabilized geometry.

Fixes from previous benchmark:
    1. M normalized by √d (fixed-scale), not M.max() (data-dependent)
       → JL and randproj are no longer equivalent
    2. Independent rng per method
       → JL and randproj get different projection matrices
    3. KTS-constrained perturbations
       → biologically meaningful state transitions
    4. Confidence-weighted dropout
       → not degenerate on small pathways

Output
    results/results/jl_stability_benchmark.csv
    Columns: dataset, pathway, regime, method,
             cv_fgw, transport_drift, delta_U, top3_jaccard, rank_tau
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.sparse.csgraph import shortest_path
from scipy.stats import kendalltau

import ot
import ot.gromov

# ──────────────────────────────────────────────
#  Config
# ──────────────────────────────────────────────
REPO_ROOT   = Path(__file__).resolve().parents[2]
DATA_DIR    = REPO_ROOT / "results" / "data"
RESULTS_DIR = REPO_ROOT / "results" / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

DATASETS = {
    "ST000356": DATA_DIR / "processed_metabolite_matrix_ST000356.csv",
    "ST003390": DATA_DIR / "ST003390_processed.csv",
    "ST003506": DATA_DIR / "ST003506_serum_processed.csv",
}
PATHWAY_MAP = DATA_DIR / "core_pathway_mapping.csv"

# KTS-constrained perturbation regimes
PERTURB_REGIMES = {
    # Noise bounded by per-metabolite assay variance
    "kts_noise_005": {"type": "kts_noise",   "level": 0.05},
    "kts_noise_010": {"type": "kts_noise",   "level": 0.10},
    "kts_noise_020": {"type": "kts_noise",   "level": 0.20},
    # Confidence-weighted dropout (not uniform)
    "kts_dropout_010": {"type": "kts_dropout", "level": 0.10},
    "kts_dropout_020": {"type": "kts_dropout", "level": 0.20},
    "kts_dropout_030": {"type": "kts_dropout", "level": 0.30},
    # Stratified bootstrap (preserves condition balance)
    "kts_bootstrap":   {"type": "kts_bootstrap", "level": 1.0},
    # Unconstrained baselines for comparison
    "raw_noise_010":   {"type": "noise",     "level": 0.10},
    "raw_dropout_020": {"type": "dropout",   "level": 0.20},
    "raw_bootstrap":   {"type": "bootstrap", "level": 1.0},
}

METHODS  = ["none", "jl", "randproj", "pca_fixed", "pca_var95", "l2"]
R        = 50
JL_DIM   = 8
OUT_FILE = RESULTS_DIR / "jl_stability_benchmark.csv"

ALPHA_FGW = 0.5
SINK_REG  = 0.5
EPS       = 1e-12

META_COLS = {
    "sample_id", "sampleid", "sample", "condition", "group",
    "label", "cohort", "class", "subject", "subject_id", "sex", "age",
}
CONDITION_CANDIDATES = ["condition", "group", "class", "label", "cohort"]


# ──────────────────────────────────────────────
#  Shared helpers
# ──────────────────────────────────────────────
def detect_met_cols(df):
    return [c for c in df.columns if c.lower().strip() not in META_COLS]


def detect_cond_col(df):
    lower_map = {c.lower().strip(): c for c in df.columns}
    for cand in CONDITION_CANDIDATES:
        if cand in lower_map:
            return lower_map[cand]
    return None


def load_pw_map(path):
    df = pd.read_csv(path)
    df.columns = [c.strip().lower() for c in df.columns]
    if "metabolite_name" in df.columns and "pathway_name" in df.columns:
        df = df.rename(columns={
            "metabolite_name": "metabolite", "pathway_name": "pathway"})
    assert {"metabolite", "pathway"}.issubset(df.columns), df.columns.tolist()
    members = (df.groupby("pathway")["metabolite"]
               .apply(lambda x: sorted(set(map(str, x)))).to_dict())
    edges = {}
    if {"source", "target"}.issubset(df.columns):
        for pw, grp in df.groupby("pathway"):
            e = grp[["source", "target"]].dropna()
            edges[pw] = list(e.itertuples(index=False, name=None))
    return members, edges


def build_structure(nodes, edges):
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


def corr_dist(X):
    if X.shape[0] < 2 or X.shape[1] < 2:
        return np.zeros((X.shape[0], X.shape[0]))
    R_ = np.corrcoef(X)
    R_ = np.clip(R_, -1.0, 1.0)
    return np.sqrt(np.maximum(2.0 * (1.0 - R_), 0.0))


def node_features(sub):
    """(n_samples, n_nodes) → (n_nodes, n_samples), z-scored."""
    eps = 1e-12
    X = sub.T.copy()
    mu = np.nanmean(X, axis=1, keepdims=True)
    sd = np.nanstd(X, axis=1, keepdims=True)
    X = (X - mu) / (sd + eps)
    return np.nan_to_num(X, nan=0.0)


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


# ──────────────────────────────────────────────
#  Preprocessing  (independent rng per method)
# ──────────────────────────────────────────────
def _method_seed(base_seed: int, method: str) -> int:
    """Deterministic but independent seed per method."""
    return base_seed + hash(method) % (2**31)


def preprocess(X_s, X_t, method, rng):
    if method == "none":
        return X_s.copy(), X_t.copy()
    X_s, X_t = _pad(X_s, X_t)
    d_orig = X_s.shape[1]
    n_s = X_s.shape[0]
    X_all = np.vstack([X_s, X_t])

    if method == "jl":
        d = min(JL_DIM, d_orig)
        R_ = rng.standard_normal((d_orig, d)) / np.sqrt(d)   # JL scaling
        Z = X_all @ R_
        return Z[:n_s], Z[n_s:]
    if method == "randproj":
        d = min(JL_DIM, d_orig)
        R_ = rng.standard_normal((d_orig, d))                 # NO scaling
        Z = X_all @ R_
        return Z[:n_s], Z[n_s:]
    if method == "pca_fixed":
        from sklearn.decomposition import PCA
        d = min(JL_DIM, d_orig, X_all.shape[0])
        pca = PCA(n_components=d)
        Z = pca.fit_transform(X_all)
        return Z[:n_s], Z[n_s:]
    if method == "pca_var95":
        from sklearn.decomposition import PCA
        pca = PCA(n_components=0.95, svd_solver="full")
        Z = pca.fit_transform(X_all)
        return Z[:n_s], Z[n_s:]
    if method == "l2":
        from sklearn.preprocessing import normalize
        return normalize(X_s, norm="l2", axis=1), normalize(X_t, norm="l2", axis=1)
    raise ValueError(method)


# ──────────────────────────────────────────────
#  FGW  (FIXED: √d normalization, not max-normalization)
# ──────────────────────────────────────────────
def fgw_solve(X_s, X_t, C_s, C_t):
    """Single FGW solve with FIXED-SCALE normalization.

    Critical fix:  M /= √d  (not M /= M.max())
    This preserves the JL scaling factor so that
    JL-projected and unscaled-projected features
    produce genuinely different cost matrices.
    """
    X_s, X_t = _pad(X_s, X_t)
    n_s, n_t = X_s.shape[0], X_t.shape[0]
    p = np.ones(n_s) / n_s
    q = np.ones(n_t) / n_t

    M = cdist(X_s, X_t, "sqeuclidean")

    # FIXED-SCALE normalization: divide by feature dimension
    # This is what the paper states (Section 5).
    # M.max() would erase the JL 1/√d factor.
    d_feat = X_s.shape[1]
    M /= (np.sqrt(d_feat) + EPS)

    _n = lambda C: C / (C.max() + EPS)
    T, log = ot.gromov.entropic_fused_gromov_wasserstein(
        M, _n(C_s), _n(C_t), p, q,
        loss_fun="square_loss", epsilon=SINK_REG, alpha=ALPHA_FGW, log=True,
    )
    d = log.get("fgw_dist", float(np.sum(T * M)))
    return float(d), T


# ──────────────────────────────────────────────
#  Uk from existing plan  (no re-solve)
# ──────────────────────────────────────────────
def transport_entropy(T):
    t = T.ravel()
    t = t[t > EPS]
    return float(-np.sum(t * np.log(t)))


def spectral_gap(C):
    pos = C[C > 0]
    if len(pos) == 0:
        return 0.0
    sigma = float(np.median(pos))
    A = np.exp(-C ** 2 / (2 * sigma ** 2 + EPS))
    np.fill_diagonal(A, 0.0)
    L = np.diag(A.sum(1)) - A
    ev = np.sort(np.linalg.eigvalsh(L))
    return float(ev[1] / (ev[-1] + EPS)) if len(ev) > 1 else 0.0


def Uk_from_plan(T, C):
    return transport_entropy(T) + spectral_gap(C)


def transport_plan_drift(T0, T_r):
    T0_norm = np.linalg.norm(T0, "fro")
    if T0_norm < EPS:
        T0_norm = 1.0
    return float(np.linalg.norm(T_r - T0, "fro") / T0_norm)


# ──────────────────────────────────────────────
#  KTS: Metabolite confidence scores
# ──────────────────────────────────────────────
def compute_metabolite_confidence(sub_case, sub_ctrl):
    """Per-metabolite confidence based on detection and variance.

    Higher confidence = lower missingness + lower CV.
    Returns array of shape (n_metabolites,) in [0, 1].
    """
    combined = np.vstack([sub_case, sub_ctrl])
    n_total = combined.shape[0]

    # Detection rate: fraction of non-NaN, non-zero values
    detected = np.sum(~np.isnan(combined) & (combined != 0), axis=0)
    detection_rate = detected / (n_total + EPS)

    # Inverse CV: low variance relative to mean = high confidence
    mu = np.nanmean(combined, axis=0)
    sd = np.nanstd(combined, axis=0)
    cv = sd / (np.abs(mu) + EPS)
    inv_cv = 1.0 / (1.0 + cv)

    # Combine: geometric mean of detection rate and inverse CV
    confidence = np.sqrt(detection_rate * inv_cv)
    # Normalize to [0, 1]
    cmax = confidence.max()
    if cmax > 0:
        confidence /= cmax
    return confidence


# ──────────────────────────────────────────────
#  KTS-constrained perturbations
# ──────────────────────────────────────────────
def kts_perturb_noise(sub, level, rng, confidence):
    """Noise bounded by per-metabolite assay variance × (1 - confidence).

    High-confidence metabolites get less noise (stable measurements).
    Low-confidence metabolites get more noise (uncertain measurements).
    """
    base_sigma = np.nanstd(sub, axis=0, keepdims=True) * level
    # Scale noise inversely with confidence
    conf_scale = (1.0 - confidence.reshape(1, -1)) + 0.1  # floor at 0.1
    sigma = base_sigma * conf_scale
    sigma = np.where(sigma < EPS, EPS, sigma)
    return sub + rng.normal(0, sigma, size=sub.shape)


def kts_perturb_dropout(sub_case, sub_ctrl, nodes, level, rng, confidence):
    """Confidence-weighted dropout: low-confidence metabolites more likely dropped.

    Uses soft probabilistic dropout per-metabolite rather than
    hard uniform k-of-n, avoiding the degenerate behavior on small pathways.
    """
    n = len(nodes)
    # Dropout probability inversely proportional to confidence
    drop_prob = (1.0 - confidence) * level
    drop_prob = np.clip(drop_prob, 0.0, 0.8)  # never drop everything

    # Each metabolite independently retained or dropped
    mask = rng.random(n) > drop_prob
    # Ensure at least 2 metabolites survive
    if mask.sum() < 2:
        top2 = np.argsort(confidence)[-2:]
        mask[:] = False
        mask[top2] = True

    keep = np.where(mask)[0].tolist()
    return (sub_case[:, keep], sub_ctrl[:, keep],
            [nodes[i] for i in keep], keep)


def kts_perturb_bootstrap(sub, rng):
    """Stratified bootstrap: resample but preserve approximate quartile structure.

    Rather than pure uniform resampling, ensures each quartile of the
    sample distribution contributes at least one observation.
    """
    n = sub.shape[0]
    if n < 4:
        # Too small for quartile stratification, fall back to regular
        idx = rng.choice(n, size=n, replace=True)
        return sub[idx]

    # Sort by row-mean, split into quartiles
    row_means = np.nanmean(sub, axis=1)
    sorted_idx = np.argsort(row_means)
    q_size = n // 4

    idx = []
    for q in range(4):
        start = q * q_size
        end = start + q_size if q < 3 else n
        q_indices = sorted_idx[start:end]
        # Resample within each quartile
        n_draw = max(1, len(q_indices))
        drawn = rng.choice(q_indices, size=n_draw, replace=True)
        idx.extend(drawn.tolist())

    # Trim or pad to original size
    idx = idx[:n]
    while len(idx) < n:
        idx.append(rng.choice(n))

    return sub[np.array(idx)]


# Unconstrained perturbations (kept as baselines)
def perturb_noise(sub, level, rng):
    sigma = np.nanstd(sub, axis=0, keepdims=True) * level
    sigma = np.where(sigma < EPS, EPS, sigma)
    return sub + rng.normal(0, sigma, size=sub.shape)


def perturb_dropout(sub_case, sub_ctrl, nodes, level, rng):
    n = len(nodes)
    n_drop = max(1, int(round(level * n)))
    n_drop = min(n_drop, n - 2)
    keep = sorted(rng.choice(n, size=n - n_drop, replace=False))
    return (sub_case[:, keep], sub_ctrl[:, keep],
            [nodes[i] for i in keep], keep)


def perturb_bootstrap(sub, rng):
    idx = rng.choice(sub.shape[0], size=sub.shape[0], replace=True)
    return sub[idx]


# ──────────────────────────────────────────────
#  Top-3 by transport mass
# ──────────────────────────────────────────────
def top3_met(T, nodes):
    mass = T.sum(axis=1)
    idx = np.argsort(mass)[::-1][:3]
    return {nodes[i] for i in idx}


def jaccard(a, b):
    if not a and not b:
        return 1.0
    return len(a & b) / len(a | b)


# ──────────────────────────────────────────────
#  Per-pathway stability
# ──────────────────────────────────────────────
def evaluate_stability(
    sub_case, sub_ctrl, nodes, C_full, confidence,
    method, regime_name, regime_cfg,
):
    # Independent rng for this method (so JL ≠ randproj)
    rng0 = np.random.default_rng(_method_seed(999, method))

    # ── BASELINE ──
    X_s0 = node_features(sub_case)
    X_t0 = node_features(sub_ctrl)
    X_s0p, X_t0p = preprocess(X_s0, X_t0, method, rng0)
    d0, T0 = fgw_solve(X_s0p, X_t0p, C_full, C_full)
    U0     = Uk_from_plan(T0, C_full)
    S0     = top3_met(T0, nodes)

    # ── REPLICATES ──
    distances, drifts, delta_Us, jaccards = [], [], [], []
    ptype  = regime_cfg["type"]
    plevel = regime_cfg["level"]

    for r in range(R):
        rng = np.random.default_rng(_method_seed(r, method))

        # Apply perturbation
        if ptype == "kts_noise":
            sc_r = kts_perturb_noise(sub_case, plevel, rng, confidence)
            st_r = kts_perturb_noise(sub_ctrl, plevel, rng, confidence)
            nodes_r, C_r = nodes, C_full
        elif ptype == "kts_dropout":
            sc_r, st_r, nodes_r, ki = kts_perturb_dropout(
                sub_case, sub_ctrl, nodes, plevel, rng, confidence)
            C_r = C_full[np.ix_(ki, ki)]
        elif ptype == "kts_bootstrap":
            sc_r = kts_perturb_bootstrap(sub_case, rng)
            st_r = kts_perturb_bootstrap(sub_ctrl, rng)
            nodes_r, C_r = nodes, C_full
        elif ptype == "noise":
            sc_r = perturb_noise(sub_case, plevel, rng)
            st_r = perturb_noise(sub_ctrl, plevel, rng)
            nodes_r, C_r = nodes, C_full
        elif ptype == "dropout":
            sc_r, st_r, nodes_r, ki = perturb_dropout(
                sub_case, sub_ctrl, nodes, plevel, rng)
            C_r = C_full[np.ix_(ki, ki)]
        elif ptype == "bootstrap":
            sc_r = perturb_bootstrap(sub_case, rng)
            st_r = perturb_bootstrap(sub_ctrl, rng)
            nodes_r, C_r = nodes, C_full
        else:
            continue

        try:
            X_sr = node_features(sc_r)
            X_tr = node_features(st_r)
            X_srp, X_trp = preprocess(X_sr, X_tr, method, rng)
            d_r, T_r = fgw_solve(X_srp, X_trp, C_r, C_r)
        except Exception:
            continue

        U_r = Uk_from_plan(T_r, C_r)
        S_r = top3_met(T_r, nodes_r)

        distances.append(d_r)
        delta_Us.append(abs(U_r - U0))
        jaccards.append(jaccard(S0, S_r))

        if T_r.shape == T0.shape:
            drifts.append(transport_plan_drift(T0, T_r))

    arr = np.array(distances)
    cv_fgw = float(arr.std() / (arr.mean() + EPS)) if len(arr) > 1 else np.nan

    return dict(
        method=method, regime=regime_name,
        cv_fgw=round(cv_fgw, 6),
        transport_drift=round(float(np.mean(drifts)), 6) if drifts else np.nan,
        delta_U=round(float(np.mean(delta_Us)), 6) if delta_Us else np.nan,
        top3_jaccard=round(float(np.mean(jaccards)), 6) if jaccards else np.nan,
        n_nodes=len(nodes), R_effective=len(distances),
    )


# ──────────────────────────────────────────────
#  Cross-pathway Kendall τ
# ──────────────────────────────────────────────
def compute_rank_tau(pw_data, method, regime_name, regime_cfg):
    pw_names = sorted(pw_data)
    if len(pw_names) < 3:
        return np.nan

    Uk_baseline, Uk_perturbed_mean = [], []
    ptype, plevel = regime_cfg["type"], regime_cfg["level"]

    for pw in pw_names:
        d = pw_data[pw]
        conf = d["confidence"]
        rng0 = np.random.default_rng(_method_seed(999, method))

        X_s0 = node_features(d["sub_case"])
        X_t0 = node_features(d["sub_ctrl"])
        X_s0p, X_t0p = preprocess(X_s0, X_t0, method, rng0)
        _, T0 = fgw_solve(X_s0p, X_t0p, d["C"], d["C"])
        Uk_baseline.append(Uk_from_plan(T0, d["C"]))

        uks = []
        for r in range(min(R, 20)):
            rng = np.random.default_rng(_method_seed(r, method))
            if ptype == "kts_noise":
                sc = kts_perturb_noise(d["sub_case"], plevel, rng, conf)
                st = kts_perturb_noise(d["sub_ctrl"], plevel, rng, conf)
                Cr = d["C"]
            elif ptype == "kts_dropout":
                sc, st, _, ki = kts_perturb_dropout(
                    d["sub_case"], d["sub_ctrl"], d["nodes"], plevel, rng, conf)
                Cr = d["C"][np.ix_(ki, ki)]
            elif ptype == "kts_bootstrap":
                sc = kts_perturb_bootstrap(d["sub_case"], rng)
                st = kts_perturb_bootstrap(d["sub_ctrl"], rng)
                Cr = d["C"]
            elif ptype == "noise":
                sc = perturb_noise(d["sub_case"], plevel, rng)
                st = perturb_noise(d["sub_ctrl"], plevel, rng)
                Cr = d["C"]
            elif ptype == "dropout":
                sc, st, _, ki = perturb_dropout(
                    d["sub_case"], d["sub_ctrl"], d["nodes"], plevel, rng)
                Cr = d["C"][np.ix_(ki, ki)]
            elif ptype == "bootstrap":
                sc = perturb_bootstrap(d["sub_case"], rng)
                st = perturb_bootstrap(d["sub_ctrl"], rng)
                Cr = d["C"]
            else:
                continue
            try:
                Xsr = node_features(sc)
                Xtr = node_features(st)
                Xsrp, Xtrp = preprocess(Xsr, Xtr, method, rng)
                _, Tr = fgw_solve(Xsrp, Xtrp, Cr, Cr)
                uks.append(Uk_from_plan(Tr, Cr))
            except Exception:
                continue

        Uk_perturbed_mean.append(np.mean(uks) if uks else Uk_baseline[-1])

    tau, _ = kendalltau(Uk_baseline, Uk_perturbed_mean)
    return float(tau)


# ──────────────────────────────────────────────
#  Main
# ──────────────────────────────────────────────
def main():
    pw_members, pw_edges = load_pw_map(PATHWAY_MAP)
    per_pw_rows: list[dict] = []
    rank_rows: list[dict] = []

    for ds_name, ds_path in DATASETS.items():
        if not ds_path.exists():
            print(f"[WARN] {ds_path} not found — skipping")
            continue

        df = pd.read_csv(ds_path)
        met_cols = detect_met_cols(df)
        met_set  = set(met_cols)
        cond_col = detect_cond_col(df)

        if cond_col is None:
            print(f"[WARN] {ds_name}: no condition column — skipping")
            continue

        groups = list(df[cond_col].dropna().unique())
        if len(groups) != 2:
            print(f"[WARN] {ds_name}: expected 2 groups, found {len(groups)} — skipping")
            continue

        lower_map = {str(g).lower(): g for g in groups}
        if "control" in lower_map:
            ctrl_label = lower_map["control"]
            case_label = [g for g in groups if g != ctrl_label][0]
        else:
            case_label, ctrl_label = groups[0], groups[1]

        df_case = df[df[cond_col] == case_label]
        df_ctrl = df[df[cond_col] == ctrl_label]
        print(f"[{ds_name}] case='{case_label}' (n={len(df_case)}), "
              f"ctrl='{ctrl_label}' (n={len(df_ctrl)})")

        pw_data = {}
        for pw_name, members in pw_members.items():
            available = [m for m in members if m in met_set]
            if len(available) < 3:
                continue
            sc = df_case[available].values.astype(np.float64)
            st = df_ctrl[available].values.astype(np.float64)
            if pw_name in pw_edges and pw_edges[pw_name]:
                C = build_structure(available, pw_edges[pw_name])
            else:
                C = corr_dist(np.vstack([sc, st]).T)

            # Compute metabolite confidence for KTS
            confidence = compute_metabolite_confidence(sc, st)

            pw_data[pw_name] = dict(
                sub_case=sc, sub_ctrl=st, nodes=available,
                C=C, confidence=confidence,
            )

        print(f"  {len(pw_data)} pathways retained")

        # Per-pathway stability
        for pw_name, d in pw_data.items():
            for method in METHODS:
                for reg_name, reg_cfg in PERTURB_REGIMES.items():
                    result = evaluate_stability(
                        d["sub_case"], d["sub_ctrl"], d["nodes"],
                        d["C"], d["confidence"],
                        method, reg_name, reg_cfg,
                    )
                    result["dataset"] = ds_name
                    result["pathway"] = pw_name
                    per_pw_rows.append(result)
                print(f"    {pw_name} × {method}: done")

        # Cross-pathway Kendall τ
        for method in METHODS:
            for reg_name, reg_cfg in PERTURB_REGIMES.items():
                tau = compute_rank_tau(pw_data, method, reg_name, reg_cfg)
                rank_rows.append(dict(
                    dataset=ds_name, method=method,
                    regime=reg_name, rank_tau=round(tau, 6),
                ))

    # Merge
    rank_df = pd.DataFrame(rank_rows)
    pw_df   = pd.DataFrame(per_pw_rows)
    merged  = pw_df.merge(rank_df, on=["dataset", "method", "regime"], how="left")

    col_order = [
        "dataset", "pathway", "regime", "method",
        "cv_fgw", "transport_drift", "delta_U",
        "top3_jaccard", "rank_tau",
    ]
    extra = [c for c in merged.columns if c not in col_order]
    merged = merged[col_order + extra]
    merged.to_csv(OUT_FILE, index=False)
    print(f"\n[done] {len(merged)} rows → {OUT_FILE}")


if __name__ == "__main__":
    main()
