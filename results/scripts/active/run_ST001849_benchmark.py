#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Real multi-pathway identifiability benchmark.

Pipeline:
1. Load metabolomics dataset
2. Parse case/control conditions
3. Restrict to pathway metabolites
4. Build node features per condition
5. Apply selected preprocessing operator
6. Compute FGW-based U_k
7. Run masking-based oracle/reveal benchmark
8. Save regret / top-k outputs

This script is intentionally self-contained.
"""

from __future__ import annotations

from pathlib import Path
import warnings
import numpy as np
import pandas as pd

from scipy.stats import kendalltau
from scipy.spatial.distance import cdist
from sklearn.decomposition import PCA

import ot

# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------

ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "results" / "data"
OUT_DIR = ROOT / "results" / "results"
OUT_DIR.mkdir(parents=True, exist_ok=True)

DATASETS = {
    "ST001849": DATA_DIR / "ST001849_benchmark_ready.csv",
}

# explicit case/control for datasets not matched by choose_case_ctrl() aliases
CASE_CTRL_OVERRIDE = {
    "ST001849": ("severe", "mild"),   # severe = case, mild = control
}

PATHWAY_MAP = DATA_DIR / "ST001849_pathway_mapping.csv"

OUT_FILE = OUT_DIR / "ST001849_benchmark_results.csv"

MASK_RATES = [0.1, 0.2, 0.3, 0.4, 0.5]
TRIALS_PER_RHO = 50
MIN_MET_PER_PATHWAY = 3

# choose the operator selected by your KTS stability benchmark
METHOD_DEFAULT = "none"   # recommended from your current results
# METHOD_DEFAULT = "l2"   # reasonable alternative

JL_DIM = 8
ALPHA_FGW = 0.5
SINK_REG = 0.5
SINK_MAX_ITER = 5000
EPS = 1e-12

warnings.filterwarnings("ignore", message="Sinkhorn did not converge")
warnings.filterwarnings("ignore", message="Warning: numerical errors at iteration")


# ---------------------------------------------------------------------
# Dataset parsing
# ---------------------------------------------------------------------

def fix_st000356(df: pd.DataFrame) -> pd.DataFrame:
    """
    ST000356 arrives in transposed / metadata-row format.

    Converts it to:
        rows = samples
        columns = metabolites
        + condition column
    """
    cols = list(df.columns)
    if "Metabolite_name" not in cols:
        return df

    met_idx = cols.index("Metabolite_name")
    sample_cols = cols[met_idx + 1:]

    # first row contains condition metadata for the sample columns
    cond_row = df.iloc[0][sample_cols]

    def extract_label(x):
        x = str(x).strip().lower()
        if any(k in x for k in ["control", "healthy", "normal"]):
            return "control"
        if any(k in x for k in ["cancer", "tumor", "tumour", "breast cancer"]):
            return "cancer"
        return "unknown"

    conditions = cond_row.apply(extract_label)

    bad = cond_row[conditions == "unknown"]
    if len(bad) > 0:
        print("\n[ST000356] Unknown condition labels detected:")
        print(bad.to_string())

    # remove metadata row, keep only sample columns after Metabolite_name
    df_data = df.iloc[1:].copy()
    df_data = df_data[["Metabolite_name"] + sample_cols]
    df_data = df_data.set_index("Metabolite_name")[sample_cols]

    # transpose so samples are rows
    df_t = df_data.T.reset_index(drop=True)
    df_t = df_t.apply(pd.to_numeric, errors="coerce")

    # attach condition and drop unknown samples
    df_t["condition"] = conditions.values
    df_t = df_t[df_t["condition"] != "unknown"].reset_index(drop=True)

    return df_t


def detect_cond_col(df: pd.DataFrame) -> str | None:
    candidates = [
        "condition", "group", "label",
        "phenotype", "class", "diagnosis", "status"
    ]
    lower_map = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand in lower_map:
            return lower_map[cand]
    return None


def choose_case_ctrl(groups: list) -> tuple:
    groups = list(groups)
    lower_map = {str(g).lower(): g for g in groups}
    aliases = ["control", "ctrl", "normal", "healthy"]

    matched_ctrl = None
    for a in aliases:
        if a in lower_map:
            matched_ctrl = lower_map[a]
            break

    if matched_ctrl is not None:
        ctrl_label = matched_ctrl
        case_label = [g for g in groups if g != ctrl_label][0]
    else:
        case_label, ctrl_label = groups[0], groups[1]

    return case_label, ctrl_label


def detect_met_cols(df: pd.DataFrame) -> list[str]:
    meta_cols = {
        "sample_id", "sampleid", "condition", "group", "label",
        "class", "diagnosis", "status", "phenotype", "refmet_name"
    }
    out = []
    for c in df.columns:
        if c.lower() in meta_cols:
            continue
        # keep numeric-like columns only
        try:
            pd.to_numeric(df[c], errors="raise")
            out.append(c)
        except Exception:
            continue
    return out


# ---------------------------------------------------------------------
# Pathway mapping
# ---------------------------------------------------------------------

def load_pathway_mapping(path: Path):
    df = pd.read_csv(path)
    df.columns = [c.strip().lower() for c in df.columns]

    if "metabolite_name" in df.columns and "pathway_name" in df.columns:
        df = df.rename(columns={
            "metabolite_name": "metabolite",
            "pathway_name": "pathway",
        })

    assert {"metabolite", "pathway"}.issubset(df.columns), df.columns.tolist()

    pw_members = (
        df.groupby("pathway")["metabolite"]
        .apply(lambda x: sorted(set(map(str, x))))
        .to_dict()
    )

    # membership table only; no explicit edge list in current mapping
    pw_edges = {pw: [] for pw in pw_members}
    return pw_members, pw_edges


# ---------------------------------------------------------------------
# Feature construction
# ---------------------------------------------------------------------

def node_features(sub: np.ndarray) -> np.ndarray:
    """
    Condition-specific node features.
    Uses sample-profile embedding (rows=nodes, cols=samples).
    """
    X = sub.T.copy()  # (n_nodes, n_samples)
    X = np.asarray(X, dtype=np.float64)
    X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)

    # per-node standardization for stability
    mu = np.mean(X, axis=1, keepdims=True)
    sd = np.std(X, axis=1, keepdims=True)
    sd[sd < EPS] = 1.0
    X = (X - mu) / sd
    return X


# ---------------------------------------------------------------------
# Structure matrices
# ---------------------------------------------------------------------

def sanitize_square(M: np.ndarray) -> np.ndarray:
    M = np.asarray(M, dtype=np.float64)
    M = np.nan_to_num(M, nan=0.0, posinf=0.0, neginf=0.0)
    M = 0.5 * (M + M.T)
    np.fill_diagonal(M, 0.0)
    return M


def corr_dist(X: np.ndarray) -> np.ndarray:
    """
    X expected as (n_nodes, n_samples) for node-wise correlation.
    """
    X = np.asarray(X, dtype=np.float64)
    X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)

    if X.shape[0] < 2 or X.shape[1] < 2:
        return np.zeros((X.shape[0], X.shape[0]), dtype=np.float64)

    R = np.corrcoef(X)
    R = np.nan_to_num(R, nan=0.0, posinf=0.0, neginf=0.0)
    R = np.clip(R, -1.0, 1.0)
    C = np.sqrt(np.maximum(2.0 * (1.0 - R), 0.0))
    return sanitize_square(C)


# ---------------------------------------------------------------------
# Preprocessing
# ---------------------------------------------------------------------

def preprocess(X_s: np.ndarray, X_t: np.ndarray, method, rng) -> tuple[np.ndarray, np.ndarray]:
    X_s = np.asarray(X_s, dtype=np.float64)
    X_t = np.asarray(X_t, dtype=np.float64)

    d_s = X_s.shape[1]
    d_t = X_t.shape[1]
    d_max = max(d_s, d_t)

    # zero-pad smaller ambient dimension so source and target share a common space
    if d_s < d_max:
        X_s = np.hstack([X_s, np.zeros((X_s.shape[0], d_max - d_s))])
    if d_t < d_max:
        X_t = np.hstack([X_t, np.zeros((X_t.shape[0], d_max - d_t))])

    X_all = np.vstack([X_s, X_t])
    n_s = X_s.shape[0]

    if method == "none":
        Z = X_all

    elif method == "l2":
        norms = np.linalg.norm(X_all, axis=1, keepdims=True)
        norms[norms < EPS] = 1.0
        Z = X_all / norms

    elif method == "jl":
        rng = np.random.default_rng(101)
        d_proj = min(JL_DIM, X_all.shape[1])
        P = rng.normal(0.0, 1.0 / np.sqrt(d_proj), size=(X_all.shape[1], d_proj))
        Z = X_all @ P

    elif method == "randproj":
        rng = np.random.default_rng(202)
        d_proj = min(JL_DIM, X_all.shape[1])
        P = rng.normal(0.0, 1.0, size=(X_all.shape[1], d_proj))
        Z = X_all @ P

    elif method == "pca_fixed":
        d_proj = min(JL_DIM, X_all.shape[0], X_all.shape[1])
        if d_proj < 1:
            Z = X_all
        else:
            pca = PCA(n_components=d_proj, svd_solver="full")
            Z = pca.fit_transform(X_all)

    elif method == "pca_var95":
        pca = PCA(n_components=0.95, svd_solver="full")
        Z = pca.fit_transform(X_all)

    else:
        raise ValueError(f"Unknown method: {method}")

    Z = np.nan_to_num(Z, nan=0.0, posinf=0.0, neginf=0.0)
    return Z[:n_s], Z[n_s:]


# ---------------------------------------------------------------------
# FGW and score
# ---------------------------------------------------------------------

def _n(M: np.ndarray) -> np.ndarray:
    M = sanitize_square(M)
    s = np.mean(M[M > 0]) if np.any(M > 0) else 1.0
    if s < EPS:
        s = 1.0
    return M / s


def fgw_cross(X_s: np.ndarray, X_t: np.ndarray, C_s: np.ndarray, C_t: np.ndarray):
    X_s = np.nan_to_num(np.asarray(X_s, dtype=np.float64), nan=0.0, posinf=0.0, neginf=0.0)
    X_t = np.nan_to_num(np.asarray(X_t, dtype=np.float64), nan=0.0, posinf=0.0, neginf=0.0)
    C_s = _n(C_s)
    C_t = _n(C_t)

    p = np.ones(X_s.shape[0], dtype=np.float64) / max(X_s.shape[0], 1)
    q = np.ones(X_t.shape[0], dtype=np.float64) / max(X_t.shape[0], 1)

    M = cdist(X_s, X_t, metric="euclidean")
    d = X_s.shape[1] if X_s.ndim == 2 else 1
    M = M / (np.sqrt(d) + EPS)
    M = np.nan_to_num(M, nan=0.0, posinf=0.0, neginf=0.0)

    try:
        T, log = ot.gromov.entropic_fused_gromov_wasserstein(
            M, C_s, C_t, p, q,
            loss_fun="square_loss",
            epsilon=SINK_REG,
            alpha=ALPHA_FGW,
            log=True,
            numItermax=SINK_MAX_ITER,
        )
    except TypeError:
        # POT version fallback
        T, log = ot.gromov.entropic_fused_gromov_wasserstein(
            M, C_s, C_t, p, q,
            loss_fun="square_loss",
            epsilon=SINK_REG,
            alpha=ALPHA_FGW,
            log=True,
            max_iter=SINK_MAX_ITER,
        )
    except Exception:
        return np.nan, None

    T = np.asarray(T, dtype=np.float64)
    T = np.nan_to_num(T, nan=0.0, posinf=0.0, neginf=0.0)

    if T.sum() <= 0 or not np.all(np.isfinite(T)):
        return np.nan, None

    dval = float(log.get("fgw_dist", np.nan)) if isinstance(log, dict) else np.nan
    if not np.isfinite(dval):
        dval = float(np.sum(T * M))
    return dval, T


def transport_entropy(T: np.ndarray) -> float:
    T = np.asarray(T, dtype=np.float64)
    T = np.maximum(T, EPS)
    T = T / (T.sum() + EPS)
    H = -np.sum(T * np.log(T + EPS))
    return float(H)


def spectral_gap(C: np.ndarray) -> float:
    C = sanitize_square(C)

    pos = C[C > 0]
    if len(pos) == 0:
        return 0.0

    sigma = float(np.median(pos))
    if not np.isfinite(sigma) or sigma < EPS:
        return 0.0

    A = np.exp(-(C ** 2) / (2.0 * sigma ** 2 + EPS))
    A = np.nan_to_num(A, nan=0.0, posinf=0.0, neginf=0.0)
    A = 0.5 * (A + A.T)
    np.fill_diagonal(A, 0.0)
    A[A < 0] = 0.0

    d = A.sum(axis=1)
    L = np.diag(d) - A
    L = np.nan_to_num(L, nan=0.0, posinf=0.0, neginf=0.0)
    L = 0.5 * (L + L.T)

    try:
        ev = np.sort(np.linalg.eigvalsh(L))
    except np.linalg.LinAlgError:
        return 0.0

    if len(ev) < 2:
        return 0.0

    denom = float(ev[-1]) + EPS
    if not np.isfinite(denom) or abs(denom) < EPS:
        return 0.0

    gap = float(ev[1] / denom)
    if not np.isfinite(gap):
        return 0.0
    return max(gap, 0.0)


def compute_Uk(X_s: np.ndarray, X_t: np.ndarray, C_s: np.ndarray, C_t: np.ndarray) -> float:
    """
    Cross-condition U_k.
    Uses transport entropy + structural gap.
    No instability re-solves here for speed.
    """
    _, T = fgw_cross(X_s, X_t, C_s, C_t)
    if T is None:
        raise ValueError("FGW failed to produce a valid transport plan")
    H = transport_entropy(T)
    sgi = spectral_gap(C_s)
    return float(H + sgi)


# ---------------------------------------------------------------------
# Masking benchmark
# ---------------------------------------------------------------------

def mask_nodes(n_nodes: int, rho: float, rng: np.random.Generator):
    n_hide = max(1, int(round(rho * n_nodes)))
    n_hide = min(n_hide, n_nodes - 2)
    perm = rng.permutation(n_nodes).tolist()
    hidden_idx = sorted(perm[:n_hide])
    kept_idx = sorted(perm[n_hide:])
    return kept_idx, hidden_idx


def _safe_z(v):
    v = np.asarray(v, dtype=np.float64)
    if len(v) == 0:
        return v
    mu = np.mean(v)
    sd = np.std(v)
    if sd < EPS:
        return np.zeros_like(v)
    return (v - mu) / (sd + EPS)


def degree_scores(C_full, hidden_idx):
    C = np.asarray(C_full, dtype=np.float64)
    C = np.nan_to_num(C, nan=0.0, posinf=0.0, neginf=0.0)
    deg = np.sum(C > 0, axis=1).astype(np.float64)
    return {hi: float(deg[hi]) for hi in hidden_idx}


def variance_scores(sub_case, sub_ctrl, hidden_idx):
    pooled = np.vstack([sub_case, sub_ctrl])
    vals = np.nanvar(pooled[:, hidden_idx], axis=0)
    return {hi: float(vals[j]) for j, hi in enumerate(hidden_idx)}


def diff_scores(sub_case, sub_ctrl, hidden_idx):
    mu_case = np.nanmean(sub_case[:, hidden_idx], axis=0)
    mu_ctrl = np.nanmean(sub_ctrl[:, hidden_idx], axis=0)
    vals = np.abs(mu_case - mu_ctrl)
    return {hi: float(vals[j]) for j, hi in enumerate(hidden_idx)}


def surrogate_scores(sub_case, sub_ctrl, C_full, hidden_idx):
    """
    Composite non-oracle surrogate:
    |mean difference| + pooled variance + structural degree
    """
    diff_d = diff_scores(sub_case, sub_ctrl, hidden_idx)
    var_d = variance_scores(sub_case, sub_ctrl, hidden_idx)
    deg_d = degree_scores(C_full, hidden_idx)

    diff_v = np.array([diff_d[i] for i in hidden_idx], dtype=np.float64)
    var_v  = np.array([var_d[i]  for i in hidden_idx], dtype=np.float64)
    deg_v  = np.array([deg_d[i]  for i in hidden_idx], dtype=np.float64)

    diff_z = _safe_z(diff_v)
    var_z  = _safe_z(var_v)
    deg_z  = _safe_z(deg_v)

    comp = diff_z + var_z + deg_z
    return {hi: float(comp[j]) for j, hi in enumerate(hidden_idx)}


def random_scores(hidden_idx, rng):
    vals = rng.random(len(hidden_idx))
    return {hi: float(vals[j]) for j, hi in enumerate(hidden_idx)}


def rank_metrics_from_scores(deltas, pred_scores, nodes):
    """
    deltas: oracle exact ΔUk keyed by hidden node index
    pred_scores: predictor score keyed by hidden node index
    nodes: full node-name list
    """
    if not deltas or not pred_scores:
        return None

    common = [hi for hi in deltas.keys() if hi in pred_scores]
    if len(common) == 0:
        return None

    oracle_sorted = sorted(common, key=lambda hi: deltas[hi], reverse=True)
    pred_sorted   = sorted(common, key=lambda hi: pred_scores[hi], reverse=True)

    hi_star = oracle_sorted[0]
    hi_hat  = pred_sorted[0]

    delta_star = float(deltas[hi_star])
    delta_hat  = float(deltas[hi_hat])

    regret  = delta_star - delta_hat
    nregret = regret / (abs(delta_star) + EPS)

    oracle_top3 = oracle_sorted[:3]
    pred_top3   = pred_sorted[:3]

    top1 = int(hi_hat == hi_star)
    top3 = len(set(oracle_top3) & set(pred_top3)) / 3.0

    oracle_rank_vals = [deltas[hi] for hi in common]
    pred_rank_vals   = [pred_scores[hi] for hi in common]
    tau, _ = kendalltau(oracle_rank_vals, pred_rank_vals)
    tau = float(tau) if np.isfinite(tau) else np.nan

    return {
        "delta_star": round(delta_star, 8),
        "delta_hat": round(delta_hat, 8),
        "regret": round(regret, 8),
        "nregret": round(nregret, 8),
        "top1": top1,
        "top3": round(top3, 4),
        "rank_tau": round(tau, 6) if np.isfinite(tau) else np.nan,
        "m_star": nodes[hi_star],
        "m_hat": nodes[hi_hat],
    }
def masked_transport_proxy(T_mask, C_full, kept_idx, hidden_idx):
    """
    Propagate masked-state transport mass from kept metabolites to hidden candidates.

    Parameters
    ----------
    T_mask : np.ndarray
        FGW transport plan on the masked pathway state, shape (n_kept, n_kept)
        or generally (n_kept, n_target). We only use row mass.
    C_full : np.ndarray
        Full pathway structure matrix, shape (n_nodes, n_nodes).
    kept_idx : list[int]
        Indices of observed/kept metabolites in the full pathway.
    hidden_idx : list[int]
        Indices of hidden candidate metabolites in the full pathway.

    Returns
    -------
    dict[int, float]
        Hidden-node index -> propagated transport proxy
    """
    T_mask = np.asarray(T_mask, dtype=np.float64)
    C_full = np.asarray(C_full, dtype=np.float64)

    # row transport mass on the kept masked state
    tau_kept = np.sum(T_mask, axis=1)
    tau_kept = np.nan_to_num(tau_kept, nan=0.0, posinf=0.0, neginf=0.0)

    out = {}
    if len(kept_idx) == 0:
        return {hi: 0.0 for hi in hidden_idx}

    # bandwidth from positive distances among kept/hidden interactions
    sub = C_full[np.ix_(hidden_idx, kept_idx)] if len(hidden_idx) > 0 else np.array([])
    pos = sub[sub > 0]
    lam = float(np.median(pos)) if len(pos) > 0 else 1.0
    if not np.isfinite(lam) or lam <= EPS:
        lam = 1.0

    for hi in hidden_idx:
        d = C_full[hi, kept_idx]
        d = np.nan_to_num(d, nan=np.inf, posinf=np.inf, neginf=np.inf)

        # similarity weights from hidden candidate to kept nodes
        w = np.exp(-(d ** 2) / (lam ** 2 + EPS))
        w = np.nan_to_num(w, nan=0.0, posinf=0.0, neginf=0.0)

        s = w.sum()
        if s <= EPS:
            out[hi] = 0.0
        else:
            w = w / s
            out[hi] = float(np.dot(w, tau_kept))

    return out


def mb2d_transport_scores(
    sub_case,
    sub_ctrl,
    C_full,
    kept_idx,
    hidden_idx,
    method=METHOD_DEFAULT,
):
    """
    Transport-aware 2D Maxwell-Boltzmann / Rayleigh surrogate.

    For each hidden metabolite m, define a 2D vector:
        u_m = [ z_diff(m), z_tau(m) ]
    where:
        z_diff = standardized |mu_case - mu_ctrl|
        z_tau  = standardized propagated masked-state transport proxy

    Then define speed:
        v_m = ||u_m||_2

    Under an isotropic Gaussian assumption in 2D, v follows the 2D
    Maxwell-Boltzmann speed law (Rayleigh distribution):
        f(v; sigma) = (v / sigma^2) * exp(-v^2 / (2 sigma^2))

    We score candidates by negative log-density:
        score(m) = -log f(v_m; sigma_hat)

    Larger score = more transport-aware / more surprising candidate.
    """
    if len(hidden_idx) == 0 or len(kept_idx) < 2:
        return {}

    # 1) masked baseline state
    X_s_mask = node_features(sub_case[:, kept_idx])
    X_t_mask = node_features(sub_ctrl[:, kept_idx])
    X_s_mask, X_t_mask = preprocess(
        X_s_mask, X_t_mask, method, np.random.default_rng(0)
    )
    C_mask = C_full[np.ix_(kept_idx, kept_idx)]

    # 2) masked FGW transport plan
    _, T_mask = fgw_cross(X_s_mask, X_t_mask, C_mask, C_mask)
    if T_mask is None:
        return {hi: 0.0 for hi in hidden_idx}

    # 3) diff-abundance component on hidden candidates
    mu_case = np.nanmean(sub_case[:, hidden_idx], axis=0)
    mu_ctrl = np.nanmean(sub_ctrl[:, hidden_idx], axis=0)
    diff_v = np.abs(mu_case - mu_ctrl)
    diff_v = np.nan_to_num(diff_v, nan=0.0, posinf=0.0, neginf=0.0)

    # 4) transport-propagated component on hidden candidates
    tau_d = masked_transport_proxy(T_mask, C_full, kept_idx, hidden_idx)
    tau_v = np.array([tau_d[hi] for hi in hidden_idx], dtype=np.float64)

    # 5) z-score both components across hidden candidates
    z_diff = _safe_z(diff_v)
    z_tau = _safe_z(tau_v)

    # 6) 2D speed
    v = np.sqrt(z_diff ** 2 + z_tau ** 2)
    v = np.nan_to_num(v, nan=0.0, posinf=0.0, neginf=0.0)

    # 7) Rayleigh / MB-2D scale estimate
    sigma_hat = np.sqrt(np.mean(v ** 2) / 2.0 + EPS)
    if not np.isfinite(sigma_hat) or sigma_hat <= EPS:
        sigma_hat = 1.0

    # 8) negative log-density score
    # f(v;sigma) = v/sigma^2 * exp(-v^2/(2 sigma^2))
    # -log f = -log(v) + 2log(sigma) + v^2/(2 sigma^2)
    scores = {}
    for j, hi in enumerate(hidden_idx):
        vv = max(float(v[j]), EPS)
        score = (
            -np.log(vv)
            + 2.0 * np.log(sigma_hat + EPS)
            + (vv ** 2) / (2.0 * sigma_hat ** 2 + EPS)
        )
        if not np.isfinite(score):
            score = 0.0
        scores[hi] = float(score)

    return scores
def gnc_commutator_scores(
    sub_case,
    sub_ctrl,
    C_full,
    kept_idx,
    hidden_idx,
):
    """
    Geometry-noncommutative masked-state surrogate.

    For each hidden metabolite m, define:
      1) masked operator L_mask on kept nodes
      2) local insertion operator E_m = a_m * (w_m w_m^T)
         where:
           - w_m encodes structural proximity of hidden m to kept nodes
           - a_m is a scalar candidate strength from diff abundance + variance

    Score:
        s_GNC(m) = || [L_mask, E_m] ||_F
                 = || L_mask E_m - E_m L_mask ||_F

    Larger score means the hidden metabolite is more noncommutatively
    incompatible with the current masked operator, and is hypothesized
    to reduce ambiguity more strongly when revealed.
    """
    if len(hidden_idx) == 0 or len(kept_idx) < 2:
        return {}

    # masked structure operator on kept nodes
    Ck = np.asarray(C_full[np.ix_(kept_idx, kept_idx)], dtype=np.float64)
    Ck = np.nan_to_num(Ck, nan=0.0, posinf=0.0, neginf=0.0)
    Ck = 0.5 * (Ck + Ck.T)
    np.fill_diagonal(Ck, 0.0)
    Ck[Ck < 0] = 0.0

    # similarity kernel from structure distances
    pos = Ck[Ck > 0]
    sigma = float(np.median(pos)) if len(pos) > 0 else 1.0
    if not np.isfinite(sigma) or sigma <= EPS:
        sigma = 1.0

    A = np.exp(-(Ck ** 2) / (2.0 * sigma ** 2 + EPS))
    A = np.nan_to_num(A, nan=0.0, posinf=0.0, neginf=0.0)
    A = 0.5 * (A + A.T)
    np.fill_diagonal(A, 0.0)
    A[A < 0] = 0.0

    D = np.diag(A.sum(axis=1))
    L_mask = D - A
    L_mask = np.nan_to_num(L_mask, nan=0.0, posinf=0.0, neginf=0.0)
    L_mask = 0.5 * (L_mask + L_mask.T)

    # scalar candidate strengths from hidden-node local statistics
    diff_d = diff_scores(sub_case, sub_ctrl, hidden_idx)
    var_d = variance_scores(sub_case, sub_ctrl, hidden_idx)

    diff_v = np.array([diff_d[hi] for hi in hidden_idx], dtype=np.float64)
    var_v = np.array([var_d[hi] for hi in hidden_idx], dtype=np.float64)

    diff_z = _safe_z(diff_v)
    var_z = _safe_z(var_v)

    out = {}

    # bandwidth for hidden-to-kept distances
    hk = C_full[np.ix_(hidden_idx, kept_idx)] if len(hidden_idx) > 0 else np.array([])
    pos_hk = hk[hk > 0]
    lam = float(np.median(pos_hk)) if len(pos_hk) > 0 else 1.0
    if not np.isfinite(lam) or lam <= EPS:
        lam = 1.0

    for j, hi in enumerate(hidden_idx):
        # local scalar strength
        a_m = float(diff_z[j] + var_z[j])

        # structural weights from hidden node to kept nodes
        d = np.asarray(C_full[hi, kept_idx], dtype=np.float64)
        d = np.nan_to_num(d, nan=np.inf, posinf=np.inf, neginf=np.inf)

        w = np.exp(-(d ** 2) / (lam ** 2 + EPS))
        w = np.nan_to_num(w, nan=0.0, posinf=0.0, neginf=0.0)

        s = w.sum()
        if s <= EPS:
            out[hi] = 0.0
            continue

        w = w / s
        w = w.reshape(-1, 1)  # column vector on kept space

        # rank-one insertion operator
        E_m = a_m * (w @ w.T)
        E_m = np.nan_to_num(E_m, nan=0.0, posinf=0.0, neginf=0.0)
        E_m = 0.5 * (E_m + E_m.T)

        # commutator score
        comm = L_mask @ E_m - E_m @ L_mask
        score = np.linalg.norm(comm, ord="fro")

        if not np.isfinite(score):
            score = 0.0

        out[hi] = float(score)

    return out

def run_trial(sub_case, sub_ctrl, nodes, C_full, rho, rng, method=METHOD_DEFAULT):
    kept_idx, hidden_idx = mask_nodes(len(nodes), rho, rng)
    if len(hidden_idx) == 0 or len(kept_idx) < 2:
        return []

    # masked baseline
    X_s_mask = node_features(sub_case[:, kept_idx])
    X_t_mask = node_features(sub_ctrl[:, kept_idx])
    X_s_mask, X_t_mask = preprocess(
        X_s_mask, X_t_mask, method, np.random.default_rng(0)
    )
    C_mask = C_full[np.ix_(kept_idx, kept_idx)]

    try:
        Uk_base = compute_Uk(X_s_mask, X_t_mask, C_mask, C_mask)
    except Exception:
        return []

    # ------------------------------------------------------------
    # 1) exact oracle deltas by full reveal
    # ------------------------------------------------------------
    deltas = {}
    for hi in hidden_idx:
        aug_idx = kept_idx + [hi]

        X_s_aug = node_features(sub_case[:, aug_idx])
        X_t_aug = node_features(sub_ctrl[:, aug_idx])
        X_s_aug, X_t_aug = preprocess(
            X_s_aug, X_t_aug, method, np.random.default_rng(0)
        )
        C_aug = C_full[np.ix_(aug_idx, aug_idx)]

        try:
            Uk_aug = compute_Uk(X_s_aug, X_t_aug, C_aug, C_aug)
        except Exception:
            continue

        delta = Uk_base - Uk_aug
        if np.isfinite(delta):
            deltas[hi] = float(delta)

    if len(deltas) == 0:
        return []

    valid_hidden = sorted(deltas.keys())

         
    # 2) predictor scores WITHOUT exact reveal
    # ------------------------------------------------------------
    model_scores = surrogate_scores(sub_case, sub_ctrl, C_full, valid_hidden)
    mb2d_scores  = mb2d_transport_scores(
        sub_case, sub_ctrl, C_full, kept_idx, valid_hidden, method=method
    )
    gnc_scores   = gnc_commutator_scores(
        sub_case, sub_ctrl, C_full, kept_idx, valid_hidden
    )
    rand_scores  = random_scores(valid_hidden, rng)
    deg_scores   = degree_scores(C_full, valid_hidden)
    var_scores   = variance_scores(sub_case, sub_ctrl, valid_hidden)
    diff_s       = diff_scores(sub_case, sub_ctrl, valid_hidden)

    predictors = {
        "surrogate": model_scores,
        "gnc_commutator": gnc_scores,
        "mb2d_transport": mb2d_scores,
        "random": rand_scores,
        "degree": deg_scores,
        "variance": var_scores,
        "diffabundance": diff_s,
    }
    

    # ------------------------------------------------------------
    # 3) compute regret / top-k / rank-correlation for each method
    # ------------------------------------------------------------
    rows = []
    for pred_name, pred_scores in predictors.items():
        metrics = rank_metrics_from_scores(deltas, pred_scores, nodes)
        if metrics is None:
            continue

        row = {
            "rho": rho,
            "n_nodes": len(nodes),
            "n_hidden": len(valid_hidden),
            "operator_method": method,      # the FGW/preprocess operator
            "predictor_method": pred_name,  # the prioritization rule
        }
        row.update(metrics)
        rows.append(row)

    return rows

# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():
    pw_members, pw_edges = load_pathway_mapping(PATHWAY_MAP)
    rows = []

    for ds_name, ds_path in DATASETS.items():
        if not ds_path.exists():
            print(f"[WARN] {ds_path} not found — skipping")
            continue

        df = pd.read_csv(ds_path)
        if ds_name == "ST000356":
            df = fix_st000356(df)

        cond_col = detect_cond_col(df)
        if cond_col is None:
            print(f"[WARN] {ds_name}: no condition column — cannot run cross-condition benchmark, skipping")
            continue

        groups = list(df[cond_col].dropna().unique())
        if len(groups) != 2:
            print(f"[WARN] {ds_name}: expected 2 groups, found {len(groups)} — skipping")
            continue

        case_label, ctrl_label = choose_case_ctrl(groups)
        if ds_name in CASE_CTRL_OVERRIDE:
            case_label, ctrl_label = CASE_CTRL_OVERRIDE[ds_name]
            print(f"  [override] case='{case_label}', ctrl='{ctrl_label}'")
        df_case = df[df[cond_col] == case_label].copy()
        df_ctrl = df[df[cond_col] == ctrl_label].copy()

        print(f"[{ds_name}] case='{case_label}' (n={len(df_case)}), ctrl='{ctrl_label}' (n={len(df_ctrl)})")

        met_cols = detect_met_cols(df)
        met_set = set(met_cols)

        for pw_name, members in pw_members.items():
            available = [m for m in members if m in met_set]
            if len(available) < MIN_MET_PER_PATHWAY:
                continue

            sub_case = df_case[available].apply(pd.to_numeric, errors="coerce").values.astype(np.float64)
            sub_ctrl = df_ctrl[available].apply(pd.to_numeric, errors="coerce").values.astype(np.float64)

            # drop metabolites that are entirely non-finite in either group
            good = []
            for j in range(len(available)):
                if np.all(~np.isfinite(sub_case[:, j])) or np.all(~np.isfinite(sub_ctrl[:, j])):
                    continue
                good.append(j)

            if len(good) < MIN_MET_PER_PATHWAY:
                continue

            available = [available[j] for j in good]
            sub_case = sub_case[:, good]
            sub_ctrl = sub_ctrl[:, good]

            if pw_name in pw_edges and pw_edges[pw_name]:
                # current mapping has no edges; fallback kept for future use
                C_full = np.zeros((len(available), len(available)), dtype=np.float64)
            else:
                C_full = corr_dist(np.vstack([sub_case, sub_ctrl]).T)

            n_trials_total = 0
            for rho in MASK_RATES:
                ok_trials = 0
                for t in range(TRIALS_PER_RHO):
                    rng = np.random.default_rng(1000 + t + int(100 * rho))
                    trial_rows = run_trial(
                        sub_case, sub_ctrl, available, C_full, rho, rng, method=METHOD_DEFAULT
                    )
                    if not trial_rows:
                        continue
                    
                    for result in trial_rows:
                        result["dataset"] = ds_name
                        result["pathway"] = pw_name
                        rows.append(result)
                    ok_trials += 1
                n_trials_total += ok_trials

            if n_trials_total > 0:
                print(f"  {pw_name}: {n_trials_total} trials")

    if not rows:
        print("[WARN] no valid benchmark rows produced")
        return

    out = pd.DataFrame(rows)
    col_order = [
          "dataset",
          "pathway",
          "rho",
          "operator_method",
          "predictor_method",
          "n_nodes",
          "n_hidden",
          "delta_star",
          "delta_hat",
          "regret",
          "nregret",
          "top1",
          "top3",
          "rank_tau",
          "m_star",
          "m_hat",
       ]
    out = out[[c for c in col_order if c in out.columns]]
    out.to_csv(OUT_FILE, index=False)
    print(f"\n[done] {len(out)} rows → {OUT_FILE}")


if __name__ == "__main__":
    main()
