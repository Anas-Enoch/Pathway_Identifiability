#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_acquisition_simulation.py  (v2 — fixed parameters)

Key changes vs v1:
  INITIAL_MASK_RATE  0.20 → 0.40  (2+ hidden nodes; strategies can now differ)
  MIN_MET            5    → 3     (more pathways retained)
  MAX_NODES          25   → 50    (larger pathways allowed)
  N_STEPS            10   → 5     (max reveals per trial; enough to see curves)
  SINK_MAX_ITER      5000 → 300   (performance)

Added metric: oracle_top1 — whether each strategy picked the metabolite
with the highest oracle ΔU_k (directly connects simulation to benchmark).

Outputs:
  results/results/simulation_delta_uk.csv
  results/results/simulation_auc.csv
  results/results/simulation_pathway_tau.csv
"""

from __future__ import annotations
import warnings

warnings.filterwarnings("ignore", message="Sinkhorn did not converge")
warnings.filterwarnings("ignore", message="n_jobs")
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", message="numerical errors")
warnings.filterwarnings("ignore", message="Solver failed")
warnings.filterwarnings("ignore", message="divide by zero")
try:
    from scipy.stats import SmallSampleWarning
    warnings.filterwarnings("ignore", category=SmallSampleWarning)
except ImportError:
    warnings.filterwarnings("ignore", message="too small")

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import kendalltau
from scipy.spatial.distance import cdist
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.exceptions import ConvergenceWarning
warnings.filterwarnings("ignore", category=ConvergenceWarning)
import ot

# ── Path resolution ────────────────────────────────────────────────────────────
def find_root(start: Path, marker="results", max_depth=6):
    p = start.resolve()
    for _ in range(max_depth):
        if (p / marker).is_dir(): return p
        p = p.parent
    return start.resolve().parents[2]

ROOT     = find_root(Path(__file__).parent)
DATA_DIR = ROOT / "results" / "data"
OUT_DIR  = ROOT / "results" / "results"
OUT_DIR.mkdir(parents=True, exist_ok=True)

DATASET_PATH = DATA_DIR / "processed_metabolite_matrix_ST000356.csv"
PATHWAY_MAP  = DATA_DIR / "core_pathway_mapping.csv"

# ── Config ─────────────────────────────────────────────────────────────────────
INITIAL_MASK_RATE = 0.40   # ← was 0.20; 40% gives 2+ hidden nodes for 5-node pathways
N_STEPS           = 5      # max sequential reveals per trial
N_TRIALS          = 50
MIN_MET           = 3      # ← was 5; allows more pathways
MAX_NODES         = 50     # ← was 25; allows medium pathways
TEST_FRAC         = 0.30
SEED_BASE         = 42

ALPHA_FGW     = 0.5
SINK_REG      = 0.5
SINK_MAX_ITER = 300
EPS           = 1e-12
LR_MAX_ITER   = 300
LR_C          = 1.0

STRATEGIES = ["commutator", "variance", "random"]

# ── Data utilities ─────────────────────────────────────────────────────────────
def fix_st000356(df):
    cols = list(df.columns)
    if "Metabolite_name" not in cols: return df
    idx = cols.index("Metabolite_name"); s_c = cols[idx+1:]
    cr  = df.iloc[0][s_c]
    def lbl(x):
        x = str(x).strip().lower()
        if any(k in x for k in ["control","healthy","normal"]): return "control"
        if any(k in x for k in ["cancer","tumor","tumour"]):    return "cancer"
        return "unknown"
    conds = cr.apply(lbl)
    df2   = df.iloc[1:].copy()[["Metabolite_name"]+s_c].set_index("Metabolite_name")[s_c]
    dft   = df2.T.reset_index(drop=True).apply(pd.to_numeric, errors="coerce")
    dft["condition"] = conds.values
    return dft[dft["condition"] != "unknown"].reset_index(drop=True)

def detect_cond_col(df):
    for c in ["condition","group","label","phenotype","class","diagnosis","status"]:
        for col in df.columns:
            if col.lower() == c: return col
    return None

def choose_case_ctrl(groups):
    groups = list(groups); lm = {str(g).lower(): g for g in groups}
    for a in ["control","ctrl","normal","healthy"]:
        if a in lm:
            ctrl = lm[a]; case = [g for g in groups if g != ctrl][0]
            return case, ctrl
    return groups[0], groups[1]

def detect_met_cols(df):
    skip = {"sample_id","sampleid","condition","group","label",
            "class","diagnosis","status","phenotype","refmet_name","-"}
    out = []
    for c in df.columns:
        if c.lower() in skip or c.strip() == "-": continue
        try: pd.to_numeric(df[c], errors="raise"); out.append(c)
        except: pass
    return out

def load_pathway_mapping(path):
    df = pd.read_csv(path); df.columns = [c.strip().lower() for c in df.columns]
    if "metabolite_name" in df.columns: df = df.rename(columns={"metabolite_name":"metabolite"})
    if "pathway_name"    in df.columns: df = df.rename(columns={"pathway_name":"pathway"})
    return df.groupby("pathway")["metabolite"].apply(
        lambda x: sorted(set(map(str, x)))).to_dict()

# ── Math utilities ─────────────────────────────────────────────────────────────
def sanitize_square(M):
    M = np.asarray(M, dtype=np.float64)
    M = np.nan_to_num(M, nan=0., posinf=0., neginf=0.)
    M = 0.5*(M+M.T); np.fill_diagonal(M, 0.); return M

def corr_distance_structure(sub_case, sub_ctrl):
    X = np.nan_to_num(np.vstack([sub_case, sub_ctrl]).astype(np.float64),
                      nan=0., posinf=0., neginf=0.)
    if X.shape[1] < 2: return np.zeros((X.shape[1], X.shape[1]))
    R = np.clip(np.corrcoef(X, rowvar=False), -1., 1.)
    R = np.nan_to_num(R, nan=0., posinf=0., neginf=0.)
    D = np.sqrt(np.maximum(2.*(1.-R), 0.))
    np.fill_diagonal(D, 0.); D = 0.5*(D+D.T)
    return D

def node_features(X):
    X = np.asarray(X, dtype=float)
    return np.column_stack([np.nanmean(X, axis=0), np.nanstd(X, axis=0)])

def fgw_cross(X_s, X_t, C_s, C_t):
    def _n(M):
        M = sanitize_square(M); s = np.mean(M[M>0]) if np.any(M>0) else 1.
        return M / (s if s > EPS else 1.)
    X_s = np.nan_to_num(np.asarray(X_s, dtype=np.float64), nan=0., posinf=0., neginf=0.)
    X_t = np.nan_to_num(np.asarray(X_t, dtype=np.float64), nan=0., posinf=0., neginf=0.)
    p = np.ones(X_s.shape[0]) / max(X_s.shape[0], 1)
    q = np.ones(X_t.shape[0]) / max(X_t.shape[0], 1)
    M = cdist(X_s, X_t, metric="euclidean")
    d = X_s.shape[1] if X_s.ndim == 2 else 1
    M = np.nan_to_num(M / (np.sqrt(d) + EPS), nan=0., posinf=0., neginf=0.)
    try:
        T, log = ot.gromov.entropic_fused_gromov_wasserstein(
            M, _n(C_s), _n(C_t), p, q, loss_fun="square_loss",
            epsilon=SINK_REG, alpha=ALPHA_FGW, log=True, numItermax=SINK_MAX_ITER)
    except TypeError:
        T, log = ot.gromov.entropic_fused_gromov_wasserstein(
            M, _n(C_s), _n(C_t), p, q, loss_fun="square_loss",
            epsilon=SINK_REG, alpha=ALPHA_FGW, log=True, max_iter=SINK_MAX_ITER)
    except Exception:
        return np.nan, None
    T = np.nan_to_num(np.asarray(T, dtype=np.float64), nan=0., posinf=0., neginf=0.)
    if T.sum() <= 0 or not np.all(np.isfinite(T)): return np.nan, None
    dval = float(log.get("fgw_dist", np.nan)) if isinstance(log, dict) else np.nan
    if not np.isfinite(dval): dval = float(np.sum(T * M))
    return dval, T

def spectral_gap(C):
    C = sanitize_square(C); pos = C[C > 0]
    if len(pos) == 0: return 0.
    sigma = float(np.median(pos))
    if not np.isfinite(sigma) or sigma < EPS: return 0.
    A = np.exp(-(C**2) / (2.*sigma**2 + EPS))
    A = np.nan_to_num(A, nan=0., posinf=0., neginf=0.)
    A = 0.5*(A+A.T); np.fill_diagonal(A, 0.); A[A < 0] = 0.
    L = np.diag(A.sum(axis=1)) - A
    L = np.nan_to_num(L, nan=0., posinf=0., neginf=0.); L = 0.5*(L+L.T)
    try: ev = np.sort(np.linalg.eigvalsh(L))
    except: return 0.
    if len(ev) < 2: return 0.
    denom = float(ev[-1]) + EPS
    return max(float(ev[1]/denom), 0.) if np.isfinite(denom) and abs(denom) > EPS else 0.

def transport_entropy(T):
    T = np.maximum(np.asarray(T, dtype=np.float64), EPS)
    T = T / (T.sum() + EPS)
    return float(-np.sum(T * np.log(T + EPS)))

def compute_Uk(X_s, X_t, C_s, C_t):
    _, T = fgw_cross(X_s, X_t, C_s, C_t)
    if T is None: raise ValueError("FGW failed")
    return float(transport_entropy(T) + spectral_gap(C_s))

def compute_Uk_from_panel(sub_case, sub_ctrl, panel_idx, C_full):
    return compute_Uk(
        node_features(sub_case[:, panel_idx]),
        node_features(sub_ctrl[:, panel_idx]),
        C_full[np.ix_(panel_idx, panel_idx)],
        C_full[np.ix_(panel_idx, panel_idx)],
    )

def mask_nodes(n_nodes, rho, rng):
    n_hide = max(1, min(int(round(rho * n_nodes)), n_nodes - 2))
    perm   = rng.permutation(n_nodes).tolist()
    return sorted(perm[n_hide:]), sorted(perm[:n_hide])   # kept, hidden

# ── Predictors ─────────────────────────────────────────────────────────────────
def gnc_commutator_score(sub_case, sub_ctrl, C_full, kept_idx, candidate_idx):
    Ck = np.asarray(C_full[np.ix_(kept_idx, kept_idx)], dtype=np.float64)
    Ck = np.nan_to_num(Ck, nan=0., posinf=0., neginf=0.)
    Ck = 0.5*(Ck+Ck.T); np.fill_diagonal(Ck, 0.); Ck[Ck < 0] = 0.
    pos = Ck[Ck > 0]
    sigma = float(np.median(pos)) if len(pos) > 0 else 1.
    if not np.isfinite(sigma) or sigma <= EPS: sigma = 1.
    A = np.exp(-(Ck**2) / (2.*sigma**2 + EPS))
    A = np.nan_to_num(A, nan=0., posinf=0., neginf=0.)
    A = 0.5*(A+A.T); np.fill_diagonal(A, 0.); A[A < 0] = 0.
    L_mask = np.diag(A.sum(axis=1)) - A
    L_mask = np.nan_to_num(L_mask, nan=0., posinf=0., neginf=0.)
    L_mask = 0.5*(L_mask+L_mask.T)
    mu_c = np.nanmean(sub_case[:, candidate_idx])
    mu_t = np.nanmean(sub_ctrl[:, candidate_idx])
    var  = np.nanvar(np.concatenate([sub_case[:, candidate_idx],
                                      sub_ctrl[:, candidate_idx]]))
    a_m  = float(abs(mu_c - mu_t) + np.sqrt(max(var, 0.)))
    d    = C_full[candidate_idx, kept_idx].astype(np.float64)
    d    = np.nan_to_num(d, nan=np.inf, posinf=np.inf, neginf=np.inf)
    hk   = C_full[np.ix_([candidate_idx], kept_idx)]
    pos_hk = hk[hk > 0]
    lam  = float(np.median(pos_hk)) if len(pos_hk) > 0 else 1.
    if not np.isfinite(lam) or lam <= EPS: lam = 1.
    w = np.exp(-(d**2) / (lam**2 + EPS))
    w = np.nan_to_num(w, nan=0., posinf=0., neginf=0.)
    s = w.sum()
    if s <= EPS: return 0.
    w = (w/s).reshape(-1, 1)
    E_m  = a_m * (w @ w.T)
    E_m  = np.nan_to_num(E_m, nan=0., posinf=0., neginf=0.)
    E_m  = 0.5*(E_m+E_m.T)
    comm = L_mask @ E_m - E_m @ L_mask
    score = np.linalg.norm(comm, ord="fro")
    return float(score) if np.isfinite(score) else 0.

def variance_score(sub_case, sub_ctrl, candidate_idx):
    pooled = np.concatenate([sub_case[:, candidate_idx], sub_ctrl[:, candidate_idx]])
    v = float(np.nanvar(pooled))
    return v if np.isfinite(v) else 0.

def select_next(strategy, sub_case, sub_ctrl, C_full, kept_idx, hidden_idx, rng):
    if strategy == "random":
        return int(rng.choice(hidden_idx))
    elif strategy == "variance":
        scores = {hi: variance_score(sub_case, sub_ctrl, hi) for hi in hidden_idx}
    elif strategy == "commutator":
        scores = {hi: gnc_commutator_score(sub_case, sub_ctrl, C_full, kept_idx, hi)
                  for hi in hidden_idx}
    else:
        raise ValueError(f"Unknown strategy: {strategy}")
    return max(scores, key=lambda hi: scores[hi])

def compute_oracle_deltas(sub_case, sub_ctrl, kept_idx, hidden_idx, C_full):
    """Exact ΔU_k for each hidden candidate. Returns dict {hi: delta}."""
    try:
        Uk_base = compute_Uk_from_panel(sub_case, sub_ctrl, kept_idx, C_full)
    except Exception:
        return {}
    deltas = {}
    for hi in hidden_idx:
        aug = sorted(kept_idx + [hi])
        try:
            Uk_aug = compute_Uk_from_panel(sub_case, sub_ctrl, aug, C_full)
            delta  = Uk_base - Uk_aug
            if np.isfinite(delta):
                deltas[hi] = float(delta)
        except Exception:
            pass
    return deltas

# ── AUC ────────────────────────────────────────────────────────────────────────
def compute_auc(sc_tr, sk_tr, sc_te, sk_te, panel_idx):
    X_tr = np.nan_to_num(np.vstack([sc_tr[:, panel_idx],
                                     sk_tr[:, panel_idx]]).astype(np.float64),
                         nan=0., posinf=0., neginf=0.)
    y_tr = np.array([1]*len(sc_tr) + [0]*len(sk_tr), dtype=int)
    X_te = np.nan_to_num(np.vstack([sc_te[:, panel_idx],
                                     sk_te[:, panel_idx]]).astype(np.float64),
                         nan=0., posinf=0., neginf=0.)
    y_te = np.array([1]*len(sc_te) + [0]*len(sk_te), dtype=int)
    if len(np.unique(y_te)) < 2 or X_tr.shape[1] == 0: return np.nan
    lr = LogisticRegression(C=LR_C, max_iter=LR_MAX_ITER, solver="lbfgs", random_state=0)
    try:
        lr.fit(X_tr, y_tr)
        return float(roc_auc_score(y_te, lr.predict_proba(X_te)[:, 1]))
    except Exception:
        return np.nan

# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    if not DATASET_PATH.exists():
        print(f"[ERROR] Dataset not found: {DATASET_PATH}"); return

    df = fix_st000356(pd.read_csv(DATASET_PATH))
    cond_col = detect_cond_col(df)
    if cond_col is None:
        print("[ERROR] No condition column found"); return
    groups = list(df[cond_col].dropna().unique())
    if len(groups) != 2:
        print(f"[ERROR] Expected 2 groups, got {groups}"); return
    case_label, ctrl_label = choose_case_ctrl(groups)
    df_case = df[df[cond_col] == case_label].copy()
    df_ctrl = df[df[cond_col] == ctrl_label].copy()
    print(f"[ST000356] case='{case_label}' n={len(df_case)}, "
          f"ctrl='{ctrl_label}' n={len(df_ctrl)}")
    print(f"Parameters: mask_rate={INITIAL_MASK_RATE}, MIN_MET={MIN_MET}, "
          f"MAX_NODES={MAX_NODES}, N_TRIALS={N_TRIALS}, N_STEPS={N_STEPS}")

    met_set    = set(detect_met_cols(df))
    pw_members = load_pathway_mapping(PATHWAY_MAP)

    # ── Build pathway data ────────────────────────────────────────────────────
    print("\nBuilding pathway data ...")
    pathway_list = []
    pathway_data = {}
    full_panel_uk = {}

    for pw_name, members in pw_members.items():
        available = [m for m in members if m in met_set]
        if len(available) < MIN_MET: continue
        if len(available) > MAX_NODES:
            continue  # silent skip for large pathways

        sc = df_case[available].apply(pd.to_numeric, errors="coerce").values.astype(np.float64)
        sk = df_ctrl[available].apply(pd.to_numeric, errors="coerce").values.astype(np.float64)
        good = [j for j in range(len(available))
                if not (np.all(~np.isfinite(sc[:, j])) or
                        np.all(~np.isfinite(sk[:, j])))]
        if len(good) < MIN_MET: continue
        available = [available[j] for j in good]
        sc = sc[:, good]; sk = sk[:, good]

        # Check masking will give >= 2 hidden nodes
        n_hide = max(1, min(int(round(INITIAL_MASK_RATE * len(available))),
                            len(available) - 2))
        if n_hide < 2:
            print(f"  [skip] {pw_name}: {len(available)} nodes → only {n_hide} hidden")
            continue

        C_full = corr_distance_structure(sc, sk)
        try:
            Uk_full = compute_Uk_from_panel(sc, sk, list(range(len(available))), C_full)
            if not np.isfinite(Uk_full): continue
        except Exception:
            continue

        full_panel_uk[pw_name] = float(Uk_full)
        pathway_list.append(pw_name)
        pathway_data[pw_name] = (available, sc, sk, C_full)
        print(f"  ✓ {pw_name} ({len(available)} nodes, {n_hide} hidden)  "
              f"Uk_full={Uk_full:.4f}")

    print(f"\nSimulation-ready pathways: {len(pathway_list)}")
    if not pathway_list:
        print("[ERROR] No pathways retained."); return

    full_panel_ranking = sorted(full_panel_uk, key=full_panel_uk.get, reverse=True)

    # ── Fixed test/train split ────────────────────────────────────────────────
    rng_split = np.random.default_rng(SEED_BASE)
    n_case, n_ctrl = len(df_case), len(df_ctrl)
    n_test_case  = max(1, int(round(TEST_FRAC * n_case)))
    n_test_ctrl  = max(1, int(round(TEST_FRAC * n_ctrl)))
    case_test_idx  = sorted(rng_split.choice(n_case, n_test_case, replace=False))
    ctrl_test_idx  = sorted(rng_split.choice(n_ctrl, n_test_ctrl, replace=False))
    case_train_idx = [i for i in range(n_case) if i not in set(case_test_idx)]
    ctrl_train_idx = [i for i in range(n_ctrl) if i not in set(ctrl_test_idx)]

    # ── Simulation ────────────────────────────────────────────────────────────
    rows_delta_uk = []
    rows_auc      = []
    rows_path_tau = []

    print("\nRunning simulation ...")
    for trial_idx in range(N_TRIALS):
        rng        = np.random.default_rng(SEED_BASE + trial_idx * 997)
        partial_uk = {}

        for pw_name in pathway_list:
            available, sc, sk, C_full = pathway_data[pw_name]
            n_nodes = len(available)

            sc_tr = sc[case_train_idx, :]; sc_te = sc[case_test_idx, :]
            sk_tr = sk[ctrl_train_idx, :]; sk_te = sk[ctrl_test_idx, :]

            init_kept, init_hidden = mask_nodes(n_nodes, INITIAL_MASK_RATE, rng)
            if len(init_hidden) < 2: continue   # enforce ≥2 hidden

            # Compute oracle deltas once for this masking (shared across strategies)
            oracle_deltas = compute_oracle_deltas(sc, sk, init_kept, init_hidden, C_full)
            if len(oracle_deltas) < 2: continue
            oracle_best = max(oracle_deltas, key=oracle_deltas.get)  # oracle top-1

            # Initial U_k for pathway ranking
            try:
                Uk_init = compute_Uk_from_panel(sc, sk, init_kept, C_full)
                if not np.isfinite(Uk_init): continue
            except Exception:
                continue
            partial_uk[pw_name] = float(Uk_init)

            # Run each strategy
            for strategy in STRATEGIES:
                kept   = list(init_kept)
                hidden = list(init_hidden)
                Uk_prev = Uk_init
                s_rng  = np.random.default_rng(
                    SEED_BASE + trial_idx*997 + hash(strategy) % 1000)

                for step in range(min(N_STEPS, len(init_hidden))):
                    if len(hidden) < 1: break

                    chosen = select_next(strategy, sc, sk, C_full, kept, hidden, s_rng)
                    is_oracle_top1 = int(chosen == oracle_best and step == 0)

                    hidden.remove(chosen)
                    kept = sorted(kept + [chosen])

                    try:
                        Uk_new = compute_Uk_from_panel(sc, sk, kept, C_full)
                        if not np.isfinite(Uk_new): continue
                    except Exception:
                        continue

                    delta_uk = float(Uk_prev - Uk_new)
                    Uk_prev  = Uk_new

                    rows_delta_uk.append({
                        "trial":         trial_idx,
                        "pathway":       pw_name,
                        "n_nodes":       n_nodes,
                        "n_hidden_init": len(init_hidden),
                        "strategy":      strategy,
                        "step":          step + 1,
                        "delta_uk":      delta_uk,
                        "uk_current":    Uk_new,
                        "oracle_top1":   is_oracle_top1,   # 1 = strategy matched oracle at step 1
                    })
                    rows_auc.append({
                        "trial":         trial_idx,
                        "pathway":       pw_name,
                        "strategy":      strategy,
                        "step":          step + 1,
                        "n_metabolites": len(kept),
                        "auc":           compute_auc(sc_tr, sk_tr, sc_te, sk_te, kept),
                    })

        # Pathway ranking stability
        common = [p for p in full_panel_ranking if p in partial_uk]
        if len(common) >= 3:
            tau, _ = kendalltau([full_panel_uk[p] for p in common],
                                [partial_uk[p]    for p in common])
            rows_path_tau.append({
                "trial":        trial_idx,
                "n_pathways":   len(common),
                "kendall_tau":  float(tau) if np.isfinite(tau) else np.nan,
                "masking_rate": INITIAL_MASK_RATE,
            })

        if trial_idx % 10 == 0:
            print(f"  Trial {trial_idx:>2}/{N_TRIALS}  "
                  f"delta_uk rows: {len(rows_delta_uk)}")

    # ── Save ──────────────────────────────────────────────────────────────────
    pd.DataFrame(rows_delta_uk).to_csv(OUT_DIR/"simulation_delta_uk.csv",   index=False)
    pd.DataFrame(rows_auc     ).to_csv(OUT_DIR/"simulation_auc.csv",        index=False)
    pd.DataFrame(rows_path_tau).to_csv(OUT_DIR/"simulation_pathway_tau.csv", index=False)

    print(f"\n[done]")
    print(f"  ΔU_k rows:      {len(rows_delta_uk)}")
    print(f"  AUC rows:       {len(rows_auc)}")
    print(f"  Pathway τ rows: {len(rows_path_tau)}")

    # ── Summary ───────────────────────────────────────────────────────────────
    if rows_delta_uk:
        df_d  = pd.DataFrame(rows_delta_uk)
        step1 = df_d[df_d["step"] == 1]

        print("\n── Step-1 summary (first reveal only) ───────────────────────────")
        print(f"Trials with ≥2 hidden nodes: {len(step1)//len(STRATEGIES)}")
        print(f"\nMedian ΔU_k at step 1 (higher = more ambiguity removed):")
        print(step1.groupby("strategy")["delta_uk"].median().round(4)
              .sort_values(ascending=False).to_string())

        print(f"\nFraction with positive ΔU_k (reveal reduced ambiguity):")
        print(step1.groupby("strategy").apply(
            lambda x: (x["delta_uk"] > 0).mean()).round(3).to_string())

        print(f"\nOracle top-1 recovery at step 1 (= surrogate fidelity in simulation):")
        print(step1.groupby("strategy")["oracle_top1"].mean().round(3)
              .sort_values(ascending=False).to_string())

    if rows_auc:
        df_a = pd.DataFrame(rows_auc)
        print("\n── Mean AUC by strategy and step ────────────────────────────────")
        sel  = df_a[df_a["step"].isin([1, 3, 5])]
        print(sel.groupby(["strategy","step"])["auc"].mean().round(4)
              .unstack("step").to_string())

    if rows_path_tau:
        df_t = pd.DataFrame(rows_path_tau)
        print(f"\n── Pathway ranking stability (Kendall τ) ────────────────────────")
        print(f"Mean τ:   {df_t['kendall_tau'].mean():.4f}")
        print(f"Std  τ:   {df_t['kendall_tau'].std():.4f}")
        print(f"(τ = 1.0: partial-panel ranking perfectly matches full-panel)")


if __name__ == "__main__":
    main()
