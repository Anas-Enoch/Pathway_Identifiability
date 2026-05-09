#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_msea_mummichog_benchmark.py  (FIXED — fast oracle)

Key changes vs. the lagging version:
  1. SINK_MAX_ITER 5000 → 300   (Sinkhorn converges well before 5000)
  2. TRIALS_PER_RHO  50 → 20   (20 is sufficient for stable estimates)
  3. MAX_NODES       inf → 25   (skip pathways with >25 metabolites)
  4. Progress printed at START of each pathway, not end
"""
from __future__ import annotations
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
from scipy.stats import kendalltau, ttest_ind
from scipy.spatial.distance import cdist
import ot

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

DATASETS = {
    "ST000356": DATA_DIR / "processed_metabolite_matrix_ST000356.csv",
    "ST003390": DATA_DIR / "ST003390_processed.csv",
    "ST003506": DATA_DIR / "ST003506_serum_processed.csv",
}
PATHWAY_MAP    = DATA_DIR / "core_pathway_mapping.csv"
OUT_FILE       = OUT_DIR  / "msea_mummichog_results.csv"

MASK_RATES     = [0.1, 0.2, 0.3, 0.4, 0.5]
TRIALS_PER_RHO = 20          # was 50
MIN_MET        = 3
MAX_NODES      = 25          # NEW: skip large pathways
METHOD_DEFAULT = "none"
ALPHA_FGW      = 0.5
SINK_REG       = 0.5
SINK_MAX_ITER  = 300         # was 5000 — this was the main bottleneck
EPS            = 1e-12
PVAL_THRESH    = 0.05

warnings.filterwarnings("ignore", message="Sinkhorn did not converge")

def fix_st000356(df):
    cols = list(df.columns)
    if "Metabolite_name" not in cols: return df
    idx = cols.index("Metabolite_name"); s_c = cols[idx+1:]
    cr  = df.iloc[0][s_c]
    def lbl(x):
        x=str(x).strip().lower()
        if any(k in x for k in ["control","healthy","normal"]): return "control"
        if any(k in x for k in ["cancer","tumor","tumour"]):    return "cancer"
        return "unknown"
    conds = cr.apply(lbl)
    df2   = df.iloc[1:].copy()[["Metabolite_name"]+s_c].set_index("Metabolite_name")[s_c]
    dft   = df2.T.reset_index(drop=True).apply(pd.to_numeric,errors="coerce")
    dft["condition"] = conds.values
    return dft[dft["condition"]!="unknown"].reset_index(drop=True)

def detect_cond_col(df):
    for c in ["condition","group","label","phenotype","class","diagnosis","status"]:
        for col in df.columns:
            if col.lower()==c: return col
    return None

def choose_case_ctrl(groups):
    groups=list(groups); lm={str(g).lower():g for g in groups}
    for a in ["control","ctrl","normal","healthy"]:
        if a in lm:
            ctrl=lm[a]; case=[g for g in groups if g!=ctrl][0]; return case,ctrl
    return groups[0],groups[1]

def detect_met_cols(df):
    skip={"sample_id","sampleid","condition","group","label",
          "class","diagnosis","status","phenotype","refmet_name"}
    out=[]
    for c in df.columns:
        if c.lower() in skip: continue
        try: pd.to_numeric(df[c],errors="raise"); out.append(c)
        except: pass
    return out

def load_pathway_mapping(path):
    df=pd.read_csv(path); df.columns=[c.strip().lower() for c in df.columns]
    if "metabolite_name" in df.columns: df=df.rename(columns={"metabolite_name":"metabolite"})
    if "pathway_name"    in df.columns: df=df.rename(columns={"pathway_name":"pathway"})
    assert {"metabolite","pathway"}.issubset(df.columns)
    pw=df.groupby("pathway")["metabolite"].apply(lambda x:sorted(set(map(str,x)))).to_dict()
    return pw,{p:[] for p in pw}

def sanitize_square(M):
    M=np.asarray(M,dtype=np.float64)
    M=np.nan_to_num(M,nan=0.,posinf=0.,neginf=0.)
    M=0.5*(M+M.T); np.fill_diagonal(M,0.); return M

def corr_dist(X):
    X=np.nan_to_num(np.asarray(X,dtype=np.float64),nan=0.,posinf=0.,neginf=0.)
    if X.shape[0]<2 or X.shape[1]<2: return np.zeros((X.shape[0],X.shape[0]))
    R=np.clip(np.corrcoef(X),-1.,1.)
    R=np.nan_to_num(R,nan=0.,posinf=0.,neginf=0.)
    return sanitize_square(np.sqrt(np.maximum(2.*(1.-R),0.)))

def node_features(X):
    """(n_samples, n_nodes) → (n_nodes, 2)  [mean, std per node]"""
    X=np.asarray(X,dtype=float)
    mu=np.nanmean(X,axis=0); sd=np.nanstd(X,axis=0)
    return np.column_stack([mu,sd])

def fgw_cross(X_s,X_t,C_s,C_t):
    def _n(M):
        M=sanitize_square(M); s=np.mean(M[M>0]) if np.any(M>0) else 1.
        return M/(s if s>EPS else 1.)
    X_s=np.nan_to_num(np.asarray(X_s,dtype=np.float64),nan=0.,posinf=0.,neginf=0.)
    X_t=np.nan_to_num(np.asarray(X_t,dtype=np.float64),nan=0.,posinf=0.,neginf=0.)
    p=np.ones(X_s.shape[0])/max(X_s.shape[0],1)
    q=np.ones(X_t.shape[0])/max(X_t.shape[0],1)
    M=cdist(X_s,X_t,metric="euclidean")
    d=X_s.shape[1] if X_s.ndim==2 else 1
    M=np.nan_to_num(M/(np.sqrt(d)+EPS),nan=0.,posinf=0.,neginf=0.)
    try:
        T,log=ot.gromov.entropic_fused_gromov_wasserstein(
            M,_n(C_s),_n(C_t),p,q,loss_fun="square_loss",
            epsilon=SINK_REG,alpha=ALPHA_FGW,log=True,numItermax=SINK_MAX_ITER)
    except TypeError:
        T,log=ot.gromov.entropic_fused_gromov_wasserstein(
            M,_n(C_s),_n(C_t),p,q,loss_fun="square_loss",
            epsilon=SINK_REG,alpha=ALPHA_FGW,log=True,max_iter=SINK_MAX_ITER)
    except Exception: return np.nan,None
    T=np.nan_to_num(np.asarray(T,dtype=np.float64),nan=0.,posinf=0.,neginf=0.)
    if T.sum()<=0 or not np.all(np.isfinite(T)): return np.nan,None
    dval=float(log.get("fgw_dist",np.nan)) if isinstance(log,dict) else np.nan
    if not np.isfinite(dval): dval=float(np.sum(T*M))
    return dval,T

def spectral_gap(C):
    C=sanitize_square(C); pos=C[C>0]
    if len(pos)==0: return 0.
    sigma=float(np.median(pos))
    if not np.isfinite(sigma) or sigma<EPS: return 0.
    A=np.exp(-(C**2)/(2.*sigma**2+EPS))
    A=np.nan_to_num(A,nan=0.,posinf=0.,neginf=0.)
    A=0.5*(A+A.T); np.fill_diagonal(A,0.); A[A<0]=0.
    d=A.sum(axis=1); L=np.diag(d)-A
    L=np.nan_to_num(L,nan=0.,posinf=0.,neginf=0.); L=0.5*(L+L.T)
    try: ev=np.sort(np.linalg.eigvalsh(L))
    except: return 0.
    if len(ev)<2: return 0.
    denom=float(ev[-1])+EPS
    return max(float(ev[1]/denom),0.) if np.isfinite(denom) and abs(denom)>EPS else 0.

def transport_entropy(T):
    T=np.maximum(np.asarray(T,dtype=np.float64),EPS)
    T=T/(T.sum()+EPS)
    return float(-np.sum(T*np.log(T+EPS)))

def compute_Uk(X_s,X_t,C_s,C_t):
    _,T=fgw_cross(X_s,X_t,C_s,C_t)
    if T is None: raise ValueError("FGW failed")
    return float(transport_entropy(T)+spectral_gap(C_s))

def mask_nodes(n_nodes,rho,rng):
    n_hide=max(1,min(int(round(rho*n_nodes)),n_nodes-2))
    perm=rng.permutation(n_nodes).tolist()
    return sorted(perm[n_hide:]),sorted(perm[:n_hide])

def rank_metrics_from_scores(deltas,pred_scores,nodes):
    if not deltas or not pred_scores: return None
    common=[hi for hi in deltas if hi in pred_scores]
    if not common: return None
    osorted=sorted(common,key=lambda h:deltas[h],reverse=True)
    psorted=sorted(common,key=lambda h:pred_scores[h],reverse=True)
    hi_star,hi_hat=osorted[0],psorted[0]
    ds=float(deltas[hi_star]); dh=float(deltas[hi_hat])
    regret=ds-dh; nregret=regret/(abs(ds)+EPS)
    top1=int(hi_hat==hi_star); top3=len(set(osorted[:3])&set(psorted[:3]))/3.
    tau,_=kendalltau([deltas[h] for h in common],[pred_scores[h] for h in common])
    return {
        "delta_star":round(ds,8),"delta_hat":round(dh,8),
        "regret":round(regret,8),"nregret":round(nregret,8),
        "top1":top1,"top3":round(top3,4),
        "rank_tau":round(float(tau),6) if np.isfinite(tau) else np.nan,
        "m_star":nodes[hi_star],"m_hat":nodes[hi_hat],
    }

def msea_proxy_scores(sub_case,sub_ctrl,kept_idx,hidden_idx):
    if not hidden_idx: return {}
    t_stats=np.zeros(sub_case.shape[1]); p_vals=np.ones(sub_case.shape[1])
    for j in range(sub_case.shape[1]):
        a=sub_case[:,j][np.isfinite(sub_case[:,j])]
        b=sub_ctrl[:,j][np.isfinite(sub_ctrl[:,j])]
        if len(a)<2 or len(b)<2: continue
        try:
            t,p=ttest_ind(a,b,equal_var=False)
            if np.isfinite(t): t_stats[j]=abs(float(t))
            if np.isfinite(p): p_vals[j]=float(p)
        except: pass
    n_hits=sum(1 for k in kept_idx if p_vals[k]<PVAL_THRESH)
    hit_frac=n_hits/max(len(kept_idx),1)
    return {hi:float(t_stats[hi]*(1.0+hit_frac)) for hi in hidden_idx}

def mummichog_proxy_scores(sub_case,sub_ctrl,C_full,hidden_idx):
    if not hidden_idx: return {}
    mu_c=np.nanmean(sub_case[:,hidden_idx],axis=0)
    mu_t=np.nanmean(sub_ctrl[:,hidden_idx],axis=0)
    diff=np.nan_to_num(np.abs(mu_c-mu_t),nan=0.,posinf=0.,neginf=0.)
    C=np.nan_to_num(np.asarray(C_full,dtype=np.float64),nan=0.,posinf=0.,neginf=0.)
    deg=np.sum(C>0,axis=1).astype(np.float64)
    return {hi:float(diff[j]*(1.0+deg[hi])) for j,hi in enumerate(hidden_idx)}

def run_trial_msea(sub_case,sub_ctrl,nodes,C_full,rho,rng):
    kept_idx,hidden_idx=mask_nodes(len(nodes),rho,rng)
    if len(hidden_idx)<2 or len(kept_idx)<2: return []
    X_s=node_features(sub_case[:,kept_idx])
    X_t=node_features(sub_ctrl[:,kept_idx])
    C_mask=C_full[np.ix_(kept_idx,kept_idx)]
    try:
        Uk_base=compute_Uk(X_s,X_t,C_mask,C_mask)
    except: return []
    deltas={}
    for hi in hidden_idx:
        aug=kept_idx+[hi]
        try:
            Uk_aug=compute_Uk(node_features(sub_case[:,aug]),
                               node_features(sub_ctrl[:,aug]),
                               C_full[np.ix_(aug,aug)],
                               C_full[np.ix_(aug,aug)])
        except: continue
        delta=Uk_base-Uk_aug
        if np.isfinite(delta): deltas[hi]=float(delta)
    if len(deltas)<2: return []
    valid=sorted(deltas.keys())
    rows=[]
    for name,scores in [("msea_proxy",    msea_proxy_scores(sub_case,sub_ctrl,kept_idx,valid)),
                         ("mummichog_proxy",mummichog_proxy_scores(sub_case,sub_ctrl,C_full,valid))]:
        m=rank_metrics_from_scores(deltas,scores,nodes)
        if m is None: continue
        rows.append({"rho":rho,"n_nodes":len(nodes),"n_hidden":len(valid),
                     "operator_method":METHOD_DEFAULT,"predictor_method":name,**m})
    return rows

def main():
    pw_members,_=load_pathway_mapping(PATHWAY_MAP)
    rows=[]
    for ds_name,ds_path in DATASETS.items():
        if not ds_path.exists():
            print(f"[WARN] {ds_path} not found — skipping"); continue
        df=pd.read_csv(ds_path)
        if ds_name=="ST000356": df=fix_st000356(df)
        cond_col=detect_cond_col(df)
        if cond_col is None:
            print(f"[WARN] {ds_name}: no condition column — skipping"); continue
        groups=list(df[cond_col].dropna().unique())
        if len(groups)!=2:
            print(f"[WARN] {ds_name}: {len(groups)} groups — skipping"); continue
        case_label,ctrl_label=choose_case_ctrl(groups)
        df_case=df[df[cond_col]==case_label].copy()
        df_ctrl=df[df[cond_col]==ctrl_label].copy()
        print(f"[{ds_name}] case='{case_label}' (n={len(df_case)}), ctrl='{ctrl_label}' (n={len(df_ctrl)})")
        met_set=set(detect_met_cols(df))
        pw_count=0
        for pw_name,members in pw_members.items():
            available=[m for m in members if m in met_set]
            if len(available)<MIN_MET: continue
            if len(available)>MAX_NODES:
                print(f"  [skip] {pw_name}: {len(available)} nodes > {MAX_NODES}"); continue
            sc=df_case[available].apply(pd.to_numeric,errors="coerce").values.astype(np.float64)
            sk=df_ctrl[available].apply(pd.to_numeric,errors="coerce").values.astype(np.float64)
            good=[j for j in range(len(available))
                  if not(np.all(~np.isfinite(sc[:,j])) or np.all(~np.isfinite(sk[:,j])))]
            if len(good)<MIN_MET: continue
            available=[available[j] for j in good]
            sc=sc[:,good]; sk=sk[:,good]
            C_full=corr_dist(np.vstack([sc,sk]).T)
            # Print at start so you see immediate activity
            print(f"  → {pw_name} ({len(available)} nodes) ...",end="",flush=True)
            n_ok=0
            for rho in MASK_RATES:
                for t in range(TRIALS_PER_RHO):
                    rng=np.random.default_rng(1000+t+int(100*rho))
                    for r in run_trial_msea(sc,sk,available,C_full,rho,rng):
                        r["dataset"]=ds_name; r["pathway"]=pw_name
                        rows.append(r); n_ok+=1
            print(f" {n_ok} OK"); pw_count+=1
        print(f"[{ds_name}] done — {pw_count} pathways")
    if not rows: print("[WARN] no rows"); return
    out=pd.DataFrame(rows)
    col_order=["dataset","pathway","rho","operator_method","predictor_method",
               "n_nodes","n_hidden","delta_star","delta_hat","regret","nregret",
               "top1","top3","rank_tau","m_star","m_hat"]
    out=out[[c for c in col_order if c in out.columns]]
    out.to_csv(OUT_FILE,index=False)
    print(f"\n[done] {len(out)} rows → {OUT_FILE}")
    hard=out[out["n_hidden"]>=2]
    if len(hard):
        print("\nHard-subset summary (n_hidden ≥ 2):")
        print(hard.groupby("predictor_method")[["regret","nregret","top1","top3","rank_tau"]]
              .mean().round(4).sort_values("rank_tau",ascending=False).to_string())

if __name__=="__main__":
    main()
