"""
Multi-pathway benchmark loop (real controlled masking on real metabolomics data)

Goal
----
For many KEGG pathways, run controlled masking on metabolites that are truly observed,
recommend next-best measurement, and score via oracle regret:
    ΔU(m) = U(masked) - U(masked + m revealed)
    regret = ΔU(m*) - ΔU(m_hat)

This skeleton is deliberately modular:
- You plug in your existing graph construction, JL-FGW, U_k computation, and impact ranking.
- You can run baselines in the same loop.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Tuple, Set, Optional, Callable
import numpy as np
import pandas as pd
import time


# =========================
# Data structures
# =========================

Condition = str  # e.g., "tumor", "normal"
Metabolite = str
PathwayID = str


@dataclass
class PathwaySpec:
    """Defines a pathway by its metabolite nodes and (optionally) reaction edges, etc."""
    pathway_id: PathwayID
    metabolites: List[Metabolite]
    # Optional: edges, stoichiometry, compartments, enzymes...
    # edges: List[Tuple[Metabolite, Metabolite]] = None


@dataclass
class ConditionStats:
    """Per-condition distributional embeddings for metabolites: mean/var on log-normalized scale."""
    mu: Dict[Metabolite, float]
    var: Dict[Metabolite, float]
    observed: Set[Metabolite]  # metabolites considered observed in this condition (<=50% missing rule)


@dataclass
class TrialResult:
    pathway_id: PathwayID
    rho: float
    trial: int
    method: str
    m_hat: Optional[Metabolite]
    m_star: Optional[Metabolite]
    delta_hat: float
    delta_star: float
    regret: float
    nregret: float
    runtime_s: float


# =========================
# Preprocessing helpers
# =========================

def compute_condition_stats(
    X: pd.DataFrame,
    meta: pd.DataFrame,
    condition_col: str,
    condition_value: str,
    epsilon: float = 1e-6,
    missing_threshold: float = 0.5,
    median_center: bool = True,
) -> ConditionStats:
    """
    X: rows=metabolites, cols=samples (raw intensities or concentrations)
    meta: rows=samples, contains condition labels
    Returns: mean/var per metabolite for samples in condition_value.
    """
    samples = meta.index[meta[condition_col] == condition_value].tolist()
    Xc = X[samples].copy()

    # Log transform
    Xc = np.log(Xc + epsilon)

    # Optional sample-wise median centering (on log scale)
    if median_center:
        med = Xc.median(axis=0, skipna=True)
        Xc = Xc.sub(med, axis=1)

    # Missingness handling: define observed set; impute only for those <= threshold
    observed = set()
    mu, var = {}, {}

    for m in Xc.index:
        col = Xc.loc[m]
        frac_missing = col.isna().mean()

        if frac_missing > missing_threshold:
            # treat as unobserved/latent in this condition
            continue

        # impute half min (within-condition)
        if col.isna().any():
            min_obs = col.min(skipna=True)
            col = col.fillna(0.5 * min_obs)

        observed.add(m)
        mu[m] = float(col.mean())
        var[m] = float(col.var(ddof=1)) if len(col) > 1 else 0.0

    return ConditionStats(mu=mu, var=var, observed=observed)


# =========================
# Your core pipeline hooks
# =========================

def build_condition_graph(
    pathway: PathwaySpec,
    stats: ConditionStats,
    latent_prior_var: float = 1.0,
) -> dict:
    """
    Build your condition-aware pathway graph with distributional node embeddings.
    Return whatever graph object your JL-FGW code expects.

    Minimal representation example:
        graph = {
            "nodes": pathway.metabolites,
            "mu": {m: mu or prior},
            "var": {m: var or prior},
            "observed": set(...)
            # edges/topology, etc.
        }
    """
    graph = {
        "nodes": list(pathway.metabolites),
        "mu": {},
        "var": {},
        "observed": set(stats.observed),
        # TODO: add edges/topology
        # "edges": pathway.edges,
    }

    for m in pathway.metabolites:
        if m in stats.mu:
            graph["mu"][m] = stats.mu[m]
            graph["var"][m] = stats.var.get(m, latent_prior_var)
        else:
            # latent node: prior mean 0, large variance
            graph["mu"][m] = 0.0
            graph["var"][m] = latent_prior_var

    return graph


def compute_Uk(
    graph_c1: dict,
    graph_c2: dict,
    weights: Tuple[float, float, float],
) -> float:
    """
    Compute U_k for a pathway across two conditions using your method:
        U = α H(T) + β Var(A) + γ SGI

    Replace this stub with your real code.
    """
    alpha, beta, gamma = weights

    # TODO:
    # - run JL-FGW alignment, get transport plan T
    # - compute entropy H(T)
    # - compute instability Var(A) via perturbations
    # - compute SGI_k
    # return alpha*H + beta*VarA + gamma*SGI

    raise NotImplementedError("Plug in your real U_k computation here.")


def recommend_metabolite_your_method(
    graph_c1: dict,
    graph_c2: dict,
    unobserved_candidates: List[Metabolite],
    weights: Tuple[float, float, float],
) -> Optional[Metabolite]:
    """
    Return m_hat using your locked estimator (gradient-based impact).
    """
    # TODO: compute impact scores for each candidate and take argmax
    # e.g., impact[m] = ||∂U/∂mu_m|| * var_m
    raise NotImplementedError("Plug in your measurement impact ranking here.")


# =========================
# Baselines (examples)
# =========================

def baseline_max_degree(pathway: PathwaySpec, candidates: List[Metabolite], graph: dict) -> Optional[Metabolite]:
    """
    Requires edges/topology. If you don't have edges yet, skip or substitute with
    max-appearance / max-neighbors from KEGG adjacency.
    """
    # TODO: implement degree from graph["edges"]
    return candidates[0] if candidates else None


def baseline_random(candidates: List[Metabolite], rng: np.random.Generator) -> Optional[Metabolite]:
    return rng.choice(candidates) if candidates else None


# =========================
# Oracle ΔU(m): reveal one metabolite and recompute U
# =========================

def reveal_metabolite(graph: dict, m: Metabolite, true_mu: float, true_var: float = 0.0) -> dict:
    """
    Return a copy of graph where metabolite m is treated as observed with true value.
    """
    g2 = {
        "nodes": graph["nodes"],
        "mu": dict(graph["mu"]),
        "var": dict(graph["var"]),
        "observed": set(graph["observed"]),
        # Copy other fields if you add them (edges, etc.)
        **{k: v for k, v in graph.items() if k not in {"nodes", "mu", "var", "observed"}}
    }
    g2["mu"][m] = float(true_mu)
    g2["var"][m] = float(true_var)
    g2["observed"].add(m)
    return g2


def deltaU_oracle(
    pathway: PathwaySpec,
    base_graph_c1: dict,
    base_graph_c2: dict,
    masked_set: Set[Metabolite],
    true_stats_c1: ConditionStats,
    true_stats_c2: ConditionStats,
    weights: Tuple[float, float, float],
    epsilon: float = 1e-9,
) -> Tuple[Optional[Metabolite], float, Dict[Metabolite, float]]:
    """
    Compute oracle m* and ΔU(m) for each masked metabolite m by revealing it and recomputing U.
    This does NOT enumerate latent completions: it only tries each candidate metabolite.
    """
    U_base = compute_Uk(base_graph_c1, base_graph_c2, weights)

    delta_map: Dict[Metabolite, float] = {}
    best_m, best_delta = None, -np.inf

    for m in masked_set:
        # Reveal using true measured values (condition-specific means/vars)
        # Here we reveal in both conditions if present; you can choose reveal only in one.
        g1 = reveal_metabolite(
            base_graph_c1, m,
            true_mu=true_stats_c1.mu.get(m, base_graph_c1["mu"][m]),
            true_var=true_stats_c1.var.get(m, 0.0),
        )
        g2 = reveal_metabolite(
            base_graph_c2, m,
            true_mu=true_stats_c2.mu.get(m, base_graph_c2["mu"][m]),
            true_var=true_stats_c2.var.get(m, 0.0),
        )

        U_reveal = compute_Uk(g1, g2, weights)
        dU = float(U_base - U_reveal)
        delta_map[m] = dU

        if dU > best_delta + epsilon:
            best_delta = dU
            best_m = m

    return best_m, float(best_delta if best_m is not None else 0.0), delta_map


# =========================
# Main benchmark loop
# =========================

def run_multi_pathway_benchmark(
    pathways: List[PathwaySpec],
    X: pd.DataFrame,
    meta: pd.DataFrame,
    condition_col: str,
    c1: Condition = "tumor",
    c2: Condition = "normal",
    rhos: List[float] = [0.1, 0.2, 0.3, 0.4, 0.5],
    trials_per_rho: int = 25,
    weights: Tuple[float, float, float] = (1.0, 1.0, 1.0),
    min_obs_overlap: int = 8,
    min_pathway_size: int = 12,
    rng_seed: int = 0,
) -> pd.DataFrame:
    rng = np.random.default_rng(rng_seed)
    results: List[TrialResult] = []

    # Compute true condition stats (used for controlled masking reveals)
    stats_c1 = compute_condition_stats(X, meta, condition_col, c1)
    stats_c2 = compute_condition_stats(X, meta, condition_col, c2)

    for pw in pathways:
        V = pw.metabolites
        if len(V) < min_pathway_size:
            continue

        overlap = list((set(V) & stats_c1.observed) & stats_c2.observed)
        if len(overlap) < min_obs_overlap:
            continue

        # Build full graphs using condition stats
        G1_full = build_condition_graph(pw, stats_c1)
        G2_full = build_condition_graph(pw, stats_c2)

        # Candidates we can mask are those truly observed (ground truth available)
        observable = overlap

        for rho in rhos:
            for t in range(trials_per_rho):
                # Controlled masking set
                k = max(1, int(round(rho * len(observable))))
                masked = set(rng.choice(observable, size=k, replace=False).tolist())

                # Create masked graphs: treat masked as latent (set prior)
                G1_mask = build_condition_graph(pw, stats_c1)
                G2_mask = build_condition_graph(pw, stats_c2)
                for m in masked:
                    # wipe observed status, replace with prior uncertainty
                    if m in G1_mask["observed"]:
                        G1_mask["observed"].remove(m)
                    if m in G2_mask["observed"]:
                        G2_mask["observed"].remove(m)
                    G1_mask["mu"][m] = 0.0
                    G2_mask["mu"][m] = 0.0
                    # keep high variance
                    G1_mask["var"][m] = max(G1_mask["var"][m], 1.0)
                    G2_mask["var"][m] = max(G2_mask["var"][m], 1.0)

                candidates = list(masked)  # evaluate recommendation among masked set only

                # ---- Your method recommendation ----
                t0 = time.perf_counter()
                try:
                    m_hat = recommend_metabolite_your_method(G1_mask, G2_mask, candidates, weights)
                except NotImplementedError:
                    m_hat = None
                t1 = time.perf_counter()

                # ---- Oracle + regret ----
                try:
                    m_star, delta_star, delta_map = deltaU_oracle(
                        pw, G1_mask, G2_mask, masked, stats_c1, stats_c2, weights
                    )
                except NotImplementedError:
                    m_star, delta_star, delta_map = None, 0.0, {}

                delta_hat = float(delta_map.get(m_hat, 0.0)) if (m_hat is not None) else 0.0
                regret = float(delta_star - delta_hat)
                nregret = float(regret / (delta_star + 1e-9)) if delta_star > 0 else 0.0

                results.append(
                    TrialResult(
                        pathway_id=pw.pathway_id,
                        rho=rho,
                        trial=t,
                        method="JL-FGW-U",
                        m_hat=m_hat,
                        m_star=m_star,
                        delta_hat=delta_hat,
                        delta_star=delta_star,
                        regret=regret,
                        nregret=nregret,
                        runtime_s=float(t1 - t0),
                    )
                )

                # ---- Baselines (optional) ----
                # Random
                m_hat_r = baseline_random(candidates, rng)
                delta_hat_r = float(delta_map.get(m_hat_r, 0.0)) if m_hat_r is not None else 0.0
                regret_r = float(delta_star - delta_hat_r)
                nregret_r = float(regret_r / (delta_star + 1e-9)) if delta_star > 0 else 0.0
                results.append(
                    TrialResult(
                        pathway_id=pw.pathway_id, rho=rho, trial=t, method="Random",
                        m_hat=m_hat_r, m_star=m_star,
                        delta_hat=delta_hat_r, delta_star=delta_star,
                        regret=regret_r, nregret=nregret_r, runtime_s=0.0
                    )
                )

                # Max-degree (if edges exist)
                # m_hat_deg = baseline_max_degree(pw, candidates, G1_mask)
                # ...

    df = pd.DataFrame([r.__dict__ for r in results])
    return df


# =========================
# Plot helpers (regret curves)
# =========================

def plot_regret_curves(df: pd.DataFrame, outfile: str = "Fig_multi_pathway_regret.pdf") -> None:
    import matplotlib.pyplot as plt

    # mean nRegret vs rho per method
    agg = df.groupby(["method", "rho"], as_index=False)["nregret"].mean()

    plt.figure(figsize=(6, 4))
    for method in sorted(agg["method"].unique()):
        sub = agg[agg["method"] == method].sort_values("rho")
        plt.plot(sub["rho"], sub["nregret"], marker="o", label=method)

    plt.xlabel("Masking rate $\\rho$")
    plt.ylabel("Mean normalized regret")
    plt.title("Multi-pathway benchmark: regret vs masking")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()
    print(f"Saved {outfile}")


# =========================
# Usage example (you replace I/O)
# =========================
if __name__ == "__main__":
    # TODO: load your metabolomics table and metadata
    # X = pd.read_csv("metabolomics_matrix.csv", index_col=0)  # metabolites x samples
    # meta = pd.read_csv("sample_metadata.csv", index_col=0)   # samples x covariates

    # TODO: build pathway list (KEGG parsing)
    # pathways = [...]
    # df = run_multi_pathway_benchmark(pathways, X, meta, condition_col="condition")
    # df.to_csv("multi_pathway_results.csv", index=False)
    # plot_regret_curves(df)

    print("Benchmark skeleton ready. Plug in X/meta/pathways + your U_k/recommender.")
