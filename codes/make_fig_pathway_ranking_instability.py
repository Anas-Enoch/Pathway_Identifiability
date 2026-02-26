#!/usr/bin/env python3
"""
Generate Fig_pathway_ranking_instability.pdf

Inputs:
  1) ranking_stability_by_method.csv
     Required columns:
       - method : string (e.g., "Mean", "Enrichment", "GSVA")
       - tau    : float  (Kendall tau per replicate)
     Optional:
       - replicate : int or string (only used for sanity checks)

  2) pathway_instability_vs_Uk.csv
     Required columns:
       - pathway : string (pathway identifier/name)
       - Uk      : float  (predicted underdetermination)
       - Instab  : float  (E[|Δrank|] per pathway)
     Optional:
       - sgi_bin : string/category for coloring (e.g., "low","mid","high")
       - label   : 0/1 flag to annotate selected pathways
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import spearmanr


DEFAULT_METHOD_ORDER = ["Mean", "Enrichment", "GSVA"]


def _assert_cols(df: pd.DataFrame, required: list[str], name: str) -> None:
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"{name} is missing required columns: {missing}. Found: {list(df.columns)}")


def make_figure(
    stability_csv: Path,
    pathway_csv: Path,
    out_pdf: Path,
    method_order: list[str] | None = None,
    rho: float = 0.3,
    R: int | None = None,
    random_seed: int = 0,
) -> None:
    rng = np.random.default_rng(random_seed)

    # -------------------------
    # Load data
    # -------------------------
    stab = pd.read_csv(stability_csv)
    _assert_cols(stab, ["method", "tau"], "ranking_stability_by_method.csv")

    pw = pd.read_csv(pathway_csv)
    _assert_cols(pw, ["pathway", "Uk", "Instab"], "pathway_instability_vs_Uk.csv")

    # Clean
    stab = stab.copy()
    stab["method"] = stab["method"].astype(str)
    stab["tau"] = pd.to_numeric(stab["tau"], errors="coerce")
    stab = stab.dropna(subset=["tau"])

    pw = pw.copy()
    pw["Uk"] = pd.to_numeric(pw["Uk"], errors="coerce")
    pw["Instab"] = pd.to_numeric(pw["Instab"], errors="coerce")
    pw = pw.dropna(subset=["Uk", "Instab"])

    # Method order
    if method_order is None:
        # Use default if present, else alphabetical
        present = list(pd.unique(stab["method"]))
        if all(m in present for m in DEFAULT_METHOD_ORDER):
            method_order = DEFAULT_METHOD_ORDER
        else:
            method_order = sorted(present)

    stab["method"] = pd.Categorical(stab["method"], categories=method_order, ordered=True)
    stab = stab.sort_values("method")

    # Spearman association
    rho_s, pval = spearmanr(pw["Uk"].values, pw["Instab"].values)

    # -------------------------
    # Plot
    # -------------------------
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2), constrained_layout=True)

    # ---- Panel A: boxplot + jitter
    ax = axes[0]

    grouped = [stab.loc[stab["method"] == m, "tau"].values for m in method_order]
    positions = np.arange(1, len(method_order) + 1)

    bp = ax.boxplot(
        grouped,
        positions=positions,
        widths=0.55,
        patch_artist=True,
        showfliers=False,  # we show replicate points ourselves
        medianprops=dict(linewidth=1.4),
        whiskerprops=dict(linewidth=1.0),
        capprops=dict(linewidth=1.0),
        boxprops=dict(linewidth=1.0),
    )
    # Do NOT set specific colors; keep default styling. (Journal-safe)
    # Matplotlib default will color patches; make them transparent to avoid color commitment.
    for box in bp["boxes"]:
        box.set_alpha(0.15)

    # Jittered replicate points
    for i, m in enumerate(method_order, start=1):
        y = stab.loc[stab["method"] == m, "tau"].values
        x = i + rng.uniform(-0.12, 0.12, size=len(y))
        ax.scatter(x, y, s=12, alpha=0.55, linewidths=0)

    ax.set_xticks(positions)
    ax.set_xticklabels(method_order)
    ax.set_ylim(0.0, 1.0)
    ax.set_ylabel("Ranking stability (Kendall τ)")
    ax.set_xlabel("Pathway activity estimator")
    ax.set_title("(A) Ranking stability under masking")

    # Small annotation
    text = f"ρ={rho:.1f}"
    if R is not None:
        text += f", R={R}"
    ax.text(
        0.02,
        0.98,
        text,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
    )

    # ---- Panel B: scatter Uk vs Instab + fit + annotation
    ax = axes[1]

    x = pw["Uk"].values
    y = pw["Instab"].values

    # Optional coloring by sgi_bin if present
    if "sgi_bin" in pw.columns:
        # Keep it minimal: scatter per category
        for cat, dfc in pw.groupby("sgi_bin"):
            ax.scatter(dfc["Uk"].values, dfc["Instab"].values, s=18, alpha=0.70, label=str(cat))
        ax.legend(frameon=False, fontsize=8, loc="best")
    else:
        ax.scatter(x, y, s=18, alpha=0.70)

    # Simple linear fit (robust enough, easy to read)
    if len(x) >= 2 and np.std(x) > 0:
        coef = np.polyfit(x, y, deg=1)
        xs = np.linspace(np.min(x), np.max(x), 200)
        ys = coef[0] * xs + coef[1]
        ax.plot(xs, ys, linewidth=1.4)

    ax.set_xlabel("Predicted underdetermination (Uₖ)")
    ax.set_ylabel("Ranking instability (E[|Δ rank|])")
    ax.set_title("(B) Uₖ predicts instability")

    # Spearman annotation
    ax.text(
        0.02,
        0.98,
        f"Spearman ρₛ={rho_s:.2f}\n$p$={pval:.2e}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
    )

    # Optional pathway labels (if column 'label' exists OR if top extremes)
    if "label" in pw.columns:
        to_label = pw[pw["label"].astype(int) == 1]
    else:
        # auto-label a few extreme points (top 3 by Uk and Instab)
        to_label = pd.concat(
            [
                pw.nlargest(3, "Uk"),
                pw.nlargest(3, "Instab"),
            ]
        ).drop_duplicates(subset=["pathway"])

    for _, row in to_label.iterrows():
        ax.annotate(
            str(row["pathway"]),
            (row["Uk"], row["Instab"]),
            textcoords="offset points",
            xytext=(5, 4),
            fontsize=7,
            alpha=0.9,
        )

    # Save
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--stability_csv", type=Path, default=Path("ranking_stability_by_method.csv"))
    p.add_argument("--pathway_csv", type=Path, default=Path("pathway_instability_vs_Uk.csv"))
    p.add_argument("--out", type=Path, default=Path("Fig_pathway_ranking_instability.pdf"))
    p.add_argument("--rho", type=float, default=0.3)
    p.add_argument("--R", type=int, default=None)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--method_order", type=str, default="Mean,Enrichment,GSVA",
                   help="Comma-separated, e.g. 'Mean,Enrichment,GSVA'")
    args = p.parse_args()

    method_order = [m.strip() for m in args.method_order.split(",") if m.strip()]
    make_figure(
        stability_csv=args.stability_csv,
        pathway_csv=args.pathway_csv,
        out_pdf=args.out,
        method_order=method_order,
        rho=args.rho,
        R=args.R,
        random_seed=args.seed,
    )


if __name__ == "__main__":
    main()
