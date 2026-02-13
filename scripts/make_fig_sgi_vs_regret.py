import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def spearman_corr(x, y):
    """Compute Spearman correlation without scipy (rank correlation)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    rx = pd.Series(x).rank(method="average").to_numpy()
    ry = pd.Series(y).rank(method="average").to_numpy()
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    denom = np.sqrt((rx**2).sum()) * np.sqrt((ry**2).sum())
    return float((rx * ry).sum() / (denom + 1e-12))


def main():
    # run from repo root:
    #   python scripts/make_fig_sgi_vs_regret.py
    df = pd.read_csv("results/tables/multi_pathway_results.csv")

    # ---- choose method ----
    # If you later add baselines, pick the proposed method only
    if "method" in df.columns and "JL-FGW-U" in set(df["method"].unique()):
        df = df[df["method"] == "JL-FGW-U"].copy()

    # ---- choose rho for difficulty score ----
    rho_eval = 0.3
    dfe = df[df["rho"].astype(float) == rho_eval].copy()
    if dfe.empty:
        raise ValueError(f"No rows at rho={rho_eval}. Available: {sorted(df['rho'].unique())}")

    # ---- complexity proxy ----
    # Use median U_pred at rho=0.1 as pathway complexity (your current proxy).
    rho0 = 0.1
    base = df[df["rho"].astype(float) == rho0].copy()
    if base.empty:
        raise ValueError(f"No rows at rho={rho0}. Available: {sorted(df['rho'].unique())}")

    complexity = (
        base.groupby("pathway_id", as_index=False)["U_pred"]
        .median()
        .rename(columns={"U_pred": "complexity"})
    )

    # ---- pathway-level difficulty ----
    # One point per pathway: mean nregret over trials at rho_eval
    difficulty = (
        dfe.groupby("pathway_id", as_index=False)["nregret"]
        .mean()
        .rename(columns={"nregret": "difficulty"})
    )

    merged = complexity.merge(difficulty, on="pathway_id", how="inner").dropna()
    x = merged["complexity"].to_numpy(dtype=float)
    y = merged["difficulty"].to_numpy(dtype=float)

    # ---- correlations ----
    rho_s = spearman_corr(x, y)

    # ---- linear fit ----
    # y = a*x + b
    a, b = np.polyfit(x, y, deg=1)
    xs = np.linspace(x.min(), x.max(), 200)
    yhat = a * xs + b

    # ---- bootstrap CI band for the fitted line ----
    rng = np.random.default_rng(0)
    B = 1000
    boot = np.empty((B, len(xs)), dtype=float)

    n = len(x)
    for i in range(B):
        idx = rng.choice(np.arange(n), size=n, replace=True)
        ai, bi = np.polyfit(x[idx], y[idx], deg=1)
        boot[i] = ai * xs + bi

    lo = np.quantile(boot, 0.025, axis=0)
    hi = np.quantile(boot, 0.975, axis=0)

    # ---- plot ----
    plt.figure(figsize=(6.2, 4.2))

    plt.scatter(x, y, alpha=0.6, s=30)
    plt.plot(xs, yhat, linewidth=2)
    plt.fill_between(xs, lo, hi, alpha=0.15)

    plt.xlabel("Pathway complexity proxy (median $U_k$ at $\\rho=0.1$)")
    plt.ylabel(f"Mean normalized regret at $\\rho={rho_eval}$")
    plt.ylim(0.0, 1.0)
    plt.title("Pathway complexity predicts identifiability difficulty")

    # Put Spearman rho in the plot (top-left)
    plt.text(
        0.02, 0.98,
        f"Spearman $\\rho_s$ = {rho_s:.2f}",
        transform=plt.gca().transAxes,
        ha="left", va="top"
    )

    plt.tight_layout()
    out = "results/figures/Fig_complexity_vs_regret.pdf"
    plt.savefig(out)
    plt.close()
    print(f"Saved {out}")


if __name__ == "__main__":
    main()
