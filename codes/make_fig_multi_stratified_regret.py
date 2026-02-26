import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def bootstrap_ci_mean(x, B=2000, alpha=0.05, seed=0):
    """
    95% bootstrap CI for mean(x). x = 1D array.
    """
    rng = np.random.default_rng(seed)
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return np.nan, np.nan, np.nan
    mu = float(np.mean(x))
    n = len(x)
    boot = np.empty(B, dtype=float)
    for b in range(B):
        samp = rng.choice(x, size=n, replace=True)
        boot[b] = np.mean(samp)
    lo = float(np.quantile(boot, alpha / 2))
    hi = float(np.quantile(boot, 1 - alpha / 2))
    return mu, lo, hi


def main():
    # IMPORTANT: run from repo root:
    #   python scripts/make_fig_multi_stratified_regret.py
    df = pd.read_csv("results/tables/multi_pathway_results.csv")

    # sanity
    required = {"pathway_id", "rho", "trial", "U_pred", "nregret", "method"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns in CSV: {missing}")

    df["rho"] = df["rho"].astype(float)

    # ---- 1) Define pathway complexity proxy ----
    # Use median U_pred at low masking (rho=0.1) as a stable complexity score.
    rho0 = 0.1
    base = df[df["rho"] == rho0].copy()
    if base.empty:
        raise ValueError(f"No rows found for rho={rho0}. Available rhos: {sorted(df['rho'].unique())}")

    complexity = (
        base.groupby("pathway_id", as_index=False)["U_pred"]
        .median()
        .rename(columns={"U_pred": "complexity"})
    )

    # Assign tertiles
    q1, q2 = complexity["complexity"].quantile([1/3, 2/3]).values
    def tier(x):
        if x <= q1:
            return "Low complexity"
        elif x <= q2:
            return "Mid complexity"
        else:
            return "High complexity"

    complexity["tier"] = complexity["complexity"].apply(tier)

    # Merge tiers back
    df = df.merge(complexity[["pathway_id", "tier"]], on="pathway_id", how="left")

    # ---- 2) Aggregate to pathway-level means (avoid treating trials as independent pathways) ----
    # For each pathway, rho, method: average nregret over trials
    agg = (
        df.groupby(["pathway_id", "rho", "method", "tier"], as_index=False)["nregret"]
        .mean()
        .rename(columns={"nregret": "nregret_pw"})
    )

    rhos = sorted(agg["rho"].unique())
    tiers = ["Low complexity", "Mid complexity", "High complexity"]

    # If later you add baselines into CSV, this will automatically support multiple methods.
    methods = sorted(agg["method"].unique())
    # Plot only JL-FGW-U if present (cleaner)
    if "JL-FGW-U" in methods:
        methods = ["JL-FGW-U"]

    plt.figure(figsize=(6.4, 4.4))

    for ti, tier_name in enumerate(tiers):
        means, lows, highs = [], [], []
        for rho in rhos:
            sub = agg[(agg["tier"] == tier_name) & (agg["rho"] == rho) & (agg["method"].isin(methods))]
            x = sub["nregret_pw"].values  # one value per pathway
            mu, lo, hi = bootstrap_ci_mean(x, B=2000, seed=1000 + 10*ti + int(100*rho))
            means.append(mu); lows.append(lo); highs.append(hi)

        means = np.array(means); lows = np.array(lows); highs = np.array(highs)
        plt.plot(rhos, means, marker="o", linewidth=2, label=tier_name)
        plt.fill_between(rhos, lows, highs, alpha=0.15)

    plt.xlabel("Masking rate $\\rho$")
    plt.ylabel("Mean normalized regret")
    plt.ylim(0.0, 1.0)
    plt.title("Stratified benchmark: regret vs masking rate (by pathway complexity)")
    plt.legend(frameon=False)
    plt.tight_layout()

    out = "results/figures/Fig_multi_stratified_regret.pdf"
    plt.savefig(out)
    plt.close()
    print(f"Saved {out}")


if __name__ == "__main__":
    main()
