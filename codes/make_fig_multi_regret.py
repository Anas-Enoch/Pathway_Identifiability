import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
CSV = Path(__file__).resolve().parent / "multi_pathway_results.csv"

def bootstrap_ci_mean(x, B=1000, alpha=0.05, seed=0):
    """95% bootstrap CI for the mean."""
    rng = np.random.default_rng(seed)
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return np.nan, np.nan, np.nan

    mu = float(np.mean(x))
    boot = []
    n = len(x)
    for _ in range(B):
        samp = rng.choice(x, size=n, replace=True)
        boot.append(np.mean(samp))
    boot = np.sort(np.array(boot))
    lo = float(np.quantile(boot, alpha / 2))
    hi = float(np.quantile(boot, 1 - alpha / 2))
    return mu, lo, hi


def main():
    # --------
    # Input
    # --------
    CSV = Path("results/tables/multi_pathway_results.csv")
    df = pd.read_csv(CSV)

    # Ensure correct types
    df["rho"] = df["rho"].astype(float)
    df["nregret"] = df["nregret"].astype(float)

    # Methods order (put yours first)
    preferred_order = ["JL-FGW-U", "FBA", "Imputation-RF", "Max-Degree", "Random"]
    methods = [m for m in preferred_order if m in df["method"].unique()]
    # include any others found
    for m in sorted(df["method"].unique()):
        if m not in methods:
            methods.append(m)

    rhos = sorted(df["rho"].unique())

    # Style (no explicit colors)
    markers = ["o", "s", "^", "D", "x", "v", "P", "*"]
    linestyles = ["-","--","--","--","--","--","--","--"]

    plt.figure(figsize=(6.2, 4.2))

    for idx, method in enumerate(methods):
        means, lows, highs = [], [], []
        for rho in rhos:
            x = df.loc[(df["method"] == method) & (df["rho"] == rho), "nregret"].values
            mu, lo, hi = bootstrap_ci_mean(x, B=1000, seed=idx + int(100*rho))
            means.append(mu); lows.append(lo); highs.append(hi)

        means = np.array(means)
        lows = np.array(lows)
        highs = np.array(highs)

        plt.plot(
            rhos, means,
            marker=markers[idx % len(markers)],
            linestyle=linestyles[idx % len(linestyles)],
            linewidth=2 if method == "JL-FGW-U" else 1.6,
            label=method
        )
        # CI band
        plt.fill_between(rhos, lows, highs, alpha=0.15)

    plt.xlabel("Masking rate $\\rho$")
    plt.ylabel("Mean normalized regret")
    plt.ylim(0.0, 1.0)
    plt.title("Multi-pathway benchmark: regret vs masking rate")

    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig("Fig_multi_pathway_regret.pdf")
    plt.close()

    print("Saved Fig_multi_pathway_regret.pdf")


if __name__ == "__main__":
    main()
