import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import binned_statistic, spearmanr

# -----------------------------
# Load results
# -----------------------------

df = pd.read_csv("results/tables/multi_pathway_results.csv")

# Keep only oracle rows (one per trial)
# If delta_star already represents oracle reduction per trial, no filtering needed

U = df["U_pred"].values
delta = df["delta_star"].values

# -----------------------------
# Bin calibration
# -----------------------------

n_bins = 8
bin_means, bin_edges, _ = binned_statistic(U, delta, statistic='mean', bins=n_bins)

# Compute bin centers
bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])

# Compute bootstrap CI per bin
ci_lower = []
ci_upper = []

for i in range(n_bins):
    mask = (U >= bin_edges[i]) & (U < bin_edges[i+1])
    values = delta[mask]
    
    if len(values) > 5:
        boot_means = []
        for _ in range(200):
            sample = np.random.choice(values, size=len(values), replace=True)
            boot_means.append(np.mean(sample))
        ci_lower.append(np.percentile(boot_means, 2.5))
        ci_upper.append(np.percentile(boot_means, 97.5))
    else:
        ci_lower.append(np.nan)
        ci_upper.append(np.nan)

# -----------------------------
# Plot
# -----------------------------

# Spearman correlation (robust, preferred by reviewers)
rho, pval = spearmanr(U, delta)

plt.figure(figsize=(6,5))

# Scatter (transparent)
plt.scatter(U, delta, alpha=0.25, s=12)

# Keep only bins with finite means (avoid broken CI bands)
finite = np.isfinite(bin_means) & np.isfinite(ci_lower) & np.isfinite(ci_upper)

# Binned mean line + CI band
plt.plot(bin_centers[finite], bin_means[finite], marker='o')
plt.fill_between(bin_centers[finite], np.array(ci_lower)[finite], np.array(ci_upper)[finite], alpha=0.2)

plt.xlabel(r"Predicted underdetermination $U_k$")
plt.ylabel(r"Observed ambiguity reduction $\Delta U_k(m^*)$")  # raw string fixes warning
plt.title("Calibration of underdetermination functional")

# Annotate correlation
plt.text(
    0.02, 0.98,
    rf"Spearman $\rho$ = {rho:.2f} (p={pval:.1e})",
    transform=plt.gca().transAxes,
    va="top"
)

plt.tight_layout()

out_path = Path("results/figures/Fig_multi_calibration.pdf")
out_path.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(out_path)
print("Saved", out_path)
