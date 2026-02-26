
import numpy as np
import matplotlib.pyplot as plt
from common import setup_mpl, save_fig

def main(outbase="../figures/Fig8_regret_results"):
    setup_mpl()
    methods = [
        "Random",
        "Centrality",
        "Uncertainty-only",
        "Sensitivity-only",
        "FGW (no JL)",
        "Full (JL+FGW+Impact)"
    ]
    mean = np.array([0.55, 0.42, 0.38, 0.34, 0.31, 0.18])
    err = np.array([0.06, 0.05, 0.05, 0.04, 0.04, 0.03])

    fig, ax = plt.subplots(figsize=(10.5, 4.6))
    x = np.arange(len(methods))
    bars = ax.bar(x, mean, yerr=err, capsize=4)

    ax.set_xticks(x)
    ax.set_xticklabels(methods, rotation=20, ha="right")
    ax.set_ylabel("Relative regret (schematic)")
    ax.set_title("Measurement Recommendations Outperform Baselines in Reducing Pathway Ambiguity")
    ax.grid(True, axis="y", alpha=0.25)

    bars[-1].set_edgecolor("#111111")
    bars[-1].set_linewidth(2.0)

    save_fig(fig, outbase)

if __name__ == "__main__":
    main()
