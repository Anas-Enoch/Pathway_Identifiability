
import numpy as np
import matplotlib.pyplot as plt
from common import setup_mpl, save_fig, OBS_COLOR, MISS_EDGE, LOD_COLOR

def normal_pdf(x, mu, sigma):
    return (1.0/(sigma*np.sqrt(2*np.pi))) * np.exp(-0.5*((x-mu)/sigma)**2)

def truncated_normal_pdf(x, mu, sigma, a):
    pdf = normal_pdf(x, mu, sigma)
    mask = x >= a
    area = np.trapezoid(pdf[mask], x[mask]) + 1e-12
    out = np.where(mask, pdf/area, 0.0)
    return out

def main(outbase="../figures/Fig3_distributional_node_features"):
    setup_mpl()
    x = np.linspace(-3, 6, 600)

    fig, axes = plt.subplots(1, 3, figsize=(12, 3.6), sharey=True)

    mu1, s1 = 1.0, 0.35
    y1 = normal_pdf(x, mu1, s1)
    axes[0].plot(x, y1, lw=2.5, color=OBS_COLOR)
    axes[0].fill_between(x, 0, y1, alpha=0.25, color=OBS_COLOR)
    axes[0].set_title("Observed metabolite")
    axes[0].text(0.02, 0.95, "$\\mu$ measured\n$\\sigma^2\\approx\\sigma^2_{assay}$",
                 transform=axes[0].transAxes, va="top")

    mu2, s2 = 1.0, 1.3
    y2 = normal_pdf(x, mu2, s2)
    axes[1].plot(x, y2, lw=2.5, color=MISS_EDGE)
    axes[1].fill_between(x, 0, y2, alpha=0.25, color=MISS_EDGE)
    axes[1].set_title("Panel-missing metabolite")
    axes[1].text(0.02, 0.95, "$\\mu$ weak prior\n$\\sigma^2$ large",
                 transform=axes[1].transAxes, va="top")

    mu3, s3 = 0.2, 0.9
    lod = 0.8
    y3 = truncated_normal_pdf(x, mu3, s3, lod)
    axes[2].plot(x, y3, lw=2.5, color=LOD_COLOR)
    axes[2].fill_between(x, 0, y3, alpha=0.25, color=LOD_COLOR)
    axes[2].axvline(lod, lw=1.8, linestyle="--", color=LOD_COLOR)
    axes[2].set_title("LOD-censored metabolite")
    axes[2].text(0.02, 0.95, "left-censored at LOD\n(truncated posterior)",
                 transform=axes[2].transAxes, va="top")

    for ax in axes:
        ax.set_xlabel("Abundance (a.u.)")
        ax.grid(True, alpha=0.25)
    axes[0].set_ylabel("Density")

    fig.suptitle("Distributional Node Features Enable Uncertainty-Aware Comparison", y=1.05)
    save_fig(fig, outbase)

if __name__ == "__main__":
    main()
