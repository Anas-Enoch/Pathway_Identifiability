import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parents[2]
input_csv = BASE_DIR / "results" / "results" / "pathway_coverage_summary.csv"
output_png = BASE_DIR / "results" / "results" / "fig_pathway_coverage_clean.png"

df = pd.read_csv(input_csv)

pivot = df.pivot(index="pathway_name", columns="dataset", values="n_metabolites").fillna(0)

# Keep a biologically cleaner order
preferred_order = [
    "Arginine and Proline Metabolism",
    "Glycine Serine Threonine Metabolism",
    "Valine Leucine Isoleucine Degradation",
    "Alanine Aspartate Glutamate Metabolism",
    "Glycolysis / Gluconeogenesis",
    "Glutathione Metabolism",
    "Purine Metabolism",
    "Pyrimidine Metabolism",
    "Tryptophan Metabolism",
    "Glycerophospholipid Metabolism",
    "Histidine Metabolism",
    "TCA Cycle",
    "Taurine and Hypotaurine Metabolism",
]

pivot = pivot.reindex([p for p in preferred_order if p in pivot.index])

ax = pivot.plot(kind="bar", figsize=(14, 7))

ax.set_ylabel("Detected metabolites")
ax.set_xlabel("")
ax.set_title("Pathway coverage across metabolomics cohorts")
plt.xticks(rotation=35, ha="right")
plt.tight_layout()

plt.savefig(output_png, dpi=300, bbox_inches="tight")
print(f"Figure saved to: {output_png}")
