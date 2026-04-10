import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parents[2]

input_csv = BASE_DIR / "results" / "results" / "retained_pathways_by_dataset.csv"
output_png = BASE_DIR / "results" / "results" / "fig_retained_pathways_heatmap.png"

df = pd.read_csv(input_csv)

# Pivot: rows = pathways, columns = datasets, values = n_detected
heat = df.pivot_table(
    index="pathway_name",
    columns="dataset",
    values="n_detected",
    aggfunc="max",
    fill_value=0
)

# Better pathway order: strongest recurring pathways first
row_order = [
    "Arginine and Proline Metabolism",
    "Glycine Serine Threonine Metabolism",
    "Alanine Aspartate Glutamate Metabolism",
    "Valine Leucine Isoleucine Degradation",
    "Pyrimidine Metabolism",
    "Purine Metabolism",
    "Tryptophan Metabolism",
    "Glutathione Metabolism",
    "Glycolysis / Gluconeogenesis",
    "Glycerophospholipid Metabolism",
    "Histidine Metabolism",
    "Taurine and Hypotaurine Metabolism",
    "TCA Cycle",
]
row_order = [r for r in row_order if r in heat.index]
heat = heat.reindex(row_order)

# Better dataset order
col_order = [c for c in ["ST000356", "ST003390", "ST003506"] if c in heat.columns]
heat = heat[col_order]

fig, ax = plt.subplots(figsize=(8, 6))
im = ax.imshow(heat.values, aspect="auto")

# ticks
ax.set_xticks(range(len(heat.columns)))
ax.set_xticklabels(heat.columns)
ax.set_yticks(range(len(heat.index)))
ax.set_yticklabels(heat.index)

ax.set_title("Cross-dataset pathway coverage")
ax.set_xlabel("Dataset")
ax.set_ylabel("Pathway")

# annotate values
for i in range(heat.shape[0]):
    for j in range(heat.shape[1]):
        ax.text(j, i, str(int(heat.iloc[i, j])), ha="center", va="center")

cbar = plt.colorbar(im, ax=ax)
cbar.set_label("Detected metabolites")

plt.tight_layout()
plt.savefig(output_png, dpi=300, bbox_inches="tight")
print(f"Saved: {output_png}")
print(heat)
