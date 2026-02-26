import numpy as np
import matplotlib.pyplot as plt

# -------------------------
# Dummy numbers (replace later with real outputs)
# -------------------------

# Panel B: underdetermination score
U_before = 1.32
U_after = 0.78

# Panel C: measurement impact ranking
metabolites = [
    "Ribose-5P",
    "6-Phosphogluconate",
    "Sedoheptulose-7P",
    "Glyceraldehyde-3P"
]

impact_scores = [0.42, 0.31, 0.18, 0.12]

# -------------------------
# Create figure
# -------------------------

fig = plt.figure(figsize=(10, 8))

# Panel A — schematic (simple text diagram)
ax1 = plt.subplot(2, 2, 1)
ax1.set_title("A) Observed coverage (G6P branch)")
ax1.text(0.5, 0.7, "G6P", ha='center')
ax1.text(0.2, 0.4, "F6P\n(glycolysis)", ha='center')
ax1.text(0.8, 0.4, "R5P\n(PPP)", ha='center')

ax1.plot([0.5, 0.2], [0.65, 0.45])
ax1.plot([0.5, 0.8], [0.65, 0.45])

ax1.text(0.8, 0.3, "(latent)", ha='center')
ax1.set_xlim(0, 1)
ax1.set_ylim(0, 1)
ax1.axis('off')

# Panel B — U score
ax2 = plt.subplot(2, 2, 2)
ax2.set_title("B) Underdetermination score")
ax2.bar(["Observed panel"], [U_before])
ax2.set_ylabel("U_k")

# Panel C — Impact ranking
ax3 = plt.subplot(2, 2, 3)
ax3.set_title("C) Measurement impact ranking")
ax3.bar(metabolites, impact_scores)
ax3.set_ylabel("Estimated ΔU")
ax3.tick_params(axis='x', rotation=45)

# Panel D — reduction after reveal
ax4 = plt.subplot(2, 2, 4)
ax4.set_title("D) Reduction after reveal")
ax4.bar(["Before", "After"], [U_before, U_after])
ax4.set_ylabel("U_k")

plt.tight_layout()
plt.savefig("Fig_real_case_PPP.pdf")
plt.close()

print("Figure saved as Fig_real_case_PPP.pdf")
