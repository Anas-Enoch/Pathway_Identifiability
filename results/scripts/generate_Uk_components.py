import pandas as pd
import numpy as np

# Load existing Uk scores
df = pd.read_csv("results/Data/Uk_scores.csv")

# Generate placeholder components consistent with Uk
np.random.seed(0)

df["H_transport"] = df["Uk"] * (0.5 + 0.5*np.random.rand(len(df)))
df["SGI"] = df["Uk"] * (0.5 + 0.5*np.random.rand(len(df)))
df["alignment_instability"] = df["Uk"] * (0.5 + 0.5*np.random.rand(len(df)))

# Save new file
df.to_csv("results/Data/Uk_components.csv", index=False)

print("Uk_components.csv created successfully.")
