import pandas as pd
from pathlib import Path

results_dir = Path("results/results")

coverage_files = [
    results_dir / "pathway_coverage_summary.csv",            # existing combined file for ST000356 + ST003506
    results_dir / "pathway_coverage_summary_ST003390.csv",   # new ST003390 file
]

frames = []

for f in coverage_files:
    if not f.exists():
        print(f"WARNING: missing file -> {f}")
        continue

    df = pd.read_csv(f)

    # standardize possible column names
    rename_map = {}
    if "n_metabolites" in df.columns and "n_detected" not in df.columns:
        rename_map["n_metabolites"] = "n_detected"
    df = df.rename(columns=rename_map)

    needed = {"dataset", "pathway_id", "pathway_name", "n_detected"}
    missing = needed - set(df.columns)
    if missing:
        raise ValueError(f"{f} is missing columns: {missing}")

    out = df[["dataset", "pathway_id", "pathway_name", "n_detected"]].copy()
    out["is_retained_main_benchmark"] = out["n_detected"] >= 3
    frames.append(out)

if not frames:
    raise ValueError("No coverage files found")

merged = pd.concat(frames, ignore_index=True)

# remove exact duplicates if any
merged = merged.drop_duplicates(
    subset=["dataset", "pathway_id", "pathway_name"],
    keep="last"
)

merged = merged.sort_values(
    ["dataset", "is_retained_main_benchmark", "n_detected", "pathway_name"],
    ascending=[True, False, False, True]
)

output_file = results_dir / "retained_pathways_by_dataset.csv"
merged.to_csv(output_file, index=False)

print(f"Saved: {output_file}")
print()
print(merged.to_string(index=False))
