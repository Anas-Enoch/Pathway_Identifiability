import pandas as pd
from pathlib import Path
import re

mapping_file = Path("results/data/core_pathway_mapping.csv")
cohort_file = Path("results/data/ST003390_processed.csv")
output_hits = Path("results/results/st003390_mapping_hits.csv")
output_cov = Path("results/results/pathway_coverage_summary_ST003390.csv")


def normalize_name(x: str) -> str:
    x = str(x).strip().lower()
    x = x.replace("α", "alpha")
    x = x.replace("β", "beta")
    x = x.replace("γ", "gamma")
    x = x.replace("^", "")
    x = x.replace("_", " ")
    x = x.replace("-", " ")
    x = re.sub(r"[(),/]", " ", x)
    x = re.sub(r"\s+", " ", x).strip()

    synonyms = {
        "l proline": "proline",
        "glutamic acid": "glutamate",
        "pyroglutamic acid": "pyroglutamate",
        "indoleacetic acid": "indoleacetic acid",
        "uric acid": "uric acid",
        "aspartic acid": "aspartate",
        "alpha ketoglutaric acid": "alpha ketoglutarate",
        "2 o methyluridine": "2 o methyluridine",
        "n n dimethylglycine": "dimethylglycine",
    }
    return synonyms.get(x, x)


# load files
mapping = pd.read_csv(mapping_file)
cohort = pd.read_csv(cohort_file)

# cohort metabolite columns
cohort_metabolites = [c for c in cohort.columns if c not in ["SampleID", "Phenotype"]]

mapping["norm_name"] = mapping["metabolite_name"].map(normalize_name)
cohort_norm = pd.DataFrame({
    "cohort_metabolite": cohort_metabolites,
    "norm_name": [normalize_name(x) for x in cohort_metabolites]
})

# exact normalized overlap
hits = mapping.merge(cohort_norm, on="norm_name", how="left")
hits["is_detected"] = hits["cohort_metabolite"].notna()

# save detailed hit table
output_hits.parent.mkdir(parents=True, exist_ok=True)
hits.to_csv(output_hits, index=False)

# pathway coverage
coverage = (
    hits.groupby(["pathway_id", "pathway_name"], as_index=False)
        .agg(
            n_mapped=("metabolite_name", "count"),
            n_detected=("is_detected", "sum")
        )
        .sort_values(["n_detected", "n_mapped"], ascending=[False, False])
)

coverage["dataset"] = "ST003390"
coverage["is_usable_ge_3"] = coverage["n_detected"] >= 3
coverage["is_strong_ge_4"] = coverage["n_detected"] >= 4

coverage.to_csv(output_cov, index=False)

print("\n=== OVERLAP SUMMARY ===")
print(f"Dataset metabolites: {len(cohort_metabolites)}")
print(f"Mapping rows: {len(mapping)}")
print(f"Matched mapping rows: {hits['is_detected'].sum()}")

print("\n=== PATHWAY COVERAGE ===")
print(coverage.to_string(index=False))

print("\n=== USABLE PATHWAYS (>=3 detected metabolites) ===")
usable = coverage[coverage["n_detected"] >= 3]
print(usable[["pathway_id", "pathway_name", "n_detected"]].to_string(index=False))

print(f"\nSaved detailed hits to: {output_hits}")
print(f"Saved pathway coverage to: {output_cov}")
