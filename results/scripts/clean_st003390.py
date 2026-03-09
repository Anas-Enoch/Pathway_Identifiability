from pathlib import Path
import re
import pandas as pd

RAW_FILE = Path("results/data/raw/ST003390_AN005559.txt")
OUT_FILE = Path("results/data/ST003390_processed.csv")


def extract_ms_block(path: Path):
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()

    start = None
    end = None
    for i, line in enumerate(lines):
        s = line.strip()
        if s == "MS_METABOLITE_DATA_START":
            start = i + 1
        elif s == "MS_METABOLITE_DATA_END":
            end = i
            break

    if start is None or end is None:
        raise ValueError("Could not find MS_METABOLITE_DATA_START / END block")

    return [ln for ln in lines[start:end] if ln.strip()]


def main():
    block = extract_ms_block(RAW_FILE)

    # --- First two rows are tab-delimited and clean ---
    sample_row = block[0].split("\t")
    factor_row = block[1].split("\t")

    if sample_row[0].strip() != "Samples":
        raise ValueError(f"Expected first row to start with 'Samples', got: {sample_row[0]!r}")
    if factor_row[0].strip() != "Factors":
        raise ValueError(f"Expected second row to start with 'Factors', got: {factor_row[0]!r}")

    sample_ids = [x.strip() for x in sample_row[1:]]
    factor_vals = [x.strip() for x in factor_row[1:]]

    n_samples = len(sample_ids)

    if len(factor_vals) != n_samples:
        factor_vals = factor_vals[:n_samples]

    # phenotype extraction
    phenotypes = []
    for x in factor_vals:
        m = re.search(r"Phenotype:(.*)$", x)
        phenotypes.append(m.group(1).strip() if m else x.strip())

    # --- Metabolite rows are malformed: extract numeric values by regex ---
    metabolite_names = []
    metabolite_matrix = []

    num_pattern = re.compile(r"-?\d+\.\d{4}")

    for line in block[2:]:
        if "\t" not in line:
            continue

        met_name, numeric_part = line.split("\t", 1)
        met_name = met_name.strip()

        values = num_pattern.findall(numeric_part)

        if len(values) != n_samples:
            print(f"WARNING: {met_name} -> expected {n_samples} values, found {len(values)}")
            continue

        metabolite_names.append(met_name)
        metabolite_matrix.append([float(v) for v in values])

    if not metabolite_matrix:
        raise ValueError("No metabolite rows were parsed successfully")

    # metabolites x samples -> samples x metabolites
    met_df = pd.DataFrame(metabolite_matrix, index=metabolite_names, columns=sample_ids)
    sample_df = met_df.T

    sample_df.insert(0, "Phenotype", phenotypes)
    sample_df.index.name = "SampleID"
    sample_df.reset_index(inplace=True)

    OUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    sample_df.to_csv(OUT_FILE, index=False)

    print(f"Saved cleaned matrix to: {OUT_FILE}")
    print(sample_df.head().to_string())
    print(sample_df.shape)
    print(sample_df["Phenotype"].value_counts(dropna=False))


if __name__ == "__main__":
    main()
