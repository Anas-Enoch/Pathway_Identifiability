#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
parse_ST002829.py

Converts the 4 ST002829 mwTab txt files into two CSVs:

  ST002829_processed_raw.csv
      rows = all samples, cols = sample_id + all metabolites
      includes every sample (pre / during / post disease)

  ST002829_benchmark_ready.csv
      rows = "during" samples only (collected at active disease)
      col  = condition (severe | mild) + all metabolites
      ready to drop into run_real_multi_pathway_benchmark.py as-is

Study: Nucleotide, phospholipid, and kynurenine metabolites robustly
       associated with COVID-19 severity (LEOCC cohort, n > 600)
Source: Metabolomics Workbench ST002829

Run from your repo root:
    python3 results/scripts/parse_ST002829.py

Place the 4 txt files in a folder called ST002829/ at repo root:
    ST002829/ST002829_AN004619.txt
    ST002829/ST002829_AN004620.txt
    ST002829/ST002829_AN004621.txt
    ST002829/ST002829_AN004622.txt
"""

from __future__ import annotations
from pathlib import Path
import re
import sys
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# ── paths (auto-resolved from script location) ────────────────────────
_HERE = Path(__file__).resolve().parent
_ROOT = _HERE
for _ in range(5):
    if (_ROOT / "results" / "data").exists():
        break
    _ROOT = _ROOT.parent

TXT_DIR  = _ROOT / "ST002829"
RAW_CSV  = _ROOT / "results" / "data" / "ST002829_processed_raw.csv"
BENCH_CSV = _ROOT / "results" / "data" / "ST002829_benchmark_ready.csv"
# ─────────────────────────────────────────────────────────────────────

BAD_ROWS = {
    "", "Factors", "Factor", "factor", "factors",
    "Samples", "Sample", "sample", "samples",
    "Metabolite", "metabolite", "Metabolite_name",
}


# ── mwTab section finder ──────────────────────────────────────────────

def find_section(lines: list[str], keys: list[str]):
    for key in keys:
        start = end = None
        for i, line in enumerate(lines):
            s = line.strip()
            if s == f"{key}_START":
                start = i + 1
            elif s == f"{key}_END":
                end = i
                break
        if start is not None and end is not None:
            return start, end, key
    return None, None, None


def normalize_block(block: list[list[str]]):
    header = block[0]
    ncols  = len(header)
    rows = []
    for r in block[1:]:
        if len(r) > ncols:
            r = r[:ncols]
        elif len(r) < ncols:
            r = r + [""] * (ncols - len(r))
        rows.append(r)
    return header, rows


# ── metabolite data parser ────────────────────────────────────────────

def parse_metabolite_data(path: Path) -> pd.DataFrame:
    """Return (n_samples, n_metabolites+1) dataframe with sample_id column."""
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        lines = fh.readlines()

    start, end, key = find_section(
        lines, ["MS_METABOLITE_DATA", "METABOLITE_DATA", "NMR_METABOLITE_DATA"]
    )
    if start is None:
        raise ValueError(f"No metabolite data section found in {path}")

    raw = []
    for line in lines[start:end]:
        line = line.rstrip("\n")
        if not line.strip():
            continue
        raw.append(line.split("\t"))

    if len(raw) < 2:
        raise ValueError(f"Metabolite block too short in {path}")

    header, rows = normalize_block(raw)
    df = pd.DataFrame(rows, columns=header)

    # standardise first column name
    df = df.rename(columns={df.columns[0]: "metabolite"})
    df["metabolite"] = df["metabolite"].astype(str).str.strip()

    # drop junk rows
    df = df[~df["metabolite"].isin(BAD_ROWS)].copy()
    df = df[~df["metabolite"].str.startswith("#", na=False)].copy()
    df = df.drop_duplicates(subset=["metabolite"], keep="first")

    # transpose → samples as rows
    df_t = df.set_index("metabolite").T.reset_index()
    df_t = df_t.rename(columns={"index": "sample_id"})

    for c in df_t.columns:
        if c != "sample_id":
            df_t[c] = pd.to_numeric(df_t[c], errors="coerce")

    return df_t


# ── SUBJECT_SAMPLE_FACTORS parser ────────────────────────────────────

def parse_factors(path: Path) -> pd.DataFrame:
    """
    Extract {sample_id, severity, time_point} from SUBJECT_SAMPLE_FACTORS.

    ST002829 factor string format (typical):
        "COVID_severity:Severe | Time_point:During"
    or
        "Severity:Mild | Collection_time:Before"

    Returns DataFrame with cols: sample_id, severity, time_point
    """
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        lines = fh.readlines()

    in_block = False
    records  = []

    for line in lines:
        s = line.strip()
        if s == "SUBJECT_SAMPLE_FACTORS_START":
            in_block = True
            continue
        if s == "SUBJECT_SAMPLE_FACTORS_END":
            in_block = False
            continue
        if not in_block or not s or s.startswith("#"):
            continue

        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue

        sample_id = parts[2].strip()
        factor_str = "\t".join(parts[3:])   # everything after sample_id

        severity   = _parse_severity(factor_str)
        time_point = _parse_timepoint(factor_str)

        if sample_id:
            records.append({
                "sample_id":  sample_id,
                "severity":   severity,
                "time_point": time_point,
            })

    return pd.DataFrame(records) if records else pd.DataFrame(
        columns=["sample_id", "severity", "time_point"]
    )


def _parse_severity(s: str) -> str | None:
    sl = s.lower()
    # explicit labels
    if re.search(r'\bsevere\b', sl):   return "severe"
    if re.search(r'\bmild\b', sl):     return "mild"
    if re.search(r'\bmoderate\b', sl): return "moderate"
    if re.search(r'\bnon-severe\b|non_severe', sl): return "mild"
    if re.search(r'\bcritical\b|\bicu\b|\bintubat', sl): return "severe"
    if re.search(r'\bhealthy\b|\bcontrol\b|\bnormal\b', sl): return "control"
    # WHO numeric
    m = re.search(r'who[:\s_-]*(\d)', sl)
    if m:
        score = int(m.group(1))
        if score <= 3: return "mild"
        if score <= 5: return "moderate"
        return "severe"
    return None


def _parse_timepoint(s: str) -> str | None:
    sl = s.lower()
    if re.search(r'\bbefore\b|\bpre[_\s-]?(disease|covid|diagnosis)\b|\bpre\b', sl):
        return "before"
    if re.search(r'\bduring\b|\bactive\b|\bacute\b|\bat[_\s]diagnosis\b', sl):
        return "during"
    if re.search(r'\bafter\b|\bpost[_\s-]?(disease|covid)\b|\brecov', sl):
        return "after"
    return None


# ── combine all txt files ─────────────────────────────────────────────

def load_all(txt_dir: Path):
    txt_files = sorted(txt_dir.glob("ST002829_AN*.txt"))
    if not txt_files:
        raise FileNotFoundError(
            f"No ST002829_AN*.txt files found in {txt_dir}\n"
            "Expected: ST002829_AN004619.txt  AN004620  AN004621  AN004622"
        )

    print(f"Found {len(txt_files)} analysis file(s):")

    met_dfs  = []
    fact_dfs = []

    for f in txt_files:
        print(f"  Parsing {f.name} ...")
        try:
            mdf = parse_metabolite_data(f)
            print(f"    metabolite block → {mdf.shape[0]} samples × {mdf.shape[1]-1} metabolites")
            met_dfs.append(mdf.set_index("sample_id"))
        except Exception as e:
            print(f"    [WARN] metabolite parse failed: {e}")

        fdf = parse_factors(f)
        if len(fdf):
            print(f"    factors block     → {len(fdf)} rows")
        fact_dfs.append(fdf)

    # combine metabolite matrices
    combined = pd.concat(met_dfs, axis=1)
    # deduplicate metabolite columns — keep first non-null value
    combined = combined.T.groupby(level=0).first().T
    combined = combined.dropna(axis=1, how="all")
    combined = combined.reset_index()   # sample_id back as column

    # combine factor metadata
    factors = pd.concat(fact_dfs, ignore_index=True).drop_duplicates("sample_id")

    return combined, factors


# ── severity normalisation for benchmark ─────────────────────────────

def normalise_severity(sev: str | None) -> str | None:
    """Map severity variants to 'severe' | 'mild' | None (drop)."""
    if sev is None:
        return None
    s = str(sev).lower().strip()
    if s in ("severe", "critical"):
        return "severe"
    if s in ("mild", "non-severe", "moderate"):
        # treat moderate as mild for a 2-class benchmark
        # change this if you want a 3-class split
        return "mild"
    if s == "control":
        return None   # exclude healthy controls from severity benchmark
    return None


# ── PCA fallback when metadata is missing ────────────────────────────

def pca_severity_fallback(bench: pd.DataFrame, met_cols: list[str]) -> pd.Series:
    print("[WARN] No severity metadata matched. Using PCA severity fallback.")
    X = bench[met_cols].values.astype(float)
    col_med = np.nanmedian(X, axis=0)
    for j in range(X.shape[1]):
        X[np.isnan(X[:, j]), j] = col_med[j]
    X_sc = StandardScaler().fit_transform(np.log1p(X))
    pc1  = PCA(n_components=1).fit_transform(X_sc)[:, 0]
    labels = pd.Series(
        np.where(pc1 > np.median(pc1), "severe", "mild"),
        index=bench.index,
    )
    print(f"  PCA split — mild: {(labels=='mild').sum()}  severe: {(labels=='severe').sum()}")
    return labels


# ── main ──────────────────────────────────────────────────────────────

def main() -> None:

    if not TXT_DIR.exists():
        sys.exit(
            f"[ERROR] TXT_DIR not found: {TXT_DIR}\n"
            "        Create a folder called ST002829/ at repo root and place\n"
            "        the 4 AN*.txt files inside it."
        )

    # 1. parse all files
    combined, factors = load_all(TXT_DIR)
    print(f"\nCombined metabolite matrix: {combined.shape}")
    print(f"Factor metadata rows:       {len(factors)}")

    # 2. save raw CSV (all samples, all timepoints)
    RAW_CSV.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(RAW_CSV, index=False)
    print(f"\nRaw CSV saved → {RAW_CSV}")

    # 3. show timepoint and severity breakdown
    merged = combined[["sample_id"]].merge(factors, on="sample_id", how="left")
    print("\nTimepoint counts:")
    print(merged["time_point"].value_counts(dropna=False).to_string())
    print("\nSeverity counts:")
    print(merged["severity"].value_counts(dropna=False).to_string())

    # 4. select "during" samples for severity benchmark
    #    (closest biological match to ST001849 d0 baseline)
    during_ids = merged[merged["time_point"] == "during"]["sample_id"]
    if len(during_ids) == 0:
        print("\n[WARN] No 'during' samples found. Using all samples with severity label.")
        bench_raw = combined.merge(
            factors[["sample_id","severity"]], on="sample_id", how="inner"
        )
    else:
        print(f"\n'During' samples: {len(during_ids)}")
        bench_raw = combined[combined["sample_id"].isin(during_ids)].merge(
            factors[["sample_id","severity"]], on="sample_id", how="left"
        )

    # 5. normalise severity labels → severe | mild
    bench_raw["condition"] = bench_raw["severity"].apply(normalise_severity)
    n_labelled = bench_raw["condition"].notna().sum()
    print(f"Labelled samples: {n_labelled}/{len(bench_raw)}")

    met_cols = [c for c in combined.columns if c != "sample_id"]

    if n_labelled < 0.6 * len(bench_raw):
        # fallback: use all "during" samples and derive severity from PCA
        bench_pca = bench_raw.dropna(subset=met_cols, how="all").reset_index(drop=True)
        bench_pca["condition"] = pca_severity_fallback(bench_pca, met_cols)
        bench_out = bench_pca
    else:
        bench_out = bench_raw.dropna(subset=["condition"]).reset_index(drop=True)

    # 6. validate two-class split
    counts = bench_out["condition"].value_counts()
    print("\nBenchmark condition split:")
    print(counts.to_string())
    if len(counts) != 2:
        sys.exit(f"[ERROR] Expected 2 conditions, got: {counts.index.tolist()}")
    if counts.min() < 10:
        print(f"[WARN] Smallest group n={counts.min()} — check metadata.")

    # 7. build output (condition first, no sample_id)
    met_cols_avail = [c for c in met_cols if c in bench_out.columns]
    out = bench_out[["condition"] + met_cols_avail].copy()

    # verify benchmark parsing compatibility
    assert "condition" in out.columns
    assert len(out["condition"].unique()) == 2
    numeric_ok = all(pd.api.types.is_numeric_dtype(out[c]) for c in met_cols_avail)
    assert numeric_ok, "Non-numeric metabolite columns found"

    # 8. save benchmark CSV
    out.to_csv(BENCH_CSV, index=False)
    print(f"\nBenchmark CSV saved → {BENCH_CSV}")
    print(f"  rows       : {len(out)}")
    print(f"  metabolites: {len(met_cols_avail)}")
    print(f"  condition  : {sorted(out['condition'].unique())}")

    print("\n── Next steps ──────────────────────────────────────────────────")
    print("1. Copy ST002829_pathway_mapping.csv to results/data/")
    print("   (or reuse ST001849_pathway_mapping.csv if metabolite overlap is sufficient)")
    print("2. Add to DATASETS in your benchmark script:")
    print(f'   "ST002829": DATA_DIR / "{BENCH_CSV.name}",')
    print("3. Add to CASE_CTRL_OVERRIDE:")
    print('   "ST002829": ("severe", "mild"),')


if __name__ == "__main__":
    main()
