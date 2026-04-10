import pandas as pd
import numpy as np
from pathlib import Path

DATA_DIR = Path("results/data")
OUT_DIR = Path("results/results")
OUT_DIR.mkdir(parents=True, exist_ok=True)

DATASETS = {
    "ST003506": DATA_DIR / "ST003506_serum_processed.csv",
    "ST000356": DATA_DIR / "processed_metabolite_matrix_ST000356.csv",
}

MAPPING_FILE = DATA_DIR / "core_pathway_mapping.csv"


def normalize_group(x: str) -> str:
    x = str(x).strip().lower()
    if "control" in x:
        return "control"
    if "lymphedema" in x:
        return "disease"
    if "cancer" in x or "breast cancer" in x:
        return "disease"
    return x


def load_dataset(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df.columns = [str(c).strip() for c in df.columns]
    if "sample_id" not in df.columns or "group" not in df.columns:
        raise ValueError(f"{path.name} must contain sample_id and group columns")
    df["group"] = df["group"].map(normalize_group)
    return df


def load_mapping(path: Path) -> pd.DataFrame:
    mp = pd.read_csv(path)
    mp.columns = [str(c).strip() for c in mp.columns]
    required = {"metabolite_name", "pathway_id", "pathway_name"}
    missing = required - set(mp.columns)
    if missing:
        raise ValueError(f"Mapping file missing columns: {missing}")
    return mp


def pathway_coverage(df: pd.DataFrame, mp: pd.DataFrame, dataset_name: str) -> pd.DataFrame:
    metabolite_cols = set(df.columns) - {"sample_id", "group"}
    sub = mp[mp["metabolite_name"].isin(metabolite_cols)].copy()
    cov = (
        sub.groupby(["pathway_id", "pathway_name"], as_index=False)
        .agg(n_metabolites=("metabolite_name", "nunique"))
    )
    cov["dataset"] = dataset_name
    return cov[["dataset", "pathway_id", "pathway_name", "n_metabolites"]]


def compute_pathway_scores(df: pd.DataFrame, mp: pd.DataFrame, dataset_name: str) -> pd.DataFrame:
    metabolite_cols = set(df.columns) - {"sample_id", "group"}
    sub = mp[mp["metabolite_name"].isin(metabolite_cols)].copy()

    out_rows = []
    for (pid, pname), g in sub.groupby(["pathway_id", "pathway_name"]):
        mets = [m for m in g["metabolite_name"].unique() if m in df.columns]
        if len(mets) < 2:
            continue

        X = df[mets].apply(pd.to_numeric, errors="coerce")
        row_mean = X.mean(axis=1, skipna=True)

        ctrl = row_mean[df["group"] == "control"]
        dis = row_mean[df["group"] == "disease"]

        if len(ctrl) < 2 or len(dis) < 2:
            continue

        score = dis.mean() - ctrl.mean()

        out_rows.append(
            {
                "dataset": dataset_name,
                "pathway_id": pid,
                "pathway_name": pname,
                "n_metabolites": len(mets),
                "mean_control": ctrl.mean(),
                "mean_disease": dis.mean(),
                "pathway_score": score,
                "abs_pathway_score": abs(score),
            }
        )

    return pd.DataFrame(out_rows)


def masking_benchmark(df: pd.DataFrame, mp: pd.DataFrame, dataset_name: str,
                      rho=0.3, n_repeats=50, min_mets=3, seed=42) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    metabolite_cols = set(df.columns) - {"sample_id", "group"}
    sub = mp[mp["metabolite_name"].isin(metabolite_cols)].copy()

    rows = []
    for (pid, pname), g in sub.groupby(["pathway_id", "pathway_name"]):
        mets = [m for m in g["metabolite_name"].unique() if m in df.columns]
        if len(mets) < min_mets:
            continue

        Xfull = df[mets].apply(pd.to_numeric, errors="coerce")

        for r in range(n_repeats):
            n_mask = max(1, int(np.floor(rho * len(mets))))
            if len(mets) - n_mask < min_mets - 1:
                n_mask = max(1, len(mets) - (min_mets - 1))
            masked = list(rng.choice(mets, size=n_mask, replace=False))
            observed = [m for m in mets if m not in masked]

            base = Xfull[observed].mean(axis=1, skipna=True)
            base_ctrl = base[df["group"] == "control"].mean()
            base_dis = base[df["group"] == "disease"].mean()
            base_score = abs(base_dis - base_ctrl)

            oracle_delta = -np.inf
            best_met = None

            for m in masked:
                score_m = Xfull[observed + [m]].mean(axis=1, skipna=True)
                ctrl_m = score_m[df["group"] == "control"].mean()
                dis_m = score_m[df["group"] == "disease"].mean()
                delta = abs(abs(dis_m - ctrl_m) - base_score)

                if delta > oracle_delta:
                    oracle_delta = delta
                    best_met = m

            var_obs = Xfull[observed].var(axis=0, skipna=True)
            pred_met = var_obs.idxmax() if len(var_obs) else observed[0]

            if pred_met in masked:
                pred_delta = oracle_delta if pred_met == best_met else 0.0
            else:
                score_p = Xfull[observed + [pred_met]].mean(axis=1, skipna=True)
                ctrl_p = score_p[df["group"] == "control"].mean()
                dis_p = score_p[df["group"] == "disease"].mean()
                pred_delta = abs(abs(dis_p - ctrl_p) - base_score)

            regret = oracle_delta - pred_delta
            nregret = regret / (oracle_delta + 1e-8) if oracle_delta > 0 else 0.0

            rows.append(
                {
                    "dataset": dataset_name,
                    "pathway_id": pid,
                    "pathway_name": pname,
                    "rho": rho,
                    "repeat": r,
                    "n_metabolites": len(mets),
                    "masked_metabolites": ";".join(masked),
                    "oracle_metabolite": best_met,
                    "predicted_metabolite": pred_met,
                    "oracle_delta": oracle_delta,
                    "pred_delta": pred_delta,
                    "regret": regret,
                    "nregret": nregret,
                }
            )

    return pd.DataFrame(rows)


def main():
    mp = load_mapping(MAPPING_FILE)

    cov_all = []
    scores_all = []
    bench_all = []

    for name, path in DATASETS.items():
        df = load_dataset(path)

        cov = pathway_coverage(df, mp, name)
        cov_all.append(cov)

        scores = compute_pathway_scores(df, mp, name)
        scores_all.append(scores)

        bench = masking_benchmark(df, mp, name, rho=0.3, n_repeats=50, min_mets=3, seed=42)
        bench_all.append(bench)

    cov_all = pd.concat(cov_all, ignore_index=True)
    scores_all = pd.concat(scores_all, ignore_index=True)
    bench_all = pd.concat(bench_all, ignore_index=True)

    cov_all.to_csv(OUT_DIR / "pathway_coverage_summary.csv", index=False)
    scores_all.to_csv(OUT_DIR / "pathway_scores_summary.csv", index=False)
    bench_all.to_csv(OUT_DIR / "core_benchmark_results.csv", index=False)

    print("Saved:")
    print(OUT_DIR / "pathway_coverage_summary.csv")
    print(OUT_DIR / "pathway_scores_summary.csv")
    print(OUT_DIR / "core_benchmark_results.csv")

    print("\nCoverage per dataset:")
    print(cov_all.sort_values(["dataset", "n_metabolites"], ascending=[True, False]).to_string(index=False))

    print("\nBenchmark summary:")
    print(
        bench_all.groupby("dataset", as_index=False)
        .agg(
            mean_regret=("regret", "mean"),
            mean_nregret=("nregret", "mean"),
            n_trials=("repeat", "count"),
        )
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
