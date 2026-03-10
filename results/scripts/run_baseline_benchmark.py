from pathlib import Path
import math
import random
import re

import numpy as np
import pandas as pd


# =========================
# Paths
# =========================
BASE_DIR = Path(__file__).resolve().parents[2]

DATA_DIR = BASE_DIR / "results" / "data"
RESULTS_DIR = BASE_DIR / "results" / "results"

MAPPING_FILE = DATA_DIR / "core_pathway_mapping.csv"
RETAINED_FILE = RESULTS_DIR / "retained_pathways_by_dataset.csv"

DATASET_FILES = {
    "ST000356": DATA_DIR / "processed_metabolite_matrix_ST000356.csv",
    "ST003390": DATA_DIR / "ST003390_processed.csv",
    "ST003506": DATA_DIR / "ST003506_serum_processed.csv",
}

TRIAL_OUT = RESULTS_DIR / "baseline_trial_results.csv"
REGRET_OUT = RESULTS_DIR / "baseline_regret_summary.csv"
TOPK_OUT = RESULTS_DIR / "baseline_topk_summary.csv"

TRIAL_OUT_NONTRIVIAL = RESULTS_DIR / "baseline_trial_results_nontrivial.csv"
REGRET_OUT_NONTRIVIAL = RESULTS_DIR / "baseline_regret_summary_nontrivial.csv"
TOPK_OUT_NONTRIVIAL = RESULTS_DIR / "baseline_topk_summary_nontrivial.csv"
GLOBAL_OUT_NONTRIVIAL = RESULTS_DIR / "baseline_global_summary_nontrivial.csv"


# =========================
# Config
# =========================
MASK_RATES = [0.4, 0.5, 0.6]
N_TRIALS = 100
METHODS = ["proposed", "random", "degree", "imputation"]
RANDOM_SEED = 42
MIN_NONTRIVIAL_MASKED = 3


# =========================
# Utilities
# =========================
def normalize_name(x: str) -> str:
    x = str(x).strip().lower()
    x = x.replace("α", "alpha").replace("β", "beta").replace("γ", "gamma")
    x = x.replace("_", " ").replace("-", " ").replace("/", " ")
    x = re.sub(r"[(),]", " ", x)
    x = re.sub(r"\s+", " ", x).strip()

    synonyms = {
        "glutamic acid": "glutamate",
        "aspartic acid": "aspartate",
        "alpha ketoglutaric acid": "alpha ketoglutarate",
        "α ketoglutaric acid": "alpha ketoglutarate",
        "pyroglutamic acid": "pyroglutamate",
        "l proline": "proline",
        "hydroxyproline": "4 hydroxy proline",
        "n n dimethylglycine": "dimethylglycine",
        "alpha aminobutyric acid": "2 aminobutyric acid",
        "uric acid": "uric acid",
        "indoleacetic acid": "indoleacetic acid",
    }
    return synonyms.get(x, x)


def load_mapping(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    needed = {"metabolite_name", "pathway_id", "pathway_name"}
    missing = needed - set(df.columns)
    if missing:
        raise ValueError(f"Mapping file missing columns: {missing}")
    df["norm_name"] = df["metabolite_name"].map(normalize_name)
    return df


def load_retained(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    needed = {"dataset", "pathway_id", "pathway_name", "is_retained_main_benchmark"}
    missing = needed - set(df.columns)
    if missing:
        raise ValueError(f"Retained file missing columns: {missing}")
    df = df[df["is_retained_main_benchmark"] == True].copy()
    return df


def load_dataset(dataset_name: str, file_path: Path) -> pd.DataFrame:
    df = pd.read_csv(file_path)
    df.columns = [str(c).strip() for c in df.columns]

    # SampleID
    if "SampleID" not in df.columns:
        sample_candidates = [
            "sampleid", "sample_id", "sample", "samples", "id", "subject", "subject_id"
        ]
        found = None
        for c in df.columns:
            if c.lower() in sample_candidates:
                found = c
                break
        if found is not None:
            df = df.rename(columns={found: "SampleID"})
        else:
            df.insert(0, "SampleID", [f"{dataset_name}_{i+1}" for i in range(len(df))])

    # Phenotype
    if "Phenotype" not in df.columns:
        pheno_candidates = ["phenotype", "group", "class", "label", "condition", "status", "diagnosis"]
        found = None
        for c in df.columns:
            if c.lower() in pheno_candidates:
                found = c
                break
        if found is not None:
            df = df.rename(columns={found: "Phenotype"})
        else:
            sid = df["SampleID"].astype(str).str.lower()
            if dataset_name == "ST000356":
                phenotype = np.where(
                    sid.str.contains("control|ctrl|healthy|normal|hc"),
                    "Control",
                    np.where(sid.str.contains("cancer|tumor|case|bc"), "Cancer", "Unknown")
                )
                df.insert(1, "Phenotype", phenotype)
            else:
                df.insert(1, "Phenotype", "Unknown")

    return df


def get_dataset_metabolite_columns(df: pd.DataFrame) -> list[str]:
    return [c for c in df.columns if c not in ["SampleID", "Phenotype"]]


def build_pathway_metabolite_map(mapping_df: pd.DataFrame, dataset_metabs: list[str]) -> dict:
    dataset_lookup = {normalize_name(m): m for m in dataset_metabs}
    out = {}
    for (pid, pname), sub in mapping_df.groupby(["pathway_id", "pathway_name"]):
        matched = []
        for _, row in sub.iterrows():
            nm = row["norm_name"]
            if nm in dataset_lookup:
                matched.append(dataset_lookup[nm])
        out[(pid, pname)] = sorted(set(matched))
    return out


def mask_metabolites(pathway_metabs: list[str], mask_rate: float, rng: random.Random):
    n = len(pathway_metabs)
    k = max(1, int(math.ceil(mask_rate * n)))
    k = min(k, n)
    masked = rng.sample(pathway_metabs, k=k)
    observed = [m for m in pathway_metabs if m not in masked]
    return observed, masked


# =========================
# Scoring functions
# =========================
def metabolite_gain_variance(df: pd.DataFrame, metabolite: str) -> float:
    vals = pd.to_numeric(df[metabolite], errors="coerce")
    score = vals.var(skipna=True)
    if pd.isna(score):
        return 0.0
    return float(score)

def metabolite_gain_proposed(df: pd.DataFrame, observed: list[str], candidate: str) -> float:
    """
    Structural gain based on partial-correlation novelty.

    Goals:
    - reward candidates that connect to the observed set
    - penalize redundancy with already observed metabolites
    - avoid collapsing into plain variance ranking
    """

    vals_c = pd.to_numeric(df[candidate], errors="coerce")
    var_c = vals_c.var(skipna=True)
    if pd.isna(var_c):
        var_c = 0.0

    # No observed metabolites yet -> weak fallback only
    if len(observed) == 0:
        return float(0.01 * var_c)

    sub = df[observed + [candidate]].apply(pd.to_numeric, errors="coerce")

    # Need at least 2 columns and some variance
    if sub.shape[1] < 2:
        return float(0.01 * var_c)

    cov = sub.cov().fillna(0.0).values

    # Small ridge for numerical stability
    ridge = 1e-6 * np.eye(cov.shape[0])
    cov_reg = cov + ridge

    try:
        precision = np.linalg.pinv(cov_reg)
    except Exception:
        return float(0.01 * var_c)

    # Candidate is last column
    j = precision.shape[0] - 1

    partial_corrs = []
    for i in range(j):
        denom = np.sqrt(abs(precision[i, i] * precision[j, j]))
        if denom <= 1e-12:
            continue
        pcorr = -precision[i, j] / denom
        if np.isfinite(pcorr):
            partial_corrs.append(abs(float(pcorr)))

    if len(partial_corrs) == 0:
        return float(0.01 * var_c)

    partial_corrs = np.array(partial_corrs, dtype=float)

    mean_pcorr = float(np.mean(partial_corrs))
    max_pcorr = float(np.max(partial_corrs))
    std_pcorr = float(np.std(partial_corrs))

    # Broad conditional connectivity is good
    broad_support = mean_pcorr

    # Penalize "one dominant edge" behavior
    redundancy_penalty = max_pcorr - mean_pcorr

    # Reward distributed conditional structure
    diversity_reward = std_pcorr

    # Tiny variance tie-breaker only
    score = (
        1.2 * broad_support
        - 1.0 * redundancy_penalty
        + 0.6 * diversity_reward
        + 0.001 * float(var_c)
    )

    return float(score)

def oracle_ranking(df: pd.DataFrame, observed: list[str], masked: list[str]) -> list[tuple[str, float]]:
    gains = [(m, metabolite_gain_proposed(df, observed, m)) for m in masked]
    gains.sort(key=lambda x: x[1], reverse=True)
    return gains


def method_predict(df: pd.DataFrame, pathway_metabs: list[str], observed: list[str], masked: list[str], method: str, rng: random.Random):
    if not masked:
        return None, []

    if method == "random":
        ranked = masked[:]
        rng.shuffle(ranked)
        return ranked[0], ranked

    elif method == "degree":
        # Proxy: rank by total non-missing connectivity to observed set via absolute correlation sum
        scores = []
        for m in masked:
            vals_m = pd.to_numeric(df[m], errors="coerce")
            score = 0.0
            for o in observed:
                vals_o = pd.to_numeric(df[o], errors="coerce")
                corr = vals_m.corr(vals_o)
                if pd.notna(corr):
                    score += abs(float(corr))
            # fallback when observed is empty or all correlations NaN
            if score == 0.0:
                score = float(vals_m.notna().sum())
            scores.append((m, score))
        scores.sort(key=lambda x: x[1], reverse=True)
        ranked = [m for m, _ in scores]
        return ranked[0], ranked

    elif method == "imputation":
        # Placeholder imputation proxy: per-metabolite variance
        scores = []
        for m in masked:
            score = metabolite_gain_variance(df, m)
            scores.append((m, score))
        scores.sort(key=lambda x: x[1], reverse=True)
        ranked = [m for m, _ in scores]
        return ranked[0], ranked

    elif method == "proposed":
        scores = []
        for m in masked:
            score = metabolite_gain_proposed(df, observed, m)
            scores.append((m, score))
        scores.sort(key=lambda x: x[1], reverse=True)
        ranked = [m for m, _ in scores]
        return ranked[0], ranked

    else:
        raise ValueError(f"Unknown method: {method}")


# =========================
# Summaries
# =========================
def summarize_results(trial_df: pd.DataFrame):
    regret_summary = (
        trial_df.groupby(
            ["dataset", "pathway_id", "pathway_name", "mask_rate", "method"],
            as_index=False
        )
        .agg(
            n_trials=("trial", "count"),
            mean_regret=("regret", "mean"),
            median_regret=("regret", "median"),
            std_regret=("regret", "std"),
            mean_normalized_regret=("normalized_regret", "mean"),
        )
    )

    topk_summary = (
        trial_df.groupby(
            ["dataset", "pathway_id", "pathway_name", "mask_rate", "method"],
            as_index=False
        )
        .agg(
            n_trials=("trial", "count"),
            top1_success=("is_top1", "mean"),
            top3_success=("is_top3", "mean"),
            top5_success=("is_top5", "mean"),
        )
    )

    return regret_summary, topk_summary


def global_summary_from(regret_summary: pd.DataFrame, topk_summary: pd.DataFrame) -> pd.DataFrame:
    g1 = regret_summary.groupby("method", as_index=False)["mean_normalized_regret"].mean()
    g2 = topk_summary.groupby("method", as_index=False)[["top1_success", "top3_success", "top5_success"]].mean()
    out = g1.merge(g2, on="method", how="outer").sort_values("mean_normalized_regret", ascending=True)
    return out


# =========================
# Main
# =========================
def main():
    random.seed(RANDOM_SEED)
    np.random.seed(RANDOM_SEED)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    mapping_df = load_mapping(MAPPING_FILE)
    retained_df = load_retained(RETAINED_FILE)

    all_trials = []

    for dataset_name, dataset_file in DATASET_FILES.items():
        if not dataset_file.exists():
            print(f"WARNING: missing dataset file for {dataset_name}: {dataset_file}")
            continue

        print(f"\n=== DATASET: {dataset_name} ===")
        data_df = load_dataset(dataset_name, dataset_file)
        dataset_metabs = get_dataset_metabolite_columns(data_df)
        pathway_map = build_pathway_metabolite_map(mapping_df, dataset_metabs)

        retained_sub = retained_df[retained_df["dataset"] == dataset_name].copy()
        if retained_sub.empty:
            print(f"WARNING: no retained pathways for {dataset_name}")
            continue

        for _, row in retained_sub.iterrows():
            pathway_id = row["pathway_id"]
            pathway_name = row["pathway_name"]
            key = (pathway_id, pathway_name)

            pathway_metabs = pathway_map.get(key, [])
            if len(pathway_metabs) < 3:
                continue

            print(f"  Pathway: {pathway_name} ({len(pathway_metabs)} metabolites)")

            for mask_rate in MASK_RATES:
                for trial in range(1, N_TRIALS + 1):
                    rng = random.Random(RANDOM_SEED + hash((dataset_name, pathway_id, mask_rate, trial)) % 10_000_000)

                    observed, masked = mask_metabolites(pathway_metabs, mask_rate, rng)
                    if len(masked) == 0:
                        continue

                    oracle_scores = oracle_ranking(data_df, observed, masked)
                    if not oracle_scores:
                        continue

                    oracle_metabolite, oracle_gain = oracle_scores[0]
                    oracle_ranked = [m for m, _ in oracle_scores]

                    for method in METHODS:
                        predicted_metabolite, ranked = method_predict(
                            data_df,
                            pathway_metabs,
                            observed,
                            masked,
                            method,
                            rng,
                        )

                        if predicted_metabolite is None:
                            continue

                        if method == "proposed":
                            predicted_gain = metabolite_gain_proposed(data_df, observed, predicted_metabolite)
                        else:
                            predicted_gain = metabolite_gain_proposed(data_df, observed, predicted_metabolite)

                        regret = float(oracle_gain - predicted_gain)
                        normalized_regret = regret / (oracle_gain + 1e-12) if oracle_gain > 0 else np.nan

                        top1 = int(predicted_metabolite == oracle_metabolite)
                        top3 = int(predicted_metabolite in oracle_ranked[:3])
                        top5 = int(predicted_metabolite in oracle_ranked[:5])

                        all_trials.append({
                            "dataset": dataset_name,
                            "pathway_id": pathway_id,
                            "pathway_name": pathway_name,
                            "mask_rate": mask_rate,
                            "trial": trial,
                            "method": method,
                            "n_pathway_metabolites": len(pathway_metabs),
                            "n_observed": len(observed),
                            "n_masked": len(masked),
                            "oracle_metabolite": oracle_metabolite,
                            "predicted_metabolite": predicted_metabolite,
                            "oracle_gain": oracle_gain,
                            "predicted_gain": predicted_gain,
                            "regret": regret,
                            "normalized_regret": normalized_regret,
                            "is_top1": top1,
                            "is_top3": top3,
                            "is_top5": top5,
                        })

    if not all_trials:
        raise ValueError("No benchmark trials were generated.")

    trial_df = pd.DataFrame(all_trials)
    regret_summary, topk_summary = summarize_results(trial_df)

    trial_df.to_csv(TRIAL_OUT, index=False)
    regret_summary.to_csv(REGRET_OUT, index=False)
    topk_summary.to_csv(TOPK_OUT, index=False)

    print(f"\nSaved trial-level results to: {TRIAL_OUT}")
    print(f"Saved regret summary to: {REGRET_OUT}")
    print(f"Saved top-k summary to: {TOPK_OUT}")

    print("\n=== REGRET SUMMARY (head) ===")
    print(regret_summary.head().to_string(index=False))

    print("\n=== TOP-K SUMMARY (head) ===")
    print(topk_summary.head().to_string(index=False))

    # Non-trivial subset
    trial_df_nontrivial = trial_df[trial_df["n_masked"] >= MIN_NONTRIVIAL_MASKED].copy()
    if len(trial_df_nontrivial) > 0:
        regret_nontrivial, topk_nontrivial = summarize_results(trial_df_nontrivial)
        global_nontrivial = global_summary_from(regret_nontrivial, topk_nontrivial)

        trial_df_nontrivial.to_csv(TRIAL_OUT_NONTRIVIAL, index=False)
        regret_nontrivial.to_csv(REGRET_OUT_NONTRIVIAL, index=False)
        topk_nontrivial.to_csv(TOPK_OUT_NONTRIVIAL, index=False)
        global_nontrivial.to_csv(GLOBAL_OUT_NONTRIVIAL, index=False)

        print("\n=== GLOBAL SUMMARY (NON-TRIVIAL TRIALS ONLY) ===")
        print(global_nontrivial.round(3).to_string(index=False))
    else:
        print("\nNo non-trivial trials found with current threshold.")


if __name__ == "__main__":
    main()
