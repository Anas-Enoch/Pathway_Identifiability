from pathlib import Path
import pandas as pd

from pathway_id.evaluation.benchmark import run_multi_pathway_benchmark
from pathway_id.io.load_metabolights import load_dataset
from pathway_id.io.load_kegg import load_kegg_pathways

def main():
    # ---------------------------
    # Load data
    # ---------------------------
    X, meta = load_dataset()
    pathways = load_kegg_pathways()

    # ---------------------------
    # Run benchmark
    # ---------------------------
    df = run_multi_pathway_benchmark(
        pathways=pathways,
        X=X,
        meta=meta,
        condition_col="condition",
        c1="tumor",
        c2="normal",
        rhos=[0.1, 0.2, 0.3, 0.4, 0.5],
        trials_per_rho=25,
        weights=(1.0, 1.0, 1.0),
    )

    # ---------------------------
    # Save results
    # ---------------------------
    out_dir = Path("results/tables")
    out_dir.mkdir(parents=True, exist_ok=True)

    out_path = out_dir / "multi_pathway_results.csv"
    df.to_csv(out_path, index=False)
    print(f"Saved {out_path.resolve()}")

if __name__ == "__main__":
    main()
