#!/bin/bash
set -e
source .venv/bin/activate
python3 scripts/make_fig1_problem_setup.py
python3 scripts/make_fig3_distributional_embeddings.py
python3 scripts/make_fig4_fgw_jl_stability_placeholder.py
python3 scripts/make_fig6_measurement_impact_placeholder.py
python3 scripts/make_fig8_regret_results_placeholder.py
python3 scripts/make_fig9_generalization_grid_placeholder.py
echo "Done. Outputs in figures/"
