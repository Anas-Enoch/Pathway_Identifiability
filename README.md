
# Starter Figure Pack (Placeholder Figures)

This pack generates placeholder versions of:
- Fig 1, Fig 3, Fig 4, Fig 6, Fig 8, Fig 9
as both **PDF** and **SVG**.

## Run locally
Install:
- numpy
- matplotlib
- networkx

Then:
```bash
python scripts/make_fig1_problem_setup.py
python scripts/make_fig3_distributional_embeddings.py
python scripts/make_fig4_fgw_jl_stability_placeholder.py
python scripts/make_fig6_measurement_impact_placeholder.py
python scripts/make_fig8_regret_results_placeholder.py
python scripts/make_fig9_generalization_grid_placeholder.py
```

Outputs appear in `figures/`.

## Color spec
- Observed: solid blue nodes
- Panel-missing: hollow grey nodes
- LOD-censored: dashed red nodes
