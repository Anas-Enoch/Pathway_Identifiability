
import matplotlib.pyplot as plt
import networkx as nx
from common import setup_mpl, save_fig, draw_graph, OBS_COLOR, MISS_EDGE, LOD_COLOR

def build_ppp_gly_graph():
    G = nx.DiGraph()
    nodes = ["Glc", "G6P", "F6P", "F1,6BP", "G3P", "R5P", "NADPH", "Pyruvate"]
    G.add_nodes_from(nodes)
    edges = [
        ("Glc", "G6P"),
        ("G6P", "F6P"),
        ("F6P", "F1,6BP"),
        ("F1,6BP", "G3P"),
        ("G3P", "Pyruvate"),
        ("G6P", "R5P"),
        ("R5P", "NADPH"),
    ]
    G.add_edges_from(edges)
    pos = {
        "Glc": (0.05, 0.55),
        "G6P": (0.25, 0.55),
        "F6P": (0.45, 0.55),
        "F1,6BP": (0.65, 0.55),
        "G3P": (0.85, 0.55),
        "Pyruvate": (1.05, 0.55),
        "R5P": (0.45, 0.25),
        "NADPH": (0.65, 0.25),
    }
    return G, pos

def main(outbase="../figures/Fig1_problem_setup"):
    setup_mpl()
    G, pos = build_ppp_gly_graph()

    node_state = {n: {"type": "obs"} for n in G.nodes}
    node_state["R5P"]["type"] = "miss"
    node_state["NADPH"]["type"] = "lod"
    node_state["F1,6BP"]["type"] = "miss"

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))

    draw_graph(axes[0], G, pos, node_state, title="(A) Interpretation A: Glycolysis-dominant")
    axes[0].annotate("Ambiguous branch point", xy=pos["G6P"], xytext=(0.12, 0.85),
                     textcoords="axes fraction",
                     arrowprops=dict(arrowstyle="->", lw=1.5),
                     fontsize=10)
    gly_edges = [("G6P","F6P"), ("F6P","F1,6BP"), ("F1,6BP","G3P"), ("G3P","Pyruvate")]
    nx.draw_networkx_edges(G, pos, ax=axes[0], edgelist=gly_edges, width=4.0,
                           edge_color="#2ca02c", arrows=True, arrowsize=14)

    draw_graph(axes[1], G, pos, node_state, title="(B) Interpretation B: PPP-dominant")
    ppp_edges = [("G6P","R5P"), ("R5P","NADPH")]
    nx.draw_networkx_edges(G, pos, ax=axes[1], edgelist=ppp_edges, width=4.0,
                           edge_color="#ff7f0e", arrows=True, arrowsize=14)

    import matplotlib.patches as mpatches
    obs_patch = mpatches.Patch(facecolor=OBS_COLOR, edgecolor="#222222", label="Observed")
    miss_patch = mpatches.Patch(facecolor="white", edgecolor=MISS_EDGE, label="Panel-missing")
    lod_patch = mpatches.Patch(facecolor="white", edgecolor=LOD_COLOR, label="LOD-censored")
    fig.legend(handles=[obs_patch, miss_patch, lod_patch], loc="lower center", ncol=3, frameon=False)

    fig.suptitle("Structural Ambiguity in Metabolic Pathways Due to Incomplete Metabolomics Panels", y=1.02)
    save_fig(fig, outbase)

if __name__ == "__main__":
    main()
