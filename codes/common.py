
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

# Color spec (as requested)
OBS_COLOR = "#1f77b4"   # solid blue
MISS_EDGE = "#9e9e9e"   # hollow grey outline
LOD_COLOR = "#d62728"   # dashed red outline
TEXT_COLOR = "#111111"

def setup_mpl():
    plt.rcParams.update({
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "font.size": 10,
        "axes.titlesize": 12,
        "axes.labelsize": 10,
        "legend.fontsize": 9,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "font.family": "DejaVu Sans",
    })

def save_fig(fig, outbase):
    fig.tight_layout()
    fig.savefig(str(outbase) + ".pdf", bbox_inches="tight")
    fig.savefig(str(outbase) + ".svg", bbox_inches="tight")
    plt.close(fig)

def draw_graph(ax, G, pos, node_state, title=None, show_labels=True):
    # edges
    nx.draw_networkx_edges(G, pos, ax=ax, width=1.6, edge_color="#444444", arrows=True, arrowsize=12)

    obs = [n for n in G.nodes if node_state.get(n, {}).get("type") == "obs"]
    miss = [n for n in G.nodes if node_state.get(n, {}).get("type") == "miss"]
    lod = [n for n in G.nodes if node_state.get(n, {}).get("type") == "lod"]

    if obs:
        nx.draw_networkx_nodes(G, pos, nodelist=obs, ax=ax,
                               node_color=OBS_COLOR, edgecolors="#222222",
                               linewidths=1.2, node_size=900)
    if miss:
        nx.draw_networkx_nodes(G, pos, nodelist=miss, ax=ax,
                               node_color="white", edgecolors=MISS_EDGE,
                               linewidths=2.0, node_size=900)
    if lod:
        nx.draw_networkx_nodes(G, pos, nodelist=lod, ax=ax,
                               node_color="white", edgecolors=LOD_COLOR,
                               linewidths=2.0, node_size=900)
        for n in lod:
            x, y = pos[n]
            circ = plt.Circle((x, y), radius=0.06, fill=False,
                              edgecolor=LOD_COLOR, linewidth=2.0, linestyle=(0,(3,2)))
            ax.add_patch(circ)

    if show_labels:
        labels = {n: str(n) for n in G.nodes}
        nx.draw_networkx_labels(G, pos, labels=labels, ax=ax, font_size=9, font_color=TEXT_COLOR)

    if title:
        ax.set_title(title)
    ax.set_axis_off()
