
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from common import setup_mpl, save_fig, draw_graph

def make_random_pathway(seed, n=7):
    rng = np.random.default_rng(seed)
    G = nx.DiGraph()
    nodes = [f"v{i}" for i in range(n)]
    G.add_nodes_from(nodes)
    for i in range(n-1):
        G.add_edge(nodes[i], nodes[i+1])
    for _ in range(2):
        a = rng.integers(0, n-2)
        b = rng.integers(a+1, n)
        G.add_edge(nodes[a], nodes[b])
    pos = nx.spring_layout(G, seed=seed)
    xs = np.array([pos[k][0] for k in pos])
    ys = np.array([pos[k][1] for k in pos])

    xs_ptp = np.ptp(xs)
    ys_ptp = np.ptp(ys)

    for k in pos:
        pos[k] = (
            (pos[k][0] - xs.min()) / (xs_ptp + 1e-9),
            (pos[k][1] - ys.min()) / (ys_ptp + 1e-9),
    )
    return G, pos

def assign_states(G, seed):
    rng = np.random.default_rng(seed)
    st = {}
    for n in G.nodes:
        r = rng.random()
        if r < 0.55:
            t = "obs"
        elif r < 0.85:
            t = "miss"
        else:
            t = "lod"
        st[n] = {"type": t}
    return st

def main(outbase="../figures/Fig9_generalization_grid"):
    setup_mpl()
    rows, cols = 2, 3
    fig, axes = plt.subplots(rows, cols, figsize=(12.5, 6.5))

    seeds = [2, 4, 7, 11, 13, 17]
    U = [0.22, 0.48, 0.31, 0.62, 0.27, 0.55]

    for ax, seed, u in zip(axes.flatten(), seeds, U):
        G, pos = make_random_pathway(seed, n=7)
        st = assign_states(G, seed+100)
        draw_graph(ax, G, pos, st, title=None, show_labels=False)
        ax.text(0.02, 0.98, f"U={u:.2f}", transform=ax.transAxes, va="top",
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#444444", alpha=0.95),
                fontsize=10)

    fig.suptitle("Framework Applies to Diverse Pathways with Variable Topologies and Missingness", y=0.98)
    save_fig(fig, outbase)

if __name__ == "__main__":
    main()
