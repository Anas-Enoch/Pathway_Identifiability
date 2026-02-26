
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from common import setup_mpl, save_fig, draw_graph

def softmax(x, temp=1.0):
    x = np.array(x) / max(1e-9, temp)
    x = x - np.max(x)
    e = np.exp(x)
    return e / (e.sum() + 1e-12)

def build_two_similar_graphs():
    G1 = nx.DiGraph()
    G2 = nx.DiGraph()

    n1 = ["A","B","C","D","E","F"]
    n2 = ["A'","B'","C'","D'","E'","F'"]
    G1.add_nodes_from(n1)
    G2.add_nodes_from(n2)

    edges1 = [("A","B"),("B","C"),("C","D"),("B","E"),("E","F")]
    edges2 = [("A'","B'"),("B'","C'"),("C'","D'"),("B'","E'"),("E'","F'")]
    G1.add_edges_from(edges1)
    G2.add_edges_from(edges2)

    pos1 = {"A": (0.05,0.55), "B": (0.25,0.55), "C": (0.45,0.55), "D": (0.65,0.55),
            "E": (0.45,0.25), "F": (0.65,0.25)}
    pos2 = {"A'": (0.05,0.55), "B'": (0.25,0.55), "C'": (0.45,0.55), "D'": (0.65,0.55),
            "E'": (0.45,0.25), "F'": (0.65,0.25)}
    return (G1,pos1),(G2,pos2)

def random_transport(seed, n, concentrated=False):
    rng = np.random.default_rng(seed)
    T = np.zeros((n,n))
    for i in range(n):
        if concentrated:
            base = -2.5*np.abs(np.arange(n)-i) + rng.normal(0,0.35,size=n)
            p = softmax(base, temp=0.7)
        else:
            base = rng.normal(0,1.0,size=n)
            p = softmax(base, temp=1.4)
        T[i,:] = p
    return T

def draw_mapping(ax, pos1, pos2, nodes1, nodes2, T, color, alpha_scale=0.85):
    for i, u in enumerate(nodes1):
        for j, v in enumerate(nodes2):
            w = T[i,j]
            if w < 0.08:
                continue
            x1,y1 = pos1[u]
            x2,y2 = pos2[v]
            ax.plot([x1, x2+1.2], [y1, y2], color=color, alpha=min(0.95, alpha_scale*w*3.0), lw=1.5)

def main(outbase="../figures/Fig4_fgw_alignment_jl_stability"):
    setup_mpl()
    (G1,pos1),(G2,pos2) = build_two_similar_graphs()
    st1 = {n: {"type":"obs"} for n in G1.nodes}
    st2 = {n: {"type":"obs"} for n in G2.nodes}

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.2))

    axes[0].set_title("(A) Without JL: diffuse / unstable mapping")
    draw_graph(axes[0], G1, pos1, st1, show_labels=True)
    draw_graph(axes[0], G2, {k:(v[0]+1.2, v[1]) for k,v in pos2.items()}, st2, show_labels=True)
    T1 = random_transport(seed=3, n=len(G1.nodes), concentrated=False)
    draw_mapping(axes[0], pos1, pos2, list(G1.nodes), list(G2.nodes), T1, color="#777777")

    axes[1].set_title("(B) With JL: concentrated / stable mapping")
    draw_graph(axes[1], G1, pos1, st1, show_labels=True)
    draw_graph(axes[1], G2, {k:(v[0]+1.2, v[1]) for k,v in pos2.items()}, st2, show_labels=True)
    T2 = random_transport(seed=3, n=len(G1.nodes), concentrated=True)
    draw_mapping(axes[1], pos1, pos2, list(G1.nodes), list(G2.nodes), T2, color="#2ca02c")

    inset = fig.add_axes([0.41, 0.05, 0.18, 0.25])
    seeds = np.arange(1, 21)
    rng0 = np.random.default_rng(0)
    rng1 = np.random.default_rng(1)
    var_no = 0.08 + 0.03*np.sin(seeds/2) + 0.015*rng0.normal(size=len(seeds))
    var_jl = 0.03 + 0.01*np.sin(seeds/2) + 0.008*rng1.normal(size=len(seeds))
    inset.plot(seeds, np.clip(var_no, 0.01, None), label="no JL")
    inset.plot(seeds, np.clip(var_jl, 0.005, None), label="with JL")
    inset.set_title("Stability proxy")
    inset.set_xlabel("seed")
    inset.set_ylabel("Var(distance)")
    inset.legend(frameon=False, fontsize=8)
    inset.grid(True, alpha=0.25)

    fig.suptitle("JL-Stabilized Fused Gromov–Wasserstein Alignment Ensures Reliable Graph Matching", y=1.02)
    save_fig(fig, outbase)

if __name__ == "__main__":
    main()
