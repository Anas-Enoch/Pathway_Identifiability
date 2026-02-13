
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from common import setup_mpl, save_fig, draw_graph

def build_graph():
    G = nx.DiGraph()
    nodes = ["M1","M2","M3","M4","M5","M6","M7"]
    G.add_nodes_from(nodes)
    edges = [("M1","M2"),("M2","M3"),("M3","M4"),("M2","M5"),("M5","M6"),("M6","M7")]
    G.add_edges_from(edges)
    pos = {"M1":(0.05,0.55),"M2":(0.25,0.55),"M3":(0.45,0.55),"M4":(0.65,0.55),
           "M5":(0.45,0.25),"M6":(0.65,0.25),"M7":(0.85,0.25)}
    return G, pos

def main(outbase="../figures/Fig6_measurement_impact"):
    setup_mpl()
    G, pos = build_graph()

    node_state = {n: {"type":"obs"} for n in G.nodes}
    unobs = ["M3","M5","M6"]
    node_state["M3"]["type"] = "miss"
    node_state["M5"]["type"] = "lod"
    node_state["M6"]["type"] = "miss"

    sigma2 = {"M3": 1.6, "M5": 1.2, "M6": 1.4}
    grad = {"M3": 0.9, "M5": 1.5, "M6": 1.1}
    impact = {m: abs(grad[m]) * sigma2[m] for m in unobs}
    m_star = max(impact, key=impact.get)

    fig = plt.figure(figsize=(12, 4.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1.0])

    ax0 = fig.add_subplot(gs[0,0])
    draw_graph(ax0, G, pos, node_state, title="(A) Pathway graph with unmeasured nodes")
    ax0.annotate("candidate measurements", xy=pos["M3"], xytext=(0.05, 0.9),
                 textcoords="axes fraction",
                 arrowprops=dict(arrowstyle="->", lw=1.5),
                 fontsize=10)

    ax1 = fig.add_subplot(gs[0,1])
    labels = unobs
    vals = [impact[m] for m in labels]
    bars = ax1.bar(labels, vals)
    ax1.set_title("(B) Estimated measurement impact")
    ax1.set_ylabel(r"$\widehat{\Delta}_m \mathcal{U}$ (placeholder)")
    ax1.grid(True, axis="y", alpha=0.25)

    for b, m in zip(bars, labels):
        if m == m_star:
            b.set_edgecolor("#111111")
            b.set_linewidth(2.0)
            ax1.annotate("recommended $m^*$", xy=(b.get_x()+b.get_width()/2, b.get_height()),
                         xytext=(0, 10), textcoords="offset points", ha="center",
                         arrowprops=dict(arrowstyle="->", lw=1.4))

    formula = r"$\widehat{\Delta}_m \mathcal{U} = \left\|\nabla_{\mu_m}\mathcal{U}\right\|\cdot\sigma_m^2$"
    ax1.text(0.02, 0.98, formula, transform=ax1.transAxes, va="top",
             bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#444444", alpha=0.95))

    fig.suptitle("Estimating Measurement Impact Without Enumerating Completions", y=1.02)
    save_fig(fig, outbase)

if __name__ == "__main__":
    main()
