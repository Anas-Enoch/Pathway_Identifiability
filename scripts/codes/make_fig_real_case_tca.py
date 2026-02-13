import numpy as np
import matplotlib.pyplot as plt


def simulate_tca_case():
    """
    Minimal reproducible TCA/anaplerosis ambiguity scenario.
    """

    # Nodes
    metabolites = [
        "Citrate",
        "α-KG",
        "Succinate",
        "Malate",
        "Glutamine",
        "Pyruvate",
        "Oxaloacetate"
    ]

    # Observed mask (simulate missing anaplerotic nodes)
    observed = {
        "Citrate": True,
        "α-KG": True,
        "Succinate": True,
        "Malate": True,
        "Glutamine": False,
        "Pyruvate": False,
        "Oxaloacetate": False
    }

    # Simulated composite U
    U_before = 0.62

    # Simulated impact ranking (descending)
    impacts = {
        "Glutamine": 0.31,
        "Oxaloacetate": 0.22,
        "Pyruvate": 0.18
    }

    top_metabolite = max(impacts, key=impacts.get)
    U_after = U_before - impacts[top_metabolite]

    return metabolites, observed, impacts, U_before, U_after


def panel_A(ax, metabolites, observed):
    """
    Draw simplified TCA + anaplerotic entry diagram.
    """

    coords = {
        "Citrate": (0.2, 0.6),
        "α-KG": (0.4, 0.75),
        "Succinate": (0.6, 0.6),
        "Malate": (0.4, 0.4),
        "Glutamine": (0.4, 0.95),
        "Pyruvate": (0.05, 0.5),
        "Oxaloacetate": (0.4, 0.25),
    }

    edges = [
        ("Citrate", "α-KG"),
        ("α-KG", "Succinate"),
        ("Succinate", "Malate"),
        ("Malate", "Citrate"),
        ("Glutamine", "α-KG"),
        ("Pyruvate", "Oxaloacetate"),
        ("Oxaloacetate", "Citrate")
    ]

    for u, v in edges:
        ax.plot(
            [coords[u][0], coords[v][0]],
            [coords[u][1], coords[v][1]]
        )

    for m in metabolites:
        x, y = coords[m]
        if observed[m]:
            ax.scatter(x, y)
        else:
            ax.scatter(x, y, marker='s')
        ax.text(x, y + 0.04, m, ha='center', fontsize=8)

    ax.set_title("(A) TCA/anaplerotic subnetwork")
    ax.set_xticks([])
    ax.set_yticks([])


def panel_B(ax, U_before):
    ax.bar(["Observed panel"], [U_before])
    ax.set_ylim(0, 1)
    ax.set_title("(B) Underdetermination score $U_k$")
    ax.set_ylabel("$U_k$")


def panel_C(ax, impacts):
    names = list(impacts.keys())
    values = list(impacts.values())
    ax.bar(names, values)
    ax.set_ylim(0, 0.4)
    ax.set_title("(C) Ranked measurement impact")
    ax.set_ylabel("$\Delta U_k$")
    ax.tick_params(axis='x', rotation=45)


def panel_D(ax, U_before, U_after):
    ax.bar(["Before", "After reveal"], [U_before, U_after])
    ax.set_ylim(0, 1)
    ax.set_title("(D) $U_k$ reduction after reveal")
    ax.set_ylabel("$U_k$")


def main():

    metabolites, observed, impacts, U_before, U_after = simulate_tca_case()

    fig, axes = plt.subplots(2, 2, figsize=(10, 8))

    panel_A(axes[0, 0], metabolites, observed)
    panel_B(axes[0, 1], U_before)
    panel_C(axes[1, 0], impacts)
    panel_D(axes[1, 1], U_before, U_after)

    plt.tight_layout()

    output_path = "results/figures/Fig_real_case_TCA.pdf"
    plt.savefig(output_path)
    plt.close()

    print(f"Saved {output_path}")


if __name__ == "__main__":
    main()
