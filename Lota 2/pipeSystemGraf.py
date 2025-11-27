import matplotlib.pyplot as plt
import networkx as nx
import numpy as np


def plot_pressure_flow_graph(edges, pressures, positions,
                             node_cmap="viridis", arrow_scale=5):
    """
    edges: list of (u, v, flow_value)
    pressures: dict {node: pressure}
    positions: dict {node: (x, y)}
    """

    G = nx.DiGraph()
    for u, v, flow in edges:
        G.add_edge(u, v, flow=flow)

    pressure_values = np.array([pressures[n] for n in G.nodes()])
    norm = plt.Normalize(vmin=pressure_values.min(),
                         vmax=pressure_values.max())

    fig, ax = plt.subplots(figsize=(10, 8))

    # Draw nodes
    nodes = nx.draw_networkx_nodes(
        G,
        pos=positions,
        node_color=pressure_values,
        cmap=node_cmap,
        node_size=1200,
        edgecolors="black",
        ax=ax
    )

    # Draw node labels
    labels = {n: f"{pressures[n]:.2f}" for n in G.nodes()}
    nx.draw_networkx_labels(G, positions, labels, ax=ax, font_size=10)

    # Draw arrows using annotate() manually
    flows = [abs(G[u][v]["flow"]) for u, v in G.edges()]
    max_flow = max(flows) if flows else 1.0

    for (u, v, flow) in edges:
        alpha = 0.2 + 0.8 * abs(flow) / max_flow

        if flow >= 0:
            start = positions[u]
            end   = positions[v]
        else:
            start = positions[v]
            end   = positions[u]

        ax.annotate(
            "",
            xy=end,
            xytext=start,
            arrowprops=dict(
                arrowstyle="->",
                linewidth=4,
                color="black",
                alpha=alpha,
                shrinkA=15,
                shrinkB=15
            )
        )

        # Place flow label
        mid = (
            ((positions[u][0] + positions[v][0]) / 2)+0.1,
            ((positions[u][1] + positions[v][1]) / 2)+0.05
        )

        ax.text(
            mid[0], mid[1],
            f"{flow:.2f}",
            fontsize=9,
            bbox=dict(facecolor="white", alpha=0.0, edgecolor="none")
        )

    # Add colorbar properly
    sm = plt.cm.ScalarMappable(cmap=node_cmap, norm=norm)
    sm.set_array([])

    fig.colorbar(sm, ax=ax, label="Pressure")

    ax.set_aspect("equal")
    ax.set_title("Pipe Network: Node Pressure + Flow Visualization")
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    edges = [ # directon defines the positive flow direction. arrow will be flipped if negative
        ("A", "B", 10),
        ("B", "C", 4),
        ("A", "C", 2),
        ("C", "D", -7)
    ]

    pressures = {
        "A": 100,
        "B": 90,
        "C": 70,
        "D": 50
    }

    positions = {
        "A": (0, 0),
        "B": (1, 1),
        "C": (2, 0),
        "D": (3, 1)
    }

    plot_pressure_flow_graph(edges, pressures, positions)
