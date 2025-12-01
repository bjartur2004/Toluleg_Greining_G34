import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import networkx as nx
import numpy as np


def plot_pressure_flow_graph(edges, pressures, positions,
                             node_cmap="viridis", arrow_scale=5):

    if node_cmap == "bright_plasma":
        n_colors = 256
        key_colors = [
            (0.0, (0.15, 0.35, 1.00)),    
            (0.30,(0.75, 0.10, 0.70)),  
            (0.75,(0.95, 0.50, 0.05)),   
            (1.0, (1.00, 0.90, 0.20))      
        ]
    node_cmap = mcolors.LinearSegmentedColormap.from_list("bright_plasma", key_colors, N=n_colors)

    
    G = nx.DiGraph()
    for u, v, flow in edges:
        G.add_edge(u, v, flow=flow)

    pressure_values = np.array([pressures[n] for n in G.nodes()])
    norm = plt.Normalize(vmin=pressure_values.min(),
                         vmax=pressure_values.max())

    fig, ax = plt.subplots(figsize=(10, 8))

    # Draw nodes using the colormap object
    nodes = nx.draw_networkx_nodes(
        G,
        pos=positions,
        node_color=pressure_values,
        cmap=node_cmap,
        node_size=1800,
        edgecolors="black",
        ax=ax
    )

    # Draw node labels
    labels = {n: f"{pressures[n]:.1f}" for n in G.nodes()}
    nx.draw_networkx_labels(G, positions, labels, ax=ax, font_size=10)

    # Draw edges with flow arrows and labels
    flows = [abs(G[u][v]["flow"]) for u, v in G.edges()]
    max_flow = max(flows) if flows else 1.0

    for (u, v, flow) in edges:
        alpha = 0.2 + 0.8 * abs(flow) / max_flow
        start, end = (positions[u], positions[v]) if flow >= 0 else (positions[v], positions[u])

        ax.annotate(
            "",
            xy=end,
            xytext=start,
            arrowprops=dict(
                arrowstyle="->",
                linewidth=4,
                color="black",
                alpha=alpha,
                shrinkA=20,
                shrinkB=20
            )
        )

        mid = ((positions[u][0]+positions[v][0])/2 - 0.06,
               (positions[u][1]+positions[v][1])/2 - 0.01)
        ax.text(mid[0], mid[1], f"{flow:.3f}",
                fontsize=9,
                bbox=dict(facecolor="white", alpha=0.8, edgecolor="none"))

    # Add colorbar with proper colormap
    sm = plt.cm.ScalarMappable(cmap=node_cmap, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, ax=ax, label="Þrýstingur [MPa]")

    ax.set_aspect("equal")
    ax.set_title("Nóðu þrýstingar [MPa] og Flæði [m^3/s]")
    ax.margins(0.15)
    plt.tight_layout()
    plt.show()

# prufu setup
if __name__ == "__main__":
    edges = [
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

    plot_pressure_flow_graph(edges, pressures, positions, node_cmap="bright_plasma")
