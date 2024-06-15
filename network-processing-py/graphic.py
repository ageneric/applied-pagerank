import networkx as nx
import numpy as np
from matplotlib import pyplot as plt


plt.style.use('ggplot')


def draw_network_pagerank(network, pagerank_dict, top_gene_names, top_pageranks, top_n=40):
    view_network = nx.subgraph(network, top_gene_names[:top_n])
    scale = 900 / max(0.01, top_pageranks[0])
    sizes = np.array([pagerank_dict[node] for node in view_network.nodes]) * scale
    nx.draw(view_network, node_size=sizes, with_labels=True)
    plt.show()

def draw_power_law(sorted_values, figure_mode=0, tex=False):
    """0 - Standard plot; 1 - Truncated plot with x^1.333 line"""
    if figure_mode % 2:
        sorted_values = sorted_values[:10000]

    plt.figure(figsize=(10, 6))
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('Gene Index')
    plt.ylabel('Gene Differential Expression')
    plt.title('Power Law Distribution Plot of Gene Differential Expressions')

    plt.plot(sorted_values)

    if figure_mode % 3:
        space = 1.1 ** np.arange(0, 98)
        plt.plot(space, space ** 1.25)

    plt.grid(True, which='both', linestyle='--', linewidth=0.25)
    plt.show()

def draw_gene_deg_plot(gene_deg, top_gene_names):
    plt.plot([gene_deg[g] for g in top_gene_names[:1500]])
    plt.show()

def draw_divergence_plot(pagerank_collection):
    plt.show()

def draw_robustness_region(pagerank_collection):
    plt.show()

def write_to_graphml(network, pagerank_dict, gene_deg, name='graph'):
    import pathlib

    nx.set_node_attributes(network, gene_deg, name='differential_expression')
    nx.set_node_attributes(network, pagerank_dict, name='pagerank')

    attrs = []
    for edge in network.edges:
        attrs.append(pagerank_dict, [edge[1]])

    nx.write_graphml(network, pathlib.Path.cwd() / (name + '.graphml'))
