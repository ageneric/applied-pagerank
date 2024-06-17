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

def draw_power_law(sorted_values, figure_mode=0, text_values=None, tex=False):
    """0 - Standard plot; 1 - Truncated plot with x^1.333 line"""
    if text_values is None:
        xl = 'Gene Index'
        yl = 'Gene Differential Expression'
        title = 'Power Law Distribution Plot of Gene Differential Expressions'
    else:
        xl, yl, title = text_values

    plt.figure(figsize=(10, 6))
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel(xl)
    plt.ylabel(yl)
    plt.title(title)

    if figure_mode % 5:
        for i in range(0, 5):
            line = plt.axvline(x=10**i, linestyle='-')
            line.set_color("white")

    if figure_mode % 2:
        sorted_values = sorted_values[:10000]

    plt.plot(sorted_values)

    if figure_mode % 3:
        space = 1.1 ** np.arange(0, 98)
        plt.plot(space, space ** 1.32 * sorted_values[0] * 1.2)
        cont_space = 1.1 ** np.arange(98, 100)
        plt.plot(cont_space, cont_space ** 1.32 * sorted_values[0] * 1.2,
                 linestyle='dotted')

    plt.grid(True, which='both')
    plt.show()

def draw_gene_deg_plot(gene_deg, top_gene_names):
    plt.plot([gene_deg[g] for g in top_gene_names[:1500]])
    plt.show()

def draw_sorted_gene_deg_plot(gene_deg):
    plt.figure(figsize=(10, 6))

    plt.xlabel('Gene Differential Expression')
    plt.ylabel('Gene Index')
    plt.title('Plot of Gene Differential Expressions')
    plt.grid(True, which='both')
    space = 1.1 ** np.arange(0, 98)
    plt.plot(space, (space ** 1.33) / 1000)
    plt.plot(sorted(gene_deg.values()))
    plt.show()

def draw_divergence_plot(pagerank_collection):
    plt.show()

def draw_robustness_region(pagerank_collection):
    plt.show()

def draw_pagerank_hist(pageranks):
    plt.hist(sorted(pageranks, key=lambda x: -x), bins=100)
    plt.xlabel('PageRank')
    plt.ylabel('Frequency Density')
    plt.title('Histogram of Gene Differential Expressions (PageRank < 0.0001)')
    plt.show()



def write_to_graphml(network, pagerank_dict, gene_deg, name='graph'):
    import pathlib

    nx.set_node_attributes(network, gene_deg, name='differential_expression')
    nx.set_node_attributes(network, pagerank_dict, name='pagerank')

    attrs = {}
    ws = nx.get_edge_attributes(network, 'weight')
    ws3 = {k: v*v for k, v in ws.items()}

    for i, edge in enumerate(network.edges):
        attrs[edge] = 1000 * pagerank_dict[edge[0]] * pagerank_dict[edge[1]] * ws3[edge]
    nx.set_edge_attributes(network, attrs, 'pagerank_weight')

    nx.set_edge_attributes(network, ws3, 'rescaled_weight')
    nx.write_graphml(network, pathlib.Path.cwd() / (name + '.graphml'))


def alpha_plot(alphas, values, text_values=None, tex=False):
    if text_values is None:
        xl = 'Alpha'
        yl = 'K-L Divergence with alpha = 0.85'
        title = 'K-L Divergence as Teleportation Probability (alpha) Varies'
    else:
        xl, yl, title = text_values

    plt.figure(figsize=(10, 6))
    plt.ylim([-0.05, 0.5])
    plt.xlabel(xl)
    plt.ylabel(yl)
    plt.title(title)

    plt.plot(alphas, values)
