import networkx as nx
import numpy as np
from matplotlib import pyplot as plt
from main import get_gene_differential_expressions, \
    get_gene_directed_neighbours, STRING
from pagerank import linear_pagerank, draw


if __name__ == '__main__':
    network = nx.DiGraph()
    gene_deg = get_gene_differential_expressions()

    threshold = 0.0075
    filtered_list = [x for x in gene_deg if abs(gene_deg[x]) > threshold]

    # draw(network)

    weighting = STRING

    for gene in filtered_list:
        neighbours = get_gene_directed_neighbours(gene, filtered_list)
        neighbour_genes, neighbour_rel = neighbours['V2'], neighbours['V3']
        network.add_weighted_edges_from((n, gene, weighting(n, gene, gene_deg)) for n in neighbour_genes)

    # Create stochastic matrix, manual method
    node_count = network.number_of_nodes()
    P = np.zeros(shape=(node_count, node_count))

    node_names = [node for node in network.nodes]

    for i, node in enumerate(network.nodes):
        edges = network.out_edges(node, data=True)
        weights = np.zeros(node_count)
        for edge in edges:
            j = node_names.index(edge[1])
            weights[j] = np.array([edge[2]['weight']])
        # not a zero row
        if np.any(weights):
            P[i, :] = weights / sum(weights)

    pagerank = linear_pagerank(P)
    print(pagerank)

    genes = [(item[0], item[1]) for item in zip(network.nodes, pagerank)]
    top_genes = sorted(genes, key=lambda x: -x[1])
    gene_names = [data[0] for data in top_genes]

    """
    nx.draw(network, with_labels=True, node_size=pagerank * 10000)
    plt.show()
    """
