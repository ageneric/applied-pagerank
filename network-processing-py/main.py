import pandas as pd
import networkx as nx
from matplotlib import pyplot as plt

from gene_data import get_gene_differential_expressions
from weighting import WeightMethod, expressions, trrust, get_STRING_subset
from pagerank import linear_pagerank
from graphic import draw_network_pagerank


def get_gene_directed_neighbours(gene_name, includes):
    return trrust[(trrust['name1'] == gene_name) & (trrust['name2'].isin(includes))]

def write_to_graphml():
    import pathlib
    nx.write_graphml(network, pathlib.Path.cwd() / 'graph.graphml')


if __name__ == '__main__':
    print('Imported modules.')

    network = nx.DiGraph()
    gene_deg = get_gene_differential_expressions(expressions)
    filter_list = [g for (g, differential) in gene_deg.items() if abs(differential) > 0.0001]
    string_df = get_STRING_subset(filter_list)
    weighting = WeightMethod(gene_deg, string_df).RMS
    print('Imported weightings.')

    for i, gene in enumerate(filter_list):
        neighbours = get_gene_directed_neighbours(gene, filter_list)
        neighbour_genes = neighbours['name2']
        network.add_weighted_edges_from((n_gene, gene, weighting(n_gene, gene))
                                        for n_gene in neighbour_genes)
        if i % 1000 == 999:
            print(f'Generating NetworkX graph: {i} complete')
    print('Generated NetworkX graph.')

    stochastic_network = nx.stochastic_graph(network)

    matrix_P = nx.adjacency_matrix(stochastic_network, stochastic_network.nodes,
                                   weight='weight').todense()

    pagerank = linear_pagerank(matrix_P)
    pagerank_dict = {k: v for k, v in zip(network.nodes, pagerank)}
    nx.set_node_attributes(network, pagerank_dict, name='pagerank')

    genes = [(item[0], item[1]) for item in zip(network.nodes, pagerank)]
    top_genes = sorted(genes, key=lambda x: -x[1])
    top_gene_names = [data[0] for data in top_genes]
    top_pageranks = [data[1] for data in top_genes]
    print('Computed PageRank list.')

    # Visualise 'top 50' genes
    draw_network_pagerank(network, top_gene_names, top_pageranks)
    plt.show()

