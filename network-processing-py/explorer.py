"""Utility script to explore a pre-built network.

Please don't commit changes to this file!

Instructions:
    Comment out either Option 1 or Option 2
    Change the paths as required
    Execute this script in interactive mode with the Python Console
        e.g. run with IDLE
"""
import numpy
import pandas as pd
import pickle
import networkx as nx
from itertools import combinations

from gene_data import get_gene_differential_expressions
from weighting import WeightVectorMethod, expressions, TF, TARGET
from pagerank import linear_system_pagerank, format_pagerank, get_personalisation_vector_by_deg
from measure import compute_rbo_less_equal_prs, compute_kl_divergence


def get_df():
    # df.to_csv('../output/df_TFLink_STRING.csv')
    df = pd.read_csv('../output/df_TFLink_STRING.csv')
    return df

def get_network():
    # pickle.dump(network, open('../output/STRING_network.pickle', 'wb'))
    network = pickle.load(open('../output/STRING_network.pickle', 'rb'))
    return network


def get_pagerank_alpha_difference(P, v):
    pagerank_collection = {
        'pagerank': [],
        'top_genes': []
    }
    alpha_list = [0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0]
    for alpha in alpha_list:
        pagerank = linear_system_pagerank(P, alpha=alpha, v=v)
        _, top_genes = format_pagerank(pagerank, network)
        # note: this requires the order of pagerank_collection keys
        # and return order of pagerank statistic to be identical
        pagerank_collection['pagerank'].append(pagerank)
        pagerank_collection['top_genes'].append(top_genes)

    rbo, kl = [], []

    for pair in combinations(range(len(alpha_list)), 2):

        rbo.append(compute_rbo_less_equal_prs(pagerank_collection['top_genes'][pair[0]],
                                              pagerank_collection['top_genes'][pair[1]]))
        kl.append(compute_kl_divergence(pagerank_collection['pagerank'][pair[0]],
                                        pagerank_collection['pagerank'][pair[1]]))

    pagerank_collection['rbo'] = rbo
    pagerank_collection['kl'] = kl
    return pagerank_collection


if __name__ == '__main__':
    gene_deg = get_gene_differential_expressions(expressions)
    df = get_df()

    # ===========================================================
    """Option 1: read network"""
    network = get_network()

    """Option 2: compute network (using a different weighting)"""
    # weighting = WeightVectorMethod(gene_deg, df).STRING
    # network = generate_networkx_graph(df, weighting)
    # ===========================================================

    # Compute pagerank using steps in main
    stochastic_network = nx.stochastic_graph(network)
    matrix_P = nx.adjacency_matrix(stochastic_network, stochastic_network.nodes,
                                   weight='weight').todense()
    personalisation = get_personalisation_vector_by_deg(network.nodes, gene_deg)
    pagerank = linear_system_pagerank(matrix_P, alpha=0.85, v=personalisation)
    pagerank_dict, top_genes = format_pagerank(pagerank, network)

    print(f'{top_genes[:25]=}')
