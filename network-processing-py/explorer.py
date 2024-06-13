"""Utility script to explore a pre-built network.

Please don't commit changes to this file!

Instructions:
    Comment out either Option 1 or Option 2
    Change the paths as required
    Execute this script in interactive mode with the Python Console
        e.g. run with IDLE
"""

import pandas as pd
import pickle
import networkx as nx

from gene_data import get_gene_differential_expressions
from weighting import WeightVectorMethod, expressions, TF, TARGET
from main import generate_networkx_graph
import pagerank
import scipy

def get_df():
    # df.to_csv('../output/df_TFLink_STRING.csv')
    df = pd.read_csv('../output/df_TFLink_STRING.csv')
    return df

def get_network():
    # pickle.dump(network, open('../output/STRING_network.pickle', 'wb'))
    network = pickle.load(open('../output/product_STRING_network.pickle', 'rb'))
    return network

def compute_kl_divergence(p, q):
    """p, q: PageRank vectors."""
    scipy.special.rel_entr()

def compute_difference(order, other_order):
    """order, other_order: String lists."""

def compute_pagerank_varying_alpha(P):
    pagerank_collection = {
        'pagerank': [],
        'pagerank_dict': [],
        'genes': [],
        'top_genes': []
    }
    for alpha in 0.5, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0:
        pagerank = pagerank.linear_pagerank(matrix_P)
        pagerank_dict = {k: v for k, v in zip(network.nodes, pagerank)}

        genes = [(item[0], item[1]) for item in zip(network.nodes, pagerank)]
        top_genes = sorted(genes, key=lambda x: -x[1])


        pagerank_collection['pagerank'].append(pagerank)
        pagerank_collection['genes'].append(genes)


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
    pagerank = pagerank.linear_pagerank(matrix_P)
    pagerank_dict = {k: v for k, v in zip(network.nodes, pagerank)}

    genes = [(item[0], item[1]) for item in zip(network.nodes, pagerank)]
    top_genes = sorted(genes, key=lambda x: -x[1])

    print(f'{top_genes[:25]=}')
