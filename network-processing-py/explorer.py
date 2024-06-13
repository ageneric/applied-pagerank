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
import rbo

from gene_data import get_gene_differential_expressions
from weighting import WeightVectorMethod, expressions, TF, TARGET
from main import generate_networkx_graph
from pagerank import linear_pagerank, get_personalisation_vector_by_deg
import scipy

def get_df():
    # df.to_csv('../output/df_TFLink_STRING.csv')
    df = pd.read_csv('../output/df_TFLink_STRING.csv')
    return df

def get_network():
    # pickle.dump(network, open('../output/STRING_network.pickle', 'wb'))
    network = pickle.load(open('../output/STRING_network.pickle', 'rb'))
    return network

def compute_kl_divergence(p, q):
    """p, q: PageRank vectors."""
    return sum(scipy.special.rel_entr(p, q))

def compute_rbo(order, other_order):
    """order, other_order: lists."""
    return rbo.RankingSimilarity(order, other_order).rbo()

def get_pagerank_statistic(P, **kwargs):
    pagerank = linear_pagerank(P, **kwargs)
    pagerank_dict = {k: v for k, v in zip(network.nodes, pagerank)}

    genes = [(item[0], item[1]) for item in zip(network.nodes, pagerank)]
    top_genes = sorted(genes, key=lambda x: -x[1])
    return pagerank, pagerank_dict, genes, top_genes

def get_pagerank_alpha_difference(P):
    pagerank_collection = {
        'pagerank': [],
        'pagerank_dict': [],
        'genes': [],
        'top_genes': []
    }
    for alpha in 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0:
        pagerank_result = get_pagerank_statistic(P, alpha=alpha)
        # note: this requires the order of pagerank_collection keys
        # and return order of pagerank statistic to be identical
        for key, val in zip(pagerank_collection.keys(), pagerank_result):
            pagerank_collection[key].append(val)

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
    pagerank, pagerank_dict, genes, top_genes = get_pagerank_statistic(matrix_P,
                                                                       alpha=0.85,
                                                                       v=get_personalisation_vector_by_deg(network.nodes, gene_deg))

    print(f'{top_genes[:25]=}')
