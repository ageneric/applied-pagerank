import networkx as nx
import numpy as np


def linear_pagerank(G, v=None, alpha=0.85):
    if not v:
        v = [1/G.number_of_nodes() for i in range(G.number_of_nodes)]
    A = nx.adjacency_matrix(G)

    np.diff(A.indptr) == 0  # Find zero rows

    # sort the rows/ relabel nodes

    # calculate P, (row-wise? we want to avoid dividing by 0)

    # implement algorithm 1 from the paper
