import networkx as nx
import numpy as np


P = np.array([[0.7, 0.7, 0. , 0.6],
              [0. , 0. , 0. , 0. ],
              [0.3, 0.3, 0. , 0.4],
              [0. , 0. , 0. , 0. ]]).transpose()


def permute_dangling_rows(P, v):
    # Find zero & non-zero rows
    non_dangling_rows = np.any(P, axis=1)  # 1d boolean array

    # Sort the rows / relabel nodes - row-wise
    P_1 = P[non_dangling_rows]

    # Column-wise
    P_11 = P_1[:, non_dangling_rows]
    P_12 = P_1[:, ~non_dangling_rows]

    # Permute v
    v_1 = v[non_dangling_rows]
    v_2 = v[~non_dangling_rows]
    return P_11, P_12, v_1, v_2, non_dangling_rows


def linear_pagerank(P, v=None, alpha=0.85):
    node_count = P.shape[0]

    # Set the default choice of v - uniform probability per outbound edge
    if not v:
        v = np.array([1/node_count for _ in range(node_count)])

    P_11, P_12, v_1, v_2, non_dangling_rows = permute_dangling_rows(P, v)
    permutation = np.arange(node_count)[non_dangling_rows]
    counter_permutation = np.arange(node_count)[~non_dangling_rows]

    # Set up the system of linear equations for PageRank
    I = np.identity(P_11.shape[0])
    X = (I - alpha*P_11)
    # Solving pi_1^T X = v_1^T <= or equivalently => X^T pi_1 = v_1
    pi_1 = np.linalg.solve(X.transpose(), v_1).transpose()

    pi_2 = alpha * pi_1.transpose() @ P_12 + v_2.transpose()

    # Invert the permutation
    ret = np.empty(node_count)
    ret[permutation] = pi_1
    ret[counter_permutation] = pi_2

    # Normalise using the 1-norm as standard
    return ret / np.linalg.norm(ret, 1)
