from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching
from graph_tool.topology import max_cardinality_matching
import numpy as np
import sys
from vc_graph import VCGraph


def get_lowerbound_max_match(g) -> int:
    """example lowerbound, basis minimum"""
    max_match = max_cardinality_matching(g, edges=True)  # heuristic=True
    return int(max_match.get_array().sum())


def get_csr_adj_matrix(edges, num_vertices) -> "csr_matrix":
    """
    Create a bipartite graph as sparse adjacency csr_matrix.
    Input: "edges", "num_vertices"
    Output: compressed sparse matrix from scipi.sparse package
    """
    num_edges = len(edges)
    data = np.ones(num_edges * 2, dtype="bool")
    edges_1, edges_2 = edges[:, 0].astype(int), edges[:, 1].astype(int)
    L_edges_1, L_edges_2 = edges_1 + num_vertices, edges_2 + num_vertices

    edges_u, edges_v = np.concatenate((edges_1, edges_2)), np.concatenate((L_edges_2, L_edges_1))
    return csr_matrix(
        (data, (edges_u, edges_v)),
        shape=(num_vertices * 2, num_vertices * 2),
        dtype=np.bool,
    )


def get_lowerbound_lp(edges, num_vertices) -> int:
    """
    Calculate the LP lower bound to the ILP by using maximum bipartite matching. Compare lecture 2 / LP bound.
    """
    adj_matrix_csr = get_csr_adj_matrix(edges, num_vertices)

    maximum_matching = maximum_bipartite_matching(adj_matrix_csr)
    return int(len(np.nonzero(maximum_matching > 0)[0]) / 2)
