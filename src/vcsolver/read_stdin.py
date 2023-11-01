import sys
import numpy as np


def read_graph_fast():
    first_line = sys.stdin.readline().strip().split(" ")
    num_vertices = int(first_line[0][1:])
    num_edges = int(first_line[1])
    adjacency_list = [set() for _ in range(num_vertices)]
    data = sys.stdin.read()
    edges = np.fromstring(data, dtype=np.intc, sep=" ")
    edges = edges - 1
    edges = edges.reshape((num_edges, 2))

    for u, w in edges:
        adjacency_list[u].add(w)
        adjacency_list[w].add(u)
    return adjacency_list, edges


if __name__ == "__main__":
    read_graph_fast()
