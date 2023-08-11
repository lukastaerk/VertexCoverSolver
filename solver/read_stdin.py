import sys
import numpy as np


def read_graph():
    first_line = sys.stdin.readline().strip().split(" ")
    num_vertices = int(first_line[0][1:])
    num_edges = int(first_line[1])
    adjacency_list = [set() for _ in range(num_vertices)]
    edges = np.zeros((num_edges, 2), dtype=int)

    edges_index = 0
    for line in sys.stdin:
        if (
            len(line.strip()) == 0 or line.strip()[0] == "#"
        ):  # skip all empty lines and starting with # includes the first line
            continue
        u, w = map(int, line.split())
        adjacency_list[u - 1].add(w - 1)
        adjacency_list[w - 1].add(u - 1)
        edges[edges_index] = (u - 1, w - 1)
        edges_index += 1
    return adjacency_list, edges


def read_graph_fast():
    first_line = sys.stdin.readline().strip().split(" ")
    num_vertices = int(first_line[0][1:])
    num_edges = int(first_line[1])
    adjacency_list = [set() for _ in range(num_vertices)]
    edges = np.zeros((num_edges, 2), dtype=int)

    edges_index = 0
    data = sys.stdin.read()
    dt = np.dtype([("e1", int), ("e2", int)])
    edges = np.fromstring(data, dtype=int, sep=" ") - 1
    edges = edges.reshape((num_edges, 2))

    for u, w in edges:
        adjacency_list[u].add(w)
        adjacency_list[w].add(u)
    return adjacency_list, edges


if __name__ == "__main__":
    read_graph_fast()
