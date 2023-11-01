import sys
import numpy as np


def read_graph_fast():
    first_line = sys.stdin.readline().strip().split(" ")
    # can look like this: #6160 40207
    # or this: p td 51795 63334 
    n_m = []
    while first_line:
        w = first_line.pop(0)
        # remove # from w 
        if w[0] == "#":
            w = w[1:]
        # check if w is a number
        if w.isdigit():
            n_m.append(int(w))
            

    num_vertices = n_m[0]
    num_edges = n_m[1]
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
