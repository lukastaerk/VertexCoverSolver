import numpy as np
from collections import defaultdict


def find_shortest_augmenting_paths(left_vertices, right_vertices, adjacency, left_pairs, right_pairs, distances):
    queue = []
    dummy_vertex = None

    for vertex in left_vertices:
        if left_pairs[vertex] is None:
            distances[vertex] = 0
            queue.append(vertex)
        else:
            distances[vertex] = float("inf")

    distances[dummy_vertex] = float("inf")

    while queue:
        current_vertex = queue.pop(0)
        if distances[current_vertex] < distances[dummy_vertex]:
            for neighbor in adjacency[current_vertex]:
                next_vertex = right_pairs[neighbor]
                if distances[next_vertex] == float("inf"):
                    distances[next_vertex] = distances[current_vertex] + 1
                    queue.append(next_vertex)

    return distances[dummy_vertex] != float("inf")


def traverse_augmenting_path(vertex, adjacency, left_pairs, right_pairs, distances):
    if vertex is not None:
        for neighbor in adjacency[vertex]:
            next_vertex = right_pairs[neighbor]
            if distances[next_vertex] == distances[vertex] + 1:
                if traverse_augmenting_path(right_pairs[neighbor], adjacency, left_pairs, right_pairs, distances):
                    right_pairs[neighbor] = vertex
                    left_pairs[vertex] = neighbor
                    return True
        distances[vertex] = float("inf")
        return False
    return True


def maximum_bipartite_matching(left_vertices, right_vertices, adjacency):
    left_pairs = defaultdict(
        lambda: None,
    )
    right_pairs = defaultdict(lambda: None)
    distances = {}

    while find_shortest_augmenting_paths(left_vertices, right_vertices, adjacency, left_pairs, right_pairs, distances):
        for vertex in left_vertices:
            if left_pairs[vertex] is None:
                traverse_augmenting_path(vertex, adjacency, left_pairs, right_pairs, distances)

    return left_pairs, right_pairs


def main():
    adjacency = {
        0: [5, 6],
        1: [5, 9],
        2: [7, 8],
        3: [5, 9],
        4: [6, 8],
        5: [0, 1, 3],
        6: [0, 4],
        7: [2],
        8: [2, 4],
        9: [1, 3],
    }
    edges = np.array(
        [
            (0, 5),
            (0, 6),
            (1, 5),
            (1, 9),
            (2, 7),
            (2, 8),
            (3, 5),
            (3, 9),
            (4, 6),
            (4, 8),
        ]
    )
    left_set = set(edges[:, 0])
    right_set = set(edges[:, 1])

    left_matching, right_matching = maximum_bipartite_matching(left_set, right_set, adjacency)
    print(dict(left_matching), dict(right_matching))


if __name__ == "__main__":
    main()
