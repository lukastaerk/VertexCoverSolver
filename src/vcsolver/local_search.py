import numpy as np
from vc_graph import VCGraph


def do_local_search(g: VCGraph, solution: np.ndarray):
    improvement = 0
    # Pick those vertices as candidates which are
    # part of the solution and whose degree is > 0 after preprocessing.
    # The solution can be modified for these vertices.
    candidates = np.where(solution & (g.degrees > 0))[0]
    candi = set()
    for v in candidates:
        v_neighbors = list(g.get_neighbors(v))

        # Get all neighbors of v which are not in the solution
        not_vc_indices = np.where(solution[v_neighbors] == False)[0]
        if len(not_vc_indices) == 0:
            # all neighbors of v are in the VC, so we can safely remove v from the VC
            solution[v] = False
            improvement += 1
        if len(not_vc_indices) == 1:
            # there is exactly one neighbor of v that is not in the VC, add that neighbor to the VC
            # and remove v
            u = v_neighbors[not_vc_indices[0]]
            solution[v] = False
            solution[u] = True
            candi.update(g.get_neighbors(u))
        # if there are more than one neighbors not in the VC, we cannot remove v

    # try to reduce the solution even further by considering the neighbors of the vertices u
    # we just added to the solution again and removing these neighbors, if they
    # can be removed safely
    for v in candi:
        if not solution[v]:
            continue
        v_neighbors = list(g.get_neighbors(v))
        # if all neighbors of v are in the VC, we can remove v from the VC
        if np.all(solution[v_neighbors]):
            solution[v] = False
            improvement += 1

    return solution, improvement
