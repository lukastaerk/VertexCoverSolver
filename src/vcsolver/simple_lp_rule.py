from collections import defaultdict
from max_bipartite_matching import maximum_bipartite_matching
from vc_graph import VCGraph


def lp_rule_get_bip_graph(g: VCGraph, current_vertices: list) -> (set, set, dict):
    # We construct the sets L and R of the bipartite graph by considering
    # the current vertices of g.
    # Since we want unique names for the vertices in L and R, we give the vertices
    # in L their original names and for the vertices in R, we use the original names
    # plus an offset.
    # We will fill the two sets in the for loop below.
    l_vertices = set()
    r_vertices = set()
    # this will be an adjacency list (i.e. a dict with vertices as keys and sets (neighbors) as values)
    g_prime_edges = dict()
    offset = len(g.adj_list)  # required for unique naming of the vertices in R
    for v in current_vertices:
        # For every edge {u,v} in the current edge set, add edges
        # {l_u, r_v} and {l_v, r_u} to the edge set of the new bipartite graph
        # Init L set
        g_prime_edges[v] = list([u + offset for u in g.adj_list[v]])
        # Init R set
        g_prime_edges[v + offset] = list(g.adj_list[v])
        l_vertices.add(v)
        r_vertices.add(v + offset)

    return l_vertices, r_vertices, g_prime_edges


# This is a modified BFS that takes as input a bipartite graph (L union R, adjlist)
# where L, R: sets of integers and adjlist: dict with vertices (int) as keys and
# lists of neighbors (int) as values. Furthermore, we require a maximum matching
# (pair_L, pair_R) [may contain unmatched vertices].
#
# Let U be the set of unmatched vertices in L.
# Then this procedure constructs a set Z that consists of vertices that are either
# in U or are connected to U by alternating paths. These are paths that start in
# a vertex in U and then alternate between edges that are not in the matching and edges
# that are in the matching.
def koenig_bfs(L: set, R: set, adjlist: dict, pair_L: dict, pair_R: dict) -> set:
    matched = dict()  # to be able to check for unmatched vertices
    U = set()
    for i in L:
        if pair_L[i] is None:
            matched[i] = True
            U.add(i)
        else:
            matched[i] = True
    for r in R:
        if pair_R[r] is not None:
            matched[r] = True
        else:
            matched[r] = False
    Z = U.copy()
    for u in U:
        discovered = defaultdict(bool)
        discovered[u] = True
        queue = list()
        queue.insert(0, u)
        while len(queue) != 0:
            v = queue.pop()
            if matched[v] and v in R:
                w = pair_R[v]
                if discovered[w]:
                    continue
                discovered[w] = True
                Z.add(w)
                queue.insert(0, w)
            else:
                for w in adjlist[v]:  # for all unmatched edges vw:
                    if discovered[w]:
                        continue
                    discovered[w] = True
                    Z.add(w)
                    queue.insert(0, w)
    return Z


# min_vc is supposed to be a minimal vertex cover for the constructed bipartite graph,
# not the original graph g!
# However, current_vertices is the list of current vertices in the original graph g.
def lp_rule_set_sol_variables(g: VCGraph, current_vertices: list, min_vc: set) -> dict:
    x = dict(zip(current_vertices, [1 / 2] * len(current_vertices)))  # LP relaxation solution
    offset = len(g.adj_list)  # offset is used for naming of R vertices
    for v in current_vertices:
        l_v = v
        r_v = v + offset
        if l_v in min_vc and r_v in min_vc:
            x[v] = 1
        if l_v not in min_vc and r_v not in min_vc:
            x[v] = 0
    return x


def lp_rule_compute_sol(g: VCGraph) -> list:
    current_vertices = g.get_current_vertices()

    # Construct the special bipartite graph g_prime from g
    L, R, adjlist = lp_rule_get_bip_graph(g, current_vertices)

    # HK requires L: set, R: set, adjlist: dict and returns maximum matching pair_L, pair_R = (dict, dict)
    pair_L, pair_R = maximum_bipartite_matching(L, R, adjlist)

    # Given the maximum matching (pair_L, pair_R) use König's Theorem to obtain a minimum
    # vertex cover for g_prime
    Z = koenig_bfs(L, R, adjlist, pair_L, pair_R)
    min_vc = L.difference(Z).union(R.intersection(Z))

    # Now, using min_vc, set the variables x_v of the LP relaxation
    x = lp_rule_set_sol_variables(g, current_vertices, min_vc)
    return x


def simple_lp_rule(g: VCGraph, k: int, improvement_mode: bool = False) -> "bool, list(int), int":
    # Compute the solution (i.e. the values for the variables x_v) of the LP relaxation
    x = lp_rule_compute_sol(g)

    # Given the optimal solution x for the LP relaxation, we define the sets (as arrays)
    # we need for the reduction
    V1 = [key for (key, val) in x.items() if val == 1]  # pick all indices (i.e. vertices) that are in the VC of g_prime
    V0 = [key for (key, val) in x.items() if val == 0]  # pick all vertices that are not in the VC of g_prime
    Vother = [key for (key, val) in x.items() if val == 0.5]  # the remaining vertices

    # We now apply the reduction
    # From the lecture, we know that the vertices in V1 are a subset of the VC of g
    # and vertices in V0 are definitely not in the VC of g
    vc_additions = V1
    flag = True
    #   We do this check here to avoid the edge removal effort, in case this LP rule
    #   application turned out to be fruitless (i.e. we already have a too big VC).
    if len(vc_additions) > k:
        flag = False
        return flag, [], 0

    # TODO The easy & slow kernelization (see lec3, slide 9) might be used here now in order to increase
    # the number of variables set to 1 (this improves the impact of the reduction rule)
    # ...
    if improvement_mode:
        # (trial and error)
        # for each vertex v check whether adding v to the vertex cover
        # increases the LP solution. If not, set x_v = 1
        pass  # TODO

    num_revert = g.remove_edges_for_array(vc_additions)
    g.num_hits_by_reduction_rule["lprule"] += len(vc_additions)
    return flag, vc_additions, num_revert
