from vc_graph import VCGraph
from branch_and_bound import update_max_deg
from utils import dom_helper
import numpy as np
from hopcroft_karp import hopcroft_karp
from collections import defaultdict

RED_MATRIX = np.array(
    [
        (True, False, True, False, False, False, False, False, False),
        (True, False, True, True, True, False, False, False, False),
        (True, True, True, True, True, False, False, False, False),
        (True, False, True, True, False, True, False, False, False),
        (True, True, True, True, True, True, False, False, False),
        (False, False, False, False, False, False, False, False, False),
        (True, False, True, True, False, False, True, False, False),
        (True, False, True, True, False, False, False, True, False),
        (True, True, True, True, True, False, False, True, False),
        (True, True, True, True, True, False, False, False, True),
    ],
    dtype=[
        ("deg_1", "bool"),
        ("deg_2", "bool"),
        ("high_deg", "bool"),
        ("buss", "bool"),
        ("dom", "bool"),
        ("crown", "bool"),
        ("lprule", "bool"),
        ("unconfined", "bool"),
        ("deg3_ind_set", "bool"),
    ],
)
RED_FREQ = np.array(
    [(1, 1, 1, 1, 1), (1, 1, 1, 5, 1), (1, 1, 1, 10, 1), (1, 10, 1, 20, 1)],
    dtype=[
        ("deg_1", "int"),
        ("deg_2", "int"),
        ("dom", "int"),
        ("crown", "int"),
        ("deg3_ind_set", "int"),
    ],
)


def is_unconfined(g: VCGraph, v: int) -> bool:
    S = set([v])
    unclear = True  # we need this loop condition in case we have to do the checks
    # multiple times (which happens if we increase S in the loop)
    while unclear:
        unclear = False
        # Find a vertex u in N(S) with |N(u) intersect S| = 1 with minimal |N(u)\N[S]|
        S_neighborhood = set()
        for s in S:
            S_neighborhood.update(g.get_neighbors(s))
        u = -1
        min_num_extra_neighbors = float("inf")
        u_extra_neighbors = set()
        for n in S_neighborhood:
            n_neighbors = g.get_neighbors(n)
            if len(n_neighbors.intersection(S)) != 1:
                continue
            # Check if |N(n)\N[S]| is minimal
            n_extra_neighbors = n_neighbors.difference(S_neighborhood.union(S))
            num_extra_neighbors = len(n_extra_neighbors)
            if num_extra_neighbors < min_num_extra_neighbors:
                min_num_extra_neighbors = num_extra_neighbors
                u = n
                u_extra_neighbors = n_extra_neighbors
        if u == -1:
            return False  # there is no such vertex u
        if len(u_extra_neighbors) == 0:
            return True
        elif len(u_extra_neighbors) == 1:
            w = min(u_extra_neighbors)  # According to StackOverflow, this seems to be
            # a good approach to get the single element of the set
            # in an efficient manner.. sets don't support indexing
            S.add(w)
            unclear = True  # do the loop again

    return False


def unconfined_rule(g: VCGraph, k: int) -> "bool, list(int), int":
    # For every vertex v in g, check if v is unconfined
    # If yes, put v into the vertex cover and remove v from g
    flag = True
    vc_additions = list()
    # We pick the *first found* unconfined vertex v, apply rule
    # and then return. Let perform_reductions handle the repetition.
    found_unconfined = False
    for bag in g.deg_bags[3:]:
        for vertex in bag:
            if is_unconfined(g, vertex):
                vc_additions.append(vertex)
                found_unconfined = True
                break
        if found_unconfined:
            break
    # saves edge removal effort by doing this check here instead of in perform_reduction
    if len(vc_additions) > k:
        flag = False
        return flag, [], 0

    num_revert = g.remove_edges_for_array(vc_additions)
    g.num_hits_by_reduction_rule["unconfined"] += len(vc_additions)
    return flag, vc_additions, num_revert


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
    pair_L, pair_R = hopcroft_karp(L, R, adjlist)

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


# TODO check if V and E is reduced graph or init graph
def buss_rule(g: VCGraph, k: int, max_deg: int) -> bool:
    flag = True
    if max_deg > k:
        raise Exception("apply high degree rule first!")
    #    num_V = g.init_num_vertices - len(g.deg_bags[0])
    num_V = g.get_num_vertices()
    # |V| > k^2 + k or |E|>k^2
    if num_V > k**2 + k or g.num_edges > k**2:
        flag = False
    return flag


def domination_rule(g: VCGraph, k: int, max_deg: int) -> "bool, list(int), int":
    flag = True
    vc_additions = list()
    num_revert = 0

    if max_deg < 3:
        return flag, vc_additions, num_revert

    relevant_vertices = g.recently_updated_vertices
    # g.recently_updated_vertices.clear() # remove all vertices because we have a copy
    while len(relevant_vertices) > 0:  # relevant_vertices: case <= 2 is handled by other reduction rules
        v = relevant_vertices.pop()
        if g.degrees[v] <= 2:
            continue
        tmp_vc_addition = dom_helper(v, g.adj_list)  # find dom vertex to cover
        if len(tmp_vc_addition) > 0:
            num_revert += g.remove_edges_for_array(tmp_vc_addition)
            vc_additions.extend(tmp_vc_addition)
        if len(vc_additions) > k:
            flag = False
            break
    g.num_hits_by_reduction_rule["dom"] += len(vc_additions)
    return flag, vc_additions, num_revert


def deg_1(g: VCGraph, k: int = float("inf")) -> "bool, list(int), int":
    flag = True
    num_revert = 0
    if len(g.deg_bags[1]) == 0:
        return flag, [], num_revert
    deg_1_vertices = g.deg_bags[1]

    deg_1_vertices_neighbors = set.union(*map(g.adj_list.__getitem__, deg_1_vertices))
    vc_additions = list(deg_1_vertices_neighbors.difference(deg_1_vertices))  # contains no vertices of degree 1

    # special case degree 1 remaining set of single edges
    # take only one vertex of each edge
    single_edges = deg_1_vertices.intersection(deg_1_vertices_neighbors)
    lesser_indices = set(
        [v1 if v1 < list(g.adj_list[v1])[0] else list(g.adj_list[v1])[0] for v1 in single_edges]
    )  # filter edges take lesser index
    vc_additions.extend(list(lesser_indices))

    # saves edge removal effort by doing this check here i nstead of in perform_reduction
    if len(vc_additions) > k:
        flag = False
        return flag, [], num_revert

    num_revert = g.remove_edges_for_array(vc_additions)
    g.num_hits_by_reduction_rule["deg_1"] += len(vc_additions)
    return flag, vc_additions, num_revert


def deg_2(g: VCGraph, k: int = float("inf"), max_deg: int = 0) -> "bool, list(int), int, int, int":
    flag = True
    num_merges = 0
    num_revert = 0
    if len(g.deg_bags[2]) == 0:
        return flag, [], num_revert, num_merges, max_deg
    # remove single triangle and rectangles
    vc_additions = list()

    v, x, y = 0, 0, 0
    pop_next_v = g.deg_bags[2].pop
    while len(g.deg_bags[2]) > 0:  # and len(g.deg_bags[1])==0 first things first TODO test this with aba
        if k < num_merges + len(vc_additions):  # check k
            flag = False
            break
        v = pop_next_v()
        (x, y) = g.adj_list[v]
        if y in g.adj_list[x]:  # triangle
            vc_additions.extend([x, y])
            num_revert += g.remove_edges_for_array([x, y])
        else:  # ' MERGE '
            num_merges += 1
            num_revert += g.merge_deg2(v, x, y)
            if g.degrees[v] > max_deg:
                max_deg = g.degrees[v]

    g.num_hits_by_reduction_rule["deg_2"] += num_merges + len(vc_additions)
    return flag, vc_additions, num_revert, num_merges, max_deg


def high_degree(g: VCGraph, k: int, max_deg: int) -> "bool, list(int), int":
    flag = True
    num_revert = 0
    if k < 0:
        flag = False
        return flag, [], num_revert

    high_degree_vertices = list(set.union(*g.deg_bags[k + 1 : max_deg + 1]))
    num_revert = 0
    if k - len(high_degree_vertices) < 0:
        flag = False
        return flag, [], num_revert

    num_revert = g.remove_edges_for_array(high_degree_vertices)

    g.num_hits_by_reduction_rule["high_deg"] += len(high_degree_vertices)
    return flag, high_degree_vertices, num_revert


def find_maximal_matching(g: VCGraph) -> "list, set":
    matched_vertices = set()
    matching = list()
    for vertices in g.deg_bags[-1:0:-1]:
        for v in vertices:
            if v in matched_vertices:
                continue
            for u in g.adj_list[v]:  # intersection would do the job too, TODO possible large neighborhoods speed-up
                if u in matched_vertices:
                    continue
                matching.append((v, u))
                matched_vertices.update({v, u})
                break
    return matching, matched_vertices


# dict list lookup for find_crown
def get_all_neighbors(graph, neighbors):
    return map(lambda key: graph[key], filter(lambda n: n in graph, neighbors))


def find_crown(g: VCGraph) -> "set, set":
    maximal_matching, A_vertices = find_maximal_matching(g)
    B_vertices = set.union(*g.deg_bags[1:])  # all vertices of the graph deg_0 left out
    B_vertices.difference_update(A_vertices)  # B = V \ A
    edges_AB = list()
    B_vertices_list = list(B_vertices)
    B_graph = dict(zip(B_vertices_list, map(g.adj_list.__getitem__, B_vertices_list)))  # B_graph init B -> A
    A_graph = dict()  # A_graph setting A -> B
    A_vertices.clear()
    for b in B_vertices:
        for a in B_graph[b]:
            if a not in A_graph:
                A_graph[a] = set()
            edges_AB.append((a, b))
            A_graph[a].add(b)
            A_vertices.add(a)

    if len(edges_AB) == 0:
        return set(), set()
    edges_AB = np.array(edges_AB, dtype=int)

    adjlist_AB = B_graph.copy()
    adjlist_AB.update(A_graph)
    pair_A, pair_B = hopcroft_karp(A_vertices, B_vertices, adjlist_AB)
    matching = np.array([(k, v) for k, v in pair_B.items() if v is not None])
    vc_candidates_B, vc_candidates_A = matching[:, 0], matching[:, 1]

    B_maching_g = dict(zip(vc_candidates_B, vc_candidates_A))
    A_maching_g = dict(zip(vc_candidates_A, vc_candidates_B))

    for key in vc_candidates_A:
        A_graph[key].remove(A_maching_g[key])

    in_A = A_vertices.difference(vc_candidates_A)
    Z = in_A.copy()
    while len(in_A) != 0:
        in_B = set.union(*list(get_all_neighbors(A_graph, in_A)))
        in_B.difference_update(Z)
        if len(in_B) == 0:
            break
        Z.update(in_B)
        in_A = set(list(get_all_neighbors(B_maching_g, in_B)))
        in_A.difference_update(Z)
        if len(in_A) == 0:
            break
        Z.update(in_A)

    X = A_vertices.difference(Z).union(B_vertices.intersection(Z))
    A_x = X.intersection(A_vertices)
    B_0 = B_vertices.difference(X)
    # print("A_x: ", A_x,"B_x: ", B_0) #  K=(L \ Z) U (R intersec Z).
    return A_x, B_0


def crown_rule(g: VCGraph, k: int) -> "bool, list(int), int":
    flag = True
    num_revert = 0
    A_x, B_0 = find_crown(g)
    len_A_x = len(A_x)
    if k < len_A_x:
        flag = False
        return flag, [], num_revert
    A_x = list(A_x)
    num_revert = g.remove_edges_for_array(A_x)
    g.num_hits_by_reduction_rule["crown"] += len_A_x
    return flag, A_x, num_revert


def packing_reduction(g: VCGraph, k: int) -> "bool, list(int), int":
    flag, num_revert = 1, 0
    neighbors = set()
    not_in_cover = g.Packing.get_reduction_vertices()
    for v in not_in_cover:
        neighbors.update(g.adj_list[v])

    in_cover = list(neighbors)
    if k < len(in_cover):
        flag = 0
        in_cover.clear()
    else:
        num_revert = g.remove_edges_for_array(in_cover)
        g.Packing.reset_reduction()
        g.num_hits_by_reduction_rule["packing"] += len(in_cover)
    return flag, in_cover, num_revert


def deg3_ind_set(g: VCGraph, max_deg: int) -> "bool, list(int), int, int":
    """
    :param g: graph
    :param max_deg: max degree of the graph
    :return: flag, list of vertices to remove, num_revert, max_deg
    """
    flag, num_revert = 1, 0
    checked_vertices = set()
    if len(g.deg_bags[3]) == 0:
        return flag, [], num_revert, max_deg
    while True:
        deg3 = g.get_vertices_by_degree(3)
        deg3_queue = deg3 - checked_vertices
        if len(deg3_queue) == 0:
            break
        v = next(iter(deg3_queue))
        neighbors_v = g.get_neighbors(v)
        a, b, c = [i for i in neighbors_v]
        neighbors_a = g.get_neighbors(a)
        neighbors_b = g.get_neighbors(b)
        # check if abc is independent set: a is not adjacent with b or c AND b is not adjacent with c
        if len(neighbors_a.intersection({b, c})) + len(neighbors_b.intersection({c})) == 0:
            g.num_hits_by_reduction_rule["deg3_ind_set"] += 1
            local_num_revert, local_max_deg = g.merge_deg3(v)
            num_revert += local_num_revert
            max_deg = max(max_deg, local_max_deg)
        checked_vertices.update({v})
    return flag, [], num_revert, max_deg


def perform_reduction(
    g: VCGraph,
    k: int,
    max_deg: int,
    red_grp: int = 2,
    red_freq: int = 0,
    rec_steps: int = -1,
    preprocessing=False,
) -> "bool, list(int), int, int":
    # CASE 1
    flag = True
    vc = list()
    num_revert = 0
    while True:
        if g.Packing and g.Packing.packing_is_violated():
            flag = 0
            break
        if g.Packing and g.Packing.has_reduction():
            flag, vc_re, re = packing_reduction(g, k)
            vc.extend(vc_re)
            num_revert += re
            k -= len(vc_re)
            if not flag:
                break  # vc and re might be non-empty
            if len(vc_re) > 0:
                continue
        if (
            RED_MATRIX[red_grp]["deg_1"] and (rec_steps + 1) % RED_FREQ[red_freq]["deg_1"] == 0
        ):  # TODO possible speed-up
            flag, vc_re, re = deg_1(g, k)
            vc.extend(vc_re)
            num_revert += re
            k -= len(vc_re)
            if not flag:
                break  # vc and re might be non-empty
            if len(vc_re) > 0:
                continue
        if RED_MATRIX[red_grp]["deg_2"] and (rec_steps + 1) % RED_FREQ[red_freq]["deg_2"] == 0:  # deg_2 rule TODO
            flag, vc_re, re, num_merges, max_deg = deg_2(g, k, max_deg)
            vc.extend(vc_re)
            num_revert += re
            k -= len(vc_re) + num_merges
            if not flag:
                break
            if len(vc_re) + num_merges > 0:
                continue
        max_deg = update_max_deg(g, max_deg)
        if RED_MATRIX[red_grp]["high_deg"] and k < max_deg and not preprocessing:
            flag, vc_re, re = high_degree(g, k, max_deg)
            vc.extend(vc_re)
            num_revert += re
            k -= len(vc_re)
            if not flag:
                break
            if len(vc_re) > 0:
                continue
        max_deg = update_max_deg(g, max_deg)
        if RED_MATRIX[red_grp]["buss"] and not preprocessing:
            flag = buss_rule(g, k, max_deg)
            if not flag:
                break
        if RED_MATRIX[red_grp]["dom"] and (rec_steps + 1) % RED_FREQ[red_freq]["dom"] == 0:
            flag, vc_re, re = domination_rule(g, k, max_deg)
            vc.extend(vc_re)
            num_revert += re
            k -= len(vc_re)
            if not flag:
                break
            if len(vc_re) > 0:
                continue
        if RED_MATRIX[red_grp]["crown"] and (rec_steps + 1) % RED_FREQ[red_freq]["crown"] == 0:
            flag, vc_re, re = crown_rule(g, k)
            vc.extend(vc_re)
            num_revert += re
            k -= len(vc_re)
            if not flag:
                break
            if len(vc_re) > 0:
                continue
        if RED_MATRIX[red_grp]["deg3_ind_set"] and (rec_steps + 1) % RED_FREQ[red_freq]["deg3_ind_set"] == 0:
            flag, vc_re, re, max_deg = deg3_ind_set(g, max_deg)
            num_revert += re
            if not flag:
                break
        if RED_MATRIX[red_grp]["lprule"]:
            flag, vc_re, re = simple_lp_rule(g, k)
            vc.extend(vc_re)
            num_revert += re
            k -= len(vc_re)
            if not flag:
                break
            if len(vc_re) > 0:
                continue
        if RED_MATRIX[red_grp]["unconfined"]:
            flag, vc_re, re = unconfined_rule(g, k)
            vc.extend(vc_re)
            num_revert += re
            k -= len(vc_re)
            if not flag:
                break
            if len(vc_re) > 0:
                continue
        break

    return flag, vc, num_revert, k, max_deg
