from vc_graph import VCGraph


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
    while g.recently_updated_vertices:
        vertex = g.recently_updated_vertices.pop()
        if is_unconfined(g, vertex):
            vc_additions.append(vertex)
            break
    # saves edge removal effort by doing this check here instead of in perform_reduction
    if len(vc_additions) > k:
        flag = False
        return flag, [], 0

    num_revert = g.remove_edges_for_array(vc_additions)
    g.reduction_hit_counter["unconfined"] += len(vc_additions)
    return flag, vc_additions, num_revert


# TODO check if V and E is reduced graph or init graph
def buss_rule(g: VCGraph, k: int) -> (bool, [], 0):
    flag = True
    if g.max_deg > k:
        raise Exception("Apply high-degree rule first!")
    #    num_V = g.init_num_vertices - len(g.deg_bags[0])
    num_V = g.get_num_vertices()
    # |V| > k^2 + k or |E|>k^2
    if num_V > k**2 + k or g.num_edges > k**2:
        flag = False
    return flag, [], 0


# domination rule
def find_domination_vertex(v: int, adj_list: list):
    v_neighbors = adj_list[v]
    v_neighbors.add(v)
    tmp_vc_addition = list()
    for u in v_neighbors:
        if u == v:
            continue
        u_neighbors = adj_list[u]
        u_neighbors.add(u)
        if v_neighbors.issubset(u_neighbors):
            tmp_vc_addition.append(u)
        u_neighbors.remove(u)

    v_neighbors.remove(v)
    return tmp_vc_addition


def domination_rule(g: VCGraph, k: int) -> "bool, list(int), int":
    flag = True
    vc_additions = list()
    num_revert = 0
    if g.max_deg < 3:
        return flag, vc_additions, num_revert
    relevant_vertices = g.recently_updated_vertices
    while len(relevant_vertices) > 0:  # relevant_vertices: case <= 2 is handled by other reduction rules
        v = relevant_vertices.pop()
        if g.degrees[v] <= 2:
            continue
        tmp_vc_addition = find_domination_vertex(v, g.adj_list)
        if len(tmp_vc_addition) > 0:
            num_revert += g.remove_edges_for_array(tmp_vc_addition)
            vc_additions.extend(tmp_vc_addition)
        if len(vc_additions) > k:
            flag = False
            break
    g.reduction_hit_counter["dom"] += len(vc_additions)
    return flag, vc_additions, num_revert


def deg_1(g: VCGraph, k: int = float("inf")) -> "bool, list(int), int":
    flag = True
    num_revert = 0
    deg_1_vertices = g.deg_bags[1]
    if not deg_1_vertices:
        return flag, [], num_revert

    deg_1_vertices_neighbors = set.union(*map(g.adj_list.__getitem__, deg_1_vertices))
    vc_additions = list(deg_1_vertices_neighbors.difference(deg_1_vertices))  # contains no vertices of degree 1

    # special case degree 1 remaining set of single edges
    # take only one vertex of each edge
    single_edges = deg_1_vertices & deg_1_vertices_neighbors
    lesser_indices = {min(v1, next(iter(g.adj_list[v1]))) for v1 in single_edges}
    vc_additions.extend(lesser_indices)
    # single_edges = deg_1_vertices.intersection(deg_1_vertices_neighbors)
    # lesser_indices = set(
    #    [v1 if v1 < list(g.adj_list[v1])[0] else list(g.adj_list[v1])[0] for v1 in single_edges]
    # )  # filter edges take lesser index
    # vc_additions.extend(list(lesser_indices))

    # saves edge removal effort by doing this check here instead of in perform_reduction
    if len(vc_additions) > k:
        flag = False
        return flag, [], num_revert

    num_revert = g.remove_edges_for_array(vc_additions)
    g.reduction_hit_counter["deg_1"] += len(vc_additions)
    return flag, vc_additions, num_revert


def deg_2(g: VCGraph, k: int = float("inf")) -> "bool, list(int), int, int":
    flag = True
    num_merges = 0
    num_revert = 0
    if len(g.deg_bags[2]) == 0:
        return flag, [], num_revert, num_merges
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
            if g.degrees[v] > g.max_deg:
                g.set_max_deg(g.degrees[v])

    g.reduction_hit_counter["deg_2"] += num_merges + len(vc_additions)
    return flag, vc_additions, num_revert, num_merges


def high_degree(g: VCGraph, k: int) -> "bool, list(int), int":
    flag = True
    num_revert = 0
    if k < 0:
        flag = False
        return flag, [], num_revert

    high_degree_vertices = list(set.union(*g.deg_bags[k + 1 : g.max_deg + 1]))
    if k - len(high_degree_vertices) < 0:
        flag = False
        return flag, [], num_revert

    num_revert = g.remove_edges_for_array(high_degree_vertices)
    g.update_max_deg()
    g.reduction_hit_counter["high_deg"] += len(high_degree_vertices)
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


def packing_reduction(g: VCGraph, k: int) -> "bool, list(int), int":
    flag, num_revert = 1, 0
    if g.Packing.packing_is_violated():
        flag = 0
        return flag, [], num_revert
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
        g.reduction_hit_counter["packing"] += len(in_cover)
    return flag, in_cover, num_revert


def deg3_ind_set(g: VCGraph, k: int) -> "bool, list(int), int, int":
    """
    :param g: graph
    :return: flag, list of vertices to remove, num_revert
    """
    flag, num_revert = 1, 0
    if len(g.deg_bags[3]) == 0:
        return flag, [], num_revert
    deg3_queue = g.get_vertices_by_degree(3).copy()
    while deg3_queue:
        v = deg3_queue.pop()
        neighbors_v = g.get_neighbors(v)
        if len(neighbors_v) != 3:
            continue
        a, b, c = neighbors_v
        neighbors_a = g.get_neighbors(a)
        neighbors_b = g.get_neighbors(b)
        # check if abc is independent set: a is not adjacent with b or c AND b is not adjacent with c
        if len(neighbors_a.intersection({b, c})) + len(neighbors_b.intersection({c})) == 0:
            g.reduction_hit_counter["deg3_ind_set"] += 1
            local_num_revert, local_max_deg = g.merge_deg3(v)
            num_revert += local_num_revert
            g.set_max_deg(max(g.max_deg, local_max_deg))
    return flag, [], num_revert


def perform_reduction(
    g: VCGraph,
    k: int,
    rec_steps: int = -1,
    reduction_mode=0,
) -> "(bool, list(int), int, int)":
    """
    g: graph
    k: current k
    rec_steps: recursion steps
    preprocessing: True if preprocessing is performed

    return: flag, vc, num_revert, k
    """
    # CASE 1
    flag = True
    vc = list()
    num_revert = 0
    # a dictionary of the form {rule_name: (condition function, function)}
    # condition: True if rule is applicable, False otherwise
    # function: function that applies the rule
    rule_dict = {
        "deg_1": (lambda rs: reduction_mode >= 0, deg_1),
        "dom": (lambda rs: reduction_mode == 0 or reduction_mode == 1, domination_rule),
        "packing": (lambda rs: reduction_mode >= 0 and g.Packing, packing_reduction),
        "deg_2": (lambda rs: reduction_mode >= 0, deg_2),
        "high_deg": (lambda rs: reduction_mode != 2 and k < g.max_deg, high_degree),
        "buss": (lambda rs: reduction_mode != 2 and k < g.max_deg, buss_rule),
        "unconfined": (lambda rs: reduction_mode >= 2, unconfined_rule),
        "deg3_ind_set": (lambda rs: reduction_mode == 2 and (rs + 1) % 1 == 0, deg3_ind_set),
    }
    while True:
        applied_rule = False
        for cond_func, rule_func in rule_dict.values():
            if not cond_func(rec_steps):
                continue

            flag, vc_re, re, *extras = rule_func(g, k)
            vc.extend(vc_re)
            num_revert += re
            num_merges = extras[0] if extras else 0
            k -= len(vc_re) + num_merges
            if not flag:
                break  # if the flag is false, it indicates that there is no VC of size k and we break the while-loop after the for-loop.
            if len(vc_re) + num_merges > 0:
                applied_rule = True
                break  # if we reduced something, we want to start over with the first rule

        if not flag or not applied_rule:
            break

    return flag, vc, num_revert, k
