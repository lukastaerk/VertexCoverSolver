import numpy as np
import copy
from utils import (
    set_remove_for_set,
    set_add_for_set,
    set_add_for_array,
    degrees_update_for_array,
    remove_edges_for_array,
    add_edges_for_array,
)
from graph_tool import Graph
from packing_constraints import Packing
from clique_cover import CliqueCover


def generate_deg_bags(degrees, indices):
    max_deg = max(degrees.max(), 3)  # deg_bags at least of size 3 to apply deg3 checks without index error
    degrees, indices = degrees.astype("int32"), indices.astype("int32")
    deg_bags = list([set() for _ in range(max_deg + 1)])
    set_add_for_array(deg_bags, degrees, indices)
    return deg_bags


class VCGraph:
    """
    ATTRIBUTES:
    adj_list: ndarray of set of adjacent vertices (i.e. for every vertex a set of adjacent vertices)
    edges: ndarray of ndarray of type int (representing the edges)
    """

    def __init__(
        self,
        adj_list: list,
        edges: np.ndarray,
        num_vertices: int,
        vc_solution=[],
        degrees=[],
        deg_bags=[],
        gt_graph=Graph(directed=False),
    ):
        self.gt_G = gt_graph  # gt_G is init in preprocessing with edges
        self.edges = edges  # edges are updated in preprocessing
        self.adj_list = adj_list
        if len(vc_solution) == 0:
            vc_solution = np.zeros(num_vertices, dtype=bool)
        self.vc_solution = vc_solution  # -1 unknown, 1 solution, 0 removed array indicating which vertex is part of the finished solution
        self.num_edges = len(edges)
        if len(degrees) == 0:
            degrees = np.array(list(map(len, adj_list)), dtype="int32")
        if len(deg_bags) == 0:
            deg_bags = generate_deg_bags(degrees, np.arange(len(degrees), dtype="int32"))
        self.degrees = degrees
        self.deg_bags = deg_bags
        self.init_num_vertices = num_vertices  # is the actual number of vertices at init time for subgraphs
        self.merge_stack = (
            list()
        )  # deg_2 rule push((u,v,w)) -> merge_deg2(u,v,w)-> u' return revert_merge  --v--u--w--x--y-- last in first out
        self.revert_stack = list()
        self.recently_updated_vertices = (
            set.union(*deg_bags[3:]) if len(deg_bags) > 3 else set()
        )  # for domination_rule only deg(v)>=3 needed
        self.num_hits_by_reduction_rule = {
            "dom": 0,
            "high_deg": 0,
            "deg_1": 0,
            "deg_2": 0,
            "crown": 0,
            "packing": 0,
            "lprule": 0,
            "unconfined": 0,
            "deg3_ind_set": 0,
        }
        # clique cover init
        self.CliqueCover = None
        self.inspect_for_clique_lb = list()
        # Packing Constraint
        self.Packing = None

    def set_packing_constraint(self):
        self.Packing = Packing(len(self.adj_list), self.adj_list, self.vc_solution)

    def set_click_cover(self):
        self.CliqueCover = CliqueCover(len(self.adj_list), self.adj_list, self.vc_solution)
        lb = self.CliqueCover.clique_lb(1, self.deg_bags)

    def copy(self):
        new_g = copy.copy(self)
        new_g.degrees = self.degrees.copy()

        adj_list = [None] * len(self.adj_list)
        for i in range(len(self.adj_list)):
            if len(self.adj_list[i]):
                adj_list[i] = self.adj_list[i].copy()

        deg_bags = [None] * (self.degrees.max() + 1)
        for i in range(len(deg_bags)):
            deg_bags[i] = self.deg_bags[i].copy()

        new_g.adj_list = adj_list
        new_g.deg_bags = deg_bags
        return new_g

    def execute_from_revert_stack(self, number_of_reverts):
        for i in range(number_of_reverts):
            self.revert_stack.pop()()  # take last element and execute

    def merge_deg2(self, v: int, x: int, y: int) -> "int":
        self.merge_stack.append((v, x, y))
        last_len_rs = len(self.revert_stack)
        if self.CliqueCover:
            lb = self.CliqueCover.update()
            # print("# lb in merge: ", lb)
            self.revert_stack.append(lambda: self.CliqueCover.revert_update())

        args = list([v, x, y])
        self.revert_stack.extend([lambda: self.merge_stack.pop(), lambda: self.update_degrees(args)])
        self.revert_stack.append(self.remove_edges(v, self.adj_list[v], update=False))
        x_y_neighbors_without_v = set.union(self.adj_list[x], self.adj_list[y])  # remove old edges

        neighbors_to_update = list(self.adj_list[x].intersection(self.adj_list[y]))
        if len(neighbors_to_update) > 0:
            self.revert_stack.append(lambda: self.update_degrees(neighbors_to_update))
        self.revert_stack.append(self.remove_edges(x, self.adj_list[x], update=False))  # delete x,y remove all edges
        self.revert_stack.append(self.remove_edges(y, self.adj_list[y], update=False))
        self.revert_stack.append(
            self.add_edges(v, x_y_neighbors_without_v, update=False)
        )  # add neighborhood union of x,y to v

        if len(x_y_neighbors_without_v) >= len(self.deg_bags):
            num_extra_deg_bags = (len(x_y_neighbors_without_v) - len(self.deg_bags)) + 1
            self.resize_deg_bags(num_extra_deg_bags)
        self.update_degrees(neighbors_to_update)
        self.update_degrees(args)
        if self.CliqueCover:
            self.CliqueCover.record_vertices(
                [], [v] + list(x_y_neighbors_without_v)
            )  # vc_size -1 so only x is in the solution and y is not
            self.CliqueCover.update_c_lb(x, y)
            lb = self.CliqueCover.update()
            # print("# lb in merge2: ", lb)
            self.revert_stack.append(lambda: self.CliqueCover.undo_update_c_lb(x, y))
            self.revert_stack.append(lambda: self.CliqueCover.revert_update())

        return len(self.revert_stack) - last_len_rs

    def merge_deg3(self, v: int) -> "int, int":
        last_len_rs = len(self.revert_stack)
        neighbors_v = self.get_neighbors(v)
        a, b, c = [i for i in neighbors_v]
        # we have to add edges for e.g. a to N(b). Because N(v) could overlap with N(a) we remove N(a). This way no edges are added  that are already there.
        # This could be an issue when we remove the edges during reverting the changes. E.g. a edge could be removed even though it was there before.
        neighbors_a_without_neighbors_c = self.get_neighbors(a) - self.get_neighbors(c)
        neighbors_b_without_neighbors_a = self.get_neighbors(b) - self.get_neighbors(a)
        neighbors_c_without_neighbors_b = self.get_neighbors(c) - self.get_neighbors(b)
        # add tuple of vertices to revert changes once solution is found and correct solution has to be derived from solution of merged graph
        self.merge_stack.append((v, a, b, c))
        self.revert_stack.append(lambda: self.merge_stack.pop())
        # add update all vertices to merge_stack first, such that this happens last when reverting from revert_stack via .pop()()
        all_vertices_to_be_updated = list(
            set.union(
                {v},
                neighbors_v,
                neighbors_a_without_neighbors_c,
                neighbors_b_without_neighbors_a,
                neighbors_c_without_neighbors_b,
            )
        )
        self.revert_stack.append(lambda: self.update_degrees(all_vertices_to_be_updated))
        # apply changes according to deg3 indep set rule
        self.revert_stack.append(self.remove_edges(v, {a, b, c}, update=False))
        self.revert_stack.append(self.add_edges(a, set({b}), update=False))
        self.revert_stack.append(self.add_edges(b, set({c}), update=False))
        self.revert_stack.append(self.add_edges(a, neighbors_b_without_neighbors_a, update=False))
        self.revert_stack.append(self.add_edges(b, neighbors_c_without_neighbors_b, update=False))
        self.revert_stack.append(self.add_edges(c, neighbors_a_without_neighbors_c, update=False))
        # degree of a, b or c could be higher than size of deg_bags -> resize if necessary
        all_vertices_without_v = set.union(
            {a, b, c},
            neighbors_b_without_neighbors_a,
            neighbors_c_without_neighbors_b,
            neighbors_a_without_neighbors_c,
        )
        all_vertices_to_be_updated = list(all_vertices_without_v)

        degrees = [len(self.adj_list[i]) for i in all_vertices_to_be_updated]
        local_max_degree = max(degrees)
        num_missing_deg_bags = local_max_degree - len(self.deg_bags) + 1
        self.resize_deg_bags(num_missing_deg_bags)
        # update degrees of all vertices involved in the merge
        all_vertices_to_be_updated.append(v)

        self.update_degrees(all_vertices_to_be_updated)
        if self.CliqueCover:
            self.CliqueCover.record_vertices([], all_vertices_to_be_updated)

        return len(self.revert_stack) - last_len_rs, local_max_degree

    def resize_deg_bags(self, num_missing_deg_bags: int):
        extra_deg_bags = [set() for _ in range(num_missing_deg_bags)]
        self.deg_bags.extend(extra_deg_bags)

    def update_solution(self, vertices, flag: bool):
        self.vc_solution[vertices] = flag

    # Returns the current number of vertices of this graph
    # (This number might not be the same as the initial number due to things like
    # preprocessing)
    def get_num_vertices(self) -> int:
        return self.init_num_vertices - len(self.deg_bags[0])

    def get_degree(self, v: int) -> int:
        return self.degrees[v]

    # Returns the current vertices (i.e. vertices *after* preprocessing and *during*
    # branching).
    # CAUTION: Costly operation
    def get_current_vertices(self) -> list:
        return list(set.union(*self.deg_bags[1:]))

    def get_vertices_by_degree(self, deg: int) -> set:
        return self.deg_bags[deg]

    def get_neighbors(self, v: int) -> set:
        return self.adj_list[v]

    def update_num_edges(self, num):
        self.num_edges -= num

    def degrees_decrement(self, vertex: int):
        deg = self.degrees[vertex]
        self.deg_bags[deg].remove(vertex)
        deg -= 1
        self.deg_bags[deg].add(vertex)
        self.degrees[vertex] = deg

    def update_to_zero_degree(self, vertex: int):
        deg = self.degrees[vertex]
        self.deg_bags[deg].remove(vertex)
        deg = 0
        self.deg_bags[deg].add(vertex)
        self.degrees[vertex] = deg

    def update_degrees(self, vertices: list):
        degrees_update_for_array(self.adj_list, self.deg_bags, self.degrees, vertices)

    # removes edges
    def remove_edges_for_array(self, vertices: list, packing_init=False, v=None) -> int:
        last_len_rs = len(self.revert_stack)
        if len(vertices) == 0:
            return 0
        if packing_init:
            self.Packing.init_satellites(v, vertices, self.adj_list)
            self.revert_stack.append(lambda: self.Packing.revert_init_satellites(vertices))
        elif self.Packing:
            self.Packing.update_constraint(vertices)
            self.revert_stack.append(lambda: self.Packing.revert_update_constraint(vertices))

        neighboorhoods = map(lambda i: self.adj_list[i], vertices)
        to_update_deg = set.union(*neighboorhoods)
        self.recently_updated_vertices.update(to_update_deg)
        to_update_deg.update(vertices)
        to_update_deg = list(to_update_deg)

        self.revert_stack.append(lambda: self.update_degrees(to_update_deg))
        old_edges, num_edges = remove_edges_for_array(self.adj_list, vertices)
        self.update_num_edges(num_edges)
        self.revert_stack.append(lambda: self.update_num_edges(-num_edges))
        self.revert_stack.append(lambda: add_edges_for_array(self.adj_list, vertices, old_edges))
        self.update_degrees(to_update_deg)

        self.update_solution(vertices, True)
        self.revert_stack.append(lambda: self.update_solution(vertices, False))

        if self.CliqueCover:
            self.CliqueCover.record_vertices(vertices, to_update_deg)

        return len(self.revert_stack) - last_len_rs

    # update
    def remove_edges(self, v: int, neighbors: set, update: bool = True, packing_init=False) -> "function":
        """
        removes all edges of v
        param: int v
               list neighbors
        return: function lambda -- to revert removing the edges
        """
        if packing_init:
            self.Packing.update_constraint_v(v, neighbors)
        # remove v from adjacency list of its neighbors
        set_remove_for_set(self.adj_list, neighbors, v)

        # clique cover update TODO
        # self.inspect_for_clique_lb.extend(self.adj_list[v])
        # self.clique_size[self.clique[neighbors]] -= 1
        self.adj_list[v] = self.adj_list[v].difference(
            neighbors
        )  # make a new set so the neighbors are not lost for revert
        if update:
            neighbors_list = list(neighbors)
            self.update_degrees(neighbors_list)
            self.update_to_zero_degree(v)
            self.recently_updated_vertices.update(
                neighbors
            )  # use the set of neighbors for update list -> set expensive, set -> list not so much!
            self.update_solution(v, True)
            if self.CliqueCover:
                self.CliqueCover.record_vertices([v], neighbors_list)
        self.num_edges -= len(neighbors)

        return lambda: self.add_edges(v, neighbors, update, packing_init)

    def add_edges(self, v: int, neighbors: set, update: bool = True, packing_init=False) -> "function":
        """add back old edges of v to its neighbors"""
        self.adj_list[v].update(neighbors)
        set_add_for_set(self.adj_list, neighbors, v)
        if packing_init:
            self.Packing.revert_update_constraint_v(v, neighbors)
        if update:
            self.update_degrees(list(neighbors))
            self.update_to_zero_degree(v)
            self.update_solution(v, False)

        self.num_edges += len(neighbors)

        return lambda: self.remove_edges(v, neighbors, update)
