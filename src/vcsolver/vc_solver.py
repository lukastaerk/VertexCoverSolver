import numpy as np
import sys
from lower_bound import get_lowerbound_lp
from vc_graph import VCGraph
from collections import Counter
import reduction
import time
import random
from local_search import do_local_search


class VCSolver:
    def __init__(
        self,
        num_vertices: int,
        time_limit=None,
        recursive_steps=0,
        preprocessing=False,
        greedy=False,
        local_search=False,
        reduction_grouping=2,
        reduction_frequency=0,
        print_lower_bound=False,
        **kwargs,
    ):
        """This is the class container for the Vertex Cover Solver
        method: run"""
        self.recursive_steps = recursive_steps  # to get this information in the results
        self.last_k = 0
        self.greedy = greedy
        (
            self.reduction_grouping,
            self.reduction_frequency,
            self.print_lower_bound,
            self.preprocessing,
            self.local_search,
        ) = (
            reduction_grouping,
            reduction_frequency,
            print_lower_bound,
            preprocessing,
            local_search,
        )
        self.reduction_hit_counter = Counter()
        self.branching_count = 0
        self.random_candidates = list()
        self.last_max_deg = -1
        self.time_limit = time_limit if time_limit else time.time() + 50
        # greedy init
        self.upperbound = num_vertices
        self.pre_solution = np.zeros(num_vertices, dtype=bool)
        self.best_solution = np.zeros(num_vertices, dtype=bool)
        self.best_solution_merge_stack = None
        print(
            "# init: reduction_grouping: %d, reduction_frequency: %d, print_lower_bound: %d, preprocessing: %d, greedy: %d, local_search: %d"
            % (
                reduction_grouping,
                reduction_frequency,
                print_lower_bound,
                preprocessing,
                greedy,
                local_search,
            )
        )

    def run(self, g: VCGraph):
        if self.print_lower_bound:
            self.print_all_lowerbounds(g)
        if self.preprocessing:
            vc_solved = self.run_preprocessing(g)
            self.reduction_hit_counter.update(g.num_hits_by_reduction_rule)
            if vc_solved:
                self.resolve_merged_vertices(g.merge_stack, g.vc_solution)
                return np.nonzero(g.vc_solution)[0]
        if self.greedy:
            self.compute_upperbound(g)
            print("# upperbound: ", self.upperbound)

        g.vc_solution = self.pre_solution.copy()
        self.best_solution_merge_stack = g.merge_stack
        g.set_packing_constraint()
        g.set_click_cover()
        g.set_max_deg(g.degrees.max())
        ub = self.constrained_branching(g, 0, self.upperbound)
        g.merge_stack = self.best_solution_merge_stack
        self.reduction_hit_counter.update(g.num_hits_by_reduction_rule)
        # return solution
        self.resolve_merged_vertices(g.merge_stack, self.best_solution)
        return np.nonzero(self.best_solution)[0]

    def run_preprocessing(self, g: VCGraph) -> bool:
        all_red_rules_activated = self.reduction_grouping  # NOTE this was equal 4 in the solver3
        flag, vc_reduction, _, _ = reduction.perform_reduction(
            g,
            float("inf"),
            all_red_rules_activated,
            self.reduction_frequency,
            preprocessing=True,
        )
        self.pre_solution[:] = g.vc_solution.copy()

        if len(g.deg_bags[0]) == len(g.adj_list):  # all vertices v deg(v) == 0 [vc] is solved
            return True
        return False

    def get_random_vertex(self, g: VCGraph):
        if self.last_max_deg != g.max_deg:
            self.last_max_deg = g.max_deg
            self.random_candidates = list(g.deg_bags[g.max_deg])
            random.shuffle(self.random_candidates)
        v_deg, v = 0, -1
        while g.max_deg != v_deg:
            v = self.random_candidates.pop()
            v_deg = g.degrees[v]
        return v

    def compute_upperbound(self, g: VCGraph):
        num_runs = min(100, max(int(10**6 / g.num_edges), 1))
        print("# number of random runs: %s" % num_runs)
        min_vc = float("inf")
        last_run_time = 0
        for i in range(num_runs):
            start_time = time.time()
            if self.time_limit < start_time + last_run_time:
                break
            g.vc_solution = self.pre_solution.copy()
            self.random_candidates = list()
            self.last_max_deg = -1
            g_copy = g.copy()
            vc_size = self.greedy_heuristic(g_copy, g_copy.degrees.max(), vc_size=0, rnd_vertex=i != 0)
            if self.local_search:  # Improve the current solution using local search
                g.vc_solution, imp = do_local_search(g, g.vc_solution)
                vc_size -= imp
            if vc_size < min_vc:
                self.best_solution = g.vc_solution.copy()
                min_vc = vc_size
            last_run_time = time.time() - start_time

        self.upperbound = np.count_nonzero(self.best_solution)

    def greedy_heuristic(self, g: VCGraph, vc_size: int = 0, rnd_vertex=False) -> np.ndarray:

        while g.max_deg != 0:

            # Iterate over all vertices with degree 1
            # and add one of their neighbors to the vertex cover
            while g.deg_bags[1]:
                # Get a vertex of degree 1, "remove" it from the graph (set degree to 0)
                first = g.deg_bags[1].pop()
                g.degrees[first] = 0
                g.deg_bags[0].add(first)

                # Pick a neighbor (neigh) of above vertex and remove above vertex from
                # the adjacency list of neigh
                neigh = g.adj_list[first].pop()
                g.adj_list[neigh].remove(first)

                # Add neigh to the vertex cover and remove neigh from the graph
                # (i.e. remove it from the adjacency lists of all neighbors of neigh
                # and update degrees)
                neighbors = g.adj_list[neigh]
                # TODO Dude this code is ugly.. why no loop?
                list(map(lambda n: g.adj_list[n].remove(neigh), neighbors))
                list(map(g.degrees_decrement, neighbors))
                neighbors.clear()
                g.update_to_zero_degree(neigh)
                g.vc_solution[neigh] = True
                vc_size += 1

            # We changed the graph (if there were degree 1 vertices), so we update
            # the max degree
            g.update_max_deg()
            if g.max_deg == 0:  # return solution size
                return vc_size

            self.recursive_steps += 1

            # If rnd_vertex is set, pick a random vertex v with maximal degree
            # otherwise pick the vertex v at the end of the max_deg degree bag.
            # The flag rnd_vertex allows for variation in the solution and may
            # (sometimes) lead to better results.
            if rnd_vertex:
                # success, vc_reduction, revert_reduction_array, k_update = reduction.perform_reduction(g, float("inf"), max_deg, reduction_grouping=self.reduction_grouping, preprocessing=True)
                # vc_size += k - k_update
                v = self.get_random_vertex(g)
                g.deg_bags[g.max_deg].remove(v)
            else:
                v = g.deg_bags[g.max_deg].pop()

            # "Remove" v from the graph and add v to the vertex cover
            g.degrees[v] = 0
            g.deg_bags[0].add(v)
            neighbors_v = g.adj_list[v]
            list(map(lambda n: g.adj_list[n].remove(v), neighbors_v))
            list(map(g.degrees_decrement, neighbors_v))
            g.adj_list[v].clear()
            g.vc_solution[v] = True
            vc_size += 1

        return vc_size

    def constrained_branching(self, g: VCGraph, vc_size: int, ub: int) -> int:
        k = ub - vc_size - 1  # be better
        if k < 0:
            g.CliqueCover.reset_records()
            return ub

        (success, vc_reduction, num_revert_red, k_update) = reduction.perform_reduction(
            g, k, reduction_grouping=self.reduction_grouping, preprocessing=True
        )
        vc_size += k - k_update
        if not success or g.Packing.packing_is_violated():
            g.CliqueCover.reset_records()
            g.execute_from_revert_stack(num_revert_red)
            g.recently_updated_vertices.clear()
            return ub

        lb = g.CliqueCover.update()

        if lb + vc_size >= ub:
            # print("####",lb, vc_size, ub)
            g.CliqueCover.revert_update()
            g.execute_from_revert_stack(num_revert_red)
            g.recently_updated_vertices.clear()
            return ub

        g.update_max_deg()
        max_deg = g.max_deg
        if g.max_deg == 0:
            if vc_size < ub:
                self.best_solution = g.vc_solution.copy()
                self.best_solution_merge_stack = g.merge_stack.copy()
                ub = vc_size
                # print("#new ub: ",ub)
            # clean_up
            g.CliqueCover.revert_update()
            g.execute_from_revert_stack(num_revert_red)
            g.recently_updated_vertices.clear()
            return ub

        self.recursive_steps += 1
        # CASE MAX_DEG FIRST
        # LEFT BRANCH ------------------------------------------
        v = g.deg_bags[g.max_deg].pop()
        g.deg_bags[g.max_deg].add(v)

        neighbors_v = list(g.adj_list[v])
        g.revert_stack.append(g.remove_edges(v, g.adj_list[v], packing_init=True))
        ub = self.constrained_branching(g, vc_size + 1, ub)
        if lb != g.CliqueCover.get_c_lb():
            print(lb, g.CliqueCover.get_c_lb())
        # clean up for v
        g.execute_from_revert_stack(1)
        g.recently_updated_vertices.clear()
        g.set_max_deg(max_deg)
        # RIGHT BRANCH --------------------------------
        # take all neighbors
        num_revert = g.remove_edges_for_array(neighbors_v, packing_init=True, v=v)

        ub = self.constrained_branching(g, vc_size + len(neighbors_v), ub)

        # if not undo changes:
        g.CliqueCover.revert_update()
        g.execute_from_revert_stack(num_revert + num_revert_red)
        g.recently_updated_vertices.clear()
        return ub

    def print_all_lowerbounds(self, g: VCGraph):
        g.set_click_cover()
        lb = g.CliqueCover.clique_lb(1, g.deg_bags)
        lb_lp = get_lowerbound_lp(g.adj_list, len(g.adj_list))

        self.last_k = max(lb, lb_lp)
        print("# clique cover update:")
        print("#lb1: %s" % lb)
        print("# lp via maximum bipartite matching:")
        print("#lb3: %s" % lb_lp)

    def get_lower_bound(self, g: VCGraph) -> int:
        lb_lp = get_lowerbound_lp(g.adj_list, len(g.adj_list))
        g.set_click_cover()
        lb = g.CliqueCover.clique_lb(1, g.deg_bags)
        return max(lb_lp, lb)

    def resolve_merged_degree_2(self, merged_vertices, solution: np.ndarray):
        (v, x, y) = merged_vertices
        if solution[v]:
            if np.any(solution[[x, y]]):
                raise Exception("merging failed for v == True!!!")
            solution[[x, y]] = True
            solution[v] = False
        else:
            solution[v] = True
            if np.any(solution[[x, y]]):
                raise Exception("merging failed!!!")

    def resolve_merged_degree_3(self, merged_vertices, solution: np.ndarray):
        (v, a, b, c) = merged_vertices
        num_of_abc_in_vc = int(solution[a] & solution[b] & solution[c])
        if num_of_abc_in_vc == 2:
            solution[v] = True
            if solution[a] & solution[b] == True:
                solution[a] = False
            if solution[b] & solution[c] == True:
                solution[b] = False
            if solution[c] & solution[a] == True:
                solution[c] = False
        if num_of_abc_in_vc == 1:
            solution[v] = True
            if solution[a] == True:
                solution[a] = False
            if solution[a] == True:
                solution[a] = False
            if solution[a] == True:
                solution[a] = False

    def resolve_merged_vertices(self, merge_stack: list, solution: np.ndarray):
        while len(merge_stack) > 0:
            merged_v = merge_stack.pop()
            if len(merged_v) == 3:
                self.resolve_merged_degree_2(merged_v, solution)
            if len(merged_v) == 4:
                self.resolve_merged_degree_3(merged_v, solution)

    def receive_signal_last_k(self, signum, frame):
        print("#last-k: %s" % (self.last_k))
        print("#recursive steps: %s" % self.recursive_steps)
        print("#", self.reduction_hit_counter)
        sys.exit()
