import numpy as np
from lower_bound import get_lowerbound_lp
from vc_graph import COVERED, DTYPE, VCGraph, get_cover_vertices, resolve_merged_vertices
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
        reduction_mode=2,
        reduction_frequency=0,
        print_lower_bound=False,
        **kwargs,
    ):
        """This is the class container for the Vertex Cover Solver
        method: run"""
        self.recursive_steps = recursive_steps  # to get this information in the results
        self.last_k = 0
        self.greedy = greedy
        (self.reduction_mode, self.reduction_frequency, self.print_lower_bound, self.preprocessing, self.local_search,) = (
            reduction_mode,
            reduction_frequency,
            print_lower_bound,
            preprocessing,
            local_search,
        )
        self.branching_count = 0
        self.random_candidates = list()
        self.last_max_deg = -1
        self.time_limit = time_limit if time_limit else time.time() + 50
        # greedy init
        self.upper_bound = num_vertices
        self.pre_solution = np.zeros(num_vertices, dtype=DTYPE)
        self.best_solution = np.zeros(num_vertices, dtype=DTYPE)
        self.best_solution_merge_stack = None
        self.print_kernel = kwargs.get("print_kernel", False)
        print(
            "# init: reduction_mode: %d, reduction_frequency: %d, print_lower_bound: %d, preprocessing: %d, greedy: %d, local_search: %d"
            % (
                reduction_mode,
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
            print("# preprocessing")
            vc_solved = self.run_preprocessing(g)
            if vc_solved:
                resolve_merged_vertices(g.merge_stack, g.vc_solution)
                return get_cover_vertices(g.vc_solution)
        if self.greedy:
            print("# greedy")
            self.compute_upper_bound(g)
            print("# upper_bound: ", self.upper_bound)

        g.vc_solution = self.pre_solution.copy()
        self.best_solution_merge_stack = g.merge_stack
        g.set_packing_constraint()
        g.set_click_cover()
        g.set_max_deg(g.degrees.max())
        print("# run constrained branching")
        _ = self.constrained_branching(g, 0, self.upper_bound)
        g.merge_stack = self.best_solution_merge_stack
        # return solution
        resolve_merged_vertices(g.merge_stack, self.best_solution)
        return get_cover_vertices(self.best_solution)

    def run_preprocessing(self, g: VCGraph) -> bool:
        flag, vc_reduction, _, _ = reduction.perform_reduction(
            g,
            float("inf"),
            rec_steps=0,
            reduction_mode=2,
        )
        self.pre_solution[:] = g.vc_solution.copy()

        if len(g.deg_bags[0]) == len(g.adj_list):  # all vertices v deg(v) == 0 [vc] is solved
            return True

        if self.print_kernel:
            with open("kernel.dimacs", "w") as f:
                f.write(f"#{len(g.adj_list)} {g.num_edges}\n")
                edges = [(u, v) for u in range(len(g.adj_list)) for v in g.adj_list[u] if u < v]
                for u, v in edges:
                    f.write(f"{u+1} {v+1}\n")
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

    def compute_upper_bound(self, g: VCGraph):
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

        self.upper_bound = get_cover_vertices(self.best_solution).size

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
                list(map(lambda n: g.adj_list[n].remove(neigh), neighbors))
                list(map(g.degrees_decrement, neighbors))
                neighbors.clear()
                g.update_to_zero_degree(neigh)
                g.update_solution(neigh, COVERED)
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
            g.update_solution(v, COVERED)
            vc_size += 1

        return vc_size

    def constrained_branching(self, g: VCGraph, vc_size: int, ub: int) -> int:
        k = ub - vc_size - 1  # be better
        if k < 0:
            g.CliqueCover.reset_records()
            return ub

        (success, vc_reduction, num_revert_red, k_update) = reduction.perform_reduction(
            g, k, rec_steps=self.recursive_steps, reduction_mode=self.reduction_mode
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
