import numpy as np
import sys, signal, resource
from graph_tool.topology import label_components
from graph_tool import GraphView
from lower_bound import get_lowerbound_max_match, get_lowerbound_lp
from vc_graph import VCGraph, generate_deg_bags
from vectorized_func import execute_all
from read_stdin import read_graph_fast
from collections import Counter
import reduction
from sat_solver import min_sat_solver
import time
import random


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


def update_max_deg(g: VCGraph, max_deg: int) -> int:
    while (
        len(g.deg_bags[max_deg]) == 0 and max_deg > 0
    ):  # while bag of vertices for given degree is empty, lower max_deg by 1
        max_deg -= 1
    return max_deg


class VC:
    def __init__(
        self,
        num_vertices: int,
        time_limit=None,
        recursive_steps=0,
        red_grp=2,
        red_freq=0,
        print_lb=0,
        pre=0,
        branching1=0,
        branching2=0,
        sat=0,
        greedy=0,
        local_search=0,
        **kwargs,
    ):
        """This is the class container for the Vertex Cover Solver
        method: vc_branch"""
        self.recursive_steps = recursive_steps  # to get this information in the results
        self.last_k = 0
        self.branching1, self.branching2, self.sat, self.greedy = (
            branching1,
            branching2,
            sat,
            greedy,
        )
        self.red_grp, self.red_freq, self.print_lb, self.pre, self.local_search = (
            red_grp,
            red_freq,
            print_lb,
            pre,
            local_search,
        )
        self.counter_num_hits_by_red_rule = Counter()
        self.branching_count = 0
        self.random_candidates = list()
        self.last_max_deg = -1
        self.crt_g = None
        self.time_limit = time_limit if time_limit else time.time() + 50
        # greedy init
        self.upperbound = num_vertices
        self.pre_solution = np.zeros(num_vertices, dtype=bool)
        self.best_solution = np.zeros(num_vertices, dtype=bool)
        self.best_solution_merge_stack = None
        if sum([self.branching1, self.branching2, self.sat]) != 1:
            raise Exception("select only one of branching1, branching2, sat")
        print(
            "# init: red_grp: %d, red_freq: %d, print_lb: %d, pre: %d, branching1: %d, branching2: %d, sat: %d, greedy: %d, local_search: %d"
            % (
                red_grp,
                red_freq,
                print_lb,
                pre,
                branching1,
                branching2,
                sat,
                greedy,
                local_search,
            )
        )

    def run(self, g: VCGraph):
        if self.print_lb:
            self.print_all_lowerbounds(g)
        if self.pre:
            vc_solved = self.preprocessing(g)
            self.counter_num_hits_by_red_rule.update(g.num_hits_by_reduction_rule)
            if vc_solved:
                self.resolve_merged_vertices(g.merge_stack, g.vc_solution)
                return np.nonzero(g.vc_solution)[0]
        if self.sat:
            self.best_solution[min_sat_solver(g.edges)] = True
        if self.greedy:
            self.compute_upperbound(g)
            print("# upperbound: ", self.upperbound)
        if self.branching1:
            self.split_graph_into_components_and_run(g)
            self.best_solution[:] = g.vc_solution
        if self.branching2:
            g.vc_solution = self.pre_solution.copy()
            self.best_solution_merge_stack = g.merge_stack
            g.set_packing_constraint()
            g.set_click_cover()
            ub = self.constrained_branching(g, g.degrees.max(), 0, self.upperbound)
            g.merge_stack = self.best_solution_merge_stack
            self.counter_num_hits_by_red_rule.update(g.num_hits_by_reduction_rule)
        # return solution
        self.resolve_merged_vertices(g.merge_stack, self.best_solution)
        return np.nonzero(self.best_solution)[0]

    def preprocessing(self, g: VCGraph) -> bool:
        all_red_rules_activated = self.red_grp  # NOTE this was equal 4 in the solver3
        flag, vc_reduction, _, _, _ = reduction.perform_reduction(
            g,
            float("inf"),
            g.degrees.max(),
            all_red_rules_activated,
            self.red_freq,
            preprocessing=True,
        )
        self.pre_solution[:] = g.vc_solution.copy()

        if len(g.deg_bags[0]) == len(g.adj_list):  # all vertices v deg(v) == 0 [vc] is solved
            return True

        edges = [
            np.array([(e1, e2) for e2 in g.adj_list[e1] if e1 < e2])
            for e1 in np.arange(len(g.adj_list))
            if len(g.adj_list[e1]) > 0
        ]
        edges = [es for es in edges if len(es) > 0]
        g.edges = np.concatenate(edges)
        return False

    def get_random_vertex(self, g: VCGraph, max_deg: int):
        if self.last_max_deg != max_deg:
            self.last_max_deg = max_deg
            self.random_candidates = list(g.deg_bags[max_deg])
            random.shuffle(self.random_candidates)
        v_deg, v = 0, -1
        while max_deg != v_deg:
            v = self.random_candidates.pop()
            v_deg = g.degrees[v]
        return v

    def compute_upperbound(self, g: VCGraph):
        num_runs = min(100, max(int(10**6 / len(g.edges)), 1))
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

    def greedy_heuristic(self, g: VCGraph, max_deg: int, vc_size: int = 0, rnd_vertex=False) -> np.ndarray:

        while max_deg != 0:

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
            max_deg = update_max_deg(g, max_deg)
            if max_deg == 0:  # return solution size
                return vc_size

            self.recursive_steps += 1

            # If rnd_vertex is set, pick a random vertex v with maximal degree
            # otherwise pick the vertex v at the end of the max_deg degree bag.
            # The flag rnd_vertex allows for variation in the solution and may
            # (sometimes) lead to better results.
            if rnd_vertex:
                # success, vc_reduction, revert_reduction_array, k_update, max_deg = reduction.perform_reduction(g, float("inf"), max_deg, red_grp=self.red_grp, preprocessing=True)
                # vc_size += k - k_update

                v = self.get_random_vertex(g, max_deg)
                g.deg_bags[max_deg].remove(v)
            else:
                v = g.deg_bags[max_deg].pop()

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

    def constrained_branching(self, g: VCGraph, max_deg: int, vc_size: int, ub: int) -> int:
        k = ub - vc_size - 1  # be better
        if k < 0:
            g.CliqueCover.reset_records()
            return ub

        (
            success,
            vc_reduction,
            num_revert_red,
            k_update,
            max_deg,
        ) = reduction.perform_reduction(g, k, max_deg, red_grp=self.red_grp, preprocessing=True)
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

        max_deg = update_max_deg(g, max_deg)
        if max_deg == 0:
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
        v = g.deg_bags[max_deg].pop()
        g.deg_bags[max_deg].add(v)

        neighbors_v = list(g.adj_list[v])
        g.revert_stack.append(g.remove_edges(v, g.adj_list[v], packing_init=True))
        ub = self.constrained_branching(g, max_deg, vc_size + 1, ub)
        if lb != g.CliqueCover.get_c_lb():
            print(lb, g.CliqueCover.get_c_lb())
        # clean up for v
        g.execute_from_revert_stack(1)
        g.recently_updated_vertices.clear()

        # RIGHT BRANCH --------------------------------
        # take all neighbors
        num_revert = g.remove_edges_for_array(neighbors_v, packing_init=True, v=v)

        ub = self.constrained_branching(g, max_deg, vc_size + len(neighbors_v), ub)

        # if not undo changes:
        g.CliqueCover.revert_update()
        g.execute_from_revert_stack(num_revert + num_revert_red)
        g.recently_updated_vertices.clear()
        return ub

    def vc_branch(self, g: VCGraph, k: int, max_deg: int) -> bool:
        """
        branching over vertices, those with highest degree first
        param: int k
               int max_deg
        return: bool vc_solved indicates whether the VC is solved for k
        """
        if k < 0:
            g.CliqueCover.reset_records()
            return False

        max_deg = update_max_deg(g, max_deg)
        if max_deg == 0:
            return True  # point of recursion end

        success, vc_reduction, num_revert_red, k, max_deg = reduction.perform_reduction(
            g, k, max_deg, self.red_grp, self.red_freq, rec_steps=self.recursive_steps
        )
        lb = g.CliqueCover.update()
        if not success or lb > k:
            # if lb>k: print("###", lb, k, max_deg)
            g.CliqueCover.revert_update()
            g.execute_from_revert_stack(num_revert_red)
            g.recently_updated_vertices.clear()
            return False

        max_deg = update_max_deg(g, max_deg)
        if max_deg == 0:
            return True  # point of recursion end

        self.recursive_steps += 1
        # CASE MAX_DEG FIRST
        # LEFT BRANCH ------------------------------------------

        v = g.deg_bags[max_deg].pop()
        g.deg_bags[max_deg].add(v)

        neighbors_v = list(g.adj_list[v])
        g.revert_stack.append(g.remove_edges(v, g.adj_list[v]))
        vc_solved = self.vc_branch(g, k - 1, max_deg)
        if vc_solved:
            return vc_solved

        # clean up for v
        g.execute_from_revert_stack(1)
        g.recently_updated_vertices.clear()

        # RIGHT BRANCH --------------------------------
        # take all neighbors
        num_revert = g.remove_edges_for_array(neighbors_v)

        vc_solved = self.vc_branch(g, k - len(neighbors_v), max_deg)
        if vc_solved:
            return vc_solved

        # if not undo changes:
        g.CliqueCover.revert_update()
        g.execute_from_revert_stack(num_revert + num_revert_red)
        g.recently_updated_vertices.clear()
        return vc_solved

    def print_all_lowerbounds(self, g: VCGraph):
        g.set_click_cover()
        lb = g.CliqueCover.clique_lb(1, g.deg_bags)
        lb_max_c_m = get_lowerbound_max_match(g.gt_G)
        lb_lp = get_lowerbound_lp(g.gt_G.get_edges(), len(g.adj_list))

        self.last_k = max(lb, lb_lp, lb_max_c_m)
        print("# clique cover update:")
        print("#lb1: %s" % lb)
        print("# max_match:")
        print("#lb2: %s" % lb_max_c_m)
        print("# lp via maximum bipartite matching:")
        print("#lb3: %s" % lb_lp)

    def get_lower_bound(self, g: VCGraph) -> int:
        lb_lp = get_lowerbound_lp(g.gt_G.get_edges(), len(g.adj_list))
        lb_max_c_m = get_lowerbound_max_match(g.gt_G)
        g.set_click_cover()
        lb = g.CliqueCover.clique_lb(1, g.deg_bags)
        return max(lb_max_c_m, lb_lp, lb)

    def get_vertex_cover(self, g: VCGraph) -> np.ndarray:
        k = self.get_lower_bound(g)
        print("# lb:", k)
        self.last_k = k + g.vc_solution.sum()  # only for printing purpose
        max_deg = update_max_deg(g, len(g.deg_bags) - 1)
        while True:
            try:
                vc_solved = self.vc_branch(g, k, max_deg)

            except RecursionError as re:
                print("#Recursion Error max depth exceded")
                print("#last-k: %s" % (self.last_k))
                sys.exit()
            if vc_solved:
                return np.nonzero(g.vc_solution)[0]  # get indices that are True and therefore part of the solution
            k += 1
            self.last_k += 1
        return None

    def split_graph_into_components_and_run(self, g: VCGraph):
        g.gt_G.add_edge_list(g.edges)
        comp, hist = label_components(g.gt_G)
        comp_index = np.where(hist > 1)[0]
        sortminfirst = hist[comp_index].argsort()
        graphs = [GraphView(g.gt_G, vfilt=comp.a == index, directed=False) for index in comp_index[sortminfirst]]
        vc = []
        self.counter_num_hits_by_red_rule.update(g.num_hits_by_reduction_rule)
        # print("#",hist[comp_index[sortminfirst]])
        for gt_g in graphs:
            vertices = gt_g.get_vertices()
            g_deg_bags = generate_deg_bags(g.degrees[vertices], indices=vertices)
            self.crt_g = VCGraph(
                g.adj_list,
                gt_g.get_edges(),
                len(vertices),
                vc_solution=g.vc_solution,
                degrees=g.degrees,
                deg_bags=g_deg_bags,
                gt_graph=gt_g,
            )
            vc = self.get_vertex_cover(self.crt_g)
            g.merge_stack.extend(self.crt_g.merge_stack)
            self.counter_num_hits_by_red_rule.update(self.crt_g.num_hits_by_reduction_rule)

        self.resolve_merged_vertices(
            g.merge_stack, g.vc_solution
        )  # resolves merges for all graphs and preprocessing at the very end

    def resolve_merged_vertices(self, merge_stack: list, solution: np.ndarray):
        while len(merge_stack) > 0:
            temp = merge_stack.pop()

            # resolve merging deg2
            if len(temp) == 3:
                (v, x, y) = temp
                if solution[v]:
                    if np.any(solution[[x, y]]):
                        raise Exception("mergeing failed for v == True!!!")
                    solution[[x, y]] = True
                    solution[v] = False
                else:
                    solution[v] = True
                    if np.any(solution[[x, y]]):
                        raise Exception("mergeing failed!!!")
            # resolve merging deg3
            if len(temp) == 4:
                (v, a, b, c) = temp
                num_of_abc_in_vc = int(solution[a] & solution[b] & solution[c])
                # if num_of_abc_in_vc == 3: # do S' = S
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

    def receive_signal_last_k(self, signum, frame):
        print("#last-k: %s" % (self.last_k))
        print("#recursive steps: %s" % self.recursive_steps)
        if self.crt_g != None:
            self.counter_num_hits_by_red_rule.update(self.crt_g.num_hits_by_reduction_rule)
        print("#", self.counter_num_hits_by_red_rule)
        sys.exit()


def main(args):
    start_time = time.time()
    kwargs = {
        "pre": 1,
        "branching1": 0,
        "branching2": 0,
        "sat": 0,
        "greedy": 0,
        "local_search": 0,
        "red_grp": 2,
        "red_freq": 0,
        "print_lb": 0,
        "ignore_limit": 0,
    }

    for arg in args:
        k = arg.split("=")[0]
        v = arg.split("=")[1]
        if k in kwargs:
            kwargs[k] = int(v)
        else:
            raise LookupError("keyword %s is not valid!" % k)

    newslimit = 65536000
    slimit, hlimit = resource.getrlimit(resource.RLIMIT_STACK)
    if kwargs["ignore_limit"] and slimit < newslimit:
        raise Exception("stack limit is too small, run ulimit -Ss 64000")
    sys.setrecursionlimit(10**6)  # careful segfault

    adjacency_list, edges = read_graph_fast()  # read from stdin
    num_vertices = len(adjacency_list)

    if len(edges) == 0:  # if graph has no edges
        print("#recursive steps: 0")
        return None  # ugly workaround, to avoid this unnecessary in vs_branch,

    print("# finished reading input after %s sec" % (time.time() - start_time))
    vertexC = VC(num_vertices, time_limit=50 + start_time, **kwargs)
    graph = VCGraph(adjacency_list, edges, num_vertices)
    signal.signal(signal.SIGINT, vertexC.receive_signal_last_k)
    print("# finished init after %s sec" % (time.time() - start_time))
    vc = vertexC.run(graph)
    print("# solution after %s sec" % (time.time() - start_time))

    print("#recursive steps: %s" % vertexC.recursive_steps)
    print("# solsize: %s" % vc.size)
    print("#", vertexC.counter_num_hits_by_red_rule)
    # for v in vc:
    #    print(v+1) # nameing back +1 {1,...,n} = Vertices
    sys.stdout.write("\n".join(map(str, vc + 1)))
    print("\n# output after %s sec" % (time.time() - start_time))


if __name__ == "__main__":
    main(sys.argv[1:])
