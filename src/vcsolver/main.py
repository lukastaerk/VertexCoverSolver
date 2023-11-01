import argparse
import resource
import time
import signal
import sys

from vc_graph import VCGraph
from vc_solver import VCSolver
from read_stdin import read_graph_fast


def parse_arguments():
    parser = argparse.ArgumentParser(description="Vertex Cover vc_solver")

    parser.add_argument(
        "--preprocessing", "-p", action="store_true", help="Enable or disable preprocessing (default: disabled)"
    )
    parser.add_argument(
        "--greedy", "-g", action="store_true", help="Enable or disable greedy algorithm (default: disabled)"
    )
    parser.add_argument(
        "--local_search", "-ls", action="store_true", help="Enable or disable local search (default: disabled)"
    )
    parser.add_argument(
        "--print_lower_bound",
        "-lb",
        action="store_true",
        help="Enable or disable printing of lower bound (default: disabled)",
    )
    parser.add_argument("--ignore_limit", "-il", action="store_true", help="Ignore limit flag (default: disabled)")

    return parser.parse_args()


def main():
    start_time = time.time()
    args = parse_arguments()
    args_dict = vars(args)
    newslimit = 65536000
    slimit, hlimit = resource.getrlimit(resource.RLIMIT_STACK)
    if not args_dict["ignore_limit"] and slimit < newslimit:
        raise Exception("Stack limit is too small, run ulimit -Ss 64000 or use --ignore_limit flag")
    sys.setrecursionlimit(10**6)  # careful segfault

    adjacency_list, edges = read_graph_fast()  # read from stdin
    num_vertices = len(adjacency_list)

    if len(edges) == 0:  # if graph has no edges
        print("#recursive steps: 0")
        # ugly workaround, to avoid this unnecessary in vs_branch,
        return None

    print("# finished reading input after %s sec" % (time.time() - start_time))
    graph = VCGraph(adjacency_list)
    vc_solver = VCSolver(num_vertices, time_limit=50 + start_time, **args_dict)

    signal.signal(signal.SIGINT, vc_solver.receive_signal_last_k)
    print("# finished init after %s sec" % (time.time() - start_time))
    vc = vc_solver.run(graph)
    print("# solution after %s sec" % (time.time() - start_time))

    print("#recursive steps: %s" % vc_solver.recursive_steps)
    print("# solsize: %s" % vc.size)
    print("#", vc_solver.reduction_hit_counter)
    # for v in vc:
    #    print(v+1) # nameing back +1 {1,...,n} = Vertices
    sys.stdout.write("\n".join(map(str, vc + 1)))
    print("\n# output after %s sec" % (time.time() - start_time))


if __name__ == "__main__":
    main()
