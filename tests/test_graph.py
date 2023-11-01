import unittest
import sys
import io
import time
from vcsolver.read_stdin import read_graph_fast
from vcsolver.vc_graph import VCGraph
from vcsolver.vc_solver import VCSolver

class TestGraph(unittest.TestCase):

    def setUp(self):
        instance = "../data/4-interst/kernel"
        file = open(instance+".dimacs", "r")
        self.mock_stdin(file.read())
        file.close()
        file_sol = open(instance+".solution", "r")
        self.solution = int(file_sol.readline().strip())
        file_sol.close()
        adjacency, edges = read_graph_fast()
        self.start_time = time.time()
        self.graph = VCGraph(adjacency)
        self.solver = VCSolver(len(adjacency), time_limit=50 + self.start_time, print_lower_bound=True)

    def mock_stdin(self, mock_input):
        sys.stdin = io.StringIO(mock_input)

    def test_get_current_vertices(self):
        start_time = time.time()
        self.assertEqual(len(self.graph.get_current_vertices()), 192)
        end_time = time.time()
        print("Time taken: ", end_time - start_time)
        self.assertGreaterEqual(1e-4, end_time - start_time)

    def test_solver(self):
        start_time = time.time()
        vc = self.solver.run(self.graph)
        end_time = time.time()
        self.assertEqual(vc.size, self.solution)
        print("Time taken: ", end_time - start_time)
        self.assertGreaterEqual(20, end_time - start_time)
        print("# ",self.solver.reduction_hit_counter)

if __name__ == "__main__":
    unittest.main()