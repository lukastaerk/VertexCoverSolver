import unittest
import sys
import io
import time
from vcsolver.read_stdin import read_graph_fast
from vcsolver.vc_graph import VCGraph
from vcsolver.vc_solver import VCSolver
from vcsolver.verifier import verifier

class TestGraph(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super(TestGraph, self).__init__(*args, **kwargs)
        self.start_time = time.time()
        instance = "../data/4-interst/kernel"
        self.adjacency, self.edges, self.solution = self.read_graph(instance)
        self.graph = VCGraph(self.adjacency)
        self.solver = VCSolver(len(self.adjacency), time_limit=50 + self.start_time, print_lower_bound=True, reduction_mode=1)

    def read_graph(self, instance):
        file = open(instance+".dimacs", "r")
        self.mock_stdin(file.read())
        file.close()
        file_sol = open(instance+".solution", "r")
        solution = int(file_sol.readline().strip())
        file_sol.close()
        adj, edges = read_graph_fast()
        return adj, edges, solution
        

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
        self.assertTrue(verifier(vc, self.edges))

    def test_deg3(self):
        instance = "../data/3-medium-sized/vc019"
        adjacency, edges, solution = self.read_graph(instance)
        graph = VCGraph(adjacency)
        solver = VCSolver(len(adjacency), time_limit=50 + self.start_time, reduction_mode=1)
        vc = solver.run(graph)
        self.assertTrue(verifier(vc, edges))
        self.assertEqual(vc.size, solution)
        



if __name__ == "__main__":
    unittest.main()