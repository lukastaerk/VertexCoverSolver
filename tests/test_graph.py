import unittest
import sys
import io
import time
from vcsolver.read_stdin import read_graph_fast
from vcsolver.vc_graph import VCGraph

class TestGraph(unittest.TestCase):

    def setUp(self):
        file = open("../data/2-social-networks/09-email.graph.dimacs", "r")
        self.mock_stdin(file.read())
        file.close()
        adjacency, edges = read_graph_fast()
        self.graph = VCGraph(adjacency)

    def mock_stdin(self, mock_input):
        sys.stdin = io.StringIO(mock_input)

    def test_get_current_vertices(self):
        start_time = time.time()
        self.assertEqual(len(self.graph.get_current_vertices()), len(self.graph.adj_list))
        end_time = time.time()
        print("Time taken: ", end_time - start_time)
        self.assertGreaterEqual(1e-4, end_time - start_time)

if __name__ == "__main__":
    unittest.main()