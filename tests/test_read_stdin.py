import unittest
import sys
import io
import time
from vcsolver.read_stdin import read_graph_fast

class TestReadGraph(unittest.TestCase):

    def mock_stdin(self, mock_input):
        sys.stdin = io.StringIO(mock_input)

    def test_read_graph_fast(self):
        # Mock the standard input
        # read input for file 
        open_file = open("tests/data/kernel.dimacs", "r")
        self.mock_stdin(open_file.read())
        open_file.close()
        start_time = time.time()
        adjacency_list, edges = read_graph_fast()
        end_time = time.time()
        print("Time taken: ", end_time - start_time)
        self.assertEqual(len(adjacency_list), 200)
        self.assertEqual(len(edges), 837)
        self.assertGreaterEqual(1.0, end_time - start_time)

if __name__ == "__main__":
    unittest.main()
