import unittest
from vcsolver.reduction import is_unconfined
from vcsolver.vc_graph import VCGraph

class TestUnconfinedRule(unittest.TestCase):

    def test_case_1(self):
        adj_list = [{1, 2, 3}, {0, 2}, {0, 1}, {0}]
        g = VCGraph(adj_list)

        self.assertTrue(is_unconfined(g, 0))
        self.assertTrue(is_unconfined(g, 1))
        self.assertTrue(is_unconfined(g, 2))
        self.assertFalse(is_unconfined(g, 3))

    def test_case_2(self):
        adj_list = [{1, 2}, {0}, {0}]
        g = VCGraph(adj_list)

        self.assertTrue(is_unconfined(g, 0))
        self.assertFalse(is_unconfined(g, 1))
        self.assertFalse(is_unconfined(g, 2))

    def test_case_3(self):
        adj_list = [{1}, {0, 2, 3, 4}, {1}, {1}, {1}]
        g = VCGraph(adj_list)

        self.assertFalse(is_unconfined(g, 0))
        self.assertTrue(is_unconfined(g, 1))
        self.assertFalse(is_unconfined(g, 2))
        self.assertFalse(is_unconfined(g, 3))
        self.assertFalse(is_unconfined(g, 4))

    def test_case_4(self):
        adj_list = [{1, 2, 4}, {0, 2, 3}, {0, 1}, {1, 4}, {0, 3}]
        g = VCGraph(adj_list)

        self.assertTrue(is_unconfined(g, 0))
        self.assertTrue(is_unconfined(g, 1))
        self.assertTrue(is_unconfined(g, 2))
        self.assertTrue(is_unconfined(g, 3))
        self.assertTrue(is_unconfined(g, 4))

if __name__ == '__main__':
    unittest.main()
