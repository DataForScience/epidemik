import unittest
from  epidemik import EpiModel, NetworkEpiModel
import networkx as nx

class NetworkEpiModelTestCase(unittest.TestCase):
    def setUp(self):
        self.N = 300
        self.G_full = nx.erdos_renyi_graph(self.N, p=1.)
        self.beta = 0.05
        self.SI_full = NetworkEpiModel(self.G_full)
    
    def test_named_parameters(self):
        self.SI_full.add_interaction("S", "I", "I", beta=self.beta)
        self.assertIn("beta",
            self.SI_full.params,
            "The parameter beta should be in the params dictionary")

    def test_unnamed_parameters(self):
        self.SI_full.add_interaction("S", "I", "I", self.beta)
        self.assertIn("rate1",
            self.SI_full.params,
            "The parameter rate1 should be in the params dictionary")

