import unittest
from epidemik import EpiModel
from sklearn.linear_model import LinearRegression
import numpy as np

class EpiModelTestCase(unittest.TestCase):
    def setUp(self):

        self.beta = 0.3
        self.mu = 0.1
        self.birth = 0.3
        self.death = 0.3
        self.fixed_birth = 10

        self.SIR = EpiModel()
        self.SIR.add_interaction("S", "I", "I", beta=self.beta)
        self.SIR.add_spontaneous("I", "R", mu=self.mu)

        self.birth_test = EpiModel()
        self.birth_test.add_interaction("S", "I", "I", beta=self.beta)
        self.birth_test.add_birth_rate("S", b=self.birth)

        self.fixed_birth_test = EpiModel()
        self.fixed_birth_test.add_interaction("S", "I", "I", beta=self.beta)
        self.fixed_birth_test.add_birth_rate("S", b=self.fixed_birth, fixed=True, global_rate=False)

        self.death_test = EpiModel()
        self.death_test.add_interaction("S", "I", "I", beta=self.beta)
        self.death_test.add_death_rate(d=self.death)

    def test_R0(self):
        self.assertEqual(self.SIR.R0(), 3.0, "incorrect R0")

    def test_edges(self):
        edges = list(self.SIR.transitions.edges(data=True))
        self.assertEqual(len(edges), 2)

        for edge in edges:
            data = edge[2]

            if "agent" in data:
                self.assertEqual(data["agent"], "I")
                self.assertEqual(data["rate"], "beta")
                self.assertEqual(self.SIR.params["beta"], self.beta)
                self.assertEqual(edge[0], "S")
                self.assertEqual(edge[1], "I")
            else:
                self.assertEqual(data["rate"], "mu")
                self.assertEqual(self.SIR.params["mu"], self.mu)
                self.assertEqual(edge[0], "I")
                self.assertEqual(edge[1], "R")

    def test_birth(self):
        self.assertEqual(self.birth_test.transitions.nodes['S']["birth"], "b")
        self.assertIn("b", self.birth_test.params.keys())
        self.assertEqual(self.birth_test.params["b"], self.birth)

    def test_birth_rate(self):
        self.birth_test.integrate(10, S=990, I=10)
        values = self.birth_test.values_
        values['total'] = values.sum(axis=1)
        values = values.reset_index()

        lm = LinearRegression()
        lm.fit(values['index'].values.reshape(-1, 1), np.log(values['total']))

        self.assertAlmostEqual(lm.coef_[0], self.birth, delta=0.01)

    def test_fixed_birth_rate(self):
        self.fixed_birth_test.integrate(10, S=990, I=10)
        values = self.fixed_birth_test.values_
        values['total'] = values.sum(axis=1)
        values = values.reset_index()

        lm = LinearRegression()
        lm.fit(values['index'].values.reshape(-1, 1), values['total'])

        self.assertAlmostEqual(lm.coef_[0], 10, delta=0.01)

    def test_death(self):
        self.assertEqual(self.death_test.transitions.nodes['S']["death"], "d")
        self.assertEqual(self.death_test.transitions.nodes['I']["death"], "d")

        self.assertIn("d", self.death_test.params.keys())
        self.assertEqual(self.death_test.params["d"], self.death)

    def test_death_rate(self):
        self.death_test.integrate(10, S=990, I=10)
        values = self.death_test.values_
        values['total'] = values.sum(axis=1)
        values = values.reset_index()

        lm = LinearRegression()
        lm.fit(values['index'].values.reshape(-1, 1), np.log(values['total']))

        self.assertAlmostEqual(lm.coef_[0], -self.death, delta=0.01)