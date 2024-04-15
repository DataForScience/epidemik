import unittest
import pandas as pd
from epidemik import MetaEpiModel
from epidemik.utils import NotInitialized

class MetaEpiModelTestCase(unittest.TestCase):
	def setUp(self):
		self.travel = pd.DataFrame({'A': [0.99, 0.1], 'B':[0.01, 0.9]}, index=["A", "B"])
		self.population = pd.DataFrame({'Population':[100_000, 10_000]}, index=["A", "B"])

		self.SIR = MetaEpiModel(self.travel, self.population)
		self.beta = 0.3
		self.mu = 0.1
		self.SIR.add_interaction('S', 'I', 'I', self.beta)
		self.SIR.add_spontaneous('I', 'R', self.mu)

	def test_number_populations(self):
		self.assertEqual(self.SIR.travel_graph.shape[0], 2)

	def test_R0(self):
		self.assertEqual(self.SIR.R0(), 3.0, 'incorrect R0')

	def test_simulate(self):
		with self.assertRaises(NotInitialized) as _:
			self.SIR.simulate(10)

	def test_initialize_populations(self):
		self.SIR._initialize_populations('S')
		self.assertEqual(self.SIR.compartments_['S'].sum(), self.population['Population'].sum())

	def test_travel(self):
		self.SIR._initialize_populations('S')
		new_compartments = self.SIR._run_travel(self.SIR.compartments_, self.travel)

		self.assertEqual(new_compartments.sum().sum(), self.population['Population'].sum())

	def test_integrate(self):
		with self.assertRaises(NotImplementedError) as _:
			self.SIR.integrate()