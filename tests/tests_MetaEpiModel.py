import unittest

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pandas as pd

from epidemik import MetaEpiModel
from epidemik.utils import NotInitialized


class MetaEpiModelTestCase(unittest.TestCase):
	def setUp(self):
		np.random.seed(42)
		self.travel = pd.DataFrame({'A': [0.99, 0.1], 'B': [0.01, 0.9]}, index=["A", "B"])
		self.population = pd.DataFrame({'Population': [100_000, 10_000]}, index=["A", "B"])

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

	def test_validation(self):
		bad_rows = pd.DataFrame({'A': [0.5, 0.1], 'B': [0.01, 0.9]}, index=["A", "B"])
		with self.assertRaises(ValueError):
			MetaEpiModel(bad_rows, self.population)

		not_square = pd.DataFrame({'A': [1.0, 1.0]}, index=["A", "B"])
		with self.assertRaises(ValueError):
			MetaEpiModel(not_square, self.population)

		with self.assertRaises(ValueError):
			MetaEpiModel(self.travel, self.population.loc[["A"]])

		with self.assertRaises(ValueError):
			MetaEpiModel(self.travel, self.population, population='Pop')

	def test_simulate_run(self):
		self.SIR.simulate(20, seed_state='A', I=10)

		total = self.population['Population'].sum()
		self.assertEqual(self.SIR.compartments_.sum().sum(), total)
		self.assertEqual(list(self.SIR.compartments_.index), ['A', 'B'])

		state = self.SIR.get_state('A')
		# initial conditions + 20 steps
		self.assertEqual(state.values_.shape, (21, 3))
		self.assertEqual(list(state.values_.index[:2]), [1, 2])
		self.assertEqual(state.values_.iloc[0].tolist(), [99_990, 10, 0])

		with self.assertRaises(ValueError):
			self.SIR.simulate(5, seed_state='Z', I=10)

	def test_simulate_deprecated_timestamp(self):
		with self.assertWarns(DeprecationWarning):
			self.SIR.simulate(timestamp=5, seed_state='A', I=10)
		self.assertEqual(self.SIR.get_state('A').values_.shape, (6, 3))

	def test_seasonality_is_forwarded(self):
		season = np.ones(365)
		self.SIR.simulate(5, seasonality=season, seed_state='A', I=10)
		self.assertIs(self.SIR.get_state('B').seasonality, season)

	def test_vaccination_and_demographics(self):
		self.SIR.add_vaccination('S', 'V', 0.01, start=2)
		self.SIR.add_birth_rate(0.001)
		self.SIR.add_death_rate(0.001)
		self.SIR.simulate(10, seed_state='A', I=10)
		self.assertGreater(self.SIR.compartments_['V'].sum(), 0)

	def test_groups(self):
		model = MetaEpiModel(self.travel, self.population)
		model.add_interaction("Sh", "Ih", "Iv", 0.3)
		model.add_interaction("Sv", "Iv", "Ih", 0.3)
		model.add_spontaneous("Ih", "Rh", 0.1)
		model.add_groups({"host": ["Sh", "Ih", "Rh"], "vector": ["Sv", "Iv"]})
		model.simulate(5, seed_state='A', susceptible='Sh', Ih=10, Sv=100)

		self.assertFalse(model.compartments_.isnull().any().any())
		self.assertEqual(model.compartments_.sum().sum(), self.population['Population'].sum())

	def test_plots(self):
		with self.assertRaises(NotInitialized):
			self.SIR.plot()

		self.SIR.simulate(10, seed_state='A', I=10)

		ax = self.SIR.plot(title='test')
		self.assertIsNotNone(ax)
		matplotlib.pyplot.close('all')

		ax = self.SIR.plot_peaks()
		self.assertIsNotNone(ax)
		matplotlib.pyplot.close('all')

		with self.assertRaises(ValueError):
			self.SIR.plot(compartment='Q')

		fig, ax = matplotlib.pyplot.subplots()
		self.SIR.draw_model(ax=ax, show=False)
		matplotlib.pyplot.close(fig)


if __name__ == '__main__':
	unittest.main()
