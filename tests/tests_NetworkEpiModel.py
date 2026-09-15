import unittest

import matplotlib
matplotlib.use("Agg")

import networkx as nx
import numpy as np

from epidemik import NetworkEpiModel


class NetworkEpiModelTestCase(unittest.TestCase):
	def setUp(self):
		np.random.seed(42)
		self.G = nx.erdos_renyi_graph(n=200, p=0.05, seed=42)
		self.beta = 0.1
		self.mu = 0.1
		self.SIR = NetworkEpiModel(self.G)
		self.SIR.add_interaction('S', 'I', 'I', self.beta)
		self.SIR.add_spontaneous('I', 'R', self.mu)

	def test_kavg(self):
		self.assertAlmostEqual(self.SIR.kavg_, 2*self.G.number_of_edges()/self.G.number_of_nodes())

	def test_R0(self):
		self.assertAlmostEqual(self.SIR.R0(), np.round(self.beta/self.mu*self.SIR.kavg_, 2))

	def test_R0_without_recovered_compartment(self):
		SIS = NetworkEpiModel(self.G)
		SIS.add_interaction('S', 'I', 'I', self.beta)
		SIS.add_spontaneous('I', 'S', self.mu)
		self.assertIsNotNone(SIS.R0())

	def test_rescale(self):
		model = NetworkEpiModel(self.G)
		model.add_interaction('S', 'I', 'I', self.beta, rescale=True)
		model.add_spontaneous('I', 'R', self.mu)
		self.assertAlmostEqual(model.R0(), np.round(self.beta/self.mu, 2))

	def test_integrate(self):
		with self.assertRaises(NotImplementedError):
			self.SIR.integrate(10)

	def test_simulate(self):
		seeds = {0: 'I', 1: 'I', 2: 'I'}
		self.SIR.simulate(50, seeds=seeds)

		self.assertEqual(self.SIR.population_.shape, (50, 200))
		self.assertEqual(list(self.SIR.population_.index[:2]), [1, 2])
		self.assertEqual(list(self.SIR.values_.columns), ['S', 'I', 'R'])
		self.assertEqual(self.SIR.values_.iloc[0].tolist(), [197, 3, 0])
		self.assertTrue((self.SIR.values_.sum(axis=1) == 200).all())
		self.assertGreater(self.SIR.values_['R'].iloc[-1], 0)

	def test_values_contain_every_compartment(self):
		# Even compartments that are never reached must be present
		SEIR = NetworkEpiModel(self.G)
		SEIR.add_interaction('S', 'E', 'I', 0.0)
		SEIR.add_spontaneous('E', 'I', 0.3)
		SEIR.add_spontaneous('I', 'R', self.mu)
		SEIR.simulate(5, seeds={0: 'I'})
		self.assertEqual(list(SEIR.values_.columns), ['S', 'E', 'I', 'R'])
		self.assertEqual(SEIR.values_['E'].sum(), 0)

	def test_multiple_spontaneous_transitions(self):
		model = NetworkEpiModel(self.G)
		model.add_interaction('S', 'I', 'I', self.beta)
		model.add_spontaneous('I', 'R', self.mu)
		model.add_spontaneous('I', 'D', 0.05)

		self.assertEqual(model.spontaneous['I'], {'R': self.mu, 'D': 0.05})

		model.simulate(100, seeds={0: 'I', 1: 'I', 2: 'I', 3: 'I', 4: 'I'})
		self.assertGreater(model.values_['R'].iloc[-1], 0)
		self.assertGreater(model.values_['D'].iloc[-1], 0)

	def test_long_compartment_names(self):
		model = NetworkEpiModel(self.G)
		model.add_interaction('Susceptible', 'Infected', 'Infected', self.beta)
		model.add_spontaneous('Infected', 'Recovered', self.mu)
		model.simulate(50, seeds={0: 'Infected', 1: 'Infected'})

		self.assertEqual(list(model.values_.columns), ['Susceptible', 'Infected', 'Recovered'])
		self.assertEqual(model.population_.iloc[0, 0], 'Infected')
		self.assertGreater(model.values_['Recovered'].iloc[-1], 0)

	def test_arbitrary_node_labels(self):
		G = nx.complete_graph(['a', 'b', 'c', 'd'])
		model = NetworkEpiModel(G)
		model.add_interaction('S', 'I', 'I', 1.0)
		model.add_spontaneous('I', 'R', 0.0)
		model.simulate(3, seeds={'a': 'I'})

		self.assertEqual(list(model.population_.columns), ['a', 'b', 'c', 'd'])
		self.assertEqual(model.values_.iloc[0].tolist(), [3, 1, 0])
		self.assertEqual(model.values_.iloc[1].tolist(), [0, 4, 0])

	def test_invalid_seeds(self):
		with self.assertRaises(ValueError):
			self.SIR.simulate(5, seeds={10_000: 'I'})
		with self.assertRaises(ValueError):
			self.SIR.simulate(5, seeds={0: 'Q'})

	def test_explicit_susceptible(self):
		self.SIR.simulate(5, seeds={0: 'I'}, susceptible='S')
		self.assertEqual(self.SIR.values_['S'].iloc[0], 199)
		with self.assertRaises(ValueError):
			self.SIR.simulate(5, seeds={0: 'I'}, susceptible='Q')

	def test_plot(self):
		self.SIR.simulate(10, seeds={0: 'I'})
		fig, ax = matplotlib.pyplot.subplots()
		self.SIR.plot(show=False, ax=ax)
		self.assertEqual(len(ax.lines), 3)
		matplotlib.pyplot.close(fig)

	def test_timesteps_validation(self):
		with self.assertRaises(ValueError):
			self.SIR.simulate(0, seeds={0: 'I'})

	def test_R0_none_without_interactions(self):
		model = NetworkEpiModel(self.G)
		model.add_spontaneous('I', 'R', self.mu)
		self.assertIsNone(model.R0())

	def test_large_spontaneous_rates_warn(self):
		model = NetworkEpiModel(self.G)
		model.add_interaction('S', 'I', 'I', self.beta)
		model.add_spontaneous('I', 'R', 0.8)
		model.add_spontaneous('I', 'D', 0.8)
		with self.assertWarns(UserWarning):
			model.simulate(10, seeds={0: 'I', 1: 'I'})
		self.assertTrue((model.values_.sum(axis=1) == 200).all())
		self.assertEqual(model.values_['I'].iloc[-1], 0)

	def test_spontaneous_into_active_state(self):
		# E -> I keeps the node active so it can go on to infect others
		SEIR = NetworkEpiModel(nx.complete_graph(20))
		SEIR.add_interaction('S', 'E', 'I', 1.0)
		SEIR.add_spontaneous('E', 'I', 1.0)
		SEIR.add_spontaneous('I', 'R', 0.0)
		SEIR.simulate(4, seeds={0: 'E'})
		self.assertEqual(SEIR.values_['I'].iloc[1], 1)
		self.assertEqual(SEIR.values_['E'].iloc[2], 19)
		self.assertEqual(SEIR.values_['I'].iloc[3], 20)

	def test_susceptible_fallback_for_SIS(self):
		SIS = NetworkEpiModel(self.G)
		SIS.add_interaction('S', 'I', 'I', self.beta)
		SIS.add_spontaneous('I', 'S', self.mu)
		SIS.simulate(20, seeds={0: 'I'})
		self.assertEqual(SIS.values_['S'].iloc[0], 199)
		self.assertTrue((SIS.values_.sum(axis=1) == 200).all())

	def test_no_seeds_and_undetectable_susceptible(self):
		self.SIR.simulate(5, seeds={})
		self.assertEqual(self.SIR.values_['S'].iloc[-1], 200)

		cyclic = NetworkEpiModel(self.G)
		cyclic.add_spontaneous('I', 'R', 0.1)
		cyclic.add_spontaneous('R', 'I', 0.1)
		with self.assertRaises(ValueError):
			cyclic.simulate(5, seeds={0: 'I'})

	def test_interaction_into_inactive_state(self):
		# S + I -> R + I: the target is not infectious, so it must not stay active
		model = NetworkEpiModel(nx.complete_graph(5))
		model.add_interaction('S', 'R', 'I', 1.0)
		model.add_spontaneous('I', 'R', 0.0)
		model.simulate(3, seeds={0: 'I'})
		# columns are S, R, I; "I" has no inflow either but is not a source
		self.assertEqual(list(model.values_.columns), ['S', 'R', 'I'])
		self.assertEqual(model.values_.iloc[0].tolist(), [4, 0, 1])
		self.assertEqual(model.values_.iloc[-1].tolist(), [0, 4, 1])


if __name__ == '__main__':
	unittest.main()
