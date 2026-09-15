import unittest

import matplotlib
matplotlib.use("Agg")

import numpy as np
from scipy.integrate import odeint

from epidemik import EpiModel
from epidemik.utils import NotInitialized


class EpiModelTestCase(unittest.TestCase):
	def setUp(self):
		np.random.seed(42)
		self.SIR = EpiModel()
		self.beta = 0.3
		self.mu = 0.1
		self.SIR.add_interaction('S', 'I', 'I', self.beta)
		self.SIR.add_spontaneous('I', 'R', self.mu)

	def test_R0(self):
		self.assertEqual(self.SIR.R0(), 3.0, 'incorrect R0')

	def test_R0_SEIR(self):
		SEIR = EpiModel()
		SEIR.add_interaction('S', 'E', 'I', self.beta)
		SEIR.add_spontaneous('E', 'I', 0.2)
		SEIR.add_spontaneous('I', 'R', self.mu)
		self.assertAlmostEqual(SEIR.R0(), self.beta/self.mu)

	def test_R0_SIRS(self):
		# Waning immunity R -> S must not break the R0 computation
		self.SIR.add_spontaneous('R', 'S', 0.01)
		self.assertAlmostEqual(self.SIR.R0(), self.beta/self.mu)

	def test_R0_without_interactions(self):
		model = EpiModel()
		model.add_spontaneous('I', 'R', self.mu)
		self.assertIsNone(model.R0())

	def test_edges(self):
		edges = list(self.SIR.transitions.edges(data=True))
		self.assertEqual(len(edges), 2)

		for edge in edges:
			data = edge[2]

			if "agent" in data:
				self.assertEqual(data['agent'], 'I')
				self.assertEqual(data['rate'], self.beta)
				self.assertEqual(edge[0], 'S')
				self.assertEqual(edge[1], 'I')
			else:
				self.assertEqual(data['rate'], self.mu)
				self.assertEqual(edge[0], 'I')
				self.assertEqual(edge[1], 'R')

	def test_integrate(self):
		self.SIR.integrate(100, S=990, I=10, R=0)
		values = self.SIR.values_

		self.assertEqual(values.shape, (100, 3))
		self.assertEqual(list(values.index[:2]), [1, 2])
		self.assertEqual(values.iloc[0].tolist(), [990, 10, 0])
		self.assertTrue(np.allclose(values.sum(axis=1), 1000))
		self.assertGreater(values['R'].iloc[-1], values['R'].iloc[0])

	def test_simulate(self):
		self.SIR.simulate(100, S=990, I=10, R=0)
		values = self.SIR.values_

		self.assertEqual(values.shape, (100, 3))
		self.assertEqual(list(values.index[:2]), [1, 2])
		# The first row holds the initial conditions, as in integrate()
		self.assertEqual(values.iloc[0].tolist(), [990, 10, 0])
		self.assertTrue((values.sum(axis=1) == 1000).all())
		self.assertGreater(values['R'].iloc[-1], 0)

	def test_t_min(self):
		self.SIR.integrate(10, t_min=5, S=990, I=10, R=0)
		self.assertEqual(list(self.SIR.values_.index), list(range(5, 15)))
		self.SIR.simulate(10, t_min=5, S=990, I=10, R=0)
		self.assertEqual(list(self.SIR.values_.index), list(range(5, 15)))

	def test_unknown_compartment_warns(self):
		with self.assertWarns(UserWarning):
			self.SIR.integrate(10, S=990, I=10, R=0, Q=5)
		self.assertEqual(list(self.SIR.values_.columns), ['S', 'I', 'R'])

		with self.assertWarns(UserWarning):
			self.SIR.simulate(10, S=990, I=10, R=0, Q=5)
		self.assertEqual(list(self.SIR.values_.columns), ['S', 'I', 'R'])

	def test_simulate_large_rates_do_not_crash(self):
		model = EpiModel()
		model.add_interaction('S', 'I', 'I', 5.0)
		model.add_spontaneous('I', 'R', 0.1)

		with self.assertWarns(UserWarning):
			model.simulate(20, S=990, I=10, R=0)

		self.assertTrue((model.values_.sum(axis=1) == 1000).all())

	def test_simulate_multiple_edges_same_target(self):
		# Two interactions S -> I with different agents must both count
		model = EpiModel()
		model.add_interaction('S', 'I', 'I', 0.2)
		model.add_interaction('S', 'I', 'A', 0.2)
		model.add_spontaneous('I', 'A', 0.5)
		model.add_spontaneous('A', 'R', 0.1)
		model.simulate(50, S=990, I=10, A=0, R=0)
		self.assertTrue((model.values_.sum(axis=1) == 1000).all())
		self.assertGreater(model.values_['R'].iloc[-1], 0)

	def test_single_step(self):
		self.SIR.single_step(S=990, I=10, R=0)
		self.assertEqual(list(self.SIR.values_.index), [1, 2])
		self.assertEqual(self.SIR.values_.iloc[0].tolist(), [990, 10, 0])

		# Without keyword arguments the simulation continues from the last state
		self.SIR.single_step()
		self.SIR.single_step()
		self.assertEqual(list(self.SIR.values_.index), [1, 2, 3, 4])
		self.assertTrue((self.SIR.values_.sum(axis=1) == 1000).all())

		# With keyword arguments the state is overridden before the step
		self.SIR.single_step(S=0, I=0, R=1000)
		self.assertEqual(self.SIR.values_.iloc[-1].tolist(), [0, 0, 1000])

	def test_getattr_and_getitem(self):
		with self.assertRaises(AttributeError):
			self.SIR.I

		with self.assertRaises(NotInitialized):
			self.SIR['I']

		self.SIR.integrate(10, S=990, I=10, R=0)

		self.assertEqual(self.SIR.I.shape, (10,))
		self.assertEqual(self.SIR['I'].shape, (10,))
		self.assertEqual(self.SIR[['S', 'R']].shape, (10, 2))
		self.assertIsNone(self.SIR['Q'])

		with self.assertRaises(AttributeError):
			self.SIR.Q

	def test_plot(self):
		with self.assertRaises(NotInitialized):
			self.SIR.plot(show=False)

		self.SIR.integrate(10, S=990, I=10, R=0)

		fig, ax = matplotlib.pyplot.subplots()
		result = self.SIR.plot(show=False, ax=ax, title='SIR')
		self.assertIs(result, ax)
		self.assertEqual(len(ax.lines), 3)
		matplotlib.pyplot.close(fig)

	def test_draw_model(self):
		fig, ax = matplotlib.pyplot.subplots()
		self.SIR.add_birth_rate(0.01)
		self.SIR.add_death_rate(0.01)
		self.SIR.draw_model(ax=ax, show=False)
		matplotlib.pyplot.close(fig)

	def test_vaccination(self):
		self.SIR.add_vaccination('S', 'V', 0.01, start=10)

		self.SIR.integrate(30, S=990, I=10, R=0, V=0)
		self.assertAlmostEqual(self.SIR.V.loc[10], 0, places=6)
		self.assertGreater(self.SIR.V.loc[12], 0)

		self.SIR.simulate(30, S=990, I=10, R=0, V=0)
		self.assertAlmostEqual(self.SIR.V.loc[10], 0, places=6)
		self.assertGreater(self.SIR.V.loc[30], 0)

	def test_demographics(self):
		self.SIR.add_birth_rate(0.01)
		self.SIR.add_death_rate(0.01)

		self.SIR.integrate(50, S=990, I=10, R=0)
		self.assertTrue(np.allclose(self.SIR.values_.sum(axis=1), 1000))

		self.SIR.simulate(50, S=990, I=10, R=0)
		self.assertEqual(self.SIR.values_.shape, (50, 3))

	def test_seasonality(self):
		season = 1 + 0.3*np.cos(2*np.pi*np.arange(365)/365)
		self.SIR.integrate(400, seasonality=season, S=990, I=10, R=0)
		self.assertTrue(np.allclose(self.SIR.values_.sum(axis=1), 1000))
		self.SIR.simulate(400, seasonality=season, S=990, I=10, R=0)
		self.assertTrue((self.SIR.values_.sum(axis=1) == 1000).all())

	def _age_model(self):
		contact = np.array([[10, 4], [4, 12]])
		population = np.array([1000, 2000])

		model = EpiModel()
		model.add_interaction('S', 'I', 'I', 0.02)
		model.add_spontaneous('I', 'R', self.mu)
		model.add_death_rate(0.001, comps=['R'])
		model.add_age_structure(contact, population)
		return model

	def test_age_structure_integrate(self):
		model = self._age_model()

		self.assertEqual(
			list(model.transitions.nodes()),
			['S_a', 'S_b', 'I_a', 'I_b', 'R_a', 'R_b']
		)
		# Death rate attached before the expansion is copied to every stratum
		self.assertEqual(model.transitions.nodes['R_a']['death'], 0.001)

		model.integrate(50, S=2990, I=10, R=0)
		first = model.values_.copy()

		self.assertEqual(list(model.values_.columns), ['S', 'I', 'R'])
		self.assertEqual(model.values_ages_.shape, (50, 6))
		self.assertEqual(model['S_a'].shape, (50,))
		self.assertEqual(model[['S_a', 'S_b']].shape, (50, 2))
		self.assertAlmostEqual(model.values_ages_['S_a'].iloc[0], 2990/3)

		# Integration is deterministic
		model.integrate(50, S=2990, I=10, R=0)
		self.assertTrue(np.allclose(first.values, model.values_.values))

	def test_age_structure_simulate(self):
		model = self._age_model()
		model.simulate(50, S=2990, I=10, R=0)

		self.assertEqual(list(model.values_.columns), ['S', 'I', 'R'])
		self.assertEqual(model.values_ages_.shape, (50, 6))
		self.assertEqual(model.values_.iloc[0].tolist(), [2990, 10, 0])

	def test_age_structure_validation(self):
		with self.assertRaises(ValueError):
			self.SIR.add_age_structure([[1, 2, 3], [4, 5, 6]], [10, 20])
		with self.assertRaises(ValueError):
			self.SIR.add_age_structure([[1, 2], [3, 4]], [10, 20, 30])

		self.SIR.add_age_structure([[1, 2], [3, 4]], [10, 20])
		with self.assertRaises(ValueError):
			self.SIR.add_groups({"all": list(self.SIR.transitions.nodes())})

	def _vector_model(self):
		bh, bv, gh, gv, mh, mv = 0.3, 0.3, 0.1, 0.1, 0.01, 0.1

		model = EpiModel()
		model.add_interaction("Sh", "Ih", "Iv", bh)
		model.add_interaction("Sv", "Iv", "Ih", bv)
		model.add_spontaneous("Ih", "Rh", gh)
		model.add_groups({"host": ["Sh", "Ih", "Rh"], "vector": ["Sv", "Iv"]})
		model.add_death_rate(mh, comps=["Sh", "Ih", "Rh"])
		model.add_death_rate(gv, comps=["Sv", "Iv"])
		model.add_birth_rate(mv, comps=["Sv"], from_group=True)

		return model, (bh, bv, gh, gv, mh, mv)

	def test_groups_match_reference_ode(self):
		model, params = self._vector_model()
		bh, bv, gh, gv, mh, mv = params

		def reference(state, t):
			Sh, Ih, Rh, Sv, Iv = state
			Nh = Sh + Ih + Rh
			Nv = Sv + Iv
			return [
				-bh*Sh*Iv/Nv - mh*Sh,
				bh*Sh*Iv/Nv - gh*Ih - mh*Ih,
				gh*Ih - mh*Rh,
				mv*Nv - bv*Sv*Ih/Nh - gv*Sv,
				bv*Sv*Ih/Nh - gv*Iv,
			]

		T = 100
		t = np.arange(1, T+1)
		ref = odeint(reference, [999, 1, 0, 1000, 10], t)

		model.integrate(T, Sh=999, Ih=1, Rh=0, Sv=1000, Iv=10)
		got = model.values_[["Sh", "Ih", "Rh", "Sv", "Iv"]].values

		self.assertTrue(np.allclose(got, ref, atol=1e-6))
		self.assertIsNone(model.R0())

	def test_groups_simulate(self):
		model, _ = self._vector_model()
		model.simulate(50, Sh=999, Ih=1, Rh=0, Sv=1000, Iv=10)
		self.assertEqual(model.values_.shape, (50, 5))
		self.assertFalse(model.values_.isnull().any().any())

	def test_groups_empty_group_does_not_produce_nan(self):
		model, _ = self._vector_model()
		# No vectors at all: the force of infection must be zero, not NaN
		model.integrate(20, Sh=999, Ih=1, Rh=0, Sv=0, Iv=0)
		self.assertFalse(model.values_.isnull().any().any())
		model.simulate(20, Sh=999, Ih=1, Rh=0, Sv=0, Iv=0)
		self.assertFalse(model.values_.isnull().any().any())

	def test_groups_validation(self):
		with self.assertRaises(ValueError):
			self.SIR.add_groups({"host": ["S", "I"]})  # R unassigned

		with self.assertRaises(ValueError):
			self.SIR.add_groups({"a": ["S", "I", "R"], "b": ["R"]})

		self.SIR.add_groups({"host": ["S", "I", "R"]})

		with self.assertRaises(ValueError):
			self.SIR.add_age_structure([[1]], [10])

		# Compartments added after add_groups must be assigned before running
		self.SIR.add_spontaneous('R', 'D', 0.01)
		with self.assertRaises(ValueError):
			self.SIR.integrate(10, S=990, I=10, R=0, D=0)

	def test_birth_from_group_requires_groups(self):
		with self.assertRaises(ValueError):
			self.SIR.add_birth_rate(0.1, comps=['S'], from_group=True)

	def test_repr(self):
		text = repr(self.SIR)
		self.assertIn('3 compartments and 2 transitions', text)
		self.assertIn('S + I = I', text)
		self.assertIn('R0=3.00', text)

	def test_constructor_compartments(self):
		model = EpiModel(['S', 'I', 'R'])
		self.assertEqual(list(model.transitions.nodes()), ['S', 'I', 'R'])
		model.add_interaction('S', 'I', 'I', self.beta)
		model.add_spontaneous('I', 'R', self.mu)
		model.integrate(5, S=990, I=10, R=0)
		self.assertEqual(list(model.values_.columns), ['S', 'I', 'R'])

	def test_timesteps_validation(self):
		with self.assertRaises(ValueError):
			self.SIR.integrate(0, S=990, I=10, R=0)
		with self.assertRaises(ValueError):
			self.SIR.simulate(0, S=990, I=10, R=0)

	def test_simulate_deaths_never_negative(self):
		# High death rate together with fast transitions used to drive
		# compartments negative because deaths were sampled independently.
		model = EpiModel()
		model.add_interaction('S', 'I', 'I', 0.9)
		model.add_spontaneous('I', 'R', 0.9)
		model.add_death_rate(0.5)
		model.add_birth_rate(0.01, comps=['S'])
		model.simulate(60, S=500, I=500, R=0)
		self.assertTrue((model.values_ >= 0).all().all())
		self.assertLess(model.values_.sum(axis=1).iloc[-1], 1000)

	def test_seasonality_length(self):
		season = np.array([1.0, 0.0])   # transmission only on odd time steps
		self.SIR.integrate(10, seasonality=season, S=990, I=10, R=0)
		self.assertTrue(np.allclose(self.SIR.values_.sum(axis=1), 1000))
		self.SIR.simulate(10, seasonality=season, S=990, I=10, R=0)
		self.assertTrue((self.SIR.values_.sum(axis=1) == 1000).all())

	def test_R0_singular(self):
		# SI model: V is singular, so R0 is undefined
		SI = EpiModel()
		SI.add_interaction('S', 'I', 'I', self.beta)
		self.assertIsNone(SI.R0())

	def test_R0_ignores_vaccination_and_demographics(self):
		self.SIR.add_vaccination('S', 'V', 0.5, start=1)
		self.SIR.add_death_rate(0.5)
		self.assertAlmostEqual(self.SIR.R0(), self.beta/self.mu)

	def test_repr_extras(self):
		self.SIR.add_vaccination('S', 'V', 0.01, start=10)
		self.SIR.add_birth_rate(0.001, comps=['S'])
		self.SIR.add_death_rate(0.002)
		text = repr(self.SIR)
		self.assertIn('S -> V 0.010000 starting at 10 days', text)
		self.assertIn('birth -> S 0.001000', text)
		self.assertIn('R -> death 0.002000', text)

		model, _ = self._vector_model()
		text = repr(model)
		self.assertIn('Groups:', text)
		self.assertIn('host: Sh, Ih, Rh', text)
		self.assertIn('birth -> Sv 0.100000 (per group total)', text)
		self.assertNotIn('R0=', text)

	def test_age_structure_vaccination_and_single_step(self):
		contact = np.array([[10, 4], [4, 12]])
		population = np.array([1000, 2000])

		model = EpiModel()
		model.add_interaction('S', 'I', 'I', 0.02)
		model.add_spontaneous('I', 'R', self.mu)
		model.add_vaccination('S', 'V', 0.01, start=2)
		model.add_age_structure(contact, population)
		# A compartment added after the expansion must still be aggregated
		model.add_spontaneous('R_a', 'D_a', 0.01)
		model.add_spontaneous('R_b', 'D_b', 0.01)

		self.assertEqual(model.transitions.edges['S_a', 'V_a', 0]['start'], 2)

		model.integrate(10, S=2990, I=10, R=0)
		self.assertEqual(list(model.values_.columns), ['S', 'I', 'R', 'V', 'D'])
		self.assertGreater(model.V.iloc[-1], 0)
		self.assertGreater(model.D.iloc[-1], 0)

		# single_step continues from the existing trajectory (10 rows so far)
		model.single_step(S=2990, I=10, R=0, V=0, D=0)
		model.single_step()
		self.assertEqual(model.values_.shape, (12, 5))
		self.assertEqual(model.values_ages_.shape, (12, 10))
		self.assertEqual(list(model.values_.index[-3:]), [10, 11, 12])
		self.assertTrue((model.values_.iloc[-2:].sum(axis=1) == 3000).all())

		# Expanded names can be used directly as initial conditions
		model.integrate(5, S_a=1000, S_b=1990, I_b=10)
		self.assertEqual(model.values_ages_['S_a'].iloc[0], 1000)
		self.assertEqual(model.values_['I'].iloc[0], 10)

		with self.assertWarns(UserWarning):
			model.integrate(5, S=2990, I=10, Q=1)

	def test_age_structure_too_many_groups(self):
		n = 27
		with self.assertRaises(ValueError):
			self.SIR.add_age_structure(np.ones((n, n)), np.ones(n))

	def test_plot_defaults(self):
		self.SIR.integrate(10, S=990, I=10, R=0)
		ax = self.SIR.plot(show=False, normed=False)
		self.assertEqual(len(ax.lines), 3)
		matplotlib.pyplot.close('all')

	def test_draw_model_default_axes(self):
		model = EpiModel()
		model.add_interaction('S', 'I', 'I', 0.2)
		model.add_interaction('S', 'I', 'A', 0.1)
		model.add_spontaneous('I', 'A', 0.1)
		model.add_spontaneous('A', 'R', 0.1)
		model.draw_model(show=False)
		matplotlib.pyplot.close('all')

	def test_dunder_lookup_does_not_hit_values(self):
		import copy
		self.SIR.integrate(5, S=990, I=10, R=0)
		clone = copy.deepcopy(self.SIR)
		self.assertTrue(clone.values_.equals(self.SIR.values_))

	def test_plot_and_draw_with_show(self):
		import warnings
		self.SIR.integrate(5)   # no initial conditions: all-zero first row
		with warnings.catch_warnings():
			warnings.simplefilter('ignore')
			self.SIR.plot(show=True)
			self.SIR.draw_model(show=True)
		matplotlib.pyplot.close('all')


if __name__ == '__main__':
	unittest.main()
