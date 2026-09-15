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
				self.assertEqual(data['rate'], 'rate1')
				self.assertEqual(self.SIR.params['rate1'], self.beta)
				self.assertEqual(edge[0], 'S')
				self.assertEqual(edge[1], 'I')
			else:
				self.assertEqual(data['rate'], 'rate2')
				self.assertEqual(self.SIR.params['rate2'], self.mu)
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
		# births proportional to the total population, all into S; deaths everywhere
		self.SIR.add_birth_rate(0.01, comps=['S'])
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
		death_key = model.transitions.nodes['R_a']['death']
		self.assertEqual(model.params[death_key], 0.001)
		# Rescaled interaction rates stay tied to the original parameter
		self.assertEqual(model.params.raw('rate1_ab'), '(rate1)*4.0')
		self.assertAlmostEqual(model.params['rate1_ab'], 0.08)

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
		self.assertIn('# Epidemic Model with 3 compartments and 2 transitions', text)
		self.assertIn('  rate1: 0.3', text)
		self.assertIn('  - S + I = I rate1', text)
		self.assertIn('  - I -> R rate2', text)
		self.assertIn('# R0=3.00', text)

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
		self.assertIn('  - S -> V rate3 starting at 10 days', text)
		self.assertIn('  - {compartment: S, birth: rate4}', text)
		self.assertIn('  - {compartment: R, death: rate5}', text)

		model, _ = self._vector_model()
		text = repr(model)
		self.assertIn('Groups:', text)
		self.assertIn('  host: [Sh, Ih, Rh]', text)
		self.assertIn('  - {compartment: Sv, birth: rate6, from_group: true}', text)
		self.assertIn('Compartments: [Sh, Ih, Iv, Sv, Rh]', text)
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

	# ------------------------------------------------------------------
	# Named parameters
	# ------------------------------------------------------------------
	def test_named_parameters_and_expressions(self):
		model = EpiModel()
		model.add_interaction('S', 'I', 'I', beta=0.3)
		model.add_spontaneous('I', 'R', mu='beta/3')
		self.assertEqual(model.params.raw('mu'), 'beta/3')
		self.assertAlmostEqual(model.R0(), 3.0)

		# Redefining a parameter changes every rate derived from it
		model.define_parameters(beta=0.6)
		self.assertAlmostEqual(model.R0(), 3.0)
		model.define_parameters(mu=0.1)
		self.assertAlmostEqual(model.R0(), 6.0)

		# Positional strings: existing parameter name or expression
		model.add_vaccination('S', 'V', 'mu', start=5)
		self.assertEqual(model.transitions.edges['S', 'V', 0]['rate'], 'mu')
		model.add_spontaneous('R', 'S', '2*mu')
		self.assertEqual(model.params.raw('rate1'), '2*mu')
		self.assertAlmostEqual(model.params['rate1'], 0.2)

		model.integrate(10, S=990, I=10)
		model.simulate(10, S=990, I=10)

	def test_rate_argument_validation(self):
		with self.assertRaises(TypeError):
			self.SIR.add_spontaneous('R', 'S', 0.1, omega=0.1)
		with self.assertRaises(TypeError):
			self.SIR.add_spontaneous('R', 'S')
		with self.assertRaises(TypeError):
			self.SIR.add_spontaneous('R', 'S', a=0.1, b=0.2)
		with self.assertRaises(TypeError):
			self.SIR.add_vaccination('S', 'V', 0.1)
		with self.assertRaises(TypeError):
			self.SIR.add_viral_generation('I', 'V', 0.1)
		with self.assertRaises(TypeError):
			self.SIR.add_viral_generation('I', 'V', a=0.1)

	def test_numpy_scalar_rate_is_not_aliased(self):
		model = EpiModel()
		model.add_interaction('S', 'I', 'I', np.array(0.3))
		model.add_spontaneous('I', 'R', np.float64(0.1))
		model.integrate(5, S=990, I=10)
		self.assertEqual(model.params.raw('rate1'), 0.3)
		self.assertEqual(model.params.raw('rate2'), 0.1)

	def test_undefined_parameter_in_expression(self):
		self.SIR.add_spontaneous('R', 'S', 'omega/2')
		with self.assertRaises(KeyError):
			self.SIR.integrate(5, S=990, I=10)
		self.assertIsNone(self.SIR.R0())

	# ------------------------------------------------------------------
	# Random number generation
	# ------------------------------------------------------------------
	def test_seed_reproducibility(self):
		runs = []
		for _ in range(2):
			model = EpiModel(seed=123)
			model.add_interaction('S', 'I', 'I', self.beta)
			model.add_spontaneous('I', 'R', self.mu)
			model.simulate(50, S=990, I=10)
			runs.append(model.values_)
		self.assertTrue(runs[0].equals(runs[1]))

		rng = np.random.default_rng(7)
		model = EpiModel(rng=rng)
		self.assertIs(model.rng, rng)

	def test_reset(self):
		self.SIR.integrate(5, S=990, I=10)
		self.SIR.reset()
		with self.assertRaises(NotInitialized):
			self.SIR['I']

	# ------------------------------------------------------------------
	# Demographics modes
	# ------------------------------------------------------------------
	def test_birth_modes(self):
		fixed = EpiModel(['S'])
		fixed.add_birth_rate(10, comps=['S'], fixed=True)
		fixed.integrate(5, S=100)
		self.assertTrue(np.allclose(fixed.S.values, 100 + 10*np.arange(5)))

		per_capita = EpiModel(['S'])
		per_capita.add_birth_rate(b=0.1, comps=['S'], global_rate=False)
		per_capita.integrate(5, S=100)
		self.assertTrue(np.allclose(per_capita.S.values, 100*np.exp(0.1*np.arange(5)), rtol=1e-4))
		self.assertIn('  - {compartment: S, birth: b, global: false}', repr(per_capita))
		self.assertIn('  - {compartment: S, birth: rate1, fixed: true}', repr(fixed))

		# Births proportional to the total population, all into S
		total = EpiModel()
		total.add_spontaneous('I', 'R', 0.0)
		total.add_birth_rate(0.1, comps=['S'])
		total.integrate(5, S=100, I=100)
		self.assertTrue(np.allclose(total.values_.sum(axis=1).values, 200*np.exp(0.1*np.arange(5)), rtol=1e-4))
		self.assertTrue((total.I.values == 100).all())

		# Old calling convention: compartments first, named rate
		compat = EpiModel()
		compat.add_interaction('S', 'I', 'I', beta=0.3)
		compat.add_birth_rate('S', b=0.3)
		compat.add_birth_rate(['I'], b2=0.1)
		self.assertEqual(compat.transitions.nodes['S']['birth'], 'b')
		self.assertEqual(compat.transitions.nodes['I']['birth'], 'b2')
		self.assertEqual(compat.params['b'], 0.3)

		for model in (fixed, per_capita, total):
			model.simulate(5, S=100, I=100)
			self.assertTrue((model.values_ >= 0).all().all())

	def test_death_modes(self):
		fixed = EpiModel(['S'])
		fixed.add_death_rate(5, comps=['S'], fixed=True)
		fixed.integrate(5, S=100)
		self.assertTrue(np.allclose(fixed.S.values, 100 - 5*np.arange(5)))

		per_capita = EpiModel(['S'])
		per_capita.add_death_rate(d=0.1)
		per_capita.integrate(5, S=100)
		self.assertTrue(np.allclose(per_capita.S.values, 100*np.exp(-0.1*np.arange(5)), rtol=1e-4))

		total = EpiModel(['S'])
		total.add_spontaneous('I', 'R', 0.0)
		total.add_death_rate(0.1, comps=['I'], global_rate=True)
		total.integrate(2, S=100, I=100)
		# dI/dt = -0.1 * (S + I) with S = 100: I(1) = 200 exp(-0.1) - 100
		self.assertAlmostEqual(total.I.iloc[1], 200*np.exp(-0.1) - 100, places=3)
		self.assertIn('  - {compartment: I, death: rate2, global: true}', repr(total))
		self.assertIn('  - {compartment: S, death: rate1, fixed: true}', repr(fixed))

		for model in (fixed, per_capita, total):
			model.simulate(20, S=100, I=100)
			self.assertTrue((model.values_ >= 0).all().all())

	# ------------------------------------------------------------------
	# Within-host models
	# ------------------------------------------------------------------
	def _viral_model(self):
		model = EpiModel()
		model.define_parameters(s=1e5, d=0.01, beta=2e-7, delta=0.5, p=100, c=2)
		model.add_birth_rate('s', comps=['T'], fixed=True)
		model.add_death_rate('d', comps=['T'])
		model.add_interaction('T', 'I', 'V', 'beta', norm=False)
		model.add_viral_generation('I', 'V', 'delta', 'p')
		model.add_death_rate('c', comps=['V'])
		return model

	def test_viral_generation_matches_reference(self):
		model = self._viral_model()
		self.assertIsNone(model.R0())

		def reference(state, t, s=1e5, d=0.01, beta=2e-7, delta=0.5, p=100, c=2):
			T, I, V = state
			return [s - d*T - beta*T*V, beta*T*V - delta*I, p*I - c*V]

		t = np.arange(1, 61)
		ref = odeint(reference, [1e7, 0, 10], t)
		model.integrate(60, T=1e7, I=0, V=10)
		self.assertTrue(np.allclose(model.values_[['T', 'I', 'V']].values, ref, rtol=1e-5, atol=1e-3))

		self.assertIn('  - T + V = I beta unnormalized', repr(model))
		self.assertIn('  - I => V delta p', repr(model))

		model.simulate(10, T=1000, I=10, V=100)
		self.assertTrue((model.values_ >= 0).all().all())

		named = EpiModel()
		named.add_viral_generation('I', 'V', delta=0.5, p=100)
		self.assertEqual(named.params.raw('delta'), 0.5)
		self.assertEqual(named.params.raw('p'), 100)

	# ------------------------------------------------------------------
	# Persistence
	# ------------------------------------------------------------------
	def _models_dir(self):
		import os
		return os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'models')

	def test_save_load_roundtrip(self):
		import tempfile

		model = EpiModel(name='Test model')
		model.add_interaction('S', 'E', 'I', beta=0.3)
		model.add_spontaneous('E', 'I', sigma='beta/2')
		model.add_spontaneous('I', 'R', mu=0.1)
		model.add_vaccination('S', 'V', nu=0.01, start=30)
		model.add_birth_rate(b=0.001, comps=['S'])
		model.add_death_rate(d=0.001)
		model.add_death_rate(5, comps=['R'], fixed=True)

		with tempfile.TemporaryDirectory() as tmp:
			path = tmp + '/model.yaml'
			model.save_model(path)
			loaded = EpiModel.load_model(path, seed=1)

		self.assertEqual(loaded.name, 'Test model')
		self.assertEqual(repr(loaded), repr(model))
		self.assertEqual(loaded.params.raw('sigma'), 'beta/2')

		model.integrate(50, S=990, I=10)
		loaded.integrate(50, S=990, I=10)
		self.assertTrue(np.allclose(model.values_.values, loaded.values_.values))

		vector, _ = self._vector_model()
		viral = self._viral_model()
		for original in (vector, viral):
			with tempfile.TemporaryDirectory() as tmp:
				original.save_model(tmp + '/m.yaml')
				loaded = EpiModel.load_model(tmp + '/m.yaml')
			self.assertEqual(repr(loaded), repr(original))

	def test_load_model_library(self):
		import os

		names = sorted(f for f in os.listdir(self._models_dir()) if f.endswith('.yaml'))
		self.assertEqual(len(names), 6)

		for name in names:
			model = EpiModel.load_model(os.path.join(self._models_dir(), name))
			self.assertGreater(model.transitions.number_of_edges(), 0)
			model.integrate(5, **{list(model.transitions.nodes())[0]: 990, 'I': 10, 'Is': 10})

		sir = EpiModel.load_model(os.path.join(self._models_dir(), 'SIR.yaml'))
		self.assertEqual(sir.name, 'SIR')
		self.assertAlmostEqual(sir.R0(), 2.0)

	def test_load_model_errors_and_literals(self):
		import tempfile

		with tempfile.TemporaryDirectory() as tmp:
			path = tmp + '/bad.yaml'

			with open(path, 'wt') as f:
				f.write('Transitions:\n  - S ++ I\n')
			with self.assertRaises(ValueError):
				EpiModel.load_model(path)

			with open(path, 'wt') as f:
				f.write('Parameters:\n  beta: 0.3\nTransitions:\n  - I -> R mu\n')
			with self.assertRaises(ValueError):
				EpiModel.load_model(path)

			with open(path, 'wt') as f:
				f.write('transitions:\n  - S + I = I 0.3\n  - I -> R 0.1\n  - I => V 0.5 100\n')
			model = EpiModel.load_model(path)
			self.assertEqual(model.params.raw('rate1'), 0.3)
			self.assertEqual(model.params.raw('rate2'), 0.1)
			self.assertEqual(model.params.raw('rate4'), 100.0)

			with open(path, 'wt') as f:
				f.write('')
			self.assertEqual(EpiModel.load_model(path).transitions.number_of_nodes(), 0)

	def test_list_and_download_models(self):
		import os
		import tempfile

		with tempfile.TemporaryDirectory() as cache:
			os.environ['EPIDEMIK_CACHE_DIR'] = cache
			try:
				listed = EpiModel.list_models(repo=self._models_dir())
				self.assertIn('SIR.yaml', listed)
				self.assertEqual(len(listed), 6)

				path = EpiModel.download_model('SIR.yaml', repo=self._models_dir(), load_model=False)
				self.assertEqual(path, os.path.join(cache, 'SIR.yaml'))
				self.assertTrue(os.path.exists(path))

				# The cached copy is reused unless force=True
				with open(path, 'wt') as f:
					f.write('Transitions:\n  - S + I = I 0.5\n')
				cached = EpiModel.download_model('SIR.yaml', repo=self._models_dir())
				self.assertEqual(cached.params['rate1'], 0.5)
				fresh = EpiModel.download_model('SIR.yaml', repo=self._models_dir(), force=True, seed=1)
				self.assertEqual(fresh.params['beta'], 0.2)

				# Files in the cache are listed too (no duplicates)
				fresh.name = 'copy'
				fresh.save_model(os.path.join(cache, 'MyModel.yaml'))
				listed = EpiModel.list_models(repo=self._models_dir())
				self.assertIn('MyModel.yaml', listed)
				self.assertEqual(listed.count('SIR.yaml'), 1)

				# GitHub repositories are fetched with urlretrieve (mocked here)
				import shutil
				import sys
				EpiModelModule = sys.modules['epidemik.EpiModel']
				fetched = []

				def fake_urlretrieve(url, path):
					fetched.append(url)
					shutil.copyfile(os.path.join(self._models_dir(), os.path.basename(url)), path)

				original = EpiModelModule.urlretrieve
				EpiModelModule.urlretrieve = fake_urlretrieve
				try:
					listed = EpiModel.list_models()
					seir = EpiModel.download_model('SEIR.yaml')
				finally:
					EpiModelModule.urlretrieve = original

				self.assertIn('SEIR.yaml', listed)
				self.assertAlmostEqual(seir.R0(), 2.0)
				self.assertEqual(fetched, [
					'https://raw.githubusercontent.com/DataForScience/epidemik/refs/heads/main/models/model_list.txt',
					'https://raw.githubusercontent.com/DataForScience/epidemik/refs/heads/main/models/SEIR.yaml',
				])
			finally:
				del os.environ['EPIDEMIK_CACHE_DIR']

	# ------------------------------------------------------------------
	# Review findings
	# ------------------------------------------------------------------
	def test_agent_only_compartment(self):
		model = EpiModel()
		model.add_interaction('S', 'I', 'W', 0.3)   # environmental reservoir
		model.add_spontaneous('I', 'R', 0.1)
		self.assertEqual(list(model.transitions.nodes()), ['S', 'I', 'W', 'R'])
		model.integrate(5, S=990, I=10, W=100)
		model.simulate(5, S=990, I=10, W=100)
		self.assertEqual(model.values_['W'].iloc[-1], 100)

	def test_R0_leaky_vaccine_is_consistent(self):
		def build(waning):
			model = EpiModel()
			model.add_interaction('S', 'I', 'I', 0.3)
			model.add_spontaneous('I', 'R', 0.1)
			model.add_interaction('V', 'I', 'I', 0.15)
			model.add_vaccination('S', 'V', 0.01, start=10)
			if waning:
				model.add_spontaneous('R', 'S', 0.05)
			return model

		self.assertAlmostEqual(build(False).R0(), 3.0)
		self.assertAlmostEqual(build(True).R0(), 3.0)
		self.assertEqual(build(True)._get_susceptible(), {'S'})
