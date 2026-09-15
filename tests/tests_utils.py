import os
import tempfile
import unittest

from epidemik.utils import Parameters, get_cache_directory, get_remote_path, OFFICIAL_REPO, epi_colors, EPI_COLORS


class ParametersTestCase(unittest.TestCase):
	def setUp(self):
		self.params = Parameters(beta=0.2)

	def test_numeric_and_string_values(self):
		self.params['mu'] = '.1'
		self.assertEqual(self.params['beta'], 0.2)
		self.assertEqual(self.params['mu'], 0.1)
		self.assertEqual(self.params.raw('mu'), 0.1)
		self.assertEqual(self.params.get('mu'), 0.1)
		self.assertIsNone(self.params.get('missing'))

	def test_expressions(self):
		self.params['mu'] = 'beta/2'
		self.params['nu'] = ' 2 * mu '
		self.params['k'] = 'exp(0) + sqrt(4) + pi - pi'
		self.assertAlmostEqual(self.params['mu'], 0.1)
		self.assertAlmostEqual(self.params['nu'], 0.2)
		self.assertAlmostEqual(self.params['k'], 3.0)
		self.assertEqual(self.params.raw('mu'), 'beta/2')
		self.assertEqual(self.params.resolve(), {'beta': 0.2, 'mu': 0.1, 'nu': 0.2, 'k': 3.0})

		# Redefining beta propagates to the derived parameters
		self.params['beta'] = 0.4
		self.assertAlmostEqual(self.params['mu'], 0.2)

	def test_errors(self):
		with self.assertRaises(KeyError):
			self.params['missing']

		self.params['bad'] = 'str(1.3)'   # builtins are not available
		with self.assertRaises(KeyError):
			self.params['bad']

		self.params['broken'] = '1 +'
		with self.assertRaises(ValueError):
			self.params['broken']

		self.params['undefined'] = 'gamma * 2'
		with self.assertRaises(KeyError):
			self.params['undefined']

		self.params['a'] = 'b'
		self.params['b'] = 'a'
		with self.assertRaises(ValueError):
			self.params['a']

		self.params['text'] = '"abc"'
		with self.assertRaises(ValueError):
			self.params['text']

		with self.assertRaises(TypeError):
			self.params['obj'] = object()

	def test_copy_and_repr(self):
		self.params['mu'] = 'beta/2'
		clone = self.params.copy()
		clone['beta'] = 1.0
		self.assertEqual(self.params['mu'], 0.1)
		self.assertEqual(clone['mu'], 0.5)
		self.assertEqual(repr(self.params), "Parameters(beta=0.2, mu='beta/2')")


class HelpersTestCase(unittest.TestCase):
	def test_cache_directory(self):
		with tempfile.TemporaryDirectory() as tmp:
			old = os.environ.get('EPIDEMIK_CACHE_DIR')
			os.environ['EPIDEMIK_CACHE_DIR'] = tmp
			try:
				self.assertEqual(get_cache_directory(), tmp)
			finally:
				if old is None:
					del os.environ['EPIDEMIK_CACHE_DIR']
				else:
					os.environ['EPIDEMIK_CACHE_DIR'] = old

		self.assertTrue(get_cache_directory().endswith('epidemik'))

	def test_remote_path(self):
		self.assertEqual(
			get_remote_path(),
			'https://raw.githubusercontent.com/DataForScience/epidemik/refs/heads/main/models/'
		)
		self.assertEqual(
			get_remote_path('https://github.com/someone/models'),
			'https://raw.githubusercontent.com/someone/models/refs/heads/main/models/'
		)
		self.assertTrue(OFFICIAL_REPO.startswith('https://github.com/'))

		with tempfile.TemporaryDirectory() as tmp:
			self.assertEqual(get_remote_path(tmp), tmp)

		with self.assertRaises(ValueError):
			get_remote_path('https://example.com/not-github')

	def test_colors(self):
		self.assertIs(epi_colors, EPI_COLORS)
		self.assertEqual(epi_colors['S'], '#51a7f9')
		self.assertEqual(epi_colors['Q'], '#f39019')


if __name__ == '__main__':
	unittest.main()
