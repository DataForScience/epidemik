import unittest
from epidemik import EpiModel

class EpiModelTestCase(unittest.TestCase):
	def setUp(self):
		self.SIR = EpiModel()
		self.beta = 0.3
		self.mu = 0.1
		self.SIR.add_interaction('S', 'I', 'I', self.beta)
		self.SIR.add_spontaneous('I', 'R', self.mu)

	def test_R0(self):
		self.assertEqual(self.SIR.R0(), 3.0, 'incorrect R0')

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
