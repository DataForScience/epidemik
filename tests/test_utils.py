import unittest
import pandas as pd
from epidemik.utils import *
import logging


class ParametersTestCase(unittest.TestCase):
    def setUp(self):
        self.params = Parameters()
    
    def test_define_parameters(self):
        self.params['beta'] = '.2'
        self.assertEqual(self.params['beta'], 0.2)

    def test_compute_parameters(self):
        self.params['beta'] = '.2'
        self.assertEqual(self.params['beta'], 0.2)
    
    def test_math(self):
        self.params['beta'] = 0.2
        self.params['mu'] = 'beta/2'
        logging.warning(self.params.items())
        self.assertEqual(self.params['mu'], 0.1)

    def test_builtin_override(self):
        self.params['beta'] = 'str(1.3)'
        self.assertIsNone(self.params['beta'])