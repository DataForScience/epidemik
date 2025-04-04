"""
epidemik
========

epidemik is a Python package for the simulation of Compartmental Epidemic Models.

See https://www.github.com/DataForScience/epidemik for complete documentation.
"""

from .EpiModel import EpiModel
from .NetworkEpiModel import NetworkEpiModel
from .MetaEpiModel import MetaEpiModel

# Don't forget to update pyproject.toml
__version__ = "0.1.3.3"
