"""
epidemik
========

:mod:`epidemik` is a Python package for the simulation of compartmental
epidemic models.

The package exposes three main model classes:

* :class:`EpiModel` -- compartmental model in a homogeneously mixed
  population, supporting numerical integration and stochastic simulation.
* :class:`NetworkEpiModel` -- compartmental model on top of a contact
  network represented by a :class:`networkx.Graph`.
* :class:`MetaEpiModel` -- metapopulation model where several sub-populations
  are coupled together by a travel matrix.

The :mod:`epidemik.utils` submodule provides a small number of shared
helpers, including the :class:`~epidemik.utils.NotInitialized` exception.

See https://www.github.com/DataForScience/epidemik for the source code and
https://epidemik.readthedocs.io/ for the full documentation.
"""


from .EpiModel import EpiModel
from .NetworkEpiModel import NetworkEpiModel
from .MetaEpiModel import MetaEpiModel
from . import utils

__all__ = ["EpiModel", "NetworkEpiModel", "MetaEpiModel", "utils"]

__version__ = "0.1"
