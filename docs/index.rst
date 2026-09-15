epidemik documentation
======================

``epidemik`` is a Python package for building, exploring and visualising
compartmental epidemic models. It supports three flavours of model:

* :class:`~epidemik.EpiModel` -- a classical homogeneously mixed
  compartmental model that can be either integrated as a system of ODEs or
  simulated stochastically.
* :class:`~epidemik.NetworkEpiModel` -- a stochastic compartmental model
  that lives on top of a contact network supplied as a
  :class:`networkx.Graph`.
* :class:`~epidemik.MetaEpiModel` -- a metapopulation model that couples
  several :class:`~epidemik.EpiModel` instances together through a travel
  matrix.

All three classes share the same fluent API for declaring compartments and
transitions, so moving from a homogeneous mixing model to a network or
metapopulation variant is mostly a matter of swapping the class.

.. code-block:: python

   from epidemik import EpiModel

   SIR = EpiModel()
   SIR.add_interaction("S", "I", "I", 0.2)
   SIR.add_spontaneous("I", "R", 0.1)
   SIR.integrate(365, S=9_990, I=10, R=0)
   SIR.plot()

.. toctree::
   :maxdepth: 2
   :caption: User guide

   installation
   quickstart
   usage/epimodel
   usage/networkepimodel
   usage/metaepimodel

.. toctree::
   :maxdepth: 2
   :caption: Reference

   api

.. toctree::
   :maxdepth: 1
   :caption: Project

   changelog


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
