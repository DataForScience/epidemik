Quick start
===========

This page walks through the canonical SIR example using
:class:`~epidemik.EpiModel`. By the end of it you will have built a model,
inspected its structure, computed :math:`R_0`, and explored its dynamics
both deterministically and stochastically.

Building an SIR model
---------------------

We start by creating an empty :class:`~epidemik.EpiModel` and adding the
two transitions that define the SIR dynamics. The interaction
``S + I -> I + I`` describes infections at rate :math:`\beta`, and the
spontaneous transition ``I -> R`` describes recoveries at rate
:math:`\mu`.

.. code-block:: python

   from epidemik import EpiModel

   beta = 0.2
   mu = 0.1

   SIR = EpiModel(seed=42)
   SIR.add_interaction("S", "I", "I", beta=beta)
   SIR.add_spontaneous("I", "R", mu=mu)

Rates are stored as named parameters (``beta`` and ``mu`` here; a plain
number would create anonymous ``rate1``, ``rate2``, ... parameters) and
can be expressions of one another, e.g. ``mu="beta/2"``. A printable,
YAML-formatted description of the model is available via ``print``:

.. code-block:: pycon

   >>> print(SIR)
   # Epidemic Model with 3 compartments and 2 transitions:

   Compartments: [S, I, R]

   Parameters:
     beta: 0.2
     mu: 0.1

   Transitions:
     - S + I = I beta
     - I -> R mu

   # R0=2.00

The same text is what :meth:`~epidemik.EpiModel.save_model` writes and
:meth:`~epidemik.EpiModel.load_model` reads; ready-made models can be
fetched with :meth:`~epidemik.EpiModel.download_model`.

The corresponding compartment graph can be drawn with
:meth:`~epidemik.EpiModel.draw_model`:

.. code-block:: python

   SIR.draw_model()

The basic reproduction number :math:`R_0` is computed from the
next-generation matrix:

.. code-block:: python

   SIR.R0()  # -> 2.0

Deterministic dynamics
----------------------

To follow the mean-field evolution of the model, call
:meth:`~epidemik.EpiModel.integrate`. The first argument is the number of
time steps, and the remaining keyword arguments give the initial size of
each compartment.

.. code-block:: python

   N = 10_000
   I0 = 10
   SIR.integrate(365, S=N - I0, I=I0, R=0)

The full trajectory is available in :attr:`~epidemik.EpiModel.values_` as a
:class:`pandas.DataFrame` with one row per time step (the first row holds
the initial conditions), and individual compartments can be obtained
either as attributes or by indexing:

.. code-block:: python

   SIR.values_.head()
   SIR.I                # equivalent to SIR.values_["I"]
   SIR[["S", "R"]]      # multiple compartments at once

A convenience plotting helper is provided as well:

.. code-block:: python

   SIR.plot(title="Deterministic SIR")

Stochastic dynamics
-------------------

To sample one realization of the corresponding stochastic process, use
:meth:`~epidemik.EpiModel.simulate` with the same signature. The output
has the same shape and time index as the one produced by
:meth:`~epidemik.EpiModel.integrate`, so the two can be compared directly:

.. code-block:: python

   SIR.simulate(365, S=N - I0, I=I0, R=0)
   SIR.plot(title="Stochastic SIR")

You can also advance a stochastic simulation one step at a time using
:meth:`~epidemik.EpiModel.single_step`: the first call takes the initial
conditions and every later call appends one more row to
:attr:`~epidemik.EpiModel.values_`. This is the building block used
internally by :class:`~epidemik.MetaEpiModel` to couple sub-populations
together.

Where to go next
----------------

* :doc:`usage/epimodel` documents the full :class:`~epidemik.EpiModel`
  API, including named parameters, demographics, vaccination, age
  structure, seasonality, multi-group (vector-borne) and within-host
  models, and saving/loading model files.
* :doc:`usage/networkepimodel` describes how to run the same kind of
  model on top of an arbitrary contact network.
* :doc:`usage/metaepimodel` shows how to spread an epidemic across many
  coupled sub-populations.
* :doc:`api` is the auto-generated API reference for the whole package.
