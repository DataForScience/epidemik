MetaEpiModel: metapopulation models
===================================

:class:`~epidemik.MetaEpiModel` simulates an epidemic spreading across
several coupled sub-populations -- think US states, cities or districts.
Each sub-population is driven by its own :class:`~epidemik.EpiModel`, and
at every time step individuals are exchanged between sub-populations
according to a travel matrix.

Building a metapopulation model
-------------------------------

The constructor takes two :class:`pandas.DataFrame` arguments:

* ``travel_graph`` -- a square data frame whose rows and columns are
  labelled by sub-population identifiers. Entry ``(i, j)`` is the
  probability that an individual currently in sub-population ``i`` travels
  to sub-population ``j`` during a single time step. Each row must sum
  to 1.
* ``populations`` -- a data frame indexed by sub-population identifier
  containing at least one column with the total population of each
  sub-population (named ``"Population"`` by default; use the
  ``population`` argument to pick another column).

.. code-block:: python

   import pandas as pd
   from epidemik import MetaEpiModel

   travel = pd.DataFrame(
       {"A": [0.99, 0.10], "B": [0.01, 0.90]}, index=["A", "B"]
   )
   populations = pd.DataFrame(
       {"Population": [100_000, 10_000]}, index=["A", "B"]
   )

   SIR = MetaEpiModel(travel, populations)
   SIR.add_interaction("S", "I", "I", 0.3)
   SIR.add_spontaneous("I", "R", 0.1)

The constructor validates its inputs and raises :class:`ValueError` if the
travel matrix is not square, if its rows do not sum to 1, or if a
sub-population is missing from ``populations``.

Every transition added at the metapopulation level --
:meth:`~epidemik.MetaEpiModel.add_interaction`,
:meth:`~epidemik.MetaEpiModel.add_spontaneous`,
:meth:`~epidemik.MetaEpiModel.add_vaccination`,
:meth:`~epidemik.MetaEpiModel.add_birth_rate`,
:meth:`~epidemik.MetaEpiModel.add_death_rate` and
:meth:`~epidemik.MetaEpiModel.add_groups` -- is propagated to every
sub-population, so the disease is assumed to be the same everywhere.

Seeding and simulating
----------------------

The initial epidemic seeds are placed in a designated *seed sub-population*
via the ``seed_state`` argument of
:meth:`~epidemik.MetaEpiModel.simulate`. The remaining keyword arguments
give the initial size of each non-susceptible compartment in that
sub-population; the susceptible compartment is adjusted automatically so
that the total population of the seed sub-population is preserved.

.. code-block:: python

   SIR.simulate(timesteps=60, seed_state="A", I=10)

Internally, every time step alternates between a *spread* phase (every
sub-population's :class:`~epidemik.EpiModel` advances by one stochastic
step via :meth:`~epidemik.EpiModel.single_step`) and a *travel* phase
(individuals are redistributed across sub-populations according to
``travel_graph``). A ``seasonality`` array, if given, is forwarded to
every sub-population.

If ``seed_state`` is not provided, a :class:`~epidemik.utils.NotInitialized`
exception is raised. The first positional argument used to be called
``timestamp``; that spelling still works but emits a
:class:`DeprecationWarning`.

Inspecting the result
---------------------

After a simulation:

* ``SIR.compartments_`` holds the current compartment counts per
  sub-population (one row per sub-population, one column per
  compartment).
* ``SIR.get_state("A")`` returns the :class:`~epidemik.EpiModel` that
  drove sub-population ``"A"``, so you can inspect or plot its trajectory
  using the usual API. Its :attr:`~epidemik.EpiModel.values_` has
  ``timesteps + 1`` rows: the initial conditions followed by the state
  after each step.

  .. code-block:: python

     SIR.get_state("A").plot()

* :meth:`~epidemik.MetaEpiModel.plot` lays out one mini-panel per
  sub-population showing the trajectory of a compartment (``"I"`` by
  default), marks the time of the local peak with a vertical line and
  colours each panel by how many steps after the start that peak occurs.
* :meth:`~epidemik.MetaEpiModel.plot_peaks` produces a heat-map of the
  same normalized trajectories, which highlights how the epidemic peak
  diffuses across sub-populations.

Calling either plotting method before
:meth:`~epidemik.MetaEpiModel.simulate` raises
:class:`~epidemik.utils.NotInitialized`.

What is *not* supported
-----------------------

* :meth:`~epidemik.MetaEpiModel.integrate` raises
  :class:`NotImplementedError`: the metapopulation dynamics only make sense
  as a stochastic process.
* Age structure is not available at the metapopulation level.
* With :meth:`~epidemik.MetaEpiModel.add_groups`, ``populations`` only
  seeds the ``susceptible`` compartment of every sub-population. The
  initial size of the other groups (e.g. the vectors) can only be set in
  ``seed_state`` through the keyword arguments of
  :meth:`~epidemik.MetaEpiModel.simulate`.
