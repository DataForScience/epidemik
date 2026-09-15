EpiModel: homogeneously mixed models
====================================

:class:`~epidemik.EpiModel` is the workhorse of the ``epidemik`` package.
This page documents every feature it currently supports, with short
examples for each.

Defining the dynamics
---------------------

A model is built up from two basic types of transition:

* **Interaction** -- ``source + agent -> target + agent`` at a given rate.
  This is the canonical mass-action contagion term: its contribution to
  the dynamics is ``rate * source * agent / N``.
* **Spontaneous** -- ``source -> target`` at a given rate, independent of
  any other compartment (recoveries, latency, waning immunity, ...).

.. code-block:: python

   from epidemik import EpiModel

   SEIR = EpiModel()
   SEIR.add_interaction("S", "E", "I", 0.3)   # S + I -> E + I
   SEIR.add_spontaneous("E", "I", 0.2)         # E -> I
   SEIR.add_spontaneous("I", "R", 0.1)         # I -> R

You can declare compartments explicitly when constructing the model, but
it is also fine to rely on the implicit creation triggered by the calls
above. Several transitions may share the same source and target (for
example two interactions driven by different agents); their contributions
are added together.

``print(model)`` gives a textual summary of every transition, of any
birth/death rates and groups, and of :math:`R_0`;
:meth:`~epidemik.EpiModel.draw_model` draws the compartment graph.

Computing :math:`R_0`
---------------------

The basic reproduction number is computed using the next-generation matrix
approach, so it works for any model topology -- SIR, SEIR, SIS, SIRS,
models with several infectious compartments, and so on:

.. code-block:: python

   SEIR.R0()

Vaccination transitions and demographics are ignored in the computation.
If the model has no infectious dynamics (or :math:`R_0` cannot be derived
for some other reason), :meth:`~epidemik.EpiModel.R0` returns ``None``.
Multi-group (vector-borne) models are not supported yet and also return
``None``.

Running the dynamics
--------------------

There are two ways to explore the dynamics of an :class:`~epidemik.EpiModel`:

* :meth:`~epidemik.EpiModel.integrate` solves the corresponding ODE system
  numerically using :func:`scipy.integrate.odeint`.
* :meth:`~epidemik.EpiModel.simulate` runs a discrete-time stochastic
  simulation: at every step, the individuals in each compartment are
  redistributed with a multinomial draw whose probabilities are the
  transition rates.

Both methods take the number of time steps as a positional argument and
the initial size of each compartment as keyword arguments:

.. code-block:: python

   SEIR.integrate(365, S=9_990, I=10, R=0)
   SEIR.simulate(365, S=9_990, I=10, R=0)

Results are stored as a :class:`pandas.DataFrame` in
:attr:`~epidemik.EpiModel.values_`, indexed by time. Both methods produce
exactly ``timesteps`` rows, the first of which (at ``t_min``, which
defaults to 1) holds the initial conditions, so deterministic and
stochastic trajectories can be compared directly. Compartments that were
not mentioned start at zero; initial conditions for names that are not
compartments of the model are ignored with a warning.

Because :meth:`~epidemik.EpiModel.simulate` uses the rates as per-step
transition probabilities, they should be small compared to 1. If the
probabilities out of some compartment add up to more than 1, they are
renormalized and a warning is emitted.

Stepping manually
-----------------

:meth:`~epidemik.EpiModel.single_step` advances a stochastic simulation
one step at a time. The first call takes the initial conditions; later
calls continue from the last row of :attr:`~epidemik.EpiModel.values_`
unless new compartment sizes are given, in which case they replace the
current state:

.. code-block:: python

   SEIR.single_step(S=9_990, I=10, R=0)   # rows at t=1 (initial) and t=2
   SEIR.single_step()                     # row at t=3
   SEIR.single_step(S=9_000, I=1_000)     # override the state, then step

This is the building block :class:`~epidemik.MetaEpiModel` uses to couple
sub-populations together.

Vaccination
-----------

To model a time-gated vaccination campaign, use
:meth:`~epidemik.EpiModel.add_vaccination`. The transition behaves like a
spontaneous transition that only becomes active for ``time > start``:

.. code-block:: python

   # Vaccinate 1% of susceptibles per day, starting on day 60
   SEIR.add_vaccination("S", "V", rate=0.01, start=60)
   SEIR.integrate(365, S=9_990, I=10, R=0, V=0)

Demographics: births and deaths
-------------------------------

Per-capita birth and death rates can be attached to one or all
compartments. Enabling either puts the model in *demographics mode*:

.. code-block:: python

   SEIR.add_birth_rate(0.0001)              # applied to every compartment
   SEIR.add_death_rate(0.0001, comps=["R"])  # applied only to "R"

Both ODE integration and stochastic simulation are aware of demographics:
in :meth:`~epidemik.EpiModel.integrate` births and deaths enter the
equations as ``rate * X`` terms, while in
:meth:`~epidemik.EpiModel.simulate` births are sampled from a binomial
distribution and deaths compete with the other transitions out of each
compartment in the same multinomial draw. By default the inflow into compartment ``X`` is
proportional to ``X`` itself; pass ``from_group=True`` to
:meth:`~epidemik.EpiModel.add_birth_rate` to make it proportional to the
total population of the group that contains ``X`` instead (see
`Multi-group and vector-borne models`_ below).

Multi-group and vector-borne models
-----------------------------------

By default every interaction is normalized by the total population
``N``. For models with several interacting populations -- hosts and
vectors, for instance -- the force of infection acting on the hosts should
instead be normalized by the number of *vectors*, and vice versa.
:meth:`~epidemik.EpiModel.add_groups` declares such a partition of the
compartments: after the call, the interaction ``X + A -> Y + A``
contributes ``rate * X * A / N_group(A)``, where ``N_group(A)`` is the
total population of the group containing the agent ``A``.

.. code-block:: python

   bh, bv, gh, gv, mh, mv = 0.3, 0.3, 0.1, 0.1, 0.01, 0.1

   vector_model = EpiModel()
   vector_model.add_interaction("Sh", "Ih", "Iv", bh)   # Sh + Iv -> Ih + Iv
   vector_model.add_interaction("Sv", "Iv", "Ih", bv)   # Sv + Ih -> Iv + Ih
   vector_model.add_spontaneous("Ih", "Rh", gh)

   vector_model.add_groups({
       "host":   ["Sh", "Ih", "Rh"],
       "vector": ["Sv", "Iv"],
   })

   # hosts die at rate mh, vectors die at rate gv, and mv * N_vector new
   # susceptible vectors are born per unit time
   vector_model.add_death_rate(mh, comps=["Sh", "Ih", "Rh"])
   vector_model.add_death_rate(gv, comps=["Sv", "Iv"])
   vector_model.add_birth_rate(mv, comps=["Sv"], from_group=True)

   vector_model.integrate(100, Sh=999, Ih=1, Rh=0, Sv=1_000, Iv=10)

This integrates the classical host-vector system

.. math::

   \begin{aligned}
   \dot S_h &= -b_h S_h I_v / N_v - m_h S_h \\
   \dot I_h &= b_h S_h I_v / N_v - g_h I_h - m_h I_h \\
   \dot R_h &= g_h I_h - m_h R_h \\
   \dot S_v &= m_v N_v - b_v S_v I_h / N_h - g_v S_v \\
   \dot I_v &= b_v S_v I_h / N_h - g_v I_v
   \end{aligned}

Every compartment of the model must belong to exactly one group;
:meth:`~epidemik.EpiModel.integrate` and
:meth:`~epidemik.EpiModel.simulate` raise :class:`ValueError` if a
compartment added after :meth:`~epidemik.EpiModel.add_groups` was left
unassigned. If a group is empty its force of infection is simply zero.
Groups cannot be combined with age structure, and
:meth:`~epidemik.EpiModel.R0` returns ``None`` for multi-group models.

Age structure
-------------

Age structure is layered on top of an existing model by calling
:meth:`~epidemik.EpiModel.add_age_structure` with a contact matrix and a
vector of per-age-group populations. Every compartment is replicated into
one compartment per age group (``S_a``, ``S_b``, ...), every interaction
rate is rescaled by the corresponding contact matrix entry, and
interactions are normalized by the population of the *agent's* age group.
Birth and death rates attached before the call are copied to every age
stratum.

.. code-block:: python

   import numpy as np

   contact = np.array([
       [10, 4, 1],
       [4, 12, 4],
       [1, 4, 8],
   ])
   population = np.array([1_000, 2_000, 1_500])

   age_SEIR = EpiModel()
   age_SEIR.add_interaction("S", "E", "I", 0.05)
   age_SEIR.add_spontaneous("E", "I", 0.2)
   age_SEIR.add_spontaneous("I", "R", 0.1)
   age_SEIR.add_age_structure(contact, population)

   age_SEIR.integrate(365, S=4_490, I=10, R=0)

Initial conditions given with the original compartment names are split
across age groups proportionally to ``population`` (exactly in
:meth:`~epidemik.EpiModel.integrate`, with a multinomial draw in
:meth:`~epidemik.EpiModel.simulate`). Expanded names such as ``S_a=500``
can also be passed directly. After the run, the per-age-group time series
are available in :attr:`~epidemik.EpiModel.values_ages_` and the totals
(across age groups) in :attr:`~epidemik.EpiModel.values_`. Individual age
strata can be retrieved with ``model["S_a"]`` or ``model[["S_a", "S_b"]]``.

Seasonality
-----------

Both :meth:`~epidemik.EpiModel.integrate` and
:meth:`~epidemik.EpiModel.simulate` accept a ``seasonality`` argument: a
array of multiplicative factors (typically one per day of the year) that
is applied to every
interaction rate. At time ``t`` the effective rate is multiplied by
``seasonality[t % len(seasonality)]``.

.. code-block:: python

   import numpy as np

   t = np.arange(365)
   season = 1 + 0.3 * np.cos(2 * np.pi * t / 365)
   SEIR.integrate(365, seasonality=season, S=9_990, I=10, R=0)

Accessing and visualising results
---------------------------------

After running either :meth:`~epidemik.EpiModel.integrate` or
:meth:`~epidemik.EpiModel.simulate`:

* :attr:`~epidemik.EpiModel.values_` is the full trajectory as a
  :class:`pandas.DataFrame` indexed by time.
* ``model.X`` returns the time series of compartment ``X`` (it is shorthand
  for ``model.values_["X"]``).
* ``model["X"]`` and ``model[["X", "Y"]]`` return the time series of one or
  several compartments; ``None`` is returned when the name is unknown.
* :meth:`~epidemik.EpiModel.plot` produces a quick line plot with
  reasonable defaults (using the colours defined in
  :data:`epidemik.utils.epi_colors`), optionally on a caller-supplied
  ``ax``.
* :meth:`~epidemik.EpiModel.draw_model` draws the compartmental structure
  itself.

Using ``model["X"]`` or :meth:`~epidemik.EpiModel.plot` before running the
model raises :class:`~epidemik.utils.NotInitialized`.
