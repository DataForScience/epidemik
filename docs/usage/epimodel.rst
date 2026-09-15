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
above (agents are created too, so an environmental reservoir that only
ever appears as an agent is a valid compartment). Several transitions may
share the same source and target (for example two interactions driven by
different agents); their contributions are added together.

``print(model)`` gives a YAML summary of the parameters, every transition,
any birth/death rates and groups, and :math:`R_0`;
:meth:`~epidemik.EpiModel.draw_model` draws the compartment graph.

Named parameters
----------------

Every rate is stored as a named parameter in
:attr:`~epidemik.EpiModel.params`, an instance of
:class:`~epidemik.utils.Parameters`. A rate can be given in three ways:

* as a plain number, in which case an anonymous parameter ``rate1``,
  ``rate2``, ... is created;
* as a single keyword argument, ``beta=0.3``, which names the parameter;
* as a string: the name of an existing parameter, or an *expression* in
  terms of other parameters such as ``"beta/2"`` (the functions ``exp``,
  ``log``, ``sqrt``, ``sin``, ``cos``, ``abs``, ``min``, ``max`` and the
  constants ``pi`` and ``e`` are available).

.. code-block:: python

   SIR = EpiModel()
   SIR.define_parameters(beta=0.3)
   SIR.add_interaction("S", "I", "I", "beta")
   SIR.add_spontaneous("I", "R", mu="beta/3")   # mu is tied to beta

   SIR.params["mu"]        # -> 0.1 (evaluated)
   SIR.params.raw("mu")    # -> "beta/3"
   SIR.define_parameters(beta=0.6)
   SIR.R0()                # -> 3.0, mu followed beta

Expressions are evaluated when the model is run (or when :math:`R_0` is
computed); an expression that refers to an unknown name raises
:class:`KeyError` at that point, and one that cannot be evaluated raises
:class:`ValueError`.

Computing :math:`R_0`
---------------------

The basic reproduction number is computed using the next-generation matrix
approach, so it works for any model topology -- SIR, SEIR, SIS, SIRS,
models with several infectious compartments, and so on:

.. code-block:: python

   SEIR.R0()

Vaccination transitions and demographics are ignored in the computation,
and every susceptible class is assumed to hold the whole population at the
disease-free equilibrium. Susceptible classes are the sources of
interactions that are not infected themselves and are only fed by infected
or post-infection compartments -- so ``S`` in an SIRS model qualifies,
while a class ``V`` fed by vaccinating ``S`` does not.

If the model has no infectious dynamics (or :math:`R_0` cannot be derived
for some other reason), :meth:`~epidemik.EpiModel.R0` returns ``None``.
Multi-group (vector-borne) models and within-host models (unnormalized
interactions or viral generation) are not supported and also return
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

All random numbers used by :meth:`~epidemik.EpiModel.simulate` come from
the model's :attr:`~epidemik.EpiModel.rng`, a
:class:`numpy.random.Generator` created from the ``seed`` argument of the
constructor (or passed in directly as ``rng``), so simulations are
reproducible:

.. code-block:: python

   SEIR = EpiModel(seed=42)

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

Birth and death rates can be attached to one or all compartments.
Enabling either puts the model in *demographics mode*:

.. code-block:: python

   SEIR.add_birth_rate(0.0001, comps=["S"])   # newborns are susceptible
   SEIR.add_death_rate(0.0001)               # applied to every compartment

By default births into a compartment are proportional to the **total**
population (``rate * N``, the classical setup where every newborn is
susceptible) and deaths are proportional to the compartment itself
(``rate * X``). Both methods accept further modes:

* ``fixed=True`` -- a constant number of individuals per unit time
  (e.g. a constant supply of target cells);
* ``global_rate=True/False`` -- proportional to the total population or
  to the compartment itself;
* ``from_group=True`` (births only) -- proportional to the population of
  the group containing the compartment (see `Multi-group and
  vector-borne models`_ below).

Both ODE integration and stochastic simulation are aware of demographics:
in :meth:`~epidemik.EpiModel.integrate` births and deaths enter the
equations as flow terms, while in :meth:`~epidemik.EpiModel.simulate`
births are sampled from a binomial (or Poisson) distribution and deaths
compete with the other transitions out of each compartment in the same
multinomial draw, so compartments never become negative.

For backwards compatibility the list of compartments may also be given
first together with a named rate: ``add_birth_rate("S", b=0.01)``.

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

Within-host models
------------------

Two additional building blocks make it possible to write models of
infection *within* a host, such as the target-cell-limited model of viral
dynamics:

* ``norm=False`` on :meth:`~epidemik.EpiModel.add_interaction` turns the
  term into the plain mass-action product ``rate * source * agent``,
  without dividing by the population;
* :meth:`~epidemik.EpiModel.add_viral_generation` encodes the asymmetric
  reaction ``source => target``: infected cells die at one rate while, per
  infected cell, virions are produced at another.

.. code-block:: python

   TCL = EpiModel()
   TCL.define_parameters(s=1e5, d=0.01, beta=2e-7, delta=0.5, p=100, c=2)
   TCL.add_birth_rate("s", comps=["T"], fixed=True)   # constant supply of cells
   TCL.add_death_rate("d", comps=["T"])
   TCL.add_interaction("T", "I", "V", "beta", norm=False)  # T + V -> I + V
   TCL.add_viral_generation("I", "V", "delta", "p")        # I dies, V is produced
   TCL.add_death_rate("c", comps=["V"])                    # viral clearance
   TCL.integrate(60, T=1e7, I=0, V=10)

This integrates :math:`\dot T = s - dT - \beta T V`,
:math:`\dot I = \beta T V - \delta I` and :math:`\dot V = pI - cV`.
In :meth:`~epidemik.EpiModel.simulate`, cell death competes with the other
transitions and virion production is drawn from a Poisson distribution.
:meth:`~epidemik.EpiModel.R0` returns ``None`` for such models, since the
within-host reproduction number depends on the initial number of target
cells.

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
can also be passed directly. The rescaled interaction rates become new
parameters named ``<rate>_<ij>`` (e.g. ``beta_ab``) defined as expressions
of the original ones, so redefining ``beta`` afterwards still works. After the run, the per-age-group time series
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
:meth:`~epidemik.EpiModel.reset` discards the results of a previous run.

Saving, loading and sharing models
----------------------------------

A model definition can be written to a YAML file and read back:

.. code-block:: python

   SIR.name = "SIR"
   SIR.save_model("SIR.yaml")
   SIR2 = EpiModel.load_model("SIR.yaml", seed=1)

The file is exactly the output of ``print(SIR)``:

.. code-block:: yaml

   # Epidemic Model with 3 compartments and 2 transitions:
   Name: SIR

   Compartments: [S, I, R]

   Parameters:
     beta: 0.3
     mu: "beta/3"

   Transitions:
     - S + I = I beta
     - I -> R mu

   # R0=3.00

Transitions are written as ``S + I = I beta`` (interaction, with a
trailing ``unnormalized`` when ``norm=False``), ``I -> R mu``
(spontaneous), ``S -> V nu starting at 60 days`` (vaccination) and
``I => V delta p`` (viral generation). Rates may be parameter names or
numeric literals. Optional ``Demographics`` (a list of
``{compartment: S, birth: b}`` / ``{compartment: I, death: d}`` mappings
with ``fixed``, ``global`` and ``from_group`` flags) and ``Groups``
sections round-trip the corresponding features; ``Compartments`` fixes
the order of the compartments and may be omitted.

A small library of ready-made models (SI, SIR, SEIR, SEIRS, SEIIR,
SEIIRD) is hosted in the ``epidemik`` GitHub repository:

.. code-block:: python

   EpiModel.list_models()                     # -> ["SEIIR.yaml", "SEIIRD.yaml", ...]
   SEIR = EpiModel.download_model("SEIR.yaml")

Downloaded files are cached in a per-user directory (see
:func:`epidemik.utils.get_cache_directory`; override it with the
``EPIDEMIK_CACHE_DIR`` environment variable). Both functions accept a
``repo`` argument that can be another GitHub repository or a local
directory of model files.
