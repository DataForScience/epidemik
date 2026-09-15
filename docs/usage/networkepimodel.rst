NetworkEpiModel: epidemics on a contact network
===============================================

:class:`~epidemik.NetworkEpiModel` lifts an :class:`~epidemik.EpiModel`
onto an arbitrary contact network supplied as a :class:`networkx.Graph`.
Every individual is represented by a single network node and can only
infect its direct neighbours. Only stochastic simulation is supported --
the corresponding ODE system has no closed-form analogue on a generic
graph.

Building a network model
------------------------

The constructor takes the contact network as its first argument; the
disease transitions are then added in exactly the same way as for an
:class:`~epidemik.EpiModel`. Node labels can be of any hashable type
(integers, strings, tuples, ...):

.. code-block:: python

   import networkx as nx
   from epidemik import NetworkEpiModel

   G = nx.erdos_renyi_graph(n=1_000, p=0.01, seed=42)

   SIR = NetworkEpiModel(G, seed=42)
   SIR.add_interaction("S", "I", "I", beta=0.05)
   SIR.add_spontaneous("I", "R", mu=0.1)

Rates are named parameters exactly as in :class:`~epidemik.EpiModel`
(numbers, keywords or expressions). On a network, a rate is the
probability that a single contact transmits the disease during one time
step (for interactions) or that an individual undergoes the transition
during one time step (for spontaneous transitions). If you want to keep
:math:`R_0` comparable to a mean-field model defined with the same rates,
you can ask :meth:`~epidemik.NetworkEpiModel.add_interaction` to rescale
the rate by the average degree of the network; the rescaled value is
stored as a new parameter ``<name>_rescaled``:

.. code-block:: python

   SIR.add_interaction("S", "I", "I", beta=0.05, rescale=True)

Vaccination campaigns are supported through
:meth:`~epidemik.NetworkEpiModel.add_vaccination` and behave like
time-gated spontaneous transitions.

Running a simulation
--------------------

Instead of supplying initial compartment sizes, simulations are seeded by
passing a ``{node_label: compartment_name}`` dictionary. All non-seeded
nodes start in the susceptible compartment, which is detected
automatically (the compartment without incoming transitions) or can be
given explicitly with ``susceptible=``:

.. code-block:: python

   seeds = {0: "I", 1: "I", 2: "I"}
   SIR.simulate(timesteps=60, seeds=seeds)

After the run:

* ``SIR.population_`` is a :class:`pandas.DataFrame` of shape
  ``(timesteps, n_nodes)`` containing the compartment of every node at
  every time step. Columns are labelled by node, rows by time.
* :attr:`~epidemik.EpiModel.values_` contains the aggregate compartment
  counts across the network over time, with one column per compartment
  of the model (compartments that are never reached contain zeros).
* :meth:`~epidemik.EpiModel.plot` produces the usual aggregate plot.

As for :class:`~epidemik.EpiModel`, the first row of both data frames
(at ``t_min``, which defaults to 1) holds the initial conditions, and all
random numbers come from the model's ``rng`` (see the ``seed`` argument
of the constructor). Seeding a node that is not part of the network, or
with a state that is not a compartment, raises :class:`ValueError`. If
the transition probabilities out of a compartment (or from one contact)
add up to more than 1 they are renormalized and a warning is emitted.
When several interactions share the same source and agent, one of the
targets is chosen at random with probabilities proportional to the rates.

R\ :sub:`0` on a network
------------------------

The :meth:`~epidemik.NetworkEpiModel.R0` method returns the network-corrected
basic reproduction number. The mean-field :math:`R_0` of the underlying
disease is multiplied by the average degree of the contact network (stored
in ``kavg_``) to account for the limited number of contacts each
individual has, and rounded to two decimals.

.. code-block:: python

   SIR.R0()

The method returns ``None`` whenever the mean-field :math:`R_0` cannot be
computed (see :meth:`EpiModel.R0 <epidemik.EpiModel.R0>`).

What is *not* supported
-----------------------

Some features available to :class:`~epidemik.EpiModel` are intentionally
omitted on networks:

* :meth:`~epidemik.NetworkEpiModel.integrate` raises
  :class:`NotImplementedError` -- explore the dynamics with
  :meth:`~epidemik.NetworkEpiModel.simulate` instead.
* Births, deaths, viral generation, groups and age structure do not make
  sense when every node is one individual of a fixed population;
  the corresponding ``add_*`` methods raise :class:`NotImplementedError`.
* Seasonality is not taken into account by
  :meth:`~epidemik.NetworkEpiModel.simulate`.
