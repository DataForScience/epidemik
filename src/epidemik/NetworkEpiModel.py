### −∗− mode : python ; −∗−
# @file NetworkEpiModel.py
# @author Bruno Goncalves
######################################################

from typing import Dict, Union
import warnings

import numpy as np
import pandas as pd

from .EpiModel import EpiModel


class NetworkEpiModel(EpiModel):
    """Compartmental epidemic model that lives on a contact network.

    Unlike :class:`~epidemik.EpiModel`, which assumes a homogeneously mixed
    population, :class:`NetworkEpiModel` simulates individual nodes that can
    only infect their neighbours along the edges of the provided
    :class:`networkx.Graph`. Only stochastic simulation is supported --
    numerical integration of the corresponding ODEs is not meaningful in
    this setting.

    Attributes
    ----------
    network : networkx.Graph
        The contact network.
    kavg_ : float
        Average degree of the contact network.
    population_ : pandas.DataFrame
        State of every node at every time step (one column per node),
        available after :meth:`simulate`.
    values_ : pandas.DataFrame
        Number of nodes in each compartment at every time step, available
        after :meth:`simulate`.
    """

    def __init__(self, network, compartments=None, seed=None, rng=None, name=None):
        """Initialize a network-based epidemic model.

        Parameters
        ----------
        network : networkx.Graph
            Contact network over which the epidemic will spread. Nodes
            represent individuals and edges represent contacts that can
            transmit the disease. Node labels can be of any hashable type.
        compartments : list of str, optional
            Optional list of compartment names. Just like in
            :class:`~epidemik.EpiModel`, compartments are also added
            implicitly the first time they appear in a transition.
        seed : int, optional
            Seed for the random number generator used by :meth:`simulate`.
        rng : numpy.random.Generator, optional
            Random number generator to use instead of a new seeded one.
        name : str, optional
            Human-readable name of the model.
        """
        super(NetworkEpiModel, self).__init__(compartments, seed=seed, rng=rng, name=name)
        self.network = network

        n_nodes = network.number_of_nodes()
        self.kavg_ = 2*network.number_of_edges()/n_nodes if n_nodes > 0 else 0.0

        self.spontaneous = {}
        self.interactions = {}

    def integrate(self, *args, **kwargs):
        """Not supported for network models.

        Raises
        ------
        NotImplementedError
            Always. Network models can only be explored via
            :meth:`simulate`.
        """
        raise NotImplementedError("Network Models don't support numerical integration")

    def add_interaction(self, source: str, target: str, agent: str, rate=None, rescale: bool = False, **rates) -> None:
        """Add a contagious interaction between two compartments.

        Parameters
        ----------
        source : str
            Name of the source compartment.
        target : str
            Name of the target compartment.
        agent : str
            Name of the agent compartment driving the interaction.
        rate : float or str, optional
            Per-contact transmission probability per time step (number,
            parameter name or expression).
        rescale : bool, default=False
            If True, divide the rate by the average degree of the underlying
            network. This is useful when porting a mean-field rate to a
            network setting while keeping :math:`R_0` comparable. The
            rescaled value is stored as a new parameter ``<name>_rescaled``
            defined in terms of the original one.
        **rates
            Alternatively, the rate as a single named parameter, e.g.
            ``beta=0.05``.
        """
        key = self._register_rate(rate, rates, "add_interaction")

        if rescale:
            rescaled = key + "_rescaled"
            self.params.define_parameters(**{rescaled: "(%s)/%r" % (key, float(self.kavg_))})
            key = rescaled

        super(NetworkEpiModel, self).add_interaction(source, target, agent, key)

        self.interactions.setdefault(source, {}).setdefault(agent, []).append({'target': target, 'rate': key})

    def add_spontaneous(self, source: str, target: str, rate=None, **rates) -> None:
        """Add a spontaneous transition between two compartments.

        Parameters
        ----------
        source : str
            Name of the source compartment.
        target : str
            Name of the target compartment.
        rate : float or str, optional
            Per-individual transition probability per time step (number,
            parameter name or expression).
        **rates
            Alternatively, the rate as a single named parameter, e.g.
            ``mu=0.1``.
        """
        key = self._register_rate(rate, rates, "add_spontaneous")
        super(NetworkEpiModel, self).add_spontaneous(source, target, key)

        self.spontaneous.setdefault(source, {})[target] = {'rate': key, 'start': None}

    def add_vaccination(self, source: str, target: str, rate=None, start: Union[int, None] = None, **rates) -> None:
        """Add a time-gated vaccination transition between two compartments.

        Behaves like :meth:`add_spontaneous` but is only active for
        ``time > start``.

        Parameters
        ----------
        source : str
            Name of the source compartment.
        target : str
            Name of the target compartment.
        rate : float or str, optional
            Per-individual vaccination probability per time step.
        start : int
            Time step at which the campaign begins.
        **rates
            Alternatively, the rate as a single named parameter.
        """
        if start is None:
            raise TypeError("add_vaccination requires the 'start' time of the campaign.")

        key = self._register_rate(rate, rates, "add_vaccination")
        super(NetworkEpiModel, self).add_vaccination(source, target, key, start=start)

        self.spontaneous.setdefault(source, {})[target] = {'rate': key, 'start': start}

    def _unsupported(self, what: str):
        raise NotImplementedError(
            "%s is not supported by NetworkEpiModel: every node is a single individual "
            "of a fixed population." % what
        )

    def add_birth_rate(self, *args, **kwargs):
        """Not supported for network models (raises :class:`NotImplementedError`)."""
        self._unsupported("add_birth_rate")

    def add_death_rate(self, *args, **kwargs):
        """Not supported for network models (raises :class:`NotImplementedError`)."""
        self._unsupported("add_death_rate")

    def add_viral_generation(self, *args, **kwargs):
        """Not supported for network models (raises :class:`NotImplementedError`)."""
        self._unsupported("add_viral_generation")

    def add_age_structure(self, *args, **kwargs):
        """Not supported for network models (raises :class:`NotImplementedError`)."""
        self._unsupported("add_age_structure")

    def add_groups(self, *args, **kwargs):
        """Not supported for network models (raises :class:`NotImplementedError`)."""
        self._unsupported("add_groups")

    def simulate(self, timesteps: int, seeds: Dict, t_min: int = 1, susceptible: Union[str, None] = None, **kwargs) -> None:
        """Run a discrete-time stochastic simulation on the network.

        At every time step each active node attempts to infect its
        neighbours and may also undergo spontaneous (or, once started,
        vaccination) transitions. The per-node state trajectory is stored
        in :attr:`population_` and the aggregate compartment counts over
        time are stored in :attr:`values_`. In both data frames the first
        row (``t_min``) holds the initial conditions. All random numbers
        come from :attr:`rng`. Seasonality is not supported on networks.

        Parameters
        ----------
        timesteps : int
            Number of time steps (rows) to produce, including the initial
            condition.
        seeds : dict
            Mapping ``{node_label: compartment_name}`` that sets the initial
            state of a subset of network nodes. All other nodes are placed
            in the susceptible compartment.
        t_min : int, default=1
            Time label of the first row of :attr:`population_` and
            :attr:`values_`.
        susceptible : str, optional
            Name of the susceptible compartment. If omitted, the compartment
            without incoming transitions is used (as for
            :class:`~epidemik.EpiModel`).
        **kwargs
            Reserved for future use; currently ignored.

        Raises
        ------
        ValueError
            If a seed node is not part of the network, if a seed state is
            not a compartment of the model, or if the susceptible
            compartment cannot be determined.

        Warns
        -----
        UserWarning
            If the spontaneous transition probabilities out of a
            compartment sum to more than 1 (they are renormalized).
        """
        if timesteps < 1:
            raise ValueError("timesteps must be a positive integer.")

        pos = {comp: i for i, comp in enumerate(self.transitions.nodes())}
        comps = list(pos)
        rates = self.params.resolve()

        nodes = list(self.network.nodes())
        node_index = {node: i for i, node in enumerate(nodes)}
        N = len(nodes)

        if susceptible is None:
            candidates = self._get_susceptible()

            if len(candidates) == 0:
                raise ValueError(
                    "Unable to determine the susceptible compartment; "
                    "pass it explicitly with susceptible=..."
                )

            susceptible = sorted(candidates)[0]
        elif susceptible not in pos:
            raise ValueError("'%s' is not a compartment of the model." % susceptible)

        population = np.empty((timesteps, N), dtype=object)
        population[0, :] = susceptible

        active_states = self._get_active()

        for node, state in seeds.items():
            if node not in node_index:
                raise ValueError("Seed node %r is not part of the network." % (node,))

            if state not in pos:
                raise ValueError("Seed state '%s' is not a compartment of the model." % state)

            population[0, node_index[node]] = state

        # Nodes that can do something on their own: infectious ones and, if
        # the susceptible state has spontaneous transitions (e.g.
        # vaccination), every susceptible node too.
        active_nodes = set(i for i in range(N) if population[0, i] in active_states)

        # {agent: {source: ([targets], [cumulative probabilities])}}
        infections = {}

        for agent, by_source in self._get_infections().items():
            infections[agent] = {}

            for source, options in by_source.items():
                targets = [option['target'] for option in options]
                probs = np.array([rates[option['rate']] for option in options], dtype='float')

                if probs.sum() > 1:
                    warnings.warn(
                        "The transmission probabilities from '%s' to '%s' sum to more than 1; "
                        "they were renormalized." % (agent, source), stacklevel=2
                    )
                    probs /= probs.sum()

                infections[agent][source] = (targets, np.cumsum(probs))

        warned = False

        for t in range(t_min+1, t_min+timesteps):
            row = t - t_min
            population[row] = population[row-1]

            if len(active_nodes) == 0:
                continue

            current_active = list(active_nodes)
            self.rng.shuffle(current_active)

            for node_i in current_active:
                state_i = population[row-1, node_i]

                if state_i in infections:
                    # contact each neighbour to see if we infect them
                    NN = [node_index[n] for n in self.network.neighbors(nodes[node_i])]
                    self.rng.shuffle(NN)

                    for node_j in NN:
                        state_j = population[row-1, node_j]

                        if state_j in infections[state_i]:
                            targets, cumulative = infections[state_i][state_j]
                            draw = self.rng.random()
                            choice = int(np.searchsorted(cumulative, draw, side='right'))

                            if choice < len(targets):
                                new_state = targets[choice]
                                population[row, node_j] = new_state

                                if new_state in active_states:
                                    active_nodes.add(node_j)
                                else:
                                    active_nodes.discard(node_j)

                if state_i in self.spontaneous:
                    prob = np.zeros(len(pos))

                    for target, option in self.spontaneous[state_i].items():
                        if option['start'] is not None and option['start'] >= t:
                            continue

                        prob[pos[target]] += rates[option['rate']]

                    total = prob.sum()

                    if total > 1:
                        if not warned:
                            warnings.warn(
                                "The spontaneous transition probabilities out of '%s' sum "
                                "to more than 1; they were renormalized." % state_i,
                                stacklevel=2
                            )
                            warned = True

                        prob /= total
                        total = 1.0

                    prob[pos[state_i]] = max(0.0, 1.0-total)

                    new_state = comps[np.argmax(self.rng.multinomial(1, prob))]

                    if new_state != state_i:
                        population[row, node_i] = new_state

                        if new_state in active_states:
                            active_nodes.add(node_i)
                        else:
                            active_nodes.discard(node_i)

        time = np.arange(t_min, t_min+timesteps, 1, dtype='int')

        self.population_ = pd.DataFrame(population, index=time, columns=nodes)

        counts = np.zeros((timesteps, len(comps)), dtype='int')

        for k, comp in enumerate(comps):
            counts[:, k] = (population == comp).sum(axis=1)

        self.values_ = pd.DataFrame(counts, index=time, columns=comps)

    def R0(self) -> Union[float, None]:
        """Compute the network-corrected basic reproduction number.

        The underlying mean-field :math:`R_0` (computed by
        :meth:`EpiModel.R0`) is multiplied by the average degree of the
        contact network, since each infected node only has :math:`\\langle k
        \\rangle` contacts on average.

        Returns
        -------
        float or None
            Network :math:`R_0` rounded to two decimals, or ``None`` if the
            mean-field :math:`R_0` cannot be computed.
        """
        r0 = super(NetworkEpiModel, self).R0()

        if r0 is None:
            return None

        return float(np.round(r0*self.kavg_, 2))
