### −∗− mode : python ; −∗−
# @file NetworkEpiModel.py
# @author Bruno Goncalves
######################################################

from typing import Dict, Union
import warnings

import numpy as np
from numpy import random
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

    def __init__(self, network, compartments=None):
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
        """
        super(NetworkEpiModel, self).__init__(compartments)
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

    def add_interaction(self, source: str, target: str, agent: str, rate: float, rescale: bool = False) -> None:
        """Add a contagious interaction between two compartments.

        Parameters
        ----------
        source : str
            Name of the source compartment.
        target : str
            Name of the target compartment.
        agent : str
            Name of the agent compartment driving the interaction.
        rate : float
            Per-contact transmission probability per time step.
        rescale : bool, default=False
            If True, divide ``rate`` by the average degree of the underlying
            network. This is useful when porting a mean-field rate to a
            network setting while keeping :math:`R_0` comparable.
        """
        if rescale:
            rate /= self.kavg_

        super(NetworkEpiModel, self).add_interaction(source, target, agent=agent, rate=rate)

        self.interactions.setdefault(source, {})[agent] = {'target': target, 'rate': rate}

    def add_spontaneous(self, source: str, target: str, rate: float) -> None:
        """Add a spontaneous transition between two compartments.

        Parameters
        ----------
        source : str
            Name of the source compartment.
        target : str
            Name of the target compartment.
        rate : float
            Per-individual transition probability per time step.
        """
        super(NetworkEpiModel, self).add_spontaneous(source, target, rate=rate)

        self.spontaneous.setdefault(source, {})[target] = rate

    def simulate(self, timesteps: int, seeds: Dict, t_min: int = 1, susceptible: Union[str, None] = None, **kwargs) -> None:
        """Run a discrete-time stochastic simulation on the network.

        At every time step each active node attempts to infect its
        neighbours and may also undergo spontaneous transitions. The
        per-node state trajectory is stored in :attr:`population_` and the
        aggregate compartment counts over time are stored in
        :attr:`values_`. In both data frames the first row (``t_min``)
        holds the initial conditions.

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
        active_nodes = set()

        for node, state in seeds.items():
            if node not in node_index:
                raise ValueError("Seed node %r is not part of the network." % (node,))

            if state not in pos:
                raise ValueError("Seed state '%s' is not a compartment of the model." % state)

            i = node_index[node]
            population[0, i] = state

            if state in active_states:
                active_nodes.add(i)

        infections = self._get_infections()
        warned = False

        for t in range(1, timesteps):
            population[t] = population[t-1]

            if len(active_nodes) == 0:
                continue

            current_active = list(active_nodes)
            np.random.shuffle(current_active)

            for node_i in current_active:
                state_i = population[t-1, node_i]

                if state_i in infections:
                    # contact each neighbour to see if we infect them
                    NN = [node_index[n] for n in self.network.neighbors(nodes[node_i])]
                    np.random.shuffle(NN)

                    for node_j in NN:
                        state_j = population[t-1, node_j]

                        if state_j in infections[state_i]:
                            if np.random.random() < infections[state_i][state_j]['rate']:
                                new_state = infections[state_i][state_j]['target']
                                population[t, node_j] = new_state

                                if new_state in active_states:
                                    active_nodes.add(node_j)
                                else:
                                    active_nodes.discard(node_j)

                if state_i in self.spontaneous:
                    prob = np.zeros(len(pos))

                    for target, rate in self.spontaneous[state_i].items():
                        prob[pos[target]] += rate

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

                    new_state = comps[np.argmax(random.multinomial(1, prob))]

                    if new_state != state_i:
                        population[t, node_i] = new_state

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
