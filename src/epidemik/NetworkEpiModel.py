### −∗− mode : python ; −∗−
# @file NetworkEpiModel.py
# @author Bruno Goncalves
######################################################

from typing import Union
import networkx as nx
import numpy as np
from numpy import linalg
from numpy import random
import pandas as pd
import matplotlib.pyplot as plt
from .EpiModel import EpiModel
from collections import Counter
from .utils import *

class NetworkEpiModel(EpiModel):
    """Compartmental epidemic model that lives on a contact network.

    Unlike :class:`~epidemik.EpiModel`, which assumes a homogeneously mixed
    population, :class:`NetworkEpiModel` simulates individual nodes that can
    only infect their neighbours along the edges of the provided
    :class:`networkx.Graph`. Only stochastic simulation is supported --
    numerical integration of the corresponding ODEs is not meaningful in
    this setting.
    """

    def __init__(self, network, compartments=None):
        """Initialize a network-based epidemic model.

        Parameters
        ----------
        network : networkx.Graph
            Contact network over which the epidemic will spread. Nodes
            represent individuals and edges represent contacts that can
            transmit the disease.
        compartments : list of str, optional
            Optional list of compartment names. Just like in
            :class:`~epidemik.EpiModel`, compartments are also added
            implicitly the first time they appear in a transition.
        """
        super(NetworkEpiModel, self).__init__(compartments)
        self.network = network
        self.kavg_ = 2*network.number_of_edges()/network.number_of_nodes() 
        self.spontaneous = {}
        self.interactions = {}

    def integrate(self, timesteps, **kwargs):
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
            Per-contact transmission rate.
        rescale : bool, default=False
            If True, divide ``rate`` by the average degree of the underlying
            network. This is useful when porting a mean-field rate to a
            network setting while keeping :math:`R_0` comparable.
        """
        if rescale:
            rate /= self.kavg_

        super(NetworkEpiModel, self).add_interaction(source, target, agent=agent, rate=rate)

        if source not in self.interactions:
            self.interactions[source] = {}

        if target not in self.interactions[source]:
            self.interactions[source] = {}

        self.interactions[source][agent] = {'target': target, 'rate': rate}

    def add_spontaneous(self, source: str, target: str, rate: float) -> None:
        """Add a spontaneous transition between two compartments.

        Parameters
        ----------
        source : str
            Name of the source compartment.
        target : str
            Name of the target compartment.
        rate : float
            Per-individual transition rate.
        """
        super(NetworkEpiModel, self).add_spontaneous(source, target, rate=rate)
        if source not in self.spontaneous:
            self.spontaneous[source] = {}

        if target not in self.spontaneous[source]:
            self.spontaneous[source] = {}

        self.spontaneous[source][target] = rate


    def simulate(self, timesteps: int, seeds, **kwargs) -> None:
        """Run a discrete-time stochastic simulation on the network.

        At every time step each infected node attempts to infect its
        neighbours and may also undergo spontaneous transitions. The
        per-node state trajectory is stored in :attr:`population_` and the
        aggregate compartment counts over time are stored in
        :attr:`values_`.

        Parameters
        ----------
        timesteps : int
            Number of time steps to simulate.
        seeds : dict
            Mapping ``{node_id: compartment_name}`` that sets the initial
            state of a subset of network nodes. All other nodes are placed
            in the susceptible compartment.
        **kwargs
            Reserved for future use; currently ignored.
        """
        pos = {comp: i for i, comp in enumerate(self.transitions.nodes())}
        N = self.network.number_of_nodes()

        population = np.zeros((timesteps, N), dtype='str')

        comps = list(self.transitions.nodes)
        time = np.arange(1, timesteps, 1, dtype='int')

        susceptible = self._get_susceptible().pop()

        active_nodes = set()
        current_active = set()
        active_states = self._get_active()

        for node in range(N):
            if node in seeds:
                population[0, node] = seeds[node]
                active_nodes.add(node)
            else:
                population[0, node] = susceptible

        infections = self._get_infections()

        for t in time:
            population[t] = np.copy(population[t-1])

            if len(active_nodes) == 0:
                continue

            current_active = list(active_nodes)
            np.random.shuffle(current_active)

            for node_i in current_active:
                state_i = population[t-1, node_i]

                if state_i in infections:
                    # contact each neighbour to see if we infect them
                    NN = list(self.network.neighbors(node_i))
                    np.random.shuffle(NN)

                    for node_j in NN:
                        state_j = population[t-1, node_j]

                        if state_j in infections[state_i]:
                            prob = np.random.random()

                            if prob < infections[state_i][state_j]['rate']:
                                new_state = infections[state_i][state_j]['target']
                                population[t, node_j] = new_state

                                if new_state in active_states:
                                    active_nodes.add(node_j)
 
                if state_i in self.spontaneous:
                    n_trans = len(self.spontaneous[state_i])

                    prob = np.zeros(len(pos))

                    for target in self.spontaneous[state_i]:
                        prob[pos[target]] = self.spontaneous[state_i][target]

                    prob[pos[state_i]] = 1-np.sum(prob)

                    new_state = comps[np.argmax(random.multinomial(1, prob))]

                    if new_state != state_i:
                        population[t, node_i] = new_state
                        
                        active_nodes.add(node_i)

                        if new_state not in active_states:
                            active_nodes.remove(node_i)
                        
                        continue


        self.population_ = pd.DataFrame(population)
        self.values_ = pd.DataFrame.from_records(self.population_.apply(lambda x: Counter(x), axis=1)).fillna(0).astype('int')

    def R0(self) -> Union[float, None]:
        """Compute the network-corrected basic reproduction number.

        The underlying mean-field :math:`R_0` (computed by the parent class)
        is multiplied by the average degree of the contact network to
        account for the limited number of contacts each individual has.

        Returns
        -------
        float or None
            Network :math:`R_0`, or ``None`` if no recovered compartment
            (``"R"``) is defined.
        """
        if 'R' not in set(self.transitions.nodes):
            return None
        return np.round(super(NetworkEpiModel, self).R0()*self.kavg_, 2)
