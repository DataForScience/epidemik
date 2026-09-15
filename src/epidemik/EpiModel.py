### −∗− mode : python ; −∗−
# @file EpiModel.py
# @author Bruno Goncalves
######################################################

from typing import Dict, List, Set, Union
import warnings
import string

import networkx as nx
import numpy as np
from numpy import linalg
from numpy import random
import scipy.integrate
import pandas as pd
import matplotlib.pyplot as plt

from .utils import *

class EpiModel(object):
    """Compartmental epidemic model on a homogeneously mixed population.

    The model is internally represented as a directed multigraph
    (:class:`networkx.MultiDiGraph`) whose nodes are compartments and whose
    edges are transitions. Each transition can be either *spontaneous*
    (``X -> Y`` at some rate) or an *interaction* mediated by an agent
    compartment (``X + A -> Y + A`` at some rate). Once defined, the model
    dynamics can be explored deterministically with :meth:`integrate` (which
    numerically solves the corresponding ODEs) or stochastically with
    :meth:`simulate`.
    """

    def __init__(self, compartments=None):
        """Initialize an empty :class:`EpiModel` instance.

        Parameters
        ----------
        compartments : list of str, optional
            Optional list of compartment names to seed the model with. New
            compartments are also added automatically the first time they
            appear in a call to :meth:`add_interaction`,
            :meth:`add_spontaneous` or :meth:`add_vaccination`.
        """
        self.transitions = nx.MultiDiGraph()
        self.seasonality = None
        self.population = None
        self.orig_comps = None
        self.demographics = False
        self.groups = None
        self.comp_to_group = None

        if compartments is not None:
            self.transitions.add_nodes_from([comp for comp in compartments])
    
    def add_interaction(self, source: str, target: str, agent: str, rate: float) -> None:
        """Add an interaction transition between two compartments.

        Encodes a reaction of the form ``source + agent -> target + agent``
        that occurs with the given per-capita rate.

        Parameters
        ----------
        source : str
            Name of the source compartment.
        target : str
            Name of the target compartment.
        agent : str
            Name of the agent compartment driving the interaction.
        rate : float
            Per-capita rate at which the interaction takes place.
        """
        self.transitions.add_edge(source, target, agent=agent, rate=rate)

    def add_spontaneous(self, source: str, target: str, rate: float) -> None:
        """Add a spontaneous transition between two compartments.

        Encodes a reaction of the form ``source -> target`` that occurs at
        the given rate, independently of any other compartment.

        Parameters
        ----------
        source : str
            Name of the source compartment.
        target : str
            Name of the target compartment.
        rate : float
            Per-capita rate at which the transition takes place.
        """
        self.transitions.add_edge(source, target, rate=rate)

    def add_birth_rate(self, rate: float, comps: Union[List, None] = None, from_group: bool = False) -> None:
        """Attach a birth rate to one or more compartments.

        Enabling births (or deaths) switches the model into demographics
        mode: per-capita inflows are added to the chosen compartments during
        both :meth:`integrate` and :meth:`simulate`.

        Parameters
        ----------
        rate : float
            Birth rate. Interpreted per capita by default, or per
            *group* total when ``from_group=True`` -- see below.
        comps : list, optional
            Compartments to which to assign this birth rate. If ``None``
            (the default) the rate is applied to every compartment.
        from_group : bool, default=False
            If ``False`` (default), each listed compartment ``X`` receives
            an inflow of ``rate * X`` per unit time (the classical
            per-capita birth term). If ``True``, each listed compartment
            ``X`` receives an inflow of ``rate * N_{group(X)}``, where
            ``N_{group(X)}`` is the total population of the group that
            contains ``X``. This is the canonical setup for vector-borne
            models, where all newborn vectors are susceptible and the
            total vector birth flow is ``mv * Nv`` channelled into
            ``Sv``. Requires :meth:`add_groups` to have been called.

        Raises
        ------
        ValueError
            If ``from_group=True`` but :meth:`add_groups` has not been
            called yet.
        """
        if from_group and self.groups is None:
            raise ValueError(
                "add_birth_rate(..., from_group=True) requires add_groups "
                "to have been called first."
            )

        self.demographics = True

        if comps is None:
            comps = self.transitions.nodes()

        for comp in comps:
            self.transitions.nodes[comp]['birth'] = rate
            self.transitions.nodes[comp]['birth_from_group'] = from_group

    def add_death_rate(self, rate: float, comps: Union[List, None] = None) -> None:
        """Attach a death rate to one or more compartments.

        Parameters
        ----------
        rate : float
            Per-capita death rate.
        comps : list, optional
            Compartments to which to assign this death rate. If ``None``
            (the default) the rate is applied to every compartment.
        """
        self.demographics = True

        if comps is None:
            comps = self.transitions.nodes()

        for comp in comps:
            self.transitions.nodes[comp]['death']=rate

    def add_vaccination(self, source: str, target: str, rate: float, start: int) -> None:
        """Add a time-gated vaccination transition between two compartments.

        Identical to :meth:`add_spontaneous` except that the transition is
        only active for ``time > start``.

        Parameters
        ----------
        source : str
            Name of the source compartment (e.g. ``"S"``).
        target : str
            Name of the target compartment (e.g. ``"V"``).
        rate : float
            Per-capita vaccination rate.
        start : int
            Time step at which the vaccination campaign begins.
        """
        self.transitions.add_edge(source, target, rate=rate, start=start)

    def add_age_structure(self, matrix: List, population: List) -> List[List]:
        """Expand the model into an age-structured version.

        Every compartment ``X`` is replicated into one compartment per age
        group ``X_a``, ``X_b``, ... and the interactions are rescaled by the
        provided contact matrix. The ``population`` argument provides the
        absolute size of each age group, which is used by :meth:`integrate`
        to split the initial population across age strata.

        Parameters
        ----------
        matrix : list of list of float
            Square contact matrix. Entry ``matrix[i][j]`` is the average
            number of contacts a member of age group ``i`` has with members
            of age group ``j``.
        population : list of int
            Number of individuals in each age group, in the same order as
            the rows of ``matrix``.

        Raises
        ------
        ValueError
            If the model is already in multi-group (vector-borne) mode via
            :meth:`add_groups`.
        """
        if self.groups is not None:
            raise ValueError(
                "add_age_structure is incompatible with add_groups. "
                "Use one or the other."
            )

        self.contact = np.asarray(matrix)
        self.population = np.asarray(population).flatten()

        assert self.contact.shape[0] == self.contact.shape[1], "The contact matrix must be square." 

        age_groups = list(string.ascii_lowercase[:len(matrix)])
        n_ages = len(age_groups)

        model = EpiModel()
        self.orig_comps = list(self.transitions.nodes())

        for node_i, node_j, data in self.transitions.edges(data=True):
            # Interacting transition
            if "agent" in data:
                for i, age_i in enumerate(age_groups):
                    node_age_i = node_i + '_' + age_i
                    node_age_j = node_j + '_' + age_i

                    for j, age_j in enumerate(age_groups):
                        agent_age = data["agent"] + '_' + age_j

                        model.add_interaction(node_age_i, node_age_j, agent_age, data["rate"]*self.contact[i][j])

            # Spontaneous transition
            else:
                for age_i in age_groups:
                    node_age_i = node_i + '_' + age_i
                    node_age_j = node_j + '_' + age_i

                    if "start" not in data:
                        model.add_spontaneous(node_age_i, node_age_j, data["rate"])
                    else:
                        # vaccination
                        model.add_vaccination(node_age_i, node_age_j, data["rate"], data["start"])

        self.transitions = model.transitions

    def add_groups(self, groups: Dict[str, List[str]]) -> None:
        """Partition compartments into named sub-populations (groups).

        When groups are declared the model becomes a *multi-group* (e.g.
        vector-borne) model: cross-group interactions normalize the force
        of infection by the **agent's** group total rather than by the
        global total population. Concretely, an interaction

        .. math::

           X + A \\to Y + A

        contributes

        .. math::

           \\text{rate} \\cdot X \\cdot \\frac{A}{N_{\\text{group}(A)}}

        to the dynamics, where :math:`N_{\\text{group}(A)}` is the total
        population of the group that contains the agent compartment ``A``.
        This recovers the standard host-vector SIR dynamics where, for
        example, susceptible hosts become infected at rate
        :math:`b_h S_h I_v / N_v` and susceptible vectors at rate
        :math:`b_v S_v I_h / N_h`.

        Parameters
        ----------
        groups : dict of str -> list of str
            Mapping from group name (e.g. ``"host"``, ``"vector"``) to the
            list of compartment names that belong to that group. Every
            compartment present in the model must appear in exactly one
            group.

        Raises
        ------
        ValueError
            If the model is already in age-structured mode (see
            :meth:`add_age_structure`), if a compartment is listed in more
            than one group, or if any compartment in the model is left
            unassigned.

        Examples
        --------
        >>> model = EpiModel()
        >>> model.add_interaction("Sh", "Ih", "Iv", 0.3)
        >>> model.add_interaction("Sv", "Iv", "Ih", 0.3)
        >>> model.add_spontaneous("Ih", "Rh", 0.1)
        >>> model.add_groups({
        ...     "host":   ["Sh", "Ih", "Rh"],
        ...     "vector": ["Sv", "Iv"],
        ... })
        """
        if self.population is not None:
            raise ValueError(
                "add_groups is incompatible with add_age_structure. "
                "Use one or the other."
            )

        comp_to_group = {}
        for group_name, comp_list in groups.items():
            for comp in comp_list:
                if comp in comp_to_group:
                    raise ValueError(
                        "Compartment '%s' is listed in more than one group." % comp
                    )
                comp_to_group[comp] = group_name

        model_comps = set(self.transitions.nodes())
        missing = model_comps - set(comp_to_group)
        if missing:
            raise ValueError(
                "The following compartments are not assigned to any group: %s"
                % sorted(missing)
            )

        self.groups = {name: list(comps) for name, comps in groups.items()}
        self.comp_to_group = comp_to_group

    def _group_populations(self, population: np.ndarray, pos: Dict) -> Dict[str, float]:
        """Compute the per-compartment effective population.

        Internal helper used by :meth:`_new_cases` and :meth:`simulate` when
        the model has been put in multi-group mode via :meth:`add_groups`.
        Returns a mapping ``comp -> N_group(comp)``: the total population of
        the group that contains ``comp``.
        """
        group_totals = {
            name: sum(population[pos[c]] for c in comps if c in pos)
            for name, comps in self.groups.items()
        }
        return {comp: group_totals[grp] for comp, grp in self.comp_to_group.items()}

    def _new_cases(self, population: np.ndarray, time: float, pos: Dict) -> np.ndarray:
        """Compute the right-hand side of the model ODE system.

        Internal helper used by :meth:`integrate`.

        Parameters
        ----------
        population : numpy.ndarray
            Current population of each compartment.
        time : float
            Current time.
        pos : dict
            Mapping from compartment name to its index in ``population``.

        Returns
        -------
        numpy.ndarray
            Array of derivatives -- one entry per compartment.
        """
        diff = np.zeros(len(pos))
        N = np.sum(population)

        if self.population is not None:
            N = {}

            for comp_i in self.transitions.nodes():
                age_group = comp_i.split('_')[-1]

                for comp_j in pos:
                    if comp_j.endswith(age_group):
                        N[comp_i] = N.get(comp_i, 0) + population[pos[comp_j]]
        elif self.groups is not None:
            N = self._group_populations(population, pos)

        for edge in self.transitions.edges(data=True):
            source = edge[0]
            target = edge[1]
            trans = edge[2]
            
            rate = trans['rate']*population[pos[source]]
            
            if 'start' in trans and trans['start'] >= time:
                continue

            if 'agent' in trans:
                agent = trans['agent']

                if isinstance(N, dict):
                    rate *= population[pos[agent]]/N[agent]
                else:
                    rate *= population[pos[agent]]/N

                if self.seasonality is not None:
                    curr_t = int(time)%365
                    season = float(self.seasonality[curr_t])
                    rate *= season
                
            diff[pos[source]] -= rate
            diff[pos[target]] += rate

        # Population dynamics (once per RHS evaluation, not per edge)
        if self.demographics:
            for comp, data in self.transitions.nodes(data=True):
                comp_id = pos[comp]

                if "birth" in data:
                    if data.get("birth_from_group", False):
                        source_pop = N[comp]
                    else:
                        source_pop = population[comp_id]
                    births = source_pop*data["birth"]
                    diff[comp_id] += births

                if "death" in data:
                    deaths = population[comp_id]*data["death"]
                    diff[comp_id] -= deaths

        return diff
    
    def plot(self, title: Union[str, None]= None, normed: bool = True, show: bool = True, ax: Union[plt.Axes, None] = None, **kwargs):
        """Plot the time evolution of every compartment.

        Must be called after :meth:`integrate` or :meth:`simulate`.

        Parameters
        ----------
        title : str, optional
            Title to set on the resulting axes.
        normed : bool, default=True
            If True, divide every series by the total initial population so
            that values represent fractions instead of absolute counts.
        show : bool, default=True
            If True, call :func:`matplotlib.pyplot.show` after plotting.
        ax : matplotlib.axes.Axes, optional
            Axes to draw onto. If ``None``, the current axes are used.
        **kwargs
            Extra keyword arguments forwarded to
            :meth:`pandas.Series.plot`.

        Returns
        -------
        matplotlib.axes.Axes
            The axes that were drawn on.
        """
        try:
            if normed:
                N = self.values_.iloc[0].sum()
            else:
                N = 1

            if ax is None:
                ax = plt.gca()

            for comp in self.values_.columns:
                (self.values_[comp]/N).plot(c=epi_colors[comp[0]], **kwargs)

            ax.legend(self.values_.columns)
            ax.set_xlabel('Time')
            ax.set_ylabel('Population')
            
            if title is not None:
                ax.set_title(title)
            
            if show:
                plt.show()

            return ax
        except Exception as e:
            print(e)
            raise NotInitialized('You must call integrate() or simulate() first')
    
    def __getattr__(self, name: str) -> pd.Series:
        """Return the time series of a compartment as an attribute.

        After :meth:`integrate` or :meth:`simulate`, the time series of any
        compartment ``X`` can be accessed via ``model.X``.

        Parameters
        ----------
        name : str
            Name of the compartment.

        Returns
        -------
        pandas.Series
            Time series of the requested compartment.
        """

        if 'values_' in self.__dict__:
            return self.values_[name]
        else:
            raise AttributeError("'EpiModel' object has no attribute '%s'" % name)

    def simulate(self, timesteps: int, t_min: int = 1, seasonality: Union[np.ndarray, None] = None, **kwargs) -> None:
        """Run a discrete-time stochastic simulation of the model.

        Each compartment transitions according to a multinomial draw at every
        time step. Results are stored in the :attr:`values_` attribute as a
        :class:`pandas.DataFrame` indexed by time.

        Parameters
        ----------
        timesteps : int
            Number of time steps to simulate.
        t_min : int, default=1
            Starting time (used to align the index of :attr:`values_` and to
            evaluate any time-gated transitions, e.g. vaccination).
        seasonality : numpy.ndarray, optional
            Length-365 array of seasonal forcing factors. Interaction rates
            at time ``t`` are multiplied by ``seasonality[t % 365]``.
        **kwargs
            Initial size of each compartment, passed as keyword arguments
            (e.g. ``S=9990, I=10, R=0``).
        """
        pos = {comp: i for i, comp in enumerate(self.transitions.nodes())}
        population=np.zeros(len(pos), dtype='int')

        for comp in kwargs:
            population[pos[comp]] = kwargs[comp]

        values = []
        values.append(population)

        comps = list(self.transitions.nodes)
        time = np.arange(t_min, t_min+timesteps, 1, dtype='int')

        self.seasonality = seasonality

        for t in time:
            pop = values[-1]
            new_pop = values[-1].copy()

            if self.groups is None:
                N = np.sum(pop)
            else:
                N = self._group_populations(pop, pos)

            for comp in comps:
                trans = list(self.transitions.edges(comp, data=True))             

                prob = np.zeros(len(comps), dtype='float')

                for _, node_j, data in trans:
                    source = pos[comp]
                    target = pos[node_j]

                    rate = data['rate']

                    if 'start' in data and data['start'] >= t:
                        continue

                    if 'agent' in data:
                        agent_name = data['agent']
                        agent = pos[agent_name]

                        if isinstance(N, dict):
                            rate *= pop[agent]/N[agent_name]
                        else:
                            rate *= pop[agent]/N

                        if self.seasonality is not None:
                            curr_t = int(t)%365
                            season = float(self.seasonality[curr_t])
                            rate *= season

                    prob[target] = rate

                prob[source] = 1-np.sum(prob)

                delta = random.multinomial(pop[source], prob)
                delta[source] = 0

                changes = np.sum(delta)

                if changes == 0:
                    continue

                new_pop[source] -= changes

                for i in range(len(delta)):
                    new_pop[i] += delta[i]

            # Population dynamics
            if self.demographics:
                for comp, data in self.transitions.nodes(data=True):
                    comp_id = pos[comp]

                    if "birth" in data:
                        if data.get("birth_from_group", False):
                            source_pop = int(N[comp])
                        else:
                            source_pop = pop[comp_id]
                        births = np.random.binomial(source_pop, data["birth"])
                        new_pop[comp_id] += births

                    if "death" in data:
                        deaths = np.random.binomial(pop[comp_id], data["death"])
                        new_pop[comp_id] -= deaths

            values.append(new_pop)

        values = np.array(values)
        self.values_ = pd.DataFrame(values[1:], columns=comps, index=time)
    
    def integrate(self, timesteps: int , t_min: int = 1, seasonality: Union[np.ndarray, None] = None, **kwargs) -> None:
        """Numerically integrate the model's ODE system.

        Uses :func:`scipy.integrate.odeint` to advance the deterministic
        mean-field equations of the model. Results are stored in
        :attr:`values_`; when an age structure has been declared via
        :meth:`add_age_structure`, the per-age-group trajectories are also
        stored in ``values_ages_``.

        Parameters
        ----------
        timesteps : int
            Number of time steps to integrate over.
        t_min : int, default=1
            Starting time of the integration.
        seasonality : numpy.ndarray, optional
            Length-365 array of seasonal forcing factors. See
            :meth:`simulate` for details.
        **kwargs
            Initial size of each compartment, passed as keyword arguments
            (e.g. ``S=9990, I=10, R=0``).
        """

        pos = {comp: i for i, comp in enumerate(self.transitions.nodes())}
        population=np.zeros(len(pos))

        for comp in kwargs:
            if self.population is None:
                if comp not in pos:
                    continue

                population[pos[comp]] = kwargs[comp]
            else:
                total_pop = self.population.sum()
                p = np.copy(self.population)/total_pop
                n = np.random.multinomial(kwargs[comp], p, 1)[0]

                for i, age in enumerate(string.ascii_lowercase[:len(p)]):
                    comp_age = comp + '_' + age
                    if comp_age not in pos:
                        continue

                    population[pos[comp_age]] = n[i]
        
        time = np.arange(t_min, t_min+timesteps, 1)

        self.seasonality = seasonality
        values = pd.DataFrame(
                scipy.integrate.odeint(
                    self._new_cases, 
                    population, 
                    time, 
                    args=(pos,)
                ), columns=pos.keys(), index=time
        )

        if self.population is None:
            self.values_ = values
        else:
            self.values_ages_ = values

            totals = values.T.copy()
            totals['key'] = totals.index.map(lambda x: '_'.join(x.split('_')[:-1]))
            totals = totals.groupby('key').sum().T
            totals.columns.name = None
            self.values_ = totals[self.orig_comps].copy()

    def single_step(self, seasonality=None, **kwargs):
        """Advance the stochastic simulation by a single time step.

        On the first call, this is equivalent to ``simulate(2, 1, ...)``.
        On subsequent calls, the simulation is continued from the end of
        :attr:`values_` -- this is the building block used by
        :class:`MetaEpiModel` to couple multiple sub-populations together.

        Parameters
        ----------
        seasonality : numpy.ndarray, optional
            Length-365 array of seasonal forcing factors.
        **kwargs
            Initial compartment populations, forwarded to :meth:`simulate`
            on the first call only.
        """
        if hasattr(self, 'values_') is False:
            self.simulate(2, 1, seasonality=seasonality, **kwargs)
        else:
            old_values = self.values_.copy()
            t_curr = self.values_.index.max()
            self.simulate(2, t_curr, seasonality=seasonality, **kwargs)
            new_values = pd.concat([old_values, self.values_.iloc[[-1]]])
            self.values_ = new_values

    def __repr__(self) -> str:
        """Return a human-readable description of the model.

        Lists every transition (interaction, spontaneous or vaccination)
        and, when applicable, the basic reproduction number :math:`R_0`.

        Returns
        -------
        str
            Multi-line description of the model.
        """
        text = 'Epidemic Model with %u compartments and %u transitions:\n\n' % \
              (self.transitions.number_of_nodes(), 
               self.transitions.number_of_edges())
        
        for edge in self.transitions.edges(data=True):
            source = edge[0]
            target = edge[1]
            trans = edge[2]
            
            rate = trans['rate']

            if 'agent' in trans:
                agent = trans['agent']
                text += "%s + %s = %s %f\n" % (source, agent, target, rate)
            elif 'start' in trans:
                start = trans['start']
                text+="%s -> %s %f starting at %s days\n" % (source, target, rate, start)
            else:
                text+="%s -> %s %f\n" % (source, target, rate)
        
        if self.groups is not None:
            text += "\nGroups:\n"
            for group_name, comp_list in self.groups.items():
                text += "  %s: %s\n" % (group_name, ", ".join(comp_list))

        R0 = self.R0()

        if R0 is not None:
            text += "\nR0=%1.2f" % R0

        return text

    def _get_active(self) -> Set:
        active = set()

        for node_i, node_j, data in self.transitions.edges(data=True):
            if "agent" in data:
                active.add(data['agent'])
            else:
                active.add(node_i)

        return active

    def _get_susceptible(self) -> Set:
        susceptible = set([node for node, deg in self.transitions.in_degree() if deg==0])

        if len(susceptible) == 0:
            for node_i, node_j, data in self.transitions.edges(data=True):
                if "agent" in data:
                    susceptible.add(node_i)

        return susceptible

    def _get_infections(self) -> Dict:
        inf = {}

        for node_i, node_j, data in self.transitions.edges(data=True):
            if "agent" in data:
                agent = data['agent']

                if agent not in inf:
                    inf[agent] = {}

                if node_i not in inf[agent]:
                    inf[agent][node_i] = {}

                inf[agent][node_i]['target'] = node_j
                inf[agent][node_i]['rate'] = data['rate']

        return inf

    def draw_model(self, ax: Union[plt.Axes, None] = None, show: bool = True) -> None:
        """Render the compartmental structure of the model as a graph.

        Each compartment is drawn as a node and each transition as a
        directed edge labelled by the agent compartment (for interactions)
        or left unlabelled (for spontaneous transitions). When births or
        deaths have been enabled, extra invisible source/sink nodes are
        added at the top/bottom of every compartment.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Axes to draw onto. If ``None`` a new figure is created.
        show : bool, default=True
            Whether to call :func:`matplotlib.pyplot.show` after drawing.
        """

        trans = self.transitions.copy()

        try:
            from networkx.drawing.nx_agraph import graphviz_layout
            pos=graphviz_layout(trans, prog='dot', args='-Grankdir="LR"')
        except:
            pos=nx.layout.spectral_layout(trans)

        pos = nx.layout.rescale_layout_dict(pos)
        pos2 = pos.copy()

        if self.demographics:
            for comp, data in self.transitions.nodes(data=True):
                if "birth" in data:
                    trans.add_edge("_b_" + comp, comp, rate=data["birth"])
                if "death" in data:
                    trans.add_edge(comp, "_d_" + comp, rate=data["death"])

        node_colors = []

        for node in trans.nodes():
            if node.startswith("_b_"):
                node_colors.append('white')
                orig_pos = pos[node[3:]]
                pos[node] = [orig_pos[0], orig_pos[1]+1]
            elif node.startswith("_d_"):
                node_colors.append('white')
                orig_pos = pos[node[3:]]
                pos[node] = [orig_pos[0], orig_pos[1]-1]
            else:
                node_colors.append(epi_colors[node[0]])

        edge_labels = {}

        for node_i, node_j, data in trans.edges(data=True):
            edge = (node_i, node_j)

            if "agent" in data:
                if edge not in edge_labels:
                    edge_labels[edge] = data["agent"]
                else:
                    edge_labels[edge] = edge_labels[edge] + "+" + data["agent"]
            else:
                edge_labels[edge] = ""

        if ax is None:
            fig, ax = plt.subplots(1, figsize=(10, 2))

        nx.draw(trans, pos, with_labels=False, arrows=True, node_shape='H', 
        node_color=node_colors, node_size=1000, ax=ax)
        nx.draw_networkx_edge_labels(trans, pos, edge_labels=edge_labels, ax=ax)
        nx.draw_networkx_labels(self.transitions, pos2, font_color='k')

        if show:
            plt.show()

    def R0(self) -> Union[float, None]:
        r"""Compute the basic reproduction number :math:`R_0` of the model.

        The calculation is fully generic and uses the next-generation matrix
        approach described in
        *Diekmann, Heesterbeek and Roberts, J. R. Soc. Interface* **7**, 873
        (2010). :math:`R_0` is returned as the spectral radius of
        :math:`F V^{-1}`, where :math:`F` collects the rates of new
        infections in each infected compartment and :math:`V` collects the
        rates of all other transitions among infected compartments.

        Returns
        -------
        float or None
            The basic reproduction number, or ``None`` if it cannot be
            computed (for example when the model has no infectious dynamics
            defined, or when the model has been put in multi-group mode via
            :meth:`add_groups` -- vector-borne :math:`R_0` is not yet
            supported).
        """
        if self.groups is not None:
            return None

        infected = set()

        susceptible = self._get_susceptible()

        for node_i, node_j, data in self.transitions.edges(data=True):
            if "agent" in data:
                infected.add(data['agent'])
                infected.add(node_j)


        infected = sorted(infected)
        N_infected = len(infected)

        F = np.zeros((N_infected, N_infected), dtype='float')
        V = np.zeros((N_infected, N_infected), dtype='float')

        pos = dict(zip(infected, np.arange(N_infected)))

        try:
            for node_i, node_j, data in self.transitions.edges(data=True):
                rate = data['rate']

                if "agent" in data:
                    target = pos[node_j]
                    agent = pos[data['agent']]

                    if node_i in susceptible:
                        F[target, agent] = rate
                elif "start" in data:
                    continue
                else:
                    source = pos[node_i]

                    V[source, source] += rate

                    if node_j in pos:
                        target = pos[node_j]
                        V[target, source] -= rate
        
            eig, v = linalg.eig(np.dot(F, linalg.inv(V)))
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                return np.real(eig.max())
        except:
            return None

    def __getitem__(self, key):
        """Return one or several compartment trajectories by name.

        Parameters
        ----------
        key : str or list of str
            Name(s) of the compartment(s) of interest. When the model has
            age structure, age-stratified compartments such as ``"S_a"`` are
            also accepted.

        Returns
        -------
        pandas.Series or pandas.DataFrame or None
            Time series of the requested compartment(s), or ``None`` if no
            matching compartment can be found.
        """
        if type(key) != type([]):
            key_check = set([key])
        else:
            key_check = set(key)

        if len(key_check & set(self.values_.columns)) > 0:
            return self.values_[key]
        elif len(key_check & set(self.values_ages_.columns)) > 0:
            return self.values_ages_[key]
        else:
            return None
