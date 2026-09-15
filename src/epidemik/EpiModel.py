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

from .utils import NotInitialized, epi_colors


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

    Attributes
    ----------
    transitions : networkx.MultiDiGraph
        Graph of compartments (nodes) and transitions (edges).
    values_ : pandas.DataFrame
        Trajectory of every compartment, indexed by time. Only available
        after :meth:`integrate` or :meth:`simulate` has been called.
    values_ages_ : pandas.DataFrame
        Age-stratified trajectory (one column per ``compartment_age``
        pair). Only available after :meth:`add_age_structure` and a call
        to :meth:`integrate` or :meth:`simulate`.
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
        self.contact = None
        self.orig_comps = None
        self.demographics = False
        self.groups = None
        self.comp_to_group = None

        if compartments is not None:
            self.transitions.add_nodes_from([comp for comp in compartments])

    # ------------------------------------------------------------------
    # Model definition
    # ------------------------------------------------------------------
    def add_interaction(self, source: str, target: str, agent: str, rate: float) -> None:
        """Add an interaction transition between two compartments.

        Encodes a reaction of the form ``source + agent -> target + agent``
        that occurs with the given per-capita rate. In the mean-field
        equations the contribution of this transition is
        ``rate * source * agent / N`` where ``N`` is the total population
        (or the population of the agent's group/age stratum -- see
        :meth:`add_groups` and :meth:`add_age_structure`).

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

        Notes
        -----
        When combined with :meth:`add_age_structure`, birth and death
        rates attached *before* the age structure is added are copied to
        every age stratum of the compartment. Rates added *after* must use
        the expanded compartment names (``"S_a"``, ``"S_b"``, ...).
        """
        if from_group and self.groups is None:
            raise ValueError(
                "add_birth_rate(..., from_group=True) requires add_groups "
                "to have been called first."
            )

        self.demographics = True

        if comps is None:
            comps = list(self.transitions.nodes())

        for comp in comps:
            self.transitions.add_node(comp)
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
            comps = list(self.transitions.nodes())

        for comp in comps:
            self.transitions.add_node(comp)
            self.transitions.nodes[comp]['death'] = rate

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

    def add_age_structure(self, matrix: List, population: List) -> None:
        """Expand the model into an age-structured version.

        Every compartment ``X`` is replicated into one compartment per age
        group ``X_a``, ``X_b``, ... and the interactions are rescaled by the
        provided contact matrix. The ``population`` argument provides the
        absolute size of each age group, which is used by :meth:`integrate`
        and :meth:`simulate` to split the initial population across age
        strata. Any birth or death rates already attached to a compartment
        are copied to every one of its age strata.

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
            :meth:`add_groups`, if ``matrix`` is not square, or if
            ``population`` does not have one entry per age group.
        """
        if self.groups is not None:
            raise ValueError(
                "add_age_structure is incompatible with add_groups. "
                "Use one or the other."
            )

        contact = np.asarray(matrix, dtype=float)
        population = np.asarray(population, dtype=float).flatten()

        if contact.ndim != 2 or contact.shape[0] != contact.shape[1]:
            raise ValueError("The contact matrix must be square.")

        if population.shape[0] != contact.shape[0]:
            raise ValueError(
                "population must have one entry per age group "
                "(%u given, %u expected)." % (population.shape[0], contact.shape[0])
            )

        if contact.shape[0] > len(string.ascii_lowercase):
            raise ValueError("At most %u age groups are supported." % len(string.ascii_lowercase))

        self.contact = contact
        self.population = population

        age_groups = list(string.ascii_lowercase[:contact.shape[0]])

        model = EpiModel()
        self.orig_comps = list(self.transitions.nodes())

        # Create the expanded compartments first so that the ordering is
        # predictable (S_a, S_b, ..., I_a, I_b, ...) and copy over any
        # demographic attributes attached to the original compartments.
        for comp, data in self.transitions.nodes(data=True):
            for age in age_groups:
                model.transitions.add_node(comp + '_' + age, **data)

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
        r"""Partition compartments into named sub-populations (groups).

        When groups are declared the model becomes a *multi-group* (e.g.
        vector-borne) model: cross-group interactions normalize the force
        of infection by the **agent's** group total rather than by the
        global total population. Concretely, an interaction

        .. math::

           X + A \to Y + A

        contributes

        .. math::

           \text{rate} \cdot X \cdot \frac{A}{N_{\text{group}(A)}}

        to the dynamics, where :math:`N_{\text{group}(A)}` is the total
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
            group. Compartments added *after* this call must also belong
            to a group; :meth:`integrate` and :meth:`simulate` check this
            and raise :class:`ValueError` otherwise.

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

    # ------------------------------------------------------------------
    # Internal helpers shared by integrate() and simulate()
    # ------------------------------------------------------------------
    def _check_groups(self, pos: Dict) -> None:
        """Verify that every compartment belongs to a group (if groups are used)."""
        if self.groups is None:
            return

        missing = set(pos) - set(self.comp_to_group)
        if missing:
            raise ValueError(
                "The following compartments are not assigned to any group: %s. "
                "Call add_groups again after adding new compartments."
                % sorted(missing)
            )

    def _group_populations(self, population: np.ndarray, pos: Dict) -> Dict[str, float]:
        """Compute the per-compartment effective population.

        Internal helper used when the model has been put in multi-group
        mode via :meth:`add_groups`. Returns a mapping
        ``comp -> N_group(comp)``: the total population of the group that
        contains ``comp``.
        """
        group_totals = {
            name: sum(population[pos[c]] for c in comps if c in pos)
            for name, comps in self.groups.items()
        }
        return {comp: group_totals[grp] for comp, grp in self.comp_to_group.items()}

    def _effective_populations(self, population: np.ndarray, pos: Dict):
        """Return the population used to normalize interaction terms.

        Returns a scalar (the total population) for plain models, or a
        ``comp -> N`` mapping for age-structured and multi-group models,
        where ``N`` is the population of the age stratum / group that
        contains ``comp``.
        """
        if self.population is not None:
            totals = {}

            for comp in pos:
                age = comp.split('_')[-1]
                totals[age] = totals.get(age, 0) + population[pos[comp]]

            return {comp: totals[comp.split('_')[-1]] for comp in pos}
        elif self.groups is not None:
            return self._group_populations(population, pos)

        return np.sum(population)

    def _transition_rate(self, data: Dict, population: np.ndarray, pos: Dict, N, time: float):
        """Return the effective per-capita rate of a transition at ``time``.

        Returns ``None`` when the transition is not active yet (e.g. a
        vaccination campaign that has not started).
        """
        if 'start' in data and data['start'] >= time:
            return None

        rate = data['rate']

        if 'agent' in data:
            agent = data['agent']
            n = N[agent] if isinstance(N, dict) else N

            if n <= 0:
                return 0.0

            rate *= population[pos[agent]]/n

            if self.seasonality is not None:
                curr_t = int(time) % len(self.seasonality)
                rate *= float(self.seasonality[curr_t])

        return rate

    def _initial_population(self, pos: Dict, kwargs: Dict, stochastic: bool) -> np.ndarray:
        """Build the initial population vector from ``**kwargs``.

        Unknown compartment names are ignored with a warning. In
        age-structured mode, base compartment names (``"S"``) are split
        across age strata proportionally to the population vector given to
        :meth:`add_age_structure` (using a multinomial draw when
        ``stochastic`` is True); expanded names (``"S_a"``) are used
        directly.
        """
        population = np.zeros(len(pos))
        unknown = []

        if self.population is not None:
            p = self.population/self.population.sum()
            ages = string.ascii_lowercase[:len(p)]

        for comp, value in kwargs.items():
            if comp in pos:
                population[pos[comp]] = value
                continue

            if self.population is None:
                unknown.append(comp)
                continue

            comp_ages = [comp + '_' + age for age in ages]

            if not any(comp_age in pos for comp_age in comp_ages):
                unknown.append(comp)
                continue

            if stochastic:
                n = np.random.multinomial(int(value), p)
            else:
                n = value*p

            for i, comp_age in enumerate(comp_ages):
                if comp_age in pos:
                    population[pos[comp_age]] = n[i]

        if unknown:
            warnings.warn(
                "Ignoring initial conditions for unknown compartment(s): %s"
                % sorted(unknown), stacklevel=3
            )

        return population

    def _store_results(self, values: pd.DataFrame) -> None:
        """Store a trajectory in :attr:`values_` (and :attr:`values_ages_`)."""
        if self.population is None:
            self.values_ = values
        else:
            self.values_ages_ = values

            totals = values.T.copy()
            totals['key'] = ['_'.join(c.split('_')[:-1]) for c in totals.index]
            totals = totals.groupby('key', sort=False).sum().T
            totals.columns.name = None

            # Original compartments first, then any added after the expansion
            columns = list(self.orig_comps)
            columns += [c for c in totals.columns if c not in columns]
            self.values_ = totals[columns].copy()

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
        N = self._effective_populations(population, pos)

        for source, target, trans in self.transitions.edges(data=True):
            rate = self._transition_rate(trans, population, pos, N, time)

            if rate is None:
                continue

            flow = rate*population[pos[source]]

            diff[pos[source]] -= flow
            diff[pos[target]] += flow

        # Population dynamics (once per RHS evaluation, not per edge)
        if self.demographics:
            for comp, data in self.transitions.nodes(data=True):
                comp_id = pos[comp]

                if "birth" in data:
                    if data.get("birth_from_group", False):
                        source_pop = N[comp]
                    else:
                        source_pop = population[comp_id]
                    diff[comp_id] += source_pop*data["birth"]

                if "death" in data:
                    diff[comp_id] -= population[comp_id]*data["death"]

        return diff

    # ------------------------------------------------------------------
    # Results access
    # ------------------------------------------------------------------
    def plot(self, title: Union[str, None] = None, normed: bool = True, show: bool = True, ax: Union[plt.Axes, None] = None, **kwargs):
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

        Raises
        ------
        NotInitialized
            If neither :meth:`integrate` nor :meth:`simulate` has been
            called yet.
        """
        if 'values_' not in self.__dict__:
            raise NotInitialized('You must call integrate() or simulate() first')

        N = self.values_.iloc[0].sum() if normed else 1

        if N == 0:
            N = 1

        if ax is None:
            ax = plt.gca()

        for comp in self.values_.columns:
            (self.values_[comp]/N).plot(ax=ax, c=epi_colors[comp[0]], label=comp, **kwargs)

        ax.legend()
        ax.set_xlabel('Time')
        ax.set_ylabel('Population')

        if title is not None:
            ax.set_title(title)

        if show:
            plt.show()

        return ax

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

        Raises
        ------
        AttributeError
            If the model has not been run yet, or if ``name`` is not a
            compartment of the model.
        """
        if name.startswith('__'):
            raise AttributeError(name)

        if 'values_' in self.__dict__ and name in self.__dict__['values_'].columns:
            return self.__dict__['values_'][name]

        raise AttributeError("'%s' object has no attribute '%s'" % (type(self).__name__, name))

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

        Raises
        ------
        NotInitialized
            If neither :meth:`integrate` nor :meth:`simulate` has been
            called yet.
        """
        if 'values_' not in self.__dict__:
            raise NotInitialized('You must call integrate() or simulate() first')

        if isinstance(key, str):
            keys = [key]
        else:
            keys = list(key)

        if all(k in self.values_.columns for k in keys):
            return self.values_[key]

        ages = self.__dict__.get('values_ages_')

        if ages is not None and all(k in ages.columns for k in keys):
            return ages[key]

        return None

    # ------------------------------------------------------------------
    # Dynamics
    # ------------------------------------------------------------------
    def simulate(self, timesteps: int, t_min: int = 1, seasonality: Union[np.ndarray, None] = None, **kwargs) -> None:
        """Run a discrete-time stochastic simulation of the model.

        Each compartment transitions according to a multinomial draw at every
        time step, with per-step transition probabilities equal to the
        corresponding rates. Results are stored in the :attr:`values_`
        attribute as a :class:`pandas.DataFrame` indexed by time: the first
        row (``t_min``) holds the initial conditions and the remaining
        ``timesteps - 1`` rows hold the state after each step, mirroring
        :meth:`integrate`.

        Parameters
        ----------
        timesteps : int
            Number of time steps (rows of :attr:`values_`) to produce,
            including the initial condition.
        t_min : int, default=1
            Starting time (used to align the index of :attr:`values_` and to
            evaluate any time-gated transitions, e.g. vaccination).
        seasonality : numpy.ndarray, optional
            Array of seasonal forcing factors, typically of length 365.
            Interaction rates at time ``t`` are multiplied by
            ``seasonality[t % len(seasonality)]``.
        **kwargs
            Initial size of each compartment, passed as keyword arguments
            (e.g. ``S=9990, I=10, R=0``). Unknown compartment names are
            ignored with a warning. For age-structured models, values are
            split across age strata with a multinomial draw.

        Warns
        -----
        UserWarning
            If, for some compartment, the sum of the outgoing transition
            probabilities exceeds 1 at some time step. The probabilities are
            renormalized so that every individual leaves the compartment;
            consider using smaller rates (i.e. a finer time step).
        """
        if timesteps < 1:
            raise ValueError("timesteps must be a positive integer.")

        pos = {comp: i for i, comp in enumerate(self.transitions.nodes())}
        self._check_groups(pos)

        population = self._initial_population(pos, kwargs, stochastic=True).astype('int')

        comps = list(pos)
        time = np.arange(t_min, t_min+timesteps, 1, dtype='int')

        self.seasonality = seasonality

        values = [population]
        warned = False

        for t in time[1:]:
            pop = values[-1]
            new_pop = pop.copy()

            N = self._effective_populations(pop, pos)

            for comp in comps:
                source = pos[comp]

                if pop[source] == 0:
                    continue

                # One slot per compartment plus a final "death" sink so that
                # deaths compete with the other transitions out of ``comp``
                # and the population can never become negative.
                prob = np.zeros(len(comps)+1, dtype='float')

                for _, node_j, data in self.transitions.edges(comp, data=True):
                    rate = self._transition_rate(data, pop, pos, N, t)

                    if rate is None:
                        continue

                    prob[pos[node_j]] += rate

                if self.demographics:
                    prob[-1] = self.transitions.nodes[comp].get('death', 0.0)

                total = prob.sum()

                if total > 1:
                    if not warned:
                        warnings.warn(
                            "The total transition probability out of '%s' exceeds 1 at "
                            "t=%u; probabilities were renormalized. Consider using "
                            "smaller rates." % (comp, t), stacklevel=2
                        )
                        warned = True

                    prob /= total
                    total = 1.0

                prob[source] = max(0.0, 1.0-total)

                delta = random.multinomial(pop[source], prob)
                deaths = delta[-1]
                delta = delta[:-1]
                delta[source] = 0

                changes = np.sum(delta) + deaths

                if changes == 0:
                    continue

                new_pop[source] -= changes
                new_pop += delta

            # Births (deaths are handled above, as competing transitions)
            if self.demographics:
                for comp, data in self.transitions.nodes(data=True):
                    if "birth" not in data:
                        continue

                    comp_id = pos[comp]

                    if data.get("birth_from_group", False):
                        source_pop = int(N[comp])
                    else:
                        source_pop = pop[comp_id]

                    new_pop[comp_id] += np.random.binomial(source_pop, data["birth"])

            values.append(new_pop)

        values = np.array(values)
        self._store_results(pd.DataFrame(values, columns=comps, index=time))

    def integrate(self, timesteps: int, t_min: int = 1, seasonality: Union[np.ndarray, None] = None, **kwargs) -> None:
        """Numerically integrate the model's ODE system.

        Uses :func:`scipy.integrate.odeint` to advance the deterministic
        mean-field equations of the model. Results are stored in
        :attr:`values_` (the first row, at ``t_min``, holds the initial
        conditions); when an age structure has been declared via
        :meth:`add_age_structure`, the per-age-group trajectories are also
        stored in :attr:`values_ages_`.

        Parameters
        ----------
        timesteps : int
            Number of time steps to integrate over.
        t_min : int, default=1
            Starting time of the integration.
        seasonality : numpy.ndarray, optional
            Array of seasonal forcing factors, typically of length 365. See
            :meth:`simulate` for details.
        **kwargs
            Initial size of each compartment, passed as keyword arguments
            (e.g. ``S=9990, I=10, R=0``). Unknown compartment names are
            ignored with a warning. For age-structured models, values are
            split across age strata proportionally to the age group sizes.
        """
        if timesteps < 1:
            raise ValueError("timesteps must be a positive integer.")

        pos = {comp: i for i, comp in enumerate(self.transitions.nodes())}
        self._check_groups(pos)

        population = self._initial_population(pos, kwargs, stochastic=False)

        time = np.arange(t_min, t_min+timesteps, 1)

        self.seasonality = seasonality
        values = pd.DataFrame(
                scipy.integrate.odeint(
                    self._new_cases,
                    population,
                    time,
                    args=(pos,)
                ), columns=list(pos), index=time
        )

        self._store_results(values)

    def single_step(self, seasonality=None, t_min: int = 1, **kwargs):
        """Advance the stochastic simulation by a single time step.

        On the first call the initial conditions given as keyword arguments
        are stored at ``t_min`` and one stochastic step is taken. On
        subsequent calls the simulation is continued from the last row of
        :attr:`values_`, unless new compartment sizes are provided as
        keyword arguments, in which case they replace the current state
        before the step is taken -- this is the building block used by
        :class:`MetaEpiModel` to couple multiple sub-populations together.

        Parameters
        ----------
        seasonality : numpy.ndarray, optional
            Array of seasonal forcing factors, typically of length 365.
        t_min : int, default=1
            Starting time, used on the first call only.
        **kwargs
            Compartment populations. Required on the first call; optional
            afterwards.
        """
        if 'values_' not in self.__dict__:
            self.simulate(2, t_min, seasonality=seasonality, **kwargs)
            return

        old_values = self.values_
        old_ages = self.__dict__.get('values_ages_')
        t_curr = int(old_values.index.max())

        if not kwargs:
            current = old_ages if old_ages is not None else old_values
            kwargs = current.iloc[-1].to_dict()

        self.simulate(2, t_curr, seasonality=seasonality, **kwargs)

        self.values_ = pd.concat([old_values, self.values_.iloc[[-1]]])

        if old_ages is not None:
            self.values_ages_ = pd.concat([old_ages, self.values_ages_.iloc[[-1]]])

    # ------------------------------------------------------------------
    # Introspection
    # ------------------------------------------------------------------
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

        for source, target, trans in self.transitions.edges(data=True):
            rate = trans['rate']

            if 'agent' in trans:
                agent = trans['agent']
                text += "%s + %s = %s %f\n" % (source, agent, target, rate)
            elif 'start' in trans:
                start = trans['start']
                text += "%s -> %s %f starting at %s days\n" % (source, target, rate, start)
            else:
                text += "%s -> %s %f\n" % (source, target, rate)

        if self.demographics:
            for comp, data in self.transitions.nodes(data=True):
                if "birth" in data:
                    if data.get("birth_from_group", False):
                        text += "birth -> %s %f (per group total)\n" % (comp, data["birth"])
                    else:
                        text += "birth -> %s %f\n" % (comp, data["birth"])
                if "death" in data:
                    text += "%s -> death %f\n" % (comp, data["death"])

        if self.groups is not None:
            text += "\nGroups:\n"
            for group_name, comp_list in self.groups.items():
                text += "  %s: %s\n" % (group_name, ", ".join(comp_list))

        R0 = self.R0()

        if R0 is not None:
            text += "\nR0=%1.2f" % R0

        return text

    def _get_active(self) -> Set:
        """Return the set of compartments that can drive a transition."""
        active = set()

        for node_i, node_j, data in self.transitions.edges(data=True):
            if "agent" in data:
                active.add(data['agent'])
            else:
                active.add(node_i)

        return active

    def _get_susceptible(self) -> Set:
        """Return the set of compartments that behave as susceptibles.

        Compartments with no incoming transition that are the source of an
        interaction are considered susceptible. If there are none, every
        compartment without incoming transitions is used; failing that (e.g.
        in an SIS model), the source compartments of every interaction.
        """
        no_inflow = set([node for node, deg in self.transitions.in_degree() if deg == 0])
        sources = set([node_i for node_i, _, data in self.transitions.edges(data=True) if "agent" in data])

        susceptible = no_inflow & sources

        if len(susceptible) == 0:
            susceptible = no_inflow

        if len(susceptible) == 0:
            susceptible = sources

        return susceptible

    def _get_infections(self) -> Dict:
        """Return interactions indexed as ``{agent: {source: {target, rate}}}``."""
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
            pos = graphviz_layout(trans, prog='dot', args='-Grankdir="LR"')
        except Exception:
            pos = nx.layout.spectral_layout(trans)

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
                if edge not in edge_labels or edge_labels[edge] == "":
                    edge_labels[edge] = data["agent"]
                else:
                    edge_labels[edge] = edge_labels[edge] + "+" + data["agent"]
            elif edge not in edge_labels:
                edge_labels[edge] = ""

        if ax is None:
            fig, ax = plt.subplots(1, figsize=(10, 2))

        nx.draw(trans, pos, with_labels=False, arrows=True, node_shape='H',
                node_color=node_colors, node_size=1000, ax=ax)
        nx.draw_networkx_edge_labels(trans, pos, edge_labels=edge_labels, ax=ax)
        nx.draw_networkx_labels(self.transitions, pos2, font_color='k', ax=ax)

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
        Vaccination transitions and demographics are ignored.

        Returns
        -------
        float or None
            The basic reproduction number, or ``None`` if it cannot be
            computed (for example when the model has no infectious dynamics
            defined, when :math:`V` is singular, or when the model has been
            put in multi-group mode via :meth:`add_groups` -- vector-borne
            :math:`R_0` is not yet supported).
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

        if N_infected == 0:
            return None

        F = np.zeros((N_infected, N_infected), dtype='float')
        V = np.zeros((N_infected, N_infected), dtype='float')

        pos = dict(zip(infected, np.arange(N_infected)))

        for node_i, node_j, data in self.transitions.edges(data=True):
            rate = data['rate']

            if "agent" in data:
                target = pos[node_j]
                agent = pos[data['agent']]

                if node_i in susceptible:
                    F[target, agent] += rate
            elif "start" in data:
                continue
            elif node_i in pos:
                source = pos[node_i]

                V[source, source] += rate

                if node_j in pos:
                    target = pos[node_j]
                    V[target, source] -= rate

        try:
            eig = linalg.eigvals(np.dot(F, linalg.inv(V)))
        except linalg.LinAlgError:
            return None

        return float(np.max(np.abs(eig)))
