### −∗− mode : python ; −∗−
# @file EpiModel.py
# @author Bruno Goncalves
######################################################

from typing import Dict, List, Set, Union
import json
import os
import re
import shutil
import warnings
import string
from urllib.request import urlretrieve

import networkx as nx
import numpy as np
from numpy import linalg
import scipy.integrate
import pandas as pd
import matplotlib.pyplot as plt
import yaml

from . import utils
from .utils import NotInitialized, Parameters, epi_colors

# Grammar of the textual transition format used by __repr__ / save_model /
# load_model. Compartments and parameters are identifiers; a rate may also
# be written as a numeric literal.
_NAME = r"[A-Za-z_][A-Za-z0-9_]*"
_TOKEN = r"[A-Za-z_][A-Za-z0-9_]*|[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"
_RE_INTERACTION = re.compile(
    r"^\s*(%s)\s*\+\s*(%s)\s*=\s*(%s)\s+(%s)(\s+unnormalized)?\s*$" % (_NAME, _NAME, _NAME, _TOKEN)
)
_RE_VACCINATION = re.compile(
    r"^\s*(%s)\s*->\s*(%s)\s+(%s)\s+starting\s+at\s+(\d+)\s+days?\s*$" % (_NAME, _NAME, _TOKEN)
)
_RE_VIRAL = re.compile(r"^\s*(%s)\s*=>\s*(%s)\s+(%s)\s+(%s)\s*$" % (_NAME, _NAME, _TOKEN, _TOKEN))
_RE_SPONTANEOUS = re.compile(r"^\s*(%s)\s*->\s*(%s)\s+(%s)\s*$" % (_NAME, _NAME, _TOKEN))


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

    Every rate is stored as a named parameter in :attr:`params`. Rates can
    be given as plain numbers (an anonymous parameter ``rate1``, ``rate2``,
    ... is created), as a named keyword argument (``beta=0.3``), or as an
    expression in terms of other parameters (``mu="beta/2"``).

    Attributes
    ----------
    name : str or None
        Optional name of the model, stored by :meth:`save_model`.
    transitions : networkx.MultiDiGraph
        Graph of compartments (nodes) and transitions (edges).
    params : epidemik.utils.Parameters
        Named parameters of the model.
    rng : numpy.random.Generator
        Random number generator used by :meth:`simulate`.
    values_ : pandas.DataFrame
        Trajectory of every compartment, indexed by time. Only available
        after :meth:`integrate` or :meth:`simulate` has been called.
    values_ages_ : pandas.DataFrame
        Age-stratified trajectory (one column per ``compartment_age``
        pair). Only available after :meth:`add_age_structure` and a call
        to :meth:`integrate` or :meth:`simulate`.
    """

    def __init__(self, compartments=None, seed=None, rng=None, name=None):
        """Initialize an empty :class:`EpiModel` instance.

        Parameters
        ----------
        compartments : list of str, optional
            Optional list of compartment names to seed the model with. New
            compartments are also added automatically the first time they
            appear in a call to :meth:`add_interaction`,
            :meth:`add_spontaneous` or :meth:`add_vaccination`.
        seed : int, optional
            Seed for the random number generator used by :meth:`simulate`.
            Ignored when ``rng`` is given.
        rng : numpy.random.Generator, optional
            Random number generator to use. If ``None`` a new generator is
            created with :func:`numpy.random.default_rng`.
        name : str, optional
            Human-readable name of the model.
        """
        self.name = name
        self.transitions = nx.MultiDiGraph()
        self.params = Parameters()
        self.seasonality = None
        self.population = None
        self.contact = None
        self.orig_comps = None
        self.demographics = False
        self.groups = None
        self.comp_to_group = None

        if rng is None:
            rng = np.random.default_rng(seed)

        self.rng = rng

        if compartments is not None:
            self.transitions.add_nodes_from([comp for comp in compartments])

    # ------------------------------------------------------------------
    # Parameters
    # ------------------------------------------------------------------
    def define_parameters(self, **kwargs) -> None:
        """Define (or redefine) named parameters.

        Parameters
        ----------
        **kwargs
            ``name=value`` pairs where ``value`` is a number or an
            expression string in terms of other parameters (``"beta/2"``).

        Examples
        --------
        >>> model = EpiModel()
        >>> model.define_parameters(beta=0.3, mu="beta/3")
        >>> model.add_interaction("S", "I", "I", "beta")
        >>> model.add_spontaneous("I", "R", "mu")
        """
        self.params.define_parameters(**kwargs)

    def _register_rate(self, rate, rates: Dict, what: str) -> str:
        """Store a rate given positionally or as a named keyword.

        Returns the name of the parameter holding the rate. A positional
        number creates an anonymous parameter ``rate<n>``; a positional
        string naming an existing parameter is used as is; any other
        positional string is stored as an anonymous expression.
        """
        if rate is not None and rates:
            raise TypeError("%s: give the rate either positionally or as a keyword, not both." % what)

        if rate is None:
            if len(rates) != 1:
                raise TypeError(
                    "%s: exactly one rate is required, got %s." % (what, sorted(rates) or "none")
                )

            key, value = next(iter(rates.items()))
            self.params.define_parameters(**{key: value})
            return key

        if isinstance(rate, str) and rate in self.params:
            return rate

        count = 1
        while "rate%u" % count in self.params:
            count += 1

        key = "rate%u" % count
        self.params.define_parameters(**{key: rate})
        return key

    # ------------------------------------------------------------------
    # Model definition
    # ------------------------------------------------------------------
    def add_interaction(self, source: str, target: str, agent: str, rate=None, norm: bool = True, **rates) -> None:
        """Add an interaction transition between two compartments.

        Encodes a reaction of the form ``source + agent -> target + agent``.
        In the mean-field equations the contribution of this transition is
        ``rate * source * agent / N`` where ``N`` is the total population
        (or the population of the agent's group/age stratum -- see
        :meth:`add_groups` and :meth:`add_age_structure`). With
        ``norm=False`` the term is the plain mass-action product
        ``rate * source * agent``, which is the natural choice for
        within-host models (cells and virions).

        Parameters
        ----------
        source : str
            Name of the source compartment.
        target : str
            Name of the target compartment.
        agent : str
            Name of the agent compartment driving the interaction.
        rate : float or str, optional
            Per-capita rate, given positionally: a number, the name of an
            existing parameter, or an expression string.
        norm : bool, default=True
            Whether to normalize the interaction by the population size.
        **rates
            Alternatively, the rate as a single named parameter, e.g.
            ``beta=0.3``.
        """
        key = self._register_rate(rate, rates, "add_interaction")
        # Keep the natural compartment order (source, target, then agent)
        self.transitions.add_nodes_from([source, target, agent])
        self.transitions.add_edge(source, target, agent=agent, rate=key, norm=norm)

    def add_spontaneous(self, source: str, target: str, rate=None, **rates) -> None:
        """Add a spontaneous transition between two compartments.

        Encodes a reaction of the form ``source -> target`` that occurs at
        the given rate, independently of any other compartment.

        Parameters
        ----------
        source : str
            Name of the source compartment.
        target : str
            Name of the target compartment.
        rate : float or str, optional
            Per-capita rate, given positionally (number, parameter name or
            expression).
        **rates
            Alternatively, the rate as a single named parameter, e.g.
            ``mu=0.1``.
        """
        key = self._register_rate(rate, rates, "add_spontaneous")
        self.transitions.add_edge(source, target, rate=key)

    def add_viral_generation(self, source: str, target: str, source_rate=None, target_rate=None, **rates) -> None:
        """Add an asymmetric production transition (e.g. virions from infected cells).

        Encodes ``source => target``: individuals in ``source`` are removed
        at rate ``source_rate`` (cell death) while, independently, new
        individuals appear in ``target`` at rate ``target_rate`` per
        ``source`` individual (viral production). Unlike a spontaneous
        transition the two flows are not equal.

        Parameters
        ----------
        source : str
            Producing compartment (e.g. infected cells).
        target : str
            Produced compartment (e.g. free virus).
        source_rate : float or str, optional
            Removal rate of ``source`` individuals.
        target_rate : float or str, optional
            Production rate of ``target`` individuals per ``source``
            individual.
        **rates
            Alternatively, exactly two named parameters, the first being the
            source rate and the second the target rate, e.g.
            ``delta=0.5, p=100``.
        """
        if source_rate is None and target_rate is None:
            if len(rates) != 2:
                raise TypeError("add_viral_generation requires exactly two rates.")
            (source_name, source_value), (target_name, target_value) = rates.items()
            source_key = self._register_rate(None, {source_name: source_value}, "add_viral_generation")
            target_key = self._register_rate(None, {target_name: target_value}, "add_viral_generation")
        elif source_rate is None or target_rate is None or rates:
            raise TypeError("add_viral_generation requires exactly two rates.")
        else:
            source_key = self._register_rate(source_rate, {}, "add_viral_generation")
            target_key = self._register_rate(target_rate, {}, "add_viral_generation")

        self.transitions.add_edge(source, target, rate=source_key, viral_source=True, viral_target=False,
                                  viral_pair=(source_key, target_key))
        self.transitions.add_edge(source, target, rate=target_key, viral_source=False, viral_target=True,
                                  viral_pair=(source_key, target_key))

    def _split_rate_and_comps(self, rate, comps, rates):
        """Accept ``add_birth_rate("S", b=0.1)``-style calls (compartments first)."""
        if isinstance(rate, (list, tuple, set)):
            return None, list(rate)

        if isinstance(rate, str) and rates and rate in self.transitions.nodes and comps is None:
            return None, [rate]

        return rate, comps

    def add_birth_rate(self, rate=None, comps: Union[List, None] = None, fixed: bool = False,
                       global_rate: bool = True, from_group: bool = False, **rates) -> None:
        """Attach a birth rate to one or more compartments.

        Enabling births (or deaths) switches the model into demographics
        mode: inflows are added to the chosen compartments during both
        :meth:`integrate` and :meth:`simulate`. The inflow into each listed
        compartment ``X`` is, per unit time:

        * ``rate * N`` with ``global_rate=True`` (the default), where ``N``
          is the total population -- the classical setup where every
          newborn enters the susceptible compartment;
        * ``rate * X`` with ``global_rate=False``;
        * ``rate * N_group(X)`` with ``from_group=True``, where
          ``N_group(X)`` is the population of the group containing ``X``
          (see :meth:`add_groups`) -- the canonical setup for vector-borne
          models where all newborn vectors are susceptible;
        * exactly ``rate`` individuals with ``fixed=True``.

        Parameters
        ----------
        rate : float or str, optional
            Birth rate (number, parameter name or expression). For backwards
            compatibility, a list of compartments may also be given here
            together with a named rate (``add_birth_rate(["S"], b=0.1)``).
        comps : list, optional
            Compartments to which to assign this birth rate. If ``None``
            (the default) the rate is applied to every compartment.
        fixed : bool, default=False
            If True, ``rate`` is a constant number of births per unit time.
        global_rate : bool, default=True
            If True, births are proportional to the total population.
        from_group : bool, default=False
            If True, births are proportional to the group total. Requires
            :meth:`add_groups`.
        **rates
            Alternatively, the rate as a single named parameter, e.g.
            ``b=0.01``.

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
        rate, comps = self._split_rate_and_comps(rate, comps, rates)

        if from_group and self.groups is None:
            raise ValueError(
                "add_birth_rate(..., from_group=True) requires add_groups "
                "to have been called first."
            )

        key = self._register_rate(rate, rates, "add_birth_rate")
        self.demographics = True

        if comps is None:
            comps = list(self.transitions.nodes())

        for comp in comps:
            self.transitions.add_node(comp)
            self.transitions.nodes[comp]['birth'] = key
            self.transitions.nodes[comp]['fixed_birth'] = fixed
            self.transitions.nodes[comp]['global_birth'] = bool(global_rate) and not from_group
            self.transitions.nodes[comp]['birth_from_group'] = from_group

    def add_death_rate(self, rate=None, comps: Union[List, None] = None, fixed: bool = False,
                       global_rate: bool = False, **rates) -> None:
        """Attach a death rate to one or more compartments.

        The outflow from each listed compartment ``X`` is, per unit time,
        ``rate * X`` by default, ``rate * N`` (with ``N`` the total
        population) with ``global_rate=True``, or exactly ``rate``
        individuals with ``fixed=True``.

        Parameters
        ----------
        rate : float or str, optional
            Death rate (number, parameter name or expression). A list of
            compartments may also be given here together with a named rate.
        comps : list, optional
            Compartments to which to assign this death rate. If ``None``
            (the default) the rate is applied to every compartment.
        fixed : bool, default=False
            If True, ``rate`` is a constant number of deaths per unit time.
        global_rate : bool, default=False
            If True, deaths are proportional to the total population.
        **rates
            Alternatively, the rate as a single named parameter, e.g.
            ``d=0.01``.
        """
        rate, comps = self._split_rate_and_comps(rate, comps, rates)

        key = self._register_rate(rate, rates, "add_death_rate")
        self.demographics = True

        if comps is None:
            comps = list(self.transitions.nodes())

        for comp in comps:
            self.transitions.add_node(comp)
            self.transitions.nodes[comp]['death'] = key
            self.transitions.nodes[comp]['fixed_death'] = fixed
            self.transitions.nodes[comp]['global_death'] = bool(global_rate)

    def add_vaccination(self, source: str, target: str, rate=None, start: Union[int, None] = None, **rates) -> None:
        """Add a time-gated vaccination transition between two compartments.

        Identical to :meth:`add_spontaneous` except that the transition is
        only active for ``time > start``.

        Parameters
        ----------
        source : str
            Name of the source compartment (e.g. ``"S"``).
        target : str
            Name of the target compartment (e.g. ``"V"``).
        rate : float or str, optional
            Per-capita vaccination rate (number, parameter name or
            expression).
        start : int
            Time step at which the vaccination campaign begins.
        **rates
            Alternatively, the rate as a single named parameter, e.g.
            ``nu=0.01``.
        """
        if start is None:
            raise TypeError("add_vaccination requires the 'start' time of the campaign.")

        key = self._register_rate(rate, rates, "add_vaccination")
        self.transitions.add_edge(source, target, rate=key, start=start)

    def add_age_structure(self, matrix: List, population: List) -> None:
        """Expand the model into an age-structured version.

        Every compartment ``X`` is replicated into one compartment per age
        group ``X_a``, ``X_b``, ... and the interactions are rescaled by the
        provided contact matrix. The ``population`` argument provides the
        absolute size of each age group, which is used by :meth:`integrate`
        and :meth:`simulate` to split the initial population across age
        strata. Any birth or death rates already attached to a compartment
        are copied to every one of its age strata. The rescaled interaction
        rates are stored as new parameters named ``<rate>_<i><j>`` defined
        in terms of the original ones, so the original parameters can still
        be changed afterwards.

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

        model = EpiModel(rng=self.rng)
        model.params = self.params.copy()
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
                        key = "%s_%s%s" % (data["rate"], age_i, age_j)
                        model.params.define_parameters(**{key: "(%s)*%r" % (data["rate"], float(contact[i][j]))})
                        model.transitions.add_edge(node_age_i, node_age_j, agent=agent_age, rate=key,
                                                   norm=data.get("norm", True))

            # Spontaneous, vaccination or viral transition: same parameters
            else:
                for age_i in age_groups:
                    node_age_i = node_i + '_' + age_i
                    node_age_j = node_j + '_' + age_i
                    model.transitions.add_edge(node_age_i, node_age_j, **data)

        self.transitions = model.transitions
        self.params = model.params

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
        """Return ``comp -> N_group(comp)`` for multi-group models."""
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

    def _transition_rate(self, data: Dict, population: np.ndarray, pos: Dict, N, time: float, rates: Dict):
        """Return the effective per-capita rate of a transition at ``time``.

        Returns ``None`` when the transition is not active yet (e.g. a
        vaccination campaign that has not started).
        """
        if 'start' in data and data['start'] >= time:
            return None

        rate = rates[data['rate']]

        if 'agent' in data:
            agent = data['agent']
            amount = population[pos[agent]]

            if data.get('norm', True):
                n = N[agent] if isinstance(N, dict) else N

                if n <= 0:
                    return 0.0

                amount = amount/n

            rate *= amount

            if self.seasonality is not None:
                curr_t = int(time) % len(self.seasonality)
                rate *= float(self.seasonality[curr_t])

        return rate

    def _expected_births(self, comp: str, data: Dict, population: np.ndarray, pos: Dict, N, rates: Dict) -> float:
        """Expected number of births into ``comp`` per unit time."""
        rate = rates[data['birth']]

        if data.get('fixed_birth', False):
            return rate
        if data.get('birth_from_group', False):
            return rate*N[comp]
        if data.get('global_birth', False):
            return rate*np.sum(population)

        return rate*population[pos[comp]]

    def _expected_deaths(self, comp: str, data: Dict, population: np.ndarray, pos: Dict, rates: Dict) -> float:
        """Expected number of deaths out of ``comp`` per unit time."""
        rate = rates[data['death']]

        if data.get('fixed_death', False):
            return rate
        if data.get('global_death', False):
            return rate*np.sum(population)

        return rate*population[pos[comp]]

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
                n = self.rng.multinomial(int(value), p)
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

    def _new_cases(self, population: np.ndarray, time: float, pos: Dict, rates: Dict) -> np.ndarray:
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
        rates : dict
            Evaluated parameters (see :meth:`Parameters.resolve
            <epidemik.utils.Parameters.resolve>`).

        Returns
        -------
        numpy.ndarray
            Array of derivatives -- one entry per compartment.
        """
        diff = np.zeros(len(pos))
        N = self._effective_populations(population, pos)

        for source, target, trans in self.transitions.edges(data=True):
            rate = self._transition_rate(trans, population, pos, N, time, rates)

            if rate is None:
                continue

            flow = rate*population[pos[source]]

            # Viral generation edges are asymmetric: the "source" edge only
            # removes from the source, the "target" edge only produces.
            if not trans.get('viral_target', False):
                diff[pos[source]] -= flow

            if not trans.get('viral_source', False):
                diff[pos[target]] += flow

        # Population dynamics (once per RHS evaluation, not per edge)
        if self.demographics:
            for comp, data in self.transitions.nodes(data=True):
                comp_id = pos[comp]

                if "birth" in data:
                    diff[comp_id] += self._expected_births(comp, data, population, pos, N, rates)

                if "death" in data:
                    diff[comp_id] -= self._expected_deaths(comp, data, population, pos, rates)

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

        At every time step the individuals of each compartment are
        redistributed with a multinomial draw whose probabilities are the
        outgoing transition rates (deaths included); births are drawn from a
        binomial (or Poisson, for fixed rates) distribution and viral
        production from a Poisson distribution. Results are stored in the
        :attr:`values_` attribute as a :class:`pandas.DataFrame` indexed by
        time: the first row (``t_min``) holds the initial conditions and the
        remaining ``timesteps - 1`` rows hold the state after each step,
        mirroring :meth:`integrate`. All random numbers come from
        :attr:`rng`.

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
        rates = self.params.resolve()

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
                production = []

                for _, node_j, data in self.transitions.edges(comp, data=True):
                    rate = self._transition_rate(data, pop, pos, N, t, rates)

                    if rate is None:
                        continue

                    if data.get('viral_target', False):
                        production.append((pos[node_j], rate*pop[source]))
                    elif data.get('viral_source', False):
                        prob[-1] += rate
                    else:
                        prob[pos[node_j]] += rate

                if self.demographics and 'death' in self.transitions.nodes[comp]:
                    expected = self._expected_deaths(comp, self.transitions.nodes[comp], pop, pos, rates)
                    prob[-1] += expected/pop[source]

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

                delta = self.rng.multinomial(pop[source], prob)
                deaths = delta[-1]
                delta = delta[:-1]
                delta[source] = 0

                changes = np.sum(delta) + deaths

                if changes > 0:
                    new_pop[source] -= changes
                    new_pop += delta

                for target, expected in production:
                    new_pop[target] += self.rng.poisson(expected)

            # Births (deaths are handled above, as competing transitions)
            if self.demographics:
                for comp, data in self.transitions.nodes(data=True):
                    if "birth" not in data:
                        continue

                    comp_id = pos[comp]
                    rate = rates[data['birth']]

                    if data.get('fixed_birth', False):
                        births = self.rng.poisson(rate)
                    elif data.get('birth_from_group', False):
                        births = self.rng.binomial(int(N[comp]), min(rate, 1.0))
                    elif data.get('global_birth', False):
                        births = self.rng.binomial(int(np.sum(pop)), min(rate, 1.0))
                    else:
                        births = self.rng.binomial(int(pop[comp_id]), min(rate, 1.0))

                    new_pop[comp_id] += births

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
        rates = self.params.resolve()

        population = self._initial_population(pos, kwargs, stochastic=False)

        time = np.arange(t_min, t_min+timesteps, 1)

        self.seasonality = seasonality
        values = pd.DataFrame(
                scipy.integrate.odeint(
                    self._new_cases,
                    population,
                    time,
                    args=(pos, rates)
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

    def reset(self) -> None:
        """Discard the results of any previous :meth:`integrate` or :meth:`simulate` call."""
        self.__dict__.pop('values_', None)
        self.__dict__.pop('values_ages_', None)

    # ------------------------------------------------------------------
    # Introspection
    # ------------------------------------------------------------------
    @staticmethod
    def _format_value(value) -> str:
        if isinstance(value, str):
            return json.dumps(value)
        return repr(float(value))

    def __repr__(self) -> str:
        """Return a human-readable, YAML-compatible description of the model.

        The text lists the parameters, every transition (interaction,
        spontaneous, vaccination or viral generation), the demographic
        rates and groups, and -- as a comment -- the basic reproduction
        number :math:`R_0`. It is exactly what :meth:`save_model` writes
        and what :meth:`load_model` reads back.

        Returns
        -------
        str
            Multi-line description of the model.
        """
        lines = ['# Epidemic Model with %u compartments and %u transitions:' % (
            self.transitions.number_of_nodes(), self.transitions.number_of_edges())]

        if self.name is not None:
            lines.append('Name: %s' % self.name)

        lines += ['', 'Compartments: [%s]' % ', '.join(self.transitions.nodes())]

        lines += ['', 'Parameters:']
        for key in self.params.keys():
            lines.append('  %s: %s' % (key, self._format_value(self.params.raw(key))))

        lines += ['', 'Transitions:']
        for source, target, trans in self.transitions.edges(data=True):
            if 'agent' in trans:
                line = '  - %s + %s = %s %s' % (source, trans['agent'], target, trans['rate'])
                if not trans.get('norm', True):
                    line += ' unnormalized'
                lines.append(line)
            elif 'start' in trans:
                lines.append('  - %s -> %s %s starting at %s days' % (source, target, trans['rate'], trans['start']))
            elif trans.get('viral_source', False):
                source_key, target_key = trans['viral_pair']
                lines.append('  - %s => %s %s %s' % (source, target, source_key, target_key))
            elif trans.get('viral_target', False):
                continue
            else:
                lines.append('  - %s -> %s %s' % (source, target, trans['rate']))

        if self.demographics:
            lines += ['', 'Demographics:']
            for comp, data in self.transitions.nodes(data=True):
                if 'birth' in data:
                    entry = 'compartment: %s, birth: %s' % (comp, data['birth'])
                    if data.get('fixed_birth', False):
                        entry += ', fixed: true'
                    if data.get('birth_from_group', False):
                        entry += ', from_group: true'
                    elif not data.get('global_birth', False):
                        entry += ', global: false'
                    lines.append('  - {%s}' % entry)
                if 'death' in data:
                    entry = 'compartment: %s, death: %s' % (comp, data['death'])
                    if data.get('fixed_death', False):
                        entry += ', fixed: true'
                    if data.get('global_death', False):
                        entry += ', global: true'
                    lines.append('  - {%s}' % entry)

        if self.groups is not None:
            lines += ['', 'Groups:']
            for group_name, comp_list in self.groups.items():
                lines.append('  %s: [%s]' % (group_name, ', '.join(comp_list)))

        R0 = self.R0()

        if R0 is not None:
            lines += ['', '# R0=%1.2f' % R0]

        return '\n'.join(lines) + '\n'

    def _get_active(self) -> Set:
        """Return the set of compartments that can drive a transition."""
        active = set()

        for node_i, node_j, data in self.transitions.edges(data=True):
            if "agent" in data:
                active.add(data['agent'])
            else:
                active.add(node_i)

        return active

    def _get_infected(self) -> Set:
        """Return the compartments that take part in transmission (agents and targets)."""
        infected = set()

        for node_i, node_j, data in self.transitions.edges(data=True):
            if "agent" in data:
                infected.add(data['agent'])
                infected.add(node_j)

        return infected

    def _get_susceptible(self) -> Set:
        """Return the set of compartments that behave as susceptibles.

        A compartment is susceptible if it is the source of an interaction,
        is not itself infected (an agent or the target of an interaction),
        and only receives individuals from infected compartments or from
        compartments downstream of the infection (e.g. ``R -> S``). In
        particular a class that is fed by vaccination of the susceptibles
        (``S -> V``) is *not* susceptible, so that :math:`R_0` describes
        the disease-free state before vaccination. If no compartment
        qualifies, every non-infected interaction source is used.
        """
        infected = self._get_infected()
        sources = set([node_i for node_i, _, data in self.transitions.edges(data=True) if "agent" in data])

        # Compartments reachable from the infected ones through spontaneous
        # transitions without passing through an interaction source.
        removed = set()
        stack = list(infected)

        while stack:
            node = stack.pop()

            for _, nxt, data in self.transitions.out_edges(node, data=True):
                if "agent" in data or nxt in infected or nxt in sources or nxt in removed:
                    continue

                removed.add(nxt)
                stack.append(nxt)

        susceptible = set()

        for comp in sources - infected:
            upstream = set([pred for pred, _, _ in self.transitions.in_edges(comp, data=True)])

            if upstream <= (infected | removed):
                susceptible.add(comp)

        if len(susceptible) == 0:
            susceptible = sources - infected

        if len(susceptible) == 0:
            susceptible = sources

        return susceptible

    def _get_infections(self) -> Dict:
        """Return interactions indexed as ``{agent: {source: [{target, rate}, ...]}}``."""
        inf = {}

        for node_i, node_j, data in self.transitions.edges(data=True):
            if "agent" in data:
                agent = data['agent']
                inf.setdefault(agent, {}).setdefault(node_i, []).append(
                    {'target': node_j, 'rate': data['rate']}
                )

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
        rates of all other transitions among infected compartments. Every
        susceptible class (see :meth:`_get_susceptible`) is assumed to hold
        the whole population at the disease-free equilibrium. Vaccination
        transitions and demographics are ignored.

        Returns
        -------
        float or None
            The basic reproduction number, or ``None`` if it cannot be
            computed: when the model has no infectious dynamics defined,
            when :math:`V` is singular, when the model has been put in
            multi-group mode via :meth:`add_groups`, or when it contains
            unnormalized interactions or viral generation transitions
            (within-host models), whose :math:`R_0` depends on the initial
            conditions.
        """
        if self.groups is not None:
            return None

        for _, _, data in self.transitions.edges(data=True):
            if ('agent' in data and not data.get('norm', True)) or data.get('viral_source') or data.get('viral_target'):
                return None

        try:
            rates = self.params.resolve()
        except (KeyError, ValueError):
            return None

        infected = sorted(self._get_infected())
        N_infected = len(infected)

        if N_infected == 0:
            return None

        susceptible = self._get_susceptible()

        F = np.zeros((N_infected, N_infected), dtype='float')
        V = np.zeros((N_infected, N_infected), dtype='float')

        pos = dict(zip(infected, np.arange(N_infected)))

        for node_i, node_j, data in self.transitions.edges(data=True):
            rate = rates[data['rate']]

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

    # ------------------------------------------------------------------
    # Persistence
    # ------------------------------------------------------------------
    def save_model(self, filename: str) -> None:
        """Save the model definition to a YAML file.

        The file contains the output of ``repr(model)``: the name,
        parameters, transitions, demographics and groups. It can be read
        back with :meth:`load_model`.

        Parameters
        ----------
        filename : str
            Path of the file to write.
        """
        with open(filename, "wt") as f:
            f.write(repr(self))

    @staticmethod
    def _rate_argument(token: str, params: Dict) -> Dict:
        """Translate a rate token of the text format into ``add_*`` arguments."""
        if re.match(r"^%s$" % _NAME, token):
            if token not in params:
                raise ValueError("Transition refers to undefined parameter '%s'." % token)
            return {token: params[token]}

        return {"rate": float(token)}

    @staticmethod
    def load_model(filename: str, seed=None, rng=None) -> "EpiModel":
        """Load a model definition from a YAML file written by :meth:`save_model`.

        The file is a mapping with the (case-insensitive) keys ``Name``,
        ``Compartments`` (optional list fixing the compartment order),
        ``Parameters`` (``name: value`` pairs, values being numbers or
        expressions), ``Transitions`` (a list of strings such as
        ``S + I = I beta``, ``I -> R mu``, ``S -> V nu starting at 60 days``
        or ``I => V delta p``), and optionally ``Demographics`` (a list of
        ``{compartment: S, birth: b}`` / ``{compartment: I, death: d}``
        mappings with optional ``fixed``, ``global`` and ``from_group``
        flags) and ``Groups`` (a mapping from group name to list of
        compartments).

        Parameters
        ----------
        filename : str
            Path of the file to read.
        seed : int, optional
            Seed for the random number generator of the new model.
        rng : numpy.random.Generator, optional
            Random number generator for the new model.

        Returns
        -------
        EpiModel
            The loaded model.

        Raises
        ------
        ValueError
            If a transition cannot be parsed or refers to an undefined
            parameter.
        """
        with open(filename, "rt") as f:
            data = yaml.safe_load(f) or {}

        data = {str(key).lower(): value for key, value in data.items()}

        model = EpiModel(seed=seed, rng=rng, name=data.get("name"))
        model.transitions.add_nodes_from([str(c) for c in data.get("compartments") or []])

        params = {str(k): v for k, v in (data.get("parameters") or {}).items()}
        model.define_parameters(**params)

        for trans in data.get("transitions") or []:
            trans = str(trans)

            match = _RE_INTERACTION.match(trans)
            if match:
                source, agent, target, rate, unnormalized = match.groups()
                model.add_interaction(source, target, agent, norm=unnormalized is None,
                                      **EpiModel._rate_argument(rate, params))
                continue

            match = _RE_VACCINATION.match(trans)
            if match:
                source, target, rate, start = match.groups()
                model.add_vaccination(source, target, start=int(start), **EpiModel._rate_argument(rate, params))
                continue

            match = _RE_VIRAL.match(trans)
            if match:
                source, target, source_rate, target_rate = match.groups()
                source_arg = EpiModel._rate_argument(source_rate, params)
                target_arg = EpiModel._rate_argument(target_rate, params)
                if "rate" in source_arg or "rate" in target_arg:
                    model.add_viral_generation(source, target, float(source_rate), float(target_rate))
                else:
                    model.add_viral_generation(source, target, **source_arg, **target_arg)
                continue

            match = _RE_SPONTANEOUS.match(trans)
            if match:
                source, target, rate = match.groups()
                model.add_spontaneous(source, target, **EpiModel._rate_argument(rate, params))
                continue

            raise ValueError("Unable to parse transition %r." % trans)

        groups = data.get("groups")
        if groups:
            model.add_groups({str(name): [str(c) for c in comps] for name, comps in groups.items()})

        for entry in data.get("demographics") or []:
            entry = {str(k).lower(): v for k, v in entry.items()}
            comps = [str(entry["compartment"])]

            if "birth" in entry:
                model.add_birth_rate(comps=comps, fixed=bool(entry.get("fixed", False)),
                                     global_rate=bool(entry.get("global", True)),
                                     from_group=bool(entry.get("from_group", False)),
                                     **EpiModel._rate_argument(str(entry["birth"]), params))
            if "death" in entry:
                model.add_death_rate(comps=comps, fixed=bool(entry.get("fixed", False)),
                                     global_rate=bool(entry.get("global", False)),
                                     **EpiModel._rate_argument(str(entry["death"]), params))

        return model

    @staticmethod
    def _fetch(remote_path: str, filename: str, force: bool = False) -> str:
        """Copy ``filename`` from a repository into the cache and return its local path."""
        cache_dir = utils.get_cache_directory()
        os.makedirs(cache_dir, exist_ok=True)
        local_path = os.path.join(cache_dir, filename)

        if force or not os.path.exists(local_path):
            if os.path.isdir(remote_path):
                shutil.copyfile(os.path.join(remote_path, filename), local_path)
            else:
                urlretrieve(remote_path + filename, local_path)

        return local_path

    @staticmethod
    def list_models(repo: Union[str, None] = None) -> List[str]:
        """List the model files available in a repository.

        Parameters
        ----------
        repo : str, optional
            GitHub repository URL (defaults to the official ``epidemik``
            repository) or local directory. A GitHub repository must
            provide a ``model_list.txt`` index next to the model files.

        Returns
        -------
        list of str
            Sorted file names, including any ``.yaml`` files previously
            saved in the local cache directory (see
            :func:`epidemik.utils.get_cache_directory`).
        """
        remote_path = utils.get_remote_path(repo)

        if os.path.isdir(remote_path):
            models = [f for f in os.listdir(remote_path) if f.endswith(".yaml")]
        else:
            local_path = EpiModel._fetch(remote_path, "model_list.txt", force=True)
            with open(local_path, "rt") as f:
                models = [line.strip() for line in f if line.strip()]

        cache_dir = utils.get_cache_directory()

        if os.path.isdir(cache_dir):
            models += [f for f in os.listdir(cache_dir) if f.endswith(".yaml")]

        return sorted(set(models))

    @staticmethod
    def download_model(filename: str, repo: Union[str, None] = None, load_model: bool = True,
                       force: bool = False, seed=None, rng=None):
        """Fetch a model file from a repository (caching it locally).

        Parameters
        ----------
        filename : str
            Name of the model file, e.g. ``"SIR.yaml"`` (see
            :meth:`list_models`).
        repo : str, optional
            GitHub repository URL (defaults to the official ``epidemik``
            repository) or local directory.
        load_model : bool, default=True
            If True, return the loaded :class:`EpiModel`; otherwise return
            the path of the cached file.
        force : bool, default=False
            Re-download the file even if it is already cached.
        seed : int, optional
            Seed for the random number generator of the loaded model.
        rng : numpy.random.Generator, optional
            Random number generator for the loaded model.

        Returns
        -------
        EpiModel or str
            The loaded model, or the local path of the cached file.
        """
        remote_path = utils.get_remote_path(repo)
        local_path = EpiModel._fetch(remote_path, filename, force=force)

        if load_model:
            return EpiModel.load_model(local_path, seed=seed, rng=rng)

        return local_path
