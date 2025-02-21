### −∗− mode : python ; −∗−
# @file EpiModel.py
# @author Bruno Goncalves
######################################################

from typing import Dict, List, Set, Union, Self
import warnings
import string
import time
import os
import re
from urllib.parse import urlparse
from urllib.request import urlretrieve

import networkx as nx
import numpy as np
from numpy import linalg
import scipy.integrate
import pandas as pd
import matplotlib.pyplot as plt
import yaml

from typing import Union

from . import utils


class EpiModel(object):
    """Simple Epidemic Model Implementation

    Provides a way to implement and numerically integrate
    """

    def __init__(self, compartments=None, seed=None, rng=None):
        """
        Initialize the EpiModel object

        Parameters:
        - compartments: list of strings, optional
            List of compartment names

        Returns:
        None
        """
        self.name = None
        self.transitions = nx.MultiDiGraph()
        self.seasonality = None
        self.population = None
        self.orig_comps = None
        self.demographics = False
        self.params = utils.Parameters()

        if seed is None:
            seed = int(time.time()) + os.getpid()

        if rng is None:
            self.rng = np.random.default_rng(seed=seed)
        else:
            self.rng = rng

        if compartments is not None:
            self.transitions.add_nodes_from([comp for comp in compartments])

    def add_interaction(
        self,
        source: str,
        target: str,
        agent: str,
        rate: Union[None, float, str] = None,
        norm=True,
        **rates,
    ) -> None:
        """
        Add an interaction between two compartments

        Parameters:
        - source: string
            Name of the source compartment
        - target: string
            Name of the target compartment
        - agent: string
            Name of the agent
        - rate: float, str, None
            Rate of the interaction
        - norm: bool
            Whether to normalize the transition rate or not
        - rates:
            Named parameters for the interaction

        Returns:
        None
        """

        if rate is not None:
            count = len(self.params) + 1
            rate_key = "rate" + str(count)
            self.params[rate_key]=rate
        else:
            self.params.define_parameters(**rates)
            rates = list(rates.keys())
            rate_key = rates[0]

        if agent not in self.transitions.nodes:
            self.transitions.add_node(agent)
    
        self.transitions.add_edge(source, target, agent=agent, rate=rate_key, norm=norm)

    def add_spontaneous(
        self, source: str, target: str, rate: Union[None, float, str] = None, **rates
    ) -> None:
        """
        Add a spontaneous transition between two compartments

        Parameters:
        - source: string
            Name of the source compartment
        - target: string
            Name of the target compartment
        - rate: float
            Rate of the transition

        Returns:
        None
        """

        if rate is not None:
            count = len(self.params) + 1
            rate_key = "rate" + str(count)
            self.params[rate_key, rate]
        else:
            self.params.define_parameters(**rates)
            rates = list(rates.keys())
            rate_key = rates[0]

        self.transitions.add_edge(source, target, rate=rate_key)

    def add_viral_generation(
        self,
        source: str,
        target: str,
        source_rate: Union[float, str, None] = None,
        target_rate: Union[float, str, None] = None,
        **rates,
    ) -> None:
        """
        Add a viral generation transition

        Parameters:
        - source: string
            Name of the source compartment
        - target: string
            Name of the target compartment
        - source_rate: float
            Rate of destruction of infected cells
        - target_rate: float
            Rate of creation of viral particles
        """

        if source_rate is not None and target_rate is not None:
            count = len(self.params) + 1
            rate_key = "rate" + str(count)
            self.params[rate_key] = source_rate

            rate_key = "rate" + str(count + 1)
            self.params[rate_key] = target_rate
        else:
            self.params.define_parameters(**rates)
            rates = list(rates.keys())
            source_rate = rates[0]
            target_rate = rates[1]

        self.transitions.add_edge(
            source, target, rate=source_rate, viral_source=True, viral_target=False,
        )
        self.transitions.add_edge(
            source, target, rate=target_rate, viral_source=False, viral_target=True
        )

    def add_birth_rate(
        self, 
        comps: Union[List, None] = None, 
        rate: Union[float, None] = None, 
        fixed=False, 
        global_rate=True,
        **rates
    ) -> None:
        """
        Add a birth rate to one or more compartments

        Parameters:
        - rate: float
            Birth rate
        - comps: list, optional, default=None
            List of compartments to which to assign this birth rate.
            If None, apply to all compartments
        """
        self.demographics = True

        if rate is not None:
            count = len(self.params) + 1
            rate_key = "rate" + str(count)
            self.params[rate_key] = rate
        else:
            self.params.update(rates)
            rate_key = list(rates.keys())[0]

        if comps is None:
            comps = self.transitions.nodes()

        for comp in comps:
            if comp not in self.transitions.nodes:
                self.transitions.add_node(comp)

            self.transitions.nodes[comp]["birth"] = rate_key
            self.transitions.nodes[comp]["fixed"] = fixed
            self.transitions.nodes[comp]["global"] = global_rate

    def add_death_rate(
        self, 
        comps: Union[List, None] = None, 
        rate: Union[None, float] = None, 
        **rates
    ) -> None:
        """
        Add a birth rate to one or more compartments

        Parameters:
        - rate: float
            Death rate
        - comps: list, optional, default=None
            List of compartments to which to assign this death rate.
            If None, apply to all compartments
        """
        self.demographics = True

        if rate is not None:
            count = len(self.params) + 1
            rate_key = "rate" + str(count)
            self.params[rate_key] = rate
        else:
            self.params.update(rates)
            rate_key = list(rates.keys())[0]

        if comps is None:
            comps = self.transitions.nodes()

        for comp in comps:
            if comp not in self.transitions.nodes:
                self.transitions.add_node(comp)

            self.transitions.nodes[comp]["death"] = rate_key

    def add_vaccination(
        self,
        source: str,
        target: str,
        start: int,
        rate: Union[None, float, str],
        **rates,
    ) -> None:
        """
        Add a vaccination transition between two compartments

        Parameters:
        - source: string
            Name of the source compartment
        - target: string
            Name of the target compartment
        - rate: float
            Rate of the vaccination
        - start: int
            Start time of the vaccination

        Returns:
        None
        """

        if rate is not None:
            count = len(self.params) + 1
            rate_key = "rate" + str(count)
            self.params[rate_key] = rate
        else:
            self.params.update(rates)
            rate_key = list(rates.keys())[0]

        self.transitions.add_edge(source, target, rate=rate_key, start=start)

    def add_age_structure(self, matrix: List, population: List) -> List[List]:
        """
        Add a vaccination transition between two compartments

        Parameters:
        - matrix: List
        - population: List
        """
        self.contact = np.asarray(matrix)
        self.population = np.asarray(population).flatten()

        assert self.contact.shape[0] == self.contact.shape[1], (
            "The contact matrix must be square."
        )

        age_groups = list(string.ascii_lowercase[: len(matrix)])
        n_ages = len(age_groups)

        model = EpiModel()
        self.orig_comps = list(self.transitions.nodes())

        for node_i, node_j, data in self.transitions.edges(data=True):
            # Interacting transition
            if "agent" in data:
                for i, age_i in enumerate(age_groups):
                    node_age_i = node_i + "_" + age_i
                    node_age_j = node_j + "_" + age_i

                    for j, age_j in enumerate(age_groups):
                        agent_age = data["agent"] + "_" + age_j

                        model.add_interaction(
                            node_age_i,
                            node_age_j,
                            agent_age,
                            data["rate"] * self.contact[i][j],
                        )

            # Spontaneous transition
            else:
                for age_i in age_groups:
                    node_age_i = node_i + "_" + age_i
                    node_age_j = node_j + "_" + age_i

                    if "start" not in data:
                        model.add_spontaneous(node_age_i, node_age_j, data["rate"])
                    else:
                        # vaccination
                        model.add_vaccination(
                            node_age_i, node_age_j, data["rate"], data["start"]
                        )

        self.transitions = model.transitions

    def _new_cases(self, time: float, population: np.ndarray, pos: Dict) -> np.ndarray:
        """
        Internal function used by integration routine

        Parameters:
        - time: float
            Current time
        - population: numpy array
            Current population of each compartment
        - pos: dict
            Dictionary mapping compartment names to indices

        Returns:
        numpy array
            Array of new cases for each compartment
        """
        diff = np.zeros(len(pos))
        N = np.sum(population)

        if self.population is not None:
            N = {}

            for comp_i in self.transitions.nodes():
                age_group = comp_i.split("_")[-1]

                for comp_j in pos:
                    if comp_j.endswith(age_group):
                        N[comp_i] = N.get(comp_i, 0) + population[pos[comp_j]]

        for edge in self.transitions.edges(data=True):
            source = edge[0]
            target = edge[1]
            trans = edge[2]

            rate_val = self.params[trans["rate"]]
            rate = rate_val * population[pos[source]]

            if "start" in trans and trans["start"] >= time:
                continue

            if "agent" in trans:
                agent = trans["agent"]
                rate *= population[pos[agent]]

                if trans["norm"]:
                    if self.population is None:
                        rate /=  N
                    else:
                        rate /= N[agent]
                
                if self.seasonality is not None:
                    curr_t = int(time) % 365
                    season = float(self.seasonality[curr_t])
                    rate *= season

            if "viral_source" not in trans or trans["viral_source"]:
                diff[pos[source]] -= rate
            # Make sure viral generations are asymetric
            if "viral_target" not in trans or trans["viral_target"]:
                diff[pos[target]] += rate

        # Population dynamics
        if self.demographics:
            for comp, data in self.transitions.nodes(data=True):
                comp_id = pos[comp]

                if "birth" in data:
                    if "fixed" in data and data["fixed"]:
                        births = self.params[data["birth"]]
                    else:
                        if data["global"]:
                            total_population = population.sum()
                            births = total_population * self.params[data["birth"]]
                        else:
                            births = population[comp_id] * self.params[data["birth"]]

                    diff[comp_id] += births

                if "death" in data:
                    deaths = population[comp_id] * self.params[data["death"]]
                    diff[comp_id] -= deaths

        return diff

    def plot(
        self,
        title: Union[str, None] = None,
        normed: bool = True,
        show: bool = True,
        ax: Union[plt.Axes, None] = None,
        **kwargs,
    ):
        """
        Convenience function for plotting

        Parameters:
        - title: string, optional, default=None
            Title of the plot
        - normed: bool, default=True
            Whether to normalize the values or not
        - ax: matplotlib Axes object, default=None
            The Axes object to plot to. If None, a new figure is created.
        - show: bool, default=True
            Whether to call plt.show() or not
        - kwargs: keyword arguments
            Additional arguments to pass to the plot function

        Returns:
        matplotlib.axes._subplots.AxesSubplot
            The plot object
        """
        try:
            if normed:
                N = self.values_.iloc[0].sum()
            else:
                N = 1

            if ax is None:
                ax = plt.gca()

            for comp in self.values_.columns:
                (self.values_[comp] / N).plot(c=utils.EPI_COLORS[comp[0]], **kwargs)

            ax.legend(self.values_.columns)
            ax.set_xlabel("Time")
            ax.set_ylabel("Population")

            if title is not None:
                ax.set_title(title)

            if show:
                plt.show()

            return ax
        except Exception as e:
            print(e)
            raise utils.NotInitialized("You must call integrate() or simulate() first")

    def __getattr__(self, name: str) -> pd.Series:
        """
        Dynamic method to return the individual compartment values

        Parameters:
        - name: string
            Name of the compartment

        Returns:
        pandas.Series
            The values of the specified compartment
        """
        if "values_" in self.__dict__:
            return self.values_[name]
        else:
            raise AttributeError("'EpiModel' object has no attribute '%s'" % name)

    def simulate(
        self,
        timesteps: int,
        t_min: int = 1,
        seasonality: Union[np.ndarray, None] = None,
        **kwargs,
    ) -> None:
        """
        Stochastically simulate the epidemic model

        Parameters:
        - timesteps: int
            Number of time steps to simulate
        - t_min: int, optional
            Starting time
        - seasonality: numpy array, optional
            Array of seasonal factors
        - kwargs: keyword arguments
            Initial population of each compartment

        Returns:
        None
        """
        pos = {comp: i for i, comp in enumerate(self.transitions.nodes())}
        population = np.zeros(len(pos), dtype="int")

        for comp in kwargs:
            population[pos[comp]] = kwargs[comp]

        values = []
        values.append(population)

        comps = list(self.transitions.nodes)
        time = np.arange(t_min, t_min + timesteps, 1, dtype="int")

        self.seasonality = seasonality

        for t in time:
            pop = values[-1]
            new_pop = values[-1].copy()
            N = np.sum(pop)

            for comp in comps:
                trans = list(self.transitions.edges(comp, data=True))

                prob = np.zeros(len(comps), dtype="float")

                for _, node_j, data in trans:
                    source = pos[comp]
                    target = pos[node_j]

                    rate = self.params[data["rate"]]  # self.params[data["rate"]]

                    if "start" in data and data["start"] >= t:
                        continue

                    if "agent" in data:
                        agent = pos[data["agent"]]
                        rate *= pop[agent] / N

                        if self.seasonality is not None:
                            curr_t = int(t) % 365
                            season = float(self.seasonality[curr_t])
                            rate *= season

                    prob[target] = rate

                prob[source] = 1 - np.sum(prob)

                delta = self.rng.multinomial(pop[source], prob)
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
                        births = self.rng.binomial(
                            pop[comp_id], self.params[data["birth"]]
                        )
                        new_pop[comp_id] += births

                    if "death" in data:
                        deaths = self.rng.binomial(
                            pop[comp_id], self.params[data["death"]]
                        )
                        new_pop[comp_id] -= deaths

            values.append(new_pop)

        values = np.array(values)
        self.values_ = pd.DataFrame(values[1:], columns=comps, index=time)

    def integrate(
        self,
        timesteps: int,
        t_min: int = 1,
        seasonality: Union[np.ndarray, None] = None,
        **kwargs,
    ) -> None:
        """
        Numerically integrate the epidemic model

        Parameters:
        - timesteps: int
            Number of time steps to integrate
        - t_min: int, optional
            Starting time
        - seasonality: numpy array, optional
            Array of seasonality values
        - kwargs: keyword arguments
            Initial population of each compartment

        Returns:
        None
        """
        pos = {comp: i for i, comp in enumerate(self.transitions.nodes())}
        population = np.zeros(len(pos))

        for comp in kwargs:
            if self.population is None:
                if comp not in pos:
                    continue

                population[pos[comp]] = kwargs[comp]
            else:
                total_pop = self.population.sum()
                p = np.copy(self.population) / total_pop
                n = self.rng.multinomial(kwargs[comp], p, 1)[0]

                for i, age in enumerate(string.ascii_lowercase[: len(p)]):
                    comp_age = comp + "_" + age
                    if comp_age not in pos:
                        continue

                    population[pos[comp_age]] = n[i]

        time = np.arange(t_min, t_min + timesteps)

        self.seasonality = seasonality
        values = pd.DataFrame(
            scipy.integrate.solve_ivp(
                fun=self._new_cases,
                t_span=(time[0], time[-1]),
                y0=population,
                t_eval=time,
                args=(pos,),
                method="LSODA",
            ).y.T,
            columns=pos.keys(),
            index=time,
        )

        if self.population is None:
            self.values_ = values
        else:
            self.values_ages_ = values

            totals = values.T.copy()
            totals["key"] = totals.index.map(lambda x: "_".join(x.split("_")[:-1]))
            totals = totals.groupby("key").sum().T
            totals.columns.name = None
            self.values_ = totals[self.orig_comps].copy()

    def single_step(self, seasonality=None, **kwargs):
        if hasattr(self, "values_") is False:
            self.simulate(2, 1, seasonality=seasonality, **kwargs)
        else:
            old_values = self.values_.copy()
            t_curr = self.values_.index.max()
            self.simulate(2, t_curr, seasonality=seasonality, **kwargs)
            new_values = pd.concat([old_values, self.values_.iloc[[-1]]])
            self.values_ = new_values

    def __repr__(self) -> str:
        """
        Return a string representation of the EpiModel object

        Returns:
        string
            String representation of the EpiModel object
        """
        text = "# Epidemic Model with %u compartments and %u transitions:" % (
            self.transitions.number_of_nodes(),
            self.transitions.number_of_edges(),
        )
        if self.name is not None:
            text += "\nName: %s\n\n" % self.name
        else:
            text += "\n\n"

        text += "Parameters:\n"
        for rate, value in self.params.items():
            text += "  %s : %s\n" % (rate, value)
        text += "\n\nTransitions:\n"

        for edge in self.transitions.edges(data=True):
            source = edge[0]
            target = edge[1]
            trans = edge[2]

            # Interaction
            if "agent" in trans:
                agent = trans["agent"]
                text += "  - %s + %s = %s %s\n" % (source, agent, target, trans["rate"])
            # Vaccination
            elif "start" in trans:
                start = trans["start"]
                text += "  - %s -> %s %s starting at %s days\n" % (
                    source,
                    target,
                    rate,
                    start,
                )
            # Viral transition
            elif "source_rate" in trans:
                text += "  - %s => %s %s %s" % (
                    source,
                    target,
                    trans["source_rate"],
                    trans["target_rate"],
                )
            # Spontaneous
            else:
                text += "  - %s -> %s %s\n" % (source, target, rate)

        if self.demographics:
            text += "\n\nDemographics:\n"
            for comp, data in self.transitions.nodes(data=True):
                if "birth" in data:
                    text += "  - -> %s: %s # birth rate\n" % (comp, data["birth"])
                if "death" in data:
                    text += "  - %s ->: %s # death rate\n" % (comp, data["death"])

        R0 = self.R0()

        if R0 is not None:
            text += "\n# R0=%1.2f" % R0

        return text

    def _get_active(self) -> Set:
        active = set()

        for node_i, node_j, data in self.transitions.edges(data=True):
            if "agent" in data:
                active.add(data["agent"])
            else:
                active.add(node_i)

        return active

    def _get_susceptible(self) -> Set:
        susceptible = set(
            [node for node, deg in self.transitions.in_degree() if deg == 0]
        )

        if len(susceptible) == 0:
            for node_i, node_j, data in self.transitions.edges(data=True):
                if "agent" in data:
                    susceptible.add(node_i)

        return susceptible

    def _get_infections(self) -> Dict:
        inf = {}

        for node_i, node_j, data in self.transitions.edges(data=True):
            if "agent" in data:
                agent = data["agent"]

                if agent not in inf:
                    inf[agent] = {}

                if node_i not in inf[agent]:
                    inf[agent][node_i] = {}

                inf[agent][node_i]["target"] = node_j
                inf[agent][node_i]["rate"] = data["rate"]

        return inf

    def draw_model(self, ax: Union[plt.Axes, None] = None, show: bool = True) -> None:
        """
        Plot the model structure

        - ax: matplotlib Axes object, default=None
            The Axes object to plot to. If None, a new figure is created.
        - show: bool, default=True
            Whether to call plt.show() or not
        """

        trans = self.transitions.copy()

        try:
            from networkx.drawing.nx_agraph import graphviz_layout

            pos = graphviz_layout(trans, prog="dot", args='-Grankdir="LR"')
        except:
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
                node_colors.append("white")
                orig_pos = pos[node[3:]]
                pos[node] = [orig_pos[0], orig_pos[1] + 1]
            elif node.startswith("_d_"):
                node_colors.append("white")
                orig_pos = pos[node[3:]]
                pos[node] = [orig_pos[0], orig_pos[1] - 1]
            else:
                node_colors.append(utils.EPI_COLORS[node[0]])

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

        nx.draw(
            trans,
            pos,
            with_labels=False,
            arrows=True,
            node_shape="H",
            node_color=node_colors,
            node_size=1000,
            ax=ax,
        )
        nx.draw_networkx_edge_labels(trans, pos, edge_labels=edge_labels, ax=ax)
        nx.draw_networkx_labels(self.transitions, pos2, font_color="k")

        if show:
            plt.show()

    def R0(self) -> Union[float, None]:
        """
        Return the value of the basic reproductive ratio, $R_0$, for the model as defined

        The calculation is completely generic as it uses the Next-Generation matrix approach
        defined in J. R. Soc Interface 7, 873 (2010)

        Returns:
        R0 - the value of the largest eigenvalue of the next generation matrix
        """

        infected = set()

        susceptible = self._get_susceptible()

        for node_i, node_j, data in self.transitions.edges(data=True):
            if "agent" in data:
                infected.add(data["agent"])
                infected.add(node_j)

        infected = sorted(infected)
        N_infected = len(infected)

        F = np.zeros((N_infected, N_infected), dtype="float")
        V = np.zeros((N_infected, N_infected), dtype="float")

        pos = dict(zip(infected, np.arange(N_infected)))

        try:
            for node_i, node_j, data in self.transitions.edges(data=True):
                rate = self.params[data["rate"]]  # self.params[data["rate"]]

                if "agent" in data:
                    target = pos[node_j]
                    agent = pos[data["agent"]]

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

    def save_model(self, filename: str) -> None:
        """
        Save the model to a file

        Parameters:
        - filename: string
            Name of the file to save the model to

        Returns:
        None
        """
        with open(filename, "wt") as f:
            f.write(self.__repr__())

    def list_models() -> List[str]:
        """
        List the models available in the official repository
        """
        remote_path = utils.get_remote_path()
        remote_path = os.path.join(remote_path, "model_list.txt")

        cache_dir = utils.get_cache_directory()

        if not os.path.exists(cache_dir):
            os.makedirs(cache_dir)

        local_path = os.path.join(cache_dir, "model_list.txt")

        # Always get the latest version
        urlretrieve(remote_path, local_path)
        models = [line.strip() for line in open(local_path, "rt").readlines()]

        # Add any user models
        if os.path.exists(cache_dir):
            for file in os.listdir(cache_dir):
                if file.endswith(".yaml"):
                    models.append(file)

        # Make sure to remove duplicates
        return sorted(set(models))

    def load_model(filename: str) -> None:
        """
        Load the model from a file

        Parameters:
        - filename: string
            Name of the file to load the model from

        Returns:
        None
        """
        data = yaml.load(open(filename, "rt"), Loader=yaml.FullLoader)
        model = EpiModel()

        i_re = r"([a-zA-Z]+)\s*\+\s*([a-zA-Z]+)\s*=\s*([a-zA-Z]+)\s*([a-zA-Z]+)"
        s_re = r"([a-zA-Z]+)\s*->\s*([a-zA-Z]+)\s+([a-zA-Z]+)$"

        for key in data:
            if key.lower() == "parameters":
                params = data[key]
            elif key.lower() == "transitions":
                for trans in data[key]:
                    if re.match(i_re, trans):
                        source, agent, target, rate = re.search(i_re, trans).groups()
                        args = {rate: params[rate]}
                        model.add_interaction(source, target, agent, **args)
                    elif re.match(s_re, trans):
                        source, target, rate = re.search(s_re, trans).groups()
                        args = {rate: params[rate]}
                        model.add_spontaneous(source, target, **args)
            elif key.lower() == "name":
                model.name = data[key]

        return model

    def download_model(
        filename: str, repo: Union[str, None] = None, load_model: bool = True
    ) -> Union[None, Self]:
        """
        Download model from offical repository
        """
        remote_path = utils.get_remote_path(repo)
        remote_path = os.path.join(remote_path, filename)

        cache_dir = utils.get_cache_directory()

        if not os.path.exists(cache_dir):
            os.makedirs(cache_dir)

        local_path = os.path.join(cache_dir, filename)

        if not os.path.exists(local_path):
            urlretrieve(remote_path, local_path)

        if load_model:
            return EpiModel.load_model(local_path)
