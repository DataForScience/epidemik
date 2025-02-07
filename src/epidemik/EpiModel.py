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
    """Simple Epidemic Model Implementation
    
        Provides a way to implement and numerically integrate 
    """
    def __init__(self, compartments=None):
        """
        Initialize the EpiModel object
        
        Parameters:
        - compartments: list of strings, optional
            List of compartment names
        
        Returns:
        None
        """
        self.transitions = nx.MultiDiGraph()
        self.seasonality = None
        self.population = None
        self.orig_comps = None
        self.demographics = False
        
        if compartments is not None:
            self.transitions.add_nodes_from([comp for comp in compartments])
    
    def add_interaction(self, source: str, target: str, agent: str, rate: float) -> None:  
        """
        Add an interaction between two compartments
        
        Parameters:
        - source: string
            Name of the source compartment
        - target: string
            Name of the target compartment
        - agent: string
            Name of the agent
        - rate: float
            Rate of the interaction
        
        Returns:
        None
        """      
        self.transitions.add_edge(source, target, agent=agent, rate=rate)        
        
    def add_spontaneous(self, source: str, target: str, rate: float) -> None:
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
        self.transitions.add_edge(source, target, rate=rate)

    def add_birth_rate(self, rate: float, comps: Union[List, None] = None) -> None:
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

        if comps is None:
            comps = self.transitions.nodes()

        for comp in comps:
            self.transitions.nodes[comp]['birth']=rate

    def add_death_rate(self, rate: float, comps: Union[List, None] = None) -> None:
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

        if comps is None:
            comps = self.transitions.nodes()

        for comp in comps:
            self.transitions.nodes[comp]['death']=rate

    def add_vaccination(self, source: str, target: str, rate: float, start: int) -> None:
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
        self.transitions.add_edge(source, target, rate=rate, start=start)

    def add_age_structure(self, matrix: List, population: List) -> List[List]:
        """ 
        Add a vaccination transition between two compartments
        
        Parameters:
        - matrix: List
        - population: List
        """
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
        
    def _new_cases(self, time: float, population: np.ndarray,  pos: Dict) -> np.ndarray:
        """
        Internal function used by integration routine
        
        Parameters:
        - population: numpy array
            Current population of each compartment
        - time: float
            Current time
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
                age_group = comp_i.split('_')[-1]

                for comp_j in pos:
                    if comp_j.endswith(age_group):
                        N[comp_i] = N.get(comp_i, 0) + population[pos[comp_j]]
        
        for edge in self.transitions.edges(data=True):
            source = edge[0]
            target = edge[1]
            trans = edge[2]
            
            rate = trans['rate']*population[pos[source]]
            
            if 'start' in trans and trans['start'] >= time:
                continue

            if 'agent' in trans:
                agent = trans['agent']

                if self.population is None:
                    rate *= population[pos[agent]]/N
                else:
                    rate *= population[pos[agent]]/N[agent]

                if self.seasonality is not None:
                    curr_t = int(time)%365
                    season = float(self.seasonality[curr_t])
                    rate *= season
                
            diff[pos[source]] -= rate
            diff[pos[target]] += rate
            
            # Population dynamics
            if self.demographics:
                for comp, data in self.transitions.nodes(data=True):
                    comp_id = pos[comp]

                    if "birth" in data:
                        births = population[comp_id]*data["birth"]
                        diff[comp_id] += births

                    if "death" in data:
                        deaths = population[comp_id]*data["death"]
                        diff[comp_id] -= deaths

        return diff
    
    def plot(self, title: Union[str, None]= None, normed: bool = True, show: bool = True, ax: Union[plt.Axes, None] = None, **kwargs):
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
        """
        Dynamic method to return the individual compartment values
        
        Parameters:
        - name: string
            Name of the compartment
        
        Returns:
        pandas.Series
            The values of the specified compartment
        """        
        if 'values_' in self.__dict__:
            return self.values_[name]
        else:
            raise AttributeError("'EpiModel' object has no attribute '%s'" % name)

    def simulate(self, timesteps: int, t_min: int = 1, seasonality: Union[np.ndarray, None] = None, **kwargs) -> None:
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
            N = np.sum(pop)


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
                        agent = pos[data['agent']]
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
                        births = np.random.binomial(pop[comp_id], data["birth"])
                        new_pop[comp_id] += births

                    if "death" in data:
                        deaths = np.random.binomial(pop[comp_id], data["death"])
                        new_pop[comp_id] -= deaths

            values.append(new_pop)

        values = np.array(values)
        self.values_ = pd.DataFrame(values[1:], columns=comps, index=time)
    
    def integrate(self, timesteps: int , t_min: int = 1, seasonality: Union[np.ndarray, None] = None, **kwargs) -> None:
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
                    scipy.integrate.solve_ivp(
                        fun=self._new_cases,
                        t_span=(time[0], time[-1]), 
                        y0=population,
                        t_eval=time,
                        args=(pos,),
                        method='LSODA',
                    ).y.T, columns=pos.keys(), index=time
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
        if hasattr(self, 'values_') is False:
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
