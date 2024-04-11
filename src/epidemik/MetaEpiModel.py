### −∗− mode : python ; −∗−
# @file MetaEpiModel.py
# @author Bruno Goncalves
######################################################

import networkx as nx
import numpy as np
from numpy import linalg
from numpy import random
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection

from .EpiModel import *

from tqdm import tqdm
tqdm.pandas()

class MetaEpiModel:
    """Simple Epidemic Model Implementation
    
        Provides a way to implement and numerically integrate 
    """
    def __init__(self, travel_graph, populations):
        """
        Initialize the EpiModel object
        
        Parameters:
        - compartments_: list of strings, optional
            List of compartment names
        
        Returns:
        None
        """
        self.travel_graph = travel_graph
        self.populations = populations

        models = {}

        self.transitions = None
        self.prototype = None
        self.seasonality = None

        for state in travel_graph.index:
            models[state] = EpiModel() 
            if self.transitions is None:
                self.transitions = models[state].transitions
                self.prototype = models[state] 

        self.models = models

    def __repr__(self):
        """
        Return a string representation of the EpiModel object
        
        Returns:
        string
            String representation of the EpiModel object
        """

        model_text = self.models[self.travel_graph.index[0]]

        text = "Metapopulation model with %u populations\n\nThe disease is defined by an %s" % (self.travel_graph.shape[0], model_text)
        return text 

    def add_interaction(self, source, target, agent, rate):  
        """
        Add an interaction between two compartments_
        
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
        for state in self.models:  
            self.models[state].add_interaction(source, target, agent, rate)       
        
    def add_spontaneous(self, source, target, rate):
        """
        Add a spontaneous transition between two compartments_
        
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
        for state in self.models:  
            self.models[state].add_spontaneous(source, target, rate)      

    def add_vaccination(self, source, target, rate, start):
        """
        Add a vaccination transition between two compartments_
        
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
        for state in self.models:  
            self.models[state].add_vaccination(source, target, rate, start)
    
    def get_state(self, state):
        """
        Return a reference to a state EpiModel object

        Parameters:
        - state: string
            Name of the state to return
        """

        return self.models[state]

    def _initialize_populations(self, susceptible):
        columns = list(self.transitions.nodes())
        self.compartments_ = pd.DataFrame(np.zeros((self.travel_graph.shape[0], len(columns)), dtype='int'), columns=columns)
        self.compartments_.index = self.populations.index

        susceptible = list(self.prototype._get_susceptible())[0]

        self.compartments_.loc[:, susceptible] = self.populations['Population']

    def _run_travel(self, compartments_, travel):
        def travel_step(x, populations):
            n = populations.loc[x.name]
            p = travel.loc[x.name].values.tolist()
            output = np.random.multinomial(n, p)

            return pd.Series(output, index=travel.columns)
        
        
        new_compartments = compartments_.copy()
        
        # Travel occurs independently for each compartment
        # since we don't allow in-flight transitions
        for comp in compartments_.columns:
            new_compartments[comp] = travel.apply(travel_step, populations=compartments_[comp]).sum(axis=1)
            
        return new_compartments
    
    def _run_spread(self, models, compartments_, seasonality):
        for state in self.compartments_.index:
            pop = dict(self.compartments_.loc[state].to_dict())
            self.models[state].single_step(**pop)
            self.compartments_.loc[state] = self.models[state].values_.iloc[[-1]].values[0]

    def simulate(self, timestamp, t_min=1, seasonality=None, seed_state=None, susceptible='S', **kwargs):
        if seed_state is None:
            raise NotInitialized("You have to specify the seed_state")

        self._initialize_populations(susceptible)

        pos = {comp: i for i, comp in enumerate(self.models[seed_state].transitions.nodes())}

        for comp in kwargs:
            if comp != susceptible and comp in pos:
                self.compartments_.loc[seed_state, comp] += kwargs[comp]
                self.compartments_.loc[seed_state, susceptible] -= kwargs[comp]

        for t in tqdm(range(t_min, timestamp+1), total=timestamp):
            self._run_spread(self.models, self.compartments_, self.seasonality)
            self.compartments_ = self._run_travel(self.compartments_, self.travel_graph)
    
    def plot(self, title=None, normed=True, layout=None, **kwargs):
        if layout is None:
            n_pop = self.travel_graph.shape[0]
            N = int(np.round(np.sqrt(n_pop), 0)+1)

            coords = [[i%N, i//N] for i in range(n_pop)]
            coords = pd.DataFrame(coords, columns=['x', 'y'])
            coords.index = self.travel_graph.index
        else:
            coords = layout

        fig, ax = plt.subplots(1, figsize=(16, 22))
        ax.set_aspect(1.)
        ax.invert_yaxis()

        patches = []
        color_list = []
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        x = np.linspace(0., 0.75, self.models['NY'].values_.shape[0])



        for state in self.models:
            daily = self.models[state]['I'].values
            
            timeline = (daily/daily.max())
            color = colors[3]
            
            peak = self.models[state]['I'][self.models[state]['I']>0].index.min()
            
            if peak <= 10:
                color = colors[0]
            elif peak <= 20:
                color = colors[1]
            elif peak <= 30:
                color = colors[2]
            elif peak <= 40:
                color = colors[3]
            elif peak > 40:
                color = colors[4]
                
            fancybox = mpatches.FancyBboxPatch([coords.x[state]-0.5, coords.y[state]-0.5], 0.8, 0.8,
                                                boxstyle=mpatches.BoxStyle("Round", pad=0.06))
            patches.append(fancybox)
            color_list.append(color)
            
            ax.text(y=coords.y[state]-0.4, x=coords.x[state]-0.4, 
                    s=state, horizontalalignment='center', verticalalignment='center', fontsize=15)
            ax.plot(coords.x[state]+x-0.5, coords.y[state]-timeline/1.4+0.25, lw=1, color=colors[1])
            
            ax.vlines(x=coords.x[state]+x[peak]-0.5, 
                    ymin=coords.y[state]-timeline.min()/1.4+0.25, 
                    ymax=coords.y[state]-timeline.max()/1.4+0.25, lw=1,
                    color='darkgray')

        collection = PatchCollection(patches, facecolors=color_list, alpha=0.3)
        ax.add_collection(collection)

        x_min, x_max = ax.get_xlim()
        y_min, y_max = ax.get_ylim()

        x_delta = (x_max-x_min+1)/7

        ax.text(x=1, y=y_min, s='[0-10]',  color=colors[0], fontsize=15, ha='center')
        ax.text(x=2, y=y_min, s='(10-20]', color=colors[1], fontsize=15, ha='center')
        ax.text(x=3, y=y_min, s='(20-30]', color=colors[2], fontsize=15, ha='center')
        ax.text(x=4, y=y_min, s='(30-40]', color=colors[3], fontsize=15, ha='center')
        ax.text(x=5, y=y_min, s='40+',     color=colors[4], fontsize=15, ha='center')


        ax.axis('off')
        fig.patch.set_facecolor('#FFFFFF')
        fig.tight_layout()

    def plot_peaks(self):
        peaks = None

        for state in self.models.values():
            if peaks is None:
                peaks = state.values_[['I']].copy()
            else:
                peaks = pd.concat([peaks, state.values_[['I']]], axis=1)

        peaks = peaks.T
        peaks.index = self.travel_graph.index
        peaks = peaks.apply(lambda x:x/x.max(), axis=1)

        n, m = peaks.shape
        ratio = int(np.round(m/n, 0)+1)
        fig, ax = plt.subplots(1, figsize=(15, 10*ratio))
        import matplotlib as mpl

        norm = mpl.colors.Normalize(vmin=0,vmax=1)
        sm = plt.cm.ScalarMappable(cmap=plt.cm.rainbow, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, orientation='horizontal')
        cbar.ax.tick_params(labelsize=10)
        cbar.set_label('Peak Fraction', fontsize=14)

        ax.imshow(peaks.loc[self.travel_graph.index].fillna(0).values, cmap=plt.cm.rainbow)
        ax.set_yticks(np.arange(0, len(self.travel_graph.index)))
        ax.set_yticklabels(self.travel_graph.index, fontsize=10)
        ax.set_xticks(np.arange(0, peaks.shape[1], 3))
        ax.set_xticklabels(np.arange(0, peaks.shape[1], 3), fontsize=10)
        # ax.set_aspect(1)
        fig.patch.set_facecolor('#FFFFFF')
    
if __name__ == '__main__':

    Nk_uk = pd.read_csv("data/United Kingdom-2020.csv", index_col=0)
    Nk_ke = pd.read_csv("data/Kenya-2020.csv", index_col=0)

    contacts_uk = pd.read_excel("data/MUestimates_all_locations_2.xlsx", sheet_name="United Kingdom of Great Britain", header=None)
    contacts_ke = pd.read_excel("data/MUestimates_all_locations_1.xlsx", sheet_name="Kenya")

    beta = 0.05
    mu = 0.1

    SIR_uk = EpiModel()
    SIR_uk.add_interaction('S', 'I', 'I', beta)
    SIR_uk.add_spontaneous('I', 'R', mu)


    SIR_ke = EpiModel()
    SIR_ke.add_interaction('S', 'I', 'I', beta)
    SIR_ke.add_spontaneous('I', 'R', mu)

    N_uk = int(Nk_uk.sum())
    N_ke = int(Nk_ke.sum())


    SIR_uk.add_age_structure(contacts_uk, Nk_uk)
    SIR_ke.add_age_structure(contacts_ke, Nk_ke)

    SIR_uk.integrate(100, S=N_uk*.99, I=N_uk*.01, R=0)
    SIR_ke.integrate(100, S=N_ke*.99, I=N_ke*.01, R=0)

    fig, ax = plt.subplots(1)
    SIR_uk.draw_model(ax)
    fig.savefig('SIR_model.png', dpi=300, facecolor='white')

    fig, ax = plt.subplots(1)

    (SIR_uk['I']*100/N_uk).plot(ax=ax)
    (SIR_ke['I']*100/N_ke).plot(ax=ax)
    ax.legend(['UK', 'Kenya'])
    ax.set_xlabel('Time')
    ax.set_ylabel('Population (%)')

    fig.savefig('SIR_age.png', dpi=300, facecolor='white')