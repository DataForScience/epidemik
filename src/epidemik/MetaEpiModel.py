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
    """Metapopulation compartmental epidemic model.

    A :class:`MetaEpiModel` is a collection of independent
    :class:`~epidemik.EpiModel` instances -- one per sub-population --
    coupled together by a travel matrix. At every simulation step each
    sub-population is advanced by one stochastic time step, after which
    individuals move between sub-populations according to the travel
    probabilities. The same disease model is shared across every
    sub-population.
    """

    def __init__(self, travel_graph: pd.DataFrame, populations: pd.DataFrame, population: str ='Population'):
        """Initialize a metapopulation model.

        Parameters
        ----------
        travel_graph : pandas.DataFrame
            Square data frame indexed and columned by sub-population labels.
            Entry ``(i, j)`` is the probability that an individual currently
            in sub-population ``i`` travels to sub-population ``j`` during
            a single time step. Each row should sum to 1.
        populations : pandas.DataFrame
            Data frame indexed by sub-population label containing at least
            one column with the total population of each sub-population.
        population : str, default='Population'
            Name of the column in ``populations`` that holds the population
            counts.
        """
        self.travel_graph = travel_graph
        self.populations = populations
        self.population = population

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

    def __repr__(self) -> str:
        """Return a human-readable description of the metapopulation model.

        Returns
        -------
        str
            Multi-line string describing the number of sub-populations and
            the underlying disease model.
        """

        model_text = self.models[self.travel_graph.index[0]]

        text = "Metapopulation model with %u populations\n\nThe disease is defined by an %s" % (self.travel_graph.shape[0], model_text)
        return text 

    def add_interaction(self, source: str, target: str, agent: str, rate: float) -> None:
        """Add an interaction transition to every sub-population.

        The transition is added to every :class:`~epidemik.EpiModel` managed
        by this metapopulation model -- the disease is assumed to be the
        same everywhere.

        Parameters
        ----------
        source : str
            Name of the source compartment.
        target : str
            Name of the target compartment.
        agent : str
            Name of the agent compartment.
        rate : float
            Per-capita rate at which the interaction takes place.
        """
        for state in self.models:  
            self.models[state].add_interaction(source, target, agent, rate)       

    def add_spontaneous(self, source: str, target: str, rate: float) -> None:
        """Add a spontaneous transition to every sub-population.

        Parameters
        ----------
        source : str
            Name of the source compartment.
        target : str
            Name of the target compartment.
        rate : float
            Per-capita rate at which the transition takes place.
        """
        for state in self.models:  
            self.models[state].add_spontaneous(source, target, rate)      

    def add_groups(self, groups: Dict[str, List[str]]) -> None:
        """Declare compartment groups (e.g. host/vector) in every sub-population.

        Forwards to :meth:`EpiModel.add_groups` on every underlying
        :class:`~epidemik.EpiModel`, so the same multi-group structure is
        used everywhere. See :meth:`EpiModel.add_groups` for the semantics.

        Parameters
        ----------
        groups : dict of str -> list of str
            Mapping from group name to the list of compartments in that
            group.
        """
        for state in self.models:
            self.models[state].add_groups(groups)

    def add_vaccination(self, source: str, target: str, rate: float, start: int) -> None:
        """Add a time-gated vaccination transition to every sub-population.

        Parameters
        ----------
        source : str
            Name of the source compartment.
        target : str
            Name of the target compartment.
        rate : float
            Per-capita vaccination rate.
        start : int
            Time step at which the vaccination campaign begins.
        """
        for state in self.models:  
            self.models[state].add_vaccination(source, target, rate, start)

    def R0(self) -> Union[float, None]:
        """Compute the basic reproduction number of the underlying disease.

        Since every sub-population shares the same disease model, the
        :math:`R_0` of the metapopulation model is simply the :math:`R_0`
        of any individual sub-population.

        Returns
        -------
        float or None
            Basic reproduction number, or ``None`` if it cannot be computed.
        """
        key = list(self.models.keys())[0]
        return self.models[key].R0()

    def get_state(self, state: str) -> EpiModel:
        """Return the :class:`~epidemik.EpiModel` for a given sub-population.

        Parameters
        ----------
        state : str
            Label of the sub-population.

        Returns
        -------
        EpiModel
            The underlying :class:`~epidemik.EpiModel` driving ``state``.
        """
        return self.models[state]

    def _initialize_populations(self, susceptible: str, population: Union[pd.DataFrame, None] =None) -> None:
        columns = list(self.transitions.nodes())
        self.compartments_ = pd.DataFrame(np.zeros((self.travel_graph.shape[0], len(columns)), dtype='int'), columns=columns)
        self.compartments_.index = self.populations.index

        if population is None:
            population = self.population

        for state in self.compartments_.index:
            self.compartments_.loc[state, susceptible]  = self.populations.loc[state, population]

    def _run_travel(self, compartments_: pd.DataFrame, travel: pd.DataFrame) -> pd.DataFrame:
        def travel_step(x, populations: pd.DataFrame) -> pd.Series:
            n = populations.loc[x.name]
            p = travel.loc[x.name].values.tolist()
            output = np.random.multinomial(n, p)

            return pd.Series(output, index=travel.columns)
        
        
        new_compartments = compartments_.copy()
        
        # Travel occurs independently for each compartment
        # since we don't allow in-flight transitions
        for comp in compartments_.columns:
            new_compartments[comp] = travel.apply(
                    travel_step, 
                    populations=compartments_[comp]
            ).sum(axis=1)
            
        return new_compartments
    
    def _run_spread(self) -> None:
        for state in self.compartments_.index:
            pop = self.compartments_.loc[state].to_dict()
            self.models[state].single_step(**pop)
            self.compartments_.loc[state] = self.models[state].values_.iloc[[-1]].values[0]

    def simulate(
            self, timestamp: int, t_min: int = 1, 
            seasonality=None, seed_state: [str, None] = None, 
            susceptible: str ='S', **kwargs
    ) -> None:
        """Run a coupled stochastic simulation across every sub-population.

        Each time step alternates between an *epidemic* phase (every
        sub-population advances by one step of its own
        :class:`~epidemik.EpiModel`) and a *travel* phase (individuals are
        redistributed across sub-populations according to
        :attr:`travel_graph`). The initial epidemic seed is placed in
        ``seed_state``.

        Parameters
        ----------
        timestamp : int
            Number of time steps to simulate.
        t_min : int, default=1
            Starting time of the simulation.
        seasonality : numpy.ndarray, optional
            Optional seasonality vector, forwarded to every sub-population.
        seed_state : str
            Label of the sub-population where the epidemic seeds are placed.
            Must be provided.
        susceptible : str, default='S'
            Name of the compartment representing susceptible individuals.
        **kwargs
            Initial size of each non-susceptible compartment in
            ``seed_state`` (e.g. ``I=10``). The susceptible compartment is
            adjusted automatically so that the total population of
            ``seed_state`` remains unchanged.

        Raises
        ------
        NotInitialized
            If ``seed_state`` is not provided.
        """
        if seed_state is None:
            raise NotInitialized("You have to specify the seed_state")

        self._initialize_populations(susceptible)

        pos = {comp: i for i, comp in enumerate(self.models[seed_state].transitions.nodes())}

        for comp in kwargs:
            if comp != susceptible and comp in pos:
                self.compartments_.loc[seed_state, comp] += kwargs[comp]
                self.compartments_.loc[seed_state, susceptible] -= kwargs[comp]

        for t in tqdm(range(t_min, timestamp+1), total=timestamp):
            self._run_spread()
            self.compartments_ = self._run_travel(self.compartments_, self.travel_graph)
    
    def integrate(self, **kwargs):
        """Not supported for metapopulation models.

        Raises
        ------
        NotImplementedError
            Always. Metapopulation models can only be explored via
            :meth:`simulate`.
        """
        raise NotImplementedError("MetaEpiModel doesn't support direct integration of the ODE")

    def draw_model(self) -> None:
        """Draw the underlying disease model.

        Delegates to :meth:`EpiModel.draw_model` on the first
        sub-population's model, since every sub-population shares the same
        disease.
        """
        return self.models.iloc[0].draw_model()

    def plot(self, title: Union[str, None] = None, normed: bool = True, layout=None, **kwargs) -> None:
        """Plot the infection time series for every sub-population.

        Each sub-population is rendered as a small panel whose position is
        given by ``layout`` (or by an automatic square layout if
        ``layout`` is omitted). The colour of each panel encodes the time
        at which the epidemic peaks locally.

        Parameters
        ----------
        title : str, optional
            Currently unused, kept for API symmetry with
            :meth:`EpiModel.plot`.
        normed : bool, default=True
            Currently unused, kept for API symmetry with
            :meth:`EpiModel.plot`.
        layout : pandas.DataFrame, optional
            Optional layout with columns ``x`` and ``y`` indexed by
            sub-population label, providing the on-screen position of every
            panel. If omitted, an automatic grid layout is generated.
        **kwargs
            Reserved for future use.
        """
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

    def plot_peaks(self) -> None:
        """Render a heat-map of the infection peak across sub-populations.

        For every sub-population, the time series of the ``I`` compartment
        is normalized by its maximum value and stacked vertically into a
        heat-map. The result highlights how the epidemic peak diffuses
        across sub-populations over time.
        """
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
