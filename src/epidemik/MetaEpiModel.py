### −∗− mode : python ; −∗−
# @file MetaEpiModel.py
# @author Bruno Goncalves
######################################################

from typing import Dict, List, Union
import warnings

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
from tqdm import tqdm

from .EpiModel import EpiModel
from .utils import NotInitialized


class MetaEpiModel:
    """Metapopulation compartmental epidemic model.

    A :class:`MetaEpiModel` is a collection of independent
    :class:`~epidemik.EpiModel` instances -- one per sub-population --
    coupled together by a travel matrix. At every simulation step each
    sub-population is advanced by one stochastic time step, after which
    individuals move between sub-populations according to the travel
    probabilities. The same disease model is shared across every
    sub-population.

    Attributes
    ----------
    travel_graph : pandas.DataFrame
        Row-stochastic travel matrix, indexed and columned by
        sub-population label.
    populations : pandas.DataFrame
        Data frame with the total population of each sub-population.
    models : dict
        Mapping from sub-population label to its :class:`~epidemik.EpiModel`.
    compartments_ : pandas.DataFrame
        Current number of individuals in each compartment (columns) of
        each sub-population (rows). Available after :meth:`simulate`.
    """

    def __init__(self, travel_graph: pd.DataFrame, populations: pd.DataFrame, population: str = 'Population'):
        """Initialize a metapopulation model.

        Parameters
        ----------
        travel_graph : pandas.DataFrame
            Square data frame indexed and columned by sub-population labels.
            Entry ``(i, j)`` is the probability that an individual currently
            in sub-population ``i`` travels to sub-population ``j`` during
            a single time step. Each row must sum to 1.
        populations : pandas.DataFrame
            Data frame indexed by sub-population label containing at least
            one column with the total population of each sub-population.
        population : str, default='Population'
            Name of the column in ``populations`` that holds the population
            counts.

        Raises
        ------
        ValueError
            If ``travel_graph`` is not square, if its rows and columns are
            not labelled by the same sub-populations, if its rows do not sum
            to 1, or if ``populations`` is missing a sub-population or the
            ``population`` column.
        """
        if travel_graph.shape[0] != travel_graph.shape[1] or set(travel_graph.index) != set(travel_graph.columns):
            raise ValueError(
                "travel_graph must be a square data frame whose rows and "
                "columns are labelled by the same sub-populations."
            )

        # Make sure columns are in the same order as the rows
        travel_graph = travel_graph.loc[travel_graph.index, travel_graph.index].astype(float)

        row_sums = travel_graph.sum(axis=1)

        if not np.allclose(row_sums, 1.0, atol=1e-6):
            bad = row_sums[~np.isclose(row_sums, 1.0, atol=1e-6)].index.tolist()
            raise ValueError(
                "Every row of travel_graph must sum to 1 (offending rows: %s)." % bad
            )

        missing = set(travel_graph.index) - set(populations.index)

        if missing:
            raise ValueError(
                "populations is missing the following sub-populations: %s" % sorted(missing)
            )

        if population not in populations.columns:
            raise ValueError("populations has no column named '%s'." % population)

        self.travel_graph = travel_graph
        self.populations = populations
        self.population = population
        self.seasonality = None

        models = {}

        for state in travel_graph.index:
            models[state] = EpiModel()

        self.models = models
        self.prototype = models[travel_graph.index[0]]
        self.transitions = self.prototype.transitions

    def __repr__(self) -> str:
        """Return a human-readable description of the metapopulation model.

        Returns
        -------
        str
            Multi-line string describing the number of sub-populations and
            the underlying disease model.
        """
        text = "Metapopulation model with %u populations\n\nThe disease is defined by an %s" % (
            self.travel_graph.shape[0], self.prototype
        )
        return text

    # ------------------------------------------------------------------
    # Model definition -- every call is forwarded to each sub-population
    # ------------------------------------------------------------------
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
        for model in self.models.values():
            model.add_interaction(source, target, agent, rate)

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
        for model in self.models.values():
            model.add_spontaneous(source, target, rate)

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
        for model in self.models.values():
            model.add_vaccination(source, target, rate, start)

    def add_birth_rate(self, rate: float, comps: Union[List, None] = None, from_group: bool = False) -> None:
        """Attach a birth rate to one or more compartments in every sub-population.

        See :meth:`EpiModel.add_birth_rate` for the semantics.

        Parameters
        ----------
        rate : float
            Birth rate.
        comps : list, optional
            Compartments to which to assign this birth rate. If ``None``
            (the default) the rate is applied to every compartment.
        from_group : bool, default=False
            Whether the inflow is proportional to the group total rather
            than to the compartment itself (requires :meth:`add_groups`).
        """
        for model in self.models.values():
            model.add_birth_rate(rate, comps=comps, from_group=from_group)

    def add_death_rate(self, rate: float, comps: Union[List, None] = None) -> None:
        """Attach a death rate to one or more compartments in every sub-population.

        See :meth:`EpiModel.add_death_rate` for the semantics.

        Parameters
        ----------
        rate : float
            Per-capita death rate.
        comps : list, optional
            Compartments to which to assign this death rate. If ``None``
            (the default) the rate is applied to every compartment.
        """
        for model in self.models.values():
            model.add_death_rate(rate, comps=comps)

    def add_groups(self, groups: Dict[str, List[str]]) -> None:
        """Declare compartment groups (e.g. host/vector) in every sub-population.

        Forwards to :meth:`EpiModel.add_groups` on every underlying
        :class:`~epidemik.EpiModel`, so the same multi-group structure is
        used everywhere. See :meth:`EpiModel.add_groups` for the semantics.

        Note that :meth:`simulate` only seeds the *susceptible* compartment
        of every sub-population from ``populations``; the initial size of
        the other groups can only be set in ``seed_state`` through the
        keyword arguments of :meth:`simulate`.

        Parameters
        ----------
        groups : dict of str -> list of str
            Mapping from group name to the list of compartments in that
            group.
        """
        for model in self.models.values():
            model.add_groups(groups)

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
        return self.prototype.R0()

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

    # ------------------------------------------------------------------
    # Simulation
    # ------------------------------------------------------------------
    def _initialize_populations(self, susceptible: str, population: Union[str, None] = None) -> None:
        """Place the whole population of every sub-population in ``susceptible``."""
        columns = list(self.transitions.nodes())
        states = self.travel_graph.index

        if population is None:
            population = self.population

        self.compartments_ = pd.DataFrame(0, index=states, columns=columns, dtype='int')
        self.compartments_[susceptible] = self.populations.loc[states, population].astype('int').values

    def _run_travel(self, compartments_: pd.DataFrame, travel: pd.DataFrame) -> pd.DataFrame:
        """Redistribute individuals across sub-populations.

        Travel occurs independently for each compartment, since in-flight
        transitions are not allowed. Each row of ``travel`` is used as the
        probability vector of a multinomial draw.
        """
        states = list(travel.index)
        probabilities = travel.loc[states, states].values

        new_compartments = pd.DataFrame(0, index=states, columns=compartments_.columns, dtype='int')

        for comp in compartments_.columns:
            counts = compartments_.loc[states, comp].values.astype('int')
            arrivals = np.zeros(len(states), dtype='int')

            for i in range(len(states)):
                if counts[i] > 0:
                    arrivals += np.random.multinomial(counts[i], probabilities[i])

            new_compartments[comp] = arrivals

        return new_compartments

    def _run_spread(self, t_min: int = 1) -> None:
        """Advance every sub-population by one stochastic time step."""
        for state in self.compartments_.index:
            pop = self.compartments_.loc[state].to_dict()
            self.models[state].single_step(seasonality=self.seasonality, t_min=t_min, **pop)
            self.compartments_.loc[state] = self.models[state].values_.iloc[-1].values

    def simulate(
            self, timesteps: Union[int, None] = None, t_min: int = 1,
            seasonality: Union[np.ndarray, None] = None, seed_state: Union[str, None] = None,
            susceptible: str = 'S', **kwargs
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
        timesteps : int
            Number of time steps to simulate. (The old spelling
            ``timestamp`` is still accepted but deprecated.)
        t_min : int, default=1
            Starting time of the simulation.
        seasonality : numpy.ndarray, optional
            Array of seasonal forcing factors (typically of length 365),
            forwarded to every sub-population.
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
        ValueError
            If ``seed_state`` is not a known sub-population, or if
            ``susceptible`` is not a compartment of the model.
        """
        if 'timestamp' in kwargs:
            warnings.warn(
                "The 'timestamp' argument of MetaEpiModel.simulate is deprecated; "
                "use 'timesteps' instead.", DeprecationWarning, stacklevel=2
            )
            timestamp = kwargs.pop('timestamp')

            if timesteps is None:
                timesteps = timestamp

        if timesteps is None:
            raise TypeError("simulate() missing required argument: 'timesteps'")

        if timesteps < 1:
            raise ValueError("timesteps must be a positive integer.")

        if seed_state is None:
            raise NotInitialized("You have to specify the seed_state")

        if seed_state not in self.models:
            raise ValueError("Unknown seed_state '%s'." % seed_state)

        if susceptible not in self.transitions.nodes():
            raise ValueError("'%s' is not a compartment of the model." % susceptible)

        self.seasonality = seasonality
        self._initialize_populations(susceptible)

        for comp, value in kwargs.items():
            if comp == susceptible:
                continue

            if comp not in self.compartments_.columns:
                warnings.warn(
                    "Ignoring initial condition for unknown compartment '%s'." % comp,
                    stacklevel=2
                )
                continue

            self.compartments_.loc[seed_state, comp] += value
            self.compartments_.loc[seed_state, susceptible] -= value

        for t in tqdm(range(t_min, t_min+timesteps), total=timesteps):
            self._run_spread(t_min)
            self.compartments_ = self._run_travel(self.compartments_, self.travel_graph)

    def integrate(self, *args, **kwargs):
        """Not supported for metapopulation models.

        Raises
        ------
        NotImplementedError
            Always. Metapopulation models can only be explored via
            :meth:`simulate`.
        """
        raise NotImplementedError("MetaEpiModel doesn't support direct integration of the ODE")

    # ------------------------------------------------------------------
    # Visualisation
    # ------------------------------------------------------------------
    def draw_model(self, ax: Union[plt.Axes, None] = None, show: bool = True) -> None:
        """Draw the underlying disease model.

        Delegates to :meth:`EpiModel.draw_model` on the first
        sub-population's model, since every sub-population shares the same
        disease.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Axes to draw onto. If ``None`` a new figure is created.
        show : bool, default=True
            Whether to call :func:`matplotlib.pyplot.show` after drawing.
        """
        return self.prototype.draw_model(ax=ax, show=show)

    def _check_simulated(self) -> None:
        if 'compartments_' not in self.__dict__ or 'values_' not in self.prototype.__dict__:
            raise NotInitialized('You must call simulate() first')

    def plot(self, title: Union[str, None] = None, normed: bool = True, layout=None, compartment: str = 'I', figsize=(16, 22), **kwargs) -> plt.Axes:
        """Plot the infection time series for every sub-population.

        Each sub-population is rendered as a small panel whose position is
        given by ``layout`` (or by an automatic square layout if
        ``layout`` is omitted). Inside each panel the trajectory of
        ``compartment`` is drawn normalized by its maximum, a vertical line
        marks the time of the local peak, and the panel colour encodes the
        number of time steps between the start of the simulation and that
        peak.

        Parameters
        ----------
        title : str, optional
            Title to set on the resulting axes.
        normed : bool, default=True
            Unused: every curve is always normalized by its own maximum.
            Kept for API symmetry with :meth:`EpiModel.plot`.
        layout : pandas.DataFrame, optional
            Optional layout with columns ``x`` and ``y`` indexed by
            sub-population label, providing the on-screen position of every
            panel. If omitted, an automatic grid layout is generated.
        compartment : str, default='I'
            Compartment whose trajectory is drawn in every panel.
        figsize : tuple, default=(16, 22)
            Size of the figure created for the plot.
        **kwargs
            Reserved for future use.

        Returns
        -------
        matplotlib.axes.Axes
            The axes that were drawn on.

        Raises
        ------
        NotInitialized
            If :meth:`simulate` has not been called yet.
        ValueError
            If ``compartment`` is not a compartment of the model.
        """
        self._check_simulated()

        if compartment not in self.transitions.nodes():
            raise ValueError("'%s' is not a compartment of the model." % compartment)

        if layout is None:
            n_pop = self.travel_graph.shape[0]
            N = int(np.round(np.sqrt(n_pop), 0)+1)

            coords = [[i % N, i//N] for i in range(n_pop)]
            coords = pd.DataFrame(coords, columns=['x', 'y'])
            coords.index = self.travel_graph.index
        else:
            coords = layout

        fig, ax = plt.subplots(1, figsize=figsize)
        ax.set_aspect(1.)
        ax.invert_yaxis()

        patches = []
        color_list = []
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

        n_steps = self.prototype.values_.shape[0]
        x = np.linspace(0., 0.75, n_steps)

        for state, model in self.models.items():
            series = model.values_[compartment]
            daily = series.values.astype(float)
            peak_value = daily.max()

            if peak_value > 0:
                timeline = daily/peak_value
            else:
                timeline = daily

            peak_pos = int(np.argmax(daily))
            peak = series.index[peak_pos] - series.index[0]
            color = self._peak_color(peak, colors)

            fancybox = mpatches.FancyBboxPatch([coords.x[state]-0.5, coords.y[state]-0.5], 0.8, 0.8,
                                               boxstyle=mpatches.BoxStyle("Round", pad=0.06))
            patches.append(fancybox)
            color_list.append(color)

            ax.text(y=coords.y[state]-0.4, x=coords.x[state]-0.4,
                    s=state, horizontalalignment='center', verticalalignment='center', fontsize=15)
            ax.plot(coords.x[state]+x[:len(timeline)]-0.5, coords.y[state]-timeline/1.4+0.25, lw=1, color=colors[1])

            ax.vlines(x=coords.x[state]+x[peak_pos]-0.5,
                      ymin=coords.y[state]-timeline.min()/1.4+0.25,
                      ymax=coords.y[state]-timeline.max()/1.4+0.25, lw=1,
                      color='darkgray')

        collection = PatchCollection(patches, facecolors=color_list, alpha=0.3)
        ax.add_collection(collection)

        x_min, x_max = ax.get_xlim()
        y_min, y_max = ax.get_ylim()

        ax.text(x=1, y=y_min, s='[0-10]',  color=colors[0], fontsize=15, ha='center')
        ax.text(x=2, y=y_min, s='(10-20]', color=colors[1], fontsize=15, ha='center')
        ax.text(x=3, y=y_min, s='(20-30]', color=colors[2], fontsize=15, ha='center')
        ax.text(x=4, y=y_min, s='(30-40]', color=colors[3], fontsize=15, ha='center')
        ax.text(x=5, y=y_min, s='40+',     color=colors[4], fontsize=15, ha='center')

        if title is not None:
            ax.set_title(title)

        ax.axis('off')
        fig.patch.set_facecolor('#FFFFFF')

        with warnings.catch_warnings():
            # tight_layout cannot always honour the fixed aspect ratio
            warnings.simplefilter('ignore', UserWarning)
            fig.tight_layout()

        return ax

    @staticmethod
    def _peak_color(peak: float, colors) -> str:
        """Colour encoding the delay (in time steps) until the local peak."""
        if peak <= 10:
            return colors[0]
        elif peak <= 20:
            return colors[1]
        elif peak <= 30:
            return colors[2]
        elif peak <= 40:
            return colors[3]
        return colors[4]

    def plot_peaks(self, compartment: str = 'I') -> plt.Axes:
        """Render a heat-map of the infection peak across sub-populations.

        For every sub-population, the time series of ``compartment`` is
        normalized by its maximum value and stacked vertically into a
        heat-map. The result highlights how the epidemic peak diffuses
        across sub-populations over time.

        Parameters
        ----------
        compartment : str, default='I'
            Compartment whose trajectory is used.

        Returns
        -------
        matplotlib.axes.Axes
            The axes that were drawn on.

        Raises
        ------
        NotInitialized
            If :meth:`simulate` has not been called yet.
        ValueError
            If ``compartment`` is not a compartment of the model.
        """
        self._check_simulated()

        if compartment not in self.transitions.nodes():
            raise ValueError("'%s' is not a compartment of the model." % compartment)

        peaks = pd.concat(
            [model.values_[[compartment]] for model in self.models.values()], axis=1
        ).T
        peaks.index = self.travel_graph.index
        peaks = peaks.apply(lambda x: x/x.max() if x.max() > 0 else x, axis=1)

        n, m = peaks.shape
        ratio = int(np.round(m/n, 0)+1)
        fig, ax = plt.subplots(1, figsize=(15, 10*ratio))

        norm = mpl.colors.Normalize(vmin=0, vmax=1)
        sm = plt.cm.ScalarMappable(cmap=plt.cm.rainbow, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, orientation='horizontal')
        cbar.ax.tick_params(labelsize=10)
        cbar.set_label('Peak Fraction', fontsize=14)

        ax.imshow(peaks.loc[self.travel_graph.index].fillna(0).values, cmap=plt.cm.rainbow)
        ax.set_yticks(np.arange(0, len(self.travel_graph.index)))
        ax.set_yticklabels(self.travel_graph.index, fontsize=10)
        ax.set_xticks(np.arange(0, peaks.shape[1], 3))
        ax.set_xticklabels(peaks.columns[::3], fontsize=10)
        fig.patch.set_facecolor('#FFFFFF')

        return ax
