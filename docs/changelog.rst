Changelog
=========

Notable changes between releases of ``epidemik``. Older releases are listed
on the `GitHub releases page
<https://github.com/DataForScience/epidemik/releases>`_ and on
`PyPI <https://pypi.org/project/epidemik/#history>`_.

0.2.0 (unreleased)
------------------

New features
~~~~~~~~~~~~

* :meth:`EpiModel.add_groups <epidemik.EpiModel.add_groups>` partitions the
  compartments into named groups (e.g. hosts and vectors). Interactions are
  then normalized by the population of the agent's group, which yields the
  standard host-vector dynamics. :meth:`EpiModel.add_birth_rate
  <epidemik.EpiModel.add_birth_rate>` gained a ``from_group`` flag for
  births proportional to a group total.
* :class:`~epidemik.MetaEpiModel` forwards
  :meth:`~epidemik.MetaEpiModel.add_groups`,
  :meth:`~epidemik.MetaEpiModel.add_birth_rate` and
  :meth:`~epidemik.MetaEpiModel.add_death_rate` to every sub-population.
* :meth:`MetaEpiModel.plot <epidemik.MetaEpiModel.plot>` and
  :meth:`MetaEpiModel.plot_peaks <epidemik.MetaEpiModel.plot_peaks>` accept
  a ``compartment`` argument instead of always using ``"I"``; ``plot`` also
  honours ``title`` and accepts ``figsize``.
* :meth:`NetworkEpiModel.simulate <epidemik.NetworkEpiModel.simulate>`
  accepts ``t_min`` and ``susceptible`` arguments and supports networks with
  arbitrary (non-integer) node labels.
* :meth:`EpiModel.single_step <epidemik.EpiModel.single_step>` continues
  from the current state when called without keyword arguments.
* Age-structured models accept initial conditions for individual strata
  (``S_a=...``), and birth/death rates attached before
  :meth:`~epidemik.EpiModel.add_age_structure` are propagated to every
  stratum.

Behaviour changes
~~~~~~~~~~~~~~~~~

* :meth:`EpiModel.simulate <epidemik.EpiModel.simulate>` now stores the
  initial conditions in the first row of ``values_`` (at ``t_min``) and
  produces exactly ``timesteps`` rows, exactly like
  :meth:`~epidemik.EpiModel.integrate`. Previously the first row was the
  state *after* one step, so stochastic and deterministic trajectories were
  offset by one time step.
* :meth:`NetworkEpiModel.simulate <epidemik.NetworkEpiModel.simulate>`
  indexes ``population_`` and ``values_`` by time starting at ``t_min``
  (default 1) rather than 0, labels the columns of ``population_`` by node,
  and always includes every compartment of the model in ``values_`` (in
  model order).
* The initial conditions of age-structured models are now split across age
  groups deterministically in :meth:`~epidemik.EpiModel.integrate`
  (previously a random multinomial draw made integration
  non-reproducible).
* The first argument of :meth:`MetaEpiModel.simulate
  <epidemik.MetaEpiModel.simulate>` is now called ``timesteps``; the old
  ``timestamp`` spelling still works but emits a ``DeprecationWarning``.
  The simulation always runs ``timesteps`` steps regardless of ``t_min``.
* :meth:`NetworkEpiModel.R0 <epidemik.NetworkEpiModel.R0>` no longer
  requires a compartment named ``"R"`` to exist.
* Unknown compartment names passed as initial conditions to
  :meth:`~epidemik.EpiModel.integrate` or
  :meth:`~epidemik.EpiModel.simulate` are ignored with a warning (they were
  silently ignored by ``integrate`` and crashed ``simulate``).
* ``model["X"]`` and :meth:`~epidemik.EpiModel.plot` raise
  :class:`~epidemik.utils.NotInitialized` when the model has not been run
  yet; the constructor of :class:`~epidemik.MetaEpiModel` validates the
  travel matrix and populations and raises :class:`ValueError` on bad
  input.

Bug fixes
~~~~~~~~~

* :meth:`MetaEpiModel.plot <epidemik.MetaEpiModel.plot>` crashed unless a
  sub-population was called ``"NY"``, and used the time label of the peak
  as a positional index. :meth:`MetaEpiModel.draw_model
  <epidemik.MetaEpiModel.draw_model>` always crashed.
* :meth:`MetaEpiModel.simulate <epidemik.MetaEpiModel.simulate>` ignored
  its ``seasonality`` argument.
* :class:`~epidemik.NetworkEpiModel` truncated compartment names to a single
  character during simulation, so multi-letter compartments never
  transitioned; only the last spontaneous transition out of a compartment
  was kept; and non-integer node labels were ignored.
* :meth:`EpiModel.simulate <epidemik.EpiModel.simulate>` crashed with
  ``pvals < 0`` when the transition probabilities out of a compartment
  exceeded 1, crashed on age-structured models, and overwrote (instead of
  adding) the contributions of several transitions with the same source
  and target.
* Multi-group models produced ``NaN`` when a group was empty.
* :meth:`EpiModel.R0 <epidemik.EpiModel.R0>` returned ``None`` for models
  with spontaneous transitions out of non-infectious compartments (e.g.
  waning immunity ``R -> S``) and hid every error behind a bare ``except``.
* :meth:`EpiModel.plot <epidemik.EpiModel.plot>` ignored the ``ax``
  argument, and :meth:`EpiModel.draw_model <epidemik.EpiModel.draw_model>`
  drew the labels on the current axes instead of ``ax``.
* ``model["X"]`` raised ``KeyError`` for unknown names on models without
  age structure instead of returning ``None``.
* :meth:`EpiModel.add_age_structure <epidemik.EpiModel.add_age_structure>`
  validates its inputs with :class:`ValueError` instead of ``assert``.

Documentation
~~~~~~~~~~~~~

* Every public class and method has a complete NumPy-style docstring, the
  Sphinx configuration uses :mod:`sphinx.ext.napoleon` and
  :mod:`sphinx.ext.intersphinx`, and the user guide covers
  :class:`~epidemik.EpiModel`, :class:`~epidemik.NetworkEpiModel` and
  :class:`~epidemik.MetaEpiModel` in detail.

0.1 (2024-06-10)
----------------

* Type hints across the public API.
* Sphinx documentation with autodoc support.

0.0.17 -- 0.0.22 (2024-04 -- 2024-05)
-------------------------------------

* Added :class:`~epidemik.MetaEpiModel`.
* Added the test suite and continuous integration.
* Improved the plotting functions of :class:`~epidemik.EpiModel`.
