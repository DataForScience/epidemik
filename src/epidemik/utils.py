"""Shared utilities used across the :mod:`epidemik` package.

This module exposes a small custom exception used to signal that an
:class:`~epidemik.EpiModel` or :class:`~epidemik.MetaEpiModel` has not been
fully initialised yet, as well as the default colour palette used by the
plotting helpers (:meth:`~epidemik.EpiModel.plot` and
:meth:`~epidemik.EpiModel.draw_model`).
"""

from collections import defaultdict


class NotInitialized(Exception):
    """Raised when a model is used before it has been fully initialised.

    Typical situations include calling :meth:`~epidemik.EpiModel.plot`
    before either :meth:`~epidemik.EpiModel.integrate` or
    :meth:`~epidemik.EpiModel.simulate`, or calling
    :meth:`~epidemik.MetaEpiModel.simulate` without specifying a
    ``seed_state``.
    """


#: Default colours used by the plotting helpers, keyed by the first letter
#: of the compartment name (``"S"`` -> blue, ``"I"`` -> magenta, ...). Any
#: unknown compartment falls back to a sensible default orange.
epi_colors = defaultdict(lambda: '#f39019')
epi_colors['S'] = '#51a7f9'
epi_colors['E'] = '#f9e351'
epi_colors['I'] = '#cf51f9'
epi_colors['R'] = '#70bf41'
epi_colors['D'] = '#8b8b8b'