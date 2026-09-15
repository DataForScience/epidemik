Installation
============

``epidemik`` is published on PyPI and supports Python 3.8 and above. The
latest release requires NumPy 1.20+, SciPy 1.10+, pandas 2.0+ and
NetworkX 3+.

From PyPI
---------

The recommended way to install the package is via :program:`pip`:

.. code-block:: bash

   pip install epidemik

From source
-----------

The repository is managed with `uv <https://docs.astral.sh/uv/>`_. To set
up a development environment with the package installed in editable mode
and the test tools pinned in ``uv.lock``:

.. code-block:: bash

   git clone https://github.com/DataForScience/epidemik.git
   cd epidemik
   uv sync                 # add --group docs to also install Sphinx
   uv run pytest           # run the test suite
   uv build                # build the wheel and sdist

A plain ``pip install -e .`` also works if you prefer not to use uv.

Dependencies
------------

``epidemik`` relies on a handful of well-known scientific Python libraries:

* `NumPy <https://numpy.org/>`_ for array maths and random number generation.
* `SciPy <https://scipy.org/>`_ for ODE integration.
* `NetworkX <https://networkx.org/>`_ as the underlying graph representation
  of the compartmental structure and contact networks.
* `pandas <https://pandas.pydata.org/>`_ for tabular outputs.
* `Matplotlib <https://matplotlib.org/>`_ for plotting helpers.
* `tqdm <https://tqdm.github.io/>`_ for progress bars in
  :class:`~epidemik.MetaEpiModel` simulations.

These dependencies are installed automatically by :program:`pip`. For nicer
:meth:`~epidemik.EpiModel.draw_model` output, you may additionally install
`pygraphviz <https://pygraphviz.github.io/>`_ -- when present, it is used
to lay out the compartment graph; otherwise a spectral fallback is used.
