API reference
=============

This section documents the public API of every module in the
:mod:`epidemik` package. The content is generated automatically from the
docstrings in the source code.

.. automodule:: epidemik
   :no-members:

.. currentmodule:: epidemik

Package overview
----------------

.. autosummary::
   :nosignatures:

   EpiModel
   NetworkEpiModel
   MetaEpiModel
   utils

EpiModel
--------

.. autoclass:: epidemik.EpiModel
   :members:
   :exclude-members: __init__
   :special-members: __init__, __repr__, __getitem__

NetworkEpiModel
---------------

.. autoclass:: epidemik.NetworkEpiModel
   :members:
   :exclude-members: __init__
   :special-members: __init__
   :show-inheritance:

MetaEpiModel
------------

.. autoclass:: epidemik.MetaEpiModel
   :members:
   :exclude-members: __init__
   :special-members: __init__, __repr__

epidemik.utils
--------------

.. automodule:: epidemik.utils
   :members:
   :exclude-members: epi_colors

.. autodata:: epidemik.utils.epi_colors
   :annotation: = <colour palette keyed by compartment initial>
