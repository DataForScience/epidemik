"""Sphinx configuration for the epidemik package.

See https://www.sphinx-doc.org/en/master/usage/configuration.html for the full
list of options.
"""

from __future__ import annotations

import os
import sys
from datetime import date

# Make the in-tree source importable so autodoc works without an install
# (e.g. when building the docs locally before ``pip install -e .``).
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "src")))

import epidemik  # noqa: E402  (import after sys.path tweak)

# -- Project information -----------------------------------------------------

project = "epidemik"
author = "Bruno Gonçalves"
copyright = f"2024-{date.today().year}, {author}"
release = epidemik.__version__
version = release

# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.duration",
    "sphinx.ext.doctest",
    "sphinx.ext.mathjax",
    "sphinx_rtd_theme",
]

templates_path = ["source/_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# Generate stub pages for autosummary directives automatically.
autosummary_generate = True

autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "show-inheritance": True,
    "member-order": "bysource",
}

# Napoleon settings -- the in-source docstrings follow NumPy conventions.
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_param = True
napoleon_use_rtype = True

# Cross-references to external project documentation. Set the
# ``EPIDEMIK_DOCS_OFFLINE`` environment variable (to any truthy value) to
# skip the network fetches when building the documentation locally.
if os.environ.get("EPIDEMIK_DOCS_OFFLINE"):
    intersphinx_mapping = {}
else:
    intersphinx_mapping = {
        "python": ("https://docs.python.org/3", None),
        "numpy": ("https://numpy.org/doc/stable/", None),
        "scipy": ("https://docs.scipy.org/doc/scipy/", None),
        "pandas": ("https://pandas.pydata.org/docs/", None),
        "networkx": ("https://networkx.org/documentation/stable/", None),
        "matplotlib": ("https://matplotlib.org/stable/", None),
    }

# Don't hang the build forever if a remote inventory is slow.
intersphinx_timeout = 10

# -- Options for HTML output -------------------------------------------------

html_theme = "sphinx_rtd_theme"
html_static_path = ["source/_static"]
html_title = f"epidemik {release} documentation"
html_theme_options = {
    "logo_only": False,
    "prev_next_buttons_location": "bottom",
    "style_external_links": True,
    "style_nav_header_background": "white",
    # Toc options
    "collapse_navigation": True,
    "sticky_navigation": True,
    "navigation_depth": 4,
    "includehidden": True,
    "titles_only": False,
}

# Object signatures look nicer when split across lines for long parameter lists.
maximum_signature_line_length = 88
