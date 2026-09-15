"""Shared utilities used across the :mod:`epidemik` package.

This module exposes:

* :class:`NotInitialized`, raised when a model is used before it has been
  run or fully set up;
* :class:`Parameters`, the dictionary of named model parameters (which
  may be numbers or expressions in terms of other parameters);
* :data:`epi_colors`, the default colour palette used by the plotting
  helpers;
* helpers used by :meth:`EpiModel.download_model
  <epidemik.EpiModel.download_model>` and :meth:`EpiModel.list_models
  <epidemik.EpiModel.list_models>` to locate the model repository and the
  local cache.
"""

from collections import defaultdict
import math
import os
import sys
from typing import Dict, Union
from urllib.parse import urlparse

__all__ = [
    "NotInitialized",
    "Parameters",
    "epi_colors",
    "EPI_COLORS",
    "OFFICIAL_REPO",
    "get_cache_directory",
    "get_remote_path",
]


class NotInitialized(Exception):
    """Raised when a model is used before it has been fully initialised.

    Typical situations include calling :meth:`~epidemik.EpiModel.plot`
    before either :meth:`~epidemik.EpiModel.integrate` or
    :meth:`~epidemik.EpiModel.simulate`, or calling
    :meth:`~epidemik.MetaEpiModel.simulate` without specifying a
    ``seed_state``.
    """


#: Names that expressions may use besides the other parameters.
_SAFE_FUNCTIONS = {
    name: getattr(math, name)
    for name in ("exp", "log", "log10", "sqrt", "sin", "cos", "tan", "pi", "e")
}
_SAFE_FUNCTIONS.update({"abs": abs, "min": min, "max": max})


class Parameters(dict):
    """Named parameters of a model.

    Behaves like a dictionary whose values are either numbers or *string
    expressions* written in terms of other parameters, e.g. ``"beta/2"``
    or ``"2*mu"``. Reading a key (``params["mu"]``) always returns the
    evaluated numeric value, while :meth:`raw` returns what was stored.

    Examples
    --------
    >>> params = Parameters(beta=0.2)
    >>> params["mu"] = "beta/2"
    >>> params["mu"]
    0.1
    >>> params.raw("mu")
    'beta/2'
    """

    def __init__(self, **kwargs):
        super().__init__()
        self.define_parameters(**kwargs)

    def __setitem__(self, key, value):
        self.define_parameters(**{key: value})

    def __getitem__(self, key):
        return self.compute_parameter(key)

    def get(self, key, default=None):
        if key not in self:
            return default
        return self.compute_parameter(key)

    def define_parameters(self, **kwargs) -> None:
        """Define one or more parameters.

        Parameters
        ----------
        **kwargs
            ``name=value`` pairs. Numeric strings (``"0.2"``) are converted
            to floats; any other string is kept as an expression to be
            evaluated lazily against the other parameters.

        Raises
        ------
        TypeError
            If a value is neither a number nor a string.
        """
        for key, value in kwargs.items():
            if isinstance(value, str):
                try:
                    value = float(value.strip())
                except ValueError:
                    value = value.strip()
            elif isinstance(value, bool) or not isinstance(value, (int, float)):
                try:
                    value = float(value)
                except (TypeError, ValueError):
                    raise TypeError(
                        "Parameter '%s' must be a number or an expression string, "
                        "got %r." % (key, value)
                    )

            super().__setitem__(key, value)

    def raw(self, key):
        """Return the stored value of ``key`` (number or expression string)."""
        return super().__getitem__(key)

    def compute_parameter(self, param: str, _stack=()) -> float:
        """Evaluate a parameter.

        Parameters
        ----------
        param : str
            Name of the parameter.

        Returns
        -------
        float
            The numeric value of the parameter, evaluating any expression
            (recursively) in terms of the other parameters.

        Raises
        ------
        KeyError
            If ``param`` (or a parameter it refers to) is not defined.
        ValueError
            If the expression cannot be evaluated, or if parameters refer
            to each other in a cycle.
        """
        if param not in self.keys():
            raise KeyError("Unknown parameter '%s'." % param)

        value = super().__getitem__(param)

        if not isinstance(value, str):
            return float(value)

        if param in _stack:
            raise ValueError(
                "Circular definition while evaluating parameter '%s'." % param
            )

        namespace = _ParameterNamespace(self, _stack + (param,))

        try:
            result = eval(value, {"__builtins__": {}}, namespace)
        except (KeyError, ValueError):
            raise
        except NameError as e:
            raise KeyError("Parameter '%s' = %r refers to an unknown name: %s" % (param, value, e))
        except Exception as e:
            raise ValueError(
                "Unable to evaluate parameter '%s' = %r: %s" % (param, value, e)
            )

        try:
            return float(result)
        except (TypeError, ValueError):
            raise ValueError(
                "Parameter '%s' = %r does not evaluate to a number." % (param, value)
            )

    def resolve(self) -> Dict[str, float]:
        """Return a plain ``{name: float}`` dictionary with every parameter evaluated."""
        return {key: self.compute_parameter(key) for key in self.keys()}

    def copy(self) -> "Parameters":
        new = Parameters()
        for key in self.keys():
            new.define_parameters(**{key: super().__getitem__(key)})
        return new

    def __repr__(self) -> str:
        return "Parameters(%s)" % ", ".join(
            "%s=%r" % (key, super(Parameters, self).__getitem__(key)) for key in self.keys()
        )


class _ParameterNamespace(dict):
    """Lookup table used when evaluating parameter expressions."""

    def __init__(self, params: Parameters, stack):
        super().__init__(_SAFE_FUNCTIONS)
        self._params = params
        self._stack = stack

    def __missing__(self, key):
        return self._params.compute_parameter(key, self._stack)


#: Default colours used by the plotting helpers, keyed by the first letter
#: of the compartment name (``"S"`` -> blue, ``"E"`` -> yellow, ``"I"`` ->
#: magenta, ``"R"`` -> green, ``"D"`` -> grey). Any other compartment falls
#: back to orange.
epi_colors = defaultdict(lambda: '#f39019')
epi_colors['S'] = '#51a7f9'
epi_colors['E'] = '#f9e351'
epi_colors['I'] = '#cf51f9'
epi_colors['R'] = '#70bf41'
epi_colors['D'] = '#8b8b8b'

#: Alias of :data:`epi_colors` kept for backwards compatibility.
EPI_COLORS = epi_colors

#: GitHub repository hosting the official library of model files.
OFFICIAL_REPO = "https://github.com/DataForScience/epidemik/"

#: Location of the model files inside :data:`OFFICIAL_REPO` (the ``models/``
#: directory of the ``main`` branch).
_MODELS_LOCATION = "refs/heads/main/models/"


def get_cache_directory() -> str:
    """Return the directory where downloaded model files are cached.

    The location can be overridden with the ``EPIDEMIK_CACHE_DIR``
    environment variable. Otherwise ``~/Library/Caches/epidemik`` is used
    on macOS and ``$XDG_CACHE_HOME/epidemik`` (defaulting to
    ``~/.cache/epidemik``) elsewhere.
    """
    override = os.environ.get("EPIDEMIK_CACHE_DIR")

    if override:
        return override

    if sys.platform == "darwin":
        base = os.path.expanduser("~/Library/Caches")
    else:
        base = os.getenv("XDG_CACHE_HOME", os.path.expanduser("~/.cache"))

    return os.path.join(base, "epidemik")


def get_remote_path(repo: Union[str, None] = None) -> str:
    """Return the location of the model files in a repository.

    Parameters
    ----------
    repo : str, optional
        Either the URL of a GitHub repository (defaults to
        :data:`OFFICIAL_REPO`) or the path of a local directory containing
        model files.

    Returns
    -------
    str
        Base URL (ending with ``/``) from which ``<name>.yaml`` files can
        be fetched, or the local directory path.

    Raises
    ------
    ValueError
        If ``repo`` is neither a GitHub URL nor an existing directory.
    """
    if repo is None:
        repo = OFFICIAL_REPO

    if os.path.isdir(repo):
        return repo

    parsed_repo = urlparse(repo)

    if parsed_repo.netloc == "github.com":
        repo = repo.replace("github.com", "raw.githubusercontent.com")

        if not repo.endswith("/"):
            repo += "/"

        return repo + _MODELS_LOCATION

    raise ValueError(
        "repo must be a GitHub repository URL or an existing local directory, got %r." % repo
    )
