from collections import defaultdict
import os
import sys
from urllib.parse import urlparse
from typing import Union

__all__ = [
    "OFFICIAL_REPO",
    "EPI_COLORS",
    "NotInitialized",
    "get_cache_directory",
    "NotImplementedError",
    "Parameters",
]


class NotInitialized(Exception):
    pass


class NotImplementedError(Exception):
    pass


APPNAME = os.path.basename(os.path.realpath("."))
OFFICIAL_REPO = "https://github.com/DataForScience/epidemik/"

EPI_COLORS = defaultdict(lambda: "#f39019")
EPI_COLORS["S"] = "#51a7f9"
EPI_COLORS["E"] = "#f9e351"
EPI_COLORS["I"] = "#cf51f9"
EPI_COLORS["R"] = "#70bf41"
EPI_COLORS["D"] = "#8b8b8b"


class Parameters(dict):
    def __init__(self):
        super().__init__()
        self.globals = {"__builtins__": None}

    def __setitem__(self, key, value):
        self.define_parameters(**{key: value})

    def __getitem__(self, key):
        return self.compute_parameter(key)

    def define_parameters(self, **kwargs) -> None:
        """
        Define one or more parameter for the model

        Parameters:
        - kwargs: keyword arguments
            Named parameters for the model

        Returns:
        None
        """
        for key, value in kwargs.items():
            if isinstance(value, str):
                try:
                    # Convert floats written as strings to float
                    value = float(value.strip())
                except:
                    pass
            
            super().__setitem__(key, value)

    def compute_parameter(self, param: Union[str]) -> float:
        """
        Compute the rate from a string

        Parameters:
        - rate: string
            Rate of the transition

        Returns:
        float
            The computed rate
        """
        import logging

        if param in self.keys():
            if isinstance(super().__getitem__(param), (int, float)):
                return super().__getitem__(param)
            else:
                try:
                    return eval(super().__getitem__(param), self.globals, self)
                except Exception as e:
                    logging.error(f"Error computing parameter {param}: {e}")
                    return None

def get_cache_directory():
    """
    Return the location of the cache directory for the current platform.
    """
    system = sys.platform

    if system == "darwin":
        path = os.path.join(os.path.expanduser("~/Library/Caches"), APPNAME)
    else:
        path = os.getenv("XDG_CACHE_HOME", os.path.expanduser("~/.cache"))
        path = os.path.join(path, APPNAME)

    return path


def get_remote_path(repo: Union[str, None] = None) -> str:
    """
    Return the location of the remote cache directory.
    """
    if repo is None:
        repo = OFFICIAL_REPO

    parsed_repo = urlparse(repo)

    if parsed_repo.netloc == "github.com":
        repo = repo.replace("github.com", "raw.githubusercontent.com")
        remote_path = repo + os.path.join("refs/heads/models/models/")
    else:
        None

    return remote_path
