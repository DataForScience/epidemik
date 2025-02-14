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
