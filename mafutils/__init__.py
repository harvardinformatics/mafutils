"""Reusable MAF toolkit package for indexing, fetching, and summarizing MAFs."""

from importlib.metadata import PackageNotFoundError, version as _installed_version

try:
    __version__ = _installed_version("mafutils")
except PackageNotFoundError:
    # Running from a source tree that was never installed (no distribution
    # metadata to read). setuptools-scm derives the real version from git
    # tags at build/install time, so there's nothing to fall back to here.
    __version__ = "unknown"

__all__ = ["cli", "fetch", "index", "stats", "__version__"]
