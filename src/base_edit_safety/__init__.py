"""Review evidence for base-editing guides at the HBG1/HBG2 promoters."""
from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("base-edit-safety")
except PackageNotFoundError:
    __version__ = "0+source-tree"
