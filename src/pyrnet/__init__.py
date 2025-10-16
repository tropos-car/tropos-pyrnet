
from importlib.metadata import version, PackageNotFoundError
try:
    __version__ = version("pyrnet")
except PackageNotFoundError:
    # package is not installed
    pass