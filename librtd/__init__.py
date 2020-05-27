import pkg_resources

__version__ = pkg_resources.get_distribution("mechwolf").version

import nimporter
from .librtd import returnTimeDistribution
