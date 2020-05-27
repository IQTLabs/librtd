# flake8: noqa

import pkg_resources

__version__ = pkg_resources.get_distribution("librtd").version

import nimporter
from .librtd import returnTimeDistribution
