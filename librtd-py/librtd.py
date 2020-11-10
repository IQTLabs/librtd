# flake8: noqa

import pkg_resources

__version__ = pkg_resources.get_distribution("librtd").version

import nimporter
from librtdpy import (
    same_kmer_return_times,
    pairwise_kmer_return_times,
    reverse_complement_return_times,
    return_time_distribution,
)
