# ----------------------------------------------------------------------------
# Copyright (c) 2024, Bokulich Lab.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from ._version import get_versions
from .build_index import build_index
from .classify_consensus import _find_consensus_annotation, classify_consensus_minimap2
from .extract_reads import extract_reads
from .filter_reads import filter_reads
from .minimap2_search import minimap2_search
from .nanoplot_stats import stats
from .trim_long_reads import trim

__version__ = get_versions()["version"]
del get_versions

__all__ = [
    "filter_reads",
    "build_index",
    "minimap2_search",
    "classify_consensus_minimap2",
    "_find_consensus_annotation",
    "extract_reads",
    "stats",
    "trim",
]
