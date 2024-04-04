# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from ._version import get_versions
from .build_index import build_index
from .classify_consensus import _find_consensus_annotation, classify_consensus_minimap2
from .extract_seqs import extract_seqs
from .filter_paired_end_reads import filter_paired_end_reads
from .filter_single_end_reads import filter_single_end_reads
from .minimap2_search import minimap2_search

__version__ = get_versions()["version"]
del get_versions

__all__ = [
    "filter_single_end_reads",
    "filter_paired_end_reads",
    "build_index",
    "minimap2_search",
    "classify_consensus_minimap2",
    "_find_consensus_annotation",
    "extract_seqs",
]
