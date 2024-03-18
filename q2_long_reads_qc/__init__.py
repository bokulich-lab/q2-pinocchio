# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from ._version import get_versions
from .build_index import build_index
from .classify_consensus import classify_consensus, find_consensus_annotation
from .extract_seqs import extract_seqs
from .filter_reads import filter_reads
from .minimap2 import minimap2

__version__ = get_versions()["version"]
del get_versions

__all__ = [
    "filter_reads",
    "build_index",
    "minimap2",
    "classify_consensus",
    "find_consensus_annotation",
    "extract_seqs",
]
