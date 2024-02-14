# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from ._version import get_versions
from .longqc import evaluate
from .minimap2 import build_minimap2_index  # , filter_reads

__version__ = get_versions()["version"]
del get_versions

__all__ = ["evaluate", "filter_reads", "build_minimap2_index"]
