# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from ._version import get_versions
from .filtering import filter_reads
from .minimap2 import minimap2_build

__version__ = get_versions()["version"]
del get_versions

__all__ = ["filter_reads", "minimap2_build"]
