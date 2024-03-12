# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os

import pandas as pd

from ..plugin_setup import plugin
from ._format import PAFDirectoryFormat, PAFFormat


@plugin.register_transformer
def _1(data: pd.DataFrame) -> PAFDirectoryFormat:
    # Create a new PAFDirectoryFormat object
    ff = PAFDirectoryFormat()
    paf_file_path = os.path.join(str(ff), "output.paf")
    data.to_csv(paf_file_path, sep="\t", index=False, header=False)
    return ff


@plugin.register_transformer
def _2(data: PAFFormat) -> pd.DataFrame:
    file = pd.read_csv(str(data), sep="\t", header=None)
    return file
