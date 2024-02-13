# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import model


class Minimap2IndexFileFormat(model.BinaryFileFormat):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def _validate_(self, level):
        pass


Minimap2IndexFileDirFmt = model.SingleFileDirectoryFormat(
    "Minimap2IndexFileDirFmt", "index.mmi", Minimap2IndexFileFormat
)
