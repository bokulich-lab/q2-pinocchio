# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from qiime2 import Artifact

from .test_long_reads_qc import LongReadsQCTestsBase


class TestMinimap2Build(LongReadsQCTestsBase):

    # This test just makes sure this runs without error, which will include
    # validating the contents.
    def test_build(self):
        genome = Artifact.load(self.get_data_path("sars2-genome.qza"))
        self.plugin.methods["minimap2_build"](genome)
