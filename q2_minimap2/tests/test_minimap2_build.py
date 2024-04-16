# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from q2_types.feature_data import FeatureData, Sequence
from qiime2 import Artifact

from .test_minimap2 import Minimap2TestsBase


class TestMinimap2Build(Minimap2TestsBase):
    def test_build(self):
        fasta_path = self.get_data_path("dna-sequences.fasta")
        genome = Artifact.import_data(FeatureData[Sequence], fasta_path)
        self.plugin.methods["build_index"](genome)
