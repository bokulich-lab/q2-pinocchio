# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import subprocess
from unittest.mock import patch

from q2_types.feature_data import FeatureData, Sequence
from qiime2 import Artifact

from .test_minimap2 import Minimap2TestsBase


class TestMinimap2Build(Minimap2TestsBase):
    def setUp(self):
        super().setUp()
        # Retrieve the file path for the test DNA sequences in FASTA format
        fasta_path = self.get_data_path("minimap2_build/dna-sequences.fasta")
        # Import the DNA sequences from the FASTA file into a QIIME 2 Artifact
        self.genome = Artifact.import_data(FeatureData[Sequence], fasta_path)

    def test_build(self):
        self.plugin.methods["build_index"](self.genome)

    def test_build_index_subprocess_failure(self):
        # Patch the subprocess.run or the internal function that calls subprocess.run
        with patch("subprocess.run") as mocked_run:
            # Configuring the mock to raise a CalledProcessError
            mocked_run.side_effect = subprocess.CalledProcessError(
                returncode=1, cmd="minimap2 -d ..."
            )

            # Calling the method under test
            with self.assertRaises(Exception) as context:
                self.plugin.methods["build_index"](self.genome)

            # Checking if the raised exception is as expected
            self.assertTrue(
                "An error was encountered while running main.nf"
                in str(context.exception)
            )
            self.assertTrue("(return code 1)" in str(context.exception))
