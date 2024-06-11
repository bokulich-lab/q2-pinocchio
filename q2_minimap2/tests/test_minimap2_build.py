# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import subprocess
from unittest.mock import patch

from q2_types.feature_data import DNAFASTAFormat

from q2_minimap2.build_index import build_index

from .test_minimap2 import Minimap2TestsBase


class TestMinimap2Build(Minimap2TestsBase):
    def setUp(self):
        super().setUp()
        # Retrieve the file path for the test DNA sequences
        self.genome = DNAFASTAFormat(
            self.get_data_path("minimap2_build/dna-sequences.fasta"), mode="r"
        )

    def test_build(self):
        build_index(self.genome)

    def test_build_index_subprocess_failure(self):
        # Patch the subprocess.run or the internal function that calls subprocess.run
        with patch("subprocess.run") as mocked_run:
            # Configuring the mock to raise a CalledProcessError
            mocked_run.side_effect = subprocess.CalledProcessError(
                returncode=1, cmd="minimap2 -d ..."
            )

            # Calling the method under test
            with self.assertRaises(Exception) as context:
                build_index(self.genome)

            # Checking if the raised exception is as expected
            self.assertTrue(
                "An error was encountered while running main.nf"
                in str(context.exception)
            )
            self.assertTrue("(return code 1)" in str(context.exception))
