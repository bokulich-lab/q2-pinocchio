# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import shutil
import tempfile
import unittest

from q2_long_reads_qc.minimap2.minimap2_search import (
    filter_by_maxaccepts,
    filter_by_perc_identity,
)

from .test_long_reads_qc import LongReadsQCTestsBase


class TestFilterByMaxAccepts(LongReadsQCTestsBase):
    def setUp(self):
        super().setUp()
        # Retrieve the paths to test data files
        self.initial_PairwiseAlignmentMN2_file = self.get_data_path(
            "initial_PairwiseAlignmentMN2_file.PairwiseAlignmentMN2"
        )
        self.exp_PairwiseAlignmentMN2_f_max1 = self.get_data_path(
            "expected_PairwiseAlignmentMN2_file_max1.PairwiseAlignmentMN2"
        )
        self.exp_PairwiseAlignmentMN2_f_max2 = self.get_data_path(
            "expected_PairwiseAlignmentMN2_file_max2.PairwiseAlignmentMN2"
        )
        self.exp_PairwiseAlignmentMN2_f_max3 = self.get_data_path(
            "expected_PairwiseAlignmentMN2_file_max3.PairwiseAlignmentMN2"
        )

    def test_filter_by_maxaccepts_max1(self):
        self._test_filter_by_maxaccepts(1, self.exp_PairwiseAlignmentMN2_f_max1)

    def test_filter_by_maxaccepts_max2(self):
        self._test_filter_by_maxaccepts(2, self.exp_PairwiseAlignmentMN2_f_max2)

    def test_filter_by_maxaccepts_max3(self):
        self._test_filter_by_maxaccepts(3, self.exp_PairwiseAlignmentMN2_f_max3)

    def _test_filter_by_maxaccepts(self, max_accepts, expected_file_path):
        # Create a temporary file to use as the output destination
        with tempfile.NamedTemporaryFile(mode="w+", delete=False) as tmp_file:
            temp_file_path = tmp_file.name

        # Copy the content of the initial PairwiseAlignmentMN2 file to the temp file
        shutil.copy(self.initial_PairwiseAlignmentMN2_file, temp_file_path)

        # Call the function with the temp file as input and output, and
        # specified maxaccepts
        filter_by_maxaccepts(temp_file_path, max_accepts)

        # Compare the content of the temp file with the expected
        # PairwiseAlignmentMN2 file for specified maxaccepts
        with open(temp_file_path, "r") as temp_file, open(
            expected_file_path, "r"
        ) as expected_file:
            temp_content = temp_file.read()
            expected_content = expected_file.read()
            self.assertEqual(
                temp_content,
                expected_content,
                "Filtered PairwiseAlignmentMN2 file content does not match expected "
                f"output for maxaccepts {max_accepts}",
            )

        # Clean up the temporary file
        os.remove(temp_file_path)


class TestFilterByPercIdentity(LongReadsQCTestsBase):
    def setUp(self):
        super().setUp()
        self.PairwiseAlignmentMN2_file = self.get_data_path(
            "initial_PairwiseAlignmentMN2_file.PairwiseAlignmentMN2"
        )
        self.expected_PairwiseAlignmentMN2_file_perc_85 = self.get_data_path(
            "expected_PairwiseAlignmentMN2_file_perc_85.PairwiseAlignmentMN2"
        )
        self.expected_PairwiseAlignmentMN2_file_perc_80 = self.get_data_path(
            "expected_PairwiseAlignmentMN2_file_perc_80.PairwiseAlignmentMN2"
        )

    def test_filter_by_perc_identity_85(self):
        self._test_filter_by_perc_identity(
            0.85, self.expected_PairwiseAlignmentMN2_file_perc_85
        )

    def test_filter_by_perc_identity_80(self):
        self._test_filter_by_perc_identity(
            0.8, self.expected_PairwiseAlignmentMN2_file_perc_80
        )

    def _test_filter_by_perc_identity(self, perc_identity, expected_file_path):
        with tempfile.NamedTemporaryFile(mode="w+", delete=False) as tmp_file:
            temp_file_path = tmp_file.name

        shutil.copy(self.PairwiseAlignmentMN2_file, temp_file_path)

        filter_by_perc_identity(temp_file_path, perc_identity)

        with open(temp_file_path, "r") as temp_file, open(
            expected_file_path, "r"
        ) as expected_file:
            temp_content = temp_file.read()
            expected_content = expected_file.read()
            self.assertEqual(
                temp_content,
                expected_content,
                "Filtered PairwiseAlignmentMN2 file content does not match expected "
                f"output for perc_identity {perc_identity}",
            )

        os.remove(temp_file_path)


if __name__ == "__main__":
    unittest.main()
