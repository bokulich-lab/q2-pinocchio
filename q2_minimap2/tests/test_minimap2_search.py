# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tempfile
import unittest

import pandas as pd
import pandas.testing as pdt
from pandas.testing import assert_frame_equal
from qiime2 import Artifact

from q2_minimap2.minimap2_search import (
    filter_by_maxaccepts,
    filter_by_perc_identity,
    minimap2_search,
)

from .test_minimap2 import Minimap2TestsBase


class TestFilterByMaxAccepts(Minimap2TestsBase):
    def setUp(self):
        super().setUp()
        # Retrieve the paths to test data files
        self.initial_PairwiseAlignmentMN2_file = self.get_data_path(
            "minimap2_search/initial_paf_file.paf"
        )
        self.exp_PairwiseAlignmentMN2_f_max1 = self.get_data_path(
            "minimap2_search/expected_paf_file_max1.paf"
        )
        self.exp_PairwiseAlignmentMN2_f_max2 = self.get_data_path(
            "minimap2_search/expected_paf_file_max2.paf"
        )
        self.exp_PairwiseAlignmentMN2_f_max3 = self.get_data_path(
            "minimap2_search/expected_paf_file_max3.paf"
        )

    def test_filter_by_maxaccepts_max1(self):
        self._test_filter_by_maxaccepts(1, self.exp_PairwiseAlignmentMN2_f_max1)

    def test_filter_by_maxaccepts_max2(self):
        self._test_filter_by_maxaccepts(2, self.exp_PairwiseAlignmentMN2_f_max2)

    def test_filter_by_maxaccepts_max3(self):
        self._test_filter_by_maxaccepts(3, self.exp_PairwiseAlignmentMN2_f_max3)

    def _test_filter_by_maxaccepts(self, max_accepts, expected_file_path):
        # Load input and expected data
        input_df = pd.read_csv(
            self.initial_PairwiseAlignmentMN2_file, sep="\t", header=None
        )
        expected_df = pd.read_csv(expected_file_path, sep="\t", header=None)

        # Generate results from function
        result_df = filter_by_maxaccepts(input_df, max_accepts)
        result_df.reset_index(drop=True, inplace=True)

        # Assert that the two data frames are equal
        assert_frame_equal(result_df, expected_df)


class TestFilterByPercIdentity(Minimap2TestsBase):
    def setUp(self):
        super().setUp()
        self.PairwiseAlignmentMN2_file = self.get_data_path(
            "minimap2_search/initial_paf_file.paf"
        )
        self.expected_PairwiseAlignmentMN2_file_perc_85 = self.get_data_path(
            "minimap2_search/expected_paf_file_perc_85.paf"
        )
        self.expected_PairwiseAlignmentMN2_file_perc_80 = self.get_data_path(
            "minimap2_search/expected_paf_file_perc_80.paf"
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
        # Load input and expected data
        input_df = pd.read_csv(self.PairwiseAlignmentMN2_file, sep="\t", header=None)
        expected_df = pd.read_csv(expected_file_path, sep="\t", header=None)

        # Generate results from function
        result_df = filter_by_perc_identity(input_df, perc_identity, True)
        result_df.reset_index(drop=True, inplace=True)

        # Assert that the two data frames are equal
        assert_frame_equal(result_df, expected_df)


class TestMinimap2(Minimap2TestsBase):
    def setUp(self):
        super().setUp()

        self.query_reads = Artifact.load(
            self.get_data_path("minimap2_search/query-seqs.qza")
        )
        self.index_database = Artifact.load(
            self.get_data_path("minimap2_search/minimap2_test_index.qza")
        )
        self.ref = Artifact.load(
            self.get_data_path("minimap2_search/se-dna-sequences.qza")
        )
        self.paf_true = Artifact.load(
            self.get_data_path("minimap2_search/minimap2_test_paf.qza")
        )
        self.paf_true_only_hits = Artifact.load(
            self.get_data_path("minimap2_search/minimap2_only_hits_test_paf.qza")
        )

    def test_minimap2(self):
        (output_artifact,) = self.plugin.methods["minimap2_search"](
            self.query_reads, self.index_database
        )

        # Create a temporary directory to hold the extracted sequences
        with tempfile.TemporaryDirectory() as temp_dir:
            output_artifact.export_data(temp_dir)
            output_fp = os.path.join(temp_dir, "output.paf")

            with tempfile.TemporaryDirectory() as temp_dir2:

                self.paf_true.export_data(temp_dir2)
                true_paf_fp = os.path.join(temp_dir2, "output.paf")

                # Compare the contents of the files
                with open(true_paf_fp, "r") as file1, open(output_fp, "r") as file2:
                    true_paf_content = file1.read()
                    output_paf_content = file2.read()

                    self.assertEqual(true_paf_content, output_paf_content)

    def test_minimap2_only_hits(self):
        (output_artifact,) = self.plugin.methods["minimap2_search"](
            self.query_reads,
            self.index_database,
            output_no_hits=False,
        )

        # Create a temporary directory to hold the extracted sequences
        with tempfile.TemporaryDirectory() as temp_dir:
            output_artifact.export_data(temp_dir)
            output_fp = os.path.join(temp_dir, "output.paf")

            with tempfile.TemporaryDirectory() as temp_dir2:

                self.paf_true_only_hits.export_data(temp_dir2)
                true_paf_fp = os.path.join(temp_dir2, "output.paf")

                # Compare the contents of the files
                with open(true_paf_fp, "r") as file1, open(output_fp, "r") as file2:
                    true_paf_content = file1.read()
                    output_paf_content = file2.read()

                    self.assertEqual(true_paf_content, output_paf_content)

    def test_minimap2_output_consistency(self):
        (result1,) = self.plugin.methods["minimap2_search"](
            self.query_reads, index_database=self.index_database
        )

        (result2,) = self.plugin.methods["minimap2_search"](
            self.query_reads, reference_reads=self.ref
        )
        pdt.assert_frame_equal(result1.view(pd.DataFrame), result2.view(pd.DataFrame))

    def test_minimap2_both_ref_and_index(self):
        with self.assertRaisesRegex(ValueError, "Only one.*can be provided.*"):
            minimap2_search(
                self.query_reads,
                reference_reads=self.ref,
                index_database=self.index_database,
            )
        with self.assertRaisesRegex(ValueError, "Either.*must be provided.*"):
            minimap2_search(self.query_reads)


if __name__ == "__main__":
    unittest.main()
