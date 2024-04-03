# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


import tempfile

import pandas as pd
import pkg_resources
from qiime2.plugin.testing import TestPluginBase

from q2_minimap2.types._format import (
    Minimap2IndexDBFmt,
    PairwiseAlignmentMN2DirectoryFormat,
    PairwiseAlignmentMN2Format,
)


class Minimap2TypesTestPluginBase(TestPluginBase):
    package = "q2_minimap2.types.tests"

    def setUp(self):
        super().setUp()
        self.temp_dir = tempfile.TemporaryDirectory(
            prefix="q2-long-reads-qc-test-temp-"
        )

    def tearDown(self):
        self.temp_dir.cleanup()

    def get_data_path(self, filename):
        return pkg_resources.resource_filename(self.package, "data/%s" % filename)


class TestMinimap2IndexDB(Minimap2TypesTestPluginBase):
    def test_Minimap2IndexDB_format_validate_positive(self):
        filepath = self.get_data_path("index.mmi")
        format = Minimap2IndexDBFmt(filepath, mode="r")
        format.validate()


class TestPairwiseAlignmentMN2(Minimap2TypesTestPluginBase):
    def test_Minimap2IndexDB_format_validate_positive(self):
        filepath = self.get_data_path("paf_file.paf")
        format = PairwiseAlignmentMN2Format(filepath, mode="r")
        format.validate()

    def test_PairwiseAlignmentMN2_format_to_dataframe_transformer(self):
        filepath = self.get_data_path("paf_file.paf")
        transformer = self.get_transformer(
            pd.DataFrame, PairwiseAlignmentMN2DirectoryFormat
        )
        df = pd.read_csv(filepath, sep="\t")
        obs = transformer(df)
        self.assertIsInstance(obs, PairwiseAlignmentMN2DirectoryFormat)

    def test_PairwiseAlignmentMN2_to_df_transformer(self):
        filepath = self.get_data_path("")
        transformer = self.get_transformer(
            PairwiseAlignmentMN2DirectoryFormat, pd.DataFrame
        )
        paf = PairwiseAlignmentMN2DirectoryFormat(filepath, mode="r")
        obs = transformer(paf)
        self.assertIsInstance(obs, pd.DataFrame)
