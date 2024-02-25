# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from qiime2.plugin.testing import TestPluginBase


class LongReadsQCTestsBase(TestPluginBase):
    package = "q2_long_reads_qc.tests"

    def setUp(self):
        super().setUp()
