# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from qiime2.plugin.testing import TestPluginBase


class Minimap2TestsBase(TestPluginBase):
    package = "q2_minimap2.tests"

    def setUp(self):
        super().setUp()
