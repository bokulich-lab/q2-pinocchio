# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


import pandas as pd
import pandas.testing as pdt
from qiime2 import Artifact

from q2_long_reads_qc.classify_consensus import (
    _PairwiseAlignmentMN2_format_df_to_series_of_lists,
)

from .test_long_reads_qc import LongReadsQCTestsBase


def series_is_subset(expected, observed):
    joined = pd.concat([expected, observed], axis=1, join="inner")
    compared = joined.apply(lambda x: x[0].startswith(x[1]), axis=1)
    return len(compared[compared]) >= len(compared[~compared])


class TestConsensusAssignment(LongReadsQCTestsBase):
    def setUp(self):
        super().setUp()
        self.taxonomy = Artifact.load(self.get_data_path("taxonomy.qza"))
        self.reads = Artifact.load(self.get_data_path("se-dna-sequences.qza"))
        self.paf = Artifact.load(self.get_data_path("search_results.qza"))
        self.exp = self.taxonomy.view(pd.Series)
        self.query_reads = Artifact.load(self.get_data_path("query-seqs.qza"))

    def test_classify_consensus_minimap2(self):
        (paf, taxonomy) = self.plugin.pipelines["classify_consensus_minimap2"](
            query=self.reads,
            reference_reads=self.reads,
            reference_taxonomy=self.taxonomy,
        )

        self.assertTrue(series_is_subset(self.exp, taxonomy.view(pd.Series)))

    def test_PairwiseAlignmentMN2_format_df_to_series_of_lists(self):
        paf = self.paf.view(pd.DataFrame)
        taxonomy = self.taxonomy.view(pd.Series)

        obs = _PairwiseAlignmentMN2_format_df_to_series_of_lists(paf, taxonomy)
        exp = pd.Series(
            {
                "1111561": [
                    "k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; "
                    "o__Legionellales; f__; g__; s__",
                    "k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; "
                    "o__; f__; g__; s__",
                ],
                "835097": [
                    "k__Bacteria; p__Chloroflexi; c__SAR202; o__; f__; g__; s__"
                ],
                "junk": ["Unassigned"],
            }
        )

        # Use assert to compare the sorted lists
        assert exp.tolist() == obs.tolist(), "Contents of the Series are not equal."

    # should fail when hit IDs are missing from reference taxonomy
    def test_PairwiseAlignmentMN2_format_df_to_series_of_lists_fail_on_missing_ids(
        self,
    ):
        paf = self.paf.view(pd.DataFrame)
        taxonomy = self.taxonomy.view(pd.Series)
        paf.loc[3] = ["junk", "", "", "", "", "lost-id"] + [""] * 17
        with self.assertRaisesRegex(KeyError, "results do not match.*lost-id"):
            _PairwiseAlignmentMN2_format_df_to_series_of_lists(paf, taxonomy)

    def test_find_consensus_annotation(self):
        (consensus,) = self.plugin.methods["find_consensus_annotation"](
            search_results=self.paf, reference_taxonomy=self.taxonomy
        )
        obs = consensus.view(pd.DataFrame)

        exp = pd.DataFrame(
            {
                "Taxon": {
                    "1111561": "k__Bacteria; p__Proteobacteria; "
                    "c__Gammaproteobacteria",
                    "835097": "k__Bacteria; p__Chloroflexi; c__SAR202; o__; f__; "
                    "g__; s__",
                    "junk": "Unassigned",
                },
                "Consensus": {"1111561": "1.0", "835097": "1.0", "junk": "1.0"},
            }
        )
        pdt.assert_frame_equal(exp, obs, check_names=False)
