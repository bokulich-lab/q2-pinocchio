# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


import pandas as pd
import pandas.testing as pdt
from q2_types.feature_data import FeatureData, Sequence, Taxonomy
from qiime2 import Artifact

from q2_minimap2.classify_consensus import _PAF_format_df_to_series_of_lists

from .test_minimap2 import Minimap2TestsBase


def series_is_subset(expected, observed):
    joined = pd.concat([expected, observed], axis=1, join="inner")
    compared = joined.apply(lambda x: x[0].startswith(x[1]), axis=1)
    return len(compared[compared]) >= len(compared[~compared])


class TestConsensusAssignment(Minimap2TestsBase):
    def setUp(self):
        super().setUp()

        self.taxonomy_tsv_path = self.get_data_path("consensus/taxonomy.tsv")
        self.reads_path = self.get_data_path("consensus/se-dna-sequences.fasta")
        self.paf = Artifact.load(self.get_data_path("search_results.qza"))

    def test_classify_consensus_minimap2(self):
        taxonomy_input = Artifact.import_data(
            FeatureData[Taxonomy], self.taxonomy_tsv_path
        )
        reads = Artifact.import_data(FeatureData[Sequence], self.reads_path)

        (paf, taxonomy) = self.plugin.pipelines["classify_consensus_minimap2"](
            query=reads,
            reference_reads=reads,
            reference_taxonomy=taxonomy_input,
        )

        expected = taxonomy_input.view(pd.Series)
        self.assertTrue(series_is_subset(expected, taxonomy.view(pd.Series)))

    def test_PAF_format_df_to_series_of_lists(self):
        paf = self.paf.view(pd.DataFrame)
        taxonomy_input = Artifact.import_data(
            FeatureData[Taxonomy], self.taxonomy_tsv_path
        )
        taxonomy = taxonomy_input.view(pd.Series)

        obs = _PAF_format_df_to_series_of_lists(paf, taxonomy)
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
    def test_PAF_format_df_to_series_of_lists_fail_on_missing_ids(
        self,
    ):
        paf = self.paf.view(pd.DataFrame)
        taxonomy_input = Artifact.import_data(
            FeatureData[Taxonomy], self.taxonomy_tsv_path
        )
        taxonomy = taxonomy_input.view(pd.Series)
        paf.loc[3] = ["junk", "", "", "", "", "lost-id"] + [""] * 17
        with self.assertRaisesRegex(KeyError, "results do not match.*lost-id"):
            _PAF_format_df_to_series_of_lists(paf, taxonomy)

    def test_find_consensus_annotation(self):
        taxonomy_input = Artifact.import_data(
            FeatureData[Taxonomy], self.taxonomy_tsv_path
        )
        (consensus,) = self.plugin.methods["_find_consensus_annotation"](
            search_results=self.paf, reference_taxonomy=taxonomy_input
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
