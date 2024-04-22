# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import gzip
import itertools
import unittest

from q2_types.feature_data import FeatureData, Sequence
from q2_types.per_sample_sequences import (
    FastqGzFormat,
    SequencesWithQuality,
    SingleLanePerSampleSingleEndFastqDirFmt,
)
from q2_types.sample_data import SampleData
from qiime2 import Artifact

from q2_minimap2.types._type import Minimap2IndexDB

from .test_minimap2 import Minimap2TestsBase

seq_ids_mapped = [
    "SARS2:6:73:941:1973#",
    "SARS2:6:73:231:3321#",
    "SARS2:6:73:233:3421#",
    "SARS2:6:73:552:2457#",
    "SARS2:6:73:567:7631#",
]
seq_ids_unmapped = ["SARS2:6:73:356:9806#"]

perc_id_mapped = [
    "SARS2:6:73:231:3321#",
    "SARS2:6:73:233:3421#",
    "SARS2:6:73:552:2457#",
    "SARS2:6:73:567:7631#",
]

perc_id_unmapped = ["SARS2:6:73:941:1973#", "SARS2:6:73:356:9806#"]


class TestFilterSingleEndReads(Minimap2TestsBase):
    def setUp(self):
        super().setUp()

        minimap2_index_path = self.get_data_path("filter_reads/index.mmi")
        query_reads_path = self.get_data_path("filter_reads/single_end/")
        reference_reads_path = self.get_data_path("filter_reads/dna-sequences.fasta")

        self.query_reads = Artifact.import_data(
            SampleData[SequencesWithQuality], query_reads_path
        )
        self.minimap2_index = Artifact.import_data(Minimap2IndexDB, minimap2_index_path)
        self.reference_reads = Artifact.import_data(
            FeatureData[Sequence], reference_reads_path
        )

    # Keep unmapped
    def test_filter_single_end_keep_unmapped(self):
        (obs_art,) = self.plugin.methods["filter_single_end_reads"](
            self.query_reads, self.minimap2_index, keep="unmapped"
        )

        obs = obs_art.view(SingleLanePerSampleSingleEndFastqDirFmt)
        obs_seqs = obs.sequences.iter_views(FastqGzFormat)

        for _, obs_fp in obs_seqs:
            with gzip.open(str(obs_fp), "rt") as obs_fh:
                self.assertNotEqual(len(obs_fh.readlines()), 0)
                obs_fh.seek(0)
                # Iterate over expected and observed reads, side-by-side
                for records in itertools.zip_longest(*[obs_fh] * 4):
                    (obs_seq_h, obs_seq, _, obs_qual) = records
                    # Make sure seqs that map to genome were removed
                    obs_id = obs_seq_h.strip("@/012\n")
                    self.assertTrue(obs_id not in seq_ids_mapped)
                    self.assertTrue(obs_id in seq_ids_unmapped)

    # Keep mapped
    def test_filter_single_end_keep_mapped(self):
        (obs_art,) = self.plugin.methods["filter_single_end_reads"](
            self.query_reads, self.minimap2_index
        )

        obs = obs_art.view(SingleLanePerSampleSingleEndFastqDirFmt)
        obs_seqs = obs.sequences.iter_views(FastqGzFormat)

        for _, obs_fp in obs_seqs:
            with gzip.open(str(obs_fp), "rt") as obs_fh:
                self.assertNotEqual(len(obs_fh.readlines()), 0)
                obs_fh.seek(0)
                # Iterate over expected and observed reads, side-by-side
                for records in itertools.zip_longest(*[obs_fh] * 4):
                    (obs_seq_h, obs_seq, _, obs_qual) = records
                    # Make sure seqs that do not map to genome were removed
                    obs_id = obs_seq_h.strip("@/012\n")
                    self.assertTrue(obs_id in seq_ids_mapped)
                    self.assertTrue(obs_id not in seq_ids_unmapped)

    # Keep mapped using reference
    def test_filter_single_end_keep_mapped_using_ref(self):
        (obs_art,) = self.plugin.methods["filter_single_end_reads"](
            self.query_reads, reference_reads=self.reference_reads
        )

        obs = obs_art.view(SingleLanePerSampleSingleEndFastqDirFmt)
        obs_seqs = obs.sequences.iter_views(FastqGzFormat)

        for _, obs_fp in obs_seqs:
            with gzip.open(str(obs_fp), "rt") as obs_fh:
                self.assertNotEqual(len(obs_fh.readlines()), 0)
                obs_fh.seek(0)
                # Iterate over expected and observed reads, side-by-side
                for records in itertools.zip_longest(*[obs_fh] * 4):
                    (obs_seq_h, obs_seq, _, obs_qual) = records
                    # Make sure seqs that do not map to genome were removed
                    obs_id = obs_seq_h.strip("@/012\n")
                    self.assertTrue(obs_id in seq_ids_mapped)
                    self.assertTrue(obs_id not in seq_ids_unmapped)

    def test_filter_single_end_keep_unmapped_with_perc_id(self):
        (obs_art,) = self.plugin.methods["filter_single_end_reads"](
            self.query_reads,
            self.minimap2_index,
            keep="unmapped",
            min_per_identity=0.99,
        )

        obs = obs_art.view(SingleLanePerSampleSingleEndFastqDirFmt)
        obs_seqs = obs.sequences.iter_views(FastqGzFormat)

        for _, obs_fp in obs_seqs:
            with gzip.open(str(obs_fp), "rt") as obs_fh:
                self.assertNotEqual(len(obs_fh.readlines()), 0)
                obs_fh.seek(0)
                # Iterate over expected and observed reads, side-by-side
                for records in itertools.zip_longest(*[obs_fh] * 4):
                    (obs_seq_h, obs_seq, _, obs_qual) = records
                    # Make sure seqs that do not map to genome were removed
                    obs_id = obs_seq_h.strip("@/012\n")
                    self.assertTrue(obs_id in perc_id_unmapped)
                    self.assertTrue(obs_id not in perc_id_mapped)

    def test_filter_single_end_keep_mapped_with_perc_id(self):
        (obs_art,) = self.plugin.methods["filter_single_end_reads"](
            self.query_reads,
            self.minimap2_index,
            keep="mapped",
            min_per_identity=0.99,
        )

        obs = obs_art.view(SingleLanePerSampleSingleEndFastqDirFmt)
        obs_seqs = obs.sequences.iter_views(FastqGzFormat)

        for _, obs_fp in obs_seqs:
            with gzip.open(str(obs_fp), "rt") as obs_fh:
                self.assertNotEqual(len(obs_fh.readlines()), 0)
                obs_fh.seek(0)
                # Iterate over expected and observed reads, side-by-side
                for records in itertools.zip_longest(*[obs_fh] * 4):
                    (obs_seq_h, obs_seq, _, obs_qual) = records
                    # Make sure seqs that do not map to genome were removed
                    obs_id = obs_seq_h.strip("@/012\n")
                    self.assertTrue(obs_id in perc_id_mapped)
                    self.assertTrue(obs_id not in perc_id_unmapped)

    def test_both_reference_and_index_provided(self):
        with self.assertRaises(ValueError) as context:
            self.plugin.methods["filter_single_end_reads"](
                self.query_reads,
                index_database=self.minimap2_index,
                reference_reads=self.reference_reads,
            )
        self.assertIn(
            "Only one of reference_reads or index_database can be provided",
            str(context.exception),
        )

    def test_neither_reference_nor_index_provided(self):
        with self.assertRaises(ValueError) as context:
            self.plugin.methods["filter_single_end_reads"](
                self.query_reads, index_database=None, reference_reads=None
            )
        self.assertIn(
            "Either reference_reads or index_database must be provided",
            str(context.exception),
        )


if __name__ == "__main__":
    unittest.main()
