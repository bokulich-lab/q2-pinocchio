# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import gzip
import itertools
import os
import unittest

from q2_types.feature_data import DNAFASTAFormat
from q2_types.per_sample_sequences import CasavaOneEightSingleLanePerSampleDirFmt

from q2_minimap2.filter_reads import filter_reads
from q2_minimap2.types._format import Minimap2IndexDBDirFmt

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

        self.query_single_reads = CasavaOneEightSingleLanePerSampleDirFmt(
            self.get_data_path("filter_reads/single_end/"), mode="r"
        )
        self.query_paired_reads = CasavaOneEightSingleLanePerSampleDirFmt(
            self.get_data_path("filter_reads/paired_end/"), mode="r"
        )
        self.minimap2_index = Minimap2IndexDBDirFmt(
            self.get_data_path("filter_reads/index/"), mode="r"
        )
        self.reference_reads = DNAFASTAFormat(
            self.get_data_path("filter_reads/dna-sequences.fasta"), mode="r"
        )

    def _check_ids(self, obs_seqs, included_ids, excluded_ids):

        fastq_files = [f for f in os.listdir(str(obs_seqs)) if f.endswith(".fastq.gz")]

        # Process each FASTQ.GZ file
        for obs_fp in fastq_files:
            file_path = os.path.join(str(obs_seqs), obs_fp)
            with gzip.open(file_path, "rt") as obs_fh:
                # Ensure the file is not empty
                self.assertNotEqual(len(obs_fh.readlines()), 0)
                obs_fh.seek(0)
                # Iterate over expected and observed reads, side-by-side
                for records in itertools.zip_longest(*[obs_fh] * 4):
                    (obs_seq_h, obs_seq, _, obs_qual) = records
                    # Make sure seqs that map to genome were removed
                    obs_id = obs_seq_h.strip("@/012\n")
                    self.assertTrue(obs_id in included_ids)
                    self.assertTrue(obs_id not in excluded_ids)

    def test_filter_single_end_keep_unmapped(self):
        obs_seqs = filter_reads(
            query_reads=self.query_single_reads,
            index_database=self.minimap2_index,
            keep="unmapped",
        )

        self._check_ids(obs_seqs, seq_ids_unmapped, seq_ids_mapped)

    def test_filter_single_end_keep_mapped(self):
        obs_seqs = filter_reads(
            query_reads=self.query_single_reads,
            index_database=self.minimap2_index,
        )
        self._check_ids(obs_seqs, seq_ids_mapped, seq_ids_unmapped)

    def test_filter_single_end_keep_mapped_sr(self):
        (obs_art,) = self.plugin.methods["filter_reads"](
            self.query_reads_single, self.minimap2_index, mapping_preset="sr"
        )

        obs = obs_art.view(SingleLanePerSampleSingleEndFastqDirFmt)
        obs_seqs = obs.sequences.iter_views(FastqGzFormat)
        self._check_ids(obs_seqs, seq_ids_mapped, seq_ids_unmapped)

    def test_filter_single_end_keep_mapped_using_ref(self):
        obs_seqs = filter_reads(
            query_reads=self.query_single_reads,
            reference_reads=self.reference_reads,
        )
        self._check_ids(obs_seqs, seq_ids_mapped, seq_ids_unmapped)

    def test_filter_single_end_keep_unmapped_with_perc_id(self):
        obs_seqs = filter_reads(
            query_reads=self.query_single_reads,
            index_database=self.minimap2_index,
            keep="unmapped",
            min_per_identity=0.99,
        )
        self._check_ids(obs_seqs, perc_id_unmapped, perc_id_mapped)

    def test_filter_single_end_keep_mapped_with_perc_id(self):
        obs_seqs = filter_reads(
            query_reads=self.query_single_reads,
            index_database=self.minimap2_index,
            keep="mapped",
            min_per_identity=0.99,
        )
        self._check_ids(obs_seqs, perc_id_mapped, perc_id_unmapped)

    def test_both_reference_and_index_provided(self):
        with self.assertRaises(ValueError) as context:
            filter_reads(
                query_reads=self.query_single_reads,
                index_database=self.minimap2_index,
                reference_reads=self.reference_reads,
            )
        self.assertIn(
            "Only one of reference_reads or index_database can be provided",
            str(context.exception),
        )

    def test_neither_reference_nor_index_provided(self):
        with self.assertRaises(ValueError) as context:
            filter_reads(
                query_reads=self.query_single_reads,
                index_database=None,
                reference_reads=None,
            )
        self.assertIn(
            "Either reference_reads or index_database must be provided",
            str(context.exception),
        )

    def test_filter_paired_end_keep_unmapped(self):
        obs_seqs = filter_reads(
            query_reads=self.query_paired_reads,
            index_database=self.minimap2_index,
            keep="unmapped",
        )
        self._check_ids(obs_seqs, seq_ids_unmapped, seq_ids_mapped)

    def test_filter_paired_end_keep_mapped(self):
        obs_seqs = filter_reads(
            query_reads=self.query_paired_reads,
            index_database=self.minimap2_index,
        )
        self._check_ids(obs_seqs, seq_ids_mapped, seq_ids_unmapped)


if __name__ == "__main__":
    unittest.main()
