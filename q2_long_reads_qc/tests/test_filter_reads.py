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

from q2_types.per_sample_sequences import (
    FastqGzFormat,
    SingleLanePerSampleSingleEndFastqDirFmt,
)
from qiime2 import Artifact

from .test_long_reads_qc import LongReadsQCTestsBase

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


class TestFilterReads(LongReadsQCTestsBase):
    def setUp(self):
        super().setUp()

        self.query_reads = Artifact.load(self.get_data_path("single-end.qza"))
        self.minimap2_index = Artifact.load(self.get_data_path("index.qza"))

    # Exclude mapped
    def test_filter_exclude_mapped(self):
        (obs_art,) = self.plugin.methods["filter_reads"](
            self.query_reads, self.minimap2_index, exclude_mapped=True
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
    def test_filter_keep_mapped(self):
        (obs_art,) = self.plugin.methods["filter_reads"](
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

    def test_exclude_mapped_with_perc_id(self):
        (obs_art,) = self.plugin.methods["filter_reads"](
            self.query_reads,
            self.minimap2_index,
            exclude_mapped=True,
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

    def test_include_unmapped_with_perc_id(self):
        (obs_art,) = self.plugin.methods["filter_reads"](
            self.query_reads,
            self.minimap2_index,
            exclude_mapped=False,
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


if __name__ == "__main__":
    unittest.main()
