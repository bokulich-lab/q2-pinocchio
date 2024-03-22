# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import itertools
import os
import tempfile
import unittest

from qiime2 import Artifact

from q2_minimap2.extract_seqs import extract_seqs

from .test_long_reads import Minimap2TestsBase

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


class TestExtractSeqs(Minimap2TestsBase):
    def setUp(self):
        super().setUp()

        self.query_reads = Artifact.load(
            self.get_data_path("extract_seqs_test_input.qza")
        )
        self.minimap2_index = Artifact.load(self.get_data_path("index.qza"))

        self.query2 = Artifact.load(self.get_data_path("query-seqs.qza"))
        self.index_database2 = Artifact.load(
            self.get_data_path("minimap2_test_index.qza")
        )
        self.ref = Artifact.load(self.get_data_path("se-dna-sequences.qza"))

    # Exclude mapped
    def test_extract_unmapped(self):
        (obs_art,) = self.plugin.methods["extract_seqs"](
            self.query_reads, self.minimap2_index, extract="unmapped"
        )

        # Create a temporary directory to hold the extracted sequences
        with tempfile.TemporaryDirectory() as temp_dir:
            obs_art.export_data(temp_dir)
            obs_fp = os.path.join(temp_dir, "dna-sequences.fasta")

            # Open and read the FASTA file
            with open(obs_fp, "rt") as obs_fh:
                # Ensure the file is not empty
                self.assertNotEqual(len(obs_fh.readlines()), 0)
                obs_fh.seek(0)

                # Iterate over expected and observed reads, side-by-side
                for records in itertools.zip_longest(*[obs_fh] * 4):
                    (obs_seq_h, obs_seq, _, obs_qual) = records
                    # Extract sequence ID from header
                    obs_id = obs_seq_h.strip(">/012\n")
                    # Ensure sequences that map to the genome were removed
                    self.assertTrue(obs_id not in seq_ids_mapped)
                    self.assertTrue(obs_id in seq_ids_unmapped)

    # Keep mapped
    def test_extract_mapped(self):
        (obs_art,) = self.plugin.methods["extract_seqs"](
            self.query_reads, self.minimap2_index
        )

        # Create a temporary directory to hold the extracted sequences
        with tempfile.TemporaryDirectory() as temp_dir:
            obs_art.export_data(temp_dir)
            obs_fp = os.path.join(temp_dir, "dna-sequences.fasta")

            # Open and read the FASTA file
            with open(obs_fp, "rt") as obs_fh:
                # Ensure the file is not empty
                self.assertNotEqual(len(obs_fh.readlines()), 0)
                obs_fh.seek(0)
                # Iterate over expected and observed reads, side-by-side
                for records in itertools.zip_longest(*[obs_fh] * 4):
                    (obs_seq_h, obs_seq, _, obs_qual) = records
                    # Make sure seqs that do not map to genome were removed
                    obs_id = obs_seq_h.strip(">/012\n")
                    self.assertTrue(obs_id in seq_ids_mapped)
                    self.assertTrue(obs_id not in seq_ids_unmapped)

    def test_extract_unmapped_with_perc_id(self):
        (obs_art,) = self.plugin.methods["extract_seqs"](
            self.query_reads,
            self.minimap2_index,
            extract="unmapped",
            min_per_identity=0.99,
        )

        # Create a temporary directory to hold the extracted sequences
        with tempfile.TemporaryDirectory() as temp_dir:
            obs_art.export_data(temp_dir)
            obs_fp = os.path.join(temp_dir, "dna-sequences.fasta")

            # Open and read the FASTA file
            with open(obs_fp, "rt") as obs_fh:
                # Ensure the file is not empty
                self.assertNotEqual(len(obs_fh.readlines()), 0)
                obs_fh.seek(0)
                # Iterate over expected and observed reads, side-by-side
                for records in itertools.zip_longest(*[obs_fh] * 4):
                    (obs_seq_h, obs_seq, _, obs_qual) = records
                    # Make sure seqs that do not map to genome were removed
                    obs_id = obs_seq_h.strip(">/012\n")
                    self.assertTrue(obs_id in perc_id_unmapped)
                    self.assertTrue(obs_id not in perc_id_mapped)

    def test_extract_mapped_with_perc_id(self):
        (obs_art,) = self.plugin.methods["extract_seqs"](
            self.query_reads,
            self.minimap2_index,
            extract="mapped",
            min_per_identity=0.99,
        )

        # Create a temporary directory to hold the extracted sequences
        with tempfile.TemporaryDirectory() as temp_dir:
            obs_art.export_data(temp_dir)
            obs_fp = os.path.join(temp_dir, "dna-sequences.fasta")

            # Open and read the FASTA file
            with open(obs_fp, "rt") as obs_fh:
                # Ensure the file is not empty
                self.assertNotEqual(len(obs_fh.readlines()), 0)
                obs_fh.seek(0)
                # Iterate over expected and observed reads, side-by-side
                for records in itertools.zip_longest(*[obs_fh] * 4):
                    (obs_seq_h, obs_seq, _, obs_qual) = records
                    # Make sure seqs that do not map to genome were removed
                    obs_id = obs_seq_h.strip(">/012\n")
                    self.assertTrue(obs_id in perc_id_mapped)
                    self.assertTrue(obs_id not in perc_id_unmapped)

    def test_extract_mapped_using_reference(self):
        (obs_art,) = self.plugin.methods["extract_seqs"](
            self.query2, reference_reads=self.ref
        )

        # Create a temporary directory to hold the extracted sequences
        with tempfile.TemporaryDirectory() as temp_dir:
            obs_art.export_data(temp_dir)
            obs_fp = os.path.join(temp_dir, "dna-sequences.fasta")
            correct_output_fp = self.get_data_path("extracted_mapped.fasta")

            with open(correct_output_fp, "r") as file1, open(obs_fp, "r") as file2:
                true_fasta_content = file1.read()
                output_fasta_content = file2.read()
                self.assertEqual(true_fasta_content, output_fasta_content)

    def test_extract_seqs_both_ref_and_index(self):
        with self.assertRaisesRegex(ValueError, "Only one.*can be provided.*"):
            extract_seqs(
                self.query2,
                reference_reads=self.ref,
                index_database=self.index_database2,
            )
        with self.assertRaisesRegex(ValueError, "Either.*must be provided.*"):
            extract_seqs(self.query2)


if __name__ == "__main__":
    unittest.main()
