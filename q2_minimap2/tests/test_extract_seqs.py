# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import itertools
import unittest

from q2_types.feature_data import DNAFASTAFormat

from q2_minimap2.extract_seqs import extract_seqs
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


class TestExtractSeqs(Minimap2TestsBase):
    def setUp(self):
        super().setUp()
        self.query1_reads = DNAFASTAFormat(
            self.get_data_path("extract_seqs/extract_seqs_test_input.fasta"), mode="r"
        )
        self.query2_reads = DNAFASTAFormat(
            self.get_data_path("extract_seqs/query-seqs.fasta"), mode="r"
        )
        self.minimap2_index = Minimap2IndexDBDirFmt(
            self.get_data_path("extract_seqs/index/"), mode="r"
        )
        self.ref = DNAFASTAFormat(
            self.get_data_path("extract_seqs/se-dna-sequences.fasta"), mode="r"
        )

    # Exclude mapped
    def test_extract_unmapped(self):
        extracted_sequences = extract_seqs(
            self.query1_reads, self.minimap2_index, extract="unmapped"
        )
        # Open and read the FASTA file
        with open(str(extracted_sequences), "rt") as obs_fh:
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
        extracted_sequences = extract_seqs(self.query1_reads, self.minimap2_index)
        # Open and read the FASTA file
        with open(str(extracted_sequences), "rt") as obs_fh:
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
        extracted_sequences = extract_seqs(
            self.query1_reads,
            self.minimap2_index,
            extract="unmapped",
            min_per_identity=0.99,
        )
        # Open and read the FASTA file
        with open(str(extracted_sequences), "rt") as obs_fh:
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
        extracted_sequences = extract_seqs(
            self.query1_reads,
            self.minimap2_index,
            extract="mapped",
            min_per_identity=0.99,
        )
        # Open and read the FASTA file
        with open(str(extracted_sequences), "rt") as obs_fh:
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
        extracted_sequences = extract_seqs(self.query2_reads, reference_reads=self.ref)
        obs_fp = str(extracted_sequences)
        correct_output_fp = self.get_data_path("extract_seqs/extracted_mapped.fasta")
        with open(correct_output_fp, "r") as file1, open(obs_fp, "r") as file2:
            true_fasta_content = file1.read()
            output_fasta_content = file2.read()
            self.assertEqual(true_fasta_content, output_fasta_content)

    def test_extract_seqs_both_ref_and_index(self):
        with self.assertRaisesRegex(ValueError, "Only one.*can be provided.*"):
            extract_seqs(
                self.query2_reads,
                reference_reads=self.ref,
                index_database=self.minimap2_index,
            )
        with self.assertRaisesRegex(ValueError, "Either.*must be provided.*"):
            extract_seqs(self.query2_reads)


if __name__ == "__main__":
    unittest.main()
