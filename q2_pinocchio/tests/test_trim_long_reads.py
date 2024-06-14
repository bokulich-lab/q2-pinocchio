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
import subprocess
import unittest
from unittest.mock import patch

from q2_types.per_sample_sequences import CasavaOneEightSingleLanePerSampleDirFmt

from q2_pinocchio.tests.test_pinocchio import PinocchioTestsBase
from q2_pinocchio.trim_long_reads import (
    construct_chopper_command,
    process_and_rezip,
    trim,
)

seq_ids_maxlen10000 = [
    "ONT1_3",
    "ONT1_2",
    "ONT1_5",
    "ONT1_4",
    "ONT1_7",
]

seq_ids_minlen10000 = [
    "ONT1_1",
    "ONT1_6",
    "ONT1_8",
    "ONT1_9",
    "ONT1_10",
]

seq_ids_pe_minq_20 = [
    "@SRR10138346.5",
    "@SRR10138346.15",
    "@SRR10138346.26",
    "@SRR10138346.10",
    "@SRR10138346.23",
    "@SRR10138346.5",
    "@SRR10138346.10",
    "@SRR10138346.7",
    "@SRR10138346.6",
]


class TestChopperUtilities(PinocchioTestsBase):
    def test_construct_chopper_command(self):
        expected_command = [
            "chopper",
            "--quality",
            "20",
            "--maxqual",
            "35",
            "--minlength",
            "100",
            "--maxlength",
            "500",
            "--headcrop",
            "5",
            "--tailcrop",
            "5",
            "--threads",
            "4",
        ]
        result = construct_chopper_command(
            quality=20,
            maxqual=35,
            minlength=100,
            maxlength=500,
            headcrop=5,
            tailcrop=5,
            threads=4,
        )
        self.assertEqual(result, expected_command)


class TestTrim(PinocchioTestsBase):
    def setUp(self):
        super().setUp()

        self.source_dir_se = self.get_data_path("trim/single_end/")
        self.source_dir_pe = self.get_data_path("trim/paired_end/")

    def test_trimmed_se_maxlen_10000(self):
        # Initialize the query reads from the temporary directory
        query_reads = CasavaOneEightSingleLanePerSampleDirFmt(
            self.source_dir_se, mode="r"
        )

        trimmed_se_maxlen_10000 = trim(query_reads, maxlength=10000)
        fastq_files = [
            os.path.join(str(trimmed_se_maxlen_10000), f)
            for f in os.listdir(str(trimmed_se_maxlen_10000))
            if f.endswith(".fastq.gz")
        ]

        # Process each FASTQ.GZ file
        for file_path in fastq_files:
            with gzip.open(file_path, "rt") as obs_fh:
                self.assertNotEqual(len(obs_fh.readlines()), 0)
                obs_fh.seek(0)
                # Iterate over expected and observed reads, side-by-side
                for records in itertools.zip_longest(*[obs_fh] * 4):
                    (obs_seq_h, obs_seq, _, obs_qual) = records
                    # Make sure seqs that do not map to genome were removed
                    obs_id = obs_seq_h.strip("@\n")
                    self.assertTrue(obs_id in seq_ids_maxlen10000)

    def test_trimmed_se_minlen_10000(self):
        query_reads = CasavaOneEightSingleLanePerSampleDirFmt(
            self.source_dir_se, mode="r"
        )
        trimmed_se_minlen_10000 = trim(query_reads, minlength=10000)

        fastq_files = [
            os.path.join(str(trimmed_se_minlen_10000), f)
            for f in os.listdir(str(trimmed_se_minlen_10000))
            if f.endswith(".fastq.gz")
        ]

        # Process each FASTQ.GZ file
        for file_path in fastq_files:
            with gzip.open(file_path, "rt") as obs_fh:
                self.assertNotEqual(len(obs_fh.readlines()), 0)
                obs_fh.seek(0)
                # Iterate over expected and observed reads, side-by-side
                for records in itertools.zip_longest(*[obs_fh] * 4):
                    (obs_seq_h, obs_seq, _, obs_qual) = records
                    # Make sure seqs that do not map to genome were removed
                    obs_id = obs_seq_h.strip("@\n")
                    self.assertTrue(obs_id in seq_ids_minlen10000)

    def test_trimmed_pe_minq_20(self):
        # Initialize the query reads from the temporary directory
        query_reads = CasavaOneEightSingleLanePerSampleDirFmt(
            self.source_dir_pe, mode="r"
        )
        trimmed_pe_minlen_10000 = trim(query_reads, quality=20)

        fastq_files = [
            os.path.join(str(trimmed_pe_minlen_10000), f)
            for f in os.listdir(str(trimmed_pe_minlen_10000))
            if f.endswith(".fastq.gz")
        ]

        # Process each FASTQ.GZ file
        for file_path in fastq_files:
            with gzip.open(file_path, "rt") as obs_fh:
                # Iterate over expected and observed reads, side-by-side
                for records in itertools.zip_longest(*[obs_fh] * 4):
                    (obs_seq_h, obs_seq, _, obs_qual) = records
                    # Make sure seqs that do not map to genome were removed
                    obs_id = obs_seq_h.split()[0]
                    self.assertTrue(obs_id in seq_ids_pe_minq_20)


class TestConstructChopperCommand(PinocchioTestsBase):
    def test_construct_chopper_command(self):
        # Define the inputs
        quality = 20
        maxqual = 40
        minlength = 100
        maxlength = 1000
        headcrop = 10
        tailcrop = 20
        threads = 4

        # Expected output
        expected_command = [
            "chopper",
            "--quality",
            "20",
            "--maxqual",
            "40",
            "--minlength",
            "100",
            "--maxlength",
            "1000",
            "--headcrop",
            "10",
            "--tailcrop",
            "20",
            "--threads",
            "4",
        ]

        # Call the function
        result = construct_chopper_command(
            quality, maxqual, minlength, maxlength, headcrop, tailcrop, threads
        )

        # Assert the result matches the expected output
        self.assertEqual(result, expected_command)


class TestProcessAndRezip(PinocchioTestsBase):
    @patch("q2_minimap2.trim_long_reads.run_commands_with_pipe")
    def test_process_and_rezip_success(self, mock_run_commands_with_pipe):
        """Test that process_and_rezip runs successfully."""
        input_file = "/fake/input/file.fastq.gz"
        chopper_cmd = ["chopper", "filter"]
        filtered_seqs_path = "/fake/output/file.fastq.gz"

        # Mock run_commands_with_pipe to simulate successful execution
        mock_run_commands_with_pipe.return_value = None  # Simulate no exception raised

        # Call the function
        process_and_rezip(input_file, chopper_cmd, filtered_seqs_path)

        # Check that run_commands_with_pipe was called with the correct arguments
        mock_run_commands_with_pipe.assert_called_once_with(
            ["gunzip", "-c", str(input_file)],
            chopper_cmd,
            ["gzip"],
            str(filtered_seqs_path),
        )

    @patch("q2_minimap2.trim_long_reads.run_commands_with_pipe")
    def test_process_and_rezip_exception(self, mock_run_commands_with_pipe):
        """Test that process_and_rezip raises an exception when chopper fails."""
        input_file = "/fake/input/file.fastq.gz"
        chopper_cmd = ["chopper", "filter"]
        filtered_seqs_path = "/fake/output/file.fastq.gz"

        # Mock run_commands_with_pipe to raise CalledProcessError
        mock_run_commands_with_pipe.side_effect = subprocess.CalledProcessError(
            returncode=1, cmd="chopper"
        )

        with self.assertRaisesRegex(
            Exception, "An error was encountered while using chopper"
        ):
            process_and_rezip(input_file, chopper_cmd, filtered_seqs_path)


if __name__ == "__main__":
    unittest.main()
