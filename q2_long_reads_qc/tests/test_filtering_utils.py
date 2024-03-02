# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import shutil
import tempfile
import unittest

from q2_long_reads_qc.minimap2._filtering_utils import (
    calculate_identity,
    get_alignment_length,
    process_sam_file,
    set_penalties,
)

from .test_long_reads_qc import LongReadsQCTestsBase


class TestSetPenalties(LongReadsQCTestsBase):

    # Test with all parameters provided to ensure correct option string formation
    def test_all_arguments_provided(self):
        self.assertEqual(
            set_penalties(1, -1, 4, 1), ["-A", "1", "-B", "-1", "-O", "4", "-E", "1"]
        )

    # Test with some parameters as None to check function's handling of optional
    # arguments
    def test_some_arguments_none(self):
        # Expecting options only for mismatch and gap_e since others are None
        self.assertEqual(set_penalties(None, -1, None, 2), ["-B", "-1", "-E", "2"])
        # Expecting options only for match and gap_o since others are None
        self.assertEqual(set_penalties(1, None, 4, None), ["-A", "1", "-O", "4"])

    # Test with all parameters as None to verify it returns an empty list, indicating
    # no options
    def test_all_arguments_none(self):
        self.assertEqual(set_penalties(None, None, None, None), [])


class TestCalculateIdentity(LongReadsQCTestsBase):
    # Test case with NM tag present
    def test_with_nm_tag(self):
        aln = "SOME\tDATA\tNM:i:5"
        total_length = 100
        expected_identity = 0.95  # (100 - 5) / 100
        self.assertAlmostEqual(calculate_identity(aln, total_length), expected_identity)

    # Test case without NM tag, should default to 0 mismatches
    def test_without_nm_tag(self):
        aln = "SOME\tDATA"
        total_length = 100
        expected_identity = 1.0  # (100 - 0) / 100
        self.assertAlmostEqual(calculate_identity(aln, total_length), expected_identity)

    # Test case with NM tag present but 0 mismatches
    def test_with_nm_tag_zero_mismatches(self):
        aln = "SOME\tDATA\tNM:i:0"
        total_length = 100
        expected_identity = 1.0  # (100 - 0) / 100
        self.assertAlmostEqual(calculate_identity(aln, total_length), expected_identity)

    # Test case with total length 0 to check for division by zero handling
    def test_total_length_zero(self):
        aln = "SOME\tDATA\tNM:i:5"
        total_length = 0
        with self.assertRaises(ZeroDivisionError):
            calculate_identity(aln, total_length)

    # Test case with negative mismatches, which shouldn't happen but let's test
    def test_negative_mismatches(self):
        aln = "SOME\tDATA\tNM:i:-5"
        total_length = 100
        expected_identity = 1.05  # (100 - (-5)) / 100, though practically impossible
        self.assertAlmostEqual(calculate_identity(aln, total_length), expected_identity)

    # Tests the function's response to a non-integer mismatch count in the NM tag
    def test_malformed_nm_tag(self):
        aln = "SOME\tDATA\tNM:i:notanumber"
        total_length = 100
        # Expecting a ValueError due to the inability to parse "notanumber" as
        # an integer
        with self.assertRaises(ValueError):
            calculate_identity(aln, total_length)

    # Verifies the function's behavior when multiple NM tags are present in
    # the alignment
    def test_multiple_nm_tags(self):
        aln = "SOME\tDATA\tNM:i:5\tANOTHER\tTAG\tNM:i:10"
        total_length = 100
        # Assuming the function picks the first NM tag
        expected_identity = 0.95  # (100 - 5) / 100
        self.assertAlmostEqual(calculate_identity(aln, total_length), expected_identity)

    # Checks how the function calculates identity with mismatches exceeding the
    # total length
    def test_large_number_of_mismatches(self):
        aln = "SOME\tDATA\tNM:i:150"
        total_length = 100
        # More mismatches than total length, resulting in negative identity
        expected_identity = -0.5  # (100 - 150) / 100
        self.assertAlmostEqual(calculate_identity(aln, total_length), expected_identity)

    # Assesses the function's ability to handle special characters in the
    # alignment string
    def test_special_characters_in_alignment(self):
        aln = "SOME\t$%^&*()\tDATA\tNM:i:5"
        total_length = 100
        expected_identity = 0.95  # (100 - 5) / 100
        self.assertAlmostEqual(calculate_identity(aln, total_length), expected_identity)

    # Evaluates the function's handling of a floating-point number for the total
    # alignment length
    def test_floating_point_total_length(self):
        aln = "SOME\tDATA\tNM:i:5"
        total_length = 100.5
        # Test handling of floating point total lengths
        expected_identity = (100.5 - 5) / 100.5
        self.assertAlmostEqual(calculate_identity(aln, total_length), expected_identity)


class TestGetAlignmentLength(LongReadsQCTestsBase):
    # Test with the CIGAR string indicating no alignment
    def test_no_alignment(self):
        self.assertEqual(get_alignment_length("*"), 0)

    # Test with a CIGAR string that has only matches
    def test_only_matches(self):
        self.assertEqual(get_alignment_length("10M"), 10)

    # Test with a CIGAR string that has only insertions
    def test_only_insertions(self):
        self.assertEqual(get_alignment_length("5I"), 5)

    # Test with a CIGAR string that has only deletions
    def test_only_deletions(self):
        self.assertEqual(get_alignment_length("3D"), 3)

    # Test with a CIGAR string that combines matches, insertions, and deletions
    def test_combination(self):
        self.assertEqual(get_alignment_length("5M2I3D"), 10)

    # Test with a more complex CIGAR string
    def test_complex_cigar(self):
        self.assertEqual(get_alignment_length("6M1I4M1D5M"), 17)

    # Test with large numbers in the CIGAR string
    def test_large_numbers(self):
        self.assertEqual(get_alignment_length("100M50I25D"), 175)

    # Test with a CIGAR string that includes operations not contributing to alignment
    #  length (e.g., 'S' for soft clipping)
    def test_with_non_contributing_operations(self):
        self.assertEqual(get_alignment_length("5M3S2I4D"), 11)

    # Test with an empty CIGAR string
    def test_empty_cigar(self):
        self.assertEqual(get_alignment_length(""), 0)

    # Test with invalid characters in the CIGAR string
    def test_invalid_characters(self):
        self.assertEqual(
            get_alignment_length("10M2X3M"), 13
        )  # Assuming 'X' is treated as invalid and ignored

    # Test CIGAR strings with missing numbers before operation characters
    def test_missing_numbers(self):
        self.assertEqual(
            get_alignment_length("MID"), 0
        )  # Assuming it results in 0 as no valid operations are detected

    # Test handling of leading and trailing non-operation characters
    def test_leading_trailing_non_operations(self):
        self.assertEqual(
            get_alignment_length("A10M3DZ"), 13
        )  # Assuming non-operations are ignored

    # Test sequential operations without numbers in-between
    def test_sequential_operations_without_numbers(self):
        self.assertEqual(get_alignment_length("2M3I2D2M"), 9)

    # Test operations specified with a length of 0
    def test_zero_length_operations(self):
        # Assuming 0-length operations are valid and result in 0 contribution
        self.assertEqual(get_alignment_length("0M10I0D"), 10)

    # Test with non-standard operations that don't contribute to alignment length
    def test_non_standard_operations(self):
        self.assertEqual(
            get_alignment_length("10M3S2H"), 10
        )  # Assuming 'S' and 'H' are ignored for alignment length


class TestProcessSamFile(LongReadsQCTestsBase):
    def setUp(self):
        super().setUp()

        self.samfile = self.get_data_path("initial.sam")
        self.only_primary = self.get_data_path("mapped.sam")
        self.only_unmapped = self.get_data_path("unmapped.sam")

    def test_keep_primary_mapped(self):
        # Create a temporary directory
        with tempfile.TemporaryDirectory() as tmpdir:
            # Concatenate the temporary directory path and the original file name
            tmpfile = os.path.join(tmpdir, "initial_tmp.sam")

            # Copy the SAM file to the temporary location
            shutil.copy(self.samfile, tmpfile)

            process_sam_file(tmpfile, False, None)

            # Check each line for equality using splitlines()
            with open(tmpfile, "r") as tmpfile_content, open(
                self.only_primary, "r"
            ) as primary_content:
                tmp_lines = tmpfile_content.read().splitlines()
                primary_lines = primary_content.read().splitlines()

                tmp_index = 0
                primary_index = 0

                while tmp_index < len(tmp_lines) and primary_index < len(primary_lines):
                    tmp_line = tmp_lines[tmp_index]
                    primary_line = primary_lines[primary_index]

                    if tmp_line.startswith("@PG"):
                        # Skip the line starting with '@PG' in tmp_lines
                        tmp_index += 1
                        primary_index += 1
                        continue

                    # Explicitly print the line
                    self.assertEqual(tmp_line, primary_line)

                    tmp_index += 1
                    primary_index += 1

    def test_keep_unmapped(self):
        # Create a temporary directory
        with tempfile.TemporaryDirectory() as tmpdir:
            # Concatenate the temporary directory path and the original file name
            tmpfile = os.path.join(tmpdir, "initial_tmp.sam")

            # Copy the SAM file to the temporary location
            shutil.copy(self.samfile, tmpfile)

            process_sam_file(tmpfile, True, None)

            # Check each line for equality using splitlines()
            with open(tmpfile, "r") as tmpfile_content, open(
                self.only_unmapped, "r"
            ) as unmapped_content:
                tmp_lines = tmpfile_content.read().splitlines()
                unmapped_lines = unmapped_content.read().splitlines()

                tmp_index = 0
                unmapped_index = 0

                while tmp_index < len(tmp_lines) and unmapped_index < len(
                    unmapped_lines
                ):
                    tmp_line = tmp_lines[tmp_index]
                    unmapped_line = unmapped_lines[unmapped_index]

                    if tmp_line.startswith("@PG"):
                        # Skip the line starting with '@PG' in tmp_lines
                        tmp_index += 1
                        unmapped_index += 1
                        continue

                    # Explicitly print the line
                    self.assertEqual(tmp_line, unmapped_line)

                    tmp_index += 1
                    unmapped_index += 1


if __name__ == "__main__":
    unittest.main()
