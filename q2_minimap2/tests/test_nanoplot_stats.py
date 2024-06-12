# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from q2_types.per_sample_sequences import CasavaOneEightSingleLanePerSampleDirFmt

from q2_minimap2.nanoplot_stats import _create_visualization, _run_nanoplot, stats
from q2_minimap2.tests.test_minimap2 import Minimap2TestsBase


class TestRunNanoPlot(Minimap2TestsBase):
    @patch("q2_minimap2.nanoplot_stats.run_command")
    def test_run_nanoplot_success(self, mock_run_command):
        """Test that _run_nanoplot runs successfully."""
        sequences_path = self.get_data_path("nanoplot/")
        with tempfile.TemporaryDirectory() as output_dir:

            # Mock run_command to simulate successful execution
            mock_run_command.return_value = None  # Simulate no exception raised

            # Call the function
            _run_nanoplot(sequences_path, output_dir)

            # Check that run_command was called correctly
            fastq_files = [
                os.path.join(sequences_path, f)
                for f in os.listdir(sequences_path)
                if f.endswith(".fastq.gz")
            ]
            nanoplot_cmd = ["NanoPlot", "--fastq", *fastq_files, "-o", output_dir]
            mock_run_command.assert_called_once_with(nanoplot_cmd, verbose=True)

    @patch("q2_minimap2.nanoplot_stats.run_command")
    def test_run_nanoplot_exception(self, mock_run_command):
        """Test that _run_nanoplot raises an exception when NanoPlot fails."""
        sequences_path = self.get_data_path("nanoplot/")
        with tempfile.TemporaryDirectory() as output_dir:
            # Mock run_command to raise CalledProcessError
            mock_run_command.side_effect = subprocess.CalledProcessError(
                returncode=1, cmd="NanoPlot"
            )
            with self.assertRaises(Exception) as context:
                _run_nanoplot(sequences_path, output_dir)
            self.assertIn(
                "An error was encountered while running nanoplot",
                str(context.exception),
            )


class TestCreateVisualization(Minimap2TestsBase):
    @patch("q2_minimap2.nanoplot_stats.q2templates.render")
    @patch("q2_minimap2.nanoplot_stats.copy_tree")
    @patch("q2_minimap2.nanoplot_stats.pkg_resources.resource_filename")
    def test_create_visualization(
        self, mock_resource_filename, mock_copy_tree, mock_render
    ):
        """Test that copies templates and data, and renders the index.html."""
        output_dir = "/fake/output/dir"
        nanoplot_output = "/fake/nanoplot/output"

        # Mock the resource filename to return a fake templates directory
        mock_resource_filename.return_value = "/fake/templates/dir"

        # Call the function
        _create_visualization(output_dir, nanoplot_output)

        # Check that resource_filename was called correctly
        mock_resource_filename.assert_called_once_with("q2_minimap2", "assets")

        # Check that copy_tree was called correctly for templates and nanoplot data
        mock_copy_tree.assert_any_call("/fake/templates/dir/nanoplot", output_dir)
        mock_copy_tree.assert_any_call(
            nanoplot_output, os.path.join(output_dir, "nanoplot_data")
        )

        # Check that q2templates.render was called correctly
        expected_context = {"tabs": [{"title": "Nanoplot", "url": "index.html"}]}
        expected_index = os.path.join("/fake/templates/dir/nanoplot", "index.html")
        mock_render.assert_called_once_with(
            [expected_index], output_dir, context=expected_context
        )


class TestStats(Minimap2TestsBase):
    @patch("q2_minimap2.nanoplot_stats._create_visualization")
    @patch("q2_minimap2.nanoplot_stats._run_nanoplot")
    def test_stats(self, mock_run_nanoplot, mock_create_visualization):
        """Test the stats function."""
        # Create the necessary objects
        output_dir = "/fake/output/dir"
        sequences = CasavaOneEightSingleLanePerSampleDirFmt()
        sequences.path = "/fake/sequences/path"  # Set a fake path attribute

        # Call the stats function
        stats(output_dir, sequences)

        # Check that _run_nanoplot was called with the correct first argument
        args, _ = mock_run_nanoplot.call_args
        self.assertEqual(args[0], sequences.path)

        # Check that _create_visualization was called with the correct first argument
        args, _ = mock_create_visualization.call_args
        self.assertEqual(args[0], output_dir)


if __name__ == "__main__":
    unittest.main()
