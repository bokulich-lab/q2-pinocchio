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
from distutils.dir_util import copy_tree

import pkg_resources
import q2templates
from q2_types.per_sample_sequences import CasavaOneEightSingleLanePerSampleDirFmt

from q2_minimap2._utils import run_command


# Run NanoPlot on sequence files in the specified directory
def _run_nanoplot(sequences_path, nanoplot_output):
    fastq_files = [
        os.path.join(sequences_path, f)
        for f in os.listdir(sequences_path)
        if f.endswith(".fastq.gz")
    ]
    nanoplot_cmd = ["NanoPlot", "--fastq", *fastq_files, "-o", nanoplot_output]
    try:
        # Run Nanoplot command
        run_command(nanoplot_cmd, verbose=True)
    except subprocess.CalledProcessError as e:
        raise Exception(
            "An error was encountered while running nanoplot, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )


def _create_visualization(output_dir, nanoplot_output):
    # Copy Nanoplot templates to the output directory
    TEMPLATES = pkg_resources.resource_filename("q2_minimap2", "assets")
    copy_tree(os.path.join(TEMPLATES, "nanoplot"), output_dir)

    # Copy Nanoplot data from the temporary directory to the output directory
    copy_tree(nanoplot_output, os.path.join(output_dir, "nanoplot_data"))

    # Generate an index.html file for Nanoplot in the output directory
    context = {"tabs": [{"title": "Nanoplot", "url": "index.html"}]}
    index = os.path.join(TEMPLATES, "nanoplot", "index.html")
    q2templates.render([index], output_dir, context=context)


def stats(
    output_dir: str,
    sequences: CasavaOneEightSingleLanePerSampleDirFmt,
):
    with tempfile.TemporaryDirectory() as nanoplot_output:
        _run_nanoplot(sequences.path, nanoplot_output)
        _create_visualization(output_dir, nanoplot_output)
