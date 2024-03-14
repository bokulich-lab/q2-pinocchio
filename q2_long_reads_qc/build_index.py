# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import subprocess

from q2_types.feature_data import DNAFASTAFormat

from q2_long_reads_qc._utils import run_command
from q2_long_reads_qc.types._format import Minimap2IndexDBDirFmt


def build_index(
    sequences: DNAFASTAFormat, kmer_length: int = 15
) -> Minimap2IndexDBDirFmt:
    """
    Generates a QIIME 2 artifact representing a Minimap2 index database.

    Args:
    sequences: Reference sequences to be indexed.
    kmer_length: The k-mer size for indexing.

    Returns:
    A QIIME 2 artifact containing the Minimap2 index.
    """

    # Initialize a directory format object to store the Minimap2 index.
    database = Minimap2IndexDBDirFmt()

    # Command to build the Minimap2 index file
    build_cmd = [
        "minimap2",
        "-d",
        str(database.path / "index.mmi"),
        str(sequences),
        "-k",
        str(kmer_length),
    ]

    try:
        # Execute the command to create the index
        run_command(build_cmd)
    except subprocess.CalledProcessError as e:
        raise Exception(
            "An error was encountered while running main.nf, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )

    return database
