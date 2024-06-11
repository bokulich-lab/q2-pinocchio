# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import subprocess

from q2_types.feature_data import DNAFASTAFormat

from q2_minimap2._utils import run_command
from q2_minimap2.types._format import Minimap2IndexDBDirFmt


def create_idx_build_cmd(database, sequences, preset, kmer_length):
    build_cmd = [
        "minimap2",
        "-x",
        preset,
        "-d",  # Option for generating index
        str(database.path / "index.mmi"),  # Output index path
        str(sequences),  # Input sequences file path
        "-k",  # Option to specify k-mer length
        str(kmer_length),  # k-mer length
    ]

    return build_cmd


def build_index(
    sequences: DNAFASTAFormat,
    preset: str = "map-ont",  # Preset to apply multiple options
    kmer_length: int = 15,  # Minimap2 uses this value as the default value
) -> Minimap2IndexDBDirFmt:
    # Initialize a directory format object to store the Minimap2 index
    database = Minimap2IndexDBDirFmt()

    # Construct the command to build the Minimap2 index file
    build_cmd = create_idx_build_cmd(database, sequences, preset, kmer_length)

    try:
        # Execute the command to create the Minimap2 index database
        run_command(build_cmd)
    except subprocess.CalledProcessError as e:
        raise Exception(
            "An error was encountered while running main.nf, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )

    return database
