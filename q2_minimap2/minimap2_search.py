# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
from q2_types.feature_data import DNAFASTAFormat

from q2_minimap2._filtering_utils import run_cmd
from q2_minimap2.types._format import Minimap2IndexDBDirFmt, PairwiseAlignmentMN2Format


# Filter a PAF file to keep only a certain number of entries
# for each read, as defined by "maxaccepts"
def filter_by_maxaccepts(input_PairwiseAlignmentMN2_path, maxaccepts):
    # Dictionary to keep track of the counts of occurancies/hits for each read ID
    counts = {}
    # List to hold the lines that meet the filtering criteria
    filtered_lines = []

    # Open the input file in read mode
    with open(input_PairwiseAlignmentMN2_path, "r") as infile:
        # Iterate over each line in the input file
        for line in infile:
            # Extract the read ID from the current line
            read_id = line.split("\t")[0]

            # Check and update the count for the specific read ID
            # We need this in order to not exceed the maximum
            # number of hits allowed per read ID
            if counts.get(read_id, 0) < maxaccepts:
                counts[read_id] = counts.get(read_id, 0) + 1
                filtered_lines.append(line)

    # Write the filtered lines back to the file
    with open(input_PairwiseAlignmentMN2_path, "w") as outfile:
        outfile.writelines(filtered_lines)


# Filter PAF entries based on a threshold of percentage identity
def filter_by_perc_identity(PairwiseAlignmentMN2_path, perc_identity, output_no_hits):
    # Open and read all lines from the input file
    with open(PairwiseAlignmentMN2_path, "r") as file:
        lines = file.readlines()

    # Initialize list to hold lines that pass the filter
    filtered_lines = []
    for line in lines:
        # Split each line of the PAF file into its components (columns)
        parts = line.strip().split("\t")

        # Add unmapped entry if it exists
        if (int(parts[10]) == 0) and (output_no_hits is True):
            filtered_lines.append(line)
            continue

        # Calculate the BLAST-like alignment identity
        # https://lh3.github.io/minimap2/minimap2.html - OUTPUT FORMAT
        identity_score = int(parts[9]) / int(parts[10])

        # If the identity score meets or exceeds the threshold
        # include the line in the filtered output
        if identity_score >= perc_identity:
            filtered_lines.append(line)
        else:
            # Modify the line to mark it as unmapped while keeping the
            # first two columns intact
            if output_no_hits is True:
                modified_line = (
                    f"{parts[0]}\t{parts[1]}\t0\t0\t*\t*\t0\t0\t0\t0\t0\t0\n"
                )
                filtered_lines.append(modified_line)

    # Overwrite the input file with only the lines that met the filtering criteria
    with open(PairwiseAlignmentMN2_path, "w") as file:
        file.writelines(filtered_lines)


# Construct the command list for the Minimap2 alignment search
def construct_command(
    idx_ref_path, query_reads, n_threads, paf_file_fp, output_no_hits
):
    cmd = [
        "minimap2",
        "-c",
        str(idx_ref_path),
        str(query_reads),
        "-t",
        str(n_threads),
        "-o",
        str(paf_file_fp),
    ]
    if output_no_hits:
        cmd.append("--paf-no-hit")
    return cmd


# Performs sequence alignment using Minimap2 and outputs results in
# PairwiseAlignmentMN2 format.
def minimap2_search(
    query_reads: DNAFASTAFormat,
    index_database: Minimap2IndexDBDirFmt = None,
    reference_reads: DNAFASTAFormat = None,
    n_threads: int = 3,
    maxaccepts: int = 1,
    perc_identity: float = None,
    output_no_hits: bool = True,
) -> pd.DataFrame:

    # Ensure that only one of reference_reads or index_database is provided
    if reference_reads and index_database:
        raise ValueError(
            "Only one of reference_reads or index_database can be provided as input. "
            "Choose one and try again."
        )

    # Ensure that at least one of reference_reads and index_database is provided
    if not reference_reads and not index_database:
        raise ValueError(
            "Either reference_reads or index_database must be provided as input."
        )

    # Determine the reference or index path based on input
    idx_ref_path = (
        str(index_database.path / "index.mmi")
        if index_database
        else str(reference_reads.path)
    )

    # Create a reference to a file with PAF format
    paf_file_fp = PairwiseAlignmentMN2Format()

    # Construct the
    cmd = construct_command(
        idx_ref_path, query_reads, n_threads, paf_file_fp, output_no_hits
    )

    # Execute the Minimap2 alignment command
    run_cmd(cmd, "Minimap2")

    # Filter the PAF file by maxaccepts (default = 1)
    filter_by_maxaccepts(str(paf_file_fp), maxaccepts)

    # Optionally filter by perc_identity
    if perc_identity is not None:
        filter_by_perc_identity(str(paf_file_fp), perc_identity, output_no_hits)

    # Read the PAF file as a pandas DataFrame
    df = pd.read_csv(str(paf_file_fp), sep="\t", header=None)

    return df
