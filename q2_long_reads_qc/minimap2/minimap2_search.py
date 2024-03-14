# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import shutil

import pandas as pd
from q2_types.feature_data import DNAFASTAFormat

from q2_long_reads_qc.filtering._filtering_utils import run_cmd
from q2_long_reads_qc.types._format import (
    Minimap2IndexDBDirFmt,
    PairwiseAlignmentMN2DirectoryFormat,
    PairwiseAlignmentMN2Format,
)


# Filter a PairwiseAlignmentMN2 file to keep up to maxaccepts entries for each read.
def filter_by_maxaccepts(input_PairwiseAlignmentMN2_path, maxaccepts):
    counts = {}
    filtered_lines = []
    with open(input_PairwiseAlignmentMN2_path, "r") as infile:
        for line in infile:
            read_id = line.split("\t")[0]
            if counts.get(read_id, 0) < maxaccepts:
                counts[read_id] = counts.get(read_id, 0) + 1
                filtered_lines.append(line)
    # Write the filtered lines back to the file
    with open(input_PairwiseAlignmentMN2_path, "w") as outfile:
        outfile.writelines(filtered_lines)


# Filter PairwiseAlignmentMN2 entries by percentage identity.
def filter_by_perc_identity(PairwiseAlignmentMN2_path, perc_identity):
    with open(PairwiseAlignmentMN2_path, "r") as file:
        lines = file.readlines()

    filtered_lines = []
    for line in lines:
        parts = line.strip().split("\t")
        if int(parts[10]) == 0:
            filtered_lines.append(line)
            continue
        identity_score = int(parts[9]) / int(
            parts[10]
        )  # Calculate the BLAST-like alignment identity
        if identity_score >= perc_identity:
            filtered_lines.append(line)

    with open(PairwiseAlignmentMN2_path, "w") as file:
        file.writelines(filtered_lines)


# Performs sequence alignment using Minimap2 and outputs results in
# PairwiseAlignmentMN2 format.
def minimap2_search(
    query_reads: DNAFASTAFormat,
    minimap2_index: Minimap2IndexDBDirFmt = None,
    reference_reads: DNAFASTAFormat = None,
    n_threads: int = 3,
    maxaccepts: int = 1,
    perc_identity: float = None,
    output_no_hits: bool = True,
) -> pd.DataFrame:

    # Ensure that only one of reference_reads or minimap2_index is provided
    if reference_reads and minimap2_index:
        raise ValueError(
            "Only one of reference_reads or minimap2_index can be provided as input. "
            "Choose one and try again."
        )

    # Ensure that at least one reference type is provided
    if not reference_reads and not minimap2_index:
        raise ValueError(
            "Either reference_reads or minimap2_index must be provided as input."
        )

    # Determine the reference or index path based on input
    idx_ref_path = (
        str(minimap2_index.path / "index.mmi")
        if minimap2_index
        else str(reference_reads.path)
    )

    # Construct output file
    paf_file_fp = PairwiseAlignmentMN2Format()

    # Build the Minimap2 command
    cmd = [
        "minimap2",
        "-c",
        idx_ref_path,
        str(query_reads),
        "-t",
        str(n_threads),
        "-o",
        str(paf_file_fp),
    ]

    if output_no_hits:
        cmd += ["--paf-no-hit"]

    # Execute the Minimap2 alignment command
    run_cmd(cmd, "Minimap2")

    # Filter the PairwiseAlignmentMN2 file by maxaccepts
    filter_by_maxaccepts(str(paf_file_fp), maxaccepts)

    # Optionally filter by perc_identity
    if perc_identity is not None:
        filter_by_perc_identity(str(paf_file_fp), perc_identity)

    # Initialize result dictionary to store the output paf file
    result = PairwiseAlignmentMN2DirectoryFormat()
    destination_path = os.path.join(result.path, "output.paf")
    shutil.copy(str(paf_file_fp), destination_path)
    df = pd.read_csv(destination_path, sep="\t", header=None)

    return df
