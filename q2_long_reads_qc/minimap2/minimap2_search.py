# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import re
import shutil

import pandas as pd
from q2_types.feature_data import DNAFASTAFormat
from q2_types.per_sample_sequences import SingleLanePerSampleSingleEndFastqDirFmt

from q2_long_reads_qc.filtering._filtering_utils import run_cmd
from q2_long_reads_qc.types._format import (
    Minimap2IndexDBDirFmt,
    PAFDirectoryFormat,
    PAFFormat,
)


# Filter a PAF file to keep up to maxaccepts entries for each read.
def filter_by_maxaccepts(input_paf_path, maxaccepts):
    counts = {}
    filtered_lines = []
    with open(input_paf_path, "r") as infile:
        for line in infile:
            read_id = line.split("\t")[0]
            if counts.get(read_id, 0) < maxaccepts:
                counts[read_id] = counts.get(read_id, 0) + 1
                filtered_lines.append(line)
    # Write the filtered lines back to the file
    with open(input_paf_path, "w") as outfile:
        outfile.writelines(filtered_lines)


# Filter PAF entries by percentage identity.
def filter_by_perc_identity(paf_path, perc_identity):
    with open(paf_path, "r") as file:
        lines = file.readlines()

    filtered_lines = []
    for line in lines:
        parts = line.strip().split("\t")
        identity_score = int(parts[9]) / int(
            parts[10]
        )  # Calculate the BLAST-like alignment identity
        if identity_score >= perc_identity:
            filtered_lines.append(line)

    with open(paf_path, "w") as file:
        file.writelines(filtered_lines)


# Performs sequence alignment using Minimap2 and outputs results in PAF format.
def minimap2_search(
    query_reads: SingleLanePerSampleSingleEndFastqDirFmt,
    minimap2_index: Minimap2IndexDBDirFmt = None,
    reference_reads: DNAFASTAFormat = None,
    maxaccepts: int = 1,
    perc_identity: float = None,
) -> PAFDirectoryFormat:

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

    # Prepare a DataFrame from the query reads manifest for iteration
    input_df = query_reads.manifest.view(pd.DataFrame)

    # Initialize result dictionary to store output paths
    result = PAFDirectoryFormat()
    print("result:", result)

    # Iterate over each sample to perform alignment
    for _, fwd in input_df.itertuples():
        # Extract sample name from the file path
        sample_name_match = re.search(r"/([^/]+)\.fastq", fwd)
        if sample_name_match:
            sample_name = sample_name_match.group(1)

            # Construct output file
            output_path = PAFFormat()

            # Build the Minimap2 command
            cmd = ["minimap2", "-c", idx_ref_path, fwd, "-o", str(output_path)]

            # Execute the Minimap2 alignment command
            run_cmd(cmd, "Minimap2")

            # Filter the PAF file by maxaccepts
            filter_by_maxaccepts(str(output_path), maxaccepts)

            # Optionally filter by perc_identity
            if perc_identity is not None:
                filter_by_perc_identity(str(output_path), perc_identity)

            # Copy the PAF file to the destination
            # Assuming result.path gives the directory path where you want to store
            # the PAF files
            destination_path = os.path.join(result.path, f"{sample_name}.paf")
            shutil.copy(str(output_path), destination_path)

            # Store the result with the sample name as key
            # result[sample_name] = output_path

    return result
