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
    # Read the input file into a DataFrame
    df = pd.read_csv(input_PairwiseAlignmentMN2_path, sep="\t", header=None)

    # Group by read_id and count occurrences
    counts = df.groupby(0).cumcount() + 1

    # Filter the DataFrame based on maxaccepts
    filtered_df = df[counts <= maxaccepts]

    # Write the output DataFrame back to a TSV file
    filtered_df.to_csv(
        input_PairwiseAlignmentMN2_path, sep="\t", header=False, index=False
    )


# Filter PAF entries based on a threshold of percentage identity
def filter_by_perc_identity(PairwiseAlignmentMN2_path, perc_identity, output_no_hits):
    # Read the input file into a DataFrame
    df = pd.read_csv(PairwiseAlignmentMN2_path, sep="\t", header=None)

    # Add a new column containing the row numbers
    df["row_num"] = df.index

    # Calculate the BLAST-like alignment identity only if there are mapped bases
    df["identity_score"] = df.apply(
        lambda row: row[9] / row[10] if row[10] > 0 else 0, axis=1
    )

    # Filter based on identity score
    mapped_df = df[df["identity_score"] >= perc_identity]
    filtered_out = df[df["identity_score"] < perc_identity]

    # Keep only the first entry/row for each unique value in column 1 (query)
    # in filtered_out
    filtered_out = filtered_out.drop_duplicates(subset=1)

    # Set columns 2 to 22 to '*'
    filtered_out.iloc[:, 2:23] = "*"

    # Set all columns after the second to 0
    filtered_out.iloc[:, 2:23] = 0

    # Set columns 4 and 5 to '*'
    filtered_out.iloc[:, 4:6] = "*"

    # Merging the rows of the two DataFrames based on row number
    merged_df = pd.concat([mapped_df, filtered_out], axis=0)
    merged_df = merged_df.sort_values(by="row_num")
    merged_df.reset_index(drop=True, inplace=True)

    # Drop unnecessary columns
    merged_df.drop(columns=["row_num", "identity_score"], inplace=True)

    # Write the filtered DataFrame back to the input file
    merged_df.to_csv(PairwiseAlignmentMN2_path, sep="\t", header=False, index=False)


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

    # Construct the command
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
