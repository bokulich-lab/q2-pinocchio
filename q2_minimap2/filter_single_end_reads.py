# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tempfile

import pandas as pd
from q2_types.feature_data import DNAFASTAFormat
from q2_types.per_sample_sequences import SingleLanePerSampleSingleEndFastqDirFmt

from q2_minimap2._filtering_utils import (
    build_filtered_out_dir,
    convert_to_fastq,
    make_mn2_cmd,
    make_samt_cmd,
    process_sam_file,
    run_cmd,
    set_penalties,
)
from q2_minimap2.types._format import Minimap2IndexDBDirFmt


# This function aligns reads to a reference using Minimap2,
# filters the alignments, and converts the filtered alignments
# to FASTQ format, storing the output in a given directory
def _minimap2_filter_reads(
    reads,  # Input reads file path
    outdir,  # Output directory for filtered reads in FASTQ format
    idx_path,  # Path to the Minimap2 index file
    n_threads,  # Number of threads for Minimap2 and samtools
    preset,  # Preset options for Minimap2 alignment
    exclude_mapped,  # Flag to exclude mapped reads, keeping only unmapped
    min_per_identity,  # Minimum percentage identity for a read to be included
    penalties,  # Scoring penalties for Minimap2 alignment
):
    # Create a temporary file for SAM output from Minimap2
    with tempfile.NamedTemporaryFile() as sam_f:
        samf_fp = sam_f.name
        # Create a temporary file for BAM output from samtools view
        with tempfile.NamedTemporaryFile() as bam_f:
            bamf_fp = bam_f.name

            # Construct and execute Minimap2 command for alignment
            mn2_cmd = make_mn2_cmd(
                preset, idx_path, n_threads, penalties, reads, samf_fp
            )
            run_cmd(mn2_cmd, "Minimap2")

            # Process the SAM file to filter alignments based on criteria
            process_sam_file(samf_fp, exclude_mapped, min_per_identity)
            # Construct and execute samtools view command to convert SAM to BAM
            samtools_view_cmd = make_samt_cmd(samf_fp, bamf_fp, n_threads)
            run_cmd(samtools_view_cmd, "samtools view")

            # Construct and execute command to convert BAM to FASTQ
            # using samtools fastq, directing output to the specified output directory
            fwd = str(outdir.path / os.path.basename(reads))
            _reads = ["-0", fwd]
            convert_to_fastq_cmd = convert_to_fastq(_reads, n_threads, bamf_fp)
            run_cmd(convert_to_fastq_cmd, "samtools fastq")


def filter_single_end_reads(
    query_reads: SingleLanePerSampleSingleEndFastqDirFmt,
    index_database: Minimap2IndexDBDirFmt = None,  # Optional pre-built Minimap2 index
    reference_reads: DNAFASTAFormat = None,  # Optional reference sequences
    n_threads: int = 1,  # Number of threads for Minimap2
    mapping_preset: str = "map-ont",  # Minimap2 mapping preset
    keep: str = "mapped",  # Keep 'mapped' or 'unmapped' reads
    min_per_identity: float = None,  # Minimum percentage identity to keep a read
    matching_score: int = None,  # Score for matching bases
    mismatching_penalty: int = None,  # Penalty for mismatched bases
    gap_open_penalty: int = None,  # Penalty for opening a gap
    gap_extension_penalty: int = None,  # Penalty for extending a gap
) -> SingleLanePerSampleSingleEndFastqDirFmt:

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

    # Determine the reference path from the provided index_database or reference_reads
    if index_database:
        idx_ref_path = str(index_database.path) + "/index.mmi"
    elif reference_reads:
        idx_ref_path = str(reference_reads.path)
    else:
        raise ValueError(
            "Either reference_reads or a minimap2_index must be provided as input."
        )

    # Initialize directory format for filtered sequences
    filtered_seqs = SingleLanePerSampleSingleEndFastqDirFmt()

    # Import data from the manifest file to a df
    input_df = query_reads.manifest.view(pd.DataFrame)

    # Set penalties for alignment based on function parameters
    penalties = set_penalties(
        matching_score, mismatching_penalty, gap_open_penalty, gap_extension_penalty
    )

    # Determine whether to exclude mapped reads based on the 'keep' parameter
    exclude_mapped = keep != "mapped"

    # Process each read, filtering according to the specified parameters
    for _, fwd in input_df.itertuples():
        _minimap2_filter_reads(
            fwd,
            filtered_seqs,
            idx_ref_path,
            n_threads,
            mapping_preset,
            exclude_mapped,
            min_per_identity,
            penalties,
        )

    # Build and return the directory of filtered output reads
    result = build_filtered_out_dir(query_reads, filtered_seqs)

    return result
