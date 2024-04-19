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
from q2_types.per_sample_sequences import SingleLanePerSamplePairedEndFastqDirFmt

from q2_minimap2._filtering_utils import (
    add_mapped_paired_read_flags,
    build_filtered_paired_end_out_dir,
    collate_bam_inplace,
    convert_to_fastq_paired,
    make_mn2_paired_end_cmd,
    make_samt_cmd,
    process_paired_unmapped_flags,
    process_sam_file,
    run_cmd,
    set_penalties,
)
from q2_minimap2.types._format import Minimap2IndexDBDirFmt


# This function aligns paired-end reads to a reference using Minimap2, filters
# the alignments based on criteria such as mapped/unmapped status and minimum
# percentage identity, and converts the filtered alignments to FASTQ format,
# saving the output in the specified directory.
def _minimap2_filter_paired_end_reads(
    reads1,  # First read of the paired-end sequences
    reads2,  # Second read of the paired-end sequences
    outdir,  # Output directory for the filtered FASTQ files
    idx_path,  # Path to the Minimap2 index database
    n_threads,  # Number of threads for Minimap2 and samtools
    preset,  # Minimap2 alignment preset
    exclude_mapped,  # Flag indicating whether to exclude mapped reads
    min_per_identity,  # Minimum percentage identity for an alignment to be kept
    penalties,  # Alignment penalties
):
    # Temporary files are used for SAM and BAM outputs
    with tempfile.NamedTemporaryFile() as sam_f:
        samf_fp = sam_f.name  # Temporary SAM file path
        with tempfile.NamedTemporaryFile() as bam_f:
            bamf_fp = bam_f.name  # Temporary BAM file path

            # Align reads to the reference using Minimap2 and generate a SAM file
            mn2_cmd = make_mn2_paired_end_cmd(
                preset, idx_path, n_threads, penalties, reads1, reads2, samf_fp
            )
            run_cmd(mn2_cmd, "Minimap2 paired-end")

            # Filter the SAM file using samtools based on the include/exclude criteria
            process_sam_file(samf_fp, exclude_mapped, min_per_identity)
            # Modify flags for paired and unmapped reads
            # making them suitable for samtools fastq command
            process_paired_unmapped_flags(samf_fp)
            add_mapped_paired_read_flags(samf_fp)

            # Convert the filtered SAM file to a BAM file using samtools view
            samtools_view_cmd = make_samt_cmd(samf_fp, bamf_fp, n_threads)
            run_cmd(samtools_view_cmd, "samtools view")

            # Collate the BAM file in place, ensuring proper read grouping
            collate_bam_inplace(bamf_fp)

            # Convert the filtered BAM file to FASTQ format, generating separate files
            # for each read in the pair
            file1 = str(outdir.path / os.path.basename(reads1))
            file2 = str(outdir.path / os.path.basename(reads2))
            _reads = ["-1", file1, "-2", file2]
            convert_to_fastq_cmd = convert_to_fastq_paired(_reads, n_threads, bamf_fp)
            run_cmd(convert_to_fastq_cmd, "samtools fastq")


def filter_paired_end_reads(
    query_reads: SingleLanePerSamplePairedEndFastqDirFmt,  # Input paired-end reads
    index_database: Minimap2IndexDBDirFmt = None,  # Optional pre-built Minimap2 index
    reference_reads: DNAFASTAFormat = None,  # Optional reference sequences for indexing
    n_threads: int = 1,  # Number of threads for Minimap2
    mapping_preset: str = "map-ont",  # Preset options for Minimap2 alignment
    keep: str = "mapped",  # Option to keep 'mapped' or 'unmapped' reads
    min_per_identity: float = None,  # Minimum percentage identity to keep a read
    matching_score: int = None,  # Score for matching bases
    mismatching_penalty: int = None,  # Penalty for mismatched bases
    gap_open_penalty: int = None,  # Penalty for opening a gap
    gap_extension_penalty: int = None,  # Penalty for extending a gap
) -> SingleLanePerSamplePairedEndFastqDirFmt:

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

    # Get the path of index database or reference sequences
    if index_database:
        idx_ref_path = str(index_database.path) + "/index.mmi"
    elif reference_reads:
        idx_ref_path = str(reference_reads.path)

    # Initialize directory format for filtered sequences
    filtered_seqs = SingleLanePerSamplePairedEndFastqDirFmt()

    # Import data from the manifest file to a dataframe
    input_df = query_reads.manifest.view(pd.DataFrame)

    # Configure alignment penalties and scoring
    penalties = set_penalties(
        matching_score, mismatching_penalty, gap_open_penalty, gap_extension_penalty
    )

    # Decide whether to exclude mapped reads based on 'keep' parameter
    exclude_mapped = keep == "unmapped"

    # Process each pair of reads for filtering
    for _, fwd, rev in input_df.itertuples():
        # Apply Minimap2 to filter paired-end reads according to alignment criteria
        _minimap2_filter_paired_end_reads(
            fwd,
            rev,
            filtered_seqs,
            idx_ref_path,
            n_threads,
            mapping_preset,
            exclude_mapped,
            min_per_identity,
            penalties,
        )

    # Compile filtered reads into a new output directory
    result = build_filtered_paired_end_out_dir(query_reads, filtered_seqs)

    return result
