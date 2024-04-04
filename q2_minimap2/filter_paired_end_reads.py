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
    build_filtered_paired_end_out_dir,
    convert_to_fastq_paired,
    make_mn2_paired_end_cmd,
    make_samt_cmd,
    process_sam_file,
    run_cmd,
    set_penalties,
)
from q2_minimap2.types._format import Minimap2IndexDBDirFmt


# This function uses Minimap2 to align reads to a reference, filters the
# resulting SAM file based on mapping criteria using samtools view, and
# finally converts the filtered BAM file to FASTA format using samtools fasta,
# saving the output in the specified directory.
def _minimap2_filter_paired_end_reads(
    reads1,
    reads2,
    outdir,
    idx_path,
    n_threads,
    preset,
    exclude_mapped,
    min_per_identity,
    penalties,
):

    with tempfile.NamedTemporaryFile() as sam_f:
        samf_fp = sam_f.name
        with tempfile.NamedTemporaryFile() as bam_f:
            bamf_fp = bam_f.name

            # Use Minimap2 to find mapped and unmapped reads
            mn2_cmd = make_mn2_paired_end_cmd(
                preset, idx_path, n_threads, penalties, reads1, reads2, samf_fp
            )
            run_cmd(mn2_cmd, "Minimap2 paired-end")

            # Filter sam file using samtools view
            process_sam_file(samf_fp, exclude_mapped, min_per_identity)
            samtools_view_cmd = make_samt_cmd(samf_fp, bamf_fp, n_threads)
            run_cmd(samtools_view_cmd, "samtools view")

            # Convert to FASTQ with samtools
            file1 = str(outdir.path / os.path.basename(reads1))
            file2 = str(outdir.path / os.path.basename(reads2))
            _reads = ["-1", file1, "-2", file2]
            convert_to_fastq_cmd = convert_to_fastq_paired(_reads, n_threads, bamf_fp)
            run_cmd(convert_to_fastq_cmd, "samtools fastq")


def filter_paired_end_reads(
    query_reads: SingleLanePerSamplePairedEndFastqDirFmt,
    index_database: Minimap2IndexDBDirFmt = None,
    reference_reads: DNAFASTAFormat = None,
    n_threads: int = 1,
    mapping_preset: str = "map-ont",
    keep: str = "mapped",
    min_per_identity: float = None,
    matching_score: int = None,
    mismatching_penalty: int = None,
    gap_open_penalty: int = None,
    gap_extension_penalty: int = None,
) -> SingleLanePerSamplePairedEndFastqDirFmt:

    if reference_reads and index_database:
        raise ValueError(
            "Only one reference_reads or index_database artifact "
            "can be provided as input. Choose one and try again."
        )

    # Get the path of index database or reference
    if index_database:
        idx_ref_path = str(index_database.path) + "/index.mmi"
    elif reference_reads:
        idx_ref_path = str(reference_reads.path)
    else:
        raise ValueError(
            "Either reference_reads or a minimap2_index must be provided " "as input."
        )

    # Initialize directory format for filtered sequences
    filtered_seqs = SingleLanePerSamplePairedEndFastqDirFmt()

    # Import data from the manifest file to a df
    input_df = query_reads.manifest.view(pd.DataFrame)

    penalties = set_penalties(
        matching_score, mismatching_penalty, gap_open_penalty, gap_extension_penalty
    )

    if keep == "mapped":
        exclude_mapped = False
    else:
        exclude_mapped = True

    # Iterate over each forward read in the DataFrame
    for _, fwd, rev in input_df.itertuples():
        # Filter the read using minimap2 according to the specified parameters
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

    result = build_filtered_paired_end_out_dir(query_reads, filtered_seqs)

    return result
