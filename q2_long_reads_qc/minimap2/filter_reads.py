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
from q2_types.per_sample_sequences import SingleLanePerSampleSingleEndFastqDirFmt

from q2_long_reads_qc.minimap2._filtering_utils import (
    build_filtered_out_dir,
    make_convert_cmd,
    make_mn2_cmd,
    make_samt_cmd,
    process_sam_file,
    run_cmd,
    set_penalties,
)
from q2_long_reads_qc.types._format import Minimap2IndexDBDirFmt


# This function uses Minimap2 to align reads to a reference, filters the
# resulting SAM file based on mapping criteria using samtools view, and
# finally converts the filtered BAM file to FASTQ format using samtools fastq,
# saving the output in the specified directory.
def _minimap2_filter(
    reads,
    outdir,
    idx,
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
            mn2_cmd = make_mn2_cmd(preset, idx, n_threads, penalties, reads, samf_fp)
            run_cmd(mn2_cmd, "Minimap2")

            # Filter sam file using samtools view
            process_sam_file(samf_fp, exclude_mapped, min_per_identity)
            samtools_view_cmd = make_samt_cmd(samf_fp, bamf_fp, n_threads)
            run_cmd(samtools_view_cmd, "samtools view")

            # Convert to FASTQ with samtools
            fwd = str(outdir.path / os.path.basename(reads))
            _reads = ["-0", fwd]
            convert_to_fastq_cmd = make_convert_cmd(_reads, n_threads, bamf_fp)
            run_cmd(convert_to_fastq_cmd, "samtools fastq")


def filter_reads(
    sequences: SingleLanePerSampleSingleEndFastqDirFmt,
    minimap2_index: Minimap2IndexDBDirFmt,
    n_threads: int = 1,
    mapping_preset: str = "map-ont",
    exclude_mapped: str = False,
    min_per_identity: float = None,
    matching_score: int = None,
    mismatching_penalty: int = None,
    gap_open_penalty: int = None,
    gap_extension_penalty: int = None,
) -> SingleLanePerSampleSingleEndFastqDirFmt:
    # Initialize directory format for filtered sequences
    filtered_seqs = SingleLanePerSampleSingleEndFastqDirFmt()
    # Import data from the manifest file to a df
    df = sequences.manifest.view(pd.DataFrame)

    penalties = set_penalties(
        matching_score, mismatching_penalty, gap_open_penalty, gap_extension_penalty
    )

    # Iterate over each forward read in the DataFrame
    for _, fwd in df.itertuples():
        # Filter the read using minimap2 according to the specified parameters
        _minimap2_filter(
            fwd,
            filtered_seqs,
            minimap2_index,
            n_threads,
            mapping_preset,
            exclude_mapped,
            min_per_identity,
            penalties,
        )

    result = build_filtered_out_dir(filtered_seqs, sequences)

    return result
