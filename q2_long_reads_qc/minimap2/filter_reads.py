# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import subprocess
import tempfile

# import pysam
from q2_types.per_sample_sequences import CasavaOneEightSingleLanePerSampleDirFmt

from q2_long_reads_qc._utils import run_command
from q2_long_reads_qc.types._format import Minimap2IndexDBDirFmt

# samtools flags
# -f 4 keeps only single alignments that are unmapped
# -f 12 keeps only paired alignments with both reads unmapped
# -F 256 removes reads that are not primary alignment
# -F 260 removes reads that are not primary alignment or unmapped
# -F 268 removes reads that are not primary alignment or unmapped
# or pair is unmapped.
KEEP_UNMAPPED_SINGLE = "4"
KEEP_UNMAPPED_PAIRED = "12"
REMOVE_SECONDARY_ALIGNMENTS = "256"
REMOVE_SECONDARY_OR_UNMAPPED_SINGLE = "260"
REMOVE_SECONDARY_OR_UNMAPPED_PAIRED = "268"


def set_penalties(match, mismatch, gap_o, gap_e):
    options = []
    if match != 2:
        options += ["-A", str(match)]
    if mismatch != 4:
        options += ["-B", str(mismatch)]
    if gap_o != 4:
        options += ["-O", str(gap_o)]
    if gap_e != 2:
        options += ["-E", str(gap_e)]

    return options


def filter_reads(
    reads: CasavaOneEightSingleLanePerSampleDirFmt,
    minimap2_index: Minimap2IndexDBDirFmt,
    n_threads: int = 1,
    mapping_preset: str = "map-ont",
    exclude_mapped: str = False,
    matching_score: int = 2,
    mismatching_penalty: int = 4,
    gap_open_penalty: int = 4,
    gap_extension_penalty: int = 2,
) -> CasavaOneEightSingleLanePerSampleDirFmt:

    filtered_seqs = CasavaOneEightSingleLanePerSampleDirFmt()
    df = reads.manifest
    fastq_paths = [record[1:] for record in df.itertuples()]

    for fwd, rev in fastq_paths:
        _minimap2_filter(
            fwd,
            filtered_seqs,
            minimap2_index,
            n_threads,
            mapping_preset,
            exclude_mapped,
            matching_score,
            mismatching_penalty,
            gap_open_penalty,
            gap_extension_penalty,
        )

    return filtered_seqs


def _minimap2_filter(
    reads,
    outdir,
    database,
    n_threads,
    mapping_preset,
    exclude_mapped,
    match,
    mismatch,
    gap_o,
    gap_e,
):
    alignment_options = set_penalties(match, mismatch, gap_o, gap_e)

    with tempfile.NamedTemporaryFile() as sam_f:
        samfile_output_path = sam_f.name
        with tempfile.NamedTemporaryFile() as bam_f:
            bamfile_output_path = bam_f.name

            # align to reference with Minimap2
            minimap2_cmd = [
                "minimap2",
                "-a",
                "-x",
                mapping_preset,
                str(database.path / "index.mmi"),
                "-t",
                str(n_threads),
            ] + alignment_options

            minimap2_cmd += [reads]
            minimap2_cmd += ["-o", samfile_output_path]

            try:
                # Execute Minimap2
                run_command(minimap2_cmd)
            except subprocess.CalledProcessError as e:
                raise Exception(
                    "An error was encountered while using Minimap2, "
                    f"(return code {e.returncode}), please inspect "
                    "stdout and stderr to learn more."
                )

            # Filter alignment and convert to BAM with samtools
            if exclude_mapped:
                # Keep only the unmapped
                sam_flags = ["-F", "256", "-f", "4"]
            else:
                # Keep only the primary alignments
                sam_flags = ["-F", "260"]

            samtools_command = [
                "samtools",
                "view",
                "-b",
                samfile_output_path,
                "-o",
                bamfile_output_path,
                *sam_flags,
                "-@",
                str(n_threads - 1),
            ]

            try:
                # Execute samtools view
                run_command(samtools_command)
            except subprocess.CalledProcessError as e:
                raise Exception(
                    "An error was encountered while using samtools view, "
                    f"(return code {e.returncode}), please inspect "
                    "stdout and stderr to learn more."
                )

            # Convert to FASTQ with samtools
            fwd = str(outdir.path / os.path.basename(reads))
            _reads = ["-0", fwd]

            # -s /dev/null excludes singletons
            # -n keeps samtools from altering header IDs!
            convert_command = [
                "samtools",
                "fastq",
                *_reads,
                "-s",
                "/dev/null",
                "-@",
                str(n_threads - 1),
                "-n",
                bamfile_output_path,
            ]

            try:
                # Execute samtools fastq
                run_command(convert_command)
            except subprocess.CalledProcessError as e:
                raise Exception(
                    "An error was encountered while using samtools fastq, "
                    f"(return code {e.returncode}), please inspect "
                    "stdout and stderr to learn more."
                )
