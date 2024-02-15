# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tempfile

from q2_types.per_sample_sequences import CasavaOneEightSingleLanePerSampleDirFmt

from q2_long_reads_qc._utils import run_command
from q2_long_reads_qc.types._format import Minimap2IndexFileDirFmt

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


def filter_reads(
    reads: CasavaOneEightSingleLanePerSampleDirFmt,
    minimap2_index: Minimap2IndexFileDirFmt,
    n_threads: int = 1,
    exclude_mapped: str = False,
) -> CasavaOneEightSingleLanePerSampleDirFmt:

    filtered_seqs = CasavaOneEightSingleLanePerSampleDirFmt()
    df = reads.manifest
    fastq_paths = [record[1:] for record in df.itertuples()]

    for fwd, rev in fastq_paths:
        _minimap2_filter(fwd, filtered_seqs, minimap2_index, n_threads, exclude_mapped)

    return filtered_seqs


def _minimap2_filter(f_read, outdir, database, n_threads, exclude_mapped):
    with tempfile.NamedTemporaryFile() as sam_f:
        samfile_output_path = sam_f.name
        with tempfile.NamedTemporaryFile() as bam_f:
            bamfile_output_path = bam_f.name

            # align to reference with Minimap2
            minimap2_cmd = [
                "minimap2",
                "-a",
                "-x",
                "map-ont",
                str(database.path / "index.mmi"),
                "-t",
                str(n_threads),
            ]

            minimap2_cmd += [f_read]
            minimap2_cmd += ["-o", samfile_output_path]
            run_command(minimap2_cmd)

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
            run_command(samtools_command)

            # Convert to FASTQ with samtools
            fwd = str(outdir.path / os.path.basename(f_read))
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
            run_command(convert_command)
