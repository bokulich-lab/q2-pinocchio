# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import shutil
import subprocess

from q2_types.per_sample_sequences import CasavaOneEightSingleLanePerSampleDirFmt

from q2_minimap2._utils import run_commands_with_pipe


# Generates a command list for the 'chopper' tool
# incorporating provided parameters for sequence quality trimming
def construct_chopper_command(
    quality: int,
    maxqual: int,
    minlength: int,
    maxlength: int,
    headcrop: int,
    tailcrop: int,
    threads: int,
) -> list:
    return [
        "chopper",
        "--quality",
        str(quality),
        "--maxqual",
        str(maxqual),
        "--minlength",
        str(minlength),
        "--maxlength",
        str(maxlength),
        "--headcrop",
        str(headcrop),
        "--tailcrop",
        str(tailcrop),
        "--threads",
        str(threads),
    ]


# Executes a pipeline that unzips FASTQ files, processes them with 'chopper'
# and then rezips them
def process_and_rezip(input_file, chopper_cmd, filtered_seqs_path):
    try:
        run_commands_with_pipe(
            ["gunzip", "-c", str(input_file)],
            chopper_cmd,
            ["gzip"],
            str(filtered_seqs_path),
        )
    except subprocess.CalledProcessError as e:
        raise Exception(
            f"An error was encountered while using chopper, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )


# Trims paired-end read FASTQ files using specified quality control parameter
def trim(
    query_reads: CasavaOneEightSingleLanePerSampleDirFmt,
    threads: int = 4,
    quality: int = 0,
    maxqual: int = 1000,
    minlength: int = 1,
    maxlength: int = 2147483647,
    headcrop: int = 0,
    tailcrop: int = 0,
) -> CasavaOneEightSingleLanePerSampleDirFmt:

    # Initialize directory format for filtered sequences
    filtered_seqs = CasavaOneEightSingleLanePerSampleDirFmt()
    shutil.copy(
        os.path.join(query_reads.path, "MANIFEST"),
        os.path.join(filtered_seqs.path, "MANIFEST"),
    )
    shutil.copy(
        os.path.join(query_reads.path, "metadata.yml"),
        os.path.join(filtered_seqs.path, "metadata.yml"),
    )

    # Iterate over the FASTQ paired-end files in the DataFrame
    # and exececute chopper command
    for _, fwd, rev in query_reads.manifest.itertuples():
        filtered_seqs_path_fwd = filtered_seqs.path / os.path.basename(fwd)
        if rev:
            filtered_seqs_path_rev = filtered_seqs.path / os.path.basename(rev)

        chopper_cmd = construct_chopper_command(
            quality, maxqual, minlength, maxlength, headcrop, tailcrop, threads
        )

        process_and_rezip(fwd, chopper_cmd, str(filtered_seqs_path_fwd))
        if rev:
            process_and_rezip(rev, chopper_cmd, str(filtered_seqs_path_rev))

    return filtered_seqs
