# ----------------------------------------------------------------------------
# Copyright (c) 2024, Bokulich Lab.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import subprocess

from q2_types.per_sample_sequences import CasavaOneEightSingleLanePerSampleDirFmt

from q2_pinocchio._utils import run_commands_with_pipe


# Generates a command list for the 'chopper' tool
# incorporating provided parameters for sequence quality trimming
def construct_chopper_command(
    min_quality: int,
    max_quality: int,
    min_length: int,
    max_length: int,
    headcrop: int,
    tailcrop: int,
    threads: int,
) -> list:
    return [
        "chopper",
        "--quality",
        str(min_quality),
        "--maxqual",
        str(max_quality),
        "--minlength",
        str(min_length),
        "--maxlength",
        str(max_length),
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
    query: CasavaOneEightSingleLanePerSampleDirFmt,
    n_threads: int = 4,
    min_quality: int = 0,
    max_quality: int = 1000,
    min_length: int = 1,
    max_length: int = 2147483647,
    headcrop: int = 0,
    tailcrop: int = 0,
) -> CasavaOneEightSingleLanePerSampleDirFmt:

    # Initialize directory format for filtered sequences
    filtered_seqs = CasavaOneEightSingleLanePerSampleDirFmt()

    # Iterate over the FASTQ paired-end files in the DataFrame
    # and exececute chopper command
    for _, fwd, rev in query.manifest.itertuples():
        filtered_seqs_path_fwd = filtered_seqs.path / os.path.basename(fwd)
        chopper_cmd = construct_chopper_command(
            min_quality,
            max_quality,
            min_length,
            max_length,
            headcrop,
            tailcrop,
            n_threads,
        )
        process_and_rezip(fwd, chopper_cmd, str(filtered_seqs_path_fwd))

        if rev:
            filtered_seqs_path_rev = filtered_seqs.path / os.path.basename(rev)
            process_and_rezip(rev, chopper_cmd, str(filtered_seqs_path_rev))

    return filtered_seqs
