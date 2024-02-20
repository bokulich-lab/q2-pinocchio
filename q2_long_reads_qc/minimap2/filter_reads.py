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
import subprocess
import tempfile

import pandas as pd
import qiime2.util
import yaml
from q2_types.per_sample_sequences import (
    FastqManifestFormat,
    SingleLanePerSampleSingleEndFastqDirFmt,
    YamlFormat,
)
from q2_types.per_sample_sequences._transformer import (
    _parse_and_validate_manifest_partial,
)

from q2_long_reads_qc._utils import run_command
from q2_long_reads_qc.types._format import Minimap2IndexDBDirFmt


def set_penalties(match, mismatch, gap_o, gap_e):
    options = []
    if match is not None:
        options += ["-A", str(match)]
    if mismatch is not None:
        options += ["-B", str(mismatch)]
    if gap_o is not None:
        options += ["-O", str(gap_o)]
    if gap_e is not None:
        options += ["-E", str(gap_e)]

    return options


# Function to calculate the identity percentage of an alignment.
def calculate_identity(aln, total_length):
    try:
        # Extracts the number of mismatches (NM tag) from the SAM file alignment line.
        nm = int([x for x in aln.split("\t") if x.startswith("NM:i:")][0].split(":")[2])
    except IndexError:
        # Defaults to 0 mismatches if the NM tag is not found.
        nm = 0

    # Calculates matches by subtracting mismatches from total length.
    matches = total_length - nm

    # Calculates identity percentage as the ratio of matches to total alignment length.
    identity_percentage = matches / total_length

    return identity_percentage


# Function to get the alignment length from a CIGAR string.
def get_alignment_length(cigar):
    if cigar == "*":
        # Returns 0 if the CIGAR string is '*', indicating no alignment.
        return 0

    # Extracts all match, insertion, and deletion operations from the CIGAR string.
    matches = re.findall(r"(\d+)([MID])", cigar)

    # Sums the lengths of matches, insertions, and deletions to get total
    # alignment length.
    total_length = sum(int(length) for length, op in matches if op in ["M", "D", "I"])

    return total_length


# Function to process a SAM file, filter based on mappings and identity percentage.
def process_sam_file(input_sam_file, exclude_mapped, min_per_identity):
    # Creates a temporary file to write filtered alignments to.
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as tmp_file:
        temp_file_path = tmp_file.name

    # Opens the input SAM file and the temporary file for writing.
    with open(input_sam_file, "r") as infile, open(temp_file_path, "w") as outfile:
        for line in infile:
            # Writes header lines directly to the output file.
            if line.startswith("@"):
                outfile.write(line)
                continue

            parts = line.split("\t")
            flag = int(parts[1])
            cigar = parts[5]

            # Calculates identity percentage for alignments with a valid CIGAR string.
            if min_per_identity is not None and cigar != "*":
                total_length = get_alignment_length(cigar)
                identity_percentage = calculate_identity(line, total_length)
            else:
                # Defaults identity percentage to 1 (100%) if no CIGAR string or no
                # min_per_identity is specified.
                identity_percentage = 1

            # Logic for excluding or including reads based on mappings and
            # identity percentage.
            if exclude_mapped:
                # Condition for keeping unmapped reads or mapped reads below the
                # identity threshold.
                keep_this_mapped = (
                    min_per_identity is not None
                    and identity_percentage < min_per_identity
                )
                if (flag & 0x4) or keep_this_mapped:
                    outfile.write(line)
            else:
                # Includes reads that are not unmapped and not secondary, based on the
                # identity threshold.
                if not (flag & 0x4) and not (flag & 0x100):
                    if (
                        min_per_identity is not None
                        and identity_percentage >= min_per_identity
                    ) or (min_per_identity is None):
                        outfile.write(line)

    # Replaces the original SAM file with the filtered temporary file.
    shutil.move(temp_file_path, input_sam_file)


def _minimap2_filter(
    reads,
    outdir,
    database,
    n_threads,
    mapping_preset,
    exclude_mapped,
    min_per_identity,
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

            process_sam_file(samfile_output_path, exclude_mapped, min_per_identity)

            samtools_command = [
                "samtools",
                "view",
                "-bS",
                str(samfile_output_path),
                "-o",
                str(bamfile_output_path),
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


def build_filtered_reads_output_directory(filtered_seqs, reads):
    # Parse the input manifest to get a DataFrame of reads
    with reads.manifest.view(FastqManifestFormat).open() as fh:
        input_manifest = _parse_and_validate_manifest_partial(
            fh, single_end=True, absolute=False
        )

    # Filter the input manifest DataFrame for forward reads
    output_df = input_manifest[input_manifest.direction == "forward"]

    # Initialize the output manifest
    output_manifest = FastqManifestFormat()
    # Write the filtered manifest to the output manifest file
    with output_manifest.open() as fh:
        output_df.to_csv(fh, index=False)

    # Initialize the result object to store filtered reads
    result = SingleLanePerSampleSingleEndFastqDirFmt()
    # Write the output manifest to the result object
    result.manifest.write_data(output_manifest, FastqManifestFormat)
    # Duplicate each filtered sequence file to the result object's directory
    for _, _, filename, _ in output_df.itertuples():
        qiime2.util.duplicate(
            str(filtered_seqs.path / filename), str(result.path / filename)
        )

    # Create metadata about the phred offset
    metadata = YamlFormat()
    metadata.path.write_text(yaml.dump({"phred-offset": 33}))
    # Attach metadata to the result
    result.metadata.write_data(metadata, YamlFormat)

    return result


def filter_reads(
    reads: SingleLanePerSampleSingleEndFastqDirFmt,
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

    df = reads.manifest.view(pd.DataFrame)
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
            matching_score,
            mismatching_penalty,
            gap_open_penalty,
            gap_extension_penalty,
        )

    result = build_filtered_reads_output_directory(filtered_seqs, reads)

    return result
