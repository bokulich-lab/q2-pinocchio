# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import re
import shutil
import subprocess
import tempfile

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

from q2_minimap2._utils import run_command


# Set Minimap2 alignment penalties based on provided parameters
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


# Function to calculate the identity percentage of an alignment
def calculate_identity(aln, total_length):
    try:
        # Extracts the number of mismatches (NM tag) from the SAM file alignment line
        nm = int([x for x in aln.split("\t") if x.startswith("NM:i:")][0].split(":")[2])
    except IndexError:
        # Defaults to 0 mismatches if the NM tag is not found
        nm = 0

    # Calculates matches by subtracting mismatches from total length
    matches = total_length - nm

    # Calculates identity percentage as the ratio of matches to total alignment length.
    identity_percentage = matches / total_length

    return identity_percentage


# Function to get the alignment length from a CIGAR string
def get_alignment_length(cigar):
    if cigar == "*":
        # Returns 0 if the CIGAR string is '*', indicating no alignment
        return 0

    # Extracts all match, insertion, and deletion operations from the CIGAR string
    matches = re.findall(r"(\d+)([MID])", cigar)

    # Sums the lengths of matches, insertions, and deletions to get total
    # alignment length
    total_length = sum(int(length) for length, op in matches if op in ["M", "D", "I"])

    return total_length


# Function to process a SAM file, filter based on mappings and identity percentage
def process_sam_file(input_sam_file, exclude_mapped, min_per_identity):
    # Creates a temporary file to write filtered alignments to
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as tmp_file:
        temp_file_path = tmp_file.name

    # Opens the input SAM file and the temporary file for writing
    with open(input_sam_file, "r") as infile, open(temp_file_path, "w") as outfile:
        for line in infile:
            # Writes header lines directly to the output file
            if line.startswith("@"):
                outfile.write(line)
                continue

            parts = line.split("\t")
            flag = int(parts[1])
            cigar = parts[5]

            # Calculates identity percentage for alignments with a valid CIGAR string
            if min_per_identity is not None and cigar != "*":
                total_length = get_alignment_length(cigar)
                identity_percentage = calculate_identity(line, total_length)
            else:
                # Defaults identity percentage to 1 (100%) if no CIGAR string or no
                # min_per_identity is specified
                identity_percentage = 1

            # Logic for excluding or including reads based on mappings and
            # identity percentage
            if exclude_mapped:
                # Condition for keeping unmapped reads or mapped reads below the
                # identity threshold
                keep_this_mapped = (
                    min_per_identity is not None
                    and identity_percentage < min_per_identity
                )
                if (flag & 0x4) or keep_this_mapped:
                    outfile.write(line)
            else:
                # Includes reads that are not unmapped and not secondary, based on the
                # identity threshold
                if not (flag & 0x4) and not (flag & 0x100):
                    if (
                        min_per_identity is not None
                        and identity_percentage >= min_per_identity
                    ) or (min_per_identity is None):
                        outfile.write(line)

    # Replaces the original SAM file with the filtered temporary file
    shutil.move(temp_file_path, input_sam_file)


# Generate samtools fasta convert command
def convert_to_fasta(_reads, n_threads, bamfile_filepath):
    # -s /dev/null excludes singletons
    # -n keeps samtools from altering header IDs
    convert_cmd = [
        "samtools",
        "fasta",
        "-0",
        str(_reads),
        "-s",
        "/dev/null",
        "-@",
        str(n_threads - 1),
        "-n",
        str(bamfile_filepath),
    ]

    return convert_cmd


def convert_to_fastq(_reads, n_threads, bamfile_filepath):
    convert_cmd = [
        "samtools",
        "fastq",
        *_reads,
        "-s",
        "/dev/null",
        "-@",
        str(n_threads - 1),
        "-n",
        str(bamfile_filepath),
    ]

    return convert_cmd


# Generate samtools view filtering command
def make_samt_cmd(samfile_filepath, bamfile_filepath, n_threads):
    samtools_cmd = [
        "samtools",
        "view",
        "-bS",
        str(samfile_filepath),
        "-o",
        str(bamfile_filepath),
        "-@",
        str(n_threads - 1),
    ]

    return samtools_cmd


# Generate Minimap2 mapping command
def make_mn2_cmd(mapping_preset, index, n_threads, penalties, reads, samf_fp):
    # align to reference with Minimap2
    minimap2_cmd = [
        "minimap2",
        "-a",
        "-x",
        mapping_preset,
        str(index),
        "-t",
        str(n_threads),
    ] + penalties

    minimap2_cmd += [reads]
    minimap2_cmd += ["-o", samf_fp]

    return minimap2_cmd


# Helper function for command execution
def run_cmd(cmd, str):
    try:
        # Execute samtools fastq
        run_command(cmd)
    except subprocess.CalledProcessError as e:
        raise Exception(
            f"An error was encountered while using {str}, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )


def build_filtered_out_dir(input_reads, filtered_seqs):
    # Parse the input manifest to get a DataFrame of reads
    with input_reads.manifest.view(FastqManifestFormat).open() as fh:
        input_manifest = _parse_and_validate_manifest_partial(
            fh, single_end=True, absolute=False
        )
        # Filter the input manifest DataFrame for forward reads
        output_df = input_manifest[input_manifest.direction == "forward"]

    # Initialize the output manifest
    output_manifest = FastqManifestFormat()
    # Copy input manifest to output manifest
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
