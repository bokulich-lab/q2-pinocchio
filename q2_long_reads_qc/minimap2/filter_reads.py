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
from q2_long_reads_qc.minimap2._filtering_utils import process_sam_file, set_penalties
from q2_long_reads_qc.types._format import Minimap2IndexDBDirFmt


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


def build_filtered_out_dir(filtered_seqs, reads):
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

    result = build_filtered_out_dir(filtered_seqs, reads)

    return result
