# ----------------------------------------------------------------------------
# Copyright (c) 2022, <developer name>.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

# from q2_types.feature_data import FeatureData, Sequence
from q2_types.per_sample_sequences import (
    PairedEndSequencesWithQuality,
    SequencesWithQuality,
)
from q2_types.sample_data import SampleData
from qiime2.plugin import Citations, Plugin, TypeMatch

import q2_long_reads_qc
from q2_long_reads_qc import __version__
from q2_long_reads_qc.types._format import (
    Minimap2IndexFileDirFmt,
    Minimap2IndexFileFormat,
)
from q2_long_reads_qc.types._type import Minimap2Index

citations = Citations.load("citations.bib", package="q2_long_reads_qc")

plugin = Plugin(
    name="long-reads-qc",
    version=__version__,
    website="https://github.com/bokulich-lab/q2-plugin-name",
    package="q2_long_reads_qc",
    description="QIIME 2 Plugin for quality control of long read sequences.",
    short_description="",
)

plugin.register_formats(Minimap2IndexFileDirFmt, Minimap2IndexFileFormat)
plugin.register_semantic_types(Minimap2Index)
plugin.register_semantic_type_to_format(
    Minimap2Index, artifact_format=Minimap2IndexFileDirFmt
)

plugin.visualizers.register_function(
    function=q2_long_reads_qc.longqc.evaluate.evaluate_long_reads,
    inputs={
        "sequences": SampleData[SequencesWithQuality],
    },
    parameters="",
    input_descriptions={
        "sequences": "Long reads to be analyzed.",
    },
    parameter_descriptions="",
    name="Evaluate long sequences using LongQC.",
    description="This method uses LongQC "
    " to assess the quality of long reads providing"
    " visualizations summarizing the results.",
    citations=[citations["Fukasawa_LongQC"]],
)


T = TypeMatch([SequencesWithQuality, PairedEndSequencesWithQuality])
plugin.methods.register_function(
    function=q2_long_reads_qc.minimap2.filter_reads,
    inputs={"demultiplexed_sequences": SampleData[T]},
    parameters="",
    outputs=[("filtered_sequences", SampleData[T])],
    input_descriptions={
        "demultiplexed_sequences": "The sequences to be trimmed.",
    },
    parameter_descriptions="",
    output_descriptions={"filtered_sequences": "The resulting filtered sequences."},
    name="Filter demultiplexed sequences by alignment to reference database.",
    description=(
        "Filter out (or keep) demultiplexed single- or paired-end "
        "sequences that align to a reference database, using minimap2 "
        "and samtools. This method can be used to filter out long "
        "sequences or alternatively (when exclude_seqs is False) to "
        "only keep sequences that do align to the reference."
    ),
)

"""
plugin.methods.register_function(
    function=q2_long_reads_qc.minimap2.minimap2_build,
    inputs={'sequences': FeatureData[Sequence]},
    parameters={'n_threads': Int % Range(1, None)},
    outputs=[('database', Bowtie2Index)],
    input_descriptions={
        'sequences': 'Reference sequences used to build Minimap2 index.'},
    parameter_descriptions={'n_threads': 'Number of threads to launch'},
    output_descriptions={'database': 'Minimap2 index.'},
    name='Build Minimap2 index from reference sequences.',
    description='Build Minimap2 index from reference sequences.',
)

"""
