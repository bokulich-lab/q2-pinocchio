# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from q2_types.feature_data import FeatureData, Sequence
from q2_types.per_sample_sequences import SequencesWithQuality
from q2_types.sample_data import SampleData
from qiime2.plugin import Bool, Choices, Citations, Float, Int, Plugin, Range, Str

import q2_long_reads_qc
from q2_long_reads_qc import __version__
from q2_long_reads_qc.types._format import Minimap2IndexDBDirFmt, Minimap2IndexDBFmt
from q2_long_reads_qc.types._type import Minimap2IndexDB

citations = Citations.load("citations.bib", package="q2_long_reads_qc")

plugin = Plugin(
    name="long-reads-qc",
    version=__version__,
    website="https://github.com/bokulich-lab/q2-plugin-name",
    package="q2_long_reads_qc",
    description="QIIME 2 Plugin for quality control of long read sequences.",
    short_description="",
)

plugin.register_formats(Minimap2IndexDBDirFmt, Minimap2IndexDBFmt)
plugin.register_semantic_types(Minimap2IndexDB)
plugin.register_semantic_type_to_format(
    Minimap2IndexDB, artifact_format=Minimap2IndexDBDirFmt
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

plugin.methods.register_function(
    function=q2_long_reads_qc.minimap2.filter_reads,
    inputs={
        "reads": SampleData[SequencesWithQuality],
        "minimap2_index": Minimap2IndexDB,
    },
    parameters={
        "n_threads": Int % Range(1, None),
        "mapping_preset": Str % Choices(["map-ont", "map-hifi", "map-pb"]),
        "exclude_mapped": Bool,
        "min_per_identity": Float % Range(0.0, 1.0, inclusive_end=True),
        "matching_score": Int,
        "mismatching_penalty": Int,
        "gap_open_penalty": Int % Range(1, None),
        "gap_extension_penalty": Int % Range(1, None),
    },
    outputs=[("filtered_sequences", SampleData[SequencesWithQuality])],
    input_descriptions={
        "reads": "The sequences to be trimmed.",
        "minimap2_index": "Minimap2 index file.",
    },
    parameter_descriptions={
        "n_threads": "Number of threads to use.",
        "mapping_preset": "Specifies the type of input sequences that will be "
        "used during the mapping process. 1) map-ont: Align noisy long reads "
        "of ~10% error rate to a reference genome. 2) map-hifi: Align PacBio "
        "high-fidelity (HiFi) reads to a reference genome. 3) map-pb: Align "
        "older PacBio continuous long (CLR) reads to a reference genome.",
        "exclude_mapped": "Exclude sequences that align to reference. When "
        "set to False it excludes sequences that do not align to the reference "
        "database.",
        "min_per_identity": "After the alignment step, mapped reads will be "
        "reclassified as unmapped if they identity percentages falls below this "
        "value. If not set there is no reclassification.",
        "matching_score": "Matching score.",
        "mismatching_penalty": "Mismatching penalty.",
        "gap_open_penalty": "Gap open penalty.",
        "gap_extension_penalty": "Gap extension penalty.",
    },
    output_descriptions={
        "filtered_sequences": "The resulting filtered sequences.",
    },
    name="Filter demultiplexed sequences by alignment "
    "to reference database using Minimap2.",
    description=(
        "Filter out (or keep) demultiplexed single- or paired-end "
        "sequences that align to a reference database, using Minimap2 "
        "and samtools. This method can be used to filter out long "
        "sequences or alternatively (when exclude_seqs is False) to "
        "only keep sequences that do align to the reference."
    ),
)


plugin.methods.register_function(
    function=q2_long_reads_qc.minimap2.minimap2_build,
    inputs={"sequences": FeatureData[Sequence]},
    parameters={"kmer_length": Int % Range(1, 28)},
    outputs=[("database", Minimap2IndexDB)],
    input_descriptions={
        "sequences": "Reference sequences used to build Minimap2 index."
    },
    parameter_descriptions={"kmer_length": "Minimizer k-mer length."},
    output_descriptions={"database": "Minimap2 index."},
    name="Build Minimap2 index from reference sequences.",
    description="Build Minimap2 index from reference sequences.",
)
