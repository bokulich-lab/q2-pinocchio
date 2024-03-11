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
from qiime2.plugin import Citations, Float, Int, Plugin, Range

import q2_long_reads_qc
from q2_long_reads_qc import __version__
from q2_long_reads_qc._action_params import (
    filter_reads_param_descriptions,
    filter_reads_params,
    minimap2_build_param_descriptions,
    minimap2_build_params,
)
from q2_long_reads_qc.types._format import (
    Minimap2IndexDBDirFmt,
    Minimap2IndexDBFmt,
    PAFDirectoryFormat,
    PAFFormat,
)
from q2_long_reads_qc.types._type import PAF, Minimap2IndexDB

citations = Citations.load("citations.bib", package="q2_long_reads_qc")

plugin = Plugin(
    name="long-reads-qc",
    version=__version__,
    website="https://github.com/bokulich-lab/q2-plugin-name",
    package="q2_long_reads_qc",
    description="QIIME 2 Plugin for quality control of long read sequences.",
    short_description="",
)

plugin.register_formats(
    Minimap2IndexDBDirFmt, Minimap2IndexDBFmt, PAFFormat, PAFDirectoryFormat
)
plugin.register_semantic_types(Minimap2IndexDB, PAF)
plugin.register_semantic_type_to_format(
    Minimap2IndexDB, artifact_format=Minimap2IndexDBDirFmt
)
# plugin.register_semantic_type_to_format(PAF, artifact_format=PAFDirectoryFormat)

plugin.register_semantic_type_to_format(SampleData[PAF], PAFDirectoryFormat)

plugin.methods.register_function(
    function=q2_long_reads_qc.filtering.filter_reads,
    inputs={
        "query_reads": SampleData[SequencesWithQuality],
        "minimap2_index": Minimap2IndexDB,
        "reference_reads": FeatureData[Sequence],
    },
    parameters=filter_reads_params,
    outputs=[("filtered_query_reads", SampleData[SequencesWithQuality])],
    input_descriptions={
        "query_reads": "The sequences to be filtered.",
        "minimap2_index": "Minimap2 index file. Incompatible with reference_reads.",
        "reference_reads": "Reference sequences. Incompatible with minimap2_index.",
    },
    parameter_descriptions=filter_reads_param_descriptions,
    output_descriptions={
        "filtered_query_reads": "The resulting filtered sequences.",
    },
    name="Filter long sequences using Minimap2.",
    description=(
        "Filter demultiplexed single- or paired-end sequences based "
        "on alignment to a reference database using Minimap2 and samtools. "
        "This versatile command allows for the exclusion of long sequences "
        "or, when exclude_mapped is set to False, selectively retains only "
        "those sequences aligning to the reference."
    ),
)


plugin.methods.register_function(
    function=q2_long_reads_qc.minimap2.minimap2_build,
    inputs={"sequences": FeatureData[Sequence]},
    parameters=minimap2_build_params,
    outputs=[("database", Minimap2IndexDB)],
    input_descriptions={
        "sequences": "Reference sequences used to build Minimap2 index."
    },
    parameter_descriptions=minimap2_build_param_descriptions,
    output_descriptions={"database": "Minimap2 index."},
    name="Build Minimap2 index from reference sequences.",
    description="Build Minimap2 index from reference sequences.",
)


plugin.methods.register_function(
    function=q2_long_reads_qc.minimap2.minimap2_search,
    inputs={
        "query_reads": SampleData[SequencesWithQuality],
        "minimap2_index": Minimap2IndexDB,
        "reference_reads": FeatureData[Sequence],
    },
    parameters={
        "maxaccepts": Int % Range(1, None),
        "perc_identity": Float % Range(0.0, 1.0, inclusive_end=True),
    },
    outputs=[("search_results", SampleData[PAF])],
    input_descriptions={
        "query_reads": "Query sequences.",
        "minimap2_index": "Minimap2 index file. Incompatible with reference_reads.",
        "reference_reads": "Reference sequences. Incompatible with minimap2_index.",
    },
    parameter_descriptions={
        "maxaccepts": "Maximum number of hits to keep for each query. Minimap2 will "
        "choose the first N hits in the reference database. NOTE: the database "
        "is not sorted by similarity to query, so these are the "
        "first N hits that pass the threshold, not necessarily the top N hits.",
        "perc_identity": "Optionally reject match if percent identity to query is "
        "lower.",
    },
    output_descriptions={
        "search_results": "Top hits for each query.",
    },
    name="Minimap2 alignment search.",
    description=(
        "Search for top hits in a reference database using alignment between the "
        "query sequences and reference database sequences using Minimap2. Returns a "
        "report of the top M hits for each query (where M=maxaccepts)."
    ),
)
