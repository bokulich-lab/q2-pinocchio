# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import importlib

from q2_types.feature_data import FeatureData
from qiime2.plugin import Citations, Plugin

import q2_minimap2
from q2_minimap2 import __version__
from q2_minimap2._action_params import (
    build_index_dsc,
    build_index_inputs,
    build_index_inputs_dsc,
    build_index_outputs,
    build_index_outputs_dsc,
    build_index_param_dsc,
    build_index_params,
    classify_consensus_minimap2_dsc,
    classify_consensus_minimap2_inputs,
    classify_consensus_minimap2_inputs_dsc,
    classify_consensus_minimap2_outputs,
    classify_consensus_minimap2_outputs_dsc,
    classify_consensus_minimap2_param_dsc,
    classify_consensus_minimap2_params,
    extract_seqs_dsc,
    extract_seqs_inputs,
    extract_seqs_inputs_dsc,
    extract_seqs_outputs,
    extract_seqs_outputs_dsc,
    extract_seqs_param_dsc,
    extract_seqs_params,
    filter_reads_dsc,
    filter_reads_inputs,
    filter_reads_inputs_dsc,
    filter_reads_outputs,
    filter_reads_outputs_dsc,
    filter_reads_param_dsc,
    filter_reads_params,
    find_consensus_annotation_dsc,
    find_consensus_annotation_inputs,
    find_consensus_annotation_inputs_dsc,
    find_consensus_annotation_outputs,
    find_consensus_annotation_outputs_dsc,
    find_consensus_annotation_params,
    find_consensus_annotation_params_dsc,
    minimap2_dsc,
    minimap2_inputs,
    minimap2_inputs_dsc,
    minimap2_outputs,
    minimap2_outputs_dsc,
    minimap2_param_dsc,
    minimap2_params,
)
from q2_minimap2.types._format import (
    Minimap2IndexDBDirFmt,
    Minimap2IndexDBFmt,
    PairwiseAlignmentMN2DirectoryFormat,
    PairwiseAlignmentMN2Format,
)
from q2_minimap2.types._type import Minimap2IndexDB, PairwiseAlignmentMN2

citations = Citations.load("citations.bib", package="q2_minimap2")

plugin = Plugin(
    name="minimap2",
    version=__version__,
    website="https://github.com/bokulich-lab/q2-plugin-name",
    package="q2_minimap2",
    description="QIIME 2 Plugin for quality control and taxonomic "
    "classification of long read sequences using Minimap2",
    short_description="",
)

plugin.register_formats(
    Minimap2IndexDBDirFmt,
    Minimap2IndexDBFmt,
    PairwiseAlignmentMN2Format,
    PairwiseAlignmentMN2DirectoryFormat,
)
plugin.register_semantic_types(Minimap2IndexDB, PairwiseAlignmentMN2)
plugin.register_semantic_type_to_format(
    Minimap2IndexDB, artifact_format=Minimap2IndexDBDirFmt
)

plugin.register_semantic_type_to_format(
    FeatureData[PairwiseAlignmentMN2], PairwiseAlignmentMN2DirectoryFormat
)

plugin.methods.register_function(
    function=q2_minimap2.extract_seqs,
    inputs=extract_seqs_inputs,
    parameters=extract_seqs_params,
    outputs=extract_seqs_outputs,
    input_descriptions=extract_seqs_inputs_dsc,
    parameter_descriptions=extract_seqs_param_dsc,
    output_descriptions=extract_seqs_outputs_dsc,
    name="Extract feature sequences using Minimap2.",
    description=extract_seqs_dsc,
    citations=[citations["Minimap2"]],
)

plugin.methods.register_function(
    function=q2_minimap2.filter_reads,
    inputs=filter_reads_inputs,
    parameters=filter_reads_params,
    outputs=filter_reads_outputs,
    input_descriptions=filter_reads_inputs_dsc,
    parameter_descriptions=filter_reads_param_dsc,
    output_descriptions=filter_reads_outputs_dsc,
    name="Filter demultiplexed single- or paired-end sequences long sequences "
    "using Minimap2.",
    description=filter_reads_dsc,
    citations=[citations["Minimap2"]],
)

plugin.methods.register_function(
    function=q2_minimap2.build_index,
    inputs=build_index_inputs,
    parameters=build_index_params,
    outputs=build_index_outputs,
    input_descriptions=build_index_inputs_dsc,
    parameter_descriptions=build_index_param_dsc,
    output_descriptions=build_index_outputs_dsc,
    name="Build Minimap2 index database from reference sequences.",
    description=build_index_dsc,
    citations=[citations["Minimap2"]],
)

plugin.methods.register_function(
    function=q2_minimap2.minimap2,
    inputs=minimap2_inputs,
    parameters=minimap2_params,
    outputs=minimap2_outputs,
    input_descriptions=minimap2_inputs_dsc,
    parameter_descriptions=minimap2_param_dsc,
    output_descriptions=minimap2_outputs_dsc,
    name="Minimap2 alignment search.",
    description=minimap2_dsc,
    citations=[citations["Minimap2"], citations["li2009sequence"]],
)

plugin.pipelines.register_function(
    function=q2_minimap2.classify_consensus_minimap2,
    inputs=classify_consensus_minimap2_inputs,
    parameters=classify_consensus_minimap2_params,
    outputs=classify_consensus_minimap2_outputs,
    input_descriptions=classify_consensus_minimap2_inputs_dsc,
    parameter_descriptions=classify_consensus_minimap2_param_dsc,
    output_descriptions=classify_consensus_minimap2_outputs_dsc,
    name="Minimap2 consensus taxonomy classifier.",
    description=classify_consensus_minimap2_dsc,
    citations=[citations["Minimap2"], citations["li2009sequence"]],
)

plugin.methods.register_function(
    function=q2_minimap2.find_consensus_annotation,
    inputs=find_consensus_annotation_inputs,
    parameters=find_consensus_annotation_params,
    outputs=find_consensus_annotation_outputs,
    input_descriptions=find_consensus_annotation_inputs_dsc,
    parameter_descriptions=find_consensus_annotation_params_dsc,
    output_descriptions=find_consensus_annotation_outputs_dsc,
    name="Find consensus among multiple annotations.",
    description=find_consensus_annotation_dsc,
)

importlib.import_module("q2_minimap2.types._transformer")
