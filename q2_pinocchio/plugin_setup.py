# ----------------------------------------------------------------------------
# Copyright (c) 2024, Bokulich Lab.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import importlib

from q2_types.feature_data import FeatureData
from qiime2.plugin import Citations, Plugin

import q2_pinocchio
from q2_pinocchio import __version__
from q2_pinocchio._action_params import (
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
    extract_reads_dsc,
    extract_reads_inputs,
    extract_reads_inputs_dsc,
    extract_reads_outputs,
    extract_reads_outputs_dsc,
    extract_reads_param_dsc,
    extract_reads_params,
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
    minimap2_search_dsc,
    minimap2_search_inputs,
    minimap2_search_inputs_dsc,
    minimap2_search_outputs,
    minimap2_search_outputs_dsc,
    minimap2_search_param_dsc,
    minimap2_search_params,
    stats_dsc,
    stats_input_descriptions,
    stats_inputs,
    trim_dsc,
    trim_input_descriptions,
    trim_inputs,
    trim_output_descriptions,
    trim_outputs,
    trim_parameter_descriptions,
    trim_parameters,
)
from q2_pinocchio.types._format import (
    Minimap2IndexDBDirFmt,
    Minimap2IndexDBFmt,
    PairwiseAlignmentMN2DirectoryFormat,
    PairwiseAlignmentMN2Format,
)
from q2_pinocchio.types._type import Minimap2IndexDB, PairwiseAlignmentMN2

citations = Citations.load("citations.bib", package="q2_pinocchio")

plugin = Plugin(
    name="pinocchio",
    version=__version__,
    website="https://github.com/bokulich-lab/q2-pinocchio",
    package="q2_pinocchio",
    description="QIIME 2 Plugin for quality control and taxonomic "
    "classification of long read sequences.",
    short_description="",
)

plugin.register_formats(
    Minimap2IndexDBDirFmt,
    Minimap2IndexDBFmt,
    PairwiseAlignmentMN2Format,
    PairwiseAlignmentMN2DirectoryFormat,
)
plugin.register_semantic_types(Minimap2IndexDB, PairwiseAlignmentMN2)

plugin.register_artifact_class(
    Minimap2IndexDB,
    directory_format=Minimap2IndexDBDirFmt,
    description=(
        "Represents a Minimap2 index database used for efficient sequence "
        "alignment by storing pre-processed genome or reference data. "
        "Ideal when we have repeated use of the same reference file."
    ),
)

plugin.register_artifact_class(
    FeatureData[PairwiseAlignmentMN2],
    directory_format=PairwiseAlignmentMN2DirectoryFormat,
    description=(
        "Represents the PAF (Pairwise mApping Format) file generated by Minimap2 "
        "genome mapper. This format is used to represent the alignment of query "
        "sequences against a reference sequence, detailing regions of similarity and "
        "providing insights into structural variations, sequence identity, "
        "and other genomic features."
    ),
)

plugin.methods.register_function(
    function=q2_pinocchio.extract_reads,
    inputs=extract_reads_inputs,
    parameters=extract_reads_params,
    outputs=extract_reads_outputs,
    input_descriptions=extract_reads_inputs_dsc,
    parameter_descriptions=extract_reads_param_dsc,
    output_descriptions=extract_reads_outputs_dsc,
    name="Extract long feature sequences using Minimap2.",
    description=extract_reads_dsc,
    citations=[citations["Minimap2"]],
)

plugin.methods.register_function(
    function=q2_pinocchio.filter_reads,
    inputs=filter_reads_inputs,
    parameters=filter_reads_params,
    outputs=filter_reads_outputs,
    input_descriptions=filter_reads_inputs_dsc,
    parameter_descriptions=filter_reads_param_dsc,
    output_descriptions=filter_reads_outputs_dsc,
    name="Filter long demultiplexed sequences using Minimap2.",
    description=filter_reads_dsc,
    citations=[citations["Minimap2"]],
)

plugin.methods.register_function(
    function=q2_pinocchio.build_index,
    inputs=build_index_inputs,
    parameters=build_index_params,
    outputs=build_index_outputs,
    input_descriptions=build_index_inputs_dsc,
    parameter_descriptions=build_index_param_dsc,
    output_descriptions=build_index_outputs_dsc,
    name="Build Minimap2 index database.",
    description=build_index_dsc,
    citations=[citations["Minimap2"]],
)

plugin.methods.register_function(
    function=q2_pinocchio.minimap2_search,
    inputs=minimap2_search_inputs,
    parameters=minimap2_search_params,
    outputs=minimap2_search_outputs,
    input_descriptions=minimap2_search_inputs_dsc,
    parameter_descriptions=minimap2_search_param_dsc,
    output_descriptions=minimap2_search_outputs_dsc,
    name="Minimap2 alignment search.",
    description=minimap2_search_dsc,
    citations=[citations["Minimap2"], citations["li2009sequence"]],
)

plugin.pipelines.register_function(
    function=q2_pinocchio.classify_consensus_minimap2,
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
    function=q2_pinocchio._find_consensus_annotation,
    inputs=find_consensus_annotation_inputs,
    parameters=find_consensus_annotation_params,
    outputs=find_consensus_annotation_outputs,
    input_descriptions=find_consensus_annotation_inputs_dsc,
    parameter_descriptions=find_consensus_annotation_params_dsc,
    output_descriptions=find_consensus_annotation_outputs_dsc,
    name="Find consensus among multiple annotations.",
    description=find_consensus_annotation_dsc,
)

plugin.visualizers.register_function(
    function=q2_pinocchio.stats,
    inputs=stats_inputs,
    parameters="",
    input_descriptions=stats_input_descriptions,
    parameter_descriptions={},
    name="Quality control statistics for long sequences using NanoPlot.",
    description=stats_dsc,
    citations=[citations["Nanopack2"]],
)

plugin.methods.register_function(
    function=q2_pinocchio.trim,
    inputs=trim_inputs,
    outputs=trim_outputs,
    parameters=trim_parameters,
    input_descriptions=trim_input_descriptions,
    output_descriptions=trim_output_descriptions,
    parameter_descriptions=trim_parameter_descriptions,
    name="Trim long sequences using Chopper.",
    description=trim_dsc,
    citations=[citations["Nanopack2"]],
)

importlib.import_module("q2_pinocchio.types._transformer")
