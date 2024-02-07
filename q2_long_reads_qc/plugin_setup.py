# ----------------------------------------------------------------------------
# Copyright (c) 2022, <developer name>.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from q2_types.per_sample_sequences import (
    PairedEndSequencesWithQuality,
    SequencesWithQuality,
)
from q2_types.sample_data import SampleData
from qiime2.plugin import Citations, Plugin

import q2_long_reads_qc
from q2_long_reads_qc import __version__

citations = Citations.load("citations.bib", package="q2_long_reads_qc")

plugin = Plugin(
    name="long-reads-qc",
    version=__version__,
    website="https://github.com/bokulich-lab/q2-plugin-name",
    package="q2_long_reads_qc",
    description="QIIME 2 Plugin for quality control of long read sequences.",
    short_description="",
)


plugin.visualizers.register_function(
    function=q2_long_reads_qc.longqc.evaluate.evaluate_long_reads,
    inputs={
        "seqs": SampleData[SequencesWithQuality | PairedEndSequencesWithQuality],
    },
    parameters="",
    input_descriptions={
        "seqs": "Long reads to be analyzed.",
    },
    parameter_descriptions="",
    name="Evaluate long sequences using LongQC.",
    description="This method uses LongQC "
    " to assess the quality of long reads providing"
    " visualizations summarizing the results.",
    citations=[citations["Fukasawa_LongQC"]],
)
