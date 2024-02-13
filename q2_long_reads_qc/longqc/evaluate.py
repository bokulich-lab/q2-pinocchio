# ----------------------------------------------------------------------------
# Copyright (c) 2024, Christos Matzoros.
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


def evaluate_long_reads(
    output_dir: str,
    sequences: SampleData[SequencesWithQuality | PairedEndSequencesWithQuality],
) -> None:

    pass
