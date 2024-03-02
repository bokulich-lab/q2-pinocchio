# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import Bool, Choices, Float, Int, Range, Str

filter_reads_params = {
    "n_threads": Int % Range(1, None),
    "mapping_preset": Str % Choices(["map-ont", "map-hifi", "map-pb"]),
    "exclude_mapped": Bool,
    "min_per_identity": Float % Range(0.0, 1.0, inclusive_end=True),
    "matching_score": Int,
    "mismatching_penalty": Int,
    "gap_open_penalty": Int % Range(1, None),
    "gap_extension_penalty": Int % Range(1, None),
}

filter_reads_param_descriptions = {
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
    "reclassified as unmapped if their identity percentage falls below this "
    "value. If not set there is no reclassification.",
    "matching_score": "Matching score.",
    "mismatching_penalty": "Mismatching penalty.",
    "gap_open_penalty": "Gap open penalty.",
    "gap_extension_penalty": "Gap extension penalty.",
}

minimap2_build_params = {"kmer_length": Int % Range(1, 28)}
minimap2_build_param_descriptions = {"kmer_length": "Minimizer k-mer length."}
