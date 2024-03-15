# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
from q2_feature_classifier._consensus_assignment import _compute_consensus_annotations


def _PairwiseAlignmentMN2_format_df_to_series_of_lists(
    assignments: pd.DataFrame,
    ref_taxa: pd.Series,
    unassignable_label: str = "Unassigned",
) -> pd.Series:

    missing_ids = set(assignments[5].values) - set(ref_taxa.index) - {"*", ""}
    if len(missing_ids) > 0:
        raise KeyError(
            "Reference taxonomy and search results do not match. "
            "The following identifiers were reported in the search "
            "results but are not present in the reference taxonomy:"
            " {0}".format(", ".join(str(i) for i in missing_ids))
        )

    # if vsearch fails to find assignment, it reports '*' as the
    # accession ID, so we will add this mapping to the reference taxonomy.
    ref_taxa["*"] = unassignable_label
    assignments_copy = assignments.copy(deep=True)
    for index, value in assignments_copy.iterrows():
        sseqid = assignments_copy.iloc[index][5]
        assignments_copy.at[index, 5] = ref_taxa.at[sseqid]
    # convert to dict of {accession_id: [annotations]}
    taxa_hits: pd.Series = assignments_copy.set_index(0)[5]
    taxa_hits = taxa_hits.groupby(taxa_hits.index).apply(list)

    return taxa_hits


def find_consensus_annotation(
    search_results: pd.DataFrame,
    reference_taxonomy: pd.Series,
    min_consensus: int = 0.51,
    unassignable_label: str = "Unassigned",
) -> pd.DataFrame:

    # load and convert PairwiseAlignmentFormat (PAF) results to dict of taxa hits
    obs_taxa = _PairwiseAlignmentMN2_format_df_to_series_of_lists(
        search_results, reference_taxonomy, unassignable_label=unassignable_label
    )

    # compute consensus annotations
    result = _compute_consensus_annotations(
        obs_taxa, min_consensus=min_consensus, unassignable_label=unassignable_label
    )
    result.index.name = "Feature ID"

    return result


def classify_consensus(
    ctx,
    query,
    reference_taxonomy,
    index_database=None,
    reference_reads=None,
    maxaccepts=10,
    perc_identity=0.7,
    output_no_hits=True,
    min_consensus=0.51,
    unassignable_label="Unassigned",
    num_threads=1,
):
    # Retrieve the actions for mapping and consensus annotation
    search_db = ctx.get_action("long_reads_qc", "minimap2")
    lca = ctx.get_action("long_reads_qc", "find_consensus_annotation")

    # Search for top hits in a reference database
    (result,) = search_db(
        query_reads=query,
        index_database=index_database,
        reference_reads=reference_reads,
        maxaccepts=maxaccepts,
        perc_identity=perc_identity,
        output_no_hits=output_no_hits,
        n_threads=num_threads,
    )

    # Find consensus annotation for each query searched against a
    # reference database
    (consensus,) = lca(
        search_results=result,
        reference_taxonomy=reference_taxonomy,
        min_consensus=min_consensus,
        unassignable_label=unassignable_label,
    )

    return result, consensus
