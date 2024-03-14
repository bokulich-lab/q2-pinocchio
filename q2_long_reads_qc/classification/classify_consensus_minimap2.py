# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import time

import pandas as pd
from qiime2.plugin import (  # Bool,; Choices,; Float,; Int,; Range,; Str,; Threads,
    get_available_cores,
)

# Specify default settings for various functions
DEFAULTMAXACCEPTS = 10
DEFAULTPERCENTID = 0.8
DEFAULTQUERYCOV = 0.8
DEFAULTSTRAND = "both"
DEFAULTEVALUE = 0.001
DEFAULTMINCONSENSUS = 0.51
DEFAULTOUTPUTNOHITS = True
DEFAULTNUMTHREADS = 1
DEFAULTUNASSIGNABLELABEL = "Unassigned"


def _PairwiseAlignmentMN2_format_df_to_series_of_lists(
    assignments: pd.DataFrame,
    ref_taxa: pd.Series,
    unassignable_label: str = DEFAULTUNASSIGNABLELABEL,
) -> pd.Series:
    """import observed assignments in blast6 format to series of lists.

    assignments: pd.DataFrame
        Taxonomy observation map in blast format 6. Each line consists of
        taxonomy assignments of a query sequence in tab-delimited format:
            <query_id>    <subject-seq-id>   <...other columns are ignored>

    ref_taxa: pd.Series
        Reference taxonomies in tab-delimited format:
            <accession ID>  Annotation
        The accession IDs in this taxonomy should match the subject-seq-ids in
        the "assignment" input.
    """
    # validate that assignments are present in reference taxonomy
    # (i.e., that the correct reference taxonomy was used).
    # Note that we drop unassigned labels from this set.

    print(assignments["sseqid"].values)
    time.sleep(1000)
    missing_ids = set(assignments["sseqid"].values) - set(ref_taxa.index) - {"*", ""}
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
        sseqid = assignments_copy.iloc[index]["sseqid"]
        assignments_copy.at[index, "sseqid"] = ref_taxa.at[sseqid]
    # convert to dict of {accession_id: [annotations]}
    taxa_hits: pd.Series = assignments_copy.set_index("qseqid")["sseqid"]
    taxa_hits = taxa_hits.groupby(taxa_hits.index).apply(list)

    return taxa_hits


def find_consensus_annotation(
    search_results: pd.DataFrame,
    reference_taxonomy: pd.Series,
    min_consensus: int = 0.51,
    unassignable_label: str = DEFAULTUNASSIGNABLELABEL,
) -> pd.Series:

    # load and convert blast6format results to dict of taxa hits
    # obs_taxa = _PairwiseAlignmentMN2_format_df_to_series_of_lists(
    #    search_results, reference_taxonomy, unassignable_label=unassignable_label
    # )

    # compute consensus annotations

    """
    result = _compute_consensus_annotations(
        obs_taxa, min_consensus=min_consensus,
        unassignable_label=unassignable_label)
    result.index.name = 'Feature ID'
    """

    return search_results


def classify_consensus_minimap2(
    ctx,
    query,
    reference_taxonomy,
    minimap2_index=None,
    reference_reads=None,
    maxaccepts=DEFAULTMAXACCEPTS,
    perc_identity=DEFAULTPERCENTID,
    query_cov=DEFAULTQUERYCOV,
    strand=DEFAULTSTRAND,
    evalue=DEFAULTEVALUE,
    output_no_hits=DEFAULTOUTPUTNOHITS,
    min_consensus=DEFAULTMINCONSENSUS,
    unassignable_label=DEFAULTUNASSIGNABLELABEL,
    num_threads=DEFAULTNUMTHREADS,
):
    if num_threads == 0:
        num_threads = get_available_cores()

    search_db = ctx.get_action("long_reads_qc", "minimap2_search")
    # lca = ctx.get_action('long_reads_qc', 'find_consensus_annotation')
    (result,) = search_db(
        query_reads=query,
        minimap2_index=minimap2_index,
        reference_reads=reference_reads,
        maxaccepts=maxaccepts,
        perc_identity=perc_identity,
        output_no_hits=output_no_hits,
        n_threads=num_threads,
    )

    # consensus, = lca(search_results=result,
    #                 reference_taxonomy=reference_taxonomy,
    #                 min_consensus=min_consensus,
    #                 unassignable_label=unassignable_label)

    # New: add BLAST6Format result as an output. This could just as well be a
    # visualizer generated from these results (using q2-metadata tabulate).
    # Would that be more useful to the user that the QZA?
    return result  # , consensus
