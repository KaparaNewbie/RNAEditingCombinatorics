import numpy as np
import pandas as pd
from Bio import Align
import parasail
from icecream import ic


def umi_seqs_overlap_early_filter_old(
    u_umi_seq,
    v_umi_seq,
    max_len_to_alignment_len_abs_diff: int = 1,
    max_gaps: int = 0,
    max_mismatches: int = 1,
):
    """Check if the UMI sequences can potentially overlap based on their lengths or an estimate of differences."""
    u_len, v_len = len(u_umi_seq), len(v_umi_seq)
    abs_len_diff = np.abs(u_len - v_len)

    # even if the shortest UMI seq of the two is fully aligned to the other,
    # we won't be able to satisfy max_len_to_alignment_len_abs_diff
    # (e.g., u_umi_seq_len = 3, v_umi_seq_len = 6, max_len_to_alignment_len_abs_diff = 1)
    if abs_len_diff / 2 > max_len_to_alignment_len_abs_diff:
        return False

    # check the worst case scenario of indels and mismatches arising from comparing
    # the two UMI sequences from the start of the longer one
    min_len = min(u_len, v_len)
    char_diff = sum(1 for i in range(min_len) if u_umi_seq[i] != v_umi_seq[i])
    if char_diff > abs_len_diff + max_gaps + max_mismatches:
        return False

    return True


def umi_seqs_overlap_early_filter(u_umi_seq, v_umi_seq, max_errors: int = 1):
    """Check if the UMI sequences can potentially overlap based on their lengths or an estimate of differences."""
    u_len, v_len = len(u_umi_seq), len(v_umi_seq)
    abs_len_diff = np.abs(u_len - v_len)

    # even if the shortest UMI seq of the two is fully aligned to the other,
    # we won't be able to satisfy max_errors
    # (e.g., u_umi_seq_len = 3, v_umi_seq_len = 6, max_errors = 1)
    if abs_len_diff / 2 > max_errors:
        return False

    # # check the worst case scenario of indels and mismatches arising from comparing
    # # the two UMI sequences from the start of the longer one
    # min_len = min(u_len, v_len)
    # char_diff = sum(1 for i in range(min_len) if u_umi_seq[i] != v_umi_seq[i])
    # if char_diff > max_errors:
    #     return False

    return True


# u_umi_seq = "AAA"
# v_umi_seq = "AAAAAA"
# assert not umi_seqs_overlap_early_filter(u_umi_seq, v_umi_seq)


# u_umi_seq = "CCAAA"
# v_umi_seq = "AAG"
# assert umi_seqs_overlap_early_filter(u_umi_seq, v_umi_seq)

# u_umi_seq = "AACCCAA"
# v_umi_seq = "ACCCA"
# assert umi_seqs_overlap_early_filter(u_umi_seq, v_umi_seq)

# u_umi_seq = "AAGGG"
# v_umi_seq = "GGG"
# assert umi_seqs_overlap_early_filter(u_umi_seq, v_umi_seq)


def umi_seqs_overlap_old(
    u,
    v,
    u_umi_seq,
    v_umi_seq,
    max_len_to_alignment_len_abs_diff: int = 1,
    max_gaps: int = 0,
    max_mismatches: int = 1,
    max_alignments: int | None = 1,
):
    umi_seqs = [u_umi_seq, v_umi_seq]

    for x_umi_seq, y_umi_seq in [umi_seqs, list(reversed(umi_seqs))]:

        x_umi_seq_len = len(x_umi_seq)
        y_umi_seq_len = len(y_umi_seq)

        aligner = Align.PairwiseAligner(
            mode="local",
            scoring="blastn",
        )

        alignments = aligner.align(x_umi_seq, y_umi_seq, strand="+")

        for alignment_i, alignment in enumerate(alignments, start=1):

            # score = alignment.score
            gaps, _, mismatches = alignment.counts()  # _ is identities
            alignment_length = alignment.length

            if (
                (gaps <= max_gaps)
                and (mismatches <= max_mismatches)
                and (
                    np.abs(x_umi_seq_len - alignment_length)
                    <= max_len_to_alignment_len_abs_diff
                )
                and (
                    np.abs(y_umi_seq_len - alignment_length)
                    <= max_len_to_alignment_len_abs_diff
                )
            ):
                return (
                    u,
                    v,
                    True,
                )  # there is a reliable overlap between the two UMI sequences

            if max_alignments is not None and alignment_i == max_alignments:
                break

    return u, v, False  # no reliable overlap found between the two UMI sequences


def umi_seqs_overlap_old_2(
    u,
    v,
    u_umi_seq,
    v_umi_seq,
    maximum_errors: int = 1,
    max_alignments: int | None = 5,
    mode: str = "local",
    debug: bool = False,
):

    u_umi_seq_len = len(u_umi_seq)
    v_umi_seq_len = len(v_umi_seq)

    if debug:
        print(f"{u_umi_seq_len=}, {v_umi_seq_len=}")

    if mode == "global":
        aligner = Align.PairwiseAligner(
            mode="global",
        )
    elif mode == "local":
        aligner = Align.PairwiseAligner(
            mode="local",
            scoring="blastn",
        )
    else:
        raise ValueError(f"Unknown mode: {mode=}")

    # minimum possible errors between the two UMI sequences observed so far
    # (we'll return this value if no reliable overlap is found, but the last alignment's errors are even
    # worse than a previous alignment's errors)
    minimal_errors = None

    alignments = aligner.align(u_umi_seq, v_umi_seq, strand="+")

    for alignment_i, alignment in enumerate(alignments, start=1):

        # score = alignment.score
        gaps, _, mismatches = alignment.counts()  # _ is identities
        # alignment_length = alignment.length

        # u_len_to_alignment_len_abs_diff = np.abs(u_umi_seq_len - alignment_length)
        # v_len_to_alignment_len_abs_diff = np.abs(v_umi_seq_len - alignment_length)

        # errors = gaps + mismatches + u_len_to_alignment_len_abs_diff + v_len_to_alignment_len_abs_diff

        errors = gaps + mismatches
        if mode == "local":
            alignment_length = alignment.length
            u_len_to_alignment_len_abs_diff = np.abs(u_umi_seq_len - alignment_length)
            v_len_to_alignment_len_abs_diff = np.abs(v_umi_seq_len - alignment_length)
            errors += u_len_to_alignment_len_abs_diff + v_len_to_alignment_len_abs_diff

        if minimal_errors is None or errors < minimal_errors:
            minimal_errors = errors

        # if debug:
        #     print(f"Alignment {alignment_i}: score = {alignment.score}, errors = {errors}\n{alignment}")

        if debug:
            print(f"Alignment {alignment_i}:")
            if mode == "global":
                print(f"{gaps=}, {mismatches=}, {errors=}")
            else:
                print(
                    f"{gaps=}, {mismatches=}, {alignment_length=}, {u_len_to_alignment_len_abs_diff=}, {v_len_to_alignment_len_abs_diff=}, {errors=}"
                )
            print(alignment)

        if errors <= maximum_errors:
            return (
                u,
                v,
                True,
                minimal_errors,
            )  # there is a reliable overlap between the two UMI sequences

        if max_alignments is not None and alignment_i == max_alignments:
            break

    return (
        u,
        v,
        False,
        minimal_errors,
    )  # no reliable overlap found between the two UMI sequences


def umi_seqs_overlap(
    u,
    v,
    u_umi_seq,
    v_umi_seq,
    # u_rtg_gene_start,
    # v_rtg_gene_start,
    # u_rtg_gene_end,
    # v_rtg_gene_end,
    # max_dist_between_two_reads_start_or_end,
    maximum_errors: int = 1,
    # matrix=parasail.matrix_create("ACGT", 1, 0),
    matrix=parasail.nuc44,
    gap_open=10,
    gap_extend=1,
    debug: bool = False,
):

    # if (
    #     np.abs(u_rtg_gene_start - v_rtg_gene_start) > max_dist_between_two_reads_start_or_end
    #     or np.abs(u_rtg_gene_end - v_rtg_gene_end) > max_dist_between_two_reads_start_or_end
    # ):
    #     return (u, v, False, np.nan) # (u, v, umis_sufficiently_similar, errors)

    u_umi_seq_len = len(u_umi_seq)
    v_umi_seq_len = len(v_umi_seq)

    # result = parasail.sw_trace(u_umi_seq, v_umi_seq, gap_open, gap_extend, matrix)
    result = parasail.nw_trace(v_umi_seq, u_umi_seq, gap_open, gap_extend, matrix)

    aligned_ref = result.traceback.ref
    aligned_query = result.traceback.query
    aligned_comp = result.traceback.comp.replace(" ", "-")

    alignment_length = len(result.traceback.ref)
    gaps = aligned_comp.count("-")
    mismatches = aligned_comp.count(".")

    # u_len_to_alignment_len_abs_diff = np.abs(u_umi_seq_len - alignment_length)
    # v_len_to_alignment_len_abs_diff = np.abs(v_umi_seq_len - alignment_length)

    errors = (
        gaps
        + mismatches
        # + u_len_to_alignment_len_abs_diff
        # + v_len_to_alignment_len_abs_diff
    )

    if debug:
        for a in [aligned_ref, aligned_comp, aligned_query]:
            print(a)

        ic(
            u_umi_seq_len,
            v_umi_seq_len,
            gaps,
            mismatches,
            alignment_length,
            # u_len_to_alignment_len_abs_diff,
            # v_len_to_alignment_len_abs_diff,
            errors,
        )

    umis_sufficiently_similar = errors <= maximum_errors

    return (u, v, umis_sufficiently_similar, errors)


# # case 1

# u_umi_seq = "AGCTCGCTAGAAACTTAG"
# v_umi_seq = "AGCTCGCTAGAAACTTA"

# umi_seqs_overlap("u", "v", u_umi_seq, v_umi_seq, debug=True, mode="global")
# umi_seqs_overlap("u", "v", u_umi_seq, v_umi_seq, debug=True, mode="local")

# assert umi_seqs_overlap_early_filter(u_umi_seq, v_umi_seq)

# # u_umi_seq_len = len(u_umi_seq)
# # v_umi_seq_len = len(v_umi_seq)

# # aligner = Align.PairwiseAligner(
# #         mode="global",
# #         # substitution_matrix = Align.substitution_matrices.load("NUC.4.4"), # https://rosalind.info/glossary/dnafull/
# #         # open_gap_score = gap_open,
# #         # extend_gap_score = gap_extend,
# #         # end_open_gap_score = end_gap_open if end_gap else 0,
# #         # end_extend_gap_score = end_gap_extend if end_gap else 0,
# #     )
# # alignments = list(aligner.align(u_umi_seq, v_umi_seq, strand="+"))
# # for i, alignment in enumerate(alignments, start=1):
# #     gaps, _, mismatches = alignment.counts()
# #     alignment_length = alignment.length
# #     u_len_to_alignment_len_abs_diff = np.abs(u_umi_seq_len - alignment_length)
# #     v_len_to_alignment_len_abs_diff = np.abs(v_umi_seq_len - alignment_length)
# #     print(f"Alignment {i}:")
# #     print(f"{gaps=}, {mismatches=}, {alignment_length=}, {u_len_to_alignment_len_abs_diff=}, {v_len_to_alignment_len_abs_diff=}")
# #     print(alignment)

# # aligner = Align.PairwiseAligner(
# #     mode="local",
# #     scoring="blastn",
# # )
# # alignments = list(aligner.align(u_umi_seq, v_umi_seq, strand="+"))
# # for i, alignment in enumerate(alignments, start=1):
# #     gaps, _, mismatches = alignment.counts()
# #     alignment_length = alignment.length
# #     u_len_to_alignment_len_abs_diff = np.abs(u_umi_seq_len - alignment_length)
# #     v_len_to_alignment_len_abs_diff = np.abs(v_umi_seq_len - alignment_length)
# #     print(f"Alignment {i}:")
# #     print(f"{gaps=}, {mismatches=}, {alignment_length=}, {u_len_to_alignment_len_abs_diff=}, {v_len_to_alignment_len_abs_diff=}")
# #     print(alignment)


# # case 2

# u_umi_seq = "GCTCGCTAGAAAGTTAG"
# v_umi_seq = "AGCTCGCTAGAAACTTA"

# umi_seqs_overlap_early_filter(u_umi_seq, v_umi_seq)

# umi_seqs_overlap("u", "v", u_umi_seq, v_umi_seq, debug=True, mode="global")
# umi_seqs_overlap("u", "v", u_umi_seq, v_umi_seq, debug=True, mode="local")


# # case 3

# u_umi_seq = "AGGTA"
# v_umi_seq = "TAGCT"

# umi_seqs_overlap_early_filter(u_umi_seq, v_umi_seq)
# umi_seqs_overlap_early_filter(u_umi_seq, v_umi_seq, max_errors=2)

# umi_seqs_overlap("u", "v", u_umi_seq, v_umi_seq, debug=True, mode="global")
# umi_seqs_overlap("u", "v", u_umi_seq, v_umi_seq, debug=True, mode="local")

# # case 3.5

# u_umi_seq = "AGCT"
# v_umi_seq = "TAGCT"

# umi_seqs_overlap("u", "v", u_umi_seq, v_umi_seq, debug=True, mode="local")
# umi_seqs_overlap("u", "v", u_umi_seq, v_umi_seq, debug=True, mode="global")

# umi_seqs_overlap_early_filter(u_umi_seq, v_umi_seq)
# umi_seqs_overlap_early_filter(u_umi_seq, v_umi_seq, max_errors=2)


# # case 4
# u_umi_seq = "AGCTCGCTAG"
# v_umi_seq = "TAGCTAGCTA"

# umi_seqs_overlap("u", "v", u_umi_seq, v_umi_seq, debug=True, mode="global")
# umi_seqs_overlap("u", "v", u_umi_seq, v_umi_seq, debug=True, mode="local")

# # case 5
# u_umi_seq = "AGCTCGCTAG"
# v_umi_seq = "TAGCTGCTA"

# umi_seqs_overlap("u", "v", u_umi_seq, v_umi_seq, debug=True, mode="global")
# umi_seqs_overlap("u", "v", u_umi_seq, v_umi_seq, debug=True, mode="local")


# gap_open = 10 # pairwise alignment score for the first residue in a gap
# gap_extend = 0.5 # pairwise alignment score for each residue in a gap after the first
# end_gap = False # apply gap penalties to the ends of sequences
# end_gap_open = 10 # score taken away when an end gap is created
# end_gap_extend = 0.5 # penalty is added to the end gap penalty for each residue in the gap (after the first?)


# some tests

# u_umi_seq = "AAAAAAAAAA"
# y_umi_seq = "AAAAAAAAGG"
# assert not umi_seqs_overlap("u", "v", u_umi_seq, y_umi_seq)[2]

# u_umi_seq = "AAAAAAAAAA"
# y_umi_seq = "AAAAAAAAAG"
# assert umi_seqs_overlap("u", "v", u_umi_seq, y_umi_seq)[2]

# u_umi_seq = "AAAAAAAAA"
# y_umi_seq = "AAAAAAAAAG"
# assert umi_seqs_overlap("u", "v", u_umi_seq, y_umi_seq)[2]

# u_umi_seq = "AACAAAAAA"
# y_umi_seq = "AAAAAAAAAG"
# assert umi_seqs_overlap("u", "v", u_umi_seq, y_umi_seq)[2]

# u_umi_seq = "ACCAAAAAAA"
# y_umi_seq = "AAAAAAAAAG"
# assert not umi_seqs_overlap("u", "v", u_umi_seq, y_umi_seq)[2]

# u_umi_seq = "ACCAAAAAA"
# y_umi_seq = "AAAAAAAAAG"
# assert not umi_seqs_overlap("u", "v", u_umi_seq, y_umi_seq)[2]


def one_batch_umi_seqs_overlap(one_batch_umi_seqs_overlap_inputs):
    return [umi_seqs_overlap(*inputs) for inputs in one_batch_umi_seqs_overlap_inputs]


def agg_one_u_overlap_vs_results_df(
    one_u_overlap_results_df: pd.DataFrame,
) -> pd.Series:
    """
    Aggregate the results for one read U from the umi_seqs_overlap_results_df DataFrame.

    `one_u_overlap_results_df` contains the columns:
    * U - the read of interest (same for all rows)
    * V - other reads that their UMI may overlap with U's
    * Overlap - boolean indicating if U and V have overlapping UMI sequences
    * MinimalErrors - minimal number of alignment errors between UMI sequences of U and V

    Return a series with:
    * Read - the read U
    * OtherReadswithIndistinguishableUMIs - list of reads V that have overlapping UMI sequences with U
    * NumOfOtherReadswithIndistinguishableUMIs - number of such reads V
    * MinimalErrorsToOtherReadswithIndistinguishableUMIs - list of minimal alignment errors to each such read V's UMI sequence
    """
    u = one_u_overlap_results_df.iloc[0]["U"]
    one_u_overlap_results_df = one_u_overlap_results_df.loc[
        one_u_overlap_results_df["Overlap"]
    ]
    return pd.Series(
        {
            "Read": u,
            "OtherReadswithIndistinguishableUMIs": one_u_overlap_results_df[
                "V"
            ].tolist(),
            "NumOfOtherReadswithIndistinguishableUMIs": one_u_overlap_results_df.shape[
                0
            ],
            "MinimalErrorsWithOtherReadswithIndistinguishableUMIs": one_u_overlap_results_df[
                "MinimalErrors"
            ].tolist(),
        }
    )


def process_one_u_overlap_vs_group(item):
    """A convenience wrapper to unpack the item tuple, when we do parallel processing over an iterable."""
    _, g = item  # _ is u
    return agg_one_u_overlap_vs_results_df(g)


def compare_u_v_editing_statuses_light(
    gene: str,
    repeat: str,
    num_of_editing_sites_in_gene: int,
    u: str,
    v: str,
    minimal_errors: int,
    used_reads_df: pd.DataFrame,
    expected_disagreements_per_position_series: pd.Series,
    used_reads_first_col_pos: int = 6,
):
    """
    Compare editing statuses of two reads u and v.

    Parameters:
    - u, v: read identifiers
    - minimal_errors: minimal number of alignment errors between u and v UMIs
    - used_reads_df: DataFrame containing editing status information for reads
    - used_reads_first_col_pos: position of the first editing status column (0-based index)

    Returns:
    A series with counts of various editing status comparisons.
    """
    u_v_reads_df = (
        used_reads_df.loc[used_reads_df["Read"].isin([u, v])]
        .set_index("Read")
        .rename(columns={"AmbigousPositions": "AmbiguousPositions"})
        .iloc[:, used_reads_first_col_pos - 1 :]
    )
    u_v_reads_df.columns = u_v_reads_df.columns.astype(int)

    u_ambiguous_positions = u_v_reads_df.loc[u, :].eq(-1)
    v_ambiguous_positions = u_v_reads_df.loc[v, :].eq(-1)

    all_ambiguous_positions = u_ambiguous_positions | v_ambiguous_positions
    shared_ambiguous_positions = u_ambiguous_positions & v_ambiguous_positions
    u_unique_ambiguous_positions = u_ambiguous_positions & ~v_ambiguous_positions
    v_unique_ambiguous_positions = ~u_ambiguous_positions & v_ambiguous_positions

    shared_unambiguous_positions = ~all_ambiguous_positions

    all_ambiguous_positions_count = all_ambiguous_positions.sum()
    shared_ambiguous_positions_count = shared_ambiguous_positions.sum()
    u_unique_ambiguous_positions_count = u_unique_ambiguous_positions.sum()
    v_unique_ambiguous_positions_count = v_unique_ambiguous_positions.sum()

    total_u_unique_ambiguous_positions_prct = (
        100 * u_unique_ambiguous_positions_count / num_of_editing_sites_in_gene
    )
    total_v_unique_ambiguous_positions_prct = (
        100 * v_unique_ambiguous_positions_count / num_of_editing_sites_in_gene
    )
    relative_u_unique_ambiguous_positions_prct = (
        100 * u_unique_ambiguous_positions_count / all_ambiguous_positions_count
    )
    relative_v_unique_ambiguous_positions_prct = (
        100 * v_unique_ambiguous_positions_count / all_ambiguous_positions_count
    )

    shared_unambiguous_positions_count = shared_unambiguous_positions.sum()

    strongly_disagreeing_positions_count = (
        u_v_reads_df.loc[:, shared_unambiguous_positions]
        .apply(lambda x: x.nunique() == 2)
        .sum()
    )
    expected_strongly_disagreeing_positions = (
        expected_disagreements_per_position_series.loc[
            shared_unambiguous_positions
        ].sum()
    )

    weakly_disagreeing_positions_count = (
        strongly_disagreeing_positions_count
        + u_v_reads_df.loc[:, all_ambiguous_positions]
        .apply(lambda x: x.nunique() == 2)
        .sum()
    )

    total_strongly_disagreeing_positions_prct = (
        100 * strongly_disagreeing_positions_count / num_of_editing_sites_in_gene
    )
    total_weakly_disagreeing_positions_prct = (
        100 * weakly_disagreeing_positions_count / num_of_editing_sites_in_gene
    )
    # % of strongly disagreeing positions out of unambiguous positions in u-v
    if shared_unambiguous_positions_count > 0:
        relative_strongly_disagreeing_positions_prct = (
            100
            * strongly_disagreeing_positions_count
            / shared_unambiguous_positions_count
        )
    else:
        relative_strongly_disagreeing_positions_prct = np.nan

    result = pd.Series(
        {
            "Gene": gene,
            "Repeat": repeat,
            "EditingSitesInGene": num_of_editing_sites_in_gene,
            "U": u,
            "V": v,
            "MinimalErrors": minimal_errors,
            "AllAmbiguousPositions": all_ambiguous_positions_count,
            "SharedAmbiguousPositions": shared_ambiguous_positions_count,
            "UUniqueAmbiguousPositions": u_unique_ambiguous_positions_count,
            "VUniqueAmbiguousPositions": v_unique_ambiguous_positions_count,
            "%TotalUUniqueAmbiguousPositions": total_u_unique_ambiguous_positions_prct,
            "%TotalVUniqueAmbiguousPositions": total_v_unique_ambiguous_positions_prct,
            "%RelativeUUniqueAmbiguousPositions": relative_u_unique_ambiguous_positions_prct,
            "%RelativeVUniqueAmbiguousPositions": relative_v_unique_ambiguous_positions_prct,
            "SharedUnambiguousPositions": shared_unambiguous_positions_count,
            "StronglyDisagreeingPositions": strongly_disagreeing_positions_count,
            "ExpectedStronglyDisagreeingPositions": expected_strongly_disagreeing_positions,
            "WeaklyDisagreeingPositions": weakly_disagreeing_positions_count,
            "%TotalStronglyDisagreeingPositions": total_strongly_disagreeing_positions_prct,
            "%TotalWeaklyDisagreeingPositions": total_weakly_disagreeing_positions_prct,
            "%RelativeStronglyDisagreeingPositions": relative_strongly_disagreeing_positions_prct,
        }
    )
    return result



def compare_u_v_noise_statuses_light(
    gene: str,
    u: str,
    v: str,
    reads_noise_df: pd.DataFrame,
    used_reads_first_col_pos: int = 1,
):
    """
    Compare editing statuses of two reads u and v.

    Parameters:
    - u, v: read identifiers
    - used_reads_df: DataFrame containing editing status information for reads
    - used_reads_first_col_pos: position of the first editing status column (0-based index)

    Returns:
    A series with counts of various editing status comparisons.
    """
    u_v_reads_df = (
        reads_noise_df.loc[reads_noise_df["Read"].isin([u, v])]
        .set_index("Read")
        .iloc[:, used_reads_first_col_pos - 1 :]
    )
    u_v_reads_df.columns = u_v_reads_df.columns.astype(int)

    u_ambiguous_positions = u_v_reads_df.loc[u, :].eq(-1)
    v_ambiguous_positions = u_v_reads_df.loc[v, :].eq(-1)

    all_ambiguous_positions = u_ambiguous_positions | v_ambiguous_positions
    shared_ambiguous_positions = u_ambiguous_positions & v_ambiguous_positions
    # u_unique_ambiguous_positions = u_ambiguous_positions & ~v_ambiguous_positions
    # v_unique_ambiguous_positions = ~u_ambiguous_positions & v_ambiguous_positions

    shared_unambiguous_positions = ~all_ambiguous_positions

    # all_ambiguous_positions_count = all_ambiguous_positions.sum()
    shared_ambiguous_positions_count = shared_ambiguous_positions.sum()
    # u_unique_ambiguous_positions_count = u_unique_ambiguous_positions.sum()
    # v_unique_ambiguous_positions_count = v_unique_ambiguous_positions.sum()

    # total_u_unique_ambiguous_positions_prct = (
    #     100 * u_unique_ambiguous_positions_count / num_of_editing_sites_in_gene
    # )
    # total_v_unique_ambiguous_positions_prct = (
    #     100 * v_unique_ambiguous_positions_count / num_of_editing_sites_in_gene
    # )
    # relative_u_unique_ambiguous_positions_prct = (
    #     100 * u_unique_ambiguous_positions_count / all_ambiguous_positions_count
    # )
    # relative_v_unique_ambiguous_positions_prct = (
    #     100 * v_unique_ambiguous_positions_count / all_ambiguous_positions_count
    # )

    shared_unambiguous_positions_count = shared_unambiguous_positions.sum()

    strongly_disagreeing_positions_count = (
        u_v_reads_df.loc[:, shared_unambiguous_positions]
        .apply(lambda x: x.nunique() == 2)
        .sum()
    )
    # expected_strongly_disagreeing_positions = (
    #     expected_disagreements_per_position_series.loc[
    #         shared_unambiguous_positions
    #     ].sum()
    # )

    # weakly_disagreeing_positions_count = (
    #     strongly_disagreeing_positions_count
    #     + u_v_reads_df.loc[:, all_ambiguous_positions]
    #     .apply(lambda x: x.nunique() == 2)
    #     .sum()
    # )

    # total_strongly_disagreeing_positions_prct = (
    #     100 * strongly_disagreeing_positions_count / num_of_editing_sites_in_gene
    # )
    # total_weakly_disagreeing_positions_prct = (
    #     100 * weakly_disagreeing_positions_count / num_of_editing_sites_in_gene
    # )
    # % of strongly disagreeing positions out of unambiguous positions in u-v
    if shared_unambiguous_positions_count > 0:
        relative_strongly_disagreeing_positions_prct = (
            100
            * strongly_disagreeing_positions_count
            / shared_unambiguous_positions_count
        )
    else:
        relative_strongly_disagreeing_positions_prct = np.nan

    result = pd.Series(
        {
            "Gene": gene,
            # "Repeat": repeat,
            # "EditingSitesInGene": num_of_editing_sites_in_gene,
            "U": u,
            "V": v,
            # "MinimalErrors": minimal_errors,
            # "AllAmbiguousPositions": all_ambiguous_positions_count,
            "SharedAmbiguousNoisePositions": shared_ambiguous_positions_count,
            # "UUniqueAmbiguousPositions": u_unique_ambiguous_positions_count,
            # "VUniqueAmbiguousPositions": v_unique_ambiguous_positions_count,
            # "%TotalUUniqueAmbiguousPositions": total_u_unique_ambiguous_positions_prct,
            # "%TotalVUniqueAmbiguousPositions": total_v_unique_ambiguous_positions_prct,
            # "%RelativeUUniqueAmbiguousPositions": relative_u_unique_ambiguous_positions_prct,
            # "%RelativeVUniqueAmbiguousPositions": relative_v_unique_ambiguous_positions_prct,
            "SharedUnambiguousNoisePositions": shared_unambiguous_positions_count,
            "StronglyDisagreeingNoisePositions": strongly_disagreeing_positions_count,
            # "ExpectedStronglyDisagreeingPositions": expected_strongly_disagreeing_positions,
            # "WeaklyDisagreeingPositions": weakly_disagreeing_positions_count,
            # "%TotalStronglyDisagreeingPositions": total_strongly_disagreeing_positions_prct,
            # "%TotalWeaklyDisagreeingPositions": total_weakly_disagreeing_positions_prct,
            "%RelativeStronglyDisagreeingNoisePositions": relative_strongly_disagreeing_positions_prct,
        }
    )
    return result