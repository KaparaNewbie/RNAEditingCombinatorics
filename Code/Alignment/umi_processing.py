import numpy as np
from Bio import Align


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


def umi_seqs_overlap(
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
