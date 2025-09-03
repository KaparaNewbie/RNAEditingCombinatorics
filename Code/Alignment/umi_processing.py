import numpy as np
from Bio import Align


def umi_seqs_overlap_early_filter(
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


# def one_batch_umi_seqs_overlap_early_filter(one_batch_umi_seqs_overlap_inputs):
#     return [
#         input
#         for input in one_batch_umi_seqs_overlap_inputs
#         if umi_seqs_overlap_early_filter(
#             input[2],  # u_umi_seq
#             input[3],  # v_umi_seq
#             input[4],  # max_len_to_alignment_len_abs_diff
#             input[5],  # max_gaps
#             input[6],  # max_mismatches
#         )
#     ]


u_umi_seq = "AAA"
v_umi_seq = "AAAAAA"
assert not umi_seqs_overlap_early_filter(u_umi_seq, v_umi_seq)


u_umi_seq = "CCAAA"
v_umi_seq = "AAG"
assert umi_seqs_overlap_early_filter(u_umi_seq, v_umi_seq)

u_umi_seq = "AACCCAA"
v_umi_seq = "ACCCA"
assert umi_seqs_overlap_early_filter(u_umi_seq, v_umi_seq)

u_umi_seq = "AAGGG"
v_umi_seq = "GGG"
assert umi_seqs_overlap_early_filter(u_umi_seq, v_umi_seq)


def umi_seqs_overlap(
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


# some tests

x_umi_seq = "AAAAAAAAAA"
y_umi_seq = "AAAAAAAAGG"
assert not umi_seqs_overlap("u", "v", x_umi_seq, y_umi_seq)[2]

x_umi_seq = "AAAAAAAAAA"
y_umi_seq = "AAAAAAAAAG"
assert umi_seqs_overlap("u", "v", x_umi_seq, y_umi_seq)[2]

x_umi_seq = "AAAAAAAAA"
y_umi_seq = "AAAAAAAAAG"
assert umi_seqs_overlap("u", "v", x_umi_seq, y_umi_seq)[2]

x_umi_seq = "AACAAAAAA"
y_umi_seq = "AAAAAAAAAG"
assert umi_seqs_overlap("u", "v", x_umi_seq, y_umi_seq)[2]

x_umi_seq = "ACCAAAAAAA"
y_umi_seq = "AAAAAAAAAG"
assert not umi_seqs_overlap("u", "v", x_umi_seq, y_umi_seq)[2]

x_umi_seq = "ACCAAAAAA"
y_umi_seq = "AAAAAAAAAG"
assert not umi_seqs_overlap("u", "v", x_umi_seq, y_umi_seq)[2]


def one_batch_umi_seqs_overlap(one_batch_umi_seqs_overlap_inputs):
    return [umi_seqs_overlap(*inputs) for inputs in one_batch_umi_seqs_overlap_inputs]
