# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.3
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Setup

# %%
import random

from Bio import Align
import parasail
from icecream import ic
import numpy as np


# %% [markdown]
# # Simulate data

# %%
def random_nucleotide_seq(min_len=10, max_len=14):
    length = random.randint(min_len, max_len)
    return ''.join(random.choices('ACGT', k=length))


# %%
def mutate_nucleotide_sequence(seq, x):
    """Mutate up to x random nucleotides in the given DNA/RNA sequence."""
    original_nucleotides = list(seq)
    indices = random.sample(range(len(original_nucleotides)), x)
    for idx in indices:
        new_nucleotide = random.choice(list("ATCG")) # Choose a new nucleotide, possibly the same
        original_nucleotides[idx] = new_nucleotide

    mutated_seq = ''.join(original_nucleotides)
    return mutated_seq

# Example usage:
seq1 = random_nucleotide_seq()
# ic(seq1)
x = 3  # number of nucleotides to mutate
mutated_seq2 = mutate_nucleotide_sequence(seq1, x)
print("Original:", seq1)
print("Mutated :", mutated_seq2)

# %%
# n_base_seqs = 100
n_base_seqs = 600
max_diverging_seqs_per_base_seq = int(n_base_seqs ** 0.5)
max_mutations_per_diverging_seq = 5
max_comparisons_per_seq = 50

# %%
base_seqs = [random_nucleotide_seq() for _ in range(n_base_seqs)]
base_seqs

# %%
diverging_seqs = [
    mutate_nucleotide_sequence(base_seq, max_mutations_per_diverging_seq)
    for base_seq in base_seqs
    for _ in range(random.randint(1, max_diverging_seqs_per_base_seq))
]
ic(len(diverging_seqs));

# %%
seqs = base_seqs + diverging_seqs
ic(len(seqs));
random.shuffle(seqs)
seqs

# %%
N = len(seqs)

# %%
comparisons_indices = [
    (i, j)
    for i in range(N)
    for j in random.choices(range(0, N), k=random.randint(1, max_comparisons_per_seq))
]
ic(len(comparisons_indices))
comparisons_indices

# %%
comparisons = [
    (seqs[i], seqs[j])
    for i, j in comparisons_indices
]
comparisons

# %%
len(comparisons)


# %% [markdown]
# # Alignment functions

# %%
def python_pairwise_align(seq1, seq2):
    aligner = Align.PairwiseAligner(
        mode="local",
        # substitution_matrix = Align.substitution_matrices.load("NUC.4.4")
    )
    alignments = aligner.align(seq1, seq2, "+")
    return alignments


# %%
type(parasail.matrix_create("ACGT", 1, 0))


# %%
def parasail_pairwise_align(
    seq1, seq2, 
    # matrix=parasail.matrix_create("ACGT", 1, 0),
    # gap_open=0,
    # gap_extend=0,
    matrix=parasail.nuc44,
    gap_open=-10,
    gap_extend=-0.5
):
    # matrix = parasail.nuc44
    # matrix = parasail.matrix_create("ACGT", 1, 0)
    result = parasail.sw_trace(seq2, seq1, gap_open, gap_extend, matrix)
    
    return result


# %% [markdown]
# # Parsing tests

# %% [markdown]
# ## Parsing test data

# %%
# def get_u_and_v_with_at_least_one_mismatch(trials=1000):
    
#     for _ in range(trials):
    
#         u_umi_seq = v_umi_seq = ""

#         while u_umi_seq == v_umi_seq:

#             u_umi_seq = random_nucleotide_seq()
#             v_umi_seq = mutate_nucleotide_sequence(u_umi_seq, 3)

#             # u_umi_seq_len = len(u_umi_seq)
#             # v_umi_seq_len = len(v_umi_seq)

#         # ic(u_umi_seq, v_umi_seq);
        
#         alignments = list(python_pairwise_align(u_umi_seq, v_umi_seq))
#         # ic(len(alignments));
        
#         for alignment in alignments:
#             _, _, mismatches = alignment.counts()
#             if mismatches > 0:
#                 # return u_umi_seq, v_umi_seq
#                 ic(u_umi_seq, v_umi_seq);
#                 return

# %%
# # u_umi_seq, v_umi_seq = comparisons[4]
# u_umi_seq, v_umi_seq = comparisons[6]

# u_umi_seq, v_umi_seq

# %%
# # u_umi_seq, v_umi_seq = comparisons[4]
# u_umi_seq, v_umi_seq = comparisons[8]

# u_umi_seq, v_umi_seq

# %%
# get_u_and_v_with_at_least_one_mismatch()

# %%
# # u_umi_seq = "GCTCGCTAGAAAGTTAG"
# # v_umi_seq = "AGCTCGCTAGAAACTTA"

# u_umi_seq = "GCTCCCTAGAAAGTTAG"
# v_umi_seq = "AGCTCGCTAGAAACTTA"

# ic(
#     u_umi_seq, 
#     v_umi_seq
# );

# %%
u_umi_seq = "AAAAGGGAGGGAAAAC"
v_umi_seq = "AAAAGGGCGGGAAAA"

# %%
# u_umi_seq = 'ACCCCCCTACAT'
# v_umi_seq = 'ACGAGTTTACAT'

# u_umi_seq = 'ACCCCCCTACATG'
# v_umi_seq = 'ACGAGTTTACAT'

# %%
u_umi_seq_len = len(u_umi_seq)
v_umi_seq_len = len(v_umi_seq)

ic(
    u_umi_seq_len,
    v_umi_seq_len
);

# %% [markdown]
# ## Biopython

# %%
alignments = list(python_pairwise_align(u_umi_seq, v_umi_seq))
ic(len(alignments));

# %%
for i, alignment in enumerate(alignments):
    gaps, _, mismatches = alignment.counts()
    alignment_length = alignment.length
    print(f"Alignment {i}: {alignment_length = }, {gaps = }, {mismatches = }")
    print(alignment)
    print()

# %%
alignment = alignments[2]
print(alignment)

# %%
gaps, _, mismatches = alignment.counts()  # _ is identities
alignment_length = alignment.length

u_len_to_alignment_len_abs_diff = np.abs(u_umi_seq_len - alignment_length)
v_len_to_alignment_len_abs_diff = np.abs(v_umi_seq_len - alignment_length)

errors = gaps + mismatches + u_len_to_alignment_len_abs_diff + v_len_to_alignment_len_abs_diff

ic(
    u_umi_seq_len,
    v_umi_seq_len,
    gaps,
    mismatches,
    alignment_length,
    u_len_to_alignment_len_abs_diff,
    v_len_to_alignment_len_abs_diff,
    errors,
);

# %% [markdown]
# ## Parasail

# %%
# for i, alignment in enumerate(alignments):
#     gaps, _, mismatches = alignment.counts()
#     alignment_length = alignment.length
#     print(f"Alignment {i}: {alignment_length = }, {gaps = }, {mismatches = }")
#     print(alignment)
#     print()

# %%
# gaps, _, mismatches = alignment.counts()  # _ is identities
# alignment_length = alignment.length

# u_len_to_alignment_len_abs_diff = np.abs(u_umi_seq_len - alignment_length)
# v_len_to_alignment_len_abs_diff = np.abs(v_umi_seq_len - alignment_length)

# errors = gaps + mismatches + u_len_to_alignment_len_abs_diff + v_len_to_alignment_len_abs_diff

# ic(
#     u_umi_seq_len,
#     v_umi_seq_len,
#     gaps,
#     mismatches,
#     alignment_length,
#     u_len_to_alignment_len_abs_diff,
#     v_len_to_alignment_len_abs_diff,
#     errors,
# );

# %%
matrix = parasail.nuc44
gap_open = 10
gap_extend = 1

# %%
u_umi_seq, v_umi_seq

# %%
# # result = parasail.sw_trace(v_umi_seq, u_umi_seq, 0, 0, parasail.nuc44)
# # result = parasail_pairwise_align(u_umi_seq, v_umi_seq)
# result = parasail.sw_trace(v_umi_seq, u_umi_seq, gap_open, gap_extend, matrix)

# aligned_ref = result.traceback.ref
# aligned_query = result.traceback.query
# aligned_comp = result.traceback.comp.replace(" ", "-")

# alignment_length = len(result.traceback.ref)
# gaps = aligned_comp.count("-")
# mismatches = aligned_comp.count(".")

# u_len_to_alignment_len_abs_diff = np.abs(u_umi_seq_len - alignment_length)
# v_len_to_alignment_len_abs_diff = np.abs(v_umi_seq_len - alignment_length)

# errors = gaps + mismatches + u_len_to_alignment_len_abs_diff + v_len_to_alignment_len_abs_diff

# print(f"{alignment_length = }")
# print()

# for a in [aligned_ref, aligned_comp, aligned_query]:
#     print(a)
    
# ic(
#     u_umi_seq_len,
#     v_umi_seq_len,
#     gaps,
#     mismatches,
#     alignment_length,
#     u_len_to_alignment_len_abs_diff,
#     v_len_to_alignment_len_abs_diff,
#     errors,
# );

# %%
# dir(result)

# %%
# result.ref

# %%
# result.query

# %%
# u_end = result.end_ref + 1
# v_end = result.end_query + 1

# u_end, v_end

# %%
# u_start = u_end - len(aligned_ref.replace("-", ""))
# v_start = v_end - len(aligned_query.replace("-", ""))

# u_start, v_start

# %%
seq_couples = [
    ["GCTCATCAATA", "GCTCATCAATA"],
    ["GCTCATCAATA", "GCTCGCAATAC"],
    ["GCTCATCAATA", "GCCCCTCAATAC"],
    ["GCTCATCAATA", "CAAGGCTTCAATA"],
    ["GCTCATCAATA", "TCCTCCAATAAA"],
]

# %%

# %%
for i, (u_umi_seq, v_umi_seq) in enumerate(seq_couples):

    u_umi_seq_len = len(u_umi_seq)
    v_umi_seq_len = len(v_umi_seq)
    
    result = parasail.nw_trace(v_umi_seq, u_umi_seq, gap_open, gap_extend, matrix)

    aligned_ref = result.traceback.ref
    aligned_query = result.traceback.query
    aligned_comp = result.traceback.comp.replace(" ", "-")

    alignment_length = len(result.traceback.ref)
    gaps = aligned_comp.count("-")
    mismatches = aligned_comp.count(".")

    # u_len_to_alignment_len_abs_diff = np.abs(u_umi_seq_len - alignment_length)
    # v_len_to_alignment_len_abs_diff = np.abs(v_umi_seq_len - alignment_length)

    # errors = gaps + mismatches + u_len_to_alignment_len_abs_diff + v_len_to_alignment_len_abs_diff

    errors = gaps + mismatches

    print(f"{i = }")
    print(f"{u_umi_seq = }\n{v_umi_seq = }")
    for a in [aligned_ref, aligned_comp, aligned_query]:
        print(a)
        
    # ic(
    #     u_umi_seq_len,
    #     v_umi_seq_len,
    #     gaps,
    #     mismatches,
    #     alignment_length,
    #     # u_len_to_alignment_len_abs_diff,
    #     # v_len_to_alignment_len_abs_diff,
    #     errors,
    # );
    
    print(
        f"Length of UMI sequences: {u_umi_seq_len}, {v_umi_seq_len}\n"
        f"Alignment length: {alignment_length}\n"
        f"Gaps: {gaps}\n"
        f"Mismatches: {mismatches}\n"
        f"Total errors: {errors}"
    )
    print()

# %%

# %%

# %%

# %% [markdown]
# # Alignment benchmarks

# %%
len(comparisons)

# %%
# %%timeit

for seq1, seq2 in comparisons:
    python_pairwise_align(seq1, seq2)

# %%
# %%timeit

for seq1, seq2 in comparisons:
    parasail_pairwise_align(seq1, seq2)

# %%
test_input_size = 197_240
test_time_secs = 2.14 # seconds
test_time_secs_per_input = test_time_secs / test_input_size
test_time_secs_per_input

# %%
real_input_size = 53_591_452
estimated_real_time_secs = real_input_size * test_time_secs_per_input
estimated_real_time_hrs = estimated_real_time_secs / 3600
estimated_real_time_hrs
