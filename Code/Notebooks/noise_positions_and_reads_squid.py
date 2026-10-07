# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.4
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Imports

# %%
# code_dir = "/private7/projects/Combinatorics/Code"

# %%
# from itertools import chain, product
from itertools import chain
from pathlib import Path
# import sys
from multiprocessing import Pool

from icecream import ic
import plotly.express as px
# import plotly.graph_objects as go
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact, chi2_contingency
from statsmodels.stats.multitest import fdrcorrection

# sys.path.append(str(Path(code_dir).absolute()))
# import 
# from Pileup.positions import read_positions_file

# %% [markdown]
# # Inputs

# %%
seed = 1892
edited_col = "Edited"
new_noise_positions_in_use = True
snp_noise_level = 0.1
final_noise_col: str | None = None
multisample = False
top_x_noisy_positions = 3
assurance_factor = 1.5

# %%
# old_reads_first_pos_loc = 6

# %%
# # if noise positions were re-created where non ref base positions with noise >= 0.1 weren't discarded.
# # in that case - the new reads files may have different new reads names,
# # so to connect between the reads with editing statuses and reads with new noise statuses,
# # we need to use to old reads col which is another col in the new reads df.
# new_noise_positions_in_use = True

# %%
# pacbio_samples = ["GRIA2", "PCLO", "ADAR1", "IQEC1"]

# pacbio_new_positions_files = [
#     "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/CompleteNoisePositions/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.positions.csv.gz",
#     "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/CompleteNoisePositions/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.positions.csv.gz",
#     "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/CompleteNoisePositions/ADAR1.Merged.r64296e203404D01.aligned.sorted.MinRQ998.positions.csv.gz",
#     "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/CompleteNoisePositions/IQEC.Merged.r64296e203404D01.aligned.sorted.MinRQ998.positions.csv.gz",
# ]
# pacbio_new_old_to_new_reads_files = [
#     "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/CompleteNoisePositions/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.OldToNewReads.csv.gz",
#     "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/CompleteNoisePositions/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.OldToNewReads.csv.gz",
#     "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/CompleteNoisePositions/ADAR1.Merged.r64296e203404D01.aligned.sorted.MinRQ998.OldToNewReads.csv.gz",
#     "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/CompleteNoisePositions/IQEC.Merged.r64296e203404D01.aligned.sorted.MinRQ998.OldToNewReads.csv.gz",
# ]

# pacbio_old_reads_files = [
#     "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.reads.csv.gz",
#     "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.reads.csv.gz",
#     "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples/ADAR1.Merged.r64296e203404D01.aligned.sorted.MinRQ998.reads.csv.gz",
#     "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples/IQEC.Merged.r64296e203404D01.aligned.sorted.MinRQ998.reads.csv.gz",
# ]
# pacbio_old_unique_reads_files = [
#     "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_reads.csv.gz",
#     "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_reads.csv.gz",
#     "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples/ADAR1.Merged.r64296e203404D01.aligned.sorted.MinRQ998.unique_reads.csv.gz",
#     "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples/IQEC.Merged.r64296e203404D01.aligned.sorted.MinRQ998.unique_reads.csv.gz",
# ]
# pacbio_old_old_to_new_reads_files = [
#     "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.OldToNewReads.csv.gz",
#     "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.OldToNewReads.csv.gz",
#     "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples/ADAR1.Merged.r64296e203404D01.aligned.sorted.MinRQ998.OldToNewReads.csv.gz",
#     "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples/IQEC.Merged.r64296e203404D01.aligned.sorted.MinRQ998.OldToNewReads.csv.gz",
# ]

# pacbio_snps_positions_files = [
#     "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.positions.snps.csv.gz",
#     "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.positions.snps.csv.gz",
#     "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples/ADAR1.Merged.r64296e203404D01.aligned.sorted.MinRQ998.positions.snps.csv.gz",
#     "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples/IQEC.Merged.r64296e203404D01.aligned.sorted.MinRQ998.positions.snps.csv.gz",
# ]
# pacbio_snps_reads_files = [
#     "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.reads.snps.csv.gz",
#     "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.reads.snps.csv.gz",
#     "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples/ADAR1.Merged.r64296e203404D01.aligned.sorted.MinRQ998.reads.snps.csv.gz",
#     "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples/IQEC.Merged.r64296e203404D01.aligned.sorted.MinRQ998.reads.snps.csv.gz",
# ]

# %%
pacbio_samples = [
    "GRIA2", "PCLO", 
    # "GRIA2", "ADAR1", "IQEC1"
]

pacbio_new_positions_files = [
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/CompleteNoisePositions/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.positions.csv.gz",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/CompleteNoisePositions/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.positions.csv.gz",
    # "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/CompleteNoisePositions/comp141693_c0_seq1.merged.MinRQ998.positions.csv.gz",
    # "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/CompleteNoisePositions/comp134400_c0_seq1_extended.merged.MinRQ998.positions.csv.gz",
    # "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/CompleteNoisePositions/comp141565_c6_seq3.merged.MinRQ998.positions.csv.gz"
]
pacbio_new_old_to_new_reads_files = [
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/CompleteNoisePositions/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.OldToNewReads.csv.gz",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/CompleteNoisePositions/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.OldToNewReads.csv.gz",
    # "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/CompleteNoisePositions/comp141693_c0_seq1.merged.MinRQ998.OldToNewReads.csv.gz",
    # "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/CompleteNoisePositions/comp134400_c0_seq1_extended.merged.MinRQ998.OldToNewReads.csv.gz",
    # "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/CompleteNoisePositions/comp141565_c6_seq3.merged.MinRQ998.OldToNewReads.csv.gz"
]

pacbio_old_reads_files = [
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.reads.csv.gz",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.reads.csv.gz",
    # "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads/comp141693_c0_seq1.merged.MinRQ998.reads.csv.gz",
    # "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads/comp134400_c0_seq1_extended.merged.MinRQ998.reads.csv.gz",
    # "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads/comp141565_c6_seq3.merged.MinRQ998.reads.csv.gz"
]
pacbio_old_unique_reads_files = [
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_reads.csv.gz",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_reads.csv.gz",
    # "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads/comp141693_c0_seq1.merged.MinRQ998.unique_reads.csv.gz",
    # "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads/comp134400_c0_seq1_extended.merged.MinRQ998.unique_reads.csv.gz",
    # "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads/comp141565_c6_seq3.merged.MinRQ998.unique_reads.csv.gz"
]
pacbio_old_old_to_new_reads_files = [
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.OldToNewReads.csv.gz",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.OldToNewReads.csv.gz",
    # "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads/comp141693_c0_seq1.merged.MinRQ998.OldToNewReads.csv.gz",
    # "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads/comp134400_c0_seq1_extended.merged.MinRQ998.OldToNewReads.csv.gz",
    # "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads/comp141565_c6_seq3.merged.MinRQ998.OldToNewReads.csv.gz"
]

# pacbio_snps_positions_files = [
#     "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.positions.snps.csv.gz",
#     "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.positions.snps.csv.gz",
#     # "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads/comp141693_c0_seq1.merged.MinRQ998.positions.snps.csv.gz",
#     # "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads/comp134400_c0_seq1_extended.merged.MinRQ998.positions.snps.csv.gz",
#     # "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads/comp141565_c6_seq3.merged.MinRQ998.positions.snps.csv.gz"
# ]
pacbio_noise_positions_files = [
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.positions.noise.csv.gz",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.positions.noise.csv.gz",
]
pacbio_snps_reads_files = [
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.reads.snps.csv.gz",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.reads.snps.csv.gz",
    # "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads/comp141693_c0_seq1.merged.MinRQ998.reads.snps.csv.gz",
    # "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads/comp134400_c0_seq1_extended.merged.MinRQ998.reads.snps.csv.gz",
    # "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads/comp141565_c6_seq3.merged.MinRQ998.reads.snps.csv.gz"
]

# %%
illumina_samples = [
    "RUSC2",
    "TRIM2",
    "CA2D3",  # here's the problem, line 43
    "ABL",
    "DGLA",
    "K0513",
    "KCNAS",
    "ACHA4",
    "ANR17",
    "TWK7",
    "SCN1",
    "CACB2",
    "RIMS2",
    "PCLO",
    "DOP1",
    "IQEC1",
    "CSKI1",
    "MTUS2",
    "ROBO2",
]
illumina_chroms = [
    "comp141881_c0_seq3",
    "comp141044_c0_seq2",
    "comp140439_c0_seq1",
    "comp126362_c0_seq1",
    "comp141517_c0_seq1",
    "comp141840_c0_seq2",
    "comp141640_c0_seq1",
    "comp140987_c3_seq1",
    "comp140910_c2_seq1",
    "comp136058_c0_seq1",
    "comp141378_c0_seq7",
    "comp141158_c1_seq2",
    "comp140712_c0_seq3",
    "comp141882_c0_seq14",
    "comp141880_c1_seq3",
    "comp141565_c6_seq3",
    "comp141684_c0_seq1",
    "comp141532_c3_seq11",
    "comp141574_c0_seq3",
]

illumina_new_positions_files = [
    f"/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/CompleteNoisePositions/reads.sorted.aligned.filtered.{chrom}.positions.csv.gz"
    for chrom in illumina_chroms
]
illumina_new_old_to_new_reads_files = [
    f"/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/CompleteNoisePositions/reads.sorted.aligned.filtered.{chrom}.OldToNewReads.csv.gz"
    for chrom in illumina_chroms
]

illumina_old_reads_files = [
    f"/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.{chrom}.reads.csv"
    for chrom in illumina_chroms
]
illumina_old_unique_reads_files = [
    f"/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.{chrom}.unique_reads.csv"
    for chrom in illumina_chroms
]
illumina_old_old_to_new_reads_files = [
    f"/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.{chrom}.OldToNewReads.csv"
    for chrom in illumina_chroms
]

# illumina_snps_positions_files = [
#     f"/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.{chrom}.positions.snps.csv.gz"
#     for chrom in illumina_chroms
# ]
illumina_noise_positions_files = [
    f"/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.{chrom}.positions.noise.csv.gz"
    for chrom in illumina_chroms
]
illumina_snps_reads_files = [
    f"/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.{chrom}.reads.snps.csv.gz"
    for chrom in illumina_chroms
]

# %%
first_illumina_sample_index = 0
last_illumina_sample_index = len(illumina_samples)

# using only a few samples for now to test the code and make sure it runs in a reasonable time.
# first_illumina_sample_index = 13
# last_illumina_sample_index = 16
# first_illumina_sample_index = 0
# last_illumina_sample_index = 2

num_illumina_samples = last_illumina_sample_index - first_illumina_sample_index

illumina_samples[first_illumina_sample_index:last_illumina_sample_index]

# %%
# platforms = ["PacBio"] * len(pacbio_samples) + ["Illumina"] * len(illumina_samples)
# parities = ["SE"] * len(pacbio_samples) + ["PE"] * len(illumina_samples)
# samples = pacbio_samples + illumina_samples
# # new positions files made for complete noise estimates, 
# # and their corresponding old_to_new_reads_files
# new_positions_files = pacbio_new_positions_files + illumina_new_positions_files
# new_old_to_new_reads_files = pacbio_new_old_to_new_reads_files + illumina_new_old_to_new_reads_files
# # original reads - used for estimating editing
# old_reads_files = pacbio_old_reads_files + illumina_old_reads_files
# old_unique_reads_files = pacbio_old_unique_reads_files + illumina_old_unique_reads_files
# old_old_to_new_reads_files = pacbio_old_old_to_new_reads_files + illumina_old_old_to_new_reads_files



platforms = ["PacBio"] * len(pacbio_samples) + ["Illumina"] * num_illumina_samples
parities = ["SE"] * len(pacbio_samples) + ["PE"] * num_illumina_samples
samples = pacbio_samples + illumina_samples[first_illumina_sample_index:last_illumina_sample_index]
# new positions files made for complete noise estimates, 
# and their corresponding old_to_new_reads_files
new_positions_files = pacbio_new_positions_files + illumina_new_positions_files[first_illumina_sample_index:last_illumina_sample_index]
new_old_to_new_reads_files = pacbio_new_old_to_new_reads_files + illumina_new_old_to_new_reads_files[first_illumina_sample_index:last_illumina_sample_index]
# original reads - used for estimating editing
old_reads_files = pacbio_old_reads_files + illumina_old_reads_files[first_illumina_sample_index:last_illumina_sample_index]
old_unique_reads_files = pacbio_old_unique_reads_files + illumina_old_unique_reads_files[first_illumina_sample_index:last_illumina_sample_index]
old_old_to_new_reads_files = pacbio_old_old_to_new_reads_files + illumina_old_old_to_new_reads_files[first_illumina_sample_index:last_illumina_sample_index]

# snps_positions_files = pacbio_snps_positions_files + illumina_snps_positions_files[first_illumina_sample_index:last_illumina_sample_index]
noise_positions_files = pacbio_noise_positions_files + illumina_noise_positions_files[first_illumina_sample_index:last_illumina_sample_index]
snps_reads_files = pacbio_snps_reads_files + illumina_snps_reads_files[first_illumina_sample_index:last_illumina_sample_index]

# %%
samples

# %%
# parities

# %%
len(samples)


# %% [markdown]
# # Functions

# %% [markdown]
# ## Positions

# %%
def find_alt_base(ref_base, a_count, t_count, c_count, g_count, seed):

    bases = list("ATCG")
    base_counts = [a_count, t_count, c_count, g_count]

    alt_base_counts = {
        base: count for base, count in zip(bases, base_counts) if base != ref_base
    }

    max_alt_base_count = max(alt_base_counts.values())

    max_alt_bases = [
        base for base, count in alt_base_counts.items() if count == max_alt_base_count
    ]

    # randomly choose one of the equally most frequent alt bases
    return np.random.default_rng(seed).choice(max_alt_bases)


# %%
def calc_noise(ref_base_count, alt_base_count):
    try:
        noise = alt_base_count / (alt_base_count + ref_base_count)
    except ZeroDivisionError:
        noise = 0  # if there are no mapped alt bases
    return noise


# %%
# def make_snps_positions(
#     positions_file: Path | str,
#     snps_positions_file: Path | str,
#     top_x_noisy_positions: int,
#     assurance_factor: float,
#     seed: int,
#     edited_col: str,
#     snp_noise_level: float,
#     final_noise_col: str | None,
#     new_old_to_new_reads_file: Path | str | None = None,
#     old_old_to_new_reads_file: Path | str | None = None,
# ):
#     noise_positions_df = pd.read_csv(positions_file, sep="\t")
#     noise_positions_df = noise_positions_df.loc[
#         (noise_positions_df["CDS"])
#         & (~noise_positions_df[edited_col])
#         & (noise_positions_df["TotalCoverage"].gt(0))
#     ]
#     if final_noise_col is not None:
#         noise_positions_df = noise_positions_df.loc[
#             noise_positions_df[final_noise_col]
#         ]
#     noise_positions_df = noise_positions_df.drop(
#         columns=[
#             "CDS",
#             edited_col,
#             "InProbRegion",
#             "EditingFrequency",
#             "KnownEditing",
#         ]
#     )
    
#     # finalize empty df with expected cols, save it and return
#     if noise_positions_df.empty:
#         snps_positions_df = noise_positions_df.assign(
#             AltBase=None,
#             AboveEditingThreshold=None,
#         )
#         snps_positions_df.to_csv(snps_positions_file, sep="\t", index=False)
#         return snps_positions_df
    
#     noise_positions_df["AltBase"] = noise_positions_df.apply(
#         lambda x: find_alt_base(
#             x["RefBase"],
#             x["A"],
#             x["T"],
#             x["C"],
#             x["G"],
#             seed,
#         ),
#         axis=1,
#     )
#     noise_positions_df = noise_positions_df.loc[
#         ~(
#             (noise_positions_df["RefBase"].eq("A"))
#             & (noise_positions_df["AltBase"].eq("G"))
#         )
#     ]
#     noise_positions_df["Noise"] = noise_positions_df.apply(
#         lambda x: calc_noise(x[x["RefBase"]], x[x["AltBase"]]), axis=1
#     )
    
#     # determine editing threshold to annotate which noise positions have high noise and are thus annotated as "strong" SNPs
#     if final_noise_col is not None:
#         noise_levels = (
#             noise_positions_df.loc[
#                 (noise_positions_df["Noise"] < snp_noise_level)
#                 & (noise_positions_df[final_noise_col]),
#                 "Noise",
#             ]
#             .sort_values(ascending=False)[:top_x_noisy_positions]
#             .tolist()
#         )
#     else:
#         noise_levels = (
#             noise_positions_df.loc[
#                 (noise_positions_df["Noise"] < snp_noise_level),
#                 "Noise",
#             ]
#             .sort_values(ascending=False)[:top_x_noisy_positions]
#             .tolist()
#         )
#     # if there are less noisy positions than `top_x_noisy_positions`, add zeros accordingly
#     noise_levels = pd.Series(
#         noise_levels + [0 for _ in range(top_x_noisy_positions - len(noise_levels))]
#     )
#     noise_threshold = noise_levels.mean()
#     if pd.isna(noise_threshold):
#         noise_threshold = 0
#     # finalize the editing threshold
#     noise_threshold *= assurance_factor
    
#     snps_positions_df = noise_positions_df.loc[
#         noise_positions_df["Noise"].ge(snp_noise_level)
#     ].copy()

#     # now annotate the strong SNPs
#     snps_positions_df["AboveEditingThreshold"] = snps_positions_df["Noise"].ge(noise_threshold)

#     if (old_old_to_new_reads_file is not None) and (new_old_to_new_reads_file is not None):

#         # map old and new compressed reads names
#         new_old_to_new_reads_df = pd.read_table(new_old_to_new_reads_file).rename(
#             columns={"NewRead": "NewCompressedReadName"}
#         )
#         old_old_to_new_reads_df = pd.read_table(old_old_to_new_reads_file).rename(
#             columns={"NewRead": "OldCompressedReadName"}
#         )
#         new_to_old_shorted_reads_df = new_old_to_new_reads_df.merge(old_old_to_new_reads_df, on="OldRead", how="outer")
#         assert not new_to_old_shorted_reads_df.isna().any().any(), f"{positions_file = }, {snps_positions_file = }"
#         new_to_old_shorted_reads_dict = dict(
#             zip(
#                 new_to_old_shorted_reads_df['NewCompressedReadName'], 
#                 new_to_old_shorted_reads_df['OldCompressedReadName']
#             )
#         )
        
#         # use the mapping to replace the new compressed reads names with the old ones (the originals)
#         snps_positions_df["SplitNewCompressedReads"] = snps_positions_df["Reads"].str.split(",")
#         snps_positions_df["SplitOldCompressedReads"] = snps_positions_df["SplitNewCompressedReads"].apply(
#             lambda x: [new_to_old_shorted_reads_dict[y] for y in x],
#         )
#         snps_positions_df["Reads"] = snps_positions_df["SplitOldCompressedReads"].apply(lambda x: ",".join(x))
#         del snps_positions_df["SplitNewCompressedReads"]
#         del snps_positions_df["SplitOldCompressedReads"]
    
#     snps_positions_df.to_csv(snps_positions_file, sep="\t", index=False)

#     return snps_positions_df

# %%
def make_noise_positions(
    positions_file: Path | str,
    noise_positions_file: Path | str,
    top_x_noisy_positions: int,
    assurance_factor: float,
    seed: int,
    edited_col: str,
    snp_noise_level: float,
    final_noise_col: str | None,
    new_old_to_new_reads_file: Path | str | None = None,
    old_old_to_new_reads_file: Path | str | None = None,
):
    noise_positions_df = pd.read_csv(positions_file, sep="\t")
    noise_positions_df = noise_positions_df.loc[
        (noise_positions_df["CDS"])
        & (~noise_positions_df[edited_col])
        & (noise_positions_df["TotalCoverage"].gt(0))
    ]
    if final_noise_col is not None:
        noise_positions_df = noise_positions_df.loc[
            noise_positions_df[final_noise_col]
        ]
    noise_positions_df = noise_positions_df.drop(
        columns=[
            "CDS",
            edited_col,
            "InProbRegion",
            "EditingFrequency",
            "KnownEditing",
        ]
    )
    
    # finalize empty df with expected cols, save it and return
    if noise_positions_df.empty:
        noise_positions_df = noise_positions_df.assign(
            AltBase=None,
            EditingThreshold=None,
            AboveEditingThreshold=None,
            SNP=None
        )
        noise_positions_df.to_csv(
            noise_positions_file, sep="\t", index=False
        )
        return noise_positions_df
    
    noise_positions_df["AltBase"] = noise_positions_df.apply(
        lambda x: find_alt_base(
            x["RefBase"],
            x["A"],
            x["T"],
            x["C"],
            x["G"],
            seed,
        ),
        axis=1,
    )
    # noise_positions_df = noise_positions_df.loc[
    #     ~(
    #         (noise_positions_df["RefBase"].eq("A"))
    #         & (noise_positions_df["AltBase"].eq("G"))
    #     )
    # ]
    noise_positions_df["Noise"] = noise_positions_df.apply(
        lambda x: calc_noise(x[x["RefBase"]], x[x["AltBase"]]), axis=1
    )

    # determine editing threshold to annotate which noise positions have high noise and are thus possible false-positives
    if final_noise_col is not None:
        noise_levels = (
            noise_positions_df.loc[
                (noise_positions_df["Noise"] < snp_noise_level)
                & (noise_positions_df["RefBase"].ne("A"))
                & (noise_positions_df[final_noise_col]),
                "Noise",
            ]
            .sort_values(ascending=False)[:top_x_noisy_positions]
            .tolist()
        )
    else:
        noise_levels = (
            noise_positions_df.loc[
                (noise_positions_df["Noise"] < snp_noise_level)
                & (noise_positions_df["RefBase"].ne("A")),
                "Noise",
            ]
            .sort_values(ascending=False)[:top_x_noisy_positions]
            .tolist()
        )
    # if there are less noisy positions than `top_x_noisy_positions`, add zeros accordingly
    noise_levels = pd.Series(
        noise_levels + [0 for _ in range(top_x_noisy_positions - len(noise_levels))]
    )
    editing_threshold = noise_levels.mean()
    if pd.isna(editing_threshold):
        editing_threshold = 0
    # finalize the editing threshold
    editing_threshold *= assurance_factor

    noise_positions_df["EditingThreshold"] = editing_threshold
    noise_positions_df["AboveEditingThreshold"] = noise_positions_df["Noise"].ge(editing_threshold)
    
    # annotate which positions are considered SNPs based on the snp_noise_level
    noise_positions_df["SNP"] = (
        (noise_positions_df["Noise"].ge(snp_noise_level))
        & (noise_positions_df["RefBase"].ne("A"))
        & (noise_positions_df["AltBase"].ne("G"))
    )

    # since the reads names in the noise positions file are the new compressed reads names, 
    # we need to map them back to the old compressed reads names (the originals) using the old_to_new_reads files
    if (old_old_to_new_reads_file is not None) and (new_old_to_new_reads_file is not None):

        # map old and new compressed reads names
        new_old_to_new_reads_df = pd.read_table(new_old_to_new_reads_file).rename(
            columns={"NewRead": "NewCompressedReadName"}
        )
        old_old_to_new_reads_df = pd.read_table(old_old_to_new_reads_file).rename(
            columns={"NewRead": "OldCompressedReadName"}
        )
        new_to_old_shorted_reads_df = new_old_to_new_reads_df.merge(old_old_to_new_reads_df, on="OldRead", how="outer")
        assert not new_to_old_shorted_reads_df.isna().any().any(), f"{positions_file = }, {noise_positions_file = }"
        new_to_old_shorted_reads_dict = dict(
            zip(
                new_to_old_shorted_reads_df['NewCompressedReadName'], 
                new_to_old_shorted_reads_df['OldCompressedReadName']
            )
        )
        
        # use the mapping to replace the new compressed reads names with the old ones (the originals)
        noise_positions_df["SplitNewCompressedReads"] = noise_positions_df["Reads"].str.split(",")
        noise_positions_df["SplitOldCompressedReads"] = noise_positions_df["SplitNewCompressedReads"].apply(
            lambda x: [new_to_old_shorted_reads_dict[y] for y in x],
        )
        noise_positions_df["Reads"] = noise_positions_df["SplitOldCompressedReads"].apply(lambda x: ",".join(x))
        del noise_positions_df["SplitNewCompressedReads"]
        del noise_positions_df["SplitOldCompressedReads"]
        
    noise_positions_df.to_csv(noise_positions_file, sep="\t", index=False)

    return noise_positions_df



# %% [markdown]
# ## Reads

# %%
def add_noise_se_reads_in_positions(noise_positions_df, sample):
    unique_reads = set(chain.from_iterable(noise_positions_df["Reads"].str.split(",")))
    positions_per_read = {read: [] for read in unique_reads}

    # alt_base = "G" if strand == "+" else "C"
    for x, position_row in enumerate(noise_positions_df.itertuples()):
        mapped_bases = position_row.MappedBases
        mapped_reads = position_row.Reads.split(",")
        alt_base = position_row.AltBase
        # deal with reads that were mapped to this pos
        for mapped_base, mapped_read in zip(mapped_bases, mapped_reads):
            if mapped_base == alt_base:
                pos_edited_in_read = 1  # mismatch
            elif mapped_base == ".":
                pos_edited_in_read = 0  # match
            else:
                pos_edited_in_read = (
                    -1
                )  # we ignore mismatches that aren't RefBase2AltBase by marking them as NaN
            positions_per_read[mapped_read].append(pos_edited_in_read)
        # deal with reads that weren't mapped to this pos
        unmapped_reads_in_pos = unique_reads - set(mapped_reads)
        for mapped_read in unmapped_reads_in_pos:
            pos_edited_in_read = -1
            positions_per_read[mapped_read].append(pos_edited_in_read)
        # check all reads got values for pos
        mapped_pos_dist = {len(bases) for bases in positions_per_read.values()}
        if len(mapped_pos_dist) != 1:
            if sample is not None:
                raise Exception(
                    f"Problem at line {x} in sample {sample}: not all reads are mapped."
                )
            else:
                raise Exception(f"Problem at line {x}: not all reads are mapped.")

    return positions_per_read


# %%
def add_noise_pe_reads_in_positions(noise_positions_df, sample):
    unique_reads = set(chain.from_iterable(noise_positions_df["Reads"].str.split(",")))
    positions_per_read = {read: [] for read in unique_reads}

    overlapping_bases = 0

    # alt_base = "G" if strand == "+" else "C"
    for x, position_row in enumerate(noise_positions_df.itertuples()):
        mapped_bases = position_row.MappedBases
        mapped_reads = position_row.Reads.split(",")
        phred_scores = position_row.Phred
        alt_base = position_row.AltBase

        last_phred_of_read_in_position = {}

        # deal with reads that were mapped to this pos
        for phred_score, mapped_base, mapped_read in zip(
            phred_scores, mapped_bases, mapped_reads
        ):
            if mapped_base == alt_base:
                pos_edited_in_read = 1  # mismatch
            elif mapped_base == ".":
                pos_edited_in_read = 0  # match
            else:
                pos_edited_in_read = (
                    -1
                )  # we ignore mismatches that aren't RefBase2AltBase by marking them as NaN

            # deal with overlapping paired-ended reads mapped to this position
            if mapped_read not in last_phred_of_read_in_position:
                positions_per_read[mapped_read].append(pos_edited_in_read)
                last_phred_of_read_in_position[mapped_read] = phred_score
            else:
                # for a pair of PE reads, we only keep the one with the highest phred score
                if ord(phred_score) > ord(last_phred_of_read_in_position[mapped_read]):
                    positions_per_read[mapped_read][
                        -1
                    ] = pos_edited_in_read  # overwrite the previous editing status of the read in this position
                    last_phred_of_read_in_position[mapped_read] = phred_score
                overlapping_bases += 1

        # deal with reads that weren't mapped to this pos
        unmapped_reads_in_pos = unique_reads - set(mapped_reads)
        for mapped_read in unmapped_reads_in_pos:
            pos_edited_in_read = -1
            positions_per_read[mapped_read].append(pos_edited_in_read)

        # check all reads got values for pos
        mapped_pos_dist = {len(bases) for bases in positions_per_read.values()}
        if len(mapped_pos_dist) != 1:
            if sample is not None:
                raise Exception(
                    f"Problem at line {x} in sample {sample}: not all reads are mapped."
                )
            else:
                raise Exception(f"Problem at line {x}: not all reads are mapped.")

    # if sample is not None:
    #     ic(sample, overlapping_bases)

    return positions_per_read


# %%
def add_noise_multisample_se_reads_in_positions(noise_positions_df, sample):
    unique_reads = set(chain.from_iterable(noise_positions_df["Reads"].str.split(",")))
    positions_per_read = {read: [] for read in unique_reads}
    sample_per_read = {}

    # alt_base = "G" if strand == "+" else "C"
    for x, position_row in enumerate(noise_positions_df.itertuples()):
        mapped_bases = position_row.MappedBases
        mapped_reads = position_row.Reads.split(",")
        mapped_samples = position_row.Samples.split(",")
        alt_base = position_row.AltBase
        # deal with reads that were mapped to this pos
        for mapped_base, mapped_read, mapped_sample in zip(mapped_bases, mapped_reads, mapped_samples):
            
            sample_per_read[mapped_read] = mapped_sample
            
            if mapped_base == alt_base:
                pos_edited_in_read = 1  # mismatch
            elif mapped_base == ".":
                pos_edited_in_read = 0  # match
            else:
                pos_edited_in_read = (
                    -1
                )  # we ignore mismatches that aren't RefBase2AltBase by marking them as NaN
            positions_per_read[mapped_read].append(pos_edited_in_read)
        # deal with reads that weren't mapped to this pos
        unmapped_reads_in_pos = unique_reads - set(mapped_reads)
        for mapped_read in unmapped_reads_in_pos:
            pos_edited_in_read = -1
            positions_per_read[mapped_read].append(pos_edited_in_read)
        # check all reads got values for pos
        mapped_pos_dist = {len(bases) for bases in positions_per_read.values()}
        if len(mapped_pos_dist) != 1:
            raise Exception(
                    f"Problem at line {x} in sample {sample}: not all reads are mapped."
                )

    return positions_per_read


# %%
def make_snps_reads(
    snps_positions_df,
    snps_reads_file,
    platform, 
    sample, 
    parity, 
    multisample, 
):
    # If there are no SNP positions, there are no reads to encode.
    if snps_positions_df is None or snps_positions_df.empty:
        reads_snps_df = pd.DataFrame(columns=["Platform", "Sample", "Read"])
        reads_snps_df.to_csv(snps_reads_file, index=False, sep="\t")
        return reads_snps_df
    
    # defensive: ensure Reads/Samples exist and are strings
    if "Reads" not in snps_positions_df.columns:
        raise KeyError("snps_positions_df is missing required column 'Reads'")
    if multisample and "Samples" not in snps_positions_df.columns:
        raise KeyError("snps_positions_df is missing required column 'Samples'")
    
    if parity == "SE":
        if multisample:
            positions_per_read = add_noise_multisample_se_reads_in_positions(
                snps_positions_df, sample=sample
            )
        else:
            positions_per_read = add_noise_se_reads_in_positions(
                snps_positions_df, sample=sample
            )
        
    elif parity == "PE":
        if multisample:
            raise NotImplementedError("Multisample PE reads in positions is not implemented yet.")
        else:
            positions_per_read = add_noise_pe_reads_in_positions(
                snps_positions_df, sample=sample
            )
    else:
        raise ValueError("parity must be either 'SE' or 'PE'")

    reads_snps_df = pd.DataFrame(positions_per_read)
    reads_snps_df = reads_snps_df.T.reset_index().rename(
        {"index": "Read"}, axis="columns"
    )
    reads_snps_df = reads_snps_df.rename(
        columns={
            old_col: pos
            for old_col, pos in zip(
                reads_snps_df.columns[1:], snps_positions_df["Position"]
            )
        }
    )
            
    reads_snps_df.insert(0, "Platform", platform)
    reads_snps_df.insert(1, "Sample", sample)
    
    reads_snps_df.to_csv(snps_reads_file, index=False, sep="\t")

    return reads_snps_df


# %% [markdown]
# ## Positions & reads

# %%
# def make_positions_and_reads_snps_dfs(
#     positions_file,
#     snps_positions_file,
#     top_x_noisy_positions,
#     assurance_factor,
#     seed,
#     edited_col,
#     snp_noise_level,
#     final_noise_col,
#     snps_reads_file,
#     platform,
#     sample,
#     parity,
#     multisample,
#     new_old_to_new_reads_file=None,
#     old_old_to_new_reads_file=None,
# ):
#     try:
#         snps_positions_df = make_snps_positions(
#             positions_file, 
#             snps_positions_file, 
#             top_x_noisy_positions,
#             assurance_factor,
#             seed, 
#             edited_col, 
#             snp_noise_level, 
#             final_noise_col,
#             new_old_to_new_reads_file=new_old_to_new_reads_file,
#             old_old_to_new_reads_file=old_old_to_new_reads_file
#         )
#         snps_reads_df = make_snps_reads(
#             snps_positions_df, snps_reads_file, platform, sample, parity, multisample
#         )
#         return snps_positions_df, snps_reads_df
#     except Exception as e:
#         print(f"Error occurred while processing {positions_file}: {e}")
#         raise

# %%
def make_noise_positions_and_snps_reads_dfs(
    positions_file,
    noise_positions_file,
    top_x_noisy_positions,
    assurance_factor,
    seed,
    edited_col,
    snp_noise_level,
    final_noise_col,
    snps_reads_file,
    platform,
    sample,
    parity,
    multisample,
    new_old_to_new_reads_file=None,
    old_old_to_new_reads_file=None,
):
    try:
        noise_positions_df = make_noise_positions(
            positions_file, 
            noise_positions_file, 
            top_x_noisy_positions,
            assurance_factor,
            seed, 
            edited_col, 
            snp_noise_level, 
            final_noise_col,
            new_old_to_new_reads_file=new_old_to_new_reads_file,
            old_old_to_new_reads_file=old_old_to_new_reads_file
        )
        snps_positions_df = noise_positions_df.loc[
            noise_positions_df["SNP"]
        ].copy()
        snps_reads_df = make_snps_reads(
            snps_positions_df, snps_reads_file, platform, sample, parity, multisample
        )
        return noise_positions_df, snps_reads_df
    except Exception as e:
        print(f"Error occurred while processing {positions_file}: {e}")
        raise

# %% [markdown]
# # Finding SNPs

# %%
# # ?make_positions_and_reads_snps_dfs

# %%
# samples

# %%
# # i = 2
# i = 0
# # i = 3
# positions_file = new_positions_files[i]
# noise_positions_file = noise_positions_files[i]
# # top_x_noisy_positions
# # assurance_factor
# # seed
# # edited_col
# # snp_noise_level
# # final_noise_col
# snps_reads_file = snps_reads_files[i]
# platform = platforms[i]
# sample = samples[i]
# parity = parities[i]
# # multisample
# new_old_to_new_reads_file = new_old_to_new_reads_files[i]
# old_old_to_new_reads_file = old_old_to_new_reads_files[i]

# %%
# positions_file

# %%
# # testing new version of make_snps_positions,
# # where all noise positions are kept - not only snps,
# # and each position is annotated with whether it is above the editing threshold or not.

# noise_positions_df = pd.read_csv(positions_file, sep="\t")
# noise_positions_df = noise_positions_df.loc[
#     (noise_positions_df["CDS"])
#     & (~noise_positions_df[edited_col])
#     & (noise_positions_df["TotalCoverage"].gt(0))
# ]
# if final_noise_col is not None:
#     noise_positions_df = noise_positions_df.loc[
#         noise_positions_df[final_noise_col]
#     ]
# noise_positions_df = noise_positions_df.drop(
#     columns=[
#         "CDS",
#         edited_col,
#         "InProbRegion",
#         "EditingFrequency",
#         "KnownEditing",
#     ]
# )
# noise_positions_df

# %%
# # finalize empty df with expected cols, save it and return
# if noise_positions_df.empty:
#     noise_positions_df = noise_positions_df.assign(
#         AltBase=None,
#         EditingThreshold=None,
#         AboveEditingThreshold=None,
#         SNP=None
#     )
#     noise_positions_df.to_csv(noise_positions_file, sep="\t", index=False)
#     # return noise_positions_df  # TODO uncomment this return statement when running the full pipeline, for now we want to see the empty df

# %%
# noise_positions_df["AltBase"] = noise_positions_df.apply(
#     lambda x: find_alt_base(
#         x["RefBase"],
#         x["A"],
#         x["T"],
#         x["C"],
#         x["G"],
#         seed,
#     ),
#     axis=1,
# )
# noise_positions_df = noise_positions_df.loc[
#     ~(
#         (noise_positions_df["RefBase"].eq("A"))
#         & (noise_positions_df["AltBase"].eq("G"))
#     )
# ]
# noise_positions_df["Noise"] = noise_positions_df.apply(
#     lambda x: calc_noise(x[x["RefBase"]], x[x["AltBase"]]), axis=1
# )

# # determine editing threshold to annotate which noise positions have high noise and are thus annotated as "strong" SNPs
# if final_noise_col is not None:
#     noise_levels = (
#         noise_positions_df.loc[
#             (noise_positions_df["Noise"] < snp_noise_level)
#             & (noise_positions_df[final_noise_col]),
#             "Noise",
#         ]
#         .sort_values(ascending=False)[:top_x_noisy_positions]
#         .tolist()
#     )
# else:
#     noise_levels = (
#         noise_positions_df.loc[
#             (noise_positions_df["Noise"] < snp_noise_level),
#             "Noise",
#         ]
#         .sort_values(ascending=False)[:top_x_noisy_positions]
#         .tolist()
#     )
# # if there are less noisy positions than `top_x_noisy_positions`, add zeros accordingly
# noise_levels = pd.Series(
#     noise_levels + [0 for _ in range(top_x_noisy_positions - len(noise_levels))]
# )
# noise_threshold = noise_levels.mean()
# if pd.isna(noise_threshold):
#     noise_threshold = 0
# # finalize the editing threshold
# noise_threshold *= assurance_factor

# ic(noise_threshold)

# noise_positions_df["EditingThreshold"] = noise_threshold
# noise_positions_df["AboveEditingThreshold"] = noise_positions_df["Noise"].ge(noise_threshold)
# noise_positions_df["SNP"] = noise_positions_df["Noise"].ge(snp_noise_level)

# noise_positions_df

# %%
# if (old_old_to_new_reads_file is not None) and (new_old_to_new_reads_file is not None):

#     # map old and new compressed reads names
#     new_old_to_new_reads_df = pd.read_table(new_old_to_new_reads_file).rename(
#         columns={"NewRead": "NewCompressedReadName"}
#     )
#     old_old_to_new_reads_df = pd.read_table(old_old_to_new_reads_file).rename(
#         columns={"NewRead": "OldCompressedReadName"}
#     )
#     new_to_old_shorted_reads_df = new_old_to_new_reads_df.merge(old_old_to_new_reads_df, on="OldRead", how="outer")
#     assert not new_to_old_shorted_reads_df.isna().any().any(), f"{positions_file = }, {noise_positions_file = }"
#     new_to_old_shorted_reads_dict = dict(
#         zip(
#             new_to_old_shorted_reads_df['NewCompressedReadName'], 
#             new_to_old_shorted_reads_df['OldCompressedReadName']
#         )
#     )
    
#     # use the mapping to replace the new compressed reads names with the old ones (the originals)
#     noise_positions_df["SplitNewCompressedReads"] = noise_positions_df["Reads"].str.split(",")
#     noise_positions_df["SplitOldCompressedReads"] = noise_positions_df["SplitNewCompressedReads"].apply(
#         lambda x: [new_to_old_shorted_reads_dict[y] for y in x],
#     )
#     noise_positions_df["Reads"] = noise_positions_df["SplitOldCompressedReads"].apply(lambda x: ",".join(x))
#     del noise_positions_df["SplitNewCompressedReads"]
#     del noise_positions_df["SplitOldCompressedReads"]
    
# noise_positions_df.to_csv(noise_positions_file, sep="\t", index=False)
# # return noise_positions_df  # TODO uncomment this return statement when running the full pipeline, for now we want to see the empty df

# noise_positions_df


# %%
# noise_positions_df = make_noise_positions(
#     positions_file,
#     noise_positions_file,
#     top_x_noisy_positions,
#     assurance_factor,
#     seed,
#     edited_col,
#     snp_noise_level,
#     final_noise_col,
#     new_old_to_new_reads_file,
#     old_old_to_new_reads_file,
# )
# noise_positions_df

# %%
# # old threshold obtained by using samtools 2

# threshold_candidates_df = noise_positions_df.loc[
#     noise_positions_df["Noise"].lt(snp_noise_level)
#     & noise_positions_df["RefBase"].ne("A")
# ]

# top3_df = threshold_candidates_df.nlargest(
#     top_x_noisy_positions, "Noise"
# )

# display(
#     top3_df[
#         ["Position", "RefBase", "AltBase", "Noise", "TotalCoverage"]
#     ]
# )

# top3_df["Noise"].sum() / top_x_noisy_positions * assurance_factor

# %%
# # new threshold obtained by using samtools 1.18

# threshold_candidates_df = noise_positions_df.loc[
#     noise_positions_df["Noise"].lt(snp_noise_level)
#     & noise_positions_df["RefBase"].ne("A")
# ]

# top3_df = threshold_candidates_df.nlargest(
#     top_x_noisy_positions, "Noise"
# )

# display(
#     top3_df[
#         ["Position", "RefBase", "AltBase", "Noise", "TotalCoverage"]
#     ]
# )

# top3_df["Noise"].sum() / top_x_noisy_positions * assurance_factor

# %%
# snps_positions_df = noise_positions_df.loc[
#     noise_positions_df["SNP"]
# ].copy()

# snps_positions_df

# %%
# snps_reads_df = make_snps_reads(
#         snps_positions_df, snps_reads_file, platform, sample, parity, multisample
#     )
# snps_reads_df

# %%
# make_positions_and_reads_snps_dfs(
#     positions_file,
#     snps_positions_file,
#     top_x_noisy_positions,
#     assurance_factor,
#     seed,
#     edited_col,
#     snp_noise_level,
#     final_noise_col,
#     snps_reads_file,
#     platform,
#     sample,
#     parity,
#     multisample,
#     new_old_to_new_reads_file,
#     old_old_to_new_reads_file
# )

# %%
# make_noise_positions_and_snps_reads_dfs(
#     positions_file,
#     noise_positions_file,
#     top_x_noisy_positions,
#     assurance_factor,
#     seed,
#     edited_col,
#     snp_noise_level,
#     final_noise_col,
#     snps_reads_file,
#     platform,
#     sample,
#     parity,
#     multisample,
#     new_old_to_new_reads_file,
#     old_old_to_new_reads_file,
# )

# %%
# with Pool(processes=len(samples)) as pool:
#     positions_and_reads_snps_dfs = pool.starmap(
#         func=make_positions_and_reads_snps_dfs,
#         iterable=[
#             (
#                 positions_file,
#                 snps_positions_file,
#                 top_x_noisy_positions,
#                 assurance_factor,
#                 seed,
#                 edited_col,
#                 snp_noise_level,
#                 final_noise_col,
#                 snps_reads_file,
#                 platform,
#                 sample,
#                 parity,
#                 multisample,
#                 new_old_to_new_reads_file,
#                 old_old_to_new_reads_file
#             )
#             for (
#                 positions_file, 
#                 snps_positions_file, 
#                 snps_reads_file, 
#                 platform, 
#                 sample, 
#                 parity, 
#                 new_old_to_new_reads_file, 
#                 old_old_to_new_reads_file
#             ) in zip(
#                 new_positions_files,
#                 snps_positions_files,
#                 snps_reads_files,
#                 platforms,
#                 samples,
#                 parities,
#                 new_old_to_new_reads_files,
#                 old_old_to_new_reads_files
#             )
            
#         ],
#     )

# positions_snps_dfs = [dfs[0] for dfs in positions_and_reads_snps_dfs]
# reads_snps_dfs = [dfs[1] for dfs in positions_and_reads_snps_dfs]

# %%
with Pool(processes=len(samples)) as pool:
    noise_positions_and_snps_reads_dfs = pool.starmap(
        func=make_noise_positions_and_snps_reads_dfs,
        iterable=[
            (
                positions_file,
                noise_positions_file,
                top_x_noisy_positions,
                assurance_factor,
                seed,
                edited_col,
                snp_noise_level,
                final_noise_col,
                snps_reads_file,
                platform,
                sample,
                parity,
                multisample,
                new_old_to_new_reads_file,
                old_old_to_new_reads_file
            )
            for (
                positions_file, 
                noise_positions_file, 
                snps_reads_file, 
                platform, 
                sample, 
                parity, 
                new_old_to_new_reads_file, 
                old_old_to_new_reads_file
            ) in zip(
                new_positions_files,
                noise_positions_files,
                snps_reads_files,
                platforms,
                samples,
                parities,
                new_old_to_new_reads_files,
                old_old_to_new_reads_files
            )
        #     ) in zip(
        #         new_positions_files[:2],
        #         noise_positions_files[:2],
        #         snps_reads_files[:2],
        #         platforms[:2],
        #         samples[:2],
        #         parities[:2],
        #         new_old_to_new_reads_files[:2],
        #         old_old_to_new_reads_files[:2]
        #     )  
        ],
    )

noise_positions_dfs = [dfs[0] for dfs in noise_positions_and_snps_reads_dfs]
snps_reads_dfs = [dfs[1] for dfs in noise_positions_and_snps_reads_dfs]

# %%
noise_positions_df = noise_positions_dfs[0]
noise_positions_df

# %%
for sample, noise_positions_df in zip(samples, noise_positions_dfs):
    ic(
        sample,
        noise_positions_df.loc[
            (noise_positions_df["RefBase"].eq("A"))
            & (noise_positions_df["AltBase"].eq("G")),
            & (noise_positions_df["Noise"].ge(noise_positions_df["EditingThreshold"]))
        ].shape[0]
    )

# %%
noise_positions_df.loc[
    (noise_positions_df["RefBase"].eq("A"))
    & (noise_positions_df["AltBase"].eq("G")),
    # & (noise_positions_df["Noise"].ge(noise_positions_df["EditingThreshold"]))
    "Noise"
].describe()

# %%
reads_noise_first_pos_loc = 3
if new_noise_positions_in_use:
    reads_noise_first_pos_loc += 1
    
reads_noise_first_pos_loc

# %%
noise_per_position_dfs = []
for reads_snps_df, sample, platform in zip(reads_snps_dfs, samples, platforms):
    noise_per_position_df = (
        reads_snps_df.iloc[:, reads_noise_first_pos_loc:]
        # .apply(lambda x: x.eq(1).sum() / x.ne(-1).sum())
        # .sort_values(ascending=False)
        # .reset_index()
        # .rename(columns={"index": "Position", 0: "Noise"})
        .apply(lambda x: 100 * x.eq(1).sum() / x.ne(-1).sum())
        .sort_values(ascending=False)
        .reset_index()
        .rename(columns={"index": "Position", 0: "%Noise"})
    )
    noise_per_position_df.insert(0, "Platform", platform)
    noise_per_position_df.insert(1, "Sample", sample)
    noise_per_position_dfs.append(noise_per_position_df)
concat_noise_per_position_df = pd.concat(noise_per_position_dfs, ignore_index=True).sort_values(
    by=["Platform", "Sample", "%Noise"], ascending=[True, True, False], ignore_index=True
)
# concat_noise_per_position_df["%Noise"] = concat_noise_per_position_df["Noise"] * 100
concat_noise_per_position_df

# %%
concat_noise_per_position_df.groupby("Sample")["%Noise"].describe().round(2)

# %%
concat_noise_per_position_df.loc[
    (concat_noise_per_position_df["%Noise"].ge(snp_noise_level)) 
    # & (concat_noise_per_position_df["%Noise"].le(max_noise_prct_threshold))
].round(2)

# %%
concat_noise_per_position_df.loc[
    (concat_noise_per_position_df["%Noise"].ge(snp_noise_level)) 
    # & (concat_noise_per_position_df["%Noise"].le(max_noise_prct_threshold))
].groupby(["Platform", "Sample"]).size().reset_index(name="NumOfSNPs")

# %%
# thresholds_step = 5
# min_min_threshold = 30
# max_min_threshold = 50
# min_max_threshold = 50
# max_max_threshold = 70
# snp_noise_levels = range(min_min_threshold, max_min_threshold+thresholds_step, thresholds_step)
# max_noise_prct_thresholds = range(min_max_threshold, max_max_threshold+thresholds_step, thresholds_step)
# # min_thresholds_diff = 20

# num_of_snps_per_sample_dfs = []
# for min_thr, max_thr in product(snp_noise_levels, max_noise_prct_thresholds):
#     # ic(min_thr, max_thr)
#     # if max_thr - min_thr < min_thresholds_diff:
#     #     continue
#     if max_thr <= min_thr:
#         continue
#     # ic(min_thr, max_thr)
#     num_of_snps_per_sample_df = concat_noise_per_position_df.loc[
#         (concat_noise_per_position_df["%Noise"].ge(min_thr)) &
#         (concat_noise_per_position_df["%Noise"].le(max_thr))
#     ].groupby(["Platform", "Sample"]).size().reset_index(name="NumOfSNPs")
#     num_of_snps_per_sample_df.insert(2, "MinNoisePrct", min_thr)
#     num_of_snps_per_sample_df.insert(3, "MaxNoisePrct", max_thr)
#     num_of_snps_per_sample_dfs.append(num_of_snps_per_sample_df)
# concat_num_of_snps_per_sample_df = pd.concat(
#     num_of_snps_per_sample_dfs, ignore_index=True
# ).sort_values(
#     by=["Platform", "Sample", "MinNoisePrct", "MaxNoisePrct"],
#     ascending=[True, True, True, True],
#     ignore_index=True,
# )
# concat_num_of_snps_per_sample_df


# %%
thresholds_step = 10
min_min_threshold = 10
max_min_threshold = 100
# min_max_threshold = 50
# max_max_threshold = 70
snp_noise_levels = range(min_min_threshold, max_min_threshold+thresholds_step, thresholds_step)
# max_noise_prct_thresholds = range(min_max_threshold, max_max_threshold+thresholds_step, thresholds_step)
# min_thresholds_diff = 20

num_of_snps_per_sample_dfs = []
for min_thr in snp_noise_levels:
    # print(f"{min_thr = }")
    num_of_snps_per_sample_df = concat_noise_per_position_df.loc[
        concat_noise_per_position_df["%Noise"].ge(min_thr)
    ].groupby(["Platform", "Sample"]).size().reset_index(name="NumOfSNPs")
    num_of_snps_per_sample_df.insert(2, "MinNoisePrct", min_thr)
    # num_of_snps_per_sample_df.insert(3, "MaxNoisePrct", max_thr)
    num_of_snps_per_sample_dfs.append(num_of_snps_per_sample_df)
concat_num_of_snps_per_sample_df = pd.concat(
    num_of_snps_per_sample_dfs, ignore_index=True
).sort_values(
    by=["Platform", "Sample", "MinNoisePrct"],
    ascending=[True, True, True],
    ignore_index=True,
)
concat_num_of_snps_per_sample_df

# %%
concat_num_of_snps_per_sample_df["MinNoisePrct"].nunique()

# %%
int(9/2)

# %%
fig = px.bar(
    concat_num_of_snps_per_sample_df,
    x="Sample",
    y="NumOfSNPs",
    facet_col="MinNoisePrct",
    facet_col_wrap=int(concat_num_of_snps_per_sample_df["MinNoisePrct"].nunique() / 2),
    facet_row_spacing=0.1,
    # facet_row="MaxNoisePrct",
    color="Platform",
    title="Number of SNPs per Sample at Different Noise Percentage Thresholds",
)
fig.update_yaxes(dtick=1)
fig.update_layout(height=500, width=1000, template="simple_white")
fig.show()

# %% [markdown]
# # Using SNPs to separate editing by allele

# %% [markdown]
# Next, only for genes in which we found SNPs, we'll try to separate editing by allele.

# %%
# snp_noise_level = 40
snp_noise_level = 10
# max_noise_prct_threshold = 60

snps_df = concat_noise_per_position_df.loc[
    (concat_noise_per_position_df["%Noise"].ge(snp_noise_level)) 
    # & (concat_noise_per_position_df["%Noise"].le(max_noise_prct_threshold))
# ].drop(columns="Noise").round(2)
].round(2)

snps_df

# %%
snps_platform_and_sample = {(sample, platform) for sample, platform in snps_df.loc[:, ["Platform", "Sample"]].values}
snps_platform_and_sample

# %%
snps_platforms = []
snps_parities = []
snps_samples = []
snps_old_reads_files = []
snps_old_unique_reads_files = []
snps_old_old_to_new_reads_files = []
snps_reads_snps_dfs = []
snps_positions_lists = []

for platform, pairity, sample, old_reads_file, old_unique_reads_file, old_old_to_new_reads_file, reads_snps_df in zip(
    platforms, parities, samples, old_reads_files, old_unique_reads_files, old_old_to_new_reads_files, reads_snps_dfs
):
    if (platform, sample) not in snps_platform_and_sample:
        continue
    print(platform, sample)
    snps_platforms.append(platform)
    snps_parities.append(pairity)
    snps_samples.append(sample)
    snps_old_reads_files.append(old_reads_file)
    snps_old_unique_reads_files.append(old_unique_reads_file)
    snps_old_old_to_new_reads_files.append(old_old_to_new_reads_file)
    snps_reads_snps_dfs.append(reads_snps_df)
    snps_positions = snps_df.loc[
        (snps_df["Platform"].eq(platform))
        & (snps_df["Sample"].eq(sample)),
        "Position"
    ].values
    snps_positions_lists.append(snps_positions)

# snps_platforms, snps_samples
ic(
    snps_platforms, snps_parities, snps_samples, snps_old_reads_files, 
    snps_old_unique_reads_files, snps_old_old_to_new_reads_files,
    snps_positions_lists
)

# %%
# snps_df

# %%
# for chrom, sample in zip(illumina_chroms, illumina_samples):
#     print(f"{chrom} - {sample}")

# %%
snps_reads_snps_dfs[0]

# %%
# snps_reads_snps_dfs[1]

# %%
snps_old_reads_dfs = [
    pd.read_table(old_reads_file)
    for old_reads_file in snps_old_reads_files
]

snps_old_reads_dfs[0]

# %%
# snps_old_reads_dfs[1]

# %%

# %%
snps_df

# %%
snps_positions = snps_positions_lists[0]
snps_reads_noise_df = snps_reads_snps_dfs[0]
snps_old_reads_df = snps_old_reads_dfs[0]
snps_old_old_to_new_reads_file = snps_old_old_to_new_reads_files[0]

# %%
snps_positions

# %%
# snps_old_old_to_new_reads_file

# %%
# snps_old_reads_df

# %%
# old_reads_first_pos_loc = 6

# %%
# incorporate OldRead column into snps_old_reads_df

snps_old_old_to_new_reads_df = pd.read_table(snps_old_old_to_new_reads_file)
# snps_old_old_to_new_reads_df

snps_old_reads_df = snps_old_reads_df.merge(
    snps_old_old_to_new_reads_df,
    left_on="Read",
    right_on="NewRead",
    how="left"
).drop(
    # retain "OldRead" and "Read" (which is the same as "NewRead")
    columns="NewRead"
)
snps_old_reads_df.insert(
    snps_old_reads_df.columns.get_loc("Read") + 1,
    "OldRead2",
    snps_old_reads_df["OldRead"]
)
snps_old_reads_df = snps_old_reads_df.drop(
    columns="OldRead"
).rename(
    columns={"OldRead2": "OldRead"}
)

snps_old_reads_df = snps_old_reads_df.sort_values(
    "OldRead",
    ignore_index=True
)

current_old_reads_first_pos_loc = old_reads_first_pos_loc + 1
ic(current_old_reads_first_pos_loc)

snps_old_reads_df

# %%
# snps_reads_noise_df

# %%
# keep only SNPs from snps_reads_noise_df
snps_reads_noise_df = pd.concat(
    [
        # reads_noise_df.iloc[:, :reads_noise_first_pos_loc],
        snps_reads_noise_df.loc[:, "OldRead"],
        snps_reads_noise_df.loc[:, snps_positions]
    ],
    axis=1
)

# keep SNPs info only for reads for which we have editing info in snps_old_reads_df
snps_reads_noise_df = snps_reads_noise_df.merge(
    snps_old_reads_df.loc[:, ["OldRead"]],
    on="OldRead",
    how="right"
)

snps_reads_noise_df = snps_reads_noise_df.sort_values(
    "OldRead",
    ignore_index=True
)

snps_reads_noise_df

# %%
assert snps_old_reads_df.shape[0] == snps_reads_noise_df.shape[0], "Number of reads mismatch between old reads df and reads noise df"

for snp_position in snps_positions:
    print(snp_position)
    # we can directly insert this col from snps_reads_noise_df into snps_old_reads_df 
    # because both dataframes are sorted by OldRead
    snps_old_reads_df.insert(
        current_old_reads_first_pos_loc,
        f"SNP_{snp_position}",
        snps_reads_noise_df[snp_position]
    )
    current_old_reads_first_pos_loc += 1
    ic(current_old_reads_first_pos_loc);
    
snps_old_reads_df

# %%

# %%
# Identify the editing site columns (starting from the current offset onwards)
editing_cols = snps_old_reads_df.columns[current_old_reads_first_pos_loc:]
editing_cols

# %%
snp_cols = [f"SNP_{pos}" for pos in snps_positions]
snp_cols

# %%
# Filter to keep only reads where the SNP positions have a valid base (0 or 1, not -1)
# Using .all(axis=1) ensures we have a complete haplotype for the selected SNPs
valid_reads_mask = snps_old_reads_df[snp_cols].ne(-1).all(axis=1)
temp_df = snps_old_reads_df.loc[valid_reads_mask].copy()

# Replace -1 (unmapped/ignored) with NaN in the editing columns
# so means are computed on covered reads only
temp_df[editing_cols] = temp_df[editing_cols].replace(-1, np.nan)

# # 3) define a single haplotype label (useful for counting / plotting / downstream merges)
# # e.g. "SNP_123=0|SNP_456=1|SNP_789=0"
# # haplotype_col = "|".join(snp_cols)
# # temp_df[haplotype_col] = temp_df[snp_cols].astype(str).agg("|".join, axis=1)
# temp_df["Haplotype"] = temp_df[snp_cols].astype(str).agg("|".join, axis=1)
temp_df.insert(
    current_old_reads_first_pos_loc,
    "Haplotype",
    temp_df[snp_cols].astype(str).agg("|".join, axis=1)
)
current_old_reads_first_pos_loc += 1
ic(current_old_reads_first_pos_loc)

temp_df

# %%
temp_df.groupby("Haplotype")["287"].value_counts().reset_index()

# %% [markdown]
# ## Basic plots

# %%
# temp_df.groupby("Haplotype")[editing_cols].agg(["mean", "std"])
editing_long_df = (
    temp_df.groupby("Haplotype")[editing_cols]
    .agg(["mean", "std"])
    .stack(level=0, future_stack=True)  # silence FutureWarning (pandas>=2.1)
    .reset_index()
    .rename(columns={"level_1": "EditingSite"})
)

# # optional: nicer dtypes / sorting
editing_long_df["EditingSite"] = editing_long_df["EditingSite"].astype(int)
editing_long_df = editing_long_df.sort_values(["EditingSite", "Haplotype"], ignore_index=True)
editing_long_df["EditingSite"] = editing_long_df["EditingSite"].astype(str)

editing_long_df

# %%
fig = px.scatter(
    editing_long_df,
    x="EditingSite",
    y="mean",
    color="Haplotype",
    # error_y="std",
    # facet_row="Haplotype",
    opacity=0.8,
)
fig.update_layout(
    height=500,
    width=1000,
    template="simple_white",
    # title=f"Editing Status Distribution per Editing Site<br>for Reads with Complete SNP Haplotypes"
)
fig.show()

# %%
# Group by the SNP alleles (haplotypes) and calculate the editing frequency for each site
# (rows=editing sites, cols=haplotypes)
# editing_profiles = temp_df.groupby("Haplotype")[editing_cols].mean().T
editing_profiles = temp_df.groupby("Haplotype")[editing_cols].mean().mul(100).T
editing_profiles.dropna(how="all", inplace=True)
# editing_profiles["MeanOfMeans"] = editing_profiles.mean(axis=1)
editing_profiles = editing_profiles.merge(
    temp_df.loc[:, editing_cols].mean().mul(100).T.rename("GlobalMeanEditing"),
    left_index=True,
    right_index=True
)
editing_profiles["MeanEditingDiff(1-0)"] = editing_profiles.apply(
    lambda x: x.iloc[1] - x.iloc[0], axis=1
)

editing_profiles

# %%
editing_profiles["MeanEditingDiff(1-0)"].sort_values(ascending=False).round(2)

# %%
fig = px.histogram(
    editing_profiles,
    x="MeanEditingDiff(1-0)",
    labels={
        "MeanEditingDiff(1-0)": "Mean % editing difference per site (Haplotype 1 - Haplotype 0)"
    }
)
fig.update_xaxes(dtick=1)
fig.update_layout(
    height=500,
    width=1000,
    template="simple_white",
)
fig.show()

# %%
# fig = px.histogram(
#     editing_profiles,
#     x="MeanOfMeans",
#     y="MeanEditingDiff(1-0)",
#     labels={
#         "MeanEditingDiff(1-0)": "Mean % editing difference per site (Haplotype 1 - Haplotype 0)",
#         "MeanOfMeans": "Mean of means % editing per site"
#     },
#     histfunc="avg",
# )
# # fig.update_xaxes(dtick=1)
# fig.update_layout(
#     height=500,
#     width=1000,
#     template="simple_white",
# )
# fig.show()

# %%
# fig = px.density_heatmap(
#     editing_profiles,
#     y="GlobalMeanEditing",
#     x="MeanEditingDiff(1-0)",
#     labels={
#         "MeanEditingDiff(1-0)": "Mean % editing difference per site (Haplotype 1 - Haplotype 0)",
#         "GlobalMeanEditing": "Global mean % editing per site"
#     },
#     # histfunc="avg",
# )
# # fig.update_xaxes(dtick=1)
# fig.update_layout(
#     height=500,
#     width=500,
#     template="simple_white",
# )
# fig.show()

# %%
fig = px.scatter(
    editing_profiles,
    x="GlobalMeanEditing",
    y="MeanEditingDiff(1-0)",
    labels={
        "MeanEditingDiff(1-0)": "Mean % editing difference per site<br>(Haplotype 1 - Haplotype 0)",
        "GlobalMeanEditing": "Global mean % editing per site"
    },
    marginal_x="histogram", marginal_y="histogram",
    trendline="ols", trendline_color_override="black"
)


results = px.get_trendline_results(fig)
print(results.iloc[0].px_fit_results.summary())

# fig.update_xaxes(dtick=1)
fig.update_layout(
    height=700,
    width=700,
    template="simple_white",
)
fig.show()

# %%

# %%
temp_df.groupby("Haplotype")["287"].value_counts().reset_index()


# %% [markdown]
# ## Global enrichment with chi-square permutations

# %%
def _per_site_counts(df, site_col, haplotype_col, haplotypes):
    """
    Returns a 2x2 (or rxc where r=len(haplotypes) and c=2 (unedited/edited)) counts:
      [[h1_edited0, h1_edited1],
       [h2_edited0, h2_edited1]]
    using df rows with non-NaN at site_col and haplotype_col in haplotypes (e.g. {h1,h2}).
    """
    sub = df.loc[df[haplotype_col].isin(haplotypes), [haplotype_col, site_col]].dropna()
    # keep only 0/1
    sub = sub.loc[sub[site_col].isin([0, 1])]
    if sub.empty:
        return None

    cros = pd.crosstab(sub[haplotype_col], sub[site_col])
    # ensure full shape
    cros = cros.reindex(index=haplotypes, columns=[0, 1], fill_value=0)
    return cros.values


# %%
def _fisher_two_sided(contingency_table):
    _, p = fisher_exact(contingency_table, alternative="two-sided")
    return p


# %%
def _chi_square_two_sided(contingency_table):
    res = chi2_contingency(contingency_table)
    return res.pvalue


# %%
def count_diff_sites(
    df, haplotype_col, haplotypes, editing_sites_cols, 
    alpha=0.05, min_reads_per_cell=5, 
    test=_chi_square_two_sided
):
    pvals = []
    tested_cols = []
    for site_col in editing_sites_cols:
        contingency_table = _per_site_counts(df, site_col, haplotype_col, haplotypes)
        if contingency_table is None:
            continue
        if any(contingency_table.flatten() < min_reads_per_cell):
            print(f"Skipping site {site_col} due to low total reads ({contingency_table.sum()})")
            continue
        # p = _fisher_two_sided(contingency_table)
        p = test(contingency_table)
        pvals.append(p)
        tested_cols.append(site_col)

    pvals = np.asarray(pvals)
    rejected, pvalue_corrected = fdrcorrection(
        pvals,
        alpha=alpha,
        # general correlated tests, see:
        # https://www.statsmodels.org/stable/generated/statsmodels.stats.multitest.fdrcorrection.html#statsmodels.stats.multitest.fdrcorrection.method
        method="n"
    )
    tests_df = pd.DataFrame(
        {
            "Site": tested_cols,
            "P": pvals,
            "P_FDR": pvalue_corrected,
            "Rejected_H0": rejected,
        }
    ).sort_values("Site")
    # num of significant sites after FDR correction
    n_sig = sum(rejected)
    return n_sig, len(tested_cols), tests_df

# n_sig, n_tested_cols, tests_df = count_diff_sites(temp_df, "Haplotype", editing_cols, "0", "1", alpha=0.05, min_reads_per_cell=5)


# %%
haplotype_col = "Haplotype"
site_col = "287"
haplotypes = ["0", "1"]
df = temp_df

# %%
_per_site_counts(temp_df, site_col, haplotype_col, haplotypes)

# %%
n_sig, n_tested_cols, tests_df = count_diff_sites(df, haplotype_col, haplotypes, editing_cols)

# %%
n_sig, n_tested_cols

# %%
tests_df

# %%

# %%
df = temp_df
haplotype_col = "Haplotype"
haplotypes = ["0", "1"]
n_perm = 200
alpha = 0.05
min_reads_per_cell = 5
seed = 1892
test = _chi_square_two_sided

# %%
rng = np.random.default_rng(seed)
# # copy relevant subset to allow shuffling haplotype labels w/o modifying original df
sub = df.loc[df[haplotype_col].isin(haplotypes)].copy()
# calculate observed number of significant sites in the original data
observed, m_tested, per_site = count_diff_sites(
    sub, haplotype_col, haplotypes, editing_cols, alpha=alpha, min_reads_per_cell=min_reads_per_cell, test=test
)
observed, m_tested

# %%
null_counts = []
hap_values = sub[haplotype_col].to_numpy()
for _ in range(n_perm):
    # Randomly shuffle (“permute”) the haplotype labels across reads
    perm = hap_values.copy()
    rng.shuffle(perm)
    sub[haplotype_col] = perm
    # Recalculate the number of significant sites for this permutation
    c, _, _ = count_diff_sites(
        sub, haplotype_col, haplotypes, editing_cols, alpha=alpha, min_reads_per_cell=min_reads_per_cell, test=test
    )
    null_counts.append(c)

# %%
null_counts

# %%
pd.Series(null_counts).value_counts()

# %%
pd.Series(null_counts).describe()

# %%
pd.Series(null_counts).ge(observed).sum()

# %%
# empirical p-value
p_emp = (1 + pd.Series(null_counts).ge(observed).sum()) / (1 + n_perm)
p_emp

# %%
np.round(p_emp, 5)


# %%

# %%

# %%
def permutation_null(df, hap_col, editing_cols, h1, h2, n_perm=200, alpha=0.05, min_total_reads=20, seed=1892):
    rng = np.random.default_rng(seed)
    # copy relevant subset to allow shuffling haplotype labels w/o modifying original df
    sub = df.loc[df[hap_col].isin([h1, h2])].copy()
    # calculate observed number of significant sites in the original data
    observed, m_tested, per_site = count_diff_sites(sub, hap_col, editing_cols, h1, h2, alpha=alpha, min_total_reads=min_total_reads)

    null_counts = []
    hap_values = sub[hap_col].to_numpy()
    for _ in range(n_perm):
        # Randomly shuffle (“permute”) the haplotype labels across reads
        perm = hap_values.copy()
        rng.shuffle(perm)
        sub[hap_col] = perm
        # Recalculate the number of significant sites for this permutation
        c, _, _ = count_diff_sites(sub, hap_col, editing_cols, h1, h2, alpha=alpha, min_total_reads=min_total_reads)
        null_counts.append(c)

    null_counts = np.asarray(null_counts)
    # empirical p-value
    p_emp = (1 + (null_counts >= observed).sum()) / (1 + n_perm)
    return observed, null_counts, p_emp, m_tested, per_site

# ---- run for your temp_df ----
hap_counts = temp_df["Haplotype"].value_counts()
top2 = hap_counts.index[:2].tolist()
h1, h2 = top2[0], top2[1]
print("Comparing:", h1, "vs", h2)
print("Reads per haplotype:\n", hap_counts.loc[[h1, h2]])

obs, null_counts, p_emp, m_tested, per_site = permutation_null(
    temp_df, "Haplotype", editing_cols, h1, h2,
    n_perm=200, alpha=0.05, min_total_reads=20, seed=1
)
print(f"Tested sites: {m_tested}")
print(f"Observed #diff sites (p<0.05): {obs}")
print(f"Null mean: {null_counts.mean():.2f}, null 95%: [{np.quantile(null_counts, 0.025)}, {np.quantile(null_counts, 0.975)}]")
print(f"Empirical p-value: {p_emp:.4g}")

per_site.head(10)

# %%

# %%

# %% [markdown]
# ## Logistic regression: EditingStatus ~ C(EditingSite) * C(Haplotype)

# %%
temp_df

# %%
reg_df = temp_df.loc[
    :,
    ["Read", "Haplotype"] + editing_cols.tolist()
].melt(
    id_vars=["Read", "Haplotype"],
    var_name="EditingSite",
    value_name="EditingStatus"
)

reg_df = reg_df.loc[
    reg_df["EditingStatus"].notna()
].reset_index(drop=True)

reg_df["EditingStatus"] = reg_df["EditingStatus"].astype(int)

reg_df

# %%
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats

# %%
# Null model — site effect only
# logitP(edited=1)=αEditingSite​
# Under this model, two reads at the same site are expected to have the
# same editing probability, regardless of haplotype.

m0 = smf.glm(
    formula="EditingStatus ~ C(EditingSite)",
    data=reg_df,
    family=sm.families.Binomial()
).fit()

print(m0.summary())

# %%
8.80e+06 / 8799382

# %%
# Alternative model — add a global haplotype shift
# logitP(edited=1)=αEditingSite​+βHaplotype​
# Interpretation:
# we still allow each site to have its own baseline rate
# but now haplotypes can differ by a constant multiplicative (log-odds) shift

m1 = smf.glm(
    formula="EditingStatus ~ C(EditingSite) + C(Haplotype)",
    data=reg_df,
    family=sm.families.Binomial()
).fit()

print(m1.summary())

# %%
m1.llf > m0.llf

# %%
np.abs(m1.llf - m0.llf)

# %%
#m0 vs m1 likelihood ratio test: Does haplotype affect editing globally?

LLR = 2 * (m1.llf - m0.llf)
df_diff = m1.df_model - m0.df_model
p_global_haplotype = stats.chi2.sf(LLR, df_diff)

print("Global haplotype effect LRT:")
print(f"  LLR = {LLR:.3f}, df = {df_diff}, p = {p_global_haplotype:.3e}")

# %%
# Model 2
# EditingStatus ~ C(EditingSite) * C(Haplotype)
# The * expands to
# EditingStatus ~ C(EditingSite) + C(Haplotype) + C(EditingSite):C(Haplotype)
# The interaction term γ allows:
# haplotype effects to differ per site
# direction & magnitude to vary
# This model asks: “Do haplotypes show distinct editing profiles, not just a uniform shift?”

m2 = smf.glm(
    formula="EditingStatus ~ C(EditingSite) * C(Haplotype)",
    data=reg_df,
    family=sm.families.Binomial()
).fit()

print(m2.summary())

# %%
m2.llf > m1.llf

# %%
np.abs(m2.llf - m1.llf)

# %%
LLR_int = 2 * (m2.llf - m1.llf)
df_int = m2.df_model - m1.df_model
p_int = stats.chi2.sf(LLR_int, df_int)

print("Haplotype × Site interaction LRT:")
print(f"  LLR = {LLR_int:.3f}, df = {df_int}, p = {p_int:.3e}")


# %%
m1.llf

# %%
m2.llf

# %%
m1.df_model

# %%
m2.df_model

# %%

# %% [markdown]
# ## Multinomial logistic regression: Haplotype  ~  EditingProfile

# %%
temp_df

# %%
# X_wide shape is (n_reads, n_sites)
X_wide = temp_df.set_index("Read").loc[:, editing_cols]
X_wide

# %%
# % missing data stats per site
X_wide.isna().sum().div(X_wide.shape[0]).mul(100).describe().round(2)

# %%
# Haplotype label per read
y = temp_df["Haplotype"]
y

# %% [markdown]
# ### Toy example

# %%
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report, confusion_matrix

# %%
# Train–test split
X_train, X_test, y_train, y_test = train_test_split(
    X_wide,
    y,
    test_size=0.2,
    stratify=y,        # keep haplotype proportions
    random_state=1892
)

# %%
# Pipeline: imputer + multinomial logistic regression
clf = Pipeline(
    [
        (
            "imputer", 
            SimpleImputer(
                # strategy="most_frequent",
                strategy="mean"  # default
            )
        ),
        (
            "logreg", 
            LogisticRegression(
                max_iter=1000,
                # multi_class="auto",   # supports >2 haplotypes
                # n_jobs=20             # use all cores
            )
        )
    ]
)

clf.fit(X_train, y_train)

# %%
# Classification performance
y_pred = clf.predict(X_test)

# %%
print(classification_report(y_test, y_pred))

# %%
conf_mat = confusion_matrix(y_test, y_pred)
conf_df = pd.DataFrame(
    conf_mat,
    index=[f"True_{hap}" for hap in clf.classes_],
    columns=[f"Pred_{hap}" for hap in clf.classes_]
)
print(conf_df)

# %%
y_proba = clf.predict_proba(X_test)  # shape: (n_test_reads, n_haplotypes)
haplotypes = clf.named_steps["logreg"].classes_

# Example: probabilities for first test read
print("Haplotypes:", haplotypes)
print("Probabilities:", y_proba[0])


# %%
y_proba

# %%
y_proba.mean(axis=0)

# %%
pd.Series(y_proba.mean(axis=1)).describe().r

# %% [markdown]
# ### Unified benchmarking pipeline w/ different classifiers

# %%
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.metrics import (
    classification_report,
    accuracy_score,
    f1_score
)

X_train, X_test, y_train, y_test = train_test_split(
    X_wide,
    y,
    test_size=0.2,
    stratify=y,
    random_state=42
)

# %%
# Define a model zoo

from sklearn.linear_model import LogisticRegression
from sklearn.svm import LinearSVC
from sklearn.naive_bayes import BernoulliNB
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.ensemble import HistGradientBoostingClassifier

# %%
# Pipelines (some need imputation, some don't)

models = {}

# strategy = "most_frequent"
strategy = "mean"

# 1) Logistic Regression (L2, class-balanced)
models["logreg_L2_balanced"] = Pipeline([
    ("imputer", SimpleImputer(strategy=strategy)),
    ("clf", LogisticRegression(
        max_iter=2000,
        class_weight="balanced"
    ))
])

# 2) Logistic Regression (L1 sparse model)
models["logreg_L1_sparse"] = Pipeline([
    ("imputer", SimpleImputer(strategy=strategy)),
    ("clf", LogisticRegression(
        penalty="elasticnet",
        l1_ratio=1.0,          # pure L1
        solver="saga",
        max_iter=4000,
        class_weight="balanced"
    ))
])

# 3) Linear SVM
models["linear_svm"] = Pipeline([
    ("imputer", SimpleImputer(strategy=strategy)),
    ("clf", LinearSVC(class_weight="balanced"))
])

# 4) Bernoulli Naive Bayes (natural fit for 0/1 features)
models["bernoulli_nb"] = Pipeline([
    ("imputer", SimpleImputer(strategy=strategy)),
    ("clf", BernoulliNB())
])

# 5) Random Forest
models["random_forest"] = Pipeline([
    ("imputer", SimpleImputer(strategy=strategy)),
    ("clf", RandomForestClassifier(
        n_estimators=300,
        max_depth=None,
        n_jobs=-1,
        class_weight="balanced_subsample"
    ))
])

# 6) Gradient Boosting
models["gradient_boosting"] = Pipeline([
    ("imputer", SimpleImputer(strategy=strategy)),
    ("clf", GradientBoostingClassifier())
])

# 7) HistGradientBoosting (handles NaN -> NO imputer)
models["hist_gb_nan_aware"] = HistGradientBoostingClassifier()


# %%
# Train + evaluate automatically

results = []

for name, model in models.items():
    print(f"\n🔹 Training {name} ...")
    model.fit(X_train, y_train)

    y_pred = model.predict(X_test)

    acc = accuracy_score(y_test, y_pred)
    f1m = f1_score(y_test, y_pred, average="macro")

    print(f"  accuracy = {acc:.3f}, macro-F1 = {f1m:.3f}")
    print(classification_report(y_test, y_pred))

    results.append((name, acc, f1m))

results_df = pd.DataFrame(results, columns=["model", "accuracy", "macro_f1"])
print(results_df.sort_values("macro_f1", ascending=False))

# %%
# Rank sites by permutation importance (HistGradientBoosting)

from sklearn.inspection import permutation_importance

# Use your trained model
hgb = models["hist_gb_nan_aware"]

r = permutation_importance(
    hgb,
    X_test,
    y_test,
    n_repeats=5,
    random_state=42,
    scoring="f1_macro",
    n_jobs=20
)

importance = pd.Series(r.importances_mean, index=X_wide.columns)\
              .sort_values(ascending=False)

# top_sites = importance.head(40).reset_index()
top_sites = importance.head(40).reset_index().rename(columns={"index": "EditingSite", 0: "Importance"})
top_sites

# %%
top_sites = importance.head(40).reset_index().rename(columns={"index": "EditingSite", 0: "Importance"})
top_sites

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %% [markdown]
# ## Reads --> unique reads --> distinct proteins ?

# %%
snps_old_unique_reads_dfs = [
    pd.read_table(old_unique_reads_file).rename(
        columns={"Transcript": "UniqueRead"}
    )
    for old_unique_reads_file in snps_old_unique_reads_files
]
snps_old_unique_reads_dfs[0]

# %%
snps_old_unique_reads_dfs[1]

# %% [markdown]
# # 12-mismatch analysis
#
# The two PacBio and 19 Illumina inputs follow cells `22e6a692` / `d99b4678`. Recovered counts and CDS positions come from `CompleteNoisePositions`, with its existing BQ30 and relative-coverage filtering. Original positions provide separate annotations and the historical threshold; no coverage filter is reapplied. Positions are 0-based and the current targets are on `+`.
#
# For a fresh session, run this section in order; no reads preparation or other squid notebook is needed. If the prepared `mm12_positions_df`, original identity table and successful threshold check are already in memory, run only `mm12-squid-decisions`, `mm12-squid-scheme-comparison`, then the revised Export cells. Do not rerun alternative selection for an export.
#
# The source decisions are carried into `12mm/v2`. Threshold reconstruction and its original-site identity check are unchanged.
#

# %% [markdown]
# ## Inputs

# %%
from pathlib import Path
import gzip
from IPython.display import display

import numpy as np
import pandas as pd

mm12_root = Path("/private6/projects/Combinatorics")
mm12_squid_root = mm12_root / "D.pealeii/MpileupAndTranscripts"
mm12_output_root = mm12_root / "Code/Notebooks/12mm/v2/positions"
mm12_seed = 1892  # Existing project seed, used only by this section's local generator.
mm12_snp_threshold = 0.10
mm12_top_noise_positions = 3
mm12_assurance_factor = 1.5
mm12_unit_cols = ["Dataset", "Platform", "Sample", "Gene", "Chrom"]
mm12_site_cols = mm12_unit_cols + ["Position"]
mm12_source_cols = [
    "Chrom", "Position", "RefBase", "TotalCoverage", "A", "T", "C", "G",
    "CDS", "Edited", "KnownEditing", "InProbRegion", "Noise", "EditingFrequency",
]

# %%
mm12_inputs_df = pd.DataFrame(
    [
        ("SquidLongReads", "GRIA2", "comp141693_c0_seq1", "GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998"),
        ("SquidLongReads", "PCLO", "comp141882_c0_seq14", "PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998"),
        ("SquidShortReads", "RUSC2", "comp141881_c0_seq3", "reads.sorted.aligned.filtered.comp141881_c0_seq3"),
        ("SquidShortReads", "TRIM2", "comp141044_c0_seq2", "reads.sorted.aligned.filtered.comp141044_c0_seq2"),
        ("SquidShortReads", "CA2D3", "comp140439_c0_seq1", "reads.sorted.aligned.filtered.comp140439_c0_seq1"),
        ("SquidShortReads", "ABL", "comp126362_c0_seq1", "reads.sorted.aligned.filtered.comp126362_c0_seq1"),
        ("SquidShortReads", "DGLA", "comp141517_c0_seq1", "reads.sorted.aligned.filtered.comp141517_c0_seq1"),
        ("SquidShortReads", "K0513", "comp141840_c0_seq2", "reads.sorted.aligned.filtered.comp141840_c0_seq2"),
        ("SquidShortReads", "KCNAS", "comp141640_c0_seq1", "reads.sorted.aligned.filtered.comp141640_c0_seq1"),
        ("SquidShortReads", "ACHA4", "comp140987_c3_seq1", "reads.sorted.aligned.filtered.comp140987_c3_seq1"),
        ("SquidShortReads", "ANR17", "comp140910_c2_seq1", "reads.sorted.aligned.filtered.comp140910_c2_seq1"),
        ("SquidShortReads", "TWK7", "comp136058_c0_seq1", "reads.sorted.aligned.filtered.comp136058_c0_seq1"),
        ("SquidShortReads", "SCN1", "comp141378_c0_seq7", "reads.sorted.aligned.filtered.comp141378_c0_seq7"),
        ("SquidShortReads", "CACB2", "comp141158_c1_seq2", "reads.sorted.aligned.filtered.comp141158_c1_seq2"),
        ("SquidShortReads", "RIMS2", "comp140712_c0_seq3", "reads.sorted.aligned.filtered.comp140712_c0_seq3"),
        ("SquidShortReads", "PCLO", "comp141882_c0_seq14", "reads.sorted.aligned.filtered.comp141882_c0_seq14"),
        ("SquidShortReads", "DOP1", "comp141880_c1_seq3", "reads.sorted.aligned.filtered.comp141880_c1_seq3"),
        ("SquidShortReads", "IQEC1", "comp141565_c6_seq3", "reads.sorted.aligned.filtered.comp141565_c6_seq3"),
        ("SquidShortReads", "CSKI1", "comp141684_c0_seq1", "reads.sorted.aligned.filtered.comp141684_c0_seq1"),
        ("SquidShortReads", "MTUS2", "comp141532_c3_seq11", "reads.sorted.aligned.filtered.comp141532_c3_seq11"),
        ("SquidShortReads", "ROBO2", "comp141574_c0_seq3", "reads.sorted.aligned.filtered.comp141574_c0_seq3"),
    ],
    columns=["Dataset", "Sample", "Chrom", "Stem"],
)
mm12_inputs_df["Platform"] = mm12_inputs_df["Dataset"].map(
    {"SquidLongReads": "PacBio", "SquidShortReads": "Illumina"}
)
mm12_inputs_df["Gene"] = mm12_inputs_df["Sample"]
mm12_inputs_df["Strand"] = "+"
mm12_inputs_df

# %%
mm12_inputs_df["SourceFile"] = [
    mm12_squid_root / "CompleteNoisePositions" / f"{stem}.positions.csv.gz"
    for stem in mm12_inputs_df["Stem"]
]
mm12_inputs_df["OriginalSourceFile"] = [
    mm12_squid_root / "RQ998.TopNoisyPositions3.BQ30" / f"{row.Stem}.positions.csv.gz"
    if row.Platform == "PacBio"
    else mm12_squid_root / "Illumina" / f"{row.Stem}.positions.csv"
    for row in mm12_inputs_df.itertuples(index=False)
]
mm12_inputs_df["OutputFile"] = [
    mm12_output_root / ("squid_long_reads" if row.Platform == "PacBio" else "squid_short_reads")
    / f"{row.Stem}.12mm.positions.csv.gz"
    for row in mm12_inputs_df.itertuples(index=False)
]
assert not mm12_inputs_df.duplicated(mm12_unit_cols).any(), "Repeated input unit"
assert mm12_inputs_df["SourceFile"].is_unique and mm12_inputs_df["OutputFile"].is_unique
mm12_missing_files = [
    file for col in ["SourceFile", "OriginalSourceFile"]
    for file in mm12_inputs_df[col] if not file.is_file()
]
if mm12_missing_files:
    raise FileNotFoundError(f"Missing explicit 12mm input files: {mm12_missing_files}")
mm12_existing_outputs = [file for file in mm12_inputs_df["OutputFile"] if file.exists()]
if mm12_existing_outputs:
    raise FileExistsError(f"12mm positions already exist; refusing to overwrite: {mm12_existing_outputs}")

# %%
# Load source columns only; never load reads, Phred strings or mapped-base strings.
mm12_original_dfs = []
mm12_recovered_dfs = []
for mm12_input in mm12_inputs_df.itertuples(index=False):
    for mm12_source_file, mm12_destination in [
        (mm12_input.OriginalSourceFile, mm12_original_dfs),
        (mm12_input.SourceFile, mm12_recovered_dfs),
    ]:
        mm12_source_df = pd.read_csv(mm12_source_file, sep="\t", usecols=mm12_source_cols)
        assert not mm12_source_df.empty, f"Empty input: {mm12_source_file}"
        assert not mm12_source_df[["Chrom", "Position"]].isna().any().any()
        assert not mm12_source_df.duplicated(["Chrom", "Position"]).any(), str(mm12_source_file)
        assert mm12_source_df["Chrom"].eq(mm12_input.Chrom).all(), str(mm12_source_file)
        assert mm12_source_df["RefBase"].isin(list("ATCG")).all(), str(mm12_source_file)
        assert mm12_source_df[["CDS", "Edited", "KnownEditing", "InProbRegion"]].isin([True, False]).all().all()
        mm12_counts = mm12_source_df[["A", "T", "C", "G", "TotalCoverage"]].to_numpy()
        assert np.isfinite(mm12_counts).all() and (mm12_counts >= 0).all()
        assert (mm12_counts == np.floor(mm12_counts)).all(), str(mm12_source_file)
        assert mm12_source_df["Position"].ge(0).all()
        assert mm12_source_df["Position"].eq(mm12_source_df["Position"].astype(int)).all()
        for mm12_col in ["Dataset", "Platform", "Sample", "Gene", "Strand"]:
            mm12_source_df[mm12_col] = getattr(mm12_input, mm12_col)
        mm12_destination.append(mm12_source_df.sort_values(["Chrom", "Position"]).reset_index(drop=True))

# %% [markdown]
# ## Original thresholds and site identity
#
# The directed-sequencing algorithm in `Code/Pileup/positions.py`, `pileup_to_positions` (also retained in `positions.old.py`), first applies its coverage filter and removes non-principal-base sites with `Noise >= 0.10`. The original positions files already contain that population. Its threshold is the mean of the first three descending **original `Noise` values**, multiplied by 1.5 once. This precedes any downstream CDS selection. Principal-base noise is NaN; pandas skips NaN in the mean. With fewer than three finite values, only those values contribute: no zero padding. A missing threshold is an error here.
#
# We do not use this notebook's later `b457f41f` threshold: that calculation restricts CDS/non-edited rows and pads with zeros, so equivalence is not assumed. The directed call uses de novo detection: `A>G` on `+` (`T>C` on `-`), frequency `alt / (ref + alt)`, strictly greater than the threshold. Its zero-denominator editing frequency is zero. There is no known-site, significance, problematic-region or gene-SNP exclusion in this original squid call. The two current input cohorts are all `+`.
#
# The acceptance check below reconstructs the original rule on **every original row**, using the original counts and no original editing flag in the reconstruction mask. It compares complete site identities with `Edited` and also checks the CDS subset used in the squid notebooks. Both directions of difference must be empty before dependent export. This check has been written, not run.

# %%
mm12_threshold_rows = []
mm12_original_identity_dfs = []
for mm12_original_df, mm12_input in zip(mm12_original_dfs, mm12_inputs_df.itertuples(index=False)):
    mm12_principal_ref = "A" if mm12_input.Strand == "+" else "T"
    mm12_principal_alt = "G" if mm12_input.Strand == "+" else "C"
    mm12_is_principal = mm12_original_df["RefBase"].eq(mm12_principal_ref)
    assert mm12_original_df.loc[mm12_is_principal, "Noise"].isna().all()
    assert mm12_original_df.loc[~mm12_is_principal, "Noise"].lt(0.10).all(), "Original noise filter differs"
    mm12_threshold = mm12_original_df["Noise"].sort_values(ascending=False).iloc[:mm12_top_noise_positions].mean()
    mm12_threshold *= mm12_assurance_factor
    if not np.isfinite(mm12_threshold):
        raise ValueError(f"Original editing threshold is missing for {mm12_input.Dataset}/{mm12_input.Sample}")
    mm12_threshold_rows.append({
        **{col: getattr(mm12_input, col) for col in mm12_unit_cols},
        "OldEditingThreshold": mm12_threshold,
        "ThresholdSource": f"{mm12_input.OriginalSourceFile}; original Noise top-3 mean x1.5; no CDS filter or zero padding",
    })
    mm12_old_denominator = mm12_original_df[mm12_principal_ref] + mm12_original_df[mm12_principal_alt]
    mm12_old_frequency = mm12_original_df[mm12_principal_alt].div(mm12_old_denominator)
    mm12_old_frequency = mm12_old_frequency.where(mm12_old_denominator.ne(0), 0.0).where(mm12_is_principal)
    assert np.allclose(mm12_old_frequency, mm12_original_df["EditingFrequency"], equal_nan=True), "Original editing-frequency definition differs"
    mm12_reconstructed = mm12_is_principal & mm12_old_frequency.gt(mm12_threshold)
    mm12_identity_df = mm12_original_df[mm12_site_cols + ["CDS"]].copy()
    mm12_identity_df["OriginalSite"] = mm12_original_df["Edited"]
    mm12_identity_df["ReconstructedSite"] = mm12_reconstructed
    mm12_identity_df["MissingOriginalSite"] = mm12_identity_df["OriginalSite"] & ~mm12_reconstructed
    mm12_identity_df["UnexpectedReconstructedSite"] = mm12_reconstructed & ~mm12_identity_df["OriginalSite"]
    mm12_original_identity_dfs.append(mm12_identity_df)
mm12_thresholds_df = pd.DataFrame(mm12_threshold_rows)
mm12_original_identity_df = pd.concat(mm12_original_identity_dfs, ignore_index=True)
assert not mm12_thresholds_df.duplicated(mm12_unit_cols).any()
assert not mm12_original_identity_df.duplicated(mm12_site_cols).any()

# %%
mm12_identity_summary_dfs = []
for mm12_scope, mm12_scope_mask in [
    ("All original positions", pd.Series(True, index=mm12_original_identity_df.index)),
    ("Original CDS positions", mm12_original_identity_df["CDS"]),
]:
    mm12_scope_df = mm12_original_identity_df.loc[mm12_scope_mask]
    mm12_summary_df = mm12_scope_df.groupby(mm12_unit_cols, as_index=False)[
        ["OriginalSite", "ReconstructedSite", "MissingOriginalSite", "UnexpectedReconstructedSite"]
    ].sum()
    mm12_summary_df["Scope"] = mm12_scope
    mm12_identity_summary_dfs.append(mm12_summary_df)
mm12_identity_summary_df = pd.concat(mm12_identity_summary_dfs, ignore_index=True)
mm12_identity_differences_df = mm12_original_identity_df.loc[
    mm12_original_identity_df["MissingOriginalSite"] | mm12_original_identity_df["UnexpectedReconstructedSite"]
].copy()
display(mm12_identity_summary_df)
if not mm12_identity_differences_df.empty:
    display(mm12_identity_differences_df.head(30))
    raise ValueError(
        "Original squid site identity failed; first 30 differences shown, all retained in "
        "mm12_identity_differences_df. Stop before exporting threshold-dependent positions."
    )
# Mark validation only after both full original and CDS site identities have passed.
mm12_thresholds_df["OriginalSchemeValidated"] = True

# %% [markdown]
# ## Annotations
#
# The merge retains every recovered row, using a 1:1 **left** join on `Chrom, Position` within the explicit input. Original-only positions are counted separately and their identities retained in `mm12_original_only_sites_df`; they are not restored or assigned zero coverage. Each source keeps its existing coverage population. Recovered base counts are never replaced by original counts. Historical flags and frequencies are prefixed `Original`; unmatched recovered rows retain missing historical annotations plus `OriginalAnnotationPresent=False`. The recovered detection fields are prefixed `Recovered`. Current CDS/known-site/problematic-region annotations remain from the recovered input, alongside their historical versions.
#
# The exported population is every recovered CDS position after the existing input coverage filter. `MismatchFrequency = AltBaseCount / (RefBaseCount + AltBaseCount)`, as in the octopus analysis, not the sum of four counts. Alternative bases are selected without A>G preference. Positive ties use a local RNG seeded with 1892; choices are made once here and exported. No observed alternative means no mismatch label; absent canonical coverage remains undefined, not an invented observation.

# %%
mm12_positions_dfs = []
mm12_merge_rows = []
mm12_original_only_sites_dfs = []
mm12_original_annotation_cols = [
    "RefBase", "TotalCoverage", "CDS", "Edited", "KnownEditing", "InProbRegion", "Noise", "EditingFrequency",
]
for mm12_recovered_df, mm12_original_df, mm12_input in zip(
    mm12_recovered_dfs, mm12_original_dfs, mm12_inputs_df.itertuples(index=False)
):
    mm12_recovered_keys = pd.MultiIndex.from_frame(mm12_recovered_df[["Chrom", "Position"]])
    mm12_original_keys = pd.MultiIndex.from_frame(mm12_original_df[["Chrom", "Position"]])
    mm12_original_only = ~mm12_original_keys.isin(mm12_recovered_keys)
    mm12_original_only_sites_dfs.append(mm12_original_df.loc[mm12_original_only, mm12_site_cols].copy())
    mm12_annotations_df = mm12_original_df[["Chrom", "Position"] + mm12_original_annotation_cols].rename(
        columns={col: f"Original{col}" for col in mm12_original_annotation_cols}
    )
    mm12_positions_df = mm12_recovered_df.rename(columns={
        "Edited": "RecoveredEdited", "Noise": "RecoveredNoise", "EditingFrequency": "RecoveredEditingFrequency",
    }).merge(mm12_annotations_df, on=["Chrom", "Position"], how="left", validate="one_to_one", indicator=True)
    assert len(mm12_positions_df) == len(mm12_recovered_df), "Recovered rows lost or duplicated"
    assert np.array_equal(mm12_positions_df[["A", "T", "C", "G"]], mm12_recovered_df[["A", "T", "C", "G"]])
    mm12_positions_df["OriginalAnnotationPresent"] = mm12_positions_df.pop("_merge").eq("both")
    mm12_matched = mm12_positions_df["OriginalAnnotationPresent"]
    assert mm12_positions_df.loc[mm12_matched, "RefBase"].eq(mm12_positions_df.loc[mm12_matched, "OriginalRefBase"]).all()
    mm12_merge_rows.append({
        **{col: getattr(mm12_input, col) for col in mm12_unit_cols},
        "RecoveredPositions": len(mm12_positions_df), "OriginalPositions": len(mm12_original_df),
        "RecoveredOnlyPositions": int((~mm12_matched).sum()),
        "OriginalOnlyPositions": int(mm12_original_only.sum()),
        "RecoveredCDSPositions": int(mm12_positions_df["CDS"].sum()),
    })
    mm12_positions_df["SourceFile"] = str(mm12_input.SourceFile)
    mm12_positions_df["OriginalSourceFile"] = str(mm12_input.OriginalSourceFile)
    mm12_positions_dfs.append(mm12_positions_df.loc[mm12_positions_df["CDS"]].copy())
mm12_positions_df = pd.concat(mm12_positions_dfs, ignore_index=True)
mm12_positions_df = mm12_positions_df.merge(mm12_thresholds_df, on=mm12_unit_cols, how="left", validate="many_to_one")
assert mm12_positions_df["OldEditingThreshold"].notna().all()
assert mm12_positions_df["OriginalSchemeValidated"].eq(True).all()
mm12_merge_summary_df = pd.DataFrame(mm12_merge_rows)
mm12_original_only_sites_df = pd.concat(mm12_original_only_sites_dfs, ignore_index=True)

# %%
# One stable input/site ordering and one local draw for each positive tie.
mm12_positions_df = mm12_positions_df.sort_values(mm12_site_cols).reset_index(drop=True)
mm12_bases = np.array(list("ATCG"))
mm12_base_counts = mm12_positions_df[list(mm12_bases)].to_numpy()
mm12_ref_indices = pd.Index(mm12_bases).get_indexer(mm12_positions_df["RefBase"])
mm12_row_indices = np.arange(len(mm12_positions_df))
mm12_ref_counts = mm12_base_counts[mm12_row_indices, mm12_ref_indices]
mm12_alternative_counts = mm12_base_counts.copy()
mm12_alternative_counts[mm12_row_indices, mm12_ref_indices] = -1
mm12_alt_counts = mm12_alternative_counts.max(axis=1)
mm12_alt_indices = mm12_alternative_counts.argmax(axis=1)
mm12_candidates = mm12_alternative_counts == mm12_alt_counts[:, None]
mm12_positive_ties = (mm12_alt_counts > 0) & (mm12_candidates.sum(axis=1) > 1)
mm12_rng = np.random.default_rng(mm12_seed)
for mm12_row in np.flatnonzero(mm12_positive_ties):
    mm12_alt_indices[mm12_row] = mm12_rng.choice(np.flatnonzero(mm12_candidates[mm12_row]))
mm12_positions_df["RefBaseCount"] = mm12_ref_counts
mm12_positions_df["AltBaseCount"] = mm12_alt_counts
mm12_positions_df["AltBase"] = pd.Series(mm12_bases[mm12_alt_indices], dtype="string").where(mm12_alt_counts > 0)
mm12_positions_df["AlternativeTie"] = mm12_positive_ties
mm12_positions_df["AlternativeSeed"] = mm12_seed
mm12_denominator = mm12_ref_counts + mm12_alt_counts
mm12_positions_df["MismatchFrequency"] = np.divide(
    mm12_alt_counts, mm12_denominator,
    out=np.full(len(mm12_positions_df), np.nan), where=mm12_denominator > 0,
)
mm12_positions_df["ObservedMismatch"] = (mm12_alt_counts > 0) & (mm12_denominator > 0)
mm12_positions_df["Mismatch"] = (
    mm12_positions_df["RefBase"].astype("string") + ">" + mm12_positions_df["AltBase"]
).where(mm12_positions_df["ObservedMismatch"])
mm12_positions_df["MismatchStatus"] = np.select(
    [mm12_denominator == 0, mm12_alt_counts == 0],
    ["NoCanonicalCoverage", "NoAlternativeObserved"], default="Observed",
)

# %%
# Decide once from the prepared full-precision values; export and plots reuse these flags.
mm12_positions_df["AboveOldEditingThreshold"] = mm12_positions_df["MismatchFrequency"].gt(
    mm12_positions_df["OldEditingThreshold"]
)
mm12_positions_df["AtOrAboveSuspectedSNPLevel"] = mm12_positions_df["MismatchFrequency"].ge(0.10)
mm12_positions_df["MismatchFrequency1"] = mm12_positions_df["MismatchFrequency"].eq(1.0)


# %%
mm12_positions_df["Species"] = "D. pealeii"
mm12_positions_df["SNPThreshold"] = mm12_snp_threshold
mm12_positions_df["AnalysisEligible"] = mm12_positions_df["CDS"]
mm12_observed = mm12_positions_df["ObservedMismatch"]
mm12_is_ag = mm12_positions_df["Mismatch"].eq("A>G").fillna(False)
mm12_positions_df["DefinitiveSNP"] = mm12_observed & mm12_positions_df["MismatchFrequency1"]
mm12_positions_df["SuspectedSNP"] = (
    mm12_observed & ~mm12_is_ag & mm12_positions_df["AtOrAboveSuspectedSNPLevel"]
    & ~mm12_positions_df["MismatchFrequency1"]
)
mm12_positions_df["SuspectedSNPWoACOrAT"] = (
    mm12_positions_df["SuspectedSNP"] & ~mm12_positions_df["Mismatch"].isin(["A>C", "A>T"])
)
mm12_snp_units_df = mm12_positions_df.groupby(mm12_unit_cols, as_index=False).agg(
    NumOfSuspectedSNPsPerChrom=("SuspectedSNP", "sum"),
    NumOfSuspectedSNPsWithoutACOrATPerChrom=("SuspectedSNPWoACOrAT", "sum"),
)
mm12_snp_units_df["EditingDetectionDisabledInChromDueToSuspectedSNPsWithoutACOrAT"] = (
    mm12_snp_units_df["NumOfSuspectedSNPsWithoutACOrATPerChrom"].gt(3)
)
mm12_positions_df = mm12_positions_df.merge(mm12_snp_units_df, on=mm12_unit_cols, how="left", validate="many_to_one")
mm12_positions_df["NoiseSite"] = (
    mm12_observed & ~mm12_is_ag & ~mm12_positions_df["DefinitiveSNP"] & ~mm12_positions_df["SuspectedSNP"]
)
mm12_positions_df["EditingSite"] = (
    mm12_observed & mm12_is_ag & ~mm12_positions_df["MismatchFrequency1"]
    & ~mm12_positions_df["EditingDetectionDisabledInChromDueToSuspectedSNPsWithoutACOrAT"]
    & mm12_positions_df["AboveOldEditingThreshold"]
)


# %%
mm12_positions_df

# %% [markdown]
# ## Checks
#
# The gene decision is per dataset/platform/sample/gene, before selecting editing sites. Suspected SNP counts include AC/AT; the disqualification count excludes AG/AC/AT and excludes frequency 1. A passing gene need not contain an editing site. KS1/EE2/EE3 are questions tested below, with no expected-zero assertion. The summary and affected-unit table retain both SNP counts; the full site table remains available for inspection.

# %%
assert not mm12_positions_df.duplicated(mm12_site_cols).any(), "Duplicate exported site identity"
assert mm12_positions_df["CDS"].eq(True).all() and mm12_positions_df["AnalysisEligible"].eq(True).all()
assert mm12_positions_df["Strand"].eq("+").all(), "Current manifest is explicitly plus-strand"
assert mm12_positions_df.groupby(mm12_unit_cols)["OldEditingThreshold"].nunique().eq(1).all()
assert mm12_positions_df["OldEditingThreshold"].between(0, 1).all()
assert mm12_positions_df["SNPThreshold"].eq(0.10).all()
mm12_defined = mm12_positions_df["MismatchFrequency"].notna()
assert mm12_positions_df.loc[mm12_defined, "MismatchFrequency"].between(0, 1).all()
assert mm12_positions_df.loc[~mm12_defined, "MismatchStatus"].eq("NoCanonicalCoverage").all()
assert mm12_positions_df.loc[~mm12_positions_df["ObservedMismatch"], "Mismatch"].isna().all()
assert mm12_positions_df.loc[mm12_positions_df["ObservedMismatch"], "AltBaseCount"].gt(0).all()
assert mm12_positions_df.loc[mm12_positions_df["ObservedMismatch"], "RefBase"].ne(
    mm12_positions_df.loc[mm12_positions_df["ObservedMismatch"], "AltBase"]
).all()
assert not (mm12_positions_df["DefinitiveSNP"] & mm12_positions_df["SuspectedSNP"]).any()
assert not mm12_positions_df.loc[mm12_positions_df["SuspectedSNPWoACOrAT"], "Mismatch"].isin(["A>G", "A>C", "A>T"]).any()
assert mm12_positions_df.loc[mm12_positions_df["DefinitiveSNP"], "MismatchFrequency"].eq(1).all()
assert mm12_positions_df.loc[mm12_positions_df["SuspectedSNP"], "MismatchFrequency"].ge(0.10).all()
assert mm12_positions_df.loc[mm12_positions_df["SuspectedSNP"], "MismatchFrequency"].lt(1).all()
mm12_unit_coverage_df = mm12_inputs_df[mm12_unit_cols].merge(
    mm12_positions_df.groupby(mm12_unit_cols).size().rename("ExportedPositions").reset_index(),
    on=mm12_unit_cols, how="left", validate="one_to_one",
)
assert mm12_unit_coverage_df["ExportedPositions"].gt(0).all(), "A current input has no recovered CDS positions"
assert mm12_unit_coverage_df["ExportedPositions"].sum() == mm12_merge_summary_df["RecoveredCDSPositions"].sum()
display(mm12_merge_summary_df)

# %%
mm12_snp_summary_df = mm12_snp_units_df.groupby(["Dataset", "Platform"], as_index=False).agg(
    SuspectedSNPSites=("NumOfSuspectedSNPsPerChrom", "sum"),
    DisqualifyingSNPSites=("NumOfSuspectedSNPsWithoutACOrATPerChrom", "sum"),
    DecisionUnits=("Chrom", "size"),
    DisqualifiedUnits=("EditingDetectionDisabledInChromDueToSuspectedSNPsWithoutACOrAT", "sum"),
)
mm12_units_with_snps_df = mm12_snp_units_df.assign(
    HasSuspectedSNP=mm12_snp_units_df["NumOfSuspectedSNPsPerChrom"].gt(0),
    HasDisqualifyingSNP=mm12_snp_units_df["NumOfSuspectedSNPsWithoutACOrATPerChrom"].gt(0),
).groupby(["Dataset", "Platform"], as_index=False)[["HasSuspectedSNP", "HasDisqualifyingSNP"]].sum()
mm12_snp_summary_df = mm12_snp_summary_df.merge(
    mm12_units_with_snps_df, on=["Dataset", "Platform"], validate="one_to_one",
)
mm12_snp_cases_df = mm12_snp_units_df.loc[
    mm12_snp_units_df["NumOfSuspectedSNPsPerChrom"].gt(0)
].sort_values(
    ["EditingDetectionDisabledInChromDueToSuspectedSNPsWithoutACOrAT", "Dataset", "Sample"],
    ascending=[False, True, True],
)
mm12_suspected_sites_df = mm12_positions_df.loc[
    mm12_positions_df["SuspectedSNP"],
    mm12_site_cols + ["Mismatch", "MismatchFrequency", "SuspectedSNPWoACOrAT"],
].copy()
display(mm12_snp_summary_df)
display(mm12_snp_cases_df)  # At most 21 current input units; both SNP counts stay visible.

# %% [markdown]
# ## Original and updated sites
#
# Compare the original CDS editing sites with `EditingSite` in the updated scheme using the old thresholds. Both sets retain their own source coverage selection. Differences are reported, not forced to zero; this is separate from the unchanged original-threshold acceptance check.
#

# %%
mm12_comparison_keys = ["Dataset", "Platform", "Sample", "Chrom", "Position"]
mm12_original_sites_df = mm12_original_identity_df.loc[
    mm12_original_identity_df["OriginalSite"] & mm12_original_identity_df["CDS"], mm12_comparison_keys
]
mm12_updated_sites_df = mm12_positions_df.loc[mm12_positions_df["EditingSite"], mm12_comparison_keys]
assert not mm12_original_sites_df.duplicated(mm12_comparison_keys).any()
assert not mm12_updated_sites_df.duplicated(mm12_comparison_keys).any()
mm12_site_comparison_df = mm12_original_sites_df.merge(
    mm12_updated_sites_df, on=mm12_comparison_keys, how="outer", validate="one_to_one", indicator=True
)
mm12_site_comparison_df = mm12_site_comparison_df.assign(
    Original=mm12_site_comparison_df["_merge"].ne("right_only"),
    Updated=mm12_site_comparison_df["_merge"].ne("left_only"),
    Shared=mm12_site_comparison_df["_merge"].eq("both"),
    OriginalOnly=mm12_site_comparison_df["_merge"].eq("left_only"),
    UpdatedOnly=mm12_site_comparison_df["_merge"].eq("right_only"),
)
mm12_comparison_counts = ["Original", "Updated", "Shared", "OriginalOnly", "UpdatedOnly"]
mm12_comparison_units = mm12_comparison_keys[:-1]
mm12_scheme_summary_df = mm12_inputs_df[mm12_comparison_units + ["Gene"]].merge(
    mm12_site_comparison_df.groupby(mm12_comparison_units, as_index=False)[mm12_comparison_counts].sum(),
    on=mm12_comparison_units, how="left", validate="one_to_one",
)
mm12_scheme_summary_df[mm12_comparison_counts] = mm12_scheme_summary_df[mm12_comparison_counts].fillna(0).astype(int)
mm12_scheme_differences_df = mm12_site_comparison_df.loc[
    mm12_site_comparison_df["OriginalOnly"] | mm12_site_comparison_df["UpdatedOnly"]
].drop(columns="_merge")
display(mm12_scheme_summary_df)
if not mm12_scheme_differences_df.empty:
    display(mm12_scheme_differences_df.head(30))  # Full identities remain in this table.


# %% [markdown]
# ## Export
#
# Export all prepared CDS rows to `Code/Notebooks/12mm/v2/positions/squid_long_reads` and `squid_short_reads`. The source decisions are `AboveOldEditingThreshold`, `AtOrAboveSuspectedSNPLevel`, `MismatchFrequency1` and the gene-exclusion flag. Original annotations and the successful squid threshold check remain separate.
#
# The cells below check all destinations, write gzip files exclusively as `.pending`, read them back, and compare full site identities, decisions and all twelve counts in each panel. Numeric comparisons are diagnostic only. Only a successful transfer is published under the final names; failed pending files remain for inspection and are never loaded by the plotting notebook.
#

# %%
import os

# Use the already prepared annotations; do not choose alternatives or rebuild thresholds here.
assert mm12_identity_differences_df.empty, "The original squid threshold check must pass."
assert mm12_positions_df["OriginalSchemeValidated"].eq(True).all()
mm12_source_df = mm12_positions_df.copy()
mm12_output_root = Path("/private6/projects/Combinatorics/Code/Notebooks/12mm/v2/positions")
mm12_inputs_df["OutputFile"] = [
    mm12_output_root / ("squid_long_reads" if row.Platform == "PacBio" else "squid_short_reads")
    / f"{row.Stem}.12mm.positions.csv.gz"
    for row in mm12_inputs_df.itertuples(index=False)
]
mm12_output_paths = mm12_inputs_df["OutputFile"].tolist()
mm12_pending_paths = [Path(str(path) + ".pending") for path in mm12_output_paths]
mm12_existing_outputs = [path for path in mm12_output_paths + mm12_pending_paths if path.exists()]
if mm12_existing_outputs:
    raise FileExistsError(f"Refusing to overwrite: {mm12_existing_outputs}")


# %%
def mm12_panel_masks(positions_df):
    """Use producer decisions for the six panel populations."""
    base = positions_df["AnalysisEligible"] & positions_df["ObservedMismatch"]
    passed = ~positions_df["EditingDetectionDisabledInChromDueToSuspectedSNPsWithoutACOrAT"]
    above = positions_df["AboveOldEditingThreshold"]
    snp_level = positions_df["AtOrAboveSuspectedSNPLevel"]
    complete = positions_df["MismatchFrequency1"]
    return {
        "1": base & complete,
        "2 Passed": base & snp_level & ~complete & passed,
        "2 Excluded": base & snp_level & ~complete & ~passed,
        "3": base & passed & ~above,
        "4": base & passed & above & ~snp_level,
        "5": base & passed & above & ~complete,
    }



# %%
mm12_transfer_keys = ["Dataset", "Platform", "Sample", "Chrom", "Position"]
mm12_transfer_flags = [
    "AboveOldEditingThreshold", "AtOrAboveSuspectedSNPLevel", "MismatchFrequency1",
    "EditingDetectionDisabledInChromDueToSuspectedSNPsWithoutACOrAT",
    "ObservedMismatch", "AnalysisEligible", "EditingSite",
]
mm12_transfer_flags += [
    "CDS", "RecoveredEdited", "KnownEditing", "InProbRegion",
    "OriginalSchemeValidated", "DefinitiveSNP", "SuspectedSNP",
    "SuspectedSNPWoACOrAT", "NoiseSite", "OriginalAnnotationPresent", "AlternativeTie",
]
# Missing historical annotations remain missing; decision flags above cannot be missing.
mm12_nullable_flags = ["OriginalCDS", "OriginalEdited", "OriginalKnownEditing", "OriginalInProbRegion"]
mm12_categories = [
    "Species", "Gene", "Strand", "RefBase", "AltBase", "Mismatch", "MismatchStatus",
    "OriginalRefBase", "SourceFile", "OriginalSourceFile", "ThresholdSource",
]
mm12_count_cols = [
    "A", "T", "C", "G", "TotalCoverage", "RefBaseCount", "AltBaseCount",
    "OriginalTotalCoverage", "AlternativeSeed", "NumOfSuspectedSNPsPerChrom",
    "NumOfSuspectedSNPsWithoutACOrATPerChrom",
]
mm12_mismatches = sorted(f"{ref}>{alt}" for ref in "ATCG" for alt in "ATCG" if ref != alt)
if mm12_source_df[mm12_transfer_keys].isna().any().any() or mm12_source_df.duplicated(mm12_transfer_keys).any():
    raise ValueError("The export source has missing or duplicate site keys.")
for mm12_col in mm12_transfer_flags:
    if not pd.api.types.is_bool_dtype(mm12_source_df[mm12_col]) or mm12_source_df[mm12_col].isna().any():
        raise ValueError(f"The source requires explicit, nonmissing Boolean decisions: {mm12_col}")
mm12_source_df = mm12_source_df.set_index(mm12_transfer_keys).sort_index()
mm12_source_masks = mm12_panel_masks(mm12_source_df)
mm12_ag_sites = mm12_source_masks["5"] & mm12_source_df["Mismatch"].eq("A>G")
mm12_editing_sites = (
    mm12_source_df["AnalysisEligible"] & mm12_source_df["ObservedMismatch"]
    & mm12_source_df["EditingSite"]
)
if not mm12_source_df.index[mm12_ag_sites].equals(mm12_source_df.index[mm12_editing_sites]):
    display(mm12_source_df.loc[mm12_ag_sites ^ mm12_editing_sites].reset_index().head(20))
    raise ValueError("Graph 5 A>G identities disagree with EditingSite in the source scheme.")


# %%
mm12_transfer_verified = False
mm12_loaded_dfs, mm12_default_dfs = [], []
for mm12_input, mm12_pending_path in zip(mm12_inputs_df.itertuples(index=False), mm12_pending_paths):
    mm12_export_mask = (
        (mm12_source_df.index.get_level_values("Dataset") == mm12_input.Dataset)
        & (mm12_source_df.index.get_level_values("Platform") == mm12_input.Platform)
        & (mm12_source_df.index.get_level_values("Sample") == mm12_input.Sample)
        & (mm12_source_df.index.get_level_values("Chrom") == mm12_input.Chrom)
    )
    mm12_export_df = mm12_source_df.loc[mm12_export_mask].reset_index()
    assert not mm12_export_df.empty, f"No positions for {mm12_input.Sample}"
    mm12_pending_path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(mm12_pending_path, "xt", encoding="utf-8", newline="") as mm12_handle:
        mm12_export_df.to_csv(mm12_handle, sep="\t", index=False)
    mm12_loaded_dfs.append(pd.read_csv(
        mm12_pending_path, sep="\t", compression="gzip", float_precision="round_trip",
        dtype={col: "boolean" for col in mm12_transfer_flags + mm12_nullable_flags},
    ))
    mm12_default_dfs.append(pd.read_csv(
        mm12_pending_path, sep="\t", compression="gzip",
        usecols=mm12_transfer_keys + ["MismatchFrequency", "OldEditingThreshold"],
    ))
mm12_loaded_df = pd.concat(mm12_loaded_dfs, ignore_index=True)
mm12_default_df = pd.concat(mm12_default_dfs, ignore_index=True)


# %%
# Compare discrete values exactly; numerical threshold diagnostics are separate below.
if mm12_loaded_df[mm12_transfer_keys].isna().any().any() or mm12_loaded_df.duplicated(mm12_transfer_keys).any():
    raise ValueError("Read-back has missing or duplicate site keys; no files are published.")
mm12_loaded_df = mm12_loaded_df.set_index(mm12_transfer_keys).sort_index()
mm12_missing_sites = mm12_source_df.index.difference(mm12_loaded_df.index)
mm12_extra_sites = mm12_loaded_df.index.difference(mm12_source_df.index)
if len(mm12_missing_sites) or len(mm12_extra_sites):
    display(mm12_missing_sites.to_frame(index=False).head(20), mm12_extra_sites.to_frame(index=False).head(20))
    raise ValueError("Source/read-back populations differ; inspect inputs or run state before plotting.")
mm12_compare_cols = mm12_categories + mm12_count_cols + mm12_transfer_flags + mm12_nullable_flags
mm12_same = (
    mm12_source_df[mm12_compare_cols].eq(mm12_loaded_df[mm12_compare_cols])
    | (mm12_source_df[mm12_compare_cols].isna() & mm12_loaded_df[mm12_compare_cols].isna())
).fillna(False)
mm12_transfer_differences_df = mm12_source_df.loc[~mm12_same.all(axis=1), mm12_compare_cols].join(
    mm12_loaded_df.loc[~mm12_same.all(axis=1), mm12_compare_cols], lsuffix="_source", rsuffix="_read"
)
if not mm12_transfer_differences_df.empty:
    display(mm12_transfer_differences_df.reset_index().head(20))
    raise ValueError("Source decisions, categories or counts changed in transfer; no files are published.")


# %%
mm12_loaded_masks = mm12_panel_masks(mm12_loaded_df)
mm12_count_index = pd.MultiIndex.from_product(
    [mm12_source_df.index.get_level_values("Dataset").unique(), mm12_mismatches],
    names=["Dataset", "Mismatch"],
)
mm12_transfer_counts_dfs = []
for mm12_panel, mm12_mask in mm12_source_masks.items():
    mm12_read_mask = mm12_loaded_masks[mm12_panel]
    mm12_source_sites = mm12_source_df.index[mm12_mask]
    mm12_read_sites = mm12_loaded_df.index[mm12_read_mask]
    if not mm12_source_sites.equals(mm12_read_sites):
        display(mm12_source_sites.symmetric_difference(mm12_read_sites).to_frame(index=False).head(20))
        raise ValueError(f"Panel {mm12_panel}: source/read-back site identities differ.")
    mm12_counts_df = pd.DataFrame({
        "Source": mm12_source_df.loc[mm12_mask].groupby(["Dataset", "Mismatch"]).size().reindex(mm12_count_index, fill_value=0),
        "ReadBack": mm12_loaded_df.loc[mm12_read_mask].groupby(["Dataset", "Mismatch"]).size().reindex(mm12_count_index, fill_value=0),
    })
    assert mm12_counts_df["Source"].eq(mm12_counts_df["ReadBack"]).all(), mm12_panel
    mm12_counts_df["Panel"] = mm12_panel
    mm12_transfer_counts_dfs.append(mm12_counts_df.reset_index())
mm12_transfer_counts_df = pd.concat(mm12_transfer_counts_dfs, ignore_index=True)
display(mm12_transfer_counts_df.groupby(["Dataset", "Panel"], sort=False)[["Source", "ReadBack"]].sum())
# The full twelve-category integer count table remains available above.


# %%
# Diagnose parser changes separately from source population/decision changes.
# The default-parser table is never used for plotting or filter decisions.
mm12_default_df = mm12_default_df.set_index(mm12_transfer_keys).sort_index()
mm12_numeric_cols = ["MismatchFrequency", "OldEditingThreshold"]
mm12_saved_decision = mm12_source_df["AboveOldEditingThreshold"]
mm12_source_numeric = mm12_source_df["MismatchFrequency"].gt(mm12_source_df["OldEditingThreshold"])
mm12_read_numeric = mm12_loaded_df["MismatchFrequency"].gt(mm12_loaded_df["OldEditingThreshold"])
mm12_default_numeric = mm12_default_df["MismatchFrequency"].gt(mm12_default_df["OldEditingThreshold"])
mm12_numeric_same = (
    mm12_source_df[mm12_numeric_cols].eq(mm12_loaded_df[mm12_numeric_cols])
    | (mm12_source_df[mm12_numeric_cols].isna() & mm12_loaded_df[mm12_numeric_cols].isna())
).all(axis=1)
mm12_diagnostic_mask = (
    ~mm12_numeric_same | mm12_saved_decision.ne(mm12_read_numeric)
    | mm12_saved_decision.ne(mm12_default_numeric) | mm12_saved_decision.ne(mm12_source_numeric)
)
mm12_float_differences_df = mm12_source_df.loc[mm12_diagnostic_mask, ["Mismatch"]].copy()
for mm12_label, mm12_frame in [
    ("Source", mm12_source_df), ("RoundTripRead", mm12_loaded_df), ("DefaultRead", mm12_default_df),
]:
    for mm12_col in mm12_numeric_cols:
        mm12_float_differences_df[f"{mm12_label}_{mm12_col}"] = mm12_frame.loc[mm12_diagnostic_mask, mm12_col]
mm12_float_differences_df["SavedDecision"] = mm12_saved_decision.loc[mm12_diagnostic_mask]
mm12_float_differences_df["SourceNumericDecision"] = mm12_source_numeric.loc[mm12_diagnostic_mask]
mm12_float_differences_df["RoundTripNumericDecision"] = mm12_read_numeric.loc[mm12_diagnostic_mask]
mm12_float_differences_df["DefaultNumericDecision"] = mm12_default_numeric.loc[mm12_diagnostic_mask]
mm12_float_differences_df["RoundTripValuesChanged"] = ~mm12_numeric_same.loc[mm12_diagnostic_mask]
if not mm12_float_differences_df.empty:
    with pd.option_context("display.float_format", lambda value: format(value, ".17g")):
        display(mm12_float_differences_df.reset_index().head(20))
    print("Float/source-state differences are diagnostic only; the saved decisions remain authoritative.")

mm12_transfer_verified = True


# %%
assert mm12_transfer_verified, "Complete all transfer checks before publishing."
for mm12_pending_path, mm12_output_path in zip(mm12_pending_paths, mm12_output_paths):
    os.link(mm12_pending_path, mm12_output_path)  # Exclusive even if a target appeared meanwhile.
    mm12_pending_path.unlink()  # Only this successful export's temporary name.

