# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% papermill={"duration": 0.071769, "end_time": "2022-02-01T09:42:43.049672", "exception": false, "start_time": "2022-02-01T09:42:42.977903", "status": "completed"} tags=["parameters"]
condition_col = "Gene"
conditions = [
    "RUSC2_MOUSE",
    "TRIM2_BOVIN",
    "CA2D3_MOUSE",
    "ABL_DROME",
    "DGLA_HUMAN",
    "K0513_MOUSE",
    "KCNAS_DROME",
    "ACHA4_MOUSE",
    "ANR17_HUMAN",
    "TWK7_CAEEL",
    "SCN1_HETBL",
    "CACB2_RABIT",
    "RIMS2_RAT",
    "PCLO_CHICK",
    "DOP1_HUMAN",
    "IQEC1_HUMAN",
    "CSKI1_MOUSE",
    "MTUS2_HUMAN",
    "ROBO2_HUMAN",
]
chroms = [
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
starts = [
    192,
    352,
    1,
    2,
    171,
    400,
    1,
    385,
    2,
    2,
    443,
    288,
    311,
    1,
    256,
    989,
    71,
    2,
    451,
]
ends = [
    6296,
    2766,
    3495,
    2377,
    3359,
    3732,
    1617,
    2190,
    3460,
    1246,
    7381,
    1808,
    2116,
    6294,
    7722,
    4195,
    2857,
    4558,
    5283,
]
strands = [
    "+",
    "+",
    "+",
    "+",
    "+",
    "+",
    "+",
    "+",
    "+",
    "+",
    "+",
    "+",
    "+",
    "+",
    "+",
    "+",
    "+",
    "+",
    "+",
]
# unaligned_bam_files = [
#     "/private7/projects/Combinatorics/D.pealeii/Data/CCS/BasicCCS/GRIA-CNS-RESUB.C0x1291.ccs.bam",
#     "/private7/projects/Combinatorics/D.pealeii/Data/CCS/BasicCCS/PCLO-CNS-RESUB.C0x1291.ccs.bam",
# ]
reads_type = "miseq"  # something like CCS / miseq / etc.
aligned_bam_files = [
    f"/private7/projects/Combinatorics/D.pealeii/Alignment/Illumina/reads.ByChrom/reads.sorted.aligned.filtered.{chrom}.bam"
    for chrom in chroms
]
# filtered_aligned_bam_files = [
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.bam",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.bam",
# ]
include_flags = None
exclude_flags = "2304"  # remove secondary and supplementary (chimeric) alignments
sep = "\t"
positions_files = [
    f"/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.{chrom}.positions.csv"
    for chrom in chroms
]
reads_files = [
    f"/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.{chrom}.reads.csv"
    for chrom in chroms
]
unique_reads_files = [
    f"/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.{chrom}.unique_reads.csv"
    for chrom in chroms
]
proteins_files = [
    f"/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.{chrom}.proteins.csv"
    for chrom in chroms
]
unique_proteins_files = [
    f"/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.{chrom}.unique_proteins.csv"
    for chrom in chroms
]
reads_first_col_pos = 6
unique_reads_first_col_pos = 8
proteins_first_col_pos = 12
unique_proteins_first_col_pos = 14
reads_editing_col = "EditingFrequency"
proteins_editing_col = "MinNonSyns"

# todo update & uncomment
# distinct_unique_reads_files = [
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA.AllRows.DistinctUniqueReads.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO.AllRows.DistinctUniqueReads.csv",
# ]

distinct_unique_proteins_files = [
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141881_c0_seq3.DistinctUniqueProteins.12.07.2022-20:54:38.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141044_c0_seq2.DistinctUniqueProteins.13.07.2022-06:33:23.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp140439_c0_seq1.DistinctUniqueProteins.12.07.2022-22:51:22.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp126362_c0_seq1.DistinctUniqueProteins.15.07.2022-06:11:18.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141517_c0_seq1.DistinctUniqueProteins.14.07.2022-07:43:15.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141840_c0_seq2.DistinctUniqueProteins.13.07.2022-20:30:25.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141640_c0_seq1.DistinctUniqueProteins.12.07.2022-19:44:02.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp140987_c3_seq1.DistinctUniqueProteins.18.07.2022-07:50:43.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp140910_c2_seq1.DistinctUniqueProteins.13.07.2022-16:15:35.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp136058_c0_seq1.DistinctUniqueProteins.21.07.2022-07:57:53.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141378_c0_seq7.DistinctUniqueProteins.19.07.2022-08:12:24.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141158_c1_seq2.DistinctUniqueProteins.13.07.2022-01:54:59.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp140712_c0_seq3.DistinctUniqueProteins.05.11.2022-02_11_33.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141882_c0_seq14.DistinctUniqueProteins.11.11.2022-22_23_54.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141880_c1_seq3.DistinctUniqueProteins.14.11.2022-20_58_20.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141565_c6_seq3.DistinctUniqueProteins.16.11.2022-12_42_59.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141684_c0_seq1.DistinctUniqueProteins.19.11.2022-01_51_21.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141532_c3_seq11.DistinctUniqueProteins.05.12.2022-20_38_20.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141574_c0_seq3.DistinctUniqueProteins.12.12.2022-19_31_01.csv",
]

# proteins_jaccard_files = [
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141881_c0_seq3.JaccardMatrixProteins.12.07.2022-20:54:38.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141044_c0_seq2.JaccardMatrixProteins.13.07.2022-06:33:23.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp140439_c0_seq1.JaccardMatrixProteins.12.07.2022-22:51:22.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp126362_c0_seq1.JaccardMatrixProteins.15.07.2022-06:11:18.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141517_c0_seq1.JaccardMatrixProteins.14.07.2022-07:43:15.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141840_c0_seq2.JaccardMatrixProteins.13.07.2022-20:30:25.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141640_c0_seq1.JaccardMatrixProteins.12.07.2022-19:44:02.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp140987_c3_seq1.JaccardMatrixProteins.18.07.2022-07:50:43.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp140910_c2_seq1.JaccardMatrixProteins.13.07.2022-16:15:35.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp136058_c0_seq1.JaccardMatrixProteins.21.07.2022-07:57:53.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141378_c0_seq7.JaccardMatrixProteins.19.07.2022-08:12:24.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141158_c1_seq2.JaccardMatrixProteins.13.07.2022-01:54:59.csv",
# ]

expression_files = [
    f"/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/{condition}.DistinctUniqueProteins.ExpressionLevels.csv"
    for condition in conditions
]

# todo update & uncomment
# alg_repetitions = 5

known_sites_file = (
    "/private7/projects/Combinatorics/D.pealeii/Annotations/D.pea.EditingSites.csv"
)
samtools_path = "/home/alu/kobish/anaconda3/envs/combinatorics/bin/samtools"
threads = 20
code_dir = "/private7/projects/Combinatorics/Code"
seed = 1892

# %%
len(chroms)

# %%
len(positions_files)

# %% [markdown] papermill={"duration": 0.029907, "end_time": "2022-02-01T09:42:43.198426", "exception": false, "start_time": "2022-02-01T09:42:43.168519", "status": "completed"}
# # Imports

# %%
# %load_ext autoreload
# %autoreload 2

# %% papermill={"duration": 2.901153, "end_time": "2022-02-01T09:42:46.125355", "exception": false, "start_time": "2022-02-01T09:42:43.224202", "status": "completed"}
import sys
from functools import reduce
from itertools import chain, combinations, product
from math import ceil
from multiprocessing import Pool
from pathlib import Path

from scipy import interpolate  # todo unimport this later?
import matplotlib.pyplot as plt
import numpy as np

# from numba import jit, njit, prange
import pandas as pd
import plotly.colors as pc
import plotly.express as px
import plotly.graph_objects as go
import scipy.stats
import seaborn as sns
from icecream import ic
from matplotlib_venn import venn2, venn3
from plotly.subplots import make_subplots
from sklearn import linear_model
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.metrics import mean_squared_error, r2_score

# from numpy.random import RandomState


sys.path.append(str(Path(code_dir).absolute()))
from Alignment.alignment_utils import (
    count_reads,
    count_reads_in_fastq,
    count_reads_in_unaligned_bam,
    count_unique_filtered_aligned_reads,
    count_unique_reads,
)


# %% [markdown] papermill={"duration": 0.040192, "end_time": "2022-02-01T09:42:46.214429", "exception": false, "start_time": "2022-02-01T09:42:46.174237", "status": "completed"}
# # Ploting utils

# %%
def rgb_change(r, g, b, d_r, d_g, d_b, scale):
    # todo: allow both changes to be in the same direction by modifying the given scale?
    values = [r, g, b]
    deltas = [int(d_v * scale) for d_v in (d_r, d_g, d_b)]
    legitimate_changes = {
        "+": [min(v + d_v, 255) for v, d_v in zip(values, deltas)],
        "-": [max(v - d_v, 0) for v, d_v in zip(values, deltas)],
    }
    complete_changes = {
        "+": sum(
            v + d_v == new_v
            for v, d_v, new_v in zip(values, deltas, legitimate_changes["+"])
        ),
        "-": sum(
            v - d_v == new_v
            for v, d_v, new_v in zip(values, deltas, legitimate_changes["-"])
        ),
    }
    if complete_changes["+"] >= complete_changes["-"]:
        r, g, b = legitimate_changes["+"]
    else:
        r, g, b = legitimate_changes["-"]
    return r, g, b


def two_subcolors_from_hex(hex_color, d_r=4, d_g=20, d_b=22, scale_1=1, scale_2=4):
    r, g, b = pc.hex_to_rgb(hex_color)
    subcolor_1 = rgb_change(r, g, b, d_r, d_g, d_b, scale_1)
    subcolor_2 = rgb_change(r, g, b, d_r, d_g, d_b, scale_2)
    subcolor_1 = pc.label_rgb(subcolor_1)
    subcolor_2 = pc.label_rgb(subcolor_2)
    return subcolor_1, subcolor_2


# %%
ic(conditions);

# %% papermill={"duration": 0.054755, "end_time": "2022-02-01T09:42:46.304499", "exception": false, "start_time": "2022-02-01T09:42:46.249744", "status": "completed"}
# plotly consts
# color_sequence = px.colors.qualitative.Pastel
# color_sequence = px.colors.qualitative.D3
# color_sequence = px.colors.qualitative.G10
color_sequence = px.colors.qualitative.Dark24
color_discrete_map = {
    condition: color for condition, color in zip(conditions, color_sequence)
}
ic(color_discrete_map)
subcolors_discrete_map = {
    condition: two_subcolors_from_hex(color_discrete_map[condition])
    for condition in conditions
}
ic(subcolors_discrete_map)
category_orders = {condition_col: conditions}
horizontal_category_orders = {
    category: list(reversed(category_orders[category])) for category in category_orders
}
# valid_shapes = ['', '/', '\\', 'x', '-', '|', '+', '.']
# pattern_shape_map = {
#     condition: shape for condition, shape in zip(conditions, cycle(valid_shapes))
# }
facet_col_spacing = 0.05
facet_col_wrap = 6
facet_row_spacing = facet_col_spacing * 6
template = "plotly_white"
zerolinewidth = 4


# %% [markdown]
#
# plotly.colors.n_colors(lowcolor, highcolor, n_colors, colortype='tuple')
#
#     Splits a low and high color into a list of n_colors colors in it
#
#     Accepts two color tuples and returns a list of n_colors colors which form the intermediate colors between lowcolor and highcolor from linearly interpolating through RGB space. If colortype is ‘rgb’ the function will return a list of colors in the same form.
#

# %%
def n_repetitions_colormap(subcolors_discrete_map, condition, n_repetitions):
    lowcolor, highcolor = subcolors_discrete_map[condition]
    colors = pc.n_colors(lowcolor, highcolor, n_repetitions, colortype="rgb")
    return {i: color for i, color in enumerate(colors, start=1)}


# %%
n_repetitions_colormap(subcolors_discrete_map, conditions[0], 10)

# %% [markdown] papermill={"duration": 0.040192, "end_time": "2022-02-01T09:42:46.214429", "exception": false, "start_time": "2022-02-01T09:42:46.174237", "status": "completed"}
# # Data

# %% [markdown] papermill={"duration": 0.02598, "end_time": "2022-02-01T09:42:46.438342", "exception": false, "start_time": "2022-02-01T09:42:46.412362", "status": "completed"}
# ## Known sites

# %%
known_sites_df = pd.read_csv(known_sites_file)
# new_known_sites_cols = ["Chrom", "SwissProt", "Position", "OriginalAA", "NewAA", "RefBase", "Editing", "%Editing", "DNA"]
new_known_sites_cols = [
    "Chrom",
    "SwissProt",
    "Position",
    "OriginalAA",
    "NewAA",
    "RefBase",
    "Editing",
    "EditingFrequency",
    "DNA",
]
known_sites_df = known_sites_df.set_axis(new_known_sites_cols, axis="columns")
known_sites_df.insert(
    5, "Coverage", known_sites_df["RefBase"] + known_sites_df["Editing"]
)
known_sites_df = known_sites_df.drop(["RefBase", "DNA"], axis="columns")
known_sites_df["Position"] = known_sites_df["Position"] - 1
# known_sites_df["%Editing"] = known_sites_df["%Editing"] * 100
known_sites_df


# %% [markdown] papermill={"duration": 0.041741, "end_time": "2022-02-01T09:42:47.760215", "exception": false, "start_time": "2022-02-01T09:42:47.718474", "status": "completed"}
# ## Positions

# %%
positions_dfs = [
    pd.read_csv(position_file, sep=sep) for position_file in positions_files
]
for positions_df, condition in zip(positions_dfs, conditions):
    positions_df.insert(0, condition_col, condition)
positions_dfs[0]


# %%
editing_positions_per_sample = [
    len(df.loc[(df["Edited"]) & (df["CDS"])]) for df in positions_dfs
]
for x in editing_positions_per_sample:
    print(x)

# %%
print(
    f"Average of {sum(editing_positions_per_sample)/len(positions_dfs)} editing sites per sample"
)

# %% [markdown] papermill={"duration": 0.02598, "end_time": "2022-02-01T09:42:46.438342", "exception": false, "start_time": "2022-02-01T09:42:46.412362", "status": "completed"}
# ## Reads

# %% [markdown] papermill={"duration": 0.02598, "end_time": "2022-02-01T09:42:46.438342", "exception": false, "start_time": "2022-02-01T09:42:46.412362", "status": "completed"}
# ### All

# %% [markdown]
# That is, all filtered reads.

# %% papermill={"duration": 1.204258, "end_time": "2022-02-01T09:42:47.668206", "exception": false, "start_time": "2022-02-01T09:42:46.463948", "status": "completed"}
reads_dfs = [pd.read_csv(reads_file, sep=sep) for reads_file in reads_files]
reads_dfs[0]


# %%
for reads_df in reads_dfs:
    mean_edited_positions = reads_df["EditedPositions"].mean()
    print(mean_edited_positions)

# %%
ambigous_positions_in_reads_df = reads_df = pd.concat(
    [reads_df["AmbigousPositions"] for reads_df in reads_dfs], ignore_index=True
)
ambigous_positions_in_reads_df

# %%
100 * ambigous_positions_in_reads_df.gt(0).sum() / len(ambigous_positions_in_reads_df)

# %%
ambigous_positions_in_reads_df.mean()

# %%
# edited_reads_dfs = [
#     reads_df.loc[reads_df[reads_editing_col] > 0] for reads_df in reads_dfs
# ]
# edited_reads_dfs[0]


# %% [markdown] papermill={"duration": 0.041741, "end_time": "2022-02-01T09:42:47.760215", "exception": false, "start_time": "2022-02-01T09:42:47.718474", "status": "completed"}
# ### Unique

# %% [markdown]
# All unique reads

# %% papermill={"duration": 0.126539, "end_time": "2022-02-01T09:42:47.923363", "exception": false, "start_time": "2022-02-01T09:42:47.796824", "status": "completed"}
unique_reads_dfs = [
    pd.read_csv(unique_reads_file, sep=sep) for unique_reads_file in unique_reads_files
]
for unique_reads_df in unique_reads_dfs:
    if "Transcript" in unique_reads_df.columns:
        unique_reads_df.rename(columns={"Transcript": "UniqueRead"}, inplace=True)
unique_reads_dfs[0]


# %%
unique_reads_dfs[1]

# %% [markdown]
# Each of these unique reads is edited at least once

# %%
# edited_unique_reads_dfs = [
#     unique_reads_df.loc[unique_reads_df[reads_editing_col] > 0]
#     for unique_reads_df in unique_reads_dfs
# ]
# edited_unique_reads_dfs[0]


# %% [markdown]
# ### Distinct unique reads

# %%
# todo update & uncomment

# distinct_unique_reads_dfs = []
# for condition, distinct_unique_reads_file, unique_reads_df in zip(
#     conditions, distinct_unique_reads_files, unique_reads_dfs
# ):
#     distinct_unique_reads_df = pd.read_csv(distinct_unique_reads_file, sep=sep)
#     distinct_unique_reads_df.insert(0, condition_col, condition)
#     distinct_unique_reads_df.insert(
#         1,
#         "NumOfReads",
#         (
#             distinct_unique_reads_df["Fraction"]
#             * unique_reads_df["NumOfReads"].sum()
#         ).astype(int),
#     )
#     distinct_unique_reads_dfs.append(distinct_unique_reads_df)

# distinct_unique_reads_df = (
#     pd.concat(distinct_unique_reads_dfs)
#     .reset_index(drop=True)
#     .rename(
#         columns={"NumUniqueSamples": "NumUniqueReads", "UniqueSamples": "UniqueReads"}
#     )
# )

# distinct_unique_reads_df


# %%
# todo update & uncomment

# expanded_distinct_unique_reads_df = (
#     distinct_unique_reads_df.copy()
#     .assign(UniqueReads2=lambda x: x.UniqueReads.str.split(","))
#     .drop("UniqueReads", axis=1)
#     .rename(columns={"UniqueReads2": "UniqueReads"})
#     .explode("UniqueReads")
#     .rename(columns={"UniqueReads": "UniqueRead", "NumOfReads": "NumOfReadsInFraction"})
#     .drop(["NumUniqueReads"], axis=1)
#     .merge(
#         pd.concat(
#             [
#                 df.loc[
#                     :,
#                     [
#                         condition_col,
#                         "UniqueRead",
#                         "NumOfReads",
#                         "Reads",
#                         "EditingFrequency",
#                         "EditedPositions",
#                     ],
#                 ]
#                 for df in unique_reads_dfs
#             ]
#         ),
#         on=[condition_col, "UniqueRead"],
#     )
#     .assign(Reads2=lambda x: x.Reads.str.split(","))
#     .drop("Reads", axis=1)
#     .rename(columns={"Reads2": "Reads"})
# )

# expanded_distinct_unique_reads_df


# %% [markdown]
# ## Proteins

# %% [markdown]
# ### All proteins

# %%
proteins_dfs = [pd.read_csv(proteins_file, sep=sep) for proteins_file in proteins_files]
for proteins_df in proteins_dfs:
    if "Transcript" in proteins_df.columns:
        proteins_df.rename(columns={"Transcript": "UniqueRead"}, inplace=True)
proteins_dfs[0]


# %%
proteins_dfs[1]

# %%
# edited_proteins_dfs = [
#     proteins_df.loc[proteins_df[proteins_editing_col] > 0]
#     for proteins_df in proteins_dfs
# ]
# edited_proteins_dfs[0]


# %% [markdown]
# ### Unique proteins

# %%
unique_proteins_dfs = [
    pd.read_csv(unique_proteins_file, sep=sep)
    for unique_proteins_file in unique_proteins_files
]
for unique_proteins_df in unique_proteins_dfs:
    unique_proteins_df.rename(
        columns={
            col: col.replace("Transcripts", "UniqueReads")
            for col in unique_proteins_df.columns[:unique_proteins_first_col_pos]
            if "Transcripts" in col
        },
        inplace=True,
    )
unique_proteins_dfs[0]


# %%
for unique_proteins_df in unique_proteins_dfs:
    mean_edited_aas = unique_proteins_df["MinNonSyns"].mean()
    print(mean_edited_aas)

# %%
unique_proteins_dfs[-1].head()

# %%
editable_aas_per_sample = [
    df.iloc[:, unique_proteins_first_col_pos:].shape[1] for df in unique_proteins_dfs
]

editable_aas_per_sample

# %%
# sorted(editable_aas_per_sample)

# %%
avg_editables_aas_per_sample = sum(editable_aas_per_sample) / len(unique_proteins_dfs)

print(f"Average of {avg_editables_aas_per_sample:.0f} editable AAs per sample")

# %%
unique_proteins_dfs[3].iloc[:, unique_proteins_first_col_pos:]

# %%
pd.DataFrame(
    {
        condition_col: conditions,
        "EditableAAs": [
            unique_proteins_df.iloc[:, unique_proteins_first_col_pos:].shape[1]
            for unique_proteins_df in unique_proteins_dfs
        ],
    }
)

# %%
unique_proteins_dfs[0].columns[:unique_proteins_first_col_pos]

# %%
# unique_edited_proteins_dfs = [
#     unique_proteins_df.loc[unique_proteins_df[proteins_editing_col] > 0]
#     for unique_proteins_df in unique_proteins_dfs
# ]
# unique_edited_proteins_dfs[0]


# %% [markdown]
# ### Distinct unique proteins

# %%
distinct_unique_proteins_dfs = []
for condition, distinct_unique_proteins_file, unique_reads_df in zip(
    conditions, distinct_unique_proteins_files, unique_reads_dfs
):
    distinct_unique_proteins_df = pd.read_csv(distinct_unique_proteins_file, sep=sep)
    distinct_unique_proteins_df.insert(0, condition_col, condition)
    distinct_unique_proteins_df.insert(
        1,
        "NumOfReads",
        (
            distinct_unique_proteins_df["Fraction"]
            * unique_reads_df["NumOfReads"].sum()
        ).astype(int),
    )
    distinct_unique_proteins_dfs.append(distinct_unique_proteins_df)

distinct_unique_proteins_df = (
    pd.concat(distinct_unique_proteins_dfs)
    .reset_index(drop=True)
    .rename(columns={"NumUniqueSamples": "NumOfProteins", "UniqueSamples": "Proteins"})
)

distinct_unique_proteins_df = distinct_unique_proteins_df.sort_values(
    [
        condition_col,
        "Fraction",
        "FractionRepetition",
        "Algorithm",
        "AlgorithmRepetition",
    ]
).reset_index(drop=True)

distinct_unique_proteins_df


# %%
distinct_unique_proteins_df.loc[
    distinct_unique_proteins_df["NumOfProteins"]
    == distinct_unique_proteins_df["NumOfProteins"].max()
]

# %%

# %%
# unique_edited_proteins_dfs[0].columns[:unique_proteins_first_col_pos]


# %%
expanded_distinct_unique_proteins_df = (
    distinct_unique_proteins_df.copy()
    .assign(Proteins2=lambda x: x.Proteins.str.split(","))
    .drop("Proteins", axis=1)
    .rename(columns={"Proteins2": "Proteins"})
    .explode("Proteins")
    .rename(columns={"Proteins": "Protein", "NumOfReads": "NumOfReadsInFraction"})
    .drop(["NumOfProteins"], axis=1)
    .merge(
        pd.concat(
            [df.iloc[:, :unique_proteins_first_col_pos] for df in unique_proteins_dfs]
        ),
        on=[condition_col, "Protein"],
    )
)

expanded_distinct_unique_proteins_df


# %%
distinct_unique_proteins_df2 = (
    expanded_distinct_unique_proteins_df.groupby(
        [
            condition_col,
            "Fraction",
            "FractionRepetition",
            "Algorithm",
            "AlgorithmRepetition",
        ]
    )["NumOfReads"]
    .sum()
    .reset_index()
    .rename(columns={"NumOfReads": "NumOfSupportingReads"})
    .merge(
        distinct_unique_proteins_df,
        on=[
            condition_col,
            "Fraction",
            "FractionRepetition",
            "Algorithm",
            "AlgorithmRepetition",
        ],
    )
    .assign(
        SupportingReadsPerProtein=lambda x: x["NumOfSupportingReads"]
        / x["NumOfProteins"],
        PercentSupportedReads=lambda x: 100
        * x["NumOfSupportingReads"]
        / x["NumOfReads"],
    )
    .rename(columns={"PercentSupportedReads": "%SupportedReads"})
)
distinct_unique_proteins_df2

# %% [markdown]
# #### Jaccards

# %%
# todo uncomment

# proteins_jaccard_dfs = []
# for jaccard_file in proteins_jaccard_files:
#     df = pd.read_csv(jaccard_file, sep=sep)
#     proteins_jaccard_dfs.append(df)
# proteins_jaccard_dfs[0]

# %%
# todo uncomment

# annotated_proteins_jaccard_dfs = []

# for condition, proteins_jaccard_df in zip(conditions, proteins_jaccard_dfs):

#     df = distinct_unique_proteins_df.loc[
#         distinct_unique_proteins_df[condition_col] == condition
#     ].reset_index(drop=True)
#     df = distinct_unique_proteins_df.loc[
#         distinct_unique_proteins_df[condition_col] == condition,
#         [
#             condition_col,
#             "Fraction",
#             "FractionRepetition",
#             "Algorithm",
#             "AlgorithmRepetition",
#         ],
#     ].reset_index(drop=True)
#     df = pd.concat([df, proteins_jaccard_df], axis=1)
#     index_details_dict = {
#         str(
#             i + 1
#         ): f'{row["Fraction"]}-{row["FractionRepetition"]}-{row["Algorithm"]}-{row["AlgorithmRepetition"]}'
#         for i, row in df.iterrows()
#     }
#     df = df.rename(columns=index_details_dict)

#     annotated_proteins_jaccard_dfs.append(df)

# annotated_proteins_jaccard_dfs[0]

# %% [markdown]
# ## Summary of data loss

# %%
total_unaligned_reads_fasta = "/private7/projects/Combinatorics/D.pealeii/Data/RutisReads/TrimmedWoDup/reads_1.fastq.gz"
total_aligned_reads_bam = (
    "/private7/projects/Combinatorics/D.pealeii/Alignment/Illumina/reads.sorted.bam"
)
total_aligned_filtered_reads_bam = "/private7/projects/Combinatorics/D.pealeii/Alignment/Illumina/reads.sorted.aligned.filtered.bam"
# per_chrom_aligned_bams = [
#     f"/private7/projects/Combinatorics/D.pealeii/Alignment/Illumina/reads.ByChrom/reads.sorted.aligned.filtered.{chrom}.bam"
#     for chrom in chroms
# ]

# %%

# %%
# todo update & uncomment

# unaligned_reads_counts = [
#     count_reads_in_unaligned_bam(samtools_path, bam, threads)
#     for bam in unaligned_bam_files
# ]

# aligned_reads_counts = [
#     count_reads(samtools_path, bam, f"{chrom}:{start+1}-{end}", include_flags, exclude_flags, threads)
#     for bam, chrom, start, end in zip(aligned_bam_files, chroms, starts, ends)
# ]

with Pool(processes=5) as pool:
    aligned_reads_counts = pool.starmap(
        func=count_unique_filtered_aligned_reads,
        iterable=[
            (bam, f"{chrom}:{start+1}-{end}", include_flags, exclude_flags, threads)
            for bam, chrom, start, end in zip(aligned_bam_files, chroms, starts, ends)
        ],
    )


# %%
conditions

# %%
aligned_reads_counts

# %%
sum(aligned_reads_counts)

# %%
total_aligned_filtered_reads - sum(aligned_reads_counts)

# %%
# filtered_aligned_reads_counts = [
#     count_reads(samtools_path, bam, f"{chrom}:{start+1}-{end}", include_flags, exclude_flags, threads)
#     for bam, chrom, start, end in zip(filtered_aligned_bam_files, chroms, starts, ends)
# ]

# pileup_reads_counts = [
#     len(set(chain.from_iterable(positions_df["Reads"].str.split(","))))
#     for positions_df in positions_dfs
# ]


def count_pileup_reads(positions_df):
    return len(set(chain.from_iterable(positions_df["Reads"].str.split(","))))


with Pool(processes=3) as pool:
    pileup_reads_counts = pool.map(count_pileup_reads, positions_dfs)

pileup_reads_counts


# %%
unique_reads_counts = [len(unique_reads_df) for unique_reads_df in unique_reads_dfs]

# max_fraction = distinct_unique_reads_df["Fraction"].max()
max_fraction = 1.0


# distinct_unique_reads_counts = (
#     distinct_unique_reads_df.loc[
#         distinct_unique_reads_df["Fraction"] == max_fraction
#     ]
#     .groupby(condition_col)["NumUniqueReads"]
#     .mean()
#     .round()
#     .astype(int)
# )


distinct_unique_proteins_counts = (
    distinct_unique_proteins_df.loc[
        distinct_unique_proteins_df["Fraction"] == max_fraction
    ]
    .groupby(condition_col)["NumOfProteins"]
    .mean()
    .round()
    .astype(int)
)

# todo update & uncomment

data_loss_df = pd.DataFrame(
    {
        # f"Unaligned {reads_type} reads": unaligned_reads_counts,
        "Aligned reads (within ORF)": aligned_reads_counts,
        # "Filtered aligned reads (within ORF)": filtered_aligned_reads_counts,
        "Pileup reads": pileup_reads_counts,
        "Unique reads": unique_reads_counts,
        # "Distinct unique reads (mean)": distinct_unique_reads_counts,
        "Distinct unique proteins (mean)": distinct_unique_proteins_counts,
    },
    index=conditions,
)
data_loss_df


# %%
total_unaligned_reads = count_reads_in_fastq(
    total_unaligned_reads_fasta, cat_cmd="zcat"
)
ic(total_unaligned_reads)

# total_aligned_reads = count_unique_reads(total_aligned_reads_bam, 10)
total_aligned_reads = count_unique_filtered_aligned_reads(
    total_aligned_reads_bam, None, None, None, threads
)
ic(total_aligned_reads)

ic(total_unaligned_reads - total_aligned_reads)

total_aligned_filtered_reads = count_unique_reads(total_aligned_filtered_reads_bam, 10)
ic(total_aligned_filtered_reads)

ic(total_aligned_reads - total_aligned_filtered_reads)

total_aligned_filtered_reads_within_targeted_transcripts = sum(aligned_reads_counts)
ic(total_aligned_filtered_reads_within_targeted_transcripts)

ic(
    total_aligned_filtered_reads
    - total_aligned_filtered_reads_within_targeted_transcripts
)

# %% [markdown] papermill={"duration": 0.045853, "end_time": "2022-02-01T09:42:48.953594", "exception": false, "start_time": "2022-02-01T09:42:48.907741", "status": "completed"}
# # Results

# %% [markdown] papermill={"duration": 0.149848, "end_time": "2022-02-01T09:43:12.800733", "exception": false, "start_time": "2022-02-01T09:43:12.650885", "status": "completed"}
# ## Data loss

# %%
total_long_data_loss_df = pd.DataFrame(
    {
        "Stage": [
            "Total unaligned reads",
            "Total aligned reads",
            "Total aligned filtered reads",
            "Total aligned filtered reads<br>within targeted transcripts"
            ],
        "Count": [
            total_unaligned_reads,
            total_aligned_reads,
            total_aligned_filtered_reads,
            total_aligned_filtered_reads_within_targeted_transcripts,
        ]
    }
)
# total_long_data_loss_df

max_total_data_count = total_long_data_loss_df["Count"].max()

def percent_count_relative_to_first_stage(count):
    percent = 100 * count / max_total_data_count
    # percent = round(percent, 1)
    if int(percent) == percent:
        percent = f"{percent:.0f}%"
    else:
        percent = f"{percent:.1f}%"
    return percent

total_long_data_loss_df["%InitialCount"] = total_long_data_loss_df.apply(
    lambda x: percent_count_relative_to_first_stage(x["Count"]),
    axis=1,
)

total_long_data_loss_df

# %%
fig = go.Figure()

x = total_long_data_loss_df["Stage"]
y = total_long_data_loss_df["Count"]
text = total_long_data_loss_df["%InitialCount"]

fig.add_trace(
    go.Bar(
        y=x,
        x=y,
        text=text,
        marker=dict(
            color="white",
            line=dict(
                color="black",
                width=2
            )
        ),
        # marker_color=color_discrete_map[condition],
        orientation="h",
        # marker_pattern_shape="x",
        # marker_pattern_shape="/",
        # name="Data",
        # name=condition,
        # showlegend=False,
        # width=0.3
    )
)

# fig.update_traces(width=0.3)

fig.update_yaxes(
    categoryorder="array",
    categoryarray=list(
        reversed(
            [
                "Total unaligned reads",
                "Total aligned reads",
                "Total aligned filtered reads",
                "Total aligned filtered reads<br>within targeted transcripts"
            ]
        )
    ),
)

fig.update_layout(
    # title="Data loss through the different processing stages",
    width=700,
    height=400,
    # showlegend=False,
    template=template,
    barmode="group",
    # legend_traceorder="reversed",
    # bargap=0.2,  # gap between bars of adjacent location coordinates.
    # bargroupgap=0,  # gap between bars of the same location coordinate. q
)

fig.write_image("Total data loss through the different preprocessing stages - Illumina.svg", width=700, height=400)

fig.show()

# %%
# todo update & uncomment

long_data_loss_df = pd.melt(
    data_loss_df.reset_index().rename(columns={"index": condition_col}),
    id_vars=[condition_col],
    var_name="Stage",
    value_name="Count",
)
long_data_loss_df = long_data_loss_df.sort_values(
    by=[condition_col, "Count"], ascending=False
).reset_index(drop=True)

max_count_per_condition = long_data_loss_df.groupby(condition_col)["Count"].max()


def percent_count_relative_to_first_stage(condition, count):
    max_count = max_count_per_condition[condition]
    percent = 100 * count / max_count
    # percent = round(percent, 1)
    if int(percent) == percent:
        percent = f"{percent:.0f}%"
    else:
        percent = f"{percent:.1f}%"
    return percent


long_data_loss_df["%InitialCount"] = long_data_loss_df.apply(
    lambda x: percent_count_relative_to_first_stage(x[condition_col], x["Count"]),
    axis=1,
)

long_data_loss_df["Stage"] = long_data_loss_df["Stage"].apply(
    lambda x: x.replace(" (", "<br>(") if "(" in x else x
)

long_data_loss_df["StageType"] = long_data_loss_df.apply(
    lambda x: "DataProcessing"
    if x["Stage"] in ["Unique reads", "Distinct unique proteins<br>(mean)"]
    else "Results",
    axis=1,
)


long_data_loss_df


# %%
fig = go.Figure()

# colors = [color_discrete_map[condition] for condition in conditions]

# non_syns_per_max_sol_exp_df = [
#     str(max_sol_exp_df.iloc[:, ML_INPUT_FIRST_COL_POS:].shape[1])
#     for max_sol_exp_df in max_sol_exp_dfs
# ]

for condition in list(reversed(conditions)):
    df = long_data_loss_df.loc[
        (long_data_loss_df[condition_col] == condition) 
        # & (long_data_loss_df["Stage"] != "Pileup reads")
    ]
    preprocessing_df = df.loc[
        ~df["Stage"].isin(["Unique reads", "Distinct unique proteins<br>(mean)"])
    ]
    results_df = df.loc[
        df["Stage"].isin(["Unique reads", "Distinct unique proteins<br>(mean)"])
    ]

    preprocessing_x = preprocessing_df["Stage"]
    preprocessing_y = preprocessing_df["Count"]
    preprocessing_text = preprocessing_df["%InitialCount"]

    fig.add_trace(
        go.Bar(
            y=preprocessing_x,
            x=preprocessing_y,
            text=preprocessing_text,
            marker_color=color_discrete_map[condition],
            orientation="h",
            # marker_pattern_shape="x",
            # marker_pattern_shape="/",
            # name="Data",
            name=condition,
            # showlegend=False,
            # width=0.3
        )
    )

    results_x = results_df["Stage"]
    results_y = results_df["Count"]
    results_text = results_df["%InitialCount"]

    fig.add_trace(
        go.Bar(
            y=results_x,
            x=results_y,
            text=results_text,
            marker_color=color_discrete_map[condition],
            orientation="h",
            # marker_pattern_shape="x",
            marker_pattern_shape="/",
            # name="Data",
            showlegend=False,
            # width=0.3
        )
    )

fig.update_traces(width=0.05)
fig.update_yaxes(
    categoryorder="array",
    categoryarray=list(
        reversed(
            [
                "Unaligned CCS reads",
                "Aligned reads<br>(within ORF)",
                "Pileup reads",
                "Unique reads",
                "Distinct unique proteins<br>(mean)",
            ]
        )
    ),
)

fig.update_layout(
    # title="Data loss through the different processing stages",
    width=1200,
    height=1400,
    # showlegend=False,
    template=template,
    barmode="group",
    legend_traceorder="reversed",
    bargap=0.1,  # gap between bars of adjacent location coordinates.
    bargroupgap=0,  # gap between bars of the same location coordinate. q
)

fig.write_image("Data loss through the different processing stages - Illumina.svg", width=1200, height=1400)

fig.show()

# %%
# todo update & uncomment

# fig = px.bar(
#     long_data_loss_df.sort_values("Count"),
#     y="Stage",
#     x="Count",
#     text="%InitialCount",
#     title="Data loss through the different processing stages",
#     labels={"Count": "", "Stage": ""},
#     color=condition_col,
#     color_discrete_map=color_discrete_map,
#     category_orders=horizontal_category_orders,
#     template=template,
#     barmode="group",
#     orientation="h",
#     width=1000,
#     height=1400,
#     log_x=True
# )

# # line_x = long_data_loss_df["Count"].max() * 1.1
# # line_y0 = "Aligned reads<br>(within ORF)"
# # line_y1 = "Compatible unique edited reads<br>(mean)"
# # fig.add_shape(
# #     type="line",
# #     x0=line_x,
# #     x1=line_x,
# #     y0=line_y0,
# #     y1=line_y1,
# #     line=dict(color="black", width=13),
# #     opacity=0.5
# # )

# # fig.update_traces(bargroupgap=0.01)

# fig.show()


# %% [markdown] papermill={"duration": 0.124528, "end_time": "2022-02-01T09:43:10.054394", "exception": false, "start_time": "2022-02-01T09:43:09.929866", "status": "completed"}
# ## Positions

# %% [markdown] papermill={"duration": 0.149848, "end_time": "2022-02-01T09:43:12.800733", "exception": false, "start_time": "2022-02-01T09:43:12.650885", "status": "completed"}
# ### Correlation matrix

# %%
reads_w_nan_dfs = [reads_df.replace({-1: np.NaN}) for reads_df in reads_dfs]
reads_w_nan_dfs[0]


# %%
# Compute the correlation matrix
corrs = [
    reads_w_nan_df.iloc[:, reads_first_col_pos:].corr()
    for reads_w_nan_df in reads_w_nan_dfs
]

vmin = min(corr.min().min() for corr in corrs)
vmax = max(corr.max().max() for corr in corrs)

# Generate a mask for the upper triangle
masks = [np.triu(np.ones_like(corr, dtype=bool)) for corr in corrs]

# Generate a custom diverging colormap
# cmap = sns.diverging_palette(230, 20, as_cmap=True)


# %%
cmap = sns.diverging_palette(230, 20, as_cmap=True)
# cmap = sns.cubehelix_palette(start=2, rot=0, dark=0, light=.95, reverse=True, as_cmap=True)

for condition, corr, mask in zip(conditions, corrs, masks):

    fig, ax = plt.subplots(figsize=(11, 9))

    # Draw the heatmap with the mask and correct aspect ratio
    sns.heatmap(
        corr,
        mask=mask,
        cmap=cmap,
        square=True,
        linewidths=0.5,
        xticklabels=False,
        yticklabels=False,
        center=0,
        cbar_kws={"shrink": 0.5},
        # vmin=vmin, vmax=vmax,
    )

    title = f"Pearson correlation coefficient between editing sites in {condition} reads"
    
    plt.title(title)
    
    if condition == "PCLO_CHICK":
        fig.savefig(f"{title}.svg", dpi=300);

    # break


# %% [markdown]
# ### Mutual information

# %%
df = reads_w_nan_dfs[0].iloc[:, reads_first_col_pos:]
df.head()


# %%
def calc_normalized_mi(df, pos_1, pos_2):

    two_positions = [pos_1, pos_2]
    complete_rows = df.loc[:, two_positions].notna().all(axis=1)
    df_1_2 = df.loc[complete_rows, two_positions].reset_index(drop=True)

    val_counts = df_1_2.apply(lambda x: x.value_counts())
    val_freqs = val_counts.apply(lambda x: x / x.sum())

    p_1_a = val_freqs.loc[0.0, pos_1]
    p_1_g = val_freqs.loc[1.0, pos_1]
    p_2_a = val_freqs.loc[0.0, pos_2]
    p_2_g = val_freqs.loc[1.0, pos_2]

    df_1_2 = df_1_2.astype(int).astype(str)
    mutual_val_counts = pd.Series(df_1_2[pos_1] + df_1_2[pos_2]).value_counts()
    mutual_val_counts_sum = mutual_val_counts.sum()

    def get_count_or_zero(val):
        try:
            return mutual_val_counts.loc[val]
        except:
            return 0.0

    p_aa = get_count_or_zero("00") / mutual_val_counts_sum
    p_ag = get_count_or_zero("01") / mutual_val_counts_sum
    p_ga = get_count_or_zero("10") / mutual_val_counts_sum
    p_gg = get_count_or_zero("11") / mutual_val_counts_sum

    # ic(p_aa, p_ag_ga, p_gg);

    def calc_partial_mi(common_prob, prob_1, prob_2):
        if common_prob / (prob_1 * prob_2) != 0:
            return common_prob * np.log(common_prob / (prob_1 * prob_2))
        return 0.0

    mi_aa = calc_partial_mi(p_aa, p_1_a, p_2_a)
    mi_ag = calc_partial_mi(p_ag, p_1_a, p_2_g)
    mi_ga = calc_partial_mi(p_ga, p_1_g, p_2_a)
    mi_gg = calc_partial_mi(p_gg, p_1_g, p_2_g)

    mi = mi_aa + mi_ag + mi_ga + mi_gg

    h_1 = -p_1_a * np.log(p_1_a) - p_1_g * np.log(p_1_g)
    h_2 = -p_2_a * np.log(p_2_a) - p_2_g * np.log(p_2_g)

    mi_normalized = mi / min(h_1, h_2)

    return mi_normalized


# %%
mi_dfs = []

for df, condition in zip(reads_w_nan_dfs, conditions):
    df = df.iloc[:, reads_first_col_pos:]

    positions_couples = list(combinations(df.columns, 2))

    with Pool(processes=min(threads, 10)) as pool:
        mis = pool.starmap(
            func=calc_normalized_mi,
            iterable=[(df, pos_1, pos_2) for pos_1, pos_2 in positions_couples],
        )

    mi_df = pd.DataFrame(
        {
            "Pos1": [positions[0] for positions in positions_couples],
            "Pos2": [positions[1] for positions in positions_couples],
            "NormalizedMutualInformation": mis,
        }
    )
    mi_df["Pos1"] = mi_df["Pos1"].astype(int)
    mi_df["Pos2"] = mi_df["Pos2"].astype(int)
    mi_df["Distance"] = np.abs(mi_df["Pos2"] - mi_df["Pos1"])
    mi_df[condition_col] = condition

    mi_dfs.append(mi_df)


merged_mi_df = pd.concat(mi_dfs)
merged_mi_df


# %%
def get_symmetric_mi_df(mi_df, positions):

    mi_df = mi_df.set_index(["Pos1", "Pos2"])

    # positions_couples = list(combinations(positions, 2))
    # positions_couples_set = set(positions_couples)
    positions_couples_set = set(list(combinations(positions, 2)))

    complementary_positions_couples = [
        (int(pos_1), int(pos_2))
        for pos_1 in positions
        for pos_2 in positions
        if (pos_1, pos_2) not in positions_couples_set
    ]

    complementary_mis = []

    for pos_1, pos_2 in complementary_positions_couples:
        if pos_1 == pos_2:
            # altought, technically, the mutual information of a position with itself is 1, is shouldn't be ploted on a heatmap
            mi_normalized = np.nan
        else:
            # get the already-computed mutual information, which appear in reverse order of positions in mi_df
            mi_normalized = mi_df.at[(pos_2, pos_1), "NormalizedMutualInformation"]
        complementary_mis.append(mi_normalized)

    complematary_mis_df = pd.DataFrame(
        data={"NormalizedMutualInformation": complementary_mis},
        index=pd.MultiIndex.from_tuples(
            complementary_positions_couples, names=["Pos1", "Pos2"]
        ),
    )

    complete_mi_df = pd.concat(
        [mi_df.loc[:, ["NormalizedMutualInformation"]], complematary_mis_df]
    ).reset_index()

    symmetric_mi_df = complete_mi_df.pivot(index="Pos1", columns="Pos2")
    symmetric_mi_df = symmetric_mi_df.rename_axis("", axis=0)
    symmetric_mi_df = symmetric_mi_df.set_axis(
        symmetric_mi_df.columns.droplevel(), axis=1
    ).rename_axis("", axis=1)

    return symmetric_mi_df


# %%
conditions

# %%
pclo_index = conditions.index("PCLO_CHICK")

symmetric_pclo_mi_df = get_symmetric_mi_df(
    mi_df=mi_dfs[pclo_index],
    positions=reads_w_nan_dfs[pclo_index].iloc[:, reads_first_col_pos:].columns,
)

# %%
sns.set_theme(style="white")

# Generate a mask for the upper triangle
mask = np.triu(np.ones_like(symmetric_pclo_mi_df, dtype=bool))

# Generate a custom diverging colormap
# cmap = sns.diverging_palette(230, 20, as_cmap=True)
cmap = sns.color_palette("YlOrBr", as_cmap=True)
# cmap = sns.light_palette("seagreen", as_cmap=True)

# vmin = symmetric_pclo_mi_df.min().min()
# vmax = symmetric_pclo_mi_df.max().max()

fig, ax = plt.subplots(figsize=(12, 12))

# Draw the heatmap with the mask and correct aspect ratio
sns.heatmap(
    symmetric_pclo_mi_df,
    mask=mask,
    cmap=cmap,
    square=True,
    linewidths=0.5,
    xticklabels=False,
    yticklabels=False,
    center=0.5,
    cbar_kws={"shrink": 0.5},
    # vmin=0, vmax=1
    # vmin=vmin, vmax=vmax,
)

title = f"Normalized mutual information between editing sites in {conditions[pclo_index]}"

plt.title(title);

fig.savefig(f"{title}.svg", dpi=300);

# %%
fig = px.histogram(
    merged_mi_df,
    x="NormalizedMutualInformation",
    color=condition_col,
    color_discrete_map=color_discrete_map,
    category_orders=category_orders,
    template=template,
    title="Mutual information between each two editing positions",
    log_y=True,
)
# fig.for_each_annotation(
#     lambda a: a.update(text=a.text.replace(f"{condition_col}=", ""))
# )
fig.update_layout(
    barmode="overlay",
    #     # showlegend=False,
)
fig.update_traces(opacity=0.7)
fig.show()

# %%
fig = px.histogram(
    merged_mi_df,
    x="NormalizedMutualInformation",
    y="Distance",
    histfunc="avg",
    # marginal="rug", # or box, violin, rug
    marginal="box",  # or box, violin, rug
    color=condition_col,
    color_discrete_map=color_discrete_map,
    category_orders=category_orders,
    template=template,
    title="Mutual information between each two editing positions",
    # log_y=True
    # log_x=True
)
# fig.for_each_annotation(
#     lambda a: a.update(text=a.text.replace(f"{condition_col}=", ""))
# )
fig.update_layout(
    barmode="overlay",
    width=1000,
    height=450
    #     # showlegend=False,
)
fig.update_traces(opacity=0.7)
fig.show()

# %%
fig = px.scatter(
    merged_mi_df,
    x="Distance",
    y="NormalizedMutualInformation",
    # color_discrete_sequence=[color_discrete_map[conditions[0]]],
    color=condition_col,
    color_discrete_map=color_discrete_map,
    # facet_col=condition_col,
    # facet_col_spacing=facet_col_spacing,
    category_orders=category_orders,
    template=template,
    title="Mutual information between each two editing positions",
    log_y=True,
    log_x=True,
)
# fig.for_each_annotation(
#     lambda a: a.update(text=a.text.replace(f"{condition_col}=", ""))
# )
# fig.update_layout(showlegend=False)
fig.update_traces(
    # opacity=0.7
    marker=dict(opacity=0.7, size=4)
)
fig.show()

# %% [markdown]
# ### Correlation between positions to sites edited in reads

# %%
reads_w_nan_dfs[0]

# %%
edited_positions_per_read_serieses = [
    reads_w_nan_df.loc[:, "EditedPositions"] for reads_w_nan_df in reads_w_nan_dfs
]
edited_positions_per_read_serieses[0]

# %%
editing_in_reads_dfs = [
    reads_w_nan_df.iloc[:, reads_first_col_pos:] for reads_w_nan_df in reads_w_nan_dfs
]
editing_in_reads_dfs[0]

# %%
positions_to_edited_sites_in_reads_correleations = []
for editing_in_reads, edited_positions_per_read in zip(
    editing_in_reads_dfs, edited_positions_per_read_serieses
):
    s = [
        editing_in_reads.iloc[:, x].corr(edited_positions_per_read)
        for x in range(editing_in_reads.shape[1])
    ]
    s = pd.Series(s, index=editing_in_reads.columns)
    positions_to_edited_sites_in_reads_correleations.append(s)
positions_to_edited_sites_in_reads_correleations[0]

# %%
_conditions = []
_conditions_corrs = []
_positions = []

for condition, condition_corrs in zip(
    conditions, positions_to_edited_sites_in_reads_correleations
):
    _conditions.extend([condition] * len(condition_corrs))
    _conditions_corrs.extend(condition_corrs)
    _positions.extend(condition_corrs.index)

positions_to_edited_sites_in_reads_correleations_df = pd.DataFrame(
    {
        condition_col: _conditions,
        "Position": _positions,
        "Pearson": _conditions_corrs,
    }
)

positions_to_edited_sites_in_reads_correleations_df

# %%
for condition in conditions:
    df = positions_to_edited_sites_in_reads_correleations_df.loc[
        positions_to_edited_sites_in_reads_correleations_df[condition_col] == condition
    ]
    fig = px.bar(
        df,
        x="Position",
        y="Pearson",
        # facet_col=condition_col,
        # facet_col_spacing=facet_col_spacing,
        color=condition_col,
        color_discrete_map=color_discrete_map,
        category_orders=category_orders,
        template=template,
        title="Pearson correlation between editing sites to number of sites edited in each read",
    )
    # fig.update_layout(showlegend=False)
    fig.show()

# %%
fig = px.histogram(
    positions_to_edited_sites_in_reads_correleations_df,
    x="Pearson",
    color=condition_col,
    color_discrete_map=color_discrete_map,
    facet_col=condition_col,
    facet_col_spacing=facet_col_spacing,
    facet_col_wrap=facet_col_wrap,
    facet_row_spacing=facet_row_spacing,
    category_orders=category_orders,
    template=template,
    title="Pearson correlation between editing sites to number of sites edited in each read",
)
fig.for_each_annotation(
    lambda a: a.update(text=a.text.replace(f"{condition_col}=", ""))
)
fig.update_layout(showlegend=False)
fig.show()

# %%
_df = positions_to_edited_sites_in_reads_correleations_df

cols = min(facet_col_wrap, len(conditions))
rows = ceil(len(conditions) / cols)
row_col_iter = list(product(range(1, rows + 1), range(1, cols + 1)))[: len(conditions)]

fig = make_subplots(
    rows=rows,
    # cols=cols,
    cols=1,
    # subplot_titles=conditions,
    shared_yaxes=True,
    y_title="Pearson",
)

for condition, (row, col) in zip(conditions, row_col_iter):

    fig.add_trace(
        go.Violin(
            x=_df.loc[_df[condition_col] == condition, condition_col],
            y=_df.loc[_df[condition_col] == condition, "Pearson"],
            line_color=color_discrete_map[condition],
            fillcolor="white",
            box_visible=True,
            meanline_visible=True,
            points="all",
            # pointpos=0,
            # jitter=0.05
        ),
        row=row,
        # col=col
        col=1,
    )

fig.update_layout(
    template=template,
    showlegend=False,
    title_text="Pearson correlation between editing sites to number of sites edited in each read",
    height=250 * rows,
)
fig.update_yaxes(zerolinewidth=zerolinewidth, tickmode="linear", tick0=0, dtick=0.2)

fig.show()

# %% [markdown]
# ### Noise in positions

# %%
noise_dfs = []
for positions_df, condition, strand in zip(positions_dfs, conditions, strands):
    ref_base = "A" if strand == "+" else "T"
    df = positions_df.loc[positions_df["RefBase"] != ref_base]
    df = df.assign(Noise2=df["Noise"] * 100).rename(columns={"Noise2": "%Noise"})
    noise_dfs.append(df)
merged_noise_df = pd.concat(noise_dfs)
merged_noise_df.iloc[[0, 1, -2, -1]]


# %%
fig = px.violin(
    merged_noise_df,
    x=condition_col,
    y="%Noise",
    color=condition_col,
    color_discrete_map=color_discrete_map,
    category_orders=category_orders,
    template=template,
    title="Noise levels",
    labels={"Gene": "Transcript"}
)

fig.update_yaxes(title="% noise", tickmode="linear", tick0=0, dtick=2)
fig.update_layout(showlegend=False, width=1200, height=350)

fig.write_image(
    "Per chrom noise levels - Illumina.svg",
    width=1200,
    height=350,
)

fig.show()


# %%
for positions_df, condition, strand in zip(positions_dfs, conditions, strands):
    ref_base = "A" if strand == "+" else "T"
    df = positions_df.loc[positions_df["RefBase"] != ref_base]
    df = df.assign(Noise2=df["Noise"] * 100).rename(columns={"Noise2": "%Noise"})
    fig = px.bar(
        df,
        x="Position",
        y="%Noise",
        color=condition_col,
        color_discrete_map=color_discrete_map,
        category_orders=category_orders,
        # template=template,
        template="simple_white",
        title=f"Noise per position in {condition}",
    )
    fig.update_layout(showlegend=False)
    fig.update_yaxes(title="% noise")
    fig.show()


# %% [markdown]
# ### Known & new editing sites

# %%
conditions_labels = {condition: ["Edited", "KnownEditing"] for condition in conditions}

conditions_sets = {
    condition: [
        set(positions_df.loc[positions_df[label], "Position"])
        for label in conditions_labels[condition]
    ]
    for positions_df, condition in zip(positions_dfs, conditions)
}

cols = min(facet_col_wrap, len(conditions), 5)
rows = ceil(len(conditions) / cols)

fig, axs = plt.subplots(
    nrows=rows,
    ncols=cols,
    figsize=(3.5 * cols, 2.5 * rows),
    constrained_layout=True,
    gridspec_kw=dict(hspace=0.2, wspace=0.3),
)

for condition, ax in zip(conditions, axs.flat):
    labels = conditions_labels[condition]
    sets = conditions_sets[condition]
    labels[0] = f"Edited\n({len(sets[0])})"
    labels[1] = f"Known editing\n({len(sets[1])})"
    labels = labels[:2]
    sets = sets[:2]
    v_func = venn2
    v_func(sets, set_labels=labels, ax=ax)
    ax.set_title(condition, fontdict=dict(fontsize=16))

# fig.suptitle(title)

if (num_of_axes_to_delete := cols * rows - len(conditions)) > 0:
    cols_to_delete = range(cols - 1, cols - 1 - num_of_axes_to_delete, -1)
    for col_to_delete in cols_to_delete:
        fig.delaxes(axs[rows - 1][col_to_delete])

plt.savefig("Known vs new editing sites - Illumina.svg", format="svg", dpi=300)
        
plt.show()


# %% [markdown]
# ### Editing vs. coverage

# %%
known_sites_df.head()


# %%
ref_base_positions_dfs = []

for positions_df, strand in zip(positions_dfs, strands):
    ref_base = "A" if strand == "+" else "T"
    alt_base = "G" if strand == "+" else "C"
    positions_df = positions_df.loc[positions_df["RefBase"] == ref_base]
    positions_df.insert(
        positions_df.columns.get_loc("EditingFrequency"),
        "EditedReads",
        positions_df[alt_base],
    )
    positions_df = positions_df.drop(
        [
            "MappedBases",
            "Reads",
            "Noise",
            "RefBase",
            "A",
            "T",
            "C",
            "G",
            "InProbRegion",
        ],
        axis=1,
    )
    ref_base_positions_dfs.append(positions_df)

merged_ref_base_positions_df = pd.concat(ref_base_positions_dfs).reset_index(drop=True)

merged_ref_base_positions_df = merged_ref_base_positions_df.merge(
    known_sites_df.loc[:, ["Chrom", "Position", "EditingFrequency"]],
    how="left",
    on=["Chrom", "Position"],
    suffixes=["", "Known"],
)

merged_ref_base_positions_df["Position"] = merged_ref_base_positions_df[
    "Position"
].astype("str")

# merged_ref_base_positions_df = merged_ref_base_positions_df.fillna(0)

merged_ref_base_positions_df.insert(
    merged_ref_base_positions_df.columns.get_loc("EditingFrequency"),
    "%Editing",
    merged_ref_base_positions_df["EditingFrequency"] * 100,
)
merged_ref_base_positions_df.insert(
    merged_ref_base_positions_df.columns.get_loc("EditingFrequencyKnown"),
    "%EditingKnown",
    merged_ref_base_positions_df["EditingFrequencyKnown"] * 100,
)
merged_ref_base_positions_df = merged_ref_base_positions_df.drop(
    ["EditingFrequency", "EditingFrequencyKnown"], axis=1
)

melt_merged_ref_base_positions_df = merged_ref_base_positions_df.melt(
    id_vars=[condition_col, "Chrom", "Position", "Edited", "KnownEditing"],
    # var_name="",
    # value_name=""
)

melt_merged_ref_base_positions_df


# %%
for condition in conditions:
    df = melt_merged_ref_base_positions_df.loc[
        melt_merged_ref_base_positions_df[condition_col] == condition
    ]
    # df = df.loc[~df["variable"].str.contains("%Editing")]
    df = df.loc[df["variable"].isin(["TotalCoverage", "EditedReads"])]
    df = df.loc[df["Edited"]]

    fig = px.bar(
        df,
        x="Position",
        y="value",
        # facet_col="KnownEditing",
        # facet_row="Edited",
        title=(
            f"Editing vs. coverage in {condition}"
            "<br>"
            "<sub>(only edited positions are presented)</sub>"
        ),
        color_discrete_sequence=[
            color_discrete_map[condition],
            color_discrete_map[condition],
        ],
        pattern_shape="variable",
        pattern_shape_map={"TotalCoverage": "", "EditedReads": "\\"},
        # pattern_shape_sequence=["", '\\', ],
        barmode="overlay",
    )
    fig.update_layout(
        legend_title="",
        yaxis_title="Reads",
        template=template,
    )
    fig.show()


# %% [markdown]
# ### Current vs. known editing levels

# %% papermill={"duration": 4.052404, "end_time": "2022-02-01T09:42:53.176715", "exception": false, "start_time": "2022-02-01T09:42:49.124311", "status": "completed"}
# todo retain nan rows and turn nans to 0?

# cols = len(conditions)

cols = min(facet_col_wrap, len(conditions), 4)
rows = ceil(len(conditions) / cols)
row_col_iter = list(product(range(1, rows + 1), range(1, cols + 1)))[: len(conditions)]

fig = make_subplots(
    rows=rows,
    cols=cols,
    subplot_titles=conditions,
    shared_yaxes=True,
    shared_xaxes=True,
    x_title="% editing",
    y_title="% known editing",
    # vertical_spacing=facet_row_spacing
)

df = merged_ref_base_positions_df.copy()
df = df.loc[df["Edited"] & df["KnownEditing"]]

for (row, col), condition, unique_reads_df in zip(
    row_col_iter, conditions, unique_reads_dfs
):
    condition_df = df.loc[df[condition_col] == condition]

    x = condition_df["%Editing"]
    y = condition_df["%EditingKnown"]

    r, pv = scipy.stats.pearsonr(x, y)

    fig.add_trace(
        go.Scatter(
            x=x,
            y=y,
            name=condition,
            mode="markers",
            marker_color=color_discrete_map[condition],
        ),
        row=row,
        col=col,
    )

    fig.add_annotation(
        row=row,
        col=col,
        x=30,
        y=70,
        xref="x",
        yref="y",
        text=f"<b>Pearson's r</b><br>p-val = {pv:.2e}<br>ρ = {r:.2g}",
        bgcolor="white",
        borderpad=4,
        font=dict(size=10),
        opacity=0.8,
        showarrow=False,
    )

fig.update_layout(
    title_text="Correlation between current & previously-reported editing levels",
    showlegend=False,
    template=template,
    height=200 * rows,
)

fig.update_xaxes(range=[0, 100])
fig.update_yaxes(range=[0, 100])

fig.show()


# %%
# todo turn to fig and drop nan (here they're treated as zeros)

for condition in conditions:
    df = melt_merged_ref_base_positions_df.loc[
        melt_merged_ref_base_positions_df[condition_col] == condition
    ]
    df = df.loc[df["variable"].isin(["%Editing", "%EditingKnown"])]
    df["variable"] = df["variable"].apply(
        lambda x: "Current" if x == "%Editing" else "Previous"
    )
    df = df.loc[df["Edited"] | df["KnownEditing"]]
    df["Position"] = df["Position"].astype(int)
    df = df.sort_values("Position")

    fig = px.line(
        df,
        x="Position",
        y="value",
        facet_row="KnownEditing",
        category_orders={"KnownEditing": [True, False]},
        facet_row_spacing=0.1,
        color="variable",
        color_discrete_sequence=subcolors_discrete_map[condition],
        symbol="variable",
        symbol_sequence=["circle", "square"],
        title=f"Current vs. previously-reported editing levels in {condition}",
        template=template,
        height=500,
    )
    fig.update_layout(legend_title="Study")
    fig.update_traces(marker_size=5)
    fig.for_each_yaxis(lambda y: y.update(title="% editing"))  # for facet_row
    fig.update_yaxes(nticks=6)
    # https://community.plotly.com/t/changing-label-of-plotly-express-facet-categories/28066/5
    facet_labels_dict = {
        "KnownEditing=False": "New<br>positions",
        "KnownEditing=True": "Known<br>positions",
    }
    fig.for_each_annotation(lambda a: a.update(text=facet_labels_dict[a.text]))
    fig.show()


# %% [markdown]
# ### Positions' editing levels distribution

# %%
df = merged_ref_base_positions_df.loc[
    merged_ref_base_positions_df["Edited"]
]  # todo or maybe all positions - not just the edited ones?
fig = px.histogram(
    df,
    x="%Editing",
    facet_col=condition_col,
    facet_col_spacing=facet_col_spacing,
    facet_col_wrap=facet_col_wrap,
    facet_row_spacing=facet_row_spacing,
    # histnorm="percent",
    # cumulative=True,
    # marginal="histogram",
    title="Distribution of editing levels in positions",
    labels={"%Editing": "% editing"},
    color=condition_col,
    color_discrete_map=color_discrete_map,
    category_orders=category_orders,
    template=template,
)

# # https://stackoverflow.com/questions/58167028/single-axis-caption-in-plotly-express-facet-plot
# for axis in fig.layout:
#     if type(fig.layout[axis]) == go.layout.YAxis:
#         fig.layout[axis].title.text = ""
# fig.update_layout(showlegend=False, yaxis_title="Positions")
fig.update_layout(showlegend=False)
fig.update_yaxes(title="Positions")
fig.for_each_annotation(lambda a: a.update(text=a.text.replace(f"{condition_col}=", "")))

fig.show()


# %% [markdown] papermill={"duration": 0.030615, "end_time": "2022-02-01T09:42:49.024262", "exception": false, "start_time": "2022-02-01T09:42:48.993647", "status": "completed"}
# ## Num of distinct unique proteins

# %%
fraction_1_gdf = distinct_unique_proteins_df.loc[
    (distinct_unique_proteins_df["Fraction"] == 1.0)
    # & (distinct_unique_proteins_df["Algorithm"] == "Descending")
].groupby(condition_col)

means = (
    fraction_1_gdf["NumOfProteins"]
    .mean()
    .reset_index()
    .rename(columns={"NumOfProteins": "Mean"})
)
stds = (
    fraction_1_gdf["NumOfProteins"]
    .std()
    .reset_index()
    .rename(columns={"NumOfProteins": "STD"})
)

mean_distinct_proteins_df = means.merge(stds, on=condition_col)

mean_distinct_proteins_df

# %%
fig = px.bar(
    mean_distinct_proteins_df,
    x=condition_col,
    y="Mean",
    error_y="STD",
    # points="all",
    labels={"Mean": "Distinct proteins (avg)"},
    color=condition_col,
    color_discrete_map=color_discrete_map,
    category_orders=category_orders,
    template=template,
)

fig.update_yaxes(range=[0, mean_distinct_proteins_df["Mean"].max()])

fig.update_layout(
    showlegend=False, xaxis_title="", width=max(70 * len(conditions), 300), height=400
)

fig.write_image(
    "Distinct proteins (avg) - Illumina.svg",
    width=max(70 * len(conditions), 300),
    height=400
)

fig.show()


# %%
# distinct_proteins_per_mapped_reads_df = distinct_unique_proteins_df.copy()
# distinct_proteins_per_mapped_reads_df = distinct_proteins_per_mapped_reads_df.loc[
#     distinct_proteins_per_mapped_reads_df["Algorithm"] == "Descending"
# ]
# distinct_proteins_per_mapped_reads_df = distinct_proteins_per_mapped_reads_df.loc[
#     distinct_proteins_per_mapped_reads_df["Fraction"] == 1.0
# ]
# distinct_proteins_per_mapped_reads_df["NumOfProteins/MappedReads"] = (
#     distinct_proteins_per_mapped_reads_df["NumOfProteins"]
#     / distinct_proteins_per_mapped_reads_df["NumOfReads"]
# )
# editable_aas_per_sample_dict = {
#     condition: editable_aas
#     for condition, editable_aas in zip(conditions, editable_aas_per_sample)
# }
# distinct_proteins_per_mapped_reads_df["EditableAAs"] = [
#     editable_aas_per_sample_dict[condition]
#     for condition in distinct_proteins_per_mapped_reads_df[condition_col]
# ]
# # df["NumOfProteins/Read/EditableAmicoAcids"] = df["NumOfProteins"] / df["NumOfReads"] / df["EditableAminoAcidsPerSample"]
# distinct_proteins_per_mapped_reads_df


# %%
# mean_distinct_proteins_per_mapped_reads_df = pd.DataFrame(
#     {
#         condition_col: conditions,
#         "EditableAAS": editable_aas_per_sample,
#     }
# )

# fraction_1_gdf = distinct_proteins_per_mapped_reads_df.loc[
#     distinct_proteins_per_mapped_reads_df["Fraction"] == 1.0
# ].groupby(condition_col)

# means = (
#     fraction_1_gdf["NumOfProteins/MappedReads"]
#     .mean()
#     .reset_index()
#     .rename(columns={"NumOfProteins/MappedReads": "Mean"})
# )
# stds = (
#     fraction_1_gdf["NumOfProteins/MappedReads"]
#     .std()
#     .reset_index()
#     .rename(columns={"NumOfProteins/MappedReads": "SD"})
# )

# mean_distinct_proteins_per_mapped_reads_df = (
#     mean_distinct_proteins_per_mapped_reads_df.merge(means).merge(stds)
# )

# mean_distinct_proteins_per_mapped_reads_df

# %%
# fig = go.Figure()

# for condition in conditions:
#     condition_df = distinct_proteins_per_mapped_reads_df.loc[
#         distinct_proteins_per_mapped_reads_df[condition_col] == condition
#     ]

#     x = condition_df["EditableAAs"]
#     y = condition_df["NumOfProteins/MappedReads"]

#     fig.add_trace(
#         go.Scatter(
#             x=x,
#             y=y,
#             name=condition,
#             mode="markers",
#             marker_color=color_discrete_map[condition],
#         )
#     )

# # correlate all x and y values
# x = distinct_proteins_per_mapped_reads_df["EditableAAs"]
# y = distinct_proteins_per_mapped_reads_df["NumOfProteins/MappedReads"]
# r, pv = scipy.stats.pearsonr(x, y)

# fig.add_annotation(
#     x=13,
#     y=0.09,
#     xref="x",
#     yref="y",
#     text=f"<b>Pearson's r</b><br>p-val = {pv:.2e}<br>ρ = {r:.2g}",
#     bgcolor="white",
#     borderpad=4,
#     font=dict(size=12),
#     opacity=0.8,
#     showarrow=False,
# )

# fig.update_xaxes(
#     range=[0, distinct_proteins_per_mapped_reads_df["EditableAAs"].max() * 1.1]
# )
# fig.update_yaxes(
#     range=[
#         0,
#         distinct_proteins_per_mapped_reads_df["NumOfProteins/MappedReads"].max() * 1.1,
#     ]
# )

# fig.update_layout(
#     # showlegend=False,
#     template=template,
#     xaxis_title="Editable amino acids",
#     yaxis_title="Distinct unique proteins /<br>mapped reads",
#     width=800,
#     height=600,
# )

# fig.show()


# %%
# fig = px.bar(
#     mean_distinct_proteins_per_mapped_reads_df,
#     x=condition_col,
#     y="Mean",
#     error_y="SD",
#     # points="all",
#     labels={"Mean": "Mean distinct unique proteins /<br>mapped reads"},
#     color=condition_col,
#     color_discrete_map=color_discrete_map,
#     category_orders=category_orders,
#     template=template,
# )

# fig.update_yaxes(range=[0, mean_distinct_proteins_per_mapped_reads_df["Mean"].max()])

# fig.update_layout(
#     showlegend=False, xaxis_title="", width=max(70 * len(conditions), 300), height=400
# )

# fig.show()


# %%
# fig = px.violin(
#     distinct_proteins_per_mapped_reads_df,
#     x=condition_col,
#     y="NumOfProteins/MappedReads",
#     # points="all",
#     labels={"NumOfProteins/MappedReads": "Distinct unique proteins /<br>mapped reads"},
#     color=condition_col,
#     color_discrete_map=color_discrete_map,
#     category_orders=category_orders,
#     template=template,
# )

# fig.update_yaxes(
#     range=[0, distinct_proteins_per_mapped_reads_df["NumOfProteins/MappedReads"].max()]
# )

# fig.update_layout(
#     showlegend=False, xaxis_title="", width=max(70 * len(conditions), 300), height=400
# )

# fig.show()


# %%
# distinct_unique_proteins_df["NumOfProteins"].max()


# %%
# condition = conditions[0]
# df = distinct_unique_proteins_df.loc[
#     distinct_unique_proteins_df[condition_col] == condition
# ]
# grouped_df = df.groupby("Fraction")
# x_fraction_mean = grouped_df["NumOfReads"].mean().reset_index()
# x_fraction_mean


# %%
# y_fraction_mean = grouped_df["NumOfProteins"].mean().reset_index()
# y_fraction_mean


# %%
# x_fraction_mean.merge(y_fraction_mean, on="Fraction")


# %%
# x_axis_name = "Reads"
# y_axis_name = "Distinct unique proteins"
# head_title = (
#     "Distinct unique proteins vs. sequencing depth"
#     # "<br>"
#     # # f"<sub>({alg_repetitions * 2} repetitions over each fraction of data)</sub>"
#     # "<sub>(100 repetitions over each fraction of data)</sub>"
# )
# _marker_size = 7
# maximal_x = 0

# # Initialize figure with subplots
# fig = make_subplots(
#     rows=1, cols=1, print_grid=False, x_title=x_axis_name, y_title=y_axis_name
# )

# first_data_trace = True

# # Add traces
# for col, condition in enumerate(conditions, start=1):

#     df = distinct_unique_proteins_df.loc[
#         distinct_unique_proteins_df[condition_col] == condition
#     ]
#     x_measured = df["NumOfReads"]
#     y_measured = df["NumOfProteins"]

#     if first_data_trace:
#         fig.add_trace(
#             go.Scatter(
#                 x=x_measured,
#                 y=y_measured,
#                 mode="markers",
#                 marker=dict(
#                     color=subcolors_discrete_map[condition][0], size=_marker_size
#                 ),
#                 legendgroup="Partial-CDS, Illumina",  # this can be any string
#                 legendgrouptitle_text="Partial-CDS, Illumina",
#                 # name="Measured",
#                 name=condition,
#             ),
#         )
#         first_data_trace = False
#     else:
#         fig.add_trace(
#             go.Scatter(
#                 x=x_measured,
#                 y=y_measured,
#                 mode="markers",
#                 marker=dict(
#                     color=subcolors_discrete_map[condition][0], size=_marker_size
#                 ),
#                 legendgroup="Partial-CDS, Illumina",  # this can be any string
#                 name=condition,
#             ),
#         )

#     grouped_df = df.groupby("Fraction")
#     x_fraction_mean = grouped_df["NumOfReads"].mean().reset_index()
#     y_fraction_mean = grouped_df["NumOfProteins"].mean().reset_index()
#     mean_fraction_df = x_fraction_mean.merge(y_fraction_mean, on="Fraction")

#     fig.add_trace(
#         go.Scatter(
#             x=mean_fraction_df["NumOfReads"],
#             y=mean_fraction_df["NumOfProteins"],
#             mode="lines",
#             marker=dict(color=subcolors_discrete_map[condition][0], size=_marker_size),
#             showlegend=False,
#         ),
#     )

#     maximal_x = max(maximal_x, x_measured.max())

# dscam_ys = [
#     19_008,
#     18_496,
# ]
# dscam_legend_names = [
#     "Theoretical maximum",
#     "Measured",
# ]
# # dscam_legend_names = ["measured", "theoretical maximum"]
# dscam_colors = ["grey", "black"]
# fig.add_trace(
#     go.Scatter(
#         x=[0.05 * maximal_x, 1.05 * maximal_x],
#         y=[dscam_ys[0], dscam_ys[0]],
#         mode="lines",
#         line=dict(
#             color=dscam_colors[0],
#             dash="dash",
#             # width=3
#         ),
#         legendgroup="DSCAM",  # this can be any string
#         legendgrouptitle_text="DSCAM",
#         name=dscam_legend_names[0],
#         # name=f"DSCAM {dscam_legend_names[1]}",
#     ),
# )
# fig.add_trace(
#     go.Scatter(
#         x=[0.05 * maximal_x, 1.05 * maximal_x],
#         y=[dscam_ys[1], dscam_ys[1]],
#         mode="lines",
#         line=dict(
#             color=dscam_colors[1],
#             dash="dash",
#             # width=3
#         ),
#         legendgroup="DSCAM",  # this can be any string
#         name=dscam_legend_names[1],
#         # name=f"DSCAM {dscam_legend_names[0]}",
#     ),
# )

# # fig.update_xaxes(type="log")
# fig.update_yaxes(type="log")

# fig.update_layout(
#     title_text=head_title,
#     # legend_title_text=condition_col,
#     # legend_title_text="Partial-CDS, Illumina",
#     template=template,
#     # legend_font=dict(size=8),
#     # legend_grouptitlefont=dict(size=8),
#     # legend_tracegroupgap=4,
#     width=900,
#     height=650,
# )
# # fig.write_image("Distinct unique proteins vs. sequencing depth - Illumina.png", format='png',engine='kaleido')
# fig.write_image(
#     "Distinct unique proteins vs. sequencing depth - Illumina.svg",
#     width=650,
#     height=500,
# )
# fig.show()


# %%
x_axis_name = "Mapped reads"
y_axis_name = "Distinct proteins"
head_title = (
    "Distinct proteins vs. sequencing depth"
    # "<br>"
    # # f"<sub>({alg_repetitions * 2} repetitions over each fraction of data)</sub>"
    # "<sub>(100 repetitions over each fraction of data)</sub>"
)
_marker_size = 7
maximal_x = 0

# Initialize figure with subplots
fig = make_subplots(
    rows=1, cols=1, print_grid=False, x_title=x_axis_name, y_title=y_axis_name
)

first_data_trace = True

# Add traces
for col, condition in enumerate(conditions, start=1):

    df = distinct_unique_proteins_df.loc[
        distinct_unique_proteins_df[condition_col] == condition
    ]
    x_measured = df["NumOfReads"]
    y_measured = df["NumOfProteins"]

    if first_data_trace:
        fig.add_trace(
            go.Scatter(
                x=x_measured,
                y=y_measured,
                mode="markers",
                marker=dict(
                    color=subcolors_discrete_map[condition][0], size=_marker_size
                ),
                legendgroup="Partial-CDS, Illumina",  # this can be any string
                legendgrouptitle_text="Partial-CDS, Illumina",
                # name="Measured",
                name=condition,
            ),
        )
        first_data_trace = False
    else:
        fig.add_trace(
            go.Scatter(
                x=x_measured,
                y=y_measured,
                mode="markers",
                marker=dict(
                    color=subcolors_discrete_map[condition][0], size=_marker_size
                ),
                legendgroup="Partial-CDS, Illumina",  # this can be any string
                name=condition,
            ),
        )

    grouped_df = df.groupby("Fraction")
    x_fraction_mean = grouped_df["NumOfReads"].mean().reset_index()
    y_fraction_mean = grouped_df["NumOfProteins"].mean().reset_index()
    mean_fraction_df = x_fraction_mean.merge(y_fraction_mean, on="Fraction")

    fig.add_trace(
        go.Scatter(
            x=mean_fraction_df["NumOfReads"],
            y=mean_fraction_df["NumOfProteins"],
            mode="lines",
            marker=dict(color=subcolors_discrete_map[condition][0], size=_marker_size),
            showlegend=False,
        ),
    )

    maximal_x = max(maximal_x, x_measured.max())

# dscam_ys = [
#     36_016,
#     18_496,
# ]
# dscam_legend_names = [
#     "Theoretical maximum",
#     "Measured",
# ]
# # dscam_legend_names = ["measured", "theoretical maximum"]
# dscam_colors = ["grey", "black"]
# fig.add_trace(
#     go.Scatter(
#         x=[0.05 * maximal_x, 1.05 * maximal_x],
#         y=[dscam_ys[0], dscam_ys[0]],
#         mode="lines",
#         line=dict(
#             color=dscam_colors[0],
#             dash="dash",
#             # width=3
#         ),
#         legendgroup="DSCAM",  # this can be any string
#         legendgrouptitle_text="DSCAM",
#         name=dscam_legend_names[0],
#         # name=f"DSCAM {dscam_legend_names[1]}",
#     ),
# )
# fig.add_trace(
#     go.Scatter(
#         x=[0.05 * maximal_x, 1.05 * maximal_x],
#         y=[dscam_ys[1], dscam_ys[1]],
#         mode="lines",
#         line=dict(
#             color=dscam_colors[1],
#             dash="dash",
#             # width=3
#         ),
#         legendgroup="DSCAM",  # this can be any string
#         name=dscam_legend_names[1],
#         # name=f"DSCAM {dscam_legend_names[0]}",
#     ),
# )

dscam_dashed_lined_width = 3.5
dscam_ys = [36_016, 18_496, ]
dscam_legend_names = ["Theoretical maximum", "Measured", ]
# dscam_legend_names = ["measured", "theoretical maximum"]
dscam_colors = ["grey", "black"]
# dscam_colors = ["black", "grey"]
# dscam_dashes = ["dash", "dot"]
dscam_dashes = ["dash", "dash"]
fig.add_trace(
    go.Scatter(
        x=[0.05 * maximal_x, 1.05 * maximal_x],
        y=[dscam_ys[0], dscam_ys[0]],
        mode="lines",
        line=dict(
            color=dscam_colors[0],
            dash=dscam_dashes[0],
            width=dscam_dashed_lined_width
        ),
        opacity=0.6,
        legendgroup="DSCAM",  # this can be any string
        legendgrouptitle_text="DSCAM",
        name=dscam_legend_names[0],
        # name=f"DSCAM {dscam_legend_names[1]}",
    ),
)
fig.add_trace(
    go.Scatter(
        x=[0.05 * maximal_x, 1.05 * maximal_x],
        y=[dscam_ys[1], dscam_ys[1]],
        mode="lines",
        line=dict(
            color=dscam_colors[1],
            dash=dscam_dashes[1],
            width=dscam_dashed_lined_width
        ),
        opacity=0.6,
        legendgroup="DSCAM",  # this can be any string
        name=dscam_legend_names[1],
        # name=f"DSCAM {dscam_legend_names[0]}",
    ),
)

# fig.update_xaxes(type="log")
# fig.update_yaxes(type="log")

fig.update_layout(
    # title_text=head_title,
    # legend_title_text=condition_col,
    # legend_title_text="Partial-CDS, Illumina",
    template=template,
    legend_font=dict(size=12),
    legend_grouptitlefont=dict(size=14),
    # legend_font=dict(size=14),
    # legend_grouptitlefont=dict(size=16),
    # legend_tracegroupgap=4,
    width=900,
    height=650,
    # width=650,
    # height=500,
)
# fig.write_image("Distinct unique proteins vs. sequencing depth - Illumina.png", format='png',engine='kaleido')
fig.write_image(
    "Distinct proteins vs. sequencing depth - Illumina.svg",
    # width=1100,
    # height=1000,
    width=900,
    height=650,
)
fig.show()


# %% [markdown]
# ### NaNs distribution

# %%
expanded_distinct_unique_proteins_df

# %%
ambigous_positions_df = (
    # expanded_distinct_unique_proteins_df.loc[
    #     expanded_distinct_unique_proteins_df["Algorithm"] == "Descending"
    # ]
    expanded_distinct_unique_proteins_df
    .groupby(
        [
            condition_col,
            "Fraction",
            "FractionRepetition",
            "Algorithm",
            "AlgorithmRepetition",
        ]
    )
    .agg(
        # NumOfSupportingReads=("NumOfReads", sum),
        SumOfAmbigousPositions=("AmbigousPositions", sum),
        MeanOfAmbigousPositions=("AmbigousPositions", np.mean),
    )
    .reset_index()
    # .rename(columns={"NumOfReads": "NumOfSupportingReads"})
    # .merge(distinct_unique_proteins_df, on=[condition_col, "Fraction", "FractionRepetition", "Algorithm", "AlgorithmRepetition"])
    # .assign(
    #     SupportingReadsPerProtein = lambda x: x["NumOfSupportingReads"] / x["NumOfProteins"],
    #     PercentSupportedReads = lambda x: 100 * x["NumOfSupportingReads"] / x["NumOfReads"]
    # )
    # .rename(columns={"PercentSupportedReads": "%SupportedReads"})
)
ambigous_positions_df

# %%
fig = px.histogram(
    ambigous_positions_df.loc[ambigous_positions_df["Fraction"] == 1.0],
    x="MeanOfAmbigousPositions",
    facet_col=condition_col,
    facet_col_wrap=facet_col_wrap,
    facet_col_spacing=facet_col_spacing,
    facet_row_spacing=facet_row_spacing,
    # labels={"MeanOfAmbigousPositions": "Mean ambigous positions in proteins of a solution"},
    title="Distribution of ambiguous positions in different solutions",
    color=condition_col,
    color_discrete_map=color_discrete_map,
    category_orders=category_orders,
    template=template,
)

# https://stackoverflow.com/questions/58167028/single-axis-caption-in-plotly-express-facet-plot
for axis in fig.layout:
    if type(fig.layout[axis]) == go.layout.YAxis:
        fig.layout[axis].title.text = ""
fig.update_layout(showlegend=False, yaxis_title="# solutions")
# fig.for_each_yaxis(lambda y: y.update(title="Transcripts"))   # for facet_row

fig.show()


# %%
df1 = (
    expanded_distinct_unique_proteins_df.loc[
        (expanded_distinct_unique_proteins_df["Fraction"] == 1.0)
        # & (expanded_distinct_unique_proteins_df["Algorithm"] == "Descending")
    ]
    .groupby([condition_col, "Protein"])
    .agg(NumOfSolutions=("Protein", "size"))
    .reset_index()
    .rename(columns={"NumOfSolutions": "#SolutionIncluded"})
)

df1

# %%
dfs2 = []

for condition, unique_proteins_df in zip(conditions, unique_proteins_dfs):
    df = df1.loc[df1[condition_col] == condition]
    df2 = (
        unique_proteins_df.loc[:, [condition_col, "Protein", "AmbigousPositions"]]
        .merge(df.loc[:, ["Protein", "#SolutionIncluded"]], how="left", on="Protein")
        .fillna(0)
    )
    df2["#SolutionIncluded"] = df2["#SolutionIncluded"].astype(int)
    dfs2.append(df2)
# df2

dfs2[0]

# %%
df3 = pd.concat(dfs2)
df3

# %%
cols = min(facet_col_wrap, len(conditions), 5)
rows = ceil(len(conditions) / cols)
row_col_iter = list(product(range(1, rows + 1), range(1, cols + 1)))[: len(conditions)]

x_title = "Ambiguous positions in a protein"
y_title = "Num solutions a protein is included at (avg)"
# title_text = "Distribution of min & max estimates of non-syn mutations per protein"

fig = make_subplots(
    rows=rows,
    cols=cols,
    subplot_titles=conditions,
    shared_yaxes=True,
    shared_xaxes=True,
    x_title=x_title,
    y_title=y_title,
)

# plot averaged histograms
for (
    (row, col),
    condition,
) in zip(row_col_iter, conditions):
    _df3 = df3.loc[df3[condition_col] == condition]
    x = _df3["AmbigousPositions"]
    y = _df3["#SolutionIncluded"]
    fig.add_trace(
        go.Histogram(
            x=x,
            y=y,
            histfunc="avg",
            marker_color=color_discrete_map[condition],
            name=condition,
        ),
        row=row,
        col=col,
    )

# find max_x and max_y
f = fig.full_figure_for_development(warn=False)
data_traces = {}
for condition in conditions:
    for data in f.data:
        if data.name == condition:
            data_traces[condition] = data
            continue
max_x = max([max(data_traces[condition].x) for condition in conditions])
max_y = max([max(data_traces[condition].y) for condition in conditions])

fig.update_layout(
    template=template,
    barmode="overlay",  # Overlay both histograms
    # title_text=title_text,
    # title_y=0.95,
    showlegend=False,
    height=max(170 * rows, 300),
    width=max(220 * cols, 500),
)

fig.update_xaxes(range=[0, max_x * 1.1])
fig.update_yaxes(range=[0, max_y * 1.1])

fig.show()


# %% [markdown]
# ### Comparing algorithms

# %% [markdown]
# #### Coupled diff comparisons

# %%
distinct_unique_proteins_df

# %%
max_distinct_proteins_per_transcript_and_alg_df = distinct_unique_proteins_df.loc[distinct_unique_proteins_df["Fraction"] == 1.0].copy()

max_distinct_proteins_per_transcript_and_alg_df[
    "MaxNumOfProteins"
] = max_distinct_proteins_per_transcript_and_alg_df.groupby(
    [condition_col, "Algorithm"]
)[
    "NumOfProteins"
].transform(
    max
)
max_distinct_proteins_per_transcript_and_alg_df["IsMaxNumOfProteins"] = (
    max_distinct_proteins_per_transcript_and_alg_df["NumOfProteins"]
    == max_distinct_proteins_per_transcript_and_alg_df["MaxNumOfProteins"]
)

max_distinct_proteins_per_transcript_and_alg_df = (
    max_distinct_proteins_per_transcript_and_alg_df.loc[
        max_distinct_proteins_per_transcript_and_alg_df["IsMaxNumOfProteins"]
    ]
)
max_distinct_proteins_per_transcript_and_alg_df = (
    max_distinct_proteins_per_transcript_and_alg_df.drop_duplicates(
        subset=[condition_col, "Algorithm"], ignore_index=True
    )
)

max_distinct_proteins_per_transcript_and_alg_df

# %%
max_distinct_proteins_per_transcript_and_alg_df.loc[max_distinct_proteins_per_transcript_and_alg_df[condition_col] == "ROBO2_HUMAN"]

# %%
25228 - 25178

# %%
asc_df = max_distinct_proteins_per_transcript_and_alg_df.loc[max_distinct_proteins_per_transcript_and_alg_df["Algorithm"] == "Ascending"].reset_index(drop=True)
desc_df = max_distinct_proteins_per_transcript_and_alg_df.loc[max_distinct_proteins_per_transcript_and_alg_df["Algorithm"] != "Ascending"].reset_index(drop=True)

ic(asc_df["NumOfProteins"].eq(desc_df["NumOfProteins"]).sum()) # ==
ic(asc_df["NumOfProteins"].gt(desc_df["NumOfProteins"]).sum())  # >
ic(asc_df["NumOfProteins"].lt(desc_df["NumOfProteins"]).sum());  # <

# %%
greater_asc_transcripts = asc_df.loc[asc_df["NumOfProteins"].gt(desc_df["NumOfProteins"]), condition_col]
greater_asc_transcripts

# %%
distinct_unique_proteins_df.loc[
    (distinct_unique_proteins_df["Fraction"] == 1.0) &
    (distinct_unique_proteins_df[condition_col].isin(greater_asc_transcripts))
].groupby([condition_col, "Algorithm"])["NumOfProteins"].value_counts()

# %%
max_distinct_proteins_per_transcript_and_alg_df.loc[max_distinct_proteins_per_transcript_and_alg_df[condition_col].isin(greater_asc_transcripts)]

# %%
max_distinct_proteins_per_transcript_and_alg_and_fracrepetition_df = (
    distinct_unique_proteins_df.copy()
)

max_distinct_proteins_per_transcript_and_alg_and_fracrepetition_df[
    "MaxNumOfProteins"
] = max_distinct_proteins_per_transcript_and_alg_and_fracrepetition_df.groupby(
    [condition_col, "Fraction", "Algorithm", "FractionRepetition"]
)[
    "NumOfProteins"
].transform(
    max
)
max_distinct_proteins_per_transcript_and_alg_and_fracrepetition_df[
    "IsMaxNumOfProteins"
] = (
    max_distinct_proteins_per_transcript_and_alg_and_fracrepetition_df["NumOfProteins"]
    == max_distinct_proteins_per_transcript_and_alg_and_fracrepetition_df[
        "MaxNumOfProteins"
    ]
)

max_distinct_proteins_per_transcript_and_alg_and_fracrepetition_df = (
    max_distinct_proteins_per_transcript_and_alg_and_fracrepetition_df.loc[
        max_distinct_proteins_per_transcript_and_alg_and_fracrepetition_df[
            "IsMaxNumOfProteins"
        ]
    ]
).reset_index(drop=True)
# max_distinct_proteins_per_transcript_and_alg_and_fracrepetition_df["Duplicated"] = max_distinct_proteins_per_transcript_and_alg_and_fracrepetition_df.duplicated(subset=[condition_col, "Fraction", "Algorithm", "FractionRepetition"])
max_distinct_proteins_per_transcript_and_alg_and_fracrepetition_df = (
    max_distinct_proteins_per_transcript_and_alg_and_fracrepetition_df.drop_duplicates(
        subset=[condition_col, "Fraction", "Algorithm", "FractionRepetition"],
        ignore_index=True,
    )
)
max_distinct_proteins_per_transcript_and_alg_and_fracrepetition_df = (
    max_distinct_proteins_per_transcript_and_alg_and_fracrepetition_df.sort_values(
        [condition_col, "Fraction", "FractionRepetition", "Algorithm"],
        ignore_index=True,
    )
)

max_distinct_proteins_per_transcript_and_alg_and_fracrepetition_df

# %%
asc_df = max_distinct_proteins_per_transcript_and_alg_and_fracrepetition_df.loc[
    max_distinct_proteins_per_transcript_and_alg_and_fracrepetition_df["Algorithm"]
    == "Ascending"
].reset_index(drop=True)
desc_df = max_distinct_proteins_per_transcript_and_alg_and_fracrepetition_df.loc[
    max_distinct_proteins_per_transcript_and_alg_and_fracrepetition_df["Algorithm"]
    != "Ascending"
].reset_index(drop=True)

assert len(asc_df) == len(desc_df)

ic(asc_df["NumOfProteins"].eq(desc_df["NumOfProteins"]).sum())  # ==
ic(asc_df["NumOfProteins"].gt(desc_df["NumOfProteins"]).sum())  # >
ic(asc_df["NumOfProteins"].lt(desc_df["NumOfProteins"]).sum()); # <

# %%
greater_asc_transcripts = asc_df.loc[
    asc_df["NumOfProteins"].gt(desc_df["NumOfProteins"]), condition_col
].unique()

ic(len(greater_asc_transcripts))

greater_asc_transcripts

# %%
df = max_distinct_proteins_per_transcript_and_alg_and_fracrepetition_df.loc[
    max_distinct_proteins_per_transcript_and_alg_and_fracrepetition_df[
        condition_col
    ].isin(greater_asc_transcripts)
]
df = df.drop(
    ["AlgorithmRepetition", "IsMaxNumOfProteins", "MaxNumOfProteins", "Proteins"],
    axis=1,
).pivot(
    index=[condition_col, "NumOfReads", "Fraction", "FractionRepetition"],
    columns="Algorithm",
)
df = df.set_axis(df.columns.get_level_values(1).values, axis=1)
df = df.reset_index()
df["Desc - Asc"] = df["Descending"] - df["Ascending"]
df = df.loc[df["Desc - Asc"] < 0].reset_index(drop=True)
df

# %%
# fig = px.histogram(
#     df,
#     x="Desc - Asc",
#     template=template,
#     # text_auto=True,
#     facet_col="Fraction",
#     facet_row=condition_col,
#     # facet_col_wrap=facet_col_wrap,
#     # color="Fraction",
#     color=condition_col,
#     color_discrete_map=color_discrete_map,
#     # category_orders={"Fraction": [0.2, 0.4, 0.6, 0.8, 1.0]},
#     # barmode="group"
# )

# # Reduce opacity to see both histograms
# # fig.update_traces(opacity=0.75)

# fig.update_layout(
#     showlegend=False,
#     # height=800,
#     # width=500,
#     width=800,
#     # # height=300
#     # # width=800,
#     height=1500
#     # barmode='overlay' # Overlay both histograms
# )



# fig.show()

# %%
fig = px.histogram(
    df,
    x="Desc - Asc",
    template=template,
    # text_auto=True,
    facet_col="Fraction",
    # facet_row="Fraction",
    # color="Fraction",
    color=condition_col,
    color_discrete_map=color_discrete_map,
    category_orders={"Fraction": [0.2, 0.4, 0.6, 0.8, 1.0]},
    # barmode="group"
)

# Reduce opacity to see both histograms
fig.update_traces(opacity=0.75)

fig.update_layout(
    # showlegend=False,
    # height=800,
    # width=500,
    width=1350,
    # # height=300
    # # width=800,
    height=430
    # barmode='overlay' # Overlay both histograms
)



fig.show()

# %%

# %%

# %%

# %%
# asc_df = distinct_unique_proteins_df.loc[
#     distinct_unique_proteins_df["Algorithm"] == "Ascending"
# ].reset_index(drop=True)
# desc_df = distinct_unique_proteins_df.loc[
#     distinct_unique_proteins_df["Algorithm"] == "Descending"
# ].reset_index(drop=True)
# desc_df

# %%
# desc_df.groupby([condition_col, "Fraction", "FractionRepetition"])[
#     "NumOfProteins"
# ].mean().reset_index()

# %%
# desc_df["NumOfProteins"]

# %%
# algs_diff_df = desc_df.loc[
#     :, [condition_col, "Fraction", "FractionRepetition", "AlgorithmRepetition"]
# ]
# algs_diff_df["Desc - Asc"] = desc_df["NumOfProteins"] - asc_df["NumOfProteins"]
# algs_diff_df["Diff"] = algs_diff_df["Desc - Asc"].apply(
#     lambda x: "> 0" if x > 0 else "= 0" if x == 0 else "< 0"
# )

# # algs_diff_df = algs_diff_df.loc[algs_diff_df["Fraction"] == 1.0].reset_index(drop=True)

# algs_diff_df

# %%
# algs_diff_df.loc[algs_diff_df["Diff"] == "> 0"]

# %%
# algs_diff_df.loc[algs_diff_df["Diff"] == "< 0"]

# %%
# df = algs_diff_df.rename(columns={"Desc - Asc": "Desc <span>&#8722;</span>  Asc"})

# fig = px.scatter(
#     df,
#     x="Fraction",
#     y="Desc <span>&#8722;</span>  Asc",
#     color=condition_col,
#     color_discrete_map=color_discrete_map,
#     facet_row=condition_col,
#     facet_row_spacing=facet_col_spacing,
#     facet_col="Diff",
#     facet_col_spacing=facet_col_spacing,
#     symbol="Diff",
#     category_orders=category_orders,
#     template=template,
#     title="Differences between coupled solutions' sizes",
#     # shared_
# )
# fig.update_layout(showlegend=False)
# # https://www.toptal.com/designers/htmlarrows/math/minus-sign/
# # fig.update_yaxes(title="Desc <span>&#8722;</span>  Asc")
# # https://community.plotly.com/t/changing-label-of-plotly-express-facet-categories/28066/5?u=kaparanewbie
# fig.for_each_annotation(lambda a: a.update(text=a.text.replace("Diff=", "Diff ")))
# fig.for_each_annotation(lambda a: a.update(text=a.text.replace(f"{condition_col}=", "")))

# # fig.update_xaxes(
# #     tickmode = 'array',
# #     tickvals = list(range(1, 11)),
# #     ticktext = [str(0.1 + 0.1 * x)[:3] for x in range(10)]
# # )

# # fig.update_xaxes(visible=True)
# # fig.update_xaxes(zeroline=True)
# # fig.update_xaxes(zerolinewidth=4)
# fig.update_yaxes(zerolinewidth=zerolinewidth)


# fig.show()

# %%
# df = algs_diff_df.rename(columns={"Desc - Asc": "Desc <span>&#8722;</span>  Asc"})

# fig = px.scatter(
#     df,
#     x="Fraction",
#     y="Desc <span>&#8722;</span>  Asc",
#     color=condition_col,
#     color_discrete_map=color_discrete_map,
#     # facet_row=condition_col,
#     # facet_row_spacing=facet_col_spacing*2,
#     facet_col="Diff",
#     facet_col_spacing=facet_col_spacing,
#     symbol="Diff",
#     category_orders=category_orders,
#     template=template,
#     title="Differences between coupled solutions' sizes",
#     # shared_
# )
# # fig.update_layout(showlegend=False)
# # https://www.toptal.com/designers/htmlarrows/math/minus-sign/
# # fig.update_yaxes(title="Desc <span>&#8722;</span>  Asc")
# # https://community.plotly.com/t/changing-label-of-plotly-express-facet-categories/28066/5?u=kaparanewbie
# # fig.for_each_annotation(lambda a: a.update(text=a.text.replace("Diff=", "Diff ")))

# # fig.update_xaxes(
# #     tickmode = 'array',
# #     tickvals = list(range(1, 11)),
# #     ticktext = [str(0.1 + 0.1 * x)[:3] for x in range(10)]
# # )

# fig.update_yaxes(zerolinewidth=zerolinewidth)

# fig.show()

# %%
# x_axis_name = "Fraction"
# y_axis_name = "#Desc <span>&#8722;</span> #Asc"
# head_title = "Differences between coupled solutions' sizes"

# # basic_diff_names = ["Diff < 0", "Diff = 0", "Diff > 0"]
# # subplot_titles = []
# # first_row = True
# # for condition in conditions:
# #     for diff in basic_diff_names:
# #         st = f"<sub>{condition}</sub>"
# #         if first_row:
# #             st = f"{diff}<br>" + st
# #         subplot_titles.append(st)
# #     first_row = False

# # rows = len(conditions)
# rows = 1
# cols = 3

# fig = make_subplots(
#     rows=rows,
#     cols=cols,
#     y_title=y_axis_name,
#     x_title=x_axis_name,
#     subplot_titles=["Diff < 0", "Diff = 0", "0 < Diff"],
#     # subplot_titles=subplot_titles,
#     shared_yaxes=True,
#     shared_xaxes=True,
# )

# algorithms = ["Ascending", "Descending"]
# diffs = ["< 0", "= 0", "> 0"]
# symbols = ["cross", "diamond", "circle"]

# for col, (diff, symbol) in enumerate(zip(diffs, symbols), start=1):
#     # for row, condition in enumerate(conditions, start=1):
#     for condition in conditions:
#         df = algs_diff_df.loc[
#             (algs_diff_df["Diff"] == diff) & (algs_diff_df[condition_col] == condition)
#         ]
#         x = df["Fraction"]
#         y = df["Desc - Asc"]

#         fig.add_trace(
#             go.Violin(
#                 x=x,
#                 y=y,
#                 name=condition,
#                 line_color=color_discrete_map[condition],
#                 # mode="markers",
#                 opacity=0.7,
#             ),
#             row=1,
#             col=col,
#         )

# min_y = algs_diff_df["Desc - Asc"].min()
# max_y = algs_diff_df["Desc - Asc"].max()
# abs_y_limit = max(abs(min_y), abs(max_y))
# # min_y = min_y * 1.1 if abs(abs_y_limit // min_y // 10) <= 1 else min_y - abs(abs_y_limit // min_y // 10)
# min_y = (
#     min_y * 1.1
#     if abs(abs_y_limit // min_y // 10) <= 1
#     else min_y - abs(abs_y_limit // min_y)
# )
# max_y = (
#     max_y * 1.1
#     if abs(abs_y_limit // max_y // 10) <= 1
#     else max_y + abs(abs_y_limit // max_y)
# )
# fig.update_yaxes(zerolinewidth=zerolinewidth, range=[min_y, max_y])

# fig.update_xaxes(tick0=0, dtick=0.2, tickangle=-60, matches="x", range=[0, 1.1])

# # fig.update_traces(orientation='h', side='positive', width=3, points=False)

# fig.update_layout(
#     title_text=head_title,
#     # title_y=0.95,
#     template=template,
#     # showlegend=False,
#     height=700,
#     # violinmode='group'
# )

# fig.show()


# %% [markdown]
# #### Solutions' sizes

# %%
min_max_fraction_1_distinct_prots_df = (
    distinct_unique_proteins_df.loc[distinct_unique_proteins_df["Fraction"] == 1.0]
    .groupby(condition_col)["NumOfProteins"]
    .agg(["min", "max"])
    .reset_index()
)
min_max_fraction_1_distinct_prots_df["%SolutionsDispersion"] = min_max_fraction_1_distinct_prots_df.apply(
    lambda x: 100 * (x["max"] - x["min"]) / x["max"],
    axis=1
)
min_max_fraction_1_distinct_prots_df

# %%
fig = px.bar(
    min_max_fraction_1_distinct_prots_df,
    x=condition_col,
    y="%SolutionsDispersion",
    color=condition_col,
    color_discrete_map=color_discrete_map,
    category_orders=category_orders,
    template=template,
)

fig.update_xaxes(title="")
fig.update_yaxes(title="% dispersion<br><sub>100 * (max - min) / max</sub>", range=[0, min_max_fraction_1_distinct_prots_df["%SolutionsDispersion"].max()])

fig.update_layout(
    showlegend=False, 
    # xaxis_title="", 
    width=max(70 * len(conditions), 300), height=400
)

fig.write_image(
    "%SolutionsDispersion - Illumina.svg",
    width=max(70 * len(conditions), 300),
    height=400
)


fig.show()


# %%
# fig = go.Figure(data=[go.Histogram(x=min_max_fraction_1_distinct_prots_df["%SolutionsDispersion"])])

# fig.update_xaxes(title="% dispersion<br><sub>100 * (max - min) / max</sub>")
# fig.update_yaxes(title="Transcripts")

# fig.update_layout(
#     # showlegend=False,
#     width=600,
#     height=400,
#     template=template
# )

# fig.write_image(
#     "%SolutionsDispersion - Illumina.svg",
#     width=600,
#     height=400
# )

# fig.show()

# %%
# y_axis_name = "Distinct unique proteins"
# head_title = "Distinct unique proteins vs. heuristic method"

# cols = min(facet_col_wrap, len(conditions), 5)
# rows = ceil(len(conditions) / cols)
# row_col_iter = list(product(range(1, rows + 1), range(1, cols + 1)))[: len(conditions)]

# fig = make_subplots(
#     rows=rows,
#     cols=cols,
#     print_grid=False,
#     y_title=y_axis_name,
#     subplot_titles=conditions,
#     # shared_yaxes=True,
# )

# algorithms = ["Ascending", "Descending"]

# max_y = 0

# for (row, col), condition in zip(row_col_iter, conditions):

#     df = distinct_unique_proteins_df.loc[
#         (distinct_unique_proteins_df[condition_col] == condition)
#         & (distinct_unique_proteins_df["Fraction"] == 1.0)
#     ]

#     xs = [
#         df.loc[df["Algorithm"] == algorithm, condition_col] for algorithm in algorithms
#     ]
#     ys = [
#         df.loc[df["Algorithm"] == algorithm, "NumOfProteins"]
#         for algorithm in algorithms
#     ]

#     max_y = max(max_y, max(ys[0]), max(ys[1]))

#     fig.add_trace(
#         go.Violin(
#             x=xs[0],
#             y=ys[0],
#             legendgrouptitle_text=condition,
#             legendgroup=condition,
#             name="Ascending",
#             # scalegroup="Ascending",
#             side="negative",
#             line_color=subcolors_discrete_map[condition][0],
#             # points="all"
#         ),
#         row=row,
#         col=col,
#     )

#     fig.add_trace(
#         go.Violin(
#             x=xs[1],
#             y=ys[1],
#             legendgroup=condition,
#             name="Descending",
#             # scalegroup="Descending",
#             side="positive",
#             line_color=subcolors_discrete_map[condition][1],
#             # points="all"
#         ),
#         row=row,
#         col=col,
#     )

# # https://stackoverflow.com/a/63221694/10249633
# for ax in fig["layout"]:
#     if ax[:5] == "xaxis":
#         fig["layout"][ax]["tickmode"] = "array"
#         fig["layout"][ax]["tickvals"] = [1]
#         fig["layout"][ax]["ticktext"] = [""]

# fig.update_layout(
#     title_text=head_title,
#     template=template,
#     legend_title_text=f"{condition_col}, Algorithm",
#     legend_font=dict(size=12),
#     legend_grouptitlefont=dict(size=9),
#     legend_tracegroupgap=4,
#     # violingap=0,
#     # violinmode='overlay'
#     height=200 * rows,
#     width=250 * cols,
# )
# # fig.update_yaxes(range=[0, max_y * 1.05])
# fig.show()


# %% [markdown]
# #### Jaccard (overlap of solutions)

# %%
# proteins_jaccard_dfs = []
# for condition in conditions:
#     df = distinct_unique_proteins_df.loc[distinct_unique_proteins_df[condition_col] == condition].reset_index(drop=True)
#     df["Proteins"] = df["Proteins"].apply(lambda x: set(x.split(",")))
#     jaccard_df = calc_jaccard_df(df, "Proteins")
#     proteins_jaccard_dfs.append(jaccard_df)
# proteins_jaccard_dfs[0]

# %%
# annotated_proteins_jaccard_dfs = []

# for condition, proteins_jaccard_df in zip(conditions, proteins_jaccard_dfs):

#     df = distinct_unique_proteins_df.loc[distinct_unique_proteins_df[condition_col] == condition].reset_index(drop=True)
#     df = (
#         distinct_unique_proteins_df.loc[
#             distinct_unique_proteins_df[condition_col] == condition,
#             [condition_col, "Fraction", "FractionRepetition", "Algorithm", "AlgorithmRepetition"]
#         ]
#         .reset_index(drop=True)
#     )
#     df = pd.concat([df, proteins_jaccard_df], axis=1)
#     index_details_dict = {
#         i: f'{row["Fraction"]}-{row["FractionRepetition"]}-{row["Algorithm"]}-{row["AlgorithmRepetition"]}'
#         for i, row in df.iterrows()
#     }
#     df = df.rename(columns=index_details_dict)

#     annotated_proteins_jaccard_dfs.append(df)

# annotated_proteins_jaccard_dfs[0]

# %% [markdown]
# ##### Heatmap

# %%
# todo uncomment

# df = annotated_proteins_jaccard_dfs[0]
# # df = df.drop("Gene", axis=1).set_index(["Fraction", "FractionRepetition", "Algorithm", "AlgorithmRepetition"])
# # df = df.set_axis(df.index, axis=1)
# df

# %%
# todo uncomment

# fractions = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
# total_len = len(df)
# fraction_len = total_len / len(fractions)
# middles = [int((fraction_len / 2) + (fraction_len * x)) for x in range(len(fractions))]
# middles

# %%
# todo uncomment

# ticks = []
# for x, fraction in enumerate(
#     [fraction for fraction in fractions for _ in range(int(fraction_len))], start=1
# ):
#     if x > total_len:
#         break
#     tick = str(fraction) if x in middles else ""
#     ticks.append(tick)

# %%
#

# for x, tick in enumerate(ticks, start=1):
#     if tick != "":
#         ic(x, tick)

# %%
#

# data = df.iloc[:, 5:].to_numpy()

# %%
#

# data.shape

# %%
# todo uncomment

# fig = px.imshow(
#     # data,
#     df.iloc[:, 5:],
#     labels=dict(x="Fraction", y="Fraction", color="Jaccard index"),
#     x=df["Fraction"],
#     y=df["Fraction"],
# )
# fig.update_xaxes(side="top")
# fig.show()

# %%
# import dash_bio
# https://plotly.com/python/clustergram/

# %%

# %%
#

# import seaborn as sns

# sns.set_theme()
# uniform_data = df.iloc[:, 5:].to_numpy()
# ax = sns.heatmap(uniform_data)

# %% [markdown]
# ##### Distribution

# %%
# todo uncomment

# melted_proteins_jaccard_dfs = []

# for annotated_proteins_jaccard_df in annotated_proteins_jaccard_dfs:

#     df = pd.melt(
#         annotated_proteins_jaccard_df,
#         id_vars=[
#             condition_col,
#             "Fraction",
#             "FractionRepetition",
#             "Algorithm",
#             "AlgorithmRepetition",
#         ],
#         var_name="ComparedAgainst",
#         value_name="JaccardIndex",
#     )

#     melted_proteins_jaccard_dfs.append(df)

# melted_proteins_jaccard_dfs[0]

# %%
# TODO - rerun this heavy cell after the notebook finished its running

# cols = len(conditions)

# x_title = "Jaccard index"
# y_title = "# solutions"
# title_text = "Distribution of Jaccard index between different solutions"

# fig = make_subplots(
#     rows=1,
#     # rows=2,
#     cols=cols,
#     subplot_titles=conditions,
#     shared_yaxes=True,
#     x_title=x_title,
#     y_title=y_title,
# )

# min_x = None
# max_x = 0

# algorithms = ["Ascending", "Descending"]

# for col, (condition, melted_proteins_jaccard_df) in enumerate(zip(conditions, melted_proteins_jaccard_dfs), start=1):

#     for i, alg in enumerate(algorithms):

#         x = melted_proteins_jaccard_df.loc[
#             melted_proteins_jaccard_df["Algorithm"] == alg, "JaccardIndex"
#         ]

#         fig.add_trace(
#             go.Histogram(
#                 x=x,
#                 marker_color=subcolors_discrete_map[condition][i],
#                 name=f"{condition}, {alg}",
#             ),
#             row=1,
#             col=col,
#         )

#         min_x = min(min_x, x.min()) if min_x else x.min()
#         max_x = max(max_x, x.max())

# fig.update_layout(
#     template=template,
#     barmode="overlay",  # Overlay both histograms
#     title_text=title_text,
#     legend_title_text=f"{condition_col}, Algorithm",
# )

# # fig.update_traces(opacity=0.75)  # Reduce opacity to see both histograms
# fig.update_xaxes(range=[min_x * 0.9, max_x * 1.1])

# # fig.show()
# fig.show(config={'staticPlot': True, 'responsive': False})

# %%
# TODO - rerun this heavy cell after the notebook finished its running

# fractions = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
# comparisons = ["Same<br>fraction", "Other<br>fractions"]
# algorithms = ["Ascending", "Descending"]

# x_title = "Fraction"
# y_title = "Jaccard index"
# title_text = "Intra- & inter-fraction Jaccard index comparisons"

# cols = len(conditions)
# rows = 2

# fig = make_subplots(
#     rows=rows,
#     cols=cols,
#     column_titles=conditions,
#     row_titles=comparisons,
#     # horizontal_spacing=0./rows,
#     shared_xaxes=True,
#     shared_yaxes=True,
#     x_title=x_title,
#     y_title=y_title,
# )

# for col, (condition, melted_proteins_jaccard_df) in enumerate(zip(conditions, melted_proteins_jaccard_dfs), start=1):

#     # melted_proteins_jaccard_df = melted_proteins_jaccard_df.sample(200) # todo temp

#     intra_fractions_lines = melted_proteins_jaccard_df.apply(lambda x: x["ComparedAgainst"].startswith(str(x["Fraction"])), axis=1)

#     for row in range(1, rows + 1):

#         if row == 1:
#             inter_or_intra_df = melted_proteins_jaccard_df.loc[intra_fractions_lines]
#             self_lines = inter_or_intra_df.apply(lambda x: f"{x['Fraction']}-{x['FractionRepetition']}-{x['Algorithm']}-{x['AlgorithmRepetition']}" == x["ComparedAgainst"], axis=1)
#             inter_or_intra_df = inter_or_intra_df.loc[~self_lines]
#         else:
#             inter_or_intra_df = melted_proteins_jaccard_df.loc[~intra_fractions_lines]

#         for i, alg in enumerate(algorithms):

#             alg_df = inter_or_intra_df.loc[
#                 inter_or_intra_df["Algorithm"] == alg
#             ]
#             x = alg_df["Fraction"]
#             y = alg_df["JaccardIndex"]

#             fig.add_trace(
#                     go.Box(
#                         x=x,
#                         y=y,
#                         marker_color=subcolors_discrete_map[condition][i],
#                         name=f"{condition}, {alg}",
#                         showlegend=row==1
#                     ),
#                     row=row,
#                     col=col,
#                 )

# fig.update_layout(
#     template=template,
#     boxmode="group",
#     title_text=title_text,
#     legend_title_text=f"{condition_col}, Algorithm",
#     height=600,
#     # width=800
# )

# fig.update_xaxes(
#     tick0 = 0,
#     dtick = 0.1,
#     tickangle = -60,
#     matches='x'
# )

# fig.update_yaxes(range=[0, 1])

# # fig.show()
# fig.show(config={'staticPlot': True, 'responsive': False})

# %%

# %% [markdown] papermill={"duration": 0.030615, "end_time": "2022-02-01T09:42:49.024262", "exception": false, "start_time": "2022-02-01T09:42:48.993647", "status": "completed"}
# ## Supporting reads' coverage

# %% [markdown] papermill={"duration": 0.030615, "end_time": "2022-02-01T09:42:49.024262", "exception": false, "start_time": "2022-02-01T09:42:48.993647", "status": "completed"}
# ### Unique reads

# %%
unique_reads_dfs[0]


# %%
# s = unique_reads_dfs[0]["NumOfReads"].sum()
# s

# %%
# unique_reads_dfs[0]["NumOfReads"].cumsum()

# %%
# 100 * unique_reads_dfs[0]["NumOfReads"].cumsum() / s

# %%
# max(df["NumOfReads"].max() for df in unique_reads_dfs)


# %%
# df = unique_reads_dfs[0].assign(
#     CumNumOfReads=unique_reads_dfs[0]["NumOfReads"].cumsum()
# )
# df["%CumNumOfReads"] = 100 * df["CumNumOfReads"] / df["NumOfReads"].sum()
# df["UniqueReadWithDecreasingSupport"] = df.index
# df

# %%
# cummulative_supporting_reads_dfs = []
# for unique_reads_df in unique_reads_dfs:
#     df = (
#         unique_reads_df.loc[:, ["UniqueRead", "NumOfReads"]]
#         .sort_values("NumOfReads", ascending=False)
#         .reset_index(drop=True)
#     )
#     df = df.assign(CumNumOfReads=df["NumOfReads"].cumsum())
#     df["%CumNumOfReads"] = 100 * df["CumNumOfReads"] / df["NumOfReads"].sum()
#     df["NumUniqueReads"] = df.index + 1
#     cummulative_supporting_reads_dfs.append(df)
# cummulative_supporting_reads_dfs[0]

# %%
# df = cummulative_supporting_reads_dfs[0]
# df

# %%
# df.iloc[:100]

# %%
# ax = sns.barplot(
#     # data=df.iloc[:10000],
#     data=df,
#     x="NumUniqueReads",
#     y="%CumNumOfReads",
#     color="blue",
#     # dodge=False
# )
# ax.set(xscale="log")

# %% papermill={"duration": 4.052404, "end_time": "2022-02-01T09:42:53.176715", "exception": false, "start_time": "2022-02-01T09:42:49.124311", "status": "completed"}
# # cummulative_supporting_reads_dfs = []
# # for unique_reads_df in unique_reads_dfs:
# #     df = unique_reads_df.loc[:, ["UniqueRead", "NumOfReads"]].sort_values("NumOfReads", ascending=False).reset_index(drop=True)
# #     df = df.assign(CumNumOfReads = df["NumOfReads"].cumsum())
# #     df["%CumNumOfReads"] = 100 * df["CumNumOfReads"] / df["NumOfReads"].sum()
# #     df["NumUniqueReads"] = df.index + 1
# #     cummulative_supporting_reads_dfs.append(df)


# cols = len(conditions)

# x_title = "Cummulative unique reads"
# y_title = "Cummulative<br>supporting reads (%)"
# title_text = "Unique reads' abundance"

# fig = make_subplots(
#     rows=1,
#     cols=cols,
#     subplot_titles=conditions,
#     shared_yaxes=True,
#     shared_xaxes=True,
#     x_title=x_title,
#     y_title=y_title,
# )

# for col, condition, unique_reads_df in zip(
#     range(1, cols + 1), conditions, unique_reads_dfs
# ):

#     df = unique_reads_df.loc[:, ["UniqueRead", "NumOfReads"]].sort_values("NumOfReads", ascending=False).reset_index(drop=True)
#     df = df.assign(CumNumOfReads = df["NumOfReads"].cumsum())
#     df["%CumNumOfReads"] = 100 * df["CumNumOfReads"] / df["NumOfReads"].sum()
#     df["NumUniqueReads"] = df.index + 1

#     # ax = sns.barplot(x="day", y="tip", data=tips, ci=68)
#     ax = sns.barplot(
#         data=df,
#         x="NumUniqueReads",
#         y="%CumNumOfReads",
#         color=color_discrete_map[condition],
#         # dodge=False
#     )
#     ax.set(xscale="log")

#     # fig.add_trace(
#     #     go.Bar(
#     #         x=df["NumUniqueReads"],
#     #         y=df["%CumNumOfReads"],
#     #         marker_color=color_discrete_map[condition]
#     #     ),
#     #     row=1,
#     #     col=col,
#     # )

# # fig.update_layout(
# #     title_text=title_text,
# #     showlegend=False,
# #     template=template
# # )

# # fig.update_xaxes(type="log")
# # fig.update_yaxes(type="log")

# # fig.show(config={'staticPlot': True, 'responsive': False})


# %% papermill={"duration": 4.052404, "end_time": "2022-02-01T09:42:53.176715", "exception": false, "start_time": "2022-02-01T09:42:49.124311", "status": "completed"}
# # cummulative_supporting_reads_dfs = []
# # for unique_reads_df in unique_reads_dfs:
# #     df = unique_reads_df.loc[:, ["UniqueRead", "NumOfReads"]].sort_values("NumOfReads", ascending=False).reset_index(drop=True)
# #     df = df.assign(CumNumOfReads = df["NumOfReads"].cumsum())
# #     df["%CumNumOfReads"] = 100 * df["CumNumOfReads"] / df["NumOfReads"].sum()
# #     df["NumUniqueReads"] = df.index + 1
# #     cummulative_supporting_reads_dfs.append(df)


# cols = len(conditions)

# x_title = "Cummulative unique reads"
# y_title = "Cummulative<br>supporting reads (%)"
# title_text = "Unique reads' abundance"

# fig = make_subplots(
#     rows=1,
#     cols=cols,
#     subplot_titles=conditions,
#     shared_yaxes=True,
#     shared_xaxes=True,
#     x_title=x_title,
#     y_title=y_title,
# )

# for col, condition, unique_reads_df in zip(
#     range(1, cols + 1), conditions, unique_reads_dfs
# ):

#     df = unique_reads_df.loc[:, ["UniqueRead", "NumOfReads"]].sort_values("NumOfReads", ascending=False).reset_index(drop=True)
#     df = df.assign(CumNumOfReads = df["NumOfReads"].cumsum())
#     df["%CumNumOfReads"] = 100 * df["CumNumOfReads"] / df["NumOfReads"].sum()
#     df["NumUniqueReads"] = df.index + 1

#     fig.add_trace(
#         go.Bar(
#             x=df["NumUniqueReads"],
#             y=df["%CumNumOfReads"],
#             marker_color=color_discrete_map[condition]
#         ),
#         row=1,
#         col=col,
#     )

# fig.update_layout(
#     title_text=title_text,
#     showlegend=False,
#     template=template
# )

# fig.update_xaxes(type="log")
# fig.update_yaxes(type="log")

# fig.show(config={'staticPlot': True, 'responsive': False})


# %% papermill={"duration": 4.052404, "end_time": "2022-02-01T09:42:53.176715", "exception": false, "start_time": "2022-02-01T09:42:49.124311", "status": "completed"}
# # cummulative_supporting_reads_dfs = []
# # for unique_reads_df in unique_reads_dfs:
# #     df = unique_reads_df.loc[:, ["UniqueRead", "NumOfReads"]].sort_values("NumOfReads", ascending=False).reset_index(drop=True)
# #     df = df.assign(CumNumOfReads = df["NumOfReads"].cumsum())
# #     df["%CumNumOfReads"] = 100 * df["CumNumOfReads"] / df["NumOfReads"].sum()
# #     df["NumUniqueReads"] = df.index + 1
# #     cummulative_supporting_reads_dfs.append(df)


# cols = len(conditions)

# x_title = "Cummulative unique reads"
# y_title = "Cummulative<br>supporting reads (%)"
# title_text = "Unique reads' abundance"

# fig = make_subplots(
#     rows=1,
#     cols=cols,
#     subplot_titles=conditions,
#     shared_yaxes=True,
#     shared_xaxes=True,
#     x_title=x_title,
#     y_title=y_title,
# )

# for col, condition, unique_reads_df in zip(
#     range(1, cols + 1), conditions, unique_reads_dfs
# ):

#     df = unique_reads_df.loc[:, ["UniqueRead", "NumOfReads"]].sort_values("NumOfReads", ascending=False).reset_index(drop=True)
#     df = df.assign(CumNumOfReads = df["NumOfReads"].cumsum())
#     df["%CumNumOfReads"] = 100 * df["CumNumOfReads"] / df["NumOfReads"].sum()
#     df["NumUniqueReads"] = df.index + 1

#     fig.add_trace(
#         go.Bar(
#             x=df["NumUniqueReads"],
#             y=df["%CumNumOfReads"],
#             marker_color=color_discrete_map[condition],
#         ),
#         row=1,
#         col=col,
#     )

# fig.update_layout(
#     title_text=title_text,
#     showlegend=False,
#     template=template
# )

# fig.update_xaxes(type="log")
# fig.update_yaxes(type="log")

# fig.update_traces(marker_line_width=0)

# fig.show(config={'staticPlot': True, 'responsive': False})


# %%

# %%

# %% [markdown] papermill={"duration": 0.030615, "end_time": "2022-02-01T09:42:49.024262", "exception": false, "start_time": "2022-02-01T09:42:48.993647", "status": "completed"}
# ### Unique proteins

# %%
# unique_proteins_dfs[0].head()


# %%
# max(df["NumOfReads"].max() for df in unique_proteins_dfs)


# %% papermill={"duration": 4.052404, "end_time": "2022-02-01T09:42:53.176715", "exception": false, "start_time": "2022-02-01T09:42:49.124311", "status": "completed"}
# cols = len(conditions)

# x_title = "Reads"
# y_title = "log<sub>10</sub> unique proteins"
# title_text = "Distribution of supporting reads per unique protein"

# fig = make_subplots(
#     rows=1,
#     cols=cols,
#     subplot_titles=conditions,
#     shared_yaxes=True,
#     shared_xaxes=True,
#     x_title=x_title,
#     y_title=y_title,
# )

# for col, condition, unique_proteins_df in zip(
#     range(1, cols + 1), conditions, unique_proteins_dfs
# ):

#     fig.add_trace(
#         go.Histogram(
#             x=unique_proteins_df["NumOfReads"],
#             marker_color=color_discrete_map[condition],
#         ),
#         row=1,
#         col=col,
#     )

# fig.update_layout(
#     title_text=title_text,
#     showlegend=False,
#     template=template,
# )

# # min_x = min(0, min(df["NumOfReads"].min() for df in unique_proteins_dfs))
# min_x = 0
# max_x = max(df["NumOfReads"].max() for df in unique_proteins_dfs) * 1.05
# fig.update_xaxes(range=[min_x, max_x])
# fig.update_yaxes(type="log")

# fig.show()


# %% [markdown] papermill={"duration": 0.030615, "end_time": "2022-02-01T09:42:49.024262", "exception": false, "start_time": "2022-02-01T09:42:48.993647", "status": "completed"}
# ### Distinct unique proteins

# %% [markdown] toc-hr-collapsed=true jp-MarkdownHeadingCollapsed=true
# #### Jaccard - TODO - erase?

# %%
# gdf = (
#     expanded_distinct_unique_proteins_df
#     .groupby([condition_col, "Fraction", "FractionRepetition", "Algorithm", "AlgorithmRepetition"])
# )

# %%
# fraction_1_dfs = []

# for condition in conditions:
#     df = expanded_distinct_unique_proteins_df.loc[
#         (expanded_distinct_unique_proteins_df[condition_col] == condition) &
#         (expanded_distinct_unique_proteins_df["Fraction"] == 1.0)
#     ]
#     df["Reads"] = df["Reads"].apply(lambda x: x.split(","))
#     df = df.groupby(["FractionRepetition", "Algorithm", "AlgorithmRepetition"])["Reads"].apply(lambda x: list(chain(x))).reset_index()
#     df["Reads"] = df["Reads"].apply(lambda x: set(chain.from_iterable(x)))
#     df["NumOfReads"] = df["Reads"].apply(lambda x: len(x))
#     fraction_1_dfs.append(df)

# fraction_1_dfs[0]

# %%
# def jaccard_index(s1: set, s2: set):
#     return len(s1 & s2) / len(s1 | s2)

# %%
# def calc_jaccard_df(df, col="Proteins"):
#     return df[col].apply(
#         lambda x: df[col].apply(lambda y: jaccard_index(x, y))
#     ).reset_index(drop=True)

# %% [markdown]
# ##### By sample

# %%
# condition = conditions[0]

# %%
# df = distinct_unique_proteins_df.loc[distinct_unique_proteins_df[condition_col] == condition]
# df = df.loc[df["Fraction"] == 1.0].reset_index(drop=True)
# df["Proteins"] = df["Proteins"].apply(lambda x: x.split(","))
# # df["Proteins"] = df["Proteins"].apply(lambda x: set(chain.from_iterable(x)))
# df["Proteins"] = df["Proteins"].apply(lambda x: set(x))

# df

# %%
# calc_jaccard_df(df, "Proteins")

# %%

# %%
# fraction_1_jaccard_dfs = []
# for condition in conditions:
#     df = distinct_unique_proteins_df.loc[distinct_unique_proteins_df[condition_col] == condition]
#     df = df.loc[df["Fraction"] == 1.0].reset_index(drop=True)
#     df["Proteins"] = df["Proteins"].apply(lambda x: x.split(","))
#     df["Proteins"] = df["Proteins"].apply(lambda x: set(x))
#     fraction_1_jaccard_dfs.append(calc_jaccard_df(df, "Proteins"))
# fraction_1_jaccard_dfs[0]

# %%
# fraction_1_jaccard_dfs[0].melt()

# %%
# cols = len(conditions)

# fig = make_subplots(
#     rows=1,
#     cols=cols,
#     subplot_titles=conditions,
#     shared_yaxes=True,
#     shared_xaxes=True,
#     # x_title="PC1",
#     # y_title="PC2",
# )

# min_x = None
# max_y = None
# max_x = None
# min_y = None


# for col, condition, jaccard_df in zip(
#     range(1, cols + 1), conditions, fraction_1_jaccard_dfs
# ):

#     f1_df = distinct_unique_proteins_df.loc[distinct_unique_proteins_df[condition_col] == condition]
#     f1_df = f1_df.loc[f1_df["Fraction"] == 1.0].reset_index(drop=True)

#     for i, alg in enumerate(f1_df["Algorithm"].unique()):

#         alg_components = components.loc[components["Algorithm"] == alg]
#         x = alg_components["PC1"]
#         y = alg_components["PC2"]

#         fig.add_trace(
#             go.Histogram(
#                 x=x,
#                 y=y,
#                 mode="markers",
#                 marker_color=subcolors_discrete_map[condition][i],
#                 name=f"{condition}, {alg}",
#             ),
#             row=1,
#             col=col,
#         )

#     fig.add_annotation(
#         row=1,
#         col=col,
#         # align="right",
#         x=min_x+abs(min_x)*0.4,
#         y=max_y-abs(max_y)*0.34,
#         xref="x",
#         yref="y",
#         text=(
#             "<b>Variance explaind</b>"
#             "<br>"
#             f"PC1 = {pca.explained_variance_ratio_[0] * 100:.2f}%"
#             "<br>"
#             f"PC2 = {pca.explained_variance_ratio_[1] * 100:.2f}%"
#         ),
#         bgcolor="white",
#         opacity=0.7,
#         borderpad=4,
#     )

# fig.update_layout(
#     title_text="Jaccard Index PCA on reads supporting different sets of distinct unique proteins",
#     legend_title_text=f"{condition_col}, Algorithm",
#     template=template,
# )

# fig.update_xaxes(range=[min_x-abs(min_x)*0.2, max_x+abs(max_x)*0.2])
# # fig.update_yaxes(range=[min_y, max_y])

# fig.show()


# %%

# %% papermill={"duration": 4.052404, "end_time": "2022-02-01T09:42:53.176715", "exception": false, "start_time": "2022-02-01T09:42:49.124311", "status": "completed"}
# cols = len(conditions)

# fig = make_subplots(
#     rows=1,
#     cols=cols,
#     subplot_titles=conditions,
#     shared_yaxes=True,
#     shared_xaxes=True,
#     x_title="PC1",
#     y_title="PC2",
# )

# min_x = None
# max_y = None
# max_x = None
# min_y = None

# components_dfs = []
# pcas = []
# for condition, f1_jaccard_df in zip(
#     conditions, fraction_1_jaccard_dfs
# ):
#     f1_df = distinct_unique_proteins_df.loc[distinct_unique_proteins_df[condition_col] == condition]
#     f1_df = f1_df.loc[f1_df["Fraction"] == 1.0].reset_index(drop=True)
#     X = f1_jaccard_df
#     features_names = X.columns
#     pca = PCA(n_components=2, random_state=1892)  # our lord & savior J.R.R.T was born in 1892
#     components = pca.fit_transform(X)
#     components = pd.concat([f1_df.loc[:, ["FractionRepetition", "Algorithm", "AlgorithmRepetition"]], pd.DataFrame(components)], axis=1).rename(columns={0: "PC1", 1: "PC2"})
#     min_x = min(min_x, components["PC1"].min()) if min_x else components["PC1"].min()
#     max_y = max(max_y, components["PC2"].max()) if max_y else components["PC2"].max()
#     max_x = max(max_x, components["PC1"].max()) if max_x else components["PC1"].max()
#     min_y = min(min_y, components["PC2"].min()) if min_y else components["PC2"].min()
#     components_dfs.append(components)
#     pcas.append(pca)

# for col, condition, components, f1_df, pca in zip(
#     range(1, cols + 1), conditions, components_dfs, fraction_1_dfs, pcas
# ):

#     for i, alg in enumerate(f1_df["Algorithm"].unique()):

#         alg_components = components.loc[components["Algorithm"] == alg]
#         x = alg_components["PC1"]
#         y = alg_components["PC2"]

#         fig.add_trace(
#             go.Scatter(
#                 x=x,
#                 y=y,
#                 mode="markers",
#                 marker_color=subcolors_discrete_map[condition][i],
#                 name=f"{condition}, {alg}",
#             ),
#             row=1,
#             col=col,
#         )

#     fig.add_annotation(
#         row=1,
#         col=col,
#         # align="right",
#         x=min_x+abs(min_x)*0.4,
#         y=max_y-abs(max_y)*0.34,
#         xref="x",
#         yref="y",
#         text=(
#             "<b>Variance explaind</b>"
#             "<br>"
#             f"PC1 = {pca.explained_variance_ratio_[0] * 100:.2f}%"
#             "<br>"
#             f"PC2 = {pca.explained_variance_ratio_[1] * 100:.2f}%"
#         ),
#         bgcolor="white",
#         opacity=0.7,
#         borderpad=4,
#     )

# fig.update_layout(
#     title_text="Jaccard Index PCA on reads supporting different sets of distinct unique proteins",
#     legend_title_text=f"{condition_col}, Algorithm",
#     template=template,
# )

# fig.update_xaxes(range=[min_x-abs(min_x)*0.2, max_x+abs(max_x)*0.2])
# # fig.update_yaxes(range=[min_y, max_y])

# fig.show()


# %% [markdown]
# ##### By sample & algorithm

# %%
# algs = ["Ascending", "Descending"]

# %%
# fraction_1_alg_condition_jaccard_dfs = {}
# for df, condition in zip(fraction_1_dfs, conditions):
#     for alg in algs:
#         jaccard_df = calc_jaccard_df(df.loc[df["Algorithm"] == alg], "Reads")
#         fraction_1_alg_condition_jaccard_dfs[(condition, alg)] = jaccard_df

# %% papermill={"duration": 4.052404, "end_time": "2022-02-01T09:42:53.176715", "exception": false, "start_time": "2022-02-01T09:42:49.124311", "status": "completed"}
# cols = len(conditions)
# rows = 2

# # subplot_titles = [f"{condition}, {alg}" for condition in conditions for alg in algs]
# subplot_titles = [f"{condition}, {alg}" for alg in algs for condition in conditions ]

# fig = make_subplots(
#     rows=rows,
#     cols=cols,
#     subplot_titles=subplot_titles,
#     shared_yaxes=True,
#     shared_xaxes=True,
#     x_title="PC1",
#     y_title="PC2",
# )

# min_x = None
# max_y = None
# max_x = None
# min_y = None

# components_dfs = {}
# pcas = {}
# for col, (f1_df, condition) in enumerate(zip(fraction_1_dfs, conditions), start=1):
#     for row, alg in enumerate(algs, start=1):
#         f1_alg_df = f1_df.loc[f1_df["Algorithm"] == alg]
#         X = fraction_1_alg_condition_jaccard_dfs[(condition, alg)]
#         features_names = X.columns
#         pca = PCA(n_components=2, random_state=1892)  # our lord & savior J.R.R.T was born in 1892
#         components = pca.fit_transform(X)
#         components = pd.concat([f1_alg_df.loc[:, ["FractionRepetition", "Algorithm", "AlgorithmRepetition"]], pd.DataFrame(components)], axis=1).rename(columns={0: "PC1", 1: "PC2"})
#         min_x = min(min_x, components["PC1"].min()) if min_x else components["PC1"].min()
#         max_y = max(max_y, components["PC2"].max()) if max_y else components["PC2"].max()
#         max_x = max(max_x, components["PC1"].max()) if max_x else components["PC1"].max()
#         min_y = min(min_y, components["PC2"].min()) if min_y else components["PC2"].min()
#         components_dfs[(condition, alg)] = components
#         pcas[(condition, alg)] = pca

# ic(min_x, max_x)
# ic(min_y, max_y)

# for col, (condition, components, f1_df, pca) in enumerate(
#     zip(conditions, components_dfs, fraction_1_dfs, pcas),
#     start=1
# ):
#     for row, alg in enumerate(algs, start=1):
#         components = components_dfs[(condition, alg)]
#         x = components["PC1"]
#         y = components["PC2"]
#         fig.add_trace(
#             go.Scatter(
#                 x=x,
#                 y=y,
#                 mode="markers",
#                 marker_color=subcolors_discrete_map[condition][row - 1],
#                 # name=f"{condition}, {alg}",
#             ),
#             row=row,
#             col=col,
#         )
#         pca = pcas[(condition, alg)]
#         fig.add_annotation(
#             row=row,
#             col=col,
#             # align="right",
#             font=dict(
#                 size=10,
#             ),
#             showarrow=False,
#             # x=min_x+abs(min_x)*0.4,
#             x=min_x+abs(min_x)*1.5,
#             y=max_y+2.25,
#             # y=3,
#             xref="x",
#             yref="y",
#             text=(
#                 "<b>Variance explaind</b>"
#                 "<br>"
#                 f"PC1 = {pca.explained_variance_ratio_[0] * 100:.2f}%"
#                 "<br>"
#                 f"PC2 = {pca.explained_variance_ratio_[1] * 100:.2f}%"
#             ),
#             bgcolor="white",
#             opacity=0.7,
#             borderpad=4,

#         )

# fig.update_layout(
#     title_text="Jaccard Index PCA on reads supporting different sets of distinct unique proteins",
#     # legend_title_text=f"{condition_col}, Algorithm",
#     showlegend=False,
#     template=template,
# )

# fig.update_xaxes(range=[min_x-abs(min_x)*0.2, max_x+abs(max_x)*0.2])
# # fig.update_yaxes(range=[min_y, max_y])

# fig.show()


# %% [markdown]
# ##### By sample, algorithm & fraction repetition

# %%
# # jaccards per fraction repetition

# f_r_jaccard_dfs = {}
# for fraction_1_df, condition in zip(fraction_1_dfs, conditions):
#     f_r_condition_jaccard_dfs = {}
#     for f_r in np.sort(fraction_1_df["FractionRepetition"].unique()):
#         f_r_df = fraction_1_df.loc[fraction_1_df["FractionRepetition"] == f_r].reset_index(drop=True)
#         f_r_condition_jaccard_dfs[f_r] = f_r_df["Reads"].apply(
#             lambda x: f_r_df["Reads"].apply(lambda y: jaccard_index(x, y))
#         )
#     f_r_jaccard_dfs[condition] = f_r_condition_jaccard_dfs
# f_r_jaccard_dfs[conditions[0]][1]

# %%
# max_cols = 5

# for condition, f1_df in zip(conditions, fraction_1_dfs):

#     f_r_condition_jaccard_dfs = f_r_jaccard_dfs[condition]

#     frs = sorted(list(f_r_condition_jaccard_dfs.keys()))
#     plots = len(frs)
#     cols = min(plots, max_cols)
#     rows = ceil(plots / cols)
#     rows_and_cols = [
#         (row, col)
#         for row in range(1, rows + 1)
#         for col in range(1, cols + 1)
#         if row * col <= plots
#     ]
#     # frs_subcolors_discrete_map = n_repetitions_colormap(subcolors_discrete_map, condition, plots)

#     fig = make_subplots(
#         rows=rows,
#         cols=cols,
#         subplot_titles=[str(fr) for fr in frs],
#         shared_yaxes=True,
#         shared_xaxes=True,
#         x_title="PC1",
#         y_title="PC2",
#     )

#     min_x = None
#     max_y = None
#     max_x = None
#     min_y = None

#     components_dfs = []
#     pcas = []

#     for fr in f_r_condition_jaccard_dfs:

#         fr_f1_df = f1_df.loc[f1_df["FractionRepetition"] == fr].reset_index(drop=True)
#         X = f_r_condition_jaccard_dfs[fr]
#         features_names = X.columns
#         pca = PCA(n_components=2, random_state=1892)  # our lord & savior J.R.R.T was born in 1892
#         components = pca.fit_transform(X)
#         components = (
#             pd.concat(
#                 [
#                     fr_f1_df.loc[:, ["FractionRepetition", "Algorithm", "AlgorithmRepetition"]],
#                     pd.DataFrame(components)
#                 ],
#                 axis=1
#             )
#             .rename(columns={0: "PC1", 1: "PC2"})
#         )
#         min_x = min(min_x, components["PC1"].min()) if min_x else components["PC1"].min()
#         max_y = max(max_y, components["PC2"].max()) if max_y else components["PC2"].max()
#         max_x = max(max_x, components["PC1"].max()) if max_x else components["PC1"].max()
#         min_y = min(min_y, components["PC2"].min()) if min_y else components["PC2"].min()
#         components_dfs.append(components)
#         pcas.append(pca)

#     for (row, col), components, pca, fr in zip(
#         rows_and_cols, components_dfs, pcas, frs
#     ):

#         for i, alg in enumerate(components["Algorithm"].unique()):

#             alg_components = components.loc[components["Algorithm"] == alg]
#             x = alg_components["PC1"]
#             y = alg_components["PC2"]

#             fig.add_trace(
#                 go.Scatter(
#                     x=x,
#                     y=y,
#                     mode="markers",
#                     marker_color=subcolors_discrete_map[condition][i],
#                     name=f"{alg}, {fr}",
#                     # name=alg,
#                 ),
#                 row=row,
#                 col=col,
#             )

#         # fig.add_annotation(
#         #     row=row,
#         #     col=col,
#         #     # align="right",
#         #     x=min_x+abs(min_x)*0.4,
#         #     y=max_y-abs(max_y)*0.34,
#         #     xref="x",
#         #     yref="y",
#         #     text=(
#         #         "<b>Variance explaind</b>"
#         #         "<br>"
#         #         f"PC1 = {pca.explained_variance_ratio_[0] * 100:.2f}%"
#         #         "<br>"
#         #         f"PC2 = {pca.explained_variance_ratio_[1] * 100:.2f}%"
#         #     ),
#         #     bgcolor="white",
#         #     opacity=0.7,
#         #     borderpad=4,
#         # )

#     fig.update_layout(
#         title_text=f"Jaccard Index PCA on reads supporting different sets of distinct unique proteins in {condition}",
#         # legend_title_text=f"{condition_col}, Algorithm",
#         legend_title_text="Algorithm",
#         template=template,
#     )

#     fig.update_xaxes(range=[min_x-abs(min_x)*0.2, max_x+abs(max_x)*0.2])
#     # fig.update_yaxes(range=[min_y, max_y])

#     fig.show()

# %%

# %%

# %%

# %%

# %% [markdown]
# #### Expression levels

# %% [markdown]
# ##### Relative expression of isoforms

# %%
expression_dfs = []
for expression_file in expression_files:
    expression_df = pd.read_csv(expression_file, sep=sep)
    expression_df["#Solution"] = expression_df["#Solution"].astype(str)
    expression_df["Diff5+"] = (
        abs(
            expression_df["TotalEqualSupportingReads"]
            - expression_df["TotalWeightedSupportingReads"]
        )
        >= 0.05
        * (
            expression_df["TotalEqualSupportingReads"]
            + expression_df["TotalWeightedSupportingReads"]
        )
        / 2
    )
    expression_dfs.append(expression_df)
expression_dfs[0]


# %%
def find_rand_maximal_solution(
    expression_df, seed, allowed_algorithms=["Ascending", "Descending"]
):
    df = (
        expression_df.loc[expression_df["Algorithm"].isin(allowed_algorithms)]
        .groupby("#Solution")
        .agg("size")
        .reset_index()
        .rename(columns={0: "Size"})
    )
    # rand_maximal_solution = df.loc[df["Size"] == df["Size"].max(), "#Solution"].sample(random_state=seed).reset_index(drop=True)
    rand_maximal_solution = (
        df.loc[df["Size"] == df["Size"].max(), "#Solution"]
        .sample(random_state=seed)
        .values[0]
    )
    return rand_maximal_solution


# %%
def make_percentile_df(
    expression_df,
    first_percentile=10,
    inclusive_last_percentile=110,
    percentile_step=10,
    allowed_algorithms=["Ascending", "Descending"],
):

    gb = expression_df.loc[expression_df["Algorithm"].isin(allowed_algorithms)].groupby(
        "#Solution"
    )
    solutions_expression_dfs = [gb.get_group(x) for x in gb.groups]

    equal_supp_reads_dfs = []
    weighted_supp_reads_dfs = []
    for df in solutions_expression_dfs:
        equal_df = df.sort_values(
            "TotalEqualSupportingReads", ascending=False
        ).reset_index(drop=True)
        equal_df["CumTotalEqualSupportingReads"] = equal_df[
            "TotalEqualSupportingReads"
        ].cumsum()
        equal_df["%CumTotalEqualSupportingReads"] = (
            100
            * equal_df["CumTotalEqualSupportingReads"]
            / equal_df["TotalEqualSupportingReads"].sum()
        )
        equal_supp_reads_dfs.append(equal_df)
        weighted_df = df.sort_values(
            "TotalWeightedSupportingReads", ascending=False
        ).reset_index(drop=True)
        weighted_df["CumTotalWeightedSupportingReads"] = weighted_df[
            "TotalWeightedSupportingReads"
        ].cumsum()
        weighted_df["%CumTotalWeightedSupportingReads"] = (
            100
            * weighted_df["CumTotalWeightedSupportingReads"]
            / weighted_df["TotalWeightedSupportingReads"].sum()
        )
        weighted_supp_reads_dfs.append(weighted_df)

    # equal_supp_reads_dfs, weighted_supp_reads_dfs = make_supp_reads_dfs(expression_df)

    solutions = []
    assignment_methods = []
    percentiles = []
    required_proteins = []
    algorithms = []
    _conditions = []

    for dfs, method, col in zip(
        [equal_supp_reads_dfs, weighted_supp_reads_dfs],
        ["Equal", "Weighted"],
        ["%CumTotalEqualSupportingReads", "%CumTotalWeightedSupportingReads"],
    ):
        for df in dfs:
            # ic(df.iloc[:1, :3], method, col)
            # break
            solution = df.loc[0, "#Solution"]
            algorithm = df.loc[0, "Algorithm"]
            _condition = df.loc[0, condition_col]
            a = df[col].to_numpy()
            # for percentile in range(50, 100, 10):
            # for percentile in range(10, 110, 10):
            for percentile in range(
                first_percentile, inclusive_last_percentile, percentile_step
            ):
                idx = (np.abs(a - percentile)).argmin()
                if a[idx] < percentile:
                    idx += 1
                solutions.append(solution)
                assignment_methods.append(method)
                percentiles.append(percentile)
                required_proteins.append(idx)
                algorithms.append(algorithm)
                _conditions.append(_condition)

    percentile_df = pd.DataFrame(
        {
            "#Solution": solutions,
            "AssignmentMethod": assignment_methods,
            "Percentile": percentiles,
            "RequiredProteins": required_proteins,
            "Algorithm": algorithms,
            condition_col: _conditions,
        }
    )

    return percentile_df


# %%
def choose_sample_solutions(
    expression_df, seed, allowed_algorithms=["Ascending", "Descending"]
):
    return (
        expression_df.loc[
            expression_df["Algorithm"].isin(allowed_algorithms),
            [
                condition_col,
                "#Solution",
                "Fraction",
                "FractionRepetition",
                "Algorithm",
                "AlgorithmRepetition",
            ],
        ]
        .groupby(["Algorithm", "#Solution"])
        .sample()
        .groupby("Algorithm")
        .sample(3, random_state=seed)
        .reset_index(drop=True)["#Solution"]
    )


# %%
maximal_solutions = [
    find_rand_maximal_solution(expression_df, seed, allowed_algorithms=["Ascending", "Descending"])
    for expression_df in expression_dfs
]
maximal_solutions

# %%
percentile_dfs = [
    make_percentile_df(
        expression_df.loc[expression_df["#Solution"] == maximal_solution].reset_index(
            drop=True
        ),
        allowed_algorithms=["Ascending", "Descending"],
    )
    for expression_df, maximal_solution in zip(expression_dfs, maximal_solutions)
]
percentile_dfs[0]

# %%
# _seeds = [np.random.default_rng(seed)]
# for _ in conditions[1:]:
#     _seeds.append(np.random.default_rng(_seeds[-1]))

# all_conditions_sample_solutions = [
#     choose_sample_solutions(expression_df, _seed, allowed_algorithms=["Descending"]) for expression_df, _seed in zip(expression_dfs, _seeds)
# ]
# all_conditions_sample_solutions[0]

# %%
# sample_supp_reads_dfs = [
#     make_supp_reads_dfs(expression_df, sample_solutions)
#     for expression_df, sample_solutions in zip(expression_dfs, all_conditions_sample_solutions)
# ]

# ic(len(sample_supp_reads_dfs)) # condition
# ic([len(sample_supp_reads_dfs[x]) for x in range(len(sample_supp_reads_dfs))]) # assignment method per condition
# ic(len(sample_supp_reads_dfs[0][0])) # sample solutions per assignment method per condition

# sample_supp_reads_dfs[0][0][0]

# %%
# df = sample_supp_reads_dfs[0][0][0]
# 100 * df["TotalEqualSupportingReads"] / df["TotalEqualSupportingReads"].sum()

# %%
# merged_sample_supp_reads_dfs = []
# for nested_dfs in sample_supp_reads_dfs:
#     dfs = []
#     for assignment_dfs in nested_dfs:
#         for solution_df in assignment_dfs:
#             dfs.append(solution_df)
#     ic(len(dfs))
#     df = pd.concat(dfs).reset_index(drop=True)
#     merged_sample_supp_reads_dfs.append(df)
# merged_sample_supp_reads_dfs[0]

# %%
# # num of proteins with different assignment results
# gb = expression_dfs[0].groupby("#Solution")
# gb.agg({"Diff5+": ['size', 'sum']})

# %%
# sol926_exp_df = expression_dfs[0].loc[expression_dfs[0]["#Solution"] == "926"]
# sol926_exp_df

# %%
# sol926_exp_df.loc[:, ["NumOfReads", "TotalEqualSupportingReads", "TotalWeightedSupportingReads", "Diff5+"]]

# %%
# # df = sol926_exp_df.loc[
# #     (sol926_exp_df["Diff5+"]) & (sol926_exp_df["TotalWeightedSupportingReads"] < sol926_exp_df["TotalEqualSupportingReads"]),
# #     ["NumOfReads", "TotalEqualSupportingReads", "TotalWeightedSupportingReads", "Diff5+"]
# # ]
# df = sol926_exp_df.loc[
#     (sol926_exp_df["Diff5+"]) & (sol926_exp_df["TotalWeightedSupportingReads"] < sol926_exp_df["TotalEqualSupportingReads"])
# ]
# df["Equal-Weighted"] = df["TotalEqualSupportingReads"] - df["TotalWeightedSupportingReads"]
# df = df.sort_values("Equal-Weighted").reset_index(drop=True)
# df

# %%
# df.loc[481]

# %%
# sol901_exp_df = expression_dfs[0].loc[expression_dfs[0]["#Solution"] == "901"]
# sol901_exp_df

# %%
# sol901_exp_df.loc[:, ["NumOfReads", "AdditionalEqualSupportingReads", "AdditionalWeightedSupportingReads", "TotalEqualSupportingReads", "TotalWeightedSupportingReads"]].sum()

# %%
# sol901_exp_df_by_equal_supp_reads = sol901_exp_df.sort_values("TotalEqualSupportingReads").reset_index(drop=True)
# sol901_exp_df_by_equal_supp_reads

# %%
# sol901_exp_df_by_equal_supp_reads["CumTotalEqualSupportingReads"] = sol901_exp_df_by_equal_supp_reads["TotalEqualSupportingReads"].cumsum()
# sol901_exp_df_by_equal_supp_reads["%CumTotalEqualSupportingReads"] = 100 * sol901_exp_df_by_equal_supp_reads["CumTotalEqualSupportingReads"] / sol901_exp_df_by_equal_supp_reads["TotalEqualSupportingReads"].sum()
# sol901_exp_df_by_equal_supp_reads

# %%
# cummulative_supporting_reads_dfs = []
# for unique_reads_df in unique_reads_dfs:
#     df = unique_reads_df.loc[:, ["UniqueRead", "NumOfReads"]].sort_values("NumOfReads", ascending=False).reset_index(drop=True)
#     df = df.assign(CumNumOfReads = df["NumOfReads"].cumsum())
#     df["%CumNumOfReads"] = 100 * df["CumNumOfReads"] / df["NumOfReads"].sum()
#     df["NumUniqueReads"] = df.index + 1
#     cummulative_supporting_reads_dfs.append(df)
# cummulative_supporting_reads_dfs[0]

# %%
# sol901_exp_df_by_weighted_supp_reads = sol901_exp_df.sort_values("TotalWeightedSupportingReads").reset_index(drop=True)
# sol901_exp_df_by_weighted_supp_reads

# %%
# # TODO - repeat this plot with the randomly-selected solutions

# fig = px.scatter(
#     sol901_exp_df,
#     x="TotalEqualSupportingReads",
#     y="TotalWeightedSupportingReads",
#     template=template,
#     color=condition_col,
#     color_discrete_map=color_discrete_map,
#     title="Weighted vs. equal assignment of supporting reads from<br>unchosen proteins of solution 901",
#     labels={
#         "TotalEqualSupportingReads": "Total equal supporting reads",
#         "TotalWeightedSupportingReads": "Total weighted<br>supporting reads"
#     },
#     height=500,
#     width=600
# )
# fig.show()

# %%
# percentile_dfs[0]

# %%
# x_axis_name = "Distinct unique protein rank"
# y_axis_name = "Cummulative relative<br>expression (%)"
# head_title = "Cummulative expression vs. distinct unique proteins"

# cols = min(facet_col_wrap, len(conditions), 4)
# rows = ceil(len(conditions) / cols)
# row_col_iter = list(product(range(1, rows + 1), range(1, cols + 1)))[: len(conditions)]

# assignment_methods = ["Equal", "Weighted"]
# symbols = ["circle", "triangle-up"]

# fig = make_subplots(
#     rows=rows,
#     cols=cols,
#     subplot_titles=conditions,
#     shared_yaxes=True,
#     x_title=x_axis_name,
#     y_title=y_axis_name,
#     vertical_spacing=facet_row_spacing / 1.5,
#     horizontal_spacing=facet_col_spacing * 1.5,
# )

# for (row, col), percentile_df, condition in zip(
#     row_col_iter, percentile_dfs, conditions
# ):
#     legend_x = [percentile_df["RequiredProteins"].max() * 5 / 6]
#     legend_ys = [[35], [20]]

#     for color, symbol, assignment_method, legend_y in zip(
#         subcolors_discrete_map[condition], symbols, assignment_methods, legend_ys
#     ):
#         _percentile_df = percentile_df.loc[
#             percentile_df["AssignmentMethod"] == assignment_method
#         ]

#         x = _percentile_df["RequiredProteins"]
#         y = _percentile_df["Percentile"]

#         x_mean = _percentile_df.groupby("Percentile")["RequiredProteins"].apply(np.mean)
#         y_unique = x_mean.index

#         fig.add_trace(
#             go.Scatter(
#                 x=x,
#                 y=y,
#                 mode="markers",
#                 marker=dict(
#                     color=color,
#                     size=6,
#                     opacity=0.7,
#                     symbol=symbol,
#                     # line=dict(width=0),
#                 ),
#             ),
#             row=row,
#             col=col,
#         )

#         fig.add_trace(
#             go.Scatter(
#                 x=x_mean,
#                 y=y_unique,
#                 mode="lines+markers",
#                 marker=dict(
#                     color=color,
#                     size=6,
#                     opacity=0.7,
#                     symbol=symbol,
#                     # line=dict(width=0),
#                 ),
#             ),
#             row=row,
#             col=col,
#         )

#         fig.add_trace(
#             go.Scatter(
#                 x=legend_x,
#                 y=legend_y,
#                 mode="markers+text",
#                 marker=dict(
#                     color=color,
#                     size=6,
#                     opacity=0.5,
#                     symbol=symbol,
#                     # line=dict(width=0),
#                 ),
#                 text=f" {assignment_method}",
#                 textposition="middle right",
#                 textfont=dict(size=8),
#             ),
#             row=row,
#             col=col,
#         )

# fig.update_xaxes(tick0=0, dtick=5_000, matches="x")
# fig.update_layout(
#     title=head_title,
#     showlegend=False,
#     template=template,
#     height=max(300, 200 * rows),
# )
# # fig.write_image(
# #     f"{head_title} - PacBio.svg",
# #     height=max(300, 200 * rows),
# #     width=max(600, 250 * cols),
# # )
# fig.show()

# %%
x_axis_name = "Distinct protein rank"
y_axis_name = "Cummulative relative<br>expression (%)"
head_title = "Weighted cummulative expression vs. distinct protein rank"

cols = min(facet_col_wrap, len(conditions), 5)
rows = ceil(len(conditions) / cols)
row_col_iter = list(product(range(1, rows + 1), range(1, cols + 1)))[: len(conditions)]

assignment_method = "Weighted"
y_col_name = "TotalWeightedSupportingReads"

fig = make_subplots(
    rows=rows,
    cols=cols,
    subplot_titles=conditions,
    shared_yaxes=True,
    x_title=x_axis_name,
    y_title=y_axis_name,
    # vertical_spacing=facet_row_spacing / 1.5,
    # horizontal_spacing=facet_col_spacing * 1.5,
)

for (row, col), percentile_df, expression_df, maximal_solution, condition in zip(
    row_col_iter, percentile_dfs, expression_dfs, maximal_solutions, conditions
):
    percentile_df = percentile_df.loc[
        percentile_df["AssignmentMethod"] == assignment_method
    ]

    x = percentile_df["RequiredProteins"]
    y = percentile_df["Percentile"]

    # plot percentiles
    fig.add_trace(
        go.Scatter(
            x=x,
            y=y,
            mode="markers+lines",
            marker=dict(
                color=color_discrete_map[condition],
                size=6,
                opacity=0.7,
            ),
        ),
        row=row,
        col=col,
    )

    expression_df = (
        expression_df.loc[expression_df["#Solution"] == maximal_solution]
        .sort_values("TotalWeightedSupportingReads", ascending=False)
        .reset_index(drop=True)
    )
    expression_df["CummulativeRelativeWeightedExpression"] = expression_df[
        ["TotalWeightedSupportingReads"]
    ].transform(lambda x: 100 * x / x.sum())

    top_x = [10, 100, 1000]
    top_y = [
        expression_df["CummulativeRelativeWeightedExpression"][:x].sum() for x in top_x
    ]

    # plot top 10/100/1000 expressed proteins
    fig.add_trace(
        go.Scatter(
            x=top_x,
            y=top_y,
            mode="markers+text",
            marker=dict(
                size=6,
                symbol="square",
                opacity=0.5,
                color="black",
            ),
            text=[
                f"  (10, {top_y[0]:.1f})",
                f"   (100, {top_y[1]:.1f})",
                f"    (1000, {top_y[2]:.1f})",
            ],
            # ways to get better text positioning:
            # https://community.plotly.com/t/solving-the-problem-of-overlapping-text-labels-in-a-scatterplot-by-manually-assigning-the-position-of-each-label/66159/2
            # https://github.com/plotly/plotly.py/issues/925
            textposition="middle right",
            textfont=dict(size=8),
        ),
        row=row,
        col=col,
    )

# fig.update_xaxes(tick0=0, dtick=10_000, matches="x")
fig.update_layout(
    # title=head_title,
    showlegend=False,
    template=template,
    height=max(300, 260 * rows),
    width=1200,
)

fig.write_image(
    f"{head_title} - Illumina.svg",
#     height=max(300, 200 * rows),
#     width=max(600, 250 * cols),
    height=max(300, 260 * rows),
    width=1200, 
)

fig.show()

# %%
# cols = min(facet_col_wrap, len(conditions), 4)
# rows = ceil(len(conditions) / cols)
# row_col_iter = list(product(range(1, rows + 1), range(1, cols + 1)))[: len(conditions)]

# # linear_spaces = [(300, 15_000), (400, 23_000)]  # (start, end) tuples for both x and y
# # forward_transforms = [(linear_to_log10, linear_to_log10), (linear_to_log10, inverse)]  # (x, y) tuples
# # reverse_transforms = [(log10_to_linear, log10_to_linear), (log10_to_linear, inverse)]  # (x, y) tuples

# maximal_dfs = [
#     expression_df.loc[expression_df["#Solution"] == maximal_solution].reset_index(
#         drop=True
#     )
#     for expression_df, maximal_solution in zip(expression_dfs, maximal_solutions)
# ]
# maximal_algorithms = [df.loc[0, "Algorithm"] for df in maximal_dfs]

# # subplot_titles = [
# #     f"{solution} ({algorithm})"
# #     for solution, algorithm in zip(maximal_solutions, maximal_algorithms)
# # ]
# subplot_titles = [
#     # f"{condition} (#{solution}, {algorithm})"
#     f"{condition}<br><sub>(#{solution}, {algorithm})</sub>"
#     for condition, solution, algorithm in zip(
#         conditions, maximal_solutions, maximal_algorithms
#     )
# ]
# x_axis_name = "#Protein"
# y_axis_name = "Relative expression (%)"
# head_title = f"Relative expression of proteins considering a largest solution in each {str(condition_col).lower()}"

# assignment_methods = ["Equal", "Weighted"]
# y_col_names = ["TotalEqualSupportingReads", "TotalWeightedSupportingReads"]
# symbols = ["circle", "triangle-up"]

# data_marker_size = 2.5
# data_opacity = 0.2
# regression_line_width = 6

# data_scatter_type = go.Scattergl
# fit_scatter_type = go.Scatter

# fig = make_subplots(
#     rows=rows,
#     cols=cols,
#     y_title=y_axis_name,
#     x_title=x_axis_name,
#     subplot_titles=subplot_titles,
#     shared_yaxes=True,
#     shared_xaxes=True,
# )

# for (
#     (row, col),
#     condition,
#     maximal_df,
#     maximal_solution,
#     maximal_algorithm,
#     # linear_space,
#     # (forward_x_transform, forward_y_transform),
#     # (reverse_x_transform, reverse_y_transform),
# ) in zip(
#     row_col_iter,
#     conditions,
#     maximal_dfs,
#     maximal_solutions,
#     maximal_algorithms,
#     # linear_spaces,
#     # forward_transforms,
#     # reverse_transforms,
# ):

#     legend_x = [3]
#     base_legend_y = None
#     # base_legend_y = percentile_df["RequiredProteins"].max() * 5 / 6
#     # legend_ys = [[base_legend_y], [base_legend_y * 0.8]]

#     for color, assignment_method, y_col_name, symbol in zip(
#         subcolors_discrete_map[condition], assignment_methods, y_col_names, symbols
#     ):

#         assignment_df = maximal_df.sort_values(y_col_name, ascending=False).reset_index(
#             drop=True
#         )
#         assignment_df["#Protein"] = list(range(1, len(assignment_df) + 1))
#         assignment_df["AssignmentMethod"] = assignment_method

#         x = assignment_df["#Protein"]
#         y = 100 * assignment_df[y_col_name] / assignment_df[y_col_name].sum()

#         base_legend_y = (
#             y.min() if base_legend_y is None else min(base_legend_y, y.min())
#         )

#         # x = np.log10(x)
#         # y = 1 / y

#         if assignment_method == assignment_methods[0]:
#             fig.add_trace(
#                 data_scatter_type(
#                     x=x,
#                     y=y,
#                     legendgrouptitle_text=condition,
#                     legendgroup=condition,
#                     name=assignment_method,
#                     mode="markers",
#                     marker_color=color,
#                     marker_size=data_marker_size,
#                     marker=dict(
#                         opacity=data_opacity, line=dict(width=0), symbol=symbol
#                     ),
#                 ),
#                 row=row,
#                 col=col,
#             )
#         else:
#             fig.add_trace(
#                 data_scatter_type(
#                     x=x,
#                     y=y,
#                     legendgroup=condition,
#                     name=assignment_method,
#                     mode="markers",
#                     marker_color=color,
#                     marker_size=data_marker_size,
#                     marker=dict(
#                         opacity=data_opacity, line=dict(width=0), symbol=symbol
#                     ),
#                 ),
#                 row=row,
#                 col=col,
#             )

#     # legend_ys = [[base_legend_y], [base_legend_y * 0.8]]
#     legend_ys = [[base_legend_y * 10], [base_legend_y]]
#     for color, assignment_method, symbol, legend_y in zip(
#         subcolors_discrete_map[condition], assignment_methods, symbols, legend_ys
#     ):
#         fig.add_trace(
#             go.Scatter(
#                 x=legend_x,
#                 y=legend_y,
#                 mode="markers+text",
#                 marker=dict(
#                     color=color,
#                     size=6,
#                     opacity=0.5,
#                     symbol=symbol,
#                     # line=dict(width=0),
#                 ),
#                 text=f" {assignment_method}",
#                 textposition="middle right",
#                 textfont=dict(size=8),
#             ),
#             row=row,
#             col=col,
#         )

# #         train_logspace = [
# #             int(i)
# #             for i in np.logspace(
# #                 np.log10(linear_space[0]), np.log10(linear_space[1]), num=1000
# #             )
# #         ]
# #         test_logspace = [
# #             int(i)
# #             for i in np.logspace(
# #                 np.log10(linear_space[0] + 20), np.log10(linear_space[1] - 20), num=1000
# #             )
# #             if int(i) not in train_logspace
# #         ]

# #         train_x = forward_x_transform(x[train_logspace])
# #         train_y = forward_y_transform(y[train_logspace])

# #         test_x = forward_x_transform(x[test_logspace])
# #         test_y = forward_y_transform(y[test_logspace])

# #         # Create linear regression object
# #         regr = linear_model.LinearRegression(n_jobs=threads)
# #         # Train the model using the training sets
# #         regr.fit(np.array(train_x).reshape(-1, 1), train_y)
# #         # Make predictions using the testing set
# #         pred_y = regr.predict(np.array(test_x).reshape(-1, 1))

# #         # transform these variables back to original scale so they can plotted
# #         test_x = reverse_x_transform(test_x)
# #         pred_y = reverse_y_transform(pred_y)

# #         fig.add_trace(
# #             fit_scatter_type(
# #                 x=test_x,
# #                 y=pred_y,
# #                 mode="lines",
# #                 marker_color=color,
# #                 line=dict(
# #                     dash="dash",
# #                     width=regression_line_width,
# #                 ),
# #                 legendgroup=condition,
# #                 name=f"{assignment_method} - fitted",
# #                 showlegend=False
# #             ),
# #             row=1,
# #             col=col,
# #         )

# #         coef = regr.coef_[0]
# #         intercept = regr.intercept_
# #         mse = mean_squared_error(test_y, pred_y)
# #         r2 = r2_score(test_y, pred_y)
# #         if intercept >= 0:
# #             operator = "+"
# #         else:
# #             operator = "-"
# #             intercept = np.abs(intercept)

# #         if assignment_method == assignment_methods[0]:
# #             textposition = "top right"
# #             i = int(len(test_x) / 10)
# #             text_x = test_x.iloc[i] + 2000
# #             text_y = pred_y[i] + 0.02
# #         else:
# #             textposition = "bottom left"
# #             i = int(len(test_x) / 3.5)
# #             text_x = test_x.iloc[0] - int(train_logspace[0] / 2)
# #             text_y = pred_y[i] - 0.002
# #         text_x = np.log10(text_x)
# #         text_y = np.log10(text_y)

# #         text = (
# #             f"<b>y = {coef:.2f}x {operator} {intercept:.2f}</b>"
# #             "<br>"
# #             f"MSE = {mse:.2g}"  # 0 is perfect prediction
# #             "<br>"
# #             f"R2 = {r2:.2g}"  # 1 is perfect prediction
# #         )

# #         fig.add_annotation(
# #             row=1,
# #             col=col,
# #             x=text_x,
# #             y=text_y,
# #             xref="x",
# #             yref="y",
# #             text=text,
# #             align="center",
# #             font=dict(
# #                 size=9,
# #                 color=color
# #             ),
# #             showarrow=False,
# #         )

# fig.update_layout(
#     title_text=head_title,
#     title_y=0.95,
#     template=template,
#     showlegend=False,
#     # legend_itemsizing="constant",
#     height=200 * rows,
# )
# fig.update_xaxes(type="log")
# fig.update_yaxes(type="log")
# fig.write_image(f"{head_title} - Illumina.svg")
# # fig.show()
# fig.show(config={'staticPlot': True, 'responsive': False})


# %%

# %%
def linear_to_log10(arr):
    return np.log10(arr)


def log10_to_linear(log10_arr):
    return np.power([10] * len(log10_arr), log10_arr)


def inverse(arr):
    return 1 / arr


def formulate_log10_equation(coef, intercept):
    return f"y = x^{coef:.2f} * 10^{intercept:.2f}"


def formulate_semilog10_equation(coef, intercept):
    if intercept >= 0:
        operator = "+"
    else:
        operator = "-"
        intercept = np.abs(intercept)
    return f"y = 1 / ({coef:.2f}*log(x) {operator} {intercept:.2f})"


# %%
facet_row_spacing, facet_col_spacing

# %%
conditions

# %%
cols = min(facet_col_wrap, len(conditions), 5)
rows = ceil(len(conditions) / cols)
row_col_iter = list(product(range(1, rows + 1), range(1, cols + 1)))[: len(conditions)]

linear_spaces = [
    (2, 4305),
    (100, 12_386),
    (943, 5475),
    (88, 31_300),  # 3
    (103, 17_500),
    (98, 21_007),
    (140, 1142),
    (83, 5553),
    (107, 4218),
    (120, 13_000),
    (34, 11_500),
    (1044, 10_000),
    (95, 7_000),
    (192, 10_000),
    (104, 28_500),
    (99, 60_000),
    (64, 30_000),
    (254, 9354),
    (85, 19_337),
]

# forward_transforms = [
#     (linear_to_log10, linear_to_log10),
#     (linear_to_log10, linear_to_log10),
#     (linear_to_log10, linear_to_log10),
#     (linear_to_log10, linear_to_log10),  # 3
#     (linear_to_log10, linear_to_log10),
#     (linear_to_log10, linear_to_log10),
#     # (linear_to_log10, linear_to_log10),
#     (linear_to_log10, inverse),
#     (linear_to_log10, linear_to_log10),
#     (linear_to_log10, linear_to_log10),
#     (linear_to_log10, linear_to_log10),
#     (linear_to_log10, linear_to_log10),
#     (linear_to_log10, linear_to_log10),
#     (linear_to_log10, linear_to_log10),
#     (linear_to_log10, linear_to_log10),
#     (linear_to_log10, linear_to_log10),
#     (linear_to_log10, linear_to_log10),
#     (linear_to_log10, linear_to_log10),
#     (linear_to_log10, inverse),
#     (linear_to_log10, linear_to_log10),
#     # (linear_to_log10, inverse),
# ]

# reverse_transforms = [
#     (log10_to_linear, log10_to_linear),
#     (log10_to_linear, log10_to_linear),
#     (log10_to_linear, log10_to_linear),
#     (log10_to_linear, log10_to_linear),  # 3
#     (log10_to_linear, log10_to_linear),
#     (log10_to_linear, log10_to_linear),
#     # (log10_to_linear, log10_to_linear),
#     (log10_to_linear, inverse),
#     (log10_to_linear, log10_to_linear),
#     (log10_to_linear, log10_to_linear),
#     (log10_to_linear, log10_to_linear),
#     (log10_to_linear, log10_to_linear),  # 3
#     (log10_to_linear, log10_to_linear),
#     (log10_to_linear, log10_to_linear),
#     (log10_to_linear, log10_to_linear),
#     (log10_to_linear, log10_to_linear),
#     (log10_to_linear, log10_to_linear),
#     (log10_to_linear, log10_to_linear),
#     (log10_to_linear, inverse),
#     (log10_to_linear, log10_to_linear),
#     # (log10_to_linear, inverse),
# ]

forward_transforms = [(linear_to_log10, linear_to_log10)] * len(conditions)
reverse_transforms = [(log10_to_linear, log10_to_linear)] * len(conditions)
formulate_equations = [formulate_log10_equation] * len(conditions)

fit_texts = [
    "    y ~ 1 / sqrt(x)",
    "    y ~ 1 / sqrt(x)",
    "    y ~ 1 / sqrt(x)",
    "    y ~ 1 / sqrt(x)",  # 3
    "    y ~ 1 / sqrt(x)",
    "    y ~ 1 / sqrt(x)",
    # "    y ~ 1 / sqrt(x)",
    "    y ~ 1 / log(x)",
    "    y ~ 1 / sqrt(x)",
    "    y ~ 1 / sqrt(x)",
    "    y ~ 1 / sqrt(x)",
    "    y ~ 1 / sqrt(x)",
    "    y ~ 1 / sqrt(x)",
    "    y ~ 1 / sqrt(x)",
    "    y ~ 1 / sqrt(x)",
    "    y ~ 1 / sqrt(x)",
    "    y ~ 1 / sqrt(x)",
    "    y ~ 1 / sqrt(x)",
    "    y ~ 1 / log(x)",
    "    y ~ 1 / sqrt(x)",
]

maximal_dfs = [
    expression_df.loc[expression_df["#Solution"] == maximal_solution].reset_index(
        drop=True
    )
    for expression_df, maximal_solution in zip(expression_dfs, maximal_solutions)
]

assignment_dfs = [
    (
        maximal_df.sort_values("TotalWeightedSupportingReads", ascending=False)
        .reset_index(drop=True)
        .assign(ProteinRank=list(range(1, len(maximal_df) + 1)))
        .rename(columns={"ProteinRank": "#Protein"})
    )
    for maximal_df in maximal_dfs
]

subplot_titles = conditions
x_axis_name = "Distinct protein rank"
y_axis_name = "Relative expression (%)"
# head_title = f"Relative expression of proteins considering a largest solution in each {str(condition_col).lower()}"
head_title = "Relative expression vs. distinct protein rank"

data_marker_size = 2.5
data_opacity = 0.2
regression_line_width = 6

data_scatter_type = go.Scattergl
fit_scatter_type = go.Scatter

y_min = None
y_max = None

fig = make_subplots(
    rows=rows,
    cols=cols,
    y_title=y_axis_name,
    x_title=x_axis_name,
    subplot_titles=subplot_titles,
    shared_yaxes=True,
    shared_xaxes=True,
    # vertical_spacing=facet_row_spacing / 2.5,
    # horizontal_spacing=facet_col_spacing * 1.5,
    vertical_spacing=0.05,
    horizontal_spacing=0.025,
)

for (
    (row, col),
    condition,
    maximal_df,
    linear_space,
    (forward_x_transform, forward_y_transform),
    (reverse_x_transform, reverse_y_transform),
    # fit_text,
    formulate_equation,
) in zip(
    row_col_iter,
    conditions,
    maximal_dfs,
    linear_spaces,
    forward_transforms,
    reverse_transforms,
    # fit_texts,
    formulate_equations,
):

    assignment_df = maximal_df.sort_values(
        "TotalWeightedSupportingReads", ascending=False
    ).reset_index(drop=True)
    assignment_df["#Protein"] = list(range(1, len(assignment_df) + 1))

    x = assignment_df["#Protein"]
    y = (
        100
        * assignment_df["TotalWeightedSupportingReads"]
        / assignment_df["TotalWeightedSupportingReads"].sum()
    )

    y_first = y.iloc[0]
    y_last = y.iloc[-1]
    y_min = y_last if y_min is None else min(y_min, y_last)
    y_max = y_first if y_max is None else max(y_max, y_first)

    fig.add_trace(
        data_scatter_type(
            x=x,
            y=y,
            # legendgrouptitle_text=condition,
            # legendgroup=condition,
            # name=assignment_method,
            mode="markers",
            marker_color=color_discrete_map[condition],
            marker_size=data_marker_size,
            marker=dict(
                opacity=data_opacity,
                line=dict(width=0),
            ),
        ),
        row=row,
        col=col,
    )

    train_logspace = [
        int(i)
        for i in np.logspace(
            np.log10(linear_space[0]), np.log10(linear_space[1]), num=1000 if condition != conditions[6] else 400
        )
    ]
    # test_logspace = [
    #     int(i)
    #     for i in np.logspace(
    #         np.log10(linear_space[0] + 20), np.log10(linear_space[1] - 20), num=1000
    #     )
    #     if int(i) not in train_logspace
    # ]
    test_logspace = [
        int(i)
        for i in np.logspace(
            np.log10(linear_space[0]), np.log10(linear_space[1]), num=800  if condition != conditions[6] else 200
        )
        if int(i) not in train_logspace
    ]

    try:
        train_x = forward_x_transform(x[train_logspace])
        train_y = forward_y_transform(y[train_logspace])
    except:
        ic(condition)
        ic(train_logspace)
        break

    test_x = forward_x_transform(x[test_logspace])
    test_y = forward_y_transform(y[test_logspace])

    # Create linear regression object
    regr = linear_model.LinearRegression(n_jobs=threads)
    # Train the model using the training sets
    regr.fit(np.array(train_x).reshape(-1, 1), train_y)
    # Make predictions using the testing set
    pred_y = regr.predict(np.array(test_x).reshape(-1, 1))

    # transform these variables back to original scale so they can plotted
    test_x = reverse_x_transform(test_x)
    pred_y = reverse_y_transform(pred_y)

    fig.add_trace(
        fit_scatter_type(
            x=test_x,
            y=pred_y,
            mode="lines",
            marker_color="grey",
            line=dict(
                dash="dash",
                width=regression_line_width,
            ),
            showlegend=False,
        ),
        row=row,
        col=col,
    )

    i = int(len(test_x) / 10)
    text_x = test_x.iloc[i] + 2000
    text_y = pred_y[i] + 0.03
    text_x = np.log10(text_x)
    text_y = np.log10(text_y)
    
    coef = regr.coef_[0]
    intercept = regr.intercept_
    # mse = mean_squared_error(test_y, pred_y)
    # r2 = r2_score(test_y, pred_y)
    # if intercept >= 0:
    #     operator = "+"
    # else:
    #     operator = "-"
    #     intercept = np.abs(intercept)
    equation = formulate_equation(coef, intercept)
    
    # if condition == conditions[6]:
    #     ic(condition)
    #     ic(equation)
    #     ic
    
    fit_text_new = (
        # f"<b>{assignment_method}</b>"
        # "<br>"
        # f"<b>y = {coef:.2f}x {operator} {intercept:.2f}</b>"
        # f"y = {coef:.2f}x {operator} {intercept:.2f}"
        f"{equation}"
        # "<br>"
        # f"MSE = {mse:.2f}"  # 0 is perfect prediction
        # "<br>"
        # f"R2 = {r2:.2f}"  # 1 is perfect prediction
    )

    fig.add_annotation(
        row=row,
        col=col,
        x=text_x,
        y=text_y,
        xref="x",
        yref="y",
        # text=fit_text,
        text=fit_text_new,
        align="center",
        font=dict(size=9, color="grey"),
        showarrow=False,
    )

fig.update_layout(
    # title_text=head_title,
    template=template,
    showlegend=False,
    height=260 * rows,
    width=1600,
)
fig.update_xaxes(type="log")
fig.update_yaxes(type="log", range=[np.log10(y_min), np.log10(y_max)])
fig.write_image(f"{head_title} - Illumina.svg")
fig.show()
# fig.show(config={'staticPlot': True, 'responsive': False})


# %%

# %%

# %% [markdown]
# ##### Clustering

# %%
max_sol_dfs = [
    expression_df.loc[expression_df["#Solution"] == maximal_solution].reset_index(
        drop=True
    )
    for expression_df, maximal_solution in zip(expression_dfs, maximal_solutions)
]
max_sol_dfs[0]


# %%
# from kmodes.kmodes import KModes

# %%
# km = KModes(n_clusters=4, init='Huang', n_init=5, verbose=1, n_jobs=20)

# clusters = km.fit_predict(df)

# %%
# clusters

# %%
# len(clusters)

# %%
# from collections import Counter

# %%
# Counter(clusters)

# %%
# # Print the cluster centroids
# print(km.cluster_centroids_)

# %%
# # ?km

# %%
# km.labels_

# %%
# df.values.shape

# %%
# df.values.reshape(-1)

# %%
# X = df.iloc[:2, :2]
# X

# %%

# %%
# import numpy as np

# from sklearn.cluster import DBSCAN
# from sklearn import metrics
# from sklearn.datasets import make_blobs
# from sklearn.preprocessing import StandardScaler

# %%
# centers = [[1, 1], [-1, -1], [1, -1]]
# X, labels_true = make_blobs(
#     n_samples=750, centers=centers, cluster_std=0.4, random_state=0
# )

# X = StandardScaler().fit_transform(X)

# %%
# X

# %%
# labels_true

# %%
# db = DBSCAN(eps=0.3, min_samples=10).fit(X)
# core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
# core_samples_mask[db.core_sample_indices_] = True  # mark samples that were chosen as core samples
# labels = db.labels_

# %%
# labels

# %%
# # Number of clusters in labels, ignoring noise if present.
# n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
# n_noise_ = list(labels).count(-1)
# ic(n_clusters_)
# ic(n_noise_);

# %%
# unique_labels = set(labels)
# unique_labels

# %%
# k =

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%
def update_cell(
    cell,
    original_aa,
    unchanged_value,
    imputed_nan_value,
    edited_value,
    check_unchanged_func=lambda cell, original_aa: cell == original_aa,
    check_nan_func=lambda cell: "," in cell,
):
    if check_unchanged_func(cell, original_aa):
        return unchanged_value
    elif check_nan_func(cell):
        return imputed_nan_value
    else:
        return edited_value


# %%
def drop_uniformative_aa_cols(df):
    """Keep only cols with 0 and 1 values, i.e., positions where some, but not all, proteins have non-syn mutations."""

    # value_counts = df.apply(pd.unique)
    # informative_cols = value_counts.loc[value_counts.apply(lambda x: 0 in x and 1 in x)].index

    cols_editing_frequency = df.applymap(
        lambda x: x if x in [0.0, 1.0] else np.nan
    ).apply(np.mean)
    informative_cols = cols_editing_frequency.loc[cols_editing_frequency > 0].index

    return df.filter(informative_cols)


# %%
ML_INPUT_FIRST_COL_POS = 11


def prepare_ml_input_df(
    max_sol_df,
    unique_proteins_df,
    unique_proteins_first_col_pos,
    sorting_col,
    impute_nans=True,
):
    """Prepare an input df fit for some machine learning algorithms such as tSNE, PCA, etc."""

    cols_to_use_from_max_sol_df = [
        condition_col,
        "Protein",
        "#Solution",
        "Algorithm",
        "TotalEqualSupportingReads",
        "TotalWeightedSupportingReads",
        "Diff5+",
    ]
    cols_to_use_from_unique_proteins_df = [
        condition_col,
        "Protein",
    ] + unique_proteins_df.columns[unique_proteins_first_col_pos:].to_list()

    df = max_sol_df.loc[:, cols_to_use_from_max_sol_df].merge(
        unique_proteins_df.loc[:, cols_to_use_from_unique_proteins_df],
        how="left",
        on=[condition_col, "Protein"],
    )

    df_a = df.iloc[:, :7]

    original_aas = [col.split("(")[1][0] for col in df.columns[7:]]

    df_b = df.iloc[:, 7:].apply(
        lambda row: [
            update_cell(cell, original_aa, 0, np.nan, 1)
            for cell, original_aa in zip(row, original_aas)
        ],
        axis=1,
        result_type="broadcast",
    )

    if impute_nans:
        mean_col_editing_freqs_wo_nan = df_b.apply(np.mean)
        df_b = df_b.apply(
            lambda row: [
                update_cell(
                    cell,
                    original_aa,
                    0.0,
                    imputed_nan_value,
                    1.0,
                    check_unchanged_func=lambda cell, original_aa: cell == 0,
                    check_nan_func=lambda cell: np.isnan(cell),
                )
                for cell, original_aa, imputed_nan_value in zip(
                    row, original_aas, mean_col_editing_freqs_wo_nan
                )
            ],
            axis=1,
            result_type="broadcast",
        )

    df_b = drop_uniformative_aa_cols(df_b)

    df = pd.concat([df_a, df_b], axis=1)

    df.insert(
        6,
        "MeanTotalSupportingReads",
        (df["TotalEqualSupportingReads"] + df["TotalWeightedSupportingReads"]) / 2,
    )
    df.insert(
        7,
        "%EqualTotalExpression",
        100 * df["TotalEqualSupportingReads"] / df["TotalEqualSupportingReads"].sum(),
    )
    df.insert(
        8,
        "%WeightedTotalExpression",
        100
        * df["TotalWeightedSupportingReads"]
        / df["TotalWeightedSupportingReads"].sum(),
    )
    df.insert(
        9,
        "%MeanTotalExpression",
        100 * df["MeanTotalSupportingReads"] / df["MeanTotalSupportingReads"].sum(),
    )
    df = df.sort_values(sorting_col, ascending=False)

    return df


# %%
def color_highest_expressed_proteins(n, rank_cutoff, color_options=["red", "black"]):
    colors = []
    for rank in range(1, n + 1):
        if rank <= rank_cutoff:
            color = color_options[0]
        else:
            color = color_options[1]
        colors.append(color)
    return colors


# %%
def run_tsnes(
    conditions,
    tsne_input_dfs,
    seed,
    perplexities=[5, 30, 50, 100],
    n_components=2,
    learning_rate="auto",
    n_iter=300,
    init="random",
    n_jobs=20,
    first_col_pos=ML_INPUT_FIRST_COL_POS,
    top_expressed_proteins=None,
):
    conditions_tsnes = []
    conditions_Xs = []

    rng = np.random.RandomState(seed)

    for condition, tsne_input_df in zip(conditions, tsne_input_dfs):

        condition_tsnes = []

        if top_expressed_proteins is None:
            _top_expressed_proteins = len(tsne_input_df)
        else:
            _top_expressed_proteins = top_expressed_proteins
        X = tsne_input_df.iloc[:_top_expressed_proteins, first_col_pos:].values
        conditions_Xs.append(X)

        for perplexity in perplexities:

            t_sne = TSNE(
                n_components=n_components,
                learning_rate=learning_rate,
                perplexity=perplexity,
                n_iter=n_iter,
                init=init,
                random_state=rng,
                n_jobs=n_jobs,
            )

            prots_perplexity_tsne = t_sne.fit_transform(X)
            condition_tsnes.append(prots_perplexity_tsne)

        conditions_tsnes.append(condition_tsnes)

    return conditions_tsnes, conditions_Xs


# %%
def run_pcas(
    conditions,
    pca_input_dfs,
    seed,
    n_components=2,
    first_col_pos=ML_INPUT_FIRST_COL_POS,
    top_expressed_proteins=None,
):
    components_dfs = []

    rng = np.random.RandomState(seed)

    for condition, pca_input_df in zip(conditions, pca_input_dfs):

        if top_expressed_proteins is None:
            _top_expressed_proteins = len(pca_input_df)
        else:
            _top_expressed_proteins = top_expressed_proteins
        X = pca_input_df.iloc[:_top_expressed_proteins, first_col_pos:].values

        pca = PCA(n_components=n_components, random_state=rng)
        components = pca.fit_transform(X)
        components_df = pd.DataFrame(components).rename(columns={0: "PC1", 1: "PC2"})
        components_dfs.append(components_df)

    return components_dfs


# %%
# equal_exp_tsne_input_dfs = [
#     prepare_ml_input_df(
#         max_sol_df,
#         unique_proteins_df,
#         unique_proteins_first_col_pos,
#         sorting_col="%EqualTotalExpression",
#     )
#     for max_sol_df, unique_proteins_df in zip(max_sol_dfs, unique_proteins_dfs)
# ]
# equal_exp_tsne_input_dfs[0]

weighted_exp_tsne_input_dfs = [
    prepare_ml_input_df(
        max_sol_df,
        unique_proteins_df,
        unique_proteins_first_col_pos,
        sorting_col="%WeightedTotalExpression",
    )
    for max_sol_df, unique_proteins_df in zip(max_sol_dfs, unique_proteins_dfs)
]
weighted_exp_tsne_input_dfs[0]

# %%
# perplexities = [5, 30, 50, 100]
# perplexities = [5, 30, 50, 100, 150]
perplexities = [150]
n_iter = 500
n_jobs = 40

# %%
perplexities[-1]

# %%
# n_iter_500_equal_conditions_tsnes, n_iter_500_equal_conditions_Xs = run_tsnes(
#     conditions,
#     equal_exp_tsne_input_dfs,
#     seed,
#     perplexities=perplexities,
#     n_iter=n_iter,
#     n_jobs=n_jobs,
# )
n_iter_500_weighted_conditions_tsnes, n_iter_500_weighted_conditions_Xs = run_tsnes(
    conditions,
    weighted_exp_tsne_input_dfs,
    seed,
    perplexities=perplexities,
    n_iter=n_iter,
    n_jobs=n_jobs,
)

# %%
len(n_iter_500_weighted_conditions_tsnes)

# %% jupyter={"source_hidden": true}
help(fig.write_image)

# %%
cols = min(facet_col_wrap, len(conditions), 5)
rows = ceil(len(conditions) / cols)
row_col_iter = list(product(range(1, rows + 1), range(1, cols + 1)))[: len(conditions)]

fig = make_subplots(
    rows=rows,
    cols=cols,
    subplot_titles=conditions,
    # shared_yaxes=True,
    # x_title=x_title,
    # y_title=y_title,
)

# min_x = None
# max_x = 0
# max_y = 0

marker_size = 1
line_width = 0.5


for (row, col), condition, X, condition_tsnes in zip(row_col_iter, conditions, n_iter_500_weighted_conditions_Xs, n_iter_500_weighted_conditions_tsnes):

    colors = "white"

    prots_perplexity_tsne = condition_tsnes[-1]  # perplexity 150

    x, y = prots_perplexity_tsne.T

    fig.add_trace(
        go.Scattergl(
        # go.Scatter(
            x=x,
            y=y,
            mode="markers",
            marker=dict(color=colors, line_width=line_width, size=marker_size),
        ),
        row=row,
        col=col,
    )
    
    
fig.update_layout(
    # title_text=head_title,
    # title_y=0.95,
    template=template,
    showlegend=False,
    width=200*cols,
    height=200*rows,
)

fig.write_image(
    "tSNE clustering - Illumina.svg", width=200*cols, height=200*rows
)

fig.show()


# %%
# rank_cutoff = 1000

# head_title = (
#     f"t-SNEs for largest solution of each {str(condition_col).lower()} under different perplexities"
#     # "<br>"
#     # f"<sub>{rank_cutoff} highest expressed proteins are colored</sub>"
# )
# row_titles = conditions
# column_titles = [f"Perplexity = {perplexity}" for perplexity in perplexities]

# fig = make_subplots(
#     rows=len(conditions),
#     cols=len(perplexities),
#     row_titles=row_titles,
#     column_titles=column_titles,
#     # shared_yaxes=True,
#     # shared_xaxes=True
# )

# marker_size = 1
# line_width = 0.5

# for row, (condition, X, condition_tsnes) in enumerate(
#     zip(
#         conditions,
#         n_iter_500_weighted_conditions_Xs,
#         n_iter_500_weighted_conditions_tsnes,
#     ),
#     start=1,
# ):

#     # n = X.shape[0]
#     # color_options = [color_discrete_map[condition], "white"]
#     # colors = color_highest_expressed_proteins(n, rank_cutoff, color_options)
#     colors = "white"

#     for col, prots_perplexity_tsne in enumerate(condition_tsnes, start=1):

#         x, y = prots_perplexity_tsne.T

#         fig.add_trace(
#             # go.Scattergl(
#             go.Scatter(
#                 x=x,
#                 y=y,
#                 mode="markers",
#                 marker=dict(color=colors, line_width=line_width, size=marker_size),
#             ),
#             row=row,
#             col=col,
#         )

# fig.update_layout(
#     # title_text=head_title,
#     # title_y=0.95,
#     template=template,
#     showlegend=False,
#     width=1500,
#     height=300 * len(conditions),
# )

# fig.write_image(
#     "tSNE clustering - Illumina.svg", width=1500, height=300 * len(conditions)
# )

# # fig.show()

# %%
# rank_cutoff = 1000

# head_title = (
#     f"t-SNEs for largest solution of each {str(condition_col).lower()} under different perplexities"
#     # "<br>"
#     # f"<sub>{rank_cutoff} highest expressed proteins are colored</sub>"
# )
# row_titles = conditions
# column_titles = [f"Perplexity = {perplexity}" for perplexity in perplexities]

# fig = make_subplots(
#     rows=len(conditions),
#     cols=len(perplexities),
#     row_titles=row_titles,
#     column_titles=column_titles,
#     # shared_yaxes=True,
#     # shared_xaxes=True
# )

# marker_size = 1
# line_width = 0.5

# for row, (condition, X, condition_tsnes) in enumerate(
#     zip(
#         conditions,
#         n_iter_500_weighted_conditions_Xs,
#         n_iter_500_weighted_conditions_tsnes,
#     ),
#     start=1,
# ):

#     # n = X.shape[0]
#     # color_options = [color_discrete_map[condition], "white"]
#     # colors = color_highest_expressed_proteins(n, rank_cutoff, color_options)
#     colors = "white"

#     for col, prots_perplexity_tsne in enumerate(condition_tsnes, start=1):

#         x, y = prots_perplexity_tsne.T

#         fig.add_trace(
#             # go.Scattergl(
#             go.Scatter(
#                 x=x,
#                 y=y,
#                 mode="markers",
#                 marker=dict(color=colors, line_width=line_width, size=marker_size),
#             ),
#             row=row,
#             col=col,
#         )

# fig.update_layout(
#     # title_text=head_title,
#     # title_y=0.95,
#     template=template,
#     showlegend=False,
#     width=1500,
#     height=300 * len(conditions),
# )

# fig.write_image(
#     "tSNE clustering - Illumina.svg", width=1500, height=300 * len(conditions)
# )

# # fig.show()

# %%
# fig.write_image(
#     "tSNE clustering - Illumina.svg",
#     # width=1500, height=300*len(conditions)
# )

# %%
# X = weighted_exp_tsne_input_dfs[0].iloc[:1000, ML_INPUT_FIRST_COL_POS:].values

# # tsne_df = tsne_dfs[0]

# rng = np.random.RandomState(seed)

# t_sne = TSNE(
#     n_components=2,
#     learning_rate="auto",
#     perplexity=30,
#     n_iter=300,
#     init="random",
#     random_state=rng,
#     n_jobs=20
# )

# prots_perplexity_tsne = t_sne.fit_transform(X)

# x, y = prots_perplexity_tsne.T

# n = X.shape[0]
# rank_cutoff = 300
# color_options = ["red", "white"]
# colors = color_highest_expressed_proteins(n, rank_cutoff, color_options)

# fig = go.Figure(
#     data=go.Scattergl(
#         x=x,
#         y=y,
#         mode='markers',
#         marker=dict(
#             color=colors,
#             line_width=0.5
#         )
#     )
# )

# fig.update_layout(
#     width=600,
#     height=600,
#     template=template
# )

# ic(t_sne.n_iter_)
# ic(t_sne.kl_divergence_)

# fig.show()

# %% [markdown]
# > All proteins

# %%
# perplexities = [5, 30, 50, 100]
perplexities = [5, 30, 50, 100, 150]
n_iter = 500
n_jobs = 60

# %%
# TODO: uncomment

# n_iter_500_equal_conditions_tsnes, n_iter_500_equal_conditions_Xs = run_tsnes(
#     conditions,
#     equal_exp_tsne_input_dfs,
#     seed,
#     perplexities=perplexities,
#     n_iter=n_iter,
#     n_jobs=n_jobs,
# )
# n_iter_500_weighted_conditions_tsnes, n_iter_500_weighted_conditions_Xs = run_tsnes(
#     conditions,
#     weighted_exp_tsne_input_dfs,
#     seed,
#     perplexities=perplexities,
#     n_iter=n_iter,
#     n_jobs=n_jobs,
# )

# %%
# TODO: uncomment


# rank_cutoff = 300

# for conditions_tsnes, conditions_Xs, sorting_method in zip(
#     [n_iter_500_equal_conditions_tsnes, n_iter_500_weighted_conditions_tsnes],
#     [n_iter_500_equal_conditions_Xs, n_iter_500_weighted_conditions_Xs],
#     ["equal", "weighted"],
# ):

#     # head_title = (
#     #     f"t-SNEs for largest solution of each {str(condition_col).lower()} under different perplexities, sorted by % of total {sorting_method} expression"
#     #     "<br>"
#     #     f"<sub>{rank_cutoff} highest expressed proteins are colored</sub>"
#     # )
#     head_title = (
#         f"t-SNEs for largest solution of each {str(condition_col).lower()} under different perplexities"
#         "<br>"
#         f"<sub>{rank_cutoff} highest expressed proteins (according to % of total {sorting_method} expression) are colored</sub>"
#     )
#     row_titles = conditions
#     column_titles = [f"Perplexity = {perplexity}" for perplexity in perplexities]

#     fig = make_subplots(
#         rows=len(conditions),
#         cols=len(perplexities),
#         row_titles=row_titles,
#         column_titles=column_titles,
#         # shared_yaxes=True,
#         # shared_xaxes=True
#     )

#     for row, (condition, X, condition_tsnes) in enumerate(
#         zip(conditions, conditions_Xs, conditions_tsnes), start=1
#     ):

#         n = X.shape[0]
#         color_options = [color_discrete_map[condition], "white"]
#         colors = color_highest_expressed_proteins(n, rank_cutoff, color_options)

#         for col, prots_perplexity_tsne in enumerate(condition_tsnes, start=1):

#             x, y = prots_perplexity_tsne.T

#             fig.add_trace(
#                 go.Scattergl(
#                     x=x, y=y, mode="markers", marker=dict(color=colors, line_width=0.5)
#                 ),
#                 row=row,
#                 col=col,
#             )

#     fig.update_layout(
#         title_text=head_title,
#         title_y=0.95,
#         template=template,
#         showlegend=False,
#         width=1200,
#         height=300 * len(conditions),
#     )

#     fig.show()

# %%
# TODO: uncomment


# rank_cutoff = 300

# for conditions_tsnes, conditions_Xs, sorting_method in zip(
#     [n_iter_500_equal_conditions_tsnes, n_iter_500_weighted_conditions_tsnes],
#     [n_iter_500_equal_conditions_Xs, n_iter_500_weighted_conditions_Xs],
#     ["equal", "weighted"],
# ):

#     # head_title = (
#     #     f"t-SNEs for largest solution of each {str(condition_col).lower()} under different perplexities, sorted by % of total {sorting_method} expression"
#     #     "<br>"
#     #     f"<sub>{rank_cutoff} highest expressed proteins are colored</sub>"
#     # )
#     head_title = (
#         f"t-SNEs for largest solution of each {str(condition_col).lower()} under different perplexities"
#         "<br>"
#         f"<sub>{rank_cutoff} highest expressed proteins (according to % of total {sorting_method} expression) are colored</sub>"
#     )
#     row_titles = conditions
#     column_titles = [f"Perplexity = {perplexity}" for perplexity in perplexities]

#     fig = make_subplots(
#         rows=len(conditions),
#         cols=len(perplexities),
#         row_titles=row_titles,
#         column_titles=column_titles,
#         # shared_yaxes=True,
#         # shared_xaxes=True
#     )

#     for row, (condition, X, condition_tsnes) in enumerate(
#         zip(conditions, conditions_Xs, conditions_tsnes), start=1
#     ):

#         n = X.shape[0]
#         color_options = [color_discrete_map[condition], "white"]
#         colors = color_highest_expressed_proteins(n, rank_cutoff, color_options)[
#             :rank_cutoff
#         ]

#         for col, prots_perplexity_tsne in enumerate(condition_tsnes, start=1):

#             x, y = prots_perplexity_tsne.T
#             x = x[:rank_cutoff]
#             y = y[:rank_cutoff]

#             fig.add_trace(
#                 go.Scattergl(
#                     x=x, y=y, mode="markers", marker=dict(color=colors, line_width=0.5)
#                 ),
#                 row=row,
#                 col=col,
#             )

#     fig.update_layout(
#         title_text=head_title,
#         title_y=0.95,
#         template=template,
#         showlegend=False,
#         width=1200,
#         height=600,
#     )

#     fig.show()

# %%
# TODO: uncomment


# equal_conditions_pcas = run_pcas(conditions, equal_exp_tsne_input_dfs, seed)
# weighted_conditions_pcas = run_pcas(conditions, weighted_exp_tsne_input_dfs, seed)

# %%
# TODO: uncomment


# rank_cutoff = 300

# for conditions_pcas, sorting_method in zip(
#     [equal_conditions_pcas, weighted_conditions_pcas], ["equal", "weighted"]
# ):

#     head_title = (
#         f"PCAs for largest solution of each {str(condition_col).lower()}"
#         "<br>"
#         f"<sub>{rank_cutoff} highest expressed proteins (according to % of total {sorting_method} expression) are colored</sub>"
#     )
#     column_titles = conditions

#     fig = make_subplots(
#         rows=1,
#         cols=len(conditions),
#         # row_titles=row_titles,
#         column_titles=column_titles,
#         shared_yaxes=True,
#         shared_xaxes=True,
#         x_title="PC1",
#         y_title="PC2",
#     )

#     for col, (condition, condition_pca) in enumerate(
#         zip(conditions, conditions_pcas), start=1
#     ):

#         n = len(condition_pca)
#         color_options = [color_discrete_map[condition], "white"]
#         colors = color_highest_expressed_proteins(n, rank_cutoff, color_options)

#         x = condition_pca["PC1"]
#         y = condition_pca["PC2"]

#         fig.add_trace(
#             go.Scattergl(
#                 x=x, y=y, mode="markers", marker=dict(color=colors, line_width=0.5)
#             ),
#             row=1,
#             col=col,
#         )

#     fig.update_layout(
#         title_text=head_title,
#         title_y=0.95,
#         template=template,
#         showlegend=False,
#         width=800,
#         height=500,
#     )

#     fig.show()

# %% [markdown]
# > Top 1000 expressed proteins

# %%
# TODO: uncomment


# # perplexities = [5, 30, 50, 100]
# perplexities = [5, 30, 50, 100, 150, 200]

# %%
# TODO: uncomment


# top_expressed_proteins = 1000
# top_1000_equal_conditions_tsnes, top_1000_equal_conditions_Xs = run_tsnes(
#     conditions,
#     equal_exp_tsne_input_dfs,
#     seed,
#     perplexities=perplexities,
#     top_expressed_proteins=top_expressed_proteins,
# )
# top_1000_weighted_conditions_tsnes, top_1000_weighted_conditions_Xs = run_tsnes(
#     conditions,
#     weighted_exp_tsne_input_dfs,
#     seed,
#     perplexities=perplexities,
#     top_expressed_proteins=top_expressed_proteins,
# )

# %%
# TODO: uncomment


# rank_cutoff = 100

# for conditions_tsnes, conditions_Xs, sorting_method in zip(
#     [top_1000_equal_conditions_tsnes, top_1000_weighted_conditions_tsnes],
#     [top_1000_equal_conditions_Xs, top_1000_weighted_conditions_Xs],
#     ["equal", "weighted"],
# ):

#     # head_title = (
#     #     f"t-SNEs for top {top_expressed_proteins} expressed proteins in largest solution of each {str(condition_col).lower()} under different perplexities, sorted by % of total {sorting_method} expression"
#     #     "<br>"
#     #     f"<sub>{rank_cutoff} highest expressed proteins are colored</sub>"
#     # )
#     head_title = (
#         f"t-SNEs for top {top_expressed_proteins} expressed proteins in largest solution of each {str(condition_col).lower()} under different perplexities"
#         "<br>"
#         f"<sub>{rank_cutoff} highest expressed proteins (according to % of total {sorting_method} expression) are colored</sub>"
#     )

#     row_titles = conditions
#     column_titles = [f"Perplexity = {perplexity}" for perplexity in perplexities]

#     fig = make_subplots(
#         rows=len(conditions),
#         cols=len(perplexities),
#         row_titles=row_titles,
#         column_titles=column_titles,
#         # shared_yaxes=True,
#         # shared_xaxes=True
#     )

#     for row, (condition, X, condition_tsnes) in enumerate(
#         zip(conditions, conditions_Xs, conditions_tsnes), start=1
#     ):

#         n = X.shape[0]
#         color_options = [color_discrete_map[condition], "white"]
#         colors = color_highest_expressed_proteins(n, rank_cutoff, color_options)

#         for col, prots_perplexity_tsne in enumerate(condition_tsnes, start=1):

#             x, y = prots_perplexity_tsne.T

#             fig.add_trace(
#                 go.Scattergl(
#                     x=x, y=y, mode="markers", marker=dict(color=colors, line_width=0.5)
#                 ),
#                 row=row,
#                 col=col,
#             )

#     fig.update_layout(
#         title_text=head_title,
#         # title_pad_b=3,
#         title_y=0.95,
#         template=template,
#         showlegend=False,
#         width=1200,
#         height=600,
#     )

#     fig.show()

# %%
# TODO: uncomment


# top_expressed_proteins = 1000
# top_1000_equal_conditions_pcas = run_pcas(
#     conditions,
#     equal_exp_tsne_input_dfs,
#     seed,
#     top_expressed_proteins=top_expressed_proteins,
# )
# top_1000_weighted_conditions_pcas = run_pcas(
#     conditions,
#     weighted_exp_tsne_input_dfs,
#     seed,
#     top_expressed_proteins=top_expressed_proteins,
# )

# %%
# TODO: uncomment


# rank_cutoff = 100

# for conditions_pcas, sorting_method in zip(
#     [top_1000_equal_conditions_pcas, top_1000_weighted_conditions_pcas],
#     ["equal", "weighted"],
# ):

#     head_title = (
#         f"PCAs for top {top_expressed_proteins} expressed proteins in largest solution of each {str(condition_col).lower()}"
#         "<br>"
#         f"<sub>{rank_cutoff} highest expressed proteins (according to % of total {sorting_method} expression) are colored</sub>"
#     )
#     column_titles = conditions

#     fig = make_subplots(
#         rows=1,
#         cols=len(conditions),
#         # row_titles=row_titles,
#         column_titles=column_titles,
#         shared_yaxes=True,
#         shared_xaxes=True,
#         x_title="PC1",
#         y_title="PC2",
#     )

#     for col, (condition, condition_pca) in enumerate(
#         zip(conditions, conditions_pcas), start=1
#     ):

#         n = len(condition_pca)
#         color_options = [color_discrete_map[condition], "white"]
#         colors = color_highest_expressed_proteins(n, rank_cutoff, color_options)

#         x = condition_pca["PC1"]
#         y = condition_pca["PC2"]

#         fig.add_trace(
#             go.Scattergl(
#                 x=x, y=y, mode="markers", marker=dict(color=colors, line_width=0.5)
#             ),
#             row=1,
#             col=col,
#         )

#     fig.update_layout(
#         title_text=head_title,
#         title_y=0.95,
#         template=template,
#         showlegend=False,
#         width=800,
#         height=500,
#     )

#     fig.show()

# %% [markdown]
# ##### Shannon's entropy

# %%
def calc_data_entropy(max_sol_exp_df, prcnt_equal_exp_col, prcnt_weighted_exp_col):
    def _calc_data_entropy(prcnt_exp_col):
        p = max_sol_exp_df[prcnt_exp_col] / 100
        s_data = sum(-p * np.log2(p))
        return s_data

    s_data_equal_exp = _calc_data_entropy(prcnt_equal_exp_col)

    s_data_weighted_exp = _calc_data_entropy(prcnt_weighted_exp_col)

    return s_data_equal_exp, s_data_weighted_exp


# %%
def calc_hypothetical_entropy(max_sol_exp_df, first_col_pos):
    p = max_sol_exp_df.iloc[:, first_col_pos:].apply(np.mean)
    s_hypo = sum(-p * np.log2(p) - (1 - p) * np.log2((1 - p)))
    return s_hypo


# %%
def calc_entropies(
    max_sol_exp_df,
    first_col_pos,
    prcnt_equal_exp_col="%EqualTotalExpression",
    prcnt_weighted_exp_col="%WeightedTotalExpression",
):
    s_data_equal_exp, s_data_weighted_exp = calc_data_entropy(
        max_sol_exp_df, prcnt_equal_exp_col, prcnt_weighted_exp_col
    )
    s_hypo = calc_hypothetical_entropy(max_sol_exp_df, first_col_pos)
    return s_data_equal_exp, s_data_weighted_exp, s_hypo


# %%
max_sol_exp_dfs = [
    prepare_ml_input_df(
        max_sol_df,
        unique_proteins_df,
        unique_proteins_first_col_pos,
        sorting_col="%WeightedTotalExpression",
        impute_nans=False,
    )
    for max_sol_df, unique_proteins_df in zip(max_sol_dfs, unique_proteins_dfs)
]

max_sol_entropies = [
    calc_entropies(
        max_sol_exp_df,
        ML_INPUT_FIRST_COL_POS,
        "%EqualTotalExpression",
        "%WeightedTotalExpression",
    )
    for max_sol_exp_df in max_sol_exp_dfs
]

entropies_names = ["EqualExpressionData", "WeightedExpressionData", "Hypothetical"]

_conditions = []
_entropies_names = []
_entropies = []
for condition, condition_entropies in zip(conditions, max_sol_entropies):
    for entropy, entropy_name in zip(condition_entropies, entropies_names):
        _conditions.append(condition)
        _entropies_names.append(entropy_name)
        _entropies.append(entropy)

shannon_df = pd.DataFrame(
    {
        condition_col: _conditions,
        "EntropyName": _entropies_names,
        "EntropyValue": _entropies,
    }
)

shannon_df

# %%
# TODO find number of non-syn cols per mal_sol_df
max_sol_exp_dfs[0].iloc[:, ML_INPUT_FIRST_COL_POS:].shape[1]

# %%
fig = go.Figure()

colors = [color_discrete_map[condition] for condition in conditions]

non_syns_per_max_sol_exp_df = [
    str(max_sol_exp_df.iloc[:, ML_INPUT_FIRST_COL_POS:].shape[1])
    for max_sol_exp_df in max_sol_exp_dfs
]

fig.add_trace(
    go.Bar(
        x=conditions,
        y=shannon_df.loc[
            shannon_df["EntropyName"] == "WeightedExpressionData", "EntropyValue"
        ],
        marker_color=colors,
        # marker_pattern_shape="x",
        marker_pattern_shape="/",
        # name="Data",
        showlegend=False,
        # width=0.3
    )
)

fig.add_trace(
    go.Bar(
        x=conditions,
        y=shannon_df.loc[shannon_df["EntropyName"] == "Hypothetical", "EntropyValue"],
        marker_color=colors,
        marker_pattern_shape="",
        text=non_syns_per_max_sol_exp_df,
        textposition="outside",
        # name="Hypothetical",
        showlegend=False,
        # width=0.3
        textfont=dict(size=20),
    )
)


# Add single entropy traces for legend
fig.add_trace(
    go.Bar(
        x=[None],
        y=[None],
        marker=dict(
            color="white",
            pattern_shape="/",
            line=dict(
                color="grey",
            ),
        ),
        legendgroup="Observed",
        name="Observed",
        showlegend=True,
        visible='legendonly'
    )
)
fig.add_trace(
    go.Bar(
        x=[None],
        y=[None],
        marker=dict(
            color="white",
            pattern_shape="",
            line=dict(
                color="grey",
            ),
        ),
        legendgroup="Hypothetical",
        name="Hypothetical",
        showlegend=True,
        visible='legendonly'
    )
)

fig.update_yaxes(title_text="Entropy")
fig.update_xaxes(
    tickfont=dict(size=10)
)  # https://plotly.com/python/axes/#set-axis-label-rotation-and-font

# fig.update_traces(
#     marker=dict(line_color="black", line_width=0.3, pattern_fillmode="replace"),
#     # width=0.3
# )

fig.update_layout(
    # title=f"Shannon's entropy of a largest solution of each {condition_col.lower()}",
    width=1100,
    height=500,
    # showlegend=False,
    template=template,
    barmode="group",
    bargap=0.2,  # gap between bars of adjacent location coordinates.
    bargroupgap=0.15,  # gap between bars of the same location coordinate.
    legend=dict(
        yanchor="bottom",
        y=1.02,
        xanchor="right",
        x=1,
        orientation="h",
    ),
)

fig.write_image(
    "Observed vs hypothetical entropy - Illumina.svg",
    height=500,
    width=1100, 
)

fig.show()

# %% [markdown] toc-hr-collapsed=true
# ## Editing in reads

# %% [markdown]
# ### Edited sites

# %% [markdown] papermill={"duration": 0.115362, "end_time": "2022-02-01T09:42:53.418861", "exception": false, "start_time": "2022-02-01T09:42:53.303499", "status": "completed"}
# #### All reads

# %%
reads_dfs[0].head()


# %%
merged_reads_df = pd.concat(
    [
        df.loc[:, [condition_col, "Read", "EditingFrequency", "EditedPositions"]]
        for df in reads_dfs
    ]
).reset_index(drop=True)

fig = px.histogram(
    merged_reads_df,
    x="EditedPositions",
    facet_col=condition_col,
    facet_col_spacing=facet_col_spacing,
    labels={"EditedPositions": "Edited sites in read"},
    title="Distribution of edited sites in all reads",
    color=condition_col,
    color_discrete_map=color_discrete_map,
    category_orders=category_orders,
    template=template,
)

# https://stackoverflow.com/questions/58167028/single-axis-caption-in-plotly-express-facet-plot
for axis in fig.layout:
    if type(fig.layout[axis]) == go.layout.YAxis:
        fig.layout[axis].title.text = ""
fig.update_layout(showlegend=False, yaxis_title="Reads")
# fig.for_each_yaxis(lambda y: y.update(title="Transcripts"))   # for facet_row

fig.show()


# %% [markdown] papermill={"duration": 0.115362, "end_time": "2022-02-01T09:42:53.418861", "exception": false, "start_time": "2022-02-01T09:42:53.303499", "status": "completed"}
# #### Distinct unique reads

# %%
# todo update & uncomment

# cols = len(conditions)

# x_title = "Edited sites per read (mean)"
# y_title = "Distinct unique reads"
# title_text = "Distribution of edited sites per distinct unique read"

# fig = make_subplots(
#     rows=1,
#     cols=cols,
#     subplot_titles=conditions,
#     shared_yaxes=True,
#     # shared_xaxes=True,
#     x_title=x_title,
#     y_title=y_title,
# )

# x_col = "EditedPositions"
# max_x = 0

# for col, condition in zip(range(1, cols + 1), conditions):
#     df = expanded_distinct_unique_reads_df.copy()
#     df = (
#         df.loc[(df["Fraction"] == 1.0) & (df[condition_col] == condition)]
#         .groupby([condition_col, "UniqueRead"])[x_col]
#         .mean()
#         .reset_index()
#     )

#     x = df[x_col]

#     fig.add_trace(
#         go.Histogram(
#             x=x,
#             marker_color=color_discrete_map[condition],
#         ),
#         row=1,
#         col=col,
#     )

#     max_x = max(max_x, x.max())

# fig.update_layout(title_text=title_text, template=template, showlegend=False)

# min_x = 0
# max_x = max_x * 1.05
# fig.update_xaxes(range=[min_x, max_x])

# fig.show()


# %% [markdown]
# ### % editing

# %% [markdown] papermill={"duration": 0.115362, "end_time": "2022-02-01T09:42:53.418861", "exception": false, "start_time": "2022-02-01T09:42:53.303499", "status": "completed"}
# #### All reads

# %%
reads_dfs[0].head()


# %%
(
    reads_dfs[0]["EditedPositions"]
    + reads_dfs[0]["UneditedPositions"]
    + reads_dfs[0]["AmbigousPositions"]
).unique()


# %%
merged_reads_df = pd.concat(
    [
        df.loc[:, [condition_col, "Read", "EditingFrequency", "EditedPositions"]]
        for df in reads_dfs
    ]
).reset_index(drop=True)
merged_reads_df["%editing"] = merged_reads_df["EditingFrequency"] * 100

fig = px.histogram(
    merged_reads_df,
    x="%editing",
    facet_col=condition_col,
    facet_col_spacing=facet_col_spacing,
    labels={"%editing": "% editing in read"},
    title="Distribution of % editing in all reads",
    color=condition_col,
    color_discrete_map=color_discrete_map,
    category_orders=category_orders,
    template=template,
)

# https://stackoverflow.com/questions/58167028/single-axis-caption-in-plotly-express-facet-plot
for axis in fig.layout:
    if type(fig.layout[axis]) == go.layout.YAxis:
        fig.layout[axis].title.text = ""
fig.update_layout(showlegend=False, yaxis_title="Reads")
# fig.for_each_yaxis(lambda y: y.update(title="Transcripts"))   # for facet_row

fig.show()


# %% [markdown] papermill={"duration": 0.115362, "end_time": "2022-02-01T09:42:53.418861", "exception": false, "start_time": "2022-02-01T09:42:53.303499", "status": "completed"}
# #### Distinct unique edited reads

# %%
# todo update & uncomment

# expanded_distinct_unique_reads_df.head()


# %%
# todo update & uncomment

# cols = len(conditions)

# x_title = "% editing in read (mean)"
# y_title = "Distinct unique<br>edited reads"
# title_text = "Distribution of % editing per Distinct unique edited read"

# fig = make_subplots(
#     rows=1,
#     cols=cols,
#     subplot_titles=conditions,
#     shared_yaxes=True,
#     # shared_xaxes=True,
#     x_title=x_title,
#     y_title=y_title,
# )

# x_col = "EditingFrequency"
# max_x = 0

# for col, condition in zip(range(1, cols + 1), conditions):
#     df = expanded_distinct_unique_reads_df.copy()
#     df = (
#         df.loc[(df["Fraction"] == 1.0) & (df[condition_col] == condition)]
#         .groupby([condition_col, "UniqueRead"])[x_col]
#         .mean()
#         .reset_index()
#     )

#     x = df[x_col] * 100

#     fig.add_trace(
#         go.Histogram(
#             x=x,
#             marker_color=color_discrete_map[condition],
#         ),
#         row=1,
#         col=col,
#     )

#     max_x = max(max_x, x.max())

# fig.update_layout(title_text=title_text, template=template, showlegend=False)

# min_x = 0
# max_x = max_x * 1.05
# fig.update_xaxes(range=[min_x, max_x])

# fig.show()


# %% [markdown]
# ## Distribution of non-syns

# %%
cols = min(facet_col_wrap, len(conditions), 5)
rows = ceil(len(conditions) / cols)
row_col_iter = list(product(range(1, rows + 1), range(1, cols + 1)))[: len(conditions)]

x_title = "Non-syn substitutions per protein"
y_title = "Proteins"
title_text = "Distribution of min & max estimates of non-syn substitutions per protein"

fig = make_subplots(
    rows=rows,
    cols=cols,
    subplot_titles=conditions,
    shared_yaxes=True,
    x_title=x_title,
    y_title=y_title,
)

min_x = None
max_x = 0
max_y = 0

col_names = ["MinNonSyns", "MaxNonSyns"]
estimate_names = ["Min", "Max"]

for (row, col), condition, proteins_df in zip(row_col_iter, conditions, proteins_dfs):

    for i, (col_name, estimate_name) in enumerate(zip(col_names, estimate_names)):

        x = proteins_df[col_name]

        fig.add_trace(
            go.Histogram(
                x=x,
                marker_color=subcolors_discrete_map[condition][i],
                name=f"{condition}, {estimate_name}",
            ),
            row=row,
            col=col,
        )

        min_x = min(min_x, x.min()) if min_x else x.min()
        max_x = max(max_x, x.max())
        max_y = max(max_y, len(x))


for (row, col), condition in zip(row_col_iter, conditions):

    for i, (col_name, estimate_name) in enumerate(zip(col_names, estimate_names)):

        fig.add_trace(
            go.Scatter(
                x=[0.75 * max_x],
                y=[(0.13 * max_y) - (20_000 * i)],
                mode="markers+text",
                marker=dict(
                    color=subcolors_discrete_map[condition][i],
                    size=9,
                    # opacity=0.7,
                    symbol="square",
                    # line=dict(width=0),
                ),
                text=estimate_name,
                textposition="middle right",
                textfont=dict(size=9),
            ),
            row=row,
            col=col,
        )

fig.update_layout(
    template=template,
    barmode="overlay",  # Overlay both histograms
    title_text=title_text,
    title_y=0.95,
    showlegend=False,
    height=200 * rows,
)

fig.update_traces(opacity=0.75)  # Reduce opacity to see both histograms
fig.update_xaxes(range=[min_x * 0.9, max_x * 1.1])
fig.update_yaxes(range=[0, max_y * 0.2])

fig.show()


# %%
cols = min(facet_col_wrap, len(conditions), 5)
rows = ceil(len(conditions) / cols)
row_col_iter = list(product(range(1, rows + 1), range(1, cols + 1)))[: len(conditions)]

x_title = "Non-syn substitutions per protein"
y_title = "Proteins"
title_text = "Distribution of min & max estimates of non-syn substitutions per protein"

fig = make_subplots(
    rows=rows,
    cols=cols,
    subplot_titles=conditions,
    shared_yaxes=True,
    x_title=x_title,
    y_title=y_title,
)

min_x = None
max_x = 0
max_y = 0

col_names = ["MinNonSyns", "MaxNonSyns"]
estimate_names = ["Min", "Max"]

for (row, col), condition, proteins_df in zip(row_col_iter, conditions, proteins_dfs):

    for i, (col_name, estimate_name) in enumerate(zip(col_names, estimate_names)):

        x = proteins_df[col_name]

        fig.add_trace(
            go.Histogram(
                x=x,
                marker_color=subcolors_discrete_map[condition][i],
                name=f"{condition}, {estimate_name}",
            ),
            row=row,
            col=col,
        )

        min_x = min(min_x, x.min()) if min_x else x.min()
        max_x = max(max_x, x.max())
        # max_y = max(max_y, len(x))

# add mean lines + text (and also keep track of max_y)

f = fig.full_figure_for_development(warn=False)
data_traces = {}
for condition in conditions:
    for estimate_name in estimate_names:
        name = f"{condition}, {estimate_name}"
        for data in f.data:
            if data.name == name:
                data_traces[name] = data
                continue
for (row, col), condition in zip(row_col_iter, conditions):
    for estimate_name in estimate_names:
        name = f"{condition}, {estimate_name}"
        data = data_traces[name]
        x = data.x
        xbins = f.data[0].xbins
        plotbins = list(
            np.arange(
                start=xbins["start"],
                stop=xbins["end"] + xbins["size"],
                step=xbins["size"],
            )
        )
        counts, bins = np.histogram(list(x), bins=plotbins)
        max_count = max(counts)
        max_y = max(max_y, max_count)
        x_mean = np.mean(x)
        fig.add_trace(
            go.Scatter(
                x=[x_mean, x_mean, x_mean],
                y=[0, max_count, max_count * 1.1],
                mode="lines+text",
                line=dict(
                    color="white",
                    dash="dash",
                    width=4,
                ),
                text=["", "", f"{x_mean:.0f}"],
                textposition="top center",
                textfont=dict(size=10),
            ),
            row=row,
            col=col,
        )

ic(max_y)

# add legends

for (row, col), condition in zip(row_col_iter, conditions):
    for i, (col_name, estimate_name) in enumerate(zip(col_names, estimate_names)):
        fig.add_trace(
            go.Scatter(
                x=[0.75 * max_x],
                # y=[(0.13 * max_y) - (1_000 * i)],
                y=[(0.7 * max_y) - (ceil(max_y / 5) * i)],
                mode="markers+text",
                marker=dict(
                    color=subcolors_discrete_map[condition][i],
                    size=9,
                    # opacity=0.7,
                    symbol="square",
                    # line=dict(width=0),
                ),
                text=estimate_name,
                textposition="middle right",
                textfont=dict(size=9),
            ),
            row=row,
            col=col,
        )

fig.update_layout(
    template=template,
    barmode="overlay",  # Overlay both histograms
    title_text=title_text,
    title_y=0.95,
    showlegend=False,
    height=max(250 * rows, 300),
    width=max(300 * cols, 800),
)

fig.update_traces(opacity=0.75)  # Reduce opacity to see both histograms
fig.update_xaxes(range=[min_x * 0.9, max_x * 1.1])
fig.update_yaxes(range=[0, max_y * 1.3])

fig.show()


# %%
non_syns_per_read_dfs = []

for proteins_df in proteins_dfs:

    proteins_df = proteins_df.loc[
        :, [condition_col, "MinNonSyns", "MaxNonSyns", "Reads"]
    ]
    proteins_df["Read"] = proteins_df["Reads"].str.split(",")
    del proteins_df["Reads"]
    proteins_df = proteins_df.explode("Read")

    non_syns_per_read_dfs.append(proteins_df)

# %%
non_syns_per_read_df = pd.concat(non_syns_per_read_dfs, ignore_index=True)
non_syns_per_read_df

# %%
non_syns_per_read_df["MinNonSyns"].sum() / len(non_syns_per_read_df)

# %%
# Distribution of min & max estimates of non-syn substitutions per *read*


cols = min(facet_col_wrap, len(conditions), 4)
rows = ceil(len(conditions) / cols)
row_col_iter = list(product(range(1, rows + 1), range(1, cols + 1)))[: len(conditions)]

x_title = "Non-syn substitutions per read"
y_title = "Reads"
title_text = "Distribution of min & max estimates of non-syn substitutions per read"

fig = make_subplots(
    rows=rows,
    cols=cols,
    subplot_titles=conditions,
    shared_yaxes=True,
    x_title=x_title,
    y_title=y_title,
)

min_x = None
max_x = 0
max_y = 0

col_names = ["MinNonSyns", "MaxNonSyns"]
estimate_names = ["Min", "Max"]

# for (row, col), condition, proteins_df in zip(row_col_iter, conditions, proteins_dfs):
for (row, col), condition, proteins_df in zip(
    row_col_iter, conditions, non_syns_per_read_dfs
):

    # proteins_df = proteins_df.loc[:, [condition_col, "MinNonSyns", "MaxNonSyns", "Reads"]]
    # proteins_df["Reads"] = proteins_df["Reads"].str.split(",")
    # proteins_df = proteins_df.explode("Reads")

    for i, (col_name, estimate_name) in enumerate(zip(col_names, estimate_names)):

        x = proteins_df[col_name]

        fig.add_trace(
            go.Histogram(
                x=x,
                marker_color=subcolors_discrete_map[condition][i],
                name=f"{condition}, {estimate_name}",
            ),
            row=row,
            col=col,
        )

        min_x = min(min_x, x.min()) if min_x else x.min()
        max_x = max(max_x, x.max())
        # max_y = max(max_y, len(x))

# add mean lines + text (and also keep track of max_y)

f = fig.full_figure_for_development(warn=False)
data_traces = {}
for condition in conditions:
    for estimate_name in estimate_names:
        name = f"{condition}, {estimate_name}"
        for data in f.data:
            if data.name == name:
                data_traces[name] = data
                continue
for (row, col), condition in zip(row_col_iter, conditions):
    for estimate_name in estimate_names:
        name = f"{condition}, {estimate_name}"
        data = data_traces[name]
        x = data.x
        xbins = f.data[0].xbins
        plotbins = list(
            np.arange(
                start=xbins["start"],
                stop=xbins["end"] + xbins["size"],
                step=xbins["size"],
            )
        )
        counts, bins = np.histogram(list(x), bins=plotbins)
        max_count = max(counts)
        max_y = max(max_y, max_count)
        x_mean = np.mean(x)
        fig.add_trace(
            go.Scatter(
                x=[x_mean, x_mean, x_mean],
                y=[0, max_count, max_count * 1.1],
                mode="lines+text",
                line=dict(
                    color="white",
                    dash="dash",
                    width=4,
                ),
                text=["", "", f"{x_mean:.0f}"],
                textposition="top center",
                textfont=dict(size=10),
            ),
            row=row,
            col=col,
        )

ic(max_y)

# add legends

for (row, col), condition in zip(row_col_iter, conditions):
    for i, (col_name, estimate_name) in enumerate(zip(col_names, estimate_names)):
        fig.add_trace(
            go.Scatter(
                x=[0.75 * max_x],
                # y=[(0.13 * max_y) - (1_000 * i)],
                # y=[(0.7 * max_y) - (1_000 * i)],
                y=[(0.7 * max_y) - (ceil(max_y / 5) * i)],
                mode="markers+text",
                marker=dict(
                    color=subcolors_discrete_map[condition][i],
                    size=9,
                    # opacity=0.7,
                    symbol="square",
                    # line=dict(width=0),
                ),
                text=estimate_name,
                textposition="middle right",
                textfont=dict(size=9),
            ),
            row=row,
            col=col,
        )

fig.update_layout(
    template=template,
    barmode="overlay",  # Overlay both histograms
    # title_text=title_text,
    title_y=0.95,
    showlegend=False,
    height=max(200 * rows, 300),
    width=max(350 * cols, 800),
)

fig.update_traces(opacity=0.75)  # Reduce opacity to see both histograms
fig.update_xaxes(range=[min_x * 0.9, max_x * 1.1])
# fig.update_yaxes(range=[0, max_y * 1.2])
fig.update_yaxes(range=[0, max_y * 1.3])

fig.write_image(
    f"{title_text} - Illumina.svg",
    width=max(350 * cols, 800),
    height=max(200 * rows, 300),
)

# fig.show()

