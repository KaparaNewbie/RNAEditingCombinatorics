# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.0
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
    "CACB2_RABIT"
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
    "comp141158_c1_seq2"
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
]

proteins_jaccard_files = [
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141881_c0_seq3.JaccardMatrixProteins.12.07.2022-20:54:38.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141044_c0_seq2.JaccardMatrixProteins.13.07.2022-06:33:23.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp140439_c0_seq1.JaccardMatrixProteins.12.07.2022-22:51:22.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp126362_c0_seq1.JaccardMatrixProteins.15.07.2022-06:11:18.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141517_c0_seq1.JaccardMatrixProteins.14.07.2022-07:43:15.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141840_c0_seq2.JaccardMatrixProteins.13.07.2022-20:30:25.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141640_c0_seq1.JaccardMatrixProteins.12.07.2022-19:44:02.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp140987_c3_seq1.JaccardMatrixProteins.18.07.2022-07:50:43.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp140910_c2_seq1.JaccardMatrixProteins.13.07.2022-16:15:35.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp136058_c0_seq1.JaccardMatrixProteins.21.07.2022-07:57:53.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141378_c0_seq7.JaccardMatrixProteins.19.07.2022-08:12:24.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141158_c1_seq2.JaccardMatrixProteins.13.07.2022-01:54:59.csv"
]

# todo update & uncomment
# expression_files = [
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA.DistinctUniqueProteins.ExpressionLevels.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO.DistinctUniqueProteins.ExpressionLevels.csv"
# ]

# todo update & uncomment
# alg_repetitions = 5

known_sites_file = (
    "/private7/projects/Combinatorics/D.pealeii/Annotations/D.pea.EditingSites.csv"
)
samtools_path = "/home/alu/kobish/anaconda3/envs/combinatorics/bin/samtools"
threads = 20
code_dir = "/private7/projects/Combinatorics/Code"
seed = 1892

# %% [markdown] papermill={"duration": 0.029907, "end_time": "2022-02-01T09:42:43.198426", "exception": false, "start_time": "2022-02-01T09:42:43.168519", "status": "completed"} tags=[]
# # Imports

# %% papermill={"duration": 2.901153, "end_time": "2022-02-01T09:42:46.125355", "exception": false, "start_time": "2022-02-01T09:42:43.224202", "status": "completed"} tags=[]
import sys
from itertools import chain, product
from pathlib import Path
from math import ceil
from functools import reduce


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
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn import linear_model
from sklearn.metrics import mean_squared_error, r2_score
# from numpy.random import RandomState


sys.path.append(str(Path(code_dir).absolute()))
from Alignment.alignment_utils import count_reads, count_reads_in_unaligned_bam


# %% [markdown] papermill={"duration": 0.040192, "end_time": "2022-02-01T09:42:46.214429", "exception": false, "start_time": "2022-02-01T09:42:46.174237", "status": "completed"} tags=[]
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

# %% papermill={"duration": 0.054755, "end_time": "2022-02-01T09:42:46.304499", "exception": false, "start_time": "2022-02-01T09:42:46.249744", "status": "completed"} tags=[]
# plotly consts
# color_sequence = px.colors.qualitative.Pastel
# color_sequence = px.colors.qualitative.D3
# color_sequence = px.colors.qualitative.G10
color_sequence = px.colors.qualitative.Dark24
color_discrete_map = {
    condition: color for condition, color in zip(conditions, color_sequence)
}
ic(color_discrete_map);
subcolors_discrete_map = {
    condition: two_subcolors_from_hex(color_discrete_map[condition])
    for condition in conditions
}
ic(subcolors_discrete_map);
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

# %% [markdown] tags=[]
# # Jaccard functions - TODO erase?

# %%
# def jaccard_index(s1: set, s2: set):
#     """
#     Calculate Jaccard index for the sets `s1` and `s2`, defined as the size of the intersection divided by the size of the union.
#     """
#     return len(s1 & s2) / len(s1 | s2)

# %%
# def calc_jaccard_df(df: pd.DataFrame, col: str):
#     """
#     Create a Jaccard index matrix according to the different sets found at `df[col]`.

#     Notice each cell in `df[col]` actually has to be of `set` type.
#     """
#     return df[col].apply(
#         lambda x: df[col].apply(lambda y: jaccard_index(x, y))
#     ).reset_index(drop=True)

# %%
# # def np_jaccard_index(a: np.ndarray, b: np.ndarray, dropna=False):
# #     if dropna:
# #         a = a[~np.isnan(a)]
# #         b = b[~np.isnan(b)]
# #     return np.intersect1d(a, b, assume_unique=True).shape[0] / np.union1d(a, b).shape[0]

# def np_jaccard_index(a: np.ndarray, b: np.ndarray):
#     return np.intersect1d(a, b, assume_unique=True).shape[0] / np.union1d(a, b).shape[0]

# %%
# # # @jit(nopython=True)
# # @jit
# # def numba_np_jaccard_index(a: np.ndarray, b: np.ndarray, dropna=False):
# #     if dropna:
# #         a = a[~np.isnan(a)]
# #         b = b[~np.isnan(b)]
# #     return np.intersect1d(a, b, assume_unique=True).shape[0] / np.union1d(a, b).shape[0]

# # @jit(nopython=True)
# @jit
# def numba_np_jaccard_index(a: np.ndarray, b: np.ndarray):
#     return np.intersect1d(a, b, assume_unique=True).shape[0] / np.union1d(a, b).shape[0]

# %%
# # def np_calc_jaccard_df(sets_array):
# #     results = np.empty((sets_array.shape[0], sets_array.shape[0]))
# #     for i in range(sets_array.shape[0]):
# #         results[i] = np.apply_along_axis(np_jaccard_index, 1, sets_array, sets_array[i])
# #     return results

# def np_calc_jaccard_matrix(sets):
#     results = np.empty((len(sets), len(sets)))
#     for i in range(len(sets)):
#         for j in range(len(sets)):
#             if i == j:
#                 jacc_ind = 1.0
#             else:
#                 jacc_ind = np_jaccard_index(sets[i], sets[j])
#             results[i][j] = jacc_ind
#     return np.array(results)

# %%
# # @jit
# # def numba_np_calc_jaccard_df(sets_array):
# #     results = np.empty((sets_array.shape[0], sets_array.shape[0]))
# #     for i in range(sets_array.shape[0]):
# #         results[i] = np.apply_along_axis(np_jaccard_index, 1, sets_array, sets_array[i])
# #     return results


# @jit
# def numba_np_calc_jaccard_matrix(sets):
#     results = np.empty((len(sets), len(sets)))
#     for i in range(len(sets)):
#         for j in range(len(sets)):
#             if i == j:
#                 jacc_ind = 1.0
#             else:
#                 jacc_ind = np_jaccard_index(sets[i], sets[j])
#             results[i][j] = jacc_ind
#     return np.array(results)


# # @njit(parallel=True, nopython=False)
# # def parallel_numba_np_calc_jaccard_matrix(sets):
# #     results = np.empty((len(sets), len(sets)))
# #     for i in prange(len(sets)):
# #         for j in prange(len(sets)):
# #             if i == j:
# #                 jacc_ind = 1.0
# #             else:
# #                 jacc_ind = np_jaccard_index(sets[i], sets[j])
# #             results[i][j] = jacc_ind
# #     return np.array(results)

# %% tags=[]
# n = 100_000
# x = [i for i in range(0, n, 2)]
# y = [i for i in range(0, n, 7)]

# %% tags=[]
# a = np.array(x)
# b = np.array(y)

# s1 = set(x)
# s2 = set(y)

# %% tags=[]
# # %%timeit
# jaccard_index(s1, s2)

# %% tags=[]
# # %%timeit
# np_jaccard_index(a, b)

# %% tags=[]
# # %%timeit
# numba_np_jaccard_index(a, b)

# %%
# condition = conditions[0]
# df = distinct_unique_proteins_df.loc[distinct_unique_proteins_df[condition_col] == condition].reset_index(drop=True)
# n_test_rows = 10
# df = df.sample(n_test_rows).reset_index(drop=True)
# # proteins_sets_array = df["Proteins"].apply(lambda x: np.array(x.split(",")))
# proteins_sets_array = np.array(df["Proteins"].apply(lambda x: np.array(x.split(","), dtype=object)), dtype=object)
# # proteins_sets_array

# %%
# # %%timeit
# np_calc_jaccard_matrix(proteins_sets_array)

# %%
# # %%timeit
# numba_np_calc_jaccard_matrix(proteins_sets_array)

# %%

# %%
# condition = conditions[0]
# df = distinct_unique_proteins_df.loc[distinct_unique_proteins_df[condition_col] == condition].reset_index(drop=True)
# proteins_sets_array = np.array(df["Proteins"].apply(lambda x: np.array(x.split(","), dtype=object)), dtype=object)
# # proteins_sets_array

# %% [markdown] papermill={"duration": 0.040192, "end_time": "2022-02-01T09:42:46.214429", "exception": false, "start_time": "2022-02-01T09:42:46.174237", "status": "completed"} tags=[] toc-hr-collapsed=true
# # Data

# %% [markdown] papermill={"duration": 0.02598, "end_time": "2022-02-01T09:42:46.438342", "exception": false, "start_time": "2022-02-01T09:42:46.412362", "status": "completed"} tags=[]
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


# %% [markdown] papermill={"duration": 0.041741, "end_time": "2022-02-01T09:42:47.760215", "exception": false, "start_time": "2022-02-01T09:42:47.718474", "status": "completed"} tags=[]
# ## Positions

# %% tags=[]
positions_dfs = [
    pd.read_csv(position_file, sep=sep) for position_file in positions_files
]
for positions_df, condition in zip(positions_dfs, conditions):
    positions_df.insert(0, condition_col, condition)
positions_dfs[0]


# %% [markdown] papermill={"duration": 0.02598, "end_time": "2022-02-01T09:42:46.438342", "exception": false, "start_time": "2022-02-01T09:42:46.412362", "status": "completed"} tags=[]
# ## Reads

# %% [markdown] papermill={"duration": 0.02598, "end_time": "2022-02-01T09:42:46.438342", "exception": false, "start_time": "2022-02-01T09:42:46.412362", "status": "completed"} tags=[]
# ### All

# %% [markdown]
# That is, all filtered reads.

# %% papermill={"duration": 1.204258, "end_time": "2022-02-01T09:42:47.668206", "exception": false, "start_time": "2022-02-01T09:42:46.463948", "status": "completed"} tags=[]
reads_dfs = [pd.read_csv(reads_file, sep=sep) for reads_file in reads_files]
reads_dfs[0]


# %%
# edited_reads_dfs = [
#     reads_df.loc[reads_df[reads_editing_col] > 0] for reads_df in reads_dfs
# ]
# edited_reads_dfs[0]


# %% [markdown] papermill={"duration": 0.041741, "end_time": "2022-02-01T09:42:47.760215", "exception": false, "start_time": "2022-02-01T09:42:47.718474", "status": "completed"} tags=[]
# ### Unique

# %% [markdown]
# All unique reads

# %% papermill={"duration": 0.126539, "end_time": "2022-02-01T09:42:47.923363", "exception": false, "start_time": "2022-02-01T09:42:47.796824", "status": "completed"} tags=[]
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


# %% [markdown] tags=[]
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
unique_proteins_dfs[1]

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
    distinct_unique_proteins_df = pd.read_csv(
        distinct_unique_proteins_file, sep=sep
    )
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

distinct_unique_proteins_df = distinct_unique_proteins_df.sort_values([condition_col, "Fraction", "FractionRepetition", "Algorithm", "AlgorithmRepetition"]).reset_index(drop=True)

distinct_unique_proteins_df


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
            [
                df.iloc[:, :unique_proteins_first_col_pos]
                for df in unique_proteins_dfs
            ]
        ),
        on=[condition_col, "Protein"],
    )
)

expanded_distinct_unique_proteins_df


# %%
distinct_unique_proteins_df2 = (
    expanded_distinct_unique_proteins_df
    .groupby([condition_col, "Fraction", "FractionRepetition", "Algorithm", "AlgorithmRepetition"])["NumOfReads"]
    .sum()
    .reset_index()
    .rename(columns={"NumOfReads": "NumOfSupportingReads"})
    .merge(distinct_unique_proteins_df, on=[condition_col, "Fraction", "FractionRepetition", "Algorithm", "AlgorithmRepetition"])
    .assign(
        SupportingReadsPerProtein = lambda x: x["NumOfSupportingReads"] / x["NumOfProteins"],
        PercentSupportedReads = lambda x: 100 * x["NumOfSupportingReads"] / x["NumOfReads"]
    )
    .rename(columns={"PercentSupportedReads": "%SupportedReads"})
)
distinct_unique_proteins_df2

# %% [markdown]
# #### Jaccards

# %%
proteins_jaccard_dfs = []
for jaccard_file in proteins_jaccard_files:
    df = pd.read_csv(jaccard_file, sep=sep)
    proteins_jaccard_dfs.append(df)
proteins_jaccard_dfs[0]

# %% tags=[]
annotated_proteins_jaccard_dfs = []

for condition, proteins_jaccard_df in zip(conditions, proteins_jaccard_dfs):
    
    df = distinct_unique_proteins_df.loc[distinct_unique_proteins_df[condition_col] == condition].reset_index(drop=True)
    df = (
        distinct_unique_proteins_df.loc[
            distinct_unique_proteins_df[condition_col] == condition, 
            [condition_col, "Fraction", "FractionRepetition", "Algorithm", "AlgorithmRepetition"]
        ]
        .reset_index(drop=True)
    )
    df = pd.concat([df, proteins_jaccard_df], axis=1)
    index_details_dict = {
        str(i+1): f'{row["Fraction"]}-{row["FractionRepetition"]}-{row["Algorithm"]}-{row["AlgorithmRepetition"]}'
        for i, row in df.iterrows()
    }
    df = df.rename(columns=index_details_dict)
    
    annotated_proteins_jaccard_dfs.append(df)
    
annotated_proteins_jaccard_dfs[0]

# %% [markdown]
# ## Summary of data loss

# %%
# todo update & uncomment

## unaligned_reads_counts = [
##     count_reads_in_unaligned_bam(samtools_path, bam, threads)
##     for bam in unaligned_bam_files
## ]

# aligned_reads_counts = [
#     count_reads(samtools_path, bam, f"{chrom}:{start+1}-{end}", include_flags, exclude_flags, threads)
#     for bam, chrom, start, end in zip(aligned_bam_files, chroms, starts, ends)
# ]

## filtered_aligned_reads_counts = [
##     count_reads(samtools_path, bam, f"{chrom}:{start+1}-{end}", include_flags, exclude_flags, threads)
##     for bam, chrom, start, end in zip(filtered_aligned_bam_files, chroms, starts, ends)
## ]

# pileup_reads_counts = [
#     len(set(chain.from_iterable(positions_df["Reads"].str.split(","))))
#     for positions_df in positions_dfs
# ]

# unique_reads_counts = [len(unique_reads_df) for unique_reads_df in unique_reads_dfs]

# max_fraction = distinct_unique_reads_df["Fraction"].max()
# max_fraction = 1.0


# distinct_unique_reads_counts = (
#     distinct_unique_reads_df.loc[
#         distinct_unique_reads_df["Fraction"] == max_fraction
#     ]
#     .groupby(condition_col)["NumUniqueReads"]
#     .mean()
#     .round()
#     .astype(int)
# )


# distinct_unique_proteins_counts = (
#     distinct_unique_proteins_df.loc[
#         distinct_unique_proteins_df["Fraction"] == max_fraction
#     ]
#     .groupby(condition_col)["NumOfProteins"]
#     .mean()
#     .round()
#     .astype(int)
# )

# todo update & uncomment

# data_loss_df = pd.DataFrame(
#     {
#         # f"Unaligned {reads_type} reads": unaligned_reads_counts,
#         "Aligned reads (within ORF)": aligned_reads_counts,
#         # "Filtered aligned reads (within ORF)": filtered_aligned_reads_counts,
#         "Pileup reads": pileup_reads_counts,
#         "Unique reads": unique_reads_counts,
#         # "Distinct unique reads (mean)": distinct_unique_reads_counts,
#         "Distinct unique proteins (mean)": distinct_unique_proteins_counts,
#     },
#     index=conditions,
# )
# data_loss_df


# %% [markdown] papermill={"duration": 0.045853, "end_time": "2022-02-01T09:42:48.953594", "exception": false, "start_time": "2022-02-01T09:42:48.907741", "status": "completed"} tags=[]
# # Results

# %% [markdown] papermill={"duration": 0.149848, "end_time": "2022-02-01T09:43:12.800733", "exception": false, "start_time": "2022-02-01T09:43:12.650885", "status": "completed"} tags=[]
# ## Data loss

# %%
# todo update & uncomment

# long_data_loss_df = pd.melt(
#     data_loss_df.reset_index().rename(columns={"index": condition_col}),
#     id_vars=[condition_col],
#     var_name="Stage",
#     value_name="Count",
# )
# long_data_loss_df = long_data_loss_df.sort_values(
#     by=[condition_col, "Count"], ascending=False
# ).reset_index(drop=True)

# max_count_per_condition = long_data_loss_df.groupby(condition_col)["Count"].max()


# def percent_count_relative_to_first_stage(condition, count):
#     max_count = max_count_per_condition[condition]
#     percent = 100 * count / max_count
#     # percent = round(percent, 1)
#     if int(percent) == percent:
#         percent = f"{percent:.0f}%"
#     else:
#         percent = f"{percent:.1f}%"
#     return percent


# long_data_loss_df["%InitialCount"] = long_data_loss_df.apply(
#     lambda x: percent_count_relative_to_first_stage(x[condition_col], x["Count"]),
#     axis=1,
# )

# long_data_loss_df["Stage"] = long_data_loss_df["Stage"].apply(
#     lambda x: x.replace(" (", "<br>(") if "(" in x else x
# )

# long_data_loss_df


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
#     height=1000,
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

# fig.show()


# %% [markdown] papermill={"duration": 0.124528, "end_time": "2022-02-01T09:43:10.054394", "exception": false, "start_time": "2022-02-01T09:43:09.929866", "status": "completed"} tags=[]
# ## Positions

# %% [markdown] papermill={"duration": 0.149848, "end_time": "2022-02-01T09:43:12.800733", "exception": false, "start_time": "2022-02-01T09:43:12.650885", "status": "completed"} tags=[]
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
cmap = sns.diverging_palette(230, 20, as_cmap=True)


# %%
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

    plt.title(
        f"Pearson correlation coefficient between editing sites in {condition} reads"
    )


# %% [markdown]
# ### Correlation between positions to sites edited in reads

# %%
reads_w_nan_dfs[0]

# %%
edited_positions_per_read_serieses = [
    reads_w_nan_df.loc[:, "EditedPositions"]
    for reads_w_nan_df in reads_w_nan_dfs
]
edited_positions_per_read_serieses[0]

# %%
editing_in_reads_dfs = [
    reads_w_nan_df.iloc[:, reads_first_col_pos:]
    for reads_w_nan_df in reads_w_nan_dfs
]
editing_in_reads_dfs[0]

# %%
positions_to_edited_sites_in_reads_correleations = []
for editing_in_reads, edited_positions_per_read in zip(editing_in_reads_dfs, edited_positions_per_read_serieses):
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

for condition, condition_corrs in zip(conditions, positions_to_edited_sites_in_reads_correleations):
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
    df = positions_to_edited_sites_in_reads_correleations_df.loc[positions_to_edited_sites_in_reads_correleations_df[condition_col] == condition]
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
fig.for_each_annotation(lambda a: a.update(text=a.text.replace(f"{condition_col}=", "")))
fig.update_layout(showlegend=False)
fig.show()

# %%
_df = positions_to_edited_sites_in_reads_correleations_df

cols = min(facet_col_wrap, len(conditions))
rows = ceil(len(conditions) / cols)
row_col_iter = list(product(range(1, rows+1), range(1, cols+1)))[:len(conditions)]

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
        col=1
     )
    
fig.update_layout(
    template=template,
    showlegend=False,
    title_text="Pearson correlation between editing sites to number of sites edited in each read",
    height=250*rows
)
fig.update_yaxes(
    zerolinewidth=zerolinewidth, 
    tickmode='linear',
    tick0=0,
    dtick=0.2
)

fig.show()

# %% [markdown]
# ### Noise in positions

# %% tags=[]
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
)
fig.update_layout(showlegend=False)
fig.update_yaxes(
    title="% noise",
    tickmode='linear',
    tick0=0,
    dtick=2
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
conditions_labels = {
    condition: ["Edited", "KnownEditing", "InProbRegion"] for condition in conditions
}

conditions_sets = {
    condition: [
        set(positions_df.loc[positions_df[label], "Position"])
        for label in conditions_labels[condition]
    ]
    for positions_df, condition in zip(positions_dfs, conditions)
}

problamatic_regions_exist = False

fig, axs = plt.subplots(
    ncols=len(conditions), figsize=(5 * len(conditions), 2.5 * len(conditions))
)
for condition, ax in zip(conditions, axs):
    labels = conditions_labels[condition]
    sets = conditions_sets[condition]
    if len(sets[2]) == 0:
        labels = labels[:2]
        sets = sets[:2]
        v_func = venn2
    else:
        v_func = venn3
        problamatic_regions_exist = True
    v_func(sets, set_labels=labels, ax=ax)
    ax.set_title(condition)

if problamatic_regions_exist:
    title = "Positions' membership: currently edited, known editing, and probalamtic regions"
else:
    title = "Positions' membership: currently edited & known editing"

fig.suptitle(title)
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

# %% papermill={"duration": 4.052404, "end_time": "2022-02-01T09:42:53.176715", "exception": false, "start_time": "2022-02-01T09:42:49.124311", "status": "completed"} tags=[]
# todo retain nan rows and turn nans to 0?

# cols = len(conditions)

cols = min(facet_col_wrap, len(conditions), 4)
rows = ceil(len(conditions) / cols)
row_col_iter = list(product(range(1, rows+1), range(1, cols+1)))[:len(conditions)]

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
    height=200*rows
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


# %% [markdown] papermill={"duration": 0.030615, "end_time": "2022-02-01T09:42:49.024262", "exception": false, "start_time": "2022-02-01T09:42:48.993647", "status": "completed"} tags=[]
# ## Num of distinct unique reads

# %%
# todo update & uncomment

# x_axis_name = "Reads"
# y_axis_name = "Distinct unique reads"
# head_title = (
#     "Distinct unique reads vs. sequencing depth"
#     "<br>"
#     # f"<sub>({alg_repetitions * 2} repetitions over each fraction of data)</sub>"
#     "<sub>(100 repetitions over each fraction of data)</sub>"
# )

# maximal_x = 0

# # Initialize figure with subplots
# fig = make_subplots(
#     rows=1, cols=1, print_grid=False, x_title=x_axis_name, y_title=y_axis_name
# )

# # Add traces
# for col, condition in enumerate(conditions, start=1):

#     df = distinct_unique_reads_df.loc[
#         distinct_unique_reads_df[condition_col] == condition
#     ]
#     x_measured = df["NumOfReads"]
#     y_measured = df["NumUniqueReads"]

#     fig.add_trace(
#         go.Scatter(
#             x=x_measured,
#             y=y_measured,
#             # mode="markers+lines",
#             mode="markers",
#             marker=dict(color=subcolors_discrete_map[condition][0], size=9),
#             legendgroup=condition,  # this can be any string
#             legendgrouptitle_text=condition,
#             name="Measured",
#         ),
#     )

#     grouped_df = df.groupby("Fraction")
#     x_fraction_mean = grouped_df["NumOfReads"].mean().reset_index()
#     y_fraction_mean = grouped_df["NumUniqueReads"].mean().reset_index()
#     mean_fraction_df = x_fraction_mean.merge(y_fraction_mean, on="Fraction")

#     fig.add_trace(
#         go.Scatter(
#             x=mean_fraction_df["NumOfReads"],
#             y=mean_fraction_df["NumUniqueReads"],
#             mode="lines",
#             marker=dict(color=subcolors_discrete_map[condition][0], size=9),
#             showlegend=False
#         ),
#     )

#     maximal_x = max(maximal_x, x_measured.max())


# fig.update_layout(
#     title_text=head_title,
#     # legend_title_text=condition_col,
#     template=template,
#     legend_font=dict(size=8),
#     legend_grouptitlefont=dict(size=9),
#     legend_tracegroupgap=4,
# )

# fig.show()


# %%
# todo update & uncomment

# distinct_unique_reads_df["NumUniqueReads"].min()


# %%
# todo update & uncomment

# y_axis_name = "Distinct unique reads"
# head_title = "Distinct unique reads vs. heuristic method"

# fig = make_subplots(
#     rows=1,
#     cols=len(conditions),
#     print_grid=False,
#     y_title=y_axis_name,
#     subplot_titles=conditions,
#     shared_yaxes=True,
# )

# algorithms = ["Ascending", "Descending"]

# for col, condition in enumerate(conditions, start=1):

#     df = distinct_unique_reads_df.loc[
#         distinct_unique_reads_df[condition_col] == condition
#     ]

#     xs = [
#         df.loc[df["Algorithm"] == algorithm, condition_col] for algorithm in algorithms
#     ]
#     ys = [
#         df.loc[df["Algorithm"] == algorithm, "NumUniqueReads"]
#         for algorithm in algorithms
#     ]

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
#         row=1,
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
#         row=1,
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
# )
# # fig.update_yaxes(range=[0, compt_unique_edited_reads_df["NumUniqueReads"].max()*1.5])
# fig.show()


# %% tags=[]
# fig = px.violin(
#     compt_unique_edited_reads_df,
#     y="NumUniqueReads",
#     facet_col="Algorithm",
#     facet_col_spacing=facet_col_spacing,
#     title=(
#         "Compatible unique reads vs. heuristic method"
#         "<br>"
#         # f"<sub>({repetitions} repetitions for each {condition_col.lower()} and algorithm)</sub>"
#     ),
#     color=condition_col,
#     color_discrete_map=color_discrete_map,
#     template=template,
#     points="all",
#     labels={"NumUniqueReads": "Compatible unique<br>edited reads"},
# )
# fig.show()


# %% [markdown] papermill={"duration": 0.030615, "end_time": "2022-02-01T09:42:49.024262", "exception": false, "start_time": "2022-02-01T09:42:48.993647", "status": "completed"} tags=[]
# ## Num of distinct unique proteins

# %%
distinct_unique_proteins_df["NumOfProteins"].max()


# %% tags=[]
# x_axis_name = "Edited reads"
# y_axis_name = "Compatible unique edited reads"
# head_title = (
#         "Compatible unique proteins vs. sequencing depth"
#         "<br>"
#         # f"<sub>({alg_repetitions * 2} repetitions over each fraction of data)</sub>"
#         "<sub>(100 repetitions over each fraction of data)</sub>"
#     )

# frac_step = 0.1
# extraploation_increment = 3

# maximal_x = 0

# # Initialize figure with subplots
# fig = make_subplots(rows=1, cols=1, print_grid=False, x_title=x_axis_name, y_title=y_axis_name)

# # Add traces
# for col, condition in enumerate(conditions, start=1):

#     df = distinct_unique_proteins_df.loc[distinct_unique_proteins_df[condition_col] == condition]

#     x_measured = df["NumOfReads"]
#     y_measured = df["NumOfProteins"]

#     fig.add_trace(go.Scatter(x=x_measured, y=y_measured, mode='markers',
#                              marker=dict(color=subcolors_discrete_map[condition][0], size=9),
#                              legendgroup=condition,  # this can be any string
#                              legendgrouptitle_text=condition,
#                              name="Measured",
#                             ),
#                   row=1, col=1)

#     x_measured = sorted(list(set(x_measured)))
#     y_measured = list(df.groupby("Fraction")["NumOfProteins"].mean())
#     # ic(condition, x_measured, y_measured)
#     f = interpolate.interp1d(x_measured, y_measured, fill_value="extrapolate")
#     max_x_measured = max(x_measured)
#     max_x_extrapolated = extraploation_increment * max_x_measured
#     maximal_x = max(maximal_x, max_x_extrapolated)
#     x_extraploated = x_measured
#     step = 1
#     new_x = max_x_measured * (1 + (frac_step * step))
#     while new_x <= max_x_extrapolated:
#         x_extraploated.append(new_x)
#         step += 1
#         new_x = max_x_measured * (1 + (frac_step * step))
#     y_extrapolated = f(x_extraploated)

#     fig.add_trace(go.Scatter(x=x_extraploated, y=y_extrapolated,
#                              mode='lines',
#                              line=dict(color=subcolors_discrete_map[condition][1],
#                                        # width=4
#                                       ),
#                              legendgroup=condition,
#                              name="Extrapolated",
#                             ),
#                   row=1, col=1)

# dscam_ys = [18_496, 19_008]
# dscam_legend_names = ["Measured", "Theoretical maximum"]
# dscam_colors = ["black", "grey"]
# dscam_legendgrouptitle_texts = ["DSCAM", None]

# fig.add_trace(go.Scatter(x=[0.025*maximal_x, 0.975*maximal_x], y=[dscam_ys[0], dscam_ys[0]],
#                                  mode='lines',
#                                  line=dict(color=dscam_colors[0], dash="dash",
#                                            # width=3
#                                           ),
#                                  legendgroup="DSCAM1",  # this can be any string
#                                  legendgrouptitle_text="DSCAM",name=dscam_legend_names[0],
#                                 ),
#                       row=1, col=1)
# fig.add_trace(go.Scatter(x=[0.025*maximal_x, 0.975*maximal_x], y=[dscam_ys[1], dscam_ys[1]],
#                                  mode='lines',
#                                  line=dict(color=dscam_colors[1], dash="dash",
#                                            # width=3
#                                           ),
#                                  legendgroup="DSCAM1",  # this can be any string
#                                  name=dscam_legend_names[1],
#                                 ),
#                       row=1, col=1)

# fig.update_layout(title_text=head_title,
#                   # legend_title_text=condition_col,
#                   template=template,
#                   legend_font=dict(size=8),
#                   legend_grouptitlefont=dict(size=9),
#                   legend_tracegroupgap=4,
#                  )

# fig.show()


# %%
condition = conditions[0]
df = distinct_unique_proteins_df.loc[
    distinct_unique_proteins_df[condition_col] == condition
]
grouped_df = df.groupby("Fraction")
x_fraction_mean = grouped_df["NumOfReads"].mean().reset_index()
x_fraction_mean


# %%
y_fraction_mean = grouped_df["NumOfProteins"].mean().reset_index()
y_fraction_mean


# %%
x_fraction_mean.merge(y_fraction_mean, on="Fraction")


# %%
x_axis_name = "Reads"
y_axis_name = "Distinct unique proteins"
head_title = (
    "Distinct unique proteins vs. sequencing depth"
    "<br>"
    # f"<sub>({alg_repetitions * 2} repetitions over each fraction of data)</sub>"
    "<sub>(100 repetitions over each fraction of data)</sub>"
)

maximal_x = 0

# Initialize figure with subplots
fig = make_subplots(
    rows=1, cols=1, print_grid=False, x_title=x_axis_name, y_title=y_axis_name
)

# Add traces
for col, condition in enumerate(conditions, start=1):

    df = distinct_unique_proteins_df.loc[
        distinct_unique_proteins_df[condition_col] == condition
    ]
    x_measured = df["NumOfReads"]
    y_measured = df["NumOfProteins"]

    fig.add_trace(
        go.Scatter(
            x=x_measured,
            y=y_measured,
            # mode="markers+lines",
            mode="markers",
            marker=dict(color=subcolors_discrete_map[condition][0], size=9),
            legendgroup=condition,  # this can be any string
            legendgrouptitle_text=condition,
            name="Measured",
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
            marker=dict(color=subcolors_discrete_map[condition][0], size=9),
            showlegend=False
        ),
    )

    maximal_x = max(maximal_x, x_measured.max())

dscam_ys = [18_496, 19_008]
dscam_legend_names = ["Measured", "Theoretical maximum"]
dscam_colors = ["black", "grey"]
dscam_legendgrouptitle_texts = ["DSCAM", None]

fig.add_trace(
    go.Scatter(
        x=[0.05 * maximal_x, 1.05 * maximal_x],
        y=[dscam_ys[0], dscam_ys[0]],
        mode="lines",
        line=dict(
            color=dscam_colors[0],
            dash="dash",
            # width=3
        ),
        legendgroup="DSCAM1",  # this can be any string
        legendgrouptitle_text="DSCAM",
        name=dscam_legend_names[0],
    ),
)
fig.add_trace(
    go.Scatter(
        x=[0.05 * maximal_x, 1.05 * maximal_x],
        y=[dscam_ys[1], dscam_ys[1]],
        mode="lines",
        line=dict(
            color=dscam_colors[1],
            dash="dash",
            # width=3
        ),
        legendgroup="DSCAM1",  # this can be any string
        name=dscam_legend_names[1],
    ),
)

fig.update_layout(
    title_text=head_title,
    # legend_title_text=condition_col,
    template=template,
    legend_font=dict(size=8),
    legend_grouptitlefont=dict(size=9),
    legend_tracegroupgap=4,
)

fig.show()


# %% [markdown]
# ### NaNs distribution

# %%
ambigous_positions_df = (
    expanded_distinct_unique_proteins_df
    .groupby([condition_col, "Fraction", "FractionRepetition", "Algorithm", "AlgorithmRepetition"])
    .agg(
        # NumOfSupportingReads=("NumOfReads", sum),
        SumOfAmbigousPositions=("AmbigousPositions", sum),
        MeanOfAmbigousPositions=("AmbigousPositions", np.mean)   
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
    ambigous_positions_df,
    x="MeanOfAmbigousPositions",
    facet_col=condition_col,
    facet_col_spacing=facet_col_spacing,
    # labels={"MeanOfAmbigousPositions": "Mean ambigous positions in proteins of a solution"},
    title="Distribution of ambigous positions in different solutions",
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
    expanded_distinct_unique_proteins_df
    .loc[expanded_distinct_unique_proteins_df["Fraction"]==1.0]
    .groupby([condition_col, "Protein"])
    .agg(
        NumOfSolutions=("Protein", "size")
    )
    .reset_index()
    .rename(columns={"NumOfSolutions": "#SolutionIncluded"})
)

df1

# %%
dfs2 = []

for condition, unique_proteins_df in zip(conditions, unique_proteins_dfs):
    df = df1.loc[df1[condition_col] == condition]
    df2 = (
        unique_proteins_df
        .loc[:, [condition_col, "Protein", "AmbigousPositions"]]
        .merge(
            df.loc[:, ["Protein", "#SolutionIncluded"]],
            how="left",
            on="Protein"
        )
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
df3.loc[df3["#SolutionIncluded"] == 0].groupby(condition_col).size()

# %%
df3.loc[df3["#SolutionIncluded"] > 0].groupby(condition_col).size()

# %%
df3.loc[df3["#SolutionIncluded"] == 100].groupby(condition_col).size()

# %%
fig = px.scatter(
    df3,
    x="AmbigousPositions",
    y="#SolutionIncluded",
    facet_col=condition_col,
    facet_col_spacing=facet_col_spacing,
    # labels={"MeanOfAmbigousPositions": "Mean ambigous positions in proteins of a solution"},
    title="#solutions as function of ambigous positions",
    color=condition_col,
    color_discrete_map=color_discrete_map,
    category_orders=category_orders,
    template=template,
)

# https://stackoverflow.com/questions/58167028/single-axis-caption-in-plotly-express-facet-plot
# for axis in fig.layout:
#     if type(fig.layout[axis]) == go.layout.YAxis:
#         fig.layout[axis].title.text = ""
# fig.update_layout(showlegend=False, yaxis_title="# solutions")
# fig.for_each_yaxis(lambda y: y.update(title="Transcripts"))   # for facet_row

fig.show()


# %%
# df3["="]

fig = px.histogram(
    df3.loc[df3["#SolutionIncluded"]==1],
    x="AmbigousPositions",
    # y="#SolutionIncluded",
    facet_col=condition_col,
    facet_col_spacing=facet_col_spacing,
    # labels={"MeanOfAmbigousPositions": "Mean ambigous positions in proteins of a solution"},
    title="#solutions as function of ambigous positions",
    color=condition_col,
    color_discrete_map=color_discrete_map,
    category_orders=category_orders,
    template=template,
)

# https://stackoverflow.com/questions/58167028/single-axis-caption-in-plotly-express-facet-plot
# for axis in fig.layout:
#     if type(fig.layout[axis]) == go.layout.YAxis:
#         fig.layout[axis].title.text = ""
# fig.update_layout(showlegend=False, yaxis_title="# solutions")
# fig.for_each_yaxis(lambda y: y.update(title="Transcripts"))   # for facet_row

fig.show()


# %% tags=[]
# # df3["="]

# fig = px.histogram(
#     df3.loc[df3["#SolutionIncluded"]==100],
#     x="AmbigousPositions",
#     # y="#SolutionIncluded",
#     facet_col=condition_col,
#     facet_col_spacing=facet_col_spacing,
#     # labels={"MeanOfAmbigousPositions": "Mean ambigous positions in proteins of a solution"},
#     title="#solutions as function of ambigous positions",
#     color=condition_col,
#     color_discrete_map=color_discrete_map,
#     category_orders=category_orders,
#     template=template,
# )

# # https://stackoverflow.com/questions/58167028/single-axis-caption-in-plotly-express-facet-plot
# # for axis in fig.layout:
# #     if type(fig.layout[axis]) == go.layout.YAxis:
# #         fig.layout[axis].title.text = ""
# # fig.update_layout(showlegend=False, yaxis_title="# solutions")
# # fig.for_each_yaxis(lambda y: y.update(title="Transcripts"))   # for facet_row

# fig.show()


# %%

# %%

# %%

# %% [markdown]
# ### Comparing algorithms

# %% [markdown]
# #### Coupled diff comparisons

# %%
distinct_unique_proteins_df

# %%
asc_df = distinct_unique_proteins_df.loc[distinct_unique_proteins_df["Algorithm"] == "Ascending"].reset_index(drop=True)
desc_df = distinct_unique_proteins_df.loc[distinct_unique_proteins_df["Algorithm"] == "Descending"].reset_index(drop=True)
desc_df

# %%
desc_df.groupby([condition_col, "Fraction", "FractionRepetition"])["NumOfProteins"].mean().reset_index()

# %%
desc_df["NumOfProteins"]

# %%
algs_diff_df = desc_df.loc[:, [condition_col, "Fraction", "FractionRepetition", "AlgorithmRepetition"]]
algs_diff_df["Desc - Asc"] = desc_df["NumOfProteins"] - asc_df["NumOfProteins"]
algs_diff_df["Diff"] = algs_diff_df["Desc - Asc"].apply(
    lambda x: "> 0" if x > 0 else "= 0" if x == 0 else "< 0"
)

algs_diff_df

# %%
algs_diff_df.loc[algs_diff_df["Diff"] != "> 0"]

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
x_axis_name = "Fraction"
y_axis_name = "#Desc <span>&#8722;</span> #Asc"
head_title = "Differences between coupled solutions' sizes"

basic_diff_names = ["Diff < 0", "Diff = 0", "Diff > 0"]
subplot_titles = []
first_row = True
for condition in conditions:
    for diff in basic_diff_names:
        st = f"<sub>{condition}</sub>"
        if first_row:
            st = f"{diff}<br>" + st
        subplot_titles.append(st)
    first_row = False

fig = make_subplots(
    rows=len(conditions),
    cols=3,
    y_title=y_axis_name,
    x_title=x_axis_name,
    # subplot_titles=["Diff < 0", "Diff = 0", "0 < Diff"],
    subplot_titles=subplot_titles,
    shared_yaxes=True,
    shared_xaxes=True
)

algorithms = ["Ascending", "Descending"]
diffs = ["< 0", "= 0", "> 0"]
symbols = ["cross", "diamond", "circle"]

for col, (diff, symbol) in enumerate(zip(diffs, symbols), start=1):
    for row, condition in enumerate(conditions, start=1):
        df = algs_diff_df.loc[(algs_diff_df["Diff"] == diff) & (algs_diff_df[condition_col] == condition)]
        x = df["Fraction"]
        y = df["Desc - Asc"]
        if col == 1:
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=y,
                    legendgrouptitle_text=condition,
                    legendgroup=condition,
                    name=condition,
                    line_color=color_discrete_map[condition],
                    mode="markers",
                    marker={"symbol": symbol}
                ),
                row=row,
                col=col,
            )
        else:
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=y,
                    # legendgrouptitle_text=condition,
                    legendgroup=condition,
                    name=condition,
                    line_color=color_discrete_map[condition],
                    mode="markers",
                    marker={"symbol": symbol}
                ),
                row=row,
                col=col,
            )

fig.update_layout(
    title_text=head_title,
    title_y=0.95,
    template=template,
    showlegend=False
)

fig.update_xaxes(
    tick0 = 0,
    dtick = 0.1,
    tickangle = -60,
    matches='x'
)

min_y = algs_diff_df["Desc - Asc"].min()
max_y = algs_diff_df["Desc - Asc"].max()
abs_y_limit = max(abs(min_y), abs(max_y))
# min_y = min_y * 1.1 if abs(abs_y_limit // min_y // 10) <= 1 else min_y - abs(abs_y_limit // min_y // 10)
min_y = min_y * 1.1 if abs(abs_y_limit // min_y // 10) <= 1 else min_y - abs(abs_y_limit // min_y)
max_y = max_y * 1.1 if abs(abs_y_limit // max_y // 10) <= 1 else max_y + abs(abs_y_limit // max_y)

fig.update_yaxes(zerolinewidth=zerolinewidth, range=[min_y, max_y])

fig.show()


# %% [markdown]
# #### Solutions' sizes

# %%
y_axis_name = "Distinct unique proteins"
head_title = "Distinct unique proteins vs. heuristic method"

fig = make_subplots(
    rows=1,
    cols=len(conditions),
    print_grid=False,
    y_title=y_axis_name,
    subplot_titles=conditions,
    shared_yaxes=True,
)

algorithms = ["Ascending", "Descending"]

for col, condition in enumerate(conditions, start=1):

    df = distinct_unique_proteins_df.loc[
        distinct_unique_proteins_df[condition_col] == condition
    ]

    xs = [
        df.loc[df["Algorithm"] == algorithm, condition_col] for algorithm in algorithms
    ]
    ys = [
        df.loc[df["Algorithm"] == algorithm, "NumOfProteins"]
        for algorithm in algorithms
    ]

    fig.add_trace(
        go.Violin(
            x=xs[0],
            y=ys[0],
            legendgrouptitle_text=condition,
            legendgroup=condition,
            name="Ascending",
            # scalegroup="Ascending",
            side="negative",
            line_color=subcolors_discrete_map[condition][0],
            # points="all"
        ),
        row=1,
        col=col,
    )

    fig.add_trace(
        go.Violin(
            x=xs[1],
            y=ys[1],
            legendgroup=condition,
            name="Descending",
            # scalegroup="Descending",
            side="positive",
            line_color=subcolors_discrete_map[condition][1],
            # points="all"
        ),
        row=1,
        col=col,
    )

# https://stackoverflow.com/a/63221694/10249633
for ax in fig["layout"]:
    if ax[:5] == "xaxis":
        fig["layout"][ax]["tickmode"] = "array"
        fig["layout"][ax]["tickvals"] = [1]
        fig["layout"][ax]["ticktext"] = [""]

fig.update_layout(
    title_text=head_title,
    template=template,
    legend_title_text=f"{condition_col}, Algorithm",
    legend_font=dict(size=12),
    legend_grouptitlefont=dict(size=9),
    legend_tracegroupgap=4,
    # violingap=0,
    # violinmode='overlay'
)
# fig.update_yaxes(range=[0, distinct_unique_proteins_df["NumOfProteins"].max()*1.05])
fig.show()


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

# %% tags=[]
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
df = annotated_proteins_jaccard_dfs[0]
# df = df.drop("Gene", axis=1).set_index(["Fraction", "FractionRepetition", "Algorithm", "AlgorithmRepetition"])
# df = df.set_axis(df.index, axis=1)
df

# %%
fractions = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
total_len = len(df)
fraction_len = total_len / len(fractions)
middles = [int((fraction_len / 2) + (fraction_len * x)) for x in range(len(fractions))]
middles

# %%
ticks = []
for x, fraction in enumerate([fraction for fraction in fractions for _ in range(int(fraction_len))], start=1):
    if x > total_len:
        break
    tick = str(fraction) if x in middles else ""
    ticks.append(tick)

# %%
for x, tick in enumerate(ticks, start=1):
    if tick != "":
        ic(x, tick)

# %%
data = df.iloc[:, 5:].to_numpy()

# %%
data.shape

# %%
fig = px.imshow(
    # data, 
    df.iloc[:, 5:],
    labels=dict(x="Fraction", y="Fraction", color="Jaccard index"),
    x=df["Fraction"], 
    y=df["Fraction"]
)
fig.update_xaxes(side="top")
fig.show()

# %%
# import dash_bio
# https://plotly.com/python/clustergram/

# %%

# %%
import seaborn as sns; sns.set_theme()
uniform_data = df.iloc[:, 5:].to_numpy()
ax = sns.heatmap(uniform_data)

# %% [markdown]
# ##### Distribution

# %% tags=[]
melted_proteins_jaccard_dfs = []

for annotated_proteins_jaccard_df in annotated_proteins_jaccard_dfs:
    
    df = pd.melt(
        annotated_proteins_jaccard_df,
        id_vars=[condition_col, "Fraction", "FractionRepetition", "Algorithm", "AlgorithmRepetition"],
        var_name="ComparedAgainst",
        value_name="JaccardIndex"
    )
    
    melted_proteins_jaccard_dfs.append(df)
    
melted_proteins_jaccard_dfs[0]

# %% tags=[]
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

# %% tags=[]
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

# %% [markdown] papermill={"duration": 0.030615, "end_time": "2022-02-01T09:42:49.024262", "exception": false, "start_time": "2022-02-01T09:42:48.993647", "status": "completed"} tags=[]
# ## Supporting reads' coverage

# %% [markdown] papermill={"duration": 0.030615, "end_time": "2022-02-01T09:42:49.024262", "exception": false, "start_time": "2022-02-01T09:42:48.993647", "status": "completed"} tags=[]
# ### Unique reads

# %%
unique_reads_dfs[0]


# %%
s = unique_reads_dfs[0]["NumOfReads"].sum()
s

# %%
unique_reads_dfs[0]["NumOfReads"].cumsum()

# %%
100 * unique_reads_dfs[0]["NumOfReads"].cumsum() / s

# %%
max(df["NumOfReads"].max() for df in unique_reads_dfs)


# %%
df = unique_reads_dfs[0].assign(CumNumOfReads = unique_reads_dfs[0]["NumOfReads"].cumsum())
df["%CumNumOfReads"] = 100 * df["CumNumOfReads"] / df["NumOfReads"].sum()
df["UniqueReadWithDecreasingSupport"] = df.index
df

# %%
cummulative_supporting_reads_dfs = []
for unique_reads_df in unique_reads_dfs:
    df = unique_reads_df.loc[:, ["UniqueRead", "NumOfReads"]].sort_values("NumOfReads", ascending=False).reset_index(drop=True)
    df = df.assign(CumNumOfReads = df["NumOfReads"].cumsum())
    df["%CumNumOfReads"] = 100 * df["CumNumOfReads"] / df["NumOfReads"].sum()
    df["NumUniqueReads"] = df.index + 1
    cummulative_supporting_reads_dfs.append(df)
cummulative_supporting_reads_dfs[0]

# %% tags=[]
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

# %% papermill={"duration": 4.052404, "end_time": "2022-02-01T09:42:53.176715", "exception": false, "start_time": "2022-02-01T09:42:49.124311", "status": "completed"} tags=[]
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


# %% papermill={"duration": 4.052404, "end_time": "2022-02-01T09:42:53.176715", "exception": false, "start_time": "2022-02-01T09:42:49.124311", "status": "completed"} tags=[]
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


# %% papermill={"duration": 4.052404, "end_time": "2022-02-01T09:42:53.176715", "exception": false, "start_time": "2022-02-01T09:42:49.124311", "status": "completed"} tags=[]
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

# %% [markdown] papermill={"duration": 0.030615, "end_time": "2022-02-01T09:42:49.024262", "exception": false, "start_time": "2022-02-01T09:42:48.993647", "status": "completed"} tags=[]
# ### Distinct unique reads

# %% tags=[]
# todo update & uncomment

# expanded_distinct_unique_reads_df


# %%
# todo update & uncomment

# (
#     expanded_distinct_unique_reads_df.query("Fraction == 1.0")
#     .groupby([condition_col, "UniqueRead"])["NumOfReads"]
#     .mean()
#     .max()
# )


# %%
# todo update & uncomment

# cols = len(conditions)
# repetitions = 100

# x_title = "Reads (mean)"
# y_title = "log<sub>10</sub> distinct unique reads"
# title_text = "Distribution of supporting reads per distinct unique read"

# fig = make_subplots(
#     rows=1,
#     cols=cols,
#     subplot_titles=conditions,
#     shared_yaxes=True,
#     # shared_xaxes=True,
#     x_title=x_title,
#     y_title=y_title,
# )

# max_x = 0

# for col, condition in zip(range(1, cols + 1), conditions):
#     df = expanded_distinct_unique_reads_df.copy()
#     df = (
#         df.loc[(df["Fraction"] == 1.0) & (df[condition_col] == condition)]
#         .groupby([condition_col, "UniqueRead"])["NumOfReads"]
#         .mean()
#         .reset_index()
#     )

#     fig.add_trace(
#         go.Histogram(
#             x=df["NumOfReads"],
#             marker_color=color_discrete_map[condition],
#         ),
#         row=1,
#         col=col,
#     )

#     max_x = max(max_x, df["NumOfReads"].max())

# fig.update_layout(title_text=title_text, template=template, showlegend=False)

# min_x = 0
# max_x = max_x * 1.05
# fig.update_xaxes(range=[min_x, max_x])
# fig.update_yaxes(type="log")

# fig.show()


# %%

# %%

# %%
# todo update & uncomment

# cols = len(conditions)
# repetitions = 100
# filtering_algs = ["Ascending", "Descending"]

# x_title = "Reads (mean)"
# y_title = "log<sub>10</sub> compatible unique<br>edited reads"
# title_text = "Distribution of supporting reads per compatible unique edited read"

# fig = make_subplots(
#     rows=2,
#     cols=cols,
#     subplot_titles=conditions,
#     shared_yaxes=True,
#     shared_xaxes=True,
#     x_title=x_title,
#     y_title=y_title,
# )

# max_x = 0

# for col, condition in zip(range(1, cols + 1), conditions):
#     df = expanded_distinct_unique_reads_df.copy()
#     df = (
#         df.loc[(df["Fraction"] == 1.0) & (df[condition_col] == condition)]
#         .groupby([condition_col, "Algorithm", "UniqueRead"])["NumOfReads"]
#         .mean()
#         .reset_index()
#     )
#     for i, alg in enumerate(filtering_algs):

#         x = df.loc[df["Algorithm"] == alg, "NumOfReads"]

#         fig.add_trace(
#             go.Histogram(
#                 x=x,
#                 marker_color=subcolors_discrete_map[condition][i],
#                 # marker_color=color_discrete_map[condition],
#                 # histfunc="avg",
#                 name=f"{condition}, {alg}",
#             ),
#             row=i + 1,
#             col=col,
#         )

#         max_x = max(max_x, x.max())

# fig.update_layout(
#     title_text=title_text,
#     legend_title=f"{condition_col}, Algorithm",
#     template=template,
#     height=450,
# )

# min_x = 0
# max_x = max_x * 1.05
# fig.update_xaxes(range=[min_x, max_x])
# fig.update_yaxes(type="log")

# fig.show()


# %% [markdown] papermill={"duration": 0.030615, "end_time": "2022-02-01T09:42:49.024262", "exception": false, "start_time": "2022-02-01T09:42:48.993647", "status": "completed"} tags=[]
# ### Unique proteins

# %%
unique_proteins_dfs[0].head()


# %%
max(df["NumOfReads"].max() for df in unique_proteins_dfs)


# %% papermill={"duration": 4.052404, "end_time": "2022-02-01T09:42:53.176715", "exception": false, "start_time": "2022-02-01T09:42:49.124311", "status": "completed"} tags=[]
cols = len(conditions)

x_title = "Reads"
y_title = "log<sub>10</sub> unique proteins"
title_text = "Distribution of supporting reads per unique protein"

fig = make_subplots(
    rows=1,
    cols=cols,
    subplot_titles=conditions,
    shared_yaxes=True,
    shared_xaxes=True,
    x_title=x_title,
    y_title=y_title,
)

for col, condition, unique_proteins_df in zip(
    range(1, cols + 1), conditions, unique_proteins_dfs
):

    fig.add_trace(
        go.Histogram(
            x=unique_proteins_df["NumOfReads"],
            marker_color=color_discrete_map[condition],
        ),
        row=1,
        col=col,
    )

fig.update_layout(
    title_text=title_text,
    showlegend=False,
    template=template,
)

# min_x = min(0, min(df["NumOfReads"].min() for df in unique_proteins_dfs))
min_x = 0
max_x = max(df["NumOfReads"].max() for df in unique_proteins_dfs) * 1.05
fig.update_xaxes(range=[min_x, max_x])
fig.update_yaxes(type="log")

fig.show()


# %% [markdown] papermill={"duration": 0.030615, "end_time": "2022-02-01T09:42:49.024262", "exception": false, "start_time": "2022-02-01T09:42:48.993647", "status": "completed"} tags=[]
# ### Distinct unique proteins

# %% [markdown]
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

# %% papermill={"duration": 4.052404, "end_time": "2022-02-01T09:42:53.176715", "exception": false, "start_time": "2022-02-01T09:42:49.124311", "status": "completed"} tags=[]
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

# %% papermill={"duration": 4.052404, "end_time": "2022-02-01T09:42:53.176715", "exception": false, "start_time": "2022-02-01T09:42:49.124311", "status": "completed"} tags=[]
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
    expression_df["Diff5+"] = abs(expression_df["TotalEqualSupportingReads"] - expression_df["TotalWeightedSupportingReads"]) >= 0.05 * (expression_df["TotalEqualSupportingReads"] + expression_df["TotalWeightedSupportingReads"]) / 2
    expression_dfs.append(expression_df)
expression_dfs[0]


# %%
def find_rand_maximal_solution(expression_df, seed):
    df = expression_df.groupby("#Solution").agg("size").reset_index().rename(columns={0: "Size"})
    # rand_maximal_solution = df.loc[df["Size"] == df["Size"].max(), "#Solution"].sample(random_state=seed).reset_index(drop=True)
    rand_maximal_solution = df.loc[df["Size"] == df["Size"].max(), "#Solution"].sample(random_state=seed).values[0]
    return rand_maximal_solution


# %%
def make_percentile_df(expression_df):
    
    gb = expression_df.groupby("#Solution")    
    solutions_expression_dfs = [gb.get_group(x) for x in gb.groups]
    
    equal_supp_reads_dfs = []
    weighted_supp_reads_dfs = []
    for df in solutions_expression_dfs:
        equal_df = df.sort_values("TotalEqualSupportingReads", ascending=False).reset_index(drop=True)
        equal_df["CumTotalEqualSupportingReads"] = equal_df["TotalEqualSupportingReads"].cumsum()
        equal_df["%CumTotalEqualSupportingReads"] = 100 * equal_df["CumTotalEqualSupportingReads"] / equal_df["TotalEqualSupportingReads"].sum()
        equal_supp_reads_dfs.append(equal_df)
        weighted_df = df.sort_values("TotalWeightedSupportingReads", ascending=False).reset_index(drop=True)
        weighted_df["CumTotalWeightedSupportingReads"] = weighted_df["TotalWeightedSupportingReads"].cumsum()
        weighted_df["%CumTotalWeightedSupportingReads"] = 100 * weighted_df["CumTotalWeightedSupportingReads"] / weighted_df["TotalWeightedSupportingReads"].sum()
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
        ["%CumTotalEqualSupportingReads", "%CumTotalWeightedSupportingReads"]
    ):
        for df in dfs:
            # ic(df.iloc[:1, :3], method, col)
            # break
            solution = df.loc[0, "#Solution"]
            algorithm = df.loc[0, "Algorithm"]
            _condition = df.loc[0, condition_col]
            a = df[col].to_numpy()
            # for percentile in range(50, 100, 10):
            for percentile in range(10, 110, 10):
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
            condition_col: _conditions
        }
    )
    
    return percentile_df
    


# %%
def choose_sample_solutions(expression_df, seed):
    return (
        expression_df
        .loc[:, [condition_col, "#Solution", "Fraction", "FractionRepetition", "Algorithm", "AlgorithmRepetition"]]
        .groupby(["Algorithm", "#Solution"])
        .sample()
        .groupby("Algorithm")
        .sample(3, random_state=seed)
        .reset_index(drop=True)
        ["#Solution"]
    )


# %%
maximal_solutions = [find_rand_maximal_solution(expression_df, seed) for expression_df in expression_dfs]
maximal_solutions

# %%
percentile_dfs = [
    make_percentile_df(expression_df)
    for expression_df in expression_dfs
]
percentile_dfs[0]

# %%
all_conditions_sample_solutions = [choose_sample_solutions(expression_df, seed) for expression_df in expression_dfs]
all_conditions_sample_solutions

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
# num of proteins with different assignment results
gb = expression_dfs[0].groupby("#Solution")    
gb.agg({"Diff5+": ['size', 'sum']})

# %% tags=[]
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
sol901_exp_df = expression_dfs[0].loc[expression_dfs[0]["#Solution"] == "901"]
sol901_exp_df

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
# TODO - repeat this plot with the randomly-selected solutions

fig = px.scatter(
    sol901_exp_df,
    x="TotalEqualSupportingReads",
    y="TotalWeightedSupportingReads",
    template=template,
    color=condition_col,
    color_discrete_map=color_discrete_map,
    title="Weighted vs. equal assignment of supporting reads from<br>unchosen proteins of solution 901",
    labels={
        "TotalEqualSupportingReads": "Total equal supporting reads",
        "TotalWeightedSupportingReads": "Total weighted<br>supporting reads"
    },
    height=500,
    width=600
)
fig.show()

# %%
# percentile_dfs[0]

# %%
for percentile_df, condition in zip(percentile_dfs, conditions):
    
    # num_of_solutions = expression_df["#Solution"].nunique()

    fig = px.box(
        percentile_df,
        x="Percentile",
        y="RequiredProteins",
        template=template,
        title=(
            f"Distinct unique proteins vs. cummulative expression in {condition}"
            # "<br>"
            # "<sub>(proteins are sorted in descending order of supporting reads)</sub>"
        ),
        labels={
            "Percentile": "Cummulative relative expression (%)",
            "RequiredProteins": "Distinct unique<br>proteins"
        },
        color="AssignmentMethod",
        color_discrete_sequence=subcolors_discrete_map[condition],
        facet_col="Algorithm",
        # facet_row_spacing=facet_col_spacing*3,
        facet_col_spacing=facet_col_spacing*2,
        # height=600,
        # width=1000
    )

    fig.update_xaxes(
        tick0 = 10,
        dtick = 10,
        tickangle = -60,
        matches='x'
    )

    min_y = 0
    max_y = percentile_df["RequiredProteins"].max()
    fig.update_yaxes(range=[min_y, max_y])

    # fig.update_layout(showlegend=False)

    fig.show()

# %% tags=[]
linear_spaces = [(300, 15_000), (100, 23_000)]
log_axes = [(True, True), (True, True)]

x_axis_name = "#Protein"
y_axis_name = "Relative expression (%)"
head_title = f"Relative expression of proteins considering the largest solution in each {str(condition_col).lower()}"

maximal_dfs = [
    expression_df.loc[expression_df["#Solution"] == maximal_solution].reset_index(
        drop=True
    )
    for expression_df, maximal_solution in zip(expression_dfs, maximal_solutions)
]

maximal_algorithms = [df.loc[0, "Algorithm"] for df in maximal_dfs]

subplot_titles = [
    f"{solution} ({algorithm})"
    for solution, algorithm in zip(maximal_solutions, maximal_algorithms)
]

fig = make_subplots(
    rows=1,
    cols=len(conditions),
    y_title=y_axis_name,
    x_title=x_axis_name,
    subplot_titles=subplot_titles,
    shared_yaxes=True,
    shared_xaxes=True,
)

assignment_methods = ["Equal", "Weighted"]
y_col_names = ["TotalEqualSupportingReads", "TotalWeightedSupportingReads"]

marker_size = 2.5
# marker_size = 3

for col, (
    condition,
    maximal_df,
    maximal_solution,
    maximal_algorithm,
    linear_space,
    (log_x_axis, log_y_axis),
) in enumerate(
    zip(
        conditions,
        maximal_dfs,
        maximal_solutions,
        maximal_algorithms,
        linear_spaces,
        log_axes,
    ),
    start=1,
):

    for color, assignment_method, y_col_name in zip(
        subcolors_discrete_map[condition], assignment_methods, y_col_names
    ):

        assignment_df = maximal_df.sort_values(y_col_name, ascending=False).reset_index(
            drop=True
        )
        # assignment_df["#Protein"] = [str(x) for x in range(1, len(assignment_df) + 1)]
        assignment_df["#Protein"] = list(range(1, len(assignment_df) + 1))
        assignment_df["AssignmentMethod"] = assignment_method

        x = assignment_df["#Protein"]
        y = 100 * assignment_df[y_col_name] / assignment_df[y_col_name].sum()

        if assignment_method == assignment_methods[0]:
            fig.add_trace(
                go.Scattergl(
                    x=x,
                    y=y,
                    legendgrouptitle_text=condition,
                    legendgroup=condition,
                    name=assignment_method,
                    mode="markers",
                    marker_color=color,
                    marker_size=marker_size,
                    marker=dict(
                        opacity=0.2,
                        line=dict(
                            width=0,
                        )
                    ),
                ),
                row=1,
                col=col,
            )
        else:
            fig.add_trace(
                go.Scattergl(
                    x=x,
                    y=y,
                    legendgroup=condition,
                    name=assignment_method,
                    mode="markers",
                    marker_color=color,
                    marker_size=marker_size,
                    marker=dict(
                        opacity=0.2,
                        line=dict(
                            width=0,
                        )
                    ),
                ),
                row=1,
                col=col,
            )

        train_logspace = [
            int(i)
            for i in np.logspace(
                np.log10(linear_space[0]), np.log10(linear_space[1]), num=1000
            )
        ]
        test_logspace = [
            int(i)
            for i in np.logspace(
                np.log10(linear_space[0] + 20), np.log10(linear_space[1] - 20), num=1000
            )
            if int(i) not in train_logspace
        ]

        train_x = x[train_logspace]
        test_x = x[test_logspace]
        train_y = y[train_logspace]
        test_y = y[test_logspace]
        
        if log_x_axis:
            train_x = np.log10(train_x)
            test_x = np.log10(test_x)
        if log_y_axis:
            train_y = np.log10(train_y)
            test_y = np.log10(test_y)

        # Create linear regression object
        regr = linear_model.LinearRegression(n_jobs=threads)
        # Train the model using the training sets
        regr.fit(np.array(train_x).reshape(-1, 1), train_y)
        # Make predictions using the testing set
        pred_y = regr.predict(np.array(test_x).reshape(-1, 1))

        # transform these variables back to original scale so they can plotted
        if log_x_axis:
            test_x = np.power([10] * len(test_x), test_x)
        if log_y_axis:
            pred_y = np.power([10] * len(pred_y), pred_y)

        fig.add_trace(
            go.Scatter(
                x=test_x,
                y=pred_y,
                mode="lines",
                marker_color=color,
                line=dict(
                    dash="dash",
                    width=5,
                ),
                # legendgroup=condition,
                # name=f"{assignment_method} - fitted",
                showlegend=False
            ),
            row=1,
            col=col,
        )
        
        # todo uncomnnet / incorporate into the plot?
        # print("Coefficient of determination: %.2f" % r2_score(test_y, pred_y)) # 1 is perfect prediction
        
        if assignment_method == assignment_methods[0]:
            textposition = "top right"
            i = int(len(test_x) / 6)
            text_x = test_x.iloc[i] + 400
            text_y = pred_y[i] + 0.001
        else:
            textposition = "bottom left"
            i = int(2 * len(test_x) / 3)
            text_x = test_x.iloc[i] - 400
            text_y = pred_y[i] - 0.001
        
        coef = regr.coef_[0]
        intercept = regr.intercept_
        if intercept >= 0:
            operator = "+" 
        else:
            operator = "-"
            intercept = np.abs(intercept)
            
        text = (
            f"y = {coef:.2g}x {operator} {intercept:.2g}"
            "<br>"
            f"Mean sq. err. = {mean_squared_error(test_y, pred_y):.2g}"
        )
        
        fig.add_trace(
            go.Scatter(
                x=[text_x],
                y=[text_y],
                mode="text",
                text=text, 
                textfont=dict(
                    size=10, 
                    color=color
                ), 
                textposition=textposition,
                showlegend=False
            ),
            row=1,
            col=col,
        )

fig.update_layout(
    title_text=head_title,
    template=template,
    legend_itemsizing="constant",
)
fig.update_xaxes(
    type="log",
)
fig.update_yaxes(
    type="log",
)
fig.show()
# fig.show(config={'staticPlot': True, 'responsive': False})


# %%
y_col_name = "TotalEqualSupportingReads"
assignment_df = maximal_dfs[0].sort_values(y_col_name, ascending=False).reset_index(drop=True)
assignment_df["#Protein"] = list(range(1, len(assignment_df) + 1))
assignment_df["AssignmentMethod"] = assignment_method
assignment_df

# %%
x = assignment_df["#Protein"]
y = 100 * assignment_df[y_col_name] / assignment_df[y_col_name].sum()
y.iloc[10_000:]

# %%

# %%
x_axis_name = "#Protein"
y_axis_name = "Relative expression (%)"
head_title = "TBD ZIPF PLOT"

subplot_titles = [
    solution 
    for condition_solutions in all_conditions_sample_solutions
    for solution in condition_solutions
]


fig = make_subplots(
    rows=len(conditions),
    cols=len(all_conditions_sample_solutions[0]),
    y_title=y_axis_name,
    x_title=x_axis_name,
    subplot_titles=subplot_titles,
    shared_yaxes=True,
    shared_xaxes=True
)


assignment_methods = ["Equal", "Weighted"]
y_col_names = ["TotalEqualSupportingReads", "TotalWeightedSupportingReads"]


for row, (condition, expression_df, sample_solutions) in enumerate(zip(conditions, expression_dfs, all_conditions_sample_solutions), start=1):
    ic(row, condition)
    gb = expression_df.groupby("#Solution")
    solutions_expression_dfs = [
        gb.get_group(x).reset_index(drop=True) 
        for x in gb.groups
        if sample_solutions.str.fullmatch(x).any()
    ]
    ic(len(solutions_expression_dfs))
    for col, df in enumerate(solutions_expression_dfs, start=1):
        
        # ic(col)
        
        solution = df.loc[0, "#Solution"]
        algorithm = df.loc[0, "Algorithm"]
        
        # ic(solution, algorithm)
        
        equal_df = df.sort_values("TotalEqualSupportingReads", ascending=False).reset_index(drop=True)        
        equal_df["#Protein"] = [str(x) for x in range(1, len(equal_df) + 1)]
        equal_df["AssignmentMethod"] = "Equal"
        weighted_df = df.sort_values("TotalWeightedSupportingReads", ascending=False).reset_index(drop=True)
        weighted_df["#Protein"] = [str(x) for x in range(1, len(weighted_df) + 1)]
        weighted_df["AssignmentMethod"] = "Weighted"
        
        assignment_dfs = [equal_df, weighted_df]
          
        x = assignment_dfs[0]["#Protein"]
        ys = [
            100 * assignment_df[y_col_name] / assignment_df[y_col_name].sum()
            for assignment_df, y_col_name in zip(assignment_dfs, y_col_names)
        ]
        y = reduce(lambda a, b: (a + b), ys) / len(ys)
        
        fig.add_trace(
            go.Scattergl(
                x=x,
                y=y,
                opacity=0.7,
                legendgrouptitle_text=condition,
                legendgroup=condition,
                name=solution,
                # name=f"{solution}, {assignment_method}",
                marker_color=subcolors_discrete_map[condition][0] if algorithm == algorithms[0] else subcolors_discrete_map[condition][1],
#                         # mode="markers",
#                         # marker={"symbol": symbol}
            ),
            row=row,
            col=col,
        )
        
#         for i, (assignment_method, y_col_name, assignment_df) in enumerate(zip(assignment_methods, y_col_names, assignment_dfs)):
        
#             # ic(assignment_method, y_col_name)
#             # ic(assignment_df)
            
#             x = assignment_df["#Protein"]
#             y = 100 * assignment_df[y_col_name] / assignment_df[y_col_name].sum()
            
#             # Xs.append(x)
#             # Ys.append(y)
            
#             # ic(x)
#             # ic(y)
#             # ic(y.sum())
            
#             if col == 1 and assignment_method == assignment_methods[0]:
#                 fig.add_trace(
#                     go.Scattergl(
#                         x=x,
#                         y=y,
#                         opacity=0.7,
#                         legendgrouptitle_text=condition,
#                         legendgroup=condition,
#                         name=f"{solution}, {assignment_method}",
#                         marker_color=subcolors_discrete_map[condition][i],
# #                         # mode="markers",
# #                         # marker={"symbol": symbol}
#                     ),
#                     row=row,
#                     col=col,
#                 )
#             else:
#                 fig.add_trace(
#                     go.Scattergl(
#                         x=x,
#                         y=y,
#                         opacity=0.7,
# #                         # legendgrouptitle_text=condition,
#                         legendgroup=condition,
#                         name=f"{solution}, {assignment_method}",
#                         marker_color=subcolors_discrete_map[condition][i],
# #                         # mode="markers",
# #                         # marker={"symbol": symbol}
#                     ),
#                     row=row,
#                     col=col,
#                 )


fig.update_layout(
    title_text=head_title,
    # title_y=0.95,
    template=template,
    # showlegend=False
)

fig.update_xaxes(
    type="log",
    # tick0 = 0,
    # dtick = 0.1,
    # tickangle = -60,
    # matches='x'
)

# min_y = algs_diff_df["Desc - Asc"].min()
# max_y = algs_diff_df["Desc - Asc"].max()
# abs_y_limit = max(abs(min_y), abs(max_y))
# # min_y = min_y * 1.1 if abs(abs_y_limit // min_y // 10) <= 1 else min_y - abs(abs_y_limit // min_y // 10)
# min_y = min_y * 1.1 if abs(abs_y_limit // min_y // 10) <= 1 else min_y - abs(abs_y_limit // min_y)
# max_y = max_y * 1.1 if abs(abs_y_limit // max_y // 10) <= 1 else max_y + abs(abs_y_limit // max_y)

fig.update_yaxes(
    type="log",
    # zerolinewidth=zerolinewidth, 
    # range=[min_y, max_y]
)

fig.show()


# %% [markdown]
# ##### Clustering

# %%
max_sol_dfs = [
    expression_df.loc[expression_df["#Solution"] == maximal_solution].reset_index(drop=True)
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
    check_unchanged_func = lambda cell, original_aa: cell == original_aa,
    check_nan_func = lambda cell: "," in cell
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
    
    cols_editing_frequency = (
        df
        .applymap(lambda x: x if x in [0.0, 1.0] else np.nan)
        .apply(np.mean)
    )
    informative_cols = cols_editing_frequency.loc[cols_editing_frequency > 0].index
    
    return df.filter(informative_cols)


# %%
ML_INPUT_FIRST_COL_POS = 11


def prepare_ml_input_df(max_sol_df, unique_proteins_df, unique_proteins_first_col_pos, sorting_col, impute_nans=True):
    """Prepare an input df fit for some machine learning algorithms such as tSNE, PCA, etc."""
    
    cols_to_use_from_max_sol_df = [condition_col, "Protein", "#Solution", "Algorithm", "TotalEqualSupportingReads", "TotalWeightedSupportingReads", "Diff5+"]
    cols_to_use_from_unique_proteins_df = [condition_col, "Protein"] + unique_proteins_df.columns[unique_proteins_first_col_pos:].to_list()
    
    df = (
        max_sol_df
        .loc[:, cols_to_use_from_max_sol_df]
        .merge(
            unique_proteins_df.loc[:, cols_to_use_from_unique_proteins_df],
            how="left",
            on=[condition_col, "Protein"]
        )
    )

    df_a = df.iloc[:, :7]

    original_aas = [col.split("(")[1][0] for col in df.columns[7:]]
    
    df_b = df.iloc[:, 7:].apply(
        lambda row: [
            update_cell(
                cell, 
                original_aa, 
                0, 
                np.nan, 
                1
            )
            for cell, original_aa in zip(row, original_aas)
        ],
        axis=1,
        result_type="broadcast"
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
                    check_nan_func=lambda cell: np.isnan(cell)
                )
                for cell, original_aa, imputed_nan_value in zip(row, original_aas, mean_col_editing_freqs_wo_nan)
            ],
            axis=1,
            result_type="broadcast"
        )
        
    df_b = drop_uniformative_aa_cols(df_b)

    df = pd.concat([df_a, df_b], axis=1)

    df.insert(6, "MeanTotalSupportingReads", (df["TotalEqualSupportingReads"] + df["TotalWeightedSupportingReads"]) / 2)
    df.insert(7, "%EqualTotalExpression", 100 * df["TotalEqualSupportingReads"] / df["TotalEqualSupportingReads"].sum())
    df.insert(8, "%WeightedTotalExpression", 100 * df["TotalWeightedSupportingReads"] / df["TotalWeightedSupportingReads"].sum())
    df.insert(9, "%MeanTotalExpression", 100 * df["MeanTotalSupportingReads"] / df["MeanTotalSupportingReads"].sum())
    df = df.sort_values(sorting_col, ascending=False)
    
    return df


# %%
equal_exp_tsne_input_dfs = [
    prepare_ml_input_df(max_sol_df, unique_proteins_df, unique_proteins_first_col_pos, sorting_col="%EqualTotalExpression")
    for max_sol_df, unique_proteins_df in zip(max_sol_dfs, unique_proteins_dfs)
]
# equal_exp_tsne_input_dfs[0]

weighted_exp_tsne_input_dfs = [
    prepare_ml_input_df(max_sol_df, unique_proteins_df, unique_proteins_first_col_pos, sorting_col="%WeightedTotalExpression")
    for max_sol_df, unique_proteins_df in zip(max_sol_dfs, unique_proteins_dfs)
]
weighted_exp_tsne_input_dfs[0]


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


# %% tags=[]
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
    top_expressed_proteins=None
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
                n_jobs=n_jobs
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
    top_expressed_proteins=None
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


# %% [markdown]
# > All proteins

# %%
# perplexities = [5, 30, 50, 100]
perplexities = [5, 30, 50, 100, 150]
n_iter = 500
n_jobs = 40

# %%
n_iter_500_equal_conditions_tsnes, n_iter_500_equal_conditions_Xs = run_tsnes(conditions, equal_exp_tsne_input_dfs, seed, perplexities=perplexities, n_iter=n_iter, n_jobs=n_jobs)
n_iter_500_weighted_conditions_tsnes, n_iter_500_weighted_conditions_Xs = run_tsnes(conditions, weighted_exp_tsne_input_dfs, seed, perplexities=perplexities, n_iter=n_iter, n_jobs=n_jobs)

# %% tags=[]
rank_cutoff = 300

for conditions_tsnes, conditions_Xs, sorting_method in zip(
    [n_iter_500_equal_conditions_tsnes, n_iter_500_weighted_conditions_tsnes],
    [n_iter_500_equal_conditions_Xs, n_iter_500_weighted_conditions_Xs],
    ["equal", "weighted"]
):

    # head_title = (
    #     f"t-SNEs for largest solution of each {str(condition_col).lower()} under different perplexities, sorted by % of total {sorting_method} expression"
    #     "<br>"
    #     f"<sub>{rank_cutoff} highest expressed proteins are colored</sub>"
    # )
    head_title = (
        f"t-SNEs for largest solution of each {str(condition_col).lower()} under different perplexities"
        "<br>"
        f"<sub>{rank_cutoff} highest expressed proteins (according to % of total {sorting_method} expression) are colored</sub>"
    )
    row_titles = conditions
    column_titles = [f"Perplexity = {perplexity}" for perplexity in perplexities]

    fig = make_subplots(
        rows=len(conditions),
        cols=len(perplexities),
        row_titles=row_titles,
        column_titles=column_titles,
        # shared_yaxes=True,
        # shared_xaxes=True
    )

    for row, (condition, X, condition_tsnes) in enumerate(zip(conditions, conditions_Xs, conditions_tsnes), start=1):

        n = X.shape[0]
        color_options = [color_discrete_map[condition], "white"]
        colors = color_highest_expressed_proteins(n, rank_cutoff, color_options)

        for col, prots_perplexity_tsne in enumerate(condition_tsnes, start=1):

            x, y = prots_perplexity_tsne.T

            fig.add_trace(
                go.Scattergl(
                    x=x,
                    y=y,
                    mode="markers",
                    marker=dict(
                        color=colors,
                        line_width=0.5
                    )
                ),
                row=row,
                col=col,
            )

    fig.update_layout(
        title_text=head_title,
        title_y=0.95,
        template=template,
        showlegend=False,
        width=1200,
        height=600,
    )

    fig.show()

# %% tags=[]
rank_cutoff = 300

for conditions_tsnes, conditions_Xs, sorting_method in zip(
    [n_iter_500_equal_conditions_tsnes, n_iter_500_weighted_conditions_tsnes],
    [n_iter_500_equal_conditions_Xs, n_iter_500_weighted_conditions_Xs],
    ["equal", "weighted"]
):

    # head_title = (
    #     f"t-SNEs for largest solution of each {str(condition_col).lower()} under different perplexities, sorted by % of total {sorting_method} expression"
    #     "<br>"
    #     f"<sub>{rank_cutoff} highest expressed proteins are colored</sub>"
    # )
    head_title = (
        f"t-SNEs for largest solution of each {str(condition_col).lower()} under different perplexities"
        "<br>"
        f"<sub>{rank_cutoff} highest expressed proteins (according to % of total {sorting_method} expression) are colored</sub>"
    )
    row_titles = conditions
    column_titles = [f"Perplexity = {perplexity}" for perplexity in perplexities]

    fig = make_subplots(
        rows=len(conditions),
        cols=len(perplexities),
        row_titles=row_titles,
        column_titles=column_titles,
        # shared_yaxes=True,
        # shared_xaxes=True
    )

    for row, (condition, X, condition_tsnes) in enumerate(zip(conditions, conditions_Xs, conditions_tsnes), start=1):

        n = X.shape[0]
        color_options = [color_discrete_map[condition], "white"]
        colors = color_highest_expressed_proteins(n, rank_cutoff, color_options)[:rank_cutoff]

        for col, prots_perplexity_tsne in enumerate(condition_tsnes, start=1):

            x, y = prots_perplexity_tsne.T
            x = x[:rank_cutoff]
            y = y[:rank_cutoff]

            fig.add_trace(
                go.Scattergl(
                    x=x,
                    y=y,
                    mode="markers",
                    marker=dict(
                        color=colors,
                        line_width=0.5
                    )
                ),
                row=row,
                col=col,
            )

    fig.update_layout(
        title_text=head_title,
        title_y=0.95,
        template=template,
        showlegend=False,
        width=1200,
        height=600,
    )

    fig.show()

# %%
equal_conditions_pcas = run_pcas(conditions, equal_exp_tsne_input_dfs, seed)
weighted_conditions_pcas = run_pcas(conditions, weighted_exp_tsne_input_dfs, seed)

# %% tags=[]
rank_cutoff = 300

for conditions_pcas, sorting_method in zip(
    [equal_conditions_pcas, weighted_conditions_pcas],
    ["equal", "weighted"]
):

    head_title = (
        f"PCAs for largest solution of each {str(condition_col).lower()}"
        "<br>"
        f"<sub>{rank_cutoff} highest expressed proteins (according to % of total {sorting_method} expression) are colored</sub>"
    )
    column_titles = conditions

    fig = make_subplots(
        rows=1,
        cols=len(conditions),
        # row_titles=row_titles,
        column_titles=column_titles,
        shared_yaxes=True,
        shared_xaxes=True,
        x_title="PC1",
        y_title="PC2",
    )

    for col, (condition, condition_pca) in enumerate(zip(conditions, conditions_pcas), start=1):
        
        n = len(condition_pca)
        color_options = [color_discrete_map[condition], "white"]
        colors = color_highest_expressed_proteins(n, rank_cutoff, color_options)
        
        x = condition_pca["PC1"]
        y = condition_pca["PC2"]
        
        fig.add_trace(
            go.Scattergl(
                x=x,
                y=y,
                mode="markers",
                marker=dict(
                    color=colors,
                    line_width=0.5
                )
            ),
            row=1,
            col=col,
        )

    fig.update_layout(
        title_text=head_title,
        title_y=0.95,
        template=template,
        showlegend=False,
        width=800,
        height=500,
    )

    fig.show()

# %% [markdown]
# > Top 1000 expressed proteins

# %%
# perplexities = [5, 30, 50, 100]
perplexities = [5, 30, 50, 100, 150, 200]

# %%
top_expressed_proteins = 1000
top_1000_equal_conditions_tsnes, top_1000_equal_conditions_Xs = run_tsnes(conditions, equal_exp_tsne_input_dfs, seed, perplexities=perplexities, top_expressed_proteins=top_expressed_proteins)
top_1000_weighted_conditions_tsnes, top_1000_weighted_conditions_Xs = run_tsnes(conditions, weighted_exp_tsne_input_dfs, seed, perplexities=perplexities, top_expressed_proteins=top_expressed_proteins)

# %% tags=[]
rank_cutoff = 100

for conditions_tsnes, conditions_Xs, sorting_method in zip(
    [top_1000_equal_conditions_tsnes, top_1000_weighted_conditions_tsnes],
    [top_1000_equal_conditions_Xs, top_1000_weighted_conditions_Xs],
    ["equal", "weighted"]
):

    # head_title = (
    #     f"t-SNEs for top {top_expressed_proteins} expressed proteins in largest solution of each {str(condition_col).lower()} under different perplexities, sorted by % of total {sorting_method} expression"
    #     "<br>"
    #     f"<sub>{rank_cutoff} highest expressed proteins are colored</sub>"
    # )
    head_title = (
        f"t-SNEs for top {top_expressed_proteins} expressed proteins in largest solution of each {str(condition_col).lower()} under different perplexities"
        "<br>"
        f"<sub>{rank_cutoff} highest expressed proteins (according to % of total {sorting_method} expression) are colored</sub>"
    )
    
    row_titles = conditions
    column_titles = [f"Perplexity = {perplexity}" for perplexity in perplexities]

    fig = make_subplots(
        rows=len(conditions),
        cols=len(perplexities),
        row_titles=row_titles,
        column_titles=column_titles,
        # shared_yaxes=True,
        # shared_xaxes=True
    )

    for row, (condition, X, condition_tsnes) in enumerate(zip(conditions, conditions_Xs, conditions_tsnes), start=1):

        n = X.shape[0]
        color_options = [color_discrete_map[condition], "white"]
        colors = color_highest_expressed_proteins(n, rank_cutoff, color_options)

        for col, prots_perplexity_tsne in enumerate(condition_tsnes, start=1):

            x, y = prots_perplexity_tsne.T

            fig.add_trace(
                go.Scattergl(
                    x=x,
                    y=y,
                    mode="markers",
                    marker=dict(
                        color=colors,
                        line_width=0.5
                    )
                ),
                row=row,
                col=col,
            )

    fig.update_layout(
        title_text=head_title,
        # title_pad_b=3,
        title_y=0.95,
        template=template,
        showlegend=False,
        width=1200,
        height=600,
    )

    fig.show()

# %%
top_expressed_proteins = 1000
top_1000_equal_conditions_pcas = run_pcas(conditions, equal_exp_tsne_input_dfs, seed, top_expressed_proteins=top_expressed_proteins)
top_1000_weighted_conditions_pcas = run_pcas(conditions, weighted_exp_tsne_input_dfs, seed, top_expressed_proteins=top_expressed_proteins)

# %% tags=[]
rank_cutoff = 100

for conditions_pcas, sorting_method in zip(
    [top_1000_equal_conditions_pcas, top_1000_weighted_conditions_pcas],
    ["equal", "weighted"]
):

    head_title = (
        f"PCAs for top {top_expressed_proteins} expressed proteins in largest solution of each {str(condition_col).lower()}"
        "<br>"
        f"<sub>{rank_cutoff} highest expressed proteins (according to % of total {sorting_method} expression) are colored</sub>"
    )
    column_titles = conditions

    fig = make_subplots(
        rows=1,
        cols=len(conditions),
        # row_titles=row_titles,
        column_titles=column_titles,
        shared_yaxes=True,
        shared_xaxes=True,
        x_title="PC1",
        y_title="PC2",
    )

    for col, (condition, condition_pca) in enumerate(zip(conditions, conditions_pcas), start=1):
        
        n = len(condition_pca)
        color_options = [color_discrete_map[condition], "white"]
        colors = color_highest_expressed_proteins(n, rank_cutoff, color_options)
        
        x = condition_pca["PC1"]
        y = condition_pca["PC2"]
        
        fig.add_trace(
            go.Scattergl(
                x=x,
                y=y,
                mode="markers",
                marker=dict(
                    color=colors,
                    line_width=0.5
                )
            ),
            row=1,
            col=col,
        )

    fig.update_layout(
        title_text=head_title,
        title_y=0.95,
        template=template,
        showlegend=False,
        width=800,
        height=500,
    )

    fig.show()


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
    p = (
        max_sol_exp_df
        .iloc[:, first_col_pos:]
        .apply(np.mean)
    )
    s_hypo = sum(-p * np.log2(p) - (1 - p) * np.log2((1 - p)))
    return s_hypo


# %%
def calc_entropies(max_sol_exp_df, first_col_pos, prcnt_equal_exp_col="%EqualTotalExpression", prcnt_weighted_exp_col="%WeightedTotalExpression"):
    s_data_equal_exp, s_data_weighted_exp = calc_data_entropy(max_sol_exp_df, prcnt_equal_exp_col, prcnt_weighted_exp_col)
    s_hypo = calc_hypothetical_entropy(max_sol_exp_df, first_col_pos)
    return s_data_equal_exp, s_data_weighted_exp, s_hypo


# %%
max_sol_exp_dfs = [
    prepare_ml_input_df(max_sol_df, unique_proteins_df, unique_proteins_first_col_pos, sorting_col="%WeightedTotalExpression", impute_nans=False)
    for max_sol_df, unique_proteins_df in zip(max_sol_dfs, unique_proteins_dfs)
]

max_sol_entropies = [
    calc_entropies(max_sol_exp_df, ML_INPUT_FIRST_COL_POS, "%EqualTotalExpression", "%WeightedTotalExpression")
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
        "EntropyValue": _entropies
    }
)

shannon_df

# %% tags=[]
fig = px.bar(
    shannon_df,
    x="EntropyName",
    y="EntropyValue",
    facet_col=condition_col,
    facet_col_spacing=facet_col_spacing,
    # labels={"EntropyName": "", "EntropyValue": "Entropy"},
    labels={"EntropyName": "Entropy", "EntropyValue": ""},
    title=f"Shannon's entropy of a largest solution of each {condition_col.lower()}",
    color=condition_col,
    color_discrete_map=color_discrete_map,
    pattern_shape="EntropyName",
    pattern_shape_sequence=["+", "x", ""],
    # barmode='group',
    category_orders=category_orders,
    template=template,
)
fig.update_layout(showlegend=False)
fig.show()

# %% [markdown] tags=[] toc-hr-collapsed=true
# ## Editing in reads

# %% [markdown]
# ### Edited sites

# %% [markdown] papermill={"duration": 0.115362, "end_time": "2022-02-01T09:42:53.418861", "exception": false, "start_time": "2022-02-01T09:42:53.303499", "status": "completed"} tags=[]
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


# %% [markdown] papermill={"duration": 0.115362, "end_time": "2022-02-01T09:42:53.418861", "exception": false, "start_time": "2022-02-01T09:42:53.303499", "status": "completed"} tags=[]
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

# %% [markdown] papermill={"duration": 0.115362, "end_time": "2022-02-01T09:42:53.418861", "exception": false, "start_time": "2022-02-01T09:42:53.303499", "status": "completed"} tags=[]
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


# %% [markdown] papermill={"duration": 0.115362, "end_time": "2022-02-01T09:42:53.418861", "exception": false, "start_time": "2022-02-01T09:42:53.303499", "status": "completed"} tags=[]
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

# %% tags=[]
cols = len(conditions)

x_title = "Non-syn mutations per protein"
y_title = "Proteins"
title_text = "Distribution of min & max estimates of non-syn mutations per protein"

fig = make_subplots(
    rows=1,
    # rows=2,
    cols=cols,
    subplot_titles=conditions,
    shared_yaxes=True,
    x_title=x_title,
    y_title=y_title,
)

min_x = None
max_x = 0

for col, condition, proteins_df in zip(range(1, cols + 1), conditions, proteins_dfs):
    col_names = ["MinNonSyns", "MaxNonSyns"]
    estimate_names = ["Min", "Max"]

    for i, (col_name, estimate_name) in enumerate(zip(col_names, estimate_names)):

        x = proteins_df[col_name]

        fig.add_trace(
            go.Histogram(
                x=x,
                marker_color=subcolors_discrete_map[condition][i],
                name=f"{condition}, {estimate_name}",
            ),
            row=1,
            # row=i + 1,
            col=col,
        )

        min_x = min(min_x, x.min()) if min_x else x.min()
        max_x = max(max_x, x.max())

fig.update_layout(
    template=template,
    barmode="overlay",  # Overlay both histograms
    title_text=title_text,
    legend_title_text=f"{condition_col}, Estimate",
)

fig.update_traces(opacity=0.75)  # Reduce opacity to see both histograms
fig.update_xaxes(range=[min_x * 0.9, max_x * 1.1])

fig.show()

