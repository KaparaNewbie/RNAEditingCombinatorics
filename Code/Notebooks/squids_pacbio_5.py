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

# %% [markdown] papermill={"duration": 0.029907, "end_time": "2022-02-01T09:42:43.198426", "exception": false, "start_time": "2022-02-01T09:42:43.168519", "status": "completed"}
# # Imports

# %%
code_dir = "/private7/projects/Combinatorics/Code"

# %%
# %load_ext autoreload
# %autoreload 2

# %% papermill={"duration": 2.901153, "end_time": "2022-02-01T09:42:46.125355", "exception": false, "start_time": "2022-02-01T09:42:43.224202", "status": "completed"}
import sys
from collections import defaultdict
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
)

# %% papermill={"duration": 0.071769, "end_time": "2022-02-01T09:42:43.049672", "exception": false, "start_time": "2022-02-01T09:42:42.977903", "status": "completed"} tags=["parameters"]
condition_col = "Gene"
conditions = ["GRIA", "PCLO"]
chroms = ["comp141693_c0_seq1", "comp141882_c0_seq14"]
starts = [170, 0]
ends = [2999, 6294]
strands = ["+", "+"]
unaligned_bam_files = [
    "/private7/projects/Combinatorics/D.pealeii/Data/CCS/BasicCCS/GRIA-CNS-RESUB.C0x1291.ccs.bam",
    "/private7/projects/Combinatorics/D.pealeii/Data/CCS/BasicCCS/PCLO-CNS-RESUB.C0x1291.ccs.bam",
]
reads_type = "CCS"  # something like CCS / miseq / etc.
aligned_bam_files = [
    "/private7/projects/Combinatorics/D.pealeii/Alignment/BestN1/GRIA-CNS-RESUB.C0x1291.aligned.sorted.bam",
    "/private7/projects/Combinatorics/D.pealeii/Alignment/BestN1/PCLO-CNS-RESUB.C0x1291.aligned.sorted.bam",
]
filtered_aligned_bam_files = [
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.bam",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.bam",
]
include_flags = None
exclude_flags = "2304"  # remove secondary and supplementary (chimeric) alignments
sep = "\t"
positions_files = [
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.positions.csv.gz",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.positions.csv.gz",
]
reads_files = [
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.reads.csv.gz",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.reads.csv.gz",
]
unique_reads_files = [
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_reads.csv.gz",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_reads.csv.gz",
]
proteins_files = [
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.proteins.csv.gz",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.proteins.csv.gz",
]
unique_proteins_files = [
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz",
]
reads_first_col_pos = 6
unique_reads_first_col_pos = 8
proteins_first_col_pos = 12
unique_proteins_first_col_pos = 14
reads_editing_col = "EditingFrequency"
proteins_editing_col = "MinNonSyns"
# distinct_unique_reads_files = [
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA.AllRows.DistinctUniqueReads.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO.AllRows.DistinctUniqueReads.csv",
# ]

# # old distinct files
# distinct_unique_proteins_files = [
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/GRIA-CNS-RESUB.DistinctUniqueProteins.03.03.2023-15:36:38.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/PCLO-CNS-RESUB.DistinctUniqueProteins.03.03.2023-15:53:01.csv",
# ]

# distinct files with available supporting reads
distinct_unique_proteins_files = [
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/GRIA-CNS-RESUB.DistinctUniqueProteins.Regular.AvailableReads.22.05.2023-11:08:27.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/PCLO-CNS-RESUB.DistinctUniqueProteins.Regular.AvailableReads.22.05.2023-11:26:31.csv",
]

distinct_unique_proteins_fraction01_files = [
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/GRIA-CNS-RESUB.DistinctUniqueProteins.Fraction0_1.11.05.2023-17:03:49.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/PCLO-CNS-RESUB.DistinctUniqueProteins.Fraction0_1.11.05.2023-17:11:00.csv",
]


# proteins_jaccard_files = [
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA.AllRows.JaccardMatrixProteins.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO.AllRows.JaccardMatrixProteins.csv",
# ]

# # old expression files
# expression_files = [
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/GRIA.DistinctUniqueProteins.ExpressionLevels.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/PCLO.DistinctUniqueProteins.ExpressionLevels.csv",
# ]
# fraction01_expression_files = [
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/GRIA.DistinctUniqueProteins.ExpressionLevels.Fraction0_1.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/PCLO.DistinctUniqueProteins.ExpressionLevels.Fraction0_1.csv"
# ]

# new expression files with supporting reads and unique proteins
expression_files = [
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/GRIA.DistinctUniqueProteins.ExpressionLevels.Regular.SupportingReadsAndProteins.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/PCLO.DistinctUniqueProteins.ExpressionLevels.Regular.SupportingReadsAndProteins.csv",
]
fraction01_expression_files = [
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/GRIA.DistinctUniqueProteins.ExpressionLevels.Fraction0_1.SupportingReadsAndProteins.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/PCLO.DistinctUniqueProteins.ExpressionLevels.Fraction0_1.SupportingReadsAndProteins.csv",
]

distinct_dissimilar_miyata_proteins_files = [
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/GRIA-CNS-RESUB.DistinctUniqueProteins.AAgroupsMiyata1979.03.03.2023-15:49:55.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/PCLO-CNS-RESUB.DistinctUniqueProteins.AAgroupsMiyata1979.03.03.2023-16:09:24.csv",
]
miyata_expression_files = [
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/GRIA.DistinctUniqueProteins.ExpressionLevels.AAgroupsMiyata1979.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/PCLO.DistinctUniqueProteins.ExpressionLevels.AAgroupsMiyata1979.csv",
]
grantham_cutoff_scores = [
    # 50, 75,
    100,
    # 125, 150
]
distinct_dissimilar_grantham_proteins_files = [
    [
        "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/GRIA-CNS-RESUB.DistinctUniqueProteins.GRANTHAM1974-100.08.03.2023-12:37:40.csv",
    ],
    [
        "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/PCLO-CNS-RESUB.DistinctUniqueProteins.GRANTHAM1974-100.08.03.2023-13:53:05.csv",
    ],
]
grantham_expression_files = [
    [
        "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/GRIA.DistinctUniqueProteins.ExpressionLevels.GRANTHAM1974-100.csv",
    ],
    [
        "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/PCLO.DistinctUniqueProteins.ExpressionLevels.GRANTHAM1974-100.csv",
    ],
]
alg_repetitions = 5
known_sites_file = (
    "/private7/projects/Combinatorics/D.pealeii/Annotations/D.pea.EditingSites.csv"
)
samtools_path = "/home/alu/kobish/anaconda3/envs/combinatorics/bin/samtools"
threads = 20
code_dir = "/private7/projects/Combinatorics/Code"
seed = 1892


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


# %% papermill={"duration": 0.054755, "end_time": "2022-02-01T09:42:46.304499", "exception": false, "start_time": "2022-02-01T09:42:46.249744", "status": "completed"}
# plotly consts
# color_sequence = px.colors.qualitative.Pastel
# color_sequence = px.colors.qualitative.D3
color_sequence = px.colors.qualitative.G10
color_discrete_map = {
    condition: color for condition, color in zip(conditions, color_sequence)
}
subcolors_discrete_map = {
    condition: two_subcolors_from_hex(color_discrete_map[condition])
    for condition in conditions
}
ic(color_discrete_map)
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
template = "plotly_white"
facet_col_wrap = 6
facet_row_spacing = facet_col_spacing * 6
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
n_repetitions_colormap(subcolors_discrete_map, "GRIA", 10)

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

# %%
25_000 / (750 * 1.17)

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


# %%
known_non_syns_df = known_sites_df.assign(
    NonSyn=known_sites_df["OriginalAA"] != known_sites_df["NewAA"]
)
known_non_syns_df

# %%
known_non_syns_per_chrom_df = (
    known_non_syns_df.groupby("Chrom")["NonSyn"]
    .count()
    .reset_index()
    .rename(columns={"NonSyn": "NonSyns"})
    .sort_values("NonSyns", ascending=False)
)
known_non_syns_per_chrom_df

# %%
known_non_syns_per_chrom_df.loc[known_non_syns_per_chrom_df["NonSyns"] >= 102]

# %%
fig = px.histogram(
    known_non_syns_per_chrom_df,
    x="NonSyns",
    log_y=True,
    # cumulative=True,
    template=template,
)
# fig['layout']['xaxis']['autorange'] = "reversed" # reverse the x-axis
fig.show()

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
reads_dfs[0]["EditedPositions"].mean()

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
unique_proteins_dfs[0]

# %%
unique_proteins_first_col_pos

# %%
editable_aas_per_sample = [
    df.iloc[:, unique_proteins_first_col_pos:].shape[1] for df in unique_proteins_dfs
]
editable_aas_per_sample

# %%
avg_editables_aas_per_sample = sum(editable_aas_per_sample) / len(unique_proteins_dfs)

print(f"Average of {avg_editables_aas_per_sample:.0f} editable AAs per sample")

# %%
unique_proteins_dfs[0].iloc[:, unique_proteins_first_col_pos:]

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
len(unique_proteins_dfs[0].columns[unique_proteins_first_col_pos:])

# %%
# unique_edited_proteins_dfs = [
#     unique_proteins_df.loc[unique_proteins_df[proteins_editing_col] > 0]
#     for unique_proteins_df in unique_proteins_dfs
# ]
# unique_edited_proteins_dfs[0]


# %%
unique_proteins_dfs[1].iloc[:, unique_proteins_first_col_pos:]

# %%
possible_nonsyns = [
    "DG",
    "EG",
    "HR",
    "IM",
    "IV",
    "KE",
    "KR",
    "MV",
    "ND",
    "NS",
    "QR",
    "RG",
    "SG",
    "TA",
    "YC",
    "*W",
]

# %%
below_0_blosum62_nonsyns = ["DG", "EG", "RG", "YC", "*W"]  # -1  # -2  # -2  # -2


# %%
def contains_blosum62_nonsyn(col: pd.Series, below_0_blosum62_nonsyns):
    original_aa = col.name.split("(")[1][0]
    for cell in col:
        aas = cell.split(",")
        for aa in aas:
            if (
                f"{original_aa}{aa}" in below_0_blosum62_nonsyns
                or f"{aa}{original_aa}" in below_0_blosum62_nonsyns
            ):
                return True
    return False


# %%
def contains_blosum62_nonsyn_2(col: pd.Series, below_0_blosum62_nonsyns):
    cells = [cell.split(",") for cell in set(col)]
    for x in range(len(cells) - 1):
        aas_x = cells[x]
        for y in range(x + 1, len(cells)):
            aas_y = cells[y]
            for aa_x in aas_x:
                for aa_y in aas_y:
                    if (
                        f"{aa_x}{aa_y}" in below_0_blosum62_nonsyns
                        or f"{aa_y}{aa_x}" in below_0_blosum62_nonsyns
                    ):
                        return True
    return False


# %%
cols_which_contain_blosum62_nonsyn = [
    col
    for col in unique_proteins_dfs[1].iloc[:, unique_proteins_first_col_pos:].columns
    if contains_blosum62_nonsyn(unique_proteins_dfs[1][col], below_0_blosum62_nonsyns)
]
unique_proteins_dfs[1].loc[:, cols_which_contain_blosum62_nonsyn]

# %%
cols_which_contain_blosum62_nonsyn_2 = [
    col
    for col in unique_proteins_dfs[1].iloc[:, unique_proteins_first_col_pos:].columns
    if contains_blosum62_nonsyn_2(unique_proteins_dfs[1][col], below_0_blosum62_nonsyns)
]
unique_proteins_dfs[1].loc[:, cols_which_contain_blosum62_nonsyn_2]

# %%

# %% [markdown]
# ### Distinct unique proteins

# %%
unique_reads_dfs[0]

# %%
distinct_unique_proteins_dfs = []
for condition, distinct_unique_proteins_file, unique_reads_df in zip(
    conditions, distinct_unique_proteins_files, unique_reads_dfs
):
    distinct_unique_proteins_df = pd.read_csv(
        distinct_unique_proteins_file, sep=sep
    ).drop("AvailableReads", axis=1)
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

# %%
distinct_unique_proteins_fraction01_dfs = []
for condition, distinct_unique_proteins_fraction01_file, unique_reads_df in zip(
    conditions, distinct_unique_proteins_fraction01_files, unique_reads_dfs
):
    distinct_unique_proteins_fraction01_df = pd.read_csv(
        distinct_unique_proteins_fraction01_file, sep=sep
    ).drop("AvailableReads", axis=1)
    distinct_unique_proteins_fraction01_df.insert(0, condition_col, condition)
    distinct_unique_proteins_fraction01_df.insert(
        1,
        "NumOfReads",
        (
            distinct_unique_proteins_fraction01_df["Fraction"]
            * unique_reads_df["NumOfReads"].sum()
        ).astype(int),
    )
    distinct_unique_proteins_fraction01_dfs.append(
        distinct_unique_proteins_fraction01_df
    )

distinct_unique_proteins_fraction01_df = (
    pd.concat(distinct_unique_proteins_fraction01_dfs)
    .reset_index(drop=True)
    .rename(columns={"NumUniqueSamples": "NumOfProteins", "UniqueSamples": "Proteins"})
)

distinct_unique_proteins_fraction01_df = (
    distinct_unique_proteins_fraction01_df.sort_values(
        [
            condition_col,
            "Fraction",
            "FractionRepetition",
            "Algorithm",
            "AlgorithmRepetition",
        ]
    ).reset_index(drop=True)
)

distinct_unique_proteins_fraction01_df


# %%
sol12 = distinct_unique_proteins_fraction01_dfs[0].iloc[11]
sol12

# %%

# %% [markdown]
# #### Distinct dissimilar

# %%
distinct_dissimilar_miyata_proteins_dfs = []
for condition, distinct_dissimilar_miyata_file, unique_reads_df in zip(
    conditions, distinct_dissimilar_miyata_proteins_files, unique_reads_dfs
):
    distinct_dissimilar_miyata_proteins_df = pd.read_csv(
        distinct_dissimilar_miyata_file, sep=sep
    )
    distinct_dissimilar_miyata_proteins_df.insert(0, condition_col, condition)
    distinct_dissimilar_miyata_proteins_df.insert(
        1,
        "NumOfReads",
        (
            distinct_dissimilar_miyata_proteins_df["Fraction"]
            * unique_reads_df["NumOfReads"].sum()
        ).astype(int),
    )
    distinct_dissimilar_miyata_proteins_dfs.append(
        distinct_dissimilar_miyata_proteins_df
    )

distinct_dissimilar_miyata_proteins_df = (
    pd.concat(distinct_dissimilar_miyata_proteins_dfs)
    .reset_index(drop=True)
    .rename(columns={"NumUniqueSamples": "NumOfProteins", "UniqueSamples": "Proteins"})
)

distinct_dissimilar_miyata_proteins_df = (
    distinct_dissimilar_miyata_proteins_df.sort_values(
        [
            condition_col,
            "Fraction",
            "FractionRepetition",
            "Algorithm",
            "AlgorithmRepetition",
        ]
    ).reset_index(drop=True)
)

distinct_dissimilar_miyata_proteins_df

# %%
distinct_dissimilar_miyata_proteins_df.groupby(condition_col)["NumOfProteins"].max()

# %%
distinct_dissimilar_grantham_proteins_dfs = []
for (
    condition,
    condition_distinct_dissimilar_grantham_proteins_files,
    unique_reads_df,
) in zip(conditions, distinct_dissimilar_grantham_proteins_files, unique_reads_dfs):
    for cutoff_score, distinct_dissimilar_grantham_file in zip(
        grantham_cutoff_scores, condition_distinct_dissimilar_grantham_proteins_files
    ):
        distinct_dissimilar_grantham_proteins_df = pd.read_csv(
            distinct_dissimilar_grantham_file, sep=sep
        )
        distinct_dissimilar_grantham_proteins_df.insert(0, condition_col, condition)
        distinct_dissimilar_grantham_proteins_df.insert(1, "CutoffScore", cutoff_score)
        distinct_dissimilar_grantham_proteins_df.insert(
            2,
            "NumOfReads",
            (
                distinct_dissimilar_grantham_proteins_df["Fraction"]
                * unique_reads_df["NumOfReads"].sum()
            ).astype(int),
        )
        distinct_dissimilar_grantham_proteins_dfs.append(
            distinct_dissimilar_grantham_proteins_df
        )

distinct_dissimilar_grantham_proteins_df = (
    pd.concat(distinct_dissimilar_grantham_proteins_dfs)
    .reset_index(drop=True)
    .rename(columns={"NumUniqueSamples": "NumOfProteins", "UniqueSamples": "Proteins"})
)

distinct_dissimilar_grantham_proteins_df = (
    distinct_dissimilar_grantham_proteins_df.sort_values(
        [
            condition_col,
            "CutoffScore",
            "Fraction",
            "FractionRepetition",
            "Algorithm",
            "AlgorithmRepetition",
        ]
    ).reset_index(drop=True)
)

distinct_dissimilar_grantham_proteins_df

# %%
distinct_dissimilar_grantham_proteins_df.groupby([condition_col, "CutoffScore"])[
    "NumOfProteins"
].max()

# %% [markdown]
# #### Jaccards

# %%
# proteins_jaccard_dfs = []
# for jaccard_file in proteins_jaccard_files:
#     df = pd.read_csv(jaccard_file, sep=sep)
#     proteins_jaccard_dfs.append(df)
# proteins_jaccard_dfs[0]

# %%
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
unaligned_reads_counts = [
    count_reads_in_unaligned_bam(samtools_path, bam, threads)
    for bam in unaligned_bam_files
]

# aligned_reads_counts = [
#     count_reads(
#         samtools_path,
#         bam,
#         f"{chrom}:{start+1}-{end}",
#         include_flags,
#         exclude_flags,
#         threads,
#     )
#     for bam, chrom, start, end in zip(aligned_bam_files, chroms, starts, ends)
# ]
aligned_reads_counts = [
    count_unique_filtered_aligned_reads(
        bam,
        f"{chrom}:{start+1}-{end}",
        include_flags,
        exclude_flags,
        threads,
    )
    for bam, chrom, start, end in zip(aligned_bam_files, chroms, starts, ends)
]

# filtered_aligned_reads_counts = [
#     count_reads(
#         samtools_path,
#         bam,
#         f"{chrom}:{start+1}-{end}",
#         include_flags,
#         exclude_flags,
#         threads,
#     )
#     for bam, chrom, start, end in zip(filtered_aligned_bam_files, chroms, starts, ends)
# ]
filtered_aligned_reads_counts = [
    count_unique_filtered_aligned_reads(
        bam,
        f"{chrom}:{start+1}-{end}",
        include_flags,
        exclude_flags,
        threads,
    )
    for bam, chrom, start, end in zip(filtered_aligned_bam_files, chroms, starts, ends)
]

# %%
pileup_reads_counts = [
    len(set(chain.from_iterable(positions_df["Reads"].str.split(","))))
    for positions_df in positions_dfs
]


# %%
unique_reads_counts = [len(unique_reads_df) for unique_reads_df in unique_reads_dfs]

# %%
# max_fraction = distinct_unique_reads_df["Fraction"].max()


# %%
# distinct_unique_reads_counts = (
#     distinct_unique_reads_df.loc[distinct_unique_reads_df["Fraction"] == max_fraction]
#     .groupby(condition_col)["NumUniqueReads"]
#     .mean()
#     .round()
#     .astype(int)
# )


# %%
distinct_unique_proteins_counts = (
    distinct_unique_proteins_df.loc[
        # distinct_unique_proteins_df["Fraction"] == max_fraction
        distinct_unique_proteins_df["Fraction"]
        == 1.0
    ]
    .groupby(condition_col)["NumOfProteins"]
    # .mean()
    # .round()
    .max()
    .astype(int)
)


# %%
data_loss_df = pd.DataFrame(
    {
        f"Unaligned {reads_type} reads": unaligned_reads_counts,
        "Aligned reads (within ORF)": aligned_reads_counts,
        "Filtered aligned reads (within ORF)": filtered_aligned_reads_counts,
        "Pileup reads": pileup_reads_counts,
        "Unique reads": unique_reads_counts,
        # "Distinct unique reads (mean)": distinct_unique_reads_counts,
        "Distinct unique proteins (max)": distinct_unique_proteins_counts,
    },
    index=conditions,
)
data_loss_df


# %% [markdown] papermill={"duration": 0.045853, "end_time": "2022-02-01T09:42:48.953594", "exception": false, "start_time": "2022-02-01T09:42:48.907741", "status": "completed"}
# # Results

# %% [markdown] papermill={"duration": 0.149848, "end_time": "2022-02-01T09:43:12.800733", "exception": false, "start_time": "2022-02-01T09:43:12.650885", "status": "completed"}
# ## Data loss

# %%
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

# long_data_loss_df = long_data_loss_df.sort_values("Count")

long_data_loss_df["StageType"] = long_data_loss_df.apply(
    lambda x: "DataProcessing"
    if x["Stage"] in ["Unique reads", "Distinct unique proteins<br>(mean)"]
    else "Results",
    axis=1,
)

long_data_loss_df


# %%
# df = long_data_loss_df.loc[long_data_loss_df["Stage"] != "Distinct unique reads<br>(mean)"].sort_values("Count")

# # category_orders = category_orders | {
# #     "Stage": ['Unaligned CCS reads', 'Aligned reads<br>(within ORF)',
# #        'Pileup reads', 'Filtered aligned reads<br>(within ORF)',
# #        'Unique reads', 'Distinct unique proteins<br>(mean)']
# # }

# # pattern_shape_map = {
# #         "Unaligned CCS reads": "",
# #         'Aligned reads<br>(within ORF)': "",
# #         'Pileup reads': "",
# #         'Filtered aligned reads<br>(within ORF)': "",
# #         'Unique reads': "/",
# #         'Distinct unique proteins<br>(mean)': "/"
# #     }

# pattern_shape_map = {
#     "DataProcessing": "/",
#     "Results": ""
# }

# category_orders = horizontal_category_orders

# # category_orders = horizontal_category_orders | {
# #     # "StageType": ["DataProcessing", "Results"]
# #     "StageType": ["Results", "DataProcessing"]
# # }


# fig = px.bar(
#     df,
#     y="Stage",
#     x="Count",
#     text="%InitialCount",
#     title="Data loss through the different processing stages",
#     labels={"Count": "", "Stage": ""},
#     color=condition_col,
#     color_discrete_map=color_discrete_map,
#     category_orders=category_orders,

#     # pattern_shape="Stage",
#     # pattern_shape_map=pattern_shape_map,

#     pattern_shape="StageType",
#     pattern_shape_map=pattern_shape_map,

#     template=template,
#     barmode="group",
#     orientation="h",
#     width=800,
#     height=550,
# )

# fig.update_layout(
#     # showlegend=False,
#     bargap=0.25, # gap between bars of adjacent location coordinates.
#     bargroupgap=0.1 # gap between bars of the same location coordinate.
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

# fig.update_traces(width=0.4)

# fig.show()


# %%
fig = go.Figure()

# colors = [color_discrete_map[condition] for condition in conditions]

# non_syns_per_max_sol_exp_df = [
#     str(max_sol_exp_df.iloc[:, ML_INPUT_FIRST_COL_POS:].shape[1])
#     for max_sol_exp_df in max_sol_exp_dfs
# ]

for condition in list(reversed(conditions)):
    df = long_data_loss_df.loc[long_data_loss_df[condition_col] == condition]
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

fig.update_traces(width=0.35)
fig.update_yaxes(
    categoryorder="array",
    # categoryarray=['Unaligned CCS reads', 'Aligned reads<br>(within ORF)', 'Pileup reads', 'Filtered aligned reads<br>(within ORF)', 'Unique reads', 'Distinct unique proteins<br>(mean)']
    categoryarray=list(
        reversed(
            [
                "Unaligned CCS reads",
                "Aligned reads<br>(within ORF)",
                "Filtered aligned reads<br>(within ORF)",
                "Pileup reads",
                "Unique reads",
                "Distinct unique proteins<br>(mean)",
            ]
        )
    ),
)

fig.update_layout(
    # title="Data loss through the different processing stages",
    # orientation="h",
    width=800,
    height=600,
    # showlegend=False,
    template=template,
    barmode="group",
    legend_traceorder="reversed",
    bargap=0.3,  # gap between bars of adjacent location coordinates.
    bargroupgap=0,  # gap between bars of the same location coordinate.
)

fig.write_image(
    "Data loss through the different processing stages - PacBio.svg",
    width=800,
    height=600,
)

fig.show()

# %%
list(
    reversed(
        [
            "Unaligned CCS reads",
            "Aligned reads<br>(within ORF)",
            "Pileup reads",
            "Filtered aligned reads<br>(within ORF)",
            "Unique reads",
            "Distinct unique proteins<br>(mean)",
        ]
    )
)

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

    title = f"Pearson correlation coefficient between editing sites in {condition} reads"
    
    plt.title(title)
    
    if condition == "PCLO":
        fig.savefig(f"{title}.svg", dpi=300);


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
symmetric_pclo_mi_df = get_symmetric_mi_df(
    mi_df=mi_dfs[1], positions=reads_w_nan_dfs[1].iloc[:, reads_first_col_pos:].columns
)

symmetric_pclo_mi_df

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

title = f"Normalized mutual information between editing sites in {conditions[1]}"

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

# %%

# %%

# %%

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
fig.update_yaxes(
    title="% noise",
    # tickmode='linear',
    # tick0=0,
    # dtick=2
)
fig.update_layout(
    showlegend=False, 
    width=300, height=350,
    title_text="Squid's PacBio",
    title_x=0.27
)
fig.write_image(
    "Per chrom noise levels - PacBio.svg",
    width=300,
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

[len(x) for x in conditions_sets[conditions[0]]]

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


cols = min(facet_col_wrap, len(conditions), 4)
rows = ceil(len(conditions) / cols)

fig, axs = plt.subplots(
    nrows=rows,
    ncols=cols,
    figsize=(3.5 * cols, 2.5 * rows),
    constrained_layout=True,
    gridspec_kw=dict(hspace=0.2, wspace=0.03),
)

for condition, ax in zip(conditions, axs.flat):
    labels = conditions_labels[condition]
    sets = conditions_sets[condition]
    labels[0] = f"Edited\n({len(sets[0])})"
    labels[1] = f"Known editing\n({len(sets[1])})"
    if len(sets[2]) == 0:
        labels = labels[:2]
        sets = sets[:2]
        v_func = venn2
    else:
        v_func = venn3
        problamatic_regions_exist = True
    v_func(sets, set_labels=labels, ax=ax)
    ax.set_title(condition, fontdict=dict(fontsize=12))

# if problamatic_regions_exist:
#     title = "Positions' membership: currently edited, known editing, and probalamtic regions"
# else:
#     title = "Positions' membership: currently edited & known editing"

fig.suptitle("Squid's PacBio", fontsize="xx-large", 
             # y=1.2
            )

plt.savefig("Known vs new editing sites - PacBio.svg", format="svg", dpi=300)

plt.show()


# %%
ic(38 + 65)
ic(65 + 8);

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

# %% papermill={"duration": 4.052404, "end_time": "2022-02-01T09:42:53.176715", "exception": false, "start_time": "2022-02-01T09:42:49.124311", "status": "completed"}
# todo retain nan rows and turn nans to 0?

cols = len(conditions)

fig = make_subplots(
    rows=1,
    cols=cols,
    subplot_titles=conditions,
    shared_yaxes=True,
    shared_xaxes=True,
    x_title="% editing",
    y_title="% known editing",
)

df = merged_ref_base_positions_df.copy()
df = df.loc[df["Edited"] & df["KnownEditing"]]

for col, condition, unique_reads_df in zip(
    range(1, cols + 1), conditions, unique_reads_dfs
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
        row=1,
        col=col,
    )

    fig.add_annotation(
        row=1,
        col=col,
        x=20,
        y=75,
        xref="x",
        yref="y",
        text=f"<b>Pearson's r</b><br>p-val = {pv:.2e}<br>ρ = {r:.2g}",
        bgcolor="white",
        borderpad=4,
        opacity=0.8,
        showarrow=False,
    )

fig.update_layout(
    title_text="Correlation between current & previously-reported editing levels",
    showlegend=False,
    template=template,
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

# https://stackoverflow.com/questions/58167028/single-axis-caption-in-plotly-express-facet-plot
for axis in fig.layout:
    if type(fig.layout[axis]) == go.layout.YAxis:
        fig.layout[axis].title.text = ""
fig.update_layout(showlegend=False, yaxis_title="Positions")
# fig.for_each_yaxis(lambda y: y.update(title="Transcripts"))   # for facet_row

fig.show()


# %% [markdown] papermill={"duration": 0.030615, "end_time": "2022-02-01T09:42:49.024262", "exception": false, "start_time": "2022-02-01T09:42:48.993647", "status": "completed"}
# ## Num of distinct unique proteins

# %%
fraction_1_gdf = distinct_unique_proteins_df.loc[
    (distinct_unique_proteins_df["Fraction"] == 1.0)
    # & (distinct_unique_proteins_df["Algorithm"] == "Descending")
].groupby(condition_col)

# means = (
#     fraction_1_gdf["NumOfProteins"]
#     .mean()
#     .reset_index()
#     .rename(columns={"NumOfProteins": "Mean"})
# )

maxes = (
    fraction_1_gdf["NumOfProteins"]
    .max()
    .reset_index()
    .rename(columns={"NumOfProteins": "Max"})
)

stds = (
    fraction_1_gdf["NumOfProteins"]
    .std()
    .reset_index()
    .rename(columns={"NumOfProteins": "STD"})
)

# mean_distinct_proteins_df = means.merge(stds, on=condition_col)

# mean_distinct_proteins_df

max_distinct_proteins_df = maxes.merge(stds, on=condition_col)

max_distinct_proteins_df

# %%
fig = px.bar(
    max_distinct_proteins_df,
    x=condition_col,
    y="Max",
    # error_y="STD",
    # points="all",
    labels={"Max": "Distinct proteins"},
    color=condition_col,
    color_discrete_map=color_discrete_map,
    category_orders=category_orders,
    template=template,
)

fig.update_yaxes(range=[0, max_distinct_proteins_df["Max"].max()])

# width=max(70 * len(conditions), 300)
width = max(50 * len(conditions), 300)
height = 400

fig.update_layout(
    showlegend=False, 
    xaxis_title="Transcript", 
    width=width, 
    height=height,
    title_text="Squid's PacBio",
    # title_x=0.3
    title_x=0.27
)

fig.write_image(
    "Max distinct proteins per transcript - PacBio.svg",
    width=width,
    height=height,
)

fig.show()


# %% jupyter={"source_hidden": true}
# fig = px.bar(
#     mean_distinct_proteins_df,
#     x=condition_col,
#     y="Mean",
#     error_y="STD",
#     # points="all",
#     labels={"Mean": "Distinct proteins (avg)"},
#     color=condition_col,
#     color_discrete_map=color_discrete_map,
#     category_orders=category_orders,
#     template=template,
# )

# fig.update_yaxes(range=[0, mean_distinct_proteins_df["Mean"].max()])

# fig.update_layout(
#     showlegend=False, xaxis_title="", width=max(70 * len(conditions), 300), height=400
# )

# fig.write_image(
#     "Distinct proteins (avg) - PacBio.svg",
#     width=max(70 * len(conditions), 300),
#     height=400,
# )

# fig.show()


# %% jupyter={"source_hidden": true}
# distinct_proteins_per_editable_aas_df = distinct_unique_proteins_df.copy()
# # distinct_proteins_per_editable_aas_df = distinct_proteins_per_editable_aas_df.loc[distinct_proteins_per_editable_aas_df["Algorithm"] == "Descending"]
# distinct_proteins_per_editable_aas_df = distinct_proteins_per_editable_aas_df.loc[distinct_proteins_per_editable_aas_df["Fraction"] == 1.0]
# # distinct_proteins_per_mapped_reads_df["NumOfProteins/MappedReads"] = distinct_proteins_per_mapped_reads_df["NumOfProteins"] / distinct_proteins_per_mapped_reads_df["NumOfReads"]
# editable_aas_per_sample_dict = {condition: editable_aas for condition, editable_aas in zip(conditions, editable_aas_per_sample)}
# distinct_proteins_per_editable_aas_df["EditableAAs"] = [editable_aas_per_sample_dict[condition] for condition in distinct_proteins_per_editable_aas_df[condition_col]]
# # df["NumOfProteins/Read/EditableAmicoAcids"] = df["NumOfProteins"] / df["NumOfReads"] / df["EditableAminoAcidsPerSample"]
# distinct_proteins_per_editable_aas_df


# %% jupyter={"source_hidden": true}
# mean_distinct_proteins_per_editable_aas_df = pd.DataFrame(
#     {
#         condition_col: conditions,
#         "EditableAAS": editable_aas_per_sample,
#     }
# )

# fraction_1_gdf = distinct_proteins_per_editable_aas_df.loc[distinct_proteins_per_editable_aas_df["Fraction"] == 1.0].groupby(condition_col)

# means = fraction_1_gdf["NumOfProteins"].mean().reset_index().rename(columns={"NumOfProteins": "Mean"})
# stds = fraction_1_gdf["NumOfProteins"].std().reset_index().rename(columns={"NumOfProteins": "SD"})

# mean_distinct_proteins_per_editable_aas_df = mean_distinct_proteins_per_editable_aas_df.merge(means).merge(stds)

# mean_distinct_proteins_per_editable_aas_df

# %% jupyter={"source_hidden": true}
# fig = go.Figure()

# for condition in conditions:
#     condition_df = distinct_proteins_per_editable_aas_df.loc[distinct_proteins_per_editable_aas_df[condition_col] == condition]

#     x = condition_df["EditableAAs"]
#     y = condition_df["NumOfProteins"]

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
# x = distinct_proteins_per_editable_aas_df["EditableAAs"]
# y = distinct_proteins_per_editable_aas_df["NumOfProteins"]
# r, pv = scipy.stats.pearsonr(x, y)

# fig.add_annotation(
#     x=30,
#     y=25_000,
#     xref="x",
#     yref="y",
#     text=f"<b>Pearson's r</b><br>p-val = {pv:.2e}<br>ρ = {r:.2g}",
#     bgcolor="white",
#     borderpad=4,
#     font=dict(size=12),
#     opacity=0.8,
#     showarrow=False,
# )

# fig.update_xaxes(range=[0, distinct_proteins_per_editable_aas_df["EditableAAs"].max() * 1.1])
# fig.update_yaxes(range=[0, distinct_proteins_per_editable_aas_df["NumOfProteins"].max() * 1.1])

# fig.update_layout(
#     # showlegend=False,
#     template=template,
#     xaxis_title="Editable amino acids",
#     yaxis_title="Distinct unique proteins",
#     width=600,
#     height=400
# )

# fig.show()


# %% jupyter={"source_hidden": true}
# distinct_proteins_per_mapped_reads_df = distinct_unique_proteins_df.copy()
# distinct_proteins_per_mapped_reads_df = distinct_proteins_per_mapped_reads_df.loc[distinct_proteins_per_mapped_reads_df["Algorithm"] == "Descending"]
# distinct_proteins_per_mapped_reads_df = distinct_proteins_per_mapped_reads_df.loc[distinct_proteins_per_mapped_reads_df["Fraction"] == 1.0]
# distinct_proteins_per_mapped_reads_df["NumOfProteins/MappedReads"] = distinct_proteins_per_mapped_reads_df["NumOfProteins"] / distinct_proteins_per_mapped_reads_df["NumOfReads"]
# editable_aas_per_sample_dict = {condition: editable_aas for condition, editable_aas in zip(conditions, editable_aas_per_sample)}
# distinct_proteins_per_mapped_reads_df["EditableAAs"] = [editable_aas_per_sample_dict[condition] for condition in distinct_proteins_per_mapped_reads_df[condition_col]]
# df["NumOfProteins/Read/EditableAmicoAcids"] = df["NumOfProteins"] / df["NumOfReads"] / df["EditableAminoAcidsPerSample"]
# distinct_proteins_per_mapped_reads_df


# %% jupyter={"source_hidden": true}
# mean_distinct_proteins_per_mapped_reads_df = pd.DataFrame(
#     {
#         condition_col: conditions,
#         "EditableAAS": editable_aas_per_sample,
#     }
# )

# fraction_1_gdf = distinct_proteins_per_mapped_reads_df.loc[distinct_proteins_per_mapped_reads_df["Fraction"] == 1.0].groupby(condition_col)

# means = fraction_1_gdf["NumOfProteins/MappedReads"].mean().reset_index().rename(columns={"NumOfProteins/MappedReads": "Mean"})
# stds = fraction_1_gdf["NumOfProteins/MappedReads"].std().reset_index().rename(columns={"NumOfProteins/MappedReads": "SD"})

# mean_distinct_proteins_per_mapped_reads_df = mean_distinct_proteins_per_mapped_reads_df.merge(means).merge(stds)

# mean_distinct_proteins_per_mapped_reads_df

# %% jupyter={"source_hidden": true}
# fig = go.Figure()

# for condition in conditions:
#     condition_df = distinct_proteins_per_mapped_reads_df.loc[distinct_proteins_per_mapped_reads_df[condition_col] == condition]

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
#     x=30,
#     y=0.3,
#     xref="x",
#     yref="y",
#     text=f"<b>Pearson's r</b><br>p-val = {pv:.2e}<br>ρ = {r:.2g}",
#     bgcolor="white",
#     borderpad=4,
#     font=dict(size=12),
#     opacity=0.8,
#     showarrow=False,
# )

# fig.update_xaxes(range=[0, distinct_proteins_per_mapped_reads_df["EditableAAs"].max() * 1.1])
# fig.update_yaxes(range=[0, distinct_proteins_per_mapped_reads_df["NumOfProteins/MappedReads"].max() * 1.1])

# fig.update_layout(
#     # showlegend=False,
#     template=template,
#     xaxis_title="Editable amino acids",
#     yaxis_title="Distinct unique proteins /<br>mapped reads",
#     width=600,
#     height=400
# )

# fig.show()


# %% jupyter={"source_hidden": true}
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
#     showlegend=False,
#     xaxis_title="",
#     width=max(70*len(conditions), 300),
#     height=400
# )

# fig.show()


# %% jupyter={"source_hidden": true}
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

# fig.update_yaxes(range=[0, distinct_proteins_per_mapped_reads_df["NumOfProteins/MappedReads"].max()])

# fig.update_layout(
#     showlegend=False,
#     xaxis_title="",
#     width=max(70*len(conditions), 300),
#     height=400
# )

# fig.show()


# %% jupyter={"source_hidden": true}

# %%
editable_aas_per_sample

# %%
distinct_unique_proteins_df

# %%
distinct_unique_proteins_df["NumOfProteins"].max()


# %% jupyter={"source_hidden": true}
# condition = conditions[0]
# df = distinct_unique_proteins_df.loc[
#     distinct_unique_proteins_df[condition_col] == condition
# ]
# grouped_df = df.groupby("Fraction")
# x_fraction_mean = grouped_df["NumOfReads"].mean().reset_index()
# x_fraction_mean


# %% jupyter={"source_hidden": true}
# y_fraction_mean = grouped_df["NumOfProteins"].mean().reset_index()
# y_fraction_mean


# %% jupyter={"source_hidden": true}
# x_fraction_mean.merge(y_fraction_mean, on="Fraction")


# %% jupyter={"source_hidden": true}
# grantham_cutoff_scores

# %% jupyter={"source_hidden": true}
# distinct_dissimilar_miyata_proteins_df.groupby(condition_col)["NumOfProteins"].max()

# %% jupyter={"source_hidden": true}
# distinct_dissimilar_grantham_proteins_df.loc[distinct_dissimilar_grantham_proteins_df["CutoffScore"] == 100].groupby(condition_col)["NumOfProteins"].max()

# %% jupyter={"source_hidden": true}
# x_axis_name = "Reads"
# y_axis_name = "Distinct unique proteins"
# head_title = (
#     "Distinct unique proteins vs. sequencing depth"
#     # "<br>"
#     # # f"<sub>({alg_repetitions * 2} repetitions over each fraction of data)</sub>"
#     # "<sub>(100 repetitions over each fraction of data)</sub>"
# )
# _marker_size = 5
# maximal_x = 0

# # Initialize figure with subplots
# fig = make_subplots(
#     rows=1, cols=1, print_grid=False, x_title=x_axis_name, y_title=y_axis_name
# )

# _distinct_dfs = {
#     "Regular": distinct_unique_proteins_df,
#     "Miyata": distinct_dissimilar_miyata_proteins_df,
#     # "Grantham": distinct_dissimilar_grantham_proteins_df
# }
# for cutoff_score in grantham_cutoff_scores:
#     _distinct_dfs[f"Grantham {cutoff_score}"] = distinct_dissimilar_grantham_proteins_df.loc[distinct_dissimilar_grantham_proteins_df["CutoffScore"] == cutoff_score]

# _dissimilar_colormaps = {
#     condition: n_repetitions_colormap(subcolors_discrete_map, condition, 3)
#     for condition in conditions
# }
# _distinctions = {"Regular": 1, "Miyata": 2, "Grantham": 3}

# # basic_symbols = ["x", "star", "circle", "square", "diamond", "star-square", "triangle-down", ]
# # symbol_modifiers = [ "", "", "-dot"]

# # symbols = ["circle", "square-dot", "star", "circle", "diamond", "star-square", "triangle-down"]
# symbols = ["circle", "square-dot", "diamond", "circle", "star", "star-square", "triangle-down"]
# n_similarities = 2
# fill_color_similarities = [True] * n_similarities + [False] * (len(symbols) - n_similarities)
# # fill_color_similarities = [False, True, False, False, False, False, False]

# first_data_trace = True

# # Add traces
# for col, condition in enumerate(conditions, start=1):

#     for _distinction, _distinct_df in _distinct_dfs.items():

#         if _distinction.startswith("Grantham"):
#             i = _distinctions[_distinction.split(" ")[0]]   # same as _distinctions["Grantham"]
#             cutoff_score = int(_distinction.split(" ")[1])
#             j = i - 1 + grantham_cutoff_scores.index(cutoff_score)
#         else:
#             i =_distinctions[_distinction]
#             j = i - 1

#         line_color = _dissimilar_colormaps[condition][i]
#         fill_color = line_color if fill_color_similarities[j] else "white"

#         symbol = symbols[j]

#         name = f"{condition} - {_distinction}"

#         df = _distinct_df.loc[
#             (_distinct_df[condition_col] == condition) &
#             (_distinct_df["Algorithm"] == "Descending")
#         ]

#         x_measured = df["NumOfReads"]
#         y_measured = df["NumOfProteins"]

#         if first_data_trace:
#             fig.add_trace(
#                 go.Scatter(
#                     x=x_measured,
#                     y=y_measured,
#                     mode="markers",
#                     # marker=dict(color=subcolors_discrete_map[condition][0], size=_marker_size),
#                     marker=dict(color=fill_color, size=_marker_size, symbol=symbol, line=dict(
#                         width=2,
#                         color=line_color
#                     )),
#                     legendgroup="Full-CDS, PacBio",  # this can be any string
#                     legendgrouptitle_text="Full-CDS, PacBio",
#                     name=name
#                 ),
#             )
#             first_data_trace = False
#         else:
#             fig.add_trace(
#                 go.Scatter(
#                     x=x_measured,
#                     y=y_measured,
#                     mode="markers",
#                     # marker=dict(color=subcolors_discrete_map[condition][0], size=_marker_size),
#                     marker=dict(color=fill_color, size=_marker_size, symbol=symbol, line=dict(
#                         width=2,
#                         color=line_color
#                     )),
#                     legendgroup="Full-CDS, PacBio",  # this can be any string
#                     name=name
#                 ),
#             )

#         grouped_df = df.groupby("Fraction")
#         x_fraction_mean = grouped_df["NumOfReads"].mean().reset_index()
#         y_fraction_mean = grouped_df["NumOfProteins"].mean().reset_index()
#         mean_fraction_df = x_fraction_mean.merge(y_fraction_mean, on="Fraction")

#         maximal_x = max(maximal_x, x_measured.max())

#         fig.add_trace(
#             go.Scatter(
#                 x=mean_fraction_df["NumOfReads"],
#                 y=mean_fraction_df["NumOfProteins"],
#                 mode="lines",
#                 # marker=dict(color=subcolors_discrete_map[condition][0], size=_marker_size),
#                 line=dict(
#                     color=line_color,
#                     width=_marker_size*.2,
#                 ),
#                 opacity=0.5,
#                 showlegend=False,
#             ),
#         )

# dscam_ys = [19_008, 18_496, ]
# dscam_legend_names = ["Theoretical maximum", "Measured", ]
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

# # fig.update_yaxes(type="log")

# fig.update_layout(
#     title_text=head_title,
#     template=template,
#     # legend_font=dict(size=8),
#     # legend_grouptitlefont=dict(size=8),
#     # legend_tracegroupgap=4,
#     # width=100*maximal_x/10
#     height=600,
#     width=1000
# )
# fig.write_image("Distinct unique proteins vs. sequencing depth - PacBio.svg", width=650, height=500)
# fig.show()


# %% jupyter={"source_hidden": true}
# distinct_unique_proteins_df["NumOfProteins"].max()

# %% jupyter={"source_hidden": true}
# np.log10(distinct_unique_proteins_df["NumOfProteins"].max())

# %% jupyter={"source_hidden": true}
# from math import floor

# %% jupyter={"source_hidden": true}
# floor(np.log10(distinct_unique_proteins_df["NumOfProteins"].max()))

# %% jupyter={"source_hidden": true}
# floor(np.log10(x_measured.max()))

# %% jupyter={"source_hidden": true}
# x_measured.max()

# %% jupyter={"source_hidden": true}
# x_axis_name = "Mapped reads"
# y_axis_name = "Distinct proteins"
# head_title = (
#     "Distinct proteins vs. sequencing depth"
#     # "<br>"
#     # # f"<sub>({alg_repetitions * 2} repetitions over each fraction of data)</sub>"
#     # "<sub>(100 repetitions over each fraction of data)</sub>"
# )
# _marker_size = 5
# maximal_x = 0

# # Initialize figure with subplots
# fig = make_subplots(
#     rows=1, cols=1, print_grid=False, x_title=x_axis_name, y_title=y_axis_name
# )

# max_y = 0
# first_data_trace = True

# # Add traces
# for condition in conditions:

#     df = distinct_unique_proteins_df.loc[
#         (distinct_unique_proteins_df[condition_col] == condition)
#         # & (distinct_unique_proteins_df["Algorithm"] == "Descending")
#     ]
#     color = color_discrete_map[condition]
#     name = condition

#     x_measured = df["NumOfReads"]
#     y_measured = df["NumOfProteins"]

#     max_y = max(max_y, y_measured.max())

#     if first_data_trace:
#         fig.add_trace(
#             go.Scatter(
#                 x=x_measured,
#                 y=y_measured,
#                 mode="markers",
#                 # mode="lines+markers",
#                 # marker=dict(color=subcolors_discrete_map[condition][0], size=_marker_size),
#                 marker=dict(
#                     color=color,
#                     size=_marker_size,
#                     # symbol=symbol,
#                     #             line=dict(
#                     #     width=2,
#                     #     # color=line_color
#                     # )
#                 ),
#                 legendgroup="Full-CDS, PacBio",  # this can be any string
#                 legendgrouptitle_text="Full-CDS, PacBio",
#                 name=name,
#             ),
#         )
#         first_data_trace = False
#     else:
#         fig.add_trace(
#             go.Scatter(
#                 x=x_measured,
#                 y=y_measured,
#                 mode="markers",
#                 # mode="lines+markers",
#                 # marker=dict(color=subcolors_discrete_map[condition][0], size=_marker_size),
#                 marker=dict(
#                     color=color,
#                     size=_marker_size,
#                     # symbol=symbol,
#                     #             line=dict(
#                     #     width=2,
#                     #     # color=line_color
#                     # )
#                 ),
#                 legendgroup="Full-CDS, PacBio",  # this can be any string
#                 name=name,
#             ),
#         )

#     grouped_df = df.groupby("Fraction")
#     x_fraction_mean = grouped_df["NumOfReads"].mean().reset_index()
#     y_fraction_mean = grouped_df["NumOfProteins"].mean().reset_index()
#     mean_fraction_df = x_fraction_mean.merge(y_fraction_mean, on="Fraction")

#     maximal_x = max(maximal_x, x_measured.max())

#     fig.add_trace(
#         go.Scatter(
#             x=mean_fraction_df["NumOfReads"],
#             y=mean_fraction_df["NumOfProteins"],
#             mode="lines",
#             # marker=dict(color=subcolors_discrete_map[condition][0], size=_marker_size),
#             line=dict(
#                 color=color,
#                 width=_marker_size * 0.2,
#             ),
#             opacity=0.5,
#             showlegend=False,
#         ),
#     )

# dscam_dashed_lined_width = 3.5
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
# # dscam_colors = ["black", "grey"]
# # dscam_dashes = ["dash", "dot"]
# dscam_dashes = ["dash", "dash"]
# fig.add_trace(
#     go.Scatter(
#         x=[0.05 * maximal_x, 1.05 * maximal_x],
#         y=[dscam_ys[0], dscam_ys[0]],
#         mode="lines",
#         line=dict(
#             color=dscam_colors[0], dash=dscam_dashes[0], width=dscam_dashed_lined_width
#         ),
#         opacity=0.6,
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
#             color=dscam_colors[1], dash=dscam_dashes[1], width=dscam_dashed_lined_width
#         ),
#         opacity=0.6,
#         legendgroup="DSCAM",  # this can be any string
#         name=dscam_legend_names[1],
#         # name=f"DSCAM {dscam_legend_names[0]}",
#     ),
# )

# fig.update_yaxes(
#     # type="log",
#     # # tick0=0,
#     # dtick="D2",
#     # # exponentformat="power",
#     # showexponent='all',
#     # range=[0, (floor(np.log10(max_y)) + ceil(np.log10(max_y))) / 2],
#     range=[0, max(max_y * 1.2, 1.05 * max(dscam_ys))],
#     # zeroline=True
# )
# fig.update_xaxes(range=[0, maximal_x * 1.1])
# fig.update_layout(
#     # title_text=head_title,
#     template=template,
#     # legend_font=dict(size=8),
#     # legend_grouptitlefont=dict(size=10),
#     # legend_font=dict(size=12),
#     # legend_grouptitlefont=dict(size=14),
#     # legend_font=dict(size=8),
#     # legend_grouptitlefont=dict(size=8),
#     # legend_tracegroupgap=4,
#     # width=100*maximal_x/10
#     height=600,
#     width=1000,
# )
# fig.write_image(
#     "Distinct proteins vs. sequencing depth - PacBio.svg", width=1000, height=600
# )
# fig.show()


# %% jupyter={"source_hidden": true}
# x_axis_name = "Mapped reads"
# y_axis_name = "Distinct proteins"
# head_title = (
#     "Distinct proteins vs. sequencing depth"
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

# max_y = 0
# first_data_trace = True

# # Add traces
# for condition in conditions:

#     df = distinct_unique_proteins_df.loc[distinct_unique_proteins_df[condition_col] == condition]
#     df = df.sort_values(["Fraction", "NumOfProteins"], ascending=False).drop_duplicates("Fraction", ignore_index=True)
    
#     color = color_discrete_map[condition]
#     name = condition

#     x_measured = df["NumOfReads"]
#     y_measured = df["NumOfProteins"]

#     max_y = max(max_y, y_measured.max())

#     if first_data_trace:
#         fig.add_trace(
#             go.Scatter(
#                 x=x_measured,
#                 y=y_measured,
#                 mode="lines+markers",
#                 marker=dict(
#                     color=color,
#                     size=_marker_size,
#                 ),
#                 line=dict(
#                     color=color,
#                     width=_marker_size * 0.2,
#                 ),
#                 legendgroup="Full-CDS, PacBio",  # this can be any string
#                 legendgrouptitle_text="Full-CDS, PacBio",
#                 name=name,
#             ),
#         )
#         first_data_trace = False
#     else:
#         fig.add_trace(
#             go.Scatter(
#                 x=x_measured,
#                 y=y_measured,
#                 mode="lines+markers",
#                 marker=dict(
#                     color=color,
#                     size=_marker_size,
#                 ),
#                 line=dict(
#                     color=color,
#                     width=_marker_size * 0.2,
#                 ),
#                 legendgroup="Full-CDS, PacBio",  # this can be any string
#                 name=name,
#             ),
#         )

#     maximal_x = max(maximal_x, x_measured.max())

# dscam_dashed_lined_width = 3.5
# dscam_ys = [
#     36_016,
#     18_496,
# ]
# dscam_legend_names = [
#     "Theoretical maximum",
#     "Measured",
# ]
# dscam_colors = ["grey", "black"]
# dscam_dashes = ["dash", "dash"]
# fig.add_trace(
#     go.Scatter(
#         x=[0.05 * maximal_x, 1.05 * maximal_x],
#         y=[dscam_ys[0], dscam_ys[0]],
#         mode="lines",
#         line=dict(
#             color=dscam_colors[0], dash=dscam_dashes[0], width=dscam_dashed_lined_width
#         ),
#         opacity=0.6,
#         legendgroup="DSCAM",  # this can be any string
#         legendgrouptitle_text="DSCAM",
#         name=dscam_legend_names[0],
#     ),
# )
# fig.add_trace(
#     go.Scatter(
#         x=[0.05 * maximal_x, 1.05 * maximal_x],
#         y=[dscam_ys[1], dscam_ys[1]],
#         mode="lines",
#         line=dict(
#             color=dscam_colors[1], dash=dscam_dashes[1], width=dscam_dashed_lined_width
#         ),
#         opacity=0.6,
#         legendgroup="DSCAM",  # this can be any string
#         name=dscam_legend_names[1],
#         # name=f"DSCAM {dscam_legend_names[0]}",
#     ),
# )

# fig.update_yaxes(
#     # type="log",
#     # # tick0=0,
#     # dtick="D2",
#     # # exponentformat="power",
#     # showexponent='all',
#     # range=[0, (floor(np.log10(max_y)) + ceil(np.log10(max_y))) / 2],
#     range=[0, max(max_y * 1.2, 1.05 * max(dscam_ys))],
#     # zeroline=True
# )
# fig.update_xaxes(range=[0, maximal_x * 1.1])
# fig.update_layout(
#     # title_text=head_title,
#     template=template,
#     # legend_font=dict(size=8),
#     # legend_grouptitlefont=dict(size=10),
#     # legend_font=dict(size=12),
#     # legend_grouptitlefont=dict(size=14),
#     # legend_font=dict(size=8),
#     # legend_grouptitlefont=dict(size=8),
#     # legend_tracegroupgap=4,
#     # width=100*maximal_x/10
#     height=600,
#     width=1000,
# )
# fig.write_image(
#     "Distinct proteins vs. sequencing depth - PacBio.svg", width=1000, height=600
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

max_y = 0
first_data_trace = True

# Add traces
for condition in conditions:

    df = distinct_unique_proteins_df.loc[distinct_unique_proteins_df[condition_col] == condition]
    df = df.sort_values(["Fraction", "NumOfProteins"], ascending=False).drop_duplicates("Fraction", ignore_index=True)
    
    color = color_discrete_map[condition]
    name = condition

    x_measured = df["NumOfReads"]
    y_measured = df["NumOfProteins"]

    max_y = max(max_y, y_measured.max())

    if first_data_trace:
        fig.add_trace(
            go.Scatter(
                x=x_measured,
                y=y_measured,
                mode="lines+markers",
                marker=dict(
                    color=color,
                    size=_marker_size,
                ),
                line=dict(
                    color=color,
                    width=_marker_size * 0.2,
                ),
                # legendgroup="Full-CDS, PacBio",  # this can be any string
                # legendgrouptitle_text="Full-CDS, PacBio",
                # legendgroup="Full-CDS, squid's PacBio",  # this can be any string
                # legendgrouptitle_text="Full-CDS, squid's PacBio",
                legendgroup="Squid's PacBio, full-CDS",  # this can be any string
                legendgrouptitle_text="Squid's PacBio, full-CDS",
                name=name,
            ),
        )
        first_data_trace = False
    else:
        fig.add_trace(
            go.Scatter(
                x=x_measured,
                y=y_measured,
                mode="lines+markers",
                marker=dict(
                    color=color,
                    size=_marker_size,
                ),
                line=dict(
                    color=color,
                    width=_marker_size * 0.2,
                ),
                # legendgroup="Full-CDS, PacBio",  # this can be any string
                # legendgroup="Full-CDS, squid's PacBio",  # this can be any string
                legendgroup="Squid's PacBio, full-CDS",  # this can be any string
                name=name,
            ),
        )

    maximal_x = max(maximal_x, x_measured.max())

dscam_dashed_lined_width = 3.5
dscam_ys = [
    36_016,
    18_496,
]
dscam_legend_names = [
    "Theoretical maximum",
    "Measured",
]
dscam_colors = ["grey", "black"]
dscam_dashes = ["dash", "dash"]
fig.add_trace(
    go.Scatter(
        x=[0.05 * maximal_x, 1.05 * maximal_x],
        y=[dscam_ys[0], dscam_ys[0]],
        mode="lines",
        line=dict(
            color=dscam_colors[0], dash=dscam_dashes[0], width=dscam_dashed_lined_width
        ),
        opacity=0.6,
        legendgroup="Drosophila’s Dscam",  # this can be any string
        legendgrouptitle_text="Drosophila’s Dscam",
        name=dscam_legend_names[0],
    ),
)
fig.add_trace(
    go.Scatter(
        x=[0.05 * maximal_x, 1.05 * maximal_x],
        y=[dscam_ys[1], dscam_ys[1]],
        mode="lines",
        line=dict(
            color=dscam_colors[1], dash=dscam_dashes[1], width=dscam_dashed_lined_width
        ),
        opacity=0.6,
        legendgroup="Drosophila’s Dscam",  # this can be any string
        name=dscam_legend_names[1],
        # name=f"DSCAM {dscam_legend_names[0]}",
    ),
)

fig.update_yaxes(
    # type="log",
    # # tick0=0,
    # dtick="D2",
    # # exponentformat="power",
    # showexponent='all',
    # range=[0, (floor(np.log10(max_y)) + ceil(np.log10(max_y))) / 2],
    range=[0, max(max_y * 1.2, 1.05 * max(dscam_ys))],
    # zeroline=True
)
fig.update_xaxes(range=[0, maximal_x * 1.1])
fig.update_layout(
    # title_text=head_title,
    template=template,
    legend_font=dict(size=10),
    legend_grouptitlefont=dict(size=12),
    
    # legend_font=dict(size=8),
    # legend_grouptitlefont=dict(size=10),
    # legend_font=dict(size=12),
    # legend_grouptitlefont=dict(size=14),
    # legend_font=dict(size=8),
    # legend_grouptitlefont=dict(size=8),
    # legend_tracegroupgap=4,
    # width=100*maximal_x/10
    height=600,
    width=900,
)
fig.write_image(
    "Distinct proteins vs. sequencing depth - PacBio.svg", width=900, height=600
)
fig.show()


# %% jupyter={"source_hidden": true}
# fig = make_subplots(
#     rows=1,
#     cols=1,
#     print_grid=False,
# )
# _distinct_dissimilar_dfs = {
#     # "Regular": distinct_unique_proteins_df,
#     "Miyata": distinct_dissimilar_miyata_proteins_df,
#     "Grantham 100": distinct_dissimilar_grantham_proteins_df.loc[
#         distinct_dissimilar_grantham_proteins_df["CutoffScore"] == 100
#     ],
# }

# _distinctions = {"Regular": 1, "Miyata": 2, "Grantham": 3}
# _dissimilar_colormaps = {
#     condition: n_repetitions_colormap(subcolors_discrete_map, condition, 3)
#     for condition in conditions
# }

# symbols = [
#     "circle",
#     "square-dot",
#     "diamond",
#     "circle",
#     "star",
#     "star-square",
#     "triangle-down",
# ]
# n_similarities = 2
# fill_color_similarities = [True] * n_similarities + [False] * (
#     len(symbols) - n_similarities
# )
# _marker_size = 3
# maximal_x = 0

# # Add traces
# for condition in conditions:
#     for _distinction, _distinct_df in _distinct_dissimilar_dfs.items():
#         if _distinction.startswith("Grantham"):
#             i = _distinctions[
#                 _distinction.split(" ")[0]
#             ]  # same as _distinctions["Grantham"]
#             cutoff_score = int(_distinction.split(" ")[1])
#             j = i - 1 + grantham_cutoff_scores.index(cutoff_score)
#         else:
#             i = _distinctions[_distinction]
#             j = i - 1
#         line_color = _dissimilar_colormaps[condition][i]
#         fill_color = line_color if fill_color_similarities[j] else "white"
#         symbol = symbols[j]
#         name = f"{condition} - {_distinction}"
#         df = _distinct_df.loc[
#             (_distinct_df[condition_col] == condition)
#             # & (_distinct_df["Algorithm"] == "Descending")
#         ]
#         x_measured = df["NumOfReads"]
#         y_measured = df["NumOfProteins"]
#         fig.add_trace(
#             go.Scatter(
#                 x=x_measured,
#                 y=y_measured,
#                 mode="markers",
#                 # marker=dict(color=subcolors_discrete_map[condition][0], size=_marker_size),
#                 marker=dict(
#                     color=fill_color,
#                     size=_marker_size,
#                     symbol=symbol,
#                     line=dict(width=2, color=line_color),
#                 ),
#                 name=name,
#             )
#         )
#         grouped_df = df.groupby("Fraction")
#         x_fraction_mean = grouped_df["NumOfReads"].mean().reset_index()
#         y_fraction_mean = grouped_df["NumOfProteins"].mean().reset_index()
#         mean_fraction_df = x_fraction_mean.merge(y_fraction_mean, on="Fraction")
#         maximal_x = max(maximal_x, x_measured.max())
#         fig.add_trace(
#             go.Scatter(
#                 x=mean_fraction_df["NumOfReads"],
#                 y=mean_fraction_df["NumOfProteins"],
#                 mode="lines",
#                 # marker=dict(color=subcolors_discrete_map[condition][0], size=_marker_size),
#                 line=dict(
#                     color=line_color,
#                     width=_marker_size * 0.2,
#                 ),
#                 opacity=0.5,
#                 showlegend=False,
#             ),
#         )

# # dscam_ys = [19_008, 18_496, ]
# # dscam_legend_names = ["Theoretical maximum", "Measured", ]
# # # dscam_legend_names = ["measured", "theoretical maximum"]
# # dscam_colors = ["grey", "black"]
# # fig.add_trace(
# #     go.Scatter(
# #         x=[0.05 * maximal_x, 1.05 * maximal_x],
# #         y=[dscam_ys[0], dscam_ys[0]],
# #         mode="lines",
# #         line=dict(
# #             color=dscam_colors[0],
# #             dash="dash",
# #             # width=3
# #         ),
# #         showlegend=False,
# #         # legendgroup="DSCAM",  # this can be any string
# #         # legendgrouptitle_text="DSCAM",
# #         # name=dscam_legend_names[0],
# #         # # name=f"DSCAM {dscam_legend_names[1]}",
# #     ),
# # )
# # fig.add_trace(
# #     go.Scatter(
# #         x=[0.05 * maximal_x, 1.05 * maximal_x],
# #         y=[dscam_ys[1], dscam_ys[1]],
# #         mode="lines",
# #         line=dict(
# #             color=dscam_colors[1],
# #             dash="dash",
# #             # width=3
# #         ),
# #         showlegend=False,
# #         # legendgroup="DSCAM",  # this can be any string
# #         # name=dscam_legend_names[1],
# #         # name=f"DSCAM {dscam_legend_names[0]}",
# #     ),
# # )

# fig.update_yaxes(
#     title="Distinct dissimilar proteins",
#     type="log",
#     # exponentformat="power"
# )
# fig.update_xaxes(range=[0, maximal_x * 1.1], title="Mapped reads")

# fig.update_layout(
#     # title_text=head_title,
#     template=template,
#     # legend=dict(
#     #     orientation="h",
#     #     yanchor="bottom",
#     #     # y=1.02,
#     #     # xanchor="right",
#     #     # x=1
#     # ),
#     # legend_font=dict(size=8),
#     # legend_grouptitlefont=dict(size=8),
#     # legend_tracegroupgap=4,
#     # width=100*maximal_x/10
#     height=500,
#     width=700,
# )
# fig.write_image(
#     "Distinct dissimilar proteins vs. sequencing depth - PacBio.svg",
#     width=700,
#     height=500,
# )
# fig.show()

# %%
fig = make_subplots(
    rows=1,
    cols=1,
    print_grid=False,
)
_distinct_dissimilar_dfs = {
    # "Regular": distinct_unique_proteins_df,
    "Miyata": distinct_dissimilar_miyata_proteins_df,
    "Grantham 100": distinct_dissimilar_grantham_proteins_df.loc[
        distinct_dissimilar_grantham_proteins_df["CutoffScore"] == 100
    ],
}

_distinctions = {"Regular": 1, "Miyata": 2, "Grantham": 3}
_dissimilar_colormaps = {
    condition: n_repetitions_colormap(subcolors_discrete_map, condition, 3)
    for condition in conditions
}

symbols = [
    "circle",
    "square-dot",
    "diamond",
    "circle",
    "star",
    "star-square",
    "triangle-down",
]
n_similarities = 2
fill_color_similarities = [True] * n_similarities + [False] * (
    len(symbols) - n_similarities
)
_marker_size = 3
maximal_x = 0
maximal_y = 0

# Add traces
for condition in conditions:
    for _distinction, _distinct_df in _distinct_dissimilar_dfs.items():
        if _distinction.startswith("Grantham"):
            i = _distinctions[
                _distinction.split(" ")[0]
            ]  # same as _distinctions["Grantham"]
            cutoff_score = int(_distinction.split(" ")[1])
            j = i - 1 + grantham_cutoff_scores.index(cutoff_score)
        else:
            i = _distinctions[_distinction]
            j = i - 1
        line_color = _dissimilar_colormaps[condition][i]
        fill_color = line_color if fill_color_similarities[j] else "white"
        symbol = symbols[j]
        
        # name = f"{condition} - {_distinction}"
        legendgroup = condition
        name = _distinction
        
        df = _distinct_df.loc[
            (_distinct_df[condition_col] == condition)
            # & (_distinct_df["Algorithm"] == "Descending")
        ]
        df = df.sort_values(["Fraction", "NumOfProteins"], ascending=False).drop_duplicates("Fraction", ignore_index=True)
        x_measured = df["NumOfReads"]
        y_measured = df["NumOfProteins"]
     
        fig.add_trace(
            go.Scatter(
                x=x_measured,
                y=y_measured,
                mode="markers",
                marker=dict(
                    color=fill_color,
                    size=_marker_size,
                    symbol=symbol,
                    line=dict(width=2, color=line_color),  
                ),
                name=name,
                legendgroup=legendgroup,  # this can be any string
                legendgrouptitle_text=legendgroup,
            )
        )
        fig.add_trace(
            go.Scatter(
                x=x_measured,
                y=y_measured,
                mode="lines",
                 line=dict(
                    color=line_color,
                    width=_marker_size * 0.2,
                ),
                # name=name,
                showlegend=False
            )
        )
        
        maximal_x = max(maximal_x, x_measured.max())
        maximal_y = max(maximal_y, y_measured.max())

fig.update_xaxes(
    range=[0, maximal_x * 1.1], 
    title="Mapped reads"
)
fig.update_yaxes(
    title="Distinct dissimilar proteins",
    type="log",
    nticks=5,
    range=[0, 0.3+np.log(maximal_y)/np.log(10)]
)

height, width = 450, 500

fig.update_layout(
    template=template,
    title="Squid's PacBio",
    title_x=0.15,
    autosize=False,
    # margin_r=200,
    legend=dict(
        orientation="h",
        entrywidth=70,
        yanchor="bottom",
        # y=1.02,
        y=0.15,
        xanchor="right",
        x=1
        
        # # orientation="h",
        # yanchor="bottom",
        # y=0.1,
        # xanchor="left",
        # # x=0.95,
        # x=0.1,
        # # itemsizing='constant'
    ),
    legend_grouptitlefont_size=8,
    legend_font=dict(size=8),
    legend_tracegroupgap=0,
    # legend_itemwidth=140,
    # legend_entrywidth=0,
    # legend_bordercolor="grey",
    # legend_borderwidth=1,
    height=height, width=width
)
fig.write_image(
    "Distinct dissimilar proteins vs. sequencing depth - PacBio.svg",
    height=height, width=width
)
fig.show()

# %%

# %% jupyter={"source_hidden": true}
# x_axis_name = "Reads"
# y_axis_name = "Distinct unique proteins"
# head_title = (
#     "Distinct unique proteins vs. sequencing depth"
#     # "<br>"
#     # # f"<sub>({alg_repetitions * 2} repetitions over each fraction of data)</sub>"
#     # "<sub>(100 repetitions over each fraction of data)</sub>"
# )
# _marker_size = 5
# maximal_x = 0

# # Initialize figure with subplots
# fig = make_subplots(
#     rows=1, cols=1, print_grid=False, x_title=x_axis_name, y_title=y_axis_name
# )

# _distinct_dfs = {
#     # "Regular": distinct_unique_proteins_df,
#     "Miyata": distinct_dissimilar_miyata_proteins_df,
#     # "Grantham": distinct_dissimilar_grantham_proteins_df
# }
# for cutoff_score in grantham_cutoff_scores:
#     _distinct_dfs[f"Grantham {cutoff_score}"] = distinct_dissimilar_grantham_proteins_df.loc[distinct_dissimilar_grantham_proteins_df["CutoffScore"] == cutoff_score]

# _dissimilar_colormaps = {
#     condition: n_repetitions_colormap(subcolors_discrete_map, condition, 3)
#     for condition in conditions
# }
# _distinctions = {
#     # "Regular": 1,
#     "Miyata": 2, "Grantham": 3
# }

# # basic_symbols = ["x", "star", "circle", "square", "diamond", "star-square", "triangle-down", ]
# # symbol_modifiers = [ "", "", "-dot"]

# # symbols = ["circle", "square-dot", "star", "circle", "diamond", "star-square", "triangle-down"]
# symbols = ["square-dot", "diamond",]
# n_similarities = 2
# fill_color_similarities = [True] * n_similarities + [False] * (len(symbols) - n_similarities)
# # fill_color_similarities = [False, True, False, False, False, False, False]

# first_data_trace = True

# # Add traces
# for col, condition in enumerate(conditions, start=1):

#     for _distinction, _distinct_df in _distinct_dfs.items():

#         if _distinction.startswith("Grantham"):
#             i = _distinctions[_distinction.split(" ")[0]]   # same as _distinctions["Grantham"]
#             cutoff_score = int(_distinction.split(" ")[1])
#             j = i - 1 + grantham_cutoff_scores.index(cutoff_score)
#         else:
#             i =_distinctions[_distinction]
#             j = i - 1

#         line_color = _dissimilar_colormaps[condition][i]
#         fill_color = line_color if fill_color_similarities[j] else "white"

#         symbol = symbols[j]

#         name = f"{condition} - {_distinction}"

#         df = _distinct_df.loc[
#             (_distinct_df[condition_col] == condition) &
#             (_distinct_df["Algorithm"] == "Descending")
#         ]

#         x_measured = df["NumOfReads"]
#         y_measured = df["NumOfProteins"]

#         if first_data_trace:
#             fig.add_trace(
#                 go.Scatter(
#                     x=x_measured,
#                     y=y_measured,
#                     mode="markers",
#                     # marker=dict(color=subcolors_discrete_map[condition][0], size=_marker_size),
#                     marker=dict(color=fill_color, size=_marker_size, symbol=symbol, line=dict(
#                         width=2,
#                         color=line_color
#                     )),
#                     legendgroup="Full-CDS, PacBio",  # this can be any string
#                     legendgrouptitle_text="Full-CDS, PacBio",
#                     name=name
#                 ),
#             )
#             first_data_trace = False
#         else:
#             fig.add_trace(
#                 go.Scatter(
#                     x=x_measured,
#                     y=y_measured,
#                     mode="markers",
#                     # marker=dict(color=subcolors_discrete_map[condition][0], size=_marker_size),
#                     marker=dict(color=fill_color, size=_marker_size, symbol=symbol, line=dict(
#                         width=2,
#                         color=line_color
#                     )),
#                     legendgroup="Full-CDS, PacBio",  # this can be any string
#                     name=name
#                 ),
#             )

#         grouped_df = df.groupby("Fraction")
#         x_fraction_mean = grouped_df["NumOfReads"].mean().reset_index()
#         y_fraction_mean = grouped_df["NumOfProteins"].mean().reset_index()
#         mean_fraction_df = x_fraction_mean.merge(y_fraction_mean, on="Fraction")

#         maximal_x = max(maximal_x, x_measured.max())

#         fig.add_trace(
#             go.Scatter(
#                 x=mean_fraction_df["NumOfReads"],
#                 y=mean_fraction_df["NumOfProteins"],
#                 mode="lines",
#                 # marker=dict(color=subcolors_discrete_map[condition][0], size=_marker_size),
#                 line=dict(
#                     color=line_color,
#                     width=_marker_size*.2,
#                 ),
#                 opacity=0.5,
#                 showlegend=False,
#             ),
#         )

# dscam_ys = [19_008, 18_496, ]
# dscam_legend_names = ["Theoretical maximum", "Measured", ]
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

# # fig.update_yaxes(type="log")

# fig.update_layout(
#     title_text=head_title,
#     template=template,
#     # legend_font=dict(size=8),
#     # legend_grouptitlefont=dict(size=8),
#     # legend_tracegroupgap=4,
#     # width=100*maximal_x/10
#     height=600,
#     width=1000
# )
# fig.write_image("Distinct unique proteins vs. sequencing depth - PacBio.svg", width=650, height=500)
# fig.show()


# %% jupyter={"source_hidden": true, "outputs_hidden": true}
# np.log10(max_y * 2)

# %% jupyter={"source_hidden": true}
# floor(np.log10(max_y))

# %% jupyter={"source_hidden": true}
# ceil(np.log10(max_y))

# %% jupyter={"source_hidden": true}
# (floor(np.log10(max_y)) + ceil(np.log10(max_y))) / 2

# %% [markdown]
# ### NaNs distribution

# %%
ambigous_positions_df = (
    # expanded_distinct_unique_proteins_df.loc[expanded_distinct_unique_proteins_df["Algorithm"] == "Descending"]
    expanded_distinct_unique_proteins_df.groupby(
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
    ambigous_positions_df,
    x="MeanOfAmbigousPositions",
    facet_col=condition_col,
    facet_col_spacing=facet_col_spacing,
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
y_title = "Num solutions a protein<br>is included at (avg)"
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

# %% [markdown] toc-hr-collapsed=true
# #### Coupled diff comparisons

# %%
distinct_unique_proteins_df

# %%
max_distinct_proteins_per_transcript_and_alg_df = distinct_unique_proteins_df.loc[
    distinct_unique_proteins_df["Fraction"] == 1.0
].copy()

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
asc_df = max_distinct_proteins_per_transcript_and_alg_df.loc[
    max_distinct_proteins_per_transcript_and_alg_df["Algorithm"] == "Ascending"
].reset_index(drop=True)
desc_df = max_distinct_proteins_per_transcript_and_alg_df.loc[
    max_distinct_proteins_per_transcript_and_alg_df["Algorithm"] != "Ascending"
].reset_index(drop=True)

ic(asc_df["NumOfProteins"].eq(desc_df["NumOfProteins"]).sum())  # ==
ic(asc_df["NumOfProteins"].gt(desc_df["NumOfProteins"]).sum())  # >
ic(asc_df["NumOfProteins"].lt(desc_df["NumOfProteins"]).sum())
# <

# %%
greater_asc_transcripts = asc_df.loc[
    asc_df["NumOfProteins"].gt(desc_df["NumOfProteins"]), condition_col
]
greater_asc_transcripts

# %%
distinct_unique_proteins_df.loc[
    (distinct_unique_proteins_df["Fraction"] == 1.0)
    & (distinct_unique_proteins_df[condition_col].isin(greater_asc_transcripts))
].groupby([condition_col, "Algorithm"])["NumOfProteins"].value_counts()

# %%
max_distinct_proteins_per_transcript_and_alg_df.loc[
    max_distinct_proteins_per_transcript_and_alg_df[condition_col].isin(
        greater_asc_transcripts
    )
]

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
ic(asc_df["NumOfProteins"].lt(desc_df["NumOfProteins"]).sum())
# <

# %%
greater_asc_transcripts = asc_df.loc[
    asc_df["NumOfProteins"].gt(desc_df["NumOfProteins"]), condition_col
].unique()

ic(len(greater_asc_transcripts))

greater_asc_transcripts

# %%
# df = max_distinct_proteins_per_transcript_and_alg_and_fracrepetition_df.loc[
#     max_distinct_proteins_per_transcript_and_alg_and_fracrepetition_df[
#         condition_col
#     ].isin(greater_asc_transcripts)
# ]
# df = df.drop(
#     ["AlgorithmRepetition", "IsMaxNumOfProteins", "MaxNumOfProteins", "Proteins"],
#     axis=1,
# ).pivot(
#     index=[condition_col, "NumOfReads", "Fraction", "FractionRepetition"],
#     columns="Algorithm",
# )
# df = df.set_axis(df.columns.get_level_values(1).values, axis=1)
# df = df.reset_index()
# df["Desc - Asc"] = df["Descending"] - df["Ascending"]
# df = df.loc[df["Desc - Asc"] < 0].reset_index(drop=True)
# df

# %%
# fig = px.histogram(
#     df,
#     x="Desc - Asc",
#     template=template,
#     # text_auto=True,
#     facet_col="Fraction",
#     # color="Fraction",
#     category_orders={"Fraction": [0.2, 0.4, 0.6, 0.8, 1.0]},
#     # barmode="group"
# )

# # Reduce opacity to see both histograms
# fig.update_traces(opacity=0.75)

# fig.update_layout(
#     # showlegend=False,
#     width=1300,
#     # height=300
#     # width=800,
#     height=400,
#     # barmode='overlay' # Overlay both histograms
# )


# fig.show()

# %%

# %%

# %%
asc_df = distinct_unique_proteins_df.loc[
    distinct_unique_proteins_df["Algorithm"] == "Ascending"
].reset_index(drop=True)
desc_df = distinct_unique_proteins_df.loc[
    distinct_unique_proteins_df["Algorithm"] == "Descending"
].reset_index(drop=True)
desc_df

# %%
desc_df.groupby([condition_col, "Fraction", "FractionRepetition"])[
    "NumOfProteins"
].mean().reset_index()

# %%
desc_df["NumOfProteins"]

# %%
algs_diff_df = desc_df.loc[
    :, [condition_col, "Fraction", "FractionRepetition", "AlgorithmRepetition"]
]
algs_diff_df["Desc - Asc"] = desc_df["NumOfProteins"] - asc_df["NumOfProteins"]


algs_diff_df["Diff"] = algs_diff_df["Desc - Asc"].apply(
    lambda x: "> 0" if x > 0 else "= 0" if x == 0 else "< 0"
)

# algs_diff_df = algs_diff_df.loc[algs_diff_df["Fraction"] == 1.0].reset_index(drop=True)

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
    shared_xaxes=True,
)

algorithms = ["Ascending", "Descending"]
diffs = ["< 0", "= 0", "> 0"]
symbols = ["cross", "diamond", "circle"]

for col, (diff, symbol) in enumerate(zip(diffs, symbols), start=1):
    for row, condition in enumerate(conditions, start=1):
        df = algs_diff_df.loc[
            (algs_diff_df["Diff"] == diff) & (algs_diff_df[condition_col] == condition)
        ]
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
                    marker={"symbol": symbol},
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
                    marker={"symbol": symbol},
                ),
                row=row,
                col=col,
            )

fig.update_layout(
    title_text=head_title, title_y=0.95, template=template, showlegend=False
)

fig.update_xaxes(tick0=0, dtick=0.1, tickangle=-60, matches="x")

min_y = algs_diff_df["Desc - Asc"].min()
max_y = algs_diff_df["Desc - Asc"].max()
abs_y_limit = max(abs(min_y), abs(max_y))
# min_y = min_y * 1.1 if abs(abs_y_limit // min_y // 10) <= 1 else min_y - abs(abs_y_limit // min_y // 10)
min_y = (
    min_y * 1.1
    if abs(abs_y_limit // min_y // 10) <= 1
    else min_y - abs(abs_y_limit // min_y)
)
max_y = (
    max_y * 1.1
    if abs(abs_y_limit // max_y // 10) <= 1
    else max_y + abs(abs_y_limit // max_y)
)

fig.update_yaxes(zerolinewidth=zerolinewidth, range=[min_y, max_y])

fig.show()


# %% [markdown]
# #### Solutions' sizes

# %%
min_max_fraction_1_distinct_prots_df = (
    distinct_unique_proteins_df.loc[distinct_unique_proteins_df["Fraction"] == 1.0]
    .groupby(condition_col)["NumOfProteins"]
    .agg(["min", "max"])
    .reset_index()
)
min_max_fraction_1_distinct_prots_df[
    "%SolutionsDispersion"
] = min_max_fraction_1_distinct_prots_df.apply(
    lambda x: 100 * (x["max"] - x["min"]) / x["max"], axis=1
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
fig.update_yaxes(
    # title="% dispersion<br><sub>100 * (max - min) / max</sub>",
    # title="% of distinct proteins<br>dispersion", 
    title="% dispersion of distinct<br>proteins sets' sizes", 
    range=[0, min_max_fraction_1_distinct_prots_df["%SolutionsDispersion"].max()],
)

fig.update_layout(
    showlegend=False,
    title="Squid's PacBio",
    title_x=0.27,
    # xaxis_title="",
    width=max(70 * len(conditions), 300),
    height=400,
)

fig.write_image(
    "%SolutionsDispersion - PacBio.svg",
    width=max(70 * len(conditions), 300),
    height=400,
)


fig.show()


# %%
# fig = go.Figure(
#     data=[go.Histogram(x=min_max_fraction_1_distinct_prots_df["%SolutionsDispersion"])]
# )

# fig.update_xaxes(title="% dispersion<br><sub>100 * (max - min) / max</sub>")
# fig.update_yaxes(title="Transcripts")

# fig.update_layout(
#     # showlegend=False,
#     width=600,
#     height=400,
#     template=template,
# )

# # fig.write_image("%SolutionsDispersion - PacBio.svg", width=600, height=400)

# fig.show()

# %%
# y_axis_name = "Distinct unique proteins"
# head_title = "Distinct unique proteins vs. heuristic method"

# fig = make_subplots(
#     rows=1,
#     cols=len(conditions),
#     print_grid=False,
#     y_title=y_axis_name,
#     subplot_titles=conditions,
#     # shared_yaxes=True,
# )

# algorithms = ["Ascending", "Descending"]

# for col, condition in enumerate(conditions, start=1):

#     df = distinct_unique_proteins_df.loc[
#         (distinct_unique_proteins_df[condition_col] == condition) &
#         (distinct_unique_proteins_df["Fraction"] == 1.0)
#     ]

#     xs = [
#         df.loc[df["Algorithm"] == algorithm, condition_col] for algorithm in algorithms
#     ]
#     ys = [
#         df.loc[df["Algorithm"] == algorithm, "NumOfProteins"]
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
# # fig.update_yaxes(range=[0, distinct_unique_proteins_df["NumOfProteins"].max()*1.05])
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
# df = annotated_proteins_jaccard_dfs[0]
# # df = df.drop("Gene", axis=1).set_index(["Fraction", "FractionRepetition", "Algorithm", "AlgorithmRepetition"])
# # df = df.set_axis(df.index, axis=1)
# df

# %%
# fractions = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
# total_len = len(df)
# fraction_len = total_len / len(fractions)
# middles = [int((fraction_len / 2) + (fraction_len * x)) for x in range(len(fractions))]
# middles

# %%
# ticks = []
# for x, fraction in enumerate(
#     [fraction for fraction in fractions for _ in range(int(fraction_len))], start=1
# ):
#     if x > total_len:
#         break
#     tick = str(fraction) if x in middles else ""
#     ticks.append(tick)

# %%
# for x, tick in enumerate(ticks, start=1):
#     if tick != "":
#         ic(x, tick)

# %%
# data = df.iloc[:, 5:].to_numpy()

# %%
# data.shape

# %%
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
# import seaborn as sns

# sns.set_theme()
# uniform_data = df.iloc[:, 5:].to_numpy()
# ax = sns.heatmap(uniform_data)

# %% [markdown]
# ##### Distribution

# %%
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

# %% [markdown] papermill={"duration": 0.030615, "end_time": "2022-02-01T09:42:49.024262", "exception": false, "start_time": "2022-02-01T09:42:48.993647", "status": "completed"}
# ## Supporting reads' coverage

# %% [markdown] papermill={"duration": 0.030615, "end_time": "2022-02-01T09:42:49.024262", "exception": false, "start_time": "2022-02-01T09:42:48.993647", "status": "completed"}
# ### Unique reads

# %%
unique_reads_dfs[1]

# %% jupyter={"source_hidden": true}
# missing_values_reads = unique_reads_dfs[1].iloc[:7]["Reads"].str.split(",").explode().reset_index(drop=True)
# # missing_values_reads = unique_reads_dfs[1]["Reads"].str.split(",").explode().reset_index(drop=True)
# # missing_values_reads

# old_to_new_reads_file = "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.OldToNewReads.csv.gz"
# old_to_new_reads_df = pd.read_csv(old_to_new_reads_file, sep="\t")
# old_to_new_reads_df = old_to_new_reads_df.set_index("NewRead")
# # old_to_new_reads_df

# missing_values_reads = [
#     old_to_new_reads_df.at[read, "OldRead"]
#     for read in missing_values_reads
# ]

# missing_values_dir = Path("/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/TestMissingValuesInReads")
# missing_values_dir.mkdir(exist_ok=True)

# missing_values_reads_file = Path(missing_values_dir, "read.txt")
# bam = Path("/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.bam")
# missing_values_reads_bam = Path(missing_values_dir, bam.name)

# with missing_values_reads_file.open("w") as mvrf:
#     for read in missing_values_reads:
#         mvrf.write(read + "\n")

# import subprocess

# filter_cmd = (
#         "samtools "
#         "view "
#         f"-N {missing_values_reads_file} "
#         "--with-header "
#         "--bam "
#         f"--output {missing_values_reads_bam} "
#         f"--threads 10 "
#         f"{bam} "
# )
# subprocess.run(filter_cmd, shell=True)

# index_cmd = f"samtools index -@ 10 {missing_values_reads_bam} "
# subprocess.run(index_cmd, shell=True)

# %%
top_100_unique_reads_combinatorics_df = (
    unique_reads_dfs[1]
    # .sort_values("%WeightedTotalExpression", ascending=False)
    # .reset_index(drop=True)
    .iloc[:100, 8:]
    # .fillna(-1)
    # # .T
)
top_100_unique_reads_combinatorics_df

# %% jupyter={"source_hidden": true}
# top_100_unique_reads_combinatorics_df.iloc[:7, 42:]

# %% jupyter={"source_hidden": true}
# from collections import Counter
# Counter(top_100_unique_reads_combinatorics_df.iloc[:7].values.reshape(-1))

# %% jupyter={"source_hidden": true}
# top_100_unique_reads_combinatorics_df.columns[
#     int(len(top_100_unique_reads_combinatorics_df))
# ]

# %%
data = top_100_unique_reads_combinatorics_df.replace(0, 0.5).replace(-1, 0).values

fig = go.Figure(
    go.Heatmap(
        z=data,
        y=top_100_unique_reads_combinatorics_df.index + 1,
        colorscale=[
            [0, "rgb(192,192,192)"],
            [0.3, "rgb(192,192,192)"],
            [0.375, "black"],
            [0.625, "black"],
            [0.7, "rgb(74,246,38)"],
            [1.0, "rgb(74,246,38)"],
        ],
        colorbar=dict(
            title="Edited?",
            # tick0=0,
            # dtick=1,
            tickmode="array",
            tickvals=[0.15, 0.5, 0.85],
            # tickvals=[1, 2, 3],
            ticktext=["Missing", "No", "Yes"],
            len=0.3,
        ),
    )
)

fig.update_xaxes(
    title="Editing sites",
    side="top",
    # ticks="",
    # ticks="inside",
    # tickcolor="white",
    showticklabels=False,
    # rangeselector_y=1
    # automargin=True,
    title_standoff=0,
)

# fig.add_annotation(
#     xref="x domain",
#     yref="y domain",
#     showarrow=False,
#     # The arrow head will be 25% along the x axis, starting from the left
#     x=top_100_unique_reads_combinatorics_df.columns[int(len(top_100_unique_reads_combinatorics_df))],
#     # The arrow head will be 40% along the y axis, starting from the bottom
#     y=0.4,
#     text="Editing sites",
#     # arrowhead=2,
# )

fig.update_yaxes(
    title="Top 100 expressed unique reads<br>in squid's PCLO (PacBio)",
    autorange="reversed",
    tick0=1,
    range=[1, 100],
)
fig.update_layout(height=700, width=650, template=template, font_size=16)

fig.write_image(
    "Combinatorics of top 100 expressed unique reads in PCLO - PacBio.svg",
    width=650,
    height=700,
)

fig.show()

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
df = unique_reads_dfs[0].assign(
    CumNumOfReads=unique_reads_dfs[0]["NumOfReads"].cumsum()
)
df["%CumNumOfReads"] = 100 * df["CumNumOfReads"] / df["NumOfReads"].sum()
df["UniqueReadWithDecreasingSupport"] = df.index
df

# %%
cummulative_supporting_reads_dfs = []
for unique_reads_df in unique_reads_dfs:
    df = (
        unique_reads_df.loc[:, ["UniqueRead", "NumOfReads"]]
        .sort_values("NumOfReads", ascending=False)
        .reset_index(drop=True)
    )
    df = df.assign(CumNumOfReads=df["NumOfReads"].cumsum())
    df["%CumNumOfReads"] = 100 * df["CumNumOfReads"] / df["NumOfReads"].sum()
    df["NumUniqueReads"] = df.index + 1
    cummulative_supporting_reads_dfs.append(df)
cummulative_supporting_reads_dfs[0]

# %%
df = cummulative_supporting_reads_dfs[0]
df

# %%
df.iloc[:100]

# %%
ax = sns.barplot(
    # data=df.iloc[:10000],
    data=df,
    x="NumUniqueReads",
    y="%CumNumOfReads",
    color="blue",
    # dodge=False
)
ax.set(xscale="log")

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

# %%

# %%
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

# %% [markdown] toc-hr-collapsed=true
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
    # expression_df["Diff5+"] = (
    #     abs(
    #         expression_df["TotalEqualSupportingReads"]
    #         - expression_df["TotalWeightedSupportingReads"]
    #     )
    #     >= 0.05
    #     * (
    #         expression_df["TotalEqualSupportingReads"]
    #         + expression_df["TotalWeightedSupportingReads"]
    #     )
    #     / 2
    # )

    expression_df["AdditionalSupportingReadsIDs"] = expression_df[
        "AdditionalSupportingReadsIDs"
    ].apply(lambda x: "" if pd.isna(x) else [y.split(",") for y in x.split(";")])
    expression_df["AdditionalSupportingProteinsIDs"] = expression_df[
        "AdditionalSupportingProteinsIDs"
    ].apply(lambda x: "" if pd.isna(x) else x.split(","))

    expression_dfs.append(expression_df)
expression_dfs[0]

# %%
expression_dfs[0].iloc[[0, 1]]


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


# %% jupyter={"source_hidden": true}
# def make_percentile_df(
#     expression_df,
#     first_percentile=10,
#     inclusive_last_percentile=110,
#     percentile_step=10,
#     allowed_algorithms=["Ascending", "Descending"],
# ):

#     gb = expression_df.loc[expression_df["Algorithm"].isin(allowed_algorithms)].groupby(
#         "#Solution"
#     )
#     solutions_expression_dfs = [gb.get_group(x) for x in gb.groups]

#     equal_supp_reads_dfs = []
#     weighted_supp_reads_dfs = []
#     for df in solutions_expression_dfs:
#         equal_df = df.sort_values(
#             "TotalEqualSupportingReads", ascending=False
#         ).reset_index(drop=True)
#         equal_df["CumTotalEqualSupportingReads"] = equal_df[
#             "TotalEqualSupportingReads"
#         ].cumsum()
#         equal_df["%CumTotalEqualSupportingReads"] = (
#             100
#             * equal_df["CumTotalEqualSupportingReads"]
#             / equal_df["TotalEqualSupportingReads"].sum()
#         )
#         equal_supp_reads_dfs.append(equal_df)
#         weighted_df = df.sort_values(
#             "TotalWeightedSupportingReads", ascending=False
#         ).reset_index(drop=True)
#         weighted_df["CumTotalWeightedSupportingReads"] = weighted_df[
#             "TotalWeightedSupportingReads"
#         ].cumsum()
#         weighted_df["%CumTotalWeightedSupportingReads"] = (
#             100
#             * weighted_df["CumTotalWeightedSupportingReads"]
#             / weighted_df["TotalWeightedSupportingReads"].sum()
#         )
#         weighted_supp_reads_dfs.append(weighted_df)

#     # equal_supp_reads_dfs, weighted_supp_reads_dfs = make_supp_reads_dfs(expression_df)

#     solutions = []
#     assignment_methods = []
#     percentiles = []
#     required_proteins = []
#     algorithms = []
#     _conditions = []

#     for dfs, method, col in zip(
#         [equal_supp_reads_dfs, weighted_supp_reads_dfs],
#         ["Equal", "Weighted"],
#         ["%CumTotalEqualSupportingReads", "%CumTotalWeightedSupportingReads"],
#     ):
#         for df in dfs:
#             # ic(df.iloc[:1, :3], method, col)
#             # break
#             solution = df.loc[0, "#Solution"]
#             algorithm = df.loc[0, "Algorithm"]
#             _condition = df.loc[0, condition_col]
#             a = df[col].to_numpy()
#             # for percentile in range(50, 100, 10):
#             # for percentile in range(10, 110, 10):
#             for percentile in range(
#                 first_percentile, inclusive_last_percentile, percentile_step
#             ):
#                 idx = (np.abs(a - percentile)).argmin()
#                 if a[idx] < percentile:
#                     idx += 1
#                 solutions.append(solution)
#                 assignment_methods.append(method)
#                 percentiles.append(percentile)
#                 required_proteins.append(idx)
#                 algorithms.append(algorithm)
#                 _conditions.append(_condition)

#     percentile_df = pd.DataFrame(
#         {
#             "#Solution": solutions,
#             "AssignmentMethod": assignment_methods,
#             "Percentile": percentiles,
#             "RequiredProteins": required_proteins,
#             "Algorithm": algorithms,
#             condition_col: _conditions,
#         }
#     )

#     return percentile_df

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
    find_rand_maximal_solution(
        expression_df, seed, allowed_algorithms=["Ascending", "Descending"]
    )
    for expression_df in expression_dfs
]
maximal_solutions

# %% jupyter={"source_hidden": true}
# percentile_dfs = [
#     make_percentile_df(
#         expression_df.loc[expression_df["#Solution"] == maximal_solution].reset_index(
#             drop=True
#         ),
#         allowed_algorithms=["Ascending", "Descending"],
#     )
#     for expression_df, maximal_solution in zip(expression_dfs, maximal_solutions)
# ]
# percentile_dfs[0]

# %%
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


for assignment_df in assignment_dfs:
    assignment_df["%RelativeExpression"] = (
        100
        * assignment_df["TotalWeightedSupportingReads"]
        / assignment_df["TotalWeightedSupportingReads"].sum()
    )
    assignment_df["%CummulativeRelativeExpression"] = assignment_df[
        "%RelativeExpression"
    ].cumsum()

assignment_dfs[0]

# %% jupyter={"source_hidden": true}
# maximal_dfs = [
#     expression_df.loc[expression_df["#Solution"] == maximal_solution].reset_index(
#         drop=True
#     )
#     for expression_df, maximal_solution in zip(expression_dfs, maximal_solutions)
# ]
# maximal_fraction01_dfs = [
#     fraction01_expression_df.loc[
#         fraction01_expression_df["#Solution"] == maximal_fraction01_solution
#     ].reset_index(drop=True)
#     for fraction01_expression_df, maximal_fraction01_solution in zip(
#         fraction01_expression_dfs, maximal_fraction01_solutions
#     )
# ]

# assignment_dfs = [
#     (
#         maximal_df.sort_values("TotalWeightedSupportingReads", ascending=False)
#         .reset_index(drop=True)
#         .assign(ProteinRank=list(range(1, len(maximal_df) + 1)))
#         .rename(columns={"ProteinRank": "#Protein"})
#     )
#     for maximal_df in maximal_dfs
# ]
# fraction01_assignment_dfs = [
#     (
#         fraction01_maximal_df.sort_values(
#             "TotalWeightedSupportingReads", ascending=False
#         )
#         .reset_index(drop=True)
#         .assign(ProteinRank=list(range(1, len(fraction01_maximal_df) + 1)))
#         .rename(columns={"ProteinRank": "#Protein"})
#     )
#     for fraction01_maximal_df in maximal_fraction01_dfs
# ]
# # fraction01_assignment_dfs[0]

# for assignment_df in assignment_dfs:
#     assignment_df["%RelativeExpression"] = (
#         100
#         * assignment_df["TotalWeightedSupportingReads"]
#         / assignment_df["TotalWeightedSupportingReads"].sum()
#     )
#     assignment_df["%CummulativeRelativeExpression"] = assignment_df[
#         "%RelativeExpression"
#     ].cumsum()
# for assignment_df in fraction01_assignment_dfs:
#     assignment_df["%RelativeExpression"] = (
#         100
#         * assignment_df["TotalWeightedSupportingReads"]
#         / assignment_df["TotalWeightedSupportingReads"].sum()
#     )
#     assignment_df["%CummulativeRelativeExpression"] = assignment_df[
#         "%RelativeExpression"
#     ].cumsum()

# assignment_dfs[0]

# %% jupyter={"source_hidden": true}
# x_axis_name = "Distinct protein rank"
# y_axis_name = "Cummulative relative<br>expression (%)"
# head_title = "Cummulative expression vs. distinct unique proteins"

# cols = min(facet_col_wrap, len(conditions), 4)
# rows = ceil(len(conditions) / cols)
# row_col_iter = list(product(range(1, rows + 1), range(1, cols + 1)))[: len(conditions)]

# # assignment_method = "Weighted"
# percentile_fractions = [0.1, 1.0]

# fig = make_subplots(
#     rows=rows,
#     cols=cols,
#     subplot_titles=conditions,
#     shared_xaxes=True,
#     shared_yaxes=True,
#     x_title=x_axis_name,
#     y_title=y_axis_name,
#     vertical_spacing=facet_row_spacing / 2.5,
#     horizontal_spacing=facet_col_spacing * 1.5,
# )

# for (row, col), assignment_df, fraction01_assignment_df, condition in zip(
#     row_col_iter, assignment_dfs, fraction01_assignment_dfs, conditions
# ):
#     # legend_x = [percentile_df["RequiredProteins"].max() * 5 / 6]
#     # legend_ys = [[25], [20]]

#     for color, df, percentile_fraction in zip(
#         subcolors_discrete_map[condition],
#         [fraction01_assignment_df, assignment_df],
#         percentile_fractions,
#     ):
#         # df = df.loc[df["AssignmentMethod"] == assignment_method]

#         x = df["#Protein"]
#         y = df["%CummulativeRelativeExpression"]

#         # x_mean = df.groupby("Percentile")["RequiredProteins"].apply(np.mean)
#         # y_unique = x_mean.index

#         fig.add_trace(
#             go.Scattergl(
#                 x=x,
#                 y=y,
#                 mode="markers",
#                 # mode="lines+markers",
#                 marker=dict(
#                     color=color,
#                     size=2,
#                     opacity=0.7,
#                     # symbol=symbol,
#                     # line=dict(width=0),
#                 ),
#             ),
#             row=row,
#             col=col,
#         )

# #         fig.add_trace(
# #             go.Scatter(
# #                 x=x_mean,
# #                 y=y_unique,
# #                 mode="lines+markers",
# #                 marker=dict(
# #                     color=color,
# #                     size=7,
# #                     opacity=0.7,
# #                     symbol=symbol,
# #                     # line=dict(width=0),
# #                 ),
# #             ),
# #             row=row,
# #             col=col,
# #         )

# #         fig.add_trace(
# #             go.Scatter(
# #                 x=legend_x,
# #                 y=legend_y,
# #                 mode="markers+text",
# #                 marker=dict(
# #                     color=color,
# #                     size=7,
# #                     opacity=0.5,
# #                     symbol=symbol,
# #                     # line=dict(width=0),
# #                 ),
# #                 text=f" Fraction = {percentile_fraction}",
# #                 textposition="middle right",
# #                 textfont=dict(size=8)
# #             ),
# #             row=row,
# #             col=col,
# #         )

# fig.update_xaxes(
#     # tick0 = -1,
#     # dtick = 5_000,
#     matches="x",
#     type="log",
#     nticks=6,
# )

# fig.update_layout(
#     # title=head_title,
#     showlegend=False,
#     template=template,
#     height=400,
#     width=900,
# )
# # fig.write_image(
# #     f"{head_title} - PacBio.svg",
# #     height=max(300, 200 * rows),
# #     width=max(600, 250 * cols),
# # )
# fig.show()

# %% jupyter={"source_hidden": true}
# x_axis_name = "Distinct protein rank"
# y_axis_name = "Cummulative relative<br>expression (%)"
# head_title = "Cummulative expression vs. distinct unique proteins"

# cols = min(facet_col_wrap, len(conditions), 4)
# rows = ceil(len(conditions) / cols)
# row_col_iter = list(product(range(1, rows + 1), range(1, cols + 1)))[: len(conditions)]

# assignment_method = "Weighted"
# percentile_fractions = [0.1, 1.0]
# symbols = ["circle", "triangle-up"]

# fig = make_subplots(
#     rows=rows,
#     cols=cols,
#     subplot_titles=conditions,
#     shared_xaxes=True,
#     shared_yaxes=True,
#     x_title=x_axis_name,
#     y_title=y_axis_name,
#     vertical_spacing=facet_row_spacing / 2.5,
#     horizontal_spacing=facet_col_spacing * 1.5,
# )

# for (row, col), percentile_df, fraction01_percentile_df, condition in zip(
#     row_col_iter, percentile_dfs, fraction01_percentile_dfs, conditions
# ):
#     legend_x = [percentile_df["RequiredProteins"].max() * 5 / 6]
#     legend_ys = [[25], [20]]

#     for color, symbol, legend_y, df, percentile_fraction in zip(
#         subcolors_discrete_map[condition],
#         symbols,
#         legend_ys,
#         [fraction01_percentile_df, percentile_df],
#         percentile_fractions,
#     ):
#         df = df.loc[df["AssignmentMethod"] == assignment_method]

#         x = df["RequiredProteins"]
#         y = df["Percentile"]

#         x_mean = df.groupby("Percentile")["RequiredProteins"].apply(np.mean)
#         y_unique = x_mean.index

#         fig.add_trace(
#             go.Scatter(
#                 x=x,
#                 y=y,
#                 mode="markers",
#                 marker=dict(
#                     color=color,
#                     size=7,
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
#                     size=7,
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
#                     size=7,
#                     opacity=0.5,
#                     symbol=symbol,
#                     # line=dict(width=0),
#                 ),
#                 text=f" Fraction = {percentile_fraction}",
#                 textposition="middle right",
#                 textfont=dict(size=8),
#             ),
#             row=row,
#             col=col,
#         )

# # fig.update_xaxes(
# #     # tick0 = -1,
# #     # dtick = 5_000,
# #     matches='x',
# #     type="log",
# #     nticks=6
# # )

# fig.update_layout(
#     # title=head_title,
#     showlegend=False,
#     template=template,
#     height=400,
#     width=900,
# )
# # fig.write_image(
# #     f"{head_title} - PacBio.svg",
# #     height=max(300, 200 * rows),
# #     width=max(600, 250 * cols),
# # )
# fig.show()

# %% jupyter={"source_hidden": true}
# x_axis_name = "Distinct protein rank"
# y_axis_name = "Cummulative relative<br>expression (%)"
# head_title = "Weighted cummulative expression vs. distinct protein rank"

# cols = min(facet_col_wrap, len(conditions), 5)
# rows = ceil(len(conditions) / cols)
# row_col_iter = list(product(range(1, rows + 1), range(1, cols + 1)))[: len(conditions)]

# assignment_method = "Weighted"
# y_col_name = "TotalWeightedSupportingReads"

# fig = make_subplots(
#     rows=rows,
#     cols=cols,
#     subplot_titles=conditions,
#     shared_yaxes=True,
#     x_title=x_axis_name,
#     y_title=y_axis_name,
#     # vertical_spacing=facet_row_spacing / 1.5,
#     # horizontal_spacing=facet_col_spacing * 1.5,
# )

# for (row, col), percentile_df, expression_df, maximal_solution, condition in zip(
#     row_col_iter, percentile_dfs, expression_dfs, maximal_solutions, conditions
# ):
#     percentile_df = percentile_df.loc[
#         percentile_df["AssignmentMethod"] == assignment_method
#     ]

#     x = percentile_df["RequiredProteins"]
#     y = percentile_df["Percentile"]

#     # plot percentiles
#     fig.add_trace(
#         go.Scatter(
#             x=x,
#             y=y,
#             mode="markers+lines",
#             marker=dict(
#                 color=color_discrete_map[condition],
#                 size=6,
#                 opacity=0.7,
#             ),
#         ),
#         row=row,
#         col=col,
#     )

#     expression_df = (
#         expression_df.loc[expression_df["#Solution"] == maximal_solution]
#         .sort_values("TotalWeightedSupportingReads", ascending=False)
#         .reset_index(drop=True)
#     )
#     expression_df["CummulativeRelativeWeightedExpression"] = expression_df[
#         ["TotalWeightedSupportingReads"]
#     ].transform(lambda x: 100 * x / x.sum())

#     top_x = [10, 100, 1000]
#     top_y = [
#         expression_df["CummulativeRelativeWeightedExpression"][:x].sum() for x in top_x
#     ]

#     # plot top 10/100/1000 expressed proteins
#     fig.add_trace(
#         go.Scatter(
#             x=top_x,
#             y=top_y,
#             mode="markers+text",
#             marker=dict(
#                 size=6,
#                 symbol="square",
#                 opacity=0.5,
#                 color="black",
#             ),
#             text=[
#                 f"  (10, {top_y[0]:.1f})",
#                 f"   (100, {top_y[1]:.1f})",
#                 f"    (1000, {top_y[2]:.1f})",
#             ],
#             # ways to get better text positioning:
#             # https://community.plotly.com/t/solving-the-problem-of-overlapping-text-labels-in-a-scatterplot-by-manually-assigning-the-position-of-each-label/66159/2
#             # https://github.com/plotly/plotly.py/issues/925
#             textposition="middle right",
#             textfont=dict(size=8),
#         ),
#         row=row,
#         col=col,
#     )

# # fig.update_xaxes(tick0=0, dtick=10_000, matches="x")
# fig.update_layout(
#     title="Squid's PacBio",
#     title_x=0.15,
#     showlegend=False,
#     template=template,
#     height=max(400, 200 * rows),
#     width=max(650, 250 * cols),
# )

# fig.write_image(
#     f"{head_title} - PacBio.svg",
#     #     height=max(300, 200 * rows),
#     #     width=max(600, 250 * cols),
#     height=max(300, 260 * rows),
#     width=max(650, 250 * cols),
# )

# fig.show()

# %%
x_axis_name = "Distinct protein rank"
y_axis_name = "Cummulative relative<br>expression (%)"
head_title = "Weighted cummulative expression vs. distinct protein rank"

cols = min(facet_col_wrap, len(conditions), 5)
rows = ceil(len(conditions) / cols)
row_col_iter = list(product(range(1, rows + 1), range(1, cols + 1)))[: len(conditions)]

fig = make_subplots(
    rows=rows,
    cols=cols,
    subplot_titles=conditions,
    shared_yaxes=True,
    x_title=x_axis_name,
    y_title=y_axis_name,
    # vertical_spacing=facet_row_spacing / 1.5,
    # horizontal_spacing=facet_col_spacing * 1.5,
    vertical_spacing=0.05,
    horizontal_spacing=0.025,
)

for (row, col), assignment_df, condition in zip(
    row_col_iter, assignment_dfs, conditions
):
    
    x = assignment_df["#Protein"]
    y = assignment_df["%CummulativeRelativeExpression"]  

    # plot the cummulative expression
    fig.add_trace(
        go.Scatter(
            x=x,
            y=y,
            # mode="markers",
            mode="lines",
            marker=dict(
                color=color_discrete_map[condition],
                size=1.5,
                opacity=0.7,
            ),
        ),
        row=row,
        col=col,
    )

    top_x = [10, 100, 1000]
    top_y = [
        assignment_df["%CummulativeRelativeExpression"][x].sum() for x in top_x
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

width = max(650, 250 * cols)
height = max(400, 200 * rows)

fig.update_layout(
    title="Squid's PacBio",
    title_x=0.15,
    showlegend=False,
    template=template,
    width=width,
    height=height,
)

fig.write_image(
    f"{head_title} - PacBio.svg",
    width=width,
    height=height,
)

fig.show()


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

# %% jupyter={"source_hidden": true}
# cols = min(facet_col_wrap, len(conditions), 3)
# rows = ceil(len(conditions) / cols)
# row_col_iter = list(product(range(1, rows + 1), range(1, cols + 1)))[: len(conditions)]

# linear_spaces = [(300, 15_000), (250, 23_000)]  # (start, end) tuples for both x and y
# forward_transforms = [
#     (linear_to_log10, linear_to_log10),
#     (linear_to_log10, inverse),
# ]  # (x, y) tuples
# reverse_transforms = [
#     (log10_to_linear, log10_to_linear),
#     (log10_to_linear, inverse),
# ]  # (x, y) tuples
# formulate_equations = [formulate_log10_equation, formulate_semilog10_equation]

# maximal_dfs = [
#     expression_df.loc[expression_df["#Solution"] == maximal_solution].reset_index(
#         drop=True
#     )
#     for expression_df, maximal_solution in zip(expression_dfs, maximal_solutions)
# ]
# maximal_algorithms = [df.loc[0, "Algorithm"] for df in maximal_dfs]

# subplot_titles = [
#     f"{condition}<br><sub>(#{solution}, {algorithm})</sub>"
#     for condition, solution, algorithm in zip(
#         conditions, maximal_solutions, maximal_algorithms
#     )
# ]
# x_axis_name = "Distinct unique protein rank"
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
#     # vertical_spacing=facet_row_spacing / 2.5,
#     horizontal_spacing=facet_col_spacing * 1.5,
# )

# for (
#     (row, col),
#     condition,
#     maximal_df,
#     maximal_solution,
#     maximal_algorithm,
#     linear_space,
#     (forward_x_transform, forward_y_transform),
#     (reverse_x_transform, reverse_y_transform),
#     formulate_equation,
# ) in zip(
#     row_col_iter,
#     conditions,
#     maximal_dfs,
#     maximal_solutions,
#     maximal_algorithms,
#     linear_spaces,
#     forward_transforms,
#     reverse_transforms,
#     formulate_equations,
# ):

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

#         train_logspace = [
#             int(i)
#             for i in np.logspace(
#                 np.log10(linear_space[0]), np.log10(linear_space[1]), num=1000
#             )
#         ]
#         test_logspace = [
#             int(i)
#             for i in np.logspace(
#                 np.log10(linear_space[0] + 20), np.log10(linear_space[1] - 20), num=1000
#             )
#             if int(i) not in train_logspace
#         ]

#         train_x = forward_x_transform(x[train_logspace])
#         train_y = forward_y_transform(y[train_logspace])

#         test_x = forward_x_transform(x[test_logspace])
#         test_y = forward_y_transform(y[test_logspace])

#         # Create linear regression object
#         regr = linear_model.LinearRegression(n_jobs=threads)
#         # Train the model using the training sets
#         regr.fit(np.array(train_x).reshape(-1, 1), train_y)
#         # Make predictions using the testing set
#         pred_y = regr.predict(np.array(test_x).reshape(-1, 1))

#         # transform these variables back to original scale so they can plotted
#         test_x = reverse_x_transform(test_x)
#         pred_y = reverse_y_transform(pred_y)

#         fig.add_trace(
#             fit_scatter_type(
#                 x=test_x,
#                 y=pred_y,
#                 mode="lines",
#                 marker_color=color,
#                 line=dict(
#                     dash="dash",
#                     width=regression_line_width,
#                 ),
#                 legendgroup=condition,
#                 name=f"{assignment_method} - fitted",
#                 showlegend=False,
#             ),
#             row=1,
#             col=col,
#         )

#         coef = regr.coef_[0]
#         intercept = regr.intercept_
#         mse = mean_squared_error(test_y, pred_y)
#         r2 = r2_score(test_y, pred_y)
#         # if intercept >= 0:
#         #     operator = "+"
#         # else:
#         #     operator = "-"
#         #     intercept = np.abs(intercept)
#         equation = formulate_equation(coef, intercept)
#         text = (
#             f"<b>{assignment_method}</b>"
#             "<br>"
#             # f"<b>y = {coef:.2f}x {operator} {intercept:.2f}</b>"
#             # f"y = {coef:.2f}x {operator} {intercept:.2f}"
#             f"{equation}"
#             "<br>"
#             f"MSE = {mse:.2f}"  # 0 is perfect prediction
#             "<br>"
#             f"R2 = {r2:.2f}"  # 1 is perfect prediction
#         )

#         if assignment_method == assignment_methods[0]:
#             # textposition = "top right"
#             i = int(len(test_x) / 10)
#             text_x = test_x.iloc[i] + 2000
#             text_y = pred_y[i] + 0.03
#         else:
#             # textposition = "bottom left"
#             i = int(len(test_x) / 3.5)
#             text_x = test_x.iloc[0] - int(3 * train_logspace[0] / 4)
#             text_y = pred_y[i] - 0.003
#         text_x = np.log10(text_x)
#         text_y = np.log10(text_y)

#         fig.add_annotation(
#             row=row,
#             col=col,
#             x=text_x,
#             y=text_y,
#             xref="x",
#             yref="y",
#             text=text,
#             align="center",
#             font=dict(size=8, color=color),
#             showarrow=False,
#         )

# fig.update_layout(
#     title_text=head_title,
#     title_y=0.95,
#     template=template,
#     showlegend=False,
#     # legend_itemsizing="constant",
#     height=max(400, 200 * rows),
#     # width=max(900, 250 * cols),
# )
# fig.update_xaxes(type="log")
# fig.update_yaxes(type="log")
# fig.write_image(
#     f"{head_title} - PacBio.svg",
#     height=max(350, 200 * rows),
#     width=max(650, 350 * cols),
# )
# fig.show()
# # fig.show(config={'staticPlot': True, 'responsive': False})


# %%
# assignment_method = "Weighted"
# y_col_name = "TotalWeightedSupportingReads"

cols = min(facet_col_wrap, len(conditions), 3)
rows = ceil(len(conditions) / cols)
row_col_iter = list(product(range(1, rows + 1), range(1, cols + 1)))[: len(conditions)]

# (start, end) tuples for both x and y
linear_spaces = [
    (300, 15_000),
    # (250, 23_000)
    (700, 15_000),
]

forward_transforms = [
    (linear_to_log10, linear_to_log10),
    # (linear_to_log10, inverse),
    (linear_to_log10, linear_to_log10),
]  # (x, y) tuples
reverse_transforms = [
    (log10_to_linear, log10_to_linear),
    # (log10_to_linear, inverse),
    (log10_to_linear, log10_to_linear),
]  # (x, y) tuples

formulate_equations = [
    formulate_log10_equation,
    # formulate_semilog10_equation
    formulate_log10_equation,
]
# fit_texts = ["    y ~ 1 / sqrt(x)", "    y ~ 1 / log(x)"]

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

fig = make_subplots(
    rows=rows,
    cols=cols,
    y_title=y_axis_name,
    x_title=x_axis_name,
    subplot_titles=subplot_titles,
    shared_yaxes=True,
    # shared_xaxes=True,
    # # vertical_spacing=facet_row_spacing / 2.5,
    # horizontal_spacing=facet_col_spacing * 1.5,
    vertical_spacing=0.05,
    horizontal_spacing=0.025,
)

for (
    (row, col),
    condition,
    # maximal_df,
    # maximal_solution,
    assignment_df,
    linear_space,
    (forward_x_transform, forward_y_transform),
    (reverse_x_transform, reverse_y_transform),
    # fit_text,
    formulate_equation,
) in zip(
    row_col_iter,
    conditions,
    # maximal_dfs,
    # maximal_solutions,
    assignment_dfs,
    linear_spaces,
    forward_transforms,
    reverse_transforms,
    # fit_texts,
    formulate_equations,
):

    x = assignment_df["#Protein"]
    y = assignment_df["%RelativeExpression"]

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

    train_x = forward_x_transform(x[train_logspace])
    train_y = forward_y_transform(y[train_logspace])

    test_x = forward_x_transform(x[test_logspace])
    test_y = forward_y_transform(y[test_logspace])

    # Create linear regression object
    regr = linear_model.LinearRegression(n_jobs=threads)
    # Train the model using the training sets
    regr.fit(np.array(train_x).reshape(-1, 1), train_y)
    # Make predictions using the testing set
    pred_y = regr.predict(np.array(test_x).reshape(-1, 1))

    # transform these variables back to original scale so that they can plotted
    test_x = reverse_x_transform(test_x)  # should be equivalent to `test_x = x[test_logspace]`?
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
            # legendgroup=condition,
            # name=f"{assignment_method} - fitted",
            showlegend=False,
        ),
        row=1,
        col=col,
    )

    i = int(len(test_x) / 10)
    text_x = test_x.iloc[i] + 2000
    text_y = pred_y[i] + 0.03
    # text_x = 1000
    # text_y = 0.05
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

    ic(condition)
    ic(equation)

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
        font=dict(size=10, color="grey"),
        showarrow=False,
    )

fig.update_layout(
    title_text="Squid's PacBio",
    title_x=0.11,
    # title_y=0.95,
    template=template,
    showlegend=False,
    # legend_itemsizing="constant",
    height=max(400, 200 * rows),
    width=max(900, 250 * cols),
)
fig.update_xaxes(type="log", nticks=6)
fig.update_yaxes(type="log")
fig.write_image(
    f"{head_title} - PacBio.svg",
    height=max(400, 200 * rows),
    # width=max(650, 250 * cols),
    width=max(830, 250 * cols),
)
fig.show()
# fig.show(config={'staticPlot': True, 'responsive': False})


# %% jupyter={"source_hidden": true}
# assignment_method = "Weighted"
# y_col_name = "TotalWeightedSupportingReads"


# maximal_dfs = [
#     expression_df.loc[expression_df["#Solution"] == maximal_solution].reset_index(
#         drop=True
#     )
#     for expression_df, maximal_solution in zip(expression_dfs, maximal_solutions)
# ]

# assignment_dfs = [
#     (
#         maximal_df.sort_values("TotalWeightedSupportingReads", ascending=False)
#         .reset_index(drop=True)
#         .assign(ProteinRank=list(range(1, len(maximal_df) + 1)))
#         .rename(columns={"ProteinRank": "#Protein"})
#     )
#     for maximal_df in maximal_dfs
# ]

# for assignment_df in assignment_dfs:
#     assignment_df["%RelativeExpression"] = 100 * assignment_df[y_col_name] / assignment_df[y_col_name].sum()

# assignment_dfs[0]

# %% [markdown]
# ##### Relative expression of isoforms - 0.1 vs 1

# %%
fraction01_expression_dfs = []
for fraction01_expression_file in fraction01_expression_files:
    fraction01_expression_df = pd.read_csv(fraction01_expression_file, sep=sep)
    fraction01_expression_df["#Solution"] = fraction01_expression_df[
        "#Solution"
    ].astype(str)
    # fraction01_expression_df["Diff5+"] = (
    #     abs(
    #         fraction01_expression_df["TotalEqualSupportingReads"]
    #         - fraction01_expression_df["TotalWeightedSupportingReads"]
    #     )
    #     >= 0.05
    #     * (
    #         fraction01_expression_df["TotalEqualSupportingReads"]
    #         + fraction01_expression_df["TotalWeightedSupportingReads"]
    #     )
    #     / 2
    # )

    fraction01_expression_df["AdditionalSupportingReadsIDs"] = fraction01_expression_df[
        "AdditionalSupportingReadsIDs"
    ].apply(lambda x: "" if pd.isna(x) else [y.split(",") for y in x.split(";")])
    fraction01_expression_df[
        "AdditionalSupportingProteinsIDs"
    ] = fraction01_expression_df["AdditionalSupportingProteinsIDs"].apply(
        lambda x: "" if pd.isna(x) else x.split(",")
    )

    fraction01_expression_dfs.append(fraction01_expression_df)
    
fraction01_expression_dfs[0]

# %%
fraction01_expression_df.iloc[[1000, 2000]]

# %%
maximal_fraction01_solutions = [
    find_rand_maximal_solution(
        expression_df, seed, allowed_algorithms=["Ascending", "Descending"]
    )
    for expression_df in fraction01_expression_dfs
]
maximal_fraction01_solutions

# %% jupyter={"source_hidden": true}
# fraction01_percentile_dfs = [
#     make_percentile_df(
#         expression_df.loc[expression_df["#Solution"] == maximal_solution].reset_index(
#             drop=True
#         ),
#         allowed_algorithms=["Ascending", "Descending"],
#     )
#     for expression_df, maximal_solution in zip(
#         fraction01_expression_dfs, maximal_fraction01_solutions
#     )
# ]
# fraction01_percentile_dfs[0]

# %%
# maximal_dfs = [
#     expression_df.loc[expression_df["#Solution"] == maximal_solution].reset_index(
#         drop=True
#     )
#     for expression_df, maximal_solution in zip(expression_dfs, maximal_solutions)
# ]
maximal_fraction01_dfs = [
    fraction01_expression_df.loc[
        fraction01_expression_df["#Solution"] == maximal_fraction01_solution
    ].reset_index(drop=True)
    for fraction01_expression_df, maximal_fraction01_solution in zip(
        fraction01_expression_dfs, maximal_fraction01_solutions
    )
]

# assignment_dfs = [
#     (
#         maximal_df.sort_values("TotalWeightedSupportingReads", ascending=False)
#         .reset_index(drop=True)
#         .assign(ProteinRank=list(range(1, len(maximal_df) + 1)))
#         .rename(columns={"ProteinRank": "#Protein"})
#     )
#     for maximal_df in maximal_dfs
# ]
fraction01_assignment_dfs = [
    (
        fraction01_maximal_df.sort_values(
            "TotalWeightedSupportingReads", ascending=False
        )
        .reset_index(drop=True)
        .assign(ProteinRank=list(range(1, len(fraction01_maximal_df) + 1)))
        .rename(columns={"ProteinRank": "#Protein"})
    )
    for fraction01_maximal_df in maximal_fraction01_dfs
]
# fraction01_assignment_dfs[0]

# for assignment_df in assignment_dfs:
#     assignment_df["%RelativeExpression"] = (
#         100
#         * assignment_df["TotalWeightedSupportingReads"]
#         / assignment_df["TotalWeightedSupportingReads"].sum()
#     )
#     assignment_df["%CummulativeRelativeExpression"] = assignment_df[
#         "%RelativeExpression"
#     ].cumsum()
for assignment_df in fraction01_assignment_dfs:
    assignment_df["%RelativeExpression"] = (
        100
        * assignment_df["TotalWeightedSupportingReads"]
        / assignment_df["TotalWeightedSupportingReads"].sum()
    )
    assignment_df["%CummulativeRelativeExpression"] = assignment_df[
        "%RelativeExpression"
    ].cumsum()

fraction01_assignment_dfs[0]

# %%
x_axis_name = "Distinct protein rank"
y_axis_name = "Cummulative relative<br>expression (%)"
head_title = "Cummulative expression vs. distinct unique proteins"

cols = min(facet_col_wrap, len(conditions), 4)
rows = ceil(len(conditions) / cols)
row_col_iter = list(product(range(1, rows + 1), range(1, cols + 1)))[: len(conditions)]

# assignment_method = "Weighted"
percentile_fractions = [0.1, 1.0]

fig = make_subplots(
    rows=rows,
    cols=cols,
    subplot_titles=conditions,
    shared_xaxes=True,
    shared_yaxes=True,
    x_title=x_axis_name,
    y_title=y_axis_name,
    vertical_spacing=facet_row_spacing / 2.5,
    horizontal_spacing=facet_col_spacing * 1.5,
)

for (row, col), assignment_df, fraction01_assignment_df, condition in zip(
    row_col_iter, assignment_dfs, fraction01_assignment_dfs, conditions
):
    # legend_x = [percentile_df["RequiredProteins"].max() * 5 / 6]
    # legend_ys = [[25], [20]]
    legend_x = [10]
    legend_ys = [[75], [65]]

    for color, df, percentile_fraction, legend_y in zip(
        subcolors_discrete_map[condition],
        [fraction01_assignment_df, assignment_df],
        percentile_fractions,
        legend_ys
    ):
        # df = df.loc[df["AssignmentMethod"] == assignment_method]

        x = df["#Protein"]
        y = df["%CummulativeRelativeExpression"]       

        # x_mean = df.groupby("Percentile")["RequiredProteins"].apply(np.mean)
        # y_unique = x_mean.index

        fig.add_trace(
            go.Scattergl(
                x=x,
                y=y,
                mode="markers",
                # mode="lines+markers",
                marker=dict(
                    color=color,
                    size=2,
                    opacity=0.7,
                    # symbol=symbol,
                    # line=dict(width=0),
                ),
            ),
            row=row,
            col=col,
        )

        fig.add_trace(
            go.Scatter(
                x=legend_x,
                y=legend_y,
                mode="markers+text",
                marker=dict(
                    color=color,
                    size=7,
                    opacity=0.7,
                    # symbol=symbol,
                    # line=dict(width=0),
                ),
                text=f"  Fraction = {percentile_fraction}",
                textposition="middle right",
                textfont=dict(size=9)
            ),
            row=row,
            col=col,
        )

fig.update_xaxes(
    # tick0 = -1,
    # dtick = 5_000,
    matches="x",
    type="log",
    nticks=6,
)

fig.update_layout(
    title="Squid's PacBio",
    title_x=0.11,
    showlegend=False,
    template=template,
    height=400,
    width=900,
)

fig.write_image(
    f"Cummulative expression vs. distinct protein rank: fractions 0.1 and 1 - PacBio.svg",
    height=400,
    width=900,
)

fig.show()

# %% jupyter={"source_hidden": true}
# x_axis_name = "Distinct protein rank"
# y_axis_name = "Cummulative relative<br>expression (%)"
# head_title = "Cummulative expression vs. distinct unique proteins"

# cols = min(facet_col_wrap, len(conditions), 4)
# rows = ceil(len(conditions) / cols)
# row_col_iter = list(product(range(1, rows + 1), range(1, cols + 1)))[: len(conditions)]

# assignment_method = "Weighted"
# percentile_fractions = [0.1, 1.0]
# symbols = ["circle", "triangle-up"]

# fig = make_subplots(
#     rows=rows,
#     cols=cols,
#     subplot_titles=conditions,
#     shared_xaxes=True,
#     shared_yaxes=True,
#     x_title=x_axis_name,
#     y_title=y_axis_name,
#     vertical_spacing=facet_row_spacing / 2.5,
#     horizontal_spacing=facet_col_spacing * 1.5,
# )

# for (row, col), percentile_df, fraction01_percentile_df, condition in zip(
#     row_col_iter, percentile_dfs, fraction01_percentile_dfs, conditions
# ):
#     legend_x = [percentile_df["RequiredProteins"].max() * 5 / 6]
#     legend_ys = [[25], [20]]

#     for color, symbol, legend_y, df, percentile_fraction in zip(
#         subcolors_discrete_map[condition],
#         symbols,
#         legend_ys,
#         [fraction01_percentile_df, percentile_df],
#         percentile_fractions,
#     ):
#         df = df.loc[df["AssignmentMethod"] == assignment_method]

#         x = df["RequiredProteins"]
#         y = df["Percentile"]

#         x_mean = df.groupby("Percentile")["RequiredProteins"].apply(np.mean)
#         y_unique = x_mean.index

#         fig.add_trace(
#             go.Scatter(
#                 x=x,
#                 y=y,
#                 mode="markers",
#                 marker=dict(
#                     color=color,
#                     size=7,
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
#                     size=7,
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
#                     size=7,
#                     opacity=0.5,
#                     symbol=symbol,
#                     # line=dict(width=0),
#                 ),
#                 text=f" Fraction = {percentile_fraction}",
#                 textposition="middle right",
#                 textfont=dict(size=8),
#             ),
#             row=row,
#             col=col,
#         )

# # fig.update_xaxes(
# #     # tick0 = -1,
# #     # dtick = 5_000,
# #     matches='x',
# #     type="log",
# #     nticks=6
# # )

# fig.update_layout(
#     # title=head_title,
#     showlegend=False,
#     template=template,
#     height=400,
#     width=900,
# )
# # fig.write_image(
# #     f"{head_title} - PacBio.svg",
# #     height=max(300, 200 * rows),
# #     width=max(600, 250 * cols),
# # )

# fig.show()

# %%
fraction01_assignment_dfs[0]

# %%
frac_1_merged_to_frac_01_assignment_dfs = []

frac_01_merged_to_frac_1_assignment_dfs = []

for fraction01_assignment_df, assignment_df in zip(
    fraction01_assignment_dfs, assignment_dfs
):

    frac_1_merged_to_frac_01_assignment_df = fraction01_assignment_df.loc[
        :,
        [
            "Protein",
            "#Protein",
            "TotalWeightedSupportingReads",
            "%RelativeExpression",
            "AdditionalSupportingReadsIDs",
            "AdditionalSupportingProteinsIDs",
            "AdditionalSupportingProteins",
        ],
    ].merge(
        assignment_df.loc[
            :,
            [
                "Protein",
                "#Protein",
                "TotalWeightedSupportingReads",
            "%RelativeExpression",
                "AdditionalSupportingReadsIDs",
                "AdditionalSupportingProteinsIDs",
                "AdditionalSupportingProteins",
            ],
        ],
        on="Protein",
        how="left",
        suffixes=["0.1", "1.0"],
    )
    # frac_1_merged_to_frac_01_assignment_df["%RelativeExpression1.0"] = (
    #     100
    #     * frac_1_merged_to_frac_01_assignment_df["TotalWeightedSupportingReads1.0"]
    #     / frac_1_merged_to_frac_01_assignment_df[
    #         "TotalWeightedSupportingReads1.0"
    #     ].sum()
    # )
    # frac_1_merged_to_frac_01_assignment_df["%RelativeExpression0.1"] = (
    #     100
    #     * frac_1_merged_to_frac_01_assignment_df["TotalWeightedSupportingReads0.1"]
    #     / frac_1_merged_to_frac_01_assignment_df[
    #         "TotalWeightedSupportingReads0.1"
    #     ].sum()
    # )

    frac_01_merged_to_frac_1_assignment_df = assignment_df.loc[
        :,
        [
            "Protein",
            "#Protein",
            "TotalWeightedSupportingReads",
            "%RelativeExpression",
            "AdditionalSupportingReadsIDs",
            "AdditionalSupportingProteinsIDs",
            "AdditionalSupportingProteins",
        ],
    ].merge(
        fraction01_assignment_df.loc[
            :,
            [
                "Protein",
                "#Protein",
                "TotalWeightedSupportingReads",
            "%RelativeExpression",
                "AdditionalSupportingReadsIDs",
                "AdditionalSupportingProteinsIDs",
                "AdditionalSupportingProteins",
            ],
        ],
        on="Protein",
        how="left",
        suffixes=["1.0", "0.1"],
    )
    # frac_01_merged_to_frac_1_assignment_df["%RelativeExpression1.0"] = (
    #     100
    #     * frac_01_merged_to_frac_1_assignment_df["TotalWeightedSupportingReads1.0"]
    #     / frac_01_merged_to_frac_1_assignment_df[
    #         "TotalWeightedSupportingReads1.0"
    #     ].sum()
    # )
    # frac_01_merged_to_frac_1_assignment_df["%RelativeExpression0.1"] = (
    #     100
    #     * frac_01_merged_to_frac_1_assignment_df["TotalWeightedSupportingReads0.1"]
    #     / frac_01_merged_to_frac_1_assignment_df[
    #         "TotalWeightedSupportingReads0.1"
    #     ].sum()
    # )
    
    frac_1_merged_to_frac_01_assignment_dfs.append(
        frac_1_merged_to_frac_01_assignment_df
    )
    frac_01_merged_to_frac_1_assignment_dfs.append(
        frac_01_merged_to_frac_1_assignment_df
    )

frac_1_merged_to_frac_01_assignment_dfs[0]

# %%
frac_01_merged_to_frac_1_assignment_dfs[0]

# %%
frac_1_merged_to_frac_01_assignment_dfs[0].loc[:, ["%RelativeExpression0.1", "%RelativeExpression1.0"]].sum()

# %%
frac_01_merged_to_frac_1_assignment_dfs[0].loc[:, ["%RelativeExpression0.1", "%RelativeExpression1.0"]].sum()

# %%

# %%
x_axis_name = "Relative expression (%): fraction 0.1"
y_axis_name = "Relative expression (%): fraction 1"

cols = len(conditions)
rows = 1

relative_assignment_dfs = [
    frac_1_merged_to_frac_01_assignment_dfs,
    frac_01_merged_to_frac_1_assignment_dfs,
]

fig = make_subplots(
    rows=rows,
    cols=cols,
    # column_titles=["Fraction 1 relative to 0.1", "Fraction 0.1 relative to 1"],
    # row_titles=conditions,
    subplot_titles=conditions,
    shared_xaxes=True,
    shared_yaxes=True,
    x_title=x_axis_name,
    y_title=y_axis_name,
    # vertical_spacing=facet_row_spacing / 2.5,
    # horizontal_spacing=facet_col_spacing * 1.5,
)

for col, condition in zip(range(cols), conditions):

    frac_1_relative_to_frac_01_df = relative_assignment_dfs[0][col]
    frac_01_relative_to_frac_1_df = relative_assignment_dfs[1][col]

    frac_1_relative_to_frac_01_x = frac_1_relative_to_frac_01_df[
        "%RelativeExpression0.1"
    ]
    frac_1_relative_to_frac_01_y = frac_1_relative_to_frac_01_df[
        "%RelativeExpression1.0"
    ]

    frac_01_relative_to_frac_1_x = frac_01_relative_to_frac_1_df[
        "%RelativeExpression0.1"
    ]
    frac_01_relative_to_frac_1_y = frac_01_relative_to_frac_1_df[
        "%RelativeExpression1.0"
    ]

    fig.add_trace(
        go.Scattergl(
            x=frac_1_relative_to_frac_01_x,
            y=frac_1_relative_to_frac_01_y,
            mode="markers",
            marker=dict(
                color=color_discrete_map[condition],
                size=3,
                # opacity=0.7,
                # symbol=symbol,
            ),
        ),
        row=1,
        col=col + 1,
    )

    frac_1_matches = frac_1_relative_to_frac_01_y.notna().sum()
    possible_frac_1_matches = len(frac_1_relative_to_frac_01_y)

    frac_01_matches = frac_01_relative_to_frac_1_x.notna().sum()
    possible_frac_01_matches = len(frac_01_relative_to_frac_1_x)

    # ic(col, row, condition, possible_matches, matches)

    fig.add_annotation(
        row=1,
        col=col + 1,
        # x=0.5,
        # y=0.35,
        x=1.5,
        y=0.05,
        xref="x",
        yref="y",
        # text=f"{frac_1_matches} / {possible_frac_1_matches} = {100 * frac_1_matches / possible_frac_1_matches:.2f}%<br>of possible matches<br>relative to fraction 0.1",
        # text=f"{frac_1_matches} / {possible_frac_1_matches} = {100 * frac_1_matches / possible_frac_1_matches:.2f}%<br>of possible matches<br>relative to fraction 1",
        text=f"{frac_1_matches} / {possible_frac_1_matches} = {100 * frac_1_matches / possible_frac_1_matches:.2f}%<br>of possible matches<br>to fraction 1 proteins",
        # text=f"{frac_1_matches} / {possible_frac_1_matches} = {100 * frac_1_matches / possible_frac_1_matches:.2f}%<br>of fraction 0.1 proteins<br>were matched to<br>fraction 1 proteins",
        align="center",
        font=dict(size=8, color="black"),
        showarrow=False,
        # bgcolor="#ff7f0e",
        bgcolor="white",
    )

    fig.add_annotation(
        row=1,
        col=col + 1,
        # x=1.5,
        # y=0.05,
        x=0.5,
        y=0.35,
        xref="x",
        yref="y",
        # text=f"{frac_01_matches} / {possible_frac_01_matches} = {100 * frac_01_matches / possible_frac_01_matches:.2f}%<br>of possible matches<br>relative to fraction 1",
        # text=f"{frac_01_matches} / {possible_frac_01_matches} = {100 * frac_01_matches / possible_frac_01_matches:.2f}%<br>of possible matches<br>relative to fraction 0.1",
        text=f"{frac_01_matches} / {possible_frac_01_matches} = {100 * frac_01_matches / possible_frac_01_matches:.2f}%<br>of possible matches<br>to fraction 0.1 proteins",
        align="center",
        # font=dict(size=10, color="grey"),
        font=dict(size=9, color="black"),
        showarrow=False,
        # bgcolor="#ff7f0e",
        bgcolor="white",
    )

width = 800
height = 400
    
fig.update_layout(
    title="Squid's PacBio",
    title_x=0.12,
    showlegend=False,
    template=template,
    width=width,
    height=height,
)

fig.write_image(
    "Fraction 0.1 expression vs. fraction 1 expression - PacBio.svg",
    width=width,
    height=height,
)

fig.show()

# %% jupyter={"source_hidden": true}
# frac_1_merged_to_frac_01_assignment_dfs[0]

# %% jupyter={"source_hidden": true}
# frac_1_merged_to_frac_01_assignment_dfs[0].loc[:, ["#Protein0.1", "#Protein1.0"]]

# %% jupyter={"source_hidden": true}
# frac_01_merged_to_frac_1_assignment_dfs[0].loc[:, ["#Protein0.1", "#Protein1.0"]]

# %% jupyter={"source_hidden": true}
# top_10_gria_assignment_df = (
#     frac_1_merged_to_frac_01_assignment_dfs[0].dropna().iloc[:10]
# )
# # top_10_gria_assignment_df = frac_1_merged_to_frac_01_assignment_dfs[0].dropna().copy()

# # top_10_gria_assignment_df

# top_10_gria_assignment_df[
#     "MutualAdditionalSupportingProteinsIDs"
# ] = top_10_gria_assignment_df.apply(
#     lambda x: list(
#         set(x["AdditionalSupportingProteinsIDs0.1"])
#         & set(x["AdditionalSupportingProteinsIDs1.0"])
#     ),
#     axis=1,
# )
# top_10_gria_assignment_df["%MutualAdditionalSupportingProteinsIDs0.1"] = (
#     100
#     * top_10_gria_assignment_df["MutualAdditionalSupportingProteinsIDs"].apply(len)
#     / top_10_gria_assignment_df["AdditionalSupportingProteinsIDs0.1"].apply(len)
# )
# top_10_gria_assignment_df["%MutualAdditionalSupportingProteinsIDs1.0"] = (
#     100
#     * top_10_gria_assignment_df["MutualAdditionalSupportingProteinsIDs"].apply(len)
#     / top_10_gria_assignment_df["AdditionalSupportingProteinsIDs1.0"].apply(len)
# )

# import itertools

# top_10_gria_assignment_df[
#     "MutualAdditionalSupportingReadsIDs"
# ] = top_10_gria_assignment_df.apply(
#     lambda x: list(
#         set(itertools.chain.from_iterable(x["AdditionalSupportingReadsIDs0.1"]))
#         & set(itertools.chain.from_iterable(x["AdditionalSupportingReadsIDs1.0"]))
#     ),
#     axis=1,
# )

# top_10_gria_assignment_df["%MutualAdditionalSupportinReadsIDs0.1"] = (
#     100
#     * top_10_gria_assignment_df["MutualAdditionalSupportingReadsIDs"].apply(len)
#     / top_10_gria_assignment_df["AdditionalSupportingReadsIDs0.1"].apply(
#         lambda x: len(list(itertools.chain.from_iterable(x)))
#     )
# )
# top_10_gria_assignment_df["%MutualAdditionalSupportingReadsIDs1.0"] = (
#     100
#     * top_10_gria_assignment_df["MutualAdditionalSupportingReadsIDs"].apply(len)
#     / top_10_gria_assignment_df["AdditionalSupportingReadsIDs1.0"].apply(
#         lambda x: len(list(itertools.chain.from_iterable(x)))
#     )
# )

# top_10_gria_assignment_df

# %% jupyter={"source_hidden": true}
# fig = px.scatter(
#     top_10_gria_assignment_df,
#     x="%MutualAdditionalSupportinReadsIDs0.1",
#     y="%MutualAdditionalSupportingReadsIDs1.0",
#     marginal_y="box",
#     marginal_x="box",
# )
# fig.update_traces(marker={"size": 5})
# fig.update_layout(width=600, height=600, template=template)
# fig.show()

# %% jupyter={"source_hidden": true}
# fig = px.scatter(
#     top_10_gria_assignment_df,
#     x="%MutualAdditionalSupportingProteinsIDs0.1",
#     y="%MutualAdditionalSupportingProteinsIDs1.0",
#     marginal_y="box",
#     marginal_x="box",
# )
# fig.update_traces(marker={"size": 5})
# fig.update_layout(width=600, height=600, template=template)
# fig.show()

# %%
x_axis_name = "Distinct protein rank: fraction 0.1"
y_axis_name = "Distinct protein rank: fraction 1"

cols = len(conditions)
rows = 1

relative_assignment_dfs = [
    frac_1_merged_to_frac_01_assignment_dfs,
    frac_01_merged_to_frac_1_assignment_dfs,
]

fig = make_subplots(
    rows=rows,
    cols=cols,
    # column_titles=["Fraction 1 relative to 0.1", "Fraction 0.1 relative to 1"],
    # row_titles=conditions,
    subplot_titles=conditions,
    shared_xaxes=True,
    shared_yaxes=True,
    x_title=x_axis_name,
    y_title=y_axis_name,
    # vertical_spacing=facet_row_spacing / 2.5,
    # horizontal_spacing=facet_col_spacing * 1.5,
)

for col, condition in zip(range(cols), conditions):

    frac_1_relative_to_frac_01_df = relative_assignment_dfs[0][col]
    frac_01_relative_to_frac_1_df = relative_assignment_dfs[1][col]

    frac_1_relative_to_frac_01_x = frac_1_relative_to_frac_01_df["#Protein0.1"]
    frac_1_relative_to_frac_01_y = frac_1_relative_to_frac_01_df["#Protein1.0"]

    frac_01_relative_to_frac_1_x = frac_01_relative_to_frac_1_df["#Protein0.1"]
    frac_01_relative_to_frac_1_y = frac_01_relative_to_frac_1_df["#Protein1.0"]

    fig.add_trace(
        go.Scattergl(
            x=frac_1_relative_to_frac_01_x,
            y=frac_1_relative_to_frac_01_y,
            mode="markers",
            marker=dict(
                color=color_discrete_map[condition],
                size=3,
                # opacity=0.7,
                # symbol=symbol,
            ),
        ),
        row=1,
        col=col + 1,
    )

    frac_1_matches = frac_1_relative_to_frac_01_y.notna().sum()
    possible_frac_1_matches = len(frac_1_relative_to_frac_01_y)

    frac_01_matches = frac_01_relative_to_frac_1_x.notna().sum()
    possible_frac_01_matches = len(frac_01_relative_to_frac_1_x)

    # ic(col, row, condition, possible_matches, matches)

    fig.add_annotation(
        row=1,
        col=col + 1,
        x=2_700,
        y=3_000,
        xref="x",
        yref="y",
        text=f"{frac_1_matches} / {possible_frac_1_matches} = {100 * frac_1_matches / possible_frac_1_matches:.2f}%<br>of possible matches<br>relative to fraction 0.1",
        align="center",
        font=dict(size=8, color="black"),
        showarrow=False,
        # bgcolor="#ff7f0e",
        bgcolor="white",
    )

    fig.add_annotation(
        row=1,
        col=col + 1,
        x=1000,
        y=23_000,
        xref="x",
        yref="y",
        text=f"{frac_01_matches} / {possible_frac_01_matches} = {100 * frac_01_matches / possible_frac_01_matches:.2f}%<br>of possible matches<br>relative to fraction 1",
        align="center",
        # font=dict(size=10, color="grey"),
        font=dict(size=8, color="grey"),
        showarrow=False,
        # bgcolor="#ff7f0e",
        bgcolor="white",
    )


# fig.update_xaxes(
#     tick0 = 0,
#     dtick = 5_000,
#     matches='x'
# )
fig.update_layout(
    # title=head_title,
    showlegend=False,
    template=template,
    height=430,
    width=700,
)
# fig.write_image(
#     f"{head_title} - PacBio.svg",
#     height=max(300, 200 * rows),
#     width=max(600, 250 * cols),
# )

fig.show()

# %% jupyter={"source_hidden": true}
# assignment_method = "Weighted"
# y_col_name = "TotalWeightedSupportingReads"

# cols = min(facet_col_wrap, len(conditions), 3)
# rows = ceil(len(conditions) / cols)
# row_col_iter = list(product(range(1, rows + 1), range(1, cols + 1)))[: len(conditions)]


# subplot_titles = conditions
# x_axis_name = "Distinct protein rank"
# y_axis_name = "Relative expression (%)"
# # head_title = f"Relative expression of proteins considering a largest solution in each {str(condition_col).lower()}"
# head_title = "Relative expression vs. distinct protein rank"

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
#     # shared_xaxes=True,
#     # # vertical_spacing=facet_row_spacing / 2.5,
#     # horizontal_spacing=facet_col_spacing * 1.5,
#     vertical_spacing=0.05,
#     horizontal_spacing=0.025,
# )

# for (
#     (row, col),
#     condition,
#     maximal_df,
#     maximal_solution,
#     linear_space,
#     (forward_x_transform, forward_y_transform),
#     (reverse_x_transform, reverse_y_transform),
#     # fit_text,
#     formulate_equation,
# ) in zip(
#     row_col_iter,
#     conditions,
#     fraction01_assignment_dfs,
#     maximal_fraction01_solutions,
#     linear_spaces,
#     forward_transforms,
#     reverse_transforms,
#     # fit_texts,
#     formulate_equations,
# ):

#     assignment_df = maximal_df.sort_values(
#         "TotalWeightedSupportingReads", ascending=False
#     ).reset_index(drop=True)
#     assignment_df["#Protein"] = list(range(1, len(assignment_df) + 1))
#     assignment_df["AssignmentMethod"] = assignment_method

#     x = assignment_df["#Protein"]
#     y = 100 * assignment_df[y_col_name] / assignment_df[y_col_name].sum()

#     fig.add_trace(
#         data_scatter_type(
#             x=x,
#             y=y,
#             # legendgrouptitle_text=condition,
#             # legendgroup=condition,
#             # name=assignment_method,
#             mode="markers",
#             marker_color=color_discrete_map[condition],
#             marker_size=data_marker_size,
#             marker=dict(
#                 opacity=data_opacity,
#                 line=dict(width=0),
#             ),
#         ),
#         row=row,
#         col=col,
#     )

# #     train_logspace = [
# #         int(i)
# #         for i in np.logspace(
# #             np.log10(linear_space[0]), np.log10(linear_space[1]), num=1000
# #         )
# #     ]
# #     test_logspace = [
# #         int(i)
# #         for i in np.logspace(
# #             np.log10(linear_space[0] + 20), np.log10(linear_space[1] - 20), num=1000
# #         )
# #         if int(i) not in train_logspace
# #     ]

# #     train_x = forward_x_transform(x[train_logspace])
# #     train_y = forward_y_transform(y[train_logspace])

# #     test_x = forward_x_transform(x[test_logspace])
# #     test_y = forward_y_transform(y[test_logspace])

# #     # Create linear regression object
# #     regr = linear_model.LinearRegression(n_jobs=threads)
# #     # Train the model using the training sets
# #     regr.fit(np.array(train_x).reshape(-1, 1), train_y)
# #     # Make predictions using the testing set
# #     pred_y = regr.predict(np.array(test_x).reshape(-1, 1))

# #     # transform these variables back to original scale so that they can plotted
# #     test_x = reverse_x_transform(test_x)  # should be equivalent to `test_x = x[test_logspace]`?
# #     pred_y = reverse_y_transform(pred_y)

# #     fig.add_trace(
# #         fit_scatter_type(
# #             x=test_x,
# #             y=pred_y,
# #             mode="lines",
# #             marker_color="grey",
# #             line=dict(
# #                 dash="dash",
# #                 width=regression_line_width,
# #             ),
# #             # legendgroup=condition,
# #             # name=f"{assignment_method} - fitted",
# #             showlegend=False,
# #         ),
# #         row=1,
# #         col=col,
# #     )

# #     i = int(len(test_x) / 10)
# #     text_x = test_x.iloc[i] + 2000
# #     text_y = pred_y[i] + 0.03
# #     # text_x = 1000
# #     # text_y = 0.05
# #     text_x = np.log10(text_x)
# #     text_y = np.log10(text_y)

# #     coef = regr.coef_[0]
# #     intercept = regr.intercept_
# #     # mse = mean_squared_error(test_y, pred_y)
# #     # r2 = r2_score(test_y, pred_y)
# #     # if intercept >= 0:
# #     #     operator = "+"
# #     # else:
# #     #     operator = "-"
# #     #     intercept = np.abs(intercept)
# #     equation = formulate_equation(coef, intercept)

# #     ic(condition)
# #     ic(equation)

# #     fit_text_new = (
# #         # f"<b>{assignment_method}</b>"
# #         # "<br>"
# #         # f"<b>y = {coef:.2f}x {operator} {intercept:.2f}</b>"
# #         # f"y = {coef:.2f}x {operator} {intercept:.2f}"
# #         f"{equation}"
# #         # "<br>"
# #         # f"MSE = {mse:.2f}"  # 0 is perfect prediction
# #         # "<br>"
# #         # f"R2 = {r2:.2f}"  # 1 is perfect prediction
# #     )

# #     fig.add_annotation(
# #         row=row,
# #         col=col,
# #         x=text_x,
# #         y=text_y,
# #         xref="x",
# #         yref="y",
# #         # text=fit_text,
# #         text=fit_text_new,
# #         align="center",
# #         font=dict(size=10, color="grey"),
# #         showarrow=False,
# #     )

# fig.update_layout(
#     # title_text=head_title,
#     # title_y=0.95,
#     template=template,
#     showlegend=False,
#     # legend_itemsizing="constant",
#     height=max(400, 200 * rows),
#     width=max(900, 250 * cols),
# )
# fig.update_xaxes(type="log", nticks=6)
# fig.update_yaxes(type="log")
# fig.write_image(
#     f"{head_title} - PacBio.svg",
#     height=max(400, 200 * rows),
#     # width=max(650, 250 * cols),
#     width=max(830, 250 * cols),
# )
# fig.show()
# # fig.show(config={'staticPlot': True, 'responsive': False})


# %% jupyter={"source_hidden": true}
# assignment_method = "Weighted"
# y_col_name = "TotalWeightedSupportingReads"

# cols = min(facet_col_wrap, len(conditions), 3)
# rows = ceil(len(conditions) / cols)
# row_col_iter = list(product(range(1, rows + 1), range(1, cols + 1)))[: len(conditions)]

# # (start, end) tuples for both x and y
# linear_spaces = [
#     (300, 15_000),
#     # (250, 23_000)
#     (700, 15_000),
# ]

# forward_transforms = [
#     (linear_to_log10, linear_to_log10),
#     # (linear_to_log10, inverse),
#     (linear_to_log10, linear_to_log10),
# ]  # (x, y) tuples
# reverse_transforms = [
#     (log10_to_linear, log10_to_linear),
#     # (log10_to_linear, inverse),
#     (log10_to_linear, log10_to_linear),
# ]  # (x, y) tuples

# formulate_equations = [
#     formulate_log10_equation,
#     # formulate_semilog10_equation
#     formulate_log10_equation,
# ]
# # fit_texts = ["    y ~ 1 / sqrt(x)", "    y ~ 1 / log(x)"]

# # fraction01_expression_dfs
# # maximal_fraction01_solutions

# maximal_fraction01_dfs = [
#     fraction01_expression_df.loc[
#         fraction01_expression_df["#Solution"] == maximal_fraction01_solution
#     ].reset_index(drop=True)
#     for fraction01_expression_df, maximal_fraction01_solution in zip(
#         fraction01_expression_dfs, maximal_fraction01_solutions
#     )
# ]

# fraction01_assignment_dfs = [
#     (
#         fraction01_maximal_df.sort_values(
#             "TotalWeightedSupportingReads", ascending=False
#         )
#         .reset_index(drop=True)
#         .assign(ProteinRank=list(range(1, len(fraction01_maximal_df) + 1)))
#         .rename(columns={"ProteinRank": "#Protein"})
#     )
#     for fraction01_maximal_df in maximal_fraction01_dfs
# ]

# subplot_titles = conditions
# x_axis_name = "Distinct protein rank"
# y_axis_name = "Relative expression (%)"
# # head_title = f"Relative expression of proteins considering a largest solution in each {str(condition_col).lower()}"
# head_title = "Relative expression vs. distinct protein rank"

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
#     # shared_xaxes=True,
#     # # vertical_spacing=facet_row_spacing / 2.5,
#     # horizontal_spacing=facet_col_spacing * 1.5,
#     vertical_spacing=0.05,
#     horizontal_spacing=0.025,
# )

# for (
#     (row, col),
#     condition,
#     maximal_df,
#     maximal_solution,
#     linear_space,
#     (forward_x_transform, forward_y_transform),
#     (reverse_x_transform, reverse_y_transform),
#     # fit_text,
#     formulate_equation,
# ) in zip(
#     row_col_iter,
#     conditions,
#     fraction01_assignment_dfs,
#     maximal_fraction01_solutions,
#     linear_spaces,
#     forward_transforms,
#     reverse_transforms,
#     # fit_texts,
#     formulate_equations,
# ):

#     assignment_df = maximal_df.sort_values(
#         "TotalWeightedSupportingReads", ascending=False
#     ).reset_index(drop=True)
#     assignment_df["#Protein"] = list(range(1, len(assignment_df) + 1))
#     assignment_df["AssignmentMethod"] = assignment_method

#     x = assignment_df["#Protein"]
#     y = 100 * assignment_df[y_col_name] / assignment_df[y_col_name].sum()

#     fig.add_trace(
#         data_scatter_type(
#             x=x,
#             y=y,
#             # legendgrouptitle_text=condition,
#             # legendgroup=condition,
#             # name=assignment_method,
#             mode="markers",
#             marker_color=color_discrete_map[condition],
#             marker_size=data_marker_size,
#             marker=dict(
#                 opacity=data_opacity,
#                 line=dict(width=0),
#             ),
#         ),
#         row=row,
#         col=col,
#     )

# #     train_logspace = [
# #         int(i)
# #         for i in np.logspace(
# #             np.log10(linear_space[0]), np.log10(linear_space[1]), num=1000
# #         )
# #     ]
# #     test_logspace = [
# #         int(i)
# #         for i in np.logspace(
# #             np.log10(linear_space[0] + 20), np.log10(linear_space[1] - 20), num=1000
# #         )
# #         if int(i) not in train_logspace
# #     ]

# #     train_x = forward_x_transform(x[train_logspace])
# #     train_y = forward_y_transform(y[train_logspace])

# #     test_x = forward_x_transform(x[test_logspace])
# #     test_y = forward_y_transform(y[test_logspace])

# #     # Create linear regression object
# #     regr = linear_model.LinearRegression(n_jobs=threads)
# #     # Train the model using the training sets
# #     regr.fit(np.array(train_x).reshape(-1, 1), train_y)
# #     # Make predictions using the testing set
# #     pred_y = regr.predict(np.array(test_x).reshape(-1, 1))

# #     # transform these variables back to original scale so that they can plotted
# #     test_x = reverse_x_transform(test_x)  # should be equivalent to `test_x = x[test_logspace]`?
# #     pred_y = reverse_y_transform(pred_y)

# #     fig.add_trace(
# #         fit_scatter_type(
# #             x=test_x,
# #             y=pred_y,
# #             mode="lines",
# #             marker_color="grey",
# #             line=dict(
# #                 dash="dash",
# #                 width=regression_line_width,
# #             ),
# #             # legendgroup=condition,
# #             # name=f"{assignment_method} - fitted",
# #             showlegend=False,
# #         ),
# #         row=1,
# #         col=col,
# #     )

# #     i = int(len(test_x) / 10)
# #     text_x = test_x.iloc[i] + 2000
# #     text_y = pred_y[i] + 0.03
# #     # text_x = 1000
# #     # text_y = 0.05
# #     text_x = np.log10(text_x)
# #     text_y = np.log10(text_y)

# #     coef = regr.coef_[0]
# #     intercept = regr.intercept_
# #     # mse = mean_squared_error(test_y, pred_y)
# #     # r2 = r2_score(test_y, pred_y)
# #     # if intercept >= 0:
# #     #     operator = "+"
# #     # else:
# #     #     operator = "-"
# #     #     intercept = np.abs(intercept)
# #     equation = formulate_equation(coef, intercept)

# #     ic(condition)
# #     ic(equation)

# #     fit_text_new = (
# #         # f"<b>{assignment_method}</b>"
# #         # "<br>"
# #         # f"<b>y = {coef:.2f}x {operator} {intercept:.2f}</b>"
# #         # f"y = {coef:.2f}x {operator} {intercept:.2f}"
# #         f"{equation}"
# #         # "<br>"
# #         # f"MSE = {mse:.2f}"  # 0 is perfect prediction
# #         # "<br>"
# #         # f"R2 = {r2:.2f}"  # 1 is perfect prediction
# #     )

# #     fig.add_annotation(
# #         row=row,
# #         col=col,
# #         x=text_x,
# #         y=text_y,
# #         xref="x",
# #         yref="y",
# #         # text=fit_text,
# #         text=fit_text_new,
# #         align="center",
# #         font=dict(size=10, color="grey"),
# #         showarrow=False,
# #     )

# fig.update_layout(
#     # title_text=head_title,
#     # title_y=0.95,
#     template=template,
#     showlegend=False,
#     # legend_itemsizing="constant",
#     height=max(400, 200 * rows),
#     width=max(900, 250 * cols),
# )
# fig.update_xaxes(type="log", nticks=6)
# fig.update_yaxes(type="log")
# fig.write_image(
#     f"{head_title} - PacBio.svg",
#     height=max(400, 200 * rows),
#     # width=max(650, 250 * cols),
#     width=max(830, 250 * cols),
# )
# fig.show()
# # fig.show(config={'staticPlot': True, 'responsive': False})


# %% [markdown]
# ##### Clustering

# %%
# raise Exception("don't run clustering right now...")

# %%
# # max_sol_dfs = [
# #     expression_df.loc[expression_df["#Solution"] == maximal_solution].reset_index(
# #         drop=True
# #     )
# #     for expression_df, maximal_solution in zip(expression_dfs, maximal_solutions)
# # ]
# # max_sol_dfs[0]

# # max_sol_dfs == maximal_dfs

# maximal_dfs[0]

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


# %% jupyter={"source_hidden": true}
# if impute_nans:
#     mean_col_editing_freqs_wo_nan = df_b.apply(np.mean)
#     df_b = df_b.apply(
#         lambda row: [
#             update_cell(
#                 cell,
#                 original_aa,
#                 0.0,
#                 imputed_nan_value,
#                 1.0,
#                 check_unchanged_func=lambda cell, original_aa: cell == 0,
#                 check_nan_func=lambda cell: np.isnan(cell),
#             )
#             for cell, original_aa, imputed_nan_value in zip(
#                 row, original_aas, mean_col_editing_freqs_wo_nan
#             )
#         ],
#         axis=1,
#         result_type="broadcast",
#     )

# df_b = drop_uniformative_aa_cols(df_b)

# df = pd.concat([df_a, df_b], axis=1)

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
        # "Diff5+",
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
ML_INPUT_FIRST_COL_POS_NEW = 5


def prepare_ml_input_df_new(
    assignment_df,
    unique_proteins_df,
    unique_proteins_first_col_pos,
    relative_expression_sorting_col="%RelativeExpression", # relative expression based on weighted supporting reads
    impute_nans=True,
):
    """Prepare an input df fit for some machine learning algorithms such as tSNE, PCA, etc."""

    cols_to_use_from_assignment_df = [
        condition_col,
        "Protein",
        "#Solution",
        "Algorithm",
        relative_expression_sorting_col
    ]
    cols_to_use_from_unique_proteins_df = [
        condition_col,
        "Protein",
    ] + unique_proteins_df.columns[unique_proteins_first_col_pos:].to_list()

    df = assignment_df.loc[:, cols_to_use_from_assignment_df].merge(
        unique_proteins_df.loc[:, cols_to_use_from_unique_proteins_df],
        how="left",
        on=[condition_col, "Protein"],
    )

    df_a = df.iloc[:, :len(cols_to_use_from_assignment_df)]

    original_aas = [col.split("(")[1][0] for col in df.columns[len(cols_to_use_from_assignment_df):]]

    df_b = df.iloc[:, len(cols_to_use_from_assignment_df):].apply(
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

    df = df.sort_values(relative_expression_sorting_col, ascending=False, ignore_index=True)

    return df


# %%
prepare_ml_input_df_new(
    assignment_dfs[0],
    unique_proteins_dfs[0],
    unique_proteins_first_col_pos,
)

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
# # equal_exp_tsne_input_dfs[0]

# weighted_exp_tsne_input_dfs = [
#     prepare_ml_input_df(
#         max_sol_df,
#         unique_proteins_df,
#         unique_proteins_first_col_pos,
#         sorting_col="%WeightedTotalExpression",
#     )
#     for max_sol_df, unique_proteins_df in zip(max_sol_dfs, unique_proteins_dfs)
# ]


weighted_exp_tsne_input_dfs = [
    prepare_ml_input_df_new(
        assignment_df,
        unique_proteins_df,
        unique_proteins_first_col_pos,
    )
    for assignment_df, unique_proteins_df in zip(assignment_dfs, unique_proteins_dfs)
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


# %% jupyter={"source_hidden": true}
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
def run_pcas(
    conditions,
    pca_input_dfs,
    seed,
    n_components=2,
    # first_col_pos=ML_INPUT_FIRST_COL_POS,
    first_col_pos=ML_INPUT_FIRST_COL_POS_NEW,
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
def run_tsnes(
    conditions,
    tsne_input_dfs,
    seed,
    perplexities=[5, 30, 50, 100, 150],
    n_components=2,
    learning_rate="auto",
    n_iter=300,
    init="random",
    n_jobs=20,
    # first_col_pos=ML_INPUT_FIRST_COL_POS,
    first_col_pos=ML_INPUT_FIRST_COL_POS_NEW,
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


# %% [markdown]
# > All proteins

# %%

# %%

# %% [markdown]
# > TEST START

# %%
# weighted_exp_tsne_input_dfs[0].iloc[:, ML_INPUT_FIRST_COL_POS:]

# %%
# # perplexities = [5, 30, 50, 100]
# # perplexities = [150]
# n_iter = 500
# n_jobs = 20
# n_iter_500_weighted_conditions_tsnes, n_iter_500_weighted_conditions_Xs = run_tsnes(
#     [conditions[0]],
#     [weighted_exp_tsne_input_dfs[0]],
#     seed,
#     perplexities=[150],
#     n_iter=n_iter,
#     n_jobs=n_jobs,
# )

# %%
# prots_perplexity_tsne = n_iter_500_weighted_conditions_tsnes[0][0]
# X_transformed = n_iter_500_weighted_conditions_Xs[0]

# %%
# X = weighted_exp_tsne_input_dfs[0].iloc[:, ML_INPUT_FIRST_COL_POS:].values
# X

# %%
# # perplexities = [5, 30, 50, 100]
# # perplexities = [150]
# n_iter = 500
# n_jobs = 20
# top_expressed_proteins = 1000
# top_1000_n_iter_500_weighted_conditions_tsnes, top_1000_n_iter_500_weighted_conditions_Xs = run_tsnes(
#     [conditions[0]],
#     [weighted_exp_tsne_input_dfs[0]],
#     seed,
#     perplexities=[150],
#     n_iter=n_iter,
#     n_jobs=n_jobs,
#     top_expressed_proteins=top_expressed_proteins
# )

# %%
# top_1000_prots_perplexity_tsne = top_1000_n_iter_500_weighted_conditions_tsnes[0][0]
# top_1000_X_transformed = top_1000_n_iter_500_weighted_conditions_Xs[0]

# %%
# top_1000_X = weighted_exp_tsne_input_dfs[0].iloc[:top_expressed_proteins, ML_INPUT_FIRST_COL_POS:].values
# top_1000_X

# %%
# # cluster_sizes = [2**x for x in range(2, 8)]
# cluster_sizes = [2**x for x in range(2, 5)]
# cluster_sizes

# %%
# kmodes_clusters = [
#     KModes(n_clusters=n_clusters, init='Huang', n_init=5, verbose=1, n_jobs=n_jobs).fit_predict(X)
#     for n_clusters in cluster_sizes
# ]

# %%
# cluster_sizes = [2**x for x in range(2, 8)]
# # cluster_sizes = [2**x for x in range(2, 5)]
# cluster_sizes

# %%
# list(range(10, 110, 10))

# %%
# # cluster_sizes = [2**x for x in range(2, 8)]
# cluster_sizes = list(range(10, 110, 10)) # 10, 20, ..., 100
# top_1000_kmodes = [
#     KModes(n_clusters=n_clusters, init='Huang', n_init=5, verbose=0, n_jobs=n_jobs).fit(top_1000_X)
#     for n_clusters in cluster_sizes
# ]
# top_1000_kmodes_sumofsq = [km.cost_ for km in top_1000_kmodes]
# top_1000_kmodes_clusters = [
#     km.predict(top_1000_X) for km in top_1000_kmodes
# ]

# %%
# cluster_sizes_2 = list(range(110, 160, 10)) # 110, 120, ..., 150
# top_1000_kmodes_2 = [
#     KModes(n_clusters=n_clusters, init='Huang', n_init=5, verbose=0, n_jobs=n_jobs).fit(top_1000_X)
#     for n_clusters in cluster_sizes_2
# ]
# top_1000_kmodes_sumofsq_2 = [km.cost_ for km in top_1000_kmodes_2]
# top_1000_kmodes_clusters_2 = [
#     km.predict(top_1000_X) for km in top_1000_kmodes_2
# ]

# %%
# fig = px.scatter(
#     x=cluster_sizes+cluster_sizes_2,
#     y=top_1000_kmodes_sumofsq+top_1000_kmodes_sumofsq_2,
#     template=template,
#     labels={
#         "x": "Number of clusters (k)",
#         "y": "Sum of square distances"
#     },
#     title="Elbow method for optimal number of clusters",
#     # log_x=True
#     # log_y=True
# )
# fig.update_layout(
#     width=700,
#     height=500
# )
# fig.update_yaxes(range=[0, max(top_1000_kmodes_sumofsq)*1.1])
# fig.show()

# %%
# len(set(top_1000_kmodes_clusters[-1]))

# %%
# fig = go.Figure()
# x, y = top_1000_prots_perplexity_tsne.T
# fig.add_trace(
#     go.Scattergl(
#         x=x, y=y,
#         mode="markers",
#         marker=dict(
#             # color="white",
#             color=top_1000_kmodes_clusters[-1],
#             line_width=line_width,
#             # size=marker_size
#         )
#     ),
#     # row=row,
#     # col=col
# )
# fig.update_layout(
#     width=700,
#     height=500,
#     template=template
# )
# # fig.update_yaxes(range=[0, max(top_1000_kmodes_sumofsq)*1.1])
# fig.show()

# %%
# cols = min(len(cluster_sizes), 5)
# rows = ceil(len(cluster_sizes) / cols)
# row_col_iter = list(product(range(1, rows + 1), range(1, cols + 1)))[:len(cluster_sizes)]

# head_title = (
#     "t-SNE for top 1000 expressed proteins in GRIA, colors by kmodes of different n_clusters"
# )
# # column_titles = cluster_sizes

# fig = make_subplots(
#     rows=rows,
#     cols=cols,
#     subplot_titles=[f"K = {k}" for k in cluster_sizes],
#     # shared_yaxes=True,
#     # shared_xaxes=True,
#     # x_title="PC1",
#     # y_title="PC2",
# )

# # color_sequence = px.colors.qualitative.Light24
# x, y = top_1000_prots_perplexity_tsne.T

# for (row, col), cluster, cluster_size in zip(row_col_iter, top_1000_kmodes_clusters, cluster_sizes):

#     # colormap = {label: color for label, color in zip(range(cluster_size), color_sequence[:cluster_size])}
#     # colors = [colormap[label] for label in cluster]

#     # rank_cutoff = 1000
#     # fig = go.Figure()
#     marker_size = 3
#     line_width = 0
#     # n = X.shape[0]
#     # color_options = [color_discrete_map[condition], "white"]
#     # colors = color_highest_expressed_proteins(n, rank_cutoff, color_options)

#     fig.add_trace(
#         go.Scattergl(
#             x=x, y=y,
#             mode="markers",
#             marker=dict(
#                 # color="white",
#                 # color=colors,
#                 color=cluster,
#                 colorscale="Electric",
#                 line_width=line_width,
#                 size=marker_size
#             )
#         ),
#         row=row,
#         col=col
#     )
# pixles = 800
# fig.update_layout(
#     title_text=head_title,
#     # title_y=0.95,
#     template=template,
#     showlegend=False,
#     # width=600,
#     # height=1800,
#     width=1400,
#     height=600,
# )

# fig.show()

# %%
# cols = min(len(clusters), 3)
# rows = ceil(len(clusters) / cols)
# row_col_iter = list(product(range(1, rows + 1), range(1, cols + 1)))[:len(clusters)]

# head_title = (
#     "t-SNE for GRIA, colors by kmodes of different n_clusters"
# )
# # column_titles = cluster_sizes

# fig = make_subplots(
#     rows=rows,
#     cols=cols,
#     subplot_titles=cluster_sizes,
#     # shared_yaxes=True,
#     # shared_xaxes=True,
#     # x_title="PC1",
#     # y_title="PC2",
# )

# color_sequence = px.colors.qualitative.Light24
# x, y = prots_perplexity_tsne.T

# for (row, col), cluster, cluster_size in zip(row_col_iter, kmodes_clusters, cluster_sizes):

#     colormap = {label: color for label, color in zip(range(cluster_size), color_sequence[:cluster_size])}
#     colors = [colormap[label] for label in cluster]

#     # rank_cutoff = 1000
#     # fig = go.Figure()
#     marker_size = 2
#     line_width = 0.5
#     # n = X.shape[0]
#     # color_options = [color_discrete_map[condition], "white"]
#     # colors = color_highest_expressed_proteins(n, rank_cutoff, color_options)

#     fig.add_trace(
#         go.Scattergl(
#             x=x, y=y,
#             mode="markers",
#             marker=dict(
#                 # color="white",
#                 color=colors,
#                 line_width=line_width, size=marker_size
#             )
#         ),
#         row=row,
#         col=col
#     )
# pixles = 800
# fig.update_layout(
#     title_text=head_title,
#     # title_y=0.95,
#     template=template,
#     showlegend=False,
#     # width=600,
#     # height=1800,
#     # width=1800,
#     # height=600,
# )

# fig.show()

# %%
# test_tsne = TSNE(
#     n_components=2,
#     learning_rate="auto",
#     perplexity=150,
#     n_iter=300,
#     init="random",
#     random_state=np.random.RandomState(seed),
#     n_jobs=20,
# )
# test_tsne = test_tsne.fit_transform(test_X)

# %%
# km = KModes(n_clusters=4, init='Huang', n_init=5, verbose=1)

# clusters = km.fit_predict(test_X)

# # Print the cluster centroids
# print(km.cluster_centroids_)

# %%
# test_colormap = {label: color for label, color in zip(range(4), ["red", "blue", "green", "yellow"])}
# colors = [test_colormap[label] for label in km.labels_]
# colors[:6]

# %%
# # prots_perplexity_tsne = n_iter_500_weighted_conditions_tsnes[0][0]
# # X = n_iter_500_weighted_conditions_Xs[0][0]

# # rank_cutoff = 1000
# fig = go.Figure()
# # marker_size = 1
# # line_width = 0.5
# # n = X.shape[0]
# # color_options = [color_discrete_map[condition], "white"]
# # colors = color_highest_expressed_proteins(n, rank_cutoff, color_options)
# x, y = test_tsne.T
# fig.add_trace(
#     go.Scattergl(
#         x=x, y=y,
#         mode="markers",
#         marker=dict(
#             color=colors,
#             # line_width=line_width,
#             # size=marker_size
#         )
#     ),
# )
# fig.update_layout(
#     # title_text=head_title,
#     # title_y=0.95,
#     template=template,
#     # showlegend=False,
#     width=550,
#     height=550,
# )

# fig.show()

# %%
# km.labels_

# %%
# len(km.labels_)

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

# %% [markdown]
# > TEST END

# %%

# %%

# %%

# %%

# %%
# perplexities = [5, 30, 50, 100]
# perplexities = [5, 30, 50, 100, 150]
perplexities = [150]
n_iter = 500
n_jobs = 40

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
cols = min(facet_col_wrap, len(conditions), 5)
rows = ceil(len(conditions) / cols)
row_col_iter = list(product(range(1, rows + 1), range(1, cols + 1)))[: len(conditions)]

fig = make_subplots(
    rows=rows,
    cols=cols,
    subplot_titles=conditions,
    # shared_yaxes=True,
    x_title="t-SNE 1",
    y_title="t-SNE 2",
)

# min_x = None
# max_x = 0
# max_y = 0

marker_size = 1
line_width = 0.5


for (row, col), condition, X, condition_tsnes in zip(
    row_col_iter,
    conditions,
    n_iter_500_weighted_conditions_Xs,
    n_iter_500_weighted_conditions_tsnes,
):

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
    title_text="Squid's PacBio",
    title_x=0.19,
    # title_y=0.95,
    template=template,
    showlegend=False,
    width=600,
    height=350,
)

fig.write_image("tSNE clustering - PacBio.svg", width=600, height=350)

fig.show()


# %%

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
#             go.Scattergl(
#                 x=x,
#                 y=y,
#                 mode="markers",
#                 marker=dict(color=colors, line_width=line_width, size=marker_size),
#             ),
#             row=row,
#             col=col,
#         )

# fig.update_layout(
#     title_text=head_title,
#     title_y=0.95,
#     template=template,
#     showlegend=False,
#     width=1200,
#     height=600,
# )

# fig.write_image("tSNE clustering - PacBio.svg", width=1200, height=600)

# fig.show()

# %%

# %%

# %%
# # equal_conditions_pcas = run_pcas(conditions, equal_exp_tsne_input_dfs, seed)
# weighted_conditions_pcas = run_pcas(conditions, weighted_exp_tsne_input_dfs, seed)

# %%
# rank_cutoff = 1000

# cols = min(facet_col_wrap, len(conditions), 4)
# rows = ceil(len(conditions) / cols)
# row_col_iter = list(product(range(1, rows + 1), range(1, cols + 1)))[: len(conditions)]

# head_title = (
#     f"PCAs for largest solution of each {str(condition_col).lower()}"
#     "<br>"
#     f"<sub>{rank_cutoff} highest expressed proteins are colored</sub>"
# )
# column_titles = conditions

# fig = make_subplots(
#     rows=rows,
#     cols=cols,
#     subplot_titles=conditions,
#     # shared_yaxes=True,
#     # shared_xaxes=True,
#     x_title="PC1",
#     y_title="PC2",
# )

# marker_size = 1

# for (row, col), condition, condition_pca in zip(
#     row_col_iter, conditions, weighted_conditions_pcas
# ):

#     n = len(condition_pca)
#     color_options = [color_discrete_map[condition], "white"]
#     colors = color_highest_expressed_proteins(n, rank_cutoff, color_options)

#     x = condition_pca["PC1"]
#     y = condition_pca["PC2"]

#     fig.add_trace(
#         go.Scattergl(
#             x=x, y=y, mode="markers", marker=dict(color=colors, line_width=0.5, size=marker_size)
#         ),
#         row=1,
#         col=col,
#     )

# fig.update_layout(
#     title_text=head_title,
#     title_y=0.95,
#     template=template,
#     showlegend=False,
#     width=800,
#     height=500,
# )

# fig.show()

# %% [markdown]
# > Top 1000 expressed proteins

# %%
# # perplexities = [5, 30, 50, 100]
# perplexities = [5, 30, 50, 100, 150, 200]

# %%
# top_expressed_proteins = 1000
# # top_1000_equal_conditions_tsnes, top_1000_equal_conditions_Xs = run_tsnes(
# #     conditions,
# #     equal_exp_tsne_input_dfs,
# #     seed,
# #     perplexities=perplexities,
# #     top_expressed_proteins=top_expressed_proteins,
# # )
# top_1000_weighted_conditions_tsnes, top_1000_weighted_conditions_Xs = run_tsnes(
#     conditions,
#     weighted_exp_tsne_input_dfs,
#     seed,
#     perplexities=perplexities,
#     top_expressed_proteins=top_expressed_proteins,
# )

# %%
# rank_cutoff = 1000

# # head_title = (
# #     f"t-SNEs for largest solution of each {str(condition_col).lower()} under different perplexities"
# #     "<br>"
# #     f"<sub>{rank_cutoff} highest expressed proteins are colored</sub>"
# # )
# head_title = (
#     f"t-SNEs for top {top_expressed_proteins} expressed proteins in largest solution of each {str(condition_col).lower()} under different perplexities"
#     "<br>"
#     f"<sub>{rank_cutoff} highest expressed proteins (according to % of total {sorting_method} expression) are colored</sub>"
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

# for row, (condition, X, condition_tsnes) in enumerate(
#     zip(conditions, top_1000_weighted_conditions_Xs, top_1000_weighted_conditions_tsnes), start=1
# ):

#     n = X.shape[0]
#     color_options = [color_discrete_map[condition], "white"]
#     colors = color_highest_expressed_proteins(n, rank_cutoff, color_options)

#     for col, prots_perplexity_tsne in enumerate(condition_tsnes, start=1):

#         x, y = prots_perplexity_tsne.T

#         fig.add_trace(
#             go.Scattergl(
#                 x=x, y=y, mode="markers", marker=dict(color=colors, line_width=0.5, size=marker_size)
#             ),
#             row=row,
#             col=col,
#         )

# fig.update_layout(
#     title_text=head_title,
#     title_y=0.95,
#     template=template,
#     showlegend=False,
#     width=1200,
#     height=600,
# )

# fig.show()

# %%
# top_expressed_proteins = 1000
# # top_1000_equal_conditions_pcas = run_pcas(
# #     conditions,
# #     equal_exp_tsne_input_dfs,
# #     seed,
# #     top_expressed_proteins=top_expressed_proteins,
# # )
# top_1000_weighted_conditions_pcas = run_pcas(
#     conditions,
#     weighted_exp_tsne_input_dfs,
#     seed,
#     top_expressed_proteins=top_expressed_proteins,
# )

# %%
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

# %% jupyter={"source_hidden": true}
# max_sol_dfs = [
#     expression_df.loc[expression_df["#Solution"] == maximal_solution].reset_index(
#         drop=True
#     )
#     for expression_df, maximal_solution in zip(expression_dfs, maximal_solutions)
# ]
# max_sol_dfs[0]

# %% jupyter={"source_hidden": true}
# def calc_data_entropy(max_sol_exp_df, prcnt_equal_exp_col, prcnt_weighted_exp_col):
#     def _calc_data_entropy(prcnt_exp_col):
#         p = max_sol_exp_df[prcnt_exp_col] / 100
#         s_data = sum(-p * np.log2(p))
#         return s_data

#     s_data_equal_exp = _calc_data_entropy(prcnt_equal_exp_col)

#     s_data_weighted_exp = _calc_data_entropy(prcnt_weighted_exp_col)

#     return s_data_equal_exp, s_data_weighted_exp

# %% jupyter={"source_hidden": true}
# def calc_entropies(
#     max_sol_exp_df,
#     first_col_pos,
#     prcnt_equal_exp_col="%EqualTotalExpression",
#     prcnt_weighted_exp_col="%WeightedTotalExpression",
# ):
#     s_data_equal_exp, s_data_weighted_exp = calc_data_entropy(
#         max_sol_exp_df, prcnt_equal_exp_col, prcnt_weighted_exp_col
#     )
#     s_hypo = calc_hypothetical_entropy(max_sol_exp_df, first_col_pos)
#     return s_data_equal_exp, s_data_weighted_exp, s_hypo

# %% jupyter={"source_hidden": true}
# def calc_hypothetical_entropy(max_sol_exp_df, first_col_pos):
#     p = max_sol_exp_df.iloc[:, first_col_pos:].apply(np.mean)
#     s_hypo = sum(-p * np.log2(p) - (1 - p) * np.log2((1 - p)))
#     return s_hypo

# %% jupyter={"source_hidden": true}
# max_sol_exp_dfs = [
#     prepare_ml_input_df(
#         max_sol_df,
#         unique_proteins_df,
#         unique_proteins_first_col_pos,
#         sorting_col="%WeightedTotalExpression",
#         impute_nans=False,
#     )
#     for max_sol_df, unique_proteins_df in zip(max_sol_dfs, unique_proteins_dfs)
# ]

# max_sol_entropies = [
#     calc_entropies(
#         max_sol_exp_df,
#         ML_INPUT_FIRST_COL_POS,
#         "%EqualTotalExpression",
#         "%WeightedTotalExpression",
#     )
#     for max_sol_exp_df in max_sol_exp_dfs
# ]

# entropies_names = ["EqualExpressionData", "WeightedExpressionData", "Hypothetical"]

# _conditions = []
# _entropies_names = []
# _entropies = []
# for condition, condition_entropies in zip(conditions, max_sol_entropies):
#     for entropy, entropy_name in zip(condition_entropies, entropies_names):
#         _conditions.append(condition)
#         _entropies_names.append(entropy_name)
#         _entropies.append(entropy)

# shannon_df = pd.DataFrame(
#     {
#         condition_col: _conditions,
#         "EntropyName": _entropies_names,
#         "EntropyValue": _entropies,
#     }
# )

# shannon_df

# %%
def calc_data_entropy_new(assignment_df, prcnt_exp_col):
    def _calc_data_entropy(prcnt_exp_col):
        p = assignment_df[prcnt_exp_col] / 100
        s_data = sum(-p * np.log2(p))
        return s_data

    s_data = _calc_data_entropy(prcnt_exp_col)

    return s_data


# %%
def calc_hypothetical_entropy_new(assignment_df, first_col_pos):
    p = assignment_df.iloc[:, first_col_pos:].apply(np.mean)
    s_hypo = sum(-p * np.log2(p) - (1 - p) * np.log2((1 - p)))
    return s_hypo


# %%
def calc_entropies_new(
    assignment_df,
    first_col_pos,
    prcnt_exp_col="%RelativeExpression",
):
    s_data = calc_data_entropy_new(assignment_df, prcnt_exp_col)
    s_hypo = calc_hypothetical_entropy_new(assignment_df, first_col_pos)
    return s_data, s_hypo


# %%
max_sol_exp_dfs_new = [
    prepare_ml_input_df_new(
        assignment_df,
        unique_proteins_df,
        unique_proteins_first_col_pos,
        relative_expression_sorting_col="%RelativeExpression", # relative expression based on weighted supporting reads
        impute_nans=False,
    )
    for assignment_df, unique_proteins_df in zip(assignment_dfs, unique_proteins_dfs)
]

max_sol_entropies = [
    calc_entropies_new(
        max_sol_exp_df_new,
        ML_INPUT_FIRST_COL_POS_NEW,
    )
    for max_sol_exp_df_new in max_sol_exp_dfs_new
]

entropies_names = [
    # "EqualExpressionData", 
    "WeightedExpressionData", "Hypothetical"
]

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
fig = go.Figure()

colors = [color_discrete_map[condition] for condition in conditions]

non_syns_per_max_sol_exp_df = [
    str(max_sol_exp_df_new.iloc[:, ML_INPUT_FIRST_COL_POS_NEW:].shape[1])
    for max_sol_exp_df_new in max_sol_exp_dfs_new
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
        visible="legendonly",
        marker=dict(
            color="white",
            pattern_shape="/",
            line=dict(
                # color="grey",
                color="black",
            ),
        ),
        legendgroup="Observed",
        name="Observed",
        showlegend=True,
    )
)
fig.add_trace(
    go.Bar(
        x=[None],
        y=[None],
        visible="legendonly",
        marker=dict(
            color="white",
            pattern_shape="",
            line=dict(
                # color="grey",
                color="black",
            ),
        ),
        legendgroup="Hypothetical",
        name="Hypothetical",
        showlegend=True,
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

width = 500
height = 700

fig.update_layout(
    # title=f"Shannon's entropy of a largest solution of each {condition_col.lower()}",
    title="Squid's PacBio",
    title_x=0.2,
    width=width,
    height=height,
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
    "Observed vs hypothetical entropy - PacBio.svg",
    width=width,
    height=height,
)

fig.show()

# %%

# %%

# %% [markdown]
# ##### Combinatorics

# %%
max_sol_exp_dfs_new[1]

# %%
top_100_combinatorics_df = (
    max_sol_exp_dfs_new[1]
    .sort_values("%RelativeExpression", ascending=False)
    .reset_index(drop=True)
    .iloc[:100, ML_INPUT_FIRST_COL_POS_NEW:]
    .fillna(-1)
    # .T
)
top_100_combinatorics_df

# %%
data = top_100_combinatorics_df.replace(0, 0.5).replace(-1, 0).values

fig = go.Figure(
    go.Heatmap(
        z=data,
        y=top_100_combinatorics_df.index + 1,
        # x=
        # colorscale=[
        #     [0, "rgb(192,192,192)"],
        #     [0.5, "black"],
        #     [1.0, "rgb(74,246,38)"]
        # ],
        # colorscale=[
        #     [0, "rgb(192,192,192)"], [0.33, "rgb(192,192,192)"],
        #     [0.33, "black"], [0.66, "black"],
        #     [0.66, "rgb(74,246,38)"], [1.0, "rgb(74,246,38)"]
        # ],
        colorscale=[
            [0, "rgb(192,192,192)"],
            [0.3, "rgb(192,192,192)"],
            [0.375, "black"],
            [0.625, "black"],
            [0.7, "rgb(74,246,38)"],
            [1.0, "rgb(74,246,38)"],
        ],
        colorbar=dict(
            title="Non-syn change?",
            # tick0=0,
            # dtick=1,
            tickmode="array",
            tickvals=[0.15, 0.5, 0.85],
            # tickvals=[1, 2, 3],
            ticktext=["Missing", "No", "Yes"],
            len=0.3,
        ),
    )
)
fig.update_xaxes(title="Editable amino acids", side="top", showticklabels=False)
fig.update_yaxes(
    title="Top 100 expressed distinct proteins<br>in squid's PCLO (PacBio)",
    autorange="reversed",
    tick0=1,
    range=[1, 100],
)
fig.update_layout(height=700, width=650, template=template, font_size=16)

fig.write_image(
    "Combinatorics of top 100 expressed proteins in PCLO - PacBio.svg",
    width=650,
    height=700,
)

fig.show()

# %%
unique_proteins_dfs[1]

# %%
unique_proteins_df.columns[unique_proteins_first_col_pos:]

# %%
max_sol_exp_dfs_new[1]

# %% jupyter={"source_hidden": true}
# max_sol_df = max_sol_dfs[1]
# unique_proteins_df = unique_proteins_dfs[1]
# sorting_col = "%WeightedTotalExpression"

# cols_to_use_from_max_sol_df = [
#     condition_col,
#     "Protein",
#     "#Solution",
#     "Algorithm",
#     "TotalEqualSupportingReads",
#     "TotalWeightedSupportingReads",
#     # "Diff5+",
# ]
# cols_to_use_from_unique_proteins_df = [
#     condition_col,
#     "Protein",
# ] + unique_proteins_df.columns[unique_proteins_first_col_pos:].to_list()

# df = max_sol_df.loc[:, cols_to_use_from_max_sol_df].merge(
#     unique_proteins_df.loc[:, cols_to_use_from_unique_proteins_df],
#     how="left",
#     on=[condition_col, "Protein"],
# )

# df

# %%
max_sol_df = max_sol_exp_dfs_new[1]
unique_proteins_df = unique_proteins_dfs[1]
sorting_col = "%RelativeExpression"

cols_to_use_from_max_sol_df = [
    condition_col,
    "Protein",
    "#Solution",
    "Algorithm",
    "%RelativeExpression",
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

df

# %%
# aa_cols

# %%
aa_cols = df.columns[5:]
original_aas = [col.split("(")[1][0] for col in aa_cols]
df_b = df.iloc[:100, 5:]

df_b_value_counts = df_b.apply(
    lambda x: x.value_counts(
        # normalize=True
    ),
    # axis=1
)
# df_b_value_counts = df_b_value_counts.fillna(0).T
# df_b_value_counts = df_b_value_counts.T
df_b_value_counts

# %%
df_b_value_counts_2 = (
    df_b_value_counts.reset_index().rename(columns={"index": "NewAA"}).fillna(0)
)

df_b_value_counts_2 = df_b_value_counts_2.melt(
    id_vars="NewAA", var_name="Position", value_name="Count"
)

df_b_value_counts_2["OriginalAA"] = df_b_value_counts_2["Position"].apply(
    lambda x: x.split("(")[1][0]
)
df_b_value_counts_2["IsRecoding"] = df_b_value_counts_2.apply(
    lambda x: x["OriginalAA"] not in set(x["NewAA"].split(",")), axis=1
)

# df_b_value_counts_2 = df_b_value_counts_2.loc[
#     (df_b_value_counts_2["Count"] > 0) &
#     (df_b_value_counts_2["IsRecoding"])
# ].reset_index(drop=True)

df_b_value_counts_2 = df_b_value_counts_2.loc[df_b_value_counts_2["IsRecoding"]]

df_b_value_counts_2

# %%
df_b_value_counts_2.loc[df_b_value_counts_2["Count"] > 0].groupby(
    "Position"
).size().reset_index()[0].value_counts()

# %%
recoding_changes = df_b_value_counts_2.loc[
    df_b_value_counts_2["Count"] > 0, "NewAA"
].unique()
all_changes = df_b_value_counts_2["NewAA"].unique()

# recoding_color_sequence = px.colors.qualitative.Plotly
recoding_color_sequence = px.colors.qualitative.Set3
recoding_color_sequence = px.colors.qualitative.Light24
recoding_color_discrete_map = {}
i = 0
for aa_change in all_changes:
    if aa_change in recoding_changes:
        color = recoding_color_sequence[i]
        i += 1
    else:
        color = "white"
    recoding_color_discrete_map[aa_change] = color

fig = px.bar(
    df_b_value_counts_2,
    x="Position",
    y="Count",
    color="NewAA",
    color_discrete_map=recoding_color_discrete_map,
    labels={"NewAA": "Recoding<br>substitution"},
    template=template,
    # text_auto=True
)

# Set the "for" loop and "if" condition to hide the legend for the first and third columns.
for trace in fig["data"]:
    # Hide legend column names not equal to 'two'
    if trace["name"] not in recoding_changes:
        trace["showlegend"] = False

# fig.for_each_annotation(
#     lambda a: a.update(text=a.text.replace(f"{condition_col}=", ""))
# )
fig.update_layout(
    title="Recoding events per top 100 expressed distinct proteins in squid's PCLO (PacBio)",
    width=1500,
    height=450,
    xaxis=dict(
        tickfont=dict(size=7),
        title="Transcriptomic positions (corresponding editable AA)",
    ),
)

fig.write_image(
    "Recoding events per top 100 expressed proteins in PCLO - PacBio.svg",
    width=1500,
    height=450,
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
expanded_distinct_unique_reads_df.head()


# %%
cols = len(conditions)

x_title = "% editing in read (mean)"
y_title = "Distinct unique<br>edited reads"
title_text = "Distribution of % editing per Distinct unique edited read"

fig = make_subplots(
    rows=1,
    cols=cols,
    subplot_titles=conditions,
    shared_yaxes=True,
    # shared_xaxes=True,
    x_title=x_title,
    y_title=y_title,
)

x_col = "EditingFrequency"
max_x = 0

for col, condition in zip(range(1, cols + 1), conditions):
    df = expanded_distinct_unique_reads_df.copy()
    df = (
        df.loc[(df["Fraction"] == 1.0) & (df[condition_col] == condition)]
        .groupby([condition_col, "UniqueRead"])[x_col]
        .mean()
        .reset_index()
    )

    x = df[x_col] * 100

    fig.add_trace(
        go.Histogram(
            x=x,
            marker_color=color_discrete_map[condition],
        ),
        row=1,
        col=col,
    )

    max_x = max(max_x, x.max())

fig.update_layout(title_text=title_text, template=template, showlegend=False)

min_x = 0
max_x = max_x * 1.05
fig.update_xaxes(range=[min_x, max_x])

fig.show()


# %% [markdown]
# ## Distribution of non-syns

# %% jupyter={"source_hidden": true}
# cols = min(facet_col_wrap, len(conditions), 4)
# rows = ceil(len(conditions) / cols)
# row_col_iter = list(product(range(1, rows + 1), range(1, cols + 1)))[: len(conditions)]

# x_title = "Non-syn mutations per protein"
# y_title = "Proteins"
# title_text = "Distribution of min & max estimates of non-syn mutations per protein"

# fig = make_subplots(
#     rows=rows,
#     cols=cols,
#     subplot_titles=conditions,
#     shared_yaxes=True,
#     x_title=x_title,
#     y_title=y_title,
# )

# min_x = None
# max_x = 0
# max_y = 0

# col_names = ["MinNonSyns", "MaxNonSyns"]
# estimate_names = ["Min", "Max"]

# for (row, col), condition, proteins_df in zip(row_col_iter, conditions, proteins_dfs):

#     for i, (col_name, estimate_name) in enumerate(zip(col_names, estimate_names)):

#         x = proteins_df[col_name]

#         fig.add_trace(
#             go.Histogram(
#                 x=x,
#                 marker_color=subcolors_discrete_map[condition][i],
#                 name=f"{condition}, {estimate_name}",
#             ),
#             row=row,
#             col=col,
#         )

#         min_x = min(min_x, x.min()) if min_x else x.min()
#         max_x = max(max_x, x.max())
#         # max_y = max(max_y, len(x))

# # add mean lines + text (and also keep track of max_y)

# f = fig.full_figure_for_development(warn=False)
# data_traces = {}
# for condition in conditions:
#     for estimate_name in estimate_names:
#         name = f"{condition}, {estimate_name}"
#         for data in f.data:
#             if data.name == name:
#                 data_traces[name] = data
#                 continue
# for (row, col), condition in zip(row_col_iter, conditions):
#     for estimate_name in estimate_names:
#         name = f"{condition}, {estimate_name}"
#         data = data_traces[name]
#         x = data.x
#         xbins = f.data[0].xbins
#         plotbins = list(
#             np.arange(
#                 start=xbins["start"],
#                 stop=xbins["end"] + xbins["size"],
#                 step=xbins["size"],
#             )
#         )
#         counts, bins = np.histogram(list(x), bins=plotbins)
#         max_count = max(counts)
#         max_y = max(max_y, max_count)
#         x_mean = np.mean(x)
#         fig.add_trace(
#             go.Scatter(
#                 x=[x_mean, x_mean, x_mean],
#                 y=[0, max_count, max_count * 1.1],
#                 mode="lines+text",
#                 line=dict(
#                     color="white",
#                     dash="dash",
#                     width=4,
#                 ),
#                 text=["", "", f"{x_mean:.0f}"],
#                 textposition="top center",
#                 textfont=dict(size=10),
#             ),
#             row=row,
#             col=col,
#         )

# ic(max_y)

# # add legends

# for (row, col), condition in zip(row_col_iter, conditions):
#     for i, (col_name, estimate_name) in enumerate(zip(col_names, estimate_names)):
#         fig.add_trace(
#             go.Scatter(
#                 x=[0.75 * max_x],
#                 # y=[(0.13 * max_y) - (1_000 * i)],
#                 y=[(0.7 * max_y) - (1_000 * i)],
#                 mode="markers+text",
#                 marker=dict(
#                     color=subcolors_discrete_map[condition][i],
#                     size=9,
#                     # opacity=0.7,
#                     symbol="square",
#                     # line=dict(width=0),
#                 ),
#                 text=estimate_name,
#                 textposition="middle right",
#                 textfont=dict(size=9),
#             ),
#             row=row,
#             col=col,
#         )

# fig.update_layout(
#     template=template,
#     barmode="overlay",  # Overlay both histograms
#     title_text=title_text,
#     title_y=0.95,
#     showlegend=False,
#     height=max(200 * rows, 300),
#     width=max(350 * cols, 800),
# )

# fig.update_traces(opacity=0.75)  # Reduce opacity to see both histograms
# fig.update_xaxes(range=[min_x * 0.9, max_x * 1.1])
# fig.update_yaxes(range=[0, max_y * 1.2])

# fig.show()


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
                y=[(0.7 * max_y) - (1_000 * i)],
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

fig.update_traces(opacity=0.75)  # Reduce opacity to see both histograms
fig.update_xaxes(range=[min_x * 0.9, max_x * 1.1])
fig.update_yaxes(range=[0, max_y * 1.2])

fig.update_layout(
    template=template,
    barmode="overlay",  # Overlay both histograms
    title_text="Squid's PacBio",
    title_x=0.1,
    # title_y=0.95,
    showlegend=False,
    height=max(200 * rows, 300),
    width=max(350 * cols, 800),
)

fig.write_image(
    f"{title_text} - PacBio.svg",
    width=max(350 * cols, 800),
    height=max(200 * rows, 300),
)

fig.show()

