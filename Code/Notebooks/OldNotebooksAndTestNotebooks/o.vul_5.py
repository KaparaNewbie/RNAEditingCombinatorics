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
from functools import reduce
from itertools import chain, combinations, product
from math import ceil
from multiprocessing import Pool
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.colors as pc
import plotly.express as px
import plotly.graph_objects as go

from scipy import interpolate  # todo unimport this later?
import scipy.stats
import seaborn as sns
from Bio import SeqIO, motifs  # biopython
from icecream import ic
from logomaker import Logo  # logomaker
from matplotlib_venn import venn2, venn3
from plotly.subplots import make_subplots
from pybedtools import BedTool
from sklearn import linear_model
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.metrics import mean_squared_error, r2_score

sys.path.append(str(Path(code_dir).absolute()))
from Alignment.alignment_utils import (
    count_reads,
    count_reads_in_fastq,
    count_reads_in_unaligned_bam,
    count_unique_filtered_aligned_reads,
)
from EditingUtils.logo import multiple_logos_from_fasta_files

# %%
condition_col = "Transcript"

orfs_bed = "/private7/projects/Combinatorics/O.vulgaris/Annotations/orfs_oct.bed"
alignment_stats_file = "/private7/projects/Combinatorics/O.vulgaris/Alignment/PRJNA791920/IsoSeq/AggregatedByChromBySampleSummary.tsv"
known_sites_file = (
    "/private7/projects/Combinatorics/O.vulgaris/Annotations/O.vul.EditingSites.csv"
)
transcriptome_file = (
    "/private7/projects/Combinatorics/O.vulgaris/Annotations/orfs_oct.fa"
)

positions_dir = Path(
    "/private7/projects/Combinatorics/O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered/PositionsFiles"
)
reads_dir = Path(
    "/private7/projects/Combinatorics/O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered/ReadsFiles/"
)
proteins_dir = Path(
    "/private7/projects/Combinatorics/O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered/ProteinsFiles/"
)

distinct_proteins_dir = Path(
    "/private7/projects/Combinatorics/O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered/DistinctProteins/"
)

neural_vs_non_neural_expression_file = Path(
    "/private7/projects/Combinatorics/O.vulgaris/Annotations/NeuralVsNonNeuralExpression.csv"
)

samples_and_tissues_file = Path(
    "/private7/projects/Combinatorics/O.vulgaris/Data/PRJNA791920/IsoSeqPolished/samples.csv"
)

reads_first_col_pos = 7
unique_reads_first_col_pos = 9
proteins_first_col_pos = 13
unique_proteins_first_col_pos = 15
reads_editing_col = "EditingFrequency"
proteins_editing_col = "MinNonSyns"

reads_type = (
    "CCS"  # something like CCS / miseq / etc. # todo fix to exact pacbio read type
)

samtools_path = "/home/alu/kobish/anaconda3/envs/combinatorics/bin/samtools"
threads = 20
seed = 1892
sep = "\t"

# %%
samples_and_tissues_df = pd.read_csv(samples_and_tissues_file)
samples_and_tissues_df["Tissue"] = samples_and_tissues_df["Tissue"].str.capitalize()
samples_and_tissues_df

# %%
samples = samples_and_tissues_df["Sample"]
tissues = samples_and_tissues_df["Tissue"]
sample_to_tissue_dict = {sample: tissue for sample, tissue in zip(samples, tissues)}
sample_to_tissue_dict

# %%
orfs_df = pd.read_csv(
    orfs_bed, sep="\t", names="Chrom Start End Name Score Strand".split()
)
orfs_df

# %%
alignment_stats_df = pd.read_csv(alignment_stats_file, sep="\t")

# alignment_stats_df["UsedForPileup"] = alignment_stats_df.apply(
#     lambda x: (x["Samples"] >= min_samples)
#     and (x["MappedReadsPerSample"] >= min_mapped_reads_per_sample)
#     and (x["KnownSites"] >= min_known_sites),
#     axis=1,
# )

alignment_stats_df

# %%
positions_files = list(positions_dir.glob("*.positions.csv.gz"))

chroms_in_positions = [
    positions_file.name.split(".")[0] for positions_file in positions_files
]

positions_data_df = pd.DataFrame(
    {
        "Chrom": chroms_in_positions,
        "PositionsFile": positions_files,
    }
)

positions_data_df

# %%
reads_files = list(reads_dir.glob("*.reads.csv.gz"))
chroms_in_reads_files = [reads_file.name.split(".")[0] for reads_file in reads_files]

unique_reads_files = list(reads_dir.glob("*.unique_reads.csv.gz"))
chroms_in_unique_reads_files = [
    unique_reads_file.name.split(".")[0] for unique_reads_file in unique_reads_files
]

reads_data_df = pd.DataFrame(
    {
        "Chrom": chroms_in_reads_files,
        "ReadsFile": reads_files,
    }
)

unique_reads_data_df = pd.DataFrame(
    {
        "Chrom": chroms_in_unique_reads_files,
        "UniqueReadsFile": unique_reads_files,
    }
)

reads_data_df = reads_data_df.merge(unique_reads_data_df, on="Chrom", how="left")

reads_data_df

# %%
proteins_files = list(proteins_dir.glob("*.proteins.csv.gz"))
chroms_in_proteins_files = [
    proteins_file.name.split(".")[0] for proteins_file in proteins_files
]

unique_proteins_files = list(proteins_dir.glob("*.unique_proteins.csv.gz"))
chroms_in_unique_proteins_files = [
    unique_proteins_file.name.split(".")[0]
    for unique_proteins_file in unique_proteins_files
]

# distinct_proteins_files = list(distinct_proteins_dir.glob("*DistinctUniqueProteins.*.csv"))
distinct_proteins_files = [
    f
    for f in distinct_proteins_dir.glob("*DistinctUniqueProteins.*.csv")
    if "expression" not in f.name.lower()
]
chroms_in_distinct_proteins_files = [
    distinct_proteins_file.name.split(".")[0]
    for distinct_proteins_file in distinct_proteins_files
]

expression_files = list(
    distinct_proteins_dir.glob("*.DistinctUniqueProteins.ExpressionLevels.csv")
)
chroms_in_expression_files = [
    expression_file.name.split(".")[0] for expression_file in expression_files
]


proteins_data_df = pd.DataFrame(
    {
        "Chrom": chroms_in_proteins_files,
        "ProteinsFile": proteins_files,
    }
)

unique_proteins_data_df = pd.DataFrame(
    {
        "Chrom": chroms_in_unique_proteins_files,
        "UniqueProteinsFile": unique_proteins_files,
    }
)

distinct_proteins_data_df = pd.DataFrame(
    {
        "Chrom": chroms_in_distinct_proteins_files,
        "DistinctProteinsFile": distinct_proteins_files,
    }
)

expression_data_df = pd.DataFrame(
    {"Chrom": chroms_in_expression_files, "ExpressionFile": expression_files}
)

proteins_data_df = (
    proteins_data_df.merge(unique_proteins_data_df, on="Chrom", how="left")
    .merge(distinct_proteins_data_df, on="Chrom", how="left")
    .merge(expression_data_df, on="Chrom", how="left")
)

proteins_data_df

# %%
data_df = (
    orfs_df.merge(alignment_stats_df, on="Chrom", how="left")
    .merge(positions_data_df, on="Chrom", how="left")
    .merge(reads_data_df, on="Chrom", how="left")
    .merge(proteins_data_df, on="Chrom", how="left")
)

data_df

# %%
complete_data_df = data_df.loc[data_df["ExpressionFile"].notna()].reset_index(drop=True)
# complete_data_df = complete_data_df.drop_duplicates(
#     "Name", keep=False, ignore_index=True
# )
complete_data_df

# %%
complete_data_df["Strand"].value_counts()

# %%
# complete_data_df.loc[complete_data_df["Name"].duplicated(keep=False)].sort_values("Name", ignore_index=True)

# %%
conditions = complete_data_df["Name"].tolist()
chroms = complete_data_df["Chrom"].tolist()
starts = complete_data_df["Start"].tolist()
ends = complete_data_df["End"].tolist()
strands = complete_data_df["Strand"].tolist()

positions_files = complete_data_df["PositionsFile"].tolist()
reads_files = complete_data_df["ReadsFile"].tolist()
unique_reads_files = complete_data_df["UniqueReadsFile"].tolist()
proteins_files = complete_data_df["ProteinsFile"].tolist()
unique_proteins_files = complete_data_df["UniqueProteinsFile"].tolist()
distinct_unique_proteins_files = complete_data_df["DistinctProteinsFile"].tolist()
expression_files = complete_data_df["ExpressionFile"].tolist()


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
sample_to_tissue_dict

# %% papermill={"duration": 0.054755, "end_time": "2022-02-01T09:42:46.304499", "exception": false, "start_time": "2022-02-01T09:42:46.249744", "status": "completed"}
# # plotly consts
# color_sequence = px.colors.qualitative.Pastel
# # color_sequence = px.colors.qualitative.D3
color_sequence = px.colors.qualitative.G10

samples_color_discrete_map = {
    sample: color for sample, color in zip(samples, color_sequence)
}
samples_subcolors_discrete_map = {
    sample: two_subcolors_from_hex(samples_color_discrete_map[sample])
    for sample in samples
}

tissues_color_discrete_map = {
    tissue: color for tissue, color in zip(tissues, color_sequence)
}
tissues_subcolors_discrete_map = {
    tissue: two_subcolors_from_hex(tissues_color_discrete_map[tissue])
    for tissue in tissues
}

# ic(color_discrete_map)
# ic(subcolors_discrete_map)
# category_orders = {condition_col: conditions}
# horizontal_category_orders = {
#     category: list(reversed(category_orders[category])) for category in category_orders
# }
# # valid_shapes = ['', '/', '\\', 'x', '-', '|', '+', '.']
# # pattern_shape_map = {
# #     condition: shape for condition, shape in zip(conditions, cycle(valid_shapes))
# # }
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
# def n_repetitions_colormap(subcolors_discrete_map, condition, n_repetitions):
#     lowcolor, highcolor = subcolors_discrete_map[condition]
#     colors = pc.n_colors(lowcolor, highcolor, n_repetitions, colortype="rgb")
#     return {i: color for i, color in enumerate(colors, start=1)}

# %%
# n_repetitions_colormap(subcolors_discrete_map, "GRIA", 10)

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
# 25_000 / (750 * 1.17)

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
positions_df = pd.concat(positions_dfs, ignore_index=True)
positions_df

# %%

# %%

# %%
editing_positions_per_sample = [len(df.loc[(df["Edited"])]) for df in positions_dfs]
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
# {
#     len(reads_df["Tissue"].unique()) for reads_df in reads_dfs
# }

# %%
# edited_reads_dfs = [
#     reads_df.loc[reads_df[reads_editing_col] > 0] for reads_df in reads_dfs
# ]
# edited_reads_dfs[0]


# %% [markdown] papermill={"duration": 0.041741, "end_time": "2022-02-01T09:42:47.760215", "exception": false, "start_time": "2022-02-01T09:42:47.718474", "status": "completed"}
# ### Unique

# %% papermill={"duration": 0.126539, "end_time": "2022-02-01T09:42:47.923363", "exception": false, "start_time": "2022-02-01T09:42:47.796824", "status": "completed"}
unique_reads_dfs = [
    pd.read_csv(unique_reads_file, sep=sep) for unique_reads_file in unique_reads_files
]
for chrom, unique_reads_df in zip(chroms, unique_reads_dfs):
    unique_reads_df.insert(0, "Chrom", chrom)
unique_reads_dfs[0]


# %%
unique_reads_dfs[0].columns

# %%
expanded_unique_reads_dfs = []
for unique_reads_df in unique_reads_dfs:
    unique_reads_df = unique_reads_df.copy()
    unique_reads_df["Samples"] = unique_reads_df["Samples"].str.split(",")
    expanded_unique_reads_df = unique_reads_df.explode("Samples").reset_index(drop=True)
    # after exploding the df by samples, it may be that some reads/unique reads only appear in certain samples,
    # so in order to get that information, one would have to merge the `expanded_unique_proteins_df`
    # with a corresponding `expanded_reads_df`/`expanded_unique_reads_df`
    # expanded_unique_reads_df = expanded_unique_reads_df.drop(
    #     [
    #         "Reads",
    #         "NumOfReads",
    #         "UniqueReads",
    #         "NumOfUniqueReads",
    #         "EditingFrequency",
    #         "EditedPositions",
    #         "UneditedPositions",
    #         "AmbigousPositions",
    #     ],
    #     axis=1,
    # )
    # expanded_unique_reads_df = expanded_unique_reads_df.rename(
    #     columns={"Samples": "Sample"}
    # )
    expanded_unique_reads_dfs.append(expanded_unique_reads_df)
    break

expanded_unique_reads_dfs[0]

# expanded_unique_proteins_df = pd.concat(expanded_unique_proteins_dfs, ignore_index=True)
# del expanded_unique_proteins_dfs
# expanded_unique_proteins_df

# %%
expanded_unique_reads_dfs[0].columns

# %% [markdown]
# ## Proteins

# %% [markdown]
# ### All proteins

# %%
proteins_dfs = [pd.read_csv(proteins_file, sep=sep) for proteins_file in proteins_files]
# for proteins_df in proteins_dfs:
#     if "Transcript" in proteins_df.columns:
#         proteins_df.rename(columns={"Transcript": "UniqueRead"}, inplace=True)
proteins_dfs[0]


# %% [markdown]
# ### Unique proteins

# %%
unique_proteins_dfs = [
    pd.read_csv(unique_proteins_file, sep=sep, dtype={"Protein": str})
    for unique_proteins_file in unique_proteins_files
]
for chrom, unique_proteins_df in zip(chroms, unique_proteins_dfs):
    unique_proteins_df.insert(0, "Chrom", chrom)
# for unique_proteins_df in unique_proteins_dfs:
#     unique_proteins_df.rename(
#         columns={
#             col: col.replace("Transcripts", "UniqueReads")
#             for col in unique_proteins_df.columns[:unique_proteins_first_col_pos]
#             if "Transcripts" in col
#         },
#         inplace=True,
#     )
unique_proteins_dfs[0]


# %%
expanded_unique_proteins_dfs = []
for unique_proteins_df in unique_proteins_dfs:
    unique_proteins_df = unique_proteins_df.copy()
    unique_proteins_df["Samples"] = unique_proteins_df["Samples"].str.split(",")
    expanded_unique_proteins_df = unique_proteins_df.explode("Samples").reset_index(
        drop=True
    )
    # after exploding the df by samples, it may be that some reads/unique reads only appear in certain samples,
    # so in order to get that information, one would have to merge the `expanded_unique_proteins_df`
    # with a corresponding `expanded_reads_df`/`expanded_unique_reads_df`
    expanded_unique_proteins_df = expanded_unique_proteins_df.drop(
        [
            "Reads",
            "NumOfReads",
            "UniqueReads",
            "NumOfUniqueReads",
            "EditingFrequency",
            "EditedPositions",
            "UneditedPositions",
            "AmbigousPositions",
        ],
        axis=1,
    )
    expanded_unique_proteins_df = expanded_unique_proteins_df.rename(
        columns={"Samples": "Sample"}
    )
    expanded_unique_proteins_dfs.append(expanded_unique_proteins_df)
    # break

expanded_unique_proteins_dfs[0]

# expanded_unique_proteins_df = pd.concat(expanded_unique_proteins_dfs, ignore_index=True)
# del expanded_unique_proteins_dfs
# expanded_unique_proteins_df

# %%
editable_aas_per_sample = [
    df.iloc[:, unique_proteins_first_col_pos:].shape[1] for df in unique_proteins_dfs
]

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

# %% [markdown]
# ### Distinct unique proteins

# %%
distinct_unique_proteins_dfs = []
for condition, chrom, distinct_unique_proteins_file, unique_reads_df in zip(
    conditions, chroms, distinct_unique_proteins_files, unique_reads_dfs
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
    distinct_unique_proteins_df.insert(0, "Chrom", chrom)
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
# complete_data_df.loc[
#     (complete_data_df["Name"].isin(realizations_count_df.loc[realizations_count_df["Count"] < 16, condition_col]), ["Chrom", "Name"])
# ].values

# %%
# num_of_reads_per_transcript_and_fraction_df = (
#     distinct_unique_proteins_df.groupby([condition_col, "Fraction"])["NumOfReads"]
#     .unique()
#     .reset_index()
# )
# # num_of_reads_per_transcript_and_fraction_df = num_of_reads_per_transcript_and_fraction_df.explode("NumOfReads", ignore_index=True)
# num_of_reads_per_transcript_and_fraction_df

# %%
# num_of_reads_per_transcript_and_fraction_df["NumOfReads"].apply(len).value_counts()

# %%
# expanded_distinct_unique_proteins_df = (
#     distinct_unique_proteins_df.copy()
#     .assign(Proteins2=lambda x: x.Proteins.str.split(","))
#     .drop("Proteins", axis=1)
#     .rename(columns={"Proteins2": "Proteins"})
#     .explode("Proteins")
#     .rename(columns={"Proteins": "Protein", "NumOfReads": "NumOfReadsInFraction"})
#     .drop(["NumOfProteins"], axis=1)
#     .merge(
#         pd.concat(
#             [df.iloc[:, :unique_proteins_first_col_pos] for df in unique_proteins_dfs]
#         ),
#         on=[condition_col, "Protein"],
#     )
# )

# expanded_distinct_unique_proteins_df


# %%
# distinct_unique_proteins_df2 = (
#     expanded_distinct_unique_proteins_df.groupby(
#         [
#             condition_col,
#             "Fraction",
#             "FractionRepetition",
#             "Algorithm",
#             "AlgorithmRepetition",
#         ]
#     )["NumOfReads"]
#     .sum()
#     .reset_index()
#     .rename(columns={"NumOfReads": "NumOfSupportingReads"})
#     .merge(
#         distinct_unique_proteins_df,
#         on=[
#             condition_col,
#             "Fraction",
#             "FractionRepetition",
#             "Algorithm",
#             "AlgorithmRepetition",
#         ],
#     )
#     .assign(
#         SupportingReadsPerProtein=lambda x: x["NumOfSupportingReads"]
#         / x["NumOfProteins"],
#         PercentSupportedReads=lambda x: 100
#         * x["NumOfSupportingReads"]
#         / x["NumOfReads"],
#     )
#     .rename(columns={"PercentSupportedReads": "%SupportedReads"})
# )
# distinct_unique_proteins_df2

# %% [markdown]
# ## Summary of data loss

# %%
alignment_stats_df = (
    pd.read_csv(alignment_stats_file, sep=sep)
    .merge(orfs_df.loc[:, ["Chrom", "Name"]])
    .rename(columns={"Name": condition_col})
)
alignment_stats_df[condition_col] = alignment_stats_df[condition_col].astype(str)
alignment_stats_df

# %%

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
    .mean()
    .round()
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
        "Distinct unique proteins (mean)": distinct_unique_proteins_counts,
    },
    index=conditions,
)
data_loss_df


# %% [markdown] papermill={"duration": 0.045853, "end_time": "2022-02-01T09:42:48.953594", "exception": false, "start_time": "2022-02-01T09:42:48.907741", "status": "completed"}
# # Results

# %% [markdown] papermill={"duration": 0.149848, "end_time": "2022-02-01T09:43:12.800733", "exception": false, "start_time": "2022-02-01T09:43:12.650885", "status": "completed"}
# ## Data loss

# %%
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
# fig = px.bar(
#     long_data_loss_df.loc[
#         long_data_loss_df["Stage"] != "Distinct unique reads<br>(mean)"
#     ].sort_values("Count"),
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
#     width=800,
#     height=550,
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


# %% [markdown] papermill={"duration": 0.124528, "end_time": "2022-02-01T09:43:10.054394", "exception": false, "start_time": "2022-02-01T09:43:09.929866", "status": "completed"}
# ## Positions

# %% [markdown] papermill={"duration": 0.149848, "end_time": "2022-02-01T09:43:12.800733", "exception": false, "start_time": "2022-02-01T09:43:12.650885", "status": "completed"}
# ### Coverage - per-transcript per-sample

# %%
positions_dfs[0]


# %%
def calc_per_transcript_per_sample_coverage(positions_df, 
                                            samples_and_tissues_df
                                            # samples
                                           ):
    expanded_positions_df = (
        positions_df.loc[~positions_df["InProbRegion"]]
        .reset_index(drop=True)
        .drop(
            [
                "Phred",
                "MappedBases",
                "Noise",
                "EditingFrequency",
                "A",
                "T",
                "C",
                "G",
                "TotalCoverage",
            ],
            axis=1,
        )
    )

    expanded_positions_df["Samples"] = expanded_positions_df["Samples"].str.split(",")
    expanded_positions_df["Reads"] = expanded_positions_df["Reads"].str.split(",")

    # now is the time the df is really expanded
    expanded_positions_df = expanded_positions_df.explode(["Samples", "Reads"])
    expanded_positions_df = expanded_positions_df.rename(
        columns={"Samples": "Sample", "Reads": "Read"}
    )

    per_sample_per_transcript_coverage_df = (
        expanded_positions_df.groupby(["Chrom", "Transcript", "Sample"])["Read"]
        .apply(lambda x: x.unique().size)
        .reset_index()
        .rename(columns={"Read": "NumOfReads"})
        .merge(samples_and_tissues_df, how="left")
    )
    
    return per_sample_per_transcript_coverage_df


# %%
# expanded_positions_df = (
#     positions_dfs[0]
#     .loc[~positions_dfs[0]["InProbRegion"]]
#     .reset_index(drop=True)
#     .drop(
#         [
#             "Phred",
#             "MappedBases",
#             "Noise",
#             "EditingFrequency",
#             "A",
#             "T",
#             "C",
#             "G",
#             "TotalCoverage",
#         ],
#         axis=1,
#     )
# )

# expanded_positions_df["Samples"] = expanded_positions_df["Samples"].str.split(",")
# expanded_positions_df["Reads"] = expanded_positions_df["Reads"].str.split(",")

# # now is the time the df is really expanded
# expanded_positions_df = expanded_positions_df.explode(["Samples", "Reads"])
# expanded_positions_df = expanded_positions_df.rename(
#     columns={"Samples": "Sample", "Reads": "Read"}
# )

# expanded_positions_df

# %%
# per_sample_coverage = expanded_positions_df.groupby("Sample")["Read"].apply(lambda x: x.unique().size)
# per_sample_coverage

# complete_per_sample_coverage = []
# for sample in samples:
#     if sample in per_sample_coverage:
#         coverage = per_sample_coverage[sample]
#     else:
#         coverage = 0
#     complete_per_sample_coverage.append(coverage)

# per_sample_per_transcript_coverage_df = pd.DataFrame({"Sample": samples, "Reads": complete_per_sample_coverage})
# # per_sample_per_transcript_coverage_df.insert(0, )
# per_sample_per_transcript_coverage_df

# %%
# per_sample_per_transcript_coverage_df = (
#     expanded_positions_df.groupby(["Chrom", "Transcript", "Sample"])["Read"]
#     .apply(lambda x: x.unique().size)
#     .reset_index()
#     .rename(columns={"Read": "NumOfReads"})
# )
# per_sample_per_transcript_coverage_df

# %%
calc_per_transcript_per_sample_coverage(positions_dfs[0])

# %%
# per_transcript_per_sample_coverage_dfs = [
#     calc_per_transcript_per_sample_coverage(positions_df)
#     for positions_df in positions_dfs
# ]

with Pool(processes=10) as pool:
    per_transcript_per_sample_coverage_dfs = pool.starmap(
        func=calc_per_transcript_per_sample_coverage, 
        iterable=[(positions_df, samples_and_tissues_df) for positions_df in positions_dfs]
    )
per_transcript_per_sample_coverage_dfs[0]

# %%
merged_per_transcript_per_sample_coverage_df = pd.concat(per_transcript_per_sample_coverage_dfs).reset_index()
# merged_per_transcript_per_sample_coverage_df = merged_per_transcript_per_sample_coverage_df.merge(samples_and_tissues_df, how="left")
merged_per_transcript_per_sample_coverage_df

# %%
fig = px.histogram(
    merged_per_transcript_per_sample_coverage_df,
    x="NumOfReads",
    color="Tissue",
    color_discrete_map=tissues_color_discrete_map,
    facet_col="Tissue",
    facet_col_wrap=3,
    facet_col_spacing=facet_col_spacing,
    facet_row_spacing=facet_row_spacing*0.5,
    log_y=True,
)

width = 900
height = 600
# width = 800
# height = 300

# # Overlay both histograms
# fig.update_layout(barmode='overlay')
# # Reduce opacity to see both histograms
# fig.update_traces(opacity=0.5)

fig.update_layout(
    template=template,
    showlegend=False,
    width=width,
    height=height
)
fig.show()


# %%

# %%

# %%

# %%

# %% [markdown] papermill={"duration": 0.149848, "end_time": "2022-02-01T09:43:12.800733", "exception": false, "start_time": "2022-02-01T09:43:12.650885", "status": "completed"}
# ### Editing index - per transcript

# %%
def editing_index_per_transcript(positions_df, strand):
    ref_base = "A" if strand == "+" else "T"
    alt_base = "G" if strand == "+" else "C"
    all_refbase_positions_df = positions_df.loc[
        (positions_df["RefBase"] == ref_base) & (~positions_df["InProbRegion"])
    ]
    num_of_all_editable_adenosines = all_refbase_positions_df["TotalCoverage"].sum()
    num_of_edited_adenosines = all_refbase_positions_df[alt_base].sum()
    editing_index = 100 * num_of_edited_adenosines / num_of_all_editable_adenosines
    return editing_index


# %%
per_transcript_editing_indices = [
    editing_index_per_transcript(positions_df, strand)
    for positions_df, strand in zip(positions_dfs, strands)
]
per_transcript_editing_indices[0]

# %%
per_transcript_editing_index_df = pd.DataFrame(
    {"Chrom": chroms, "EditingIndex": per_transcript_editing_indices}
)
per_transcript_editing_index_df

# %%
fig = px.histogram(
    per_transcript_editing_index_df,
    x="EditingIndex",
    # y="TotalCoverage",
    # color="EditingStatus",
    # log_y=True
    color_discrete_sequence=["black"],
    labels={"EditingIndex": "% editing"},
)
fig.update_layout(
    #  xaxis_title="Editing frequency",
    # title="Octopus",
    title="Pooled octopus data",
    title_x=0.15,
    yaxis_title="Transcripts",
    template=template,
    width=800 * 0.7,
    height=800 * 0.5,
    #  showlegend=False
)

# fig.write_image(
#     "Mean per chrom noise levels - Octopus.svg",
#     # width=800,
#     # height=500,
#     width=width*0.7,
#     height=width*0.5,
# )

fig.show()


# %% [markdown] papermill={"duration": 0.149848, "end_time": "2022-02-01T09:43:12.800733", "exception": false, "start_time": "2022-02-01T09:43:12.650885", "status": "completed"}
# ### Editing index - per sample

# %%
def calc_per_transcript_per_sample_a_and_g_counts(positions_df, strand, samples):

    ref_base = "A" if strand == "+" else "T"
    alt_base = "G" if strand == "+" else "C"

    expanded_all_refbase_positions_df = (
        positions_df.loc[
            (positions_df["RefBase"] == ref_base) & (~positions_df["InProbRegion"])
        ]
        .reset_index(drop=True)
        .drop(
            [
                "Phred",
                "Reads",
                "Noise",
                "EditingFrequency",
                "A",
                "T",
                "C",
                "G",
                "TotalCoverage",
            ],
            axis=1,
        )
    )

    expanded_all_refbase_positions_df["Samples"] = expanded_all_refbase_positions_df[
        "Samples"
    ].str.split(",")
    expanded_all_refbase_positions_df[
        "MappedBases"
    ] = expanded_all_refbase_positions_df["MappedBases"].apply(list)

    # now is the time the df is really expanded
    expanded_all_refbase_positions_df = expanded_all_refbase_positions_df.explode(
        ["Samples", "MappedBases"]
    )
    expanded_all_refbase_positions_df = expanded_all_refbase_positions_df.rename(
        columns={"Samples": "Sample", "MappedBases": "MappedBase"}
    )

    expanded_all_refbase_positions_df["MappedBase"] = expanded_all_refbase_positions_df[
        "MappedBase"
    ].apply(lambda x: ref_base if x == "." else x)

    per_transcript_base_counts = expanded_all_refbase_positions_df.groupby("Sample")[
        "MappedBase"
    ].value_counts()

    per_transcript_per_sample_a_and_g_counts = []
    for sample in samples:
        if sample in per_transcript_base_counts:
            # a_count represents the unedited adenosines (possibly zero)
            try:
                a_count = per_transcript_base_counts[sample][ref_base]
            except KeyError:
                a_count = 0
            # g_count represnts the edited adenosines (possibly zero as well)
            try:
                g_count = per_transcript_base_counts[sample][alt_base]
            except KeyError:
                g_count = 0
        else:
            a_count = 0
            g_count = 0
        per_transcript_per_sample_a_and_g_counts.append((a_count, g_count))

    return per_transcript_per_sample_a_and_g_counts


# %%
def calc_all_per_sample_a_and_g_counts(positions_dfs, strands, samples, processes=1):
    # all_per_transcript_per_sample_a_and_g_counts = [
    #     calc_per_transcript_per_sample_a_and_g_counts(positions_df, strand, samples)
    #     for positions_df, strand in zip(positions_dfs, strands)
    # ]
    with Pool(processes=processes) as pool:
        all_per_transcript_per_sample_a_and_g_counts = pool.starmap(
            func=calc_per_transcript_per_sample_a_and_g_counts,
            iterable=[
                (positions_df, strand, samples)
                for positions_df, strand in zip(positions_dfs, strands)
            ],
        )
    # all_per_transcript_per_sample_a_and_g_counts

    all_per_sample_a_and_g_counts = [[0, 0] for _ in samples]
    for (
        per_transcript_per_sample_a_and_g_counts
    ) in all_per_transcript_per_sample_a_and_g_counts:
        for i, (a_count, g_count) in enumerate(
            per_transcript_per_sample_a_and_g_counts
        ):
            all_per_sample_a_and_g_counts[i][0] += a_count
            all_per_sample_a_and_g_counts[i][1] += g_count

    return all_per_sample_a_and_g_counts


# %%
def calc_per_sample_editing_index(positions_dfs, strands, samples, processes=1):
    all_per_sample_a_and_g_counts = calc_all_per_sample_a_and_g_counts(
        positions_dfs, strands, samples, processes
    )
    per_sample_editing_index = [
        100 * g_count / (a_count + g_count)
        for a_count, g_count in all_per_sample_a_and_g_counts
    ]
    return per_sample_editing_index


# %%
per_sample_editing_index = calc_per_sample_editing_index(
    positions_dfs, strands, samples, processes=10
)
per_sample_editing_index

# %%
per_sample_editing_index_df = pd.DataFrame(
    {"Sample": samples, "EditingIndex": per_sample_editing_index}
)
per_sample_editing_index_df

# %% [markdown]
# ### Noise in positions

# %%
positions_df

# %%
per_chrom_mean_noise_levels = (
    positions_df
    .loc[positions_df["Noise"] <= 0.1]
    .groupby("Chrom")["Noise"].nlargest(3)
    .reset_index().drop("level_1", axis=1)
    .groupby("Chrom").mean().reset_index()
)
per_chrom_mean_noise_levels

# %%
per_chrom_mean_noise_levels["% noise"] = per_chrom_mean_noise_levels["Noise"] * 100
fig = px.histogram(
    per_chrom_mean_noise_levels,
    x="% noise",
    # y="TotalCoverage",
    # color="EditingStatus",
    # log_y=True
    color_discrete_sequence=['black'],
    labels={"% noise": "% mean noise"}
)
width = 560
height = 400

fig.update_layout(
    #  xaxis_title="Editing frequency",
    # title="Octopus",
    title="Pooled octopus data",
    title_x=0.15,
     yaxis_title="Transcripts",
     template=template,
     width=width,
     height=height,
    #  showlegend=False
     
)

fig.write_image(
    "Mean per chrom noise levels - Octopus.svg",
    width=width,
    height=height,
)

fig.show()

# %% [markdown]
# ### Known & new editing sites

# %%
cols = 1
rows = 1

fig, ax = plt.subplots(
    # nrows=rows,
    # ncols=cols,
    figsize=(3.5 * cols, 2.5 * rows),
    constrained_layout=True,
    gridspec_kw=dict(hspace=0.2, wspace=0.3),
)

labels = ["Edited", "KnownEditing"]

sets = [
    set(
        positions_df.loc[
            positions_df[label], [condition_col, "Chrom", "Position"]
        ].itertuples(index=False)
    )
    for label in labels
]

labels[0] = f"Edited\n({len(sets[0])})"
labels[1] = f"Known editing\n({len(sets[1])})"

venn2(sets, set_labels=labels, ax=ax)

fig.suptitle(
    "Pooled octopus data",
    fontsize="xx-large",
    # y=1.2
)
# plt.title("Pooled octopus data", fontsize=16,
#           # y=1.2
#          )
# ax.set_title("Pooled octopus data", fontdict=dict(fontsize=16))
# # ax.set_title("", fontdict=dict(fontsize=14))
# # fig.suptitle("Pooled octopus data", fontsize="xx-large", y=1.2)
# fig.tight_layout()

plt.savefig("Known vs new editing sites - Octopus.svg", format="svg", dpi=300)

plt.show()

# %% [markdown]
# ### ADAR motif

# %%
unique_positions_df = (
    positions_df.loc[positions_df["Edited"], ["Chrom", "Position"]]
    .drop_duplicates(["Chrom", "Position"], ignore_index=True)
    .merge(
        orfs_df.loc[:, ["Chrom", "Strand", "End"]].rename(columns={"End": "ChromEnd"}),
        how="left",
    )
    .rename(columns={"Position": "Start"})
)
unique_positions_df.insert(2, "End", unique_positions_df["Start"] + 1)
unique_positions_df.insert(3, "Name", ".")
unique_positions_df.insert(4, "Score", ".")

# extend start and end positions s.t. each region spans 3 bases, with the edited adenosine in the middle
unique_positions_df["Start"] = (
    unique_positions_df["Start"] - 1
)  # upstream base of edited adenosine (or downstream for a transcript expressed from the negative strand)
unique_positions_df["End"] = (
    unique_positions_df["End"] + 1
)  # downstream base of edited adenosine (or upstream for a transcript expressed from the negative strand)

# don't consider editing sites located at the first/last position of a transcript (if there are such editing sites, they are negligble)
unique_positions_df = unique_positions_df.loc[unique_positions_df["Start"] >= 0]
unique_positions_df = unique_positions_df.loc[
    unique_positions_df.apply(lambda x: x["End"] <= x["ChromEnd"], axis=1)
]
del unique_positions_df["ChromEnd"]

unique_positions_df

# %%
editing_sites_bedtool = (
    BedTool()
    .from_dataframe(unique_positions_df)
    .sort()
    .sequence(fi=transcriptome_file, s=True)
)


fasta_files = [editing_sites_bedtool.seqfn]
main_title = None
sub_titles = ["Pooled octopus sites"]

# %%
out_file = Path("ADAR motif of pooled editing sites - Octopus.svg")

# %%
multiple_logos_from_fasta_files(
    fasta_files,
    main_title,
    sub_titles,
    out_file,
    width=14 / 3,
    height=4,
);

# %%

# %% [markdown] papermill={"duration": 0.030615, "end_time": "2022-02-01T09:42:49.024262", "exception": false, "start_time": "2022-02-01T09:42:48.993647", "status": "completed"}
# ## Num of distinct proteins

# %%
neural_vs_non_neural_expression_df = pd.read_csv(
    "/private7/projects/Combinatorics/O.vulgaris/Annotations/NeuralVsNonNeuralExpression.BySalmonAndOrthoFinder.tsv",
    sep="\t",
)

# # the original file from Y. Shoshan's paper contained a line per editing site,
# # but the per-transcript ("Chrom") expression levels are the same for each transcript,
# # so we remove duplicates s.t. each transcript will appear only oncee
# neural_vs_non_neural_expression_df = neural_vs_non_neural_expression_df.drop_duplicates(
#     subset="Chrom", ignore_index=True
# )

# # determine whether a transcript is highly expressed in neural tissues
# neural_vs_non_neural_expression_df[
#     "IsNeural"
# ] = neural_vs_non_neural_expression_df.apply(
#     lambda x: "Yes" if x["NeuralExpression"] > 4 * x["NonNeuralExpression"] else "No",
#     axis=1,
# )
neural_vs_non_neural_expression_df

# %%
fig = px.histogram(
    neural_vs_non_neural_expression_df,
    x="NeuralObimOrthologs/ObimOrthologs",
    log_y=True,
)
fig.update_layout(width=600, height=400, template=template)
fig.show()

# %%
fig = px.histogram(neural_vs_non_neural_expression_df, x="IsNeural", log_y=True)
fig.update_layout(width=600, height=400, template=template)
fig.show()

# %%
max_distinct_proteins_df = (
    distinct_unique_proteins_df.sort_values("Fraction", ascending=False)
    .groupby("Chrom")
    .apply(pd.DataFrame.nlargest, n=1, columns="NumOfProteins")
)
max_distinct_proteins_df = (
    max_distinct_proteins_df.drop("Chrom", axis=1).reset_index().drop("level_1", axis=1)
)

# # max_distinct_proteins_df[condition_col] = max_distinct_proteins_df[
# #     condition_col
# # ].astype(str)

max_distinct_proteins_df = max_distinct_proteins_df.merge(
    alignment_stats_df,
    on="Chrom",
    # how="left",
    how="right",
)

# max_distinct_proteins_df["NumOfProteins"] = max_distinct_proteins_df[
#     "NumOfProteins"
# ].fillna(1)

# max_distinct_proteins_df["NumOfReads"] = max_distinct_proteins_df.apply(
#     lambda x: x["NumOfReads"] if not pd.isna(x["NumOfReads"]) else x["MappedReads"],
#     axis=1,
# )

max_distinct_proteins_df = max_distinct_proteins_df.dropna().reset_index(drop=True)

max_distinct_proteins_df["DistinctProteins/Reads"] = (
    max_distinct_proteins_df["NumOfProteins"] / max_distinct_proteins_df["NumOfReads"]
)

max_distinct_proteins_df = max_distinct_proteins_df.merge(
    neural_vs_non_neural_expression_df.loc[:, ["OvulChrom", "IsNeural"]].rename(
        columns={"OvulChrom": "Chrom"}
    ),
    on="Chrom",
    how="left",
)

max_distinct_proteins_df["IsNeural"] = max_distinct_proteins_df["IsNeural"].fillna(
    "Missing"
)

max_distinct_proteins_df

# %%
# max_distinct_proteins_df.loc[max_distinct_proteins_df["NumOfProteins"] == 1]

# %%
max_distinct_proteins_df["IsNeural"].value_counts()

# %%
max_distinct_proteins_df["NumOfProteins"].sort_values(ascending=False).head(8)

# %%
fig = px.histogram(max_distinct_proteins_df, x="IsNeural", log_y=True)
fig.update_layout(width=600, height=400, template=template)
fig.show()

# %%
df = (
    max_distinct_proteins_df.loc[:, ["NumOfProteins"]]
    .sort_values("NumOfProteins")
    .reset_index(drop=True)
)
df["CummulativeTranscripts"] = 100 * (df.index + 1) / len(df)
df["CummulativeTranscripts"] = df["CummulativeTranscripts"][::-1].values
df = df.drop_duplicates(subset="NumOfProteins").reset_index(drop=True)
# df

neural_conditions = ["Yes", "No", "Missing"]
neural_trace_names = ["Neural", "Non-neural", "Missing"]
neural_color_discrete_map = {
    "Yes": "red",
    "No": "rgb(0,170,255)",  # kind of azure
    "Missing": "rgb(192,192,192)",  # kind of terminal grey
}
neural_dfs = []
for neural_condition in neural_conditions:
    neural_df = (
        max_distinct_proteins_df.loc[
            max_distinct_proteins_df["IsNeural"] == neural_condition, ["NumOfProteins"]
        ]
        .sort_values("NumOfProteins")
        .reset_index(drop=True)
    )
    neural_df["CummulativeTranscripts"] = 100 * (neural_df.index + 1) / len(neural_df)
    neural_df["CummulativeTranscripts"] = neural_df["CummulativeTranscripts"][
        ::-1
    ].values
    neural_df = neural_df.drop_duplicates(subset="NumOfProteins").reset_index(drop=True)
    neural_dfs.append(neural_df)

y_min = 1    

fig = make_subplots(
    # rows=1,
    # cols=2,
    rows=2,
    cols=1,
    x_title="Distinct proteins per transcript",
    y_title="% of transcripts",
    shared_yaxes=True,
    shared_xaxes=True,
    # vertical_spacing=facet_row_spacing / 2.5,
    # horizontal_spacing=facet_col_spacing * 1.5,
    vertical_spacing=0.05,
    # horizontal_spacing=0.025,
)

x = df["NumOfProteins"]
y = df["CummulativeTranscripts"]

y_min = min(y_min, y.min())

fig.add_trace(
    go.Scatter(
        x=x,
        y=y,
        mode="lines+markers",
        marker=dict(color="black", size=5),
        line=dict(color="grey", dash="dash"),
        name="All",
    ),
    row=1,
    col=1,
)

fig.add_shape(
    type="rect",
    x0=1,
    y0=0,
    x1=5,
    y1=100,
    line=dict(
        # color="RoyalBlue",
        width=0,
    ),
    # fillcolor="LightSkyBlue",
    fillcolor="orange",
    opacity=0.2,
    row=1,
    col=1,
)

fig.add_shape(
    type="rect",
    x0=5,
    y0=0,
    x1=50,
    y1=100,
    line=dict(
        # color="RoyalBlue",
        width=0,
    ),
    # fillcolor="LightSkyBlue",
    fillcolor="red",
    opacity=0.2,
    row=1,
    col=1,
)

fig.add_trace(
    go.Scatter(
        x=[2.25, 17],
        # y=[80, 85],
        y=[0.3, 0.3],
        text=[
            "~5 isoforms<br>per gene<br>due to<br>alternative<br>splicing",
            #   "Alternative splicing:<br>an average of ~5 isoforms per gene",
            "~50 distinct<br>polypeptides<br>per gene",
        ],
        mode="text",
        showlegend=False,
    ),
    row=1,
    col=1,
)

for neural_condition, neural_trace_name, neural_df in zip(
    neural_conditions, neural_trace_names, neural_dfs
):
    x = neural_df["NumOfProteins"]
    y = neural_df["CummulativeTranscripts"]
    color = neural_color_discrete_map[neural_condition]

    fig.add_trace(
        go.Scatter(
            x=x,
            y=y,
            mode="lines+markers",
            marker=dict(color=color, size=5),
            line=dict(color=color, dash="dash", width=0.5),
            name=neural_trace_name,
        ),
        # row=1, col=2
        row=2,
        col=1,
    )
    
    y_min = min(y_min, y.min())


x = max_distinct_proteins_df.loc[
    max_distinct_proteins_df["IsNeural"] == "Yes", "NumOfProteins"
]
y = max_distinct_proteins_df.loc[
    max_distinct_proteins_df["IsNeural"] == "No", "NumOfProteins"
]
statistic, pv = scipy.stats.mannwhitneyu(x, y)

fig.add_annotation(
    x=np.log(10) / np.log(10),
    y=np.log(1) / np.log(10),
    xref="x",
    yref="y",
    text=f"<b>Mann-Whitney U between<br>neural to non-neural transcripts</b><br>p-val = {pv:.2e}<br>statistic = {statistic:.2g}",
    bgcolor="white",
    borderpad=4,
    font=dict(size=11),
    opacity=0.8,
    showarrow=False,
    row=2,
    col=1,
)

fig.update_xaxes(type="log")
fig.update_yaxes(
    type="log", 
    # range=[-2, 2.2]
    range=[np.log(y_min)*1.1 / np.log(10), 2.2]
)

width = 700
height = 800

fig.update_layout(
    # xaxis_title="Distinct proteins per transcript",
    # yaxis_title="% of transcripts",
    title="Pooled octopus data",
    title_x=0.15,
    template=template,
    width=width,
    height=height,
    # showlegend=False,
)

# fig.write_image(
#     "Distinct proteins per transcript vs. % of transcripts - Octopus.svg",
#     width=width,
#     height=height,
# )

fig.show()

# %% [markdown]
# ### Distinct isoforms per sample

# %%
max_distinct_proteins_df.sort_values("NumOfProteins", ascending=False).iloc[:9]

# %%
# 9 chroms mostly diversified due to A-to-I RNA editing
strongly_diversified_max_distinct_proteins_df = max_distinct_proteins_df.sort_values(
    "NumOfProteins", ascending=False
).iloc[:9]
strongly_diversified_chroms = strongly_diversified_max_distinct_proteins_df[
    "Chrom"
].to_list()
strongly_diversified_max_num_of_proteins = (
    strongly_diversified_max_distinct_proteins_df["NumOfProteins"].to_list()
)
strongly_diversified_transcripts = strongly_diversified_max_distinct_proteins_df[
    "Transcript"
].to_list()

# %%
expanded_unique_proteins_dfs[0]

# %%
expanded_max_distinct_proteins_df = max_distinct_proteins_df.copy()
expanded_max_distinct_proteins_df["Proteins"] = expanded_max_distinct_proteins_df[
    "Proteins"
].str.split(",")
expanded_max_distinct_proteins_df = expanded_max_distinct_proteins_df.explode(
    "Proteins"
).reset_index(drop=True)
expanded_max_distinct_proteins_df = expanded_max_distinct_proteins_df.rename(
    columns={"Proteins": "Protein"}
)
expanded_max_distinct_proteins_df = expanded_max_distinct_proteins_df.drop(
    [
        "NumOfProteins",
        "NumOfReads",
        "MappedReads",
        "Samples",
        "MappedReadsPerSample",
        "DistinctProteins/Reads",
    ],
    axis=1,
)

# expanded_max_distinct_proteins_df

expanded_max_distinct_proteins_dfs = [
    expanded_max_distinct_proteins_df.loc[
        expanded_max_distinct_proteins_df["Chrom"] == chrom
    ]
    .merge(
        expanded_unique_proteins_df, on=["Chrom", "Transcript", "Protein"], how="left"
    )
    .merge(samples_and_tissues_df, on="Sample", how="left")
    for chrom, expanded_unique_proteins_df in zip(chroms, expanded_unique_proteins_dfs)
]

del expanded_max_distinct_proteins_df

for expanded_max_distinct_proteins_df in expanded_max_distinct_proteins_dfs:
    expanded_max_distinct_proteins_df.rename(
        columns={"Sample": "Sample2", "Protein": "Protein2", "Tissue": "Tissue2"},
        inplace=True,
    )
    expanded_max_distinct_proteins_df.insert(
        2, "Sample", expanded_max_distinct_proteins_df["Sample2"]
    )
    expanded_max_distinct_proteins_df.insert(
        3, "Tissue", expanded_max_distinct_proteins_df["Tissue2"]
    )
    expanded_max_distinct_proteins_df.insert(
        4, "Protein", expanded_max_distinct_proteins_df["Protein2"]
    )
    expanded_max_distinct_proteins_df.drop(
        ["Sample2", "Tissue2", "Protein2"], axis=1, inplace=True
    )

expanded_max_distinct_proteins_dfs[1]

# %%
strongly_diversified_expanded_max_distinct_proteins_dfs = [
    expanded_max_distinct_proteins_df
    for expanded_max_distinct_proteins_df, chrom in zip(
        expanded_max_distinct_proteins_dfs, chroms
    )
    if chrom in strongly_diversified_chroms
]

strongly_diversified_expanded_max_distinct_proteins_dfs[0]

# %%
strongly_diversified_num_of_proteins_per_sample_dfs = []

for (
    strongly_diversified_expanded_max_distinct_proteins_df,
    num_of_unique_proteins,
) in zip(
    strongly_diversified_expanded_max_distinct_proteins_dfs,
    strongly_diversified_max_num_of_proteins,
):

    gb = strongly_diversified_expanded_max_distinct_proteins_df.groupby(
        ["Sample", "Tissue"]
    )

    df = gb.apply(len).reset_index().rename(columns={0: "NumOfProteins"})
    df2 = (
        gb.apply(lambda x: 100 * len(x) / num_of_unique_proteins)
        .reset_index()
        .rename(columns={0: "%RelativeNumOfProteins"})
    )
    df = df.merge(df2)

    strongly_diversified_num_of_proteins_per_sample_dfs.append(df)

strongly_diversified_num_of_proteins_per_sample_dfs[0]

# %%
cols = min(facet_col_wrap, len(strongly_diversified_num_of_proteins_per_sample_dfs), 3)
rows = ceil(len(strongly_diversified_num_of_proteins_per_sample_dfs) / cols)
row_col_iter = list(product(range(1, rows + 1), range(1, cols + 1)))[
    : len(strongly_diversified_num_of_proteins_per_sample_dfs)
]

x_title = "Tissue"
y_title = "Distinct isoforms"
# title_text = "Distribution of min & max estimates of non-syn substitutions per read"

# subplot_titles = strongly_diversified_chroms
subplot_titles = [
    f"{transcript.split('_')[0]}<br><sub>({chrom})</sub>"
    for transcript, chrom in zip(
        strongly_diversified_transcripts, strongly_diversified_chroms
    )
]

tissues_order = [
    "Axial nerve cord",
    "Frontal & vertical lobe",
    "Pedunculate & olfactory lobe",
    "Stellate g. & visceral g.",
    "Sucker",
    "Retina & optic lobe",
    "Non-neuronal tissues mix",
]
tissue_to_legendrank = {tissue: x for x, tissue in enumerate(tissues_order, start=1)}

fig = make_subplots(
    rows=rows,
    cols=cols,
    subplot_titles=subplot_titles,
    shared_yaxes=True,
    x_title=x_title,
    y_title=y_title,
    vertical_spacing=0.12,
)

max_y = 0
legend_constructed = False

for (
    (row, col),
    strongly_diversified_num_of_proteins_per_sample_df,
    strongly_diversified_max_num_of_protein,
) in zip(
    row_col_iter,
    strongly_diversified_num_of_proteins_per_sample_dfs,
    strongly_diversified_max_num_of_proteins,
):
    _tissues = strongly_diversified_num_of_proteins_per_sample_df["Tissue"]

    for tissue in _tissues:

        x = [tissue]
        y = strongly_diversified_num_of_proteins_per_sample_df.loc[
            strongly_diversified_num_of_proteins_per_sample_df["Tissue"] == tissue,
            "NumOfProteins",
        ]

        max_y = max(max_y, y.max())

        if not legend_constructed:
            fig.add_trace(
                go.Bar(
                    x=x,
                    y=y,
                    marker_color=tissues_color_discrete_map[tissue],
                    name=tissue,
                    legendrank=tissue_to_legendrank[tissue],
                ),
                row=row,
                col=col,
            )
        else:
            fig.add_trace(
                go.Bar(
                    x=x,
                    y=y,
                    marker_color=tissues_color_discrete_map[tissue],
                    # name=tissue,
                    showlegend=False,
                ),
                row=row,
                col=col,
            )

    legend_constructed = True

#         fig.add_shape(
#             type="line",
#             # xref="paper",
#             # yref="paper",
#             # x0=0.1, x1=4,
#             x0='Axial nerve cord', x1="Retina & optic lobe",
#             y0=strongly_diversified_max_num_of_protein, y1=strongly_diversified_max_num_of_protein,
#             line=dict(
#                 color="grey",
#                 dash="dash",
#                 width=2
#             ),
#             row=row,
#             col=col,
#         )

#         max_y = max(max_y, strongly_diversified_max_num_of_protein)


# fig.update_xaxes(tickangle=45, automargin=True)

fig.update_xaxes(
    showticklabels=False,  # Hide x axis ticks
    categoryorder="array",
    categoryarray=tissues_order,
)
fig.update_yaxes(
    range=[0, max_y * 1.1],
    # range=[0, 0.5+np.log(max_y)/np.log(10)],
    # type="log"
)

width = 1000
height = 800


fig.update_layout(
    template=template,
    title_text="Octopus",
    title_x=0.1,
    title_y=0.97,
    # showlegend=False,
    legend_title_text="Tissue",
    width=width,
    height=height,
)

# fig.write_image(
#     f"{title_text} - PacBio.svg",
#     width=max(350 * cols, 800),
#     height=max(200 * rows, 300),
# )

fig.show()


# %%
strongly_diversified_per_transcript_per_sample_coverage_dfs = [
    df
    for df, chrom in zip(per_transcript_per_sample_coverage_dfs, chroms)
    if chrom in strongly_diversified_chroms
]
ic(len(strongly_diversified_per_transcript_per_sample_coverage_dfs))
strongly_diversified_per_transcript_per_sample_coverage_dfs[0]

# %%
cols = min(facet_col_wrap, len(strongly_diversified_num_of_proteins_per_sample_dfs), 3)
rows = ceil(len(strongly_diversified_num_of_proteins_per_sample_dfs) / cols)
row_col_iter = list(product(range(1, rows + 1), range(1, cols + 1)))[
    : len(strongly_diversified_num_of_proteins_per_sample_dfs)
]

x_title = "Tissue"
y_title = "Mapped reads"
# title_text = "Distribution of min & max estimates of non-syn substitutions per read"

# subplot_titles = strongly_diversified_chroms
subplot_titles = [
    f"{transcript.split('_')[0]}<br><sub>({chrom})</sub>"
    for transcript, chrom in zip(
        strongly_diversified_transcripts, strongly_diversified_chroms
    )
]

tissues_order = [
    "Axial nerve cord",
    "Frontal & vertical lobe",
    "Pedunculate & olfactory lobe",
    "Stellate g. & visceral g.",
    "Sucker",
    "Retina & optic lobe",
    "Non-neuronal tissues mix",
]
tissue_to_legendrank = {tissue: x for x, tissue in enumerate(tissues_order, start=1)}

fig = make_subplots(
    rows=rows,
    cols=cols,
    subplot_titles=subplot_titles,
    shared_yaxes=True,
    x_title=x_title,
    y_title=y_title,
    vertical_spacing=0.12,
)

max_y = 0
legend_constructed = False

for (
    (row, col),
    strongly_diversified_per_transcript_per_sample_coverage_df,
    strongly_diversified_max_num_of_protein,
) in zip(
    row_col_iter,
    strongly_diversified_per_transcript_per_sample_coverage_dfs,
    strongly_diversified_max_num_of_proteins,
):
    _tissues = strongly_diversified_num_of_proteins_per_sample_df["Tissue"]

    for tissue in _tissues:

        x = [tissue]
        y = strongly_diversified_per_transcript_per_sample_coverage_df.loc[
            strongly_diversified_per_transcript_per_sample_coverage_df["Tissue"] == tissue,
            "NumOfReads",
        ]

        max_y = max(max_y, y.max())

        if not legend_constructed:
            fig.add_trace(
                go.Bar(
                    x=x,
                    y=y,
                    marker_color=tissues_color_discrete_map[tissue],
                    name=tissue,
                    legendrank=tissue_to_legendrank[tissue],
                ),
                row=row,
                col=col,
            )
        else:
            fig.add_trace(
                go.Bar(
                    x=x,
                    y=y,
                    marker_color=tissues_color_discrete_map[tissue],
                    # name=tissue,
                    showlegend=False,
                ),
                row=row,
                col=col,
            )

    legend_constructed = True


# fig.update_xaxes(tickangle=45, automargin=True)

fig.update_xaxes(
    showticklabels=False,  # Hide x axis ticks
    categoryorder="array",
    categoryarray=tissues_order,
)
fig.update_yaxes(
    range=[0, max_y * 1.1],
    # range=[0, 0.5+np.log(max_y)/np.log(10)],
    # type="log"
)

width = 1000
height = 800


fig.update_layout(
    template=template,
    title_text="Octopus",
    title_x=0.1,
    title_y=0.97,
    # showlegend=False,
    legend_title_text="Tissue",
    width=width,
    height=height,
)

# fig.write_image(
#     f"{title_text} - PacBio.svg",
#     width=max(350 * cols, 800),
#     height=max(200 * rows, 300),
# )

fig.show()


# %%
cols = min(facet_col_wrap, len(strongly_diversified_num_of_proteins_per_sample_dfs), 3)
rows = ceil(len(strongly_diversified_num_of_proteins_per_sample_dfs) / cols)
row_col_iter = list(product(range(1, rows + 1), range(1, cols + 1)))[
    : len(strongly_diversified_num_of_proteins_per_sample_dfs)
]

x_title = "Mapped reads"
y_title = "Distinct isoforms"

# title_text = "Distribution of min & max estimates of non-syn substitutions per read"

# subplot_titles = strongly_diversified_chroms
subplot_titles = [
    f"{transcript.split('_')[0]}<br><sub>({chrom})</sub>"
    for transcript, chrom in zip(
        strongly_diversified_transcripts, strongly_diversified_chroms
    )
]

tissues_order = [
    "Axial nerve cord",
    "Frontal & vertical lobe",
    "Pedunculate & olfactory lobe",
    "Stellate g. & visceral g.",
    "Sucker",
    "Retina & optic lobe",
    "Non-neuronal tissues mix",
]
tissue_to_legendrank = {tissue: x for x, tissue in enumerate(tissues_order, start=1)}

fig = make_subplots(
    rows=rows,
    cols=cols,
    subplot_titles=subplot_titles,
    shared_yaxes=True,
    x_title=x_title,
    y_title=y_title,
    vertical_spacing=0.12,
    # horizontal_spacing=0.2,
)

max_x = 0
max_y = 0
legend_constructed = False

for (
    (row, col),
    strongly_diversified_num_of_proteins_per_sample_df,
    strongly_diversified_max_num_of_protein,
    strongly_diversified_per_transcript_per_sample_coverage_df
) in zip(
    row_col_iter,
    strongly_diversified_num_of_proteins_per_sample_dfs,
    strongly_diversified_max_num_of_proteins,
    strongly_diversified_per_transcript_per_sample_coverage_dfs
):
    _tissues = strongly_diversified_num_of_proteins_per_sample_df["Tissue"]

    for tissue in _tissues:

        try:
            x = strongly_diversified_per_transcript_per_sample_coverage_df.loc[
                strongly_diversified_per_transcript_per_sample_coverage_df["Tissue"] == tissue,
                "NumOfReads",
            ]
        except KeyError:
            x = [0]
        
        y = strongly_diversified_num_of_proteins_per_sample_df.loc[
            strongly_diversified_num_of_proteins_per_sample_df["Tissue"] == tissue,
            "NumOfProteins",
        ]
        # ic(samp)

        max_x = max(max_x, x.max())
        max_y = max(max_y, y.max())

        if not legend_constructed:
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=y,
                    marker_color=tissues_color_discrete_map[tissue],
                    name=tissue,
                    legendrank=tissue_to_legendrank[tissue],
                    # marker_pattern_shape="/",
                    mode='markers'
                ),
                row=row,
                col=col,
            )
        else:
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=y,
                    marker_color=tissues_color_discrete_map[tissue],
                    # name=tissue,
                    showlegend=False,
                    # marker_pattern_shape="/",
                    mode='markers'
                ),
                row=row,
                col=col,
            )
            
       

    legend_constructed = True


# fig.update_xaxes(tickangle=45, automargin=True)

fig.update_xaxes(
    # showticklabels=False,  # Hide x axis ticks
    # categoryorder="array",
    # categoryarray=tissues_order,
    range=[0, max_x*1.1]
)
# fig.update_yaxes(
#     range=[0, max_y * 1.1],
#     # range=[0, 0.5+np.log(max_y)/np.log(10)],
#     # type="log"
# )
fig.update_yaxes(
    # title_text=primary_y_title, 
    range=[0, max_y * 1.1],
    # secondary_y=False
)
# fig.update_yaxes(
#     title_text=secondary_y_title, 
#     range=[0, max_y_2 * 1.1],
#     secondary_y=True
# )

width = 950
height = 800

fig.update_traces(opacity=0.7, marker_size=6)

fig.update_layout(
    template=template,
    title_text="Octopus",
    title_x=0.1,
    title_y=0.97,
    # showlegend=False,
    legend_title_text="Tissue",
    width=width,
    height=height,
    # barmode="overlay",
)

# fig.write_image(
#     f"{title_text} - PacBio.svg",
#     width=max(350 * cols, 800),
#     height=max(200 * rows, 300),
# )

fig.show()


# %%
cols = min(facet_col_wrap, len(strongly_diversified_num_of_proteins_per_sample_dfs), 3)
rows = ceil(len(strongly_diversified_num_of_proteins_per_sample_dfs) / cols)
row_col_iter = list(product(range(1, rows + 1), range(1, cols + 1)))[
    : len(strongly_diversified_num_of_proteins_per_sample_dfs)
]

x_title = "Tissue"
primary_y_title = "Distinct isoforms"
secondary_y_title = "Mapped reads"
# title_text = "Distribution of min & max estimates of non-syn substitutions per read"

# subplot_titles = strongly_diversified_chroms
subplot_titles = [
    f"{transcript.split('_')[0]}<br><sub>({chrom})</sub>"
    for transcript, chrom in zip(
        strongly_diversified_transcripts, strongly_diversified_chroms
    )
]

tissues_order = [
    "Axial nerve cord",
    "Frontal & vertical lobe",
    "Pedunculate & olfactory lobe",
    "Stellate g. & visceral g.",
    "Sucker",
    "Retina & optic lobe",
    "Non-neuronal tissues mix",
]
tissue_to_legendrank = {tissue: x for x, tissue in enumerate(tissues_order, start=1)}

fig = make_subplots(
    rows=rows,
    cols=cols,
    subplot_titles=subplot_titles,
    shared_yaxes=True,
    x_title=x_title,
    # y_title=y_title,
    vertical_spacing=0.12,
    horizontal_spacing=0.2,
    specs=[
        [{"secondary_y": True} for _ in range(cols)]
        for _ in range(rows)
    ]
)

max_y_1 = 0
max_y_2 = 0
legend_constructed = False

for (
    (row, col),
    strongly_diversified_num_of_proteins_per_sample_df,
    strongly_diversified_max_num_of_protein,
    strongly_diversified_per_transcript_per_sample_coverage_df
) in zip(
    row_col_iter,
    strongly_diversified_num_of_proteins_per_sample_dfs,
    strongly_diversified_max_num_of_proteins,
    strongly_diversified_per_transcript_per_sample_coverage_dfs
):
    _tissues = strongly_diversified_num_of_proteins_per_sample_df["Tissue"]

    for tissue in _tissues:

        # distinct isoforms
        
        x = [tissue]
        y = strongly_diversified_num_of_proteins_per_sample_df.loc[
            strongly_diversified_num_of_proteins_per_sample_df["Tissue"] == tissue,
            "NumOfProteins",
        ]
        # ic(samp)

        max_y_1 = max(max_y_1, y.max())

        if not legend_constructed:
            fig.add_trace(
                go.Bar(
                    x=x,
                    y=y,
                    marker_color=tissues_color_discrete_map[tissue],
                    name=tissue,
                    legendrank=tissue_to_legendrank[tissue],
                    # marker_pattern_shape="/",
                ),
                row=row,
                col=col,
                secondary_y=False,
            )
        else:
            fig.add_trace(
                go.Bar(
                    x=x,
                    y=y,
                    marker_color=tissues_color_discrete_map[tissue],
                    # name=tissue,
                    showlegend=False,
                    # marker_pattern_shape="/",
                ),
                row=row,
                col=col,
                secondary_y=False,
            )
            
       # mapped reads
    
        try:
            y = strongly_diversified_per_transcript_per_sample_coverage_df.loc[
                strongly_diversified_per_transcript_per_sample_coverage_df["Tissue"] == tissue,
                "NumOfReads",
            ]
        except KeyError:
            y = [0]
            
        max_y_2 = max(max_y_2, y.max())
    
        fig.add_trace(
            go.Bar(
            # go.Scatter(
                x=x,
                y=y,
                # marker_color=tissues_color_discrete_map[tissue],
                marker_color="black",
                showlegend=False,
                # name=tissue,
                # legendrank=tissue_to_legendrank[tissue],
                # marker_pattern_shape="",
            ),
            row=row,
            col=col,
            secondary_y=True,
        )

    legend_constructed = True


# fig.update_xaxes(tickangle=45, automargin=True)

fig.update_xaxes(
    showticklabels=False,  # Hide x axis ticks
    categoryorder="array",
    categoryarray=tissues_order,
)
# fig.update_yaxes(
#     range=[0, max_y * 1.1],
#     # range=[0, 0.5+np.log(max_y)/np.log(10)],
#     # type="log"
# )
fig.update_yaxes(
    title_text=primary_y_title, 
    range=[0, max_y_1 * 1.1],
    secondary_y=False
)
fig.update_yaxes(
    title_text=secondary_y_title, 
    range=[0, max_y_2 * 1.1],
    secondary_y=True
)

width = 1200
height = 800

fig.update_traces(opacity=0.5)

fig.update_layout(
    template=template,
    title_text="Octopus",
    title_x=0.1,
    title_y=0.97,
    # showlegend=False,
    legend_title_text="Tissue",
    width=width,
    height=height,
    barmode="overlay",
)

# fig.write_image(
#     f"{title_text} - PacBio.svg",
#     width=max(350 * cols, 800),
#     height=max(200 * rows, 300),
# )

fig.show()


# %%
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Create sample data
x = [1, 2, 3, 4, 5]
y1 = [10, 20, 30, 40, 50]
y2 = [100, 200, 300, 400, 500]

# Create subplots with shared y-axis
fig = make_subplots(rows=1, cols=2, shared_yaxes=True)

# Add traces to subplots
fig.add_trace(go.Scatter(x=x, y=y1, name='Trace 1', yaxis='y1'), row=1, col=1)
fig.add_trace(go.Scatter(x=x, y=y2, name='Trace 2', yaxis='y2'), row=1, col=2)

# Retrieve tick values from primary y-axis
primary_tickvals = fig['layout']['yaxis1']['tickvals']

# Set tick values for secondary y-axis
fig['layout']['yaxis2']['tickvals'] = primary_tickvals

# Update layout
fig.update_layout(title='Subplots with Synchronized Y-Axis Tick Values')

# Show the figure
fig.show()

# %%
# # Get the names of all y-axes
# y_axes = []
# for axis in fig['layout']['yaxis']:
#     if axis.startswith('y'):
#         y_axes.append(fig['layout']['yaxis'][axis]['title']['text'])

# print(y_axes)



# %%
# fig['layout']

# %%

# %%
# fig = make_subplots(rows=2, cols=2,
#                     # shared_yaxes=True,
#                     specs=[[{"secondary_y": True}, {"secondary_y": True}],
#                            [{"secondary_y": True}, {"secondary_y": True}]])

# # Top left
# fig.add_trace(
#     go.Scatter(x=[1, 2, 3], y=[2, 52, 62], name="yaxis data"),
#     row=1, col=1, secondary_y=False)

# fig.add_trace(
#     go.Scatter(x=[1, 2, 3], y=[40, 50, 60], name="yaxis2 data"),
#     row=1, col=1, secondary_y=True,
# )

# # Top right
# fig.add_trace(
#     go.Scatter(x=[1, 2, 3], y=[2, 52, 62], name="yaxis3 data"),
#     row=1, col=2, secondary_y=False,
# )

# fig.add_trace(
#     go.Scatter(x=[1, 2, 3], y=[40, 50, 60], name="yaxis4 data"),
#     row=1, col=2, secondary_y=True,
# )

# # Bottom left
# fig.add_trace(
#     go.Scatter(x=[1, 2, 3], y=[2, 52, 62], name="yaxis5 data"),
#     row=2, col=1, secondary_y=False,
# )

# fig.add_trace(
#     go.Scatter(x=[1, 2, 3], y=[40, 50, 60], name="yaxis6 data"),
#     row=2, col=1, secondary_y=True,
# )

# # Bottom right
# fig.add_trace(
#     go.Scatter(x=[1, 2, 3], y=[2, 52, 62], name="yaxis7 data"),
#     row=2, col=2, secondary_y=False,
# )

# fig.add_trace(
#     go.Scatter(x=[1, 2, 3], y=[40, 50, 60], name="yaxis8 data"),
#     row=2, col=2, secondary_y=True,
# )

# fig.update_yaxes(
#     title_text=primary_y_title, 
#     # range=[0, max_y * 1.1],
#     # range=[0, 0.5+np.log(max_y)/np.log(10)],
#     # type="log"
#     secondary_y=False
# )
# fig.update_yaxes(title_text=secondary_y_title, secondary_y=True)

# width = 1000
# height = 800


# fig.update_layout(
#     template=template,
#     title_text="Octopus",
#     title_x=0.1,
#     title_y=0.97,
#     # showlegend=False,
#     legend_title_text="Tissue",
#     width=width,
#     height=height,
# )

# fig.show()

# %%

# %%

# %%

# %%

# %% [markdown]
# ### NaNs distribution

# %%
ambigous_positions_df = (
    expanded_distinct_unique_proteins_df.loc[
        expanded_distinct_unique_proteins_df["Algorithm"] == "Descending"
    ]
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
        & (expanded_distinct_unique_proteins_df["Algorithm"] == "Descending")
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
df3.loc[df3["#SolutionIncluded"] == 0].groupby(condition_col).size()

# %%
df3.loc[df3["#SolutionIncluded"] > 0].groupby(condition_col).size()

# %%
df3.loc[df3["#SolutionIncluded"] == 100].groupby(condition_col).size()

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


# %% [markdown]
# ### Comparing algorithms

# %% [markdown]
# #### Coupled diff comparisons

# %%
distinct_unique_proteins_df

# %%
# distinct_unique_proteins_df.to_csv("O.vul.DistinctProteins.tsv", sep="\t", index=False)

# %%
# distinct_unique_proteins_df.loc[distinct_unique_proteins_df["Fraction"] == 1.0].sort_values("NumOfProteins", ascending=False)

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
# mean distinct proteins per transcript
max_distinct_proteins_per_transcript_and_alg_df.sort_values(
    [condition_col, "NumOfProteins"], ascending=False
).drop_duplicates(condition_col, ignore_index=True)["NumOfProteins"].mean()

# %%
# num of transcripts with at least 5 variants
max_distinct_proteins_per_transcript_and_alg_df.sort_values(
    [condition_col, "NumOfProteins"], ascending=False
).drop_duplicates(condition_col, ignore_index=True)["NumOfProteins"].ge(5).sum()

# %%
# num of transcripts with at least 50 variants
max_distinct_proteins_per_transcript_and_alg_df.sort_values(
    [condition_col, "NumOfProteins"], ascending=False
).drop_duplicates(condition_col, ignore_index=True)["NumOfProteins"].ge(50).sum()

# %%
asc_df = max_distinct_proteins_per_transcript_and_alg_df.loc[
    max_distinct_proteins_per_transcript_and_alg_df["Algorithm"] == "Ascending"
].reset_index(drop=True)
desc_df = max_distinct_proteins_per_transcript_and_alg_df.loc[
    max_distinct_proteins_per_transcript_and_alg_df["Algorithm"] != "Ascending"
].reset_index(drop=True)

ic(len(asc_df))
ic(len(desc_df))

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
df.loc[df["Fraction"] < 1.0]

# %%
fig = px.histogram(df, x="Desc - Asc", template=template, text_auto=True)

fig.update_layout(
    # showlegend=False,
    # width=1500,
    # height=300
    width=600,
    height=400,
)

fig.show()

# %%
fig = px.histogram(
    df,
    x="Desc - Asc",
    template=template,
    # text_auto=True,
    facet_col="Fraction",
    # color="Fraction",
    category_orders={"Fraction": [0.2, 0.4, 0.6, 0.8, 1.0]},
    # barmode="group"
)

# Reduce opacity to see both histograms
fig.update_traces(opacity=0.75)

fig.update_layout(
    # showlegend=False,
    width=1300,
    # height=300
    # width=800,
    height=400,
    # barmode='overlay' # Overlay both histograms
)


fig.show()

# %%

# %%

# %%

# %%

# %%

# %%
# ic(asc_df["Fraction"].eq(desc_df["Fraction"]).all())
# ic(asc_df["FractionRepetition"].eq(desc_df["FractionRepetition"]).all())
# ic(asc_df["NumOfReads"].eq(desc_df["NumOfReads"]).all())
# ic(asc_df["AlgorithmRepetition"].eq(desc_df["AlgorithmRepetition"]).all());

# %%
# algs_diff_df = asc_df.rename(
#     columns={
#         "NumOfProteins": "NumOfProteinsAscending",
#         "Proteins": "ProteinsAscending",
#         "IsMaxNumOfProteins": "IsMaxNumOfProteinsAscending",
#     }
# )
# del algs_diff_df["Algorithm"]
# algs_diff_df["NumOfProteinsDescending"] = desc_df["NumOfProteins"]
# algs_diff_df["ProteinsDescending"] = desc_df["Proteins"]
# algs_diff_df["IsMaxNumOfProteinsDescending"] = desc_df["IsMaxNumOfProteins"]

# algs_diff_df["Desc - Asc"] = (
#     algs_diff_df["NumOfProteinsDescending"] - algs_diff_df["NumOfProteinsAscending"]
# )
# algs_diff_df["Diff"] = algs_diff_df["Desc - Asc"].apply(
#     lambda x: "> 0" if x > 0 else "= 0" if x == 0 else "< 0"
# )

# # algs_diff_df["IsMaxNumOfProteins"] = algs_diff_df.apply(
# #     lambda x: "Both"
# #     if (x["IsMaxNumOfProteinsAscending"] and x["IsMaxNumOfProteinsDescending"])
# #     else "Ascending"
# #     if x["IsMaxNumOfProteinsAscending"]
# #     else "Descending"
# #     if x["IsMaxNumOfProteinsDescending"]
# #     else "None",
# #     axis=1,
# # )

# algs_diff_df["IsMaxNumOfProteins"] = algs_diff_df.apply(
#     lambda x: "Both"
#     if (x["IsMaxNumOfProteinsAscending"] and x["IsMaxNumOfProteinsDescending"])
#     else "Yes"
#     if x["IsMaxNumOfProteinsAscending"]
#     else "Yes"
#     if x["IsMaxNumOfProteinsDescending"]
#     else "No",
#     axis=1,
# )

# algs_diff_df["IsMaxNumOfProteins2"] = algs_diff_df.apply(
#     lambda x: "Both"
#     if (x["IsMaxNumOfProteinsAscending"] and x["IsMaxNumOfProteinsDescending"])
#     else "Ascending"
#     if x["IsMaxNumOfProteinsAscending"]
#     else "Descending"
#     if x["IsMaxNumOfProteinsDescending"]
#     else "None",
#     axis=1,
# )

# # algs_diff_df["IsMaxNumOfProteins"] = algs_diff_df.apply(
# #     lambda x: True
# #     if (x["IsMaxNumOfProteinsAscending"] or x["IsMaxNumOfProteinsDescending"])
# #     else False,
# #     axis=1,
# # )

# algs_diff_df

# %%
# c  # df = algs_diff_df.groupby([condition_col, "Diff"])[["IsMaxNumOfProteins"]].value_counts().reset_index().rename(columns={0: "Count"})
# # df = df.rename(columns={"Diff": "#Desc <span>&#8722;</span> #Asc"})

# fig = px.histogram(
#     algs_diff_df.loc[algs_diff_df["Desc - Asc"]],
#     x="Desc - Asc",
#     # y="Diff",
#     # color="IsMaxNumOfProteins",
#     # opacity=0.5,
#     # color_discrete_map=color_discrete_map,
#     # facet_row=condition_col,
#     # facet_row_spacing=facet_col_spacing,
#     # facet_col="#Desc <span>&#8722;</span> #Asc",
#     # facet_col_spacing=facet_col_spacing,
#     # symbol="IsMaxNumOfProteins",
#     # category_orders=category_orders,
#     template=template,
#     # title="Differences between coupled solutions' sizes",
#     # shared_
#     # log_x=True,
#     # log_y=True,
#     # marginal_y="box",
#     # marginal_x="box",
#     # text_auto=True
# )

# # https://www.toptal.com/designers/htmlarrows/math/minus-sign/
# # fig.update_yaxes(title="Desc <span>&#8722;</span>  Asc")
# # https://community.plotly.com/t/changing-label-of-plotly-express-facet-categories/28066/5?u=kaparanewbie
# # fig.for_each_annotation(lambda a: a.update(text=a.text.replace("#Desc <span>&#8722;</span> #Asc=", "#Desc <span>&#8722;</span> #Asc ")))
# # fig.for_each_annotation(lambda a: a.update(text=a.text.replace(f"{condition_col}=", "")))

# # fig.update_traces(marker={"size": 4})

# # fig.update_xaxes(visible=True)
# # fig.update_xaxes(zeroline=True)
# # fig.update_xaxes(zerolinewidth=4)
# # fig.update_yaxes(zerolinewidth=zerolinewidth)

# fig.update_layout(
#     # showlegend=False,
#     # width=1500,
#     # height=300
#     width=600,
#     height=500,
# )

# fig.show()

# %%
# df = algs_diff_df.rename(columns={"Diff": "#Desc <span>&#8722;</span> #Asc"})
# df = df.loc[df["Desc - Asc"] < 0]

# # x = df["Desc - Asc"]
# # y_asc = df["NumOfProteinsAscending"]
# # y_desc = df["NumOfProteinsDescending"]

# fig = px.scatter_3d(
#     df,
#     x="Desc - Asc",
#     y="NumOfProteinsAscending",
#     z="NumOfProteinsDescending",
#     # y="Diff",
#     color="IsMaxNumOfProteins",
#     labels={"NumOfProteinsAscending": "Asc", "NumOfProteinsDescending": "Desc"},
#     # opacity=0.5,
#     # color_discrete_map=color_discrete_map,
#     # facet_row=condition_col,
#     # facet_row_spacing=facet_col_spacing,
#     # facet_col="#Desc <span>&#8722;</span> #Asc",
#     # facet_col_spacing=facet_col_spacing,
#     # symbol="IsMaxNumOfProteins",
#     # category_orders=category_orders,
#     template=template,
#     # title="Differences between coupled solutions' sizes",
#     # shared_
#     # log_x=True,
#     # log_y=True,
#     # marginal_y="box",
#     # marginal_x="box",
#     # text_auto=True
# )

# # https://www.toptal.com/designers/htmlarrows/math/minus-sign/
# # fig.update_yaxes(title="Desc <span>&#8722;</span>  Asc")
# # https://community.plotly.com/t/changing-label-of-plotly-express-facet-categories/28066/5?u=kaparanewbie
# # fig.for_each_annotation(lambda a: a.update(text=a.text.replace("#Desc <span>&#8722;</span> #Asc=", "#Desc <span>&#8722;</span> #Asc ")))
# # fig.for_each_annotation(lambda a: a.update(text=a.text.replace(f"{condition_col}=", "")))

# # fig.update_traces(marker={"size": 4})

# # fig.update_xaxes(visible=True)
# # fig.update_xaxes(zeroline=True)
# # fig.update_xaxes(zerolinewidth=4)
# # fig.update_yaxes(zerolinewidth=zerolinewidth)

# fig.update_layout(
#     # showlegend=False,
#     # width=1500,
#     # height=300
#     # width=800,
#     height=600,
# )

# fig.show()

# %%
# # df = algs_diff_df.groupby([condition_col, "Diff"])[["IsMaxNumOfProteins"]].value_counts().reset_index().rename(columns={0: "Count"})
# df = algs_diff_df.rename(columns={"Diff": "#Desc <span>&#8722;</span> #Asc"})
# df = df.loc[df["Desc - Asc"] < 0]

# # Create figure with secondary y-axis
# fig = make_subplots(specs=[[{"secondary_y": True}]])

# x = df["Desc - Asc"]
# y_asc = df["NumOfProteinsAscending"]
# y_desc = df["NumOfProteinsDescending"]

# # # Add traces
# # fig.add_trace(
# #     go.Scatter(x=x, y=y_asc, name="NumOfProteinsAscending", showlegend=False),
# #     secondary_y=False,
# # )

# # fig.add_trace(
# #     go.Scatter(x=x, y=y_desc, name="NumOfProteinsDescending data", marker=dict(color="white"), showlegend=False),
# #     secondary_y=True,
# # )

# trace1 = go.Scatter(x=x, y=y_asc, name="NumOfProteinsAscending", showlegend=False)
# trace2 = go.Scatter(x=x, y=y_asc, name="NumOfProteinsDescending data", showlegend=False)

# # The first trace is referenced to the default xaxis, yaxis (ie. xaxis='x1', yaxis='y1')
# fig.add_trace(trace1, secondary_y=False)

# # The second trace is referenced to xaxis='x1'(i.e. 'x1' is common for the two traces)
# # and yaxis='y2' (the right side yaxis)

# fig.add_trace(trace2, secondary_y=True)

# # Set x-axis title
# fig.update_xaxes(title_text="#Desc <span>&#8722;</span> #Asc")

# # Set y-axes titles
# fig.update_yaxes(
#     title_text="NumOfProteinsAscending", secondary_y=False, range=[0, y_asc.max() * 1.1]
# )
# fig.update_yaxes(
#     title_text="NumOfProteinsDescending",
#     secondary_y=True,
#     range=[0, y_desc.max() * 1.1],
#     showgrid=False,
# )

# # fig.update_traces(marker={"size": 4})

# # fig.update_xaxes(visible=True)
# # fig.update_xaxes(zeroline=True)
# # fig.update_xaxes(zerolinewidth=4)
# # fig.update_yaxes(zerolinewidth=zerolinewidth)

# fig.update_layout(
#     # showlegend=False,
#     # width=1500,
#     # height=300
#     width=500,
#     height=500,
#     template=template,
# )

# fig.show()

# %%
# df = (
#     algs_diff_df.groupby([condition_col, "Diff"])[["IsMaxNumOfProteins"]]
#     .value_counts()
#     .reset_index()
#     .rename(columns={0: "Count"})
# )
# df = df.rename(columns={"Diff": "#Desc <span>&#8722;</span> #Asc"})

# fig = px.histogram(
#     df,
#     x="#Desc <span>&#8722;</span> #Asc",
#     # y="Diff",
#     color="IsMaxNumOfProteins",
#     # opacity=0.5,
#     # color_discrete_map=color_discrete_map,
#     # facet_row=condition_col,
#     # facet_row_spacing=facet_col_spacing,
#     # facet_col="#Desc <span>&#8722;</span> #Asc",
#     facet_col_spacing=facet_col_spacing,
#     # symbol="IsMaxNumOfProteins",
#     # category_orders=category_orders,
#     template=template,
#     # title="Differences between coupled solutions' sizes",
#     # shared_
#     # log_x=True,
#     # log_y=True,
#     # marginal_y="box",
#     # marginal_x="box",
#     text_auto=True,
# )

# # https://www.toptal.com/designers/htmlarrows/math/minus-sign/
# # fig.update_yaxes(title="Desc <span>&#8722;</span>  Asc")
# # https://community.plotly.com/t/changing-label-of-plotly-express-facet-categories/28066/5?u=kaparanewbie
# fig.for_each_annotation(
#     lambda a: a.update(
#         text=a.text.replace(
#             "#Desc <span>&#8722;</span> #Asc=", "#Desc <span>&#8722;</span> #Asc "
#         )
#     )
# )
# # fig.for_each_annotation(lambda a: a.update(text=a.text.replace(f"{condition_col}=", "")))

# # fig.update_traces(marker={"size": 4})

# # fig.update_xaxes(visible=True)
# # fig.update_xaxes(zeroline=True)
# # fig.update_xaxes(zerolinewidth=4)
# # fig.update_yaxes(zerolinewidth=zerolinewidth)

# fig.update_layout(
#     # showlegend=False,
#     # width=1500,
#     # height=300
#     width=800,
#     height=400,
# )

# fig.show()

# %%
# df = (
#     algs_diff_df.groupby([condition_col, "Diff"])[["IsMaxNumOfProteins"]]
#     .value_counts()
#     .reset_index()
#     .rename(columns={0: "Count"})
# )
# df = df.rename(columns={"Diff": "#Desc <span>&#8722;</span> #Asc"})

# fig = px.histogram(
#     df,
#     x="IsMaxNumOfProteins",
#     # y="Diff",
#     # color="#Desc <span>&#8722;</span> #Asc",
#     # opacity=0.5,
#     # color_discrete_map=color_discrete_map,
#     # facet_row=condition_col,
#     # facet_row_spacing=facet_col_spacing,
#     facet_col="#Desc <span>&#8722;</span> #Asc",
#     facet_col_spacing=facet_col_spacing,
#     # symbol="IsMaxNumOfProteins",
#     # category_orders=category_orders,
#     template=template,
#     # title="Differences between coupled solutions' sizes",
#     # shared_
#     # log_x=True,
#     # log_y=True,
#     # marginal_y="box",
#     # marginal_x="box",
#     text_auto=True,
# )

# # https://www.toptal.com/designers/htmlarrows/math/minus-sign/
# # fig.update_yaxes(title="Desc <span>&#8722;</span>  Asc")
# # https://community.plotly.com/t/changing-label-of-plotly-express-facet-categories/28066/5?u=kaparanewbie
# fig.for_each_annotation(
#     lambda a: a.update(
#         text=a.text.replace(
#             "#Desc <span>&#8722;</span> #Asc=", "#Desc <span>&#8722;</span> #Asc "
#         )
#     )
# )
# # fig.for_each_annotation(lambda a: a.update(text=a.text.replace(f"{condition_col}=", "")))

# # fig.update_traces(marker={"size": 4})

# # fig.update_xaxes(visible=True)
# # fig.update_xaxes(zeroline=True)
# # fig.update_xaxes(zerolinewidth=4)
# # fig.update_yaxes(zerolinewidth=zerolinewidth)

# fig.update_layout(
#     # showlegend=False,
#     # width=1500,
#     # height=300
#     width=800,
#     height=400,
# )

# fig.show()

# %%
# df = algs_diff_df.rename(columns={"Diff": "#Desc <span>&#8722;</span> #Asc"})
# # df["abs(Desc - Asc) div 10"] = df["Desc - Asc"].abs().div(10).add(0.1)
# df["abs(Desc - Asc)"] = df["Desc - Asc"].abs().add(1)

# fig = px.scatter(
#     df,
#     x="NumOfProteinsAscending",
#     y="NumOfProteinsDescending",
#     color="#Desc <span>&#8722;</span> #Asc",
#     size="abs(Desc - Asc)",
#     # opacity=0.5,
#     # color_discrete_map=color_discrete_map,
#     # facet_row=condition_col,
#     # facet_row_spacing=facet_col_spacing,
#     facet_col="IsMaxNumOfProteins",
#     # facet_col_spacing=facet_col_spacing,
#     # symbol="IsMaxNumOfProteins",
#     # category_orders=category_orders,
#     template=template,
#     # title="Differences between coupled solutions' sizes",
#     # shared_
#     log_x=True,
#     log_y=True,
#     marginal_y="box",
#     marginal_x="box",
#     # trendline="lowess"
# )

# # https://www.toptal.com/designers/htmlarrows/math/minus-sign/
# # fig.update_yaxes(title="Desc <span>&#8722;</span>  Asc")
# # https://community.plotly.com/t/changing-label-of-plotly-express-facet-categories/28066/5?u=kaparanewbie
# # fig.for_each_annotation(lambda a: a.update(text=a.text.replace("Diff=", "Diff ")))
# # fig.for_each_annotation(lambda a: a.update(text=a.text.replace(f"{condition_col}=", "")))

# # fig.update_traces(marker={"size": 4})

# # fig.update_xaxes(visible=True)
# # fig.update_xaxes(zeroline=True)
# # fig.update_xaxes(zerolinewidth=4)
# # fig.update_yaxes(zerolinewidth=zerolinewidth)

# fig.update_layout(
#     # showlegend=False,
#     # width=1500,
#     # height=300
#     width=1400,
#     height=600,
# )

# fig.show()

# %%
# df = algs_diff_df.rename(columns={"Diff": "#Desc <span>&#8722;</span> #Asc"})

# fig = px.scatter(
#     df,
#     x="NumOfProteinsAscending",
#     y="NumOfProteinsDescending",
#     color="#Desc <span>&#8722;</span> #Asc",
#     # color_discrete_map=color_discrete_map,
#     # facet_row=condition_col,
#     # facet_row_spacing=facet_col_spacing,
#     facet_col="Fraction",
#     # facet_col_spacing=facet_col_spacing,
#     # symbol="Diff",
#     # category_orders=category_orders,
#     template=template,
#     # title="Differences between coupled solutions' sizes",
#     # shared_
#     log_x=True,
#     log_y=True,
#     marginal_y="box",
#     marginal_x="box",
# )

# # https://www.toptal.com/designers/htmlarrows/math/minus-sign/
# # fig.update_yaxes(title="Desc <span>&#8722;</span>  Asc")
# # https://community.plotly.com/t/changing-label-of-plotly-express-facet-categories/28066/5?u=kaparanewbie
# # fig.for_each_annotation(lambda a: a.update(text=a.text.replace("Diff=", "Diff ")))
# # fig.for_each_annotation(lambda a: a.update(text=a.text.replace(f"{condition_col}=", "")))

# fig.update_traces(marker={"size": 5})

# # fig.update_xaxes(visible=True)
# # fig.update_xaxes(zeroline=True)
# # fig.update_xaxes(zerolinewidth=4)
# fig.update_yaxes(zerolinewidth=zerolinewidth)

# fig.update_layout(
#     # showlegend=False,
#     # width=1500,
#     # height=300
#     # width=550,
#     # height=500
# )

# fig.show()

# %%
# df = algs_diff_df.rename(columns={"Desc - Asc": "#Desc <span>&#8722;</span> #Asc"})

# fig = px.box(
#     df,
#     x="Fraction",
#     y="#Desc <span>&#8722;</span> #Asc",
#     # color=condition_col,
#     # color_discrete_map=color_discrete_map,
#     # facet_row=condition_col,
#     # facet_row_spacing=facet_col_spacing,
#     facet_col="Diff",
#     # facet_col_spacing=facet_col_spacing,
#     # symbol="Diff",
#     # category_orders=category_orders,
#     template=template,
#     title="Differences between coupled solutions' sizes",
#     # shared_
# )
# fig.update_layout(showlegend=False)
# # https://www.toptal.com/designers/htmlarrows/math/minus-sign/
# # fig.update_yaxes(title="Desc <span>&#8722;</span>  Asc")
# # https://community.plotly.com/t/changing-label-of-plotly-express-facet-categories/28066/5?u=kaparanewbie
# fig.for_each_annotation(lambda a: a.update(text=a.text.replace("Diff=", "Diff ")))
# # fig.for_each_annotation(lambda a: a.update(text=a.text.replace(f"{condition_col}=", "")))

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
# https://www.toptal.com/designers/htmlarrows/math/minus-sign/
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

# basic_diff_names = ["Diff < 0", "Diff = 0", "Diff > 0"]
# subplot_titles = []
# first_row = True
# for condition in conditions:
#     for diff in basic_diff_names:
#         st = f"<sub>{condition}</sub>"
#         if first_row:
#             st = f"{diff}<br>" + st
#         subplot_titles.append(st)
#     first_row = False

# fig = make_subplots(
#     rows=len(conditions),
#     cols=3,
#     y_title=y_axis_name,
#     x_title=x_axis_name,
#     # subplot_titles=["Diff < 0", "Diff = 0", "0 < Diff"],
#     subplot_titles=subplot_titles,
#     shared_yaxes=True,
#     shared_xaxes=True,
# )

# algorithms = ["Ascending", "Descending"]
# diffs = ["< 0", "= 0", "> 0"]
# symbols = ["cross", "diamond", "circle"]

# for col, (diff, symbol) in enumerate(zip(diffs, symbols), start=1):
#     for row, condition in enumerate(conditions, start=1):
#         df = algs_diff_df.loc[
#             (algs_diff_df["Diff"] == diff) & (algs_diff_df[condition_col] == condition)
#         ]
#         x = df["Fraction"]
#         y = df["Desc - Asc"]
#         if col == 1:
#             fig.add_trace(
#                 go.Scatter(
#                     x=x,
#                     y=y,
#                     legendgrouptitle_text=condition,
#                     legendgroup=condition,
#                     name=condition,
#                     line_color=color_discrete_map[condition],
#                     mode="markers",
#                     marker={"symbol": symbol},
#                 ),
#                 row=row,
#                 col=col,
#             )
#         else:
#             fig.add_trace(
#                 go.Scatter(
#                     x=x,
#                     y=y,
#                     # legendgrouptitle_text=condition,
#                     legendgroup=condition,
#                     name=condition,
#                     line_color=color_discrete_map[condition],
#                     mode="markers",
#                     marker={"symbol": symbol},
#                 ),
#                 row=row,
#                 col=col,
#             )

# fig.update_layout(
#     title_text=head_title, title_y=0.95, template=template, showlegend=False
# )

# fig.update_xaxes(tick0=0, dtick=0.1, tickangle=-60, matches="x")

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

# fig.show()


# %% [markdown]
# #### Solutions' sizes

# %%
distinct_unique_proteins_df

# %%
# min_max_fraction_1_distinct_prots_df = (
#     distinct_unique_proteins_df.loc[distinct_unique_proteins_df["Fraction"] == 1.0]
#     .groupby(condition_col)["NumOfProteins"]
#     .agg(["min", "max"])
#     .reset_index()
# )
# min_max_fraction_1_distinct_prots_df[
#     "%SolutionsDispersion"
# ] = min_max_fraction_1_distinct_prots_df.apply(
#     lambda x: 100 * (x["max"] - x["min"]) / x["max"], axis=1
# )
# min_max_fraction_1_distinct_prots_df["HighDispersion"] = (
#     min_max_fraction_1_distinct_prots_df["%SolutionsDispersion"] > 1
# )
# min_max_fraction_1_distinct_prots_df = min_max_fraction_1_distinct_prots_df.sort_values(
#     [condition_col, "min", "max", "%SolutionsDispersion", "HighDispersion"]
# ).reset_index(drop=True)
# min_max_fraction_1_distinct_prots_df

# %%
dispersion_df = distinct_unique_proteins_df.loc[
    distinct_unique_proteins_df["Fraction"] == 1.0,
    ["Chrom", condition_col, "NumOfReads", "NumOfProteins"],
].reset_index(drop=True)

gb = dispersion_df.groupby("Chrom")["NumOfProteins"]
dispersion_df["MaxNumOfProteins"] = gb.transform(max)
dispersion_df["MinNumOfProteins"] = gb.transform(min)

dispersion_df = dispersion_df.drop("NumOfProteins", axis=1)
dispersion_df = dispersion_df.drop_duplicates(ignore_index=True)

dispersion_df["%SolutionsDispersion"] = dispersion_df.apply(
    lambda x: 100
    * (x["MaxNumOfProteins"] - x["MinNumOfProteins"])
    / x["MaxNumOfProteins"],
    axis=1,
)
dispersion_df["HighDispersion"] = dispersion_df["%SolutionsDispersion"] > 1

dispersion_df

# %%
# len(dispersion_df.loc[dispersion_df["HighDispersion"]])

# %%
# fig = px.scatter(
#     dispersion_df,
#     x="MinNumOfProteins",
#     y="MaxNumOfProteins",
#     # size="%SolutionsDispersion",
#     # facet_col="HighDispersion"
#     color="HighDispersion",
#     marginal_y="box",
#     marginal_x="box",
#     # log_y=True
# )

# # fig.update_xaxes(title="% dispersion<br><sub>100 * (max - min) / max</sub>")
# # fig.update_yaxes(title="Transcripts", type="log")

# fig.update_layout(
#     # showlegend=False,+
#     width=600,
#     height=500,
#     template=template,
# )

# # fig.write_image(
# #     "%SolutionsDispersion - Octopus.svg",
# #     width=600,
# #     height=400
# # )

# fig.show()

# %%
# fig = px.scatter(
#     dispersion_df.loc[dispersion_df["%SolutionsDispersion"] > 0],
#     x="NumOfReads",
#     y="%SolutionsDispersion",
#     # size="%SolutionsDispersion",
#     # facet_col="HighDispersion",
#     color="HighDispersion",
#     # color="NumOfReads",
#     # marginal_y="box",
#     # marginal_x="box",
#     # log_x=True,
#     log_y=True,
# )

# # fig.update_xaxes(title="% dispersion<br><sub>100 * (max - min) / max</sub>")
# # fig.update_yaxes(title="Transcripts", type="log")

# fig.update_layout(
#     # showlegend=False,+
#     width=600,
#     height=500,
#     template=template,
# )

# # fig.write_image(
# #     "%SolutionsDispersion - Octopus.svg",
# #     width=600,
# #     height=400
# # )

# fig.show()

# %%
# fig = px.scatter(
#     dispersion_df.loc[dispersion_df["%SolutionsDispersion"] > 0],
#     x="MinNumOfProteins",
#     y="MaxNumOfProteins",
#     # size="%SolutionsDispersion",
#     facet_col="HighDispersion",
#     # color="HighDispersion",
#     color="NumOfReads",
#     # marginal_y="box",
#     # marginal_x="box",
#     # log_y=True
# )

# # fig.update_xaxes(title="% dispersion<br><sub>100 * (max - min) / max</sub>")
# # fig.update_yaxes(title="Transcripts", type="log")

# fig.update_layout(
#     # showlegend=False,+
#     width=1200,
#     height=500,
#     template=template,
# )

# # fig.write_image(
# #     "%SolutionsDispersion - Octopus.svg",
# #     width=600,
# #     height=400
# # )

# fig.show()

# %%
fig = go.Figure(
    go.Histogram(
        x=dispersion_df["%SolutionsDispersion"],
        marker_color="black",
    )
)

# fig.update_xaxes(title="% dispersion<br><sub>100 * (max - min) / max</sub>")
fig.update_xaxes(title="% dispersion of distinct proteins sets' sizes")
fig.update_yaxes(title="Transcripts", type="log")

fig.update_layout(
    # showlegend=False,
    title="Pooled octopus data",
    title_x=0.15,
    width=600,
    height=400,
    template=template,
)

fig.write_image("%SolutionsDispersion - Octopus.svg", width=600, height=400)

fig.show()

# %%
# y_axis_name = "Distinct unique proteins"
# head_title = "Distinct unique proteins vs. heuristic method"

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

#     df = distinct_unique_proteins_df.loc[
#         distinct_unique_proteins_df[condition_col] == condition
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


# %% [markdown] papermill={"duration": 0.030615, "end_time": "2022-02-01T09:42:49.024262", "exception": false, "start_time": "2022-02-01T09:42:48.993647", "status": "completed"} toc-hr-collapsed=true
# ## Supporting reads' coverage

# %% [markdown] papermill={"duration": 0.030615, "end_time": "2022-02-01T09:42:49.024262", "exception": false, "start_time": "2022-02-01T09:42:48.993647", "status": "completed"}
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

# %% [markdown] papermill={"duration": 0.030615, "end_time": "2022-02-01T09:42:49.024262", "exception": false, "start_time": "2022-02-01T09:42:48.993647", "status": "completed"}
# ### Distinct unique reads

# %%
# expanded_distinct_unique_reads_df


# %%
# (
#     expanded_distinct_unique_reads_df.query("Fraction == 1.0")
#     .groupby([condition_col, "UniqueRead"])["NumOfReads"]
#     .mean()
#     .max()
# )


# %%
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
unique_proteins_dfs[0].head()


# %%
max(df["NumOfReads"].max() for df in unique_proteins_dfs)


# %% papermill={"duration": 4.052404, "end_time": "2022-02-01T09:42:53.176715", "exception": false, "start_time": "2022-02-01T09:42:49.124311", "status": "completed"}
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
y_col_name = "TotalWeightedSupportingReads"
expression_df = (
    expression_dfs[0].sort_values(y_col_name, ascending=False).reset_index(drop=True)
)
expression_df["CummulativeRelativeWeightedExpression"] = expression_df.groupby(
    "#Solution"
)[[y_col_name]].transform(lambda x: 100 * x / x.sum())
expression_df = expression_df.loc[expression_df["#Solution"] == "1000"].reset_index(
    drop=True
)
expression_df

# %%
for n_top_expressed in [10, 100, 1000]:
    fig = go.Figure()
    x = [f"{n_top_expressed} top expressed"] * n_top_expressed + ["Rest"] * (
        len(expression_df) - n_top_expressed
    )
    y = expression_df["MinNonSyns"]
    fig.add_trace(
        go.Box(
            x=x,
            y=y,
            # boxpoints='all,
            # mode="markers",
            # marker=dict(
            #     # size=16,
            #     # cmax=39,
            #     # cmin=0,
            #     color=z,
            #     colorbar=dict(
            #         title="MinNonSyns"
            #     ),
            #     # colorscale="Viridis"
            # ),
        )
    )
    # fig.update_xaxes(type="log")
    # fig.update_yaxes(type="log")
    fig.update_layout(
        height=400,
        template=template,
        yaxis_title="MinNonSyns",
    )
    # fig.update_traces(boxpoints='all')
    fig.show()

# %%
fig = go.Figure()
# x = expression_df["MinNonSynsFrequency"]
x = expression_df.index + 1
y = expression_df["CummulativeRelativeWeightedExpression"]
# z = expression_df["MinNonSynsFrequency"]
z = expression_df["MinNonSyns"]
fig.add_trace(
    go.Scattergl(
        x=x,
        y=y,
        mode="markers",
        marker=dict(
            # size=16,
            # cmax=39,
            # cmin=0,
            color=z,
            colorbar=dict(title="MinNonSyns"),
            # colorscale="Viridis"
        ),
    )
)
fig.update_xaxes(type="log")
fig.update_yaxes(type="log")
fig.update_layout(height=400, template=template)
fig.show()


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
    find_rand_maximal_solution(expression_df, seed, allowed_algorithms=["Descending"])
    for expression_df in expression_dfs
]
maximal_solutions

# %%
percentile_dfs = [
    make_percentile_df(
        expression_df.loc[expression_df["#Solution"] == maximal_solution].reset_index(
            drop=True
        ),
        allowed_algorithms=["Descending"],
    )
    for expression_df, maximal_solution in zip(expression_dfs, maximal_solutions)
]
percentile_dfs[0]

# %%
_seeds = [np.random.default_rng(seed)]
for _ in conditions[1:]:
    _seeds.append(np.random.default_rng(_seeds[-1]))

all_conditions_sample_solutions = [
    choose_sample_solutions(expression_df, _seed, allowed_algorithms=["Descending"])
    for expression_df, _seed in zip(expression_dfs, _seeds)
]
all_conditions_sample_solutions[0]

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
gb.agg({"Diff5+": ["size", "sum"]})

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
#         "TotalWeightedSupportingReads": "Total weighted<br>supporting reads",
#     },
#     height=500,
#     width=600,
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
#     vertical_spacing=facet_row_spacing / 2.5,
#     horizontal_spacing=facet_col_spacing * 1.5,
# )

# for (row, col), percentile_df, condition in zip(
#     row_col_iter, percentile_dfs, conditions
# ):
#     legend_x = [percentile_df["RequiredProteins"].max() * 5 / 6]
#     legend_ys = [[25], [20]]

#     for color, symbol, assignment_method, legend_y in zip(
#         subcolors_discrete_map[condition], symbols, assignment_methods, legend_ys
#     ):
#         _percentile_df = percentile_df.loc[percentile_df["AssignmentMethod"] == assignment_method]

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
#                 textfont=dict(size=8)
#             ),
#             row=row,
#             col=col,
#         )

# fig.update_xaxes(
#     tick0 = 0,
#     dtick = 5_000,
#     matches='x'
# )
# fig.update_layout(
#     title=head_title,
#     showlegend=False,
#     template=template,
# )
# # fig.write_image(
# #     f"{head_title} - PacBio.svg",
# #     height=max(300, 200 * rows),
# #     width=max(600, 250 * cols),
# # )
# fig.show()

# %%
# # only Weighted assignment method for poster

# assignment_method = "Weighted"
# y_col_name = "TotalWeightedSupportingReads"

# x_axis_name = "Distinct unique protein rank"
# y_axis_name = "Cummulative relative<br>expression (%)"
# # head_title = f"Distinct unique proteins vs. {assignment_method.lower()} cummulative expression (POSTER)"
# head_title = f"Weighted cummulative expression vs. distinct protein rank"

# cols = min(facet_col_wrap, len(conditions), 4)
# rows = ceil(len(conditions) / cols)
# row_col_iter = list(product(range(1, rows + 1), range(1, cols + 1)))[: len(conditions)]

# fig = make_subplots(
#     rows=rows,
#     cols=cols,
#     subplot_titles=conditions,
#     shared_yaxes=True,
#     x_title=x_axis_name,
#     y_title=y_axis_name,
#     vertical_spacing=facet_row_spacing / 2.5,
#     horizontal_spacing=facet_col_spacing * 1.5,
# )

# for (row, col), percentile_df, condition, expression_df in zip(
#     row_col_iter, percentile_dfs, conditions, expression_dfs
# ):
#     _percentile_df = percentile_df.loc[percentile_df["AssignmentMethod"] == assignment_method]

#     x = _percentile_df["RequiredProteins"]
#     y = _percentile_df["Percentile"]

#     x_mean = _percentile_df.groupby("Percentile")["RequiredProteins"].apply(np.mean)
#     y_unique = x_mean.index

#     expression_df = expression_df.sort_values(y_col_name, ascending=False).reset_index(drop=True)
#     expression_df["CummulativeRelativeWeightedExpression"] = expression_df.groupby("#Solution")[[y_col_name]].transform(lambda x: 100 * x / x.sum())
#     gb = expression_df.groupby("#Solution")["CummulativeRelativeWeightedExpression"]
#     # mean_cumm_exp_10_most_frequent = expression_df.groupby("#Solution")["CummulativeRelativeWeightedExpression"].apply(lambda x: x[:10].sum()).mean()
#     mean_cumm_exp_10_most_frequent = gb.apply(lambda x: x[:10].sum()).mean()
#     mean_cumm_exp_100_most_frequent = gb.apply(lambda x: x[:100].sum()).mean()
#     mean_cumm_exp_1000_most_frequent = gb.apply(lambda x: x[:1000].sum()).mean()


#     # all data
#     fig.add_trace(
#         go.Scatter(
#             x=x,
#             y=y,
#             mode="markers",
#             marker=dict(
#                 color=color_discrete_map[condition],
#                 size=4,
#                 opacity=0.7,
#                 # line=dict(width=0),
#             ),
#         ),
#         row=row,
#         col=col,
#     )

#     # mean data
#     fig.add_trace(
#         go.Scatter(
#             x=x_mean,
#             y=y_unique,
#             mode="lines+markers",
#             marker=dict(
#                 color=color_discrete_map[condition],
#                 size=4,
#                 opacity=0.7,
#                 # line=dict(width=0),
#             ),
#         ),
#         row=row,
#         col=col,
#     )

#     # mean cummulative exp of 10 most expressed proteins in each solution
#     fig.add_trace(
#         go.Scatter(
#             x=[10, 100, 1000],
#             y=[mean_cumm_exp_10_most_frequent, mean_cumm_exp_100_most_frequent, mean_cumm_exp_1000_most_frequent],
#             mode="markers+text",
#             marker=dict(
#                 color="black",
#                 size=6,
#                 # opacity=0.7,
#                 # symbol="triangle-up",
#                 symbol="square",
#                 # line=dict(width=0),
#             ),
#             text=[
#                 f"  (10, {mean_cumm_exp_10_most_frequent:.1f})",
#                 f"   (100, {mean_cumm_exp_100_most_frequent:.1f})",
#                 f"    (1000, {mean_cumm_exp_1000_most_frequent:.1f})"
#             ],
#             textposition="middle right",
#             textfont=dict(size=12)
#         ),
#         row=row,
#         col=col,
#     )

# fig.update_layout(
#     title=head_title,
#     showlegend=False,
#     template=template,
# )
# fig.update_xaxes(
#     tick0 = 0,
#     dtick = 5_000,
#     matches='x',
#     # type="log"
# )
# fig.update_yaxes(
#     tick0 = 0,
#     dtick = 20,
#     matches='y',
#     # type="log"
# )
# # fig.write_image(
# #     f"{head_title} - PacBio.svg",
# #     height=max(400, 200 * rows),
# #     width=max(650, 250 * cols),
# # )
# fig.show()

# %%
x_axis_name = "Distinct unique protein rank"
y_axis_name = "Cummulative relative<br>expression (%)"
head_title = f"Weighted cummulative expression vs. distinct protein rank"

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
    title=head_title,
    showlegend=False,
    template=template,
    height=max(300, 260 * rows),
    width=max(650, 250 * cols),
)
# fig.write_image(
#     f"{head_title} - Illumina.svg",
#     height=max(300, 200 * rows),
#     width=max(600, 250 * cols),
# )
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

# %%
# cols = min(facet_col_wrap, len(conditions), 3)
# rows = ceil(len(conditions) / cols)
# row_col_iter = list(product(range(1, rows + 1), range(1, cols + 1)))[: len(conditions)]

# linear_spaces = [(300, 15_000), (250, 23_000)]  # (start, end) tuples for both x and y
# forward_transforms = [(linear_to_log10, linear_to_log10), (linear_to_log10, inverse)]  # (x, y) tuples
# reverse_transforms = [(log10_to_linear, log10_to_linear), (log10_to_linear, inverse)]  # (x, y) tuples
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
#     for condition, solution, algorithm in zip(conditions, maximal_solutions, maximal_algorithms)
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
#     horizontal_spacing=facet_col_spacing * 1.5
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
#     formulate_equation
# ) in zip(
#     row_col_iter,
#     conditions,
#     maximal_dfs,
#     maximal_solutions,
#     maximal_algorithms,
#     linear_spaces,
#     forward_transforms,
#     reverse_transforms,
#     formulate_equations
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
#                         opacity=data_opacity,
#                         line=dict(width=0),
#                         symbol=symbol
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
#                         opacity=data_opacity,
#                         line=dict(width=0),
#                         symbol=symbol
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
#                 showlegend=False
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
#              f"<b>{assignment_method}</b>"
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
#             font=dict(
#                 size=8,
#                 color=color
#             ),
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
assignment_method = "Weighted"
y_col_name = "TotalWeightedSupportingReads"

cols = min(facet_col_wrap, len(conditions), 3)
rows = ceil(len(conditions) / cols)
row_col_iter = list(product(range(1, rows + 1), range(1, cols + 1)))[: len(conditions)]

linear_spaces = [(300, 15_000), (250, 23_000)]  # (start, end) tuples for both x and y
forward_transforms = [
    (linear_to_log10, linear_to_log10),
    (linear_to_log10, inverse),
]  # (x, y) tuples
reverse_transforms = [
    (log10_to_linear, log10_to_linear),
    (log10_to_linear, inverse),
]  # (x, y) tuples
# formulate_equations = [formulate_log10_equation, formulate_semilog10_equation]
fit_texts = ["    y ~ 1 / sqrt(x)", "    y ~ 1 / log(x)"]

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
x_axis_name = "Distinct unique protein rank"
y_axis_name = "Relative expression (%)"
head_title = f"Relative expression of proteins considering a largest solution in each {str(condition_col).lower()}"

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
    maximal_df,
    maximal_solution,
    linear_space,
    (forward_x_transform, forward_y_transform),
    (reverse_x_transform, reverse_y_transform),
    fit_text,
) in zip(
    row_col_iter,
    conditions,
    maximal_dfs,
    maximal_solutions,
    linear_spaces,
    forward_transforms,
    reverse_transforms,
    fit_texts,
):

    assignment_df = maximal_df.sort_values(
        "TotalWeightedSupportingReads", ascending=False
    ).reset_index(drop=True)
    assignment_df["#Protein"] = list(range(1, len(assignment_df) + 1))
    assignment_df["AssignmentMethod"] = assignment_method

    x = assignment_df["#Protein"]
    y = 100 * assignment_df[y_col_name] / assignment_df[y_col_name].sum()

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

    fig.add_annotation(
        row=row,
        col=col,
        x=text_x,
        y=text_y,
        xref="x",
        yref="y",
        text=fit_text,
        align="center",
        font=dict(size=12, color="grey"),
        showarrow=False,
    )

fig.update_layout(
    title_text=head_title,
    # title_y=0.95,
    template=template,
    showlegend=False,
    # legend_itemsizing="constant",
    height=max(400, 200 * rows),
    # width=max(900, 250 * cols),
)
fig.update_xaxes(type="log", nticks=6)
fig.update_yaxes(type="log")
fig.write_image(
    f"{head_title} - PacBio.svg",
    height=max(400, 200 * rows),
    width=max(650, 250 * cols),
)
fig.show()
# fig.show(config={'staticPlot': True, 'responsive': False})


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


# %% [markdown]
# > All proteins

# %%

# %%

# %% [markdown]
# > TEST START

# %%
weighted_exp_tsne_input_dfs[0].iloc[:, ML_INPUT_FIRST_COL_POS:]

# %%
# perplexities = [5, 30, 50, 100]
# perplexities = [150]
n_iter = 500
n_jobs = 20
n_iter_500_weighted_conditions_tsnes, n_iter_500_weighted_conditions_Xs = run_tsnes(
    [conditions[0]],
    [weighted_exp_tsne_input_dfs[0]],
    seed,
    perplexities=[150],
    n_iter=n_iter,
    n_jobs=n_jobs,
)

# %%
prots_perplexity_tsne = n_iter_500_weighted_conditions_tsnes[0][0]
X_transformed = n_iter_500_weighted_conditions_Xs[0]

# %%
X = weighted_exp_tsne_input_dfs[0].iloc[:, ML_INPUT_FIRST_COL_POS:].values
X

# %%
# perplexities = [5, 30, 50, 100]
# perplexities = [150]
n_iter = 500
n_jobs = 20
top_expressed_proteins = 1000
(
    top_1000_n_iter_500_weighted_conditions_tsnes,
    top_1000_n_iter_500_weighted_conditions_Xs,
) = run_tsnes(
    [conditions[0]],
    [weighted_exp_tsne_input_dfs[0]],
    seed,
    perplexities=[150],
    n_iter=n_iter,
    n_jobs=n_jobs,
    top_expressed_proteins=top_expressed_proteins,
)

# %%
top_1000_prots_perplexity_tsne = top_1000_n_iter_500_weighted_conditions_tsnes[0][0]
top_1000_X_transformed = top_1000_n_iter_500_weighted_conditions_Xs[0]

# %%
top_1000_X = (
    weighted_exp_tsne_input_dfs[0]
    .iloc[:top_expressed_proteins, ML_INPUT_FIRST_COL_POS:]
    .values
)
top_1000_X

# %%
# cluster_sizes = [2**x for x in range(2, 8)]
cluster_sizes = [2**x for x in range(2, 5)]
cluster_sizes

# %%
kmodes_clusters = [
    KModes(
        n_clusters=n_clusters, init="Huang", n_init=5, verbose=1, n_jobs=n_jobs
    ).fit_predict(X)
    for n_clusters in cluster_sizes
]

# %%
cluster_sizes = [2**x for x in range(2, 8)]
# cluster_sizes = [2**x for x in range(2, 5)]
cluster_sizes

# %%
list(range(10, 110, 10))

# %%
# cluster_sizes = [2**x for x in range(2, 8)]
cluster_sizes = list(range(10, 110, 10))  # 10, 20, ..., 100
top_1000_kmodes = [
    KModes(n_clusters=n_clusters, init="Huang", n_init=5, verbose=0, n_jobs=n_jobs).fit(
        top_1000_X
    )
    for n_clusters in cluster_sizes
]
top_1000_kmodes_sumofsq = [km.cost_ for km in top_1000_kmodes]
top_1000_kmodes_clusters = [km.predict(top_1000_X) for km in top_1000_kmodes]

# %%
cluster_sizes_2 = list(range(110, 160, 10))  # 110, 120, ..., 150
top_1000_kmodes_2 = [
    KModes(n_clusters=n_clusters, init="Huang", n_init=5, verbose=0, n_jobs=n_jobs).fit(
        top_1000_X
    )
    for n_clusters in cluster_sizes_2
]
top_1000_kmodes_sumofsq_2 = [km.cost_ for km in top_1000_kmodes_2]
top_1000_kmodes_clusters_2 = [km.predict(top_1000_X) for km in top_1000_kmodes_2]

# %%
fig = px.scatter(
    x=cluster_sizes + cluster_sizes_2,
    y=top_1000_kmodes_sumofsq + top_1000_kmodes_sumofsq_2,
    template=template,
    labels={"x": "Number of clusters (k)", "y": "Sum of square distances"},
    title="Elbow method for optimal number of clusters",
    # log_x=True
    # log_y=True
)
fig.update_layout(width=700, height=500)
fig.update_yaxes(range=[0, max(top_1000_kmodes_sumofsq) * 1.1])
fig.show()

# %%
len(set(top_1000_kmodes_clusters[-1]))

# %%
fig = go.Figure()
x, y = top_1000_prots_perplexity_tsne.T
fig.add_trace(
    go.Scattergl(
        x=x,
        y=y,
        mode="markers",
        marker=dict(
            # color="white",
            color=top_1000_kmodes_clusters[-1],
            line_width=line_width,
            # size=marker_size
        ),
    ),
    # row=row,
    # col=col
)
fig.update_layout(width=700, height=500, template=template)
# fig.update_yaxes(range=[0, max(top_1000_kmodes_sumofsq)*1.1])
fig.show()

# %%
cols = min(len(cluster_sizes), 5)
rows = ceil(len(cluster_sizes) / cols)
row_col_iter = list(product(range(1, rows + 1), range(1, cols + 1)))[
    : len(cluster_sizes)
]

head_title = "t-SNE for top 1000 expressed proteins in GRIA, colors by kmodes of different n_clusters"
# column_titles = cluster_sizes

fig = make_subplots(
    rows=rows,
    cols=cols,
    subplot_titles=[f"K = {k}" for k in cluster_sizes],
    # shared_yaxes=True,
    # shared_xaxes=True,
    # x_title="PC1",
    # y_title="PC2",
)

# color_sequence = px.colors.qualitative.Light24
x, y = top_1000_prots_perplexity_tsne.T

for (row, col), cluster, cluster_size in zip(
    row_col_iter, top_1000_kmodes_clusters, cluster_sizes
):

    # colormap = {label: color for label, color in zip(range(cluster_size), color_sequence[:cluster_size])}
    # colors = [colormap[label] for label in cluster]

    # rank_cutoff = 1000
    # fig = go.Figure()
    marker_size = 3
    line_width = 0
    # n = X.shape[0]
    # color_options = [color_discrete_map[condition], "white"]
    # colors = color_highest_expressed_proteins(n, rank_cutoff, color_options)

    fig.add_trace(
        go.Scattergl(
            x=x,
            y=y,
            mode="markers",
            marker=dict(
                # color="white",
                # color=colors,
                color=cluster,
                colorscale="Electric",
                line_width=line_width,
                size=marker_size,
            ),
        ),
        row=row,
        col=col,
    )
pixles = 800
fig.update_layout(
    title_text=head_title,
    # title_y=0.95,
    template=template,
    showlegend=False,
    # width=600,
    # height=1800,
    width=1400,
    height=600,
)

fig.show()

# %%
cols = min(len(clusters), 3)
rows = ceil(len(clusters) / cols)
row_col_iter = list(product(range(1, rows + 1), range(1, cols + 1)))[: len(clusters)]

head_title = "t-SNE for GRIA, colors by kmodes of different n_clusters"
# column_titles = cluster_sizes

fig = make_subplots(
    rows=rows,
    cols=cols,
    subplot_titles=cluster_sizes,
    # shared_yaxes=True,
    # shared_xaxes=True,
    # x_title="PC1",
    # y_title="PC2",
)

color_sequence = px.colors.qualitative.Light24
x, y = prots_perplexity_tsne.T

for (row, col), cluster, cluster_size in zip(
    row_col_iter, kmodes_clusters, cluster_sizes
):

    colormap = {
        label: color
        for label, color in zip(range(cluster_size), color_sequence[:cluster_size])
    }
    colors = [colormap[label] for label in cluster]

    # rank_cutoff = 1000
    # fig = go.Figure()
    marker_size = 2
    line_width = 0.5
    # n = X.shape[0]
    # color_options = [color_discrete_map[condition], "white"]
    # colors = color_highest_expressed_proteins(n, rank_cutoff, color_options)

    fig.add_trace(
        go.Scattergl(
            x=x,
            y=y,
            mode="markers",
            marker=dict(
                # color="white",
                color=colors,
                line_width=line_width,
                size=marker_size,
            ),
        ),
        row=row,
        col=col,
    )
pixles = 800
fig.update_layout(
    title_text=head_title,
    # title_y=0.95,
    template=template,
    showlegend=False,
    # width=600,
    # height=1800,
    # width=1800,
    # height=600,
)

fig.show()

# %%
test_tsne = TSNE(
    n_components=2,
    learning_rate="auto",
    perplexity=150,
    n_iter=300,
    init="random",
    random_state=np.random.RandomState(seed),
    n_jobs=20,
)
test_tsne = test_tsne.fit_transform(test_X)

# %%
km = KModes(n_clusters=4, init="Huang", n_init=5, verbose=1)

clusters = km.fit_predict(test_X)

# Print the cluster centroids
print(km.cluster_centroids_)

# %%
test_colormap = {
    label: color for label, color in zip(range(4), ["red", "blue", "green", "yellow"])
}
colors = [test_colormap[label] for label in km.labels_]
colors[:6]

# %%
# prots_perplexity_tsne = n_iter_500_weighted_conditions_tsnes[0][0]
# X = n_iter_500_weighted_conditions_Xs[0][0]

# rank_cutoff = 1000
fig = go.Figure()
# marker_size = 1
# line_width = 0.5
# n = X.shape[0]
# color_options = [color_discrete_map[condition], "white"]
# colors = color_highest_expressed_proteins(n, rank_cutoff, color_options)
x, y = test_tsne.T
fig.add_trace(
    go.Scattergl(
        x=x,
        y=y,
        mode="markers",
        marker=dict(
            color=colors,
            # line_width=line_width,
            # size=marker_size
        ),
    ),
)
fig.update_layout(
    # title_text=head_title,
    # title_y=0.95,
    template=template,
    # showlegend=False,
    width=550,
    height=550,
)

fig.show()

# %%
km.labels_

# %%
len(km.labels_)

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
perplexities = [5, 30, 50, 100, 150]
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
rank_cutoff = 1000

head_title = (
    f"t-SNEs for largest solution of each {str(condition_col).lower()} under different perplexities"
    # "<br>"
    # f"<sub>{rank_cutoff} highest expressed proteins are colored</sub>"
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

marker_size = 1
line_width = 0.5

for row, (condition, X, condition_tsnes) in enumerate(
    zip(
        conditions,
        n_iter_500_weighted_conditions_Xs,
        n_iter_500_weighted_conditions_tsnes,
    ),
    start=1,
):

    # n = X.shape[0]
    # color_options = [color_discrete_map[condition], "white"]
    # colors = color_highest_expressed_proteins(n, rank_cutoff, color_options)
    colors = "white"

    for col, prots_perplexity_tsne in enumerate(condition_tsnes, start=1):

        x, y = prots_perplexity_tsne.T

        fig.add_trace(
            go.Scattergl(
                x=x,
                y=y,
                mode="markers",
                marker=dict(color=colors, line_width=line_width, size=marker_size),
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
# equal_conditions_pcas = run_pcas(conditions, equal_exp_tsne_input_dfs, seed)
weighted_conditions_pcas = run_pcas(conditions, weighted_exp_tsne_input_dfs, seed)

# %%
rank_cutoff = 1000

cols = min(facet_col_wrap, len(conditions), 4)
rows = ceil(len(conditions) / cols)
row_col_iter = list(product(range(1, rows + 1), range(1, cols + 1)))[: len(conditions)]

head_title = (
    f"PCAs for largest solution of each {str(condition_col).lower()}"
    "<br>"
    f"<sub>{rank_cutoff} highest expressed proteins are colored</sub>"
)
column_titles = conditions

fig = make_subplots(
    rows=rows,
    cols=cols,
    subplot_titles=conditions,
    # shared_yaxes=True,
    # shared_xaxes=True,
    x_title="PC1",
    y_title="PC2",
)

marker_size = 1

for (row, col), condition, condition_pca in zip(
    row_col_iter, conditions, weighted_conditions_pcas
):

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
            marker=dict(color=colors, line_width=0.5, size=marker_size),
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
# top_1000_equal_conditions_tsnes, top_1000_equal_conditions_Xs = run_tsnes(
#     conditions,
#     equal_exp_tsne_input_dfs,
#     seed,
#     perplexities=perplexities,
#     top_expressed_proteins=top_expressed_proteins,
# )
top_1000_weighted_conditions_tsnes, top_1000_weighted_conditions_Xs = run_tsnes(
    conditions,
    weighted_exp_tsne_input_dfs,
    seed,
    perplexities=perplexities,
    top_expressed_proteins=top_expressed_proteins,
)

# %%
rank_cutoff = 1000

# head_title = (
#     f"t-SNEs for largest solution of each {str(condition_col).lower()} under different perplexities"
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

marker_size = 1

for row, (condition, X, condition_tsnes) in enumerate(
    zip(
        conditions, top_1000_weighted_conditions_Xs, top_1000_weighted_conditions_tsnes
    ),
    start=1,
):

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
                marker=dict(color=colors, line_width=0.5, size=marker_size),
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
top_expressed_proteins = 1000
# top_1000_equal_conditions_pcas = run_pcas(
#     conditions,
#     equal_exp_tsne_input_dfs,
#     seed,
#     top_expressed_proteins=top_expressed_proteins,
# )
top_1000_weighted_conditions_pcas = run_pcas(
    conditions,
    weighted_exp_tsne_input_dfs,
    seed,
    top_expressed_proteins=top_expressed_proteins,
)

# %%
rank_cutoff = 100

for conditions_pcas, sorting_method in zip(
    [top_1000_equal_conditions_pcas, top_1000_weighted_conditions_pcas],
    ["equal", "weighted"],
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

    for col, (condition, condition_pca) in enumerate(
        zip(conditions, conditions_pcas), start=1
    ):

        n = len(condition_pca)
        color_options = [color_discrete_map[condition], "white"]
        colors = color_highest_expressed_proteins(n, rank_cutoff, color_options)

        x = condition_pca["PC1"]
        y = condition_pca["PC2"]

        fig.add_trace(
            go.Scattergl(
                x=x, y=y, mode="markers", marker=dict(color=colors, line_width=0.5)
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
    )
)

fig.update_layout(
    title=f"Shannon's entropy of a largest solution of each {condition_col.lower()}",
    width=500,
    height=700,
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

fig.update_yaxes(title_text="Entropy")
fig.update_xaxes(
    tickfont=dict(size=10)
)  # https://plotly.com/python/axes/#set-axis-label-rotation-and-font

# fig.update_traces(
#     marker=dict(line_color="black", line_width=0.3, pattern_fillmode="replace"),
#     # width=0.3
# )
fig.show()

# %% [markdown]
# ##### Combinatorics

# %%
top_100_combinatorics_df = (
    max_sol_exp_dfs[1]
    .sort_values("%WeightedTotalExpression", ascending=False)
    .reset_index(drop=True)
    .iloc[:100, ML_INPUT_FIRST_COL_POS:]
    .fillna(-1)
    # .T
)
top_100_combinatorics_df

# %%
print("ho")

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
    title="Top 100 expressed proteins<br>in PCLO",
    autorange="reversed",
    tick0=1,
    range=[1, 100],
)
fig.update_layout(height=700, width=650, template=template, font_size=16)

fig.write_image(
    "Combinatorics of top 100 expressed proteins in PCLO - Flash talk - PacBio.svg",
    width=650,
    height=700,
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

# %%
cols = min(facet_col_wrap, len(conditions), 4)
rows = ceil(len(conditions) / cols)
row_col_iter = list(product(range(1, rows + 1), range(1, cols + 1)))[: len(conditions)]

x_title = "Non-syn mutations per protein"
y_title = "Proteins"
title_text = "Distribution of min & max estimates of non-syn mutations per protein"

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

fig.update_layout(
    template=template,
    barmode="overlay",  # Overlay both histograms
    title_text=title_text,
    title_y=0.95,
    showlegend=False,
    height=max(200 * rows, 300),
    width=max(350 * cols, 800),
)

fig.update_traces(opacity=0.75)  # Reduce opacity to see both histograms
fig.update_xaxes(range=[min_x * 0.9, max_x * 1.1])
fig.update_yaxes(range=[0, max_y * 1.2])

fig.show()

