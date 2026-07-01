# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.15.2
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
# # %autosave 600

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
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
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
from EditingUtils.seq import make_fasta_dict

# %% [markdown]
# # Data loading

# %%
condition_col = "Transcript"

orfs_bed = "/private7/projects/Combinatorics/O.vulgaris/Annotations/orfs_oct.bed"
# alignment_stats_file = "/private7/projects/Combinatorics/O.vulgaris/Alignment/PRJNA791920/IsoSeq/AggregatedByChromBySampleSummary.tsv"
alignment_stats_file = "/private7/projects/Combinatorics/O.vulgaris/Alignment/PRJNA791920/IsoSeq.Polished.Unclustered/AggregatedByChromBySampleSummary.tsv"

known_sites_file = (
    "/private7/projects/Combinatorics/O.vulgaris/Annotations/O.vul.EditingSites.csv"
)
transcriptome_file = (
    "/private7/projects/Combinatorics/O.vulgaris/Annotations/orfs_oct.fa"
)

main_data_dir = Path("/private7/projects/Combinatorics/O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BHAfterNoise/")

positions_dir = Path(main_data_dir, "PositionsFiles")
reads_dir = Path(main_data_dir, "ReadsFiles")
proteins_dir = Path(main_data_dir, "ProteinsFiles")
distinct_proteins_dir = Path(main_data_dir, "DistinctProteins2")
# expression_dir = Path(main_data_dir, "ExpressionLevels")

neural_vs_non_neural_expression_file = Path(
    "/private7/projects/Combinatorics/O.vulgaris/Annotations/NeuralVsNonNeuralExpression.csv"
)

samples_and_tissues_file = Path(
    "/private7/projects/Combinatorics/O.vulgaris/Data/PRJNA791920/IsoSeqPolished/samples.csv"
)

# reads_first_col_pos = 7
# unique_reads_first_col_pos = 9
# proteins_first_col_pos = 13
unique_proteins_first_col_pos = 15
reads_editing_col = "EditingFrequency"
proteins_editing_col = "MinNonSyns"

robo2_chrom = "comp182237_c0_seq56"

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
tmr50_alignment_stats_df = alignment_stats_df.loc[alignment_stats_df["MappedReads"] >= 50]
tmr50_alignment_stats_df

# %%
tmr1000_alignment_stats_df = alignment_stats_df.loc[alignment_stats_df["MappedReads"] >= 1000]
tmr1000_alignment_stats_df

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
tmr50_alignment_stats_df.loc[~tmr50_alignment_stats_df["Chrom"].isin(positions_data_df["Chrom"])]

# %%

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
tmr50_alignment_stats_df.loc[~tmr50_alignment_stats_df["Chrom"].isin(reads_data_df["Chrom"])]

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

# expression_files = list(
#     expression_dir.glob("*.DistinctUniqueProteins.ExpressionLevels.csv")
# )
# chroms_in_expression_files = [
#     expression_file.name.split(".")[0] for expression_file in expression_files
# ]


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

# expression_data_df = pd.DataFrame(
#     {"Chrom": chroms_in_expression_files, "ExpressionFile": expression_files}
# )

proteins_data_df = (
    proteins_data_df.merge(unique_proteins_data_df, on="Chrom", how="left")
    .merge(distinct_proteins_data_df, on="Chrom", how="left")
    # .merge(expression_data_df, on="Chrom", how="left")
)

proteins_data_df

# %%
data_df = (
    # orfs_df.merge(alignment_stats_df, on="Chrom", how="left")
    orfs_df.merge(tmr50_alignment_stats_df, on="Chrom", how="right")
    .merge(positions_data_df, on="Chrom", how="left")
    .merge(reads_data_df, on="Chrom", how="left")
    .merge(proteins_data_df, on="Chrom", how="left")
)

data_df

# %%
possibly_na_conditions = data_df["Name"].tolist()
possibly_na_chroms = data_df["Chrom"].tolist()
possibly_na_starts = data_df["Start"].tolist()
possibly_na_ends = data_df["End"].tolist()
possibly_na_strands = data_df["Strand"].tolist()

possibly_na_positions_files = data_df["PositionsFile"].tolist()
possibly_na_reads_files = data_df["ReadsFile"].tolist()
possibly_na_unique_reads_files = data_df["UniqueReadsFile"].tolist()
possibly_na_proteins_files = data_df["ProteinsFile"].tolist()
possibly_na_unique_proteins_files = data_df["UniqueProteinsFile"].tolist()
possibly_na_distinct_unique_proteins_files = data_df["DistinctProteinsFile"].tolist()
# expression_files = complete_data_df["ExpressionFile"].tolist()

# %%
# assert (
#     data_df.loc[data_df["ExpressionFile"].notna()].reset_index(drop=True).shape
#     == data_df.loc[data_df["DistinctProteinsFile"].notna()].reset_index(drop=True).shape
# ), "some distinct proteins don't have expression levels"

# %%
# complete_data_df = data_df.loc[data_df["ExpressionFile"].notna()].reset_index(drop=True)
complete_data_df = data_df.loc[data_df["DistinctProteinsFile"].notna()].reset_index(drop=True)

# complete_data_df = complete_data_df.drop_duplicates(
#     "Name", keep=False, ignore_index=True
# )
complete_data_df

# %%
# complete_data_df[["Chrom"]].to_csv("TMR50.CompleteData.Chroms.tsv", sep="\t", index=False)

# %%
robo2_index = complete_data_df.loc[complete_data_df["Chrom"] == robo2_chrom].index[0]
robo2_index

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
# expression_files = complete_data_df["ExpressionFile"].tolist()

# %%
len(data_df["UniqueReadsFile"])

# %%
len(complete_data_df["UniqueReadsFile"])

# %% [markdown]
# # Data loading - TMR 1000

# %%
tmr1000_main_data_dir = Path("/private7/projects/Combinatorics/O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000.BQ30.AHL.BHAfterNoise/")

tmr1000_positions_dir = Path(tmr1000_main_data_dir, "PositionsFiles")
tmr1000_reads_dir = Path(tmr1000_main_data_dir, "ReadsFiles")
tmr1000_proteins_dir = Path(tmr1000_main_data_dir, "ProteinsFiles")
tmr1000_distinct_proteins_dir = Path(tmr1000_main_data_dir, "DistinctProteins2")

# %%
tmr1000_positions_files = list(tmr1000_positions_dir.glob("*.positions.csv.gz"))

tmr1000_chroms_in_positions = [
    positions_file.name.split(".")[0] for positions_file in tmr1000_positions_files
]

tmr1000_positions_data_df = pd.DataFrame(
    {
        "Chrom": tmr1000_chroms_in_positions,
        "PositionsFile": tmr1000_positions_files,
    }
)

tmr1000_positions_data_df

# %%
tmr1000_reads_files = list(tmr1000_reads_dir.glob("*.reads.csv.gz"))
tmr1000_chroms_in_reads_files = [
    reads_file.name.split(".")[0] for reads_file in tmr1000_reads_files
]

tmr1000_unique_reads_files = list(tmr1000_reads_dir.glob("*.unique_reads.csv.gz"))
tmr1000_chroms_in_unique_reads_files = [
    unique_reads_file.name.split(".")[0]
    for unique_reads_file in tmr1000_unique_reads_files
]

tmr1000_reads_data_df = pd.DataFrame(
    {
        "Chrom": tmr1000_chroms_in_reads_files,
        "ReadsFile": tmr1000_reads_files,
    }
)

tmr1000_unique_reads_data_df = pd.DataFrame(
    {
        "Chrom": tmr1000_chroms_in_unique_reads_files,
        "UniqueReadsFile": tmr1000_unique_reads_files,
    }
)

tmr1000_reads_data_df = tmr1000_reads_data_df.merge(
    tmr1000_unique_reads_data_df, on="Chrom", how="left"
)

tmr1000_reads_data_df

# %%
tmr1000_proteins_files = list(tmr1000_proteins_dir.glob("*.proteins.csv.gz"))
tmr1000_chroms_in_proteins_files = [
    proteins_file.name.split(".")[0] for proteins_file in tmr1000_proteins_files
]

tmr1000_unique_proteins_files = list(
    tmr1000_proteins_dir.glob("*.unique_proteins.csv.gz")
)
tmr1000_chroms_in_unique_proteins_files = [
    unique_proteins_file.name.split(".")[0]
    for unique_proteins_file in tmr1000_unique_proteins_files
]

# distinct_proteins_files = list(distinct_proteins_dir.glob("*DistinctUniqueProteins.*.csv"))
tmr1000_distinct_proteins_files = [
    f
    for f in tmr1000_distinct_proteins_dir.glob("*DistinctUniqueProteins.*.csv")
    if "expression" not in f.name.lower()
]
tmr1000_chroms_in_distinct_proteins_files = [
    distinct_proteins_file.name.split(".")[0]
    for distinct_proteins_file in tmr1000_distinct_proteins_files
]

# expression_files = list(
#     distinct_proteins_dir.glob("*.DistinctUniqueProteins.ExpressionLevels.csv")
# )
# chroms_in_expression_files = [
#     expression_file.name.split(".")[0] for expression_file in expression_files
# ]


tmr1000_proteins_data_df = pd.DataFrame(
    {
        "Chrom": tmr1000_chroms_in_proteins_files,
        "ProteinsFile": tmr1000_proteins_files,
    }
)

tmr1000_unique_proteins_data_df = pd.DataFrame(
    {
        "Chrom": tmr1000_chroms_in_unique_proteins_files,
        "UniqueProteinsFile": tmr1000_unique_proteins_files,
    }
)

tmr1000_distinct_proteins_data_df = pd.DataFrame(
    {
        "Chrom": tmr1000_chroms_in_distinct_proteins_files,
        "DistinctProteinsFile": tmr1000_distinct_proteins_files,
    }
)

# expression_data_df = pd.DataFrame(
#     {"Chrom": chroms_in_expression_files, "ExpressionFile": expression_files}
# )

tmr1000_proteins_data_df = (
    tmr1000_proteins_data_df.merge(
        tmr1000_unique_proteins_data_df, on="Chrom", how="left"
    ).merge(tmr1000_distinct_proteins_data_df, on="Chrom", how="left")
    # .merge(expression_data_df, on="Chrom", how="left")
)

tmr1000_proteins_data_df

# %%
tmr1000_data_df = (
    orfs_df.merge(alignment_stats_df, on="Chrom", how="left")
    .merge(tmr1000_positions_data_df, on="Chrom", how="left")
    .merge(tmr1000_reads_data_df, on="Chrom", how="left")
    .merge(tmr1000_proteins_data_df, on="Chrom", how="left")
)

tmr1000_data_df

# %%
# complete_data_df = data_df.loc[data_df["ExpressionFile"].notna()].reset_index(drop=True)
tmr1000_complete_data_df = tmr1000_data_df.loc[
    tmr1000_data_df["DistinctProteinsFile"].notna()
].reset_index(drop=True)

# complete_data_df = complete_data_df.drop_duplicates(
#     "Name", keep=False, ignore_index=True
# )
tmr1000_complete_data_df

# %%
# complete_data_df.loc[complete_data_df["Name"].duplicated(keep=False)].sort_values("Name", ignore_index=True)

# %%
tmr1000_possibly_na_conditions = tmr1000_data_df["Name"].tolist()
tmr1000_possibly_na_chroms = tmr1000_data_df["Chrom"].tolist()
tmr1000_possibly_na_starts = tmr1000_data_df["Start"].tolist()
tmr1000_possibly_na_ends = tmr1000_data_df["End"].tolist()
tmr1000_possibly_na_strands = tmr1000_data_df["Strand"].tolist()

tmr1000_possibly_na_positions_files = tmr1000_data_df["PositionsFile"].tolist()
tmr1000_possibly_na_reads_files = tmr1000_data_df["ReadsFile"].tolist()
tmr1000_possibly_na_unique_reads_files = tmr1000_data_df["UniqueReadsFile"].tolist()
tmr1000_possibly_na_proteins_files = tmr1000_data_df["ProteinsFile"].tolist()
tmr1000_possibly_na_unique_proteins_files = tmr1000_data_df["UniqueProteinsFile"].tolist()
tmr1000_possibly_na_distinct_unique_proteins_files = tmr1000_data_df["DistinctProteinsFile"].tolist()
# expression_files = complete_data_df["ExpressionFile"].tolist()

# %%
tmr1000_conditions = tmr1000_complete_data_df["Name"].tolist()
tmr1000_chroms = tmr1000_complete_data_df["Chrom"].tolist()
tmr1000_starts = tmr1000_complete_data_df["Start"].tolist()
tmr1000_ends = tmr1000_complete_data_df["End"].tolist()
tmr1000_strands = tmr1000_complete_data_df["Strand"].tolist()

tmr1000_positions_files = tmr1000_complete_data_df["PositionsFile"].tolist()
tmr1000_reads_files = tmr1000_complete_data_df["ReadsFile"].tolist()
tmr1000_unique_reads_files = tmr1000_complete_data_df["UniqueReadsFile"].tolist()
tmr1000_proteins_files = tmr1000_complete_data_df["ProteinsFile"].tolist()
tmr1000_unique_proteins_files = tmr1000_complete_data_df["UniqueProteinsFile"].tolist()
tmr1000_distinct_unique_proteins_files = complete_data_df[
    "DistinctProteinsFile"
].tolist()
# expression_files = complete_data_df["ExpressionFile"].tolist()

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
# # Data preprocessing

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
fig.update_layout(width=800, height=400)
fig.show()

# %% [markdown] papermill={"duration": 0.041741, "end_time": "2022-02-01T09:42:47.760215", "exception": false, "start_time": "2022-02-01T09:42:47.718474", "status": "completed"}
# ## Positions

# %% jupyter={"source_hidden": true}
# positions_dfs = [
#     pd.read_csv(position_file, sep=sep, dtype={"Reads": str}) for position_file in positions_files
# ]
# for positions_df, condition in zip(positions_dfs, conditions):
#     positions_df.insert(0, condition_col, condition)
# positions_dfs[0]


# %% jupyter={"source_hidden": true}
# positions_df = pd.concat(positions_dfs, ignore_index=True)
# # positions_df.insert(positions_df.columns.get_loc("G")+1, "ATCGs", positions_df.loc[:, ["A", "T", "C", "G"]].sum(axis=1))
# positions_df

# %% jupyter={"source_hidden": true}
# def make_concat_positions_df(positions_files, condition_col, conditions):
#     positions_dfs = [
#         pd.read_csv(position_file, sep=sep, dtype={"Reads": str}) for position_file in positions_files
#     ]
#     for positions_df, condition in zip(positions_dfs, conditions):
#         positions_df.insert(0, condition_col, condition)
#     concat_positions_df = pd.concat(positions_dfs, ignore_index=True)
#     return concat_positions_df

# %% jupyter={"source_hidden": true}
# concat_positions_df = make_concat_positions_df(positions_files, condition_col, conditions)
# concat_positions_df

# %% jupyter={"source_hidden": true}
# data_df.loc[data_df["PositionsFile"].isna()]

# %% jupyter={"source_hidden": true}
# possibly_na_positions_files[108]

# %% jupyter={"source_hidden": true}
# pd.notna(possibly_na_positions_files[108])

# %% jupyter={"source_hidden": true}
# all_positions_dfs = [
#     pd.read_csv(position_file, sep=sep, dtype={"Reads": str}) for position_file in possibly_na_positions_files if pd.notna(position_file)
# ]
# for positions_df, condition in zip(all_positions_dfs, possibly_na_conditions):
#     positions_df.insert(0, condition_col, condition)
# all_positions_dfs[0]

# %%
# all_positions_df = pd.concat(all_positions_dfs, ignore_index=True)
# # positions_df.insert(positions_df.columns.get_loc("G")+1, "ATCGs", positions_df.loc[:, ["A", "T", "C", "G"]].sum(axis=1))
# all_positions_df

# %%
def make_concat_all_positions_df(possibly_na_positions_files, condition_col, conditions):
    all_positions_dfs = [
        pd.read_csv(position_file, sep=sep, dtype={"Reads": str}) 
        for position_file in possibly_na_positions_files 
        if pd.notna(position_file)
    ]
    for positions_df, condition in zip(all_positions_dfs, possibly_na_conditions):
        positions_df.insert(0, condition_col, condition)
    concat_all_positions_df = pd.concat(all_positions_dfs, ignore_index=True)
    return concat_all_positions_df


# %%
concat_all_positions_df = make_concat_all_positions_df(possibly_na_positions_files, condition_col, conditions)
concat_all_positions_df

# %%
# # previoisly edited positions, which are now not considered edited
# positions_df.loc[
#     (positions_df["Edited"]) & (~positions_df["EditedFinal"])
# ]

# %%
concat_all_positions_df.loc[
    (concat_all_positions_df["NoisyCorrected"].fillna(False))
    & (concat_all_positions_df["Noise"] <= 0.1)
]

# %%
x = concat_all_positions_df.loc[
    (concat_all_positions_df["NoisyCorrected"].fillna(False)),
    "Noise"
] * 100
                                   
fig = go.Figure()
fig.add_trace(go.Histogram(x=x, cumulative_enabled=True, 
                           # histnorm='percent'
                          ))

fig.update_xaxes(title="Noise [%]")
fig.update_yaxes(type="log", title="Positions")
fig.update_layout(width=700, height=500, template=template)
fig.show()

# %%
concat_all_positions_df.loc[
    concat_all_positions_df["EditedFinal"],
    [
        condition_col,
        "Chrom",
        "Position",
        "RefBase",
        "TotalCoverage",
        "A",
        "T",
        "C",
        "G",
        "EditingFrequency",
    ],
]

# %%
concat_all_positions_df.loc[
    # all edited positions in transcripts whose pooled noise levels is < 6%
    (concat_all_positions_df["EditedFinal"]) & (concat_all_positions_df["Chrom"].isin(chroms)),
    [
        condition_col,
        "Chrom",
        "Position",
        "RefBase",
        "TotalCoverage",
        "A",
        "T",
        "C",
        "G",
        "EditingFrequency",
    ],
]

# %%
# edited_positions_df = positions_df.loc[
#     positions_df["EditedFinal"],
#     [
#         condition_col,
#         "Chrom",
#         "Position",
#         "RefBase",
#         "TotalCoverage",
#         "A",
#         "T",
#         "C",
#         "G",
#         "EditingFrequency",
#     ],
# ]
# edited_positions_df

# %%
# edited_positions_df.loc[edited_positions_df["RefBase"]=="A"]

# %%
# edited_positions_df["TotalCoverage"].describe()

# %%
# edited_positions_df["TotalCoverage"].value_counts()

# %%
# totcov_vs_editfreq_df = edited_positions_df.groupby(["TotalCoverage", "EditingFrequency"]).size().reset_index().rename(columns={0: "DenovoSites"})
# totcov_vs_editfreq_df

# %%
# totcov_vs_editing_df = edited_positions_df.groupby(["TotalCoverage", "G"]).size().reset_index().rename(columns={0: "DenovoSites"})
# totcov_vs_editing_df

# %%
# totcov_vs_editing_df.head(30)

# %%
# totcov_vs_editing_df.groupby("G").size().reset_index()

# %%
# totcov_vs_editing_df.loc[totcov_vs_editing_df["G"] == 1]

# %%
# fig = px.scatter(
#     totcov_vs_editfreq_df,
#     x="TotalCoverage",
#     y="EditingFrequency",
#     color="DenovoSites",
#     template=template,
#     log_x=True
# )
# fig.update_layout(
#     width=600,
#     height=400
# )
# fig.show()

# %%
# fig = px.scatter(
#     totcov_vs_editing_df,
#     x="TotalCoverage",
#     y="G",
#     color="DenovoSites",
#     template=template,
#     log_x=True,
#     marginal_x="rug",
#     marginal_y="histogram"
# )
# fig.update_layout(
#     width=600,
#     height=400
# )
# fig.show()

# %%
# fig = px.histogram(
#     edited_positions_df,
#     x="TotalCoverage",
#     # y="EditingFrequency",
#     # histfunc="avg",
#     marginal="box",
#     template=template,
#     # log_x=True,
#     log_y=True
# )
# fig.update_layout(
#     width=600,
#     height=400
# )
# fig.show()

# %%
# fig = px.histogram(
#     edited_positions_df,
#     x="TotalCoverage",
#     y="EditingFrequency",
#     histfunc="avg",
#     marginal="box",
#     template=template,
#     # log_x=True
# )
# fig.update_layout(
#     width=600,
#     height=400
# )
# fig.show()

# %%
# fig = px.histogram(
#     totcov_vs_editing_df,
#     x="TotalCoverage",
#     y="G",
#     histfunc="avg",
#     marginal="box",
#     template=template,
#     # log_x=True
# )
# fig.update_layout(
#     width=600,
#     height=400
# )
# fig.show()

# %%
# fig = px.histogram(
#     edited_positions_df.loc[edited_positions_df["TotalCoverage"]<=1000],
#     x="TotalCoverage",
#     y="EditingFrequency",
#     histfunc="avg",
#     marginal="box",
#     template=template
# )
# fig.update_layout(
#     width=600,
#     height=400
# )
# fig.show()

# %%

# %%

# %%

# %%

# %%
# editing_positions_per_sample = [len(df.loc[(df["EditedFinal"])]) for df in positions_dfs]
# print(
#     f"Average of {sum(editing_positions_per_sample)/len(positions_dfs)} editing sites per sample"
# )

# %%
# avg editing positions per transcript, considering transcripts whose pooled noise levels is < 6%
(
    concat_all_positions_df.loc[(concat_all_positions_df["EditedFinal"]) & (concat_all_positions_df["Chrom"].isin(chroms))]
    .groupby("Chrom")
    .size()
    .mean()
    .round(2)
    
)


# %% [markdown] papermill={"duration": 0.02598, "end_time": "2022-02-01T09:42:46.438342", "exception": false, "start_time": "2022-02-01T09:42:46.412362", "status": "completed"}
# ## Reads

# %% [markdown] papermill={"duration": 0.02598, "end_time": "2022-02-01T09:42:46.438342", "exception": false, "start_time": "2022-02-01T09:42:46.412362", "status": "completed"} jp-MarkdownHeadingCollapsed=true
# ### All

# %% [markdown]
# That is, all filtered reads.

# %% papermill={"duration": 1.204258, "end_time": "2022-02-01T09:42:47.668206", "exception": false, "start_time": "2022-02-01T09:42:46.463948", "status": "completed"}
# reads_dfs = [pd.read_csv(reads_file, sep=sep, dtype={"Read": str}) for reads_file in reads_files]
# reads_dfs[0]


# %%
# reads_dfs[1]

# %%
# len(reads_dfs[0])

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
    pd.read_csv(unique_reads_file, sep=sep, dtype={"UniqueRead": str, "Reads": str}) for unique_reads_file in unique_reads_files
]
for chrom, unique_reads_df in zip(chroms, unique_reads_dfs):
    unique_reads_df.insert(0, "Chrom", chrom)
unique_reads_dfs[0]


# %%
unique_reads_dfs[0].columns

# %%
# def make_concat_unique_reads_df(unique_reads_files, chroms):
#     unique_reads_dfs = [
#         pd.read_csv(unique_reads_file, sep=sep, dtype={"UniqueRead": str, "Reads": str}) 
#         for unique_reads_file in unique_reads_files
#     ]
#     for chrom, unique_reads_df in zip(chroms, unique_reads_dfs):
#         unique_reads_df.insert(0, "Chrom", chrom)
#     concat_unique_reads_df = pd.concat(unique_reads_dfs, ignore_index=True)
#     return concat_unique_reads_df

# %%
# concat_unique_reads_df = make_concat_unique_reads_df(unique_reads_files, chroms)
# concat_unique_reads_df

# %%
# expanded_unique_reads_dfs = []
# for unique_reads_df in unique_reads_dfs:
#     unique_reads_df = unique_reads_df.copy()
#     unique_reads_df["Samples"] = unique_reads_df["Samples"].str.split(",")
#     expanded_unique_reads_df = unique_reads_df.explode("Samples").reset_index(drop=True)
#     # after exploding the df by samples, it may be that some reads/unique reads only appear in certain samples,
#     # so in order to get that information, one would have to merge the `expanded_unique_proteins_df`
#     # with a corresponding `expanded_reads_df`/`expanded_unique_reads_df`
#     # expanded_unique_reads_df = expanded_unique_reads_df.drop(
#     #     [
#     #         "Reads",
#     #         "NumOfReads",
#     #         "UniqueReads",
#     #         "NumOfUniqueReads",
#     #         "EditingFrequency",
#     #         "EditedPositions",
#     #         "UneditedPositions",
#     #         "AmbigousPositions",
#     #     ],
#     #     axis=1,
#     # )
#     # expanded_unique_reads_df = expanded_unique_reads_df.rename(
#     #     columns={"Samples": "Sample"}
#     # )
#     expanded_unique_reads_dfs.append(expanded_unique_reads_df)
#     break

# expanded_unique_reads_dfs[0]

# # expanded_unique_proteins_df = pd.concat(expanded_unique_proteins_dfs, ignore_index=True)
# # del expanded_unique_proteins_dfs
# # expanded_unique_proteins_df

# %%
# expanded_unique_reads_dfs[0].columns

# %% [markdown] papermill={"duration": 0.041741, "end_time": "2022-02-01T09:42:47.760215", "exception": false, "start_time": "2022-02-01T09:42:47.718474", "status": "completed"}
# ### Unique - TMR 1000

# %% papermill={"duration": 0.126539, "end_time": "2022-02-01T09:42:47.923363", "exception": false, "start_time": "2022-02-01T09:42:47.796824", "status": "completed"}
tmr1000_unique_reads_dfs = [
    pd.read_csv(unique_reads_file, sep=sep, dtype={"UniqueRead": str, "Reads": str})
    for unique_reads_file in tmr1000_unique_reads_files
]
for chrom, unique_reads_df in zip(tmr1000_chroms, tmr1000_unique_reads_dfs):
    unique_reads_df.insert(0, "Chrom", chrom)
tmr1000_unique_reads_dfs[0]


# %%
len(tmr1000_unique_reads_dfs)

# %%
# tmr1000_concat_unique_reads_df = make_concat_unique_reads_df(tmr1000_unique_reads_files, tmr1000_chroms)
# tmr1000_concat_unique_reads_df

# %% [markdown]
# ## Proteins

# %% [markdown] jp-MarkdownHeadingCollapsed=true
# ### All proteins

# %%
# proteins_dfs = [pd.read_csv(proteins_file, sep=sep, dtype={"UniqueRead": str, "Reads": str}) for proteins_file in proteins_files]
# # for proteins_df in proteins_dfs:
# #     if "Transcript" in proteins_df.columns:
# #         proteins_df.rename(columns={"Transcript": "UniqueRead"}, inplace=True)
# proteins_dfs[0]


# %% [markdown]
# ### Unique proteins

# %%
unique_proteins_dfs = [
    pd.read_csv(unique_proteins_file, sep=sep, dtype={"Protein": str, "Reads": str})
    for unique_proteins_file in unique_proteins_files
]
for chrom, unique_proteins_df in zip(chroms, unique_proteins_dfs):
    unique_proteins_df.insert(0, "Chrom", chrom)
# the position of the first column needs to be updated due to the Chrom col insertion
unique_proteins_first_col_pos += 1
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
unique_proteins_dfs[0].columns

# %%
unique_proteins_dfs[0].iloc[:, unique_proteins_first_col_pos:]

# %%
# def make_concat_unique_proteins_df(unique_proteins_files, chroms):
#     unique_proteins_dfs = [
#         pd.read_csv(unique_proteins_file, sep=sep, dtype={"Protein": str, "Reads": str})
#         for unique_proteins_file in unique_proteins_files
#     ]
#     for chrom, unique_proteins_df in zip(chroms, unique_proteins_dfs):
#         unique_proteins_df.insert(0, "Chrom", chrom)
#     unique_proteins_first_col_pos += 1
#     concat_unique_proteins_df = pd.concat(unique_proteins_dfs, ignore_index=True)
#     return concat_unique_proteins_df

# %%
# concat_unique_proteins_df = make_concat_unique_proteins_df(unique_proteins_files, chroms)
# concat_unique_proteins_df

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
# statistics for num of editable AAs per transcript

editable_aas_per_sample = pd.Series([
    df.iloc[:, unique_proteins_first_col_pos:].shape[1] for df in unique_proteins_dfs
])

editable_aas_per_sample.describe()

# %%
unique_proteins_dfs[0].iloc[:, unique_proteins_first_col_pos:]

# %%
# pd.DataFrame(
#     {
#         condition_col: conditions,
#         "EditableAAs": [
#             unique_proteins_df.iloc[:, unique_proteins_first_col_pos:].shape[1]
#             for unique_proteins_df in unique_proteins_dfs
#         ],
#     }
# )

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
assert len(conditions) == len(chroms) == len(distinct_unique_proteins_files) == len(unique_reads_dfs)

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
# ### Distinct unique proteins - TMR 1000

# %%
assert len(tmr1000_conditions) == len(tmr1000_chroms) == len(tmr1000_distinct_proteins_files) == len(tmr1000_unique_reads_dfs)

tmr1000_distinct_unique_proteins_dfs = []
for condition, chrom, distinct_unique_proteins_file, unique_reads_df in zip(
    tmr1000_conditions,
    tmr1000_chroms,
    tmr1000_distinct_proteins_files,
    tmr1000_unique_reads_dfs,
):
    tmr1000_distinct_unique_proteins_df = pd.read_csv(
        distinct_unique_proteins_file, sep=sep
    )
    tmr1000_distinct_unique_proteins_df.insert(0, condition_col, condition)
    tmr1000_distinct_unique_proteins_df.insert(
        1,
        "NumOfReads",
        (
            tmr1000_distinct_unique_proteins_df["Fraction"]
            * unique_reads_df["NumOfReads"].sum()
        ).astype(int),
    )
    tmr1000_distinct_unique_proteins_df.insert(0, "Chrom", chrom)
    tmr1000_distinct_unique_proteins_dfs.append(tmr1000_distinct_unique_proteins_df)

ic(len(tmr1000_distinct_unique_proteins_dfs))

tmr1000_distinct_unique_proteins_df = (
    pd.concat(tmr1000_distinct_unique_proteins_dfs)
    .reset_index(drop=True)
    .rename(columns={"NumUniqueSamples": "NumOfProteins", "UniqueSamples": "Proteins"})
)

tmr1000_distinct_unique_proteins_df = tmr1000_distinct_unique_proteins_df.sort_values(
    [
        condition_col,
        "Fraction",
        "FractionRepetition",
        "Algorithm",
        "AlgorithmRepetition",
    ]
).reset_index(drop=True)

tmr1000_distinct_unique_proteins_df


# %% [markdown] papermill={"duration": 0.045853, "end_time": "2022-02-01T09:42:48.953594", "exception": false, "start_time": "2022-02-01T09:42:48.907741", "status": "completed"}
# # Results

# %% [markdown] papermill={"duration": 0.124528, "end_time": "2022-02-01T09:43:10.054394", "exception": false, "start_time": "2022-02-01T09:43:09.929866", "status": "completed"}
# ## Positions

# %% [markdown] papermill={"duration": 0.149848, "end_time": "2022-02-01T09:43:12.800733", "exception": false, "start_time": "2022-02-01T09:43:12.650885", "status": "completed"}
# ### Coverage - per-transcript per-sample

# %%
concat_all_positions_df


# %%
def calc_per_transcript_per_sample_coverage(
    positions_df,
    samples_and_tissues_df
    # samples
):
    expanded_positions_df = (
        positions_df.loc[(~positions_df["InProbRegion"]) & (positions_df["CDS"])]
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
# # per_transcript_per_sample_coverage_dfs = [
# #     calc_per_transcript_per_sample_coverage(positions_df)
# #     for positions_df in positions_dfs
# # ]

# with Pool(processes=4) as pool:
#     per_transcript_per_sample_coverage_dfs = pool.starmap(
#         func=calc_per_transcript_per_sample_coverage,
#         iterable=[
#             (positions_df, samples_and_tissues_df) for positions_df in positions_dfs
#         ],
#     )
# per_transcript_per_sample_coverage_dfs[0]

# %%
def calc_per_transcript_per_sample_coverage_dfs(concat_all_positions_df, possibly_na_positions_files, possibly_na_chroms, samples_and_tissues_df, processes=4):
    positions_dfs = [
        concat_all_positions_df.loc[concat_all_positions_df["Chrom"] == chrom]
        for position_file, chrom in zip(possibly_na_positions_files, possibly_na_chroms)
        if pd.notna(position_file)
        
    ]
    with Pool(processes=processes) as pool:
        per_transcript_per_sample_coverage_dfs = pool.starmap(
            func=calc_per_transcript_per_sample_coverage,
            iterable=[
                (positions_df, samples_and_tissues_df) for positions_df in positions_dfs
            ],
        )
    return per_transcript_per_sample_coverage_dfs
    


# %%
per_transcript_per_sample_coverage_dfs = calc_per_transcript_per_sample_coverage_dfs(
    concat_all_positions_df, possibly_na_positions_files, possibly_na_chroms, samples_and_tissues_df
)
per_transcript_per_sample_coverage_dfs[0]

# %%
len(per_transcript_per_sample_coverage_dfs)

# %%
tmr50_alignment_stats_df.shape

# %%
# merged_per_transcript_per_sample_coverage_df = pd.concat(
#     per_transcript_per_sample_coverage_dfs
# ).reset_index()
# # merged_per_transcript_per_sample_coverage_df = merged_per_transcript_per_sample_coverage_df.merge(samples_and_tissues_df, how="left")
# merged_per_transcript_per_sample_coverage_df

# %%
fig = px.histogram(
    pd.concat(per_transcript_per_sample_coverage_dfs),
    x="NumOfReads",
    color="Tissue",
    color_discrete_map=tissues_color_discrete_map,
    facet_col="Tissue",
    facet_col_wrap=3,
    facet_col_spacing=facet_col_spacing,
    facet_row_spacing=facet_row_spacing * 0.5,
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

fig.update_layout(template=template, showlegend=False, width=width, height=height)
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
        (positions_df["RefBase"] == ref_base) & (~positions_df["InProbRegion"]) & (positions_df["CDS"])
    ]
    num_of_all_editable_adenosines = all_refbase_positions_df["TotalCoverage"].sum()
    num_of_edited_adenosines = all_refbase_positions_df[alt_base].sum()
    editing_index = 100 * num_of_edited_adenosines / num_of_all_editable_adenosines
    return editing_index


def make_one_input(concat_all_positions_df, chrom, strand): 
        return concat_all_positions_df.loc[concat_all_positions_df["Chrom"] == chrom], strand

    
def make_editing_index_per_transcript_inputs(concat_all_positions_df, chroms, strands):
    
    concat_all_positions_df = concat_all_positions_df.loc[concat_all_positions_df["Chrom"].isin(chroms)]
    
    editing_index_per_transcript_inputs = [
        (concat_all_positions_df.loc[concat_all_positions_df["Chrom"] == chrom], strand)
        for chrom, strand in zip(chroms, strands)
    ]
    
    # editing_index_per_transcript_inputs = [
    #     x
    #     for x in editing_index_per_transcript_inputs
    #     if pd.notna(x)
    # ]
    
    return editing_index_per_transcript_inputs


def calc_per_transcript_editing_index_df(concat_all_positions_df, chroms, strands, processes=4):
    
#     # todo comment out - this is for testing
#     chroms, strands = chroms[:100], strands[:100]
    
    editing_index_per_transcript_inputs = make_editing_index_per_transcript_inputs(
        concat_all_positions_df, chroms, strands
    )
    
    with Pool(processes=processes) as pool:
        per_transcript_editing_indices = pool.starmap(
            func=editing_index_per_transcript,
            iterable=editing_index_per_transcript_inputs,
        )
    per_transcript_editing_index_df = pd.DataFrame(
        {
            "Chrom": chroms,
            "EditingIndex": per_transcript_editing_indices
        }
    )
    
    return per_transcript_editing_index_df


# %%
# %%time

per_transcript_editing_index_df = calc_per_transcript_editing_index_df(
    concat_all_positions_df, chroms, strands, 12
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
            (positions_df["RefBase"] == ref_base) & (~positions_df["InProbRegion"]) & (positions_df["CDS"])
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
def calc_per_sample_editing_index_df(concat_all_positions_df, chroms, strands, samples, processes=1):
    
    concat_all_positions_df = concat_all_positions_df.loc[concat_all_positions_df["Chrom"].isin(chroms)]
    
    positions_dfs = [
        concat_all_positions_df.loc[concat_all_positions_df["Chrom"] == chrom]
        for chrom in chroms
    ]
    
    all_per_sample_a_and_g_counts = calc_all_per_sample_a_and_g_counts(
        positions_dfs, strands, samples, processes
    )
    per_sample_editing_index = [
        100 * g_count / (a_count + g_count)
        for a_count, g_count in all_per_sample_a_and_g_counts
    ]
    
    per_sample_editing_index_df = pd.DataFrame(
        {"Sample": samples, "EditingIndex": per_sample_editing_index}
    )
    
    return per_sample_editing_index_df


# %%
# %%time 

per_sample_editing_index_df = calc_per_sample_editing_index_df(
    concat_all_positions_df, chroms, strands, samples, processes=8
)
per_sample_editing_index_df

# %% [markdown]
# ### Noise in positions

# %%
(
    concat_all_positions_df.loc[
        (concat_all_positions_df["Noise"] <= 0.1)
            & (concat_all_positions_df["NoisyCorrected"]),
    ]
    .groupby("Chrom")
    .size()
    .reset_index()
    .rename(columns={0: "NoisePositions"})
    ["NoisePositions"].describe()
)


# %%
def mean_noise_levels(positions_df, top_x_noisy_positions=3, snp_noise_level=0.1):
    # if positions_df.empty:
    #     return 0.0
    noise_levels = (
        positions_df.loc[
            (positions_df["Noise"] <= snp_noise_level)
            & (positions_df["NoisyCorrected"]),
            "Noise",
        ].sort_values(ascending=False)[:top_x_noisy_positions].tolist()
    )
    # if there are less noisy positions than `top_x_noisy_positions`, add zeros accordingly
    noise_levels = pd.Series(
        noise_levels + [0 for _ in range(top_x_noisy_positions - len(noise_levels))]
    )
    return noise_levels.mean()


# %%
tmr50_alignment_stats_df

# %%
all_per_chrom_mean_noise_levels =(
    concat_all_positions_df
    # .groupby(["Transcript", "Chrom"])
    .groupby("Chrom")
    .apply(mean_noise_levels)
    .reset_index()
    .rename(columns={0: "Noise"})
    .merge(
        tmr50_alignment_stats_df.loc[:, ["Chrom"]],
        on="Chrom",
        how="right"
    )
    .fillna(0.0)
    .sort_values(["Chrom", "Noise"])
)
all_per_chrom_mean_noise_levels["%Noise"] = 100 * all_per_chrom_mean_noise_levels["Noise"]
all_per_chrom_mean_noise_levels

# %%
saved_per_chrom_mean_noise_levels_df = all_per_chrom_mean_noise_levels.merge(orfs_df[["Chrom", "Name"]], how="left")
saved_per_chrom_mean_noise_levels_df.insert(1, "Gene", saved_per_chrom_mean_noise_levels_df["Name"])
del saved_per_chrom_mean_noise_levels_df["Name"]
saved_per_chrom_mean_noise_levels_df.insert(0, "Platform", "Whole-transcriptome octopus data")
saved_per_chrom_mean_noise_levels_df.to_csv("Noise.Octopus.tsv", sep="\t", index=False)
saved_per_chrom_mean_noise_levels_df

# %%
# per_chrom_mean_noise_levels["%Noise"].describe()

# %%
all_per_chrom_mean_noise_levels["%Noise"].describe()

# %%
# described_noise_df = per_chrom_mean_noise_levels["%Noise"].describe()
# quartiles = ["25%", "50%", "75%"]
# noise_quartiles = described_noise_df.loc[quartiles].values
# noise_quartiles

# %%
all_described_noise_df = all_per_chrom_mean_noise_levels["%Noise"].describe()
quartiles = ["25%", "50%", "75%"]
all_noise_quartiles = all_described_noise_df.loc[quartiles].values
all_noise_quartiles

# %%
# new version with the noise quartiles

fig = px.histogram(
    all_per_chrom_mean_noise_levels,
    x="%Noise",
    color_discrete_sequence=["black"],
    labels={"%Noise": "Per-gene noise level [%]"},
    log_y=True,
    opacity=0.4
)

f = fig.full_figure_for_development(warn=False)
x = f.data[0]["x"]
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
max_noise_quartiles_y = max_count * 0.8

for i, (quartile, noise_quartile) in enumerate(zip(quartiles, all_noise_quartiles)):
    try:
        # no point in plotting two lines very close to each other - the next one should suffice
        ic(all_noise_quartiles[i], all_noise_quartiles[i+1])
        if np.isclose(all_noise_quartiles[i], all_noise_quartiles[i+1]):
            continue
    except IndexError:
        # if this is the last line to plot
        pass
    fig.add_shape(
        type="line",
        x0=noise_quartile,
        x1=noise_quartile,
        y0=0,
        y1=max_noise_quartiles_y - (i * 0.3 * max_noise_quartiles_y),
        line=dict(
            color="red",
            width=5,
            dash="dash",
        ),
        opacity=0.5,
        label=dict(
            text=f"{quartile}<br>     of genes",
            # text=f"{quartile} of genes",
            textposition="end",
            textangle=45,
            # font_size=14,
        ),
    )

# fig.update_xaxes(dtick=1)
# fig.update_yaxes(dtick=50)


width = 650
height = 650 * 400 / 560

fig.update_layout(
    #  xaxis_title="Editing frequency",
    # title="Octopus",
    # title="Pooled octopus data",
    title="Noise level<br><sub>Whole-transcriptome long-reads (octopus)</sub>",
    title_x=0.15,
    yaxis_title="Genes",
    template=template,
    width=width,
    height=height,
    #  showlegend=False
)


fig.write_image(
    "Mean per chrom noise levels - with quartiles - Octopus.svg",
    width=width,
    height=height,
)

fig.show()

# %% [markdown]
# ### Machine noise

# %%
all_machine_noise_df = concat_all_positions_df.loc[:, [condition_col, "Chrom", "RefBase", "TotalCoverage", "A", "T", "C", "G"]]
all_machine_noise_df["ATCGs"] = all_machine_noise_df.loc[:, ["A", "T", "C", "G"]].sum(axis=1)
all_machine_noise_df["Matches"] = all_machine_noise_df.apply(
    lambda x: x[x["RefBase"]],
    axis=1
)
all_machine_noise_df["Mismatches"] = all_machine_noise_df.apply(
    lambda x: x["ATCGs"] - x["Matches"],
    axis=1
)
all_machine_noise_df

# %%
all_pooled_per_chrom_machine_noise_df = all_machine_noise_df.groupby("Chrom")[["Matches", "Mismatches"]].sum().reset_index()
all_pooled_per_chrom_machine_noise_df["%PooledMachineNoise"] = all_pooled_per_chrom_machine_noise_df.apply(lambda x: 100 * x["Mismatches"] / x["Matches"], axis=1)
all_pooled_per_chrom_machine_noise_df = all_pooled_per_chrom_machine_noise_df.merge(
        tmr50_alignment_stats_df.loc[:, ["Chrom"]],
        on="Chrom",
        how="right"
    ).fillna(0.0)
all_pooled_per_chrom_machine_noise_df

# %%
fig = px.histogram(
    all_pooled_per_chrom_machine_noise_df,
    x="%PooledMachineNoise",
    color_discrete_sequence=["black"],
    # labels={"% noise": "Per-gene noise level [%]"},
    log_y=True
)

# fig.update_xaxes(dtick=1)
# fig.update_yaxes(dtick=50)


width = 500
height = 650 * 400 / 560

fig.update_layout(
    #  xaxis_title="Editing frequency",
    # title="Octopus",
    # title="Pooled octopus data",
    # title="Noise level<br><sub>Whole-transcriptome long-reads (octopus)</sub>",
    title_x=0.15,
    yaxis_title="Genes",
    template=template,
    width=width,
    height=height,
    #  showlegend=False
)


# fig.write_image(
#     "Mean per chrom noise levels - with quartiles - Octopus.svg",
#     width=width,
#     height=height,
# )

fig.show()

# %% [markdown]
# ### Known & new editing sites

# %%
cols = 1
rows = 1

fig, ax = plt.subplots(
    # nrows=rows,
    # ncols=cols,
    # figsize=(3.5 * cols, 2.5 * rows),
    figsize=(4.5 * cols, 2.5 * rows),
    constrained_layout=True,
    gridspec_kw=dict(hspace=0.2, wspace=0.3),
)

labels = ["EditedFinal", "KnownEditing"]

sets = [
    set(
        concat_all_positions_df.loc[
            (concat_all_positions_df[label]) 
            & (concat_all_positions_df["CDS"]) 
            # only transcripts whose pooled noise levels is < 6%
            & (concat_all_positions_df["Chrom"].isin(chroms)), 
            [condition_col, "Chrom", "Position"]
        ].itertuples(index=False)
    )
    for label in labels
]
# labels[0] = f"Edited\n({len(sets[0])})"
# labels[1] = f"Known editing\n({len(sets[1])})"
labels[0] = f"De-novo\n({len(sets[0])})"
labels[1] = f"Known\n({len(sets[1])})"

venn2(sets, set_labels=labels, ax=ax)

fig.suptitle(
    # "Pooled octopus data",
    "Whole-transcriptome octopus data",
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
# Finding the number of genes/transcripts with no known editing sites

# %%
all_chroms = concat_all_positions_df["Chrom"].unique().size
chroms_with_new_sites = (
    concat_all_positions_df.loc[concat_all_positions_df["EditedFinal"] & ~concat_all_positions_df["KnownEditing"], "Chrom"]
    .unique()
    .size
)
chroms_without_new_sites = all_chroms - chroms_with_new_sites
chroms_without_new_sites

# %%
f"{100 * chroms_without_new_sites / all_chroms:.2f}%"

# %% [markdown]
# ### ADAR motif

# %%
unique_positions_df = (
    concat_all_positions_df.loc[
        # only edited positions in transcripts whose pooled noise levels is < 6%
        (concat_all_positions_df["EditedFinal"]) & (concat_all_positions_df["Chrom"].isin(chroms)),
        ["Chrom", "Position"]
    ]
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
# main_title = None
# sub_titles = ["Pooled octopus sites"]
main_title = "Whole-transcriptome octopus data"
sub_titles = [""]

# %%
out_file = Path("ADAR motif of pooled editing sites - Octopus.svg")

# %%
multiple_logos_from_fasta_files(
    fasta_files,
    main_title,
    sub_titles,
    out_file,
    width=0.33*14,
    height=4,
    dpi=300
);

# %% [markdown] papermill={"duration": 0.030615, "end_time": "2022-02-01T09:42:49.024262", "exception": false, "start_time": "2022-02-01T09:42:48.993647", "status": "completed"}
# ## Num of distinct proteins

# %% [markdown]
# ### Pooled

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

# %% jupyter={"source_hidden": true}
# max_distinct_proteins_df = (
#     distinct_unique_proteins_df.sort_values("Fraction", ascending=False)
#     .groupby("Chrom")
#     .apply(pd.DataFrame.nlargest, n=1, columns="NumOfProteins")
# )
# max_distinct_proteins_df = (
#     max_distinct_proteins_df.drop("Chrom", axis=1).reset_index().drop("level_1", axis=1)
# )

# # # max_distinct_proteins_df[condition_col] = max_distinct_proteins_df[
# # #     condition_col
# # # ].astype(str)

# max_distinct_proteins_df = max_distinct_proteins_df.merge(
#     alignment_stats_df,
#     on="Chrom",
#     # how="left",
#     how="right",
# )

# # max_distinct_proteins_df["NumOfProteins"] = max_distinct_proteins_df[
# #     "NumOfProteins"
# # ].fillna(1)

# # max_distinct_proteins_df["NumOfReads"] = max_distinct_proteins_df.apply(
# #     lambda x: x["NumOfReads"] if not pd.isna(x["NumOfReads"]) else x["MappedReads"],
# #     axis=1,
# # )

# max_distinct_proteins_df = max_distinct_proteins_df.dropna().reset_index(drop=True)

# max_distinct_proteins_df["DistinctProteins/Reads"] = (
#     max_distinct_proteins_df["NumOfProteins"] / max_distinct_proteins_df["NumOfReads"]
# )

# max_distinct_proteins_df = max_distinct_proteins_df.merge(
#     neural_vs_non_neural_expression_df.loc[:, ["OvulChrom", "IsNeural"]].rename(
#         columns={"OvulChrom": "Chrom"}
#     ),
#     on="Chrom",
#     how="left",
# )

# max_distinct_proteins_df["IsNeural"] = max_distinct_proteins_df["IsNeural"].fillna(
#     "Missing"
# )

# max_distinct_proteins_df = max_distinct_proteins_df.sort_values(
#     "NumOfProteins", ascending=False, ignore_index=True
# )
# # max_distinct_proteins_df["CummulativeTranscripts"] = 100 * (max_distinct_proteins_df.index + 1) / len(max_distinct_proteins_df)
# # max_distinct_proteins_df["CummulativeTranscripts"] = max_distinct_proteins_df["CummulativeTranscripts"][::-1].values

# max_distinct_proteins_df

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
    tmr50_alignment_stats_df,
    on="Chrom",
    # how="left",
    how="right",
)

max_distinct_proteins_df["NumOfProteins"] = max_distinct_proteins_df[
    "NumOfProteins"
].fillna(1)

max_distinct_proteins_df["NumOfReads"] = max_distinct_proteins_df.apply(
    lambda x: x["NumOfReads"] if not pd.isna(x["NumOfReads"]) else x["MappedReads"],
    axis=1,
)

# max_distinct_proteins_df = max_distinct_proteins_df.dropna().reset_index(drop=True)

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

max_distinct_proteins_df = max_distinct_proteins_df.sort_values(
    "NumOfProteins", ascending=False, ignore_index=True
)
# max_distinct_proteins_df["CummulativeTranscripts"] = 100 * (max_distinct_proteins_df.index + 1) / len(max_distinct_proteins_df)
# max_distinct_proteins_df["CummulativeTranscripts"] = max_distinct_proteins_df["CummulativeTranscripts"][::-1].values

max_distinct_proteins_df

# %%
max_distinct_proteins_df.loc[max_distinct_proteins_df["Chrom"] == robo2_chrom]

# %%
max_distinct_proteins_df["NumOfProteins"].mean()

# %%
# number of chroms with only 1 distinct protein
max_distinct_proteins_df.loc[max_distinct_proteins_df["NumOfProteins"] == 1].shape[0]

# %%
# % of chroms with only 1 distinct protein
100 * max_distinct_proteins_df.loc[
    max_distinct_proteins_df["NumOfProteins"] == 1
].shape[0] / max_distinct_proteins_df.shape[0]

# %%
# number of chroms with at least 5 distinct proteins
max_distinct_proteins_df.loc[max_distinct_proteins_df["NumOfProteins"] >= 5].shape[0]

# %%
# % of chroms with at least 5 distinct proteins
100 * max_distinct_proteins_df.loc[
    max_distinct_proteins_df["NumOfProteins"] >= 5
].shape[0] / max_distinct_proteins_df.shape[0]

# %%
# number of chroms with at least 50 distinct proteins
max_distinct_proteins_df.loc[max_distinct_proteins_df["NumOfProteins"] >= 50].shape[0]

# %%
# % of chroms with at least 50 distinct proteins
100 * max_distinct_proteins_df.loc[
    max_distinct_proteins_df["NumOfProteins"] >= 50
].shape[0] / max_distinct_proteins_df.shape[0]

# %%
max_distinct_proteins_df["IsNeural"].value_counts()

# %%
max_distinct_proteins_df["IsNeural"].value_counts(normalize=True).mul(100).round(2)

# %%
fig = px.histogram(max_distinct_proteins_df, x="IsNeural", log_y=True)
fig.update_layout(width=600, height=400, template=template)
fig.show()

# %% jupyter={"source_hidden": true}
# tmr1000_max_distinct_proteins_df = (
#     tmr1000_distinct_unique_proteins_df.sort_values("Fraction", ascending=False)
#     .groupby("Chrom")
#     .apply(pd.DataFrame.nlargest, n=1, columns="NumOfProteins")
# )
# tmr1000_max_distinct_proteins_df = (
#     tmr1000_max_distinct_proteins_df.drop("Chrom", axis=1)
#     .reset_index()
#     .drop("level_1", axis=1)
# )

# # # max_distinct_proteins_df[condition_col] = max_distinct_proteins_df[
# # #     condition_col
# # # ].astype(str)

# tmr1000_max_distinct_proteins_df = tmr1000_max_distinct_proteins_df.merge(
#     alignment_stats_df,
#     on="Chrom",
#     # how="left",
#     how="right",
# )

# # max_distinct_proteins_df["NumOfProteins"] = max_distinct_proteins_df[
# #     "NumOfProteins"
# # ].fillna(1)

# # max_distinct_proteins_df["NumOfReads"] = max_distinct_proteins_df.apply(
# #     lambda x: x["NumOfReads"] if not pd.isna(x["NumOfReads"]) else x["MappedReads"],
# #     axis=1,
# # )

# tmr1000_max_distinct_proteins_df = (
#     tmr1000_max_distinct_proteins_df.dropna().reset_index(drop=True)
# )

# tmr1000_max_distinct_proteins_df["DistinctProteins/Reads"] = (
#     tmr1000_max_distinct_proteins_df["NumOfProteins"]
#     / tmr1000_max_distinct_proteins_df["NumOfReads"]
# )

# tmr1000_max_distinct_proteins_df = tmr1000_max_distinct_proteins_df.merge(
#     neural_vs_non_neural_expression_df.loc[:, ["OvulChrom", "IsNeural"]].rename(
#         columns={"OvulChrom": "Chrom"}
#     ),
#     on="Chrom",
#     how="left",
# )

# tmr1000_max_distinct_proteins_df["IsNeural"] = tmr1000_max_distinct_proteins_df[
#     "IsNeural"
# ].fillna("Missing")

# tmr1000_max_distinct_proteins_df = tmr1000_max_distinct_proteins_df.sort_values(
#     "NumOfProteins", ascending=False, ignore_index=True
# )

# tmr1000_max_distinct_proteins_df

# %%
tmr1000_max_distinct_proteins_df = (
    tmr1000_distinct_unique_proteins_df.sort_values("Fraction", ascending=False)
    .groupby("Chrom")
    .apply(pd.DataFrame.nlargest, n=1, columns="NumOfProteins")
)
tmr1000_max_distinct_proteins_df = (
    tmr1000_max_distinct_proteins_df.drop("Chrom", axis=1)
    .reset_index()
    .drop("level_1", axis=1)
)

# # max_distinct_proteins_df[condition_col] = max_distinct_proteins_df[
# #     condition_col
# # ].astype(str)

tmr1000_max_distinct_proteins_df = tmr1000_max_distinct_proteins_df.merge(
    tmr50_alignment_stats_df,
    on="Chrom",
    # how="left",
    how="right",
)

max_distinct_proteins_df["NumOfProteins"] = max_distinct_proteins_df[
    "NumOfProteins"
].fillna(1)

max_distinct_proteins_df["NumOfReads"] = max_distinct_proteins_df.apply(
    lambda x: x["NumOfReads"] if not pd.isna(x["NumOfReads"]) else x["MappedReads"],
    axis=1,
)

tmr1000_max_distinct_proteins_df = (
    tmr1000_max_distinct_proteins_df.dropna().reset_index(drop=True)
)

tmr1000_max_distinct_proteins_df["DistinctProteins/Reads"] = (
    tmr1000_max_distinct_proteins_df["NumOfProteins"]
    / tmr1000_max_distinct_proteins_df["NumOfReads"]
)

tmr1000_max_distinct_proteins_df = tmr1000_max_distinct_proteins_df.merge(
    neural_vs_non_neural_expression_df.loc[:, ["OvulChrom", "IsNeural"]].rename(
        columns={"OvulChrom": "Chrom"}
    ),
    on="Chrom",
    how="left",
)

tmr1000_max_distinct_proteins_df["IsNeural"] = tmr1000_max_distinct_proteins_df[
    "IsNeural"
].fillna("Missing")

tmr1000_max_distinct_proteins_df = tmr1000_max_distinct_proteins_df.sort_values(
    "NumOfProteins", ascending=False, ignore_index=True
)

tmr1000_max_distinct_proteins_df

# %%
# number of chroms with at least 5 distinct proteins
tmr1000_max_distinct_proteins_df.loc[
    tmr1000_max_distinct_proteins_df["NumOfProteins"] >= 5
].shape[0]

# %%
# % of chroms with at least 5 distinct proteins
100 * tmr1000_max_distinct_proteins_df.loc[
    tmr1000_max_distinct_proteins_df["NumOfProteins"] >= 5
].shape[0] / tmr1000_max_distinct_proteins_df.shape[0]

# %%
# number of chroms with at least 50 distinct proteins
tmr1000_max_distinct_proteins_df.loc[
    tmr1000_max_distinct_proteins_df["NumOfProteins"] >= 50
].shape[0]

# %%
# % of chroms with at least 50 distinct proteins
100 * tmr1000_max_distinct_proteins_df.loc[
    tmr1000_max_distinct_proteins_df["NumOfProteins"] >= 50
].shape[0] / tmr1000_max_distinct_proteins_df.shape[0]

# %%
df = (
    max_distinct_proteins_df.loc[:, ["NumOfProteins"]]
    .sort_values("NumOfProteins")
    .reset_index(drop=True)
)
df["CummulativeTranscripts"] = 100 * (df.index + 1) / len(df)
df["CummulativeTranscripts"] = df["CummulativeTranscripts"][::-1].values
x = df["NumOfProteins"]
y = df["CummulativeTranscripts"]
x_mean = x.mean()
x_mean_closest = x.iloc[(x - x_mean).abs().argsort()[:1]]
x_mean_closest_k = x_mean_closest.index.values[0]
if x_mean == x_mean_closest.values[0]:
    y_mean = y.iloc[x_mean_closest_k]
else:
    if x_mean < x_mean_closest.values[0]:
        i = x_mean_closest_k - 1
        j = x_mean_closest_k + 1
    else:
        i = x_mean_closest_k
        j = x_mean_closest_k + 2
    y_mean = np.interp(x_mean, x.iloc[i:j], y.iloc[i:j])
df = df.drop_duplicates(subset="NumOfProteins").reset_index(drop=True)

tmr1000_df = (
    tmr1000_max_distinct_proteins_df.loc[:, ["NumOfProteins"]]
    .sort_values("NumOfProteins")
    .reset_index(drop=True)
)
tmr1000_df["CummulativeTranscripts"] = 100 * (tmr1000_df.index + 1) / len(tmr1000_df)
tmr1000_df["CummulativeTranscripts"] = tmr1000_df["CummulativeTranscripts"][::-1].values
tmr1000_x = tmr1000_df["NumOfProteins"]
tmr1000_y = tmr1000_df["CummulativeTranscripts"]
tmr1000_x_mean = tmr1000_x.mean()
tmr1000_x_mean_closest = tmr1000_x.iloc[
    (tmr1000_x - tmr1000_x_mean).abs().argsort()[:1]
]
tmr1000_x_mean_closest_k = tmr1000_x_mean_closest.index.values[0]
if tmr1000_x_mean == tmr1000_x_mean_closest.values[0]:
    tmr1000_y_mean = tmr1000_y.iloc[tmr1000_x_mean_closest_k]
else:
    if tmr1000_x_mean < tmr1000_x_mean_closest.values[0]:
        i = tmr1000_x_mean_closest_k - 1
        j = tmr1000_x_mean_closest_k + 1
    else:
        i = tmr1000_x_mean_closest_k
        j = tmr1000_x_mean_closest_k + 2
    tmr1000_y_mean = np.interp(tmr1000_x_mean, tmr1000_x.iloc[i:j], tmr1000_y.iloc[i:j])
tmr1000_df = tmr1000_df.drop_duplicates(subset="NumOfProteins").reset_index(drop=True)

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

tmr50_legendtitle = "50 reads"
tmr1000_legendtitle = "1000 reads"
legend_title_text = "Minimum coverage per gene      "

marker_size = 4

fig = make_subplots(
    rows=2,
    cols=1,
    # x_title="Distinct proteins per gene",
    x_title="Distinct protein isoforms per gene",
    y_title="% of genes",
    shared_yaxes=True,
    shared_xaxes=True,
    # vertical_spacing=facet_row_spacing / 2.5,
    # horizontal_spacing=facet_col_spacing * 1.5,
    vertical_spacing=0.05,
    # horizontal_spacing=0.025,
)

# tmr50

x = df["NumOfProteins"]
y = df["CummulativeTranscripts"]

y_min = min(y_min, y.min())

tmr50_all_color = "purple"

fig.add_trace(
    go.Scatter(
        x=x,
        y=y,
        mode="lines+markers",
        marker=dict(color=tmr50_all_color, size=marker_size),
        line=dict(color=tmr50_all_color, dash="dash"),
        # name="All",
        # legendgroup=tmr50_legendtitle,  # this can be any string
        # legendgrouptitle_text=tmr50_legendtitle,
        name=">= 50 reads",
        legend="legend",
    ),
    row=1,
    col=1,
)


fig.add_trace(
    go.Scatter(
        x=[x_mean],
        y=[y_mean],
        mode="markers+text",
        marker=dict(
            color=tmr50_all_color,
            size=marker_size * 2.5,
            # line=dict(
            #     color="yellow",
            #     width=3
            # )
        ),
        showlegend=False,
        # text=f"{x_mean:.0f} distinct proteins<br>(avg)",
        text=f"{x_mean:.0f} distinct<br>proteins<br>(avg)",
        textposition="bottom left",
        textfont=dict(color=tmr50_all_color, size=11),
    ),
    row=1,
    col=1,
)

# tmr1000

x = tmr1000_df["NumOfProteins"]
y = tmr1000_df["CummulativeTranscripts"]

y_min = min(y_min, y.min())

fig.add_trace(
    go.Scatter(
        x=x,
        y=y,
        mode="lines+markers",
        marker=dict(color="green", size=marker_size),
        line=dict(color="green", dash="dash"),
        # name="All",
        # legendgroup=tmr1000_legendtitle,  # this can be any string
        # legendgrouptitle_text=tmr1000_legendtitle,
        name=">= 1000 reads   ",
        legend="legend",
    ),
    row=1,
    col=1,
)


fig.add_trace(
    go.Scatter(
        x=[tmr1000_x_mean],
        y=[tmr1000_y_mean],
        mode="markers+text",
        marker=dict(
            color="green",
            size=marker_size * 2.5,
            # line=dict(
            #     color="yellow",
            #     width=3
            # )
        ),
        showlegend=False,
        text=f"{tmr1000_x_mean:.0f} distinct proteins<br>(avg)",
        textposition="top right",
        textfont=dict(color="green", size=11),
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
    fillcolor="LightSkyBlue",
    # fillcolor="orange",
    # fillcolor="red",
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
        textfont=dict(size=11),
        showlegend=False,
    ),
    row=1,
    col=1,
)

for neural_condition, neural_trace_name, neural_df in zip(
    neural_conditions, neural_trace_names, neural_dfs
):
    if neural_condition == "Missing":
        continue
    x = neural_df["NumOfProteins"]
    y = neural_df["CummulativeTranscripts"]
    color = neural_color_discrete_map[neural_condition]

    fig.add_trace(
        go.Scatter(
            x=x,
            y=y,
            mode="lines+markers",
            marker=dict(color=color, size=marker_size),
            line=dict(color=color, dash="dash", width=0.5),
            name=neural_trace_name,
            # legendgroup=tmr50_legendtitle,  # this can be any string
            # legendgrouptitle_text=tmr50_legendtitle,
            legend="legend1",
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
    text=f"<b>Mann-Whitney U between<br>neural to non-neural genes</b><br>p-val = {pv:.2e}<br>statistic = {statistic:.2g}",
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
    range=[np.log(y_min) * 1.1 / np.log(10), 2.2],
)

width = 800
height = 800

fig.update_layout(
    # xaxis_title="Distinct proteins per transcript",
    # yaxis_title="% of transcripts",
    # title="Pooled octopus data",
    title="Whole-transcriptome octopus data",
    title_x=0.15,
    template=template,
    width=width,
    height=height,
    # legend_title_text=legend_title_text,
    # legend_font=dict(size=10),
    # legend_grouptitlefont=dict(size=12),
    # showlegend=False,
#     legend={
#             # "title": "By country",
#             # "xref": "container",
#             # "yref": "container",
#          "xref": "paper",
#             "yref": "paper",
#             "y": 0.9,
#             # "bgcolor": "Orange",
#         },
#         legend1={
#             # "title": "By continent",
#             # "xref": "container",
#             # "yref": "container",
#              "xref": "paper",
#             "yref": "paper",
#             "y": 0.5,
#             # "bgcolor": "Gold",

#         },
)

fig.write_image(
    "Distinct proteins per gene vs. % of genes - log(y) - Octopus.svg",
    width=width,
    height=height,
)

fig.show()

# %%
df = (
    max_distinct_proteins_df.loc[:, ["NumOfProteins"]]
    .sort_values("NumOfProteins")
    .reset_index(drop=True)
)
df["CummulativeTranscripts"] = 100 * (df.index + 1) / len(df)
df["CummulativeTranscripts"] = df["CummulativeTranscripts"][::-1].values
x = df["NumOfProteins"]
y = df["CummulativeTranscripts"]
x_mean = x.mean()
x_mean_closest = x.iloc[(x - x_mean).abs().argsort()[:1]]
x_mean_closest_k = x_mean_closest.index.values[0]
if x_mean == x_mean_closest.values[0]:
    y_mean = y.iloc[x_mean_closest_k]
else:
    if x_mean < x_mean_closest.values[0]:
        i = x_mean_closest_k - 1
        j = x_mean_closest_k + 1
    else:
        i = x_mean_closest_k
        j = x_mean_closest_k + 2
    y_mean = np.interp(x_mean, x.iloc[i:j], y.iloc[i:j])
df = df.drop_duplicates(subset="NumOfProteins").reset_index(drop=True)

tmr1000_df = (
    tmr1000_max_distinct_proteins_df.loc[:, ["NumOfProteins"]]
    .sort_values("NumOfProteins")
    .reset_index(drop=True)
)
tmr1000_df["CummulativeTranscripts"] = 100 * (tmr1000_df.index + 1) / len(tmr1000_df)
tmr1000_df["CummulativeTranscripts"] = tmr1000_df["CummulativeTranscripts"][::-1].values
tmr1000_x = tmr1000_df["NumOfProteins"]
tmr1000_y = tmr1000_df["CummulativeTranscripts"]
tmr1000_x_mean = tmr1000_x.mean()
tmr1000_x_mean_closest = tmr1000_x.iloc[
    (tmr1000_x - tmr1000_x_mean).abs().argsort()[:1]
]
tmr1000_x_mean_closest_k = tmr1000_x_mean_closest.index.values[0]
if tmr1000_x_mean == tmr1000_x_mean_closest.values[0]:
    tmr1000_y_mean = tmr1000_y.iloc[tmr1000_x_mean_closest_k]
else:
    if tmr1000_x_mean < tmr1000_x_mean_closest.values[0]:
        i = tmr1000_x_mean_closest_k - 1
        j = tmr1000_x_mean_closest_k + 1
    else:
        i = tmr1000_x_mean_closest_k
        j = tmr1000_x_mean_closest_k + 2
    tmr1000_y_mean = np.interp(tmr1000_x_mean, tmr1000_x.iloc[i:j], tmr1000_y.iloc[i:j])
tmr1000_df = tmr1000_df.drop_duplicates(subset="NumOfProteins").reset_index(drop=True)

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

tmr50_legendtitle = "50 reads"
tmr1000_legendtitle = "1000 reads"
legend_title_text = "Minimum coverage per gene      "

marker_size = 4

fig = make_subplots(
    rows=2,
    cols=1,
    x_title="Distinct proteins per gene",
    y_title="% of genes",
    shared_yaxes=True,
    shared_xaxes=True,
    # vertical_spacing=facet_row_spacing / 2.5,
    # horizontal_spacing=facet_col_spacing * 1.5,
    vertical_spacing=0.05,
    # horizontal_spacing=0.025,
)

# tmr50

x = df["NumOfProteins"]
y = df["CummulativeTranscripts"]

y_min = min(y_min, y.min())

tmr50_all_color = "purple"

fig.add_trace(
    go.Scatter(
        x=x,
        y=y,
        mode="lines+markers",
        marker=dict(color=tmr50_all_color, size=marker_size),
        line=dict(color=tmr50_all_color, dash="dash"),
        name="All",
        legendgroup=tmr50_legendtitle,  # this can be any string
        legendgrouptitle_text=tmr50_legendtitle,
    ),
    row=1,
    col=1,
)


fig.add_trace(
    go.Scatter(
        x=[x_mean],
        y=[y_mean],
        mode="markers+text",
        marker=dict(
            color=tmr50_all_color,
            size=marker_size * 2.5,
            # line=dict(
            #     color="yellow",
            #     width=3
            # )
        ),
        showlegend=False,
        # text=f"{x_mean:.0f} distinct proteins<br>(avg)",
        text=f"{x_mean:.0f} distinct<br>proteins<br>(avg)",
        textposition="bottom left",
        textfont=dict(color=tmr50_all_color, size=11),
    ),
    row=1,
    col=1,
)

# tmr1000

x = tmr1000_df["NumOfProteins"]
y = tmr1000_df["CummulativeTranscripts"]

y_min = min(y_min, y.min())

fig.add_trace(
    go.Scatter(
        x=x,
        y=y,
        mode="lines+markers",
        marker=dict(color="green", size=marker_size),
        line=dict(color="green", dash="dash"),
        name="All",
        legendgroup=tmr1000_legendtitle,  # this can be any string
        legendgrouptitle_text=tmr1000_legendtitle,
    ),
    row=1,
    col=1,
)


fig.add_trace(
    go.Scatter(
        x=[tmr1000_x_mean],
        y=[tmr1000_y_mean],
        mode="markers+text",
        marker=dict(
            color="green",
            size=marker_size * 2.5,
            # line=dict(
            #     color="yellow",
            #     width=3
            # )
        ),
        showlegend=False,
        text=f"{tmr1000_x_mean:.0f} distinct proteins<br>(avg)",
        textposition="top right",
        textfont=dict(color="green", size=11),
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
    fillcolor="LightSkyBlue",
    # fillcolor="orange",
    # fillcolor="red",
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
        textfont=dict(size=11),
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
            marker=dict(color=color, size=marker_size),
            line=dict(color=color, dash="dash", width=0.5),
            name=neural_trace_name,
            legendgroup=tmr50_legendtitle,  # this can be any string
            legendgrouptitle_text=tmr50_legendtitle,
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
    text=f"<b>Mann-Whitney U between<br>neural to non-neural genes</b><br>p-val = {pv:.2e}<br>statistic = {statistic:.2g}",
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
    range=[np.log(y_min) * 1.1 / np.log(10), 2.2],
)

width = 800
height = 800

fig.update_layout(
    # xaxis_title="Distinct proteins per transcript",
    # yaxis_title="% of transcripts",
    title="Pooled octopus data",
    title_x=0.15,
    template=template,
    width=width,
    height=height,
    legend_title_text=legend_title_text,
    # legend_font=dict(size=10),
    # legend_grouptitlefont=dict(size=12),
    # showlegend=False,
)

# fig.write_image(
#     "Distinct proteins per gene vs. % of genes - log(y) - Octopus.svg",
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
len(per_transcript_per_sample_coverage_dfs)

# %%
# strongly_diversified_per_transcript_per_sample_coverage_dfs = [
#     df
#     for df, chrom in zip(per_transcript_per_sample_coverage_dfs, chroms)
#     if chrom in strongly_diversified_chroms
# ]
# ic(len(strongly_diversified_per_transcript_per_sample_coverage_dfs))
# strongly_diversified_per_transcript_per_sample_coverage_dfs[0]

# %%
strongly_diversified_per_transcript_per_sample_coverage_dfs = [
    df
    for df, chrom, position_file in zip(per_transcript_per_sample_coverage_dfs, possibly_na_chroms, possibly_na_positions_files)
    if chrom in strongly_diversified_chroms and pd.notna(position_file)
]
ic(len(strongly_diversified_per_transcript_per_sample_coverage_dfs))
strongly_diversified_per_transcript_per_sample_coverage_dfs[0]

# %% jupyter={"source_hidden": true}
# cols = min(facet_col_wrap, len(strongly_diversified_num_of_proteins_per_sample_dfs), 3)
# rows = ceil(len(strongly_diversified_num_of_proteins_per_sample_dfs) / cols)
# row_col_iter = list(product(range(1, rows + 1), range(1, cols + 1)))[
#     : len(strongly_diversified_num_of_proteins_per_sample_dfs)
# ]

# x_title = "Tissue"
# y_title = "Mapped reads"
# # title_text = "Distribution of min & max estimates of non-syn substitutions per read"

# # subplot_titles = strongly_diversified_chroms
# subplot_titles = [
#     f"{transcript.split('_')[0]}<br><sub>({chrom})</sub>"
#     for transcript, chrom in zip(
#         strongly_diversified_transcripts, strongly_diversified_chroms
#     )
# ]

# tissues_order = [
#     "Axial nerve cord",
#     "Frontal & vertical lobe",
#     "Pedunculate & olfactory lobe",
#     "Stellate g. & visceral g.",
#     "Sucker",
#     "Retina & optic lobe",
#     "Non-neuronal tissues mix",
# ]
# tissue_to_legendrank = {tissue: x for x, tissue in enumerate(tissues_order, start=1)}

# fig = make_subplots(
#     rows=rows,
#     cols=cols,
#     subplot_titles=subplot_titles,
#     shared_yaxes=True,
#     x_title=x_title,
#     y_title=y_title,
#     vertical_spacing=0.12,
# )

# max_y = 0
# legend_constructed = False

# for (
#     (row, col),
#     strongly_diversified_per_transcript_per_sample_coverage_df,
#     strongly_diversified_max_num_of_protein,
# ) in zip(
#     row_col_iter,
#     strongly_diversified_per_transcript_per_sample_coverage_dfs,
#     strongly_diversified_max_num_of_proteins,
# ):
#     _tissues = strongly_diversified_num_of_proteins_per_sample_df["Tissue"]

#     for tissue in _tissues:

#         x = [tissue]
#         y = strongly_diversified_per_transcript_per_sample_coverage_df.loc[
#             strongly_diversified_per_transcript_per_sample_coverage_df["Tissue"] == tissue,
#             "NumOfReads",
#         ]

#         max_y = max(max_y, y.max())

#         if not legend_constructed:
#             fig.add_trace(
#                 go.Bar(
#                     x=x,
#                     y=y,
#                     marker_color=tissues_color_discrete_map[tissue],
#                     name=tissue,
#                     legendrank=tissue_to_legendrank[tissue],
#                 ),
#                 row=row,
#                 col=col,
#             )
#         else:
#             fig.add_trace(
#                 go.Bar(
#                     x=x,
#                     y=y,
#                     marker_color=tissues_color_discrete_map[tissue],
#                     # name=tissue,
#                     showlegend=False,
#                 ),
#                 row=row,
#                 col=col,
#             )

#     legend_constructed = True


# # fig.update_xaxes(tickangle=45, automargin=True)

# fig.update_xaxes(
#     showticklabels=False,  # Hide x axis ticks
#     categoryorder="array",
#     categoryarray=tissues_order,
# )
# fig.update_yaxes(
#     range=[0, max_y * 1.1],
#     # range=[0, 0.5+np.log(max_y)/np.log(10)],
#     # type="log"
# )

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

# # fig.write_image(
# #     f"{title_text} - PacBio.svg",
# #     width=max(350 * cols, 800),
# #     height=max(200 * rows, 300),
# # )

# fig.show()


# %%
strongly_diversified_max_num_of_proteins

# %%
cols = min(facet_col_wrap, len(strongly_diversified_num_of_proteins_per_sample_dfs), 3)
rows = ceil(len(strongly_diversified_num_of_proteins_per_sample_dfs) / cols)
row_col_iter = list(product(range(1, rows + 1), range(1, cols + 1)))[
    : len(strongly_diversified_num_of_proteins_per_sample_dfs)
]

x_title = "Coverage"
y_title = "Distinct protein isoforms"

# title_text = "Distribution of min & max estimates of non-syn substitutions per read"


# subplot_titles = [
#     f"{transcript.split('_')[0]}<br><sub>Pooled distinct proteins = {int(total_pooled_isoforms)}</sub>"
#     for transcript, total_pooled_isoforms in zip(
#         strongly_diversified_transcripts, strongly_diversified_max_num_of_proteins
#     )
# ]

subplot_titles = [
    f"{transcript.split('_')[0]} ({int(total_pooled_isoforms)})"
    for transcript, total_pooled_isoforms in zip(
        strongly_diversified_transcripts, strongly_diversified_max_num_of_proteins
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
    strongly_diversified_per_transcript_per_sample_coverage_df,
) in zip(
    row_col_iter,
    strongly_diversified_num_of_proteins_per_sample_dfs,
    strongly_diversified_max_num_of_proteins,
    strongly_diversified_per_transcript_per_sample_coverage_dfs,
):
    _tissues = strongly_diversified_num_of_proteins_per_sample_df["Tissue"]

    for tissue in _tissues:
        try:
            x = strongly_diversified_per_transcript_per_sample_coverage_df.loc[
                strongly_diversified_per_transcript_per_sample_coverage_df["Tissue"]
                == tissue,
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
                    mode="markers",
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
                    mode="markers",
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
    range=[0, max_x * 1.1]
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
    # title_text="Octopus",
    # title_x=0.1,
    # title_y=0.97,
    # showlegend=False,
    legend_title_text="Tissue",
    width=width,
    height=height,
    # barmode="overlay",
)

fig.write_image(
    "Distinct proteins vs coverage per 9 strong transcripts - Octopus.svg",
    width=width,
    height=height,
)

fig.show()


# %% [markdown]
# ### Comparing algorithms

# %% [markdown]
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
# df = max_distinct_proteins_per_transcript_and_alg_and_fracrepetition_df.loc[
#     max_distinct_proteins_per_transcript_and_alg_and_fracrepetition_df[
#         condition_col
#     ].isin(greater_asc_transcripts)
# ]
# df = df.drop(
#     ["AlgorithmRepetition", "IsMaxNumOfProteins", "MaxNumOfProteins", "Proteins"],
#     axis=1,
# )
# df = df.pivot(
#     index=[condition_col, "NumOfReads", "Fraction", "FractionRepetition"],
#     columns="Algorithm",
# )
# # df = df.set_axis(df.columns.get_level_values(1).values, axis=1)
# # df = df.reset_index()
# # df["Desc - Asc"] = df["Descending"] - df["Ascending"]
# # df = df.loc[df["Desc - Asc"] < 0].reset_index(drop=True)
# df

# %%
# df.loc[df["Fraction"] < 1.0]

# %%
# fig = px.histogram(df, x="Desc - Asc", template=template, text_auto=True)

# fig.update_layout(
#     # showlegend=False,
#     # width=1500,
#     # height=300
#     width=600,
#     height=400,
# )

# fig.show()

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

# %% [markdown]
# #### Solutions' sizes

# %%
distinct_unique_proteins_df

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

# %% jupyter={"source_hidden": true}
# len(dispersion_df.loc[dispersion_df["HighDispersion"]])

# %% jupyter={"source_hidden": true}
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

# %% jupyter={"source_hidden": true}
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

# %% jupyter={"source_hidden": true}
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
fig.update_yaxes(title="Genes", type="log")

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
fig = go.Figure(
    go.Histogram(
        y=dispersion_df["%SolutionsDispersion"],
        marker_color="black",
    )
)

# fig.update_xaxes(title="% dispersion<br><sub>100 * (max - min) / max</sub>")
fig.update_yaxes(title="Dispersion [%]", 
                 # type="log"
                )
fig.update_xaxes(title="Genes", 
                 type="log"
                )

fig.update_layout(
    # showlegend=False,
    title="Pooled octopus data",
    title_x=0.15,
    width=600,
    height=400,
    template=template,
)

# fig.write_image("%SolutionsDispersion - Octopus.svg", width=600, height=400)

fig.show()

# %%
grouped_dispersion_df = dispersion_df.groupby("%SolutionsDispersion").size().reset_index().rename(columns={0: "Genes"})
grouped_dispersion_df

# %%
fig = go.Figure(
    go.Bar(
        x=grouped_dispersion_df["Genes"],
        y=grouped_dispersion_df["%SolutionsDispersion"],
        marker_color="black",
    )
)

# fig.update_xaxes(title="% dispersion<br><sub>100 * (max - min) / max</sub>")
fig.update_xaxes(title="Genes", type="log")
fig.update_yaxes(
    title="Dispersion [%]", 
    # type="log"
                )

fig.update_layout(
    # showlegend=False,
    title="Pooled octopus data",
    title_x=0.15,
    width=600,
    height=400,
    template=template,
)

# fig.write_image("%SolutionsDispersion - Octopus.svg", width=600, height=400)

fig.show()

# %% [markdown] jp-MarkdownHeadingCollapsed=true
# Saving dispersion df

# %%
saved_dispersion_df = dispersion_df.rename(columns={"Transcript": "Gene"})
saved_dispersion_df.insert(0, "Platform", "Whole-transcriptome octopus data")
saved_dispersion_df.to_csv("Dispersion.Octopus.tsv", sep="\t", index=False)
saved_dispersion_df

# %% jupyter={"source_hidden": true}
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


# %% [markdown]
# ## Expression levels

# %% [markdown]
# ### Relative expression of isoforms

# %%
# expression_dfs = []
# for expression_file in expression_files:
#     # expression_df = pd.read_csv(expression_file, sep=sep)
#     # expression_df["#Solution"] = expression_df["#Solution"].astype(str)
#     expression_df = pd.read_csv(
#         expression_file,
#         sep=sep,
#         dtype={
#             "#Solution": str,
#             "AdditionalSupportingReadsIDs": str,
#             "AdditionalSupportingProteinsIDs": str,
#         },
#     )

#     # expression_df["Diff5+"] = (
#     #     abs(
#     #         expression_df["TotalEqualSupportingReads"]
#     #         - expression_df["TotalWeightedSupportingReads"]
#     #     )
#     #     >= 0.05
#     #     * (
#     #         expression_df["TotalEqualSupportingReads"]
#     #         + expression_df["TotalWeightedSupportingReads"]
#     #     )
#     #     / 2
#     # )

#     expression_df["AdditionalSupportingReadsIDs"] = expression_df[
#         "AdditionalSupportingReadsIDs"
#     ].apply(lambda x: "" if pd.isna(x) else [y.split(",") for y in x.split(";")])
#     expression_df["AdditionalSupportingProteinsIDs"] = expression_df[
#         "AdditionalSupportingProteinsIDs"
#     ].apply(lambda x: "" if pd.isna(x) else x.split(","))

#     expression_dfs.append(expression_df)
# expression_dfs[0]

# %%
# def find_rand_maximal_solution(
#     expression_df, seed, allowed_algorithms=["Ascending", "Descending"]
# ):
#     df = (
#         expression_df.loc[expression_df["Algorithm"].isin(allowed_algorithms)]
#         .groupby("#Solution")
#         .agg("size")
#         .reset_index()
#         .rename(columns={0: "Size"})
#     )
#     # rand_maximal_solution = df.loc[df["Size"] == df["Size"].max(), "#Solution"].sample(random_state=seed).reset_index(drop=True)
#     rand_maximal_solution = (
#         df.loc[df["Size"] == df["Size"].max(), "#Solution"]
#         .sample(random_state=seed)
#         .values[0]
#     )
#     return rand_maximal_solution

# %%
# def choose_sample_solutions(
#     expression_df, seed, allowed_algorithms=["Ascending", "Descending"]
# ):
#     return (
#         expression_df.loc[
#             expression_df["Algorithm"].isin(allowed_algorithms),
#             [
#                 condition_col,
#                 "#Solution",
#                 "Fraction",
#                 "FractionRepetition",
#                 "Algorithm",
#                 "AlgorithmRepetition",
#             ],
#         ]
#         .groupby(["Algorithm", "#Solution"])
#         .sample()
#         .groupby("Algorithm")
#         .sample(3, random_state=seed)
#         .reset_index(drop=True)["#Solution"]
#     )

# %%
# maximal_solutions = [
#     find_rand_maximal_solution(
#         expression_df, seed, allowed_algorithms=["Ascending", "Descending"]
#     )
#     for expression_df in expression_dfs
# ]
# maximal_solutions[:5]

# %%
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
#     assignment_df["%RelativeExpression"] = (
#         100
#         * assignment_df["TotalWeightedSupportingReads"]
#         / assignment_df["TotalWeightedSupportingReads"].sum()
#     )
#     assignment_df["%CummulativeRelativeExpression"] = assignment_df[
#         "%RelativeExpression"
#     ].cumsum()

# assignment_dfs[0]

# %%
# y_col_name = "TotalWeightedSupportingReads"
# expression_df = (
#     expression_dfs[0].sort_values(y_col_name, ascending=False).reset_index(drop=True)
# )
# expression_df["CummulativeRelativeWeightedExpression"] = expression_df.groupby(
#     "#Solution"
# )[[y_col_name]].transform(lambda x: 100 * x / x.sum())
# expression_df = expression_df.loc[expression_df["#Solution"] == "1000"].reset_index(
#     drop=True
# )
# expression_df

# %%
# for n_top_expressed in [10, 100, 1000]:
#     fig = go.Figure()
#     x = [f"{n_top_expressed} top expressed"] * n_top_expressed + ["Rest"] * (
#         len(expression_df) - n_top_expressed
#     )
#     y = expression_df["MinNonSyns"]
#     fig.add_trace(
#         go.Box(
#             x=x,
#             y=y,
#             # boxpoints='all,
#             # mode="markers",
#             # marker=dict(
#             #     # size=16,
#             #     # cmax=39,
#             #     # cmin=0,
#             #     color=z,
#             #     colorbar=dict(
#             #         title="MinNonSyns"
#             #     ),
#             #     # colorscale="Viridis"
#             # ),
#         )
#     )
#     # fig.update_xaxes(type="log")
#     # fig.update_yaxes(type="log")
#     fig.update_layout(
#         height=400,
#         template=template,
#         yaxis_title="MinNonSyns",
#     )
#     # fig.update_traces(boxpoints='all')
#     fig.show()

# %%
# fig = go.Figure()
# # x = expression_df["MinNonSynsFrequency"]
# x = expression_df.index + 1
# y = expression_df["CummulativeRelativeWeightedExpression"]
# # z = expression_df["MinNonSynsFrequency"]
# z = expression_df["MinNonSyns"]
# fig.add_trace(
#     go.Scattergl(
#         x=x,
#         y=y,
#         mode="markers",
#         marker=dict(
#             # size=16,
#             # cmax=39,
#             # cmin=0,
#             color=z,
#             colorbar=dict(title="MinNonSyns"),
#             # colorscale="Viridis"
#         ),
#     )
# )
# fig.update_xaxes(type="log")
# fig.update_yaxes(type="log")
# fig.update_layout(height=400, template=template)
# fig.show()

# %%
# def find_rand_maximal_solution(
#     expression_df, seed, allowed_algorithms=["Ascending", "Descending"]
# ):
#     df = (
#         expression_df.loc[expression_df["Algorithm"].isin(allowed_algorithms)]
#         .groupby("#Solution")
#         .agg("size")
#         .reset_index()
#         .rename(columns={0: "Size"})
#     )
#     # rand_maximal_solution = df.loc[df["Size"] == df["Size"].max(), "#Solution"].sample(random_state=seed).reset_index(drop=True)
#     rand_maximal_solution = (
#         df.loc[df["Size"] == df["Size"].max(), "#Solution"]
#         .sample(random_state=seed)
#         .values[0]
#     )
#     return rand_maximal_solution

# %%
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
# def choose_sample_solutions(
#     expression_df, seed, allowed_algorithms=["Ascending", "Descending"]
# ):
#     return (
#         expression_df.loc[
#             expression_df["Algorithm"].isin(allowed_algorithms),
#             [
#                 condition_col,
#                 "#Solution",
#                 "Fraction",
#                 "FractionRepetition",
#                 "Algorithm",
#                 "AlgorithmRepetition",
#             ],
#         ]
#         .groupby(["Algorithm", "#Solution"])
#         .sample()
#         .groupby("Algorithm")
#         .sample(3, random_state=seed)
#         .reset_index(drop=True)["#Solution"]
#     )

# %%
# maximal_solutions = [
#     find_rand_maximal_solution(expression_df, seed, allowed_algorithms=["Descending"])
#     for expression_df in expression_dfs
# ]
# maximal_solutions

# %%
# percentile_dfs = [
#     make_percentile_df(
#         expression_df.loc[expression_df["#Solution"] == maximal_solution].reset_index(
#             drop=True
#         ),
#         allowed_algorithms=["Descending"],
#     )
#     for expression_df, maximal_solution in zip(expression_dfs, maximal_solutions)
# ]
# percentile_dfs[0]

# %%
# _seeds = [np.random.default_rng(seed)]
# for _ in conditions[1:]:
#     _seeds.append(np.random.default_rng(_seeds[-1]))

# all_conditions_sample_solutions = [
#     choose_sample_solutions(expression_df, _seed, allowed_algorithms=["Descending"])
#     for expression_df, _seed in zip(expression_dfs, _seeds)
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
# assignment_method = "Weighted"
# y_col_name = "TotalWeightedSupportingReads"

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
# # formulate_equations = [formulate_log10_equation, formulate_semilog10_equation]
# fit_texts = ["    y ~ 1 / sqrt(x)", "    y ~ 1 / log(x)"]

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

# subplot_titles = conditions
# x_axis_name = "Distinct unique protein rank"
# y_axis_name = "Relative expression (%)"
# head_title = f"Relative expression of proteins considering a largest solution in each {str(condition_col).lower()}"

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
#     fit_text,
# ) in zip(
#     row_col_iter,
#     conditions,
#     maximal_dfs,
#     maximal_solutions,
#     linear_spaces,
#     forward_transforms,
#     reverse_transforms,
#     fit_texts,
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

#     train_logspace = [
#         int(i)
#         for i in np.logspace(
#             np.log10(linear_space[0]), np.log10(linear_space[1]), num=1000
#         )
#     ]
#     test_logspace = [
#         int(i)
#         for i in np.logspace(
#             np.log10(linear_space[0] + 20), np.log10(linear_space[1] - 20), num=1000
#         )
#         if int(i) not in train_logspace
#     ]

#     train_x = forward_x_transform(x[train_logspace])
#     train_y = forward_y_transform(y[train_logspace])

#     test_x = forward_x_transform(x[test_logspace])
#     test_y = forward_y_transform(y[test_logspace])

#     # Create linear regression object
#     regr = linear_model.LinearRegression(n_jobs=threads)
#     # Train the model using the training sets
#     regr.fit(np.array(train_x).reshape(-1, 1), train_y)
#     # Make predictions using the testing set
#     pred_y = regr.predict(np.array(test_x).reshape(-1, 1))

#     # transform these variables back to original scale so they can plotted
#     test_x = reverse_x_transform(test_x)
#     pred_y = reverse_y_transform(pred_y)

#     fig.add_trace(
#         fit_scatter_type(
#             x=test_x,
#             y=pred_y,
#             mode="lines",
#             marker_color="grey",
#             line=dict(
#                 dash="dash",
#                 width=regression_line_width,
#             ),
#             # legendgroup=condition,
#             # name=f"{assignment_method} - fitted",
#             showlegend=False,
#         ),
#         row=1,
#         col=col,
#     )

#     i = int(len(test_x) / 10)
#     text_x = test_x.iloc[i] + 2000
#     text_y = pred_y[i] + 0.03
#     # text_x = 1000
#     # text_y = 0.05
#     text_x = np.log10(text_x)
#     text_y = np.log10(text_y)

#     fig.add_annotation(
#         row=row,
#         col=col,
#         x=text_x,
#         y=text_y,
#         xref="x",
#         yref="y",
#         text=fit_text,
#         align="center",
#         font=dict(size=12, color="grey"),
#         showarrow=False,
#     )

# fig.update_layout(
#     title_text=head_title,
#     # title_y=0.95,
#     template=template,
#     showlegend=False,
#     # legend_itemsizing="constant",
#     height=max(400, 200 * rows),
#     # width=max(900, 250 * cols),
# )
# fig.update_xaxes(type="log", nticks=6)
# fig.update_yaxes(type="log")
# fig.write_image(
#     f"{head_title} - PacBio.svg",
#     height=max(400, 200 * rows),
#     width=max(650, 250 * cols),
# )
# fig.show()
# # fig.show(config={'staticPlot': True, 'responsive': False})


# %% [markdown]
# ### Sequneces of ROBO2's 20 most-common proteins

# %% [markdown]
# This section is commented-out not because it's not working (it does), but because it was meant for a specific collaboration at a specific time.

# %%
# robo2_assignment_df = (
#     [
#         assignment_df
#         for assignment_df in assignment_dfs
#         if assignment_df["Transcript"].str.contains("ROBO2").any()
#     ][0]
#     .sort_values("%RelativeExpression", ascending=False)
#     .reset_index(drop=True)
#     .iloc[:20]
# )
# robo2_assignment_df

# %%
# robo2_unique_proteins_df = [
#     unique_proteins_df
#     for unique_proteins_df in unique_proteins_dfs
#     if unique_proteins_df["Transcript"].str.contains("ROBO2").any()
# ][0]
# robo2_unique_proteins_df

# %%
# robo2_unique_proteins_df["NumOfReads"].sum()

# %%
# robo2_20_recoding_sites = (
#     robo2_unique_proteins_df.loc[
#         robo2_unique_proteins_df["Protein"].isin(robo2_assignment_df["Protein"])
#     ]
#     .iloc[:, unique_proteins_first_col_pos:]
#     .reset_index(drop=True)
# )
# robo2_20_recoding_sites

# %%
# top_x_robo2_proteins = min(robo2_20_recoding_sites.shape[0], 20)
# top_x_robo2_proteins

# %%
# robo2_20_recoding_sites.set_axis(
#     [f"ROBO2_{i}" for i in range(top_x_robo2_proteins)], axis="index"
# ).to_csv(f"ROBO2_{top_x_robo2_proteins}_recoding_substitutions.O.vul.tsv", sep="\t")

# %%
# transcriptome_dict = make_fasta_dict(transcriptome_file)
# robo2_mrna_seq = transcriptome_dict[robo2_chrom]
# robo2_mrna_seq

# %%
# robo2_start = starts[robo2_index]
# robo2_end = ends[robo2_index]
# robo2_start, robo2_end

# %%

# %%
# robo2_mrna_seq[robo2_start:robo2_end].translate()

# %%
# robo2_mrna_seq[robo2_start : robo2_start + 3]

# %%
# robo2_mrna_seq[robo2_start : robo2_start + 3].translate()

# %%
# robo2_20_proteins_df = pd.DataFrame(
#     {
#         # f"{x}:{x+3}({robo2_mrna_seq[x:x+3].translate()})": [str(robo2_mrna_seq[x:x+3].translate())] * 20
#         f"{x}:{x+3}({robo2_mrna_seq[x:x+3].translate()})": [
#             str(robo2_mrna_seq[x : x + 3].translate())
#         ]
#         * top_x_robo2_proteins
#         for x in range(robo2_start, robo2_end, 3)
#     }
# )
# robo2_20_proteins_df.update(robo2_20_recoding_sites)
# robo2_20_proteins_df = robo2_20_proteins_df.set_axis(
#     [f"ROBO2_{i}" for i in range(top_x_robo2_proteins)], axis="index"
# )
# robo2_20_proteins_df.to_csv(f"ROBO2_{top_x_robo2_proteins}_seqs.O.vul.tsv", sep="\t")
# robo2_20_proteins_df

# %%
# robo2_20_proteins_simplified_df = robo2_20_proteins_df.applymap(
#     lambda aa: aa if "," not in aa else "X"
# )
# robo2_20_proteins_simplified_df

# %%
# robo2_20_proteins_simplified_seq_records = [
#     SeqRecord(
#         seq=Seq("".join(aa for aa in robo2_20_proteins_simplified_df.iloc[i])),
#         id=f"ROBO2_{i}",
#         description="",
#     )
#     for i in range(top_x_robo2_proteins)
# ]
# simplified_seq_records_output_file = (
#     f"ROBO2_{top_x_robo2_proteins}_simplified_seqs.O.vul.fasta"
# )
# SeqIO.write(
#     robo2_20_proteins_simplified_seq_records,
#     simplified_seq_records_output_file,
#     "fasta",
# )
