# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
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

# %% [markdown] jp-MarkdownHeadingCollapsed=true papermill={"duration": 0.029907, "end_time": "2022-02-01T09:42:43.198426", "exception": false, "start_time": "2022-02-01T09:42:43.168519", "status": "completed"}
# # Imports

# %%
code_dir = "/private7/projects/Combinatorics/Code"

# %%

# %%
# %load_ext autoreload
# %autoreload 2
# # %autosave 600

# %% papermill={"duration": 2.901153, "end_time": "2022-02-01T09:42:46.125355", "exception": false, "start_time": "2022-02-01T09:42:43.224202", "status": "completed"}
import subprocess
import sys
import time
import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import reduce
from itertools import chain, combinations, product
from math import ceil
from multiprocessing import Pool
from pathlib import Path
from random import choice

warnings.filterwarnings("error")

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.colors as pc
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
import scipy.stats
import seaborn as sns
import seaborn.objects as so
from Bio import SeqIO, motifs  # biopython
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from icecream import ic
from logomaker import Logo  # logomaker
from matplotlib_venn import venn2, venn3
from plotly.subplots import make_subplots
from pybedtools import BedTool
from scipy import interpolate  # todo unimport this later?
from seaborn import axes_style
from sklearn import linear_model
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.metrics import mean_squared_error, r2_score
from statsmodels.stats.multitest import fdrcorrection

sys.path.append(str(Path(code_dir).absolute()))
from Alignment.alignment_utils import (
    count_reads,
    count_reads_in_fastq,
    count_reads_in_unaligned_bam,
    count_unique_filtered_aligned_reads,
)
from EditingUtils.logo import multiple_logos_from_fasta_files
from EditingUtils.seq import make_fasta_dict

# %%
pd.set_option("display.max_columns", 500)

# %% [markdown] jp-MarkdownHeadingCollapsed=true
# # Data loading

# %%
condition_col = "Transcript"

orfs_bed = "/private7/projects/Combinatorics/O.vulgaris/Annotations/orfs_oct.bed"
# alignment_stats_file = "/private7/projects/Combinatorics/O.vulgaris/Alignment/PRJNA791920/IsoSeq/AggregatedByChromBySampleSummary.tsv"
alignment_stats_file = "/private10/Projects/David/Kobi_octupus/splitSites/Whole_collapsed_sample.regatedByChromBySampleSummary.tsv"

known_sites_file = (
    "/private7/projects/Combinatorics/O.vulgaris/Annotations/O.vul.EditingSites.csv"
)
transcriptome_file = (
    "/private7/projects/Combinatorics/O.vulgaris/Annotations/orfs_oct.fa"
)

# main_data_dir = Path(
#     "/private7/projects/Combinatorics/O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BHAfterNoise.3/"
# )
main_data_dir = Path(
    "/private7/projects/Combinatorics/O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage50"
)

positions_dir = Path(main_data_dir, "PositionsFiles")
reads_dir = Path(main_data_dir, "ReadsFiles")
proteins_dir = Path(main_data_dir, "ProteinsFiles")
distinct_proteins_dir = Path(main_data_dir, "DistinctProteins")
expression_dir = Path(main_data_dir, "ExpressionLevels")
max_expression_dir = Path(main_data_dir, "MaxExpressionLevels")
max_expression_dir.mkdir(exist_ok=True)

neural_vs_non_neural_expression_file = Path(
    "/private7/projects/Combinatorics/O.vulgaris/Annotations/NeuralVsNonNeuralExpression.csv"
)

samples_and_tissues_file = Path(
    "/private7/projects/Combinatorics/O.vulgaris/Data/PRJNA791920/IsoSeqPolished/samples.csv"
)

raw_reads_info_file = "/private10/Projects/David/Kobi_octupus/splitSites/Whole_collapsed_sample_extracted_read_data.withAnnotations.csv"

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
# samples_and_tissues_df = pd.read_csv(samples_and_tissues_file)
# samples_and_tissues_df["Tissue"] = samples_and_tissues_df["Tissue"].str.capitalize()
# samples_and_tissues_df

# %%
# samples = samples_and_tissues_df["Sample"]
# tissues = samples_and_tissues_df["Tissue"]
# sample_to_tissue_dict = {sample: tissue for sample, tissue in zip(samples, tissues)}
# sample_to_tissue_dict

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
tmr50_alignment_stats_df = alignment_stats_df.loc[
    alignment_stats_df["MappedReads"] >= 50
]
tmr50_alignment_stats_df

# %%
tmr1000_alignment_stats_df = alignment_stats_df.loc[
    alignment_stats_df["MappedReads"] >= 1000
]
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
tmr50_alignment_stats_df.loc[
    ~tmr50_alignment_stats_df["Chrom"].isin(positions_data_df["Chrom"])
]

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
# tmr50_alignment_stats_df.loc[
#     ~tmr50_alignment_stats_df["Chrom"].isin(reads_data_df["Chrom"])
# ]

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
    expression_dir.glob("*.DistinctUniqueProteins.ExpressionLevels.csv")
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
# proteins_data_df["DistinctProteinsFile"].dropna()

# %%
data_df = (
    # orfs_df.merge(alignment_stats_df, on="Chrom", how="left")
    orfs_df.merge(tmr50_alignment_stats_df, on="Chrom", how="right")
    .merge(positions_data_df, on="Chrom", how="left")
    .merge(reads_data_df, on="Chrom", how="left")
    .merge(proteins_data_df, on="Chrom", how="left")
)
# data_df = data_df.loc[data_df["Chrom"].isin(tmr50_alignment_stats_df["Chrom"])].reset_index(drop=True)
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
possibly_na_expression_files = data_df["ExpressionFile"].tolist()

# %%
assert (
    data_df.loc[data_df["ExpressionFile"].notna()].reset_index(drop=True).shape
    == data_df.loc[data_df["DistinctProteinsFile"].notna()].reset_index(drop=True).shape
), "some distinct proteins don't have expression levels"

# %%
complete_data_df = data_df.loc[data_df["ExpressionFile"].notna()].reset_index(drop=True)
# complete_data_df = data_df.loc[data_df["DistinctProteinsFile"].notna()].reset_index(
#     drop=True
# )

# complete_data_df = complete_data_df.drop_duplicates(
#     "Name", keep=False, ignore_index=True
# )
complete_data_df

# %%
# complete_data_df[["Chrom"]].to_csv("TMR50.CompleteData.Chroms.tsv", sep="\t", index=False)

# %%
# complete_data_df.loc[complete_data_df["Chrom"] == robo2_chrom]

# %%
# robo2_index = complete_data_df.loc[complete_data_df["Chrom"] == robo2_chrom].index[0]
# robo2_index

# %%
# complete_data_df["Strand"].value_counts()

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

# %%
# len(data_df["UniqueReadsFile"])

# %%
# len(complete_data_df["UniqueReadsFile"])

# %%
# 100 * len(chroms) / len(possibly_na_chroms)

# %% [markdown] jp-MarkdownHeadingCollapsed=true
# # Data loading - TMR 1000

# %%
# tmr1000_main_data_dir = Path(
#     "/private7/projects/Combinatorics/O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000.BQ30.AHL.BHAfterNoise.3/"
# )
tmr1000_main_data_dir = Path(
    "/private7/projects/Combinatorics/O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage1000"
)
tmr1000_positions_dir = Path(tmr1000_main_data_dir, "PositionsFiles")
tmr1000_reads_dir = Path(tmr1000_main_data_dir, "ReadsFiles")
tmr1000_proteins_dir = Path(tmr1000_main_data_dir, "ProteinsFiles")
tmr1000_distinct_proteins_dir = Path(tmr1000_main_data_dir, "DistinctProteins")
tmr1000_expression_dir = Path(tmr1000_main_data_dir, "ExpressionLevels")

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

tmr1000_expression_files = list(
    tmr1000_expression_dir.glob("*.DistinctUniqueProteins.ExpressionLevels.csv")
)
tmr1000_chroms_in_expression_files = [
    expression_file.name.split(".")[0] for expression_file in tmr1000_expression_files
]


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

tmr1000_expression_data_df = pd.DataFrame(
    {
        "Chrom": tmr1000_chroms_in_expression_files,
        "ExpressionFile": tmr1000_expression_files,
    }
)

tmr1000_proteins_data_df = (
    tmr1000_proteins_data_df.merge(
        tmr1000_unique_proteins_data_df, on="Chrom", how="left"
    )
    .merge(tmr1000_distinct_proteins_data_df, on="Chrom", how="left")
    .merge(tmr1000_expression_data_df, on="Chrom", how="left")
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
tmr1000_complete_data_df = tmr1000_data_df.loc[
    tmr1000_data_df["ExpressionFile"].notna()
].reset_index(drop=True)

# tmr1000_complete_data_df = tmr1000_data_df.loc[
#     tmr1000_data_df["DistinctProteinsFile"].notna()
# ].reset_index(drop=True)

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
tmr1000_possibly_na_unique_proteins_files = tmr1000_data_df[
    "UniqueProteinsFile"
].tolist()
tmr1000_possibly_na_distinct_unique_proteins_files = tmr1000_data_df[
    "DistinctProteinsFile"
].tolist()
tmr1000_possibly_na_expression_files = tmr1000_data_df["ExpressionFile"].tolist()

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
tmr1000_distinct_unique_proteins_files = tmr1000_complete_data_df[
    "DistinctProteinsFile"
].tolist()
tmr1000_expression_files = tmr1000_complete_data_df["ExpressionFile"].tolist()

# %% [markdown]
# # Data loading - neuronal

# %%
neuronal_alignment_stats_file = "/private10/Projects/David/Kobi_octupus/splitSites/neuronal_reads/neuronal_summary.tsv"


neuronal_main_data_dir = Path(
    "/private7/projects/Combinatorics/O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage50.Neuronal/"
)

neuronal_positions_dir = Path(neuronal_main_data_dir, "PositionsFiles")
neuronal_reads_dir = Path(neuronal_main_data_dir, "ReadsFiles")
neuronal_proteins_dir = Path(neuronal_main_data_dir, "ProteinsFiles")
neuronal_distinct_proteins_dir = Path(neuronal_main_data_dir, "DistinctProteins")
neuronal_expression_dir = Path(neuronal_main_data_dir, "ExpressionLevels")
neuronal_max_expression_dir = Path(neuronal_main_data_dir, "MaxExpressionLevels")
neuronal_max_expression_dir.mkdir(exist_ok=True)

neural_vs_non_neural_expression_file = Path(
    "/private7/projects/Combinatorics/O.vulgaris/Annotations/NeuralVsNonNeuralExpression.csv"
)

# samples_and_tissues_file = Path(
#     "/private7/projects/Combinatorics/O.vulgaris/Data/PRJNA791920/IsoSeqPolished/samples.csv"
# )

# %%
# samples_and_tissues_df = pd.read_csv(samples_and_tissues_file)
# samples_and_tissues_df["Tissue"] = samples_and_tissues_df["Tissue"].str.capitalize()
# samples_and_tissues_df

# %%
# samples = samples_and_tissues_df["Sample"]
# tissues = samples_and_tissues_df["Tissue"]
# sample_to_tissue_dict = {sample: tissue for sample, tissue in zip(samples, tissues)}
# sample_to_tissue_dict

# %%
neuronal_alignment_stats_df = pd.read_csv(neuronal_alignment_stats_file, sep="\t")

# alignment_stats_df["UsedForPileup"] = alignment_stats_df.apply(
#     lambda x: (x["Samples"] >= min_samples)
#     and (x["MappedReadsPerSample"] >= min_mapped_reads_per_sample)
#     and (x["KnownSites"] >= min_known_sites),
#     axis=1,
# )

neuronal_alignment_stats_df

# %%
neuronal_tmr50_alignment_stats_df = neuronal_alignment_stats_df.loc[
    neuronal_alignment_stats_df["MappedReads"] >= 50
]
neuronal_tmr50_alignment_stats_df

# %%
# tmr1000_alignment_stats_df = alignment_stats_df.loc[
#     alignment_stats_df["MappedReads"] >= 1000
# ]
# tmr1000_alignment_stats_df

# %%
neuronal_positions_files = list(neuronal_positions_dir.glob("*.positions.csv.gz"))

chroms_in_neuronal_positions = [
    positions_file.name.split(".")[0] for positions_file in neuronal_positions_files
]

neuronal_positions_data_df = pd.DataFrame(
    {
        "Chrom": chroms_in_neuronal_positions,
        "PositionsFile": neuronal_positions_files,
    }
)

neuronal_positions_data_df

# %%
neuronal_tmr50_alignment_stats_df.loc[
    ~neuronal_tmr50_alignment_stats_df["Chrom"].isin(
        neuronal_positions_data_df["Chrom"]
    )
]

# %%
neuronal_reads_files = list(neuronal_reads_dir.glob("*.reads.csv.gz"))
neuronal_chroms_in_reads_files = [
    reads_file.name.split(".")[0] for reads_file in neuronal_reads_files
]

neuronal_unique_reads_files = list(neuronal_reads_dir.glob("*.unique_reads.csv.gz"))
neuronal_chroms_in_unique_reads_files = [
    unique_reads_file.name.split(".")[0]
    for unique_reads_file in neuronal_unique_reads_files
]

neuronal_reads_data_df = pd.DataFrame(
    {
        "Chrom": neuronal_chroms_in_reads_files,
        "ReadsFile": neuronal_reads_files,
    }
)

neuronal_unique_reads_data_df = pd.DataFrame(
    {
        "Chrom": neuronal_chroms_in_unique_reads_files,
        "UniqueReadsFile": neuronal_unique_reads_files,
    }
)

neuronal_reads_data_df = neuronal_reads_data_df.merge(
    neuronal_unique_reads_data_df, on="Chrom", how="left"
)

neuronal_reads_data_df

# %%
# tmr50_alignment_stats_df.loc[
#     ~tmr50_alignment_stats_df["Chrom"].isin(reads_data_df["Chrom"])
# ]

# %%
neuronal_proteins_files = list(neuronal_proteins_dir.glob("*.proteins.csv.gz"))
neuronal_chroms_in_proteins_files = [
    proteins_file.name.split(".")[0] for proteins_file in neuronal_proteins_files
]

neuronal_unique_proteins_files = list(
    neuronal_proteins_dir.glob("*.unique_proteins.csv.gz")
)
neuronal_chroms_in_unique_proteins_files = [
    unique_proteins_file.name.split(".")[0]
    for unique_proteins_file in neuronal_unique_proteins_files
]

# distinct_proteins_files = list(distinct_proteins_dir.glob("*DistinctUniqueProteins.*.csv"))
neuronal_distinct_proteins_files = [
    f
    for f in neuronal_distinct_proteins_dir.glob("*DistinctUniqueProteins.*.csv")
    if "expression" not in f.name.lower()
]
neuronal_chroms_in_distinct_proteins_files = [
    distinct_proteins_file.name.split(".")[0]
    for distinct_proteins_file in neuronal_distinct_proteins_files
]

neuronal_expression_files = list(
    neuronal_expression_dir.glob("*.DistinctUniqueProteins.ExpressionLevels.csv")
)
neuronal_chroms_in_expression_files = [
    expression_file.name.split(".")[0] for expression_file in neuronal_expression_files
]


neuronal_proteins_data_df = pd.DataFrame(
    {
        "Chrom": neuronal_chroms_in_proteins_files,
        "ProteinsFile": neuronal_proteins_files,
    }
)

neuronal_unique_proteins_data_df = pd.DataFrame(
    {
        "Chrom": neuronal_chroms_in_unique_proteins_files,
        "UniqueProteinsFile": neuronal_unique_proteins_files,
    }
)

neuronal_distinct_proteins_data_df = pd.DataFrame(
    {
        "Chrom": neuronal_chroms_in_distinct_proteins_files,
        "DistinctProteinsFile": neuronal_distinct_proteins_files,
    }
)

neuronal_expression_data_df = pd.DataFrame(
    {
        "Chrom": neuronal_chroms_in_expression_files,
        "ExpressionFile": neuronal_expression_files,
    }
)

neuronal_proteins_data_df = (
    neuronal_proteins_data_df.merge(
        neuronal_unique_proteins_data_df, on="Chrom", how="left"
    )
    .merge(neuronal_distinct_proteins_data_df, on="Chrom", how="left")
    .merge(neuronal_expression_data_df, on="Chrom", how="left")
)

neuronal_proteins_data_df

# %%
# proteins_data_df["DistinctProteinsFile"].dropna()

# %%
neuronal_data_df = (
    # orfs_df.merge(alignment_stats_df, on="Chrom", how="left")
    orfs_df.merge(neuronal_tmr50_alignment_stats_df, on="Chrom", how="right")
    .merge(neuronal_positions_data_df, on="Chrom", how="left")
    .merge(neuronal_reads_data_df, on="Chrom", how="left")
    .merge(neuronal_proteins_data_df, on="Chrom", how="left")
)
neuronal_data_df

# %%
neuronal_possibly_na_conditions = neuronal_data_df["Name"].tolist()
neuronal_possibly_na_chroms = neuronal_data_df["Chrom"].tolist()
neuronal_possibly_na_starts = neuronal_data_df["Start"].tolist()
neuronal_possibly_na_ends = neuronal_data_df["End"].tolist()
neuronal_possibly_na_strands = neuronal_data_df["Strand"].tolist()

neuronal_possibly_na_positions_files = neuronal_data_df["PositionsFile"].tolist()
neuronal_possibly_na_reads_files = neuronal_data_df["ReadsFile"].tolist()
neuronal_possibly_na_unique_reads_files = neuronal_data_df["UniqueReadsFile"].tolist()
neuronal_possibly_na_proteins_files = neuronal_data_df["ProteinsFile"].tolist()
neuronal_possibly_na_unique_proteins_files = neuronal_data_df[
    "UniqueProteinsFile"
].tolist()
neuronal_possibly_na_distinct_unique_proteins_files = neuronal_data_df[
    "DistinctProteinsFile"
].tolist()
# expression_files = complete_data_df["ExpressionFile"].tolist()
neuronal_possibly_na_expression_files = neuronal_data_df["ExpressionFile"].tolist()

# %%
assert (
    neuronal_data_df.loc[neuronal_data_df["ExpressionFile"].notna()]
    .reset_index(drop=True)
    .shape
    == neuronal_data_df.loc[neuronal_data_df["DistinctProteinsFile"].notna()]
    .reset_index(drop=True)
    .shape
), "some distinct proteins don't have expression levels"

# %%
neuronal_complete_data_df = neuronal_data_df.loc[
    neuronal_data_df["ExpressionFile"].notna()
].reset_index(drop=True)
# complete_data_df = data_df.loc[data_df["DistinctProteinsFile"].notna()].reset_index(
#     drop=True
# )

# complete_data_df = complete_data_df.drop_duplicates(
#     "Name", keep=False, ignore_index=True
# )
neuronal_complete_data_df

# %%
# complete_data_df[["Chrom"]].to_csv("TMR50.CompleteData.Chroms.tsv", sep="\t", index=False)

# %%
# complete_data_df.loc[complete_data_df["Chrom"] == robo2_chrom]

# %%
# robo2_index = complete_data_df.loc[complete_data_df["Chrom"] == robo2_chrom].index[0]
# robo2_index

# %%
# complete_data_df["Strand"].value_counts()

# %%
# complete_data_df.loc[complete_data_df["Name"].duplicated(keep=False)].sort_values("Name", ignore_index=True)

# %%
neuronal_conditions = neuronal_complete_data_df["Name"].tolist()
neuronal_chroms = neuronal_complete_data_df["Chrom"].tolist()
neuronal_starts = neuronal_complete_data_df["Start"].tolist()
neuronal_ends = neuronal_complete_data_df["End"].tolist()
neuronal_strands = neuronal_complete_data_df["Strand"].tolist()

neuronal_positions_files = neuronal_complete_data_df["PositionsFile"].tolist()
neuronal_reads_files = neuronal_complete_data_df["ReadsFile"].tolist()
neuronal_unique_reads_files = neuronal_complete_data_df["UniqueReadsFile"].tolist()
neuronal_proteins_files = neuronal_complete_data_df["ProteinsFile"].tolist()
neuronal_unique_proteins_files = neuronal_complete_data_df[
    "UniqueProteinsFile"
].tolist()
neuronal_distinct_unique_proteins_files = neuronal_complete_data_df[
    "DistinctProteinsFile"
].tolist()
neuronal_expression_files = neuronal_complete_data_df["ExpressionFile"].tolist()

# %% [markdown]
# # Data loading - non-neuronal

# %%
non_neuronal_alignment_stats_file = "/private10/Projects/David/Kobi_octupus/splitSites/non_neuronal_reads/non_neuronal_summary.tsv"


non_neuronal_main_data_dir = Path(
    "/private7/projects/Combinatorics/O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage50.NonNeuronal/"
)

non_neuronal_positions_dir = Path(non_neuronal_main_data_dir, "PositionsFiles")
non_neuronal_reads_dir = Path(non_neuronal_main_data_dir, "ReadsFiles")
non_neuronal_proteins_dir = Path(non_neuronal_main_data_dir, "ProteinsFiles")
non_neuronal_distinct_proteins_dir = Path(
    non_neuronal_main_data_dir, "DistinctProteins"
)
non_neuronal_expression_dir = Path(non_neuronal_main_data_dir, "ExpressionLevels")
non_neuronal_max_expression_dir = Path(
    non_neuronal_main_data_dir, "MaxExpressionLevels"
)
non_neuronal_max_expression_dir.mkdir(exist_ok=True)

# %%
non_neuronal_alignment_stats_df = pd.read_csv(
    non_neuronal_alignment_stats_file, sep="\t"
)

non_neuronal_alignment_stats_df

# %%
non_neuronal_tmr50_alignment_stats_df = non_neuronal_alignment_stats_df.loc[
    non_neuronal_alignment_stats_df["MappedReads"] >= 50
]
non_neuronal_tmr50_alignment_stats_df

# %%
non_neuronal_positions_files = list(
    non_neuronal_positions_dir.glob("*.positions.csv.gz")
)

chroms_in_non_neuronal_positions = [
    positions_file.name.split(".")[0] for positions_file in non_neuronal_positions_files
]

non_neuronal_positions_data_df = pd.DataFrame(
    {
        "Chrom": chroms_in_non_neuronal_positions,
        "PositionsFile": non_neuronal_positions_files,
    }
)

non_neuronal_positions_data_df

# %%
non_neuronal_tmr50_alignment_stats_df.loc[
    ~non_neuronal_tmr50_alignment_stats_df["Chrom"].isin(
        non_neuronal_positions_data_df["Chrom"]
    )
]

# %%
non_neuronal_reads_files = list(non_neuronal_reads_dir.glob("*.reads.csv.gz"))
non_neuronal_chroms_in_reads_files = [
    reads_file.name.split(".")[0] for reads_file in non_neuronal_reads_files
]

non_neuronal_unique_reads_files = list(
    non_neuronal_reads_dir.glob("*.unique_reads.csv.gz")
)
non_neuronal_chroms_in_unique_reads_files = [
    unique_reads_file.name.split(".")[0]
    for unique_reads_file in non_neuronal_unique_reads_files
]

non_neuronal_reads_data_df = pd.DataFrame(
    {
        "Chrom": non_neuronal_chroms_in_reads_files,
        "ReadsFile": non_neuronal_reads_files,
    }
)

non_neuronal_unique_reads_data_df = pd.DataFrame(
    {
        "Chrom": non_neuronal_chroms_in_unique_reads_files,
        "UniqueReadsFile": non_neuronal_unique_reads_files,
    }
)

non_neuronal_reads_data_df = non_neuronal_reads_data_df.merge(
    non_neuronal_unique_reads_data_df, on="Chrom", how="left"
)

non_neuronal_reads_data_df

# %%
non_neuronal_proteins_files = list(non_neuronal_proteins_dir.glob("*.proteins.csv.gz"))
non_neuronal_chroms_in_proteins_files = [
    proteins_file.name.split(".")[0] for proteins_file in non_neuronal_proteins_files
]

non_neuronal_unique_proteins_files = list(
    neuronal_proteins_dir.glob("*.unique_proteins.csv.gz")
)
non_neuronal_chroms_in_unique_proteins_files = [
    unique_proteins_file.name.split(".")[0]
    for unique_proteins_file in non_neuronal_unique_proteins_files
]

non_neuronal_distinct_proteins_files = [
    f
    for f in non_neuronal_distinct_proteins_dir.glob("*DistinctUniqueProteins.*.csv")
    if "expression" not in f.name.lower()
]
non_neuronal_chroms_in_distinct_proteins_files = [
    distinct_proteins_file.name.split(".")[0]
    for distinct_proteins_file in non_neuronal_distinct_proteins_files
]

non_neuronal_expression_files = list(
    non_neuronal_expression_dir.glob("*.DistinctUniqueProteins.ExpressionLevels.csv")
)
non_neuronal_chroms_in_expression_files = [
    expression_file.name.split(".")[0]
    for expression_file in non_neuronal_expression_files
]


non_neuronal_proteins_data_df = pd.DataFrame(
    {
        "Chrom": non_neuronal_chroms_in_proteins_files,
        "ProteinsFile": non_neuronal_proteins_files,
    }
)

non_neuronal_unique_proteins_data_df = pd.DataFrame(
    {
        "Chrom": non_neuronal_chroms_in_unique_proteins_files,
        "UniqueProteinsFile": non_neuronal_unique_proteins_files,
    }
)

non_neuronal_distinct_proteins_data_df = pd.DataFrame(
    {
        "Chrom": non_neuronal_chroms_in_distinct_proteins_files,
        "DistinctProteinsFile": non_neuronal_distinct_proteins_files,
    }
)

non_neuronal_expression_data_df = pd.DataFrame(
    {
        "Chrom": non_neuronal_chroms_in_expression_files,
        "ExpressionFile": non_neuronal_expression_files,
    }
)

non_neuronal_proteins_data_df = (
    non_neuronal_proteins_data_df.merge(
        non_neuronal_unique_proteins_data_df, on="Chrom", how="left"
    )
    .merge(non_neuronal_distinct_proteins_data_df, on="Chrom", how="left")
    .merge(non_neuronal_expression_data_df, on="Chrom", how="left")
)

non_neuronal_proteins_data_df

# %%
# proteins_data_df["DistinctProteinsFile"].dropna()

# %%
non_neuronal_data_df = (
    # orfs_df.merge(alignment_stats_df, on="Chrom", how="left")
    orfs_df.merge(non_neuronal_tmr50_alignment_stats_df, on="Chrom", how="right")
    .merge(non_neuronal_positions_data_df, on="Chrom", how="left")
    .merge(non_neuronal_reads_data_df, on="Chrom", how="left")
    .merge(non_neuronal_proteins_data_df, on="Chrom", how="left")
)
non_neuronal_data_df

# %%
non_neuronal_possibly_na_conditions = non_neuronal_data_df["Name"].tolist()
non_neuronal_possibly_na_chroms = non_neuronal_data_df["Chrom"].tolist()
non_neuronal_possibly_na_starts = non_neuronal_data_df["Start"].tolist()
non_neuronal_possibly_na_ends = non_neuronal_data_df["End"].tolist()
non_neuronal_possibly_na_strands = non_neuronal_data_df["Strand"].tolist()

non_neuronal_possibly_na_positions_files = non_neuronal_data_df[
    "PositionsFile"
].tolist()
non_neuronal_possibly_na_reads_files = non_neuronal_data_df["ReadsFile"].tolist()
non_neuronal_possibly_na_unique_reads_files = non_neuronal_data_df[
    "UniqueReadsFile"
].tolist()
non_neuronal_possibly_na_proteins_files = non_neuronal_data_df["ProteinsFile"].tolist()
non_neuronal_possibly_na_unique_proteins_files = non_neuronal_data_df[
    "UniqueProteinsFile"
].tolist()
non_neuronal_possibly_na_distinct_unique_proteins_files = non_neuronal_data_df[
    "DistinctProteinsFile"
].tolist()
# expression_files = complete_data_df["ExpressionFile"].tolist()
non_neuronal_possibly_na_expression_files = non_neuronal_data_df[
    "ExpressionFile"
].tolist()

# %%
assert (
    non_neuronal_data_df.loc[non_neuronal_data_df["ExpressionFile"].notna()]
    .reset_index(drop=True)
    .shape
    == non_neuronal_data_df.loc[non_neuronal_data_df["DistinctProteinsFile"].notna()]
    .reset_index(drop=True)
    .shape
), "some distinct proteins don't have expression levels"

# %%
non_neuronal_complete_data_df = non_neuronal_data_df.loc[
    non_neuronal_data_df["ExpressionFile"].notna()
].reset_index(drop=True)

non_neuronal_complete_data_df

# %%
non_neuronal_conditions = non_neuronal_complete_data_df["Name"].tolist()
non_neuronal_chroms = non_neuronal_complete_data_df["Chrom"].tolist()
non_neuronal_starts = non_neuronal_complete_data_df["Start"].tolist()
non_neuronal_ends = non_neuronal_complete_data_df["End"].tolist()
non_neuronal_strands = non_neuronal_complete_data_df["Strand"].tolist()

non_neuronal_positions_files = non_neuronal_complete_data_df["PositionsFile"].tolist()
non_neuronal_reads_files = non_neuronal_complete_data_df["ReadsFile"].tolist()
non_neuronal_unique_reads_files = non_neuronal_complete_data_df[
    "UniqueReadsFile"
].tolist()
non_neuronal_proteins_files = non_neuronal_complete_data_df["ProteinsFile"].tolist()
non_neuronal_unique_proteins_files = non_neuronal_complete_data_df[
    "UniqueProteinsFile"
].tolist()
non_neuronal_distinct_unique_proteins_files = non_neuronal_complete_data_df[
    "DistinctProteinsFile"
].tolist()
non_neuronal_expression_files = non_neuronal_complete_data_df["ExpressionFile"].tolist()


# %% [markdown] jp-MarkdownHeadingCollapsed=true papermill={"duration": 0.040192, "end_time": "2022-02-01T09:42:46.214429", "exception": false, "start_time": "2022-02-01T09:42:46.174237", "status": "completed"}
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
# sample_to_tissue_dict

# %% jupyter={"source_hidden": true}
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

# %% jupyter={"source_hidden": true} papermill={"duration": 0.054755, "end_time": "2022-02-01T09:42:46.304499", "exception": false, "start_time": "2022-02-01T09:42:46.249744", "status": "completed"}
# # # plotly consts
# # color_sequence = px.colors.qualitative.Pastel
# # # color_sequence = px.colors.qualitative.D3
# color_sequence = px.colors.qualitative.G10

# samples_color_discrete_map = {
#     sample: color for sample, color in zip(samples, color_sequence)
# }
# samples_subcolors_discrete_map = {
#     sample: two_subcolors_from_hex(samples_color_discrete_map[sample])
#     for sample in samples
# }

# tissues_color_discrete_map = {
#     tissue: color for tissue, color in zip(tissues, color_sequence)
# }
# tissues_subcolors_discrete_map = {
#     tissue: two_subcolors_from_hex(tissues_color_discrete_map[tissue])
#     for tissue in tissues
# }

# # ic(color_discrete_map)
# # ic(subcolors_discrete_map)
# # category_orders = {condition_col: conditions}
# # horizontal_category_orders = {
# #     category: list(reversed(category_orders[category])) for category in category_orders
# # }
# # # valid_shapes = ['', '/', '\\', 'x', '-', '|', '+', '.']
# # # pattern_shape_map = {
# # #     condition: shape for condition, shape in zip(conditions, cycle(valid_shapes))
# # # }
# facet_col_spacing = 0.05
# template = "plotly_white"
# facet_col_wrap = 6
# facet_row_spacing = facet_col_spacing * 6
# zerolinewidth = 4

# %%
template = "plotly_white"

# %%
pio.templates.default = template

# %%
so.Plot.config.theme.update(axes_style("whitegrid"))

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

# %% [markdown] jp-MarkdownHeadingCollapsed=true papermill={"duration": 0.02598, "end_time": "2022-02-01T09:42:46.438342", "exception": false, "start_time": "2022-02-01T09:42:46.412362", "status": "completed"}
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

# %% [markdown] jp-MarkdownHeadingCollapsed=true papermill={"duration": 0.041741, "end_time": "2022-02-01T09:42:47.760215", "exception": false, "start_time": "2022-02-01T09:42:47.718474", "status": "completed"}
# ## Positions

# %%
# positions_dfs = [
#     pd.read_csv(position_file, sep=sep, dtype={"Reads": str}) for position_file in positions_files
# ]
# for positions_df, condition in zip(positions_dfs, conditions):
#     positions_df.insert(0, condition_col, condition)
# positions_dfs[0]


# %%
# positions_df = pd.concat(positions_dfs, ignore_index=True)
# # positions_df.insert(positions_df.columns.get_loc("G")+1, "ATCGs", positions_df.loc[:, ["A", "T", "C", "G"]].sum(axis=1))
# positions_df

# %%
# def make_concat_positions_df(positions_files, condition_col, conditions):
#     positions_dfs = [
#         pd.read_csv(position_file, sep=sep, dtype={"Reads": str}) for position_file in positions_files
#     ]
#     for positions_df, condition in zip(positions_dfs, conditions):
#         positions_df.insert(0, condition_col, condition)
#     concat_positions_df = pd.concat(positions_dfs, ignore_index=True)
#     return concat_positions_df

# %%
# concat_positions_df = make_concat_positions_df(positions_files, condition_col, conditions)
# concat_positions_df

# %%
# data_df.loc[data_df["PositionsFile"].isna()]

# %%
# possibly_na_positions_files[108]

# %%
# pd.notna(possibly_na_positions_files[108])

# %%
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
def make_concat_all_positions_df(
    possibly_na_positions_files, condition_col, possibly_na_conditions
):
    # all_positions_dfs = [
    #     pd.read_csv(position_file, sep=sep, dtype={"Reads": str})
    #     for position_file in possibly_na_positions_files
    #     if pd.notna(position_file)
    # ]
    # for positions_df, condition, position_file in zip(all_positions_dfs, possibly_na_conditions, possibly_na_positions_files):
    #     if pd.notna(position_file):
    #     positions_df.insert(0, condition_col, condition)

    all_positions_dfs = []
    for position_file, condition in zip(
        possibly_na_positions_files, possibly_na_conditions
    ):
        if pd.isna(position_file):
            continue
        positions_df = pd.read_csv(position_file, sep=sep, dtype={"Reads": str})
        positions_df.insert(0, condition_col, condition)
        all_positions_dfs.append(positions_df)

    concat_all_positions_df = pd.concat(all_positions_dfs, ignore_index=True)
    return concat_all_positions_df


# %%
concat_all_positions_df = make_concat_all_positions_df(
    possibly_na_positions_files, condition_col, possibly_na_conditions
)
concat_all_positions_df

# %%
transcript_and_chrom_from_positions_df = concat_all_positions_df.loc[
    :, ["Chrom", "Transcript"]
].drop_duplicates()
transcript_and_chrom_from_positions_df = transcript_and_chrom_from_positions_df.merge(
    orfs_df.loc[:, ["Chrom", "Name"]].rename(columns={"Name": "RealTranscript"}),
    on="Chrom",
    how="left",
)
transcript_and_chrom_from_positions_df

# %%
transcript_and_chrom_from_positions_df.drop_duplicates("Chrom")

# %%
transcript_and_chrom_from_positions_df.loc[
    transcript_and_chrom_from_positions_df["Transcript"]
    != transcript_and_chrom_from_positions_df["RealTranscript"]
]

# %%
test_cols = [
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
    "Edited",
    "EditingCorrectedPVal",
    "EditedCorrected",
    "EditedFinal",
    "Noise",
    "NoisyCorrected",
]

# %%
concat_all_positions_df.loc[
    (concat_all_positions_df["NoisyCorrected"].fillna(False))
    & (concat_all_positions_df["Noise"] <= 0.1),
    test_cols,
]

# %% editable=true slideshow={"slide_type": ""}
x = (
    concat_all_positions_df.loc[
        (concat_all_positions_df["NoisyCorrected"].fillna(False)), "Noise"
    ]
    * 100
)

fig = go.Figure()
fig.add_trace(
    go.Histogram(
        x=x,
        cumulative_enabled=True,
        # histnorm='percent'
    )
)

fig.update_xaxes(title="Noise [%]")
fig.update_yaxes(type="log", title="Positions")
fig.update_layout(width=700, height=500, template=template)
fig.show()

# %%

# %%

# %%
concat_all_positions_df.loc[concat_all_positions_df["Edited"], test_cols].shape

# %%
concat_all_positions_df.loc[concat_all_positions_df["EditedCorrected"], test_cols].shape

# %%
concat_all_positions_df.loc[
    # all edited positions in all transcripts - including ones whose pooled noise levels is >= 6%
    (concat_all_positions_df["EditedFinal"]),
    test_cols,
].shape

# %%

# %%
# no overlap between noise and editing positions
concat_all_positions_df.loc[
    (
        (concat_all_positions_df["Edited"])
        | (concat_all_positions_df["EditedCorrected"])
        | (concat_all_positions_df["EditedFinal"])
    )
    & (concat_all_positions_df["NoisyCorrected"]),
    test_cols,
].shape

# %%

# %%
concat_all_positions_df.loc[
    concat_all_positions_df["Edited"] & (concat_all_positions_df["Chrom"].isin(chroms)),
    test_cols,
].shape

# %%
concat_all_positions_df.loc[
    concat_all_positions_df["EditedCorrected"]
    & (concat_all_positions_df["Chrom"].isin(chroms)),
    test_cols,
].shape

# %%
len(chroms)

# %%
concat_all_positions_df.loc[
    (concat_all_positions_df["Edited"])
    & (~concat_all_positions_df["EditedCorrected"])
    & (concat_all_positions_df["Chrom"].isin(chroms)),
    test_cols,
].shape

# %%
concat_all_positions_df.loc[
    (~concat_all_positions_df["Edited"])
    & (concat_all_positions_df["EditedCorrected"])
    & (concat_all_positions_df["Chrom"].isin(chroms)),
    test_cols,
].shape

# %%
concat_all_positions_df.loc[
    (concat_all_positions_df["Edited"])
    & (concat_all_positions_df["EditedCorrected"])
    & (concat_all_positions_df["Chrom"].isin(chroms)),
    test_cols,
].shape

# %%
concat_all_positions_df.loc[
    (concat_all_positions_df["EditedFinal"])
    & (concat_all_positions_df["Chrom"].isin(chroms)),
    test_cols,
].shape

# %%
concat_all_positions_df.loc[
    (concat_all_positions_df["EditedFinal"])
    & (
        (~concat_all_positions_df["Edited"])
        | (~concat_all_positions_df["EditedCorrected"])
    )
    & (concat_all_positions_df["Chrom"].isin(chroms)),
    test_cols,
].shape

# %%

# %%

# %%

# %%

# %%

# %%

# %%

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
    concat_all_positions_df.loc[
        (concat_all_positions_df["EditedFinal"])
        & (concat_all_positions_df["Chrom"].isin(chroms))
    ]
    .groupby("Chrom")
    .size()
    .mean()
    .round(2)
)


# %% [markdown]
# ## Raw reads

# %%
neuronality_of_annotaion_df = pd.DataFrame(
    {
        "Annotation": [
            "ACH1",
            "ACH2",
            "ACH3",
            "CCAP",
            "DOP1",
            "DOP2",
            "DOP3",
            "EC",
            "FBL",
            "GABA",
            "GLIA1",
            "GLIA2",
            "GLIA3",
            "GLUT1",
            "GLUT2",
            "GLUT3",
            "GLUT4",
            "HC",
            "IGL1-OA",
            "IGL2-GLUT/DOP",
            # "IGL2-GLUT/DOP",
            "IGL3",
            "IGL4-L11",
            "OA",
            "OGL1",
            "OGL2-DOP",
            "OGL3-OA",
            "PEP-APWG",
            "PEP-Burs",
            "PEP-Fmrfa1",
            "PEP-Fmrfa3",
            "PREC",
            "SERT",
            "SUB",
            "VL",
        ],
        "NeuronalStrRep": [
            "Yes",
            "Yes",
            "Yes",
            "No",
            "Yes",
            "Yes",
            "Yes",
            "No",
            "No",
            "Yes",
            "No",
            "No",
            "No",
            "Yes",
            "Yes",
            "Yes",
            "Yes",
            "No",
            "Yes",
            "Yes",
            # "Yes",
            "Yes",
            "Yes",
            "Yes",
            "Yes",
            "Yes",
            "Yes",
            "Yes",
            "Yes",
            "Yes",
            "Yes",
            "Yes",
            "Yes",
            "Yes",
            "Yes",
        ],
    }
)
# neuronality_of_annotaion_df["NeuronalStrRep"] = neuronality_of_annotaion_df["Neuronal"]
neuronality_of_annotaion_df["Neuronal"] = neuronality_of_annotaion_df[
    "NeuronalStrRep"
].apply(lambda x: True if x == "Yes" else False)
# neuronality_of_annotaion_df = neuronality_of_annotaion_df.drop_duplicates(
#     ignore_index=True
# )
# # .replace({"Yes": True, "No": False})
neuronality_of_annotaion_df

# %%
neuronality_of_annotaion_df.to_csv(
    "/private7/projects/Combinatorics/Code/Data/NeuronalityOfAnnotation.tsv",
    sep=sep,
    index=False,
)

# %%
neuronality_of_annotaion_df["NeuronalStrRep"].value_counts(dropna=False)

# %%
neuronality_of_annotaion_df["Neuronal"].value_counts(dropna=False)

# %%
raw_reads_info_df = pd.read_csv(
    raw_reads_info_file,
    usecols=[
        # "FileName",
        # "FilePath",
        "Gene",
        "Library",
        "ReadID",
        "CB",
        "UB",
        "seurat_clusters",
        "Annotation",
    ],
).rename(
    columns={"Gene": "Chrom", "Library": "Sample", "seurat_clusters": "SeuratCluster"}
)
raw_reads_info_df = raw_reads_info_df.merge(neuronality_of_annotaion_df, how="left")
raw_reads_info_df["NeuronalStrRep"] = raw_reads_info_df["NeuronalStrRep"].fillna("NA")

raw_reads_info_df

# %%
raw_reads_info_df.groupby(["Sample", "CB"]).size()

# %%
(
    so.Plot(reads_per_cell_df, x="Reads", color="Sample").add(so.Area(), so.Hist())
    # .add(so.Bars(), so.Hist(), color="EditingDetectedInChrom")
    # .scale(color=so.Nominal(order=[True, False]))
    # .limit(x=(0, 100))
    .label(
        x="Reads per cell",
        y="Cells",
        # color="Editing detected in gene",
    )
)

# %%
(
    so.Plot(
        raw_reads_info_df.loc[raw_reads_info_df["Chrom"].isin(chroms)], "%Yes/Total"
    )
    .add(so.Area(), so.Hist(), color="EditingDetectedInChrom")
    # .add(so.Bars(), so.Hist(), color="EditingDetectedInChrom")
    .scale(color=so.Nominal(order=[True, False]))
    .limit(x=(0, 100))
    .label(
        x="Reads from neuronal cells / all reads [%]",
        y="Genes",
        color="Editing detected in gene",
    )
)

# %%
raw_reads_info_df["Chrom"].nunique()

# %%
raw_reads_info_df.loc[raw_reads_info_df["Chrom"].isin(chroms), "Chrom"].nunique()

# %%
raw_reads_info_df.loc[
    raw_reads_info_df["Chrom"].isin(chroms), "NeuronalStrRep"
].value_counts(normalize=True).round(2).mul(100)

# %%
raw_reads_info_df.loc[
    (raw_reads_info_df["NeuronalStrRep"].ne("NA")),
    ["Chrom", "NeuronalStrRep"],
]

# %%
raw_reads_info_df.loc[
    (raw_reads_info_df["NeuronalStrRep"].ne("NA")),
    ["Chrom", "NeuronalStrRep"],
].groupby("Chrom")["NeuronalStrRep"].value_counts()

# %%
raw_reads_info_df.loc[
    (raw_reads_info_df["Chrom"].isin(chroms))
    & (raw_reads_info_df["NeuronalStrRep"].ne("NA")),
    "NeuronalStrRep",
].value_counts(normalize=True).round(2).mul(100)

# %%
raw_reads_info_df.loc[(raw_reads_info_df["Chrom"].isin(chroms))].groupby(
    ["Sample", "Neuronal"], dropna=False
)["CB"].value_counts(dropna=False).r

# %%
raw_reads_info_df.loc[(raw_reads_info_df["Chrom"].isin(chroms))].groupby(
    ["Sample", "Neuronal"], dropna=False
)["CB"].nunique(dropna=False)

# %%
raw_reads_info_df.loc[(raw_reads_info_df["Chrom"].isin(chroms))].groupby(
    ["Sample", "Neuronal", "Annotation"], dropna=False
)["CB"].nunique(dropna=False)

# %%
df = Out[890].reset_index()
df

# %%
df.groupby(["Sample", "Neuronal"], dropna=False)["CB"].mean()

# %%
raw_reads_info_df.loc[
    (raw_reads_info_df["Chrom"].isin(chroms))
    & (raw_reads_info_df["NeuronalStrRep"].ne("NA"))
].groupby("Chrom")["NeuronalStrRep"].value_counts()


# %% jupyter={"source_hidden": true}
# neuronal_reads_per_chrom_df = (
#     raw_reads_info_df
#     .groupby("Chrom")["NeuronalStrRep"]
#     .value_counts()
#     .reset_index(name="Reads")
# )

# neuronal_reads_per_chrom_df["MaxReadsPerChrom"] = (
#     neuronal_reads_per_chrom_df.groupby("Chrom")["Reads"].transform("max")
# )
# neuronal_reads_per_chrom_df["IsMaxReadsPerChrom"] = (
#     neuronal_reads_per_chrom_df["Reads"].eq(
#         neuronal_reads_per_chrom_df["MaxReadsPerChrom"]
#     )
# )
# neuronal_reads_per_chrom_df = neuronal_reads_per_chrom_df.merge(
#     (
#         neuronal_reads_per_chrom_df.loc[
#             neuronal_reads_per_chrom_df["IsMaxReadsPerChrom"]
#         ]
#         .groupby("Chrom")
#         .apply(
#             lambda x: x["NeuronalStrRep"].sample(n=1, random_state=seed),
#             include_groups=False,
#         )
#         .reset_index(name="ChromNeuronality")
#         .drop(columns="level_1")
#     ),
#     how="left",
# )

# neuronal_reads_per_chrom_df = neuronal_reads_per_chrom_df.merge(
#     neuronal_reads_per_chrom_df.groupby("Chrom")
#     .apply(
#         lambda x: x.loc[x["NeuronalStrRep"].ne("NA"), "Reads"].max(),
#         include_groups=False,
#     )
#     .reset_index(name="MaxReadsPerChromWoNa"),
#     how="left",
# )
# neuronal_reads_per_chrom_df["IsMaxReadsPerChromWoNa"] = (
#     neuronal_reads_per_chrom_df.apply(
#         lambda x: (
#             x["Reads"] == x["MaxReadsPerChromWoNa"]
#             if x["MaxReadsPerChromWoNa"] != "NA"
#             else np.nan
#         ),
#         axis=1,
#     )
# )

# neuronal_reads_per_chrom_df

# neuronal_reads_per_chrom_df = neuronal_reads_per_chrom_df.merge(
#     (
#         neuronal_reads_per_chrom_df.groupby("Chrom")
#         .apply(
#             lambda x: x.loc[x["IsMaxReadsPerChromWoNa"], "NeuronalStrRep"].sample(
#                 n=1, random_state=seed
#             ),
#             include_groups=False,
#         )
#         .reset_index(name="ChromNeuronalityWoNa")
#         .drop(columns="level_1")
#     ),
#     how="left",
# )

# neuronal_reads_per_chrom_df

# %%
# def decide_chrom_neurality(rng, num_of_neural_reads, num_of_non_neural_reads, num_of_na_reads=None):
#     if num_of_na_reads is not None:
#         nums_of_reads = [num_of_neural_reads, num_of_non_neural_reads, num_of_na_reads]
#         neurality_ids = ["Yes", "No", "NA"]
#     else:
#         nums_of_reads = [num_of_neural_reads, num_of_non_neural_reads]
#         neurality_ids = ["Yes", "No"]
#     # ic(nums_of_reads, neurality_ids)
#     max_num_of_reads = max(nums_of_reads)
#     max_ids = [neurality_id for neurality_id, num_of_reads in zip(neurality_ids, nums_of_reads) if num_of_reads == max_num_of_reads]
#     # ic(nums_of_reads, neurality_ids, max_num_of_reads, max_ids)
#     return rng.choice(max_ids)

# %%
# def decide_complete_chrom_neurality_by_expected(observed_neural_reads, observed_non_neural_reads, expected_neural_reads, observed_non_neural_reads):
#     if num_of_na_reads is not None:
#         nums_of_reads = [num_of_neural_reads, num_of_non_neural_reads, num_of_na_reads]
#         neurality_ids = ["Yes", "No", "NA"]
#     else:
#         nums_of_reads = [num_of_neural_reads, num_of_non_neural_reads]
#         neurality_ids = ["Yes", "No"]
#     # ic(nums_of_reads, neurality_ids)
#     max_num_of_reads = max(nums_of_reads)
#     max_ids = [neurality_id for neurality_id, num_of_reads in zip(neurality_ids, nums_of_reads) if num_of_reads == max_num_of_reads]
#     # ic(nums_of_reads, neurality_ids, max_num_of_reads, max_ids)
#     return rng.choice(max_ids)

# %%
def decide_chrom_neuronality_2(
    neural_reads,
    non_neural_reads,
    na_reads,
    mean_neural_to_non_neural_reads_ratio,
    std_neural_to_non_neural_reads_ratio,
    mean_neural_to_total_reads_prct,
):
    if neural_reads > 0 and non_neural_reads > 0:
        if (
            neural_reads / non_neural_reads
            >= mean_neural_to_non_neural_reads_ratio
            + 3 * std_neural_to_non_neural_reads_ratio
        ):
            return "Yes"
        return "No"
    if neural_reads > 0 and non_neural_reads == 0:
        if (100 * neural_reads) / (
            neural_reads + na_reads
        ) >= mean_neural_to_total_reads_prct:
            return "Yes"
        return "No"
    if neural_reads == 0 and non_neural_reads > 0:
        return "No"
    if neural_reads == non_neural_reads == 0:
        return "NA"
    raise ValueError("something weird")


# %%
neuronal_reads_per_chrom_df = (
    raw_reads_info_df.loc[
        raw_reads_info_df["Chrom"].isin(possibly_na_chroms)
    ]  # use only chroms with sufficient mapped reads
    .groupby("Chrom")["NeuronalStrRep"]
    .value_counts()
    .reset_index(name="Reads")
    .pivot(index=["Chrom"], columns=["NeuronalStrRep"])
    .fillna(0)
)
neuronal_reads_per_chrom_df = neuronal_reads_per_chrom_df.set_axis(
    neuronal_reads_per_chrom_df.columns.get_level_values(1).values, axis=1
).reset_index()
neuronal_reads_per_chrom_df = neuronal_reads_per_chrom_df.loc[
    :, ["Chrom", "Yes", "No", "NA"]
]

neuronal_reads_per_chrom_df["Total"] = neuronal_reads_per_chrom_df.loc[
    :, ["Yes", "No", "NA"]
].sum(axis=1)
# neuronal_reads_per_chrom_df["AllReadsWoNA"] = neuronal_reads_per_chrom_df.loc[:, ["Yes", "No"]].sum(axis=1)

# rng = np.random.default_rng(seed)
# neuronal_reads_per_chrom_df["ChromNeuronality"] = neuronal_reads_per_chrom_df.apply(
#     lambda x: decide_chrom_neurality(rng, x["Yes"], x["No"], x["NA"]),
#     axis=1
# )
# neuronal_reads_per_chrom_df["ChromNeuronalityWoNa"] = neuronal_reads_per_chrom_df.apply(
#     lambda x: decide_chrom_neurality(rng, x["Yes"], x["No"]),
#     axis=1
# )

neuronal_reads_per_chrom_df["%Yes/Total"] = (
    neuronal_reads_per_chrom_df["Yes"]
    .div(neuronal_reads_per_chrom_df["Total"])
    .mul(100)
)

neuronal_reads_per_chrom_df["ChromNeuronality"] = neuronal_reads_per_chrom_df[
    "%Yes/Total"
].ge(
    neuronal_reads_per_chrom_df["%Yes/Total"].mean()
    + neuronal_reads_per_chrom_df["%Yes/Total"].std() * 3
)
neuronal_reads_per_chrom_df["ChromNeuronality"] = neuronal_reads_per_chrom_df[
    "ChromNeuronality"
].apply(lambda x: "Yes" if x else "No")

neuronal_reads_per_chrom_df["Yes/No"] = (
    neuronal_reads_per_chrom_df["Yes"]
    .div(neuronal_reads_per_chrom_df["No"])
    .replace([np.inf, -np.inf], np.nan)
)
neuronal_reads_per_chrom_df["ChromNeuronality2"] = neuronal_reads_per_chrom_df.apply(
    lambda x: decide_chrom_neuronality_2(
        x["Yes"],
        x["No"],
        x["NA"],
        neuronal_reads_per_chrom_df["Yes/No"].mean(),
        neuronal_reads_per_chrom_df["Yes/No"].std(),
        neuronal_reads_per_chrom_df["%Yes/Total"].mean(),
    ),
    axis=1,
)

neuronal_reads_per_chrom_df["EditingDetectedInChrom"] = neuronal_reads_per_chrom_df[
    "Chrom"
].isin(chroms)

neuronal_reads_per_chrom_df

# %%
neuronal_reads_per_chrom_df["%Yes/Total"].describe().round(3)

# %%
neuronal_reads_per_chrom_df["Yes/No"].describe().round(2)

# %%
(
    so.Plot(neuronal_reads_per_chrom_df, "%Yes/Total")
    .add(so.Area(), so.Hist(), color="EditingDetectedInChrom")
    # .add(so.Bars(), so.Hist(), color="EditingDetectedInChrom")
    .scale(color=so.Nominal(order=[True, False]))
    .limit(x=(0, 100))
    .label(
        x="Reads from neuronal cells / all reads [%]",
        y="Genes",
        color="Editing detected in gene",
    )
)

# %%
(
    so.Plot(neuronal_reads_per_chrom_df, "Yes/No")
    .add(so.Area(), so.Hist(), color="EditingDetectedInChrom")
    # .add(so.Bars(), so.Hist(), color="EditingDetectedInChrom")
    .scale(color=so.Nominal(order=[True, False]), x="log")
    .limit(x=(0, None))
    .label(
        x="Reads from neuronal cells / reads from non-neural cells",
        y="Genes",
        color="Editing detected in gene",
    )
)

# %%
# (
#     so.Plot(neuronal_reads_per_chrom_df, y="Total", x="%Yes/Total")
#     .add(so.Dot(), color="EditingDetectedInChrom")
#     # .add(so.Bars(), so.Hist(), color="EditingDetectedInChrom")
#     .scale(color=so.Nominal())
#     .limit(x=(0, 100), y=(0, None))
#     .label(x="Reads from neuronal cells / all reads [%]", y="All reads", color="Editing detected in gene")
# )

# %%
# neuronal_reads_per_edited_chrom_df = (
#     neuronal_reads_per_chrom_df.loc[(neuronal_reads_per_chrom_df["Chrom"].isin(chroms))]
#     .reset_index(drop=True)
# )
neuronal_reads_per_edited_chrom_df = neuronal_reads_per_chrom_df.loc[
    neuronal_reads_per_chrom_df["EditingDetectedInChrom"]
].reset_index(drop=True)

neuronal_reads_per_edited_chrom_df

# %%
neuronal_reads_per_chrom_df["ChromNeuronality"].value_counts()

# %%
neuronal_reads_per_edited_chrom_df["ChromNeuronality"].value_counts()

# %%
neuronal_reads_per_chrom_df["ChromNeuronality2"].value_counts()

# %%
neuronal_reads_per_edited_chrom_df["ChromNeuronality2"].value_counts()

# %%
# neuronal_reads_per_edited_chrom_df["ChromNeuronalityWoNa"].value_counts()

# %%
# neuronal_reads_per_edited_chrom_df.loc[
#     neuronal_reads_per_edited_chrom_df["Chrom"].isin(chroms_with_at_least_5_isoforms)
# ]

# %%
# neuronal_reads_per_edited_chrom_df.loc[
#     neuronal_reads_per_edited_chrom_df["ChromNeuronality"].eq("No"),
#     "ChromNeuronalityWoNa",
# ].value_counts()

# %%
# raw_reads_info_df.loc[
#     raw_reads_info_df["Chrom"].isin(
#         neuronal_reads_per_edited_chrom_df.loc[
#             neuronal_reads_per_edited_chrom_df["ChromNeuronality"].eq("No"), "Chrom"
#         ].unique()
#     ),
# ].groupby("NeuronalStrRep")["Annotation"].value_counts(dropna=False)

# %%

# %%
# num of unique cells
raw_reads_info_df["CB"].nunique()

# %%
# num of unique cells
raw_reads_info_df.groupby("Sample")["CB"].nunique()

# %%
# num of unique cells
raw_reads_info_df.groupby("Sample")["CB"].nunique().sum()

# %%
# num of unique molecules
raw_reads_info_df.drop_duplicates(["CB", "UB"]).shape[0]

# %%
# num of unique molecules per cell
raw_reads_info_df.groupby("CB")["UB"].nunique().describe()

# %%
raw_reads_info_df["Annotation"].nunique()

# %%
set(raw_reads_info_df["Annotation"].unique()) - set(
    neuronality_of_annotaion_df["Annotation"].unique()
)

# %%
raw_reads_info_df["Annotation"].str.startswith("TBA")

# %%
raw_reads_info_df["Annotation"].value_counts(dropna=False)

# %% jupyter={"source_hidden": true}
# expanded_annotation_labels_string_2

# %% jupyter={"source_hidden": true}
# expanded_annotation_labels_string = """A anterior, ACH cholinergic neurons, AcTub acetylated tubulin, CCAP cardioactive peptide cells, DOP dopaminergic neurons, D dorsal, EC endothelial cells, es esophagus, FBL fibroblasts, fu funnel, GABA GABAergic neurons, GLUT glutamatergic neurons, HC hemocytes, igl inner granular layer, IGL inner granular layer cells, me medulla, OA octopaminergic neurons, ogl outer granular layer, OGL outer granular layer cells, ol optic lobe, P posterior, PEP peptidergic neurons, PREC precursor cells, plx plexiform layer, sem supraesophageal mass, SERT serotonergic neurons, sub subesophageal mass, SUB subesophageal neurons, st statocysts, TBA to be annotated, V ventral, vl vertical lobe, VL vertical lobe cells"""

# expanded_annotation_labels_string_2 = """ACH, cholinergic neurons; CCAP, cardioactive peptide cells; DOP, dopaminergic neurons; EC, endothelial cells; FBL, fibroblasts; GABA, GABAergic neurons; GLUT, glutamatergic neurons; HC, hemocytes; IGL, inner granular layer cells; OA, octopaminergic neurons; OGL, outer granular layer cells; PEP, peptidergic neurons; PREC, precursors; SERT, serotonergic neurons; SUB, subesophageal neurons; TBA, to be annotated; VL, vertical lobe cells"""
# expanded_annotation_labels_string_2 = expanded_annotation_labels_string_2.replace(
#     ",", ""
# ).replace(";", ",")

# expanded_annotation_labels_string = (
#     expanded_annotation_labels_string + ", " + expanded_annotation_labels_string_2
# )

# expanded_annotation_labels_dict = {
#     x.split(" ", maxsplit=1)[0]: x.split(" ", maxsplit=1)[1]
#     for x in expanded_annotation_labels_string.split(", ")
# }
# expanded_annotation_labels_dict

# expanded_cell_annotation_df = pd.DataFrame(
#     {
#         "Annotation": expanded_annotation_labels_dict.keys(),
#         "ExpandedAnnotation": expanded_annotation_labels_dict.values(),
#     }
# )
# expanded_cell_annotation_df

# %% jupyter={"source_hidden": true}
# well_annotated_cells_df = (
#     raw_reads_info_df.loc[
#         ~raw_reads_info_df["Annotation"].isin(
#             [
#                 "unstable",
#                 "not separated",
#                 "TBA1",
#                 "TBA2",
#                 "TBA3",
#                 "TBA4",
#                 "TBA5",
#                 "TBA6",
#                 "TBA7",
#                 "TBA8",
#             ]
#         )
#     ]
#     .drop(columns=["Chrom", "ReadID", "UB"])
#     .dropna(axis=0)
#     .drop_duplicates(["Sample", "CB"])
#     # .merge(
#     #     expanded_cell_annotation_df,
#     #     on="Annotation",
#     #     how="left",
#     #     # how="outer",
#     #     indicator="indicator",
#     # )
#     .merge(neuronality_of_annotaion_df)
# )
# # well_annotated_cells_df.insert(
# #     well_annotated_cells_df.columns.get_loc("ExpandedAnnotation") + 1,
# #     "MissingExpandedAnnotation",
# #     well_annotated_cells_df["ExpandedAnnotation"].isna(),
# # )
# # well_annotated_cells_df["MissingExpandedAnnotation"] = well_annotated_cells_df["ExpandedAnnotation"].isna()
# well_annotated_cells_df

# %% jupyter={"source_hidden": true}
# expanded_cell_annotation_with_indicator_df = expanded_cell_annotation_df.merge(
#     well_annotated_cells_df.loc[:, ["Annotation"]].drop_duplicates(),
#     how="left",
#     indicator="indicator",
# )
# expanded_cell_annotation_with_indicator_df.loc[
#     expanded_cell_annotation_with_indicator_df["indicator"] == "left_only"
# ]
# # expanded_cell_annotation_with_indicator_df

# %% jupyter={"source_hidden": true}
# well_annotated_cells_df["indicator"].value_counts()

# %% jupyter={"source_hidden": true}
# # how many cells are missing an expanded annotation
# well_annotated_cells_df.groupby(["Sample", "MissingExpandedAnnotation"]).size()

# %% [markdown] papermill={"duration": 0.041741, "end_time": "2022-02-01T09:42:47.760215", "exception": false, "start_time": "2022-02-01T09:42:47.718474", "status": "completed"}
# ## Neural expression

# %% jupyter={"source_hidden": true}
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
# fig = px.histogram(
#     neural_vs_non_neural_expression_df,
#     x="NeuralObimOrthologs/ObimOrthologs",
#     log_y=True,
# )
# fig.update_layout(width=600, height=400, template=template)
# fig.show()

# %%
fig = px.histogram(neural_vs_non_neural_expression_df, x="IsNeural", log_y=True)
fig.update_layout(width=600, height=400, template=template)
fig.show()

# %%

# %% [markdown] jp-MarkdownHeadingCollapsed=true papermill={"duration": 0.041741, "end_time": "2022-02-01T09:42:47.760215", "exception": false, "start_time": "2022-02-01T09:42:47.718474", "status": "completed"}
# ## Combined neural expression

# %%
combined_per_chrom_neurality_df = (
    neural_vs_non_neural_expression_df.loc[:, ["OvulChrom", "IsNeural"]]
    .rename(columns={"OvulChrom": "Chrom", "IsNeural": "ChromNeuronalityOld"})
    .merge(
        # neuronal_reads_per_chrom_df.drop(columns=["Yes", "No", "NA"]).rename(columns={"ChromNeuronality": "ChromNeuronalityNew"}),
        neuronal_reads_per_chrom_df.rename(
            columns={
                "ChromNeuronality": "ChromNeuronalityNew",
                "ChromNeuronality2": "ChromNeuronalityNew2",
            }
        ),
        # how="outer"
        how="right",
    )
)
# combined_per_chrom_neurality_df["OldAndNewMethodsAgree"] = (
#     combined_per_chrom_neurality_df["ChromNeuronalityOld"].eq(
#         combined_per_chrom_neurality_df["ChromNeuronalityNew"]
#     )
# )
# combined_per_chrom_neurality_df["OldAndNewMethods2Agree"] = (
#     combined_per_chrom_neurality_df["ChromNeuronalityOld"].eq(
#         combined_per_chrom_neurality_df["ChromNeuronalityNew2"]
#     )
# )
combined_per_chrom_neurality_df

# %%
complete_combined_per_chrom_neurality_df = combined_per_chrom_neurality_df.loc[
    ~combined_per_chrom_neurality_df[
        ["ChromNeuronalityOld", "ChromNeuronalityNew", "ChromNeuronalityNew2"]
    ]
    .isna()
    .any(axis=1)
].copy()
# complete_combined_per_chrom_neurality_df["NeuralityAgreementDetails"] = (
#     complete_combined_per_chrom_neurality_df.apply(
#         lambda x: (
#             "Yes"
#             if x["ChromNeuronalityOld"] == "Yes" and x["ChromNeuronalityNew"] == "No"
#             else (
#                 "Old - yes, New - no"
#                 if x["ChromNeuronalityOld"] == "Yes"
#                 else (
#                     "Old - no, New - yes" if x["ChromNeuronalityNew"] == "Yes" else "No"
#                 )
#             )
#         ),
#         axis=1,
#     )
# )
# complete_combined_per_chrom_neurality_df["NeuralityAgreementDetails2"] = (
#     complete_combined_per_chrom_neurality_df.apply(
#         lambda x: (
#             "NA" if x["ChromNeuronalityNew2"] == "NA"
#             else "Yes" if x["ChromNeuronalityOld"] == "Yes" and x["ChromNeuronalityNew2"] == "No"
#             else (
#                 "Old - yes, New - no"
#                 if x["ChromNeuronalityOld"] == "Yes"
#                 else (
#                     "Old - no, New - yes" if x["ChromNeuronalityNew2"] == "Yes" else "No"
#                 )
#             )
#         ),
#         axis=1,
#     )
# )
complete_combined_per_chrom_neurality_df

# %%
complete_combined_per_chrom_neurality_df[
    ["ChromNeuronalityOld", "ChromNeuronalityNew", "ChromNeuronalityNew2"]
].value_counts(dropna=False)

# %%
# combined_per_chrom_neurality_df["OldAndNewMethodsAgree"].value_counts(dropna=False)

# %%
# complete_combined_per_chrom_neurality_df["OldAndNewMethodsAgree"].value_counts()

# %%
# complete_combined_per_chrom_neurality_df["OldAndNewMethods2Agree"].value_counts()

# %%
complete_combined_per_chrom_neurality_df.loc[
    (
        complete_combined_per_chrom_neurality_df["ChromNeuronalityOld"].eq(
            complete_combined_per_chrom_neurality_df["ChromNeuronalityNew"]
        )
    )
    & (
        complete_combined_per_chrom_neurality_df["ChromNeuronalityOld"].eq(
            complete_combined_per_chrom_neurality_df["ChromNeuronalityNew2"]
        )
    )
]

# %%
complete_combined_per_chrom_neurality_df.loc[
    (
        complete_combined_per_chrom_neurality_df["ChromNeuronalityOld"].eq(
            complete_combined_per_chrom_neurality_df["ChromNeuronalityNew"]
        )
    )
    & (
        complete_combined_per_chrom_neurality_df["ChromNeuronalityOld"].eq(
            complete_combined_per_chrom_neurality_df["ChromNeuronalityNew2"]
        )
    )
    & (complete_combined_per_chrom_neurality_df["EditingDetectedInChrom"])
]

# %%
# complete_combined_per_chrom_neurality_df["NeuralityAgreementDetails"].value_counts()

# %%
# complete_combined_per_chrom_neurality_df["NeuralityAgreementDetails2"].value_counts()

# %%
# so.Plot(complete_combined_per_chrom_neurality_df, "%Yes/Total").add(so.Bars(), so.Hist())

# %%
# (
#     so.Plot(complete_combined_per_chrom_neurality_df, "%Yes/Total")
#     .add(so.Area(), so.Hist(), color="NeuralityAgreementDetails")
#     # .add(so.Bars(), so.Hist(), color="NeuralityAgreementDetails")
#     # .scale(color=so.Nominal(order=['Yes', 'Old - yes, New - no', 'Old - no, New - yes', 'No']))
#     .limit(x=(0, 100))
#     .layout(size=(6, 4.5))
#     .label(x="Reads from neuronal cells / all reads [%]", y="Genes")
# )

# %%
# (
#     so.Plot(complete_combined_per_chrom_neurality_df, "%Yes/Total")
#     .facet(col="OldAndNewMethodsAgree")
#     .add(so.Bars(), so.Hist(), color="NeuralityAgreementDetails")
#     # .scale(y="log")
#     .scale(marker=so.Nominal(order=['Yes', 'Old - yes, New - no', 'Old - no, New - yes', "No"]))
#     .limit(x=(0, 100))
# )

# %%
# (
#     so.Plot(complete_combined_per_chrom_neurality_df, "%Yes/Total")
#     .facet(col="NeuralityAgreementDetails", order=['Yes', 'Old - yes, New - no', 'Old - no, New - yes', "No"])
#     # .add(so.Bars(), so.Hist())
#     .add(so.Area(), so.Hist())
#     # .add(so.Area(), so.Hist(), color="NeuralityAgreementDetails")
#     # .scale(y="log")
#     .limit(x=(0, 100))
#     .layout(size=(10, 4))
#     # .label(color=None, marker=None, col=None)
# )

# %% [markdown] papermill={"duration": 0.02598, "end_time": "2022-02-01T09:42:46.438342", "exception": false, "start_time": "2022-02-01T09:42:46.412362", "status": "completed"}
# ## Reads

# %% [markdown] jp-MarkdownHeadingCollapsed=true papermill={"duration": 0.02598, "end_time": "2022-02-01T09:42:46.438342", "exception": false, "start_time": "2022-02-01T09:42:46.412362", "status": "completed"}
# ### All

# %% [markdown]
# That is, all filtered reads.

# %%
len(possibly_na_reads_files)

# %% papermill={"duration": 1.204258, "end_time": "2022-02-01T09:42:47.668206", "exception": false, "start_time": "2022-02-01T09:42:46.463948", "status": "completed"}
# dfs of reads mapped to edited positions
reads_dfs = [
    pd.read_csv(reads_file, sep=sep, dtype={"Read": str}) for reads_file in reads_files
]
for reads_df, chrom in zip(reads_dfs, chroms):
    reads_df.insert(0, "Chrom", chrom)
reads_dfs[0]


# %%
reads_mapped_to_edited_positions_info_dfs = []

for chrom, reads_df in zip(chroms, reads_dfs):
    one_chrom_raw_reads_info_df = raw_reads_info_df.loc[
        raw_reads_info_df["Chrom"] == chrom
    ]
    # one_chrom_raw_reads_info_df

    one_chrom_old_to_new_reads_df = pd.read_table(
        Path(positions_dir, f"{chrom}.OldToNewReads.csv.gz"),
        dtype={"OldRead": str, "NewRead": str},
    )
    # one_chrom_old_to_new_reads_df

    one_chrom_new_reads_info_df = (
        one_chrom_raw_reads_info_df.merge(
            one_chrom_old_to_new_reads_df,
            how="left",
            left_on="ReadID",
            right_on="OldRead",
        )
        .drop(columns=["ReadID", "OldRead"])
        .rename(columns={"NewRead": "Read"})
    )
    # one_chrom_new_reads_info_df

    # one_chrom_edited_positions_df = concat_all_positions_df.loc[
    #     (concat_all_positions_df["Chrom"] == chrom)
    #     & (concat_all_positions_df["EditedFinal"])
    # ]
    # # one_chrom_edited_positions_df
    # one_chrom_unique_reads_mapped_to_edited_positions = list(
    #     set(chain.from_iterable(one_chrom_edited_positions_df["Reads"].str.split(",")))
    # )
    one_chrom_unique_reads_mapped_to_edited_positions = reads_df["Read"].tolist()
    # ic(len(one_chrom_unique_reads_mapped_to_edited_positions))
    # one_chrom_unique_reads_mapped_to_edited_positions[:3]

    reads_mapped_to_edited_positions_info_df = one_chrom_new_reads_info_df.loc[
        one_chrom_new_reads_info_df["Read"].isin(
            one_chrom_unique_reads_mapped_to_edited_positions
        )
    ]

    reads_mapped_to_edited_positions_info_dfs.append(
        reads_mapped_to_edited_positions_info_df
    )
    # break

concat_reads_mapped_to_edited_positions_info_dfs = pd.concat(
    reads_mapped_to_edited_positions_info_dfs, ignore_index=True
)
concat_reads_mapped_to_edited_positions_info_dfs

# %%
neuronal_order = ["Yes", "No"]

# %%
well_annotated_reads_mapped_to_edited_positions_per_cell_df = (
    concat_reads_mapped_to_edited_positions_info_dfs.loc[
        ~concat_reads_mapped_to_edited_positions_info_dfs["Neuronal"].isna()
    ]
    .groupby(["Sample", "Neuronal", "NeuronalStrRep", "Annotation", "CB"], dropna=False)
    .size()
    .reset_index(name="Reads")
)
# Convert NeuronalStrRep to a categorical variable with the desired order
well_annotated_reads_mapped_to_edited_positions_per_cell_df["CatNeuronalStrRep"] = (
    pd.Categorical(
        well_annotated_reads_mapped_to_edited_positions_per_cell_df["NeuronalStrRep"],
        categories=neuronal_order,
        ordered=True,
    )
)
well_annotated_reads_mapped_to_edited_positions_per_cell_df = (
    well_annotated_reads_mapped_to_edited_positions_per_cell_df.sort_values(
        ["CatNeuronalStrRep", "Annotation"], ignore_index=True
    )
)
well_annotated_reads_mapped_to_edited_positions_per_cell_df

# %%
well_annotated_reads_mapped_to_edited_positions_per_annotation_df = (
    concat_reads_mapped_to_edited_positions_info_dfs.loc[
        ~concat_reads_mapped_to_edited_positions_info_dfs["Neuronal"].isna()
    ]
    .groupby(["Sample", "Neuronal", "NeuronalStrRep", "Annotation"], dropna=False)
    .size()
    .reset_index(name="Reads")
)
# Convert NeuronalStrRep to a categorical variable with the desired order
well_annotated_reads_mapped_to_edited_positions_per_annotation_df[
    "CatNeuronalStrRep"
] = pd.Categorical(
    well_annotated_reads_mapped_to_edited_positions_per_annotation_df["NeuronalStrRep"],
    categories=neuronal_order,
    ordered=True,
)
well_annotated_reads_mapped_to_edited_positions_per_annotation_df = (
    well_annotated_reads_mapped_to_edited_positions_per_annotation_df.sort_values(
        ["CatNeuronalStrRep", "Annotation"], ignore_index=True
    )
)
well_annotated_reads_mapped_to_edited_positions_per_annotation_df

# %%
well_annotated_reads_mapped_to_edited_positions_per_cell_df.groupby("Sample")[
    "Reads"
].sum()

# %%
well_annotated_reads_mapped_to_edited_positions_per_cell_df.groupby("Sample")[
    "Reads"
].describe()

# %%
well_annotated_reads_mapped_to_edited_positions_per_cell_df.groupby(
    ["Sample", "Neuronal"]
)["Reads"].describe().round(2)

# %%
well_annotated_reads_mapped_to_edited_positions_per_cell_df.groupby(
    ["Sample", "Annotation"]
)["Reads"].describe()

# %%
well_annotated_reads_mapped_to_edited_positions_per_annotation_df.groupby(
    ["Sample", "Neuronal"]
)["Reads"].describe().round(2)

# %%
width = 4
height = 4

p = (
    so.Plot(
        well_annotated_reads_mapped_to_edited_positions_per_cell_df,
        x="Reads",
        color="NeuronalStrRep",
    )
    .add(so.Area(), so.Hist())
    # .add(so.Bars(), so.Hist())
    .facet(col="Sample")
    .scale(
        # x="log",
        # y="log",
        # y="symlog",
        # color=so.Nominal(order=["Yes", "No", "NA"])
        color=so.Nominal(order=["Yes", "No"])
    )
    .limit(x=(0, None), y=(0, None))
    .label(
        x="Reads per cell",
        y="Cells",
        color="Neuronal",
    )
    .layout(size=(width * 2, height))
)
p

# %%
neuronal_order = ["Yes", "No"]

g = sns.catplot(
    well_annotated_reads_mapped_to_edited_positions_per_cell_df,
    x="Reads",
    y="Annotation",
    col="Sample",
    hue="CatNeuronalStrRep",
    hue_order=neuronal_order,
    kind="boxen",
    # kind="violin",
    # fill=False,
    height=6,
    aspect=0.6,
)

g.set_axis_labels(x_var="Reads per cell", y_var="Cluster")
# Remove duplicate legend entries (keep only the first one)
handles, labels = g.axes.flat[
    0
].get_legend_handles_labels()  # Get handles from first facet
g._legend.set_title("Neuronal")  # Alternative method if needed
g._legend.legendHandles = handles  # Ensure only unique handles appear

g

# %%
well_annotated_reads_mapped_to_edited_positions_per_annotation_df.head()

# %%
# g.set(xticks=[10, 30, 50])

# %%
neuronal_order = ["Yes", "No"]

g = sns.catplot(
    well_annotated_reads_mapped_to_edited_positions_per_annotation_df,
    x="Reads",
    y="Annotation",
    col="Sample",
    hue="CatNeuronalStrRep",
    hue_order=neuronal_order,
    kind="bar",
    # kind="violin",
    # fill=False,
    height=6,
    aspect=0.8,
    # log_scale=True
    # log=True
)

# # Rotate x-tick labels
# for ax in g.axes.flat:
#     ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")

g.set_axis_labels(x_var="Reads from all cells within cluster", y_var="Cluster")
# Remove duplicate legend entries (keep only the first one)
handles, labels = g.axes.flat[
    0
].get_legend_handles_labels()  # Get handles from first facet
g._legend.set_title("Neuronal")  # Alternative method if needed
g._legend.legendHandles = handles  # Ensure only unique handles appear

g

# %%
len(reads_mapped_to_edited_positions_info_dfs)


# %%

# %%
# reads_dfs[1]

# %%
def count_lines_in_file(file, cat_cmd="zcat"):
    cmd = f"{cat_cmd} {file} | wc -l"
    clean_str_output = (
        subprocess.run(cmd, shell=True, capture_output=True)
        .stdout.decode()
        .removesuffix("\n")
    )
    return int(clean_str_output)


# %%
num_of_reads_in_reads_files = pd.Series(
    [count_lines_in_file(reads_file) - 1 for reads_file in reads_files]
)
num_of_reads_in_reads_files

# %%
num_of_reads_in_reads_files.describe()

# %%
num_of_reads_in_reads_files.sum()

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


# %% [markdown] jp-MarkdownHeadingCollapsed=true papermill={"duration": 0.041741, "end_time": "2022-02-01T09:42:47.760215", "exception": false, "start_time": "2022-02-01T09:42:47.718474", "status": "completed"}
# ### Unique

# %% papermill={"duration": 0.126539, "end_time": "2022-02-01T09:42:47.923363", "exception": false, "start_time": "2022-02-01T09:42:47.796824", "status": "completed"}
unique_reads_dfs = [
    pd.read_csv(unique_reads_file, sep=sep, dtype={"UniqueRead": str, "Reads": str})
    for unique_reads_file in unique_reads_files
]
for chrom, unique_reads_df in zip(chroms, unique_reads_dfs):
    unique_reads_df.insert(0, "Chrom", chrom)
unique_reads_dfs[0]


# %%
unique_reads_dfs[0].columns

# %%

# %%

# %%
num_of_reads = []  # summed from unique reads
num_of_mapped_reads = []  # according to the alignment stats
for chrom, unique_reads_df in zip(chroms, unique_reads_dfs):
    num_of_reads.append(unique_reads_df["NumOfReads"].sum())
    num_of_mapped_reads.append(
        tmr50_alignment_stats_df.loc[
            tmr50_alignment_stats_df["Chrom"] == chrom, "MappedReads"
        ].values[0]
    )
num_of_reads_comparison_df = pd.DataFrame(
    {
        "AllMappedReads": num_of_mapped_reads,
        "ReadsInReadsFiles": num_of_reads_in_reads_files,
        "ReadsInUniqueReadsFiles": num_of_reads,
    }
)
num_of_reads_comparison_df["Diff"] = num_of_reads_comparison_df["AllMappedReads"].sub(
    num_of_reads_comparison_df["ReadsInUniqueReadsFiles"]
)
num_of_reads_comparison_df["%Decrease"] = (
    100
    * num_of_reads_comparison_df["Diff"]
    / num_of_reads_comparison_df["AllMappedReads"]
)
num_of_reads_comparison_df

# %%
num_of_reads_comparison_df["ReadsInReadsFiles"].sub(
    num_of_reads_comparison_df["ReadsInUniqueReadsFiles"]
).describe()

# %%
num_of_reads_comparison_df["Diff"].describe()

# %%
num_of_reads_comparison_df["%Decrease"].describe()

# %%

# %%

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

# %% [markdown] jp-MarkdownHeadingCollapsed=true papermill={"duration": 0.041741, "end_time": "2022-02-01T09:42:47.760215", "exception": false, "start_time": "2022-02-01T09:42:47.718474", "status": "completed"}
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
# ### Unique - neuronal

# %%
neuronal_unique_reads_dfs = [
    pd.read_csv(unique_reads_file, sep=sep, dtype={"UniqueRead": str, "Reads": str})
    for unique_reads_file in neuronal_unique_reads_files
]
for chrom, unique_reads_df in zip(neuronal_chroms, neuronal_unique_reads_dfs):
    unique_reads_df.insert(0, "Chrom", chrom)
neuronal_unique_reads_dfs[0]

# %%
len(neuronal_unique_reads_dfs)

# %% [markdown]
# ### Unique - non-neuronal

# %%
non_neuronal_unique_reads_dfs = [
    pd.read_csv(unique_reads_file, sep=sep, dtype={"UniqueRead": str, "Reads": str})
    for unique_reads_file in non_neuronal_unique_reads_files
]
for chrom, unique_reads_df in zip(non_neuronal_chroms, non_neuronal_unique_reads_dfs):
    unique_reads_df.insert(0, "Chrom", chrom)
non_neuronal_unique_reads_dfs[0]

# %%
len(non_neuronal_unique_reads_dfs)

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


# %%
# proteins_df = pd.read_csv(proteins_files[0], sep=sep, dtype={"UniqueRead": str, "Reads": str})
# proteins_df.insert(0, "Chrom", chroms[0])
# proteins_df

# %%
# unique_proteins_df = proteins_df.copy()

# proteins_grouped_by_positions = unique_proteins_df.groupby(
#     unique_proteins_df.columns[14:].tolist()
# )
# # join unique reads
# unique_proteins_df.insert(
#     unique_proteins_df.columns.get_loc("UniqueRead"),
#     "UniqueReads",
#     proteins_grouped_by_positions["UniqueRead"].transform(lambda x: ",".join(x)),
# )
# # drop the col of individual unique reads
# unique_proteins_df = unique_proteins_df.drop(["UniqueRead"], axis=1)
# # count the num of unique reads
# unique_proteins_df.insert(
#     unique_proteins_df.columns.get_loc("UniqueReads") + 1,
#     "NumOfUniqueReads",
#     unique_proteins_df["UniqueReads"].apply(lambda x: len(x.split(","))),
# )
# # first_pos_loc += 1

# # join samples
# unique_proteins_df["Samples"] = proteins_grouped_by_positions["Samples"].transform(
#     lambda x: ",".join(x)
# )

# # join reads
# unique_proteins_df["Reads"] = proteins_grouped_by_positions["Reads"].transform(
#     lambda x: ",".join(x)
# )
# # sum joined reads
# unique_proteins_df["NumOfReads"] = proteins_grouped_by_positions[
#     "NumOfReads"
# ].transform(sum)

# unique_proteins_df

# %%
# unique_proteins_df = (
#     unique_proteins_df.drop_duplicates(
#         subset=unique_proteins_df.columns[first_pos_loc:]
#     )
#     .sort_values("NumOfReads", ascending=False)
#     .reset_index(drop=True)
# )

# assert proteins_df["NumOfReads"].sum() == unique_proteins_df["NumOfReads"].sum()
# assert sum(proteins_df["NumOfReads"] * proteins_df["MinNonSyns"]) == sum(
#     unique_proteins_df["NumOfReads"] * unique_proteins_df["MinNonSyns"]
# )
# assert sum(proteins_df["NumOfReads"] * proteins_df["MaxNonSyns"]) == sum(
#     unique_proteins_df["NumOfReads"] * unique_proteins_df["MaxNonSyns"]
# )

# unique_proteins_df.insert(
#     1, "Protein", serialize_compressed_ids(len(unique_proteins_df))
# )

# %% [markdown] jp-MarkdownHeadingCollapsed=true
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
unique_proteins_dfs[0]

# %%
unique_proteins_dfs[0]["Reads"].str.split(",").apply(len)

# %%
unique_proteins_dfs[0]["Samples"].str.split(",").apply(len)

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
assert (
    len(conditions)
    == len(chroms)
    == len(distinct_unique_proteins_files)
    == len(unique_reads_dfs)
)

distinct_unique_proteins_dfs = []

for condition, chrom, distinct_unique_proteins_file, unique_reads_df in zip(
    conditions, chroms, distinct_unique_proteins_files, unique_reads_dfs
):
    distinct_unique_proteins_df = pd.read_csv(
        distinct_unique_proteins_file, sep=sep, dtype={"AvailableReads": str}
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
    distinct_unique_proteins_df.insert(0, "Chrom", chrom)
    distinct_unique_proteins_df.insert(
        distinct_unique_proteins_df.columns.get_loc("AvailableReads") + 1,
        "NumOfAvailableReads",
        distinct_unique_proteins_df["AvailableReads"].str.split(",").str.len(),
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
    distinct_unique_proteins_df["NumOfReads"]
    != distinct_unique_proteins_df["NumOfAvailableReads"]
]

# %%
distinct_unique_proteins_df["NumOfReads"].sub(
    distinct_unique_proteins_df["NumOfAvailableReads"]
).describe()

# %% jupyter={"source_hidden": true}
# complete_data_df.loc[
#     (complete_data_df["Name"].isin(realizations_count_df.loc[realizations_count_df["Count"] < 16, condition_col]), ["Chrom", "Name"])
# ].values

# %% jupyter={"source_hidden": true}
# num_of_reads_per_transcript_and_fraction_df = (
#     distinct_unique_proteins_df.groupby([condition_col, "Fraction"])["NumOfReads"]
#     .unique()
#     .reset_index()
# )
# # num_of_reads_per_transcript_and_fraction_df = num_of_reads_per_transcript_and_fraction_df.explode("NumOfReads", ignore_index=True)
# num_of_reads_per_transcript_and_fraction_df

# %% jupyter={"source_hidden": true}
# num_of_reads_per_transcript_and_fraction_df["NumOfReads"].apply(len).value_counts()

# %% jupyter={"source_hidden": true}
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


# %% jupyter={"source_hidden": true}
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

# %%

# %% [markdown]
# #### Max distinct per fraction

# %%
# create a temp col which will soon be deleted
distinct_unique_proteins_df["TempIndex"] = distinct_unique_proteins_df.index

# first, create a df with all the largset solutions per chrom and fraction
# (this happens due to n=1, keep="all" in the nlargst function)
max_distinct_per_fraction_df = (
    distinct_unique_proteins_df.groupby(["Chrom", "Fraction"]).apply(
        pd.DataFrame.nlargest,
        n=1,
        keep="all",
        columns="NumOfProteins",
        include_groups=False,
    )
    # .reset_index(drop=True)
)
# use merge to also get the chrom and fraction cols lost due to include_groups=False
max_distinct_per_fraction_df = distinct_unique_proteins_df.loc[
    distinct_unique_proteins_df["TempIndex"].isin(
        max_distinct_per_fraction_df["TempIndex"]
    )
].drop(columns="TempIndex")

del distinct_unique_proteins_df["TempIndex"]

# then, sample a single largest solution (a set of distinct proteins)
# per each fraction in each chrom
max_distinct_per_fraction_df = (
    max_distinct_per_fraction_df.groupby(["Chrom", "Fraction"])
    .sample(n=1, random_state=seed)
    .reset_index(drop=True)
)

max_distinct_per_fraction_df

# %%
max_distinct_per_fraction_df.groupby(["Chrom", "Fraction"]).size().unique()

# %% [markdown]
# #### Max distinct

# %%
# max_distinct_proteins_df = (
#     distinct_unique_proteins_df.sort_values("Fraction", ascending=False)
#     .groupby("Chrom")
#     .apply(pd.DataFrame.nlargest, n=1, columns="NumOfProteins")
# )
# max_distinct_proteins_df = (
#     max_distinct_proteins_df.drop("Chrom", axis=1).reset_index().drop("level_1", axis=1)
# )

max_distinct_proteins_df = max_distinct_per_fraction_df.loc[
    max_distinct_per_fraction_df["Fraction"] == 1
]

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

# max_distinct_proteins_df["NumOfReads"] = max_distinct_proteins_df.apply(
#     lambda x: x["NumOfReads"] if not pd.isna(x["NumOfReads"]) else x["MappedReads"],
#     axis=1,
# )

# # max_distinct_proteins_df = max_distinct_proteins_df.dropna().reset_index(drop=True)

# max_distinct_proteins_df["DistinctProteins/Reads"] = (
#     max_distinct_proteins_df["NumOfProteins"] / max_distinct_proteins_df["NumOfReads"]
# )

max_distinct_proteins_df = max_distinct_proteins_df.merge(
    neural_vs_non_neural_expression_df.loc[:, ["OvulChrom", "IsNeural"]].rename(
        columns={"OvulChrom": "Chrom"}
    ),
    on="Chrom",
    how="left",
)
# max_distinct_proteins_df = max_distinct_proteins_df.merge(
#     neural_vs_non_neural_expression_df.loc[:, ["OvulChrom", "IsNeural"]].rename(
#         columns={"OvulChrom": "Chrom", "IsNeural": "IsNeuralOld"}
#     ),
#     on="Chrom",
#     how="left",
# )
max_distinct_proteins_df["IsNeural"] = max_distinct_proteins_df["IsNeural"].fillna(
    "Missing"
)

# max_distinct_proteins_df = max_distinct_proteins_df.merge(
#     neuronal_reads_per_chrom_df.loc[
#         :, ["Chrom", "ChromNeuronality", "ChromNeuronality2"]
#     ].rename(
#         columns={"ChromNeuronality": "IsNeuralNew", "ChromNeuronality2": "IsNeuralNew2"}
#     ),
#     on="Chrom",
#     how="left",
# )
# max_distinct_proteins_df["IsNeuralNew2"] = max_distinct_proteins_df[
#     "IsNeuralNew2"
# ].replace("NA", "Missing")

max_distinct_proteins_df = max_distinct_proteins_df.sort_values(
    "NumOfProteins", ascending=False, ignore_index=True
)
# max_distinct_proteins_df["CummulativeTranscripts"] = 100 * (max_distinct_proteins_df.index + 1) / len(max_distinct_proteins_df)
# max_distinct_proteins_df["CummulativeTranscripts"] = max_distinct_proteins_df["CummulativeTranscripts"][::-1].values

max_distinct_proteins_df

# %%
max_distinct_proteins_df[["IsNeural", "IsNeuralNew", "IsNeuralNew2"]].value_counts()

# %%
max_distinct_proteins_df[["IsNeural", "IsNeuralNew2"]].value_counts()

# %%
assert max_distinct_proteins_df.loc[
    max_distinct_proteins_df["MappedReads"]
    < max_distinct_proteins_df["NumOfAvailableReads"]
].empty

# %%
max_distinct_proteins_df["MappedReads"].describe()

# %%
max_distinct_proteins_df.loc[
    max_distinct_proteins_df["MappedReads"]
    > max_distinct_proteins_df["NumOfAvailableReads"]
]

# %%
max_distinct_proteins_df["MappedReads"].sub(
    max_distinct_proteins_df["NumOfAvailableReads"]
).describe()

# %%
# max_distinct_proteins_df.loc[max_distinct_proteins_df["Chrom"] == robo2_chrom]

# %%
# mean distinct isoforms per gene
max_distinct_proteins_df["NumOfProteins"].mean().round(2)

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

# %%
# number of chroms with only 1 distinct protein, considering chroms with editing (below 6% noise)
max_distinct_proteins_df.loc[
    (max_distinct_proteins_df["NumOfProteins"] == 1)
    & (max_distinct_proteins_df["Chrom"].isin(chroms))
].shape[0]

# %%
# % of chroms with only 1 distinct protein, considering chroms with editing (below 6% noise)
100 * max_distinct_proteins_df.loc[
    (max_distinct_proteins_df["NumOfProteins"] == 1)
    & (max_distinct_proteins_df["Chrom"].isin(chroms))
].shape[0] / len(chroms)

# %%
# number of chroms with at least 5 distinct proteins, considering chroms with editing (below 6% noise)
max_distinct_proteins_df.loc[
    (max_distinct_proteins_df["NumOfProteins"] >= 5)
    & (max_distinct_proteins_df["Chrom"].isin(chroms))
].shape[0]

# %%
# % of chroms with at least 5 distinct proteins, considering chroms with editing (below 6% noise)
100 * max_distinct_proteins_df.loc[
    (max_distinct_proteins_df["NumOfProteins"] >= 5)
    & (max_distinct_proteins_df["Chrom"].isin(chroms))
].shape[0] / len(chroms)

# %%
# number of chroms with at least 50 distinct proteins, considering chroms with editing (below 6% noise)
max_distinct_proteins_df.loc[
    (max_distinct_proteins_df["NumOfProteins"] >= 50)
    & (max_distinct_proteins_df["Chrom"].isin(chroms))
].shape[0]

# %%
# % of chroms with at least 50 distinct proteins, considering chroms with editing (below 6% noise)
100 * max_distinct_proteins_df.loc[
    (max_distinct_proteins_df["NumOfProteins"] >= 50)
    & (max_distinct_proteins_df["Chrom"].isin(chroms))
].shape[0] / len(chroms)

# %%

# %%
max_distinct_proteins_df["IsNeural"].value_counts()

# %%
max_distinct_proteins_df["IsNeural"].value_counts(normalize=True).mul(100).round(1)

# %%
fig = px.histogram(max_distinct_proteins_df, x="IsNeural", log_y=True)
fig.update_layout(width=600, height=400, template=template)
fig.show()


# %% [markdown]
# ### Expression levels (max, f1, optionally 5+ isoforms)

# %%
def remove_wrapping_quote_marks_from_elements(elements):
    return [x[1:-1] for x in elements]


# %%
# number of chroms with at least 5 distinct proteins
chroms_with_at_least_5_isoforms = max_distinct_proteins_df.loc[
    max_distinct_proteins_df["NumOfProteins"] >= 5, "Chrom"
].tolist()
ic(len(chroms_with_at_least_5_isoforms));


# %%
# chroms_and_exp_files_of_chroms_with_at_least_5_isoforms_sorted_as_all_chroms = [
#     (chrom, expression_file)
#     for chrom, expression_file in zip(chroms, expression_files)
#     if chrom in chroms_with_at_least_5_isoforms
# ]

# %%
# chrom

# %%
# concat_f1_5plus_max_expression_df

# %%
# concat_f1_5plus_max_expression_df["Chrom"].nunique()

# %% jupyter={"source_hidden": true}
# def get_f1_max_expression_df(
#     expression_file, chrom, max_distinct_proteins_df, sep, condition_col
# ):
#     expression_df = pd.read_csv(
#         expression_file,
#         sep=sep,
#         dtype={
#             "#Solution": str,
#             "AdditionalSupportingReadsIDs": str,
#             "AdditionalSupportingProteinsIDs": str,
#         },
#     )
#     expression_df = expression_df.loc[expression_df["Fraction"] == 1].reset_index(
#         drop=True
#     )

#     expression_df = expression_df.drop(
#         columns=[
#             "#Solution",
#             "AmbigousPositions",
#             "EditedPositions",
#             "EditingFrequency",
#             "Index",
#             "NumOfUniqueReads",
#             "Samples",
#             "UneditedPositions",
#             "UniqueReads",
#             "AdditionalEqualSupportingReads",
#             "TotalEqualSupportingReads",
#             "MinNonSyns",
#             "MaxNonSyns",
#             "MinNonSynsFrequency",
#             "MaxNonSynsFrequency",
#             "AdditionalWeightedSupportingReads",
#             "TotalWeightedSupportingReads",
#             # "TotalAdditionalSupportingProteins"
#         ]
#     )

#     expression_df.insert(0, "Chrom", chrom)
#     expression_df = expression_df.merge(
#         max_distinct_proteins_df.loc[
#             max_distinct_proteins_df["Chrom"] == chrom,
#             [
#                 "Chrom",
#                 condition_col,
#                 "Fraction",
#                 "FractionRepetition",
#                 "Algorithm",
#                 "AlgorithmRepetition",
#             ],
#         ],
#         how="right",
#     )

#     expression_df = expression_df.drop(
#         columns=["Fraction", "FractionRepetition", "Algorithm", "AlgorithmRepetition"]
#     )

#     expression_df["Reads"] = (
#         expression_df["Reads"]
#         .str.removeprefix("SubString{String}[")
#         .str.removesuffix("]")
#         .str.split(", ")
#         .apply(remove_wrapping_quote_marks_from_elements)
#     )
#     expression_df["AdditionalSupportingReadsIDs"] = expression_df[
#         "AdditionalSupportingReadsIDs"
#     ].apply(lambda x: "" if pd.isna(x) else [y.split(",") for y in x.split(";")])
#     expression_df["AdditionalSupportingProteinsIDs"] = expression_df[
#         "AdditionalSupportingProteinsIDs"
#     ].apply(lambda x: "" if pd.isna(x) else x.split(","))

#     return expression_df

# %% jupyter={"source_hidden": true}
# def get_f1_5plus_exapnded_max_expression_df(
#     chrom,
#     expression_file,
#     positions_dir,
#     one_chrom_raw_reads_info_df,
#     max_distinct_proteins_df,
#     sep,
#     condition_col,
#     out_file=None,
#     try_using_existing_out_file=True,
# ):
#     if out_file is not None and try_using_existing_out_file:
#         try:
#             expanded_max_expression_df = pd.read_csv(
#                 out_file, sep=sep, dtype={"Protein": str, "Read": str}
#             )
#             if not expanded_max_expression_df.empty:
#                 return expanded_max_expression_df
#         except (FileNotFoundError, EOFError):
#             pass

#     one_chrom_old_to_new_reads_file = Path(
#         positions_dir, f"{chrom}.OldToNewReads.csv.gz"
#     )
#     one_chrom_old_to_new_reads_df = pd.read_table(
#         one_chrom_old_to_new_reads_file, dtype={"OldRead": str, "NewRead": str}
#     )
#     one_chrom_new_reads_info_df = (
#         one_chrom_raw_reads_info_df.merge(
#             one_chrom_old_to_new_reads_df,
#             how="left",
#             left_on="ReadID",
#             right_on="OldRead",
#         )
#         .drop(columns=["ReadID", "OldRead"])
#         .rename(columns={"NewRead": "Read"})
#     )

#     max_expression_df = get_f1_max_expression_df(
#         expression_file, chrom, max_distinct_proteins_df, sep, condition_col
#     )

#     # explode per additional supporting protein
#     supposed_num_of_rows_after_explosion = (
#         max_expression_df["AdditionalSupportingProteins"].sum()
#         + max_expression_df.loc[
#             max_expression_df["AdditionalSupportingProteins"] == 0
#         ].shape[0]
#     )
#     expanded_max_expression_df = (
#         max_expression_df.explode(
#             ["AdditionalSupportingReadsIDs", "AdditionalSupportingProteinsIDs"],
#             ignore_index=True,
#         ).drop(
#             columns=["AdditionalSupportingProteins", "AdditionalSupportingProteinsIDs"]
#         )
#         # .rename(
#         #     columns={
#         #         "AdditionalSupportingProteinsIDs": "AdditionalSupportingProteinID",
#         #         # "AdditionalSupportingProteins": "TotalAdditionalSupportingProteins",
#         #     }
#         # )
#     )
#     assert expanded_max_expression_df.shape[0] == supposed_num_of_rows_after_explosion

#     # further explode the df s.t. for each protein,
#     # there are (newly_supporting_proteins x newly_supproting_reads) rows
#     supposed_num_of_rows_after_explosion = (
#         expanded_max_expression_df["AdditionalSupportingReadsIDs"].apply(len).sum()
#     ) + expanded_max_expression_df.loc[
#         expanded_max_expression_df["AdditionalSupportingReadsIDs"].apply(len).eq(0)
#     ].shape[
#         0
#     ]
#     expanded_max_expression_df = expanded_max_expression_df.explode(
#         ["AdditionalSupportingReadsIDs"], ignore_index=True
#     ).rename(columns={"AdditionalSupportingReadsIDs": "AdditionalSupportingReadID"})
#     assert expanded_max_expression_df.shape[0] == supposed_num_of_rows_after_explosion

#     # further explode the df by theoriginally supporting reads
#     # s.t. there's a row for each original read of
#     # each protein and each of that protein's newly supporting reads
#     supposed_num_of_rows_after_explosion = expanded_max_expression_df[
#         "NumOfReads"
#     ].sum()
#     expanded_max_expression_df = (
#         expanded_max_expression_df.explode(["Reads"], ignore_index=True)
#         .rename(columns={"Reads": "Read"})
#         .drop(columns=["NumOfReads"])
#     )
#     assert expanded_max_expression_df.shape[0] == supposed_num_of_rows_after_explosion
#     assert (
#         expanded_max_expression_df.groupby(["Protein", "Read"])
#         .size()
#         .eq(
#             expanded_max_expression_df.groupby(["Protein", "Read"])[
#                 "AdditionalSupportingReadID"
#             ].nunique()
#         )
#         .all()
#     )

#     # now melt the df s.t. for each protein,
#     # each row consists of one read - an original or an additional one
#     supposed_num_of_rows_after_melt = expanded_max_expression_df.shape[0] * 2
#     expanded_max_expression_df = expanded_max_expression_df.melt(
#         id_vars=[
#             "Chrom",
#             condition_col,
#             "Protein",
#         ],
#         value_vars=["Read", "AdditionalSupportingReadID"],
#         var_name="ReadStatus",
#     )
#     expanded_max_expression_df["ReadStatus"] = expanded_max_expression_df[
#         "ReadStatus"
#     ].replace({"Read": "Original", "AdditionalSupportingReadID": "Additional"})
#     expanded_max_expression_df = expanded_max_expression_df.rename(
#         columns={"value": "Read"}
#     )
#     assert expanded_max_expression_df.shape[0] == supposed_num_of_rows_after_melt
#     # some proteins are not supported by additional reads
#     expanded_max_expression_df = expanded_max_expression_df.loc[
#         ~(
#             (expanded_max_expression_df["Read"].eq(""))
#             & (expanded_max_expression_df["ReadStatus"].eq("Additional"))
#         )
#     ]
#     # due to the melt, some original/additional reads shows more than once
#     expanded_max_expression_df = expanded_max_expression_df.drop_duplicates(
#         ignore_index=True
#     )

#     # merge info per read - new and additional ones alike
#     expanded_max_expression_df = expanded_max_expression_df.merge(
#         one_chrom_new_reads_info_df, how="left"
#     )

#     if expanded_max_expression_df.empty:
#         raise Exception(f"the resulting df of {chrom = } is empty!")

#     if out_file is not None:
#         expanded_max_expression_df.to_csv(out_file, sep=sep, index=False)

#     return expanded_max_expression_df

# %%
def simplified_get_f1_exapnded_max_expression_df(
    chrom,
    expression_file,
    positions_dir,
    one_chrom_raw_reads_info_df,
    max_distinct_proteins_df,
    sep,
    condition_col,
    out_file=None,
    try_using_existing_out_file=True,
):
    if out_file is not None and try_using_existing_out_file:
        try:
            expanded_max_expression_df = pd.read_csv(
                out_file, sep=sep, dtype={"Protein": str, "Read": str}
            )
            if not expanded_max_expression_df.empty:
                return expanded_max_expression_df
        except (FileNotFoundError, EOFError):
            pass

    one_chrom_old_to_new_reads_file = Path(
        positions_dir, f"{chrom}.OldToNewReads.csv.gz"
    )
    one_chrom_old_to_new_reads_df = pd.read_table(
        one_chrom_old_to_new_reads_file, dtype={"OldRead": str, "NewRead": str}
    )
    one_chrom_new_reads_info_df = (
        one_chrom_raw_reads_info_df.merge(
            one_chrom_old_to_new_reads_df,
            how="left",
            left_on="ReadID",
            right_on="OldRead",
        )
        .drop(columns=["ReadID", "OldRead"])
        .rename(columns={"NewRead": "Read"})
    )

    max_expression_df = get_f1_max_expression_df(
        expression_file, chrom, max_distinct_proteins_df, sep, condition_col
    )
    max_expression_df = max_expression_df.drop(
        columns=[
            "NumOfReads",
            "AdditionalSupportingProteinsIDs",
            "AdditionalSupportingProteins",
        ]
    )
    max_expression_df["AdditionalSupportingReadsIDs"] = max_expression_df[
        "AdditionalSupportingReadsIDs"
    ].apply(lambda x: sorted(set(chain.from_iterable(x))))
    max_expression_df["AllReads"] = max_expression_df.apply(
        lambda x: x["Reads"] + x["AdditionalSupportingReadsIDs"], axis=1
    )
    max_expression_df["AllReadsStatuses"] = max_expression_df.apply(
        lambda x: ["Original"] * len(x["Reads"])
        + ["Additional"] * len(x["AdditionalSupportingReadsIDs"]),
        axis=1,
    )
    max_expression_df = max_expression_df.drop(
        columns=[
            "Reads",
            "AdditionalSupportingReadsIDs",
        ]
    )
    assert (
        max_expression_df["AllReads"]
        .apply(len)
        .eq(max_expression_df["AllReadsStatuses"].apply(len))
        .all()
    )

    expanded_max_expression_df = max_expression_df.explode(
        ["AllReads", "AllReadsStatuses"],
        ignore_index=True,
    ).rename(
        columns={
            "AllReads": "Read",
            "AllReadsStatuses": "ReadStatus",
        }
    )
    # some proteins are not supported by additional reads
    expanded_max_expression_df = expanded_max_expression_df.loc[
        ~(
            (expanded_max_expression_df["Read"].eq(""))
            & (expanded_max_expression_df["ReadStatus"].eq("Additional"))
        )
    ]
    assert (
        expanded_max_expression_df.shape[0]
        == expanded_max_expression_df.drop_duplicates().shape[0]
    )
    # merge info per read - new and additional ones alike
    expanded_max_expression_df = expanded_max_expression_df.merge(
        one_chrom_new_reads_info_df, how="left"
    )

    if expanded_max_expression_df.empty:
        raise Exception(f"the resulting df of {chrom = } is empty!")

    if out_file is not None:
        expanded_max_expression_df.to_csv(out_file, sep=sep, index=False)

    return expanded_max_expression_df


# %%
# chroms_with_at_least_5_isoforms

# %% jupyter={"source_hidden": true}
# # x = 1130
# x = 1129

# (
#     chrom,
#     expression_file,
# ) = chroms_and_exp_files_of_chroms_with_at_least_5_isoforms_sorted_as_all_chroms[x]
# chrom, expression_file

# %% jupyter={"source_hidden": true}
# one_chrom_raw_reads_info_df = raw_reads_info_df.loc[raw_reads_info_df["Chrom"] == chrom]
# one_chrom_raw_reads_info_df

# %% jupyter={"source_hidden": true}
# one_chrom_raw_reads_info_df["Neuronal"].value_counts(dropna=False)

# %% jupyter={"source_hidden": true}
# copies_df

# %% jupyter={"source_hidden": true}
# out_file = Path(max_expression_dir, f"{chrom}.gz")
# ic(out_file)
# # assert not out_file.exists()
# assert out_file.exists()

# %% jupyter={"source_hidden": true}
# one_chrom_old_to_new_reads_file = Path(positions_dir, f"{chrom}.OldToNewReads.csv.gz")
# one_chrom_old_to_new_reads_df = pd.read_table(
#     one_chrom_old_to_new_reads_file, dtype={"OldRead": str, "NewRead": str}
# )
# one_chrom_new_reads_info_df = (
#     one_chrom_raw_reads_info_df.merge(
#         one_chrom_old_to_new_reads_df,
#         how="left",
#         left_on="ReadID",
#         right_on="OldRead",
#     )
#     .drop(columns=["ReadID", "OldRead"])
#     .rename(columns={"NewRead": "Read"})
# )

# %% jupyter={"source_hidden": true}
# max_expression_df = get_f1_max_expression_df(
#     expression_file, chrom, max_distinct_proteins_df, sep, condition_col
# )
# max_expression_df = max_expression_df.drop(
#     columns=[
#         "NumOfReads",
#         "AdditionalSupportingProteinsIDs",
#         "AdditionalSupportingProteins",
#     ]
# )
# max_expression_df["AdditionalSupportingReadsIDs"] = max_expression_df[
#     "AdditionalSupportingReadsIDs"
# ].apply(lambda x: sorted(set(chain.from_iterable(x))))
# max_expression_df["AllReads"] = max_expression_df.apply(
#     lambda x: x["Reads"] + x["AdditionalSupportingReadsIDs"], axis=1
# )
# max_expression_df["AllReadsStatuses"] = max_expression_df.apply(
#     lambda x: ["Original"] * len(x["Reads"])
#     + ["Additional"] * len(x["AdditionalSupportingReadsIDs"]),
#     axis=1,
# )
# max_expression_df = max_expression_df.drop(
#     columns=[
#         "Reads",
#         "AdditionalSupportingReadsIDs",
#     ]
# )
# assert max_expression_df["AllReads"].apply(len).eq(max_expression_df["AllReadsStatuses"].apply(len)).all()
# max_expression_df

# %% jupyter={"source_hidden": true}
# expanded_max_expression_df = (
#     max_expression_df.explode(
#         ["AllReads", "AllReadsStatuses"],
#         ignore_index=True,
#     )
#     # .drop(columns=["AdditionalSupportingProteins", "AdditionalSupportingProteinsIDs"])
#     .rename(
#         columns={
#             "AllReads": "Read",
#             "AllReadsStatuses": "ReadStatus",
#         }
#     )
# )

# # some proteins are not supported by additional reads
# expanded_max_expression_df = expanded_max_expression_df.loc[
#     ~(
#         (expanded_max_expression_df["Read"].eq(""))
#         & (expanded_max_expression_df["ReadStatus"].eq("Additional"))
#     )
# ]

# assert expanded_max_expression_df.shape[0] == expanded_max_expression_df.drop_duplicates().shape[0]

# # merge info per read - new and additional ones alike
# expanded_max_expression_df = expanded_max_expression_df.merge(
#     one_chrom_new_reads_info_df, how="left"
# )
# expanded_max_expression_df

# %% jupyter={"source_hidden": true}
# f1_5plus_exapnded_max_expression_dfs = []
# # f1_5plus_exapnded_max_expression_dfs_2 = []

# start_time = time.time()  # Start time of the whole cell
# loop_times = []

# # i = 0
# # i = 1119
# # j = 1

# try_using_previous_out_file = True

# for i, (chrom, expression_file) in enumerate(
#     chroms_and_exp_files_of_chroms_with_at_least_5_isoforms_sorted_as_all_chroms,
#     start=1,
# ):
#     loop_start = time.time()  # Start time of each loop iteration

#     # i += 1
#     state = "start"
#     # if i >= 1120:
#     #     ic(i, state, chrom)
#     # elif i % 50 == 0:
#     #     ic(i, state)
#     if i % 50 == 0 or i == len(
#         chroms_and_exp_files_of_chroms_with_at_least_5_isoforms_sorted_as_all_chroms
#     ):
#         ic(i, state)

#     one_chrom_raw_reads_info_df = raw_reads_info_df.loc[
#         raw_reads_info_df["Chrom"] == chrom
#     ]

#     out_file = Path(max_expression_dir, f"{chrom}.gz")

#     # exapnded_max_expression_df = get_f1_5plus_exapnded_max_expression_df(
#     exapnded_max_expression_df = simplified_get_f1_exapnded_max_expression_df(
#         chrom,
#         expression_file,
#         positions_dir,
#         one_chrom_raw_reads_info_df,
#         max_distinct_proteins_df,
#         sep,
#         condition_col,
#         out_file,
#         try_using_previous_out_file,
#     )

#     f1_5plus_exapnded_max_expression_dfs.append(exapnded_max_expression_df)
#     # f1_5plus_exapnded_max_expression_dfs_2.append(expanded_max_expression_df)

#     state = "end"
#     # if i >= 1140:
#     #     ic(i, state, chrom)
#     # elif i % 50 == 0:
#     #     ic(i, state)
#     if i % 50 == 0 or i == len(
#         chroms_and_exp_files_of_chroms_with_at_least_5_isoforms_sorted_as_all_chroms
#     ):
#         ic(i, state)
#     # i += 1
#     # j += 1

#     loop_end = time.time()  # End time of each loop iteration
#     loop_times.append(loop_end - loop_start)
#     # if i == 10:
#     #     break
#     # if j == 30:
#     #     break

#     # break

# end_time = time.time()  # End time of the whole cell
# loop_times = pd.Series(loop_times)
# print(f"Total execution time: {end_time - start_time:.2f} seconds")
# print(f"Mean execution time: {loop_times.mean():.2f} seconds")
# print(f"Median execution time: {loop_times.median():.2f} seconds")

# ic(
#     len(f1_5plus_exapnded_max_expression_dfs),
#     len(chroms_with_at_least_5_isoforms),
#     len(f1_5plus_exapnded_max_expression_dfs) == len(chroms_with_at_least_5_isoforms),
# )

# f1_5plus_exapnded_max_expression_dfs[0]
# # f1_5plus_exapnded_max_expression_dfs_2[0]

# %%
f1_exapnded_max_expression_dfs = []

start_time = time.time()  # Start time of the whole cell
loop_times = []

try_using_previous_out_file = True
strictly_use_previous_out_file_wo_verification = True

for i, (chrom, expression_file) in enumerate(
    zip(chroms, expression_files),
    start=1,
):
    loop_start = time.time()  # Start time of each loop iteration

    state = "start"
    if i % 50 == 0 or i == len(chroms):
        ic(i, state)

    one_chrom_raw_reads_info_df = raw_reads_info_df.loc[
        raw_reads_info_df["Chrom"] == chrom
    ]

    out_file = Path(max_expression_dir, f"{chrom}.gz")

    if strictly_use_previous_out_file_wo_verification:
        
        expanded_max_expression_df = pd.read_csv(
            out_file, sep=sep, dtype={"Protein": str, "Read": str}
        )

    else:

        expanded_max_expression_df = simplified_get_f1_exapnded_max_expression_df(
            chrom,
            expression_file,
            positions_dir,
            one_chrom_raw_reads_info_df,
            max_distinct_proteins_df,
            sep,
            condition_col,
            out_file,
            try_using_previous_out_file,
        )

    f1_exapnded_max_expression_dfs.append(expanded_max_expression_df)

    state = "end"
    if i % 50 == 0 or i == len(chroms):
        ic(i, state)

    loop_end = time.time()  # End time of each loop iteration
    loop_times.append(loop_end - loop_start)


end_time = time.time()  # End time of the whole cell
loop_times = pd.Series(loop_times)
print(f"Total execution time: {end_time - start_time:.2f} seconds")
print(f"Mean execution time: {loop_times.mean():.2f} seconds")
print(f"Median execution time: {loop_times.median():.2f} seconds")

ic(
    len(f1_exapnded_max_expression_dfs),
    len(chroms),
    len(f1_exapnded_max_expression_dfs) == len(chroms),
)

f1_exapnded_max_expression_dfs[0]

# %%
len(f1_exapnded_max_expression_dfs)

# %%

# %%
# try_using_previous_out_file = True

# parallel_f1_5plus_exapnded_max_expression_inputs_2 = (
#     (
#         chrom,
#         expression_file,
#         positions_dir,
#         raw_reads_info_df.loc[raw_reads_info_df["Chrom"] == chrom],
#         max_distinct_proteins_df,
#         sep,
#         condition_col,
#         Path(max_expression_dir, f"{chrom}.gz"),
#         try_using_previous_out_file,
#     )
#     for chrom, expression_file in zip(chroms, expression_files)
#     if chrom in chroms_with_at_least_5_isoforms[530:]
# )

# # f1_5plus_exapnded_max_expression_dfs = []

# start_time = time.time()  # Start time of the whole cell

# with ThreadPoolExecutor(max_workers=8) as executor:
#     futures = [
#         executor.submit(get_f1_5plus_exapnded_max_expression_df, *args)
#         for args in parallel_f1_5plus_exapnded_max_expression_inputs_2
#     ]

#     # for future in as_completed(futures):
#     #     f1_5plus_exapnded_max_expression_dfs.append(future.result())

# f1_5plus_exapnded_max_expression_dfs_2 = [
#     future.result() for future in as_completed(futures)
# ]

# end_time = time.time()  # End time of the whole cell
# print(f"Total execution time: {end_time - start_time:.4f} seconds")
# print(
#     f"Mean execution time: {(end_time - start_time) / len(f1_5plus_exapnded_max_expression_dfs_2):.4f} seconds"
# )

# ic(len(f1_5plus_exapnded_max_expression_dfs_2))
# f1_5plus_exapnded_max_expression_dfs_2[0]

# %% [markdown]
# ### Distinct unique proteins - TMR 1000

# %%
assert (
    len(tmr1000_conditions)
    == len(tmr1000_chroms)
    == len(tmr1000_distinct_proteins_files)
    == len(tmr1000_unique_reads_dfs)
)

tmr1000_distinct_unique_proteins_dfs = []
for condition, chrom, distinct_unique_proteins_file, unique_reads_df in zip(
    tmr1000_conditions,
    tmr1000_chroms,
    tmr1000_distinct_proteins_files,
    tmr1000_unique_reads_dfs,
):
    tmr1000_distinct_unique_proteins_df = pd.read_csv(
        distinct_unique_proteins_file, sep=sep, dtype={"AvailableReads": str}
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


# %% [markdown]
# #### Max distinct per fraction - TMR 1000

# %%
# create a temp col which will soon be deleted
tmr1000_distinct_unique_proteins_df["TempIndex"] = (
    tmr1000_distinct_unique_proteins_df.index
)

# first, create a df with all the largset solutions per chrom and fraction
# (this happens due to n=1, keep="all" in the nlargst function)
tmr1000_max_distinct_per_fraction_df = (
    tmr1000_distinct_unique_proteins_df.groupby(["Chrom", "Fraction"]).apply(
        pd.DataFrame.nlargest,
        n=1,
        keep="all",
        columns="NumOfProteins",
        include_groups=False,
    )
    # .reset_index(drop=True)
)
# use merge to also get the chrom and fraction cols lost due to include_groups=False
tmr1000_max_distinct_per_fraction_df = tmr1000_distinct_unique_proteins_df.loc[
    tmr1000_distinct_unique_proteins_df["TempIndex"].isin(
        tmr1000_max_distinct_per_fraction_df["TempIndex"]
    )
].drop(columns="TempIndex")

del tmr1000_distinct_unique_proteins_df["TempIndex"]

# then, sample a single largest solution (a set of distinct proteins)
# per each fraction in each chrom
tmr1000_max_distinct_per_fraction_df = (
    tmr1000_max_distinct_per_fraction_df.groupby(["Chrom", "Fraction"])
    .sample(n=1, random_state=seed)
    .reset_index(drop=True)
)

tmr1000_max_distinct_per_fraction_df

# %%
tmr1000_max_distinct_per_fraction_df.groupby(["Chrom", "Fraction"]).size().unique()

# %% [markdown]
# #### Max distinct - TMR 1000

# %%
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
tmr1000_max_distinct_proteins_df = tmr1000_max_distinct_per_fraction_df.loc[
    tmr1000_max_distinct_per_fraction_df["Fraction"] == 1
]

# # max_distinct_proteins_df[condition_col] = max_distinct_proteins_df[
# #     condition_col
# # ].astype(str)

tmr1000_max_distinct_proteins_df = tmr1000_max_distinct_proteins_df.merge(
    tmr1000_alignment_stats_df,
    on="Chrom",
    # how="left",
    how="right",
)

tmr1000_max_distinct_proteins_df["NumOfProteins"] = tmr1000_max_distinct_proteins_df[
    "NumOfProteins"
].fillna(1)

# tmr1000_max_distinct_proteins_df["NumOfReads"] = tmr1000_max_distinct_proteins_df.apply(
#     lambda x: x["NumOfReads"] if not pd.isna(x["NumOfReads"]) else x["MappedReads"],
#     axis=1,
# )

# # tmr1000_max_distinct_proteins_df = (
# #     tmr1000_max_distinct_proteins_df.dropna().reset_index(drop=True)
# # )

# tmr1000_max_distinct_proteins_df["DistinctProteins/Reads"] = (
#     tmr1000_max_distinct_proteins_df["NumOfProteins"]
#     / tmr1000_max_distinct_proteins_df["NumOfReads"]
# )

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

# tmr1000_max_distinct_proteins_df = tmr1000_max_distinct_proteins_df.merge(
#     neuronal_reads_per_chrom_df.loc[
#         :, ["Chrom", "ChromNeuronality", "ChromNeuronality2"]
#     ].rename(
#         columns={"ChromNeuronality": "IsNeuralNew", "ChromNeuronality2": "IsNeuralNew2"}
#     ),
#     on="Chrom",
#     how="left",
# )


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
# number of chroms with at least 5 distinct proteins, considering chroms with editing (below 6% noise)
tmr1000_max_distinct_proteins_df.loc[
    (tmr1000_max_distinct_proteins_df["NumOfProteins"] >= 5)
    & (tmr1000_max_distinct_proteins_df["Chrom"].isin(tmr1000_chroms))
].shape[0]

# %%
# % of chroms with at least 5 distinct proteins, considering chroms with editing (below 6% noise)
100 * tmr1000_max_distinct_proteins_df.loc[
    (tmr1000_max_distinct_proteins_df["NumOfProteins"] >= 5)
    & (tmr1000_max_distinct_proteins_df["Chrom"].isin(tmr1000_chroms))
].shape[0] / len(tmr1000_chroms)

# %%
# number of chroms with at least 50 distinct proteins, considering chroms with editing (below 6% noise)
tmr1000_max_distinct_proteins_df.loc[
    (tmr1000_max_distinct_proteins_df["NumOfProteins"] >= 50)
    & (tmr1000_max_distinct_proteins_df["Chrom"].isin(tmr1000_chroms))
].shape[0]

# %%
# % of chroms with at least 50 distinct proteins, considering chroms with editing (below 6% noise)
100 * tmr1000_max_distinct_proteins_df.loc[
    (tmr1000_max_distinct_proteins_df["NumOfProteins"] >= 50)
    & (tmr1000_max_distinct_proteins_df["Chrom"].isin(tmr1000_chroms))
].shape[0] / len(tmr1000_chroms)

# %%

# %%
raise Error("don't run the notebook from this point")

# %% [markdown]
# ### Distinct unique proteins - neuronal

# %%
assert (
    len(neuronal_conditions)
    == len(neuronal_chroms)
    == len(neuronal_distinct_unique_proteins_files)
    == len(neuronal_unique_reads_dfs)
)

neuronal_distinct_unique_proteins_dfs = []
for condition, chrom, distinct_unique_proteins_file, unique_reads_df in zip(
    neuronal_conditions,
    neuronal_chroms,
    neuronal_distinct_proteins_files,
    neuronal_unique_reads_dfs,
):
    neuronal_distinct_unique_proteins_df = pd.read_csv(
        distinct_unique_proteins_file, sep=sep, dtype={"AvailableReads": str}
    )
    neuronal_distinct_unique_proteins_df.insert(0, condition_col, condition)
    neuronal_distinct_unique_proteins_df.insert(
        1,
        "NumOfReads",
        (
            neuronal_distinct_unique_proteins_df["Fraction"]
            * unique_reads_df["NumOfReads"].sum()
        ).astype(int),
    )
    neuronal_distinct_unique_proteins_df.insert(0, "Chrom", chrom)
    neuronal_distinct_unique_proteins_dfs.append(neuronal_distinct_unique_proteins_df)

ic(len(neuronal_distinct_unique_proteins_dfs))

neuronal_distinct_unique_proteins_df = (
    pd.concat(neuronal_distinct_unique_proteins_dfs)
    .reset_index(drop=True)
    .rename(columns={"NumUniqueSamples": "NumOfProteins", "UniqueSamples": "Proteins"})
)

neuronal_distinct_unique_proteins_df = neuronal_distinct_unique_proteins_df.sort_values(
    [
        condition_col,
        "Fraction",
        "FractionRepetition",
        "Algorithm",
        "AlgorithmRepetition",
    ]
).reset_index(drop=True)

neuronal_distinct_unique_proteins_df

# %% [markdown]
# #### Max distinct per fraction - neuronal

# %%
# create a temp col which will soon be deleted
neuronal_distinct_unique_proteins_df["TempIndex"] = (
    neuronal_distinct_unique_proteins_df.index
)

# first, create a df with all the largset solutions per chrom and fraction
# (this happens due to n=1, keep="all" in the nlargst function)
neuronal_max_distinct_per_fraction_df = (
    neuronal_distinct_unique_proteins_df.groupby(["Chrom", "Fraction"]).apply(
        pd.DataFrame.nlargest,
        n=1,
        keep="all",
        columns="NumOfProteins",
        include_groups=False,
    )
    # .reset_index(drop=True)
)
# use merge to also get the chrom and fraction cols lost due to include_groups=False
neuronal_max_distinct_per_fraction_df = neuronal_distinct_unique_proteins_df.loc[
    neuronal_distinct_unique_proteins_df["TempIndex"].isin(
        neuronal_max_distinct_per_fraction_df["TempIndex"]
    )
].drop(columns="TempIndex")

del neuronal_distinct_unique_proteins_df["TempIndex"]

# then, sample a single largest solution (a set of distinct proteins)
# per each fraction in each chrom
neuronal_max_distinct_per_fraction_df = (
    neuronal_max_distinct_per_fraction_df.groupby(["Chrom", "Fraction"])
    .sample(n=1, random_state=seed)
    .reset_index(drop=True)
)

neuronal_max_distinct_per_fraction_df

# %%
neuronal_max_distinct_per_fraction_df.groupby(["Chrom", "Fraction"]).size().unique()

# %% [markdown]
# #### Max distinct - neuronal

# %%
neuronal_max_distinct_proteins_df = neuronal_max_distinct_per_fraction_df.loc[
    neuronal_max_distinct_per_fraction_df["Fraction"] == 1
]


neuronal_max_distinct_proteins_df = neuronal_max_distinct_proteins_df.merge(
    neuronal_alignment_stats_df,
    on="Chrom",
    # how="left",
    how="right",
)

neuronal_max_distinct_proteins_df["NumOfProteins"] = neuronal_max_distinct_proteins_df[
    "NumOfProteins"
].fillna(1)


# neuronal_max_distinct_proteins_df = neuronal_max_distinct_proteins_df.merge(
#     neural_vs_non_neural_expression_df.loc[:, ["OvulChrom", "IsNeural"]].rename(
#         columns={"OvulChrom": "Chrom"}
#     ),
#     on="Chrom",
#     how="left",
# )
# tmr1000_max_distinct_proteins_df["IsNeural"] = tmr1000_max_distinct_proteins_df[
#     "IsNeural"
# ].fillna("Missing")

# tmr1000_max_distinct_proteins_df = tmr1000_max_distinct_proteins_df.merge(
#     neuronal_reads_per_chrom_df.loc[
#         :, ["Chrom", "ChromNeuronality", "ChromNeuronality2"]
#     ].rename(
#         columns={"ChromNeuronality": "IsNeuralNew", "ChromNeuronality2": "IsNeuralNew2"}
#     ),
#     on="Chrom",
#     how="left",
# )


neuronal_max_distinct_proteins_df = neuronal_max_distinct_proteins_df.sort_values(
    "NumOfProteins", ascending=False, ignore_index=True
)

neuronal_max_distinct_proteins_df

# %%
# number of chroms with at least 5 distinct proteins
neuronal_max_distinct_proteins_df.loc[
    neuronal_max_distinct_proteins_df["NumOfProteins"] >= 5
].shape[0]

# %%
# % of chroms with at least 5 distinct proteins
100 * neuronal_max_distinct_proteins_df.loc[
    neuronal_max_distinct_proteins_df["NumOfProteins"] >= 5
].shape[0] / neuronal_max_distinct_proteins_df.shape[0]

# %%
# number of chroms with at least 50 distinct proteins
neuronal_max_distinct_proteins_df.loc[
    neuronal_max_distinct_proteins_df["NumOfProteins"] >= 50
].shape[0]

# %%
# % of chroms with at least 50 distinct proteins
100 * neuronal_max_distinct_proteins_df.loc[
    neuronal_max_distinct_proteins_df["NumOfProteins"] >= 50
].shape[0] / neuronal_max_distinct_proteins_df.shape[0]

# %%
# number of chroms with at least 5 distinct proteins, considering chroms with editing (below 6% noise)
neuronal_max_distinct_proteins_df.loc[
    (neuronal_max_distinct_proteins_df["NumOfProteins"] >= 5)
    & (neuronal_max_distinct_proteins_df["Chrom"].isin(neuronal_chroms))
].shape[0]

# %%
# % of chroms with at least 5 distinct proteins, considering chroms with editing (below 6% noise)
100 * neuronal_max_distinct_proteins_df.loc[
    (neuronal_max_distinct_proteins_df["NumOfProteins"] >= 5)
    & (neuronal_max_distinct_proteins_df["Chrom"].isin(neuronal_chroms))
].shape[0] / len(neuronal_chroms)

# %%
# number of chroms with at least 50 distinct proteins, considering chroms with editing (below 6% noise)
neuronal_max_distinct_proteins_df.loc[
    (neuronal_max_distinct_proteins_df["NumOfProteins"] >= 50)
    & (neuronal_max_distinct_proteins_df["Chrom"].isin(neuronal_chroms))
].shape[0]

# %%
# % of chroms with at least 50 distinct proteins, considering chroms with editing (below 6% noise)
100 * neuronal_max_distinct_proteins_df.loc[
    (neuronal_max_distinct_proteins_df["NumOfProteins"] >= 50)
    & (neuronal_max_distinct_proteins_df["Chrom"].isin(neuronal_chroms))
].shape[0] / len(neuronal_chroms)

# %% [markdown]
# ### Distinct unique proteins - non-neuronal

# %%
assert (
    len(non_neuronal_conditions)
    == len(non_neuronal_chroms)
    == len(non_neuronal_distinct_unique_proteins_files)
    == len(non_neuronal_unique_reads_dfs)
)

non_neuronal_distinct_unique_proteins_dfs = []
for condition, chrom, distinct_unique_proteins_file, unique_reads_df in zip(
    non_neuronal_conditions,
    non_neuronal_chroms,
    non_neuronal_distinct_proteins_files,
    non_neuronal_unique_reads_dfs,
):
    non_neuronal_distinct_unique_proteins_df = pd.read_csv(
        distinct_unique_proteins_file, sep=sep, dtype={"AvailableReads": str}
    )
    non_neuronal_distinct_unique_proteins_df.insert(0, condition_col, condition)
    non_neuronal_distinct_unique_proteins_df.insert(
        1,
        "NumOfReads",
        (
            non_neuronal_distinct_unique_proteins_df["Fraction"]
            * unique_reads_df["NumOfReads"].sum()
        ).astype(int),
    )
    non_neuronal_distinct_unique_proteins_df.insert(0, "Chrom", chrom)
    non_neuronal_distinct_unique_proteins_dfs.append(
        non_neuronal_distinct_unique_proteins_df
    )

ic(len(non_neuronal_distinct_unique_proteins_dfs))

non_neuronal_distinct_unique_proteins_df = (
    pd.concat(non_neuronal_distinct_unique_proteins_dfs)
    .reset_index(drop=True)
    .rename(columns={"NumUniqueSamples": "NumOfProteins", "UniqueSamples": "Proteins"})
)

non_neuronal_distinct_unique_proteins_df = (
    non_neuronal_distinct_unique_proteins_df.sort_values(
        [
            condition_col,
            "Fraction",
            "FractionRepetition",
            "Algorithm",
            "AlgorithmRepetition",
        ]
    ).reset_index(drop=True)
)

non_neuronal_distinct_unique_proteins_df

# %% [markdown]
# #### Max distinct per fraction - non-neuronal

# %%
# create a temp col which will soon be deleted
non_neuronal_distinct_unique_proteins_df["TempIndex"] = (
    non_neuronal_distinct_unique_proteins_df.index
)

# first, create a df with all the largset solutions per chrom and fraction
# (this happens due to n=1, keep="all" in the nlargst function)
non_neuronal_max_distinct_per_fraction_df = (
    non_neuronal_distinct_unique_proteins_df.groupby(["Chrom", "Fraction"]).apply(
        pd.DataFrame.nlargest,
        n=1,
        keep="all",
        columns="NumOfProteins",
        include_groups=False,
    )
    # .reset_index(drop=True)
)
# use merge to also get the chrom and fraction cols lost due to include_groups=False
non_neuronal_max_distinct_per_fraction_df = (
    non_neuronal_distinct_unique_proteins_df.loc[
        non_neuronal_distinct_unique_proteins_df["TempIndex"].isin(
            non_neuronal_max_distinct_per_fraction_df["TempIndex"]
        )
    ].drop(columns="TempIndex")
)

del non_neuronal_distinct_unique_proteins_df["TempIndex"]

# then, sample a single largest solution (a set of distinct proteins)
# per each fraction in each chrom
non_neuronal_max_distinct_per_fraction_df = (
    non_neuronal_max_distinct_per_fraction_df.groupby(["Chrom", "Fraction"])
    .sample(n=1, random_state=seed)
    .reset_index(drop=True)
)

non_neuronal_max_distinct_per_fraction_df

# %%
non_neuronal_max_distinct_per_fraction_df.groupby(["Chrom", "Fraction"]).size().unique()

# %% [markdown]
# #### Max distinct - non-neuronal

# %%
non_neuronal_max_distinct_proteins_df = non_neuronal_max_distinct_per_fraction_df.loc[
    non_neuronal_max_distinct_per_fraction_df["Fraction"] == 1
]


non_neuronal_max_distinct_proteins_df = non_neuronal_max_distinct_proteins_df.merge(
    non_neuronal_alignment_stats_df,
    on="Chrom",
    # how="left",
    how="right",
)

non_neuronal_max_distinct_proteins_df["NumOfProteins"] = (
    non_neuronal_max_distinct_proteins_df["NumOfProteins"].fillna(1)
)


# neuronal_max_distinct_proteins_df = neuronal_max_distinct_proteins_df.merge(
#     neural_vs_non_neural_expression_df.loc[:, ["OvulChrom", "IsNeural"]].rename(
#         columns={"OvulChrom": "Chrom"}
#     ),
#     on="Chrom",
#     how="left",
# )
# tmr1000_max_distinct_proteins_df["IsNeural"] = tmr1000_max_distinct_proteins_df[
#     "IsNeural"
# ].fillna("Missing")

# tmr1000_max_distinct_proteins_df = tmr1000_max_distinct_proteins_df.merge(
#     neuronal_reads_per_chrom_df.loc[
#         :, ["Chrom", "ChromNeuronality", "ChromNeuronality2"]
#     ].rename(
#         columns={"ChromNeuronality": "IsNeuralNew", "ChromNeuronality2": "IsNeuralNew2"}
#     ),
#     on="Chrom",
#     how="left",
# )


non_neuronal_max_distinct_proteins_df = (
    non_neuronal_max_distinct_proteins_df.sort_values(
        "NumOfProteins", ascending=False, ignore_index=True
    )
)

non_neuronal_max_distinct_proteins_df

# %%
# number of chroms with at least 5 distinct proteins
non_neuronal_max_distinct_proteins_df.loc[
    non_neuronal_max_distinct_proteins_df["NumOfProteins"] >= 5
].shape[0]

# %%
# % of chroms with at least 5 distinct proteins
100 * non_neuronal_max_distinct_proteins_df.loc[
    non_neuronal_max_distinct_proteins_df["NumOfProteins"] >= 5
].shape[0] / non_neuronal_max_distinct_proteins_df.shape[0]

# %%
# number of chroms with at least 50 distinct proteins
non_neuronal_max_distinct_proteins_df.loc[
    non_neuronal_max_distinct_proteins_df["NumOfProteins"] >= 50
].shape[0]

# %%
# % of chroms with at least 50 distinct proteins
100 * non_neuronal_max_distinct_proteins_df.loc[
    non_neuronal_max_distinct_proteins_df["NumOfProteins"] >= 50
].shape[0] / non_neuronal_max_distinct_proteins_df.shape[0]

# %%
# number of chroms with at least 5 distinct proteins, considering chroms with editing (below 6% noise)
non_neuronal_max_distinct_proteins_df.loc[
    (non_neuronal_max_distinct_proteins_df["NumOfProteins"] >= 5)
    & (non_neuronal_max_distinct_proteins_df["Chrom"].isin(non_neuronal_chroms))
].shape[0]

# %%
# % of chroms with at least 5 distinct proteins, considering chroms with editing (below 6% noise)
100 * non_neuronal_max_distinct_proteins_df.loc[
    (non_neuronal_max_distinct_proteins_df["NumOfProteins"] >= 5)
    & (non_neuronal_max_distinct_proteins_df["Chrom"].isin(non_neuronal_chroms))
].shape[0] / len(non_neuronal_chroms)

# %%
# number of chroms with at least 50 distinct proteins, considering chroms with editing (below 6% noise)
non_neuronal_max_distinct_proteins_df.loc[
    (non_neuronal_max_distinct_proteins_df["NumOfProteins"] >= 50)
    & (non_neuronal_max_distinct_proteins_df["Chrom"].isin(non_neuronal_chroms))
].shape[0]

# %%
# % of chroms with at least 50 distinct proteins, considering chroms with editing (below 6% noise)
100 * non_neuronal_max_distinct_proteins_df.loc[
    (non_neuronal_max_distinct_proteins_df["NumOfProteins"] >= 50)
    & (non_neuronal_max_distinct_proteins_df["Chrom"].isin(non_neuronal_chroms))
].shape[0] / len(non_neuronal_chroms)

# %% [markdown] papermill={"duration": 0.045853, "end_time": "2022-02-01T09:42:48.953594", "exception": false, "start_time": "2022-02-01T09:42:48.907741", "status": "completed"}
# # Results

# %% [markdown] jp-MarkdownHeadingCollapsed=true papermill={"duration": 0.124528, "end_time": "2022-02-01T09:43:10.054394", "exception": false, "start_time": "2022-02-01T09:43:09.929866", "status": "completed"}
# ## Positions

# %% [markdown] papermill={"duration": 0.149848, "end_time": "2022-02-01T09:43:12.800733", "exception": false, "start_time": "2022-02-01T09:43:12.650885", "status": "completed"}
# ### Coverage - per-transcript per-sample

# %%
concat_all_positions_df


# %%
def calc_per_transcript_per_sample_coverage(
    positions_df,
    samples_and_tissues_df,
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
def calc_per_transcript_per_sample_coverage_dfs(
    concat_all_positions_df,
    possibly_na_positions_files,
    possibly_na_chroms,
    samples_and_tissues_df,
    processes=4,
):
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
    concat_all_positions_df,
    possibly_na_positions_files,
    possibly_na_chroms,
    samples_and_tissues_df,
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


# %% [markdown] papermill={"duration": 0.149848, "end_time": "2022-02-01T09:43:12.800733", "exception": false, "start_time": "2022-02-01T09:43:12.650885", "status": "completed"}
# ### Editing index - per transcript

# %%
def editing_index_per_transcript(positions_df, strand):
    ref_base = "A" if strand == "+" else "T"
    alt_base = "G" if strand == "+" else "C"
    all_refbase_positions_df = positions_df.loc[
        (positions_df["RefBase"] == ref_base)
        & (~positions_df["InProbRegion"])
        & (positions_df["CDS"])
    ]
    num_of_all_editable_adenosines = all_refbase_positions_df["TotalCoverage"].sum()
    num_of_edited_adenosines = all_refbase_positions_df[alt_base].sum()
    editing_index = 100 * num_of_edited_adenosines / num_of_all_editable_adenosines
    return editing_index


def make_one_input(concat_all_positions_df, chrom, strand):
    return (
        concat_all_positions_df.loc[concat_all_positions_df["Chrom"] == chrom],
        strand,
    )


def make_editing_index_per_transcript_inputs(concat_all_positions_df, chroms, strands):
    concat_all_positions_df = concat_all_positions_df.loc[
        concat_all_positions_df["Chrom"].isin(chroms)
    ]

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


def calc_per_transcript_editing_index_df(
    concat_all_positions_df, chroms, strands, processes=4
):
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
        {"Chrom": chroms, "EditingIndex": per_transcript_editing_indices}
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
            (positions_df["RefBase"] == ref_base)
            & (~positions_df["InProbRegion"])
            & (positions_df["CDS"])
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
    expanded_all_refbase_positions_df["MappedBases"] = (
        expanded_all_refbase_positions_df["MappedBases"].apply(list)
    )

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
def calc_per_sample_editing_index_df(
    concat_all_positions_df, chroms, strands, samples, processes=1
):
    concat_all_positions_df = concat_all_positions_df.loc[
        concat_all_positions_df["Chrom"].isin(chroms)
    ]

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
        & (concat_all_positions_df["NoisyFinal"]),
    ]
    .groupby("Chrom")
    .size()
    .reset_index()
    .rename(columns={0: "NoisePositions"})["NoisePositions"]
    .describe()
)


# %%
def mean_noise_levels(positions_df, top_x_noisy_positions=3, snp_noise_level=0.1):
    # if positions_df.empty:
    #     return 0.0
    noise_levels = (
        positions_df.loc[
            (positions_df["Noise"] <= snp_noise_level) & (positions_df["NoisyFinal"]),
            "Noise",
        ]
        .sort_values(ascending=False)[:top_x_noisy_positions]
        .tolist()
    )
    # if there are less noisy positions than `top_x_noisy_positions`, add zeros accordingly
    noise_levels = pd.Series(
        noise_levels + [0 for _ in range(top_x_noisy_positions - len(noise_levels))]
    )
    return noise_levels.mean()


# %%
all_per_chrom_mean_noise_levels = (
    concat_all_positions_df
    # .groupby(["Transcript", "Chrom"])
    .groupby("Chrom")
    .apply(mean_noise_levels)
    .reset_index()
    .rename(columns={0: "Noise"})
    .merge(tmr50_alignment_stats_df.loc[:, ["Chrom"]], on="Chrom", how="right")
    .fillna(0.0)
    .sort_values(["Chrom", "Noise"])
)
all_per_chrom_mean_noise_levels["%Noise"] = (
    100 * all_per_chrom_mean_noise_levels["Noise"]
)
all_per_chrom_mean_noise_levels

# %%
saved_per_chrom_mean_noise_levels_df = all_per_chrom_mean_noise_levels.merge(
    orfs_df[["Chrom", "Name"]], how="left"
)
saved_per_chrom_mean_noise_levels_df.insert(
    1, "Gene", saved_per_chrom_mean_noise_levels_df["Name"]
)
del saved_per_chrom_mean_noise_levels_df["Name"]
saved_per_chrom_mean_noise_levels_df.insert(
    0, "Platform", "Whole-transcriptome octopus data"
)
saved_per_chrom_mean_noise_levels_df.to_csv("Noise.Octopus.tsv", sep="\t", index=False)
saved_per_chrom_mean_noise_levels_df

# %%
# per_chrom_mean_noise_levels["%Noise"].describe()

# %%
# described_noise_df = per_chrom_mean_noise_levels["%Noise"].describe()
# quartiles = ["25%", "50%", "75%"]
# noise_quartiles = described_noise_df.loc[quartiles].values
# noise_quartiles

# %%
all_described_noise_df = all_per_chrom_mean_noise_levels["%Noise"].describe()
all_described_noise_df.round(1)

# %%
round(1.2)

# %%
round(scipy.stats.iqr(all_per_chrom_mean_noise_levels["%Noise"]), 1)

# %%
np.percentile(all_per_chrom_mean_noise_levels["%Noise"], [25, 75])

# %%
np.round(np.percentile(all_per_chrom_mean_noise_levels["%Noise"], [25, 75]), 1)

# %%
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
    opacity=0.5,
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
        ic(all_noise_quartiles[i], all_noise_quartiles[i + 1])
        if np.isclose(all_noise_quartiles[i], all_noise_quartiles[i + 1]):
            continue
    except IndexError:
        # if this is the last line to plot
        pass
    ic(noise_quartile)
    fig.add_shape(
        type="line",
        x0=noise_quartile,
        x1=noise_quartile,
        y0=0,
        # y1=max_noise_quartiles_y - (i * 0.3 * max_noise_quartiles_y),
        y1=ic(max_noise_quartiles_y - (i * 0.3 * max_noise_quartiles_y)),
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


# fig.write_image(
#     "Mean per chrom noise levels - with quartiles - Octopus.svg",
#     width=width,
#     height=height,
# )

fig.show()


# %% [markdown]
# ### Mismatches by type

# %%
def find_alt_base(
    ref_base: str, a_count: int, t_count: int, c_count: int, g_count: int
):
    """
    Find the base with most supporting reads other than `ref_base`.
    If there are two or more such bases, the function picks one at random.
    """
    base_counts_dict = {"A": a_count, "T": t_count, "C": c_count, "G": g_count}
    alt_bases = set(base_counts_dict) - {ref_base}
    alt_base_counts_dict = {
        base: base_count
        for base, base_count in base_counts_dict.items()
        if base in alt_bases
    }
    max_alt_base_count = max(alt_base_counts_dict.values())
    max_alt_bases = [
        base
        for base, base_count in alt_base_counts_dict.items()
        if base_count == max_alt_base_count
    ]
    alt_base = choice(max_alt_bases)
    return alt_base


# %%
def define_mismatch_type(ref_base: str, alt_base: str, strand: str):
    if strand == "-":
        opposing_bases = {"A": "T", "T": "A", "C": "G", "G": "C"}
        ref_base = opposing_bases[ref_base]
        alt_base = opposing_bases[alt_base]
    mismatch_type = f"{ref_base}>{alt_base}"
    return mismatch_type


# %%
def mismatch_frequency(ref_base_count, alt_base_count):
    return alt_base_count / (ref_base_count + alt_base_count)


# %%
mismatches_df = concat_all_positions_df.loc[
    :,
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
        "Edited",
        # "EditingBinomPVal",
        # "EditingCorrectedPVal",
        "EditedCorrected",
        "EditedFinal",
        "Noise",
        "NoisyCorrected",
        "NoisyFinal",
    ],
]

mismatches_df.insert(
    # mismatches_df.columns.get_loc("G") + 1,
    mismatches_df.columns.get_loc("RefBase") + 1,
    "AltBase",
    mismatches_df.apply(
        lambda x: find_alt_base(
            x["RefBase"],
            x["A"],
            x["T"],
            x["C"],
            x["G"],
        ),
        axis=1,
    ),
)

# no point in keeping positions without mismatches
mismatches_df.insert(
    # mismatches_df.columns.get_loc("AltBase") + 1,
    mismatches_df.columns.get_loc("G") + 1,
    "AltBaseCount",
    mismatches_df.apply(
        lambda x: x[x["AltBase"]],
        axis=1,
    ),
)
mismatches_df = mismatches_df.loc[mismatches_df["AltBaseCount"] > 0]
# del mismatches_df["AltBaseCount"]

# annotate mismatch frequency as an alternative and general annotation to `EditingFrequency` and `Noise`
mismatches_df.insert(
    mismatches_df.columns.get_loc("G") + 1,
    "MismatchFrequency",
    mismatches_df.apply(
        lambda x: mismatch_frequency(x[x["RefBase"]], x[x["AltBase"]]), axis=1
    ),
)

# get strand information to define mismatch type on the minus strand
mismatches_df = mismatches_df.merge(data_df.loc[:, ["Chrom", "Strand"]], how="left")
mismatches_df.insert(
    mismatches_df.columns.get_loc("Chrom") + 1, "Strand2", mismatches_df["Strand"]
)
del mismatches_df["Strand"]
mismatches_df = mismatches_df.rename(columns={"Strand2": "Strand"})

mismatches_df.insert(
    mismatches_df.columns.get_loc("AltBase") + 1,
    "Mismatch",
    mismatches_df.apply(
        lambda x: define_mismatch_type(x["RefBase"], x["AltBase"], x["Strand"]), axis=1
    ),
)

mismatches_df

# %%
true_a2g_mismatches = mismatches_df.loc[
    (mismatches_df["EditedFinal"])
    & (mismatches_df["Mismatch"] == "A>G")
    & (mismatches_df["Chrom"].isin(chroms))
]
true_a2g_mismatches

# %%
wrong_a2g_mismatches = mismatches_df.loc[
    (mismatches_df["EditedFinal"])
    & (mismatches_df["Mismatch"] != "A>G")
    & (mismatches_df["Chrom"].isin(chroms))
]
wrong_a2g_mismatches

# %%
len(wrong_a2g_mismatches) / len(true_a2g_mismatches)

# %%
wrong_a2g_mismatches_per_chrom = (
    wrong_a2g_mismatches.groupby(["Chrom", condition_col])
    .size()
    .reset_index()
    .rename(columns={0: "WrongEditingSitesPerChrom"})
    .sort_values("WrongEditingSitesPerChrom", ascending=False)
)
wrong_a2g_mismatches_per_chrom

# %%
wrong_a2g_mismatches_per_chrom["WrongEditingSitesPerChrom"].describe()

# %%
significant_mismatches_df = mismatches_df.loc[
    (mismatches_df["Chrom"].isin(chroms))
    & ((mismatches_df["EditedFinal"]) | (mismatches_df["NoisyFinal"]))
]
# significant_mismatches_df["MismatchFrequency>10%"] = significant_mismatches_df["MismatchFrequency"] > 0.1
significant_mismatches_df

# %%
significant_mismatches_df.loc[
    (significant_mismatches_df["EditingFrequency"] == 1)
    | (significant_mismatches_df["Noise"] == 1)
]

# %%
# mismatches_color_sequence = px.colors.qualitative.Set3
mismatches_color_sequence = px.colors.qualitative.Dark24
mismatch_dolor_map = {
    mismatch: color
    for mismatch, color in zip(
        significant_mismatches_df["Mismatch"].unique(), mismatches_color_sequence
    )
}
# mismatch_dolor_map

# %%
fig = px.histogram(
    significant_mismatches_df,
    x="Mismatch",
    # x="MismatchFrequency",
    # facet_col="EditedFinal",
    color="Mismatch",
    color_discrete_map=mismatch_dolor_map,
    # facet_col_wrap=4,
    # log_y=True,
    template=template,
    title="Number of significant mismatches in transcripts with pooled noise level < 6%",
)

# Reduce opacity to see both histograms
# fig.update_traces(opacity=0.75)
fig.update_layout(
    width=600,
    height=500,
    showlegend=False,
    # barmode='overlay' # Overlay both histograms
)

fig.show()

# %%
fig = px.histogram(
    significant_mismatches_df,
    x="MismatchFrequency",
    facet_col="Mismatch",
    color="Mismatch",
    color_discrete_map=mismatch_dolor_map,
    facet_col_wrap=4,
    # log_y=True,
    template=template,
    title="Cummulative relative mismatch frequency distribution of significant mismatches in transcripts with pooled noise level < 6%",
    cumulative=True,
    histnorm="percent",
)

fig.update_layout(width=1600, height=800, showlegend=False)

fig.show()

# %%
fig = px.scatter(
    significant_mismatches_df.loc[
        (significant_mismatches_df["EditedFinal"])
        & (significant_mismatches_df["Mismatch"] != "A>G")
    ],
    x="EditingFrequency",
    y="MismatchFrequency",
    facet_col="Mismatch",
    color="Mismatch",
    color_discrete_map=mismatch_dolor_map,
    # facet_col_wrap=4,
    # log_y=True,
    template=template,
    title="Editing frequency vs mismatch frequency of `wrong` significantly editing sites<br>in transcripts with pooled noise level < 6%",
)

fig.update_traces(marker_size=7)
fig.update_layout(width=800, height=400, showlegend=False)

fig.show()

# %%
fig = px.histogram(
    wrong_a2g_mismatches_per_chrom,
    x="WrongEditingSitesPerChrom",
    # x="MismatchFrequency",
    # facet_col="EditedFinal",
    # color="Mismatch",
    # color_discrete_map=mismatch_dolor_map,
    # facet_col_wrap=4,
    log_y=True,
    template=template,
    title="Number of wrong A>G mismatches per transcripts with pooled noise<br>level < 6%",
)
fig.update_xaxes(dtick=1)
fig.update_layout(
    width=700,
    height=400,
    showlegend=False,
    # barmode='overlay' # Overlay both histograms
)

fig.show()

# %% [markdown]
# ### Machine noise

# %%
all_machine_noise_df = concat_all_positions_df.loc[
    :, [condition_col, "Chrom", "RefBase", "TotalCoverage", "A", "T", "C", "G"]
]
all_machine_noise_df["ATCGs"] = all_machine_noise_df.loc[:, ["A", "T", "C", "G"]].sum(
    axis=1
)
all_machine_noise_df["Matches"] = all_machine_noise_df.apply(
    lambda x: x[x["RefBase"]], axis=1
)
all_machine_noise_df["Mismatches"] = all_machine_noise_df.apply(
    lambda x: x["ATCGs"] - x["Matches"], axis=1
)
all_machine_noise_df

# %%
all_pooled_per_chrom_machine_noise_df = (
    all_machine_noise_df.groupby("Chrom")[["Matches", "Mismatches"]].sum().reset_index()
)
all_pooled_per_chrom_machine_noise_df["%PooledMachineNoise"] = (
    all_pooled_per_chrom_machine_noise_df.apply(
        lambda x: 100 * x["Mismatches"] / x["Matches"], axis=1
    )
)
all_pooled_per_chrom_machine_noise_df = all_pooled_per_chrom_machine_noise_df.merge(
    tmr50_alignment_stats_df.loc[:, ["Chrom"]], on="Chrom", how="right"
).fillna(0.0)
all_pooled_per_chrom_machine_noise_df

# %%
fig = px.histogram(
    all_pooled_per_chrom_machine_noise_df,
    x="%PooledMachineNoise",
    color_discrete_sequence=["black"],
    # labels={"% noise": "Per-gene noise level [%]"},
    log_y=True,
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
            (concat_all_positions_df[label]) & (concat_all_positions_df["CDS"])
            # only transcripts whose pooled noise levels is < 6%
            & (concat_all_positions_df["Chrom"].isin(chroms)),
            [condition_col, "Chrom", "Position"],
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
    "Single-cell octopus data",
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

# plt.savefig("Known vs new editing sites - Octopus.svg", format="svg", dpi=300)

plt.show()

# %% [markdown]
# Finding the number of genes/transcripts with no known editing sites

# %%
all_chroms = concat_all_positions_df["Chrom"].unique().size
chroms_with_new_sites = (
    concat_all_positions_df.loc[
        (concat_all_positions_df["EditedFinal"])
        & (~concat_all_positions_df["KnownEditing"]),
        "Chrom",
    ]
    .unique()
    .size
)
chroms_without_new_sites = all_chroms - chroms_with_new_sites
chroms_without_new_sites

# %%
f"{100 * chroms_without_new_sites / all_chroms:.2f}%"

# %%
(
    concat_all_positions_df.loc[
        (concat_all_positions_df["EditedFinal"])
        & (~concat_all_positions_df["KnownEditing"])
        & (concat_all_positions_df["Chrom"].isin(chroms)),
        "Chrom",
    ]
    .unique()
    .size
)

# %%
num_of_chroms_with_editing_sites = len(chroms)

num_of_chroms_with_new_sites = (
    concat_all_positions_df.loc[
        (concat_all_positions_df["EditedFinal"])
        & (~concat_all_positions_df["KnownEditing"])
        & (concat_all_positions_df["Chrom"].isin(chroms)),
        "Chrom",
    ]
    .unique()
    .size
)

num_of_chroms_without_new_sites = (
    num_of_chroms_with_editing_sites - num_of_chroms_with_new_sites
)
num_of_chroms_without_new_sites

# %%
100 * num_of_chroms_without_new_sites / num_of_chroms_with_editing_sites

# %% [markdown]
# ### ADAR motif

# %%
unique_positions_df = (
    concat_all_positions_df.loc[
        # only edited positions in transcripts whose pooled noise levels is < 6%
        (concat_all_positions_df["EditedFinal"])
        & (concat_all_positions_df["Chrom"].isin(chroms)),
        ["Chrom", "Position"],
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
main_title = "Single-cell octopus data"
sub_titles = [""]

# out_file = Path("ADAR motif of pooled editing sites - Octopus - Single-cell.svg")
out_file = None

# %%
multiple_logos_from_fasta_files(
    fasta_files, main_title, sub_titles, out_file, width=0.33 * 14, height=4, dpi=300
);

# %% [markdown] papermill={"duration": 0.030615, "end_time": "2022-02-01T09:42:49.024262", "exception": false, "start_time": "2022-02-01T09:42:48.993647", "status": "completed"}
# ## Num of distinct proteins

# %% [markdown]
# ### Pooled

# %%
# # % of chroms with at least 50 distinct proteins, considering chroms with editing (below 6% noise)
# 100 * tmr1000_max_distinct_proteins_df.loc[
#     (tmr1000_max_distinct_proteins_df["NumOfProteins"] >= 50)
#     & (tmr1000_max_distinct_proteins_df["Chrom"].isin(tmr1000_chroms))
# ].shape[0] / len(tmr1000_chroms)

# %%
x = max_distinct_proteins_df.loc[
    max_distinct_proteins_df["IsNeural"] == "Yes", "NumOfProteins"
]
y = max_distinct_proteins_df.loc[
    max_distinct_proteins_df["IsNeural"] == "No", "NumOfProteins"
]
statistic, pv = scipy.stats.mannwhitneyu(x, y)
statistic, pv

# %%
max_distinct_proteins_df

# %%
neuronal_max_distinct_proteins_df

# %%
non_neuronal_max_distinct_proteins_df

# %%
x = neuronal_max_distinct_proteins_df["NumOfProteins"]
y = non_neuronal_max_distinct_proteins_df["NumOfProteins"]
statistic, pv = scipy.stats.mannwhitneyu(x, y)
statistic, pv

# %%

# %%
max_distinct_proteins_df.loc[
    (max_distinct_proteins_df["IsNeuralNew"] == "Yes")
    & (max_distinct_proteins_df["IsNeural"] == "Yes")
]

# %%
x = max_distinct_proteins_df.loc[
    max_distinct_proteins_df["IsNeuralNew"] == "Yes", "NumOfProteins"
]
y = max_distinct_proteins_df.loc[
    max_distinct_proteins_df["IsNeuralNew"] == "No", "NumOfProteins"
]
statistic, pv = scipy.stats.mannwhitneyu(x, y)
statistic, pv

# %%

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

# %%
x

# %%
x_mean = x.mean()
x_mean

# %%
x_std = x.std()
x_std

# %%
# rank the closeness of x values to the mean
(x - x_mean).abs()

# %%
# rank the closeness of x values to the mean
(x - x_mean).abs().sort_values()

# %%
(x - x_mean).abs().argsort()

# %%
(x - x_mean).abs().argsort()[:1]

# %%
x_mean_closest = x.iloc[(x - x_mean).abs().argsort()[:1]]
x_mean_closest

# %%
x_mean_closest_k = x_mean_closest.index.values[0]
x_mean_closest_k

# %%
x_mean_closest.values[0]

# %%
x_mean == x_mean_closest.values[0]

# %%

# %%
x_mean < x_mean_closest.values[0]

# %%
if x_mean == x_mean_closest.values[0]:
    y_mean = y.iloc[x_mean_closest_k]
else:
    if x_mean < x_mean_closest.values[0]:
        i = x_mean_closest_k - 1
        j = x_mean_closest_k + 1
    else:
        i = x_mean_closest_k
        j = x_mean_closest_k + 2
    ic(i, j)
    y_mean = np.interp(x_mean, x.iloc[i:j], y.iloc[i:j])

y_mean

# %%
x.iloc[i:j]

# %%
y.iloc[i:j]

# %%
np.interp(x_mean, x, y)

# %%
# a version of the above plot,
# but the neural and non-neural data are comming from separtion of cells
# and not genes, as in the above plot

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
x_std = x.std()
ic(x_mean, x_std)
x_mean_closest = x.iloc[(x - x_mean).abs().argsort()[:1]]
x_mean_closest_k = x_mean_closest.index.values[0]
if x_mean == x_mean_closest.values[0]:
    y_mean = y.iloc[x_mean_closest_k]
else:
    # if x_mean < x_mean_closest.values[0]:
    #     i = x_mean_closest_k - 1
    #     j = x_mean_closest_k + 1
    # else:
    #     i = x_mean_closest_k
    #     j = x_mean_closest_k + 2
    # y_mean = np.interp(x_mean, x.iloc[i:j], y.iloc[i:j])
    y_mean = np.interp(x_mean, x, y)
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
tmr1000_x_std = tmr1000_x.std()
ic(tmr1000_x_mean, tmr1000_x_std)
tmr1000_x_mean_closest = tmr1000_x.iloc[
    (tmr1000_x - tmr1000_x_mean).abs().argsort()[:1]
]
tmr1000_x_mean_closest_k = tmr1000_x_mean_closest.index.values[0]
if tmr1000_x_mean == tmr1000_x_mean_closest.values[0]:
    tmr1000_y_mean = tmr1000_y.iloc[tmr1000_x_mean_closest_k]
else:
    # if tmr1000_x_mean < tmr1000_x_mean_closest.values[0]:
    #     i = tmr1000_x_mean_closest_k - 1
    #     j = tmr1000_x_mean_closest_k + 1
    # else:
    #     i = tmr1000_x_mean_closest_k
    #     j = tmr1000_x_mean_closest_k + 2
    # tmr1000_y_mean = np.interp(tmr1000_x_mean, tmr1000_x.iloc[i:j], tmr1000_y.iloc[i:j])
    tmr1000_y_mean = np.interp(tmr1000_x_mean, tmr1000_x, tmr1000_y)
tmr1000_df = tmr1000_df.drop_duplicates(subset="NumOfProteins").reset_index(drop=True)

neural_conditions = ["Yes", "No"]
neural_trace_names = ["Neural", "Non-neural"]
neural_color_discrete_map = {
    "Yes": "red",
    "No": "rgb(0,170,255)",  # kind of azure
}
neural_dfs = []
for neural_df in [
    neuronal_max_distinct_proteins_df,
    non_neuronal_max_distinct_proteins_df,
]:
    neural_df = (
        neural_df.loc[:, ["NumOfProteins"]]
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
        legend="legend",
        showlegend=False,
        # text=f"{x_mean:.0f} distinct proteins<br>(avg)",
        text=f"{x_mean:.0f} distinct<br>proteins<br>(avg)",
        # text=f"{x_mean:.0f} distinct<br>proteins<br>(avg, STD = {x_std:.0f})",
        # text=f"{x_mean:.0f} ± {x_std:.0f}<br>distinct proteins",
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
        legend="legend",
        text=f"{tmr1000_x_mean:.0f} distinct proteins<br>(avg)",
        # text=f"{tmr1000_x_mean:.0f} distinct proteins<br>(avg, STD = {tmr1000_x_std:.0f})",
        # text=f"{tmr1000_x_mean:.0f} ± {tmr1000_x_std:.0f}<br>distinct proteins",
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
    legend="legend",
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
    legend="legend",
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
        legend="legend",
    ),
    row=1,
    col=1,
)

# neural_conditions = ["Yes", "No"]
# neural_trace_names = ["Neural", "Non-neural"]

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
            # legendgroup=tmr50_legendtitle,  # this can be any string
            # legendgrouptitle_text=tmr50_legendtitle,
            legend="legend1",
        ),
        # row=1, col=2
        row=2,
        col=1,
    )

    y_min = min(y_min, y.min())


x = neuronal_max_distinct_proteins_df["NumOfProteins"]
y = non_neuronal_max_distinct_proteins_df["NumOfProteins"]
statistic, pv = scipy.stats.mannwhitneyu(x, y)
if pv < 10**-22:
    neural_vs_non_text = f"<b>Mann-Whitney U between<br>neural to non-neural cells</b><br>p-val < 1E-22<br>statistic = {statistic:.2g}"
else:
    neural_vs_non_text = f"<b>Mann-Whitney U between<br>neural to non-neural cells</b><br>p-val = {pv:.2e}<br>statistic = {statistic:.2g}"

fig.add_annotation(
    x=np.log(10) / np.log(10),
    # y=np.log(1) / np.log(10),
    y=np.log(10) / np.log(10),
    xref="x",
    yref="y",
    text=neural_vs_non_text,
    bgcolor="white",
    borderpad=4,
    font=dict(size=11),
    opacity=0.8,
    showarrow=False,
    row=2,
    col=1,
)

fig.update_xaxes(type="log", range=[0, None])
fig.update_yaxes(
    type="log",
    # range=[-2, 2.2]
    range=[np.log(y_min) * 1.1 / np.log(10), 2.2],
)

width = 800
# width = 650
height = 800

fig.update_layout(
    # xaxis_title="Distinct proteins per transcript",
    # yaxis_title="% of transcripts",
    # title="Pooled octopus data",
    title="Single-cell octopus data",
    title_x=0.15,
    template=template,
    width=width,
    height=height,
    # legend_title_text=legend_title_text,
    # legend_font=dict(size=10),
    # legend_grouptitlefont=dict(size=12),
    # showlegend=False,
    #     legend={
    #         "title": "Coverage per gene",
    #         # "xref": "container",
    #         # "yref": "container",
    #         # "xref": "paper",
    #         # "yref": "paper",
    #         "y": 0.9,
    #         # "bgcolor": "Orange",
    #     },
    #     legend1={
    #         "title": "Cells",
    #         # "xref": "container",
    #         # "yref": "container",
    #         # "xref": "paper",
    #         # "yref": "paper",
    #         # "y": 0.5,
    #         # "bgcolor": "Gold",
    #     },
)

fig.write_image(
    "Distinct proteins per gene vs. prct of genes - log(y) - Octopus - Single-cell.svg",
    width=width,
    height=height,
)

fig.show()

# %%
font_size = 24

neural_conditions = ["Yes", "No"]
neural_trace_names = ["Neural", "Non-neural"]
neural_color_discrete_map = {
    "Yes": "red",
    "No": "rgb(0,170,255)",  # kind of azure
}
neural_dfs = []
for neural_df in [
    neuronal_max_distinct_proteins_df,
    non_neuronal_max_distinct_proteins_df,
]:
    neural_df = (
        neural_df.loc[:, ["NumOfProteins"]]
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

# tmr50_legendtitle = "50 reads"
# tmr1000_legendtitle = "1000 reads"
# legend_title_text = "Minimum coverage per gene      "

marker_size = 4

fig = make_subplots(
    rows=1,
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


neural_conditions = ["Yes", "No"]
neural_trace_names = ["Neural", "Non-neural"]

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
            # legendgroup=tmr50_legendtitle,  # this can be any string
            # legendgrouptitle_text=tmr50_legendtitle,
            legend="legend1",
        ),
        # row=1, col=2
        row=1,
        col=1,
    )

    y_min = min(y_min, y.min())


x = neuronal_max_distinct_proteins_df["NumOfProteins"]
y = non_neuronal_max_distinct_proteins_df["NumOfProteins"]
statistic, pv = scipy.stats.mannwhitneyu(x, y)
if pv < 10**-22:
    # neural_vs_non_text = f"<b>Mann-Whitney U between<br>neural to non-neural cells</b><br>p-val < 1E-22<br>statistic = {statistic:.2g}"
    neural_vs_non_text = "<b>p-val < 1E-22</b><br>(Mann-Whitney U)"
else:
    # neural_vs_non_text = f"<b>Mann-Whitney U between<br>neural to non-neural cells</b><br>p-val = {pv:.2e}<br>statistic = {statistic:.2g}"
    neural_vs_non_text = f"p-val = {pv:.2e} (Mann-Whitney U)"

fig.update_annotations(font_size=font_size)

fig.add_annotation(
    x=np.log(10) / np.log(10),
    # y=np.log(1) / np.log(10),
    y=np.log(10) / np.log(10),
    xref="x",
    yref="y",
    text=neural_vs_non_text,
    bgcolor="white",
    borderpad=4,
    # font=dict(size=11),
    font=dict(size=0.7 * font_size),
    opacity=0.8,
    showarrow=False,
    row=1,
    col=1,
)

fig.update_xaxes(type="log", range=[0, None], tickfont=dict(size=0.7 * font_size))
fig.update_yaxes(
    type="log",
    # range=[-2, 2.2]
    range=[np.log(y_min) * 1.1 / np.log(10), 2.2],
    tickfont=dict(size=0.7 * font_size),
)

# fig.update_annotations(font_size=font_size)

# width = 800
width = 900
# width = 650
# height = 800
height = 600

fig.update_layout(
    # xaxis_title="Distinct proteins per transcript",
    # yaxis_title="% of transcripts",
    # title="Pooled octopus data",
    # legend_font=dict(size=0.7 * font_size),
    legend_font=dict(size=font_size),
    # title="Single-cell octopus data",
    title=dict(text="Single-cell octopus data", font=dict(size=font_size * 1.5)),
    # title_x=0.15,
    title_y=0.93,
    template=template,
    width=width,
    height=height,
    legend_title_text="Cells               ",
    # legend_font=dict(size=10),
    # legend_grouptitlefont=dict(size=12),
    # showlegend=False,
    #     legend={
    #         "title": "Coverage per gene",
    #         # "xref": "container",
    #         # "yref": "container",
    #         # "xref": "paper",
    #         # "yref": "paper",
    #         "y": 0.9,
    #         # "bgcolor": "Orange",
    #     },
    #     legend1={
    #         "title": "Cells",
    #         # "xref": "container",
    #         # "yref": "container",
    #         # "xref": "paper",
    #         # "yref": "paper",
    #         # "y": 0.5,
    #         # "bgcolor": "Gold",
    #     },
)

fig.write_image(
    "Distinct proteins per gene vs. prct of genes - log(y) - Octopus - Single-cell - bottom panel.svg",
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
x_std = x.std()
ic(x_std)
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
tmr1000_x_std = tmr1000_x.std()
ic(tmr1000_x_std)
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
        # text=f"{x_mean:.0f} distinct<br>proteins<br>(avg, STD = {x_std:.0f})",
        # text=f"{x_mean:.0f} ± {x_std:.0f}<br>distinct proteins",
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
        # text=f"{tmr1000_x_mean:.0f} distinct proteins<br>(avg, STD = {tmr1000_x_std:.0f})",
        # text=f"{tmr1000_x_mean:.0f} ± {tmr1000_x_std:.0f}<br>distinct proteins",
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
    title="Single-cell octopus data",
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

# fig.write_image(
#     "Distinct proteins per gene vs. % of genes - log(y) - Octopus.svg",
#     width=width,
#     height=height,
# )

fig.show()

# %%
# version with only top panel for a conferenct

font_size = 24

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
x_std = x.std()
ic(x_std)
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
tmr1000_x_std = tmr1000_x.std()
ic(tmr1000_x_std)
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
    # rows=2,
    rows=1,
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
        # text=f"{x_mean:.0f} distinct<br>proteins<br>(avg, STD = {x_std:.0f})",
        # text=f"{x_mean:.0f} ± {x_std:.0f}<br>distinct proteins",
        textposition="bottom left",
        textfont=dict(
            color=tmr50_all_color,
            # size=11
            size=0.7 * font_size,
        ),
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
        # text=f"{tmr1000_x_mean:.0f} distinct proteins<br>(avg, STD = {tmr1000_x_std:.0f})",
        # text=f"{tmr1000_x_mean:.0f} ± {tmr1000_x_std:.0f}<br>distinct proteins",
        textposition="top right",
        textfont=dict(
            color="green",
            # size=11
            size=0.7 * font_size,
        ),
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
        # textfont=dict(size=11),
        textfont=dict(size=0.7 * font_size),
        showlegend=False,
    ),
    row=1,
    col=1,
)

# for neural_condition, neural_trace_name, neural_df in zip(
#     neural_conditions, neural_trace_names, neural_dfs
# ):
#     if neural_condition == "Missing":
#         continue
#     x = neural_df["NumOfProteins"]
#     y = neural_df["CummulativeTranscripts"]
#     color = neural_color_discrete_map[neural_condition]

#     fig.add_trace(
#         go.Scatter(
#             x=x,
#             y=y,
#             mode="lines+markers",
#             marker=dict(color=color, size=marker_size),
#             line=dict(color=color, dash="dash", width=0.5),
#             name=neural_trace_name,
#             # legendgroup=tmr50_legendtitle,  # this can be any string
#             # legendgrouptitle_text=tmr50_legendtitle,
#             legend="legend1",
#         ),
#         # row=1, col=2
#         row=2,
#         col=1,
#     )

#     y_min = min(y_min, y.min())


# x = max_distinct_proteins_df.loc[
#     max_distinct_proteins_df["IsNeural"] == "Yes", "NumOfProteins"
# ]
# y = max_distinct_proteins_df.loc[
#     max_distinct_proteins_df["IsNeural"] == "No", "NumOfProteins"
# ]
# statistic, pv = scipy.stats.mannwhitneyu(x, y)

# fig.add_annotation(
#     x=np.log(10) / np.log(10),
#     y=np.log(1) / np.log(10),
#     xref="x",
#     yref="y",
#     text=f"<b>Mann-Whitney U between<br>neural to non-neural genes</b><br>p-val = {pv:.2e}<br>statistic = {statistic:.2g}",
#     bgcolor="white",
#     borderpad=4,
#     font=dict(size=11),
#     opacity=0.8,
#     showarrow=False,
#     row=2,
#     col=1,
# )

fig.update_xaxes(
    type="log",
    tickfont=dict(size=0.7 * font_size),
)
fig.update_yaxes(
    type="log",
    # range=[-2, 2.2]
    range=[np.log(y_min) * 1.1 / np.log(10), 2.2],
    tickfont=dict(size=0.7 * font_size),
)

fig.update_annotations(font_size=font_size)

width = 800
height = 600

fig.update_layout(
    # xaxis_title="Distinct proteins per transcript",
    # yaxis_title="% of transcripts",
    # title="Pooled octopus data",
    legend_font=dict(size=0.7 * font_size),
    # title="Whole-transcriptome octopus data",
    title=dict(text="Whole-transcriptome octopus data", font=dict(size=font_size)),
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

# fig.write_image(
#     "Distinct proteins per gene vs. % of genes - log(y) - Octopus.svg",
#     width=width,
#     height=height,
# )

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

# %%

# %%

# %%

# %% [markdown]
# ### Distinct isoforms per cell

# %%
# f1_5plus_exapnded_max_expression_dfs[0]

# %%
# len(f1_5plus_exapnded_max_expression_dfs)

# %%
f1_exapnded_max_expression_dfs[0]

# %%
len(f1_exapnded_max_expression_dfs)

# %%
copies_df = (
    # pd.concat(f1_5plus_exapnded_max_expression_dfs)
    pd.concat(f1_exapnded_max_expression_dfs)
    .groupby(
        [
            "Chrom",
            condition_col,
            "Sample",
            "CB",
            "SeuratCluster",
            "Annotation",
            "Protein",
        ],
        dropna=False,
    )
    .size()
    .reset_index(name="Copies")
    .merge(neuronality_of_annotaion_df, on="Annotation", how="left")
)

copies_df

# %%
copies_df

# %%
copies_df["Copies"].sum() <= raw_reads_info_df.shape[0]

# %%
raw_reads_info_df

# %%
reads_mapped_to_edited_positions_per_cell_df = (
    copies_df.groupby(["Sample", "Neuronal", "Annotation", "CB"], dropna=False)[
        "Copies"
    ]
    .sum()
    .reset_index(name="Reads")
)
reads_mapped_to_edited_positions_per_cell_df

# %%
reads_mapped_to_edited_positions_per_annotation_df = (
    copies_df.groupby(["Sample", "Neuronal", "Annotation"], dropna=False)["Copies"]
    .sum()
    .reset_index(name="Reads")
)
reads_mapped_to_edited_positions_per_annotation_df

# %%
reads_mapped_to_edited_positions_per_cell_df.groupby("Sample")["Reads"].describe()

# %%
reads_mapped_to_edited_positions_per_cell_df.groupby(["Sample", "Annotation"])[
    "Reads"
].describe()

# %%
p = (
    so.Plot(
        reads_per_cell_df,
        x="Reads",
        # color="Sample"
    ).add(so.Area(), so.Hist())
    # .add(so.Bars(), so.Hist())
    .scale(
        # x="log",
        y="log"
    )
    # .limit(x=(0, 100))
    .label(
        x="Reads per cell",
        y="Cells",
        # color="Editing detected in gene",
    )
)
p

# %%
copies_df.groupby([])

# %%
ge5_isoforms_per_gene_copies_df = copies_df.loc[
    copies_df["Chrom"].isin(chroms_with_at_least_5_isoforms)
]
ge5_isoforms_per_gene_copies_df

# %%
copies_df["Copies"].value_counts(dropna=False)

# %%
ge5_isoforms_per_gene_copies_df["Copies"].value_counts(dropna=False)

# %%
copies_df["Neuronal"].value_counts(dropna=False)

# %%
ge5_isoforms_per_gene_copies_df["Neuronal"].value_counts(dropna=False)

# %%
# annotations where it's unknown if the cell is neural
copies_df.loc[copies_df["Neuronal"].isna(), "Annotation"].value_counts(dropna=False)

# %%
# annotations where it's unknown if the cell is neural
ge5_isoforms_per_gene_copies_df.loc[
    ge5_isoforms_per_gene_copies_df["Neuronal"].isna(), "Annotation"
].value_counts(dropna=False)

# %%
# copies_df["Annotation"].value_counts(dropna=False)

# %%
# num of unique cells per sample
copies_df.groupby("Sample")["CB"].nunique()

# %%
# num of unique cells per sample
ge5_isoforms_per_gene_copies_df.groupby("Sample")["CB"].nunique()

# %%
# %time

ge5_per_sample_and_gene_diversity_of_isoforms_per_cell_dfs = {}

for sample in samples:

    per_sample_copies_df = ge5_isoforms_per_gene_copies_df.loc[
        ge5_isoforms_per_gene_copies_df["Sample"] == sample
    ].drop(columns=["SeuratCluster", "Annotation"])
    unique_chroms_per_sample = per_sample_copies_df["Chrom"].unique().tolist()

    for chrom in unique_chroms_per_sample:

        one_gene_diversity_of_isoforms_per_cell_df = (
            per_sample_copies_df.loc[per_sample_copies_df["Chrom"] == chrom]
            .pivot(
                index=["Sample", "Chrom", condition_col, "Protein"],
                columns="CB",
                values="Copies",
            )
            .rename_axis(None, axis=1)
        )

        ge5_per_sample_and_gene_diversity_of_isoforms_per_cell_dfs[(sample, chrom)] = (
            one_gene_diversity_of_isoforms_per_cell_df
        )

ic(len(ge5_per_sample_and_gene_diversity_of_isoforms_per_cell_dfs));

# %%
# df = list(ge5_per_sample_and_gene_diversity_of_isoforms_per_cell_dfs.values())[
#     0
# ].fillna(0)
# df

# %%
# def report_chi_2_assumptions(sample, chrom, df):
#     df = df.fillna(0)

#     row_sums = df.sum(axis=1).values.reshape(-1, 1)  # Sum of each row
#     col_sums = df.sum(axis=0).values.reshape(1, -1)  # Sum of each column
#     total = df.values.sum()  # Total sum of all elements

#     expected = (row_sums @ col_sums) / total  # Compute expected values matrix
#     expected_df = pd.DataFrame(expected, index=df.index, columns=df.columns)

#     num_of_all_cells = expected_df.shape[0] * expected_df.shape[1]
#     num_of_cells_with_expected_ge_5 = expected_df.map(lambda x: x >= 5).sum().sum()
#     num_of_cells_with_expected_lt_1 = expected_df.map(lambda x: x < 1).sum().sum()
#     prct_of_cells_with_expected_ge_5 = (
#         100 * num_of_cells_with_expected_ge_5 / num_of_all_cells
#     )
#     prct_of_cells_with_expected_lt_1 = (
#         100 * num_of_cells_with_expected_lt_1 / num_of_all_cells
#     )

#     return (
#         sample,
#         chrom,
#         num_of_all_cells,
#         prct_of_cells_with_expected_ge_5,
#         prct_of_cells_with_expected_lt_1,
#     )

# %%
# with Pool(processes=3) as pool:
#     per_sample_and_gene_chi_2_assumptions = pool.starmap(
#         func=report_chi_2_assumptions,
#         iterable=[
#             (sample, chrom, df)
#             for (
#                 sample,
#                 chrom,
#             ), df in ge5_per_sample_and_gene_diversity_of_isoforms_per_cell_dfs.items()
#             # list(lt5_per_sample_and_gene_diversity_of_isoforms_per_cell_dfs.items())[:30]
#         ],
#     )
# ge5_chi_2_assumptions = pd.DataFrame(
#     per_sample_and_gene_chi_2_assumptions,
#     columns=[
#         "Sample",
#         "Chrom",
#         "NumOfCells",
#         "%CellsWithExpected>=5",
#         "%CellsWithExpected<1",
#     ],
# )
# ge5_chi_2_assumptions

# %%
# ge5_chi_2_assumptions.loc[ge5_chi_2_assumptions["NumOfCells"] <= 20]

# %%
# ge5_chi_2_assumptions.loc[chi_2_assumptions["Chrom"] == "comp178306_c1_seq1"]

# %%
# ge5_chi_2_assumptions.loc[chi_2_assumptions["%CellsWithExpected>=5"] >= 20]

# %%
# ge5_chi_2_assumptions.groupby("Sample")["%CellsWithExpected<1"].describe()

# %% jupyter={"source_hidden": true}
# df = list(lt5_per_sample_and_gene_diversity_of_isoforms_per_cell_dfs.values())[0]
# res = scipy.stats.fisher_exact(df.fillna(0).values)
# statistic, pval = res.statistic, res.pvalue
# statistic, pval

# %% jupyter={"source_hidden": true}
# df = list(lt5_per_sample_and_gene_diversity_of_isoforms_per_cell_dfs.values())[0]
# res = scipy.stats.fisher_exact(df.T.fillna(0).values)
# statistic, pval = res.statistic, res.pvalue
# statistic, pval

# %% jupyter={"source_hidden": true}
# df = list(lt5_per_sample_and_gene_diversity_of_isoforms_per_cell_dfs.values())[0]
# res = scipy.stats.chi2_contingency(df.fillna(0).values)
# statistic, pval = res.statistic, res.pvalue
# statistic, pval

# %% jupyter={"source_hidden": true}
# df = list(lt5_per_sample_and_gene_diversity_of_isoforms_per_cell_dfs.values())[0]
# res = scipy.stats.chi2_contingency(df.T.fillna(0).values)
# statistic, pval = res.statistic, res.pvalue
# statistic, pval

# %% jupyter={"source_hidden": true}
# # %time


# # per_sample_and_gene_kruskal_wallis_isoforms_per_cell_results = {}
# per_sample_and_gene_kruskal_wallis_isoforms_per_cell_results = []

# for (
#     sample,
#     chrom,
# ), df in per_sample_and_gene_kruskal_wallis_isoforms_per_cell_dfs.items():
#     # ic(sample, chrom)
#     try:
#         # # statistic, pval = scipy.stats.kruskal(*df.fillna(0).values, nan_policy="omit")
#         # test_completed = True

#         statistic, pval = scipy.stats.kruskal(*df.values, nan_policy="omit")
#         test_completed = True

#         # try:
#         #     statistic, pval = scipy.stats.kruskal(
#         #         *df.dropna(axis=1).values, nan_policy="omit"
#         #     )
#         #     test_completed = True
#         # except RuntimeWarning as e:
#         #     # happens if df.dropna(axis=1).shape[1] < 5
#         #     if str(e) == "One or more sample arguments is too small; all returned values will be NaN. See documentation for sample size requirements.":
#         #         statistic = np.nan
#         #         pval = np.nan
#         #         test_completed = False
#         #     else:
#         #         raise

#     except ValueError as e:
#         if str(e) == "All numbers are identical in kruskal":
#             # Handle the ValueError if it has the specific message
#             statistic = np.nan
#             pval = 1
#             test_completed = False
#             # test_completed = True
#         else:
#             # Re-raise the error if the message does not match
#             raise
#     # per_sample_and_gene_kruskal_wallis_isoforms_per_cell_results[(sample, chrom)] = (test_completed, statistic, pval)
#     per_sample_and_gene_kruskal_wallis_isoforms_per_cell_results.append(
#         (sample, chrom, test_completed, statistic, pval)
#     )

# per_sample_and_gene_kruskal_wallis_isoforms_per_cell_results_df = pd.DataFrame(
#     per_sample_and_gene_kruskal_wallis_isoforms_per_cell_results,
#     columns=["Sample", "Chrom", "TestCompleted", "Statistic", "PVal"],
# )
# # ).sort_values(["Sample", "Chrom"], ignore_index=True)

# corrected_dfs = []
# for sample in samples:
#     df = per_sample_and_gene_kruskal_wallis_isoforms_per_cell_results_df.loc[
#         per_sample_and_gene_kruskal_wallis_isoforms_per_cell_results_df["Sample"]
#         == sample
#     ].copy()
#     rejected, corrected_pval = fdrcorrection(df["PVal"])
#     df["RejectedAfterCorrection"] = rejected
#     df["CorrectedPVal"] = corrected_pval
#     corrected_dfs.append(df)

# per_sample_and_gene_kruskal_wallis_isoforms_per_cell_results_df = pd.concat(
#     corrected_dfs, ignore_index=True
# )

# per_sample_and_gene_kruskal_wallis_isoforms_per_cell_results_df = (
#     per_sample_and_gene_kruskal_wallis_isoforms_per_cell_results_df.merge(
#         copies_df.loc[:, ["Sample", "Chrom", condition_col]].drop_duplicates(),
#         how="left",
#     ).loc[
#         :,
#         [
#             "Sample",
#             "Chrom",
#             condition_col,
#             "TestCompleted",
#             "Statistic",
#             "PVal",
#             "RejectedAfterCorrection",
#             "CorrectedPVal",
#         ],
#     ]
# )

# per_sample_and_gene_kruskal_wallis_isoforms_per_cell_results_df

# %%
# len(per_sample_and_gene_diversity_of_isoforms_per_cell_results)

# %%
def calc_exact_fisher_for_chrom_and_sample(
    sample, chrom, per_sample_and_gene_diversity_of_isoforms_per_cell_df
):
    try:
        res = scipy.stats.fisher_exact(
            per_sample_and_gene_diversity_of_isoforms_per_cell_df.fillna(0).values
        )
        statistic, pval = res.statistic, res.pvalue
        if (
            type(statistic) == type(pval) == np.ndarray
            and len(statistic) == len(pval) == 1
        ):
            statistic, pval = statistic[0], pval[0]
        test_completed = True
        return (sample, chrom, test_completed, statistic, pval)
    except:
        ic(sample, chrom)
        raise

# %%
# len(ge5_per_sample_and_gene_diversity_of_isoforms_per_cell_results)

# %% jupyter={"source_hidden": true}
# # %time


# with Pool(processes=3) as pool:
#     ge5_per_sample_and_gene_diversity_of_isoforms_per_cell_results = pool.starmap(
#         func=calc_exact_fisher_for_chrom_and_sample,
#         iterable=[
#             (sample, chrom, df)
#             for (
#                 sample,
#                 chrom,
#             ), df in lt5_per_sample_and_gene_diversity_of_isoforms_per_cell_dfs.items()
#             # list(lt5_per_sample_and_gene_diversity_of_isoforms_per_cell_dfs.items())[:30]
#         ],
#     )

# # # per_sample_and_gene_kruskal_wallis_isoforms_per_cell_results = {}
# # per_sample_and_gene_diversity_of_isoforms_per_cell_results = []

# # i = 0
# # for (
# #     sample,
# #     chrom,
# # ), df in lt5_per_sample_and_gene_diversity_of_isoforms_per_cell_dfs.items():

# #     # ic(sample, chrom)
# #     # res = scipy.stats.chi2_contingency(df.fillna(0).values)
# #     res = scipy.stats.fisher_exact(df.fillna(0).values)
# #     # statistic, pval = res.statistic, res.pvalue
# #     statistic, pval = res.statistic[0], res.pvalue[0]
# #     test_completed = True

# #     per_sample_and_gene_diversity_of_isoforms_per_cell_results.append(
# #         (sample, chrom, test_completed, statistic, pval)
# #     )

# #     i += 1
# #     if i == 10:
# #         break

# ge5_per_sample_and_gene_diversity_of_isoforms_per_cell_results_df = pd.DataFrame(
#     ge5_per_sample_and_gene_diversity_of_isoforms_per_cell_results,
#     columns=["Sample", "Chrom", "TestCompleted", "Statistic", "PVal"],
# )

# corrected_dfs = []
# for sample in samples:
#     df = per_sample_and_gene_diversity_of_isoforms_per_cell_results_df.loc[
#         per_sample_and_gene_diversity_of_isoforms_per_cell_results_df["Sample"]
#         == sample
#     ].copy()
#     rejected, corrected_pval = fdrcorrection(df["PVal"])
#     df["RejectedAfterCorrection"] = rejected
#     df["CorrectedPVal"] = corrected_pval
#     corrected_dfs.append(df)
# per_sample_and_gene_diversity_of_isoforms_per_cell_results_df = pd.concat(
#     corrected_dfs, ignore_index=True
# )

# ge5_per_sample_and_gene_diversity_of_isoforms_per_cell_results_df = (
#     per_sample_and_gene_diversity_of_isoforms_per_cell_results_df.merge(
#         copies_df.loc[:, ["Sample", "Chrom", condition_col]].drop_duplicates(),
#         how="left",
#     ).loc[
#         :,
#         [
#             "Sample",
#             "Chrom",
#             condition_col,
#             "TestCompleted",
#             "Statistic",
#             "PVal",
#             "RejectedAfterCorrection",
#             "CorrectedPVal",
#         ],
#     ]
# )

# ge5_per_sample_and_gene_diversity_of_isoforms_per_cell_results_df

# %%
# %time

ge5_per_sample_and_gene_diversity_of_isoforms_per_cell_results = []

i = 0
for (
    sample,
    chrom,
), df in ge5_per_sample_and_gene_diversity_of_isoforms_per_cell_dfs.items():

    # ic(sample, chrom)
    res = scipy.stats.chi2_contingency(df.fillna(0).values)
    statistic, pval = res.statistic, res.pvalue
    # res = scipy.stats.fisher_exact(df.fillna(0).values)
    # statistic, pval = res.statistic[0], res.pvalue[0]
    test_completed = True

    ge5_per_sample_and_gene_diversity_of_isoforms_per_cell_results.append(
        (sample, chrom, test_completed, statistic, pval)
    )

    i += 1
    # if i == 10:
    #     break

ge5_per_sample_and_gene_diversity_of_isoforms_per_cell_results_df = pd.DataFrame(
    ge5_per_sample_and_gene_diversity_of_isoforms_per_cell_results,
    columns=["Sample", "Chrom", "TestCompleted", "Statistic", "PVal"],
)

corrected_dfs = []
for sample in samples:
    df = ge5_per_sample_and_gene_diversity_of_isoforms_per_cell_results_df.loc[
        ge5_per_sample_and_gene_diversity_of_isoforms_per_cell_results_df["Sample"]
        == sample
    ].copy()
    rejected, corrected_pval = fdrcorrection(df["PVal"])
    df["RejectedAfterCorrection"] = rejected
    df["CorrectedPVal"] = corrected_pval
    corrected_dfs.append(df)
ge5_per_sample_and_gene_diversity_of_isoforms_per_cell_results_df = pd.concat(
    corrected_dfs, ignore_index=True
)

ge5_per_sample_and_gene_diversity_of_isoforms_per_cell_results_df = (
    ge5_per_sample_and_gene_diversity_of_isoforms_per_cell_results_df.merge(
        copies_df.loc[:, ["Sample", "Chrom", condition_col]].drop_duplicates(),
        how="left",
    ).loc[
        :,
        [
            "Sample",
            "Chrom",
            condition_col,
            "TestCompleted",
            "Statistic",
            "PVal",
            "RejectedAfterCorrection",
            "CorrectedPVal",
        ],
    ]
)

ge5_per_sample_and_gene_diversity_of_isoforms_per_cell_results_df

# %%
ge5_per_sample_and_gene_diversity_of_isoforms_per_cell_results_df.loc[
    ge5_per_sample_and_gene_diversity_of_isoforms_per_cell_results_df[
        "RejectedAfterCorrection"
    ]
]

# %%
# ge5_chi_2_assumptions.loc[ge5_chi_2_assumptions["NumOfCells"] <= 20]

# %%
ge5_per_sample_and_gene_diversity_of_isoforms_per_cell_results_df.loc[
    ge5_per_sample_and_gene_diversity_of_isoforms_per_cell_results_df["PVal"] < 0.1
]

# %% jupyter={"source_hidden": true}
# rejected_per_sample_and_gene_diversity_of_isoforms_per_cell_results_df = (
#     per_sample_and_gene_diversity_of_isoforms_per_cell_results_df.loc[
#         per_sample_and_gene_diversity_of_isoforms_per_cell_results_df[
#             "RejectedAfterCorrection"
#         ]
#     ]
# )

# assert rejected_per_sample_and_gene_diversity_of_isoforms_per_cell_results_df.loc[
#     ~rejected_per_sample_and_gene_diversity_of_isoforms_per_cell_results_df[
#         "TestCompleted"
#     ]
# ].empty

# rejected_per_sample_and_gene_diversity_of_isoforms_per_cell_results_df

# %% jupyter={"source_hidden": true}
# rejected_per_sample_and_gene_diversity_of_isoforms_per_cell_results_df.groupby(
#     "Sample"
# ).size()

# %%

# %%

# %%

# %%

# %%
ge5_well_annotated_copies_df = ge5_isoforms_per_gene_copies_df.loc[
    ~ge5_isoforms_per_gene_copies_df["Annotation"].isin(
        [
            "unstable",
            "not separated",
            "TBA1",
            "TBA2",
            "TBA3",
            "TBA4",
            "TBA5",
            "TBA6",
            "TBA7",
            "TBA8",
        ]
    )
].dropna(axis=0)
# ge5_well_annotated_copies_df = ge5_well_annotated_copies_df.merge(
#     neuronality_of_annotaion_df, on="Annotation", how="left"
# )
assert ge5_well_annotated_copies_df.loc[
    ge5_well_annotated_copies_df.isna().any(axis=1)
].empty
ge5_well_annotated_copies_df

# %%
ge5_well_annotated_copies_df.groupby("Sample")["Chrom"].nunique()

# %%
ge5_annotation_df = (
    ge5_well_annotated_copies_df.groupby(["Sample", "Chrom", condition_col, "Protein"])[
        "Annotation"
    ]
    .unique()
    .reset_index(name="UniqueAnnotations")
)
ge5_annotation_df["NumOfUniqueAnnotations"] = ge5_annotation_df[
    "UniqueAnnotations"
].apply(len)
ge5_annotation_df

# %%
ge5_annotation_df.groupby("Sample")["Chrom"].nunique()

# %%
ge5_annotation_df.groupby("Sample")["Chrom"].nunique().sum()


# %%
# annotation_df_2 = well_annotated_copies_df.groupby(["Sample", "Chrom", condition_col, "Protein", "Copies"])["Annotation"].value_counts().reset_index(name="Cells").rename(columns={"Copies": "CopiesPerCell"})
# annotation_df_2

# %%
def mean_jaccard_similarity(s):
    n = len(s)
    total_score = 0
    count = 0

    for i in range(n):
        set_i = set(s.iloc[i])
        # set_i = set(s[i])
        for j in range(i + 1, n):  # Only unique comparisons
            set_j = set(s.iloc[j])
            # set_j = set(s[j])
            intersection = len(set_i & set_j)
            union = len(set_i | set_j)
            score = intersection / union if union != 0 else 0
            total_score += score
            count += 1

    return total_score / count if count > 0 else 0


# Example usage:
s = pd.Series([["a", "b"], ["a", "b", "c"], ["c", "e"]])
mean_jaccard = mean_jaccard_similarity(s)
print(mean_jaccard)

# %%
ge5_jaccard_annotation_df = (
    ge5_annotation_df.groupby(["Sample", "Chrom", condition_col])["UniqueAnnotations"]
    .apply(mean_jaccard_similarity)
    .reset_index(name="AnnotationsJaccardIndex")
)
ge5_jaccard_annotation_df

# %%
ge5_jaccard_annotation_df.groupby("Sample")["AnnotationsJaccardIndex"].describe().round(
    2
)

# %%
fig = px.histogram(
    ge5_jaccard_annotation_df,
    x="AnnotationsJaccardIndex",
    color="Sample",
    # facet_col="Sample",
    # color_discrete_map=color_discrete_map,
    # labels={"AnnotationsJaccardIndex": "Unique genes per cell"},
    histnorm="percent",
    cumulative=True,
    # log_x=True,
    # log_y=True
)
fig.update_traces(opacity=0.75)
# fig.update_xaxes(title="Mean Jaccard Index over each<br>isoform's unique annotations", dtick=0.2)
fig.update_xaxes(
    title="Min avg Jaccard Index over each<br>gene's isoforms' unique annotations",
    dtick=0.2,
)
# fig.update_yaxes(title="Genes", dtick=25)
fig.update_yaxes(title="Genes [%]", dtick=25)
fig.update_layout(
    width=500,
    height=350,
    template=template,
    barmode="overlay",
    # title="Num of unique distinct isoforms per sample-gene-cell",
    # showlegend=False,
)
fig.show()

# %% jupyter={"source_hidden": true}
fig = px.histogram(
    ge5_jaccard_annotation_df,
    x="AnnotationsJaccardIndex",
    color="Sample",
    # facet_col="Sample",
    # color_discrete_map=color_discrete_map,
    # labels={"AnnotationsJaccardIndex": "Unique genes per cell"},
    # histnorm="percent",
    # cumulative=True,
    # log_x=True,
    # log_y=True
)
fig.update_traces(opacity=0.75)
fig.update_xaxes(
    title="Mean Jaccard Index over each<br>gene's isoforms' unique annotations",
    dtick=0.2,
)
# fig.update_xaxes(title="Min avg Jaccard Index over each<br>gene's isoforms' unique annotations", dtick=0.2)
fig.update_yaxes(title="Genes", dtick=25)
# fig.update_yaxes(title="Genes [%]", dtick=25)
fig.update_layout(
    width=500,
    height=350,
    template=template,
    barmode="overlay",
    # title="Num of unique distinct isoforms per sample-gene-cell",
    # showlegend=False,
)
fig.show()

# %% jupyter={"source_hidden": true}
fig = px.histogram(
    ge5_jaccard_annotation_df,
    x="AnnotationsJaccardIndex",
    color="Sample",
    # facet_col="Sample",
    # color_discrete_map=color_discrete_map,
    # labels={"AnnotationsJaccardIndex": "Unique genes per cell"},
    histnorm="percent",
    # cumulative=True,
    # log_x=True,
    # log_y=True
)
fig.update_traces(opacity=0.75)
fig.update_xaxes(
    title="Mean Jaccard Index over each<br>gene's isoforms' unique annotations",
    dtick=0.2,
)
# fig.update_xaxes(title="Min avg Jaccard Index over each<br>gene's isoforms' unique annotations", dtick=0.2)
# fig.update_yaxes(title="Genes", dtick=25)
fig.update_yaxes(title="Genes [%]", dtick=2)
fig.update_layout(
    width=500,
    height=350,
    template=template,
    barmode="overlay",
    # title="Num of unique distinct isoforms per sample-gene-cell",
    # showlegend=False,
)
fig.show()


# %%
def unique_elements_where_nan_equals_nan(elements):
    other_than_nan_elements = [x for x in elements if not pd.isna(x)]
    unique_elements_other_than_nan = set(other_than_nan_elements)
    # there's at least one nan value
    if len(other_than_nan_elements) < len(elements):
        unique_elements = unique_elements_other_than_nan | {np.nan}
    else:
        unique_elements = unique_elements_other_than_nan
    return sorted(unique_elements)


# %%
# ge5_isoforms_per_gene_copies_df

# %%
grouped_copies_df = copies_df.groupby(["Chrom", condition_col, "Sample", "CB"])

df1 = grouped_copies_df.agg(
    UniqueProteins=("Protein", "unique"),
    NumOfUniqueProteins=("Protein", "nunique"),
)
# df1

df2 = grouped_copies_df.agg(
    {
        "Protein": list,
        "SeuratCluster": lambda x: list(set(x))[0],
        "Annotation": lambda x: list(set(x))[0],
        "Copies": list,
        "Neuronal": lambda x: list(set(x))[0],
    }
).rename(
    columns={
        "Protein": "Proteins",
        "Copies": "CopiesPerProtein",
    }
)
# df2

per_sample_per_cb_copies_df = df1.join(df2).reset_index()

per_sample_per_cb_copies_df.insert(
    per_sample_per_cb_copies_df.columns.get_loc("CopiesPerProtein") + 1,
    "Copies",
    per_sample_per_cb_copies_df["CopiesPerProtein"].apply(sum),
)
# per_sample_per_cb_copies_df["Copies-UniqueIsoforms"] = (
#     per_sample_per_cb_copies_df["Copies"]
#     - per_sample_per_cb_copies_df["NumOfUniqueProteins"]
# )
per_sample_per_cb_copies_df["Copies/UniqueIsoforms"] = (
    per_sample_per_cb_copies_df["Copies"]
    / per_sample_per_cb_copies_df["NumOfUniqueProteins"]
)
per_sample_per_cb_copies_df["EachReadYieldsUniqueIsoform"] = (
    per_sample_per_cb_copies_df.apply(
        lambda x: x["NumOfUniqueProteins"] == x["Copies"], axis=1
    )
)

per_sample_per_cb_copies_df.insert(
    per_sample_per_cb_copies_df.columns.get_loc("Neuronal") + 1,
    "NeuronalStrRep",
    per_sample_per_cb_copies_df["Neuronal"].apply(
        lambda x: "NA" if pd.isna(x) else "Neuronal" if x else "Non-neuronal"
    ),
)

per_sample_per_cb_copies_df

# %%
ge5_per_sample_per_cb_copies_df = per_sample_per_cb_copies_df.loc[
    per_sample_per_cb_copies_df["Chrom"].isin(chroms_with_at_least_5_isoforms)
]
ge5_per_sample_per_cb_copies_df

# %%
per_sample_per_cb_copies_df["Neuronal"].value_counts(dropna=False)

# %%
ge5_per_sample_per_cb_copies_df["Neuronal"].value_counts(dropna=False)

# %% jupyter={"source_hidden": true}
# per_sample_per_cb_copies_df.groupby(
#     ["Sample", "Copies-UniqueIsoforms"]
# ).size().reset_index(name="GeneCellCouples").sort_values(
#     [
#         "Copies-UniqueIsoforms",
#         "Sample",
#     ],
#     ignore_index=True,
# )

# %%
# per_sample_per_cb_copies_df.groupby(
#     ["Sample", "Copies/UniqueIsoforms"]
# ).size().reset_index(name="GeneCellCouples").sort_values(
#     [
#         "Copies/UniqueIsoforms",
#         "Sample",
#     ],
#     ignore_index=True,
# )

# %%
# num of genes found for each cell in each sample
unique_genes_per_cell_df = (
    per_sample_per_cb_copies_df.groupby(["Sample", "CB"])["Chrom"]
    .nunique()
    .reset_index(name="NumOfUniqueGenes")
)
unique_genes_per_cell_df

# %%
# num of genes found for each cell in each sample
ge5_unique_genes_per_cell_df = (
    ge5_per_sample_per_cb_copies_df.groupby(["Sample", "CB"])["Chrom"]
    .nunique()
    .reset_index(name="NumOfUniqueGenes")
)
ge5_unique_genes_per_cell_df

# %%
all_and_ge5_genes_per_cell_df = pd.concat(
    [
        unique_genes_per_cell_df.assign(Genes="All genes"),
        ge5_unique_genes_per_cell_df.assign(Genes="Genes with 5+ isoforms"),
    ]
)
all_and_ge5_genes_per_cell_df

# %%
# num of cells found for each gene in each sample
unique_cells_per_gene_df = (
    per_sample_per_cb_copies_df.groupby(["Sample", "Chrom"])["CB"]
    .nunique()
    .reset_index(name="NumOfUniqueCells")
)

unique_cbs_per_sample = per_sample_per_cb_copies_df.groupby(["Sample"])["CB"].unique()
# unique_cbs_per_sample
unique_cells_per_gene_df["%UniqueCellsOfAllUniqueCells"] = (
    unique_cells_per_gene_df.apply(
        lambda x: 100
        * x["NumOfUniqueCells"]
        / unique_cbs_per_sample.loc[x["Sample"]].size,
        axis=1,
    )
)

unique_cells_per_gene_df

# %%
# num of cells found for each gene in each sample
ge5_unique_cells_per_gene_df = (
    ge5_per_sample_per_cb_copies_df.groupby(["Sample", "Chrom"])["CB"]
    .nunique()
    .reset_index(name="NumOfUniqueCells")
)

ge5_unique_cbs_per_sample = ge5_per_sample_per_cb_copies_df.groupby(["Sample"])[
    "CB"
].unique()
# unique_cbs_per_sample
ge5_unique_cells_per_gene_df["%UniqueCellsOfAllUniqueCells"] = (
    ge5_unique_cells_per_gene_df.apply(
        lambda x: 100
        * x["NumOfUniqueCells"]
        / ge5_unique_cbs_per_sample.loc[x["Sample"]].size,
        axis=1,
    )
)

ge5_unique_cells_per_gene_df

# %%
unique_genes_per_cell_df.groupby("Sample")["NumOfUniqueGenes"].describe().round(2)

# %%
ge5_unique_genes_per_cell_df.groupby("Sample")["NumOfUniqueGenes"].describe().round(2)

# %%
unique_cells_per_gene_df.groupby("Sample")["NumOfUniqueCells"].describe().round(1)

# %%
ge5_unique_cells_per_gene_df.groupby("Sample")["NumOfUniqueCells"].describe().round(1)

# %%
# unique cells-and-annotation per gene

# %%
unique_cells_per_gene_df.groupby("Sample")[
    "%UniqueCellsOfAllUniqueCells"
].describe().round(2)

# %%
ge5_unique_cells_per_gene_df.groupby("Sample")[
    "%UniqueCellsOfAllUniqueCells"
].describe().round(2)

# %%
# # width = 4.5
# width = 5
# height = 3

# %%
# p = (
#     so.Plot(unique_genes_per_cell_df, "NumOfUniqueGenes")
#     .layout(size=(width, height))
#     .add(so.Area(), so.Hist(), color="Sample")
#     .scale(x="log")
#     .label(x="Unique genes per cell", y="Cells", title="All genes")
# )
# p

# %%
# p = (
#     so.Plot(ge5_unique_genes_per_cell_df, "NumOfUniqueGenes")
#     .layout(size=(width, height))
#     .add(so.Area(), so.Hist(), color="Sample")
#     .scale(x="log")
#     .label(x="Unique genes per cell", y="Cells", title="Genes with 5+ isoforms")
# )
# p

# %%
# width = 4.5
width = 5
height = 3
p = (
    so.Plot(all_and_ge5_genes_per_cell_df, "NumOfUniqueGenes")
    .facet(row="Genes", wrap=3)
    .layout(size=(width, height * 2))
    .add(so.Area(), so.Hist(), color="Sample")
    .scale(
        x="log",
        y=so.Continuous().tick(every=20_000).label(unit=""),
    )
    .label(x="Unique genes per cell", y="Cells")
)
p

# %%

# %%
fig = px.histogram(
    unique_genes_per_cell_df,
    x="NumOfUniqueGenes",
    color="Sample",
    # facet_col="Sample",
    # color_discrete_map=color_discrete_map,
    labels={"NumOfUniqueGenes": "Unique genes per cell"},
    # histnorm="percent",
    # cumulative=True,
    log_x=True,
    # log_y=True
)
fig.update_traces(opacity=0.75)
fig.update_yaxes(title="Cells")
fig.update_layout(
    width=500,
    height=350,
    template=template,
    barmode="overlay",
    # title="Num of unique distinct isoforms per sample-gene-cell",
    # showlegend=False,
)
fig.show()

# %%
fig = px.histogram(
    unique_cells_per_gene_df,
    x="NumOfUniqueCells",
    color="Sample",
    # facet_col="Sample",
    # color_discrete_map=color_discrete_map,
    labels={"NumOfUniqueCells": "Unique cells per gene"},
    # histnorm="percent",
    # cumulative=True,
    # log_x=True,
    log_y=True,
)
fig.update_traces(opacity=0.75)
fig.update_yaxes(title="Genes")
fig.update_layout(
    width=500,
    height=350,
    template=template,
    barmode="overlay",
    # title="Num of unique distinct isoforms per sample-gene-cell",
    # showlegend=False,
)
fig.show()

# %%
fig = px.histogram(
    unique_cells_per_gene_df,
    x="%UniqueCellsOfAllUniqueCells",
    color="Sample",
    # facet_col="Sample",
    # color_discrete_map=color_discrete_map,
    labels={
        "%UniqueCellsOfAllUniqueCells": "Unique cells per gene /<br>all unique cells [%]"
    },
    histnorm="percent",
    # cumulative=True,
    # log_x=True,
    log_y=True,
)
fig.update_traces(opacity=0.75)
# fig.update_yaxes(title="Genes")
fig.update_yaxes(title="Genes [%]")
fig.update_layout(
    width=500,
    height=350,
    template=template,
    barmode="overlay",
    # title="Num of unique distinct isoforms per sample-gene-cell",
    # showlegend=False,
)
fig.show()

# %%
per_sample_per_cb_copies_df.groupby("Sample")["NumOfUniqueProteins"].describe()

# %%
# mean num of unique isoforms per each cell's genes
mean_unique_isoforms_per_cell_genes_df = (
    per_sample_per_cb_copies_df.groupby(["Sample", "CB"])["NumOfUniqueProteins"]
    .mean()
    .reset_index(name="MeanNumOfUniqueProteinsPerGene")
)
mean_unique_isoforms_per_cell_genes_df

# %%
mean_unique_isoforms_per_cell_genes_df.groupby("Sample")[
    "MeanNumOfUniqueProteinsPerGene"
].describe().round(2)

# %%
# mean num of unique isoforms per each gene's cells
mean_unique_isoforms_per_gene_cells_df = (
    per_sample_per_cb_copies_df.groupby(["Sample", "Chrom"])["NumOfUniqueProteins"]
    .mean()
    .reset_index(name="MeanNumOfUniqueProteinsPerCell")
)
mean_unique_isoforms_per_gene_cells_df

# %%
# mean_unique_isoforms_per_gene_cells_df["MeanNumOfUniqueProteinsPerCell"].describe()

# %%
mean_unique_isoforms_per_gene_cells_df.groupby("Sample")[
    "MeanNumOfUniqueProteinsPerCell"
].describe().round(2)

# %%
fig = px.histogram(
    mean_unique_isoforms_per_gene_cells_df,
    x="MeanNumOfUniqueProteinsPerCell",
    color="Sample",
    # facet_col="Sample",
    # color_discrete_map=color_discrete_map,
    labels={"MeanNumOfUniqueProteinsPerCell": "Unique isoforms per cell (avg)"},
    histnorm="percent",
    # cumulative=True
)
fig.update_traces(opacity=0.75)
fig.update_xaxes(tick0=0, dtick=10)
fig.update_yaxes(title="Genes [%]")
fig.update_layout(
    width=500,
    height=350,
    template=template,
    barmode="overlay",
    # title="Num of unique distinct isoforms per sample-gene-cell",
    # showlegend=False,
)
fig.show()

# %%
fig = px.histogram(
    mean_unique_isoforms_per_gene_cells_df,
    x="MeanNumOfUniqueProteinsPerCell",
    color="Sample",
    # facet_col="Sample",
    # color_discrete_map=color_discrete_map,
    labels={"MeanNumOfUniqueProteinsPerCell": "Min unique isoforms per cell (avg)"},
    histnorm="percent",
    cumulative=True,
)
fig.update_traces(opacity=0.75)
fig.update_traces(cumulative_direction="decreasing", selector=dict(type="histogram"))
fig.update_xaxes(tick0=0, dtick=10)
fig.update_yaxes(title="Genes [%]")
fig.update_layout(
    width=500,
    height=350,
    template=template,
    barmode="overlay",
    # title="Num of unique distinct isoforms per sample-gene-cell",
    # showlegend=False,
)
fig.show()

# %%
per_sample_per_cb_copies_df.groupby("Sample")["NumOfUniqueProteins"].describe().round(2)

# %%
# fig = px.histogram(
#     per_sample_per_cb_copies_df,
#     x="NumOfUniqueProteins",
#     color="Sample",
#     # facet_col="Sample",
#     # color_discrete_map=color_discrete_map,
#     labels={"NumOfUniqueProteins": "Unique isoforms per gene-cell couple"},
#     # histnorm="percent",
#     # cumulative=True
# )
# fig.update_traces(opacity=0.75)
# fig.update_xaxes(tick0=0, dtick=10)
# fig.update_yaxes(title="Gene-cell couples")
# fig.update_layout(
#     width=500,
#     height=350,
#     template=template,
#     barmode="overlay",
#     # title="Num of unique distinct isoforms per sample-gene-cell",
#     # showlegend=False,
# )
# fig.show()

# %%
# fig = px.histogram(
#     per_sample_per_cb_copies_df,
#     x="NumOfUniqueProteins",
#     color="Sample",
#     # facet_col="Sample",
#     # color_discrete_map=color_discrete_map,
#     labels={"NumOfUniqueProteins": "Min unique isoforms per gene-cell couple"},
#     histnorm="percent",
#     cumulative=True
# )
# fig.update_traces(opacity=0.75)
# fig.update_traces(cumulative_direction="decreasing", selector=dict(type='histogram'))
# fig.update_xaxes(tick0=0, dtick=10)
# fig.update_yaxes(title="Gene-cell couples [%]")
# fig.update_layout(
#     width=500,
#     height=350,
#     template=template,
#     barmode="overlay",
#     # title="Num of unique distinct isoforms per sample-gene-cell",
#     # showlegend=False,
# )
# fig.show()

# %%
fig = make_subplots(rows=1, cols=1)

for sample, color in zip(samples, colors):
    x = per_sample_per_cb_copies_df.loc[
        per_sample_per_cb_copies_df["Sample"] == sample, "Copies/UniqueIsoforms"
    ]
    fig.add_trace(
        go.Histogram(
            x=x,
            marker_color=color,
            histnorm="percent",
        ),
        row=1,
        col=1,
    )


fig.update_traces(opacity=0.75)
fig.update_xaxes(title="Reads per unique isoform (avg)")
fig.update_yaxes(title="Gene-cell couples [%]")
fig.update_layout(
    width=500,
    height=350,
    # template=template,
    barmode="overlay",
    # title="Num of unique distinct isoforms per sample-gene-cell",
    # showlegend=False,
)
fig.show()

# %%
fig = make_subplots(rows=1, cols=1)

for sample, color in zip(samples, colors):
    x = ge5_per_sample_per_cb_copies_df.loc[
        ge5_per_sample_per_cb_copies_df["Sample"] == sample, "Copies/UniqueIsoforms"
    ]
    fig.add_trace(
        go.Histogram(
            x=x,
            marker_color=color,
            histnorm="percent",
        ),
        row=1,
        col=1,
    )


fig.update_traces(opacity=0.75)
fig.update_xaxes(title="Reads per unique isoform (avg)")
fig.update_yaxes(title="Gene-cell couples [%]")
fig.update_layout(
    width=500,
    height=350,
    # template=template,
    barmode="overlay",
    # title="Num of unique distinct isoforms per sample-gene-cell",
    # showlegend=False,
)
fig.show()

# %%
# fig = px.histogram(
#     per_sample_per_cb_copies_df,
#     x="Copies/UniqueIsoforms",
#     color="Sample",
#     # facet_col="Sample",
#     # color_discrete_map=color_discrete_map,
#     # labels={"NumOfUniqueGenes": "Unique genes per cell"},
#     histnorm="percent",
#     # cumulative=True,
#     # log_x=True,
#     # log_y=True
# )
# fig.update_traces(opacity=0.75)
# fig.update_xaxes(title="Reads / unique isoforms")
# fig.update_yaxes(title="Gene-cell couples")
# fig.update_layout(
#     width=500,
#     height=350,
#     # template=template,
#     barmode="overlay",
#     # title="Num of unique distinct isoforms per sample-gene-cell",
#     # showlegend=False,
# )
# fig.show()

# %%
# fig = px.histogram(
#     per_sample_per_cb_copies_df,
#     x="Copies-UniqueIsoforms",
#     color="Sample",
#     # facet_col="Sample",
#     # color_discrete_map=color_discrete_map,
#     # labels={"NumOfUniqueGenes": "Unique genes per cell"},
#     # histnorm="percent",
#     # cumulative=True,
#     # log_x=True,
#     log_y=True,
# )
# fig.update_traces(opacity=0.75)
# fig.update_yaxes(title="Gene-cell Couples")
# fig.update_layout(
#     width=500,
#     height=350,
#     template=template,
#     barmode="overlay",
#     # title="Num of unique distinct isoforms per sample-gene-cell",
#     # showlegend=False,
# )
# fig.show()

# %% editable=true slideshow={"slide_type": ""}
per_sample_per_copies_stats_df = (
    per_sample_per_cb_copies_df.groupby(["Sample", "Copies"])
    .apply(
        lambda x: 100 * x["EachReadYieldsUniqueIsoform"].sum() / x.size,
        include_groups=False,
    )
    # .round(2)
    .reset_index(name="RelPrctOfGeneCellCouplesWithXReadsAndXUniqueIsoforms")
    .merge(
        per_sample_per_cb_copies_df.groupby(["Sample", "Copies"])
        .apply(lambda x: x["NumOfUniqueProteins"].mean(), include_groups=False)
        .reset_index(name="MeanNumOfUniqueProteinsPerGeneCellCoupleWithXReads")
    )
)
per_sample_per_copies_stats_df[
    "MeanNumOfUniqueProteinsPerGeneCellCoupleWithXReads/Reads"
] = (
    per_sample_per_copies_stats_df["MeanNumOfUniqueProteinsPerGeneCellCoupleWithXReads"]
    / per_sample_per_copies_stats_df["Copies"]
)
per_sample_per_copies_stats_df[
    "Reads/MeanNumOfUniqueProteinsPerGeneCellCoupleWithXReads"
] = (
    per_sample_per_copies_stats_df["Copies"]
    / per_sample_per_copies_stats_df[
        "MeanNumOfUniqueProteinsPerGeneCellCoupleWithXReads"
    ]
)
per_sample_per_copies_stats_df

# %% jupyter={"source_hidden": true}
# fig = px.scatter(
#     per_sample_per_copies_stats_df,
#     x="Copies",
#     y="RelPrctOfGeneCellCouplesWithXReadsAndXUniqueIsoforms",
#     color="Sample",
#     # log_x=True,
#     # log_y=True
#     # marginal_x="box",
#     # marginal_y="box",
# )
# fig.update_traces(opacity=0.75)
# fig.update_xaxes(title="Reads")
# # fig.update_yaxes(title="% of gene-cell couples with X reads<br>whose num of unique isoforms<br>is exactly X")
# # fig.update_yaxes(title="Relative % of gene-cell couples with<br>X reads and exactly X unique isoforms")
# fig.update_yaxes(
#     title="Relative % of gene-cell couples<br>with X reads and exactly<br>X unique isoforms"
# )
# fig.update_layout(
#     width=500,
#     height=350,
#     template=template,
#     # barmode="overlay",
#     # title="Num of unique distinct isoforms per sample-gene-cell",
#     # showlegend=False,
# )
# fig.show()

# %% jupyter={"source_hidden": true}
# fig = px.scatter(
#     per_sample_per_copies_stats_df,
#     x="Copies",
#     y="MeanNumOfUniqueProteinsPerGeneCellCoupleWithXReads",
#     color="Sample",
#     # log_x=True,
#     # log_y=True
#     # marginal_x="box",
#     # marginal_y="box",
# )
# fig.update_traces(opacity=0.75)
# fig.update_xaxes(title="Reads")
# fig.update_yaxes(title="Unique isoforms per gene-cell<br>couple with X reads (avg)")
# fig.update_layout(
#     width=500,
#     height=350,
#     template=template,
#     # barmode="overlay",
#     # title="Num of unique distinct isoforms per sample-gene-cell",
#     # showlegend=False,
# )
# fig.show()

# %% jupyter={"source_hidden": true}
fig = px.scatter(
    per_sample_per_copies_stats_df,
    x="Copies",
    y="RelPrctOfGeneCellCouplesWithXReadsAndXUniqueIsoforms",
    color="Sample",
    log_x=True,
    # log_y=True
    # marginal_x="box",
    # marginal_y="box",
)
fig.update_traces(opacity=0.75)
fig.update_xaxes(title="Reads")
# fig.update_yaxes(title="% of gene-cell couples with X reads<br>whose num of unique isoforms<br>is exactly X")
# fig.update_yaxes(title="Relative % of gene-cell couples with<br>X reads and exactly X unique isoforms")
fig.update_yaxes(
    title="Relative % of gene-cell couples<br>with X reads and exactly<br>X unique isoforms"
)
fig.update_layout(
    width=500,
    height=350,
    template=template,
    # barmode="overlay",
    # title="Num of unique distinct isoforms per sample-gene-cell",
    # showlegend=False,
)
fig.show()

# %%
fig = px.scatter(
    per_sample_per_copies_stats_df,
    x="Copies",
    y="MeanNumOfUniqueProteinsPerGeneCellCoupleWithXReads",
    color="Sample",
    log_x=True,
    # log_y=True
    # marginal_x="box",
    # marginal_y="box",
)
fig.update_traces(opacity=0.75)
fig.update_xaxes(title="Reads")
fig.update_yaxes(title="Unique isoforms per gene-cell<br>couple with X reads (avg)")
fig.update_layout(
    width=500,
    height=350,
    template=template,
    # barmode="overlay",
    # title="Num of unique distinct isoforms per sample-gene-cell",
    # showlegend=False,
)
fig.show()

# %%
fig = px.scatter(
    per_sample_per_copies_stats_df,
    x="Copies",
    y="MeanNumOfUniqueProteinsPerGeneCellCoupleWithXReads",
    facet_col="Sample",
    color="Neuronal",
    log_x=True,
    # log_y=True
    # marginal_x="box",
    # marginal_y="box",
)
fig.update_traces(opacity=0.75)
fig.update_xaxes(title="Reads")
fig.update_yaxes(title="Unique isoforms per gene-cell<br>couple with X reads (avg)")
fig.update_layout(
    width=500,
    height=350,
    template=template,
    # barmode="overlay",
    # title="Num of unique distinct isoforms per sample-gene-cell",
    # showlegend=False,
)
fig.show()

# %% jupyter={"source_hidden": true}
fig = px.scatter(
    per_sample_per_copies_stats_df,
    x="Copies",
    y="MeanNumOfUniqueProteinsPerGeneCellCoupleWithXReads/Reads",
    color="Sample",
    log_x=True,
    # log_y=True
    # marginal_x="box",
    # marginal_y="box",
    # labels={"MeanNumOfUniqueProteinsPerGeneCellCoupleWithXReads/Reads": ""}
)
fig.update_traces(opacity=0.75)
fig.update_xaxes(title="Reads")
fig.update_yaxes(
    title="Unique isoforms per gene-cell<br>couple with X reads (avg)<br>/ reads"
)
fig.update_layout(
    width=500,
    height=350,
    template=template,
    # barmode="overlay",
    # title="Num of unique distinct isoforms per sample-gene-cell",
    # showlegend=False,
)
fig.show()

# %%
# fig = px.scatter(
#     per_sample_per_copies_stats_df,
#     x="Copies",
#     y="Reads/MeanNumOfUniqueProteinsPerGeneCellCoupleWithXReads",
#     color="Sample",
#     # log_x=True,
#     # log_y=True
#     # marginal_x="box",
#     # marginal_y="box",
#     # labels={"MeanNumOfUniqueProteinsPerGeneCellCoupleWithXReads/Reads": ""}
# )
# fig.update_traces(opacity=0.75)
# fig.update_xaxes(title="Reads")
# fig.update_yaxes(
#     title="Reads /<br>unique isoforms per gene-cell<br>couple with X reads (avg)<br>"
# )
# fig.update_layout(
#     width=500,
#     height=350,
#     template=template,
#     # barmode="overlay",
#     # title="Num of unique distinct isoforms per sample-gene-cell",
#     # showlegend=False,
# )
# fig.show()

# %%
rows = 3

fig = make_subplots(
    rows=rows,
    # rows=3,
    cols=1,
    # subplot_titles=("Plot 1", "Plot 2", "Plot 3", "Plot 4")
    # shared_xaxes="all",
    x_title="Reads",
    vertical_spacing=0.05,
)

y_cols = [
    "RelPrctOfGeneCellCouplesWithXReadsAndXUniqueIsoforms",
    "MeanNumOfUniqueProteinsPerGeneCellCoupleWithXReads",
    "MeanNumOfUniqueProteinsPerGeneCellCoupleWithXReads/Reads",
]
y_titles = [
    "Relative % of gene-cell couples<br>with X reads and exactly<br>X unique isoforms",
    "Unique isoforms per gene-cell<br>couple with X reads (avg)",
    "Unique isoforms per gene-cell<br>couple with X reads (avg)<br>/ reads",
]
colors = ["blue", "red"]

mode = "markers"
opacity = 0.75

for row, (y_col, y_title) in enumerate(zip(y_cols, y_titles), start=1):
    # add trace
    for sample, color in zip(samples, colors):
        x = per_sample_per_copies_stats_df.loc[
            per_sample_per_copies_stats_df["Sample"] == sample, "Copies"
        ]
        y = per_sample_per_copies_stats_df.loc[
            per_sample_per_copies_stats_df["Sample"] == sample, y_col
        ]
        if row == 1:
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=y,
                    mode=mode,
                    marker_color=color,
                    opacity=opacity,
                    name=sample,
                    legendgroup=sample,
                ),
                row=row,
                col=1,
            )
        else:
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=y,
                    mode=mode,
                    marker_color=color,
                    opacity=opacity,
                    showlegend=False,
                ),
                row=row,
                col=1,
            )
    # Update yaxis properties
    fig.update_yaxes(title_text=y_title, row=row, col=1)

# Update xaxis properties
fig.update_xaxes(
    type="log",
)

# Update title and height
fig.update_layout(
    template=template,
    # title_text="Customizing Subplot Axes",
    width=500,
    height=300 * rows,
)

fig.show()

# %% jupyter={"source_hidden": true}
# base_rows = 3
# # base_rows = 1

# fig = make_subplots(
#     rows=base_rows*2,
#     # rows=3,
#     cols=1,
#     row_heights=[0.2, 0.8] * base_rows,
#     # subplot_titles=("Plot 1", "Plot 2", "Plot 3", "Plot 4")
#     # shared_xaxes="all",
#     x_title="Reads",
#     vertical_spacing=0.05,
# )

# y_cols = [
#     "RelPrctOfGeneCellCouplesWithXReadsAndXUniqueIsoforms",
#     "MeanNumOfUniqueProteinsPerGeneCellCoupleWithXReads",
#     "MeanNumOfUniqueProteinsPerGeneCellCoupleWithXReads/Reads",
# ]
# y_titles = [
#     "Relative % of gene-cell couples<br>with X reads and exactly<br>X unique isoforms",
#     "Unique isoforms per gene-cell<br>couple with X reads (avg)",
#     "Unique isoforms per gene-cell<br>couple with X reads (avg)<br>/ reads",
# ]
# colors = ["blue", "red"]

# mode = "markers"

# # main_cols = list(range(1, base_rows * 2, 2))
# main_cols = list(range(2, base_rows * 2 + 1, 2))

# for y_col, y_title, main_row in zip(y_cols, y_titles, main_cols):
#     ic(y_col, y_title, main_row)
#     # add trace
#     for sample, color in zip(samples, colors):
#         x = per_sample_per_copies_stats_df.loc[
#             per_sample_per_copies_stats_df["Sample"] == sample, "Copies"
#         ]
#         y = per_sample_per_copies_stats_df.loc[
#             per_sample_per_copies_stats_df["Sample"] == sample, y_col
#         ]
#         if main_row == 2:
#             fig.add_trace(
#                 go.Scatter(
#                     x=x,
#                     y=y,
#                     mode=mode,
#                     marker_color=color,
#                     name=sample,
#                     legendgroup=sample,
#                 ),
#                 row=main_row,
#                 col=1,
#             )
#         else:
#             fig.add_trace(
#                 go.Scatter(x=x, y=y, mode=mode, marker_color=color, showlegend=False),
#                 row=main_row,
#                 col=1,
#             )
#         fig.add_trace(
#                 go.Histogram(x=x, marker_color=color, showlegend=False, opacity=0.75),
#                 row=main_row-1,
#                 col=1,
#             )
#     # Update yaxis properties
#     fig.update_yaxes(title_text=y_title, row=main_row, col=1)
#     # break


# # Update xaxis properties
# fig.update_xaxes(
#     # title_text="Reads",
#     tick0=0,
#     dtick=1000,
#     # range=[0, x.max()],
#     # type="log",
# )

# # Update title and height
# fig.update_layout(
#     template=template,
#     # title_text="Customizing Subplot Axes",
#     width=500,
#     height=350*base_rows,
#     barmode="overlay"
# )

# fig.show()

# %%
# per_sample_per_cb_copies_df.loc[per_sample_per_cb_copies_df["NumOfUniqueProteins"] == 1].shape[0]

# %%
# per_sample_per_cb_copies_df.loc[per_sample_per_cb_copies_df["NumOfUniqueProteins"] == 2].shape[0]

# %%
# per_sample_per_cb_copies_df.loc[per_sample_per_cb_copies_df["NumOfUniqueProteins"] == 3].shape[0]

# %%
# per_sample_per_cb_copies_df.loc[per_sample_per_cb_copies_df["NumOfUniqueProteins"] == 4].shape[0]

# %%
# per_sample_per_cb_copies_df.loc[per_sample_per_cb_copies_df["NumOfUniqueProteins"] == 5].shape[0]

# %%
per_sample_per_cb_copies_df

# %%
agg_gene_cell_couples_df_1 = (
    per_sample_per_cb_copies_df.groupby(["Sample", "NumOfUniqueProteins"])
    .agg(
        NumOfGeneCellCouples=("Chrom", "size"),
        NumOfUniqueGenes=("Chrom", "nunique"),
        NumOfUniqueCells=("CB", "nunique"),
    )
    .reset_index()
    .rename(columns={"NumOfUniqueProteins": "NumOfUniqueProteinsPerGeneCell"})
)
# agg_gene_cell_couples_df_1

agg_gene_cell_couples_df_2 = (
    per_sample_per_cb_copies_df.groupby(["Sample", "NumOfUniqueProteins", "Chrom"])[
        "UniqueProteins"
    ]
    .apply(lambda x: len(set(chain.from_iterable(x))))
    .reset_index(name="NumOfUniqueProteinsPerChrom")
    .groupby(["Sample", "NumOfUniqueProteins"])["NumOfUniqueProteinsPerChrom"]
    .apply(list)
    .reset_index(name="NumOfUniqueProteinsPerChroms")
    .rename(columns={"NumOfUniqueProteins": "NumOfUniqueProteinsPerGeneCell"})
)
agg_gene_cell_couples_df_2["AvgNumOfUniqueProteinsPerGene"] = (
    agg_gene_cell_couples_df_2["NumOfUniqueProteinsPerChroms"].apply(
        lambda x: sum(x) / len(x)
    )
)
agg_gene_cell_couples_df_2["NumOfUniqueProteinsInAllGenes"] = (
    agg_gene_cell_couples_df_2["NumOfUniqueProteinsPerChroms"].apply(sum)
)
del agg_gene_cell_couples_df_2["NumOfUniqueProteinsPerChroms"]
# agg_gene_cell_couples_df_2

agg_gene_cell_couples_df = agg_gene_cell_couples_df_1.merge(agg_gene_cell_couples_df_2)
del agg_gene_cell_couples_df_1, agg_gene_cell_couples_df_2

dfs = []
for sample in samples:
    df = agg_gene_cell_couples_df.loc[
        agg_gene_cell_couples_df["Sample"] == sample
    ].sort_values("NumOfUniqueProteinsPerGeneCell", ascending=False)
    df["DecreasingCummlativeGeneCellCouples%"] = (
        df["NumOfGeneCellCouples"]
        .mul(100)
        .div(df["NumOfGeneCellCouples"].sum())
        .cumsum()
    )
    dfs.append(df)
agg_gene_cell_couples_df = pd.concat(dfs).sort_values(
    ["Sample", "NumOfUniqueProteinsPerGeneCell"], ignore_index=True
)

agg_gene_cell_couples_df

# %%
# ge5_per_sample_per_cb_copies_df

# %%
ge5_agg_gene_cell_couples_df_1 = (
    ge5_per_sample_per_cb_copies_df.groupby(["Sample", "NumOfUniqueProteins"])
    .agg(
        NumOfGeneCellCouples=("Chrom", "size"),
        NumOfUniqueGenes=("Chrom", "nunique"),
        NumOfUniqueCells=("CB", "nunique"),
    )
    .reset_index()
    .rename(columns={"NumOfUniqueProteins": "NumOfUniqueProteinsPerGeneCell"})
)
# agg_gene_cell_couples_df_1

ge5_agg_gene_cell_couples_df_2 = (
    ge5_per_sample_per_cb_copies_df.groupby(["Sample", "NumOfUniqueProteins", "Chrom"])[
        "UniqueProteins"
    ]
    .apply(lambda x: len(set(chain.from_iterable(x))))
    .reset_index(name="NumOfUniqueProteinsPerChrom")
    .groupby(["Sample", "NumOfUniqueProteins"])["NumOfUniqueProteinsPerChrom"]
    .apply(list)
    .reset_index(name="NumOfUniqueProteinsPerChroms")
    .rename(columns={"NumOfUniqueProteins": "NumOfUniqueProteinsPerGeneCell"})
)
ge5_agg_gene_cell_couples_df_2["AvgNumOfUniqueProteinsPerGene"] = (
    ge5_agg_gene_cell_couples_df_2["NumOfUniqueProteinsPerChroms"].apply(
        lambda x: sum(x) / len(x)
    )
)
ge5_agg_gene_cell_couples_df_2["NumOfUniqueProteinsInAllGenes"] = (
    ge5_agg_gene_cell_couples_df_2["NumOfUniqueProteinsPerChroms"].apply(sum)
)
del ge5_agg_gene_cell_couples_df_2["NumOfUniqueProteinsPerChroms"]
# agg_gene_cell_couples_df_2

ge5_agg_gene_cell_couples_df = ge5_agg_gene_cell_couples_df_1.merge(
    ge5_agg_gene_cell_couples_df_2
)
del ge5_agg_gene_cell_couples_df_1, ge5_agg_gene_cell_couples_df_2

dfs = []
for sample in samples:
    df = ge5_agg_gene_cell_couples_df.loc[
        ge5_agg_gene_cell_couples_df["Sample"] == sample
    ].sort_values("NumOfUniqueProteinsPerGeneCell", ascending=False)
    df["DecreasingCummlativeGeneCellCouples%"] = (
        df["NumOfGeneCellCouples"]
        .mul(100)
        .div(df["NumOfGeneCellCouples"].sum())
        .cumsum()
    )
    dfs.append(df)
ge5_agg_gene_cell_couples_df = pd.concat(dfs).sort_values(
    ["Sample", "NumOfUniqueProteinsPerGeneCell"], ignore_index=True
)

ge5_agg_gene_cell_couples_df

# %%
ge5_neuronal_agg_gene_cell_couples_df_1 = (
    ge5_per_sample_per_cb_copies_df.groupby(
        ["Sample", "NumOfUniqueProteins", "Neuronal"], dropna=False
    )
    .agg(
        NumOfGeneCellCouples=("Chrom", "size"),
        NumOfUniqueGenes=("Chrom", "nunique"),
        NumOfUniqueCells=("CB", "nunique"),
    )
    .reset_index()
    .rename(columns={"NumOfUniqueProteins": "NumOfUniqueProteinsPerGeneCell"})
)
# ge5_neuronal_agg_gene_cell_couples_df_1

ge5_neuronal_agg_gene_cell_couples_df_2 = (
    ge5_per_sample_per_cb_copies_df.groupby(
        ["Sample", "NumOfUniqueProteins", "Neuronal", "Chrom"], dropna=False
    )["UniqueProteins"]
    .apply(lambda x: len(set(chain.from_iterable(x))))
    .reset_index(name="NumOfUniqueProteinsPerChrom")
    .groupby(["Sample", "NumOfUniqueProteins", "Neuronal"], dropna=False)[
        "NumOfUniqueProteinsPerChrom"
    ]
    .apply(list)
    .reset_index(name="NumOfUniqueProteinsPerChroms")
    .rename(columns={"NumOfUniqueProteins": "NumOfUniqueProteinsPerGeneCell"})
)
ge5_neuronal_agg_gene_cell_couples_df_2["AvgNumOfUniqueProteinsPerGene"] = (
    ge5_neuronal_agg_gene_cell_couples_df_2["NumOfUniqueProteinsPerChroms"].apply(
        lambda x: sum(x) / len(x)
    )
)
ge5_neuronal_agg_gene_cell_couples_df_2["NumOfUniqueProteinsInAllGenes"] = (
    ge5_neuronal_agg_gene_cell_couples_df_2["NumOfUniqueProteinsPerChroms"].apply(sum)
)
del ge5_neuronal_agg_gene_cell_couples_df_2["NumOfUniqueProteinsPerChroms"]
# ge5_neuronal_agg_gene_cell_couples_df_2

ge5_neuronal_agg_gene_cell_couples_df = ge5_neuronal_agg_gene_cell_couples_df_1.merge(
    ge5_neuronal_agg_gene_cell_couples_df_2
)
del ge5_neuronal_agg_gene_cell_couples_df_1, ge5_neuronal_agg_gene_cell_couples_df_2
# ge5_neuronal_agg_gene_cell_couples_df

samples = ["nuclearRNA", "totalRNA"]
dfs = []
for sample in samples:
    df = ge5_neuronal_agg_gene_cell_couples_df.loc[
        ge5_neuronal_agg_gene_cell_couples_df["Sample"] == sample
    ].sort_values("NumOfUniqueProteinsPerGeneCell", ascending=False)
    # per neuronal/non-neuronal/NA separately
    unique_neuronalities_per_sample = df["Neuronal"].unique()
    for neuronality in unique_neuronalities_per_sample:
        if pd.isna(neuronality):
            df_of_one_neuronality = df.loc[df["Neuronal"].isna()].copy()
        else:
            df_of_one_neuronality = df.loc[df["Neuronal"].eq(neuronality)].copy()
        df_of_one_neuronality["DecreasingCummlativeGeneCellCouples%"] = (
            df_of_one_neuronality["NumOfGeneCellCouples"]
            .mul(100)
            .div(df_of_one_neuronality["NumOfGeneCellCouples"].sum())
            .cumsum()
        )
        dfs.append(df_of_one_neuronality)
ge5_neuronal_agg_gene_cell_couples_df = pd.concat(dfs).sort_values(
    ["Sample", "NumOfUniqueProteinsPerGeneCell"], ignore_index=True
)

ge5_neuronal_agg_gene_cell_couples_df["NeuronalStrRep"] = (
    ge5_neuronal_agg_gene_cell_couples_df["Neuronal"].apply(
        lambda x: "NA" if pd.isna(x) else "Neuronal" if x else "Non-neuronal"
    )
)

ge5_neuronal_agg_gene_cell_couples_df

# %%
agg_gene_cell_couples_df.groupby("Sample")[
    [
        "NumOfGeneCellCouples",
        "NumOfUniqueGenes",
        "NumOfUniqueCells",
        "AvgNumOfUniqueProteinsPerGene",
    ]
].agg(["count", "mean", "std", "min", "median", "max"]).round(2)

# %%
ge5_agg_gene_cell_couples_df.groupby("Sample")[
    [
        "NumOfGeneCellCouples",
        "NumOfUniqueGenes",
        "NumOfUniqueCells",
        "AvgNumOfUniqueProteinsPerGene",
    ]
].agg(["count", "mean", "std", "min", "median", "max"]).round(2)

# %%
ge5_neuronal_agg_gene_cell_couples_df.groupby(["Sample", "Neuronal"], dropna=False)[
    [
        "NumOfGeneCellCouples",
        "NumOfUniqueGenes",
        "NumOfUniqueCells",
        "AvgNumOfUniqueProteinsPerGene",
    ]
].agg(["count", "mean", "std", "min", "median", "max"]).round(2)

# %%
fig = px.line(
    agg_gene_cell_couples_df,
    x="NumOfUniqueProteinsPerGeneCell",
    y="NumOfGeneCellCouples",
    color="Sample",
    markers=True,
    # facet_col="Sample",
    # color_discrete_map=color_discrete_map,
    # labels={"NumOfUniqueProteins": "Min unique isoforms per gene-cell couple"},
    # histnorm="percent",
    # cumulative=True
    # log_x=True,
    log_y=True,
)
fig.update_traces(opacity=0.75, marker_size=4)
# fig.update_traces(cumulative_direction="decreasing", selector=dict(type='histogram'))
fig.update_xaxes(tick0=0, dtick=10, title="Unique isoforms per gene-cell couple")
fig.update_yaxes(title="Gene-cell couples", range=[0, 6])
fig.update_layout(
    width=500,
    height=350,
    template=template,
    # barmode="overlay",
    # title="Num of unique distinct isoforms per sample-gene-cell",
    # showlegend=False,
)
fig.show()

# %%
# width = 4.5
# height = 4.5
# p = (
#     so.Plot(ge5_agg_gene_cell_couples_df, x="NumOfUniqueProteinsPerGeneCell", y="NumOfGeneCellCouples")
#     .layout(size=(width, height))
#     .add(so.Line(), color="Sample")
#     .scale(
#         y="log",
#         x=so.Continuous().tick(every=10)
#     )
#     .label(x="Unique isoforms per gene-cell couple", y="Gene-cell couples")
#     .limit(x=(0, None), y=(0, None))
# )
# p

# %%
width = 4.5
height = 4.5
p = (
    so.Plot(
        ge5_agg_gene_cell_couples_df,
        x="NumOfUniqueProteinsPerGeneCell",
        y="NumOfGeneCellCouples",
        # color="Sample",
        # linestyle="NeuronalStrRep"
    )
    .facet(col="Sample")
    .layout(size=(width * 2, height))
    .add(so.Line())
    .scale(
        y="log",
        x=so.Continuous().tick(every=10),
        # linestyle=so.Nominal(order=["Neuronal", "Non-neuronal", "NA"])
    )
    .label(
        x="Unique isoforms per gene-cell couple",
        y="Gene-cell couples",
        # linestyle="Neuronal",
    )
    .limit(x=(0, None), y=(0, None))
)
p

# %%
width = 4.5
height = 4.5
p = (
    so.Plot(
        # ge5_neuronal_agg_gene_cell_couples_df.loc[
        #     ~ge5_neuronal_agg_gene_cell_couples_df["Neuronal"].isna()
        # ]
        ge5_neuronal_agg_gene_cell_couples_df,
        x="NumOfUniqueProteinsPerGeneCell",
        y="NumOfGeneCellCouples",
        # color="Sample",
        linestyle="NeuronalStrRep",
    )
    .facet(col="Sample")
    .layout(size=(width * 2, height))
    .add(so.Line())
    .scale(
        y="log",
        x=so.Continuous().tick(every=10),
        linestyle=so.Nominal(order=["Neuronal", "Non-neuronal", "NA"]),
    )
    .label(
        x="Unique isoforms per gene-cell couple",
        y="Gene-cell couples",
        # linestyle="Neuronal",
        linestyle="",
    )
    .limit(x=(0, None), y=(0, None))
)


p

# %% jupyter={"source_hidden": true}
# g = sns.catplot(
#     # data=ge5_per_sample_per_cb_copies_df.sort_values(["Neuronal", "Annotation"]),
#     data=ge5_per_sample_per_cb_copies_df.loc[~ge5_per_sample_per_cb_copies_df["Neuronal"].isna()],
#     x="Sample",
#     y="NumOfUniqueProteins",
#     # hue="NeuronalStrRep",
#     # hue_order=["Neuronal", "Non-neuronal"],
#     # kind="boxen",
#     kind="violin",
#     fill=False,
#     # split=True,
#     # sharex=False
#     height=3.5,
#     # height=8,
#     aspect=1.4,
#     # legend=False
# )
# g.set_axis_labels(y_var="Gene-cell couples with X unique isoforms")

# # # Remove duplicate legend entries (keep only the first one)
# # handles, labels = g.axes.flat[0].get_legend_handles_labels()  # Get handles from first facet
# # g._legend.set_title("")  # Alternative method if needed
# # g._legend.legendHandles = handles  # Ensure only unique handles appear

# g

# %% jupyter={"source_hidden": true}
# g = sns.catplot(
#     # data=ge5_per_sample_per_cb_copies_df.sort_values(["Neuronal", "Annotation"]),
#     data=ge5_per_sample_per_cb_copies_df.loc[~ge5_per_sample_per_cb_copies_df["Neuronal"].isna()],
#     x="Sample",
#     y="NumOfUniqueProteins",
#     hue="NeuronalStrRep",
#     hue_order=["Neuronal", "Non-neuronal"],
#     # kind="boxen",
#     kind="violin",
#     split=True,
#     # fill=False,
#     height=3.5,
#     # height=8,
#     aspect=1.4,
#     # legend=False
# )
# g.set_axis_labels(y_var="Gene-cell couples with X unique isoforms")

# # Remove duplicate legend entries (keep only the first one)
# handles, labels = g.axes.flat[0].get_legend_handles_labels()  # Get handles from first facet
# g._legend.set_title("")  # Alternative method if needed
# g._legend.legendHandles = handles  # Ensure only unique handles appear

# g

# %%
neuronal_order = ["Neuronal", "Non-neuronal", "NA"]

# Convert NeuronalStrRep to a categorical variable with the desired order
ge5_per_sample_per_cb_copies_df.loc[:, "CatNeuronalStrRep"] = pd.Categorical(
    ge5_per_sample_per_cb_copies_df["NeuronalStrRep"],
    categories=neuronal_order,
    ordered=True,
)

g = sns.catplot(
    data=ge5_per_sample_per_cb_copies_df.sort_values(
        ["CatNeuronalStrRep", "Annotation"]
    ),  # Sort the dataframe by NeuronalStrRep and Annotation
    x="NumOfUniqueProteins",
    y="Annotation",
    col="Sample",
    hue="CatNeuronalStrRep",
    hue_order=neuronal_order,
    kind="boxen",
    height=8,
    aspect=0.6,
)

del ge5_per_sample_per_cb_copies_df["CatNeuronalStrRep"]

g.set_axis_labels(x_var="Gene-cell couples with X unique isoforms", y_var="Cluster")
# Remove duplicate legend entries (keep only the first one)
handles, labels = g.axes.flat[
    0
].get_legend_handles_labels()  # Get handles from first facet
g._legend.set_title("")  # Alternative method if needed
g._legend.legendHandles = handles  # Ensure only unique handles appear

g

# %% jupyter={"source_hidden": true}
# # Count the number of observations per category
# col_order = ["Neuronal", "Non-neuronal", "NA"]
# # obs_counts = ge5_per_sample_per_cb_copies_df["NeuronalStrRep"].value_counts()
# obs_counts = ge5_per_sample_per_cb_copies_df.groupby("NeuronalStrRep")["Annotation"].nunique()
# col_widths = [obs_counts.get(col, 1) for col in col_order]  # Get count, default to 1

# # Normalize widths
# col_widths = np.array(col_widths) / sum(col_widths)  # Normalize to sum to 1

# # Create subplots with adjusted column widths
# fig, axes = plt.subplots(ncols=len(col_order), figsize=(12, 5),
#                          gridspec_kw={"width_ratios": col_widths}, sharey=True)

# # Store legend handles
# legend_handles_dict = {}

# # Loop through each category and plot separately
# for ax, col in zip(axes, col_order):
#     subset = ge5_per_sample_per_cb_copies_df[ge5_per_sample_per_cb_copies_df["NeuronalStrRep"] == col].sort_values("Annotation")
#     sns.violinplot(data=subset, x="Annotation", y="NumOfUniqueProteins", hue="Sample", split=True, ax=ax, legend=False, fill=False)  # Disable individual legends
#     ax.set_title(col)
#     # Remove individual x-axis labels
#     ax.set_xlabel("")
#     # Remove individual y-axis labels
#     ax.set_ylabel("")
#     # Despine each subplot
#     sns.despine(ax=ax)
#     # Rotate x-axis labels safely
#     plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha="right")
#     # Capture legend handles only once (from the first plot)
#     if not legend_handles_dict:
#         handles, labels = ax.get_legend_handles_labels()
#         legend_handles_dict = {label: h for h, label in zip(handles, labels)}

# # Reduce gaps between facets
# plt.subplots_adjust(wspace=0.03)  # Adjust horizontal spacing

# # Create a single legend
# if legend_handles_dict:
#     legend_patches = [mpatches.Patch(color=h.get_facecolor()[0], label=label) for label, h in legend_handles_dict.items()]
#     fig.legend(legend_patches, legend_handles_dict.keys(), title="Sample", loc="upper center", bbox_to_anchor=(0.5, 1.1), ncol=len(legend_handles_dict))

# # Add a single x-axis title
# fig.supxlabel("Cluster")
# # Add a single y-axis title
# fig.supylabel("Gene-cell couples with X unique isoforms")

# plt.tight_layout()
# plt.show()

# %%

# %% jupyter={"source_hidden": true}
fig = px.line(
    agg_gene_cell_couples_df,
    x="NumOfUniqueProteinsPerGeneCell",
    y="NumOfUniqueCells",
    color="Sample",
    markers=True,
    # facet_col="Sample",
    # color_discrete_map=color_discrete_map,
    # labels={"NumOfUniqueProteins": "Min unique isoforms per gene-cell couple"},
    # histnorm="percent",
    # cumulative=True
    # log_x=True,
    log_y=True,
)
fig.update_traces(opacity=0.75, marker_size=4)
# fig.update_traces(cumulative_direction="decreasing", selector=dict(type='histogram'))
fig.update_xaxes(tick0=0, dtick=10, title="Unique isoforms per gene-cell couple")
fig.update_yaxes(
    title="Unique cells in gene-cell couples",
    # range=[0, 6]
)
fig.update_layout(
    width=500,
    height=350,
    template=template,
    # barmode="overlay",
    # title="Num of unique distinct isoforms per sample-gene-cell",
    # showlegend=False,
)
fig.show()

# %% jupyter={"source_hidden": true}
fig = px.line(
    agg_gene_cell_couples_df,
    x="NumOfUniqueProteinsPerGeneCell",
    y="NumOfUniqueGenes",
    color="Sample",
    markers=True,
    # facet_col="Sample",
    # color_discrete_map=color_discrete_map,
    # labels={"NumOfUniqueProteins": "Min unique isoforms per gene-cell couple"},
    # histnorm="percent",
    # cumulative=True
    # log_x=True,
    log_y=True,
)
fig.update_traces(opacity=0.75, marker_size=4)
# fig.update_traces(cumulative_direction="decreasing", selector=dict(type='histogram'))
fig.update_xaxes(tick0=0, dtick=10, title="Unique isoforms per gene-cell couple")
fig.update_yaxes(
    title="Unique genes in gene-cell couples",
    # range=[0, 6]
)
fig.update_layout(
    width=500,
    height=350,
    template=template,
    # barmode="overlay",
    # title="Num of unique distinct isoforms per sample-gene-cell",
    # showlegend=False,
)
fig.show()

# %% jupyter={"source_hidden": true}
fig = px.line(
    agg_gene_cell_couples_df,
    x="NumOfUniqueProteinsPerGeneCell",
    y="AvgNumOfUniqueProteinsPerGene",
    color="Sample",
    markers=True,
    # facet_col="Sample",
    # color_discrete_map=color_discrete_map,
    # labels={"NumOfUniqueProteins": "Min unique isoforms per gene-cell couple"},
    # histnorm="percent",
    # cumulative=True
    # log_x=True,
    # log_y=True
)
fig.update_traces(opacity=0.75, marker_size=4)
# fig.update_traces(cumulative_direction="decreasing", selector=dict(type='histogram'))
fig.update_xaxes(tick0=0, dtick=10, title="Unique isoforms per gene-cell couple")
fig.update_yaxes(
    title="Unique proteins per gene<br>in gene-cell couples (avg)",
    # range=[0, 6]
    # tick0=0,
    dtick=10,
)
fig.update_layout(
    width=500,
    height=350,
    template=template,
    # barmode="overlay",
    # title="Num of unique distinct isoforms per sample-gene-cell",
    # showlegend=False,
)
fig.show()

# %% jupyter={"source_hidden": true}
fig = px.line(
    agg_gene_cell_couples_df,
    x="NumOfUniqueProteinsPerGeneCell",
    y="DecreasingCummlativeGeneCellCouples%",
    color="Sample",
    markers=True,
    # facet_col="Sample",
    # color_discrete_map=color_discrete_map,
    # labels={"NumOfUniqueProteins": "Min unique isoforms per gene-cell couple"},
    # histnorm="percent",
    # cumulative=True
    # log_y=True
)
fig.update_traces(opacity=0.75, marker_size=4)
# fig.update_traces(cumulative_direction="decreasing", selector=dict(type='histogram'))
fig.update_xaxes(tick0=0, dtick=10, title="Min unique isoforms per gene-cell couple")
fig.update_yaxes(
    title="Gene-cell couples [%]",
    # range=[0, 6]
)
fig.update_layout(
    width=500,
    height=350,
    template=template,
    # barmode="overlay",
    # title="Num of unique distinct isoforms per sample-gene-cell",
    # showlegend=False,
)
fig.show()

# %% jupyter={"source_hidden": true}
# # width = 4.5
# width = 4
# height = 4
# p = (
#     so.Plot(ge5_agg_gene_cell_couples_df, x="NumOfUniqueProteinsPerGeneCell", y="DecreasingCummlativeGeneCellCouples%")
#     .layout(size=(width, height))
#     .add(so.Line(), color="Sample")
#     # .scale(
#     #     x="log",
#     #     y=so.Continuous().tick(every=20_000).label(unit=""),
#     # )
#     .label(x="Min unique isoforms per gene-cell couple", y="Gene-cell couples [%]")
#     .limit(x=(0, None), y=(0, 100))
# )
# p

# %%
# width = 4.5
width = 4.5
height = 4.5
p = (
    so.Plot(
        ge5_agg_gene_cell_couples_df,
        x="NumOfUniqueProteinsPerGeneCell",
        y="DecreasingCummlativeGeneCellCouples%",
    )
    .layout(size=(width * 2, height))
    .add(so.Line())
    .facet(col="Sample")
    # .scale(
    #     x="log",
    #     y=so.Continuous().tick(every=20_000).label(unit=""),
    # )
    .label(x="Min unique isoforms per gene-cell couple", y="Gene-cell couples [%]")
    .limit(x=(0, None), y=(0, 100))
)
p

# %%
width = 4.5
height = 4.5
p = (
    so.Plot(
        # ge5_neuronal_agg_gene_cell_couples_df.loc[
        #     ~ge5_neuronal_agg_gene_cell_couples_df["Neuronal"].isna()
        # ],
        ge5_neuronal_agg_gene_cell_couples_df,
        x="NumOfUniqueProteinsPerGeneCell",
        y="DecreasingCummlativeGeneCellCouples%",
        # color="Sample",
        linestyle="NeuronalStrRep",
    )
    .facet(col="Sample")
    .layout(size=(width * 2, height))
    .add(so.Line())
    .scale(linestyle=so.Nominal(order=["Neuronal", "Non-neuronal", "NA"]))
    .label(
        x="Min unique isoforms per gene-cell couple",
        y="Gene-cell couples [%]",
        linestyle="",
    )
    .limit(x=(0, None), y=(0, 100))
)
p


# %%
# # width = 4.5
# width = 4
# height = 3.5
# p = (
#     so.Plot(
#         ge5_neuronal_agg_gene_cell_couples_df,
#         x="NumOfUniqueProteinsPerGeneCell",
#         y="DecreasingCummlativeGeneCellCouples%",
#         color="Sample",
#         linestyle="NeuronalStrRep"
#     )
#     # .facet(col="Sample", row="NeuronalStrRep")
#     .layout(size=(width, height))
#     .add(so.Line())
#     # .scale(marker=so.Nominal(order=["Neuronal", "Non-neuronal"]))
#     .label(
#         x="Min unique isoforms per gene-cell couple",
#         y="Gene-cell couples [%]",
#         linestyle="Neuronal",
#     )
#     .limit(x=(0, None), y=(0, 100))
# )
# p

# %%
# # width = 4.5
# width = 4
# height = 3.5
# p = (
#     so.Plot(
#         ge5_neuronal_agg_gene_cell_couples_df.loc[
#             ~ge5_neuronal_agg_gene_cell_couples_df["Neuronal"].isna()
#         ],
#         x="NumOfUniqueProteinsPerGeneCell",
#         y="DecreasingCummlativeGeneCellCouples%",
#         color="Sample",
#         linestyle="NeuronalStrRep"
#     )
#     .facet(col="Sample", row="NeuronalStrRep")
#     .layout(size=(width * 2, height * 2))
#     .add(so.Line(), legend=False)
#     .scale(marker=so.Nominal(order=["Neuronal", "Non-neuronal"]))
#     .label(
#         x="Min unique isoforms per gene-cell couple",
#         y="Gene-cell couples [%]",
#         linestyle="Neuronal",
#     )
#     .limit(x=(0, None), y=(0, 100))
# )
# p

# %%
# # width = 4.5
# width = 4
# height = 3.5
# p = (
#     so.Plot(
#         ge5_neuronal_agg_gene_cell_couples_df.loc[
#             ~ge5_neuronal_agg_gene_cell_couples_df["Neuronal"].isna()
#         ],
#         x="NumOfUniqueProteinsPerGeneCell",
#         y="DecreasingCummlativeGeneCellCouples%",
#         color="Sample",
#         # linestyle="NeuronalStrRep"
#     )
#     .facet(col="Sample", row="NeuronalStrRep")
#     .layout(size=(width * 2, height * 2))
#     # .add(so.Bars())
#     .add(so.Area(), legend=False)
#     .scale(marker=so.Nominal(order=["Neuronal", "Non-neuronal"]))
#     .label(
#         x="Min unique isoforms per gene-cell couple",
#         y="Gene-cell couples [%]",
#         linestyle="Neuronal",
#     )
#     .limit(x=(0, None), y=(0, 100))
# )
# p

# %%
# # width = 4.5
# width = 4
# height = 4
# p = (
#     so.Plot(ge5_neuronal_agg_gene_cell_couples_df, x="NumOfUniqueProteinsPerGeneCell", y="DecreasingCummlativeGeneCellCouples%")
#     # .facet(row="Neuronal")
#     .facet(row="NeuronalStrRep", order={"NeuronalStrRep": ["Neuronal", "Non-neuronal", "NA"]})
#     .layout(size=(width, height*3))
#     .add(so.Line(), color="Sample")
#     # .scale(
#     #     x="log",
#     #     y=so.Continuous().tick(every=20_000).label(unit=""),
#     # )
#     .label(x="Min unique isoforms per gene-cell couple", y="Gene-cell couples [%]")
#     .limit(x=(0, None), y=(0, 100))
# )
# p

# %%

# %%

# %% [markdown]
# #### 2+ reads per cell

# %%
def make_agg_gene_cell_couples_x_plus_copies_df(
    per_sample_per_cb_copies_df, x_plus_copies
):
    per_sample_per_cb_x_plus_copies_df = per_sample_per_cb_copies_df.loc[
        per_sample_per_cb_copies_df["Copies"] >= x_plus_copies
    ]

    agg_gene_cell_couples_df_1 = (
        per_sample_per_cb_x_plus_copies_df.groupby(["Sample", "NumOfUniqueProteins"])
        .agg(
            NumOfGeneCellCouples=("Chrom", "size"),
            NumOfUniqueGenes=("Chrom", "nunique"),
            NumOfUniqueCells=("CB", "nunique"),
        )
        .reset_index()
        .rename(columns={"NumOfUniqueProteins": "NumOfUniqueProteinsPerGeneCell"})
    )
    # agg_gene_cell_couples_df_1

    agg_gene_cell_couples_df_2 = (
        per_sample_per_cb_x_plus_copies_df.groupby(
            ["Sample", "NumOfUniqueProteins", "Chrom"]
        )["UniqueProteins"]
        .apply(lambda x: len(set(chain.from_iterable(x))))
        .reset_index(name="NumOfUniqueProteinsPerChrom")
        .groupby(["Sample", "NumOfUniqueProteins"])["NumOfUniqueProteinsPerChrom"]
        .apply(list)
        .reset_index(name="NumOfUniqueProteinsPerChroms")
        .rename(columns={"NumOfUniqueProteins": "NumOfUniqueProteinsPerGeneCell"})
    )
    agg_gene_cell_couples_df_2["AvgNumOfUniqueProteinsPerGene"] = (
        agg_gene_cell_couples_df_2["NumOfUniqueProteinsPerChroms"].apply(
            lambda x: sum(x) / len(x)
        )
    )
    agg_gene_cell_couples_df_2["NumOfUniqueProteinsInAllGenes"] = (
        agg_gene_cell_couples_df_2["NumOfUniqueProteinsPerChroms"].apply(sum)
    )
    del agg_gene_cell_couples_df_2["NumOfUniqueProteinsPerChroms"]
    # agg_gene_cell_couples_df_2

    agg_gene_cell_couples_df = agg_gene_cell_couples_df_1.merge(
        agg_gene_cell_couples_df_2
    )
    del agg_gene_cell_couples_df_1, agg_gene_cell_couples_df_2

    samples = ["nuclearRNA", "totalRNA"]
    dfs = []
    for sample in samples:
        df = agg_gene_cell_couples_df.loc[
            agg_gene_cell_couples_df["Sample"] == sample
        ].sort_values("NumOfUniqueProteinsPerGeneCell", ascending=False)
        # df["DecreasingCummlativeGeneCellCouples%"] = df["NumOfGeneCellCouples"].mul(100).div(df["NumOfGeneCellCouples"].sum()).cumsum()
        df["GeneCellCouples%"] = (
            df["NumOfGeneCellCouples"].mul(100).div(df["NumOfGeneCellCouples"].sum())
        )
        df["DecreasingCummlativeGeneCellCouples%"] = df["GeneCellCouples%"].cumsum()
        dfs.append(df)
    agg_gene_cell_couples_df = pd.concat(dfs).sort_values(
        ["Sample", "NumOfUniqueProteinsPerGeneCell"], ignore_index=True
    )

    agg_gene_cell_couples_df.insert(
        agg_gene_cell_couples_df.columns.get_loc("Sample") + 1,
        # 0,
        "MinNumOfReadsPerGeneCell",
        x_plus_copies,
    )

    return agg_gene_cell_couples_df


# %%
agg_gene_cell_couples_x_plus_copies_dfs = [
    make_agg_gene_cell_couples_x_plus_copies_df(
        per_sample_per_cb_copies_df, x_plus_copies
    )
    for x_plus_copies in range(1, 6)
]
concat_agg_gene_cell_couples_x_plus_copies_df = (
    pd.concat(agg_gene_cell_couples_x_plus_copies_dfs).sort_values(
        ["Sample", "MinNumOfReadsPerGeneCell", "NumOfUniqueProteinsPerGeneCell"],
        ignore_index=True,
    )
    # .sort_values(["MinNumOfReadsPerGeneCell", "Sample", "NumOfUniqueProteinsPerGeneCell"], ignore_index=True)
)
concat_agg_gene_cell_couples_x_plus_copies_df

# %%
len(agg_gene_cell_couples_x_plus_copies_dfs)

# %%
agg_gene_cell_couples_x_plus_copies_dfs[1]

# %%
agg_gene_cell_couples_x_plus_copies_dfs[4]

# %%
concat_agg_gene_cell_couples_x_plus_copies_df.drop(
    columns="DecreasingCummlativeGeneCellCouples%"
).to_csv(
    Path(max_expression_dir, "ConcatAggGeneCellCouplesXPlusCopies.tsv"),
    sep="\t",
    index=False,
)

# %%
(
    concat_agg_gene_cell_couples_x_plus_copies_df.groupby(
        ["Sample", "MinNumOfReadsPerGeneCell"]
    )
    # .agg({
    #     "NumOfGeneCellCouples": "sum"
    # })
    .agg(TotalNumOfGeneCellCouples=("NumOfGeneCellCouples", "sum"))
)

# %%
concat_agg_gene_cell_couples_x_plus_copies_df.loc[
    concat_agg_gene_cell_couples_x_plus_copies_df["NumOfUniqueProteinsPerGeneCell"] <= 5
].drop(columns=["DecreasingCummlativeGeneCellCouples%"])

# %%

# %%
# todo retain only True and use this with merge to filter the original table
per_sample_per_cb_copies_df.loc[(per_sample_per_cb_copies_df["Copies"] >= 2)].groupby(
    ["Chrom", condition_col, "Sample"]
)["NumOfUniqueProteins"].apply(lambda x: x.eq(1).all()).reset_index()

# %%
# todo retain only True and use this with merge to filter the original table
per_sample_per_cb_copies_df.loc[(per_sample_per_cb_copies_df["Copies"] >= 2)].groupby(
    ["Chrom", condition_col, "Sample"]
)["NumOfUniqueProteins"].apply(lambda x: x.eq(1).all()).reset_index(name="")

# %% [markdown] jp-MarkdownHeadingCollapsed=true
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
expanded_max_distinct_proteins_df = max_distinct_proteins_df.copy()
expanded_max_distinct_proteins_df["Proteins"] = expanded_max_distinct_proteins_df[
    "Proteins"
].str.split(",")
expanded_max_distinct_proteins_df = (
    expanded_max_distinct_proteins_df.explode("Proteins")
    .reset_index(drop=True)
    .rename(columns={"Proteins": "Protein"})
)
expanded_max_distinct_proteins_df = expanded_max_distinct_proteins_df.drop(
    [
        "NumOfProteins",
        "NumOfReads",
        "MappedReads",
        "Samples",
        "MappedReadsPerSample",
        "AvailableReads",
        "NumOfAvailableReads",
        # "DistinctProteins/Reads",
    ],
    axis=1,
)

# expanded_max_distinct_proteins_df

strongly_diversified_expanded_max_distinct_proteins_dfs = [
    (
        expanded_max_distinct_proteins_df.loc[
            expanded_max_distinct_proteins_df["Chrom"] == chrom
        ].merge(
            expanded_unique_proteins_df,
            on=["Chrom", "Transcript", "Protein"],
            how="left",
        )
        # .merge(samples_and_tissues_df, on="Sample", how="left")
    )
    for chrom, expanded_unique_proteins_df in zip(chroms, expanded_unique_proteins_dfs)
    if chrom in strongly_diversified_chroms
]

del expanded_max_distinct_proteins_df

# for (
#     strongly_diversified_expanded_max_distinct_proteins_df
# ) in strongly_diversified_expanded_max_distinct_proteins_dfs:
#     strongly_diversified_expanded_max_distinct_proteins_df.rename(
#         columns={"Sample": "Sample2", "Protein": "Protein2", "Tissue": "Tissue2"},
#         inplace=True,
#     )
#     strongly_diversified_expanded_max_distinct_proteins_df.insert(
#         2, "Sample", strongly_diversified_expanded_max_distinct_proteins_df["Sample2"]
#     )
#     strongly_diversified_expanded_max_distinct_proteins_df.insert(
#         3, "Tissue", strongly_diversified_expanded_max_distinct_proteins_df["Tissue2"]
#     )
#     strongly_diversified_expanded_max_distinct_proteins_df.insert(
#         4, "Protein", strongly_diversified_expanded_max_distinct_proteins_df["Protein2"]
#     )
#     strongly_diversified_expanded_max_distinct_proteins_df.drop(
#         ["Sample2", "Tissue2", "Protein2"], axis=1, inplace=True
#     )

# expanded_max_distinct_proteins_dfs[1]
strongly_diversified_expanded_max_distinct_proteins_dfs[0].loc[
    :, :"MaxNonSynsFrequency"
]

# %%
samples_and_tissues_df

# %%
strongly_diversified_expanded_max_distinct_proteins_dfs[0]["Protein"].unique().size

# %%
strongly_diversified_transcripts

# %%
# for strongly_diversified_expanded_max_distinct_proteins_df in strongly_diversified_expanded_max_distinct_proteins_dfs:
#     print(strongly_diversified_expanded_max_distinct_proteins_df.loc[:3, :"KnownSites"])

# %%
strongly_diversified_num_of_proteins_per_sample_dfs = []

for (
    strongly_diversified_expanded_max_distinct_proteins_df,
    # num_of_unique_proteins,
) in zip(
    strongly_diversified_expanded_max_distinct_proteins_dfs,
    # strongly_diversified_max_num_of_proteins,
):
    gb = strongly_diversified_expanded_max_distinct_proteins_df.groupby(
        ["Chrom", "Transcript", "Sample", "Tissue"]
    )

    df = gb.apply(len).reset_index().rename(columns={0: "NumOfProteins"})

    num_of_unique_proteins = (
        strongly_diversified_expanded_max_distinct_proteins_df["Protein"].unique().size
    )

    df2 = (
        gb.apply(lambda x: 100 * len(x) / num_of_unique_proteins)
        .reset_index()
        .rename(columns={0: "%RelativeNumOfProteins"})
    )
    df = df.merge(df2)

    strongly_diversified_num_of_proteins_per_sample_dfs.append(df)

strongly_diversified_num_of_proteins_per_sample_dfs[0]

# %%
# total_non_unique_proteins_per_strongly_diversified_chroms = [
#     len(strongly_diversified_expanded_max_distinct_proteins_df)
#     for strongly_diversified_expanded_max_distinct_proteins_df in strongly_diversified_expanded_max_distinct_proteins_dfs
# ]
# total_non_unique_proteins_per_strongly_diversified_chroms

# %%
# total_non_unique_proteins_per_strongly_diversified_chroms_2 = [
#     strongly_diversified_num_of_proteins_per_sample_df["NumOfProteins"].sum()
#     for strongly_diversified_num_of_proteins_per_sample_df in strongly_diversified_num_of_proteins_per_sample_dfs
# ]
# total_non_unique_proteins_per_strongly_diversified_chroms_2

# %%
# strongly_diversified_max_num_of_proteins

# %%
len(per_transcript_per_sample_coverage_dfs)

# %%
len(possibly_na_chroms)

# %%
len(possibly_na_positions_files)

# %%
strongly_diversified_per_transcript_per_sample_coverage_dfs[0]["Chrom"].iloc[0]

# %%
strongly_diversified_per_transcript_per_sample_coverage_dfs = [
    df
    for df in per_transcript_per_sample_coverage_dfs
    if df["Chrom"].iloc[0] in strongly_diversified_chroms
]
ic(len(strongly_diversified_per_transcript_per_sample_coverage_dfs))
strongly_diversified_per_transcript_per_sample_coverage_dfs[0]

# %%
concat_strongly_diversified_per_transcript_per_sample_coverage_df = pd.concat(
    strongly_diversified_per_transcript_per_sample_coverage_dfs
)
concat_strongly_diversified_per_transcript_per_sample_coverage_df

# %%
strongly_diversified_num_of_reads_and_proteins_per_sample_dfs = []
for (
    strongly_diversified_num_of_proteins_per_sample_df
) in strongly_diversified_num_of_proteins_per_sample_dfs:
    strongly_diversified_num_of_reads_and_proteins_per_sample_df = (
        strongly_diversified_num_of_proteins_per_sample_df.merge(
            concat_strongly_diversified_per_transcript_per_sample_coverage_df,
            on=["Chrom", "Transcript", "Sample", "Tissue"],
            how="left",
        )
    )
    assert (
        strongly_diversified_num_of_reads_and_proteins_per_sample_df["NumOfReads"]
        .ge(
            strongly_diversified_num_of_reads_and_proteins_per_sample_df[
                "NumOfProteins"
            ]
        )
        .all()
    )
    strongly_diversified_num_of_reads_and_proteins_per_sample_df[
        "NumOfReads/NumOfProteins"
    ] = (
        strongly_diversified_num_of_reads_and_proteins_per_sample_df["NumOfReads"]
        / strongly_diversified_num_of_reads_and_proteins_per_sample_df["NumOfProteins"]
    )
    strongly_diversified_num_of_reads_and_proteins_per_sample_dfs.append(
        strongly_diversified_num_of_reads_and_proteins_per_sample_df
    )
concat_strongly_diversified_num_of_reads_and_proteins_per_sample_df = pd.concat(
    strongly_diversified_num_of_reads_and_proteins_per_sample_dfs
)
concat_strongly_diversified_num_of_reads_and_proteins_per_sample_df

# %%
fig = px.box(
    concat_strongly_diversified_num_of_reads_and_proteins_per_sample_df,
    x="Tissue",
    y="NumOfProteins",
    color="Tissue",
    color_discrete_map=tissues_color_discrete_map,
    points="all",
    category_orders={"Tissue": tissues_order},
)
fig.update_layout(width=900, height=500, template=template, showlegend=False)
fig.show()

# %%
fig = px.box(
    concat_strongly_diversified_num_of_reads_and_proteins_per_sample_df,
    x="Tissue",
    y="NumOfReads/NumOfProteins",
    color="Tissue",
    color_discrete_map=tissues_color_discrete_map,
    points="all",
    category_orders={"Tissue": tissues_order},
)
fig.update_yaxes(
    range=[
        0,
        concat_strongly_diversified_num_of_reads_and_proteins_per_sample_df[
            "NumOfReads/NumOfProteins"
        ].max()
        * 1.4,
    ],
    tick0=0,
    dtick=1,
    zeroline=True,
    zerolinewidth=2,
)
fig.update_layout(width=900, height=500, template=template, showlegend=False)

fig.show()

# %%
concat_strongly_diversified_num_of_reads_and_proteins_per_sample_df["NumOfReads"].max()

# %%
concat_strongly_diversified_num_of_reads_and_proteins_per_sample_df[
    "NumOfProteins"
].max()

# %%
strongly_diversified_chroms

# %%
strongly_diversified_max_num_of_proteins

# %%
strongly_diversified_transcripts

# %%
ceil(300)

# %%
cols = min(facet_col_wrap, len(strongly_diversified_num_of_proteins_per_sample_dfs), 3)
rows = ceil(len(strongly_diversified_num_of_proteins_per_sample_dfs) / cols)
row_col_iter = list(product(range(1, rows + 1), range(1, cols + 1)))[
    : len(strongly_diversified_num_of_proteins_per_sample_dfs)
]

x_title = "Coverage"
y_title = "Distinct protein isoforms"


subplot_titles = [
    f"{transcript.split('_')[0]} ({int(total_pooled_isoforms)})"
    for transcript, total_pooled_isoforms in zip(
        strongly_diversified_transcripts, strongly_diversified_max_num_of_proteins
    )
]


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
    strongly_diversified_chrom,
    strongly_diversified_max_num_of_protein,
) in zip(
    row_col_iter,
    strongly_diversified_chroms,
    strongly_diversified_max_num_of_proteins,
):
    strongly_diversified_num_of_reads_and_proteins_per_sample_df = (
        concat_strongly_diversified_num_of_reads_and_proteins_per_sample_df.loc[
            concat_strongly_diversified_num_of_reads_and_proteins_per_sample_df["Chrom"]
            == strongly_diversified_chrom
        ]
    )

    _tissues = strongly_diversified_num_of_reads_and_proteins_per_sample_df["Tissue"]

    for tissue in _tissues:
        try:
            x = strongly_diversified_num_of_reads_and_proteins_per_sample_df.loc[
                strongly_diversified_num_of_reads_and_proteins_per_sample_df["Tissue"]
                == tissue,
                "NumOfReads",
            ]
        except KeyError:
            x = [0]

        y = strongly_diversified_num_of_reads_and_proteins_per_sample_df.loc[
            strongly_diversified_num_of_reads_and_proteins_per_sample_df["Tissue"]
            == tissue,
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

max_x_y = ceil(max(max_x, max_y)) * 1.05

fig.update_xaxes(
    # showticklabels=False,  # Hide x axis ticks
    # categoryorder="array",
    # categoryarray=tissues_order,
    # range=[0, max_x * 1.1],
    range=[0, max_x_y],
    tick0=0,
    # dtick=50
)
# fig.update_yaxes(
#     range=[0, max_y * 1.1],
#     # range=[0, 0.5+np.log(max_y)/np.log(10)],
#     # type="log"
# )
fig.update_yaxes(
    # title_text=primary_y_title,
    # range=[0, max_y * 1.1],
    range=[0, max_x_y],
    tick0=0,
    # dtick=50
    # secondary_y=False
)
# fig.update_yaxes(
#     title_text=secondary_y_title,
#     range=[0, max_y_2 * 1.1],
#     secondary_y=True
# )

width = 950
height = 800

fig.update_traces(opacity=0.7, marker_size=8)

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


# %%

# %%

# %%

# %%

# %%

# %%
# cols = min(facet_col_wrap, len(strongly_diversified_num_of_proteins_per_sample_dfs), 3)
# rows = ceil(len(strongly_diversified_num_of_proteins_per_sample_dfs) / cols)
# row_col_iter = list(product(range(1, rows + 1), range(1, cols + 1)))[
#     : len(strongly_diversified_num_of_proteins_per_sample_dfs)
# ]

# x_title = "Coverage"
# y_title = "Distinct protein isoforms"

# # title_text = "Distribution of min & max estimates of non-syn substitutions per read"


# # subplot_titles = [
# #     f"{transcript.split('_')[0]}<br><sub>Pooled distinct proteins = {int(total_pooled_isoforms)}</sub>"
# #     for transcript, total_pooled_isoforms in zip(
# #         strongly_diversified_transcripts, strongly_diversified_max_num_of_proteins
# #     )
# # ]

# subplot_titles = [
#     f"{transcript.split('_')[0]} ({int(total_pooled_isoforms)})"
#     for transcript, total_pooled_isoforms in zip(
#         strongly_diversified_transcripts, strongly_diversified_max_num_of_proteins
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
#     # horizontal_spacing=0.2,
# )

# max_x = 0
# max_y = 0
# legend_constructed = False

# for (
#     (row, col),
#     strongly_diversified_num_of_proteins_per_sample_df,
#     strongly_diversified_max_num_of_protein,
#     strongly_diversified_per_transcript_per_sample_coverage_df,
# ) in zip(
#     row_col_iter,
#     strongly_diversified_num_of_proteins_per_sample_dfs,
#     strongly_diversified_max_num_of_proteins,
#     strongly_diversified_per_transcript_per_sample_coverage_dfs,
# ):
#     _tissues = strongly_diversified_num_of_proteins_per_sample_df["Tissue"]

#     for tissue in _tissues:
#         try:
#             x = strongly_diversified_per_transcript_per_sample_coverage_df.loc[
#                 strongly_diversified_per_transcript_per_sample_coverage_df["Tissue"]
#                 == tissue,
#                 "NumOfReads",
#             ]
#         except KeyError:
#             x = [0]

#         y = strongly_diversified_num_of_proteins_per_sample_df.loc[
#             strongly_diversified_num_of_proteins_per_sample_df["Tissue"] == tissue,
#             "NumOfProteins",
#         ]
#         # ic(samp)

#         max_x = max(max_x, x.max())
#         max_y = max(max_y, y.max())

#         if not legend_constructed:
#             fig.add_trace(
#                 go.Scatter(
#                     x=x,
#                     y=y,
#                     marker_color=tissues_color_discrete_map[tissue],
#                     name=tissue,
#                     legendrank=tissue_to_legendrank[tissue],
#                     # marker_pattern_shape="/",
#                     mode="markers",
#                 ),
#                 row=row,
#                 col=col,
#             )
#         else:
#             fig.add_trace(
#                 go.Scatter(
#                     x=x,
#                     y=y,
#                     marker_color=tissues_color_discrete_map[tissue],
#                     # name=tissue,
#                     showlegend=False,
#                     # marker_pattern_shape="/",
#                     mode="markers",
#                 ),
#                 row=row,
#                 col=col,
#             )

#     legend_constructed = True


# # fig.update_xaxes(tickangle=45, automargin=True)

# fig.update_xaxes(
#     # showticklabels=False,  # Hide x axis ticks
#     # categoryorder="array",
#     # categoryarray=tissues_order,
#     range=[0, max_x * 1.1]
# )
# # fig.update_yaxes(
# #     range=[0, max_y * 1.1],
# #     # range=[0, 0.5+np.log(max_y)/np.log(10)],
# #     # type="log"
# # )
# fig.update_yaxes(
#     # title_text=primary_y_title,
#     range=[0, max_y * 1.1],
#     # secondary_y=False
# )
# # fig.update_yaxes(
# #     title_text=secondary_y_title,
# #     range=[0, max_y_2 * 1.1],
# #     secondary_y=True
# # )

# width = 950
# height = 800

# fig.update_traces(opacity=0.7, marker_size=6)

# fig.update_layout(
#     template=template,
#     # title_text="Octopus",
#     # title_x=0.1,
#     # title_y=0.97,
#     # showlegend=False,
#     legend_title_text="Tissue",
#     width=width,
#     height=height,
#     # barmode="overlay",
# )

# fig.write_image(
#     "Distinct proteins vs coverage per 9 strong transcripts - Octopus.svg",
#     width=width,
#     height=height,
# )

# fig.show()


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

max_distinct_proteins_per_transcript_and_alg_df["MaxNumOfProteins"] = (
    max_distinct_proteins_per_transcript_and_alg_df.groupby(
        [condition_col, "Algorithm"]
    )["NumOfProteins"].transform(max)
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

# %%
dispersion_df["%SolutionsDispersion"].describe().round(1)

# %%
scipy.stats.iqr(dispersion_df["%SolutionsDispersion"])

# %%
round(scipy.stats.iqr(dispersion_df["%SolutionsDispersion"]), 2)

# %%
np.percentile(dispersion_df["%SolutionsDispersion"], [25, 75])

# %%
np.round(np.percentile(dispersion_df["%SolutionsDispersion"], [25, 75]), 1)

# %%
dispersion_df["%SolutionsDispersion"].size

# %%
# percent_of_octopus_genes_with_nonzero_dispersion
np.round(
    100
    * len(dispersion_df.loc[dispersion_df["%SolutionsDispersion"] > 0])
    / len(dispersion_df),
    1,
)

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
fig.update_yaxes(title="Genes", type="log")

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
fig = go.Figure(
    go.Histogram(
        y=dispersion_df["%SolutionsDispersion"],
        marker_color="black",
    )
)

# fig.update_xaxes(title="% dispersion<br><sub>100 * (max - min) / max</sub>")
fig.update_yaxes(
    title="Dispersion [%]",
    # type="log"
)
fig.update_xaxes(title="Genes", type="log")

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
grouped_dispersion_df = (
    dispersion_df.groupby("%SolutionsDispersion")
    .size()
    .reset_index()
    .rename(columns={0: "Genes"})
)
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


# %% [markdown] toc-hr-collapsed=true
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
