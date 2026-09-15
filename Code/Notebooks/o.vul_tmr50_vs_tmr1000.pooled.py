# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
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

# %% [markdown] papermill={"duration": 0.029907, "end_time": "2022-02-01T09:42:43.198426", "exception": false, "start_time": "2022-02-01T09:42:43.168519", "status": "completed"}
# # Imports

# %%
code_dir = "/private7/projects/Combinatorics/Code"

out_dir = "/private7/projects/Combinatorics/Code/Notebooks"

# %%
# %load_ext autoreload
# %autoreload 2
# # %autosave 600

# %% papermill={"duration": 2.901153, "end_time": "2022-02-01T09:42:46.125355", "exception": false, "start_time": "2022-02-01T09:42:43.224202", "status": "completed"}
import sys
from functools import reduce
from itertools import chain, combinations, product
import math
from math import ceil
from multiprocessing import Pool
from pathlib import Path
from random import choice
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.colors as pc
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio

from scipy import interpolate  # todo unimport this later?
from scipy.stats import fisher_exact, chi2_contingency
from statsmodels.stats.multitest import fdrcorrection, multipletests
from statsmodels.stats.proportion import binom_test
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

# %%
pd.set_option('display.max_columns', 50)

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

main_data_dir = Path(
    "/private6/projects/Combinatorics/O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.PooledSamples"
)

positions_dir = Path(main_data_dir, "PositionsFiles")
reads_dir = Path(main_data_dir, "ReadsFiles")
proteins_dir = Path(main_data_dir, "ProteinsFiles")
distinct_proteins_dir = Path(main_data_dir, "DistinctProteins")
expression_dir = Path(main_data_dir, "ExpressionLevels")

neural_vs_non_neural_expression_file = Path(
    "/private7/projects/Combinatorics/O.vulgaris/Annotations/NeuralVsNonNeuralExpression.csv"
)

samples_and_tissues_file = Path(
    "/private7/projects/Combinatorics/O.vulgaris/Data/PRJNA791920/IsoSeqPolished/samples.csv"
)

reads_first_col_pos = 7
# unique_reads_first_col_pos = 9
# proteins_first_col_pos = 13
unique_proteins_first_col_pos = 15
reads_editing_col = "EditingFrequency"
proteins_editing_col = "MinNonSyns"
max_snps_per_gene_to_allow_editing_detection = 3
snp_noise_level = 0.05


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
tmr50_alignment_stats_df = alignment_stats_df.loc[
    alignment_stats_df["MappedReads"] >= 50
]
tmr50_alignment_stats_df

# %%
tmr50_alignment_stats_df.loc[
    tmr50_alignment_stats_df["MappedReads"].ge(100)
]

# %%
tmr50_alignment_stats_df.loc[
    tmr50_alignment_stats_df["MappedReads"].ge(500)
]

# %%
tmr50_alignment_stats_df.loc[
    tmr50_alignment_stats_df["MappedReads"].ge(1000)
]

# %%
fig = px.histogram(
    tmr50_alignment_stats_df,
    x="MappedReads"
)

fig.show()

# %%
# tmr1000_alignment_stats_df = alignment_stats_df.loc[
#     alignment_stats_df["MappedReads"] >= 500
# ]
# tmr1000_alignment_stats_df

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

# snps_positions_files = list(positions_dir.glob("*.positions.snps.csv.gz"))
# chroms_in_snps_positions = [
#     positions_file.name.split(".")[0] for positions_file in snps_positions_files
# ]
# snps_positions_data_df = pd.DataFrame(
#     {
#         "Chrom": chroms_in_snps_positions,
#         "SNPsPositionsFile": snps_positions_files,
#     }
# )

# positions_data_df = positions_data_df.merge(
#     snps_positions_data_df, on="Chrom", how="outer"
# )
# assert not positions_data_df.isna().any().any(), "There are some chroms that are missing either positions file or snps positions file. Please check the dataframes to see which ones are missing and fix the issue."


noise_positions_files = list(positions_dir.glob("*.positions.noise.csv.gz"))
chroms_in_noise_positions = [
    positions_file.name.split(".")[0] for positions_file in noise_positions_files
]
noise_positions_data_df = pd.DataFrame(
    {
        "Chrom": chroms_in_noise_positions,
        "NoisePositionsFile": noise_positions_files,
    }
)

positions_data_df = positions_data_df.merge(
    noise_positions_data_df, on="Chrom", how="outer"
)
assert not positions_data_df.isna().any().any(), "There are some chroms that are missing either positions file or noise positions file. Please check the dataframes to see which ones are missing and fix the issue."



positions_data_df

# %%
mismatches_files = list(positions_dir.glob("*.Mismatches.csv.gz"))

chroms_in_mismatches_files = [
    mismatches_file.name.split(".")[0] for mismatches_file in mismatches_files
]

mismatches_data_df = pd.DataFrame(
    {
        "Chrom": chroms_in_mismatches_files,
        "PositionsFile": mismatches_files,
    }
)

mismatches_data_df

# %%
noise_threshold_files = list(positions_dir.glob("*.NoiseThreshold.csv"))

chroms_in_noise_threshold_files = [
    noise_threshold_file.name.split(".")[0] for noise_threshold_file in noise_threshold_files
]

noise_threshold_files_df = pd.DataFrame(
    {
        "Chrom": chroms_in_noise_threshold_files,
        "NoiseThresholdFile": noise_threshold_files,
    }
)

noise_threshold_files_df

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

reads_data_df = reads_data_df.merge(unique_reads_data_df, on="Chrom", how="outer")
assert not positions_data_df.isna().any().any(), "There should be as many as reads as unique reads files."


snps_reads_files = list(reads_dir.glob("*.reads.snps.csv.gz"))
chroms_in_snps_reads_files = [reads_file.name.split(".")[0] for reads_file in snps_reads_files]
snps_reads_df = pd.DataFrame(
    {
        "Chrom": chroms_in_snps_reads_files,
        "SNPsReadsFile": snps_reads_files,
    }
)

reads_data_df = reads_data_df.merge(
    snps_reads_df, on="Chrom", how="left"
)

reads_data_df

# %%
tmr50_alignment_stats_df.loc[
    ~tmr50_alignment_stats_df["Chrom"].isin(reads_data_df["Chrom"])
]

# %%
# snps_reads_files = list(reads_dir.glob("*.reads.snps.csv.gz"))
# chroms_in_snps_reads_files = [reads_file.name.split(".")[0] for reads_file in snps_reads_files]
# snps_reads_df = pd.DataFrame(
#     {
#         "Chrom": chroms_in_snps_reads_files,
#         "SNPsReadsFile": snps_reads_files,
#     }
# )
# snps_reads_df

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

distinct_proteins_files = [
    f
    # for f in distinct_proteins_dir.glob("*DistinctUniqueProteins.*.csv")
    # for f in distinct_proteins_dir.glob("*DistinctUniqueProteins.*.Updated.csv")
    for f in distinct_proteins_dir.glob("*DistinctUniqueProteins.*.Updated.*csv")
    if "expression" not in f.name.lower()
]
chroms_in_distinct_proteins_files = [
    distinct_proteins_file.name.split(".")[0]
    for distinct_proteins_file in distinct_proteins_files
]

expression_files = list(
    expression_dir.glob("*.DistinctUniqueProteins.ExpressionLevels.*csv")
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
    orfs_df.merge(tmr50_alignment_stats_df, on="Chrom", how="right")
    .merge(positions_data_df, on="Chrom", how="left")
    .merge(reads_data_df, on="Chrom", how="left")
    .merge(proteins_data_df, on="Chrom", how="left")
)
data_df

# %%
# complete_data_df = data_df.loc[data_df["PositionsFile"].notna()].reset_index(drop=True)
# complete_data_df = data_df.loc[data_df["ReadsFile"].notna()].reset_index(drop=True)
# complete_data_df = data_df.loc[data_df["DistinctProteinsFile"].notna()].reset_index(
#     drop=True
# )

# this way we can be sure we filter `data_df` to by expression files
# in order to create a complete dataset without losing any valid distinct files
assert proteins_data_df.isna().sum().sum() == 0

complete_data_df = data_df.loc[data_df["ExpressionFile"].notna()].reset_index(drop=True)

# complete_data_df = complete_data_df.drop_duplicates(
#     "Name", keep=False, ignore_index=True
# )
complete_data_df

# %%
possibly_na_conditions = data_df["Name"].tolist()
possibly_na_chroms = data_df["Chrom"].tolist()
possibly_na_starts = data_df["Start"].tolist()
possibly_na_ends = data_df["End"].tolist()
possibly_na_strands = data_df["Strand"].tolist()

possibly_na_positions_files = data_df["PositionsFile"].tolist()
possibly_na_reads_files = data_df["ReadsFile"].tolist()
possibly_na_snps_reads_files = data_df["SNPsReadsFile"].tolist()
possibly_na_unique_reads_files = data_df["UniqueReadsFile"].tolist()
possibly_na_proteins_files = data_df["ProteinsFile"].tolist()
possibly_na_unique_proteins_files = data_df["UniqueProteinsFile"].tolist()
possibly_na_distinct_unique_proteins_files = data_df["DistinctProteinsFile"].tolist()

expression_files = complete_data_df["ExpressionFile"].tolist()

# %%
# assert (
#     data_df.loc[data_df["ExpressionFile"].notna()].reset_index(drop=True).shape
#     == data_df.loc[data_df["DistinctProteinsFile"].notna()].reset_index(drop=True).shape
# ), "some distinct proteins don't have expression levels"

if not complete_data_df["Chrom"].is_unique:
    raise ValueError("Multiple metadata rows per Chrom; resolve file versions first")

required_for_counts = ["Chrom", "Name", "DistinctProteinsFile", "UniqueReadsFile"]
if complete_data_df[required_for_counts].isna().any().any():
    raise ValueError("Missing metadata or required input for an existing distinct result")

# %%
# complete_data_df[["Chrom"]].to_csv("TMR50.CompleteData.Chroms.tsv", sep="\t", index=False)

# %%
# complete_data_df.loc[complete_data_df["Chrom"] == robo2_chrom]

# %%
# robo2_index = complete_data_df.loc[complete_data_df["Chrom"] == robo2_chrom].index[0]
# robo2_index

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
snps_reads_files = complete_data_df["SNPsReadsFile"].tolist()
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
100 * len(chroms) / len(possibly_na_chroms)

# %% [markdown]
# # Data loading - TMR 1000

# %%
tmr1000_main_data_dir = Path(
    "/private6/projects/Combinatorics/O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage1000.PooledSamples"
)

tmr1000_positions_dir = Path(tmr1000_main_data_dir, "PositionsFiles")
tmr1000_reads_dir = Path(tmr1000_main_data_dir, "ReadsFiles")
tmr1000_proteins_dir = Path(tmr1000_main_data_dir, "ProteinsFiles")
tmr1000_distinct_proteins_dir = Path(tmr1000_main_data_dir, "DistinctProteins")

# %%
tmr1000_alignment_stats_df = alignment_stats_df.loc[
    alignment_stats_df["MappedReads"] >= 1000
]
tmr1000_alignment_stats_df

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
tmr1000_mismatches_files = list(tmr1000_positions_dir.glob("*.Mismatches.csv.gz"))

tmr1000_chroms_in_mismatches_files = [
    mismatches_file.name.split(".")[0] for mismatches_file in tmr1000_mismatches_files
]

tmr1000_mismatches_data_df = pd.DataFrame(
    {
        "Chrom": tmr1000_chroms_in_mismatches_files,
        "PositionsFile": tmr1000_mismatches_files,
    }
)

tmr1000_mismatches_data_df

# %%
tmr1000_noise_threshold_files = list(tmr1000_positions_dir.glob("*.NoiseThreshold.csv"))

tmr1000_chroms_in_noise_threshold_files = [
    noise_threshold_file.name.split(".")[0] for noise_threshold_file in tmr1000_noise_threshold_files
]

tmr1000_noise_threshold_files_df = pd.DataFrame(
    {
        "Chrom": tmr1000_chroms_in_noise_threshold_files,
        "NoiseThresholdFile": tmr1000_noise_threshold_files,
    }
)

tmr1000_noise_threshold_files_df

# %%
tmr1000_alignment_stats_df.loc[
    ~tmr1000_alignment_stats_df["Chrom"].isin(tmr1000_positions_data_df["Chrom"])
]

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
tmr1000_alignment_stats_df.loc[
    ~tmr1000_alignment_stats_df["Chrom"].isin(tmr1000_reads_data_df["Chrom"])
]

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

# TODO uncomment later
# expression_data_df = pd.DataFrame(
#     {"Chrom": chroms_in_expression_files, "ExpressionFile": expression_files}
# )

tmr1000_proteins_data_df = (
    tmr1000_proteins_data_df
    .merge(
        tmr1000_unique_proteins_data_df, on="Chrom", how="left"
    )
    .merge(tmr1000_distinct_proteins_data_df, on="Chrom", how="left") # TODO uncomment later
    # .merge(expression_data_df, on="Chrom", how="left")
)

tmr1000_proteins_data_df

# %%
tmr1000_data_df = (
    orfs_df.merge(tmr1000_alignment_stats_df, on="Chrom", how="right")
    .merge(tmr1000_positions_data_df, on="Chrom", how="left")
    .merge(tmr1000_reads_data_df, on="Chrom", how="left")
    .merge(tmr1000_proteins_data_df, on="Chrom", how="left")
)

tmr1000_data_df

# %%
# tmr1000_complete_data_df = tmr1000_data_df.loc[
#     tmr1000_data_df["ReadsFile"].notna()
# ].reset_index(drop=True)
tmr1000_complete_data_df = tmr1000_data_df.loc[
    tmr1000_data_df["DistinctProteinsFile"].notna()
].reset_index(drop=True)
# complete_data_df = data_df.loc[data_df["ExpressionFile"].notna()].reset_index(drop=True)

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
tmr1000_distinct_unique_proteins_files = tmr1000_complete_data_df[
    "DistinctProteinsFile"
].tolist()
# expression_files = complete_data_df["ExpressionFile"].tolist()


if not tmr1000_complete_data_df["Chrom"].is_unique:
    raise ValueError("TMR1000 metadata contain duplicate Chrom values")

for chrom, distinct_file, unique_reads_file in tmr1000_complete_data_df[
    ["Chrom", "DistinctProteinsFile", "UniqueReadsFile"]
].itertuples(index=False, name=None):
    for path in (distinct_file, unique_reads_file):
        if pd.isna(path) or Path(path).name.split(".")[0] != chrom:
            raise ValueError(f"File/Chrom mismatch: {chrom!r}, {path!r}")


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

# %%
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
pio.templates.default = "plotly_white"
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

# %%
mismatches = sorted(
    list(
        {
            f"{ref_base}>{alt_base}"
            for ref_base in "ATCG"
            for alt_base in "ATCG"
            if ref_base != alt_base
        }
    )
)
# mismatches
mismatches_color_sequence = px.colors.qualitative.Dark24
mismatch_dolor_map = {
    mismatch: color
    for mismatch, color in zip(
        mismatches, mismatches_color_sequence
    )
}
# mismatch_dolor_map

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
    # .count()
    .sum() # fixed bug on 7.9.2026
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

# %% [markdown]
# ### Base positions

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
concat_all_edited_positions_df = concat_all_positions_df.loc[
    concat_all_positions_df["EditedFinal"]
]
concat_all_edited_positions_df

# %%
concat_all_positions_df.columns

# %%
concat_all_positions_df["NoisyFinal"].value_counts(dropna=False)

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

# %%
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
# ### SNPs [to be deleted]

# %%
# positions_data_df

# %%
# snps_positions_files_of_edited_genes = positions_data_df.loc[
#     positions_data_df["Chrom"].isin(chroms),
#     "SNPsPositionsFile"
# ].values
# snps_positions_files_of_edited_genes

# %%
# concat_snps_positions_df = pd.concat(
#     [
#         pd.read_table(
#             snps_positions_file, sep="\t", dtype={"Reads": str}
#         ).drop(
#             columns=[
#                 "EditingBinomPVal", "EditingCorrectedPVal", "EditedCorrected", "BelowEditingFreq1"
#             ]
#         )
#         for snps_positions_file in positions_data_df["SNPsPositionsFile"].values
#     ]
# )
# concat_snps_positions_df

# %%
# concat_snps_positions_df["Chrom"].nunique()

# %%
# concat_snps_positions_df["Chrom"].value_counts().describe().round(2)

# %% [markdown]
# ### Mismatches

# %%
def find_alt_base(
    ref_base: str, a_count: int, t_count: int, c_count: int, g_count: int,
    edited: bool
):
    """
    Find the base with most supporting reads other than `ref_base`.
    If there are two or more such bases, the function picks one at random.
    If position is `edited` (and thus its `ref_base` is `A`), its alt base is G by definition.
    """
    if edited:
        if ref_base != "A":
            raise ValueError(f"Edited position should have ref_base 'A', but got {ref_base}")
        return "G"
        
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
noise_threshold_df = pd.concat(
    [
        pd.read_csv(
            noise_threshold_file, 
            sep="\t",
            names=["Chrom", "NoiseThreshold"]
        )
        for noise_threshold_file in noise_threshold_files_df["NoiseThresholdFile"].values
    ],
    ignore_index=True
)
noise_threshold_df

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
            x["EditedFinal"]
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

mismatches_df = mismatches_df.merge(
    noise_threshold_df,
    on="Chrom",
    how="left"
)

mismatches_df["SNP"] = (
    (mismatches_df["NoisyFinal"])
    & (mismatches_df["Noise"] >= snp_noise_level)
)

mismatches_df

# %%
# # so it seems that the same data about snps can be extracted from the general positions files and from the 
# # snps-specific positions files - which is reassuring. the mismatch annotation also allows to get the mismatch 
# # type and frequency, which can be useful for downstream analyses and for setting a noise threshold based on the 
# # distribution of mismatch frequencies at known snps.

# df = concat_snps_positions_df.loc[
#     :,
#     ["Chrom", "Position", "Noise", "AboveEditingThreshold"]
# ].merge(
#     mismatches_df.loc[
#         (mismatches_df["NoisyFinal"])
#         & (mismatches_df["MismatchFrequency"].ge(snp_noise_level)),
#         ["Chrom", "Position", "MismatchFrequency", "NoiseThreshold"]
#     ],
#     how="outer"
# )

# assert df.apply(
#     lambda x: np.isclose(x["MismatchFrequency"], x["Noise"]),
#     axis=1
# ).all()

# assert df.loc[
#     (df["AboveEditingThreshold"])
#     & (df["Noise"].lt(df["NoiseThreshold"])),
# ].empty

# del df

# # following that, i can use whatever is convinient for me - 
# # but it's still good i made ahead the snps reads files

# %%
tweleve_mismatches_to_save_df = (
    mismatches_df
    .loc[
        :,
        [
            'Transcript', 'Chrom', 'Position', 'Mismatch', 'TotalCoverage',
            'MismatchFrequency', 'NoiseThreshold', 'NoisyFinal', 'EditedFinal',
            "SNP"
        ]
    ]
    .rename(
        columns={
            "NoiseThreshold": "EditingThreshold", 
            "NoisyFinal": "Noisy", 
            "EditedFinal": "Edited",
            "Transcript": "Gene"
        }
    )
    .assign(
        AboveEditingThreshold=lambda x: x["MismatchFrequency"] >= x["EditingThreshold"]
    )
)

# only consider octopus genes in which we detected editing
edited_octopus_chroms = tweleve_mismatches_to_save_df.loc[
    tweleve_mismatches_to_save_df["Edited"],
    "Chrom"
].unique()
ic(len(edited_octopus_chroms))
tweleve_mismatches_to_save_df = tweleve_mismatches_to_save_df.loc[
    tweleve_mismatches_to_save_df["Chrom"].isin(edited_octopus_chroms)
].reset_index(drop=True)

tweleve_mismatches_to_save_df

# %%
# noise_threshold_df.loc[
#     (noise_threshold_df["NoiseThreshold"].eq(0))
#     & (noise_threshold_df["Chrom"].isin(edited_octopus_chroms))
# ]

# %%
# tweleve_mismatches_to_save_df.loc[
#     (tweleve_mismatches_to_save_df["Mismatch"].eq("A>G"))
#     & (tweleve_mismatches_to_save_df["MismatchFrequency"].ge(tweleve_mismatches_to_save_df["EditingThreshold"]))
#     & (tweleve_mismatches_to_save_df["MismatchFrequency"].gt(0))
#     # & (~tweleve_mismatches_to_save_df["Edited"])
# ]

# %%
# tweleve_mismatches_to_save_df.loc[
#     (tweleve_mismatches_to_save_df["Mismatch"].eq("A>G"))
#     & (tweleve_mismatches_to_save_df["MismatchFrequency"].ge(tweleve_mismatches_to_save_df["EditingThreshold"]))
#     & (tweleve_mismatches_to_save_df["MismatchFrequency"].gt(0))
#     # & (~tweleve_mismatches_to_save_df["Edited"])
# ].groupby("Chrom")["Edited"].sum().describe()

# %%
tweleve_mismatches_to_save_df.to_csv(
    Path(out_dir, "12Mismatches.Octopus.PacBio.csv"),
    sep="\t",
    index=False
)

# %%

# %% [markdown]
# ### SNPs and coverage per gene

# %%
snps_and_coverage_per_gene_df = mismatches_df.loc[
    (mismatches_df["NoisyFinal"])
    & (mismatches_df["MismatchFrequency"].ge(snp_noise_level))
].groupby("Chrom").size().reset_index(name="SNPs").merge(
    alignment_stats_df.loc[:, ["Chrom", "MappedReads"]],
    how="outer"
).fillna(0)
snps_and_coverage_per_gene_df["SNPs"] = (
    snps_and_coverage_per_gene_df["SNPs"].astype(int)
)
snps_and_coverage_per_gene_df

# %%
tmr50_alignment_stats_df

# %% [markdown]
# The next two dfs allow to easily filter for genes which were intially considered as "editable" (i.e. genes with a low number of SNPs).

# %%
tmr50_alignment_stats_and_snps_df = tmr50_alignment_stats_df.merge(
    snps_and_coverage_per_gene_df.loc[
        snps_and_coverage_per_gene_df["SNPs"].le(max_snps_per_gene_to_allow_editing_detection)
    ],
    how="inner"
)
tmr50_alignment_stats_and_snps_df

# %%
tmr1000_alignment_stats_and_snps_df = tmr1000_alignment_stats_df.merge(
    snps_and_coverage_per_gene_df.loc[
        snps_and_coverage_per_gene_df["SNPs"].le(max_snps_per_gene_to_allow_editing_detection)
    ],
    how="inner"
)
tmr1000_alignment_stats_and_snps_df


# %% [markdown]
# ### Mismatches - theoritical fix

# %% [markdown]
# How the mismatches table would've looked if we positions with A in the reference wouldn't automatically be considered as
# A-to-G but rather A-to-X based on the highest covered base other than A.

# %%
def find_alt_base_theoretical_fix(
    ref_base: str, 
    a_count: int, t_count: int, c_count: int, g_count: int,
):
    """
    Find the base with most supporting reads other than `ref_base`.
    If there are two or more such bases, the function picks one at random.
    If all base counts other than `ref_base` are zero, the function returns NaN - signaling that there is no alternative base at this position.
    """
        
    base_counts_dict = {"A": a_count, "T": t_count, "C": c_count, "G": g_count}
    alt_bases = set(base_counts_dict) - {ref_base}
    alt_base_counts_dict = {
        base: base_count
        for base, base_count in base_counts_dict.items()
        if base in alt_bases
    }
    max_alt_base_count = max(alt_base_counts_dict.values())
    if max_alt_base_count == 0:
        return np.nan
    max_alt_bases = [
        base
        for base, base_count in alt_base_counts_dict.items()
        if base_count == max_alt_base_count
    ]
    alt_base = choice(max_alt_bases)
    return alt_base


# %%
def define_mismatch_type_theoretical_fix(
    ref_base: str, alt_base: str, 
    # strand: str, 
    ):
    # this is not a real mismatch, so we cannot define a mismatch type
    if pd.isna(alt_base):
        return np.nan
    
    # if strand == "-":
    #     opposing_bases = {"A": "T", "T": "A", "C": "G", "G": "C"}
    #     ref_base = opposing_bases[ref_base]
    #     alt_base = opposing_bases[alt_base]
    
    mismatch_type = f"{ref_base}>{alt_base}"
    
    return mismatch_type


# %%
def mismatch_frequency_theoretical_fix(ref_base_count, alt_base_count):
    if alt_base_count == 0:
        return np.nan
    return alt_base_count / (ref_base_count + alt_base_count)


# %%
def make_rows_of_all_positions_in_orf(
    transcript: str,
    chrom: str, 
    start: int, 
    end: int
) -> pd.DataFrame:
    # inclusive_end = end - 1
    df = pd.DataFrame(
        {
            "Transcript": [transcript] * (end - start),
            "Chrom": [chrom] * (end - start),
            "Position": list(range(start, end)),
        }
    )
    assert df.shape[0] % 3 == 0, f"ORF length is not a multiple of 3: {chrom}:{start}-{end}"
    return df


# %%
def define_editing_threshold(
    one_chrom_positions_df: pd.DataFrame,
    top_x_noisy_positions: int = 3,
    assurance_factor: float = 1.5
):
    noise_levels = (
        one_chrom_positions_df
        .loc[
            one_chrom_positions_df["NoiseSite"],
            "MismatchFrequency"
        ]
        .sort_values(ascending=False)
        [:top_x_noisy_positions]
        .tolist()
    )
    # if there are less noisy positions than `top_x_noisy_positions`, add zeros accordingly
    noise_levels = pd.Series(
        noise_levels + [0 for _ in range(top_x_noisy_positions - len(noise_levels))]
    )
    editing_threshold = noise_levels.mean()
    if pd.isna(editing_threshold):
        editing_threshold = 0
    # anyway, we finalize the editing threshold
    editing_threshold *= assurance_factor
    return editing_threshold


# %%
noise_threshold_df = pd.concat(
    [
        pd.read_csv(
            noise_threshold_file, 
            sep="\t",
            names=["Chrom", "NoiseThreshold"]
        )
        for noise_threshold_file in noise_threshold_files_df["NoiseThresholdFile"].values
    ],
    ignore_index=True
)
noise_threshold_df

# %%
# load the positions data as found by our pipeline
mismatches_theoretical_fix_df = concat_all_positions_df.loc[
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
        # "EditingFrequency",
        # "Edited",
        # "EditedCorrected",
        "EditedFinal",
        # "Noise",
        # "NoisyCorrected",
        "NoisyFinal",
    ],
]

# mismatches_theoretical_fix_df = mismatches_theoretical_fix_df.merge(data_df.loc[:, ["Chrom", "Strand"]], how="left")
# mismatches_theoretical_fix_df.insert(
#     mismatches_theoretical_fix_df.columns.get_loc("Chrom") + 1, "Strand2", mismatches_theoretical_fix_df["Strand"]
# )
# del mismatches_theoretical_fix_df["Strand"]
# mismatches_theoretical_fix_df = mismatches_theoretical_fix_df.rename(columns={"Strand2": "Strand"})

# all orfs are on the positive strand so we don't have to deal with positive strand/coding strand normalizations
assert set(orfs_df["Strand"].values) == {"+"}

# find the "real" alternative base for each position
mismatches_theoretical_fix_df.insert(
    mismatches_theoretical_fix_df.columns.get_loc("RefBase") + 1,
    "AltBase",
    mismatches_theoretical_fix_df.apply(
        lambda x: find_alt_base_theoretical_fix(
            x["RefBase"],
            x["A"],
            x["T"],
            x["C"],
            x["G"],
        ),
        axis=1,
    ),
)
# annotathe mismatch type according to alt- and ref-base
mismatches_theoretical_fix_df.insert(
    mismatches_theoretical_fix_df.columns.get_loc("AltBase") + 1,
    "Mismatch",
    mismatches_theoretical_fix_df.apply(
        lambda x: define_mismatch_type_theoretical_fix(x["RefBase"], x["AltBase"]), axis=1
    ),
)

# annotate Alt- and RefBaseCount as dedicated cols
mismatches_theoretical_fix_df.insert(
    mismatches_theoretical_fix_df.columns.get_loc("G") + 1,
    "RefBaseCount",
    mismatches_theoretical_fix_df.apply(
        lambda x: x[x["RefBase"]],
        axis=1
    )
)
mismatches_theoretical_fix_df.insert(
    mismatches_theoretical_fix_df.columns.get_loc("G") + 2,
    "AltBaseCount",
    mismatches_theoretical_fix_df.apply(
        lambda x: x[x["AltBase"]] if pd.notna(x["AltBase"]) else np.nan,
        axis=1
    )
)

# annotate mismatch frequency as an alternative and general annotation to `EditingFrequency` and `Noise`
mismatches_theoretical_fix_df.insert(
    mismatches_theoretical_fix_df.columns.get_loc("AltBaseCount") + 1,
    "MismatchFrequency",
    mismatches_theoretical_fix_df.apply(
        lambda x: mismatch_frequency(x["RefBaseCount"], x["AltBaseCount"]), axis=1
    ),
)

# add to mismatches_theoretical_fix_df all positions in the ORFs, even if they are not covered by reads and thus 
# not present in the positions files
cocnat_all_cds_positions_df = pd.concat(
    [
        make_rows_of_all_positions_in_orf(transcript, chrom, start, end)
        for transcript, chrom, start, end in orfs_df.loc[:, ["Name", "Chrom", "Start", "End"]].itertuples(index=False, name=None)
    ]
)
mismatches_theoretical_fix_df_2 = mismatches_theoretical_fix_df.merge(
    cocnat_all_cds_positions_df, 
    on=["Chrom", "Transcript", "Position"], 
    how="right"
)
assert mismatches_theoretical_fix_df_2.shape[0] == cocnat_all_cds_positions_df.shape[0]
assert mismatches_theoretical_fix_df_2.loc[
    mismatches_theoretical_fix_df_2["RefBase"].notna()
].shape[0] == mismatches_theoretical_fix_df.shape[0]
mismatches_theoretical_fix_df = mismatches_theoretical_fix_df_2
# fill coverage of positions that are not covered in the original positions files with 0
uncovered_positions = mismatches_theoretical_fix_df["RefBase"].isna()
mismatches_theoretical_fix_df.loc[
    uncovered_positions,
    "TotalCoverage"
] = 0

# add original noise threshold per chrom
mismatches_theoretical_fix_df = mismatches_theoretical_fix_df.merge(
    noise_threshold_df,
    on="Chrom",
    how="left"
)

# keep only chroms that have at least one position with coverage > 0, as otherwise the "original" noise threshold is meaningless
# mismatches_theoretical_fix_df = mismatches_theoretical_fix_df.loc[
#     mismatches_theoretical_fix_df.groupby("Chrom")["TotalCoverage"].transform("sum").gt(0)
# ]
# keep only chroms that have at least 50 mapped reads
mismatches_theoretical_fix_df = mismatches_theoretical_fix_df.loc[
    mismatches_theoretical_fix_df["Chrom"].isin(
        tmr50_alignment_stats_df["Chrom"]
    )
]
# verify that all chroms in the resulting df indeed linked to the "original" noise threshold
assert mismatches_theoretical_fix_df["NoiseThreshold"].isna().sum() == 0

# intialy set all pvals to 1, which represnts the pvalue of the perfect null hypothesis - no mismatch in this position
mismatches_theoretical_fix_df["BinomPVal"] = 1.0
# then, for positions that have an alternative base, calculate the binomial test p-value for the observed number of 
# alternative base reads given the total coverage and a null hypothesis probability of 0.001
positions_with_alt_base = mismatches_theoretical_fix_df["AltBase"].notna()
mismatches_theoretical_fix_df.loc[
    positions_with_alt_base,
    "BinomPVal"
] = mismatches_theoretical_fix_df.loc[positions_with_alt_base].apply(
    lambda x: binom_test(
        x["AltBaseCount"],
        x["RefBaseCount"] + x["AltBaseCount"],
        0.001,
        alternative="larger",
    ),
    axis=1,
)
# perform the BH correction
bh_rejections, bh_corrected_pvals = fdrcorrection(
    mismatches_theoretical_fix_df["BinomPVal"]
)
# add the corrected p-values and rejections to the concatenated dataframe
mismatches_theoretical_fix_df["BHCorrectedPVal"] = bh_corrected_pvals
mismatches_theoretical_fix_df["BHRejection"] = bh_rejections

# %%
mismatches_theoretical_fix_df

# %%
mismatches_theoretical_fix_df["Chrom"].nunique()

# %%
significant_mismatches_theoretical_fix_df = mismatches_theoretical_fix_df.loc[
    mismatches_theoretical_fix_df["BHRejection"]
].reset_index(drop=True)

# SNPs (definitive & suspected SNPs)

significant_mismatches_theoretical_fix_df["AtOrAboveSuspectedSNPLevel"] = (
    significant_mismatches_theoretical_fix_df["MismatchFrequency"].ge(snp_noise_level)
)
significant_mismatches_theoretical_fix_df["MismatchFrequency1"] = (
    significant_mismatches_theoretical_fix_df["MismatchFrequency"].eq(1)
)

significant_mismatches_theoretical_fix_df["SuspectedSNP"] = (
    (
        significant_mismatches_theoretical_fix_df["AtOrAboveSuspectedSNPLevel"]
        & ~significant_mismatches_theoretical_fix_df["MismatchFrequency1"]
        & significant_mismatches_theoretical_fix_df["Mismatch"].ne("A>G")
    )
)
significant_mismatches_theoretical_fix_df["DefinitiveSNP"] = (
    significant_mismatches_theoretical_fix_df["MismatchFrequency1"]
)
significant_mismatches_theoretical_fix_df["SNP"] = (
    significant_mismatches_theoretical_fix_df["SuspectedSNP"]
    | significant_mismatches_theoretical_fix_df["DefinitiveSNP"]
)

significant_mismatches_theoretical_fix_df["NumOfSuspectedSNPsPerChrom"] = (
    significant_mismatches_theoretical_fix_df.groupby("Chrom")["SuspectedSNP"].transform("sum")
)
significant_mismatches_theoretical_fix_df["NumOfDefinitiveSNPsPerChrom"] = (
    significant_mismatches_theoretical_fix_df.groupby("Chrom")["DefinitiveSNP"].transform("sum")
)
significant_mismatches_theoretical_fix_df["NumOfSNPsPerChrom"] = (
    significant_mismatches_theoretical_fix_df.groupby("Chrom")["SNP"].transform("sum")
)

# don't allow editing detection in genes where the number of suspected SNPs is above a certain threshold, 
# as this may indicate that the gene is highly polymorphic and thus not suitable for editing detection
# as it's hard to distinguish between editing and suspected SNPs
significant_mismatches_theoretical_fix_df["EditingDetectionDisabledInChromDueToSuspectedSNPs"] = (
    significant_mismatches_theoretical_fix_df["NumOfSuspectedSNPsPerChrom"].gt(
        max_snps_per_gene_to_allow_editing_detection
    )
)

significant_mismatches_theoretical_fix_df["NoiseSite"] = (
    significant_mismatches_theoretical_fix_df["Mismatch"].ne("A>G")
    & ~significant_mismatches_theoretical_fix_df["SNP"]
)

new_editing_thresholds_df = (
    significant_mismatches_theoretical_fix_df
    .groupby("Chrom")
    .apply(define_editing_threshold)
    .reset_index(name="NewEditingThreshold")
)
significant_mismatches_theoretical_fix_df = significant_mismatches_theoretical_fix_df.merge(
    new_editing_thresholds_df,
    on="Chrom",
    how="left"
)

significant_mismatches_theoretical_fix_df["AboveNewEditingThreshold"] = (
    significant_mismatches_theoretical_fix_df["MismatchFrequency"].gt(
        significant_mismatches_theoretical_fix_df["NewEditingThreshold"]
    )
)

significant_mismatches_theoretical_fix_df["EditingSite"] = (
    significant_mismatches_theoretical_fix_df["Mismatch"].eq("A>G")
    & ~significant_mismatches_theoretical_fix_df["SNP"]
    & ~significant_mismatches_theoretical_fix_df["EditingDetectionDisabledInChromDueToSuspectedSNPs"]
    & significant_mismatches_theoretical_fix_df["AboveNewEditingThreshold"]
)

significant_mismatches_theoretical_fix_df

# %%
significant_mismatches_theoretical_fix_df["Chrom"].nunique()


# %%
def is_x_gt_or_eq_or_lt_y(x, y):
    if pd.isna(x) or pd.isna(y):
        return pd.NA
    elif x > y:
        return ">"
    elif x == y:
        return "=="
    else:
        return "<"


# %%
(
    significant_mismatches_theoretical_fix_df
    .drop_duplicates("Chrom")
    .reset_index(drop=True)
    .loc[
        :,
        [
            # "Chrom", 
         "NoiseThreshold", "NewEditingThreshold"]
    ]
    .apply(
        lambda x: is_x_gt_or_eq_or_lt_y(x["NoiseThreshold"], x["NewEditingThreshold"]),
        axis=1
    )
    .value_counts()
    .reset_index()
    .rename(
        columns={
            "index": "NoiseThreshold vs NewEditingThreshold", 
        }
    )
)

# %%
(
    significant_mismatches_theoretical_fix_df
    .groupby("Chrom")
    [["SuspectedSNP", "DefinitiveSNP", "SNP"]]
    .sum().describe().round(2)
)

# %%
(
    significant_mismatches_theoretical_fix_df
    ["Mismatch"].value_counts(dropna=False)
)

# %%
significant_edited_positions_theoretical_fix_df = significant_mismatches_theoretical_fix_df.loc[
    significant_mismatches_theoretical_fix_df["EditingSite"]
]
significant_edited_positions_theoretical_fix_df

# %%
concat_all_edited_positions_df

# %%
significant_edited_positions_theoretical_fix_df.columns

# %%
significant_mismatches_theoretical_fix_df.loc[
    significant_mismatches_theoretical_fix_df["NumOfSuspectedSNPsPerChrom"].gt(max_snps_per_gene_to_allow_editing_detection)
]

# %%
significant_mismatches_theoretical_fix_df["NumOfSuspectedSNPsPerChrom"].describe().round(2)

# %%
significant_edited_positions_theoretical_fix_df["NumOfSuspectedSNPsPerChrom"].describe().round(2)

# %%
concat_old_vs_new_editing_positions_df = (
    significant_edited_positions_theoretical_fix_df
    .loc[
        :,
        [
            "Chrom", "Position", "MismatchFrequency", "NoiseThreshold", "NewEditingThreshold",
            "NumOfSuspectedSNPsPerChrom"
         ]
    ]
    .rename(
        columns={
            "MismatchFrequency": "NewMismatchFrequency",
            "NoiseThreshold": "OldEditingThreshold",
            "NumOfSuspectedSNPsPerChrom": "NumOfNewSuspectedSNPsPerChrom",
        }
    )
    .merge(
        (
            concat_all_edited_positions_df
            .loc[
                :,
                ["Chrom", "Position", "EditingFrequency"]
            ]
            .rename(
                columns={
                    "EditingFrequency": "OldMismatchFrequency",
                }
            )
        ),
        on=["Chrom", "Position"],
        how="outer",
        suffixes=("_New", "_Old"),
        indicator=True
    )
    .loc[
        :,
        [
            "Chrom", "Position", "OldMismatchFrequency", "NewMismatchFrequency",
            "OldEditingThreshold", "NewEditingThreshold", "NumOfNewSuspectedSNPsPerChrom", "_merge"
        ]
    ]
)
concat_old_vs_new_editing_positions_df["_merge"] = (
    concat_old_vs_new_editing_positions_df["_merge"]
    .astype(str)
    .replace(
        {
            "left_only": "New",
            "right_only": "Old",
            "both": "Both"
        }
    )
)
concat_old_vs_new_editing_positions_df

# %%
concat_old_vs_new_editing_positions_df["_merge"].value_counts()

# %%
concat_old_vs_new_editing_positions_only_new_df = (
    concat_old_vs_new_editing_positions_df
    .loc[
        concat_old_vs_new_editing_positions_df["_merge"].eq("New")
    ]
)
concat_old_vs_new_editing_positions_only_old_df = (
    concat_old_vs_new_editing_positions_df
    .loc[
        concat_old_vs_new_editing_positions_df["_merge"].eq("Old")
    ]
)


# %%
concat_old_vs_new_editing_positions_only_new_df

# %%
concat_old_vs_new_editing_positions_only_old_df

# %%
concat_old_vs_new_editing_positions_only_old_df["Chrom"].nunique()

# %%
sites_sets_per_chrom = (
    concat_old_vs_new_editing_positions_df.groupby("Chrom")["_merge"].value_counts()
    .reset_index()
    .pivot(
        index="Chrom", columns="_merge", values="count"
    )
    .fillna(0)
    .astype(int)
)
sites_sets_per_chrom

# %%
# sites_sets_per_chrom.loc[
#     (
#         sites_sets_per_chrom["Old"].gt(0)
#     )
# ]

# %%
chroms_with_only_old_sites = sites_sets_per_chrom.loc[
    (
        sites_sets_per_chrom["Both"].eq(0)
        & sites_sets_per_chrom["Old"].gt(0)
        & sites_sets_per_chrom["New"].eq(0) 
    )
].index.tolist()

ic(len(chroms_with_only_old_sites))

chroms_with_only_old_sites

# %%
sites_sets_per_chrom.loc[
    chroms_with_only_old_sites,
    "Old"
].sum()

# %%
# concat_old_vs_new_editing_positions_only_old_df.loc[
#     concat_old_vs_new_editing_positions_only_old_df["Chrom"].isin(
#         chroms_with_only_old_sites
#     )
# ]

# %%
# (
#     significant_mismatches_theoretical_fix_df
#     .loc[
#         significant_mismatches_theoretical_fix_df["Chrom"].isin(
#             chroms_with_only_old_sites
#         ),
#         ["Chrom", "NumOfSuspectedSNPsPerChrom"]
#     ]
#     .drop_duplicates("Chrom")
# )

# %%
significant_mismatches_theoretical_fix_chroms_with_too_much_suspected_snps_df = (
    significant_mismatches_theoretical_fix_df
    .loc[
        significant_mismatches_theoretical_fix_df["Chrom"].isin(
            chroms_with_only_old_sites
        ),
        # ["Chrom", "NumOfSuspectedSNPsPerChrom"]
    ]
    # .drop_duplicates("Chrom")
)
significant_mismatches_theoretical_fix_chroms_with_too_much_suspected_snps_df

# %%
significant_mismatches_theoretical_fix_chroms_with_too_much_suspected_snps_df["Chrom"].nunique()

# %%
significant_mismatches_theoretical_fix_df.loc[
    (
        significant_mismatches_theoretical_fix_df["Mismatch"].eq("A>G")
        & ~significant_mismatches_theoretical_fix_df["SNP"]
        & ~significant_mismatches_theoretical_fix_df["EditingDetectionDisabledInChromDueToSuspectedSNPs"]
        & significant_mismatches_theoretical_fix_df["AboveNewEditingThreshold"]
    )
]


# %%
# this table contains A>G sites that are not considered edited only because
# the gene they are in has too many suspected SNPs, which may indicate that the gene is highly polymorphic and thus not suitable for editing detection
significant_mismatches_sites_disabled_by_suspected_snps_theoretical_fix_df = significant_mismatches_theoretical_fix_df.loc[
    (
        significant_mismatches_theoretical_fix_df["Mismatch"].eq("A>G")
        & ~significant_mismatches_theoretical_fix_df["SNP"]
        & significant_mismatches_theoretical_fix_df["EditingDetectionDisabledInChromDueToSuspectedSNPs"]
        & significant_mismatches_theoretical_fix_df["AboveNewEditingThreshold"]
    )
]
significant_mismatches_sites_disabled_by_suspected_snps_theoretical_fix_df

# %%
(
    significant_mismatches_sites_disabled_by_suspected_snps_theoretical_fix_df
    ["Chrom"]
    .nunique()
)

# %%
(
    significant_mismatches_sites_disabled_by_suspected_snps_theoretical_fix_df
    .drop_duplicates(["Chrom"])
    ["NumOfSuspectedSNPsPerChrom"].describe().round(2)
)

# %%
(
    significant_mismatches_sites_disabled_by_suspected_snps_theoretical_fix_df
    .groupby("NumOfSuspectedSNPsPerChrom")
    .size()
    .reset_index(name="NumOfRejectedEditingSites")
    .assign(
        CumulativeRejectedEditingSites=lambda x: x["NumOfRejectedEditingSites"].cumsum()
    )
    # .describe().round(2)
)

# %%
fig = px.scatter(
    (
        significant_mismatches_sites_disabled_by_suspected_snps_theoretical_fix_df
        .groupby("NumOfSuspectedSNPsPerChrom")
        .size()
        .reset_index(name="NumOfRejectedEditingSites")
        .assign(
            CumulativeRejectedEditingSites=lambda x: x["NumOfRejectedEditingSites"].cumsum()
        )
        # .describe().round(2)
    ),
    x="NumOfSuspectedSNPsPerChrom",
    # y="NumOfRejectedEditingSites",
    y="CumulativeRejectedEditingSites",
    log_x=True,
    labels={
        "NumOfSuspectedSNPsPerChrom": "Number of suspected SNPs per gene",
        # "NumOfRejectedEditingSites": "Number of editing sites rejected<br>due to suspected SNPs",
        "CumulativeRejectedEditingSites": "Cumulative number of editing sites<br>rejected due to suspected SNPs",
    }
)
fig.update_yaxes(tick0=0, rangemode="tozero", dtick=5000)
fig.update_layout(
    width=600,
    height=400
)
fig.show()

# %%
(
    significant_mismatches_sites_disabled_by_suspected_snps_theoretical_fix_df
    .groupby("NumOfSuspectedSNPsPerChrom")
    ["EditedFinal"]
    .apply(
        lambda x: x.sum()
    )
    .reset_index(name="NumOfRejectedEditingSitesPreviouslyFound")
    .assign(
        CumulativeNumOfRejectedEditingSitesPreviouslyFound=lambda x: x["NumOfRejectedEditingSitesPreviouslyFound"].cumsum()
    )
)

# %%
# # how many SNPs in 4-6 SPNs per gene are new AC/AT?
# (
#     significant_mismatches_theoretical_fix_df
#     .loc[
#         (
#             significant_mismatches_theoretical_fix_df["NumOfSuspectedSNPsPerChrom"].between(4, 6)
#             & significant_mismatches_theoretical_fix_df["SuspectedSNP"]
#         )
#     ]
#     .groupby("NumOfSuspectedSNPsPerChrom")
#     ["Mismatch"]
#     # .value_counts(normalize=True)
#     # .mul(100)
#     # .round(2)
#     .value_counts()
#     .unstack(fill_value=0)
# )

# %%
significant_mismatches_theoretical_fix_df

# %%
# sites_sets_per_chrom.loc[chroms_with_only_old_sites, ["Old"]].reset_index()

# %%
# how many SNPs in 4-6 SPNs per gene are new AC/AT? 
new_ac_at_snps_at_genes_with_only_old_sites_df = (
    significant_mismatches_theoretical_fix_df
    .loc[
        (
            significant_mismatches_theoretical_fix_df["NumOfSuspectedSNPsPerChrom"].between(4, 6)
            & significant_mismatches_theoretical_fix_df["SuspectedSNP"]
            & significant_mismatches_theoretical_fix_df["Chrom"].isin(chroms_with_only_old_sites)
            & significant_mismatches_theoretical_fix_df["Mismatch"].isin(["A>C", "A>T"])
        ),
    ]
    .groupby(["Chrom", "NumOfSuspectedSNPsPerChrom"])
    .size()
    .reset_index(name="NumOfACOrATSuspectedSNPsPerChrom")
    .assign(
        NumOfSuspectedSNPsPerChromWoACOrAT = lambda x: x["NumOfSuspectedSNPsPerChrom"] - x["NumOfACOrATSuspectedSNPsPerChrom"]
    )
    .assign(
        EditingSitesRescuedByRemovingACOrATSuspectedSNPs = lambda x: (
            x["NumOfSuspectedSNPsPerChromWoACOrAT"] <= max_snps_per_gene_to_allow_editing_detection
        )
    )
    .merge(
        (
            sites_sets_per_chrom.loc[chroms_with_only_old_sites, ["Old"]]
            .reset_index()
            .rename(columns={"Old": "NumOfLostEditingSitesPerChrom"})
        ),
        how="left"
    )
)
new_ac_at_snps_at_genes_with_only_old_sites_df

# %%
# (
#     new_ac_at_snps_at_genes_with_only_old_sites_df
#     .groupby("NumOfSuspectedSNPsPerChrom")
#     ["NumOfSuspectedSNPsPerChromWoACOrAT"].describe().round(2)D
# )

# %%
(
    new_ac_at_snps_at_genes_with_only_old_sites_df
    [["NumOfSuspectedSNPsPerChrom", "NumOfSuspectedSNPsPerChromWoACOrAT"]].value_counts()
    .reset_index(name="NumOfChroms")
)

# %%

# %%
(
    significant_mismatches_sites_disabled_by_suspected_snps_theoretical_fix_df
    .groupby("NumOfSuspectedSNPsPerChrom")
    ["EditedFinal"]
    .apply(
        lambda x: 100 * x.sum() / x.size
    )
    .reset_index(name="%OfRejectedEditingSitesPreviouslyFound")
)

# %%
fig = px.line(
    (
        significant_mismatches_sites_disabled_by_suspected_snps_theoretical_fix_df
        .groupby("NumOfSuspectedSNPsPerChrom")
        ["EditedFinal"]
        .apply(
            lambda x: 100 * x.sum() / x.size
        )
        .reset_index(name="%OfRejectedEditingSitesPreviouslyFound")
    ),
    x="NumOfSuspectedSNPsPerChrom",
    y="%OfRejectedEditingSitesPreviouslyFound",
    log_x=True,
    markers=True,
    # line_shape="linear",
    labels={
        "NumOfSuspectedSNPsPerChrom": "Number of suspected SNPs per gene",
        "%OfRejectedEditingSitesPreviouslyFound": "% of rejected<br>editing sites previously found"
    }
)
fig.update_layout(
    width=600,
    height=400
)
fig.show()

# %%

# %%

# %%
# (
#     significant_mismatches_theoretical_fix_chroms_with_too_much_suspected_snps_df.loc[
#         (
#             significant_mismatches_theoretical_fix_chroms_with_too_much_suspected_snps_df["Mismatch"]
#             (
#             significant_mismatches_theoretical_fix_chroms_with_too_much_suspected_snps_df["Mismatch"].eq("A>G")
#             & ~significant_mismatches_theoretical_fix_chroms_with_too_much_suspected_snps_df["SNP"]
#             & significant_mismatches_theoretical_fix_chroms_with_too_much_suspected_snps_df["AboveNewEditingThreshold"]
#         )
#         )
#     ]
# )

# %%

# %%

# %%
(
    concat_old_vs_new_editing_positions_df
    .loc[
        concat_old_vs_new_editing_positions_df["_merge"].eq("both"),
        ["Chrom", "Position", "OldMismatchFrequency", "NewMismatchFrequency"]
    ]
    .apply(
        lambda x: is_x_gt_or_eq_or_lt_y(x["OldMismatchFrequency"], x["NewMismatchFrequency"]),
        axis=1
    )
    .value_counts()
    .reset_index()
    .rename(
        columns={
            "index": "OldMismatchFrequency vs NewMismatchFrequency", 
        }
    )
)

# %%

# %%
# (
#     significant_mismatches_theoretical_fix_df.loc[
#         significant_mismatches_theoretical_fix_df["AboveNewEditingThreshold"]
#     ]
#     ["Mismatch"].value_counts(dropna=False)
# )

# %%
# (
#     significant_mismatches_theoretical_fix_df
#     .groupby(["EditedFinal", "AboveNoiseThreshold"])
#     ["Mismatch"].value_counts(dropna=False)
# )

# %%
significant_mismatches_theoretical_fix_df.head()

# %%
(
    significant_mismatches_theoretical_fix_df.loc[
        ~significant_mismatches_theoretical_fix_df["AboveNewEditingThreshold"],
        "MismatchFrequency"
    ]
    .describe().round(2)
)

# %%
(
    significant_mismatches_theoretical_fix_df.loc[
        (
            significant_mismatches_theoretical_fix_df["AboveNewEditingThreshold"]
            & ~significant_mismatches_theoretical_fix_df["SNP"]
        ),
        "MismatchFrequency"
    ]
    .describe().round(2)
)

# %%
(
    significant_mismatches_theoretical_fix_df.loc[
        (
            significant_mismatches_theoretical_fix_df["AboveNewEditingThreshold"]
            & ~significant_mismatches_theoretical_fix_df["AtOrAboveSuspectedSNPLevel"]
            # & significant_mismatches_theoretical_fix_df["Mismatch"].n("A>G")
        ),
        "MismatchFrequency"
    ]
)

# %%
df = (
        significant_mismatches_theoretical_fix_df.loc[
            ~significant_mismatches_theoretical_fix_df["AboveNewEditingThreshold"]
        ]
    )
fig = px.histogram(
    df,
    x="Mismatch",
    color="Mismatch",
    color_discrete_map=mismatch_dolor_map,
    # facet_col="Platform",
    # facet_col_spacing=0.04,
    log_y=True,
    template=template,
    category_orders={"Mismatch": mismatches},
    title="Mismatches at or below editing threshold",
)

width = 700
height = 500

# Use for_each_annotation to customize each title (i.e., remove the "Platform=" prefix)
fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))

fig.update_xaxes(tickangle=35)
# fig.update_yaxes(dtick=10)

fig.update_layout(
    width=width,
    height=height,
    showlegend=False
)

# fig.write_image(
#     Path(out_dir, "12 npn-SNP mismatches distribution - absolute - combined.svg"),
#     width=width,
#     height=height,
# )

display(
    (
        df.groupby("Mismatch").size()
        .reset_index(name="Count")
        .assign(
            Percentage=lambda x: np.round(
                100 * x["Count"] / x["Count"].sum(),
                0
            )
        )
    )
)

fig.show()

# %%
df = (
        significant_mismatches_theoretical_fix_df.loc[
            (
                significant_mismatches_theoretical_fix_df["AboveNewEditingThreshold"]
                & ~significant_mismatches_theoretical_fix_df["AtOrAboveSuspectedSNPLevel"]
            )
        ]
    )
fig = px.histogram(
    df,
    x="Mismatch",
    color="Mismatch",
    color_discrete_map=mismatch_dolor_map,
    # facet_col="Platform",
    # facet_col_spacing=0.04,
    log_y=True,
    template=template,
    category_orders={"Mismatch": mismatches},
    title="Mismatches above editing threshold, but below suspected SNP threshold",
)

width = 700
height = 500

# Use for_each_annotation to customize each title (i.e., remove the "Platform=" prefix)
fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))

fig.update_xaxes(tickangle=35)
# fig.update_yaxes(dtick=10)

fig.update_layout(
    width=width,
    height=height,
    showlegend=False
)

display(
    (
        df.groupby("Mismatch").size()
        .reset_index(name="Count")
        .assign(
            Percentage=lambda x: np.round(
                100 * x["Count"] / x["Count"].sum(),
                0
            )
        )
    )
)

# fig.write_image(
#     Path(out_dir, "12 npn-SNP mismatches distribution - absolute - combined.svg"),
#     width=width,
#     height=height,
# )



fig.show()

# %%
df = (
        significant_mismatches_theoretical_fix_df.loc[
            (
                significant_mismatches_theoretical_fix_df["AtOrAboveSuspectedSNPLevel"]
            )
        ]
    )
fig = px.histogram(
    df,
    x="Mismatch",
    color="Mismatch",
    color_discrete_map=mismatch_dolor_map,
    # facet_col="Platform",
    facet_col_spacing=0.04,
    facet_col="MismatchFrequency1",
    log_y=True,
    template=template,
    category_orders={"Mismatch": mismatches},
    title="Mismatches at or above suspected SNP threshold",
)

width = 1200
height = 500

# Use for_each_annotation to customize each title (i.e., remove the "Platform=" prefix)
fig.for_each_annotation(
    lambda a: a.update(
        text="Mismatch frequency = 100%" if a.text == "MismatchFrequency1=True" else
        "Mismatch frequency < 100%"
    )
)

fig.update_xaxes(tickangle=35)
# fig.update_yaxes(dtick=10)

fig.update_layout(
    width=width,
    height=height,
    showlegend=False
)

# fig.write_image(
#     Path(out_dir, "12 npn-SNP mismatches distribution - absolute - combined.svg"),
#     width=width,
#     height=height,
# )

display(
    (
        df.groupby(["MismatchFrequency1", "Mismatch"]).size()
        .reset_index(name="Count")
        .assign(
            Percentage=lambda x: np.round(
                100 * x["Count"] / x.groupby("MismatchFrequency1")["Count"].transform("sum"),
                0
            )
        )
    )
)

fig.show()

# %% [markdown]
# #### Original-threshold sensitivity analysis
#
# **O** is the set of original final editing sites (`EditedFinal=True`); **N** is the existing fully updated `EditingSite` set; **H** preserves all updated decisions and substitutes the original gene threshold for the final strict `>` comparison.
# A gene **permitted for detection** can have intermediate `Edited=True` positions and no final editing sites, so permission evidence and `OriginalHasFinalEditing` are reported separately.
# This revision corrects reporting populations and readability without changing the detection rules, parameters, frozen BH results, alternate-base choices, or mismatch-plot populations.
# Overlap fractions measure retention of calls, not biological accuracy.
#
# **Audit of original flag uses.** The uses below classify membership/counting, intermediate-stage diagnostics, and gene-permission inference; no global replacement of `Edited` is appropriate.
#
# | Use in this subsection | Classification | Interpretation |
# |---|---|---|
# | `ots_O_unscoped`, `ots_O`, original subset assertions: `EditedFinal` | Final membership/counting | Exactly the unique original final-site keys; the generic `ots_site_set` helper receives `EditedFinal` for O. |
# | Audit `OriginalEditingSite = OriginalEditedFinal.eq(True)`; all original loss, recovery, gene and scope counts derived from O or this flag | Final membership/counting | `OriginalHasFinalEditing` and `OriginalEditingSites` use final sites only. |
# | Copy original `Edited`, `EditedCorrected`, `EditedFinal` into the audit; discard the theoretical table's duplicated `EditedFinal` before joining | Diagnostic storage and final-status provenance | Preserve the original flags independently; attaching them never filters N or H. |
# | `OriginalEdited` → `OriginalNaiveEditedStatus` in gained-call origins | Intermediate-stage diagnostic | Passing the naive threshold/permission gate is not a final editing call. |
# | `OriginalEditedCorrected` → `OriginalBHStatus` and BH transition/p-value tables | Intermediate-stage diagnostic | Passing original editing BH alone is not a final editing call. |
# | `Edited`/`EditedFinal` in `ots_old_positive_genes`; original subset `Edited` assertion | Inference about original gene-filter permission | Positive intermediate `Edited` implies permission in `annotate_edited_sites`; the assertion checks that final calls in this dataset also pass that intermediate gate, not that it defines final calls. |
#
# Only this subsection is rerun using the prepared DataFrames. The paired Python file receives only the corresponding subsection revision.

# %% [markdown]
# **Prepared state.** This checks that the original and theoretical DataFrames, SNP settings, and plotting settings already exist. Missing prerequisites stop the subsection rather than launching the pipeline; numerical diagnostics use `atol=1e-12, rtol=1e-10`, which never alter strict calling.

# %%
import numpy as ots_np
import pandas as ots_pd
import plotly.express as ots_px
from IPython.display import display as ots_display, Markdown as ots_Markdown, HTML as ots_HTML

ots_required = [
    "concat_all_positions_df", "concat_all_edited_positions_df",
    "mismatches_theoretical_fix_df", "significant_mismatches_theoretical_fix_df",
    "significant_edited_positions_theoretical_fix_df", "noise_threshold_df",
    "new_editing_thresholds_df", "tmr50_alignment_stats_df",
    "max_snps_per_gene_to_allow_editing_detection", "snp_noise_level",
    "mismatches", "mismatch_dolor_map", "template",
]
ots_missing = [ots_name for ots_name in ots_required if ots_name not in globals()]
if ots_missing:
    raise RuntimeError(
        "Original-threshold analysis needs existing prepared objects: "
        + ", ".join(ots_missing)
        + ". Restore the prepared kernel/data; do not run the whole notebook. "
        "Prerequisites: original positions/calls, TMR50 alignment metadata, original "
        "thresholds, frozen full theoretical/BH table, significant annotations/new "
        "thresholds/calls, and the existing SNP and plot settings."
    )
ots_keys = ["Chrom", "Position"]
ots_atol, ots_rtol = 1e-12, 1e-10

# %% [markdown]
# **Frozen inputs and site keys.** The input population is the full updated significant mismatch table in the TMR50 gene universe. A copy preserves its mismatch selections and annotations; small existing helpers define unique site keys, explicitly denominated ratios, and tolerance comparisons.

# %%
ots_hybrid_df = significant_mismatches_theoretical_fix_df.copy(deep=True)
ots_gene_universe = ots_pd.Index(tmr50_alignment_stats_df["Chrom"].unique(), name="Chrom")
assert ots_hybrid_df["Chrom"].isin(ots_gene_universe).all()
assert mismatches_theoretical_fix_df["Chrom"].isin(ots_gene_universe).all()


def ots_site_set(ots_frame, ots_flag=None):
    ots_rows = ots_frame if ots_flag is None else ots_frame.loc[ots_frame[ots_flag].eq(True)]
    return set(ots_rows[ots_keys].itertuples(index=False, name=None))


def ots_ratio(ots_numerator, ots_denominator):
    return ots_numerator / ots_denominator if ots_denominator else ots_np.nan


def ots_close(ots_first, ots_second):
    return ots_np.isclose(ots_first, ots_second, atol=ots_atol, rtol=ots_rtol, equal_nan=False)


# %% [markdown]
# **Independent threshold annotations.** Original thresholds are checked against one row per gene in `noise_threshold_df`, and updated thresholds against `new_editing_thresholds_df`. Original `NoiseThreshold` already contains the safety factor, so it is neither multiplied again nor filled with zero when missing.

# %%
# Duplicate sites would make both audit joins and SNP counts ambiguous: fail early.
ots_validation_df = ots_pd.DataFrame([
    {"Object": ots_name, "Rows": len(globals()[ots_name]),
     "Duplicate site keys": int(globals()[ots_name].duplicated(ots_keys).sum())}
    for ots_name in ots_required[:5]
])
assert ots_validation_df["Duplicate site keys"].eq(0).all(), ots_validation_df
ots_original_thresholds = noise_threshold_df[["Chrom", "NoiseThreshold"]].copy()
ots_updated_thresholds = new_editing_thresholds_df[["Chrom", "NewEditingThreshold"]].copy()
for ots_threshold_table, ots_threshold_col in [
    (ots_original_thresholds, "NoiseThreshold"),
    (ots_updated_thresholds, "NewEditingThreshold"),
]:
    assert not ots_threshold_table["Chrom"].duplicated().any()
    assert ots_threshold_table[ots_threshold_col].notna().all()
    assert ots_np.isfinite(ots_threshold_table[ots_threshold_col]).all()
    assert ots_threshold_table[ots_threshold_col].ge(0).all()
ots_thresholds_df = (
    ots_pd.DataFrame(index=ots_gene_universe)
    .join(ots_original_thresholds.set_index("Chrom"))
    .join(ots_updated_thresholds.set_index("Chrom"))
)
assert ots_thresholds_df["NoiseThreshold"].notna().all(), "Missing original threshold; no zero imputation"
for ots_threshold_col in ["NoiseThreshold", "NewEditingThreshold"]:
    ots_mapped_threshold = ots_hybrid_df["Chrom"].map(ots_thresholds_df[ots_threshold_col])
    assert ots_mapped_threshold.notna().all()
    assert ots_hybrid_df[ots_threshold_col].eq(ots_mapped_threshold).all(), (
        "Threshold disagrees with independent per-Chrom table", ots_threshold_col
    )
assert mismatches_theoretical_fix_df["NoiseThreshold"].eq(
    mismatches_theoretical_fix_df["Chrom"].map(ots_thresholds_df["NoiseThreshold"])
).all()

# %% [markdown]
# **Hybrid calls with frozen updated decisions.** This verifies the existing BH population and mismatch identities before changing only the final threshold comparison on the copy. H uses strict `MismatchFrequency > NoiseThreshold`; neither H nor N requires original `EditedFinal`, allowing novel final calls.

# %%
assert ots_hybrid_df["BHRejection"].eq(True).all()
assert ots_site_set(mismatches_theoretical_fix_df, "BHRejection") == ots_site_set(ots_hybrid_df)
ots_frozen_cols = ["Mismatch", "MismatchFrequency", "BinomPVal", "BHCorrectedPVal", "BHRejection"]
ots_pd.testing.assert_frame_equal(
    mismatches_theoretical_fix_df.loc[mismatches_theoretical_fix_df["BHRejection"], ots_keys + ots_frozen_cols]
        .set_index(ots_keys).sort_index(),
    ots_hybrid_df[ots_keys + ots_frozen_cols].set_index(ots_keys).sort_index(),
    check_dtype=False, check_exact=True,
)
ots_hybrid_df["AboveOriginalEditingThreshold"] = (
    ots_hybrid_df["MismatchFrequency"] > ots_hybrid_df["NoiseThreshold"]
)
ots_hybrid_df["EditingSiteUsingOriginalThreshold"] = (
    ots_hybrid_df["Mismatch"].eq("A>G")
    & ~ots_hybrid_df["SNP"]
    & ~ots_hybrid_df["EditingDetectionDisabledInChromDueToSuspectedSNPs"]
    & ots_hybrid_df["AboveOriginalEditingThreshold"]
)
assert ots_hybrid_df["AboveNewEditingThreshold"].equals(
    ots_hybrid_df["MismatchFrequency"] > ots_hybrid_df["NewEditingThreshold"]
)
assert ots_hybrid_df["EditingSite"].equals(
    ots_hybrid_df["Mismatch"].eq("A>G") & ~ots_hybrid_df["SNP"]
    & ~ots_hybrid_df["EditingDetectionDisabledInChromDueToSuspectedSNPs"]
    & ots_hybrid_df["AboveNewEditingThreshold"]
)

# %% [markdown]
# **Original final-site membership.** O is constructed only from `EditedFinal=True`, and direct assertions compare both O and `concat_all_edited_positions_df` with independently extracted unique final-site keys. A separate `Edited` assertion checks permission evidence for those final calls; it does not select or count them.

# %%
ots_O_unscoped = ots_site_set(concat_all_positions_df, "EditedFinal")
ots_O = {ots_key for ots_key in ots_O_unscoped if ots_key[0] in ots_gene_universe}
ots_N = ots_site_set(ots_hybrid_df, "EditingSite")
ots_H = ots_site_set(ots_hybrid_df, "EditingSiteUsingOriginalThreshold")
# Independent final-membership checks; intermediate Edited is not the selector.
ots_original_final_keys = set(concat_all_positions_df.loc[
    concat_all_positions_df["EditedFinal"].eq(True), ots_keys
].drop_duplicates().itertuples(index=False, name=None))
assert ots_O_unscoped == ots_original_final_keys
assert ots_O == {ots_key for ots_key in ots_original_final_keys if ots_key[0] in ots_gene_universe}
assert concat_all_edited_positions_df["EditedFinal"].eq(True).all()
assert ots_site_set(concat_all_edited_positions_df) == ots_original_final_keys
# Permission evidence only: this assertion never supplies final-site membership.
assert concat_all_edited_positions_df["Edited"].eq(True).all(), "Final calls lack expected intermediate permission evidence"
assert ots_O_unscoped == ots_site_set(concat_all_edited_positions_df)
assert ots_N == ots_site_set(significant_edited_positions_theoretical_fix_df)
ots_sets = {"O": ots_O, "N": ots_N, "H": ots_H}

# %% [markdown]
# **All-TMR50 baseline.** This compares measured original/updated site counts and overlap with the saved baseline, and reports H's size in the same gene universe. Differences are displayed as regression observations rather than corrected by changing data.

# %%
ots_baseline_df = ots_pd.DataFrame({
    "Metric": ["Original sites", "Updated sites", "Shared O & N", "O minus N", "N minus O"],
    "Saved baseline": [11711, 10714, 10403, 1308, 311],
    "Measured TMR50": [len(ots_O), len(ots_N), len(ots_O & ots_N), len(ots_O - ots_N), len(ots_N - ots_O)],
})
ots_baseline_df["Difference from saved"] = ots_baseline_df["Measured TMR50"] - ots_baseline_df["Saved baseline"]
ots_display(ots_baseline_df)
print("Original calls outside the fixed TMR50 scope:", len(ots_O_unscoped - ots_O))
print("Hybrid calls:", len(ots_H), "; unique TMR50 genes:", len(ots_gene_universe))
ots_display(ots_validation_df)

# %% [markdown]
# **Genes with final editing versus updated permission.** For every TMR50 gene, `OriginalHasFinalEditing` is true exactly when O contains at least one of its sites, and `OriginalEditingSites` counts those final sites. Independently, `UpdatedEligible` applies the unchanged suspected-SNP limit to the full significant population, including zero counts; no updated final call is needed.

# %%
ots_gene_df = ots_thresholds_df.copy()
ots_gene_df["UpdatedSuspectedSNPs"] = (
    ots_hybrid_df.groupby("Chrom")["SuspectedSNP"].sum().reindex(ots_gene_universe, fill_value=0).astype(int)
)
ots_gene_df["UpdatedEligible"] = ots_gene_df["UpdatedSuspectedSNPs"].le(max_snps_per_gene_to_allow_editing_detection)
assert (~ots_hybrid_df["EditingDetectionDisabledInChromDueToSuspectedSNPs"]).eq(
    ots_hybrid_df["Chrom"].map(ots_gene_df["UpdatedEligible"])
).all()
ots_gene_df["OriginalHasFinalEditing"] = ots_gene_df.index.isin({ots_key[0] for ots_key in ots_O})
ots_gene_df["OriginalEditingSites"] = ots_pd.Series([ots_key[0] for ots_key in ots_O]).value_counts().reindex(ots_gene_universe, fill_value=0)
assert ots_gene_df["OriginalHasFinalEditing"].equals(ots_gene_df["OriginalEditingSites"].gt(0))

# %% [markdown]
# **Evidence of original detection permission.** In [annotate_edited_sites](../Pileup/positions.py), `editing_detection_possible=False` forces every intermediate `Edited` flag to false, so a positive `Edited` flag is evidence that the gene filter permitted detection. A retained original SNP count above the limit proves exclusion, but a lower count is only a coverage-filtered lower bound; genes with no positive evidence remain unknown, irrespective of final-call absence.

# %%
ots_old_positive_genes = set(concat_all_positions_df.loc[
    concat_all_positions_df["Edited"].eq(True) | concat_all_positions_df["EditedFinal"].eq(True), "Chrom"
])
ots_old_snp_lower_bound = concat_all_positions_df.loc[
    concat_all_positions_df["NoisyFinal"].eq(True) & concat_all_positions_df["Noise"].ge(snp_noise_level)
].groupby("Chrom").size().reindex(ots_gene_universe, fill_value=0)
ots_gene_df["OriginalSNPCountLowerBound"] = ots_old_snp_lower_bound
ots_gene_df["OriginalEligibility"] = "unknown: only coverage-filtered evidence"
ots_gene_df.loc[ots_gene_df.index.isin(ots_old_positive_genes), "OriginalEligibility"] = "permitted: positive original intermediate Edited flag"
ots_old_proven_excluded = ots_old_snp_lower_bound.gt(max_snps_per_gene_to_allow_editing_detection)
assert not (ots_old_proven_excluded & ots_gene_df.index.isin(ots_old_positive_genes)).any()
ots_gene_df.loc[ots_old_proven_excluded, "OriginalEligibility"] = "ineligible: retained SNP count exceeds limit"
ots_gene_df["OriginalEligible"] = ots_pd.Series(ots_pd.NA, index=ots_gene_df.index, dtype="boolean")
ots_gene_df.loc[ots_gene_df.index.isin(ots_old_positive_genes), "OriginalEligible"] = True
ots_gene_df.loc[ots_old_proven_excluded, "OriginalEligible"] = False

# %% [markdown]
# **Primary and secondary reporting scopes.** The primary restriction is exactly `OriginalHasFinalEditing & UpdatedEligible`; genes that lose all updated final calls remain in its denominator. The broader “Established eligible in both” scope is retained as a secondary diagnostic and can include genes with intermediate `Edited=True` but no `EditedFinal=True`; `ots_scope_regression_df` reports both denominators.

# %%
ots_common_genes = set(ots_gene_df.index[
    ots_gene_df["OriginalEligible"].eq(True).fillna(False) & ots_gene_df["UpdatedEligible"]
])
ots_gene_eligibility_counts = ots_gene_df.groupby(["OriginalEligibility", "UpdatedEligible"]).size().rename("Genes").reset_index()
ots_display(ots_gene_eligibility_counts)
ots_primary_scope = "Genes with original final editing that pass the updated SNP gene filter"
ots_primary_genes = set(ots_gene_df.index[
    ots_gene_df["OriginalHasFinalEditing"] & ots_gene_df["UpdatedEligible"]
])
ots_scope_definitions = [
    ("All TMR50", set(ots_gene_universe), "All genes"),
    (ots_primary_scope, ots_primary_genes, "Primary restricted comparison"),
    ("Established eligible in both", ots_common_genes, "Secondary permission-evidence diagnostic"),
]
ots_scope_regression_df = ots_pd.DataFrame({
    "Population": ["Genes with original final editing", ots_primary_scope, "Established eligible in both"],
    "Saved expectation": [1262, 1159, 1532],
    "Measured genes": [int(ots_gene_df["OriginalHasFinalEditing"].sum()), len(ots_primary_genes), len(ots_common_genes)],
})
ots_scope_regression_df["Difference"] = ots_scope_regression_df["Measured genes"] - ots_scope_regression_df["Saved expectation"]
ots_display(ots_scope_regression_df)
print("Secondary permission-evidence genes without original final editing:", len(ots_common_genes - ots_primary_genes))


# %% [markdown]
# **Pairwise denominators.** For each fixed reporting scope, the following helper reports both set sizes, their intersection, each directional difference, and overlap fractions with explicit denominators. Jaccard uses the union of called sites, without adding non-editing positions to the denominator.

# %%
def ots_pairwise(ots_first, ots_second, ots_first_label, ots_second_label, ots_scope):
    ots_shared = len(ots_first & ots_second)
    return {
        "Scope": ots_scope, "First": ots_first_label, "Second": ots_second_label,
        "First sites": len(ots_first), "Second sites": len(ots_second),
        "Intersection": ots_shared, "Union sites": len(ots_first | ots_second), "First minus second": len(ots_first - ots_second),
        "Second minus first": len(ots_second - ots_first),
        "Intersection / first sites": ots_ratio(ots_shared, len(ots_first)),
        "Intersection / second sites": ots_ratio(ots_shared, len(ots_second)),
        "Jaccard: intersection / union": ots_ratio(ots_shared, len(ots_first | ots_second)),
    }


# %% [markdown]
# **O/N/H comparisons within each scope.** The same three site sets are restricted only by the stated gene population, and all three pairwise comparisons are recomputed. `ots_pairwise_df` reports population genes, genes with final calls, site counts, and overlaps, while `ots_gene_counts_df` gives one row per scheme and scope.

# %%
ots_pair_rows, ots_gene_count_rows = [], []
for ots_scope, ots_scope_genes, ots_role in ots_scope_definitions:
    ots_scope_sets = {ots_label: {ots_key for ots_key in ots_sites if ots_key[0] in ots_scope_genes}
                      for ots_label, ots_sites in ots_sets.items()}
    for ots_first_label, ots_second_label in [("O", "N"), ("O", "H"), ("N", "H")]:
        ots_pair_rows.append({
            **ots_pairwise(ots_scope_sets[ots_first_label], ots_scope_sets[ots_second_label],
                           ots_first_label, ots_second_label, ots_scope),
            "Reporting role": ots_role, "Population genes": len(ots_scope_genes),
            "First genes with final calls": len({ots_key[0] for ots_key in ots_scope_sets[ots_first_label]}),
            "Second genes with final calls": len({ots_key[0] for ots_key in ots_scope_sets[ots_second_label]}),
        })
    for ots_label, ots_sites in ots_scope_sets.items():
        ots_gene_count_rows.append({"Scope": ots_scope, "Reporting role": ots_role, "Scheme": ots_label,
                                   "Universe genes": len(ots_scope_genes), "Editing sites": len(ots_sites),
                                   "Genes with final editing": len({ots_key[0] for ots_key in ots_sites})})
ots_pairwise_df = ots_pd.DataFrame(ots_pair_rows)
ots_gene_counts_df = ots_pd.DataFrame(ots_gene_count_rows)
ots_display(ots_HTML(ots_pairwise_df.to_html(index=False)))
ots_display(ots_HTML(ots_gene_counts_df.to_html(index=False)))

# %% [markdown]
# **Complete original-versus-theoretical audit.** This starts with all theoretical candidates, including BH failures, and outer-joins retained original positions; it never starts only with final calls. Both gene thresholds are attached independently of call status, with missing updated thresholds left explicit in genes without significant mismatches; `ots_audit_df` therefore retains original-only sites and their available evidence.

# %%
ots_original_cols = [ots_col for ots_col in [
    "Transcript", "RefBase", "TotalCoverage", "A", "T", "C", "G", "EditingFrequency",
    "Edited", "EditedCorrected", "EditedFinal", "EditingBinomPVal", "EditingCorrectedPVal",
    "Noise", "NoisyCorrected", "NoisyFinal", "NoiseBinomPVal", "NoiseCorrectedPVal",
    "BelowNoiseFreq1", "BelowEditingFreq1", "CDS", "KnownEditing", "InProbRegion",
] if ots_col in concat_all_positions_df.columns]
ots_original_audit = concat_all_positions_df.loc[
    concat_all_positions_df["Chrom"].isin(ots_gene_universe), ots_keys + ots_original_cols
].rename(columns={ots_col: "Original" + ots_col for ots_col in ots_original_cols}).set_index(ots_keys)
ots_original_audit["PresentInOriginalPositions"] = True
ots_audit_df = mismatches_theoretical_fix_df.drop(columns=["NoiseThreshold", "EditedFinal", "NoisyFinal"]).copy().set_index(ots_keys)
ots_audit_df["PresentInTheoreticalUniverse"] = True
ots_audit_df = ots_audit_df.join(ots_original_audit, how="outer", validate="one_to_one").reset_index()
for ots_col in ["PresentInOriginalPositions", "PresentInTheoreticalUniverse"]:
    ots_audit_df[ots_col] = ots_audit_df[ots_col].eq(True)
ots_audit_df["NoiseThreshold"] = ots_audit_df["Chrom"].map(ots_thresholds_df["NoiseThreshold"])
ots_audit_df["NewEditingThreshold"] = ots_audit_df["Chrom"].map(ots_thresholds_df["NewEditingThreshold"])
ots_audit_df["UpdatedThresholdAvailable"] = ots_audit_df["NewEditingThreshold"].notna()
assert ots_audit_df["NoiseThreshold"].notna().all()

# %% [markdown]
# **Final status and gene evidence on each site.** Original final membership is exactly `OriginalEditedFinal.eq(True)`, while the unchanged N/H flags come from the significant table. Intermediate SNP annotations remain unknown outside that population, and final-editing presence, original permission evidence, and updated permission remain separate gene columns.

# %%
ots_audit_annotation_cols = [
    "AtOrAboveSuspectedSNPLevel", "MismatchFrequency1", "SuspectedSNP", "DefinitiveSNP", "SNP",
    "NoiseSite", "EditingSite", "EditingSiteUsingOriginalThreshold",
]
ots_audit_df = ots_audit_df.merge(ots_hybrid_df[ots_keys + ots_audit_annotation_cols], on=ots_keys, how="left", validate="one_to_one")
for ots_col in ["AtOrAboveSuspectedSNPLevel", "MismatchFrequency1", "SuspectedSNP", "DefinitiveSNP", "SNP", "NoiseSite"]:
    ots_audit_df[ots_col] = ots_audit_df[ots_col].astype("boolean")
ots_audit_df["OriginalEditingSite"] = ots_audit_df["OriginalEditedFinal"].eq(True)
for ots_col in ["EditingSite", "EditingSiteUsingOriginalThreshold"]:
    ots_audit_df[ots_col] = ots_audit_df[ots_col].eq(True)
ots_audit_df["OriginalHasFinalEditing"] = ots_audit_df["Chrom"].map(ots_gene_df["OriginalHasFinalEditing"])
ots_audit_df["OriginalEligible"] = ots_audit_df["Chrom"].map(ots_gene_df["OriginalEligible"])
ots_audit_df["OriginalEligibility"] = ots_audit_df["Chrom"].map(ots_gene_df["OriginalEligibility"])
ots_audit_df["UpdatedEligible"] = ots_audit_df["Chrom"].map(ots_gene_df["UpdatedEligible"])
ots_audit_df["EditingDetectionDisabledInChromDueToSuspectedSNPs"] = ~ots_audit_df["UpdatedEligible"]
ots_audit_df["AbsentOriginalPositionIsNotProofOfZeroCoverage"] = ~ots_audit_df["PresentInOriginalPositions"]
ots_audit_df["AboveOriginalEditingThreshold"] = ots_audit_df["MismatchFrequency"].gt(ots_audit_df["NoiseThreshold"])
ots_audit_df["AboveNewEditingThreshold"] = ots_audit_df["MismatchFrequency"].gt(ots_audit_df["NewEditingThreshold"])

# %% [markdown]
# **Independent failure flags.** For each audited candidate, mismatch identity, updated BH failure, SNP exclusion, gene exclusion, and threshold failure are marked separately wherever observable. Missing candidates or annotations have explicit flags; absence from retained original data is not proof of zero pre-filter coverage, and numerical failure is not inferred from a missing threshold.

# %%
ots_audit_df["AbsentTheoreticalCandidate"] = ~ots_audit_df["PresentInTheoreticalUniverse"]
ots_audit_df["LeadingMismatchNotAG"] = ots_audit_df["Mismatch"].notna() & ots_audit_df["Mismatch"].ne("A>G")
ots_audit_df["NewBHFailure"] = ots_audit_df["PresentInTheoreticalUniverse"] & ots_audit_df["BHRejection"].eq(False)
ots_audit_df["SNPExclusion"] = ots_audit_df["SNP"].eq(True).fillna(False)
ots_audit_df["UpdatedGeneExcluded"] = ~ots_audit_df["UpdatedEligible"]
ots_audit_df["UpdatedThresholdFailure"] = (
    ots_audit_df["MismatchFrequency"].notna() & ots_audit_df["NewEditingThreshold"].notna()
    & ~ots_audit_df["AboveNewEditingThreshold"]
)
ots_audit_df["OriginalThresholdFailure"] = (
    ots_audit_df["MismatchFrequency"].notna() & ots_audit_df["NoiseThreshold"].notna()
    & ~ots_audit_df["AboveOriginalEditingThreshold"]
)
ots_audit_df["MismatchUnavailable"] = ots_audit_df["PresentInTheoreticalUniverse"] & ots_audit_df["Mismatch"].isna()
ots_audit_df["UpdatedThresholdUnavailable"] = ~ots_audit_df["UpdatedThresholdAvailable"]

# %% [markdown]
# **Lost sites and exhaustive combinations.** The input populations are O minus N, O minus H, recovered original losses, and N minus H. Overlapping reason flags are retained alongside mutually exclusive combinations, and assertions require every lost site to appear exactly once with no unexplained combination.

# %%
ots_lost_df = ots_audit_df.loc[ots_audit_df["OriginalEditingSite"] & ~ots_audit_df["EditingSite"]].copy()
ots_remaining_df = ots_audit_df.loc[ots_audit_df["OriginalEditingSite"] & ~ots_audit_df["EditingSiteUsingOriginalThreshold"]].copy()
ots_recovered_df = ots_lost_df.loc[ots_lost_df["EditingSiteUsingOriginalThreshold"]].copy()
ots_n_to_h_lost_df = ots_audit_df.loc[ots_audit_df["EditingSite"] & ~ots_audit_df["EditingSiteUsingOriginalThreshold"]].copy()
ots_reason_cols = ["AbsentTheoreticalCandidate", "LeadingMismatchNotAG", "NewBHFailure", "SNPExclusion",
                   "UpdatedGeneExcluded", "UpdatedThresholdFailure", "MismatchUnavailable", "UpdatedThresholdUnavailable"]
ots_h_reason_cols = [ots_col for ots_col in ots_reason_cols if ots_col not in ["UpdatedThresholdFailure", "UpdatedThresholdUnavailable"]] + ["OriginalThresholdFailure"]


def ots_reason_summary(ots_frame, ots_columns):
    if ots_frame.empty:
        return ots_pd.DataFrame(columns=["Mutually exclusive reason combination", "Sites"])
    ots_combinations = ots_frame[ots_columns].apply(
        lambda ots_row: " + ".join(ots_col for ots_col in ots_columns if ots_row[ots_col]) or "Unexplained (inspect)", axis=1
    ).value_counts().rename_axis("Mutually exclusive reason combination").reset_index(name="Sites")
    assert ots_combinations["Sites"].sum() == len(ots_frame)
    assert not ots_combinations["Mutually exclusive reason combination"].eq("Unexplained (inspect)").any(), ots_combinations
    return ots_combinations


ots_lost_reasons_df = ots_lost_df[ots_reason_cols].sum().rename("Lost original sites (overlapping)").rename_axis("Reason").reset_index()
ots_lost_combinations_df = ots_reason_summary(ots_lost_df, ots_reason_cols)
ots_remaining_combinations_df = ots_reason_summary(ots_remaining_df, ots_h_reason_cols)
assert len(ots_lost_df) == len(ots_O - ots_N)
assert len(ots_remaining_df) == len(ots_O - ots_H)
assert len(ots_recovered_df) == len((ots_O - ots_N) & ots_H)
assert ots_site_set(ots_audit_df, "OriginalEditingSite") == ots_O
assert ots_site_set(ots_audit_df, "EditingSite") == ots_N
assert ots_site_set(ots_audit_df, "EditingSiteUsingOriginalThreshold") == ots_H

# %% [markdown]
# **Loss and recovery summary.** `ots_loss_recovery_df` counts losses, recoveries, residual losses, and any N-to-H losses using final-site keys. The overlapping reason totals and exhaustive combination tables explain those counts; the audit preview confirms independent threshold annotations on original losses.

# %%
ots_loss_recovery_df = ots_pd.DataFrame([
    {"Measure": "Original losses O minus N", "Sites": len(ots_lost_df)},
    {"Measure": "Recovered original losses (O minus N) & H", "Sites": len(ots_recovered_df)},
    {"Measure": "Still missing from H: O minus H", "Sites": len(ots_remaining_df)},
    {"Measure": "Lost from N to H", "Sites": len(ots_n_to_h_lost_df)},
])
ots_display(ots_loss_recovery_df)
print("Recovered / original losses:", ots_ratio(len(ots_recovered_df), len(ots_lost_df)))
ots_display(ots_lost_reasons_df)
ots_display(ots_HTML(ots_lost_combinations_df.to_html(index=False)))
ots_display(ots_HTML(ots_remaining_combinations_df.to_html(index=False)))
ots_display(ots_reason_summary(ots_n_to_h_lost_df, ots_h_reason_cols))
print("Full audit rows:", len(ots_audit_df), "; both thresholds available on lost originals:",
      int(ots_lost_df[["NoiseThreshold", "NewEditingThreshold"]].notna().all(axis=1).sum()))
ots_display(ots_lost_df[ots_keys + ["OriginalEditingFrequency", "Mismatch", "MismatchFrequency", "BHRejection",
                                 "NoiseThreshold", "NewEditingThreshold", "EditingSiteUsingOriginalThreshold"] + ots_reason_cols].head(12))

# %% [markdown]
# **Origins of gained final calls.** This retains updated-only, hybrid-only, recovered, and N-to-H lost calls and groups them by available original intermediate-stage evidence. `OriginalNaiveEditedStatus` describes `OriginalEdited`, `OriginalBHStatus` describes `OriginalEditedCorrected`, and only `OriginalEditingSite` indicates original final editing; `ots_gain_origins_df` is a reporting table, not a reconstructed original counterfactual.

# %%
ots_gain_sites_df = ots_audit_df.loc[
    (ots_audit_df["EditingSite"] | ots_audit_df["EditingSiteUsingOriginalThreshold"])
    & ~(ots_audit_df["OriginalEditingSite"] & ots_audit_df["EditingSite"] & ots_audit_df["EditingSiteUsingOriginalThreshold"])
].copy()
ots_gain_sites_df["OriginalBHStatus"] = ots_gain_sites_df["OriginalEditedCorrected"].map({True: "pass", False: "fail"}).fillna("unavailable")
ots_gain_sites_df["OriginalNaiveEditedStatus"] = ots_gain_sites_df["OriginalEdited"].map({True: "pass", False: "fail"}).fillna("unavailable")
ots_gain_sites_df["OriginalFrequencyAboveThreshold"] = ots_gain_sites_df["OriginalEditingFrequency"].gt(ots_gain_sites_df["NoiseThreshold"])
ots_gain_sites_df["FrequencyAboveOriginalThresholdChanged"] = (
    ots_gain_sites_df["OriginalEditingFrequency"].notna()
    & ots_gain_sites_df["OriginalFrequencyAboveThreshold"].ne(ots_gain_sites_df["AboveOriginalEditingThreshold"])
)
ots_gain_summary_parts = []
for ots_label, ots_gain_set in [("N minus O", ots_N - ots_O), ("H minus O", ots_H - ots_O),
                               ("H minus N", ots_H - ots_N), ("N minus H", ots_N - ots_H)]:
    ots_gain_part = ots_gain_sites_df.loc[ots_pd.MultiIndex.from_frame(ots_gain_sites_df[ots_keys]).isin(ots_gain_set)].copy()
    ots_gain_part["Comparison"] = ots_label
    ots_gain_summary_parts.append(ots_gain_part.groupby([
        "Comparison", "OriginalEditingSite", "PresentInOriginalPositions", "OriginalBHStatus",
        "OriginalNaiveEditedStatus", "OriginalEligibility", "FrequencyAboveOriginalThresholdChanged",
    ], dropna=False).size().reset_index(name="Sites"))
    assert len(ots_gain_part) == len(ots_gain_set)
ots_gain_origins_df = ots_pd.concat(ots_gain_summary_parts, ignore_index=True)
ots_display(ots_gain_origins_df)

# %% [markdown]
# **AC/AT diagnostic 1: threshold candidates only.** The population is the same updated significant `NoiseSite` candidates, first unchanged and then excluding A>C/A>T only from threshold candidacy. The top-three calculation preserves zero padding and factor 1.5; position breaks equal-frequency contributor ordering without reselecting any frozen alternate base, and tied contributors are flagged.

# %%
ots_noise_candidates = ots_hybrid_df.loc[ots_hybrid_df["NoiseSite"], ots_keys + ["Mismatch", "MismatchFrequency"]].copy()
ots_noise_candidates["ACorAT"] = ots_noise_candidates["Mismatch"].isin(["A>C", "A>T"])
ots_represented_genes = ots_pd.Index(ots_hybrid_df["Chrom"].unique(), name="Chrom")


def ots_top_three(ots_candidates):
    ots_ranked = ots_candidates.sort_values(["Chrom", "MismatchFrequency", "Position"], ascending=[True, False, True], kind="stable").copy()
    ots_ranked["TiedFrequencyCandidates"] = ots_ranked.groupby(["Chrom", "MismatchFrequency"])["Position"].transform("size")
    ots_ranked["Rank"] = ots_ranked.groupby("Chrom").cumcount() + 1
    ots_grid = ots_pd.MultiIndex.from_product([ots_represented_genes, [1, 2, 3]], names=["Chrom", "Rank"])
    ots_top = ots_ranked.loc[ots_ranked["Rank"].le(3)].set_index(["Chrom", "Rank"]).reindex(ots_grid).reset_index()
    ots_top["ZeroPadding"] = ots_top["Position"].isna()
    ots_top["ThresholdContributionFrequency"] = ots_top["MismatchFrequency"].fillna(0.0)
    ots_top["ACorAT"] = ots_top["ACorAT"].eq(True)
    ots_top["TopThreeFrequencyTie"] = ots_top["TiedFrequencyCandidates"].gt(1)
    return ots_top


# %% [markdown]
# **Threshold contributions and the saved comparison.** Per-gene thresholds are computed from the two candidate lists and checked against the existing updated thresholds. `ots_gene_df` stores their difference as the direct AC/AT threshold contribution conditional on updated BH, while `ots_threshold_baseline_df` compares exact old/new inequalities with the saved counts.

# %%
ots_top_three_df = ots_top_three(ots_noise_candidates)
ots_top_three_without_acat_df = ots_top_three(ots_noise_candidates.loc[~ots_noise_candidates["ACorAT"]])
ots_gene_df["RecomputedUpdatedThreshold"] = ots_top_three_df.groupby("Chrom")["ThresholdContributionFrequency"].mean().mul(1.5)
ots_gene_df["ThresholdWithoutACATConditionalOnUpdatedBH"] = ots_top_three_without_acat_df.groupby("Chrom")["ThresholdContributionFrequency"].mean().mul(1.5)
ots_gene_df["TopThreeACATContributors"] = ots_top_three_df.groupby("Chrom")["ACorAT"].sum()
ots_gene_df["ACATThresholdContribution"] = ots_gene_df["RecomputedUpdatedThreshold"] - ots_gene_df["ThresholdWithoutACATConditionalOnUpdatedBH"]
ots_gene_df["NewMinusOriginalThreshold"] = ots_gene_df["NewEditingThreshold"] - ots_gene_df["NoiseThreshold"]
ots_gene_df["LostOriginalSites"] = ots_lost_df.groupby("Chrom").size().reindex(ots_gene_universe, fill_value=0)
assert ots_close(ots_gene_df.loc[ots_represented_genes, "RecomputedUpdatedThreshold"],
                 ots_gene_df.loc[ots_represented_genes, "NewEditingThreshold"]).all()
assert ots_gene_df.loc[ots_represented_genes, "ACATThresholdContribution"].ge(-ots_atol).all()
ots_gene_df["ThresholdsClose"] = ots_close(ots_gene_df["NoiseThreshold"], ots_gene_df["NewEditingThreshold"])
ots_threshold_baseline_df = ots_pd.DataFrame({
    "Comparison": ["equal (exact)", "new higher (exact)", "new lower (exact)"],
    "Saved": [2594, 1153, 6],
    "Measured": [int(ots_gene_df["NewEditingThreshold"].eq(ots_gene_df["NoiseThreshold"]).sum()),
                 int(ots_gene_df["NewEditingThreshold"].gt(ots_gene_df["NoiseThreshold"]).sum()),
                 int(ots_gene_df["NewEditingThreshold"].lt(ots_gene_df["NoiseThreshold"]).sum())],
})
ots_threshold_baseline_df["Difference from saved"] = ots_threshold_baseline_df["Measured"] - ots_threshold_baseline_df["Saved"]
ots_display(ots_threshold_baseline_df)

# %% [markdown]
# **Threshold changes in three gene populations.** The existing summaries cover all represented genes, genes with original final editing, and genes with lost original final sites. Exact changes and tolerance-based changes remain separate, and `ots_threshold_summary_df` reports the prevalence and size of direct AC/AT contributions within each population.

# %%
ots_threshold_summary_rows = []
for ots_scope, ots_scope_mask in [
    ("All represented genes", ots_gene_df.index.isin(ots_represented_genes)),
    ("Genes with original editing", ots_gene_df["OriginalEditingSites"].gt(0)),
    ("Genes with lost original sites", ots_gene_df["LostOriginalSites"].gt(0)),
]:
    ots_threshold_part = ots_gene_df.loc[ots_scope_mask & ots_gene_df["NewEditingThreshold"].notna()]
    ots_threshold_summary_rows.append({
        "Scope": ots_scope, "Genes": len(ots_threshold_part),
        "Threshold unavailable": int((ots_scope_mask & ots_gene_df["NewEditingThreshold"].isna()).sum()),
        "Equal within tolerance": int(ots_threshold_part["ThresholdsClose"].sum()),
        "Meaningfully higher": int((~ots_threshold_part["ThresholdsClose"] & ots_threshold_part["NewMinusOriginalThreshold"].gt(0)).sum()),
        "Meaningfully lower": int((~ots_threshold_part["ThresholdsClose"] & ots_threshold_part["NewMinusOriginalThreshold"].lt(0)).sum()),
        "Genes with top-three AC/AT": int(ots_threshold_part["TopThreeACATContributors"].gt(0).sum()),
        "Threshold reduced excluding AC/AT": int((~ots_close(ots_threshold_part["RecomputedUpdatedThreshold"], ots_threshold_part["ThresholdWithoutACATConditionalOnUpdatedBH"])).sum()),
        "Mean new minus original": ots_threshold_part["NewMinusOriginalThreshold"].mean(),
        "Mean AC/AT contribution": ots_threshold_part["ACATThresholdContribution"].mean(),
        "Max AC/AT contribution": ots_threshold_part["ACATThresholdContribution"].max(),
    })
ots_threshold_summary_df = ots_pd.DataFrame(ots_threshold_summary_rows)

# %% [markdown]
# **Lower thresholds and relevant contributors.** This lists genes whose updated threshold is strictly lower and the top-three contributors in genes with original final editing or losses. The tables distinguish meaningful differences from roundoff; excluding AC/AT from the updated candidates is not automatically the historical threshold calculation.

# %%
ots_lower_thresholds_df = ots_gene_df.loc[ots_gene_df["NewEditingThreshold"].lt(ots_gene_df["NoiseThreshold"]),
    ["NoiseThreshold", "NewEditingThreshold", "NewMinusOriginalThreshold", "ThresholdsClose",
     "ThresholdWithoutACATConditionalOnUpdatedBH", "OriginalEditingSites", "LostOriginalSites"]].copy()
ots_relevant_top_three_df = ots_top_three_df.merge(
    ots_gene_df[["OriginalEditingSites", "LostOriginalSites", "NoiseThreshold", "NewEditingThreshold"]], on="Chrom", how="left", validate="many_to_one"
).loc[lambda ots_frame: ots_frame["OriginalEditingSites"].gt(0) | ots_frame["LostOriginalSites"].gt(0)].copy()
ots_display(ots_threshold_summary_df)
ots_display(ots_lower_thresholds_df.style.format(precision=17))
ots_display(ots_top_three_df.loc[ots_top_three_df["Chrom"].isin(ots_lower_thresholds_df.index)])
ots_display(ots_relevant_top_three_df.head(18))
print("Complete top-three table rows:", len(ots_top_three_df), "; relevant rows:", len(ots_relevant_top_three_df))

# %% [markdown]
# **Threshold-only recovery diagnostic.** For the full significant population, this counterfactual removes AC/AT only from threshold candidates while retaining BH, SNP and gene-filter decisions. The printed recovery count is the intersection with original final losses; it is separate from the suspected-SNP gene-count diagnostic below.

# %%
# Direct AC/AT threshold counterfactual retains all other updated decisions.
ots_hybrid_df["AboveThresholdWithoutACATDiagnostic"] = ots_hybrid_df["MismatchFrequency"].gt(
    ots_hybrid_df["Chrom"].map(ots_gene_df["ThresholdWithoutACATConditionalOnUpdatedBH"])
)
ots_acat_threshold_sites = ots_site_set(ots_hybrid_df.loc[
    ots_hybrid_df["Mismatch"].eq("A>G") & ~ots_hybrid_df["SNP"]
    & ~ots_hybrid_df["EditingDetectionDisabledInChromDueToSuspectedSNPs"]
    & ots_hybrid_df["AboveThresholdWithoutACATDiagnostic"]
])
print("Original losses recovered by removing only AC/AT threshold candidates:", len((ots_O - ots_N) & ots_acat_threshold_sites))

# %% [markdown]
# **Complete contributor table.** The input is every relevant gene's three ranked contributors, including padded zeros. The expandable table lists all positions, mismatch types, frequencies and AC/AT flags so that each gene can be inspected without relying on a truncated preview.

# %%
from IPython.display import HTML as ots_HTML

# Keep every relevant gene/rank inspectable without a many-page expanded table.
ots_contributor_columns = ["Chrom", "Rank", "Position", "Mismatch", "MismatchFrequency", "ACorAT",
                           "ZeroPadding", "TopThreeFrequencyTie", "OriginalEditingSites", "LostOriginalSites",
                           "NoiseThreshold", "NewEditingThreshold"]
ots_display(ots_HTML(
    "<details><summary>Complete top-three contributors: all genes with original editing or original losses "
    f"({ots_relevant_top_three_df['Chrom'].nunique():,} genes, {len(ots_relevant_top_three_df):,} rows)</summary>"
    '<div style="max-height:500px;overflow:auto">'
    + ots_relevant_top_three_df[ots_contributor_columns].to_html(index=False, float_format=lambda ots_value: f"{ots_value:.12g}")
    + "</div></details>"
))

# %% [markdown]
# **Check the lower-threshold cases using retained data.** For the genes with lower updated thresholds, the same top-three rule is also applied to retained original noise candidates and compared with both stored thresholds. The detailed table shows whether those contributors pass updated BH; it does not reconstruct removed pre-coverage rows or infer that missing positions were uncovered.

# %%
# Inspect the six numerically lower thresholds against retained original noise evidence.
# This is a check of available original rows, not a pre-coverage reconstruction.
ots_lower_original_noise_df = ots_audit_df.loc[
    ots_audit_df["Chrom"].isin(ots_lower_thresholds_df.index)
    & ots_audit_df["OriginalNoisyFinal"].eq(True)
    & ots_audit_df["OriginalNoise"].lt(snp_noise_level),
    ots_keys + ["OriginalRefBase", "OriginalNoise", "OriginalNoiseBinomPVal", "OriginalNoiseCorrectedPVal",
                "Mismatch", "MismatchFrequency", "BinomPVal", "BHCorrectedPVal", "BHRejection",
                "SuspectedSNP", "DefinitiveSNP", "SNP", "NoiseSite"],
].sort_values(["Chrom", "OriginalNoise", "Position"], ascending=[True, False, True]).copy()
ots_lower_original_noise_df["OriginalRetainedRank"] = ots_lower_original_noise_df.groupby("Chrom").cumcount() + 1
ots_lower_original_top_three_df = ots_lower_original_noise_df.loc[ots_lower_original_noise_df["OriginalRetainedRank"].le(3)].copy()
ots_lower_original_threshold_check_df = ots_lower_thresholds_df.copy()
ots_lower_original_threshold_check_df["ThresholdFromRetainedOriginalNoise"] = (
    ots_lower_original_top_three_df.groupby("Chrom")["OriginalNoise"].sum()
    .reindex(ots_lower_thresholds_df.index, fill_value=0).div(3).mul(1.5)
)
ots_lower_original_threshold_check_df["RetainedReconstructionMatchesOriginal"] = ots_close(
    ots_lower_original_threshold_check_df["NoiseThreshold"],
    ots_lower_original_threshold_check_df["ThresholdFromRetainedOriginalNoise"],
)
ots_lower_original_threshold_check_df["RetainedTopThreeFailUpdatedBH"] = (
    ots_lower_original_top_three_df.loc[ots_lower_original_top_three_df["BHRejection"].eq(False)]
    .groupby("Chrom").size().reindex(ots_lower_thresholds_df.index, fill_value=0)
)
ots_lower_original_threshold_check_df["RetainedReconstructionMatchesUpdated"] = ots_close(
    ots_lower_original_threshold_check_df["NewEditingThreshold"],
    ots_lower_original_threshold_check_df["ThresholdFromRetainedOriginalNoise"],
)
ots_display(ots_lower_original_threshold_check_df.style.format(precision=17))
ots_display(ots_lower_original_top_three_df)
print("Lower thresholds beyond tolerance:", int((~ots_lower_thresholds_df["ThresholdsClose"]).sum()),
      "; original editing sites in these genes:", int(ots_lower_thresholds_df["OriginalEditingSites"].sum()))

# %% [markdown]
# In the saved baseline all six decreases were beyond tolerance and none occurred in a gene with original final editing. Retained original top-three noise reproduced each updated threshold, not its stored original threshold, with the retained contributors still passing updated BH. This is consistent with original threshold calculation preceding coverage filtering; those retained tables cannot identify the missing historical contributors.

# %% [markdown]
# **AC/AT diagnostic 2: suspected-SNP gene counts only.** Starting from all updated significant `SuspectedSNP` sites, count the A>C and A>T contributions per TMR50 gene and subtract only those counts. `PushedOverLimitByACAT` identifies genes that cross the limit in this conditional count comparison; `SuspectedSNPsWithoutACAT` is an updated-count counterfactual, not a historical SNP count.

# %%
ots_acat_suspected_df = ots_hybrid_df.loc[ots_hybrid_df["SuspectedSNP"] & ots_hybrid_df["Mismatch"].isin(["A>C", "A>T"])]
ots_gene_df["ACATSuspectedSNPs"] = ots_acat_suspected_df.groupby("Chrom").size().reindex(ots_gene_universe, fill_value=0)
for ots_mismatch_type, ots_count_column in [("A>C", "ACSuspectedSNPs"), ("A>T", "ATSuspectedSNPs")]:
    ots_gene_df[ots_count_column] = ots_acat_suspected_df.loc[ots_acat_suspected_df["Mismatch"].eq(ots_mismatch_type)].groupby("Chrom").size().reindex(ots_gene_universe, fill_value=0)
ots_gene_df["SuspectedSNPsWithoutACAT"] = ots_gene_df["UpdatedSuspectedSNPs"] - ots_gene_df["ACATSuspectedSNPs"]
ots_gene_df["PushedOverLimitByACAT"] = (
    ots_gene_df["UpdatedSuspectedSNPs"].gt(max_snps_per_gene_to_allow_editing_detection)
    & ots_gene_df["SuspectedSNPsWithoutACAT"].le(max_snps_per_gene_to_allow_editing_detection)
)

# %% [markdown]
# **Detailed affected-gene table.** `ots_acat_pushed_genes_df` separates original final-editing presence from intermediate permission evidence for every pushed gene. It includes `OriginalSNPCountLowerBound`, explicitly a retained-data lower bound, beside the separate updated `SuspectedSNPsWithoutACAT` counterfactual; the copied loss subsets define the populations used in the next summary.

# %%
ots_acat_pushed_genes_df = ots_gene_df.loc[ots_gene_df["PushedOverLimitByACAT"], [
    "UpdatedSuspectedSNPs", "ACSuspectedSNPs", "ATSuspectedSNPs", "ACATSuspectedSNPs",
    "SuspectedSNPsWithoutACAT", "OriginalSNPCountLowerBound", "OriginalHasFinalEditing",
    "OriginalEditingSites", "LostOriginalSites", "OriginalEligibility",
]].copy()
ots_acat_pushed_original_genes_df = ots_acat_pushed_genes_df.loc[
    ots_acat_pushed_genes_df["OriginalHasFinalEditing"]
].copy()
ots_pushed_original_losses_df = ots_lost_df.loc[
    ots_lost_df["Chrom"].isin(ots_acat_pushed_genes_df.index)
].copy()
ots_gene_excluded_original_losses_df = ots_lost_df.loc[ots_lost_df["UpdatedGeneExcluded"]].copy()
ots_display(ots_HTML(
    "<details><summary>Detailed AC/AT-pushed genes; OriginalSNPCountLowerBound is a retained-data lower bound</summary>"
    + ots_acat_pushed_genes_df.to_html() + "</details>"
))

# %% [markdown]
# **Corrected AC/AT gene-filter summary.** `ots_gene_filter_summary_df` now labels every row with its population and unit, separating all pushed genes, pushed genes with original final editing, and pushed genes without it. BH failures and reclassification are counted independently among all original final losses, updated-gene-excluded original losses, and original losses in pushed genes, so the global and affected-gene counts cannot be conflated.

# %%
ots_gene_filter_summary_rows = []
for ots_population, ots_gene_part in [
    ("All AC/AT-pushed genes", ots_acat_pushed_genes_df),
    ("AC/AT-pushed genes with original final editing", ots_acat_pushed_original_genes_df),
    ("AC/AT-pushed genes without original final editing",
     ots_acat_pushed_genes_df.loc[~ots_acat_pushed_genes_df["OriginalHasFinalEditing"]]),
]:
    ots_gene_filter_summary_rows.append({"Population": ots_population, "Measure": "Genes", "Unit": "genes", "Count": len(ots_gene_part)})
ots_gene_filter_summary_rows.append({
    "Population": "All AC/AT-pushed genes", "Measure": "Original final editing sites", "Unit": "sites",
    "Count": int(ots_acat_pushed_genes_df["OriginalEditingSites"].sum()),
})
for ots_population, ots_loss_part in [
    ("All original final editing losses (O minus N)", ots_lost_df),
    ("Original final losses in updated-excluded genes", ots_gene_excluded_original_losses_df),
    ("Original final losses in AC/AT-pushed genes", ots_pushed_original_losses_df),
]:
    for ots_measure, ots_count in [
        ("Lost original final sites", len(ots_loss_part)),
        ("New BH failures (overlapping)", int(ots_loss_part["NewBHFailure"].sum())),
        ("Leading mismatch not A>G (overlapping)", int(ots_loss_part["LeadingMismatchNotAG"].sum())),
    ]:
        ots_gene_filter_summary_rows.append({"Population": ots_population, "Measure": ots_measure, "Unit": "sites", "Count": ots_count})
ots_gene_filter_summary_df = ots_pd.DataFrame(ots_gene_filter_summary_rows)
ots_display(ots_HTML(ots_gene_filter_summary_df.to_html(index=False)))
ots_display(ots_lost_df.loc[ots_lost_df["LeadingMismatchNotAG"]].groupby(["Mismatch", "BHRejection"], dropna=False).size().reset_index(name="Lost original final sites"))

# %% [markdown]
# **Frequency and strict-boundary diagnostics.** The audit supplies both original and updated frequencies and base counts, with frozen A>G sites distinguished from reclassified original calls. Absolute differences, direct count equality, and proximity to each threshold are measured using the stated tolerance; every actual call continues to use strict `>`.

# %%
ots_frequency_df = ots_audit_df.loc[
    ots_audit_df["OriginalEditingFrequency"].notna() & ots_audit_df["MismatchFrequency"].notna(),
    ots_keys + ["OriginalEditingSite", "EditingSite", "EditingSiteUsingOriginalThreshold", "Mismatch",
                "OriginalEditingFrequency", "MismatchFrequency", "NoiseThreshold", "NewEditingThreshold",
                "OriginalA", "OriginalT", "OriginalC", "OriginalG", "A", "T", "C", "G",
                "OriginalEditingBinomPVal", "BinomPVal", "OriginalEditingCorrectedPVal", "BHCorrectedPVal",
                "OriginalEditedCorrected", "BHRejection"],
].copy()
ots_frequency_df["FrequencyDelta"] = ots_frequency_df["MismatchFrequency"] - ots_frequency_df["OriginalEditingFrequency"]
ots_frequency_df["AbsFrequencyDelta"] = ots_frequency_df["FrequencyDelta"].abs()
ots_frequency_df["FrequenciesClose"] = ots_close(ots_frequency_df["OriginalEditingFrequency"], ots_frequency_df["MismatchFrequency"])
ots_frequency_df["CountsIdentical"] = ots_np.column_stack([
    ots_frequency_df["Original" + ots_base].eq(ots_frequency_df[ots_base]).to_numpy() for ots_base in "ATCG"
]).all(axis=1)
for ots_threshold_col, ots_suffix in [("NoiseThreshold", "Original"), ("NewEditingThreshold", "Updated")]:
    ots_frequency_df["Near" + ots_suffix + "Threshold"] = (
        ots_close(ots_frequency_df["MismatchFrequency"], ots_frequency_df[ots_threshold_col])
        | ots_close(ots_frequency_df["OriginalEditingFrequency"], ots_frequency_df[ots_threshold_col])
    )
    ots_frequency_df["StrictDecisionChangesAt" + ots_suffix + "Threshold"] = (
        ots_frequency_df["MismatchFrequency"].gt(ots_frequency_df[ots_threshold_col])
        != ots_frequency_df["OriginalEditingFrequency"].gt(ots_frequency_df[ots_threshold_col])
    ) & ots_frequency_df[ots_threshold_col].notna()

# %% [markdown]
# **Frequency comparisons in explicit site populations.** The existing summaries separately cover shared O/N sites, original final calls with frozen A>G, all comparable A>G, and reclassified original final calls. Exact inequalities and tolerance-based agreement are reported together, so roundoff is not mistaken for changed base counts.

# %%
ots_frequency_summary_rows = []
for ots_scope, ots_freq_mask in [
    ("Shared O & N (saved exact comparison)", ots_frequency_df["OriginalEditingSite"] & ots_frequency_df["EditingSite"]),
    ("Original calls, frozen A>G", ots_frequency_df["OriginalEditingSite"] & ots_frequency_df["Mismatch"].eq("A>G")),
    ("All comparable frozen A>G", ots_frequency_df["Mismatch"].eq("A>G")),
    ("Original calls reclassified", ots_frequency_df["OriginalEditingSite"] & ots_frequency_df["Mismatch"].ne("A>G")),
]:
    ots_freq_part = ots_frequency_df.loc[ots_freq_mask]
    ots_frequency_summary_rows.append({
        "Scope": ots_scope, "Sites": len(ots_freq_part),
        "Old < new (exact)": int(ots_freq_part["FrequencyDelta"].gt(0).sum()),
        "Old = new (exact)": int(ots_freq_part["FrequencyDelta"].eq(0).sum()),
        "Old > new (exact)": int(ots_freq_part["FrequencyDelta"].lt(0).sum()),
        "Equal within tolerance": int(ots_freq_part["FrequenciesClose"].sum()),
        "Counts differ": int((~ots_freq_part["CountsIdentical"]).sum()),
        "Median absolute delta": ots_freq_part["AbsFrequencyDelta"].median(),
        "99th percentile absolute delta": ots_freq_part["AbsFrequencyDelta"].quantile(.99),
        "Max absolute delta": ots_freq_part["AbsFrequencyDelta"].max(),
        "Strict decision changes at old threshold": int(ots_freq_part["StrictDecisionChangesAtOriginalThreshold"].sum()),
        "Strict decision changes at new threshold": int(ots_freq_part["StrictDecisionChangesAtUpdatedThreshold"].sum()),
    })

# %% [markdown]
# **Boundary sites and available testing-universe evidence.** `ots_boundary_sites_df` retains comparable A>G rows near a threshold or with a strict comparison change, while `ots_testing_universe_df` counts the already available theoretical and retained-original populations. The original pipeline performs separate noise/editing BH before coverage filtering, whereas the preceding theoretical section uses all ORF positions built from final retained data; missing retained rows therefore do not prove zero pre-filter coverage.

# %%
ots_frequency_summary_df = ots_pd.DataFrame(ots_frequency_summary_rows)
ots_boundary_sites_df = ots_frequency_df.loc[
    ots_frequency_df["Mismatch"].eq("A>G")
    & (ots_frequency_df["NearOriginalThreshold"] | ots_frequency_df["NearUpdatedThreshold"]
       | ots_frequency_df["StrictDecisionChangesAtOriginalThreshold"]
       | ots_frequency_df["StrictDecisionChangesAtUpdatedThreshold"])
].copy()
ots_recovered_boundary_sites = ots_site_set(ots_boundary_sites_df) & ots_site_set(ots_recovered_df)
ots_testing_universe_df = ots_pd.DataFrame([
    {"Population": "Full theoretical BH universe", "Positions": len(mismatches_theoretical_fix_df)},
    {"Population": "Theoretical BH rejections", "Positions": len(ots_hybrid_df)},
    {"Population": "Theoretical rows absent from retained original positions (coverage unknown)",
     "Positions": int((ots_audit_df["PresentInTheoreticalUniverse"] & ~ots_audit_df["PresentInOriginalPositions"]).sum())},
    {"Population": "Original editing-test p-values available after coverage filtering (not full BH denominator)",
     "Positions": int(ots_audit_df["OriginalEditingBinomPVal"].notna().sum())},
    {"Population": "Original noise-test p-values available after coverage filtering (not full BH denominator)",
     "Positions": int(ots_audit_df["OriginalNoiseBinomPVal"].notna().sum())},
])

# %% [markdown]
# **Reading the frequency and BH tables.** These outputs summarize the same frozen-frequency diagnostics and existing original versus updated BH annotations. `OriginalEditedCorrected` remains an intermediate significance flag throughout; no BH rerun, missing-count reconstruction, or new coverage investigation is performed.

# %%
ots_display(ots_frequency_summary_df.style.format(precision=17))
print("Smallest single-read fraction among comparable A>G denominators:",
      1.0 / (ots_frequency_df.loc[ots_frequency_df["Mismatch"].eq("A>G"), "A"]
             + ots_frequency_df.loc[ots_frequency_df["Mismatch"].eq("A>G"), "G"]).max())
print("Boundary-sensitive comparable sites:", len(ots_boundary_sites_df), "; recovered originals among them:", len(ots_recovered_boundary_sites))
ots_display(ots_boundary_sites_df.head(12))
ots_display(ots_testing_universe_df)
ots_bh_transition_df = ots_frequency_df.groupby(
    ["OriginalEditingSite", "Mismatch", "OriginalEditedCorrected", "BHRejection"], dropna=False
).size().reset_index(name="Sites with both frequencies")
ots_display(ots_bh_transition_df)
ots_display(ots_remaining_df.loc[ots_remaining_df["NewBHFailure"], ots_keys + [
    "OriginalEditingBinomPVal", "BinomPVal", "OriginalEditingCorrectedPVal", "BHCorrectedPVal",
    "OriginalEditedCorrected", "BHRejection", "Mismatch", "OriginalEditingFrequency", "MismatchFrequency",
]].head(12))

# %% [markdown]
# **Mismatch-distribution count tables.** All three plots use the full updated significant mismatch table, including non-editing sites and genes excluded from final editing detection; they are not filtered to original `EditedFinal`. Count tables retain zero categories and unchanged frequency windows, assert plot 3's threshold independence, and report overlap between plots 1 and 3 when thresholds reach/exceed the suspected-SNP level.

# %%
ots_plot_masks = {
    "Updated threshold": [~ots_hybrid_df["AboveNewEditingThreshold"],
                          ots_hybrid_df["AboveNewEditingThreshold"] & ~ots_hybrid_df["AtOrAboveSuspectedSNPLevel"],
                          ots_hybrid_df["AtOrAboveSuspectedSNPLevel"]],
    "Original threshold": [~ots_hybrid_df["AboveOriginalEditingThreshold"],
                           ots_hybrid_df["AboveOriginalEditingThreshold"] & ~ots_hybrid_df["AtOrAboveSuspectedSNPLevel"],
                           ots_hybrid_df["AtOrAboveSuspectedSNPLevel"]],
}
ots_plot_count_parts = []
for ots_scheme, ots_masks in ots_plot_masks.items():
    for ots_plot_number, ots_plot_mask in enumerate(ots_masks, 1):
        if ots_plot_number == 3:
            ots_count_index = ots_pd.MultiIndex.from_product([mismatches, [False, True]], names=["Mismatch", "MismatchFrequency1"])
            ots_plot_counts = ots_hybrid_df.loc[ots_plot_mask].groupby(["Mismatch", "MismatchFrequency1"]).size().reindex(ots_count_index, fill_value=0).reset_index(name="Count")
        else:
            ots_plot_counts = ots_hybrid_df.loc[ots_plot_mask, "Mismatch"].value_counts().reindex(mismatches, fill_value=0).rename_axis("Mismatch").reset_index(name="Count")
            ots_plot_counts["MismatchFrequency1"] = "all"
        ots_plot_counts["Scheme"], ots_plot_counts["Plot"] = ots_scheme, ots_plot_number
        assert ots_plot_counts["Count"].sum() == ots_plot_mask.sum()
        ots_plot_count_parts.append(ots_plot_counts)
ots_plot_counts_df = ots_pd.concat(ots_plot_count_parts, ignore_index=True)
ots_plot_comparison_df = ots_plot_counts_df.pivot(index=["Plot", "Mismatch", "MismatchFrequency1"], columns="Scheme", values="Count").reset_index()
ots_plot_comparison_df["Original minus updated"] = ots_plot_comparison_df["Original threshold"] - ots_plot_comparison_df["Updated threshold"]
assert ots_plot_masks["Original threshold"][2].equals(ots_plot_masks["Updated threshold"][2])
assert ots_plot_comparison_df.loc[ots_plot_comparison_df["Plot"].eq(3), "Original minus updated"].eq(0).all()
ots_plot2_ag_summary_df = ots_plot_comparison_df.loc[ots_plot_comparison_df["Plot"].eq(2)].assign(
    MismatchGroup=lambda ots_frame: ots_np.where(ots_frame["Mismatch"].eq("A>G"), "A>G", "non-A>G")
).groupby("MismatchGroup")[["Updated threshold", "Original threshold", "Original minus updated"]].sum()
ots_plot_overlap_df = ots_pd.DataFrame([
    {"Scheme": ots_scheme, "Threshold > suspected-SNP level genes": int(ots_thresholds_df[ots_col].gt(snp_noise_level).sum()),
     "Threshold = suspected-SNP level genes": int(ots_thresholds_df[ots_col].eq(snp_noise_level).sum()),
     "Sites shared by plots 1 and 3": int((ots_plot_masks[ots_scheme][0] & ots_plot_masks[ots_scheme][2]).sum())}
    for ots_scheme, ots_col in [("Updated threshold", "NewEditingThreshold"), ("Original threshold", "NoiseThreshold")]
])
ots_display(ots_plot2_ag_summary_df)
ots_display(ots_plot_overlap_df)

# %% [markdown]
# **Plot 1: below the original threshold.** The population is the complement of strict `AboveOriginalEditingThreshold`, so equality is included without an extra SNP-frequency restriction. The count table and Plotly histogram retain the existing colors, mismatch order, 700×500 dimensions, and log-y scale.

# %%
ots_display(ots_plot_comparison_df.loc[ots_plot_comparison_df["Plot"].eq(1)])
ots_fig_below_original = ots_px.histogram(
    ots_hybrid_df.loc[ots_plot_masks["Original threshold"][0]],
    x="Mismatch", color="Mismatch", color_discrete_map=mismatch_dolor_map,
    log_y=True, template=template, category_orders={"Mismatch": mismatches},
    title="Mismatches below editing threshold<br><sup>Original-threshold sensitivity analysis</sup>",
)
ots_fig_below_original.for_each_annotation(lambda ots_annotation: ots_annotation.update(text=ots_annotation.text.split("=")[-1]))
ots_fig_below_original.update_xaxes(tickangle=35)
ots_fig_below_original.update_layout(width=700, height=500, showlegend=False)
ots_fig_below_original.show()

# %% [markdown]
# **Plot 2: above the original threshold and below the suspected-SNP level.** The population uses `AboveOriginalEditingThreshold & ~AtOrAboveSuspectedSNPLevel` on all significant mismatches. Its count table and histogram retain the exact previous window and formatting, including non-final mismatches and genes excluded from detection.

# %%
ots_display(ots_plot_comparison_df.loc[ots_plot_comparison_df["Plot"].eq(2)])
ots_fig_above_original = ots_px.histogram(
    ots_hybrid_df.loc[ots_plot_masks["Original threshold"][1]],
    x="Mismatch", color="Mismatch", color_discrete_map=mismatch_dolor_map,
    log_y=True, template=template, category_orders={"Mismatch": mismatches},
    title="Mismatches above editing threshold, but below suspected SNP threshold<br><sup>Original-threshold sensitivity analysis</sup>",
)
ots_fig_above_original.for_each_annotation(lambda ots_annotation: ots_annotation.update(text=ots_annotation.text.split("=")[-1]))
ots_fig_above_original.update_xaxes(tickangle=35)
ots_fig_above_original.update_layout(width=700, height=500, showlegend=False)
ots_fig_above_original.show()

# %% [markdown]
# **Plot 3: at or above the suspected-SNP level.** This uses `AtOrAboveSuspectedSNPLevel`, faceted by `MismatchFrequency1`, independently of either editing threshold. The 1000×500 log-y plot and its counts remain unchanged, and its population may overlap plot 1 rather than partitioning all sites into exclusive windows.

# %%
ots_display(ots_plot_comparison_df.loc[ots_plot_comparison_df["Plot"].eq(3)])
ots_fig_snp_original = ots_px.histogram(
    ots_hybrid_df.loc[ots_plot_masks["Original threshold"][2]],
    x="Mismatch", color="Mismatch", color_discrete_map=mismatch_dolor_map,
    facet_col_spacing=0.04, facet_col="MismatchFrequency1",
    log_y=True, template=template, category_orders={"Mismatch": mismatches},
    title="Mismatches at or above suspected SNP threshold<br><sup>Original-threshold sensitivity analysis</sup>",
)
ots_fig_snp_original.update_xaxes(tickangle=35)
ots_fig_snp_original.update_layout(width=1000, height=500, showlegend=False)
ots_fig_snp_original.show()

# %% [markdown]
# **Regression checks for the reporting revision.** The frozen significant columns are checked against the input, and loss combinations must account for every lost final site. `ots_regression_df` compares measured global counts with the saved expectations, while `ots_affected_gene_regression_df` checks the explicitly separated gene populations; differences are reported and must be traced rather than removed by tuning.

# %%
ots_pd.testing.assert_frame_equal(
    ots_hybrid_df.loc[:, significant_mismatches_theoretical_fix_df.columns],
    significant_mismatches_theoretical_fix_df, check_exact=True,
)
assert ots_lost_combinations_df["Sites"].sum() == len(ots_O - ots_N)
assert ots_remaining_combinations_df["Sites"].sum() == len(ots_O - ots_H)
ots_regression_df = ots_pd.DataFrame({
    "Global metric": ["O", "N", "H", "O intersection N", "O intersection H", "O minus N", "Original losses in updated-excluded genes"],
    "Saved expectation": [11711, 10714, 10816, 10403, 10504, 1308, 1205],
    "Measured": [len(ots_O), len(ots_N), len(ots_H), len(ots_O & ots_N), len(ots_O & ots_H), len(ots_O - ots_N), len(ots_gene_excluded_original_losses_df)],
})
ots_regression_df["Difference"] = ots_regression_df["Measured"] - ots_regression_df["Saved expectation"]
ots_affected_gene_regression_df = ots_pd.DataFrame({
    "Population / measure": ["All pushed genes", "Pushed genes with original final editing", "Pushed genes without original final editing", "Original final losses in pushed genes", "Reclassified among all original final losses", "Reclassified within updated-gene-excluded original losses"],
    "Saved expectation": [111, 103, 8, 1205, 5, 3],
    "Measured": [len(ots_acat_pushed_genes_df), len(ots_acat_pushed_original_genes_df), int((~ots_acat_pushed_genes_df["OriginalHasFinalEditing"]).sum()), len(ots_pushed_original_losses_df), int(ots_lost_df["LeadingMismatchNotAG"].sum()), int(ots_gene_excluded_original_losses_df["LeadingMismatchNotAG"].sum())],
})
ots_affected_gene_regression_df["Difference"] = ots_affected_gene_regression_df["Measured"] - ots_affected_gene_regression_df["Saved expectation"]
ots_display(ots_regression_df)
ots_display(ots_affected_gene_regression_df)
if ots_regression_df["Difference"].ne(0).any():
    print("Global regression differs: inspect changed site keys and upstream state before interpreting this reporting revision.")

# %% [markdown]
# **Before/after report.** The report below states whether final-site membership changed and distinguishes reporting-scope corrections from detection rules. It uses measured global and restricted-scope results and points to the corrected population summary and detailed gene table.
#
# File verification for this revision: all 719 cells outside this subsection retain their latest source, IDs, metadata, execution counts, and saved outputs exactly, including concurrent edits in the preceding analysis. All 37 revised code cells executed in the prepared kernel; original O/N/H site keys, existing comparison scopes, all plot count tables, and the three complete Plotly figure specifications match the captured baseline.

# %%
ots_old_retention = ots_ratio(len(ots_O & ots_N), len(ots_O))
ots_hybrid_retention = ots_ratio(len(ots_O & ots_H), len(ots_O))
ots_primary_pairs_df = ots_pairwise_df.loc[ots_pairwise_df["Scope"].eq(ots_primary_scope)].copy()
ots_primary_oh = ots_primary_pairs_df.loc[ots_primary_pairs_df["First"].eq("O") & ots_primary_pairs_df["Second"].eq("H")].iloc[0]
ots_primary_on = ots_primary_pairs_df.loc[ots_primary_pairs_df["First"].eq("O") & ots_primary_pairs_df["Second"].eq("N")].iloc[0]
ots_secondary_counts = ots_gene_counts_df.loc[ots_gene_counts_df["Scope"].eq("Established eligible in both")].set_index("Scheme")["Editing sites"]
ots_conclusion = (
    "**Before/after.** No original final editing-site count was based on intermediate `Edited`: "
    "O and the original audit flag already used `EditedFinal`. Direct key assertions now make that explicit. "
    "The logic fix is to label and consistently calculate AC/AT summary populations; the added primary restriction "
    "and the separation of final editing from permission evidence are reporting-scope corrections and clarifications, "
    "not detection-algorithm changes.\n\n"
    f"{int(ots_regression_df['Difference'].ne(0).sum())} global regression counts differ from the saved baseline: "
    f"O={len(ots_O):,}, N={len(ots_N):,}, H={len(ots_H):,}. Retention remains "
    f"{ots_old_retention:.2%} for N and {ots_hybrid_retention:.2%} for H, with "
    f"{len(ots_recovered_df):,} original losses recovered. The primary restriction contains "
    f"{len(ots_primary_genes):,} of {int(ots_gene_df['OriginalHasFinalEditing'].sum()):,} original-final-editing genes; "
    f"its O/N/H site counts are {int(ots_primary_on['First sites']):,}/"
    f"{int(ots_primary_on['Second sites']):,}/{int(ots_primary_oh['Second sites']):,}. "
    f"The previous broader permission-evidence scope remains secondary ({len(ots_common_genes):,} genes, "
    f"including {len(ots_common_genes - ots_primary_genes):,} without original final editing), with O/N/H counts "
    f"{ots_secondary_counts['O']:,}/{ots_secondary_counts['N']:,}/{ots_secondary_counts['H']:,}.\n\n"
    f"Of {len(ots_acat_pushed_genes_df):,} AC/AT-pushed genes, "
    f"{len(ots_acat_pushed_original_genes_df):,} have original final editing and "
    f"{int((~ots_acat_pushed_genes_df['OriginalHasFinalEditing']).sum()):,} do not; "
    f"{len(ots_pushed_original_losses_df):,} original final sites are lost in those genes. "
    f"Reclassification affects {int(ots_lost_df['LeadingMismatchNotAG'].sum())} sites among all original losses, "
    f"versus {int(ots_gene_excluded_original_losses_df['LeadingMismatchNotAG'].sum())} within gene-excluded losses. "
    "See `ots_pairwise_df` for corrected scope denominators, `ots_gene_filter_summary_df` for the population-labelled "
    "summary, and `ots_acat_pushed_genes_df` for every affected gene and its retained-data SNP lower bound. "
    "The mismatch plot populations and strict threshold windows are unchanged."
)
ots_display(ots_Markdown(ots_conclusion))

# %% [markdown] papermill={"duration": 0.02598, "end_time": "2022-02-01T09:42:46.438342", "exception": false, "start_time": "2022-02-01T09:42:46.412362", "status": "completed"}
# ## Reads

# %% [markdown] papermill={"duration": 0.02598, "end_time": "2022-02-01T09:42:46.438342", "exception": false, "start_time": "2022-02-01T09:42:46.412362", "status": "completed"}
# ### All

# %% [markdown]
# That is, all filtered reads.

# %% papermill={"duration": 1.204258, "end_time": "2022-02-01T09:42:47.668206", "exception": false, "start_time": "2022-02-01T09:42:46.463948", "status": "completed"}
reads_dfs = [
    pd.read_csv(reads_file, sep=sep, dtype={"Read": str}) for reads_file in reads_files
]
for chrom, reads_df in zip(chroms, reads_dfs):
    reads_df.insert(0, "Chrom", chrom)
    
reads_first_col_pos += 1
reads_dfs[0]


# %%

# %%
# reads_dfs[1]

# %%
import subprocess


def count_lines_in_file(file, cat_cmd="zcat"):
    cmd = f"{cat_cmd} {file} | wc -l"
    clean_str_output = (
        subprocess.run(cmd, shell=True, capture_output=True)
        .stdout.decode()
        .removesuffix("\n")
    )
    return int(clean_str_output)


# %%
num_of_reads_in_reads_files = [
    count_lines_in_file(reads_file) - 1 for reads_file in reads_files
]
num_of_reads_in_reads_files[:3]

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


# %% [markdown]
# ### Reads with full haplotype information (SNPs and editing positions covered)

# %%
def make_editing_and_snps_reads_df(
    reads_snps_file,
    reads_df,
    reads_first_col_pos,
    sep
):
    # this df contain snps status per position
    reads_snps_df = pd.read_csv(
        reads_snps_file, 
        sep=sep,
        dtype={"Read": str}
    ).rename(
        columns={
            "Sample": "Chrom"
        }
    ).drop(
        columns="Platform"
    )
    
    if reads_snps_df.empty:
        return None
    
    # add all reads (even those without a single SNP covered) to reads_snps_df
    # to allow direct insertion of data from reads_snps_df into editing_and_snps_reads_df (see below)
    # based on ordering of the two dfs by the Read col
    reads_snps_df = reads_snps_df.merge(
        reads_df.loc[:, ["Chrom", "Read"]],
        how="right"
    )
    reads_snps_df = reads_snps_df.sort_values("Read", ignore_index=True)

    # fill missing SNPs information with -1 and change dtype to int
    reads_snps_df.iloc[:, 2:] = reads_snps_df.iloc[:, 2:].fillna(-1)
    for col in reads_snps_df.columns[2:]:
        reads_snps_df[col] = reads_snps_df[col].astype(int)

    snps_positions = reads_snps_df.columns[2:]
    # snps_positions
    
    editing_and_snps_reads_df = (
        reads_df
        .merge(
            reads_snps_df.loc[:, ["Chrom", "Read"]],
            how="left"
        )
        .sort_values("Read", ignore_index=True)
    )

    # insert each indivdual SNP to reads with editing status df
    reads_first_col_pos_in_editing_and_snps_reads_df = reads_first_col_pos
    # ic(reads_first_col_pos_in_editing_and_snps_reads_df)

    for snp_position in snps_positions:
        # ic(snp_position)
        # we can directly insert this col from reads_snps_df into editing_and_snps_reads_df 
        # because both dataframes are sorted by Read
        editing_and_snps_reads_df.insert(
            reads_first_col_pos_in_editing_and_snps_reads_df,
            f"SNP_{snp_position}",
            reads_snps_df[snp_position]
        )
        reads_first_col_pos_in_editing_and_snps_reads_df += 1
        # ic(reads_first_col_pos_in_editing_and_snps_reads_df);


    # Identify the editing site columns (starting from the current offset onwards)
    editing_cols = editing_and_snps_reads_df.columns[reads_first_col_pos_in_editing_and_snps_reads_df:]
    # editing_cols

    snp_cols = [f"SNP_{pos}" for pos in snps_positions]
    # snp_cols

    # Filter to keep only reads where the SNP positions have a valid base (0 or 1, not -1)
    # Using .all(axis=1) ensures we have a complete haplotype for the selected SNPs
    valid_reads_mask = editing_and_snps_reads_df[snp_cols].ne(-1).all(axis=1)
    editing_and_snps_reads_df = editing_and_snps_reads_df.loc[valid_reads_mask].copy()
    
    # if no reads hold complete haplotype information, terminate early
    if editing_and_snps_reads_df.empty:
        return None

    # Replace -1 (unmapped/ignored) with NaN in the editing columns
    # so means are computed on covered reads only
    editing_and_snps_reads_df[editing_cols] = editing_and_snps_reads_df[editing_cols].replace(-1, np.nan)

    # define a single haplotype label (useful for counting / plotting / downstream merges)
    # e.g. "SNP_123=0|SNP_456=1|SNP_789=0"
    editing_and_snps_reads_df.insert(
        reads_first_col_pos_in_editing_and_snps_reads_df,
        "Haplotype",
        editing_and_snps_reads_df[snp_cols].astype(str).agg("|".join, axis=1)
    )
    reads_first_col_pos_in_editing_and_snps_reads_df += 1
    # ic(reads_first_col_pos_in_editing_and_snps_reads_df)

    # return editing_and_snps_reads_df, reads_first_col_pos_in_editing_and_snps_reads_df
    return editing_and_snps_reads_df

# %%
editing_and_snps_reads_dfs = [
        make_editing_and_snps_reads_df(
        reads_snps_file,
        reads_df,
        reads_first_col_pos,
        sep
    )
    for reads_snps_file, reads_df in zip(
        snps_reads_files, reads_dfs,
        strict=True
    )
]

# keep only results where the editing and snps reads df is not None (i.e. there were reads with complete haplotype information)
editing_and_snps_reads_dfs = [
    res
    for res in editing_and_snps_reads_dfs
    if type(res) == pd.core.frame.DataFrame
]

ic(len(editing_and_snps_reads_dfs));

# %%
editing_and_snps_reads_dfs[0]

# %%

# %% [markdown] papermill={"duration": 0.041741, "end_time": "2022-02-01T09:42:47.760215", "exception": false, "start_time": "2022-02-01T09:42:47.718474", "status": "completed"}
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
        distinct_unique_proteins_file, 
        sep=sep,
        dtype={"UniqueSamples": str, "AvailableReads": str, "MissingUniqueSamples": str}
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

if "NumMissingUniqueSamples" in distinct_unique_proteins_df.columns:
    distinct_unique_proteins_df["NumMissingUniqueSamples"] = distinct_unique_proteins_df["NumMissingUniqueSamples"].fillna(0)
    assert distinct_unique_proteins_df["NumMissingUniqueSamples"].eq(0).all()
    del distinct_unique_proteins_df["MissingUniqueSamples"]
    del distinct_unique_proteins_df["NumMissingUniqueSamples"]

distinct_unique_proteins_df


# %%
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
unique_proteins_dfs[0]

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
assert (
    len(tmr1000_conditions)
    == len(tmr1000_chroms)
    == len(tmr1000_distinct_unique_proteins_files)
    == len(tmr1000_unique_reads_dfs)
)

tmr1000_distinct_unique_proteins_dfs = []
for condition, chrom, distinct_unique_proteins_file, unique_reads_df in zip(
    tmr1000_conditions,
    tmr1000_chroms,
    # tmr1000_distinct_proteins_files,
    tmr1000_distinct_unique_proteins_files, # 7.9.2026 update
    tmr1000_unique_reads_dfs,
):
    tmr1000_distinct_unique_proteins_df = pd.read_csv(
        distinct_unique_proteins_file, 
        sep=sep,
        dtype={"UniqueSamples": str, "AvailableReads": str}
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
# ## Expression

# %%
expression_dfs = []
for chrom, expression_file in zip(chroms, expression_files):
    # expression_df = pd.read_csv(expression_file, sep=sep)
    # expression_df["#Solution"] = expression_df["#Solution"].astype(str)
    expression_df = pd.read_csv(
        expression_file,
        sep=sep,
        dtype={
            "#Solution": str,
            "AdditionalSupportingReadsIDs": str,
            "AdditionalSupportingProteinsIDs": str,
        },
    )
    expression_df.insert(0, "Chrom", chrom)

    expression_df["AdditionalSupportingReadsIDs"] = expression_df[
        "AdditionalSupportingReadsIDs"
    ].apply(lambda x: "" if pd.isna(x) else [y.split(",") for y in x.split(";")])
    expression_df["AdditionalSupportingProteinsIDs"] = expression_df[
        "AdditionalSupportingProteinsIDs"
    ].apply(lambda x: "" if pd.isna(x) else x.split(","))

    expression_dfs.append(expression_df)
expression_dfs[0]

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
concat_per_transcript_per_sample_coverage_df = pd.concat(
    per_transcript_per_sample_coverage_dfs
)
concat_per_transcript_per_sample_coverage_df

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
def calc_per_sample_editing_index_df(
    concat_all_positions_df, chroms, strands, samples, processes=1
):
    """Calc A-to-I RNA editing index per sample, considering edited genes."""
    
    # keep only chroms in which we detected editing (at least 1 editing site)
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

# %%
per_sample_editing_index_df.round(2)

# %% [markdown]
# ### Editing - per-site per-sample

# %% [markdown]
# For each site, how many samples support the MM? Maybe cut-off by coverage per sample [50].

# %%
concat_all_positions_df

# %%
concat_edited_positions_df = concat_all_positions_df.loc[
    concat_all_positions_df["EditedFinal"]
]
concat_edited_positions_df

# %%
concat_edited_positions_df["Chrom"].nunique()

# %%
concat_edited_positions_df.groupby("Chrom").size().describe().round(2)

# %%
expanded_concat_edited_positions_df = (
        concat_edited_positions_df
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
                # "TotalCoverage",
                "Edited",
                "CDS",
                "KnownEditing",
                "InProbRegion",
                "RefBase",
                "NoiseBinomPVal",
                "NoiseCorrectedPVal",
                "NoisyCorrected",
                "BelowNoiseFreq1",
                "NoisyFinal",
                "EditingBinomPVal",
                "EditingCorrectedPVal",
                "EditedCorrected",
                "EditedFinal",
                "BelowEditingFreq1"
            ],
            axis=1,
        )
    )

expanded_concat_edited_positions_df["Samples"] = expanded_concat_edited_positions_df[
    "Samples"
].str.split(",")
expanded_concat_edited_positions_df[
    "MappedBases"
] = expanded_concat_edited_positions_df["MappedBases"].apply(list)

# now is the time the df is really expanded
expanded_concat_edited_positions_df = expanded_concat_edited_positions_df.explode(
    ["Samples", "MappedBases"]
).rename(
    columns={"Samples": "Sample", "MappedBases": "MappedBase"}
)

expanded_concat_edited_positions_df["MappedBase"] = expanded_concat_edited_positions_df["MappedBase"].apply(
    lambda x: "A" if x == "."
    else "G" if x == "G"
    else np.nan
)

expanded_concat_edited_positions_df = expanded_concat_edited_positions_df.loc[
    expanded_concat_edited_positions_df["MappedBase"].notna()
]

expanded_concat_edited_positions_df = expanded_concat_edited_positions_df.groupby(
    ["Chrom", "Transcript", "Sample", "Position"]
)["MappedBase"].value_counts(dropna=False).reset_index()

expanded_concat_edited_positions_df

# %%
per_sample_agged_expanded_concat_edited_positions_df = expanded_concat_edited_positions_df.pivot(
    index=["Chrom", "Transcript", "Sample", "Position"],
    columns="MappedBase",
    values="count",
    # fill_value=0,
).reset_index().fillna(0)

# per_sample_agged_expanded_concat_edited_positions_df["TotalCoverage"] = per_sample_agged_expanded_concat_edited_positions_df.loc[:, ["A", "G", np.nan]].sum(axis=1)
per_sample_agged_expanded_concat_edited_positions_df["TotalCoverage"] = (
    per_sample_agged_expanded_concat_edited_positions_df.loc[:, ["A", "G"]].sum(axis=1)
)
per_sample_agged_expanded_concat_edited_positions_df["EditingFrequency"] = (
    per_sample_agged_expanded_concat_edited_positions_df["G"] / per_sample_agged_expanded_concat_edited_positions_df["TotalCoverage"]
)

per_sample_agged_expanded_concat_edited_positions_df

# %%
per_sample_agged_expanded_concat_edited_positions_df.groupby(
    ["Chrom", "Position"]
)["TotalCoverage"].mean().describe()

# %%
fig = px.histogram(
    per_sample_agged_expanded_concat_edited_positions_df.groupby(
        ["Chrom", "Position"]
    )["TotalCoverage"].mean().reset_index(name="MeanTotalCoveragePerSampleInEditingSite"),
    x="MeanTotalCoveragePerSampleInEditingSite",
    # log_x=True,
    # log_y=True,
    histnorm="percent",
    # cumulative=True,
)
fig.update_layout(
    template=template,
    width=600,
    height=400
)
fig.show()

# %%
min_tot_covs = [0, 10, 50, 250]
min_editing_freqs = [0.001, 0.003, 0.005, 0.01]

# %%
num_of_sites_per_sample_per_min_cov_and_editing_freq_dfs = []
for min_tot_cov, min_editing_freq in product(min_tot_covs, min_editing_freqs):
    # min_tot_cov = min_tot_covs[1]
    # min_editing_freq = min_editing_freqs[0]
    print(f"Min total coverage: {min_tot_cov}, min editing frequency: {min_editing_freq}")
    num_of_sites_per_sample_per_min_cov_and_editing_freq_df = per_sample_agged_expanded_concat_edited_positions_df.loc[
        (per_sample_agged_expanded_concat_edited_positions_df["TotalCoverage"].ge(min_tot_cov))
        & (per_sample_agged_expanded_concat_edited_positions_df["EditingFrequency"].ge(min_editing_freq))
    ].groupby(
        ["Chrom", "Transcript", "Position"]
    )["Sample"].apply(list).reset_index(name="Samples")
    num_of_sites_per_sample_per_min_cov_and_editing_freq_df.insert(0, "MinTotalCoverage", min_tot_cov)
    num_of_sites_per_sample_per_min_cov_and_editing_freq_df.insert(1, "MinEditingFrequency", min_editing_freq)
    num_of_sites_per_sample_per_min_cov_and_editing_freq_df = num_of_sites_per_sample_per_min_cov_and_editing_freq_df.explode("Samples").rename(
        columns={"Samples": "Sample"}
    ).groupby(
        #    ["MinTotalCoverage", "MinEditingFrequency", "NumOfSitesCoveredByAtLeastOneSample", "Sample",]
        ["MinTotalCoverage", "MinEditingFrequency", "Sample",]
    ).size().reset_index(name="NumOfEditingSitesPerSample")
    # num_of_all_sites_covered_by_at_least_one_sample = per_sample_agged_expanded_concat_edited_positions_df.loc[
    #     (per_sample_agged_expanded_concat_edited_positions_df["TotalCoverage"].ge(min_tot_cov))
    # ].drop_duplicates(
    #     ["Chrom", "Transcript", "Position"]
    # ).shape[0]
    # num_of_sites_per_sample_per_min_cov_and_editing_freq_df["NumOfSitesCoveredByAtLeastOneSample"] = num_of_all_sites_covered_by_at_least_one_sample
    
    sites_covered_per_sample_df = per_sample_agged_expanded_concat_edited_positions_df.loc[
        (per_sample_agged_expanded_concat_edited_positions_df["TotalCoverage"].ge(min_tot_cov))
    ].groupby("Sample").size().reset_index(name="NumOfSitesCoveredPerSample")

    num_of_sites_per_sample_per_min_cov_and_editing_freq_df = num_of_sites_per_sample_per_min_cov_and_editing_freq_df.merge(
        sites_covered_per_sample_df, on="Sample"
    )
    num_of_sites_per_sample_per_min_cov_and_editing_freq_df["%EditedOfCoveredSites"] = (
        100
        * num_of_sites_per_sample_per_min_cov_and_editing_freq_df["NumOfEditingSitesPerSample"]
        / num_of_sites_per_sample_per_min_cov_and_editing_freq_df["NumOfSitesCoveredPerSample"]
    )
    
    num_of_sites_per_sample_per_min_cov_and_editing_freq_dfs.append(num_of_sites_per_sample_per_min_cov_and_editing_freq_df)
    
    # break
    
concat_num_of_sites_per_sample_per_min_cov_and_editing_freq_df = pd.concat(
    num_of_sites_per_sample_per_min_cov_and_editing_freq_dfs, ignore_index=True
)

# concat_num_of_sites_per_sample_per_min_cov_and_editing_freq_df["NumOfSamples"] = concat_num_of_sites_per_sample_per_min_cov_and_editing_freq_df["Samples"].apply(len)

# # concat_num_of_sites_per_sample_per_min_cov_and_editing_freq_df = concat_num_of_sites_per_sample_per_min_cov_and_editing_freq_df.rename(
# #     columns={
# #         "count": "NumOfEditingSites",
# #         "Sample": "NumOfSamples"
# #     }
# # )


# concat_num_of_sites_per_sample_per_min_cov_and_editing_freq_df = concat_num_of_sites_per_sample_per_min_cov_and_editing_freq_df.explode("Samples").rename(
#     columns={"Samples": "Sample"}
# ).groupby(
# #    ["MinTotalCoverage", "MinEditingFrequency", "NumOfSitesCoveredByAtLeastOneSample", "Sample",]
#    ["MinTotalCoverage", "MinEditingFrequency", "Sample",]
# ).size().reset_index(name="NumOfEditingSitesPerSample")

# concat_num_of_sites_per_sample_per_min_cov_and_editing_freq_df["%EditedSitesOfAllSitesCoveredByAtLeastOneSample"] = (
#     100 
#     * concat_num_of_sites_per_sample_per_min_cov_and_editing_freq_df["NumOfEditingSitesPerSample"] 
#     / concat_num_of_sites_per_sample_per_min_cov_and_editing_freq_df["NumOfSitesCoveredByAtLeastOneSample"]
# )

concat_num_of_sites_per_sample_per_min_cov_and_editing_freq_df

# %%

# %%
per_sample_agged_expanded_concat_edited_positions_df

# %%

# %%
# # min_tot_cov = min_tot_covs[0]
# # min_editing_freq = min_editing_freqs[0]

# min_tot_cov = min_tot_covs[1]
# min_editing_freq = min_editing_freqs[1]

# %%
# per_sample_agged_expanded_concat_edited_positions_df

# %%
# # create a df with 5 cols (MinTotalCoverage, MinEditingFrequency, condition_col, Position, Replicates)
# # denoting for each gene-position the replicates that have at least min_tot_cov coverage and at least min_editing_freq editing frequency
# replicates_per_gene_per_position_per_min_cov_and_editing_freq_df = (
#     per_sample_agged_expanded_concat_edited_positions_df
#     .loc[
#         (per_sample_agged_expanded_concat_edited_positions_df["TotalCoverage"].ge(min_tot_cov))
#         & (per_sample_agged_expanded_concat_edited_positions_df["EditingFrequency"].ge(min_editing_freq))
#     ]
#     .groupby(["Chrom", condition_col, "Position"])
#     ["Sample"].apply(list)
#     .reset_index(name="Samples")
# )
# replicates_per_gene_per_position_per_min_cov_and_editing_freq_df["Samples"] = (
#         replicates_per_gene_per_position_per_min_cov_and_editing_freq_df["Samples"].apply(
#         lambda x: x + ["All"] if len(x) == 7 else x
#     )
# )
# replicates_per_gene_per_position_per_min_cov_and_editing_freq_df.insert(0, "MinTotalCoverage", min_tot_cov)
# replicates_per_gene_per_position_per_min_cov_and_editing_freq_df.insert(1, "MinEditingFrequency", min_editing_freq)
# replicates_per_gene_per_position_per_min_cov_and_editing_freq_df

# %%
# replicates_per_gene_per_position_per_min_cov_and_editing_freq_df["Chrom"].nunique()

# %%
# replicates_per_gene_per_position_per_min_cov_and_editing_freq_df["Samples"].apply(
#     len
# ).sum()

# %%
# # create a df with 6 cols (MinTotalCoverage, MinEditingFrequency, condition_col, Replicate, NumOfEditingSitesPerReplicate)
# # denoting for each gene-replicate the number of editing sites that have at least min_tot_cov coverage and at least min_editing_freq editing frequency
# num_of_sites_per_gene_per_sample_per_min_cov_and_editing_freq_df = (
#     replicates_per_gene_per_position_per_min_cov_and_editing_freq_df
#     .explode("Samples")
#     .rename(
#         columns={"Samples": "Sample"}
#     )
#     .groupby(
#         ["MinTotalCoverage", "MinEditingFrequency", "Chrom", condition_col, "Sample",]
#     )
#     .size()
#     .reset_index(name="NumOfEditingSitesPerSample")
# )
# num_of_sites_per_gene_per_sample_per_min_cov_and_editing_freq_df

# %%
# num_of_sites_per_gene_per_sample_per_min_cov_and_editing_freq_df["Chrom"].nunique()

# %%
# num_of_sites_per_gene_per_sample_per_min_cov_and_editing_freq_df["NumOfEditingSitesPerSample"].sum()

# %%
# # get all sites, across all genes, that have at least min_tot_cov coverage, per each individual replicate
# # (regardless of their editing frequency)
# sites_covered_per_gene_per_sample_per_individual_replicate_df = (
#     per_sample_agged_expanded_concat_edited_positions_df
#     .loc[
#         (per_sample_agged_expanded_concat_edited_positions_df["TotalCoverage"].ge(min_tot_cov))
#     ]
#     .groupby(["Chrom", condition_col, "Sample"])
#     .size()
#     .reset_index(name="NumOfSitesCoveredPerSample")
# )
# sites_covered_per_gene_per_sample_per_individual_replicate_df

# %%
# sites_covered_per_gene_per_sample_per_individual_replicate_df["Chrom"].nunique()

# %%
# sites_covered_per_gene_per_sample_in_all_replicates_df = (
#     per_sample_agged_expanded_concat_edited_positions_df
#     .loc[
#         (per_sample_agged_expanded_concat_edited_positions_df["TotalCoverage"].ge(min_tot_cov))
#     ]
#     .groupby(["Chrom", condition_col, "Position"])
#     .size()
#     .reset_index(name="NumOfSitesCoveredPerSample")
#     .loc[lambda x: x["NumOfSitesCoveredPerSample"].eq(7)]
#     .assign(Sample="All")
#     .groupby(["Chrom", condition_col, "Sample"])
#     .size()
#     .reset_index(name="NumOfSitesCoveredPerSample")
# )
# chroms_with_at_least_one_site_covered_in_all_replicates = set(sites_covered_per_gene_per_sample_in_all_replicates_df["Chrom"].unique())
# chroms = set(complete_data_df["Chrom"].tolist())
# chroms_with_not_even_one_site_covered_in_all_replicates = chroms - chroms_with_at_least_one_site_covered_in_all_replicates
# no_sites_covered_per_gene_per_sample_in_all_replicates_df = (
#     pd.DataFrame(
#         {
#             "Chrom": list(chroms_with_not_even_one_site_covered_in_all_replicates),
#         }
#     )
#     .merge(
#         orfs_df.loc[:, ["Chrom", "Name"]], 
#         on="Chrom", 
#         how="left"
#     )
#     .rename(
#         columns={"Name": condition_col}
#     )
#     .assign(
#         Sample="All",
#         NumOfSitesCoveredPerSample=0
#     )
# )
# sites_covered_per_gene_per_sample_in_all_replicates_df = (
#     pd.concat(
#         [
#             sites_covered_per_gene_per_sample_in_all_replicates_df,
#             no_sites_covered_per_gene_per_sample_in_all_replicates_df
#         ]
#     )
#     .sort_values("Chrom", ignore_index=True)
# )
# sites_covered_per_gene_per_sample_in_all_replicates_df

# %%
# sites_covered_per_gene_per_sample_in_all_replicates_df["Chrom"].nunique()

# %%
# # finally, we have a df with num of covered sites per gene per sample, 
# # including the counts for each individual replicate and the counts for all replicates together
# sites_covered_per_gene_per_sample_df = (
#     pd.concat(
#         [
#             sites_covered_per_gene_per_sample_per_individual_replicate_df,
#             sites_covered_per_gene_per_sample_in_all_replicates_df,
#         ]
#     )
#     .sort_values(
#         ["Chrom", condition_col, "Sample"],
#         ignore_index=True
#     )
# )
# sites_covered_per_gene_per_sample_df

# %%
# sites_covered_per_gene_per_sample_df["Chrom"].nunique()

# %%
# sites_covered_per_gene_per_sample_df.loc[
#     sites_covered_per_gene_per_sample_df["Sample"].eq("All"),
#     "NumOfSitesCoveredPerSample"
# ].describe().round(2)

# %%
# # merge the two dfs to get a df with 7 cols 
# # (MinTotalCoverage, MinEditingFrequency, condition_col, Replicate, NumOfEditingSitesPerReplicate, NumOfSitesCoveredPerReplicate, %EditedOfCoveredSites)
# num_of_sites_per_gene_per_sample_per_min_cov_and_editing_freq_df = num_of_sites_per_gene_per_sample_per_min_cov_and_editing_freq_df.merge(
#     sites_covered_per_gene_per_sample_df, 
#     on=["Chrom", condition_col, "Sample"]
# )
# # normalize the number of editing sites by the number of covered sites to get the percentage of edited sites out of the covered sites
# num_of_sites_per_gene_per_sample_per_min_cov_and_editing_freq_df["%EditedOfCoveredSites"] = (
#     100
#     * num_of_sites_per_gene_per_sample_per_min_cov_and_editing_freq_df["NumOfEditingSitesPerSample"]
#     / num_of_sites_per_gene_per_sample_per_min_cov_and_editing_freq_df["NumOfEditingSitesPerSample"]
# )
# num_of_sites_per_gene_per_sample_per_min_cov_and_editing_freq_df

# %%
# num_of_sites_per_gene_per_sample_per_min_cov_and_editing_freq_df.groupby("Sample")["%EditedOfCoveredSites"].describe().round(2)

# %%

# %%
per_sample_agged_expanded_concat_edited_positions_df

# %%
# min_tot_cov = min_tot_covs[0]
# min_editing_freq = min_editing_freqs[0]

min_tot_cov = min_tot_covs[1]
min_editing_freq = min_editing_freqs[1]

# %%
# create a df with 5 cols (MinTotalCoverage, MinEditingFrequency, condition_col, Position, Replicates)
# denoting for each gene-position the replicates that have at least min_tot_cov coverage and at least min_editing_freq editing frequency
replicates_per_gene_per_position_per_min_cov_and_editing_freq_df = (
    per_sample_agged_expanded_concat_edited_positions_df
    .loc[
        (per_sample_agged_expanded_concat_edited_positions_df["TotalCoverage"].ge(min_tot_cov))
        & (per_sample_agged_expanded_concat_edited_positions_df["EditingFrequency"].ge(min_editing_freq))
    ]
    .groupby(["Chrom", condition_col, "Position"])
    ["Sample"].apply(list)
    .reset_index(name="Samples")
)
replicates_per_gene_per_position_per_min_cov_and_editing_freq_df["Samples"] = (
        replicates_per_gene_per_position_per_min_cov_and_editing_freq_df["Samples"].apply(
        lambda x: x + ["All"] if len(x) == 7 else x
    )
)
replicates_per_gene_per_position_per_min_cov_and_editing_freq_df.insert(0, "MinTotalCoverage", min_tot_cov)
replicates_per_gene_per_position_per_min_cov_and_editing_freq_df.insert(1, "MinEditingFrequency", min_editing_freq)
replicates_per_gene_per_position_per_min_cov_and_editing_freq_df

# %%
replicates_per_gene_per_position_per_min_cov_and_editing_freq_df["Chrom"].nunique()

# %%
replicates_per_gene_per_position_per_min_cov_and_editing_freq_df["Samples"].apply(len).sum()

# %%
# create a df with 5 cols (MinTotalCoverage, MinEditingFrequency, Replicate, NumOfEditingSitesPerReplicate)
# denoting for each replicate the number of editing sites that have at least min_tot_cov coverage and at least min_editing_freq editing frequency
# across all 3 genes
num_of_sites_per_sample_per_min_cov_and_editing_freq_df = (
    replicates_per_gene_per_position_per_min_cov_and_editing_freq_df
    .explode("Samples")
    .rename(
        columns={"Samples": "Sample"}
    )
    .groupby(
        ["MinTotalCoverage", "MinEditingFrequency", "Sample",]
    )
    .size()
    .reset_index(name="NumOfEditingSitesPerSample")
)
num_of_sites_per_sample_per_min_cov_and_editing_freq_df

# %%
num_of_sites_per_sample_per_min_cov_and_editing_freq_df["NumOfEditingSitesPerSample"].sum()

# %%
# get all sites, across all genes, that have at least min_tot_cov coverage, per each individual replicate
# (regardless of their editing frequency)
sites_covered_per_individual_replicate_df = (
    per_sample_agged_expanded_concat_edited_positions_df
    .loc[
        (per_sample_agged_expanded_concat_edited_positions_df["TotalCoverage"].ge(min_tot_cov))
    ]
    .groupby("Sample")
    .size()
    .reset_index(name="NumOfSitesCoveredPerSample")
)
sites_covered_per_individual_replicate_df

# %%
# sites_covered_per_sample_in_all_replicates_df = (
#     per_sample_agged_expanded_concat_edited_positions_df
#     .loc[
#         (per_sample_agged_expanded_concat_edited_positions_df["TotalCoverage"].ge(min_tot_cov))
#     ]
#     # .groupby(["Chrom", condition_col, "Position"])
#     # .size()
#     # .reset_index(name="NumOfSitesCoveredPerSample")
#     # .loc[lambda x: x["NumOfSitesCoveredPerSample"].eq(7)]
#     # .assign(Sample="All")
#     # .groupby(["Chrom", condition_col, "Sample"])
#     # .size()
#     # .reset_index(name="NumOfSitesCoveredPerSample")
# )
# sites_covered_per_sample_in_all_replicates_df

# %%
sites_covered_per_gene_per_sample_in_all_replicates_df = (
    per_sample_agged_expanded_concat_edited_positions_df
    .loc[
        (per_sample_agged_expanded_concat_edited_positions_df["TotalCoverage"].ge(min_tot_cov))
    ]
    .groupby(["Chrom", condition_col, "Position"])
    .size()
    .reset_index(name="NumOfSitesCoveredPerSample")
    .loc[lambda x: x["NumOfSitesCoveredPerSample"].eq(7)]
    .assign(Sample="All")
    .groupby(["Chrom", condition_col, "Sample"])
    .size()
    .reset_index(name="NumOfSitesCoveredPerSample")
)
sites_covered_per_gene_per_sample_in_all_replicates_df

# %%
chroms_with_at_least_one_site_covered_in_all_replicates = set(sites_covered_per_gene_per_sample_in_all_replicates_df["Chrom"].unique())
chroms = set(complete_data_df["Chrom"].tolist())
chroms_with_not_even_one_site_covered_in_all_replicates = chroms - chroms_with_at_least_one_site_covered_in_all_replicates


# %%

no_sites_covered_per_gene_per_sample_in_all_replicates_df = (
    pd.DataFrame(
        {
            "Chrom": list(chroms_with_not_even_one_site_covered_in_all_replicates),
        }
    )
    .merge(
        orfs_df.loc[:, ["Chrom", "Name"]], 
        on="Chrom", 
        how="left"
    )
    .rename(
        columns={"Name": condition_col}
    )
    .assign(
        Sample="All",
        NumOfSitesCoveredPerSample=0
    )
)

# %%

# %%

# %%
site_covered_by_all_replicates = (
    per_sample_agged_expanded_concat_edited_positions_df
    .loc[
        (per_sample_agged_expanded_concat_edited_positions_df["TotalCoverage"].ge(min_tot_cov))
    ]
    .groupby(["Chrom", condition_col, "Position"])["Sample"]
    .nunique()
    .eq(7)
    .sum()
)
site_covered_by_all_replicates

# %%
sites_covered_per_sample_in_all_replicates_df = pd.DataFrame(
    {
        "Sample": ["All"],
        "NumOfSitesCoveredPerSample": [site_covered_by_all_replicates.sum()]
    }
)
# finally, we have a df with num of covered sites per sample, 
# including the counts for each individual replicate and the counts for all replicates together
sites_covered_per_sample_df = (
    pd.concat(
        [
            sites_covered_per_individual_replicate_df,
            sites_covered_per_sample_in_all_replicates_df,
        ]
    )
    .sort_values(
        "Sample",
        ignore_index=True
    )
)
sites_covered_per_sample_df

# %%
# merge the two dfs to get a df with 6 cols 
# (MinTotalCoverage, MinEditingFrequency, Replicate, NumOfEditingSitesPerReplicate, NumOfSitesCoveredPerReplicate, %EditedOfCoveredSites)
num_of_sites_per_sample_per_min_cov_and_editing_freq_df = num_of_sites_per_sample_per_min_cov_and_editing_freq_df.merge(
    sites_covered_per_sample_df, on="Sample"
)
# normalize the number of editing sites by the number of covered sites to get the percentage of edited sites out of the covered sites
num_of_sites_per_sample_per_min_cov_and_editing_freq_df["%EditedOfCoveredSites"] = (
    100
    * num_of_sites_per_sample_per_min_cov_and_editing_freq_df["NumOfEditingSitesPerSample"]
    / num_of_sites_per_sample_per_min_cov_and_editing_freq_df["NumOfSitesCoveredPerSample"]
)
num_of_sites_per_sample_per_min_cov_and_editing_freq_df

# %%

# %%

# %%

# %%
# num_of_sites_per_gene_per_sample_per_min_cov_and_editing_freq_dfs = []
num_of_sites_per_sample_per_min_cov_and_editing_freq_dfs = []

for min_tot_cov, min_editing_freq in product(min_tot_covs, min_editing_freqs):
    
    # print(f"Min total coverage: {min_tot_cov}, min editing frequency: {min_editing_freq}")
    
    # create a df with 5 cols (MinTotalCoverage, MinEditingFrequency, condition_col, Position, Replicates)
    # denoting for each gene-position the replicates that have at least min_tot_cov coverage and at least min_editing_freq editing frequency
    replicates_per_gene_per_position_per_min_cov_and_editing_freq_df = (
        per_sample_agged_expanded_concat_edited_positions_df
        .loc[
            (per_sample_agged_expanded_concat_edited_positions_df["TotalCoverage"].ge(min_tot_cov))
            & (per_sample_agged_expanded_concat_edited_positions_df["EditingFrequency"].ge(min_editing_freq))
        ]
        .groupby(["Chrom", condition_col, "Position"])
        ["Sample"].apply(list)
        .reset_index(name="Samples")
    )
    replicates_per_gene_per_position_per_min_cov_and_editing_freq_df["Samples"] = (
            replicates_per_gene_per_position_per_min_cov_and_editing_freq_df["Samples"].apply(
            lambda x: x + ["All"] if len(x) == 7 else x
        )
    )
    replicates_per_gene_per_position_per_min_cov_and_editing_freq_df.insert(0, "MinTotalCoverage", min_tot_cov)
    replicates_per_gene_per_position_per_min_cov_and_editing_freq_df.insert(1, "MinEditingFrequency", min_editing_freq)

    # # create a df with 6 cols (MinTotalCoverage, MinEditingFrequency, condition_col, Replicate, NumOfEditingSitesPerReplicate)
    # # denoting for each gene-replicate the number of editing sites that have at least min_tot_cov coverage and at least min_editing_freq editing frequency
    # num_of_sites_per_gene_per_sample_per_min_cov_and_editing_freq_df = (
    #     replicates_per_gene_per_position_per_min_cov_and_editing_freq_df
    #     .explode("Samples")
    #     .rename(
    #         columns={"Samples": "Sample"}
    #     )
    #     .groupby(
    #         ["MinTotalCoverage", "MinEditingFrequency", "Chrom", condition_col, "Sample",]
    #     )
    #     .size()
    #     .reset_index(name="NumOfEditingSitesPerSample")
    # )
    # # get all sites, across all genes, that have at least min_tot_cov coverage, per each individual replicate
    # # (regardless of their editing frequency)
    # sites_covered_per_gene_per_sample_per_individual_replicate_df = (
    #     per_sample_agged_expanded_concat_edited_positions_df
    #     .loc[
    #         (per_sample_agged_expanded_concat_edited_positions_df["TotalCoverage"].ge(min_tot_cov))
    #     ]
    #     .groupby(["Chrom", condition_col, "Sample"])
    #     .size()
    #     .reset_index(name="NumOfSitesCoveredPerSample")
    # )
    # sites_covered_per_gene_per_sample_in_all_replicates_df = (
    #     per_sample_agged_expanded_concat_edited_positions_df
    #     .loc[
    #         (per_sample_agged_expanded_concat_edited_positions_df["TotalCoverage"].ge(min_tot_cov))
    #     ]
    #     .groupby(["Chrom", condition_col, "Position"])
    #     .size()
    #     .reset_index(name="NumOfSitesCoveredPerSample")
    #     .loc[lambda x: x["NumOfSitesCoveredPerSample"].eq(7)]
    #     .assign(Sample="All")
    #     .groupby(["Chrom", condition_col, "Sample"])
    #     .size()
    #     .reset_index(name="NumOfSitesCoveredPerSample")
    # )
    # chroms_with_at_least_one_site_covered_in_all_replicates = set(sites_covered_per_gene_per_sample_in_all_replicates_df["Chrom"].unique())
    # chroms = set(complete_data_df["Chrom"].tolist())
    # chroms_with_not_even_one_site_covered_in_all_replicates = chroms - chroms_with_at_least_one_site_covered_in_all_replicates
    # no_sites_covered_per_gene_per_sample_in_all_replicates_df = (
    #     pd.DataFrame(
    #         {
    #             "Chrom": list(chroms_with_not_even_one_site_covered_in_all_replicates),
    #         }
    #     )
    #     .merge(
    #         orfs_df.loc[:, ["Chrom", "Name"]], 
    #         on="Chrom", 
    #         how="left"
    #     )
    #     .rename(
    #         columns={"Name": condition_col}
    #     )
    #     .assign(
    #         Sample="All",
    #         NumOfSitesCoveredPerSample=0
    #     )
    # )
    # sites_covered_per_gene_per_sample_in_all_replicates_df = (
    #     pd.concat(
    #         [
    #             sites_covered_per_gene_per_sample_in_all_replicates_df,
    #             no_sites_covered_per_gene_per_sample_in_all_replicates_df
    #         ]
    #     )
    #     .sort_values("Chrom", ignore_index=True)
    # )
    # # finally, we have a df with num of covered sites per gene per sample, 
    # # including the counts for each individual replicate and the counts for all replicates together
    # sites_covered_per_gene_per_sample_df = (
    #     pd.concat(
    #         [
    #             sites_covered_per_gene_per_sample_per_individual_replicate_df,
    #             sites_covered_per_gene_per_sample_in_all_replicates_df,
    #         ]
    #     )
    #     .sort_values(
    #         ["Chrom", condition_col, "Sample"],
    #         ignore_index=True
    #     )
    # )
    # # merge the two dfs to get a df with 7 cols 
    # # (MinTotalCoverage, MinEditingFrequency, condition_col, Replicate, NumOfEditingSitesPerReplicate, NumOfSitesCoveredPerReplicate, %EditedOfCoveredSites)
    # num_of_sites_per_gene_per_sample_per_min_cov_and_editing_freq_df = num_of_sites_per_gene_per_sample_per_min_cov_and_editing_freq_df.merge(
    #     sites_covered_per_gene_per_sample_df, 
    #     on=["Chrom", condition_col, "Sample"]
    # )
    # # normalize the number of editing sites by the number of covered sites to get the percentage of edited sites out of the covered sites
    # num_of_sites_per_gene_per_sample_per_min_cov_and_editing_freq_df["%EditedOfCoveredSites"] = (
    #     100
    #     * num_of_sites_per_gene_per_sample_per_min_cov_and_editing_freq_df["NumOfEditingSitesPerSample"]
    #     / num_of_sites_per_gene_per_sample_per_min_cov_and_editing_freq_df["NumOfEditingSitesPerSample"]
    # )

    # create a df with 5 cols (MinTotalCoverage, MinEditingFrequency, Replicate, NumOfEditingSitesPerReplicate)
    # denoting for each replicate the number of editing sites that have at least min_tot_cov coverage and at least min_editing_freq editing frequency
    # across all 3 genes
    num_of_sites_per_sample_per_min_cov_and_editing_freq_df = (
        replicates_per_gene_per_position_per_min_cov_and_editing_freq_df
        .explode("Samples")
        .rename(
            columns={"Samples": "Sample"}
        )
        .groupby(
            ["MinTotalCoverage", "MinEditingFrequency", "Sample",]
        )
        .size()
        .reset_index(name="NumOfEditingSitesPerSample")
    )
    # get all sites, across all genes, that have at least min_tot_cov coverage, per each individual replicate
    # (regardless of their editing frequency)
    sites_covered_per_individual_replicate_df = (
        per_sample_agged_expanded_concat_edited_positions_df
        .loc[
            (per_sample_agged_expanded_concat_edited_positions_df["TotalCoverage"].ge(min_tot_cov))
        ]
        .groupby("Sample")
        .size()
        .reset_index(name="NumOfSitesCoveredPerSample")
    )
    site_covered_by_all_replicates = (
        per_sample_agged_expanded_concat_edited_positions_df
        .loc[
            (per_sample_agged_expanded_concat_edited_positions_df["TotalCoverage"].ge(min_tot_cov))
        ]
        .groupby(["Chrom", condition_col, "Position"])["Sample"]
        .nunique()
        .eq(7)
        .sum()
    )
    sites_covered_per_sample_in_all_replicates_df = pd.DataFrame(
        {
            "Sample": ["All"],
            "NumOfSitesCoveredPerSample": [site_covered_by_all_replicates.sum()]
        }
    )
    # finally, we have a df with num of covered sites per sample, 
    # including the counts for each individual replicate and the counts for all replicates together
    sites_covered_per_sample_df = (
        pd.concat(
            [
                sites_covered_per_individual_replicate_df,
                sites_covered_per_sample_in_all_replicates_df,
            ]
        )
        .sort_values(
            "Sample",
            ignore_index=True
        )
    )
    # merge the two dfs to get a df with 6 cols 
    # (MinTotalCoverage, MinEditingFrequency, Replicate, NumOfEditingSitesPerReplicate, NumOfSitesCoveredPerReplicate, %EditedOfCoveredSites)
    num_of_sites_per_sample_per_min_cov_and_editing_freq_df = num_of_sites_per_sample_per_min_cov_and_editing_freq_df.merge(
        sites_covered_per_sample_df, on="Sample"
    )
    # normalize the number of editing sites by the number of covered sites to get the percentage of edited sites out of the covered sites
    num_of_sites_per_sample_per_min_cov_and_editing_freq_df["%EditedOfCoveredSites"] = (
        100
        * num_of_sites_per_sample_per_min_cov_and_editing_freq_df["NumOfEditingSitesPerSample"]
        / num_of_sites_per_sample_per_min_cov_and_editing_freq_df["NumOfSitesCoveredPerSample"]
    )
    
    # num_of_sites_per_gene_per_sample_per_min_cov_and_editing_freq_dfs.append(num_of_sites_per_gene_per_sample_per_min_cov_and_editing_freq_df)
    num_of_sites_per_sample_per_min_cov_and_editing_freq_dfs.append(num_of_sites_per_sample_per_min_cov_and_editing_freq_df)
    
    # break 

# concat_num_of_sites_per_gene_per_sample_per_min_cov_and_editing_freq_df = pd.concat(
#     num_of_sites_per_gene_per_sample_per_min_cov_and_editing_freq_dfs, 
#     ignore_index=True
# )
concat_num_of_sites_per_sample_per_min_cov_and_editing_freq_df = pd.concat(
    num_of_sites_per_sample_per_min_cov_and_editing_freq_dfs, 
    ignore_index=True
)

# %%
# concat_num_of_sites_per_gene_per_sample_per_min_cov_and_editing_freq_df

# %%
# concat_num_of_sites_per_gene_per_sample_per_min_cov_and_editing_freq_df.groupby("Sample")["%EditedOfCoveredSites"].describe().round(2)

# %%
concat_num_of_sites_per_sample_per_min_cov_and_editing_freq_df

# %%
# concat_num_of_sites_per_sample_per_min_cov_and_editing_freq_df.groupby("Sample")["%EditedOfCoveredSites"].describe().round(2)

# %%
# per_sample_agged_expanded_concat_edited_positions_df

# %%
# # aggreated measure with num of sites supported by numer of samples

# num_of_samples_per_sites_per_min_cov_and_editing_freq_dfs = []
# for min_tot_cov, min_editing_freq in product(min_tot_covs, min_editing_freqs):
#     print(f"Min total coverage: {min_tot_cov}, min editing frequency: {min_editing_freq}")
#     num_of_samples_per_sites_per_min_cov_and_editing_freq_df = per_sample_agged_expanded_concat_edited_positions_df.loc[
#         (per_sample_agged_expanded_concat_edited_positions_df["TotalCoverage"].ge(min_tot_cov))
#         & (per_sample_agged_expanded_concat_edited_positions_df["EditingFrequency"].ge(min_editing_freq))
#     ].groupby(
#         ["Chrom", "Transcript", "Position"]
#     )["Sample"].nunique().value_counts().reset_index()
#     num_of_samples_per_sites_per_min_cov_and_editing_freq_df.insert(0, "MinTotalCoverage", min_tot_cov)
#     num_of_samples_per_sites_per_min_cov_and_editing_freq_df.insert(1, "MinEditingFrequency", min_editing_freq)
    
#     # num_of_all_sites_covered_by_at_least_one_sample = per_sample_agged_expanded_concat_edited_positions_df.loc[
#     #     (per_sample_agged_expanded_concat_edited_positions_df["TotalCoverage"].ge(min_tot_cov))
#     # ].drop_duplicates(
#     #     ["Chrom", "Transcript", "Position"]
#     # ).shape[0]
#     # num_of_samples_per_sites_per_min_cov_and_editing_freq_df["NumOfSitesCoveredByAtLeastOneSample"] = num_of_all_sites_covered_by_at_least_one_sample
    
#     num_of_samples_per_sites_per_min_cov_and_editing_freq_df = num_of_samples_per_sites_per_min_cov_and_editing_freq_df.rename(
#         columns={
#             "count": "NumOfEditingSites",
#             "Sample": "NumOfSamples"
#         }
#     )
#     num_of_samples_per_sites_per_min_cov_and_editing_freq_df = num_of_samples_per_sites_per_min_cov_and_editing_freq_df.sort_values(
#         "NumOfSamples", ascending=False, ignore_index=True
#     )
#     num_of_samples_per_sites_per_min_cov_and_editing_freq_df["ReverseCumulativeNumOfEditingSites"] = (
#         num_of_samples_per_sites_per_min_cov_and_editing_freq_df["NumOfEditingSites"].cumsum()
#     )
#     num_of_samples_per_sites_per_min_cov_and_editing_freq_df["%ReverseCumulativeNumOfEditingSites"] = (
#         100
#         * num_of_samples_per_sites_per_min_cov_and_editing_freq_df["ReverseCumulativeNumOfEditingSites"]
#         / num_of_samples_per_sites_per_min_cov_and_editing_freq_df["NumOfEditingSites"].sum()
#     )
    
#     num_of_samples_per_sites_per_min_cov_and_editing_freq_dfs.append(num_of_samples_per_sites_per_min_cov_and_editing_freq_df)
    
# concat_num_of_samples_per_sites_per_min_cov_and_editing_freq_df = pd.concat(
#     num_of_samples_per_sites_per_min_cov_and_editing_freq_dfs, ignore_index=True
# )

# # concat_num_of_samples_per_sites_per_min_cov_and_editing_freq_df["%EditedSitesOfAllSitesCoveredByAtLeastOneSample"] = (
# #         100 
# #         * concat_num_of_samples_per_sites_per_min_cov_and_editing_freq_df["NumOfEditingSites"] 
# #         / concat_num_of_samples_per_sites_per_min_cov_and_editing_freq_df["NumOfSitesCoveredByAtLeastOneSample"]
# #     )

# concat_num_of_samples_per_sites_per_min_cov_and_editing_freq_df

# %%
fig = px.bar(
    concat_num_of_sites_per_sample_per_min_cov_and_editing_freq_df,
    x="Sample",
    y="%EditedOfCoveredSites",
    facet_row="MinTotalCoverage",
    facet_col="MinEditingFrequency",
    labels={
        # "NumOfEditingSitesPerSample": "Editing sites",
        "%EditedOfCoveredSites": "Edited /<br>covered sites [%]",
    },
    color="Sample",
    color_discrete_map=samples_color_discrete_map,
    category_orders={
        "Sample": samples.tolist() + ["All"],
    }
    # log_y=True
)

# Loop through all facet annotations and modify their text
for annotation in fig.layout.annotations:
    if "MinTotalCoverage=" in annotation.text:
        annotation.text = annotation.text.replace("MinTotalCoverage=", "Coverage / sample ≥ ")
    if "MinEditingFrequency=" in annotation.text:
        annotation.text = annotation.text.replace("MinEditingFrequency=", "Editing frequency / sample ≥ ")
        
# fig.update_yaxes(dtick=5000)

width = 1200
height = 900

fig.update_layout(
    template=template,
    width=width,
        height=height,
    title="Number of sites edited per sample",
    showlegend=False
)

fig.write_image(
    Path(
        out_dir,
        "Number of sites edited per sample - Octopus - pooled.svg",
    ),
    width=width,
    height=height,
)

fig.show()

# %%
# fig = px.line(
#     concat_num_of_samples_per_sites_per_min_cov_and_editing_freq_df,
#     x="NumOfSamples",
#     y="%ReverseCumulativeNumOfEditingSites",
#     facet_row="MinTotalCoverage",
#     facet_col="MinEditingFrequency",
#     facet_row_spacing=0.05,
#     labels={
#         "NumOfSamples": "Samples",
#         "%ReverseCumulativeNumOfEditingSites": "% of edited sites",
#     },
#     markers=True
# )

# # Loop through all facet annotations and modify their text
# for annotation in fig.layout.annotations:
#     if "MinTotalCoverage=" in annotation.text:
#         annotation.text = annotation.text.replace("MinTotalCoverage=", "Coverage / sample ≥ ")
#     if "MinEditingFrequency=" in annotation.text:
#         annotation.text = annotation.text.replace("MinEditingFrequency=", "Editing freq / sample ≥ ")

# # fig.update_xaxes(dtick=1)
# fig.update_xaxes(dtick=1, autorange="reversed")
# fig.update_yaxes(range=[0, 105], dtick=25)

# fig.update_layout(
#     template=template,
#     width=800,
#     height=800,
#     title="Cumulative % of edited sites supported by X samples or less"
# )
# fig.show()

# %%
# fig = px.line(
#     concat_num_of_samples_per_sites_per_min_cov_and_editing_freq_df.loc[
#         (concat_num_of_samples_per_sites_per_min_cov_and_editing_freq_df["MinTotalCoverage"] == 0)
#         & (concat_num_of_samples_per_sites_per_min_cov_and_editing_freq_df["MinEditingFrequency"] == 0.01)
#     ],
#     x="NumOfSamples",
#     y="%ReverseCumulativeNumOfEditingSites",
#     # facet_row="MinTotalCoverage",
#     # facet_col="MinEditingFrequency",
#     # facet_row_spacing=0.05,
#     labels={
#         "NumOfSamples": "Samples",
#         "%ReverseCumulativeNumOfEditingSites": "% of edited sites",
#     },
#     markers=True
# )

# # # Loop through all facet annotations and modify their text
# # for annotation in fig.layout.annotations:
# #     if "MinTotalCoverage=" in annotation.text:
# #         annotation.text = annotation.text.replace("MinTotalCoverage=", "Coverage / sample ≥ ")
# #     if "MinEditingFrequency=" in annotation.text:
# #         annotation.text = annotation.text.replace("MinEditingFrequency=", "Editing freq / sample ≥ ")

# # fig.update_xaxes(dtick=1)
# fig.update_xaxes(dtick=1, autorange="reversed")
# fig.update_yaxes(range=[0, 105], dtick=25)

# width = 400
# height = 400

# fig.update_layout(
#     template=template,
#     width=width,
#     height=height,
#     title="Cumulative % of edited sites supported<br>by X samples or less"
# )

# fig.write_image(
#     Path(
#         out_dir,
#         "Cumulative % of edited sites supported by X samples or less - Octopus - pooled.svg",
#     ),
#     width=width,
#     height=height,
# )

# fig.show()

# %%
# raise Exception("Stop notebook execution here")

# %% [markdown]
# ### Noise in positions

# %% [markdown]
# Only show estimated noise levels for genes with max 3 SNPs.

# %%
(
    concat_all_positions_df.loc[
        (concat_all_positions_df["Noise"] <= snp_noise_level)
        & (concat_all_positions_df["NoisyFinal"])
        & (concat_all_positions_df["Chrom"].isin(tmr50_alignment_stats_and_snps_df["Chrom"].values)) # only consider genes with <= 3 SNPs
    ]
    .groupby("Chrom")
    .size()
    .reset_index()
    .rename(columns={0: "NoisePositions"})["NoisePositions"]
    .describe()
)


# %%
def mean_noise_levels(positions_df, top_x_noisy_positions=3, snp_noise_level=0.05):
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
    .apply(mean_noise_levels, 3, snp_noise_level, include_groups=False)
    .reset_index()
    .rename(columns={0: "Noise"})
    # .merge(tmr50_alignment_stats_df.loc[:, ["Chrom"]], on="Chrom", how="right")
    .merge(tmr50_alignment_stats_and_snps_df.loc[:, ["Chrom"]], on="Chrom", how="right")
    .fillna(0.0)
    .sort_values(["Chrom", "Noise"])
)
all_per_chrom_mean_noise_levels["%Noise"] = (
    100 * all_per_chrom_mean_noise_levels["Noise"]
)
all_per_chrom_mean_noise_levels

# %%
all_per_chrom_mean_noise_levels["Noise"].describe()

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
saved_per_chrom_mean_noise_levels_df.to_csv(
    Path(out_dir, "Noise.Octopus.WholeTranscriptome.Pooled.tsv"), 
    sep="\t", 
    index=False
)
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
# ### Mismatches by type (old, redistribute to other sections)

# %%
significant_mismatches_df = mismatches_df.loc[
    ((mismatches_df["EditedFinal"]) | (mismatches_df["NoisyFinal"]))
    & (mismatches_df["MismatchFrequency"].ge(mismatches_df["NoiseThreshold"]))
    # & (mismatches_df["Chrom"].isin(chroms))
]
significant_mismatches_df

# %%
# significant_mismatches_above_noise_threshold_df = significant_mismatches_df.loc[
#     significant_mismatches_df["MismatchFrequency"].ge(significant_mismatches_df["NoiseThreshold"])

# ]
# # significant_mismatches_df["MismatchFrequency>10%"] = significant_mismatches_df["MismatchFrequency"] > 0.1
# significant_mismatches_above_noise_threshold_df

# %%
assert significant_mismatches_df.loc[
    significant_mismatches_df["EditedFinal"],
    "Chrom"
].nunique() == complete_data_df.shape[0]

# %%
significant_mismatches_df.loc[
    significant_mismatches_df["EditedFinal"],
    "Chrom"
].nunique()

# %%
true_a2g_mismatches = significant_mismatches_df.loc[
    (significant_mismatches_df["EditedFinal"])
    & (significant_mismatches_df["Mismatch"] == "A>G")
    # & (significant_mismatches_df["Chrom"].isin(chroms))
]
true_a2g_mismatches

# %%
wrong_a2g_mismatches = significant_mismatches_df.loc[
    (significant_mismatches_df["EditedFinal"])
    & (significant_mismatches_df["Mismatch"] != "A>G")
    # & (significant_mismatches_df["Chrom"].isin(chroms))
]
wrong_a2g_mismatches

# %%
100 * len(wrong_a2g_mismatches) / len(true_a2g_mismatches)

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
significant_mismatches_df.loc[
    (significant_mismatches_df["EditedFinal"])
].groupby("Chrom").size().describe().round(2)

# %%
significant_mismatches_df.loc[
    (significant_mismatches_df["NoisyFinal"])
    # & (significant_mismatches_df["Chrom"].isin(chroms))
].groupby("Chrom").size().describe().round(2)

# %%
# significant_mismatches_above_noise_threshold_df.loc[
#     (significant_mismatches_above_noise_threshold_df["NoisyFinal"])
#     # & (significant_mismatches_df["Chrom"].isin(chroms))
# ].groupby("Chrom").size().describe().round(2)

# %%

# %%
# todo what defines a SNP? above snp_noise_level or above noise threshold? update after deciding
snps_per_gene_df = significant_mismatches_df.loc[
    (significant_mismatches_df["NoisyFinal"])
    & (significant_mismatches_df["MismatchFrequency"] >= snp_noise_level)
].groupby("Chrom").size().reset_index(name="SNPs")
snps_per_gene_df

# %%
editing_sites_per_gene_df = significant_mismatches_df.loc[
    significant_mismatches_df["EditedFinal"]
].groupby("Chrom").size().reset_index(name="EditingSites")
editing_sites_per_gene_df

# %%
tmr50_alignment_stats_df

# %%
editing_sites_and_snps_and_coverage_per_gene_df = (
    tmr50_alignment_stats_df.loc[:, ["Chrom", "MappedReads", "KnownSites"]]
    .merge(
        snps_per_gene_df,
        # how="left"
        how="outer"
    )
    .merge(
        editing_sites_per_gene_df,
        # how="left",
        how="outer"
    )
    .fillna(0)
    .sort_values("Chrom", ignore_index=True)
)
editing_sites_and_snps_and_coverage_per_gene_df["SNPs"] = (
    editing_sites_and_snps_and_coverage_per_gene_df["SNPs"].astype(int)
)
editing_sites_and_snps_and_coverage_per_gene_df

# %%
editing_sites_and_snps_and_coverage_per_edited_gene_df = editing_sites_and_snps_and_coverage_per_gene_df.loc[
    editing_sites_and_snps_and_coverage_per_gene_df["EditingSites"].gt(0)
]

editing_sites_and_snps_and_coverage_per_edited_gene_df

# %%
max_known_sites = editing_sites_and_snps_and_coverage_per_gene_df["KnownSites"].max()
max_current_sites = editing_sites_and_snps_and_coverage_per_gene_df["EditingSites"].max()
fig = px.scatter(
    editing_sites_and_snps_and_coverage_per_gene_df,
    x="KnownSites",
    y="EditingSites",
    marginal_x="histogram",
    marginal_y="histogram",
    trendline="ols",
    trendline_color_override="black"
)
fig.update_xaxes(dtick=50)
fig.update_yaxes(dtick=50)
height = 400
width = height * max_known_sites / max_current_sites
fig.update_layout(
    height=height,
    width=width, 
    title="Known vs. current editing sites per gene"
)
fig.show()

# %%
max_known_sites = editing_sites_and_snps_and_coverage_per_edited_gene_df["KnownSites"].max()
max_current_sites = editing_sites_and_snps_and_coverage_per_edited_gene_df["EditingSites"].max()
fig = px.scatter(
    editing_sites_and_snps_and_coverage_per_edited_gene_df,
    x="KnownSites",
    y="EditingSites",
    marginal_x="histogram",
    marginal_y="histogram",
    trendline="ols",
    trendline_color_override="black"
)
fig.update_xaxes(dtick=50)
fig.update_yaxes(dtick=50)
height = 400
width = height * max_known_sites / max_current_sites
fig.update_layout(
    height=height,
    width=width, 
    title="Known vs. current editing sites per gene (currently-edited genes only)"
)
fig.show()

# %%
editing_sites_and_snps_and_coverage_per_edited_gene_df["EditingSites"].sum()

# %%
editing_sites_and_snps_and_coverage_per_edited_gene_df["SNPs"].sum()

# %%
editing_sites_and_snps_and_coverage_per_edited_gene_df["EditingSites"].describe().round(2)

# %%
editing_sites_and_snps_and_coverage_per_edited_gene_df["SNPs"].describe().round(2)

# %%
assert editing_sites_and_snps_and_coverage_per_edited_gene_df.loc[
    editing_sites_and_snps_and_coverage_per_edited_gene_df["SNPs"] > max_snps_per_gene_to_allow_editing_detection,
].empty

# %%
fig = px.histogram(
    editing_sites_and_snps_and_coverage_per_edited_gene_df,
    x="SNPs",
    # log_y=True,
    # histnorm="percent",
    cumulative=True,
)
fig.update_xaxes(dtick=25)
fig.update_layout(width=600, height=450)
fig.show()

# %%
# fig = px.histogram(
#     # noise_sites_per_gene,
#     editing_sites_and_snps_and_coverage_per_edited_gene_df.loc[
#         editing_sites_and_snps_and_coverage_per_edited_gene_df["SNPs"].le(10)
#     ],
#     x="SNPs",
#     # log_y=True
#     cumulative=True,
# )
# fig.update_xaxes(dtick=1)
# # fig.update_yaxes(dtick=500, range=[0, 2500])
# fig.update_layout(width=600, height=450)
# fig.show()

# %%
# min_mapped = 1000
# editing_sites_and_snps_and_coverage_per_edited_gene_df.loc[
#     (editing_sites_and_snps_and_coverage_per_edited_gene_df["SNPs"].le(3))
#     & (editing_sites_and_snps_and_coverage_per_edited_gene_df["MappedReads"].ge(min_mapped)),
#     ["SNPs"]
# ].value_counts().reset_index(name=f"GenesWith{min_mapped}+MappedReads").sort_values("SNPs", ignore_index=True)

# %%
# min_mapped = 800
# editing_sites_and_snps_and_coverage_per_edited_gene_df.loc[
#     (editing_sites_and_snps_and_coverage_per_edited_gene_df["SNPs"].le(3))
#     & (editing_sites_and_snps_and_coverage_per_edited_gene_df["MappedReads"].ge(min_mapped)),
#     ["SNPs"]
# ].value_counts().reset_index(name=f"GenesWith{min_mapped}+MappedReads").sort_values("SNPs", ignore_index=True)

# %%
# fig = px.histogram(
#     # noise_sites_per_gene,
#     editing_sites_and_snps_and_coverage_per_edited_gene_df.loc[
#         editing_sites_and_snps_and_coverage_per_edited_gene_df["SNPs"].le(3)
#     ],
#     x="MappedReads",
#     facet_col="SNPs",
#     facet_col_spacing=0.06,
#     nbins=20,
#     # facet
#     # log_y=True
#     # cumulative=True,
# )
# fig.update_xaxes(tickangle=45, dtick=1000)
# fig.update_layout(width=1200, height=450)
# fig.show()

# %%
# fig = px.histogram(
#     # noise_sites_per_gene,
#     editing_sites_and_snps_and_coverage_per_edited_gene_df.loc[
#         # editing_sites_and_snps_and_coverage_per_edited_gene_df["SNPs"].le(10)
#         editing_sites_and_snps_and_coverage_per_edited_gene_df["SNPs"].le(19)
#     ],
#     x="EditingSites",
#     facet_col="SNPs",
#     facet_col_spacing=0.06,
#     facet_col_wrap=10,
#     facet_row_spacing=0.12,
#     nbins=20,
#     # facet
#     # log_y=True
#     cumulative=True,
# )
# # fig.update_xaxes(tickangle=45, dtick=25)
# fig.update_layout(width=1400, height=600)
# fig.show()

# %%
# editing_sites_and_snps_and_coverage_per_edited_gene_df.loc[
#     editing_sites_and_snps_and_coverage_per_edited_gene_df["SNPs"].le(10)
# ].groupby(
#     "SNPs"
# ).agg(
#     # MeanMappedReads=("MappedReads", "mean"),
#     # STDMappedReads=("MappedReads", "std"),
#     MeanEditingSites=("EditingSites", "mean"),
#     STDEditingSites=("EditingSites", "std"),
#     SumEditingSites=("EditingSites", "sum"),
# ).round(2).reset_index()

# %%
fig = px.box(
    editing_sites_and_snps_and_coverage_per_edited_gene_df.loc[
        editing_sites_and_snps_and_coverage_per_edited_gene_df["SNPs"].le(10)
    ],
    x="SNPs",
    y="EditingSites",
    # facet_col="SNPs",
    # facet_col_spacing=0.06,
    # facet_col_wrap=10,
    # facet_row_spacing=0.12,
    # nbins=20,
    # # facet
    log_y=True,
    # cumulative=True,
)
# fig.update_xaxes(tickangle=45, dtick=25)
fig.update_layout(width=800, height=450)
fig.show()

# %%
# fig = px.scatter(
#     editing_sites_and_snps_and_coverage_per_edited_gene_df.loc[
#         # editing_sites_and_snps_and_coverage_per_edited_gene_df["SNPs"].le(10)
#         :
#     ].groupby(
#         "SNPs"
#     ).agg(
#         # MeanMappedReads=("MappedReads", "mean"),
#         # STDMappedReads=("MappedReads", "std"),
#         MeanEditingSites=("EditingSites", "mean"),
#         STDEditingSites=("EditingSites", "std"),
#         SumEditingSites=("EditingSites", "sum"),
#     ).round(2).reset_index(),
#     x="SNPs",
#     y="MeanEditingSites",
#     error_y="STDEditingSites",
#     # markers=True
# )
# fig.update_xaxes(dtick=50)
# fig.update_yaxes(dtick=50)
# fig.update_layout(width=650, height=400)
# # fig.update_layout(width=1000, height=600)
# fig.show()

# %%
fig = px.scatter(
    editing_sites_and_snps_and_coverage_per_edited_gene_df,
    x="SNPs",
    y="EditingSites",
    # error_y="STDEditingSites",
    trendline="ols",
    trendline_color_override="black"
)
fig.update_xaxes(dtick=50)
fig.update_yaxes(dtick=50)
fig.update_layout(width=650, height=400)
fig.show()

# %%
fig = px.line(
    editing_sites_and_snps_and_coverage_per_edited_gene_df.loc[
        editing_sites_and_snps_and_coverage_per_edited_gene_df["SNPs"].le(10)
        # :
    ].groupby(
        "SNPs"
    ).agg(
        # MeanMappedReads=("MappedReads", "mean"),
        # STDMappedReads=("MappedReads", "std"),
        MeanEditingSites=("EditingSites", "mean"),
        STDEditingSites=("EditingSites", "std"),
        SumEditingSites=("EditingSites", "sum"),
    ).round(2).reset_index(),
    x="SNPs",
    y="MeanEditingSites",
    error_y="STDEditingSites",
    markers=True
    # ols
)
fig.update_xaxes(dtick=1)
fig.update_yaxes(dtick=5)
fig.update_layout(width=800, height=500)
fig.show()

# %%

# %%

# %%

# %%

# %% [markdown]
# ### 12 non-SNP mismatches distribution

# %% [markdown]
# For this analysis, we only consider genes in which we detected editing.  
# This excludes genes with too many SNPs where we didn't even look for editing, or genes we looked at but didn't find any editing.  
# Also, one needs to remember that by our definition of editing, a site is considered an A->G if its mismatch frequency
# is above the noise threshold. Therefore, not many sites of other mismatches can possibly be above the noise threshold.

# %%
mismatches_df

# %%
significant_non_snp_mismatches_df = mismatches_df.loc[
    (
        (mismatches_df["EditedFinal"])
        | (
            (mismatches_df["NoisyFinal"])
            & (mismatches_df["MismatchFrequency"].ge(mismatches_df["NoiseThreshold"]))
            & (mismatches_df["MismatchFrequency"].lt(snp_noise_level))
        )
    )
    & (mismatches_df["Chrom"].isin(chroms))
]
significant_non_snp_mismatches_df = significant_non_snp_mismatches_df.sort_values("Mismatch", ignore_index=True)
significant_non_snp_mismatches_df

# %%
significant_non_snp_mismatches_df.to_csv(
    Path(out_dir, "12NonSNPMismatchesAboveNoiseThreshold.Octopus.WholeTranscriptome.Pooled.csv"),
    sep="\t",
    index=False
)

# %% [markdown]
# ### 12 mismatches distribution

# %% [markdown]
# For this analysis, we only consider genes in which we detected editing.  
# This excludes genes with too many SNPs where we didn't even look for editing, or genes we looked at but didn't find any editing.  
# Also, one needs to remember that by our definition of editing, a site is considered an A->G if its mismatch frequency
# is above the noise threshold. Therefore, not many sites of other mismatches can possibly be above the noise threshold.

# %%
significant_mismatches_df = mismatches_df.loc[
    ((mismatches_df["EditedFinal"]) | (mismatches_df["NoisyFinal"]))
    & (mismatches_df["MismatchFrequency"].ge(mismatches_df["NoiseThreshold"]))
    & (mismatches_df["Chrom"].isin(chroms))
]
significant_mismatches_df = significant_mismatches_df.sort_values("Mismatch", ignore_index=True)
significant_mismatches_df

# %%
significant_mismatches_df.to_csv(
    Path(out_dir, "12MismatchsAboveNoiseThreshold.Octopus.WholeTranscriptome.Pooled.csv"),
    sep="\t",
    index=False
)

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
    title="Total number of significant sites per mismatch across<br>all octopus genes",
)

width = 600
height = 500

# Reduce opacity to see both histograms
# fig.update_traces(opacity=0.75)
fig.update_layout(
    width=width,
    height=height,
    showlegend=False
    # barmode='overlay' # Overlay both histograms
)

fig.show()

# %%
# Per gene: (A>G sites / all sites in gene) divided by (top other mismatch sites / all sites in gene)
# == A>G_count / top_other_mismatch_count (also returns the underlying fractions)

per_gene_mismatch_ratio_df = (
    significant_mismatches_df.loc[
        significant_mismatches_df["Chrom"].isin(chroms), 
        ["Chrom", "Mismatch"]
    ]
    .value_counts()
    .reset_index(name="NumSites")
    .merge(
        (
            significant_mismatches_df.loc[
                significant_mismatches_df["Chrom"].isin(chroms)
            ]
           .groupby("Chrom")
           .size()
           .reset_index(name="NumSitesInGene")
        ),
        on="Chrom",
        how="left",
    )
    .assign(FractionInGene=lambda x: x["NumSites"] / x["NumSitesInGene"])
)

a2g_per_gene_df = (
    per_gene_mismatch_ratio_df.loc[
        per_gene_mismatch_ratio_df["Mismatch"] == "A>G", 
        ["Chrom", "NumSites", "FractionInGene"]
    ]
    .rename(
        columns={
            "NumSites": "NumA2GSites",
            "FractionInGene": "FracA2G",
        }
    )
)

top_other_per_gene_df = (
    per_gene_mismatch_ratio_df.loc[
        per_gene_mismatch_ratio_df["Mismatch"] != "A>G", 
        ["Chrom", "Mismatch", "NumSites", "FractionInGene"]
    ]
    .sort_values(["Chrom", "NumSites"], ascending=[True, False])
    .drop_duplicates("Chrom", keep="first", ignore_index=True)
    .rename(
        columns={
            "Mismatch": "TopOtherMismatch",
            "NumSites": "NumTopOtherSites",
            "FractionInGene": "FracTopOther",
        }
    )
)

# per_gene_a2g_over_top_other_ratio_df = (
#     significant_mismatches_df.loc[
#         significant_mismatches_df["Chrom"].isin(chroms)
#     ]
#     .groupby("Chrom").size().reset_index(name="NumSitesInGene")
#     .merge(a2g_per_gene_df, on="Chrom", how="left")
#     .merge(top_other_per_gene_df, on="Chrom", how="left")
#     .fillna({"NumA2GSites": 0, "FracA2G": 0.0})
#     # .fillna({"NumA2GSites": 0, "FracA2G": 0.0, "NumTopOtherSites": 0, "FracTopOther": 0.0})
#     .assign(
#         A2G_over_TopOther_FracRatio=lambda x: x["FracA2G"] / x["FracTopOther"],
#         A2G_over_TopOther_CountRatio=lambda x: x["NumA2GSites"] / x["NumTopOtherSites"],
#         SNR=lambda x: x["NumA2GSites"] / (x["NumA2GSites"] + x["NumTopOtherSites"])
#     )
#     # # .sort_values("A2G_over_TopOther_FracRatio", ascending=False, ignore_index=True)
#     .rename(columns={"SNR": "AG / (AG + top other)"})
#     .sort_values("AG / (AG + top other)", ascending=False, ignore_index=True)
# )

per_gene_a2g_over_top_other_ratio_df = (
    significant_mismatches_df.loc[
        significant_mismatches_df["Chrom"].isin(chroms)
    ]
    .groupby("Chrom").size().reset_index(name="NumSitesInGene")
    .merge(a2g_per_gene_df, on="Chrom", how="left")
    .merge(top_other_per_gene_df, on="Chrom", how="left")
    # .fillna({"NumA2GSites": 0, "FracA2G": 0.0})
    .fillna({"NumA2GSites": 0, "FracA2G": 0.0, "NumTopOtherSites": 0, "FracTopOther": 0.0})
    # .assign(
    #     A2G_over_TopOther_FracRatio=lambda x: x["FracA2G"] / x["FracTopOther"],
    #     A2G_over_TopOther_CountRatio=lambda x: x["NumA2GSites"] / x["NumTopOtherSites"],
    #     SNR=lambda x: x["NumA2GSites"] / (x["NumA2GSites"] + x["NumTopOtherSites"])
    # )
    # # # .sort_values("A2G_over_TopOther_FracRatio", ascending=False, ignore_index=True)
    # .rename(columns={"SNR": "AG / (AG + top other)"})
    # .sort_values("AG / (AG + top other)", ascending=False, ignore_index=True)
)
assert per_gene_a2g_over_top_other_ratio_df["NumA2GSites"].gt(0).all()
per_gene_a2g_over_top_other_ratio_df["A2G_over_TopOther_FracRatio"] = per_gene_a2g_over_top_other_ratio_df.apply(
    lambda x: 
        1 
        if pd.isna(x["TopOtherMismatch"]) 
        else x["FracA2G"] / x["FracTopOther"],
    axis=1
)
per_gene_a2g_over_top_other_ratio_df["AG / (AG + top other)"] = per_gene_a2g_over_top_other_ratio_df.apply(
    lambda x: 
        x["NumA2GSites"] / (x["NumA2GSites"] + x["NumTopOtherSites"]),
    axis=1
)


per_gene_a2g_over_top_other_ratio_df

# %%
# num of edited sites w/o any snp
per_gene_a2g_over_top_other_ratio_df.loc[
    per_gene_a2g_over_top_other_ratio_df["TopOtherMismatch"].isna()
].shape[0]

# %%
fig = px.histogram(
    per_gene_a2g_over_top_other_ratio_df,
    x="AG / (AG + top other)",
    title="Cumulative % of genes with SNR below given value",
    histnorm="percent",
    cumulative=True,
    nbins=20
)
fig.update_xaxes(dtick=0.1)
fig.update_yaxes(dtick=10)
fig.update_layout(
    width=600,
    height=500,
    showlegend=False
)

fig.show()

# %%
fig = px.scatter(
    per_gene_a2g_over_top_other_ratio_df,
    x="AG / (AG + top other)",
    y="NumA2GSites",
    # title="Cumulative % of genes with SNR below given value",
    # histnorm="percent",
    # cumulative=True,
    # nbins=20
)
fig.update_xaxes(dtick=0.1)
# fig.update_yaxes(dtick=10)
fig.update_layout(
    width=600,
    height=500,
    showlegend=False
)

fig.show()

# %%

# %% [markdown]
# ### Haplotype-based editing

# %% [markdown]
# #### Genes enriched with SNPs

# %%
# at_least_x_snps_per_gene_options = [6, 8, 10]
at_least_x_mapped_reads_per_gene_options = [200, 600, 1000]
# min_profound_snp_freq = 0.2
# max_profound_snp_freq = 0.8
min_profound_snp_freq = 0.3
max_profound_snp_freq = 0.7

# take SNPs with frequency between min_profound_snp_freq and max_profound_snp_freq in genes with significant 
# editing and noise, and annotate number of SNPs per gene and at least x mapped reads per gene
enriched_snps_df = significant_mismatches_df.loc[
    (significant_mismatches_df["Chrom"].isin(chroms))
    & (significant_mismatches_df["NoisyFinal"])
    & (significant_mismatches_df["MismatchFrequency"].ge(min_profound_snp_freq))
    & (significant_mismatches_df["MismatchFrequency"].le(max_profound_snp_freq)),
]

enriched_snps_df = enriched_snps_df.drop(
    columns=[
        'EditingFrequency', 'Edited', 'EditedCorrected',
        'EditedFinal', 'Noise', 'NoisyCorrected', 'NoisyFinal',
    ]
)

enriched_snps_df["SNPsPerGene"] = enriched_snps_df.groupby("Chrom").transform("size")

# enriched_snps_df = enriched_snps_df.loc[
#     enriched_snps_df["SNPsPerGene"].eq(enriched_snps_df["SNPsPerGene"].max())
# ]

# # annotate at least x SNPs in gene
# for at_least_x_snps in at_least_x_snps_per_gene_options:
#     snps_per_gene_flag_df = (
#         enriched_snps_df
#         .groupby("Chrom")
#         .size()
#         .reset_index(name="NumSNPsInGene")
#         .assign(**{f"AtLeast{at_least_x_snps}SNPsInGene": lambda x: x["NumSNPsInGene"] >= at_least_x_snps})
#         .loc[:, ["Chrom", f"AtLeast{at_least_x_snps}SNPsInGene"]]
#     )
#     enriched_snps_df = enriched_snps_df.merge(
#         snps_per_gene_flag_df,
#         on="Chrom",
#         how="left",
    # )
    
# # annotate num of SNPs per gene
# enriched_snps_df = enriched_snps_df.merge(
#     enriched_snps_df.groupby("Chrom").size().reset_index(name="SNPsPerGene"),
#     on="Chrom",
#     how="left",
# )

# annotate at least x mapped reads in gene
for at_least_x_mapped_reads in at_least_x_mapped_reads_per_gene_options:
    mapped_reads_per_gene_flag_df = (
        alignment_stats_df.loc[:, ["Chrom", "MappedReads"]]
        .assign(**{f"AtLeast{at_least_x_mapped_reads}MappedReadsPerGene": lambda x: x["MappedReads"] >= at_least_x_mapped_reads})
        .loc[:, ["Chrom", f"AtLeast{at_least_x_mapped_reads}MappedReadsPerGene"]]
    )
    enriched_snps_df = enriched_snps_df.merge(
        mapped_reads_per_gene_flag_df,
        on="Chrom",
        how="left",
    )

enriched_snps_df

# %%
enriched_snps_df["Chrom"].value_counts().describe().round(2)

# %%
enriched_snps_df.drop_duplicates("Chrom")["SNPsPerGene"].describe().round(2)

# %%
fig = px.histogram(
    enriched_snps_df.drop_duplicates("Chrom"),
    x="SNPsPerGene",
    # facet_col="Mismatch",
    # nbins=30,
    # log_y=True,
    # template=template,
    # title="Mismatch frequency distribution of enriched SNPs<br>in transcripts with pooled noise level < 6%",
)
fig.update_xaxes(dtick=1)
fig.update_layout(width=600, height=400, showlegend=False)
fig.show()

# %%
alignment_stats_df.loc[
    alignment_stats_df["MappedReads"].ge(200), 
    ["Chrom", "MappedReads",]
]


# %%
def find_highest_val_x_is_ge(x, vals):
    for val in sorted(vals, reverse=True):
        if x >= val:
            return val
    return np.nan

def find_lowest_val_x_is_le(x, vals):
    for val in sorted(vals):
        if x <= val:
            return val
    return np.nan


# %%
enriched_snps_df

# %%
# enriched_snps_per_gene_df = (
#     enriched_snps_df.drop_duplicates("Chrom").loc[:, ["Chrom", "SNPsPerGene"]]
#     .merge(
#         alignment_stats_df.loc[:, ["Chrom", "MappedReads"]],
#         on="Chrom",
#         how="inner",
#     )
# )

# # at_least_x_snps_per_gene_options = [6, 8, 10]
# # at_most_x_snps_per_gene_options = [16, 18, 20]
# at_least_x_snps_per_gene_options = [3]
# at_most_x_snps_per_gene_options = [3]
# at_least_x_mapped_reads_per_gene_options = [200, 600, 1000]

# enriched_snps_per_gene_df["AtLeastXSNPs"] = enriched_snps_per_gene_df["SNPs"].apply(
#     lambda x: find_highest_val_x_is_ge(x, at_least_x_snps_per_gene_options)
# )
# enriched_snps_per_gene_df["AtMostXSNPs"] = enriched_snps_per_gene_df["SNPs"].apply(
#     lambda x: find_lowest_val_x_is_le(x, at_most_x_snps_per_gene_options)
# )
# enriched_snps_per_gene_df["AtLeastXMappedReads"] = enriched_snps_per_gene_df["MappedReads"].apply(
#     lambda x: find_highest_val_x_is_ge(x, at_least_x_mapped_reads_per_gene_options)
# )

# enriched_snps_per_gene_df = enriched_snps_per_gene_df.loc[
#     (enriched_snps_per_gene_df["AtLeastXSNPs"].eq(min(at_least_x_snps_per_gene_options)))
#     & (enriched_snps_per_gene_df["AtLeastXMappedReads"].eq(min(at_least_x_mapped_reads_per_gene_options)))
# ]

# enriched_snps_per_gene_df

# %%
# enriched_snps_per_gene_df = (
#     enriched_snps_df.drop_duplicates("Chrom").loc[:, ["Chrom", "SNPsPerGene"]]
#     .merge(
#         alignment_stats_df.loc[:, ["Chrom", "MappedReads"]],
#         on="Chrom",
#         how="inner",
#     )
# )

enriched_snps_per_gene_df = (
    enriched_snps_df.drop_duplicates("Chrom").loc[
        enriched_snps_df["SNPsPerGene"].eq(enriched_snps_df["SNPsPerGene"].max()), 
        ["Chrom", "SNPsPerGene"]
    ]
    .rename(columns={"SNPsPerGene": "SNPs"})
    .merge(
        alignment_stats_df.loc[:, ["Chrom", "MappedReads"]],
        on="Chrom",
        how="inner",
    )
)

# at_least_x_snps_per_gene_options = [6, 8, 10]
# at_most_x_snps_per_gene_options = [16, 18, 20]
at_least_x_mapped_reads_per_gene_options = [200, 600, 1000]

# enriched_snps_per_gene_df["AtLeastXSNPs"] = enriched_snps_per_gene_df["SNPs"].apply(
#     lambda x: find_highest_val_x_is_ge(x, at_least_x_snps_per_gene_options)
# )
# enriched_snps_per_gene_df["AtMostXSNPs"] = enriched_snps_per_gene_df["SNPs"].apply(
#     lambda x: find_lowest_val_x_is_le(x, at_most_x_snps_per_gene_options)
# )
enriched_snps_per_gene_df["AtLeastXMappedReads"] = enriched_snps_per_gene_df["MappedReads"].apply(
    lambda x: find_highest_val_x_is_ge(x, at_least_x_mapped_reads_per_gene_options)
)

enriched_snps_per_gene_df = enriched_snps_per_gene_df.loc[
    (enriched_snps_per_gene_df["AtLeastXMappedReads"].eq(min(at_least_x_mapped_reads_per_gene_options)))
    # & (enriched_snps_per_gene_df["AtLeastXSNPs"].eq(min(at_least_x_snps_per_gene_options)))
]

enriched_snps_per_gene_df

# %%
enriched_snps_per_gene_df["SNPs"].value_counts()

# %%
enriched_snps_per_gene_df["AtLeastXSNPs"].value_counts()

# %%
enriched_snps_per_gene_df["AtLeastXMappedReads"].value_counts()

# %%
enriched_snps_per_gene_df["MappedReads"].describe().round(2)

# %%
enriched_snps_df.merge(
    enriched_snps_per_gene_df.loc[:, ["Chrom"]],
    on="Chrom",
    how="inner"
).groupby("Chrom")["MismatchFrequency"].agg(["mean", "std"])

# %%
fig = px.histogram(
    enriched_snps_df.merge(
        enriched_snps_per_gene_df.loc[:, ["Chrom"]],
        on="Chrom",
        how="inner"
    ),
    # x="MappedReads",
    x="MismatchFrequency",
    facet_col="Chrom",
    facet_col_wrap=10,
    facet_col_spacing=0.1,
    nbins=10,
    # log_y=True,
    # template=template,
    # title="Mismatch frequency distribution of enriched SNPs<br>in transcripts with pooled noise level < 6%",
)
fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
fig.update_xaxes(dtick=0.1)
fig.update_yaxes(dtick=1)
fig.update_layout(
    # width=1800, 
    # height=1000, 
    width=max(600, 280 * enriched_snps_per_gene_df.shape[0]), 
    height=400, 
    # showlegend=False
)
fig.show()

# %%

# %%

# %% [markdown]
# ### Total mismatches

# %%
mismatches_dfs = []
for mismatches_file in mismatches_files:
    if Path(mismatches_file).exists():
        mismatches_df = pd.read_csv(mismatches_file, sep="\t")
        mismatches_dfs.append(mismatches_df)
mismatches_df = pd.concat(
    mismatches_dfs,
    ignore_index=True
)
del mismatches_dfs

mismatches_df

# %%
concat_edited_positions_df = concat_all_positions_df.loc[
    concat_all_positions_df["EditedFinal"]
]
# concat_edited_positions_df.loc[:, ["Chrom"]]

# %%
# avg editing positions per transcript, considering transcripts whose pooled noise levels is < 6%
num_edited_positions_per_chrom_df = (
    concat_all_positions_df.loc[
        (concat_all_positions_df["EditedFinal"])
        # & (concat_all_positions_df["Chrom"].isin(chroms))
    ]
    .groupby("Chrom")
    .size()
    .reset_index(name="NumOfEditingSites")
)
num_edited_positions_per_chrom_df

# %%
alignment_stats_df

# %%
# Build two summary dfs from mismatches_df:
# 1) AG vs all other mismatches combined
# 2) AG vs the highest-count other mismatch (per Chrom)

# -- if "mismatches_df" not in globals():
# --     raise NameError("mismatches_df is not defined. Run the cell that creates mismatches_df first.")

req_cols = {"Chrom", "Mismatch", "Count"}
missing = req_cols - set(mismatches_df.columns)
if missing:
    raise ValueError(f"mismatches_df is missing required columns: {sorted(missing)}")

mm = mismatches_df.copy()
mm["Count"] = pd.to_numeric(mm["Count"], errors="coerce").fillna(0).astype(int)

# --- 1) AG vs all others together ---
ag_vs_all_others_df = (
    mm.assign(IsAG=mm["Mismatch"].eq("AG"))
      .groupby(["Chrom", "IsAG"], as_index=False)["Count"].sum()
      .pivot(index="Chrom", columns="IsAG", values="Count")
      .rename(columns={True: "AG", False: "AllOthers"})
      .fillna(0)
      .astype(int)
      .reset_index()
)
ag_vs_all_others_df["Total"] = ag_vs_all_others_df["AG"] + ag_vs_all_others_df["AllOthers"]
ag_vs_all_others_df["AG_Frac"] = np.where(
    ag_vs_all_others_df["Total"] > 0,
    ag_vs_all_others_df["AG"] / ag_vs_all_others_df["Total"],
    np.nan,
)

ag_vs_all_others_df = ag_vs_all_others_df.merge(
    per_transcript_editing_index_df,
    how="outer",
    on="Chrom",
)

ag_vs_all_others_df = ag_vs_all_others_df.merge(
    num_edited_positions_per_chrom_df,
    how="outer",
    on="Chrom",
)
ag_vs_all_others_df["NumOfEditingSites"] = ag_vs_all_others_df["NumOfEditingSites"].fillna(0).astype(int)

ag_vs_all_others_df = ag_vs_all_others_df.merge(
    alignment_stats_df.loc[:, ["Chrom", "MappedReads", "MappedReadsPerSample"]],
    how="left",
    on="Chrom",
)

# --- 2) AG vs highest other mismatch ---
ag_df = mm.loc[mm["Mismatch"].eq("AG"), ["Chrom", "Count"]].rename(columns={"Count": "AG"})

non_ag = mm.loc[~mm["Mismatch"].eq("AG")].copy()
highest_other = (
    non_ag.sort_values(["Chrom", "Count", "Mismatch"], ascending=[True, False, True])
          .groupby("Chrom", as_index=False)
          .first()
          .rename(columns={"Mismatch": "TopOtherMismatch", "Count": "TopOtherCount"})
)

ag_vs_top_other_df = (
    ag_df.merge(highest_other, on="Chrom", how="outer")
         .fillna({"AG": 0, "TopOtherCount": 0, "TopOtherMismatch": ""})
)
ag_vs_top_other_df["AG"] = ag_vs_top_other_df["AG"].astype(int)
ag_vs_top_other_df["TopOtherCount"] = ag_vs_top_other_df["TopOtherCount"].astype(int)
ag_vs_top_other_df["Total_AG_TopOther"] = ag_vs_top_other_df["AG"] + ag_vs_top_other_df["TopOtherCount"]
ag_vs_top_other_df["AG_Frac_AG_TopOther"] = np.where(
    ag_vs_top_other_df["Total_AG_TopOther"] > 0,
    ag_vs_top_other_df["AG"] / ag_vs_top_other_df["Total_AG_TopOther"],
    np.nan,
)

ag_vs_top_other_df = ag_vs_top_other_df.merge(
    per_transcript_editing_index_df,
    how="outer",
    on="Chrom",
)

ag_vs_top_other_df = ag_vs_top_other_df.merge(
    num_edited_positions_per_chrom_df,
    how="outer",
    on="Chrom",
)
ag_vs_top_other_df["NumOfEditingSites"] = ag_vs_top_other_df["NumOfEditingSites"].fillna(0).astype(int)

ag_vs_top_other_df = ag_vs_top_other_df.merge(
    alignment_stats_df.loc[:, ["Chrom", "MappedReads", "MappedReadsPerSample"]],
    how="left",
    on="Chrom",
)

# ag_vs_all_others_df, ag_vs_top_other_df

# %%
per_transcript_editing_index_df

# %%
ag_vs_all_others_df

# %%
ag_vs_top_other_df

# %%
ag_vs_top_other_df["AG"] / ag_vs_top_other_df["TopOtherCount"]

# %%
fig = px.histogram(
    ag_vs_top_other_df,
    x="AG_Frac_AG_TopOther",
    labels={"AG_Frac_AG_TopOther": "AG mismatch / (AG mismatch + top other mismatch)"},
)
fig.update_layout(
    width=600,
    height=400,
)
fig.show()

# %%
fig = px.histogram(
    ag_vs_all_others_df,
    x="AG_Frac",
    labels={"AG_Frac": "AG mismatch / all 12 mismatches"},
)
fig.update_layout(
    width=600,
    height=400,
)
fig.show()

# %%
fig = px.density_heatmap(
    ag_vs_top_other_df,
    x="EditingIndex",
    y="AG_Frac_AG_TopOther",
    labels={
        "AG_Frac_AG_TopOther": "AG mismatch / (AG mismatch + top other mismatch)",
        "EditingIndex": "Editing index",
    },
    # color_continuous_scale='YlGnBu',
    # color_continuous_scale=px.colors.sequential.YlGnBu, 
    marginal_x="histogram",
    marginal_y="histogram"
)
fig.update_xaxes(dtick=1)
fig.update_layout(
    width=600,
    height=550,
)
fig.show()

# %%
fig = px.density_heatmap(
    ag_vs_all_others_df,
    x="EditingIndex",
    y="AG_Frac",
    labels={
        "AG_Frac": "AG mismatch / all 12 mismatches",
        "EditingIndex": "Editing index",
    },
    # color_continuous_scale='YlGnBu',
    marginal_x="histogram",
    marginal_y="histogram"
)
fig.update_xaxes(dtick=1)
fig.update_layout(
    width=600,
    height=550,
)
fig.show()

# %%
fig = px.scatter(
    ag_vs_top_other_df,
    x="EditingIndex",
    y="AG_Frac_AG_TopOther",
    labels={
        "AG_Frac_AG_TopOther": "AG mismatch / (AG mismatch + top other mismatch)",
        "EditingIndex": "Editing index",
        "NumOfEditingSites": "Editing sites",
    },
    trendline="ols",
    trendline_color_override="black",
    marginal_x="histogram", 
    marginal_y="histogram",
    color="NumOfEditingSites",
    # color_continuous_scale=px.colors.sequential.YlGnBu,
    color_continuous_scale=px.colors.sequential.Turbo,
    # color_continuous_scale=px.colors.sequential.Blackbody,
    # color_continuous_scale=px.colors.sequential.Bluered,
    # color_continuous_scale=px.colors.sequential.RdBu,
)

tr_line=[]
for k, trace in enumerate(fig.data):
    # ic(k, trace)
    try:
        if trace.mode is not None and trace.mode == 'lines':
            tr_line.append(k)
    except AttributeError:
        pass # for Histogram which has no attribute 'mode' 
for id in tr_line:
    fig.data[id].update(line_width=4)

fig.update_xaxes(dtick=1)
fig.update_layout(
    width=600,
    height=550,
)

# results = px.get_trendline_results(fig)
# print(results)

# fig.show()

# %%
fig = px.scatter(
    ag_vs_all_others_df,
    x="EditingIndex",
    y="AG_Frac",
    labels={
        "AG_Frac": "AG mismatch / all 12 mismatches",
        "EditingIndex": "Editing index",
        "NumOfEditingSites": "Editing sites",
    },
    trendline="ols",
    trendline_color_override="black",
    marginal_x="histogram", 
    marginal_y="histogram",
    color="NumOfEditingSites",
    # color_continuous_scale=px.colors.sequential.YlGnBu,
    color_continuous_scale=px.colors.sequential.Turbo,
    # color_continuous_scale=px.colors.sequential.Blackbody,
    # color_continuous_scale=px.colors.sequential.Bluered,
    # color_continuous_scale=px.colors.sequential.RdBu,
)

tr_line=[]
for k, trace in enumerate(fig.data):
    # ic(k, trace)
    try:
        if trace.mode is not None and trace.mode == 'lines':
            tr_line.append(k)
    except AttributeError:
        pass # for Histogram which has no attribute 'mode' 
for id in tr_line:
    fig.data[id].update(line_width=4)

fig.update_xaxes(dtick=1)
fig.update_layout(
    width=600,
    height=550,
)
fig.show()

# %%
fig = px.histogram(
    ag_vs_top_other_df,
    x="EditingIndex",
    y="AG_Frac_AG_TopOther",
    labels={
        "AG_Frac_AG_TopOther": "AG mismatch / (AG mismatch + top other mismatch)",
        "EditingIndex": "Editing index",
    },
    histfunc="avg",
)
fig.update_xaxes(dtick=1)
fig.update_layout(
    width=600,
    height=550,
)
fig.show()

# %%
fig = px.histogram(
    ag_vs_all_others_df,
    x="EditingIndex",
    y="AG_Frac",
    labels={
        "AG_Frac": "AG mismatch / all 12 mismatches",
        "EditingIndex": "Editing index",
    },
    histfunc="avg",
)
fig.update_xaxes(dtick=1)
fig.update_layout(
    width=600,
    height=550,
)
fig.show()

# %%
fig = px.histogram(
    ag_vs_top_other_df,
    x="NumOfEditingSites",
    y="AG_Frac_AG_TopOther",
    labels={
        "AG_Frac_AG_TopOther": "AG mismatch / (AG mismatch + top other mismatch)",
        # "EditingIndex": "Editing index",
        "NumOfEditingSites": "Editing sites",
    },
    histfunc="avg",
)
# fig.update_xaxes(dtick=1)
fig.update_layout(
    width=600,
    height=550,
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
all_pooled_per_chrom_machine_noise_df[
    "%PooledMachineNoise"
] = all_pooled_per_chrom_machine_noise_df.apply(
    lambda x: 100 * x["Mismatches"] / x["Matches"], axis=1
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
# ### Noise sites above noise threshold

# %%
noise_threshold_files_df

# %%
noise_threshold_df = pd.concat(
    [
        pd.read_csv(
            noise_threshold_file, 
            sep="\t",
            names=["Chrom", "NoiseThreshold"]
        )
        for noise_threshold_file in noise_threshold_files_df["NoiseThresholdFile"].values
    ],
    ignore_index=True
)
noise_threshold_df

# %%
concat_all_positions_df

# %%
noise_positions_df = concat_all_positions_df.loc[
    concat_all_positions_df["NoisyFinal"]
]

noise_positions_df = noise_positions_df.merge(
    noise_threshold_df,
    on="Chrom",
    how="left"
)

noise_positions_df["NoiseAboveNoiseThreshold"] = (
    noise_positions_df["Noise"].ge(noise_positions_df["NoiseThreshold"])
)

noise_positions_df

# %%
noise_positions_df.loc[
    noise_positions_df["NoiseAboveNoiseThreshold"]
]

# %%
noise_positions_df.loc[
    (noise_positions_df["NoiseAboveNoiseThreshold"])
    & (noise_positions_df["Noise"] == 1)
]

# %%
concat_edited_positions_df = concat_all_positions_df.loc[
    concat_all_positions_df["EditedFinal"]
]
concat_edited_positions_df

# %%
fig = px.histogram(
    concat_edited_positions_df,
    x="EditingFrequency",
    histnorm="percent",
)
fig.update_layout(
    width=600,
    height=450,
    title="Editing freq. of final editing sites"
)
fig.show()

# %%
fig = px.histogram(
    noise_positions_df.loc[
        noise_positions_df["NoiseAboveNoiseThreshold"]
    ],
    x="Noise",
    histnorm="percent",
)
fig.update_layout(
    width=600,
    height=450,
    title="Noise freq. of noise sites >= editing threshold"
)
fig.show()

# %% [markdown]
# ### Known & new editing sites

# %%
cols = 1
rows = 1

fig, ax = plt.subplots(
    # figsize=(3.2 * cols, 2.5 * rows),
    figsize=(5 * cols, 2.5 * rows),
    constrained_layout=True,
    facecolor="white",
)

ax.set_facecolor("white")

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

labels[0] = f"De-novo\n({len(sets[0])})"
labels[1] = f"Known\n({len(sets[1])})"

venn = venn2(sets, set_labels=labels, ax=ax)

# Make all Venn text black
for text in ax.texts:
    text.set_color("black")

fig.suptitle(
    "Whole-transcriptome octopus data",
    fontsize="xx-large",
    color="black",
)

plt.savefig(
    Path(out_dir, "Known vs new editing sites - Octopus - pooled.svg"),
    format="svg",
    dpi=300,
    facecolor="white",
)

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
main_title = "Whole-transcriptome octopus data"
sub_titles = [""]

out_file = Path(out_dir, "ADAR motif of pooled editing sites - Octopus - pooled.svg")

# %%
multiple_logos_from_fasta_files(
    fasta_files, main_title, sub_titles, out_file, 
    # width=0.33 * 14, 
    width=0.33 * 14, 
    height=4, 
    dpi=300,
    # tighthen_layout=True
);


# %% [markdown]
# ## Editing by haplotype

# %% [markdown]
# ### Chi-square permutations enrichment functions

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
            # print(f"Skipping site {site_col} due to low total reads ({contingency_table.sum()})")
            print(f"Skipping site {site_col} due to low reads in one or more cells ({contingency_table})")
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
def permutation_null(
    df, 
    haplotypes, 
    n_perm=200, 
    alpha=0.05, 
    min_reads_per_cell=5, 
    seed=seed, 
    hap_col="Haplotype",
    test=_chi_square_two_sided
):
    # set up random generator with fixed seed for reproducibility
    rng = np.random.default_rng(seed)
    # the editing columns are all columns after the haplotype column
    editing_cols = df.iloc[:, df.columns.get_loc(hap_col) + 1:].columns
    # copy relevant subset to allow shuffling haplotype labels w/o modifying original df
    sub = df.loc[df[hap_col].isin(haplotypes)].copy()
    # calculate observed number of significant sites in the original data
    observed, m_tested, per_site = count_diff_sites(sub, hap_col, haplotypes, editing_cols, alpha=alpha, min_reads_per_cell=min_reads_per_cell, test=test)

    null_counts = []
    hap_values = sub[hap_col].to_numpy()
    for _ in range(n_perm):
        # Randomly shuffle (“permute”) the haplotype labels across reads
        perm = hap_values.copy()
        rng.shuffle(perm)
        sub[hap_col] = perm
        # Recalculate the number of significant sites for this permutation
        c, _, _ = count_diff_sites(sub, hap_col, haplotypes, editing_cols, alpha=alpha, min_reads_per_cell=min_reads_per_cell, test=test)
        null_counts.append(c)

    null_counts = np.asarray(null_counts)
    # empirical p-value
    p_emp = (1 + (null_counts >= observed).sum()) / (1 + n_perm)
    return observed, null_counts, p_emp, m_tested, per_site



# %%
def all_haplotypes_permutation_null(
    editing_and_snps_reads_df,
    n_perm=200, 
    alpha=0.05, 
    seed=seed,
    # min_reads_per_cell=5, 
    # test=_chi_square_two_sided
    min_reads_per_cell=0, 
    test=_fisher_two_sided,
    hap_col="Haplotype"
):
    chrom = editing_and_snps_reads_df["Chrom"].iloc[0]
    
    hap_counts = editing_and_snps_reads_df[hap_col].value_counts()
    haplotypes_couples = list(combinations(hap_counts.index, 2))

    results_dict = defaultdict(list)

    for h1, h2 in haplotypes_couples:
        haplotypes = [h1, h2]
        # print(haplotypes)
        # print("Comparing:", h1, "vs", h2)
        # print("Reads per haplotype:\n", hap_counts.loc[haplotypes])
        obs, null_counts, p_emp, m_tested, _ = permutation_null(
            editing_and_snps_reads_df, 
            haplotypes,
            n_perm=n_perm, 
            alpha=alpha,
            min_reads_per_cell=min_reads_per_cell,  
            seed=seed,
            hap_col=hap_col,
            test=test
        )
        # print(f"Tested sites: {m_tested}")
        # print(f"Observed #diff sites (p<0.05): {obs}")
        # print(f"Null mean: {null_counts.mean():.2f}, null 95%: [{np.quantile(null_counts, 0.025)}, {np.quantile(null_counts, 0.975)}]")
        # print(f"Empirical p-value: {p_emp:.4g}")
        null_counts_ge_observed = (null_counts >= obs).sum()
        results_dict["Haplotypes"].append(haplotypes)
        results_dict["Tested sites"].append(m_tested)
        results_dict["Observed"].append(obs)
        results_dict["Null >= Observed"].append(null_counts_ge_observed)
        # results_dict["Null Mean"].append(null_counts.mean())
        # results_dict["Null 95%"].append([np.quantile(null_counts, 0.025), np.quantile(null_counts, 0.975)])
        results_dict["Empirical p-value"].append(p_emp)

    results_df = pd.DataFrame(results_dict)
    results_df.insert(0, "Chrom", chrom)
    
    # also correct for multiple testing across the different haplotype pairs (e.g. 3 pairs for 3 haplotypes) 
    # using Benjamini-Hochberg FDR correction
    reject, pvals_corr, *_ = multipletests(results_df["Empirical p-value"].tolist(), method='fdr_bh')
    results_df["Empirical p-value (FDR)"] = pvals_corr
    results_df["Reject H0 (FDR)"] = reject
    
    return results_df


# %%
# editing_and_snps_reads_df

# %%
# n_perm=200
# alpha=0.05
# seed=seed
# # min_reads_per_cell=5, 
# # test=_chi_square_two_sided
# min_reads_per_cell=0
# test=_fisher_two_sided
# hap_col="Haplotype"

# %%
# chrom = editing_and_snps_reads_df["Chrom"].iloc[0]
# chrom

# %%
# hap_counts = editing_and_snps_reads_df[hap_col].value_counts()
# hap_counts

# %%
# haplotypes_couples = list(combinations(hap_counts.index, 2))
# haplotypes_couples

# %%
# h1, h2 = haplotypes_couples[0]
# haplotypes = [h1, h2]
# h1, h2, haplotypes

# %%
# obs, null_counts, p_emp, m_tested, _ = permutation_null(
#     editing_and_snps_reads_df, 
#     haplotypes,
#     n_perm=n_perm, 
#     alpha=alpha,
#     min_reads_per_cell=min_reads_per_cell,  
#     seed=seed,
#     hap_col=hap_col,
#     test=test
# )
# obs, null_counts, p_emp, m_tested

# %%
# # work on df copy
# df = editing_and_snps_reads_df.copy()
# # df

# %%
# # set up random generator with fixed seed for reproducibility
# rng = np.random.default_rng(seed)

# # the editing columns are all columns after the haplotype column
# editing_cols = df.iloc[:, df.columns.get_loc(hap_col) + 1:].columns
# editing_cols

# %%
# # # copy relevant subset to allow shuffling haplotype labels w/o modifying original df
# sub = df.loc[df[hap_col].isin(haplotypes)].copy()
# sub

# %%
# # calculate observed number of significant sites in the original data
# observed, m_tested, per_site = count_diff_sites(
#     sub, hap_col, haplotypes, editing_cols, 
#     alpha=alpha, min_reads_per_cell=min_reads_per_cell, test=test
# )

# null_counts = []
# hap_values = sub[hap_col].to_numpy()
# for _ in range(n_perm):
#     # Randomly shuffle (“permute”) the haplotype labels across reads
#     perm = hap_values.copy()
#     rng.shuffle(perm)
#     sub[hap_col] = perm
#     # Recalculate the number of significant sites for this permutation
#     c, _, _ = count_diff_sites(sub, hap_col, haplotypes, editing_cols, alpha=alpha, min_reads_per_cell=min_reads_per_cell, test=test)
#     null_counts.append(c)

# null_counts = np.asarray(null_counts)
# # empirical p-value
# p_emp = (1 + (null_counts >= observed).sum()) / (1 + n_perm)

# %%

# %% [markdown]
# ### Selected diverse haplotypes tests

# %%
def validate_haplotypes_diversity(
    editing_and_snps_reads_df,
    min_haplotypes,
    min_reads_per_haplotype,
    min_editing_sites=0
):
    
    num_editing_sites = len(editing_and_snps_reads_df.loc[:, "Haplotype":].columns) - 1
    enough_editing_sites = num_editing_sites >= min_editing_sites
    
    haplotypes_value_counts = editing_and_snps_reads_df["Haplotype"].value_counts()
    enough_well_covered_haplotypes = haplotypes_value_counts.ge(min_reads_per_haplotype).sum() >= min_haplotypes
    
    if enough_editing_sites and enough_well_covered_haplotypes:
        return True
    return False


# %%
# editing_and_snps_reads_df = editing_and_snps_reads_dfs[2]
# editing_and_snps_reads_df

# %%

# %%
diverse_haplotypes_editing_and_snps_reads_dfs = [
    df
    for df in editing_and_snps_reads_dfs
    # if validate_haplotypes_diversity(df, 2, 50)
    # if validate_haplotypes_diversity(df, 2, 300) or validate_haplotypes_diversity(df, 3, 100)
    # if validate_haplotypes_diversity(df, 3, 100)
    # if validate_haplotypes_diversity(df, 3, 70, 5)
    if validate_haplotypes_diversity(df, 2, 200, 15)
]
ic(len(diverse_haplotypes_editing_and_snps_reads_dfs));

# %%
diverse_haplotypes_editing_and_snps_reads_dfs[0]

# %%
editing_long_dfs = []
editing_profiles_dfs = []

for editing_and_snps_reads_df in diverse_haplotypes_editing_and_snps_reads_dfs:
    
    chrom = editing_and_snps_reads_df["Chrom"].iloc[0]
    editing_cols = editing_and_snps_reads_df.loc[:, "Haplotype":].iloc[:, 1:].columns
    
    editing_long_df = (
        editing_and_snps_reads_df.groupby("Haplotype")[editing_cols]
        .agg(["mean", "std"])
        .stack(level=0, future_stack=True)  # silence FutureWarning (pandas>=2.1)
        .reset_index()
        .rename(columns={"level_1": "EditingSite"})
    )
    # # optional: nicer dtypes / sorting
    editing_long_df["EditingSite"] = editing_long_df["EditingSite"].astype(int)
    editing_long_df = editing_long_df.sort_values(["EditingSite", "Haplotype"], ignore_index=True)
    editing_long_df["EditingSite"] = editing_long_df["EditingSite"].astype(str)
    editing_long_df.insert(0, "Chrom", chrom)
    editing_long_dfs.append(editing_long_df)
    
    editing_profiles_df = editing_and_snps_reads_df.groupby("Haplotype")[editing_cols].mean().mul(100).T
    editing_profiles_df.dropna(how="all", inplace=True)
    editing_profiles_df.insert(0, "Chrom", chrom)
    editing_profiles_dfs.append(editing_profiles_df)
    
concat_editing_long_df = pd.concat(editing_long_dfs)
concat_editing_profiles_dfs = pd.concat(editing_profiles_dfs)

concat_editing_long_df["EditingSite"] = concat_editing_long_df["EditingSite"].astype(int)
concat_editing_long_df["%mean"] = concat_editing_long_df["mean"].mul(100)

# del editing_long_dfs, editing_profiles_dfs
# del editing_long_dfs

# %%
concat_editing_long_df 

# %%
concat_editing_profiles_dfs

# %%
concat_editing_profiles_dfs.iloc[:, 1:].max().max()

# %%
editing_profiles_dfs[0]


# %%
def make_haplotype_diffs_df(editing_profiles_df):
    haplotype_diffs_df = editing_profiles_df.copy()

    haplotypes = haplotype_diffs_df.columns[1:]

    for h1, h2 in combinations(haplotypes, 2):
        haplotype_diffs_df[f"{h1} - {h2}"] = haplotype_diffs_df[h1] - haplotype_diffs_df[h2]
        
    haplotype_diffs_df = haplotype_diffs_df.drop(columns=haplotypes)
    
    return haplotype_diffs_df


# %%
haplotype_diffs_dfs = [
    make_haplotype_diffs_df(editing_profiles_df)
    for editing_profiles_df in editing_profiles_dfs
]
haplotype_diffs_dfs[0]

# %%
# -----------------------
# Chrom -> subplot title
# -----------------------
chrom_to_title = (
    orfs_df.drop_duplicates("Chrom")
    .set_index("Chrom")["Name"]
    .to_dict()
)

def chrom_to_plot_title(chrom):
    name = chrom_to_title.get(chrom, chrom)
    return str(name).split("_")[0]


# %%

df = concat_editing_long_df  # expects columns: Chrom, EditingSite, %mean, Haplotype

facet_col = "Chrom"
x_col = "EditingSite"
y_col = "%mean"
color_col = "Haplotype"

facet_col_wrap = 3
opacity = 0.8

width = 1000
height = 600

# Global axis title styling
global_axis_font_size = 20
# global_x_title_y = -0.085
# global_y_title_x = -0.075
# global_x_title_y = -0.04
# global_y_title_x = -0.04
global_x_title_y = -0.055
global_y_title_x = -0.040


# ---------------------------
# 1) Decide facet order + grid
# ---------------------------
chroms = list(pd.unique(df[facet_col]))
n_panels = len(chroms)
ncols = facet_col_wrap
nrows = math.ceil(n_panels / ncols)

# -----------------------------------------
# 2) Stable color map for haplotypes (global)
# -----------------------------------------
haplotypes = list(pd.unique(df[color_col]))

palette = px.colors.qualitative.D3
if len(haplotypes) > len(palette):
    palette = (palette * (len(haplotypes) // len(palette) + 1))[: len(haplotypes)]

color_map = {h: c for h, c in zip(haplotypes, palette)}

# ---------------------------
# 3) Create subplots
#    Titles should be mapped gene names (not Chrom)
# ---------------------------
subplot_titles = [chrom_to_plot_title(c) for c in chroms]

fig = make_subplots(
    rows=nrows,
    cols=ncols,
    subplot_titles=subplot_titles,
    horizontal_spacing=0.06,
    vertical_spacing=0.14,
)

panel_legends = {}  # (row, col) -> list[(haplotype, color)]

# ---------------------------
# 4) Add traces per panel
# ---------------------------
for i, chrom in enumerate(chroms):
    r = i // ncols + 1
    c = i % ncols + 1

    sub = df[df[facet_col] == chrom]

    present_haps = list(pd.unique(sub[color_col]))
    panel_legends[(r, c)] = [(h, color_map[h]) for h in present_haps]

    for h in present_haps:
        sh = sub[sub[color_col] == h]
        fig.add_trace(
            go.Scattergl(
                x=sh[x_col],
                y=sh[y_col],
                mode="markers",
                marker=dict(color=color_map[h], size=7, opacity=opacity),
                name=str(h),
                showlegend=False,  # legend via annotations
            ),
            row=r,
            col=c,
        )

    # # X axis label only on bottom row
    # fig.update_xaxes(
    #     title_text=("Editing site" if r == nrows else None),
    #     row=r,
    #     col=c,
    # )

    # # Y axis label only on left-most column
    # fig.update_yaxes(
    #     title_text=("Mean editing [%]" if c == 1 else None),
    #     row=r,
    #     col=c,
    # )

# ---------------------------
# 5) Layout + styling
# ---------------------------
fig.update_layout(
    width=width,
    height=height,
    # template="simple_white",
    margin=dict(l=80, r=20, t=60, b=60),
)

fig.update_annotations(font=dict(size=12))


# ---------------------------
# 6) Global/shared axis titles (single title per axis)
# ---------------------------
fig.add_annotation(
    x=0.5,
    y=global_x_title_y,
    xref="paper",
    yref="paper",
    text="Editing site",
    showarrow=False,
    xanchor="center",
    yanchor="top",
    font=dict(size=global_axis_font_size),
)

fig.add_annotation(
    x=global_y_title_x,
    y=0.5,
    xref="paper",
    yref="paper",
    text="Mean editing [%]",
    showarrow=False,
    xanchor="right",
    yanchor="middle",
    textangle=-90,
    font=dict(size=global_axis_font_size),
)

# ---------------------------
# 7) "Legend inside each subplot" via annotations
# ---------------------------
layout = fig.layout

for i, chrom in enumerate(chroms):
    r = i // ncols + 1
    c = i % ncols + 1

    axis_index = i + 1
    xaxis_name = "xaxis" if axis_index == 1 else f"xaxis{axis_index}"
    yaxis_name = "yaxis" if axis_index == 1 else f"yaxis{axis_index}"

    xdom = getattr(layout, xaxis_name).domain
    ydom = getattr(layout, yaxis_name).domain

    items = panel_legends[(r, c)]
    legend_lines = ["<span style='font-size:12px;'><b>Haplotype</b></span>"]
    for h, colr in items:
        legend_lines.append(
            f"<span style='color:{colr};'>●</span> <span style='font-size:11px;'>{h}</span>"
        )
    legend_html = "<br>".join(legend_lines)

    fig.add_annotation(
        x=xdom[1] - 0.01,
        y=ydom[1] - 0.01,
        xref="paper",
        yref="paper",
        xanchor="right",
        yanchor="top",
        text=legend_html,
        showarrow=False,
        align="left",
        bgcolor="rgba(255,255,255,0.7)",
        bordercolor="rgba(0,0,0,0.15)",
        borderwidth=1,
        borderpad=4,
    )
    


fig.show()

# %%
len(diverse_haplotypes_editing_and_snps_reads_dfs)

# %%
with Pool(processes=6) as pool:
    selected_examples_fisher_haplotypes_dfs = pool.map(
        all_haplotypes_permutation_null,
        diverse_haplotypes_editing_and_snps_reads_dfs
    )
concat_selected_examples_fisher_haplotypes_df = pd.concat(selected_examples_fisher_haplotypes_dfs, ignore_index=True)
concat_selected_examples_fisher_haplotypes_df.insert(
    1, 
    condition_col,
    concat_selected_examples_fisher_haplotypes_df["Chrom"].apply(chrom_to_plot_title)
)
concat_selected_examples_fisher_haplotypes_df

# %%
concat_selected_examples_fisher_haplotypes_df.round(5)

# %%
editing_and_snps_reads_df = diverse_haplotypes_editing_and_snps_reads_dfs[0].copy()
editing_and_snps_reads_df

# %%
# COMPLETE CODE REPLACEMENT (Plotly) — fixes:
# 1) y=x legend entry: shows a dashed sample in the legend (like the line), while the plotted line stays normal dash
# 2) axis titles: truly bigger (and we avoid shrinking them via update_annotations at the end)

import math
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from itertools import combinations
import ast

# -----------------------
# Config
# -----------------------
facet_col_wrap = 3
width = 1100
height = 820

axis_max = 100
dtick = 10
x_grid = np.linspace(0, axis_max, 200)

show_points = True
points_opacity = 0.5
marker_size = 4
band_alpha = 0.15
line_width = 2

# y=x line in plot (keep normal dashed)
base_line_width = line_width + 1
base_line_dash = "dash"   # normal dashed

# Axis titles (make them obviously bigger)
# global_axis_font_size = 38   # bumped up so it's unmistakably larger
global_axis_font_size = 20   # bumped up so it's unmistakably larger
global_x_title_y = -0.055
global_y_title_x = -0.040

# spacing
horizontal_spacing = 0.035
vertical_spacing = 0.105

pair_palette = [
    "#EF553B", "#00CC96", "#AB63FA", "#FFA15A", "#19D3F3",
    "#FF6692", "#B6E880", "#FF97FF", "#FECB52", "#636EFA"
]

# Legend placement: left edge corresponds to x=3 (data units)
legend_data_x = 3.0
legend_top_inset = 0.002

# Legend fonts (smaller as requested)
# panel_header_size = 11
panel_header_size = 12
panel_item_size = 9
panel_section_size = 10
panel_footer_size = 9

# ---------- helpers ----------
def ols_line_and_ci(x, y, x_grid, alpha=0.05):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]

    n = x.size
    if n < 3 or np.allclose(x.var(), 0):
        return None

    X = np.column_stack([np.ones(n), x])
    XtX = X.T @ X
    beta = np.linalg.solve(XtX, X.T @ y)

    y_fit = X @ beta
    resid = y - y_fit
    dof = n - 2
    s2 = (resid @ resid) / dof
    cov_beta = s2 * np.linalg.inv(XtX)

    Xg = np.column_stack([np.ones_like(x_grid), x_grid])
    y_hat = Xg @ beta

    se_mean = np.sqrt(np.sum(Xg @ cov_beta * Xg, axis=1))
    tcrit = 1.96

    lo = y_hat - tcrit * se_mean
    hi = y_hat + tcrit * se_mean
    return y_hat, lo, hi, beta


def rgba_from_hex(hex_color, alpha):
    hex_color = hex_color.lstrip("#")
    r = int(hex_color[0:2], 16)
    g = int(hex_color[2:4], 16)
    b = int(hex_color[4:6], 16)
    return f"rgba({r},{g},{b},{alpha})"


def parse_hap_pair(v):
    if isinstance(v, (list, tuple)) and len(v) == 2:
        return frozenset(map(str, v))
    if isinstance(v, str):
        try:
            obj = ast.literal_eval(v)
            if isinstance(obj, (list, tuple)) and len(obj) == 2:
                return frozenset(map(str, obj))
        except Exception:
            return None
    return None


def legend_xpaper_for_data_x(layout, axis_index, data_x):
    xaxis_name = "xaxis" if axis_index == 1 else f"xaxis{axis_index}"
    xdom = getattr(layout, xaxis_name).domain
    frac = (data_x - 0.0) / (axis_max - 0.0)
    frac = float(np.clip(frac, 0.0, 1.0))
    return xdom[0] + frac * (xdom[1] - xdom[0])


# -----------------------
# Chrom -> plot title (first token of Name)
# -----------------------
chrom_to_title = (
    orfs_df.drop_duplicates("Chrom")
    .set_index("Chrom")["Name"]
    .to_dict()
)

def chrom_to_plot_title(chrom):
    name = chrom_to_title.get(chrom, chrom)
    return str(name).split("_")[0]


# -----------------------
# Enrichment map
# -----------------------
enrich_df = concat_selected_examples_fisher_haplotypes_df.copy()
reject_col = "Reject H0 (FDR)"
if reject_col not in enrich_df.columns:
    raise KeyError(f"Expected column '{reject_col}' in concat_selected_examples_fisher_haplotypes_df.")

enrich_df["_pairset"] = enrich_df["Haplotypes"].apply(parse_hap_pair)

enriched_pairs_by_chrom = {}
for chrom, sub in enrich_df.groupby("Chrom", dropna=False):
    enriched_pairs_by_chrom[str(chrom)] = set(
        sub.loc[(sub[reject_col] == True) & sub["_pairset"].notna(), "_pairset"].tolist()
    )


# -----------------------
# Prepare panels
# -----------------------
chroms = [df["Chrom"].iloc[0] for df in editing_profiles_dfs]
titles = [chrom_to_plot_title(c) for c in chroms]

n_panels = len(editing_profiles_dfs)
ncols = facet_col_wrap
nrows = math.ceil(n_panels / ncols)

fig = make_subplots(
    rows=nrows,
    cols=ncols,
    subplot_titles=titles,
    horizontal_spacing=horizontal_spacing,
    vertical_spacing=vertical_spacing,
)

panel_legends = {}

# -----------------------
# Add traces
# -----------------------
for i, editing_profiles in enumerate(editing_profiles_dfs):
    chrom = str(editing_profiles["Chrom"].iloc[0])
    hap_cols = list(editing_profiles.iloc[:, 1:].columns)

    r = i // ncols + 1
    c = i % ncols + 1

    # y=x baseline (normal dash)
    fig.add_trace(
        go.Scatter(
            x=[0, axis_max],
            y=[0, axis_max],
            mode="lines",
            line=dict(color="black", width=base_line_width, dash=base_line_dash),
            showlegend=False,
            hoverinfo="skip",
        ),
        row=r,
        col=c,
    )

    enriched_set = enriched_pairs_by_chrom.get(chrom, set())
    enriched_items = []
    not_enriched_items = []

    for k, (h1, h2) in enumerate(combinations(hap_cols, 2)):
        color = pair_palette[k % len(pair_palette)]
        fill_color = rgba_from_hex(color, band_alpha)

        x = pd.to_numeric(editing_profiles[h1], errors="coerce").to_numpy()
        y = pd.to_numeric(editing_profiles[h2], errors="coerce").to_numpy()

        res = ols_line_and_ci(x, y, x_grid)
        if res is None:
            continue

        y_hat, lo, hi, beta = res

        label = f"{h1} vs. {h2}"
        pairset = frozenset([str(h1), str(h2)])
        (enriched_items if pairset in enriched_set else not_enriched_items).append((label, color))

        if show_points:
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=y,
                    mode="markers",
                    marker=dict(color=color, size=marker_size, opacity=points_opacity),
                    showlegend=False,
                    hoverinfo="skip",
                ),
                row=r,
                col=c,
            )

        # CI band
        fig.add_trace(
            go.Scatter(x=x_grid, y=hi, mode="lines", line=dict(width=0),
                       showlegend=False, hoverinfo="skip"),
            row=r, col=c
        )
        fig.add_trace(
            go.Scatter(x=x_grid, y=lo, mode="lines", line=dict(width=0),
                       fill="tonexty", fillcolor=fill_color,
                       showlegend=False, hoverinfo="skip"),
            row=r, col=c
        )

        # regression line
        fig.add_trace(
            go.Scatter(
                x=x_grid,
                y=y_hat,
                mode="lines",
                line=dict(color=color, width=line_width),
                showlegend=False,
            ),
            row=r,
            col=c,
        )

    total_pairs = len(enriched_items) + len(not_enriched_items)
    pct_enriched = (100.0 * len(enriched_items) / total_pairs) if total_pairs else 0.0

    panel_legends[(r, c)] = dict(
        enriched=enriched_items,
        not_enriched=not_enriched_items,
        n_enriched=len(enriched_items),
        n_total=total_pairs,
        pct=pct_enriched,
    )

    fig.update_xaxes(range=[0, axis_max], tick0=0, dtick=dtick, row=r, col=c)
    fig.update_yaxes(range=[0, axis_max], tick0=0, dtick=dtick,
                     scaleanchor=f"x{i+1}", scaleratio=1, row=r, col=c)

    fig.update_yaxes(showticklabels=(c == 1), row=r, col=c)
    fig.update_xaxes(showticklabels=(r == nrows), row=r, col=c)


# -----------------------
# Panel legends (HTML annotations)
# -----------------------
layout = fig.layout

for i in range(n_panels):
    r = i // ncols + 1
    c = i % ncols + 1
    axis_index = i + 1

    yaxis_name = "yaxis" if axis_index == 1 else f"yaxis{axis_index}"
    ydom = getattr(layout, yaxis_name).domain

    info = panel_legends.get((r, c), {})
    enriched_items = info.get("enriched", [])
    not_items = info.get("not_enriched", [])
    n_enriched = info.get("n_enriched", 0)
    n_total = info.get("n_total", 0)
    pct = info.get("pct", 0.0)

    # left edge at x=3 (data)
    legend_x_paper = legend_xpaper_for_data_x(layout, axis_index, legend_data_x)

    # (1) legend entry for y=x should LOOK dashed.
    # We simulate a dashed sample using a monospace string of dashes with gaps.
    # dashed_sample = "<span style='font-family:monospace;'>- - -</span>"
    # dashed_sample = "<span style='font-family:monospace;'>- -</span>"
    dashed_sample = "<span style='font-family:monospace;'>— —</span>"

    legend_lines = [
        f"<span style='font-size:{panel_header_size}px;'><b>Pairs</b></span>",
        f"{dashed_sample} <span style='font-size:{panel_item_size}px; color:black;'>y = x</span>",
    ]

    # (your enriched / not enriched structure kept; if you still want the headers, keep them)
    if enriched_items:
        legend_lines.append(f"<span style='font-size:{panel_section_size}px;'><b>Enriched</b></span>")
        for label, colr in enriched_items:
            legend_lines.append(
                f"<span style='color:{colr};'>—</span> "
                f"<span style='font-size:{panel_item_size}px;'>{label}</span>"
            )

    if not_items or (not enriched_items):
        legend_lines.append(f"<span style='font-size:{panel_section_size}px;'><b>Not enriched</b></span>")
        if not_items:
            for label, colr in not_items:
                legend_lines.append(
                    f"<span style='color:{colr};'>—</span> "
                    f"<span style='font-size:{panel_item_size}px;'>{label}</span>"
                )
        else:
            legend_lines.append(
                f"<span style='font-size:{panel_item_size}px; color:rgba(0,0,0,0.55);'>None</span>"
            )

    legend_lines.append(
        f"<span style='font-size:{panel_footer_size}px; color:rgba(0,0,0,0.75);'>"
        f"({n_enriched}/{n_total} enriched; {pct:.0f}%)</span>"
    )

    fig.add_annotation(
        x=legend_x_paper,
        y=ydom[1] - legend_top_inset,
        xref="paper",
        yref="paper",
        xanchor="left",
        yanchor="top",
        text="<br>".join(legend_lines),
        showarrow=False,
        align="left",
        bgcolor="rgba(255,255,255,0.75)",
        bordercolor="rgba(0,0,0,0.15)",
        borderwidth=1,
        borderpad=4,
    )


# -----------------------
# Global axis titles (BIGGER)
# -----------------------
fig.add_annotation(
    x=0.5,
    y=global_x_title_y,
    xref="paper",
    yref="paper",
    text="Haplotype X mean editing [%]",
    showarrow=False,
    xanchor="center",
    yanchor="top",
    font=dict(size=global_axis_font_size),
)

fig.add_annotation(
    x=global_y_title_x,
    y=0.5,
    xref="paper",
    yref="paper",
    text="Haplotype Y mean editing [%]",
    showarrow=False,
    xanchor="right",
    yanchor="middle",
    textangle=-90,
    font=dict(size=global_axis_font_size),
)

# -----------------------
# Layout
# -----------------------
fig.update_layout(
    width=width,
    height=height,
    margin=dict(l=85, r=25, t=70, b=105),
)

# IMPORTANT: do NOT shrink all annotations globally anymore (this was making titles look small)
# fig.update_annotations(font=dict(size=12))

fig.show()


# %% [markdown]
# ### Global enrichment with chi-square permutations

# %%
diverse_haplotypes_editing_and_snps_reads_dfs_for_fisher_tests = [
    df
    for df in editing_and_snps_reads_dfs
    if validate_haplotypes_diversity(df, 2, 50)
    # if validate_haplotypes_diversity(df, 2, 300) or validate_haplotypes_diversity(df, 3, 100)
    # if validate_haplotypes_diversity(df, 3, 100)
]
ic(len(diverse_haplotypes_editing_and_snps_reads_dfs_for_fisher_tests));

# %%
with Pool(processes=6) as pool:
    fisher_haplotypes_dfs = pool.map(
        all_haplotypes_permutation_null,
        diverse_haplotypes_editing_and_snps_reads_dfs_for_fisher_tests
    )
concat_fisher_haplotypes_df = pd.concat(fisher_haplotypes_dfs, ignore_index=True)
concat_fisher_haplotypes_df

# %%
concat_fisher_haplotypes_df.sort_values("Chrom")

# %%
summarized_fisher_df = (
    concat_fisher_haplotypes_df
    .groupby("Chrom")["Reject H0 (FDR)"]
    .agg(
        NumHaplotypePairs=lambda x: x.size,
        NumSignificantPairs=lambda x: x.sum(),
        PrctSignificantPairs=lambda x: 100 * x.sum() / x.size
    )
    .reset_index()
    .rename(columns={"PrctSignificantPairs": "%SignificantPairs"})
)
summarized_fisher_df

# %%
summarized_fisher_df.loc[
    summarized_fisher_df["%SignificantPairs"].eq(0)
]

# %%
summarized_fisher_df["NumHaplotypePairs"].max()

# %%
fig = px.histogram(
    summarized_fisher_df,
    x="NumHaplotypePairs",
    # y="%SignificantPairs",
    labels={
        "NumHaplotypePairs": "Number of haplotype pairs",
        # "%SignificantPairs": "% of haplotype pairs with<br>significant editing differences",
        # "%SignificantPairs": "% of differentially edited haplotype pairs",
    },
    nbins=int(summarized_fisher_df["NumHaplotypePairs"].max()),
    histnorm="percent",
)
fig.update_xaxes(dtick=1)
fig.update_yaxes(title="% of tested genes")
fig.update_layout(
    width=600,
    height=450,
)
fig.show()

# %%
fig = px.histogram(
    summarized_fisher_df,
    x="%SignificantPairs",
    labels={
        "%SignificantPairs": "% of differentially edited haplotype pairs",
    },
    nbins=20,
    histnorm="percent",
)
fig.update_xaxes(dtick=10, range=[0, 100])
fig.update_yaxes(title="% of tested genes")
fig.update_layout(
    width=600,
    height=450,
)
fig.show()

# %%
fig = px.ecdf(
    summarized_fisher_df,
    x="%SignificantPairs",
    labels={
        "%SignificantPairs": "% of differentially edited haplotype pairs",
    },
    ecdfmode="reversed"
)
fig.update_xaxes(dtick=10, range=[0, 100])
fig.update_yaxes(title="Fraction of tested genes")
fig.update_layout(
    width=600,
    height=450,
    # ti
)
fig.show()

# %%
fig = px.scatter(
    summarized_fisher_df,
    x="NumHaplotypePairs",
    y="%SignificantPairs",
    labels={
        "NumHaplotypePairs": "Haplotype pairs",
        # "%SignificantPairs": "% of haplotype pairs with<br>significant editing differences",
        "%SignificantPairs": "% of differentially edited pairs",
    },
    trendline="ols"
)
fig.update_xaxes(dtick=1)
# fig.update_yaxes(title="% of tested genes")
fig.update_layout(
    width=600,
    height=450,
)
fig.show()

# %%
fig = px.scatter(
    summarized_fisher_df.groupby(
        "NumHaplotypePairs"
    )["%SignificantPairs"].agg(["mean", "std"]).reset_index(),
    x="NumHaplotypePairs",
    y="mean",
    error_y="std",
    labels={
        "NumHaplotypePairs": "Haplotype pairs",
        # "%SignificantPairs": "% of haplotype pairs with<br>significant editing differences",
        # "%SignificantPairs": "% of differentially edited pairs",
        "mean": "Mean % of differentially edited pairs",
    },
    # trendline="ols"
)
fig.update_xaxes(dtick=1)
# fig.update_yaxes(title="% of tested genes")
fig.update_layout(
    width=600,
    height=450,
)
fig.show()

# %%
fig = px.box(
    summarized_fisher_df,
    x="NumHaplotypePairs",
    y="%SignificantPairs",
    labels={
        "NumHaplotypePairs": "Haplotype pairs",
        "%SignificantPairs": "% of differentially edited pairs",
    },
    category_orders={
        "NumHaplotypePairs": summarized_fisher_df["NumHaplotypePairs"].sort_values().drop_duplicates().tolist()
    },
    points="all"
    # opacity=0.5,
    # markers=True, lines=False
)
fig.update_xaxes(dtick=1)
# fig.update_yaxes(title="Fraction of tested genes")
fig.update_layout(
    width=600,
    height=450,
    # ti
)
fig.show()

# %%
fig = px.histogram(
    summarized_fisher_df,
    x="%SignificantPairs",
    labels={
        "%SignificantPairs": "% of differentially edited haplotype pairs",
    },
    # nbins=101,
    histnorm="percent",
    facet_col="NumHaplotypePairs",
    category_orders={
        "NumHaplotypePairs": summarized_fisher_df["NumHaplotypePairs"].sort_values().drop_duplicates().tolist()
    }
)
fig.update_xaxes(dtick=20)
# fig.update_yaxes(title="% of tested genes")
fig.update_layout(
    # width=600,
    height=450,
)
fig.show()

# %%
fig = px.ecdf(
    summarized_fisher_df,
    x="%SignificantPairs",
    labels={
        "%SignificantPairs": "% of differentially edited haplotype pairs",
    },
    # ecdfmode="reversed",
    facet_col="NumHaplotypePairs",
    category_orders={
        "NumHaplotypePairs": summarized_fisher_df["NumHaplotypePairs"].sort_values().drop_duplicates().tolist()
    },
    # markers=True, lines=False
)
fig.update_xaxes(dtick=10, range=[0, 100])
# fig.update_yaxes(title="Fraction of tested genes")
fig.update_layout(
    # width=600,
    height=450,
    # ti
)
fig.show()

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

# %%
snps_and_coverage_per_gene_df = mismatches_df.loc[
    (mismatches_df["NoisyFinal"])
    & (mismatches_df["MismatchFrequency"].ge(snp_noise_level))
].groupby("Chrom").size().reset_index(name="SNPs").merge(
    alignment_stats_df.loc[:, ["Chrom", "MappedReads"]],
    how="outer"
).fillna(0)
snps_and_coverage_per_gene_df["SNPs"] = (
    snps_and_coverage_per_gene_df["SNPs"].astype(int)
)
snps_and_coverage_per_gene_df

# %%
num_genes_with_at_most_3_snps_and_at_least_50_reads = (
    snps_and_coverage_per_gene_df.loc[
        (snps_and_coverage_per_gene_df["SNPs"].le(max_snps_per_gene_to_allow_editing_detection))
        & (snps_and_coverage_per_gene_df["MappedReads"].ge(50))
    ].shape[0]
)
num_genes_with_at_most_3_snps_and_at_least_1000_reads = (
    snps_and_coverage_per_gene_df.loc[
        (snps_and_coverage_per_gene_df["SNPs"].le(max_snps_per_gene_to_allow_editing_detection))
        & (snps_and_coverage_per_gene_df["MappedReads"].ge(1000))
    ].shape[0]
)
num_genes_with_at_most_3_snps_and_at_least_50_reads, num_genes_with_at_most_3_snps_and_at_least_1000_reads

# %%
max_distinct_proteins_df = (
    distinct_unique_proteins_df.sort_values("Fraction", ascending=False)
    .groupby("Chrom")
    .apply(pd.DataFrame.nlargest, n=1, columns="NumOfProteins")
)
# max_distinct_proteins_df = (
#     max_distinct_proteins_df.drop("Chrom", axis=1).reset_index().drop("level_1", axis=1)
# )

# # max_distinct_proteins_df[condition_col] = max_distinct_proteins_df[
# #     condition_col
# # ].astype(str)

# max_distinct_proteins_df = max_distinct_proteins_df.merge(
#     tmr50_alignment_stats_df,
#     on="Chrom",
#     # how="left",
#     how="right",
# )

max_distinct_proteins_df = max_distinct_proteins_df.merge(
    snps_and_coverage_per_gene_df.loc[
        (snps_and_coverage_per_gene_df["SNPs"].le(max_snps_per_gene_to_allow_editing_detection))
        & (snps_and_coverage_per_gene_df["MappedReads"].ge(50))
    ],
    on="Chrom",
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
fig = px.scatter(
    max_distinct_proteins_df,
    x="NumOfReads",
    y="NumOfProteins"
    # log
)
fig.update_layout(
    width=500,
    height=500
)
fig.show()

# %%
max_distinct_proteins_df.loc[
    max_distinct_proteins_df["MappedReads"]
    < max_distinct_proteins_df["NumOfAvailableReads"]
]

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
max_distinct_proteins_df["NumOfProteins"].mean()

# %%
alignment_stats_df


# %%
def report_isoforms_per_gene(
    max_distinct_prots_df, 
    chroms_of_edited_positions=None,
    sig_digits=0, 
):    
    max_possibly_edited_genes = max_distinct_prots_df.shape[0]
    
    ic(max_possibly_edited_genes)
    
    if chroms_of_edited_positions is not None:
        num_genes_with_edited_positions = len(chroms_of_edited_positions)
        prct_of_genes_with_edited_positions = np.round(
            100 * num_genes_with_edited_positions / max_possibly_edited_genes,
            sig_digits
        )
        ic(num_genes_with_edited_positions, prct_of_genes_with_edited_positions)
    
    exactly_x_isoforms_per_gene = [1]
    
    at_least_x_isoforms_per_gene = [5, 50]
    
    for exactly_x_isoforms in exactly_x_isoforms_per_gene:
        genes_with_exactly_x_isoforms = max_distinct_prots_df.loc[
            max_distinct_prots_df["NumOfProteins"] == exactly_x_isoforms
        ].shape[0]
        prct_of_genes_with_exactly_x_isoforms = np.round(
            100 * genes_with_exactly_x_isoforms / max_possibly_edited_genes,
            sig_digits
        )
        ic(exactly_x_isoforms, genes_with_exactly_x_isoforms, prct_of_genes_with_exactly_x_isoforms)

    for at_least_x_isoforms in at_least_x_isoforms_per_gene:
        genes_with_at_least_x_isoforms = max_distinct_prots_df.loc[
            max_distinct_prots_df["NumOfProteins"] >= at_least_x_isoforms
        ].shape[0]
        prct_of_genes_with_at_least_x_isoforms = np.round(
            100 * genes_with_at_least_x_isoforms / max_possibly_edited_genes,
            sig_digits
        )
        ic(at_least_x_isoforms, genes_with_at_least_x_isoforms, prct_of_genes_with_at_least_x_isoforms)


# %%
report_isoforms_per_gene(max_distinct_proteins_df, chroms_of_edited_positions=chroms)

# %%

# %%
max_distinct_proteins_df["IsNeural"].value_counts()

# %%
max_distinct_proteins_df["IsNeural"].value_counts(normalize=True).mul(100).round(1)

# %%
fig = px.histogram(max_distinct_proteins_df, x="IsNeural", log_y=True)
fig.update_layout(width=600, height=400, template=template)
fig.show()

# %%
tmr1000_max_distinct_proteins_df = (
    tmr1000_distinct_unique_proteins_df.sort_values("Fraction", ascending=False)
    .groupby("Chrom")
    .apply(pd.DataFrame.nlargest, n=1, columns="NumOfProteins")
)
# tmr1000_max_distinct_proteins_df = (
#     tmr1000_max_distinct_proteins_df.drop("Chrom", axis=1)
#     .reset_index()
#     .drop("level_1", axis=1)
# )

# # max_distinct_proteins_df[condition_col] = max_distinct_proteins_df[
# #     condition_col
# # ].astype(str)

# tmr1000_max_distinct_proteins_df = tmr1000_max_distinct_proteins_df.merge(
#     tmr1000_alignment_stats_df,
#     on="Chrom",
#     # how="left",
#     how="right",
# )
tmr1000_max_distinct_proteins_df = tmr1000_max_distinct_proteins_df.merge(
    snps_and_coverage_per_gene_df.loc[
        (snps_and_coverage_per_gene_df["SNPs"].le(max_snps_per_gene_to_allow_editing_detection))
        & (snps_and_coverage_per_gene_df["MappedReads"].ge(1000))
    ],
    on="Chrom",
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

tmr1000_max_distinct_proteins_df = tmr1000_max_distinct_proteins_df.sort_values(
    "NumOfProteins", ascending=False, ignore_index=True
)

tmr1000_max_distinct_proteins_df

# %%
report_isoforms_per_gene(
    tmr1000_max_distinct_proteins_df, 
    # chroms_of_edited_positions=chroms
)

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

# fig.write_image(
#     "Distinct proteins per gene vs. % of genes - log(y) - Octopus.svg",
#     width=width,
#     height=height,
# )

fig.show()

# %%
# top panel - for main figures

# font_size = 24

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

y_min = 1

tmr50_legendtitle = "50 reads"
tmr1000_legendtitle = "1000 reads"
# legend_title_text = "Minimum coverage per gene      "

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
            # size=0.7 * font_size,
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
            # size=0.7 * font_size,
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
        # textfont=dict(size=0.7 * font_size),
        showlegend=False,
    ),
    row=1,
    col=1,
)


fig.update_xaxes(
    type="log",
    # tickfont=dict(size=0.7 * font_size),
)
fig.update_yaxes(
    type="log",
    # range=[-2, 2.2]
    range=[np.log(y_min) * 1.1 / np.log(10), 2.2],
    # tickfont=dict(size=0.7 * font_size),
)

# fig.update_annotations(font_size=font_size)

width = 800
# width = 900
height = 500

fig.update_layout(
    # xaxis_title="Distinct proteins per transcript",
    # yaxis_title="% of transcripts",
    # title="Pooled octopus data",
    # legend_font=dict(size=font_size),
    # title="Whole-transcriptome octopus data",
    title=dict(
        text="Whole-transcriptome octopus data", 
        # font=dict(size=font_size * 1.5)
    ),
    # title_x=0.15,
    title_y=0.93,
    template=template,
    width=width,
    height=height,
    # legend_title_text=legend_title_text,
    # legend_title_text="Min coverage per gene       ",
    # legend_title_text="Coverage per gene       ",
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
    "Distinct proteins per gene vs. % of genes - log(y) - Octopus - Pooled - top panel.svg",
    width=width,
    height=height,
)

fig.show()

# %%

# %%

# %%
# variables needed for main fig. top panel - 
# to be saved and combined with SC as a bottom panel

max_distinct_proteins_df.to_csv(
    Path(out_dir, "MaxDistinctProtsForFig6.TMR50.Octopus.Pooled.csv"),
    sep="\t",
    index=False
)
tmr1000_max_distinct_proteins_df.to_csv(
    Path(out_dir, "MaxDistinctProtsForFig6.TMR1000.Octopus.Pooled.csv"),
    sep="\t",
    index=False
)

# %%
# top panel - for main figures

# font_size = 24

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

y_min = 1

tmr50_legendtitle = "50 reads"
tmr1000_legendtitle = "1000 reads"
# legend_title_text = "Minimum coverage per gene      "

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
            # size=0.7 * font_size,
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
            # size=0.7 * font_size,
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
        # textfont=dict(size=0.7 * font_size),
        showlegend=False,
    ),
    row=1,
    col=1,
)


fig.update_xaxes(
    type="log",
    # tickfont=dict(size=0.7 * font_size),
)
fig.update_yaxes(
    type="log",
    # range=[-2, 2.2]
    range=[np.log(y_min) * 1.1 / np.log(10), 2.2],
    # tickfont=dict(size=0.7 * font_size),
)

# fig.update_annotations(font_size=font_size)

width = 800
# width = 900
height = 500

fig.update_layout(
    # xaxis_title="Distinct proteins per transcript",
    # yaxis_title="% of transcripts",
    # title="Pooled octopus data",
    # legend_font=dict(size=font_size),
    # title="Whole-transcriptome octopus data",
    title=dict(
        text="Whole-transcriptome octopus data", 
        # font=dict(size=font_size * 1.5)
    ),
    # title_x=0.15,
    title_y=0.93,
    template=template,
    width=width,
    height=height,
    # legend_title_text=legend_title_text,
    # legend_title_text="Min coverage per gene       ",
    # legend_title_text="Coverage per gene       ",
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
    "Distinct proteins per gene vs. % of genes - log(y) - Octopus - Pooled - top panel.svg",
    width=width,
    height=height,
)

fig.show()

# %%
# bottom panel - for supp. figures

# font_size = 24

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


# legend_title_text = "Minimum coverage per gene      "

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
        row=1,
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


# neural_vs_non_text = f"<b>Mann-Whitney U</b><br>p-val = {pv:.2e}<br>statistic = {statistic:.2g}",

if pv < 10**-22:
    # neural_vs_non_text = f"<b>Mann-Whitney U between<br>neural to non-neural cells</b><br>p-val < 1E-22<br>statistic = {statistic:.2g}"
    neural_vs_non_text = "<b>p-val < 1E-22</b><br>(Mann-Whitney U)"
else:
    # neural_vs_non_text = f"<b>Mann-Whitney U between<br>neural to non-neural cells</b><br>p-val = {pv:.2e}<br>statistic = {statistic:.2g}"
    neural_vs_non_text = f"<b>p-val = {pv:.2e}</b><br>(Mann-Whitney U)"

fig.add_annotation(
    x=np.log(10) / np.log(10),
    y=np.log(1) / np.log(10),
    xref="x",
    yref="y",
    # text=f"<b>Mann-Whitney U between<br>neural to non-neural genes</b><br>p-val = {pv:.2e}<br>statistic = {statistic:.2g}",
    text=neural_vs_non_text,
    bgcolor="white",
    borderpad=4,
    font=dict(size=13),
    opacity=0.8,
    showarrow=False,
    row=1,
    col=1,
)





fig.update_xaxes(
    type="log",
    # tickfont=dict(size=0.7 * font_size),
)
fig.update_yaxes(
    type="log",
    # range=[-2, 2.2]
    range=[np.log(y_min) * 1.1 / np.log(10), 2.2],
    # tickfont=dict(size=0.7 * font_size),
)

# fig.update_annotations(font_size=font_size)

width = 800
# width = 900
height = 500

fig.update_layout(
    # xaxis_title="Distinct proteins per transcript",
    # yaxis_title="% of transcripts",
    # title="Pooled octopus data",
    # legend_font=dict(size=font_size),
    # title="Whole-transcriptome octopus data",
    title=dict(
        text="Whole-transcriptome octopus data", 
        # font=dict(size=font_size * 1.5)
    ),
    # title_x=0.15,
    title_y=0.93,
    template=template,
    width=width,
    height=height,
    # legend_title_text=legend_title_text,
    # legend_title_text="Min coverage per gene       ",
    legend_title_text="Genes",
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
    Path(out_dir, "Distinct proteins per gene vs. % of genes - log(y) - Octopus - Pooled - bottom panel.svg"),
    width=width,
    height=height,
)

fig.show()

# %%
# version with only top panel for a conference

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

# width = 800
width = 900
height = 600

fig.update_layout(
    # xaxis_title="Distinct proteins per transcript",
    # yaxis_title="% of transcripts",
    # title="Pooled octopus data",
    legend_font=dict(size=font_size),
    # title="Whole-transcriptome octopus data",
    title=dict(
        text="Whole-transcriptome octopus data", font=dict(size=font_size * 1.5)
    ),
    # title_x=0.15,
    title_y=0.93,
    template=template,
    width=width,
    height=height,
    # legend_title_text=legend_title_text,
    # legend_title_text="Min coverage per gene       ",
    legend_title_text="Coverage per gene       ",
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
#     "Distinct proteins per gene vs. % of genes - log(y) - Octopus - top panel.svg",
#     width=width,
#     height=height,
# )

fig.show()

# %%

# %%
# version with only top panel for a conference

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

# width = 800
width = 900
height = 600

fig.update_layout(
    # xaxis_title="Distinct proteins per transcript",
    # yaxis_title="% of transcripts",
    # title="Pooled octopus data",
    legend_font=dict(size=font_size),
    # title="Whole-transcriptome octopus data",
    title=dict(
        text="Whole-transcriptome octopus data", font=dict(size=font_size * 1.5)
    ),
    # title_x=0.15,
    title_y=0.93,
    template=template,
    width=width,
    height=height,
    # legend_title_text=legend_title_text,
    # legend_title_text="Min coverage per gene       ",
    legend_title_text="Coverage per gene       ",
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
#     "Distinct proteins per gene vs. % of genes - log(y) - Octopus - top panel.svg",
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

expanded_max_distinct_proteins_df

# %%
(
    expanded_max_distinct_proteins_df.loc[
        expanded_max_distinct_proteins_df["Chrom"] == chrom
    ]
    .merge(
        expanded_unique_proteins_df,
        on=["Chrom", "Transcript", "Protein"],
        how="left",
    )
    .merge(samples_and_tissues_df, on="Sample", how="left")
)

# %%

strongly_diversified_expanded_max_distinct_proteins_dfs = [
    (
        expanded_max_distinct_proteins_df.loc[
            expanded_max_distinct_proteins_df["Chrom"] == chrom
        ]
        .merge(
            expanded_unique_proteins_df,
            on=["Chrom", "Transcript", "Protein"],
            how="left",
        )
        .merge(samples_and_tissues_df, on="Sample", how="left")
    )
    for chrom, expanded_unique_proteins_df in zip(chroms, expanded_unique_proteins_dfs)
    if chrom in strongly_diversified_chroms
]

del expanded_max_distinct_proteins_df

for (
    strongly_diversified_expanded_max_distinct_proteins_df
) in strongly_diversified_expanded_max_distinct_proteins_dfs:
    strongly_diversified_expanded_max_distinct_proteins_df.rename(
        columns={"Sample": "Sample2", "Protein": "Protein2", "Tissue": "Tissue2"},
        inplace=True,
    )
    strongly_diversified_expanded_max_distinct_proteins_df.insert(
        2, "Sample", strongly_diversified_expanded_max_distinct_proteins_df["Sample2"]
    )
    strongly_diversified_expanded_max_distinct_proteins_df.insert(
        3, "Tissue", strongly_diversified_expanded_max_distinct_proteins_df["Tissue2"]
    )
    strongly_diversified_expanded_max_distinct_proteins_df.insert(
        4, "Protein", strongly_diversified_expanded_max_distinct_proteins_df["Protein2"]
    )
    strongly_diversified_expanded_max_distinct_proteins_df.drop(
        ["Sample2", "Tissue2", "Protein2"], axis=1, inplace=True
    )

# expanded_max_distinct_proteins_dfs[1]
strongly_diversified_expanded_max_distinct_proteins_dfs[0].loc[
    :, :"MaxNonSynsFrequency"
]

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
    
    num_of_unique_proteins = strongly_diversified_expanded_max_distinct_proteins_df["Protein"].unique().size
    
    df2 = (
        gb.apply(lambda x: 100 * len(x) / num_of_unique_proteins)
        .reset_index()
        .rename(columns={0: "%RelativeNumOfProteins"})
    )
    df = df.merge(df2)

    strongly_diversified_num_of_proteins_per_sample_dfs.append(df)

strongly_diversified_num_of_proteins_per_sample_dfs[0]

# %% jupyter={"source_hidden": true}
# total_non_unique_proteins_per_strongly_diversified_chroms = [
#     len(strongly_diversified_expanded_max_distinct_proteins_df)
#     for strongly_diversified_expanded_max_distinct_proteins_df in strongly_diversified_expanded_max_distinct_proteins_dfs
# ]
# total_non_unique_proteins_per_strongly_diversified_chroms

# %% jupyter={"source_hidden": true}
# total_non_unique_proteins_per_strongly_diversified_chroms_2 = [
#     strongly_diversified_num_of_proteins_per_sample_df["NumOfProteins"].sum()
#     for strongly_diversified_num_of_proteins_per_sample_df in strongly_diversified_num_of_proteins_per_sample_dfs
# ]
# total_non_unique_proteins_per_strongly_diversified_chroms_2

# %% jupyter={"source_hidden": true}
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
concat_strongly_diversified_per_transcript_per_sample_coverage_df = pd.concat(strongly_diversified_per_transcript_per_sample_coverage_dfs)
concat_strongly_diversified_per_transcript_per_sample_coverage_df

# %%
strongly_diversified_num_of_reads_and_proteins_per_sample_dfs = []
for strongly_diversified_num_of_proteins_per_sample_df in strongly_diversified_num_of_proteins_per_sample_dfs:
    strongly_diversified_num_of_reads_and_proteins_per_sample_df = strongly_diversified_num_of_proteins_per_sample_df.merge(
        concat_strongly_diversified_per_transcript_per_sample_coverage_df,
        on=["Chrom", "Transcript", "Sample", "Tissue"],
        how="left"
    )
    assert strongly_diversified_num_of_reads_and_proteins_per_sample_df["NumOfReads"].ge(strongly_diversified_num_of_reads_and_proteins_per_sample_df["NumOfProteins"]).all()
    strongly_diversified_num_of_reads_and_proteins_per_sample_df["NumOfReads/NumOfProteins"] = strongly_diversified_num_of_reads_and_proteins_per_sample_df["NumOfReads"] / strongly_diversified_num_of_reads_and_proteins_per_sample_df["NumOfProteins"]
    strongly_diversified_num_of_reads_and_proteins_per_sample_dfs.append(strongly_diversified_num_of_reads_and_proteins_per_sample_df)
concat_strongly_diversified_num_of_reads_and_proteins_per_sample_df = pd.concat(strongly_diversified_num_of_reads_and_proteins_per_sample_dfs)
concat_strongly_diversified_num_of_reads_and_proteins_per_sample_df

# %%
fig = px.box(
    concat_strongly_diversified_num_of_reads_and_proteins_per_sample_df,
    x="Tissue",
    y="NumOfProteins",
    color="Tissue",
    color_discrete_map=tissues_color_discrete_map,
    points="all",
    category_orders={"Tissue": tissues_order}
)
fig.update_layout(
    width=900,
    height=500,
    template=template,
    showlegend=False
)
fig.show()

# %%
fig = px.box(
    concat_strongly_diversified_num_of_reads_and_proteins_per_sample_df,
    x="Tissue",
    y="NumOfReads/NumOfProteins",
    color="Tissue",
    color_discrete_map=tissues_color_discrete_map,
    points="all",
    category_orders={"Tissue": tissues_order}
)
fig.update_yaxes(
    range=[0, concat_strongly_diversified_num_of_reads_and_proteins_per_sample_df["NumOfReads/NumOfProteins"].max()*1.4], 
    tick0=0, dtick=1,    
    zeroline=True, zerolinewidth=2
)
fig.update_layout(
    width=900,
    height=500,
    template=template,
    showlegend=False
)

fig.show()

# %%
concat_strongly_diversified_num_of_reads_and_proteins_per_sample_df["NumOfReads"].max()

# %%
concat_strongly_diversified_num_of_reads_and_proteins_per_sample_df["NumOfProteins"].max()

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
    strongly_diversified_num_of_reads_and_proteins_per_sample_df = concat_strongly_diversified_num_of_reads_and_proteins_per_sample_df.loc[
        concat_strongly_diversified_num_of_reads_and_proteins_per_sample_df["Chrom"] == strongly_diversified_chrom
    ]
    
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
            strongly_diversified_num_of_reads_and_proteins_per_sample_df["Tissue"] == tissue,
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
# ### Distinct isoforms per sample 2

# %%
expression_dfs[0]

# %%
f1_concat_expression_df = pd.concat(
    [
        expression_df.loc[
            expression_df["Fraction"] == 1,
            [
                "Chrom",
                condition_col,
                "Protein",
                "#Solution",
                "Fraction",
                "FractionRepetition",
                "Algorithm",
                "AlgorithmRepetition",
                "Samples",
                # "Reads",
                # "AdditionalSupportingReadsIDs",
                "TotalWeightedSupportingReads",
            ],
        ]
        for expression_df in expression_dfs
    ]
)

f1_concat_expression_df

# %%
max_distinct_proteins_df

# %%
max_f1_concat_expression_df = f1_concat_expression_df.merge(
    max_distinct_proteins_df.loc[
        :,
        [
            "Chrom",
            condition_col,
            #  'NumOfReads',
            "Fraction",
            "FractionRepetition",
            "Algorithm",
            "AlgorithmRepetition",
            "NumOfProteins",
            #    'Proteins',
            #    'AvailableReads', 'NumOfAvailableReads', 'MappedReads',
            # "Samples",
            # "MappedReadsPerSample",
            #    'KnownSites', 'IsNeural'
        ],
    ],
    # on=[
    #     "Chrom",
    #     condition_col,
    #     "Fraction",
    #     "FractionRepetition",
    #     "Algorithm",
    #     "AlgorithmRepetition",
    # ],
).drop(
    columns=[
        "Fraction",
        "FractionRepetition",
        "Algorithm",
        "AlgorithmRepetition",
        "#Solution",
    ]
)

max_f1_concat_expression_df.insert(1, "MockChrom", max_f1_concat_expression_df["Chrom"])

max_f1_concat_expression_df


# %%
def per_sample_contribution_stats(df, digits=3):
    chrom = df["MockChrom"].iloc[0]
    condition = df[condition_col].iloc[0]
    num_of_proteins = df["NumOfProteins"].iloc[0]

    # test_chrom, test_condition
    supported_isoforms_per_sample = (
        (df["Samples"].str.split(",").apply(set).explode().value_counts())
        .reset_index(name="NumOfProteinsPerGenePerSample")
        .rename(columns={"Samples": "Sample"})
    )

    supported_isoforms_per_sample.insert(0, "Chrom", chrom)
    supported_isoforms_per_sample.insert(1, condition_col, condition)
    supported_isoforms_per_sample.insert(3, "NumOfProteinsPerGene", num_of_proteins)

    supported_isoforms_per_sample["%OfProteinsPerGenePerSample"] = (
        supported_isoforms_per_sample["NumOfProteinsPerGenePerSample"]
        .div(num_of_proteins)
        .mul(100)
        .round(digits)
    )

    return supported_isoforms_per_sample


# %%
# chrom = "comp106951_c0_seq1"
# df = max_f1_concat_expression_df.loc[max_f1_concat_expression_df["Chrom"] == chrom]
# df

# %%
# per_sample_contribution_stats(df)

# %%
# concat_per_transcript_per_sample_coverage_dfs.loc[:, ["Chrom", "Sample", "NumOfReads"]]

# %%
per_sample_contribution_df = (
    max_f1_concat_expression_df.groupby("Chrom")
    .apply(per_sample_contribution_stats, include_groups=False)
    .reset_index(drop=True)
)

temp_samples_and_tissues_df = samples_and_tissues_df.set_index("Sample")
per_sample_contribution_df.insert(
    per_sample_contribution_df.columns.get_loc("Sample") + 1,
    "Tissue",
    per_sample_contribution_df.apply(
        lambda x: temp_samples_and_tissues_df.at[x["Sample"], "Tissue"], axis=1
    ),
)
del temp_samples_and_tissues_df

per_sample_contribution_df = per_sample_contribution_df.merge(
    concat_per_transcript_per_sample_coverage_df.loc[
        :, ["Chrom", "Sample", "NumOfReads"]
    ],
    how="left",
)

per_sample_contribution_df

# %%

# %%
# 9 chroms mostly diversified due to A-to-I RNA editing
strongly_diversified_max_num_of_proteins = (
    max_distinct_proteins_df.sort_values("NumOfProteins", ascending=False)
    .iloc[:9]["NumOfProteins"]
    .to_list()
)

strongly_diversified_chroms = (
    max_distinct_proteins_df.sort_values("NumOfProteins", ascending=False)
    .iloc[:9]["Chrom"]
    .to_list()
)

strongly_diversified_transcripts = (
    max_distinct_proteins_df.sort_values("NumOfProteins", ascending=False)
    .iloc[:9][condition_col]
    .to_list()
)

strongly_diversified_chroms

# %%
concat_strongly_diversified_num_of_reads_and_proteins_per_sample_df = (
    per_sample_contribution_df.loc[
        per_sample_contribution_df["Chrom"].isin(strongly_diversified_chroms)
    ]
)

concat_strongly_diversified_num_of_reads_and_proteins_per_sample_df

# %%
cols = min(facet_col_wrap, len(strongly_diversified_chroms), 3)
rows = ceil(len(strongly_diversified_chroms) / cols)
row_col_iter = list(product(range(1, rows + 1), range(1, cols + 1)))[
    : len(strongly_diversified_chroms)
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
    # vertical_spacing=0.12,
    # horizontal_spacing=0.12,
)

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
            "NumOfProteinsPerGenePerSample",
        ]
        # ic(samp)

        # max_x = max(max_x, x.max())
        # max_y = max(max_y, y.max())

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

max_x = concat_strongly_diversified_num_of_reads_and_proteins_per_sample_df[
    "NumOfReads"
].max()
max_y = concat_strongly_diversified_num_of_reads_and_proteins_per_sample_df[
    "NumOfProteinsPerGenePerSample"
].max()
max_x_y = ceil(max(max_x, max_y)) * 1.05
# max_x_y = max(max_x, max_y) * 1.05

fig.update_xaxes(
    # showticklabels=False,  # Hide x axis ticks
    # categoryorder="array",
    # categoryarray=tissues_order,
    # range=[0, max_x * 1.1],
    range=[0, max_x_y],
    tick0=0,
    dtick=100,
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
    dtick=100,
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
    Path(out_dir, "Distinct proteins vs coverage per 9 strong transcripts - Octopus.svg"),
    width=width,
    height=height,
)

fig.show()

# %%
per_sample_contribution_df

# %%
fig = px.histogram(
    per_sample_contribution_df,
    x="%OfProteinsPerGenePerSample",
    color="Tissue",
    color_discrete_map=tissues_color_discrete_map,
    # opacity=0.7,
    facet_row="Tissue",
)
fig.update_layout(
    template=template,
    height=1000,
    width=600,
    # barmode="overlay",
    showlegend=False,
)
fig.show()

# %%
ceil(max(max_x, max_y)) * 1.05

# %%
per_sample_contribution_df["%OfProteinsPerGenePerSample"].describe()

# %%
per_sample_contribution_df.groupby("Tissue")["%OfProteinsPerGenePerSample"].describe()

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
    ["Chrom", condition_col, "Algorithm"]
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
        subset=["Chrom", condition_col, "Algorithm"], ignore_index=True
    )
)

max_distinct_proteins_per_transcript_and_alg_df

# %%
# mean distinct proteins per transcript
max_distinct_proteins_per_transcript_and_alg_df.sort_values(
    ["Chrom", condition_col, "NumOfProteins"], ascending=False
).drop_duplicates(["Chrom", condition_col], ignore_index=True)["NumOfProteins"].mean()

# %%
# num of transcripts with at least 5 variants
max_distinct_proteins_per_transcript_and_alg_df.sort_values(
    ["Chrom", condition_col, "NumOfProteins"], ascending=False
).drop_duplicates(["Chrom", condition_col], ignore_index=True)["NumOfProteins"].ge(5).sum()

# %%
# num of transcripts with at least 50 variants
max_distinct_proteins_per_transcript_and_alg_df.sort_values(
   ["Chrom", condition_col, "NumOfProteins"], ascending=False
).drop_duplicates(["Chrom", condition_col], ignore_index=True)["NumOfProteins"].ge(50).sum()


# %%
def align_algorithm_results(df, keys):
    unexpected = set(df["Algorithm"].dropna()) - {"Ascending", "Descending"}
    if unexpected:
        raise ValueError(f"Unexpected algorithms: {sorted(unexpected)}")
    asc = df.loc[df["Algorithm"].eq("Ascending")].set_index(keys)
    desc = df.loc[df["Algorithm"].eq("Descending")].set_index(keys)
    assert asc.index.is_unique
    assert desc.index.is_unique
    if len(asc.index.difference(desc.index)) or len(desc.index.difference(asc.index)):
        raise ValueError("Ascending/Descending have different gene/subsample keys")
    asc = asc.sort_index()
    desc = desc.reindex(asc.index)
    return asc.reset_index(), desc.reset_index()

# asc_df = max_distinct_proteins_per_transcript_and_alg_df.loc[
#     max_distinct_proteins_per_transcript_and_alg_df["Algorithm"] == "Ascending"
# ].reset_index(drop=True)
# desc_df = max_distinct_proteins_per_transcript_and_alg_df.loc[
#     max_distinct_proteins_per_transcript_and_alg_df["Algorithm"] != "Ascending"
# ].reset_index(drop=True)

asc_df, desc_df = align_algorithm_results(
    max_distinct_proteins_per_transcript_and_alg_df,
    ["Chrom"],
)

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

greater_asc_chroms = asc_df.loc[
    asc_df["NumOfProteins"].gt(desc_df["NumOfProteins"]), "Chrom"
].unique()
greater_asc_chroms

# %%
distinct_unique_proteins_df.loc[
    (distinct_unique_proteins_df["Fraction"] == 1.0)
    # & (distinct_unique_proteins_df[condition_col].isin(greater_asc_transcripts))
    & (distinct_unique_proteins_df["Chrom"].isin(greater_asc_chroms))
].groupby(["Chrom", condition_col, "Algorithm"])["NumOfProteins"].value_counts()

# %%
max_distinct_proteins_per_transcript_and_alg_df.loc[
    max_distinct_proteins_per_transcript_and_alg_df["Chrom"].isin(
        greater_asc_chroms
    )
]

# %%
max_distinct_proteins_per_transcript_and_alg_and_fracrepetition_df = (
    distinct_unique_proteins_df.copy()
)

max_distinct_proteins_per_transcript_and_alg_and_fracrepetition_df[
    "MaxNumOfProteins"
] = max_distinct_proteins_per_transcript_and_alg_and_fracrepetition_df.groupby(
    ["Chrom", condition_col, "Fraction", "Algorithm", "FractionRepetition"]
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
        subset=["Chrom", condition_col, "Fraction", "Algorithm", "FractionRepetition"],
        ignore_index=True,
    )
)
max_distinct_proteins_per_transcript_and_alg_and_fracrepetition_df = (
    max_distinct_proteins_per_transcript_and_alg_and_fracrepetition_df.sort_values(
        ["Chrom", condition_col, "Fraction", "FractionRepetition", "Algorithm"],
        ignore_index=True,
    )
)

max_distinct_proteins_per_transcript_and_alg_and_fracrepetition_df

# %%
# asc_df = max_distinct_proteins_per_transcript_and_alg_and_fracrepetition_df.loc[
#     max_distinct_proteins_per_transcript_and_alg_and_fracrepetition_df["Algorithm"]
#     == "Ascending"
# ].reset_index(drop=True)
# desc_df = max_distinct_proteins_per_transcript_and_alg_and_fracrepetition_df.loc[
#     max_distinct_proteins_per_transcript_and_alg_and_fracrepetition_df["Algorithm"]
#     != "Ascending"
# ].reset_index(drop=True)

asc_df, desc_df = align_algorithm_results(
    max_distinct_proteins_per_transcript_and_alg_and_fracrepetition_df,
    ["Chrom", "Fraction", "FractionRepetition"],
)

assert len(asc_df) == len(desc_df)

ic(asc_df["NumOfProteins"].eq(desc_df["NumOfProteins"]).sum())  # ==
ic(asc_df["NumOfProteins"].gt(desc_df["NumOfProteins"]).sum())  # >
ic(asc_df["NumOfProteins"].lt(desc_df["NumOfProteins"]).sum())
# <

# %%
# greater_asc_transcripts = asc_df.loc[
#     asc_df["NumOfProteins"].gt(desc_df["NumOfProteins"]), condition_col
# ].unique()

# ic(len(greater_asc_transcripts))

# greater_asc_transcripts

greater_asc_chroms = asc_df.loc[
    asc_df["NumOfProteins"].gt(desc_df["NumOfProteins"]), "Chrom"
].unique()
greater_asc_chroms

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
np.round(100 * len(dispersion_df.loc[dispersion_df["%SolutionsDispersion"] > 0]) / len(dispersion_df), 1)

# %%
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
# ##### Saving dispersion df

# %%
saved_dispersion_df = dispersion_df.rename(columns={"Transcript": "Gene"})
saved_dispersion_df.insert(0, "Platform", "Whole-transcriptome octopus data")
saved_dispersion_df.to_csv(
    Path(out_dir, "Dispersion.Octopus.WholeTranscriptome.Pooled.tsv"), 
    sep="\t", index=False
)
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

# %% [markdown]
# ### Distinct isoforms vs. coverage, compared to targeted squid genes

# %%
homologs_dir = Path("/private6/projects/Combinatorics/SquidVsOctopusHomologs")


# %%
subprocess.run(
    f"python /private7/projects/Combinatorics/Code/squid_vs_oct_orfs.py --output-dir {homologs_dir}", 
    shell=True
)

# %%
squid_to_oct_homologs_file = Path(homologs_dir, "comprehensive_homologs_clean.csv")
squid_to_oct_homologs_df = pd.read_csv(squid_to_oct_homologs_file)
squid_to_oct_homologs_df = squid_to_oct_homologs_df.loc[
    squid_to_oct_homologs_df["MatchDirection"] != "No match"
]
squid_to_oct_homologs_df = squid_to_oct_homologs_df.sort_values("IsReciprocal", ascending=False, ignore_index=True)
squid_to_oct_homologs_df = squid_to_oct_homologs_df.loc[
    squid_to_oct_homologs_df["OctopusToSquidCoverage"] >= 0.3
]
squid_to_oct_homologs_df = squid_to_oct_homologs_df.merge(
    tmr50_alignment_stats_df.loc[:, ["Chrom", "MappedReads"]].rename(columns={"Chrom": "OctopusGene", "MappedReads": "OctopusMappedReads"}),
    how="left"
)
squid_to_oct_homologs_df = squid_to_oct_homologs_df.sort_values("TargetGene", ignore_index=True)
squid_to_oct_homologs_df = squid_to_oct_homologs_df.loc[
    squid_to_oct_homologs_df["OctopusMappedReads"] > 0
].reset_index(drop=True)
squid_to_oct_homologs_df = squid_to_oct_homologs_df.drop(columns=["IsTargetGene", "MatchedTarget", "SquidUniProt", "MatchDirection"])
squid_to_oct_homologs_df

# %%
homologs_of_tageted_squid_genes_distinct_unique_proteins_df = distinct_unique_proteins_df.merge(
    # squid_to_oct_homologs_df.loc[:, ["OctopusGene", "SquidGene", "SquidUniProt"]],
    squid_to_oct_homologs_df.loc[:, ["OctopusGene", "SquidGene", "TargetGene"]].rename(columns={"TargetGene": "SquidShortUniProt"}),
    how="inner",
    left_on="Chrom",
    right_on="OctopusGene"
)
homologs_of_tageted_squid_genes_distinct_unique_proteins_df

# %%
hsd = Path("/private6/projects/Combinatorics/D.pealeii/OctopusHomologs/MpileupAndTranscripts")

# %%
homologous_squid_genes = [
    "K0513",
    "KCNAS",
    "PCLO",
    "SCN1",
    "GRIA2",
    "PCLO",
    "ADAR1",
]
homologous_squid_chroms = ['comp141840_c0_seq2',
 'comp141640_c0_seq1',
 'comp141882_c0_seq14',
 'comp141378_c0_seq7',
 'comp141693_c0_seq1',
 'comp141882_c0_seq14',
 'comp134400_c0_seq1_extended']
# homologous_squid_chroms
homologous_squid_platforms = ["Short-reads"] * 4 + ["Long-reads"] * 3


# %%
def dashed_platform_rep_to_upper_camel_case(platform):
    a, b = platform.split("-")
    b = b.title()
    return a + b


# %%
homologous_squid_unique_reads_files = [
    list(hsd.glob(f"*{gene}*{dashed_platform_rep_to_upper_camel_case(platform)}*.unique_reads.csv.gz"))[0]
    for gene, platform in zip(
        homologous_squid_genes, homologous_squid_platforms
    )
]
ic(homologous_squid_unique_reads_files);

homologous_squid_distinct_proteins_files = [
    list(hsd.glob(f"*{gene}*{dashed_platform_rep_to_upper_camel_case(platform)}*.DistinctUniqueProteins*.csv"))[0]
    for gene, platform in zip(
        homologous_squid_genes, homologous_squid_platforms
    )
]
ic(homologous_squid_distinct_proteins_files);

# %%
homologous_squid_distinct_unique_proteins_dfs = []

for gene, chrom, platform, distinct_unique_proteins_file, unique_reads_file in zip(
    homologous_squid_genes, homologous_squid_chroms, homologous_squid_platforms, homologous_squid_distinct_proteins_files, homologous_squid_unique_reads_files
):
    
    unique_reads_df = pd.read_csv(unique_reads_file, sep=sep, dtype={"UniqueRead": str, "Reads": str})
    unique_reads_df.insert(0, "Chrom", chrom)
    
    _distinct_unique_proteins_df = pd.read_csv(distinct_unique_proteins_file, sep=sep)
    _distinct_unique_proteins_df.insert(0, "Gene", gene)
    _distinct_unique_proteins_df.insert(
        1,
        "NumOfReads",
        (
            _distinct_unique_proteins_df["Fraction"]
            * unique_reads_df["NumOfReads"].sum()
        ).astype(int),
    )
    _distinct_unique_proteins_df.insert(0, "Chrom", chrom)
    _distinct_unique_proteins_df.insert(1, "Platform", platform)
    _distinct_unique_proteins_df.insert(
        _distinct_unique_proteins_df.columns.get_loc("AvailableReads") + 1,
        "NumOfAvailableReads",
        _distinct_unique_proteins_df["AvailableReads"].str.split(",").str.len(),
    )
    homologous_squid_distinct_unique_proteins_dfs.append(_distinct_unique_proteins_df)

homologous_squid_distinct_unique_proteins_df = (
    pd.concat(homologous_squid_distinct_unique_proteins_dfs)
    .reset_index(drop=True)
    .rename(columns={"NumUniqueSamples": "NumOfProteins", "UniqueSamples": "Proteins"})
)

homologous_squid_distinct_unique_proteins_df = homologous_squid_distinct_unique_proteins_df.sort_values(
    [
        "Gene",
        "Fraction",
        "FractionRepetition",
        "Algorithm",
        "AlgorithmRepetition",
    ]
).reset_index(drop=True)

homologous_squid_distinct_unique_proteins_df

# %%
max_per_fraction_homologs_of_tageted_squid_genes_distinct_unique_proteins_df = homologs_of_tageted_squid_genes_distinct_unique_proteins_df.sort_values(
    ["SquidShortUniProt", "Fraction", "NumOfProteins"],
    ascending=False
).drop_duplicates(
    ["SquidShortUniProt", "Fraction"], 
    ignore_index=True
)
max_per_fraction_homologs_of_tageted_squid_genes_distinct_unique_proteins_df

# %%
max_per_fraction_homologous_squid_distinct_unique_proteins_df = homologous_squid_distinct_unique_proteins_df.sort_values(
    ["Gene", "Platform", "Fraction", "NumOfProteins"],
    ascending=False
).drop_duplicates(
    ["Gene", "Platform", "Fraction"], 
    ignore_index=True
)
max_per_fraction_homologous_squid_distinct_unique_proteins_df

# %%
max_per_fraction_homologs_of_tageted_squid_genes_distinct_unique_proteins_df.groupby(
    "SquidShortUniProt"
)[["NumOfAvailableReads", "NumOfProteins"]].max().reset_index().sort_values("SquidShortUniProt")

# %%
max_per_fraction_homologous_squid_distinct_unique_proteins_df.groupby(
    ["Gene", "Platform"]
)[["NumOfAvailableReads", "NumOfProteins"]].max().reset_index().sort_values(["Gene", "Platform"])


# %%
def add_trace_to_squid_vs_oct(fig, row, col, legend_constructed, x, y, platform, color_discrete_map):
    if not legend_constructed:
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=y,
                    marker_color=color_discrete_map[platform],
                    name=platform,
                    # legendrank=tissue_to_legendrank[tissue],
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
                marker_color=color_discrete_map[platform],
                # name=tissue,
                showlegend=False,
                # marker_pattern_shape="/",
                mode="markers",
            ),
            row=row,
            col=col,
        )

    legend_constructed = True
    
    return legend_constructed


cols = min(facet_col_wrap, len(set(homologous_squid_genes)), 3)
rows = ceil(len(set(homologous_squid_genes)) / cols)
row_col_iter = list(product(range(1, rows + 1), range(1, cols + 1)))[
    : len(set(homologous_squid_genes))
]

gene_row_col_dict = {(row, col): gene for (row, col), gene in zip(row_col_iter, set(homologous_squid_genes))}
ic(gene_row_col_dict)

x_title = "Coverage"
y_title = "Distinct protein isoforms"


subplot_titles = [
    gene_row_col_dict[(row, col)]
    for (row, col) in row_col_iter
]
ic(subplot_titles);

homologs_color_seq = px.colors.qualitative.G10[:3]

homologs_color_discrete_map = {
    platform: color
    for platform, color in zip(
        ["Octopus long-reads", "Squid long-reads", "Squid short-reads"],
        # homologs_color_seq
        ["black"] + homologs_color_seq[:2]
    )
}

# homologs_color_discrete_map

fig = make_subplots(
    rows=rows,
    cols=cols,
    subplot_titles=subplot_titles,
    shared_xaxes="all",
    shared_yaxes="all",
    x_title=x_title,
    y_title=y_title,
    vertical_spacing=0.1,
    # horizontal_spacing=0.12,
)

homologs_legend_constructed = {
    platform: False
    for platform in homologs_color_discrete_map
}

for row, col in row_col_iter:
    
    squid_gene = gene_row_col_dict[(row, col)]
    
    
    octopus_df = max_per_fraction_homologs_of_tageted_squid_genes_distinct_unique_proteins_df.loc[
        max_per_fraction_homologs_of_tageted_squid_genes_distinct_unique_proteins_df["SquidShortUniProt"]
        == squid_gene
    ].copy()
    x = octopus_df["NumOfAvailableReads"]
    y = octopus_df["NumOfProteins"]
    platform = "Octopus long-reads"
    legend_constructed = homologs_legend_constructed[platform]
    homologs_legend_constructed[platform] = add_trace_to_squid_vs_oct(fig, row, col, legend_constructed, x, y, platform, homologs_color_discrete_map)
    
    squid_df = max_per_fraction_homologous_squid_distinct_unique_proteins_df.loc[
        max_per_fraction_homologous_squid_distinct_unique_proteins_df["Gene"]
        == squid_gene
    ].copy()
    for squid_platform in squid_df["Platform"].unique():
        _squid_df = squid_df.loc[
            squid_df["Platform"] == squid_platform
        ]
        x = _squid_df["NumOfAvailableReads"]
        y = _squid_df["NumOfProteins"]
        platform = f"Squid {squid_platform.lower()}"
        legend_constructed = homologs_legend_constructed[platform]
        homologs_legend_constructed[platform] = add_trace_to_squid_vs_oct(fig, row, col, legend_constructed, x, y, platform, homologs_color_discrete_map)

# fig.update_xaxes(tickangle=45, automargin=True)

# max_x = max(
#     max_per_fraction_homologs_of_tageted_squid_genes_distinct_unique_proteins_df["NumOfAvailableReads"].max(),
#     max_per_fraction_homologous_squid_distinct_unique_proteins_df["NumOfAvailableReads"].max(),
# )
# max_y = max(
#     max_per_fraction_homologs_of_tageted_squid_genes_distinct_unique_proteins_df["NumOfProteins"].max(),
#     max_per_fraction_homologous_squid_distinct_unique_proteins_df["NumOfProteins"].max(),
# )
# max_x_y = ceil(max(max_x, max_y)) * 1.05
# # max_x_y = max(max_x, max_y) * 1.05

fig.update_traces(opacity=0.7, marker_size=8)

width = 1000
height = 600

fig.update_layout(
    template=template,
    width=width,
    height=height,
)

fig.show()

# now also plot the fig as log-log

fig.update_xaxes(
    type="log"
)
fig.update_yaxes(
    type="log"
)
# fig.write_image(
#     "Distinct proteins vs coverage per 9 strong transcripts - Octopus.svg",
#     width=width,
#     height=height,
# )

fig.show()

# %%
raise Error()

# %%

# %%
homologous_color_sequence = px.colors.qualitative.G10
homologous_conditions = target_squid_gene_names
homologous_color_discrete_map = {
    condition: color for condition, color in zip(homologous_conditions, homologous_color_sequence)
}

# %%
fig = px.scatter(
    homologs_of_tageted_squid_genes_distinct_unique_proteins_df.sort_values(
        ["SquidUniProt", "Fraction", "NumOfProteins"],
        ascending=False
    ).drop_duplicates(
        ["SquidUniProt", "Fraction"], 
        ignore_index=True
    ),
    x="NumOfReads",
    y="NumOfProteins",
    facet_col="SquidUniProt",
    facet_col_wrap=3,
    color="SquidUniProt",
    log_x=True,
    log_y=True,
)
fig.update_layout(
    width=800,
    height=600,
    showlegend=False,
    template=template,
)
fig.show()

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
for condition in homologous_conditions:
    ic(condition)
    df = homologs_of_tageted_squid_genes_distinct_unique_proteins_df.loc[
        homologs_of_tageted_squid_genes_distinct_unique_proteins_df["SquidUniProt"] == condition
    ]
    df = df.sort_values(["Fraction", "NumOfProteins"], ascending=False).drop_duplicates(
        "Fraction", ignore_index=True
    )

    color = homologous_color_discrete_map[condition]
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
                legendgroup="Squid's full-CDS,<br>long-reads",  # this can be any string
                # legendgrouptitle_text="Squid's full-CDS,<br>long-reads",
                legendgrouptitle_text="Octopus homologs of<br>targeted squid genes",
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
                legendgroup="Squid's full-CDS,<br>long-reads",  # this can be any string
                name=name,
            ),
        )

    maximal_x = max(maximal_x, x_measured.max())

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
#         legendgroup="Drosophila’s Dscam",  # this can be any string
#         legendgrouptitle_text="Drosophila’s Dscam",
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
#         legendgroup="Drosophila’s Dscam",  # this can be any string
#         name=dscam_legend_names[1],
#         # name=f"DSCAM {dscam_legend_names[0]}",
#     ),
# )

fig.update_yaxes(
    # type="log",
    # # tick0=0,
    # dtick="D2",
    # # exponentformat="power",
    # showexponent='all',
    # range=[0, (floor(np.log10(max_y)) + ceil(np.log10(max_y))) / 2],
    # range=[0, max(max_y * 1.2, 1.05 * max(dscam_ys))],
    range=[0, None],
    # type='log'
    # zeroline=True
)
fig.update_xaxes(range=[0, maximal_x * 1.1])
fig.update_layout(
    # title_text=head_title,
    title_text="Distinct isoforms in pooled octopus samples of targeted squid homologs",
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
# fig.write_image(
#     "Distinct proteins vs. sequencing depth - PacBio.svg", width=900, height=600
# )
fig.show()

# %%

# %%

# %%

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
