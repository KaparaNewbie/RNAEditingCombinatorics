# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.18.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Imports

# %%
code_dir = "/private7/projects/Combinatorics/Code"

out_dir = "/private7/projects/Combinatorics/Code/Notebooks"

# %% papermill={"duration": 2.901153, "end_time": "2022-02-01T09:42:46.125355", "exception": false, "start_time": "2022-02-01T09:42:43.224202", "status": "completed"}
import sys
from functools import reduce
from itertools import chain, combinations, product
from math import ceil, floor
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
from pybedtools import BedTool
from sklearn import linear_model
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.metrics import mean_squared_error, r2_score

# from numpy.random import RandomState


sys.path.append(str(Path(code_dir).absolute()))
from Alignment.alignment_utils import count_reads, count_reads_in_unaligned_bam
from EditingUtils.logo import make_freq_df, multiple_logos_from_fasta_files
from EditingUtils.seq import make_fasta_dict

# %%
px.colors

# %%
# %load_ext autoreload
# %autoreload 2

# %% jupyter={"source_hidden": true} papermill={"duration": 0.071769, "end_time": "2022-02-01T09:42:43.049672", "exception": false, "start_time": "2022-02-01T09:42:42.977903", "status": "completed"} tags=["parameters"]
# condition_col = "Gene"
# conditions = ["GRIA", "PCLO"]
# # conditions = ["GRIA2_HUMAN", "PCLO_CHICK"]
# chroms = ["comp141693_c0_seq1", "comp141882_c0_seq14"]
# starts = [170, 0]
# ends = [2999, 6294]
# strands = ["+", "+"]
# unaligned_bam_files = [
#     "/private7/projects/Combinatorics/D.pealeii/Data/CCS/BasicCCS/GRIA-CNS-RESUB.C0x1291.ccs.bam",
#     "/private7/projects/Combinatorics/D.pealeii/Data/CCS/BasicCCS/PCLO-CNS-RESUB.C0x1291.ccs.bam",
# ]
# reads_type = "CCS"  # something like CCS / miseq / etc.
# aligned_bam_files = [
#     "/private7/projects/Combinatorics/D.pealeii/Alignment/BestN1/GRIA-CNS-RESUB.C0x1291.aligned.sorted.bam",
#     "/private7/projects/Combinatorics/D.pealeii/Alignment/BestN1/PCLO-CNS-RESUB.C0x1291.aligned.sorted.bam",
# ]
# filtered_aligned_bam_files = [
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.bam",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.bam",
# ]
# include_flags = None
# exclude_flags = "2304"  # remove secondary and supplementary (chimeric) alignments
# sep = "\t"
# positions_files = [
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.positions.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.positions.csv",
# ]
# reads_files = [
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.reads.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.reads.csv",
# ]
# unique_reads_files = [
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_reads.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_reads.csv",
# ]
# proteins_files = [
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.proteins.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.proteins.csv",
# ]
# unique_proteins_files = [
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv",
# ]
# reads_first_col_pos = 6
# unique_reads_first_col_pos = 8
# proteins_first_col_pos = 12
# unique_proteins_first_col_pos = 14
# reads_editing_col = "EditingFrequency"
# proteins_editing_col = "MinNonSyns"
# distinct_unique_reads_files = [
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA.AllRows.DistinctUniqueReads.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO.AllRows.DistinctUniqueReads.csv",
# ]
# distinct_unique_proteins_files = [
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA.AllRows.DistinctUniqueProteins.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO.AllRows.DistinctUniqueProteins.csv",
# ]
# proteins_jaccard_files = [
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA.AllRows.JaccardMatrixProteins.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO.AllRows.JaccardMatrixProteins.csv",
# ]
# expression_files = [
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA.DistinctUniqueProteins.ExpressionLevels.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO.DistinctUniqueProteins.ExpressionLevels.csv",
# ]
# distinct_dissimilar_miyata_proteins_files = [
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.DistinctUniqueProteins.AAgroupsMiyata1979.06.12.2022-20:29:59.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.DistinctUniqueProteins.AAgroupsMiyata1979.06.12.2022-20:48:08.csv"
# ]
# miyata_expression_files = [
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA.DistinctUniqueProteins.ExpressionLevels.AAgroupsMiyata1979.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO.DistinctUniqueProteins.ExpressionLevels.AAgroupsMiyata1979.csv"
# ]
# grantham_cutoff_scores = [
#     50, 75, 100,
#     # 125, 150
# ]
# distinct_dissimilar_grantham_proteins_files = [
#     [
#         "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.DistinctUniqueProteins.GRANTHAM1974-50.06.12.2022-20:24:15.csv",
#         "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.DistinctUniqueProteins.GRANTHAM1974-75.06.12.2022-21:09:37.csv",
#         "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.DistinctUniqueProteins.GRANTHAM1974-100.07.12.2022-08:25:55.csv",
#         # "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.DistinctUniqueProteins.GRANTHAM1974-125.06.12.2022-22:44:40.csv",
#         # "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.DistinctUniqueProteins.GRANTHAM1974-150.07.12.2022-00:03:50.csv"
#     ],
#     [
#         "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.DistinctUniqueProteins.GRANTHAM1974-50.06.12.2022-20:38:37.csv",
#         "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.DistinctUniqueProteins.GRANTHAM1974-75.06.12.2022-22:01:57.csv",
#         "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.DistinctUniqueProteins.GRANTHAM1974-100.07.12.2022-09:37:48.csv",
#         # "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.DistinctUniqueProteins.GRANTHAM1974-125.07.12.2022-00:06:52.csv",
#         # "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.DistinctUniqueProteins.GRANTHAM1974-150.07.12.2022-02:17:58.csv"
#     ]
# ]
# grantham_expression_files = [
#     [
#         "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA.DistinctUniqueProteins.ExpressionLevels.GRANTHAM1974-50.csv",
#         "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA.DistinctUniqueProteins.ExpressionLevels.GRANTHAM1974-75.csv",
#         "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA.DistinctUniqueProteins.ExpressionLevels.GRANTHAM1974-100.csv",
#         # "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA.DistinctUniqueProteins.ExpressionLevels.GRANTHAM1974-125.csv",
#         # "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA.DistinctUniqueProteins.ExpressionLevels.GRANTHAM1974-150.csv"
#     ],
#     [
#         "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO.DistinctUniqueProteins.ExpressionLevels.GRANTHAM1974-50.csv",
#         "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO.DistinctUniqueProteins.ExpressionLevels.GRANTHAM1974-75.csv",
#         "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO.DistinctUniqueProteins.ExpressionLevels.GRANTHAM1974-100.csv",
#         # "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO.DistinctUniqueProteins.ExpressionLevels.GRANTHAM1974-125.csv",
#         # "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO.DistinctUniqueProteins.ExpressionLevels.GRANTHAM1974-150.csv"
#     ]
# ]
# alg_repetitions = 5
# known_sites_file = (
#     "/private7/projects/Combinatorics/D.pealeii/Annotations/D.pea.EditingSites.csv"
# )
# samtools_path = "/home/alu/kobish/anaconda3/envs/combinatorics/bin/samtools"
# threads = 20
# code_dir = "/private7/projects/Combinatorics/Code"
# seed = 1892

# %% jupyter={"source_hidden": true} papermill={"duration": 0.071769, "end_time": "2022-02-01T09:42:43.049672", "exception": false, "start_time": "2022-02-01T09:42:42.977903", "status": "completed"} tags=["parameters"]
# condition_col = "Gene"
# conditions = [
#     "RUSC2_MOUSE",
#     "TRIM2_BOVIN",
#     "CA2D3_MOUSE",
#     "ABL_DROME",
#     "DGLA_HUMAN",
#     "K0513_MOUSE",
#     "KCNAS_DROME",
#     "ACHA4_MOUSE",
#     "ANR17_HUMAN",
#     "TWK7_CAEEL",
#     "SCN1_HETBL",
#     "CACB2_RABIT",
#     "RIMS2_RAT",
#     "PCLO_CHICK",
#     "DOP1_HUMAN",
#     "IQEC1_HUMAN",
#     "CSKI1_MOUSE",
# ]
# chroms = [
#     "comp141881_c0_seq3",
#     "comp141044_c0_seq2",
#     "comp140439_c0_seq1",
#     "comp126362_c0_seq1",
#     "comp141517_c0_seq1",
#     "comp141840_c0_seq2",
#     "comp141640_c0_seq1",
#     "comp140987_c3_seq1",
#     "comp140910_c2_seq1",
#     "comp136058_c0_seq1",
#     "comp141378_c0_seq7",
#     "comp141158_c1_seq2",
#     "comp140712_c0_seq3",
#     "comp141882_c0_seq14",
#     "comp141880_c1_seq3",
#     "comp141565_c6_seq3",
#     "comp141684_c0_seq1",
# ]
# starts = [192, 352, 1, 2, 171, 400, 1, 385, 2, 2, 443, 288, 311, 1, 256, 989, 71]
# ends = [
#     6296,
#     2766,
#     3495,
#     2377,
#     3359,
#     3732,
#     1617,
#     2190,
#     3460,
#     1246,
#     7381,
#     1808,
#     2116,
#     6294,
#     7722,
#     4195,
#     2857,
# ]
# strands = [
#     "+",
#     "+",
#     "+",
#     "+",
#     "+",
#     "+",
#     "+",
#     "+",
#     "+",
#     "+",
#     "+",
#     "+",
#     "+",
#     "+",
#     "+",
#     "+",
#     "+",
# ]
# # unaligned_bam_files = [
# #     "/private7/projects/Combinatorics/D.pealeii/Data/CCS/BasicCCS/GRIA-CNS-RESUB.C0x1291.ccs.bam",
# #     "/private7/projects/Combinatorics/D.pealeii/Data/CCS/BasicCCS/PCLO-CNS-RESUB.C0x1291.ccs.bam",
# # ]
# reads_type = "miseq"  # something like CCS / miseq / etc.
# aligned_bam_files = [
#     f"/private7/projects/Combinatorics/D.pealeii/Alignment/Illumina/reads.ByChrom/reads.sorted.aligned.filtered.{chrom}.bam"
#     for chrom in chroms
# ]
# # filtered_aligned_bam_files = [
# #     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.bam",
# #     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.bam",
# # ]
# include_flags = None
# exclude_flags = "2304"  # remove secondary and supplementary (chimeric) alignments
# sep = "\t"
# positions_files = [
#     f"/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.{chrom}.positions.csv"
#     for chrom in chroms
# ]
# reads_files = [
#     f"/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.{chrom}.reads.csv"
#     for chrom in chroms
# ]
# unique_reads_files = [
#     f"/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.{chrom}.unique_reads.csv"
#     for chrom in chroms
# ]
# proteins_files = [
#     f"/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.{chrom}.proteins.csv"
#     for chrom in chroms
# ]
# unique_proteins_files = [
#     f"/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.{chrom}.unique_proteins.csv"
#     for chrom in chroms
# ]
# reads_first_col_pos = 6
# unique_reads_first_col_pos = 8
# proteins_first_col_pos = 12
# unique_proteins_first_col_pos = 14
# reads_editing_col = "EditingFrequency"
# proteins_editing_col = "MinNonSyns"

# # todo update & uncomment
# # distinct_unique_reads_files = [
# #     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA.AllRows.DistinctUniqueReads.csv",
# #     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO.AllRows.DistinctUniqueReads.csv",
# # ]

# distinct_unique_proteins_files = [
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141881_c0_seq3.DistinctUniqueProteins.12.07.2022-20:54:38.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141044_c0_seq2.DistinctUniqueProteins.13.07.2022-06:33:23.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp140439_c0_seq1.DistinctUniqueProteins.12.07.2022-22:51:22.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp126362_c0_seq1.DistinctUniqueProteins.15.07.2022-06:11:18.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141517_c0_seq1.DistinctUniqueProteins.14.07.2022-07:43:15.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141840_c0_seq2.DistinctUniqueProteins.13.07.2022-20:30:25.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141640_c0_seq1.DistinctUniqueProteins.12.07.2022-19:44:02.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp140987_c3_seq1.DistinctUniqueProteins.18.07.2022-07:50:43.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp140910_c2_seq1.DistinctUniqueProteins.13.07.2022-16:15:35.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp136058_c0_seq1.DistinctUniqueProteins.21.07.2022-07:57:53.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141378_c0_seq7.DistinctUniqueProteins.19.07.2022-08:12:24.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141158_c1_seq2.DistinctUniqueProteins.13.07.2022-01:54:59.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp140712_c0_seq3.DistinctUniqueProteins.05.11.2022-02_11_33.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141882_c0_seq14.DistinctUniqueProteins.11.11.2022-22_23_54.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141880_c1_seq3.DistinctUniqueProteins.14.11.2022-20_58_20.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141565_c6_seq3.DistinctUniqueProteins.16.11.2022-12_42_59.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141684_c0_seq1.DistinctUniqueProteins.19.11.2022-01_51_21.csv",
# ]

# # proteins_jaccard_files = [
# #     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141881_c0_seq3.JaccardMatrixProteins.12.07.2022-20:54:38.csv",
# #     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141044_c0_seq2.JaccardMatrixProteins.13.07.2022-06:33:23.csv",
# #     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp140439_c0_seq1.JaccardMatrixProteins.12.07.2022-22:51:22.csv",
# #     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp126362_c0_seq1.JaccardMatrixProteins.15.07.2022-06:11:18.csv",
# #     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141517_c0_seq1.JaccardMatrixProteins.14.07.2022-07:43:15.csv",
# #     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141840_c0_seq2.JaccardMatrixProteins.13.07.2022-20:30:25.csv",
# #     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141640_c0_seq1.JaccardMatrixProteins.12.07.2022-19:44:02.csv",
# #     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp140987_c3_seq1.JaccardMatrixProteins.18.07.2022-07:50:43.csv",
# #     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp140910_c2_seq1.JaccardMatrixProteins.13.07.2022-16:15:35.csv",
# #     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp136058_c0_seq1.JaccardMatrixProteins.21.07.2022-07:57:53.csv",
# #     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141378_c0_seq7.JaccardMatrixProteins.19.07.2022-08:12:24.csv",
# #     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141158_c1_seq2.JaccardMatrixProteins.13.07.2022-01:54:59.csv",
# # ]

# expression_files = [
#     f"/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/{condition}.DistinctUniqueProteins.ExpressionLevels.csv"
#     for condition in conditions
# ]

# # todo update & uncomment
# # alg_repetitions = 5

# known_sites_file = (
#     "/private7/projects/Combinatorics/D.pealeii/Annotations/D.pea.EditingSites.csv"
# )
# samtools_path = "/home/alu/kobish/anaconda3/envs/combinatorics/bin/samtools"
# threads = 20
# code_dir = "/private7/projects/Combinatorics/Code"
# seed = 1892

# %% [markdown] papermill={"duration": 0.040192, "end_time": "2022-02-01T09:42:46.214429", "exception": false, "start_time": "2022-02-01T09:42:46.174237", "status": "completed"}
# # Ploting utils

# %%
facet_col_spacing = 0.05
template = "plotly_white"
facet_col_wrap = 6
facet_row_spacing = facet_col_spacing * 6
zerolinewidth = 4

# %% jupyter={"source_hidden": true} papermill={"duration": 0.054755, "end_time": "2022-02-01T09:42:46.304499", "exception": false, "start_time": "2022-02-01T09:42:46.249744", "status": "completed"}
# # plotly consts
# # color_sequence = px.colors.qualitative.Pastel
# # color_sequence = px.colors.qualitative.D3
# color_sequence = px.colors.qualitative.G10
# color_discrete_map = {
#     condition: color for condition, color in zip(conditions, color_sequence)
# }
# subcolors_discrete_map = {
#     condition: two_subcolors_from_hex(color_discrete_map[condition])
#     for condition in conditions
# }
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


# %% [markdown]
#
# plotly.colors.n_colors(lowcolor, highcolor, n_colors, colortype='tuple')
#
#     Splits a low and high color into a list of n_colors colors in it
#
#     Accepts two color tuples and returns a list of n_colors colors which form the intermediate colors between lowcolor and highcolor from linearly interpolating through RGB space. If colortype is ‘rgb’ the function will return a list of colors in the same form.
#

# %% jupyter={"source_hidden": true}
# def n_repetitions_colormap(subcolors_discrete_map, condition, n_repetitions):
#     lowcolor, highcolor = subcolors_discrete_map[condition]
#     colors = pc.n_colors(lowcolor, highcolor, n_repetitions, colortype="rgb")
#     return {i: color for i, color in enumerate(colors, start=1)}

# %% jupyter={"source_hidden": true}
# n_repetitions_colormap(subcolors_discrete_map, "GRIA", 10)

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


# %% [markdown] papermill={"duration": 0.041741, "end_time": "2022-02-01T09:42:47.760215", "exception": false, "start_time": "2022-02-01T09:42:47.718474", "status": "completed"}
# # Editing levels in PCLO - squid

# %%
px.colors

# %%

# %%
known_sites_file = (
    "/private7/projects/Combinatorics/D.pealeii/Annotations/D.pea.EditingSites.csv"
)
positions_files = [
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.positions.csv.gz",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.comp141882_c0_seq14.positions.csv",
]
sep = "\t"
strands = ["+", "+"]
condition_col = "Platform"
conditions = ["PacBio", "Illumina"]
chrom = "comp141882_c0_seq14"
colors = [
    px.colors.qualitative.G10[1],
    px.colors.qualitative.Dark24[13],
    # px.colors.qualitative.Set3[11]
    px.colors.qualitative.Vivid[5],
]
start = 0
end = 6294

# %%
known_sites_df = pd.read_csv(known_sites_file)

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

# known_sites_df["Position"] = known_sites_df["Position"].astype(int)
known_sites_df = known_sites_df.loc[
    (known_sites_df["Chrom"] == chrom)
    & (start <= known_sites_df["Position"])
    & (known_sites_df["Position"] < end)
]

known_sites_df


# %%
positions_dfs = []
for positions_file, condition in zip(positions_files, conditions):
    positions_df = pd.read_csv(positions_file, sep=sep).drop(
        ["CDS", "Phred"], axis=1, errors="ignore"
    )  # ignore error from droping these cols which only exist in Illumina data
    positions_df.insert(0, condition_col, condition)
    positions_dfs.append(positions_df)
positions_dfs[0]

# %%
ref_base_positions_dfs = []

for positions_df, strand in zip(positions_dfs, strands):
    ref_base = "A" if strand == "+" else "T"
    alt_base = "G" if strand == "+" else "C"
    positions_df = (
        positions_df.loc[positions_df["RefBase"] == ref_base]
        .sort_values(["Chrom", "Position"])
        .reset_index(drop=True)
    )
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

ref_base_positions_dfs[0]


# %%
col_suffixes = conditions + ["Known"]
col_suffixes

# %%
editing_freq_cols = [f"EditingFrequency_{col_suffix}" for col_suffix in col_suffixes]
editing_percent_cols = [f"%Editing_{col_suffix}" for col_suffix in col_suffixes]
edited_cols = [f"Edited_{col_suffix}" for col_suffix in col_suffixes]


# %%
def current_study_min_max_positions(merged_ref_base_positions_df):
    currently_edited_positions = merged_ref_base_positions_df.loc[
        (merged_ref_base_positions_df["Edited_PacBio"]) | (merged_ref_base_positions_df["Edited_Illumina"]),
    ].sort_values("Position")["Position"]
    min_position = currently_edited_positions.iloc[0]
    max_position = currently_edited_positions.iloc[-1]
    return min_position, max_position


# %%
merged_ref_base_positions_df = (
    ref_base_positions_dfs[0]
    .loc[:, ["Chrom", "Position", "EditingFrequency"]]
    .merge(
        ref_base_positions_dfs[1].loc[:, ["Chrom", "Position", "EditingFrequency"]],
        how="outer",
        on=["Chrom", "Position"],
        suffixes=[f"_{condition}" for condition in conditions],
    )
)

merged_ref_base_positions_df = merged_ref_base_positions_df.merge(
    known_sites_df.loc[
        known_sites_df["Chrom"] == chrom, ["Chrom", "Position", "EditingFrequency"]
    ].rename(columns={"EditingFrequency": "EditingFrequency_Known"}),
    how="outer",
    on=["Chrom", "Position"],
)

for editing_freq_col, editing_percent_col, edited_col in zip(
    editing_freq_cols, editing_percent_cols, edited_cols
):
    merged_ref_base_positions_df[editing_percent_col] = (
        merged_ref_base_positions_df[editing_freq_col].fillna(0) * 100
    )
    merged_ref_base_positions_df[edited_col] = (
        merged_ref_base_positions_df[editing_percent_col] > 0
    )

merged_ref_base_positions_df = merged_ref_base_positions_df.drop(
    editing_freq_cols, axis=1
)

# merged_ref_base_positions_df["Position"] = merged_ref_base_positions_df["Position"].astype("str")
merged_ref_base_positions_df = merged_ref_base_positions_df.sort_values(
    ["Chrom", "Position"]
).reset_index(drop=True)

current_study_min_x, current_study_max_x = current_study_min_max_positions(merged_ref_base_positions_df)
merged_ref_base_positions_df = merged_ref_base_positions_df.loc[
    (merged_ref_base_positions_df["Position"] >= current_study_min_x) &
    (merged_ref_base_positions_df["Position"] <= current_study_max_x)
]
merged_ref_base_positions_df = merged_ref_base_positions_df.loc[
    (merged_ref_base_positions_df["Edited_PacBio"]) |
    (merged_ref_base_positions_df["Edited_Illumina"]) |
    (merged_ref_base_positions_df["Edited_Known"])
]

merged_ref_base_positions_df

# %%
col_suffixes_updated_legend = {
    "PacBio": "Long-reads",
    "Illumina": "Short-reads",
    "Known": "Known (Short-reads)    ",
}

# %%
fig = make_subplots(
    rows=2,
    cols=1,
    shared_xaxes=True,
    # subplot_titles=["Known positions", "New positions"],
    row_titles=["Known<br>positions", "New<br>positions"],
    x_title="Position",
    y_title="% editing in squid's PCLO",
    # vertical_spacing=0.15,
    vertical_spacing=0.13,
)

symbols = ["circle", "square"]

for row in [1, 2]:
    if row == 1:
        df = merged_ref_base_positions_df.loc[
            merged_ref_base_positions_df["Edited_Known"]
        ]
    else:
        df = merged_ref_base_positions_df.loc[
            ~merged_ref_base_positions_df["Edited_Known"]
        ]

    for color, editing_percent_col, edited_col, col_suffix in zip(
        px.colors.qualitative.Pastel[:3],
        editing_percent_cols,
        edited_cols,
        col_suffixes,
    ):
        if row == 2 and col_suffix == "Known":
            continue

        x = df["Position"]
        y = df[editing_percent_col]

        if row == 1:
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=y,
                    mode="lines+markers",
                    marker=dict(
                        color=color,
                        opacity=0.7,
                        size=4,
                        # line=dict(
                        #     width=2,
                        #     color='DarkSlateGrey'
                        # )
                    ),
                    # name=col_suffix,
                    name=col_suffixes_updated_legend[col_suffix],
                ),
                row=row,
                col=1,
            )
        else:
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=y,
                    mode="lines+markers",
                    marker=dict(color=color, opacity=0.7, size=4),
                    showlegend=False,
                ),
                row=row,
                col=1,
            )

# fig.update_xaxes(
#     range=[current_study_min_x, current_study_max_x],
# )

# custom_tick_labels = ['0', '1', '2', '3', '4', '5']
# fig.update_xaxes(ticktext=custom_tick_labels)
fig.update_yaxes(
    range=[0, 100],
    tick0=0,
    dtick=20
    # range=[0, 100], tick0=0, dtick=20,
    # tick0=0, dtick=20
    # range=[0, 2],
    # range=[-2.8, 2],
    # type="log",
    # nticks=3, tick0=0,
    # ticktext=['0', '1', '2', '3', '4', '5']
)

width = 1100
height = 600

fig.update_layout(
    template=template,
    width=width,
    height=height,
    legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=0.9),
)

# fig.write_image(
#     "Comparison of editing levels across platforms and studies in squid’s PCLO.svg",
#     width=width,
#     height=height,
# )

fig.show()

# %%
fig = make_subplots(
    rows=2,
    cols=1,
    shared_xaxes=True,
    # subplot_titles=["Known positions", "New positions"],
    row_titles=["Known<br>positions", "De-novo<br>positions"],
    x_title="Position",
    y_title="Editing level in squid's PCLO [%]",
    # vertical_spacing=0.15,
    vertical_spacing=0.13,
)

symbols = ["circle", "square"]

# current_study_min_x = end
# current_study_max_x = start

# current_study_min_x, current_study_max_x = current_study_min_max_positions(merged_ref_base_positions_df)

# colors = [
#     px.colors.qualitative.G10[1],
#     # px.colors.qualitative.Dark24[13],
#     # px.colors.qualitative.Set3[11]
#     px.colors.qualitative.G10[4],
#     px.colors.qualitative.Vivid[5],
# ]
# colors = px.colors.qualitative.Pastel[:3]
colors = px.colors.qualitative.Pastel[:2] + [px.colors.qualitative.Pastel[4]]

for row in [1, 2]:
    if row == 1:
        df = merged_ref_base_positions_df.loc[
            merged_ref_base_positions_df["Edited_Known"]
        ]
    else:
        df = merged_ref_base_positions_df.loc[
            ~merged_ref_base_positions_df["Edited_Known"]
        ]

    for color, editing_percent_col, edited_col, col_suffix in zip(
        colors,
        editing_percent_cols,
        edited_cols,
        col_suffixes,
    ):
        if row == 2 and col_suffix == "Known":
            continue

        x = df["Position"]
        y = df[editing_percent_col]
        
        if row == 1:
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=y,
                    # mode="lines+markers",
                    mode="markers",
                    marker=dict(
                        color=color,
                        opacity=0.7,
                        # size=4,
                        size=6,
                        # line=dict(
                        #     width=2,
                        #     color='DarkSlateGrey'
                        # )
                    ),
                    # name=col_suffix,
                    name=col_suffixes_updated_legend[col_suffix],
                ),
                row=row,
                col=1,
            )
        else:
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=y,
                    # mode="lines+markers",
                    mode="markers",
                    marker=dict(color=color, 
                                opacity=0.7, 
                                # size=4,
                                size=6
                               ),
                    showlegend=False,
                ),
                row=row,
                col=1,
            )

# fig.update_xaxes(
#     range=[current_study_min_x, current_study_max_x],
# )

lowest_y_greater_than_0 = (
    pd.Series(
        merged_ref_base_positions_df.loc[
            :, ["%Editing_PacBio", "%Editing_Illumina", "%Editing_Known"]
        ].values.reshape(-1)
    )
    .replace(0, np.nan)
    .min()
)

fig["layout"]["yaxis"].update(
    range=[0, 100],
    tick0=0,
)

fig["layout"]["yaxis2"].update(
    range=[np.log(lowest_y_greater_than_0) / np.log(10), 2],
    type="log",
    tick0=0,
)

width = 1100
height = 600

fig.update_layout(
    template=template,
    width=width,
    height=height,
    legend=dict(
        orientation="h", 
        x=0.83,
        y=0.8,
        xref="container",
        yref="container",
        xanchor="right",
    ),
)

fig.write_image(
    "Comparison of editing levels across platforms and studies in squid’s PCLO - log(y) for new positions.svg",
    width=width,
    height=height,
)

fig.show()

# %%
cat_merged_ref_base_positions_df = merged_ref_base_positions_df.copy()
cat_merged_ref_base_positions_df["Position"] = cat_merged_ref_base_positions_df["Position"].astype(str)
cat_merged_ref_base_positions_df

# %%
# cat_merged_ref_base_positions_df["Position"].astype(int).sort_values().astype(str).tolist()

# %%
fig = make_subplots(
    rows=2,
    cols=1,
    shared_xaxes="all",
    # subplot_titles=["Known positions", "New positions"],
    row_titles=["Known<br>positions", "De-novo<br>positions"],
    x_title="Position",
    y_title="Editing level in squid's PCLO [%]",
    # vertical_spacing=0.15,
    vertical_spacing=0.13,
)

# current_study_min_x = end
# current_study_max_x = start



# colors = [
#     px.colors.qualitative.G10[1],
#     # px.colors.qualitative.Dark24[13],
#     # px.colors.qualitative.Set3[11]
#     px.colors.qualitative.G10[4],
#     px.colors.qualitative.Vivid[5],
# ]
# colors = px.colors.qualitative.Pastel[:3]
colors = px.colors.qualitative.Pastel[:2] + [px.colors.qualitative.Pastel[4]]

for row in [1, 2]:
    if row == 1:
        df = cat_merged_ref_base_positions_df.loc[
            cat_merged_ref_base_positions_df["Edited_Known"]
        ]
    else:
        df = cat_merged_ref_base_positions_df.loc[
            ~cat_merged_ref_base_positions_df["Edited_Known"]
        ]

    for color, editing_percent_col, edited_col, col_suffix in zip(
        colors,
        editing_percent_cols,
        edited_cols,
        col_suffixes,
    ):
        if row == 2 and col_suffix == "Known":
            continue

        x = df["Position"]
        y = df[editing_percent_col]

        if row == 1:
            fig.add_trace(
                go.Bar(
                    x=x,
                    y=y,
                    marker=dict(
                        color=color,
                        opacity=0.7,
                        # line=dict(
                        #     width=2,
                        #     color='DarkSlateGrey'
                        # )
                    ),
                    # name=col_suffix,
                    name=col_suffixes_updated_legend[col_suffix],
                ),
                row=row,
                col=1,
            )
        else:
            fig.add_trace(
                go.Bar(
                    x=x,
                    y=y,
                    marker=dict(color=color, 
                                opacity=0.7,
                               ),
                    showlegend=False,
                ),
                row=row,
                col=1,
            )

# ic(current_study_min_x, current_study_max_x)

fig.update_xaxes(
    type='category',
    title_standoff = 70,
    automargin=False
    # range=[current_study_min_x, current_study_max_x],
)

lowest_y_greater_than_0 = (
    pd.Series(
        cat_merged_ref_base_positions_df.loc[
            :, ["%Editing_PacBio", "%Editing_Illumina", "%Editing_Known"]
        ].values.reshape(-1)
    )
    .replace(0, np.nan)
    .min()
)

fig["layout"]["yaxis"].update(
    range=[0, 100],
    tick0=0,
    dtick=20
)

fig["layout"]["yaxis2"].update(
    range=[np.log(lowest_y_greater_than_0) / np.log(10), 2],
    type="log",
    tick0=0,
)

width = 1400
height = width * 500 / 1100

fig.update_annotations(
    font_size=14,
    # yshift=-50, 
)

fig.update_layout(
    template=template,
    width=width,
    height=height,
    legend=dict(
        orientation="h", 
        x=0.85,
        y=0.8,
        xref="container",
        yref="container",
        xanchor="right",
    ),
    barmode='overlay', 
    xaxis={
        'categoryorder':'array', 
        'categoryarray':cat_merged_ref_base_positions_df["Position"].astype(int).sort_values().astype(str).tolist()
    },
    bargap=0
    
)

# fig.write_image(
#     "Comparison of editing levels across platforms and studies in squid’s PCLO - log(y) for new positions.svg",
#     width=width,
#     height=height,
# )

fig.show()

# %%
editing_percent_cols, edited_cols, col_suffixes,

# %%
current_study_min_x, current_study_max_x

# %% jupyter={"source_hidden": true}
# fig = go.Figure()

# x = merged_ref_base_positions_df["Position"]

# # for color, editing_percent_col, col_suffix in zip(colors, editing_percent_cols, col_suffixes):
# for color, editing_percent_col, col_suffix in zip(
#     px.colors.qualitative.Pastel[:3], editing_percent_cols, col_suffixes
# ):

#     y = merged_ref_base_positions_df[editing_percent_col]

#     fig.add_trace(
#         go.Scatter(
#             x=x,
#             y=y,
#             mode="lines+markers",
#             marker=dict(color=color, opacity=0.5, size=6),
#             name=col_suffix,
#         )
#     )

# fig.update_layout(template=template, height=600)

# fig.show()

# %% [markdown]
# # ADAR motif - squid

# %%
condition_col = "Platform"
conditions = ["PacBio", "Illumina"]
sep = "\t"
# notice a transcriptome is different from a genome when considering up- and downstream bases of edited adenosines in order to plot ADAR's motif
# due to splicing, so it may be better to take the bases from the reads
transcriptome_file = (
    "/private7/projects/Combinatorics/D.pealeii/Annotations/orfs_squ.fa"
)

# %%
pacbio_data_df = pd.DataFrame(
    {
        # condition_col: ["GRIA", "PCLO"],
        "Chrom": ["comp141693_c0_seq1", "comp141882_c0_seq14"],
        "Start": [170, 0],
        "End": [2999, 6294],
        "Strand": ["+", "+"],
        # "PositionFile": [
        #     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.positions.csv",
        #     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.positions.csv",
        # ],
        "PositionFile": [
            "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.positions.csv.gz",
            "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.positions.csv.gz",
        ],
    }
)

pacbio_data_df.insert(0, condition_col, conditions[0])

pacbio_data_df

# %% papermill={"duration": 0.071769, "end_time": "2022-02-01T09:42:43.049672", "exception": false, "start_time": "2022-02-01T09:42:42.977903", "status": "completed"} tags=["parameters"]
illumina_data_df = pd.DataFrame(
    {
        # condition_col: [
        #     "RUSC2_MOUSE",
        #     "TRIM2_BOVIN",
        #     "CA2D3_MOUSE",
        #     "ABL_DROME",
        #     "DGLA_HUMAN",
        #     "K0513_MOUSE",
        #     "KCNAS_DROME",
        #     "ACHA4_MOUSE",
        #     "ANR17_HUMAN",
        #     "TWK7_CAEEL",
        #     "SCN1_HETBL",
        #     "CACB2_RABIT",
        #     "RIMS2_RAT",
        #     "PCLO_CHICK",
        #     "DOP1_HUMAN",
        #     "IQEC1_HUMAN",
        #     "CSKI1_MOUSE",
        #     "MTUS2_HUMAN",
        #     "ROBO2_HUMAN",
        # ],
        "Chrom": [
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
        ],
        "Start": [
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
        ],
        "End": [
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
        ],
        "Strand": [
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
        ],
    }
)

illumina_data_df.insert(0, condition_col, conditions[1])

illumina_data_df["PositionFile"] = [
    f"/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.{chrom}.positions.csv"
    for chrom in illumina_data_df["Chrom"]
]

illumina_data_df

# %% jupyter={"source_hidden": true}
# condition = "GRIA"
# strand = "+"
# positions_file = "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.positions.csv"

# positions_df = pd.read_csv(positions_file, sep=sep).drop(
#     ["CDS", "Phred"], axis=1, errors="ignore"
# )  # ignore error from droping these cols which only exist in Illumina data
# positions_df.insert(0, condition_col, condition)

# positions_df

# %% jupyter={"source_hidden": true}
# ref_base = "A" if strand == "+" else "T"

# editing_sites_df = (
#     positions_df.loc[
#         (positions_df["RefBase"] == ref_base) & (positions_df["Edited"] > 0),
#         [condition_col, "Chrom", "Position"]
#     ]
#     .sort_values(["Chrom", "Position"])
#     .reset_index(drop=True)
#     .rename(columns={"Position": "Start"})
# )

# editing_sites_df

# %% jupyter={"source_hidden": true}
# editing_sites_df.insert(
#     editing_sites_df.columns.get_loc("Start") + 1,
#     "End",
#     editing_sites_df["Start"] + 1
# )
# editing_sites_df

# %% jupyter={"source_hidden": true}
# editing_sites_df.insert(
#     editing_sites_df.columns.get_loc("End") + 1,
#     "Score",
#     "."
# )
# editing_sites_df.insert(
#     editing_sites_df.columns.get_loc("Score") + 1,
#     "Strand",
#     strand
# )
# editing_sites_df.insert(
#     editing_sites_df.columns.get_loc("End") + 1,
#     "Name",
#     (
#         editing_sites_df[condition_col]
#         + ":"
#         + editing_sites_df["Chrom"]
#         + ":"
#         + editing_sites_df["Start"].astype(str)
#         + ":"
#         + editing_sites_df["End"].astype(str)
#         + ":"
#         + editing_sites_df["Strand"].astype(str)
#     )
# )
# editing_sites_df

# %%
fasta_dict = make_fasta_dict(transcriptome_file)
chromosomal_lengths = {chrom: len(seq) for chrom, seq in fasta_dict.items()}

data_dfs = [pacbio_data_df, illumina_data_df]

merged_platforms_editing_sites_dfs = []

for data_df in data_dfs:
    platform_editing_sites_dfs = []

    for _, row in data_df.iterrows():  # _ is the row's index
        condition, _, _, _, strand, positions_file = row  # unpack all cols at once

        positions_df = pd.read_csv(positions_file, sep=sep).drop(
            ["CDS", "Phred"], axis=1, errors="ignore"
        )  # ignore error from droping these cols which only exist in Illumina data
        positions_df.insert(0, condition_col, condition)

        ref_base = "A" if strand == "+" else "T"

        editing_sites_df = (
            positions_df.loc[
                (positions_df["RefBase"] == ref_base) & (positions_df["Edited"] > 0),
                [condition_col, "Chrom", "Position"],
            ]
            .sort_values(["Chrom", "Position"])
            .reset_index(drop=True)
            .rename(columns={"Position": "Start"})
        )
        editing_sites_df.insert(
            editing_sites_df.columns.get_loc("Start") + 1,
            "End",
            editing_sites_df["Start"] + 1,
        )
        editing_sites_df.insert(
            editing_sites_df.columns.get_loc("End") + 1, "Score", "."
        )
        editing_sites_df.insert(
            editing_sites_df.columns.get_loc("Score") + 1, "Strand", strand
        )
        editing_sites_df.insert(
            editing_sites_df.columns.get_loc("End") + 1,
            "Name",
            (
                editing_sites_df[condition_col]
                + ":"
                + editing_sites_df["Chrom"]
                + ":"
                + editing_sites_df["Start"].astype(str)
                + ":"
                + editing_sites_df["End"].astype(str)
                + ":"
                + editing_sites_df["Strand"].astype(str)
            ),
        )

        # editing_sites_df = (
        #     positions_df.loc[
        #         (positions_df["RefBase"] == ref_base) & (positions_df["Edited"] > 0),
        #         [condition_col, "Chrom", "Position"]
        #     ]
        #     .sort_values(["Chrom", "Position"])
        #     .reset_index(drop=True)
        #     .rename(columns={"Position": "End"})
        # )
        # editing_sites_df.insert(
        #     editing_sites_df.columns.get_loc("Chrom") + 1,
        #     "Start",
        #     editing_sites_df["End"] - 1
        # )
        # editing_sites_df.insert(
        #     editing_sites_df.columns.get_loc("End") + 1,
        #     "Score",
        #     "."
        # )
        # editing_sites_df.insert(
        #     editing_sites_df.columns.get_loc("Score") + 1,
        #     "Strand",
        #     strand
        # )
        # editing_sites_df.insert(
        #     editing_sites_df.columns.get_loc("End") + 1,
        #     "Name",
        #     (
        #         editing_sites_df[condition_col]
        #         + ":"
        #         + editing_sites_df["Chrom"]
        #         + ":"
        #         + editing_sites_df["Start"].astype(str)
        #         + ":"
        #         + editing_sites_df["End"].astype(str)
        #         + ":"
        #         + editing_sites_df["Strand"].astype(str)
        #     )
        # )

        platform_editing_sites_dfs.append(editing_sites_df)

    merged_platform_editing_sites_df = pd.concat(
        platform_editing_sites_dfs, ignore_index=True
    )

    # extend start and end positions s.t. each region spans 3 bases, with the edited adenosine in the middle
    merged_platform_editing_sites_df["Start"] = (
        merged_platform_editing_sites_df["Start"] - 1
    )  # upstream base of edited adenosine (or downstream for a transcript expressed from the negative strand)
    merged_platform_editing_sites_df["End"] = (
        merged_platform_editing_sites_df["End"] + 1
    )  # downstream base of edited adenosine (or upstream for a transcript expressed from the negative strand)

    # don't consider editing sites located at the first/last position of a transcript (if there are such editing sites, they are negligble)
    merged_platform_editing_sites_df = merged_platform_editing_sites_df.loc[
        merged_platform_editing_sites_df["Start"] >= 0
    ]
    merged_platform_editing_sites_df = merged_platform_editing_sites_df.loc[
        merged_platform_editing_sites_df.apply(
            lambda x: x["End"] <= chromosomal_lengths[x["Chrom"]], axis=1
        )
    ]

    merged_platforms_editing_sites_dfs.append(merged_platform_editing_sites_df)

merged_platforms_editing_sites_dfs[0]

# %% jupyter={"source_hidden": true}
# merged_platforms_editing_sites_dfs[1]

# %%
merged_editing_sites_df = pd.concat(
    merged_platforms_editing_sites_dfs, ignore_index=True
)

# remove duplicates editing sites from transcripts present in both platforms (namely PCLO)
merged_editing_sites_df = merged_editing_sites_df.drop_duplicates(
    ["Chrom", "Start", "End", "Strand"], ignore_index=True
)

merged_editing_sites_df

# %%
editing_sites_dfs = merged_platforms_editing_sites_dfs + [merged_editing_sites_df]

# %%
editing_sites_bedtools = [
    (
        BedTool()
        .from_dataframe(editing_sites_df.iloc[:, 1:])
        .sort()
        .sequence(fi=transcriptome_file, s=True)
    )
    for editing_sites_df in editing_sites_dfs
]

# %%
fasta_files = [
    editing_sites_bedtool.seqfn for editing_sites_bedtool in editing_sites_bedtools
]
main_title = None
sub_titles = conditions + ["All"]

# %%
out_file = "ADAR motif of all currently-edited sites - Squid.svg"

# %%

# %%
import matplotlib.pyplot as plt  # matplotlib
import pandas as pd  # pandas
from Bio import SeqIO, motifs  # biopython
from logomaker import Logo  # logomaker

# %%
records = SeqIO.parse(fasta_files[0], "fasta")
seqs = [record.seq for record in records]

# %%
len(seqs)

# %%

# %%

# %%
# _sub_titles = [f"Squid's {sub_title}" for sub_title in sub_titles]
# multiple_logos_from_fasta_files(
#     fasta_files, main_title, _sub_titles, out_file, width=14, height=4, dpi=300
# );

# %%
_sub_titles = [f"Squid's {sub_title}" for sub_title in ["Long-reads", "Short-reads", "pooled data"]]
multiple_logos_from_fasta_files(
    fasta_files, main_title, _sub_titles, out_file, width=14, height=4, dpi=300
);

# %%
_sub_titles_2 = ["Long-reads", "Short-reads"]
out_file_2 = "ADAR motif of all currently-edited sites -  PacBio vs Illumina - Squid.svg"
fasta_files_2 = fasta_files[:2]
main_title_2 = "Squid"
multiple_logos_from_fasta_files(
    fasta_files_2, main_title_2, _sub_titles_2, out_file_2, width=0.66*14, height=4, dpi=300
);

# %% [markdown]
# # Pooled editing levels comparisons - squid

# %%
known_sites_file = (
    "/private7/projects/Combinatorics/D.pealeii/Annotations/D.pea.EditingSites.csv"
)

# pacbio_positions_files = [
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.positions.csv.gz",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.positions.csv.gz",
# ]
pacbio_positions_files = [
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.positions.csv.gz",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.positions.csv.gz",
]
illumina_positions_files = [
    f"/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.{chrom}.positions.csv"
    for chrom in [
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
]
positions_files = pacbio_positions_files + illumina_positions_files

sep = "\t"

pacbio_strands = [
    "+",
    "+",
]
illumina_strands = [
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
strands = pacbio_strands + illumina_strands

# condition_cols = ["Platform", "Transcript"]

condition_col = "Transcript-Platform"

pacbio_transcripts = ["GRIA", "PCLO"]
illumina_transcripts = [
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
transcripts = pacbio_transcripts + illumina_transcripts

platforms = ["PacBio"] * len(pacbio_transcripts) + ["Illumina"] * len(
    illumina_transcripts
)

conditions = [
    f"{transcript}-{platform}" for transcript, platform in zip(transcripts, platforms)
]

# conditions = ["PacBio", "Illumina"]

# chrom = "comp141882_c0_seq14"
# colors = [
#     px.colors.qualitative.G10[1],
#     px.colors.qualitative.Dark24[13],
#     # px.colors.qualitative.Set3[11]
#     px.colors.qualitative.Vivid[5],
# ]
# start = 0
# end = 6294

# %%
positions_dfs = [
    pd.read_csv(position_file, sep=sep) for position_file in positions_files
]
for positions_df, condition in zip(positions_dfs, conditions):
    positions_df.insert(0, condition_col, condition)
positions_dfs[0]


# %%
positions_dfs[5]

# %%
len(positions_dfs)

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

merged_ref_base_positions_df["Transcript-Platform"] = merged_ref_base_positions_df["Transcript-Platform"].apply(
    lambda x: "GRIA2-PacBio" if x == "GRIA-PacBio" else x
)

merged_ref_base_positions_df


# %% [markdown]
# ## Current vs. known editing levels

# %%
both_df = merged_ref_base_positions_df.loc[
    merged_ref_base_positions_df["Edited"]
    & merged_ref_base_positions_df["KnownEditing"]
]

all_df = merged_ref_base_positions_df.loc[
    merged_ref_base_positions_df["Edited"]
    | merged_ref_base_positions_df["KnownEditing"]
]
all_df


# %%
def editing_status_color(
    edited: bool,
    known_editing: bool,
    both_color="black",
    edited_color="red",
    known_editing_color="green",
):
    if edited and known_editing:
        return both_color
    elif edited:
        return edited_color
    else:
        return known_editing_color

# %% jupyter={"source_hidden": true} papermill={"duration": 4.052404, "end_time": "2022-02-01T09:42:53.176715", "exception": false, "start_time": "2022-02-01T09:42:49.124311", "status": "completed"}
# # todo retain nan rows and turn nans to 0?

# both_df = merged_ref_base_positions_df.loc[
#     merged_ref_base_positions_df["Edited"]
#     & merged_ref_base_positions_df["KnownEditing"]
# ]
# edited_df = merged_ref_base_positions_df.loc[
#     merged_ref_base_positions_df["Edited"]
#     & ~merged_ref_base_positions_df["KnownEditing"]
# ]
# known_editing_df = merged_ref_base_positions_df.loc[
#     ~merged_ref_base_positions_df["Edited"]
#     & merged_ref_base_positions_df["KnownEditing"]
# ]

# # all_df = merged_ref_base_positions_df.loc[merged_ref_base_positions_df["Edited"] | merged_ref_base_positions_df["KnownEditing"]]

# fig = make_subplots(
#     rows=1,
#     cols=1,
#     # subplot_titles=conditions,
#     # shared_yaxes=True,
#     # shared_xaxes=True,
#     x_title="% editing",
#     y_title="% known editing",
# )

# # df = merged_ref_base_positions_df.copy()
# # df = merged_ref_base_positions_df.loc[merged_ref_base_positions_df["Edited"] & merged_ref_base_positions_df["KnownEditing"]]
# x = all_df["%Editing"].fillna(0)
# y = all_df["%EditingKnown"].fillna(0)

# all_colors = all_df.apply(
#     lambda x: editing_status_color(x["Edited"], x["KnownEditing"]), axis=1
# )

# x_both = both_df["%Editing"].fillna(0)
# y_both = both_df["%EditingKnown"].fillna(0)
# r, pv = scipy.stats.pearsonr(x_both, y_both)

# fig.add_trace(
#     go.Scatter(
#         x=x,
#         y=y,
#         # name=condition,
#         mode="markers",
#         # marker_color="black",
#         marker_color=all_colors,
#         # marker_line_color=colors,
#         # marker_line=dict(
#         #     color=colors
#         # ),
#         # marker_color="white",
#         marker_size=6,
#         # marker_symbol="circle-open"
#     ),
#     row=1,
#     col=1,
# )

# fig.add_annotation(
#     row=1,
#     col=1,
#     x=20,
#     y=90,
#     xref="x",
#     yref="y",
#     text=f"<b>Pearson's r</b><br>p-val = {pv:.2e}<br>ρ = {r:.2g}",
#     bgcolor="white",
#     borderpad=4,
#     # opacity=0.7,
#     showarrow=False,
# )

# fig.update_layout(
#     # title_text="Correlation between current & previously-reported editing levels",
#     title_text="Correlation between current & previously-reported editing levels",
#     showlegend=False,
#     template=template,
#     width=600,
#     height=500,
# )

# fig.update_xaxes(range=[0, 100])
# fig.update_yaxes(range=[0, 100])

# fig.show()


# %% papermill={"duration": 4.052404, "end_time": "2022-02-01T09:42:53.176715", "exception": false, "start_time": "2022-02-01T09:42:49.124311", "status": "completed"}
# todo retain nan rows and turn nans to 0?

both_df = merged_ref_base_positions_df.loc[
    merged_ref_base_positions_df["Edited"]
    & merged_ref_base_positions_df["KnownEditing"]
]
edited_df = merged_ref_base_positions_df.loc[
    (merged_ref_base_positions_df["Edited"])
    & (~merged_ref_base_positions_df["KnownEditing"])
]
known_editing_df = merged_ref_base_positions_df.loc[
    (~merged_ref_base_positions_df["Edited"])
    & (merged_ref_base_positions_df["KnownEditing"])
]
all_df = merged_ref_base_positions_df.loc[
    merged_ref_base_positions_df["Edited"]
    | merged_ref_base_positions_df["KnownEditing"]
]

assert len(both_df) + len(edited_df) + len(known_editing_df) == len(all_df)


# %%
len(both_df)

# %%
len(edited_df)

# %%
len(known_editing_df)

# %%
# edited_df["Transcript-Platform"] = edited_df["Transcript-Platform"].apply(
#     lambda x: "GRIA2-PacBio" if x == "GRIA-PacBio" else x
# )
edited_df

# %%

# %% papermill={"duration": 4.052404, "end_time": "2022-02-01T09:42:53.176715", "exception": false, "start_time": "2022-02-01T09:42:49.124311", "status": "completed"}
fig = make_subplots(
    rows=1,
    cols=1,
    # x_title="% current editing",
    # y_title="% known editing",
    x_title="Editing level (targeted sequencing) [%]",
    y_title="Previously reported editing level [%]"
)

# all_colors = all_df.apply(
#     lambda x: editing_status_color(x["Edited"], x["KnownEditing"]),
#     axis=1
# )

x_both = both_df["%Editing"].fillna(0)
y_both = both_df["%EditingKnown"].fillna(0)

x_edited = edited_df["%Editing"].fillna(0)
y_edited = edited_df["%EditingKnown"].fillna(0)

x_known_editing = known_editing_df["%Editing"].fillna(0)
y_known_editing = known_editing_df["%EditingKnown"].fillna(0)

xs_ys = [(x_both, y_both), (x_edited, y_edited), (x_known_editing, y_known_editing)]
# names = ["Both", "Currently edited", "Known editing"]
names = ["Both", "De-novo", "Known"]
colors = ["black", "green", "red"]

for (x, y), name, color in zip(xs_ys, names, colors):
    fig.add_trace(
        go.Scatter(
            x=x,
            y=y,
            name=name,
            mode="markers",
            marker_color=color,
            marker_size=6,
        ),
        row=1,
        col=1,
    )

# x_all = all_df["%Editing"].fillna(0)
# y_all = all_df["%EditingKnown"].fillna(0)

pearson_df = all_df.loc[:, ["%Editing", "%EditingKnown"]].reset_index(drop=True).dropna()
x_all = pearson_df["%Editing"]
y_all = pearson_df["%EditingKnown"]

r, pv = scipy.stats.pearsonr(x_all, y_all)

if pv > 0:
    pearson_text = text=f"Pearson's r = {r:.2g}<br>p-value = {pv:.2g}"
else:
    # guarnteed as explained here: https://stackoverflow.com/questions/45914221/minimal-p-value-for-scipy-stats-pearsonr
    pearson_text=f"Pearson's r = {r:.2g}<br>p-value < 1E-22"


fig.add_annotation(
    row=1,
    col=1,
    x=20,
    # y=85,
    y=90,
    xref="x",
    yref="y",
    # text=f"<b>Pearson's r</b><br>p-val = {pv:.2e}<br>ρ = {r:.2g}",
    text=pearson_text,
    bgcolor="white",
    borderpad=4,
    # opacity=0.7,
    showarrow=False,
)

fig.update_layout(
    # title_text="Correlation between current & previously-reported editing levels",
    title_text="Pooled editing levels in squid",
    title_x=0.15,
    # legend_title_text="Editing status",
    legend_title_text="Detected",
    # showlegend=False,
    template=template,
    width=600,
    height=500,
)

fig.update_xaxes(range=[0, 100])
fig.update_yaxes(range=[0, 100])

fig.write_image(
    "Correlation between current vs previously-reported editing levels - Squid.svg",
    width=600,
    height=500,
)

fig.show()


# %% papermill={"duration": 4.052404, "end_time": "2022-02-01T09:42:53.176715", "exception": false, "start_time": "2022-02-01T09:42:49.124311", "status": "completed"}
fig = make_subplots(
    rows=1,
    cols=1,
    x_title="Detected",
    y_title="Editing level [%]",
)


y_both = both_df["%Editing"].fillna(0)
y_edited = edited_df["%Editing"].fillna(0)
y_known_editing = known_editing_df["%EditingKnown"].fillna(0)
ys = [y_both, y_edited, y_known_editing]

# names = ["Both", "Currently edited", "Known editing"]
names = ["Both", "De-novo", "Known"]
colors = ["black", "green", "red"]

for name, y, color in zip(names, ys, colors):
    fig.add_trace(
        go.Box(
            x=[name] * len(y),
            y=y,
            boxpoints="all",
            # jitter=0.8,
            whiskerwidth=0.2,
            marker_size=2,
            line_width=1,
            # fillcolor=color,
            line_color=color,
            # points="all"
            # name=name,
            # mode="markers",
            marker_color=color,
            # marker_size=6,
        ),
        row=1,
        col=1,
    )


fig.update_layout(
    # title_text="Correlation between current & previously-reported editing levels",
    title_text="Pooled editing levels in squid",
    title_x=0.15,
    showlegend=False,
    template=template,
    width=600,
    height=500,
)
# fig.update_traces(
#     boxpoints='all',
#     # jitter=0
# )

# fig.update_xaxes(range=[0, 100])
fig.update_yaxes(range=[0, 100])

fig.write_image(
    "Distribution of pooled editing levels according to their editing status in current vs. previous study - Squid.svg",
    width=600,
    height=500,
)

fig.show()


# %% [markdown]
# # Isoforms in 80K coverage - squid

# %%
condition_col = "Gene"
sep = "\t"
unique_proteins_first_col_pos = 14
# platforms = ["PacBio", "Illumina"]
platforms = ["Long-reads", "Short-reads"]

# %%
# pacbio_distinct_unique_proteins_files = [
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/GRIA-CNS-RESUB.DistinctUniqueProteins.03.03.2023-15:36:38.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/PCLO-CNS-RESUB.DistinctUniqueProteins.03.03.2023-15:53:01.csv",
# ]
pacbio_distinct_unique_proteins_files = [
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/GRIA-CNS-RESUB.DistinctUniqueProteins.06.02.2024-09:29:20.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO-CNS-RESUB.DistinctUniqueProteins.06.02.2024-09:46:24.csv",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples.80K/ADAR1.Merged.r64296e203404D01.aligned.sorted.Sampled80000.unique_proteins.csv.gz.DistinctUniqueProteins.04.02.2026-15:16:16.csv",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples.80K/IQEC.Merged.r64296e203404D01.aligned.sorted.Sampled80000.unique_proteins.csv.gz.DistinctUniqueProteins.04.02.2026-14:00:24.csv"
]
pacbio_unique_reads_files = [
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_reads.csv.gz",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_reads.csv.gz",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples.80K/ADAR1.Merged.r64296e203404D01.aligned.sorted.Sampled80000.unique_reads.csv.gz",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples.80K/IQEC.Merged.r64296e203404D01.aligned.sorted.Sampled80000.unique_reads.csv.gz"
]
pacbio_unique_proteins_files = [
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples.80K/ADAR1.Merged.r64296e203404D01.aligned.sorted.Sampled80000.unique_proteins.csv.gz",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples.80K/IQEC.Merged.r64296e203404D01.aligned.sorted.Sampled80000.unique_proteins.csv.gz"
]

# pacbio_conditions = ["GRIA", "PCLO"]
pacbio_conditions = [
    "GRIA2", 
    "PCLO",
    "ADAR1",
    "IQEC1"
]

pacbio_color_sequence = px.colors.qualitative.G10
pacbio_color_discrete_map = {
    condition: color
    for condition, color in zip(pacbio_conditions, pacbio_color_sequence)
}

# %%
# old sampling where there are 40K reads pairs instead of 80K

# illumina_distinct_unique_proteins_files = [
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80K/comp141881_c0_seq3.Sampled80000.DistinctUniqueProteins.11.04.2023-16:36:40.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80K/comp141044_c0_seq2.Sampled80000.DistinctUniqueProteins.11.04.2023-17:21:47.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80K/comp140439_c0_seq1.Sampled80000.DistinctUniqueProteins.11.04.2023-16:34:07.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80K/comp126362_c0_seq1.Sampled80000.DistinctUniqueProteins.11.04.2023-17:27:11.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80K/comp141517_c0_seq1.Sampled80000.DistinctUniqueProteins.11.04.2023-18:09:19.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80K/comp141840_c0_seq2.Sampled80000.DistinctUniqueProteins.11.04.2023-16:40:11.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80K/comp141640_c0_seq1.Sampled80000.DistinctUniqueProteins.11.04.2023-16:35:41.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80K/comp140987_c3_seq1.Sampled80000.DistinctUniqueProteins.11.04.2023-16:51:06.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80K/comp140910_c2_seq1.Sampled80000.DistinctUniqueProteins.11.04.2023-16:46:22.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80K/comp136058_c0_seq1.Sampled80000.DistinctUniqueProteins.11.04.2023-17:00:54.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80K/comp141378_c0_seq7.Sampled80000.DistinctUniqueProteins.11.04.2023-18:02:54.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80K/comp141158_c1_seq2.Sampled80000.DistinctUniqueProteins.11.04.2023-16:38:37.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80K/comp140712_c0_seq3.Sampled80000.DistinctUniqueProteins.11.04.2023-16:55:55.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80K/comp141882_c0_seq14.Sampled80000.DistinctUniqueProteins.11.04.2023-17:08:14.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80K/comp141880_c1_seq3.Sampled80000.DistinctUniqueProteins.11.04.2023-17:47:57.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80K/comp141565_c6_seq3.Sampled80000.DistinctUniqueProteins.11.04.2023-17:55:20.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80K/comp141684_c0_seq1.Sampled80000.DistinctUniqueProteins.11.04.2023-17:41:03.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80K/comp141532_c3_seq11.Sampled80000.DistinctUniqueProteins.11.04.2023-17:37:06.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80K/comp141574_c0_seq3.Sampled80000.DistinctUniqueProteins.11.04.2023-17:13:57.csv",
# ]

# illumina_chroms = [
#     "comp141881_c0_seq3",
#     "comp141044_c0_seq2",
#     "comp140439_c0_seq1",
#     "comp126362_c0_seq1",
#     "comp141517_c0_seq1",
#     "comp141840_c0_seq2",
#     "comp141640_c0_seq1",
#     "comp140987_c3_seq1",
#     "comp140910_c2_seq1",
#     "comp136058_c0_seq1",
#     "comp141378_c0_seq7",
#     "comp141158_c1_seq2",
#     "comp140712_c0_seq3",
#     "comp141882_c0_seq14",
#     "comp141880_c1_seq3",
#     "comp141565_c6_seq3",
#     "comp141684_c0_seq1",
#     "comp141532_c3_seq11",
#     "comp141574_c0_seq3",
# ]
# illumina_unique_reads_files = [
#     f"/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80K/reads.sorted.aligned.filtered.{chrom}.Sampled80000.unique_reads.csv.gz"
#     for chrom in illumina_chroms
# ]
# illumina_unique_proteins_files = [
#     f"/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80K/reads.sorted.aligned.filtered.{chrom}.Sampled80000.unique_proteins.csv.gz"
#     for chrom in illumina_chroms
# ]

# # illumina_conditions = [
# #     "RUSC2_MOUSE",
# #     "TRIM2_BOVIN",
# #     "CA2D3_MOUSE",
# #     "ABL_DROME",
# #     "DGLA_HUMAN",
# #     "K0513_MOUSE",
# #     "KCNAS_DROME",
# #     "ACHA4_MOUSE",
# #     "ANR17_HUMAN",
# #     "TWK7_CAEEL",
# #     "SCN1_HETBL",
# #     "CACB2_RABIT",
# #     "RIMS2_RAT",
# #     "PCLO_CHICK",
# #     "DOP1_HUMAN",
# #     "IQEC1_HUMAN",
# #     "CSKI1_MOUSE",
# #     "MTUS2_HUMAN",
# #     "ROBO2_HUMAN"
# # ]

# illumina_conditions = [
#     "RUSC2",
#     "TRIM2",
#     "CA2D3",
#     "ABL",
#     "DGLA",
#     "K0513",
#     "KCNAS",
#     "ACHA4",
#     "ANR17",
#     "TWK7",
#     "SCN1",
#     "CACB2",
#     "RIMS2",
#     "PCLO",
#     "DOP1",
#     "IQEC1",
#     "CSKI1",
#     "MTUS2",
#     "ROBO2",
# ]

# illumina_color_sequence = px.colors.qualitative.Dark24
# illumina_color_discrete_map = {
#     condition: color
#     for condition, color in zip(illumina_conditions, illumina_color_sequence)
# }

# %%
illumina_distinct_unique_proteins_files = [
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80KPairs/comp141881_c0_seq3.Sampled160000.DistinctUniqueProteins.04.01.2026-13:20:40.csv",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80KPairs/comp141044_c0_seq2.Sampled160000.DistinctUniqueProteins.04.01.2026-14:51:11.csv",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80KPairs/comp140439_c0_seq1.Sampled160000.DistinctUniqueProteins.04.01.2026-13:15:56.csv",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80KPairs/comp126362_c0_seq1.Sampled160000.DistinctUniqueProteins.04.01.2026-15:01:44.csv",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80KPairs/comp141517_c0_seq1.Sampled160000.DistinctUniqueProteins.04.01.2026-16:23:30.csv",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80KPairs/comp141840_c0_seq2.Sampled160000.DistinctUniqueProteins.04.01.2026-13:26:32.csv",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80KPairs/comp141640_c0_seq1.Sampled160000.DistinctUniqueProteins.04.01.2026-13:19:05.csv",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80KPairs/comp140987_c3_seq1.Sampled160000.DistinctUniqueProteins.04.01.2026-13:46:32.csv",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80KPairs/comp140910_c2_seq1.Sampled160000.DistinctUniqueProteins.04.01.2026-13:38:31.csv",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80KPairs/comp136058_c0_seq1.Sampled160000.DistinctUniqueProteins.04.01.2026-14:05:46.csv",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80KPairs/comp141378_c0_seq7.Sampled160000.DistinctUniqueProteins.04.01.2026-16:10:40.csv",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80KPairs/comp141158_c1_seq2.Sampled160000.DistinctUniqueProteins.04.01.2026-13:23:51.csv",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80KPairs/comp140712_c0_seq3.Sampled160000.DistinctUniqueProteins.04.01.2026-13:55:22.csv",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80KPairs/comp141882_c0_seq14.Sampled160000.DistinctUniqueProteins.04.01.2026-14:21:47.csv",    
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80KPairs/comp141880_c1_seq3.Sampled160000.DistinctUniqueProteins.04.01.2026-15:42:45.csv",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80KPairs/comp141565_c6_seq3.Sampled160000.DistinctUniqueProteins.04.01.2026-15:55:45.csv",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80KPairs/comp141684_c0_seq1.Sampled160000.DistinctUniqueProteins.04.01.2026-15:29:17.csv",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80KPairs/comp141532_c3_seq11.Sampled160000.DistinctUniqueProteins.04.01.2026-15:21:31.csv",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80KPairs/comp141574_c0_seq3.Sampled160000.DistinctUniqueProteins.04.01.2026-14:32:56.csv",
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

illumina_unique_reads_files = [
    f"/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80KPairs/reads.sorted.aligned.filtered.{chrom}.Sampled160000.unique_reads.csv.gz"
    for chrom in illumina_chroms
]
illumina_unique_proteins_files = [
    f"/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80KPairs/reads.sorted.aligned.filtered.{chrom}.Sampled160000.unique_proteins.csv.gz"
    for chrom in illumina_chroms
]

illumina_conditions = [
    "RUSC2",
    "TRIM2",
    "CA2D3",
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

illumina_color_sequence = px.colors.qualitative.Dark24
illumina_color_discrete_map = {
    condition: color
    for condition, color in zip(illumina_conditions, illumina_color_sequence)
}

# %%
platforms_distinct_proteins_files = [
    pacbio_distinct_unique_proteins_files,
    illumina_distinct_unique_proteins_files,
]
platforms_unique_reads_files = [pacbio_unique_reads_files, illumina_unique_reads_files]
platforms_unique_proteins_files = [
    pacbio_unique_proteins_files,
    illumina_unique_proteins_files,
]
platforms_conditions = [pacbio_conditions, illumina_conditions]

editable_aas_per_platform_and_condition = {}

platforms_distinct_proteins_dfs = []

for (
    platform,
    platform_distinct_proteins_files,
    platform_unique_reads_files,
    platform_unique_proteins_files,
    platform_conditions,
) in zip(
    platforms,
    platforms_distinct_proteins_files,
    platforms_unique_reads_files,
    platforms_unique_proteins_files,
    platforms_conditions,
):
    distinct_unique_proteins_dfs = []
    for (
        condition,
        distinct_unique_proteins_file,
        unique_reads_file,
        unique_proteins_file,
    ) in zip(
        platform_conditions,
        platform_distinct_proteins_files,
        platform_unique_reads_files,
        platform_unique_proteins_files,
    ):
        # unique_reads_df
        unique_reads_df = pd.read_csv(unique_reads_file, sep=sep)
        if "Transcript" in unique_reads_df.columns:
            unique_reads_df.rename(columns={"Transcript": "UniqueRead"}, inplace=True)

        # unique proteins df --> editable AAs
        unique_proteins_df = pd.read_csv(unique_proteins_file, sep=sep)
        unique_proteins_df.rename(
            columns={
                col: col.replace("Transcripts", "UniqueReads")
                for col in unique_proteins_df.columns[:unique_proteins_first_col_pos]
                if "Transcripts" in col
            },
            inplace=True,
        )
        editable_aas_per_platform_and_condition[
            (platform, condition)
        ] = unique_proteins_df.iloc[:, unique_proteins_first_col_pos:].shape[1]

        distinct_unique_proteins_df = pd.read_csv(
            distinct_unique_proteins_file, sep=sep
        )
        # if platform == "PacBio":
        #     distinct_unique_proteins_df = distinct_unique_proteins_df.drop("AvailableReads", axis=1)
        if platform == "Long-reads":
            distinct_unique_proteins_df = distinct_unique_proteins_df.drop(
                "AvailableReads", axis=1
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
        distinct_unique_proteins_df.insert(0, "Platform", platform)
        distinct_unique_proteins_dfs.append(distinct_unique_proteins_df)

    distinct_unique_proteins_df = (
        pd.concat(distinct_unique_proteins_dfs)
        .reset_index(drop=True)
        .rename(
            columns={"NumUniqueSamples": "NumOfProteins", "UniqueSamples": "Proteins"}
        )
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
    platforms_distinct_proteins_dfs.append(distinct_unique_proteins_df)

distinct_proteins_df = pd.concat(platforms_distinct_proteins_dfs, ignore_index=True)
distinct_proteins_df


# %%
editable_aas_per_platform_and_condition

# %%
distinct_proteins_per_editable_aas_df = distinct_proteins_df.copy()

# distinct_proteins_per_editable_aas_df = distinct_proteins_per_editable_aas_df.loc[
#     distinct_proteins_per_editable_aas_df["Algorithm"] == "Descending"
# ]

distinct_proteins_per_editable_aas_df = distinct_proteins_per_editable_aas_df.loc[
    distinct_proteins_per_editable_aas_df["Fraction"] == 1.0
]

distinct_proteins_per_editable_aas_df["EditableAAs"] = [
    editable_aas_per_platform_and_condition[(platform, condition)]
    for platform, condition in zip(
        distinct_proteins_per_editable_aas_df["Platform"],
        distinct_proteins_per_editable_aas_df[condition_col],
    )
]

distinct_proteins_per_editable_aas_df = (
    distinct_proteins_per_editable_aas_df.sort_values(
        ["NumOfProteins", "EditableAAs"], ascending=False
    ).reset_index(drop=True)
)

distinct_proteins_per_editable_aas_df


# %%
(
    distinct_proteins_per_editable_aas_df.sort_values(
        ["Platform", "Gene", "NumOfProteins"], ascending=False
    ).drop_duplicates(["Platform", "Gene", "NumOfProteins"], ignore_index=True)
)

# %% jupyter={"source_hidden": true}
# platforms_color_map = {
#     platform: color_map
#     for platform, color_map in zip(
#         platforms, [pacbio_color_discrete_map, illumina_color_discrete_map]
#     )
# }

# # ["circle", "square-dot", "diamond", "circle", "star", "star-square", "triangle-down"]
# # symbols = ["diamond", "square"]
# # symbols = ["diamond", "circle"]
# symbols = ["star", "triangle-up"]
# platforms_symbols = {
#     platform: symbol for platform, symbol in zip(platforms, symbols)
# }


# fig = go.Figure()

# for platform, condition in editable_aas_per_platform_and_condition:

#     platform_and_condition_df = distinct_proteins_per_editable_aas_df.loc[
#         (distinct_proteins_per_editable_aas_df["Platform"] == platform)
#         & (distinct_proteins_per_editable_aas_df[condition_col] == condition)
#     ]

#     x = platform_and_condition_df["EditableAAs"]
#     y = platform_and_condition_df["NumOfProteins"]

#     x_mean = [x.iloc[0]]
#     y_mean = [y.mean()]

#     fig.add_trace(
#         go.Scatter(
#             x=x_mean,
#             y=y_mean,
#             mode="markers",
#             marker=dict(
#                 color=platforms_color_map[platform][condition],
#                 symbol=platforms_symbols[platform],
#                 size=8,
#             ),
#             legendgrouptitle_text=platform,
#             legendgroup=platform,  # this can be any string
#             name=condition,
#         )
#     )

# # correlate all x and y values
# x = distinct_proteins_per_editable_aas_df["EditableAAs"]
# y = distinct_proteins_per_editable_aas_df["NumOfProteins"]
# r, pv = scipy.stats.pearsonr(x, y)

# fig.add_annotation(
#     x=40,
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

# fig.update_xaxes(
#     range=[0, distinct_proteins_per_editable_aas_df["EditableAAs"].max() * 1.1]
# )
# fig.update_yaxes(
#     range=[0, distinct_proteins_per_editable_aas_df["NumOfProteins"].max() * 1.1]
# )

# fig.update_layout(
#     # legend_font=dict(size=8),
#     # legend_grouptitlefont=dict(size=10),
#     template=template,
#     xaxis_title="Editable amino acids",
#     yaxis_title="Distinct proteins (avg)",
#     width=700,
#     height=650,
# )

# fig.write_image(
#     "Distinct proteins vs. editable AAs (80K) - PacBio vs. Illumina.svg",
#     width=700,
#     height=650
# )

# fig.show()


# %%
platform = "PacBio"
condition = "PCLO"
platform_and_condition_df = distinct_proteins_per_editable_aas_df.loc[
    (distinct_proteins_per_editable_aas_df["Platform"] == platform)
    & (distinct_proteins_per_editable_aas_df[condition_col] == condition)
]
platform_and_condition_df

# %%
x = platform_and_condition_df["EditableAAs"]
y = platform_and_condition_df["NumOfProteins"]

# %%
y.max()

# %%
# y.argmax()

# %%
# y_max = [y.max()]
# x_corresponding_y_max = [x.iloc[y.argmax()]]

# %%
# x_mean = [x.iloc[0]]
# y_mean = [y.mean()]

# %%

# %%
platforms_color_map = {
    platform: color_map
    for platform, color_map in zip(
        platforms, [pacbio_color_discrete_map, illumina_color_discrete_map]
    )
}

# ["circle", "square-dot", "diamond", "circle", "star", "star-square", "triangle-down"]
# symbols = ["diamond", "square"]
# symbols = ["diamond", "circle"]
# symbols = ["star", "triangle-up"]
symbols = ["star", "circle"]
symbols_sizes = {"star": 12, "circle": 10}
platforms_symbols = {platform: symbol for platform, symbol in zip(platforms, symbols)}


fig = go.Figure()

for platform, condition in editable_aas_per_platform_and_condition:
    platform_and_condition_df = distinct_proteins_per_editable_aas_df.loc[
        (distinct_proteins_per_editable_aas_df["Platform"] == platform)
        & (distinct_proteins_per_editable_aas_df[condition_col] == condition)
    ].sort_values("NumOfProteins", ascending=False)

    x = platform_and_condition_df["EditableAAs"]
    y = platform_and_condition_df["NumOfProteins"]

    x_mean = [x.iloc[0]]
    y_mean = [y.mean()]

    y_max = [y.max()]
    x_corresponding_y_max = [x.iloc[y.argmax()]]

    name = "GRIA2" if condition == "GRIA" and platform == "Long-reads" else condition
    
    fig.add_trace(
        go.Scatter(
            x=x_corresponding_y_max,
            y=y_max,
            mode="markers",
            marker=dict(
                color=platforms_color_map[platform][condition],
                symbol=platforms_symbols[platform],
                # size=8,
                # size=12,
                size=symbols_sizes[platforms_symbols[platform]],
            ),
            legendgrouptitle_text=platform,
            legendgroup=platform,  # this can be any string
            name=name,
        )
    )

# correlate all x and y values
x = distinct_proteins_per_editable_aas_df["EditableAAs"]
y = distinct_proteins_per_editable_aas_df["NumOfProteins"]
r, pv = scipy.stats.pearsonr(x, y)

# pearson_text = f"<b>Pearson's r</b><br>p-val = {pv:.2e}<br>ρ = {r:.2g}"
if pv > 0:
    pearson_text = text=f"Pearson's r = {r:.2g}<br>p-value = {pv:.2g}"
else:
    # guarnteed as explained here: https://stackoverflow.com/questions/45914221/minimal-p-value-for-scipy-stats-pearsonr
    pearson_text=f"Pearson's r = {r:.2g}<br>p-value < 1E-22"

fig.add_annotation(
    # x=27,
    # y=18_000,
    x=110,
    y=2_500,
    xref="x",
    yref="y",
    text=pearson_text,
    bgcolor="white",
    borderpad=10,
    font=dict(size=14),
    opacity=0.8,
    showarrow=False,
)

fig.update_xaxes(
    range=[0, distinct_proteins_per_editable_aas_df["EditableAAs"].max() * 1.1]
)
fig.update_yaxes(
    range=[0, distinct_proteins_per_editable_aas_df["NumOfProteins"].max() * 1.1]
)

width = 750
# height = width * 650 / 700
height = 700


fig.update_layout(
    # legend_font=dict(size=8),
    # legend_grouptitlefont=dict(size=10),
    # title="Squid",
    # title_x=0.15,
    template=template,
    xaxis_title="Editable codons",
    yaxis_title="Distinct protein isoforms",
    width=width,
    height=height,
)

fig.write_image(
    Path(out_dir, "Distinct proteins vs. editable AAs (80K) - PacBio vs. Illumina.svg"),
    width=width,
    height=height,
)

fig.show()


# %% [markdown]
# # Combined expression plots for squid

# %%
condition_col = "Gene"
platforms = ["Long-reads", "Short-reads"]

pacbio_conditions = ["GRIA2", "PCLO", "ADAR1", "IQEC1"]

illumina_conditions = [
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
illumina_conditions = [condition.split("_")[0] for condition in illumina_conditions]

pacbio_color_sequence = px.colors.qualitative.G10
pacbio_color_discrete_map = {
    condition: color
    for condition, color in zip(pacbio_conditions, pacbio_color_sequence)
}

illumina_color_sequence = px.colors.qualitative.Dark24
illumina_color_discrete_map = {
    condition: color
    for condition, color in zip(illumina_conditions, illumina_color_sequence)
}

# %%
in_dir = Path("/private7/projects/Combinatorics/Code/Notebooks")
out_dir = in_dir

illumina_merged_assignment_file = Path(in_dir, "AssignedExpression.Illumina.tsv")
pacbio_merged_assignment_file = Path(in_dir, "AssignedExpression.PacBio.tsv")
pacbio_umi_merged_assignment_file = Path(
    in_dir, "AssignedExpression.PacBio.WithUMIs.tsv"
)

# %%
expression_cols = [
    "Gene",
    "Platform",
    "#Protein",
    "%RelativeExpression",
    "%CummulativeRelativeExpression",
]

# %%
illumina_merged_assignment_df = pd.read_table(
    illumina_merged_assignment_file, usecols=expression_cols
)
illumina_merged_assignment_df[condition_col] = (
    illumina_merged_assignment_df[condition_col].str.split("_").str[0]
)
# illumina_merged_assignment_df

# pacbio_merged_assignment_df = pd.read_table(pacbio_merged_assignment_file)
pacbio_merged_assignment_df = pd.concat(
    [
        pd.read_table(pacbio_merged_assignment_file, usecols=expression_cols),
        pd.read_table(pacbio_umi_merged_assignment_file, usecols=expression_cols),
    ]
)
pacbio_merged_assignment_df[condition_col] = pacbio_merged_assignment_df[
    condition_col
].apply(lambda x: "GRIA2" if x == "GRIA" else x)
pacbio_merged_assignment_df

# %%
illumina_assignment_dfs = [
    illumina_merged_assignment_df.loc[
        illumina_merged_assignment_df[condition_col] == condition
    ].reset_index(drop=True)
    for condition in illumina_conditions
]
ic(len(illumina_assignment_dfs))
illumina_assignment_dfs[0]

# %%
pacbio_assignment_dfs = [
    pacbio_merged_assignment_df.loc[
        pacbio_merged_assignment_df[condition_col] == condition
    ].reset_index(drop=True)
    for condition in pacbio_conditions
]
ic(len(pacbio_assignment_dfs))
pacbio_assignment_dfs[0]

# %%
platforms_color_map = {
    platform: color_map
    for platform, color_map in zip(
        platforms, [pacbio_color_discrete_map, illumina_color_discrete_map]
    )
}

# dashes = ["dash", "dot"]
dashes = ["dash", "solid"]
platforms_dashes = {platform: dash for platform, dash in zip(platforms, dashes)}

symbols = ["star", "triangle-up"]
platforms_symbols = {platform: symbol for platform, symbol in zip(platforms, symbols)}


joined_platforms = [platforms[0]] * len(pacbio_conditions) + [platforms[1]] * len(
    illumina_conditions
)

joined_conditions = pacbio_conditions + illumina_conditions

joined_assignments_dfs = pacbio_assignment_dfs + illumina_assignment_dfs

# %% [markdown]
# ## All squid genes

# %%
fig = make_subplots(
    rows=1,
    cols=2,
    y_title="Relative expression [%]",
    x_title="Isoform rank",
    subplot_titles=platforms,
    shared_yaxes="all",
    shared_xaxes="all",
    # vertical_spacing=facet_row_spacing / 2.5,
    # horizontal_spacing=facet_col_spacing * 1.5,
    # vertical_spacing=0.05,
    vertical_spacing=0.03,
    horizontal_spacing=0.015,
)

for platform, condition, assignment_df in zip(
    joined_platforms, joined_conditions, joined_assignments_dfs
):
    x = assignment_df["#Protein"]
    y = assignment_df["%RelativeExpression"]

    color = platforms_color_map[platform][condition]
    col = 1 if platform == "Long-reads" else 2

    fig.add_trace(
        go.Scattergl(
            x=x,
            y=y,
            mode="markers",
            marker=dict(
                color=color,
                # symbol=platforms_symbols[platform],
                size=3,
                opacity=0.3,
            ),
            # mode="lines",
            # line=dict(
            #     color=platforms_color_map[platform][condition],
            #     dash=platforms_dashes[platform],
            #     # size=8,
            # ),
            # opacity=0.3,
            legendgrouptitle_text=platform,
            legendgroup=platform,  # this can be any string
            name=condition,
            showlegend=False,
        ),
        row=1,
        col=col,
    )

    # Add trace for legend so the legened wont be opaque
    fig.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            mode="markers",
            marker=dict(
                color=color,
                # symbol=platforms_symbols[platform],
                size=5,
                # size=3,
            ),
            legendgrouptitle_text=platform,
            legendgroup=platform,  # this can be any string
            name=condition,
            # name=f"{condition}  ",
            showlegend=True,
            # visible="legendonly",
        ),
        row=1,
        col=col,
    )

fig.update_xaxes(type="log")
fig.update_yaxes(
    type="log",
    # range=[np.log10(y_min), np.log10(y_max)], nticks=6
)
width = 1200
# height = 650
height = 700

fig.update_layout(
    legend_font=dict(size=10),
    # legend_grouptitlefont=dict(size=8),
    # legend_tracegroupgap=4,
    title="Squid",
    title_x=0.1,
    template=template,
    width=width,
    height=height,
)

fig.write_image(
    Path(
        out_dir,
        "Relative expression vs. distinct protein rank - PacBio vs. Illumina.svg",
    ),
    width=width,
    height=height,
)

fig.show()


# %%
Path(
    out_dir,
    "Relative expression vs. distinct protein rank - PacBio vs. Illumina.svg",
)

# %%
font_size = 28

fig = make_subplots(
    rows=1,
    cols=2,
    y_title="Cumulative relative<br>expression [%]",
    x_title="Isoform rank",
    subplot_titles=platforms,
    shared_yaxes="all",
    shared_xaxes="all",
    vertical_spacing=0.03,
    # horizontal_spacing=0.015,
    horizontal_spacing=0.03,
)

for platform, condition, assignment_df in zip(
    joined_platforms, joined_conditions, joined_assignments_dfs
):
    x = assignment_df["#Protein"]
    y = assignment_df["%CummulativeRelativeExpression"]

    color = platforms_color_map[platform][condition]
    col = 1 if platform == "Long-reads" else 2

    fig.add_trace(
        go.Scattergl(
            x=x,
            y=y,
            mode="markers",
            marker=dict(
                color=color,
                # symbol=platforms_symbols[platform],
                size=3,
                opacity=0.3,
            ),
            # mode="lines",
            # line=dict(
            #     color=platforms_color_map[platform][condition],
            #     dash=platforms_dashes[platform],
            #     # size=8,
            # ),
            # opacity=0.3,
            legendgrouptitle_text=platform,
            legendgroup=platform,  # this can be any string
            name=condition,
            showlegend=False,
        ),
        row=1,
        col=col,
    )

    # Add trace for legend so the legened wont be opaque
    fig.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            mode="markers",
            marker=dict(
                color=color,
                # symbol=platforms_symbols[platform],
                size=5,
            ),
            legendgrouptitle_text=platform,
            legendgroup=platform,  # this can be any string
            name=condition,
            showlegend=True,
            # visible="legendonly",
        ),
        row=1,
        col=col,
    )

fig.update_xaxes(
    type="log",
    tickfont=dict(size=0.7 * font_size),
    tickangle=20,
)
fig.update_yaxes(
    # type="log",
    tickfont=dict(size=0.7 * font_size),
    # range=[np.log10(y_min), np.log10(y_max)], nticks=6
)

fig.update_annotations(font_size=font_size)

width = 1200
height = 650

fig.update_layout(
    # legend_font=dict(size=0.5 * font_size),
    legend_font=dict(size=10),
    legend_grouptitlefont=dict(size=0.5 * font_size),
    # legend_font=dict(size=6),
    # legend_grouptitlefont=dict(size=8),
    # legend_tracegroupgap=4,
    # title="Squid",
    # title_x=0.1,
    template=template,
    width=width,
    height=height,
    margin=dict(l=120, r=40, t=60, b=80),
)

fig.write_image(
    Path(
        out_dir,
        "Weighted cumulative expression vs. distinct protein rank - PacBio vs. Illumina.svg",
    ),
    width=width,
    height=height,
)

fig.show()


# %%
font_size = 28

fig = make_subplots(
    rows=1,
    cols=2,
    y_title="Cumulative relative<br>expression [%]",
    x_title="Isoform rank",
    subplot_titles=platforms,
    shared_yaxes="all",
    shared_xaxes="all",
    vertical_spacing=0.03,
    # horizontal_spacing=0.015,
    horizontal_spacing=0.03,
)

for platform, condition, assignment_df in zip(
    joined_platforms, joined_conditions, joined_assignments_dfs
):
    x = assignment_df["#Protein"]
    y = assignment_df["%CummulativeRelativeExpression"]

    color = platforms_color_map[platform][condition]
    col = 1 if platform == "Long-reads" else 2

    fig.add_trace(
        go.Scattergl(
            x=x,
            y=y,
            mode="markers",
            marker=dict(
                color=color,
                # symbol=platforms_symbols[platform],
                size=3,
                opacity=0.3,
            ),
            # mode="lines",
            # line=dict(
            #     color=platforms_color_map[platform][condition],
            #     dash=platforms_dashes[platform],
            #     # size=8,
            # ),
            # opacity=0.3,
            legendgrouptitle_text=platform,
            legendgroup=platform,  # this can be any string
            name=condition,
            showlegend=False,
        ),
        row=1,
        col=col,
    )

    # Add trace for legend so the legened wont be opaque
    fig.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            mode="markers",
            marker=dict(
                color=color,
                # symbol=platforms_symbols[platform],
                size=5,
            ),
            legendgrouptitle_text=platform,
            legendgroup=platform,  # this can be any string
            name=condition,
            showlegend=True,
            # visible="legendonly",
        ),
        row=1,
        col=col,
    )

fig.update_xaxes(
    type="log",
    tickfont=dict(size=0.7 * font_size),
    tickangle=20,
)
fig.update_yaxes(
    # type="log",
    tickfont=dict(size=0.7 * font_size),
    # range=[np.log10(y_min), np.log10(y_max)], nticks=6
)

fig.update_annotations(font_size=font_size)

width = 1200
height = 650

fig.update_layout(
    # legend_font=dict(size=0.5 * font_size),
    legend_font=dict(size=10),
    legend_grouptitlefont=dict(size=0.5 * font_size),
    # legend_font=dict(size=6),
    # legend_grouptitlefont=dict(size=8),
    # legend_tracegroupgap=4,
    # title="Squid",
    # title_x=0.1,
    template=template,
    width=width,
    height=height,
    margin=dict(l=120, r=40, t=60, b=80),
)

# fig.write_image(
#     Path(
#         out_dir,
#         "Weighted cumulative expression vs. distinct protein rank - PacBio vs. Illumina.svg",
#     ),
#     width=width,
#     height=height,
# )

fig.show()

# %%
font_size = 36

fig = make_subplots(
    rows=2,
    cols=1,
    y_title="Cumulative relative expression [%]",
    x_title="Isoform rank",
    subplot_titles=platforms,
    shared_yaxes="all",
    shared_xaxes="all",
    vertical_spacing=0.03,
    # horizontal_spacing=0.015,
    horizontal_spacing=0.03,
)

for platform, condition, assignment_df in zip(
    joined_platforms, joined_conditions, joined_assignments_dfs
):
    x = assignment_df["#Protein"]
    y = assignment_df["%CummulativeRelativeExpression"]

    color = platforms_color_map[platform][condition]
    # col = 1 if platform == "Long-reads" else 2
    row = 1 if platform == "Long-reads" else 2

    fig.add_trace(
        go.Scattergl(
            x=x,
            y=y,
            mode="markers",
            marker=dict(
                color=color,
                # symbol=platforms_symbols[platform],
                size=3,
                opacity=0.3,
            ),
            # mode="lines",
            # line=dict(
            #     color=platforms_color_map[platform][condition],
            #     dash=platforms_dashes[platform],
            #     # size=8,
            # ),
            # opacity=0.3,
            # legendgrouptitle_text=platform,
            legendgrouptitle_text=platform + "  ",
            legendgroup=platform,  # this can be any string
            name=condition,
            showlegend=False,
        ),
        row=row,
        col=1,
    )

    # Add trace for legend so the legened wont be opaque
    fig.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            mode="markers",
            marker=dict(
                color=color,
                # symbol=platforms_symbols[platform],
                size=5,
            ),
            # legendgrouptitle_text=platform,
            legendgrouptitle_text=platform + "  ",
            legendgroup=platform,  # this can be any string
            name=condition,
            showlegend=True,
            # visible="legendonly",
        ),
        row=row,
        col=1,
    )

fig.update_xaxes(
    type="log",
    tickfont=dict(size=0.7 * font_size),
    # tickangle=20,
)
fig.update_yaxes(
    # type="log",
    tickfont=dict(size=0.7 * font_size),
    # range=[np.log10(y_min), np.log10(y_max)], nticks=6
)

fig.update_annotations(font_size=font_size)

# width = 650
width = 1200
height = 1200


fig.update_layout(
    # legend_font=dict(size=0.5 * font_size),
    legend_font=dict(size=18),
    # legend_grouptitlefont=dict(size=0.5 * font_size),
    legend_grouptitlefont=dict(size=22),
    # legend_font=dict(size=6),
    # legend_grouptitlefont=dict(size=8),
    # legend_tracegroupgap=4,
    # title="Squid",
    # title_x=0.1,
    template=template,
    width=width,
    height=height,
    margin=dict(l=120, r=60, t=60, b=80),
)

fig.write_image(
    Path(
        out_dir,
        "Weighted cumulative expression vs. distinct protein rank - vertical - PacBio vs. Illumina.svg",
    ),
    width=width,
    height=height,
)

fig.show()

# %% [markdown]
# ## Combined expression plots for long-reads squid genes

# %%
# TODO copy the code for the individual expression plots of previous long reads genes (GRIA2 & PCLO)
# and extend them with ADAR1 and IQEC1

# %% [markdown]
# ## Duped vs de-duped expression

# %%
# TODO reads dedopped assignment dfs and compare them to the non-dedopped ones (ADAR1 and IQEC1)

# %%
pacbio_umi_conditions = ["ADAR1", "IQEC1"]

umis_subcolors_discrete_map = {
    condition: two_subcolors_from_hex(pacbio_color_discrete_map[condition])
    for condition in pacbio_umi_conditions
}
umis_subcolors_discrete_map

# %%
pacbio_umi_unique_merged_assignment_file = Path(
    in_dir, "AssignedExpression.PacBio.WithUMIs.Unique.tsv"
)

# %%
merged_umi_duped_assignment_df = pd.read_table(
    pacbio_umi_merged_assignment_file, usecols=expression_cols
)
merged_umi_duped_assignment_df.insert(1, "Deduped", False)
# merged_umi_dooped_assignment_df

# %%
merged_umi_deduped_assignment_df = pd.read_table(
    pacbio_umi_unique_merged_assignment_file, usecols=expression_cols
)
merged_umi_deduped_assignment_df.insert(1, "Deduped", True)
# merged_umi_unqiue_assignment_df

# %%
merged_umi_assignment_df = pd.concat(
    [
        merged_umi_duped_assignment_df,
        merged_umi_deduped_assignment_df,
    ]
)
merged_umi_assignment_df

# %%
merged_umi_assignment_df.loc[
    :,
    [
        "Gene",
        "Deduped",
    ],
].value_counts().reset_index().sort_values("Gene")

# %%
x_axis_name = "Isoform rank"
y_axis_name = "Cumulative relative<br>expression [%]"
head_title = "Cumulative expression vs. distinct unique proteins"

cols = min(facet_col_wrap, len(pacbio_umi_conditions), 4)
rows = ceil(len(pacbio_umi_conditions) / cols)
row_col_iter = list(product(range(1, rows + 1), range(1, cols + 1)))[
    : len(pacbio_umi_conditions)
]

# assignment_method = "Weighted"
# percentile_fractions = [0.1, 1.0]

fig = make_subplots(
    rows=rows,
    cols=cols,
    subplot_titles=pacbio_umi_conditions,
    shared_xaxes=True,
    shared_yaxes=True,
    x_title=x_axis_name,
    y_title=y_axis_name,
    # vertical_spacing=0.01,
    horizontal_spacing=facet_col_spacing,
)

# legend_x = [5]
# legend_ys = [[75], [65]]
legend_x = [3]
# legend_ys = [[85], [75]]
legend_ys = [[95], [85]]
inner_texts = ["Full data", "Deduplicated data (UMIs)"]

for (row, col), condition in zip(row_col_iter, pacbio_umi_conditions):

    duped_aassignment_df = merged_umi_assignment_df.loc[
        (~merged_umi_assignment_df["Deduped"])
        & (merged_umi_assignment_df["Gene"] == condition)
    ]
    deduped_aassignment_df = merged_umi_assignment_df.loc[
        (merged_umi_assignment_df["Deduped"])
        & (merged_umi_assignment_df["Gene"] == condition)
    ]

    for color, df, legend_y, inner_text in zip(
        umis_subcolors_discrete_map[condition],
        [duped_aassignment_df, deduped_aassignment_df],
        legend_ys,
        inner_texts,
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
                    # opacity=0.7,
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
                    # opacity=0.7,
                    # symbol=symbol,
                    # line=dict(width=0),
                ),
                text=inner_text,
                textposition="middle right",
                textfont=dict(size=10),
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

width = 900
height = 450

fig.update_layout(
    # title="Squid's Long-reads",
    # title_x=0.11,
    showlegend=False,
    template=template,
    width=width,
    height=height,
)

fig.write_image(
    Path(
        out_dir,
        "Cumulative expression vs. distinct protein rank - PacBio - full data vs. deduped.svg",
    ),
    width=width,
    height=height,
)

fig.show()

# %%

# %% [markdown]
# # Combined noise plots

# %%
condition_col = "Gene"
platforms = ["Long-reads", "Short-reads"]

pacbio_conditions = ["GRIA2", "PCLO"]

illumina_conditions = [
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
illumina_conditions = [condition.split("_")[0] for condition in illumina_conditions]

pacbio_color_sequence = px.colors.qualitative.G10
pacbio_color_discrete_map = {
    condition: color
    for condition, color in zip(pacbio_conditions, pacbio_color_sequence)
}

illumina_color_sequence = px.colors.qualitative.Dark24
illumina_color_discrete_map = {
    condition: color
    for condition, color in zip(illumina_conditions, illumina_color_sequence)
}

# %%
illumina_merged_noise_file = "NoiseLevels.Illumina.tsv"
pacbio_merged_noise_file = "NoiseLevels.PacBio.tsv"
octopus_noise_file = "Noise.Octopus.tsv"

# %%
illumina_merged_noise_df = pd.read_table(illumina_merged_noise_file)
illumina_merged_noise_df[condition_col] = (
    illumina_merged_noise_df[condition_col].str.split("_").str[0]
)
# illumina_merged_assignment_df

pacbio_merged_noise_df = pd.read_table(pacbio_merged_noise_file)
# pacbio_merged_assignment_df

# %%
illumina_noise_dfs = [
    illumina_merged_noise_df.loc[
        illumina_merged_noise_df[condition_col] == condition
    ].reset_index(drop=True)
    for condition in illumina_conditions
]
illumina_noise_dfs[0]

# %%
pacbio_noise_dfs = [
    pacbio_merged_noise_df.loc[
        pacbio_merged_noise_df[condition_col] == condition
    ].reset_index(drop=True)
    for condition in pacbio_conditions
]
pacbio_noise_dfs[0]

# %%
octopus_noise_df = pd.read_table(octopus_noise_file)
octopus_noise_df

# %%
octopus_noise_df.loc[octopus_noise_df["%Noise"] > 0]

# %%
100 * len(octopus_noise_df.loc[octopus_noise_df["%Noise"] > 0]) / len(octopus_noise_df)

# %%
octopus_noise_df["%Noise"].describe()

# %%
platforms_color_map = {
    platform: color_map
    for platform, color_map in zip(
        platforms, [pacbio_color_discrete_map, illumina_color_discrete_map]
    )
}

# dashes = ["dash", "dot"]
# dashes = ["dash", "solid"]
# platforms_dashes = {platform: dash for platform, dash in zip(platforms, dashes)}

# symbols = ["star", "triangle-up"]
# platforms_symbols = {platform: symbol for platform, symbol in zip(platforms, symbols)}


joined_platforms = [platforms[0]] * len(pacbio_conditions) + [platforms[1]] * len(
    illumina_conditions
)

joined_conditions = pacbio_conditions + illumina_conditions

joined_noise_dfs = pacbio_noise_dfs + illumina_noise_dfs

# %%
len(joined_noise_dfs)

# %%
max_squid_noise = ceil(max([df["%Noise"].max() for df in joined_noise_dfs]))
max_octopus_noise = ceil(octopus_noise_df["%Noise"].max())
max_noise = ic(max([max_squid_noise, max_octopus_noise]))

# %%
column_widths = [len(pacbio_conditions), len(illumina_conditions), len(illumina_conditions)/2]
# subplot_titles = ["Squid's Long-reads", "Squids' Short-reads", "Whole-transcriptome octopus data"]
subplot_titles = ["Squids'<br>Long-reads", "Squids'<br>Short-reads", "Whole-transcriptome<br>octopus data"]

fig = make_subplots(
    rows=1,
    cols=3,
    # y_title="Dispersion [%]",
    # x_title="Gene",
    subplot_titles=subplot_titles,
    # shared_yaxes="all",
    # shared_xaxes="all",
    # vertical_spacing=0.04,
    horizontal_spacing=0.06,
    column_widths=column_widths,
)

for platform, condition, noise_df in zip(
    joined_platforms, joined_conditions, joined_noise_dfs
):
    x = noise_df[condition_col]
    y = noise_df["%Noise"]
    
    color = platforms_color_map[platform][condition]
    col = 1 if platform == "Long-reads" else 2

    fig.add_trace(
        go.Violin(
            x=x,
            y=y,
            # mode="markers",
            # fillcolor=color,
            marker=dict(
                color=color,
                # symbol=platforms_symbols[platform],
                # size=3,
                # opacity=0.3,
            ),
            # mode="lines",
            # line=dict(
            #     color=platforms_color_map[platform][condition],
            #     dash=platforms_dashes[platform],
            #     # size=8,
            # ),
            # opacity=0.3,
            name=condition,
        ),
        row=1,
        col=col,
    )
    
    fig.update_xaxes(row=1, col=col, title_text="Gene", tickangle = 30, tickfont=dict(size=10))
    fig.update_yaxes(row=1, col=col, range=[0, max_noise])
    
    # ic(platform, condition)
    # ic(x)
    # ic(y)
    # break

fig.update_yaxes(row=1, col=1, title_text="Per-site noise level [%]")
    
fig.add_trace(
        go.Histogram(
        y=octopus_noise_df["%Noise"],
        marker_color="black",
        # opacity=0.5
    ),
    row=1,
    col=3,
)
fig.update_xaxes(row=1, col=3, title_text="Genes", type="log")
fig.update_yaxes(row=1, col=3, title_text="Per-gene noise level [%]", range=[0, max_noise])

    
# fig.update_xaxes(
#         tickangle = 30, 
#     # title_standoff = 10,
#     tickfont=dict(size=10)
# )    

# fig.update_yaxes(type="log",
#                  # range=[np.log10(y_min), np.log10(y_max)], nticks=6
#                 )
width = 1200
height = 450

fig.update_annotations(
    # yshift=20, 
    font_size=14
)

fig.update_layout(
    title="Noise levels",
    title_x=0.07,
    title_y=0.95,
    template=template,
    width=width,
    height=height,
    showlegend=False
)

fig.write_image(
    "Per chrom noise levels - PacBio vs. Illumina vs. octopus.svg",
    width=width,
    height=height,
)

fig.show()


# %% jupyter={"source_hidden": true}
# column_widths = [len(conditions) for conditions in [pacbio_conditions, illumina_conditions]]

# fig = make_subplots(
#     rows=1,
#     cols=2,
#     y_title="Per-site noise levels [%]",
#     x_title="Gene",
#     subplot_titles=platforms,
#     shared_yaxes="all",
#     # shared_xaxes="all",
#     vertical_spacing=0.03,
#     # horizontal_spacing=0.015,
#     horizontal_spacing=0.03,
#     column_widths=column_widths,
# )

# for platform, condition, noise_df in zip(
#     joined_platforms, joined_conditions, joined_noise_dfs
# ):
#     x = noise_df[condition_col]
#     y = noise_df["%Noise"]
    
#     color = platforms_color_map[platform][condition]
#     col = 1 if platform == "Long-reads" else 2

#     fig.add_trace(
#         go.Violin(
#             x=x,
#             y=y,
#             # mode="markers",
#             # fillcolor=color,
#             marker=dict(
#                 color=color,
#                 # symbol=platforms_symbols[platform],
#                 # size=3,
#                 # opacity=0.3,
#             ),
#             # mode="lines",
#             # line=dict(
#             #     color=platforms_color_map[platform][condition],
#             #     dash=platforms_dashes[platform],
#             #     # size=8,
#             # ),
#             # opacity=0.3,
#             name=condition,
#         ),
#         row=1,
#         col=col,
#     )
    
    
#     # ic(platform, condition)
#     # ic(x)
#     # ic(y)
#     # break

# fig.update_xaxes(
#         tickangle = 30, 
#     # title_standoff = 10,
#     tickfont=dict(size=10)
# )    

# # fig.update_yaxes(type="log",
# #                  # range=[np.log10(y_min), np.log10(y_max)], nticks=6
# #                 )
# width = 1000
# height = 450

# fig.update_layout(
#     title="Noise levels",
#     title_x=0.07,
#     template=template,
#     width=width,
#     height=height,
#     showlegend=False
# )

# fig.write_image(
#     "Per chrom noise levels - PacBio vs. Illumina.svg",
#     width=width,
#     height=height,
# )

# fig.show()


# %% [markdown]
# # Combined distinct proteins plots for squid

# %%
condition_col = "Gene"
platforms = ["Long-reads", "Short-reads"]

pacbio_conditions = ["GRIA2", "PCLO", "ADAR1", "IQEC1"]

illumina_conditions = [
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
illumina_conditions = [condition.split("_")[0] for condition in illumina_conditions]

pacbio_color_sequence = px.colors.qualitative.G10
pacbio_color_discrete_map = {
    condition: color
    for condition, color in zip(pacbio_conditions, pacbio_color_sequence)
}

illumina_color_sequence = px.colors.qualitative.Dark24
illumina_color_discrete_map = {
    condition: color
    for condition, color in zip(illumina_conditions, illumina_color_sequence)
}

# %%
in_dir = Path("/private7/projects/Combinatorics/Code/Notebooks")

illumina_merged_distinct_proteins_file = Path(in_dir, "DistinctProteins.Illumina.tsv")
pacbio_merged_distinct_proteins_file = Path(in_dir, "DistinctProteins.PacBio.tsv")
pacbio_umi_merged_distinct_proteins_file = Path(
    in_dir, "DistinctProteins.PacBio.UMI.tsv"
)

illumina_merged_max_distinct_proteins_file = Path(
    in_dir, "MaxDistinctProteinsF1.Illumina.tsv"
)
pacbio_merged_max_distinct_proteins_file = Path(
    in_dir, "MaxDistinctProteinsF1.PacBio.tsv"
)
pacbio_umi_merged_max_distinct_proteins_file = Path(
    in_dir, "MaxDistinctProteinsF1.PacBio.UMI.tsv"
)

# %%
platforms_color_map = {
    platform: color_map
    for platform, color_map in zip(
        platforms, [pacbio_color_discrete_map, illumina_color_discrete_map]
    )
}

# %%
platforms_conditions = [pacbio_conditions, illumina_conditions]

# %% [markdown]
# ## Max distinct

# %%
# pacbio_max_distinct_proteins_df = pd.read_table(pacbio_merged_max_distinct_proteins_file)
pacbio_max_distinct_proteins_df = pd.concat(
    [
        pd.read_table(pacbio_merged_max_distinct_proteins_file),
        pd.read_table(pacbio_umi_merged_max_distinct_proteins_file),
    ]
)
illumina_max_distinct_proteins_df = pd.read_table(
    illumina_merged_max_distinct_proteins_file
)

merged_max_distinct_df = pd.concat(
    [pacbio_max_distinct_proteins_df, illumina_max_distinct_proteins_df]
)
merged_max_distinct_df["Color"] = merged_max_distinct_df.apply(
    lambda x: platforms_color_map[x["Platform"]][x[condition_col]], axis=1
)
merged_max_distinct_df

# %%
column_widths = [
    len(conditions) for conditions in [pacbio_conditions, illumina_conditions]
]

font_size = 28

fig = make_subplots(
    rows=1,
    cols=2,
    # y_title="Max distinct protein isoforms observed",
    y_title="Max distinct protein<br>isoforms observed",
    x_title="Gene",
    subplot_titles=platforms,
    shared_yaxes="all",
    # shared_xaxes="all",
    vertical_spacing=0.04,
    # horizontal_spacing=0.015,
    horizontal_spacing=0.04,
    column_widths=column_widths,
)

for col, (platform, platform_conditions) in enumerate(
    zip(platforms, platforms_conditions), start=1
):
    platform_df = merged_max_distinct_df.loc[
        merged_max_distinct_df["Platform"] == platform
    ]
    for condition in platform_conditions:
        condition_df = platform_df.loc[platform_df[condition_col] == condition]
        x = condition_df[condition_col]
        y = condition_df["Max"]
        colors = condition_df["Color"]
        fig.add_trace(
            go.Bar(
                x=x,
                y=y,
                marker=dict(
                    color=colors,
                ),
                name=condition,
            ),
            row=1,
            col=col,
        )


fig.update_xaxes(
    tickangle=30,
    # title_standoff = 10,
    # tickfont=dict(size=10),
    title_font=dict(size=font_size),
    title_standoff=100,
)

fig.update_yaxes(
    title_font=dict(size=font_size),
    tickfont=dict(size=0.7 * font_size),
    title_standoff=100,
)

fig.update_annotations(font_size=font_size)

# width = 1050
width = 1200
height = 550

fig.update_layout(
    # title="Noise levels",
    title_x=0.07,
    template=template,
    width=width,
    height=height,
    showlegend=False,
    margin=dict(l=120, r=20, t=60, b=80),
)

fig.write_image(
    Path(out_dir, "Max distinct proteins per gene - PacBio vs Illumina.svg"),
    width=width,
    height=height,
)

fig.show()


# %% [markdown]
# ## All fractions

# %%
# illumina_merged_distinct_proteins_file = "DistinctProteins.Illumina.tsv"
# pacbio_merged_distinct_proteins_file = "DistinctProteins.PacBio.tsv"

pacbio_distinct_proteins_df = pd.concat(
    [
        pd.read_table(pacbio_merged_distinct_proteins_file),
        pd.read_table(pacbio_umi_merged_distinct_proteins_file),
    ]
)
illumina_distinct_proteins_df = pd.read_table(illumina_merged_distinct_proteins_file)
illumina_distinct_proteins_df[condition_col] = (
    illumina_distinct_proteins_df[condition_col].str.split("_", expand=True).iloc[:, 0]
)
illumina_distinct_proteins_df = illumina_distinct_proteins_df.sort_values(condition_col)

merged_distinct_df = pd.concat(
    [pacbio_distinct_proteins_df, illumina_distinct_proteins_df]
)
# merged_distinct_df["Color"] = merged_distinct_df.apply(
#     lambda x: platforms_color_map[x["Platform"]][x[condition_col]],
#     axis=1
# )
merged_distinct_df

# %%
illumina_distinct_proteins_df

# %%
list(merged_distinct_df.groupby("Platform")["NumOfReads"].max())

# %%
x_axis_name = "Coverage"
y_axis_name = "Distinct protein isoforms observed"
# column_widths = list(merged_distinct_df.groupby("Platform")["NumOfReads"].max())
# column_widths = [300_000, 1_500_000]
column_widths = [500_000, 1_500_000]

_marker_size = 7

font_size = 28

# Initialize figure with subplots
fig = make_subplots(
    rows=1,
    cols=2,
    print_grid=False,
    x_title=x_axis_name,
    y_title=y_axis_name,
    shared_yaxes="all",
    subplot_titles=platforms,
    horizontal_spacing=0.04,
    column_widths=column_widths,
)

# first_data_trace = True

dscam_dashed_lined_width = 3.5
dscam_colors = ["grey", "black"]
dscam_dashes = ["dash", "dash"]
dscam_ys = [
    36_016,
    18_496,
]
dscam_legend_names = [
    "Theoretical maximum    ",
    "Observed",
]


for col, (platform, platform_conditions) in enumerate(
    zip(platforms, platforms_conditions), start=1
):

    # add all conditions under a sequencing platform

    platform_df = merged_distinct_df.loc[merged_distinct_df["Platform"] == platform]

    for condition in platform_conditions:

        condition_df = (
            platform_df.loc[platform_df[condition_col] == condition]
            .sort_values(["Fraction", "NumOfProteins"], ascending=False)
            .drop_duplicates("Fraction", ignore_index=True)
        )

        x_measured = condition_df["NumOfReads"]
        y_measured = condition_df["NumOfProteins"]
        color = platforms_color_map[platform][condition]
        name = condition

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
                legendgroup=platform,  # this can be any string
                legendgrouptitle_text=platform,
                name=name,
            ),
            row=1,
            col=col,
        )

    # add dscam data to compare it with the platform

    platform_maximal_x = platform_df["NumOfReads"].max()

    showlegend = col == 2

    fig.add_trace(
        go.Scatter(
            x=[0.05 * platform_maximal_x, 1.05 * platform_maximal_x],
            y=[dscam_ys[0], dscam_ys[0]],
            mode="lines",
            line=dict(
                color=dscam_colors[0],
                dash=dscam_dashes[0],
                width=dscam_dashed_lined_width,
            ),
            opacity=0.6,
            legendgroup="Drosophila’s DSCAM",  # this can be any string
            legendgrouptitle_text="Drosophila’s DSCAM",
            name=dscam_legend_names[0],
            showlegend=showlegend,
        ),
        row=1,
        col=col,
    )
    fig.add_trace(
        go.Scatter(
            x=[0.05 * platform_maximal_x, 1.05 * platform_maximal_x],
            y=[dscam_ys[1], dscam_ys[1]],
            mode="lines",
            line=dict(
                color=dscam_colors[1],
                dash=dscam_dashes[1],
                width=dscam_dashed_lined_width,
            ),
            opacity=0.6,
            legendgroup="Drosophila’s DSCAM",  # this can be any string
            legendgrouptitle_text="Drosophila’s DSCAM",
            name=dscam_legend_names[1],
            showlegend=showlegend,
        ),
        row=1,
        col=col,
    )


fig.update_xaxes(
    # tickfont=dict(size=0.7 * font_size),
    tickfont=dict(size=0.6 * font_size),
    # tickangle=20,
    tickangle=25,
    # tickangle=10,
    # type="log"
    title_standoff=100,
    automargin=True,
)


fig.update_yaxes(
    # tickfont=dict(size=0.7 * font_size),
    tickfont=dict(size=0.6 * font_size),
    # type="log"
)

fig.update_annotations(font_size=font_size)

# width = 1100
# width = 1300
# height = width*650/900

width = 1250
height = width * 700 / 900

fig.update_layout(
    template=template,
    legend_font=dict(size=0.6 * font_size),
    legend_grouptitlefont=dict(size=0.6 * font_size),
    # legend_font=dict(size=14),
    # legend_grouptitlefont=dict(size=16),
    # legend_tracegroupgap=4,
    width=width,
    height=height,
    margin=dict(l=120, r=40, t=60, b=120),
    # autosize=False,
)
# fig.write_image("Distinct unique proteins vs. sequencing depth - Illumina.png", format='png',engine='kaleido')
fig.write_image(
    Path(out_dir, "Distinct proteins vs. sequencing depth - PacBio vs Illumina.svg"),
    width=width,
    height=height,
)
fig.show()


# %%
import plotly.graph_objects as go
from plotly.subplots import make_subplots

import numpy as np

# Generate synthetic data
x1 = np.linspace(0, 100, 20)
y1 = np.random.randint(10, 100, size=20)
x2 = 
y2 = np.random.randint(10, 100, size=20)

font_size = 28

fig = make_subplots(
    rows=1,
    cols=2,
    x_title="X-axis title",
    y_title="Y-axis title",
    subplot_titles=["Group A", "Group B"],
    shared_yaxes=True,
    horizontal_spacing=0.04,
)

fig.add_trace(
    go.Scatter(
        x=x,
        y=y1,
        mode="lines+markers",
        name="A",
        marker=dict(size=8, color="blue"),
    ),
    row=1,
    col=1,
)
fig.add_trace(
    go.Scatter(
        x=x,
        y=y2,
        mode="lines+markers",
        name="B",
        marker=dict(size=8, color="red"),
    ),
    row=1,
    col=2,
)

fig.update_xaxes(
    tickfont=dict(size=0.7 * font_size),
    tickangle=30,
    title_standoff=100,  # <-- This has no effect!
    automargin=True,
)
fig.update_xaxes(title_standoff=100, automargin=True, row=1, col=1)
fig.update_xaxes(title_standoff=100, automargin=True, row=1, col=2)

fig.update_yaxes(
    tickfont=dict(size=0.7 * font_size),
)

fig.update_annotations(font_size=font_size)

width = 900
height = 500

fig.update_layout(
    template="plotly_white",
    width=width,
    height=height,
    margin=dict(l=120, r=40, t=60, b=120),
)

fig.show()

# %%
height

# %% [markdown]
# # Combined editing in motifs plots for squid long-reads

# %%
condition_col = "Gene"
platforms = ["Long-reads", "Short-reads"]

pacbio_conditions = ["GRIA2", "PCLO", "ADAR1", "IQEC1"]

# illumina_conditions = [
#     "RUSC2_MOUSE",
#     "TRIM2_BOVIN",
#     "CA2D3_MOUSE",
#     "ABL_DROME",
#     "DGLA_HUMAN",
#     "K0513_MOUSE",
#     "KCNAS_DROME",
#     "ACHA4_MOUSE",
#     "ANR17_HUMAN",
#     "TWK7_CAEEL",
#     "SCN1_HETBL",
#     "CACB2_RABIT",
#     "RIMS2_RAT",
#     "PCLO_CHICK",
#     "DOP1_HUMAN",
#     "IQEC1_HUMAN",
#     "CSKI1_MOUSE",
#     "MTUS2_HUMAN",
#     "ROBO2_HUMAN",
# ]
# illumina_conditions = [condition.split("_")[0] for condition in illumina_conditions]

pacbio_color_sequence = px.colors.qualitative.G10
pacbio_color_discrete_map = {
    condition: color
    for condition, color in zip(pacbio_conditions, pacbio_color_sequence)
}

# %%
in_dir = Path("/private7/projects/Combinatorics/Code/Notebooks")

# illumina_merged_distinct_proteins_file = Path(in_dir, "DistinctProteins.Illumina.tsv")
pacbio_editing_in_motifs_file = Path(in_dir, "EditingInMotifs.PacBio.tsv")
pacbio_umi_editing_in_motifs_file = Path(
    in_dir,
    "EditingInMotifs.PacBio.UMI.tsv",
)

# illumina_merged_max_distinct_proteins_file = Path(
#     in_dir, "MaxDistinctProteinsF1.Illumina.tsv"
# )

# %%
platforms_color_map = {
    platform: color_map
    for platform, color_map in zip(
        platforms,
        [
            pacbio_color_discrete_map,
            # illumina_color_discrete_map
        ],
    )
}
platforms_color_map

# %%
platforms_conditions = [
    pacbio_conditions,
    # illumina_conditions
]

# %%
# pacbio_max_distinct_proteins_df = pd.read_table(pacbio_merged_max_distinct_proteins_file)
pacbio_concat_signature_positions_df = pd.concat(
    [
        pd.read_table(
            pacbio_editing_in_motifs_file,
            dtype={"Edited": bool, "PositionInSignature": str},
            na_filter=False,
        ),
        pd.read_table(
            pacbio_umi_editing_in_motifs_file,
            dtype={"Edited": bool, "PositionInSignature": str},
            na_filter=False,
        ),
    ]
)
# illumina_max_distinct_proteins_df = pd.read_table(
#     illumina_merged_max_distinct_proteins_file
# )

# merged_max_distinct_df = pd.concat(
#     [pacbio_max_distinct_proteins_df, illumina_max_distinct_proteins_df]
# )
# merged_max_distinct_df["Color"] = merged_max_distinct_df.apply(
#     lambda x: platforms_color_map[x["Platform"]][x[condition_col]], axis=1
# )
# merged_max_distinct_df

pacbio_concat_signature_positions_df["Color"] = (
    pacbio_concat_signature_positions_df.apply(
        lambda x: platforms_color_map[x["Platform"]][x[condition_col]],
        axis=1,
    )
)

pacbio_concat_signature_positions_df

# %%
pacbio_concat_signature_positions_df["PositionInSignature"].value_counts(dropna=False)

# %%
pacbio_concat_signature_positions_proportion_df = (
    pacbio_concat_signature_positions_df.groupby([condition_col, "Edited"])[
        "PositionInSignature"
    ]
    .value_counts(normalize=True)
    .mul(100)
    .reset_index()
    .rename(columns={"proportion": "%OfSites"})
)

pacbio_concat_signature_positions_proportion_df[
    "Edited"
] = pacbio_concat_signature_positions_proportion_df["Edited"].apply(
    # lambda x: "Edited" if x else "Not edited"
    lambda x: "Edited" if x else "Not<br>edited"
)

pacbio_concat_signature_positions_proportion_df

# %%
df = pacbio_concat_signature_positions_proportion_df.loc[
    pacbio_concat_signature_positions_proportion_df["PositionInSignature"].eq("No")
]

df

# %%
df[condition_col].apply(lambda x: pacbio_color_discrete_map[x])

# %%
pacbio_color_discrete_map

# %%
x = [df[condition_col].values, df["Edited"].values]
y = df["%OfSites"].values

x

# %%
y


# %%
def signature_editing_test(signature_positions_df):
    edited_in_signature = signature_positions_df.loc[
        signature_positions_df["Edited"]
        & signature_positions_df["PositionInSignature"].eq("Yes")
    ].shape[0]
    not_edited_in_signature = signature_positions_df.loc[
        ~signature_positions_df["Edited"]
        & signature_positions_df["PositionInSignature"].eq("Yes")
    ].shape[0]
    edited_outside_signature = signature_positions_df.loc[
        signature_positions_df["Edited"]
        & signature_positions_df["PositionInSignature"].eq("No")
    ].shape[0]
    not_edited_outside_signature = signature_positions_df.loc[
        ~signature_positions_df["Edited"]
        & signature_positions_df["PositionInSignature"].eq("No")
    ].shape[0]

    if edited_in_signature + not_edited_in_signature == 0:
        return np.nan

    table = [
        [edited_in_signature, edited_outside_signature],
        [not_edited_in_signature, not_edited_outside_signature],
    ]

    oddsratio, pvalue = scipy.stats.fisher_exact(
        table,
        # alternative="greater"
    )  # "greater" tests for enrichment

    edited_to_non_edited_in_signature_ratio = (
        edited_in_signature / not_edited_in_signature
    )
    edited_to_non_edited_outside_signature_ratio = (
        edited_outside_signature / not_edited_outside_signature
    )

    # return (
    #     oddsratio,
    #     pvalue,
    #     edited_to_non_edited_in_signature_ratio,
    #     edited_to_non_edited_outside_signature_ratio,
    # )

    return pvalue


# %%
# Prepare to collect p-values for annotation
fisher_pvals = {}

# We'll compare "Yes" vs "No" for each condition (Gene)
for condition in pacbio_concat_signature_positions_df[condition_col].unique():
    df = pacbio_concat_signature_positions_df[
        pacbio_concat_signature_positions_df[condition_col] == condition
    ]
    pvalue = signature_editing_test(df)
    fisher_pvals[condition] = pvalue

fisher_pvals

# %%
fig = go.Figure()

font_size = 14

for edited_in_siganture, marker_pattern_shape in zip(
    ["No", "Yes", "NA"], ["", "/", "."]
):
    df = pacbio_concat_signature_positions_proportion_df.loc[
        pacbio_concat_signature_positions_proportion_df["PositionInSignature"]
        == edited_in_siganture
    ]
    x = [df[condition_col].values, df["Edited"].values]
    y = df["%OfSites"].values

    colors = df[condition_col].apply(lambda x: pacbio_color_discrete_map[x])

    fig.add_trace(
        go.Bar(
            x=x,
            y=y,
            name=edited_in_siganture,
            marker_pattern_shape=marker_pattern_shape,
            marker_color=colors,
            showlegend=False,
        )
    )

# Add a black and white legend for pattern shapes: ["No", "Yes", "NA"] and patterns ["", "/", "."]
legend_labels = ["No", "Yes", "NA"]
legend_patterns = ["", "/", "."]
legend_colors = ["black"] * 3  # All black

for label, pattern in zip(legend_labels, legend_patterns):
    fig.add_trace(
        go.Bar(
            x=[None],
            y=[None],
            marker=dict(
                color="black",
                line_color="black",
                pattern_fillmode="replace",
                pattern_shape=pattern,
            ),
            name=label,
            showlegend=True,
        )
    )

# --- In your plotting code, after the bars are added, annotate p-values ---

for condition, pval in fisher_pvals.items():
    if np.isnan(pval):
        continue
    fig.add_annotation(
        x=(condition, "Edited"),  # <-- tuple for grouped bar x-axis
        xshift=-35,
        # x=[
        #     [condition] * 2,
        #     ["Not<br>edited", "Edited"],
        # ],  # <-- tuple for grouped bar x-axis
        y=101,  # adjust as needed
        # text=f"p = {pval:.2e}" if not np.isnan(pval) else "n/a",
        text=f"p-val = {pval:.3g}",
        showarrow=False,
        font=dict(color="black", size=font_size),
        yanchor="bottom",
    )

fig.update_xaxes(tickangle=30)

fig.update_yaxes(
    title="% of coding adenosines",
    # title_font=dict(size=font_size),
    #     tickfont=dict(size=0.7 * font_size)
)

# fig.update_yaxes(
#     title_font=dict(size=font_size),
#     tickfont=dict(size=0.7 * font_size),
#     title_standoff=100,
# )

# fig.update_annotations(font_size=font_size)

width = 800
height = 450

fig.update_layout(
    barmode="relative",
    width=width,
    height=height,
    template=template,
    legend=dict(
        title_text="In signature   ",
    ),
    # barcornerradius=15,
)

fig.write_image(
    Path(out_dir, "Editing in domains - PacBio.svg"),
    width=width,
    height=height,
)


fig.show()

# %%
fig = px.bar(
    pacbio_concat_signature_positions_proportion_df,
    x="Edited",
    y="%OfSites",
    color=condition_col,
    color_discrete_map=pacbio_color_discrete_map,
    pattern_shape="PositionInSignature",
    #  pattern_shape_map={
    #  "China": ".", "Canada": "/", "South Korea": "+"
    #  },
    pattern_shape_sequence=["", "/", "."],
    facet_col=condition_col,
)

width = 700
height = 450

fig.update_layout(
    width=width,
    height=height,
    template=template,
)


fig.show()

# %% [markdown]
# # Combined dispersion plots

# %%
condition_col = "Gene"
platforms = ["Long-reads", "Short-reads"]

pacbio_conditions = ["GRIA2", "PCLO"]

illumina_conditions = [
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
illumina_conditions = [condition.split("_")[0] for condition in illumina_conditions]

pacbio_color_sequence = px.colors.qualitative.G10
pacbio_color_discrete_map = {
    condition: color
    for condition, color in zip(pacbio_conditions, pacbio_color_sequence)
}

illumina_color_sequence = px.colors.qualitative.Dark24
illumina_color_discrete_map = {
    condition: color
    for condition, color in zip(illumina_conditions, illumina_color_sequence)
}

# %%
platforms_color_map = {
    platform: color_map
    for platform, color_map in zip(
        platforms, [pacbio_color_discrete_map, illumina_color_discrete_map]
    )
}

# %%
platforms_color_map

# %%
pacbio_dispersion_file = "Dispersion.PacBio.tsv"
illumina_dispersion_file = "Dispersion.Illumina.tsv"
octopus_dispersion_file = "Dispersion.Octopus.tsv"

# %%
pacbio_dispersion_df = pd.read_table(pacbio_dispersion_file)
pacbio_dispersion_df[condition_col] = pacbio_dispersion_df[condition_col].apply(
    lambda x: "GRIA2" if x == "GRIA" else x
)

illumina_dispersion_df = pd.read_table(illumina_dispersion_file)
# illumina_dispersion_df

octopus_dispersion_df = pd.read_table(octopus_dispersion_file)
octopus_dispersion_df

# %%
octopus_dispersion_df["%SolutionsDispersion"].describe()

# %%
octopus_dispersion_df.loc[octopus_dispersion_df["%SolutionsDispersion"] > 0].sort_values("%SolutionsDispersion")

# %%
# described_octopus_dispersion_df = octopus_dispersion_df["%SolutionsDispersion"].describe()
# quartiles = ["25%", "50%", "75%"]
# # dispersion_quartiles = described_octopus_dispersion_df.loc[quartiles].mul(100).values
# dispersion_quartiles = described_octopus_dispersion_df.loc[quartiles].values
# dispersion_quartiles

# %%
percent_of_octopus_genes_with_nonzero_dispersion = 100 * len(octopus_dispersion_df.loc[octopus_dispersion_df["%SolutionsDispersion"] > 0]) / len(octopus_dispersion_df)
percent_of_octopus_genes_with_nonzero_dispersion

# %%
column_widths = [len(pacbio_conditions), len(illumina_conditions), len(illumina_conditions)/2]
# subplot_titles = ["Squid's Long-reads", "Squids' Short-reads", "Whole-transcriptome octopus data"]
subplot_titles = ["Squids'<br>Long-reads", "Squids'<br>Short-reads", "Whole-transcriptome<br>octopus data"]

fig = make_subplots(
    rows=1,
    cols=3,
    y_title="Dispersion [%]",
    # x_title="Gene",
    subplot_titles=subplot_titles,
    # shared_yaxes="all",
    # shared_xaxes="all",
    # vertical_spacing=0.04,
    horizontal_spacing=0.06,
    column_widths=column_widths,
)
                 
dispersion_dfs = [pacbio_dispersion_df, illumina_dispersion_df, octopus_dispersion_df]
x_titles = ["Gene", "Gene", "Genes"]
platforms_conditions = [pacbio_conditions, illumina_conditions, None]
platforms = ["Long-reads", "Short-reads", None]

max_squid_dispersion = ceil(max([pacbio_dispersion_df["%SolutionsDispersion"].max(), illumina_dispersion_df["%SolutionsDispersion"].max()]))

for col, (dispersion_df, x_title, platform, platform_conditions) in enumerate(zip(dispersion_dfs, x_titles, platforms, platforms_conditions), start=1):
    
    if col != 3:
    
        for condition in platform_conditions:
            condition_df = dispersion_df.loc[dispersion_df[condition_col] == condition]
            x = condition_df[condition_col]
            y = condition_df["%SolutionsDispersion"]
            # colors = condition_df["Color"]
            # ic(platform, platforms_color_map[platform], condition)
            
            color = platforms_color_map[platform][condition]
            fig.add_trace(
                go.Bar(
                    x=x,
                    y=y,
                    marker=dict(
                        # color=colors,
                        color=color
                    ),
                    name=condition,
                ),
                row=1,
                col=col,
            )
        fig.update_xaxes(row=1, col=col, title_text=x_title, tickangle = 30, tickfont=dict(size=10))
        fig.update_yaxes(row=1, col=col, range=[0, max_squid_dispersion])
    
    else:
        fig.add_trace(
            go.Histogram(
                y=dispersion_df["%SolutionsDispersion"],
                marker_color="black",
                # opacity=0.5
            ),
            row=1,
            col=col,
        )
        fig.update_xaxes(row=1, col=col, title_text=x_title, type="log", )

# f = fig.full_figure_for_development(warn=False)
# f_data = f.data[-1]
# y = f_data["y"]
# ybins = f_data.ybins
# plotbins = list(
#     np.arange(
#         start=ybins["start"],
#         stop=ybins["end"] + ybins["size"],
#         step=ybins["size"],
#     )
# )
# counts, bins = np.histogram(list(y), bins=plotbins)
# # max_count = max(counts)
# # percent_of_octopus_genes_with_nonzero_dispersion_x = max_count * 0.6 # 1762 - too high
# counts.flatten().sort()
# # percent_of_octopus_genes_with_nonzero_dispersion_x = counts[:2].mean() # 1474 - too high
# max_count_1, max_count_2 = counts[:2]
# percent_of_octopus_genes_with_nonzero_dispersion_x = max_count_2 * 3.5
# # percent_of_octopus_genes_with_nonzero_dispersion_x = np.log10(counts[:2].mean())
# # percent_of_octopus_genes_with_nonzero_dispersion_x = np.log10(max_count_1)
# ic(percent_of_octopus_genes_with_nonzero_dispersion_x)

# text = f"%{100-percent_of_octopus_genes_with_nonzero_dispersion:.1f} of<br>genes<br>"
# # text = f"                         %{100-percent_of_octopus_genes_with_nonzero_dispersion:.1f} of genes"
# fig.add_shape(
#     type="line",
#     x0=0,
#     x1=percent_of_octopus_genes_with_nonzero_dispersion_x,
#     y0=0.5,
#     y1=0.5,
#     line=dict(
#         color="red",
#         width=4,
#         dash="dash",
#     ),
#     opacity=0.5,
#     label=dict(
#         text=text,
#         textposition="end",
#         textangle=12,
#         # font_size=14,
#     ),
#     row=1,
#     col=3
# )
        
fig.update_annotations(yshift=20, font_size=14)

width = 1200
height = 450

fig.update_layout(
    title="Dispersion of estimates for the number of distinct protein isoforms<br><br>",
    title_x=0.05,
    template=template,
    width=width,
    height=height,
    showlegend=False,
    title_y=0.95,
    # automargin=True, 
    # yref='paper'
)

fig.write_image(
    "SolutionsDispersion - PacBio vs Illumina vs octopus.svg",
    width=width,
    height=height,
)

fig.show()


# %% [markdown]
# # Combined correlation and MI plots

# %%
corr_colorscales = [
    "RdBu",
    [
        [0.00, "#67001f"],
        [0.20, "#b2182b"],
        [0.35, "#ef8a62"],
        [0.48, "#f7f7f7"],  # narrow neutral zone
        [0.52, "#f7f7f7"],
        [0.65, "#67a9cf"],
        [0.80, "#2166ac"],
        [1.00, "#053061"],
    ],
    [
        [0.00, "#3b0f70"],
        [0.20, "#8c2981"],
        [0.40, "#de4968"],
        [0.50, "#f7f7f7"],
        [0.60, "#9bd3c6"],
        [0.80, "#4aa6c3"],
        [1.00, "#1f3c88"],
    ]
]
corrs_final_colorscale = corr_colorscales[0]

# %%
mi_colorscales = [
    "Blues",
    # # sns_colorscale
    # [
    #     [0.0, '#3f7f93'], #cmap = sns.diverging_palette(220, 10, as_cmap = True)
    #     [0.071, '#5890a1'],
    #     [0.143, '#72a1b0'],
    #     [0.214, '#8cb3bf'],
    #     [0.286, '#a7c5cf'],
    #     [0.357, '#c0d6dd'],
    #     [0.429, '#dae8ec'],
    #     [0.5, '#f2f2f2'],
    #     [0.571, '#f7d7d9'],
    #     [0.643, '#f2bcc0'],
    #     [0.714, '#eda3a9'],
    #     [0.786, '#e8888f'],
    #     [0.857, '#e36e76'],
    #     [0.929, '#de535e'],
    #     [1.0, '#d93a46']
    # ],
    # upper part of corrs_final_colorscale
    [
        [0.00, "#f7f7f7"],  # MI ≈ 0 → neutral / background
        [0.25, "#d1e5f0"],
        [0.50, "#67a9cf"],
        [0.75, "#2166ac"],
        [1.00, "#053061"],  # max MI
    ]
]
mi_final_colorscale = mi_colorscales[-1]

# %% [markdown]
# ## PCLO - long- vs short-reads

# %%
platforms = ["Long-reads", "Short-reads"]

# %%
pacbio_corr_file = Path(out_dir, "PCLOMaskedBonferroniCorr.PacBio.tsv")
illumina_corr_file = Path(out_dir, "PCLOMaskedBonferroniCorr.Illumina.tsv")

pacbio_mi_file = Path(out_dir, "PCLOMaskedMI.PacBio.tsv")
illumina_mi_file = Path(out_dir, "PCLOMaskedMI.Illumina.tsv")

# %%
pacbio_corr = pd.read_table(pacbio_corr_file).values
illumina_corr = pd.read_table(illumina_corr_file).values

pacbio_mi = pd.read_table(pacbio_mi_file).values
illumina_mi = pd.read_table(illumina_mi_file).values

# %%
pacbio_corr

# %%
pclo_corrs_cmin = min(np.nanmin(pacbio_corr), np.nanmin(illumina_corr))
pclo_corrs_cmax = max(np.nanmax(pacbio_corr), np.nanmax(illumina_corr))
pclo_corrs_abs_max = np.round(np.max([abs(pclo_corrs_cmin), abs(pclo_corrs_cmax)]), 1)
ic(pclo_corrs_cmin, pclo_corrs_cmax, pclo_corrs_abs_max);

# %%
fig = make_subplots(
    rows=1,
    cols=2,
    subplot_titles=platforms,
    horizontal_spacing=0.06,
)


fig.add_trace(
    go.Heatmap(
        z=pacbio_corr,
        xgap=1, 
        ygap=1,
        coloraxis="coloraxis",
    ),
    row=1,
    col=1
)
fig.add_trace(
    go.Heatmap(
        z=illumina_corr,
        xgap=1, 
        ygap=1,
        coloraxis="coloraxis",
    ),
    row=1,
    col=2
)

fig.update_xaxes(
    showticklabels=False,
    showgrid=False,
    zeroline=False
)
fig.update_yaxes(
    showgrid=False,
    zeroline=False,
    autorange="reversed",
    showticklabels=False,
)

width = 1100
height = 600

fig.update_layout(
    width=width,
    height=height,
    template=template, 
    coloraxis=dict(
        colorscale=corrs_final_colorscale,
        colorbar_thickness=20,
        colorbar_ticklen=3,
        cmid=0,
        # cmax=pclo_corrs_abs_max,
        # cmin=-pclo_corrs_abs_max,
        cmax=1,
        cmin=-1,
        colorbar=dict(
            dtick=0.25,
            len=0.6,
        ),
    ),
    title=dict(
        x=0.1,
        text="Site-Site correlations (PCLO)",
  ),
)

fig.write_image(
    Path(out_dir, "Pearson's r between editing sites in squid's PCLO - pacbio vs illumina.svg"),
    width=width,
    height=height,
)

fig.show()

# %%
pclo_mis_cmin = min(np.nanmin(pacbio_mi), np.nanmin(illumina_mi))
pclo_mis_cmax = max(np.nanmax(pacbio_mi), np.nanmax(illumina_mi))

# %%
fig = make_subplots(
    rows=1,
    cols=2,
    subplot_titles=platforms,
    horizontal_spacing=0.06,
)


fig.add_trace(
    go.Heatmap(
        z=pacbio_mi,
        xgap=1, 
        ygap=1,
        coloraxis="coloraxis",
    ),
    row=1,
    col=1
)

fig.add_trace(
    go.Heatmap(
        z=illumina_mi,
        xgap=1, 
        ygap=1,
        coloraxis="coloraxis",
    ),
    row=1,
    col=2
)

fig.update_xaxes(
    showticklabels=False,
    showgrid=False,
    zeroline=False
)

fig.update_yaxes(
    showgrid=False,
    zeroline=False,
    autorange="reversed",
    showticklabels=False,
    # tick0=1,
    # range=[1, 100],
)

width = 1100
height = 600


fig.update_layout(
    width=width,
    height=height,
    template=template, 
    coloraxis=dict(
        colorscale=mi_final_colorscale,
        # cmin=pclo_mis_cmin,
        # cmax=pclo_mis_cmax,
        cmin=0,
        cmax=1,
        colorbar_thickness=20,
        colorbar_ticklen=3,
        colorbar=dict(
            # dtick=0.25,
            dtick=0.2,
            len=0.6,
        ),
    ),
    title=dict(
        x=0.1,
        text="Mutual information (PCLO)",
        # font=dict(size=18)
    ),
    # font_size=16,
)

fig.write_image(
    Path(out_dir, "Normalized mutual information between editing sites in squid's PCLO - pacbio vs illumina.svg"),
    width=width,
    height=height,
)

fig.show()

# %% [markdown]
# ## Other 3 long-reads genes

# %%
platforms = ["Long-reads", "Short-reads"]

# %%
conditions = ["GRIA2","ADAR1", "IQEC1"]

# %%
corr_files = [
    Path(out_dir, f"{condition}MaskedBonferroniCorr.PacBio.tsv")
    for condition in conditions
]
corr_files

# %%
mi_files = [
    Path(out_dir, f"{condition}MaskedMI.PacBio.tsv")
    for condition in conditions
]
mi_files

# %%
corrs = [
    pd.read_table(corr_file).values
    for corr_file in corr_files
]

mis = [
    pd.read_table(mi_file).values
    for mi_file in mi_files
]

# %%
not_na_min_corrs = [
    np.nanmin(corr)
    for corr in corrs
]
not_na_max_corrs = [
    np.nanmax(corr)
    for corr in corrs
]
not_na_min_corrs + not_na_max_corrs
max_abs_cor_val = np.max(np.abs(not_na_min_corrs + not_na_max_corrs))
max_abs_cor_val = np.round(max_abs_cor_val, 1)
max_abs_cor_val

# %%
fig = make_subplots(
    rows=1,
    cols=3,
    subplot_titles=conditions,
    horizontal_spacing=0.06,
)

for col, corr in enumerate(corrs, start=1):
    fig.add_trace(
        go.Heatmap(
            z=corr,
            xgap=1, 
            ygap=1,
            coloraxis="coloraxis",
        ),
        row=1,
        col=col
    )



fig.update_xaxes(
    showticklabels=False,
    showgrid=False,
    zeroline=False
)

fig.update_yaxes(
    showgrid=False,
    zeroline=False,
    autorange="reversed",
    showticklabels=False,
)

width = 1500
height = 600


fig.update_layout(
    width=width,
    height=height,
    template=template, 
    coloraxis=dict(
        colorscale=corrs_final_colorscale,
        colorbar_thickness=20,
        colorbar_ticklen=3,
        cmid=0,
        cmin=-1,
        cmax=1,
        # cmax=max_abs_cor_val,
        # cmin=-max_abs_cor_val,
        colorbar=dict(
            dtick=0.25,
            len=0.6,
        ),
    ),
    title=dict(
        x=0.1,
        # text="Site-Site correlations (PCLO)",
        text="Site-Site correlations (Long-reads)",
    ),
)

fig.write_image(
    Path(out_dir, "Pearson's r between editing sites in squid except PCLO - pacbio.svg"),
    width=width,
    height=height,
)

fig.show()

# %%
fig = make_subplots(
    rows=1,
    cols=3,
    subplot_titles=conditions,
    horizontal_spacing=0.06,
)

for col, mi in enumerate(mis, start=1):
    fig.add_trace(
        go.Heatmap(
            z=mi,
            xgap=1, 
            ygap=1,
            coloraxis="coloraxis",
        ),
        row=1,
        col=col
    )



fig.update_xaxes(
    showticklabels=False,
    showgrid=False,
    zeroline=False
)

fig.update_yaxes(
    showgrid=False,
    zeroline=False,
    autorange="reversed",
    showticklabels=False,
)

width = 1500
height = 600


fig.update_layout(
    width=width,
    height=height,
    template=template, 
    coloraxis=dict(
        colorscale=mi_final_colorscale,
        cmin=0,
        cmax=1,
        colorbar_thickness=20,
        colorbar_ticklen=3,
        colorbar=dict(
            # dtick=0.25,
            dtick=0.2,
            len=0.6,
        ),
    ),
    title=dict(
        x=0.1,
        # text="Site-Site correlations (PCLO)",
        text="Mutual information (Long-reads)",
    ),
)

fig.write_image(
    Path(out_dir, "Normalized mutual information between editing sites in squid's except PCLO - pacbio .svg"),
    width=width,
    height=height,
)

fig.show()

# %% [markdown]
# # Combined entropy plots

# %%
condition_col = "Gene"
platforms = ["Long-reads", "Short-reads"]

pacbio_conditions = ["GRIA2", "PCLO"]

illumina_conditions = [
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
illumina_conditions = [condition.split("_")[0] for condition in illumina_conditions]

platforms_conditions = [pacbio_conditions, illumina_conditions]

pacbio_color_sequence = px.colors.qualitative.G10
pacbio_color_discrete_map = {
    condition: color
    for condition, color in zip(pacbio_conditions, pacbio_color_sequence)
}

illumina_color_sequence = px.colors.qualitative.Dark24
illumina_color_discrete_map = {
    condition: color
    for condition, color in zip(illumina_conditions, illumina_color_sequence)
}

platforms_color_map = {
    platform: color_map
    for platform, color_map in zip(
        platforms, [pacbio_color_discrete_map, illumina_color_discrete_map]
    )
}

# %%
pacbio_file = "ShanonEntropy.PacBio.tsv"
illumina_file = "ShanonEntropy.Illumina.tsv"

# %%
illumina_df = pd.read_table(illumina_file)
illumina_df

pacbio_df = pd.read_table(pacbio_file)
pacbio_df[condition_col] = pacbio_df[condition_col].apply(lambda x: x if x != "GRIA" else "GRIA2")
# pacbio_df.insert(0, "Platform")
pacbio_df

# %%
column_widths = [len(conditions) for conditions in [pacbio_conditions, illumina_conditions]]

fig = make_subplots(
    rows=1,
    cols=2,
    y_title="Entropy",
    x_title="Gene",
    subplot_titles=platforms,
    shared_yaxes="all",
    # shared_xaxes="all",
    vertical_spacing=0.03,
    # horizontal_spacing=0.015,
    horizontal_spacing=0.03,
    column_widths=column_widths,
)

shannon_dfs = [pacbio_df, illumina_df]

for col, (shannon_df, platform, conditions) in enumerate(zip(shannon_dfs, platforms, platforms_conditions), start=1):
    
    colors = [platforms_color_map[platform][condition] for condition in conditions]

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
        ),
        row=1,
        col=col
    )

    fig.add_trace(
        go.Bar(
            x=conditions,
            y=shannon_df.loc[shannon_df["EntropyName"] == "Hypothetical", "EntropyValue"],
            marker_color=colors,
            marker_pattern_shape="",
            # text=non_syns_per_max_sol_exp_df,
            # textposition="outside",
            # name="Hypothetical",
            showlegend=False,
            # width=0.3
            textfont=dict(size=20),
        ),
        row=1,
        col=col
    )
    
# Add single entropy traces for legend
fig.add_trace(
    go.Bar(
        x=[None],
        y=[None],
        # visible="legendonly",
        marker=dict(
            # color="white",
            color="black",
            # size=8,
            pattern_shape="/",
            line=dict(
                # color="grey",
                color="black",
                # width=8,
            ),
        ),
        legendgroup="Observed",
        name="Observed",
        showlegend=True,
        visible="legendonly",
    )
)
fig.add_trace(
    go.Bar(
        x=[None],
        y=[None],
        # visible="legendonly",
        marker=dict(
            # color="white",
            color="black",
            # size=8,
            pattern_shape="",
            line=dict(
                # color="grey",
                color="black",
                # width=4,
            ),
        ),
        legendgroup="Hypothetical",
        name="Hypothetical",
        showlegend=True,
        visible="legendonly",
    )
)


# fig.update_yaxes(title_text="Entropy")
fig.update_xaxes(
    # title="Gene",
    tickfont=dict(size=10)
)  # https://plotly.com/python/axes/#set-axis-label-rotation-and-font

# fig.update_traces(
#     marker=dict(line_color="black", line_width=0.3, pattern_fillmode="replace"),
#     # width=0.3
# )

width = 1100
height = 500

fig.update_layout(
    # title=f"Shannon's entropy of a largest solution of each {condition_col.lower()}",
    # title="Squid's Long-reads",
    # title_x=0.2,
    width=width,
    height=height,
    # showlegend=False,
    template=template,
    barmode="group",
    bargap=0.2,  # gap between bars of adjacent location coordinates.
    bargroupgap=0.15,  # gap between bars of the same location coordinate.
    legend=dict(
        yanchor="bottom",
        # y=0.91,
        y=0.8,
        xanchor="right",
        x=0.99,
        orientation="h",
        # itemwidth=40,
        
    ),
)

fig.write_image(
    "Observed vs hypothetical entropy - PacBio vs Illumina.svg",
    width=width,
    height=height,
)

fig.show()

# %% [markdown]
# # Combined kmeans silhouette scores for squid long-reads

# %%
condition_col = "Gene"
platforms = ["Long-reads", "Short-reads"]

pacbio_conditions = ["GRIA2", "PCLO", "ADAR1", "IQEC1"]

# illumina_conditions = [
#     "RUSC2_MOUSE",
#     "TRIM2_BOVIN",
#     "CA2D3_MOUSE",
#     "ABL_DROME",
#     "DGLA_HUMAN",
#     "K0513_MOUSE",
#     "KCNAS_DROME",
#     "ACHA4_MOUSE",
#     "ANR17_HUMAN",
#     "TWK7_CAEEL",
#     "SCN1_HETBL",
#     "CACB2_RABIT",
#     "RIMS2_RAT",
#     "PCLO_CHICK",
#     "DOP1_HUMAN",
#     "IQEC1_HUMAN",
#     "CSKI1_MOUSE",
#     "MTUS2_HUMAN",
#     "ROBO2_HUMAN",
# ]
# illumina_conditions = [condition.split("_")[0] for condition in illumina_conditions]

pacbio_color_sequence = px.colors.qualitative.G10
pacbio_color_discrete_map = {
    condition: color
    for condition, color in zip(pacbio_conditions, pacbio_color_sequence)
}

# %%
in_dir = Path("/private7/projects/Combinatorics/Code/Notebooks")

# illumina_merged_distinct_proteins_file = Path(in_dir, "DistinctProteins.Illumina.tsv")
pacbio_kmeans_silhouette_scores_file = Path(in_dir, "KMeansMeanSilhouetteScores.PacBio.tsv")
pacbio_umi_kmeans_silhouette_scores_file = Path(
    in_dir,
    "KMeansMeanSilhouetteScores.PacBio.WithUMIs.tsv",
)

# illumina_merged_max_distinct_proteins_file = Path(
#     in_dir, "MaxDistinctProteinsF1.Illumina.tsv"
# )

# %%
platforms_color_map = {
    platform: color_map
    for platform, color_map in zip(
        platforms,
        [
            pacbio_color_discrete_map,
            # illumina_color_discrete_map
        ],
    )
}
platforms_color_map

# %%
platforms_conditions = [
    pacbio_conditions,
    # illumina_conditions
]

# %%
# pacbio_max_distinct_proteins_df = pd.read_table(pacbio_merged_max_distinct_proteins_file)
pacbio_concat_kmeans_silhouette_scores_df = pd.concat(
    [
        pd.read_table(
            pacbio_kmeans_silhouette_scores_file,
            # dtype={"Edited": bool, "PositionInSignature": str},
            # na_filter=False,
        ),
        pd.read_table(
            pacbio_umi_kmeans_silhouette_scores_file,
            # dtype={"Edited": bool, "PositionInSignature": str},
            # na_filter=False,
        ),
    ]
)
# illumina_max_distinct_proteins_df = pd.read_table(
#     illumina_merged_max_distinct_proteins_file
# )

# merged_max_distinct_df = pd.concat(
#     [pacbio_max_distinct_proteins_df, illumina_max_distinct_proteins_df]
# )
# merged_max_distinct_df["Color"] = merged_max_distinct_df.apply(
#     lambda x: platforms_color_map[x["Platform"]][x[condition_col]], axis=1
# )
# merged_max_distinct_df

pacbio_concat_kmeans_silhouette_scores_df[condition_col] = pacbio_concat_kmeans_silhouette_scores_df[condition_col].apply(
    lambda x: x if x != "GRIA" else "GRIA2"
)

pacbio_concat_kmeans_silhouette_scores_df["Color"] = (
    pacbio_concat_kmeans_silhouette_scores_df.apply(
        lambda x: platforms_color_map[x["Platform"]][x[condition_col]],
        axis=1,
    )
)

pacbio_concat_kmeans_silhouette_scores_df

# %%
title_font_size = 24

fig = px.scatter(
    pacbio_concat_kmeans_silhouette_scores_df,
    x="ClusterSize",
    y="MeanSilhouetteScore",
    color=condition_col,
    color_discrete_map=platforms_color_map[platforms[0]],
    opacity=0.7,
    labels={
        # "ClusterSize": "Number of clusters (k)",
        # "MeanSilhouetteScore": "Mean silhouette score<br>of MiniBatchKMeans"
        "ClusterSize": "Number of MiniBatchKMeans clusters (k)",
        "MeanSilhouetteScore": "Mean silhouette score"
    }
)

fig.update_xaxes(dtick=100, title_font=dict(size=title_font_size))

fig.update_yaxes(rangemode="tozero", tick0=0.0, dtick=0.01, title_font=dict(size=title_font_size))

# width = 500
width = 700
# height = 350
# height = 300
height = 0.6 * width

fig.update_annotations(font_size=title_font_size)

fig.update_layout(
    title_text="Long-reads",
    title_x=0.15,
    # title_y=0.95,
    template=template,
    # showlegend=False,
    width=width,
    height=height,
    title_font=dict(size=title_font_size+4)
    # legend_font=dict(size=18),
    # legend_grouptitlefont=dict(size=22),
)

# width = 650
# width = 1200
# height = 1200
    
# fig.write_image(
#     "Mean silhouette score of MiniBatchKMeans vs. K size - PacBio combined.svg",
#     width=600,
#     height=350,
# )

fig.write_image(
    Path(out_dir, "Mean silhouette score of MiniBatchKMeans vs. K size - PacBio - dense.svg",),
    width=width,
    height=height,
)

fig.show()

# %%
cols = min(facet_col_wrap, len(conditions), 5)
rows = ceil(len(conditions) / cols)
row_col_iter = list(product(range(1, rows + 1), range(1, cols + 1)))[: len(conditions)]

fig = make_subplots(
    rows=rows,
    cols=cols,
    subplot_titles=fixed_conditions,
    shared_xaxes=True,
    shared_yaxes=True,
    # x_title="Number of clusters (k)",
    # y_title="Mean silhouette score<br>of MiniBatchKMeans",
    x_title="Number of MiniBatchKMeans clusters (k)",
    y_title="Mean silhouette score",
    # column_titles=["Mean silhouette score", "Inertia score"],
    # row_titles=conditions
)

# min_x = None
# max_x = 0
# max_y = 0

marker_size = 3
line_width = 0.5

x = cluster_sizes

for (
    (row, col),
    condition,
    kmeans,
    kmeans_silhouette_scores,
) in zip(
    row_col_iter,
    conditions,
    conditions_kmeans,
    conditions_kmeans_silhouette_scores,
):
    color = color_discrete_map[condition]

    fig.add_trace(
        go.Scattergl(
            # go.Scatter(
            x=x,
            y=kmeans_silhouette_scores,
            mode="markers",
            marker=dict(
                color=color,
                # color=labels[clustered],
                # colorscale="Blackbody",
                line_width=line_width,
                size=marker_size,
            ),
        ),
        row=row,
        col=col,
    )

fig.update_yaxes(rangemode="tozero", tick0=0.0, dtick=0.01)

fig.update_layout(
    title_text="Long-reads",
    title_x=0.15,
    # title_y=0.95,
    template=template,
    showlegend=False,
    width=600,
    height=350,
)

# fig.write_image(
#     "Mean silhouette score of MiniBatchKMeans vs. K size - PacBio.svg",
#     width=600,
#     height=350,
# )

fig.show()
