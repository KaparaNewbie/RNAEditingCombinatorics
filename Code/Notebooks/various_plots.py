# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.4
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Imports

# %%
code_dir = "/private7/projects/Combinatorics/Code"

# %% papermill={"duration": 2.901153, "end_time": "2022-02-01T09:42:46.125355", "exception": false, "start_time": "2022-02-01T09:42:43.224202", "status": "completed"} tags=[]
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

# %% papermill={"duration": 0.071769, "end_time": "2022-02-01T09:42:43.049672", "exception": false, "start_time": "2022-02-01T09:42:42.977903", "status": "completed"} tags=["parameters"]
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

# %% papermill={"duration": 0.071769, "end_time": "2022-02-01T09:42:43.049672", "exception": false, "start_time": "2022-02-01T09:42:42.977903", "status": "completed"} tags=["parameters"]
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

# %% [markdown] papermill={"duration": 0.040192, "end_time": "2022-02-01T09:42:46.214429", "exception": false, "start_time": "2022-02-01T09:42:46.174237", "status": "completed"} tags=[]
# # Ploting utils

# %%
facet_col_spacing = 0.05
template = "plotly_white"
facet_col_wrap = 6
facet_row_spacing = facet_col_spacing * 6
zerolinewidth = 4

# %% papermill={"duration": 0.054755, "end_time": "2022-02-01T09:42:46.304499", "exception": false, "start_time": "2022-02-01T09:42:46.249744", "status": "completed"} tags=[]
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

# %%
# def n_repetitions_colormap(subcolors_discrete_map, condition, n_repetitions):
#     lowcolor, highcolor = subcolors_discrete_map[condition]
#     colors = pc.n_colors(lowcolor, highcolor, n_repetitions, colortype="rgb")
#     return {i: color for i, color in enumerate(colors, start=1)}

# %%
# n_repetitions_colormap(subcolors_discrete_map, "GRIA", 10)

# %% [markdown] papermill={"duration": 0.041741, "end_time": "2022-02-01T09:42:47.760215", "exception": false, "start_time": "2022-02-01T09:42:47.718474", "status": "completed"} tags=[]
# # Editing levels in PCLO

# %%
px.colors

# %%

# %%
known_sites_file = (
    "/private7/projects/Combinatorics/D.pealeii/Annotations/D.pea.EditingSites.csv"
)
# positions_files = [
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.positions.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.comp141882_c0_seq14.positions.csv",
# ]
positions_files = [
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.positions.csv.gz",
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

# %% tags=[]
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

# %% tags=[]
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

merged_ref_base_positions_df

# %%
fig = make_subplots(
    rows=2,
    cols=1,
    shared_xaxes=True,
    subplot_titles=["Known positions", "New positions"],
    x_title="Position",
    y_title="% editing",
    vertical_spacing=0.15,
)

symbols = ["circle", "square"]

current_study_min_x = end
current_study_max_x = start

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

        if col_suffix != "Known":
            _min_x = min([i for i, j in zip(x, y) if j != 0])
            current_study_min_x = min(_min_x, current_study_min_x)
            current_study_max_x = max(max(x), current_study_max_x)

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
                    name=col_suffix,
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

fig.update_xaxes(
    range=[current_study_min_x, current_study_max_x],
)

fig.update_yaxes(range=[0, 100], tick0=0, dtick=20)

fig.update_layout(template=template, height=600)

fig.show()

# %%
fig = go.Figure()

x = merged_ref_base_positions_df["Position"]

# for color, editing_percent_col, col_suffix in zip(colors, editing_percent_cols, col_suffixes):
for color, editing_percent_col, col_suffix in zip(
    px.colors.qualitative.Pastel[:3], editing_percent_cols, col_suffixes
):

    y = merged_ref_base_positions_df[editing_percent_col]

    fig.add_trace(
        go.Scatter(
            x=x,
            y=y,
            mode="lines+markers",
            marker=dict(color=color, opacity=0.5, size=6),
            name=col_suffix,
        )
    )

fig.update_layout(template=template, height=600)

fig.show()

# %% [markdown]
# # ADAR motif

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
            "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.positions.csv.gz",
            "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.positions.csv.gz",
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

# %%
# condition = "GRIA"
# strand = "+"
# positions_file = "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.positions.csv"

# positions_df = pd.read_csv(positions_file, sep=sep).drop(
#     ["CDS", "Phred"], axis=1, errors="ignore"
# )  # ignore error from droping these cols which only exist in Illumina data
# positions_df.insert(0, condition_col, condition)

# positions_df

# %%
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

# %%
# editing_sites_df.insert(
#     editing_sites_df.columns.get_loc("Start") + 1,
#     "End",
#     editing_sites_df["Start"] + 1
# )
# editing_sites_df

# %%
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

# %% tags=[]
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

# %%
merged_platforms_editing_sites_dfs[1]

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
out_file = None

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
multiple_logos_from_fasta_files(
    fasta_files, main_title, sub_titles, out_file, width=14, height=4
);

# %% [markdown] tags=[]
# # Pooled editing levels comparisons

# %%
known_sites_file = (
    "/private7/projects/Combinatorics/D.pealeii/Annotations/D.pea.EditingSites.csv"
)

# pacbio_positions_files = [
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.positions.csv.gz",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.positions.csv.gz",
# ]
pacbio_positions_files = [
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.positions.csv.gz",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.positions.csv.gz",
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

# %% tags=[]
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


# %% papermill={"duration": 4.052404, "end_time": "2022-02-01T09:42:53.176715", "exception": false, "start_time": "2022-02-01T09:42:49.124311", "status": "completed"} tags=[]
# todo retain nan rows and turn nans to 0?

both_df = merged_ref_base_positions_df.loc[
    merged_ref_base_positions_df["Edited"]
    & merged_ref_base_positions_df["KnownEditing"]
]
edited_df = merged_ref_base_positions_df.loc[
    merged_ref_base_positions_df["Edited"]
    & ~merged_ref_base_positions_df["KnownEditing"]
]
known_editing_df = merged_ref_base_positions_df.loc[
    ~merged_ref_base_positions_df["Edited"]
    & merged_ref_base_positions_df["KnownEditing"]
]

# all_df = merged_ref_base_positions_df.loc[merged_ref_base_positions_df["Edited"] | merged_ref_base_positions_df["KnownEditing"]]

fig = make_subplots(
    rows=1,
    cols=1,
    # subplot_titles=conditions,
    # shared_yaxes=True,
    # shared_xaxes=True,
    x_title="% editing",
    y_title="% known editing",
)

# df = merged_ref_base_positions_df.copy()
# df = merged_ref_base_positions_df.loc[merged_ref_base_positions_df["Edited"] & merged_ref_base_positions_df["KnownEditing"]]
x = all_df["%Editing"].fillna(0)
y = all_df["%EditingKnown"].fillna(0)

all_colors = all_df.apply(
    lambda x: editing_status_color(x["Edited"], x["KnownEditing"]), axis=1
)

x_both = both_df["%Editing"].fillna(0)
y_both = both_df["%EditingKnown"].fillna(0)
r, pv = scipy.stats.pearsonr(x_both, y_both)

fig.add_trace(
    go.Scatter(
        x=x,
        y=y,
        # name=condition,
        mode="markers",
        # marker_color="black",
        marker_color=all_colors,
        # marker_line_color=colors,
        # marker_line=dict(
        #     color=colors
        # ),
        # marker_color="white",
        marker_size=6,
        # marker_symbol="circle-open"
    ),
    row=1,
    col=1,
)

fig.add_annotation(
    row=1,
    col=1,
    x=20,
    y=85,
    xref="x",
    yref="y",
    text=f"<b>Pearson's r</b><br>p-val = {pv:.2e}<br>ρ = {r:.2g}",
    bgcolor="white",
    borderpad=4,
    opacity=0.7,
    showarrow=False,
)

fig.update_layout(
    title_text="Correlation between current & previously-reported editing levels",
    showlegend=False,
    template=template,
    width=600,
    height=500,
)

fig.update_xaxes(range=[0, 100])
fig.update_yaxes(range=[0, 100])

fig.show()


# %% papermill={"duration": 4.052404, "end_time": "2022-02-01T09:42:53.176715", "exception": false, "start_time": "2022-02-01T09:42:49.124311", "status": "completed"} tags=[]
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
edited_df

# %% papermill={"duration": 4.052404, "end_time": "2022-02-01T09:42:53.176715", "exception": false, "start_time": "2022-02-01T09:42:49.124311", "status": "completed"} tags=[]
fig = make_subplots(
    rows=1,
    cols=1,
    x_title="% current editing",
    y_title="% known editing",
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
names = ["Both", "Currently edited", "Known editing"]
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

x_all = all_df["%Editing"].fillna(0)
y_all = all_df["%EditingKnown"].fillna(0)
r, pv = scipy.stats.pearsonr(x_all, y_all)

fig.add_annotation(
    row=1,
    col=1,
    x=20,
    y=85,
    xref="x",
    yref="y",
    text=f"<b>Pearson's r</b><br>p-val = {pv:.2e}<br>ρ = {r:.2g}",
    bgcolor="white",
    borderpad=4,
    opacity=0.7,
    showarrow=False,
)

fig.update_layout(
    title_text="Correlation between current & previously-reported editing levels",
    legend_title_text="Editing status",
    # showlegend=False,
    template=template,
    width=600,
    height=500,
)

fig.update_xaxes(range=[0, 100])
fig.update_yaxes(range=[0, 100])

fig.show()


# %% papermill={"duration": 4.052404, "end_time": "2022-02-01T09:42:53.176715", "exception": false, "start_time": "2022-02-01T09:42:49.124311", "status": "completed"} tags=[]
fig = make_subplots(
    rows=1,
    cols=1,
    x_title="Editing status",
    y_title="% editing",
)

y_both = both_df["%Editing"].fillna(0)
y_edited = edited_df["%Editing"].fillna(0)
y_known_editing = known_editing_df["%EditingKnown"].fillna(0)
ys = [y_both, y_edited, y_known_editing]

names = ["Both", "Currently edited", "Known editing"]
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

fig.show()


# %% [markdown]
# # Isoforms in 80K coverage

# %%
condition_col = "Gene"
sep = "\t"
unique_proteins_first_col_pos = 14
platforms = ["PacBio", "Illumina"]

# %%
pacbio_distinct_unique_proteins_files = [
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/GRIA-CNS-RESUB.DistinctUniqueProteins.03.03.2023-15:36:38.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/PCLO-CNS-RESUB.DistinctUniqueProteins.03.03.2023-15:53:01.csv",
]

pacbio_unique_reads_files = [
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_reads.csv.gz",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_reads.csv.gz",
]
pacbio_unique_proteins_files = [
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz",
]

pacbio_conditions = ["GRIA", "PCLO"]

pacbio_color_sequence = px.colors.qualitative.G10
pacbio_color_discrete_map = {
    condition: color
    for condition, color in zip(pacbio_conditions, pacbio_color_sequence)
}

# %%
illumina_distinct_unique_proteins_files = [
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80K/comp141881_c0_seq3.Sampled80000.DistinctUniqueProteins.11.04.2023-16:36:40.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80K/comp141044_c0_seq2.Sampled80000.DistinctUniqueProteins.11.04.2023-17:21:47.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80K/comp140439_c0_seq1.Sampled80000.DistinctUniqueProteins.11.04.2023-16:34:07.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80K/comp126362_c0_seq1.Sampled80000.DistinctUniqueProteins.11.04.2023-17:27:11.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80K/comp141517_c0_seq1.Sampled80000.DistinctUniqueProteins.11.04.2023-18:09:19.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80K/comp141840_c0_seq2.Sampled80000.DistinctUniqueProteins.11.04.2023-16:40:11.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80K/comp141640_c0_seq1.Sampled80000.DistinctUniqueProteins.11.04.2023-16:35:41.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80K/comp140987_c3_seq1.Sampled80000.DistinctUniqueProteins.11.04.2023-16:51:06.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80K/comp140910_c2_seq1.Sampled80000.DistinctUniqueProteins.11.04.2023-16:46:22.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80K/comp136058_c0_seq1.Sampled80000.DistinctUniqueProteins.11.04.2023-17:00:54.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80K/comp141378_c0_seq7.Sampled80000.DistinctUniqueProteins.11.04.2023-18:02:54.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80K/comp141158_c1_seq2.Sampled80000.DistinctUniqueProteins.11.04.2023-16:38:37.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80K/comp140712_c0_seq3.Sampled80000.DistinctUniqueProteins.11.04.2023-16:55:55.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80K/comp141882_c0_seq14.Sampled80000.DistinctUniqueProteins.11.04.2023-17:08:14.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80K/comp141880_c1_seq3.Sampled80000.DistinctUniqueProteins.11.04.2023-17:47:57.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80K/comp141565_c6_seq3.Sampled80000.DistinctUniqueProteins.11.04.2023-17:55:20.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80K/comp141684_c0_seq1.Sampled80000.DistinctUniqueProteins.11.04.2023-17:41:03.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80K/comp141532_c3_seq11.Sampled80000.DistinctUniqueProteins.11.04.2023-17:37:06.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80K/comp141574_c0_seq3.Sampled80000.DistinctUniqueProteins.11.04.2023-17:13:57.csv",
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
    f"/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80K/reads.sorted.aligned.filtered.{chrom}.Sampled80000.unique_reads.csv.gz"
    for chrom in illumina_chroms
]
illumina_unique_proteins_files = [
    f"/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina80K/reads.sorted.aligned.filtered.{chrom}.Sampled80000.unique_proteins.csv.gz"
    for chrom in illumina_chroms
]

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
#     "ROBO2_HUMAN"
# ]

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
platforms_color_map = {
    platform: color_map
    for platform, color_map in zip(
        platforms, [pacbio_color_discrete_map, illumina_color_discrete_map]
    )
}

# ["circle", "square-dot", "diamond", "circle", "star", "star-square", "triangle-down"]
# symbols = ["diamond", "square"]
# symbols = ["diamond", "circle"]
symbols = ["star", "triangle-up"]
platforms_symbols = {
    platform: symbol for platform, symbol in zip(platforms, symbols)
}


fig = go.Figure()

for platform, condition in editable_aas_per_platform_and_condition:

    platform_and_condition_df = distinct_proteins_per_editable_aas_df.loc[
        (distinct_proteins_per_editable_aas_df["Platform"] == platform)
        & (distinct_proteins_per_editable_aas_df[condition_col] == condition)
    ]

    x = platform_and_condition_df["EditableAAs"]
    y = platform_and_condition_df["NumOfProteins"]
    
    x_mean = [x.iloc[0]]
    y_mean = [y.mean()]

    fig.add_trace(
        go.Scatter(
            x=x_mean,
            y=y_mean,
            mode="markers",
            marker=dict(
                color=platforms_color_map[platform][condition],
                symbol=platforms_symbols[platform],
                size=8,
            ),
            legendgrouptitle_text=platform,
            legendgroup=platform,  # this can be any string
            name=condition,
        )
    )

# correlate all x and y values
x = distinct_proteins_per_editable_aas_df["EditableAAs"]
y = distinct_proteins_per_editable_aas_df["NumOfProteins"]
r, pv = scipy.stats.pearsonr(x, y)

fig.add_annotation(
    x=40,
    y=25_000,
    xref="x",
    yref="y",
    text=f"<b>Pearson's r</b><br>p-val = {pv:.2e}<br>ρ = {r:.2g}",
    bgcolor="white",
    borderpad=4,
    font=dict(size=12),
    opacity=0.8,
    showarrow=False,
)

fig.update_xaxes(
    range=[0, distinct_proteins_per_editable_aas_df["EditableAAs"].max() * 1.1]
)
fig.update_yaxes(
    range=[0, distinct_proteins_per_editable_aas_df["NumOfProteins"].max() * 1.1]
)

fig.update_layout(
    # legend_font=dict(size=8),
    # legend_grouptitlefont=dict(size=10),
    template=template,
    xaxis_title="Editable amino acids",
    yaxis_title="Distinct proteins (avg)",
    width=700,
    height=650,
)

fig.write_image(
    "Distinct proteins vs. editable AAs (80K) - PacBio vs. Illumina.svg",
    width=700,
    height=650
)

fig.show()


# %%
