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

# %% [markdown]
# # Imports

# %%
code_dir = "/private7/projects/Combinatorics/Code"

# %% papermill={"duration": 2.901153, "end_time": "2022-02-01T09:42:46.125355", "exception": false, "start_time": "2022-02-01T09:42:43.224202", "status": "completed"} tags=[]
import sys
from functools import reduce
from itertools import chain, product, combinations
from math import ceil
from pathlib import Path
from multiprocessing import Pool

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
from Alignment.alignment_utils import count_reads, count_reads_in_unaligned_bam

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
template = "plotly_white"

# %% papermill={"duration": 0.054755, "end_time": "2022-02-01T09:42:46.304499", "exception": false, "start_time": "2022-02-01T09:42:46.249744", "status": "completed"} tags=[]
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

# %% [markdown] papermill={"duration": 0.041741, "end_time": "2022-02-01T09:42:47.760215", "exception": false, "start_time": "2022-02-01T09:42:47.718474", "status": "completed"} tags=[]
# # Editing levels in PCLO

# %%
px.colors

# %%
known_sites_file = "/private7/projects/Combinatorics/D.pealeii/Annotations/D.pea.EditingSites.csv"
positions_files = [
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.positions.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.comp141882_c0_seq14.positions.csv"
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
    px.colors.qualitative.Vivid[5]
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
    positions_df = (
        pd.read_csv(positions_file, sep=sep)
        .drop(["CDS", "Phred"], axis=1, errors="ignore")  # ignore error from droping these cols which only exist in Illumina data
    ) 
    positions_df.insert(0, condition_col, condition)
    positions_dfs.append(positions_df)
positions_dfs[0]

# %% tags=[]
ref_base_positions_dfs = []

for positions_df, strand in zip(positions_dfs, strands):
    ref_base = "A" if strand == "+" else "T"
    alt_base = "G" if strand == "+" else "C"
    positions_df = positions_df.loc[positions_df["RefBase"] == ref_base].sort_values(["Chrom", "Position"]).reset_index(drop=True)
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
merged_ref_base_positions_df = ref_base_positions_dfs[0].loc[:, ["Chrom", "Position", "EditingFrequency"]].merge(
    ref_base_positions_dfs[1].loc[:, ["Chrom", "Position", "EditingFrequency"]],
    how="outer",
    on=["Chrom", "Position"],
    suffixes=[f"_{condition}" for condition in conditions],
)

merged_ref_base_positions_df = merged_ref_base_positions_df.merge(
    known_sites_df.loc[known_sites_df["Chrom"] == chrom, ["Chrom", "Position", "EditingFrequency"]].rename(columns={"EditingFrequency": "EditingFrequency_Known"}),
    how="outer",
    on=["Chrom", "Position"],
)

for editing_freq_col, editing_percent_col, edited_col in zip(editing_freq_cols, editing_percent_cols, edited_cols):
    merged_ref_base_positions_df[editing_percent_col] = merged_ref_base_positions_df[editing_freq_col].fillna(0) * 100
    merged_ref_base_positions_df[edited_col] = merged_ref_base_positions_df[editing_percent_col] > 0
    
merged_ref_base_positions_df = merged_ref_base_positions_df.drop(editing_freq_cols, axis=1)

# merged_ref_base_positions_df["Position"] = merged_ref_base_positions_df["Position"].astype("str")
merged_ref_base_positions_df = merged_ref_base_positions_df.sort_values(["Chrom", "Position"]).reset_index(drop=True)

merged_ref_base_positions_df

# %%
fig = make_subplots(
    rows=2, cols=1, 
    shared_xaxes=True,
    subplot_titles=["Known positions", "New positions"],
    x_title="Position",
    y_title="% editing",
    vertical_spacing=0.15
)

symbols = ["circle", "square"]

current_study_min_x = end
current_study_max_x = start

for row in [1, 2]:

    if row == 1:
        df = merged_ref_base_positions_df.loc[merged_ref_base_positions_df["Edited_Known"]]
    else:
        df = merged_ref_base_positions_df.loc[~merged_ref_base_positions_df["Edited_Known"]]
   
    for color, editing_percent_col, edited_col, col_suffix in zip(px.colors.qualitative.Pastel[:3], editing_percent_cols, edited_cols, col_suffixes):
        
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
                    name=col_suffix
                ),
                row=row,
                col=1
            )
        else:
              fig.add_trace(
                go.Scatter(
                    x=x,
                    y=y,
                    mode="lines+markers",
                    marker=dict(
                        color=color,
                        opacity=0.7,
                        size=4
                    ),
                    showlegend=False
                ),
                row=row,
                col=1
            )
                
fig.update_xaxes(
    range=[current_study_min_x, current_study_max_x],
)
                
fig.update_yaxes(
    range=[0, 100],
    tick0 = 0,
    dtick = 20
)
                
fig.update_layout(
    template=template,
    height=600
)
    
fig.show()

# %%
fig = go.Figure()

x = merged_ref_base_positions_df["Position"]

# for color, editing_percent_col, col_suffix in zip(colors, editing_percent_cols, col_suffixes):
for color, editing_percent_col, col_suffix in zip(px.colors.qualitative.Pastel[:3], editing_percent_cols, col_suffixes):
    
    y = merged_ref_base_positions_df[editing_percent_col]
    
    fig.add_trace(
        go.Scatter(
            x=x,
            y=y,
            mode="lines+markers",
            marker=dict(
                color=color,
                opacity=0.5,
                size=6
            ),
            name=col_suffix
        )
    )
    
fig.update_layout(
    template=template,
    height=600
)
    
fig.show()

# %% [markdown]
# # ADAR motif

# %%
