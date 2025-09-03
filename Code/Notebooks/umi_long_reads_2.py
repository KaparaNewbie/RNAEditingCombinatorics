# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.2
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Settings
#

# %%
code_dir = "/private7/projects/Combinatorics/Code"

# %%
# %load_ext autoreload
# %autoreload 2

# %%
from collections import Counter, defaultdict
from itertools import chain, product, repeat
import itertools
import multiprocessing as mp
from multiprocessing import Pool
from pathlib import Path
import subprocess
import copy
import sys

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.io as pio
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pysam
from Bio import Align
from Bio.Seq import Seq, reverse_complement, complement
from icecream import ic
import matplotlib.pyplot as plt
import networkx as nx
import more_itertools
from IPython.display import clear_output

sys.path.append(str(Path(code_dir).absolute()))
from Alignment import umi_processing

# %%
pio.templates.default = "plotly_white"

# %%
seed = 1892

# %%
pd.set_option("display.max_columns", 500)

# %%
unmapped_bams_dir = Path(
    "/private7/projects/Combinatorics/D.pealeii/Data/RawWithUMIs/30-1097162729/CCS"
)
mapped_bams_dir = Path(
    "/private7/projects/Combinatorics/D.pealeii/Alignment/UMILongReads"
)

# %%
unmapped_bam_files = list(unmapped_bams_dir.glob("*.bam"))
unmapped_bam_files

# %%
mapped_bam_files = list(mapped_bams_dir.glob("*.bam"))
mapped_bam_files

# %%
merged_bams_dir = Path(
    "/private7/projects/Combinatorics/D.pealeii/Alignment/UMILongReads.MergedSamples"
)

# %%
mapped_merged_bam_files = list(reversed(list(merged_bams_dir.glob("*.bam"))))
mapped_merged_bam_files

# %%
min_read_quality = 0.998

# %%
mapped_merged_filtered_bams_dir = Path(
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples"
)

# %%
merged_annotated_reads_files = [
    Path(mapped_merged_filtered_bams_dir, "ADAR1.Merged.r64296e203404D01.aligned.sorted.MinRQ998.reads.csv.gz"),
    Path(mapped_merged_filtered_bams_dir, "IQEC.Merged.r64296e203404D01.aligned.sorted.MinRQ998.reads.csv.gz")
]

# %%
merged_old_to_new_reads_files = [
    Path(mapped_merged_filtered_bams_dir, "ADAR1.Merged.r64296e203404D01.aligned.sorted.MinRQ998.OldToNewReads.csv.gz"),
    Path(mapped_merged_filtered_bams_dir, "IQEC.Merged.r64296e203404D01.aligned.sorted.MinRQ998.OldToNewReads.csv.gz")
]

# %%
reads_with_recognizable_barcodes_dir = Path(
    "/private7/projects/Combinatorics/D.pealeii/Alignment/UMILongReads.ReadsWithRecognizableBarcodes.MergedSamples"
)
reads_with_recognizable_barcodes_dir.mkdir(parents=True, exist_ok=True)

# %% [markdown]
# # Biological settings
#

# %% [markdown]
# > ADAR1_RT_Div
#
# CTAGCAATGCTTCAGGCTGTGGNNNNNNNNNNNNgcagtctggaaaggaatggccatc
#
# > IQEC_RT2_Div
#
# CTAGCAATGCTTCAGGCTGTGGNNNNNNNNNNNNCTGGGAGAAGAGCAGATGACTG
#

# %%
primers_dict = {
    "PCR": Seq("CTAGCAATGCTTCAGGCTGTGG"),
    "IQEC": Seq("CTGGGAGAAGAGCAGATGACTG"),
    # "ADAR": Seq("gcagtctggaaaggaatggccatc".upper()),
    "ADAR1": Seq("gcagtctggaaaggaatggccatc".upper()),
}
primers_dict

# %%
STRANDS = ["+", "-"]

# %%
genes = ["ADAR1", "IQEC"]
GENES = genes

# %%
chrom_per_gene_dict = {
    "ADAR1": "comp134400_c0_seq1_extended",
    "IQEC": "comp141565_c6_seq3",
}

# %%
gene_by_chrom_dict = {v: k for k, v in chrom_per_gene_dict.items()}

# %%
genes_orfs_boundries_dict = {"ADAR1": [0, 3741], "IQEC": [988, 4195]}

# %%
genes_seq_dict = {
    "ADAR1": "ATGGCGAACTCTAATCTATCTTCTTCCATGAATCAAAACATGGCTCACATGGGTGGTGGGAGTTTAGTAAATGGCTACTACAAACAAGTACCATATTCAGGTGGAAGAAGCCGAAATGCAAGTGGTAGGTCCGGTAGTCGTGGCCGCGGCAAACCTGCAGTCAGGGAAACTAGTTTGAATGTTCATCCAGAATGGGAAGAGCGTATTGTAAACTACCTTGCCCACAAAACTCATCCTGTAAAGACCATGGAGCTGGCACGCCTCGTGAATGTTCGCTCACGCAAAGAAGTGAATCCCACTTTGTACAGCATGGACAGACGAGGTCTTATCAGGAAACATGGAATGCAACCTCCGACATGGGTAATTGCAGACCCACCCCAATCTCACGGCGGATACAACCAAAATGAGACACACTATTCAAGTAGCCCAGGAATTTACCAGCATAGTCCGGTTTCGAGAACTCCTCAGAACTTTTATCCCAATAATCGAGAGAGTTATCGAGGACACAAAGCTCCAAATAGTAATTACCCGCGCTCCAAACGGACTTCATACAGGAATGACTGGCATAACTTTTGCTCCCCTCCATCTCACATGTACCCAGAGGGCAAAAATGAATCTTTGATCTATAGTCACAGTAACAAAGATAATGAGATGTTATCAATGGGAAACGCTAGTTCTCCAAACAGATTGCGGTCTGAAAGTTGTAGCCCAGATGAATGTGAAGCGCTGGAAACAGCCATTGCAAATTATCCACATAGTGGTGGGTATGGCCAGGGTTACTCAGGACATTTCCCTGTAACCCCAGAACCAAAGCAAACTGGCAAGAGAAGGAGAAATTGTGATAACTTTGGTTTACAACAACATCCATCAAATGGTACGATGCCGATGAAAAACAGTGAGAAGATCCAACAGAAAAAATTGGAATTTCAGGATGAAAGATATTATGATGCAAGCTACCCATCCTATTCTGGAAATTTTACTATGAACTATGCAAATCATTTTGATATGCCTGTCTATCATCCGATAGACCGGGAGGACCCGAAGGATACTCCACCTCCGTCACGTGTGTCAATGGACTTGATAAAAAAGGATTCCAAGGACATCTCGTCACATGAACGAATCTCTCCCAAGAGGAATTCAAACAGTAAGGGTTGTAATTCTGACGCTCATGACCCACAGGCAAGAGTTATTTCCTTCCTGGATAAAACTATGAATGGCTCTGCGAAGTCACGAGAAATCGCCAAGCATACAAGTCTTTCTCTGGAAGATACCCAGAAGATATTGCATAGTTTGTGTAAGAAGAAAATAGTCGCAACAATTGGTGAAGATATCTACATAATGGCTAAAAATGCAGCCAGTTATGAGACTGAGATTCCAGCAGGAAAAAACTCCTCATCAAACATGAATTCAAACATGGCACGCCAGTTCTCCAGTGGAAATCGGCAGCCCCCTGCGCCCCCACATGTACTATTGGCGGATAATGGCATCAATTCCGGCAGCATGAAAAACGTTTATTTCCAGGGTAATAATGCTCCCAAACAATCTGGGTCCAACTCGAGTGAATCAAAATCAGCACAGAGCCAGGTGGGCAGAAGCCCTCATCTACCCCCTTCCCCTCATGAACTATTAGCAAAGGACCCAATTTTCAAGGGAGACATTACTGCACCCAATACAAACGCTTCAAAGGACTACAACCAGTCGTCATCATCTTCGTCAGCATCCTTGTCGTCCTCAACGTCAAAGAACTCAAGGTGGAATAGCAACACTGCAGCGACAGAGAGTTCCAGAGCTCCAAACACGACCTCTGCTTCAACATCGTCAACTACATCATTTGCTCCCACTCCTAGTAAGTCTGCCTCTAATTCAAAACAGACTGCTCCTAGTCCCAAGCAACCATCTCCAAGTCCTAAGCAGAACACCCCTAAGAGTTCCAAGAGTTCCAAAAGTTCCAAGCAGAGAGCCACAAGCCCCAAACAAAACAGCACTCCTAGCTCCCAGGCGTCCTCTCAGTCAAACTCCAATACTACTACAACTGCCACCTCAAGCAGCAGCAAAAATAATAAAAATAACAACAATAACAACACCTCAGTAGAGAATTTGCAAGATGCCCTCAAAAATGTGTCTATCTCGTCCCCGACTGAGACTACTGAGAGCAAAACGCCCACATTGGCCGAGATCAAGGCGGCAGCAGTGGCGCAGGCGTTGGCCGACAAAGCGGCTGAGAAAGGAGCTGACAAGTCTGGTACTGATTCATTGGCACCAAATCTACAGATCACCTCAGAAAGTTTCGCCGCTCTCAACAAAAATCCAGTTAGCGCACTGATGGAATATGCTCAACAGCGACACTTACCCGTTGAATTTAAGCTTTTGTCACACAGAGGACCTTCTCATCGACCGTTGTTCAAATTTGCCGTGATTCTTGGTAAACGCCAGTTCCCCAGTATGGAGTGCAACAGTAAGAAGGATGGTAAGAAAGAGGCAGCCGATCTGACATTGCGCATTCTCATTGCTGAAGGACAGTATCAACTGGAGAACACCGTCTCAGCATTGAAAACAATTCCACCTGCTGAAATGACACATTTCGACAAAATGGCTGCCTTGAGTCACCAGGCATTTAACAACATTGCCTTGCAAATCCCTGAGAACCTTGCTGGGAGAAAGGTCATCGCTGCTTTGGTGATGAAGCGATCACCAACGGATACGGGAATTGTTATCAGTGTTGGAACTGGTAACCGCTGTTTAACCGGTGATCATTTGAGTTTGGAAGGCAACAGTGTCAATGACTCTCATGCTGAAATAATCACACGCCGAGGTTTTCTGAGATATCTGTACAAACATTTACTGGAGTATGATCCCGAAAAACCCCATGACCTATTTGAGAAAGGTGAACGTAGTCTTTGCCGGATAAAAACCAACATTACATTCCATCTGTATATATCAACTGCTCCTTGTGGTGATGGAGCACTTTTTTCACCCAGGGATACCGACTCCAGTAATGTGAAAGTGGATGAGGAAAATAAGCACGTCCATAATCCGACTTTTTCAAGCAGTGTTCAGGGATTGCTGAGAACCAAAGTGGAAGGAGGTGAAGGGACCATTCCAATAGATGCTGATTTCACTGAACAAACATGGGATGGAATTCAACGAGGTGAAAGATTGCGCACAATGTCATGTTCAGATAAAATATGTCGATGGAACGTTGTTGGTCTGCAAGGAGCTTTGCTTAGTCACTTTGTGGAACCAATCTACCTGGAATCTCTGACATTAGGTTATCTTTATGATCATGGCCACTTAGCACGAGCTGTTTGCTGCCGTATTGAACGGGGAGAGGCCTCTGTCAACCAACTACTACCTGAGGGCTACCGATTGAACCATCCTTGGCTTGGCAGAGTTACTGCTTGTGATCCACCTAGAGAAACCCAAAAGACGAAATCGTTAAGTATCAACTGGTGCTATGATGATGAAAAGTCTGAAGTTCTCGATGGTACAGCAGGCATCTGTTACACAGCGATTGAGAAAAATCTCTTCTCTCGCTTAACAAAGCACAGCTTATATGAAGAATTCAAAAAAGTGTGTCAGAAATTTGAACGCGAGGACTTGATGAATGTCACTTCTTACAACAAAGCCAAGATGATGGCCATTCCTTTCCAGACTGCCAAAAACGTAATGTTGAAAAAACTCAAGGAAAACAACTGCGGAACTTGGGTGTCAAAACCTATTGAGGAGGAGATGTTTACGTAG",
    "IQEC": "AAAAAAAAGAAAAGAAAAAGAAAAGAGGTTTTTTTACTTTATATTATTAATATATATTTGTAATTAACCTACCAAAGGCCTAATTGTGGGCATTTCCGGATTTTACACTGTTGCTCGCTTGTTTTCTTTCTCCCGTGACAAACCAAATGAATAGAACATTAAAAAAGGGGAAATGTTTCTTACTTGTATGTTGCCGATTACACTACGACCATGGTCACAAATAGACACACTATGGTCTACGGCCTGGACGCCCACAAAAAGCTCGCTTACCCCTTTTACACCTGTGGCTGGACGTCTTTCACCGAAAAAAGAAACGATACTGGTCATATCTTCACCCCCAAATTTTGACTAAACAAACGCAAAAATCACTAACACAATCGATAACAACTAGTACTTACTACTACTATATTCTCTTATTTAGAATTATTACTTACTACCTTCAATCATTCAGTGGTTGTCTTATCAACCAAAAAATTATAATACAAGCTGTAGTAAGATTTTTCAGTCACATAGTCTTTTTTTTGTATATTTTTTAGTAACAATTAAGTGGGTAACTACTGTTTGTTAACAACACTTACTTAACTAATGATTTACAACTAGTTATATTTACTTGCAATTAAATGACTGCTGTCGATTATATTTTTACACACACACAAACACTAAAACACATGACATAGTACCCTAATACATTATTATTACTATTACACCCACGTACTATAAACAATCTAATTTATTCTTTACTGCTGTCATTTTATTTATTTCCTTCCCCAGAACCCTTGTAATTATTAAATTGGTTAATTAATTACCCTACTGCTTATTAGCGCTGTGGACTACCCAAGGATACTAATTTCTGAAAACATCGACGCCTTTTAAAAATTCACTAGGACGGATTATTATTACCAATTATAATTTAATAACCCCAAAACGTTTTATTTCCATGTTGTGCCCTATTCAAGATATTATTTTTATCACAAATACTTGAAACTAAACCAGCACCAGCCAAATCAACACGATGGAATTACCACCTCAGGAAGGTGTCACAAGCGCCAGGTCAAGTGAGGACCTACTTCAAGGAAAGGCGCGTCCTGTTGCCAGTGGAAACGGCTGCAGCAAAGAGCTTGGCCGCCCAAAGGGAGTTAATCAAAACAGCTCGACTGCTGGTGCTGCAGTCCCCAATGGCACCACCCCCTCCCCTCTCCCCACCCCAACACCAACAGCAACCAATGCCCAGCGGGCGCAGGTTCAAAGGCAACGCTCGTTCAAAAGAACAGACGATGAACTAATAAAACGGTCGCGCGTGCAGTCCAGCTCCTATGAGCTCTCCCAAGACTTGCAAGAAAAGCAATTGGAAATGCTCGAGCGGAAATATGGCGGAACAATACGAGCTCGACGTGCAGCGCGAGTAATCCAGAGAGCGTTTCGCAAATATTGTATGAATAAGAATTTTGAGAAATTGAGAAACTCAGTCGGCGAACGTCGATTGTCCAAACGGCTATCCGAACTTGGACGGAGCAATACCGTTTGGACTGACCGCATAAGTTCGGATATGCAGTTTACAATGGGAGGTACAGCATTACAAAATACAGACTCAAAGATTGACCAATTGTCACAGGAGAATGTGCGTAAGCGAATCGCAGAGATTGAATCCCGAAAGAATTTGGATAACAAGCTCCACCACCAATTGTCTTTAGACCCTTCATTGAAGATTGACAACCTCAAGCGAACACAACGAAAGGAAAGGAGGCGTTTGGAAAGGTCGATGGAGATAGATTTAACGGACATTCCCACCGAAAACGTGAAGGAGGAGGAGCGAGAAGAGAAATGTAGTGCCCAGTTGGCTGCCTCCATCGCCGCGGAAACGAATAACAATAGAAACTCTTACCCAGAGCTTAATGACAGTAGTGCCAGTGACTCCCCACAAGCAACGCCCGTGGAGAGCACAGTTGATTTGCATAGTTTAGATTTTGAAAATTTGCTTGAGAGCAAAGAGACTGATATTTTGACAGATAGTTTCCATAGTGACAGTATCCATAGCGATGGGAGCCAGGAGACTGGGTCCCTGATAGGCCATCACAAGCCTCAAGCCCTGTCTCTGTCAATGGAACTTGAGCCCCGGTTGCACAATGTGGTTAGTAAAGATCTGCATGGGTCCTGTGATAGTATGGGTAGGTTGTACCCGGAGTCGGGGGTCTCAGAGATCCCCTCTGAGGTCCAGATCAAAGTTGATCTGGCCTCCCCCGAGGACACCGGACTGGACGATCAGAAGACCCTGAATGAAGAGACTGTAAAATATTACATGAACACAAAGGTGAAATTACGTGGCGGTGGTAGTAGTGGTGTTAATAGCGGCGGTGGAGGTGGCGGGGGCAGCGTACGTGATATTAAGGACAAGAAAGGTGTGATGCCACACCGAGGACCAGAGGCATCTCCGATTTGGAAACGAAAGAGTGCCACTAACTCAACAGTGTCAGTCGTGAATAAAGGTGAGGTGAAGCGAATGAGCAACATCTCAGAGACGAGCGAACCAGAAAGTTTGGATGGCCAGTGTTCAAGTTCCCCCAGCTCTGATAATGTGAGTTCGGAGAACATAAGTATAGGCAGTGAGACTAGTATAAGTTATCAAAGGAAAATCCGTATCAGCGTCACAAACGACCAACACATGACGAGGATTGGCGATCGGGAAAGGAAGCGACTGTATAGAATTGGGTTGAATTTATTTAACAAAAAACCAGAAAAGGGCCTCCGCTTCCTCATAGAACGTGGTTTTGTTGATCAATCTCCACCAAATATCGCCAAATTTCTGATTACACGAAAAGGACTTAGTAAACAGATGATTGGGGAGTACCTTGGGAACTTACAGATTCAACTGAACCAAGATGTCTTAGATTTCTTTGCTGAGGAAATAGATTTGTCAGGACTCCAAGTGGATGTAGCACTTCGAAAGTTCCAGTCCTTCTTCAGAATGCCGGGAGAAGCACAGAAAATTGAAAAACTAATGGAGGCATTTGCCCATCGATATTGTATATGCAACCCTGACCAAATAAAAAACTTCCGAAACCCAGACACTATCTTCCTGCTGGCCTTTGCTATCATCATGCTCAACACAGACCTGCATAATTCCAACATCAAAACAGAAAAGCGAATGAAACTAGAAGACTTTATCAAAAACCTCAGAGGAATTGATGAAACAGAAGATATTGACCGTGACTTATTACGAGGAATATATGAACGGATCCGTGCCCAAGAGTTTCAGCCAGGAGTAGATCATCAGTCTCAAGTAATGAAGGTTGAACAGACCATTGTAGGCAAAAAACCACAACTGGCTTTACCACACAGACGCCTTGTTTGCTACTGCCGCCTTTATGAAATTCATGATCCCTACAAAAAAGAGAAAATTGGTCTTCATCAGCGTGAAGTATTTCTTTTTAATGACATCATTCTTATTACAAAAATTTTTAGCAAGAAAAAAACAGGAATCACCTACTCCTTCAAGCAGTCCTTCTCCTTGTATGGCATGCAGGTGTATCTGTTTGGAACAGCTCATTACCAATATGGTATTAGGTTGTCCAGTAACGGCAAGGTGCTAATTACGTTCAACGCTAGGAATGACCATGATAGGCAGAAATTTGTCGAAGATCTCAAAGAAGCCATTTTAGAGACCACCGGCTTAGAAAGCTTACGGATCGAAGAGGAATTGCAACGACATCGTGCCACTCACAACACCATTGATGACGACTCTCGAGTCTTGATTTATGACGTTTTAAAGCCATCAGATCCTGCGGTAAATAGGCTCTCTGCACCTGAGTGTGGTTTGAAGAAGTCTGCTTTGTCCAATTCTCTCCTTGACCTCTGTGAAGGTCAAGGAATGAAGAGAGAGGACAGTGGAGGATCTCTCGACAGTGGTGTGATTCTGTTTGGTGGTGTGAAAGTACGCATATCTCAGTGCAAACTGAGTGAAATCAGTCATCTGCTCTTCTCCCAGGCGTCAGTGGGAACGGTGAGCGGTCACAGTCGTGATGATTCTGTCTCCACCATTCTAAATCCTCATCGGGGATCCATTCTTCATACGAACCCAGCCACCCATGAACGCTGCCATTCTTCAATATCTGTTTCTTCAACACGGCCCAACAATCCCAAGCGAATGACAGAACCGGCCACTCTTAGCCGAAGTAGCAGCAGCAGCAGCAGCAACAAC",
}

# %%
expected_barcode_locations_on_gene_dict = {}

for gene in genes:

    aligner = Align.PairwiseAligner(
        # match_score=1.0, mismatch_score=-2.0, gap_score = -2.5,
        mode="local",  # otherwise we'd get scattered matches of the primer across the read
        scoring="blastn",
        # scoring="megablast"
    )

    barcode_seq = primers_dict[gene]
    gene_seq = genes_seq_dict[gene]

    alignments = aligner.align(gene_seq, barcode_seq, strand="-")
    assert (
        len(alignments) == 1
    ), f"Expected exactly one alignment for {gene}, got {len(alignments)}"
    alignment = alignments[0]
    print(gene)
    print(alignment)
    btr_read_coords = alignment.aligned[0][0]
    expected_barcode_locations_on_gene_dict[gene] = btr_read_coords

expected_barcode_locations_on_gene_dict


# %% [markdown]
# # Utility functions

# %%
def get_read_object_from_bam_files(
    read_name: str, bam_files: list[Path], sample: str = None
):
    if sample is not None:
        bam_files = [bam_file for bam_file in bam_files if sample in bam_file.name]
    for bam_file in bam_files:
        with pysam.AlignmentFile(
            bam_file,
            "rb",
            threads=10,
        ) as samfile:
            for read in samfile:
                if read.query_name == read_name:
                    return read
    return None


# test the function with a specific read name
read_name = "m64296e_241222_071206/74516384/ccs"
get_read_object_from_bam_files(read_name, mapped_bam_files)

# %% [markdown]
# ![image.png](attachment:image.png)

# %% [markdown]
# # Unmapped BAMs
#

# %%
unmapped_bam_files = list(unmapped_bams_dir.glob("*.bam"))
unmapped_bam_files

# %%
# !samtools view -c --threads 10 /private7/projects/Combinatorics/D.pealeii/Data/RawWithUMIs/30-1097162729/CCS/LP1ADAR1.r64296e203404D01.hifireads.bam

# %%
# def count_unique_reads(bam: Path, threads: int = 1):
#     with pysam.AlignmentFile(bam, "rb", threads=threads) as samfile:
#         unique_reads_names = {read.query_name for read in samfile}
#     return len(unique_reads_names)


# def count_reads(bam: Path, threads: int = 1):
#     with pysam.AlignmentFile(bam, "rb", threads=threads) as samfile:
#         reads_names = [read.query_name for read in samfile]
#     return reads_names

# %%

# %%
# bam_file = unmapped_bam_files[0]
# bam_file

# %%
# # !samtools view $bam_file | head -n 1

# %%
# with pysam.AlignmentFile(
#     bam_file,
#     "rb",
#     threads=10,
#     check_sq=False,
#     require_index=False,
#     # index_filename=str(Path(bam_file.parent, f"{bam_file.name}.pbi")),
# ) as samfile:
#     reads = [read for read in samfile]
# reads[0]
# # print(reads[0]|)

# %%
# read = reads[0]

# %%
# # fn: Number of passes in PacBio sequencing (or similar, platform-dependent).
# read.get_tag("fn")

# %%
# # np: Number of passes through the template (PacBio circular consensus reads).
# read.get_tag("np")

# %%
# Seq("CTAGCAATGCTTCAGGCTGTGG")

# %%
# read.query_sequence

# %%
# read.get_tag("rx")

# %%
# read.get_tag("mi")

# %% [markdown]
# bx:B:i,16,16   
# bc:B:S,57,57   
# bq:i:100   
# cx:i:12   
# bl:Z:TACTGACGCTGAGTAT   
# bt:Z:ATACTCAGCGTCAGTA
#

# %%
# read.get_tag("RG")

# %%
# # bc: Barcode sequence (raw or uncorrected).
# read.get_tag("bc")

# %%
# # bq: Base qualities of the barcode sequence.
# read.get_tag("bq")

# %%
# # bx: Barcode sequence, often corrected (used in 10x Genomics or similar applications).
# read.get_tag("bx")

# %%
# # bt: Barcode tag quality or type (platform-specific).
# read.get_tag("bt")

# %%
# len(read.get_tag("bt"))

# %%
# # bt: Barcode tag quality or type (platform-specific).
# read.get_tag("bl")

# %%
# len(read.get_tag("bl"))

# %%
unmapped_bam_dfs = []
for bam_file in unmapped_bam_files:
    sample = bam_file.name.split(".")[0]
    gene = sample[3:]
    repeat = sample[2]
    with pysam.AlignmentFile(
        bam_file,
        "rb",
        threads=10,
        check_sq=False,
        require_index=False,
        # index_filename=str(Path(bam_file.parent, f"{bam_file.name}.pbi")),
    ) as samfile:
        reads = reads = [read for read in samfile]
        reads_names = [read.query_name for read in reads]
        raw_barcods_seqs = [read.get_tag("bc") for read in reads]
        corrected_barcode_seqs = [read.get_tag("bx") for read in reads]
        barcode_tags = [read.get_tag("bt") for read in reads]
        # unique_reads_names = set(reads_names)
        # ic(bam_file.name.split(".")[0], len(reads_names));
        reads_seqs = [read.query_sequence for read in reads]
        reads_lengths = [len(seq) for seq in reads_seqs]
        df = pd.DataFrame(
            {
                "Sample": sample,
                "Gene": gene,
                "Repeat": repeat,
                "Read": reads_names,
                "RawBarcodeSeq": raw_barcods_seqs,
                "CorrectedBarcodeSeq": corrected_barcode_seqs,
                "BarcodeTag": barcode_tags,
                "Seq": reads_seqs,
                "ReadLength": reads_lengths,
            }
        )
        unmapped_bam_dfs.append(df)
concat_unmapped_bams_df = pd.concat(unmapped_bam_dfs, ignore_index=True)
concat_unmapped_bams_df

# %%
concat_unmapped_bams_df["CorrectedBarcodeSeq"][0]

# %%
concat_unmapped_bams_df["CorrectedBarcodeSeq"].astype(str)

# %%
concat_unmapped_bams_df.groupby(
    [
        "Gene",
        "Repeat",
    ]
).size().reset_index(name="NumOfReads")

# %%
concat_reads_per_barcode_tag_df = (
    concat_unmapped_bams_df.groupby(
        [
            "Gene",
            "Repeat",
            "BarcodeTag",
        ]
    )
    .size()
    .reset_index(name="NumOfReadsPerBarcodeTag")
)
concat_reads_per_barcode_tag_df

# %%
concat_reads_per_barcode_tag_df.groupby(["Gene", "Repeat"]).size().reset_index(
    name="NumOfUniqueBarcodeTags"
)

# %%
# concat_reads_per_corrected_barcode_seq_df = (
#     concat_bams_df.groupby(
#         [
#             "Gene",
#             "Repeat",
#             "CorrectedBarcodeSeq",
#         ]
#     )
#     .size()
#     .reset_index(name="NumOfReadsPerCorrectedBarcodeSeq")
# )
# concat_reads_per_corrected_barcode_seq_df

# %% [markdown]
# # Mapped BAMs
#

# %%
# mapped_bam_files = list(mapped_bams_dir.glob("*.bam"))
mapped_bam_files


# %%
# # !samtools view -c --threads 10 /private7/projects/Combinatorics/D.pealeii/Alignment/UMILongReads/LP2IQEC.r64296e203404D01.aligned.sorted.bam

# %%
# mapped_bam_file = mapped_bam_files[0]
# mapped_bam_file

# %% jupyter={"source_hidden": true}
# with pysam.AlignmentFile(
#     mapped_bam_file,
#     "rb",
#     threads=10,
# ) as samfile:
#     mapped_reads = [read for read in samfile]
# mapped_reads[0]
# print(mapped_reads[0])

# %% jupyter={"source_hidden": true}
# seq = read.query_sequence
# len(seq)

# %% jupyter={"source_hidden": true}
# soft_clipped_seq = read.query_sequence
# len(soft_clipped_seq)

# %%
# primers_dict

# %% jupyter={"source_hidden": true}
# read = mapped_reads[0]
# print(read)

# %%
# read.get_cigar_stats()

# %%
# read.get_aligned_pairs(with_cigar=True)[0][2] == 4

# %%
# read.get_aligned_pairs(with_cigar=True)

# %% jupyter={"source_hidden": true}
# query_name = "IQEC"
# query_seq = str(primers_dict[query_name])
# target_name = read.query_name
# target_seq = read.query_sequence  # soft clipped seq of the read

# ic(query_name, query_seq)
# ic(target_name, target_seq);

# %%

# %% jupyter={"source_hidden": true}
# aligner = Align.PairwiseAligner(
#     # match_score=1.0, mismatch_score=-2.0, gap_score = -2.5,
#     # mode="local",
#     scoring="blastn"
#     # scoring="megablast"
# )

# global_df = find_best_pairwise_alignments_for_query_and_target(
#     aligner,
#     query_name,
#     target_name,
#     query_seq,
#     target_seq,
#     first_n_alignments_per_strand=10,
#     # debug_mode=True
# )
# global_df

# %% jupyter={"source_hidden": true}
# global_df.loc[
#     (global_df["%QueryIdentity"] == 100)
#     & (global_df["%QueryCoverage"] == 100)
# ].sort_values("Gaps", ascending=False)

# %%
# print(global_df.loc[0, "AlignmentObject"])

# %%
# aligner.mode

# %% jupyter={"source_hidden": true}
# aligner = Align.PairwiseAligner(
#     # match_score=1.0, mismatch_score=-2.0, gap_score = -2.5,
#     mode="local",  # otherwise we'd get scattered matches of the primer across the read
#     scoring="blastn",
#     # scoring="megablast"
# )

# local_df = find_best_pairwise_alignments_for_query_and_target(
#     aligner,
#     query_name,
#     target_name,
#     query_seq,
#     target_seq,
#     first_n_alignments_per_strand=10,
#     debug_mode=True,
# )
# local_df

# %% jupyter={"source_hidden": true}
# local_df.loc[
#     (local_df["%QueryIdentity"] == 100)
#     & (local_df["%QueryCoverage"] == 100)
#     & (local_df["Gaps"] == 0)
# ]

# %% jupyter={"source_hidden": true}
# local_df.loc[
#     (local_df["%QueryIdentity"] == 100)
#     & (local_df["%QueryCoverage"] == 100)
#     & (local_df["Gaps"] == 0),
#     "AlignmentObject",
# ].apply(lambda x: x.aligned)

# %% jupyter={"source_hidden": true}
# for i, alignment in enumerate(positive_strand_alignments, start=1):
#     # ic(i)
#     score = alignment.score
#     gaps, identities, mismatches = alignment.counts()
#     alignment_length = alignment.length
#     print(
#         f"{score = }, {gaps = }, {identities = }, {mismatches = }, {alignment_length = }\n"
#     )
#     print(alignment)
#     if i == 10:
#         break
#     print("\n")

# %% jupyter={"source_hidden": true}
# alignment = alignments[0]

# print(alignment)
# # ic(alignment)

# %% jupyter={"source_hidden": true}
# score = alignment.score
# gaps, identities, mismatches = alignment.counts()
# alignment_length = alignment.length
# ic(score, gaps, identities, mismatches, alignment_length);

# %%
# alignment.coordinates

# %%
# alignment.aligned

# %% jupyter={"source_hidden": true}
# aligned_target_indices, aligned_query_indices = alignment.aligned
# aligned_target_indices, aligned_query_indices

# %%
# target_aligned_indices = alignment.aligned[0]

# target_aligned_indices

# %% jupyter={"source_hidden": true}
# query_aligned_indices = alignment.aligned[1]
# query_aligned_indices

# %% jupyter={"source_hidden": true}
# query_aligned_indices = alignment.aligned[1]
# sum([q_end_i - q_start_i for q_start_i, q_end_i in query_aligned_indices])

# %%

# %%

# %%

# %% jupyter={"source_hidden": true}
# #     length of the aligned query sequence.

# #    This is equal to query_alignment_end - query_alignment_start

# read.query_alignment_length

# %% jupyter={"source_hidden": true}
# """
# infer query length from CIGAR alignment.

# This method deduces the query length from the CIGAR alignment but does not include hard-clipped bases.

# Returns None if CIGAR alignment is not present.

# If always is set to True, infer_read_length is used instead. This is deprecated and only present for backward compatibility.

# """

# read.infer_query_length()

# %% jupyter={"source_hidden": true}
# """
# infer read length from CIGAR alignment.

# This method deduces the read length from the CIGAR alignment including hard-clipped bases.

# Returns None if CIGAR alignment is not present.

# """

# read.infer_read_length()

# %% jupyter={"source_hidden": true}
# # fn: Number of passes in PacBio sequencing (or similar, platform-dependent).
# read.get_tag("fn")

# %% jupyter={"source_hidden": true}
# # np: Number of passes through the template (PacBio circular consensus reads).
# read.get_tag("np")

# %%
# read.query_sequence

# %% jupyter={"source_hidden": true}
# # bq: Base qualities of the barcode sequence.
# read.get_tag("rq")

# %%
# /private7/projects/Combinatorics/D.pealeii/Alignment/UMILongReads/LP2IQEC.r64296e203404D01.aligned.sorted.bam

# %%

# %%
def get_read_and_target_aligned_starts_and_ends(aligned_read_target_pairs):
    # aligned_read_target_pairs = read.get_aligned_pairs(matches_only=True)
    aligned_read_start, aligned_target_start = aligned_read_target_pairs[0]
    aligned_read_end, aligned_target_end = aligned_read_target_pairs[-1]
    aligned_read_end += 1
    aligned_target_end += 1
    return (
        aligned_read_start,
        aligned_read_end,
        aligned_target_start,
        aligned_target_end,
    )


# %%
mapped_bam_dfs = []

for bam_file in mapped_bam_files:

    sample = bam_file.name.split(".")[0]
    gene = sample[3:]
    repeat = sample[2]

    # expected_chrom = chrom_per_gene_dict[gene]

    with pysam.AlignmentFile(
        bam_file,
        "rb",
        threads=10,
        # check_sq=False,
        # require_index=False,
        # index_filename=str(Path(bam_file.parent, f"{bam_file.name}.pbi")),
    ) as samfile:
        reads = [read for read in samfile]
        reads_names = [read.query_name for read in reads]
        read_quality_tags = [read.get_tag("rq") for read in reads]
        mapped_strands = ["+" if read.is_forward else "-" for read in reads]
        mapped_chroms = [read.reference_name for read in reads]

        all_reads_and_targets_aligned_starts_and_ends = [
            get_read_and_target_aligned_starts_and_ends(
                read.get_aligned_pairs(matches_only=True)
            )
            for read in reads
        ]

        df = pd.DataFrame(
            {
                "Sample": sample,
                "Gene": gene,
                "Repeat": repeat,
                "Read": reads_names,
                "ReadQuality": read_quality_tags,
                "MappedStrand": mapped_strands,"MappedChrom": mapped_chroms,
            }
        )

        alignment_boundries_df = pd.DataFrame(
            all_reads_and_targets_aligned_starts_and_ends,
            columns=["ReadStart", "ReadEnd", "GeneStart", "GeneEnd"],
        )

        df = pd.concat([df, alignment_boundries_df], axis=1)

        mapped_bam_dfs.append(df)

concat_mapped_bams_df = pd.concat(mapped_bam_dfs, ignore_index=True)
concat_mapped_bams_df

# %%
concat_mapped_bams_df.groupby("Gene").size()

# %%

# %%
# mapped_bam_file = mapped_bam_files[0]
# mapped_bam_file

# %%
# mapped_bam_files

# %%
# # example of read mapped to the - strand

# read_name = "m64296e_241222_071206/74516384/ccs"

# with pysam.AlignmentFile(
#     mapped_bam_file,
#     "rb",
#     threads=10,
# ) as samfile:
#     read = [read for read in samfile if read.query_name == read_name][0]
# print(read)

# %%
# aligned_read_target_pairs = read.get_aligned_pairs(matches_only=True)
# aligned_read_start, aligned_target_start = aligned_read_target_pairs[0]
# aligned_read_end, aligned_target_end = aligned_read_target_pairs[-1]
# aligned_read_end += 1
# aligned_target_end += 1

# (
#     aligned_read_start,
#     aligned_read_end,
#     aligned_target_start,
#     aligned_target_end,
# )

# %%
# # example of read mapped to the + strand

# read_name = "m64296e_241222_071206/139266104/ccs"

# with pysam.AlignmentFile(
#     mapped_bam_file,
#     "rb",
#     threads=10,
# ) as samfile:
#     read = [read for read in samfile if read.query_name == read_name][0]
# print(read)

# %%
# aligned_read_target_pairs = read.get_aligned_pairs(matches_only=True)
# aligned_read_start, aligned_target_start = aligned_read_target_pairs[0]
# aligned_read_end, aligned_target_end = aligned_read_target_pairs[-1]
# aligned_read_end += 1
# aligned_target_end += 1

# (
#     aligned_read_start,
#     aligned_read_end,
#     aligned_target_start,
#     aligned_target_end,
# )

# %%

# %%

# %% [markdown]
# # Concat BAMs
#

# %%
concat_bams_df = (
    concat_unmapped_bams_df.drop(
        columns=["RawBarcodeSeq", "CorrectedBarcodeSeq", "BarcodeTag"]
    )
    # .rename(columns={"Seq": "UnalignedSeq", "ReadLength": "UnalignedSeqLength"})
    # .rename(columns={"ReadLength": "SeqLength"})
    .merge(
        concat_mapped_bams_df
        # .rename(columns={"Seq": "AlignedSeq", "ReadLength": "AlignedSeqLength"})
        ,
        how="outer",
        indicator="indicator",
    )
)

# make sure all aligned reads also appear in the unaligned reads' bams
assert concat_bams_df["indicator"].value_counts()["right_only"] == 0

concat_bams_df.insert(
    concat_bams_df.columns.get_loc("Read") + 1,
    "Mapped",
    concat_bams_df["indicator"].apply(lambda x: True if x == "both" else False),
)

concat_bams_df["MappedGene"] = concat_bams_df.apply(
    lambda x: (
        gene_by_chrom_dict.get(x["MappedChrom"], "Other") if x["Mapped"] else "Unmapped"
    ),
    axis=1,
)

# after using this col for validating a correct merge,
# and for determining wether a read is edited,
# we don't need it anymore
del concat_bams_df["indicator"]

concat_bams_df.loc[concat_bams_df["MappedStrand"].eq("-"), "Seq"] = concat_bams_df.loc[
    concat_bams_df["MappedStrand"].eq("-"), "Seq"
].apply(reverse_complement)

concat_bams_df.insert(
    concat_bams_df.columns.get_loc("ReadQuality") + 1,
    "HighQualityRead",
    concat_bams_df["ReadQuality"].ge(min_read_quality),
)

concat_bams_df

# %%
concat_bams_df.loc[concat_bams_df["Mapped"]]["Gene"].value_counts()

# %%
mapped_bam_files

# %%
# # example of read mapped to the - strand

# read_name = "m64296e_241222_071206/100007945/ccs"

# with pysam.AlignmentFile(
#     mapped_bam_files[5],  # LP1ADAR1
#     "rb",
#     threads=10,
# ) as samfile:
#     read = [read for read in samfile if read.query_name == read_name][0]
# # print(read)
# concat_bams_df.loc[concat_bams_df["Read"] == read_name]

# %%
# start index of the aligned query portion of the sequence (0-based, inclusive).
# This the index of the first base in query_sequence that is not soft-clipped.

# read.query_alignment_start

# %%
# end index of the aligned query portion of the sequence (0-based, exclusive)
# This the index just past the last base in query_sequence that is not soft-clipped.

# read.query_alignment_end

# %%
# # the way the seq was stored in the original FASTQ file
# original_seq = concat_bams_df.loc[concat_bams_df["Read"] == read_name, "Seq"].values[0]

# # the way the seq was originally stored in the BAM file
# rev_comp_seq = reverse_complement(
#     (concat_bams_df.loc[concat_bams_df["Read"] == read_name, "Seq"].values[0])
# )

# %%
# # the way the seq was stored in the original FASTQ file
# original_seq = reverse_complement(
#     (concat_bams_df.loc[concat_bams_df["Read"] == read_name, "Seq"].values[0])
# )

# # the way the seq is now stored in the BAM file
# rev_comp_seq = concat_bams_df.loc[concat_bams_df["Read"] == read_name, "Seq"].values[0]

# %%
# soft_clipped_seq = read.query_sequence
# # aligned_seq = read.query_alignment_sequence

# soft_clipped_seq == original_seq

# %%
# read.query_sequence == read.get_forward_sequence()

# %%
# read.get_forward_sequence() == rev_comp_seq

# %%
# gene_seq = genes_seq_dict["ADAR1"]

# %%
# aligner = Align.PairwiseAligner(
#     # match_score=1.0, mismatch_score=-2.0, gap_score = -2.5,
#     # mode="global",
#     scoring="blastn"
#     # scoring="megablast"
# )

# alignments = aligner.align(gene_seq, original_seq, strand="+")
# # ic(len(alignments))
# alignment = alignments[0]
# ic(alignment.score)
# print(alignment)

# %%
# aligner = Align.PairwiseAligner(
#     # match_score=1.0, mismatch_score=-2.0, gap_score = -2.5,
#     # mode="global",
#     scoring="blastn"
#     # scoring="megablast"
# )

# alignments = aligner.align(gene_seq, rev_comp_seq, strand="-")
# ic(len(alignments))
# alignment = alignments[0]
# ic(alignment.score)
# print(alignment)

# %%
# aligner = Align.PairwiseAligner(
#     # match_score=1.0, mismatch_score=-2.0, gap_score = -2.5,
#     # mode="global",
#     scoring="blastn"
#     # scoring="megablast"
# )

# alignments = aligner.align(gene_seq, rev_comp_seq, strand="+")
# # ic(len(alignments))
# alignment = alignments[0]
# ic(alignment.score)
# print(alignment)

# %%
# aligner = Align.PairwiseAligner(
#     # match_score=1.0, mismatch_score=-2.0, gap_score = -2.5,
#     # mode="global",
#     scoring="blastn"
#     # scoring="megablast"
# )

# alignments = aligner.align(gene_seq, original_seq, strand="-")
# ic(len(alignments))
# alignment = alignments[0]
# print(alignment)

# %%

# %%

# %%
# aligned_read_target_pairs = read.get_aligned_pairs(matches_only=True)
# aligned_read_start, aligned_target_start = aligned_read_target_pairs[0]
# aligned_read_end, aligned_target_end = aligned_read_target_pairs[-1]
# aligned_read_end += 1
# aligned_target_end += 1

# (
#     aligned_read_start,
#     aligned_read_end,
#     aligned_target_start,
#     aligned_target_end,
# )

# %%
concat_bams_df["Read"].nunique()

# %%
concat_bams_df["MappedGene"].value_counts()

# %%
concat_bams_df.groupby(["Gene", "Repeat"])["MappedGene"].value_counts().reset_index()

# %%
mapped_reads_stats_df = (
    concat_bams_df.groupby(["Gene", "Repeat", "Mapped", "HighQualityRead"])
    .size()
    .reset_index(name="NumOfReads")
)
mapped_reads_stats_df

# %%
concat_bams_df.groupby(
    [
        "Gene",
        "Repeat",
    ]
)["Mapped"].apply(
    lambda x: 100 * x.sum() / x.size
).reset_index(name="%MappedReads/AllReads").round(2)

# %%
mapped_reads_stats_df.groupby(["Gene", "Repeat"]).apply(
    lambda x: 100
    * x.loc[x["Mapped"] & x["HighQualityRead"], "NumOfReads"]
    / x["NumOfReads"].sum(),
    include_groups=False,
).round(2).reset_index(name="%HighQualityMappedReads/AllReads").drop(columns="level_2")

# %%
mapped_reads_stats_df.loc[mapped_reads_stats_df["Mapped"]].groupby(
    ["Gene", "Repeat"]
).apply(
    lambda x: 100 * x.loc[x["HighQualityRead"], "NumOfReads"] / x["NumOfReads"].sum(),
    include_groups=False,
).round(
    2
).reset_index(
    name="%HighQualityMappedReads/AllMappedReads"
).drop(
    columns="level_2"
)

# %%
# concat_bams_df.dropna()["%AlignedSeqLength/RawSeqLength"].value_counts(dropna=False)

# %%
# fig = px.histogram(
#     concat_bams_df.dropna(),
#     # x="UMILength",
#     x="%AlignedSeqLength/RawSeqLength",
#     # color="Gene",
#     facet_col="Gene",
#     color="Repeat",
#     # facet_row="Sample",
#     # log_x=True,
#     log_y=True,
#     histnorm="percent",
#     # cummulative=True
# )
# # Reduce opacity to see both histograms
# fig.update_traces(opacity=0.5)
# # Overlay both histograms
# fig.update_layout(
#     barmode="overlay",
#     template="plotly_white",
#     width=900,
#     height=400,
#     # showlegend=False
# )

# fig.show()

# %%
# fig = px.histogram(
#     concat_bams_df.dropna(),
#     # x="UMILength",
#     x="%AlignedSeqLength/RawSeqLength",
#     # color="Gene",
#     facet_col="Gene",
#     facet_row="HighQualityRead",
#     color="Repeat",
#     # facet_row="Sample",
#     # log_x=True,
#     log_y=True,
#     histnorm="percent",
#     # cummulative=True
# )
# # Reduce opacity to see both histograms
# fig.update_traces(opacity=0.5)
# # Overlay both histograms
# fig.update_layout(
#     barmode="overlay",
#     template="plotly_white",
#     width=900,
#     height=600,
#     # showlegend=False
# )

# fig.show()

# %% [markdown]
# # Final used reads

# %%
merged_old_to_new_reads_dfs = []

for old_to_new_reads_file, gene in zip(merged_old_to_new_reads_files, genes):
    # ic(old_to_new_reads_file)
    old_to_new_reads_df = pd.read_csv(old_to_new_reads_file, sep="\t")
    old_to_new_reads_df.insert(0, "Gene", gene)
    # old_to_new_reads_df = old_to_new_reads_df.merge(
    #     concat_merged_mapped_bam_df,
    #     left_on="OldRead",
    #     right_on="Read",
    #     how="left",
    # ).drop(columns="Read")
    # ic(old_to_new_reads_df.head(3))
    merged_old_to_new_reads_dfs.append(old_to_new_reads_df)

# merged_old_to_new_reads_dfs[0]

concat_merged_old_to_new_reads_df = pd.concat(
    merged_old_to_new_reads_dfs, ignore_index=True
)

# concat_merged_old_to_new_reads_df["Gene"] = concat_merged_old_to_new_reads_df["Gene"].apply(
#     lambda x: "IQEC1" if x == "IQEC" else x
# )

concat_merged_old_to_new_reads_df

# %%
used_reads_df =  pd.concat(
    [
        pd.read_csv(reads_file, sep="\t").loc[:, ["Gene", "Read"]] 
        for reads_file in merged_annotated_reads_files
    ], 
    ignore_index=True
)

used_reads_df["Gene"] = used_reads_df["Gene"].apply(
    lambda x: "IQEC" if x == "IQEC1" else x
)

used_reads_df = used_reads_df.merge(
    concat_merged_old_to_new_reads_df,
    left_on=["Gene", "Read"],
    right_on=["Gene", "NewRead"],
    how="inner",
).drop(columns="Read")

used_reads_df

# %%
used_reads_df["Gene"].value_counts()

# %%
concat_bams_df

# %%
concat_bams_df = concat_bams_df.merge(
    used_reads_df,
    left_on=["Gene", "Read"],
    right_on=["Gene", "OldRead"],
    how="inner",
)

concat_bams_df = concat_bams_df.loc[
    concat_bams_df["MappedGene"].eq(concat_bams_df["Gene"])
]

assert concat_bams_df.shape[0] == concat_bams_df.drop_duplicates(["Gene", "Read"]).shape[0]

concat_bams_df.insert(
    3,
    "NR2",
    concat_bams_df["NewRead"],
)

concat_bams_df = concat_bams_df.drop(columns=["OldRead", "NewRead"])
concat_bams_df = concat_bams_df.rename(columns={"Read": "OldRead"})
concat_bams_df = concat_bams_df.rename(columns={"NR2": "Read"})
concat_bams_df

# %%
concat_bams_df["Gene"].value_counts()

# %%
concat_bams_df.groupby(["Gene", "Repeat"]).size().reset_index(name="NumOfReads")


# %%

# %% [markdown]
# # Find barcodes in reads
#

# %%
def find_first_n_pairwise_alignments_for_barcode_and_read(
    sample,
    gene,
    repeat,
    barcode,  # query_name
    read,  # target_name
    old_read, # old_read is the original read name in the bam file, while read is a shortend version of it
    mapped_gene,  # ADAR1, IQEC, Other, Unmapped
    rtg_strand,  # the genomic strand to which the read was mapped to
    rtg_read_start,
    rtg_read_end,
    rtg_gene_start,
    rtg_gene_end,
    barcode_seq,  # query_seq
    read_seq,  # target_seq
    first_n_alignments_per_strand=10,
    debug_mode=False,
):

    aligner = Align.PairwiseAligner(
        # match_score=1.0, mismatch_score=-2.0, gap_score = -2.5,
        mode="local",  # otherwise we'd get scattered matches of the primer across the read
        scoring="blastn",
        # scoring="megablast"
    )

    read_seq_len = len(read_seq)
    barcode_seq_len = len(barcode_seq)

    STRANDS = ["+", "-"]

    rows = []

    for strand in STRANDS:

        alignments = aligner.align(read_seq, barcode_seq, strand=strand)

        try:
            num_of_alignments = len(alignments)
        except OverflowError:
            num_of_alignments = np.inf
        if debug_mode:
            print(f"{num_of_alignments = }")

        for i, alignment in enumerate(alignments, start=1):

            score = alignment.score
            gaps, identities, mismatches = alignment.counts()
            alignment_length = alignment.length

            # prct_barcode_identity = 100 * identities / len(barcode_seq) # todo this is a bit stringent, as it takes len(barcode_seq) instead of alignment_length
            prct_barcode_identity = 100 * identities / alignment_length
            prct_barcode_cov = 100 * (identities + mismatches) / len(barcode_seq)

            btr_read_coords, btr_barcode_coords = alignment.aligned

            num_of_btr_read_coords = len(btr_read_coords)
            num_of_btr_barcode_coords = len(btr_barcode_coords)
            num_of_btr_barcode_gap_openings = num_of_btr_barcode_coords - 1

            if debug_mode:
                # print(f"aligner = {aligner}")
                print(
                    f"{score = }, {gaps = }, {identities = }, {mismatches = }, {alignment_length = }, {prct_barcode_identity = :.2f}, {prct_barcode_cov = :.2f}\n"
                )
                print(alignment)
                print("\n")

            row = [
                sample,  # Sample"
                gene,  # Gene
                repeat,  # Repeat
                read,  # Read
                old_read, # OldRead
                barcode,  # Barcode
                mapped_gene,  # MappedGene
                read_seq,  # TargetSeq
                rtg_strand,  # RTGStrand - the genomic strand to which the read was mapped to
                rtg_read_start,
                rtg_read_end,
                rtg_gene_start,
                rtg_gene_end,
                strand,  # BTRStrand
                num_of_alignments,  # BTRNumOfAllAlignmentsPerStrand
                i,  # BTR#AlignmentPerStrand
                score,  # BTRScore
                gaps,  # BTRGaps
                identities,  # BTRIdentitites
                mismatches,  # BTRMismatches
                prct_barcode_identity,  # %BTRBarcodeIdentity
                prct_barcode_cov,  # %BTRBarcodeCoverage
                read_seq_len,  # ReadSeqLength
                barcode_seq_len,  # BarcodeSeqLength
                alignment_length,  # AlignmentLength
                btr_read_coords,  # BTRReadCoords
                btr_barcode_coords,  # BTRBarcodeCoords
                num_of_btr_read_coords,  # "NumOfBTRReadCoords",
                num_of_btr_barcode_coords,  # "NumOfBTRBarcodeCoords",
                num_of_btr_barcode_gap_openings,  # "NumOfBTRBarcodeGapOpenings",
                alignment,  # BTRAlignmentObject
            ]
            rows.append(row)

            if (
                first_n_alignments_per_strand is not None
                and i == first_n_alignments_per_strand
            ):
                break

    df = pd.DataFrame(
        rows,
        columns=[
            "Sample",
            "Gene",
            "Repeat",
            "Read",  # "TargetName",
            "OldRead",
            "Barcode",  # "QueryName",
            "MappedGene",
            "ReadSeq",  # "TargetSeq",
            "RTGStrand",  # "GenomicTargetStrand",
            "RTGReadStart",
            "RTGReadEnd",
            "RTGGeneStart",
            "RTGGeneEnd",
            "BTRStrand",  # "TargetStrand",
            "BTRNumOfAllAlignmentsPerStrand",
            "BTR#AlignmentPerStrand",
            "BTRScore",
            "BTRGaps",
            "BTRIdentitites",
            "BTRMismatches",
            "%BTRBarcodeIdentity",  # "%QueryIdentity",
            "%BTRBarcodeCoverage",  # "%QueryCoverage",
            "ReadSeqLength",  # "TargetSeqLength",
            "BarcodeSeqLength",  # "QuerySeqLength",
            "BTRAlignmentLength",  # "AlignmentLength",
            "BTRReadCoords",  # "AlignedTargetCoords",
            "BTRBarcodeCoords",  # "AlignedQueryCoords",
            "NumOfBTRReadCoords",
            "NumOfBTRBarcodeCoords",
            "NumOfBTRBarcodeGapOpenings",
            "BTRAlignmentObject",  # "AlignmentObject",
        ],
    )

    return df


# %%
def align_reads_to_barcode_one_batch(
    one_gene_bams_df,
    # aligner,
    gene,
    barcode,  # primer name
    barcode_seq,  # primer seq
    first_n_alignments_per_strand,
    debug_mode,
):
    """Find the first n alignments per strand between the barcode (query) and each read in the batch."""
    alignments_dfs = one_gene_bams_df.apply(
        lambda x: find_first_n_pairwise_alignments_for_barcode_and_read(
            x["Sample"],
            gene,
            x["Repeat"],
            # aligner,
            barcode,
            x["Read"],  # target_name
            x["OldRead"],  # old_read
            x["MappedGene"],  # mapped_gene
            x[
                "MappedStrand"
            ],  # genomic_target_strand - the genomic strand to which the read was mapped to
            x["ReadStart"],
            x["ReadEnd"],
            x["GeneStart"],
            x["GeneEnd"],
            barcode_seq,
            x["Seq"],  # read_seq
            first_n_alignments_per_strand=first_n_alignments_per_strand,
            debug_mode=debug_mode,
        ),
        axis=1,
    )
    alignments_df = pd.concat([df for df in alignments_dfs])

    return alignments_df


# %%
def align_reads_to_barcode_by_parallel_batches(
    processes,
    one_gene_bams_df,
    gene,
    barcode,
    barcode_seq,
    first_n_alignments_per_strand=10,
    debug_mode=False,
    batches_per_process=10,
):
    n = processes * batches_per_process

    # Compute split indices
    split_sizes = [
        len(one_gene_bams_df) // n + (i < len(one_gene_bams_df) % n) for i in range(n)
    ]
    split_indices = [0] + list(pd.Series(split_sizes).cumsum())

    # Split DataFrame
    split_dfs = [
        one_gene_bams_df.iloc[split_indices[i] : split_indices[i + 1]] for i in range(n)
    ]

    assert sum(len(df) for df in split_dfs) == len(one_gene_bams_df)
    assert pd.concat(split_dfs).equals(one_gene_bams_df)

    with Pool(processes=processes) as pool:
        alignments_dfs = pool.starmap(
            func=align_reads_to_barcode_one_batch,
            iterable=[
                (
                    df,
                    gene,
                    barcode,
                    barcode_seq,
                    first_n_alignments_per_strand,
                    debug_mode,
                )
                for df in split_dfs
            ],
        )
    for i, df in enumerate(alignments_dfs):
        if len(df) == 0:
            print(f"Warning: alignments_df for batch {i} is empty.")
    concat_alignments_df = pd.concat(alignments_dfs, ignore_index=True)
    return concat_alignments_df


# %%
def get_genomic_coord_for_read_coord(
    aligned_pairs,
    required_read_pos,
    required_read_pos_is_exclusive_end_pos_of_some_feature_on_the_read=False,
):
    if required_read_pos_is_exclusive_end_pos_of_some_feature_on_the_read:
        required_read_pos -= 1
    for read_pos, ref_pos, _ in aligned_pairs:  # _ is for CIGAR operation
        if read_pos == required_read_pos:
            if (
                required_read_pos_is_exclusive_end_pos_of_some_feature_on_the_read
                and ref_pos is not None
            ):
                ref_pos += 1  # convert back to inclusive end position
            return ref_pos
    return None


# %%

def find_btg_gene_coords_one_sample(input_df, mapped_bam_file, threads=10):

    # required_reads = input_df["Read"].unique().tolist()
    unique_required_old_reads_names = input_df["OldRead"].unique().tolist()

    with pysam.AlignmentFile(
        mapped_bam_file,
        "rb",
        threads=threads,
    ) as samfile:
        # get all reads in this specific bam file that are in the input_df
        # reads = [read for read in samfile if read.query_name in required_reads]
        reads = [read for read in samfile if read.query_name in unique_required_old_reads_names]

    # names of reads found in this specific bam file,
    # in the order they appear in the bam file
    # reads_names = [read.query_name for read in reads]
    old_reads_names = [read.query_name for read in reads]
    # if len(reads_names) == 0:
    #     continue
    aligned_pairs = [
        read.get_aligned_pairs(matches_only=False, with_cigar=True)
        for read in reads
    ]
    # aligned_pairs_series = pd.Series(aligned_pairs, index=reads_names)
    aligned_pairs_series = pd.Series(aligned_pairs, index=old_reads_names)

    # annotated_df = input_df.loc[input_df["Read"].isin(reads_names)].copy()
    annotated_df = input_df.loc[input_df["OldRead"].isin(old_reads_names)].copy()

    # annotated_df["BTGGeneStart"] = annotated_df.apply(
    #     lambda x: get_genomic_coord_for_read_coord(
    #         aligned_pairs_series.loc[x["Read"]], x["BTRReadStart"]
    #     ),
    #     axis=1,
    # )
    # annotated_df["BTGGeneEnd"] = annotated_df.apply(
    #     lambda x: get_genomic_coord_for_read_coord(
    #         aligned_pairs_series.loc[x["Read"]], x["BTRReadEnd"], True
    #     ),
    #     axis=1,
    # )
    annotated_df["BTGGeneStart"] = annotated_df.apply(
        lambda x: get_genomic_coord_for_read_coord(
            aligned_pairs_series.loc[x["OldRead"]], x["BTRReadStart"]
        ),
        axis=1,
    )
    annotated_df["BTGGeneEnd"] = annotated_df.apply(
        lambda x: get_genomic_coord_for_read_coord(
            aligned_pairs_series.loc[x["OldRead"]], x["BTRReadEnd"], True
        ),
        axis=1,
    )

    return annotated_df


# %%
def parallel_find_btg_gene_coords(input_df, mapped_bam_files, processes=6, threads=10):
    """Find `BTGGeneStart` and `BTGGeneEnd` for each read in the input_df."""
    samples = [bam_file.name.split(".")[0] for bam_file in mapped_bam_files]
    # input_df = input_df.loc[:, ["Sample", "Read", "BTRReadStart", "BTRReadEnd"]].copy()
    input_dfs = [input_df.loc[input_df["Sample"].eq(sample)] for sample in samples]

    with Pool(processes=processes) as pool:
        annotated_dfs = pool.starmap(
            func=find_btg_gene_coords_one_sample,
            iterable=zip(input_dfs, mapped_bam_files, itertools.repeat(threads)),
        )
    concat_annotated_df = pd.concat(annotated_dfs)
    return concat_annotated_df


# %%
processes = 10
first_n_alignments_per_strand = 100
debug_mode = False

alignments_dfs = []

for gene in genes:

    # one_gene_bams_df = concat_bams_df.loc[
    #     (concat_bams_df["Gene"].eq(gene)) & (concat_bams_df["Mapped"])
    # ]

    # consider all reads from a certain gene, whether they were mapped or not,
    # and whether they were mapped to the correct gene or not
    one_gene_bams_df = concat_bams_df.loc[(concat_bams_df["Gene"].eq(gene))]

    # primers = [gene, "PCR"]

    # for primer in primers:

    # align each primer against each read, wether we expect that read to contain that primer/barcode or not
    for barcode in primers_dict:

        ic(gene, barcode)
        
        barcode_seq = primers_dict[barcode]

        alignments_df = align_reads_to_barcode_by_parallel_batches(
            processes,
            one_gene_bams_df,
            gene,
            barcode,  # primer name
            barcode_seq,  # primer seq
            first_n_alignments_per_strand=first_n_alignments_per_strand,
            debug_mode=debug_mode,
        )

        alignments_dfs.append(alignments_df)

        # break # break from 2nd primer

    # break # break from 2nd gene

concat_alignments_df = pd.concat(alignments_dfs, ignore_index=True)

concat_alignments_df.insert(
    concat_alignments_df.columns.get_loc("BTRReadCoords") + 1,
    "BTRReadStart",
    concat_alignments_df["BTRReadCoords"].apply(lambda x: x[0][0]),
)
concat_alignments_df.insert(
    concat_alignments_df.columns.get_loc("BTRReadCoords") + 2,
    "BTRReadEnd",
    concat_alignments_df["BTRReadCoords"].apply(lambda x: x[-1][-1]),
)

# # get BTGeneStart and BTGGeneEnd for each read
# concat_alignments_df = find_btg_gene_coords(concat_alignments_df)

concat_alignments_df

# %%
btg_gene_coords_input_df = concat_alignments_df.loc[
    :, ["Sample", "Read", "OldRead", "BTRReadStart", "BTRReadEnd"]
]
btg_gene_coords_input_df

# %%
concat_annotated_btg_gene_coords_df = parallel_find_btg_gene_coords(
    btg_gene_coords_input_df, mapped_bam_files
)
concat_annotated_btg_gene_coords_df

# %%
concat_alignments_df = concat_alignments_df.drop(
    columns=["BTGGeneStart", "BTGGeneEnd"], errors="ignore"
).join(
    concat_annotated_btg_gene_coords_df.drop(
        columns=["Sample", "Read", "OldRead", "BTRReadStart", "BTRReadEnd"]
    )
)
concat_alignments_df

# %%
# concat_alignments_df["Read"].value_counts().describe()

# %%
# # the local aligner outputed (at most?) only one pairwise alignment for each read-primer couple, one for each strand
# assert concat_alignments_df[
#     "NumOfAllAlignmentsPerStrand"
# ].value_counts().index.values == np.array([1])

# %%
# # the local aligner outputed (at most?) only one pairwise alignment for each read-primer couple, one for each strand
# assert concat_alignments_df[
#     "NumOfAllAlignmentsPerStrand"
# ].value_counts().index.values == np.array([1])

# # count how many times each read is present in the read-primer alignment table -
# # if each read-primer couple is present once for each strand, than each read should show up exactly 4 times,
# # one for the PCR primer and another for each of the gene-specific primers, each of which for each of the strands
# assert concat_alignments_df.shape[0] / concat_bams_df.shape[0] == 6
# # since these two tests passed successfuly, we may safely delete the following columns from the df for clarity purposes
# concat_alignments_df = concat_alignments_df.drop(
#     columns=["NumOfAllAlignmentsPerStrand", "#AlignmentPerStrand"]
# )

# concat_alignments_df["NumOfAlignedTargetCoords"] = concat_alignments_df[
#     "AlignedTargetCoords"
# ].apply(len)
# concat_alignments_df["NumOfAlignedQueryCoords"] = concat_alignments_df[
#     "AlignedQueryCoords"
# ].apply(len)
# concat_alignments_df["QueryGapOpenings"] = (
#     concat_alignments_df["NumOfAlignedQueryCoords"] - 1
# )
# assert concat_alignments_df.loc[
#     (concat_alignments_df["Gaps"] == 0) & (concat_alignments_df["QueryGapOpenings"] > 0)
# ].empty

# concat_alignments_df["NumOfBTRReadCoords"] = concat_alignments_df[
#     "BTRReadCoords"
# ].apply(len)
# concat_alignments_df["NumOfBTRBarcodeCoords"] = concat_alignments_df[
#     "BTRBarcodeCoords"
# ].apply(len)
# concat_alignments_df["BTRBarcodeGapOpenings"] = (
#     concat_alignments_df["NumOfBTRBarcodeCoords"] - 1
# )
assert concat_alignments_df.loc[
    (concat_alignments_df["BTRGaps"] == 0)
    & (concat_alignments_df["NumOfBTRBarcodeGapOpenings"] > 0)
].empty

# concat_alignments_df

# %%
# read = "m64296e_241222_071206/10027148/ccs"

# row = concat_alignments_df.loc[
#     (concat_alignments_df["Read"] == read)
#     & (concat_alignments_df["RTGStrand"] == "-")
#     & (concat_alignments_df["Barcode"].eq(concat_alignments_df["Gene"]))
#     & (concat_alignments_df["BTRStrand"] == "-")
# ].iloc[0]
# row

# %%
# alignment = row["BTRAlignmentObject"]

# print(alignment)

# %%
# row["ReadSeq"][3378:3402]

# %%
# barcode_seq = primers_dict[row["Barcode"]]
# # barcode_seq

# score = alignment.score
# gaps, identities, mismatches = alignment.counts()
# alignment_length = alignment.length

# prct_barcode_identity = 100 * identities / len(barcode_seq)
# prct_barcode_cov = 100 * (identities + mismatches) / len(barcode_seq)

# ic(score, gaps, identities, mismatches, alignment_length)
# ic(prct_barcode_identity, prct_barcode_cov)

# %%
# read_object = get_read_object_from_bam_files(read, mapped_bam_files)
# read_object

# %%
# btr_read_coords, btr_barcode_coords = alignment.aligned
# # btr_read_coords

# btr_read_start = btr_read_coords[0][0]
# btr_read_end = btr_read_coords[-1][-1]
# btr_read_start, btr_read_end

# %%
# btr_read_end - btr_read_start

# %%
# aligned_pairs_df = pd.DataFrame(
#     read_object.get_aligned_pairs(matches_only=False, with_cigar=True),
#     columns=["ReadPos", "RefPos", "Cigar"],
# )
# aligned_pairs_df

# %%
# aligned_pairs_df.loc[
#     (aligned_pairs_df["ReadPos"] >= btr_read_start)
#     & (aligned_pairs_df["ReadPos"] < btr_read_end)
# ]

# %%
# expected_start, expected_end = expected_barcode_locations_on_gene_dict[row["Gene"]]
# ic(expected_start, expected_end, expected_end - expected_start)

# %%
# observed_start = aligned_pairs_df.loc[
#     aligned_pairs_df["ReadPos"].eq(btr_read_start), "RefPos"
# ].values[0]
# observed_end = aligned_pairs_df.loc[
#     aligned_pairs_df["ReadPos"].eq(btr_read_end - 1), "RefPos"
# ].values[0]

# ic(observed_start, observed_end, observed_end - observed_start)

# %%
# coverage = min(observed_end + 1, expected_end) - max(observed_start, expected_start)
# coverage

# %%

# %% [markdown]
# # Per-position coverage

# %%
# coverage_depth_files = [
#     Path(mapped_bams_dir, f"{gene}.CoverageDepth.tsv") for gene in genes
# ]

# %%
# for (gene, chrom), coverage_depth_file in zip(
#     chrom_per_gene_dict.items(), coverage_depth_files
# ):
#     gene_mapped_bam_files = [file for file in mapped_bam_files if gene in file.name]
#     gene_mapped_bam_files = " ".join(str(file) for file in gene_mapped_bam_files)
#     cmd = f"samtools depth -a -H --min-BQ 30 -r {chrom} -o {coverage_depth_file} {gene_mapped_bam_files}"
#     print(cmd)
#     subprocess.run(cmd, shell=True, check=True)

# %%

# %%
coverage_depth_files = [
    Path(mapped_merged_filtered_bams_dir, f"{gene}.FinalCoverageDepth.tsv") for gene in genes
]

# %%
# filter the bam files according to the final reads used
temp_dir = Path(mapped_merged_filtered_bams_dir, "temp")
temp_dir.mkdir(exist_ok=True)

mapped_merged_filtered_bams = []

for in_bam_path, gene in zip(mapped_merged_bam_files, genes):
    
    # Get the set of unique read names to keep
    reads_to_keep = set(
        concat_bams_df.loc[
            concat_bams_df["Gene"] == gene,
            "OldRead"
        ]
    )
    
    out_bam_path = Path(temp_dir, in_bam_path.name)
    
    mapped_merged_filtered_bams.append(out_bam_path)

    with pysam.AlignmentFile(in_bam_path, "rb") as in_bam, pysam.AlignmentFile(
        out_bam_path, "wb", template=in_bam
    ) as out_bam:
        for read in in_bam:
            if read.query_name in reads_to_keep:
                out_bam.write(read)

    print(f"Filtered BAM written to: {out_bam_path}")
    
# !samtools index -M {temp_dir}/*.bam


# create coverage depth files using the filtered bam files
for (gene, chrom), coverage_depth_file, final_filtered_bam in zip(
    chrom_per_gene_dict.items(), coverage_depth_files, mapped_merged_filtered_bams
):
    # gene_mapped_bam_files = [file for file in mapped_bam_files if gene in file.name]
    # gene_mapped_bam_files = " ".join(str(file) for file in gene_mapped_bam_files)
    # cmd = f"samtools depth -a -H --min-BQ 30 -r {chrom} -o {coverage_depth_file} {gene_mapped_bam_files}"
    cmd = f"samtools depth -a -H --min-BQ 30 -r {chrom} -o {coverage_depth_file} {final_filtered_bam}"
    print(cmd)
    subprocess.run(cmd, shell=True, check=True)
    
    
# clean up the temporary directory
# !rm -rf {temp_dir}

# %%

# %%
gene_cov_dfs = []

for (gene, chrom), coverage_depth_file in zip(
    chrom_per_gene_dict.items(), coverage_depth_files
):
    with open(coverage_depth_file) as f:
        header_line = f.readline().lstrip("#").strip()
        colnames = header_line.split("\t")
        colnames = ["Chrom", "Position"] + [
            colname.split("/")[-1].split(".")[0] for colname in colnames[2:]
        ]

    gene_cov_df = pd.read_table(
        coverage_depth_file,
        comment=None,  # don't ignore lines starting with #
        names=colnames,
        skiprows=1,  # skip the header line as data
    )
    gene_cov_df["Position"] = gene_cov_df["Position"] - 1
    gene_cov_df = gene_cov_df.melt(
        id_vars=["Chrom", "Position"],
        var_name="Sample",
        value_name="Coverage",
    )
    gene_cov_df = gene_cov_df.groupby("Position")["Coverage"].sum().reset_index()
    gene_cov_df = gene_cov_df.rename(columns={"Coverage": "Reads"})

    max_mapped_reads_per_gene = concat_bams_df.loc[
        (concat_bams_df["Gene"].eq(gene))
        & (concat_bams_df["Gene"].eq(concat_bams_df["MappedGene"]))
        & (concat_bams_df["Mapped"])
    ].shape[0]

    gene_cov_df["%OfAllMappedReads"] = (
        100 * gene_cov_df["Reads"] / max_mapped_reads_per_gene
    ).round(2)

    gene_cov_dfs.append(gene_cov_df)

gene_cov_dfs[0]

# %%
gene_cov_dfs[1]

# %%
for gene, gene_cov_df in zip(genes, gene_cov_dfs):
    fig = px.area(
        gene_cov_df,
        x="Position",
        y="Reads",
        title=gene,
        labels={"Reads": "Mapped reads"},
    )
    max_mapped_reads_per_gene = concat_bams_df.loc[
        (concat_bams_df["Gene"].eq(gene))
        & (concat_bams_df["Gene"].eq(concat_bams_df["MappedGene"]))
        & (concat_bams_df["Mapped"])
    ].shape[0]
    last_position_per_gene = gene_cov_df["Position"].max()
    fig.add_shape(
        type="line",
        x0=0,
        x1=last_position_per_gene,
        y0=max_mapped_reads_per_gene,
        y1=max_mapped_reads_per_gene,
        line=dict(
            color="LightSeaGreen",
            width=4,
            dash="dashdot",
        ),
    )
    expected_barcode_start, expected_barcode_end = (
        expected_barcode_locations_on_gene_dict[gene]
    )
    fig.add_shape(
        type="line",
        x0=expected_barcode_start,
        x1=expected_barcode_start,
        y0=0,
        y1=max_mapped_reads_per_gene,
        line=dict(
            color="red",
            width=2,
            dash="dash",
        ),
    )
    fig.add_shape(
        type="line",
        x0=expected_barcode_end - 1,
        x1=expected_barcode_end - 1,
        y0=0,
        y1=max_mapped_reads_per_gene,
        line=dict(
            color="red",
            width=2,
            dash="dash",
        ),
    )
    fig.update_xaxes(dtick=250)
    fig.update_yaxes(dtick=25_000)
    fig.update_layout(
        width=1200,
        height=500,
    )
    fig.show()

# %%
for gene, gene_cov_df in zip(genes, gene_cov_dfs):
    fig = px.area(
        gene_cov_df,
        x="Position",
        y="%OfAllMappedReads",
        labels={"%OfAllMappedReads": "Mapped reads [%]"},
        title=gene,
    )
    expected_barcode_start, expected_barcode_end = (
        expected_barcode_locations_on_gene_dict[gene]
    )
    fig.add_shape(
        type="line",
        x0=expected_barcode_start,
        x1=expected_barcode_start,
        y0=0,
        y1=100,
        line=dict(
            color="red",
            width=2,
            dash="dash",
        ),
    )
    fig.add_shape(
        type="line",
        x0=expected_barcode_end - 1,
        x1=expected_barcode_end - 1,
        y0=0,
        y1=100,
        line=dict(
            color="red",
            width=2,
            dash="dash",
        ),
    )
    fig.update_xaxes(dtick=250)
    fig.update_yaxes(dtick=10, range=[0, 101])
    fig.update_layout(
        width=1200,
        height=500,
    )
    fig.show()

# %%
main_mapping_boundaries_per_gene = []

# for gene_cov_df in gene_cov_dfs:
#     main_mapping_boundaries = (
#         gene_cov_df.loc[gene_cov_df["Reads"].gt(25_000), "Position"]
#         .agg(["min", "max"])
#         .values
#     )
#     main_mapping_boundaries_per_gene.append(main_mapping_boundaries)

for gene in genes:

    mapping_boundaries_stats_df = (
        concat_bams_df.loc[
            (concat_bams_df["Gene"].eq(gene))
            & (concat_bams_df["Gene"].eq(concat_bams_df["MappedGene"]))
            & (concat_bams_df["Mapped"]),
            ["GeneStart", "GeneEnd"],
        ]
        .describe()
        .T
    )

    mean_gene_start = mapping_boundaries_stats_df.at["GeneStart", "25%"]
    mean_gene_end = mapping_boundaries_stats_df.at["GeneEnd", "75%"]
    main_mapping_boundaries = (mean_gene_start, mean_gene_end)
    main_mapping_boundaries_per_gene.append(main_mapping_boundaries)


main_mapping_boundaries_per_gene

# %%
for gene, main_mapping_boundaries in zip(genes, main_mapping_boundaries_per_gene):
    main_mapping_boundary_start, main_mapping_boundary_end = main_mapping_boundaries

    df = concat_bams_df.loc[
        (concat_bams_df["Gene"].eq(gene))
        & (concat_bams_df["Gene"].eq(concat_bams_df["MappedGene"]))
        & (concat_bams_df["Mapped"])
    ].copy()

    df["ReadAlignmentLength"] = df["ReadEnd"] - df["ReadStart"]
    df["GeneAlignmentLength"] = df["GeneEnd"] - df["GeneStart"]

    max_mapped_reads_per_gene = df.shape[0]

    # reads_spanning_from_main_mapping_boundaries_from_start_to_end = df.loc[
    #     (df["GeneStart"].le(main_mapping_boundary_start))
    #     & (df["GeneEnd"].sub(1).ge(main_mapping_boundary_end))
    # ].shape[0]

    # reads_spanning_from_main_mapping_boundaries_from_start_to_end = df.loc[
    #     (df["GeneStart"].le(main_mapping_boundary_start))
    #     & (df["GeneEnd"].ge(main_mapping_boundary_end))
    # ].shape[0]

    df["WithinMainMappingBoundaries"] = df["GeneStart"].le(
        main_mapping_boundary_start
    ) & df["GeneEnd"].ge(main_mapping_boundary_end)
    reads_spanning_from_main_mapping_boundaries_from_start_to_end = df.loc[
        df["WithinMainMappingBoundaries"]
    ].shape[0]

    prct_reads_spanning_from_main_mapping_boundaries_from_start_to_end = np.round(
        100
        * reads_spanning_from_main_mapping_boundaries_from_start_to_end
        / max_mapped_reads_per_gene,
        2,
    )

    ic(
        gene,
        # main_mapping_boundary_start,
        # main_mapping_boundary_end,
        main_mapping_boundaries,
        max_mapped_reads_per_gene,
        reads_spanning_from_main_mapping_boundaries_from_start_to_end,
        prct_reads_spanning_from_main_mapping_boundaries_from_start_to_end,
    )

    # break

# %%
# df

# %%
# fig = px.scatter_matrix(
#     df,
#     dimensions=["ReadLength", "ReadAlignmentLength", "GeneAlignmentLength"],
#     color="WithinMainMappingBoundaries",
#     # title=f"{gene} - Read vs Gene Alignment Lengths",
#     title=gene,
#     labels={
#         "ReadLength": "Read<br>length",
#         "ReadAlignmentLength": "Read<br>alignment<br>length",
#         "GeneAlignmentLength": "Gene<br>alignment<br>length",
#         "WithinMainMappingBoundaries": "Within<br>main<br>mapping<br>boundaries",
#     },
# )
# fig.update_traces(diagonal_visible=False)
# fig.update_layout(width=800, height=600)
# fig.show()

# %%
mapping_stats_df = concat_bams_df.loc[
    # (concat_bams_df["Gene"].eq(gene))
    # &
    (concat_bams_df["Gene"].eq(concat_bams_df["MappedGene"]))
    & (concat_bams_df["Mapped"])
].copy()
mapping_stats_df["ReadAlignmentLength"] = (
    mapping_stats_df["ReadEnd"] - mapping_stats_df["ReadStart"]
)
mapping_stats_df

# %%
mapping_stats_df.groupby("Gene")[
    ["ReadLength", "ReadAlignmentLength", "GeneStart", "GeneEnd"]
].describe().round(2)

# %%
mapping_stats_df.loc[:, ["Gene", "ReadLength", "ReadAlignmentLength"]].groupby(
    "Gene"
).describe().round(2).T

# %%
fig = px.histogram(
    mapping_stats_df,
    x="ReadLength",
    # y="ReadAlignmentLength",
    # histfunc="avg",
    color="Gene",
    histnorm="percent",
    cumulative=True,
    facet_col="Gene",
    # facet_col_spacing=0.05,
    labels={
        "ReadLength": "Read length [bp]",
        # "ReadAlignmentLength": "Read alignment length [bp]",
        # "MappedGene": "Mapped gene",
    },
    # trendline="ols",
    opacity=0.75,
)

fig.update_layout(
    width=1200,
    height=500,
    barmode="overlay",
    # title="Read alignment length vs. read length",
)
fig.show()

# %%

# %% [markdown]
# # Choose best gene-specific barcode alignment

# %% [markdown]
# ## Tests before construction

# %%
# concat_alignments_df

# %%
# test_concat_alignments_df = pd.concat(
#     [
#         concat_alignments_df.loc[
#             (concat_alignments_df["Barcode"] != "PCR")
#             & (concat_alignments_df["Gene"].eq(concat_alignments_df["MappedGene"]))
#             & (concat_alignments_df["Gene"].eq(concat_alignments_df["Barcode"]))
#             & (concat_alignments_df["BTRGaps"].gt(0))
#         ].sample(n=3, random_state=seed),
#         concat_alignments_df.loc[
#             (concat_alignments_df["Barcode"] != "PCR")
#             & (concat_alignments_df["Gene"].eq(concat_alignments_df["MappedGene"]))
#             & (concat_alignments_df["Gene"].eq(concat_alignments_df["Barcode"]))
#             & (concat_alignments_df["BTRGaps"].eq(0))
#         ].sample(n=3, random_state=seed),
#     ],
#     ignore_index=True,
# )

# test_concat_alignments_df

# %%
# test_concat_alignments_df = pd.concat(
#     [
#         concat_alignments_df.loc[
#             (concat_alignments_df["Barcode"] != "PCR")
#             & (concat_alignments_df["Gene"].eq(concat_alignments_df["MappedGene"]))
#             & (concat_alignments_df["Gene"].eq(concat_alignments_df["Barcode"]))
#             & (concat_alignments_df["BTRGaps"].gt(0))
#         ].sample(n=3, random_state=seed),
#         concat_alignments_df.loc[
#             (concat_alignments_df["Barcode"] != "PCR")
#             & (concat_alignments_df["Gene"].eq(concat_alignments_df["MappedGene"]))
#             & (concat_alignments_df["Gene"].eq(concat_alignments_df["Barcode"]))
#             & (concat_alignments_df["BTRGaps"].eq(0))
#         ].sample(n=3, random_state=seed),
#     ],
#     ignore_index=True,
# )

# test_concat_alignments_df.insert(
#     test_concat_alignments_df.columns.get_loc("BTRReadCoords") + 1,
#     "BTRReadStart",
#     test_concat_alignments_df["BTRReadCoords"].apply(lambda x: x[0][0]),
# )
# test_concat_alignments_df.insert(
#     test_concat_alignments_df.columns.get_loc("BTRReadCoords") + 2,
#     "BTRReadEnd",
#     test_concat_alignments_df["BTRReadCoords"].apply(lambda x: x[-1][-1]),
# )

# test_concat_alignments_df

# %%
# def get_genomic_coord_for_read_coord(aligned_pairs, required_read_pos):
#     for read_pos, ref_pos, cigar_op in aligned_pairs:
#         if read_pos == required_read_pos:
#             return ref_pos
#     return None

# %%
# input_df = test_concat_alignments_df


# annotated_dfs = []

# all_required_reads = input_df["Read"].values

# for bam_file in mapped_bam_files:

#     with pysam.AlignmentFile(
#         bam_file,
#         "rb",
#         threads=10,
#     ) as samfile:
#         reads = [read for read in samfile if read.query_name in all_required_reads]
#         reads_names = [
#             read.query_name for read in reads
#         ]  # names of reads found in this specific bam file
#         if len(reads_names) == 0:
#             continue
#         aligned_pairs = [
#             read.get_aligned_pairs(matches_only=False, with_cigar=True)
#             for read in reads
#         ]
#         annotated_df = input_df.loc[input_df["Read"].isin(reads_names)].copy()
#         aligned_pairs_series = pd.Series(aligned_pairs, index=reads_names)
#         annotated_df["BTGGeneStart"] = annotated_df.apply(
#             lambda x: get_genomic_coord_for_read_coord(
#                 aligned_pairs_series.loc[x["Read"]], x["BTRReadStart"]
#             ),
#             axis=1,
#         )
#         annotated_df["BTGGeneEnd"] = annotated_df.apply(
#             lambda x: get_genomic_coord_for_read_coord(
#                 aligned_pairs_series.loc[x["Read"]], x["BTRReadEnd"]
#             ),
#             axis=1,
#         )
#         annotated_dfs.append(annotated_df)

# concat_annotated_df = pd.concat(annotated_dfs)
# concat_annotated_df

# %%
# annotated_df

# %%
# aligned_pairs_series

# %%
# df = pd.DataFrame(
#     aligned_pairs_series.loc[annotated_df["Read"].values[0]],
#     columns=["ReadPos", "RefPos", "CigarOp"],
# ).set_index("ReadPos")
# df.loc[3371 - 3 :]

# %%
# df = pd.DataFrame(
#     aligned_pairs_series.loc[annotated_df["Read"].values[0]],
#     columns=["ReadPos", "RefPos", "CigarOp"],
# ).set_index("ReadPos")
# df.loc[3371 - 3 :]["CigarOp"].value_counts()

# %%
# def find_btg_gene_coords(input_df):

#     annotated_dfs = []

#     all_required_reads = input_df["Read"].values

#     for bam_file in mapped_bam_files:

#         with pysam.AlignmentFile(
#             bam_file,
#             "rb",
#             threads=10,
#         ) as samfile:
#             reads = [read for read in samfile if read.query_name in all_required_reads]
#             reads_names = [
#                 read.query_name for read in reads
#             ]  # names of reads found in this specific bam file
#             if len(reads_names) == 0:
#                 continue
#             aligned_pairs = [
#                 read.get_aligned_pairs(matches_only=False, with_cigar=True)
#                 for read in reads
#             ]
#             annotated_df = input_df.loc[input_df["Read"].isin(reads_names)].copy()
#             aligned_pairs_series = pd.Series(aligned_pairs, index=reads_names)
#             annotated_df["BTGGeneStart"] = annotated_df.apply(
#                 lambda x: get_genomic_coord_for_read_coord(
#                     aligned_pairs_series.loc[x["Read"]], x["BTRReadStart"]
#                 ),
#                 axis=1,
#             )
#             annotated_df["BTGGeneEnd"] = annotated_df.apply(
#                 lambda x: get_genomic_coord_for_read_coord(
#                     aligned_pairs_series.loc[x["Read"]], x["BTRReadEnd"]
#                 ),
#                 axis=1,
#             )
#             annotated_dfs.append(annotated_df)

#     concat_annotated_df = pd.concat(annotated_dfs)
#     return concat_annotated_df

# %%
# annotated_test_concat_alignments_df = find_btg_gene_coords(test_concat_alignments_df)
# annotated_test_concat_alignments_df

# %%

# %%
# test_aligned_pairs_dfs = []

# for bam_file in mapped_bam_files:

#     with pysam.AlignmentFile(
#         bam_file,
#         "rb",
#         threads=10,
#     ) as samfile:
#         reads = [
#             read
#             for read in samfile
#             if read.query_name in test_concat_alignments_df["Read"].values
#         ]
#         reads_names = [read.query_name for read in reads]
#         aligned_pairs = [
#             read.get_aligned_pairs(matches_only=False, with_cigar=True)
#             for read in reads
#         ]
#         df = pd.Series(aligned_pairs, index=reads_names)
#         test_aligned_pairs_dfs.append(df)

# concat_test_aligned_pairs_df = pd.concat(test_aligned_pairs_dfs)
# concat_test_aligned_pairs_df

# %%
# def get_genomic_coord_for_read_coord(aligned_pairs, required_read_pos):
#     for read_pos, ref_pos, cigar_op in aligned_pairs:
#         if read_pos == required_read_pos:
#             return ref_pos
#     return None

# %%
# read = "m64296e_241222_071206/1376592/ccs"

# test_concat_alignments_df.loc[test_concat_alignments_df["Read"].eq(read)]

# %%
# concat_test_aligned_pairs_df.loc[read]

# %%
# test_concat_alignments_df["BTGGeneStart2"] = test_concat_alignments_df.apply(
#     lambda x: get_genomic_coord_for_read_coord(
#         concat_test_aligned_pairs_df.loc[x["Read"]], x["BTRReadStart"]
#     ),
#     axis=1,
# )

# test_concat_alignments_df

# %%
# df = pd.DataFrame(
#     concat_test_aligned_pairs_df.loc["m64296e_241222_071206/32115244/ccs"],
#     columns=["ReadPos", "RefPos", "CigarOp"],
# ).set_index("ReadPos")
# df

# %%
# df.loc[571 - 2 : 571 + 2]

# %%
# 1238 + 571

# %%
# df.loc[df["CigarOp"] == 2]

# %%
# df.loc[:571, "CigarOp"].value_counts()

# %%
# df["CigarOp"].value_counts()

# %%

# %%

# %%

# %%
# single_barcode_min_umi_seq_length = 12
# single_barcode_max_umi_seq_length = 14

# %%

# %%
# concat_alignments_df.loc[
#     (concat_alignments_df["Barcode"] != "PCR")
#     & (concat_alignments_df["Gene"].eq(concat_alignments_df["MappedGene"]))
#     & (concat_alignments_df["Gene"].eq(concat_alignments_df["Barcode"]))
#     # & (concat_alignments_df["Read"].eq(read))
# ]

# %%
# concat_alignments_df.loc[
#     (concat_alignments_df["Barcode"] != "PCR")
#     & (concat_alignments_df["Gene"].eq(concat_alignments_df["MappedGene"]))
#     & (concat_alignments_df["Gene"].eq(concat_alignments_df["Barcode"]))
#     & (concat_alignments_df["BTRGaps"].gt(0))
# ]

# %%
# read = "m64296e_241222_071206/100007945/ccs"


# concat_alignments_df.loc[
#     (concat_alignments_df["Barcode"] != "PCR")
#     & (concat_alignments_df["Gene"].eq(concat_alignments_df["MappedGene"]))
#     & (concat_alignments_df["Gene"].eq(concat_alignments_df["Barcode"]))
#     & (concat_alignments_df["Read"].eq(read))
# ]

# %%
# read = "m64296e_241222_071206/100008153/ccs"

# concat_alignments_df.loc[
#     (concat_alignments_df["Barcode"] != "PCR")
#     & (concat_alignments_df["Gene"].eq(concat_alignments_df["MappedGene"]))
#     & (concat_alignments_df["Gene"].eq(concat_alignments_df["Barcode"]))
#     & (concat_alignments_df["Read"].eq(read))
# ]

# %%
# read = "m64296e_241222_071206/100073695/ccs"

# concat_alignments_df.loc[
#     (concat_alignments_df["Barcode"] != "PCR")
#     & (concat_alignments_df["Gene"].eq(concat_alignments_df["MappedGene"]))
#     & (concat_alignments_df["Gene"].eq(concat_alignments_df["Barcode"]))
#     & (concat_alignments_df["Read"].eq(read))
# ]

# %%

# %%

# %%

# %% [markdown]
# ## Construction

# %%
expected_barcode_locations_on_gene_dict


# %%
# def expected_barcode_to_gene_coverage(
#     observed_start: int,
#     observed_inclusive_end: int,
#     expected_start: int,
#     expected_exclusive_end: int,
# ):
#     """
#     Calculate the expected coverage of the barcode on the gene based on the expected start and end positions.
#     """
#     if pd.isna(
#         [observed_start, observed_inclusive_end, expected_start, expected_exclusive_end]
#     ).any():
#         return np.nan
#     observed_exclusive_end = observed_inclusive_end + 1  # Adjust for inclusive end
#     if (
#         expected_exclusive_end < observed_start
#         or observed_exclusive_end < expected_start
#     ):
#         return 0.0
#     coverage = min(observed_exclusive_end, expected_exclusive_end) - max(
#         observed_start, expected_start
#     )
#     return 100 * coverage / (expected_exclusive_end - expected_start)

# %%
def expected_barcode_to_gene_coverage(
    observed_start: int,
    observed_end: int,
    expected_start: int,
    expected_end: int,
):
    """
    Calculate the expected coverage of the barcode on the gene based on the expected start and end positions.
    """
    if pd.isna([observed_start, observed_end, expected_start, expected_end]).any():
        return np.nan
    if expected_end < observed_start or observed_end < expected_start:
        return 0.0
    coverage = min(observed_end, expected_end) - max(observed_start, expected_start)
    return 100 * coverage / (expected_end - expected_start)


# %%
gene_specific_concat_alignments_df = concat_alignments_df.loc[
    (concat_alignments_df["Barcode"] != "PCR")
    & (concat_alignments_df["Gene"].eq(concat_alignments_df["MappedGene"]))
    & (concat_alignments_df["Gene"].eq(concat_alignments_df["Barcode"]))
    # & (concat_alignments_df["BTRStrand"].eq("-"))
].reset_index(drop=True)

gene_specific_concat_alignments_df["BTGExpectedGeneStart"] = (
    gene_specific_concat_alignments_df.apply(
        lambda x: expected_barcode_locations_on_gene_dict[x["Gene"]][0], axis=1
    )
)
gene_specific_concat_alignments_df["BTGExpectedGeneEnd"] = (
    gene_specific_concat_alignments_df.apply(
        lambda x: expected_barcode_locations_on_gene_dict[x["Gene"]][1], axis=1
    )
)

gene_specific_concat_alignments_df["BTG%GeneCoverage"] = (
    gene_specific_concat_alignments_df.apply(
        lambda x: expected_barcode_to_gene_coverage(
            x["BTGGeneStart"],
            x["BTGGeneEnd"],
            x["BTGExpectedGeneStart"],
            x["BTGExpectedGeneEnd"],
        ),
        axis=1,
    )
)

gene_specific_concat_alignments_df

# %%
gene_specific_concat_alignments_df["BTG%GeneCoverage"].value_counts(dropna=False)

# %%
gene_specific_concat_alignments_df.loc[
    gene_specific_concat_alignments_df["BTGGeneStart"]
    .sub(gene_specific_concat_alignments_df["BTGExpectedGeneStart"])
    .abs()
    .le(2)
]

# %%
gene_specific_concat_alignments_df.loc[
    gene_specific_concat_alignments_df["BTGGeneEnd"]
    .sub(gene_specific_concat_alignments_df["BTGExpectedGeneEnd"])
    .abs()
    .le(1)
]

# %%
gene_specific_concat_alignments_df.loc[
    gene_specific_concat_alignments_df["BTG%GeneCoverage"].ge(30)
]

# %%
gene_specific_concat_alignments_df.groupby("Gene")["BTG%GeneCoverage"].describe()

# %%
gene_specific_concat_alignments_df.groupby("Gene")["BTGGeneStart"].describe()

# %%
# gene_specific_concat_alignments_df.insert(
#     gene_specific_concat_alignments_df.columns.get_loc("BTRReadCoords") + 3,
#     "%BTRRelLocOfReadStart",
#     gene_specific_concat_alignments_df["BTRReadStart"]
#     .mul(100)
#     .div(gene_specific_concat_alignments_df["ReadSeqLength"]),
# )

# gene_specific_concat_alignments_df.insert(
#     gene_specific_concat_alignments_df.columns.get_loc("MappedGene") + 1,
#     "ORFStart",
#     gene_specific_concat_alignments_df["Gene"].apply(
#         lambda x: genes_orfs_boundries_dict[x][0]
#     ),
# )
# gene_specific_concat_alignments_df.insert(
#     gene_specific_concat_alignments_df.columns.get_loc("MappedGene") + 2,
#     "ORFEnd",
#     gene_specific_concat_alignments_df["Gene"].apply(
#         lambda x: genes_orfs_boundries_dict[x][1]
#     ),
# )
# gene_specific_concat_alignments_df.insert(
#     gene_specific_concat_alignments_df.columns.get_loc("MappedGene") + 3,
#     "ORFLength",
#     gene_specific_concat_alignments_df.apply(
#         lambda x: x["ORFEnd"] - x["ORFStart"], axis=1
#     ),
# )


# gene_specific_concat_alignments_df.insert(
#     gene_specific_concat_alignments_df.columns.get_loc("RTGGeneEnd") + 1,
#     "%RTGORFCoverage",
#     gene_specific_concat_alignments_df.apply(
#         lambda x: 100
#         * (min(x["ORFEnd"], x["RTGGeneEnd"]) - max(x["ORFStart"], x["RTGGeneStart"]))
#         / x["ORFLength"],
#         axis=1,
#     ),
# )

# gene_specific_concat_alignments_df["UMISeq"] = gene_specific_concat_alignments_df.apply(
#     lambda x: x["ReadSeq"][
#         x["BTRReadEnd"] : min(
#             x["BTRReadEnd"] + single_barcode_max_umi_seq_length, x["ReadSeqLength"]
#         )
#     ],
#     axis=1,
# )
# gene_specific_concat_alignments_df["UMISeqLength"] = gene_specific_concat_alignments_df[
#     "UMISeq"
# ].apply(len)
# gene_specific_concat_alignments_df[
#     "UMIUniqueSubSeqs"
# ] = gene_specific_concat_alignments_df["UMISeq"].apply(
#     lambda x: split_umi_seq_to_unique_sub_seqs(x, single_barcode_min_umi_seq_length)
# )


# gene_specific_concat_alignments_df["%RTGORFCoverageUntilBarcode"] = (
#     gene_specific_concat_alignments_df.apply(
#         lambda x: 100
#         * (
#             min(x["ORFEnd"], x["RTGGeneEnd"], x["BTGGeneStart"])
#             - max(x["ORFStart"], x["RTGGeneStart"])
#         )
#         / x["ORFLength"],
#         axis=1,
#     )
# )
# gene_specific_concat_alignments_df["BTGGeneStart-ORFEnd"] = (
#     gene_specific_concat_alignments_df["BTGGeneStart"].sub(
#         gene_specific_concat_alignments_df["ORFEnd"]
#     )
# )

# gene_specific_concat_alignments_df["RTGReadTailLength"] = (
#     gene_specific_concat_alignments_df["ReadSeqLength"].sub(
#         gene_specific_concat_alignments_df["RTGReadEnd"]
#     )
# )
# gene_specific_concat_alignments_df["BTRBarcodeInTail"] = (
#     gene_specific_concat_alignments_df["BTRReadStart"].ge(
#         gene_specific_concat_alignments_df["RTGReadEnd"]
#     )
# )

# gene_specific_concat_alignments_df["BTGGeneStart-RTGGeneEnd"] = (
#     gene_specific_concat_alignments_df["BTGGeneStart"].sub(
#         gene_specific_concat_alignments_df["RTGGeneEnd"]
#     )
# )

# gene_specific_concat_alignments_df

# %% [markdown]
# ## Tests

# %%
# gene_specific_concat_alignments_df.groupby("Gene")["BTRBarcodeInTail"].value_counts(
#     normalize=True
# ).mul(100).round(2)

# %%
# gene_specific_concat_alignments_df.groupby("Gene")["BTGGeneStart-RTGGeneEnd"].describe()

# %%
# # marginal = "histogram"
# marginal = "box"

# fig = px.scatter(
#     gene_specific_concat_alignments_df.loc[
#         gene_specific_concat_alignments_df["BTGGeneStart"].ge(0)
#     ],
#     # df,
#     # x="UMILength",
#     x="BTGGeneStart",
#     y="RTGGeneEnd",
#     color="Gene",
#     marginal_x=marginal,
#     marginal_y=marginal,
#     # facet_col="Gene",
#     # facet_row="BTRBarcodeInTail",
#     # color="Repeat",
#     # facet_row="Sample",
#     # histnorm="percent",
#     # cumulative=True,
#     # labels={"UMISeqLength": "UMI length"},
#     # log_y=True,
# )
# # fig.update_xaxes(range=[0, None])
# fig.update_layout(
#     width=600,
#     #   height=350
#     height=600,
# )
# fig.show()

# %%
# fig = px.histogram(
#     gene_specific_concat_alignments_df.loc[
#         (gene_specific_concat_alignments_df["BTGGeneStart"].ge(0))
#         # & (~gene_specific_concat_alignments_df["BTRBarcodeInTail"])
#     ],
#     # df,
#     # x="UMILength",
#     x="RTGGeneEnd",
#     # color="Gene",
#     facet_col="Gene",
#     # facet_row="BTRBarcodeInTail",
#     # color="BTRBarcodeInTail",
#     # color="Repeat",
#     # facet_row="Sample",
#     # histnorm="percent",
#     # cumulative=True,
#     # labels={"UMISeqLength": "UMI length"},
#     # log_x=True,
#     log_y=True,
#     # opacity=0.7,
# )
# # fig.update_xaxes(range=[0, None], tick0=0, dtick=100)
# fig.update_layout(
#     width=800,
#     #   height=350
#     height=350,
#     barmode="overlay",
#     # barmode="stack",
# )
# fig.show()

# %%
# fig = px.histogram(
#     gene_specific_concat_alignments_df.loc[
#         (gene_specific_concat_alignments_df["BTGGeneStart"].ge(0))
#         # & (~gene_specific_concat_alignments_df["BTRBarcodeInTail"])
#     ],
#     # df,
#     # x="UMILength",
#     x="RTGGeneEnd",
#     # color="Gene",
#     facet_col="Gene",
#     # facet_row="BTRBarcodeInTail",
#     # color="BTRBarcodeInTail",
#     # color="Repeat",
#     # facet_row="Sample",
#     histnorm="percent",
#     # cumulative=True,
#     # labels={"UMISeqLength": "UMI length"},
#     # log_x=True,
#     log_y=True,
#     # opacity=0.7,
# )
# # fig.update_xaxes(range=[0, None], tick0=0, dtick=100)
# fig.update_layout(
#     width=800,
#     #   height=350
#     height=350,
#     barmode="overlay",
#     # barmode="stack",
# )
# fig.show()

# %%
# fig = px.histogram(
#     gene_specific_concat_alignments_df.loc[
#         (gene_specific_concat_alignments_df["BTGGeneStart"].ge(0))
#         & (~gene_specific_concat_alignments_df["BTRBarcodeInTail"])
#     ],
#     # df,
#     # x="UMILength",
#     x="BTGGeneStart-RTGGeneEnd",
#     # color="Gene",
#     facet_col="Gene",
#     # facet_row="BTRBarcodeInTail",
#     # color="BTRBarcodeInTail",
#     # color="Repeat",
#     # facet_row="Sample",
#     # histnorm="percent",
#     # cumulative=True,
#     # labels={"UMISeqLength": "UMI length"},
#     # log_x=True,
#     log_y=True,
#     # opacity=0.7,
# )
# fig.update_xaxes(range=[0, None], tick0=0, dtick=100)
# fig.update_layout(
#     width=800,
#     #   height=350
#     height=350,
#     barmode="overlay",
#     # barmode="stack",
# )
# fig.show()

# %%
# fig = px.histogram(
#     gene_specific_concat_alignments_df.loc[
#         (gene_specific_concat_alignments_df["BTGGeneStart"].ge(0))
#         # & (~gene_specific_concat_alignments_df["BTRBarcodeInTail"])
#     ],
#     # df,
#     # x="UMILength",
#     x="BTGGeneStart-RTGGeneEnd",
#     # color="Gene",
#     # facet_row="Gene",
#     # # facet_row="BTRBarcodeInTail",
#     # # color="BTRBarcodeInTail",
#     # facet_col="BTRBarcodeInTail",
#     facet_col="Gene",
#     # facet_row="BTRBarcodeInTail",
#     # color="BTRBarcodeInTail",
#     facet_row="BTRBarcodeInTail",
#     # color="Repeat",
#     # facet_row="Sample",
#     # histnorm="percent",
#     # cumulative=True,
#     # labels={"UMISeqLength": "UMI length"},
#     # log_x=True,
#     log_y=True,
#     # opacity=0.7,
# )
# # fig.update_xaxes(
# #     # range=[0, None], tick0=0,
# #     dtick=20
# # )
# fig.update_layout(
#     width=800,
#     #   height=350
#     height=550,
#     # barmode="overlay",
#     # barmode="stack",
# )
# fig.show()

# %%
# fig = px.histogram(
#     gene_specific_concat_alignments_df.loc[
#         (gene_specific_concat_alignments_df["BTGGeneStart"].ge(0))
#         & (~gene_specific_concat_alignments_df["BTRBarcodeInTail"])
#     ],
#     # df,
#     # x="UMILength",
#     x="BTGGeneStart-RTGGeneEnd",
#     # color="Gene",
#     facet_row="Gene",
#     # facet_row="BTRBarcodeInTail",
#     # color="BTRBarcodeInTail",
#     # color="Repeat",
#     # facet_row="Sample",
#     histnorm="percent",
#     # cumulative=True,
#     # labels={"UMISeqLength": "UMI length"},
#     # log_x=True,
#     # log_y=True,
#     # opacity=0.7,
# )
# # fig.update_xaxes(range=[0, None], tick0=0, dtick=20)
# fig.update_layout(
#     width=1400,
#     #   height=350
#     height=550,
#     # barmode="overlay",
#     # barmode="stack",
# )
# fig.show()

# %%
# fig = px.histogram(
#     gene_specific_concat_alignments_df.loc[
#         gene_specific_concat_alignments_df["BTGGeneStart"].ge(0)
#     ],
#     # df,
#     # x="UMILength",
#     x="BTGGeneStart",
#     # color="Gene",
#     facet_col="Gene",
#     # facet_row="BTRBarcodeInTail",
#     color="BTRBarcodeInTail",
#     # color="Repeat",
#     # facet_row="Sample",
#     # histnorm="percent",
#     # cumulative=True,
#     # labels={"UMISeqLength": "UMI length"},
#     log_y=True,
#     opacity=0.7,
# )
# fig.update_xaxes(range=[0, None], tick0=0, dtick=1000)
# fig.update_layout(
#     width=800,
#     #   height=350
#     height=450,
#     barmode="overlay",
#     # barmode="stack",
# )
# fig.show()

# %%
# fig = px.histogram(
#     gene_specific_concat_alignments_df.loc[
#         gene_specific_concat_alignments_df["BTGGeneStart"].ge(0)
#     ],
#     # df,
#     # x="UMILength",
#     x="BTGGeneStart",
#     # color="Gene",
#     facet_col="Gene",
#     # color="BTRBarcodeInTail",
#     # facet_row="BTRBarcodeInTail",
#     # color="Repeat",
#     # facet_row="Sample",
#     histnorm="percent",
#     # cumulative=True,
#     # labels={"UMISeqLength": "UMI length"},
#     log_y=True,
#     opacity=0.7,
# )
# fig.update_xaxes(range=[0, None], tick0=0, dtick=1000)
# fig.update_layout(
#     width=800,
#     height=350,
#     # height=600,
#     # height=450,
#     # barmode="overlay",
#     # barmode="stack",
# )
# fig.show()

# %%
# gene_specific_concat_alignments_df.groupby("Gene")["BTGGeneStart"].describe()

# %%
# gene_specific_concat_alignments_df.loc[
#     gene_specific_concat_alignments_df["BTGGeneStart"].lt(0),
#     # [
#     #     "Sample",
#     #     "Gene",
#     #     "Repeat",
#     #     "Read",
#     #     "ReadSeq",
#     #     "RTGStrand",
#     #     "RTGGeneStart",
#     #     "BTRReadStart",
#     #     "RTGReadStart",
#     #     "BTGGeneStart",
#     # ],
# ]

# %%
# Out[279]["RTGGeneStart"].describe()

# %%
# Out[279]["BTRReadStart"].describe()

# %%
# Out[279]["RTGReadStart"].describe()

# %%
# Out[279]["BTRReadStart"] - Out[279]["RTGReadStart"]

# %%
# fig = px.histogram(
#     gene_specific_concat_alignments_df,
#     # df,
#     # x="UMILength",
#     x="BTGGeneStart",
#     # color="Gene",
#     facet_col="Gene",
#     # facet_row="UMISeqLength>20",
#     # color="Repeat",
#     # facet_row="Sample",
#     # histnorm="percent",
#     # cumulative=True,
#     # labels={"UMISeqLength": "UMI length"},
# )
# # fig.update_xaxes(dtick=1)
# fig.update_layout(width=800, height=350)
# fig.show()

# %%
# gene_specific_concat_alignments_df["BTRReadStart"]

# %%
# (
#     gene_specific_concat_alignments_df["ReadSeqLength"].sub(
#         gene_specific_concat_alignments_df["BTRReadStart"]
#     )
# ).describe()

# %%
# (
#     (
#         gene_specific_concat_alignments_df["ReadSeqLength"].sub(
#             gene_specific_concat_alignments_df["BTRReadStart"]
#         )
#     )
#     .mul(100)
#     .div(gene_specific_concat_alignments_df["ReadSeqLength"])
# ).describe()

# %%
# (
#     (
#         gene_specific_concat_alignments_df["ReadSeqLength"].sub(
#             gene_specific_concat_alignments_df["BTRReadStart"]
#         )
#     )
#     .mul(100)
#     .div(gene_specific_concat_alignments_df["ReadSeqLength"])
# ).add(gene_specific_concat_alignments_df["%BTRRelLocOfReadStart"]).describe()

# %%
# gene_specific_concat_alignments_df.loc[
#     (
#         gene_specific_concat_alignments_df["Gene"].eq(
#             gene_specific_concat_alignments_df["MappedGene"]
#         )
#     )
#     & (
#         gene_specific_concat_alignments_df["Gene"].eq(
#             gene_specific_concat_alignments_df["Barcode"]
#         )
#     ),
#     # & (gene_specific_concat_alignments_df["BTRStrand"].eq("-"))
#     "RTGReadTailLength",
# ].describe().round(2)

# %%
# gene_specific_concat_alignments_df.loc[
#     (
#         gene_specific_concat_alignments_df["Gene"].eq(
#             gene_specific_concat_alignments_df["MappedGene"]
#         )
#     )
#     & (
#         gene_specific_concat_alignments_df["Gene"].eq(
#             gene_specific_concat_alignments_df["Barcode"]
#         )
#     ),
#     # & (gene_specific_concat_alignments_df["BTRStrand"].eq("-"))
#     "RTGReadTailLength",
# ].describe().round(2)

# %%
# gene_specific_concat_alignments_df.loc[
#     (
#         gene_specific_concat_alignments_df["Gene"].eq(
#             gene_specific_concat_alignments_df["MappedGene"]
#         )
#     )
#     & (
#         gene_specific_concat_alignments_df["Gene"].eq(
#             gene_specific_concat_alignments_df["Barcode"]
#         )
#     )
#     & (gene_specific_concat_alignments_df["BTRStrand"].eq("-")),
#     # "RTGReadTail",
# ].groupby("Gene")["RTGReadTailLength"].describe().round(2)

# %%
# gene_specific_concat_alignments_df.loc[
#     (
#         gene_specific_concat_alignments_df["Gene"].eq(
#             gene_specific_concat_alignments_df["MappedGene"]
#         )
#     )
#     & (
#         gene_specific_concat_alignments_df["Gene"].eq(
#             gene_specific_concat_alignments_df["Barcode"]
#         )
#     )
#     & (gene_specific_concat_alignments_df["BTRStrand"].eq("-"))
#     & (gene_specific_concat_alignments_df["RTGReadTailLength"].ge(100)),
#     # "RTGReadTail",
# ].groupby("Gene")["RTGReadTailLength"].describe().round(2)

# %%
# gene_specific_concat_alignments_df.loc[
#     (
#         gene_specific_concat_alignments_df["Gene"].eq(
#             gene_specific_concat_alignments_df["MappedGene"]
#         )
#     )
#     & (
#         gene_specific_concat_alignments_df["Gene"].eq(
#             gene_specific_concat_alignments_df["Barcode"]
#         )
#     )
#     & (gene_specific_concat_alignments_df["BTRStrand"].eq("-"))
#     & (gene_specific_concat_alignments_df["RTGReadTailLength"].ge(33)),
#     # "RTGReadTail",
# ].groupby("Gene")["RTGReadTailLength"].describe().round(2)["count"].sum()

# %%
# 142493.0 + 169132.0

# %%
# pcr_primer_len = len(primers_dict["PCR"])
# for gene in genes:
#     primers_len_without_pcr = len(primers_dict[gene]) + 10  # minimal UMI length
#     primers_len_with_pcr = primers_len_without_pcr + pcr_primer_len
#     ic(gene, primers_len_without_pcr, primers_len_with_pcr)

# %%
# gene_specific_concat_alignments_df.loc[
#     (
#         (gene_specific_concat_alignments_df["Gene"].eq("ADAR1"))
#         & (gene_specific_concat_alignments_df["RTGReadTailLength"].ge(56))
#     )
#     | (
#         (gene_specific_concat_alignments_df["Gene"].eq("IQEC"))
#         & (gene_specific_concat_alignments_df["RTGReadTailLength"].ge(54))
#     ),
#     # "RTGReadTail",
# ]

# %%
# gene_specific_concat_alignments_df.loc[
#     (
#         (gene_specific_concat_alignments_df["Gene"].eq("ADAR1"))
#         & (gene_specific_concat_alignments_df["RTGReadTailLength"].ge(56))
#     )
#     | (
#         (gene_specific_concat_alignments_df["Gene"].eq("IQEC"))
#         & (gene_specific_concat_alignments_df["RTGReadTailLength"].ge(54))
#     ),
#     # "RTGReadTail",
# ].groupby("Gene")["RTGReadTailLength"].describe().round(2)

# %%
# fig = px.histogram(
#     gene_specific_concat_alignments_df.loc[
#         (
#             gene_specific_concat_alignments_df["Gene"].eq(
#                 gene_specific_concat_alignments_df["MappedGene"]
#             )
#         )
#         & (
#             gene_specific_concat_alignments_df["Gene"].eq(
#                 gene_specific_concat_alignments_df["Barcode"]
#             )
#         )
#         & (gene_specific_concat_alignments_df["BTRStrand"].eq("-"))
#         & (gene_specific_concat_alignments_df["RTGReadTailLength"].le(100))
#     ],
#     x="RTGReadTailLength",
#     facet_col="Gene",
#     # facet_row="BTRStrand",
#     # log_x=True,
#     log_y=True,
# )
# fig.update_xaxes(dtick=5)
# fig.update_layout(width=800, height=350, template="plotly_white")
# fig.show()

# %%
# fig = px.histogram(
#     gene_specific_concat_alignments_df.loc[
#         (
#             gene_specific_concat_alignments_df["Gene"].eq(
#                 gene_specific_concat_alignments_df["MappedGene"]
#             )
#         )
#         & (
#             gene_specific_concat_alignments_df["Gene"].eq(
#                 gene_specific_concat_alignments_df["Barcode"]
#             )
#         )
#         # & (gene_specific_concat_alignments_df["BTRStrand"].eq("-"))
#     ],
#     x="BTGGeneStart-ORFEnd",
#     facet_col="Gene",
#     facet_row="BTRStrand",
#     log_y=True,
# )
# fig.update_layout(width=800, height=500, template="plotly_white")
# fig.show()

# %%
# (
#     gene_specific_concat_alignments_df.loc[
#         (
#             gene_specific_concat_alignments_df["Gene"].eq(
#                 gene_specific_concat_alignments_df["MappedGene"]
#             )
#         )
#         & (
#             gene_specific_concat_alignments_df["Gene"].eq(
#                 gene_specific_concat_alignments_df["Barcode"]
#             )
#         )
#         # & (gene_specific_concat_alignments_df["BTRStrand"].eq("-"))
#     ]
#     .groupby(
#         [
#             "Gene",
#             "BTRStrand",
#         ]
#     )["BTGGeneStart-ORFEnd"]
#     .describe()
#     .round(2)
# )

# %%
# gene_specific_concat_alignments_df.loc[
#     (
#         gene_specific_concat_alignments_df["Gene"].eq(
#             gene_specific_concat_alignments_df["MappedGene"]
#         )
#     )
#     & (
#         gene_specific_concat_alignments_df["Gene"].eq(
#             gene_specific_concat_alignments_df["Barcode"]
#         )
#     )
#     & (gene_specific_concat_alignments_df["%RTGORFCoverageUntilBarcode"] < 0),
#     "%RTGORFCoverageUntilBarcode",
# ].describe()

# %%
# genes_orfs_boundries_dict

# %%
# fig = px.histogram(
#     gene_specific_concat_alignments_df.loc[
#         (
#             gene_specific_concat_alignments_df["Gene"].eq(
#                 gene_specific_concat_alignments_df["MappedGene"]
#             )
#         )
#         & (
#             gene_specific_concat_alignments_df["Gene"].eq(
#                 gene_specific_concat_alignments_df["Barcode"]
#             )
#         )
#     ],
#     x="%RTGORFCoverageUntilBarcode",
#     facet_col="Gene",
#     histnorm="percent",
#     log_y=True,
#     cumulative=True,
# )
# fig.update_xaxes(dtick=10)
# fig.update_layout(width=800, height=400, template="plotly_white")
# fig.show()

# %%
# gene_specific_concat_alignments_df.loc[
#     (
#         gene_specific_concat_alignments_df["BTGGeneStart"].ge(
#             gene_specific_concat_alignments_df["ORFEnd"]
#         )
#     )
# ]

# %%
# gene_specific_concat_alignments_df.loc[
#     (
#         gene_specific_concat_alignments_df["BTGGeneStart"].lt(
#             gene_specific_concat_alignments_df["ORFEnd"]
#         )
#     )
# ]

# %%

# %%
# gene_specific_concat_alignments_df.loc[
#     (
#         gene_specific_concat_alignments_df["BTGGeneStart"].ge(
#             gene_specific_concat_alignments_df["ORFEnd"]
#         )
#     )
#     & (
#         gene_specific_concat_alignments_df["BTRReadStart"].ge(
#             gene_specific_concat_alignments_df["RTGReadEnd"]
#         )
#     )
#     & (
#         gene_specific_concat_alignments_df["RTGGeneEnd"].ge(
#             gene_specific_concat_alignments_df["ORFEnd"]
#         )
#     )
# ]

# %%
# gene_specific_concat_alignments_df.loc[
#     (
#         gene_specific_concat_alignments_df["BTRReadStart"].ge(
#             gene_specific_concat_alignments_df["RTGReadEnd"]
#         )
#     )
#     & (
#         gene_specific_concat_alignments_df["RTGGeneEnd"].ge(
#             gene_specific_concat_alignments_df["ORFEnd"]
#         )
#     )
# ]

# %%
# gene_specific_concat_alignments_df.apply(
#     lambda x: 
# )

# %%
# genes_orfs_boundries_dict

# %%
# gene_specific_concat_alignments_df.drop_duplicates(["Gene", "Repeat", "Read"]).shape[0]

# %%
# {
#     "TargetName": "Read",
#     "QueryName": "Barcode",
#     "TargetSeq": "ReadSeq",
#     "GenomicTargetStrand": "RTGStrand",
#     "TargetStrand": "BTRStrand",
#     "%QueryIdentity": "%BTRBarcodeIdentity",
#     "%QueryCoverage": "%BTRBarcodeCoverage",
#     "TargetSeqLength": "ReadSeqLength",
#     "QuerySeqLength": "BarcodeSeqLength",
#     "AlignmentLength": "BTRAlignmentLength",
#     "AlignedTargetCoords": "BTRReadCoords",
#     "AlignedQueryCoords": "BTRBarcodeCoords",
#     "AlignmentObject": "BTRAlignmentObject",
#     "NumOfAlignedTargetCoords": "NumOfBTRReadCoords",
#     "NumOfAlignedQueryCoords": "NumOfBTRBarcodeCoords",
#     "QueryGapOpenings": "BTRBarcodeGapOpenings",
#     "Score": "BTRScore",
#     "Gaps": "BTRGaps",
#     "Identitites": "BTRIdentitites",
#     "Mismatches": "BTRMismatches",
#     "TargetStart": "BTRReadStart",
#     "TargetEnd": "BTRReadEnd",
#     "%RelLocOfTargetStart": "%BTRRelLocOfReadStart",
# }

# %% [markdown]
# ## Actually choosing

# %%
gene_specific_concat_alignments_df


# %%
def choose_best_btr_barcode_alignment_for_read(
    read_df: pd.DataFrame,
) -> pd.Series:
    """
    Choose the best barcode-to-read alignment for this read based on the criteria defined.
    """
    # if read_df.empty:
    #     return pd.Series()

    criteria_cols_and_optimization_funcs = [
        ("%BTRBarcodeCoverage", max),
        ("%BTRBarcodeIdentity", max),
        ("NumOfBTRBarcodeGapOpenings", min),
        ("BTRGaps", min),
    ]
    optimized_cols = []
    is_optimized_cols = []
    for col, optimization_func in criteria_cols_and_optimization_funcs:
        optimized_col = f"{col}_best"
        is_optimized_col = f"{col}_best_is_optimized"
        optimized_cols.append(optimized_col)
        is_optimized_cols.append(is_optimized_col)
        read_df[optimized_col] = optimization_func(read_df[col])
        read_df[is_optimized_col] = read_df[col].eq(read_df[optimized_col])

    read_df["NumOfOptimizedCols"] = read_df[is_optimized_cols].sum(axis=1)

    criteria_cols = []
    ascending_sortings = []
    for col, optimization_func in criteria_cols_and_optimization_funcs:

        criteria_cols.append(optimized_col)
        ascending_sortings.append(True if optimization_func == max else False)

    read_df = (
        read_df.loc[
            read_df["NumOfOptimizedCols"].eq(read_df["NumOfOptimizedCols"].max())
        ]
        .sort_values(by=criteria_cols, ascending=ascending_sortings)
        .drop(
            columns=criteria_cols
            + optimized_cols
            + is_optimized_cols
            + ["NumOfOptimizedCols"]
        )
        .iloc[0]
    )

    return read_df


# %%
def split_umi_seq_to_unique_sub_seqs(umi_seq, min_umi_seq_len):
    return {
        umi_seq[x : x + min_umi_seq_len]
        for x in range(0, len(umi_seq) - min_umi_seq_len + 1)
    }


# %%
# gene-repeat-read combinations with (one or two) best combined tags' alignment(s) w.r.t strand

max_abs_distance_from_expected_gene_barcode_start_or_end = 2

# min_prct_btg_gene_cov = 85

# min_prct_btr_barcode_cov = 85
# min_prct_btr_barcode_identity = 85
# min_prct_btr_barcode_cov = 95
min_prct_btr_barcode_cov = 100
# min_prct_btr_barcode_identity = 95
min_prct_btr_barcode_identity = 100

# min_rel_loc_of_target_start = 90

# max_gap_openings_per_btr_barcode = 1

# # max_total_gaps_per_btr_barcode = 4
# max_total_gaps_per_btr_barcode = 1


# possible_gap_openings_per_btr_barcode = list(
#     range(0, max_gap_openings_per_btr_barcode + 1)
# )
# possible_total_gaps_per_btr_barcode = list(range(0, max_total_gaps_per_btr_barcode + 1))


# max_umi_seq_length = 12


best_gene_specific_concat_alignments_df = (
    gene_specific_concat_alignments_df.loc[
        (
            gene_specific_concat_alignments_df["BTGGeneStart"]
            .sub(gene_specific_concat_alignments_df["BTGExpectedGeneStart"])
            .abs()
            .le(max_abs_distance_from_expected_gene_barcode_start_or_end)
        )
        & (
            gene_specific_concat_alignments_df["BTGGeneEnd"]
            .sub(gene_specific_concat_alignments_df["BTGExpectedGeneEnd"])
            .abs()
            .le(max_abs_distance_from_expected_gene_barcode_start_or_end)
        )
        & (
            gene_specific_concat_alignments_df["%BTRBarcodeCoverage"].ge(
                min_prct_btr_barcode_cov
            )
        )
        & (
            gene_specific_concat_alignments_df["%BTRBarcodeIdentity"].ge(
                min_prct_btr_barcode_identity
            )
        )
        # & (
        #     gene_specific_concat_alignments_df["BTG%GeneCoverage"].ge(
        #         min_prct_btg_gene_cov
        #     )
        # )
        # & (
        #     gene_specific_concat_alignments_df["NumOfBTRBarcodeGapOpenings"].isin(
        #         possible_gap_openings_per_btr_barcode
        #     )
        # )
        # & (
        #     gene_specific_concat_alignments_df["BTRGaps"].isin(
        #         possible_total_gaps_per_btr_barcode
        #     )
        # )
    ]
    .sort_values(
        [
            "Gene",
            "Repeat",
            "Read",
            "MappedGene",
            "BTRStrand",
            # "StrQueries",
        ]
    )
    .reset_index(drop=True)
)

# make sure there are no duplicates - only one alignment per gene-repeat-read
assert (
        best_gene_specific_concat_alignments_df.drop_duplicates(
            ["Gene", "Repeat", "Read"]
        ).shape[0]
        == best_gene_specific_concat_alignments_df.shape[0]
    )

# # drop duplicates if needed - keep only best alignment per gene-repeat-read
# if (
#     best_gene_specific_concat_alignments_df.drop_duplicates(
#         ["Gene", "Repeat", "Read"]
#     ).shape[0]
#     != best_gene_specific_concat_alignments_df.shape[0]
# ):
#     best_gene_specific_concat_alignments_df["Read2"] = (
#         best_gene_specific_concat_alignments_df["Read"]
#     )
#     best_gene_specific_concat_alignments_df = (
#         best_gene_specific_concat_alignments_df.groupby("Read2")
#         .apply(choose_best_btr_barcode_alignment_for_read, include_groups=False)
#         .reset_index(drop=True)
#     )
#     assert (
#         best_gene_specific_concat_alignments_df.drop_duplicates(
#             ["Gene", "Repeat", "Read"]
#         ).shape[0]
#         == best_gene_specific_concat_alignments_df.shape[0]
#     )

assert (
    best_gene_specific_concat_alignments_df["BTRStrand"].value_counts()["-"]
    == best_gene_specific_concat_alignments_df.shape[0]
)

assert best_gene_specific_concat_alignments_df.loc[
    :,
    [
        "%BTRBarcodeCoverage",
        "%BTRBarcodeIdentity",
        "NumOfBTRBarcodeGapOpenings",
        "BTRGaps",
    ],
].min().eq([100, 100, 0, 0]).all(), "No gaps are expected since we required 100% coverage and identity."



best_gene_specific_concat_alignments_df

# %%
best_gene_specific_concat_alignments_df.groupby("Gene").size()

# %%
principle_exact_umi_seq_length = 12

best_gene_specific_concat_alignments_df["ExactUMISeq"] = (
    best_gene_specific_concat_alignments_df.apply(
        lambda x: x["ReadSeq"][
            x["BTRReadEnd"] : min(
                x["BTRReadEnd"] + principle_exact_umi_seq_length, x["ReadSeqLength"]
            )
        ],
        axis=1,
    )
)
best_gene_specific_concat_alignments_df["ExactUMISeqLength"] = (
    best_gene_specific_concat_alignments_df["ExactUMISeq"].apply(len)
)

# %%
best_gene_specific_concat_alignments_df.loc[
    best_gene_specific_concat_alignments_df["ExactUMISeqLength"].lt(
        principle_exact_umi_seq_length
    ),
].groupby("Gene").size()

# %%
best_gene_specific_concat_alignments_df.groupby("Gene")["ExactUMISeqLength"].describe()

# %%

# %%
min_spanning_umi_seq_len = 10
max_spanning_umi_seq_len = 14

# umi_sub_seq_len = 10
umi_sub_seq_len = 12


best_gene_specific_concat_alignments_df["SpanningUMISeq"] = (
    best_gene_specific_concat_alignments_df.apply(
        lambda x: x["ReadSeq"][
            x["BTRReadEnd"] : min(
                x["BTRReadEnd"] + max_spanning_umi_seq_len, x["ReadSeqLength"]
            )
        ],
        axis=1,
    )
)
best_gene_specific_concat_alignments_df["SpanningUMISeqLength"] = (
    best_gene_specific_concat_alignments_df["SpanningUMISeq"].apply(len)
)

best_gene_specific_concat_alignments_df["UMIUniqueSubSeqs"] = (
    best_gene_specific_concat_alignments_df["SpanningUMISeq"].apply(
        lambda x: split_umi_seq_to_unique_sub_seqs(x, umi_sub_seq_len)
    )
)

# %%
best_gene_specific_concat_alignments_df.loc[
    best_gene_specific_concat_alignments_df["SpanningUMISeqLength"].lt(
        min_spanning_umi_seq_len
    ),
].groupby("Gene").size()

# %%
best_gene_specific_concat_alignments_df.groupby("Gene")[
    "SpanningUMISeqLength"
].describe()

# %%
best_gene_specific_concat_alignments_df.groupby("Gene").apply(
    lambda x: x["UMIUniqueSubSeqs"].apply(len).describe(),
    include_groups=False,
)

# %%

# %%

# %%
best_gene_specific_concat_alignments_df.loc[
    best_gene_specific_concat_alignments_df["BTGGeneEnd"].isna()
]


# %% [markdown]
# ## Further tests

# %%
# i = 0

# row = best_gene_specific_concat_alignments_df.iloc[i]
# print(row)

# %%
# sample = row["Sample"]
# mapped_bam_file = [file for file in mapped_bam_files if sample in file.name][0]

# read_name = row["Read"]

# with pysam.AlignmentFile(
#     mapped_bam_file,
#     "rb",
#     threads=10,
# ) as samfile:
#     read = [read for read in samfile if read.query_name == read_name][0]
# # print()
# print(read)

# %%
# btr_read_start = row["BTRReadStart"]
# btr_read_end = row["BTRReadEnd"]
# btg_gene_start = row["BTGGeneStart"]
# btg_gene_end = row["BTGGeneEnd"]
# btg_expected_gene_start = row["BTGExpectedGeneStart"]
# btg_expected_gene_end = row["BTGExpectedGeneEnd"]

# pd.DataFrame(
#     [
#         [btr_read_start, btr_read_end],
#         [btg_gene_start, btg_gene_end],
#         [btg_expected_gene_start, btg_expected_gene_end],
#     ],
#     index=["BTRRead", "BTGGene", "BTGExpectedGene"],
#     columns=["Start", "End"],
# )

# %%
# 3389 - 3365

# %%
# aligned_pairs = read.get_aligned_pairs(matches_only=False, with_cigar=True)
# aligned_pairs

# %%
# aligned_pairs_df.loc[
#     aligned_pairs_df["ReadPos"].ge(3365) & aligned_pairs_df["ReadPos"].le(3389)
# ]

# %%
# aligned_pairs_df = pd.DataFrame(
#     read.get_aligned_pairs(matches_only=False, with_cigar=True),
#     columns=["ReadPos", "RefPos", "Cigar"],
# )
# aligned_pairs_df

# %%
# aligned_pairs

# %%
# def get_genomic_coord_for_read_coord(aligned_pairs, required_read_pos):
#     for read_pos, ref_pos, cigar_op in aligned_pairs:
#         if read_pos == required_read_pos:
#             return ref_pos
#     return None

# %%
# best_gene_specific_concat_alignments_df.groupby("Gene")["BTG%GeneCoverage"].describe()

# %%
# best_gene_specific_concat_alignments_df.loc[
#     best_gene_specific_concat_alignments_df["BTG%GeneCoverage"].ge(80)
# ]

# %%
# fig = px.histogram(
#     best_gene_specific_concat_alignments_df,
#     x="%BTRBarcodeCoverage",
#     facet_col="Gene",
#     # facet_row="BTRStrand",
#     # log_x=True,
#     # log_y=True,
#     histnorm="percent",
#     cumulative=True,
#     color="Gene",
# )
# fig.update_layout(showlegend=False)
# fig.show()

# %% [markdown]
# # Find unique reads by UMI sub-seqs in gene-specific PCR-amplified selected barcodes
#

# %% [markdown]
# ## Find unique acceptable PCR barcodes

# %%
def is_pcr_barcode_found(
    btr_barcode_coords, # "BTRBarcodeCoords_PCR"
    btr_barcode_identity, # "%BTRBarcodeIdentity_PCR"
    btr_barcode_coverage, # "%BTRBarcodeCoverage_PCR"
    btr_alignment_length, # "BTRAlignmentLength_PCR"
    pcr_barcode_seq_length, # "BarcodeSeqLength_PCR"
    gene_btr_read_end, # "BTRReadEnd"
    pcr_btr_read_start, # "BTRReadStart_PCR"
    acceptable_range_between_gene_and_pcr=range(10, 15), # range(expected_spanning_umi_seq_len - 2, expected_spanning_umi_seq_len + 3)
    min_regular_btr_barcode_coverage=90,
    min_regular_btr_barcode_identity=90,
    min_special_btr_alignment_length=10,
    debug=False,
):

    if (
        pcr_btr_read_start - gene_btr_read_end
    ) not in acceptable_range_between_gene_and_pcr:
        return False

    if (
        btr_barcode_coverage >= min_regular_btr_barcode_coverage
        and btr_barcode_identity >= min_regular_btr_barcode_identity
    ):
        return True

    btr_barcode_coords_end, btr_barcode_coords_start = btr_barcode_coords[0]

    if debug:
        ic(
            btr_barcode_coords,
            btr_barcode_identity,
            btr_barcode_coverage,
            btr_alignment_length,
            pcr_barcode_seq_length,
            gene_btr_read_end,
            pcr_btr_read_start,
            acceptable_range_between_gene_and_pcr,
            min_regular_btr_barcode_coverage,
            min_regular_btr_barcode_identity,
            min_special_btr_alignment_length,
        )

    if (
        # (
        #     (pcr_btr_read_start - gene_btr_read_end)
        #     in acceptable_range_between_gene_and_pcr
        # )
        # and
        (len(btr_barcode_coords) == 1)
        and (btr_barcode_coords_end == pcr_barcode_seq_length)
        and (btr_barcode_coords_end - btr_barcode_coords_start == btr_alignment_length)
        and (btr_alignment_length >= min_special_btr_alignment_length)
        and np.isclose(
            100 * btr_alignment_length / pcr_barcode_seq_length, btr_barcode_coverage
        )
    ):
        return True

    return False


# %%
best_gene_specific_pcr_amplified_concat_alignments_df = concat_alignments_df.loc[
    (concat_alignments_df["Barcode"].eq("PCR"))
    & (concat_alignments_df["BTRStrand"].eq("-"))
].merge(
    best_gene_specific_concat_alignments_df,
    on=["Sample", "Gene", "Repeat", "Read", "OldRead", "MappedGene"],
    how="inner",
    suffixes=["_PCR", None],
)

# min_spanning_umi_seq_len = 10
# max_spanning_umi_seq_len = 14
# umi_sub_seq_len = 10
umi_sub_seq_len = 10
expected_spanning_umi_seq_len = 12

pcr_barcode_found = best_gene_specific_pcr_amplified_concat_alignments_df.apply(
    lambda x: is_pcr_barcode_found(
        x["BTRBarcodeCoords_PCR"],
        x["%BTRBarcodeIdentity_PCR"],
        x["%BTRBarcodeCoverage_PCR"],
        x["BTRAlignmentLength_PCR"],
        x["BarcodeSeqLength_PCR"],
        x["BTRReadEnd"],
        x["BTRReadStart_PCR"],
        # min_special_btr_alignment_length=12,
        acceptable_range_between_gene_and_pcr=range(
            expected_spanning_umi_seq_len - 2, expected_spanning_umi_seq_len + 3
        ),  # range(10, 15)
    ),
    axis=1,
)
best_gene_specific_pcr_amplified_concat_alignments_df = (
    best_gene_specific_pcr_amplified_concat_alignments_df.loc[pcr_barcode_found]
)


best_gene_specific_pcr_amplified_concat_alignments_df["SpanningUMISeq"] = (
    best_gene_specific_pcr_amplified_concat_alignments_df.apply(
        lambda x: x["ReadSeq"][x["BTRReadEnd"] : x["BTRReadStart_PCR"]],
        axis=1,
    )
)
best_gene_specific_pcr_amplified_concat_alignments_df["SpanningUMISeqLength"] = (
    best_gene_specific_pcr_amplified_concat_alignments_df["SpanningUMISeq"].apply(len)
)
best_gene_specific_pcr_amplified_concat_alignments_df["UMIUniqueSubSeqs"] = (
    best_gene_specific_pcr_amplified_concat_alignments_df["SpanningUMISeq"].apply(
        lambda x: split_umi_seq_to_unique_sub_seqs(x, umi_sub_seq_len)
    )
)

best_gene_specific_pcr_amplified_concat_alignments_df

# %%
# best_gene_specific_pcr_amplified_concat_alignments_df["Read"].nunique()

best_gene_specific_pcr_amplified_concat_alignments_df.drop_duplicates(["Gene", "Read"]).shape[0]

# %%
best_gene_specific_pcr_amplified_concat_alignments_df["SpanningUMISeqLength"].describe()

# %%
best_gene_specific_pcr_amplified_concat_alignments_df[
    "SpanningUMISeqLength"
].value_counts().sort_index()

# %%
best_gene_specific_pcr_amplified_concat_alignments_df.loc[
    :, ["%BTRBarcodeIdentity_PCR", "%BTRBarcodeCoverage_PCR"]
].describe().round(2)

# %%
best_gene_specific_pcr_amplified_concat_alignments_df.loc[
    best_gene_specific_pcr_amplified_concat_alignments_df.duplicated(
        subset=["Gene", "Read"], keep=False
    )
]

# %%
best_gene_specific_pcr_amplified_concat_alignments_df.loc[
    best_gene_specific_pcr_amplified_concat_alignments_df.duplicated(
        subset=["Gene", "Read"], keep=False
    )
].loc[:, ["Gene", "Read"]].value_counts().describe().round(2)


# %%
def choose_best_pcr_btr_barcode_alignment_for_read(
    read_df: pd.DataFrame, seed: int, ideal_spanning_umi_seq_len: int = 12
) -> pd.Series:
    """
    Choose the best barcode-to-read alignment for this read based on the criteria defined.
    """
    if read_df.shape[0] == 1:
        # return read_df.iloc[0]
        return read_df
    
    # ic(seed, ideal_spanning_umi_seq_len)
    
    read_df["DistanceFromIdealSpanningUMISeqLength"] = (
        read_df["SpanningUMISeqLength"].sub(ideal_spanning_umi_seq_len).abs()
    )

    criteria_cols_and_optimization_funcs = [
        ("DistanceFromIdealSpanningUMISeqLength", min),
        ("%BTRBarcodeCoverage_PCR", max),
        ("%BTRBarcodeIdentity_PCR", max),
        ("NumOfBTRBarcodeGapOpenings_PCR", min),
        ("BTRGaps_PCR", min),
        ("BTRAlignmentLength_PCR", max)
    ]
    optimized_cols = []
    is_optimized_cols = []
    for col, optimization_func in criteria_cols_and_optimization_funcs:
        optimized_col = f"{col}_best"
        is_optimized_col = f"{col}_best_is_optimized"
        optimized_cols.append(optimized_col)
        is_optimized_cols.append(is_optimized_col)
        read_df[optimized_col] = optimization_func(read_df[col])
        read_df[is_optimized_col] = read_df[col].eq(read_df[optimized_col])

    read_df["NumOfOptimizedCols"] = read_df[is_optimized_cols].sum(axis=1)

    # criteria_cols = []
    # ascending_sortings = []
    # for col, optimization_func in criteria_cols_and_optimization_funcs:
    #     criteria_cols.append(optimized_col)
    #     ascending_sortings.append(True if optimization_func == max else False)
    criteria_cols = [col for col, _ in criteria_cols_and_optimization_funcs]

    read_df = (
        read_df.loc[
            read_df["NumOfOptimizedCols"].eq(read_df["NumOfOptimizedCols"].max())
        ]
        # .sort_values(by=criteria_cols, ascending=ascending_sortings)
        .drop(
            columns=criteria_cols
            + optimized_cols
            + is_optimized_cols
            + ["NumOfOptimizedCols"]
        )
        # .iloc[0]
        .sample(n=1, random_state=seed)  
    )

    # del read_df["DistanceFromIdealSpanningUMISeqLength"]

    return read_df


# %%
best_gene_specific_pcr_amplified_concat_alignments_df

# %%
best_gene_specific_pcr_amplified_concat_alignments_df["Gene"] + "-" + best_gene_specific_pcr_amplified_concat_alignments_df["Read"]

# %%
# best_gene_specific_pcr_amplified_concat_alignments_df["Read2"] = (
#     best_gene_specific_pcr_amplified_concat_alignments_df["Read"]
# )
best_gene_specific_pcr_amplified_concat_alignments_df["Read2"] = (
    best_gene_specific_pcr_amplified_concat_alignments_df["Gene"] 
    + "-" 
    + best_gene_specific_pcr_amplified_concat_alignments_df["Read"]
)


best_gene_specific_pcr_amplified_concat_alignments_df = (
    best_gene_specific_pcr_amplified_concat_alignments_df.groupby("Read2")
    .apply(choose_best_pcr_btr_barcode_alignment_for_read, seed, expected_spanning_umi_seq_len, include_groups=False)
    # .apply(choose_best_pcr_btr_barcode_alignment_for_read, include_groups=False, seed=seed, ideal_spanning_umi_seq_len=expected_spanning_umi_seq_len)
    # .apply(choose_best_pcr_btr_barcode_alignment_for_read, (seed, expected_spanning_umi_seq_len), include_groups=False)
    .reset_index(drop=True)
)
best_gene_specific_pcr_amplified_concat_alignments_df

# %%
# best_gene_specific_pcr_amplified_concat_alignments_df["Read"].nunique()


best_gene_specific_pcr_amplified_concat_alignments_df.drop_duplicates(["Gene", "Read"]).shape[0]

# %% [markdown]
# ### Save reads with recognized barcodes, before de-duplication

# %%
# reads_with_recognizable_barcodes_dir = Path(
#     "/private7/projects/Combinatorics/D.pealeii/Alignment/UMILongReads.ReadsWithRecognizableBarcodes.MergedSamples"
# )
# reads_with_recognizable_barcodes_dir.mkdir(parents=True, exist_ok=True)

# %%
reads_with_recognizable_barcodes_out_file = Path(
    merged_bams_dir,
    "ReadsWithRecognizableBarcodes.tsv",
)

best_gene_specific_pcr_amplified_concat_alignments_df.loc[
    :, ["Sample", "Gene", "Repeat", "Read", "OldRead"]
].to_csv(
    reads_with_recognizable_barcodes_out_file,
    sep="\t",
    index=False,
    # na_rep="NA",
    # float_format="%.2f",
)

# %%
# # Get the set of unique read names to keep
# reads_to_keep = set(
#     best_gene_specific_pcr_amplified_concat_alignments_df["Read"]
# )

for in_bam_path, gene in zip(mapped_merged_bam_files, genes):
    
    # Get the set of unique read names to keep
    reads_to_keep = set(
        best_gene_specific_pcr_amplified_concat_alignments_df.loc[
            best_gene_specific_pcr_amplified_concat_alignments_df["Gene"] == gene,
            "OldRead"
        ]
    )
    out_bam_path = Path(reads_with_recognizable_barcodes_dir, in_bam_path.name)

    with pysam.AlignmentFile(in_bam_path, "rb") as in_bam, pysam.AlignmentFile(
        out_bam_path, "wb", template=in_bam
    ) as out_bam:
        for read in in_bam:
            if read.query_name in reads_to_keep:
                out_bam.write(read)

    print(f"Filtered BAM written to: {out_bam_path}")
    
# !samtools index -M {reads_with_recognizable_barcodes_dir}/*.bam

# %% [markdown]
# ## Find unique reads by UMI sub-seqs

# %% [markdown]
# ### Graph preleminaries

# %%
def compute_reads_with_indistinguishable_umi_subseqs(df):
    """
    For each read, find all other reads that share at least one UMI sub-sequence.
    Returns a DataFrame with new columns:
      - OtherReadswithIndistinguishableUMISubSeqs
      - NumOfOtherReadswithIndistinguishableUMISubSeqs
    """
    # Build a mapping from each sub-sequence to the set of reads containing it
    subseq_to_reads = {}
    for _, (read, subseqs) in df[
        ["Read", "UMIUniqueSubSeqs"]
    ].iterrows():  # _ is the row index
        for subseq in subseqs:
            subseq_to_reads.setdefault(subseq, set()).add(read)

    # For each read, collect all reads sharing any sub-sequence
    reads_with_overlap = []
    for _, (read, subseqs) in df[
        ["Read", "UMIUniqueSubSeqs"]
    ].iterrows():  # _ is the row index
        overlapping_reads = set()
        for subseq in subseqs:
            overlapping_reads.update(subseq_to_reads.get(subseq, set()))
        # reads_with_overlap.append(list(overlapping_reads))
        reads_with_overlap.append(list(overlapping_reads - {read}))

    df = df.copy()
    # df["ReadswithIndistinguishableUMISubSeqs"] = reads_with_overlap
    # df["OtherReadswithIndistinguishableUMISubSeqs"] = [
    #     [r for r in reads if r != read]
    #     for read, reads in zip(df["Read"], df["ReadswithIndistinguishableUMISubSeqs"])
    # ]
    df["OtherReadswithIndistinguishableUMISubSeqs"] = reads_with_overlap
    df["NumOfOtherReadswithIndistinguishableUMISubSeqs"] = df[
        "OtherReadswithIndistinguishableUMISubSeqs"
    ].apply(len)
    return df


# %%
umi_sub_seq_len

# %%
gene_specific_pcr_amplified_dfs = []

for gene, repeat in product(genes, list("123")):

    ic(gene, repeat)

    gene_and_repeat_df = best_gene_specific_pcr_amplified_concat_alignments_df.loc[
        (best_gene_specific_pcr_amplified_concat_alignments_df["Gene"] == gene)
        & (best_gene_specific_pcr_amplified_concat_alignments_df["Repeat"] == repeat)
    ].copy()
    
    gene_and_repeat_df["UMIUniqueSubSeqs"] = (
        gene_and_repeat_df["SpanningUMISeq"].apply(
            lambda x: split_umi_seq_to_unique_sub_seqs(x, umi_sub_seq_len)
        )
    )

    gene_specific_pcr_amplified_dfs.append(gene_and_repeat_df)

with Pool(processes=6) as pool:
    processed_gene_specific_pcr_amplified_dfs = pool.map(
        # func=process_gene_and_repeat_df,
        func=compute_reads_with_indistinguishable_umi_subseqs,
        iterable=gene_specific_pcr_amplified_dfs,
    )

best_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df = pd.concat(
    processed_gene_specific_pcr_amplified_dfs, ignore_index=True
)
best_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df

# %%
best_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df.groupby(
    "Sample"
).size()

# %%
best_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df.loc[
    best_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df[
        "NumOfOtherReadswithIndistinguishableUMISubSeqs"
    ].eq(0),
].groupby("Sample").size()

# %%
best_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df.loc[
    best_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df[
        "NumOfOtherReadswithIndistinguishableUMISubSeqs"
    ].eq(0),
].groupby("Sample").size().mul(100).div(
    best_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df.groupby(
        "Sample"
    ).size()
).round(2)

# %%
best_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df[
    "NumOfOtherReadswithIndistinguishableUMISubSeqs"
].eq(0).sum()

# %%
best_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df.loc[
    best_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df[
        "NumOfOtherReadswithIndistinguishableUMISubSeqs"
    ].eq(0),
].shape[0] * 100 / best_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df.shape[0]


# %% [markdown]
# ### Analyze graph structure

# %%
def create_indistinguishable_graph(gene_and_repeat_df):
    # create g by explicitly deepcopying the needed cols
    reads = gene_and_repeat_df["Read"].tolist()
    other_indistinguishable_reads = (
        gene_and_repeat_df["OtherReadswithIndistinguishableUMISubSeqs"]
        .apply(lambda x: copy.deepcopy(x))
        .tolist()
    )

    G = nx.Graph()

    G.add_nodes_from(reads)
    # ic(G.number_of_nodes())

    for read, neighbours in zip(reads, other_indistinguishable_reads):
        for neighbour in neighbours:
            G.add_edge(read, neighbour)
            
    # ic(G.number_of_edges());
    return G


# %%
# Clear any previous outputs and run the loop cleanly
clear_output()

print("Starting loop execution...")
print("Expected iterations:", len(list(product(genes, list("123")))))
print("-" * 40)

Gs = []
ccs_dfs = []
for i, (gene, repeat) in enumerate(product(genes, list("123"))):
    # i_s.append(i)
    print(f"Iteration {i}: gene={gene}, repeat={repeat}")

    gene_and_repeat_df = best_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df.loc[
        (best_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df["Gene"] == gene)
        & (best_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df["Repeat"] == repeat),
        ["Read", "OtherReadswithIndistinguishableUMISubSeqs"]
    ]

    G = create_indistinguishable_graph(gene_and_repeat_df)
    Gs.append(G)
    
    ccs_df = pd.DataFrame(
        {
                "Gene": gene,
                "Repeat": int(repeat),
                "CC":  list(nx.connected_components(G)),
        }
    )
    ccs_df["Size"] = ccs_df["CC"].apply(len)
    ccs_df["Degrees"] = ccs_df["CC"].apply(lambda x: [len(G.adj[read]) for read in x])
    ccs_df["MeanDegree"] = ccs_df.apply(lambda x: sum(x["Degrees"]) / x["Size"], axis=1)
    ccs_dfs.append(ccs_df)

print("-" * 40)
print(f"Loop completed. Processed {i} iterations, created {len(Gs)} graphs.")

concat_ccs_df = pd.concat(ccs_dfs, ignore_index=True)

assert concat_ccs_df.loc[
    concat_ccs_df["Degrees"].apply(max).ge(concat_ccs_df["Size"])
].empty, "Max degree in connected components should not exceed the size of the component."

concat_ccs_df["Edges"] = concat_ccs_df["Degrees"].apply(
    lambda x: sum(x) / 2
)
concat_ccs_df["Cliquishness"] = concat_ccs_df.apply(
    # lambda x: x["Edges"] / ((x["Size"] * (x["Size"] - 1)) / 2) if x["Size"] > 1 else np.nan,
    lambda x: x["Edges"] / ((x["Size"] * (x["Size"] - 1)) / 2) if x["Size"] > 1 else 1,
    axis=1
)

concat_ccs_df

# %%

# %%
# how many reads are in connected components of size >= 2 that are not cliquish in
concat_ccs_df.loc[
    (concat_ccs_df["Size"] >= 1)
    & (concat_ccs_df["Cliquishness"] < 1),    
].groupby(["Gene", "Repeat"])["Size"].sum()

# %%
# num of nodes per graph
concat_ccs_df.groupby(["Gene", "Repeat"])["Size"].apply(np.sum).astype(int).reset_index()

# %%
# num of edges per graph
concat_ccs_df.groupby(["Gene", "Repeat"])["Edges"].apply(np.sum).astype(int).reset_index()

# %%
# num of connected components per graph
concat_ccs_df[["Gene", "Repeat"]].value_counts().reset_index().sort_values(["Gene", "Repeat"]).reset_index(drop=True)

# %%
# connected components sizes stats
concat_ccs_df.groupby(["Gene", "Repeat"])["Size"].describe().round(2).reset_index()

# %%
# % of connected components with size 1
concat_ccs_df.groupby(["Gene", "Repeat"])["Size"].apply(lambda x: x.eq(1).sum() * 100 / x.size).round(2)

# %%
# % of nodes in connected components with size 1
concat_ccs_df.groupby(["Gene", "Repeat"]).apply(
    lambda x: x.loc[x["Size"].eq(1), "Size"].sum() * 100 / x["Size"].sum(),
    include_groups=False
).round(2)

# %%
fig = px.histogram(
    concat_ccs_df,
    x="Size",
    # y="MeanDegree",
    # facet_col="Gene",
    # facet_row="Repeat",
    facet_col="Repeat",
    facet_row="Gene",
    # color="Repeat",
    # opacity=0.6,
    # log_x=True,
    log_y=True,
    histnorm="percent",
    # cumulative=True,
    labels={"Size": "Connected component size"},
)
fig.update_layout(
    width=1200,
    height=600,
    # barmode='overlay'
    title="Connected components size distribution"
)
fig.show()

# %%
fig = px.histogram(
    concat_ccs_df.loc[concat_ccs_df["Size"].ge(2)],
    x="MeanDegree",
    # y="MeanDegree",
    # facet_col="Gene",
    # facet_row="Repeat",
    facet_col="Repeat",
    facet_row="Gene",
    # color="Repeat",
    # opacity=0.6,
    # log_x=True,
    log_y=True,
    # histnorm="percent",
    # cumulative=True,
    labels={"MeanDegree": "Connected component mean degree"},
)
fig.update_layout(
    width=1200,
    height=600,
    # barmode='overlay',
    title="Mean degree distribution per connected components with size >= 2"
)
fig.show()

# %%
fig = px.scatter(
    concat_ccs_df.groupby(["Gene", "Repeat", "Size"])["MeanDegree"].agg(["mean", "std"]).reset_index(),
    x="Size",
    y="mean",
    error_y="std",
    facet_col="Repeat",
    facet_row="Gene",
    # color="Repeat",
    # opacity=0.6,
    # log_x=True,
    # log_y=True,
    # histnorm="percent",
    # cumulative=True,
    labels={
        "Size": "Connected component size", 
        # "mean": "Connected component mean degree"
        "mean": "Mean degree"
        },
)
fig.update_xaxes(tick0=0, dtick=5)
fig.update_yaxes(tick0=0, dtick=5)
fig.update_layout(
    width=1200,
    height=400,
    # barmode='overlay',
    # title="Mean degree distribution per connected components with size >= 2"
)
fig.show()

# %%
fig = px.scatter(
    concat_ccs_df.loc[
        concat_ccs_df["Size"].ge(2)
    ].groupby(["Gene", "Repeat", "Size"])["Cliquishness"].agg(["mean", "std"]).reset_index(),
    x="Size",
    y="mean",
    error_y="std",
    facet_col="Repeat",
    facet_row="Gene",
    # color="Repeat",
    # opacity=0.6,
    # log_x=True,
    # log_y=True,
    # histnorm="percent",
    # cumulative=True,
    labels={
        "Size": "Connected component size", 
        # "mean": "Connected component mean degree"
        "mean": "Mean cliquishness<br>(|E| / |E|_clique)"
        # "mean": "Mean cliquishness<br>(|E| / |V|*(|V|-1)/2)"
        },
)
fig.update_layout(
    width=1200,
    height=600,
    # barmode='overlay',
    title="How close is a connected component of size >= 2 to be a clique?"
)
fig.show()

# %%
concat_ccs_df

# %%
ccs_df = concat_ccs_df.loc[
    (concat_ccs_df["Gene"] == "ADAR1")
    & (concat_ccs_df["Repeat"] == 1)
]
ccs_df

# %%
G = Gs[0]

# %%
max_cliques = list(nx.find_cliques(G))
len(max_cliques)

# %%
for i, (gene, repeat) in enumerate(product(genes, list("123"))):
    G = Gs[i]
    max_cliques = list(nx.find_cliques(G))
    print(f"Gene: {gene}, Repeat: {repeat}, Max cliques: {len(max_cliques)}")
    # if len(max_cliques) > 0:
    #     print("Example max clique:", max_cliques[0])
    # else:
    #     print("No max cliques found.")
    # print("-" * 40)

# %%
for i, (gene, repeat) in enumerate(product(genes, list("123"))):
    
    print(f"Iteration {i}: gene={gene}, repeat={repeat}")
    
    G = Gs[i]
    
    ccs_df = concat_ccs_df.loc[
        (concat_ccs_df["Gene"] == gene)
        & (concat_ccs_df["Repeat"] == int(repeat))
    ]
    # ccs_df

    gene_and_repeat_df = best_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df.loc[
        (best_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df["Gene"] == gene)
        & (best_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df["Repeat"] == repeat),
        # ["Read", "OtherReadswithIndistinguishableUMISubSeqs"]
    ]
    
    break

gene_and_repeat_df

# %%
ccs_df


# %%
# how can we make more stringent definition of an edge between two reads?
# for each edge (u, v), report:
# 1 - u's spanning barcode length
# 2 - v's spanning barcode length
# 3 - length of the longest common subsequence of u's and v's spanning barcodes

# %%
def make_edges_umi_df(gene, repeat, G, gene_and_repeat_df, max_alignments):
    edges_umi_data = []

    for u, v in G.edges():
        
        u_v_gene_and_repeat_df = gene_and_repeat_df.loc[
            gene_and_repeat_df["Read"].isin([u, v]),
            ["Read", "SpanningUMISeq", "SpanningUMISeqLength",]
        ].set_index("Read")
        
        # u_spanning_umi_seq = u_v_gene_and_repeat_df.loc[u, "SpanningUMISeq"]
        # v_spanning_umi_seq = u_v_gene_and_repeat_df.loc[v, "SpanningUMISeq"]
        u_spanning_umi_seq_len = u_v_gene_and_repeat_df.loc[u, "SpanningUMISeqLength"]
        v_spanning_umi_seq_len = u_v_gene_and_repeat_df.loc[v, "SpanningUMISeqLength"]
        
        if v_spanning_umi_seq_len < u_spanning_umi_seq_len:
            # swap u and v to ensure u has the shorter spanning UMI seq
            u, v = v, u
            u_spanning_umi_seq_len, v_spanning_umi_seq_len = (
                v_spanning_umi_seq_len,
                u_spanning_umi_seq_len,
            )
            
        
        for x, y in [(u, v), (v, u)]:
            x_spanning_umi_seq = u_v_gene_and_repeat_df.loc[x, "SpanningUMISeq"]
            y_spanning_umi_seq = u_v_gene_and_repeat_df.loc[y, "SpanningUMISeq"]
            # x_spanning_umi_seq_len = u_v_gene_and_repeat_df.loc[x, "SpanningUMISeqLength"]
            # y_spanning_umi_seq_len = u_v_gene_and_repeat_df.loc[y, "SpanningUMISeqLength"]
            # ic(x, y)
            
            aligner = Align.PairwiseAligner(
                # match_score=1.0, mismatch_score=-2.0, gap_score = -2.5,
                mode="local",  # otherwise we'd get scattered matches of the primer across the read
                scoring="blastn",
                # scoring="megablast"
            )
            
            alignments = aligner.align(x_spanning_umi_seq, y_spanning_umi_seq, strand="+")

            for alignment_i, alignment in enumerate(alignments, start=1):

                score = alignment.score
                gaps, identities, mismatches = alignment.counts()
                alignment_length = alignment.length

                order = "U, V" if (x, y) == (u, v) else "V, U"
                u_v_umi_data = [
                    u, v,
                    u_spanning_umi_seq_len, v_spanning_umi_seq_len,
                    order,
                    alignment_i, score, gaps, identities, mismatches, alignment_length, alignment
                ]
                
                edges_umi_data.append(u_v_umi_data)
                
                if (
                    max_alignments is not None
                    and i == max_alignments
                ):
                    break

    edges_umi_df = pd.DataFrame(
        {
            "Gene": gene,
            "Repeat": repeat,
            "U": [data[0] for data in edges_umi_data],
            "V": [data[1] for data in edges_umi_data],
            "USpanningUMILength": [data[2] for data in edges_umi_data],
            "VSpanningUMILength": [data[3] for data in edges_umi_data],
            "Order": [data[4] for data in edges_umi_data],
            "AlignmentIndex": [data[5] for data in edges_umi_data],
            "Score": [data[6] for data in edges_umi_data],
            "Gaps": [data[7] for data in edges_umi_data],
            "Identities": [data[8] for data in edges_umi_data],
            "Mismatches": [data[9] for data in edges_umi_data],
            "AlignmentLength": [data[10] for data in edges_umi_data],
            "Alignment": [data[11] for data in edges_umi_data],
        }
    )

    assert edges_umi_df.groupby(["U", "V"])["AlignmentLength"].nunique().describe()["max"] == 1, "There should be only one alignment length per edge (U, V), regardless of the order (U, V) or (V, U)."

    # following the last assertment, we can safely drop the duplicate edges,
    # as they will have the same alignment length and other properties,
    # and also the columns representing the order of the alignment
    # (as we only keep one (u, v) alignment)
    edges_umi_df = edges_umi_df.drop_duplicates(["U", "V"]).drop(columns=["Order", "AlignmentIndex", "Score"])
    
    edges_umi_df.insert(
        edges_umi_df.columns.get_loc("VSpanningUMILength") + 1,
        "SpanningUMISeqAbsDiff",
        edges_umi_df["USpanningUMILength"].sub(edges_umi_df["VSpanningUMILength"]).abs()
    )

    edges_umi_df["%AlignmentLength/USpanningUMILength"] = (
        edges_umi_df["AlignmentLength"].mul(100).div(edges_umi_df["USpanningUMILength"])
    )
    edges_umi_df["%AlignmentLength/VSpanningUMILength"] = (
        edges_umi_df["AlignmentLength"].mul(100).div(edges_umi_df["VSpanningUMILength"])
    )
        
    return edges_umi_df

# %%
max_alignments = 100

with Pool() as pool:
    edges_umi_dfs = pool.starmap(
        func=make_edges_umi_df,
        iterable=[
            (
                gene,
                repeat,
                Gs[i], 
                best_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df.loc[
                    (best_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df["Gene"] == gene)
                    & (best_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df["Repeat"] == repeat),
                    ["Read", "SpanningUMISeq", "SpanningUMISeqLength",]
                ],
                max_alignments
            )
            for i, (gene, repeat) in enumerate(product(genes, list("123")))
        ]
    )

concat_edges_umi_df = pd.concat(edges_umi_dfs, ignore_index=True)

assert concat_edges_umi_df["Gaps"].describe()["max"] == 0
del concat_edges_umi_df["Gaps"]

concat_edges_umi_df

# %%
# concat_edges_umi_df.loc[
#     (concat_edges_umi_df["USpanningUMILength"].eq(concat_edges_umi_df["VSpanningUMILength"]))
#     & (concat_edges_umi_df["USpanningUMILength"].eq(concat_edges_umi_df["AlignmentLength"]))
# ]

# %%
concat_edges_umi_df.loc[
    (concat_edges_umi_df["AlignmentLength"].sub(concat_edges_umi_df["USpanningUMILength"]).abs().le(1))
    & (concat_edges_umi_df["AlignmentLength"].sub(concat_edges_umi_df["VSpanningUMILength"]).abs().le(1))
]

# %%
concat_edges_umi_df.loc[
    (concat_edges_umi_df["AlignmentLength"].sub(concat_edges_umi_df["USpanningUMILength"]).abs().le(1))
    & (concat_edges_umi_df["AlignmentLength"].sub(concat_edges_umi_df["VSpanningUMILength"]).abs().le(1)),
    
].groupby(["Gene", "Repeat"])[["AlignmentLength", "Mismatches"]].describe()

# %%
concat_edges_umi_df["Mismatches"].describe()

# %%
concat_edges_umi_df["SpanningUMISeqAbsDiff"].describe()

# %%
# concat_edges_umi_df.loc[:, ["Gene", "Repeat", "USpanningUMILength", "VSpanningUMILength"]].value_counts()

# %%
fig = px.box(
    concat_edges_umi_df.transform(
        lambda x: x.astype({"USpanningUMILength": "str"})
    ),
    x="USpanningUMILength",
    y="VSpanningUMILength",
    facet_row_spacing=0.15,
    color="Repeat",
    # size="Edges",
    # facet_col="Repeat",
    facet_row="Gene",
    category_orders={
        "USpanningUMILength": [str(x) for x in range(10, 15)],
    },
    labels={
        "USpanningUMILength": "U's length",
        "VSpanningUMILength": "V's length",
        },
)
# fig.update_xaxes(tick0=0, dtick=2)
fig.update_layout(
    width=800,
    height=600,
    # Reduce spacing between groups and boxes
    # boxmode='group',   # group by "x" 
    # boxgap=0.1,        # space between boxes in one group (default ~0.3)
    # boxgroupgap=0.01   # space between groups (different x values)
)
fig.show()

# %%
# fig = px.box(
#     concat_edges_umi_df.transform(
#         lambda x: x.astype({"USpanningUMILength": "str"})
#     ),
#     x="USpanningUMILength",
#     y="VSpanningUMILength",
#     facet_row_spacing=0.15,
#     color="AlignmentLength",
#     # size="Edges",
#     facet_col="Repeat",
#     facet_row="Gene",
#     category_orders={
#         "USpanningUMILength": [str(x) for x in range(10, 15)],
#     },
#     # labels={
#     #     "%AlignmentLength/USpanningUMILength": "% alignment length /<br>U's length",
#     #     "%AlignmentLength/VSpanningUMILength": "% alignment length /<br>V's length",
#     #     "Edges": "% edges in sample"
#     #     },
# )
# # fig.update_xaxes(tick0=0, dtick=2)
# fig.update_layout(
#     width=1400,
#     height=600,
#     # Reduce spacing between groups and boxes
#     # boxmode='group',   # group by "x" 
#     # boxgap=0.1,        # space between boxes in one group (default ~0.3)
#     # boxgroupgap=0.01   # space between groups (different x values)
# )
# fig.show()

# %%
fig = px.box(
    concat_edges_umi_df.transform(
        lambda x: x.astype({"USpanningUMILength": "str"})
    ),
    x="USpanningUMILength",
    y="VSpanningUMILength",
    facet_row_spacing=0.15,
    # color="AlignmentLength",
    # facet_col="Repeat",
    facet_col="AlignmentLength",
    color="Repeat",
    facet_row="Gene",
    category_orders={
        "USpanningUMILength": [str(x) for x in range(10, 15)],
    },
    labels={
        "USpanningUMILength": "U's length",
        "VSpanningUMILength": "V's length",
        "AlignmentLength": "Alignment length"
        },
)
# fig.update_xaxes(tick0=0, dtick=2)
fig.update_layout(
    width=1400,
    height=600,
    # Reduce spacing between groups and boxes
    # boxmode='group',   # group by "x" 
    # boxgap=0.1,        # space between boxes in one group (default ~0.3)
    # boxgroupgap=0.01   # space between groups (different x values)
)
fig.show()

# %%
fig = px.box(
    concat_edges_umi_df.transform(
        lambda x: x.astype({"USpanningUMILength": "str"})
    ),
    x="USpanningUMILength",
    # y="%AlignmentLength/USpanningUMILength",
    y="AlignmentLength",
    # y="VSpanningUMILength",
    facet_row_spacing=0.15,
    # color="AlignmentLength",
    # facet_col="Repeat",
    # facet_col="AlignmentLength",
    # facet_col="VSpanningUMILength",
    color="Repeat",
    facet_row="Gene",
    category_orders={
        "USpanningUMILength": [str(x) for x in range(10, 15)],
    },
    labels={
        "USpanningUMILength": "U's length",
        # "%AlignmentLength/USpanningUMILength": "% alignment length /<br>U's length",
        "AlignmentLength": "Alignment length",
        # "%AlignmentLength/VSpanningUMILength": "% alignment length /<br>V's length",
        # "Edges": "% edges in sample"
        },
)
# fig.update_xaxes(tick0=0, dtick=2)
fig.update_traces(
    line=dict(width=4)  # bolder median line
)
fig.update_layout(
    width=800,
    height=600,
    # Reduce spacing between groups and boxes
    # boxmode='group',   # group by "x" 
    # boxgap=0.1,        # space between boxes in one group (default ~0.3)
    # boxgroupgap=0.01   # space between groups (different x values)
)
fig.show()

# %%
# fig = px.histogram(
#     concat_edges_umi_df,
#     x="USpanningUMILength",
#     y="VSpanningUMILength",
#     facet_col="Repeat",
#     facet_row="Gene",
#     facet_row_spacing=0.1,
#     # color="Repeat",
#     # opacity=0.6,
#     # log_x=True,
#     # log_y=True,
#     # histnorm="percent",
#     histfunc="avg",
#     # cumulative=True,
#     # labels={"MeanDegree": "Connected component mean degree"},
# )
# fig.update_xaxes(dtick=1)
# fig.update_yaxes(dtick=1, range=[10, None])
# fig.update_layout(
#     width=1200,
#     height=600,
#     # barmode='overlay',
#     # title="Mean degree distribution per connected components with size >= 2"
# )
# fig.show()

# %%
fig = px.histogram(
    concat_edges_umi_df,
    x="SpanningUMISeqAbsDiff",
    # y="MeanDegree",
    facet_col="Repeat",
    facet_row="Gene",
    # color="Repeat",
    # opacity=0.6,
    # log_x=True,
    # log_y=True,
    # histnorm="percent",
    # cumulative=True,
    # labels={"MeanDegree": "Connected component mean degree"},
)
fig.update_xaxes(dtick=1)
fig.update_layout(
    width=800,
    height=400,
    # barmode='overlay',
    # title="Mean degree distribution per connected components with size >= 2"
)
fig.show()

# %%
fig = px.histogram(
    concat_edges_umi_df,
    x="SpanningUMISeqAbsDiff",
    # y="MeanDegree",
    facet_col="Repeat",
    facet_row="Gene",
    # color="Repeat",
    # opacity=0.6,
    # log_x=True,
    # log_y=True,
    histnorm="percent",
    # cumulative=True,
    # labels={"MeanDegree": "Connected component mean degree"},
)
fig.update_xaxes(dtick=1)
fig.update_layout(
    width=800,
    height=400,
    # barmode='overlay',
    # title="Mean degree distribution per connected components with size >= 2"
)
fig.show()

# %%
fig = px.histogram(
    concat_edges_umi_df,
    x="AlignmentLength",
    # y="MeanDegree",
    facet_col="Repeat",
    facet_row="Gene",
    facet_row_spacing=0.1,
    # color="Repeat",
    # opacity=0.6,
    # log_x=True,
    # log_y=True,
    # histnorm="percent",
    # cumulative=True,
    # labels={"MeanDegree": "Connected component mean degree"},
)
fig.update_xaxes(dtick=1)
fig.update_layout(
    width=800,
    height=400,
    # barmode='overlay',
    # title="Mean degree distribution per connected components with size >= 2"
)
fig.show()

# %%
fig = px.histogram(
    concat_edges_umi_df,
    x="AlignmentLength",
    # y="MeanDegree",
    facet_col="Repeat",
    facet_row="Gene",
    facet_row_spacing=0.1,
    # color="Repeat",
    # opacity=0.6,
    # log_x=True,
    # log_y=True,
    histnorm="percent",
    # cumulative=True,
    # labels={"MeanDegree": "Connected component mean degree"},
)
fig.update_xaxes(dtick=1)
fig.update_layout(
    width=800,
    height=400,
    # barmode='overlay',
    # title="Mean degree distribution per connected components with size >= 2"
)
fig.show()

# %%
# fig = px.histogram(
#     concat_edges_umi_df,
#     x="USpanningUMILength",
#     y="AlignmentLength",
#     facet_col="Repeat",
#     facet_row="Gene",
#     facet_row_spacing=0.1,
#     # color="Repeat",
#     # opacity=0.6,
#     # log_x=True,
#     # log_y=True,
#     # histnorm="percent",
#     histfunc="avg",
#     # cumulative=True,
#     # labels={"MeanDegree": "Connected component mean degree"},
# )
# fig.update_xaxes(dtick=1)
# fig.update_yaxes(
#     # dtick=1, 
#                  range=[concat_edges_umi_df["AlignmentLength"].min(), None]
#                  )
# fig.update_layout(
#     width=1200,
#     height=600,
#     # barmode='overlay',
#     # title="Mean degree distribution per connected components with size >= 2"
# )
# fig.show()

# %%
# # APPROACH 5: Build plot trace by trace in correct order


# df = concat_edges_umi_df.assign(
#     Sample = concat_edges_umi_df["Gene"] + "-" + concat_edges_umi_df["Repeat"]
# ).loc[
#     :, ["Sample", "%AlignmentLength/USpanningUMILength", "%AlignmentLength/VSpanningUMILength", "SpanningUMISeqAbsDiff"]
# ].groupby("Sample").value_counts(normalize=True).mul(100).div(10).astype(int).apply(lambda x: f"{x*10}-{x*10+10}").reset_index(name="Edges")

# # Get unique samples and spanning diff values for subplot layout
# samples = [f"{gene}-{repeat}" for gene, repeat in product(genes, list("123"))]
# spanning_diffs = sorted(df["SpanningUMISeqAbsDiff"].unique())

# # Create subplots
# fig = make_subplots(
#     rows=len(spanning_diffs), 
#     cols=len(samples),
#     subplot_titles=samples,
#     vertical_spacing=0.02,
#     horizontal_spacing=0.02,
#     row_titles=[f"UMI abs diff = {diff}" for diff in spanning_diffs],
#     shared_xaxes="all",
#     shared_yaxes="all",
#     x_title="% alignment length / U's length",      
#     y_title="% alignment length / V's length",
    
# )

# # Define color mapping for consistent colors
# desired_order = ["0-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-90", "90-100"]
# colors = px.colors.qualitative.Set1
# # colors = px.colors.sequential.Plasma_r
# # colors = px.colors.sequential.Rainbow
# # colors = px.colors.qualitative.G10
# colors = px.colors.qualitative.Vivid
# color_map = {edge: colors[i % len(colors)] for i, edge in enumerate(desired_order)}

# # print("Color mapping:", color_map)

# # Track which edges we've added to legend
# legend_added = set()

# # Add traces in the correct order
# for edge in desired_order:
#     if edge in df["Edges"].values:
#         edge_data = df[df["Edges"] == edge]
        
#         for _, row in edge_data.iterrows():
#             sample_idx = samples.index(row["Sample"]) + 1
#             diff_idx = spanning_diffs.index(row["SpanningUMISeqAbsDiff"]) + 1
            
#             fig.add_trace(
#                 go.Scatter(
#                     x=[row["%AlignmentLength/USpanningUMILength"]],
#                     y=[row["%AlignmentLength/VSpanningUMILength"]],
#                     mode='markers',
#                     marker=dict(color=color_map[edge], size=8),
#                     name=edge,
#                     legendgroup=edge,
#                     showlegend=edge not in legend_added,
#                     # opacity=0.6
#                 ),
#                 row=diff_idx,
#                 col=sample_idx
#             )
            
#             legend_added.add(edge)

# fig.update_layout(
#     width=1200,
#     height=1200*11/14,
#     # title="Edges by Sample with Correctly Ordered Legend",
#     legend_title_text="% of edges in sample"
    
# )

# # print("Legend order should be:", [edge for edge in desired_order if edge in legend_added])
# fig.show()

# %%
# for gene, repeat in product(genes, list("123")):

#     fig = px.scatter(
#         concat_edges_umi_df.loc[
#             (concat_edges_umi_df["Gene"] == gene) 
#             & (concat_edges_umi_df["Repeat"] == repeat),
#              ["%AlignmentLength/USpanningUMILength", "%AlignmentLength/VSpanningUMILength", "SpanningUMISeqAbsDiff"]
#             ].value_counts().reset_index().rename(columns={"count": "Edges"}),
#         x="%AlignmentLength/USpanningUMILength",
#         y="%AlignmentLength/VSpanningUMILength",
#         # y="MeanDegree",
#         # facet_col="Repeat",
#         # facet_col="Repeat",
#         # facet_row="Gene",
#         # facet_row_spacing=0.1,
#         color="Edges",
#         facet_col="SpanningUMISeqAbsDiff",
#         # opacity=0.6,
#         # log_x=True,
#         # log_y=True,
#         # histnorm="percent",
#         # cumulative=True,
#         labels={
#             "%AlignmentLength/USpanningUMILength": "% alignment length /<br>U's length",
#             "%AlignmentLength/VSpanningUMILength": "% alignment length /<br>V's length",
#             },
#     )
#     # fig.update_xaxes(tick0=0, dtick=2)
#     fig.update_layout(
#         width=1600,
#         height=1600/5,
#         # barmode='overlay',
#         title=f"{gene} - {repeat}",
#     )
#     fig.show()

# %%

# %%

# %%
concat_ccs_df.loc[
    concat_ccs_df["Size"].ge(8),
    ["Gene", "Repeat"]
].value_counts()

# %%
sample_ccs_df = concat_ccs_df.loc[
    concat_ccs_df["Size"].ge(8),
].groupby(["Gene", "Repeat"]).sample(
    n=1, random_state=seed, replace=False
).set_index(["Gene", "Repeat"])
sample_ccs_df

# %%
nodes_primers_dfs = []

for i, (gene, repeat) in enumerate(product(genes, list("123"))):

    nodes = sample_ccs_df.loc[(gene, int(repeat)), "CC"]
    
    print(f"Iteration {i}: gene={gene}, repeat={repeat}, {len(nodes)} nodes")
    
    nodes_primers_df = best_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df.loc[
        (best_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df["Gene"] == gene)
        & (best_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df["Repeat"] == repeat)
        & (best_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df["Read"].isin(nodes))
    ]
    
    nodes_primers_dfs.append(nodes_primers_df)

# %%
primers_dict

# %%
nodes_primers_df = nodes_primers_dfs[0]
nodes_primers_df

# %%
for _, (read, read_seq, btr_read_start) in nodes_primers_df.loc[:, ["Read", "ReadSeq", "BTRReadStart"]].iterrows():
    print(f">{read}\n{read_seq[btr_read_start-10:]}")

# %%
for _, (read, spanning_umi_seq) in nodes_primers_df.loc[:, ["Read", "SpanningUMISeq"]].iterrows():
    print(f">{read}\n{spanning_umi_seq}")

# %%
primers_dict["ADAR1"].reverse_complement()

# %%
primers_dict["PCR"].reverse_complement()

# %% [markdown]
# ![image-2.png](attachment:image-2.png)

# %% [markdown]
# ![image.png](attachment:image.png)

# %%
nodes_primers_df = nodes_primers_dfs[1]
nodes_primers_df

# %%
for _, (read, read_seq, btr_read_start) in nodes_primers_df.loc[:, ["Read", "ReadSeq", "BTRReadStart"]].iterrows():
    print(f">{read}\n{read_seq[btr_read_start-10:]}")

# %%
for _, (read, spanning_umi_seq) in nodes_primers_df.loc[:, ["Read", "SpanningUMISeq"]].iterrows():
    print(f">{read}\n{spanning_umi_seq}")

# %%
primers_dict["ADAR1"].reverse_complement()

# %%
primers_dict["PCR"].reverse_complement()

# %% [markdown]
# ![image.png](attachment:image.png)

# %% [markdown]
# ![image.png](attachment:image.png)

# %%
nodes_primers_df = nodes_primers_dfs[2]
nodes_primers_df

# %%
for _, (read, read_seq, btr_read_start) in nodes_primers_df.loc[:, ["Read", "ReadSeq", "BTRReadStart"]].iterrows():
    print(f">{read}\n{read_seq[btr_read_start-10:]}")

# %%
for _, (read, spanning_umi_seq) in nodes_primers_df.loc[:, ["Read", "SpanningUMISeq"]].iterrows():
    print(f">{read}\n{spanning_umi_seq}")

# %%
primers_dict["ADAR1"].reverse_complement()

# %%
primers_dict["PCR"].reverse_complement()

# %% [markdown]
#

# %% [markdown]
#

# %%
nodes_primers_df = nodes_primers_dfs[3]
nodes_primers_df

# %%
for _, (read, read_seq, btr_read_start) in nodes_primers_df.loc[:, ["Read", "ReadSeq", "BTRReadStart"]].iterrows():
    print(f">{read}\n{read_seq[btr_read_start-10:]}")

# %%
for _, (read, spanning_umi_seq) in nodes_primers_df.loc[:, ["Read", "SpanningUMISeq"]].iterrows():
    print(f">{read}\n{spanning_umi_seq}")

# %%
primers_dict["IQEC"].reverse_complement()

# %%
primers_dict["PCR"].reverse_complement()

# %% [markdown]
# ![image.png](attachment:image.png)

# %% [markdown]
# ![image.png](attachment:image.png)

# %%
nodes_primers_df = nodes_primers_dfs[4]
nodes_primers_df

# %%
for _, (read, read_seq, btr_read_start) in nodes_primers_df.loc[:, ["Read", "ReadSeq", "BTRReadStart"]].iterrows():
    print(f">{read}\n{read_seq[btr_read_start-10:]}")

# %%
for _, (read, spanning_umi_seq) in nodes_primers_df.loc[:, ["Read", "SpanningUMISeq"]].iterrows():
    print(f">{read}\n{spanning_umi_seq}")

# %%
primers_dict["IQEC"].reverse_complement()

# %%
primers_dict["PCR"].reverse_complement()

# %% [markdown]
# ![image.png](attachment:image.png)

# %% [markdown]
# ![image.png](attachment:image.png)

# %%
nodes_primers_df = nodes_primers_dfs[5]
nodes_primers_df

# %%
for _, (read, read_seq, btr_read_start) in nodes_primers_df.loc[:, ["Read", "ReadSeq", "BTRReadStart"]].iterrows():
    print(f">{read}\n{read_seq[btr_read_start-10:]}")

# %%
for _, (read, spanning_umi_seq) in nodes_primers_df.loc[:, ["Read", "SpanningUMISeq"]].iterrows():
    print(f">{read}\n{spanning_umi_seq}")

# %%
primers_dict["IQEC"].reverse_complement()

# %%
primers_dict["PCR"].reverse_complement()


# %%

# %%

# %%

# %% [markdown]
# ### Obtain MIS

# %%
def choose_mis_of_reads_with_distinct_umi_sub_seqs(df, seed):
    # intialize a random number generator with the given seed
    # to select a random high-degree read at each step
    rng = np.random.default_rng(seed)

    # create g by explicitly deepcopying the needed cols
    reads = df["Read"].tolist()
    other_indistinguishable_reads = (
        df["OtherReadswithIndistinguishableUMISubSeqs"]
        .apply(lambda x: copy.deepcopy(x))
        .tolist()
    )
    g = pd.DataFrame(
        {
            "Read": reads,
            "OtherReadswithIndistinguishableUMISubSeqs": other_indistinguishable_reads,
        }
    ).set_index("Read")
    g["Degree"] = g["OtherReadswithIndistinguishableUMISubSeqs"].apply(len)

    # at each step, we choose a read with the maximum degree,
    # remove it from the graph, and update the degrees of its neighbors
    # we continue until there are no reads left with a degree greater than 0
    # this will ensure that we retain only a subset of reads that are distinct from each other
    while (max_degree := g["Degree"].max()) > 0:
        reads_with_max_degree = g.loc[g["Degree"].eq(max_degree)].index
        read = rng.choice(reads_with_max_degree)
        indistinguishable_reads_for_read = g.loc[
            read, "OtherReadswithIndistinguishableUMISubSeqs"
        ]

        # Remove the read from the indistinguishable reads' neighbors
        try:
            g.loc[
                indistinguishable_reads_for_read,
                "OtherReadswithIndistinguishableUMISubSeqs",
            ].apply(lambda x: x.remove(read))
        except KeyError as e:
            ic(reads_with_max_degree, read, indistinguishable_reads_for_read)
            break

        # Following that, decrease the degree of the indistinguishable reads
        g.loc[indistinguishable_reads_for_read, "Degree"] -= 1

        # Remove the read from the graph
        g.drop(read, inplace=True)

    return g.index


# copy.deepcopy()

# %%
distinct_umi_sub_seq_gene_specific_pcr_amplified_alignments_dfs = []

for gene, repeat in product(genes, list("123")):

    ic(gene, repeat)

    gene_and_repeat_df = (
        best_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df.loc[
            (
                best_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df[
                    "Gene"
                ]
                == gene
            )
            & (
                best_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df[
                    "Repeat"
                ]
                == repeat
            )
        ]
    )
    mis_reads = choose_mis_of_reads_with_distinct_umi_sub_seqs(gene_and_repeat_df, seed)
    gene_and_repeat_df = gene_and_repeat_df.loc[
        gene_and_repeat_df["Read"].isin(mis_reads)
    ]

    distinct_umi_sub_seq_gene_specific_pcr_amplified_alignments_dfs.append(
        gene_and_repeat_df
    )


distinct_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df = pd.concat(
    distinct_umi_sub_seq_gene_specific_pcr_amplified_alignments_dfs, ignore_index=True
)
distinct_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df

# %%
distinct_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df.groupby(
    "Sample"
).size()

# %%
distinct_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df.groupby(
    "Sample"
).size()

# %%
ccs_vs_distinct_reads_df = concat_ccs_df[["Gene", "Repeat"]].value_counts().reset_index().sort_values(
        ["Gene", "Repeat"]
    ).reset_index(drop=True).rename(columns={"count": "ConnectedComponents"}).transform(
    lambda x: x.astype({"Repeat": "str"})
    ).merge(
    distinct_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df.groupby(
        ["Gene", "Repeat"]
    ).size().reset_index().rename(columns={0: "DistinctReads"}),
    on=["Gene", "Repeat"],
    how="left",
)
    
ccs_vs_distinct_reads_df["DistinctReads/ConnectedComponent"] = ccs_vs_distinct_reads_df["DistinctReads"].div(
    ccs_vs_distinct_reads_df["ConnectedComponents"]
).round(2)

# assert ccs_vs_distinct_reads_df["DistinctReads"].ge(
#     ccs_vs_distinct_reads_df["ConnectedComponents"]
# ).all(), "Each connected component should have at least one distinct read in it."

assert ccs_vs_distinct_reads_df["DistinctReads/ConnectedComponent"].ge(1).all(), "Each connected component should have at least one distinct read in it."

ccs_vs_distinct_reads_df

# %%
ccs_vs_distinct_reads_df["DistinctReads"].div(
    ccs_vs_distinct_reads_df["ConnectedComponents"]
).round(2)

# %%
# num of connected components per graph
concat_ccs_df[["Gene", "Repeat"]].value_counts().reset_index().sort_values(
        ["Gene", "Repeat"]
    ).reset_index(drop=True).rename(columns={"count": "ConnectedComponents"}).transform(
    lambda x: x.astype({"Gene": "str", "Repeat": "str"})
    )

# %%
distinct_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df.groupby(
    "Gene"
).size()

# %% [markdown]
# ### Save reads after de-duplication

# %%
unique_reads_out_file = Path(
    merged_bams_dir,
    "UniqueReadsByUMISubSeq.tsv",
)

distinct_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df.loc[
    :, ["Sample", "Gene", "Repeat", "Read"]
].to_csv(
    unique_reads_out_file,
    sep="\t",
    index=False,
    # na_rep="NA",
    # float_format="%.2f",
)

# %%
unique_reads_dir = Path(
    "/private7/projects/Combinatorics/D.pealeii/Alignment/UMILongReads.UniqueReadsByUMISubSeq.MergedSamples"
)
unique_reads_dir.mkdir(parents=True, exist_ok=True)

# %%
# Get the set of unique read names to keep
reads_to_keep = set(
    distinct_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df["Read"]
)

for bam_path in mapped_merged_bam_files:
    out_bam_path = Path(unique_reads_dir, Path(bam_path).name)

    with pysam.AlignmentFile(bam_path, "rb") as in_bam, pysam.AlignmentFile(
        out_bam_path, "wb", template=in_bam
    ) as out_bam:
        for read in in_bam:
            if read.query_name in reads_to_keep:
                out_bam.write(read)

    print(f"Filtered BAM written to: {out_bam_path}")


# %%
# !samtools index -M {unique_reads_dir}/*.bam

# %% [markdown]
# ## Find unique reads by overlapping UMI seq

# %% [markdown]
# ### Graph preleminaries

# %%
def split_umi_seq_to_unique_sub_seqs(umi_seq, min_umi_seq_len):
    return {
        umi_seq[x : x + min_umi_seq_len]
        for x in range(0, len(umi_seq) - min_umi_seq_len + 1)
    }


# %%
def compute_reads_with_indistinguishable_umi_subseqs(df):
    """
    For each read, find all other reads that share at least one UMI sub-sequence.
    Returns a DataFrame with new columns:
      - OtherReadswithIndistinguishableUMISubSeqs
      - NumOfOtherReadswithIndistinguishableUMISubSeqs
    """
    # Build a mapping from each sub-sequence to the set of reads containing it
    subseq_to_reads = {}
    for _, (read, subseqs) in df[
        ["Read", "UMIUniqueSubSeqs"]
    ].iterrows():  # _ is the row index
        for subseq in subseqs:
            subseq_to_reads.setdefault(subseq, set()).add(read)

    # For each read, collect all reads sharing any sub-sequence
    reads_with_overlap = []
    for _, (read, subseqs) in df[
        ["Read", "UMIUniqueSubSeqs"]
    ].iterrows():  # _ is the row index
        overlapping_reads = set()
        for subseq in subseqs:
            overlapping_reads.update(subseq_to_reads.get(subseq, set()))
        # reads_with_overlap.append(list(overlapping_reads))
        reads_with_overlap.append(list(overlapping_reads - {read}))

    df = df.copy()
    # df["ReadswithIndistinguishableUMISubSeqs"] = reads_with_overlap
    # df["OtherReadswithIndistinguishableUMISubSeqs"] = [
    #     [r for r in reads if r != read]
    #     for read, reads in zip(df["Read"], df["ReadswithIndistinguishableUMISubSeqs"])
    # ]
    df["OtherReadswithIndistinguishableUMISubSeqs"] = reads_with_overlap
    df["NumOfOtherReadswithIndistinguishableUMISubSeqs"] = df[
        "OtherReadswithIndistinguishableUMISubSeqs"
    ].apply(len)
    return df


# %%
partial_umi_sub_seq_len = 5

gene_specific_pcr_amplified_dfs = []

for gene, repeat in product(genes, list("123")):

    ic(gene, repeat)

    gene_and_repeat_df = best_gene_specific_pcr_amplified_concat_alignments_df.loc[
        (best_gene_specific_pcr_amplified_concat_alignments_df["Gene"] == gene)
        & (best_gene_specific_pcr_amplified_concat_alignments_df["Repeat"] == repeat)
    ].copy()
    
    gene_and_repeat_df["UMIUniqueSubSeqs"] = (
        gene_and_repeat_df["SpanningUMISeq"].apply(
            lambda x: split_umi_seq_to_unique_sub_seqs(x, partial_umi_sub_seq_len)
        )
    )

    gene_specific_pcr_amplified_dfs.append(gene_and_repeat_df)

with Pool(processes=6) as pool:
    processed_gene_specific_pcr_amplified_dfs = pool.map(
        # func=process_gene_and_repeat_df,
        func=compute_reads_with_indistinguishable_umi_subseqs,
        iterable=gene_specific_pcr_amplified_dfs,
    )

best_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df = pd.concat(
    processed_gene_specific_pcr_amplified_dfs, ignore_index=True
)
best_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df

# %%
best_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df.groupby(
    "Sample"
)["NumOfOtherReadswithIndistinguishableUMISubSeqs"].describe().round(2)

# %%
max_gaps = 0
max_mismatches = 1
max_len_to_alignment_len_abs_diff = 1
max_alignments = 1
processes = 30

# %%
gene = genes[0]
repeat = "1"
    
one_sample_df = best_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df.loc[
    (best_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df["Gene"] == gene)
    & (best_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df["Repeat"] == repeat)
]
one_sample_df

# %%
potential_edges = set(
    (u, v)
    for u, vs in one_sample_df.loc[:, ["Read", "OtherReadswithIndistinguishableUMISubSeqs"]].values.tolist()
    for v in vs
)
ic(len(potential_edges))
assert len(potential_edges) <= one_sample_df["OtherReadswithIndistinguishableUMISubSeqs"].apply(len).sum()
# potential_edges

# %%
read_and_umi_series = one_sample_df.loc[:, ["Read", "SpanningUMISeq"]].set_index("Read").squeeze()
read_and_umi_series

# %%
umi_seqs_overlap_inputs = (
    (
        u, 
        v, 
        read_and_umi_series[u], 
        read_and_umi_series[v], 
        max_len_to_alignment_len_abs_diff,
        max_gaps, 
        max_mismatches, 
        max_alignments
    )
    for u, v in potential_edges
)
umi_seqs_overlap_inputs = [
    input 
    for input in umi_seqs_overlap_inputs
    if umi_processing.umi_seqs_overlap_early_filter(
        input[2], input[3], max_len_to_alignment_len_abs_diff, max_gaps, max_mismatches
    )
]
len(umi_seqs_overlap_inputs)

# %%
vs = one_sample_df.shape[0]

e_clique = vs * (vs - 1) / 2

ic(len(umi_seqs_overlap_inputs), e_clique, np.round(100 * len(umi_seqs_overlap_inputs) / e_clique, 2));

# %%
# umi_seqs_overlap_inputs_batches = more_itertools.divide(processes, umi_seqs_overlap_inputs)

# %%
ctx = mp.get_context("spawn")
with ctx.Pool(processes=processes) as pool:
    umi_seqs_overlap_batched_results = pool.map(
        func=umi_processing.one_batch_umi_seqs_overlap,
        iterable=more_itertools.divide(processes, umi_seqs_overlap_inputs)
    )

# %%
# umi_seqs_overlap_results = list(chain.from_iterable(umi_seqs_overlap_batched_results))

# create a symmetric df with all pairs (U, V) and (V, U)
umi_seqs_overlap_results_df = pd.DataFrame(
    # umi_seqs_overlap_results,
    list(chain.from_iterable(umi_seqs_overlap_batched_results)),
    columns=["U", "V", "Overlap"]
)
umi_seqs_overlap_results_df_2 = umi_seqs_overlap_results_df.copy()
umi_seqs_overlap_results_df_2.insert(0, "U2", umi_seqs_overlap_results_df_2["V"])
umi_seqs_overlap_results_df_2.insert(1, "V2", umi_seqs_overlap_results_df_2["U"])
umi_seqs_overlap_results_df_2 = umi_seqs_overlap_results_df_2.drop(columns=["U", "V"]).rename(columns={"U2": "U", "V2": "V"})
# now the df is symmetric
umi_seqs_overlap_results_df = pd.concat(
    [umi_seqs_overlap_results_df, umi_seqs_overlap_results_df_2],
    ignore_index=True
)

# find for each read the other reads with indistinguishable UMI sequences
reads_and_indistinguishable_reads_df = umi_seqs_overlap_results_df.groupby("U").apply(
    lambda x: x.loc[x["Overlap"], "V"].tolist(),
    include_groups=False
).reset_index(name="OtherReadswithIndistinguishableUMIs").rename(columns={"U": "Read"})
reads_and_indistinguishable_reads_df["NumOfOtherReadswithIndistinguishableUMIs"] = reads_and_indistinguishable_reads_df["OtherReadswithIndistinguishableUMIs"].apply(len)

one_sample_df = one_sample_df.merge(
    reads_and_indistinguishable_reads_df,
    how="left"
)

best_umi_overlap_seq_gene_specific_pcr_amplified_alignments_dfs = [one_sample_df]

best_umi_overlap_seq_gene_specific_pcr_amplified_concat_alignments_df = pd.concat(
    best_umi_overlap_seq_gene_specific_pcr_amplified_alignments_dfs, ignore_index=True
)
best_umi_overlap_seq_gene_specific_pcr_amplified_concat_alignments_df

# %%

# %%

# %%
best_umi_overlap_seq_gene_specific_pcr_amplified_concat_alignments_df.groupby(
    "Sample"
).size()

# %%
best_umi_overlap_seq_gene_specific_pcr_amplified_concat_alignments_df.groupby(
    "Sample"
)["NumOfOtherReadswithIndistinguishableUMIs"].describe().round(2)

# %%
best_umi_overlap_seq_gene_specific_pcr_amplified_concat_alignments_df.loc[
    best_umi_overlap_seq_gene_specific_pcr_amplified_concat_alignments_df[
        "NumOfOtherReadswithIndistinguishableUMIs"
    ].eq(0),
].groupby("Sample").size()

# %%
best_umi_overlap_seq_gene_specific_pcr_amplified_concat_alignments_df.loc[
    best_umi_overlap_seq_gene_specific_pcr_amplified_concat_alignments_df[
        "NumOfOtherReadswithIndistinguishableUMIs"
    ].eq(0),
].groupby("Sample").size().mul(100).div(
    best_umi_overlap_seq_gene_specific_pcr_amplified_concat_alignments_df.groupby(
        "Sample"
    ).size()
).round(2)

# %%
best_umi_overlap_seq_gene_specific_pcr_amplified_concat_alignments_df[
    "NumOfOtherReadswithIndistinguishableUMIs"
].eq(0).sum()

# %%
best_umi_overlap_seq_gene_specific_pcr_amplified_concat_alignments_df.loc[
    best_umi_overlap_seq_gene_specific_pcr_amplified_concat_alignments_df[
        "NumOfOtherReadswithIndistinguishableUMIs"
    ].eq(0),
].shape[0] * 100 / best_umi_overlap_seq_gene_specific_pcr_amplified_concat_alignments_df.shape[0]

# %%




# %%
def compute_reads_with_indistinguishable_umis(
    one_sample_df: pd.DataFrame,
    max_gaps: int = 2,
    max_mismatches: int = 2,
    max_len_to_alignment_len_abs_diff: int = 1, 
    processes: int = 10,
    # batch_size: int | None = None,
    max_alignments: int | None = 2,
):
    """
    For each read, find all other reads that have sufficient overlap with its UMI.
    Returns a DataFrame with new columns:
      - OtherReadswithIndistinguishableUMIIs
      - NumOfOtherReadswithIndistinguishableUMIIs
    """
    reads_and_umis = one_sample_df.loc[:, ["Read", "SpanningUMISeq"]].values.tolist()
    
    umi_seqs_overlap_inputs = [
        (u, v, u_umi_seq, v_umi_seq, max_gaps, max_mismatches, max_len_to_alignment_len_abs_diff, max_alignments)
        for (u, u_umi_seq), (v, v_umi_seq) in itertools.combinations(reads_and_umis, 2)
    ]
    # assert len(umi_seqs_overlap_inputs) == (len(reads_and_umis) * (len(reads_and_umis) - 1)) / 2
    
    # if batch_size is None:
    #     batch_size = len(umi_seqs_overlap_inputs) // processes
    # ic(len(umi_seqs_overlap_inputs), processes, batch_size)
    # umi_seqs_overlap_inputs_batches = itertools.batched(umi_seqs_overlap_inputs, n=batch_size)
    
    # TODO consider using more_itertools' distribute or divide for this purpose to explicitly divide the data according to number of processes
    # https://more-itertools.readthedocs.io/en/stable/api.html#more_itertools.distribute

    umi_seqs_overlap_inputs_batches = more_itertools.divide(processes, umi_seqs_overlap_inputs)

    # with Pool(processes=processes) as pool:
    #     umi_seqs_overlap_batched_results = pool.map(
    #         func=one_batch_umi_seqs_overlap,
    #         iterable=umi_seqs_overlap_inputs_batches
    #     )
        
    ctx = mp.get_context("spawn")
    with ctx.Pool(processes=processes) as pool:
        umi_seqs_overlap_batched_results = pool.map(
            # func=one_batch_umi_seqs_overlap,
            func=umi_processing.one_batch_umi_seqs_overlap,
            iterable=umi_seqs_overlap_inputs_batches
        )
        
    umi_seqs_overlap_results = list(chain.from_iterable(umi_seqs_overlap_batched_results))
    
    # create a symmetric df with all pairs (U, V) and (V, U)
    umi_seqs_overlap_results_df = pd.DataFrame(
        umi_seqs_overlap_results,
        columns=["U", "V", "Overlap"]
    )
    umi_seqs_overlap_results_df_2 = umi_seqs_overlap_results_df.copy()
    umi_seqs_overlap_results_df_2.insert(0, "U2", umi_seqs_overlap_results_df_2["V"])
    umi_seqs_overlap_results_df_2.insert(1, "V2", umi_seqs_overlap_results_df_2["U"])
    umi_seqs_overlap_results_df_2 = umi_seqs_overlap_results_df_2.drop(columns=["U", "V"]).rename(columns={"U2": "U", "V2": "V"})
    # now the df is symmetric
    umi_seqs_overlap_results_df = pd.concat(
        [umi_seqs_overlap_results_df, umi_seqs_overlap_results_df_2],
        ignore_index=True
    )

    # find for each read the other reads with indistinguishable UMI sequences
    reads_and_indistinguishable_reads_df = umi_seqs_overlap_results_df.groupby("U").apply(
        lambda x: x.loc[x["Overlap"], "V"].tolist(),
        include_groups=False
    ).reset_index(name="OtherReadswithIndistinguishableUMIs").rename(columns={"U": "Read"})
    reads_and_indistinguishable_reads_df["NumOfOtherReadswithIndistinguishableUMIs"] = reads_and_indistinguishable_reads_df["OtherReadswithIndistinguishableUMIs"].apply(len)

    one_sample_df = one_sample_df.merge(
        reads_and_indistinguishable_reads_df,
        how="left"
    )
    return one_sample_df


# %%
max_gaps = 2
max_mismatches = 2
max_len_to_alignment_len_abs_diff = 1
max_alignments = 2
# batch_size = 1000
# batch_size = None
processes = 8

# compute_reads_with_indistinguishable_umis_inputs = [
#     (
#         best_gene_specific_pcr_amplified_concat_alignments_df.loc[
#             (best_gene_specific_pcr_amplified_concat_alignments_df["Gene"] == gene)
#             & (best_gene_specific_pcr_amplified_concat_alignments_df["Repeat"] == repeat)
#         ],
#         max_gaps,
#         max_mismatches,
#         max_len_to_alignment_len_abs_diff,
#         max_alignments,
#         batch_size,
#         processes,
#     )
#     for gene, repeat in product(genes, list("123"))
# ]

# with Pool(processes=6) as pool:
#     processed_gene_specific_pcr_amplified_dfs_2 = pool.starmap(
#         # func=process_gene_and_repeat_df,
#         func=compute_reads_with_indistinguishable_umis,
#         iterable=compute_reads_with_indistinguishable_umis_inputs,
#     )


best_umi_overlap_seq_gene_specific_pcr_amplified_alignments_dfs = []

for gene, repeat in product(genes, list("123")):
    
    ic(gene, repeat)
    
    one_sample_df = best_gene_specific_pcr_amplified_concat_alignments_df.loc[
        (best_gene_specific_pcr_amplified_concat_alignments_df["Gene"] == gene)
        & (best_gene_specific_pcr_amplified_concat_alignments_df["Repeat"] == repeat)
    ]
    one_sample_df = compute_reads_with_indistinguishable_umis(
        one_sample_df,
        max_gaps,
        max_mismatches,
        max_len_to_alignment_len_abs_diff,
        processes,
        # batch_size,
        max_alignments,
    )
    best_umi_overlap_seq_gene_specific_pcr_amplified_alignments_dfs.append(
        one_sample_df
    )
    
    break # TODO: remove - just for testing



best_umi_overlap_seq_gene_specific_pcr_amplified_concat_alignments_df = pd.concat(
    best_umi_overlap_seq_gene_specific_pcr_amplified_alignments_dfs, ignore_index=True
)
best_umi_overlap_seq_gene_specific_pcr_amplified_concat_alignments_df

# %%
best_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df.groupby(
    "Sample"
).size()

# %%
best_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df.loc[
    best_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df[
        "NumOfOtherReadswithIndistinguishableUMISubSeqs"
    ].eq(0),
].groupby("Sample").size()

# %%
best_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df.loc[
    best_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df[
        "NumOfOtherReadswithIndistinguishableUMISubSeqs"
    ].eq(0),
].groupby("Sample").size().mul(100).div(
    best_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df.groupby(
        "Sample"
    ).size()
).round(2)

# %%
best_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df[
    "NumOfOtherReadswithIndistinguishableUMISubSeqs"
].eq(0).sum()

# %%
best_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df.loc[
    best_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df[
        "NumOfOtherReadswithIndistinguishableUMISubSeqs"
    ].eq(0),
].shape[0] * 100 / best_umi_sub_seq_gene_specific_pcr_amplified_concat_alignments_df.shape[0]


# %% [markdown]
# ### Analyze graph structure

# %%
def create_indistinguishable_graph(
    gene_and_repeat_df, 
    indistinguishable_reads_col="OtherReadswithIndistinguishableUMIs"
):
    # create g by explicitly deepcopying the needed cols
    reads = gene_and_repeat_df["Read"].tolist()
    other_indistinguishable_reads = (
        gene_and_repeat_df[indistinguishable_reads_col]
        .apply(lambda x: copy.deepcopy(x))
        .tolist()
    )

    G = nx.Graph()

    G.add_nodes_from(reads)
    # ic(G.number_of_nodes())

    for read, neighbours in zip(reads, other_indistinguishable_reads):
        for neighbour in neighbours:
            G.add_edge(read, neighbour)
            
    # ic(G.number_of_edges());
    return G


# %%
# Clear any previous outputs and run the loop cleanly
clear_output()

print("Starting loop execution...")
print("Expected iterations:", len(list(product(genes, list("123")))))
print("-" * 40)

Gs = []
ccs_dfs = []
for i, (gene, repeat) in enumerate(product(genes, list("123"))):
    # i_s.append(i)
    print(f"Iteration {i}: gene={gene}, repeat={repeat}")

    gene_and_repeat_df = best_umi_overlap_seq_gene_specific_pcr_amplified_concat_alignments_df.loc[
        (best_umi_overlap_seq_gene_specific_pcr_amplified_concat_alignments_df["Gene"] == gene)
        & (best_umi_overlap_seq_gene_specific_pcr_amplified_concat_alignments_df["Repeat"] == repeat),
        ["Read", "OtherReadswithIndistinguishableUMIs"]
    ]

    G = create_indistinguishable_graph(gene_and_repeat_df)
    Gs.append(G)
    
    ccs_df = pd.DataFrame(
        {
                "Gene": gene,
                "Repeat": int(repeat),
                "CC":  list(nx.connected_components(G)),
        }
    )
    ccs_df["Size"] = ccs_df["CC"].apply(len)
    ccs_df["Degrees"] = ccs_df["CC"].apply(lambda x: [len(G.adj[read]) for read in x])
    ccs_df["MeanDegree"] = ccs_df.apply(lambda x: sum(x["Degrees"]) / x["Size"], axis=1)
    ccs_dfs.append(ccs_df)
    
    break # TODO: remove - just for testing

print("-" * 40)
print(f"Loop completed. Processed {i} iterations, created {len(Gs)} graphs.")

concat_ccs_df = pd.concat(ccs_dfs, ignore_index=True)

assert concat_ccs_df.loc[
    concat_ccs_df["Degrees"].apply(max).ge(concat_ccs_df["Size"])
].empty, "Max degree in connected components should not exceed the size of the component."

concat_ccs_df["Edges"] = concat_ccs_df["Degrees"].apply(
    lambda x: sum(x) / 2
)
concat_ccs_df["Cliquishness"] = concat_ccs_df.apply(
    # lambda x: x["Edges"] / ((x["Size"] * (x["Size"] - 1)) / 2) if x["Size"] > 1 else np.nan,
    lambda x: x["Edges"] / ((x["Size"] * (x["Size"] - 1)) / 2) if x["Size"] > 1 else 1,
    axis=1
)

concat_ccs_df

# %%
# how many reads are in connected components of size >= 2 that are not cliquish in
concat_ccs_df.loc[
    (concat_ccs_df["Size"] >= 1)
    & (concat_ccs_df["Cliquishness"] < 1),    
].groupby(["Gene", "Repeat"])["Size"].sum()

# %%
# num of nodes per graph
concat_ccs_df.groupby(["Gene", "Repeat"])["Size"].apply(np.sum).astype(int).reset_index()

# %%
# num of edges per graph
concat_ccs_df.groupby(["Gene", "Repeat"])["Edges"].apply(np.sum).astype(int).reset_index()

# %%
# num of connected components per graph
concat_ccs_df[["Gene", "Repeat"]].value_counts().reset_index().sort_values(["Gene", "Repeat"]).reset_index(drop=True)

# %%
# connected components sizes stats
concat_ccs_df.groupby(["Gene", "Repeat"])["Size"].describe().round(2).reset_index()

# %%
# % of connected components with size 1
concat_ccs_df.groupby(["Gene", "Repeat"])["Size"].apply(lambda x: x.eq(1).sum() * 100 / x.size).round(2)

# %%
# % of nodes in connected components with size 1
concat_ccs_df.groupby(["Gene", "Repeat"]).apply(
    lambda x: x.loc[x["Size"].eq(1), "Size"].sum() * 100 / x["Size"].sum(),
    include_groups=False
).round(2)
