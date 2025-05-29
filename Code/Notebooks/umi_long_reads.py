# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.7
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Settings
#

# %%
from collections import Counter, defaultdict
from itertools import chain, product
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.io as pio
import pysam
from Bio import Align
from Bio.Seq import Seq, reverse_complement, complement
from icecream import ic

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

# %%
mapped_bams_dir = Path(
    "/private7/projects/Combinatorics/D.pealeii/Alignment/UMILongReads"
)

# %%
min_read_quality = 0.998

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

# %% [markdown]
# # Unmapped BAMs
#

# %% [markdown]
# ## Concat unmapped bams

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
# ## Finding exact primers in unmapped reads
#

# %%
# len(concat_bams_df["Read"].unique()) == concat_bams_df.shape[0]

# %% jupyter={"source_hidden": true}
# read_seq = read.query_sequence
# read_seq

# %% jupyter={"source_hidden": true}
# read_name = read.query_name
# read_name

# %% jupyter={"source_hidden": true}
# for primer_name, primer_seq in primers_dict.items():
#     lowest_index_forward = read_seq.find(str(primer_seq))
#     lowest_index_reverse = read_seq.find(str(primer_seq.reverse_complement))
#     found_forward = lowest_index_forward >= 0
#     found_reverse = lowest_index_reverse >= 0
#     if found_forward and found_reverse:
#         found = "Both"
#     elif found_forward:
#         found = "Forward"
#     elif found_reverse:
#         found = "Reverse"
#     else:
#         found = False
#     print(f"{primer_name}, {found}, {lowest_index_forward = }, {lowest_index_forward = }")

# %% jupyter={"source_hidden": true}
# primer_seq = Seq("TCGTGAATG")

# %% jupyter={"source_hidden": true}
# str(primer_seq) in read_seq

# %% jupyter={"source_hidden": true}
# read_seq.find(str(primer_seq))

# %% jupyter={"source_hidden": true}
# read_seq.index(str(primer_seq))

# %% jupyter={"source_hidden": true}
# from dataclasses import dataclass

# # @dataclass
# # class PrimersFoundInRead:
# #     Read: str
# #     Primer: str
# #     PrimerFound: str

# @dataclass
# class PrimersFoundInRead:
#     Read: str
#     Primer: str
#     # PrimerFound: str
#     ForwardLocationFound: int
#     ReverseLocationFound: int

# %% jupyter={"source_hidden": true}
# def get_primers_found_per_read_df(primers_dict, read):
#     primers_found_per_read = []

#     read_name = read.query_name
#     read_seq = read.query_sequence

#     for primer_name, primer_seq in primers_dict.items():
#         found_forward = str(primer_seq) in read_seq
#         found_reverse = str(primer_seq.reverse_complement) in read_seq
#         if found_forward and found_reverse:
#             found = "Both"
#         elif found_forward:
#             found = "Forward"
#         elif found_reverse:
#             found = "Reverse"
#         else:
#             found = False
#         primers_found_per_read.append(PrimersFoundInRead(read_name, primer_name, found))

#     df = pd.DataFrame(primers_found_per_read).pivot(columns=["Primer"], index="Read")
#     df.columns = df.columns.droplevel(0)
#     df = df.rename_axis(None, axis=1)
#     df = df.rename_axis(None).reset_index(names="Read")

#     return df

# %% jupyter={"source_hidden": true}
# primers_found_per_read = []

# for primer_name, primer_seq in primers_dict.items():
#     # found_forward = str(primer_seq) in read_seq
#     # found_reverse = str(primer_seq.reverse_complement) in read_seq
#     lowest_index_forward = read_seq.find(str(primer_seq))
#     lowest_index_reverse = read_seq.find(str(primer_seq.reverse_complement))
#     # found_forward = lowest_index_forward >= 0
#     # found_reverse = lowest_index_reverse >= 0
#     # if found_forward and found_reverse:
#     #     found = "Both"
#     # elif found_forward:
#     #     found = "Forward"
#     # elif found_reverse:
#     #     found = "Reverse"
#     # else:
#     #     found = False
#     # primers_found_per_read.append(PrimersFoundInRead(read_name, primer_name, found, lowest_index_forward, lowest_index_reverse))
#     primers_found_per_read.append(PrimersFoundInRead(read_name, primer_name, lowest_index_forward, lowest_index_reverse))
# primers_found_per_read
#     # print(f"{primer_name}, {found}")

# %% jupyter={"source_hidden": true}
# df = pd.DataFrame(primers_found_per_read)
# df

# %% jupyter={"source_hidden": true}
# df.pivot(columns=["Primer"], index="Read")

# %% jupyter={"source_hidden": true}
# df = df.pivot(columns=["Primer"], index="Read")
# df.columns = df.columns.droplevel(0)
# df = df.rename_axis(None, axis=1)
# df = df.rename_axis(None).reset_index(names="Read")
# df

# %%
def get_primers_found_per_read_df(primers_dict, read):
    read_name = read.query_name
    read_seq = read.query_sequence

    # primers_found_per_read = []
    primers_found_per_read = {"Read": read_name}

    for primer_name, primer_seq in primers_dict.items():
        # found_forward = str(primer_seq) in read_seq
        # found_reverse = str(primer_seq.reverse_complement) in read_seq
        # if found_forward and found_reverse:
        #     found = "Both"
        # elif found_forward:
        #     found = "Forward"
        # elif found_reverse:
        #     found = "Reverse"
        # else:
        #     found = False
        # primers_found_per_read.append(PrimersFoundInRead(read_name, primer_name, found))

        lowest_index_forward = read_seq.find(str(primer_seq))
        lowest_index_reverse = read_seq.find(str(primer_seq.reverse_complement()))

        found_forward = lowest_index_forward >= 0
        found_reverse = lowest_index_reverse >= 0
        if found_forward and found_reverse:
            found = "Both"
        elif found_forward:
            found = "Forward"
        elif found_reverse:
            found = "Reverse"
        else:
            found = False

        primers_found_per_read[primer_name] = found
        primers_found_per_read[f"{primer_name}_F"] = lowest_index_forward
        primers_found_per_read[f"{primer_name}_R"] = lowest_index_reverse

    # df = pd.DataFrame(primers_found_per_read).pivot(columns=["Primer"], index="Read")
    # df.columns = df.columns.droplevel(0)
    # df = df.rename_axis(None, axis=1)
    # df = df.rename_axis(None).reset_index(names="Read")

    df = pd.DataFrame([primers_found_per_read])

    return df


# %% jupyter={"source_hidden": true}
# primers_found_per_read = {"Read": read_name}

# for primer_name, primer_seq in primers_dict.items():
#     # found_forward = str(primer_seq) in read_seq
#     # found_reverse = str(primer_seq.reverse_complement) in read_seq
#     lowest_index_forward = read_seq.find(str(primer_seq))
#     lowest_index_reverse = read_seq.find(str(primer_seq.reverse_complement()))

#     found_forward = lowest_index_forward >= 0
#     found_reverse = lowest_index_reverse >= 0
#     if found_forward and found_reverse:
#         found = "Both"
#     elif found_forward:
#         found = "Forward"
#     elif found_reverse:
#         found = "Reverse"
#     else:
#         found = False
#     # primers_found_per_read.append(PrimersFoundInRead(read_name, primer_name, found, lowest_index_forward, lowest_index_reverse))
#     # primers_found_per_read.append(PrimersFoundInRead(read_name, primer_name, lowest_index_forward, lowest_index_reverse))

#     primers_found_per_read[primer_name] = found
#     primers_found_per_read[f"{primer_name}_F"] = lowest_index_forward
#     primers_found_per_read[f"{primer_name}_R"] = lowest_index_reverse

# primers_found_per_read

# %% jupyter={"source_hidden": true}
# get_primers_found_per_read_df(primers_dict, read)

# %% jupyter={"source_hidden": true}
# df = pd.DataFrame([primers_found_per_read])
# df

# %% jupyter={"source_hidden": true}
# df.pivot(columns=["Primer"], index="Read")

# %%
bam_dfs = []
for bam_file in bam_files:
    with pysam.AlignmentFile(
        bam_file,
        "rb",
        threads=10,
        check_sq=False,
        require_index=False,
        # index_filename=str(Path(bam_file.parent, f"{bam_file.name}.pbi")),
    ) as samfile:
        reads = [read for read in samfile]
        reads_names = [read.query_name for read in reads]
        reads_seqs = [read.query_sequence for read in reads]
        # raw_barcods_seqs = [read.get_tag("bc") for read in reads]
        # corrected_barcode_seqs = [read.get_tag("bx") for read in reads]
        # barcode_tags = [read.get_tag("bt") for read in reads]
        reads_lengths = [len(read.query_sequence) for read in reads]
        # unique_reads_names = set(reads_names)
        # ic(bam_file.name.split(".")[0], len(reads_names));
        sample = bam_file.name.split(".")[0]
        gene = sample[3:]
        repeat = sample[2]
        df = pd.DataFrame(
            {
                "Sample": sample,
                "Gene": gene,
                "Repeat": repeat,
                "Read": reads_names,
                "Seq": reads_seqs,
                # "RawBarcodeSeq": raw_barcods_seqs,
                # "CorrectedBarcodeSeq": corrected_barcode_seqs,
                # "BarcodeTag": barcode_tags,
                "ReadLength": reads_lengths,
            }
        )

        primers_found_per_read_dfs = [
            get_primers_found_per_read_df(primers_dict, read) for read in reads
        ]
        primers_found_per_read_df = pd.concat(primers_found_per_read_dfs)

        df = df.merge(primers_found_per_read_df, how="left")

        bam_dfs.append(df)

        # break # todo: remove
concat_bams_df = pd.concat(bam_dfs, ignore_index=True)
concat_bams_df

# %%
concat_bams_df.groupby("Sample").size()

# %%
concat_bams_df.groupby(["Gene", "Sample"]).size().reset_index(name="Reads").groupby(
    "Gene"
)["Reads"].mean()

# %%
concat_bams_df.groupby("Sample").size().mean()

# %%
concat_bams_df.groupby(["Gene", "Sample"]).size().reset_index(name="Reads")

# %%
concat_bams_df.groupby(["Gene", "Sample"]).size().reset_index(name="Reads").groupby(
    "Gene"
)["Reads"].mean()

# %%
# concat_bams_df["ADAR_F"].value_counts()

# %%
concat_bams_df["PCR"].value_counts(normalize=True).mul(100).round(1)

# %%
concat_bams_df["IQEC"].value_counts(normalize=True).mul(100).round(1)

# %%
concat_bams_df.loc[concat_bams_df["Gene"] == "IQEC", "IQEC"].value_counts(
    normalize=True
).mul(100).round(1)

# %%
concat_bams_df.loc[concat_bams_df["Gene"] == "IQEC", ["PCR", "IQEC"]].value_counts(
    normalize=True
).mul(100).round(1)

# %%
concat_bams_df.loc[
    (concat_bams_df["Gene"] == "IQEC")
    & (concat_bams_df["PCR"] != False)
    & (concat_bams_df["IQEC"] != False),
    # ["PCR", "IQEC"]
].groupby("Repeat").size()

# %%
concat_bams_df.loc[concat_bams_df["Gene"] != "IQEC", "IQEC"].value_counts(
    normalize=True
).mul(100)

# %%
concat_bams_df["ADAR1"].value_counts(normalize=True).mul(100).round(1)

# %%
concat_bams_df.loc[concat_bams_df["Gene"] == "ADAR1", "ADAR1"].value_counts(
    normalize=True
).mul(100).round(1)

# %%
concat_bams_df.loc[
    (concat_bams_df["Gene"] == "ADAR1")
    & (concat_bams_df["PCR"] != False)
    & (concat_bams_df["ADAR1"] != False),
    # ["PCR", "IQEC"]
].groupby("Repeat").size()

# %%
concat_bams_df.loc[concat_bams_df["Gene"] != "ADAR1", "ADAR1"].value_counts(
    normalize=True
).mul(100)

# %%
# concat_bams_df.loc[:, primers_dict.keys()].value_counts(normalize=True).mul(100).round(1).reset_index(name="%").sort_values(by=list(primers_dict.keys()))

# %%
# concat_bams_df.groupby("Gene")[list(primers_dict.keys())].value_counts(normalize=True).mul(100).reset_index(name="%")
concat_bams_df.groupby("Gene")[list(primers_dict.keys())].value_counts(
    normalize=True
).mul(100)

# %%
f_r_locations = (
    concat_bams_df.loc[:, ["PCR_F", "PCR_R", "IQEC_F", "IQEC_R", "ADAR1_F", "ADAR1_R"]]
    .stack()
    .reset_index()
)
del f_r_locations["level_0"]
f_r_locations = f_r_locations.rename(columns={"level_1": "Primer", 0: "Location"})
f_r_locations

# %%
f_r_locations.loc[
    (f_r_locations["Primer"].str.endswith("_R")) & (f_r_locations["Location"] != -1)
]

# %%
from random import choices

import numpy as np

# %%
primers_dict

# %% [markdown]
# ### ADAR
#

# %%
# bam_file = bam_files[0]
# ic(bam_file)

# with pysam.AlignmentFile(
#     bam_file,
#     "rb",
#     threads=10,
#     check_sq=False,
#     require_index=False,
#     # index_filename=str(Path(bam_file.parent, f"{bam_file.name}.pbi")),
# ) as samfile:
#     # adar_reads = choices([read for read in samfile], k=3)
#     # rng = np.random.RandomState(18)
#     rng = np.random.default_rng()
#     adar_reads = rng.choice([read for read in samfile], size=3, replace=False)

# adar_read_names = [read.query_name for read in adar_reads]
# adar_read_seqs = [read.query_sequence for read in adar_reads]

# df = pd.DataFrame({"Read": adar_read_names, "Seq": adar_read_seqs})

# adar_test_concat_bams_df = concat_bams_df.loc[concat_bams_df["Read"].isin(adar_read_names)]
# adar_test_concat_bams_df = adar_test_concat_bams_df.merge(df)

# adar_read_names = adar_test_concat_bams_df["Read"].tolist()
# adar_read_seqs = adar_test_concat_bams_df["Seq"].tolist()
# adar_primers_found = adar_test_concat_bams_df["ADAR"].tolist()

# adar_test_concat_bams_df

# %%
# primers_dict["ADAR"]

# %%
# primers_dict["ADAR"].reverse_complement()

# %%
# adar_read_seqs[0]

# %%
# adar_read_seqs[1]

# %%
# adar_read_seqs[2]

# %%
# for adar_read_seq, adar_primer_found in zip(adar_read_seqs, adar_primers_found):
#     if adar_primer_found == "Forward":
#         assert str(primers_dict["ADAR"]) in adar_read_seq
#     elif adar_primer_found == "Reverse":
#         assert str(primers_dict["ADAR"].reverse_complement()) in adar_read_seq
#     elif adar_primer_found == "Both":
#         assert str(primers_dict["ADAR"]) in adar_read_seq and str(primers_dict["ADAR"].reverse_complement()) in adar_read_seq
#     elif adar_primer_found == False:
#         assert str(primers_dict["ADAR"]) not in adar_read_seq and str(primers_dict["ADAR"].reverse_complement()) not in adar_read_seq
#     else:
#         raise Exception(f"{adar_primer_found = } is illegal")

#     # if adar_primer_found:
#     #     assert str(primers_dict["ADAR"]) in adar_read_seq
#     # else:
#     #     assert str(primers_dict["ADAR"]) not in adar_read_seq

# %%

# %%

# %%

# %%
adar_bams_df = concat_bams_df.loc[(concat_bams_df["Gene"] == "ADAR1")]
adar_bams_df

# %%
adar_bams_df

# %%
adar_bams_df[["PCR", "ADAR1", "IQEC"]].value_counts(normalize=True).mul(100).round(4)

# %%
adar_bams_df = adar_bams_df.drop(columns=["IQEC", "IQEC_F", "IQEC_R"], errors="ignore")

# %%
adar_bams_df.loc[
    (adar_bams_df["PCR"] != False) & (adar_bams_df["ADAR1"] != False),
    # ["PCR", "IQEC"]
]

# %%
adar_bams_df.loc[
    (adar_bams_df["PCR"] != False) & (adar_bams_df["ADAR1"] != False), ["PCR", "ADAR1"]
].value_counts(normalize=True).mul(100).round(3)

# %% [markdown]
# #### Forward ADAR
#

# %%
first_primer = "PCR"
last_primer = "ADAR1"

correct_orientation = "F"
wrong_orientation = "R"

# %%
adar_f_bams_df = adar_bams_df.loc[
    (adar_bams_df[first_primer] == "Forward")
    & (adar_bams_df[last_primer] == "Forward"),
].copy()
adar_f_bams_df["Orientation"] = correct_orientation
# adar_f_bams_df = adar_f_bams_df.drop(
#     columns=[
#         f"{first_primer}_{wrong_orientation}",
#         f"{last_primer}_{wrong_orientation}",
#     ],
#     errors="ignore",
# )
adar_f_bams_df["UMILength"] = adar_f_bams_df[
    f"{last_primer}_{correct_orientation}"
].sub(
    adar_f_bams_df[f"{first_primer}_{correct_orientation}"].add(
        len(primers_dict[first_primer])
    )
)
# adar_f_bams_df[f"UMISeq_{correct_orientation}"] = adar_f_bams_df.apply(
#     lambda x: x["Seq"][
#         x[f"{last_primer}_{correct_orientation}"]
#         - x["UMILength"] : x[f"{last_primer}_{correct_orientation}"]
#     ],
#     axis=1,
# )
# rev_comp_umi_seq = adar_f_bams_df[f"UMISeq_{correct_orientation}"].apply(
#     lambda x: str(Seq(x).reverse_complement())
# )
# if correct_orientation == "F":
#     adar_f_bams_df[f"UMISeq_{wrong_orientation}"] = rev_comp_umi_seq
# else:
#     adar_f_bams_df.insert(
#         adar_f_bams_df.columns.get_loc(f"UMISeq_{correct_orientation}"),
#         f"UMISeq_{wrong_orientation}",
#         rev_comp_umi_seq,
#     )
adar_f_bams_df["UMISeq"] = adar_f_bams_df.apply(
    lambda x: x["Seq"][
        x[f"{last_primer}_{correct_orientation}"]
        - x["UMILength"] : x[f"{last_primer}_{correct_orientation}"]
    ],
    axis=1,
)

# umi_seqs = adar_f_bams_df["UMISeq"].to_list()
# parent_umi_seqs = defaultdict(list)
# # for s1 in umi_seqs:
# #     for s2 in umi_seqs:
# #         if s1 != s2 and s1 in s2:
# #             parent_umi_seqs[s1].append(s2)
# # Counter(len(v) for v in parent_umi_seqs.values())
# parent_umi_seqs = set()
# for s1 in umi_seqs:
#     for s2 in umi_seqs:
#         if s1 != s2 and s1 in s2:
#             parent_umi_seqs.add(s1)
#             break
# adar_f_bams_df["SubstringOfOtherUMI"] = [(x in parent_umi_seqs) for x in umi_seqs]

adar_f_bams_df

# %%
adar_f_bams_df.groupby("Repeat")["UMISeq"].nunique()

# %%
adar_f_bams_df["UMILength"].describe()

# %%
100 * adar_f_bams_df.loc[adar_f_bams_df["UMILength"] == 12].shape[
    0
] / adar_f_bams_df.shape[0]

# %%
adar_f_bams_df.loc[:, ["PCR_F", "ADAR1_F"]].value_counts(normalize=True).mul(100).round(
    3
)

# %%
adar_f_bams_df.loc[:, ["PCR_F", "ADAR1_F"]].describe()

# %% [markdown]
# #### Reverse ADAR
#

# %%
first_primer = "ADAR1"
last_primer = "PCR"

correct_orientation = "R"
wrong_orientation = "F"

# %%
adar_r_bams_df = adar_bams_df.loc[
    (adar_bams_df[first_primer] == "Reverse")
    & (adar_bams_df[last_primer] == "Reverse"),
].copy()
adar_r_bams_df["Orientation"] = correct_orientation
# adar_r_bams_df = adar_r_bams_df.drop(
#     columns=[
#         f"{first_primer}_{wrong_orientation}",
#         f"{last_primer}_{wrong_orientation}",
#     ],
#     errors="ignore",
# )
adar_r_bams_df["UMILength"] = adar_r_bams_df[
    f"{last_primer}_{correct_orientation}"
].sub(
    adar_r_bams_df[f"{first_primer}_{correct_orientation}"].add(
        len(primers_dict[first_primer])
    )
)
# adar_r_bams_df[f"UMISeq_{correct_orientation}"] = adar_r_bams_df.apply(
#     lambda x: x["Seq"][
#         x[f"{last_primer}_{correct_orientation}"]
#         - x["UMILength"] : x[f"{last_primer}_{correct_orientation}"]
#     ],
#     axis=1,
# )
# rev_comp_umi_seq = adar_r_bams_df[f"UMISeq_{correct_orientation}"].apply(
#     lambda x: str(Seq(x).reverse_complement())
# )
# if correct_orientation == "F":
#     adar_r_bams_df[f"UMISeq_{wrong_orientation}"] = rev_comp_umi_seq
# else:
#     adar_r_bams_df.insert(
#         adar_r_bams_df.columns.get_loc(f"UMISeq_{correct_orientation}"),
#         f"UMISeq_{wrong_orientation}",
#         rev_comp_umi_seq,
#     )
adar_r_bams_df["UMISeq"] = adar_r_bams_df.apply(
    lambda x: x["Seq"][
        x[f"{last_primer}_{correct_orientation}"]
        - x["UMILength"] : x[f"{last_primer}_{correct_orientation}"]
    ],
    axis=1,
)
adar_r_bams_df["UMISeq"] = adar_r_bams_df["UMISeq"].apply(
    lambda x: str(Seq(x).reverse_complement())
)
adar_r_bams_df

# %%
adar_r_bams_df.groupby("Repeat")["UMISeq"].nunique()

# %%
adar_r_bams_df["UMILength"].describe()

# %%
100 * adar_r_bams_df.loc[adar_r_bams_df["UMILength"] == 12].shape[
    0
] / adar_r_bams_df.shape[0]


# %%

# %%

# %% [markdown]
# #### ADAR recombined
#

# %%
def find_umis_which_are_substrings_of_others(umi_seqs):
    parent_umi_seqs = set()
    for s1 in umi_seqs:
        for s2 in umi_seqs:
            if s1 != s2 and s1 in s2:
                parent_umi_seqs.add(s1)
                break
    # return {x: x in parent_umi_seqs for x in umi_seqs}
    # return pd.DataFrame([(x, x in parent_umi_seqs) for x in umi_seqs], columns=["UMISeq", "UMISeqIsSubstringOfOther"])
    return [x in parent_umi_seqs for x in umi_seqs]


# %%

# %%

# %%
f_r_adar_bam_dfs = pd.concat([adar_f_bams_df, adar_r_bams_df]).sort_values(
    ["Repeat"], ignore_index=True
)

per_repeat_f_r_adar_bam_dfs = [
    f_r_adar_bam_dfs.loc[f_r_adar_bam_dfs["Repeat"] == repeat].copy()
    for repeat in list("123")
]
for df in per_repeat_f_r_adar_bam_dfs:
    df["UMISeqIsSubstringOfOther"] = find_umis_which_are_substrings_of_others(
        df["UMISeq"]
    )
per_repeat_f_r_adar_bam_dfs[0]

f_r_adar_bam_dfs = pd.concat(per_repeat_f_r_adar_bam_dfs, ignore_index=True)

f_r_adar_bam_dfs

# %%
f_r_adar_bam_dfs["UMILength"].value_counts()

# %%
f_r_adar_bam_dfs["UMILength"].describe()

# %%
fig = px.histogram(
    f_r_adar_bam_dfs, x="UMILength", color="Sample", log_x=True, log_y=True
)
# Reduce opacity to see both histograms
fig.update_traces(opacity=0.75)
# Overlay both histograms
fig.update_layout(barmode="overlay", template="plotly_white", width=600, height=400)

fig.show()

# %%
f_r_adar_bam_dfs.groupby("Repeat")["UMISeqIsSubstringOfOther"].value_counts(
    normalize=True, dropna=False
).mul(100).round(2)

# %%
f_r_adar_bam_dfs.groupby("Sample").size()

# %%
f_r_adar_bam_dfs.groupby("Sample").size().mean()

# %%
# unique umis per ADAR sample, when considering together F and R reads
f_r_adar_bam_dfs.loc[~f_r_adar_bam_dfs["UMISeqIsSubstringOfOther"]].groupby("Sample")[
    "UMISeq"
].nunique().reset_index(name="UniqueUMIs")

# %%
# unique umis per ADAR sample, when considering together F and R reads
f_r_adar_bam_dfs.loc[~f_r_adar_bam_dfs["UMISeqIsSubstringOfOther"]].groupby("Sample")[
    "UMISeq"
].nunique().mean()

# %%
# unique umis per ADAR sample, when separately considering F and R reads in each sample
f_r_adar_bam_dfs.loc[~f_r_adar_bam_dfs["UMISeqIsSubstringOfOther"]].groupby(
    ["Sample", "Orientation"]
)["UMISeq"].nunique().reset_index(name="UniqueUMIs").groupby("Sample")[
    "UniqueUMIs"
].sum().reset_index()

# %%

# %% [markdown]
# ### IQEC
#

# %%
# bam_file = bam_files[3]
# ic(bam_file)

# with pysam.AlignmentFile(
#     bam_file,
#     "rb",
#     threads=10,
#     check_sq=False,
#     require_index=False,
#     # index_filename=str(Path(bam_file.parent, f"{bam_file.name}.pbi")),
# ) as samfile:
#     rng = np.random.default_rng()
#     iqec_reads = rng.choice([read for read in samfile], size=3, replace=False)

# iqec_read_names = [read.query_name for read in iqec_reads]
# iqec_read_seqs = [read.query_sequence for read in iqec_reads]

# df = pd.DataFrame({"Read": iqec_read_names, "Seq": iqec_read_seqs})

# iqec_test_concat_bams_df = concat_bams_df.loc[concat_bams_df["Read"].isin(iqec_read_names)]
# iqec_test_concat_bams_df = iqec_test_concat_bams_df.merge(df)

# iqec_read_names = iqec_test_concat_bams_df["Read"].tolist()
# iqec_read_seqs = iqec_test_concat_bams_df["Seq"].tolist()
# iqec_primers_found = iqec_test_concat_bams_df["IQEC"].tolist()

# iqec_test_concat_bams_df

# %%
# primers_dict["IQEC"]

# %%
# primers_dict["IQEC"].reverse_complement()

# %%
# iqec_read_seqs[0]

# %%
# iqec_read_seqs[1]

# %%
# iqec_read_seqs[2]

# %%
# for iqec_read_seq, iqec_primer_found in zip(iqec_read_seqs, iqec_primers_found):
#     if iqec_primer_found == "Forward":
#         assert str(primers_dict["IQEC"]) in iqec_read_seq
#     elif iqec_primer_found == "Reverse":
#         assert str(primers_dict["IQEC"].reverse_complement()) in iqec_read_seq
#     elif iqec_primer_found == "Both":
#         assert str(primers_dict["IQEC"]) in iqec_read_seq and str(primers_dict["IQEC"].reverse_complement()) in iqec_read_seq
#     elif iqec_primer_found == False:
#         assert str(primers_dict["IQEC"]) not in iqec_read_seq and str(primers_dict["IQEC"].reverse_complement()) not in iqec_read_seq
#     else:
#         raise Exception(f"{iqec_primer_found = } is illegal")

#     # if adar_primer_found:
#     #     assert str(primers_dict["ADAR"]) in adar_read_seq
#     # else:
#     #     assert str(primers_dict["ADAR"]) not in adar_read_seq

# %%
iqec_bams_df = concat_bams_df.loc[(concat_bams_df["Gene"] == "IQEC")]
iqec_bams_df

# %%
iqec_bams_df[["PCR", "IQEC", "ADAR1"]].value_counts(normalize=True).mul(100).round(4)

# %%
iqec_bams_df = iqec_bams_df.drop(
    columns=["ADAR1", "ADAR1_F", "ADAR1_R"], errors="ignore"
)

# %%
iqec_bams_df.loc[
    (iqec_bams_df["PCR"] != False) & (iqec_bams_df["IQEC"] != False),
    # ["PCR", "IQEC"]
]

# %%
primers_dict

# %%
[len(p) for p in primers_dict.values()]

# %%
iqec_bams_df.loc[
    (iqec_bams_df["PCR"] == False) & (iqec_bams_df["IQEC"] == False),
    # ["PCR", "IQEC"]
].iloc[0]

# %%
iqec_bams_df.loc[
    (iqec_bams_df["PCR"] == False) & (iqec_bams_df["IQEC"] == False),
    # ["PCR", "IQEC"]
].iloc[0]["Seq"]

# %%
iqec_bams_df.loc[
    (iqec_bams_df["PCR"] != False) & (iqec_bams_df["IQEC"] != False), ["PCR", "IQEC"]
].value_counts(normalize=True).mul(100).round(3)

# %% [markdown]
# #### Forward IQEC
#

# %%
first_primer = "PCR"
last_primer = "IQEC"

correct_orientation = "F"
wrong_orientation = "R"

# %%
iqec_f_bams_df = iqec_bams_df.loc[
    (iqec_bams_df[first_primer] == "Forward")
    & (iqec_bams_df[last_primer] == "Forward"),
].copy()
iqec_f_bams_df["Orientation"] = correct_orientation

iqec_f_bams_df["UMILength"] = iqec_f_bams_df[
    f"{last_primer}_{correct_orientation}"
].sub(
    iqec_f_bams_df[f"{first_primer}_{correct_orientation}"].add(
        len(primers_dict[first_primer])
    )
)

iqec_f_bams_df["UMISeq"] = iqec_f_bams_df.apply(
    lambda x: x["Seq"][
        x[f"{last_primer}_{correct_orientation}"]
        - x["UMILength"] : x[f"{last_primer}_{correct_orientation}"]
    ],
    axis=1,
)

iqec_f_bams_df

# %%
iqec_f_bams_df["UMILength"].describe()

# %%
100 * iqec_f_bams_df.loc[iqec_f_bams_df["UMILength"] == 12].shape[
    0
] / iqec_f_bams_df.shape[0]

# %% [markdown]
# #### Reverse IQEC
#

# %%
first_primer = "IQEC"
last_primer = "PCR"

correct_orientation = "R"
wrong_orientation = "F"

# %%
iqec_r_bams_df = iqec_bams_df.loc[
    (iqec_bams_df[first_primer] == "Reverse")
    & (iqec_bams_df[last_primer] == "Reverse"),
].copy()
iqec_r_bams_df["Orientation"] = correct_orientation

iqec_r_bams_df["UMILength"] = iqec_r_bams_df[
    f"{last_primer}_{correct_orientation}"
].sub(
    iqec_r_bams_df[f"{first_primer}_{correct_orientation}"].add(
        len(primers_dict[first_primer])
    )
)

iqec_r_bams_df["UMISeq"] = iqec_r_bams_df.apply(
    lambda x: x["Seq"][
        x[f"{last_primer}_{correct_orientation}"]
        - x["UMILength"] : x[f"{last_primer}_{correct_orientation}"]
    ],
    axis=1,
)
iqec_r_bams_df["UMISeq"] = iqec_r_bams_df["UMISeq"].apply(
    lambda x: str(Seq(x).reverse_complement())
)
iqec_r_bams_df

# %%
iqec_r_bams_df["UMILength"].describe()

# %%
100 * iqec_r_bams_df.loc[iqec_r_bams_df["UMILength"] == 12].shape[
    0
] / iqec_r_bams_df.shape[0]

# %% [markdown]
# #### IQEC recombined
#

# %%
f_r_iqec_bam_dfs = pd.concat([iqec_f_bams_df, iqec_r_bams_df]).sort_values(
    ["Repeat"], ignore_index=True
)

per_repeat_f_r_iqec_bam_dfs = [
    f_r_iqec_bam_dfs.loc[f_r_iqec_bam_dfs["Repeat"] == repeat].copy()
    for repeat in list("123")
]
for df in per_repeat_f_r_iqec_bam_dfs:
    df["UMISeqIsSubstringOfOther"] = find_umis_which_are_substrings_of_others(
        df["UMISeq"]
    )
# per_repeat_f_r_adar_bam_dfs[0]

f_r_iqec_bam_dfs = pd.concat(per_repeat_f_r_iqec_bam_dfs, ignore_index=True)

f_r_iqec_bam_dfs

# %%
fig = px.histogram(
    f_r_iqec_bam_dfs.loc[f_r_iqec_bam_dfs["UMILength"] >= 0],
    x="UMILength",
    color="Sample",
    # facet_row="Sample",
    log_x=True,
    log_y=True,
)
# Reduce opacity to see both histograms
fig.update_traces(opacity=0.5)
# Overlay both histograms
fig.update_layout(
    barmode="overlay",
    template="plotly_white",
    width=600,
    height=400,
    # showlegend=False
)

fig.show()

# %%
len(
    "GCGCAGGCGTTGGCCGACAAAGCGGCTGAGAAAGGAGCTGACAAGTCTGGTACTGATTCATTGGCACCAAATCTACAGATCACCTCAGAAAGTTTCGCCGCTCTCAACAAAAATCCAGTTAGCGCACTGATGGAATATGCTCAACAGCGACACTTACCCGTTGAATTTAAGCTTTTGTCACACAGAGGACCTTCTCATCGACCGTTGTTCAAATTTGCCGTGATTCTTGGTAAACGCCAGTTCCCCAGTATGGAGTGCAACAGTAAGAAGGATGGTAAGAAAGAGGCAGCCGATCTGACATTGCGCATTCTCATTGCTGAAGGACAGTATCAACTGGAGAACACCGTCTCAGCATTGAAAACAATTCCACCTGCTGAAATGACACATTTCGACAAAATGGCTGCCTTGAGTCACCAGGCATTTAACAACATTGCCTTGCAAATCCCTGAGAACCTTGCTGGGAGAAAGGTCATCGCTGCTTTGGTGATGAAGCGATCACCAACGGATACGGGAATTGTTATCAGTGTTGGAACTGGTAACCGCTGTTTAACCGGTGATCATTTGAGTTTGGAAGGCAACAGTGTCAATGACTCTCATGCTGAAATAATCACACGCCGAGGTTTTCTGAGATATCTGTACAAACATTTACTGGAGTATGATCCCGAAAAACCCCATGACCTATTTGAGAAAGGTGAACGTAGTCTTTGCCGGATAAAAACCAACATTACATTCCATCTGTATATATCAACTGCTCCTTGTGGTGATGGAGCACTTTTTTCACCCAGGGATACCGACTCCAGTAATGTGAAAGTGGATGAGGAAAATAAGCACGTCCATAATCCGACTTTTTCAAGCAGTGTTCAGGGATTGCTGAGAACCAAAGTGGAAGGAGGTGAAGGGACCATTCCAATAGATGCTGATTTCACTGAACAAACATGGGATGGAATTCAACGAGGTGAAAGATTGCGCACAATGTCATGTTCAGATAAAATATGTCGATGGAACGTTGTTGGTCTGCAAGGAGCTTTGCTTAGTCACTTTGTGGAACCAATCTACCTGGAATCTCTGACATTAGGTTATCTTTATGATCATGGCCACTTAGCACGAGCTGTTTGCTGCCGTATTGAACGGGGAGAGGCCTCTGTCAACCAACTACTACCTGAGGGCTACCGATTGAACCATCCTTGGCTTGGCAGAGTTACTGCTTGTGATCCACCTAGAGAAACCCAAAAGACGAAATCGTTAAGTATCAACTGGTGCTATGATGATGAAAAGTCTGAAGTTCTCGATGGTACAGCAGGCATCTGTTACACAGCGATTGAGAAAAATCTCTTCTCTCGCTTAACAAAGCACAGCTTATATGAAGAATTCAAAAAAGTGTGTCAGAAATTTGAACGCGAGGACTTGATGAATGTCACTTCTTACAACAAAGCCAAGATGATGGCCATTCCTTTCCAGACTGCCAAAAACGTAATGTTGAAAAAACTCAAGGAAAACAACTGCGGAACTTGGGTGTCAAAACCTATTGAGGAGGAGATGTTTACGTAGATGTCTGCTGCTACAACACCTATATAGTACATAAAGACAAAAACTACTATTATAAATGTTGC"
)

# %%
f_r_iqec_bam_dfs.groupby("Repeat")["UMISeqIsSubstringOfOther"].value_counts(
    normalize=True, dropna=False
).mul(100).round(2)

# %%
f_r_iqec_bam_dfs.groupby("Sample").size()

# %%
f_r_iqec_bam_dfs.groupby("Sample").size().mean()

# %%
# unique umis per IQEC sample, when considering together F and R reads
f_r_iqec_bam_dfs.loc[~f_r_iqec_bam_dfs["UMISeqIsSubstringOfOther"]].groupby("Sample")[
    "UMISeq"
].nunique().reset_index(name="UniqueUMIs")

# %%
# unique umis per IQEC sample, when separately considering F and R reads in each sample
f_r_iqec_bam_dfs.loc[~f_r_iqec_bam_dfs["UMISeqIsSubstringOfOther"]].groupby(
    ["Sample", "Orientation"]
)["UMISeq"].nunique().reset_index(name="UniqueUMIs").groupby("Sample")[
    "UniqueUMIs"
].sum().reset_index()

# %%

# %%

# %% [markdown]
# # Mapped BAMs
#

# %%
mapped_bam_files = list(mapped_bams_dir.glob("*.bam"))
mapped_bam_files

# %%
# # !samtools view -c --threads 10 /private7/projects/Combinatorics/D.pealeii/Alignment/UMILongReads/LP2IQEC.r64296e203404D01.aligned.sorted.bam

# %%
mapped_bam_file = mapped_bam_files[0]
mapped_bam_file

# %% jupyter={"source_hidden": true}
with pysam.AlignmentFile(
    mapped_bam_file,
    "rb",
    threads=10,
) as samfile:
    mapped_reads = [read for read in samfile]
mapped_reads[0]
print(mapped_reads[0])

# %% jupyter={"source_hidden": true}
# seq = read.query_sequence
# len(seq)

# %% jupyter={"source_hidden": true}
# soft_clipped_seq = read.query_sequence
# len(soft_clipped_seq)

# %%
# primers_dict

# %% jupyter={"source_hidden": true}
read = mapped_reads[0]
print(read)

# %%
read.get_cigar_stats()

# %%
read.get_aligned_pairs(with_cigar=True)[0][2] == 4

# %%
read.get_aligned_pairs(with_cigar=True)


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
        # reads_mapped_to_positive_strand = [read.is_forward for read in reads]
        mapped_strands = ["+" if read.is_forward else "-" for read in reads]
        # reads_soft_clipped_seqs = [read.query_sequence for read in reads]
        # reads_soft_clipped_lengths = [len(seq) for seq in reads_soft_clipped_seqs]
        # aligned_seqs = [read.query_alignment_sequence for read in reads]
        # aligned_seqs_lengths = [len(seq) for seq in aligned_seqs]
        mapped_chroms = [read.reference_name for read in reads]
        # mapped_genes = [gene_by_chrom_dict[read.reference_name] for read in reads]
        # aligned_pairs = [
        #     read.get_aligned_pairs(matches_only=False, with_cigar=True)
        #     for read in reads
        # ]

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
                "MappedStrand": mapped_strands,
                # "RTGStrand": mapped_strands,
                # "SoftClippedSeq": reads_soft_clipped_seqs,
                # "SoftClippedSeqLength": reads_soft_clipped_lengths,
                # "AlignedSeq": aligned_seqs,
                # "AlignedSeqLength": aligned_seqs_lengths,
                # "ExpectedChrom": expected_chrom,
                "MappedChrom": mapped_chroms,
                # "MappedGene": mapped_genes,
            }
        )

        alignment_boundries_df = pd.DataFrame(
            all_reads_and_targets_aligned_starts_and_ends,
            columns=["ReadStart", "ReadEnd", "GeneStart", "GeneEnd"],
            # columns=["RTGReadStart", "RTGReadEnd", "RTGGeneStart", "RTGGeneEnd"],
        )

        df = pd.concat([df, alignment_boundries_df], axis=1)

        # def expected_gene_mapping_status(
        #     expected_chrom, mapped_chrom, gene_by_chrom_dict
        # ):
        #     if expected_chrom == mapped_chrom:
        #         return "Correct"
        #     elif mapped_chrom in gene_by_chrom_dict:
        #         return "Incorrect"
        #     else:
        #         return "Other"

        # df["MappedToExpectedChrom"] = df["ExpectedChrom"].eq(df["MappedChrom"])
        # df["MappedToExpectedGene"] = df["Gene"].eq(df["MappedGene"])
        # df["MappedGene"] = df.apply(
        #     lambda x: gene_by_chrom_dict.get(x["MappedChrom"], "Other"), axis=1
        # )
        # df["ExpectedGeneMappingStatus"] = df.apply(
        #     lambda x: expected_gene_mapping_status(
        #         x["ExpectedChrom"], x["MappedChrom"], gene_by_chrom_dict
        #     ),
        #     axis=1,
        # )

        # # this assertion is needed in order to be sure we can rev-comp both
        # # the "SoftClippedSeq" and "AlignedSeq" cols of reads mapped to
        # # the negative strand
        # assert (
        #     df.apply(
        #         lambda x: x["AlignedSeq"] in x["SoftClippedSeq"], axis=1
        #     ).value_counts()[True]
        #     == df.shape[0]
        # )

        # for col in ["SoftClippedSeq", "AlignedSeq"]:
        #     df.loc[df["MappedStrand"].eq("-"), col] = df.loc[
        #         df["MappedStrand"].eq("-"), col
        #     ].apply(reverse_complement)

        # df.loc[df["MappedStrand"].eq("-"), col] = df.loc[
        #         df["MappedStrand"].eq("-"), col
        #     ].apply(reverse_complement)

        mapped_bam_dfs.append(df)

        # break

concat_mapped_bams_df = pd.concat(mapped_bam_dfs, ignore_index=True)
concat_mapped_bams_df

# %%
# concat_alignments_df

# %%

# %%
mapped_bam_file = mapped_bam_files[0]
mapped_bam_file

# %%
# example of read mapped to the - strand

read_name = "m64296e_241222_071206/74516384/ccs"

with pysam.AlignmentFile(
    mapped_bam_file,
    "rb",
    threads=10,
) as samfile:
    read = [read for read in samfile if read.query_name == read_name][0]
print(read)

# %%
aligned_read_target_pairs = read.get_aligned_pairs(matches_only=True)
aligned_read_start, aligned_target_start = aligned_read_target_pairs[0]
aligned_read_end, aligned_target_end = aligned_read_target_pairs[-1]
aligned_read_end += 1
aligned_target_end += 1

(
    aligned_read_start,
    aligned_read_end,
    aligned_target_start,
    aligned_target_end,
)

# %%
# example of read mapped to the + strand

read_name = "m64296e_241222_071206/139266104/ccs"

with pysam.AlignmentFile(
    mapped_bam_file,
    "rb",
    threads=10,
) as samfile:
    read = [read for read in samfile if read.query_name == read_name][0]
print(read)

# %%
aligned_read_target_pairs = read.get_aligned_pairs(matches_only=True)
aligned_read_start, aligned_target_start = aligned_read_target_pairs[0]
aligned_read_end, aligned_target_end = aligned_read_target_pairs[-1]
aligned_read_end += 1
aligned_target_end += 1

(
    aligned_read_start,
    aligned_read_end,
    aligned_target_start,
    aligned_target_end,
)

# %%

# %%

# %% [markdown]
# # Concat BAMs
#

# %% jupyter={"source_hidden": true}
# concat_bams_df = (
#     concat_unmapped_bams_df.drop(
#         columns=["RawBarcodeSeq", "CorrectedBarcodeSeq", "BarcodeTag"]
#     )
#     # .rename(columns={"Seq": "UnalignedSeq", "ReadLength": "UnalignedSeqLength"})
#     # .rename(columns={"ReadLength": "SeqLength"})
#     .rename(columns={"Seq": "RawSeq", "ReadLength": "RawSeqLength"}).merge(
#         concat_mapped_bams_df
#         # .rename(columns={"Seq": "AlignedSeq", "ReadLength": "AlignedSeqLength"})
#         ,
#         how="outer",
#         indicator="indicator",
#     )
# )

# # make sure all aligned reads also appear in the unaligned reads' bams
# assert concat_bams_df["indicator"].value_counts()["right_only"] == 0

# concat_bams_df.insert(
#     concat_bams_df.columns.get_loc("Read") + 1,
#     "Mapped",
#     concat_bams_df["indicator"].apply(lambda x: True if x == "both" else False),
# )

# # after using this col for validating a correct merge,
# # and for determining wether a read is edited,
# # we don't need it anymore
# del concat_bams_df["indicator"]

# concat_bams_df["%AlignedSeqLength/RawSeqLength"] = (
#     concat_bams_df["AlignedSeqLength"].mul(100).div(concat_bams_df["RawSeqLength"])
# )

# # make sure that for mapped reads whose soft-clipped aligned seq is of the
# # same length as the original ("Raw") seq, these sequences are also identical
# # (this is due to the rev-comp of the soft-clipped seqs before)
# assert (
#     concat_bams_df.loc[
#         (concat_bams_df["Mapped"])
#         & (concat_bams_df["RawSeqLength"] == concat_bams_df["SoftClippedSeqLength"]),
#         "RawSeq",
#     ]
#     .eq(
#         concat_bams_df.loc[
#             (concat_bams_df["Mapped"])
#             & (
#                 concat_bams_df["RawSeqLength"] == concat_bams_df["SoftClippedSeqLength"]
#             ),
#             "SoftClippedSeq",
#         ]
#     )
#     .all()
# )

# concat_bams_df.insert(
#     concat_bams_df.columns.get_loc("ReadQuality") + 1,
#     "HighQualityRead",
#     concat_bams_df["ReadQuality"].ge(min_read_quality),
# )

# concat_bams_df

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
concat_bams_df.groupby(["Gene", "Repeat"])["MappedGene"].value_counts()

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
# # Find barcodes in reads
#

# %%
# concat_bams_df.head()

# %%
def find_first_n_pairwise_alignments_for_barcode_and_read(
    sample,
    gene,
    repeat,
    barcode,  # query_name
    read,  # target_name
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

            prct_barcode_identity = 100 * identities / len(barcode_seq)
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
# alignments_df_1_2_b = align_reads_to_barcode_by_parallel_batches(
#     processes,
#     pd.concat([split_df_1, split_df_2], ignore_index=True),
#     gene,
#     barcode,
#     barcode_seq,
#     first_n_alignments_per_strand=first_n_alignments_per_strand,
#     debug_mode=False,
#     batches_per_process=10,
# )
# alignments_df_1_2_b

# %%

# %%

# %%
# processes = 10
# first_n_alignments_per_strand = 100
# debug_mode = False

# aligner = Align.PairwiseAligner(
#     # match_score=1.0, mismatch_score=-2.0, gap_score = -2.5,
#     mode="local",  # otherwise we'd get scattered matches of the primer across the read
#     scoring="blastn",
#     # scoring="megablast"
# )

# %%
# gene = genes[0]

# barcode = gene
# # query_name = barcode
# barcode_seq = primers_dict[barcode]

# %%
# one_gene_bams_df = concat_bams_df.loc[(concat_bams_df["Gene"].eq(gene))]
# one_gene_bams_df

# %%
# batches_per_process = 10
# # first_n_alignments_per_strand = 10
# debug_mode = False
# n = processes * batches_per_process

# # Compute split indices
# split_sizes = [
#     len(one_gene_bams_df) // n + (i < len(one_gene_bams_df) % n) for i in range(n)
# ]
# split_indices = [0] + list(pd.Series(split_sizes).cumsum())

# # Split DataFrame
# split_dfs = [
#     one_gene_bams_df.iloc[split_indices[i] : split_indices[i + 1]] for i in range(n)
# ]

# assert pd.concat(split_dfs).equals(one_gene_bams_df)
# assert sum(len(df) for df in split_dfs) == len(one_gene_bams_df)

# # split_df = split_dfs[0].iloc[:3]
# # split_df = split_dfs[0]
# # split_df

# %%
# split_df_1 = split_dfs[0]
# split_df_2 = split_dfs[1]

# %%
# split_df_1

# %%
# split_df_2

# %%
# 2257 + 2257

# %%
# pd.concat([split_df_1, split_df_2], ignore_index=True)

# %%
# alignments_df_1 = align_reads_to_barcode_one_batch(
#     one_gene_bams_df=split_df_1,
#     gene=gene,
#     barcode=barcode,  # primer name
#     barcode_seq=barcode_seq,  # primer seq
#     first_n_alignments_per_strand=first_n_alignments_per_strand,
#     debug_mode=debug_mode,
# )

# alignments_df_1

# %%
# alignments_df_2 = align_reads_to_barcode_one_batch(
#     one_gene_bams_df=split_df_2,
#     gene=gene,
#     barcode=barcode,  # primer name
#     barcode_seq=barcode_seq,  # primer seq
#     first_n_alignments_per_strand=first_n_alignments_per_strand,
#     debug_mode=debug_mode,
# )

# alignments_df_2

# %%
# alignments_df_1_2_a = pd.concat([alignments_df_1, alignments_df_2], ignore_index=True)
# alignments_df_1_2_a

# %%
# alignments_df_1_2_b = align_reads_to_barcode_by_parallel_batches(
#     processes,
#     pd.concat([split_df_1, split_df_2], ignore_index=True),
#     gene,
#     barcode,
#     barcode_seq,
#     first_n_alignments_per_strand=first_n_alignments_per_strand,
#     debug_mode=False,
#     batches_per_process=10,
# )
# alignments_df_1_2_b

# %%
# .shape[0] == alignments_df_1_2.shape[0]

# %%
# alignments_df["BTRNumOfAllAlignmentsPerStrand"].describe()

# %%
# alignments_df.loc[
#     (alignments_df["Gene"].eq(alignments_df["MappedGene"]))
#     & (alignments_df["BTRNumOfAllAlignmentsPerStrand"].gt(1))
# ]["BTRNumOfAllAlignmentsPerStrand"].describe()

# %%
# read = "m64296e_241222_071206/100007945/ccs"

# test_df = concat_bams_df.loc[concat_bams_df["Read"].eq(read)]
# test_df

# %%
# gene = test_df.iloc[0]["Gene"]
# barcode = gene
# barcode_seq = primers_dict[barcode]

# alignments_df = test_df.apply(
#     lambda x: find_first_n_pairwise_alignments_for_barcode_and_read(
#         x["Sample"],
#         x["Gene"],
#         x["Repeat"],
#         aligner,
#         barcode,
#         x["Read"],  # target_name
#         x["MappedGene"],  # mapped_gene
#         x[
#             "MappedStrand"
#         ],  # genomic_target_strand - the genomic strand to which the read was mapped to
#         x["ReadStart"],
#         x["ReadEnd"],
#         x["GeneStart"],
#         x["GeneEnd"],
#         barcode_seq,
#         x["Seq"],  # read_seq
#         first_n_alignments_per_strand=first_n_alignments_per_strand,
#         debug_mode=debug_mode,
#     ),
#     axis=1,
# )
# alignments_df

# %%
# sample = test_df.iloc[0]["Sample"]
# repeat = test_df.iloc[0]["Repeat"]
# read = test_df.iloc[0]["Read"]
# mapped_gene = test_df.iloc[0]["MappedGene"]
# rtg_strand = test_df.iloc[0]["MappedStrand"]
# rtg_read_start = test_df.iloc[0]["ReadStart"]
# rtg_read_end = test_df.iloc[0]["ReadEnd"]
# rtg_gene_start = test_df.iloc[0]["GeneStart"]
# rtg_gene_end = test_df.iloc[0]["GeneEnd"]
# read_seq = test_df.iloc[0]["Seq"]

# find_first_n_pairwise_alignments_for_barcode_and_read(
#     sample,
#     gene,
#     repeat,
#     aligner,
#     barcode,  # query_name
#     read,  # target_name
#     mapped_gene,  # ADAR1, IQEC, Other, Unmapped
#     rtg_strand,  # the genomic strand to which the read was mapped to
#     rtg_read_start,
#     rtg_read_end,
#     rtg_gene_start,
#     rtg_gene_end,
#     barcode_seq,  # query_seq
#     read_seq,  # target_seq
#     first_n_alignments_per_strand=first_n_alignments_per_strand,
#     # debug_mode=True,
# )

# %%

# %%
# for strand in STRANDS:
#     print(f"\n{strand = }")

#     alignments = aligner.align(read_seq, barcode_seq, strand=strand)

#     num_of_alignments = len(alignments)
#     print(f"{num_of_alignments = }")
#     for alignment in alignments:
#         print(alignment)
#         score = alignment.score
#         gaps, identities, mismatches = alignment.counts()
#         alignment_length = alignment.length

#         btr_read_coords, btr_barcode_coords = alignment.aligned

#         prct_barcode_identity = 100 * identities / len(barcode)
#         prct_barcode_cov = 100 * (identities + mismatches) / len(barcode)

#         btr_read_coords, btr_barcode_coords = alignment.aligned

#         # num_of_btr_read_coords = len(btr_read_coords)
#         # num_of_btr_barcode_coords = len(btr_barcode_coords)
#         # num_of_btr_barcode_gap_openings = num_of_btr_barcode_coords - 1

#         print(
#             # f"{score = }, {gaps = }, {identities = }, {mismatches = }, "
#             # f"{alignment_length = } "
#             f"{prct_barcode_identity = :.2f}, {prct_barcode_cov = :.2f} "
#             f"{btr_read_coords = }, {btr_barcode_coords = }\n"
#         )

#     print()

# %%
# barcode_seq[4:12]

# %%
# barcode_seq[5:12]

# %%
# barcode_seq[5:12].reverse_complement()

# %%
# i = 1

# # ic(split_df.iloc[i])

# sample = split_df.iloc[i]["Sample"]
# repeat = split_df.iloc[i]["Repeat"]
# read = split_df.iloc[i]["Read"]
# mapped_gene = split_df.iloc[i]["MappedGene"]
# rtg_gene_strand = split_df.iloc[i]["MappedStrand"]
# rtg_read_start = split_df.iloc[i]["ReadStart"]
# rtg_read_end = split_df.iloc[i]["ReadEnd"]
# rtg_gene_start = split_df.iloc[i]["GeneStart"]
# rtg_gene_end = split_df.iloc[i]["GeneEnd"]
# read_seq = split_df.iloc[i]["Seq"]


# for strand in STRANDS:
#     print(f"\n{strand = }")

#     alignments = aligner.align(read_seq, barcode_seq, strand=strand)

#     num_of_alignments = len(alignments)
#     for alignment in alignments:
#         print(alignment)
#         score = alignment.score
#         gaps, identities, mismatches = alignment.counts()
#         alignment_length = alignment.length

#         btr_read_coords, btr_barcode_coords = alignment.aligned

#         prct_barcode_identity = 100 * identities / len(barcode)
#         prct_barcode_cov = 100 * (identities + mismatches) / len(barcode)

#         btr_read_coords, btr_barcode_coords = alignment.aligned

#         # num_of_btr_read_coords = len(btr_read_coords)
#         # num_of_btr_barcode_coords = len(btr_barcode_coords)
#         # num_of_btr_barcode_gap_openings = num_of_btr_barcode_coords - 1

#         print(
#             # f"{score = }, {gaps = }, {identities = }, {mismatches = }, "
#             # f"{alignment_length = } "
#             f"{prct_barcode_identity = :.2f}, {prct_barcode_cov = :.2f} "
#             f"{btr_read_coords = }, {btr_barcode_coords = }\n"
#         )

#     print()

# %%
# find_first_n_pairwise_alignments_for_barcode_and_read(
#     sample,
#     gene,
#     repeat,
#     aligner,
#     query_name,
#     target_name,
#     mapped_gene,  # ADAR1, IQEC, Other, Unmapped
#     genomic_target_strand,  # the genomic strand to which the read was mapped to
#     rtg_read_start,
#     rtg_read_end,
#     rtg_gene_start,
#     rtg_gene_end,
#     query_seq,
#     target_seq,
#     first_n_alignments_per_strand=10,
#     debug_mode=False,
# )

# %%
# alignments_dfs = split_df.apply(
#     lambda x: find_first_n_pairwise_alignments_for_barcode_and_read(
#         x["Sample"],
#         gene,
#         x["Repeat"],
#         aligner,
#         query_name,
#         x["Read"],  # target_name
#         # x["Mapped"],  # target_mapped_to_ref
#         x["MappedGene"],  # mapped_gene
#         x[
#             "MappedStrand"
#         ],  # genomic_target_strand - the genomic strand to which the read was mapped to
#         x["ReadStart"],
#         x["ReadEnd"],
#         x["GeneStart"],
#         x["GeneEnd"],
#         query_seq,
#         # x["SoftClippedSeq"],  # target_seq
#         x["Seq"],  # target_seq
#         first_n_alignments_per_strand=first_n_alignments_per_strand,
#         debug_mode=debug_mode,
#     ),
#     axis=1,
# )
# alignments_dfs[0]

# %%

# %%

# %%
def get_genomic_coord_for_read_coord(aligned_pairs, required_read_pos):
    for read_pos, ref_pos, cigar_op in aligned_pairs:
        if read_pos == required_read_pos:
            return ref_pos
    return None


# %%
def find_btg_gene_coords(input_df):

    annotated_dfs = []

    all_required_reads = input_df["Read"].values

    for bam_file in mapped_bam_files:

        with pysam.AlignmentFile(
            bam_file,
            "rb",
            threads=10,
        ) as samfile:
            reads = [read for read in samfile if read.query_name in all_required_reads]
            reads_names = [
                read.query_name for read in reads
            ]  # names of reads found in this specific bam file
            if len(reads_names) == 0:
                continue
            aligned_pairs = [
                read.get_aligned_pairs(matches_only=False, with_cigar=True)
                for read in reads
            ]
            annotated_df = input_df.loc[input_df["Read"].isin(reads_names)].copy()
            aligned_pairs_series = pd.Series(aligned_pairs, index=reads_names)
            annotated_df["BTGGeneStart"] = annotated_df.apply(
                lambda x: get_genomic_coord_for_read_coord(
                    aligned_pairs_series.loc[x["Read"]], x["BTRReadStart"]
                ),
                axis=1,
            )
            annotated_df["BTGGeneEnd"] = annotated_df.apply(
                lambda x: get_genomic_coord_for_read_coord(
                    aligned_pairs_series.loc[x["Read"]], x["BTRReadEnd"]
                ),
                axis=1,
            )
            annotated_dfs.append(annotated_df)

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


# get BTGeneStart and BTGGeneEnd for each read
concat_alignments_df = find_btg_gene_coords(concat_alignments_df)

concat_alignments_df

# %%
concat_alignments_df["BTRNumOfAllAlignmentsPerStrand"].describe()

# %%
read = "m64296e_241222_071206/100007945/ccs"

concat_alignments_df.loc[concat_alignments_df["Read"].eq(read)]

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

concat_alignments_df

# %%
{
    "TargetName": "Read",
    "QueryName": "Barcode",
    "TargetSeq": "ReadSeq",
    "GenomicTargetStrand": "RTGStrand",
    "TargetStrand": "BTRStrand",
    "%QueryIdentity": "%BTRBarcodeIdentity",
    "%QueryCoverage": "%BTRBarcodeCoverage",
    "TargetSeqLength": "ReadSeqLength",
    "QuerySeqLength": "BarcodeSeqLength",
    "AlignmentLength": "BTRAlignmentLength",
    "AlignedTargetCoords": "BTRReadCoords",
    "AlignedQueryCoords": "BTRBarcodeCoords",
    "AlignmentObject": "BTRAlignmentObject",
    "NumOfAlignedTargetCoords": "NumOfBTRReadCoords",
    "NumOfAlignedQueryCoords": "NumOfBTRBarcodeCoords",
    "QueryGapOpenings": "BTRBarcodeGapOpenings",
    "Score": "BTRScore",
    "Gaps": "BTRGaps",
    "Identitites": "BTRIdentitites",
    "Mismatches": "BTRMismatches",
}

# %%
concat_alignments_df[["BTRBarcodeGapOpenings", "BTRGaps"]].value_counts(
    normalize=True
).mul(100).round(4)

# %%
concat_alignments_df.loc[
    concat_alignments_df["BTRBarcodeGapOpenings"] > 0,
    [
        "BTRBarcodeGapOpenings",
        "BTRGaps",
    ],
].value_counts()

# %%
# concat_alignments_df.groupby(["Gene", "Repeat", "QueryName"])[
#     ["%QueryIdentity", "%QueryCoverage"]
# ].describe()

# %%
# concat_alignments_df.groupby(["Gene", "QueryName"])[
#     ["%QueryIdentity", "%QueryCoverage"]
# ].describe()

# %%
concat_alignments_df.loc[
    concat_alignments_df["Gene"].eq(concat_alignments_df["MappedGene"])
].groupby(
    [
        "Gene",
        "Barcode",
    ]
)[
    ["%BTRBarcodeIdentity", "%BTRBarcodeCoverage"]
].describe()

# %%
concat_alignments_df.loc[
    ~concat_alignments_df["Gene"].eq(concat_alignments_df["MappedGene"])
].groupby(
    [
        "Gene",
        "Barcode",
    ]
)[
    ["%BTRBarcodeIdentity", "%BTRBarcodeCoverage"]
].describe()

# %%
concat_alignments_df.groupby(["Gene", "Barcode", "MappedGene"])[
    ["%BTRBarcodeIdentity", "%BTRBarcodeCoverage"]
].describe()

# %%
concat_alignments_df.head()

# %%
# for gene in genes:

#     density_heatmap_concat_alignments_df = (
#         concat_alignments_df.loc[concat_alignments_df["Gene"] == gene]
#         .groupby(
#             [
#                 "Gene",
#                 "MappedGene",
#                 # "Repeat",
#                 "QueryName",
#                 "GenomicTargetStrand",
#                 "TargetStrand",
#                 "%QueryCoverage",
#                 "%QueryIdentity",
#             ]
#         )
#         .size()
#         .reset_index(name="Matches")
#     )
#     density_heatmap_concat_alignments_df["Log10(Matches)"] = np.log10(
#         density_heatmap_concat_alignments_df["Matches"]
#     )
#     density_heatmap_concat_alignments_df["Log2(Matches)"] = np.log2(
#         density_heatmap_concat_alignments_df["Matches"]
#     )

#     for acceptable_mapped_genes_options in [
#         ["ADAR1", "Other", "IQEC", "Unmapped"],
#         [gene],
#     ]:

#         fig = px.scatter(
#             density_heatmap_concat_alignments_df.loc[
#                 density_heatmap_concat_alignments_df["MappedGene"].isin(
#                     acceptable_mapped_genes_options
#                 )
#             ],
#             x="%QueryCoverage",
#             y="%QueryIdentity",
#             color="QueryName",
#             size="Matches",
#             facet_col="GenomicTargetStrand",
#             facet_col_spacing=0.05,
#             facet_row="TargetStrand",
#             facet_row_spacing=0.05,
#             labels={
#                 "%QueryCoverage": "Query coverage [%]",
#                 "%QueryIdentity": "Query identity [%]",
#                 "QueryName": "Tag",
#                 "GenomicTargetStrand": "Genomic read strand",
#                 "TargetStrand": "Relative read strand",
#             },
#         )
#         fig.update_xaxes(range=[0, 105], dtick=10)
#         fig.update_yaxes(range=[0, 105], dtick=10)
#         # Reduce opacity to see both histograms
#         fig.update_traces(opacity=0.7)
#         # Overlay both histograms
#         fig.update_layout(
#             # barmode="overlay",
#             template="plotly_white",
#             width=850,
#             height=600,
#             title=f"{gene} samples, reads mapped to {", ".join(acceptable_mapped_genes_options)}",
#             # showlegend=False
#         )

#         fig.show()
# # density_heatmap_concat_alignments_df

# %% [markdown]
# # Combine barcodes' alignments per read
#

# %%
# def combine_one_read_aligned_tags_per_strand(
#     one_read_aligned_primers_df,
#     gene,
#     repeat,
#     read_name,
# ):
#     if one_read_aligned_primers_df.index.names == ["TargetStrand", "QueryName"]:
#         # the index for the df is already set correctly - so there's nothing that should be done in this aspect
#         pass
#     elif {"TargetStrand", "QueryName"}.issubset(one_read_aligned_primers_df.columns):
#         one_read_aligned_primers_df = one_read_aligned_primers_df.set_index(
#             [
#                 "TargetStrand",
#                 "QueryName",
#             ]
#         ).sort_index()
#     else:
#         raise Error(
#             f"{gene=}, {repeat=}, {read_name=}: 'TargetStrand' and 'QueryName' should be both either columns or the only index levels"
#         )

#     sample = one_read_aligned_primers_df.at[("+", gene), "Sample"]
#     genomic_read_strand = one_read_aligned_primers_df.at[
#         ("+", gene), "GenomicTargetStrand"
#     ]
#     read_seq = one_read_aligned_primers_df.at[("+", gene), "TargetSeq"]

#     gene_lexically_smaller_than_pcr = gene < "PCR"

#     rows = []

#     for strand in ["+", "-"]:

#         row = [gene, repeat, sample, read_name, genomic_read_strand, strand, read_seq]

#         one_strand_one_read_aligned_primers_df = one_read_aligned_primers_df.loc[
#             (strand,),
#         ]

#         gene_tag_coords_on_read = one_strand_one_read_aligned_primers_df.at[
#             gene, "AlignedTargetCoords"
#         ]
#         pcr_tag_coords_on_read = one_strand_one_read_aligned_primers_df.at[
#             "PCR", "AlignedTargetCoords"
#         ]
#         # ic(gene_tag_coords_on_read, pcr_tag_coords_on_read)

#         tags_coords_overlap = False

#         first_gene_tag_coord_on_read = gene_tag_coords_on_read[0][0]
#         last_gene_tag_coord_on_read = gene_tag_coords_on_read[-1][-1]
#         first_pcr_tag_coord_on_read = pcr_tag_coords_on_read[0][0]
#         last_pcr_tag_coord_on_read = pcr_tag_coords_on_read[-1][-1]

#         if (
#             first_gene_tag_coord_on_read
#             < last_gene_tag_coord_on_read
#             < first_pcr_tag_coord_on_read
#             < last_pcr_tag_coord_on_read
#         ):
#             tags = [gene, "PCR"]
#         elif (
#             first_pcr_tag_coord_on_read
#             < last_pcr_tag_coord_on_read
#             < first_gene_tag_coord_on_read
#             < last_gene_tag_coord_on_read
#         ):
#             tags = ["PCR", gene]
#         else:
#             tags_coords_overlap = True
#             if (
#                 first_gene_tag_coord_on_read < last_gene_tag_coord_on_read
#                 and first_pcr_tag_coord_on_read < last_pcr_tag_coord_on_read
#             ):
#                 if first_gene_tag_coord_on_read < first_pcr_tag_coord_on_read:
#                     tags = [gene, "PCR"]
#                 elif first_pcr_tag_coord_on_read < first_gene_tag_coord_on_read:
#                     tags = ["PCR", gene]
#                 elif last_gene_tag_coord_on_read < last_pcr_tag_coord_on_read:
#                     tags = [gene, "PCR"]
#                 elif last_pcr_tag_coord_on_read < last_gene_tag_coord_on_read:
#                     tags = ["PCR", gene]
#                 else:
#                     raise ValueError(
#                         f"{gene=}, {repeat=}, {read_name=}, {strand=}: unexpected order of tags' coordinates w.r.t their mapping to the read. "
#                         f"{first_pcr_tag_coord_on_read=}, {last_pcr_tag_coord_on_read=}, {first_gene_tag_coord_on_read=}, {last_gene_tag_coord_on_read=}."
#                     )
#             else:
#                 raise ValueError(
#                     f"{gene=}, {repeat=}, {read_name=}, {strand=}: unexpected order of tags' coordinates w.r.t their mapping to the read. "
#                     f"{first_pcr_tag_coord_on_read=}, {last_pcr_tag_coord_on_read=}, {first_gene_tag_coord_on_read=}, {last_gene_tag_coord_on_read=}."
#                 )

#         if tags == [gene, "PCR"]:
#             tags_coords_row_extension = [
#                 tags,
#                 tags_coords_overlap,
#                 first_gene_tag_coord_on_read,
#                 last_gene_tag_coord_on_read,
#                 first_pcr_tag_coord_on_read,
#                 last_pcr_tag_coord_on_read,
#             ]
#         elif tags == ["PCR", gene]:
#             tags_coords_row_extension = [
#                 tags,
#                 tags_coords_overlap,
#                 first_pcr_tag_coord_on_read,
#                 last_pcr_tag_coord_on_read,
#                 first_gene_tag_coord_on_read,
#                 last_gene_tag_coord_on_read,
#             ]
#         else:
#             raise Exception()

#         row.extend(tags_coords_row_extension)

#         # sort the df for the next extenstion of the row
#         if (tags == [gene, "PCR"] and gene_lexically_smaller_than_pcr) or (
#             tags == ["PCR", gene] and not gene_lexically_smaller_than_pcr
#         ):
#             one_strand_one_read_aligned_primers_df = (
#                 one_strand_one_read_aligned_primers_df.sort_index()
#             )
#         elif (tags == [gene, "PCR"] and not gene_lexically_smaller_than_pcr) or (
#             tags == ["PCR", gene] and gene_lexically_smaller_than_pcr
#         ):
#             one_strand_one_read_aligned_primers_df = (
#                 one_strand_one_read_aligned_primers_df.sort_index(ascending=False)
#             )
#         else:
#             raise Exception()

#         row.extend(
#             one_strand_one_read_aligned_primers_df.loc[
#                 :,
#                 [
#                     "Score",
#                     "Gaps",
#                     "Identitites",
#                     "Mismatches",
#                     "%QueryIdentity",
#                     "%QueryCoverage",
#                     "TargetSeqLength",
#                     "QuerySeqLength",
#                     "AlignmentLength",
#                     "AlignedTargetCoords",
#                     "AlignedQueryCoords",
#                     "AlignmentObject",
#                     "NumOfAlignedTargetCoords",
#                     "NumOfAlignedQueryCoords",
#                     "QueryGapOpenings",
#                 ],
#             ].T.values.tolist()
#         )
#         rows.append(row)

#     combined_one_read_aligned_tags_per_strand_df = pd.DataFrame(
#         rows,
#         columns=[
#             "Gene",
#             "Repeat",
#             "Sample",
#             "TargetName",
#             "GenomicTargetStrand",
#             "RelativeTargetStrand",
#             "TargetSeq",
#             "Queries",
#             "QueriesCoordsOverlap",
#             "Query1Start",
#             "Query1End",
#             "Query2Start",
#             "Query2End",
#             "Score",
#             "Gaps",
#             "Identitites",
#             "Mismatches",
#             "%QueryIdentity",
#             "%QueryCoverage",
#             "TargetSeqLength",
#             "QuerySeqLength",
#             "AlignmentLength",
#             "AlignedTargetCoords",
#             "AlignedQueryCoords",
#             "AlignmentObject",
#             "NumOfAlignedTargetCoords",
#             "NumOfAlignedQueryCoords",
#             "QueryGapOpenings",
#         ],
#     )

#     combined_one_read_aligned_tags_per_strand_df.insert(
#         combined_one_read_aligned_tags_per_strand_df.columns.get_loc("Query2End") + 1,
#         "UMISeq",
#         combined_one_read_aligned_tags_per_strand_df.apply(
#             lambda x: x["TargetSeq"][x["Query1End"] : x["Query2Start"]], axis=1
#         ),
#     )
#     combined_one_read_aligned_tags_per_strand_df.loc[
#         combined_one_read_aligned_tags_per_strand_df["RelativeTargetStrand"].eq("-"),
#         "UMISeq",
#     ] = combined_one_read_aligned_tags_per_strand_df.loc[
#         combined_one_read_aligned_tags_per_strand_df["RelativeTargetStrand"].eq("-"),
#         "UMISeq",
#     ].apply(
#         reverse_complement
#     )
#     combined_one_read_aligned_tags_per_strand_df.insert(
#         combined_one_read_aligned_tags_per_strand_df.columns.get_loc("UMISeq") + 1,
#         "UMISeqLength",
#         combined_one_read_aligned_tags_per_strand_df["UMISeq"].apply(len),
#     )

#     assert (
#         combined_one_read_aligned_tags_per_strand_df.loc[
#             combined_one_read_aligned_tags_per_strand_df["QueriesCoordsOverlap"],
#             "UMISeq",
#         ]
#         .eq("")
#         .all()
#     )
#     assert (
#         combined_one_read_aligned_tags_per_strand_df.loc[
#             combined_one_read_aligned_tags_per_strand_df["QueriesCoordsOverlap"],
#             "UMISeqLength",
#         ]
#         .eq(0)
#         .all()
#     )

#     combined_one_read_aligned_tags_per_strand_df.insert(
#         combined_one_read_aligned_tags_per_strand_df.columns.get_loc("UMISeqLength")
#         + 1,
#         "GoodUMISeqLength",
#         combined_one_read_aligned_tags_per_strand_df["UMISeqLength"].isin([11, 12, 13]),
#     )
#     combined_one_read_aligned_tags_per_strand_df.insert(
#         combined_one_read_aligned_tags_per_strand_df.columns.get_loc("UMISeqLength")
#         + 2,
#         "MinUMISeqLength",
#         combined_one_read_aligned_tags_per_strand_df["UMISeqLength"].eq(
#             combined_one_read_aligned_tags_per_strand_df["UMISeqLength"].min()
#         ),
#     )
#     combined_one_read_aligned_tags_per_strand_df.insert(
#         combined_one_read_aligned_tags_per_strand_df.columns.get_loc("UMISeqLength")
#         + 3,
#         "MaxUMISeqLength",
#         combined_one_read_aligned_tags_per_strand_df["UMISeqLength"].eq(
#             combined_one_read_aligned_tags_per_strand_df["UMISeqLength"].max()
#         ),
#     )

#     combined_one_read_aligned_tags_per_strand_df.insert(
#         combined_one_read_aligned_tags_per_strand_df.columns.get_loc("%QueryIdentity")
#         + 1,
#         "Mean%QueryIdentity",
#         combined_one_read_aligned_tags_per_strand_df["%QueryIdentity"].apply(
#             lambda x: sum(x) / 2
#         ),
#     )
#     combined_one_read_aligned_tags_per_strand_df.insert(
#         combined_one_read_aligned_tags_per_strand_df.columns.get_loc("%QueryIdentity")
#         + 2,
#         "MaxMean%QueryIdentity",
#         combined_one_read_aligned_tags_per_strand_df["Mean%QueryIdentity"].eq(
#             combined_one_read_aligned_tags_per_strand_df["Mean%QueryIdentity"].max()
#         ),
#     )
#     combined_one_read_aligned_tags_per_strand_df.insert(
#         combined_one_read_aligned_tags_per_strand_df.columns.get_loc("%QueryCoverage")
#         + 1,
#         "Mean%QueryCoverage",
#         combined_one_read_aligned_tags_per_strand_df["%QueryCoverage"].apply(
#             lambda x: sum(x) / 2
#         ),
#     )
#     combined_one_read_aligned_tags_per_strand_df.insert(
#         combined_one_read_aligned_tags_per_strand_df.columns.get_loc("%QueryCoverage")
#         + 2,
#         "MaxMean%QueryCoverage",
#         combined_one_read_aligned_tags_per_strand_df["Mean%QueryCoverage"].eq(
#             combined_one_read_aligned_tags_per_strand_df["Mean%QueryCoverage"].max()
#         ),
#     )

#     return combined_one_read_aligned_tags_per_strand_df

# %%
list(product(STRANDS, [["PCR", gene] for gene in GENES]))


# %%
# for tag_names in [["PCR", "ADAR1"], ["PCR", "IQEC"]]:
#     ic(tag_names)
#     tag_names.remove("PCR")
#     gene_tag_name = tag_names[0]
#     tag_names.append("PCR")
#     ic(gene_tag_name, tag_names)

# %%
def combine_one_read_aligned_tags_per_strand(
    one_read_aligned_primers_df,
    sample_gene,  # the gene from which the sample was taken from
    repeat,
    read_name,
):
    if one_read_aligned_primers_df.index.names == ["TargetStrand", "QueryName"]:
        # the index for the df is already set correctly - so there's nothing that should be done in this aspect
        pass
    elif {"TargetStrand", "QueryName"}.issubset(one_read_aligned_primers_df.columns):
        one_read_aligned_primers_df = one_read_aligned_primers_df.set_index(
            [
                "TargetStrand",
                "QueryName",
            ]
        ).sort_index()
    else:
        raise Error(
            f"{sample_gene=}, {repeat=}, {read_name=}: 'TargetStrand' and 'QueryName' should be both either columns or the only index levels"
        )

    sample = one_read_aligned_primers_df.at[("+", sample_gene), "Sample"]
    mapped_gene = one_read_aligned_primers_df.at[("+", sample_gene), "MappedGene"]
    genomic_read_strand = one_read_aligned_primers_df.at[
        ("+", sample_gene), "GenomicTargetStrand"
    ]  # could be np.nan if the read was not mapped ("Unmapped" in the "MappedGene" column)
    read_seq = one_read_aligned_primers_df.at[("+", sample_gene), "TargetSeq"]
    read_seq_len = len(read_seq)

    rows = []

    # for strand in ["+", "-"]:

    # iterate over each of the 4 combinations of strand and PCR with each of the gene-specific tags
    for strand, tag_names in product(STRANDS, [["PCR", gene] for gene in GENES]):

        tag_names.remove("PCR")
        gene = tag_names[0]
        tag_names.append("PCR")

        gene_lexically_smaller_than_pcr = gene < "PCR"

        row = [
            sample_gene,
            repeat,
            sample,
            read_name,
            mapped_gene,
            genomic_read_strand,
            strand,
            read_seq,
        ]

        one_strand_one_gene_tag_one_read_aligned_primers_df = (
            one_read_aligned_primers_df.loc[(strand,),].loc[tag_names]
        )

        gene_tag_coords_on_read = (
            one_strand_one_gene_tag_one_read_aligned_primers_df.at[
                gene, "AlignedTargetCoords"
            ]
        )
        pcr_tag_coords_on_read = one_strand_one_gene_tag_one_read_aligned_primers_df.at[
            "PCR", "AlignedTargetCoords"
        ]
        # ic(gene_tag_coords_on_read, pcr_tag_coords_on_read)

        tags_coords_overlap = False

        first_gene_tag_coord_on_read = gene_tag_coords_on_read[0][0]
        last_gene_tag_coord_on_read = gene_tag_coords_on_read[-1][-1]
        first_pcr_tag_coord_on_read = pcr_tag_coords_on_read[0][0]
        last_pcr_tag_coord_on_read = pcr_tag_coords_on_read[-1][-1]

        if (
            first_gene_tag_coord_on_read
            < last_gene_tag_coord_on_read
            < first_pcr_tag_coord_on_read
            < last_pcr_tag_coord_on_read
        ):
            tags = [gene, "PCR"]
        elif (
            first_pcr_tag_coord_on_read
            < last_pcr_tag_coord_on_read
            < first_gene_tag_coord_on_read
            < last_gene_tag_coord_on_read
        ):
            tags = ["PCR", gene]
        else:
            tags_coords_overlap = True
            if (
                first_gene_tag_coord_on_read < last_gene_tag_coord_on_read
                and first_pcr_tag_coord_on_read < last_pcr_tag_coord_on_read
            ):
                if first_gene_tag_coord_on_read < first_pcr_tag_coord_on_read:
                    tags = [gene, "PCR"]
                elif first_pcr_tag_coord_on_read < first_gene_tag_coord_on_read:
                    tags = ["PCR", gene]
                elif last_gene_tag_coord_on_read < last_pcr_tag_coord_on_read:
                    tags = [gene, "PCR"]
                elif last_pcr_tag_coord_on_read < last_gene_tag_coord_on_read:
                    tags = ["PCR", gene]
                elif (
                    first_gene_tag_coord_on_read == first_pcr_tag_coord_on_read
                    and last_gene_tag_coord_on_read == last_pcr_tag_coord_on_read
                ):
                    tags = [gene, "PCR"]
                    np.random.default_rng().shuffle(tags)  # randomly shuffle the tags
                else:
                    raise ValueError(
                        f"{sample_gene=}, {repeat=}, {read_name=}, {gene=}, {strand=}: unexpected order of tags' coordinates w.r.t their mapping to the read. "
                        f"{first_pcr_tag_coord_on_read=}, {last_pcr_tag_coord_on_read=}, {first_gene_tag_coord_on_read=}, {last_gene_tag_coord_on_read=}."
                    )
            else:
                raise ValueError(
                    f"{sample_gene=}, {repeat=}, {read_name=}, {gene=}, {strand=}: unexpected order of tags' coordinates w.r.t their mapping to the read. "
                    f"{first_pcr_tag_coord_on_read=}, {last_pcr_tag_coord_on_read=}, {first_gene_tag_coord_on_read=}, {last_gene_tag_coord_on_read=}."
                )

        if tags == [gene, "PCR"]:
            tags_coords_row_extension = [
                tags,
                tags_coords_overlap,
                first_gene_tag_coord_on_read,
                last_gene_tag_coord_on_read,
                first_pcr_tag_coord_on_read,
                last_pcr_tag_coord_on_read,
            ]
        elif tags == ["PCR", gene]:
            tags_coords_row_extension = [
                tags,
                tags_coords_overlap,
                first_pcr_tag_coord_on_read,
                last_pcr_tag_coord_on_read,
                first_gene_tag_coord_on_read,
                last_gene_tag_coord_on_read,
            ]
        else:
            raise Exception()

        row.extend(tags_coords_row_extension)

        # sort the df for the next extenstion of the row
        if (tags == [gene, "PCR"] and gene_lexically_smaller_than_pcr) or (
            tags == ["PCR", gene] and not gene_lexically_smaller_than_pcr
        ):
            one_strand_one_gene_tag_one_read_aligned_primers_df = (
                one_strand_one_gene_tag_one_read_aligned_primers_df.sort_index()
            )
        elif (tags == [gene, "PCR"] and not gene_lexically_smaller_than_pcr) or (
            tags == ["PCR", gene] and gene_lexically_smaller_than_pcr
        ):
            one_strand_one_gene_tag_one_read_aligned_primers_df = (
                one_strand_one_gene_tag_one_read_aligned_primers_df.sort_index(
                    ascending=False
                )
            )
        else:
            raise Exception()

        row.extend(
            one_strand_one_gene_tag_one_read_aligned_primers_df.loc[
                :,
                [
                    "Score",
                    "Gaps",
                    "Identitites",
                    "Mismatches",
                    "%QueryIdentity",
                    "%QueryCoverage",
                    "TargetSeqLength",
                    "QuerySeqLength",
                    "AlignmentLength",
                    "AlignedTargetCoords",
                    "AlignedQueryCoords",
                    "AlignmentObject",
                    "NumOfAlignedTargetCoords",
                    "NumOfAlignedQueryCoords",
                    "QueryGapOpenings",
                ],
            ].T.values.tolist()
        )
        rows.append(row)

    combined_one_read_aligned_tags_per_strand_df = pd.DataFrame(
        rows,
        columns=[
            "Gene",
            "Repeat",
            "Sample",
            "TargetName",
            "MappedGene",
            "GenomicTargetStrand",
            "RelativeTargetStrand",
            "TargetSeq",
            "Queries",
            "QueriesCoordsOverlap",
            "Query1Start",
            "Query1End",
            "Query2Start",
            "Query2End",
            "Score",
            "Gaps",
            "Identitites",
            "Mismatches",
            "%QueryIdentity",
            "%QueryCoverage",
            "TargetSeqLength",
            "QuerySeqLength",
            "AlignmentLength",
            "AlignedTargetCoords",
            "AlignedQueryCoords",
            "AlignmentObject",
            "NumOfAlignedTargetCoords",
            "NumOfAlignedQueryCoords",
            "QueryGapOpenings",
        ],
    )

    combined_one_read_aligned_tags_per_strand_df.insert(
        combined_one_read_aligned_tags_per_strand_df.columns.get_loc("Query2End") + 1,
        "UMISeq",
        combined_one_read_aligned_tags_per_strand_df.apply(
            lambda x: x["TargetSeq"][x["Query1End"] : x["Query2Start"]], axis=1
        ),
    )
    combined_one_read_aligned_tags_per_strand_df.loc[
        combined_one_read_aligned_tags_per_strand_df["RelativeTargetStrand"].eq("-"),
        "UMISeq",
    ] = combined_one_read_aligned_tags_per_strand_df.loc[
        combined_one_read_aligned_tags_per_strand_df["RelativeTargetStrand"].eq("-"),
        "UMISeq",
    ].apply(
        reverse_complement
    )
    combined_one_read_aligned_tags_per_strand_df.insert(
        combined_one_read_aligned_tags_per_strand_df.columns.get_loc("UMISeq") + 1,
        "UMISeqLength",
        combined_one_read_aligned_tags_per_strand_df["UMISeq"].apply(len),
    )

    assert (
        combined_one_read_aligned_tags_per_strand_df.loc[
            combined_one_read_aligned_tags_per_strand_df["QueriesCoordsOverlap"],
            "UMISeq",
        ]
        .eq("")
        .all()
    )
    assert (
        combined_one_read_aligned_tags_per_strand_df.loc[
            combined_one_read_aligned_tags_per_strand_df["QueriesCoordsOverlap"],
            "UMISeqLength",
        ]
        .eq(0)
        .all()
    )

    combined_one_read_aligned_tags_per_strand_df.insert(
        combined_one_read_aligned_tags_per_strand_df.columns.get_loc("UMISeqLength")
        + 1,
        "%RelLocOfQuery1Start",
        combined_one_read_aligned_tags_per_strand_df.apply(
            lambda x: 100 * x["Query1Start"] / read_seq_len, axis=1
        ),
    )
    combined_one_read_aligned_tags_per_strand_df.insert(
        combined_one_read_aligned_tags_per_strand_df.columns.get_loc("UMISeqLength")
        + 2,
        "%RelLocOfQuery2End",
        # combined_one_read_aligned_tags_per_strand_df.apply(
        #     lambda x: 100 * (read_seq_len - x["Query2End"]) / read_seq_len, axis=1
        # ),
        combined_one_read_aligned_tags_per_strand_df.apply(
            lambda x: 100 * x["Query2End"] / read_seq_len, axis=1
        ),
    )

    combined_one_read_aligned_tags_per_strand_df.insert(
        combined_one_read_aligned_tags_per_strand_df.columns.get_loc("%QueryIdentity")
        + 1,
        "Mean%QueryIdentity",
        combined_one_read_aligned_tags_per_strand_df["%QueryIdentity"].apply(
            lambda x: sum(x) / 2
        ),
    )
    combined_one_read_aligned_tags_per_strand_df.insert(
        combined_one_read_aligned_tags_per_strand_df.columns.get_loc("%QueryIdentity")
        + 2,
        "MaxMean%QueryIdentity",
        combined_one_read_aligned_tags_per_strand_df["Mean%QueryIdentity"].eq(
            combined_one_read_aligned_tags_per_strand_df["Mean%QueryIdentity"].max()
        ),
    )
    combined_one_read_aligned_tags_per_strand_df.insert(
        combined_one_read_aligned_tags_per_strand_df.columns.get_loc("%QueryCoverage")
        + 1,
        "Mean%QueryCoverage",
        combined_one_read_aligned_tags_per_strand_df["%QueryCoverage"].apply(
            lambda x: sum(x) / 2
        ),
    )
    combined_one_read_aligned_tags_per_strand_df.insert(
        combined_one_read_aligned_tags_per_strand_df.columns.get_loc("%QueryCoverage")
        + 2,
        "MaxMean%QueryCoverage",
        combined_one_read_aligned_tags_per_strand_df["Mean%QueryCoverage"].eq(
            combined_one_read_aligned_tags_per_strand_df["Mean%QueryCoverage"].max()
        ),
    )

    return combined_one_read_aligned_tags_per_strand_df


# %%
gene = "ADAR1"
repeat = "1"
read_name = "m64296e_241222_071206/125042831/ccs"

# %%
# concat_bams_df

# %%
concat_bams_df.loc[
    (concat_bams_df["Gene"] == gene)
    & (concat_bams_df["Repeat"] == repeat)
    & (concat_bams_df["Read"] == read_name),
]

# %%
adar_gene_seq = "ATGGCGAACTCTAATCTATCTTCTTCCATGAATCAAAACATGGCTCACATGGGTGGTGGGAGTTTAGTAAATGGCTACTACAAACAAGTACCATATTCAGGTGGAAGAAGCCGAAATGCAAGTGGTAGGTCCGGTAGTCGTGGCCGCGGCAAACCTGCAGTCAGGGAAACTAGTTTGAATGTTCATCCAGAATGGGAAGAGCGTATTGTAAACTACCTTGCCCACAAAACTCATCCTGTAAAGACCATGGAGCTGGCACGCCTCGTGAATGTTCGCTCACGCAAAGAAGTGAATCCCACTTTGTACAGCATGGACAGACGAGGTCTTATCAGGAAACATGGAATGCAACCTCCGACATGGGTAATTGCAGACCCACCCCAATCTCACGGCGGATACAACCAAAATGAGACACACTATTCAAGTAGCCCAGGAATTTACCAGCATAGTCCGGTTTCGAGAACTCCTCAGAACTTTTATCCCAATAATCGAGAGAGTTATCGAGGACACAAAGCTCCAAATAGTAATTACCCGCGCTCCAAACGGACTTCATACAGGAATGACTGGCATAACTTTTGCTCCCCTCCATCTCACATGTACCCAGAGGGCAAAAATGAATCTTTGATCTATAGTCACAGTAACAAAGATAATGAGATGTTATCAATGGGAAACGCTAGTTCTCCAAACAGATTGCGGTCTGAAAGTTGTAGCCCAGATGAATGTGAAGCGCTGGAAACAGCCATTGCAAATTATCCACATAGTGGTGGGTATGGCCAGGGTTACTCAGGACATTTCCCTGTAACCCCAGAACCAAAGCAAACTGGCAAGAGAAGGAGAAATTGTGATAACTTTGGTTTACAACAACATCCATCAAATGGTACGATGCCGATGAAAAACAGTGAGAAGATCCAACAGAAAAAATTGGAATTTCAGGATGAAAGATATTATGATGCAAGCTACCCATCCTATTCTGGAAATTTTACTATGAACTATGCAAATCATTTTGATATGCCTGTCTATCATCCGATAGACCGGGAGGACCCGAAGGATACTCCACCTCCGTCACGTGTGTCAATGGACTTGATAAAAAAGGATTCCAAGGACATCTCGTCACATGAACGAATCTCTCCCAAGAGGAATTCAAACAGTAAGGGTTGTAATTCTGACGCTCATGACCCACAGGCAAGAGTTATTTCCTTCCTGGATAAAACTATGAATGGCTCTGCGAAGTCACGAGAAATCGCCAAGCATACAAGTCTTTCTCTGGAAGATACCCAGAAGATATTGCATAGTTTGTGTAAGAAGAAAATAGTCGCAACAATTGGTGAAGATATCTACATAATGGCTAAAAATGCAGCCAGTTATGAGACTGAGATTCCAGCAGGAAAAAACTCCTCATCAAACATGAATTCAAACATGGCACGCCAGTTCTCCAGTGGAAATCGGCAGCCCCCTGCGCCCCCACATGTACTATTGGCGGATAATGGCATCAATTCCGGCAGCATGAAAAACGTTTATTTCCAGGGTAATAATGCTCCCAAACAATCTGGGTCCAACTCGAGTGAATCAAAATCAGCACAGAGCCAGGTGGGCAGAAGCCCTCATCTACCCCCTTCCCCTCATGAACTATTAGCAAAGGACCCAATTTTCAAGGGAGACATTACTGCACCCAATACAAACGCTTCAAAGGACTACAACCAGTCGTCATCATCTTCGTCAGCATCCTTGTCGTCCTCAACGTCAAAGAACTCAAGGTGGAATAGCAACACTGCAGCGACAGAGAGTTCCAGAGCTCCAAACACGACCTCTGCTTCAACATCGTCAACTACATCATTTGCTCCCACTCCTAGTAAGTCTGCCTCTAATTCAAAACAGACTGCTCCTAGTCCCAAGCAACCATCTCCAAGTCCTAAGCAGAACACCCCTAAGAGTTCCAAGAGTTCCAAAAGTTCCAAGCAGAGAGCCACAAGCCCCAAACAAAACAGCACTCCTAGCTCCCAGGCGTCCTCTCAGTCAAACTCCAATACTACTACAACTGCCACCTCAAGCAGCAGCAAAAATAATAAAAATAACAACAATAACAACACCTCAGTAGAGAATTTGCAAGATGCCCTCAAAAATGTGTCTATCTCGTCCCCGACTGAGACTACTGAGAGCAAAACGCCCACATTGGCCGAGATCAAGGCGGCAGCAGTGGCGCAGGCGTTGGCCGACAAAGCGGCTGAGAAAGGAGCTGACAAGTCTGGTACTGATTCATTGGCACCAAATCTACAGATCACCTCAGAAAGTTTCGCCGCTCTCAACAAAAATCCAGTTAGCGCACTGATGGAATATGCTCAACAGCGACACTTACCCGTTGAATTTAAGCTTTTGTCACACAGAGGACCTTCTCATCGACCGTTGTTCAAATTTGCCGTGATTCTTGGTAAACGCCAGTTCCCCAGTATGGAGTGCAACAGTAAGAAGGATGGTAAGAAAGAGGCAGCCGATCTGACATTGCGCATTCTCATTGCTGAAGGACAGTATCAACTGGAGAACACCGTCTCAGCATTGAAAACAATTCCACCTGCTGAAATGACACATTTCGACAAAATGGCTGCCTTGAGTCACCAGGCATTTAACAACATTGCCTTGCAAATCCCTGAGAACCTTGCTGGGAGAAAGGTCATCGCTGCTTTGGTGATGAAGCGATCACCAACGGATACGGGAATTGTTATCAGTGTTGGAACTGGTAACCGCTGTTTAACCGGTGATCATTTGAGTTTGGAAGGCAACAGTGTCAATGACTCTCATGCTGAAATAATCACACGCCGAGGTTTTCTGAGATATCTGTACAAACATTTACTGGAGTATGATCCCGAAAAACCCCATGACCTATTTGAGAAAGGTGAACGTAGTCTTTGCCGGATAAAAACCAACATTACATTCCATCTGTATATATCAACTGCTCCTTGTGGTGATGGAGCACTTTTTTCACCCAGGGATACCGACTCCAGTAATGTGAAAGTGGATGAGGAAAATAAGCACGTCCATAATCCGACTTTTTCAAGCAGTGTTCAGGGATTGCTGAGAACCAAAGTGGAAGGAGGTGAAGGGACCATTCCAATAGATGCTGATTTCACTGAACAAACATGGGATGGAATTCAACGAGGTGAAAGATTGCGCACAATGTCATGTTCAGATAAAATATGTCGATGGAACGTTGTTGGTCTGCAAGGAGCTTTGCTTAGTCACTTTGTGGAACCAATCTACCTGGAATCTCTGACATTAGGTTATCTTTATGATCATGGCCACTTAGCACGAGCTGTTTGCTGCCGTATTGAACGGGGAGAGGCCTCTGTCAACCAACTACTACCTGAGGGCTACCGATTGAACCATCCTTGGCTTGGCAGAGTTACTGCTTGTGATCCACCTAGAGAAACCCAAAAGACGAAATCGTTAAGTATCAACTGGTGCTATGATGATGAAAAGTCTGAAGTTCTCGATGGTACAGCAGGCATCTGTTACACAGCGATTGAGAAAAATCTCTTCTCTCGCTTAACAAAGCACAGCTTATATGAAGAATTCAAAAAAGTGTGTCAGAAATTTGAACGCGAGGACTTGATGAATGTCACTTCTTACAACAAAGCCAAGATGATGGCCATTCCTTTCCAGACTGCCAAAAACGTAATGTTGAAAAAACTCAAGGAAAACAACTGCGGAACTTGGGTGTCAAAACCTATTGAGGAGGAGATGTTTACGTAG"

# %%
raw_seq = concat_bams_df.loc[
    (concat_bams_df["Gene"] == gene)
    & (concat_bams_df["Repeat"] == repeat)
    & (concat_bams_df["Read"] == read_name),
    "Seq",
].values[0]

# %%
aligner = Align.PairwiseAligner(
    # match_score=1.0, mismatch_score=-2.0, gap_score = -2.5,
    mode="local",  # otherwise we'd get scattered matches of the primer across the read
    scoring="blastn",
    # scoring="megablast"
)

# %%
primers_dict["PCR"]

# %%
len(primers_dict["IQEC"])

# %%
alignments = aligner.align(raw_seq, primers_dict["PCR"], strand="+")
ic(len(alignments))
alignment = alignments[0]
print(alignment)

# %%
alignments = aligner.align(raw_seq, primers_dict["PCR"], strand="-")
ic(len(alignments))
alignment = alignments[0]
print(alignment)

# %%
alignments = aligner.align(adar_gene_seq, raw_seq, strand="+")

# %%
len(alignments)

# %%
alignment = alignments[0]
print(alignment)

# %%
gene = "ADAR1"
repeat = "1"
read_name = "m64296e_241222_071206/125042831/ccs"

one_read_aligned_primers_df = (
    concat_alignments_df.loc[
        (concat_alignments_df["Gene"] == gene)
        & (concat_alignments_df["Repeat"] == repeat)
        & (concat_alignments_df["TargetName"] == read_name)
    ]
    .set_index(
        [
            "TargetStrand",
            "QueryName",
        ]
    )
    .sort_index()
)
one_read_aligned_primers_df

# %%
print(one_read_aligned_primers_df.at[("+", gene), "AlignmentObject"])

# %%
print(one_read_aligned_primers_df.at[("+", "PCR"), "AlignmentObject"])

# %%
print(one_read_aligned_primers_df.at[("-", gene), "AlignmentObject"])

# %%
print(one_read_aligned_primers_df.at[("-", "PCR"), "AlignmentObject"])

# %%
combined_one_read_aligned_tags_per_strand_df = combine_one_read_aligned_tags_per_strand(
    one_read_aligned_primers_df,
    gene,
    repeat,
    read_name,
)

combined_one_read_aligned_tags_per_strand_df

# %%
# combined_one_read_aligned_tags_per_strand_df.loc[
#     (combined_one_read_aligned_tags_per_strand_df["MaxMean%QueryIdentity"])
#     & (combined_one_read_aligned_tags_per_strand_df["MaxMean%QueryCoverage"])
#     & (combined_one_read_aligned_tags_per_strand_df["GoodUMISeqLength"])
# ]

# %%
gene = "ADAR1"
repeat = "1"
sample = "LP1ADAR1"
read_name = "m64296e_241222_071206/100007961/ccs"

one_read_aligned_primers_df = (
    concat_alignments_df.loc[
        (concat_alignments_df["Gene"] == gene)
        & (concat_alignments_df["Repeat"] == repeat)
        & (concat_alignments_df["TargetName"] == read_name)
    ]
    .set_index(
        [
            "TargetStrand",
            "QueryName",
        ]
    )
    .sort_index()
)
one_read_aligned_primers_df

# %%
strand = "+"
tag_names = ["PCR", gene]
one_strand_one_gene_tag_one_read_aligned_primers_df = one_read_aligned_primers_df.loc[
    (strand,),
].loc[tag_names]
one_strand_one_gene_tag_one_read_aligned_primers_df

# %%
# print(one_read_aligned_primers_df.at[("+", "ADAR1"), "AlignmentObject"])

# %%
# print(one_read_aligned_primers_df.at[("+", "PCR"), "AlignmentObject"])

# %%
# print(one_read_aligned_primers_df.at[("-", "ADAR1"), "AlignmentObject"])

# %%
# print(one_read_aligned_primers_df.at[("-", "PCR"), "AlignmentObject"])

# %%
# "ATGGCCAT" in "GCAGTCTGGAAAGGAATGGCCATC"

# %%
combined_one_read_aligned_tags_per_strand_df = combine_one_read_aligned_tags_per_strand(
    one_read_aligned_primers_df,
    gene,
    repeat,
    read_name,
)

combined_one_read_aligned_tags_per_strand_df

# %%
# combined_one_read_aligned_tags_per_strand_df.loc[
#     (combined_one_read_aligned_tags_per_strand_df["MaxMean%QueryIdentity"])
#     & (combined_one_read_aligned_tags_per_strand_df["MaxMean%QueryCoverage"])
#     & (combined_one_read_aligned_tags_per_strand_df["GoodUMISeqLength"])
# ]

# %%
# ValueError: sample_gene='ADAR1', repeat='3', read_name='m64296e_241222_071206/84214648/ccs',
# gene='ADAR1', strand='+':
# unexpected order of tags' coordinates w.r.t their mapping to the read.
# first_pcr_tag_coord_on_read=np.int64(20),
# last_pcr_tag_coord_on_read=np.int64(26),
# first_gene_tag_coord_on_read=np.int64(20),
# last_gene_tag_coord_on_read=np.int64(26).

gene = "ADAR1"
repeat = "3"
read_name = "m64296e_241222_071206/84214648/ccs"

one_read_aligned_primers_df = (
    concat_alignments_df.loc[
        (concat_alignments_df["Gene"] == gene)
        & (concat_alignments_df["Repeat"] == repeat)
        & (concat_alignments_df["TargetName"] == read_name)
    ]
    .set_index(
        [
            "TargetStrand",
            "QueryName",
        ]
    )
    .sort_index()
)
one_read_aligned_primers_df

# %%
combined_one_read_aligned_tags_per_strand_df = combine_one_read_aligned_tags_per_strand(
    one_read_aligned_primers_df,
    gene,
    repeat,
    read_name,
)

combined_one_read_aligned_tags_per_strand_df

# %%
# ValueError: sample_gene='ADAR1', repeat='3', read_name='m64296e_241222_071206/84214648/ccs',
# gene='ADAR1', strand='+':
# unexpected order of tags' coordinates w.r.t their mapping to the read.
# first_pcr_tag_coord_on_read=np.int64(20),
# last_pcr_tag_coord_on_read=np.int64(26),
# first_gene_tag_coord_on_read=np.int64(20),
# last_gene_tag_coord_on_read=np.int64(26).

gene = "ADAR1"
repeat = "1"
read_name = "m64296e_241222_071206/100008482/ccs"

one_read_aligned_primers_df = (
    concat_alignments_df.loc[
        (concat_alignments_df["Gene"] == gene)
        & (concat_alignments_df["Repeat"] == repeat)
        & (concat_alignments_df["TargetName"] == read_name)
    ]
    .set_index(
        [
            "TargetStrand",
            "QueryName",
        ]
    )
    .sort_index()
)
one_read_aligned_primers_df

# %%
combined_one_read_aligned_tags_per_strand_df = combine_one_read_aligned_tags_per_strand(
    one_read_aligned_primers_df,
    gene,
    repeat,
    read_name,
)

combined_one_read_aligned_tags_per_strand_df

# %%

# %%
concat_alignments_df.head()

# %%
concat_alignments_df.loc[:, ["Gene", "Repeat", "TargetName"]].drop_duplicates().shape[0]

# %%
concat_alignments_df = concat_alignments_df.sort_values(
    ["Gene", "Repeat", "TargetName"], ignore_index=True
)
concat_alignments_df

# %%
unique_gene_repeat_read_combintations = (
    concat_alignments_df.loc[:, ["Gene", "Repeat", "TargetName"]]
    .drop_duplicates()
    .values.tolist()
)

ic(unique_gene_repeat_read_combintations[:3])

# %%
# concat_alignments_df.head(7)

# %%
rows_per_read = 6

one_read_aligned_primers_dfs = [
    concat_alignments_df.iloc[i : i + rows_per_read]
    for i in range(0, len(concat_alignments_df), rows_per_read)
]
ic(len(one_read_aligned_primers_dfs))
one_read_aligned_primers_dfs[0]

# %%
assert len(one_read_aligned_primers_dfs) == len(unique_gene_repeat_read_combintations)

combine_one_read_aligned_tags_per_strand_inputs = [
    (
        one_read_aligned_primers_df,
        gene,
        repeat,
        read_name,
    )
    for one_read_aligned_primers_df, (gene, repeat, read_name) in zip(
        one_read_aligned_primers_dfs, unique_gene_repeat_read_combintations, strict=True
    )
]

ic(len(combine_one_read_aligned_tags_per_strand_inputs))

# combine_one_read_aligned_tags_per_strand_inputs[0]

# %% [markdown]
#

# %%
with Pool(processes=7) as pool:
    combined_one_read_aligned_tags_per_strand_dfs = pool.starmap(
        func=combine_one_read_aligned_tags_per_strand,
        iterable=combine_one_read_aligned_tags_per_strand_inputs,
    )
concat_combined_one_read_aligned_tags_per_strand_df = pd.concat(
    combined_one_read_aligned_tags_per_strand_dfs, ignore_index=True
)

concat_combined_one_read_aligned_tags_per_strand_df.insert(
    concat_combined_one_read_aligned_tags_per_strand_df.columns.get_loc("Queries") + 1,
    "StrQueries",
    concat_combined_one_read_aligned_tags_per_strand_df["Queries"].apply(str),
)

concat_combined_one_read_aligned_tags_per_strand_df

# %%
df = concat_combined_one_read_aligned_tags_per_strand_df.copy()
# df["StrQueries"] = df["Queries"].apply(sorted).apply(lambda x: str(x))
# df["StrQueries"] = df["Queries"].apply(lambda x: str(x))
df["GenomicStrand/ReadStrand"] = (
    df["GenomicTargetStrand"] + "/" + df["RelativeTargetStrand"]
)

# %%
# df.groupby(["Gene", "GenomicTargetStrand", "RelativeTargetStrand", "StrQueries"])[
#     [
#         "%RelLocOfQuery1Start",
#         "%RelLocOfQuery2End",
#     ]
# ].describe()

# %%
df.groupby(["Gene", "StrQueries"])[
    [
        "%RelLocOfQuery1Start",
        "%RelLocOfQuery2End",
    ]
].describe()

# %%
df.loc[df["Gene"] == df["MappedGene"]].groupby(["Gene", "StrQueries"])[
    [
        "%RelLocOfQuery1Start",
        "%RelLocOfQuery2End",
    ]
].describe().round(1)

# %%
density_heatmap_concat_alignments_df = (
    df.loc[df["Gene"] == df["MappedGene"]]
    .groupby(
        [
            "Gene",
            "StrQueries",
            "GenomicStrand/ReadStrand",
            # "GenomicTargetStrand",
            # "RelativeTargetStrand",
            "Mean%QueryCoverage",
            "Mean%QueryIdentity",
            "%RelLocOfQuery1Start",
            "%RelLocOfQuery2End",
        ]
    )
    .size()
    .reset_index(name="Matches")
)
density_heatmap_concat_alignments_df["Log10(Matches)"] = np.log10(
    density_heatmap_concat_alignments_df["Matches"]
)
density_heatmap_concat_alignments_df

# %%
density_heatmap_concat_alignments_df["Matches"].describe()

# %%
density_heatmap_concat_alignments_df["Log10(Matches)"].describe()

# %%
fig = px.scatter(
    density_heatmap_concat_alignments_df,
    x="Mean%QueryCoverage",
    y="Mean%QueryIdentity",
    color="StrQueries",
    # size="Matches",
    size="Log10(Matches)",
    facet_col="GenomicStrand/ReadStrand",
    facet_col_spacing=0.05,
    # facet_col_wrap=2,
    facet_row="Gene",
    facet_row_spacing=0.05,
    labels={
        "Mean%QueryCoverage": "Mean queries' coverage [%]",
        "Mean%QueryIdentity": "Mean queries' identity [%]",
        "StrQueries": "Queries",
        # "GenomicTargetStrand": "Genomic read strand",
        # "RelativeTargetStrand": "Relative read strand",
    },
    render_mode="webgl",
)
fig.update_xaxes(range=[0, 105], dtick=10)
fig.update_yaxes(range=[0, 105], dtick=10)
# Reduce opacity to see both histograms
fig.update_traces(opacity=0.7)
# Overlay both histograms
fig.update_layout(
    # barmode="overlay",
    template="plotly_white",
    width=1200,
    height=600,
    # title=f"{gene} samples, reads mapped to {", ".join(acceptable_mapped_genes_options)}",
    # showlegend=False
)

fig.show()

# %%
fig = px.histogram(
    df.loc[df["Gene"] == df["MappedGene"]],
    x="Mean%QueryCoverage",
    y="Mean%QueryIdentity",
    color="StrQueries",
    # size="Matches",
    # size="Log10(Matches)",
    histfunc="avg",
    facet_col="GenomicStrand/ReadStrand",
    facet_col_spacing=0.05,
    # facet_col_wrap=2,
    facet_row="Gene",
    facet_row_spacing=0.05,
    labels={
        "Mean%QueryCoverage": "Mean queries' coverage [%]",
        "Mean%QueryIdentity": "Mean queries' identity [%]",
        "StrQueries": "Queries",
        # "GenomicTargetStrand": "Genomic read strand",
        # "RelativeTargetStrand": "Relative read strand",
    },
    # render_mode="webgl",
)
fig.update_xaxes(range=[0, 105], dtick=10)
fig.update_yaxes(range=[0, 105], dtick=10)
# Reduce opacity to see both histograms
fig.update_traces(opacity=0.7)
# Overlay both histograms
fig.update_layout(
    barmode="overlay",
    template="plotly_white",
    width=1200,
    height=650,
    # title=f"{gene} samples, reads mapped to {", ".join(acceptable_mapped_genes_options)}",
    # showlegend=False
)

fig.show()

# %%

# %%
# for gene in genes:

#     density_heatmap_concat_alignments_df = (
#         df.loc[df["Gene"] == gene]
#         .groupby(
#             [
#                 "Gene",
#                 "MappedGene",
#                 # "Repeat",
#                 "StrQueries",
#                 "GenomicTargetStrand",
#                 "RelativeTargetStrand",
#                 "Mean%QueryCoverage",
#                 "Mean%QueryIdentity",
#                 "%RelLocOfQuery1Start",
#                 "%RelLocOfQuery2End",
#             ]
#         )
#         .size()
#         .reset_index(name="Matches")
#     )
#     density_heatmap_concat_alignments_df["Log10(Matches)"] = np.log10(
#         density_heatmap_concat_alignments_df["Matches"]
#     )
#     density_heatmap_concat_alignments_df["Log2(Matches)"] = np.log2(
#         density_heatmap_concat_alignments_df["Matches"]
#     )

#     for acceptable_mapped_genes_options in [
#         ["ADAR1", "Other", "IQEC", "Unmapped"],
#         [gene],
#     ]:

#         fig = px.scatter(
#             density_heatmap_concat_alignments_df.loc[
#                 density_heatmap_concat_alignments_df["MappedGene"].isin(
#                     acceptable_mapped_genes_options
#                 )
#             ],
#             x="Mean%QueryCoverage",
#             y="Mean%QueryIdentity",
#             color="StrQueries",
#             size="Matches",
#             facet_col="GenomicTargetStrand",
#             facet_col_spacing=0.05,
#             facet_row="RelativeTargetStrand",
#             facet_row_spacing=0.05,
#             labels={
#                 "Mean%QueryCoverage": "Mean queries' coverage [%]",
#                 "Mean%QueryIdentity": "Mean queries' identity [%]",
#                 "StrQueries": "Queries",
#                 "GenomicTargetStrand": "Genomic read strand",
#                 "RelativeTargetStrand": "Relative read strand",
#             },
#         )
#         fig.update_xaxes(range=[0, 105], dtick=10)
#         fig.update_yaxes(range=[0, 105], dtick=10)
#         # Reduce opacity to see both histograms
#         fig.update_traces(opacity=0.7)
#         # Overlay both histograms
#         fig.update_layout(
#             # barmode="overlay",
#             template="plotly_white",
#             width=850,
#             height=600,
#             title=f"{gene} samples, reads mapped to {", ".join(acceptable_mapped_genes_options)}",
#             # showlegend=False
#         )

#         fig.show()
# # density_heatmap_concat_alignments_df

# %% [markdown]
# # Choose best combined alignment(s) per read
#

# %% [markdown]
# ## Best

# %%
concat_combined_one_read_aligned_tags_per_strand_df

# %%
# sampled_gene_repeat_and_targets = (
#     concat_combined_one_read_aligned_tags_per_strand_df.groupby(
#         [
#             "Gene",
#             "MappedGene",
#             "GenomicTargetStrand",
#         ],
#         dropna=False,
#     )
#     .sample(n=2, random_state=seed)
#     .loc[:, ["Gene", "Repeat", "TargetName"]]
#     .drop_duplicates()
#     .reset_index(drop=True)
# )

# # sampled_gene_repeat_and_targets

# concat_combined_one_read_aligned_tags_per_strand_df.merge(
#     sampled_gene_repeat_and_targets, how="inner"
# )

# %%
# gene-repeat-read combinations with (one or two) best combined tags' alignment(s) w.r.t strand

best_umi_seq_len = 12
# allowed_umi_seq_len_diff = 2
# min_umi_seq_len = best_umi_seq_len - allowed_umi_seq_len_diff
# max_umi_seq_len = best_umi_seq_len + allowed_umi_seq_len_diff
min_umi_seq_len = 10
max_umi_seq_len = 17

# min_mean_prct_query_cov_and_identity = 95
# min_mean_prct_query_cov_and_identity = 90
# min_mean_prct_query_cov_and_identity = 70
min_mean_prct_query_cov = 85
min_mean_prct_query_identity = 70

min_rel_loc_of_query_1_start = 90

max_gap_openings_per_query = 1
max_total_gaps_per_query = 2
# max_total_gaps_per_query = 4

possible_gap_openings_per_query = range(0, max_gap_openings_per_query + 1)
allowed_query_gap_openings = [
    [x, y] for x, y in product(possible_gap_openings_per_query, repeat=2)
]

possible_total_gaps_per_query = range(0, max_total_gaps_per_query + 1)
allowed_query_total_gaps = [
    [x, y] for x, y in product(possible_total_gaps_per_query, repeat=2)
]

best_concat_combined_one_read_aligned_tags_per_strand_df = (
    concat_combined_one_read_aligned_tags_per_strand_df.loc[
        # (concat_combined_one_read_aligned_tags_per_strand_df["MaxMean%QueryIdentity"])
        # & (concat_combined_one_read_aligned_tags_per_strand_df["MaxMean%QueryCoverage"]) &
        (
            concat_combined_one_read_aligned_tags_per_strand_df[
                "Mean%QueryCoverage"
            ].ge(min_mean_prct_query_cov)
        )
        & (
            concat_combined_one_read_aligned_tags_per_strand_df[
                "Mean%QueryIdentity"
            ].ge(min_mean_prct_query_identity)
        )
        & (
            concat_combined_one_read_aligned_tags_per_strand_df[
                "QueryGapOpenings"
            ].isin(allowed_query_gap_openings)
        )
        & (
            concat_combined_one_read_aligned_tags_per_strand_df["Gaps"].isin(
                allowed_query_total_gaps
            )
        )
        & (
            concat_combined_one_read_aligned_tags_per_strand_df["UMISeqLength"].ge(
                min_umi_seq_len
            )
        )
        & (
            concat_combined_one_read_aligned_tags_per_strand_df["UMISeqLength"].le(
                max_umi_seq_len
            )
        )
        & (
            concat_combined_one_read_aligned_tags_per_strand_df[
                "%RelLocOfQuery1Start"
            ].ge(min_rel_loc_of_query_1_start)
        )
        & (
            concat_combined_one_read_aligned_tags_per_strand_df["Queries"].apply(
                lambda x: x[1] == "PCR"
            )
        )
        & (
            concat_combined_one_read_aligned_tags_per_strand_df[
                "RelativeTargetStrand"
            ].eq("-")
        )
    ]
    .sort_values(
        [
            "Gene",
            "Repeat",
            "TargetName",
            "MappedGene",
            "RelativeTargetStrand",
            "StrQueries",
        ]
    )
    .reset_index(drop=True)
)
best_concat_combined_one_read_aligned_tags_per_strand_df

# %%
# how many unique gene-repeat-read combinations have a best combined tags' alignment w.r.t strand
best_concat_combined_one_read_aligned_tags_per_strand_df.drop_duplicates(
    ["Gene", "Repeat", "TargetName"]
).shape[0]

# %%
best_concat_combined_one_read_aligned_tags_per_strand_df.groupby("Gene")[
    "MappedGene"
].value_counts()

# %%
# best_concat_combined_one_read_aligned_tags_per_strand_df[
#     "%RelLocOfQuery1Start"
# ].describe()

# %%
best_concat_combined_one_read_aligned_tags_per_strand_df[
    "%RelLocOfQuery2End"
].describe()

# %% [markdown]
# ## Bad

# %%
# gene-repeat-read combinations which donesn't have a best combined tags' alignment w.r.t strand
bad_concat_combined_one_read_aligned_tags_per_strand_df = concat_combined_one_read_aligned_tags_per_strand_df.merge(
    best_concat_combined_one_read_aligned_tags_per_strand_df.loc[
        :,
        [
            "Gene",
            "Repeat",
            "Sample",
            "TargetName",
            "GenomicTargetStrand",
            "RelativeTargetStrand",
        ],
    ],
    how="left",
    # on=['Gene', 'Repeat', 'Sample', 'TargetName', 'GenomicTargetStrand',
    #    'RelativeTargetStrand', ],
    indicator="indicator",
)

assert (
    bad_concat_combined_one_read_aligned_tags_per_strand_df.shape[0]
    == concat_combined_one_read_aligned_tags_per_strand_df.shape[0]
)


# in the "bad" df - retain only reads which don't appear in the "good" df
# considering any of their strands or combined tags per strand
bad_concat_combined_one_read_aligned_tags_per_strand_df = (
    bad_concat_combined_one_read_aligned_tags_per_strand_df.loc[
        (
            bad_concat_combined_one_read_aligned_tags_per_strand_df.groupby(
                [
                    "Gene",
                    "Repeat",
                    "TargetName",
                ]
            )["indicator"]
            .transform("nunique")
            .eq(1)
        )
        & (
            bad_concat_combined_one_read_aligned_tags_per_strand_df["indicator"]
            == "left_only"
        )
    ].reset_index(drop=True)
)

assert (
    bad_concat_combined_one_read_aligned_tags_per_strand_df.shape[0] % 4 == 0
), "The number of combined matches should be a multiple of 4 - 2 strands * 2 combined tags per strand"

bad_concat_combined_one_read_aligned_tags_per_strand_df

# %%
# how many unique gene-repeat-read combinations don't have a best combined tags' alignment w.r.t strand
bad_concat_combined_one_read_aligned_tags_per_strand_df.drop_duplicates(
    ["Gene", "Repeat", "TargetName"]
).shape[0]

# %%
assert (
    concat_combined_one_read_aligned_tags_per_strand_df.drop_duplicates(
        ["Gene", "Repeat", "TargetName"]
    ).shape[0]
    == best_concat_combined_one_read_aligned_tags_per_strand_df.drop_duplicates(
        ["Gene", "Repeat", "TargetName"]
    ).shape[0]
    + bad_concat_combined_one_read_aligned_tags_per_strand_df.drop_duplicates(
        ["Gene", "Repeat", "TargetName"]
    ).shape[0]
)

# %% [markdown]
# ### Rescuable 1?

# %%
rescuable_bad_concat_combined_one_read_aligned_tags_per_strand_df = bad_concat_combined_one_read_aligned_tags_per_strand_df.loc[
    (
        bad_concat_combined_one_read_aligned_tags_per_strand_df[
            "Mean%QueryCoverage"
        ].ge(min_mean_prct_query_cov)
    )
    | (
        bad_concat_combined_one_read_aligned_tags_per_strand_df[
            "Mean%QueryIdentity"
        ].ge(min_mean_prct_query_cov_and_identity)
    )
    # | (
    #     bad_concat_combined_one_read_aligned_tags_per_strand_df[
    #         "QueryGapOpenings"
    #     ].isin(allowed_query_gap_openings)
    # )
    # | (
    #     bad_concat_combined_one_read_aligned_tags_per_strand_df["Gaps"].isin(
    #         allowed_query_total_gaps
    #     )
    # )
    | (
        (
            bad_concat_combined_one_read_aligned_tags_per_strand_df["UMISeqLength"].ge(
                min_umi_seq_len
            )
        )
        & (
            bad_concat_combined_one_read_aligned_tags_per_strand_df["UMISeqLength"].le(
                max_umi_seq_len
            )
        )
    ),
    [
        "Gene",
        "Repeat",
        "TargetName",
        "MappedGene",
        "GenomicTargetStrand",
        "RelativeTargetStrand",
        "StrQueries",
        "UMISeqLength",
        "%RelLocOfQuery1Start",
        "%RelLocOfQuery2End",
        "Mean%QueryCoverage",
        "Mean%QueryIdentity",
    ],
].sort_values(
    [
        "Gene",
        "Repeat",
        "TargetName",
    ],
    ignore_index=True,
)
rescuable_bad_concat_combined_one_read_aligned_tags_per_strand_df

# %%
rescuable_bad_concat_combined_one_read_aligned_tags_per_strand_df.loc[
    :, ["Gene", "Repeat", "TargetName"]
].drop_duplicates()

# %%
# checking if the rescuable merged alignments are the best of each group
concat_combined_one_read_aligned_tags_per_strand_df.loc[
    :,
    [
        "Gene",
        "Repeat",
        "TargetName",
        "MappedGene",
        "GenomicTargetStrand",
        "RelativeTargetStrand",
        "StrQueries",
        "UMISeqLength",
        "%RelLocOfQuery1Start",
        "%RelLocOfQuery2End",
        "Mean%QueryCoverage",
        "Mean%QueryIdentity",
    ],
].merge(
    rescuable_bad_concat_combined_one_read_aligned_tags_per_strand_df.loc[
        :, ["Gene", "Repeat", "TargetName"]
    ].drop_duplicates(),
    how="inner",
).sort_values(
    [
        "Gene",
        "Repeat",
        "TargetName",
    ],
    ignore_index=True,
)

# %% [markdown]
# ### Rescuable 2?

# %% [markdown]
# Perhaps long UMI seqs are result of shorter-than-expected alignment of the tags to the read

# %%
rescuable_2_df = bad_concat_combined_one_read_aligned_tags_per_strand_df.loc[
    (
        (
            bad_concat_combined_one_read_aligned_tags_per_strand_df[
                "Mean%QueryCoverage"
            ].ge(min_mean_prct_query_cov)
        )
        | (
            bad_concat_combined_one_read_aligned_tags_per_strand_df[
                "Mean%QueryIdentity"
            ].ge(min_mean_prct_query_identity)
        )
    )
    & (
        bad_concat_combined_one_read_aligned_tags_per_strand_df["UMISeqLength"]
        > max_umi_seq_len
    )
    & (
        bad_concat_combined_one_read_aligned_tags_per_strand_df.apply(
            lambda x: any(
                [
                    alignment_len < query_seq_len
                    for alignment_len, query_seq_len in zip(
                        x["AlignmentLength"], x["QuerySeqLength"]
                    )
                ]
            ),
            axis=1,
        )
    )
].reset_index(drop=True)

rescuable_2_df["MissingAlignmentLength"] = rescuable_2_df.apply(
    lambda x: sum(x["QuerySeqLength"]) - sum(x["AlignmentLength"]), axis=1
)
rescuable_2_df["UMISeqLength-MissingAlignmentLength"] = rescuable_2_df[
    "UMISeqLength"
].sub(rescuable_2_df["MissingAlignmentLength"])

rescuable_2_df

# %%
# how many unique gene-repeat-read combinations have a best combined tags' alignment w.r.t strand
rescuable_2_df.drop_duplicates(["Gene", "Repeat", "TargetName"]).shape[0]

# %%
rescuable_2_df["UMISeqLength-MissingAlignmentLength"].describe()

# %%
rescuable_2_df.loc[
    (rescuable_2_df["UMISeqLength-MissingAlignmentLength"].ge(min_umi_seq_len))
    & (rescuable_2_df["UMISeqLength-MissingAlignmentLength"].le(max_umi_seq_len))
]

# %% [markdown]
# ## Ambiguous best

# %%
# gene-repeat-read combinations that have a best combined tags' alignment w.r.t strand - on both strands
ambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df = (
    best_concat_combined_one_read_aligned_tags_per_strand_df.loc[
        best_concat_combined_one_read_aligned_tags_per_strand_df.duplicated(
            ["Gene", "Repeat", "TargetName"], keep=False
        ),
    ].reset_index(drop=True)
)
ambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df

# %%
# how many unique gene-repeat-read combinations have a best combined tags' alignment w.r.t strand - on both strand
ambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df.drop_duplicates(
    ["Gene", "Repeat", "TargetName"]
).shape[0]

# %%
ambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df.groupby(
    ["Gene", "MappedGene", "RelativeTargetStrand"]
)["StrQueries"].value_counts()


# %% [markdown]
# ## Unambiguous best

# %%
def split_umi_seq_to_unique_sub_seqs(umi_seq, min_umi_seq_len):
    return {
        umi_seq[x : x + min_umi_seq_len]
        for x in range(0, len(umi_seq) - min_umi_seq_len + 1)
    }


# %%
# gene-repeat-read combinations that have a best combined tags' alignment w.r.t strand - on exactly one strand
unambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df = (
    best_concat_combined_one_read_aligned_tags_per_strand_df.loc[
        ~best_concat_combined_one_read_aligned_tags_per_strand_df.duplicated(
            ["Gene", "Repeat", "TargetName"], keep=False
        ),
    ].reset_index(drop=True)
)

# this df should contain no duplicated combined matches per read
assert (
    unambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df.shape
    == unambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df.drop_duplicates(
        ["Gene", "Repeat", "TargetName"]
    ).shape
)

unambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df.insert(
    unambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df.columns.get_loc(
        "UMISeq"
    )
    + 1,
    "UMIUniqueSubSeqs",
    unambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df[
        "UMISeq"
    ].apply(lambda x: split_umi_seq_to_unique_sub_seqs(x, min_umi_seq_len)),
)

assert (
    unambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df[
        "QueriesCoordsOverlap"
    ].value_counts()[False]
    == unambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df.shape[0]
)
del unambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df[
    "QueriesCoordsOverlap"
]

unambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df

# %%
unambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df.groupby(
    ["RelativeTargetStrand", "StrQueries"]
).size()

# %%
unambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df.groupby(
    ["Gene", "MappedGene", "RelativeTargetStrand", "StrQueries"]
).size()

# %%
unambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df.groupby(
    ["Gene", "MappedGene", "GenomicTargetStrand", "RelativeTargetStrand", "StrQueries"]
).size()

# %%
assert (
    best_concat_combined_one_read_aligned_tags_per_strand_df.drop_duplicates(
        ["Gene", "Repeat", "TargetName"]
    ).shape[0]
    == ambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df.drop_duplicates(
        ["Gene", "Repeat", "TargetName"]
    ).shape[
        0
    ]
    + unambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df.shape[0]
)

# %%
# % of gene-repeat-read combinations that have a best combined tags' alignment w.r.t strand - on exactly one strand
100 * unambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df.shape[
    0
] / concat_combined_one_read_aligned_tags_per_strand_df.drop_duplicates(
    ["Gene", "Repeat", "TargetName"]
).shape[
    0
]

# %%
fig = px.histogram(
    unambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df,
    # df,
    # x="UMILength",
    x="UMISeqLength",
    # color="Gene",
    facet_col="Gene",
    # facet_row="UMISeqLength>20",
    color="Repeat",
    # facet_row="Sample",
    histnorm="percent",
    # cumulative=True,
    labels={"UMISeqLength": "UMI length"},
)

fig.update_xaxes(dtick=1)
fig.update_yaxes(dtick=10)

fig.update_traces(
    # # Reduce opacity to see both histograms
    opacity=0.7,
    cumulative_direction="decreasing",
    # selector=dict(type='histogram')
)

# Overlay both histograms
fig.update_layout(
    barmode="overlay",
    template="plotly_white",
    width=600,
    # height=600,
    height=350,
    # showlegend=False
)

fig.show()

# %% [markdown]
# ## Relevant unambiguous best

# %% [markdown]
# Considering only reads mapped to expected gene.

# %%
relevant_unambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df = (
    unambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df.loc[
        (
            unambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df["Gene"]
            == unambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df[
                "MappedGene"
            ]
        )
    ].reset_index(drop=True)
)
relevant_unambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df

# %%
relevant_unambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df.groupby(
    ["Gene", "GenomicTargetStrand", "RelativeTargetStrand", "StrQueries"]
).size()

# %% [markdown]
# # Find unique reads by UMI seqs
#

# %% [markdown] jp-MarkdownHeadingCollapsed=true
# Now, let's see how many unique UMI sequences can I get by:
#
# 1. For each read, split its UMI seq to a list of sub-seqs of length 10 (1-10, 2-11, 3-12, 4-13, 5-14).
# 2. Than, each read is considered unique if none of its sub-UMI seqs is present in any other read.
#

# %% [markdown]
# ## Test sub-UMI seqs comparisons

# %%
gene = "ADAR1"
repeat = "3"
n_rows = 5000

df = (
    unambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df.loc[
        (
            unambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df["Gene"]
            == gene
        )
        & (
            unambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df[
                "Repeat"
            ]
            == repeat
        ),
        ["TargetName", "UMISeq", "UMIUniqueSubSeqs"],
    ]
    .sample(n_rows, random_state=seed)
    .reset_index(drop=True)
)

df


# %%
def old_method(df):

    df = df.explode("UMIUniqueSubSeqs").rename(
        columns={"UMIUniqueSubSeqs": "UMIUniqueSubSeq"}
    )

    df["ReadswithIndistinguishableUMISubSeqs"] = df.apply(
        lambda x: df.loc[
            df["UMIUniqueSubSeq"].eq(x["UMIUniqueSubSeq"]),
            "TargetName",
        ]
        .drop_duplicates()
        .tolist(),
        axis=1,
    )
    df["OtherReadswithIndistinguishableUMISubSeqs"] = df.apply(
        lambda x: [
            read
            for read in x["ReadswithIndistinguishableUMISubSeqs"]
            if read != x["TargetName"]
        ],
        axis=1,
    )
    df["NumOfOtherReadswithIndistinguishableUMISubSeqs"] = df[
        "OtherReadswithIndistinguishableUMISubSeqs"
    ].apply(len)

    df = df.drop_duplicates("TargetName")

    df = df.drop(columns=["UMIUniqueSubSeq"])

    return df


# %%
def at_least_one_shared_subset(umi_1_unique_subseqs: set, umi_2_unique_subseqs: set):

    # return umi_1_unique_subseqs & umi_2_unique_subseqs != set()

    for umi_1_unique_subseq in umi_1_unique_subseqs:
        if umi_1_unique_subseq in umi_2_unique_subseqs:
            return True
    return False


# %%
def new_method(df):
    df = df.copy()
    df["ReadswithIndistinguishableUMISubSeqs"] = df.apply(
        lambda x: df.loc[
            df["UMIUniqueSubSeqs"].apply(
                lambda y: at_least_one_shared_subset(x["UMIUniqueSubSeqs"], y)
            ),
            "TargetName",
        ].tolist(),
        axis=1,
    )
    df["OtherReadswithIndistinguishableUMISubSeqs"] = df.apply(
        lambda x: [
            read
            for read in x["ReadswithIndistinguishableUMISubSeqs"]
            if read != x["TargetName"]
        ],
        axis=1,
    )
    df["NumOfOtherReadswithIndistinguishableUMISubSeqs"] = df[
        "OtherReadswithIndistinguishableUMISubSeqs"
    ].apply(len)

    # df = df.drop(columns=["UMIUniqueSubSeqs"])

    return df


# %%
df_created_by_old_method = old_method(df)
df_created_by_old_method

# %%
df_created_by_new_method = new_method(df)
df_created_by_new_method

# %%
old_with_new_df = df_created_by_new_method.loc[
    :,
    [
        "TargetName",
        "UMISeq",
        "UMIUniqueSubSeqs",
        "NumOfOtherReadswithIndistinguishableUMISubSeqs",
        "OtherReadswithIndistinguishableUMISubSeqs",
    ],
].merge(
    df_created_by_old_method.loc[
        :,
        [
            "TargetName",
            "UMISeq",
            "NumOfOtherReadswithIndistinguishableUMISubSeqs",
            "OtherReadswithIndistinguishableUMISubSeqs",
        ],
    ],
    how="outer",
    on=["TargetName", "UMISeq"],
    indicator="indicator",
    suffixes=("_new", "_old"),
)

assert old_with_new_df["indicator"].value_counts()["both"] == old_with_new_df.shape[0]
del old_with_new_df["indicator"]

old_with_new_df

# %%
disagreeing_old_with_new_df = (
    old_with_new_df.loc[
        old_with_new_df["NumOfOtherReadswithIndistinguishableUMISubSeqs_new"]
        != old_with_new_df["NumOfOtherReadswithIndistinguishableUMISubSeqs_old"],
    ]
    .reset_index(drop=True)
    .drop(columns=["UMISeq"])
)

disagreeing_old_with_new_df

# %%
old_with_new_df.loc[
    old_with_new_df["TargetName"] == "m64296e_241222_071206/100008461/ccs",
    "UMIUniqueSubSeqs",
].values[0]

# %%
umi_1_unique_subseqs = old_with_new_df.loc[
    old_with_new_df["TargetName"] == "m64296e_241222_071206/100008461/ccs",
    "UMIUniqueSubSeqs",
].values[0]
umi_2_unique_subseqs = old_with_new_df.loc[
    old_with_new_df["TargetName"] == "m64296e_241222_071206/25363407/ccs",
    "UMIUniqueSubSeqs",
].values[0]

umi_1_unique_subseqs, umi_2_unique_subseqs

# %%
at_least_one_shared_subset(umi_1_unique_subseqs, umi_2_unique_subseqs)

# %%
df2 = (
    df.loc[
        df["TargetName"].isin(
            [
                "m64296e_241222_071206/100008461/ccs",
                "m64296e_241222_071206/25363407/ccs",
            ]
        )
    ]
    .explode("UMIUniqueSubSeqs")
    .rename(columns={"UMIUniqueSubSeqs": "UMIUniqueSubSeq"})
)

df2["ReadswithIndistinguishableUMISubSeqs"] = df2.apply(
    lambda x: df2.loc[
        df2["UMIUniqueSubSeq"].eq(x["UMIUniqueSubSeq"]),
        "TargetName",
    ]
    .drop_duplicates()
    .tolist(),
    axis=1,
)

df2["OtherReadswithIndistinguishableUMISubSeqs"] = df2.apply(
    lambda x: [
        read
        for read in x["ReadswithIndistinguishableUMISubSeqs"]
        if read != x["TargetName"]
    ],
    axis=1,
)
df2["NumOfOtherReadswithIndistinguishableUMISubSeqs"] = df2[
    "OtherReadswithIndistinguishableUMISubSeqs"
].apply(len)

# df2 = df2.drop_duplicates("TargetName")

# df2 = df2.drop(columns=["UMIUniqueSubSeq"])

df2

# %% [markdown]
# ## The real deal

# %%
min_umi_seq_len, max_umi_seq_len

# %%
relevant_unambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df

# %%
relevant_unambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df[
    "UMISeq"
].nunique()

# %%
relevant_unambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df[
    "UMIUniqueSubSeqs"
]

# %%
# dfs = []

# for gene in ["ADAR1", "IQEC"]:

#     for repeat in list("123"):

#         ic(gene, repeat)

#         gene_and_repeat_df = unambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df.loc[
#             (
#                 unambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df[
#                     "Gene"
#                 ]
#                 == gene
#             )
#             & (
#                 unambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df[
#                     "Repeat"
#                 ]
#                 == repeat
#             ),
#             # ["Gene", "Repeat", "TargetName", "UMIUniqueSubSeqs"],
#         ]

#         df = gene_and_repeat_df.loc[:, ["TargetName", "UMIUniqueSubSeqs"]]
#         df = df.explode("UMIUniqueSubSeqs").rename(
#             columns={"UMIUniqueSubSeqs": "UMIUniqueSubSeq"}
#         )

#         df["ReadswithIndistinguishableUMISubSeqs"] = df.apply(
#             lambda x: df.loc[
#                 df["UMIUniqueSubSeq"].eq(x["UMIUniqueSubSeq"]),
#                 "TargetName",
#             ]
#             .drop_duplicates()
#             .tolist(),
#             axis=1,
#         )
#         df["OtherReadswithIndistinguishableUMISubSeqs"] = df.apply(
#             lambda x: [
#                 read
#                 for read in x["ReadswithIndistinguishableUMISubSeqs"]
#                 if read != x["TargetName"]
#             ],
#             axis=1,
#         )
#         df["NumOfOtherReadswithIndistinguishableUMISubSeqs"] = df[
#             "OtherReadswithIndistinguishableUMISubSeqs"
#         ].apply(len)
#         # df = df.drop_duplicates("TargetName")

#         gene_and_repeat_df = gene_and_repeat_df.merge(
#             df.drop_duplicates("TargetName").loc[
#                 :,
#                 [
#                     "TargetName",
#                     "OtherReadswithIndistinguishableUMISubSeqs",
#                     "NumOfOtherReadswithIndistinguishableUMISubSeqs",
#                 ],
#             ],
#             on="TargetName",
#             how="outer",
#             indicator="indicator",
#         )

#         dfs.append(gene_and_repeat_df)

# new_unambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df = pd.concat(
#     dfs, ignore_index=True
# )
# new_unambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df

# %%
dfs = []

for gene in ["ADAR1", "IQEC"]:

    for repeat in list("123"):

        ic(gene, repeat)

        gene_and_repeat_df = relevant_unambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df.loc[
            (
                relevant_unambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df[
                    "Gene"
                ]
                == gene
            )
            & (
                relevant_unambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df[
                    "Repeat"
                ]
                == repeat
            ),
            # ["Gene", "Repeat", "TargetName", "UMIUniqueSubSeqs"],
        ]

        df = gene_and_repeat_df.loc[:, ["TargetName", "UMIUniqueSubSeqs"]]

        # df = df.explode("UMIUniqueSubSeqs").rename(
        #     columns={"UMIUniqueSubSeqs": "UMIUniqueSubSeq"}
        # )
        # df["ReadswithIndistinguishableUMISubSeqs"] = df.apply(
        #     lambda x: df.loc[
        #         df["UMIUniqueSubSeq"].eq(x["UMIUniqueSubSeq"]),
        #         "TargetName",
        #     ]
        #     .drop_duplicates()
        #     .tolist(),
        #     axis=1,
        # )
        # df["OtherReadswithIndistinguishableUMISubSeqs"] = df.apply(
        #     lambda x: [
        #         read
        #         for read in x["ReadswithIndistinguishableUMISubSeqs"]
        #         if read != x["TargetName"]
        #     ],
        #     axis=1,
        # )
        # df["NumOfOtherReadswithIndistinguishableUMISubSeqs"] = df[
        #     "OtherReadswithIndistinguishableUMISubSeqs"
        # ].apply(len)
        # # df = df.drop_duplicates("TargetName")

        # gene_and_repeat_df = gene_and_repeat_df.merge(
        #     df.drop_duplicates("TargetName").loc[
        #         :,
        #         [
        #             "TargetName",
        #             "OtherReadswithIndistinguishableUMISubSeqs",
        #             "NumOfOtherReadswithIndistinguishableUMISubSeqs",
        #         ],
        #     ],
        #     on="TargetName",
        #     how="outer",
        #     indicator="indicator",
        # )

        df["ReadswithIndistinguishableUMISubSeqs"] = df.apply(
            lambda x: df.loc[
                df["UMIUniqueSubSeqs"].apply(
                    lambda y: at_least_one_shared_subset(x["UMIUniqueSubSeqs"], y)
                ),
                "TargetName",
            ].tolist(),
            axis=1,
        )
        df["OtherReadswithIndistinguishableUMISubSeqs"] = df.apply(
            lambda x: [
                read
                for read in x["ReadswithIndistinguishableUMISubSeqs"]
                if read != x["TargetName"]
            ],
            axis=1,
        )
        df["NumOfOtherReadswithIndistinguishableUMISubSeqs"] = df[
            "OtherReadswithIndistinguishableUMISubSeqs"
        ].apply(len)

        gene_and_repeat_df = gene_and_repeat_df.merge(
            df.loc[
                :,
                [
                    "TargetName",
                    "OtherReadswithIndistinguishableUMISubSeqs",
                    "NumOfOtherReadswithIndistinguishableUMISubSeqs",
                ],
            ],
            on="TargetName",
            how="outer",
            indicator="indicator",
        )

        assert (
            gene_and_repeat_df["indicator"].value_counts()["both"]
            == gene_and_repeat_df.shape[0]
        )
        del gene_and_repeat_df["indicator"]

        dfs.append(gene_and_repeat_df)

new_unambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df = pd.concat(
    dfs, ignore_index=True
)
new_unambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df

# %%
len(dfs)

# %%
new_unambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df

# %%
new_unambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df.loc[
    new_unambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df[
        "NumOfOtherReadswithIndistinguishableUMISubSeqs"
    ].eq(0)
].groupby(["Gene", "Repeat"]).size()

# %%
new_unambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df[
    "NumOfOtherReadswithIndistinguishableUMISubSeqs"
].eq(0).sum()

# %%
# % of reads with unique UMI relative to all reads
100 * new_unambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df[
    "NumOfOtherReadswithIndistinguishableUMISubSeqs"
].eq(0).sum() / concat_bams_df.shape[0]

# %%
# df.loc[df["StrQueries"] == "['PCR', 'ADAR1']"]

# %%
fig = px.histogram(
    df.loc[df["StrQueries"] == "['PCR', 'ADAR1']"],
    # df,
    # x="UMILength",
    x="Query1Start",
    # color="Gene",
    facet_col="Gene",
    # facet_row="UMISeqLength>20",
    color="Repeat",
    # facet_row="Sample",
    # histnorm="percent",
    # cumulative=True,
    # labels={"UMISeqLength": "UMI length"},
)

fig.update_xaxes(dtick=1)
fig.update_yaxes(dtick=10)

fig.update_traces(
    # # Reduce opacity to see both histograms
    opacity=0.7,
    # cumulative_direction="decreasing",
    # selector=dict(type='histogram')
)

# Overlay both histograms
fig.update_layout(
    barmode="overlay",
    template="plotly_white",
    width=600,
    # height=600,
    height=350,
    # showlegend=False
)

fig.show()

# %%

# %%
(
    new_unambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df.assign(
        StrQueries=new_unambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df[
            "Queries"
        ].apply(
            lambda x: str(x)
        )
    )
    .groupby(["GenomicTargetStrand", "RelativeTargetStrand", "StrQueries"])[
        ["Query1Start", "Query2End"]
    ]
    .value_counts()
)

# %%

# %%
new_unambiguous_best_concat_combined_one_read_aligned_tags_per_strand_df.groupby(
    ["Gene", "Repeat"]
)["NumOfOtherReadswithIndistinguishableUMISubSeqs"].value_counts()

# %% [markdown]
# # Choose best gene-specific barcode alignment

# %% [markdown]
# ## Tests before construction

# %%
concat_alignments_df

# %%
test_concat_alignments_df = pd.concat(
    [
        concat_alignments_df.loc[
            (concat_alignments_df["Barcode"] != "PCR")
            & (concat_alignments_df["Gene"].eq(concat_alignments_df["MappedGene"]))
            & (concat_alignments_df["Gene"].eq(concat_alignments_df["Barcode"]))
            & (concat_alignments_df["BTRGaps"].gt(0))
        ].sample(n=3, random_state=seed),
        concat_alignments_df.loc[
            (concat_alignments_df["Barcode"] != "PCR")
            & (concat_alignments_df["Gene"].eq(concat_alignments_df["MappedGene"]))
            & (concat_alignments_df["Gene"].eq(concat_alignments_df["Barcode"]))
            & (concat_alignments_df["BTRGaps"].eq(0))
        ].sample(n=3, random_state=seed),
    ],
    ignore_index=True,
)

test_concat_alignments_df.insert(
    test_concat_alignments_df.columns.get_loc("BTRReadCoords") + 1,
    "BTRReadStart",
    test_concat_alignments_df["BTRReadCoords"].apply(lambda x: x[0][0]),
)
test_concat_alignments_df.insert(
    test_concat_alignments_df.columns.get_loc("BTRReadCoords") + 2,
    "BTRReadEnd",
    test_concat_alignments_df["BTRReadCoords"].apply(lambda x: x[-1][-1]),
)
test_concat_alignments_df["BTGGeneStart"] = test_concat_alignments_df.apply(
    lambda x: x["RTGGeneStart"] + (x["BTRReadStart"] - x["RTGReadStart"]), axis=1
)

test_concat_alignments_df

# %%
test_concat_alignments_df = pd.concat(
    [
        concat_alignments_df.loc[
            (concat_alignments_df["Barcode"] != "PCR")
            & (concat_alignments_df["Gene"].eq(concat_alignments_df["MappedGene"]))
            & (concat_alignments_df["Gene"].eq(concat_alignments_df["Barcode"]))
            & (concat_alignments_df["BTRGaps"].gt(0))
        ].sample(n=3, random_state=seed),
        concat_alignments_df.loc[
            (concat_alignments_df["Barcode"] != "PCR")
            & (concat_alignments_df["Gene"].eq(concat_alignments_df["MappedGene"]))
            & (concat_alignments_df["Gene"].eq(concat_alignments_df["Barcode"]))
            & (concat_alignments_df["BTRGaps"].eq(0))
        ].sample(n=3, random_state=seed),
    ],
    ignore_index=True,
)

test_concat_alignments_df.insert(
    test_concat_alignments_df.columns.get_loc("BTRReadCoords") + 1,
    "BTRReadStart",
    test_concat_alignments_df["BTRReadCoords"].apply(lambda x: x[0][0]),
)
test_concat_alignments_df.insert(
    test_concat_alignments_df.columns.get_loc("BTRReadCoords") + 2,
    "BTRReadEnd",
    test_concat_alignments_df["BTRReadCoords"].apply(lambda x: x[-1][-1]),
)

test_concat_alignments_df


# %%
def get_genomic_coord_for_read_coord(aligned_pairs, required_read_pos):
    for read_pos, ref_pos, cigar_op in aligned_pairs:
        if read_pos == required_read_pos:
            return ref_pos
    return None


# %%
input_df = test_concat_alignments_df


annotated_dfs = []

all_required_reads = input_df["Read"].values

for bam_file in mapped_bam_files:

    with pysam.AlignmentFile(
        bam_file,
        "rb",
        threads=10,
    ) as samfile:
        reads = [read for read in samfile if read.query_name in all_required_reads]
        reads_names = [
            read.query_name for read in reads
        ]  # names of reads found in this specific bam file
        if len(reads_names) == 0:
            continue
        aligned_pairs = [
            read.get_aligned_pairs(matches_only=False, with_cigar=True)
            for read in reads
        ]
        annotated_df = input_df.loc[input_df["Read"].isin(reads_names)].copy()
        aligned_pairs_series = pd.Series(aligned_pairs, index=reads_names)
        annotated_df["BTGGeneStart"] = annotated_df.apply(
            lambda x: get_genomic_coord_for_read_coord(
                aligned_pairs_series.loc[x["Read"]], x["BTRReadStart"]
            ),
            axis=1,
        )
        annotated_df["BTGGeneEnd"] = annotated_df.apply(
            lambda x: get_genomic_coord_for_read_coord(
                aligned_pairs_series.loc[x["Read"]], x["BTRReadEnd"]
            ),
            axis=1,
        )
        annotated_dfs.append(annotated_df)

concat_annotated_df = pd.concat(annotated_dfs)
concat_annotated_df

# %%
annotated_df

# %%
aligned_pairs_series

# %%
df = pd.DataFrame(
    aligned_pairs_series.loc[annotated_df["Read"].values[0]],
    columns=["ReadPos", "RefPos", "CigarOp"],
).set_index("ReadPos")
df.loc[3371 - 3 :]

# %%
df = pd.DataFrame(
    aligned_pairs_series.loc[annotated_df["Read"].values[0]],
    columns=["ReadPos", "RefPos", "CigarOp"],
).set_index("ReadPos")
df.loc[3371 - 3 :]["CigarOp"].value_counts()


# %%
def find_btg_gene_coords(input_df):

    annotated_dfs = []

    all_required_reads = input_df["Read"].values

    for bam_file in mapped_bam_files:

        with pysam.AlignmentFile(
            bam_file,
            "rb",
            threads=10,
        ) as samfile:
            reads = [read for read in samfile if read.query_name in all_required_reads]
            reads_names = [
                read.query_name for read in reads
            ]  # names of reads found in this specific bam file
            if len(reads_names) == 0:
                continue
            aligned_pairs = [
                read.get_aligned_pairs(matches_only=False, with_cigar=True)
                for read in reads
            ]
            annotated_df = input_df.loc[input_df["Read"].isin(reads_names)].copy()
            aligned_pairs_series = pd.Series(aligned_pairs, index=reads_names)
            annotated_df["BTGGeneStart"] = annotated_df.apply(
                lambda x: get_genomic_coord_for_read_coord(
                    aligned_pairs_series.loc[x["Read"]], x["BTRReadStart"]
                ),
                axis=1,
            )
            annotated_df["BTGGeneEnd"] = annotated_df.apply(
                lambda x: get_genomic_coord_for_read_coord(
                    aligned_pairs_series.loc[x["Read"]], x["BTRReadEnd"]
                ),
                axis=1,
            )
            annotated_dfs.append(annotated_df)

    concat_annotated_df = pd.concat(annotated_dfs)
    return concat_annotated_df


# %%
annotated_test_concat_alignments_df = find_btg_gene_coords(test_concat_alignments_df)
annotated_test_concat_alignments_df

# %%

# %%
test_aligned_pairs_dfs = []

for bam_file in mapped_bam_files:

    with pysam.AlignmentFile(
        bam_file,
        "rb",
        threads=10,
    ) as samfile:
        reads = [
            read
            for read in samfile
            if read.query_name in test_concat_alignments_df["Read"].values
        ]
        reads_names = [read.query_name for read in reads]
        aligned_pairs = [
            read.get_aligned_pairs(matches_only=False, with_cigar=True)
            for read in reads
        ]
        df = pd.Series(aligned_pairs, index=reads_names)
        test_aligned_pairs_dfs.append(df)

concat_test_aligned_pairs_df = pd.concat(test_aligned_pairs_dfs)
concat_test_aligned_pairs_df


# %%
def get_genomic_coord_for_read_coord(aligned_pairs, required_read_pos):
    for read_pos, ref_pos, cigar_op in aligned_pairs:
        if read_pos == required_read_pos:
            return ref_pos
    return None


# %%
read = "m64296e_241222_071206/1376592/ccs"

test_concat_alignments_df.loc[test_concat_alignments_df["Read"].eq(read)]

# %%
# concat_test_aligned_pairs_df.loc[read]

# %%
test_concat_alignments_df["BTGGeneStart2"] = test_concat_alignments_df.apply(
    lambda x: get_genomic_coord_for_read_coord(
        concat_test_aligned_pairs_df.loc[x["Read"]], x["BTRReadStart"]
    ),
    axis=1,
)

test_concat_alignments_df

# %%
df = pd.DataFrame(
    concat_test_aligned_pairs_df.loc["m64296e_241222_071206/32115244/ccs"],
    columns=["ReadPos", "RefPos", "CigarOp"],
).set_index("ReadPos")
df

# %%
df.loc[571 - 2 : 571 + 2]

# %%
1238 + 571

# %%
df.loc[df["CigarOp"] == 2]

# %%
df.loc[:571, "CigarOp"].value_counts()

# %%
df["CigarOp"].value_counts()

# %%

# %%

# %%

# %%
single_barcode_min_umi_seq_length = 12
single_barcode_max_umi_seq_length = 14

# %%

# %%
concat_alignments_df.loc[
    (concat_alignments_df["Barcode"] != "PCR")
    & (concat_alignments_df["Gene"].eq(concat_alignments_df["MappedGene"]))
    & (concat_alignments_df["Gene"].eq(concat_alignments_df["Barcode"]))
    # & (concat_alignments_df["Read"].eq(read))
]

# %%
concat_alignments_df.loc[
    (concat_alignments_df["Barcode"] != "PCR")
    & (concat_alignments_df["Gene"].eq(concat_alignments_df["MappedGene"]))
    & (concat_alignments_df["Gene"].eq(concat_alignments_df["Barcode"]))
    & (concat_alignments_df["BTRGaps"].gt(0))
]

# %%
read = "m64296e_241222_071206/100007945/ccs"


concat_alignments_df.loc[
    (concat_alignments_df["Barcode"] != "PCR")
    & (concat_alignments_df["Gene"].eq(concat_alignments_df["MappedGene"]))
    & (concat_alignments_df["Gene"].eq(concat_alignments_df["Barcode"]))
    & (concat_alignments_df["Read"].eq(read))
]

# %%
read = "m64296e_241222_071206/100008153/ccs"

concat_alignments_df.loc[
    (concat_alignments_df["Barcode"] != "PCR")
    & (concat_alignments_df["Gene"].eq(concat_alignments_df["MappedGene"]))
    & (concat_alignments_df["Gene"].eq(concat_alignments_df["Barcode"]))
    & (concat_alignments_df["Read"].eq(read))
]

# %%
read = "m64296e_241222_071206/100073695/ccs"

concat_alignments_df.loc[
    (concat_alignments_df["Barcode"] != "PCR")
    & (concat_alignments_df["Gene"].eq(concat_alignments_df["MappedGene"]))
    & (concat_alignments_df["Gene"].eq(concat_alignments_df["Barcode"]))
    & (concat_alignments_df["Read"].eq(read))
]

# %%

# %%

# %%

# %% [markdown]
# ## Construction

# %%
gene_specific_concat_alignments_df = concat_alignments_df.loc[
    (concat_alignments_df["Barcode"] != "PCR")
    & (concat_alignments_df["Gene"].eq(concat_alignments_df["MappedGene"]))
    & (concat_alignments_df["Gene"].eq(concat_alignments_df["Barcode"]))
    # & (concat_alignments_df["BTRStrand"].eq("-"))
].reset_index(drop=True)

gene_specific_concat_alignments_df.insert(
    gene_specific_concat_alignments_df.columns.get_loc("BTRReadCoords") + 1,
    "BTRReadStart",
    gene_specific_concat_alignments_df["BTRReadCoords"].apply(lambda x: x[0][0]),
)
gene_specific_concat_alignments_df.insert(
    gene_specific_concat_alignments_df.columns.get_loc("BTRReadCoords") + 2,
    "BTRReadEnd",
    gene_specific_concat_alignments_df["BTRReadCoords"].apply(lambda x: x[-1][-1]),
)

# # gene_specific_concat_alignments_df["BTGGeneStart"] = (
# #     gene_specific_concat_alignments_df.apply(
# #         lambda x: x["RTGGeneStart"] + (x["BTRReadStart"] - x["RTGReadStart"]), axis=1
# #     )
# # )

# # get BTGeneStart and BTGGeneEnd for each read
# gene_specific_concat_alignments_df = find_btg_gene_coords(
#     gene_specific_concat_alignments_df
# )

gene_specific_concat_alignments_df

# %%
Out[291]

# %%
gene_specific_concat_alignments_df["BTGGeneStart"].describe()

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

gene_specific_concat_alignments_df

# %% [markdown]
# ## Tests

# %%
gene_specific_concat_alignments_df.groupby("Gene")["BTRBarcodeInTail"].value_counts(
    normalize=True
).mul(100).round(2)

# %%
gene_specific_concat_alignments_df.groupby("Gene")["BTGGeneStart-RTGGeneEnd"].describe()

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
fig = px.histogram(
    gene_specific_concat_alignments_df.loc[
        (gene_specific_concat_alignments_df["BTGGeneStart"].ge(0))
        # & (~gene_specific_concat_alignments_df["BTRBarcodeInTail"])
    ],
    # df,
    # x="UMILength",
    x="RTGGeneEnd",
    # color="Gene",
    facet_col="Gene",
    # facet_row="BTRBarcodeInTail",
    # color="BTRBarcodeInTail",
    # color="Repeat",
    # facet_row="Sample",
    # histnorm="percent",
    # cumulative=True,
    # labels={"UMISeqLength": "UMI length"},
    # log_x=True,
    log_y=True,
    # opacity=0.7,
)
# fig.update_xaxes(range=[0, None], tick0=0, dtick=100)
fig.update_layout(
    width=800,
    #   height=350
    height=350,
    barmode="overlay",
    # barmode="stack",
)
fig.show()

# %%
fig = px.histogram(
    gene_specific_concat_alignments_df.loc[
        (gene_specific_concat_alignments_df["BTGGeneStart"].ge(0))
        # & (~gene_specific_concat_alignments_df["BTRBarcodeInTail"])
    ],
    # df,
    # x="UMILength",
    x="RTGGeneEnd",
    # color="Gene",
    facet_col="Gene",
    # facet_row="BTRBarcodeInTail",
    # color="BTRBarcodeInTail",
    # color="Repeat",
    # facet_row="Sample",
    histnorm="percent",
    # cumulative=True,
    # labels={"UMISeqLength": "UMI length"},
    # log_x=True,
    log_y=True,
    # opacity=0.7,
)
# fig.update_xaxes(range=[0, None], tick0=0, dtick=100)
fig.update_layout(
    width=800,
    #   height=350
    height=350,
    barmode="overlay",
    # barmode="stack",
)
fig.show()

# %%
fig = px.histogram(
    gene_specific_concat_alignments_df.loc[
        (gene_specific_concat_alignments_df["BTGGeneStart"].ge(0))
        & (~gene_specific_concat_alignments_df["BTRBarcodeInTail"])
    ],
    # df,
    # x="UMILength",
    x="BTGGeneStart-RTGGeneEnd",
    # color="Gene",
    facet_col="Gene",
    # facet_row="BTRBarcodeInTail",
    # color="BTRBarcodeInTail",
    # color="Repeat",
    # facet_row="Sample",
    # histnorm="percent",
    # cumulative=True,
    # labels={"UMISeqLength": "UMI length"},
    # log_x=True,
    log_y=True,
    # opacity=0.7,
)
fig.update_xaxes(range=[0, None], tick0=0, dtick=100)
fig.update_layout(
    width=800,
    #   height=350
    height=350,
    barmode="overlay",
    # barmode="stack",
)
fig.show()

# %%
fig = px.histogram(
    gene_specific_concat_alignments_df.loc[
        (gene_specific_concat_alignments_df["BTGGeneStart"].ge(0))
        # & (~gene_specific_concat_alignments_df["BTRBarcodeInTail"])
    ],
    # df,
    # x="UMILength",
    x="BTGGeneStart-RTGGeneEnd",
    # color="Gene",
    # facet_row="Gene",
    # # facet_row="BTRBarcodeInTail",
    # # color="BTRBarcodeInTail",
    # facet_col="BTRBarcodeInTail",
    facet_col="Gene",
    # facet_row="BTRBarcodeInTail",
    # color="BTRBarcodeInTail",
    facet_row="BTRBarcodeInTail",
    # color="Repeat",
    # facet_row="Sample",
    # histnorm="percent",
    # cumulative=True,
    # labels={"UMISeqLength": "UMI length"},
    # log_x=True,
    log_y=True,
    # opacity=0.7,
)
# fig.update_xaxes(
#     # range=[0, None], tick0=0,
#     dtick=20
# )
fig.update_layout(
    width=800,
    #   height=350
    height=550,
    # barmode="overlay",
    # barmode="stack",
)
fig.show()

# %%
fig = px.histogram(
    gene_specific_concat_alignments_df.loc[
        (gene_specific_concat_alignments_df["BTGGeneStart"].ge(0))
        & (~gene_specific_concat_alignments_df["BTRBarcodeInTail"])
    ],
    # df,
    # x="UMILength",
    x="BTGGeneStart-RTGGeneEnd",
    # color="Gene",
    facet_row="Gene",
    # facet_row="BTRBarcodeInTail",
    # color="BTRBarcodeInTail",
    # color="Repeat",
    # facet_row="Sample",
    histnorm="percent",
    # cumulative=True,
    # labels={"UMISeqLength": "UMI length"},
    # log_x=True,
    # log_y=True,
    # opacity=0.7,
)
# fig.update_xaxes(range=[0, None], tick0=0, dtick=20)
fig.update_layout(
    width=1400,
    #   height=350
    height=550,
    # barmode="overlay",
    # barmode="stack",
)
fig.show()

# %%
fig = px.histogram(
    gene_specific_concat_alignments_df.loc[
        gene_specific_concat_alignments_df["BTGGeneStart"].ge(0)
    ],
    # df,
    # x="UMILength",
    x="BTGGeneStart",
    # color="Gene",
    facet_col="Gene",
    # facet_row="BTRBarcodeInTail",
    color="BTRBarcodeInTail",
    # color="Repeat",
    # facet_row="Sample",
    # histnorm="percent",
    # cumulative=True,
    # labels={"UMISeqLength": "UMI length"},
    log_y=True,
    opacity=0.7,
)
fig.update_xaxes(range=[0, None], tick0=0, dtick=1000)
fig.update_layout(
    width=800,
    #   height=350
    height=450,
    barmode="overlay",
    # barmode="stack",
)
fig.show()

# %%
fig = px.histogram(
    gene_specific_concat_alignments_df.loc[
        gene_specific_concat_alignments_df["BTGGeneStart"].ge(0)
    ],
    # df,
    # x="UMILength",
    x="BTGGeneStart",
    # color="Gene",
    facet_col="Gene",
    # color="BTRBarcodeInTail",
    # facet_row="BTRBarcodeInTail",
    # color="Repeat",
    # facet_row="Sample",
    histnorm="percent",
    # cumulative=True,
    # labels={"UMISeqLength": "UMI length"},
    log_y=True,
    opacity=0.7,
)
fig.update_xaxes(range=[0, None], tick0=0, dtick=1000)
fig.update_layout(
    width=800,
    height=350,
    # height=600,
    # height=450,
    # barmode="overlay",
    # barmode="stack",
)
fig.show()

# %%
gene_specific_concat_alignments_df.groupby("Gene")["BTGGeneStart"].describe()

# %%
gene_specific_concat_alignments_df.loc[
    gene_specific_concat_alignments_df["BTGGeneStart"].lt(0),
    # [
    #     "Sample",
    #     "Gene",
    #     "Repeat",
    #     "Read",
    #     "ReadSeq",
    #     "RTGStrand",
    #     "RTGGeneStart",
    #     "BTRReadStart",
    #     "RTGReadStart",
    #     "BTGGeneStart",
    # ],
]

# %%
Out[279]["RTGGeneStart"].describe()

# %%
Out[279]["BTRReadStart"].describe()

# %%
Out[279]["RTGReadStart"].describe()

# %%
Out[279]["BTRReadStart"] - Out[279]["RTGReadStart"]

# %%
fig = px.histogram(
    gene_specific_concat_alignments_df,
    # df,
    # x="UMILength",
    x="BTGGeneStart",
    # color="Gene",
    facet_col="Gene",
    # facet_row="UMISeqLength>20",
    # color="Repeat",
    # facet_row="Sample",
    # histnorm="percent",
    # cumulative=True,
    # labels={"UMISeqLength": "UMI length"},
)
# fig.update_xaxes(dtick=1)
fig.update_layout(width=800, height=350)
fig.show()

# %%
gene_specific_concat_alignments_df["BTRReadStart"]

# %%
(
    gene_specific_concat_alignments_df["ReadSeqLength"].sub(
        gene_specific_concat_alignments_df["BTRReadStart"]
    )
).describe()

# %%
(
    (
        gene_specific_concat_alignments_df["ReadSeqLength"].sub(
            gene_specific_concat_alignments_df["BTRReadStart"]
        )
    )
    .mul(100)
    .div(gene_specific_concat_alignments_df["ReadSeqLength"])
).describe()

# %%
(
    (
        gene_specific_concat_alignments_df["ReadSeqLength"].sub(
            gene_specific_concat_alignments_df["BTRReadStart"]
        )
    )
    .mul(100)
    .div(gene_specific_concat_alignments_df["ReadSeqLength"])
).add(gene_specific_concat_alignments_df["%BTRRelLocOfReadStart"]).describe()

# %%
gene_specific_concat_alignments_df.loc[
    (
        gene_specific_concat_alignments_df["Gene"].eq(
            gene_specific_concat_alignments_df["MappedGene"]
        )
    )
    & (
        gene_specific_concat_alignments_df["Gene"].eq(
            gene_specific_concat_alignments_df["Barcode"]
        )
    ),
    # & (gene_specific_concat_alignments_df["BTRStrand"].eq("-"))
    "RTGReadTailLength",
].describe().round(2)

# %%
gene_specific_concat_alignments_df.loc[
    (
        gene_specific_concat_alignments_df["Gene"].eq(
            gene_specific_concat_alignments_df["MappedGene"]
        )
    )
    & (
        gene_specific_concat_alignments_df["Gene"].eq(
            gene_specific_concat_alignments_df["Barcode"]
        )
    ),
    # & (gene_specific_concat_alignments_df["BTRStrand"].eq("-"))
    "RTGReadTailLength",
].describe().round(2)

# %%
gene_specific_concat_alignments_df.loc[
    (
        gene_specific_concat_alignments_df["Gene"].eq(
            gene_specific_concat_alignments_df["MappedGene"]
        )
    )
    & (
        gene_specific_concat_alignments_df["Gene"].eq(
            gene_specific_concat_alignments_df["Barcode"]
        )
    )
    & (gene_specific_concat_alignments_df["BTRStrand"].eq("-")),
    # "RTGReadTail",
].groupby("Gene")["RTGReadTailLength"].describe().round(2)

# %%
gene_specific_concat_alignments_df.loc[
    (
        gene_specific_concat_alignments_df["Gene"].eq(
            gene_specific_concat_alignments_df["MappedGene"]
        )
    )
    & (
        gene_specific_concat_alignments_df["Gene"].eq(
            gene_specific_concat_alignments_df["Barcode"]
        )
    )
    & (gene_specific_concat_alignments_df["BTRStrand"].eq("-"))
    & (gene_specific_concat_alignments_df["RTGReadTailLength"].ge(100)),
    # "RTGReadTail",
].groupby("Gene")["RTGReadTailLength"].describe().round(2)

# %%
gene_specific_concat_alignments_df.loc[
    (
        gene_specific_concat_alignments_df["Gene"].eq(
            gene_specific_concat_alignments_df["MappedGene"]
        )
    )
    & (
        gene_specific_concat_alignments_df["Gene"].eq(
            gene_specific_concat_alignments_df["Barcode"]
        )
    )
    & (gene_specific_concat_alignments_df["BTRStrand"].eq("-"))
    & (gene_specific_concat_alignments_df["RTGReadTailLength"].ge(33)),
    # "RTGReadTail",
].groupby("Gene")["RTGReadTailLength"].describe().round(2)["count"].sum()

# %%
142493.0 + 169132.0

# %%
pcr_primer_len = len(primers_dict["PCR"])
for gene in genes:
    primers_len_without_pcr = len(primers_dict[gene]) + 10  # minimal UMI length
    primers_len_with_pcr = primers_len_without_pcr + pcr_primer_len
    ic(gene, primers_len_without_pcr, primers_len_with_pcr)

# %%
gene_specific_concat_alignments_df.loc[
    (
        (gene_specific_concat_alignments_df["Gene"].eq("ADAR1"))
        & (gene_specific_concat_alignments_df["RTGReadTailLength"].ge(56))
    )
    | (
        (gene_specific_concat_alignments_df["Gene"].eq("IQEC"))
        & (gene_specific_concat_alignments_df["RTGReadTailLength"].ge(54))
    ),
    # "RTGReadTail",
]

# %%
gene_specific_concat_alignments_df.loc[
    (
        (gene_specific_concat_alignments_df["Gene"].eq("ADAR1"))
        & (gene_specific_concat_alignments_df["RTGReadTailLength"].ge(56))
    )
    | (
        (gene_specific_concat_alignments_df["Gene"].eq("IQEC"))
        & (gene_specific_concat_alignments_df["RTGReadTailLength"].ge(54))
    ),
    # "RTGReadTail",
].groupby("Gene")["RTGReadTailLength"].describe().round(2)

# %%
fig = px.histogram(
    gene_specific_concat_alignments_df.loc[
        (
            gene_specific_concat_alignments_df["Gene"].eq(
                gene_specific_concat_alignments_df["MappedGene"]
            )
        )
        & (
            gene_specific_concat_alignments_df["Gene"].eq(
                gene_specific_concat_alignments_df["Barcode"]
            )
        )
        & (gene_specific_concat_alignments_df["BTRStrand"].eq("-"))
        & (gene_specific_concat_alignments_df["RTGReadTailLength"].le(100))
    ],
    x="RTGReadTailLength",
    facet_col="Gene",
    # facet_row="BTRStrand",
    # log_x=True,
    log_y=True,
)
fig.update_xaxes(dtick=5)
fig.update_layout(width=800, height=350, template="plotly_white")
fig.show()

# %%
fig = px.histogram(
    gene_specific_concat_alignments_df.loc[
        (
            gene_specific_concat_alignments_df["Gene"].eq(
                gene_specific_concat_alignments_df["MappedGene"]
            )
        )
        & (
            gene_specific_concat_alignments_df["Gene"].eq(
                gene_specific_concat_alignments_df["Barcode"]
            )
        )
        # & (gene_specific_concat_alignments_df["BTRStrand"].eq("-"))
    ],
    x="BTGGeneStart-ORFEnd",
    facet_col="Gene",
    facet_row="BTRStrand",
    log_y=True,
)
fig.update_layout(width=800, height=500, template="plotly_white")
fig.show()

# %%
(
    gene_specific_concat_alignments_df.loc[
        (
            gene_specific_concat_alignments_df["Gene"].eq(
                gene_specific_concat_alignments_df["MappedGene"]
            )
        )
        & (
            gene_specific_concat_alignments_df["Gene"].eq(
                gene_specific_concat_alignments_df["Barcode"]
            )
        )
        # & (gene_specific_concat_alignments_df["BTRStrand"].eq("-"))
    ]
    .groupby(
        [
            "Gene",
            "BTRStrand",
        ]
    )["BTGGeneStart-ORFEnd"]
    .describe()
    .round(2)
)

# %%
gene_specific_concat_alignments_df.loc[
    (
        gene_specific_concat_alignments_df["Gene"].eq(
            gene_specific_concat_alignments_df["MappedGene"]
        )
    )
    & (
        gene_specific_concat_alignments_df["Gene"].eq(
            gene_specific_concat_alignments_df["Barcode"]
        )
    )
    & (gene_specific_concat_alignments_df["%RTGORFCoverageUntilBarcode"] < 0),
    "%RTGORFCoverageUntilBarcode",
].describe()

# %%
genes_orfs_boundries_dict

# %%
fig = px.histogram(
    gene_specific_concat_alignments_df.loc[
        (
            gene_specific_concat_alignments_df["Gene"].eq(
                gene_specific_concat_alignments_df["MappedGene"]
            )
        )
        & (
            gene_specific_concat_alignments_df["Gene"].eq(
                gene_specific_concat_alignments_df["Barcode"]
            )
        )
    ],
    x="%RTGORFCoverageUntilBarcode",
    facet_col="Gene",
    histnorm="percent",
    log_y=True,
    cumulative=True,
)
fig.update_xaxes(dtick=10)
fig.update_layout(width=800, height=400, template="plotly_white")
fig.show()

# %%
gene_specific_concat_alignments_df.loc[
    (
        gene_specific_concat_alignments_df["BTGGeneStart"].ge(
            gene_specific_concat_alignments_df["ORFEnd"]
        )
    )
]

# %%
gene_specific_concat_alignments_df.loc[
    (
        gene_specific_concat_alignments_df["BTGGeneStart"].lt(
            gene_specific_concat_alignments_df["ORFEnd"]
        )
    )
]

# %%

# %%
gene_specific_concat_alignments_df.loc[
    (
        gene_specific_concat_alignments_df["BTGGeneStart"].ge(
            gene_specific_concat_alignments_df["ORFEnd"]
        )
    )
    & (
        gene_specific_concat_alignments_df["BTRReadStart"].ge(
            gene_specific_concat_alignments_df["RTGReadEnd"]
        )
    )
    & (
        gene_specific_concat_alignments_df["RTGGeneEnd"].ge(
            gene_specific_concat_alignments_df["ORFEnd"]
        )
    )
]

# %%
gene_specific_concat_alignments_df.loc[
    (
        gene_specific_concat_alignments_df["BTRReadStart"].ge(
            gene_specific_concat_alignments_df["RTGReadEnd"]
        )
    )
    & (
        gene_specific_concat_alignments_df["RTGGeneEnd"].ge(
            gene_specific_concat_alignments_df["ORFEnd"]
        )
    )
]

# %%
gene_specific_concat_alignments_df.apply(
    lambda x: 
)

# %%
genes_orfs_boundries_dict

# %%
gene_specific_concat_alignments_df.drop_duplicates(["Gene", "Repeat", "Read"]).shape[0]

# %%
{
    "TargetName": "Read",
    "QueryName": "Barcode",
    "TargetSeq": "ReadSeq",
    "GenomicTargetStrand": "RTGStrand",
    "TargetStrand": "BTRStrand",
    "%QueryIdentity": "%BTRBarcodeIdentity",
    "%QueryCoverage": "%BTRBarcodeCoverage",
    "TargetSeqLength": "ReadSeqLength",
    "QuerySeqLength": "BarcodeSeqLength",
    "AlignmentLength": "BTRAlignmentLength",
    "AlignedTargetCoords": "BTRReadCoords",
    "AlignedQueryCoords": "BTRBarcodeCoords",
    "AlignmentObject": "BTRAlignmentObject",
    "NumOfAlignedTargetCoords": "NumOfBTRReadCoords",
    "NumOfAlignedQueryCoords": "NumOfBTRBarcodeCoords",
    "QueryGapOpenings": "BTRBarcodeGapOpenings",
    "Score": "BTRScore",
    "Gaps": "BTRGaps",
    "Identitites": "BTRIdentitites",
    "Mismatches": "BTRMismatches",
    "TargetStart": "BTRReadStart",
    "TargetEnd": "BTRReadEnd",
    "%RelLocOfTargetStart": "%BTRRelLocOfReadStart",
}

# %% [markdown]
# ## Actually choosing

# %%
# gene-repeat-read combinations with (one or two) best combined tags' alignment(s) w.r.t strand


# min_prct_query_cov = 90
min_prct_query_cov = 85
# min_prct_query_cov = 80
# min_prct_query_identity = 90
min_prct_query_identity = 85
# min_prct_query_identity = 80


min_rel_loc_of_target_start = 90

max_gap_openings_per_query = 1
# max_total_gaps_per_query = 2
max_total_gaps_per_query = 4

possible_gap_openings_per_gene_specific_query = list(
    range(0, max_gap_openings_per_query + 1)
)
possible_total_gaps_per_gene_specific_query = list(
    range(0, max_total_gaps_per_query + 1)
)


best_gene_specific_concat_alignments_df = (
    gene_specific_concat_alignments_df.loc[
        (
            gene_specific_concat_alignments_df["%BTRBarcodeCoverage"].ge(
                min_prct_query_cov
            )
        )
        & (
            gene_specific_concat_alignments_df["%BTRBarcodeIdentity"].ge(
                min_prct_query_identity
            )
        )
        & (
            gene_specific_concat_alignments_df["NumOfBTRBarcodeGapOpenings"].isin(
                possible_gap_openings_per_gene_specific_query
            )
        )
        & (
            gene_specific_concat_alignments_df["BTRGaps"].isin(
                possible_total_gaps_per_gene_specific_query
            )
        )
        # & (
        #     gene_specific_concat_alignments_df["%BTRRelLocOfReadStart"].ge(
        #         min_rel_loc_of_target_start
        #     )
        # )
        # & (gene_specific_concat_alignments_df["BTRStrand"].eq("-"))
        # & (
        #     gene_specific_concat_alignments_df["Gene"].eq(
        #         gene_specific_concat_alignments_df["MappedGene"]
        #     )
        # )
        # & (
        #     gene_specific_concat_alignments_df["Gene"].eq(
        #         gene_specific_concat_alignments_df["Barcode"]
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

best_gene_specific_concat_alignments_df["UMISeqLength"] = (
    best_gene_specific_concat_alignments_df["UMISeq"].apply(len)
)
best_gene_specific_concat_alignments_df["UMIUniqueSubSeqs"] = (
    best_gene_specific_concat_alignments_df["UMISeq"].apply(
        lambda x: split_umi_seq_to_unique_sub_seqs(x, 12)
    )
)

assert (
    best_gene_specific_concat_alignments_df.drop_duplicates(
        ["Gene", "Repeat", "Read"]
    ).shape[0]
    == best_gene_specific_concat_alignments_df.shape[0]
)


best_gene_specific_concat_alignments_df

# %%
best_gene_specific_concat_alignments_df["Gene"].value_counts()

# %% [markdown]
# ## Further tests

# %%
fig = px.histogram(
    best_gene_specific_concat_alignments_df,
    # df,
    # x="UMILength",
    x="BTGGeneStart",
    # color="Gene",
    facet_col="Gene",
    facet_row="BTRBarcodeInTail",
    # color="Repeat",
    # facet_row="Sample",
    # histnorm="percent",
    # cumulative=True,
    # labels={"UMISeqLength": "UMI length"},
    log_y=True,
)
# fig.update_xaxes(dtick=1)
fig.update_layout(
    width=800,
    #   height=350
    height=600,
)
fig.show()

# %%
fig = px.histogram(
    best_gene_specific_concat_alignments_df,
    # df,
    # x="UMILength",
    x="BTGGeneStart",
    # color="Gene",
    facet_col="Gene",
    facet_row="BTRBarcodeInTail",
    # color="Repeat",
    # facet_row="Sample",
    histnorm="percent",
    cumulative=True,
    # labels={"UMISeqLength": "UMI length"},
    # log_y=True,
)
# fig.update_xaxes(dtick=1)
fig.update_layout(
    width=800,
    #   height=350
    height=600,
)
fig.show()

# %%
best_gene_specific_concat_alignments_df.loc[
    (best_gene_specific_concat_alignments_df["Gene"].eq("IQEC"))
    & (best_gene_specific_concat_alignments_df["BTGGeneStart"] >= 3770)
    & (best_gene_specific_concat_alignments_df["BTGGeneStart"] <= 3780)
]

# %%
Out[395].iloc[1]["ReadSeq"]

# %%
print(Out[395].iloc[1]["BTRAlignmentObject"])

# %%
best_gene_specific_concat_alignments_df.loc[
    best_gene_specific_concat_alignments_df["UMISeq"].apply(len) <= 12
]

# %%
best_gene_specific_concat_alignments_df["UMISeq"].apply(len).value_counts()

# %%
best_gene_specific_concat_alignments_df["UMIUniqueSubSeqs"].apply(len).value_counts()

# %%
# best_gene_specific_concat_alignments_df.drop_duplicates(
#     ["Gene", "Repeat", "TargetName"]
# ).shape[0]

# %% [markdown]
# # Find unique reads by UMI seqs in gene-specific selected barcodes
#

# %%
gene_specific_dfs = []

for gene in ["ADAR1", "IQEC"]:

    for repeat in list("123"):

        ic(gene, repeat)

        gene_and_repeat_df = best_gene_specific_concat_alignments_df.loc[
            (best_gene_specific_concat_alignments_df["Gene"] == gene)
            & (best_gene_specific_concat_alignments_df["Repeat"] == repeat)
        ]

        df = gene_and_repeat_df.loc[:, ["TargetName", "UMIUniqueSubSeqs"]]

        df["ReadswithIndistinguishableUMISubSeqs"] = df.apply(
            lambda x: df.loc[
                df["UMIUniqueSubSeqs"].apply(
                    lambda y: at_least_one_shared_subset(x["UMIUniqueSubSeqs"], y)
                ),
                "TargetName",
            ].tolist(),
            axis=1,
        )
        df["OtherReadswithIndistinguishableUMISubSeqs"] = df.apply(
            lambda x: [
                read
                for read in x["ReadswithIndistinguishableUMISubSeqs"]
                if read != x["TargetName"]
            ],
            axis=1,
        )
        df["NumOfOtherReadswithIndistinguishableUMISubSeqs"] = df[
            "OtherReadswithIndistinguishableUMISubSeqs"
        ].apply(len)

        gene_and_repeat_df = gene_and_repeat_df.merge(
            df.loc[
                :,
                [
                    "TargetName",
                    "OtherReadswithIndistinguishableUMISubSeqs",
                    "NumOfOtherReadswithIndistinguishableUMISubSeqs",
                ],
            ],
            on="TargetName",
            how="outer",
            indicator="indicator",
        )

        assert (
            gene_and_repeat_df["indicator"].value_counts()["both"]
            == gene_and_repeat_df.shape[0]
        )
        del gene_and_repeat_df["indicator"]

        gene_specific_dfs.append(gene_and_repeat_df)

        # break  # todo remove - this is for test purposes

    # break  # todo remove - this is for test purposes

new_best_gene_specific_concat_alignments_df = pd.concat(
    gene_specific_dfs, ignore_index=True
)
new_best_gene_specific_concat_alignments_df

# %%
new_best_gene_specific_concat_alignments_df

# %%
new_best_gene_specific_concat_alignments_df[
    "NumOfOtherReadswithIndistinguishableUMISubSeqs"
].eq(0).sum()
