# ---
# jupyter:
#   jupytext:
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

# %%
from pathlib import Path

import pandas as pd
import numpy as numpy
import plotly.colors as pc
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# %%
orfs_bed_file = Path("D.pealeii/Annotations/Jan2025/orfs_squ.bed").absolute()

# %%
targeted_chroms = [
    "comp141693_c0_seq1", "comp141882_c0_seq14", "comp134400_c0_seq1_extended", "comp141565_c6_seq3"
]
targeted_gene_names = [
    "GRIA2", "PCLO", "ADAR1", "IQEC1"
]

# %%
bioprojects = ["PRJNA1216919", "PRJNA641326"]
tissues = ["SG", "OG"]

# %%
main_in_dirs = [
    Path(
        "D.pealeii/Salmon", bioproject
    ).absolute()
    for bioproject in bioprojects
]
main_in_dirs

# %%
samples_per_bioproject = [
    [
        x.name 
        for x in main_in_dir.glob("SRR*")
        if x.is_dir()
    ]
    for main_in_dir in main_in_dirs
]
samples_per_bioproject

# %%
quant_files_per_bioproject = [
    [
        Path(main_in_dir, sample, "quant.sf")
        for sample in samples
    ]
    for main_in_dir, samples in zip(main_in_dirs, samples_per_bioproject)
]
quant_files_per_bioproject

# %%
orfs_df = pd.read_table(
    orfs_bed_file,
    names="Chrom Start End Name Score Strand".split()
)
orfs_df

# %%
quant_dfs = []
for bioproject, tissue, quant_files, samples in zip(bioprojects, tissues, quant_files_per_bioproject, samples_per_bioproject):
    for quant_file, sample in zip(quant_files, samples):
        df = pd.read_csv(quant_file, sep="\t")
        df = df.rename(columns={"Name": "Chrom"})
        df.insert(0, "BioProject", bioproject)
        df.insert(1, "Tissue", tissue)
        df.insert(2, "Sample", sample)
        quant_dfs.append(df)
concat_quant_df = pd.concat(quant_dfs, ignore_index=True)
concat_quant_df = concat_quant_df.merge(
    orfs_df.loc[:, ["Chrom", "Name"]],
    on="Chrom",
    how="left"
)
concat_quant_df["Name"] = concat_quant_df["Name"].str.split("_").str[0]
# all cols taken, reordered to my liking
concat_quant_df = concat_quant_df.loc[
    :,
    ['BioProject', 'Tissue', 'Sample', 'Chrom', 'Name', 
     'Length', 'EffectiveLength', 'TPM', 'NumReads', ]
]
concat_quant_df

# %%
concat_targeted_quant_df = concat_quant_df.loc[
    concat_quant_df["Chrom"].isin(targeted_chroms)
]

concat_targeted_quant_df

# %%
['BioProject', 'Tissue', 'Sample', 'Chrom', 'Name', 
     'Length', 'EffectiveLength', 'TPM', 'NumReads', ]

# %%
(
    concat_targeted_quant_df
    .groupby(
        ['Chrom', 'Name'],
        as_index=False
    )
    ["TPM"].describe()
    .round(2)
)

# %%
(
    concat_targeted_quant_df
    .groupby(
        ['Chrom', 'Name', 'BioProject', 'Tissue'],
        as_index=False
    )
    ["TPM"].describe()
    .round(2)
)
