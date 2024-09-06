from multiprocessing import Pool
from pathlib import Path
import sys

import pandas as pd
import numpy as np
from icecream import ic

# print(Path(".").absolute())
# exit()

code_dir = "/private7/projects/Combinatorics/Code"
sys.path.append(str(Path(code_dir).absolute()))

from pileup_with_subparsers import (
    simulate_complete_and_corresponding_partially_unknown_main,
)


pd.set_option("display.max_columns", 20)

rng = np.random.RandomState(1892)

n_reads = 10_000
num_of_mock_genes = 3
min_known_sites_per_mock_gene = 50
max_known_sites_per_mock_gene = 150

processes = 6
threads = 30


# n_reads = 100
# unknown_probability = 0.05  # chance of an editing site's status being unknown
# unknown_probabilities = [
#     0.01,
#     0.03,
#     0.05,
# ]  # chance of an editing site's status being unknown
unknown_probabilities = [
    round(x / 126, 2) for x in [11, 17, 23]
]  # chance of an editing site's status being unknown
snp_noise_level = 0.1
top_x_noisy_positions = 3
assurance_factor = 1.5
denovo_detection = False
out_files_sep = "\t"


known_sites_csv_file = (
    "/private7/projects/Combinatorics/D.pealeii/Annotations/D.pea.EditingSites.csv"
)
known_sites_file = (
    "/private7/projects/Combinatorics/D.pealeii/Annotations/D.pea.EditingSites.bed"
)
orfs_bed_file = "/private7/projects/Combinatorics/D.pealeii/Annotations/orfs_squ.bed"


cols = ["Trinity id", "Editing location (base1)", "SwissProt name", "Editing Level"]
editing_sites_df = (
    pd.read_csv(known_sites_csv_file)
    .filter(cols)
    .rename(
        columns={
            "Trinity id": "Chrom",
            "Editing location (base1)": "Position",
            "SwissProt name": "SwissProtName",
            "Editing Level": "EditingLevel",
        }
    )
    .sort_values(["Chrom", "Position"])
    .reset_index(drop=True)
)
editing_sites_df["Position"] = editing_sites_df["Position"] - 1

orfs_df = (
    pd.read_table(
        orfs_bed_file, names=["Chrom", "Start", "End", "Name", "Score", "Strand"]
    )
    .sort_values("Chrom")
    .reset_index(drop=True)
)

chroms = (
    editing_sites_df["Chrom"]
    .value_counts()[
        (editing_sites_df["Chrom"].value_counts() >= min_known_sites_per_mock_gene)
        & (editing_sites_df["Chrom"].value_counts() <= max_known_sites_per_mock_gene)
    ]
    .sample(num_of_mock_genes, random_state=rng)
    .reset_index()["Chrom"]
    .sort_values()
    .values
)

selected_orfs_df = orfs_df.loc[orfs_df["Chrom"].isin(chroms)]

selected_editing_sites_df = editing_sites_df.loc[editing_sites_df["Chrom"].isin(chroms)]
# selected_editing_sites_df.groupby("Chrom").size()

per_chrom_editing_freqs = (
    selected_editing_sites_df.groupby("Chrom")["EditingLevel"].apply(list).values
)
per_chrom_adenosine_positions = (
    selected_editing_sites_df.groupby("Chrom")["Position"].apply(list).values
)


per_chrom_gene_names = selected_orfs_df["Name"].values
per_chrom_orf_start = selected_orfs_df["Start"].values
per_chrom_orf_end = selected_orfs_df["End"].values
per_chrom_strand = selected_orfs_df["Strand"].values

out_dir = Path("/private7/projects/Combinatorics/Simulations/GraphAssessment")

transcriptome = "/private7/projects/Combinatorics/D.pealeii/Annotations/orfs_squ.fa"


for unknown_probability in unknown_probabilities:
    for (
        chrom,
        gene_name,
        orf_start,
        orf_end,
        strand,
        editing_freqs,
        adenosine_positions,
    ) in zip(
        chroms,
        per_chrom_gene_names,
        per_chrom_orf_start,
        per_chrom_orf_end,
        per_chrom_strand,
        per_chrom_editing_freqs,
        per_chrom_adenosine_positions,
    ):
        ic(unknown_probability, chrom)

        prefix = f"{chrom}.{gene_name}.UP{str(unknown_probability).replace('.', '_')}."

        simulate_complete_and_corresponding_partially_unknown_main(
            transcriptome=transcriptome,
            chrom=chrom,
            orf_start=orf_start,
            orf_end=orf_end,
            strand=strand,
            n_reads=n_reads,
            known_editing_freqs=editing_freqs,
            adenosine_positions=adenosine_positions,
            unknown_probability=unknown_probability,
            snp_noise_level=snp_noise_level,
            top_x_noisy_positions=top_x_noisy_positions,
            assurance_factor=assurance_factor,
            out_dir=out_dir,
            known_sites_file=known_sites_file,
            out_files_sep=out_files_sep,
            denovo_detection=denovo_detection,
            common_output_prefix=prefix,
        )


# distinct_files_cmd = fr"""INFILES=$(echo Simulations/GraphAssessment/*.UniqueProteins.tsv);
# julia \
# --project=. \
# --threads {threads} --proc {processes} \
# Code/Simulations/maximal_independent_set_5.jl \
# --infiles $INFILES \
# --postfix_to_remove .UniqueProteins.tsv \
# --idcol Protein \
# --firstcolpos 15 \
# --datatype Proteins \
# --outdir Simulations/GraphAssessment \
# --fracstep 0.2 \
# --fracrepetitions 4 \
# --algrepetitions 2 \
# --algs Ascending Descending \
# --run_solve_threaded
# """

# import subprocess
# # subprocess.run(distinct_files_cmd, shell=True)
# subprocess.run(distinct_files_cmd)

# import plotly.express as px


distinct_files = list(out_dir.glob("*.DistinctUniqueProteins.*.csv"))
distinct_files_dict = {}
for distinct_file in distinct_files:
    file_name_parts = distinct_file.name.split(".")
    # file_name_parts = distinct_files[0].name.split(".")
    chrom = file_name_parts[0]
    unknown_probability = float(file_name_parts[1].replace("UP", "").replace("_", "."))
    data_type = file_name_parts[2]
    distinct_files_dict[(chrom, unknown_probability, data_type)] = distinct_file

dfs = []
for (
    chrom,
    unknown_probability,
    data_type,
), distinct_file in distinct_files_dict.items():
    df = pd.read_table(distinct_file)
    df.insert(0, "Chrom", chrom)
    df.insert(1, "UnknownProbability", unknown_probability)
    df.insert(2, "DataType", data_type)
    dfs.append(df)
merged_df = (
    pd.concat(dfs, ignore_index=True)
    .rename(
        columns={
            "NumUniqueSamples": "NumDistinctProteins",
            "UniqueSamples": "DistinctProteins",
        }
    )
    .merge(selected_orfs_df.loc[:, ["Chrom", "Name"]], on="Chrom")
)


# merged_distinct_df.columns

# fig = px.scatter(
#     merged_distinct_df, x="Fraction", y="NumDistinctProteins", color="Data"
# )

merged_df.groupby(["Name", "UnknownProbability", "DataType", "Fraction"])[
    "NumDistinctProteins"
].max()

merged_df.groupby(["Name", "UnknownProbability", "DataType"])[
    "NumDistinctProteins"
].max()