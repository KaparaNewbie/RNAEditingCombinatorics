from multiprocessing import Pool
from pathlib import Path
from itertools import product
import sys

import pandas as pd
import numpy as np
from icecream import ic

# print(Path(".").absolute())
# exit()

code_dir = "/private7/projects/Combinatorics/Code"
sys.path.append(str(Path(code_dir).absolute()))

from EditingUtils.seq import make_fasta_dict
from pileup_with_subparsers import (
    simulate_complete_and_corresponding_errored_partially_unknown_main,
)

# num. of columns to disply when printing a df in the terminal
pd.set_option("display.max_columns", 20)

processes = 4
threads = 30

# parameters for simulating reads
n_reads = 100_000
# n_reads = 10_000
num_of_mock_genes = 3
min_known_sites_per_mock_gene = 50
max_known_sites_per_mock_gene = 150
unknown_probabilities = [
    round(x / 126, 2) for x in [11, 17, 23]
]  # chance of a position's status being unknown 
error_probability = 0.001  # chance of a base being called wrong
rng = np.random.RandomState(1892) # seed for reproducibility - when selecting a gene and when inserting editing, sequencing error, missing information
data_creation_repetitions = list(range(1, 11)) # num of times to repeat the data creation process - for each complete and errored data set


# list(product(unknown_probabilities, data_creation_repetitions))

# parameters for denovo editing detection & formatting of output files
snp_noise_level = 0.1
top_x_noisy_positions = 3
assurance_factor = 1.5
denovo_detection = False
out_files_sep = "\t"
min_percent_of_max_coverage = 0.1
group_col = "Gene"
compression_postfix = ".gz"

# paths of input & output files
out_dir = Path("/private7/projects/Combinatorics/Simulations/GraphAssessment")
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
per_chrom_known_edited_adenosine_positions = (
    selected_editing_sites_df.groupby("Chrom")["Position"].apply(list).values
)

transcriptome = "/private7/projects/Combinatorics/D.pealeii/Annotations/orfs_squ.fa"
transcriptome_dict = make_fasta_dict(transcriptome)


per_chrom_gene_names = selected_orfs_df["Name"].values
per_chrom_orf_start = selected_orfs_df["Start"].values
per_chrom_orf_end = selected_orfs_df["End"].values
per_chrom_strand = selected_orfs_df["Strand"].values
per_chrom_seq = [transcriptome_dict[chrom] for chrom in chroms]








starmap_inputs = []
# for unknown_probability in unknown_probabilities:
for unknown_probability, data_creation_repetition in product(unknown_probabilities, data_creation_repetitions):
    for (
        chrom,
        gene_name,
        orf_start,
        orf_end,
        strand,
        seq,
        editing_freqs,
        known_edited_adenosine_positions,
    ) in zip(
        chroms,
        per_chrom_gene_names,
        per_chrom_orf_start,
        per_chrom_orf_end,
        per_chrom_strand,
        per_chrom_seq,
        per_chrom_editing_freqs,
        per_chrom_known_edited_adenosine_positions,
    ):
        # common_output_prefix = (
        #     f"{chrom}.{gene_name}.UP{str(unknown_probability).replace('.', '_')}."
        # )
        common_output_prefix = (
            f"{chrom}.{gene_name}.UP{str(unknown_probability).replace('.', '_')}.Rep{data_creation_repetition}."
        )

        starmap_inputs.append(
            (
                transcriptome,
                chrom,
                orf_start,
                orf_end,
                strand,
                seq,
                n_reads,
                editing_freqs,
                known_edited_adenosine_positions,
                unknown_probability,
                error_probability,
                snp_noise_level,
                top_x_noisy_positions,
                min_percent_of_max_coverage,
                assurance_factor,
                out_dir,
                known_sites_file,
                out_files_sep,
                denovo_detection,
                group_col,
                common_output_prefix,
                compression_postfix,
            )
        )


# ic(starmap_inputs)


with Pool(processes) as pool:
    pool.starmap(
        func=simulate_complete_and_corresponding_errored_partially_unknown_main,
        iterable=starmap_inputs,
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


# distinct_files = list(out_dir.glob("*.DistinctUniqueProteins.*.csv"))
# distinct_files_dict = {}
# for distinct_file in distinct_files:
#     file_name_parts = distinct_file.name.split(".")
#     # file_name_parts = distinct_files[0].name.split(".")
#     chrom = file_name_parts[0]
#     unknown_probability = float(file_name_parts[1].replace("UP", "").replace("_", "."))
#     data_type = file_name_parts[2]
#     distinct_files_dict[(chrom, unknown_probability, data_type)] = distinct_file

# dfs = []
# for (
#     chrom,
#     unknown_probability,
#     data_type,
# ), distinct_file in distinct_files_dict.items():
#     df = pd.read_table(distinct_file)
#     df.insert(0, "Chrom", chrom)
#     df.insert(1, "UnknownProbability", unknown_probability)
#     df.insert(2, "DataType", data_type)
#     dfs.append(df)
# merged_df = (
#     pd.concat(dfs, ignore_index=True)
#     .rename(
#         columns={
#             "NumUniqueSamples": "NumDistinctProteins",
#             "UniqueSamples": "DistinctProteins",
#         }
#     )
#     .merge(selected_orfs_df.loc[:, ["Chrom", "Name"]], on="Chrom")
# )


# # merged_distinct_df.columns

# # fig = px.scatter(
# #     merged_distinct_df, x="Fraction", y="NumDistinctProteins", color="Data"
# # )

# merged_df.groupby(["Name", "UnknownProbability", "DataType", "Fraction"])[
#     "NumDistinctProteins"
# ].max()

# merged_df.groupby(["Name", "UnknownProbability", "DataType"])[
#     "NumDistinctProteins"
# ].max()




def get_files_for_fp_tests(in_dir: Path,
    chrom: str,
    swiss_prot_name: str,
    unknown_probability: float,
    repetition: str,
    data_types: list[str] = ["Complete", "Errored.PartiallyUnknown"],) -> list[Path]:
    coupled_unique_proteins_files = [
            Path(
            in_dir,
            f"{chrom}.{swiss_prot_name}.UP{str(unknown_probability).replace('.', '_')}.Rep{repetition}.{data_type.replace('+', '.')}.UniqueProteins.tsv.gz",
        )
        for data_type in data_types
    ]
    for f in coupled_unique_proteins_files:
        if not f.exists():
            raise FileNotFoundError(f)
    out_file = Path(
        in_dir,
        f"{chrom}.{swiss_prot_name}.UP{str(unknown_probability).replace('.', '_')}.Rep{repetition}.FalsePositives.tsv.gz"
    )
    return coupled_unique_proteins_files + [out_file]
    

files_for_fp_tests = [
    get_files_for_fp_tests(
        out_dir, chrom, swiss_prot_name, unknown_probability, repetition
    )
    for chrom, swiss_prot_name in zip(
        chroms, per_chrom_gene_names
    )
    for unknown_probability, repetition in product(
        unknown_probabilities, data_creation_repetitions
    )
]

# len(files_for_fp_tests)

complete_infiles_fofn = Path(out_dir, "CompleteInfiles.fofn")
errored_na_files_fofn = Path(out_dir, "ErroredNAFiles.fofn")
false_positives_out_files_fofn = Path(out_dir, "FalsePositivesOutFiles.fofn")

complete_infiles = " ".join([str(f[0]) for f in files_for_fp_tests])
errored_na_files = " ".join([str(f[1]) for f in files_for_fp_tests])
false_positives_out_files = " ".join([str(f[2]) for f in files_for_fp_tests])

for fofn, files in zip(
    [complete_infiles_fofn, errored_na_files_fofn, false_positives_out_files_fofn],
    [complete_infiles, errored_na_files, false_positives_out_files],
):
    with open(fofn, "w") as f:
        f.write(files)


# complete_infiles = complete_infiles[:2]
# errored_na_files = errored_na_files[:2]
# out_files = out_files[:2]


# mis_5_assessment_cmd = f"""
# julia \
# --project=. \
# Code/Simulations/mis_5_assessment.jl \
# --complete_infiles {" ".join(complete_infiles)} \
# --errored_na_files {" ".join(errored_na_files)} \
# --out_files {" ".join(out_files)}
# """

# import subprocess
# subprocess.run(mis_5_assessment_cmd, shell=True)




