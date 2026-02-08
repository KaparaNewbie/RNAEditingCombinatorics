# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.18.1
#   kernelspec:
#     display_name: combinatorics2
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Imports

# %%
# %load_ext autoreload
# %autoreload 2

# %%
import argparse
from multiprocessing import Pool
from pathlib import Path
from typing import Union
from collections import defaultdict, Counter
import subprocess
import sys
from itertools import chain

import pandas as pd
import numpy as np
from pybedtools import BedTool
from icecream import ic
from Bio import Seq
import plotly.io as pio
import plotly.express as px

pio.templates.default = "plotly_white"

code_dir = "/private7/projects/Combinatorics/Code"
sys.path.append(str(Path(code_dir).absolute()))

# from General.argparse_utils import abs_path_from_str, expanded_path_from_str
from Alignment.alignment_utils import filter_bam_by_read_quality, sample_bam
from Pileup.mpileup import mpileup
import Pileup.positions as positions_module
# from Pileup.positions import (
#     pileup_to_positions,
#     # multisample_pileups_to_positions_old,
#     multisample_pileups_to_positions_all_transcripts,
#     simulate_complete_and_corresponding_errored_partially_unknown_positions_dfs,
# )
# from Pileup.reads import reads_and_unique_reads, multisample_reads_and_unique_reads
# from Pileup.proteins import (
#     proteins_and_unique_proteins,
#     multisample_proteins_and_unique_proteins,
# )
from General.consts import final_words, ic_prefix

# configure icecream to print the time of the print and the context (file, line, function)
ic.configureOutput(includeContext=True, prefix=ic_prefix)


DataFrameOrSeries = Union[pd.DataFrame, pd.Series]

# %% [markdown]
# # Setting variables and steps prios to editing detection

# %%
transcriptome = Path("/private7/projects/Combinatorics/O.vulgaris/Annotations/orfs_oct.fa")
# samples_table = "O.vulgaris/Data/PRJNA791920/IsoSeqPolished/samples.csv"
sample_col = ...
group_col = ...
known_editing_sites = Path("/private7/projects/Combinatorics/O.vulgaris/Annotations/O.vul.EditingSites.bed")
cds_regions = Path("/private7/projects/Combinatorics/O.vulgaris/Annotations/orfs_oct.bed")
min_percent_of_max_coverage = 0.1
min_mapped_reads_per_position = 0
include_flags = None
exclude_flags = "2304"
parity = "SE"
snp_noise_level = 0.05 # TODO or even 0.04 for pooled samples? 0.1 is for unpooled samples
top_x_noisy_positions = 3
assurance_factor = 1.5
pooled_transcript_noise_threshold = 0.06 # TODO probably lower in accordance with the lowered snp_noise_level?
# out_dir = "O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BHAfterNoise.3"
out_dir = Path("/private7/projects/Combinatorics/O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.PooledSamples")
samtools_path = "samtools"
processes = 20
threads = 5
min_rq = 0.998
min_bq = 30
out_files_sep = "\t"
keep_pileup_files = True
keep_bam_files = True
override_existing_pileup_files = False
override_existing_bam_files = False
gz_compression = True
alignments_stats_table = Path("/private7/projects/Combinatorics/O.vulgaris/Alignment/PRJNA791920/IsoSeq.Polished.Unclustered/AggregatedByChromBySampleSummary.tsv")
alignments_stats_table_sep = "\t"
min_samples = 0
min_mapped_reads_per_sample = 0
total_mapped_reads = 50
min_known_sites = 0
main_by_chrom_dir = Path("/private7/projects/Combinatorics/O.vulgaris/Alignment/PRJNA791920/IsoSeq.Polished.Unclustered/ByChrom")
postfix = ".bam"
interfix_start = ".aligned.sorted"
remove_non_refbase_noisy_positions = False
known_sites_only = False
pileup_dir_name: str = "PileupFiles"
positions_dir_name: str = "PositionsFiles"
reads_dir_name: str = "ReadsFiles"
proteins_dir_name: str = "ProteinsFiles"
sample_reads: bool = False
num_sampled_reads: int = 80_000
seed: int = 1892
final_editing_scheme = "BH after noise thresholding"
disregard_alt_base_freq_1 = True
alternative_hypothesis = "larger"
stop_after_making_noise_positions: bool = True # TODO just for the sake of testing the detection scheme
samples_are_pooled = True

# %%
prob_regions_bed = None  # we don't deal with individual "problamitic" sites in such large scale analysis

# denovo_detection = True
denovo_detection = not known_sites_only

out_dir.mkdir(exist_ok=True)

# %%
# 1 - get data

cds_df = pd.read_csv(
    cds_regions,
    sep="\t",
    names="Chrom Start End Name Score Strand".split(),
    comment="#",
)
# samples_df = pd.read_csv(samples_table, sep=samples_table_sep)
# sample_to_group_dict = {
#     row[sample_col]: row[group_col] for _, row in samples_df.iterrows()
# }
alignments_stats_df = pd.read_csv(
    alignments_stats_table, sep=alignments_stats_table_sep
)

# while redundant, as the first filtertation also accounts for cases where min_known_sites is 0,
# it allows us to accept an alignments_stats_df which doesn't have a "KnownSites" column to begin with
if min_known_sites > 0:
    alignments_stats_df = alignments_stats_df.loc[
        (alignments_stats_df["Samples"] >= min_samples)
        & (
            alignments_stats_df["MappedReadsPerSample"]
            >= min_mapped_reads_per_sample
        )
        & (alignments_stats_df["KnownSites"] >= min_known_sites)
        & (alignments_stats_df["MappedReads"] >= total_mapped_reads)
    ]
else:
    alignments_stats_df = alignments_stats_df.loc[
        (alignments_stats_df["Samples"] >= min_samples)
        & (
            alignments_stats_df["MappedReadsPerSample"]
            >= min_mapped_reads_per_sample
        )
        & (alignments_stats_df["MappedReads"] >= total_mapped_reads)
    ]

alignments_stats_df = alignments_stats_df.merge(cds_df, how="left")

# todo - remove - this is for debugging
# test_chroms = [
#     "comp183313_c0_seq12",
#     "comp162994_c0_seq1",
#     "comp183909_c0_seq7",
#     "comp181233_c0_seq10",
#     "comp183670_c0_seq3",
#     "comp183256_c0_seq35",
#     "comp183782_c0_seq5",
#     "comp183377_c0_seq11",
#     "comp181723_c2_seq2",
#     "comp169467_c0_seq1",
#     "comp183713_c0_seq9",
#     "comp181924_c0_seq4",
# ]
tests_chrom_with_more_than_3_snps = [
    "comp176758_c0_seq1",
    "comp73217_c0_seq1",
    "comp148544_c0_seq1",
    "comp72138_c0_seq1",
    "comp183760_c0_seq7",
]
test_chroms_with_3_or_less_snps = [
    "comp182826_c6_seq4",
    "comp72421_c0_seq1",
    "comp179489_c1_seq1",
    "comp181820_c2_seq1",
    "comp73741_c0_seq1",
]
test_chroms = tests_chrom_with_more_than_3_snps + test_chroms_with_3_or_less_snps
alignments_stats_df = alignments_stats_df.loc[
    alignments_stats_df["Chrom"].isin(test_chroms)
]

alignments_stats_df

# %%
chroms = []
swissprot_names = []
samples = []
# groups = []
orfs_starts = []
orfs_ends = []
pileup_formatted_regions = []
orfs_strands = []
bam_files = []
for _, row in alignments_stats_df.iterrows():
    chrom = row["Chrom"]
    swissprot_name = row["Name"]
    start = row["Start"]
    end = row["End"]
    strand = row["Strand"]
    region = f"{chrom}:{start+1}-{end}"
    bams_in_chrom = list(
        Path(main_by_chrom_dir, chrom).glob(f"*{postfix}")
    )  # should be something like `*.bam`
    # ic(bams_in_chrom)  # todo - remove
    for bam_in_chrom in bams_in_chrom:
        # ic(bam_in_chrom)
        sample = bam_in_chrom.stem[: bam_in_chrom.stem.index(interfix_start)]
        # group = sample_to_group_dict[sample]
        chroms.append(chrom)
        swissprot_names.append(swissprot_name)
        samples.append(sample)
        # groups.append(group)
        orfs_starts.append(start)
        orfs_ends.append(end)
        pileup_formatted_regions.append(region)
        orfs_strands.append(strand)
        bam_files.append(bam_in_chrom)

ic(chroms, swissprot_names, samples, orfs_starts, orfs_ends, pileup_formatted_regions, orfs_strands, bam_files);

# %%
# 2 - filter reads by mapping quality

if min_rq is not None:
    filterd_bams_dir = Path(out_dir, "FilterdBAMs")
    filterd_bams_dir.mkdir(exist_ok=True)
    with Pool(processes=processes) as pool:
        filtered_bam_files = pool.starmap(
            func=filter_bam_by_read_quality,
            iterable=[
                (
                    samtools_path,
                    in_bam,
                    min_rq,
                    threads,
                    filterd_bams_dir,
                    override_existing_bam_files,
                )
                for in_bam in bam_files
            ],
        )
else:
    filtered_bam_files = bam_files
    
ic(len(filtered_bam_files));
# filtered_bam_files

# %%
# 3 - sample reads to desired number of reads

if sample_reads:
    sampled_bams_dir = Path(out_dir, "SampledBAMs")
    sampled_bams_dir.mkdir(exist_ok=True)
    sampled_filtered_bam_files = [
        Path(sampled_bams_dir, f"{bam_file.stem}.Sampled{num_sampled_reads}.bam")
        for bam_file in filtered_bam_files
    ]
    with Pool(processes=processes) as pool:
        pool.starmap(
            func=sample_bam,
            iterable=[
                (samtools_path, in_bam, out_bam, num_sampled_reads, seed, threads)
                for in_bam, out_bam in zip(
                    filtered_bam_files, sampled_filtered_bam_files
                )
            ],
        )
else:
    sampled_filtered_bam_files = filtered_bam_files
    
ic(len(sampled_filtered_bam_files));


# %%
# 3 - run mpileup

# pileup_dir = Path(out_dir, "PileupFiles")
pileup_dir = Path(out_dir, pileup_dir_name)
pileup_dir.mkdir(exist_ok=True)
pileup_files = [
    Path(pileup_dir, f"{sample}.{chrom}.pileup")
    for sample, chrom in zip(samples, chroms)
]

with Pool(processes=processes) as pool:
    pool.starmap(
        func=mpileup,
        iterable=[
            (
                samtools_path,
                transcriptome,
                region,
                include_flags,
                exclude_flags,
                min_bq,
                in_bam,
                out_pileup,
                threads,
                override_existing_pileup_files,
            )
            for region, in_bam, out_pileup in zip(
                pileup_formatted_regions,
                # filtered_bam_files,
                sampled_filtered_bam_files,
                pileup_files,
            )
        ],
    )

if not keep_bam_files:
    if min_rq is not None:
        subprocess.run(f"rm -rf {filterd_bams_dir}", shell=True)
    # delete sampled reads, whether they were previously filtered by base quality or not
    if sample_reads:
        subprocess.run(f"rm -rf {sampled_bams_dir}", shell=True)

# %%
compression_postfix = ".gz" if gz_compression else ""

# 4 - pileup files -> positions dfs

group_col = "Transcript"  # todo make more visible

unique_chroms = []
unique_orfs_starts = []
unique_orfs_ends = []
unique_orfs_strands = []
unique_swissprot_names = []
for chrom, start, end, strand, swissprot_name in zip(
    chroms, orfs_starts, orfs_ends, orfs_strands, swissprot_names
):
    if chrom in unique_chroms:
        continue
    unique_chroms.append(chrom)
    unique_orfs_starts.append(start)
    unique_orfs_ends.append(end)
    unique_orfs_strands.append(strand)
    unique_swissprot_names.append(swissprot_name)

# positions_dir_name = (
#     "PositionsFiles" if denovo_detection else "PositionsFilesKnownSites"
# )
positions_dir = Path(out_dir, positions_dir_name)
positions_dir.mkdir(exist_ok=True)

reads_mapping_files = [
    Path(positions_dir, f"{chrom}.OldToNewReads.csv{compression_postfix}")
    for chrom in unique_chroms
]
mismatches_files = [
    Path(positions_dir, f"{chrom}.Mismatches.csv{compression_postfix}")
    for chrom in unique_chroms
]
positions_files = [
    Path(positions_dir, f"{chrom}.positions.csv{compression_postfix}")
    for chrom in unique_chroms
]
corrected_noise_files = [
    Path(positions_dir, f"{chrom}.CorrectedNoise.csv{compression_postfix}")
    for chrom in unique_chroms
]
corrected_editing_files = [
    Path(positions_dir, f"{chrom}.CorrectedEditing.csv{compression_postfix}")
    for chrom in unique_chroms
]

pileup_files_per_chroms = defaultdict(list)
samples_per_chroms = defaultdict(list)
for sample, chrom, pileup_file in zip(samples, chroms, pileup_files):
    pileup_files_per_chroms[chrom].append(pileup_file)
    samples_per_chroms[chrom].append(sample)

binom_noise_pval_col = "NoiseBinomPVal"
bh_noise_pval_col = "NoiseCorrectedPVal"
bh_noisy_col = "NoisyCorrected"
binom_editing_pval_col = "EditingBinomPVal"
bh_editing_pval_col = "EditingCorrectedPVal"
bh_editing_col = "EditedCorrected"

# %%
# multisample_pileups_to_positions_all_transcripts(
#     processes,
#     unique_chroms,
#     unique_orfs_strands,
#     pileup_files_per_chroms,
#     samples_per_chroms,
#     reads_mapping_files,
#     positions_files,
#     corrected_noise_files,
#     corrected_editing_files,
#     min_percent_of_max_coverage,
#     min_mapped_reads_per_position,  # min_absolute_coverage
#     snp_noise_level,
#     top_x_noisy_positions,
#     assurance_factor,
#     transcriptome,
#     final_editing_scheme,
#     prob_regions_bed,  # problamatic_regions_file
#     known_editing_sites,  # known_sites_file
#     cds_regions,  # cds_regions_file
#     out_files_sep,
#     keep_pileup_files,
#     remove_non_refbase_noisy_positions,
#     denovo_detection,
#     alternative_hypothesis,
#     binom_noise_pval_col,
#     bh_noise_pval_col,
#     bh_noisy_col,
#     binom_editing_pval_col,
#     bh_editing_pval_col,
#     bh_editing_col,
#     disregard_alt_base_freq_1,
# )

# for the sake of testing the detection scheme for pooled samples in the notebook,
# where we dont go beyond noise and editing positions,
# we change the variables names to those use in multisample_pileups_to_positions_all_transcripts
# so can copy lines from it directly

processes: int = processes
chroms: list[str] = unique_chroms
strands: list[str] = unique_orfs_strands
pileup_files_per_chroms: defaultdict[str, list[Path]] = pileup_files_per_chroms
samples_per_chroms: defaultdict[str, list[str]] = samples_per_chroms 
reads_mapping_files: list[Path] = reads_mapping_files
positions_files: list[Path] = positions_files
corrected_noise_files: list[Path] = corrected_noise_files
corrected_editing_files: list[Path] = corrected_editing_files
min_percent_of_max_coverage: float = min_percent_of_max_coverage
min_absolute_coverage: Union[int, None] = min_mapped_reads_per_position
snp_noise_level: float = snp_noise_level
top_x_noisy_positions: int = top_x_noisy_positions
assurance_factor: float = assurance_factor
transcriptome_file: Path = transcriptome
final_editing_scheme: str = final_editing_scheme
problamatic_regions_file: Union[Path, str, None] = prob_regions_bed
known_sites_file: Union[Path, str, None] = known_editing_sites
cds_regions_file: Union[Path, str, None] = cds_regions
out_files_sep: str = out_files_sep
keep_pileup_files: bool = keep_pileup_files
remove_non_refbase_noisy_positions: bool = remove_non_refbase_noisy_positions
denovo_detection: bool = denovo_detection
alternative_hypothesis: str = alternative_hypothesis
binom_noise_pval_col: str = binom_noise_pval_col
bh_noise_pval_col: str = bh_noise_pval_col
bh_noisy_col: str = bh_noisy_col
binom_editing_pval_col: str = binom_editing_pval_col
bh_editing_pval_col: str = bh_editing_pval_col
bh_editing_col: str = bh_editing_col
disregard_alt_base_freq_1: bool = disregard_alt_base_freq_1

# %%

# %% [markdown]
# # Editing detection

# %% [markdown]
# ## Main

# %%
cds_regions_df = pd.read_csv(
    cds_regions_file,
    sep="\t",
    names=["Chrom", "Start", "End", "Name", "Score", "Strand"],
)
cds_regions_df

# %%
transcriptome_dict = positions_module.make_fasta_dict(transcriptome_file)

# %% [markdown]
# ## multisample_pileups_to_positions_part_1 - test

# %%
# # get the positions dfs up until noise correction (not included)
# with Pool(processes=processes) as pool:
#     pool.starmap(
#         func=positions_module.multisample_pileups_to_positions_part_1,
#         iterable=[
#             (
#                 pileup_files_per_chroms[chrom], # pileup_files: list[Path]
#                 samples_per_chroms[chrom], # samples: list[str]
#                 strand, # strand: str
#                 positions_file, # positions_file: Path
#                 problamatic_regions_file, # problamatic_regions_file: Union[Path, str, None] = None
#                 known_sites_file, # known_sites_file: Union[Path, str, None] = None
#                 cds_regions_file, # cds_regions_file: Union[Path, str, None] = None
#                 reads_mapping_file, # reads_mapping_file: Union[Path, str, None] = None
#                 out_files_sep, # out_files_sep: str = "\t"
#             )
#             for reads_mapping_file, chrom, strand, positions_file in zip(
#                 reads_mapping_files, chroms, strands, positions_files
#             )
#         ],
#     )

# %%
chroms.index("comp181924_c0_seq4")

# %%
# take one chrom to test

# i = 0
i = chroms.index("comp181924_c0_seq4") # 7

reads_mapping_file = reads_mapping_files[i]
chrom = chroms[i]
strand = strands[i]
positions_file = positions_files[i]
mismatches_file = mismatches_files[i]

pileup_files = pileup_files_per_chroms[chrom]
samples = samples_per_chroms[chrom]

ic(reads_mapping_file, chrom, strand, positions_file, pileup_files, samples);

# %%
cols = [
    "Chrom",
    "Position",
    "RefBase",
    "TotalCoverage",
    "MappedBases",
    "Phred",
    "Reads",
]

positions_dfs = [
    pd.read_csv(pileup_file, sep="\t", names=cols) for pileup_file in pileup_files
]
for positions_df, sample in zip(positions_dfs, samples):
    positions_df.insert(0, "Sample", sample)

# merge the separate pileups tables into a multi-sample one
positions_df = pd.concat(positions_dfs, ignore_index=True)

# # filter out zero-coverage positions (probably due to deletions) - they are irrelevent, and also disturb `replace_reads_names` later
positions_df = positions_df.loc[positions_df["TotalCoverage"] > 0]
# if len(positions_df) == 0:
#     return
# todo: in future versions, also remove positions with insufficient coverage of A, T, C, or G bases
# (such positions are not informative)
# then, define the total coverage as the sum of the coverage of A, T, C, and G bases, and ignore
# unknown bases?
# positions_df["ATCGs"] = positions_df.loc[:, ["A", "T", "C", "G"]].sum(axis=1)
# positions_df = positions_df.loc[positions_df["ATCGs"] > 0]
# del positions_df["ATCGs"]
# if len(positions_df) == 0:
#     return

# del positions_df["Phred"]

# change editing position to 0-based
positions_df["Position"] = positions_df["Position"] - 1

# present all reads as if they were mapped to the positive strand
positions_df["MappedBases"] = positions_df["MappedBases"].str.upper()
positions_df["MappedBases"] = positions_df["MappedBases"].replace(
    {r"[,]": ".", r"[<]": r">"}, regex=True
)

# make sure we know how to deal with all kinds of mapped bases
unique_mapped_bases = set(chain.from_iterable(positions_df["MappedBases"]))
if not unique_mapped_bases <= {">", "*", ".", "A", "C", "G", "T", "N"}:
    raise Exception(f"{unique_mapped_bases = }")

# verify that the number of mapped bases (6th col) corresponds
# to the number of reads' names (8th col)

for row_num, row in enumerate(positions_df.itertuples()):
    bases = row.MappedBases
    reads = row.Reads.split(",")
    if len(bases) != len(reads):
        sample = row.Sample
        # get the first element of a signle-element list
        pileup_file = [
            pileup_file
            for pileup_file in pileup_files
            if sample in pileup_file.name
        ][0]
        raise ValueError(
            f"Line {row_num} in {pileup_file} contains indels and/or doesn't have the "
            "same number of mapped bases and reads."
        )

# replace reads' names with a shortened, memory-efficient version
positions_module.replace_reads_names(
    positions_df,
    "TotalCoverage",
    "Reads",
    reads_sep=",",
    mapping_out_file=reads_mapping_file,
    out_files_sep=out_files_sep,
)

positions_df

# %%
positions_df["Position"].value_counts()

# %%
# if len(positions_df) > 0:
positions_df["Sample"] = positions_df["Sample"].astype(str)
positions_df["MappedBases"] = positions_df["MappedBases"].astype(str)
positions_df["Phred"] = positions_df["Phred"].astype(str)
positions_df["Reads"] = positions_df["Reads"].astype(str)

if len(positions_df) > 0:
    positions_df.insert(
        positions_df.columns.get_loc("Sample") + 1,
        "Samples",
        positions_df.apply(
            lambda x: ",".join([x["Sample"] for _ in range(len(x["MappedBases"]))]),
            axis=1,
        ),
    )
else:
    # add an empty "Samples" column in the empty df for the next step
    positions_df.insert(positions_df.columns.get_loc("Sample") + 1, "Samples", None)
    
positions_df

# %%
# positions_df.loc[
#     positions_df["Position"] == 15717
# ]

# %%
# join multi-sample data per position
try:
    positions_df = (
        positions_df.groupby(["Chrom", "Position", "RefBase"])
        .agg(
            {
                "Samples": ",".join,
                "TotalCoverage": "sum", # TODO update to "sum" in the original function
                "MappedBases": "".join,
                "Phred": "".join,
                "Reads": ",".join,
            }
        )
        .reset_index()
    )
except TypeError as e:
    ic()
    ic(e)
    ic(pileup_files)
    
positions_df

# %%
# # remove positions with insufficient coverage
# max_coverage = positions_df["TotalCoverage"].max()
# required_coverage = max_coverage * min_percent_of_max_coverage
# positions_df = positions_df.loc[positions_df["TotalCoverage"] >= required_coverage]

# %%
positions_module.annotate_prolamatic_sites(positions_df, strand, problamatic_regions_file)
positions_module.annotate_known_sites(positions_df, strand, known_sites_file)
positions_module.annotate_coding_sites(positions_df, strand, cds_regions_file)
positions_module.annotate_base_counts(positions_df)

# todo: in future versions, explictly define mismatch per position and, based on that,
# only then calculate noise
positions_module.annotate_noise(positions_df, strand)

# %%
positions_df


# %%
def count_mismatches_per_position(
    ref_base, a_count, t_count, c_count, g_count
):
    """
    Count mismatches in one position.
    """
    bases = ["A", "T", "C", "G"]
    base_counts = [a_count, t_count, c_count, g_count]

    mismatches_counts = {}
    for alt_base, count in zip(bases, base_counts):
        # if alt_base != ref_base:
        #     mismatches_counts[f"{ref_base}{alt_base}"] = count
        mismatches_counts[f"{ref_base}{alt_base}"] = count
    
    mismatches_counts = Counter(mismatches_counts)
    return mismatches_counts


def count_mismatches(
    positions_df: pd.DataFrame,
    chrom: str,
    mismatches_file: Path | str,
    mismatches_file_sep: str,
):
    """
    Count mismatches in all positions and save them to a csv file.
    """
    per_position_mismatches_counters = positions_df.apply(
        lambda row: count_mismatches_per_position(
            row["RefBase"],
            row["A"],
            row["T"],
            row["C"],
            row["G"],
        ),
        axis=1
    ).tolist()

    mismatches_counter = per_position_mismatches_counters[0]
    for counter in per_position_mismatches_counters[1:]:
        mismatches_counter += counter

    true_mismatches_counter = Counter(
        {
        mismatch: mismatches_counter[mismatch]
            for mismatch in mismatches_counter
            if mismatch[0] != mismatch[1]  
        }
    )
    
    all_possible_mismatches = {
        f"{ref_base}{alt_base}"
        for ref_base in "ATCG"
        for alt_base in "ATCG"
        if ref_base != alt_base
    }
    for mismatch in all_possible_mismatches:
        if mismatch not in true_mismatches_counter:
            true_mismatches_counter[mismatch] = 0
        
    mismatches_df = pd.DataFrame.from_dict(
        true_mismatches_counter, orient="index", columns=["Count"]
    ).reset_index()
    mismatches_df = mismatches_df.rename(columns={"index": "Mismatch"})
    mismatches_df = mismatches_df.sort_values(
        ["Count", "Mismatch"], 
        ascending=[False, True]
    )
    mismatches_df.insert(0, "Chrom", chrom)
    
    mismatches_df.to_csv(mismatches_file, sep=mismatches_file_sep, index=False)

    # return mismatches_df


# %%
try:
    if not positions_df.empty:
        count_mismatches(
            positions_df.loc[:, ["RefBase", "A", "T", "C", "G"]],
            # positions_df["Chrom"][0],
            positions_df["Chrom"].iloc[0],
            mismatches_file,
            out_files_sep,
        )
    else:
        print(f"{pileup_file} is empty")
except e:
    print(positions_df.shape)
    print(positions_df.head())

# %%
positions_df["Chrom"][0]

# %%
count_mismatches(
    positions_df.loc[:, ["RefBase", "A", "T", "C", "G"]],
    positions_df["Chrom"][0],
    mismatches_file,
    out_files_sep,
)

# mismatches_df

# %%
positions_df["Chrom"].iloc[0]

# %%
mismatches_counter.most_common(2)

# %%
mismatches_counter.most_common(2)

# %%

# %% [markdown]
# ## multisample_pileups_to_positions_part_1 - run

# %%
# get the positions dfs up until noise correction (not included)
with Pool(processes=processes) as pool:
    pool.starmap(
        func=positions_module.multisample_pileups_to_positions_part_1,
        iterable=[
            (
                pileup_files_per_chroms[chrom],
                samples_per_chroms[chrom],
                strand,
                positions_file,
                problamatic_regions_file,
                known_sites_file,
                cds_regions_file,
                reads_mapping_file,
                mismatches_file,
                out_files_sep,
            )
            for reads_mapping_file, mismatches_file, chrom, strand, positions_file in zip(
                reads_mapping_files, mismatches_files, chroms, strands, positions_files
            )
        ],
    )

# %% [markdown]
# ## binom_and_bh_correction_for_noise_all_transcripts

# %%
positions_module.binom_and_bh_correction_for_noise_all_transcripts(
    positions_files,
    corrected_noise_files,
    out_files_sep,
    chroms,
    processes,
    cds_regions_df,
    transcriptome_dict,
    alternative_hypothesis,
    binom_noise_pval_col,
    bh_noise_pval_col,
    bh_noisy_col,
)

# %% [markdown]
# ## multisample_pileups_to_positions_part_2 - 4+ SNPs test

# %%
# tests_chrom_with_more_than_3_snps = [
#         "comp176758_c0_seq1",
#         "comp73217_c0_seq1",
#         "comp148544_c0_seq1",
#         "comp72138_c0_seq1",
#         "comp183760_c0_seq7",
#     ]

# %%

# %%
sep = "\t"
max_snps_per_gene_to_allow_editing_detection = 3

# %%
# take one chrom to test

i = 1
# i = chroms.index("comp181924_c0_seq4") # 7

positions_file = positions_files[i]
corrected_noise_file = corrected_noise_files[i]
strand = strands[i]
chrom = chroms[i]

ic(positions_file, corrected_noise_file, strand, chrom);

# %%
ref_base = "A" if strand == "+" else "T"

# read the positions df
positions_df = pd.read_csv(positions_file, sep=sep, dtype={"Reads": str})

# read the corrected noise df
corrected_noise_df = pd.read_csv(corrected_noise_file, sep=sep)

# merge the corrected noise into the positions df
if not corrected_noise_df.empty:
    positions_df = positions_df.merge(
        corrected_noise_df.loc[
            :,
            [
                "Chrom",
                "Position",
                binom_noise_pval_col,
                bh_noise_pval_col,
                bh_noisy_col,
            ],
        ],
        on=["Chrom", "Position"],
        how="left",
    )
    positions_df[bh_noisy_col] = positions_df[bh_noisy_col].fillna(
        False
    )  # because all([np.nan]) == True

else:
    positions_df[binom_noise_pval_col] = positions_df["RefBase"].apply(
        lambda x: 1.0 if x != ref_base else np.nan
    )
    positions_df[bh_noise_pval_col] = positions_df["RefBase"].apply(
        lambda x: 1.0 if x != ref_base else np.nan
    )
    # because all([np.nan]) == True
    positions_df[bh_noisy_col] = False

positions_df.insert(
    positions_df.columns.get_loc(bh_noisy_col) + 1,
    "BelowNoiseFreq1",
    positions_df["Noise"] < 1,
)

if disregard_alt_base_freq_1:
    noise_cols = [bh_noisy_col, "BelowNoiseFreq1"]
else:
    noise_cols = [bh_noisy_col]

positions_module.annotate_noisy_sites_final_decision(
    positions_df, positions_df.columns.get_loc("BelowNoiseFreq1") + 1, noise_cols
)

# %%
positions_df

# %%
#  possibly remove non refbase positions with too-high noise & define noise threshold
if remove_non_refbase_noisy_positions:
    # # remove non refbase positions with too-high noise or whose corrected noise p-value is not significant
    # positions_before = len(positions_df)
    # positions_df = positions_df.loc[
    #     (positions_df["RefBase"] == ref_base)
    #     | ((positions_df["Noise"] < snp_noise_level) & (positions_df["NoisyFinal"]))
    # ]
    # positions_after = len(positions_df)
    # removed_positions = positions_before - positions_after
    # print(
    #     f"{removed_positions} extremely- or insignificantly-noisy positions removed"
    # )
    # # define noise threshold
    # noise_levels = (
    #     positions_df["Noise"]
    #     .sort_values(ascending=False)[:top_x_noisy_positions]
    #     .tolist()
    # )
    # # if there are less noisy positions than `top_x_noisy_positions`, add zeros accordingly
    # noise_levels = pd.Series(
    #     noise_levels + [0 for _ in range(top_x_noisy_positions - len(noise_levels))]
    # )
    # noise_threshold = noise_levels.mean()
    raise NotImplementedError(
        "Removing non-refbase noisy positions is no longer an option."
    )
else:
    # define noise threshold
    noise_levels = (
        positions_df.loc[
            (positions_df["Noise"] < snp_noise_level)
            & (positions_df["NoisyFinal"]),
            "Noise",
        ]
        .sort_values(ascending=False)[:top_x_noisy_positions]
        .tolist()
    )
    # if there are less noisy positions than `top_x_noisy_positions`, add zeros accordingly
    noise_levels = pd.Series(
        noise_levels + [0 for _ in range(top_x_noisy_positions - len(noise_levels))]
    )
    noise_threshold = noise_levels.mean()
if pd.isna(noise_threshold):
    noise_threshold = 0

editing_detection_possible = True

# note 1:
# we require `>` rather than `>=`
# in order to enforce processing the positions into reads
# by setting `top_x_noisy_positions = 1`
# (or any value of `snp_noise_level` that will lead to
# `pooled_transcript_noise == pooled_transcript_noise_threshold`,
# for that matter)
# note 2:
# this pooled_transcript_noise_threshold was previously defined downstream,
# when converting positions to reads, proteins, etc.
# but it does belong already here as we are defining the noise threshold here
# and if the noise is too high, we don't want to consider this gene as edited
# and specifically - none of its positions
if noise_threshold > pooled_transcript_noise_threshold:
    # raise Warning(
    #     f"Warning: {positions_file = } has {noise_threshold = } > {pooled_transcript_noise_threshold = } "
    #     "and thus won't be considered for editing detection and further processing downstream."
    # )
    print(
        f"Warning: {positions_file = } has {noise_threshold = } > {pooled_transcript_noise_threshold = } "
        "and thus won't be considered for editing detection and further processing downstream."
    )
    editing_detection_possible = False
# also, don't consider a gene as edited if there are too many SNPs in it
num_of_snps_in_gene = positions_df.loc[
    (positions_df["NoisyFinal"])
    & (positions_df["Noise"].ge(snp_noise_level))
].shape[0]
if num_of_snps_in_gene > max_snps_per_gene_to_allow_editing_detection:
    # raise Warning(
    #     f"Warning: {positions_file = } has {num_of_snps_in_gene = } > {max_snps_per_gene_to_allow_editing_detection = } "
    #     "and thus won't be considered for editing detection and further processing downstream."
    # )
    print(
        f"Warning: {positions_file = } has {num_of_snps_in_gene = } > {max_snps_per_gene_to_allow_editing_detection = } "
        "and thus won't be considered for editing detection and further processing downstream."
    )
    editing_detection_possible = False

ic(editing_detection_possible)

# anyway, we finalize the noise threshold
# (although it won't matter if editing_detection_possible == False)
noise_threshold *= assurance_factor

noise_threshold

# %%
# write noise threshold to file
chrom = positions_file.name.split('.')[0]
with open(
    positions_file.parent / f"{chrom}.NoiseThreshold.csv",
    "w",
) as f:
    f.write(f"{chrom}\t{noise_threshold}\n")

# %%
# annotate editing frequency and "naive" editing status according to noise threshold
positions_module.annotate_editing_frequency_per_position(positions_df, strand)

# %%
positions_df

# %%
positions_df.loc[positions_df["RefBase"] == ref_base]

# %%
ref_base = "A" if strand == "+" else "T"

if editing_detection_possible:
    if denovo_detection:
        edited_positions = positions_df.loc[positions_df["RefBase"] == ref_base].apply(
            lambda x: x["EditingFrequency"] > noise_threshold, axis=1
        )
    else:
        edited_positions = positions_df.loc[positions_df["RefBase"] == ref_base].apply(
            lambda x: (x["EditingFrequency"] > noise_threshold) and (x["KnownEditing"]),
            axis=1,
        )
    edited = [
        i in edited_positions.loc[edited_positions].index for i in positions_df.index
    ]
else:
    edited = [False] * len(positions_df)
positions_df.insert(positions_df.columns.get_loc("Position") + 1, "Edited", edited)

positions_df

# %%
positions_df.loc[positions_df["RefBase"] == ref_base, "Edited"].value_counts()

# %%

# annotate_edited_sites(
#     positions_df, strand, noise_threshold, denovo_detection,
#     editing_detection_possible
# )

# %% [markdown]
# ## multisample_pileups_to_positions_part_2 - <3 SNPs test

# %%
# test_chroms_with_3_or_less_snps = [
#         "comp182826_c6_seq4",
#         "comp72421_c0_seq1",
#         "comp179489_c1_seq1",
#         "comp181820_c2_seq1",
#         "comp73741_c0_seq1",
#     ]

# %%
sep = "\t"
max_snps_per_gene_to_allow_editing_detection = 3

# %%
# take one chrom to test

i = 0
# i = chroms.index("comp181924_c0_seq4") # 7

positions_file = positions_files[i]
corrected_noise_file = corrected_noise_files[i]
strand = strands[i]
chrom = chroms[i]

ic(positions_file, corrected_noise_file, strand, chrom);

# %%
ref_base = "A" if strand == "+" else "T"

# read the positions df
positions_df = pd.read_csv(positions_file, sep=sep, dtype={"Reads": str})

# read the corrected noise df
corrected_noise_df = pd.read_csv(corrected_noise_file, sep=sep)

# merge the corrected noise into the positions df
if not corrected_noise_df.empty:
    positions_df = positions_df.merge(
        corrected_noise_df.loc[
            :,
            [
                "Chrom",
                "Position",
                binom_noise_pval_col,
                bh_noise_pval_col,
                bh_noisy_col,
            ],
        ],
        on=["Chrom", "Position"],
        how="left",
    )
    positions_df[bh_noisy_col] = positions_df[bh_noisy_col].fillna(
        False
    )  # because all([np.nan]) == True

else:
    positions_df[binom_noise_pval_col] = positions_df["RefBase"].apply(
        lambda x: 1.0 if x != ref_base else np.nan
    )
    positions_df[bh_noise_pval_col] = positions_df["RefBase"].apply(
        lambda x: 1.0 if x != ref_base else np.nan
    )
    # because all([np.nan]) == True
    positions_df[bh_noisy_col] = False

positions_df.insert(
    positions_df.columns.get_loc(bh_noisy_col) + 1,
    "BelowNoiseFreq1",
    positions_df["Noise"] < 1,
)

if disregard_alt_base_freq_1:
    noise_cols = [bh_noisy_col, "BelowNoiseFreq1"]
else:
    noise_cols = [bh_noisy_col]

positions_module.annotate_noisy_sites_final_decision(
    positions_df, positions_df.columns.get_loc("BelowNoiseFreq1") + 1, noise_cols
)

# %%
positions_df

# %%
#  possibly remove non refbase positions with too-high noise & define noise threshold
if remove_non_refbase_noisy_positions:
    # # remove non refbase positions with too-high noise or whose corrected noise p-value is not significant
    # positions_before = len(positions_df)
    # positions_df = positions_df.loc[
    #     (positions_df["RefBase"] == ref_base)
    #     | ((positions_df["Noise"] < snp_noise_level) & (positions_df["NoisyFinal"]))
    # ]
    # positions_after = len(positions_df)
    # removed_positions = positions_before - positions_after
    # print(
    #     f"{removed_positions} extremely- or insignificantly-noisy positions removed"
    # )
    # # define noise threshold
    # noise_levels = (
    #     positions_df["Noise"]
    #     .sort_values(ascending=False)[:top_x_noisy_positions]
    #     .tolist()
    # )
    # # if there are less noisy positions than `top_x_noisy_positions`, add zeros accordingly
    # noise_levels = pd.Series(
    #     noise_levels + [0 for _ in range(top_x_noisy_positions - len(noise_levels))]
    # )
    # noise_threshold = noise_levels.mean()
    raise NotImplementedError(
        "Removing non-refbase noisy positions is no longer an option."
    )
else:
    # define noise threshold
    noise_levels = (
        positions_df.loc[
            (positions_df["Noise"] < snp_noise_level)
            & (positions_df["NoisyFinal"]),
            "Noise",
        ]
        .sort_values(ascending=False)[:top_x_noisy_positions]
        .tolist()
    )
    # if there are less noisy positions than `top_x_noisy_positions`, add zeros accordingly
    noise_levels = pd.Series(
        noise_levels + [0 for _ in range(top_x_noisy_positions - len(noise_levels))]
    )
    noise_threshold = noise_levels.mean()
if pd.isna(noise_threshold):
    noise_threshold = 0

editing_detection_possible = True

# note 1:
# we require `>` rather than `>=`
# in order to enforce processing the positions into reads
# by setting `top_x_noisy_positions = 1`
# (or any value of `snp_noise_level` that will lead to
# `pooled_transcript_noise == pooled_transcript_noise_threshold`,
# for that matter)
# note 2:
# this pooled_transcript_noise_threshold was previously defined downstream,
# when converting positions to reads, proteins, etc.
# but it does belong already here as we are defining the noise threshold here
# and if the noise is too high, we don't want to consider this gene as edited
# and specifically - none of its positions
if noise_threshold > pooled_transcript_noise_threshold:
    # raise Warning(
    #     f"Warning: {positions_file = } has {noise_threshold = } > {pooled_transcript_noise_threshold = } "
    #     "and thus won't be considered for editing detection and further processing downstream."
    # )
    print(
        f"Warning: {positions_file = } has {noise_threshold = } > {pooled_transcript_noise_threshold = } "
        "and thus won't be considered for editing detection and further processing downstream."
    )
    editing_detection_possible = False
# also, don't consider a gene as edited if there are too many SNPs in it
num_of_snps_in_gene = positions_df.loc[
    (positions_df["NoisyFinal"])
    & (positions_df["Noise"].ge(snp_noise_level))
].shape[0]
if num_of_snps_in_gene > max_snps_per_gene_to_allow_editing_detection:
    # raise Warning(
    #     f"Warning: {positions_file = } has {num_of_snps_in_gene = } > {max_snps_per_gene_to_allow_editing_detection = } "
    #     "and thus won't be considered for editing detection and further processing downstream."
    # )
    print(
        f"Warning: {positions_file = } has {num_of_snps_in_gene = } > {max_snps_per_gene_to_allow_editing_detection = } "
        "and thus won't be considered for editing detection and further processing downstream."
    )
    editing_detection_possible = False

ic(editing_detection_possible)

# anyway, we finalize the noise threshold
# (although it won't matter if editing_detection_possible == False)
noise_threshold *= assurance_factor

noise_threshold

# %%
# write noise threshold to file
chrom = positions_file.name.split('.')[0]
with open(
    positions_file.parent / f"{chrom}.NoiseThreshold.csv",
    "w",
) as f:
    f.write(f"{chrom}\t{noise_threshold}\n")

# %%
# annotate editing frequency and "naive" editing status according to noise threshold
positions_module.annotate_editing_frequency_per_position(positions_df, strand)

# %%
positions_df

# %%
positions_df.loc[positions_df["RefBase"] == ref_base]

# %%
ref_base = "A" if strand == "+" else "T"

if editing_detection_possible:
    if denovo_detection:
        edited_positions = positions_df.loc[positions_df["RefBase"] == ref_base].apply(
            lambda x: x["EditingFrequency"] > noise_threshold, axis=1
        )
    else:
        edited_positions = positions_df.loc[positions_df["RefBase"] == ref_base].apply(
            lambda x: (x["EditingFrequency"] > noise_threshold) and (x["KnownEditing"]),
            axis=1,
        )
    edited = [
        i in edited_positions.loc[edited_positions].index for i in positions_df.index
    ]
else:
    edited = [False] * len(positions_df)
positions_df.insert(positions_df.columns.get_loc("Position") + 1, "Edited", edited)

positions_df

# %%
positions_df.loc[positions_df["RefBase"] == ref_base, "Edited"].value_counts()

# %%

# annotate_edited_sites(
#     positions_df, strand, noise_threshold, denovo_detection,
#     editing_detection_possible
# )
