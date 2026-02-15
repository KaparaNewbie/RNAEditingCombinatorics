# ---
# jupyter:
#   jupytext:
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
# code_dir = "/private7/projects/Combinatorics/Code"

# %%
# from itertools import chain, product
from itertools import chain
from pathlib import Path
# import sys
from multiprocessing import Pool

from icecream import ic
import plotly.express as px
# import plotly.graph_objects as go
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact, chi2_contingency
from statsmodels.stats.multitest import fdrcorrection

# sys.path.append(str(Path(code_dir).absolute()))
# import 
# from Pileup.positions import read_positions_file

# %% [markdown]
# # Inputs

# %%
PROCESSES = 8

# %%
platform = "PacBio"
parity = "SE"
seed = 1892
edited_col = "EditedFinal"
new_noise_positions_in_use = False
snp_noise_level = 0.05
final_noise_col = "NoisyFinal"
multisample = True
top_x_noisy_positions = 3
assurance_factor = 1.5

# %%
main_data_dir = Path(
    "/private6/projects/Combinatorics/O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.PooledSamples"
)

positions_dir = Path(main_data_dir, "PositionsFiles")
reads_dir = Path(main_data_dir, "ReadsFiles")

# alignment_stats_file = "/private7/projects/Combinatorics/O.vulgaris/Alignment/PRJNA791920/IsoSeq.Polished.Unclustered/AggregatedByChromBySampleSummary.tsv"

# %%
# tmr50_alignment_stats_df = pd.read_csv(alignment_stats_file, sep="\t")
# tmr50_alignment_stats_df = tmr50_alignment_stats_df.loc[
#     tmr50_alignment_stats_df["MappedReads"] >= 50
# ].reset_index(drop=True)
# tmr50_alignment_stats_df

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

positions_data_df["SNPsPositionsFile"] = positions_data_df["Chrom"].apply(
    lambda x: Path(positions_dir, f"{x}.positions.snps.csv.gz")
)
positions_data_df["SNPsReadsFile"] = positions_data_df["Chrom"].apply(
    lambda x: Path(reads_dir, f"{x}.reads.snps.csv.gz")
)

positions_data_df


# %%
# first_x_chroms = 30

# positions_data_df = positions_data_df.iloc[:first_x_chroms]

# positions_data_df

# %% [markdown]
# # Functions

# %% [markdown]
# ## Positions

# %%
def find_alt_base(ref_base, a_count, t_count, c_count, g_count, seed):

    bases = list("ATCG")
    base_counts = [a_count, t_count, c_count, g_count]

    alt_base_counts = {
        base: count for base, count in zip(bases, base_counts) if base != ref_base
    }

    max_alt_base_count = max(alt_base_counts.values())

    max_alt_bases = [
        base for base, count in alt_base_counts.items() if count == max_alt_base_count
    ]

    # randomly choose one of the equally most frequent alt bases
    return np.random.default_rng(seed).choice(max_alt_bases)


# %%
def calc_noise(ref_base_count, alt_base_count):
    try:
        noise = alt_base_count / (alt_base_count + ref_base_count)
    except ZeroDivisionError:
        noise = 0  # if there are no mapped alt bases
    return noise


# %%
# for (
#     positions_file, 
#     snps_positions_file, 
#     snps_reads_file,  
#     sample, 
# ) in positions_data_df.loc[
#     :,
#     ["PositionsFile", "SNPsPositionsFile", "SNPsReadsFile", "Chrom"]
# ].values:
#     ic(positions_file, snps_positions_file, snps_reads_file, sample);
#     break

# %%
# noise_positions_df = pd.read_csv(positions_file, sep="\t")
# noise_positions_df = noise_positions_df.loc[
#     (noise_positions_df["CDS"])
#     & (~noise_positions_df[edited_col])
#     & (noise_positions_df["TotalCoverage"].gt(0))
# ]
# if final_noise_col is not None:
#     noise_positions_df = noise_positions_df.loc[
#         noise_positions_df[final_noise_col]
#     ]
# noise_positions_df = noise_positions_df.drop(
#     columns=[
#         "CDS",
#         edited_col,
#         "InProbRegion",
#         "EditingFrequency",
#         "KnownEditing",
#     ]
# )
# noise_positions_df

# %%
def make_snps_positions(
    positions_file: Path | str,
    snps_positions_file: Path | str,
    top_x_noisy_positions: int,
    assurance_factor: float,
    seed: int,
    edited_col: str,
    snp_noise_level: float,
    final_noise_col: str | None,
    new_old_to_new_reads_file: Path | str | None = None,
    old_old_to_new_reads_file: Path | str | None = None,
):
    noise_positions_df = pd.read_csv(positions_file, sep="\t")
    noise_positions_df = noise_positions_df.loc[
        (noise_positions_df["CDS"])
        & (~noise_positions_df[edited_col])
        & (noise_positions_df["TotalCoverage"].gt(0))
    ]
    if final_noise_col is not None:
        noise_positions_df = noise_positions_df.loc[
            noise_positions_df[final_noise_col]
        ]
    noise_positions_df = noise_positions_df.drop(
        columns=[
            "CDS",
            edited_col,
            "InProbRegion",
            "EditingFrequency",
            "KnownEditing",
        ]
    )
    
    # finalize empty df with expected cols, save it and return
    if noise_positions_df.empty:
        snps_positions_df = noise_positions_df.assign(
            AltBase=None,
            AboveEditingThreshold=None,
        )
        snps_positions_df.to_csv(snps_positions_file, sep="\t", index=False)
        return snps_positions_df
    
    noise_positions_df["AltBase"] = noise_positions_df.apply(
        lambda x: find_alt_base(
            x["RefBase"],
            x["A"],
            x["T"],
            x["C"],
            x["G"],
            seed,
        ),
        axis=1,
    )
    noise_positions_df = noise_positions_df.loc[
        ~(
            (noise_positions_df["RefBase"].eq("A"))
            & (noise_positions_df["AltBase"].eq("G"))
        )
    ]
    noise_positions_df["Noise"] = noise_positions_df.apply(
        lambda x: calc_noise(x[x["RefBase"]], x[x["AltBase"]]), axis=1
    )
    
    # determine editing threshold to annotate which noise positions have high noise and are thus annotated as "strong" SNPs
    if final_noise_col is not None:
        noise_levels = (
            noise_positions_df.loc[
                (noise_positions_df["Noise"] < snp_noise_level)
                & (noise_positions_df[final_noise_col]),
                "Noise",
            ]
            .sort_values(ascending=False)[:top_x_noisy_positions]
            .tolist()
        )
    else:
        noise_levels = (
            noise_positions_df.loc[
                (noise_positions_df["Noise"] < snp_noise_level),
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
    # finalize the editing threshold
    noise_threshold *= assurance_factor
    
    snps_positions_df = noise_positions_df.loc[
        noise_positions_df["Noise"].ge(snp_noise_level)
    ].copy()

    # now annotate the strong SNPs
    snps_positions_df["AboveEditingThreshold"] = snps_positions_df["Noise"].ge(noise_threshold)

    if (old_old_to_new_reads_file is not None) and (new_old_to_new_reads_file is not None):

        # map old and new compressed reads names
        new_old_to_new_reads_df = pd.read_table(new_old_to_new_reads_file).rename(
            columns={"NewRead": "NewCompressedReadName"}
        )
        old_old_to_new_reads_df = pd.read_table(old_old_to_new_reads_file).rename(
            columns={"NewRead": "OldCompressedReadName"}
        )
        new_to_old_shorted_reads_df = new_old_to_new_reads_df.merge(old_old_to_new_reads_df, on="OldRead", how="outer")
        assert not new_to_old_shorted_reads_df.isna().any().any()
        new_to_old_shorted_reads_dict = dict(
            zip(
                new_to_old_shorted_reads_df['NewCompressedReadName'], 
                new_to_old_shorted_reads_df['OldCompressedReadName']
            )
        )
        
        # use the mapping to replace the new compressed reads names with the old ones (the originals)
        snps_positions_df["SplitNewCompressedReads"] = snps_positions_df["Reads"].str.split(",")
        snps_positions_df["SplitOldCompressedReads"] = snps_positions_df["SplitNewCompressedReads"].apply(
            lambda x: [new_to_old_shorted_reads_dict[y] for y in x],
        )
        snps_positions_df["Reads"] = snps_positions_df["SplitOldCompressedReads"].apply(lambda x: ",".join(x))
        del snps_positions_df["SplitNewCompressedReads"]
        del snps_positions_df["SplitOldCompressedReads"]
    
    snps_positions_df.to_csv(snps_positions_file, sep="\t", index=False)

    return snps_positions_df


# %% [markdown]
# ## Reads

# %%
def add_noise_se_reads_in_positions(noise_positions_df, sample):
    unique_reads = set(chain.from_iterable(noise_positions_df["Reads"].str.split(",")))
    positions_per_read = {read: [] for read in unique_reads}

    # alt_base = "G" if strand == "+" else "C"
    for x, position_row in enumerate(noise_positions_df.itertuples()):
        mapped_bases = position_row.MappedBases
        mapped_reads = position_row.Reads.split(",")
        alt_base = position_row.AltBase
        # deal with reads that were mapped to this pos
        for mapped_base, mapped_read in zip(mapped_bases, mapped_reads):
            if mapped_base == alt_base:
                pos_edited_in_read = 1  # mismatch
            elif mapped_base == ".":
                pos_edited_in_read = 0  # match
            else:
                pos_edited_in_read = (
                    -1
                )  # we ignore mismatches that aren't RefBase2AltBase by marking them as NaN
            positions_per_read[mapped_read].append(pos_edited_in_read)
        # deal with reads that weren't mapped to this pos
        unmapped_reads_in_pos = unique_reads - set(mapped_reads)
        for mapped_read in unmapped_reads_in_pos:
            pos_edited_in_read = -1
            positions_per_read[mapped_read].append(pos_edited_in_read)
        # check all reads got values for pos
        mapped_pos_dist = {len(bases) for bases in positions_per_read.values()}
        if len(mapped_pos_dist) != 1:
            if sample is not None:
                raise Exception(
                    f"Problem at line {x} in sample {sample}: not all reads are mapped."
                )
            else:
                raise Exception(f"Problem at line {x}: not all reads are mapped.")

    return positions_per_read


# %%
def add_noise_pe_reads_in_positions(noise_positions_df, sample):
    unique_reads = set(chain.from_iterable(noise_positions_df["Reads"].str.split(",")))
    positions_per_read = {read: [] for read in unique_reads}

    overlapping_bases = 0

    # alt_base = "G" if strand == "+" else "C"
    for x, position_row in enumerate(noise_positions_df.itertuples()):
        mapped_bases = position_row.MappedBases
        mapped_reads = position_row.Reads.split(",")
        phred_scores = position_row.Phred
        alt_base = position_row.AltBase

        last_phred_of_read_in_position = {}

        # deal with reads that were mapped to this pos
        for phred_score, mapped_base, mapped_read in zip(
            phred_scores, mapped_bases, mapped_reads
        ):
            if mapped_base == alt_base:
                pos_edited_in_read = 1  # mismatch
            elif mapped_base == ".":
                pos_edited_in_read = 0  # match
            else:
                pos_edited_in_read = (
                    -1
                )  # we ignore mismatches that aren't RefBase2AltBase by marking them as NaN

            # deal with overlapping paired-ended reads mapped to this position
            if mapped_read not in last_phred_of_read_in_position:
                positions_per_read[mapped_read].append(pos_edited_in_read)
                last_phred_of_read_in_position[mapped_read] = phred_score
            else:
                # for a pair of PE reads, we only keep the one with the highest phred score
                if ord(phred_score) > ord(last_phred_of_read_in_position[mapped_read]):
                    positions_per_read[mapped_read][
                        -1
                    ] = pos_edited_in_read  # overwrite the previous editing status of the read in this position
                    last_phred_of_read_in_position[mapped_read] = phred_score
                overlapping_bases += 1

        # deal with reads that weren't mapped to this pos
        unmapped_reads_in_pos = unique_reads - set(mapped_reads)
        for mapped_read in unmapped_reads_in_pos:
            pos_edited_in_read = -1
            positions_per_read[mapped_read].append(pos_edited_in_read)

        # check all reads got values for pos
        mapped_pos_dist = {len(bases) for bases in positions_per_read.values()}
        if len(mapped_pos_dist) != 1:
            if sample is not None:
                raise Exception(
                    f"Problem at line {x} in sample {sample}: not all reads are mapped."
                )
            else:
                raise Exception(f"Problem at line {x}: not all reads are mapped.")

    # if sample is not None:
    #     ic(sample, overlapping_bases)

    return positions_per_read


# %%
def add_noise_multisample_se_reads_in_positions(noise_positions_df, sample):
    unique_reads = set(chain.from_iterable(noise_positions_df["Reads"].str.split(",")))
    positions_per_read = {read: [] for read in unique_reads}
    sample_per_read = {}

    # alt_base = "G" if strand == "+" else "C"
    for x, position_row in enumerate(noise_positions_df.itertuples()):
        mapped_bases = position_row.MappedBases
        mapped_reads = position_row.Reads.split(",")
        mapped_samples = position_row.Samples.split(",")
        alt_base = position_row.AltBase
        # deal with reads that were mapped to this pos
        for mapped_base, mapped_read, mapped_sample in zip(mapped_bases, mapped_reads, mapped_samples):
            
            sample_per_read[mapped_read] = mapped_sample
            
            if mapped_base == alt_base:
                pos_edited_in_read = 1  # mismatch
            elif mapped_base == ".":
                pos_edited_in_read = 0  # match
            else:
                pos_edited_in_read = (
                    -1
                )  # we ignore mismatches that aren't RefBase2AltBase by marking them as NaN
            positions_per_read[mapped_read].append(pos_edited_in_read)
        # deal with reads that weren't mapped to this pos
        unmapped_reads_in_pos = unique_reads - set(mapped_reads)
        for mapped_read in unmapped_reads_in_pos:
            pos_edited_in_read = -1
            positions_per_read[mapped_read].append(pos_edited_in_read)
        # check all reads got values for pos
        mapped_pos_dist = {len(bases) for bases in positions_per_read.values()}
        if len(mapped_pos_dist) != 1:
            raise Exception(
                    f"Problem at line {x} in sample {sample}: not all reads are mapped."
                )

    return positions_per_read


# %%
def make_snps_reads(
    snps_positions_df,
    snps_reads_file,
    platform, 
    sample, 
    parity, 
    multisample, 
):
    # If there are no SNP positions, there are no reads to encode.
    if snps_positions_df is None or snps_positions_df.empty:
        reads_snps_df = pd.DataFrame(columns=["Platform", "Sample", "Read"])
        reads_snps_df.to_csv(snps_reads_file, index=False, sep="\t")
        return reads_snps_df
    
    # defensive: ensure Reads/Samples exist and are strings
    if "Reads" not in snps_positions_df.columns:
        raise KeyError("snps_positions_df is missing required column 'Reads'")
    if multisample and "Samples" not in snps_positions_df.columns:
        raise KeyError("snps_positions_df is missing required column 'Samples'")
    
    if parity == "SE":
        if multisample:
            positions_per_read = add_noise_multisample_se_reads_in_positions(
                snps_positions_df, sample=sample
            )
        else:
            positions_per_read = add_noise_se_reads_in_positions(
                snps_positions_df, sample=sample
            )
        
    elif parity == "PE":
        if multisample:
            raise NotImplementedError("Multisample PE reads in positions is not implemented yet.")
        else:
            positions_per_read = add_noise_pe_reads_in_positions(
                snps_positions_df, sample=sample
            )
    else:
        raise ValueError("parity must be either 'SE' or 'PE'")

    reads_snps_df = pd.DataFrame(positions_per_read)
    reads_snps_df = reads_snps_df.T.reset_index().rename(
        {"index": "Read"}, axis="columns"
    )
    reads_snps_df = reads_snps_df.rename(
        columns={
            old_col: pos
            for old_col, pos in zip(
                reads_snps_df.columns[1:], snps_positions_df["Position"]
            )
        }
    )
            
    reads_snps_df.insert(0, "Platform", platform)
    reads_snps_df.insert(1, "Sample", sample)
    
    reads_snps_df.to_csv(snps_reads_file, index=False, sep="\t")

    return reads_snps_df


# %% [markdown]
# ## Positions & reads

# %%
def make_positions_and_reads_snps_dfs(
    positions_file,
    snps_positions_file,
    top_x_noisy_positions,
    assurance_factor,
    seed,
    edited_col,
    snp_noise_level,
    final_noise_col,
    snps_reads_file,
    platform,
    sample,
    parity,
    multisample,
    new_old_to_new_reads_file=None,
    old_old_to_new_reads_file=None,
):
    snps_positions_df = make_snps_positions(
        positions_file, 
        snps_positions_file, 
        top_x_noisy_positions,
        assurance_factor,
        seed, 
        edited_col, 
        snp_noise_level, 
        final_noise_col,
        new_old_to_new_reads_file=new_old_to_new_reads_file,
        old_old_to_new_reads_file=old_old_to_new_reads_file
    )
    snps_reads_df = make_snps_reads(
        snps_positions_df, snps_reads_file, platform, sample, parity, multisample
    )
    return snps_positions_df, snps_reads_df


# %% [markdown]
# # Finding SNPs

# %%
# # ?make_positions_and_reads_snps_dfs

# %%
with Pool(processes=PROCESSES) as pool:
    positions_and_reads_snps_dfs = pool.starmap(
        func=make_positions_and_reads_snps_dfs,
        iterable=[
            (
                positions_file,
                snps_positions_file,
                top_x_noisy_positions,
                assurance_factor,
                seed,
                edited_col,
                snp_noise_level,
                final_noise_col,
                snps_reads_file,
                platform,
                sample,
                parity,
                multisample,
                None, # new_old_to_new_reads_file
                None  # old_old_to_new_reads_file
            )
            for (
                positions_file, 
                snps_positions_file, 
                snps_reads_file,  
                sample, 
            ) in positions_data_df.loc[
                :,
                ["PositionsFile", "SNPsPositionsFile", "SNPsReadsFile", "Chrom"]
            ].values

        ],
    )

positions_snps_dfs = [dfs[0] for dfs in positions_and_reads_snps_dfs]
reads_snps_dfs = [dfs[1] for dfs in positions_and_reads_snps_dfs]

# %%
# non_empty_chroms = []
# for chrom, positions_snps_df, reads_snps_df in zip(
#     positions_data_df["Chrom"].values,
#     positions_snps_dfs,
#     reads_snps_dfs
# ):
#     if not positions_snps_df.empty:
#         non_empty_chroms.append(chrom)
#     # ic(
#     #     chrom,
#     #     positions_snps_df.shape[0],
#     #     reads_snps_df.shape[0]
#     # );
# non_empty_chroms

# %%
# for chrom, positions_snps_df, reads_snps_df in zip(
#     positions_data_df["Chrom"].values,
#     positions_snps_dfs,
#     reads_snps_dfs
# ):
#     if not positions_snps_df.empty:
#         # non_empty_chroms.append(chrom)
#         ic(
#             chrom,
#             positions_snps_df.shape[0],
#             reads_snps_df.shape[0]
#         );
# # non_empty_chroms

# %%
# positions_snps_df

# %%
# reads_snps_df
