from itertools import chain

# from pathlib import Path
# import sys
from multiprocessing import Pool

from icecream import ic
import plotly.express as px

# import plotly.graph_objects as go
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact, chi2_contingency
from statsmodels.stats.multitest import fdrcorrection


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


def calc_noise(ref_base_count, alt_base_count):
    try:
        noise = alt_base_count / (alt_base_count + ref_base_count)
    except ZeroDivisionError:
        noise = 0  # if there are no mapped alt bases
    return noise


def make_noise_positions_df(positions_file, seed, edited_col="Edited"):
    noise_positions_df = pd.read_csv(positions_file, sep="\t")
    noise_positions_df = noise_positions_df.loc[
        (noise_positions_df["CDS"])
        # & (noise_positions_df["Noise"].notna())
        & (~noise_positions_df[edited_col])
        & (noise_positions_df["TotalCoverage"].gt(0))
    ].drop(
        columns=[
            "CDS",
            edited_col,
            "InProbRegion",
            "EditingFrequency",
            # "Phred",
            "KnownEditing",
        ]
    )
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

    return noise_positions_df


def add_noise_se_reads_in_positions(noise_positions_df, sample=None):
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


def add_noise_pe_reads_in_positions(noise_positions_df, sample=None):
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


def make_reads_noise_df(
    positions_file,
    seed,
    edited_col,
    platform,
    sample,
    parity,
    old_to_new_reads_file=None,
):
    noise_positions_df = make_noise_positions_df(
        positions_file, seed, edited_col=edited_col
    )
    if parity == "SE":
        positions_per_read = add_noise_se_reads_in_positions(
            noise_positions_df, sample=sample
        )
    elif parity == "PE":
        positions_per_read = add_noise_pe_reads_in_positions(
            noise_positions_df, sample=sample
        )
    else:
        raise ValueError("parity must be either 'SE' or 'PE'")

    reads_noise_df = pd.DataFrame(positions_per_read)
    reads_noise_df = reads_noise_df.T.reset_index().rename(
        {"index": "Read"}, axis="columns"
    )
    reads_noise_df = reads_noise_df.rename(
        columns={
            old_col: pos
            for old_col, pos in zip(
                reads_noise_df.columns[1:], noise_positions_df["Position"]
            )
        }
    )

    if old_to_new_reads_file is not None:
        old_to_new_reads_df = pd.read_table(old_to_new_reads_file)
        reads_noise_df = old_to_new_reads_df.merge(
            reads_noise_df, left_on="NewRead", right_on="Read"
        ).drop(
            # retain "OldRead" and "Read" (which is the same as "NewRead")
            columns="NewRead"
        )

    reads_noise_df.insert(0, "Platform", platform)
    reads_noise_df.insert(1, "Sample", sample)

    return reads_noise_df


seed = 1892
edited_col = "Edited"

old_reads_first_pos_loc = 6

# if noise positions were re-created where non ref base positions with noise >= 0.1 weren't discarded.
# in that case - the new reads files may have different new reads names,
# so to connect between the reads with editing statuses and reads with new noise statuses,
# we need to use to old reads col which is another col in the new reads df.
new_noise_positions_in_use = True


pacbio_samples = ["GRIA2", "PCLO", "ADAR1", "IQEC1"]

pacbio_new_positions_files = [
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/CompleteNoisePositions/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.positions.csv.gz",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/CompleteNoisePositions/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.positions.csv.gz",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/CompleteNoisePositions/ADAR1.Merged.r64296e203404D01.aligned.sorted.MinRQ998.positions.csv.gz",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/CompleteNoisePositions/IQEC.Merged.r64296e203404D01.aligned.sorted.MinRQ998.positions.csv.gz",
]
pacbio_new_old_to_new_reads_files = [
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/CompleteNoisePositions/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.OldToNewReads.csv.gz",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/CompleteNoisePositions/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.OldToNewReads.csv.gz",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/CompleteNoisePositions/ADAR1.Merged.r64296e203404D01.aligned.sorted.MinRQ998.OldToNewReads.csv.gz",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/CompleteNoisePositions/IQEC.Merged.r64296e203404D01.aligned.sorted.MinRQ998.OldToNewReads.csv.gz",
]

pacbio_old_reads_files = [
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.reads.csv.gz",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.reads.csv.gz",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples/ADAR1.Merged.r64296e203404D01.aligned.sorted.MinRQ998.reads.csv.gz",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples/IQEC.Merged.r64296e203404D01.aligned.sorted.MinRQ998.reads.csv.gz",
]
pacbio_old_unique_reads_files = [
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_reads.csv.gz",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_reads.csv.gz",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples/ADAR1.Merged.r64296e203404D01.aligned.sorted.MinRQ998.unique_reads.csv.gz",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples/IQEC.Merged.r64296e203404D01.aligned.sorted.MinRQ998.unique_reads.csv.gz",
]
pacbio_old_old_to_new_reads_files = [
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.OldToNewReads.csv.gz",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.OldToNewReads.csv.gz",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples/ADAR1.Merged.r64296e203404D01.aligned.sorted.MinRQ998.OldToNewReads.csv.gz",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples/IQEC.Merged.r64296e203404D01.aligned.sorted.MinRQ998.OldToNewReads.csv.gz",
]

illumina_samples = [
    "RUSC2",
    "TRIM2",
    "CA2D3",  # here's the problem, line 43
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

illumina_new_positions_files = [
    f"/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/CompleteNoisePositions/reads.sorted.aligned.filtered.{chrom}.positions.csv.gz"
    for chrom in illumina_chroms
]
illumina_new_old_to_new_reads_files = [
    f"/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/CompleteNoisePositions/reads.sorted.aligned.filtered.{chrom}.OldToNewReads.csv.gz"
    for chrom in illumina_chroms
]

illumina_old_reads_files = [
    f"/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.{chrom}.reads.csv"
    for chrom in illumina_chroms
]
illumina_old_unique_reads_files = [
    f"/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.{chrom}.unique_reads.csv"
    for chrom in illumina_chroms
]
illumina_old_old_to_new_reads_files = [
    f"/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.{chrom}.OldToNewReads.csv"
    for chrom in illumina_chroms
]

platforms = ["PacBio"] * len(pacbio_samples) + ["Illumina"] * len(illumina_samples)
parities = ["SE"] * len(pacbio_samples) + ["PE"] * len(illumina_samples)
samples = pacbio_samples + illumina_samples
# new positions files made for complete noise estimates,
# and their corresponding old_to_new_reads_files
new_positions_files = pacbio_new_positions_files + illumina_new_positions_files
new_old_to_new_reads_files = (
    pacbio_new_old_to_new_reads_files + illumina_new_old_to_new_reads_files
)
# original reads - used for estimating editing
old_reads_files = pacbio_old_reads_files + illumina_old_reads_files
old_unique_reads_files = pacbio_old_unique_reads_files + illumina_old_unique_reads_files
old_old_to_new_reads_files = (
    pacbio_old_old_to_new_reads_files + illumina_old_old_to_new_reads_files
)

with Pool(processes=len(samples)) as pool:
    reads_noise_dfs = pool.starmap(
        func=make_reads_noise_df,
        iterable=[
            (
                positions_file,
                seed,
                edited_col,
                platform,
                sample,
                parity,
                old_to_new_reads_file,
            )
            for positions_file, platform, sample, parity, old_to_new_reads_file in zip(
                new_positions_files,
                platforms,
                samples,
                parities,
                new_old_to_new_reads_files,
            )
        ],
    )

reads_noise_first_pos_loc = 3
if new_noise_positions_in_use:
    reads_noise_first_pos_loc += 1