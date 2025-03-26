from collections import defaultdict
import subprocess
from pathlib import Path
from typing import Union
from itertools import chain
from multiprocessing import Pool

from pybedtools import BedTool
import numpy as np
import pandas as pd

from icecream import ic
from statsmodels.stats.proportion import binom_test
from statsmodels.stats.multitest import fdrcorrection
from Bio import Seq

from Pileup.ids import compress_ids
from EditingUtils.seq import make_fasta_dict


def replace_reads_names_in_pos(old_reads_str, old_to_new_mapping, reads_sep):
    old_reads = old_reads_str.split(reads_sep)
    new_reads = [old_to_new_mapping[old_read] for old_read in old_reads]
    new_reads_str = reads_sep.join(new_reads)
    return new_reads_str


def replace_reads_names(
    positions_df,
    total_coverage_col,
    reads_col,
    reads_sep=",",
    mapping_out_file=None,
    out_files_sep="\t",
):
    if len(positions_df) == 0:
        return
    # replace
    old_unique_reads = set(
        chain.from_iterable(positions_df[reads_col].str.split(reads_sep))
    )
    # # the previous line is equivalent to the following, but the previous is sometimes killed by the kernel
    # old_unique_reads = set()
    # l = len(positions_df)
    # window = 100
    # for i in range(0, l, window):
    #     j = min(i + window, l)
    #     old_unique_reads_in_window = set(
    #         [
    #             read
    #             for reads_in_position in positions_df.iloc[i:j][reads_col]
    #             for read in set(reads_in_position.split(reads_sep))
    #         ]
    #     )
    #     old_unique_reads = old_unique_reads | old_unique_reads_in_window
    old_to_new_mapping = compress_ids(old_unique_reads)
    positions_df[reads_col] = positions_df[reads_col].apply(
        lambda old_reads_str: replace_reads_names_in_pos(
            old_reads_str, old_to_new_mapping, reads_sep
        )
    )

    # verify correct replacement of substrings of old to new reads, `both` seperated by `reads_sep`
    # tests 1 - proof that `reads_sep` is used in `reads_col` only to sepearate reads, i.e., reads don't have `reads_sep` in their names
    # if positions_df[total_coverage_col].ne(positions_df[reads_col].str.count(reads_sep) + 1).any():
    if (
        ~positions_df[total_coverage_col]
        .eq(positions_df[reads_col].str.count(reads_sep) + 1)
        .all()
    ):
        raise Exception()  # todo add text to exception
    # test 2
    new_unique_reads_from_df = set(
        chain.from_iterable(positions_df[reads_col].str.split(reads_sep))
    )
    new_unique_reads_from_mapping = set(old_to_new_mapping.values())
    bad_lengths_cond = (
        not len(old_unique_reads)
        == len(new_unique_reads_from_df)
        == len(new_unique_reads_from_mapping)
    )
    bad_replacing_cond = new_unique_reads_from_df != new_unique_reads_from_mapping
    if bad_lengths_cond or bad_replacing_cond:
        raise Exception()  # todo add text to exception

    # save mapping to a csv file
    if mapping_out_file is not None:
        old_keys = []
        new_keys = []
        for old_key, new_key in old_to_new_mapping.items():
            old_keys.append(old_key)
            new_keys.append(new_key)
        mapping_df = pd.DataFrame({"OldRead": old_keys, "NewRead": new_keys})
        mapping_df.to_csv(mapping_out_file, sep=out_files_sep, index=False)


def positions_df_to_bed(positions_df: pd.DataFrame, strand: str) -> BedTool:
    """Convert positions DataFrame to a 6-columns BedTool.

    Args:
        positions_df (pd.DataFrame): A DataFrame of positions created from the pileup file.
        strand (str): The ORF's strand.

    Returns:
        BedTool: A BedTool of positions.
    """
    positions_bed = BedTool.from_dataframe(
        positions_df.loc[:, ["Chrom", "Position"]].assign(
            # EndPosition=positions_df.iloc[:, 1] + 1,
            EndPosition=positions_df["Position"] + 1,
            Name=".",
            Score=".",
            Strand=strand,
        )
    )
    return positions_bed


def annotate_prolamatic_sites(positions_df, strand, problamatic_regions_file=None):
    if len(positions_df) > 0:
        if problamatic_regions_file:
            problamatic_regions_bed = BedTool(problamatic_regions_file)
            positions_bed = positions_df_to_bed(positions_df, strand)
            prob_positions_bed = positions_bed.intersect(problamatic_regions_bed)
            prob_positions = {interval.start for interval in prob_positions_bed}

            def is_prob(position: int) -> bool:
                return position in prob_positions

            in_prob_region = [
                is_prob(position) for position in positions_df["Position"]
            ]
        else:
            in_prob_region = [False for _ in range(len(positions_df))]
        positions_df.insert(
            positions_df.columns.get_loc("Position") + 1, "InProbRegion", in_prob_region
        )
    else:
        # add an empty column to the empty dataframe for comapatibility with the following steps
        positions_df.insert(
            positions_df.columns.get_loc("Position") + 1, "InProbRegion", None
        )


def annotate_known_sites(positions_df, strand, known_sites_file=None):
    if len(positions_df) > 0:
        if known_sites_file:
            known_sites_bed = BedTool(known_sites_file)
            positions_bed = positions_df_to_bed(positions_df, strand)
            known_positions_bed = positions_bed.intersect(known_sites_bed, s=True)
            known_positions = {interval.start for interval in known_positions_bed}

            def is_known(position: int) -> bool:
                return position in known_positions

            is_position_known = [
                is_known(position) for position in positions_df["Position"]
            ]
        else:
            is_position_known = [False for _ in range(len(positions_df))]
        positions_df.insert(
            positions_df.columns.get_loc("Position") + 1,
            "KnownEditing",
            is_position_known,
        )
    else:
        # add an empty column to the empty dataframe for comapatibility with the following steps
        positions_df.insert(
            positions_df.columns.get_loc("Position") + 1, "KnownEditing", None
        )


def annotate_coding_sites(positions_df, strand, cds_regions_file=None):
    if len(positions_df) > 0:
        if cds_regions_file:
            cds_regions_bed = BedTool(cds_regions_file)
            positions_bed = positions_df_to_bed(positions_df, strand)
            coding_positions_bed = positions_bed.intersect(cds_regions_bed, s=True)
            coding_positions = {interval.start for interval in coding_positions_bed}

            def is_coding(position: int) -> bool:
                return position in coding_positions

            is_position_coding = [
                is_coding(position) for position in positions_df["Position"]
            ]
        else:
            is_position_coding = [True for _ in range(len(positions_df))]
        positions_df.insert(
            positions_df.columns.get_loc("Position") + 1, "CDS", is_position_coding
        )
    else:
        # add an empty column to the empty dataframe for comapatibility with the following steps
        positions_df.insert(positions_df.columns.get_loc("Position") + 1, "CDS", None)


def count_bases(ref_base, mapped_bases):
    mapped_bases = mapped_bases.replace(".", ref_base)
    return (
        mapped_bases.count("A"),
        mapped_bases.count("T"),
        mapped_bases.count("C"),
        mapped_bases.count("G"),
    )


def annotate_base_counts(positions_df):
    if len(positions_df) > 0:
        atcgs = positions_df.apply(
            lambda x: count_bases(x["RefBase"], x["MappedBases"]),
            axis=1,
            result_type="expand",
        ).rename(columns={0: "A", 1: "T", 2: "C", 3: "G"})
        positions_df[["A", "T", "C", "G"]] = atcgs
    else:
        # add empty columns to the empty dataframe for comapatibility with the following steps
        positions_df[["A", "T", "C", "G"]] = None


def position_noise_level(
    ref_base: str, strand: str, a_count: int, t_count: int, c_count: int, g_count: int
) -> float:
    """
    Calculate noise level in a position.

    The noise is defined as number of bases of the most abundant `alt_base`,
    divided by the same number + the number of `ref_base` bases.

    Args:
        `ref_base` (str): The reference base of the position.
        `strand` (str): The ORF's strand.
        `a_count` (int): Number of As mapped to position.
        `t_count` (int): Number of Ts mapped to position.
        `c_count` (int): Number of Cs mapped to position.
        `g_count` (int): Number of Gs mapped to position.

    Returns:
        float: The noise level.
    """
    base_counts = [a_count, t_count, c_count, g_count]

    if sum(base_counts) == 0:
        noise = np.nan

    # we only measure noise in positions that don't undergo RNA editing by ADAR
    # todo: in futurue versions, do consider adenosine positions if the most abundant base is not G
    elif (strand == "+" and ref_base == "A") or (strand == "-" and ref_base == "T"):
        noise = np.nan

    else:
        # we measure noise for T positions only on the positive strand
        if strand == "+" and ref_base == "T":
            ref_base_count = t_count
        # we measure noise for A positions only on the negative strand
        elif strand == "-" and ref_base == "A":
            ref_base_count = a_count
        # we measure noise for A & G positions on both strand
        elif ref_base == "C":
            ref_base_count = c_count
        else:  # ref_base == "G"
            ref_base_count = g_count

        bases = ["A", "T", "C", "G"]
        alt_bases = set(bases) - {ref_base}
        alt_base_counts = [
            base_count
            for base_count, base in zip(base_counts, bases)
            if base in alt_bases
        ]
        max_alt_base_count = max(alt_base_counts)

        try:
            noise = max_alt_base_count / (max_alt_base_count + ref_base_count)
        except ZeroDivisionError:
            noise = 0  # if there are no mapped alt bases

    return noise


def annotate_noise(positions_df, strand):
    if len(positions_df) > 0:
        positions_df["Noise"] = positions_df.apply(
            lambda x: position_noise_level(
                x["RefBase"], strand, x["A"], x["T"], x["C"], x["G"]
            ),
            axis=1,
        )
    else:
        # add an empty column to the empty dataframe for comapatibility with the following steps
        positions_df["Noise"] = None


# def editing_frequency_per_position(
#     ref_base: str, base: str, ref_base_count: int, alt_base_count: int
# ) -> float:
#     if base == ref_base:
#         try:
#             freq = alt_base_count / (ref_base_count + alt_base_count)
#         except ZeroDivisionError:
#             freq = 0
#     else:
#         freq = np.nan
#     return freq


def editing_frequency_per_position(
    principal_ref_base: str, ref_base: str, ref_base_count: int, alt_base_count: int
) -> float:
    if ref_base == principal_ref_base:
        try:
            freq = alt_base_count / (ref_base_count + alt_base_count)
        except ZeroDivisionError:
            freq = 0
    else:
        freq = np.nan
    return freq


def annotate_editing_frequency_per_position(positions_df: pd.DataFrame, strand: str):
    principal_ref_base = "A" if strand == "+" else "T"
    principal_alt_base = "G" if strand == "+" else "C"

    positions_df.insert(
        len(positions_df.columns),
        "EditingFrequency",
        positions_df.apply(
            lambda x: editing_frequency_per_position(
                principal_ref_base,
                x["RefBase"],
                x[principal_ref_base],
                x[principal_alt_base],
            ),
            axis=1,
        ),
    )


def annotate_edited_sites(
    positions_df: pd.DataFrame,
    strand: str,
    noise_threshold: float,
    denovo_detection: bool = True,
):
    """Determine which sites are currently edited.

    For each position, do this by checking if editing_frequency > noise_threshold.

    Args:
        positions_df (pd.DataFrame): The positions dataframe.
        strand (str): The ORF's strand.
        noise_threshold (float): A site is considered edited if its editing frequency is above the noise threshold.
        denovo_detection (bool): Whether to determine both new ("denovo") and known sites as edited, or only known ones. Defaults to True.
    """
    ref_base = "A" if strand == "+" else "T"

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
    positions_df.insert(positions_df.columns.get_loc("Position") + 1, "Edited", edited)


def pileup_to_positions(
    pileup_file: Path,
    strand: str,
    min_percent_of_max_coverage: float,
    snp_noise_level: float,
    top_x_noisy_positions: int,
    assurance_factor: float,
    problamatic_regions_file: Union[Path, str, None] = None,
    known_sites_file: Union[Path, str, None] = None,
    cds_regions_file: Union[Path, str, None] = None,
    positions_out_file: Union[Path, str, None] = None,
    reads_mapping_file: Union[Path, str, None] = None,
    out_files_sep: str = "\t",
    keep_pileup_file: bool = False,
    remove_non_refbase_noisy_positions: bool = True,
    denovo_detection: bool = True,
) -> pd.DataFrame:
    """Read pileup file into a DataFrame.
    Verify format and change to 0-based coordinates.
    Possibly remove unedited positions, and possibly mask positions in problamatic regions.

    Args:
        pileup_file (Path): Path to a file creates by the `mpileup` function.
        strand (str): The ORF's strand.
        min_percent_of_max_coverage (float): Keep only positions with coverage >= min_percent_of_max_coverage * max_coverage.
        snp_noise_level (float): Treat non-refbase positions with noise >= snp_noise_level as SNPs, so they won't be considered for noise threshold.
        top_x_noisy_positions (int): Use this many positions with highest noise levels to define initial noise threshold.
        assurance_factor (float): Multiply the measured noise by this factor to define the actual noise threshold.
        problamatic_regions_file (Union[Path, str, None]): If given, annotate positions as residing in problamatic genomic regions.
        known_sites_file (Union[Path, str, None]): If given, annotate positions as known editing sites.
        cds_regions_file (Union[Path, str, None]): If given, annotate positions as coding sites. Else, *all sites are considered as coding!*.
        positions_out_file (Union[Path, str, None]): If given, write positions_df to this file.
        reads_mapping_file (Union[Path, str, None]): If given, write mapping from the original reads to their shortened versions to this file.
        out_files_sep (str): Use this char as a seperator for csv out files. Defaults to tab.
        keep_pileup_file(bool): If True, keep the pileup file from which the positions.csv file is made.
        remove_non_refbase_noisy_positions (bool): If true, remove non refbase positions with too-high noise. Even if this is False, these positions are not considered for noise_threshold determination.
        denovo_detection (bool): Whether to determine both new ("denovo") and known sites as edited, or only known ones. Defaults to True.

    Raises:
        ValueError: If the 5th col doens't contain only base match/ bash mismatch characters, and/or the number of
        reads' names in the 7th col isn't equal to the number of characters in the 5th col.

    Returns:
        pd.DataFrame: A df of the parsed pileup file.
    """
    cols = [
        "Chrom",
        "Position",
        "RefBase",
        "TotalCoverage",
        "MappedBases",
        "Phred",
        "Reads",
    ]
    positions_df = pd.read_csv(pileup_file, sep="\t", names=cols)

    # filter out zero-coverage positions (probably due to deletions) - they are irrelevent, and also disturb `replace_reads_names` later
    positions_df = positions_df.loc[positions_df["TotalCoverage"] > 0]

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

    # verify that the number of mapped base (5th col) corresponds
    # to the number of reads' names (7th col)
    for row_num, row in enumerate(positions_df.itertuples()):
        bases = row.MappedBases
        reads = row.Reads.split(",")
        if len(bases) != len(reads):
            raise ValueError(
                f"Line {row_num} in {pileup_file} contains indels and/or doesn't have the "
                "same number of mapped bases and reads."
            )

    # replace reads' names with a shortened, memory-efficient version
    ic(pileup_file, "before replacing reads' names in pileup_to_positions")
    replace_reads_names(
        positions_df,
        "TotalCoverage",
        "Reads",
        reads_sep=",",
        mapping_out_file=reads_mapping_file,
        out_files_sep=out_files_sep,
    )
    ic(pileup_file, "after replacing reads' names in pileup_to_positions")

    # remove positions with insufficient coverage
    max_coverage = positions_df["TotalCoverage"].max()
    required_coverage = max_coverage * min_percent_of_max_coverage
    positions_df = positions_df.loc[positions_df["TotalCoverage"] >= required_coverage]

    annotate_prolamatic_sites(positions_df, strand, problamatic_regions_file)
    annotate_known_sites(positions_df, strand, known_sites_file)
    annotate_coding_sites(positions_df, strand, cds_regions_file)
    annotate_base_counts(positions_df)
    annotate_noise(positions_df, strand)

    ref_base = "A" if strand == "+" else "T"

    if remove_non_refbase_noisy_positions:
        # remove non refbase positions with too-high noise
        positions_before = len(positions_df)
        positions_df = positions_df.loc[
            (positions_df["RefBase"] == ref_base)
            | (positions_df["Noise"] < snp_noise_level)
        ]
        positions_after = len(positions_df)
        removed_positions = positions_before - positions_after
        print(f"{removed_positions} extremely noisy positions removed")

        # same as `positions_df["Noise"].max()` if `top_x_noise_samples == 1`
        # noise_threshold = positions_df.loc[positions_df["RefBase"] != "A", "Noise"].sort_values(ascending=False)[:top_x_noisy_positions].mean()
        noise_threshold = (
            positions_df["Noise"]
            .sort_values(ascending=False)[:top_x_noisy_positions]
            .mean()
        )
    else:
        noise_threshold = (
            positions_df.loc[positions_df["Noise"] < snp_noise_level, "Noise"]
            .sort_values(ascending=False)[:top_x_noisy_positions]
            .mean()
        )
    noise_threshold *= assurance_factor
    annotate_editing_frequency_per_position(positions_df, strand)
    annotate_edited_sites(positions_df, strand, noise_threshold, denovo_detection)

    # verify that noise is only applied to non ref_base positions, and vice versa for editing frequency
    non_ref_base_edit_freq = positions_df.loc[
        positions_df["RefBase"] != ref_base, "EditingFrequency"
    ].unique()
    assert len(non_ref_base_edit_freq) == 1 and np.isnan(non_ref_base_edit_freq[0])
    ref_base_noise = positions_df.loc[
        positions_df["RefBase"] == ref_base, "Noise"
    ].unique()
    assert len(ref_base_noise) == 1 and np.isnan(ref_base_noise[0])

    if positions_out_file:
        positions_df.to_csv(
            positions_out_file, sep=out_files_sep, index=False, na_rep=np.nan
        )

    if not keep_pileup_file:
        subprocess.run(f"rm {pileup_file}", shell=True)

    # return positions_df


# pileup_to_positions(*pileup_to_positions_inputs)


def multisample_pileups_to_positions_old(
    pileup_files: list[Path],
    samples: list[str],
    strand: str,
    min_percent_of_max_coverage: float,
    snp_noise_level: float,
    top_x_noisy_positions: int,
    assurance_factor: float,
    problamatic_regions_file: Union[Path, str, None] = None,
    known_sites_file: Union[Path, str, None] = None,
    cds_regions_file: Union[Path, str, None] = None,
    positions_out_file: Union[Path, str, None] = None,
    reads_mapping_file: Union[Path, str, None] = None,
    out_files_sep: str = "\t",
    keep_pileup_file: bool = False,
    remove_non_refbase_noisy_positions: bool = True,
    denovo_detection: bool = True,
):
    """Read pileup file into a DataFrame.
    Verify format and change to 0-based coordinates.
    Possibly remove unedited positions, and possibly mask positions in problamatic regions.

    Args:
        `pileup_files` (Path): A list of paths to files creates by the `mpileup` function, each belonging to a different sample.
        `samples`: list[str]: A list of samples' names corresponding to `pileup_files`.
        `strand` (str): The ORF's strand.
        `min_percent_of_max_coverage` (float): Keep only positions with coverage >= min_percent_of_max_coverage * max_coverage.
        `snp_noise_level` (float): Treat non-refbase positions with noise >= snp_noise_level as SNPs, so they won't be considered for noise threshold.
        `top_x_noisy_positions` (int): Use this many positions with highest noise levels to define initial noise threshold.
        `assurance_factor` (float): Multiply the measured noise by this factor to define the actual noise threshold.
        `problamatic_regions_file` (Union[Path, str, None]): If given, annotate positions as residing in problamatic genomic regions.
        `known_sites_file` (Union[Path, str, None]): If given, annotate positions as known editing sites.
        `cds_regions_file` (Union[Path, str, None]): If given, annotate positions as coding sites. Else, *all sites are considered as coding!*.
        `positions_out_file` (Union[Path, str, None]): If given, write positions_df to this file.
        `reads_mapping_file` (Union[Path, str, None]): If given, write mapping from the original reads to their shortened versions to this file.
        `out_files_sep` (str): Use this char as a seperator for csv out files. Defaults to tab.
        `keep_pileup_file`(bool): If True, keep the pileup file from which the positions.csv file is made.
        `remove_non_refbase_noisy_positions` (bool): If true, remove non refbase positions with too-high noise. Even if this is False, these positions are not considered for noise_threshold determination.
        `denovo_detection` (bool): Whether to determine both new ("denovo") and known sites as edited, or only known ones. Defaults to True.

    Raises:
        ValueError: If the 5th col doens't contain only base match/ bash mismatch characters, and/or the number of
        reads' names in the 7th col isn't equal to the number of characters in the 5th col.

    Returns:
        pd.DataFrame: A df of the parsed pileup file.
    """
    cols = [
        "Chrom",
        "Position",
        "RefBase",
        "TotalCoverage",
        "MappedBases",
        "Phred",
        "Reads",
    ]
    # ic()
    # ic(pileup_files)
    # ic(samples)

    positions_dfs = [
        pd.read_csv(pileup_file, sep="\t", names=cols) for pileup_file in pileup_files
    ]
    for positions_df, sample in zip(positions_dfs, samples):
        positions_df.insert(0, "Sample", sample)

    # merge the separate pileups tables into a multi-sample one
    positions_df = pd.concat(positions_dfs, ignore_index=True)

    # filter out zero-coverage positions (probably due to deletions) - they are irrelevent, and also disturb `replace_reads_names` later
    positions_df = positions_df.loc[positions_df["TotalCoverage"] > 0]
    if len(positions_df) == 0:
        return

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
    replace_reads_names(
        positions_df,
        "TotalCoverage",
        "Reads",
        reads_sep=",",
        mapping_out_file=reads_mapping_file,
        out_files_sep=out_files_sep,
    )

    positions_df["Sample"] = positions_df["Sample"].astype(str)
    positions_df["MappedBases"] = positions_df["MappedBases"].astype(str)
    positions_df["Phred"] = positions_df["Phred"].astype(str)
    positions_df["Reads"] = positions_df["Reads"].astype(str)

    positions_df.insert(
        positions_df.columns.get_loc("Sample") + 1,
        "Samples",
        positions_df.apply(
            lambda x: ",".join([x["Sample"] for _ in range(len(x["MappedBases"]))]),
            axis=1,
        ),
    )

    # join multi-sample data per position
    try:
        positions_df = (
            positions_df.groupby(["Chrom", "Position", "RefBase"])
            .agg(
                {
                    "Samples": ",".join,
                    "TotalCoverage": sum,
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

    # remove positions with insufficient coverage
    max_coverage = positions_df["TotalCoverage"].max()
    required_coverage = max_coverage * min_percent_of_max_coverage
    positions_df = positions_df.loc[positions_df["TotalCoverage"] >= required_coverage]

    annotate_prolamatic_sites(positions_df, strand, problamatic_regions_file)
    annotate_known_sites(positions_df, strand, known_sites_file)
    annotate_coding_sites(positions_df, strand, cds_regions_file)
    annotate_base_counts(positions_df)
    annotate_noise(positions_df, strand)

    # positions_df.loc[~positions_df["KnownEditing"], "Noise"].describe()

    ref_base = "A" if strand == "+" else "T"

    if remove_non_refbase_noisy_positions:
        # remove non refbase positions with too-high noise
        positions_before = len(positions_df)
        positions_df = positions_df.loc[
            (positions_df["RefBase"] == ref_base)
            | (positions_df["Noise"] < snp_noise_level)
        ]
        positions_after = len(positions_df)
        removed_positions = positions_before - positions_after
        print(f"{removed_positions} extremely noisy positions removed")

        # same as `positions_df["Noise"].max()` if `top_x_noise_samples == 1`
        # noise_threshold = positions_df.loc[positions_df["RefBase"] != "A", "Noise"].sort_values(ascending=False)[:top_x_noisy_positions].mean()
        noise_threshold = (
            positions_df["Noise"]
            .sort_values(ascending=False)[:top_x_noisy_positions]
            .mean()
        )
    else:
        noise_threshold = (
            positions_df.loc[positions_df["Noise"] < snp_noise_level, "Noise"]
            .sort_values(ascending=False)[:top_x_noisy_positions]
            .mean()
        )
    if pd.isna(noise_threshold):
        noise_threshold = 0
    noise_threshold *= assurance_factor
    annotate_editing_frequency_per_position(positions_df, strand)
    annotate_edited_sites(positions_df, strand, noise_threshold, denovo_detection)

    # verify that noise is only applied to non ref_base positions, and vice versa for editing frequency
    non_ref_base_edit_freq = positions_df.loc[
        positions_df["RefBase"] != ref_base, "EditingFrequency"
    ].unique()
    assert len(non_ref_base_edit_freq) == 1 and np.isnan(non_ref_base_edit_freq[0])
    ref_base_noise = positions_df.loc[
        positions_df["RefBase"] == ref_base, "Noise"
    ].unique()
    assert len(ref_base_noise) == 1 and np.isnan(ref_base_noise[0])

    # verify that the number of mapped bases ==
    # number of reads' names ==
    # number of samples
    for row_num, row in enumerate(positions_df.itertuples()):
        mapped_bases = row.MappedBases
        mapped_reads = row.Reads.split(",")
        mapped_samples = row.Samples.split(",")
        if not (len(mapped_bases) == len(mapped_reads) == len(mapped_samples)):
            raise ValueError

    if positions_out_file:
        positions_df.to_csv(
            positions_out_file, sep=out_files_sep, index=False, na_rep=np.nan
        )

    if not keep_pileup_file:
        for pileup_file in pileup_files:
            subprocess.run(f"rm {pileup_file}", shell=True)

    # return positions_df


def multisample_pileups_to_positions_all_transcripts(
    processes: int,
    chroms: list[str],
    strands: list[str],
    pileup_files_per_chroms: defaultdict[str, list[Path]],
    samples_per_chroms: defaultdict[str, list[str]],
    reads_mapping_files: list[Path],
    positions_files: list[Path],
    min_percent_of_max_coverage: float,
    min_absolute_coverage: Union[int, None],
    snp_noise_level: float,
    top_x_noisy_positions: int,
    assurance_factor: float,
    transcriptome_file: Path,
    final_editing_scheme: str,
    problamatic_regions_file: Union[Path, str, None] = None,
    known_sites_file: Union[Path, str, None] = None,
    cds_regions_file: Union[Path, str, None] = None,
    out_files_sep: str = "\t",
    keep_pileup_files: bool = False,
    remove_non_refbase_noisy_positions: bool = True,
    denovo_detection: bool = True,
    alternative_hypothesis: str = "greater",
    binom_noise_pval_col: str = "NoiseBinomPVal",
    bh_noise_pval_col: str = "NoiseCorrectedPVal",
    bh_noisy_col: str = "NoisyCorrected",
    binom_editing_pval_col: str = "EditingBinomPVal",
    bh_editing_pval_col: str = "EditingCorrectedPVal",
    bh_editing_col: str = "EditedCorrected",
):
    """
    Read pileup files of reads mapped to a certain transcript from one or more samples into a DataFrame.
    Verify format and change to 0-based coordinates.
    Possibly remove unedited positions, and possibly mask positions in problamatic regions.
    Replace reads' names with a shortened, memory-efficient version, and possibly write the mapping from the original
    reads to their shortened versions to a file.
    Annotate base counts and noise.
    Perform statistical tests for noise, and use these to annotate the basic/naivr editing status of A-to-G positions.
    Then, apply some coverage based filtering, and annotate the final editing status according to the final editing
    scheme.
    Then, also perform statistical tests for editing, and use these to annotate the final editing status.
    Possibly writes the final positions_df to a file, and possibly delete the original pileup files.
    Along the way, perform some sanity checks.

    Args:
        `processes` (int): The number of processes to use in parallel.
        `chroms` (list[str]): The chromosomes of transcripts to process.
        `strands` (list[str]): The strands of transcripts to process.
        `pileup_files_per_chroms` (defaultdict[str, list[Path]]): A dictionary with chroms as keys and lists of pileup files from multiple samples as values.
        `samples_per_chroms` (defaultdict[str, list[str]]): A dictionary with chroms as keys and lists of samples' names (corresponding to `pileup_files_per_chroms`) as values.
        `reads_mapping_files` (list[Path]): Write mapping from the original reads to their shortened versions to these files. There's one file per chrom.
        `positions_files` (list[Path]): Write the positions dataframes, after parsing the pileup and applying statistical tests to the noise and editing levels, to these files.
        `min_percent_of_max_coverage` (float): In each transcript, keep only positions with `coverage` >= `min_percent_of_max_coverage` * `max_coverage`, where `max_coverage` is the maximum coverage in that transcript.
        `min_absolute_coverage` (Union[int, None]): In each transcript, keep only positions with `coverage` >= `min_absolute_coverage`. If None, `min_percent_of_max_coverage` is used instead.
        `snp_noise_level` (float): Treat non-refbase positions with `noise` >= `snp_noise_level` as SNPs, so they won't be considered for noise threshold.
        `top_x_noisy_positions` (int): Use this many positions with highest noise levels to define initial noise threshold.
        `assurance_factor` (float): Multiply the measured noise by this factor to define the actual noise threshold.
        `transcriptome_file` (Path): The transcriptome fasta file.
        `final_editing_scheme` (str): The final editing scheme to use. Can be "BH only" or "BH after noise thresholding".
        `problamatic_regions_file` (Union[Path, str, None], optional): If given, annotate positions as residing in problamatic genomic regions. Defaults to None.
        `known_sites_file` (Union[Path, str, None], optional): If given, annotate positions as known editing sites. Defaults to None.
        `cds_regions_file` (Union[Path, str, None], optional): If given, annotate positions as coding sites. Else, all sites are considered as coding! Defaults to None.
        `out_files_sep` (str, optional): Use this char as a seperator for csv out files. Defaults to "\t".
        `keep_pileup_files` (bool, optional): If True, keep the pileup files from which the positions.csv files are made. Defaults to False.
        `remove_non_refbase_noisy_positions` (bool, optional): If true, remove non refbase positions with too-high noise. Even if False, these positions are not considered for `noise_threshold` determination. Defaults to True.
        `denovo_detection` (bool, optional): Whether to determine both new ("denovo") and known sites as naively edited (not the final decision), or only known ones. Defaults to True.
        `alternative_hypothesis` (str, optional): The alternative hypothesis for the binomial tests, either "larger", "smaller", or "two-sided". Defaults to "greater".
        `binom_noise_pval_col` (str, optional): The name of the column in the returned DataFrame that will contain the p-value of the binomial test for noise. Defaults to "NoiseBinomPVal".
        `bh_noise_pval_col` (str, optional): The name of the column in the returned DataFrame that will contain the corrected p-value for noise after applying the Benjamini-Hochberg procedure. Defaults to "NoiseCorrectedPVal".
        `bh_noisy_col` (str, optional): The name of the column in the returned DataFrame that will contain the information on whether the corrected noise p-value after applying the Benjamini-Hochberg procedure is significant (i.e., "rejected"). Defaults to "NoisyCorrected".
        `binom_editing_pval_col` (str, optional): The name of the column in the returned DataFrame that will contain the p-value of the binomial test for editing. Defaults to "EditingBinomPVal".
        `bh_editing_pval_col` (str, optional): The name of the column in the returned DataFrame that will contain the corrected p-value for editing after applying the Benjamini-Hochberg procedure. Defaults to "EditingCorrectedPVal".
        `bh_editing_col` (str, optional): The name of the column in the returned DataFrame that will contain the information on whether the corrected editing p-value after applying the Benjamini-Hochberg procedure is significant (i.e., "rejected"). Defaults to "EditedCorrected".
    """

    cds_regions_df = pd.read_csv(
        cds_regions_file,
        sep="\t",
        names=["Chrom", "Start", "End", "Name", "Score", "Strand"],
    )

    transcriptome_dict = make_fasta_dict(transcriptome_file)

    # get the positions dfs up until noise correction (not included)
    with Pool(processes=processes) as pool:
        positions_dfs = pool.starmap(
            func=multisample_pileups_to_positions_part_1,
            iterable=[
                (
                    pileup_files_per_chroms[chrom],
                    samples_per_chroms[chrom],
                    strand,
                    problamatic_regions_file,
                    known_sites_file,
                    cds_regions_file,
                    reads_mapping_file,
                    out_files_sep,
                )
                for reads_mapping_file, chrom, strand in zip(
                    reads_mapping_files,
                    chroms,
                    strands,
                )
            ],
        )

    # binom # BH correcttion for noise
    positions_dfs = binom_and_bh_correction_for_noise_all_transcripts(
        positions_dfs,
        processes,
        cds_regions_df,
        transcriptome_dict,
        alternative_hypothesis,
        binom_noise_pval_col,
        bh_noise_pval_col,
        bh_noisy_col,
    )

    # determine "naive" editing status
    with Pool(processes=processes) as pool:
        positions_dfs = pool.starmap(
            func=multisample_pileups_to_positions_part_2,
            iterable=[
                (
                    positions_df,
                    strand,
                    snp_noise_level,
                    top_x_noisy_positions,
                    assurance_factor,
                    remove_non_refbase_noisy_positions,
                    denovo_detection,
                    bh_noisy_col,
                )
                for positions_df, strand in zip(
                    positions_dfs,
                    strands,
                )
            ],
        )

    # binom # BH correcttion for editing
    positions_dfs = binom_and_bh_correction_for_editing_all_transcripts(
        positions_dfs,
        processes,
        cds_regions_df,
        transcriptome_dict,
        alternative_hypothesis,
        binom_editing_pval_col,
        bh_editing_pval_col,
        bh_editing_col,
    )

    # final editing decision, some cleaning, and writing to file
    with Pool(processes=processes) as pool:
        positions_dfs = pool.starmap(
            func=multisample_pileups_to_positions_part_3,
            iterable=[
                (
                    positions_df,
                    pileup_files_per_chroms[chrom],
                    strand,
                    min_percent_of_max_coverage,
                    min_absolute_coverage,
                    final_editing_scheme,
                    positions_file,
                    out_files_sep,
                    keep_pileup_files,
                    bh_editing_col,
                )
                for positions_df, chrom, strand, positions_file in zip(
                    positions_dfs,
                    chroms,
                    strands,
                    positions_files,
                )
            ],
        )


def multisample_pileups_to_positions_part_1(
    pileup_files: list[Path],
    samples: list[str],
    strand: str,
    problamatic_regions_file: Union[Path, str, None] = None,
    known_sites_file: Union[Path, str, None] = None,
    cds_regions_file: Union[Path, str, None] = None,
    reads_mapping_file: Union[Path, str, None] = None,
    out_files_sep: str = "\t",
):
    """
    Read pileup files of reads mapped to a certain transcript from one or more samples into a DataFrame.
    Verify format and change to 0-based coordinates.
    Possibly remove unedited positions, and possibly mask positions in problamatic regions.
    Replace reads' names with a shortened, memory-efficient version, and possibly write the mapping from the original
    reads to their shortened versions to a file.
    Annotate base counts and noise.

    Args:
        `pileup_files` (Path): A list of paths to files creates by the `mpileup` function, each belonging to a different sample.
        `samples`: list[str]: A list of samples' names corresponding to `pileup_files`.
        `strand` (str): The ORF's strand.
        `problamatic_regions_file` (Union[Path, str, None]): If given, annotate positions as residing in problamatic genomic regions.
        `known_sites_file` (Union[Path, str, None]): If given, annotate positions as known editing sites.
        `cds_regions_file` (Union[Path, str, None]): If given, annotate positions as coding sites. Else, *all sites are considered as coding!*.
        `reads_mapping_file` (Union[Path, str, None]): If given, write mapping from the original reads to their shortened versions to this file.
        `out_files_sep` (str): Use this char as a seperator for csv out files. Defaults to tab.


    Returns:
        pd.DataFrame: A df of the parsed pileup files.
    """
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
    # todo uncomment!!!
    # todo uncomment!!!
    # # todo uncomment!!!
    # # todo uncomment!!!
    # # todo uncomment!!!
    # # todo uncomment!!!
    # # todo uncomment!!!
    # # todo uncomment!!!
    # # todo uncomment!!!
    # # todo uncomment!!!
    #
    # replace_reads_names(
    #     positions_df,
    #     "TotalCoverage",
    #     "Reads",
    #     reads_sep=",",
    #     mapping_out_file=reads_mapping_file,
    #     out_files_sep=out_files_sep,
    # )
    # todo uncomment!!!
    # todo uncomment!!!
    # # todo uncomment!!!
    # # todo uncomment!!!
    # # todo uncomment!!!
    # # todo uncomment!!!
    # # todo uncomment!!!
    # # todo uncomment!!!
    # # todo uncomment!!!
    # # todo uncomment!!!

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

    # join multi-sample data per position
    try:
        positions_df = (
            positions_df.groupby(["Chrom", "Position", "RefBase"])
            .agg(
                {
                    "Samples": ",".join,
                    "TotalCoverage": sum,
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

    # # remove positions with insufficient coverage
    # max_coverage = positions_df["TotalCoverage"].max()
    # required_coverage = max_coverage * min_percent_of_max_coverage
    # positions_df = positions_df.loc[positions_df["TotalCoverage"] >= required_coverage]

    annotate_prolamatic_sites(positions_df, strand, problamatic_regions_file)
    annotate_known_sites(positions_df, strand, known_sites_file)
    annotate_coding_sites(positions_df, strand, cds_regions_file)
    annotate_base_counts(positions_df)

    # todo: in future versions, explictly define mismatch per position and, based on that,
    # only then calculate noise
    annotate_noise(positions_df, strand)

    return positions_df


def multisample_pileups_to_positions_part_2(
    positions_df: pd.DataFrame,
    strand: str,
    snp_noise_level: float,
    top_x_noisy_positions: int,
    assurance_factor: float,
    remove_non_refbase_noisy_positions: bool,
    denovo_detection: bool,
    bh_noisy_col: str,
):
    """
    The input positions_df is given after BH correction for noise, and is used to calculate the noise threshold and
    basic editing status according to it.

    Args:
        `positions_df` (pd.DataFrame): The positions df after BH correction for noise.
        `strand` (str): The ORF's strand.
        `snp_noise_level` (float): Treat non-refbase positions with noise >= snp_noise_level as SNPs, so they won't be considered for noise threshold.
        `top_x_noisy_positions` (int): Use this many positions with highest noise levels to define initial noise threshold.
        `assurance_factor` (float): Multiply the measured noise by this factor to define the actual noise threshold.
        `remove_non_refbase_noisy_positions` (bool): If true, remove non refbase positions with too-high noise. Even if this is False, these positions are not considered for noise_threshold determination.
        `denovo_detection` (bool): Whether to determine both new ("denovo") and known sites as naively edited (not the final decision), or only known ones. Defaults to True.
        `bh_noisy_col` (str): The name of the column in the returned DataFrame that will contain the information on whether the corrected noise p-value after applying the Benjamini-Hochberg procedure is significant (i.e., "rejected").

    Returns:
        pd.DataFrame: A df after basic/naive editing status annotation, and before BH correction for editing.
    """

    ref_base = "A" if strand == "+" else "T"

    #  possibly remove non refbase positions with too-high noise & define noise threshold
    if remove_non_refbase_noisy_positions:
        # remove non refbase positions with too-high noise or whose corrected noise p-value is not significant
        positions_before = len(positions_df)
        positions_df = positions_df.loc[
            (positions_df["RefBase"] == ref_base)
            | ((positions_df["Noise"] < snp_noise_level) & (positions_df[bh_noisy_col]))
        ]
        positions_after = len(positions_df)
        removed_positions = positions_before - positions_after
        print(
            f"{removed_positions} extremely- or insignificantly-noisy positions removed"
        )
        # define noise threshold
        noise_threshold = (
            positions_df["Noise"]
            .sort_values(ascending=False)[:top_x_noisy_positions]
            .mean()
        )
    else:
        # define noise threshold
        noise_threshold = (
            positions_df.loc[
                (positions_df["Noise"] < snp_noise_level)
                & (positions_df[bh_noisy_col]),
                "Noise",
            ]
            .sort_values(ascending=False)[:top_x_noisy_positions]
            .mean()
        )
    if pd.isna(noise_threshold):
        noise_threshold = 0
    noise_threshold *= assurance_factor

    # annotate editing frequency and "naive" editing status according to noise threshold
    annotate_editing_frequency_per_position(positions_df, strand)
    annotate_edited_sites(positions_df, strand, noise_threshold, denovo_detection)

    return positions_df


def multisample_pileups_to_positions_part_3(
    positions_df: pd.DataFrame,
    pileup_files: list[Path],
    strand: str,
    min_percent_of_max_coverage: float,
    min_absolute_coverage: Union[int, None],
    final_editing_scheme: str,
    positions_out_file: Union[Path, str, None] = None,
    out_files_sep: str = "\t",
    keep_pileup_files: bool = False,
    bh_editing_col: str = "EditedCorrected",
):
    """
    The input positions_df is given after BH correction for noise & editing.
    Then, apply some coverage based filtering, and annotate the final editing status according to the final editing scheme.
    In addition, some verifications are made.
    Possibly writes the final positions_df to a file, and possibly delete the original pileup files.

    Args:
        `positions_df` (pd.DataFrame): The positions df after BH correction for noise.
        `pileup_files` (Path): A list of paths to files creates by the `mpileup` function from reads mapped to the transcript, each belonging to a different sample.
        `strand` (str): The transcript's strand.
        `min_percent_of_max_coverage` (float): In each transcript, keep only positions with `coverage` >= `min_percent_of_max_coverage` * `max_coverage`, where `max_coverage` is the maximum coverage in that transcript.
        `min_absolute_coverage` (Union[int, None]): In each transcript, keep only positions with `coverage` >= `min_absolute_coverage`. If None, `min_percent_of_max_coverage` is used instead.
        `final_editing_scheme` (str): The final editing scheme to use. Can be "BH only" or "BH after noise thresholding".
        `positions_out_file` (Union[Path, str, None]): If given, write positions_df to this file.
        `out_files_sep` (str): Use this char as a seperator for csv out files. Defaults to tab.
        `bh_editing_col` (str): The name of the column in the returned DataFrame that will contain the information on whether the corrected editing p-value after applying the Benjamini-Hochberg procedure is significant (i.e., "rejected"). Defaults to "EditedCorrected".

    Returns:
        pd.DataFrame: The final positions df after BH correction for noise & editing of that transcript.
    """

    ref_base = "A" if strand == "+" else "T"

    # remove positions with insufficient coverage defined by absolute or relative coverage criteria
    # (we can finally do this after BH correction for noise & editing)
    if min_absolute_coverage is not None:
        required_coverage = min_absolute_coverage
    else:
        max_coverage = positions_df["TotalCoverage"].max()
        required_coverage = max_coverage * min_percent_of_max_coverage
    positions_df = positions_df.loc[positions_df["TotalCoverage"] >= required_coverage]

    # verify that noise is only applied to non ref_base positions, and vice versa for editing frequency
    non_ref_base_edit_freq = positions_df.loc[
        positions_df["RefBase"] != ref_base, "EditingFrequency"
    ].unique()
    assert len(non_ref_base_edit_freq) == 1 and np.isnan(non_ref_base_edit_freq[0])
    ref_base_noise = positions_df.loc[
        positions_df["RefBase"] == ref_base, "Noise"
    ].unique()
    assert len(ref_base_noise) == 1 and np.isnan(ref_base_noise[0])

    # verify that the number of mapped bases ==
    # number of reads' names ==
    # number of samples
    for row_num, row in enumerate(positions_df.itertuples()):
        mapped_bases = row.MappedBases
        mapped_reads = row.Reads.split(",")
        mapped_samples = row.Samples.split(",")
        if not (len(mapped_bases) == len(mapped_reads) == len(mapped_samples)):
            raise ValueError

    # finally, annotate editing based on one of two possible schemes:
    # 1. based on the corrected editing p-values alone
    # 2. based on the corrected editing p-values, and also on the noise threshold (which is also corrected), as defined
    # by the "naive" editing status
    if final_editing_scheme == "BH only":
        editing_cols = [bh_editing_col]
    elif final_editing_scheme == "BH after noise thresholding":
        editing_cols = [bh_editing_col, "Edited"]
    else:
        raise ValueError(
            f"{final_editing_scheme = } is not a valid final_editing_scheme, must be either 'BH only' or 'BH after noise thresholding'"
        )
    annotate_edited_sites_final_decision(positions_df, bh_editing_col, editing_cols)

    positions_df = positions_df.loc[positions_df["TotalCoverage"] > 0]
    if len(positions_df) == 0:
        return

    if positions_out_file:
        positions_df.to_csv(
            positions_out_file, sep=out_files_sep, index=False, na_rep=np.nan
        )

    if not keep_pileup_files:
        for pileup_file in pileup_files:
            subprocess.run(f"rm {pileup_file}", shell=True)


def binom_and_bh_correction_for_noise_all_transcripts(
    positions_dfs: list[pd.DataFrame],
    processes: int,
    cds_regions_df: pd.DataFrame,
    transcriptome_dict: dict[str, Seq.Seq],
    alternative_hypothesis: str,
    binom_pval_col: str,
    bh_pval_col: str,
    bh_rejection_col: str,
    out_file: Path | None = None,
) -> list[pd.DataFrame]:
    with Pool(processes=processes) as pool:
        per_transcript_coding_non_adenosines_dfs = pool.starmap(
            func=get_covered_and_uncovered_coding_non_adenosines_in_transcript,
            iterable=[
                (
                    positions_df,
                    chrom,
                    start,
                    end,
                    strand,
                    transcriptome_dict,
                    alternative_hypothesis,
                    binom_pval_col,
                )
                for positions_df, (chrom, start, end, strand) in zip(
                    positions_dfs,
                    cds_regions_df[["Chrom", "Start", "End", "Strand"]].itertuples(
                        index=False
                    ),
                )
            ],
        )

    all_transcripts_coding_non_adenosines_df = pd.concat(
        per_transcript_coding_non_adenosines_dfs
    ).reset_index(drop=True)

    bh_rejections, bh_corrected_pvals = fdrcorrection(
        all_transcripts_coding_non_adenosines_df[binom_pval_col]
    )

    all_transcripts_coding_non_adenosines_df[bh_pval_col] = bh_corrected_pvals
    all_transcripts_coding_non_adenosines_df[bh_rejection_col] = bh_rejections

    # merge back the corrected p-values and rejections to the original positions_dfs
    corrected_positions_dfs = []
    for positions_df in positions_dfs:
        # if len(positions_df) == 0:
        #     continue
        corrected_positions_df = positions_df.merge(
            all_transcripts_coding_non_adenosines_df.loc[
                :,
                [
                    "Chrom",
                    "Position",
                    binom_pval_col,
                    bh_pval_col,
                    bh_rejection_col,
                ],
            ],
            on=["Chrom", "Position"],
            how="left",
        )
        corrected_positions_dfs.append(corrected_positions_df)

    if out_file is not None:
        all_transcripts_coding_non_adenosines_df.to_csv(out_file, sep="\t", index=False)

    return corrected_positions_dfs


def binom_and_bh_correction_for_editing_all_transcripts(
    positions_dfs: list[pd.DataFrame],
    processes: int,
    cds_regions_df: pd.DataFrame,
    transcriptome_dict: dict[str, Seq.Seq],
    alternative_hypothesis: str,
    binom_pval_col: str,
    bh_pval_col: str,
    bh_rejection_col: str,
    out_file: Path | None = None,
) -> list[pd.DataFrame]:
    """
    Perform binomial test for editing followed by Benjamini-Hochberg correction for all adenosines in the
    transcriptome.

    Args:
        `positions_dfs` (list[pd.DataFrame]): A list of DataFrames, where each DataFrame contains the positions of adenosines in a transcript.
        `processes` (int): The number of processes to use in parallel.
        `cds_regions_df` (Path): A DataFrame represent a BED-6 coordinates of the transcriptome's CDS.
        `transcriptome_dict` (dict[str, Seq.Seq]): A dictionary of the transcriptome, where the keys are the chromosomes and the values are the sequences.
        `alternative_hypothesis` (str, optional): The alternative hypothesis for the binomial test, either "larger", "smaller", or "two-sided".
        `binom_pval_col` (str, optional): The name of the column in the returned DataFrame that will contain the p-value of the binomial test.
        `bh_pval_col` (str, optional): The name of the column in the returned DataFrame that will contain the corrected p-value after applying the Benjamini-Hochberg procedure.
        `bh_rejection_col` (str, optional): The name of the column in the returned DataFrame that will contain the information on whether the corrected p-value after applying the Benjamini-Hochberg procedure is significant (i.e., "rejected").
        `out_file` (Path, optional): If given, write all corrected positions from all transcripts to this file.

    Returns:
        list[pd.DataFrame]: The same list of DataFrames as `positions_dfs`, but with the binomial p-values, corrected
        p-values and corresponding rejections.
    """

    with Pool(processes=processes) as pool:
        per_transcript_coding_adenosines_dfs = pool.starmap(
            func=get_covered_and_uncovered_coding_adenosines_in_transcript,
            iterable=[
                (
                    positions_df,
                    chrom,
                    start,
                    end,
                    strand,
                    transcriptome_dict,
                    alternative_hypothesis,
                    binom_pval_col,
                )
                for positions_df, (chrom, start, end, strand) in zip(
                    positions_dfs,
                    cds_regions_df[["Chrom", "Start", "End", "Strand"]].itertuples(
                        index=False
                    ),
                )
            ],
        )

    all_transcripts_coding_adenosines_df = pd.concat(
        per_transcript_coding_adenosines_dfs
    ).reset_index(drop=True)

    bh_rejections, bh_corrected_pvals = fdrcorrection(
        all_transcripts_coding_adenosines_df[binom_pval_col]
    )

    all_transcripts_coding_adenosines_df[bh_pval_col] = bh_corrected_pvals
    all_transcripts_coding_adenosines_df[bh_rejection_col] = bh_rejections

    # merge back the corrected p-values and rejections to the original positions_dfs
    corrected_positions_dfs = []
    for positions_df in positions_dfs:
        corrected_positions_df = positions_df.merge(
            all_transcripts_coding_adenosines_df.loc[
                :,
                [
                    "Chrom",
                    "Position",
                    binom_pval_col,
                    bh_pval_col,
                    bh_rejection_col,
                ],
            ],
            on=["Chrom", "Position"],
            how="left",
        )
        corrected_positions_dfs.append(corrected_positions_df)

    if out_file is not None:
        all_transcripts_coding_adenosines_df.to_csv(out_file, sep="\t", index=False)

    return corrected_positions_dfs


def get_covered_and_uncovered_coding_non_adenosines_in_transcript(
    positions_df: pd.DataFrame,
    chrom: str,
    start: int,
    end: int,
    strand: str,
    transcriptome_dict: dict[str, Seq.Seq],
    alternative_hypothesis: str,
    binom_pval_col: str,
) -> pd.DataFrame:
    """
    Get a DataFrame of coding non-adenosine bases positions in a transcript with the p-value of binomial test for the
    noise in each position.

    The positions are taken from a two sources: the `positions_df` which contains covered non-adenosine bases (i.e.,
    positions with coverage > 0), and all the other coding non-adenosine bases in the transcript (their binomial P-value
    is 1.0 by definition).

    Note! This function assumes that the coverage of bases not in `positions_df` is exactly 0
    (however, it's possible that the entire coverage in a position is due to unkown bases (N)).
    Any other coverage-based filter should be applied only after this function is called.

    Args:
        `positions_df` (pd.DataFrame): A DataFrame of mapped positions in the transcript, where each row is a position.
        `chrom` (str): The chromosome of the transcript.
        `start` (int): The start position of the transcript, 0-based.
        `end` (int): The end position of the transcript, exclusive.
        `strand` (str): The strand of the transcript.
        `transciptome_dict` (dict[str, Seq.Seq]): A dictionary of the transcriptome, where the keys are the chromosomes and the values are the sequences.
        `alternative_hypothesis` (str): The alternative hypothesis for the binomial test, either "larger", "smaller", or "two-sided".
        `binom_pval_col` (str): The name of the column in the returned DataFrame that will contain the p-value of the binomial test.

    Returns:
        pd.DataFrame: _description_
    """

    # retain only coding non-adenosine positions with total coverage > 0
    negative_ref_base = "A" if strand == "+" else "T"

    covered_coding_non_adenosines_positions_df = positions_df.loc[
        (positions_df["CDS"]) & (positions_df["RefBase"] != negative_ref_base),
        [
            "Chrom",
            "Position",
            "RefBase",
            "A",
            "T",
            "C",
            "G",
            "TotalCoverage",
            "Noise",
        ],
    ]

    # get all coding non-adenosines in the transcript
    # (including those with total coverage > 0, which are already in `positions_df`)
    # this df's columns are (1) Chrom, (2) Position and (3) RefBase
    all_coding_non_adenosines_df = get_positions_of_desired_bases_in_transcript(
        chrom, start, end, strand, transcriptome_dict, {"T", "C", "G"}
    )

    # add some columns to all_coding_non_adenosines_df before returning it
    # (we don't have any covered positions with desired bases in this transcript,
    # so we can't calculate the binomial p-values for them)
    if len(covered_coding_non_adenosines_positions_df) == 0:
        all_coding_non_adenosines_df[binom_pval_col] = 1.0
        all_coding_non_adenosines_df["Noise"] = 0.0
        for col in [
            "TotalCoverage",
            "A",
            "T",
            "C",
            "G",
            "RefBaseCount",
            "AltBaseCount",
        ]:
            all_coding_non_adenosines_df[col] = 0
        return all_coding_non_adenosines_df

    # perform binomial test on each covered position
    ref_and_alt_base_counts_df = covered_coding_non_adenosines_positions_df.apply(
        lambda x: ref_and_alt_base_count_in_noise_position(
            x["RefBase"],
            strand,
            x["A"],
            x["T"],
            x["C"],
            x["G"],
        ),
        axis=1,
        result_type="expand",
    ).rename(columns={0: "RefBaseCount", 1: "AltBaseCount"})
    covered_coding_non_adenosines_positions_df[["RefBaseCount", "AltBaseCount"]] = (
        ref_and_alt_base_counts_df
    )
    covered_coding_non_adenosines_positions_df[binom_pval_col] = (
        covered_coding_non_adenosines_positions_df.apply(
            lambda x: binom_test(
                x["AltBaseCount"],
                x["RefBaseCount"] + x["AltBaseCount"],
                0.001,
                alternative=alternative_hypothesis,
            ),
            axis=1,
        )
    )

    # merge the two DataFrames, such that the positions in `positions_df` will have the P-value from the binomial test,
    # as well as original editing status (if available), the total coverage, and the number of Gs
    all_coding_non_adenosines_df = all_coding_non_adenosines_df.merge(
        covered_coding_non_adenosines_positions_df,
        on=["Chrom", "Position", "RefBase"],
        how="left",
    ).sort_values(by=["Chrom", "Position"], ignore_index=True)

    # fill in the missing values for the uncovered positions
    all_coding_non_adenosines_df[binom_pval_col] = all_coding_non_adenosines_df[
        binom_pval_col
    ].fillna(1.0)
    for col in [
        "TotalCoverage",
        "A",
        "T",
        "C",
        "G",
        "Noise",
        "RefBaseCount",
        "AltBaseCount",
    ]:
        all_coding_non_adenosines_df[col] = all_coding_non_adenosines_df[col].fillna(0)

    # return covered_coding_non_adenosines_positions_df
    return all_coding_non_adenosines_df


def get_covered_and_uncovered_coding_adenosines_in_transcript(
    positions_df: pd.DataFrame,
    chrom: str,
    start: int,
    end: int,
    strand: str,
    transcriptome_dict: dict[str, Seq.Seq],
    alternative_hypothesis: str,
    binom_editing_pval_col: str,
) -> pd.DataFrame:
    """
    Get a DataFrame of coding adenosines positions in a transcript with the p-value of binomial test for the editing in
    each position.

    The positions are taken from a two sources: the `positions_df` which contains covered adenosiens (i.e.,
    positions with coverage > 0), and all the other coding adenosiens in the transcript (their binomial P-value is 1.0
    by definition).

    Note! This function assumes that the coverage of adenosines not in `positions_df` is exactly 0
    (however, it's possible that the entire coverage in a position is due to unkown bases (N)).
    Any other coverage-based filter should be applied only after this function is called.

    Args:
        positions_df (pd.DataFrame): A DataFrame of mapped positions in the transcript, where each row is a position.
        chrom (str): The chromosome of the transcript.
        start (int): The start position of the transcript, 0-based.
        end (int): The end position of the transcript, exclusive.
        strand (str): The strand of the transcript.
        transciptome_dict (dict[str, Seq.Seq]): A dictionary of the transcriptome, where the keys are the chromosomes
        and the values are the sequences.
        alternative (str): The alternative hypothesis for the binomial test, either "larger", "smaller", or "two-sided".
        binom_editing_pval_col (str): The name of the column in the returned DataFrame that will contain the p-value of
        the binomial test.

    Returns:
        pd.DataFrame: _description_
    """
    # retain only coding adenosines positions with total coverage > 0
    # (however, it's possible that the entire coverage in a position is due to unkown bases (N))
    ref_base = "A" if strand == "+" else "T"
    covered_coding_adenosines_positions_df = positions_df.loc[
        (positions_df["CDS"]) & (positions_df["RefBase"] == ref_base),
        [
            "Chrom",
            "Position",
            "RefBase",
            "A",
            "T",
            "C",
            "G",
            "G",
            "TotalCoverage",
            "Edited",
        ],
    ]

    # perform binomial test on each covered position
    ref_and_alt_base_counts_df = covered_coding_adenosines_positions_df.apply(
        lambda x: ref_and_alt_base_count_in_editing_position(
            strand,
            x["A"],
            x["T"],
            x["C"],
            x["G"],
        ),
        axis=1,
        result_type="expand",
    ).rename(columns={0: "RefBaseCount", 1: "AltBaseCount"})
    covered_coding_adenosines_positions_df[["RefBaseCount", "AltBaseCount"]] = (
        ref_and_alt_base_counts_df
    )
    covered_coding_adenosines_positions_df[binom_editing_pval_col] = (
        covered_coding_adenosines_positions_df.apply(
            lambda x: binom_test(
                x["AltBaseCount"],
                x["RefBaseCount"] + x["AltBaseCount"],
                0.001,
                alternative=alternative_hypothesis,
            ),
            axis=1,
        )
    )

    # get all coding adenosines in the transcript (including the coverage > 0, which are already in `positions_df`)
    all_coding_adenosines_df = get_positions_of_desired_bases_in_transcript(
        chrom, start, end, strand, transcriptome_dict, {"A"}
    )

    # merge the two DataFrames, such that the positions in `positions_df` will have the P-value from the binomial test,
    # as well as original editing status (if available), the total coverage, and the number of Gs
    all_coding_adenosines_df = all_coding_adenosines_df.merge(
        covered_coding_adenosines_positions_df, on=["Chrom", "Position"], how="left"
    ).sort_values(by=["Chrom", "Position"], ignore_index=True)

    # fill in the missing values for the uncovered positions
    all_coding_adenosines_df["Edited"] = all_coding_adenosines_df["Edited"].fillna(
        False
    )
    all_coding_adenosines_df[binom_editing_pval_col] = all_coding_adenosines_df[
        binom_editing_pval_col
    ].fillna(1.0)
    for col in ["TotalCoverage", "A", "T", "C", "G", "RefBaseCount", "AltBaseCount"]:
        all_coding_adenosines_df[col] = all_coding_adenosines_df[col].fillna(0)

    # return covered_coding_adenosines_positions_df
    return all_coding_adenosines_df


def get_positions_of_desired_bases_in_transcript(
    chrom: str,
    start: int,
    end: int,
    strand: str,
    transcriptome_dict: dict[str, Seq.Seq],
    desired_bases: set[str],
) -> pd.DataFrame:
    """
    Get all positions in a transcript whose reference base is in `desired_bases`.

    Args:
        `chrom` (str): The chromosome of the transcript.
        `start` (int): The first position in the transcript to be considered, 0-based.
        `end` (int): The last position in the transcript to be considered, exclusive.
        `strand` (str): The strand of the transcript.
        `transciptome_dict` (dict[str, Seq.Seq]): A dictionary of the transcriptome, where the keys are the chromosomes and the values are the sequences.
        `desired_bases` (set[str]): A set of bases to be considered. They are given w.r.t the biological data, so if the strand is "-", the desired bases are the complements of the given ones.

    Returns:
        pd.DataFrame: A dataframe, where each row contains the chromosome, the (start) position, and the reference of a
        desired base in the transcript.
    """
    if strand == "-":
        rev_comp_dict = {"A": "T", "T": "A", "C": "G", "G": "C"}
        desired_bases = {rev_comp_dict[base] for base in desired_bases}

    positions_of_desired_bases = []

    seq = transcriptome_dict[chrom]

    for x in range(start, end):
        base = seq[x : x + 1]
        # if strand == "-":
        #     base = base.reverse_complement()
        if base in desired_bases:
            positions_of_desired_bases.append((chrom, x, base))

    positions_of_desired_bases_df = pd.DataFrame(
        positions_of_desired_bases, columns=["Chrom", "Position", "RefBase"]
    )

    return positions_of_desired_bases_df


def ref_and_alt_base_count_in_editing_position(
    strand: str, a_count: int, t_count: int, c_count: int, g_count: int
) -> tuple[int, int]:
    if strand == "+":
        ref_base_count, alt_base_count = a_count, g_count
    else:
        ref_base_count, alt_base_count = t_count, c_count
    return ref_base_count, alt_base_count


def ref_and_alt_base_count_in_noise_position(
    ref_base: str, strand: str, a_count: int, t_count: int, c_count: int, g_count: int
) -> tuple[int, int]:
    """
    Return the number of reads supporting `ref_base` and `alt_base`, where `alt_base` is the most abundant base other
    than `ref_base` in the position.

    Note: This function assumes that the reference base is not an adenosine, and that the total coverage > 0.
    However, it's possible that the total coverage in the position is made entirely of Ns, and thus,
    `sum([a_count, t_count, c_count, g_count]) == 0`.


    Args:
        `ref_base` (str): The reference base of the position.
        `strand` (str): The transcript's strand.
        `a_count` (int): Number of As mapped to position.
        `t_count` (int): Number of Ts mapped to position.
        `c_count` (int): Number of Cs mapped to position.
        `g_count` (int): Number of Gs mapped to position.

    Returns:
        tuple[int, int]: The numbers of reads supporting `ref_base` and `alt_base`, respectively.
    """
    base_counts = [a_count, t_count, c_count, g_count]

    # if (
    #     (strand == "+" and ref_base == "A") or (strand == "-" and ref_base == "T")
    # ) and sum(base_counts) == 0:
    #     raise ValueError(
    #         "This function is only for non-adenosine positions with coverage > 0."
    #     )

    # elif (strand == "+" and ref_base == "A") or (strand == "-" and ref_base == "T"):
    #     raise ValueError("This function is only for non-adenosine positions.")
    # elif sum(base_counts) == 0:
    #     raise ValueError("This function is only for positions with coverage > 0.")

    if (strand == "+" and ref_base == "A") or (strand == "-" and ref_base == "T"):
        raise ValueError("This function is only for non-adenosine positions.")

    # we measure noise for T positions only on the positive strand
    if strand == "+" and ref_base == "T":
        ref_base_count = t_count
    # we measure noise for A positions only on the negative strand
    elif strand == "-" and ref_base == "A":
        ref_base_count = a_count
    # we measure noise for A & G positions on both strand
    elif ref_base == "C":
        ref_base_count = c_count
    elif ref_base == "G":
        ref_base_count = g_count
    else:
        raise ValueError(
            "Something weird with alt_base_count w.r.t. ref_base and strand - check it out!"
        )

    bases = ["A", "T", "C", "G"]
    alt_bases = set(bases) - {ref_base}
    alt_base_counts = [
        base_count for base_count, base in zip(base_counts, bases) if base in alt_bases
    ]
    max_alt_base_count = max(alt_base_counts)

    return ref_base_count, max_alt_base_count


def annotate_edited_sites_final_decision(
    positions_df: pd.DataFrame,
    bh_editing_col: str,
    editing_cols: list[str],
) -> None:
    """
    Determine the final editing status of each position, based on the agreement between the editing columns
    `editing_cols`.

    One of the editing cols is the BH-corrected editing column `bh_editing_col`. The other (if any) is the "Edited" col,
    which is the result of the noise thresholding.

    Args:
        positions_df (pd.DataFrame): The positions dataframe.
        bh_editing_col (str): The BH-corrected editing column.
        editing_cols (list[str]): The editing columns to consider.
    """
    final_editing_status = positions_df.apply(
        lambda x: all([x[col] for col in editing_cols]), axis=1
    )
    positions_df.insert(
        positions_df.columns.get_loc(bh_editing_col) + 1,
        "EditedFinal",
        final_editing_status,
    )


# positions_files_sep = "\t"
# processes = 20
# positions_files_dir = Path(
#     "/private7/projects/Combinatorics/O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50/PositionsFiles"
# )
# positions_files = list(positions_files_dir.glob("*.positions.csv.gz"))

# covered_coding_adenosines_df = get_covered_coding_adenosines_positions_df(
#     positions_files, positions_files_sep, processes
# )


# cds_regions_file = Path("O.vulgaris/Annotations/orfs_oct.bed").absolute()
# transcriptome_file = Path("O.vulgaris/Annotations/orfs_oct.fa").absolute()

# cds_regions_df = pd.read_csv(
#     cds_regions_file,
#     sep="\t",
#     names=["Chrom", "Start", "End", "Name", "Score", "Strand"],
# )

# transciptome_dict = make_fasta_dict(transcriptome_file)

# with Pool(processes=processes) as pool:
#     per_transcript_coding_adenosines = pool.starmap(
#         func=get_all_coding_adenosines_in_transcript,
#         iterable=[
#             (chrom, start, end, strand, transciptome_dict)
#             for chrom, start, end, strand in cds_regions_df[
#                 ["Chrom", "Start", "End", "Strand"]
#             ].itertuples(index=False)
#         ],
#     )

# all_coding_adenosines = [
#     coding_adenosine
#     for coding_adenosines in per_transcript_coding_adenosines
#     for coding_adenosine in coding_adenosines
# ]

# all_coding_adenosines_df = pd.DataFrame(
#     all_coding_adenosines, columns=["Chrom", "Position"]
# )


# merged_adenosines_df = all_coding_adenosines_df.merge(
#     covered_coding_adenosines_df, on=["Chrom", "Position"], how="left"
# ).sort_values(by=["Chrom", "Position"], ignore_index=True)

# merged_adenosines_df["Edited"] = merged_adenosines_df["Edited"].fillna(False)
# merged_adenosines_df["PVal"] = merged_adenosines_df["PVal"].fillna(1.0)

# bh_rejections, bh_corrected_pvals = fdrcorrection(merged_adenosines_df["PVal"])

# merged_adenosines_df["BHRejection"] = bh_rejections
# merged_adenosines_df["BHCorrectedPVal"] = bh_corrected_pvals


# tmr50_complete_data_chroms_df = pd.read_table(
#     "/private7/projects/Combinatorics/Code/Notebooks/TMR50.CompleteData.Chroms.tsv"
# )
# tmr50_complete_data_chroms_df

# tmr50_complete_data_merged_adenosines_df = merged_adenosines_df.loc[
#     merged_adenosines_df["Chrom"].isin(tmr50_complete_data_chroms_df["Chrom"])
# ]
# tmr50_complete_data_merged_adenosines_df


# tmr50_complete_data_merged_adenosines_df[
#     "Edited"
# ].sum()  # 56626 - exactly as in Fig. S2a
# tmr50_complete_data_merged_adenosines_df.loc[
#     tmr50_complete_data_merged_adenosines_df["BHRejection"]
#     & tmr50_complete_data_merged_adenosines_df["Edited"]
# ]  # 47860 sites: 8766 less than those we report in Fig. S2a


# tmr50_complete_data_merged_adenosines_df.loc[
#     tmr50_complete_data_merged_adenosines_df["BHRejection"]
# ]


# newly_edited_df = tmr50_complete_data_merged_adenosines_df.loc[
#     tmr50_complete_data_merged_adenosines_df["BHRejection"]
#     & ~tmr50_complete_data_merged_adenosines_df["Edited"]
# ]
# newly_edited_df["%EditingFrequency"] = (
#     100 * newly_edited_df["G"] / newly_edited_df["TotalCoverage"]
# )

# percentiles = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]

# newly_edited_df["TotalCoverage"].describe(percentiles=percentiles)
# newly_edited_df["G"].describe(percentiles=percentiles)

# newly_edited_df.loc[newly_edited_df["TotalCoverage"] <= 123, "G"].describe(
#     percentiles=percentiles
# )
# newly_edited_df.loc[
#     newly_edited_df["TotalCoverage"] <= 123, "%EditingFrequency"
# ].describe(percentiles=percentiles)


# newly_edited_df["%EditingFrequency"].describe()
# newly_edited_df["%EditingFrequency"].describe(percentiles=percentiles)
