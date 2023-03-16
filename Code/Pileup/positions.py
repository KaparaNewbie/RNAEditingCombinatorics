import subprocess
from pathlib import Path
from typing import Union
from itertools import chain

from pybedtools import BedTool
import numpy as np
import pandas as pd
from icecream import ic

from Pileup.ids import compress_ids


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

    # replace
    old_unique_reads = set(
        chain.from_iterable(positions_df[reads_col].str.split(reads_sep))
    )
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
    if problamatic_regions_file:
        problamatic_regions_bed = BedTool(problamatic_regions_file)
        positions_bed = positions_df_to_bed(positions_df, strand)
        prob_positions_bed = positions_bed.intersect(problamatic_regions_bed)
        prob_positions = {interval.start for interval in prob_positions_bed}

        def is_prob(position: int) -> bool:
            return position in prob_positions

        in_prob_region = [is_prob(position) for position in positions_df["Position"]]
    else:
        in_prob_region = [False for _ in range(len(positions_df))]
    positions_df.insert(
        positions_df.columns.get_loc("Position") + 1, "InProbRegion", in_prob_region
    )


def annotate_known_sites(positions_df, strand, known_sites_file=None):
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
        positions_df.columns.get_loc("Position") + 1, "KnownEditing", is_position_known
    )


def annotate_coding_sites(positions_df, strand, cds_regions_file=None):
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


def count_bases(ref_base, mapped_bases):
    mapped_bases = mapped_bases.replace(".", ref_base)
    return (
        mapped_bases.count("A"),
        mapped_bases.count("T"),
        mapped_bases.count("C"),
        mapped_bases.count("G"),
    )


def annotate_base_counts(positions_df):
    atcgs = positions_df.apply(
        lambda x: count_bases(x["RefBase"], x["MappedBases"]),
        axis=1,
        result_type="expand",
    ).rename(columns={0: "A", 1: "T", 2: "C", 3: "G"})
    positions_df[["A", "T", "C", "G"]] = atcgs


def position_noise_level(
    ref_base: str, strand: str, a_count: int, t_count: int, c_count: int, g_count: int
) -> float:
    """
    Calculate noise level in a position.

    The noise is defined as number of bases of the must abundant alt_base, divided by the same number +
    the number of ref_base bases.

    Args:
        ref_base (str): The reference base of the position.
        strand (str): The ORF's strand.
        a_count (int): Number of A's mapped to position.
        t_count (int): Number of T's mapped to position.
        c_count (int): Number of C's mapped to position.
        g_count (int): Number of G's mapped to position.

    Returns:
        float: The noise level.
    """
    base_counts = [a_count, t_count, c_count, g_count]

    if sum(base_counts) == 0:
        noise = np.NaN

    # we only measure noise in positions that don't undergo RNA editing by ADAR
    elif (strand == "+" and ref_base == "A") or (strand == "-" and ref_base == "T"):
        noise = np.NaN

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
    positions_df["Noise"] = positions_df.apply(
        lambda x: position_noise_level(
            x["RefBase"], strand, x["A"], x["T"], x["C"], x["G"]
        ),
        axis=1,
    )


def editing_frequency_per_position(
    ref_base: str, base: str, ref_base_count: int, alt_base_count: int
) -> float:
    if base == ref_base:
        try:
            freq = alt_base_count / (ref_base_count + alt_base_count)
        except ZeroDivisionError:
            freq = 0
    else:
        freq = np.NaN
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
    replace_reads_names(
        positions_df,
        "TotalCoverage",
        "Reads",
        reads_sep=",",
        mapping_out_file=reads_mapping_file,
        out_files_sep=out_files_sep,
    )

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
            positions_out_file, sep=out_files_sep, index=False, na_rep=np.NaN
        )

    if not keep_pileup_file:
        subprocess.run(f"rm {pileup_file}", shell=True)

    # return positions_df


def multisample_pileups_to_positions(
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
        if not(len(mapped_bases) == len(mapped_reads) == len(mapped_samples)):
            raise ValueError

    if positions_out_file:
        positions_df.to_csv(
            positions_out_file, sep=out_files_sep, index=False, na_rep=np.NaN
        )

    if not keep_pileup_file:
        subprocess.run(f"rm {pileup_file}", shell=True)

    # return positions_df
