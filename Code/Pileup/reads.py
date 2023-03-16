from itertools import chain
from typing import Union
from pathlib import Path

import pandas as pd
import numpy as np
from icecream import ic

from Pileup.ids import serialize_compressed_ids


def count_values(row, first_pos_loc, value):
    return len(row.iloc[first_pos_loc:].loc[row.iloc[first_pos_loc:] == value])


def annotate_edited_positions(df, first_pos_loc):
    df.insert(
        first_pos_loc,
        "EditedPositions",
        df.apply(lambda row: count_values(row, first_pos_loc, 1), axis="columns"),
    )


def annotate_unedited_positions(df, first_pos_loc):
    df.insert(
        first_pos_loc,
        "UneditedPositions",
        df.apply(lambda row: count_values(row, first_pos_loc, 0), axis="columns"),
    )


def annotate_ambigous_positions(df, first_pos_loc):
    df.insert(
        first_pos_loc,
        "AmbigousPositions",
        df.apply(lambda row: count_values(row, first_pos_loc, -1), axis="columns"),
    )


def annotate_editing_frequency_per_read(reads_df, new_col_pos, positions_df, strand):
    ref_base = "A" if strand == "+" else "T"
    all_refbase_positions_df = positions_df.loc[
        (positions_df["RefBase"] == ref_base) & (~positions_df["InProbRegion"])
    ]
    all_refbase_positions = len(all_refbase_positions_df)
    reads_df.insert(
        new_col_pos,
        "EditingFrequency",
        reads_df.apply(
            lambda row: row["EditedPositions"] / all_refbase_positions, axis="columns"
        ),
    )


def add_pacbio_se_reads_in_positions(edited_positions_df, strand, group):
    unique_reads = set(chain.from_iterable(edited_positions_df["Reads"].str.split(",")))
    positions_per_read = {read: [] for read in unique_reads}

    alt_base = "G" if strand == "+" else "C"
    for x, position_row in enumerate(edited_positions_df.itertuples()):
        mapped_bases = position_row.MappedBases
        mapped_reads = position_row.Reads.split(",")
        # deal with reads that were mapped to this pos
        for mapped_base, mapped_read in zip(mapped_bases, mapped_reads):
            # 1 = A2G/T2C editing, 0 = match
            if mapped_base == alt_base:
                # pos_edited_in_read = "1" # A2G/T2C editing
                pos_edited_in_read = 1  # A2G/T2C editing
            elif mapped_base == ".":
                # pos_edited_in_read = "0" # match
                pos_edited_in_read = 0  # match
            else:
                # pos_edited_in_read = "-1" # we ignore mismatches that aren't RNA editing by marking them as NaN
                pos_edited_in_read = (
                    -1
                )  # we ignore mismatches that aren't RNA editing by marking them as NaN
            positions_per_read[mapped_read].append(pos_edited_in_read)
        # deal with reads that weren't mapped to this pos
        unmapped_reads_in_pos = unique_reads - set(mapped_reads)
        for mapped_read in unmapped_reads_in_pos:
            # pos_edited_in_read = "-1"
            pos_edited_in_read = -1
            positions_per_read[mapped_read].append(pos_edited_in_read)
        # check all reads got values for pos
        mapped_pos_dist = {len(bases) for bases in positions_per_read.values()}
        if len(mapped_pos_dist) != 1:
            raise Exception(
                f"Problem at line {x} in {group}'s positions_df: not all reads are mapped."
            )

    return positions_per_read


def add_illumina_pe_reads_in_positions(edited_positions_df, strand, group):

    unique_reads = set(chain.from_iterable(edited_positions_df["Reads"].str.split(",")))
    positions_per_read = {read: [] for read in unique_reads}

    overlapping_bases = 0

    alt_base = "G" if strand == "+" else "C"
    for x, position_row in enumerate(edited_positions_df.itertuples(), start=1):

        # deal with reads that were mapped to this pos
        mapped_bases = position_row.MappedBases
        mapped_reads = position_row.Reads.split(",")
        phred_scores = position_row.Phred

        last_phred_of_read_in_position = {}

        for phred_score, mapped_base, mapped_read in zip(
            phred_scores, mapped_bases, mapped_reads
        ):

            if mapped_base == alt_base:
                pos_edited_in_read = 1  # A2G/T2C editing
            elif mapped_base == ".":
                pos_edited_in_read = 0  # match
            else:
                # we ignore mismatches that aren't RNA editing by marking them as NaN
                pos_edited_in_read = -1

            # deal with overlapping paired-ended reads mapped to this position

            # if mapped_read not in last_phred_of_read_in_position:
            #     positions_per_read[mapped_read].append(pos_edited_in_read)
            #     last_phred_of_read_in_position[mapped_read] = phred_score
            # elif ord(phred_score) > ord(last_phred_of_read_in_position[mapped_read]):
            #     positions_per_read[mapped_read][-1] = pos_edited_in_read
            #     last_phred_of_read_in_position[mapped_read] = phred_score
            #     overlapping_bases += 1
            # else:
            #     raise Exception(
            #         f"{group = }: {x = }, {len(positions_per_read[mapped_read]) = }, {last_phred_of_read_in_position[mapped_read] = }"
            #     )
            if mapped_read not in last_phred_of_read_in_position:
                positions_per_read[mapped_read].append(pos_edited_in_read)
                last_phred_of_read_in_position[mapped_read] = phred_score
            else:
                if ord(phred_score) > ord(last_phred_of_read_in_position[mapped_read]):
                    positions_per_read[mapped_read][-1] = pos_edited_in_read
                    last_phred_of_read_in_position[mapped_read] = phred_score
                overlapping_bases += 1

        # deal with reads that weren't mapped to this pos
        unmapped_reads_in_pos = unique_reads - set(mapped_reads)
        for unmapped_read in unmapped_reads_in_pos:
            pos_edited_in_read = -1  #  i.e., we use this value as nan
            positions_per_read[unmapped_read].append(pos_edited_in_read)

        # check all reads got values for pos
        mapped_pos_dist = {len(bases) for bases in positions_per_read.values()}
        if len(mapped_pos_dist) != 1:
            raise Exception(
                f"Problem at line {x} in {group}'s positions_df: not all reads are mapped."
            )

    ic(group, overlapping_bases)

    return positions_per_read


def add_pacbio_se_multisample_reads_in_positions(edited_positions_df, strand, group):
    unique_reads = set(chain.from_iterable(edited_positions_df["Reads"].str.split(",")))
    positions_per_read = {read: [] for read in unique_reads}
    sample_per_read = {}

    alt_base = "G" if strand == "+" else "C"
    for x, position_row in enumerate(edited_positions_df.itertuples()):
        mapped_bases = position_row.MappedBases
        mapped_reads = position_row.Reads.split(",")
        mapped_samples = position_row.Samples.split(",")
        # deal with reads that were mapped to this pos
        for mapped_base, mapped_read, mapped_sample in zip(
            mapped_bases, mapped_reads, mapped_samples
        ):

            sample_per_read[mapped_read] = mapped_sample

            # 1 = A2G/T2C editing, 0 = match
            if mapped_base == alt_base:
                # pos_edited_in_read = "1" # A2G/T2C editing
                pos_edited_in_read = 1  # A2G/T2C editing
            elif mapped_base == ".":
                # pos_edited_in_read = "0" # match
                pos_edited_in_read = 0  # match
            else:
                # pos_edited_in_read = "-1" # we ignore mismatches that aren't RNA editing by marking them as NaN
                pos_edited_in_read = (
                    -1
                )  # we ignore mismatches that aren't RNA editing by marking them as NaN
            positions_per_read[mapped_read].append(pos_edited_in_read)
        # deal with reads that weren't mapped to this pos
        unmapped_reads_in_pos = unique_reads - set(mapped_reads)
        for mapped_read in unmapped_reads_in_pos:
            # pos_edited_in_read = "-1"
            pos_edited_in_read = -1
            positions_per_read[mapped_read].append(pos_edited_in_read)
        # check all reads got values for pos
        mapped_pos_dist = {len(bases) for bases in positions_per_read.values()}
        if len(mapped_pos_dist) != 1:
            raise Exception(
                f"Problem at line {x} in {group}'s positions_df: not all reads are mapped."
            )

    return positions_per_read, sample_per_read


def positions_to_reads(
    positions_file: Union[str, Path],
    sep: str,
    strand: str,
    group_col,
    group,
    parity: str,
):

    positions_df = pd.read_csv(positions_file, sep=sep)

    # retain only edited & reliable positions within coding regions
    edited_positions_df = positions_df.loc[
        positions_df["CDS"] & positions_df["Edited"] & ~positions_df["InProbRegion"]
    ]

    if parity in ["se", "SE"]:
        positions_per_read = add_pacbio_se_reads_in_positions(
            edited_positions_df, strand, group
        )
    elif parity in ["pe", "PE"]:
        positions_per_read = add_illumina_pe_reads_in_positions(
            edited_positions_df, strand, group
        )
    else:
        raise Exception(f"{group = }: Parity {parity} not recognized.")

    reads_df = pd.DataFrame(positions_per_read)
    reads_df = reads_df.T.reset_index().rename({"index": "Read"}, axis="columns")
    reads_df = reads_df.rename(
        columns={
            old_col: pos
            for old_col, pos in zip(
                reads_df.columns[1:], edited_positions_df["Position"]
            )
        }
    )
    reads_df.insert(0, group_col, group)
    annotate_edited_positions(reads_df, 2)
    annotate_unedited_positions(reads_df, 3)
    annotate_ambigous_positions(reads_df, 4)
    annotate_editing_frequency_per_read(reads_df, 2, positions_df, strand)
    return reads_df


def multisample_positions_to_reads(
    positions_file: Union[str, Path],
    sep: str,
    strand: str,
    group_col,
    group,
    parity: str,
    snp_noise_level: float,
    top_x_noisy_positions: int,
    pooled_transcript_noise_threshold: float,
):
    positions_df = pd.read_csv(positions_file, sep=sep)

    pooled_transcript_noise = (
        positions_df.loc[positions_df["Noise"] < snp_noise_level, "Noise"]
        .sort_values(ascending=False)[:top_x_noisy_positions]
        .mean()
    )
    if pooled_transcript_noise >= pooled_transcript_noise_threshold:
        raise ValueError(
            f"ERROR: {positions_file = } has {pooled_transcript_noise = } >= {pooled_transcript_noise_threshold = } "
            "and will not undergo fuether processing into reads, proteins, etc."
        )

    # retain only edited & reliable positions within coding regions
    edited_positions_df = positions_df.loc[
        positions_df["CDS"] & positions_df["Edited"] & ~positions_df["InProbRegion"]
    ]

    if parity in ["se", "SE"]:
        (
            positions_per_read,
            sample_per_read,
        ) = add_pacbio_se_multisample_reads_in_positions(
            edited_positions_df, strand, group
        )
    elif parity in ["pe", "PE"]:
        raise NotImplemented
    else:
        raise Exception(f"{group = }: Parity {parity} not recognized.")

    reads_df = pd.DataFrame(positions_per_read)
    reads_df = reads_df.T.reset_index().rename({"index": "Read"}, axis="columns")
    reads_df = reads_df.rename(
        columns={
            old_col: pos
            for old_col, pos in zip(
                reads_df.columns[1:], edited_positions_df["Position"]
            )
        }
    )
    reads_df.insert(0, group_col, group)
    reads_df.insert(1, "Sample", [sample_per_read[read] for read in reads_df["Read"]])
    annotate_edited_positions(reads_df, 3)
    annotate_unedited_positions(reads_df, 4)
    annotate_ambigous_positions(reads_df, 5)
    annotate_editing_frequency_per_read(reads_df, 3, positions_df, strand)

    return reads_df


def reads_to_unique_reads(reads_df, group_col, first_pos_loc=6):
    # data per read -> data per unique reads, and their supporting individual reads
    transcripts_df = reads_df.copy()  # todo no need to copy

    transcripts_df.insert(
        1,
        "Reads",
        (
            transcripts_df.groupby(transcripts_df.columns[first_pos_loc:].tolist())[
                "Read"
            ].transform(lambda x: ",".join(x))
        ),
    )
    del transcripts_df["Read"]

    transcripts_df.insert(
        2, "NumOfReads", transcripts_df["Reads"].apply(lambda x: len(x.split(",")))
    )
    first_pos_loc += 1

    transcripts_df = (
        transcripts_df.loc[
            ~transcripts_df.duplicated(subset=transcripts_df.columns[first_pos_loc:])
        ]
        .sort_values("NumOfReads", ascending=False)
        .reset_index(drop=True)
    )

    # transcripts_df.insert(
    #     1,
    #     "Transcript",
    #     pd.Series(
    #         transcripts_df[group_col] + "-" + transcripts_df.index.values.astype(str)
    #     ),
    # )
    transcripts_df.insert(
        1, "Transcript", serialize_compressed_ids(len(transcripts_df))
    )

    return transcripts_df


def multisample_reads_to_unique_reads(reads_df, first_pos_loc=7):
    # data per read and sample -> data per unique reads, and their supporting individual reads
    unique_reads_df = reads_df.copy()

    reads_grouped_by_positions = unique_reads_df.groupby(
        unique_reads_df.columns[first_pos_loc:].tolist()
    )

    unique_reads_df.insert(
        unique_reads_df.columns.get_loc("Read"),
        "Reads",
        reads_grouped_by_positions["Read"].transform(lambda x: ",".join(x)),
    )
    del unique_reads_df["Read"]

    unique_reads_df.insert(
        unique_reads_df.columns.get_loc("Sample"),
        "Samples",
        reads_grouped_by_positions["Sample"].transform(lambda x: ",".join(x)),
    )
    del unique_reads_df["Sample"]

    unique_reads_df.insert(
        unique_reads_df.columns.get_loc("Reads") + 1,
        "NumOfReads",
        unique_reads_df["Reads"].apply(lambda x: len(x.split(","))),
    )
    first_pos_loc += 1

    unique_reads_df = (
        unique_reads_df.drop_duplicates(subset=unique_reads_df.columns[first_pos_loc:])
        .sort_values("NumOfReads", ascending=False)
        .reset_index(drop=True)
    )

    assert unique_reads_df["NumOfReads"].sum() == len(reads_df)

    unique_reads_df.insert(
        1, "UniqueRead", serialize_compressed_ids(len(unique_reads_df))
    )

    return unique_reads_df


# def reads_and_unique_reads(
#     positions_df,
#     strand,
#     group_col,
#     group,
#     parity: str,
#     reads_out_file: Union[Path, str, None] = None,
#     unique_reads_out_file: Union[Path, str, None] = None,
#     out_files_sep: str = "\t",
# ) -> tuple[pd.DataFrame, pd.DataFrame]:

#     reads_df = positions_to_reads(positions_df, strand, group_col, group, parity)
#     unique_reads_df = reads_to_unique_reads(reads_df, group_col)

#     if reads_out_file:
#         reads_df.to_csv(reads_out_file, sep=out_files_sep, index=False, na_rep=np.NaN)

#     if unique_reads_out_file:
#         unique_reads_df.to_csv(
#             unique_reads_out_file, sep=out_files_sep, index=False, na_rep=np.NaN
#         )

#     return reads_df, unique_reads_df


def reads_and_unique_reads(
    positions_file,
    strand,
    group_col,
    group,
    parity: str,
    reads_out_file: Union[Path, str, None] = None,
    unique_reads_out_file: Union[Path, str, None] = None,
    sep: str = "\t",
) -> tuple[pd.DataFrame, pd.DataFrame]:

    if not Path(positions_file).exists():
        return

    try:
        reads_df = multisample_positions_to_reads(
            positions_file, sep, strand, group_col, group, parity
        )
        unique_reads_df = reads_to_unique_reads(reads_df, group_col)

        if reads_out_file:
            reads_df.to_csv(reads_out_file, sep=sep, index=False, na_rep=np.NaN)

        if unique_reads_out_file:
            unique_reads_df.to_csv(
                unique_reads_out_file, sep=sep, index=False, na_rep=np.NaN
            )
    except ValueError:  # will be raised if len(edited_positions_df) == 0
        pass

    # return reads_df, unique_reads_df


def multisample_reads_and_unique_reads(
    positions_file,
    strand,
    group_col,
    group,
    parity: str,
    snp_noise_level: float,
    top_x_noisy_positions: int,
    pooled_transcript_noise_threshold: float,
    reads_out_file: Union[Path, str, None] = None,
    unique_reads_out_file: Union[Path, str, None] = None,
    sep: str = "\t",
) -> tuple[pd.DataFrame, pd.DataFrame]:

    if not Path(positions_file).exists():
        return

    try:
        reads_df = multisample_positions_to_reads(
            positions_file,
            sep,
            strand,
            group_col,
            group,
            parity,
            snp_noise_level,
            top_x_noisy_positions,
            pooled_transcript_noise_threshold,
        )

        unique_reads_df = multisample_reads_to_unique_reads(reads_df)

        if reads_out_file:
            reads_df.to_csv(reads_out_file, sep=sep, index=False, na_rep=np.NaN)

        if unique_reads_out_file:
            unique_reads_df.to_csv(
                unique_reads_out_file, sep=sep, index=False, na_rep=np.NaN
            )
    except ValueError as e:  # will be raised if pooled_transcript_noise_threshold is met
        print(e)
