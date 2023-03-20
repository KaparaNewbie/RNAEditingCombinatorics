from typing import Union
from collections import defaultdict
from pathlib import Path

import pandas as pd
import numpy as np

from EditingUtils.seq import make_fasta_dict, AA_BY_TRIPLET
from Pileup.ids import serialize_compressed_ids

StrOrNoneOrInt = Union[str, None, int]


def filter_triplets(possible_triplets, base, pos):
    if base == -1:
        return possible_triplets
    else:
        return [triplet for triplet in possible_triplets if triplet[pos] == base]


# def guided_translation(
#     x1: str,
#     y1: str,
#     z1: str,
#     x2: StrOrNoneOrInt,
#     y2: StrOrNoneOrInt,
#     z2: StrOrNoneOrInt,
# ):
#     """
#     Translate a triplet (x2, y2, z2) to an amino acid.

#     The (x1, y1, z1) triplet is the original triplet from the transcriptome.

#     The (x2, y2, z2) triplet is a corresponding triplet of a read.

#     Use the original triplet as a template in case any of the bases in the new triplet is None.
#     Return a set of all possible translations of amino acids in case any of the new bases equals -1.
#     """
#     if x2 is None:
#         x2 = x1
#     if y2 is None:
#         y2 = y1
#     if z2 is None:
#         z2 = z1

#     possible_triplets = list(AA_BY_TRIPLET.keys())

#     for base, pos in zip([x2, y2, z2], [0, 1, 2]):
#         possible_triplets = filter_triplets(possible_triplets, base, pos)

#     amino_acids = {AA_BY_TRIPLET[triplet] for triplet in possible_triplets}
#     amino_acids = sorted(list(amino_acids))
#     amino_acids = ",".join(amino_acids)

#     return amino_acids


def guided_translation(
    x1: str,
    y1: str,
    z1: str,
    x2: StrOrNoneOrInt,
    y2: StrOrNoneOrInt,
    z2: StrOrNoneOrInt,
):
    """
    Return all possible translations of the triplet (x2, y2, z2).

    The (x1, y1, z1) triplet is the original reference triplet.
    The (x2, y2, z2) triplet is a corresponding triplet of a read.

    Use the original triplet as a template in case any of the bases in the new triplet is None.
    This is the case for each reference base other than A.

    When an A base is marked as -1, it could be either A or G in the read.
    """
    if x2 is None:
        x2 = x1
    if y2 is None:
        y2 = y1
    if z2 is None:
        z2 = z1

    possible_triplets = list(AA_BY_TRIPLET.keys())

    # filter possible triplets according to the base in each position
    for base, pos in zip([x2, y2, z2], [0, 1, 2]):
        if base == -1:
            possible_bases_in_pos = ["A", "G"]
        else:
            possible_bases_in_pos = [base]
        possible_triplets = [
            triplet
            for triplet in possible_triplets
            if triplet[pos] in possible_bases_in_pos
        ]

    amino_acids = {AA_BY_TRIPLET[triplet] for triplet in possible_triplets}
    amino_acids = sorted(list(amino_acids))
    amino_acids = ",".join(amino_acids)

    return amino_acids


def annotate_min_max_non_syns(
    df, transalted_orf_length, original_aas, first_aa_pos, new_cols_first_pos=None
):
    def min_non_syns_in_row(row):
        return sum(
            [1 for aa, original_aa in zip(row, original_aas) if original_aa not in aa]
        )

    def max_non_syns_in_row(row):
        return sum(
            [1 for aa, original_aa in zip(row, original_aas) if aa != original_aa]
        )

    min_non_syns = df.iloc[:, first_aa_pos:].apply(min_non_syns_in_row, axis="columns")
    max_non_syns = df.iloc[:, first_aa_pos:].apply(max_non_syns_in_row, axis="columns")

    assert max_non_syns.ge(min_non_syns).all()

    min_non_syns_frequency = min_non_syns / transalted_orf_length
    max_non_syns_frequency = max_non_syns / transalted_orf_length

    new_cols_first_pos = new_cols_first_pos if new_cols_first_pos else first_aa_pos
    df.insert(new_cols_first_pos, "MinNonSyns", min_non_syns)
    df.insert(new_cols_first_pos + 1, "MaxNonSyns", max_non_syns)
    df.insert(new_cols_first_pos + 2, "MinNonSynsFrequency", min_non_syns_frequency)
    df.insert(new_cols_first_pos + 3, "MaxNonSynsFrequency", max_non_syns_frequency)


def unique_reads_to_proteins(
    unique_reads_file: Union[str, Path],
    sep: str,
    region: str,
    start: int,
    end: int,
    strand: str,
    transcriptome: Union[str, Path],
    first_pos_loc: int = 8,
):
    unique_reads_df = pd.read_csv(unique_reads_file, sep=sep)

    # really shouldn't happen, that bug has been fixed by not creating ids with "NA*"
    if len(unique_reads_df.loc[unique_reads_df["Transcript"].isna()]) != 0:
        raise Exception("NaN values in `Transcripts` col of `unique_reads_df`!")

    unique_reads_df = unique_reads_df.rename(
        columns={col: int(col) for col in unique_reads_df.columns[first_pos_loc:]}
    )

    positions_cols = set(unique_reads_df.columns[first_pos_loc:])

    # if strand == "+":
    #     ref_base, alt_base = "A", "G"
    # else:
    #     ref_base, alt_base = "T", "C"  # get complementary bases
    # position_to_base = {-1: -1, 0: ref_base, 1: alt_base}
    position_to_base = {-1: -1, 0: "A", 1: "G"}

    transcriptome_dict = make_fasta_dict(transcriptome)
    m_rna = transcriptome_dict[region][start:end]
    original_triplets = [m_rna[x : x + 3] for x in range(0, len(m_rna), 3)]
    if strand == "-":
        original_triplets = [
            triplet.reverse_complement() for triplet in original_triplets
        ]
    assert len(original_triplets[0]) == len(original_triplets[-1]) == 3
    transalted_orf_length = len(original_triplets)

    new_triplets_positions = [(x, x + 1, x + 2) for x in range(start, end, 3)]

    transcripts_dict = unique_reads_df.to_dict(orient="records")

    aas_per_unique_read = defaultdict(list)  # amino acids per unique read
    triplets_cols = []
    original_aas = []

    for original_triplet, new_triplet_positions in zip(
        original_triplets, new_triplets_positions
    ):
        # the original triplet's bases according to the reference
        x1, y1, z1 = list(original_triplet)
        # corresponding indices - to be used by the reads
        pos_x2, pos_y2, pos_z2 = new_triplet_positions
        pos_x2_exist = pos_x2 in positions_cols
        pos_y2_exist = pos_y2 in positions_cols
        pos_z2_exist = pos_z2 in positions_cols
        # check wether there are reads with any of theses positions edited
        # if not - continue to the next triplet
        if not any([pos_x2_exist, pos_y2_exist, pos_z2_exist]):
            continue
        original_aa = AA_BY_TRIPLET[original_triplet]
        triplets_cols.append(f"{pos_x2}:{pos_x2 + 3}({original_aa})")
        original_aas.append(original_aa)
        for row in transcripts_dict:
            row_id = row["Transcript"]
            x2 = position_to_base[row[pos_x2]] if pos_x2_exist else None
            y2 = position_to_base[row[pos_y2]] if pos_y2_exist else None
            z2 = position_to_base[row[pos_z2]] if pos_z2_exist else None
            if strand == "-":
                # reverse the bases (they're already complementary due to position_to_base)
                x2, z2 = z2, x2
            # aa = guided_translation(x1, y1, z1, x2, y2, z2)
            aa = guided_translation(x1, y1, z1, x2, y2, z2)
            aas_per_unique_read[row_id].append(aa)

    proteins_df = pd.DataFrame(aas_per_unique_read)
    proteins_df = proteins_df.T.reset_index().rename(
        {"index": "Transcript"}, axis="columns"
    )
    proteins_df = proteins_df.set_axis(["Transcript"] + triplets_cols, axis="columns")
    proteins_df = unique_reads_df.iloc[:, :first_pos_loc].merge(
        proteins_df, on="Transcript"
    )
    annotate_min_max_non_syns(
        proteins_df, transalted_orf_length, original_aas, first_pos_loc, 1
    )
    return proteins_df


def multisample_unique_reads_to_proteins(
    unique_reads_file: Union[str, Path],
    sep: str,
    region: str,
    start: int,
    end: int,
    strand: str,
    transcriptome: Union[str, Path],
    first_pos_loc: int = 9,
):
    unique_reads_df = pd.read_csv(unique_reads_file, sep=sep, dtype={"UniqueRead": str})

    # really shouldn't happen, that bug has been fixed by not creating ids with "NA*"
    if len(unique_reads_df.loc[unique_reads_df["UniqueRead"].isna()]) != 0:
        raise Exception("NaN values in `UniqueRead` col of `unique_reads_df`!")

    unique_reads_df = unique_reads_df.rename(
        columns={col: int(col) for col in unique_reads_df.columns[first_pos_loc:]}
    )

    positions_cols = set(unique_reads_df.columns[first_pos_loc:])

    position_to_base = {-1: -1, 0: "A", 1: "G"}

    transcriptome_dict = make_fasta_dict(transcriptome)
    m_rna = transcriptome_dict[region][start:end]
    original_triplets = [m_rna[x : x + 3] for x in range(0, len(m_rna), 3)]
    if strand == "-":
        original_triplets = [
            triplet.reverse_complement() for triplet in original_triplets
        ]
    assert len(original_triplets[0]) == len(original_triplets[-1]) == 3
    transalted_orf_length = len(original_triplets)

    new_triplets_positions = [(x, x + 1, x + 2) for x in range(start, end, 3)]

    unique_reads_dict = unique_reads_df.to_dict(orient="records")

    aas_per_unique_read = defaultdict(list)  # amino acids per unique read
    triplets_cols = []
    original_aas = []

    for original_triplet, new_triplet_positions in zip(
        original_triplets, new_triplets_positions
    ):
        # the original triplet's bases according to the reference
        x1, y1, z1 = list(original_triplet)
        # corresponding indices - to be used by the reads
        pos_x2, pos_y2, pos_z2 = new_triplet_positions
        pos_x2_exist = pos_x2 in positions_cols
        pos_y2_exist = pos_y2 in positions_cols
        pos_z2_exist = pos_z2 in positions_cols
        # check wether there are reads with any of theses positions edited
        # if not - continue to the next triplet
        if not any([pos_x2_exist, pos_y2_exist, pos_z2_exist]):
            continue
        original_aa = AA_BY_TRIPLET[original_triplet]
        triplets_cols.append(f"{pos_x2}:{pos_x2 + 3}({original_aa})")
        original_aas.append(original_aa)
        for row in unique_reads_dict:
            row_id = row["UniqueRead"]
            x2 = position_to_base[row[pos_x2]] if pos_x2_exist else None
            y2 = position_to_base[row[pos_y2]] if pos_y2_exist else None
            z2 = position_to_base[row[pos_z2]] if pos_z2_exist else None
            if strand == "-":
                # reverse the bases (they're already complementary due to position_to_base)
                x2, z2 = z2, x2
            # aa = guided_translation(x1, y1, z1, x2, y2, z2)
            aa = guided_translation(x1, y1, z1, x2, y2, z2)
            aas_per_unique_read[row_id].append(aa)

    proteins_df = pd.DataFrame(aas_per_unique_read)
    proteins_df = proteins_df.T.reset_index().rename(
        {"index": "UniqueRead"}, axis="columns"
    )
    proteins_df = proteins_df.set_axis(["UniqueRead"] + triplets_cols, axis="columns")
    proteins_df = unique_reads_df.iloc[:, :first_pos_loc].merge(
        proteins_df, on="UniqueRead"
    )
    annotate_min_max_non_syns(
        proteins_df, transalted_orf_length, original_aas, first_pos_loc, 1
    )

    return proteins_df


def proteins_to_unique_proteins(
    proteins_df: pd.DataFrame, group_col: str, first_pos_loc: int = 12
):
    unique_proteins_df = proteins_df.drop("NumOfReads", axis=1)
    first_pos_loc -= 1

    unique_proteins_df.insert(
        unique_proteins_df.columns.get_loc("Transcript"),
        "Transcripts",
        (
            unique_proteins_df.groupby(
                unique_proteins_df.columns[first_pos_loc:].tolist()
            )["Transcript"].transform(lambda x: ",".join(x))
        ),
    )
    del unique_proteins_df["Transcript"]

    unique_proteins_df.insert(
        unique_proteins_df.columns.get_loc("Transcripts") + 1,
        "NumOfTranscripts",
        unique_proteins_df["Transcripts"].apply(lambda x: len(x.split(","))),
    )
    first_pos_loc += 1
    unique_proteins_df.insert(
        unique_proteins_df.columns.get_loc("Reads"),
        "Reads2",
        (
            unique_proteins_df.groupby(
                unique_proteins_df.columns[first_pos_loc:].tolist()
            )["Reads"].transform(lambda x: ",".join(x))
        ),
    )
    unique_proteins_df = unique_proteins_df.drop("Reads", axis="columns").rename(
        columns={"Reads2": "Reads"}
    )
    unique_proteins_df.insert(
        unique_proteins_df.columns.get_loc("Reads") + 1,
        "NumOfReads",
        unique_proteins_df["Reads"].apply(lambda x: len(x.split(","))),
    )
    first_pos_loc += 1

    unique_proteins_df = (
        unique_proteins_df.loc[
            ~unique_proteins_df.duplicated(
                subset=unique_proteins_df.columns[first_pos_loc:]
            )
        ]
        .sort_values("NumOfReads", ascending=False)
        .reset_index(drop=True)
    )
    assert proteins_df["NumOfReads"].sum() == unique_proteins_df["NumOfReads"].sum()
    assert sum(proteins_df["NumOfReads"] * proteins_df["MinNonSyns"]) == sum(
        unique_proteins_df["NumOfReads"] * unique_proteins_df["MinNonSyns"]
    )
    assert sum(proteins_df["NumOfReads"] * proteins_df["MaxNonSyns"]) == sum(
        unique_proteins_df["NumOfReads"] * unique_proteins_df["MaxNonSyns"]
    )

    unique_proteins_df.insert(
        1, "Protein", serialize_compressed_ids(len(unique_proteins_df))
    )

    return unique_proteins_df


def multisample_proteins_to_unique_proteins(
    proteins_df: pd.DataFrame, group_col: str, first_pos_loc: int = 13
):
    unique_proteins_df = proteins_df.copy()

    proteins_grouped_by_positions = unique_proteins_df.groupby(
        unique_proteins_df.columns[first_pos_loc:].tolist()
    )

    # join unique reads
    unique_proteins_df.insert(
        unique_proteins_df.columns.get_loc("UniqueRead"),
        "UniqueReads",
        proteins_grouped_by_positions["UniqueRead"].transform(lambda x: ",".join(x)),
    )
    # drop the col of individual unique reads
    unique_proteins_df = unique_proteins_df.drop(["UniqueRead"], axis=1)
    # count the num of unique reads
    unique_proteins_df.insert(
        unique_proteins_df.columns.get_loc("UniqueReads") + 1,
        "NumOfUniqueReads",
        unique_proteins_df["UniqueReads"].apply(lambda x: len(x.split(","))),
    )
    first_pos_loc += 1

    # join samples
    unique_proteins_df["Samples"] = proteins_grouped_by_positions["Samples"].transform(
        lambda x: ",".join(x)
    )

    # join reads
    unique_proteins_df["Reads"] = proteins_grouped_by_positions["Reads"].transform(
        lambda x: ",".join(x)
    )
    # sum joined reads
    unique_proteins_df["NumOfReads"] = proteins_grouped_by_positions[
        "NumOfReads"
    ].transform(sum)

    unique_proteins_df = (
        unique_proteins_df.drop_duplicates(
            subset=unique_proteins_df.columns[first_pos_loc:]
        )
        .sort_values("NumOfReads", ascending=False)
        .reset_index(drop=True)
    )

    assert proteins_df["NumOfReads"].sum() == unique_proteins_df["NumOfReads"].sum()
    assert sum(proteins_df["NumOfReads"] * proteins_df["MinNonSyns"]) == sum(
        unique_proteins_df["NumOfReads"] * unique_proteins_df["MinNonSyns"]
    )
    assert sum(proteins_df["NumOfReads"] * proteins_df["MaxNonSyns"]) == sum(
        unique_proteins_df["NumOfReads"] * unique_proteins_df["MaxNonSyns"]
    )

    unique_proteins_df.insert(
        1, "Protein", serialize_compressed_ids(len(unique_proteins_df))
    )

    return unique_proteins_df


def proteins_and_unique_proteins(
    unique_reads_file: Union[str, Path],
    region: str,
    start: int,
    end: int,
    strand: str,
    transcriptome: Union[str, Path],
    group_col: str,
    proteins_out_file: Union[Path, str, None] = None,
    unique_proteins_out_file: Union[Path, str, None] = None,
    sep: str = "\t",
):
    if not Path(unique_reads_file).exists():
        return

    proteins_df = unique_reads_to_proteins(
        unique_reads_file, sep, region, start, end, strand, transcriptome
    )
    unique_proteins_df = proteins_to_unique_proteins(proteins_df, group_col)

    if proteins_out_file:
        proteins_df.to_csv(proteins_out_file, sep=sep, index=False, na_rep=np.NaN)
    if unique_proteins_out_file:
        unique_proteins_df.to_csv(
            unique_proteins_out_file, sep=sep, index=False, na_rep=np.NaN
        )

    return proteins_df, unique_proteins_df


def multisample_proteins_and_unique_proteins(
    unique_reads_file: Union[str, Path],
    region: str,
    start: int,
    end: int,
    strand: str,
    transcriptome: Union[str, Path],
    group_col: str,
    proteins_out_file: Union[Path, str, None] = None,
    unique_proteins_out_file: Union[Path, str, None] = None,
    sep: str = "\t",
):
    if not Path(unique_reads_file).exists():
        return

    proteins_df = multisample_unique_reads_to_proteins(
        unique_reads_file, sep, region, start, end, strand, transcriptome
    )
    unique_proteins_df = multisample_proteins_to_unique_proteins(proteins_df, group_col)

    if proteins_out_file:
        proteins_df.to_csv(proteins_out_file, sep=sep, index=False, na_rep=np.NaN)
    if unique_proteins_out_file:
        unique_proteins_df.to_csv(
            unique_proteins_out_file, sep=sep, index=False, na_rep=np.NaN
        )

    return proteins_df, unique_proteins_df
