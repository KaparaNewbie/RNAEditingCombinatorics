# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 14:55:24 2020

@author: shosh
"""
import argparse

import pandas as pd
from Bio import SeqIO

# from .General.argparse_utils import abs_path_from_str

# trinity_file = "/private7/projects/Combinatorics/O.bim/Annotations/orfs_bim.fa"


def read_trinity_mrna_files(trinity_file):
    """
    read relevat information from transcriptome file
    """
    data = []
    col_names = [
        "id",
        "protein",
        "protein_name",
        "orfs_start",
        "orfs_end",
        "strand",
        "sequence",
    ]

    # records = list(SeqIO.parse(open(trinity_file, "r"), "fasta"))
    # record = records[0]

    for record in SeqIO.parse(open(trinity_file, "r"), "fasta"):
        rec_data = record.description.split("\t")
        if (
            rec_data[0][-1] == " "
        ):  # some fastafiles have spaces after each id, so fixing it here.
            rec_data[0] = rec_data[0].replace(" ", "")
        protein = (
            rec_data[-1].split("|")[2].split(" ")[0]
        )  # reading proteing from description assuming it was added to header using the transcriptome built pipeline we have for trinity
        rec_data = (
            rec_data[0],
            protein,
            protein.split("_")[0],
            int(rec_data[2]),
            int(rec_data[4]),
            rec_data[6],
            record.seq,
        )
        data.append(rec_data)

    df = pd.DataFrame(data=data, columns=col_names)

    return df


def create_gtf_table(trinity_df):
    data = []
    columns = [
        "chr",
        "protein",
        "region",
        "start",
        "end",
        "score",
        "strand",
        "reading_frame",
        "gene_id_to_transcript_id",
    ]

    for i, row in trinity_df.iterrows():
        if row["orfs_start"] != 1:
            data.append(
                (
                    row["id"],
                    row["protein"],
                    "UTR",
                    1,
                    row["orfs_start"] - 1,
                    0,
                    row["strand"],
                    0,
                    'gene_id "' + row["id"] + '"; transcript_id "' + row["id"] + '";',
                )
            )
            data.append(
                (
                    row["id"],
                    row["protein"],
                    "CDS",
                    row["orfs_start"],
                    row["orfs_end"],
                    0,
                    row["strand"],
                    0,
                    'gene_id "' + row["id"] + '"; transcript_id "' + row["id"] + '";',
                )
            )
            if row["orfs_end"] != len(row["sequence"]):
                data.append(
                    (
                        row["id"],
                        row["protein"],
                        "UTR",
                        row["orfs_end"] + 1,
                        len(row["sequence"]),
                        0,
                        row["strand"],
                        0,
                        'gene_id "'
                        + row["id"]
                        + '"; transcript_id "'
                        + row["id"]
                        + '";',
                    )
                )
        else:
            data.append(
                (
                    row["id"],
                    row["protein"],
                    "CDS",
                    row["orfs_start"],
                    row["orfs_end"],
                    0,
                    row["strand"],
                    0,
                    'gene_id "' + row["id"] + '"; transcript_id "' + row["id"] + '";',
                )
            )
            if row["orfs_end"] != len(row["sequence"]):
                data.append(
                    (
                        row["id"],
                        row["protein"],
                        "UTR",
                        row["orfs_end"] + 1,
                        len(row["sequence"]),
                        0,
                        row["strand"],
                        0,
                        'gene_id "'
                        + row["id"]
                        + '"; transcript_id "'
                        + row["id"]
                        + '";',
                    )
                )

    return pd.DataFrame(data=data, columns=columns)


def main(trinity_file, output_path):
    trinity_df = read_trinity_mrna_files(trinity_file)
    gtf_df = create_gtf_table(trinity_df)
    gtf_df.to_csv(output_path, header=False, index=False, sep="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="prepare salmon gtf input from transcriptome"
    )
    parser.add_argument(
        "--trinity_file",
        help="path to trinity transcriptome file",
        # type=abs_path_from_str,
        required=True,
    )
    parser.add_argument(
        "--output_path",
        help="path to output gtf file",
        # type=abs_path_from_str,
        required=True,
    )
    args = parser.parse_args()
    main(args.trinity_file, args.output_path)