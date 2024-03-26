from pathlib import Path
import sys
from multiprocessing import Pool
import gzip
import argparse

import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from icecream import ic

sys.path.append(str(Path("Code").absolute()))
from EditingUtils.seq import make_fasta_dict
from General.argparse_utils import abs_path_from_str
from General.consts import final_words, ic_prefix

# configure icecream to print the time of the print and the context (file, line, function)
ic.configureOutput(includeContext=True, prefix=ic_prefix)


def decode_degenerate_base(
    degenerate_code: int, original_base_location: int, transcript_seq: Seq, chrom: str
):
    try:
        original_base = transcript_seq[original_base_location]
    except IndexError as e:
        raise IndexError(
            f"IndexError: index {original_base_location} out of range in {chrom} "
        )
    if original_base != "A":
        return original_base
    if degenerate_code == 0:
        return "A"
    elif degenerate_code == 1:
        return "G"
    elif degenerate_code == -1:
        return "X"
    else:
        raise ValueError(f"Invalid degenerate code {degenerate_code}")


def degenerate_row_to_seqrecord(
    degenerate_row: pd.Series,
    transcript_seq: Seq,
    chrom: str,
    start: int,
    end: int,
    strand: str,
    name: str,
):
    read_id = degenerate_row.name
    decoded_degenrate_read_bases = [
        decode_degenerate_base(
            degenerate_code, original_base_location, transcript_seq, chrom
        )
        for original_base_location, degenerate_code in degenerate_row.items()
    ]
    read_seq = Seq("".join(decoded_degenrate_read_bases))
    read_record = SeqRecord(
        read_seq, id=read_id, description=f"{chrom}:{start}-{end}({strand}) {name}"
    )
    return read_record


def degenerate_reads_in_sample(
    reads_file: Path,
    chrom: str,
    start: int,
    end: int,
    strand: str,
    name: str,
    degenerate_reads_out_file: Path,
    sep: str,
    reads_first_col_pos: int,
    processes: int,
    transcriptome_dict: dict[str, Seq],
):
    ic(
        reads_file,
        chrom,
        start,
        end,
        strand,
        name,
        degenerate_reads_out_file,
        sep,
        reads_first_col_pos,
        processes,
    )
    reads_df = pd.read_csv(reads_file, sep=sep, dtype={"Read": str})
    reads_df = reads_df.set_index("Read")
    reads_first_col_pos -= 1  # since we used the "Read" column as the index, we need to subtract 1 from reads_first_col_pos
    reads_df = reads_df.iloc[:, reads_first_col_pos:]
    reads_df = reads_df.rename(columns=lambda x: int(x))

    rows = len(reads_df)
    cols = end - start

    degenerate_reads_df = pd.DataFrame(
        np.zeros((rows, cols), dtype=int),
        columns=range(start, end),
        index=reads_df.index,
    )
    degenerate_reads_df.update(reads_df)

    # transcript_seq = transcriptome_dict[chrom][start:end]
    transcript_seq = transcriptome_dict[chrom]

    starmap_input = []
    for x in range(rows):
        try:
            degenerate_row = degenerate_reads_df.iloc[x, :]
            starmap_input.append(
                (degenerate_row, transcript_seq, chrom, start, end, strand, name)
            )
        except IndexError as e:
            ic(x, rows, degenerate_reads_df.shape)
            raise e

    with Pool(processes=processes) as pool:
        read_records = pool.starmap(
            func=degenerate_row_to_seqrecord, iterable=starmap_input
        )

    if degenerate_reads_out_file.name.endswith(".gz"):
        # Open the gzipped FASTA file in write mode
        with gzip.open(degenerate_reads_out_file, "wt") as handle:
            # Write the sequence records to the gzipped FASTA file
            SeqIO.write(read_records, handle, "fasta")
    else:
        SeqIO.write(read_records, degenerate_reads_out_file, "fasta")


def directly_given_variables_main(
    *,
    reads_files,
    chroms,
    starts,
    ends,
    strands,
    names,
    out_dir,
    sep,
    reads_first_col_pos,
    transcriptome_file,
    processes,
    postfix_to_remove,
    postfix_to_add,
    **kwargs,
):
    ic(
        reads_files,
        chroms,
        starts,
        ends,
        strands,
        names,
        out_dir,
        sep,
        reads_first_col_pos,
        transcriptome_file,
        processes,
        postfix_to_remove,
        postfix_to_add,
    )
    transcriptome_dict = make_fasta_dict(transcriptome_file)
    # out_dir = degenerate_reads_out_file.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    for reads_file, chrom, start, end, strand, name in zip(
        reads_files, chroms, starts, ends, strands, names
    ):

        # postfix_to_remove = ".csv.gz"
        # postfix_to_add = ".degenerate.fa.gz"
        degenerate_reads_out_file = Path(
            out_dir, Path(reads_file).name.replace(postfix_to_remove, postfix_to_add)
        )

        try:
            degenerate_reads_in_sample(
                reads_file,
                chrom,
                start,
                end,
                strand,
                name,
                degenerate_reads_out_file,
                sep,
                reads_first_col_pos,
                processes,
                transcriptome_dict,
            )
        except IndexError as e:
            ic(e)
            ic(reads_file, chrom, start, end, strand, name)
            continue


def undirectly_given_variables_main(
    *,
    cds_regions,
    sep,
    in_dir,
    postfix_to_find,
    prefix_to_chrom,
    postfix_to_remove,
    postfix_to_add,
    reads_first_col_pos,
    transcriptome_file,
    processes,
    out_dir,
    **kwargs,
):
    ic(
        cds_regions,
        sep,
        in_dir,
        postfix_to_find,
        prefix_to_chrom,
        postfix_to_remove,
        postfix_to_add,
        reads_first_col_pos,
        transcriptome_file,
        processes,
        out_dir,
    )
    reads_files = list(in_dir.glob(f"*{postfix_to_find}"))
    # ic(reads_files)
    # len(reads_files)
    chroms = [
        reads_file.name.replace(prefix_to_chrom, "").split(".")[0]
        for reads_file in reads_files
    ]
    cds_df = pd.read_csv(
        cds_regions,
        sep=sep,
        names="Chrom Start End Name Score Strand".split(),
        comment="#",
    )
    cds_df = cds_df.loc[cds_df["Chrom"].isin(chroms)]

    reads_file_and_chroms_df = pd.DataFrame({"ReadsFile": reads_files, "Chrom": chroms})
    cds_df = cds_df.merge(reads_file_and_chroms_df, on="Chrom", how="left")

    reads_files = cds_df["ReadsFile"].tolist()
    chroms = cds_df["Chrom"].tolist()
    starts = cds_df["Start"].tolist()
    ends = cds_df["End"].tolist()
    strands = cds_df["Strand"].tolist()
    names = cds_df["Name"].tolist()

    directly_given_variables_main(
        reads_files=reads_files,
        chroms=chroms,
        starts=starts,
        ends=ends,
        strands=strands,
        names=names,
        out_dir=out_dir,
        sep=sep,
        reads_first_col_pos=reads_first_col_pos,
        transcriptome_file=transcriptome_file,
        processes=processes,
        postfix_to_remove=postfix_to_remove,
        postfix_to_add=postfix_to_add,
    )


def define_args() -> argparse.Namespace:
    # create parser
    # description = "TBD description"
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        # description=description
    )
    subparsers = parser.add_subparsers(help="sub-command help")

    # define args
    parser.add_argument(
        "--transcriptome_file",
        required=True,
        type=abs_path_from_str,
        help="A transcritpome reference fasta file. Should have no introns.",
    )
    parser.add_argument("--out_dir", type=abs_path_from_str, required=True)
    parser.add_argument(
        "--processes",
        type=int,
        default=6,
        help="Maximal number of processes to run in parallel.",
    )
    parser.add_argument(
        "--sep", default="\t", help="Delimiter for input tabular files."
    )
    parser.add_argument(
        "--reads_first_col_pos",
        required=True,
        type=int,
        help="The first column of the reads file that contains editing information, 0-based.",
    )
    parser.add_argument(
        "--postfix_to_remove",
        required=True,
        help="Remove this part from the reads file name and replace it with `--postfix_to_add`.",
    )
    parser.add_argument(
        "--postfix_to_add",
        default=".degenerate.fa.gz",
        help="Add this postfix to the degenerate out fasta files.",
    )

    directly_given_variables_subparser = subparsers.add_parser(
        "directly_given_variables",
        # help="Original use of this program by suplying a data table",
    )
    directly_given_variables_subparser.set_defaults(func=directly_given_variables_main)
    directly_given_variables_subparser.add_argument(
        "--reads_files",
        required=True,
        nargs="+",
        help="Paths to the reads files.",
        type=abs_path_from_str,
    )
    directly_given_variables_subparser.add_argument(
        "--chroms", required=True, nargs="+", help="Chromosomes of the reads files."
    )
    directly_given_variables_subparser.add_argument(
        "--starts",
        required=True,
        nargs="+",
        help="Start coordinates of the reads files' transcripts, 0-based, inclusive.",
        type=int,
    )
    directly_given_variables_subparser.add_argument(
        "--ends",
        required=True,
        nargs="+",
        help="End coordinates of the reads files' transcripts, 1-based, exclusive.",
        type=int,
    )
    directly_given_variables_subparser.add_argument(
        "--strands",
        required=True,
        nargs="+",
        help="Strands of the reads files' transcripts.",
    )
    directly_given_variables_subparser.add_argument(
        "--names",
        required=True,
        nargs="+",
        help="Protein names of the reads files' transcripts.",
    )

    undirectly_given_variables_subparser = subparsers.add_parser(
        "undirectly_given_variables",
        # help="Original use of this program by suplying a data table",
    )
    undirectly_given_variables_subparser.set_defaults(
        func=undirectly_given_variables_main
    )
    undirectly_given_variables_subparser.add_argument(
        "--in_dir",
        type=abs_path_from_str,
        help="Search for reads files in this directory.",
    )
    undirectly_given_variables_subparser.add_argument(
        "--cds_regions",
        type=abs_path_from_str,
        help="6-col bed file of coding regions in the transcriptome file. ",
    )
    undirectly_given_variables_subparser.add_argument(
        "--postfix_to_find",
        required=True,
        help="Find reads files with this postfix in `in_dir`.",
    )
    undirectly_given_variables_subparser.add_argument(
        "--prefix_to_chrom",
        required=True,
        help="In each reads file name, the chromosome name is between this prefix to its left, and the first dot to its right.",
    )

    return parser


if __name__ == "__main__":
    parser = define_args()
    args = parser.parse_args()
    ic(args)
    args.func(
        **vars(args)
    )  # https://stackoverflow.com/a/35824590/10249633 argparse.Namespace -> Dict
    final_words()


# in_dir = Path(
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina"
# )
# postfix_to_find = ".reads.csv"
# prefix_to_chrom = "reads.sorted.aligned.filtered."

# reads_dir = Path(
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30"
# )
# reads_postfix = "C0x1291.aligned.sorted.MinRQ998.reads.csv.gz"
# cds_regions = "/private7/projects/Combinatorics/D.pealeii/Annotations/orfs_squ.bed"


# reads_file = Path(
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.reads.csv.gz"
# )
# chrom = "comp141882_c0_seq14"
# start = 0
# end = 6294
# strand = "+"
# name = "PCLO_CHICK"

# degenerate_reads_out_file = Path(
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.reads.degenerate.fa.gz"
# )

# sep = "\t"
# reads_first_col_pos = 6

# processes = 10

# transcriptome_file = Path(
#     "/private7/projects/Combinatorics/D.pealeii/Annotations/orfs_squ.fa"
# )
