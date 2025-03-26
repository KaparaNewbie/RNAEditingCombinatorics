"""
This script is an updated version of create_transcript_pileup.py.
It should perform pileup on the complete ORF of each gene, and then create 4 matrices:
1. Reads (0/1/NaN per base)
2. Unique reads ("Transcripts")
3. Proteins (derived from Reads by translation)
4. Unique proteins
"""

import argparse
from multiprocessing import Pool
from pathlib import Path
from typing import Union

import pandas as pd
import numpy as np
from pybedtools import BedTool

from General.consts import final_words
from General.argparse_utils import abs_path_from_str, expanded_path_from_str
from EditingUtils.summary_utils import execute_notebook
from Alignment.alignment_utils import filter_bam_by_read_quality
from Pileup.mpileup import mpileup
from Pileup.positions import pileup_to_positions
from Pileup.reads import reads_and_unique_reads
from Pileup.proteins import proteins_and_unique_proteins


DataFrameOrSeries = Union[pd.DataFrame, pd.Series]


# def summarize(
#     *,
#     out_dir: Path,
#     template_transcripts_notebook: Path,
#     executed_transcripts_notebook_base_name: str,
#     reads_files: list[Path],
#     transcripts_files: list[Path],
#     condition_col: str,
#     conditions: list[str],
# ):

#     execute_notebook(
#         template_transcripts_notebook,
#         out_dir,
#         executed_transcripts_notebook_base_name,
#         reads_files=[str(reads_file) for reads_file in reads_files],
#         transcripts_files=[
#             str(transcripts_file) for transcripts_file in transcripts_files
#         ],
#         condition_col=condition_col,
#         conditions=conditions,
#     )


def main(
    *,
    transcriptome: Path,
    data_table: Path,
    data_table_sep: str,
    sample_col: str,
    group_col: str,
    region_col: str,
    start_col: str,
    end_col: str,
    strand_col: str,
    path_col: str,
    prob_regions_bed_col: str,
    known_editing_sites: Path,
    cds_regions: Path,
    min_percent_of_max_coverage: float,
    snp_noise_level: float,
    top_x_noisy_positions: int,
    assurance_factor: float,
    include_flags: str,
    exclude_flags: str,
    parity: str,
    out_dir: Path,
    samtools_path: Path,
    processes: int,
    threads: int,
    # transcripts_notebook_template: Path,
    min_rq: Union[float, None],
    min_bq: int,
    out_files_sep: str,
    keep_pileup_files: bool,
    gz_compression: bool,
):
    """
    1. get input

        1. data_table.csv which contains
            * sample name
            * group (here - tissue)
            * region, written in a samtools-acceptable notation (`chr`:`start`-`end`)
            (only look for known editing sites within that region in that sample)
            * path to bam file
            * prob regions bed [optional] - exclude these regions from the analysis, if they are present in the region column
        it must have a header with the supplied cols' names

        2. bed file with known editing sites

    (# todo have `prepare_data.py` and `align.py` accept such a table too, and change `prepare_data.py` s.t. the
    update table it'll produce will use `align.py`; in the same manner, the updated table produced by `align.py` will
    be used for this script)

    2. what to do with reads that are too short?
        1 - discard completly
        2 - define another sub minimal span within the designated region
        3 - just introduce a "partial" annotation for those reads

    3. run mpileup for each set of reads
    (following the previous steps, there may be a number of subsets of reads for each sample)

    4. parse the pileup files to sets of transcripts according the editing sites they consist of
    """

    out_dir.mkdir(exist_ok=True)

    # 1 - get data

    data_table = pd.read_csv(data_table, sep=data_table_sep)
    samples = data_table[
        sample_col
    ].tolist()  # needed if there are number of samples in each group
    groups = [str(group) for group in data_table[group_col]]
    regions = data_table[region_col].tolist()
    starts = data_table[start_col].tolist()
    ends = data_table[end_col].tolist()
    strands = data_table[strand_col].tolist()
    bam_files = [Path(path).absolute() for path in data_table[path_col]]
    prob_regions_beds = list(
        map(
            lambda path: BedTool(path) if path else None,
            data_table[prob_regions_bed_col].replace({np.nan: ""}),
        )
    )
    # prob_regions_beds = [BedTool(path) for path in data_table[prob_regions_bed_col]]

    # 2 - filter reads by mapping quality

    if min_rq is not None:
        with Pool(processes=processes) as pool:
            filtered_bam_files = pool.starmap(
                func=filter_bam_by_read_quality,
                iterable=[
                    (samtools_path, in_bam, min_rq, threads, out_dir)
                    for in_bam in bam_files
                ],
            )
    else:
        filtered_bam_files = bam_files

    # 3 - run mpileup

    # pileup_files = [Path(out_dir, f"{sample}.pileup") for sample in samples]
    pileup_files = [
        Path(out_dir, f"{filtered_bam_file.stem}.pileup")
        for filtered_bam_file in filtered_bam_files
    ]
    pileup_formatted_regions = [
        f"{region}:{start+1}-{end}" for region, start, end in zip(regions, starts, ends)
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
                )
                for region, in_bam, out_pileup in zip(
                    pileup_formatted_regions,
                    filtered_bam_files,
                    pileup_files,
                )
            ],
        )

    compression_postfix = ".gz" if gz_compression else ""

    # 4 - pileup files -> positions dfs

    positions_files = [
        Path(out_dir, f"{pileup_file.stem}.positions.csv{compression_postfix}")
        for pileup_file in pileup_files
    ]
    reads_mapping_files = [
        Path(out_dir, f"{pileup_file.stem}.OldToNewReads.csv{compression_postfix}")
        for pileup_file in pileup_files
    ]

    with Pool(processes=int(processes / 2)) as pool:
        pool.starmap(
            func=pileup_to_positions,
            iterable=[
                (
                    pileup_file,
                    strand,
                    min_percent_of_max_coverage,
                    snp_noise_level,
                    top_x_noisy_positions,
                    assurance_factor,
                    prob_regions_bed,
                    known_editing_sites,
                    cds_regions,
                    positions_file,
                    reads_mapping_file,
                    out_files_sep,
                    keep_pileup_files,
                )
                for pileup_file, reads_mapping_file, strand, prob_regions_bed, positions_file in zip(
                    pileup_files,
                    reads_mapping_files,
                    strands,
                    prob_regions_beds,
                    positions_files,
                )
            ],
        )

    # 5 - positions dfs -> reads & unique reads dfs

    reads_files = [
        Path(out_dir, f"{pileup_file.stem}.reads.csv{compression_postfix}")
        for pileup_file in pileup_files
    ]
    unique_reads_files = [
        Path(out_dir, f"{pileup_file.stem}.unique_reads.csv{compression_postfix}")
        for pileup_file in pileup_files
    ]

    with Pool(processes=processes) as pool:
        pool.starmap(
            func=reads_and_unique_reads,
            iterable=[
                (
                    positions_file,
                    strand,
                    group_col,
                    group,
                    parity,
                    reads_file,
                    unique_reads_file,
                    out_files_sep,
                )
                for positions_file, strand, group, reads_file, unique_reads_file in zip(
                    positions_files, strands, groups, reads_files, unique_reads_files
                )
            ],
        )

    # 6 - reads & unique reads dfs -> proteins & unique proteins dfs

    proteins_files = [
        Path(out_dir, f"{pileup_file.stem}.proteins.csv{compression_postfix}")
        for pileup_file in pileup_files
    ]
    unique_proteins_files = [
        Path(out_dir, f"{pileup_file.stem}.unique_proteins.csv{compression_postfix}")
        for pileup_file in pileup_files
    ]

    with Pool(processes=processes) as pool:
        pool.starmap(
            func=proteins_and_unique_proteins,
            iterable=[
                (
                    unique_reads_file,
                    region,
                    start,
                    end,
                    strand,
                    transcriptome,
                    group_col,
                    proteins_file,
                    unique_proteins_file,
                    out_files_sep,
                )
                for unique_reads_file, region, start, end, strand, proteins_file, unique_proteins_file in zip(
                    unique_reads_files,
                    regions,
                    starts,
                    ends,
                    strands,
                    proteins_files,
                    unique_proteins_files,
                )
            ],
        )


def define_args() -> argparse.Namespace:

    # create parser
    # description = "TBD description"
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        # description=description
    )
    # parser = argparse.ArgumentParser()

    # define args
    parser.add_argument(
        "--transcriptome",
        required=True,
        type=abs_path_from_str,
        help=("A transcritpome reference fasta file. Should have no introns."),
    )
    parser.add_argument("--data_table", type=abs_path_from_str, required=True)
    parser.add_argument(
        "--data_table_sep", default=",", help="Delimiter used in `data_table`."
    )
    parser.add_argument(
        "--sample_col", default="Sample", help="Sample col label in `data_table`."
    )
    parser.add_argument(
        "--group_col", default="Gene", help="Group col label in `data_table`."
    )
    parser.add_argument(
        "--region_col",
        default="Region",
        help=(
            "Region col label in `data_table`."
            "It's values should be given in a samtools-compatible format, i.e., either `chrom` of `chrom:start-end`, "
            "where `start` is 1-based and `end` - inclusive. "
            "Should represent an ORF (or part of it) in the transcriptome, so its length should be a multiple of 3."
        ),
    )
    parser.add_argument(
        "--start_col",
        default="Start",
        help=(
            "Start col label in `data_table`. Its values should be 0-based (as in BED format).",
        ),
    )
    parser.add_argument(
        "--end_col",
        default="End",
        help=(
            "End col label in `data_table`. Its values should be inclusive (as in BED format)."
        ),
    )
    parser.add_argument(
        "--strand_col", default="Strand", help="Strand col label in `data_table`."
    )
    parser.add_argument(
        "--path_col", default="Path", help="Path col label in `data_table`."
    )
    parser.add_argument(
        "--prob_regions_bed_col",
        default="ProbRegionsBED",
        help="Problematic regions to exclude col label in `data_table`. Its values shpuld be paths to BED files.",
    )
    parser.add_argument(
        "--known_editing_sites",
        type=abs_path_from_str,
        required=True,
        help=(
            "6-col bed file with known A-to-I editing sites. "
            "Currently, only editing sites on the positive strand are supported."
        ),
    )
    parser.add_argument(
        "--cds_regions",
        type=abs_path_from_str,
        help=(
            "6-col bed file of coding regions. "
            "Currently, only positions insided those regions are used for isoforms analysis. "
            "Pay attention! In its abscence, all positions are assumed to be within coding regions!"
        ),
    )
    parser.add_argument(
        "--min_percent_of_max_coverage",
        type=float,
        default=0.1,
        help="Keep only positions with coverage >= min_percent_of_max_coverage * max_coverage.",
    )
    parser.add_argument(
        "--snp_noise_level",
        type=float,
        default=0.1,
        help="Treat non-refbase positions with noise >= snp_noise_level as SNPs, so they won't be considered for noise threshold.",
    )
    parser.add_argument(
        "--top_x_noisy_positions",
        type=int,
        default=1,
        help="Consider the top X noisy positions to define initial noise threshold.",
    )
    parser.add_argument(
        "--assurance_factor",
        type=float,
        default=1.5,
        help=(
            "The noise in a certain position is defined as number of bases of the must abundant alt_base, "
            "divided by the same number + the number of ref_base bases. "
            "This is done for all ref_bases that aren't A on the positive strand or T on the negative strand."
            "Then, the actual noise level is determined by the "
            "initial noise threshold (see `top_x_noisy_positions`) * assurance_factor. "
            "That actual noise level is used to decide wether a position undergoes RNA editing."
        ),
    )
    parser.add_argument(
        "--include_flags",
        help="Require reads with this flag. 3 is rcommended for PE reads.",
    )
    parser.add_argument(
        "--exclude_flags",
        help=(
            "Exclude reads with this flag. "
            "The recommended 2304 remove secondary and supplementary (chimeric) alignments.",
        ),
    )
    parser.add_argument(
        "--parity",
        help="Parity of the reads.",
        choices=["SE", "se", "PE", "pe"],
        required=True,
    )
    parser.add_argument("--out_dir", type=abs_path_from_str, required=True)
    parser.add_argument(
        "--samtools_path",
        type=expanded_path_from_str,
        default=Path("~/anaconda3/envs/combinatorics/bin/samtools").expanduser(),
        help="Samtools executable.",
    )
    parser.add_argument(
        "--processes",
        type=int,
        default=8,
        help="Maximal number of processes to run in parallel. Some memory consuming may use fewer.",
    )
    parser.add_argument(
        "--threads", type=int, default=10, help="Threads used in each process."
    )
    # parser.add_argument(
    #     "--transcripts_notebook_template",
    #     default=Path("Code/Notebooks/transcripts_w_cond_w_nan.ipynb").absolute(),
    #     type=abs_path_from_str,
    #     help="Summary notebook template file.",
    # )
    parser.add_argument(
        "--min_rq",
        # default=0.99,
        type=float,
        help="Minimum predicted accuracy in [0, 1]. Use this flag for PacBio reads.",
    )
    parser.add_argument(
        "--min_bq",
        default=13,
        type=int,
        help=(
            "Minimum base quality for a base to be considered. "
            "Note base-quality 0 is used as a filtering mechanism for overlap removal. "
            "Hence using --min-BQ 0 will disable the overlap removal code and act "
            "as if the --ignore-overlaps option has been set. "
            "Use this flag for Illumina reads."
        ),
    )
    parser.add_argument(
        "--out_files_sep", default="\t", help="Delimiter for csv out files."
    )
    parser.add_argument(
        "--keep_pileup_files",
        action="store_true",
        help="Keep the pileup files from which the $sample.positions.csv files are parsed.",
    )
    parser.add_argument("--gz_compression", action="store_true")

    return parser

    # # parse args
    # args = parser.parse_args()
    # return args


if __name__ == "__main__":

    # run

    # main(
    #     **vars(define_and_parse_args())
    # )  # https://stackoverflow.com/a/35824590/10249633 argparse.Namespace -> Dict

    parser = define_args()
    args = parser.parse_args()
    main(
        **vars(args)
    )  # https://stackoverflow.com/a/35824590/10249633 argparse.Namespace -> Dict

    # end
    final_words()
