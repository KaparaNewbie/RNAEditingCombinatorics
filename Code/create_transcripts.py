import argparse
import subprocess
from multiprocessing import Pool
from pathlib import Path
from random import choice, seed
from typing import Union

import pandas as pd
from numpy import NaN
from icecream import ic
from frozendict import frozendict
from more_itertools import powerset
from pybedtools import BedTool

from General.consts import final_words
from General.argparse_utils import abs_path_from_str, expanded_path_from_str
from EditingUtils.summary_utils import execute_notebook
from Alignment.alignment_utils import filter_bam_by_read_quality

DataFrameOrSeries = Union[pd.DataFrame, pd.Series]


def count_reads(
    samtools_path: Path, in_bam: Path, region: str, exclude_flags: str, threads: int
):
    cmd = (
        f"{samtools_path} "
        "view "
        f"-F {exclude_flags} "
        "-c "
        f"--threads {threads} "
        f"{in_bam} "
        f"{region} "
    )
    reads = int(subprocess.run(cmd, shell=True, capture_output=True).stdout.decode())
    return reads


def define_max_depth(
    samtools_path: Path,
    in_bam: Path,
    region: str,
    exclude_flags: str,
    threads: int,
    max_depth_padding: int = 10_000,
):
    legal_reads_in_region = count_reads(
        samtools_path, in_bam, region, exclude_flags, threads
    )
    return legal_reads_in_region + max_depth_padding


def mpileup(
    samtools_path: Path,
    genome: Path,
    known_editing_sites: Path,
    region: str,
    prob_regions_bed: Union[BedTool, None],
    exclude_flags: str,
    in_bam: Path,
    out_file: Path,
    threads: int,
):
    """
    samtools mpileup \
    --fasta-ref D.pealeii/Annotations/orfs_squ.fa \
    --positions D.pealeii/Annotations/D.pea.EditingSites.bed \
    --region comp141882_c0_seq14 \
    --no-BAQ \
    --no-output-ins --no-output-ins \
    --no-output-del --no-output-del \
    --no-output-ends \
    --output-QNAME \
    D.pealeii/Alignment/Naive/PCLO-CNS-RESUB.C0x1291.aligned.sorted.bam \
    --output D.pealeii/Alignment.Test/PCLO.mpileup
    """
    # filter known editing sites by problamtic regions
    known_editing_sites = BedTool(known_editing_sites)
    filtered_known_editing_sites_file = Path(
        out_file.parent, f"{out_file.stem}.FilteredKnownSites.bed"
    )
    if prob_regions_bed is not None:
        known_editing_sites = known_editing_sites.intersect(prob_regions_bed, v=True)
    known_editing_sites.saveas(filtered_known_editing_sites_file)

    max_depth = define_max_depth(samtools_path, in_bam, region, exclude_flags, threads)
    mpileup_cmd = (
        f"{samtools_path} "
        "mpileup "
        f"--fasta-ref {genome} "
        f"--positions {filtered_known_editing_sites_file} "
        f"--region {region} "
        # "--no-BAQ "  # disable base alignment quality (BAQ) computation
        "--no-output-ins --no-output-ins "
        "--no-output-del --no-output-del "
        "--no-output-ends "
        "--output-QNAME "
        f"--excl-flags {exclude_flags} "  # 2304 (default) = remove secondary and supplementary (chimeric) alignments
        f"--max-depth {max_depth} "
        f"{in_bam} "
        f"--output {out_file} "
    )
    subprocess.run(mpileup_cmd, shell=True)


def read_pileup(pileup_file: Path, remove_non_edited_positions: bool) -> pd.DataFrame:
    """Read pileup file into a DataFrame.
    Verify format, change to 0-based coordinates, and possibly remove non edited positions.

    Args:
        pileup_file (Path): Path to a file creates by the `mpileup` function.
        remove_non_edited_positions (bool, optional): Remove positions of known editing sites if they weren't edited
        in any of the reads. Defaults to True.

    Raises:
        ValueError: If the 5th col doens't contain only base match/ bash mismatch characters, and/or the number of
        reads' names in the 7th col isn't equal to the number of characters in the 5th col.

    Returns:
        pd.DataFrame: The parsed pileup file.
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
    df = pd.read_csv(pileup_file, sep="\t", names=cols)
    # verify that the 5th col contains only base match/mismatch, and correspoding number
    # of reads' names in the 7th col
    for row_num, row in enumerate(df.itertuples()):
        bases = row.MappedBases
        reads = row.Reads.split(",")
        if len(bases) != len(reads):
            raise ValueError(
                f"Line {row_num} in {pileup_file} contains indels and/or doesn't have the "
                "same number of mapped bases and reads."
            )
    # change editing position to 0-based
    df["Position"] = df["Position"] - 1
    # keep only positions where at least one of the reads had editing
    if remove_non_edited_positions:
        editing_in_position = []
        for row in df.itertuples():
            bases = row.MappedBases
            if "G" in bases or "g" in bases:
                editing_in_position.append(True)
            else:
                editing_in_position.append(False)
        df = df.loc[editing_in_position]
    # return
    return df


def pileup_to_transcripts(
    pileup_file: Path,
    processes: int,
    remove_non_edited_positions: bool,
    unmapped_bases_as_nan: bool,
    group_col: str,
    group: str,
):

    # data per position

    positions_df = read_pileup(pileup_file, remove_non_edited_positions)

    # data per read -> data per read

    unique_reads = {
        read for reads in positions_df["Reads"] for read in reads.split(",")
    }

    edited_positions_in_reads = {read: [] for read in unique_reads}

    for x, row in enumerate(positions_df.itertuples()):
        bases = row.MappedBases
        reads = row.Reads.split(",")
        # deal with reads that were mapped to this pos
        for base, read in zip(bases, reads):

            # edited_in_pos = True if base in ["G", "g"] else False  # todo check
            edited_in_pos = "1" if base in ["G", "g"] else "0"

            edited_positions_in_reads[read].append(edited_in_pos)
        # deal with reads that weren't mapped to this pos
        unmapped_reads_in_pos = unique_reads - set(reads)
        for read in unmapped_reads_in_pos:

            # edited_in_pos = False
            edited_in_pos = "-1" if unmapped_bases_as_nan else "0"

            edited_positions_in_reads[read].append(edited_in_pos)
        # check all reads got values for pos
        mapped_pos_dist = {len(bases) for bases in edited_positions_in_reads.values()}
        if len(mapped_pos_dist) != 1:
            raise Exception(
                f"Problem at line {x} in positions_df: not all reads are mapped."
            )

    reads_df = pd.DataFrame(edited_positions_in_reads)
    reads_df = reads_df.T.reset_index().rename({"index": "Read"}, axis="columns")

    reads_df = reads_df.rename(
        columns={
            old_col: pos
            for old_col, pos in zip(reads_df.columns[1:], positions_df["Position"])
        }
    )

    reads_df.insert(0, group_col, group)

    # data per read -> data per transcript

    transcripts_df = reads_df.copy()

    transcripts_df.insert(
        1,
        "Reads",
        (
            transcripts_df.groupby(transcripts_df.columns[2:].tolist())[
                "Read"
            ].transform(lambda x: ",".join(x))
        ),
    )

    del transcripts_df["Read"]

    transcripts_df.insert(
        2, "NumOfReads", transcripts_df["Reads"].apply(lambda x: len(x.split(",")))
    )

    transcripts_df = (
        transcripts_df.loc[~transcripts_df.iloc[:, 3:].duplicated()]
        .sort_values("NumOfReads", ascending=False)
        .reset_index(drop=True)
    )

    transcripts_df.insert(
        1,
        "Transcript",
        pd.Series(
            transcripts_df[group_col] + "-" + transcripts_df.index.values.astype(str)
        ),
    )

    # transcripts_df = (
    #     reads_df.groupby(reads_df.columns[1:].tolist(), dropna=False)
    #     .size()
    #     .to_frame("TempNumOfReads")
    #     .reset_index()
    # )
    # transcripts_df.insert(0, "NumOfReads", transcripts_df["TempNumOfReads"])
    # del transcripts_df["TempNumOfReads"]
    # transcripts_df = transcripts_df.sort_values(
    #     "NumOfReads", ascending=False
    # ).reset_index(drop=True)

    # with Pool(processes=processes) as pool:
    #     all_supporting_reads = pool.starmap(
    #         func=get_supporting_reads,
    #         iterable=[
    #             (transcripts_df.iloc[transcripts_df_line_num, 1:], reads_df)
    #             for transcripts_df_line_num in range(len(transcripts_df))
    #         ],
    #     )

    # transcripts_df.insert(1, "Reads", all_supporting_reads)

    # end

    return reads_df, transcripts_df


# def get_supporting_reads(transcripts_df_line_wo_count_col, reads_df):
#     reads_df_wo_read_col = reads_df.iloc[:, 1:]
#     lines_with_supporting_reads = reads_df.loc[
#         transcripts_df_line_wo_count_col.eq(reads_df_wo_read_col).all(axis="columns")
#     ]
#     return ",".join(read for read in lines_with_supporting_reads["Read"])


def no_na_cells(df: DataFrameOrSeries, first_col=0) -> bool:
    """Return `True` if there are no `NaN` cells."""
    df = pd.DataFrame(df)  # in case df is actually a series
    df = df.iloc[:, first_col:]
    any_na_cells = (
        df.isna().any().any()
    )  # 2nd any is needed for df and doesn't harm a series
    not_even_one_na_cell = not any_na_cells
    return not_even_one_na_cell


def no_valid_cells(df: DataFrameOrSeries, first_col=0) -> bool:
    """Return `True` if all cells are `NaN`."""
    df = pd.DataFrame(df)  # in case df is actually a series
    df = df.iloc[:, first_col:]
    only_na_cells = df.isna().all().all()
    return only_na_cells


def rows_opt_3(
    df: pd.DataFrame, temp_seed: float = 1892.0, first_col: int = 0, debug: bool = False
):
    """Find a sub-`df` with maximal subgroup of rows and cols that are all valid (i.e., different than `NaN`)."""

    df = pd.DataFrame(df)

    if no_na_cells(df, first_col):  # no optimization is needed :)
        return df

    if no_valid_cells(df, first_col):
        raise Exception("All cells are NaN.")

    seed(temp_seed)

    row_names = df.index.values
    col_names = df.columns[first_col:]

    # k: row_name, v: col_names to delete in order to satisfy row_name
    delete_cols_to_satisfy_row = {}
    for row_name in row_names:  # todo parallelize this loop over rows
        row = df.iloc[row_name, first_col:]
        delete_cols_to_satisfy_row[row_name] = frozenset(
            col_name
            for col_name, col_is_na_in_row in zip(row.index, row.isna())
            if col_is_na_in_row
        )
    if debug:
        ic(delete_cols_to_satisfy_row)

    # k: col_name, v: row_names satisfied if col_name isn't deleted
    possible_rows_satisfied_by_col_retention = {
        col_name: frozenset(
            row_name
            for row_name in delete_cols_to_satisfy_row
            if col_name not in delete_cols_to_satisfy_row[row_name]
        )
        for col_name in col_names
    }
    if debug:
        ic(possible_rows_satisfied_by_col_retention)

    disjoint_subsets_of_possible_satisfied_rows = set()
    for (
        col_name_1,
        satisfied_rows_1,
    ) in possible_rows_satisfied_by_col_retention.items():
        # the elements of disjoint_subsets_of_possible_satisfied_rows are, of course, disjoint...
        satisfied_rows_1_are_contained = False
        for (
            col_name_2,
            satisfied_rows_2,
        ) in possible_rows_satisfied_by_col_retention.items():
            if col_name_1 != col_name_2 and satisfied_rows_1 < satisfied_rows_2:
                satisfied_rows_1_are_contained = True
                break
        if satisfied_rows_1_are_contained:
            continue
        # find all cols whose satisfied rows are a subset of satisfied_rows_1
        possible_rows_satisfied_by_col_retention_subset = {}
        for (
            col_name_2,
            satisfied_rows_2,
        ) in possible_rows_satisfied_by_col_retention.items():
            if satisfied_rows_2 <= satisfied_rows_1:
                possible_rows_satisfied_by_col_retention_subset[
                    col_name_2
                ] = satisfied_rows_2
        disjoint_subsets_of_possible_satisfied_rows.add(
            frozendict(possible_rows_satisfied_by_col_retention_subset)
        )
    if debug:
        ic(disjoint_subsets_of_possible_satisfied_rows)

    # todo: return a separate df for each disjoint subset of possible satisfied rows?

    max_valid_cells = 0
    final_rows = None
    final_cols = None

    for (
        possible_rows_satisfied_by_col_retention_subset
    ) in disjoint_subsets_of_possible_satisfied_rows:

        # find the main_col and its corresponding main_rows - the rows that could possibly by satisfied by that col retention
        # (recall that all other rows in possible_rows_satisfied_by_col_retention_subset are subsets of main_rows)

        subsets_of_satisfied_rows = set(
            possible_rows_satisfied_by_col_retention_subset.values()
        )  # there may be two or more cols satisfying the same rows
        max_rows_in_subset = max(len(rows) for rows in subsets_of_satisfied_rows)
        main_col = None
        for col_name in possible_rows_satisfied_by_col_retention_subset:
            # there may be two or more cols with max rows (i.e., with the same rows) but that's ok
            if (
                len(possible_rows_satisfied_by_col_retention_subset[col_name])
                == max_rows_in_subset
            ):
                main_col = col_name
                break
        if main_col is None:
            raise Exception("main_col is None")  # shouldn't really happen
        else:
            if debug:
                ic(main_col)
        main_rows = possible_rows_satisfied_by_col_retention_subset[main_col]

        # iterate over all possible subsets of main_rows and find other rows that can agree with them
        # this should work for all the following cases:
        # 1 - all cols satisfy the same rows
        # (which is also the case if possible_rows_satisfied_by_col_retention_subset contains only one col)
        # 2 - there are two or more cols satisfying different rows, and one of these rows is a superset of all the others

        main_rows_powerset = [
            set(main_rows_subset)
            for main_rows_subset in powerset(main_rows)
            if 0 < len(main_rows_subset)  # don't include empty sets
        ]
        if debug:
            ic(main_rows)
            ic(len(main_rows_powerset))
        for main_rows_subset in main_rows_powerset:
            cols_matching_main_rows_subset = tuple(
                col_name
                for col_name, satisfied_rows in possible_rows_satisfied_by_col_retention_subset.items()
                if main_rows_subset <= satisfied_rows
            )
            valid_cells = len(main_rows_subset) * len(cols_matching_main_rows_subset)
            if (valid_cells > max_valid_cells) or (
                valid_cells == max_valid_cells and choice([0, 1]) == 0
            ):
                max_valid_cells = valid_cells
                final_rows = main_rows_subset
                final_cols = cols_matching_main_rows_subset

    final_rows = list(final_rows)
    final_cols = list(final_cols)
    df = df.loc[final_rows, final_cols]
    return df


def repeat_rand_opt():
    # base_seed = 1892
    # seed(base_seed)
    # first_state = getstate()
    # num_of_temp_seeds = 10
    # temp_seeds = [random() for _ in range(num_of_temp_seeds)]
    # for temp_seed in temp_seeds:
    #     change_seed(temp_seed)
    # seed(base_seed)
    # new_temp_seeds = [random() for _ in range(num_of_temp_seeds)]
    # print(f"{temp_seeds = }")
    # print(f"{new_temp_seeds = }")
    # assert temp_seeds == new_temp_seeds
    pass


def summarize(
    *,
    out_dir: Path,
    template_transcripts_notebook: Path,
    executed_transcripts_notebook_base_name: str,
    reads_files: list[Path],
    transcripts_files: list[Path],
    condition_col: str,
    conditions: list[str],
):

    execute_notebook(
        template_transcripts_notebook,
        out_dir,
        executed_transcripts_notebook_base_name,
        reads_files=[str(reads_file) for reads_file in reads_files],
        transcripts_files=[
            str(transcripts_file) for transcripts_file in transcripts_files
        ],
        condition_col=condition_col,
        conditions=conditions,
    )


def main(
    *,
    genome: Path,
    data_table: Path,
    data_table_sep: str,
    sample_col: str,
    group_col: str,
    region_col: str,
    path_col: str,
    prob_regions_bed_col: str,
    known_editing_sites: Path,
    keep_non_edited_positions: bool,
    unmapped_bases_as_false: bool,
    exclude_flags: str,
    out_dir: Path,
    samtools_path: Path,
    processes: int,
    threads: int,
    transcripts_notebook_template: Path,
    min_rq: float,
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

        2. bed file with known editing sites  (# todo deal with editing sites on the negative strand)

    (# todo have `prepare_data.py` and `align.py` accept such a table too, and change `prepare_data.py` s.t. the
    update table it'll produce will use `align.py`; in the same manner, the updated table produced by `align.py` will
    be used for this script)

    2. what to do with reads that are too short?
        1 - discard completly
        2 - define another sub minimal span within the designated region
        3 - just introduce a "partial" annotation for those reads

    3. run mpileup for each set of reads
    (following the previous steps, there may be a number of subsets of reads for each sample)
    (# todo split mpileup to subregions in case of large data)

    4. parse the pileup files to sets of transcripts according the editing sites they consist of
    """

    out_dir.mkdir(exist_ok=True)

    # 1 - get data

    data_table = pd.read_csv(data_table, sep=data_table_sep)
    samples = data_table[sample_col].tolist()
    groups = [str(group) for group in data_table[group_col]]
    regions = data_table[region_col].tolist()
    bam_files = [Path(path).absolute() for path in data_table[path_col]]
    prob_regions_beds = list(
        map(
            lambda path: BedTool(path) if path else None,
            data_table[prob_regions_bed_col].replace({NaN: ""}),
        )
    )
    # prob_regions_beds = [BedTool(path) for path in data_table[prob_regions_bed_col]]

    # 2 - filter reads by mapping quality
    # todo: any other filtering? (e.g. by read length)
    with Pool(processes=processes) as pool:
        filtered_bam_files = pool.starmap(
            func=filter_bam_by_read_quality,
            iterable=[
                (samtools_path, in_bam, min_rq, threads, out_dir)
                for in_bam in bam_files
            ],
        )

    # 3 - run mpileup

    # pileup_files = [Path(out_dir, f"{sample}.pileup") for sample in samples]
    pileup_files = [
        Path(out_dir, f"{filtered_bam_file.stem}.pileup")
        for filtered_bam_file in filtered_bam_files
    ]
    with Pool(processes=processes) as pool:
        pool.starmap(
            func=mpileup,
            iterable=[
                (
                    samtools_path,
                    genome,
                    known_editing_sites,
                    region,
                    prob_regions_bed,
                    exclude_flags,
                    in_bam,
                    out_pileup,
                    threads,
                )
                for region, prob_regions_bed, in_bam, out_pileup in zip(
                    regions, prob_regions_beds, filtered_bam_files, pileup_files
                )
            ],
        )

    # 4 - pileup files -> sets of transcripts

    remove_non_edited_positions = not keep_non_edited_positions
    reads_files = []
    transcripts_files = []
    unmapped_bases_as_nan = not unmapped_bases_as_false
    for sample, group, pileup_file in zip(samples, groups, pileup_files):

        reads_df, transcripts_df = pileup_to_transcripts(
            pileup_file,
            processes,
            remove_non_edited_positions,
            unmapped_bases_as_nan,
            group_col,
            group,
        )

        # reads_file = Path(out_dir, f"{sample}.reads.csv")
        # transcripts_file = Path(out_dir, f"{sample}.transcripts.csv")
        reads_file = Path(out_dir, f"{pileup_file.stem}.reads.csv")
        transcripts_file = Path(out_dir, f"{pileup_file.stem}.transcripts.csv")

        reads_files.append(reads_file)
        transcripts_files.append(transcripts_file)

        reads_df.to_csv(reads_file, index=False, na_rep=NaN)
        transcripts_df.to_csv(transcripts_file, index=False, na_rep=NaN)

    summarize(
        out_dir=out_dir,
        template_transcripts_notebook=transcripts_notebook_template,
        executed_transcripts_notebook_base_name=f"ReadsAndTranscripts",
        reads_files=reads_files,
        transcripts_files=transcripts_files,
        condition_col=group_col,
        conditions=groups,
    )

    # 4.1 - sampled fractions of reads: pileup files -> sets of transcripts


def define_and_parse_args() -> argparse.Namespace:

    # create parser
    # description = "TBD description"
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        # description=description
    )
    # parser = argparse.ArgumentParser()

    # define args
    parser.add_argument(
        "--genome",
        required=True,
        type=abs_path_from_str,
        help=(
            "Genome/transcritpome (mostly transcrtiptome as we don't deal with strand) "
            "reference fasta file."
        ),
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
        "--region_col", default="Region", help="Region col label in `data_table`."
    )
    parser.add_argument(
        "--path_col", default="Path", help="Path col label in `data_table`."
    )
    parser.add_argument(
        "--prob_regions_bed_col",
        default="ProbRegionsBED",
        help="Problematic regions to exclude col label in `data_table`.",
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
        "--keep_non_edited_positions",
        action="store_true",
        help=(
            "Keep known positions (defined by `known_editing_sites`) even if they "
            "weren't edited in any of the reads in a sample."  # todo change the condition from sample to group
        ),
    )
    # parser.add_argument(
    #     "--unmapped_bases_as_nan",
    #     action="store_true",
    #     help=(
    #         "By default, in each read, edited positions are marked as `True`, while both unedited and "
    #         "unmapped positions are marked as `False`. Setting this flag will mark unmapped positions "
    #         "with `NaN` instead of `False`."
    #     ),
    # )
    parser.add_argument(
        "--unmapped_bases_as_false",
        action="store_true",
        help=(
            "By default, in each read, edited positions are marked as `1` (`True`), unedited ones as "
            "`0` (`False`), and unmapped positions as `-1` (`NaN`). Setting this flag will mark unmapped "
            "positions with `0` (False) instead of `-1` (`NaN`)."
        ),
    )
    parser.add_argument(
        "--exclude_flags",
        default="2304",
        help="Exclude reads with this flag. The default 2304 remove secondary and supplementary (chimeric) alignments.",
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
        default=6,
        help="Maximal number of processes to run in parallel.",
    )
    parser.add_argument(
        "--threads", type=int, default=15, help="Threads used in each process."
    )
    parser.add_argument(
        "--transcripts_notebook_template",
        default=Path("Code/Notebooks/transcripts_w_cond_w_nan.ipynb").absolute(),
        type=abs_path_from_str,
        help="Summary notebook template file.",
    )
    parser.add_argument(
        "--min_rq",
        default=0.99,
        type=float,
        help="Minimum predicted accuracy in [0, 1].",
    )

    # parse args
    args = parser.parse_args()
    return args


if __name__ == "__main__":

    # run
    main(
        **vars(define_and_parse_args())
    )  # https://stackoverflow.com/a/35824590/10249633 argparse.Namespace -> Dict

    # end
    final_words()
