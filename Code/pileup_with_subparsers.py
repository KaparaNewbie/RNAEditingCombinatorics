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
from collections import defaultdict
import subprocess

import pandas as pd
import numpy as np
from pybedtools import BedTool
from icecream import ic

from General.argparse_utils import abs_path_from_str, expanded_path_from_str
from Alignment.alignment_utils import filter_bam_by_read_quality, sample_bam
from Pileup.mpileup import mpileup
from Pileup.positions import (
    pileup_to_positions,
    # multisample_pileups_to_positions_old,
    multisample_pileups_to_positions_all_transcripts,
)
from Pileup.reads import reads_and_unique_reads, multisample_reads_and_unique_reads
from Pileup.proteins import (
    proteins_and_unique_proteins,
    multisample_proteins_and_unique_proteins,
)
from General.consts import final_words, ic_prefix

# configure icecream to print the time of the print and the context (file, line, function)
ic.configureOutput(includeContext=True, prefix=ic_prefix)


DataFrameOrSeries = Union[pd.DataFrame, pd.Series]


def undirected_sequencing_main(
    *,
    transcriptome: Path,
    samples_table: Path,
    samples_table_sep: str,
    sample_col: str,
    group_col: str,
    known_editing_sites: Path,
    cds_regions: Path,
    min_percent_of_max_coverage: float,
    min_mapped_reads_per_position: Union[int, None],
    include_flags: str,
    exclude_flags: str,
    parity: str,
    snp_noise_level: float,
    top_x_noisy_positions: int,
    assurance_factor: float,
    pooled_transcript_noise_threshold: float,
    out_dir: Path,
    samtools_path: Path,
    processes: int,
    threads: int,
    min_rq: Union[float, None],
    min_bq: int,
    out_files_sep: str,
    keep_pileup_files: bool,
    keep_bam_files: bool,
    override_existing_pileup_files: bool,
    override_existing_bam_files: bool,
    gz_compression: bool,
    alignments_stats_table: Path,
    alignments_stats_table_sep: str,
    min_samples: int,
    min_mapped_reads_per_sample: float,
    total_mapped_reads: int,
    min_known_sites: int,
    main_by_chrom_dir: Path,
    postfix: str,
    interfix_start: str,
    remove_non_refbase_noisy_positions: bool,
    known_sites_only: bool,
    pileup_dir_name: str,  # "PileupFiles"
    positions_dir_name: str,  # "PositionsFiles"
    reads_dir_name: str,  # "ReadsFiles"
    proteins_dir_name: str,  # "ProteinsFiles"
    sample_reads: bool,
    num_sampled_reads: int,
    seed: int,
    final_editing_scheme: str,
    alternative_hypothesis: str,
    # binom_noise_pval_col: str,
    # bh_noise_pval_col: str,
    # bh_noisy_col: str,
    # binom_editing_pval_col: str,
    # bh_editing_pval_col: str,
    # bh_editing_col: str,
    **kwargs,
):
    prob_regions_bed = None  # we don't deal with individual "problamitic" sites in such large scale analysis

    # denovo_detection = True
    denovo_detection = not known_sites_only

    out_dir.mkdir(exist_ok=True)

    # 1 - get data

    cds_df = pd.read_csv(
        cds_regions,
        sep="\t",
        names="Chrom Start End Name Score Strand".split(),
        comment="#",
    )
    # samples_df = pd.read_csv(samples_table, sep=samples_table_sep)
    # sample_to_group_dict = {
    #     row[sample_col]: row[group_col] for _, row in samples_df.iterrows()
    # }
    alignments_stats_df = pd.read_csv(
        alignments_stats_table, sep=alignments_stats_table_sep
    )

    alignments_stats_df = alignments_stats_df.loc[
        (alignments_stats_df["Samples"] >= min_samples)
        & (alignments_stats_df["MappedReadsPerSample"] >= min_mapped_reads_per_sample)
        & (alignments_stats_df["KnownSites"] >= min_known_sites)
        & (alignments_stats_df["MappedReads"] >= total_mapped_reads)
    ]

    alignments_stats_df = alignments_stats_df.merge(cds_df, how="left")

    # todo - remove - this is for debugging
    # test_chroms = [
    #     "comp183313_c0_seq12",
    #     "comp162994_c0_seq1",
    #     "comp183909_c0_seq7",
    #     "comp181233_c0_seq10",
    #     "comp183670_c0_seq3",
    #     "comp183256_c0_seq35",
    #     "comp183782_c0_seq5",
    #     "comp183377_c0_seq11",
    #     "comp181723_c2_seq2",
    #     "comp169467_c0_seq1",
    #     "comp183713_c0_seq9",
    # ]
    # alignments_stats_df = alignments_stats_df.loc[
    #     alignments_stats_df["Chrom"].isin(test_chroms)
    # ]

    chroms = []
    swissprot_names = []
    samples = []
    # groups = []
    orfs_starts = []
    orfs_ends = []
    pileup_formatted_regions = []
    orfs_strands = []
    bam_files = []
    for _, row in alignments_stats_df.iterrows():
        chrom = row["Chrom"]

        # # todo - remove - this is for debugging
        # if chrom not in test_chroms:
        #     continue

        swissprot_name = row["Name"]
        start = row["Start"]
        end = row["End"]
        strand = row["Strand"]
        region = f"{chrom}:{start+1}-{end}"
        bams_in_chrom = list(
            Path(main_by_chrom_dir, chrom).glob(f"*{postfix}")
        )  # should be something like `*.bam`
        ic(bams_in_chrom)  # todo - remove
        for bam_in_chrom in bams_in_chrom:
            # ic(bam_in_chrom)
            sample = bam_in_chrom.stem[: bam_in_chrom.stem.index(interfix_start)]
            # group = sample_to_group_dict[sample]
            chroms.append(chrom)
            swissprot_names.append(swissprot_name)
            samples.append(sample)
            # groups.append(group)
            orfs_starts.append(start)
            orfs_ends.append(end)
            pileup_formatted_regions.append(region)
            orfs_strands.append(strand)
            bam_files.append(bam_in_chrom)

    # 2 - filter reads by mapping quality

    if min_rq is not None:
        filterd_bams_dir = Path(out_dir, "FilterdBAMs")
        filterd_bams_dir.mkdir(exist_ok=True)
        with Pool(processes=processes) as pool:
            filtered_bam_files = pool.starmap(
                func=filter_bam_by_read_quality,
                iterable=[
                    (
                        samtools_path,
                        in_bam,
                        min_rq,
                        threads,
                        filterd_bams_dir,
                        override_existing_bam_files,
                    )
                    for in_bam in bam_files
                ],
            )
    else:
        filtered_bam_files = bam_files

    # 3 - sample reads to desired number of reads

    if sample_reads:
        sampled_bams_dir = Path(out_dir, "SampledBAMs")
        sampled_bams_dir.mkdir(exist_ok=True)
        sampled_filtered_bam_files = [
            Path(sampled_bams_dir, f"{bam_file.stem}.Sampled{num_sampled_reads}.bam")
            for bam_file in filtered_bam_files
        ]
        with Pool(processes=processes) as pool:
            pool.starmap(
                func=sample_bam,
                iterable=[
                    (samtools_path, in_bam, out_bam, num_sampled_reads, seed, threads)
                    for in_bam, out_bam in zip(
                        filtered_bam_files, sampled_filtered_bam_files
                    )
                ],
            )
    else:
        sampled_filtered_bam_files = filtered_bam_files

    # 3 - run mpileup

    # pileup_dir = Path(out_dir, "PileupFiles")
    pileup_dir = Path(out_dir, pileup_dir_name)
    pileup_dir.mkdir(exist_ok=True)
    pileup_files = [
        Path(pileup_dir, f"{sample}.{chrom}.pileup")
        for sample, chrom in zip(samples, chroms)
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
                    override_existing_pileup_files,
                )
                for region, in_bam, out_pileup in zip(
                    pileup_formatted_regions,
                    # filtered_bam_files,
                    sampled_filtered_bam_files,
                    pileup_files,
                )
            ],
        )

    if not keep_bam_files:
        if min_rq is not None:
            subprocess.run(f"rm -rf {filterd_bams_dir}", shell=True)
        # delete sampled reads, whether they were previously filtered by base quality or not
        if sample_reads:
            subprocess.run(f"rm -rf {sampled_bams_dir}", shell=True)

    compression_postfix = ".gz" if gz_compression else ""

    # 4 - pileup files -> positions dfs

    group_col = "Transcript"  # todo make more visible

    unique_chroms = []
    unique_orfs_starts = []
    unique_orfs_ends = []
    unique_orfs_strands = []
    unique_swissprot_names = []
    for chrom, start, end, strand, swissprot_name in zip(
        chroms, orfs_starts, orfs_ends, orfs_strands, swissprot_names
    ):
        if chrom in unique_chroms:
            continue
        unique_chroms.append(chrom)
        unique_orfs_starts.append(start)
        unique_orfs_ends.append(end)
        unique_orfs_strands.append(strand)
        unique_swissprot_names.append(swissprot_name)

    # positions_dir_name = (
    #     "PositionsFiles" if denovo_detection else "PositionsFilesKnownSites"
    # )
    positions_dir = Path(out_dir, positions_dir_name)
    positions_dir.mkdir(exist_ok=True)

    reads_mapping_files = [
        Path(positions_dir, f"{chrom}.OldToNewReads.csv{compression_postfix}")
        for chrom in unique_chroms
    ]
    positions_files = [
        Path(positions_dir, f"{chrom}.positions.csv{compression_postfix}")
        for chrom in unique_chroms
    ]
    corrected_noise_files = [
        Path(positions_dir, f"{chrom}.CorrectedNoise.csv{compression_postfix}")
        for chrom in unique_chroms
    ]
    corrected_editing_files = [
        Path(positions_dir, f"{chrom}.CorrectedEditing.csv{compression_postfix}")
        for chrom in unique_chroms
    ]

    pileup_files_per_chroms = defaultdict(list)
    samples_per_chroms = defaultdict(list)
    for sample, chrom, pileup_file in zip(samples, chroms, pileup_files):
        pileup_files_per_chroms[chrom].append(pileup_file)
        samples_per_chroms[chrom].append(sample)

    binom_noise_pval_col = "NoiseBinomPVal"
    bh_noise_pval_col = "NoiseCorrectedPVal"
    bh_noisy_col = "NoisyCorrected"
    binom_editing_pval_col = "EditingBinomPVal"
    bh_editing_pval_col = "EditingCorrectedPVal"
    bh_editing_col = "EditedCorrected"

    multisample_pileups_to_positions_all_transcripts(
        processes,
        unique_chroms,
        unique_orfs_strands,
        pileup_files_per_chroms,
        samples_per_chroms,
        reads_mapping_files,
        positions_files,
        corrected_noise_files,
        corrected_editing_files,
        min_percent_of_max_coverage,
        min_mapped_reads_per_position,  # min_absolute_coverage
        snp_noise_level,
        top_x_noisy_positions,
        assurance_factor,
        transcriptome,
        final_editing_scheme,
        prob_regions_bed,  # problamatic_regions_file
        known_editing_sites,  # known_sites_file
        cds_regions,  # cds_regions_file
        out_files_sep,
        keep_pileup_files,
        remove_non_refbase_noisy_positions,
        denovo_detection,
        alternative_hypothesis,
        binom_noise_pval_col,
        bh_noise_pval_col,
        bh_noisy_col,
        binom_editing_pval_col,
        bh_editing_pval_col,
        bh_editing_col,
    )

    # the pileup files themsevles were deleted by each processes turning pileup into positions,
    # so pileup_dir is empty now
    if not keep_pileup_files:
        subprocess.run(f"rm -rf {pileup_dir}", shell=True)

    # # todo comment out - this is for debugging
    # return

    # 5 - positions dfs -> reads & unique reads dfs

    # reads_dir_name = "ReadsFiles" if denovo_detection else "ReadsFilesKnownSites"
    reads_dir = Path(out_dir, reads_dir_name)
    reads_dir.mkdir(exist_ok=True)

    reads_files = [
        Path(reads_dir, f"{chrom}.reads.csv{compression_postfix}")
        for chrom in unique_chroms
    ]
    unique_reads_files = [
        Path(reads_dir, f"{chrom}.unique_reads.csv{compression_postfix}")
        for chrom in unique_chroms
    ]

    with Pool(processes=processes) as pool:
        pool.starmap(
            func=multisample_reads_and_unique_reads,
            iterable=[
                (
                    positions_file,
                    strand,
                    group_col,
                    group,
                    parity,
                    snp_noise_level,
                    top_x_noisy_positions,
                    pooled_transcript_noise_threshold,
                    bh_noisy_col,
                    reads_file,
                    unique_reads_file,
                    out_files_sep,
                )
                for positions_file, strand, group, reads_file, unique_reads_file in zip(
                    positions_files,
                    unique_orfs_strands,
                    unique_swissprot_names,
                    reads_files,
                    unique_reads_files,
                )
            ],
        )

    # 6 - reads & unique reads dfs -> proteins & unique proteins dfs

    # proteins_dir_name = (
    #     "ProteinsFiles" if denovo_detection else "ProteinsFilesKnownSites"
    # )
    proteins_dir = Path(out_dir, proteins_dir_name)
    proteins_dir.mkdir(exist_ok=True)

    Path(reads_dir, f"{chrom}.reads.csv{compression_postfix}")

    proteins_files = [
        Path(proteins_dir, f"{chrom}.proteins.csv{compression_postfix}")
        for chrom in unique_chroms
    ]
    unique_proteins_files = [
        Path(proteins_dir, f"{chrom}.unique_proteins.csv{compression_postfix}")
        for chrom in unique_chroms
    ]

    with Pool(processes=processes) as pool:
        pool.starmap(
            func=multisample_proteins_and_unique_proteins,
            iterable=[
                (
                    unique_reads_file,
                    chrom,
                    start,
                    end,
                    strand,
                    transcriptome,
                    group_col,
                    proteins_file,
                    unique_proteins_file,
                    out_files_sep,
                )
                for unique_reads_file, chrom, start, end, strand, proteins_file, unique_proteins_file in zip(
                    unique_reads_files,
                    unique_chroms,
                    unique_orfs_starts,
                    unique_orfs_ends,
                    unique_orfs_strands,
                    proteins_files,
                    unique_proteins_files,
                )
            ],
        )


def directed_sequencing_main(
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
    min_rq: Union[float, None],
    min_bq: int,
    out_files_sep: str,
    keep_pileup_files: bool,
    override_existing_pileup_files: bool,
    gz_compression: bool,
    sample_reads: bool,
    num_sampled_reads: int,
    seed: int,
    **kwargs,
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

    ic()

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
            data_table[prob_regions_bed_col].replace({np.NaN: ""}),
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

    # 3 - sample reads to desired number of reads

    if sample_reads:
        sampled_filtered_bam_files = [
            Path(out_dir, f"{bam_file.stem}.Sampled{num_sampled_reads}.bam")
            for bam_file in filtered_bam_files
        ]
        with Pool(processes=processes) as pool:
            pool.starmap(
                func=sample_bam,
                iterable=[
                    (samtools_path, in_bam, out_bam, num_sampled_reads, seed, threads)
                    for in_bam, out_bam in zip(
                        filtered_bam_files, sampled_filtered_bam_files
                    )
                ],
            )
    else:
        sampled_filtered_bam_files = filtered_bam_files

    # 4 - run mpileup

    # pileup_files = [Path(out_dir, f"{sample}.pileup") for sample in samples]
    pileup_files = [
        Path(out_dir, f"{bam_file.stem}.pileup")
        for bam_file in sampled_filtered_bam_files
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
                    override_existing_pileup_files,
                )
                for region, in_bam, out_pileup in zip(
                    pileup_formatted_regions,
                    sampled_filtered_bam_files,
                    pileup_files,
                )
            ],
        )

    compression_postfix = ".gz" if gz_compression else ""

    # 5 - pileup files -> positions dfs

    positions_files = [
        Path(out_dir, f"{pileup_file.stem}.positions.csv{compression_postfix}")
        for pileup_file in pileup_files
    ]
    reads_mapping_files = [
        Path(out_dir, f"{pileup_file.stem}.OldToNewReads.csv{compression_postfix}")
        for pileup_file in pileup_files
    ]

    remove_non_refbase_noisy_positions = True

    # pileup_to_positions_inputs = [
    #     (
    #         pileup_file,
    #         strand,
    #         min_percent_of_max_coverage,
    #         snp_noise_level,
    #         top_x_noisy_positions,
    #         assurance_factor,
    #         prob_regions_bed,
    #         known_editing_sites,
    #         cds_regions,
    #         positions_file,
    #         reads_mapping_file,
    #         out_files_sep,
    #         keep_pileup_files,
    #         remove_non_refbase_noisy_positions,
    #     )
    #     for pileup_file, reads_mapping_file, strand, prob_regions_bed, positions_file in zip(
    #         pileup_files,
    #         reads_mapping_files,
    #         strands,
    #         prob_regions_beds,
    #         positions_files,
    #     )
    # ]

    # ic(pileup_to_positions_inputs)

    # for inputs in pileup_to_positions_inputs:
    #     pileup_to_positions(*inputs)

    # pileup_to_positions(*pileup_to_positions_inputs[0])

    if processes >= 2:
        pileup_to_positions_processes = max(2, int(processes / 2))
    else:
        pileup_to_positions_processes = 1

    with Pool(processes=pileup_to_positions_processes) as pool:
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
                    remove_non_refbase_noisy_positions,
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

    # 6 - positions dfs -> reads & unique reads dfs

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

    # 7 - reads & unique reads dfs -> proteins & unique proteins dfs

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
    subparsers = parser.add_subparsers(help="sub-command help")

    # define args
    parser.add_argument(
        "--transcriptome",
        required=True,
        type=abs_path_from_str,
        help="A transcritpome reference fasta file. Should have no introns.",
    )
    parser.add_argument(
        "--known_editing_sites",
        type=abs_path_from_str,
        required=True,
        help=(
            "6-col bed file with known A-to-I editing sites. "
            + "Currently, only editing sites on the positive strand are supported."
        ),
    )
    parser.add_argument(
        "--min_percent_of_max_coverage",
        type=float,
        default=0.1,
        help="Keep only positions with coverage >= min_percent_of_max_coverage * max_coverage.",
    )
    parser.add_argument(
        "--include_flags",
        help="Require reads with this flag. 3 is rcommended for PE reads.",
    )
    parser.add_argument(
        "--exclude_flags",
        help=(
            "Exclude reads with this flag. "
            + "The recommended 2304 remove secondary and supplementary (chimeric) alignments."
        ),
    )
    parser.add_argument(
        "--parity", help="Parity of the reads.", choices=["SE", "se", "PE", "pe"]
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
            + "Note base-quality 0 is used as a filtering mechanism for overlap removal. "
            + "Hence using --min-BQ 0 will disable the overlap removal code and act "
            + "as if the --ignore-overlaps option has been set. "
            + "Use this flag for Illumina reads."
        ),
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
        default=3,
        help="Consider the top X noisy positions to define initial noise threshold.",
    )
    parser.add_argument(
        "--assurance_factor",
        type=float,
        default=1.5,
        help=(
            "The noise in a certain position is defined as number of bases of the most abundant alt_base, "
            + "divided by the same number + the number of ref_base bases. "
            + "This is done for all ref_bases that aren't A on the positive strand or T on the negative strand. "
            + "Then, the actual noise level is determined by the "
            + "initial noise threshold (see `top_x_noisy_positions`) * assurance_factor. "
            + "That actual noise level is used to decide wether a position undergoes RNA editing."
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
    parser.add_argument(
        "--override_existing_pileup_files",
        action="store_true",
        help="Recreate existing pileup files from previous runs.",
    )
    parser.add_argument("--gz_compression", action="store_true")
    parser.add_argument(
        "--sample_reads",
        action="store_true",
        help="Whether to use only a sample of reads from the (filtered) reads in the BAM file.",
    )
    parser.add_argument(
        "--num_sampled_reads",
        default=80_000,
        type=int,
        help="Sample `num_sampled_reads` from the (filtered) reads in the BAM file.",
    )
    parser.add_argument(
        "--seed",
        default=1892,
        type=int,
        help="Subsampling seed used to influence which subset of reads is kept.",
    )

    directed_sequencing_subparser = subparsers.add_parser(
        "directed_sequencing_data",
        help="Original use of this program by suplying a data table",
    )
    directed_sequencing_subparser.set_defaults(func=directed_sequencing_main)
    directed_sequencing_subparser.add_argument(
        "--data_table", type=abs_path_from_str, required=True
    )
    directed_sequencing_subparser.add_argument(
        "--data_table_sep", default=",", help="Delimiter used in `data_table`."
    )
    directed_sequencing_subparser.add_argument(
        "--sample_col", default="Sample", help="Sample col label in `data_table`."
    )
    directed_sequencing_subparser.add_argument(
        "--group_col", default="Gene", help="Group col label in `data_table`."
    )
    directed_sequencing_subparser.add_argument(
        "--region_col",
        default="Region",
        help=(
            "Region col label in `data_table`."
            "It's values should be given in a samtools-compatible format, i.e., either `chrom` of `chrom:start-end`, "
            "where `start` is 1-based and `end` - inclusive. "
            "Should represent an ORF (or part of it) in the transcriptome, so its length should be a multiple of 3."
        ),
    )
    directed_sequencing_subparser.add_argument(
        "--start_col",
        default="Start",
        help="Start col label in `data_table`. Its values should be 0-based (as in BED format).",
    )
    directed_sequencing_subparser.add_argument(
        "--end_col",
        default="End",
        help="End col label in `data_table`. Its values should be inclusive (as in BED format).",
    )
    directed_sequencing_subparser.add_argument(
        "--strand_col", default="Strand", help="Strand col label in `data_table`."
    )
    directed_sequencing_subparser.add_argument(
        "--path_col", default="Path", help="Path col label in `data_table`."
    )
    directed_sequencing_subparser.add_argument(
        "--cds_regions",
        type=abs_path_from_str,
        help=(
            "6-col bed file of coding regions. "
            "Currently, only positions insided those regions are used for isoforms analysis. "
            "Pay attention! In its abscence, all positions are assumed to be within coding regions!"
        ),
    )
    directed_sequencing_subparser.add_argument(
        "--prob_regions_bed_col",
        default="ProbRegionsBED",
        help="Problematic regions to exclude col label in `data_table`. Its values should be paths to BED files.",
    )

    undirected_sequencing_subparser = subparsers.add_parser(
        "undirected_sequencing_data",
        help=(
            "Use multiple BAMs, each aligned to different trasncript, of multiple samples, "
            "using either predefined editing sites and/or denovo detected ones."
        ),
    )
    undirected_sequencing_subparser.set_defaults(func=undirected_sequencing_main)
    undirected_sequencing_subparser.add_argument(
        "--alignments_stats_table",
        help="By-chrom-by-sample alignments stats file. This table also includes the number of known editing sites per chrom.",
        required=True,
        type=abs_path_from_str,
    )
    undirected_sequencing_subparser.add_argument(
        "--alignments_stats_table_sep",
        default="\t",
        help="Delimiter of the by-chrom-by-sample alignments stats table.",
    )
    undirected_sequencing_subparser.add_argument(
        "--min_samples",
        default=0,
        help="Use only chroms with at least that many samples mapped to them.",
        type=int,
    )
    undirected_sequencing_subparser.add_argument(
        "--min_mapped_reads_per_sample",
        default=0,
        help="Use only chroms with whose mean coverage is at least that.",
        type=float,
    )
    undirected_sequencing_subparser.add_argument(
        "--total_mapped_reads",
        default=50,
        help="Use only chroms with with total (pooled) coverage of at least that.",
        type=int,
    )
    undirected_sequencing_subparser.add_argument(
        "--min_mapped_reads_per_position",
        type=int,
        default=None,
        help=(
            "For each position pooled from multiple samples, keep only positions with "
            "such absolute coverage. "
            "(To be effective, it should probably be lower than `total_mapped_reads`.) "
            "By default, the minimal coverage is relative, and determined by the "
            "`min_percent_of_max_coverage` parameter."
        ),
    )
    undirected_sequencing_subparser.add_argument(
        "--min_known_sites",
        default=0,
        help="Use only chroms with at least that many known editing sites in them.",
        type=int,
    )
    undirected_sequencing_subparser.add_argument(
        "--pooled_transcript_noise_threshold",
        default=0.06,
        help=(
            "Use only transcripts whose pooled_transcript_noise "
            "(same as noise level before the multiplication by assurance factor) "
            "is below `pooled_transcript_noise_threshold`."
        ),
        type=float,
    )
    undirected_sequencing_subparser.add_argument(
        "--main_by_chrom_dir",
        type=abs_path_from_str,
        help="Folder whose each of its subfolders contain BAMs of different samples aligned to a certain chrom.",
    )
    undirected_sequencing_subparser.add_argument("--postfix", default=".bam")
    undirected_sequencing_subparser.add_argument(
        "--interfix_start",
        default=".aligned.sorted",
        help=(
            "For a BAM named `SRR17321899.aligned.sorted.comp22_c1_seq1`, with a given postfix (see above), "
            "the interfix_start allows to get the sample's name."
        ),
    )
    undirected_sequencing_subparser.add_argument(
        "--cds_regions",
        type=abs_path_from_str,
        help="6-col bed file of coding regions - used to define ORFs' boundries.",
        required=True,
    )
    undirected_sequencing_subparser.add_argument(
        "--samples_table", type=abs_path_from_str, required=True
    )
    undirected_sequencing_subparser.add_argument(
        "--samples_table_sep", default=",", help="Delimiter used in `samples_table`."
    )
    undirected_sequencing_subparser.add_argument(
        "--sample_col", default="Sample", help="Sample col label in `samples_table`."
    )
    undirected_sequencing_subparser.add_argument(
        "--group_col", default="Tissue", help="Group col label in `samples_table`."
    )
    undirected_sequencing_subparser.add_argument(
        "--remove_non_refbase_noisy_positions",
        action="store_true",
        help="Remove non ref-base positions whose noise are high.",
    )
    undirected_sequencing_subparser.add_argument(
        "--known_sites_only",
        action="store_true",
        help="Allow only known sites to be considered as edited, rather than denovo ones too.",
    )
    undirected_sequencing_subparser.add_argument(
        "--final_editing_scheme",
        choices=["BH only", "BH after noise thresholding"],
        required=True,
        help=(
            "The final editing scheme to use. "
            "It can either depend only on corrected editing p-values across the CDS of the whole transcriptome, "
            "or also depend on the noise thresholding (which is corrected for multiple testing in the same manner)."
        ),
    )
    undirected_sequencing_subparser.add_argument(
        "--alternative_hypothesis",
        choices=["larger", "smaller", "two-sided"],
        required=True,
        help="Alternative hypothesis for the binomial test for noise and editing.",
    )
    undirected_sequencing_subparser.add_argument(
        "--pileup_dir_name", default="PileupFiles"
    )
    undirected_sequencing_subparser.add_argument(
        "--positions_dir_name", default="PositionsFiles"
    )
    undirected_sequencing_subparser.add_argument(
        "--reads_dir_name", default="ReadsFiles"
    )
    undirected_sequencing_subparser.add_argument(
        "--proteins_dir_name", default="ProteinsFiles"
    )
    undirected_sequencing_subparser.add_argument(
        "--keep_bam_files",
        action="store_true",
        help=(
            "Keep the filtered and/or sampled bam files. "
            "Altough this option perhaps could be suitable to the directed data as well, "
            "in practice, the number of directed samples is relatively small."
        ),
    )
    undirected_sequencing_subparser.add_argument(
        "--override_existing_bam_files",
        action="store_true",
        help="Recreate existing filtered and/or sampled bam files from previous runs.",
    )

    # subparsers.default = "default"

    return parser


if __name__ == "__main__":
    # run

    # main(
    #     **vars(define_and_parse_args())
    # )  # https://stackoverflow.com/a/35824590/10249633 argparse.Namespace -> Dict

    parser = define_args()
    args = parser.parse_args()
    ic(args)
    # ic(args)
    args.func(
        **vars(args)
    )  # https://stackoverflow.com/a/35824590/10249633 argparse.Namespace -> Dict
    # main(
    #     **vars(args)
    # )  # https://stackoverflow.com/a/35824590/10249633 argparse.Namespace -> Dict

    # end
    final_words()
