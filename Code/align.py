import argparse
import subprocess
from multiprocessing import Pool
from pathlib import Path
from typing import Sequence
from typing import Union

import pysam
import pandas as pd
from pybedtools import BedTool

from General.multiqc import multiqc
from General.consts import final_words
from General.os_utils import (
    extract_sample_name,
    find_files,
    group_pe_fastq_files,
    decompress,
    delete_folder_with_files,
    copy_bytes,
)
from General.argparse_utils import abs_path_from_str, expanded_path_from_str
from Alignment.alignment_utils import count_reads
from EditingUtils.gff3 import read_gff
from General.pandas_utils import reorder_df_by_wanted_cols

# from Alignment.alignment_utils import samtools_statistics


def samtools_statistics(
    samtools_path: Path,
    alignment_file: Path,
    sample_name: str = None,
    out_dir: Path = None,
    programs: Sequence[str] = ("flagstat", "idxstats", "stats"),
):
    """Generate statistics for alignment file."""
    sample_name = (
        sample_name
        if sample_name
        else extract_sample_name(sample=alignment_file, postfix=".bam")
    )
    out_dir = out_dir if out_dir else alignment_file.parent
    for program in programs:
        cmd = f"{samtools_path} {program} {alignment_file} > {Path(out_dir, f'{sample_name}.{program}')}"
        subprocess.run(cmd, shell=True)


def pacbio_index_genome(
    base_conda_env_dir: Path,
    pb_conda_env_name: str,
    pbmm2_path: Path,
    genome: Path,
    genome_index_file: Path,
    preset: str,
    threads: int,
):
    index_cmd = (
        f". {Path(base_conda_env_dir, 'etc/profile.d/conda.sh')} && "
        f"conda activate {pb_conda_env_name} && "
        f"{pbmm2_path} "
        "index "
        f"--preset {preset} "
        f"--num-threads {threads} "
        f"{genome} "
        f"{genome_index_file} "
    )
    subprocess.run(index_cmd, shell=True, executable="/bin/bash")


def bwa_index_genome(
    samtools_path: Path,
    bwa_path: Path,
    genome: Path,
    index_dir: Path,
):
    if index_dir.exists():
        delete_folder_with_files(index_dir)
    index_dir.mkdir()
    linked_genome = Path(index_dir, genome.name)
    linked_genome.symlink_to(genome)
    samtools_faidx_cmd = f"{samtools_path} faidx {linked_genome}"
    subprocess.run(samtools_faidx_cmd, shell=True, cwd=index_dir)
    bwa_index_cmd = f"{bwa_path} index {linked_genome}"
    subprocess.run(bwa_index_cmd, shell=True, cwd=index_dir)
    genome_index_file = Path(index_dir, genome.name)
    return genome_index_file


def bwa_align_pe_reads(
    bwa_path: Path,
    threads: int,
    genome_index_file: Path,
    left_in_file: Path,
    right_in_file: Path,
    out_dir: Path,
    sample_name: Path,
    samtools_path: Path,
    require_flags: int = 3,
    exclude_flags: int = 2304,
    separate_by_chrom: bool = True,
):
    sam = Path(out_dir, f"{sample_name}.sam")
    align_cmd = (
        f"{bwa_path} "
        "mem "
        f"-t {threads} "
        "-M "
        f"{genome_index_file} "
        f"{left_in_file} "
        f"{right_in_file} "
        f"> {sam} "
    )
    subprocess.run(align_cmd, shell=True)
    bam = Path(out_dir, f"{sample_name}.bam")
    sam_to_bam_cmd = f"{samtools_path} view -@ {threads} -h -o {bam} {sam}"
    subprocess.run(sam_to_bam_cmd, shell=True)
    sam.unlink()
    sorted_bam = Path(out_dir, f"{sample_name}.sorted.bam")
    sort_cmd = f"{samtools_path} sort -@ {threads} -o {sorted_bam} {bam}"
    subprocess.run(sort_cmd, shell=True)
    bam.unlink()
    aligned_filtered_sorted_bam = Path(
        out_dir, f"{sample_name}.sorted.aligned.filtered.bam"
    )
    filter_cmd = (
        f"{samtools_path} view "
        f"-@ {threads} "
        "-h "
        f"-f {require_flags} "
        f"-F {exclude_flags} "
        f"-o {aligned_filtered_sorted_bam} "
        f"{sorted_bam}"
    )
    subprocess.run(filter_cmd, shell=True)
    index_cmd = f"{samtools_path} index {aligned_filtered_sorted_bam}"
    subprocess.run(index_cmd, shell=True, cwd=out_dir)

    if separate_by_chrom:
        separte_bam_to_chroms(
            samtools_path,
            out_dir,
            threads,
            sample_name,
            bam=aligned_filtered_sorted_bam,
            interfix=".sorted.aligned.filtered",
        )

    # if not separate_by_chrom:
    #     return
    # by_chrom_dir = Path(out_dir, f"{sample_name}.ByChrom")
    # by_chrom_dir.mkdir(exist_ok=True)
    # chroms = {
    #     record.reference_name
    #     for record in pysam.AlignmentFile(aligned_filtered_sorted_bam, "rb")
    # }
    # for chrom in chroms:
    #     bam_in_region = Path(
    #         by_chrom_dir, f"{sample_name}.sorted.aligned.filtered.{chrom}.bam"
    #     )
    #     filter_cmd = f"{samtools_path} view -@ {threads} -h -o {bam_in_region} {aligned_filtered_sorted_bam} {chrom}"
    #     subprocess.run(filter_cmd, shell=True)
    #     index_cmd = f"{samtools_path} index {bam_in_region}"
    #     subprocess.run(index_cmd, shell=True, cwd=by_chrom_dir)

    # # generate statistics per sample
    # samtools_statistics(samtools_path, out_file)


def separte_bam_to_chroms(
    samtools_path: Path,
    out_dir: Path,
    threads: int,
    sample_name: str,
    bam: Path,
    interfix: str,  # https://www.wikiwand.com/en/Interfix
):
    by_chrom_dir = Path(out_dir, f"{sample_name}.ByChrom")
    by_chrom_dir.mkdir(exist_ok=True)
    chroms = {record.reference_name for record in pysam.AlignmentFile(bam, "rb")}
    for chrom in chroms:
        bam_in_region = Path(by_chrom_dir, f"{sample_name}{interfix}.{chrom}.bam")
        filter_cmd = (
            f"{samtools_path} view -@ {threads} -h -o {bam_in_region} {bam} {chrom}"
        )
        subprocess.run(filter_cmd, shell=True)
        # it is somehow possible (why?) to get chrom with no reads aligned to it --> the bam_in_region is not created,
        # and, naturally, we would like to index only existing bam files...
        if bam_in_region.exists():
            index_cmd = f"{samtools_path} index {bam_in_region}"
            subprocess.run(index_cmd, shell=True, cwd=by_chrom_dir)


def get_nonoverlapping_gene_regions(gff_file: Path) -> list[str]:
    """
    Return a list of nonoverlapping gene regions, whether they overlap on the same or on opposing strands,
    formatted according to samtools.
    """
    gff_df = read_gff(gff_file)
    genes_gff_df = gff_df.loc[gff_df["Type"] == "gene"]
    genes_gff_df

    genes_bed_df = genes_gff_df.loc[:, ["Chrom", "Start", "End", "Strand"]]
    genes_bed_df["Start"] = genes_bed_df["Start"].astype(int)
    genes_bed_df["Start"] = genes_bed_df["Start"] - 1
    genes_bed_df["End"] = genes_bed_df["End"].astype(int)
    genes_bed_df.insert(genes_bed_df.columns.get_loc("End") + 1, "Name", ".")
    genes_bed_df.insert(genes_bed_df.columns.get_loc("Name") + 1, "Score", ".")

    genes_bedtool = BedTool().from_dataframe(genes_bed_df).sort()
    merged_genes_bedtool = genes_bedtool.merge(c=1, o="count")
    overlapping_merged_genes_bedtool = BedTool(
        [feature[:3] for feature in merged_genes_bedtool if int(feature[3]) > 1]
    )
    non_overlapping_genes_bedtool = genes_bedtool.intersect(
        overlapping_merged_genes_bedtool, v=True
    )

    non_overlapping_gene_regions = [
        f"{gene.chrom}:{gene.start+1}-{gene.end}"
        for gene in non_overlapping_genes_bedtool
    ]

    return non_overlapping_gene_regions


# todo merge separte_bam_to_chroms + separate_bam_to_genes --> separate_bam_to_regions
def separate_bam_to_genes(
    samtools_path: Path,
    out_dir: Path,
    threads: int,
    sample_name: str,
    bam: Path,
    interfix: str,  # https://www.wikiwand.com/en/Interfix
    gene_regions: list[str],  # preferably nonoverlapping
    include_flags: Union[int, str, None],
    exclude_flags: Union[int, str, None],
):
    by_gene_dir = Path(out_dir, f"{sample_name}.ByGene")
    by_gene_dir.mkdir(exist_ok=True)

    for gene_region in gene_regions:

        mapped_reads_count = count_reads(
            samtools_path, bam, gene_region, include_flags, exclude_flags, threads
        )
        if mapped_reads_count == 0:
            continue

        bam_in_region = Path(by_gene_dir, f"{sample_name}{interfix}.{gene_region}.bam")
        filter_cmd = f"{samtools_path} view -@ {threads} -h -o {bam_in_region} {bam} {gene_region}"
        subprocess.run(filter_cmd, shell=True)
        # it is somehow possible (why?) to get chrom with no reads aligned to it --> the bam_in_region is not created,
        # and, naturally, we would like to index only existing bam files...
        if bam_in_region.exists():
            index_cmd = f"{samtools_path} index {bam_in_region}"
            subprocess.run(index_cmd, shell=True, cwd=by_gene_dir)


def pacbio_align(
    base_conda_env_dir: Path,
    pb_conda_env_name: str,
    pbmm2_path: Path,
    preset: str,
    best_n_alignments_per_read: int,
    threads: int,
    genome_index_file: Path,
    in_file: Path,
    out_file: Path,
    samtools_path: Path,
    out_dir: Path,
    sample_name: str,
    interfix: str,
    gene_regions: Union[list[str], None],  # preferably nonoverlapping
    include_flags: Union[int, str, None],
    exclude_flags: Union[int, str, None],
    separate_by_chrom: bool = False,
    seperate_by_gene: bool = False,
):
    # align
    align_cmd = (
        f". {Path(base_conda_env_dir, 'etc/profile.d/conda.sh')} && "
        f"conda activate {pb_conda_env_name} && "
        f"{pbmm2_path} "
        "align "
        "--sort "
        f"--preset {preset} "
        f"--best-n {best_n_alignments_per_read} "
        f"--num-threads {threads} "
        f"{genome_index_file} "
        f"{in_file} "
        f"{out_file} "
    )
    subprocess.run(align_cmd, shell=True, executable="/bin/bash")
    # generate statistics per whole sample
    samtools_statistics(samtools_path, out_file)

    if separate_by_chrom:
        separte_bam_to_chroms(
            samtools_path,
            out_dir,
            threads,
            sample_name,
            bam=out_file,
            interfix=interfix,
        )

    if seperate_by_gene:
        if gene_regions is None:
            raise TypeError(
                "gene_regions is None but should be list[str] when seperate_by_gene is True"
            )
        separate_bam_to_genes(
            samtools_path,
            out_dir,
            threads,
            sample_name,
            bam=out_file,
            interfix=interfix,
            gene_regions=gene_regions,
            include_flags=include_flags,
            exclude_flags=exclude_flags,
        )


def illumina_main(
    *,
    genome: Path,
    in_dir: Path,
    postfix: str,
    mate_prefix: str,
    recursive: bool,
    out_dir: Path,
    processes: int,
    threads: int,
    samtools_path: Path,
    decompress_cmd: str,
    bwa_path: Path,
    require_flags: int,
    exclude_flags: int,
    separate_by_chrom: bool,
    **kwargs,
):
    out_dir.mkdir(exist_ok=True)

    index_dir = Path(out_dir, "GenomeIndex")
    genome_index_file = bwa_index_genome(samtools_path, bwa_path, genome, index_dir)

    compressed_in_files = find_files(in_dir, postfix, recursive)

    decompressed_temp_dir = Path(out_dir, "Decompressed.Temp")
    decompressed_temp_dir.mkdir(exist_ok=True)

    decompressed_in_files = [
        Path(decompressed_temp_dir, compressed_in_file.stem)
        for compressed_in_file in compressed_in_files
    ]

    with Pool(processes=processes) as pool:
        pool.starmap(
            func=decompress,
            iterable=[
                (decompress_cmd, compressed_file, decompressed_file)
                for compressed_file, decompressed_file in zip(
                    compressed_in_files, decompressed_in_files
                )
            ],
        )

    final_postfix = postfix.split(".")[-1]
    main_postfix = postfix.removesuffix(f".{final_postfix}")

    paired_decompressed_in_files = group_pe_fastq_files(
        decompressed_in_files, main_postfix, mate_prefix
    )

    paired_samples_names = [
        extract_sample_name(left_file, f"{mate_prefix}1{main_postfix}")
        for left_file, _ in paired_decompressed_in_files
    ]

    with Pool(processes=processes) as pool:
        pool.starmap(
            func=bwa_align_pe_reads,
            iterable=[
                (
                    bwa_path,
                    threads,
                    genome_index_file,
                    left_file,
                    right_file,
                    out_dir,
                    sample_name,
                    samtools_path,
                    require_flags,
                    exclude_flags,
                    separate_by_chrom,
                )
                for (left_file, right_file), sample_name in zip(
                    paired_decompressed_in_files, paired_samples_names
                )
            ],
        )

    delete_folder_with_files(decompressed_temp_dir)


def whole_transcriptome_isoseq_main(
    *,
    genome: Path,
    # gff: Path,
    known_sites_bed_file: Path,
    in_dir: Path,
    postfix: str,
    recursive: bool,
    out_dir: Path,
    base_conda_env_dir: Path,
    pb_conda_env_name: str,
    pbmm2_path: Path,
    best_n_alignments_per_read: int,
    processes: int,
    threads: int,
    samtools_path: Path,
    multiqc_path: Path,
    # include_flags: Union[int, str, None],
    # exclude_flags: Union[int, str, None],
    **kwargs,
):
    preset = "ISOSEQ"
    interfix = ".aligned.sorted"

    out_dir.mkdir(exist_ok=True)

    in_files = find_files(in_dir, postfix, recursive)
    samples_names = [extract_sample_name(in_file, postfix) for in_file in in_files]

    genome_index_file = Path(out_dir, f"{genome.name}.mmi")

    # index genome
    pacbio_index_genome(
        base_conda_env_dir,
        pb_conda_env_name,
        pbmm2_path,
        genome,
        genome_index_file,
        preset,
        threads,
    )

    # gene_regions = get_nonoverlapping_gene_regions(gff)
    gene_regions = None
    separate_by_chrom = True
    seperate_by_gene = False
    include_flags = None
    exclude_flags = None

    # align
    with Pool(processes=processes) as pool:
        pool.starmap(
            func=pacbio_align,
            iterable=[
                (
                    base_conda_env_dir,
                    pb_conda_env_name,
                    pbmm2_path,
                    preset,
                    best_n_alignments_per_read,
                    threads,
                    genome_index_file,
                    in_file,
                    Path(out_dir, f"{sample_name}{interfix}.bam"),
                    samtools_path,
                    out_dir,
                    sample_name,
                    interfix,
                    gene_regions,
                    include_flags,
                    exclude_flags,
                    separate_by_chrom,
                    seperate_by_gene,
                )
                for in_file, sample_name in zip(in_files, samples_names)
            ],
        )

    # gather separated-by-chrom bams
    # first to main sub-dir
    # and then subdir for chrom, containing bams from different samples
    by_chrom_samples_dirs = [
        Path(out_dir, f"{sample_name}.ByChrom") for sample_name in samples_names
    ]
    chroms = {
        f.suffixes[-2].removeprefix(".")
        for by_chrom_sample_dir in by_chrom_samples_dirs
        for f in by_chrom_sample_dir.iterdir()
        if f.suffix.endswith("bam")
    }

    main_by_chrom_dir = Path(out_dir, "ByChrom")
    main_by_chrom_dir.mkdir(exist_ok=True)

    with Pool(processes=processes) as pool:
        by_chrom_dfs = pool.starmap(
            func=gather_by_chrom_bams_and_collect_stats,
            iterable=[
                (
                    samtools_path,
                    main_by_chrom_dir,
                    chrom,
                    samples_names,
                    by_chrom_samples_dirs,
                    interfix,
                    threads,
                    include_flags,
                    exclude_flags,
                )
                for chrom in chroms
            ],
        )

    # write mapping stats
    df = pd.concat(by_chrom_dfs, ignore_index=True)

    agg_df = (
        df.groupby("Chrom")
        .agg({"MappedReads": sum, "Sample": len})
        .reset_index()
        .rename(columns={"Sample": "Samples"})
    )
    agg_df["MappedReadsPerSample"] = agg_df["MappedReads"] / agg_df["Samples"]

    # known_sites_bed_file = "/private7/projects/Combinatorics/O.vulgaris/Annotations/O.vul.EditingSites.bed"
    # min_known_sites = 5
    known_sites_df = pd.read_csv(
        known_sites_bed_file,
        sep="\t",
        names="Chrom Start End Name Score Strand".split(),
        comment="#",
    )
    known_sites_per_chrom_df = (
        known_sites_df.groupby("Chrom")
        .size()
        .reset_index()
        .rename(columns={0: "KnownSites"})
        .sort_values("KnownSites", ascending=False)
    )

    agg_df = agg_df.merge(known_sites_per_chrom_df, how="left")
    agg_df["KnownSites"] = agg_df["KnownSites"].fillna(0)

    agg_df = agg_df.sort_values(
        ["KnownSites", "MappedReadsPerSample", "MappedReads", "Samples"],
        ascending=False,
    ).reset_index(drop=True)

    df.to_csv(Path(out_dir, "ByChromBySampleSummary.tsv"), sep="\t", index=False)
    agg_df.to_csv(
        Path(out_dir, "AggregatedByChromBySampleSummary.tsv"), sep="\t", index=False
    )

    # run MultiQC on the reports created for each whole aligned file
    multiqc(multiqc_path=multiqc_path, data_dir=out_dir)

    # delete the genome's index (it's quite heavy)
    genome_index_file.unlink(missing_ok=True)


# def pacbio_preprocessed_isoseq_main(
#     *,
#     genome: Path,
#     gff: Path,
#     in_dir: Path,
#     postfix: str,
#     recursive: bool,
#     out_dir: Path,
#     base_conda_env_dir: Path,
#     pb_conda_env_name: str,
#     pbmm2_path: Path,
#     best_n_alignments_per_read: int,
#     processes: int,
#     threads: int,
#     samtools_path: Path,
#     multiqc_path: Path,
#     # include_flags: Union[int, str, None],
#     # exclude_flags: Union[int, str, None],
#     **kwargs,
# ):
#     preset = "ISOSEQ"
#     interfix = ".aligned.sorted"

#     out_dir.mkdir(exist_ok=True)

#     in_files = find_files(in_dir, postfix, recursive)
#     samples_names = [extract_sample_name(in_file, postfix) for in_file in in_files]

#     genome_index_file = Path(out_dir, f"{genome.name}.mmi")

#     # index genome
#     pacbio_index_genome(
#         base_conda_env_dir,
#         pb_conda_env_name,
#         pbmm2_path,
#         genome,
#         genome_index_file,
#         preset,
#         threads,
#     )

#     gene_regions = get_nonoverlapping_gene_regions(gff)
#     separate_by_chrom = False
#     seperate_by_gene = True
#     include_flags = None
#     exclude_flags = None

#     # align
#     with Pool(processes=processes) as pool:
#         pool.starmap(
#             func=pacbio_align,
#             iterable=[
#                 (
#                     base_conda_env_dir,
#                     pb_conda_env_name,
#                     pbmm2_path,
#                     preset,
#                     best_n_alignments_per_read,
#                     threads,
#                     genome_index_file,
#                     in_file,
#                     Path(out_dir, f"{sample_name}{interfix}.bam"),
#                     samtools_path,
#                     out_dir,
#                     sample_name,
#                     interfix,
#                     gene_regions,
#                     include_flags,
#                     exclude_flags,
#                     separate_by_chrom,
#                     seperate_by_gene,
#                 )
#                 for in_file, sample_name in zip(in_files, samples_names)
#             ],
#         )

#     be_gene_samples_dirs = [
#         Path(out_dir, f"{sample_name}.ByGene") for sample_name in samples_names
#     ]

#     main_by_gene_dir = Path(out_dir, "ByGene")
#     main_by_gene_dir.mkdir(exist_ok=True)

#     with Pool(processes=processes) as pool:
#         by_gene_dfs = pool.starmap(
#             func=gather_by_gene_bams_and_collect_stats,
#             iterable=[
#                 (
#                     samtools_path,
#                     main_by_gene_dir,
#                     gene_region,
#                     samples_names,
#                     be_gene_samples_dirs,
#                     interfix,
#                     threads,
#                     include_flags,
#                     exclude_flags,
#                 )
#                 for gene_region in gene_regions
#             ],
#         )

#     # write mapping stats
#     by_gene_df = pd.concat(by_gene_dfs, ignore_index=True)

#     chrom_startend = by_gene_df["GeneRegion"].str.split(":", expand=True)
#     chrom = chrom_startend.iloc[:, 0]
#     start_end = chrom_startend.iloc[:, 1].str.split("-", expand=True)
#     start = start_end.iloc[:, 0].astype(int)
#     end = start_end.iloc[:, 1].astype(int)

#     by_gene_df.insert(by_gene_df.columns.get_loc("GeneRegion") + 1, "Chrom", chrom)
#     by_gene_df.insert(by_gene_df.columns.get_loc("GeneRegion") + 2, "Start", start)
#     by_gene_df.insert(by_gene_df.columns.get_loc("GeneRegion") + 3, "End", end)

#     by_gene_df = by_gene_df.sort_values(
#         ["Chrom", "Start", "End", "Sample"]
#     ).reset_index(drop=True)

#     mapped_genes_bed3_df = by_gene_df.loc[:, ["Chrom", "Start", "End"]].drop_duplicates(
#         ignore_index=True
#     )
#     mapped_genes_bed3_df["Start"] = mapped_genes_bed3_df["Start"] - 1
#     mapped_genes_bed3_bedtool = BedTool().from_dataframe(mapped_genes_bed3_df).sort()

#     gff_df = read_gff(gff)

#     exons_gff_df = gff_df.loc[gff_df["Type"] == "exon"]
#     exons_bed6_df = exons_gff_df.loc[:, ["Chrom", "Start", "End", "Score", "Strand"]]
#     exons_bed6_df["Start"] = exons_bed6_df["Start"].astype(int) - 1
#     exons_bed6_df["End"] = exons_bed6_df["End"].astype(int)
#     exons_bed6_df.insert(exons_bed6_df.columns.get_loc("End") + 1, "Name", ".")
#     exons_bed6_bedtools = BedTool().from_dataframe(exons_bed6_df)

#     genes_gff_df = gff_df.loc[gff_df["Type"] == "gene"]
#     genes_bed_df = genes_gff_df.loc[:, ["Chrom", "Start", "End", "Strand"]]
#     genes_bed_df["Start"] = genes_bed_df["Start"].astype(int)
#     genes_bed_df["Start"] = genes_bed_df["Start"] - 1
#     genes_bed_df["End"] = genes_bed_df["End"].astype(int)
#     genes_bed_df.insert(genes_bed_df.columns.get_loc("End") + 1, "Name", ".")
#     genes_bed_df.insert(genes_bed_df.columns.get_loc("Name") + 1, "Score", ".")
#     genes_bedtool = BedTool().from_dataframe(genes_bed_df).sort()

#     mapped_genes_bed6_bedtool = genes_bedtool.intersect(
#         mapped_genes_bed3_bedtool
#     ).intersect(  # retain only mapped genes
#         exons_bed6_bedtools, s=True, c=True
#     )  # count number in exons in each such gene

#     mapped_genes_df = (
#         mapped_genes_bed6_bedtool.to_dataframe()
#         .drop(["name", "score"], axis=1)
#         .rename(
#             columns={
#                 "chrom": "Chrom",
#                 "start": "Start",
#                 "end": "End",
#                 "strand": "Strand",
#                 "thickStart": "ExonsInGene",
#             }
#         )
#     )
#     mapped_genes_df["Start"] = mapped_genes_df["Start"] + 1

#     by_gene_df = by_gene_df.merge(mapped_genes_df, on=["Chrom", "Start", "End"])
#     by_gene_df = reorder_df_by_wanted_cols(
#         by_gene_df, wanted_first_cols=["GeneRegion", "Chrom", "Start", "End", "Strand"]
#     )

#     agg_df = (
#         by_gene_df.groupby(
#             ["GeneRegion", "Chrom", "Start", "End", "Strand", "ExonsInGene"]
#         )
#         .agg({"MappedReads": sum, "Sample": len})
#         .reset_index()
#         .rename(columns={"Sample": "Samples"})
#     )
#     agg_df["MappedReadsPerSample"] = agg_df["MappedReads"] / agg_df["Samples"]
#     agg_df["CodingGene"] = agg_df["ExonsInGene"] > 0
#     agg_df = agg_df.sort_values(
#         ["CodingGene", "MappedReadsPerSample", "MappedReads", "Samples"],
#         ascending=False,
#     ).reset_index(drop=True)

#     by_gene_df.to_csv(
#         Path(out_dir, "ByGeneRegionBySampleSummary.tsv"), sep="\t", index=False
#     )
#     agg_df.to_csv(
#         Path(out_dir, "AggregatedByGeneBySampleSummary.tsv"), sep="\t", index=False
#     )

#     # # gather separated-by-chrom bams
#     # # first to main sub-dir
#     # # and then subdir for chrom, containing bams from different samples
#     # by_chrom_samples_dirs = [
#     #     Path(out_dir, f"{sample_name}.ByChrom") for sample_name in samples_names
#     # ]
#     # chroms = {
#     #     f.suffixes[-2].removeprefix(".")
#     #     for by_chrom_sample_dir in by_chrom_samples_dirs
#     #     for f in by_chrom_sample_dir.iterdir()
#     #     if f.suffix.endswith("bam")
#     # }

#     # main_by_chrom_dir = Path(out_dir, "ByChrom")
#     # main_by_chrom_dir.mkdir(exist_ok=True)

#     # # include_flags = None
#     # # exclude_flags = None
#     # with Pool(processes=processes) as pool:
#     #     by_chrom_dfs = pool.starmap(
#     #         func=gather_by_chrom_bams_and_collect_stats,
#     #         iterable=[
#     #             (
#     #                 samtools_path,
#     #                 main_by_chrom_dir,
#     #                 chrom,
#     #                 samples_names,
#     #                 by_chrom_samples_dirs,
#     #                 interfix,
#     #                 threads,
#     #                 include_flags,
#     #                 exclude_flags,
#     #             )
#     #             for chrom in chroms
#     #         ],
#     #     )

#     # # write mapping stats
#     # df = pd.concat(by_chrom_dfs, ignore_index=True)

#     # agg_df = (
#     #     df.groupby("Chrom")
#     #     .agg({"MappedReads": sum, "Sample": len})
#     #     .reset_index()
#     #     .rename(columns={"Sample": "Samples"})
#     # )
#     # agg_df["MappedReadsPerSample"] = agg_df["MappedReads"] / agg_df["Samples"]

#     # gff_df = read_gff(gff)
#     # exon_conatining_chroms = gff_df.loc[gff_df["Type"] == "exon"][
#     #     "Chrom"
#     # ].drop_duplicates()

#     # agg_df = pd.merge(
#     #     agg_df, exon_conatining_chroms, on="Chrom", how="left", indicator=True
#     # ).rename(columns={"_merge": "InExon"})
#     # agg_df["InExon"] = (
#     #     agg_df["InExon"].replace({"both": True, "left_only": False}).astype(bool)
#     # )
#     # agg_df = agg_df.sort_values(
#     #     ["InExon", "MappedReadsPerSample", "MappedReads", "Samples"], ascending=False
#     # ).reset_index(drop=True)

#     # df.to_csv(Path(out_dir, "ByChromBySampleSummary.tsv"), sep="\t", index=False)
#     # agg_df.to_csv(
#     #     Path(out_dir, "AggregatedByChromBySampleSummary.tsv"), sep="\t", index=False
#     # )

#     # run MultiQC on the reports created for each whole aligned file
#     multiqc(multiqc_path=multiqc_path, data_dir=out_dir)

#     # delete the genome's index (it's quite heavy)
#     genome_index_file.unlink(missing_ok=True)


def gather_by_gene_bams_and_collect_stats(
    samtools_path: Path,
    main_by_gene_dir: Path,
    gene_region: str,
    samples_names: list[str],
    by_gene_samples_dirs: list[Path],
    interfix: str,
    threads: int,
    include_flags: Union[int, str, None],
    exclude_flags: Union[int, str, None],
):
    # stats containers
    _samples_names = []
    _mapped_reads = []
    # the per-gene dir that will hold the per-sample-per-gene bams
    gene_dir = Path(main_by_gene_dir, gene_region)
    gene_dir_created = False
    for sample_name, by_gene_sample_dir in zip(samples_names, by_gene_samples_dirs):
        bam_in_region = Path(
            by_gene_sample_dir, f"{sample_name}{interfix}.{gene_region}.bam"
        )
        # not each sample has reads mapped to this specific gene region
        if bam_in_region.exists():
            # create this dir only if necessary - maybe no sample has reads mapped to that region
            if not gene_dir_created:
                gene_dir.mkdir(exist_ok=True)
                gene_dir_created = True

            # # copy the per-sample-per-gene bam to the per-gene dir
            # bam_in_region_copy = Path(chrom_dir, bam_in_region.name)
            # copy_bytes(bam_in_region, bam_in_region_copy)

            # soft-link the per-sample-per-gene bam to the per-gene dir
            bam_in_region_link = Path(gene_dir, bam_in_region.name)
            bam_in_region_link.symlink_to(bam_in_region)
            # also link the .bai index file
            bam_in_region_index = Path(
                bam_in_region.parent, f"{bam_in_region.name}.bai"
            )
            bam_in_region_index_link = Path(gene_dir, bam_in_region_index.name)
            bam_in_region_index_link.symlink_to(bam_in_region_index)

            # update stats
            _samples_names.append(sample_name)
            _mapped_reads.append(
                count_reads(
                    samtools_path,
                    bam_in_region,
                    None,
                    include_flags,
                    exclude_flags,
                    threads,
                )
            )
    df = pd.DataFrame(
        {
            "GeneRegion": gene_region,
            "Sample": _samples_names,
            "MappedReads": _mapped_reads,
        }
    )
    return df


def gather_by_chrom_bams_and_collect_stats(
    samtools_path: Path,
    main_by_chrom_dir: Path,
    chrom: str,
    samples_names: list[str],
    by_chrom_samples_dirs: list[Path],
    interfix: str,
    threads: int,
    include_flags: Union[int, str, None],
    exclude_flags: Union[int, str, None],
):
    # stats containers
    _samples_names = []
    _mapped_reads = []
    # the per-chrom dir that will hold the per-sample-per-chrom bams
    chrom_dir = Path(main_by_chrom_dir, chrom)
    chrom_dir_created = False
    for sample_name, by_chrom_sample_dir in zip(samples_names, by_chrom_samples_dirs):
        bam_in_region = Path(
            by_chrom_sample_dir, f"{sample_name}{interfix}.{chrom}.bam"
        )
        # not each sample has reads mapped to this specific chrom
        if bam_in_region.exists():
            # create this dir only if necessary - maybe no sample has reads mapped to that chrom
            if not chrom_dir_created:
                chrom_dir.mkdir(exist_ok=True)
                chrom_dir_created = True

            # # copy the per-sample-per-chrom bam to the per-chrom dir
            # bam_in_region_copy = Path(chrom_dir, bam_in_region.name)
            # copy_bytes(bam_in_region, bam_in_region_copy)

            # soft-link the per-sample-per-chrom bam to the per-chrom dir
            bam_in_region_link = Path(chrom_dir, bam_in_region.name)
            bam_in_region_link.symlink_to(bam_in_region)
            # also link the .bai index file
            bam_in_region_index = Path(
                bam_in_region.parent, f"{bam_in_region.name}.bai"
            )
            bam_in_region_index_link = Path(chrom_dir, bam_in_region_index.name)
            bam_in_region_index_link.symlink_to(bam_in_region_index)

            # update stats
            _samples_names.append(sample_name)
            _mapped_reads.append(
                count_reads(
                    samtools_path,
                    bam_in_region,
                    None,
                    include_flags,
                    exclude_flags,
                    threads,
                )
            )
    df = pd.DataFrame(
        {"Chrom": chrom, "Sample": _samples_names, "MappedReads": _mapped_reads}
    )
    return df


def pacbio_main(
    *,
    genome: Path,
    in_dir: Path,
    postfix: str,
    recursive: bool,
    out_dir: Path,
    base_conda_env_dir: Path,
    pb_conda_env_name: str,
    pbmm2_path: Path,
    best_n_alignments_per_read: int,
    processes: int,
    threads: int,
    samtools_path: Path,
    multiqc_path: Path,
    include_flags: Union[int, str, None],
    exclude_flags: Union[int, str, None],
    **kwargs,
):
    preset = "CCS"

    out_dir.mkdir(exist_ok=True)

    in_files = find_files(in_dir, postfix, recursive)
    samples_names = [extract_sample_name(in_file, postfix) for in_file in in_files]
    final_postfix = postfix.split(".")[-1]
    assert final_postfix in ["fastq", "fq", "bam"]
    if final_postfix in ["fastq", "fq"]:
        raise NotImplementedError  # todo: 1 - implement, 2 - add directly to arg parser

    genome_index_file = Path(out_dir, f"{genome.name}.mmi")

    pacbio_index_genome(
        base_conda_env_dir,
        pb_conda_env_name,
        pbmm2_path,
        genome,
        genome_index_file,
        preset,
        threads,
    )

    interfix = ".aligned.sorted"
    gene_regions = None
    separate_by_chrom = False
    seperate_by_gene = False

    with Pool(processes=processes) as pool:
        pool.starmap(
            func=pacbio_align,
            iterable=[
                (
                    base_conda_env_dir,
                    pb_conda_env_name,
                    pbmm2_path,
                    preset,
                    best_n_alignments_per_read,
                    threads,
                    genome_index_file,
                    in_file,
                    Path(out_dir, f"{sample_name}{interfix}.{final_postfix}"),
                    samtools_path,
                    out_dir,
                    sample_name,
                    interfix,
                    gene_regions,
                    include_flags,
                    exclude_flags,
                    separate_by_chrom,
                    seperate_by_gene,
                )
                for in_file, sample_name in zip(in_files, samples_names)
            ],
        )

    # run MultiQC on the reports created for each aligned file
    multiqc(multiqc_path=multiqc_path, data_dir=out_dir)


def define_args() -> argparse.Namespace:
    # create parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    subparsers = parser.add_subparsers(help="sub-command help")

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
    parser.add_argument(
        "--in_dir",
        required=True,
        type=abs_path_from_str,
        help=(
            "A folder with input (bam) files, where each sample's name format is something like "
            "`$sample.ccs.bam`. See `postfix` below."
        ),
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Whether to search recursively in subdirectories of `in_dir` for input files.",
    )
    parser.add_argument("--out_dir", type=abs_path_from_str, required=True)
    parser.add_argument(
        "--samtools_path",
        type=expanded_path_from_str,
        # default=Path("~/anaconda3/envs/combinatorics/bin/samtools").expanduser(),
        default=Path("samtools"),
        help="Samtools executable.",
    )
    parser.add_argument(
        "--multiqc_path",
        type=expanded_path_from_str,
        # default=Path("~/anaconda3/envs/combinatorics/bin/multiqc").expanduser(),
        default=Path("multiqc"),
        help="MultiQC executable.",
    )
    parser.add_argument(
        "--processes",
        type=int,
        default=2,
        help="Maximal number of processes to run in parallel.",
    )
    parser.add_argument(
        "--threads", type=int, default=25, help="Threads used in each process."
    )

    # pacbio args

    pacbio_parser = subparsers.add_parser("pacbio", help="Pacbio reads help")
    pacbio_parser.set_defaults(func=pacbio_main)

    pacbio_parser.add_argument(
        "--postfix",
        default=".ccs.bam",
        help="Postfix of wanted files in `in_dir`. Should be the *full* postfix.",
    )
    pacbio_parser.add_argument(
        "--base_conda_env_dir",
        default=Path("~/anaconda3").expanduser(),
        type=expanded_path_from_str,
    )
    pacbio_parser.add_argument(
        "--pb_conda_env_name",
        default="pacbiocomb",
        help="Contains all PacBio's software packages (seperate env due to python 2.7 requirement).",
    )
    pacbio_parser.add_argument(
        "--pbmm2_path",
        type=expanded_path_from_str,
        # default=Path("~/anaconda3/envs/pacbiocomb/bin/pbmm2").expanduser(),
        default=Path("pbmm2"),
        help="Note it should match `--pb_conda_env`.",
    )
    pacbio_parser.add_argument(
        "--best_n_alignments_per_read",
        type=int,
        default=1,
        help="Output at maximum N alignments for each read, 0 means no maximum.",
    )
    pacbio_parser.add_argument(
        "--include_flags",
        default=None,
        type=int,
        help="Use only reads with this flag.",
    )
    pacbio_parser.add_argument(
        "--exclude_flags",
        default=None,
        type=int,
        help="Exclude reads with this flag.",
    )

    # pacbio isoseq args

    pacbio_isoseq_undirected_seq_parser = subparsers.add_parser(
        "whole_transcriptome_isoseq", help="Whole-transcriptome IsoSeq reads help"
    )
    pacbio_isoseq_undirected_seq_parser.set_defaults(
        func=whole_transcriptome_isoseq_main
    )

    pacbio_isoseq_undirected_seq_parser.add_argument(
        "--postfix",
        default=".fastq.gz",
        help="Postfix of wanted files in `in_dir`. Should be the *full* postfix.",
    )
    pacbio_isoseq_undirected_seq_parser.add_argument(
        "--base_conda_env_dir",
        default=Path("~/anaconda3").expanduser(),
        type=expanded_path_from_str,
    )
    pacbio_isoseq_undirected_seq_parser.add_argument(
        "--pb_conda_env_name",
        default="pacbiocomb",
        help="Contains all PacBio's software packages (seperate env due to python 2.7 requirement).",
    )
    pacbio_isoseq_undirected_seq_parser.add_argument(
        "--pbmm2_path",
        type=expanded_path_from_str,
        # default=Path("~/anaconda3/envs/pacbiocomb/bin/pbmm2").expanduser(),
        default=Path("pbmm2"),
        help="Note it should match `--pb_conda_env`.",
    )
    pacbio_isoseq_undirected_seq_parser.add_argument(
        "--best_n_alignments_per_read",
        type=int,
        default=1,
        help="Output at maximum N alignments for each read, 0 means no maximum.",
    )
    # pacbio_preprocessed_isoseq_parser.add_argument(
    #     "--gff",
    #     # required=True,
    #     type=abs_path_from_str,
    #     help=("GFF3 annotation file for by-chrom separation of BAM files."),
    # )
    pacbio_isoseq_undirected_seq_parser.add_argument(
        "--known_sites_bed_file",
        required=True,
        type=abs_path_from_str,
        help=(
            "BED file of known editing sites used for by-chrom separation of BAM files."
        ),
    )

    # illumina args

    illumina_parser = subparsers.add_parser("illumina", help="Illumina reads help")
    illumina_parser.set_defaults(func=illumina_main)

    illumina_parser.add_argument(
        "--postfix",
        default=".fastq.gz",
        help="Postfix of wanted files in `in_dir`. Should be the *full* postfix.",
    )
    illumina_parser.add_argument(
        "--mate_prefix",
        default="_",
        help="Mate prefix, e.g., `_` for `reads_1.fastq` and `reads_2.fastq`.",
    )
    illumina_parser.add_argument(
        "--decompress_cmd",
        help=(
            "Command to decompress files, e.g., `gunzip -c`. "
            "Should match the form `$decompress_cmd $in_file > $out_file`."
        ),
        default="gunzip -c",
    )
    illumina_parser.add_argument(
        "--bwa_path",
        type=expanded_path_from_str,
        # default=Path("~/anaconda3/envs/combinatorics/bin/bwa-mem2").expanduser(),
        default=Path("bwa-mem2"),
        help="BWA mem2 executable.",
    )
    illumina_parser.add_argument(
        "--require_flags",
        default="3",
        type=int,
        help="Use only reads with this flag. The default 2304 required mapped propely-paired PE reads.",
    )
    illumina_parser.add_argument(
        "--exclude_flags",
        default="2304",
        type=int,
        help="Exclude reads with this flag. The default 2304 remove secondary and supplementary (chimeric) alignments.",
    )
    illumina_parser.add_argument(
        "--separate_by_chrom",
        action="store_true",
        help="Separate aligned reads by the chromosome they were mapped to.",
    )

    return parser

    # # parse args
    # args = parser.parse_args()
    # return args


if __name__ == "__main__":

    # run
    # pacbio_main(
    #     **vars(define_and_parse_args())
    # )  # https://stackoverflow.com/a/35824590/10249633 argparse.Namespace -> Dict
    parser = define_args()
    args = parser.parse_args()
    args.func(**vars(args))

    # end
    final_words()
