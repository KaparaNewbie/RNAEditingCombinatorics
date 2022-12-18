import argparse
from concurrent.futures import thread
import subprocess
from multiprocessing import Pool
from pathlib import Path
from typing import Sequence

import pysam
import pandas as pd

from General.multiqc import multiqc
from General.consts import final_words
from General.os_utils import (
    extract_sample_name,
    find_files,
    group_pe_fastq_files,
    decompress,
    delete_folder_with_files,
    copy,
)
from General.argparse_utils import abs_path_from_str, expanded_path_from_str
from Alignment.alignment_utils import count_reads

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
        separte_to_chroms(
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


def separte_to_chroms(
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
        index_cmd = f"{samtools_path} index {bam_in_region}"
        subprocess.run(index_cmd, shell=True, cwd=by_chrom_dir)


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
    separate_by_chrom: bool = False,
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
        separte_to_chroms(
            samtools_path,
            out_dir,
            threads,
            sample_name,
            bam=out_file,
            interfix=interfix,
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


def pacbio_preprocessed_isoseq_main(
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
    **kwargs,
):
    preset = "ISOSEQ"

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

    # align

    interfix = ".aligned.sorted"

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
                    True,
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
    _chroms = []  # for stats
    _samples = []  # for stats
    _reads = []  # for stats
    main_by_chrom_dir = Path(out_dir, "ByChrom")
    for chrom in chroms:
        chrom_dir = Path(main_by_chrom_dir, chrom)
        chrom_dir_created = False
        for sample_name, by_chrom_sample_dir in zip(
            samples_names, by_chrom_samples_dirs
        ):
            bam_in_region = Path(
                by_chrom_sample_dir, f"{sample_name}{interfix}.{chrom}.bam"
            )
            if bam_in_region.exists():
                if not chrom_dir_created:
                    chrom_dir.mkdir(exist_ok=True)
                    chrom_dir_created = True
                bam_in_region_copy = Path(chrom_dir, bam_in_region.name)
                copy(bam_in_region, bam_in_region_copy)
                # update stats
                _chroms.append(chrom)
                _samples.append(sample_name)
                _reads.append(
                    count_reads(samtools_path, bam_in_region, None, None, None, threads)
                )
    # write by=chrom stats
    df = pd.DataFrame({"Chrom": _chroms, "Sample": _samples, "MappedReads": _reads})
    df.to_csv(Path(out_dir, "ByChromSummary.tsv"), sep="\t", index=False)

    # run MultiQC on the reports created for each whole aligned file
    multiqc(multiqc_path=multiqc_path, data_dir=out_dir)

    # delete the genome's index (it's quite heavy)
    genome_index_file.unlink(missing_ok=True)


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
                    False,
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
        default=Path("~/anaconda3/envs/combinatorics/bin/samtools").expanduser(),
        help="Samtools executable.",
    )
    parser.add_argument(
        "--multiqc_path",
        type=expanded_path_from_str,
        default=Path("~/anaconda3/envs/combinatorics/bin/multiqc").expanduser(),
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
        default=Path("~/anaconda3/envs/pacbiocomb/bin/pbmm2").expanduser(),
        help="Note it should match `--pb_conda_env`.",
    )
    pacbio_parser.add_argument(
        "--best_n_alignments_per_read",
        type=int,
        default=1,
        help="Output at maximum N alignments for each read, 0 means no maximum.",
    )

    # pacbio isoseq args

    pacbio_preprocessed_isoseq_parser = subparsers.add_parser(
        "pacbio_preprocessed_isoseq", help="Preprocessed IsoSeq reads help"
    )
    pacbio_preprocessed_isoseq_parser.set_defaults(func=pacbio_preprocessed_isoseq_main)

    pacbio_preprocessed_isoseq_parser.add_argument(
        "--postfix",
        default=".fastq.gz",
        help="Postfix of wanted files in `in_dir`. Should be the *full* postfix.",
    )
    pacbio_preprocessed_isoseq_parser.add_argument(
        "--base_conda_env_dir",
        default=Path("~/anaconda3").expanduser(),
        type=expanded_path_from_str,
    )
    pacbio_preprocessed_isoseq_parser.add_argument(
        "--pb_conda_env_name",
        default="pacbiocomb",
        help="Contains all PacBio's software packages (seperate env due to python 2.7 requirement).",
    )
    pacbio_preprocessed_isoseq_parser.add_argument(
        "--pbmm2_path",
        type=expanded_path_from_str,
        default=Path("~/anaconda3/envs/pacbiocomb/bin/pbmm2").expanduser(),
        help="Note it should match `--pb_conda_env`.",
    )
    pacbio_preprocessed_isoseq_parser.add_argument(
        "--best_n_alignments_per_read",
        type=int,
        default=1,
        help="Output at maximum N alignments for each read, 0 means no maximum.",
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
        default=Path("~/anaconda3/envs/combinatorics/bin/bwa-mem2").expanduser(),
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
