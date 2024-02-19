from pathlib import Path
from typing import Union
import subprocess

# import os

import pysam
from icecream import ic


def sample_bam(
    samtools_path: Path,
    in_bam: Path,
    out_bam: Path,
    num_reads: float,
    seed: int,
    threads: int,
    override_existing_out_file: bool = True
):
    """Sample a BAM file.

    Args:
        samtools_path (Path): Path to the samtools executable.
        in_bam (Path): Path to the input BAM file.
        out_bam (Path): Path to the output BAM file.
        num_reads (int): Number of reads to sample.
        seed (int): Seed for the random number generator.
        threads (int): Number of threads to use.
        override_existing_out_file (bool, optional): If True, override the existing output file. Defaults to True.
    """
    # check if the output file already exists
    if out_bam.exists() and not override_existing_out_file:
        return
    # calculate fraction of reads to sample by dividing the number of total number in the number of reads to sample
    total_reads = count_reads(samtools_path, in_bam, None, None, None, threads)
    fraction = num_reads / total_reads
    # run samtools view with the fraction of reads to sample
    cmd = (
        f"{samtools_path} "
        "view "
        "--with-header "
        f"--subsample {fraction} "
        f"--subsample-seed {seed} "
        f"--threads {threads} "
        "--bam "
        f"--output {out_bam} "
        f"{in_bam} "
    )
    subprocess.run(cmd, shell=True)
    # index the new bam file
    index_cmd = f"{samtools_path} index -@ {threads} {out_bam} "
    print(index_cmd)
    subprocess.run(index_cmd, shell=True)


def filter_bam_by_read_quality(
    samtools_path: Path, in_bam: Path, min_rq: float, threads: int, out_dir: Path, override_existing_out_file: bool = True
) -> Path:
    """Filter a BAM file by read quality.

    Args:
        samtools_path (Path): Path to the samtools executable.
        in_bam (Path): Path to the input BAM file.
        min_rq (float): Minimum read quality to keep.
        threads (int): Number of threads to use.
        out_dir (Path): Write the filtered BAM file to this directory.
        override_existing_out_file (bool, optional): If True, override the existing output file. Defaults to True.

    Returns:
        Path: Path of the filtered BAM file.
    """
    # 0 - check if the output file already exists
    filtered_bam = Path(out_dir, f"{in_bam.stem}.MinRQ{str(min_rq)[2:]}.bam")
    if filtered_bam.exists() and not override_existing_out_file:
        return filtered_bam
    # 1 - find reads with quality equal to or greater than the threshold min_rq
    with pysam.AlignmentFile(in_bam, "rb") as samfile:
        high_quality_reads_names = [
            read.query_name for read in samfile if read.get_tag("rq") >= min_rq
        ]
    # 2 - save their names to a file
    high_quality_reads_file = Path(
        out_dir, f"{in_bam.stem}.ReadsNames.MinRQ{str(min_rq)[2:]}.txt"
    )
    with high_quality_reads_file.open("w") as f:
        f.write("\n".join(high_quality_reads_names))
    # 3 - create a new bam by using only those high quality reads
    # filtered_bam = Path(out_dir, f"{in_bam.stem}.MinRQ{str(min_rq)[2:]}.bam")
    filter_cmd = (
        f"{samtools_path} "
        "view "
        f"-N {high_quality_reads_file} "
        "--with-header "
        "--bam "
        f"--output {filtered_bam} "
        f"--threads {threads} "
        f"{in_bam} "
    )
    subprocess.run(filter_cmd, shell=True)
    # 4 - index the new bam file
    index_cmd = f"{samtools_path} index -@ {threads} {filtered_bam} "
    print(index_cmd)
    subprocess.run(index_cmd, shell=True)
    # 5 - return path to the new bam file
    return filtered_bam


def count_reads(
    samtools_path: Path,
    in_bam: Path,
    region: Union[str, None],
    include_flags: Union[int, str, None],
    exclude_flags: Union[int, str, None],
    threads: int,
):
    """Count number of reads in `in_bam`. Optionally, count only reads satisfying certain conditions.

    Args:
        samtools_path (Path): path to executable
        in_bam (Path): input BAM file
        region (Union[str, None]): if not None, count only reads mapped to `region`
        include_flags (Union[int, str, None]): if not None, count only reads with these flags
        exclude_flags (Union[int, str, None]): if not None, count only reads *without* these flags
        threads (int): number of `threads` to use

    Returns:
        int: number of reads in `in_bam` satisfying certain conditions (or all reads if no conditions were specified)
    """
    cmd = f"{samtools_path} view -c --threads {threads} "
    if include_flags:
        cmd += f"--require-flags {include_flags} "
    if exclude_flags:
        cmd += f"--excl-flags {exclude_flags} "
    cmd += f"{in_bam} "
    if region:
        cmd += f"{region}"

    # reads = int(subprocess.run(cmd, shell=True, capture_output=True).stdout.decode())

    stdout = subprocess.run(cmd, shell=True, capture_output=True).stdout.decode()
    if stdout == "":
        reads = 0
    else:
        reads = int(stdout)

    return reads


def count_unique_reads(bam: Path, threads: int = 1):
    with pysam.AlignmentFile(bam, "rb", threads=threads) as samfile:
        unique_reads_names = {read.query_name for read in samfile}
    return len(unique_reads_names)


def count_unique_filtered_aligned_reads(
    bam: Path,
    region: Union[str, None],
    include_flags: Union[int, str, None],
    exclude_flags: Union[int, str, None],
    threads: int = 1,
):
    """Count number of unique reads in `bam` satisfying certain conditions.

    Args:
        bam (Path): input BAM file.
        region (Union[str, None]): if not None, count only reads mapped to `region`.
        include_flags (Union[int, str, None]): if not None, count only reads with these flags.
        exclude_flags (Union[int, str, None]): if not None, count only reads *without* these flags.
        threads (int, optional): Number of `threads` to use. Defaults to 1.

    Returns:
        int: number of unique reads in `bam` satisfying certain conditions (or all reads if no conditions were specified).
    """
    cmd = f"samtools view --threads {threads} "
    if include_flags:
        cmd += f"--require-flags {include_flags} "
    if exclude_flags:
        cmd += f"--excl-flags {exclude_flags} "
    cmd += f"{bam} "
    if region:
        cmd += f"{region}"
    cmd += " | cut -f 1 | sort | uniq | wc -l"
    # print(cmd)
    reads = int(subprocess.run(cmd, shell=True, capture_output=True).stdout.decode())
    return reads


def count_reads_in_unaligned_bam(samtools_path: Path, in_bam: Path, threads: int):
    cmd = f"{samtools_path} view -c --threads {threads} {in_bam}"
    reads = int(subprocess.run(cmd, shell=True, capture_output=True).stdout.decode())
    return reads


def count_reads_in_fastq(fatsq_path: Path, cat_cmd: str = "cat"):
    cmd = f"echo $({cat_cmd} {fatsq_path} | wc -l) / 4 | bc"
    reads = int(subprocess.run(cmd, shell=True, capture_output=True).stdout.decode())
    return reads
