from pathlib import Path
from typing import Union
import subprocess

import pysam


def filter_bam_by_read_quality(
    samtools_path: Path, in_bam: Path, min_rq: float, threads: int, out_dir: Path
) -> Path:
    """Filter a BAM file by read quality.

    Args:
        samtools_path (Path): Path to the samtools executable.
        in_bam (Path): Path to the input BAM file.
        min_rq (float): Minimum read quality to keep.
        threads (int): Number of threads to use.
        out_dir (Path): Write the filtered BAM file to this directory.

    Returns:
        Path: Path of the filtered BAM file.
    """
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
    filtered_bam = Path(out_dir, f"{in_bam.stem}.MinRQ{str(min_rq)[2:]}.bam")
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
    cmd = f"{samtools_path} view -c --threads {threads} "
    if include_flags:
        cmd += f"--require-flags {include_flags} "
    if exclude_flags:
        cmd += f"--excl-flags {exclude_flags} "
    cmd += f"{in_bam} " 
    if region:
        cmd += f"{region}"
    reads = int(subprocess.run(cmd, shell=True, capture_output=True).stdout.decode())
    return reads


def count_reads_in_unaligned_bam(samtools_path: Path, in_bam: Path, threads: int):
    cmd = f"{samtools_path} view -c --threads {threads} {in_bam}"
    reads = int(subprocess.run(cmd, shell=True, capture_output=True).stdout.decode())
    return reads
