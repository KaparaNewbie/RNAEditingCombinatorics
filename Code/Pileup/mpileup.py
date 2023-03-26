from pathlib import Path
from typing import Union
import subprocess

from icecream import ic

from Alignment.alignment_utils import count_reads


def define_max_depth(
    samtools_path: Path,
    in_bam: Path,
    region: str,
    include_flags: Union[int, str, None],
    exclude_flags: Union[int, str, None],
    threads: int,
    # max_depth_padding: int = 10_000,
    max_depth_padding: float = 3.0,
):
    legal_reads_in_region = count_reads(
        samtools_path, in_bam, region, include_flags, exclude_flags, threads
    )
    legal_reads_in_region = count_reads(
        samtools_path, in_bam, None, None, None, threads
    )
    # return legal_reads_in_region + max_depth_padding
    return int(legal_reads_in_region * max_depth_padding)


def mpileup(
    samtools_path: Path,
    genome: Path,
    region: str,
    include_flags: Union[int, str, None],
    exclude_flags: Union[int, str, None],
    min_bq: Union[int, None],
    in_bam: Path,
    out_file: Path,
    threads: int,
    override_existing_pileup_files: bool, 
):
    # only run mpileup if `out_file` doesn't exist or if it does, 
    # but the user whishes to recreate it
    if out_file.exists() and not override_existing_pileup_files:
        return
    
    max_depth = define_max_depth(
        samtools_path, in_bam, region, include_flags, exclude_flags, threads
    )
    mpileup_cmd = (
        f"{samtools_path} "
        "mpileup "
        f"--fasta-ref {genome} "
        f"--region {region} "
        "--no-output-ins --no-output-ins "
        "--no-output-del --no-output-del "
        "--no-output-ends "
        "--output-QNAME "
        f"--max-depth {max_depth} "
    )
    if include_flags:
        mpileup_cmd += f"--incl-flags {include_flags} "
    if exclude_flags:
        mpileup_cmd += f"--excl-flags {exclude_flags} "  # 2304 (default) = remove secondary and supplementary (chimeric) alignments
    if min_bq:
        mpileup_cmd += f"--min-BQ {min_bq} "
    # mpileup_cmd += f"{in_bam} --output {out_file} "
    mpileup_cmd += f"--output {out_file} {in_bam}"

    ic(mpileup_cmd)

    subprocess.run(mpileup_cmd, shell=True)
