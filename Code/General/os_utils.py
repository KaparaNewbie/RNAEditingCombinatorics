from pathlib import Path
import subprocess

from General.type_hints import StrOrPath


def find_files(in_dir: Path, postfix: str, recursive: bool) -> list[Path]:
    if recursive:
        files = list(in_dir.glob(f"**/*{postfix}"))
    else:
        files = list(in_dir.glob(f"*{postfix}"))
    return files


def extract_sample_name(sample: Path, postfix: str) -> str:
    name = sample.name.removesuffix(postfix)
    return name


def delete_folder_with_files(folder: StrOrPath):
    """
    Delete folder and everything in it.
    """
    subprocess.run(f"rm -rf {folder}", shell=True)


def group_pe_fastq_files(fastq_files: list[Path], postfix: str, mate_prefix: str):
    """Sort files in fastq_files into tuples of paired-ended files.

    Args:
        fastq_files (list[Path]): paths of all fastq files
        postfix (str): e.g., postfix of `sample_1.fastq.gz` is `.fastq.gz`
        mate_prefix (str): e.g., mate_prefix of `sample_1.fastq.gz` is `_`

    Returns:
        list[tuple[Path, Path]]: list of tuples of paired-ended fastq files
    """
    sorted_fastq_files = sorted(
        fastq_files,
        key=lambda fastq_file: fastq_file.name.rstrip(postfix).split(mate_prefix),
    )
    paired_fastq_files = [
        (sorted_fastq_files[x], sorted_fastq_files[x + 1])
        for x in range(0, len(sorted_fastq_files), 2)
    ]
    return paired_fastq_files


def decompress(decompress_cmd: str, in_file: Path, out_file: Path):
    cmd = f"{decompress_cmd} {in_file} > {out_file}"
    subprocess.run(cmd, shell=True)


def compress(compress_cmd: str, in_file: Path, out_file: Path):
    cmd = f"{compress_cmd} {in_file} > {out_file}"
    subprocess.run(cmd, shell=True)


def copy_text(in_file: Path, out_file: Path):
    out_file.write_text(in_file.read_text())


def copy_bytes(in_file: Path, out_file: Path):
    out_file.write_bytes(in_file.read_bytes())
