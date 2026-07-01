import argparse
import subprocess
from multiprocessing import Pool
from pathlib import Path
import sys
# from itertools import chain
# import sys
# import os
# import shutil

from icecream import ic

code_dir = "/private7/projects/Combinatorics/Code"
sys.path.append(str(Path(code_dir).absolute()))

# from Alignment.alignment_utils import pacbio_index
from General.consts import final_words
from General.os_utils import (
    extract_sample_name,
    find_files,
    delete_folder_with_files,
    group_pe_fastq_files,
)
from General.argparse_utils import abs_path_from_str, expanded_path_from_str


def create_salmon_index(trinity_file: Path, index_dir: Path, threads: int):
    index_cmd = f"salmon index -t {trinity_file} -i {index_dir} --threads {threads}"
    subprocess.run(index_cmd, shell=True)


def run_salmon(
    index_dir: Path, fastq_1: Path, fastq_2: Path, threads: int, out_dir: Path
):
    quant_cmd = f"salmon quant -i {index_dir} -l A -1 {fastq_1} -2 {fastq_2} -p {threads} -o {out_dir}"
    subprocess.run(quant_cmd, shell=True)
    
    
def define_args() -> argparse.ArgumentParser:
    # create common & sub parsers

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )


    parser.add_argument(
        "--transcriptome_file",
        help="Path to trinity transcriptome file",
        type=abs_path_from_str,
        required=True,
    )
    parser.add_argument(
        "--processes",
        type=int,
        default=6,
        help="Maximal number of processes to run in parallel.",
    )
    parser.add_argument(
        "--threads", type=int, default=10, help="Threads used in each process."
    )
    parser.add_argument("--in_dir", type=abs_path_from_str, required=True)
    parser.add_argument(
        "--postfix",
        default=".fastq.gz",
        help="Postfix of wanted files in `in_dir`. Should be the *full* postfix.",
    )
    parser.add_argument(
        "--mate_prefix",
        default="_",
        help="Mate prefix, e.g., `_` for `$sample_1.fastq.gz` and `$sample_2.fastq.gz`.",
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Whether to search recursively in subdirectories of `in_dir` for input files.",
    )
    parser.add_argument("--out_dir", type=abs_path_from_str, required=True)

    return parser


if __name__ == "__main__":
    parser = define_args()
    args = parser.parse_args()
    transcriptome_fasta_file = args.transcriptome_file
    processes = args.processes
    threads = args.threads
    in_dir = args.in_dir
    postfix = args.postfix
    mate_prefix = args.mate_prefix
    recursive = args.recursive
    out_dir = args.out_dir
    
else:
    transcriptome_fasta_file = Path("/private7/projects/Combinatorics/D.pealeii/Annotations/Jan2025/orfs_squ.fa")
    in_dir = Path("/private6/projects/Combinatorics/D.pealeii/Data/BulkRNA/PRJNA1216919/TrimmedWoDup")
    out_dir = Path("/private6/projects/Combinatorics/D.pealeii/Salmon/PRJNA1216919")
    postfix = ".fastq.gz"
    mate_prefix = "_"
    recursive = False
    processes = 6
    threads = 10
    

max_threads = processes * threads

out_dir.mkdir(exist_ok=True, parents=True)

in_files = find_files(in_dir, postfix, recursive)

paired_in_files = group_pe_fastq_files(
    in_files, postfix, mate_prefix
)

sample_names = [
    extract_sample_name(fastq_1, f"_1{postfix}") for fastq_1, _ in paired_in_files
]
samples_out_dirs = [
    Path(out_dir, sample_name)
    for sample_name in sample_names
]

salmon_index_dir = Path(out_dir, "SalmonIndex")

create_salmon_index(transcriptome_fasta_file, salmon_index_dir, max_threads)

index_cmd = f"salmon index -t {transcriptome_fasta_file} -i {salmon_index_dir} --threads {max_threads}"
subprocess.run(index_cmd, shell=True)


with Pool(processes) as pool:
    pool.starmap(
        func=run_salmon,
        iterable=[
            (
                salmon_index_dir,
                fastq_1,
                fastq_2,
                threads,
                samples_out_dir,
            )
            for (fastq_1, fastq_2), samples_out_dir in zip(paired_in_files, samples_out_dirs)
        ],
    )

delete_folder_with_files(salmon_index_dir)

final_words()