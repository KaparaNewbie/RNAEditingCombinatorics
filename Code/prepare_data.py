import argparse
import subprocess
from multiprocessing import Pool
from pathlib import Path
from itertools import chain
import sys

from icecream import ic

# from Alignment.alignment_utils import pacbio_index
from General.consts import final_words
from General.os_utils import (
    extract_sample_name,
    find_files,
    delete_folder_with_files,
    group_pe_fastq_files,
    decompress,
    compress,
)
from General.argparse_utils import abs_path_from_str, expanded_path_from_str
from General.multiqc import multiqc


def pacbio_index(
    base_conda_env_dir: Path,
    pb_conda_env_name: str,
    pbindex_path: Path,
    alignment_file: Path,
):
    index_cmd = (
        f". {Path(base_conda_env_dir, 'etc/profile.d/conda.sh')} && "
        f"conda activate {pb_conda_env_name} && "
        f"{pbindex_path} "
        f"{alignment_file} "
    )
    subprocess.run(index_cmd, shell=True, executable="/bin/bash")


def basic_ccs(
    base_conda_env_dir: Path,
    pb_conda_env_name: str,
    ccs_path: Path,
    in_file: Path,
    out_file: Path,
    threads: int,
    fastqc_path: Path,
    min_passes: int,
    min_rq: float,
):
    ccs_cmd = (
        f". {Path(base_conda_env_dir, 'etc/profile.d/conda.sh')} && "
        f"conda activate {pb_conda_env_name} && "
        f"{ccs_path} "
        f"--num-threads {threads} "
        f"--min-passes {min_passes} "
        f"--min-rq {min_rq} "
        f"{in_file} "
        f"{out_file} "
    )
    subprocess.run(ccs_cmd, shell=True, executable="/bin/bash")
    fastqc(
        fastqc_path=fastqc_path,
        seq_file=out_file,
        out_dir=out_file.parent,
        threads=threads,
    )


def extract_hifi(
    base_conda_env_dir: Path,
    pb_conda_env_name: str,
    extracthifi_path: Path,
    in_file: Path,
    out_file: Path,
    threads: int,
    fastqc_path: Path,
):
    extract_hifi_cmd = (
        f". {Path(base_conda_env_dir, 'etc/profile.d/conda.sh')} && "
        f"conda activate {pb_conda_env_name} && "
        f"{extracthifi_path} "
        f"{in_file} "
        f"{out_file} "
    )
    subprocess.run(extract_hifi_cmd, shell=True, executable="/bin/bash")
    fastqc(
        fastqc_path=fastqc_path,
        seq_file=out_file,
        out_dir=out_file.parent,
        threads=threads,
    )


def remove_pacbio_duplicates(
    base_conda_env_dir: Path,
    pb_conda_env_name: str,
    pbmarkdup_path: Path,
    in_file: Path,
    out_file: Path,
    threads: int,
    fastqc_path: Path,
):
    rm_dup_cmd = (
        f". {Path(base_conda_env_dir, 'etc/profile.d/conda.sh')} && "
        f"conda activate {pb_conda_env_name} && "
        f"{pbmarkdup_path} "
        "--rmdup "  # Exclude duplicates from OUTFILE
        "--clobber "  # Overwrite OUTFILE if it exists
        f"--num-threads {threads} "
        f"{in_file} "
        f"{out_file} "
    )
    subprocess.run(rm_dup_cmd, shell=True, executable="/bin/bash")
    fastqc(
        fastqc_path=fastqc_path,
        seq_file=out_file,
        out_dir=out_file.parent,
        threads=threads,
    )


def remove_illumina_duplicates(
    prinseq_lite_path: Path,
    left_in_file: Path,
    right_in_file: Path,
    out_file: Path,
    min_qual_mean: int,
    trim_to_len: int,
    min_len: int,
    trim_left: int,
    trim_right: int,
    fastqc_path: Path,
    threads: int,
):
    # out_file = Path(out_dir, sample_name)
    rm_dup_cmd = (
        f"{prinseq_lite_path} "
        f"-fastq {left_in_file} "
        f"-fastq2 {right_in_file} "
        f"-out_bad null "
        f"-out_good {out_file} "
        f"-derep 14 "  # Type of duplicates to filter; 1 = exact duplicates, 4 = reverse complement exact duplicate
        f"-min_qual_mean {min_qual_mean} "
        f"-trim_to_len {trim_to_len} "
        f"-min_len {min_len} "
        f"-trim_left {trim_left} "
        f"-trim_right {trim_right} "
    )
    subprocess.run(rm_dup_cmd, shell=True)
    left_out_file = Path(out_file.parent, f"{out_file.name}_1.fastq")
    right_out_file = Path(out_file.parent, f"{out_file.name}_2.fastq")
    fastqc(
        fastqc_path=fastqc_path,
        seq_file=left_out_file,
        out_dir=left_out_file.parent,
        threads=threads,
    )
    fastqc(
        fastqc_path=fastqc_path,
        seq_file=right_out_file,
        out_dir=right_out_file.parent,
        threads=threads,
    )


def trimmomatic(
    trimmomatic_path: Path,
    trimmomatic_adapters_file: Path,
    threads: int,
    phred_offset: int,
    left_in_file: Path,
    right_in_file: Path,
    left_paired_out_file: Path,
    left_unpaired_out_file: Path,
    right_paired_out_file: Path,
    right_unpaired_out_file: Path,
    fastqc_path: Path,
):
    # 1 - trim
    # http://www.usadellab.org/cms/?page=trimmomatic
    cmd = (
        f"{trimmomatic_path} "
        "PE "
        f"-threads {threads} "
        f"-phred{phred_offset} "  # 33 or 64 (for older data)
        f"{left_in_file} "
        f"{right_in_file} "
        f"{left_paired_out_file} "
        f"{left_unpaired_out_file} "
        f"{right_paired_out_file} "
        f"{right_unpaired_out_file} "
        f"ILLUMINACLIP:{trimmomatic_adapters_file}:2:30:10 "  # Remove adapters (ILLUMINACLIP:TruSeq3-PE.fa:2:30:10)
        # f"LEADING:3 "  # Remove leading low quality or N bases (below quality 3) (LEADING:3)
        # f"TRAILING:3 "  # Remove trailing low quality or N bases (below quality 3) (TRAILING:3)
        # f"SLIDINGWINDOW:4:15 "  # Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 15 (SLIDINGWINDOW:4:15)
        # f"MINLEN:36 "  # Drop reads below the 36 bases long (MINLEN:36)
    )
    subprocess.run(cmd, shell=True)
    # 2 - delete unpaired files
    subprocess.run(f"rm -rf {left_unpaired_out_file}", shell=True)
    subprocess.run(f"rm -rf {right_unpaired_out_file}", shell=True)
    # 3 - fastqc for good trimmed files
    fastqc(
        fastqc_path=fastqc_path,
        seq_file=left_paired_out_file,
        out_dir=left_paired_out_file.parent,
        threads=threads,
    )
    fastqc(
        fastqc_path=fastqc_path,
        seq_file=right_paired_out_file,
        out_dir=right_paired_out_file.parent,
        threads=threads,
    )


def fastqc(
    fastqc_path: Path,
    seq_file: Path,
    out_dir: Path,
    threads: int,
):
    # fastqc [-o output dir] [--(no)extract] [-f fastq|bam|sam]
    #    [-c contaminant file] seqfile1 .. seqfileN
    out_dir.mkdir(exist_ok=True)
    temp_dir = Path(out_dir, f"{seq_file.name}.FastQCTempDir")
    temp_dir.mkdir(exist_ok=True)
    fastqc_cmd = (
        f"{fastqc_path} "
        f"--outdir {out_dir} "
        f"--threads {threads} "
        f"--dir {temp_dir} "
        f"{seq_file} "
    )
    subprocess.run(fastqc_cmd, shell=True)
    delete_folder_with_files(temp_dir)


def illumina_main(
    *,
    in_dir: Path,
    postfix: str,
    mate_prefix: str,
    recursive: bool,
    decompress_cmd: str,
    compress_cmd: str,
    out_dir: Path,
    fastqc_path: Path,
    multiqc_path: Path,
    prinseq_lite_path: Path,
    processes: int,
    threads: int,
    min_qual_mean: int,
    trim_to_len: int,
    min_len: int,
    trim_left: int,
    trim_right: int,
    trimmomatic_path: Path,
    trimmomatic_adapters_file: Path,
    phred_offset: int,
    **kwargs,
):
    out_dir.mkdir(exist_ok=True)
    raw_files = find_files(in_dir, postfix, recursive)
    # ic(raw_files)
    # samples_names = [extract_sample_name(in_file, postfix) for in_file in raw_files]
    final_postfix = postfix.split(".")[
        -1
    ]  # should be gz or lzma or something like that

    # 1 - QC for raw fastq files

    decompressed_raw_files = [
        Path(fastq_file.parent, fastq_file.stem) for fastq_file in raw_files
    ]
    # ic(decompressed_raw_files)

    with Pool(processes=processes) as pool:
        pool.starmap(
            func=decompress,
            iterable=[
                (decompress_cmd, raw_fastq_file, decompressed_raw_fastq_file)
                for raw_fastq_file, decompressed_raw_fastq_file in zip(
                    raw_files, decompressed_raw_files
                )
            ],
        )

    with Pool(processes=processes) as pool:
        pool.starmap(
            func=fastqc,
            iterable=[
                (fastqc_path, decompressed_raw_fastq_file, in_dir, threads)
                for decompressed_raw_fastq_file in decompressed_raw_files
            ],
        )

    multiqc(multiqc_path, in_dir)

    # 2 - treatments

    paired_decompressed_raw_files = group_pe_fastq_files(
        decompressed_raw_files, postfix, mate_prefix
    )
    # ic(paired_decompressed_raw_files)

    paired_samples_names = [
        extract_sample_name(
            left_file, f"{mate_prefix}1{postfix.removesuffix(f'.{final_postfix}')}"
        )
        for left_file, _ in paired_decompressed_raw_files
    ]
    # ic(paired_samples_names)

    # 2.1 - trim adapters

    trimmed_dir = Path(out_dir, "Trimmed")
    trimmed_dir.mkdir(exist_ok=True)

    paired_decompressed_good_trimmed_files = [
        (
            Path(trimmed_dir, f"{paired_sample_name}{mate_prefix}1.good.fastq"),
            Path(trimmed_dir, f"{paired_sample_name}{mate_prefix}2.good.fastq"),
        )
        for paired_sample_name in paired_samples_names
    ]
    # ic(paired_decompressed_good_trimmed_files)
    paired_decompressed_bad_trimmed_files = [
        (
            Path(trimmed_dir, f"{paired_sample_name}{mate_prefix}1.bad.fastq"),
            Path(trimmed_dir, f"{paired_sample_name}{mate_prefix}2.bad.fastq"),
        )
        for paired_sample_name in paired_samples_names
    ]
    # ic(paired_decompressed_bad_trimmed_files)

    with Pool(processes=processes) as pool:
        pool.starmap(
            func=trimmomatic,
            iterable=[
                (
                    trimmomatic_path,
                    trimmomatic_adapters_file,
                    threads,
                    phred_offset,
                    left_in_file,
                    right_in_file,
                    left_paired_out_file,
                    left_unpaired_out_file,
                    right_paired_out_file,
                    right_unpaired_out_file,
                    fastqc_path,
                )
                for (
                    (left_in_file, right_in_file),
                    (left_paired_out_file, right_paired_out_file),
                    (left_unpaired_out_file, right_unpaired_out_file),
                ) in zip(
                    paired_decompressed_raw_files,
                    paired_decompressed_good_trimmed_files,
                    paired_decompressed_bad_trimmed_files,
                )
            ],
        )

    multiqc(multiqc_path, trimmed_dir)

    # 2.2 - remove duplicates

    trimmed_wo_dup_dir = Path(out_dir, "TrimmedWoDup")
    trimmed_wo_dup_dir.mkdir(exist_ok=True)

    prinseq_out_files = [
        Path(trimmed_wo_dup_dir, sample_name) for sample_name in paired_samples_names
    ]

    with Pool(processes=processes) as pool:
        pool.starmap(
            func=remove_illumina_duplicates,
            iterable=[
                (
                    prinseq_lite_path,
                    left_in_file,
                    right_in_file,
                    prinseq_out_file,
                    min_qual_mean,
                    trim_to_len,
                    min_len,
                    trim_left,
                    trim_right,
                    fastqc_path,
                    threads,
                )
                for (left_in_file, right_in_file), prinseq_out_file in zip(
                    paired_decompressed_good_trimmed_files, prinseq_out_files
                )
            ],
        )

    multiqc(multiqc_path, trimmed_wo_dup_dir)

    # 3 - finishes

    # 3.1 - compress needed decompressed files

    paired_compressed_trimmed_files = [
        (
            Path(trimmed_dir, f"{decompressed_left_file.name}.{final_postfix}"),
            Path(trimmed_dir, f"{decompressed_right_file.name}.{final_postfix}"),
        )
        for decompressed_left_file, decompressed_right_file in paired_decompressed_good_trimmed_files
    ]
    with Pool(processes=processes) as pool:
        pool.starmap(
            func=compress,
            iterable=[
                (compress_cmd, decompressed_in_file, compressed_out_file)
                for decompressed_paired_files, compressed_paired_files in zip(
                    paired_decompressed_good_trimmed_files,
                    paired_compressed_trimmed_files,
                )
                for decompressed_in_file, compressed_out_file in zip(
                    decompressed_paired_files, compressed_paired_files
                )
            ],
        )

    decompressed_trimmed_wo_dup_files = [
        Path(trimmed_wo_dup_dir, f"{sample_name}{mate_prefix}{mate}.fastq")
        for sample_name in paired_samples_names
        for mate in [1, 2]
    ]
    compressed_trimmed_wo_dup_files = [
        Path(trimmed_wo_dup_dir, f"{decompressed_file.name}.{final_postfix}")
        for decompressed_file in decompressed_trimmed_wo_dup_files
    ]
    with Pool(processes=processes) as pool:
        pool.starmap(
            func=compress,
            iterable=[
                (compress_cmd, decompressed_in_file, compressed_out_file)
                for decompressed_in_file, compressed_out_file in zip(
                    decompressed_trimmed_wo_dup_files, compressed_trimmed_wo_dup_files
                )
            ],
        )

    # 3.2 - delete decompressed & other unneeded files

    files_to_delete = (
        decompressed_raw_files
        + list(chain.from_iterable(paired_decompressed_good_trimmed_files))
        + list(chain.from_iterable(paired_decompressed_bad_trimmed_files))
        + decompressed_trimmed_wo_dup_files
        + [f for f in trimmed_wo_dup_dir.iterdir() if "singletons" in f.name]
    )
    for file_to_delete in files_to_delete:
        subprocess.run(f"rm -rf {file_to_delete}", shell=True)


def pacbio_main(
    *,
    in_dir: Path,
    postfix: str,
    recursive: bool,
    out_dir: Path,
    base_conda_env_dir: Path,
    pb_conda_env_name: str,
    ccs_path: Path,
    extracthifi_path: Path,
    pbindex_path: Path,
    pbmarkdup_path: Path,
    fastqc_path: Path,
    multiqc_path: Path,
    processes: int,
    threads: int,
    min_passes: int,
    min_rq: float,
    **kwargs,
):
    out_dir.mkdir(exist_ok=True)
    in_files = find_files(in_dir, postfix, recursive)
    samples_names = [extract_sample_name(in_file, postfix) for in_file in in_files]
    final_postfix = postfix.split(".")[-1]
    assert final_postfix in ["fastq", "bam"]
    if final_postfix == "fastq":
        raise NotImplementedError  # todo: 1 - implement, 2 - add directly to arg parser

    # 1 - basic ccs files

    basic_ccs_dir = Path(out_dir, "BasicCCS")
    basic_ccs_dir.mkdir(exist_ok=True)

    basic_ccs_files = {
        sample_name: Path(basic_ccs_dir, f"{sample_name}.ccs.{final_postfix}")
        for sample_name in samples_names
    }

    with Pool(processes=processes) as pool:
        pool.starmap(
            func=basic_ccs,
            iterable=[
                (
                    base_conda_env_dir,
                    pb_conda_env_name,
                    ccs_path,
                    in_file,
                    basic_ccs_files[sample_name],
                    threads,
                    fastqc_path,
                    min_passes,
                    min_rq,
                )
                for in_file, sample_name in zip(in_files, samples_names)
            ],
        )

    multiqc(multiqc_path, basic_ccs_dir)

    # 2 - extract hifi reads

    # (a) extract

    hifi_dir = Path(out_dir, "HiFi")
    hifi_dir.mkdir(exist_ok=True)

    hifi_files = {
        sample_name: Path(hifi_dir, f"{sample_name}.ccs.{final_postfix}")
        for sample_name in samples_names
    }

    with Pool(processes=processes) as pool:
        pool.starmap(
            func=extract_hifi,
            iterable=[
                (
                    base_conda_env_dir,
                    pb_conda_env_name,
                    extracthifi_path,
                    basic_ccs_files[sample_name],
                    hifi_files[sample_name],
                    threads,
                    fastqc_path,
                )
                for sample_name in samples_names
            ],
        )

    multiqc(multiqc_path, hifi_dir)

    # (b) index

    with Pool(processes=processes) as pool:
        pool.starmap(
            func=pacbio_index,
            iterable=[
                (
                    base_conda_env_dir,
                    pb_conda_env_name,
                    pbindex_path,
                    hifi_files[sample_name],
                )
                for sample_name in samples_names
            ],
        )

    # 3 - remove duplicates

    wo_dup_dir = Path(out_dir, "WoDup")
    wo_dup_dir.mkdir(exist_ok=True)

    wo_dup_files = {
        sample_name: Path(wo_dup_dir, f"{sample_name}.ccs.{final_postfix}")
        for sample_name in samples_names
    }

    with Pool(processes=processes) as pool:
        pool.starmap(
            func=remove_pacbio_duplicates,
            iterable=[
                (
                    base_conda_env_dir,
                    pb_conda_env_name,
                    pbmarkdup_path,
                    hifi_files[sample_name],
                    wo_dup_files[sample_name],
                    threads,
                    fastqc_path,
                )
                for sample_name in samples_names
            ],
        )

    multiqc(multiqc_path, wo_dup_dir)


# def define_and_parse_args() -> argparse.Namespace:
def define_args() -> argparse.ArgumentParser:

    # create common & sub parsers

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    subparsers = parser.add_subparsers(help="sub-command help")

    # common args

    parser.add_argument(
        "--in_dir",
        required=True,
        type=abs_path_from_str,
        help=(
            "A folder with input (bam) files, where each sample's name format is something like "
            "`$sample.subreads.bam`. See `postfix` below."
        ),
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Whether to search recursively in subdirectories of `in_dir` for input files.",
    )
    parser.add_argument("--out_dir", type=abs_path_from_str, required=True)
    parser.add_argument(
        "--fastqc_path",
        type=expanded_path_from_str,
        default=Path("~/anaconda3/envs/combinatorics/bin/fastqc").expanduser(),
        help="FastQC executable.",
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
        default=6,
        help="Maximal number of processes to run in parallel.",
    )
    parser.add_argument(
        "--threads", type=int, default=10, help="Threads used in each process."
    )

    # pacbio args

    pacbio_parser = subparsers.add_parser("pacbio", help="Pacbio reads help")
    pacbio_parser.set_defaults(func=pacbio_main)

    pacbio_parser.add_argument(
        "--postfix",
        default=".subreads.bam",
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
        "--ccs_path",
        type=expanded_path_from_str,
        default=Path("~/anaconda3/envs/pacbiocomb/bin/ccs").expanduser(),
        help="Note it should match `--pb_conda_env`.",
    )
    pacbio_parser.add_argument(
        "--extracthifi_path",
        type=expanded_path_from_str,
        default=Path("~/anaconda3/envs/pacbiocomb/bin/extracthifi").expanduser(),
        help="Note it should match `--pb_conda_env`.",
    )
    pacbio_parser.add_argument(
        "--pbindex_path",
        type=expanded_path_from_str,
        default=Path("~/anaconda3/envs/pacbiocomb/bin/pbindex").expanduser(),
        help="Note it should match `--pb_conda_env`.",
    )
    pacbio_parser.add_argument(
        "--pbmarkdup_path",
        type=expanded_path_from_str,
        default=Path("~/anaconda3/envs/pacbiocomb/bin/pbmarkdup").expanduser(),
        help="Note it should match `--pb_conda_env`.",
    )
    pacbio_parser.add_argument(
        "--min_passes",
        default=3,
        type=int,
        help="Minimum number of full-length subreads required to generate CCS for a ZMW.",
    )
    pacbio_parser.add_argument(
        "--min_rq",
        default=0.99,
        type=float,
        help="Minimum predicted accuracy in [0, 1] of PacBio's CCS.",
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
        "--compress_cmd",
        help=(
            "Command to compress files, e.g., `gzip`. "
            "Should match the form `$compress_cmd $in_file > $out_file`."
        ),
        default="gzip -c",
    )
    illumina_parser.add_argument(
        "--prinseq_lite_path",
        type=expanded_path_from_str,
        default=Path("~/anaconda3/envs/combinatorics/bin/prinseq-lite.pl").expanduser(),
        help="Prinseq-lite executable.",
    )
    illumina_parser.add_argument(
        "--min_qual_mean", default=30, type=int, help="Minimum mean quality for a read."
    )
    illumina_parser.add_argument(
        "--trim_to_len",
        default=sys.maxsize,
        type=int,
        help="Trim all sequence from the 3'-end to result in sequence with this length.",
    )
    illumina_parser.add_argument(
        "--min_len", default=0, type=int, help="Filter sequence shorter than min_len."
    )
    illumina_parser.add_argument(
        "--trim_left",
        default=0,
        type=int,
        help="Trim sequence at the 5'-end by trim_left positions.",
    )
    illumina_parser.add_argument(
        "--trim_right",
        default=0,
        type=int,
        help="Trim sequence at the 3'-end by trim_left positions.",
    )

    illumina_parser.add_argument(
        "--trimmomatic_path",
        type=expanded_path_from_str,
        default=Path("~/anaconda3/envs/combinatorics/bin/trimmomatic").expanduser(),
        help="Trimmomatic executable.",
    )
    illumina_parser.add_argument(
        "--trimmomatic_adapters_file",
        type=expanded_path_from_str,
        default=Path(
            "~/anaconda3/envs/combinatorics/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa"
        ).expanduser(),
        help="A fasta file with adapters to trim.",
    )
    illumina_parser.add_argument(
        "--phred_offset",
        default=33,
        type=int,
        help="Phred offset. Either 33 (usually) or 64 (old Illumina machines).",
    )

    return parser

    # parse args

    # args = parser.parse_args()
    # return args


if __name__ == "__main__":

    # run
    # main(
    #     **vars(defineargs())
    # )  # https://stackoverflow.com/a/35824590/10249633 argparse.Namespace -> Dict
    parser = define_args()
    args = parser.parse_args()
    args.func(**vars(args))

    # end
    final_words()
