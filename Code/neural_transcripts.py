from pathlib import Path
import argparse
from typing import Union
import subprocess
from multiprocessing import Pool

import pandas as pd

from General.consts import final_words
from General.argparse_utils import abs_path_from_str

# from NeuralTranscripts.prepare_salmon_gtf_input_from_transcriptome import (
#     main as prepare_gtf_for_salmon,
# )
from NeuralTranscripts.orthofinder import main as orthofinder
from General.os_utils import (
    extract_sample_name,
    find_files,
    group_pe_fastq_files,
)


def create_salmon_index(trinity_file: Path, index_dir: Path, threads: int):
    index_cmd = f"salmon index -t {trinity_file} -i {index_dir} --threads {threads}"
    subprocess.run(index_cmd, shell=True)


def run_salmon(
    index_dir: Path, fastq_1: Path, fastq_2: Path, threads: int, out_dir: Path
):
    quant_cmd = f"salmon quant -i {index_dir} -l A -1 {fastq_1} -2 {fastq_2} -p {threads} -o {out_dir}"
    subprocess.run(quant_cmd, shell=True)


def determine_obim_neural_transcripts(
    sra_run_table_file: Path, salmon_quant_dirs: Path, out_table: Path
):
    sra_df = pd.read_csv(sra_run_table_file, usecols=["Run", "Tissue"]).rename(
        columns={"Run": "Sample"}
    )
    neural_tissues = [
        "optic_lobe",
        "axial nerve cord",
        "subesophageal brain",
        "supraesophageal brain",
    ]
    sra_df["IsNeuralTissue"] = sra_df["Tissue"].isin(neural_tissues)
    # sra_df

    # out_dir = Path("/private7/projects/Combinatorics/O.bim/Expression")
    # salmon_quant_dirs = list(out_dir.glob("SRR*"))
    quant_files = [
        Path(salmon_quant_dir, "quant.sf") for salmon_quant_dir in salmon_quant_dirs
    ]
    quant_dfs = [
        pd.read_csv(quant_file, sep="\t", usecols=["Name", "TPM"]).rename(
            columns={"Name": "Chrom"}
        )
        for quant_file in quant_files
    ]
    for salmon_quant_dir, quant_df in zip(salmon_quant_dirs, quant_dfs):
        sample_name = salmon_quant_dir.name
        quant_df.insert(0, "Sample", sample_name)
    # quant_dfs[0]
    quant_df = pd.concat(quant_dfs, ignore_index=True)
    quant_df = quant_df.merge(sra_df, on="Sample", how="left")

    # quant_df

    mean_tpm_per_chrom_df = (
        quant_df.groupby(["IsNeuralTissue", "Chrom"])["TPM"].mean().reset_index()
    )
    mean_tpm_per_chrom_df = (
        mean_tpm_per_chrom_df.pivot(
            index="Chrom", columns="IsNeuralTissue", values="TPM"
        )
        .rename(columns={True: "MeanNeuralTPM", False: "MeanNonNeuralTPM"})
        .reset_index()
        .rename_axis(None, axis=1)
        .rename_axis(None)
    )
    mean_tpm_per_chrom_df["IsNeuralTranscript"] = mean_tpm_per_chrom_df.apply(
        lambda x: (x["MeanNeuralTPM"] >= 1)
        and (x["MeanNeuralTPM"] >= x["MeanNonNeuralTPM"] * 4),
        axis=1,
    )
    mean_tpm_per_chrom_df.to_csv(out_table, index=False, sep="\t")


def neural_obim_to_ovul_transcripts(
    obim_neural_transcripts_file,
    obim_to_ovul_orthologs_file,
    ovul_neural_transcripts_file,
):
    """Determine neural transcripts in O. vulgaris based on their O. bimaculoides' orthologs."""

    # obim_neural_transcripts_file = (
    #     "/private7/projects/Combinatorics/O.bim/Expression/NeuralTranscripts.tsv"
    # )
    # obim_to_ovul_orthologs_file = "/private7/projects/Combinatorics/O.vulgaris/Annotations/OrthoFinderAgainstObim/OrthoFinderResults/Results_May18/Orthologues/Orthologues_O.bim/O.bim__v__O.vul.tsv"

    obim_neural_transcripts_df = pd.read_table(obim_neural_transcripts_file).rename(
        columns={"Chrom": "ObimChrom"}
    )
    # obim_neural_transcripts_df

    obim_to_ovul_orthologs_df = (
        pd.read_table(obim_to_ovul_orthologs_file)
        .rename(columns={"O.bim": "ObimChrom", "O.vul": "OvulChrom"})
        .drop("Orthogroup", axis=1)
    )
    obim_to_ovul_orthologs_df["ObimChrom"] = obim_to_ovul_orthologs_df[
        "ObimChrom"
    ].str.split(",")
    obim_to_ovul_orthologs_df["OvulChrom"] = obim_to_ovul_orthologs_df[
        "OvulChrom"
    ].str.split(",")
    obim_to_ovul_orthologs_df = obim_to_ovul_orthologs_df.explode("ObimChrom").explode(
        "OvulChrom", ignore_index=True
    )
    obim_to_ovul_orthologs_df = obim_to_ovul_orthologs_df.merge(
        obim_neural_transcripts_df.loc[:, ["ObimChrom", "IsNeuralTranscript"]],
        on="ObimChrom",
        how="left",
    )

    # obim_to_ovul_orthologs_df = obim_to_ovul_orthologs_df.drop_duplicates(
    #     ["OvulChrom", "IsNeuralTranscript"], ignore_index=True
    # )
    # obim_to_ovul_orthologs_df[
    #     "IsOvulNeuralAnnotationConsistent"
    # ] = obim_to_ovul_orthologs_df.groupby("OvulChrom")["IsNeuralTranscript"].transform(
    #     lambda x: len(x.unique()) == 1
    # )

    # ovul_neural_transcripts_df = (
    #     obim_to_ovul_orthologs_df.loc[
    #         obim_to_ovul_orthologs_df["IsOvulNeuralAnnotationConsistent"],
    #         ["OvulChrom", "IsNeuralTranscript"],
    #     ]
    #     .reset_index(drop=True)
    #     .rename(columns={"OvulChrom": "Chrom"})
    # )

    ovul_neural_transcripts_df = (
        obim_to_ovul_orthologs_df.replace({"IsNeuralTranscript": {True: 1, False: 0}})
        .groupby("OvulChrom")
        .agg({"IsNeuralTranscript": ["sum", "size"]})
        .reset_index()
    )

    ovul_neural_transcripts_df = ovul_neural_transcripts_df.set_axis(
        ["OvulChrom", "NeuralObimOrthologs", "ObimOrthologs"], axis=1
    )
    ovul_neural_transcripts_df["NeuralObimOrthologs/ObimOrthologs"] = (
        ovul_neural_transcripts_df["NeuralObimOrthologs"]
        / ovul_neural_transcripts_df["ObimOrthologs"]
    )

    # def is_neural(neural_to_total_orthologs_ratio):
    #     if neural_to_total_orthologs_ratio == 1:
    #         return "Yes"
    #     elif neural_to_total_orthologs_ratio >= 0.5:
    #         return "Maybe"
    #     else:
    #         return "No"

    def is_neural(neural_to_total_orthologs_ratio):
        if neural_to_total_orthologs_ratio >= 0.5:
            return "Yes"
        return "No"

    ovul_neural_transcripts_df["IsNeural"] = ovul_neural_transcripts_df[
        "NeuralObimOrthologs/ObimOrthologs"
    ].apply(is_neural)

    ovul_neural_transcripts_df.to_csv(
        ovul_neural_transcripts_file, sep="\t", index=False
    )


def find_obim_to_ovul_orthologs_file(orthofineder_out_dir):
    # orthofineder_out_dir = Path(
    #     "/private7/projects/Combinatorics/O.vulgaris/Annotations/OrthoFinderAgainstObim/OrthoFinderResults/"
    # )
    obim_to_ovul_orthologs_file = Path(
        list(Path(orthofineder_out_dir, "OrthoFinderResults").iterdir())[0],
        "Orthologues/Orthologues_O.bim/O.bim__v__O.vul.tsv",
    )
    return obim_to_ovul_orthologs_file


def main(
    *,
    obim_trinity_file: Path,
    ovul_trinity_file: Path,
    # gtf_file: Union[Path, None],
    processes: int,
    threads: int,
    obim_salmon_out_dir: Path,
    obim_in_fastq_dir: Path,
    postfix: str,
    mate_prefix: str,
    recursive: bool,
    sra_run_table: Path,
    orthofinder_out_dir: Path,
    ovul_neural_transcripts_file: Path,
    **kwargs,
):
    # # 1 - prepare salmon gtf input from transcriptome

    # if gtf_file is None:
    #     gtf_file = Path(trinity_file.parent, f"{trinity_file.stem}.gtf")
    # prepare_gtf_for_salmon(trinity_file, gtf_file)

    # 2 - index the transcriptome
    obim_salmon_out_dir.mkdir(exist_ok=True)
    index_dir = Path(obim_salmon_out_dir, f"{obim_trinity_file.stem}_salmon_index")
    # we use threads*processes because this is the number of threads used in total,
    # while this is a single process
    max_threads = threads * processes
    create_salmon_index(obim_trinity_file, index_dir, max_threads)

    # 3 - run salmon

    # find fastq files and group paired reads by sample name
    fastq_files = find_files(obim_in_fastq_dir, postfix, recursive)
    paired_fastq_files = group_pe_fastq_files(fastq_files, postfix, mate_prefix)
    # run salmon quantification for each sample
    per_sample_out_dirs = [
        Path(
            obim_salmon_out_dir,
            extract_sample_name(fastq_1, postfix).split(mate_prefix)[0],
        )
        for fastq_1, _ in paired_fastq_files
    ]
    with Pool(processes=processes) as pool:
        pool.starmap(
            func=run_salmon,
            iterable=[
                (index_dir, fastq_1, fastq_2, threads, per_sample_out_dir)
                for (fastq_1, fastq_2), per_sample_out_dir in zip(
                    paired_fastq_files, per_sample_out_dirs
                )
            ],
        )

    # 4 - determine neural transcripts in O.bim

    obim_neural_transcripts_file = Path(obim_salmon_out_dir, "NeuralTranscripts.tsv")
    determine_obim_neural_transcripts(
        sra_run_table, per_sample_out_dirs, obim_neural_transcripts_file
    )

    # 5 - determine neural transcripts in O.vul by finding orthologs to neural transcripts in O.bim
    # 5.1 - find orthologs
    orthofinder(
        ovul_trinity_file,
        obim_trinity_file,
        "O.vul",
        "O.bim",
        orthofinder_out_dir,
        parallel_sequence_search_threads=max_threads,
        parallel_analysis_threads=None,
    )
    # 5.2 - get required orthologs file
    obim_to_ovul_orthologs_file = find_obim_to_ovul_orthologs_file(orthofinder_out_dir)
    # 5.3 - determine neural transcripts
    neural_obim_to_ovul_transcripts(
        obim_neural_transcripts_file,
        obim_to_ovul_orthologs_file,
        ovul_neural_transcripts_file,
    )


def define_args() -> argparse.ArgumentParser:
    # create common & sub parsers

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "--obim_trinity_file",
        help="Path to trinity transcriptome file of O.bim",
        type=abs_path_from_str,
        required=True,
    )
    parser.add_argument(
        "--ovul_trinity_file",
        help="Path to trinity transcriptome file of O.vul",
        type=abs_path_from_str,
        required=True,
    )
    # parser.add_argument(
    #     "--gtf_file",
    #     help="Path to output gtf file based on `trinity_file`. Defaults to `Path(trinity_file.parent, f'{trinity_file.stem}.gtf')`",
    #     type=abs_path_from_str,
    #     # required=True,
    # )
    parser.add_argument(
        "--processes",
        type=int,
        default=8,
        help="Maximal number of processes to run in parallel.",
    )
    parser.add_argument(
        "--threads", type=int, default=6, help="Threads used in each process."
    )
    parser.add_argument(
        "--obim_in_fastq_dir",
        required=True,
        type=abs_path_from_str,
        help=(
            "A folder with input fastq files, where each sample's name format is something like "
            "`$sample_{1, 2}.fastq.gz`. See `postfix` below."
        ),
    )
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
    parser.add_argument("--obim_salmon_out_dir", type=abs_path_from_str, required=True)
    parser.add_argument(
        "--sra_run_table",
        type=abs_path_from_str,
        default=Path("Code/NeuralTranscripts/SraRunTable.txt").absolute(),
    )
    parser.add_argument("--orthofinder_out_dir", type=abs_path_from_str, required=True)
    parser.add_argument(
        "--ovul_neural_transcripts_file", type=abs_path_from_str, required=True
    )

    return parser


if __name__ == "__main__":
    # run
    # https://stackoverflow.com/a/35824590/10249633 argparse.Namespace -> Dict
    parser = define_args()
    args = parser.parse_args()
    main(**vars(args))

    # end
    final_words()
