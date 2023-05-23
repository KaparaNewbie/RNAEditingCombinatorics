import argparse
from pathlib import Path
import subprocess
from typing import Union
import sys

from Bio import SeqIO, SeqRecord

grandparent_dir = Path(__file__).absolute().parent.parent
sys.path.append(str(grandparent_dir))

from General.argparse_utils import abs_path_from_str
from General.consts import final_words


def transcriptome_to_proteome(transcriptome: Path, proteome: Path):
    protein_records = []
    transcript_records = SeqIO.parse(transcriptome, "fasta")
    for transcript_record in transcript_records:
        description = transcript_record.description.split("\t")

        start_index = description.index("OrfStart") + 1
        end_index = description.index("OrfEnd") + 1
        strand_index = description.index("Strand") + 1
        orf_start = int(description[start_index]) - 1
        orf_end = int(description[end_index])
        strand = description[strand_index]

        transcript_name = transcript_record.id
        protein_name = description[-1].split()[0].split("|")[-1]

        transcript_seq = transcript_record[orf_start:orf_end]
        if strand == "-":
            transcript_seq = transcript_seq.reverse_complement()
        protein_seq = transcript_seq.translate().seq

        protein_record = SeqRecord.SeqRecord(
            protein_seq, id=transcript_name, description=protein_name
        )
        protein_records.append(protein_record)

    with open(proteome, "w") as output_file:
        SeqIO.write(protein_records, output_file, "fasta")


def run_orthofinder(
    proteomes_dir: Path,
    out_dir: Path,
    parallel_sequence_search_threads: int,
    parallel_analysis_threads: int,
):
    orthofinder_cmd = f"orthofinder -f {proteomes_dir} -t {parallel_sequence_search_threads} -a {parallel_analysis_threads} -o {out_dir}"
    subprocess.run(orthofinder_cmd, shell=True)


def main(
    transcriptome_1: Path,
    transcriptome_2: Path,
    species_1: str,
    species_2: str,
    out_dir: Path,
    parallel_sequence_search_threads: int,
    parallel_analysis_threads: Union[int, None],
):
    out_dir.mkdir(exist_ok=True)
    proteomes_dir = Path(out_dir, "Proteomes")
    proteomes_dir.mkdir(exist_ok=True)
    transcriptome_to_proteome(transcriptome_1, Path(proteomes_dir, f"{species_1}.fa"))
    transcriptome_to_proteome(transcriptome_2, Path(proteomes_dir, f"{species_2}.fa"))
    orthofineder_out_dir = Path(out_dir, "OrthoFinderResults")
    if parallel_analysis_threads is None:
        parallel_analysis_threads = int(parallel_sequence_search_threads / 8)
    run_orthofinder(
        proteomes_dir,
        orthofineder_out_dir,
        parallel_sequence_search_threads,
        parallel_analysis_threads,
    )


def define_args() -> argparse.ArgumentParser:
    # create common & sub parsers

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "--transcriptome_1",
        help="Path to a trinity transcriptome file",
        type=abs_path_from_str,
        required=True,
    )
    parser.add_argument(
        "--transcriptome_2",
        help="Path to another trinity transcriptome file (from a different species)",
        type=abs_path_from_str,
        required=True,
    )
    parser.add_argument(
        "--species_1",
        help="Name of the first species",
        required=True,
    )
    parser.add_argument(
        "--species_2",
        help="Name of the second species",
        required=True,
    )
    parser.add_argument(
        "--out_dir",
        help="Path to output directory",
        type=abs_path_from_str,
        required=True,
    )
    parser.add_argument(
        "-t",
        "--parallel_sequence_search_threads",
        type=int,
        default=30,
        help=(
            "This option should always be used. It specifies the number of parallel processes for the "
            "BLAST/DIAMOND searches and tree inference steps. These steps represent most of the runtime and are "
            "highly-parallelisable and so you should typically use as many threads as there are cores available on "
            "your computer. This is the value it will default to if not specified by the user."
        ),
    )
    parser.add_argument(
        "-a",
        "--parallel_analysis_threads",
        type=int,
        help=(
            "In addition to the above, all of the critical internal steps of the OrthoFinder algorithm have been "
            "parallelised. The number of threads for these steps is controlled using the '-a' option. These steps "
            "typically have larger RAM requirements and so using a value 4-8x smaller than that used for the '-t' "
            "option is usually a good choice. Since these steps are a small component of the overall runtime it is not "
            "important to set '-a' as high as possible in order to get good performance. Not running out of RAM is a "
            "more important consideration. If the '-a' parameter is not set it will default to one eighth of the "
            "'-t' parameter."
        ),
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
