from dataclasses import dataclass
from pathlib import Path
from typing import Union
from collections import Counter

import Bio
from Bio import SeqRecord, SeqIO, Seq


def make_fasta_dict(fasta_file: Union[Path, str]) -> dict[str, Seq.Seq]:
    records = SeqIO.parse(Path(fasta_file).absolute(), "fasta")
    fasta_dict = {record.id: record.seq for record in records}
    return fasta_dict

AA_BY_TRIPLET = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "AGT": "S",
    "AGC": "S",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGA": "*",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "TGT": "C",
    "TGC": "C",
    "TGG": "W",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGA": "R",
    "AGG": "R",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


def extract_orfs_to_bed(genome: Union[Path, str], bed: Union[Path, str]) -> None:
    records = Bio.SeqIO.parse(Path(genome).absolute(), "fasta")
    orfs = []
    for record in records:
        description = record.description.split("\t")
        start_index = description.index("OrfStart") + 1
        end_index = description.index("OrfEnd") + 1
        strand_index = description.index("Strand") + 1
        start = str(int(description[start_index]) - 1)
        end = description[end_index]
        strand = description[strand_index]
        orf = [record.id, start, end, ".", ".", strand]
        orfs.append(orf)
        # seq = record.seq[orf_start:orf_end]
    with Path(bed).absolute().open("w") as bed:
        for orf in orfs:
            bed.write("\t".join(orf) + "\n")


# @dataclass
# class Result:
#     record: SeqRecord.SeqRecord
#     start: int
#     end: int
#     strand: str
#     seq: Seq.Seq
#     protein = None


# # the following code shows that the orfs are problamatic

# genome = (
#     "/private7/projects/Combinatorics/D.pealeii/Annotations/December2017/orfs_squ.fa"
# )
# # records = list(Bio.SeqIO.parse(Path(genome).absolute(), "fasta"))
# records = list(SeqIO.parse(Path(genome).absolute(), "fasta"))
# good_orfs = []
# bad_orfs = []
# mod_3_lengths = Counter()
# all_proteins = []
# for record in records:
#     description = record.description.split("\t")
#     start_index = description.index("OrfStart") + 1
#     end_index = description.index("OrfEnd") + 1
#     strand_index = description.index("Strand") + 1
#     start = int(description[start_index]) - 1
#     end = int(description[end_index])
#     strand = description[strand_index]
#     seq = record.seq[start:end]
#     mod_3_lengths[len(seq) % 3] += 1
#     seq = seq.reverse_complement() if strand == "-" else seq
#     try:
#         protein = seq.translate(to_stop=True, cds=True)
#         result = Result(record, start, end, strand, seq, protein)
#         good_orfs.append(result)
#         all_proteins.append(protein)
#     except:
#         # protein = seq.translate(to_stop=True)
#         # result = Result(record, start, end, strand, protein)
#         result = Result(record, start, end, strand, seq)
#         bad_orfs.append(result)
#         protein = seq.translate()
#         all_proteins.append(protein)


# bad_first_aa = Counter(orf.seq[:3].translate() for orf in bad_orfs)
