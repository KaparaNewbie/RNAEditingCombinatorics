"""
Input:

x_fasta: a genome/transcriptome of an organism x whose RNA editing profile is unknown
y_fasta: a genome/transcriptome of an another organism y whose RNA editing profile is known
y_known_sites_bed: known editing sites of y

Output:

x_predicted_sites_bed: editing sites of x, predicted by aliging x_fasta and y_fasta
"""

from pathlib import Path
import subprocess

from icecream import ic

from General.os_utils import copy_text

x_fasta = Path("/private7/projects/Combinatorics/O.vulgaris/Annotations/chromosomer.fa")
y_fasta = Path("/private7/projects/Combinatorics/D.pealeii/Annotations/orfs_squ.fa")
out_dir = Path(
    "/private7/projects/Combinatorics/O.vulgaris/Annotations/tblastn.chromosomer.D.pea"
)
out_dir.mkdir(exist_ok=True)
x_fasta_copy = Path(out_dir, x_fasta.name)
y_fasta_copy = Path(out_dir, y_fasta.name)
x_y_alignment = Path(out_dir, "TBlastXResults.tsv")
# y_known_sites_bed = Path("/private7/projects/Combinatorics/D.pealeii/Annotations/D.pea.EditingSites.bed")
# x_predicted_sites_bed = ("/private7/projects/Combinatorics/O.vulgaris/Annotations/chromosomer.PredictedEditingSites.bed")
threads = 20

copy_text(x_fasta, x_fasta_copy)
copy_text(y_fasta, y_fasta_copy)

faidx_cmds = [f"samtools faidx {fasta}" for fasta in [x_fasta_copy, y_fasta_copy]]
ic(faidx_cmds)
for faidx_cmd in faidx_cmds:
    subprocess.run(faidx_cmd, shell=True, cwd=out_dir)

makeblastdb_cmd = f"makeblastdb -in {y_fasta_copy} -dbtype nucl -parse_seqids"

ic(makeblastdb_cmd)

subprocess.run(makeblastdb_cmd, shell=True)

tblastx_cmd = (
    "tblastx "
    f"-query {x_fasta_copy} "  # Input file name
    f"-db {y_fasta_copy} "  # BLAST database name
    f"-out {x_y_alignment} "
    '-outfmt "6 qaccver saccver pident length qframe mismatch gapopen qstart qend sstart send evalue bitscore" '
    f"-num_threads {threads} "  # Default = `1'
    "-evalue 1e-6"
)

ic(tblastx_cmd)

subprocess.run(tblastx_cmd, shell=True)

# we align seqs from y to x, and only consider


#  -max_target_seqs <Integer, >=1>
#    Maximum number of aligned sequences to keep
#    (value of 5 or more is recommended)
#    Default = `500'
#     * Incompatible with:  num_descriptions, num_alignments
#  -num_descriptions <Integer, >=0>
#    Number of database sequences to show one-line descriptions for
#    Not applicable for outfmt > 4
#    Default = `500'
#     * Incompatible with:  max_target_seqs
#  -num_alignments <Integer, >=0>
#    Number of database sequences to show alignments for
#    Default = `250'
#     * Incompatible with:  max_target_seqs
