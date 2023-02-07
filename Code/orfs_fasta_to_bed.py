import argparse

import pandas as pd
from Bio import SeqIO

from General.argparse_utils import abs_path_from_str


parser = argparse.ArgumentParser()
parser.add_argument("--in_fasta", required=True, type=abs_path_from_str)
parser.add_argument("--out_bed", required=True, type=abs_path_from_str)
args = parser.parse_args()

records = SeqIO.parse(args.in_fasta, "fasta")
transcripts = []
starts = []
ends = []
names = []
strands = []
for record in records:
    description = record.description.split("\t")
    start_index = description.index("OrfStart") + 1
    end_index = description.index("OrfEnd") + 1
    strand_index = description.index("Strand") + 1
    orf_start = int(description[start_index]) - 1
    orf_end = int(description[end_index])
    orf_name = description[-1].split()[0].split("|")[-1]
    strand = description[strand_index]
    transcripts.append(record.id)
    starts.append(orf_start)
    ends.append(orf_end)
    names.append(orf_name)
    strands.append(strand)
#
transcripts_strands_df = pd.DataFrame(
    {
        "Chrom": transcripts,
        "Start": starts,
        "End": ends,
        "Name": names,
        "Strand": strands,
    }
)
transcripts_strands_df.insert(
    transcripts_strands_df.columns.get_loc("Strand"), "Score", "."
)
transcripts_strands_df = transcripts_strands_df.sort_values(["Chrom", "Start"])

transcripts_strands_df.to_csv(args.out_bed, sep="\t", index=False, header=False)
