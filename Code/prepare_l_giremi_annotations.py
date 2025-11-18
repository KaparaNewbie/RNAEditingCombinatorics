## %%
from pathlib import Path
import os
import shutil
import subprocess

import pandas as pd
import Bio
from Bio import SeqIO, SeqRecord
from icecream import ic


## %%
threads = 30

## %%
pwd = Path.cwd()
if pwd.name == "Code":
    os.chdir(pwd.parent)
pwd = Path.cwd()
print(pwd)

# make sure we are in the right directory
paths_in_pwd = list(pwd.glob("*"))
assert any(list(map(lambda x: "Code" in x.name, paths_in_pwd)))


## %%
original_annotations_dir = Path("D.pealeii/Annotations/Jan2025/").absolute()
l_giremi_annotations_dir = Path("D.pealeii/L-GIREMI/Annotations/").absolute()


# original_annotations_dir, l_giremi_annotations_dir


## %%
if l_giremi_annotations_dir.exists():
    rm_cmd = f"rm -rf {l_giremi_annotations_dir}"
    ic(rm_cmd)
    subprocess.run(rm_cmd, shell=True)

shutil.copytree(original_annotations_dir, l_giremi_annotations_dir)

## %%

# Shorten transcripts IDs - needed to be <= 13 for EDTA

original_fasta = Path(l_giremi_annotations_dir, "orfs_squ.fa")
short_ids_fasta = Path(l_giremi_annotations_dir, "orfs_squ.short_ids.fa")
long_to_short_ids_file = Path(l_giremi_annotations_dir, "long_to_short_ids.csv")

records = list(SeqIO.parse(original_fasta, "fasta"))

original_ids = []
shortened_ids = []
new_records = []

# i = 0
# record = records[i]

# description = record.description.split("\t")


for i, record in enumerate(records):
    original_id = record.id
    seq = record.seq
    # transcript_end = len(record.seq)
    shortened_id = f"chr_{i}"
    assert (
        len(shortened_id) <= 13
    ), f"Shortened ID {shortened_id} ({original_id}) is too long!"
    new_record = SeqRecord.SeqRecord(seq, id=shortened_id, description="")
    original_ids.append(original_id)
    shortened_ids.append(shortened_id)
    new_records.append(new_record)

with open(short_ids_fasta, "w") as output_file:
    SeqIO.write(new_records, output_file, "fasta-2line")

long_to_short_ids_df = pd.DataFrame({"LongID": original_ids, "ShortID": shortened_ids})

long_to_short_ids_df.to_csv(long_to_short_ids_file, index=False)


len(original_ids)
len(shortened_ids)
original_to_short_ids = {
    original_id: shortened_id
    for original_id, shortened_id in zip(original_ids, shortened_ids)
}

# short_ids_fasta = Path(l_giremi_annotations_dir, "orfs_squ.short_ids.fa")
short_ids_cds_fasta = Path(l_giremi_annotations_dir, "orfs_squ.short_ids.cds.fa")
long_ids_cds_bed = Path(l_giremi_annotations_dir, "orfs_squ.bed")

long_ids_cds_df = pd.read_csv(
    long_ids_cds_bed,
    sep="\t",
    names=["Chrom", "Start", "End", "Name", "Score", "Strand"],
)
long_ids_cds_df["ShortChrom"] = long_ids_cds_df["Chrom"].apply(
    lambda x: original_to_short_ids[x]
)
long_ids_cds_df = long_ids_cds_df.set_index("ShortChrom")


new_cds_records = []

for i, record in enumerate(new_records):
    shortened_id = record.id
    seq = record.seq
    orf_start, orf_end = long_ids_cds_df.loc[shortened_id, ["Start", "End"]]
    cds_seq = seq[orf_start:orf_end]
    new_record = SeqRecord.SeqRecord(cds_seq, id=shortened_id, description="")
    new_cds_records.append(new_record)

with open(short_ids_cds_fasta, "w") as output_file:
    SeqIO.write(new_cds_records, output_file, "fasta-2line")


# create orfs bed -> gtf -> gff

transcripts_bed_12_file = Path(
    l_giremi_annotations_dir, "transcripts_squ.bed"
)  # a bed file of the complete transcriptome

long_ids_cds_df = pd.read_csv(
    long_ids_cds_bed,
    sep="\t",
    names=["Chrom", "Start", "End", "Name", "Score", "Strand"],
).set_index("Chrom")

transcript_records = list(SeqIO.parse(original_fasta, "fasta"))

transcripts = []
transcripts_starts = []
transcripts_ends = []
names = []
strands = []
orfs_starts = []
orfs_ends = []

for record in transcript_records:
    # description = record.description.split("\t")
    # transcript_start = 0
    # transcript_end = len(record.seq)
    # transcript_name = description[-1].split()[0].split("|")[-1]
    # strand_index = description.index("Strand") + 1
    # strand = description[strand_index]
    # orf_start_index = description.index("OrfStart") + 1
    # orf_end_index = description.index("OrfEnd") + 1
    # orf_start = int(description[orf_start_index]) - 1
    # orf_end = int(description[orf_end_index])
    chrom = record.id
    transcript_start = 0
    transcript_end = len(record.seq)
    transcript_name = long_ids_cds_df.loc[chrom, "Name"]
    strand = long_ids_cds_df.loc[chrom, "Strand"]
    orf_start = long_ids_cds_df.loc[chrom, "Start"]
    orf_end = long_ids_cds_df.loc[chrom, "End"]

    # transcripts.append(record.id)
    transcripts.append(chrom)
    transcripts_starts.append(transcript_start)
    transcripts_ends.append(transcript_end)
    names.append(transcript_name)
    strands.append(strand)
    orfs_starts.append(orf_start)
    orfs_ends.append(orf_end)

transcripts_strands_df = pd.DataFrame(
    {
        "Chrom": transcripts,
        "Start": transcripts_starts,
        "End": transcripts_ends,
        "Name": names,
        "Strand": strands,
        "ThickStart": orfs_starts,
        "ThickEnd": orfs_ends,
        "itemRgb": ".",
        "blockCount": 1,
        "blockStarts": 0,
    }
)
transcripts_strands_df = transcripts_strands_df.sort_values(["Chrom", "Start"])
transcripts_strands_df.insert(
    transcripts_strands_df.columns.get_loc("Strand"), "Score", "."
)
transcripts_strands_df.insert(
    transcripts_strands_df.columns.get_loc("blockCount") + 1,
    "BlockSizes",
    transcripts_strands_df["ThickEnd"] - transcripts_strands_df["ThickStart"],
)
transcripts_strands_df.to_csv(
    transcripts_bed_12_file, sep="\t", index=False, header=False
)


## %%

# edta_cmd = "EDTA.pl " f"--genome {short_ids_fasta} " f"--threads {threads}"
# ic(edta_cmd)

# # just for edta execution, than return to pwd
# os.chdir(l_giremi_annotations_dir)

# perl = os.path.expanduser("~/anaconda3/envs/l_giremi/bin/perl")
# edta = os.path.expanduser("~/anaconda3/envs/l_giremi/share/EDTA/EDTA.pl")

# result = subprocess.run(
#     [
#         perl,
#         edta,
#         f"--genome={short_ids_fasta}",
#         f"--threads={threads}",
#         "--overwrite 1 --sensitive 1 --anno 1",
#     ],
#     # shell=False,
#     # check=True,
#     stdout=subprocess.PIPE,
#     stderr=subprocess.PIPE,
#     text=True,
# )


# print("STDOUT:\n", result.stdout)
# print("STDERR:\n", result.stderr)
# print("Exit code:", result.returncode)


# ## %%
