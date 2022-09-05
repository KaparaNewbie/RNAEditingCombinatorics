# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.0
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# This notebook is an attempt to parse the mpileup file created by the command:
#
# ```
# cd /private7/projects/Combinatorics
#
# samtools \
# mpileup \
# --fasta-ref D.pealeii/Annotations/orfs_squ.fa \
# --positions D.pealeii/Annotations/D.pea.EditingSites.bed \
# --region comp141882_c0_seq14 \
# --no-BAQ \
# --no-output-ins --no-output-ins \
# --no-output-del --no-output-del \
# --no-output-ends \
# --output-QNAME \
# --excl-flags 2304 \
# D.pealeii/Alignment/Naive/PCLO-CNS-RESUB.C0x1291.aligned.sorted.bam \
# --output D.pealeii/Alignment.Test/PCLO.primary.linear.mpileup
# ```
# If we can retain only base match / base mismatch in the 5th col, and corresponding number of reads' names in the last col, then we are done.

# %%
import pandas as pd

# %% [markdown]
# ## positions

# %%
# pileup_file = "/private7/projects/Combinatorics/D.pealeii/Alignment.Test/PCLO.primary.linear.mpileup"
# pileup_file = "/private7/projects/Combinatorics/D.pealeii/Alignment.Test.N1/PCLO.primary.linear.mpileup"
pileup_file = "/private7/projects/Combinatorics/D.pealeii/Alignment.Test.N1/PCLO.primary.linear.max_depth_90k.mpileup"

# %%
cols = ["Chrom", "Position", "RefBase", "TotalCoverage", "MappedBases", "Phred", "Reads"]
df = pd.read_csv(pileup_file, sep="\t", names=cols)
df["Position"] = df["Position"] - 1
df

# %%
unique_reads_in_raw_df = {
    read 
    for reads in df["Reads"]
    for read in reads.split(",")
}
len(unique_reads_in_raw_df)

# %%
for row_num, row in enumerate(df.itertuples()):
    bases = row.MappedBases
    reads = row.Reads.split(",")
    if len(bases) != len(reads):
        print(row_num)
        break

# %% [markdown]
# OK! That's working!  
# So, as long I add this sanity check to the actual script, I'm good.  Also, it should be documented if someone uses the  
# script with different data.

# %% [markdown]
# Now, let's aggregate editing sites to transcripts.

# %%
# keep only positions where at least one of the reads had editing
editing_in_position = []
for row in df.itertuples():
    position = row.Position
    bases = row.MappedBases
    if "G" in bases or "g" in bases:
        editing_in_position.append(True)
    else:
        editing_in_position.append(False)
df = df.loc[editing_in_position]
df = df.reset_index(drop=True)
df

# %% [markdown]
# ## reads

# %%
unique_reads = {
    read 
    for reads in df["Reads"]
    for read in reads.split(",")
}
len(unique_reads)

# %%
edited_positions_in_reads: dict[str, list[bool]] = {read: [] for read in unique_reads}

for x, row in enumerate(df.itertuples()):
    
    # position = row.Position
    bases = row.MappedBases
    reads = row.Reads.split(",")
    
    # deal with reads that were mapped to this pos
    for base, read in zip(bases, reads):
        edited_in_pos = True if base in ["G", "g"] else False
        edited_positions_in_reads[read].append(edited_in_pos)
    
    # deal with reads that weren't mapped to this pos
    unmapped_reads_in_pos = unique_reads - set(reads)
    for read in unmapped_reads_in_pos:
        edited_in_pos = False
        edited_positions_in_reads[read].append(edited_in_pos)
        
    # check all reads got values for pos
    mapped_pos_dist = {len(bases) for bases in edited_positions_in_reads.values()}
    if len(mapped_pos_dist) != 1:
        print(x)
        break

# %%
df2 = pd.DataFrame(edited_positions_in_reads)
df2 = df2.T.reset_index().rename({'index':'Read'}, axis = 'columns')

df2 = df2.rename(
    columns = {
        old_col: f"{chrom}:{pos}"
        for old_col, chrom, pos 
        in zip(df2.columns[1:], df["Chrom"], df["Position"])
    }
)

df2

# %% [markdown]
# ## transcripts

# %% [markdown]
# So every read has its edited sites annotated.   
# Now, let's finaly aggregate reads into transcripts.

# %%
df3 = df2.groupby(df2.columns[1:].tolist()).size().to_frame('count').reset_index()
df3.insert(0, "Count", df3["count"])
del df3["count"]
df3 = df3.sort_values("Count", ascending=False).reset_index(drop=True)
df3

# %% [markdown]
# Adding the original reads supporting each transcript:

# %%
df3.iloc[0, 1:].isin(df3.iloc[0, 1:]).any()

# %%
df3.iloc[0, 1:].eq(df3.iloc[0, 1:]).all()

# %%
# the first line contain transcripts without any editing activity in their known editing sites
line_from_df3 = df3.iloc[0, 1:]
line_from_df3.any()

# %%
line_from_df3 = df3.iloc[1, 1:]
line_from_df3.any()

# %%
line_from_df3 = df3.iloc[1]
line_from_df3.Count

# %%
df3_line_wo_count_col = df3.iloc[0, 1:]
df3_line_wo_count_col

# %%
df2_wo_read_col = df2.iloc[:, 1:]
df2_wo_read_col

# %%
df3_line_wo_count_col.eq(df2_wo_read_col).all(axis="columns")

# %%
len([x for x in df3_line_wo_count_col.eq(df2_wo_read_col).all(axis="columns") if x])

# %%
df2.loc[df3_line_wo_count_col.eq(df2_wo_read_col).all(axis="columns")]

# %%
# lines_with_supporting_reads = df2.loc[df3_line_wo_count_col.eq(df2_wo_read_col).all(axis="columns")]
# ",".join(read for read in lines_with_supporting_reads["Read"])

# %%
# def get_supporting_reads(df3_line_wo_count_col, df2):
#     df2_wo_read_col = df2.iloc[:, 1:]
#     lines_with_supporting_reads = df2.loc[df3_line_wo_count_col.eq(df2_wo_read_col).all(axis="columns")]
#     # return [read for read in lines_with_supporting_reads["Read"]]  # todo return str representation (the line below)
#     return ",".join(read for read in lines_with_supporting_reads["Read"])

# %%
# processes = 6
# df3_lines_by_proccesses = list(range(0, len(df3), int(len(df3) / processes)))
# df3_lines_by_proccesses

# %%
# processes = 6
# from multiprocessing import Pool
# with Pool(processes=processes) as pool:
#     all_supporting_reads = pool.startmap(
#         func = get_supporting_reads,
#         Iterable = [(df3.iloc[df3_line_num, 1:], df2)
#                     for df3_line_num in range(len(df3))]
#     )  

# %%

# %%
# all_supporting_reads = []
# for i, df3_line_num in enumerate(range(len(df3))):
#     df3_line_wo_count_col = df3.iloc[df3_line_num, 1:]
    
#     line_supporting_reads = get_supporting_reads(df3_line_wo_count_col, df2)
    
#     # line_supporting_reads = []
#     # for df2_line_num in range(len(df2)):
#     #     df2_line = df2.iloc[df2_line_num]
#     #     if df3_line_wo_count_col.eq(df2_line.iloc[1:]).all():
#     #         line_supporting_reads.append(df2_line.Read)
    
#     all_supporting_reads.append(line_supporting_reads)
#     # if i == 2:
#     #     break

# %%
# print([len(line_supporting_reads) for line_supporting_reads in all_supporting_reads])

# %%

# %%
# all_supporting_reads = {}
# for i, df3_line_num in enumerate(range(len(df3))):
#     df3_line = df3.iloc[df3_line_num, 1:]
#     line_supporting_reads = []
#     for df2_line_num in range(len(df2)):
#         df2_line = df2.iloc[df2_line_num]
#         if df3_line.eq(df2_line.iloc[1:]).all():
#             line_supporting_reads.append(df2_line.Read)
#     all_supporting_reads[df3_line_num] = line_supporting_reads
#     if i == 2:
#         break

# %%
# all_supporting_reads

# %%
# [len(line_supporting_reads) for line_supporting_reads in all_supporting_reads.values()]

# %%

# %%

# %%
