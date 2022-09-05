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

# %%
mapped_bams = [
    # "/private7/projects/Combinatorics/D.pealeii/Alignment/BestN1/GRIA-CNS-RESUB.C0x1291.aligned.sorted.bam",
    "/private7/projects/Combinatorics/D.pealeii/Alignment/BestN1/PCLO-CNS-RESUB.C0x1291.aligned.sorted.bam",
]
condition_col = "Gene"
conditions = [
    # "GRIA", 
    "PCLO"
]

# %%
# from pathlib import Path
import pandas as pd
import pysam
from icecream import ic

# %%
for mapped_bam in mapped_bams:
    with pysam.AlignmentFile(mapped_bam, "rb") as samfile:
        for read in samfile:
            if read.query_name == "m54278_210226_220806/20972305/ccs":
            # if 2100 <= read.query_alignment_start <= 2200:
                ic(read.query_name)
                # ic(read.query_alignment_start)
                # ic(read.query_sequence)
                ic(read.get_aligned_pairs(matches_only=False, with_seq=True)[:10])
                break
                # ic(read.to_string())
            # rq_tag = read.get_tag("rq")
            # print(f"{rq_tag = }, {type(rq_tag) = }")
                
    break

# %% [markdown]
# <!-- ![image.png](attachment:b1be2604-01be-447f-9493-f917862ca841.png)
# ![image.png](attachment:ba62ebeb-aa91-40f2-bfb7-edec7697bdfd.png) -->

# %% [markdown]
# ![image.png](attachment:4eb6eeaf-8c8d-4c41-88af-b18ffe7dd363.png)  
# ![image.png](attachment:3c3d5f27-8c6d-4d80-b07f-1d098be9fa40.png)  
# ![image.png](attachment:7f57ea4c-48de-4df7-960b-0a173fefec5b.png)

# %%
ref = "AGCATGTTAGATAA**GATAGCTGTGCTAGTAGGCAGTCAGCGCCAT"
ic(len(ref));

# %%
sam = pd.DataFrame(
    {
        "QNAME": ["r001", "r002", "r003", "r004", "r003", "r001"],
        "FLAG": [99, 0, 0, 0, 2064, 147],
        "RNAME": ["ref", "ref", "ref", "ref", "ref", "ref"],
        "POS": [7, 9, 9, 16, 29, 37],
        "MAPQ": [30, 30, 30, 30, 17, 30],
        # "CIGAR": ["8M2I4M1D3M", "3S6M1P1I4M", "5S6M", "6M14N5M", "6H5M", "9M"],
        "CIGAR": [
            [(8, "M"), (2, "I"), (4, "M"), (1, "D"), (3, "M")], 
            [(3, "S"), (6, "M"), (1, "P"), (1, "I"), (4, "M")], 
            [(5, "S"), (6, "M")], 
            [(6, "M"), (14, "N"), (5, "M")], 
            [(6, "H"), (5, "M")], 
            [(9, "M")]
        ],
        "RNEXT": ["=", "*", "*", "*", "*", "="],
        "PNEXT": [37, 0, 0, 0, 0, 7],
        "TLEN": [39, 0, 0, 0, 0, -39],
        "SEQ": ["TTAGATAAAGGATACTG", "AAAAGATAAGGATA", "GCCTAAGCTAA", "ATAGCTTCAGC", "TAGGC", "CAGCGGCAT"],
        "QUAL": ["*", "*", "*", "*", "*", "*"],
        "OPTIONAL": ["", "", "SA:Z:ref,29,-,6H5M,17,0;", "", "SA:Z:ref,9,+,5S6M,30,1;", "NM:i:1"],
        "ALNSEQ": ["TTAGATAAAGGATA*CTG", "aaaAGATAA*GGATA", "gcctaAGCTAA", "ATAGCT..............TCAGC", "ttagctTAGGC", "CAGCGGCAT"]
    }
)
sam

# %%
# proof that the cigar operations describe the entire aligned seq - lengthwise
for row in sam.itertuples(index=False):
    ops_len = sum(op_len for op_len, _ in row.CIGAR)
    assert ops_len == len(row.ALNSEQ)

# %%
for row in sam.itertuples(index=False):
    seq = row.SEQ
    alnseq = row.ALNSEQ
    parsed_alnseq = ""
    x = 0
    parsed_alnseq_addition = ""
    x_addition = 0
    for op_len, op in row.CIGAR:
        # Op, Description, Consumes query, Consumes reference
        if op in ["M", "=", "X"]: # alignment match (can be a sequence match or mismatch) / sequence match / sequence mismatch, yes, yes
            parsed_alnseq_addition = seq[x:x + op_len]
            x_addition = op_len
        elif op == "I": # insertion to the reference, yes, no
            parsed_alnseq_addition = seq[x:x + op_len]
            x_addition = op_len
        elif op == "D": # eletion from the reference, no, yes
            parsed_alnseq_addition = "*" * op_len
        elif op == "N": # skipped region from the reference, no, yes
            parsed_alnseq_addition = "." * op_len
        elif op == "S": # soft clipping (clipped sequences present in SEQ), yes, no
            parsed_alnseq_addition = seq[x:x + op_len].lower()
            x_addition = op_len
        elif op == "H": # hard clipping (clipped sequences NOT present in SEQ), no, no
            if parsed_alnseq == "":
                alnseq = alnseq[op_len:]
            else: # we should have consumed all query up to this point
                alnseq = alnseq[:op_len]
        elif op == "P": # padding (silent deletion from padded reference), no, no
            parsed_alnseq_addition = "*" * op_len
        else:
            raise Exception("Unknown CIGAR operation")
        parsed_alnseq += parsed_alnseq_addition
        x += x_addition
        parsed_alnseq_addition = ""
        x_addition = 0
    ic(alnseq, parsed_alnseq)
    assert alnseq == parsed_alnseq
        
    # break
