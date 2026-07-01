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
# unmapped_bams = [
#     "/private7/projects/Combinatorics/D.pealeii/Data/CCS/WoDup/PCLO-CNS-RESUB.C0x1291.ccs.sam",
#     "/private7/projects/Combinatorics/D.pealeii/Data/CCS/WoDup/GRIA-CNS-RESUB.C0x1291.ccs.sam"
# ]
# condition_col = "Tissue"
# conditions = ["GRIA", "PCLO"]

# %%
unmapped_bams = [
    "/private7/projects/Combinatorics/D.pealeii/Data/CCS/BasicCCS/GRIA-CNS-RESUB.C0x1291.ccs.bam",
    "/private7/projects/Combinatorics/D.pealeii/Data/CCS/BasicCCS/PCLO-CNS-RESUB.C0x1291.ccs.bam",
]
condition_col = "Tissue"
conditions = ["GRIA", "PCLO"]

# %%
# from pathlib import Path
import pandas as pd
import plotly.express as px
import pysam

# %%
color_discrete_map = {
    condition: color for condition, color in zip(conditions, px.colors.qualitative.G10)
}
category_orders = {condition_col: conditions}

# %% [markdown]
# * `np:i:34` - NumPasses (1 for subreads, variable for CCS---encodes number of complete passes of the insert)
# * `rq:f:1` - Float in [0, 1] encoding expected accuracy

# %%
all_read_lens = []
all_nps = []
all_rqs = []
for unmapped_bam in unmapped_bams:
    with pysam.AlignmentFile(unmapped_bams[0], "rb", check_sq=False) as samfile:
        read_lens = []
        nps = []
        rqs = []
        for read in samfile:
            read_lens.append(len(read.seq))
            tags = dict(read.tags)
            nps.append(int(tags["np"]))
            rqs.append(float(tags["rq"]))
        all_read_lens.append(read_lens)
        all_nps.append(nps)
        all_rqs.append(rqs)

# %%
dfs = [
    (
        pd.DataFrame({"ReadLen": read_lens, "NumOfPasses": nps, "ReadQuality": rqs})
        .assign(condition_col=condition)
        .rename(columns={"condition_col": condition_col})
    )
    for read_lens, nps, rqs, condition in zip(
        all_read_lens, all_nps, all_rqs, conditions
    )
]
dfs[0]

# %%
df = pd.concat(dfs)
df

# %%
fig = px.histogram(
    df,
    x="ReadLen",
    facet_col=condition_col,
    labels={"ReadLen": "CCS read length"},
    title="CCS read length",
    color=condition_col,
    color_discrete_map=color_discrete_map,
    category_orders=category_orders,
)
fig.update_layout(showlegend=False)
fig.show()

# %%
fig = px.histogram(
    df,
    x="ReadQuality",
    facet_col=condition_col,
    labels={"ReadQuality": "CCS read quality"},
    title="CCS read quality",
    color=condition_col,
    color_discrete_map=color_discrete_map,
    category_orders=category_orders,
)
fig.update_layout(showlegend=False)
fig.show()

# %%
fig = px.histogram(
    df,
    x="NumOfPasses",
    facet_col=condition_col,
    labels={"NumOfPasses": "Number of passes"},
    title="Number of passes",
    color=condition_col,
    color_discrete_map=color_discrete_map,
    category_orders=category_orders,
)
fig.update_layout(showlegend=False)
fig.show()

# %%
fig = px.scatter(
    df,
    x="NumOfPasses",
    y="ReadQuality",
    facet_col=condition_col,
    labels={"NumOfPasses": "Number of passes", "ReadQuality": "Read quality"},
    title="Number of passes vs. read quality",
    color=condition_col,
    color_discrete_map=color_discrete_map,
    category_orders=category_orders,
)
fig.update_layout(showlegend=False)
fig.show()
