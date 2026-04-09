# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.18.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %%
mapped_bams = [
    "/private6/projects/Combinatorics/D.pealeii/Alignment/AdditionalUMILongReads/JR1.aligned.sorted.bam",
    "/private6/projects/Combinatorics/D.pealeii/Alignment/AdditionalUMILongReads/JR2.aligned.sorted.bam",
    "/private6/projects/Combinatorics/D.pealeii/Alignment/AdditionalUMILongReads/JR3.aligned.sorted.bam"
]
condition_col = "Gene"
conditions = ["JR1", "JR2", "JR3"]

# %%
allowed_chroms = ["comp141693_c0_seq1", "comp134400_c0_seq1_extended", "comp141565_c6_seq3"]

# %%
# from pathlib import Path
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pysam
from icecream import ic
import scipy.stats

# %%
# plotly consts
# color_sequence = px.colors.qualitative.Pastel
# color_sequence = px.colors.qualitative.D3
color_sequence = px.colors.qualitative.G10
color_discrete_map = {
    condition: color for condition, color in zip(conditions, color_sequence)
}
category_orders = {condition_col: conditions}
facet_col_spacing = 0.05
template = "plotly_white"


# %% [markdown]
# * `np:i:34` - NumPasses (1 for subreads, variable for CCS---encodes number of complete passes of the insert)
# * `rq:f:1` - Float in [0, 1] encoding expected accuracy

# %%
def count_deletions(read: pysam.AlignedSegment):
    cigar = read.cigarstring
    return sum(1 for x in cigar if x == "D")


# %%
all_read_lens = []
all_nps = []
all_rqs = []
all_deletion_events = []
all_mapped_chroms = []
for mapped_bam in mapped_bams:
    with pysam.AlignmentFile(mapped_bam, "rb") as samfile:
        read_lens = []
        nps = []
        rqs = []
        deletion_events = []
        mapped_chroms = []
        for read in samfile:
            # print(read.query_name)
            # rq_tag = read.get_tag("rq")
            # print(f"{rq_tag = }, {type(rq_tag) = }")
            # break
            read_lens.append(len(read.seq))
            deletion_events.append(count_deletions(read))
            tags = dict(read.tags)
            nps.append(int(tags["np"]))
            rqs.append(float(tags["rq"]))
            mapped_chroms.append(read.reference_name)
        all_read_lens.append(read_lens)
        all_nps.append(nps)
        all_rqs.append(rqs)
        all_deletion_events.append(deletion_events)
        all_mapped_chroms.append(mapped_chroms)

# %%
dfs = [
    (
        pd.DataFrame(
            {
                "MappedChrom": mapped_chroms,
                "ReadLen": read_lens,
                "NumOfPasses": nps,
                "ReadQuality": rqs,
                "Deletions": deletions,
            }
        )
        .assign(condition_col=condition)
        .rename(columns={"condition_col": condition_col})
    )
    for read_lens, nps, rqs, condition, deletions, mapped_chroms in zip(
        all_read_lens, all_nps, all_rqs, conditions, all_deletion_events, all_mapped_chroms
    )
]

df = pd.concat(dfs).reset_index(drop=True)

df

# %%
df.groupby(condition_col)["MappedChrom"].value_counts().reset_index()

# %%
(
    df.groupby(condition_col)["MappedChrom"].value_counts().reset_index()
    .groupby(condition_col).head(3)
    .sort_values([condition_col, "count"], ascending=[True, False])
)

# %%
df["MappedChrom"].value_counts()

# %%
df["ReadQuality"].mean().round(4)

# %%
df["ReadQuality"].median().round(4)

# %%

# %%
np.round(np.percentile(df["ReadQuality"], [25, 75]), 4)

# %%
0.9999 - 0.9986

# %%
pd.Series([scipy.stats.iqr(df["ReadQuality"])]).round(4)

# %%
df.groupby(condition_col)["ReadQuality"].mean().round(4)

# %%
df.groupby(condition_col)["ReadQuality"].median().round(4)

# %%
# IQR of reads with quality score >= 0.998
df.groupby(condition_col)["ReadQuality"].apply(lambda x: scipy.stats.iqr(x)).round(4)

# %%
# IQR of reads with quality score >= 0.998
df.groupby(condition_col).apply(lambda x: np.round(np.percentile(x["ReadQuality"], [25, 75]), 4))

# %%
# % of reads with quality score >= 0.998
100 - df.groupby(condition_col)["ReadQuality"].apply(lambda x: scipy.stats.percentileofscore(x, 0.998, kind="strict"))

# %%
# reads w/o any deletion event
(df.loc[df["Deletions"] == 0].groupby(condition_col).size())

# %%
fig = px.histogram(
    df,
    x="ReadQuality",
    facet_col=condition_col,
    facet_col_spacing=facet_col_spacing,
    histnorm="percent",
    cumulative=True,
    marginal="histogram",
    # opacity=0.5,
    # barmode="group",
    labels={"ReadQuality": "CCS read quality"},
    title="CCS read quality",
    color=condition_col,
    color_discrete_map=color_discrete_map,
    category_orders=category_orders,
    template=template,
)
fig.update_layout(showlegend=False)
fig.show()

# %%
fig = px.histogram(
    df,
    x="NumOfPasses",
    facet_col=condition_col,
    facet_col_spacing=facet_col_spacing,
    histnorm="percent",
    cumulative=True,
    marginal="histogram",
    labels={"NumOfPasses": "Number of passes"},
    title="Number of passes",
    color=condition_col,
    color_discrete_map=color_discrete_map,
    category_orders=category_orders,
    template=template,
)

# fig.update_layout(showlegend=False)

# https://stackoverflow.com/questions/58167028/single-axis-caption-in-plotly-express-facet-plot
for axis in fig.layout:
    if type(fig.layout[axis]) == go.layout.YAxis:
        fig.layout[axis].title.text = ""
fig.update_layout(showlegend=False, yaxis_title="Accumulated <br>% of reads")

fig.show()

# %%
fig = px.histogram(
    df,
    x="Deletions",
    facet_col=condition_col,
    facet_col_spacing=facet_col_spacing,
    histnorm="percent",
    cumulative=True,
    marginal="histogram",
    labels={"Deletions": "Deletion events"},
    title="Occurrence of deletion events (regardless of their length)",
    color=condition_col,
    color_discrete_map=color_discrete_map,
    category_orders=category_orders,
    template=template,
)
fig.update_layout(showlegend=False)
fig.show()

# %%
fig = px.scatter(
    df,
    x="NumOfPasses",
    y="ReadQuality",
    facet_col=condition_col,
    facet_col_spacing=facet_col_spacing,
    labels={"NumOfPasses": "Number of passes", "ReadQuality": "Read quality"},
    title="Number of passes vs. read quality",
    color=condition_col,
    color_discrete_map=color_discrete_map,
    category_orders=category_orders,
    template=template,
)
fig.update_layout(showlegend=False)
fig.show()

# %%
df

# %%
df.groupby(condition_col)["ReadLen"].median().round(1)

# %%
df.groupby(condition_col)["ReadLen"].apply(lambda x: scipy.stats.iqr(x)).round(1)

# %%
# IQR of reads with quality score >= 0.998
df.groupby(condition_col).apply(lambda x: np.round(np.percentile(x["ReadLen"], [25, 75]), 1))

# %%

# %%
100 * df["ReadLen"][df["ReadLen"] >= 3000].size / df["ReadLen"].size

# %%
df["ReadLen"].mean()

# %%
df["ReadLen"].std()

# %%
df["ReadLen"].median()

# %%
scipy.stats.iqr(df["ReadLen"])

# %%
fig = px.histogram(
    df,
    x="ReadLen",
    facet_col=condition_col,
    facet_col_spacing=facet_col_spacing,
    histnorm="percent",
    cumulative=True,
    marginal="histogram",
    # opacity=0.75,
    # barmode="group",
    labels={
    "ReadLen": "Read length", 
            # "Gene": "Transcript"
           },
    # title="Squid's PacBio",
    # title="Squid's Long-reads",
    title="Long-reads coverage",
    color=condition_col,
    color_discrete_map=color_discrete_map,
    category_orders=category_orders,
    template=template,
)

# https://stackoverflow.com/questions/58167028/single-axis-caption-in-plotly-express-facet-plot
for axis in fig.layout:
    if type(fig.layout[axis]) == go.layout.YAxis:
        fig.layout[axis].title.text = ""

fig.for_each_annotation(lambda a: a.update(text=a.text.replace("Gene=", "")))

fig.update_xaxes(tick0=0, dtick=1000, tickangle=45)

width=1000
height=400

fig.update_layout(
    showlegend=False, 
    # yaxis_title="Accumulated <br>% of reads", 
    # yaxis_title="Accumulated <br>reads [%]", 
    yaxis_title="Reads [%]",
    width=width, height=height, title_x=0.13,
)

# fig.write_image(
#     "Accumulated % of reads length - PacBio.svg",
#     width=width, height=height,
# )

fig.show()

# %%
fig = px.ecdf(
    df,
    x="ReadLen",
    # facet_col=condition_col,
    # facet_col_spacing=facet_col_spacing,
    ecdfnorm='percent',
    # cumulative=True,
    # marginal="histogram",
    # opacity=0.75,
    # barmode="group",
    labels={
    "ReadLen": "Read length", 
            # "Gene": "Transcript"
           },
    # title="Squid's PacBio",
    # title="Squid's Long-reads",
    title="Long-reads coverage",
    color=condition_col,
    color_discrete_map=color_discrete_map,
    category_orders=category_orders,
    template=template,
)

# # https://stackoverflow.com/questions/58167028/single-axis-caption-in-plotly-express-facet-plot
# for axis in fig.layout:
#     if type(fig.layout[axis]) == go.layout.YAxis:
#         fig.layout[axis].title.text = ""

# fig.for_each_annotation(lambda a: a.update(text=a.text.replace("Gene=", "")))

fig.update_xaxes(tick0=0, dtick=1000, tickangle=45)

width = 600
height = 500

fig.update_layout(
    # showlegend=False, 
    # yaxis_title="Accumulated <br>% of reads", 
    # yaxis_title="Accumulated <br>reads [%]", 
    yaxis_title="Reads [%]",
    width=width, height=height, title_x=0.13,
)

# fig.write_image(
#     "Accumulated % of reads length - PacBio.svg",
#     width=width, height=height,
# )

fig.show()

# %%
fig = px.ecdf(
    df.loc[df["MappedChrom"].isin(allowed_chroms)],
    x="ReadLen",
    # facet_col=condition_col,
    # facet_col_spacing=facet_col_spacing,
    ecdfnorm='percent',
    # cumulative=True,
    # marginal="histogram",
    # opacity=0.75,
    # barmode="group",
    labels={
    "ReadLen": "Read length", 
            # "Gene": "Transcript"
           },
    # title="Squid's PacBio",
    # title="Squid's Long-reads",
    title="Long-reads coverage in targeted chroms only",
    color=condition_col,
    color_discrete_map=color_discrete_map,
    category_orders=category_orders,
    template=template,
)

# # https://stackoverflow.com/questions/58167028/single-axis-caption-in-plotly-express-facet-plot
# for axis in fig.layout:
#     if type(fig.layout[axis]) == go.layout.YAxis:
#         fig.layout[axis].title.text = ""

# fig.for_each_annotation(lambda a: a.update(text=a.text.replace("Gene=", "")))

fig.update_xaxes(tick0=0, dtick=1000, tickangle=45)

width = 600
height = 500

fig.update_layout(
    # showlegend=False, 
    # yaxis_title="Accumulated <br>% of reads", 
    # yaxis_title="Accumulated <br>reads [%]", 
    yaxis_title="Reads [%]",
    width=width, height=height, title_x=0.13,
)

# fig.write_image(
#     "Accumulated % of reads length - PacBio.svg",
#     width=width, height=height,
# )

fig.show()

# %%
df.loc[
    (df["ReadLen"] >= 1000)
    & (df["MappedChrom"].isin(allowed_chroms))
].groupby(
    condition_col
).size()

# %%
df.loc[
    (df["ReadLen"] >= 1700)
    & (df["MappedChrom"].isin(allowed_chroms))
].groupby(
    condition_col
).size()

# %%
df.loc[
    (df["ReadLen"] >= 1700)
    & (df["MappedChrom"].isin(allowed_chroms))
].groupby(
    condition_col
)["MappedChrom"].value_counts()

# %%
df.loc[
    (df["ReadLen"] >= 1800)
    & (df["MappedChrom"].isin(allowed_chroms))
].groupby(
    condition_col
).size()

# %%
df.loc[
    (df["ReadLen"] >= 2000)
    & (df["MappedChrom"].isin(allowed_chroms))
].groupby(
    condition_col
).size()

# %%
df.loc[
    (df["ReadLen"] >= 3000)
    & (df["MappedChrom"].isin(allowed_chroms))
].groupby(
    condition_col
).size()

# %%
# fig = make_subplots(
#     rows=1,
#     cols=2,
#     subplot_titles=conditions,
#     shared_yaxes="all",
#     x_title="Read quality",
#     y_title="Deletion events (avg)",
# )

# fig = go.Figure()
fig = make_subplots(specs=[[{"secondary_y": True}]])

for condition in conditions:
    
    x = df.loc[df[condition_col] == condition, "ReadQuality"]
    y = df.loc[df[condition_col] == condition, "Deletions"]
    
    fig.add_trace(
            go.Histogram(
                x=x,
                y=y,
                marker_color=color_discrete_map[condition],
                name=condition,
                bingroup=1,
                histfunc="avg",
            ),
        )
    
    cum_reads_df = (
        df.loc[df[condition_col] == condition, ["ReadQuality"]]
        .sort_values("ReadQuality", ascending=False)
        .reset_index(drop=True)
    )
    cum_reads_df["%CummulativeReads"] = 100 * (cum_reads_df.index + 1) / len(cum_reads_df)
    # cum_reads_df["%CummulativeReads"] = cum_reads_df["%CummulativeReads"][::-1].values
    # cum_reads_df["%CummulativeReads"] = 100 - cum_reads_df["%CummulativeReads"]
    x = cum_reads_df["ReadQuality"]
    y = cum_reads_df["%CummulativeReads"]
    
    fig.add_trace(
            go.Scatter(
                x=x,
                y=y,
                mode="lines",
                marker_color=color_discrete_map[condition],
                # name=condition,
                # bingroup=1,
                # histfunc="avg",
                 # histnorm='percent',
                showlegend=False,
            ),
        secondary_y=True,
        )


fig.update_xaxes(title="Read quality")
fig.update_yaxes(title="Deletion events<br>(avg)", gridcolor='black',linewidth=5
                 # tick0=0, dtick=1
                )
fig.update_yaxes(
    title_text="Reads [%]", 
                 secondary_y=True, 
                 # tickmode="sync", 
                 range=[0, 100], 
                 tick0=0, dtick=25,
                 showgrid=True, 
                 gridcolor='LightPink',
                 griddash='dash',
                 linewidth=5
                )

fig.update_traces(opacity=0.5)

width = 550
height = 350

fig.update_layout(
    template=template,
    width=width,
    height=height,
    barmode="overlay",
    bargap=0.1,
    title="Long-reads quality",
    title_x=0.17,
    legend=dict(
        orientation="h", 
        x=0.85,
        y=0.7,
        xref="container",
        yref="container",
        xanchor="right",)
)

# fig.write_image(
#     "Avg deletion events vs read quality - PacBio.svg",
#     width=width,
#     height=height,
# )

fig.show()

# %%
healthy = 50
sick = 6

np.round(100 * sick / (sick + healthy), 2)

# %%
healthy = 60
sick = 6

np.round(100 * sick / (sick + healthy), 2)

# %%
healthy = 70
sick = 6

np.round(100 * sick / (sick + healthy), 2)

# %%
healthy = 80
sick = 6

np.round(100 * sick / (sick + healthy), 2)

# %%
healthy = 100
sick = 6

np.round(100 * sick / (sick + healthy), 2)

# %%
