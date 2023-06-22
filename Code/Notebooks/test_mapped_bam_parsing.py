# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %%
mapped_bams = [
    "/private7/projects/Combinatorics/D.pealeii/Alignment/BestN1/GRIA-CNS-RESUB.C0x1291.aligned.sorted.bam",
    "/private7/projects/Combinatorics/D.pealeii/Alignment/BestN1/PCLO-CNS-RESUB.C0x1291.aligned.sorted.bam",
]
condition_col = "Gene"
conditions = ["GRIA", "PCLO"]

# %%
# from pathlib import Path
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import pysam
from icecream import ic

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
for mapped_bam in mapped_bams:
    with pysam.AlignmentFile(mapped_bam, "rb") as samfile:
        read_lens = []
        nps = []
        rqs = []
        deletion_events = []
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
        all_read_lens.append(read_lens)
        all_nps.append(nps)
        all_rqs.append(rqs)
        all_deletion_events.append(deletion_events)

# %%
dfs = [
    (
        pd.DataFrame(
            {
                "ReadLen": read_lens,
                "NumOfPasses": nps,
                "ReadQuality": rqs,
                "Deletions": deletions,
            }
        )
        .assign(condition_col=condition)
        .rename(columns={"condition_col": condition_col})
    )
    for read_lens, nps, rqs, condition, deletions in zip(
        all_read_lens, all_nps, all_rqs, conditions, all_deletion_events
    )
]

df = pd.concat(dfs).reset_index(drop=True)

df

# %%
# reads w/o any deletion event
(df.loc[df["Deletions"] == 0].groupby(condition_col).size())

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
    labels={"ReadLen": "CCS read length", "Gene": "Transcript"},
    title="Squid's PacBio",
    color=condition_col,
    color_discrete_map=color_discrete_map,
    category_orders=category_orders,
    template=template,
)

# https://stackoverflow.com/questions/58167028/single-axis-caption-in-plotly-express-facet-plot
for axis in fig.layout:
    if type(fig.layout[axis]) == go.layout.YAxis:
        fig.layout[axis].title.text = ""
fig.update_layout(
    showlegend=False, yaxis_title="Accumulated <br>% of reads", width=600, height=350, title_x=0.13,
)

fig.write_image(
    "Accumulated % of reads length - PacBio.svg",
    width=600,
    height=350,
)

fig.show()

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
fig = px.histogram(
    df,
    x="ReadQuality",
    y="Deletions",
    # facet_col=condition_col,
    # facet_col_spacing=facet_col_spacing,
    opacity=0.5,
    barmode="overlay",
    marginal="box",
    histfunc="avg",
    labels={"Deletions": "deletion events", "ReadQuality": "Read quality", "Gene": "Transcript"},
    # title="Occurrence of deletion events (regardless of their length) vs. read quality",
    title="Squid's PacBio",
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
fig.update_layout(
    title_x=0.17,
    yaxis_title="Deletion events (avg)",
    width=600,
    height=350,
)

fig.write_image(
    "Avg deletion events vs read quality - PacBio.svg",
    width=600,
    height=350,
)

fig.show()

# %%

# %%

# %%
n = 20

sites = list(range(1, n + 1))
isoforms = [2**x for x in sites]

df = pd.DataFrame({"Editing sites": sites, "Max possible isoforms": isoforms})
df[">=15"] = df["Editing sites"] >= 15

fig = px.scatter(
    df,
    x="Editing sites",
    y="Max possible isoforms",
    log_y=True,
    color=">=15",
    color_discrete_sequence=["green", "orange"],
    symbol=">=15",
    template=template,
    # height=400
)
fig.update_traces(
    marker=dict(size=8),
    selector=dict(mode="markers"),
)
fig.update_layout(showlegend=False)
fig.show()
