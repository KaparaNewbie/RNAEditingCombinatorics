# ---
# jupyter:
#   jupytext:
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
import pandas as pd
import numpy as np
import plotly.io as pio
import plotly.express as px

pio.templates.default = "plotly_white"

# %%
pidstat_file = "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/TestFixedExpression/expressionlevels.RUSC2.21.01.2026-10:51:31.pidstat"

# %%
df = pd.read_csv(pidstat_file, sep='\s+', skiprows=2)

# the titles show every time a sample is taken -
# we want to remove these rows
df = df.loc[
    ~df["#"].apply(lambda x: x.startswith("#"))
].reset_index(drop=True)


# the 1st column (#) is time in 12-hour format, 
# the 2nd column (Time) is AM/PM indicator
# we want to combine these two columns into a single datetime column
df = df.rename(
    columns={
        "#": "Time_12",
        "Time": "AM/PM"
    }
)

df.insert(
    0,
    "Time",
    df["Time_12"] + " " + df["AM/PM"]
)

df = df.drop(
    columns=["Time_12", "AM/PM"]
)

# now we can change the Time column to datetime
df["Time"] = pd.to_datetime(df["Time"], format="%I:%M:%S %p").dt.time

df["%CPU"] = pd.to_numeric(df["%CPU"], errors='coerce')
df["%MEM"] = pd.to_numeric(df["%MEM"], errors='coerce')
df["RSS"] = pd.to_numeric(df["RSS"], errors='coerce')

df.insert(
    df.columns.get_loc("RSS") + 1,
    "RSS_MB",
    df["RSS"] / 1024
)

df

# %%
df_tg = df.loc[
    df["TGID"].ne("-")
].reset_index(drop=True)
df_tg

# %%
df_tg["TGID"].value_counts()

# %%
# fig = px.box(
#     df, 
#     x="Time", 
#     y="%CPU"
# )
# fig.update_layout(
#     title="CPU usage over time",
#     xaxis_title="Time",
#     yaxis_title="% CPU",
#     width=800,
#     height=500
# )
# fig.show()

# %%
# Average effective cores
np.round(
    df_tg["%CPU"].mean() / 100,
    2
)

# %%
fig = px.area(
    df_tg, 
    x="Time", 
    y="%CPU"
)
fig.update_layout(
    title="CPU usage over time",
    xaxis_title="Time",
    yaxis_title="% CPU",
    width=800,
    height=500
)
fig.show()

# %%
fig = px.area(
    df_tg, 
    x="Time", 
    y="%MEM"
)
fig.update_layout(
    title="Memory usage over time",
    xaxis_title="Time",
    yaxis_title="% Memory",
    width=800,
    height=500
)
fig.show()

# %%
fig = px.area(
    df_tg, 
    x="Time", 
    y="RSS_MB"
)
fig.update_layout(
    title="Memory usage over time",
    xaxis_title="Time",
    yaxis_title="RSS (MB)",
    width=800,
    height=500
)
fig.show()
