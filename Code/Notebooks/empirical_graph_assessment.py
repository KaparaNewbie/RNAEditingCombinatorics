# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.15.2
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Imports, settings, etc.

# %%
import re
import warnings
from collections import defaultdict, namedtuple
from itertools import product
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
from icecream import ic
from plotly.subplots import make_subplots

# from math impo

# %%
pio.templates.default = "plotly_white"

pd.set_option("display.max_columns", 500)

# %%
out_dir = Path("/private7/projects/Combinatorics/Simulations/GraphAssessment")
distinct_files_dir = Path(
    "/private7/projects/Combinatorics/Simulations/GraphAssessment/SameFrac"
)
false_positives_files_fofn = Path(out_dir, "FalsePositivesOutFiles.fofn")

n_reads = 100_000

unique_proteins_first_pos_col = 14

seed = 1892

fractions = [0.2, 0.4, 0.6, 0.8, 1.0]

# %% [markdown]
# # Recovered isoforms (positives)

# %% [markdown]
# ## Read results from disk into dfs

# %%
beginning_of_postfix = ".DistinctUniqueProteins."

# %%
# distinct_files = list(out_dir.glob("*.DistinctUniqueProteins.*.csv"))
# distinct_files = list(out_dir.glob(f"*{beginning_of_postfix}*.csv"))
distinct_files = list(distinct_files_dir.glob(f"*{beginning_of_postfix}*.csv"))
# distinct_files

# %%
distinct_files[0]

# %%
len(distinct_files)

# %%
unknown_probabilities = set()
swiss_prot_name_by_chrom = dict()
repetitions = set()

distinct_files_dict = {}
for distinct_file in distinct_files:
    distinct_file_name = re.sub(rf"{beginning_of_postfix}*.csv", "", distinct_file.name)

    # ic(distinct_file_name)

    file_name_parts = distinct_file_name.split(".")
    chrom = file_name_parts[0]
    # swiss_prot_name = file_name_parts[1].split("_")[0]
    swiss_prot_name = file_name_parts[1]
    unknown_probability = float(file_name_parts[2].replace("UP", "").replace("_", "."))
    repetition = file_name_parts[3]
    data_type = "+".join(file_name_parts[4::])
    distinct_files_dict[
        (chrom, swiss_prot_name, unknown_probability, repetition, data_type)
    ] = distinct_file

    unknown_probabilities.add(unknown_probability)
    swiss_prot_name_by_chrom[chrom] = swiss_prot_name
    repetitions.add(repetition)

    # break

unknown_probabilities = sorted(list(unknown_probabilities))
swiss_prot_name_by_chrom = dict(
    sorted(swiss_prot_name_by_chrom.items(), key=lambda x: x[1])
)

distinct_dfs = []
for (
    chrom,
    swiss_prot_name,
    unknown_probability,
    repetition,
    data_type,
), distinct_file in distinct_files_dict.items():
    df = pd.read_table(distinct_file)
    df.insert(0, "Chrom", chrom)
    df.insert(1, "SwissProtName", swiss_prot_name)
    df.insert(2, "UnknownProbability", unknown_probability)
    df.insert(3, "DataCreationRepetition", repetition)
    df.insert(4, "DataType", data_type)
    distinct_dfs.append(df)

# %%
swiss_prot_name_by_chrom

# %%
swiss_prot_name_by_chrom.values()

# %%
short_swiss_prot_names = [x.split("_")[0] for x in swiss_prot_name_by_chrom.values()]
short_swiss_prot_names

# %%
len(distinct_dfs)

# %%
distinct_dfs[0]

# %%
merged_distinct_df = pd.concat(distinct_dfs, ignore_index=True).rename(
    columns={
        "NumUniqueSamples": "NumDistinctProteins",
        "UniqueSamples": "DistinctProteins",
    }
)
merged_distinct_df.insert(
    6, "NumOfReads", merged_distinct_df["Fraction"].mul(n_reads).astype(int)
)
merged_distinct_df["SwissProtName"] = (
    merged_distinct_df["SwissProtName"].str.split("_").str[0]
)
merged_distinct_df.insert(
    merged_distinct_df.columns.get_loc("NumDistinctProteins") + 1,
    "IsMaxDistinctProteinsPerGroup",
    merged_distinct_df.groupby(
        [
            "Chrom",
            "SwissProtName",
            "UnknownProbability",
            "DataCreationRepetition",
            "DataType",
            "Fraction",
        ]
    )["NumDistinctProteins"]
    .transform(max)
    .eq(merged_distinct_df["NumDistinctProteins"]),
)
merged_distinct_df

# %%
# max_per_fraction_merged_distinct_df = (
#     merged_distinct_df.loc[merged_distinct_df["IsMaxDistinctProteinsPerGroup"],]
#     .groupby(["Chrom", "SwissProtName", "UnknownProbability", "DataType", "Fraction"])
#     .sample(n=1, random_state=seed)
#     .reset_index(drop=True)
#     .rename(columns={"NumDistinctProteins": "MaxNumDistinctProteins"})
# )
# max_per_fraction_merged_distinct_df

# %%
# # stats over DataCreationRepetition

# stats_of_merged_distinct_df = (
#     merged_distinct_df.groupby(
#         [
#             "Chrom",
#             "SwissProtName",  # same as Chrom
#             "UnknownProbability",
#             "DataType",
#             "Fraction",
#             "NumOfReads",  # same as Fraction
#         ]
#     )
#     .agg(
#         MeanNumDistinctProteins=pd.NamedAgg(
#             column="NumDistinctProteins", aggfunc="mean"
#         ),
#         StdNumDistinctProteins=pd.NamedAgg(column="NumDistinctProteins", aggfunc="std"),
#         MaxNumDistinctProteins=pd.NamedAgg(column="NumDistinctProteins", aggfunc="max"),
#     )
#     .reset_index()
# )

# stats_of_merged_distinct_df

# %%
max_per_fraction_per_data_creation_merged_distinct_df = (
    merged_distinct_df.loc[merged_distinct_df["IsMaxDistinctProteinsPerGroup"],]
    .groupby(
        [
            "Chrom",
            "SwissProtName",
            "UnknownProbability",
            "DataCreationRepetition",
            "DataType",
            "Fraction",
        ]
    )
    .sample(n=1, random_state=seed)
    .reset_index(drop=True)
    .rename(columns={"NumDistinctProteins": "MaxNumDistinctProteins"})
)
max_per_fraction_per_data_creation_merged_distinct_df

# %%
# stats over DataCreationRepetition

stats_of_max_per_fraction_per_data_creation_merged_distinct_df = (
    max_per_fraction_per_data_creation_merged_distinct_df.groupby(
        [
            "Chrom",
            "SwissProtName",  # same as Chrom
            "UnknownProbability",
            "DataType",
            # "DataCreationRepetition",
            "Fraction",
            "NumOfReads",  # same as Fraction
        ]
    )
    .agg(
        MeanOfMaxDistinctProteins=pd.NamedAgg(
            column="MaxNumDistinctProteins", aggfunc="mean"
        ),
        StdOfMaxDistinctProteins=pd.NamedAgg(
            column="MaxNumDistinctProteins", aggfunc="std"
        ),
        # MaxNumDistinctProteins=pd.NamedAgg(column="NumDistinctProteins", aggfunc="max"),
    )
    .reset_index()
)

stats_of_max_per_fraction_per_data_creation_merged_distinct_df

# %% [markdown]
# ## Absolute results

# %%
# dashes = ["solid", "dash", "dot"]
# dashes = ["dash", "dot"]
data_types = ["Complete", "Errored+PartiallyUnknown"]
# colors = ["red", "blue", "green"]
colors = ["black", "grey"]


# unknown_probabilities_dash_dict = {unknown_probability: color for unknown_probability, color in zip(unknown_probabilities, ["red", "blue", "green"])}
# chrom_color_dict = {chrom: dash for chrom, dash in zip(swiss_prot_name_by_chrom, ['dash', 'dot', 'solid'])}

data_types_color_dict = {
    data_type: color for data_type, color in zip(data_types, colors)
}
# chrom_color_dict = {
#     chrom: color for chrom, color in zip(swiss_prot_name_by_chrom, colors)
# }

rows = len(unknown_probabilities)
cols = len(swiss_prot_name_by_chrom)

fig = make_subplots(
    rows=rows,
    cols=cols,
    print_grid=False,
    # row_titles=unknown_probabilities,
    row_titles=[
        f"P[NA]={unknown_probability}" for unknown_probability in unknown_probabilities
    ],
    column_titles=short_swiss_prot_names,
    x_title="Simulated reads",
    y_title="Distinct proteins",
    shared_xaxes="all",
    shared_yaxes="all",
)

for row, unknown_probability in enumerate(unknown_probabilities, start=1):
    for col, (chrom, swiss_prot_name) in enumerate(
        swiss_prot_name_by_chrom.items(), start=1
    ):
        df = stats_of_max_per_fraction_per_data_creation_merged_distinct_df.loc[
            (
                stats_of_max_per_fraction_per_data_creation_merged_distinct_df["Chrom"]
                == chrom
            )
            & (
                stats_of_max_per_fraction_per_data_creation_merged_distinct_df[
                    "UnknownProbability"
                ]
                == unknown_probability
            )
        ]
        for data_type in data_types:
            data_type_df = df.loc[df["DataType"] == data_type]

            x = data_type_df["NumOfReads"]
            y = data_type_df["MeanOfMaxDistinctProteins"]
            error_y = data_type_df["StdOfMaxDistinctProteins"]

            color = data_types_color_dict[data_type]
            # dash = unknown_probabilities_dash_dict[unknown_probability]
            # dash = data_types_dash_dict[data_type]
            swiss_prot_name = swiss_prot_name.split("_")[0]
            # ic(color, dash)
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=y,
                    error_y=dict(
                        type="data",  # value of error bar given in data coordinates
                        array=error_y,
                        visible=True,
                    ),
                    # legendgroup=swiss_prot_name,  # this can be any string
                    # legendgrouptitle_text=swiss_prot_name,
                    # name=unknown_probability,
                    # legendgroup=unknown_probability,  # this can be any string
                    # legendgrouptitle_text=unknown_probability,
                    # name=swiss_prot_name,
                    showlegend=False,
                    line=dict(
                        color=color,
                        # dash=dash
                    ),
                    opacity=0.7,
                ),
                row=row,
                col=col,
            )

# Custom Legend
for data_type, color in data_types_color_dict.items():
    if data_type == "Errored+PartiallyUnknown":
        data_type = "Errored + partially NA      "
    fig.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            # mode="lines+markers",
            mode="lines",
            name=data_type,
            marker=dict(size=7, color=color),
        )
    )

fig.update_xaxes(
    # title_standoff = 50
    ticksuffix="    "
)

width = 900
# height = width * 600 / 1000
height = 600

# width = 1000
# # height = width * 600 / 1000
# height = 650

fig.update_layout(
    width=width,
    height=height,
    legend=dict(
        # title="Swiss prot name, unknown probability",
        # title="Unknown probability, gene",
        # title="Data type",
        # title_side="top",
        orientation="h",
        yanchor="bottom",
        y=1.1,
        xanchor="right",
        # x=0.9,
        x=1,
        itemwidth=30,
    ),
)

fig.write_image(
    "Graph assessment - distinct proteins recovered vs num of reads.svg",
    width=width,
    height=height,
)

fig.show()

# %%
# fig = px.line(
#     max_per_fraction_merged_distinct_df,
#     x="NumOfReads",
#     y="MaxNumDistinctProteins",
#     facet_col="SwissProtName",
#     facet_row="UnknownProbability",
#     color="DataType",
#     markers=True,
# )
# fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
# # fig.update_xaxes(range=[0, merged_distinct_df["NumOfReads"].max()*1.1])
# fig.update_layout(width=1000, height=700)
# fig.show()

# %%
# fig = px.scatter(
#     stats_of_merged_distinct_df,
#     x="NumOfReads",
#     y="MeanNumDistinctProteins",
#     error_y="StdNumDistinctProteins",
#     facet_col="SwissProtName",
#     facet_row="UnknownProbability",
#     color="DataType",
#     # markers=True,
# )
# # fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
# fig.for_each_annotation(
#     lambda a: a.update(
#         # text=a.text.split("=")[-1] if "UnknownProbability" not in a.text else f"UP = {a.text.split('=')[-1]}"
#         text=a.text.split("=")[-1]
#         if "UnknownProbability" not in a.text
#         else f"P[NA] = {a.text.split('=')[-1]}"
#     )
# )
# # fig.update_xaxes(range=[0, merged_distinct_df["NumOfReads"].max()*1.1])
# fig.update_layout(width=1000, height=700)
# fig.show()

# %%
# fig = px.box(
#     merged_distinct_df,
#     x="NumOfReads",
#     y="NumDistinctProteins",
#     facet_col="SwissProtName",
#     facet_row="UnknownProbability",
#     color="DataType",
#     # markers=True,
# )
# # fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
# fig.for_each_annotation(
#     lambda a: a.update(
#         # text=a.text.split("=")[-1] if "UnknownProbability" not in a.text else f"UP = {a.text.split('=')[-1]}"
#         text=a.text.split("=")[-1]
#         if "UnknownProbability" not in a.text
#         else f"P[NA] = {a.text.split('=')[-1]}"
#     )
# )
# # fig.update_xaxes(range=[0, merged_distinct_df["NumOfReads"].max()*1.1])
# fig.update_layout(width=1000, height=600)
# fig.show()

# %% [markdown]
# ## Relative % of recovered isoforms

# %%
max_per_fraction_per_data_creation_merged_distinct_df

# %%
complete_max_per_fraction_per_data_creation_merged_distinct_df = (
    max_per_fraction_per_data_creation_merged_distinct_df.loc[
        max_per_fraction_per_data_creation_merged_distinct_df["DataType"] == "Complete",
        [
            "Chrom",
            "SwissProtName",  # same as Chrom
            "UnknownProbability",
            "DataCreationRepetition",
            "DataType",
            "Fraction",
            "NumOfReads",  # same as Fraction
            "MaxNumDistinctProteins",
        ],
    ]
    .drop(columns=["DataType"])
    .rename(columns={"MaxNumDistinctProteins": "CompleteMaxNumOfDistinctProteins"})
)

partially_unknown_max_per_fraction_per_data_creation_merged_distinct_df = (
    max_per_fraction_per_data_creation_merged_distinct_df.loc[
        max_per_fraction_per_data_creation_merged_distinct_df["DataType"] != "Complete",
        [
            "Chrom",
            "SwissProtName",
            "UnknownProbability",
            "DataCreationRepetition",
            "DataType",
            "Fraction",
            "NumOfReads",
            "MaxNumDistinctProteins",
        ],
    ]
    .drop(columns=["DataType"])
    .rename(
        columns={"MaxNumDistinctProteins": "PartiallyUnknownMaxNumOfDistinctProteins"}
    )
)

wide_max_per_fraction_per_data_creation_merged_distinct_df = (
    complete_max_per_fraction_per_data_creation_merged_distinct_df.merge(
        partially_unknown_max_per_fraction_per_data_creation_merged_distinct_df,
    )
)

wide_max_per_fraction_per_data_creation_merged_distinct_df[
    "%DistinctProteinsRecovered"
] = (
    wide_max_per_fraction_per_data_creation_merged_distinct_df[
        "PartiallyUnknownMaxNumOfDistinctProteins"
    ]
    .mul(100)
    .div(
        wide_max_per_fraction_per_data_creation_merged_distinct_df[
            "CompleteMaxNumOfDistinctProteins"
        ]
    )
    .round(2)
)

# wide_max_per_fraction_merged_distinct_df["SwissProtName - UnknownProbability"] = wide_max_per_fraction_merged_distinct_df["SwissProtName"] + " - " + wide_max_per_fraction_merged_distinct_df["UnknownProbability"].astype(str)

wide_max_per_fraction_per_data_creation_merged_distinct_df = (
    wide_max_per_fraction_per_data_creation_merged_distinct_df.sort_values(
        [
            "Chrom",
            "SwissProtName",
            "UnknownProbability",
            "DataCreationRepetition",
            "Fraction",
            "NumOfReads",
        ]
    )
).reset_index(drop=True)

wide_max_per_fraction_per_data_creation_merged_distinct_df
# wide_max_per_fraction_per_data_creation_merged_distinct_df.head()

# %%
# wide_max_per_fraction_per_data_creation_merged_distinct_df.head(20)

# %%
# stats over data creation repetitions

stats_of_max_per_fraction_merged_distinct_df = (
    wide_max_per_fraction_per_data_creation_merged_distinct_df.groupby(
        [
            "Chrom",
            "SwissProtName",  # same as Chrom
            "UnknownProbability",
            "Fraction",
            "NumOfReads",  # same as Fraction
        ]
    )
    .agg(
        MeanPrctDistinctProteinsRecovered=pd.NamedAgg(
            column="%DistinctProteinsRecovered", aggfunc="mean"
        ),
        StdPrctDistinctProteinsRecovered=pd.NamedAgg(
            column="%DistinctProteinsRecovered", aggfunc="std"
        ),
    )
    .reset_index()
    # .rename(columns={"MeanPrctDistinctProteinsRecovered": "MeanPrctDistinctProteinsRecovered"})
)
stats_of_max_per_fraction_merged_distinct_df

# %%
# unknown_probabilities

# %%
# swiss_prot_name_by_chrom

# %%
# dashes = ['dash', 'dot', 'solid']
# colors = ["red", "blue", "green"]

# # unknown_probabilities_dash_dict = {unknown_probability: color for unknown_probability, color in zip(unknown_probabilities, ["red", "blue", "green"])}
# # chrom_color_dict = {chrom: dash for chrom, dash in zip(swiss_prot_name_by_chrom, ['dash', 'dot', 'solid'])}

# unknown_probabilities_dash_dict = {unknown_probability: dash for unknown_probability, dash in zip(unknown_probabilities, dashes)}
# chrom_color_dict = {chrom: color for chrom, color in zip(swiss_prot_name_by_chrom, colors)}

# # unknown_probabilities_dash_dict

# %%
dashes = ["solid", "dash", "dot"]
colors = ["red", "blue", "green"]

# unknown_probabilities_dash_dict = {unknown_probability: color for unknown_probability, color in zip(unknown_probabilities, ["red", "blue", "green"])}
# chrom_color_dict = {chrom: dash for chrom, dash in zip(swiss_prot_name_by_chrom, ['dash', 'dot', 'solid'])}

unknown_probabilities_dash_dict = {
    unknown_probability: dash
    for unknown_probability, dash in zip(unknown_probabilities, dashes)
}
chrom_color_dict = {
    chrom: color for chrom, color in zip(swiss_prot_name_by_chrom, colors)
}

fig = make_subplots(
    rows=1,
    cols=1,
    print_grid=False,
    x_title="Simulated reads",
    y_title="Mean % of distinct proteins recovered",
)

for unknown_probability in unknown_probabilities:
    for chrom, swiss_prot_name in swiss_prot_name_by_chrom.items():
        df = stats_of_max_per_fraction_merged_distinct_df.loc[
            (stats_of_max_per_fraction_merged_distinct_df["Chrom"] == chrom)
            & (
                stats_of_max_per_fraction_merged_distinct_df["UnknownProbability"]
                == unknown_probability
            )
        ]
        x = df["NumOfReads"]
        y = df["MeanPrctDistinctProteinsRecovered"]
        error_y = df["StdPrctDistinctProteinsRecovered"]

        color = chrom_color_dict[chrom]
        dash = unknown_probabilities_dash_dict[unknown_probability]
        swiss_prot_name = swiss_prot_name.split("_")[0]
        # ic(color, dash)
        fig.add_trace(
            go.Scatter(
                x=x,
                y=y,
                error_y=dict(
                    type="data",  # value of error bar given in data coordinates
                    array=error_y,
                    visible=True,
                ),
                # legendgroup=swiss_prot_name,  # this can be any string
                # legendgrouptitle_text=swiss_prot_name,
                # name=unknown_probability,
                legendgroup=unknown_probability,  # this can be any string
                legendgrouptitle_text=unknown_probability,
                name=swiss_prot_name,
                line=dict(color=color, dash=dash),
                opacity=0.7,
            )
        )

width = 700
# height = width * 600 / 1000
height = 550

fig.update_layout(
    width=width,
    height=height,
    legend=dict(
        # title="Swiss prot name, unknown probability",
        # title="Unknown probability, gene",
        # title="Unknown probability, mock gene",
        title="NA probability, mock gene       ",
    ),
    margin_r=140,
)

fig.write_image(
    "Graph assessment - mean % of distinct proteins recovered.svg",
    width=width,
    height=height,
)

fig.show()

# %%

# %% [markdown]
# # False positives

# %%
# !zcat /private7/projects/Combinatorics/Simulations/GraphAssessment/comp140826_c0_seq1.VPS8_HUMAN.UP0_09.Rep2.FalsePositives.tsv.gz | head

# %%
false_positives_files_fofn

# %%
with open(false_positives_files_fofn, "r") as f:
    putative_false_positives_files = f.read().split(" ")
ic(len(putative_false_positives_files))
# putative_false_positives_files[:3]

false_positives_files = []
unexisting_putative_files = 0
for putative_false_positives_file in putative_false_positives_files:
    putative_false_positives_file = Path(putative_false_positives_file)
    if putative_false_positives_file.exists():
        false_positives_files.append(putative_false_positives_file)
    else:
        warnings.warn(
            f"{putative_false_positives_file} doesn't exist - maybe just yet?"
        )
        unexisting_putative_files += 1
ic(len(false_positives_files))
ic(unexisting_putative_files)
false_positives_files[:3]

# %%
# false_positivs_dfs_dict = {}

false_positives_dfs = []

for false_positives_file in false_positives_files:
    # file_name = re.sub(rf"{beginning_of_postfix}*.csv", "", distinct_file.name)
    file_name_parts = false_positives_file.name.split(".")
    chrom = file_name_parts[0]
    swiss_prot_name = file_name_parts[1].split("_")[0]
    unknown_probability = float(file_name_parts[2].replace("UP", "").replace("_", "."))
    repetition = file_name_parts[3]

    # false_positivs_dfs_dict[(chrom, swiss_prot_name, unknown_probability, repetition)] = pd.read_table()

    false_positives_df = pd.read_table(false_positives_file)
    false_positives_df["IndistinguishableCompleteProteins"] = false_positives_df[
        "IndistinguishableCompleteProteins"
    ].str.split(",")
    false_positives_df.insert(0, "Chrom", chrom)
    false_positives_df.insert(1, "SwissProtName", swiss_prot_name)
    false_positives_df.insert(2, "UnknownProbability", unknown_probability)
    false_positives_df.insert(3, "DataCreationRepetition", repetition)

    false_positives_dfs.append(false_positives_df)

concat_false_positives_df = pd.concat(false_positives_dfs, ignore_index=True)
concat_false_positives_df

# %%
# concat_false_positives_df.isna().any()

# %%
# concat_false_positives_df.loc[
#     concat_false_positives_df["IndistinguishableCompleteProteins"].isna()
# ]

# %%
# Out[132]["NumOfIndistinguishableCompleteProteins"].value_counts(dropna=False)

# %%
# concat_false_positives_df.loc[
#     concat_false_positives_df["NumOfIndistinguishableCompleteProteins"] == 0,
# ]

# %%
# Out[135]["IndistinguishableCompleteProteins"].value_counts(dropna=False)

# %%
# false_positives_df["IndistinguishableCompleteProteins"].replace(np.nan, [])

# %%
# next I should use the false-positives and this df to test for false-positives only inside the largest MIS of each dataset

max_per_fraction_per_data_creation_merged_distinct_df

# %%
chrom = "comp140826_c0_seq1"
swiss_prot_name = "VPS8"
unknown_prob = 0.09
repetition = "Rep1"
# fraction = 0.2
fraction = 0.6
# fraction = 1

# %%
# chrom = "comp141434_c0_seq1"
# swiss_prot_name = "VPS8"
# unknown_prob = 0.09
# repetition = "Rep1"
# fraction = 0.2
# # fraction = 0.6
# # fraction = 1

# chrom: ''
#     swiss_prot_name: 'CBPC1_HUMAN'
#     unknown_probability: 0.09
#     repetition: 'Rep4'
#     fraction: 0.2

# %%
one_fraction_and_data_creation_max_distinct_df = (
    max_per_fraction_per_data_creation_merged_distinct_df.loc[
        (max_per_fraction_per_data_creation_merged_distinct_df["Chrom"] == chrom)
        & (
            max_per_fraction_per_data_creation_merged_distinct_df["SwissProtName"]
            == swiss_prot_name
        )
        & (
            max_per_fraction_per_data_creation_merged_distinct_df["UnknownProbability"]
            == unknown_prob
        )
        & (
            max_per_fraction_per_data_creation_merged_distinct_df[
                "DataCreationRepetition"
            ]
            == repetition
        )
        & (
            max_per_fraction_per_data_creation_merged_distinct_df["Fraction"]
            == fraction
        ),
        ["DataType", "DistinctProteins", "Fraction", "NumOfReads"],
    ].copy()
)
one_fraction_and_data_creation_max_distinct_df[
    "DistinctProteins"
] = one_fraction_and_data_creation_max_distinct_df["DistinctProteins"].str.split(",")
one_fraction_and_data_creation_max_distinct_df

# %%
num_of_reads = one_fraction_and_data_creation_max_distinct_df["NumOfReads"].iloc[0]
# num_of_reads = one_fraction_and_data_creation_max_distinct_df["NumOfReads"].iloc[0]
num_of_reads

# %%
# one_fraction_and_data_creation_max_distinct_df.loc[
#     one_fraction_and_data_creation_max_distinct_df["DataType"] == "Complete",
#     "DistinctProteins"
# ].values[0]

# %%
complete_distinct_proteins = one_fraction_and_data_creation_max_distinct_df.loc[
    one_fraction_and_data_creation_max_distinct_df["DataType"] == "Complete",
    "DistinctProteins",
].values[0]
errored_na_distinct_proteins = one_fraction_and_data_creation_max_distinct_df.loc[
    one_fraction_and_data_creation_max_distinct_df["DataType"]
    == "Errored+PartiallyUnknown",
    "DistinctProteins",
].values[0]
complete_distinct_proteins = set(complete_distinct_proteins)
errored_na_distinct_proteins = set(errored_na_distinct_proteins)
ic(len(complete_distinct_proteins), len(errored_na_distinct_proteins));

# %%
one_fraction_and_data_creation_false_positives_df = concat_false_positives_df.loc[
    (concat_false_positives_df["Chrom"] == chrom)
    & (concat_false_positives_df["SwissProtName"] == swiss_prot_name)
    & (concat_false_positives_df["UnknownProbability"] == unknown_prob)
    & (concat_false_positives_df["DataCreationRepetition"] == repetition)
    & (
        concat_false_positives_df["Errored+PartiallyUnknownProtein"].isin(
            errored_na_distinct_proteins
        )
    ),
].copy()

one_fraction_and_data_creation_false_positives_df.insert(
    one_fraction_and_data_creation_false_positives_df.columns.get_loc(
        "DataCreationRepetition"
    )
    + 1,
    "Fraction",
    fraction,
)
one_fraction_and_data_creation_false_positives_df.insert(
    one_fraction_and_data_creation_false_positives_df.columns.get_loc(
        "DataCreationRepetition"
    )
    + 2,
    "NumOfReads",
    num_of_reads,
)

one_fraction_and_data_creation_false_positives_df

# %%
is_false_positive = one_fraction_and_data_creation_false_positives_df.apply(
    lambda x: (x["NumOfIndistinguishableCompleteProteins"] == 0)
    or (
        x["NumOfIndistinguishableCompleteProteins"] > 0
        and set(x["IndistinguishableCompleteProteins"]) & complete_distinct_proteins
        == set()
    ),
    axis=1,
)

ic(100 * is_false_positive.sum() / len(is_false_positive))

is_false_positive

# %%
# validated_false_positives_df = one_fraction_and_data_creation_false_positives_df.loc[
#     is_false_positive
# ]
# validated_false_positives_df

# %%
# one_fraction_and_data_creation_false_positives_df.loc[
#     :, "IsFalsePositive"
# ] = is_false_positive.tolist()
# one_fraction_and_data_creation_false_positives_df

# %%
one_fraction_and_data_creation_false_positives_df[
    "IsFalsePositive"
] = is_false_positive.tolist()
one_fraction_and_data_creation_false_positives_df

# %%
one_fraction_and_data_creation_false_positives_df.loc[
    (one_fraction_and_data_creation_false_positives_df["IsFalsePositive"])
    & (
        one_fraction_and_data_creation_false_positives_df[
            "NumOfIndistinguishableCompleteProteins"
        ]
        > 1
    )
].head(100)

# %%
"5cX" in complete_distinct_proteins

# %%
"axF" in complete_distinct_proteins

# %%
"6Rp" in complete_distinct_proteins

# %%
"brB" in complete_distinct_proteins

# %%
# "ahs" in complete_distinct_proteins

# %%
# product(swiss_prot_name_by_chrom.items(), unknown_probabilities, repetitions, fractions)

# %%

# %%
validated_false_positives_dfs = []

for (chrom, swiss_prot_name), unknown_probability, repetition, fraction in product(
    swiss_prot_name_by_chrom.items(), unknown_probabilities, repetitions, fractions
):
    # ic(chrom, swiss_prot_name, unknown_probability, repetition, fraction)

    swiss_prot_name = swiss_prot_name.split("_")[0]

    one_fraction_and_data_creation_max_distinct_df = (
        max_per_fraction_per_data_creation_merged_distinct_df.loc[
            (max_per_fraction_per_data_creation_merged_distinct_df["Chrom"] == chrom)
            & (
                max_per_fraction_per_data_creation_merged_distinct_df["SwissProtName"]
                == swiss_prot_name
            )
            & (
                max_per_fraction_per_data_creation_merged_distinct_df[
                    "UnknownProbability"
                ]
                == unknown_probability
            )
            & (
                max_per_fraction_per_data_creation_merged_distinct_df[
                    "DataCreationRepetition"
                ]
                == repetition
            )
            & (
                max_per_fraction_per_data_creation_merged_distinct_df["Fraction"]
                == fraction
            ),
            ["DataType", "DistinctProteins", "Fraction", "NumOfReads"],
        ].copy()
    )
    one_fraction_and_data_creation_max_distinct_df[
        "DistinctProteins"
    ] = one_fraction_and_data_creation_max_distinct_df["DistinctProteins"].str.split(
        ","
    )

    num_of_reads = one_fraction_and_data_creation_max_distinct_df["NumOfReads"].iloc[0]

    complete_distinct_proteins = one_fraction_and_data_creation_max_distinct_df.loc[
        one_fraction_and_data_creation_max_distinct_df["DataType"] == "Complete",
        "DistinctProteins",
    ].values[0]
    errored_na_distinct_proteins = one_fraction_and_data_creation_max_distinct_df.loc[
        one_fraction_and_data_creation_max_distinct_df["DataType"]
        == "Errored+PartiallyUnknown",
        "DistinctProteins",
    ].values[0]
    complete_distinct_proteins = set(complete_distinct_proteins)
    errored_na_distinct_proteins = set(errored_na_distinct_proteins)
    # ic(len(complete_distinct_proteins), len(errored_na_distinct_proteins));

    one_fraction_and_data_creation_false_positives_df = concat_false_positives_df.loc[
        (concat_false_positives_df["Chrom"] == chrom)
        & (concat_false_positives_df["SwissProtName"] == swiss_prot_name)
        & (concat_false_positives_df["UnknownProbability"] == unknown_probability)
        & (concat_false_positives_df["DataCreationRepetition"] == repetition)
        & (
            concat_false_positives_df["Errored+PartiallyUnknownProtein"].isin(
                errored_na_distinct_proteins
            )
        ),
    ].copy()
    one_fraction_and_data_creation_false_positives_df.insert(
        one_fraction_and_data_creation_false_positives_df.columns.get_loc(
            "DataCreationRepetition"
        )
        + 1,
        "Fraction",
        fraction,
    )
    one_fraction_and_data_creation_false_positives_df.insert(
        one_fraction_and_data_creation_false_positives_df.columns.get_loc(
            "DataCreationRepetition"
        )
        + 2,
        "NumOfReads",
        num_of_reads,
    )

    is_false_positive = one_fraction_and_data_creation_false_positives_df.apply(
        lambda x: (x["NumOfIndistinguishableCompleteProteins"] == 0)
        or (
            x["NumOfIndistinguishableCompleteProteins"] > 0
            and set(x["IndistinguishableCompleteProteins"]) & complete_distinct_proteins
            == set()
        ),
        axis=1,
    )

    # validated_false_positives_df = one_fraction_and_data_creation_false_positives_df.loc[is_false_positive]
    # validated_false_positives_dfs.append(validated_false_positives_df)

    one_fraction_and_data_creation_false_positives_df.loc[
        :, "IsFalsePositive"
    ] = is_false_positive
    validated_false_positives_dfs.append(
        one_fraction_and_data_creation_false_positives_df
    )

concat_validated_false_positives_df = pd.concat(
    validated_false_positives_dfs, ignore_index=True
)
concat_validated_false_positives_df

# %%
concat_validated_false_positives_prct_df = (
    concat_validated_false_positives_df.groupby(
        [
            "Chrom",
            "SwissProtName",
            "UnknownProbability",
            "Fraction",
            "NumOfReads",
            "DataCreationRepetition",
        ]
    )["IsFalsePositive"]
    .apply(lambda x: 100 * sum(x) / len(x))
    .reset_index()
    .rename(columns={"IsFalsePositive": "%FalsePositivse"})
    .groupby(
        [
            "Chrom",
            "SwissProtName",
            "UnknownProbability",
            "Fraction",
            "NumOfReads",
        ]
    )["%FalsePositivse"]
    .agg(["mean", "std"])
    .reset_index()
    .rename(
        columns={
            "mean": "MeanPrctOfPalsePositives",
            "std": "STDOfMeanPrctOfPalsePositives",
        }
    )
)


concat_validated_false_positives_prct_df

# %%
concat_validated_false_positives_prct_df["STDOfMeanPrctOfPalsePositives"].describe()

# %%
dashes = ["solid", "dash", "dot"]
colors = ["red", "blue", "green"]

# unknown_probabilities_dash_dict = {unknown_probability: color for unknown_probability, color in zip(unknown_probabilities, ["red", "blue", "green"])}
# chrom_color_dict = {chrom: dash for chrom, dash in zip(swiss_prot_name_by_chrom, ['dash', 'dot', 'solid'])}

unknown_probabilities_dash_dict = {
    unknown_probability: dash
    for unknown_probability, dash in zip(unknown_probabilities, dashes)
}
chrom_color_dict = {
    chrom: color for chrom, color in zip(swiss_prot_name_by_chrom, colors)
}

fig = make_subplots(
    rows=1,
    cols=1,
    print_grid=False,
    x_title="Simulated reads",
    y_title="Mean % of false-positive distinct proteins",
)

for unknown_probability in unknown_probabilities:
    for chrom, swiss_prot_name in swiss_prot_name_by_chrom.items():
        df = concat_validated_false_positives_prct_df.loc[
            (concat_validated_false_positives_prct_df["Chrom"] == chrom)
            & (
                concat_validated_false_positives_prct_df["UnknownProbability"]
                == unknown_probability
            )
        ]
        x = df["NumOfReads"]
        y = df["MeanPrctOfPalsePositives"]
        error_y = df["STDOfMeanPrctOfPalsePositives"]

        color = chrom_color_dict[chrom]
        dash = unknown_probabilities_dash_dict[unknown_probability]
        swiss_prot_name = swiss_prot_name.split("_")[0]
        # ic(color, dash)
        fig.add_trace(
            go.Scatter(
                x=x,
                y=y,
                error_y=dict(
                    type="data",  # value of error bar given in data coordinates
                    array=error_y,
                    visible=True,
                ),
                # legendgroup=swiss_prot_name,  # this can be any string
                # legendgrouptitle_text=swiss_prot_name,
                # name=unknown_probability,
                legendgroup=unknown_probability,  # this can be any string
                legendgrouptitle_text=unknown_probability,
                name=swiss_prot_name,
                line=dict(color=color, dash=dash),
                opacity=0.7,
            )
        )

width = 700
# height = width * 600 / 1000
height = 550

fig.update_yaxes(
    range=[
        0,
        concat_validated_false_positives_prct_df["MeanPrctOfPalsePositives"].max()
        * 1.05,
    ]
)

fig.update_layout(
    width=width,
    height=height,
    legend=dict(
        # title="Swiss prot name, unknown probability",
        # title="Unknown probability, gene",
        # title="Unknown probability, mock gene",
        title="NA probability, mock gene       ",
    ),
)

fig.write_image(
    "Graph assessment - mean % of false-positive distinct proteins.svg",
    width=width,
    height=height,
)

fig.show()

# %%
swiss_prot_name_by_chrom

# %%
# !cat /private7/projects/Combinatorics/D.pealeii/Annotations/D.pea.EditingSites.bed | grep CBPC1_HUMAN | wc -l

# %%
# !cat /private7/projects/Combinatorics/D.pealeii/Annotations/D.pea.EditingSites.bed | grep MY18A_HUMAN | wc -l

# %%
# !cat /private7/projects/Combinatorics/D.pealeii/Annotations/D.pea.EditingSites.bed | grep VPS8_HUMAN | wc -l

# %%
# # !cat /private7/projects/Combinatorics/D.pealeii/Annotations/D.pea.EditingSites.bed | grep comp141434_c0_seq1 | wc -l

# %%
# # !cat /private7/projects/Combinatorics/D.pealeii/Annotations/D.pea.EditingSites.bed | grep comp141359_c0_seq2 | wc -l

# %%
# # !cat /private7/projects/Combinatorics/D.pealeii/Annotations/D.pea.EditingSites.bed | grep comp140826_c0_seq1 | wc -l

# %%

# %%

# %%

# %%

# %% [markdown]
# ## False positives QC

# %% [markdown]
# ### FP in fraction 1.0

# %%
# chrom = "comp141434_c0_seq1"
# unknown_prob = 0.09
# repetition = "Rep1"
# fraction = 1.0

# %%
# false_positives_df

# %%
concat_validated_false_positives_df

# %%
# concat_validated_false_positives_df["IsFalsePositive"].value_counts(dropna=False)

# %%
f_1_concat_validated_false_positives_df = concat_validated_false_positives_df.loc[
    (concat_validated_false_positives_df["Fraction"] == 1.0)
    & (concat_validated_false_positives_df["IsFalsePositive"]),
]
f_1_concat_validated_false_positives_df

# %%
f_1_concat_validated_false_positives_df[
    "NumOfIndistinguishableCompleteProteins"
].value_counts()

# %% [markdown]
# It seems that each protein that is FP at fraction 1.0 is distinct from each and every other complete protein.

# %%
f_1_concat_validated_false_positives_df.groupby(
    ["Chrom", "SwissProtName", "UnknownProbability", "DataCreationRepetition"]
).size()

# %% [markdown]
# ### FP in fraction 0.2

# %%
chrom = "comp141434_c0_seq1"
unknown_prob = 0.09
repetition = "Rep1"
fraction = 0.2

# %%
false_positives_df

# %%
concat_validated_false_positives_df

# %%
concat_validated_false_positives_df["IsFalsePositive"].value_counts(dropna=False)

# %%
fp_02_sample_df = concat_validated_false_positives_df.loc[
    (concat_validated_false_positives_df["Chrom"] == chrom)
    & (concat_validated_false_positives_df["UnknownProbability"] == unknown_prob)
    & (concat_validated_false_positives_df["DataCreationRepetition"] == repetition)
    & (concat_validated_false_positives_df["Fraction"] == fraction)
    & (concat_validated_false_positives_df["IsFalsePositive"])
    & (
        concat_validated_false_positives_df["NumOfIndistinguishableCompleteProteins"]
        > 0
    )
    # & (concat_validated_false_positives_df[""] == )
    # ].sample(3, random_state=1892)
].sample(3, random_state=193)

fp_02_sample_df

# %%
fp_02_proteins = fp_02_sample_df["Errored+PartiallyUnknownProtein"].tolist()
fp_02_proteins

# %%
neighbours_of_fp_02_proteins = [
    c_protein
    for c_proteins in fp_02_sample_df["IndistinguishableCompleteProteins"]
    .dropna()
    .tolist()
    for c_protein in c_proteins
]
neighbours_of_fp_02_proteins

# %%
fp_qc_data_types = ["Errored+PartiallyUnknown"] * len(fp_02_proteins) + [
    "Complete"
] * len(neighbours_of_fp_02_proteins)
fp_qc_data_types

# %%
fp_02_proteins + neighbours_of_fp_02_proteins

# %%
distinct_in_fraction_dict = defaultdict(list)

for protein, data_type in zip(
    fp_02_proteins + neighbours_of_fp_02_proteins, fp_qc_data_types
):
    ic(protein, data_type)

    for fraction in fractions:
        one_fraction_and_data_creation_max_distinct_df = max_per_fraction_per_data_creation_merged_distinct_df.loc[
            (max_per_fraction_per_data_creation_merged_distinct_df["Chrom"] == chrom)
            & (
                max_per_fraction_per_data_creation_merged_distinct_df[
                    "UnknownProbability"
                ]
                == unknown_prob
            )
            & (
                max_per_fraction_per_data_creation_merged_distinct_df[
                    "DataCreationRepetition"
                ]
                == repetition
            )
            & (
                max_per_fraction_per_data_creation_merged_distinct_df["Fraction"]
                == fraction
            ),
            # ["DataType", "DistinctProteins", "Fraction", "NumOfReads"],
        ].copy()
        one_fraction_and_data_creation_max_distinct_df[
            "DistinctProteins"
        ] = one_fraction_and_data_creation_max_distinct_df[
            "DistinctProteins"
        ].str.split(
            ","
        )

        distinct_in_fraction = (
            protein
            in one_fraction_and_data_creation_max_distinct_df.loc[
                one_fraction_and_data_creation_max_distinct_df["DataType"] == data_type,
                "DistinctProteins",
            ].values[0]
        )
        distinct_in_fraction_dict[protein].append(distinct_in_fraction)
        # distinct_in_fraction_dict[protein].append(distinct_in_fraction)
    # break

# %%
distinct_in_fraction_df = (
    pd.DataFrame(distinct_in_fraction_dict)
    .T.set_axis(fractions, axis="columns")
    .reset_index(names="Protein")
)
distinct_in_fraction_df.insert(0, "DataType", fp_qc_data_types)
distinct_in_fraction_df

# %%

# %%
fp_in_fraction_dict = defaultdict(list)

for protein in fp_02_proteins:
    ic(protein)

    for fraction in fractions:
        fp_in_fraction = (
            protein
            in concat_validated_false_positives_df.loc[
                (concat_validated_false_positives_df["Chrom"] == chrom)
                & (
                    concat_validated_false_positives_df["UnknownProbability"]
                    == unknown_prob
                )
                & (
                    concat_validated_false_positives_df["DataCreationRepetition"]
                    == repetition
                )
                & (concat_validated_false_positives_df["Fraction"] == fraction)
                & (concat_validated_false_positives_df["IsFalsePositive"]),
                "Errored+PartiallyUnknownProtein",
            ].tolist()
        )

        fp_in_fraction_dict[protein].append(fp_in_fraction)
        # distinct_in_fraction_dict[protein].append(distinct_in_fraction)
    # break

# %%
fp_02_sample_df.loc[
    :, ["Errored+PartiallyUnknownProtein", "IndistinguishableCompleteProteins"]
]

# %%
fp_in_fraction_df = (
    pd.DataFrame(fp_in_fraction_dict)
    .T.set_axis(fractions, axis="columns")
    .reset_index(names="Protein")
)
fp_in_fraction_df

# %%
distinct_in_fraction_df

# %% [markdown]
# #### kqT 	[2IR, 2Qv]

# %%
distinct_in_fraction_df.iloc[[0, 3, 4]]

# %%
fp_in_fraction_df.iloc[[0]]

# %% [markdown]
# #### 9i9 	[2Oq]

# %%
distinct_in_fraction_df.iloc[[1, 5]]

# %%
fp_in_fraction_df.iloc[[1]]

# %% [markdown]
# #### 8Db 	[34w]

# %%
distinct_in_fraction_df.iloc[[2, 6]]

# %%
fp_in_fraction_df.iloc[[2]]

# %% [markdown]
# # False negatives

# %%
max_per_fraction_per_data_creation_merged_distinct_df

# %%
concat_validated_false_positives_df

# %%
concat_validated_false_positives_df.drop_duplicates(
    ["Chrom", "SwissProtName", "UnknownProbability", "Errored+PartiallyUnknownProtein"]
)["NumOfIndistinguishableCompleteProteins"].describe()

# %%
concat_validated_false_positives_df.loc[
    ~concat_validated_false_positives_df["IsFalsePositive"]
].drop_duplicates(
    ["Chrom", "SwissProtName", "UnknownProbability", "Errored+PartiallyUnknownProtein"]
)[
    "NumOfIndistinguishableCompleteProteins"
].describe()

# %%
concat_validated_false_positives_df.loc[
    ~concat_validated_false_positives_df["IsFalsePositive"]
].drop_duplicates(
    [
        "Chrom",
        "SwissProtName",
        "UnknownProbability",
        "Fraction",
        "Errored+PartiallyUnknownProtein",
    ]
).groupby(
    ["Fraction"]
)[
    "NumOfIndistinguishableCompleteProteins"
].describe()

# %%
chrom = "comp141434_c0_seq1"
unknown_prob = 0.09
repetition = "Rep1"
fraction = 0.2

# %%
# one_fraction_and_data_creation_max_distinct_df = max_per_fraction_per_data_creation_merged_distinct_df.loc[
#     (max_per_fraction_per_data_creation_merged_distinct_df["Chrom"] == chrom)
#     & (max_per_fraction_per_data_creation_merged_distinct_df["UnknownProbability"] == unknown_prob)
#     & (max_per_fraction_per_data_creation_merged_distinct_df["DataCreationRepetition"] == repetition)
#     & (max_per_fraction_per_data_creation_merged_distinct_df["Fraction"] == fraction)
#     # & (concat_validated_false_positives_df["IsFalsePositive"])
#     # & (concat_validated_false_positives_df["NumOfIndistinguishableCompleteProteins"] > 0)
# ]
# one_fraction_and_data_creation_max_distinct_df

# %%
one_fraction_and_data_creation_max_distinct_df = (
    max_per_fraction_per_data_creation_merged_distinct_df.loc[
        (max_per_fraction_per_data_creation_merged_distinct_df["Chrom"] == chrom)
        & (
            max_per_fraction_per_data_creation_merged_distinct_df["UnknownProbability"]
            == unknown_prob
        )
        & (
            max_per_fraction_per_data_creation_merged_distinct_df[
                "DataCreationRepetition"
            ]
            == repetition
        )
        & (
            max_per_fraction_per_data_creation_merged_distinct_df["Fraction"]
            == fraction
        ),
        [
            "DataType",
            "MaxNumDistinctProteins",
            "DistinctProteins",
            "Fraction",
            "NumOfReads",
        ],
    ].copy()
)
one_fraction_and_data_creation_max_distinct_df[
    "DistinctProteins"
] = one_fraction_and_data_creation_max_distinct_df["DistinctProteins"].str.split(",")
one_fraction_and_data_creation_max_distinct_df

# %%
one_fraction_and_data_creation_concat_validated_false_positives_df = concat_validated_false_positives_df.loc[
    (concat_validated_false_positives_df["Chrom"] == chrom)
    & (concat_validated_false_positives_df["UnknownProbability"] == unknown_prob)
    & (concat_validated_false_positives_df["DataCreationRepetition"] == repetition)
    & (concat_validated_false_positives_df["Fraction"] == fraction)
    & (~concat_validated_false_positives_df["IsFalsePositive"])
    # & (concat_validated_false_positives_df["NumOfIndistinguishableCompleteProteins"] > 0)
]
one_fraction_and_data_creation_concat_validated_false_positives_df

# %%
all_indistinguishable_complete_proteins_of_errored_distinct_proteins = [
    c_protein
    for c_proteins in one_fraction_and_data_creation_concat_validated_false_positives_df[
        "IndistinguishableCompleteProteins"
    ]
    .dropna()
    .tolist()
    for c_protein in c_proteins
]
ic(len(all_indistinguishable_complete_proteins_of_errored_distinct_proteins))

complete_distinct_proteins = one_fraction_and_data_creation_max_distinct_df.loc[
    one_fraction_and_data_creation_max_distinct_df["DataType"] == "Complete",
    "DistinctProteins",
].values[0]
ic(len(complete_distinct_proteins));

# %%
all_indistinguishable_complete_proteins_of_errored_distinct_proteins[:3]

# %%
uncovered_complete_distinct_proteins = [
    protein
    for protein in complete_distinct_proteins
    if protein
    not in all_indistinguishable_complete_proteins_of_errored_distinct_proteins
]

ic(len(uncovered_complete_distinct_proteins));

# %%
one_fraction_and_data_creation_concat_validated_false_positives_df.shape[0]

# %%
# chrom = "comp141434_c0_seq1"
# unknown_prob = 0.09
# repetition = "Rep1"
# fraction = 0.2

# row = FNRow(
#     chrom,
#     swiss_prot_name_by_chrom[chrom],
#     unknown_prob,
#     repetition,
#     fraction,
#     len(complete_distinct_proteins),
#     len(complete_distinct_proteins) - len(uncovered_complete_distinct_proteins),
#     len(uncovered_complete_distinct_proteins),
#     one_fraction_and_data_creation_concat_validated_false_positives_df.shape[0],
# )

# row

# %%
# pd.DataFrame([row, row])

# %%
swiss_prot_name_by_chrom


# %%
# fn_df_rows = []

# i = 1

# for (chrom, swiss_prot_name), unknown_probability, repetition, fraction in product(
#     swiss_prot_name_by_chrom.items(), unknown_probabilities, repetitions, fractions
# ):
#     swiss_prot_name = swiss_prot_name.split("_")[0]

#     ic(i, chrom, swiss_prot_name, unknown_probability, repetition, fraction)

#     i += 1

#     one_fraction_and_data_creation_max_distinct_df = (
#         max_per_fraction_per_data_creation_merged_distinct_df.loc[
#             (max_per_fraction_per_data_creation_merged_distinct_df["Chrom"] == chrom)
#             & (
#                 max_per_fraction_per_data_creation_merged_distinct_df[
#                     "UnknownProbability"
#                 ]
#                 == unknown_probability
#             )
#             & (
#                 max_per_fraction_per_data_creation_merged_distinct_df[
#                     "DataCreationRepetition"
#                 ]
#                 == repetition
#             )
#             & (
#                 max_per_fraction_per_data_creation_merged_distinct_df["Fraction"]
#                 == fraction
#             ),
#             [
#                 "DataType",
#                 "MaxNumDistinctProteins",
#                 "DistinctProteins",
#                 "Fraction",
#                 "NumOfReads",
#             ],
#         ].copy()
#     )
#     one_fraction_and_data_creation_max_distinct_df[
#         "DistinctProteins"
#     ] = one_fraction_and_data_creation_max_distinct_df["DistinctProteins"].str.split(
#         ","
#     )

#     num_of_reads = one_fraction_and_data_creation_max_distinct_df["NumOfReads"].iloc[0]

#     one_fraction_and_data_creation_concat_validated_false_positives_df = (
#         concat_validated_false_positives_df.loc[
#             (concat_validated_false_positives_df["Chrom"] == chrom)
#             & (
#                 concat_validated_false_positives_df["UnknownProbability"]
#                 == unknown_prob
#             )
#             & (
#                 concat_validated_false_positives_df["DataCreationRepetition"]
#                 == repetition
#             )
#             & (concat_validated_false_positives_df["Fraction"] == fraction)
#             & (~concat_validated_false_positives_df["IsFalsePositive"])
#         ]
#     )

#     all_indistinguishable_complete_proteins_of_errored_distinct_proteins = [
#         c_protein
#         for c_proteins in one_fraction_and_data_creation_concat_validated_false_positives_df[
#             "IndistinguishableCompleteProteins"
#         ]
#         .dropna()
#         .tolist()
#         for c_protein in c_proteins
#     ]
#     # ic(len(all_indistinguishable_complete_proteins_of_errored_distinct_proteins));

#     complete_distinct_proteins = one_fraction_and_data_creation_max_distinct_df.loc[
#         one_fraction_and_data_creation_max_distinct_df["DataType"] == "Complete",
#         "DistinctProteins",
#     ].values[0]

#     uncovered_complete_distinct_proteins = [
#         protein
#         for protein in complete_distinct_proteins
#         if protein
#         not in all_indistinguishable_complete_proteins_of_errored_distinct_proteins
#     ]

#     row = FNRow(
#         chrom,
#         swiss_prot_name,
#         unknown_probability,
#         repetition,
#         fraction,
#         num_of_reads,
#         len(complete_distinct_proteins),
#         len(complete_distinct_proteins) - len(uncovered_complete_distinct_proteins),
#         len(uncovered_complete_distinct_proteins),
#         one_fraction_and_data_creation_concat_validated_false_positives_df.shape[0],
#     )

#     fn_df_rows.append(row)

#     # break

# fn_df = pd.DataFrame(fn_df_rows)
# fn_df

# %%
# FNRow = namedtuple(
#     "FNType1Row",
#     [
#         "Chrom",
#         "SwissProtName",
#         "UnknownProbability",
#         "DataCreationRepetition",
#         "Fraction",
#         "NumOfReads",
#         "NumOfDistinctCompleteProteins",
#         "NumOfCoveredDistinctCompleteProteins",
#         "NumOfUncoveredDistinctCompleteProteins",
#         "NumOfTrueCoveringDistinctErroredPartiallyUnknownProteins",
#     ],
# )

# %%
def make_one_fn_row(
    chrom,
    swiss_prot_name,
    unknown_probability,
    repetition,
    fraction,
    max_per_fraction_per_data_creation_merged_distinct_df,
    concat_validated_false_positives_df,
):
    swiss_prot_name = swiss_prot_name.split("_")[0]

    # ic(i, chrom, swiss_prot_name, unknown_probability, repetition, fraction)

    one_fraction_and_data_creation_max_distinct_df = (
        max_per_fraction_per_data_creation_merged_distinct_df.loc[
            (max_per_fraction_per_data_creation_merged_distinct_df["Chrom"] == chrom)
            & (
                max_per_fraction_per_data_creation_merged_distinct_df[
                    "UnknownProbability"
                ]
                == unknown_probability
            )
            & (
                max_per_fraction_per_data_creation_merged_distinct_df[
                    "DataCreationRepetition"
                ]
                == repetition
            )
            & (
                max_per_fraction_per_data_creation_merged_distinct_df["Fraction"]
                == fraction
            ),
            [
                "DataType",
                "MaxNumDistinctProteins",
                "DistinctProteins",
                "Fraction",
                "NumOfReads",
            ],
        ].copy()
    )
    one_fraction_and_data_creation_max_distinct_df[
        "DistinctProteins"
    ] = one_fraction_and_data_creation_max_distinct_df["DistinctProteins"].str.split(
        ","
    )

    num_of_reads = one_fraction_and_data_creation_max_distinct_df["NumOfReads"].iloc[0]

    one_fraction_and_data_creation_concat_validated_false_positives_df = (
        concat_validated_false_positives_df.loc[
            (concat_validated_false_positives_df["Chrom"] == chrom)
            & (
                concat_validated_false_positives_df["UnknownProbability"]
                == unknown_prob
            )
            & (
                concat_validated_false_positives_df["DataCreationRepetition"]
                == repetition
            )
            & (concat_validated_false_positives_df["Fraction"] == fraction)
            & (~concat_validated_false_positives_df["IsFalsePositive"])
        ]
    )

    all_indistinguishable_complete_proteins_of_errored_distinct_proteins = [
        c_protein
        for c_proteins in one_fraction_and_data_creation_concat_validated_false_positives_df[
            "IndistinguishableCompleteProteins"
        ]
        .dropna()
        .tolist()
        for c_protein in c_proteins
    ]
    # ic(len(all_indistinguishable_complete_proteins_of_errored_distinct_proteins));

    complete_distinct_proteins = one_fraction_and_data_creation_max_distinct_df.loc[
        one_fraction_and_data_creation_max_distinct_df["DataType"] == "Complete",
        "DistinctProteins",
    ].values[0]

    uncovered_complete_distinct_proteins = [
        protein
        for protein in complete_distinct_proteins
        if protein
        not in all_indistinguishable_complete_proteins_of_errored_distinct_proteins
    ]

    # row = FNRow(
    #     chrom,
    #     swiss_prot_name,
    #     unknown_probability,
    #     repetition,
    #     fraction,
    #     num_of_reads,
    #     len(complete_distinct_proteins),
    #     len(complete_distinct_proteins) - len(uncovered_complete_distinct_proteins),
    #     len(uncovered_complete_distinct_proteins),
    #     one_fraction_and_data_creation_concat_validated_false_positives_df.shape[0],
    # )
    row = (
        chrom,
        swiss_prot_name,
        unknown_probability,
        repetition,
        fraction,
        num_of_reads,
        len(complete_distinct_proteins),
        len(complete_distinct_proteins) - len(uncovered_complete_distinct_proteins),
        len(uncovered_complete_distinct_proteins),
        one_fraction_and_data_creation_concat_validated_false_positives_df.shape[0],
    )

    return row


# %%
# fn_df_rows = []

# i = 1

# for (chrom, swiss_prot_name), unknown_probability, repetition, fraction in product(
#     swiss_prot_name_by_chrom.items(), unknown_probabilities, repetitions, fractions
# ):
#     row = make_one_fn_row(
#         chrom,
#         swiss_prot_name,
#         unknown_probability,
#         repetition,
#         fraction,
#         max_per_fraction_per_data_creation_merged_distinct_df,
#         concat_validated_false_positives_df,
#     )

#     fn_df_rows.append(row)

#     i += 1

#     if i > 3:
#         break

# fn_df = pd.DataFrame(
#     fn_df_rows,
#     columns=[
#         "Chrom",
#         "SwissProtName",
#         "UnknownProbability",
#         "DataCreationRepetition",
#         "Fraction",
#         "NumOfReads",
#         "NumOfDistinctCompleteProteins",
#         "NumOfCoveredDistinctCompleteProteins",
#         "NumOfUncoveredDistinctCompleteProteins",
#         "NumOfTrueCoveringDistinctErroredPartiallyUnknownProteins",
#     ],
# )
# fn_df

# %%
make_fn_rows_starmap_input = [
    (
        chrom,
        swiss_prot_name,
        unknown_probability,
        repetition,
        fraction,
        max_per_fraction_per_data_creation_merged_distinct_df,
        concat_validated_false_positives_df,
    )
    for (chrom, swiss_prot_name), unknown_probability, repetition, fraction in product(
        swiss_prot_name_by_chrom.items(), unknown_probabilities, repetitions, fractions
    )
]
ic(len(make_fn_rows_starmap_input));

# %%
with Pool(processes=10) as pool:
    fn_df_rows = pool.starmap(func=make_one_fn_row, iterable=make_fn_rows_starmap_input)

fn_df = pd.DataFrame(
    fn_df_rows,
    columns=[
        "Chrom",
        "SwissProtName",
        "UnknownProbability",
        "DataCreationRepetition",
        "Fraction",
        "NumOfReads",
        "NumOfDistinctCompleteProteins",
        "NumOfCoveredDistinctCompleteProteins",
        "NumOfUncoveredDistinctCompleteProteins",
        "NumOfTrueCoveringDistinctErroredPartiallyUnknownProteins",
    ],
)
fn_df

# %%
fn_df["%UncoveredProteins"] = (
    fn_df["NumOfUncoveredDistinctCompleteProteins"]
    .mul(100)
    .div(fn_df["NumOfDistinctCompleteProteins"])
)
# fn_df[
#     "NumOfTrueDistinctErroredPartiallyUnknownProteins/NumOfDistinctCompleteProteins"
# ] = fn_df["NumOfTrueCoveringDistinctErroredPartiallyUnknownProteins"].div(
#     fn_df["NumOfDistinctCompleteProteins"]
# )
fn_df[
    "NumOfTrueDistinctErroredPartiallyUnknownProteins/NumOfCoveredDistinctCompleteProteins"
] = fn_df["NumOfTrueCoveringDistinctErroredPartiallyUnknownProteins"].div(
    fn_df["NumOfCoveredDistinctCompleteProteins"]
)
fn_df

# %% [markdown]
# ## FN type I - fraction of complete uncovered proteins

# %% [markdown]
# ![image.png](attachment:a9a3dbf2-4c41-4bde-96f2-3ac3302bebd2.png)

# %%
dashes = ["solid", "dash", "dot"]
colors = ["red", "blue", "green"]

# unknown_probabilities_dash_dict = {unknown_probability: color for unknown_probability, color in zip(unknown_probabilities, ["red", "blue", "green"])}
# chrom_color_dict = {chrom: dash for chrom, dash in zip(swiss_prot_name_by_chrom, ['dash', 'dot', 'solid'])}

unknown_probabilities_dash_dict = {
    unknown_probability: dash
    for unknown_probability, dash in zip(unknown_probabilities, dashes)
}
chrom_color_dict = {
    chrom: color for chrom, color in zip(swiss_prot_name_by_chrom, colors)
}

fig = make_subplots(
    rows=1,
    cols=1,
    print_grid=False,
    x_title="Simulated reads",
    y_title="Mean % of uncovered distinct proteins",
)

for unknown_probability in unknown_probabilities:
    for chrom, swiss_prot_name in swiss_prot_name_by_chrom.items():
        df = fn_df.loc[
            (fn_df["Chrom"] == chrom)
            & (fn_df["UnknownProbability"] == unknown_probability)
        ]
        df = (
            df.groupby("NumOfReads")
            .agg({"%UncoveredProteins": ["mean", "std"]})
            .reset_index()
        )
        # df

        x = df["NumOfReads"]
        y = df[("%UncoveredProteins", "mean")]
        error_y = df[("%UncoveredProteins", "std")]

        color = chrom_color_dict[chrom]
        dash = unknown_probabilities_dash_dict[unknown_probability]
        swiss_prot_name = swiss_prot_name.split("_")[0]
        # ic(color, dash)
        fig.add_trace(
            go.Scatter(
                x=x,
                y=y,
                error_y=dict(
                    type="data",  # value of error bar given in data coordinates
                    array=error_y,
                    visible=True,
                ),
                # legendgroup=swiss_prot_name,  # this can be any string
                # legendgrouptitle_text=swiss_prot_name,
                # name=unknown_probability,
                legendgroup=unknown_probability,  # this can be any string
                legendgrouptitle_text=unknown_probability,
                name=swiss_prot_name,
                line=dict(color=color, dash=dash),
                opacity=0.7,
            )
        )

width = 700
# height = width * 600 / 1000
height = 550


grouped_fn_df = (
    fn_df.groupby(["Chrom", "UnknownProbability", "Fraction"])
    .agg({"%UncoveredProteins": ["mean", "std"]})
    .reset_index()
)
max_y_axis = (
    grouped_fn_df[("%UncoveredProteins", "mean")]
    .add(grouped_fn_df[("%UncoveredProteins", "std")])
    .max()
    * 1.05
)
fig.update_yaxes(range=[0, max_y_axis])

fig.update_layout(
    width=width,
    height=height,
    legend=dict(
        # title="Swiss prot name, unknown probability",
        # title="Unknown probability, gene",
        # title="Unknown probability, mock gene",
        title="NA probability, mock gene       ",
    ),
)

fig.write_image(
    "Graph assessment - mean % of uncovered distinct proteins.svg",
    width=width,
    height=height,
)

fig.show()

# %% [markdown]
# ## FN type II - ratio of errored/complete proteins

# %% [markdown]
# ![image.png](attachment:c41ed835-e6b5-4a5e-b792-f5c0a8c2ddbc.png)

# %%
dashes = ["solid", "dash", "dot"]
colors = ["red", "blue", "green"]

# unknown_probabilities_dash_dict = {unknown_probability: color for unknown_probability, color in zip(unknown_probabilities, ["red", "blue", "green"])}
# chrom_color_dict = {chrom: dash for chrom, dash in zip(swiss_prot_name_by_chrom, ['dash', 'dot', 'solid'])}

unknown_probabilities_dash_dict = {
    unknown_probability: dash
    for unknown_probability, dash in zip(unknown_probabilities, dashes)
}
chrom_color_dict = {
    chrom: color for chrom, color in zip(swiss_prot_name_by_chrom, colors)
}

fig = make_subplots(
    rows=1,
    cols=1,
    print_grid=False,
    x_title="Simulated reads",
    y_title="Mean ratio of true errored + partially-unknown proteins<br>covering distinct complete proteins",
)

base_y_col = "NumOfTrueDistinctErroredPartiallyUnknownProteins/NumOfCoveredDistinctCompleteProteins"

for unknown_probability in unknown_probabilities:
    for chrom, swiss_prot_name in swiss_prot_name_by_chrom.items():
        df = fn_df.loc[
            (fn_df["Chrom"] == chrom)
            & (fn_df["UnknownProbability"] == unknown_probability)
        ]
        df = (
            df.groupby("NumOfReads")
            .agg(
                {
                    base_y_col: [
                        "mean",
                        "std",
                    ]
                }
            )
            .reset_index()
        )
        # df

        x = df["NumOfReads"]
        y = df[
            (
                base_y_col,
                "mean",
            )
        ]
        error_y = df[
            (
                base_y_col,
                "std",
            )
        ]

        color = chrom_color_dict[chrom]
        dash = unknown_probabilities_dash_dict[unknown_probability]
        swiss_prot_name = swiss_prot_name.split("_")[0]
        # ic(color, dash)
        fig.add_trace(
            go.Scatter(
                x=x,
                y=y,
                error_y=dict(
                    type="data",  # value of error bar given in data coordinates
                    array=error_y,
                    visible=True,
                ),
                # legendgroup=swiss_prot_name,  # this can be any string
                # legendgrouptitle_text=swiss_prot_name,
                # name=unknown_probability,
                legendgroup=unknown_probability,  # this can be any string
                legendgrouptitle_text=unknown_probability,
                name=swiss_prot_name,
                line=dict(color=color, dash=dash),
                opacity=0.7,
            )
        )

width = 700
# height = width * 600 / 1000
height = 550


grouped_fn_df = (
    fn_df.groupby(["Chrom", "UnknownProbability", "Fraction"])
    .agg({base_y_col: ["mean", "std"]})
    .reset_index()
)
max_y_axis = (
    grouped_fn_df[(base_y_col, "mean")].add(grouped_fn_df[(base_y_col, "std")]).max()
    * 1.05
)
fig.update_yaxes(
    range=[0, max_y_axis],
    # type="log",
    # range=[-1, np.log(max_y_axis)/np.log(10)]
)

fig.update_layout(
    width=width,
    height=height,
    legend=dict(
        # title="Swiss prot name, unknown probability",
        # title="Unknown probability, gene",
        # title="Unknown probability, mock gene",
        title="NA probability, mock gene       ",
    ),
)

fig.write_image(
    "Graph assessment - mean ratio of true errored proteins covering distinct complete proteins.svg",
    width=width,
    height=height,
)

fig.show()
