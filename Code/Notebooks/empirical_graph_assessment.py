# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.7
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
from itertools import chain, product
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
unique_proteins_files_dir = out_dir
false_positives_files_fofn = Path(out_dir, "FalsePositivesOutFiles.fofn")

known_sites_csv_file = (
    "/private7/projects/Combinatorics/D.pealeii/Annotations/D.pea.EditingSites.csv"
)

n_reads = 100_000

unique_proteins_first_pos_col = 14

seed = 1892

fractions = [0.2, 0.4, 0.6, 0.8, 1.0]

# %%
import matplotlib.colors as mcolors


def css_color_with_opacity(color_name: str, x: float) -> str:
    """
    Convert a named CSS color to an RGBA color with the given opacity.

    Parameters:
        color_name (str): The CSS color name (e.g., "red").
        x (float): The opacity percentage (0-100).

    Returns:
        str: The RGBA representation of the color with opacity.
    """
    if color_name.lower() not in mcolors.CSS4_COLORS:
        raise ValueError("Invalid CSS color name")

    rgb = mcolors.to_rgb(mcolors.CSS4_COLORS[color_name.lower()])
    alpha = max(0, min(1, x / 100))  # Ensure alpha is in range [0,1]

    return f"rgba({int(rgb[0] * 255)}, {int(rgb[1] * 255)}, {int(rgb[2] * 255)}, {alpha:.2f})"


# Example usage
print(css_color_with_opacity("red", 50))  # Output: 'rgba(255, 0, 0, 0.50)'

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

# %% jupyter={"source_hidden": true}
# max_per_fraction_merged_distinct_df = (
#     merged_distinct_df.loc[merged_distinct_df["IsMaxDistinctProteinsPerGroup"],]
#     .groupby(["Chrom", "SwissProtName", "UnknownProbability", "DataType", "Fraction"])
#     .sample(n=1, random_state=seed)
#     .reset_index(drop=True)
#     .rename(columns={"NumDistinctProteins": "MaxNumDistinctProteins"})
# )
# max_per_fraction_merged_distinct_df

# %% jupyter={"source_hidden": true}
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

# %% jupyter={"source_hidden": true}
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

# %% jupyter={"source_hidden": true}
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

# %% jupyter={"source_hidden": true}
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
wide_max_per_fraction_per_data_creation_merged_distinct_df

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

# %% jupyter={"source_hidden": true}
# unknown_probabilities

# %% jupyter={"source_hidden": true}
# swiss_prot_name_by_chrom

# %% jupyter={"source_hidden": true}
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

# %% [markdown]
# # Unique proteins

# %%
# base_unique_proteins_first_col_pos = 14

# %%
# # !zcat /private7/projects/Combinatorics/Simulations/GraphAssessment/comp141434_c0_seq1.CBPC1_HUMAN.UP0_18.Rep8.Complete.UniqueProteins.tsv.gz | head -n 2

# %%
unique_proteins_first_col_pos = base_unique_proteins_first_col_pos

complete_unique_proteins_dfs = []

for (chrom, swiss_prot_name), unknown_probability, repetition in product(
    swiss_prot_name_by_chrom.items(), unknown_probabilities, repetitions
):
    complete_unique_proteins_file = Path(
        unique_proteins_files_dir,
        f"{chrom}.{swiss_prot_name}.UP{str(unknown_probability).replace(".", "_")}.{repetition}.Complete.UniqueProteins.tsv.gz",
    )
    complete_unique_proteins_df = pd.read_csv(
        complete_unique_proteins_file,
        sep="\t",
        dtype={"Protein": str, "Reads": str},
        usecols=[
            "Gene",
            "Protein",
            "MinNonSyns",
            "MaxNonSyns",
            "MinNonSynsFrequency",
            "MaxNonSynsFrequency",
            "NumOfReads",
        ],
    )
    swiss_prot_name = swiss_prot_name.split("_")[0]
    complete_unique_proteins_df.insert(0, "Chrom", chrom)
    complete_unique_proteins_df.insert(1, "SwissProtName", swiss_prot_name)
    complete_unique_proteins_df.insert(2, "UnknownProbability", unknown_probability)
    complete_unique_proteins_df.insert(3, "DataCreationRepetition", repetition)
    complete_unique_proteins_dfs.append(complete_unique_proteins_df)

ic(len(complete_unique_proteins_dfs))

# unique_proteins_first_col_pos = base_unique_proteins_first_col_pos + 4

concat_complete_unique_proteins_df = pd.concat(
    complete_unique_proteins_dfs, ignore_index=True
)
concat_complete_unique_proteins_df

# %%
complete_unique_proteins_dfs = [
    pd.read_csv(unique_proteins_file, sep=sep, dtype={"Protein": str, "Reads": str})
    for unique_proteins_file in unique_proteins_files
]
for chrom, unique_proteins_df in zip(chroms, unique_proteins_dfs):
    unique_proteins_df.insert(0, "Chrom", chrom)
# the position of the first column needs to be updated due to the Chrom col insertion
unique_proteins_first_col_pos += 1
# for unique_proteins_df in unique_proteins_dfs:
#     unique_proteins_df.rename(
#         columns={
#             col: col.replace("Transcripts", "UniqueReads")
#             for col in unique_proteins_df.columns[:unique_proteins_first_col_pos]
#             if "Transcripts" in col
#         },
#         inplace=True,
#     )
unique_proteins_dfs[0]

# %% [markdown]
# # Editing sites

# %%
chroms = list(swiss_prot_name_by_chrom.keys())
chroms

# %%
cols = ["Trinity id", "Editing location (base1)", "SwissProt name", "Editing Level"]
editing_sites_df = (
    pd.read_csv(known_sites_csv_file)
    .filter(cols)
    .rename(
        columns={
            "Trinity id": "Chrom",
            "Editing location (base1)": "Position",
            "SwissProt name": "SwissProtName",
            "Editing Level": "EditingLevel",
        }
    )
    .sort_values(["Chrom", "Position"])
    .reset_index(drop=True)
)
editing_sites_df = editing_sites_df.loc[
    editing_sites_df["Chrom"].isin(chroms)
].reset_index(drop=True)
editing_sites_df["Position"] = editing_sites_df["Position"] - 1
editing_sites_df["%EditingLevel"] = editing_sites_df["EditingLevel"].mul(100)
editing_sites_df["SwissProtName"] = (
    editing_sites_df["SwissProtName"].str.split("_").str[0]
)
editing_sites_df

# %%
editing_sites_df.groupby(["Chrom", "SwissProtName"])["%EditingLevel"].describe().round(
    2
)

# %%
colors = ["red", "blue", "green"]

swiss_prot_color_dict = {
    swiss_prot_name.split("_")[0]: color
    for swiss_prot_name, color in zip(swiss_prot_name_by_chrom.values(), colors)
}


fig = make_subplots(
    rows=1,
    cols=1,
    print_grid=False,
    x_title="Gene",
    y_title="Known sites' editing levels [%]",
)

for chrom, swiss_prot_name in swiss_prot_name_by_chrom.items():
    df = editing_sites_df.loc[(editing_sites_df["Chrom"] == chrom)]
    x = df["SwissProtName"]
    y = df["%EditingLevel"]

    color = chrom_color_dict[chrom]
    swiss_prot_name = swiss_prot_name.split("_")[0]
    # ic(color, dash)
    fig.add_trace(
        # go.Box(
        go.Violin(
            x=x,
            y=y,
            name=swiss_prot_name,
            # line=dict(color=color,),
            line=dict(
                color="black",
            ),
            # line=dict(color="grey",),
            # line=dict(color=css_color_with_opacity("black", 70)),
            marker=dict(
                color=color,
            ),
            fillcolor=css_color_with_opacity(color, 30),
            # opacity=0.7,
            points="all",
            box_visible=True,
            meanline_visible=True,
        )
    )

width = 700
# height = width * 600 / 1000
height = 400

# fig.update_yaxes(tick0=0, dtick=10)

fig.update_layout(
    width=width,
    height=height,
    showlegend=False,
)

# fig.write_image(
#     "Graph assessment - mean % of distinct proteins recovered.svg",
#     width=width,
#     height=height,
# )

fig.show()

# %%
colors = ["red", "blue", "green"]

# chrom_color_dict = {
#     chrom: color for chrom, color in zip(swiss_prot_name_by_chrom, colors)
# }
swiss_prot_color_dict = {
    swiss_prot_name.split("_")[0]: color
    for swiss_prot_name, color in zip(swiss_prot_name_by_chrom.values(), colors)
}

fig = px.histogram(
    editing_sites_df,
    x="%EditingLevel",
    color="SwissProtName",
    color_discrete_map=swiss_prot_color_dict,
    # marginal="box"
    # opacity=0.5,
    # facet_col="SwissProtName"
)

width = 500
# height = width * 600 / 1000
height = 350
fig.update_xaxes(title="Known sites")
fig.update_yaxes(
    title="Editing levels [%]",
    dtick=10,
    # type="log",
)
fig.update_layout(
    width=width,
    height=height,
    # barmode="overlay",
    barmode="group",
    legend=dict(
        title="Gene       ",
    ),
    margin_r=140,
)

# fig.write_image(
#     "Graph assessment - mean % of distinct proteins recovered.svg",
#     width=width,
#     height=height,
# )

fig.show()

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

# %% jupyter={"source_hidden": true}
# Out[132]["NumOfIndistinguishableCompleteProteins"].value_counts(dropna=False)

# %% jupyter={"source_hidden": true}
# concat_false_positives_df.loc[
#     concat_false_positives_df["NumOfIndistinguishableCompleteProteins"] == 0,
# ]

# %% jupyter={"source_hidden": true}
# Out[135]["IndistinguishableCompleteProteins"].value_counts(dropna=False)

# %% jupyter={"source_hidden": true}
# false_positives_df["IndistinguishableCompleteProteins"].replace(np.nan, [])

# %%
# next I should use the false-positives and this df to test for false-positives only inside the largest MIS of each dataset

max_per_fraction_per_data_creation_merged_distinct_df

# %%
chrom = "comp141434_c0_seq1"
unknown_probability = 0.13
repetition = "Rep4"
fraction = 0.2

# %% jupyter={"source_hidden": true}
# chrom = "comp140826_c0_seq1"
# swiss_prot_name = "VPS8"
# unknown_prob = 0.09
# repetition = "Rep1"
# # fraction = 0.2
# fraction = 0.6
# # fraction = 1

# %%
one_fraction_and_data_creation_max_distinct_df = (
    max_per_fraction_per_data_creation_merged_distinct_df.loc[
        (max_per_fraction_per_data_creation_merged_distinct_df["Chrom"] == chrom)
        & (
            max_per_fraction_per_data_creation_merged_distinct_df["UnknownProbability"]
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
] = one_fraction_and_data_creation_max_distinct_df["DistinctProteins"].str.split(",")
one_fraction_and_data_creation_max_distinct_df

# %%
num_of_reads = one_fraction_and_data_creation_max_distinct_df["NumOfReads"].iloc[0]
num_of_reads

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
def subset_indistinguishable_complete_proteins_in_fraction_and_data_creation(complete_distinct_proteins_in_fraction_and_data_creation, all_indistinguishable_complete_proteins):
    """
    An errored protein `x` (in a certain data creation) is indistinguishable from 3 complete proteins `y1`, `y2`, `y3`. 
    However, only `y1`, `y2` are considered distinct in a MIS of some fraction, so return these.
    Return `np.nan` if none exist.
    """
    try:
        indistinguishable_complete_proteins_in_fraction_and_data_creation = [p for p in all_indistinguishable_complete_proteins if p in complete_distinct_proteins_in_fraction_and_data_creation]
        if len(indistinguishable_complete_proteins_in_fraction_and_data_creation) == 0:
            indistinguishable_complete_proteins_in_fraction_and_data_creation = np.nan
    except:
        indistinguishable_complete_proteins_in_fraction_and_data_creation = np.nan
    return indistinguishable_complete_proteins_in_fraction_and_data_creation


# %%
one_fraction_and_data_creation_false_positives_df = concat_false_positives_df.loc[
    (concat_false_positives_df["Chrom"] == chrom)
    # & (concat_false_positives_df["SwissProtName"] == swiss_prot_name)
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

one_fraction_and_data_creation_false_positives_df = one_fraction_and_data_creation_false_positives_df.sort_values("NumOfIndistinguishableCompleteProteins")

one_fraction_and_data_creation_false_positives_df.insert(
    one_fraction_and_data_creation_false_positives_df.columns.get_loc(
        "IndistinguishableCompleteProteins"
    )
    + 1,
    "IndistinguishableCompleteProteinsInFractionAndDataCreation",
    one_fraction_and_data_creation_false_positives_df.apply(
        lambda x: 
        subset_indistinguishable_complete_proteins_in_fraction_and_data_creation(complete_distinct_proteins, x["IndistinguishableCompleteProteins"]), 
        axis=1
    )
)
one_fraction_and_data_creation_false_positives_df.insert(
    one_fraction_and_data_creation_false_positives_df.columns.get_loc(
        "IndistinguishableCompleteProteinsInFractionAndDataCreation"
    )
    + 1,
    "NumIndistinguishableCompleteProteinsInFractionAndDataCreation",
    one_fraction_and_data_creation_false_positives_df["IndistinguishableCompleteProteinsInFractionAndDataCreation"].apply(lambda x: len(x) if type(x) == list else 0)
)
assert one_fraction_and_data_creation_false_positives_df.loc[
    (one_fraction_and_data_creation_false_positives_df["IndistinguishableCompleteProteins"].isna())
& (~one_fraction_and_data_creation_false_positives_df["IndistinguishableCompleteProteinsInFractionAndDataCreation"].isna())
].empty

one_fraction_and_data_creation_false_positives_df = one_fraction_and_data_creation_false_positives_df.sort_values("NumIndistinguishableCompleteProteinsInFractionAndDataCreation")

one_fraction_and_data_creation_false_positives_df

# %%
# is_false_positive = one_fraction_and_data_creation_false_positives_df.apply(
#     lambda x: (x["NumOfIndistinguishableCompleteProteins"] == 0)
#     or (
#         x["NumOfIndistinguishableCompleteProteins"] > 0
#         and set(x["IndistinguishableCompleteProteins"]) & complete_distinct_proteins
#         == set()
#     ),
#     axis=1,
# )

is_false_positive = one_fraction_and_data_creation_false_positives_df["NumIndistinguishableCompleteProteinsInFractionAndDataCreation"].eq(0)

# is_false_positive

ic(100 * is_false_positive.sum() / len(is_false_positive))

is_false_positive

# %% jupyter={"source_hidden": true}
# one_fraction_and_data_creation_false_positives_df[
#     "IsFalsePositive"
# ] = is_false_positive.tolist()
# one_fraction_and_data_creation_false_positives_df

# %% jupyter={"source_hidden": true}
# one_fraction_and_data_creation_false_positives_df.loc[
#     (one_fraction_and_data_creation_false_positives_df["IsFalsePositive"])
#     & (
#         one_fraction_and_data_creation_false_positives_df[
#             "NumOfIndistinguishableCompleteProteins"
#         ]
#         > 1
#     )
# ].head(100)

# %%
# "5cX" in complete_distinct_proteins

# %% jupyter={"source_hidden": true}
# "axF" in complete_distinct_proteins

# %% jupyter={"source_hidden": true}
# "6Rp" in complete_distinct_proteins

# %% jupyter={"source_hidden": true}
# "brB" in complete_distinct_proteins

# %% jupyter={"source_hidden": true}
# "ahs" in complete_distinct_proteins

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
    one_fraction_and_data_creation_max_distinct_df["DistinctProteins"] = (
        one_fraction_and_data_creation_max_distinct_df["DistinctProteins"].str.split(
            ","
        )
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

    one_fraction_and_data_creation_false_positives_df.insert(
        one_fraction_and_data_creation_false_positives_df.columns.get_loc(
            "IndistinguishableCompleteProteins"
        )
        + 1,
        "IndistinguishableCompleteProteinsInFractionAndDataCreation",
        one_fraction_and_data_creation_false_positives_df.apply(
            lambda x: 
            subset_indistinguishable_complete_proteins_in_fraction_and_data_creation(complete_distinct_proteins, x["IndistinguishableCompleteProteins"]), 
            axis=1
        )
    )
    one_fraction_and_data_creation_false_positives_df.insert(
        one_fraction_and_data_creation_false_positives_df.columns.get_loc(
            "IndistinguishableCompleteProteinsInFractionAndDataCreation"
        )
        + 1,
        "NumIndistinguishableCompleteProteinsInFractionAndDataCreation",
        one_fraction_and_data_creation_false_positives_df["IndistinguishableCompleteProteinsInFractionAndDataCreation"].apply(lambda x: len(x) if type(x) == list else 0)
    )
    assert one_fraction_and_data_creation_false_positives_df.loc[
        (one_fraction_and_data_creation_false_positives_df["IndistinguishableCompleteProteins"].isna())
    & (~one_fraction_and_data_creation_false_positives_df["IndistinguishableCompleteProteinsInFractionAndDataCreation"].isna())
    ].empty

    # is_false_positive = one_fraction_and_data_creation_false_positives_df.apply(
    #     lambda x: (x["NumOfIndistinguishableCompleteProteins"] == 0)
    #     or (
    #         x["NumOfIndistinguishableCompleteProteins"] > 0
    #         and set(x["IndistinguishableCompleteProteins"]) & complete_distinct_proteins
    #         == set()
    #     ),
    #     axis=1,
    # )
    is_false_positive = one_fraction_and_data_creation_false_positives_df["NumIndistinguishableCompleteProteinsInFractionAndDataCreation"].eq(0)

    # validated_false_positives_df = one_fraction_and_data_creation_false_positives_df.loc[is_false_positive]
    # validated_false_positives_dfs.append(validated_false_positives_df)

    one_fraction_and_data_creation_false_positives_df.loc[:, "IsFalsePositive"] = (
        is_false_positive
    )
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
        one_fraction_and_data_creation_max_distinct_df["DistinctProteins"] = (
            one_fraction_and_data_creation_max_distinct_df[
                "DistinctProteins"
            ].str.split(",")
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

# %%

# %%

# %%

# %%

# %% [markdown]
# # False negatives

# %%
max_per_fraction_per_data_creation_merged_distinct_df

# %%
concat_validated_false_positives_df

# %% jupyter={"source_hidden": true}
# concat_validated_false_positives_df.drop_duplicates(
#     ["Chrom", "SwissProtName", "UnknownProbability", "Errored+PartiallyUnknownProtein"]
# )["NumOfIndistinguishableCompleteProteins"].describe()

# %% jupyter={"source_hidden": true}
# concat_validated_false_positives_df.loc[
#     ~concat_validated_false_positives_df["IsFalsePositive"]
# ].drop_duplicates(
#     ["Chrom", "SwissProtName", "UnknownProbability", "Errored+PartiallyUnknownProtein"]
# )[
#     "NumOfIndistinguishableCompleteProteins"
# ].describe()

# %% jupyter={"source_hidden": true}
# concat_validated_false_positives_df.loc[
#     ~concat_validated_false_positives_df["IsFalsePositive"]
# ].drop_duplicates(
#     ["Chrom", "SwissProtName", "UnknownProbability", "Errored+PartiallyUnknownProtein"]
# )[
#     "NumIndistinguishableCompleteProteinsInFractionAndDataCreation"
# ].describe()

# %% jupyter={"source_hidden": true}
# concat_validated_false_positives_df.loc[
#     ~concat_validated_false_positives_df["IsFalsePositive"]
# ].drop_duplicates(
#     [
#         "Chrom",
#         "SwissProtName",
#         "UnknownProbability",
#         "Fraction",
#         "Errored+PartiallyUnknownProtein",
#     ]
# ).groupby(
#     ["Fraction"]
# )[
#     "NumOfIndistinguishableCompleteProteins"
# ].describe()

# %%
# concat_validated_false_positives_df.loc[
#     ~concat_validated_false_positives_df["IsFalsePositive"]
# ].drop_duplicates(
#     [
#         "Chrom",
#         "SwissProtName",
#         "UnknownProbability",
#         "Fraction",
#         "Errored+PartiallyUnknownProtein",
#     ]
# ).groupby(
#     ["Fraction"]
# )[
#     "NumIndistinguishableCompleteProteinsInFractionAndDataCreation"
# ].describe()

# %%
# concat_validated_false_positives_df.loc[
#     ~concat_validated_false_positives_df["IsFalsePositive"]
# ].groupby(
#     ["Chrom",
#         "SwissProtName",
#         "UnknownProbability","Fraction"]
# )[
#     "NumIndistinguishableCompleteProteinsInFractionAndDataCreation"
# ].describe()

# %% [markdown]
# ## Test

# %%
# chrom = "comp141434_c0_seq1"
# unknown_probability = 0.09
# repetition = "Rep1"
# fraction = 0.2

# %%
chrom = "comp141434_c0_seq1"
unknown_probability = 0.13
repetition = "Rep4"
fraction = 0.2

# %%
one_fraction_and_data_creation_max_distinct_df = (
    max_per_fraction_per_data_creation_merged_distinct_df.loc[
        (max_per_fraction_per_data_creation_merged_distinct_df["Chrom"] == chrom)
        & (
            max_per_fraction_per_data_creation_merged_distinct_df["UnknownProbability"]
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
            "UnknownProbability",
            "DataCreationRepetition",
            "NumOfReads",
        ],
    ].copy()
)
one_fraction_and_data_creation_max_distinct_df["DistinctProteins"] = (
    one_fraction_and_data_creation_max_distinct_df["DistinctProteins"].str.split(",")
)
one_fraction_and_data_creation_max_distinct_df

# %%
num_of_reads = one_fraction_and_data_creation_max_distinct_df["NumOfReads"].iloc[0]
ic(num_of_reads)

complete_distinct_proteins = one_fraction_and_data_creation_max_distinct_df.loc[
    one_fraction_and_data_creation_max_distinct_df["DataType"] == "Complete",
    "DistinctProteins",
].values[0]
ic(len(complete_distinct_proteins));

# %%
one_fraction_and_data_creation_concat_validated_non_false_positives_df = (
    concat_validated_false_positives_df.loc[
        (concat_validated_false_positives_df["Chrom"] == chrom)
        & (
            concat_validated_false_positives_df["UnknownProbability"]
            == unknown_probability
        )
        & (concat_validated_false_positives_df["DataCreationRepetition"] == repetition)
        & (concat_validated_false_positives_df["Fraction"] == fraction)
        & (~concat_validated_false_positives_df["IsFalsePositive"])
    ]
).copy()

# one_fraction_and_data_creation_concat_validated_non_false_positives_df[
#     "CoveredIndistinguishableCompleteProteins"
# ] = one_fraction_and_data_creation_concat_validated_non_false_positives_df[
#     "IndistinguishableCompleteProteins"
# ].apply(
#     lambda x: [y for y in x if y in complete_distinct_proteins]
# )
# one_fraction_and_data_creation_concat_validated_non_false_positives_df[
#     "NumOfCoveredIndistinguishableCompleteProteins"
# ] = one_fraction_and_data_creation_concat_validated_non_false_positives_df[
#     "CoveredIndistinguishableCompleteProteins"
# ].apply(
#     len
# )
one_fraction_and_data_creation_concat_validated_non_false_positives_df[
    "CoveredIndistinguishableCompleteProteins"
] = one_fraction_and_data_creation_concat_validated_non_false_positives_df[
    "IndistinguishableCompleteProteinsInFractionAndDataCreation"
]
one_fraction_and_data_creation_concat_validated_non_false_positives_df[
    "NumOfCoveredIndistinguishableCompleteProteins"
] = one_fraction_and_data_creation_concat_validated_non_false_positives_df[
    "NumIndistinguishableCompleteProteinsInFractionAndDataCreation"
]
assert (
    one_fraction_and_data_creation_concat_validated_non_false_positives_df[
        "NumOfCoveredIndistinguishableCompleteProteins"
    ]
    .ge(1)
    .all()
)

one_fraction_and_data_creation_concat_validated_non_false_positives_df

# %%
# one_fraction_and_data_creation_concat_validated_non_false_positives_df["NumOfCoveredIndistinguishableCompleteProteins"].describe()

# %%
# len(list(chain.from_iterable(one_fraction_and_data_creation_concat_validated_non_false_positives_df["CoveredIndistinguishableCompleteProteins"].values)))

# %%
covered_complete_distinct_proteins = set(
    chain.from_iterable(
        one_fraction_and_data_creation_concat_validated_non_false_positives_df[
            "CoveredIndistinguishableCompleteProteins"
        ]
    )
)
# ic(len(covered_complete_distinct_proteins));

uncovered_complete_distinct_proteins = (
    set(complete_distinct_proteins) - covered_complete_distinct_proteins
)

num_of_distinct_complete_proteins = len(complete_distinct_proteins)
num_of_covered_distinct_complete_proteins = len(covered_complete_distinct_proteins)
num_of_uncovered_distinct_complete_proteins = len(uncovered_complete_distinct_proteins)
num_of_true_covering_distinct_errored_partially_unknown_proteins = (
    one_fraction_and_data_creation_concat_validated_non_false_positives_df.shape[0]
)

ic(
    num_of_distinct_complete_proteins,
    num_of_covered_distinct_complete_proteins,
    num_of_uncovered_distinct_complete_proteins,
    num_of_true_covering_distinct_errored_partially_unknown_proteins,
);

# %%
num_of_covered_distinct_complete_proteins + num_of_uncovered_distinct_complete_proteins

# %%
make_one_fn_row(
    chrom,
    swiss_prot_name_by_chrom[chrom],
    unknown_probability,
    repetition,
    fraction,
    max_per_fraction_per_data_creation_merged_distinct_df,
    concat_validated_false_positives_df,
)[:9]


# %% [markdown]
# ## Run

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
    one_fraction_and_data_creation_max_distinct_df["DistinctProteins"] = (
        one_fraction_and_data_creation_max_distinct_df["DistinctProteins"].str.split(
            ","
        )
    )

    num_of_reads = one_fraction_and_data_creation_max_distinct_df["NumOfReads"].iloc[0]

    complete_distinct_proteins = one_fraction_and_data_creation_max_distinct_df.loc[
        one_fraction_and_data_creation_max_distinct_df["DataType"] == "Complete",
        "DistinctProteins",
    ].values[0]

    one_fraction_and_data_creation_concat_validated_non_false_positives_df = (
        concat_validated_false_positives_df.loc[
            (concat_validated_false_positives_df["Chrom"] == chrom)
            & (
                concat_validated_false_positives_df["UnknownProbability"]
                == unknown_probability
            )
            & (
                concat_validated_false_positives_df["DataCreationRepetition"]
                == repetition
            )
            & (concat_validated_false_positives_df["Fraction"] == fraction)
            & (~concat_validated_false_positives_df["IsFalsePositive"])
        ].copy()
    )

    # one_fraction_and_data_creation_concat_validated_non_false_positives_df[
    #     "CoveredIndistinguishableCompleteProteins"
    # ] = one_fraction_and_data_creation_concat_validated_non_false_positives_df[
    #     "IndistinguishableCompleteProteins"
    # ].apply(
    #     lambda x: [y for y in x if y in complete_distinct_proteins]
    # )
    # one_fraction_and_data_creation_concat_validated_non_false_positives_df[
    #     "NumOfCoveredIndistinguishableCompleteProteins"
    # ] = one_fraction_and_data_creation_concat_validated_non_false_positives_df[
    #     "CoveredIndistinguishableCompleteProteins"
    # ].apply(
    #     len
    # )
    one_fraction_and_data_creation_concat_validated_non_false_positives_df[
        "CoveredIndistinguishableCompleteProteins"
    ] = one_fraction_and_data_creation_concat_validated_non_false_positives_df[
        "IndistinguishableCompleteProteinsInFractionAndDataCreation"
    ]
    one_fraction_and_data_creation_concat_validated_non_false_positives_df[
        "NumOfCoveredIndistinguishableCompleteProteins"
    ] = one_fraction_and_data_creation_concat_validated_non_false_positives_df[
        "NumIndistinguishableCompleteProteinsInFractionAndDataCreation"
    ]
    assert (
        one_fraction_and_data_creation_concat_validated_non_false_positives_df[
            "NumOfCoveredIndistinguishableCompleteProteins"
        ]
        .ge(1)
        .all()
    )

    covered_complete_distinct_proteins = set(
        chain.from_iterable(
            one_fraction_and_data_creation_concat_validated_non_false_positives_df[
                "CoveredIndistinguishableCompleteProteins"
            ]
        )
    )
    uncovered_complete_distinct_proteins = (
        set(complete_distinct_proteins) - covered_complete_distinct_proteins
    )

    num_of_distinct_complete_proteins = len(complete_distinct_proteins)
    num_of_covered_distinct_complete_proteins = len(covered_complete_distinct_proteins)
    num_of_uncovered_distinct_complete_proteins = len(
        uncovered_complete_distinct_proteins
    )
    num_of_true_covering_distinct_errored_partially_unknown_proteins = (
        one_fraction_and_data_creation_concat_validated_non_false_positives_df.shape[0]
    )

    complete_distinct_proteins = list(covered_complete_distinct_proteins) + list(
        uncovered_complete_distinct_proteins
    )
    coverage_status = [True] * len(covered_complete_distinct_proteins) + [False] * len(
        uncovered_complete_distinct_proteins
    )

    row = (
        chrom,
        swiss_prot_name,
        unknown_probability,
        repetition,
        fraction,
        num_of_reads,
        num_of_distinct_complete_proteins,
        num_of_covered_distinct_complete_proteins,
        num_of_uncovered_distinct_complete_proteins,
        num_of_true_covering_distinct_errored_partially_unknown_proteins,
        complete_distinct_proteins,
        coverage_status,
    )

    return row


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
ic(len(make_fn_rows_starmap_input))

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
        "DistinctCompleteProteins",
        "AreDistinctCompleteProteinsCovered",
    ],
)

fn_df["%UncoveredProteins"] = (
    fn_df["NumOfUncoveredDistinctCompleteProteins"]
    .mul(100)
    .div(fn_df["NumOfDistinctCompleteProteins"])
)
fn_df["%CoveredProteins"] = (
    fn_df["NumOfCoveredDistinctCompleteProteins"]
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

# each errored protein covers at least one complete protein
assert fn_df["NumOfTrueDistinctErroredPartiallyUnknownProteins/NumOfCoveredDistinctCompleteProteins"].le(1).all()

fn_df

# %%

# %% [markdown]
# ## FN type I - fraction of complete uncovered proteins

# %%
chrom = "comp141434_c0_seq1"
unknown_probability = 0.13
repetition = "Rep4"
fraction = 0.2

# %%
true_positives_wide_max_per_fraction_per_data_creation_merged_distinct_df.loc[
    (true_positives_wide_max_per_fraction_per_data_creation_merged_distinct_df["Chrom"] == chrom)
    & (
        true_positives_wide_max_per_fraction_per_data_creation_merged_distinct_df["UnknownProbability"]
        == unknown_probability
    )
    & (
        true_positives_wide_max_per_fraction_per_data_creation_merged_distinct_df[
            "DataCreationRepetition"
        ]
        == repetition
    )
    & (
        true_positives_wide_max_per_fraction_per_data_creation_merged_distinct_df["Fraction"]
        == fraction
    )
]

# %% jupyter={"source_hidden": true}
fn_df.loc[
    (fn_df["Chrom"] == chrom)
    & (
        fn_df["UnknownProbability"]
        == unknown_probability
    )
    & (
        fn_df[
            "DataCreationRepetition"
        ]
        == repetition
    )
    & (
        fn_df["Fraction"]
        == fraction
    )
]

# %%

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
# ### What rules which complete proteins don't get covered?

# %%
# complete_proteins_inclusion_df = (
#     fn_df.drop(
#         columns=[
#             "NumOfReads",
#             "NumOfDistinctCompleteProteins",
#             "NumOfCoveredDistinctCompleteProteins",
#             "NumOfUncoveredDistinctCompleteProteins",
#             "NumOfTrueCoveringDistinctErroredPartiallyUnknownProteins",
#             "%UncoveredProteins",
#             "NumOfTrueDistinctErroredPartiallyUnknownProteins/NumOfCoveredDistinctCompleteProteins",
#         ]
#     )
#     .explode(["DistinctCompleteProteins", "AreDistinctCompleteProteinsCovered"])
#     .rename(
#         columns={
#             "DistinctCompleteProteins": "DistinctCompleteProtein",
#             "AreDistinctCompleteProteinsCovered": "Covered",
#         }
#     )
# )
# complete_proteins_inclusion_df

# %%
# complete_proteins_inclusion_df.loc[complete_proteins_inclusion_df["Covered"]]

# %%
# (
#     complete_proteins_inclusion_df
#     .merge(
#         concat_complete_unique_proteins_df.drop(columns=["Gene"],
#         how="left",
#         left_on=["Chrom", "SwissProtName", "UnknownProbability", "DataCreationRepetition", "DistinctCompleteProtein"],
#         right_on=["Chrom", "SwissProtName", "UnknownProbability", "DataCreationRepetition", "Protein"]
#     )
# )

# %%
complete_proteins_inclusion_df = (
    fn_df.drop(
        columns=[
            "NumOfReads",
            "NumOfDistinctCompleteProteins",
            "NumOfCoveredDistinctCompleteProteins",
            "NumOfUncoveredDistinctCompleteProteins",
            "NumOfTrueCoveringDistinctErroredPartiallyUnknownProteins",
            "%UncoveredProteins",
            "NumOfTrueDistinctErroredPartiallyUnknownProteins/NumOfCoveredDistinctCompleteProteins",
        ]
    )
    .explode(["DistinctCompleteProteins", "AreDistinctCompleteProteinsCovered"])
    .rename(
        columns={
            "DistinctCompleteProteins": "DistinctCompleteProtein",
            "AreDistinctCompleteProteinsCovered": "Covered",
        }
    )
    .merge(
        concat_complete_unique_proteins_df.drop(columns=["Gene"]),
        how="left",
        left_on=["Chrom", "SwissProtName", "UnknownProbability", "DataCreationRepetition", "DistinctCompleteProtein"],
        right_on=["Chrom", "SwissProtName", "UnknownProbability", "DataCreationRepetition", "Protein"]
    )
)

complete_proteins_inclusion_df

# %%
agg_complete_proteins_inclusion_df = complete_proteins_inclusion_df.groupby(["Chrom", "SwissProtName", "UnknownProbability", "DataCreationRepetition", "DistinctCompleteProtein", "NumOfReads"]).agg(
    CoveredInFractions=("Fraction", list)
).reset_index()
agg_complete_proteins_inclusion_df["CoveredInXFractions"] = agg_complete_proteins_inclusion_df["CoveredInFractions"].apply(len)

# so far, the agg_complete_proteins_inclusion_df contains info of complete proteins which are covered at least once - 
# but we also need to account for those never included
agg_complete_proteins_inclusion_df = agg_complete_proteins_inclusion_df.merge(
    concat_complete_unique_proteins_df.drop(columns=["Gene", "MinNonSyns", "MaxNonSyns", "MinNonSynsFrequency", "MaxNonSynsFrequency"]),
    how="outer",
    left_on=["Chrom", "SwissProtName", "UnknownProbability", "DataCreationRepetition", "DistinctCompleteProtein", "NumOfReads"],
    right_on=["Chrom", "SwissProtName", "UnknownProbability", "DataCreationRepetition", "Protein", "NumOfReads"]
)
# it turns out that each unique protein is covered at least once - 
# we use the following assertion, which will save us time in reformatting
# the df to include 0-times covered unique proteins
# (which there're non of, as stated)
assert agg_complete_proteins_inclusion_df.isna().any(axis=1).any() == False

agg_complete_proteins_inclusion_df

# %%
# stats_agg_complete_proteins_inclusion_df = agg_complete_proteins_inclusion_df.groupby(["Chrom", "SwissProtName", "UnknownProbability", "NumOfReads"]).agg(
#     CoveredInXFractionsMean=("CoveredInXFractions", "mean"),
#     CoveredInXFractionsStd=("CoveredInXFractions", "std"),
# ).reset_index()

# stats_agg_complete_proteins_inclusion_df

# %%
agg_complete_proteins_inclusion_df.groupby(["Chrom", "SwissProtName", "NumOfReads"]).agg(
    CoveredInXFractionsMean=("CoveredInXFractions", "mean"),
    CoveredInXFractionsStd=("CoveredInXFractions", "std"),
).reset_index()

# %%
colors = ["red", "blue", "green"]

swiss_prot_color_dict = {
    swiss_prot_name.split("_")[0]: color
    for swiss_prot_name, color in zip(swiss_prot_name_by_chrom.values(), colors)
}

fig = px.scatter(
    agg_complete_proteins_inclusion_df.groupby(["Chrom", "SwissProtName", "NumOfReads"]).agg(
        CoveredInXFractionsMean=("CoveredInXFractions", "mean"),
        CoveredInXFractionsStd=("CoveredInXFractions", "std"),
    ).reset_index(),
    x="NumOfReads",
    y="CoveredInXFractionsMean",
    error_y="CoveredInXFractionsStd",
    color="SwissProtName",
    color_discrete_map=swiss_prot_color_dict,
    log_x=True,
    opacity=0.5,
    facet_col="SwissProtName",
    # facet_row="UnknownProbability",
)

# width = 700
width = 900
height = 300
# height = 800

fig.update_yaxes(range=[0, 6], dtick=1)

fig.update_layout(
    width=width,
    height=height,
    showlegend=False,
)

# fig.write_image(
#     "Graph assessment - mean % of distinct proteins recovered.svg",
#     width=width,
#     height=height,
# )

fig.show()

# %%
colors = ["red", "blue", "green"]

swiss_prot_color_dict = {
    swiss_prot_name.split("_")[0]: color
    for swiss_prot_name, color in zip(swiss_prot_name_by_chrom.values(), colors)
}

fig = px.scatter(
    agg_complete_proteins_inclusion_df.groupby(["Chrom", "SwissProtName", "UnknownProbability", "NumOfReads"]).agg(
        CoveredInXFractionsMean=("CoveredInXFractions", "mean"),
        CoveredInXFractionsStd=("CoveredInXFractions", "std"),
    ).reset_index(),
    x="NumOfReads",
    y="CoveredInXFractionsMean",
    error_y="CoveredInXFractionsStd",
    color="SwissProtName",
    color_discrete_map=swiss_prot_color_dict,
    log_x=True,
    opacity=0.5,
    facet_col="SwissProtName",
    facet_row="UnknownProbability",
)

# width = 700
width = 900
# height = 400
height = 800

fig.update_yaxes(range=[0, 6], dtick=1)

fig.update_layout(
    width=width,
    height=height,
    showlegend=False,
)

# fig.write_image(
#     "Graph assessment - mean % of distinct proteins recovered.svg",
#     width=width,
#     height=height,
# )

fig.show()

# %%
colors = ["red", "blue", "green"]

swiss_prot_color_dict = {
    swiss_prot_name.split("_")[0]: color
    for swiss_prot_name, color in zip(swiss_prot_name_by_chrom.values(), colors)
}


fig = make_subplots(
    rows=1,
    cols=1,
    print_grid=False,
    x_title="Gene",
    y_title="Known sites' editing levels [%]",
)

for chrom, swiss_prot_name in swiss_prot_name_by_chrom.items():
    df = agg_complete_proteins_inclusion_df.loc[(agg_complete_proteins_inclusion_df["Chrom"] == chrom)]
    x = df["SwissProtName"]
    y = df["%EditingLevel"]

    color = chrom_color_dict[chrom]
    swiss_prot_name = swiss_prot_name.split("_")[0]
    # ic(color, dash)
    fig.add_trace(
        # go.Box(
        go.Violin(
            x=x,
            y=y,
            name=swiss_prot_name,
            # line=dict(color=color,),
            line=dict(
                color="black",
            ),
            # line=dict(color="grey",),
            # line=dict(color=css_color_with_opacity("black", 70)),
            marker=dict(
                color=color,
            ),
            fillcolor=css_color_with_opacity(color, 30),
            # opacity=0.7,
            points="all",
            box_visible=True,
            meanline_visible=True,
        )
    )

width = 700
# height = width * 600 / 1000
height = 400

# fig.update_yaxes(tick0=0, dtick=10)

fig.update_layout(
    width=width,
    height=height,
    showlegend=False,
)

# fig.write_image(
#     "Graph assessment - mean % of distinct proteins recovered.svg",
#     width=width,
#     height=height,
# )

fig.show()

# %%

# %%

# %%

# %%

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


# grouped_fn_df = (
#     fn_df.groupby(["Chrom", "UnknownProbability", "Fraction"])
#     .agg({base_y_col: ["mean", "std"]})
#     .reset_index()
# )
# max_y_axis = (
#     grouped_fn_df[(base_y_col, "mean")].add(grouped_fn_df[(base_y_col, "std")]).max()
#     * 1.05
# )
# fig.update_yaxes(
#     range=[0, max_y_axis],
#     # type="log",
#     # range=[-1, np.log(max_y_axis)/np.log(10)]
# )

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

# %% [markdown]
# # True positives

# %%
# concat_validated_false_positives_df

# %%
concat_validated_false_positives_df.loc[
    ~concat_validated_false_positives_df["IsFalsePositive"]
]

# %%
concat_validated_true_positives_df = (
    concat_validated_false_positives_df.loc[
        ~concat_validated_false_positives_df["IsFalsePositive"]
    ]
    .groupby(
        [
            "Chrom",
            "SwissProtName",
            "UnknownProbability",
            "Fraction",
            "NumOfReads",
            "DataCreationRepetition",
        ]
    )
    .size()
    .reset_index(name="TruePositivePartiallyUnknownMaxNumOfDistinctProteins")
)
concat_validated_true_positives_df

# %%
wide_max_per_fraction_per_data_creation_merged_distinct_df

# %%
true_positives_wide_max_per_fraction_per_data_creation_merged_distinct_df = (
    wide_max_per_fraction_per_data_creation_merged_distinct_df.merge(
        concat_validated_true_positives_df
    )
)
true_positives_wide_max_per_fraction_per_data_creation_merged_distinct_df[
    "%TruePositiveDistinctProteinsRecovered"
] = (
    true_positives_wide_max_per_fraction_per_data_creation_merged_distinct_df[
        "TruePositivePartiallyUnknownMaxNumOfDistinctProteins"
    ]
    .mul(100)
    .div(
        true_positives_wide_max_per_fraction_per_data_creation_merged_distinct_df[
            "CompleteMaxNumOfDistinctProteins"
        ]
    )
    .round(2)
)
true_positives_wide_max_per_fraction_per_data_creation_merged_distinct_df

# %%
(
    true_positives_wide_max_per_fraction_per_data_creation_merged_distinct_df.loc[
:, ["Chrom", "SwissProtName", "UnknownProbability", "DataCreationRepetition", "Fraction", "NumOfReads", "CompleteMaxNumOfDistinctProteins", "TruePositivePartiallyUnknownMaxNumOfDistinctProteins", "%TruePositiveDistinctProteinsRecovered"]
]
    .rename(columns={"CompleteMaxNumOfDistinctProteins": "NumOfDistinctCompleteProteins", "TruePositivePartiallyUnknownMaxNumOfDistinctProteins": "NumOfTrueCoveringDistinctErroredPartiallyUnknownProteins"})
    .sort_values(["Chrom", "SwissProtName", "UnknownProbability", "DataCreationRepetition", "Fraction",], ignore_index=True)
)["%TruePositiveDistinctProteinsRecovered"]

# %%
(
    fn_df
    .drop(columns=["DistinctCompleteProteins", "AreDistinctCompleteProteinsCovered", "NumOfTrueDistinctErroredPartiallyUnknownProteins/NumOfCoveredDistinctCompleteProteins"])
    .sort_values(["Chrom", "SwissProtName", "UnknownProbability", "DataCreationRepetition", "Fraction",], ignore_index=True)
    .loc[:, ['Chrom', 'SwissProtName', 'UnknownProbability',
       'DataCreationRepetition', 'Fraction', 'NumOfReads',
       'NumOfDistinctCompleteProteins',  'NumOfTrueCoveringDistinctErroredPartiallyUnknownProteins', 'NumOfCoveredDistinctCompleteProteins',    'NumOfUncoveredDistinctCompleteProteins', '%UncoveredProteins', '%CoveredProteins']]
)["%CoveredProteins"]

# %%
old_cov_prct = Out[685]
new_cov_prct = Out[686]

# %%
new_cov_prct.sub(old_cov_prct).describe()

# %%
# tp_fn_merged_dfs = true_positives_wide_max_per_fraction_per_data_creation_merged_distinct_df.merge(fn_df)

# assert tp_fn_merged_dfs["CompleteMaxNumOfDistinctProteins"].eq(tp_fn_merged_dfs["NumOfDistinctCompleteProteins"]).all()
# assert tp_fn_merged_dfs["TruePositivePartiallyUnknownMaxNumOfDistinctProteins"].eq(tp_fn_merged_dfs["NumOfTrueCoveringDistinctErroredPartiallyUnknownProteins"]).all()

# tp_fn_merged_dfs

# %%

# %%
# wide_max_per_fraction_per_data_creation_merged_distinct_df.head(20)

# %%
# # stats over data creation repetitions

# stats_of_true_positives_max_per_fraction_merged_distinct_df = (
#     true_positives_wide_max_per_fraction_per_data_creation_merged_distinct_df.groupby(
#         [
#             "Chrom",
#             "SwissProtName",  # same as Chrom
#             "UnknownProbability",
#             "Fraction",
#             "NumOfReads",  # same as Fraction
#         ]
#     )
#     .agg(
#         MeanPrctDistinctProteinsRecovered=pd.NamedAgg(
#             column="%TruePositiveDistinctProteinsRecovered", aggfunc="mean"
#         ),
#         StdPrctDistinctProteinsRecovered=pd.NamedAgg(
#             column="%TruePositiveDistinctProteinsRecovered", aggfunc="std"
#         ),
#     )
#     .reset_index()
#     # .rename(columns={"MeanPrctDistinctProteinsRecovered": "MeanPrctDistinctProteinsRecovered"})
# )
# stats_of_true_positives_max_per_fraction_merged_distinct_df

# %%

# %%
# stats_of_fn_df = (
#     fn_df
#     .groupby([
#             "Chrom",
#             "SwissProtName",  # same as Chrom
#             "UnknownProbability",
#             "Fraction",
#             "NumOfReads",  # same as Fraction
#         ])
#     .agg(
#         MeanPrctDistinctProteinsUncovered=pd.NamedAgg(
#             column="%UncoveredProteins", aggfunc="mean"
#         ),
#         StdPrctDistinctProteinsUncovered=pd.NamedAgg(
#             column="%UncoveredProteins", aggfunc="std"
#         ),
#     )
#     .reset_index()
# )

# tp_fn_merged_stats_df = stats_of_true_positives_max_per_fraction_merged_distinct_df.merge(
#     stats_of_fn_df
# )
# tp_fn_merged_stats_df

# %%
# tp_fn_merged_stats_df["MeanPrctDistinctProteinsRecovered"].add(tp_fn_merged_stats_df["MeanPrctDistinctProteinsUncovered"])

# %%

# %%
# stats_of_true_positives_max_per_fraction_merged_distinct_df[
#     "MeanPrctDistinctProteinsRecovered"
# ].describe()

# %%
# dashes = ["solid", "dash", "dot"]
# colors = ["red", "blue", "green"]

# # unknown_probabilities_dash_dict = {unknown_probability: color for unknown_probability, color in zip(unknown_probabilities, ["red", "blue", "green"])}
# # chrom_color_dict = {chrom: dash for chrom, dash in zip(swiss_prot_name_by_chrom, ['dash', 'dot', 'solid'])}

# unknown_probabilities_dash_dict = {
#     unknown_probability: dash
#     for unknown_probability, dash in zip(unknown_probabilities, dashes)
# }
# chrom_color_dict = {
#     chrom: color for chrom, color in zip(swiss_prot_name_by_chrom, colors)
# }

# fig = make_subplots(
#     rows=1,
#     cols=1,
#     print_grid=False,
#     x_title="Simulated reads",
#     y_title="Mean % of distinct proteins recovered",
# )

# for unknown_probability in unknown_probabilities:
#     for chrom, swiss_prot_name in swiss_prot_name_by_chrom.items():
#         df = stats_of_true_positives_max_per_fraction_merged_distinct_df.loc[
#             (
#                 stats_of_true_positives_max_per_fraction_merged_distinct_df["Chrom"]
#                 == chrom
#             )
#             & (
#                 stats_of_true_positives_max_per_fraction_merged_distinct_df[
#                     "UnknownProbability"
#                 ]
#                 == unknown_probability
#             )
#         ]
#         x = df["NumOfReads"]
#         y = df["MeanPrctDistinctProteinsRecovered"]
#         error_y = df["StdPrctDistinctProteinsRecovered"]

#         color = chrom_color_dict[chrom]
#         dash = unknown_probabilities_dash_dict[unknown_probability]
#         swiss_prot_name = swiss_prot_name.split("_")[0]
#         # ic(color, dash)
#         fig.add_trace(
#             go.Scatter(
#                 x=x,
#                 y=y,
#                 error_y=dict(
#                     type="data",  # value of error bar given in data coordinates
#                     array=error_y,
#                     visible=True,
#                 ),
#                 # legendgroup=swiss_prot_name,  # this can be any string
#                 # legendgrouptitle_text=swiss_prot_name,
#                 # name=unknown_probability,
#                 legendgroup=unknown_probability,  # this can be any string
#                 legendgrouptitle_text=unknown_probability,
#                 name=swiss_prot_name,
#                 line=dict(color=color, dash=dash),
#                 opacity=0.7,
#             )
#         )

# width = 700
# # height = width * 600 / 1000
# height = 550

# fig.update_layout(
#     width=width,
#     height=height,
#     legend=dict(
#         # title="Swiss prot name, unknown probability",
#         # title="Unknown probability, gene",
#         # title="Unknown probability, mock gene",
#         title="NA probability, mock gene       ",
#     ),
#     margin_r=140,
# )

# fig.write_image(
#     "Graph assessment - mean % of distinct proteins recovered.svg",
#     width=width,
#     height=height,
# )

# fig.show()

# %% [markdown]
# # Merged plots

# %% jupyter={"source_hidden": true}
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

# fig.write_image(
#     "Graph assessment - mean % of false-positive distinct proteins.svg",
#     width=width,
#     height=height,
# )

fig.show()

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

# fig.write_image(
#     "Graph assessment - mean % of uncovered distinct proteins.svg",
#     width=width,
#     height=height,
# )

fig.show()

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
    y_title="Mean % of covered distinct proteins",
)

for unknown_probability in unknown_probabilities:
    for chrom, swiss_prot_name in swiss_prot_name_by_chrom.items():
        df = fn_df.loc[
            (fn_df["Chrom"] == chrom)
            & (fn_df["UnknownProbability"] == unknown_probability)
        ]
        df = (
            df.groupby("NumOfReads")
            .agg({"%CoveredProteins": ["mean", "std"]})
            .reset_index()
        )
        # df

        x = df["NumOfReads"]
        y = df[("%CoveredProteins", "mean")]
        error_y = df[("%CoveredProteins", "std")]

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


# grouped_fn_df = (
#     fn_df.groupby(["Chrom", "UnknownProbability", "Fraction"])
#     .agg({"%CoveredProteins": ["mean", "std"]})
#     .reset_index()
# )
# max_y_axis = (
#     grouped_fn_df[("%CoveredProteins", "mean")]
#     .add(grouped_fn_df[("%CoveredProteins", "std")])
#     .max()
#     * 1.05
# )
# fig.update_yaxes(range=[0, max_y_axis])

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

# fig.write_image(
#     "Graph assessment - mean % of uncovered distinct proteins.svg",
#     width=width,
#     height=height,
# )

fig.show()

# %% jupyter={"source_hidden": true}
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
        df = stats_of_true_positives_max_per_fraction_merged_distinct_df.loc[
            (
                stats_of_true_positives_max_per_fraction_merged_distinct_df["Chrom"]
                == chrom
            )
            & (
                stats_of_true_positives_max_per_fraction_merged_distinct_df[
                    "UnknownProbability"
                ]
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

# fig.write_image(
#     "Graph assessment - mean % of distinct proteins recovered.svg",
#     width=width,
#     height=height,
# )

fig.show()

# %% jupyter={"source_hidden": true}
# dashes = ["solid", "dash", "dot"]
# colors = ["red", "blue", "green"]

# # unknown_probabilities_dash_dict = {unknown_probability: color for unknown_probability, color in zip(unknown_probabilities, ["red", "blue", "green"])}
# # chrom_color_dict = {chrom: dash for chrom, dash in zip(swiss_prot_name_by_chrom, ['dash', 'dot', 'solid'])}

# unknown_probabilities_dash_dict = {
#     unknown_probability: dash
#     for unknown_probability, dash in zip(unknown_probabilities, dashes)
# }
# chrom_color_dict = {
#     chrom: color for chrom, color in zip(swiss_prot_name_by_chrom, colors)
# }
# # y_titles = [
# #     "Mean % of distinct proteins recovered",
# #     "Mean % of false-positive distinct proteins",
# #     "Mean % of uncovered distinct proteins",
# # ]
# y_titles = [
#     "Mean % of distinct<br>proteins recovered",
#     "Mean % of false-positive<br>distinct proteins",
#     "Mean % of uncovered<br>distinct proteins",
# ]

# rows = 3

# fig = make_subplots(
#     rows=rows,
#     cols=1,
#     print_grid=False,
#     x_title="Simulated reads",
#     shared_xaxes="all",
#     vertical_spacing=0.03,
# )

# for row, y_title in zip(range(1, rows + 1), y_titles):

#     # max_y_axis = 0

#     for unknown_probability in unknown_probabilities:
#         for chrom, swiss_prot_name in swiss_prot_name_by_chrom.items():

#             if row == 1:
#                 df = stats_of_true_positives_max_per_fraction_merged_distinct_df.loc[
#                     (
#                         stats_of_true_positives_max_per_fraction_merged_distinct_df[
#                             "Chrom"
#                         ]
#                         == chrom
#                     )
#                     & (
#                         stats_of_true_positives_max_per_fraction_merged_distinct_df[
#                             "UnknownProbability"
#                         ]
#                         == unknown_probability
#                     )
#                 ]
#                 x = df["NumOfReads"]
#                 y = df["MeanPrctDistinctProteinsRecovered"]
#                 error_y = df["StdPrctDistinctProteinsRecovered"]

#             elif row == 2:
#                 df = concat_validated_false_positives_prct_df.loc[
#                     (concat_validated_false_positives_prct_df["Chrom"] == chrom)
#                     & (
#                         concat_validated_false_positives_prct_df["UnknownProbability"]
#                         == unknown_probability
#                     )
#                 ]
#                 x = df["NumOfReads"]
#                 y = df["MeanPrctOfPalsePositives"]
#                 error_y = df["STDOfMeanPrctOfPalsePositives"]

#             elif row == 3:
#                 df = fn_df.loc[
#                     (fn_df["Chrom"] == chrom)
#                     & (fn_df["UnknownProbability"] == unknown_probability)
#                 ]
#                 df = (
#                     df.groupby("NumOfReads")
#                     .agg({"%UncoveredProteins": ["mean", "std"]})
#                     .reset_index()
#                 )
#                 x = df["NumOfReads"]
#                 y = df[("%UncoveredProteins", "mean")]
#                 error_y = df[("%UncoveredProteins", "std")]

#             else:
#                 raise ValueError(f"{row=} not programed for this plot")

#             color = chrom_color_dict[chrom]
#             dash = unknown_probabilities_dash_dict[unknown_probability]
#             swiss_prot_name = swiss_prot_name.split("_")[0]

#             # max_y_axis = max(max_y_axis, y.add(error_y).max())

#             if row == 1:
#                 fig.add_trace(
#                     go.Scatter(
#                         x=x,
#                         y=y,
#                         error_y=dict(
#                             type="data",  # value of error bar given in data coordinates
#                             array=error_y,
#                             visible=True,
#                         ),
#                         legendgroup=unknown_probability,  # this can be any string
#                         legendgrouptitle_text=unknown_probability,
#                         name=swiss_prot_name,
#                         line=dict(color=color, dash=dash),
#                         opacity=0.7,
#                     ),
#                     row=row,
#                     col=1,
#                 )
#             else:
#                 fig.add_trace(
#                     go.Scatter(
#                         x=x,
#                         y=y,
#                         error_y=dict(
#                             type="data",  # value of error bar given in data coordinates
#                             array=error_y,
#                             visible=True,
#                         ),
#                         showlegend=False,
#                         line=dict(color=color, dash=dash),
#                         opacity=0.7,
#                     ),
#                     row=row,
#                     col=1,
#                 )

#     fig.update_yaxes(
#         title_text=y_title,
#         row=row,
#         col=1,
#         # range=[0, max_y_axis*1.05]
#     )

# width = 600
# # height = width * 600 / 1000
# height = 300 * rows

# fig.update_layout(
#     width=width,
#     height=height,
#     legend=dict(
#         # title="Swiss prot name, unknown probability",
#         # title="Unknown probability, gene",
#         # title="Unknown probability, mock gene",
#         title="NA probability, mock gene       ",
#     ),
# )

# fig.write_image(
#     "Graph assessment - merged plots.svg",
#     width=width,
#     height=height,
# )

# fig.show()

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

# y_titles = [
#     "Mean % of covered<br>distinct proteins",
#     "Mean % of false-positive<br>distinct proteins",
# ]
# y_titles = [
#     "Mean % of covered distinct proteins",
#     "Mean % of false-positive distinct proteins",
# ]
y_titles = [
    "Mean % of covered isoforms",
    "Mean % of false-positive isoforms",
]

rows = 2

fig = make_subplots(
    rows=rows,
    cols=1,
    print_grid=False,
    x_title="Simulated reads",
    shared_xaxes="all",
    vertical_spacing=0.05,
)

for row, y_title in zip(range(1, rows + 1), y_titles):

    # max_y_axis = 0

    for unknown_probability in unknown_probabilities:
        for chrom, swiss_prot_name in swiss_prot_name_by_chrom.items():

            if row == 1:
                df = fn_df.loc[
                    (fn_df["Chrom"] == chrom)
                    & (fn_df["UnknownProbability"] == unknown_probability)
                ]
                df = (
                    df.groupby("NumOfReads")
                    .agg({"%CoveredProteins": ["mean", "std"]})
                    .reset_index()
                )
                # df
        
                x = df["NumOfReads"]
                y = df[("%CoveredProteins", "mean")]
                error_y = df[("%CoveredProteins", "std")]

            elif row == 2:
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

            # elif row == 3:
            #     df = fn_df.loc[
            #         (fn_df["Chrom"] == chrom)
            #         & (fn_df["UnknownProbability"] == unknown_probability)
            #     ]
            #     df = (
            #         df.groupby("NumOfReads")
            #         .agg({"%UncoveredProteins": ["mean", "std"]})
            #         .reset_index()
            #     )
            #     x = df["NumOfReads"]
            #     y = df[("%UncoveredProteins", "mean")]
            #     error_y = df[("%UncoveredProteins", "std")]

            else:
                raise ValueError(f"{row=} not programed for this plot")

            color = chrom_color_dict[chrom]
            dash = unknown_probabilities_dash_dict[unknown_probability]
            swiss_prot_name = swiss_prot_name.split("_")[0]

            # max_y_axis = max(max_y_axis, y.add(error_y).max())

            if row == 1:
                fig.add_trace(
                    go.Scatter(
                        x=x,
                        y=y,
                        error_y=dict(
                            type="data",  # value of error bar given in data coordinates
                            array=error_y,
                            visible=True,
                        ),
                        legendgroup=unknown_probability,  # this can be any string
                        legendgrouptitle_text=unknown_probability,
                        name=swiss_prot_name,
                        line=dict(color=color, dash=dash),
                        opacity=0.7,
                    ),
                    row=row,
                    col=1,
                )
            else:
                fig.add_trace(
                    go.Scatter(
                        x=x,
                        y=y,
                        error_y=dict(
                            type="data",  # value of error bar given in data coordinates
                            array=error_y,
                            visible=True,
                        ),
                        showlegend=False,
                        line=dict(color=color, dash=dash),
                        opacity=0.7,
                    ),
                    row=row,
                    col=1,
                )

    fig.update_yaxes(
        title_text=y_title,
        row=row,
        col=1,
        # range=[0, max_y_axis*1.05]
    )

width = 650
# height = width * 600 / 1000
height = 400 * rows

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
    "Graph assessment - merged plots.svg",
    width=width,
    height=height,
)

fig.show()
