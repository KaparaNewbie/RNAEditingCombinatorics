# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.4
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Prep

# %%
from pathlib import Path
from itertools import product

import pandas as pd
import numpy as np
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from icecream import ic

# %%
pacbio_platforms = (
    ["PacBio1"] * 2 +
    ["PacBio2"] * 2 +
    ["PacBio3"] * 3
)
pacbio_samples = [
	"GRIA2", "PCLO",
	"ADAR1", "IQEC1",
	"GRIA2", "ADAR1", "IQEC1",
]

# %%
shortened_platform_names = {
    "PacBio1": "PB1",
    "PacBio2": "PB2",
    "PacBio3": "PB3",
    "Illumina1": "IL1"
}

# %%
illumina_samples = [
	"RUSC2",
	"TRIM2",
	"CA2D3",
	"ABL",
	"DGLA",
	"K0513",
	"KCNAS",
	"ACHA4",
	"ANR17",
	"TWK7",
	"SCN1",
	"CACB2",
	"RIMS2",
	"PCLO",
	"DOP1",
	"IQEC1",
	"CSKI1",
	"MTUS2",
	"ROBO2"
]
illumina_platforms = ["Illumina1"] * len(illumina_samples)

# %%
unique_proteins_first_col_pos = 14

pacbio_unique_proteins_files = [
	# PacBio1
	"D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz",
	"D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz",
	# PacBio2
	"D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples/ADAR1.Merged.r64296e203404D01.aligned.sorted.MinRQ998.unique_proteins.csv.gz",
	"D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples/IQEC.Merged.r64296e203404D01.aligned.sorted.MinRQ998.unique_proteins.csv.gz",
	# PacBio3
	"D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads/comp141693_c0_seq1.merged.MinRQ998.unique_proteins.csv.gz",
	"D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads/comp134400_c0_seq1_extended.merged.MinRQ998.unique_proteins.csv.gz",
	"D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads/comp141565_c6_seq3.merged.MinRQ998.unique_proteins.csv.gz",
]
pacbio_expression_files = [
	# PacBio1
	"D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/GRIA.DistinctUniqueProteins.ExpressionLevels.EntropyConsidered.csv",
	"D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO.DistinctUniqueProteins.ExpressionLevels.EntropyConsidered.csv",
	# PacBio2
	"D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples/ADAR1.DistinctUniqueProteins.ExpressionLevels.EntropyConsidered.csv",
	"D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples/IQEC1.DistinctUniqueProteins.ExpressionLevels.EntropyConsidered.csv",
	# PacBio3
	"D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads/GRIA2.DistinctUniqueProteins.ExpressionLevels.EntropyConsidered.csv",
	"D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads/ADAR1.DistinctUniqueProteins.ExpressionLevels.EntropyConsidered.csv",
	"D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads/IQEC1.DistinctUniqueProteins.ExpressionLevels.EntropyConsidered.csv",
]

illumina_chroms = [
    "comp141881_c0_seq3",
    "comp141044_c0_seq2",
    "comp140439_c0_seq1",
    "comp126362_c0_seq1",
    "comp141517_c0_seq1",
    "comp141840_c0_seq2",
    "comp141640_c0_seq1",
    "comp140987_c3_seq1",
    "comp140910_c2_seq1",
    "comp136058_c0_seq1",
    "comp141378_c0_seq7",
    "comp141158_c1_seq2",
    "comp140712_c0_seq3",
    "comp141882_c0_seq14",
    "comp141880_c1_seq3",
    "comp141565_c6_seq3",
    "comp141684_c0_seq1",
    "comp141532_c3_seq11",
    "comp141574_c0_seq3",
]
illumina_unique_proteins_files = [
	f"/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.{chrom}.unique_proteins.csv"
    for chrom in illumina_chroms
]

illumina_expression_files = [
    Path(
        "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/FinalFixedExpressionFromCloud",
        sample,
        f"{sample}.DistinctUniqueProteins.ExpressionLevels.EntropyConsidered.csv"
	)
	for sample in illumina_samples
]
for exp_file in illumina_expression_files:
    assert exp_file.exists(), exp_file 

# %%
platforms = pacbio_platforms + illumina_platforms
samples = pacbio_samples + illumina_samples
unique_proteins_files = pacbio_unique_proteins_files + illumina_unique_proteins_files
expression_files = pacbio_expression_files + illumina_expression_files


# %%
def make_expression_df(expression_file):
    assignment_df = pd.read_csv(
        expression_file, 
        sep="\t", 
        dtype={
            "#Solution": str,
        },
        usecols = [
            'Gene', 'Protein', '#Solution', 'Fraction', 'SufficientEntropy', 'TotalWeightedSupportingReads'
        ]
    )

    assignment_df = assignment_df.loc[
        assignment_df["Fraction"].eq(1)
    ]

    solution_sizes = (
        assignment_df
        .groupby("#Solution")
        .size()
        .reset_index(name="Size")
    )
    max_solution_size = solution_sizes["Size"].max()
    max_solutions = solution_sizes.loc[
        solution_sizes["Size"].eq(max_solution_size),
        "#Solution"
    ].values
    if (len(max_solutions) > 1):
        ic(len(max_solutions), expression_file)
    rng = np.random.default_rng(1892)
    max_solution = rng.choice(max_solutions)

    assignment_df = assignment_df.loc[
        assignment_df["#Solution"].eq(max_solution)
    ]

    assignment_df = assignment_df.drop(columns=["#Solution", "Fraction"])

    assignment_df["%RelativeExpression"] = (
        100
        * assignment_df["TotalWeightedSupportingReads"]
        / assignment_df["TotalWeightedSupportingReads"].sum()
    )

    assignment_df = (
        assignment_df
        .sort_values("TotalWeightedSupportingReads")
        # .reset_index(drop=True)
    )
    assignment_df["RareToCommonProteinIndex"] = list(range(1, assignment_df.shape[0] + 1))


    assignment_df = (
        assignment_df
        .sort_values("TotalWeightedSupportingReads", ascending=False)
        .reset_index(drop=True)
    )
    assignment_df["CommonToRareProteinIndex"] = list(range(1, assignment_df.shape[0] + 1))
    # assignment_df["%CommonToRareCummulativeRelativeExpression"] = assignment_df[
    #     "%RelativeExpression"
    # ].cumsum()

    # assignment_df.insert(0, "Chrom", chrom)

    return assignment_df

# %%
X_common_proteins = [0.01, 0.03, 0.05]
Y_rare_proteins = [0.1, 0.2, 0.3]

x_and_y_proteins_denote_fractions = True

# %%
out_dir = Path("/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/JointChimericReadsAnalysis")


# %%
def define_out_file(out_dir, platform, sample_name, x_common_proteins, y_rare_proteins):
    out_file = Path(
		out_dir,
		f"{platform}.{sample_name}.X_{x_common_proteins}.Y_{y_rare_proteins}.csv.gz"
	)
    if not out_file.exists():
        raise FileNotFoundError(f"Expected output file {out_file} does not exist.")
    return out_file

out_files = [
	define_out_file(out_dir, platform, sample_name, x_common_proteins, y_rare_proteins)
	for platform, sample_name in zip(
        platforms, samples
    )
    for x_common_proteins, y_rare_proteins in product(
        X_common_proteins, Y_rare_proteins
    )
]
len(out_files)


# %%
def parse_pairs_str(s, out_delim=";", in_delim=","):
    if s in ["", "None", np.nan]:
        return []
    pairs = s.split(out_delim)
    return [pair.split(in_delim) for pair in pairs]


# %%
def parse_int_pairs_str(s, out_delim=";", in_delim=","):
    if s in ["", "None", np.nan]:
        return []
    pairs = s.split(out_delim)
    split_pairs = [pair.split(in_delim) for pair in pairs]
    return [
        [int(x) for x in pair]
        for pair in split_pairs
    ]


# %%
num_of_aas_per_platform_and_sample = {
    (platform, sample): (
        pd.read_table(unique_proteins_file, nrows=5)
        .iloc[:, unique_proteins_first_col_pos:]
        .shape[1]
    )
    for platform, sample, unique_proteins_file in zip(platforms, samples, unique_proteins_files)
}

num_of_aas_per_platform_and_sample_df = (
    pd.DataFrame(num_of_aas_per_platform_and_sample, index=["Platform", "Sample", "TotalNumOfAAs"])
    .T
    .reset_index()
    .drop(columns=["Platform", "Sample"])
    .rename(columns={"level_0": "Platform", "level_1": "Sample"})
)
num_of_aas_per_platform_and_sample_df


# %%
def read_results_file(out_file):
    df = pd.read_csv(
        out_file, 
        sep="\t",
        dtype={
            "ChimericUniqueProteinsCombinations": str,
            "ChimericProtsIntersectingAAsIndices": str
        }
    )
    df = df.drop(columns="ChimericProtsIndices")
    df["ChimericUniqueProteinsCombinations"] = df["ChimericUniqueProteinsCombinations"].apply(parse_pairs_str)
    df["ChimericProtsIntersectingAAsIndices"] = df["ChimericProtsIntersectingAAsIndices"].apply(parse_int_pairs_str)
    return df


# %%
result_df = pd.concat(
    [read_results_file(out_file) for out_file in out_files],
    ignore_index=True
)
result_df = num_of_aas_per_platform_and_sample_df.merge(
    result_df,
    how="right"
)
result_df.insert(
    1,
    "ShortPlatform",
    result_df["Platform"].map(shortened_platform_names)
)

result_df


# %%
def make_x_y_expression_df(
    platforms,
    samples,
    expression_files,
    X_common_proteins,
    Y_rare_proteins,
    stats_df
):

    results = []

    for platform, sample, expression_file in zip(
        platforms, samples, expression_files
    ):
        expression_df = make_expression_df(expression_file)
        
        for x_common_proteins, y_rare_proteins in product(X_common_proteins, Y_rare_proteins):
            
            actual_x_common_proteins, actual_y_rare_proteins = stats_df.loc[
                (stats_df["Platform"].eq(platform))
                & (stats_df["Sample"].eq(sample))
                & (stats_df["XCommonProteins"].eq(x_common_proteins))
                & (stats_df["YRareProteins"].eq(y_rare_proteins)),
                ["ActualXCommonProteins", "ActualYRareProteins"]
            ].squeeze()
            
            x_common_proteins_cumulative_expression = expression_df.loc[
                expression_df["CommonToRareProteinIndex"].le(actual_x_common_proteins),
                "%RelativeExpression"
            ].sum()
            y_rare_proteins_cumulative_expression = expression_df.loc[
                expression_df["RareToCommonProteinIndex"].le(actual_y_rare_proteins),
                "%RelativeExpression"
            ].sum()
            
            result = [
                platform, sample, x_common_proteins, y_rare_proteins, 
                x_common_proteins_cumulative_expression, y_rare_proteins_cumulative_expression
            ]
            
            results.append(result)
            
    x_y_expression_df = pd.DataFrame(
        results,
        columns=[
            "Platform", "Sample", "XCommonProteins", "YRareProteins", 
            "XCommonProteinsCumulativeExpression", "YRareProteinsCumulativeExpression"
        ]
    )
    
    return x_y_expression_df


# %%
stats_df = (
    result_df
    .groupby(
        ["Platform", "Sample", "XYProteinsDenoteFractions", "XCommonProteins", "YRareProteins"],
        as_index=False
    )
    .agg(
        ActualXCommonProteins = ("ActualXCommonProteins", "first"),
        ActualYRareProteins = ("ActualYRareProteins", "first"),
        NumOfChimericProteinsOnAALevel = ("IsChimericOnAALevel", "sum"),
        TotalCommonUniqueReads = ("TotalCommonUniqueReads", "first"),
        TotalRareUniqueReads = ("TotalRareUniqueReads", "first"),
        TotalCommonReads = ("TotalCommonReads", "first"),
        TotalRareReads = ("TotalRareReads", "first"),
    )
)
stats_df.insert(
    stats_df.columns.get_loc("NumOfChimericProteinsOnAALevel") + 1,
    "%OfChimericProteinsOnAALevel",
    stats_df["NumOfChimericProteinsOnAALevel"].mul(100).div(stats_df["ActualYRareProteins"])
)

stats_df.insert(
    1,
    "ShortPlatform",
    stats_df["Platform"].map(shortened_platform_names)
)

stats_df.insert(
    3,
    "Platform-Sample",
    stats_df["ShortPlatform"] + "-" + stats_df["Sample"]
)

x_y_expression_df = make_x_y_expression_df(
    platforms,
    samples,
    expression_files,
    X_common_proteins,
    Y_rare_proteins,
    stats_df
)
stats_df = stats_df.merge(
    x_y_expression_df
)

stats_df

# %%
x_y_expression_df.loc[
    (x_y_expression_df["Platform"].isin(set(pacbio_platforms)))
    & (x_y_expression_df["XCommonProteins"].eq(0.01)),
    "XCommonProteinsCumulativeExpression"
].describe()

# %%
platform_colormap = {
    platform: color
    for platform, color in zip(
        sorted(set(platforms)),
        # px.colors.qualitative.Plotly
        px.colors.qualitative.D3
        # px.colors.qualitative.T10
    )
}
platform_colormap

# %% [markdown]
# # % of chimeric proteins

# %%
# def make_subplot_title(x_common_proteins, y_rare_proteins, x_and_y_proteins_denote_fractions):
#     if x_and_y_proteins_denote_fractions:
#         x_label = f"Common: {x_common_proteins * 100:.0f}%"
#         y_label = f"Rare: {y_rare_proteins * 100:.0f}%"
#     else:
#         x_label = f"Common: {x_common_proteins}"
#         y_label = f"Rare: {y_rare_proteins}"
#     return f"{x_label}<br>{y_label}"

# subplot_titles = [
#     make_subplot_title(x_common_proteins, y_rare_proteins, x_and_y_proteins_denote_fractions)
#     for x_common_proteins, y_rare_proteins in product(
#         X_common_proteins, Y_rare_proteins
#     )
# ]
# subplot_titles

# %%
row_titles = [
    f"Common: {x_common_proteins * 100:.0f}%"
    for x_common_proteins in X_common_proteins
]
col_titles = [
    f"Rare: {y_rare_proteins * 100:.0f}%"
    for y_rare_proteins in Y_rare_proteins
]
row_titles, col_titles

# %%
fig = make_subplots(
    rows=len(X_common_proteins), 
    cols=len(Y_rare_proteins), 
    shared_xaxes="all", 
    shared_yaxes="all",
    y_title="% of chimeric proteins on AA level",
    # subplot_titles=subplot_titles,
    column_titles=col_titles,
    row_titles=row_titles,
    vertical_spacing=0.2/len(X_common_proteins),
    horizontal_spacing=0.1/len(Y_rare_proteins),
)

for i, x_common_proteins in enumerate(X_common_proteins, start=1):
    for j, y_rare_proteins in enumerate(Y_rare_proteins, start=1):
        subset_df = stats_df[
            (stats_df["XCommonProteins"] == x_common_proteins) &
            (stats_df["YRareProteins"] == y_rare_proteins)
        ]
        if subset_df.shape[0] < 26:
            ic(x_common_proteins, y_rare_proteins, subset_df.shape[0])
        fig.add_trace(
            go.Bar(
                x=subset_df.loc[:, ["ShortPlatform", "Sample"]].T.values.tolist(),
                y=subset_df["%OfChimericProteinsOnAALevel"],
                marker=dict(
                    color=subset_df["Platform"].map(platform_colormap)
                ),
            ),
            row=i,
            col=j
        )
fig.update_xaxes(
    tickfont_size=10
)
fig.update_yaxes(
    dtick=20
)

fig.update_layout(
    template="plotly_white",
    showlegend=False,
    width=1400,
    height=800
)
        
fig.show()

# %%
# ...existing code...

metric_colors = {
    "XCommonProteinsCumulativeExpression": "#1f77b4",
    "YRareProteinsCumulativeExpression":   "#ff7f0e",
}

fig = make_subplots(
    rows=len(X_common_proteins),
    cols=len(Y_rare_proteins),
    shared_xaxes="all",
    shared_yaxes="all",
    y_title="Cumulative expression [%]",
    column_titles=col_titles,
    row_titles=row_titles,
    vertical_spacing=0.2 / len(X_common_proteins),
    horizontal_spacing=0.1 / len(Y_rare_proteins),
)

for i, x_common_proteins in enumerate(X_common_proteins, start=1):
    for j, y_rare_proteins in enumerate(Y_rare_proteins, start=1):
        subset_df = stats_df[
            (stats_df["XCommonProteins"] == x_common_proteins)
            & (stats_df["YRareProteins"] == y_rare_proteins)
        ].copy()

        # 1D categorical x axis
        subset_df["BoxGroup"] = subset_df["ShortPlatform"] + "-" + subset_df["Sample"]
        subset_df = subset_df.sort_values(["ShortPlatform", "Sample"])

        # Trace 1: common cumulative expression
        fig.add_trace(
            go.Bar(
                # x=subset_df["BoxGroup"],
                x=subset_df.loc[:, ["ShortPlatform", "Sample"]].T.values.tolist(),
                y=subset_df["XCommonProteinsCumulativeExpression"],
                name="Common cumulative expression",
                marker_color=metric_colors["XCommonProteinsCumulativeExpression"],
                # marker_pattern_shape="x",
                # marker=dict(
                #     color=subset_df["Platform"].map(platform_colormap)
                # ),
                offsetgroup="common",
                legendgroup="common",
                showlegend=(i == 1 and j == 1),
            ),
            row=i,
            col=j,
        )

        # Trace 2: rare cumulative expression
        fig.add_trace(
            go.Bar(
                # x=subset_df["BoxGroup"],
                x=subset_df.loc[:, ["ShortPlatform", "Sample"]].T.values.tolist(),
                y=subset_df["YRareProteinsCumulativeExpression"],
                name="Rare cumulative expression",
                marker_color=metric_colors["YRareProteinsCumulativeExpression"],
                # marker=dict(
                #     color=subset_df["Platform"].map(platform_colormap)
                # ),
                # marker_pattern_shape="/",
                offsetgroup="rare",
                legendgroup="rare",
                showlegend=(i == 1 and j == 1),
            ),
            row=i,
            col=j,
        )

fig.update_layout(
    template="plotly_white",
    width=1400,
    height=800,
    barmode="group",
    legend=dict(orientation="h", yanchor="top", y=-0.2, xanchor="left", x=0),
)

fig.update_xaxes(tickfont_size=10)
fig.update_yaxes(
    # dtick=20,
    dtick=10,
    range=[0, 100]
#     type="log"
)


fig.show()

# ...existing code...

# %%

# %%
# fig = px.scatter(
#     stats_df,
#     x="TotalCommonReads",
#     y="%OfChimericProteinsOnAALevel",
#     color="Platform",
#     hover_data=["Platform-Sample"],
#     color_discrete_map=platform_colormap,
#     # color="Platform-Sample"
#     log_x=True,
#     labels={
#         "TotalCommonReads": "Total common reads",
#         "%OfChimericProteinsOnAALevel": "% of chimeric proteins on AA level"
#     }
# )
# fig.update_layout(
#     template="plotly_white",
#     width=800,
#     height=600
# )
# fig.show()

# %%
# fig = px.scatter(
#     stats_df,
#     x="TotalCommonReads",
#     # color="XCommonProteins",
#     y="%OfChimericProteinsOnAALevel",
#     color="Platform",
#     facet_col="XCommonProteins",
#     hover_data=["Platform-Sample"],
#     # color_discrete_map=platform_colormap,
#     # color="Platform-Sample"
#     log_x=True,
#     labels={
#         "TotalCommonReads": "Total common reads",
#         "XCommonProteins": "% of common proteins",
#         "%OfChimericProteinsOnAALevel": "% of chimeric proteins on AA level"
#     }
# )
# fig.update_xaxes(dtick=1)
# fig.update_layout(
#     template="plotly_white",
#     width=1200,
#     height=600
# )
# fig.show()

# %%
# fig = px.scatter_3d(
#     stats_df,
#     x="TotalCommonReads",
#     y="TotalRareReads",
#     z="%OfChimericProteinsOnAALevel",
#     # z="Platform",
#     # color="%OfChimericProteinsOnAALevel",
#     color="Platform",
#     hover_data=["Platform", "Sample"],
#     color_discrete_map=platform_colormap,
#     # color="Platform-Sample"
#     log_x=True,
#     labels={
#         "TotalCommonReads": "Total common reads",
#         "TotalRareReads": "Total rare reads",
#         "%OfChimericProteinsOnAALevel": "% of chimeric proteins on AA level"
#     }
# )
# fig.update_layout(
#     template="plotly_white",
#     width=800,
#     height=800
# )
# fig.show()

# %%
# fig = px.scatter(
#     stats_df,
#     x="TotalCommonReads",
#     y="%OfChimericProteinsOnAALevel",
#     facet_col="Platform",
#     color="Sample",
#     color_discrete_sequence=px.colors.qualitative.Dark24,
#     log_x=True,
#     labels={
#         "TotalCommonReads": "Total common reads",
#         "%OfChimericProteinsOnAALevel": "% of chimeric proteins on AA level"
#     }
# )
# fig.update_xaxes(dtick=1)
# fig.update_layout(
#     template="plotly_white",
#     width=1500,
#     height=600
# )
# fig.show()

# %%
# fig = px.scatter(
#     stats_df,
#     x="TotalCommonReads",
#     y="TotalRareReads",
#     color="%OfChimericProteinsOnAALevel",
#     facet_col="Platform",
#     log_x=True,
#     labels={
#         "TotalCommonReads": "Total common reads",
#         "TotalRareReads": "Total rare reads",
#         "%OfChimericProteinsOnAALevel": "% of chimeric proteins<br>on AA level"
#     }
# )
# fig.update_xaxes(dtick=1)
# fig.update_layout(
#     template="plotly_white",
#     width=1400,
#     height=400
# )
# fig.show()

# %%
# fig = px.scatter(
#     stats_df,
#     x="TotalCommonReads",
#     y="TotalRareReads",
#     # color="%OfChimericProteinsOnAALevel",
#     color="Sample",
#     size="%OfChimericProteinsOnAALevel",
#     facet_col="Platform",
#     log_x=True,
#     labels={
#         "TotalCommonReads": "Total common reads",
#         "TotalRareReads": "Total rare reads",
#         "%OfChimericProteinsOnAALevel": "% of chimeric proteins<br>on AA level"
#     }
# )
# fig.update_xaxes(dtick=1)
# fig.update_layout(
#     template="plotly_white",
#     width=1600,
#     height=600
# )
# fig.show()

# %%
# # 1) Build a stable sample->symbol mapping (21 samples)
# samples_sorted = sorted(stats_df["Sample"].unique())

# # Plotly marker symbols (string names). Provide >= 21 distinct entries.
# symbol_sequence = [
#     "circle", "square", "diamond", "cross", "x",
#     "triangle-up", "triangle-down", "triangle-left", "triangle-right",
#     "pentagon", "hexagon", "hexagon2", "octagon",
#     "star", "hourglass", "bowtie",
#     "circle-open", "square-open", "diamond-open",
#     "triangle-up-open", "triangle-down-open",
# ]

# if len(samples_sorted) > len(symbol_sequence):
#     raise ValueError(
#         f"Need {len(samples_sorted)} symbols but only have {len(symbol_sequence)}. "
#         "Add more symbols to symbol_sequence."
#     )

# sample_symbol_map = dict(zip(samples_sorted, symbol_sequence))

# fig = px.scatter(
#     stats_df,
#     x="TotalCommonReads",
#     y="TotalRareReads",
#     color="%OfChimericProteinsOnAALevel",
#     symbol="Sample",
#     symbol_map=sample_symbol_map,  # <-- ensures unique symbol per sample
#     facet_col="Platform",
#     log_x=True,
#     labels={
#         "TotalCommonReads": "Total common reads",
#         "TotalRareReads": "Total rare reads",
#     },
    
# )

# fig.update_xaxes(dtick=1)

# # keep your non-overlapping legend/colorbar layout
# fig.update_layout(
#     template="plotly_white",
#     width=1600,
#     height=600,
#     legend=dict(
#         orientation="h",
#         yanchor="top",
#         y=-0.25,
#         xanchor="left",
#         x=0.0,
#         # title_text="Sample",
#         title_text="Gene",
#     ),
#     margin=dict(b=140, r=180),
# )
# fig.update_coloraxes(
#     colorbar=dict(
#         x=1.08,
#         y=0.5,
#         len=0.9,
#         title="% of chimeric proteins<br>on AA level",
#     )
# )

# fig.show()

# %% [markdown]
# # Intersection points

# %%
result_df

# %%
chimeric_result_df = (
    result_df.loc[
        result_df["IsChimericOnAALevel"],
        [
            "Platform", "ShortPlatform", "Sample", "TotalNumOfAAs", "XCommonProteins", "YRareProteins", "Protein", 
            "ChimericUniqueProteinsCombinations", "ChimericProtsIntersectingAAsIndices"
        ]
    ]
    .reset_index(drop=True)
    # .drop(columns="IsChimericOnAALevel")
)

chimeric_result_df


# %%
def compute_breakdown_points(b, e, total_aas):
    b_len = b # 1 is the first AA, as this is a julia result, so length is b, not b-1
    e_len = total_aas - e + 1

    # this breakdown point represent the midpoint of the shared region
    middle_breakdown_point = e + ((b - e) / 2)

    # this breakdown point represents the end of the longer region (or the end of the shared region if they are of equal length)
    if b_len >= e_len:
        longer_breakdown_point = b
    else:
        longer_breakdown_point = e
        
    normalized_middle_breakdown_point = 100 * middle_breakdown_point / total_aas
    normalized_longer_breakdown_point = 100 * longer_breakdown_point / total_aas

    return middle_breakdown_point, normalized_middle_breakdown_point, longer_breakdown_point, normalized_longer_breakdown_point


# %%
flat_chimeric_result_df = (
    chimeric_result_df
    # .drop(columns="NumOfChimericCombinations")
    .explode(
        ["ChimericUniqueProteinsCombinations", "ChimericProtsIntersectingAAsIndices"],
        ignore_index=True
    )
)

res = flat_chimeric_result_df.loc[
    :, ["ChimericProtsIntersectingAAsIndices", "TotalNumOfAAs"]
].apply(
    lambda row: compute_breakdown_points(*row["ChimericProtsIntersectingAAsIndices"], row["TotalNumOfAAs"]),
    axis=1,
    result_type="expand"
).rename(
    columns={
        0: "MiddleBreakdownPoint", 1: "NormalizedMiddleBreakdownPoint", 
        2: "LongerBreakdownPoint", 3: "NormalizedLongerBreakdownPoint"
    }
)

flat_chimeric_result_df = pd.concat(
    [
        flat_chimeric_result_df,
        res
    ],
    axis=1
)

flat_chimeric_result_df

# %%
flat_chimeric_result_df.sort_values(["Platform", "Sample"])

# %%
flat_chimeric_result_df.loc[:, ["ShortPlatform", "Sample"]].drop_duplicates()

# %%
flat_chimeric_result_df["ShortPlatform"].drop_duplicates().sort_values()

# %%
list(sorted(set(platforms)))

# %%
frac = 0.3
# frac = 0.5

fig = make_subplots(
    rows=len(X_common_proteins), 
    cols=len(Y_rare_proteins), 
    shared_xaxes="all", 
    shared_yaxes="all",
    y_title="Normalized middle breakdown point [%]",
    column_titles=col_titles,
    row_titles=row_titles,
    vertical_spacing=0.2/len(X_common_proteins),
    horizontal_spacing=0.1/len(Y_rare_proteins),
)

for i, x_common_proteins in enumerate(X_common_proteins, start=1):
    for j, y_rare_proteins in enumerate(Y_rare_proteins, start=1):
        for platform in list(sorted(set(platforms))):
            subset_df = flat_chimeric_result_df[
                (flat_chimeric_result_df["XCommonProteins"] == x_common_proteins) &
                (flat_chimeric_result_df["YRareProteins"] == y_rare_proteins) &
                (flat_chimeric_result_df["Platform"] == platform)
            ]
            subset_df = subset_df.sample(frac=frac, random_state=1892)  # sample 10% of the data for plotting to prevent crashes due to too many points
            subset_df = subset_df.sort_values(["Platform", "Sample"])
            box_color = platform_colormap[platform]  # first platform in that subset
            fig.add_trace(
                go.Box(
                    x=subset_df.loc[:, ["ShortPlatform", "Sample"]].T.values.tolist(),
                    y=subset_df["NormalizedMiddleBreakdownPoint"],
                    marker=dict(
                        color=box_color
                    ),
                ),
                row=i,
                col=j
            )
                
fig.update_xaxes(
    tickfont_size=10
)
fig.update_yaxes(
    dtick=20
)

fig.update_layout(
    template="plotly_white",
    showlegend=False,
    width=1400,
    height=800,
    title=f"{np.round(100*frac, 0)}% of all chimerizing pairs"
)

# fig.show()
fig.show(config={'staticPlot': True})

# %%
frac = 0.3
# frac = 0.5

fig = make_subplots(
    rows=len(X_common_proteins), 
    cols=len(Y_rare_proteins), 
    shared_xaxes="all", 
    shared_yaxes="all",
    y_title="Normalized longer breakdown point [%]",
    column_titles=col_titles,
    row_titles=row_titles,
    vertical_spacing=0.2/len(X_common_proteins),
    horizontal_spacing=0.1/len(Y_rare_proteins),
)

for i, x_common_proteins in enumerate(X_common_proteins, start=1):
    for j, y_rare_proteins in enumerate(Y_rare_proteins, start=1):
        for platform in list(sorted(set(platforms))):
            subset_df = flat_chimeric_result_df[
                (flat_chimeric_result_df["XCommonProteins"] == x_common_proteins) &
                (flat_chimeric_result_df["YRareProteins"] == y_rare_proteins) &
                (flat_chimeric_result_df["Platform"] == platform)
            ]
            subset_df = subset_df.sample(frac=frac, random_state=1892)  # sample 10% of the data for plotting to prevent crashes due to too many points
            subset_df = subset_df.sort_values(["Platform", "Sample"])
            box_color = platform_colormap[platform]  # first platform in that subset
            fig.add_trace(
                go.Box(
                    x=subset_df.loc[:, ["ShortPlatform", "Sample"]].T.values.tolist(),
                    y=subset_df["NormalizedLongerBreakdownPoint"],
                    marker=dict(
                        color=box_color
                    ),
                ),
                row=i,
                col=j
            )
                
fig.update_xaxes(
    tickfont_size=10
)
fig.update_yaxes(
    dtick=20
)

fig.update_layout(
    template="plotly_white",
    showlegend=False,
    width=1400,
    height=800,
    title=f"{np.round(100*frac, 0)}% of all chimerizing pairs"
)

# fig.show()
fig.show(config={'staticPlot': True})

# %%
flat_chimeric_result_df

# %%
flat_chimeric_result_df.columns

# %%
mean_chimeric_result_df = (
    flat_chimeric_result_df
    .groupby(
        [
            'Platform', 'ShortPlatform', 'Sample', 'TotalNumOfAAs',
            'XCommonProteins', 'YRareProteins', 'Protein',
        ]
    )
    .agg(
        NumOfChimericCombinations = ("ChimericUniqueProteinsCombinations", "count"),
        MeanNormalizedMiddleBreakdownPoint = ("NormalizedMiddleBreakdownPoint", "mean"),
        # STDNormalizedMiddleBreakdownPoint = ("NormalizedMiddleBreakdownPoint", "std"),
        MeanNormalizedLongerBreakdownPoint = ("NormalizedLongerBreakdownPoint", "mean"),
        # STDNormalizedLongerBreakdownPoint = ("NormalizedLongerBreakdownPoint", "std"),
    )
    .reset_index()
)
mean_chimeric_result_df

# %%
chimeric_result_df.shape[0]

# %%
frac = 0.3
# frac = 0.5

fig = make_subplots(
    rows=len(X_common_proteins), 
    cols=len(Y_rare_proteins), 
    shared_xaxes="all", 
    shared_yaxes="all",
    y_title="Per-protein mean normalized middle breakdown point [%]",
    column_titles=col_titles,
    row_titles=row_titles,
    vertical_spacing=0.2/len(X_common_proteins),
    horizontal_spacing=0.1/len(Y_rare_proteins),
)

for i, x_common_proteins in enumerate(X_common_proteins, start=1):
    for j, y_rare_proteins in enumerate(Y_rare_proteins, start=1):
        for platform in list(sorted(set(platforms))):
            subset_df = mean_chimeric_result_df[
                (mean_chimeric_result_df["XCommonProteins"] == x_common_proteins) &
                (mean_chimeric_result_df["YRareProteins"] == y_rare_proteins) &
                (mean_chimeric_result_df["Platform"] == platform)
            ]
            # subset_df = subset_df.sample(frac=frac, random_state=1892)  # sample 10% of the data for plotting to prevent crashes due to too many points
            subset_df = subset_df.sort_values(["Platform", "Sample"])
            box_color = platform_colormap[platform]  # first platform in that subset
            fig.add_trace(
                go.Box(
                    x=subset_df.loc[:, ["ShortPlatform", "Sample"]].T.values.tolist(),
                    y=subset_df["MeanNormalizedMiddleBreakdownPoint"],
                    marker=dict(
                        color=box_color
                    ),
                ),
                row=i,
                col=j
            )
                
fig.update_xaxes(
    tickfont_size=10
)
fig.update_yaxes(
    dtick=20
)

fig.update_layout(
    template="plotly_white",
    showlegend=False,
    width=1400,
    height=800,
    # title=f"{np.round(100*frac, 0)}% of all chimerizing pairs"
)

# fig.show()
fig.show(config={'staticPlot': True})

# %%
frac = 0.3
# frac = 0.5

fig = make_subplots(
    rows=len(X_common_proteins), 
    cols=len(Y_rare_proteins), 
    shared_xaxes="all", 
    shared_yaxes="all",
    y_title="Per-protein mean normalized longer breakdown point [%]",
    column_titles=col_titles,
    row_titles=row_titles,
    vertical_spacing=0.2/len(X_common_proteins),
    horizontal_spacing=0.1/len(Y_rare_proteins),
)

for i, x_common_proteins in enumerate(X_common_proteins, start=1):
    for j, y_rare_proteins in enumerate(Y_rare_proteins, start=1):
        for platform in list(sorted(set(platforms))):
            subset_df = mean_chimeric_result_df[
                (mean_chimeric_result_df["XCommonProteins"] == x_common_proteins) &
                (mean_chimeric_result_df["YRareProteins"] == y_rare_proteins) &
                (mean_chimeric_result_df["Platform"] == platform)
            ]
            # subset_df = subset_df.sample(frac=frac, random_state=1892)  # sample 10% of the data for plotting to prevent crashes due to too many points
            subset_df = subset_df.sort_values(["Platform", "Sample"])
            box_color = platform_colormap[platform]  # first platform in that subset
            fig.add_trace(
                go.Box(
                    x=subset_df.loc[:, ["ShortPlatform", "Sample"]].T.values.tolist(),
                    y=subset_df["MeanNormalizedLongerBreakdownPoint"],
                    marker=dict(
                        color=box_color
                    ),
                ),
                row=i,
                col=j
            )
                
fig.update_xaxes(
    tickfont_size=10
)
fig.update_yaxes(
    dtick=20
)

fig.update_layout(
    template="plotly_white",
    showlegend=False,
    width=1400,
    height=800,
    # title=f"{np.round(100*frac, 0)}% of all chimerizing pairs"
)

# fig.show()
fig.show(config={'staticPlot': True})
