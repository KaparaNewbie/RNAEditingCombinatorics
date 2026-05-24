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

# %% [markdown]
# # Prep

# %%
from pathlib import Path
from itertools import product, chain

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

# %%
pacbio_expression_files = [
	# PacBio1
	"/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/GRIA.DistinctUniqueProteins.ExpressionLevels.EntropyConsidered.csv",
	"/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO.DistinctUniqueProteins.ExpressionLevels.EntropyConsidered.csv",
	# PacBio2
	"/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples/ADAR1.DistinctUniqueProteins.ExpressionLevels.EntropyConsidered.csv",
	"/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples/IQEC1.DistinctUniqueProteins.ExpressionLevels.EntropyConsidered.csv",
	# PacBio3
	"/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads/GRIA2.DistinctUniqueProteins.ExpressionLevels.EntropyConsidered.csv",
	"/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads/ADAR1.DistinctUniqueProteins.ExpressionLevels.EntropyConsidered.csv",
	"/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads/IQEC1.DistinctUniqueProteins.ExpressionLevels.EntropyConsidered.csv",
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
unique_reads_first_col_pos = 8

pacbio_unique_reads_files = [
	# PacBio1
	"/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_reads.csv.gz",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_reads.csv.gz",
	# PacBio2
	"/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples/ADAR1.Merged.r64296e203404D01.aligned.sorted.MinRQ998.unique_reads.csv.gz",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples/IQEC.Merged.r64296e203404D01.aligned.sorted.MinRQ998.unique_reads.csv.gz",
	# PacBio3
	'/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads/comp141693_c0_seq1.merged.MinRQ998.unique_reads.csv.gz',
    '/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads/comp134400_c0_seq1_extended.merged.MinRQ998.unique_reads.csv.gz',
    '/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads/comp141565_c6_seq3.merged.MinRQ998.unique_reads.csv.gz'
]

illumina_unique_reads_files = [
	f"/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.{chrom}.unique_reads.csv"
    for chrom in illumina_chroms
]

# %%
platforms = pacbio_platforms + illumina_platforms
samples = pacbio_samples + illumina_samples
# unique_proteins_files = pacbio_unique_proteins_files + illumina_unique_proteins_files
unique_reads_files = pacbio_unique_reads_files + illumina_unique_reads_files
expression_files = pacbio_expression_files + illumina_expression_files

expression_files = [Path(x) for x in expression_files]

# %%
expression_files


# %%
def find_max_solution(expression_file, interfix="MaxSolution_"):
    max_sol_files = list(
        expression_file.parent.glob(
            f"{expression_file.name}.{interfix}*"
        )
    )
    if len(max_sol_files) == 0:
        raise FileNotFoundError(f"No max solution file found for {expression_file}")
    if len(max_sol_files) > 1:
        raise ValueError(f"Multiple max solution files for {expression_file}: {max_sol_files}")
    max_sol_file = max_sol_files[0]
    max_sol = int(max_sol_file.name.rsplit("_")[-1])
    return max_sol


# %%
def make_expression_df(expression_file):
    assignment_df = pd.read_csv(
        expression_file, 
        sep="\t", 
        dtype={
            "#Solution": int,
        },
        usecols = [
            'Gene', 'Protein', '#Solution', 'SufficientEntropy', 'TotalWeightedSupportingReads'
        ]
    )
    
    max_sol = find_max_solution(expression_file)
    assignment_df = assignment_df.loc[
        assignment_df["#Solution"].eq(max_sol)
    ]

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
X_common_proteins = [0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05]
Y_rare_proteins = [0.1, 0.2, 0.3]

x_and_y_proteins_denote_fractions = True

# %%
num_of_x_y_combinations = len(list(product(X_common_proteins, Y_rare_proteins)))
num_of_x_y_combinations

# %%
soft_comparisons = [True, False]

# %%
out_dir = Path("/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/JointChimericReadsAnalysis2")


# %%
# !ls {out_dir} | head -n 3

# %%
def define_out_file(out_dir, platform, sample_name, X_common_proteins, Y_rare_proteins):
    out_file = Path(
		out_dir,
		f"{platform}.{sample_name}.X_{'_'.join(map(str, X_common_proteins))}.Y_{'_'.join(map(str, Y_rare_proteins))}.csv.gz"
	)
    if not out_file.exists():
        raise FileNotFoundError(f"Expected output file {out_file} does not exist.")
    return out_file

out_files = [
	define_out_file(out_dir, platform, sample_name, X_common_proteins, Y_rare_proteins)
	for platform, sample_name in zip(
        platforms, samples
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
num_of_editing_sites_per_platform_and_sample = {
    (platform, sample): (
        pd.read_table(unique_reads_file, nrows=5)
        .iloc[:, unique_reads_first_col_pos:]
        .shape[1]
    )
    for platform, sample, unique_reads_file in zip(platforms, samples, unique_reads_files)
}

num_of_editing_sites_per_platform_and_sample_df = (
    pd.DataFrame(num_of_editing_sites_per_platform_and_sample, index=["Platform", "Sample", "TotalNumOfEditingSites"])
    .T
    .reset_index()
    .drop(columns=["Platform", "Sample"])
    .rename(columns={"level_0": "Platform", "level_1": "Sample"})
)
num_of_editing_sites_per_platform_and_sample_df


# %%
def read_results_file(out_file):
    df = pd.read_csv(
        out_file, 
        sep="\t",
        dtype={
            "Platform": str,
            "Sample": str,
            "Solution": int,
            "IsSoftComparison": bool,
            "XCommonProteins": float,
            "YRareProteins": float,
            "XYProteinsDenoteFractions": bool,
            "ActualXCommonProteins": int,
            "ActualYRareProteins": int,
            "Protein": str,
            "UniqueRead": str,
            "IsChimeric": bool,
            "NumOfChimericCombinations": int,
            "ChimerizingProteinPairs": str,
            "ChimerizingReadPairs": str,
            "ChimerizingReadPairsIntersectingSitesIndices": str
        }
    )
    df["ChimerizingProteinPairs"] = df["ChimerizingProteinPairs"].apply(parse_pairs_str)
    df["ChimerizingReadPairs"] = df["ChimerizingReadPairs"].apply(parse_pairs_str)
    df["ChimerizingReadPairsIntersectingSitesIndices"] = df["ChimerizingReadPairsIntersectingSitesIndices"].apply(parse_int_pairs_str)
    return df


# %%
result_df = pd.concat(
    [read_results_file(out_file) for out_file in out_files],
    ignore_index=True
)
result_df = num_of_editing_sites_per_platform_and_sample_df.merge(
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
results_categories_df = (
    result_df
    .loc[
        :,
        [
            "Platform", "Sample", "XCommonProteins", "YRareProteins", "IsSoftComparison"
        ]
    ]
    .drop_duplicates()
)
results_categories_df["ExistingCombination"] = True

for platform, sample in zip(
    platforms, samples
):
    for x_common_proteins, y_rare_proteins, soft_comparison in product(
        X_common_proteins, Y_rare_proteins, soft_comparisons
    ):
        if results_categories_df.loc[
            (results_categories_df["Platform"] == platform) &
            (results_categories_df["Sample"] == sample) &
            (results_categories_df["XCommonProteins"] == x_common_proteins) &
            (results_categories_df["YRareProteins"] == y_rare_proteins) &
            (results_categories_df["IsSoftComparison"] == soft_comparison)
        ].empty:
            results_categories_df = pd.concat([
                results_categories_df,
                pd.DataFrame({
                    "Platform": [platform],
                    "Sample": [sample],
                    "XCommonProteins": [x_common_proteins],
                    "YRareProteins": [y_rare_proteins],
                    "IsSoftComparison": [soft_comparison],
                    "ExistingCombination": [False]
                })
            ], ignore_index=True)

assert results_categories_df.shape[0] == len(platforms) * len(X_common_proteins) * len(Y_rare_proteins) * len(soft_comparisons), (
    f"Expected {len(platforms) * len(X_common_proteins) * len(Y_rare_proteins) * len(soft_comparisons)} rows in results_categories_df, "
    f"but got {results_categories_df.shape[0]} rows. Please check the logic for generating all combinations."
)

results_categories_df = results_categories_df.sort_values(
    by=["Platform", "Sample", "XCommonProteins", "YRareProteins", "IsSoftComparison"],
    ignore_index=True
)

results_categories_df

# %%
(
    results_categories_df
    .groupby(["XCommonProteins", "YRareProteins", "IsSoftComparison"])
    ["ExistingCombination"].value_counts()
    .reset_index(name="Count")
)

# %% [markdown]
# The `result_df` has 1 row per rare unique read, and since a rare protein can be composed of multiple such rare   
# unique reads (although very rarely), and each such unique read can be chimeric in multiple ways,  
# we need to agg the `result_df` to have 1 row per rare protein.

# %%
proteins_results_df = (
    result_df
    .assign(
        UniqueReadPerChimerization=lambda df: [
            [read] * int(n)
            for read, n in zip(df["UniqueRead"], df["NumOfChimericCombinations"])
        ]
    )
    .groupby(
        [
            'Platform', 'ShortPlatform', 'Sample', 'TotalNumOfEditingSites',
            'Solution', 'IsSoftComparison', 'XCommonProteins', 'YRareProteins',
            'XYProteinsDenoteFractions', 'ActualXCommonProteins',
            'ActualYRareProteins', 'Protein'
        ],
        # as_index=False
    )
    .agg(
        IsChimeric = ("IsChimeric", "any"),
        NumOfUniqueReads = ("UniqueRead", "size"),
        UniqueReads = ("UniqueRead", list),
        IsUniqueReadChimeric = ("IsChimeric", list),
        NumOfChimericCombinations = ("NumOfChimericCombinations", "sum"),
        UniqueReadPerChimerization = ("UniqueReadPerChimerization", lambda x: list(chain.from_iterable(x))),
        ChimerizingProteinPairs = ("ChimerizingProteinPairs", lambda x: list(chain.from_iterable(x))),
        ChimerizingReadPairs = ("ChimerizingReadPairs", lambda x: list(chain.from_iterable(x))),
        ChimerizingReadPairsIntersectingSitesIndices = ("ChimerizingReadPairsIntersectingSitesIndices", lambda x: list(chain.from_iterable(x)))
    )
    .reset_index()
)
proteins_results_df.insert(
    proteins_results_df.columns.get_loc("IsUniqueReadChimeric") + 1,
    "%OfUniqueReadsChimeric",
    proteins_results_df["IsUniqueReadChimeric"].apply(
        lambda x: 100 * sum(x) / len(x)
    )
)
assert (
    proteins_results_df.loc[
        :,
        ["UniqueReadPerChimerization", "ChimerizingProteinPairs", "ChimerizingReadPairs", "ChimerizingReadPairsIntersectingSitesIndices"]
    ]
    .map(len).apply(set, axis=1)
    .apply(len).eq(1).all()
), "Expected all lists in the specified columns to have the same length within each row, but found a row with differing lengths."
proteins_results_df

# %%
proteins_results_df.loc[
    proteins_results_df["NumOfUniqueReads"].gt(1)
].shape[0]

# %%
proteins_results_df.loc[
    proteins_results_df["NumOfUniqueReads"].gt(1),
    ["Platform", "Sample"]
].value_counts()

# %%
proteins_results_df.loc[
    proteins_results_df["NumOfUniqueReads"].gt(1),
    "%OfUniqueReadsChimeric"
].value_counts()

# %%
# def make_x_y_expression_df(
#     platforms,
#     samples,
#     expression_files,
#     X_common_proteins,
#     Y_rare_proteins,
#     stats_df
# ):

#     results = []

#     for platform, sample, expression_file in zip(
#         platforms, samples, expression_files
#     ):
#         expression_df = make_expression_df(expression_file)
        
#         # ic(expression_df.dtypes)
        
#         for x_common_proteins, y_rare_proteins in product(X_common_proteins, Y_rare_proteins):
            
#             # we drop duplicates here because we get 2 identical rows due to the soft comparison 
#             # (one with soft_comparison=True and one with soft_comparison=False); 
#             # however, from the X and Y proteins perspective they are the same, 
#             # so we only need one of the rows to get the actual X and Y proteins
            
#             # actual_x_common_proteins, actual_y_rare_proteins = stats_df.loc[
#             #     (stats_df["Platform"].eq(platform))
#             #     & (stats_df["Sample"].eq(sample))
#             #     & (stats_df["XCommonProteins"].eq(x_common_proteins))
#             #     & (stats_df["YRareProteins"].eq(y_rare_proteins)),
#             #     ["ActualXCommonProteins", "ActualYRareProteins"]
#             # ].drop_duplicates().squeeze()
            
#             actual_x_y_series = stats_df.loc[
#                 (stats_df["Platform"].eq(platform))
#                 & (stats_df["Sample"].eq(sample))
#                 & (stats_df["XCommonProteins"].eq(x_common_proteins))
#                 & (stats_df["YRareProteins"].eq(y_rare_proteins)),
#                 ["ActualXCommonProteins", "ActualYRareProteins"]
#             ].drop_duplicates().squeeze()
            
#             if not actual_x_y_series.empty:
            
#                 actual_x_common_proteins, actual_y_rare_proteins = actual_x_y_series
                
#                 if type(actual_x_common_proteins) is str or type(actual_y_rare_proteins) is str:
#                     ic(platform, sample, x_common_proteins, y_rare_proteins)
                
#                 x_common_proteins_cumulative_expression = expression_df.loc[
#                     expression_df["CommonToRareProteinIndex"].le(actual_x_common_proteins),
#                     "%RelativeExpression"
#                 ].sum()
#                 y_rare_proteins_cumulative_expression = expression_df.loc[
#                     expression_df["RareToCommonProteinIndex"].le(actual_y_rare_proteins),
#                     "%RelativeExpression"
#                 ].sum()
                
#             else:
#                 actual_x_common_proteins = np.nan
#                 actual_y_rare_proteins = np.nan
#                 x_common_proteins_cumulative_expression = np.nan
#                 y_rare_proteins_cumulative_expression = np.nan
            
#             # ic(type(actual_x_common_proteins), type(actual_y_rare_proteins))
            
#             for soft_comparison in [True, False]:
#                 result = [
#                     platform, sample, x_common_proteins, y_rare_proteins, soft_comparison,
#                     x_common_proteins_cumulative_expression, y_rare_proteins_cumulative_expression
#                 ]
        
#                 results.append(result)
            
            
#     x_y_expression_df = pd.DataFrame(
#         results,
#         columns=[
#             "Platform", "Sample", "XCommonProteins", "YRareProteins", "IsSoftComparison",
#             "XCommonProteinsCumulativeExpression", "YRareProteinsCumulativeExpression"
#         ]
#     )
    
#     assert len(platforms) * len(X_common_proteins) * len(Y_rare_proteins) * 2 == x_y_expression_df.shape[0], (
#         f"Expected {len(platforms) * len(X_common_proteins) * len(Y_rare_proteins) * 2} rows in the resulting "
#         "DataFrame, but got {x_y_expression_df.shape[0]} rows. "
#         "Please check the input lists and the logic of the function."
#     )
    
#     return x_y_expression_df

# %%
def make_x_y_expression_df(
    platforms,
    samples,
    expression_files,
    X_common_proteins,
    Y_rare_proteins,
    # stats_df
):
    results = []

    for platform, sample, expression_file in zip(
        platforms, samples, expression_files
    ):
        expression_df = make_expression_df(expression_file)
        
        num_of_distinct_prots = expression_df.shape[0]
        
        # ic(expression_df.dtypes)
        
        for x_common_proteins, y_rare_proteins in product(X_common_proteins, Y_rare_proteins):
            
            actual_x_common_proteins = int(np.round(x_common_proteins * num_of_distinct_prots))
            actual_y_rare_proteins = int(np.round(y_rare_proteins * num_of_distinct_prots))

            x_common_proteins_cumulative_expression = expression_df.loc[
                expression_df["CommonToRareProteinIndex"].le(actual_x_common_proteins),
                "%RelativeExpression"
            ].sum()
            y_rare_proteins_cumulative_expression = expression_df.loc[
                expression_df["RareToCommonProteinIndex"].le(actual_y_rare_proteins),
                "%RelativeExpression"
            ].sum()
        
        
            result = [
                platform, sample, 
                x_common_proteins, y_rare_proteins,
                actual_x_common_proteins, actual_y_rare_proteins,
                x_common_proteins_cumulative_expression, y_rare_proteins_cumulative_expression
            ]
        
            results.append(result)
            
            
    x_y_expression_df = pd.DataFrame(
        results,
        columns=[
            "Platform", "Sample", 
            "XCommonProteins", "YRareProteins",
            "ActualXCommonProteins", "ActualYRareProteins",
            "XCommonProteinsCumulativeExpression", "YRareProteinsCumulativeExpression"
        ]
    )
    
    assert len(platforms) * len(X_common_proteins) * len(Y_rare_proteins) == x_y_expression_df.shape[0], (
        f"Expected {len(platforms) * len(X_common_proteins) * len(Y_rare_proteins)} rows in the resulting "
        f"DataFrame, but got {x_y_expression_df.shape[0]} rows. "
        "Please check the input lists and the logic of the function."
    )
    
    return x_y_expression_df

# %%
x_y_expression_df = make_x_y_expression_df(
    platforms,
    samples,
    expression_files,
    X_common_proteins,
    Y_rare_proteins,
    # proteins_stats_df
)
x_y_expression_df

# %%
x_y_expression_df.isna().sum()

# %%
proteins_stats_df = (
    proteins_results_df
    .groupby(
        [
            "Platform", "ShortPlatform", "Sample", "XYProteinsDenoteFractions", "XCommonProteins", "YRareProteins", "IsSoftComparison",
            # 'ActualXCommonProteins', 'ActualYRareProteins'
         ],
        as_index=False
    )
    .agg(
        NumOfChimericProteinsOnEditingSitesLevel = ("IsChimeric", "sum"),
    )
)

proteins_stats_df.insert(
    proteins_stats_df.columns.get_loc("Sample") + 1,
    "Platform-Sample",
    proteins_stats_df["ShortPlatform"] + "-" + proteins_stats_df["Sample"]
)

added_empty_rows = []
for platform, sample in zip(platforms, samples):
    for x_common_proteins, y_rare_proteins in product(X_common_proteins, Y_rare_proteins):
        for soft_comparison in soft_comparisons:
            subset_df = proteins_stats_df[
                (proteins_stats_df["Platform"].eq(platform))
                & (proteins_stats_df["Sample"].eq(sample))
                & (proteins_stats_df["XCommonProteins"].eq(x_common_proteins))
                & (proteins_stats_df["YRareProteins"].eq(y_rare_proteins))
                & (proteins_stats_df["IsSoftComparison"].eq(soft_comparison))
            ]
            if subset_df.empty:
                added_empty_rows.append(
                    {
                        "Platform": platform,
                        "ShortPlatform": shortened_platform_names[platform],
                        "Sample": sample,
                        "Platform-Sample": shortened_platform_names[platform] + "-" + sample,
                        "XYProteinsDenoteFractions": x_and_y_proteins_denote_fractions,
                        "XCommonProteins": x_common_proteins,
                        "YRareProteins": y_rare_proteins,
                        "IsSoftComparison": soft_comparison,
                        # "ActualXCommonProteins": np.nan,
                        # "ActualYRareProteins": np.nan,
                        "NumOfChimericProteinsOnEditingSitesLevel": 0,
                        # "%OfChimericProteinsOnEditingSitesLevel": 0.0
                    }
                )
proteins_stats_df = pd.concat([proteins_stats_df, pd.DataFrame(added_empty_rows)], ignore_index=True)

proteins_stats_df = proteins_stats_df.merge(
    x_y_expression_df,
    how="left"
)

proteins_stats_df.insert(
    proteins_stats_df.columns.get_loc("NumOfChimericProteinsOnEditingSitesLevel") + 1,
    "%OfChimericProteinsOnEditingSitesLevel",
    proteins_stats_df["NumOfChimericProteinsOnEditingSitesLevel"].mul(100).div(proteins_stats_df["ActualYRareProteins"])
)

proteins_stats_df

# %%
assert proteins_stats_df.isna().sum().eq(0).all()

# %% [markdown]
# # % of chimeric proteins

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

# %%
# row_titles = [
#     f"Common: {x_common_proteins * 100:.0f}%"
#     for x_common_proteins in X_common_proteins
# ]
row_titles = [
    f"Common: {x_common_proteins * 100:.1f}%"
    for x_common_proteins in X_common_proteins
]
col_titles = [
    f"Rare: {y_rare_proteins * 100:.0f}%"
    for y_rare_proteins in Y_rare_proteins
]
row_titles, col_titles

# %%
for soft_comparison in [False, True]:

    fig = make_subplots(
        # rows=len(X_common_proteins), 
        # cols=len(Y_rare_proteins), 
        rows=len(Y_rare_proteins), 
        cols=len(X_common_proteins), 
        shared_xaxes="all", 
        shared_yaxes="all",
        y_title="% of chimeric proteins on editing sites level",
        # subplot_titles=subplot_titles,
        # column_titles=col_titles,
        # row_titles=row_titles,
        column_titles=row_titles,
        row_titles=col_titles,
        # vertical_spacing=0.2/len(X_common_proteins),
        # horizontal_spacing=0.1/len(Y_rare_proteins),
        vertical_spacing=0.1/len(Y_rare_proteins),
        horizontal_spacing=0.2/len(X_common_proteins),
    )
    
    max_y = 0

    # for each row i
    # for i, x_common_proteins in enumerate(X_common_proteins, start=1):
    for i, y_rare_proteins in enumerate(Y_rare_proteins, start=1):
        # for each column j
        # for j, y_rare_proteins in enumerate(Y_rare_proteins, start=1):
        for j, x_common_proteins in enumerate(X_common_proteins, start=1):
            subset_df = proteins_stats_df[
                (proteins_stats_df["XCommonProteins"] == x_common_proteins) &
                (proteins_stats_df["YRareProteins"] == y_rare_proteins) &
                (proteins_stats_df["IsSoftComparison"] == soft_comparison)
            ]
            if subset_df.shape[0] < 26:
                ic(x_common_proteins, y_rare_proteins, subset_df.shape[0])
            x = subset_df.loc[:, ["ShortPlatform", "Sample"]].T.values.tolist()
            y = subset_df["%OfChimericProteinsOnEditingSitesLevel"]
            max_y = max(max_y, y.max())
            fig.add_trace(
                go.Bar(
                    x=x,
                    y=y,
                    marker=dict(
                        color=subset_df["Platform"].map(platform_colormap)
                    ),
                ),
                row=i,
                col=j
            )
    fig.update_xaxes(
        # tickfont_size=10
        tickfont_size=5
    )
    fig.update_yaxes(
        # dtick=20
        # dtick=20 if max_y >= 80 else 10
        dtick=10
    )

    fig.update_layout(
        template="plotly_white",
        showlegend=False,
        width=1700,
        height=800,
        title=f"Soft Comparison = {soft_comparison}"
    )
            
    fig.show()

# %%
metric_colors = {
    "XCommonProteinsCumulativeExpression": "#1f77b4",
    "YRareProteinsCumulativeExpression":   "#ff7f0e",
}


fig = make_subplots(
    # rows=len(X_common_proteins),
    # cols=len(Y_rare_proteins),
    rows=len(Y_rare_proteins),
    cols=len(X_common_proteins),
    shared_xaxes="all",
    shared_yaxes="all",
    y_title="Cumulative expression [%]",
    # column_titles=col_titles,
    # row_titles=row_titles,
    # vertical_spacing=0.2 / len(X_common_proteins),
    # horizontal_spacing=0.1 / len(Y_rare_proteins),
    column_titles=row_titles,
    row_titles=col_titles,
    vertical_spacing=0.1 / len(Y_rare_proteins),
    # horizontal_spacing=0.2 / len(X_common_proteins),
    horizontal_spacing=0.1 / len(X_common_proteins),
)

# for i, x_common_proteins in enumerate(X_common_proteins, start=1):
#     for j, y_rare_proteins in enumerate(Y_rare_proteins, start=1):
for i, y_rare_proteins in enumerate(Y_rare_proteins, start=1):
    for j, x_common_proteins in enumerate(X_common_proteins, start=1):
        
        subset_df = proteins_stats_df[
            (proteins_stats_df["XCommonProteins"] == x_common_proteins)
            & (proteins_stats_df["YRareProteins"] == y_rare_proteins)
            # & (proteins_stats_df["IsSoftComparison"].eq(soft_comparison))
        ].copy()
        
        subset_df = subset_df.drop_duplicates(
            subset=["Platform", "Sample", "ShortPlatform", "XCommonProteins", "YRareProteins"]
        )

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
fig.update_xaxes(tickfont_size=6)
fig.update_yaxes(
    # dtick=20,
    dtick=10,
    range=[0, 100]
#     type="log"
)
fig.update_layout(
    template="plotly_white",
    width=1700,
    height=900,
    barmode="group",
    legend=dict(orientation="h", yanchor="top", y=-0.2, xanchor="left", x=0),
    # title_text=f"Soft Comparison = {soft_comparison}"
)
fig.show()

# %% [markdown]
# # Intersection points

# %%
result_df

# %%
# proteins_results_df

# %%
chimeric_result_df = (
    result_df.loc[
        result_df["IsChimeric"],
        [
            "Platform", "ShortPlatform", "Sample", "TotalNumOfEditingSites", "IsSoftComparison", 
            "XCommonProteins", "YRareProteins", "Protein", "UniqueRead",
            "ChimerizingProteinPairs", "ChimerizingReadPairs", "ChimerizingReadPairsIntersectingSitesIndices"
        ]
    ]
    .reset_index(drop=True)
    # .drop(columns="IsChimericOnAALevel")
)

chimeric_result_df


# %%
def compute_breakdown_points(b, e, total_editing_sites):
    b_len = b # 1 is the first editing site, as this is a julia result, so length is b, not b-1
    e_len = total_editing_sites - e + 1

    # this breakdown point represent the midpoint of the shared region
    middle_breakdown_point = e + ((b - e) / 2)

    # this breakdown point represents the end of the longer region (or the end of the shared region if they are of equal length)
    if b_len >= e_len:
        longer_breakdown_point = b
    else:
        longer_breakdown_point = e
        
    normalized_middle_breakdown_point = 100 * middle_breakdown_point / total_editing_sites
    normalized_longer_breakdown_point = 100 * longer_breakdown_point / total_editing_sites

    return middle_breakdown_point, normalized_middle_breakdown_point, longer_breakdown_point, normalized_longer_breakdown_point


# %%
flat_chimeric_result_df = (
    chimeric_result_df
    # .drop(columns="NumOfChimericCombinations")
    .explode(
        [
            "ChimerizingProteinPairs", "ChimerizingReadPairs", "ChimerizingReadPairsIntersectingSitesIndices"
        ],
        ignore_index=True
    )
)

res = flat_chimeric_result_df.loc[
    :, ["ChimerizingReadPairsIntersectingSitesIndices", "TotalNumOfEditingSites"]
].apply(
    lambda row: compute_breakdown_points(
        *row["ChimerizingReadPairsIntersectingSitesIndices"], 
        row["TotalNumOfEditingSites"]
    ),
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
# # previoiusly we had 10,000,000 unique chimerizing pairs and we took 0.3 of it - 3,000,000;
# # to take now 3,000,000 out of 79,762,758 we
# # frac * (79 * 10**6) --> 3000000
# frac = 3000000 / (79 * 10**6)
# frac

# %%
# frac = 0.05
frac = 0.01

# %%
sampled_flat_chimeric_result_df = (
    flat_chimeric_result_df
    .groupby(
        ["Platform", "ShortPlatform", "Sample", "IsSoftComparison", "XCommonProteins", "YRareProteins"]
    )
    .sample(
        frac=frac, 
        random_state=np.random.default_rng(seed=1892)
    )
    .sort_values(["ShortPlatform", "Sample"])  
)
sampled_flat_chimeric_result_df

# %%
for soft_comparison in [False, True]:
    subset_df = sampled_flat_chimeric_result_df[
        sampled_flat_chimeric_result_df["IsSoftComparison"] == soft_comparison
    ]

    fig = px.scatter(
        subset_df,
        x="NormalizedMiddleBreakdownPoint",
        y="NormalizedLongerBreakdownPoint",
        # histnorm='probability',
        facet_col="XCommonProteins",
        facet_row="YRareProteins",
        category_orders={
            "XCommonProteins": X_common_proteins,
            "YRareProteins": Y_rare_proteins
        },
        labels={
            "NormalizedMiddleBreakdownPoint": "Middle<br>point",
            "NormalizedLongerBreakdownPoint": "Longer<br>point",
        },
        title=f"Breakdown points of chimeric reads (soft comparison = {soft_comparison} | sampled {frac:.1%} of chimeric read pairs)"
    )
    fig.update_layout(
        width=1600,
        height=700,
        showlegend=False,
        template="plotly_white",
    )
    # fig.show()
    fig.show(config={'staticPlot': True})

# %%
for soft_comparison in [False, True]:
    subset_df = sampled_flat_chimeric_result_df[
        sampled_flat_chimeric_result_df["IsSoftComparison"] == soft_comparison
    ].copy()
    subset_df["XCommonProteins"] = subset_df["XCommonProteins"].astype("category")

    fig = px.scatter(
        subset_df,
        x="NormalizedMiddleBreakdownPoint",
        y="NormalizedLongerBreakdownPoint",
        color="XCommonProteins",
        color_discrete_map={
            x: color
            for x, color in zip(
                X_common_proteins,
                px.colors.qualitative.Plotly[:len(X_common_proteins)]
            )
        },
        # symbol="YRareProteins",
        opacity=0.5,
        # facet_col="XCommonProteins",
        # facet_row="YRareProteins",
        # facet_col="YRareProteins",
        category_orders={
            # "XCommonProteins": X_common_proteins,
            "XCommonProteins": [str(x) for x in X_common_proteins],
            "YRareProteins": Y_rare_proteins
        },
        labels={
            "NormalizedMiddleBreakdownPoint": "Middle point",
            "NormalizedLongerBreakdownPoint": "Longer point",
        },
        title=f"Breakdown points of chimeric reads (soft comparison = {soft_comparison} | sampled {frac:.1%} of chimeric read pairs)"
    )
    fig.update_layout(
        width=1000,
        height=700,
        # showlegend=False,
        template="plotly_white",
    )
    # fig.show()
    fig.show(config={'staticPlot': True})

# %%
df = sampled_flat_chimeric_result_df

for soft_comparison in [False, True]:

    fig = make_subplots(
        rows=len(Y_rare_proteins),
        cols=len(X_common_proteins),
        shared_xaxes="all",
        shared_yaxes="all",
        y_title="Normalized middle breakdown point [%]",
        column_titles=row_titles,
        row_titles=col_titles,
        vertical_spacing=0.1 / len(Y_rare_proteins),
        # horizontal_spacing=0.2 / len(X_common_proteins), # too big
        horizontal_spacing=0.1 / len(X_common_proteins),  # fine
    )

    for i, y_rare_proteins in enumerate(Y_rare_proteins, start=1):
        for j, x_common_proteins in enumerate(X_common_proteins, start=1):
            for platform in list(sorted(set(platforms))):
                subset_df = df[
                    (df["XCommonProteins"] == x_common_proteins) &
                    (df["YRareProteins"] == y_rare_proteins) &
                    (df["Platform"] == platform) &
                    (df["IsSoftComparison"].eq(soft_comparison))
                ]
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
        # tickfont_size=10
        tickfont_size=6
    )
    fig.update_yaxes(
        dtick=20
    )

    fig.update_layout(
        template="plotly_white",
        showlegend=False,
        width=1700,
        height=800,
        title=f"Soft comparison = {soft_comparison} ({np.round(100*frac, 2)}% of all chimerizing pairs)"
    )

    # fig.show()
    fig.show(config={'staticPlot': True})

# %%
df = sampled_flat_chimeric_result_df

for soft_comparison in [False, True]:

    fig = make_subplots(
        rows=len(Y_rare_proteins),
        cols=len(X_common_proteins),
        shared_xaxes="all",
        shared_yaxes="all",
        y_title="Normalized longer breakdown point [%]",
        column_titles=row_titles,
        row_titles=col_titles,
        vertical_spacing=0.1 / len(Y_rare_proteins),
        # horizontal_spacing=0.2 / len(X_common_proteins), # too big
        horizontal_spacing=0.1 / len(X_common_proteins),  # fine
    )

    for i, y_rare_proteins in enumerate(Y_rare_proteins, start=1):
        for j, x_common_proteins in enumerate(X_common_proteins, start=1):
            for platform in list(sorted(set(platforms))):
                subset_df = df[
                    (df["XCommonProteins"] == x_common_proteins) &
                    (df["YRareProteins"] == y_rare_proteins) &
                    (df["Platform"] == platform) &
                    (df["IsSoftComparison"].eq(soft_comparison))
                ]
                box_color = platform_colormap[platform]
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
        # tickfont_size=10
        # tickfont_size=8
        tickfont_size=6
        # tickfont_size=4
    )
    fig.update_yaxes(
        dtick=20
    )

    fig.update_layout(
        template="plotly_white",
        showlegend=False,
        width=1700,
        height=800,
        title=f"Soft comparison = {soft_comparison} ({np.round(100*frac, 2)}% of all chimerizing pairs)"
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
            'Platform', 'ShortPlatform', 'Sample', 'TotalNumOfEditingSites', 'IsSoftComparison',
            'XCommonProteins', 'YRareProteins', 'Protein',
        ]
    )
    .agg(
        NumOfChimericCombinations = ("ChimerizingProteinPairs", "count"),
        MeanNormalizedMiddleBreakdownPoint = ("NormalizedMiddleBreakdownPoint", "mean"),
        MeanNormalizedLongerBreakdownPoint = ("NormalizedLongerBreakdownPoint", "mean"),
    )
    .reset_index()
    .sort_values(["ShortPlatform", "Sample"])
)
mean_chimeric_result_df

# %%
for soft_comparison in [False, True]:

    fig = make_subplots(
        rows=len(Y_rare_proteins),
        cols=len(X_common_proteins),
        shared_xaxes="all",
        shared_yaxes="all",
        y_title="Per-protein mean normalized middle breakdown point [%]",
        column_titles=row_titles,
        row_titles=col_titles,
        vertical_spacing=0.1 / len(Y_rare_proteins),
        # horizontal_spacing=0.2 / len(X_common_proteins), # too big
        horizontal_spacing=0.1 / len(X_common_proteins),  # fine
    )

    for i, y_rare_proteins in enumerate(Y_rare_proteins, start=1):
        for j, x_common_proteins in enumerate(X_common_proteins, start=1):
            for platform in list(sorted(set(platforms))):
                subset_df = mean_chimeric_result_df[
                    (mean_chimeric_result_df["XCommonProteins"] == x_common_proteins) &
                    (mean_chimeric_result_df["YRareProteins"] == y_rare_proteins) &
                    (mean_chimeric_result_df["Platform"] == platform) &
                    (mean_chimeric_result_df["IsSoftComparison"].eq(soft_comparison))
                ]
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
        # tickfont_size=10
        tickfont_size=6
    )
    fig.update_yaxes(
        dtick=20
    )

    fig.update_layout(
        template="plotly_white",
        showlegend=False,
        width=1700,
        height=800,
        # title=f"Soft comparison = {soft_comparison} ({np.round(100*frac, 2)}% of all chimerizing pairs)"
        title=f"Soft comparison = {soft_comparison}"
    )

    # fig.show()
    fig.show(config={'staticPlot': True})

# %%
for soft_comparison in [False, True]:

    fig = make_subplots(
        rows=len(Y_rare_proteins),
        cols=len(X_common_proteins),
        shared_xaxes="all",
        shared_yaxes="all",
        y_title="Per-protein mean normalized longer breakdown point [%]",
        column_titles=row_titles,
        row_titles=col_titles,
        vertical_spacing=0.1 / len(Y_rare_proteins),
        # horizontal_spacing=0.2 / len(X_common_proteins), # too big
        horizontal_spacing=0.1 / len(X_common_proteins),  # fine
    )

    for i, y_rare_proteins in enumerate(Y_rare_proteins, start=1):
        for j, x_common_proteins in enumerate(X_common_proteins, start=1):
            for platform in list(sorted(set(platforms))):
                subset_df = mean_chimeric_result_df[
                    (mean_chimeric_result_df["XCommonProteins"] == x_common_proteins) &
                    (mean_chimeric_result_df["YRareProteins"] == y_rare_proteins) &
                    (mean_chimeric_result_df["Platform"] == platform) &
                    (mean_chimeric_result_df["IsSoftComparison"].eq(soft_comparison))
                ]
                box_color = platform_colormap[platform]
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
        # tickfont_size=10
        # tickfont_size=8
        tickfont_size=6
        # tickfont_size=4
    )
    fig.update_yaxes(
        dtick=20
    )

    fig.update_layout(
        template="plotly_white",
        showlegend=False,
        width=1700,
        height=800,
        # title=f"Soft comparison = {soft_comparison} ({np.round(100*frac, 2)}% of all chimerizing pairs)"
        title=f"Soft comparison = {soft_comparison}"
        
    )

    # fig.show()
    fig.show(config={'staticPlot': True})
