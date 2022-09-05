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
# This notebook is an attempt to create and parse the mpileup file created by considering the complete ORF.

# %%
from pathlib import Path
from typing import Union
from itertools import chain
from collections import defaultdict

import pandas as pd
import numpy as np
from pybedtools import BedTool
from icecream import ic
from Bio import SeqIO, Seq


# %%
region = "comp141882_c0_seq14"
start = 0
end = 6294
strand = "+"
transcriptome = "/private7/projects/Combinatorics/D.pealeii/Annotations/orfs_squ.fa"

group_col = "Gene"
group = "PCLO"

assurance_factor = 1.5


# %%
assert (end - start) % 3 == 0


# %%
# cd_cmd = "cd /private7/projects/Combinatorics"
# mpileup_cmd = f"""samtools \
# mpileup \
# --fasta-ref D.pealeii/Annotations/orfs_squ.fa \
# --region {region}:{start+1}-{end} \
# --no-output-ins --no-output-ins \
# --no-output-del --no-output-del \
# --no-output-ends \
# --output-QNAME \
# --excl-flags 2304 \
# --max-depth 200000 \
# D.pealeii/Alignment/BestN1/PCLO-CNS-RESUB.C0x1291.aligned.sorted.bam \
# --output D.pealeii/Alignment.Test.N1/PCLO.FullORF.mpileup"""
# # !{cd_cmd} && {mpileup_cmd}


# %%
known_sites_file = (
    "/private7/projects/Combinatorics/D.pealeii/Annotations/D.pea.EditingSites.bed"
)
# problamatic_regions_file = "/private7/projects/Combinatorics/D.pealeii/Alignment/BestN1/CLO-CNS-RESUB.C0x1291.ProbRegions.bed"
problamatic_regions_file = None

pileup_file = (
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.pileup"
)


# %% [markdown]
# # Positions

# %% [markdown]
# ## Problematic regions & known editing sites

# %%
def positions_df_to_bed(positions_df: pd.DataFrame, strand: str) -> BedTool:
    """Convert positions DataFrame to a 6-columns BedTool.

    Args:
        positions_df (pd.DataFrame): A DataFrame of positions created from the pileup file.
        strand (str): The ORF's strand.

    Returns:
        BedTool: A BedTool of positions.
    """
    positions_bed = BedTool.from_dataframe(
        positions_df.loc[:, ["Chrom", "Position"]].assign(
            EndPosition=positions_df.iloc[:, 1] + 1, Name=".", Score=".", Strand=strand
        )
    )
    return positions_bed



# %%
def annotate_prolamatic_sites(positions_df, strand, problamatic_regions_file=None):
    if problamatic_regions_file:
        problamatic_regions_bed = BedTool(problamatic_regions_file)
        positions_bed = positions_df_to_bed(positions_df, strand)
        prob_positions_bed = positions_bed.intersect(problamatic_regions_bed)
        prob_positions = {interval.start for interval in prob_positions_bed}

        def is_prob(position: int) -> bool:
            return position in prob_positions

        in_prob_region = [is_prob(position) for position in positions_df["Position"]]
    else:
        in_prob_region = [False for _ in range(len(positions_df))]
    positions_df.insert(
        positions_df.columns.get_loc("Position") + 1, "InProbRegion", in_prob_region
    )



# %%
def annotate_known_sites(positions_df, strand, known_sites_file=None):
    if known_sites_file:
        known_sites_bed = BedTool(known_sites_file)
        positions_bed = positions_df_to_bed(positions_df, strand)
        known_positions_bed = positions_bed.intersect(known_sites_bed, s=True)
        known_positions = {interval.start for interval in known_positions_bed}

        def is_known(position: int) -> bool:
            return position in known_positions

        is_position_known = [
            is_known(position) for position in positions_df["Position"]
        ]
    else:
        is_position_known = [False for _ in range(len(positions_df))]
    positions_df.insert(
        positions_df.columns.get_loc("Position") + 1, "KnownEditing", is_position_known
    )



# %% [markdown]
# ## Currently-edited sites

# %%
def count_bases(ref_base, mapped_bases):
    mapped_bases = mapped_bases.replace(".", ref_base)
    return (
        mapped_bases.count("A"),
        mapped_bases.count("T"),
        mapped_bases.count("C"),
        mapped_bases.count("G"),
    )


def annotate_base_counts(positions_df):
    atcgs = positions_df.apply(
        lambda x: count_bases(x["RefBase"], x["MappedBases"]),
        axis=1,
        result_type="expand",
    ).rename(columns={0: "A", 1: "T", 2: "C", 3: "G"})
    positions_df[["A", "T", "C", "G"]] = atcgs



# %%
def position_noise_level(
    ref_base: str, strand: str, a_count: int, t_count: int, c_count: int, g_count: int
) -> float:
    """
    Calculate noise level in a position.

    The noise is defined as number of bases of the must abundant alt_base, divided by the same number +
    the number of ref_base bases.

    Args:
        ref_base (str): The reference base of the position.
        strand (str): The ORF's strand.
        a_count (int): Number of A's mapped to position.
        t_count (int): Number of T's mapped to position.
        c_count (int): Number of C's mapped to position.
        g_count (int): Number of G's mapped to position.

    Returns:
        float: The noise level.
    """
    base_counts = [a_count, t_count, c_count, g_count]

    if sum(base_counts) == 0:
        noise = np.NaN

    # we only measure noise in positions that don't undergo RNA editing by ADAR
    elif (strand == "+" and ref_base == "A") or (strand == "-" and ref_base == "T"):
        noise = np.NaN

    else:
        # we measure noise for T positions only on the positive strand
        if strand == "+" and ref_base == "T":
            ref_base_count = t_count
        # we measure noise for A positions only on the negative strand
        elif strand == "-" and ref_base == "A":
            ref_base_count = a_count
        # we measure noise for A & G positions on both strand
        elif ref_base == "C":
            ref_base_count = c_count
        else:  # ref_base == "G"
            ref_base_count = g_count

        bases = ["A", "T", "C", "G"]
        alt_bases = set(bases) - {ref_base}
        alt_base_counts = [
            base_count
            for base_count, base in zip(base_counts, bases)
            if base in alt_bases
        ]
        max_alt_base_count = max(alt_base_counts)

        try:
            noise = max_alt_base_count / (max_alt_base_count + ref_base_count)
        except ZeroDivisionError:
            noise = 0  # if there are no mapped alt bases

    return noise


def annotate_noise(positions_df, strand):
    positions_df["Noise"] = positions_df.apply(
        lambda x: position_noise_level(
            x["RefBase"], strand, x["A"], x["T"], x["C"], x["G"]
        ),
        axis=1,
    )



# %%
def editing_frequency_per_position(
    ref_base: str, base: str, ref_base_count: int, alt_base_count: int
) -> float:
    if base == ref_base:
        try:
            freq = alt_base_count / (ref_base_count + alt_base_count)
        except ZeroDivisionError:
            freq = 0
    else:
        freq = np.NaN
    return freq


def annotate_editing_frequency_per_position(positions_df: pd.DataFrame, strand: str):
    principal_ref_base = "A" if strand == "+" else "T"
    principal_alt_base = "G" if strand == "+" else "C"

    positions_df.insert(
        len(positions_df.columns),
        "EditingFrequency",
        positions_df.apply(
            lambda x: editing_frequency_per_position(
                principal_ref_base,
                x["RefBase"],
                x[principal_ref_base],
                x[principal_alt_base],
            ),
            axis=1,
        ),
    )



# %%
def annotate_edited_sites(
    positions_df: pd.DataFrame,
    strand: str,
    max_noise_level: float,
    assurance_factor: float,
):
    """Find currently-edited sites.

    For each position, check if editing_frequency > max_noise_level * assurance_factor.

    Args:
        positions_df (pd.DataFrame): The positions dataframe.
        strand (str): The ORF's strand.
        max_noise_level (float): The measured max noise level.
        assurance_factor (float): Set max_noise_level * assurance_factor to be the actual max noise level.
    """
    max_noise_level *= assurance_factor
    ref_base = "A" if strand == "+" else "T"

    edited_positions = positions_df.loc[positions_df["RefBase"] == ref_base].apply(
        lambda x: x["EditingFrequency"] > max_noise_level, axis=1
    )
    edited = [
        i in edited_positions.loc[edited_positions].index for i in positions_df.index
    ]
    positions_df.insert(positions_df.columns.get_loc("Position") + 1, "Edited", edited)



# %% [markdown]
# ## Pileup to positions (finally)

# %%
def pileup_to_positions(
    pileup_file: Path,
    strand: str,
    assurance_factor: float,
    problamatic_regions_file: Union[Path, str, None] = None,
    known_sites_file: Union[Path, str, None] = None,
    positions_out_file: Union[Path, str, None] = None,
) -> pd.DataFrame:
    """Read pileup file into a DataFrame.
    Verify format and change to 0-based coordinates.
    Possibly remove unedited positions, and possibly mask positions in problamatic regions.

    Args:
        pileup_file (Path): Path to a file creates by the `mpileup` function.
        strand (str): The ORF's strand.
        assurance_factor (float): Multiply the measured max_noise_level by this factor to define the actual max noise level.
        problamatic_regions_file (Union[Path, str, None]): If given, annotate positions as residing in problamatic genomic regions.
        known_sites_file (Union[Path, str, None]): If given, annotate positions as known editing sites.
        positions_out_file (Union[Path, str, None]): If given, write positions_df to this file.

    Raises:
        ValueError: If the 5th col doens't contain only base match/ bash mismatch characters, and/or the number of
        reads' names in the 7th col isn't equal to the number of characters in the 5th col.

    Returns:
        pd.DataFrame: A df of the parsed pileup file.
    """
    cols = [
        "Chrom",
        "Position",
        "RefBase",
        "TotalCoverage",
        "MappedBases",
        "Phred",
        "Reads",
    ]
    positions_df = pd.read_csv(pileup_file, sep="\t", names=cols)
    del positions_df["Phred"]
    # verify that the 5th col contains only base match/mismatch, and correspoding number
    # of reads' names in the 7th col
    for row_num, row in enumerate(positions_df.itertuples()):
        bases = row.MappedBases
        reads = row.Reads.split(",")
        if len(bases) != len(reads):
            raise ValueError(
                f"Line {row_num} in {pileup_file} contains indels and/or doesn't have the "
                "same number of mapped bases and reads."
            )
    # change editing position to 0-based
    positions_df["Position"] = positions_df["Position"] - 1

    # # present all reads as if they were on their positive strand
    positions_df["MappedBases"] = positions_df["MappedBases"].str.upper()
    positions_df["MappedBases"] = positions_df["MappedBases"].replace(
        {r"[,]": "."}, regex=True
    )

    unique_mapped_bases = set(chain.from_iterable(positions_df["MappedBases"]))
    assert unique_mapped_bases <= {"*", ".", "A", "C", "G", "T", "N"}

    annotate_prolamatic_sites(positions_df, strand, problamatic_regions_file)
    annotate_known_sites(positions_df, strand, known_sites_file)
    annotate_base_counts(positions_df)
    annotate_noise(positions_df, strand)
    max_noise_level = positions_df["Noise"].max()
    annotate_editing_frequency_per_position(positions_df, strand)
    annotate_edited_sites(positions_df, strand, max_noise_level, assurance_factor)

    # verify that noise is only applied to non ref_base positions, and vice versa for editing frequency
    ref_base = "A" if strand == "+" else "T"
    non_ref_base_edit_freq = positions_df.loc[
        positions_df["RefBase"] != ref_base, "EditingFrequency"
    ].unique()
    assert len(non_ref_base_edit_freq) == 1 and np.isnan(non_ref_base_edit_freq[0])
    ref_base_noise = positions_df.loc[
        positions_df["RefBase"] == ref_base, "Noise"
    ].unique()
    assert len(ref_base_noise) == 1 and np.isnan(ref_base_noise[0])
    
    if positions_out_file:
        positions_df.to_csv(positions_out_file, sep=",", index=False, na_rep=np.NaN)  

    return positions_df



# %%
positions_df = pileup_to_positions(
    pileup_file,
    strand,
    assurance_factor,
    problamatic_regions_file=problamatic_regions_file,
    known_sites_file=known_sites_file,
)
positions_df


# %%
import plotly.express as px


# %%
fig = px.histogram(
    positions_df,
    x="EditingFrequency",
    color="Edited",
    pattern_shape="KnownEditing",
    nbins=100,
    title="Editing frequency distribution",
    log_y=True,
)

fig.update_layout(
    xaxis=dict(
        tickmode="linear",
        tick0=min(positions_df["EditingFrequency"]),
        # dtick = 0.025
        dtick=0.01,
    )
)

fig.show()


# %%
df = positions_df.loc[positions_df["RefBase"] != "A"]
# df["Position"] = df["Position"].astype("str")
fig = px.bar(
    df,
    x="Position",
    y="Noise",
    template="simple_white",
    color="RefBase",
    color_discrete_sequence=["black"],
    title="Noise per position",
)
fig.update_layout(showlegend=False)
fig.show()


# %%
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn2

labels = ["Edited", "KnownEditing", "InProbRegion"]
sets = [set(positions_df.loc[positions_df[label], "Position"]) for label in labels]

if len(sets[2]) == 0:
    labels = labels[:2]
    sets = sets[:2]
    title = "Positions' membership: currently edited & known editing"
    v_func = venn2
    figsize = (3.5, 3.5)
else:
    title = "Positions' membership: currently edited, known editing, and probalamtic regions"
    v_func = venn3
    figsize = (5, 5)

plt.figure(figsize=figsize)
ax = plt.gca()
v = v_func(sets, set_labels=labels, ax=ax)
plt.title(title)
plt.show()


# %%
positions_df.loc[positions_df["Edited"] & ~positions_df["InProbRegion"]]


# %% [markdown]
# # Reads

# %%
def count_values(row, first_pos_loc, value):
    return len(row.iloc[first_pos_loc:].loc[row.iloc[first_pos_loc:] == value])

def annotate_edited_positions(df, first_pos_loc):
    df.insert(
        first_pos_loc,
        "EditedPositions",
        df.apply(lambda row: count_values(row, first_pos_loc, 1), axis="columns")
    )

def annotate_unedited_positions(df, first_pos_loc):
    df.insert(
        first_pos_loc,
        "UneditedPositions",
        df.apply(lambda row: count_values(row, first_pos_loc, 0), axis="columns")
    )

def annotate_ambigous_positions(df, first_pos_loc):
    df.insert(
        first_pos_loc,
        "AmbigousPositions",
        df.apply(lambda row: count_values(row, first_pos_loc, -1), axis="columns")
    )


# %%
def annotate_editing_frequency_per_read(reads_df, new_col_pos, positions_df, strand):
    ref_base = "A" if strand == "+" else "T"
    all_refbase_positions_df = (
        positions_df.loc[
            (positions_df["RefBase"] == ref_base) & 
            (~positions_df["InProbRegion"])
        ]
    )
    all_refbase_positions = len(all_refbase_positions_df)
    reads_df.insert(
        new_col_pos,
        "EditingFrequency",
        reads_df.apply(lambda row: row["EditedPositions"] / all_refbase_positions, axis="columns")
    )


# %%
def positions_to_reads(
    positions_df,
    strand,
    group_col,
    group,
):
    # retain only edited & reliable positions
    edited_positions_df = positions_df.loc[
        positions_df["Edited"] & ~positions_df["InProbRegion"]
    ]

    unique_reads = set(chain.from_iterable(edited_positions_df["Reads"].str.split(",")))
    positions_per_read = {read: [] for read in unique_reads}

    alt_base = "G" if strand == "+" else "C"
    for x, row in enumerate(edited_positions_df.itertuples()):
        bases = row.MappedBases
        reads = row.Reads.split(",")
        # deal with reads that were mapped to this pos
        for mapped_base, read in zip(bases, reads):
            # 1 = A2G/T2C editing, 0 = match
            if mapped_base == alt_base:
                # pos_edited_in_read = "1" # A2G/T2C editing
                pos_edited_in_read = 1  # A2G/T2C editing
            elif mapped_base == ".":
                # pos_edited_in_read = "0" # match
                pos_edited_in_read = 0  # match
            else:
                # pos_edited_in_read = "-1" # we ignore mismatches that aren't RNA editing by marking them as NaN
                pos_edited_in_read = (
                    -1
                )  # we ignore mismatches that aren't RNA editing by marking them as NaN
            positions_per_read[read].append(pos_edited_in_read)
        # deal with reads that weren't mapped to this pos
        unmapped_reads_in_pos = unique_reads - set(reads)
        for read in unmapped_reads_in_pos:
            # pos_edited_in_read = "-1"
            pos_edited_in_read = -1
            positions_per_read[read].append(pos_edited_in_read)
        # check all reads got values for pos
        mapped_pos_dist = {len(bases) for bases in positions_per_read.values()}
        if len(mapped_pos_dist) != 1:
            raise Exception(
                f"Problem at line {x} in positions_df: not all reads are mapped."
            )

    reads_df = pd.DataFrame(positions_per_read)
    reads_df = reads_df.T.reset_index().rename({"index": "Read"}, axis="columns")
    reads_df = reads_df.rename(
        columns={
            old_col: pos
            for old_col, pos in zip(reads_df.columns[1:], edited_positions_df["Position"])
        }
    )
    reads_df.insert(0, group_col, group)
    annotate_edited_positions(reads_df, 2)
    annotate_unedited_positions(reads_df, 3)
    annotate_ambigous_positions(reads_df, 4)
    annotate_editing_frequency_per_read(reads_df, 2, positions_df, strand)
    return reads_df



# %%
positions_df.loc[positions_df["Edited"]].iloc[:10, :]

# %%
reads_df = positions_to_reads(positions_df, strand, group_col, group)
# reads_df = positions_to_reads(positions_df.loc[positions_df["Edited"]].iloc[:10, :], strand, group_col, group)
reads_df


# %%
# find_uninformative_cols(reads_df, reads_df.columns[6:])

# %%
df = reads_df.assign(EditedOrAmbigousPositions = reads_df["EditedPositions"] + reads_df["AmbigousPositions"])
df

# %% tags=[]
from plotly.subplots import make_subplots
import plotly.graph_objects as go

df = reads_df.assign(EditedOrAmbigousPositions = reads_df["EditedPositions"] + reads_df["AmbigousPositions"])

x_title = "Edited sites per read"
y_title = "Reads"
title_text = "Distribution of min & max estimates of edited sites per read"

fig = make_subplots(
    # rows=1,
    rows=2,
    cols=1,
    # subplot_titles=conditions,
    shared_xaxes=True,
    x_title=x_title,
    y_title=y_title,
)

min_x = None
max_x = 0

col_names = ["EditedPositions", "EditedOrAmbigousPositions"]
estimate_names = ["Min", "Max"]

for i, (col_name, estimate_name) in enumerate(zip(col_names, estimate_names)):

    x = df[col_name]

    fig.add_trace(
        go.Histogram(
            x=x,
            # marker_color=subcolors_discrete_map[condition][i],
            name=estimate_name,
        ),
        # row=1,
        row=i + 1,
        col=1,
    )

    min_x = min(min_x, x.min()) if min_x else x.min()
    max_x = max(max_x, x.max())

fig.update_layout(
    # template=template,
    barmode="overlay",  # Overlay both histograms
    title_text=title_text,
    legend_title_text="Estimate",
)

# fig.update_traces(opacity=0.75)  # Reduce opacity to see both histograms
fig.update_xaxes(range=[min_x * 0.9, max_x * 1.1])

fig.show()


# %%


fig = px.histogram(
    reads_df,
    x="EditedPositions",
    labels={"EditedPositions": "Edited sites in read"},
    title="Distribution of edited sites in all reads",
)

# https://stackoverflow.com/questions/58167028/single-axis-caption-in-plotly-express-facet-plot
for axis in fig.layout:
    if type(fig.layout[axis]) == go.layout.YAxis:
        fig.layout[axis].title.text = ""
fig.update_layout(showlegend=False, yaxis_title="Reads")

fig.show()


# %% [markdown]
# # Unique reads ("transcripts")

# %%
def reads_to_unique_reads(reads_df, group_col, first_pos_loc=6):
    # data per read -> data per unique reads, and their supporting individual reads
    transcripts_df = reads_df.copy()

    transcripts_df.insert(
        1,
        "Reads",
        (
            transcripts_df.groupby(transcripts_df.columns[first_pos_loc:].tolist())[
                "Read"
            ].transform(lambda x: ",".join(x))
        ),
    )
    del transcripts_df["Read"]

    transcripts_df.insert(
        2, "NumOfReads", transcripts_df["Reads"].apply(lambda x: len(x.split(",")))
    )
    first_pos_loc += 1

    transcripts_df = (
        transcripts_df.loc[
            ~transcripts_df.duplicated(subset=transcripts_df.columns[first_pos_loc:])
        ]
        .sort_values("NumOfReads", ascending=False)
        .reset_index(drop=True)
    )

    transcripts_df.insert(
        1,
        "Transcript",
        pd.Series(
            transcripts_df[group_col] + "-" + transcripts_df.index.values.astype(str)
        ),
    )

    return transcripts_df



# %%
transcripts_df = reads_to_unique_reads(reads_df, group_col)
transcripts_df


# %%
import plotly.graph_objects as go

fig = px.histogram(
    transcripts_df,
    x="EditedPositions",
    labels={"EditedPositions": "Edited sites in read"},
    title="Distribution of edited sites in unique reads",
)

# https://stackoverflow.com/questions/58167028/single-axis-caption-in-plotly-express-facet-plot
for axis in fig.layout:
    if type(fig.layout[axis]) == go.layout.YAxis:
        fig.layout[axis].title.text = ""
fig.update_layout(showlegend=False, yaxis_title="Unique reads")

fig.show()


# %% [markdown]
# # Proteins

# %% [markdown]
# The goal is to have each 1-3 positions in the reads_df belong to an original AA.

# %%
def make_fasta_dict(fasta_file: Union[Path, str]) -> dict[str, Seq.Seq]:
    records = SeqIO.parse(Path(fasta_file).absolute(), "fasta")
    fasta_dict = {record.id: record.seq for record in records}
    return fasta_dict



# %% jupyter={"source_hidden": true} tags=[]
AA_BY_TRIPLET = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "AGT": "S",
    "AGC": "S",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGA": "*",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "TGT": "C",
    "TGC": "C",
    "TGG": "W",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGA": "R",
    "AGG": "R",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


# %%
StrOrNoneOrInt = Union[str, None, int]


# %% jupyter={"source_hidden": true} tags=[]
def filter_triplets(possible_triplets, base, pos):
    if base == -1:
        return possible_triplets
    else:
        return [triplet for triplet in possible_triplets if triplet[pos] == base]
    
def guided_translation(
    x1: str,
    y1: str,
    z1: str,
    x2: StrOrNoneOrInt,
    y2: StrOrNoneOrInt,
    z2: StrOrNoneOrInt,
):
    """
    Translate a triplet (x2, y2, z2) to an amino acid.

    The (x1, y1, z1) triplet is the original triplet from the transcriptome.

    The (x2, y2, z2) triplet is a corresponding triplet of a read.

    Use the original triplet as a template in case any of the bases in the new triplet is None.
    Return a set of all possible translations of amino acids in case any of the new bases equals -1.
    """
    if x2 is None:
        x2 = x1
    if y2 is None:
        y2 = y1
    if z2 is None:
        z2 = z1

    possible_triplets = list(AA_BY_TRIPLET.keys())

    for base, pos in zip([x2, y2, z2], [0, 1, 2]):
        possible_triplets = filter_triplets(possible_triplets, base, pos)

    amino_acids = {AA_BY_TRIPLET[triplet] for triplet in possible_triplets}
    amino_acids = sorted(list(amino_acids))
    amino_acids = ",".join(amino_acids)

    return amino_acids



# %%
def guided_translation_2(
    x1: str,
    y1: str,
    z1: str,
    x2: StrOrNoneOrInt,
    y2: StrOrNoneOrInt,
    z2: StrOrNoneOrInt,
):
    """
    Return all possible translations of the triplet (x2, y2, z2).

    The (x1, y1, z1) triplet is the original reference triplet.
    The (x2, y2, z2) triplet is a corresponding triplet of a read.

    Use the original triplet as a template in case any of the bases in the new triplet is None. This is the case for each reference base other than A.
    When an A base is marked as -1, it could be either A or G in the read.
    
    Return a set of all possible translations of amino acids in case any of the new bases equals -1.
    """
    if x2 is None:
        x2 = x1
    if y2 is None:
        y2 = y1
    if z2 is None:
        z2 = z1

    possible_triplets = list(AA_BY_TRIPLET.keys())

    # filter possible triplets according to the base in each position
    for base, pos in zip([x2, y2, z2], [0, 1, 2]):
        if base == -1:
            # possible_triplets = possible_triplets
            possible_bases_in_pos = ["A", "G"]
        else:
            # possible_triplets = [triplet for triplet in possible_triplets if triplet[pos] == base]
            possible_bases_in_pos = [base]
        possible_triplets = [
            triplet 
            for triplet in possible_triplets 
            if triplet[pos] in possible_bases_in_pos
        ]

    amino_acids = {AA_BY_TRIPLET[triplet] for triplet in possible_triplets}
    amino_acids = sorted(list(amino_acids))
    amino_acids = ",".join(amino_acids)

    return amino_acids



# %% [markdown]
# Q tests

# %%
x1, y1, z1 = ["C", "A", "A"]  # Q

x2, y2, z2 = [None, None, "G"]  # Q
ic(guided_translation(x1, y1, z1, x2, y2, z2))

x2, y2, z2 = [None, -1, "G"]  # L, P, Q, R
ic(guided_translation(x1, y1, z1, x2, y2, z2))

x2, y2, z2 = [None, -1, -1]  # L, P, H, Q, R
ic(guided_translation(x1, y1, z1, x2, y2, z2));

# %%
x1, y1, z1 = ["C", "A", "A"]  # Q

x2, y2, z2 = [None, None, "G"]  # Q
ic(guided_translation_2(x1, y1, z1, x2, y2, z2))

x2, y2, z2 = [None, -1, "G"]  # Q, R
ic(guided_translation_2(x1, y1, z1, x2, y2, z2))

x2, y2, z2 = [None, -1, -1]  # Q, R
ic(guided_translation_2(x1, y1, z1, x2, y2, z2));

# %% [markdown]
# Y tests

# %%
x1, y1, z1 = ["T", "A", "C"]  # Y

x2, y2, z2 = [None, "A", None]  # Y
ic(guided_translation(x1, y1, z1, x2, y2, z2))

x2, y2, z2 = [None, "G", None]  # C
ic(guided_translation(x1, y1, z1, x2, y2, z2))

x2, y2, z2 = [None, -1, None]  # C,F,S,Y
ic(guided_translation(x1, y1, z1, x2, y2, z2))

x2, y2, z2 = [None, -1, None]  # C,Y
ic(guided_translation_2(x1, y1, z1, x2, y2, z2));


# %%
def find_uninformative_cols(df, cols):
    return [col for col in cols if len(df[col].unique()) == 1]



# %%
def annotate_min_max_non_syns(df, transalted_orf_length, original_aas, first_aa_pos, new_cols_first_pos=None):
    def min_non_syns_in_row(row):
        return sum(
            [1 for aa, original_aa in zip(row, original_aas) if original_aa not in aa]
        )

    def max_non_syns_in_row(row):
        return sum(
            [1 for aa, original_aa in zip(row, original_aas) if aa != original_aa]
        )

    min_non_syns = df.iloc[:, first_aa_pos:].apply(min_non_syns_in_row, axis="columns")
    max_non_syns = df.iloc[:, first_aa_pos:].apply(max_non_syns_in_row, axis="columns")

    assert max_non_syns.ge(min_non_syns).all()
    
    min_non_syns_frequency = min_non_syns / transalted_orf_length
    max_non_syns_frequency = max_non_syns / transalted_orf_length

    new_cols_first_pos = new_cols_first_pos if new_cols_first_pos else first_aa_pos
    df.insert(new_cols_first_pos, "MinNonSyns", min_non_syns)
    df.insert(new_cols_first_pos + 1, "MaxNonSyns", max_non_syns)
    df.insert(new_cols_first_pos + 2, "MinNonSynsFrequency", min_non_syns_frequency)
    df.insert(new_cols_first_pos + 3, "MaxNonSynsFrequency", max_non_syns_frequency)



# %%
first_pos_loc = 8

positions_cols = set(transcripts_df.columns[first_pos_loc:])
# positions_cols

# %%
if strand == "+":
    ref_base, alt_base = "A", "G"
else:
    ref_base, alt_base = "T", "C"  # get complementary bases

# %%
position_to_base = {-1: -1, 0: ref_base, 1: alt_base}

# %%
transcriptome_dict = make_fasta_dict(transcriptome)
m_rna = transcriptome_dict[region][start:end]
original_triplets = [m_rna[x : x + 3] for x in range(0, len(m_rna), 3)]
if strand == "-":
    original_triplets = [triplet.reverse_complement() for triplet in original_triplets]
assert len(original_triplets[0]) == len(original_triplets[-1]) == 3
ic(original_triplets[0])

# %%
transalted_orf_length = len(original_triplets)

# %%
new_triplets_positions = [(x, x + 1, x + 2) for x in range(start, end, 3)]
ic(new_triplets_positions[0])

# %%
transcripts_dict = transcripts_df.to_dict(orient="records")

aas_per_unique_read = defaultdict(list)  # amino acids per unique read
triplets_cols = []
original_aas = []

for original_triplet, new_triplet_positions in zip(
    original_triplets, new_triplets_positions
):
    x1, y1, z1 = list(original_triplet)
    pos_x2, pos_y2, pos_z2 = new_triplet_positions
    pos_x2_exist = pos_x2 in positions_cols
    pos_y2_exist = pos_y2 in positions_cols
    pos_z2_exist = pos_z2 in positions_cols
    if not any([pos_x2_exist, pos_y2_exist, pos_z2_exist]):
        continue
    original_aa = AA_BY_TRIPLET[original_triplet]
    triplets_cols.append(f"{pos_x2}:{pos_x2 + 3}({original_aa})")
    original_aas.append(original_aa)
    for row in transcripts_dict:
        row_id = row["Transcript"]
        x2 = position_to_base[row[pos_x2]] if pos_x2_exist else None
        y2 = position_to_base[row[pos_y2]] if pos_y2_exist else None
        z2 = position_to_base[row[pos_z2]] if pos_z2_exist else None
        if strand == "-":
            # reverse the bases (they're already complementary due to position_to_base)
            x2, z2 = z2, x2
        # aa = guided_translation(x1, y1, z1, x2, y2, z2)
        aa = guided_translation_2(x1, y1, z1, x2, y2, z2)
        aas_per_unique_read[row_id].append(aa)

proteins_df = pd.DataFrame(aas_per_unique_read)
proteins_df = proteins_df.T.reset_index().rename(
    {"index": "Transcript"}, axis="columns"
)
proteins_df = proteins_df.set_axis(["Transcript"] + triplets_cols, axis="columns")
proteins_df = transcripts_df.iloc[:, :first_pos_loc].merge(proteins_df, on="Transcript")

# uninformative_cols = find_uninformative_cols(
#     proteins_df, proteins_df.columns[first_pos_loc:]
# )  # remove cols with a single AA
# proteins_df = proteins_df.drop(uninformative_cols, axis="columns")

annotate_min_max_non_syns(proteins_df, transalted_orf_length, original_aas, first_pos_loc, 1)

# uninformative_cols = find_uninformative_cols(
#     proteins_df, proteins_df.columns[first_pos_loc:]
# )  # remove cols with a single AA
# proteins_df = proteins_df.drop(uninformative_cols, axis="columns")


proteins_df

# %%

# %%

# %%
transcripts_df.loc[(transcripts_df["EditedPositions"] == 0) & (transcripts_df["AmbigousPositions"] == 0)].iloc[:, :8]

# %%
transcripts_df.loc[(transcripts_df["EditedPositions"] == 0) & (transcripts_df["AmbigousPositions"] == 0)].iloc[:, 8:]

# %%
proteins_df.loc[(proteins_df["EditedPositions"] == 0) & (proteins_df["AmbigousPositions"] == 0)].iloc[:, :12]

# %%
proteins_df.loc[(proteins_df["EditedPositions"] == 0) & (proteins_df["AmbigousPositions"] == 0)].iloc[:, 12:]

# %%

# %%

# %%

# %%
first_pos_loc

# %%
proteins_df["2214:2217(Y)"].unique()

# %%
proteins_df.loc[proteins_df["2214:2217(Y)"] == "C,F,S,Y"]

# %%
proteins_df.loc[:, ["MinNonSyns", "MaxNonSyns", "EditedPositions"]]

# %%
all(proteins_df["EditedPositions"] <= proteins_df["MinNonSyns"])

# %% jupyter={"source_hidden": true} tags=[]
# import plotly.graph_objects as go

# fig = go.Figure()
# fig.add_trace(go.Histogram(x=proteins_df["MinNonSyns"], name="Min"))
# fig.add_trace(go.Histogram(x=proteins_df["MaxNonSyns"], name="Max"))

# # Overlay both histograms
# fig.update_layout(
#     barmode="overlay",
#     title_text="Distribution of min & max estimates of non-syn mutations per protein",
#     xaxis_title_text="Non-syn mutations per protein",  # xaxis label
#     yaxis_title_text="Proteins",  # yaxis label
#     legend_title_text="Estimate",
# )
# # Reduce opacity to see both histograms
# fig.update_traces(opacity=0.75)
# # fig.update_xaxes(range=[0, max([proteins_df["MinNonSyns"].max(), proteins_df["MaxNonSyns"].max()])*1.05])
# fig.update_xaxes(
#     range=[
#         min([proteins_df["MinNonSyns"].min(), proteins_df["MaxNonSyns"].min()]) * 0.95,
#         max([proteins_df["MinNonSyns"].max(), proteins_df["MaxNonSyns"].max()]) * 1.05,
#     ]
# )
# fig.show()


# %% tags=[]
from plotly.subplots import make_subplots

x_title = "Non-syn mutations per protein"
y_title = "Proteins"
title_text = "Distribution of min & max estimates of non-syn mutations per protein"

fig = make_subplots(
    # rows=1,
    rows=2,
    cols=1,
    # subplot_titles=conditions,
    shared_xaxes=True,
    x_title=x_title,
    y_title=y_title,
)

min_x = None
max_x = 0


col_names = ["MinNonSyns", "MaxNonSyns"]
estimate_names = ["Min", "Max"]

for i, (col_name, estimate_name) in enumerate(zip(col_names, estimate_names)):

    x = proteins_df[col_name]

    fig.add_trace(
        go.Histogram(
            x=x,
            # marker_color=subcolors_discrete_map[condition][i],
            name=estimate_name,
        ),
        # row=1,
        row=i + 1,
        col=1,
    )

    min_x = min(min_x, x.min()) if min_x else x.min()
    max_x = max(max_x, x.max())

fig.update_layout(
    # template=template,
    barmode="overlay",  # Overlay both histograms
    title_text=title_text,
    legend_title_text="Estimate",
)

# fig.update_traces(opacity=0.75)  # Reduce opacity to see both histograms
fig.update_xaxes(range=[min_x * 0.9, max_x * 1.1])

fig.show()


# %%
# number of proteins without any non-syn mutations
proteins_df.loc[proteins_df["MinNonSyns"] == 0]


# %% [markdown]
# # Unique proteins

# %% tags=[]
first_pos_loc = 12

unique_proteins_df = proteins_df.drop("NumOfReads", axis=1)
first_pos_loc -= 1

unique_proteins_df.insert(
    unique_proteins_df.columns.get_loc("Transcript"),
    "Transcripts",
    (
        unique_proteins_df.groupby(unique_proteins_df.columns[first_pos_loc:].tolist())[
            "Transcript"
        ].transform(lambda x: ",".join(x))
    ),
)
del unique_proteins_df["Transcript"]

unique_proteins_df.insert(
    unique_proteins_df.columns.get_loc("Transcripts") + 1,
    "NumOfTranscripts",
    unique_proteins_df["Transcripts"].apply(lambda x: len(x.split(","))),
)
first_pos_loc += 1
unique_proteins_df.insert(
    unique_proteins_df.columns.get_loc("Reads"),
    "Reads2",
    (
        unique_proteins_df.groupby(unique_proteins_df.columns[first_pos_loc:].tolist())[
            "Reads"
        ].transform(lambda x: ",".join(x))
    ),
)
unique_proteins_df = unique_proteins_df.drop("Reads", axis="columns").rename(columns={"Reads2": "Reads"})
unique_proteins_df.insert(
    unique_proteins_df.columns.get_loc("Reads") + 1,
    "NumOfReads",
    unique_proteins_df["Reads"].apply(lambda x: len(x.split(",")))
)
first_pos_loc += 1

unique_proteins_df = (
    unique_proteins_df.loc[
        ~unique_proteins_df.duplicated(subset=unique_proteins_df.columns[first_pos_loc:])
    ]
    .sort_values("NumOfReads", ascending=False)
    .reset_index(drop=True)
)
assert proteins_df["NumOfReads"].sum() == unique_proteins_df["NumOfReads"].sum()
assert (
    sum(proteins_df["NumOfReads"] * proteins_df["MinNonSyns"]) == 
    sum(unique_proteins_df["NumOfReads"] * unique_proteins_df["MinNonSyns"])
)
assert (
    sum(proteins_df["NumOfReads"] * proteins_df["MaxNonSyns"]) == 
    sum(unique_proteins_df["NumOfReads"] * unique_proteins_df["MaxNonSyns"])
)

unique_proteins_df.insert(
    1,
    "Protein",
    pd.Series(
        unique_proteins_df[group_col]
        + "-"
        + unique_proteins_df.index.values.astype(str)
    ),
)

unique_proteins_df


# %%
unique_proteins_df.columns[:15]

# %%
# num_of_transcripts = sorted(list(unique_proteins_df["NumOfTranscripts"].unique()))
# fig = px.scatter(
#     unique_proteins_df,
#     x="MinNonSyns",
#     y="MaxNonSyns",
#     # size="NumOfReads",
#     facet_col="NumOfTranscripts",
#     facet_col_wrap=int(len(num_of_transcripts) ** 0.5),
#     category_orders={"NumOfTranscripts": num_of_transcripts},
#     facet_row_spacing=0.1, # default is 0.07 when facet_col_wrap is used
#     # facet_col_spacing=0.1,
#     color="NumOfReads",
#     title="Min vs Max non-synonymous mutations per protein",
#     template="plotly_white",
#     labels={"MinNonSyns": "Min", "MaxNonSyns": "Max"},
# )
# # fig.update_xaxes(rangeselector_font_size=0.5)
# # fig.update_yaxes(rangeselector_font_size=0.5)
# fig.show()
