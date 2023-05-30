import argparse
from pathlib import Path
from math import ceil, floor
from typing import Union
from itertools import chain

from Bio import SeqIO, motifs  # biopython
import pandas as pd  # pandas
import matplotlib.pyplot as plt  # matplotlib
from logomaker import Logo  # logomaker

from General.argparse_utils import abs_path_from_str, expanded_path_from_str


WIDTH = 9.6
HEIGHT = 5
MAX_COLS = 4
BASES_FROM_EACH_SIDE = 1
DPI = 300


def make_freq_df(fasta_file, bases_from_each_side=BASES_FROM_EACH_SIDE):
    """
    This function creates a frequency dataframe out of a fasta file, and returns it.
    """
    # get sequences from the fasta file
    records = SeqIO.parse(fasta_file, "fasta")
    seqs = [record.seq for record in records]
    # make sure all sequences are the same length
    assert len({len(seq) for seq in seqs}) == 1
    # create a motif object out of the sequences
    motif = motifs.create(seqs)
    # create a counts dataframe out of the motif object
    counts_df = pd.DataFrame({base: motif.counts[base] for base in "ACGT"})
    # transfrom the counts to frequencies
    # (all sequences have the same length, so total coverage in each position = len(seqs)
    # return the frequency dataframe)
    freq_df = counts_df / len(seqs)
    # retain only $bases_from_each_side from the middle of the motif
    middle_position = floor(len(freq_df) / 2)
    remaining_positions = (
        [x for x in range(middle_position - bases_from_each_side, middle_position)]
        + [middle_position]
        + [
            x
            for x in range(
                middle_position + 1, middle_position + bases_from_each_side + 1
            )
        ]
    )
    freq_df = freq_df.iloc[remaining_positions].reset_index(drop=True)
    # return the frequency dataframe
    return freq_df


# Plot a figure with multiple motifs of ADAR out of the frequency dfs.


def plot_single_logo(
    freq_df: pd.DataFrame,
    # title: str,
    main_title: Union[str, None],
    width: float = WIDTH,
    height: float = HEIGHT,
    # max_cols: int = MAX_COLS,
):
    fig, ax = plt.subplots(figsize=(width, height), gridspec_kw={"wspace": 0.5})

    # fig.suptitle(main_title, fontsize="x-large")
    # ax.set_title(main_title, pad=15)
    ax.set_title(main_title)
    plot = Logo(df=freq_df, ax=ax, vpad=0.02)

    # style plot
    # style plots using Logo methods
    plot.style_spines(visible=False)
    plot.style_spines(spines=("left", "bottom"), visible=True)
    # style plots using matplotlib's Axes methods
    # plot.ax.set_ylabel("Probability", labelpad=5)
    # plot.ax.set_xlabel("Position", labelpad=5)

    xticks = list(freq_df.index)
    plot.ax.set_xticks(xticks)
    middle_xtick = floor(len(xticks) / 2)
    xticklabels = (
        [str(x) for x in range(-middle_xtick, 0)]
        + ["0"]
        + [str(x) for x in range(1, middle_xtick + 1)]
    )
    plot.ax.set_xticklabels(xticklabels)
    fig.subplots_adjust(hspace=0.4)

    return fig


def plot_multiple_logos(
    freq_dfs: list[pd.DataFrame],
    titles: list[str],
    main_title: Union[str, None],
    width: float = WIDTH,
    height: float = HEIGHT,
    max_cols: int = MAX_COLS,
):
    """
    Plot a figure with multiple motifs of ADAR out of the frequency dfs.
    """
    # define fig, sub plots and titles
    num_of_plots = len(freq_dfs)

    # ncols = max_cols if max_cols >= 4 else num_of_plots
    # nrows = ceil(num_of_plots / ncols)

    # cols = min(facet_col_wrap, len(conditions), 4)
    # rows = ceil(len(conditions) / cols)

    ncols = min(max_cols if max_cols >= 4 else num_of_plots, num_of_plots)
    nrows = ceil(num_of_plots / ncols)

    fig, axes = plt.subplots(
        nrows=nrows, ncols=ncols, figsize=(width, height), gridspec_kw={"wspace": 0.5}
    )
    # axes = list(chain.from_iterable(axes))
    axes = axes.flat
    if main_title:
        fig.subplots_adjust(top=0.8)
        # fig.suptitle(main_title, y=1.0, fontsize="x-large")
        fig.suptitle(main_title, fontsize="x-large")
    plots = []
    for ax, title, freq_df in zip(axes, titles, freq_dfs):
        ax.set_title(title, pad=15)
        ax.set_title(title)
        plot = Logo(df=freq_df, ax=ax, vpad=0.02)

        # style plot
        # style plots using Logo methods
        plot.style_spines(visible=False)
        plot.style_spines(spines=("left", "bottom"), visible=True)
        # style plots using matplotlib's Axes methods
        # plot.ax.set_ylabel("Probability", labelpad=5)
        # plot.ax.set_xlabel("Position", labelpad=5)

        xticks = list(freq_df.index)
        plot.ax.set_xticks(xticks)
        middle_xtick = floor(len(xticks) / 2)
        xticklabels = (
            [str(x) for x in range(-middle_xtick, 0)]
            + ["0"]
            + [str(x) for x in range(1, middle_xtick + 1)]
        )
        plot.ax.set_xticklabels(xticklabels)

        plots.append(plot)

    # todo check if this helps to show x-axis title
    _, top_ylim = plt.ylim()

    # fig.tight_layout(h_pad=10, w_pad=10)
    fig.subplots_adjust(hspace=0.4)

    return fig


def multiple_logos_from_fasta_dir(
    fasta_dir: Path,
    postfix: str,
    main_title: Union[str, None],
    out_file: Union[Path, None],
    bases_from_each_side: int = BASES_FROM_EACH_SIDE,
    width: float = WIDTH,
    height: float = HEIGHT,
    max_cols: int = MAX_COLS,
    dpi: int = DPI,
):
    """
    Create multiple logos out of fasta files in `fasta_dir` and plot them together, side by side.

    Optionally, by providing a nonempty string for `output_file`, the plot can be also saved to disk.
    """

    fasta_dir = Path(fasta_dir).absolute()
    fasta_files = [f for f in fasta_dir.iterdir() if f.name.endswith(postfix)]
    samples_names = [f.name.removesuffix(postfix) for f in fasta_files]
    freq_dfs = [
        make_freq_df(fasta_file, bases_from_each_side) for fasta_file in fasta_files
    ]

    # plot the motif according to the frequency file, and save it to temp_dir
    fig = plot_multiple_logos(
        freq_dfs, samples_names, main_title, width, height, max_cols
    )
    if out_file:
        fig.savefig(Path(out_file).absolute(), dpi=dpi)
    return fig


def multiple_logos_from_fasta_files(
    fasta_files: Union[list[Path], list[str]],
    main_title: Union[str, None],
    sub_titles: Union[list[str], None],
    out_file: Union[Path, None],
    bases_from_each_side: int = BASES_FROM_EACH_SIDE,
    width: float = WIDTH,
    height: float = HEIGHT,
    max_cols: int = MAX_COLS,
    dpi: int = DPI,
):
    """
    Create multiple logos out of fasta files in `fasta_dir` and plot them together, side by side.

    Optionally, by providing a nonempty string for `output_file`, the plot can be also saved to disk.
    """

    sub_titles = sub_titles if sub_titles is not None else [""] * len(fasta_files)
    freq_dfs = [
        make_freq_df(fasta_file, bases_from_each_side) for fasta_file in fasta_files
    ]

    # plot the motif according to the frequency file, and save it to temp_dir
    if len(freq_dfs) > 1:
        fig = plot_multiple_logos(
            freq_dfs, sub_titles, main_title, width, height, max_cols
        )
    else:
        fig = plot_single_logo(freq_dfs[0], sub_titles[0], width, height)
    if out_file:
        fig.savefig(Path(out_file).absolute(), dpi=dpi)
    return fig


if __name__ == "__main__":
    # create parser
    parser = argparse.ArgumentParser()
    # define args
    parser.add_argument("fasta_dir", type=abs_path_from_str())
    parser.add_argument("out_file", type=abs_path_from_str())
    parser.add_argument("--postfix", default=".fa")
    parser.add_argument("--main_title")
    parser.add_argument("--width", type=float, default=WIDTH)
    parser.add_argument("--height", type=float, default=HEIGHT)
    parser.add_argument("--max_cols", type=int, default=MAX_COLS)
    parser.add_argument(
        "--bases_from_each_side", type=int, default=BASES_FROM_EACH_SIDE
    )
    parser.add_argument("--dpi", type=int, default=DPI)

    # parse args
    args = parser.parse_args()
    # run script with parsed args
    multiple_logos_from_fasta_dir(
        args.fasta_dir,
        args.postfix,
        args.main_title,
        args.out_file,
        args.bases_from_each_side,
        args.width,
        args.height,
        args.max_cols,
        args.dpi,
    )
