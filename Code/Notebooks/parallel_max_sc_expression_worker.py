from ast import literal_eval
from gzip import BadGzipFile
from itertools import chain
from pathlib import Path
from tempfile import NamedTemporaryFile
from time import perf_counter

import pandas as pd
from pandas.errors import EmptyDataError


def parse_reads(value):
    """Return one normalized list of string read IDs."""
    if isinstance(value, (list, tuple)):
        return [str(read) for read in value]

    if pd.isna(value):
        return []

    text = str(value).strip()
    if not text:
        return []

    java_prefix = "SubString{String}["
    if text.startswith(java_prefix) and text.endswith("]"):
        text = text[len(java_prefix):-1]
        if not text:
            return []

        try:
            parsed = literal_eval(f"[{text}]")
        except (SyntaxError, ValueError):
            parsed = [
                read.strip().strip("\"'")
                for read in text.split(",")
            ]
        return [str(read) for read in parsed]

    if text.startswith("[") and text.endswith("]"):
        try:
            parsed = literal_eval(text)
        except (SyntaxError, ValueError):
            parsed = None

        if isinstance(parsed, (list, tuple)):
            return [str(read) for read in parsed]

    # Current SC expression files use either one unwrapped read ID or a
    # comma-separated sequence such as "X,U,2a".
    return [
        read.strip()
        for read in text.split(",")
        if read.strip()
    ]


def write_gzip_dataframe_atomically(dataframe, out_file, sep):
    """Publish a complete gzip file without exposing a partial write."""
    out_path = Path(out_file)

    with NamedTemporaryFile(
        dir=out_path.parent,
        prefix=f".{out_path.name}.",
        suffix=".tmp.gz",
        delete=False,
    ) as temporary_file:
        temporary_path = Path(temporary_file.name)

    try:
        dataframe.to_csv(
            temporary_path,
            sep=sep,
            index=False,
            compression="gzip",
        )
        temporary_path.replace(out_path)
    finally:
        temporary_path.unlink(missing_ok=True)


def get_f1_max_expression_df(
    expression_file,
    chrom,
    max_distinct_proteins_df,
    sep,
    condition_col,
):
    expression_df = pd.read_csv(
        expression_file,
        sep=sep,
        dtype={
            "#Solution": str,
            "Protein": str,
            "Reads": str,
            "AdditionalSupportingReadsIDs": str,
            "AdditionalSupportingProteinsIDs": str,
        },
    )

    expression_df = expression_df.loc[
        expression_df["Fraction"] == 1
    ].reset_index(drop=True)

    expression_df["Reads"] = expression_df["Reads"].apply(
        parse_reads
    )
    expected_num_reads = pd.to_numeric(
        expression_df["NumOfReads"],
        errors="coerce",
    )
    parsed_num_reads = expression_df["Reads"].apply(len)
    invalid_reads = (
        expected_num_reads.notna()
        & parsed_num_reads.ne(expected_num_reads)
    )
    if invalid_reads.any():
        invalid_index = invalid_reads[invalid_reads].index[0]
        raise ValueError(
            "Could not parse Reads consistently with NumOfReads for "
            f"{chrom=} in {expression_file}: "
            f"parsed {parsed_num_reads.loc[invalid_index]} read(s), "
            f"expected {expected_num_reads.loc[invalid_index]:g}"
        )

    # Expression files from different pipeline versions do not all contain
    # the position-summary columns. These columns are discarded when present,
    # but their absence should not prevent processing.
    optional_columns_to_drop = [
        "#Solution",
        "AmbigousPositions",
        "EditedPositions",
        "EditingFrequency",
        "Index",
        "NumOfUniqueReads",
        "Samples",
        "UneditedPositions",
        "UniqueReads",
        "AdditionalEqualSupportingReads",
        "TotalEqualSupportingReads",
        "MinNonSyns",
        "MaxNonSyns",
        "MinNonSynsFrequency",
        "MaxNonSynsFrequency",
        "AdditionalWeightedSupportingReads",
        "TotalWeightedSupportingReads",
    ]
    expression_df = expression_df.drop(
        columns=expression_df.columns.intersection(
            optional_columns_to_drop
        )
    )

    expression_df.insert(0, "Chrom", chrom)

    expression_df = expression_df.merge(
        max_distinct_proteins_df.loc[
            max_distinct_proteins_df["Chrom"] == chrom,
            [
                "Chrom",
                condition_col,
                "Fraction",
                "FractionRepetition",
                "Algorithm",
                "AlgorithmRepetition",
            ],
        ],
        how="right",
    )

    expression_df = expression_df.drop(
        columns=[
            "Fraction",
            "FractionRepetition",
            "Algorithm",
            "AlgorithmRepetition",
        ]
    )

    expression_df["AdditionalSupportingReadsIDs"] = expression_df[
        "AdditionalSupportingReadsIDs"
    ].apply(
        lambda value: (
            ""
            if pd.isna(value)
            else [reads.split(",") for reads in value.split(";")]
        )
    )

    expression_df["AdditionalSupportingProteinsIDs"] = expression_df[
        "AdditionalSupportingProteinsIDs"
    ].apply(
        lambda value: "" if pd.isna(value) else value.split(",")
    )

    return expression_df


def simplified_get_f1_exapnded_max_expression_df(
    chrom,
    expression_file,
    positions_dir,
    one_chrom_raw_reads_info_df,
    max_distinct_proteins_df,
    sep,
    condition_col,
    out_file=None,
    try_using_existing_out_file=True,
):
    if out_file is not None and try_using_existing_out_file:
        try:
            expanded_max_expression_df = pd.read_csv(
                out_file,
                sep=sep,
                dtype={"Protein": str, "Read": str},
            )

            if not expanded_max_expression_df.empty:
                return expanded_max_expression_df

        except (
            FileNotFoundError,
            EOFError,
            EmptyDataError,
            BadGzipFile,
        ):
            pass

    one_chrom_old_to_new_reads_file = Path(
        positions_dir,
        f"{chrom}.OldToNewReads.csv.gz",
    )

    one_chrom_old_to_new_reads_df = pd.read_table(
        one_chrom_old_to_new_reads_file,
        dtype={"OldRead": str, "NewRead": str},
    )

    one_chrom_new_reads_info_df = (
        one_chrom_raw_reads_info_df.merge(
            one_chrom_old_to_new_reads_df,
            how="left",
            left_on="ReadID",
            right_on="OldRead",
        )
        .drop(columns=["ReadID", "OldRead"])
        .rename(columns={"NewRead": "Read"})
    )

    max_expression_df = get_f1_max_expression_df(
        expression_file,
        chrom,
        max_distinct_proteins_df,
        sep,
        condition_col,
    )

    max_expression_df = max_expression_df.drop(
        columns=[
            "NumOfReads",
            "AdditionalSupportingProteinsIDs",
            "AdditionalSupportingProteins",
        ]
    )

    max_expression_df["AdditionalSupportingReadsIDs"] = (
        max_expression_df["AdditionalSupportingReadsIDs"].apply(
            lambda value: sorted(set(chain.from_iterable(value)))
        )
    )

    max_expression_df["AllReads"] = max_expression_df.apply(
        lambda row: (
            row["Reads"] + row["AdditionalSupportingReadsIDs"]
        ),
        axis=1,
    )

    max_expression_df["AllReadsStatuses"] = max_expression_df.apply(
        lambda row: (
            ["Original"] * len(row["Reads"])
            + ["Additional"]
            * len(row["AdditionalSupportingReadsIDs"])
        ),
        axis=1,
    )

    max_expression_df = max_expression_df.drop(
        columns=[
            "Reads",
            "AdditionalSupportingReadsIDs",
        ]
    )

    assert (
        max_expression_df["AllReads"]
        .apply(len)
        .eq(max_expression_df["AllReadsStatuses"].apply(len))
        .all()
    )

    expanded_max_expression_df = (
        max_expression_df.explode(
            ["AllReads", "AllReadsStatuses"],
            ignore_index=True,
        ).rename(
            columns={
                "AllReads": "Read",
                "AllReadsStatuses": "ReadStatus",
            }
        )
    )

    expanded_max_expression_df = expanded_max_expression_df.loc[
        ~(
            expanded_max_expression_df["Read"].eq("")
            & expanded_max_expression_df["ReadStatus"].eq(
                "Additional"
            )
        )
    ]

    assert (
        expanded_max_expression_df.shape[0]
        == expanded_max_expression_df.drop_duplicates().shape[0]
    )

    expanded_max_expression_df = (
        expanded_max_expression_df.merge(
            one_chrom_new_reads_info_df,
            how="left",
        )
    )

    if expanded_max_expression_df.empty:
        raise RuntimeError(
            f"The resulting DataFrame for {chrom=} is empty"
        )

    if out_file is not None:
        write_gzip_dataframe_atomically(
            expanded_max_expression_df,
            out_file,
            sep,
        )

    return expanded_max_expression_df


def process_chromosome(job):
    (
        index,
        chrom,
        expression_file,
        positions_dir,
        one_chrom_raw_reads_info_df,
        one_chrom_max_distinct_proteins_df,
        sep,
        condition_col,
        out_file,
        try_using_previous_out_file,
        strictly_use_previous_out_file_wo_verification,
    ) = job

    start = perf_counter()

    try:
        if strictly_use_previous_out_file_wo_verification:
            expanded_max_expression_df = pd.read_csv(
                out_file,
                sep=sep,
                dtype={"Protein": str, "Read": str},
            )
        else:
            expanded_max_expression_df = (
                simplified_get_f1_exapnded_max_expression_df(
                    chrom=chrom,
                    expression_file=expression_file,
                    positions_dir=positions_dir,
                    one_chrom_raw_reads_info_df=(
                        one_chrom_raw_reads_info_df
                    ),
                    max_distinct_proteins_df=(
                        one_chrom_max_distinct_proteins_df
                    ),
                    sep=sep,
                    condition_col=condition_col,
                    out_file=out_file,
                    try_using_existing_out_file=(
                        try_using_previous_out_file
                    ),
                )
            )
    except Exception as error:
        raise RuntimeError(
            f"Failed processing {chrom=} from {expression_file}"
        ) from error

    elapsed = perf_counter() - start
    return index, expanded_max_expression_df, elapsed
