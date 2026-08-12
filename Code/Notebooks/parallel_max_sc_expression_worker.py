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
            # compression="gzip",
            compression={
                "method": "gzip",
                "compresslevel": 1,
            }
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

    # The mapping must be one-to-one.
    if not one_chrom_old_to_new_reads_df["OldRead"].is_unique:
        raise ValueError(
            f"OldRead is not unique in the read-name mapping for {chrom=}"
        )

    if not one_chrom_old_to_new_reads_df["NewRead"].is_unique:
        raise ValueError(
            f"NewRead is not unique in the read-name mapping for {chrom=}"
        )

    # Each raw molecule/read should have one metadata record within the gene.
    if not one_chrom_raw_reads_info_df["ReadID"].is_unique:
        duplicated_raw_metadata = one_chrom_raw_reads_info_df.loc[
            one_chrom_raw_reads_info_df["ReadID"].duplicated(keep=False)
        ].sort_values("ReadID")

        raise ValueError(
            f"ReadID is not unique in the raw metadata for {chrom=}:\n"
            f"{duplicated_raw_metadata.head(20)}"
        )

    # Start from the mapping rather than from all raw reads.
    # This keeps only reads that received a shortened NewRead identifier.
    one_chrom_new_reads_info_df = (
        one_chrom_old_to_new_reads_df.merge(
            one_chrom_raw_reads_info_df,
            left_on="OldRead",
            right_on="ReadID",
            how="left",
            validate="one_to_one",
            indicator=True,
        )
    )

    # Every mapped read should have corresponding cell metadata.
    missing_metadata = one_chrom_new_reads_info_df["_merge"].ne("both")

    if missing_metadata.any():
        missing_old_reads = one_chrom_new_reads_info_df.loc[
            missing_metadata,
            "OldRead",
        ].head(20).tolist()

        raise ValueError(
            f"Mapped reads without raw metadata for {chrom=}: "
            f"{missing_old_reads}"
        )

    one_chrom_new_reads_info_df = (
        one_chrom_new_reads_info_df
        .drop(columns=["OldRead", "ReadID", "_merge"])
        .rename(columns={"NewRead": "Read"})
    )

    assert one_chrom_new_reads_info_df["Read"].notna().all()
    assert one_chrom_new_reads_info_df["Read"].is_unique

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

    # expanded_max_expression_df = (
    #     expanded_max_expression_df.merge(
    #         one_chrom_new_reads_info_df,
    #         how="left",
    #     )
    # )
    
    join_cols = ["Chrom", "Read"]

    # Every old/new read identifier should occur only once in the mapping.
    assert one_chrom_old_to_new_reads_df["OldRead"].is_unique
    assert one_chrom_old_to_new_reads_df["NewRead"].is_unique

    # The same read must not have more than one metadata record.
    duplicated_read_metadata = one_chrom_new_reads_info_df.duplicated(
        subset=join_cols,
        keep=False,
    )

    if duplicated_read_metadata.any():
        duplicated_rows = one_chrom_new_reads_info_df.loc[
            duplicated_read_metadata
        ].sort_values(join_cols)

        raise ValueError(
            f"More than one metadata record was found for a read in {chrom=}:\n"
            f"{duplicated_rows.head(20)}"
        )

    expanded_max_expression_df = expanded_max_expression_df.merge(
        one_chrom_new_reads_info_df,
        on=join_cols,
        how="left",
        validate="many_to_one",
    )
    
    # Detect expression reads for which the metadata merge found no match.
    metadata_cols = ["Sample", "CB", "UB"]

    missing_joined_metadata = (
        expanded_max_expression_df[metadata_cols]
        .isna()
        .all(axis=1)
    )

    if missing_joined_metadata.any():
        missing_reads = (
            expanded_max_expression_df.loc[
                missing_joined_metadata,
                ["Chrom", "Read", "ReadStatus"],
            ]
            .drop_duplicates()
            .head(20)
        )

        raise ValueError(
            f"Expression reads without cell metadata for {chrom=}:\n"
            f"{missing_reads}"
        )

    # The metadata merge must not create completely duplicated rows.
    assert not expanded_max_expression_df.duplicated().any()

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


def process_chromosome_to_file(job):
    """
    Process one chromosome and write its expanded DataFrame to disk.

    Return only lightweight metadata, rather than sending the complete
    DataFrame back to the parent process.
    """
    index, _expanded_df, elapsed = process_chromosome(job)
    return index, elapsed