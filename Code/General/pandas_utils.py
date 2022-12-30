import pandas as pd


def reorder_df_by_wanted_cols(df: pd.DataFrame, wanted_first_cols: list) -> pd.DataFrame:
    """
    Reorder the columns in df such that columns in wanted_first_cols will be the first ones.
    The relative order of the other columns will remain unchanged.
    @param df: data frame
    @type df: pandas.DataFrame
    @param wanted_first_cols: columns to bring upfront (to the left)
    @type wanted_first_cols: list
    @return: df
    @rtype: pandas.DataFrame
    """
    wanted_first_cols_set = set(wanted_first_cols)
    all_cols_set = set(df.columns)
    other_cols_set = all_cols_set - wanted_first_cols_set
    other_cols = [col for col in df.columns if col in other_cols_set]
    new_cols = wanted_first_cols + other_cols
    df = df.reindex(new_cols, axis=1)
    return df
