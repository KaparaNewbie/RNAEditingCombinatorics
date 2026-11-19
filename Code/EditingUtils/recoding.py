"""Binary recoding features for selected isoforms, one column per codon."""
import numpy as np
import pandas as pd

AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY*")


def binary_recoding_state(state, reference_aa):
    """0: reference only; 1: reference impossible; NA: recoding unresolved."""
    assert reference_aa in AMINO_ACIDS
    if pd.isna(state) or state == "":
        return np.nan
    possible_aas = set(state.split(","))
    assert possible_aas and possible_aas <= AMINO_ACIDS, state
    if possible_aas == {reference_aa}:
        return 0.0
    if reference_aa not in possible_aas:
        return 1.0
    return np.nan


def _binary_recoding_matrix(amino_acids, impute_nans):
    """Preserve the existing column-mean imputation and positive-frequency filter."""
    encoded = {}
    for codon in amino_acids:
        reference_aa = codon.split("(")[1][0]
        states = amino_acids[codon].unique()
        mapping = {state: binary_recoding_state(state, reference_aa) for state in states}
        encoded[codon] = amino_acids[codon].map(mapping).astype(float)
    encoded = pd.DataFrame(encoded, index=amino_acids.index)
    if impute_nans:
        encoded = encoded.fillna(encoded.mean())
    retained = encoded.where(encoded.isin([0.0, 1.0])).mean().gt(0)
    return encoded.loc[:, retained]


def prepare_recoding_input(assignment_df, unique_proteins_df, first_aa_column,
                          metadata_columns, relative_expression_column,
                          impute_nans=True, use_sklearn_standard_scaler=False):
    """Join selected isoforms and prepare their binary recoding features."""
    keys = metadata_columns[:2]
    codons = unique_proteins_df.columns[first_aa_column:].tolist()
    joined = assignment_df.loc[:, metadata_columns].merge(
        unique_proteins_df.loc[:, keys + codons], on=keys, how="left",
        validate="one_to_one", indicator=True)
    assert joined["_merge"].eq("both").all()
    matrix = _binary_recoding_matrix(joined[codons], impute_nans=impute_nans)
    if use_sklearn_standard_scaler:
        from sklearn.preprocessing import StandardScaler
        matrix.loc[:, :] = StandardScaler().fit_transform(matrix)
    prepared = pd.concat([joined[metadata_columns], matrix], axis=1)
    return prepared.sort_values(relative_expression_column, ascending=False,
                                ignore_index=True, kind="stable")
