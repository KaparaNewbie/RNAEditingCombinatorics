"""Native selected assignments and original RNA/AA inputs for manuscript producers.

The TSV names the already accepted solutions. This code does not select a new
solution or reassign expression. Run from a clone with the documented native data.
"""
from pathlib import Path
from manuscript_paths import ASSIGNMENTS_DIR
from collections import Counter
import hashlib
import os
import numpy as np
import pandas as pd

NOTEBOOKS = Path(__file__).resolve().parent
PROJECT = Path(os.environ.get('MANUSCRIPT_DATA_ROOT', NOTEBOOKS.parents[1])).resolve()
PARAMETERS = ['#Solution', 'Fraction', 'FractionRepetition', 'Algorithm', 'AlgorithmRepetition']
SUPPORT = ['NumOfReads', 'AdditionalWeightedSupportingReads', 'TotalWeightedSupportingReads']
SCALAR_FIELDS = ['Gene', 'Protein', *PARAMETERS, *SUPPORT, 'MeanEditedPositions',
                 'SufficientEntropy', 'MinNonSyns', 'MaxNonSyns', 'MinNonSynsFrequency',
                 'MaxNonSynsFrequency', 'MeanUneditedPositions', 'MeanNAPositions',
                 'AdditionalSupportingProteins']


def selections():
    """Return the manuscript's fixed gene/cohort order and native source identities."""
    table = pd.read_table(NOTEBOOKS / 'manuscript_input_data/selected_solutions.tsv', dtype={'Solution': str})
    assert not table.duplicated(['Platform', 'Gene', 'Metric']).any()
    assert table.Order.tolist() == list(range(25))
    return table


def ids_sha256(values):
    return hashlib.sha256(('\n'.join(values) + '\n').encode()).hexdigest()


def verified_file(root, name, expected):
    path = Path(root) / name
    with path.open('rb') as handle:
        assert hashlib.file_digest(handle, 'sha256').hexdigest() == expected, path
    return path


def source(path, expected_sha256=None):
    """Resolve an input, optionally checking its documented version; no audit writes."""
    path = Path(path).resolve()
    assert path.is_file(), path
    if expected_sha256 is not None:
        return verified_file(path.parent, path.name, expected_sha256)
    return path


def selection_bindings(metric='ordinary'):
    """Small compatibility view of the fixed scientific selection parameters."""
    result = {}
    for selection in selections().loc[lambda t: t.Metric.eq(metric)].itertuples(index=False):
        result.setdefault(selection.Platform, {})[selection.Gene] = {
            'parameters': {'#Solution': selection.Solution,
                **{key: getattr(selection, key) for key in PARAMETERS[1:]}},
            'solution_id': selection.Solution, 'maximum_cardinality': selection.SelectedProteins,
            'expression_realpath': str(PROJECT / selection.ExpressionFile),
            'expression_sha256': selection.ExpressionSHA256,
            'canonical_file': f'{selection.Platform}.{selection.Gene}.{selection.Metric}.tsv',
            'canonical_sha256': selection.AssignmentSHA256,
            'ordered_native_ids_sha256': selection.ProteinOrderSHA256}
    return result


def read_selected_assignment(selection, project_root=PROJECT):
    """Preserve native parsing, descending rank and within-tie order exactly."""
    path = verified_file(project_root, selection.ExpressionFile, selection.ExpressionSHA256)
    header = pd.read_table(path, nrows=0).columns.tolist()
    columns = [column for column in SCALAR_FIELDS if column in header]
    assert set(['Gene', 'Protein', *PARAMETERS, *SUPPORT, 'MeanEditedPositions']) <= set(columns)
    pieces, counts = [], Counter()
    for chunk in pd.read_table(path, usecols=columns, dtype={'Protein': str, '#Solution': str}, chunksize=10000):
        allowed = chunk.loc[chunk.Algorithm.isin(['Ascending', 'Descending'])]
        counts.update(allowed['#Solution'].value_counts().to_dict())
        pieces.append(chunk.loc[chunk['#Solution'].eq(selection.Solution)])
    selected = pd.concat(pieces, ignore_index=True)
    parameters = {'#Solution': selection.Solution, **{key: getattr(selection, key) for key in PARAMETERS[1:]}}
    assert selected[PARAMETERS].drop_duplicates().to_dict('records') == [parameters]
    assert selected.Protein.is_unique
    assert len(selected) == max(counts.values()) == selection.SelectedProteins
    # Keep the accepted producer's default sort, including its tie order.
    selected = selected.sort_values('TotalWeightedSupportingReads', ascending=False).reset_index(drop=True)
    assert ids_sha256(selected.Protein.tolist()) == selection.ProteinOrderSHA256
    np.testing.assert_allclose(selected.NumOfReads + selected.AdditionalWeightedSupportingReads,
                               selected.TotalWeightedSupportingReads, rtol=0, atol=1e-9)
    selected['Gene'] = selection.Gene
    selected.insert(0, 'Platform', 'Long-reads' if selection.Platform == 'PacBio' else 'Short-reads')
    selected['#Protein'] = np.arange(1, len(selected) + 1)
    selected['%RelativeExpression'] = 100 * selected.TotalWeightedSupportingReads / selected.TotalWeightedSupportingReads.sum()
    selected['%CummulativeRelativeExpression'] = selected['%RelativeExpression'].cumsum()
    return selected


def selected_amino_acids(selection, assignment, project_root=PROJECT):
    """One row per selected protein; preserve native codon order and missing states."""
    path = verified_file(project_root, selection.UniqueProteinsFile, selection.UniqueProteinsSHA256)
    codons = pd.read_table(path, nrows=0).columns.tolist()[14:]
    assert all(':' in c and '(' in c and c.endswith(')') for c in codons)
    pieces = []
    for chunk in pd.read_table(path, usecols=['Protein', 'NumOfReads', *codons], dtype={'Protein': str},
                               keep_default_na=False, chunksize=20000):
        pieces.append(chunk.loc[chunk.Protein.isin(assignment.Protein)])
    proteins = pd.concat(pieces, ignore_index=True).set_index('Protein')
    assert proteins.index.is_unique and set(proteins.index) == set(assignment.Protein)
    proteins = proteins.loc[assignment.Protein]
    np.testing.assert_array_equal(proteins.NumOfReads, assignment.NumOfReads)
    return proteins[codons]


def original_rna_patterns(selection, assignment, project_root=PROJECT):
    """Join original read → unique RNA → selected protein, retaining multiplicities.

    Multiple original RNA patterns can support one protein. Reassigned support is
    absent from this population; NumOfReads is original direct support only.
    """
    proteins_path = verified_file(project_root, selection.UniqueProteinsFile, selection.UniqueProteinsSHA256)
    proteins = pd.read_table(proteins_path, usecols=['Protein', 'Transcripts', 'Reads', 'NumOfReads'],
        dtype={'Protein': str, 'Transcripts': str, 'Reads': str}, keep_default_na=False).set_index('Protein')
    assert proteins.index.is_unique
    proteins = proteins.loc[assignment.Protein]
    np.testing.assert_array_equal(proteins.NumOfReads, assignment.NumOfReads)
    rna_path = verified_file(project_root, selection.UniqueReadsFile, selection.UniqueReadsSHA256)
    sites = pd.read_table(rna_path, nrows=0).columns.tolist()[8:]
    assert all(site.isdigit() for site in sites)
    rna = pd.read_table(rna_path, usecols=['Transcript', 'Reads', 'NumOfReads', *sites],
        dtype={'Transcript': str, 'Reads': str, **dict.fromkeys(sites, 'int8')},
        keep_default_na=False).set_index('Transcript')
    assert rna.index.is_unique and rna.NumOfReads.gt(0).all()
    mapping = proteins.Transcripts.str.split(',').explode()
    assert mapping.is_unique and mapping.isin(rna.index).all()
    selected = rna.loc[mapping]
    np.testing.assert_array_equal(selected.Reads.str.count(',') + 1, selected.NumOfReads)
    seen = set()
    for protein, transcripts in mapping.groupby(level=0, sort=False):
        actual = [read for ids in rna.loc[transcripts, 'Reads'] for read in ids.split(',')]
        expected = proteins.loc[protein, 'Reads'].split(',')
        assert len(actual) == len(set(actual)) == len(expected) == proteins.loc[protein, 'NumOfReads']
        assert set(actual) == set(expected) and seen.isdisjoint(actual)
        seen.update(actual)
    states = selected[sites].to_numpy(dtype=np.int8)
    assert set(np.unique(states)) <= {-1, 0, 1}
    return {'Protein': np.asarray(mapping.index, dtype=str), 'UniqueRead': np.asarray(mapping, dtype=str),
            'columns': np.asarray(sites, dtype=str), 'states': states,
            'NumOfReads': selected.NumOfReads.to_numpy(dtype=np.int64)}


def read_assignment(selection, assignments_dir):
    """Read a cache only when its version, selected IDs and order are verified."""
    stem = f'{selection.Platform}.{selection.Gene}.{selection.Metric}'
    path = verified_file(assignments_dir, stem + '.tsv', selection.AssignmentSHA256)
    table = pd.read_table(path, dtype={'Protein': str, '#Solution': str}, float_precision='round_trip')
    assert len(table) == selection.SelectedProteins and table.Protein.is_unique
    assert ids_sha256(table.Protein.tolist()) == selection.ProteinOrderSHA256
    assert table['#Solution'].eq(selection.Solution).all()
    return table


def read_ordinary(platform, assignments_dir=ASSIGNMENTS_DIR):
    """Load native-produced ordinary assignments in the publication's gene order."""
    selected = selections()
    selected = selected.loc[selected.Platform.eq(platform) & selected.Metric.eq('ordinary')]
    assert len(selected) == (2 if platform == 'PacBio' else 19)
    return pd.concat([read_assignment(s, assignments_dir) for s in selected.itertuples(index=False)], ignore_index=True)


NUMERIC = SUPPORT + ['#Protein', '%RelativeExpression', '%CummulativeRelativeExpression', 'MeanEditedPositions']

def verified_merged_dataframe(platform, assignment_dfs):
    """Build afresh only after validating in-memory sources against the frozen files."""
    expected=read_ordinary(platform)
    actual=pd.concat(assignment_dfs,ignore_index=True).copy()
    actual['Gene']=actual.Gene.str.split('_').str[0].replace({'GRIA':'GRIA2'})
    actual['Protein']=actual.Protein.astype(str)
    actual['#Solution']=actual['#Solution'].astype(str)
    keys=['Gene','Protein']
    assert not actual.duplicated(keys).any()
    assert actual[keys].values.tolist()==expected[keys].values.tolist()
    assert actual['#Solution'].tolist()==expected['#Solution'].tolist()
    for column in ['Fraction','FractionRepetition','Algorithm','AlgorithmRepetition']:
        assert actual[column].tolist()==expected[column].tolist(), column
    np.testing.assert_allclose(actual[NUMERIC],expected[NUMERIC],rtol=0,atol=1e-9,equal_nan=True)
    return expected


def read_merged(platform):
    name='AssignedExpression.PacBio.tsv' if platform=='PacBio' else 'AssignedExpression.Fixed.Illumina.tsv'
    actual=pd.read_table(ASSIGNMENTS_DIR / name,dtype={'Protein':str,'#Solution':str},float_precision='round_trip')
    expected=read_ordinary(platform)
    pd.testing.assert_frame_equal(actual,expected,check_dtype=False)
    return actual


def export_assignments(output_dir, project_root=PROJECT):
    """Generate the 25 fixed assignment tables in a new directory."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=False)
    ordinary = {'PacBio': [], 'Illumina': []}
    for selection in selections().itertuples(index=False):
        table = read_selected_assignment(selection, project_root)
        stem = f'{selection.Platform}.{selection.Gene}.{selection.Metric}'
        table.to_csv(output_dir / f'{stem}.tsv', sep='\t', index=False, na_rep='NA', mode='x')
        if selection.Metric == 'ordinary': ordinary[selection.Platform].append(table)
        print(stem, len(table), flush=True)
    for platform, parts in ordinary.items():
        name = 'AssignedExpression.PacBio.tsv' if platform == 'PacBio' else 'AssignedExpression.Fixed.Illumina.tsv'
        pd.concat(parts, ignore_index=True).to_csv(output_dir / name, sep='\t', index=False, na_rep='NA', mode='x')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True, help='New directory; existing outputs are never replaced')
    parser.add_argument('--project-root', type=Path, default=PROJECT)
    args = parser.parse_args()
    export_assignments(args.output, args.project_root)
