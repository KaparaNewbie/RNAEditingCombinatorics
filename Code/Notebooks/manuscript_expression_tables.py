"""F4a, S8/S9 and S17 exports from verified native selected assignments.

One row is a selected isoform. This does not rerun expression reassignment.
S8 preserves the original zero-based training indices and log-space operations.
"""
from pathlib import Path
import argparse
import math
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from manuscript_inputs import selections, read_assignment

S8_INTERVALS = {
    'PacBio': [(173, 8402), (615, 17356)],
    'Illumina': [(2, 4305), (100, 12386), (943, 5475), (88, 31300),
        (103, 17500), (98, 21007), (140, 1142), (83, 5553), (107, 4218),
        (120, 13000), (34, 11500), (1044, 10000), (95, 7000), (192, 10000),
        (104, 28500), (99, 60000), (64, 30000), (254, 9354), (85, 19337)]}


def expression_tables(assignments_dir, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=False)
    ordinary = selections().loc[lambda t: t.Metric.eq('ordinary')]
    assignments, fits, equations, markers, observations = [], [], [], [], []
    for platform, selected in ordinary.groupby('Platform', sort=False):
        for facet, selection in enumerate(selected.itertuples(index=False)):
            table = read_assignment(selection, assignments_dir)
            label, gene = table.Platform.iloc[0], selection.Gene
            assignments.append(table)
            lo, hi = S8_INTERVALS[platform][facet]
            n_train = 400 if platform == 'Illumina' and facet == 6 else 1000
            train = [int(v) for v in np.logspace(np.log10(lo), np.log10(hi), num=n_train)]
            bounds = (lo + 20, hi - 20) if platform == 'PacBio' else (lo, hi)
            n_test = 1000 if platform == 'PacBio' else (200 if facet == 6 else 800)
            test = [int(v) for v in np.logspace(*np.log10(bounds), num=n_test) if int(v) not in train]
            x, y = table['#Protein'].to_numpy(), table['%RelativeExpression'].to_numpy()
            model = LinearRegression(n_jobs=1).fit(np.log10(x[train]).reshape(-1, 1), np.log10(y[train]))
            prediction = 10 ** model.predict(np.log10(x[test]).reshape(-1, 1))
            fits.append(pd.DataFrame({'Platform': label, 'Gene': gene,
                'Rank': 10 ** np.log10(x[test]), 'FittedRelativeExpression': prediction}))
            equations.append({'Platform': label, 'Gene': gene,
                'Slope': model.coef_[0], 'Intercept': model.intercept_,
                'Equation': f'y = x^{model.coef_[0]:.2f} * 10^{model.intercept_:.2f}'})
            for rank in [10, 100, 1000]:
                value = table['%CummulativeRelativeExpression'].iloc[rank - 1]
                markers.append({'Platform': label, 'Gene': gene, 'RankLabel': rank,
                    'CumulativeExpression': value, 'Label': f'({rank}, {value:.1f})',
                    'SourceArrayIndex': rank - 1, 'Protein': table.Protein.iloc[rank - 1],
                    'Solution': selection.Solution})
            common = table.iloc[:math.ceil(.01 * len(table))].assign(Abundancy='1% most expressed')
            rare = table.iloc[-math.ceil(.30 * len(table)):].assign(Abundancy='30% least expressed')
            chosen = pd.concat([common, rare]).copy()
            chosen['Facet'] = facet
            chosen['PlotObservation'] = np.arange(len(chosen))
            observations.append(chosen)
    expression = pd.concat(assignments, ignore_index=True)
    tables = {
        'figure4a_expression': expression,
        'figure_s8_s9_expression': expression.rename(columns={'#Protein': 'Rank',
            '%RelativeExpression': 'RelativeExpression', '%CummulativeRelativeExpression': 'CumulativeExpression'}),
        'figure_s8_fits': pd.concat(fits, ignore_index=True),
        'figure_s8_equations': pd.DataFrame(equations),
        'figure_s9_markers': pd.DataFrame(markers),
        'figure_s17_observations': pd.concat(observations, ignore_index=True)}
    for name, table in tables.items():
        table.to_csv(output_dir / f'{name}.tsv', sep='\t', index=False, na_rep='NA', mode='x')
    return tables


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--assignments', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    expression_tables(args.assignments, args.output)
