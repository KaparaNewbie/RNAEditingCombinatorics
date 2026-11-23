"""Small final-table exports from native notebook products, without saved displays."""
from pathlib import Path
import numpy as np
import pandas as pd
from manuscript_figures.styles import S3_COLORS


def figure6_tables(notebook_dir, output_dir):
    """One observation is an eligible gene in one of the four native cohorts."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=False)
    definitions = [
        ('WholeTranscriptome_TMR50', 'TMR50.Octopus.Pooled'),
        ('WholeTranscriptome_TMR1000', 'TMR1000.Octopus.Pooled'),
        ('SingleCell_Neuronal_FixedTMR50', 'Neuronal.FixedTMR50.Octopus.SC'),
        ('SingleCell_NonNeuronal_FixedTMR50', 'NonNeuronal.FixedTMR50.Octopus.SC')]
    columns = ['Chrom', 'Transcript', 'NumOfReads', 'NumOfProteins', 'Fraction',
        'FractionRepetition', 'Algorithm', 'AlgorithmRepetition', 'MappedReads']
    cohorts, curves, markers = [], [], []
    for index, (cohort, suffix) in enumerate(definitions):
        frame = pd.read_table(Path(notebook_dir) / f'MaxDistinctProtsForFig6.{suffix}.csv', usecols=columns)
        assert frame.Chrom.is_unique and frame.Chrom.notna().all()
        assert frame.NumOfProteins.ge(1).all()
        cohorts.append(frame[columns].assign(Cohort=cohort, SourceRow=np.arange(len(frame))))
        curve = frame[['NumOfProteins']].sort_values('NumOfProteins').reset_index(drop=True)
        curve['CummulativeTranscripts'] = 100 * (curve.index + 1) / len(curve)
        curve['CummulativeTranscripts'] = curve.CummulativeTranscripts.iloc[::-1].values
        if index < 2:
            x, y = curve.NumOfProteins, curve.CummulativeTranscripts
            x_mean = x.mean()
            closest = x.iloc[(x - x_mean).abs().argsort()[:1]]
            k = closest.index.values[0]
            if x_mean == closest.values[0]:
                y_mean = y.iloc[k]
            else:
                i, j = (k - 1, k + 1) if x_mean < closest.values[0] else (k, k + 2)
                y_mean = np.interp(x_mean, x.iloc[i:j], y.iloc[i:j])
            markers.append({'Cohort': cohort, 'X': x_mean, 'Y': y_mean})
        curve = curve.drop_duplicates('NumOfProteins').reset_index(drop=True)
        curves.append(curve.assign(Cohort=cohort, DenominatorGenes=len(frame)))
    tables = {'figure6_cohorts': pd.concat(cohorts, ignore_index=True),
        'figure6_distributions': pd.concat(curves, ignore_index=True),
        'figure6_mean_markers': pd.DataFrame(markers)}
    for name, table in tables.items():
        table.to_csv(output_dir / f'{name}.tsv', sep='\t', index=False, na_rep='NA', mode='x')
    return tables


def figure_s3_tables(native_dir, output_dir):
    """Join native read distributions and their original mean/population SD."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=False)
    counts, statistics = [], []
    for platform in ['PacBio', 'Illumina']:
        prefix = Path(native_dir) / f'figure_s3_{platform.lower()}_native'
        count = pd.read_table(str(prefix) + '_counts.tsv', float_precision='round_trip')
        stats = pd.read_table(str(prefix) + '_statistics.tsv', float_precision='round_trip')
        identity = count[['Platform', 'Gene', 'GeneOrder']].drop_duplicates()
        assert identity.Gene.is_unique
        stats = identity.merge(stats, on='Gene', validate='one_to_many')
        assert not stats.duplicated(['Gene', 'Estimate']).any()
        totals = count.groupby(['Gene', 'Estimate'], sort=False).Reads.sum()
        np.testing.assert_array_equal(stats.DenominatorReads, totals.reindex(pd.MultiIndex.from_frame(stats[['Gene', 'Estimate']])))
        stats['Annotation'] = [f'{m:.0f}±{s:.0f}' for m, s in zip(stats.Mean, stats.Std)]
        stats['Color'] = [S3_COLORS[(platform, r.Gene, r.Estimate)] for r in stats.itertuples()]
        stats['YLimit'] = 7995 * 1.2 if platform == 'PacBio' else 284793 * 1.3
        stats['XLimit'] = 94 if platform == 'PacBio' else 43
        counts.append(count)
        statistics.append(stats)
    tables = {'figure_s3_counts': pd.concat(counts, ignore_index=True),
        'figure_s3_statistics': pd.concat(statistics, ignore_index=True)}
    for name, table in tables.items():
        table.to_csv(output_dir / f'{name}.tsv', sep='\t', index=False, mode='x')
    return tables


def assemble_fisher_tables(producer_dir, output_dir):
    """Combine native site counts/results; site order remains the producer order."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=False)
    tables = {}
    for kind in ['sites', 'summary']:
        parts = [pd.read_table(Path(producer_dir) / f's18_fisher_{kind}.{p}.tsv',
            dtype={'EditingSite': str}, float_precision='round_trip') for p in ['PacBio', 'Illumina']]
        tables[kind] = pd.concat(parts, ignore_index=True)
    assert len(tables['sites']) == 1186 and len(tables['summary']) == 21
    assert not tables['sites'].duplicated(['Platform', 'Gene', 'EditingSite']).any()
    tables['sites'].to_csv(output_dir / 's18_fisher_site_results.tsv', sep='\t', index=False, na_rep='NA', mode='x')
    tables['summary'].to_csv(output_dir / 's18_fisher_gene_summary.tsv', sep='\t', index=False, mode='x')
    plot = tables['sites'].copy()
    plot.Platform = plot.Platform.map({'PacBio': 'Long-reads', 'Illumina': 'Short-reads'})
    # Original S18 plotting order is long reads, then short reads, native gene/site order.
    plot.to_csv(output_dir / 'figure_s18_fisher.tsv', sep='\t', index=False, na_rep='NA', mode='x')
    return tables


def apply_by_matrices(matrix_table, producer_dir):
    """Replace only Pearson values using the native BY masks; MI remains untouched."""
    table = matrix_table.copy()
    original_mi = table.loc[table.Measure.eq('MI')].copy()
    for (platform, gene), group in table.loc[table.Measure.eq('Pearson')].groupby(['Platform', 'Gene'], sort=False):
        matrix = pd.read_table(Path(producer_dir) / f'{gene}MaskedBYCorr.{platform}.tsv', float_precision='round_trip').to_numpy()
        order = pd.read_table(Path(producer_dir) / f'{gene}.BYSiteOrder.{platform}.tsv').EditingSite.to_numpy()
        np.testing.assert_array_equal(group.loc[group.ColumnIndex.eq(0), 'RowSite'], order)
        np.testing.assert_array_equal(group.loc[group.RowIndex.eq(0), 'ColumnSite'], order)
        values = matrix[group.RowIndex, group.ColumnIndex]
        table.loc[group.index, 'Value'] = values
        table.loc[group.index, 'Masked'] = np.isnan(values)
    pd.testing.assert_frame_equal(table.loc[table.Measure.eq('MI')], original_mi, check_exact=True)
    return table


def export_native_s3(project_root, output_dir):
    """Original protein-pattern order expanded by original read multiplicity."""
    from manuscript_inputs import selections
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=False)
    selected = selections().loc[lambda t: t.Metric.eq('ordinary')]
    for platform, genes in selected.groupby('Platform', sort=False):
        counts, statistics = [], []
        for order, selection in enumerate(genes.itertuples(index=False)):
            path = Path(project_root) / selection.UniqueProteinsFile.replace('.unique_proteins.', '.proteins.')
            proteins = pd.read_table(path, usecols=['MinNonSyns', 'MaxNonSyns', 'Reads', 'NumOfReads'])
            np.testing.assert_array_equal(proteins.Reads.str.count(',') + 1, proteins.NumOfReads)
            assert proteins.NumOfReads.gt(0).all()
            for estimate in ['Min', 'Max']:
                values = np.repeat(proteins[estimate + 'NonSyns'].to_numpy(), proteins.NumOfReads.to_numpy())
                assert np.issubdtype(values.dtype, np.integer) and (values >= 0).all()
                counts.extend({'Platform': platform, 'Gene': selection.Gene, 'GeneOrder': order,
                    'Estimate': estimate, 'Substitutions': value, 'Reads': int(count)}
                    for value, count in enumerate(np.bincount(values)))
                statistics.append({'Gene': selection.Gene, 'Estimate': estimate,
                    'DenominatorReads': len(values), 'Mean': np.mean(values), 'Std': np.std(values)})
        pd.DataFrame(counts).to_csv(output_dir / f'figure_s3_{platform.lower()}_native_counts.tsv', sep='\t', index=False, mode='x')
        pd.DataFrame(statistics).to_csv(output_dir / f'figure_s3_{platform.lower()}_native_statistics.tsv', sep='\t', index=False, mode='x')


def figure_s13_tables(notebook_dir, project_root, output_dir):
    """Top nine pooled genes: original selected-protein and pileup-read membership."""
    from collections import defaultdict
    from manuscript_figures.styles import S13_SAMPLE_ORDER, S13_COLORS
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=False)
    base = Path(project_root) / 'O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.PooledSamples'
    cohort = pd.read_table(Path(notebook_dir) / 'MaxDistinctProtsForFig6.TMR50.Octopus.Pooled.csv', dtype={'Proteins': str})
    cohort = cohort.sort_values('NumOfProteins', ascending=False).iloc[:9]
    samples = pd.read_csv(Path(project_root) / 'O.vulgaris/Data/PRJNA791920/IsoSeqPolished/samples.csv')
    samples.Tissue = samples.Tissue.str.capitalize()
    assert samples.Sample.is_unique and samples.Tissue.is_unique
    parts, memberships, coverage = [], [], []
    for facet, gene in enumerate(cohort.itertuples()):
        proteins = pd.read_table(base / 'ProteinsFiles' / f'{gene.Chrom}.unique_proteins.csv.gz',
            usecols=['Transcript', 'Protein', 'Samples'], dtype=str)
        selected = gene.Proteins.split(',')
        assert len(selected) == int(gene.NumOfProteins) == len(set(selected))
        proteins = proteins.set_index('Protein').loc[selected].reset_index()
        assert proteins.Transcript.eq(gene.Transcript).all()
        membership = proteins.assign(Sample=proteins.Samples.str.split(',').map(set)).explode('Sample')[['Protein', 'Sample']]
        assert not membership.duplicated().any()
        counts = membership.Sample.value_counts().rename_axis('Sample').reset_index(name='NumOfProteinsPerGenePerSample')
        counts.insert(0, 'Transcript', gene.Transcript)
        counts.insert(0, 'Chrom', gene.Chrom)
        counts['PooledIsoforms'], counts['Facet'] = int(gene.NumOfProteins), facet
        positions = pd.read_table(base / 'PositionsFiles' / f'{gene.Chrom}.positions.csv.gz',
            usecols=['CDS', 'InProbRegion', 'Samples', 'Reads'])
        eligible = positions.loc[positions.CDS & ~positions.InProbRegion, ['Samples', 'Reads']]
        assert eligible.notna().all().all()
        seen = defaultdict(set)
        for sample_list, read_list in eligible.itertuples(index=False, name=None):
            sample_values, read_values = sample_list.split(','), read_list.split(',')
            assert len(sample_values) == len(read_values)
            for sample, read in zip(sample_values, read_values):
                seen[sample].add(read)
        counts['NumOfReads'] = counts.Sample.map({k: len(v) for k, v in seen.items()})
        counts = counts.merge(samples, on='Sample', how='left', validate='many_to_one')
        assert counts.notna().all().all()
        assert counts.NumOfReads.ge(counts.NumOfProteinsPerGenePerSample).all()
        counts['Color'] = counts.Tissue.map(S13_COLORS)
        counts = counts.set_index('Sample').loc[S13_SAMPLE_ORDER[gene.Chrom]].reset_index()
        counts['ReferenceTraceOrder'] = np.arange(len(counts))
        parts.append(counts[['Chrom', 'Transcript', 'Sample', 'NumOfProteinsPerGenePerSample',
            'PooledIsoforms', 'Facet', 'NumOfReads', 'Tissue', 'Color', 'ReferenceTraceOrder']])
        memberships.append(membership.assign(Chrom=gene.Chrom))
        coverage.extend({'Chrom': gene.Chrom, 'Sample': sample, 'Read': read}
            for sample, reads in seen.items() for read in sorted(reads))
    tables = {'figure_s13_plot': pd.concat(parts, ignore_index=True),
        'figure_s13_protein_membership': pd.concat(memberships, ignore_index=True).sort_values(['Chrom', 'Protein', 'Sample']).reset_index(drop=True),
        'figure_s13_coverage_membership': pd.DataFrame(coverage)}
    for name, table in tables.items():
        table.to_csv(output_dir / f'{name}.tsv', sep='\t', index=False, na_rep='NA', mode='x')
    return tables


def figure_s16_table(native_dir, destination):
    """Format the two native graph-assessment exports; do not rerun simulations."""
    parts = []
    for panel in ['a', 'b']:
        table = pd.read_table(Path(native_dir) / f'source_s16{panel}.tsv', float_precision='round_trip')
        assert len(table) == 90 and table.Panel.eq(panel).all()
        assert table[['Mean', 'SD']].notna().all().all() and table.SD.ge(0).all()
        table['Row'] = table.Trace // (6 if panel == 'a' else 9)
        table['Column'] = (table.Trace % 6) // 2 if panel == 'a' else 0
        table['Color'] = table.DataType.map({'Complete': 'black', 'Errored+PartiallyUnknown': 'grey'}) if panel == 'a' else table.Gene.map({'CBPC1': 'red', 'MY18A': 'blue', 'VPS8': 'green'})
        table['Dash'] = 'solid' if panel == 'a' else table.UnknownProbability.map({.09: 'solid', .13: 'dash', .18: 'dot'})
        table['Opacity'] = .7
        parts.append(table[['Panel', 'Metric', 'Chrom', 'Gene', 'UnknownProbability', 'DataType',
            'NumOfReads', 'Fraction', 'Mean', 'SD', 'Trace', 'Row', 'Column', 'Color', 'Dash', 'Opacity']])
    table = pd.concat(parts, ignore_index=True)
    table.to_csv(destination, sep='\t', index=False, na_rep='NA', mode='x')
    return table


def figure_s1_tables(notebook_dir, output_dir):
    """Native read ECDFs/deletion bins and gene-noise histogram; no KDE refit."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=False)
    raw = pd.read_table(Path(notebook_dir) / 'RawReadsStats.Squid.PacBio.csv')
    curves, bins = [], []
    edges = .989975 + np.arange(202) * .00005
    for gene in ['GRIA2', 'PCLO']:
        data = raw.loc[raw.Gene.eq(gene)]
        rank = np.arange(1, len(data) + 1)
        for panel, x, y in [('a', np.sort(data.ReadLen), rank / len(data) * 100),
                            ('b', np.sort(data.ReadQuality)[::-1], 100 * rank / len(data))]:
            curves.append(pd.DataFrame({'Panel': panel, 'Gene': gene, 'Rank': rank,
                'X': x, 'Percent': y, 'DenominatorReads': len(data)}))
        count = np.histogram(data.ReadQuality, edges)[0]
        total = np.histogram(data.ReadQuality, edges, weights=data.Deletions)[0]
        assert count.sum() == len(data) and (count > 0).all()
        bins.append(pd.DataFrame({'Gene': gene, 'Bin': np.arange(201), 'Left': edges[:-1],
            'Right': edges[1:], 'ReadCount': count, 'TotalDeletions': total, 'MeanDeletions': total / count}))
    octopus = pd.read_table(Path(notebook_dir) / 'Noise.Octopus.WholeTranscriptome.Pooled.tsv')
    edges = -.05 + np.arange(51) * .1
    counts = np.histogram(octopus['%Noise'], edges)[0]
    assert octopus.Chrom.is_unique and counts.sum() == len(octopus) == 1711
    tables = {'figure_s1_curves': pd.concat(curves, ignore_index=True),
        'figure_s1_deletions': pd.concat(bins, ignore_index=True),
        'figure_s1_octopus_histogram': pd.DataFrame({'Bin': range(50), 'Left': edges[:-1], 'Right': edges[1:], 'Genes': counts})}
    for name, table in tables.items():
        table.to_csv(output_dir / f'{name}.tsv', sep='\t', index=False, na_rep='NA', mode='x')
    return tables


def figure3_tables(notebook_dir, assignments_dir, output_dir, project_root):
    """Original Figure3 definitions; heatmap ambiguity is the legacy display model."""
    from manuscript_inputs import selections, read_assignment, selected_amino_acids
    from manuscript_figures.styles import GENE_COLORS
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=False)
    choice = selections().loc[lambda t: t.Platform.eq('PacBio') & t.Gene.eq('PCLO') & t.Metric.eq('ordinary')].iloc[0]
    assignments = read_assignment(choice, assignments_dir)
    aa = selected_amino_acids(choice, assignments, project_root)
    # Preserve the approved Figure3 model; S7 uses the separately corrected binary helper.
    matrix = pd.DataFrame({c: np.where(aa[c].eq(c.split('(')[1][0]), 0.,
        np.where(aa[c].str.contains(',', regex=False), np.nan, 1.)) for c in aa})
    matrix = matrix.loc[:, matrix.mean().gt(0)]
    a = pd.concat([assignments[['Gene', 'Protein', '%RelativeExpression']].reset_index(drop=True), matrix.reset_index(drop=True)], axis=1)
    a = a.sort_values('%RelativeExpression', ascending=False).reset_index(drop=True).iloc[:100].copy()
    a.insert(0, 'IsoformRank', range(1, 101))
    assert a.shape == (100, 98) and a.iloc[:, 4:].isna().sum().sum() == 623
    bars, curves = [], []
    for platform in ['PacBio', 'Illumina']:
        label = 'Long-reads' if platform == 'PacBio' else 'Short-reads'
        maxima = pd.read_table(Path(notebook_dir) / f'MaxDistinctProteinsF1.{platform}.tsv')
        maxima.Gene = maxima.Gene.str.split('_').str[0].replace({'GRIA': 'GRIA2'})
        maxima = maxima.sort_values('Gene')[['Platform', 'Gene', 'Max']]
        maxima['Color'] = maxima.Gene.map(GENE_COLORS)
        bars.append(maxima)
        native = pd.read_table(Path(notebook_dir) / f'DistinctProteins.{platform}.tsv', usecols=['Platform', 'Gene', 'NumOfReads', 'Fraction', 'FractionRepetition', 'Algorithm', 'AlgorithmRepetition', 'NumOfProteins'])
        native.Gene = native.Gene.str.split('_').str[0].replace({'GRIA': 'GRIA2'})
        if platform == 'Illumina': native = native.sort_values('Gene')
        for gene in maxima.Gene:
            selected = native.loc[native.Gene.eq(gene)].sort_values(['Fraction', 'NumOfProteins'], ascending=False).drop_duplicates('Fraction', ignore_index=True)
            curves.append(selected[['Platform', 'Gene', 'NumOfReads', 'Fraction', 'FractionRepetition', 'Algorithm', 'AlgorithmRepetition', 'NumOfProteins']])
    dissimilar = pd.read_table(Path(notebook_dir) / 'DistinctDissimilarProteins.PacBio.tsv', usecols=['Gene', 'NumOfReads', 'Fraction', 'FractionRepetition', 'Algorithm', 'AlgorithmRepetition', 'NumOfProteins', 'Dissimilarity', 'CutoffScore'])
    dissimilar.Gene = dissimilar.Gene.replace({'GRIA': 'GRIA2'})
    different = []
    for metric in ['Miyata', 'Grantham 100']:
        for gene in ['GRIA2', 'PCLO']:
            table = dissimilar.loc[dissimilar.Dissimilarity.eq(metric) & dissimilar.Gene.eq(gene)]
            table = table.sort_values(['Fraction', 'NumOfProteins'], ascending=False).drop_duplicates('Fraction', ignore_index=True)
            different.append(table[['Gene', 'NumOfReads', 'Fraction', 'FractionRepetition', 'Algorithm', 'AlgorithmRepetition', 'NumOfProteins', 'Dissimilarity', 'CutoffScore']])
    tables = {'a': a, 'b': pd.concat(bars, ignore_index=True), 'c': pd.concat(curves, ignore_index=True), 'd': pd.concat(different, ignore_index=True)}
    for panel, table in tables.items():
        table.to_csv(output_dir / f'squid_original_figure3{panel}_plot.tsv', sep='\t', index=False, na_rep='', mode='x')
    return tables


def figure_s12_tables(notebook_dir, output_dir):
    """Native BY-significant site pairs; original sorting and Hazen arithmetic."""
    from plotly.colors import sample_colorscale
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=False)
    pairs, boxes, outliers = [], [], []
    for platform, tag in [('Long-reads', 'PB1'), ('Short-reads', 'IL1')]:
        native = pd.read_table(Path(notebook_dir) / f'Corrs.{tag}.tsv', float_precision='round_trip')
        native = native.loc[native.fdr_by_rejection].reset_index(drop=True)
        assert native.r.notna().all() and not native.duplicated(['Gene', 'Site1', 'Site2']).any()
        np.testing.assert_array_equal(native.AbsDistance100BPsBin, ((native.Site2-native.Site1).abs()//100)*100+100)
        native['SitePairsInAbsDistance100BPsBin'] = native.groupby(['Gene', 'AbsDistance100BPsBin']).r.transform('size')
        native['bin_start'] = (native.SitePairsInAbsDistance100BPsBin // 100) * 100
        native['bin_label'] = native.bin_start.astype(str) + ' - ' + (native.bin_start + 100).astype(str)
        native = native.sort_values(['bin_start', 'Gene', 'AbsDistance100BPsBin'] if platform == 'Long-reads' else ['bin_start'])
        labels = native.bin_label.drop_duplicates().tolist()
        genes = native.Gene.drop_duplicates().tolist() if platform == 'Long-reads' else sorted(native.Gene.unique())
        colors = sample_colorscale('Viridis', [i/max(1, len(labels)-1) for i in range(len(labels))])
        trace = 0
        for label, color in zip(labels, colors):
            for gene in genes:
                data = native.loc[native.Gene.eq(gene) & native.bin_label.eq(label)].copy()
                if data.empty: continue
                data.insert(0, 'Platform', platform)
                data['Trace'], data['Color'] = trace, color
                trace += 1
                pairs.append(data)
                for distance, group in data.groupby('AbsDistance100BPsBin', sort=True):
                    values = np.sort(group.r.to_numpy())
                    def quartile(p):
                        index = len(values)*p-.5
                        lo, hi = int(np.floor(index)), int(np.ceil(index))
                        return values[max(0, lo)] + (values[min(len(values)-1, hi)] - values[max(0, lo)])*(index-lo)
                    q1, median, q3 = [quartile(p) for p in [.25, .5, .75]]
                    low, high = values[values >= q1-1.5*(q3-q1)][0], values[values <= q3+1.5*(q3-q1)][-1]
                    boxes.append({'Platform': platform, 'Gene': gene, 'DistanceBin': int(distance), 'Pairs': len(values),
                        'PairCountLabel': label, 'Color': color, 'q1': q1, 'median': median, 'q3': q3, 'whislo': low, 'whishi': high})
                    for row in group.loc[(group.r < low) | (group.r > high)].itertuples():
                        outliers.append({'Platform': platform, 'Gene': gene, 'DistanceBin': int(distance),
                            'Site1': row.Site1, 'Site2': row.Site2, 'r': row.r, 'Color': color})
    tables = {'pairs': pd.concat(pairs, ignore_index=True), 'boxes': pd.DataFrame(boxes), 'outliers': pd.DataFrame(outliers)}
    for kind, table in tables.items():
        table.to_csv(output_dir / f'figure_s12_{kind}.tsv', sep='\t', index=False, mode='x')
    return tables


def figure_s10_tables(notebook_dir, assignments_dir, output_dir):
    """PB1 fraction comparison and separate PB3 recognizable-barcode/UMI branches."""
    from sklearn.model_selection import train_test_split
    from sklearn.linear_model import LinearRegression
    from manuscript_inputs import selections, read_assignment
    notebook_dir, output_dir = Path(notebook_dir), Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=False)
    columns = ['Platform', 'Gene', 'Protein', '#Protein0.1', '#Protein1.0', '%RelativeExpression0.1', '%RelativeExpression1.0']
    paired = pd.read_table(notebook_dir / 'AssignedExpression.Fractions1Vs01.PacBio.tsv', usecols=columns, dtype={'Protein': str})
    paired = paired.dropna(subset=columns[-2:]).reset_index(drop=True)
    fits = []
    for gene in ['GRIA2', 'PCLO']:
        pair = paired.loc[paired.Gene.eq(gene)]
        x_train, x_test, y_train, _ = train_test_split(pair['%RelativeExpression0.1'], pair['%RelativeExpression1.0'], test_size=.2, random_state=1892)
        model = LinearRegression(n_jobs=1).fit(x_train.to_numpy().reshape(-1, 1), y_train)
        fits.append(pd.DataFrame({'Gene': gene, 'X': x_test, 'Y': model.predict(x_test.to_numpy().reshape(-1, 1))}))
    curves = []
    definitions = [('b', 'Subsampled data (10%)', 'AssignedExpression.Fraction01.PacBio.tsv'),
        ('b', 'Full data', None), ('c', 'Full data', 'AssignedExpression.PacBio3.ReadsWithRecognizableBarcodes.tsv'),
        ('c', 'Deduplicated data (UMIs)', 'AssignedExpression.PacBio3.DeduppedReads.tsv')]
    for panel, label, filename in definitions:
        if filename is None:
            choices = selections().loc[lambda t: t.Platform.eq('PacBio') & t.Metric.eq('ordinary')]
            native = pd.concat([read_assignment(s, assignments_dir) for s in choices.itertuples(index=False)], ignore_index=True)
        else:
            header = pd.read_table(notebook_dir / filename, nrows=0).columns
            wanted = ['Gene', '#Protein', '%CummulativeRelativeExpression'] + (['Protein'] if 'Protein' in header else [])
            native = pd.read_table(notebook_dir / filename, usecols=wanted, dtype={'Protein': str})
        native.Gene = native.Gene.replace({'GRIA': 'GRIA2'})
        for gene in (['GRIA2', 'PCLO'] if panel == 'b' else ['GRIA2', 'ADAR1', 'IQEC1']):
            columns = [c for c in ['Gene', '#Protein', '%CummulativeRelativeExpression', 'Protein'] if c in native]
            data = native.loc[native.Gene.eq(gene), columns].copy()
            data = data.rename(columns={'#Protein': 'Rank', '%CummulativeRelativeExpression': 'CumulativeExpression'})
            data.insert(0, 'Panel', panel)
            data.insert(2, 'Group', label)
            assert data.Rank.tolist() == list(range(1, len(data)+1))
            curves.append(data)
    tables = {'paired': paired, 'fits': pd.concat(fits, ignore_index=True), 'cumulative': pd.concat(curves, ignore_index=True)[['Panel', 'Gene', 'Group', 'Protein', 'Rank', 'CumulativeExpression']]}
    for kind, table in tables.items():
        table.to_csv(output_dir / f'figure_s10_{kind}.tsv', sep='\t', index=False, na_rep='NA', mode='x')
    return tables


def figure_s6_tables(notebook_dir, assignments_dir, output_dir, project_root):
    """Dispersion, PCLO top100 raw AA counts, and the distinct 80K-pair cohort."""
    from plotly.colors import qualitative
    from scipy.stats import pearsonr
    from manuscript_inputs import selections, read_assignment, selected_amino_acids, verified_file
    from manuscript_figures.styles import GENE_COLORS
    notebook_dir, output_dir = Path(notebook_dir), Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=False)
    bars = []
    for platform in ['PacBio', 'Illumina']:
        data = pd.read_table(notebook_dir / f'Dispersion.{platform}.tsv')
        data.Gene = data.Gene.str.split('_').str[0].replace({'GRIA': 'GRIA2'})
        order = ['GRIA2', 'PCLO'] if platform == 'PacBio' else ['PCLO'] + sorted(set(data.Gene)-{'PCLO'})
        for gene in order:
            row = data.loc[data.Gene.eq(gene)].iloc[0]
            bars.append({'Platform': platform, 'Gene': gene, 'MinimumSolutions': int(row['min']),
                'MaximumSolutions': int(row['max']), 'PercentDispersion': row['%SolutionsDispersion'], 'Color': GENE_COLORS[gene]})
    octopus = pd.read_table(notebook_dir / 'Dispersion.Octopus.WholeTranscriptome.Pooled.tsv')
    assert octopus.Chrom.is_unique and len(octopus) == 1262
    edges = -.25 + np.arange(102)*.5
    hist = np.histogram(octopus['%SolutionsDispersion'], edges)[0]
    assert hist.sum() == len(octopus)
    histogram = pd.DataFrame({'Bin': range(101), 'Left': edges[:-1], 'Right': edges[1:], 'Genes': hist})

    choice = selections().loc[lambda t: t.Platform.eq('PacBio') & t.Gene.eq('PCLO') & t.Metric.eq('ordinary')].iloc[0]
    assignment = read_assignment(choice, assignments_dir)
    aa = selected_amino_acids(choice, assignment, project_root).iloc[:100]
    assert aa.shape == (100, 131)
    counts = aa.apply(lambda x: x.value_counts()).reset_index().rename(columns={'index': 'NewAA'}).fillna(0)
    counts = counts.melt(id_vars='NewAA', var_name='Position', value_name='Count')
    counts['OriginalAA'] = counts.Position.str.split('(').str[1].str[0]
    counts = counts.loc[[ref not in state.split(',') for ref, state in zip(counts.OriginalAA, counts.NewAA)]]
    present = set(counts.loc[counts.Count.gt(0), 'NewAA'])
    categories = counts.NewAA.unique()
    positive_colors = iter(qualitative.Light24)
    colors = {aa: next(positive_colors) if aa in present else 'white' for aa in categories}
    codon_order = list(dict.fromkeys(c for aa in categories for c in counts.loc[counts.NewAA.eq(aa), 'Position']))
    recoding = []
    for trace, state in enumerate(categories):
        for row in counts.loc[counts.NewAA.eq(state)].itertuples():
            recoding.append({'Trace': trace, 'Codon': row.Position, 'CodonOrder': codon_order.index(row.Position),
                'NewAA': state, 'Count': int(row.Count), 'Color': colors[state], 'ShowLegend': state in present})

    definitions = pd.read_table(Path(__file__).parent / 'manuscript_input_data/s6_depth_inputs.tsv')
    runs, maxima = [], []
    for item in definitions.itertuples(index=False):
        path = verified_file(project_root, item.Solutions, item.SolutionsSHA256)
        run = pd.read_table(path, usecols=['Fraction', 'FractionRepetition', 'Algorithm', 'AlgorithmRepetition', 'NumUniqueSamples'])
        run = run.loc[run.Fraction.eq(1)].rename(columns={'NumUniqueSamples': 'NumOfProteins'})
        read_count = pd.read_table(Path(project_root) / item.UniqueReads, usecols=['NumOfReads']).NumOfReads.sum()
        codons = pd.read_table(Path(project_root) / item.UniqueProteins, nrows=0).columns[14:]
        run.insert(0, 'Platform', item.Platform)
        run.insert(1, 'Gene', item.Gene)
        run.insert(2, 'NumOfReads', (run.Fraction * read_count).astype(int))
        run['EditableAAs'] = len(codons)
        assert len(run) == 16
        runs.append(run)
        maxima.append({'Platform': item.Platform, 'Gene': item.Gene, 'EditableCodons': len(codons),
            'MaximumIsoforms': int(run.NumOfProteins.max()), 'Color': GENE_COLORS[item.Gene],
            'Marker': 'star' if item.Platform == 'Long-reads' else 'circle'})
    runs = pd.concat(runs, ignore_index=True).sort_values(['NumOfProteins', 'EditableAAs'], ascending=False).reset_index(drop=True)
    maxima = pd.DataFrame(maxima)
    correlation = pearsonr(maxima.EditableCodons, maxima.MaximumIsoforms)
    statistics = pd.DataFrame([{'N': len(maxima), 'PearsonR': correlation.statistic, 'PValue': correlation.pvalue,
        'Annotation': f'Pearson r = {correlation.statistic:.2f}\np = {correlation.pvalue:.2g}'}])
    tables = {'figure_s6a_squid': pd.DataFrame(bars), 'figure_s6a_histogram': histogram,
        'figure_s6b_recoding': pd.DataFrame(recoding), 'figure_s6c_solutions': runs,
        'figure_s6c_maxima': maxima, 'figure_s6c_statistics': statistics}
    for name, table in tables.items():
        table.to_csv(output_dir / f'{name}.tsv', sep='\t', index=False, mode='x')
    return tables


def dependency_matrix_tables(notebook_dir, by_dir, output_root):
    """Format native MI matrices and newly generated BY Pearson matrices.

    Site identities come from the native unordered-pair tables. This performs no
    new correlation inference and preserves each native matrix's flattening order.
    """
    definitions = [('Figure5', 'figure5_pclo_matrices.BY', [('PacBio', 'PCLO', 'PB1'), ('Illumina', 'PCLO', 'IL1')]),
                   ('FigureS11', 'figure_s11_matrices.BY', [('PacBio', 'GRIA2', 'PB1')])]
    for folder, name, cohorts in definitions:
        directory = Path(output_root) / folder
        directory.mkdir(parents=True, exist_ok=True)
        parts = []
        for platform, gene, tag in cohorts:
            pairs = pd.read_table(Path(notebook_dir) / f'Corrs.{tag}.tsv')
            pairs = pairs.loc[pairs.Gene.eq(gene)]
            sites = np.array(sorted(set(pairs.Site1) | set(pairs.Site2)))
            pearson = pd.read_table(Path(by_dir) / f'{gene}MaskedBYCorr.{platform}.tsv', float_precision='round_trip').to_numpy()
            mi = pd.read_table(Path(notebook_dir) / f'{gene}MaskedMI.{platform}.tsv').to_numpy()
            for measure, values in [('Pearson', pearson), ('MI', mi)]:
                assert values.shape == (len(sites), len(sites))
                rows, columns = np.indices(values.shape)
                table = pd.DataFrame({'Platform': platform, 'Gene': gene, 'Measure': measure,
                    'RowIndex': rows.ravel(), 'ColumnIndex': columns.ravel(),
                    'RowSite': sites[rows.ravel()], 'ColumnSite': sites[columns.ravel()],
                    'Masked': np.isnan(values).ravel(), 'Value': values.ravel()})
                if folder == 'FigureS11':
                    table = table[['Gene', 'Platform', *table.columns[2:]]]
                parts.append(table)
        pd.concat(parts, ignore_index=True).to_csv(directory / f'{name}.tsv', sep='\t', index=False, na_rep='NA', mode='x')


def figure2_table(project_root, destination):
    """PCLO's 100 most abundant unique RNA patterns; preserve cutoff tie order."""
    from manuscript_inputs import selections
    selection = selections().loc[lambda t: t.Platform.eq('PacBio') & t.Gene.eq('PCLO') & t.Metric.eq('ordinary')].iloc[0]
    source = pd.read_table(Path(project_root) / selection.UniqueReadsFile, dtype={'Transcript': str})
    source = source.rename(columns={'Transcript': 'UniqueRead'})
    assert source.UniqueRead.is_unique
    sites = source.columns[8:].tolist()
    assert len(sites) == 144 and all(str(c).isdigit() for c in sites)
    table = source.sort_values('NumOfReads', ascending=False, kind='stable').iloc[:100]
    table = table[['Gene', 'UniqueRead', 'NumOfReads', *sites]].copy()
    table.insert(0, 'UniqueReadRank', range(1, 101))
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(destination, sep='\t', index=False, mode='x')
    return table
