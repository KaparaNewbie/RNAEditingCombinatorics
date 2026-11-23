"""S2 native site cohorts and plotting tables, preserving the distinct motif cohort.

Squid uses its original final stored Edited flag; octopus uses EditedFinal.
The octopus Venn and motif populations intentionally differ (1262/2701 genes).
Native site projection is a separate, potentially expensive step; rendering never
reads original per-position read lists or reruns site calling.
"""
from pathlib import Path
import argparse
import csv
import gzip
import hashlib
import json
import re
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
from Bio import SeqIO
from scipy.stats import pearsonr
from manuscript_inputs import PROJECT, selections


def digest(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def read_sites(path, octopus=False):
    columns = ['Chrom', 'Position', 'Edited', 'CDS', 'KnownEditing', 'InProbRegion', 'RefBase', 'EditingFrequency']
    if octopus:
        columns += ['EditedFinal']
    opener = gzip.open if path.suffix == '.gz' else open
    # Project scalar fields before pandas; a row may have millions of read IDs.
    import io
    stream = io.StringIO()
    with opener(path, 'rt') as source:
        header = source.readline().rstrip('\n').split('\t')
        indices = [header.index(c) for c in columns]
        writer = csv.writer(stream, delimiter='\t', lineterminator='\n')
        writer.writerow(columns)
        for line in source:
            fields = line.rstrip('\n').split('\t')
            assert len(fields) == len(header)
            writer.writerow([fields[i] for i in indices])
    stream.seek(0)
    table = pd.read_table(stream)
    assert not table.duplicated(['Chrom', 'Position']).any()
    return table


def export_sites(branch, output_dir, project_root=PROJECT):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=False)
    parts = []
    if branch == 'squid':
        cohort = selections().loc[lambda t: t.Metric.eq('ordinary')]
        for row in cohort.itertuples(index=False):
            path = Path(project_root) / row.UniqueProteinsFile.replace('.unique_proteins.', '.positions.')
            table = read_sites(path)
            table.insert(0, 'Gene', row.Gene)
            table.insert(0, 'Platform', row.Platform)
            parts.append(table)
        name = 's2_squid_sites'
    else:
        motif = branch == 'octopus_motif'
        suffix = 'BQ30.AHL.BHAfterNoise.3' if motif else 'PooledSamples'
        base = Path(project_root) / f'O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.{suffix}'
        pattern = '*.DistinctUniqueProteins.ExpressionLevels.csv' if motif else '*.DistinctUniqueProteins.ExpressionLevels.*csv'
        expressions = sorted((base / 'ExpressionLevels').glob(pattern))
        chroms = [p.name.split('.')[0] for p in expressions]
        assert len(chroms) == len(set(chroms)) == (2701 if motif else 1262)
        expected_cohorts = {'octopus': '913ad1a9e4b1396104d1c92a6e85561a1c8b723b07cd4030142b385ae60ae59d', 'octopus_motif': '166b2231506ba1b0a64f749bc3182d8e41fb4fee8959ba4f66eac6bf68fbef3c'}
        assert hashlib.sha256(('\n'.join(chroms) + '\n').encode()).hexdigest() == expected_cohorts[branch]
        orfs = pd.read_table(Path(project_root) / 'O.vulgaris/Annotations/orfs_oct.bed', names=['Chrom', 'Start', 'End', 'Gene', 'Score', 'Strand'])
        cohort = pd.DataFrame({'Chrom': chroms}).merge(orfs, on='Chrom', validate='one_to_one')
        assert len(cohort) == len(chroms)
        for row in cohort.itertuples(index=False):
            table = read_sites(base / 'PositionsFiles' / f'{row.Chrom}.positions.csv.gz', octopus=True)
            keep = table.EditedFinal if motif else table.EditedFinal | (table.CDS & table.KnownEditing)
            table = table.loc[keep].copy()
            table['Strand'], table['ChromEnd'], table['Gene'] = row.Strand, row.End, row.Gene
            parts.append(table)
        name = 's2_octopus_original_motif_sites' if motif else 's2_octopus_sites'
    table = pd.concat(parts, ignore_index=True)
    table.to_csv(output_dir / f'{name}.tsv', sep='\t', index=False, na_rep='NA', mode='x')
    return table


def produce_s2(site_dir, output_dir, distribution_svg, project_root=PROJECT):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=False)
    PROJECT = Path(project_root)
    tables = {}
    def save(name, table, settings):
        tables[name] = table
        if name == 's2_pooled':
            tables['s2_statistics'] = pd.DataFrame([{k: settings[k] for k in ['pearson_n', 'pearson_r', 'pearson_p']}])
    squid = pd.read_table(Path(site_dir) / 's2_squid_sites.tsv', float_precision='round_trip')
    sites = squid
    octopus = pd.read_table(Path(site_dir) / 's2_octopus_sites.tsv', float_precision='round_trip')
    octopus_motif = pd.read_table(Path(site_dir) / 's2_octopus_original_motif_sites.tsv', float_precision='round_trip')
    known_path = PROJECT/'D.pealeii/Annotations/D.pea.EditingSites.csv'
    known = pd.read_csv(known_path)
    known.columns = ['Chrom', 'SwissProt', 'Position', 'OriginalAA', 'NewAA', 'RefBase',
                     'Editing', 'EditingFrequency', 'DNA']
    known.Position -= 1
    reference = known[['Chrom','Position','EditingFrequency']].rename(columns={'EditingFrequency':'EditingFrequencyKnown'})
    pooled = sites.loc[sites.RefBase.eq('A')].merge(reference, on=['Chrom','Position'], how='left', validate='many_to_one')
    pooled['PercentEditing'] = 100*pooled.EditingFrequency
    pooled['PercentEditingKnown'] = 100*pooled.EditingFrequencyKnown
    pooled = pooled.loc[pooled.Edited | pooled.KnownEditing].copy()
    pooled['Detected'] = np.select([pooled.Edited & pooled.KnownEditing, pooled.Edited], ['Both','De-novo'], default='Known')
    pearson_data = pooled[['PercentEditing','PercentEditingKnown']].dropna()
    correlation = pearsonr(pearson_data.PercentEditing, pearson_data.PercentEditingKnown)
    save('s2_pooled', pooled, {'known_file': str(known_path), 'known_sha256': digest(known_path),
        'source_cells': [79,80,85,89,91], 'pearson_n': len(pearson_data), 'pearson_r': float(correlation.statistic),
        'pearson_p': float(correlation.pvalue), 'missing_plot_values': 'fill with zero, native c/d definition'})
    # Current VP section "Editing levels in PCLO - squid - fixed": calls define
    # membership; measured frequencies remain present even without a platform call.
    columns = ['Chrom', 'Position', 'EditingFrequency', 'Edited']
    pclo = sites.loc[sites.Gene.eq('PCLO') & sites.RefBase.eq('A')]
    merged = pclo.loc[pclo.Platform.eq('PacBio'), columns].merge(
        pclo.loc[pclo.Platform.eq('Illumina'), columns], how='outer',
        on=['Chrom', 'Position'], suffixes=['_PacBio', '_Illumina'], validate='one_to_one')
    catalog = reference.loc[reference.Chrom.eq('comp141882_c0_seq14')].assign(Edited_Known=True)
    merged = merged.merge(catalog, how='outer', on=['Chrom', 'Position'], validate='one_to_one')
    flags = ['Edited_PacBio', 'Edited_Illumina', 'Edited_Known']
    merged[flags] = merged[flags].fillna(False).astype(bool)
    for source, destination in [('EditingFrequency_PacBio', 'PacBio'),
                                ('EditingFrequency_Illumina', 'Illumina')]:
        merged[destination] = merged.pop(source) * 100
    merged['EditingFrequencyKnown'] *= 100
    merged = merged[['Chrom', 'Position', *flags, 'PacBio', 'Illumina', 'EditingFrequencyKnown']]
    merged = merged.sort_values(['Chrom', 'Position']).reset_index(drop=True)
    current = merged.loc[merged.Edited_PacBio | merged.Edited_Illumina, 'Position']
    merged = merged.loc[merged.Position.between(current.min(), current.max()) & merged[flags].any(axis=1)]
    save('s2_pclo', merged, {})
    print('Pooled',len(pooled),'Pearson',correlation,flush=True)


    venns=[]
    for (platform,gene),data in squid.groupby(['Platform','Gene'],sort=False):
        edited=set(data.loc[data.Edited,'Position']);known=set(data.loc[data.KnownEditing,'Position'])
        venns.append({'Platform':platform,'Gene':gene,'NewOnly':len(edited-known),'KnownOnly':len(known-edited),'Both':len(edited&known)})
    edited=set(octopus.loc[octopus.CDS & octopus.EditedFinal,['Chrom','Position']].itertuples(index=False,name=None))
    known=set(octopus.loc[octopus.CDS & octopus.KnownEditing,['Chrom','Position']].itertuples(index=False,name=None))
    assert [len(edited-known),len(known-edited),len(edited&known)]==[4531,13161,7180]
    venns.append({'Platform':'Octopus','Gene':'Whole transcriptome','NewOnly':len(edited-known),'KnownOnly':len(known-edited),'Both':len(edited&known)})
    save('s2_venns', pd.DataFrame(venns), {})

    motif_rows=[];site_rows=[];frequency_frames=[]
    for organism,source,site_table in [('Squid',PROJECT/'D.pealeii/Annotations/orfs_squ.fa',squid),
                                      ('Octopus',PROJECT/'O.vulgaris/Annotations/orfs_oct.fa',octopus_motif)]:
        needed=set(site_table.Chrom)
        sequences={rec.id:rec.seq for rec in SeqIO.parse(source,'fasta') if rec.id in needed}
        assert set(sequences)==needed
        flag='Edited' if organism=='Squid' else 'EditedFinal'
        selected=site_table.loc[site_table[flag]].copy()
        if organism=='Squid':
            selected=selected.loc[selected.RefBase.eq('A')];selected['Strand']='+'
        else:selected['Platform']='Octopus'
        for platform,data in selected.groupby('Platform',sort=False):
            strings=[]
            for row in data.itertuples():
                end=len(sequences[row.Chrom]) if organism=='Squid' else row.ChromEnd
                if row.Position-1<0 or row.Position+2>end:continue
                seq=sequences[row.Chrom][row.Position-1:row.Position+2]
                if row.Strand=='-':seq=seq.reverse_complement()
                seq=str(seq).upper();assert len(seq)==3 and seq[1]=='A' and set(seq)<=set('ACGT')
                strings.append(seq)
                site_rows.append({'Platform':platform,'Gene':row.Gene,'Chrom':row.Chrom,'Position':row.Position,'Strand':row.Strand,'Sequence':seq})
            for offset in range(3):
                for base in 'ACGT':
                    count=sum(s[offset]==base for s in strings)
                    motif_rows.append({'Platform':platform,'Offset':offset-1,'Base':base,'Count':count,'Denominator':len(strings),'Probability':count/len(strings)})
            frequency_frames.append({'Platform':platform,'Fasta':str(source),'FASTA_SHA256':digest(source),'Sites':len(strings)})
    motifs=pd.DataFrame(motif_rows)
    save('s2_motif_sites',pd.DataFrame(site_rows),{'selection':'Original Edited (squid) / EditedFinal in original BQ30.AHL.BHAfterNoise.3 cohort (octopus), strand-aware ±1 base; native edge exclusion', 'fasta_bindings':frequency_frames})
    save('s2_motifs',motifs,{'source':'Code/EditingUtils/logo.py make_freq_df: nucleotide counts / number of equal-length sequences','fasta_bindings':frequency_frames})


    # Display-only horizontal jitter comes from the approved native SVG; Y is native.
    root = ET.parse(distribution_svg)
    ns = {'s': 'http://www.w3.org/2000/svg'}
    groups = [g for g in root.findall('.//s:g', ns) if g.get('class', '').startswith('trace boxes')]
    observations, boxes = [], []
    assert len(groups) == 3
    for i, (label, group) in enumerate(zip(['Both', 'De-novo', 'Known'], groups, strict=True)):
        table = pooled.loc[pooled.Detected.eq(label)].copy()
        values = table['PercentEditingKnown' if label == 'Known' else 'PercentEditing'].fillna(0).to_numpy()
        order = np.argsort(values, kind='stable')
        values, table = values[order], table.iloc[order]
        points = np.array([list(map(float, re.search(r'translate\(([^)]+)\)', p.get('transform'))[1].split(',')))
            for p in group.findall('.//s:path[@class="point"]', ns)])
        np.testing.assert_allclose(points[:, 1], 320 - values * 3.2, atol=.005001, rtol=0)
        x = (points[:, 0] - (74.17 + 146.33 * i)) / 146.33 + i
        observations.append(table.assign(Detected=label, DisplayValue=values, DisplayX=x))
        q1, median, q3 = np.quantile(values, [.25, .5, .75], method='hazen')
        low = values[values >= q1 - 1.5 * (q3 - q1)].min()
        high = values[values <= q3 + 1.5 * (q3 - q1)].max()
        boxes.append({'Detected': label, 'Category': i, 'N': len(values), 'Q1': q1, 'Median': median, 'Q3': q3, 'Low': low, 'High': high})
    tables['s2_distribution'] = pd.concat(observations, ignore_index=True)
    tables['s2_boxes'] = pd.DataFrame(boxes)
    for name, table in tables.items():
        table.to_csv(output_dir / f'{name}.tsv', sep='\t', index=False, na_rep='NA', mode='x')
    return tables


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    sites = sub.add_parser('sites')
    sites.add_argument('--branch', choices=['squid', 'octopus', 'octopus_motif'], required=True)
    sites.add_argument('--output', type=Path, required=True)
    sites.add_argument('--project-root', type=Path, default=PROJECT)
    tables = sub.add_parser('tables')
    tables.add_argument('--sites', type=Path, required=True)
    tables.add_argument('--output', type=Path, required=True)
    tables.add_argument('--distribution-svg', type=Path, required=True)
    tables.add_argument('--project-root', type=Path, default=PROJECT)
    args = parser.parse_args()
    if args.command == 'sites': export_sites(args.branch, args.output, args.project_root)
    else: produce_s2(args.sites, args.output, args.distribution_svg, args.project_root)
