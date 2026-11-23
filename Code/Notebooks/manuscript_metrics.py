"""Native Tables S1/S3 PacBio metrics, retaining stage-specific read populations.

PB1 has one library per gene; PB3 uses original JR1/JR2/JR3 read identities.
These native legacy positions tables use Edited as their final stored flag and
have no EditedFinal column. It must match the actual read-matrix site panel.
"""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import pysam
from manuscript_inputs import NOTEBOOKS, PROJECT, selections


def primary_reads(path, chrom, start, end):
    """Unique primary mapped IDs, plus the subset overlapping the half-open ORF."""
    mapped, overlap = set(), set()
    with pysam.AlignmentFile(path, 'rb', threads=1) as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.flag & 2304:
                continue
            mapped.add(read.query_name)
            if read.reference_name == chrom and read.reference_start < end and read.reference_end > start:
                overlap.add(read.query_name)
    return mapped, overlap


def replicate_map(project_root):
    """Original sample identity comes from the three delivered unmapped BAMs."""
    mapping = {}
    for replicate in [1, 2, 3]:
        path = Path(project_root) / f'D.pealeii/Data/AdditionalRawWithUMIs/CCS/JR{replicate}.r841072252411_C01.hifireads.bam'
        with pysam.AlignmentFile(path, 'rb', threads=1, check_sq=False) as bam:
            for read in bam.fetch(until_eof=True):
                assert read.query_name not in mapping
                mapping[read.query_name] = replicate
    return mapping


def produce_metrics(output_dir, project_root=PROJECT):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=False)
    inputs = pd.read_table(NOTEBOOKS / 'manuscript_input_data/pacbio_metrics.tsv', keep_default_na=False)
    sample_by_old = replicate_map(project_root)
    metrics, read_metrics = [], []
    for item in inputs.itertuples(index=False):
        prefix = str(Path(project_root) / item.Prefix)
        reads = pd.read_table(prefix + '.reads.csv.gz', dtype={'Read': str}, keep_default_na=False)
        assert reads.Gene.eq(item.OriginalGene).all() and reads.Read.is_unique
        sites = reads.columns[6:].tolist()
        assert all(s.isdigit() for s in sites)
        assert reads[sites].isin([-1, 0, 1]).all().all()
        np.testing.assert_array_equal(reads[sites].eq(1).sum(axis=1), reads.EditedPositions)
        patterns = pd.read_table(prefix + '.unique_reads.csv.gz', usecols=['Transcript', 'Reads', 'NumOfReads'], dtype={'Transcript': str})
        proteins = pd.read_table(prefix + '.proteins.csv.gz', dtype={'Transcript': str}, keep_default_na=False)
        assert patterns.Transcript.is_unique and proteins.Transcript.is_unique
        np.testing.assert_array_equal(patterns.Transcript, proteins.Transcript)
        check = patterns.merge(proteins[['Transcript', 'Reads', 'NumOfReads']], on='Transcript', validate='one_to_one', suffixes=('_rna', '_protein'))
        assert check.Reads_rna.eq(check.Reads_protein).all() and check.NumOfReads_rna.eq(check.NumOfReads_protein).all()
        aa_columns = proteins.columns[12:].tolist()
        minima = sum(proteins[c].map(lambda aas: c.split('(')[1][:-1] not in aas.split(',')) for c in aa_columns)
        np.testing.assert_array_equal(minima, proteins.MinNonSyns)
        per_read = proteins[['Transcript', 'Reads', 'MinNonSyns']].assign(Read=proteins.Reads.str.split(',')).explode('Read').drop(columns='Reads')
        assert per_read.Read.is_unique and set(per_read.Read) == set(reads.Read)
        per_read = per_read.merge(reads[['Read', 'EditedPositions']], on='Read', validate='one_to_one')
        mapping = pd.read_table(prefix + '.OldToNewReads.csv.gz', dtype=str, keep_default_na=False).set_index('NewRead')
        assert mapping.index.is_unique and mapping.OldRead.is_unique
        per_read['OldRead'] = per_read.Read.map(mapping.OldRead)
        assert per_read.OldRead.notna().all()
        per_read['Replicate'] = per_read.OldRead.map(sample_by_old) if item.Dataset == 'PB3' else 'Pooled'
        assert per_read.Replicate.notna().all()
        # S3 uses the original pandas grouped mean/sample-SD implementation.
        aggregations = dict(MeanEditedSitesPerRead=('EditedPositions', 'mean'),
            SDEditedSitesPerRead=('EditedPositions', lambda x: x.std(ddof=1)),
            MeanMinimalRecodingSubstitutionsPerRead=('MinNonSyns', 'mean'),
            SDMinimalRecodingSubstitutionsPerRead=('MinNonSyns', lambda x: x.std(ddof=1)))
        grouped_metrics = per_read.groupby('Replicate').agg(**aggregations)
        pooled_metrics = per_read.assign(Group=0).groupby('Group').agg(**aggregations).iloc[0]
        known, detected, pileup = set(), set(), set()
        columns = ['Position', 'CDS', 'KnownEditing', item.DetectedFlag, 'Reads']
        for chunk in pd.read_table(prefix + '.positions.csv.gz', usecols=columns, chunksize=200):
            assert chunk.Position.between(item.Start, item.End - 1).all() and chunk.CDS.all()
            known.update(chunk.loc[chunk.KnownEditing, 'Position'])
            detected.update(chunk.loc[chunk[item.DetectedFlag], 'Position'])
            for ids in chunk.Reads:
                pileup.update(ids.split(','))
        assert detected == set(map(int, sites)) and set(reads.Read) <= pileup
        old_pileup = set(mapping.loc[list(pileup), 'OldRead'])
        aligned, aligned_orf = primary_reads(Path(project_root) / item.RawBam, item.Chrom, item.Start, item.End)
        filtered, filtered_orf = primary_reads(prefix + '.bam', item.Chrom, item.Start, item.End)
        assert old_pileup <= filtered_orf <= filtered <= aligned
        graph = pd.read_table(Path(project_root) / item.Graph, usecols=['Fraction', 'NumUniqueSamples'])
        maximum = int(graph.loc[graph.Fraction.eq(1), 'NumUniqueSamples'].max())
        total_ccs = None
        if item.RawCCSBam:
            with pysam.AlignmentFile(Path(project_root) / item.RawCCSBam, 'rb', check_sq=False, threads=1) as bam:
                total_ccs = len({read.query_name for read in bam.fetch(until_eof=True)})
        replicates = [1, 2, 3, 'Pooled'] if item.Dataset == 'PB3' else ['Pooled']
        for replicate in replicates:
            selected = per_read if replicate == 'Pooled' else per_read.loc[per_read.Replicate.eq(replicate)]
            statistics = pooled_metrics if replicate == 'Pooled' else grouped_metrics.loc[replicate]
            def stage_count(ids):
                return len(ids) if replicate == 'Pooled' else sum(sample_by_old[r] == replicate for r in ids)
            metrics.append({'Dataset': item.Dataset, 'Gene': item.Gene, 'Replicate': replicate,
                'Chrom': item.Chrom, 'Start': item.Start, 'End': item.End, 'ORFAminoAcids': (item.End - item.Start) // 3,
                'KnownSites': len(known) if replicate == 'Pooled' and item.Gene != 'ADAR1' else np.nan,
                'DetectedSites': len(detected) if replicate == 'Pooled' else np.nan,
                'EditableAAs': len(aa_columns) if replicate == 'Pooled' else np.nan,
                'MeanEditedSitesPerRead': statistics.MeanEditedSitesPerRead,
                'SDEditedSitesPerRead': statistics.SDEditedSitesPerRead,
                'MeanMinimalRecodingSubstitutionsPerRead': statistics.MeanMinimalRecodingSubstitutionsPerRead,
                'SDMinimalRecodingSubstitutionsPerRead': statistics.SDMinimalRecodingSubstitutionsPerRead,
                'TotalCCSReads': total_ccs, 'AlignedORFReads': stage_count(aligned_orf),
                'FilteredMappedReads': stage_count(filtered), 'FilteredORFReads': stage_count(filtered_orf),
                'PileupReads': stage_count(old_pileup), 'AnalysisReads': len(selected),
                'UniqueReads': selected.Transcript.nunique(),
                'DistinctProteins': maximum if replicate == 'Pooled' else np.nan})
        read_metrics.append(per_read.assign(Dataset=item.Dataset, Gene=item.Gene))
        print(item.Dataset, item.Gene, len(reads), 'original analysis reads', flush=True)
    table = pd.DataFrame(metrics)
    table.to_csv(output_dir / 'pacbio_native_metrics.tsv', sep='\t', index=False, na_rep='NA', mode='x')
    pd.concat(read_metrics, ignore_index=True).to_csv(output_dir / 'pacbio_original_read_metrics.tsv.gz', sep='\t', index=False, mode='x')
    return table


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--project-root', type=Path, default=PROJECT)
    args = parser.parse_args()
    produce_metrics(args.output, args.project_root)
