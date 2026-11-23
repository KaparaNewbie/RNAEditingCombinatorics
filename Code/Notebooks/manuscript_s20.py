"""Original-RNA chimera classification for the accepted S20 selected isoforms.

X=.01 and Y=.30 use ties-to-even rounding. A rare protein is chimeric if
any originally supporting RNA pattern has an eligible prefix/suffix witness.
Reassigned expression contributes to ranking, never to original RNA identities.
"""
from pathlib import Path
import argparse
import hashlib
import numpy as np
import pandas as pd
from manuscript_inputs import PROJECT, selections, original_rna_patterns, ids_sha256
from manuscript_figures.styles import S20_JITTER, S20_COLORS, S20_SOURCE_ORDER, S20_DISPLAY_ORDER

def common_pattern_index(common):
    """One bit per original common RNA; bit order retains the original witness order."""
    masks=[[int.from_bytes(np.packbits(common[:,column]==state,bitorder='little').tobytes(),'little')
            for state in [-1,0,1]] for column in range(common.shape[1])]
    return masks,len(common),(1<<len(common))-1

def classify_indexed(index, rare, soft):
    """Count the same E <= B pairs through exact integer set intersections."""
    masks,n,all_reads=index;length=len(rare)
    compatible=[all_reads if soft and state==-1 else
                (masks[column][0] | masks[column][state+1]) if soft else masks[column][state+1]
                for column,state in enumerate(rare)]
    prefix=[n]+[0]*length;state=all_reads;left_bits=state;B=0
    for position,allowed in enumerate(compatible):
        state &= allowed;prefix[position+1]=state.bit_count()
        if not state:break
        left_bits=state;B=position+1
    assert prefix[length]==0,'Original producer forbids a self-pair; inspect cohort before proceeding'
    suffix=[0]*(length+2);suffix[length+1]=n;state=all_reads;right_bits=state;E=length+1
    for position in range(length-1,-1,-1):
        state &= compatible[position];suffix[position+1]=state.bit_count()
        if not state:break
        right_bits=state;E=position+1
    # prefix[b]-prefix[b+1] counts B==b; suffix[b] counts E<=b.
    count=sum((prefix[b]-prefix[b+1])*suffix[b] for b in range(1,length))
    left=(left_bits & -left_bits).bit_length()-1
    right=(right_bits & -right_bits).bit_length()-1
    if count:assert left!=right and E<=B
    return count,left,right,B,E


def classify_original_rna(selection, assignment, rna, assignment_sha):
    n_common, n_rare = round(.01 * len(assignment)), round(.30 * len(assignment))
    common_ids = assignment.Protein.iloc[:n_common].tolist()
    rare_ids = assignment.Protein.iloc[-n_rare:].tolist()
    common_indices = np.flatnonzero(np.isin(rna['Protein'], common_ids))
    rare_indices = np.flatnonzero(np.isin(rna['Protein'], rare_ids))
    assert set(rna['Protein'][rare_indices]) == set(rare_ids)
    index = common_pattern_index(rna['states'][common_indices])
    rows = []
    for position in rare_indices:
        for soft in [False, True]:
            count, left, right, prefix_end, suffix_start = classify_indexed(index, rna['states'][position], soft)
            rows.append({'Platform': selection.Platform + '1', 'Gene': selection.Gene,
                'Solution': selection.Solution, 'AssignmentSHA256': assignment_sha,
                'Protein': rna['Protein'][position], 'UniqueRead': rna['UniqueRead'][position],
                'IsSoftComparison': soft, 'IsChimeric': count > 0, 'NumOfChimericCombinations': count,
                'ActualXCommonProteins': n_common, 'ActualYRareProteins': n_rare,
                'WitnessPrefixRead': rna['UniqueRead'][common_indices[left]] if count else '',
                'WitnessSuffixRead': rna['UniqueRead'][common_indices[right]] if count else '',
                'WitnessPrefixEnd': prefix_end if count else -1,
                'WitnessSuffixStart': suffix_start if count else -1})
    return pd.DataFrame(rows)


def common_expression(selection, assignment, project_root=PROJECT):
    """Preserve the original normalization BEFORE expression-rank sorting.

    The denominator's floating-point reduction order is scientifically relevant
    to exact reproduction. Align the resulting percentages to the frozen IDs;
    never change a tied canonical selection to reproduce a display value.
    """
    from manuscript_inputs import verified_file
    path = verified_file(project_root, selection.ExpressionFile, selection.ExpressionSHA256)
    parts = []
    for chunk in pd.read_table(path, usecols=['Protein', '#Solution', 'TotalWeightedSupportingReads'],
            dtype={'Protein': str, '#Solution': str}, chunksize=10000):
        parts.append(chunk.loc[chunk['#Solution'].eq(selection.Solution)])
    native = pd.concat(parts, ignore_index=True)
    assert native.Protein.is_unique and set(native.Protein) == set(assignment.Protein)
    native['RelativeExpression'] = 100 * native.TotalWeightedSupportingReads / native.TotalWeightedSupportingReads.sum()
    ranked = native.set_index('Protein').loc[assignment.Protein]
    np.testing.assert_array_equal(ranked.TotalWeightedSupportingReads, assignment.TotalWeightedSupportingReads)
    return float(ranked.RelativeExpression.iloc[:round(.01 * len(assignment))].sum())


def produce_s20(assignments_dir, output_dir, project_root=PROJECT, *, work_dir=None):
    """Write new RNA, classification and final-table files; never replace outputs."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=False)
    work_dir = Path(work_dir) if work_dir is not None else output_dir
    if work_dir != output_dir:
        work_dir.mkdir(parents=True, exist_ok=False)
    expression_rows, counts = [], {}
    ordinary = selections().loc[lambda t: t.Metric.eq('ordinary')]
    # The publication box plot keeps Illumina first and three settings per gene.
    ordinary = pd.concat([ordinary.loc[ordinary.Platform.eq(p)] for p in ['Illumina', 'PacBio']])
    for selection in ordinary.itertuples(index=False):
        stem = f'{selection.Platform}.{selection.Gene}.ordinary'
        path = Path(assignments_dir) / f'{stem}.tsv'
        assignment_sha = hashlib.sha256(path.read_bytes()).hexdigest()
        assert assignment_sha == selection.AssignmentSHA256
        assignment = pd.read_table(path, dtype={'Protein': str, '#Solution': str}, float_precision='round_trip')
        assert ids_sha256(assignment.Protein.tolist()) == selection.ProteinOrderSHA256
        assert assignment['#Solution'].eq(selection.Solution).all()
        rna = original_rna_patterns(selection, assignment, project_root)
        np.savez_compressed(work_dir / f'{stem}.original_rna.npz', **rna)
        classified = classify_original_rna(selection, assignment, rna, assignment_sha)
        classified.to_csv(work_dir / f'{stem}.classification.tsv.gz', sep='\t', index=False, na_rep='NA')
        native_platform = selection.Platform + '1'
        n_common, n_rare = round(.01 * len(assignment)), round(.30 * len(assignment))
        value = common_expression(selection, assignment, project_root)
        for fraction in [.1, .2, .3]:
            position = sum(r['Platform'] == native_platform for r in expression_rows)
            expression_rows.append({'Platform': native_platform, 'Gene': selection.Gene,
                'YRareProteins': fraction, 'XCommonProteins': .01,
                'TotalDistinctProteins': len(assignment), 'ActualXCommonProteins': n_common,
                'PlotObservation': position, 'CumulativeExpression': value,
                'JitterOffset': S20_JITTER[native_platform][position],
                'Solution': selection.Solution, 'AssignmentSHA256': assignment_sha})
        for soft, comparison in [(False, 'Strict comparison'), (True, 'Soft comparison')]:
            by_protein = classified.loc[classified.IsSoftComparison.eq(soft)].groupby('Protein').IsChimeric.any()
            assert len(by_protein) == n_rare
            counts[(comparison, native_platform, selection.Gene)] = {
                'TotalDistinctProteins': len(assignment), 'ActualYRareProteins': n_rare,
                'RecoveredChimericCount': int(by_protein.sum()),
                'ChimericPercent': 100 * int(by_protein.sum()) / n_rare,
                'Solution': selection.Solution, 'AssignmentSHA256': assignment_sha}
        print(stem, len(rna['Protein']), 'original RNA patterns classified', flush=True)
    expression = pd.DataFrame(expression_rows)
    boxes = []
    for order, (platform, data) in enumerate(expression.groupby('Platform', sort=False)):
        values = data.CumulativeExpression.to_numpy()
        q1, median, q3 = np.quantile(values, [.25, .5, .75], method='hazen')
        inliers = values[(values >= q1 - 1.5 * (q3 - q1)) & (values <= q3 + 1.5 * (q3 - q1))]
        boxes.append({'Platform': platform, 'Order': order, 'Color': S20_COLORS[platform],
            'Q1': q1, 'Median': median, 'Q3': q3, 'WhiskerLow': inliers.min(), 'WhiskerHigh': inliers.max(),
            'NPlotRows': len(values), 'NGenes': data.Gene.nunique()})
    bars = []
    for comparison in ['Strict comparison', 'Soft comparison']:
        for source_order, (platform, gene) in enumerate(S20_SOURCE_ORDER):
            record = counts[(comparison, platform, gene)]
            bars.append({'Comparison': comparison, 'Platform': platform, 'Gene': gene,
                'SourceOrder': source_order, 'DisplayOrder': S20_DISPLAY_ORDER.index((platform, gene)),
                **{key: record[key] for key in ['TotalDistinctProteins', 'ActualYRareProteins',
                    'RecoveredChimericCount', 'ChimericPercent']}, 'Color': S20_COLORS[platform],
                'Solution': record['Solution'], 'AssignmentSHA256': record['AssignmentSHA256']})
    for name, table in [('expression', expression), ('boxes', pd.DataFrame(boxes)), ('chimeras', pd.DataFrame(bars))]:
        table.to_csv(output_dir / f'figure_s20_{name}.tsv', sep='\t', index=False, mode='x')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--assignments', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True, help='New output directory')
    parser.add_argument('--project-root', type=Path, default=PROJECT)
    parser.add_argument('--work', type=Path, required=True, help='New intermediate directory outside figure_sources')
    args = parser.parse_args()
    produce_s20(args.assignments, args.output, args.project_root, work_dir=args.work)
