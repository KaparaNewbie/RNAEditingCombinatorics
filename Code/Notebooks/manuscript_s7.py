"""Corrected binary-recoding matrices and complete S7 silhouette scores.

One observation is a canonical selected isoform. Features are unweighted binary
recoding codons, with the accepted positive-frequency filter and column means.
Each branch has one shared RandomState(1892), in the TSV's gene/cohort order.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import argparse
import sys
import numpy as np
import pandas as pd
from manuscript_inputs import PROJECT, selections, selected_amino_acids, ids_sha256

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from EditingUtils.recoding import _binary_recoding_matrix


def prepare_matrices(assignments_dir, output_dir, project_root=PROJECT):
    """Build features from native AA sets, verifying selected identity and order."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=False)
    for selection in selections().itertuples(index=False):
        stem = f'{selection.Platform}.{selection.Gene}.{selection.Metric}'
        assignment = pd.read_table(Path(assignments_dir) / f'{stem}.tsv',
            dtype={'Protein': str, '#Solution': str}, float_precision='round_trip')
        assert len(assignment) == selection.SelectedProteins
        assert assignment['#Solution'].eq(selection.Solution).all()
        assert ids_sha256(assignment.Protein.tolist()) == selection.ProteinOrderSHA256
        amino_acids = selected_amino_acids(selection, assignment, project_root)
        matrix = _binary_recoding_matrix(amino_acids, impute_nans=True)
        assert np.isfinite(matrix).all().all()
        np.savez_compressed(output_dir / f'{stem}.npz', X=matrix.to_numpy(float),
            Protein=np.asarray(assignment.Protein, dtype=str), columns=np.asarray(matrix.columns, dtype=str))
        print(stem, list(matrix.shape), flush=True)


def _reduce_silhouettes(distances, start, labels):
    # Same sklearn reduction as silhouette_samples, reused across the k fits.
    from sklearn.metrics.cluster._unsupervised import _silhouette_reduce
    counts = np.bincount(labels)
    intra, inter = _silhouette_reduce(distances, start, labels, counts)
    with np.errstate(divide='ignore', invalid='ignore'):
        intra /= (counts - 1)[labels[start:start + len(distances)]]
        values = (inter - intra) / np.maximum(intra, inter)
    return np.nan_to_num(values)


def all_observation_silhouettes(matrix, label_sets):
    """Mean over every observation; no subsampling or expression weighting."""
    from sklearn.metrics import pairwise_distances_chunked
    from sklearn.preprocessing import LabelEncoder
    from threadpoolctl import threadpool_limits
    labels = [LabelEncoder().fit_transform(values) for values in label_sets]
    scores = np.empty((len(labels), len(matrix)), dtype='float64')
    start = 0
    # Start Python with OPENBLAS_NUM_THREADS=7: the local OpenBLAS build cannot
    # safely grow a pool initialized at one thread. Fitting still uses one thread.
    with threadpool_limits(limits=7), ThreadPoolExecutor(max_workers=7) as pool:
        for distances in pairwise_distances_chunked(matrix, metric='euclidean', working_memory=256):
            stop = start + len(distances)
            scores[:, start:stop] = np.asarray(list(pool.map(
                lambda label: _reduce_silhouettes(distances, start, label), labels)))
            start = stop
    assert start == len(matrix)
    return scores.mean(axis=1)


def branch_scores(matrices_dir, branch, output_path):
    """Fit the accepted k series in order, with one fitting thread and shared RNG."""
    import sklearn
    from sklearn.cluster import MiniBatchKMeans
    from threadpoolctl import threadpool_limits
    assert sklearn.__version__ == '1.9.0', 'Accepted S7 fits require scikit-learn 1.9.0'
    output_path = Path(output_path)
    assert not output_path.exists(), output_path
    output_path.parent.mkdir(parents=True, exist_ok=True)
    platform, kind = branch.split('_')
    selected = selections()
    selected = selected.loc[selected.Platform.eq(platform) &
        (selected.Metric.eq('ordinary') if kind == 'ordinary' else selected.Metric.ne('ordinary'))]
    sizes = list(range(10, 501, 10)) if kind == 'ordinary' else list(range(3, 41))
    rng = np.random.RandomState(1892)
    records = []
    for selection in selected.itertuples(index=False):
        stem = f'{selection.Platform}.{selection.Gene}.{selection.Metric}'
        with np.load(Path(matrices_dir) / f'{stem}.npz', allow_pickle=False) as saved:
            matrix = saved['X']
            assert ids_sha256(saved['Protein'].tolist()) == selection.ProteinOrderSHA256
        assert len(matrix) == selection.SelectedProteins and np.isfinite(matrix).all()
        labels = []
        for k in sizes:
            with threadpool_limits(limits=1):
                model = MiniBatchKMeans(n_clusters=k, init='k-means++', n_init=5,
                                       verbose=0, random_state=rng).fit(matrix)
            labels.append(model.labels_)
        print(stem, 'fits complete', flush=True)
        means = all_observation_silhouettes(matrix, labels)
        for k, score in zip(sizes, means, strict=True):
            records.append({'Panel': 'a' if kind == 'ordinary' else 'b', 'Platform': platform,
                'Gene': selection.Gene, 'Dissimilarity': selection.Metric, 'ClusterSize': k,
                'MeanSilhouetteScore': float(score)})
        # Partial progress is inspectable if this bounded branch is interrupted.
        pd.DataFrame(records).to_csv(output_path, sep='\t', index=False)
        print(stem, 'all-observation scores complete', flush=True)
    assert len(records) == len(selected) * len(sizes)
    return pd.DataFrame(records)


def combine_scores(branch_files, output_path):
    """Join the three complete score branches without changing values or row order."""
    assert len(branch_files) == 3
    parts = [pd.read_table(p, float_precision='round_trip') for p in branch_files]
    assert [len(p) for p in parts] == [100, 950, 152]
    table = pd.concat(parts, ignore_index=True)
    keys = ['Panel', 'Platform', 'Gene', 'Dissimilarity', 'ClusterSize']
    assert not table.duplicated(keys).any() and np.isfinite(table.MeanSilhouetteScore).all()
    expected = []
    for selection in selections().itertuples(index=False):
        ordinary = selection.Metric == 'ordinary'
        sizes = range(10, 501, 10) if ordinary else range(3, 41)
        expected.extend([('a' if ordinary else 'b', selection.Platform, selection.Gene, selection.Metric, k) for k in sizes])
    assert set(table[keys].itertuples(index=False, name=None)) == set(expected)
    output_path = Path(output_path)
    assert not output_path.exists(), output_path
    output_path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output_path, sep='\t', index=False, mode='x')
    return table


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    prepare = sub.add_parser('prepare')
    prepare.add_argument('--assignments', type=Path, required=True)
    prepare.add_argument('--output', type=Path, required=True)
    prepare.add_argument('--project-root', type=Path, default=PROJECT)
    scores = sub.add_parser('scores')
    scores.add_argument('--matrices', type=Path, required=True)
    scores.add_argument('--output', type=Path, required=True)
    scores.add_argument('--branch', choices=['PacBio_ordinary', 'PacBio_dissimilar', 'Illumina_ordinary'], required=True)
    combine = sub.add_parser('combine')
    combine.add_argument('--branches', nargs=3, type=Path, required=True, metavar=('PB', 'IL', 'DISSIMILAR'))
    combine.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'prepare': prepare_matrices(args.assignments, args.output, args.project_root)
    elif args.command == 'scores': branch_scores(args.matrices, args.branch, args.output)
    else: combine_scores(args.branches, args.output)
