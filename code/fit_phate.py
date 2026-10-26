"""PHATE embedding adapter. Inputs are already standardized; no PCA clustering.

PHATE's internal n_pca parameter is distinct from protein-module clustering.
Choose the curve backend explicitly after this step; see PROVENANCE.md.
"""
import argparse
import json
import importlib.metadata
import numpy as np
import pandas as pd
from analysis_core import read_matrix, align_metadata, output_folder


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--matrix', required=True)
    p.add_argument('--metadata', required=True)
    p.add_argument('--features', required=True, help='CSV with feature column, in analysis order')
    p.add_argument('--output', required=True)
    p.add_argument('--seed', required=True, type=int)
    p.add_argument('--knn', type=int, default=50)
    p.add_argument('--t', type=int, default=50)
    p.add_argument('--threads', type=int, default=1)
    a = p.parse_args()
    import phate
    from threadpoolctl import threadpool_limits
    X = read_matrix(a.matrix)
    meta = align_metadata(X, a.metadata, ['age'])
    features = pd.read_csv(a.features).feature
    if features.isna().any() or not features.is_unique:
        raise ValueError('Features must be unique and nonmissing.')
    values = X.loc[:, features].to_numpy()
    op = phate.PHATE(n_components=2, knn=a.knn, t=a.t, decay=40, gamma=1,
                     n_pca=100, n_landmark=2000, knn_dist='euclidean',
                     mds_dist='euclidean', mds='metric', mds_solver='smacof',
                     n_jobs=a.threads, random_state=a.seed, verbose=0)
    with threadpool_limits(limits=a.threads):
        xy = op.fit_transform(values)
    if not np.isfinite(xy).all():
        raise ValueError('Nonfinite PHATE coordinates.')
    out = output_folder(a.output)
    pd.DataFrame({'sample_id': X.index, 'age': meta.age,
                  'PHATE1': xy[:, 0], 'PHATE2': xy[:, 1]}).to_csv(
        out / 'PRIVATE_embedding.csv', index=False)
    (out / 'parameters.json').write_text(json.dumps({'parameters': op.get_params(),
        'phate_version': importlib.metadata.version('phate'), 'seed': a.seed,
        'n_samples': len(X), 'n_features': len(features)}, indent=2, default=str))


if __name__ == '__main__':
    main()
