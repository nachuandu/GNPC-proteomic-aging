"""Recovered age-bin medians -> interpolation -> direct k-means; no PCA."""
import argparse
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from threadpoolctl import threadpool_limits
from analysis_core import read_matrix, align_metadata, output_folder
from compare_k_resolutions import make_curves, module_correlations, select_directions


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--matrix', required=True)
    p.add_argument('--metadata', required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--k', type=int, default=20)
    p.add_argument('--seed', type=int, default=0)
    p.add_argument('--clusters', help='Explicit numbered clusters, e.g. 9,19. Labels are run-specific.')
    a = p.parse_args()
    X = read_matrix(a.matrix)
    meta = align_metadata(X, a.metadata, ['age'])
    use = meta.age.between(20, 90)
    values, ages = X.loc[use].to_numpy(), meta.loc[use, 'age'].to_numpy()
    with threadpool_limits(limits=1):
        labels = KMeans(n_clusters=a.k, random_state=a.seed, n_init=10,
                        algorithm='lloyd').fit_predict(make_curves(values, ages))
    modules = module_correlations(values, ages, labels)
    if a.clusters:
        chosen = [int(v) for v in a.clusters.split(',')]
        if not set(chosen) <= set(labels):
            raise ValueError('Requested cluster absent.')
    else:
        up, down, _ = select_directions(modules)
        if up is None:
            raise ValueError('Both positive and negative modules are required.')
        chosen = [int(up.cluster), int(down.cluster)]
    assignments = pd.DataFrame({'feature': X.columns, 'cluster': labels})
    assignments['selected'] = np.isin(labels, chosen)
    out = output_folder(a.output)
    assignments.to_csv(out / 'PRIVATE_cluster_assignments.csv', index=False)
    assignments.loc[assignments.selected, ['feature']].to_csv(out / 'PRIVATE_selected_features.csv', index=False)
    modules.to_csv(out / 'module_age_correlations.csv', index=False)


if __name__ == '__main__':
    main()
