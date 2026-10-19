"""Descriptive module-age comparison under recovered historical binning.

Portable transcription, not a clustering stability or prediction benchmark.
Inputs must already be preprocessed and standardized. No study data are bundled.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.cluster import KMeans
from threadpoolctl import threadpool_limits


def historical_bins():
    """Retain historical overlap at 45-46 and exclusion of exact age 90."""
    definitions = [(a, a + 1) for a in range(20, 30)]
    definitions += [(a, a + 2) for a in range(30, 45, 2)]
    definitions += [(a, a + 1) for a in range(45, 90)]
    centers = np.array([(a + b) / 2 for a, b in definitions], dtype=float)
    return definitions, centers


def interp_nan(y, x):
    """Historical fallback and NumPy endpoint extension, unchanged."""
    y = y.copy()
    good = ~np.isnan(y)
    if good.sum() < 2:
        y[~good] = 0.0
        return y
    y[~good] = np.interp(x[~good], x[good], y[good])
    return y


def make_curves(matrix, ages):
    definitions, centers = historical_bins()
    curves = np.full((matrix.shape[1], len(definitions)), np.nan)
    for j, (lower, upper) in enumerate(definitions):
        mask = (ages >= lower) & (ages < upper)
        if mask.sum() >= 2:
            curves[:, j] = np.median(matrix[mask], axis=0)
    result = np.vstack([interp_nan(row, centers) for row in curves])
    if not np.isfinite(result).all():
        raise ValueError('Non-finite clustering curves.')
    return result


def module_correlations(matrix, ages, labels):
    rows = []
    for label in np.unique(labels):
        mask = labels == label
        score = matrix[:, mask].mean(axis=1)
        rho = float(spearmanr(ages, score).statistic)
        rows.append({'cluster': int(label), 'n_analytes': int(mask.sum()),
                     'rho_module_age': rho})
    return pd.DataFrame(rows)


def select_directions(table):
    valid = table[np.isfinite(table.rho_module_age)]
    positive = valid[valid.rho_module_age > 0]
    negative = valid[valid.rho_module_age < 0]
    if positive.empty or negative.empty:
        return None, None, np.nan
    increasing = positive.loc[positive.rho_module_age.idxmax()]
    decreasing = negative.loc[negative.rho_module_age.idxmin()]
    mean_absolute = (abs(increasing.rho_module_age) + abs(decreasing.rho_module_age)) / 2
    return increasing, decreasing, float(mean_absolute)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--matrix', type=Path, required=True)
    parser.add_argument('--ages', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--threads', type=int, default=1)
    args = parser.parse_args()
    expression = pd.read_csv(args.matrix, index_col=0)
    ages = pd.read_csv(args.ages)
    if not expression.index.is_unique or not expression.columns.is_unique:
        raise ValueError('Expression keys and assay columns must be unique.')
    if not ages.sample_id.is_unique:
        raise ValueError('Age sample keys must be unique.')
    aligned = ages.set_index('sample_id').loc[expression.index, 'age'].to_numpy(float)
    matrix = expression.to_numpy(float)
    if not np.isfinite(aligned).all() or not np.isfinite(matrix).all():
        raise ValueError('Expected finite ages and an already processed matrix.')
    eligible = (aligned >= 20) & (aligned <= 90)
    matrix, aligned = matrix[eligible], aligned[eligible]
    if len(aligned) < 3 or matrix.shape[1] < 25:
        raise ValueError('Insufficient rows or features for K=2-25.')
    curves = make_curves(matrix, aligned)
    args.output.mkdir(parents=True, exist_ok=False)
    summaries = []
    with threadpool_limits(limits=args.threads):
        for k in range(2, 26):
            labels = KMeans(n_clusters=k, random_state=0, n_init=10,
                            algorithm='lloyd').fit_predict(curves)
            table = module_correlations(matrix, aligned, labels)
            increasing, decreasing, score = select_directions(table)
            pd.DataFrame({'feature': expression.columns, 'cluster': labels}).to_csv(
                args.output / f'PRIVATE_assignments_K{k}.csv', index=False)
            table.to_csv(args.output / f'module_age_correlations_K{k}.csv', index=False)
            row = {'K': k, 'mean_absolute_module_age_rho': score,
                   'both_directions_present': increasing is not None}
            for prefix, module in [('increasing', increasing), ('decreasing', decreasing)]:
                row[prefix + '_cluster'] = None if module is None else int(module.cluster)
                row[prefix + '_n_analytes'] = None if module is None else int(module.n_analytes)
                row[prefix + '_rho'] = None if module is None else float(module.rho_module_age)
            summaries.append(row)
    pd.DataFrame(summaries).to_csv(args.output / 'K_resolution_comparison.csv', index=False)


if __name__ == '__main__':
    main()
