"""Shared numerical backend of the completed single-visit IADRC analysis.

Extracted without altering the PHATE or principal-curve operations.
This is not the historical R princurve implementation.
"""
import contextlib
import io
import time
import warnings
import numpy as np
from scipy.stats import spearmanr

def correlation(x, y):
    return float(spearmanr(x, y).statistic)


def percentile(reference_age, reference_gap, query_age, query_gap, h=5.):
    weights = np.exp(-.5 * ((query_age[:, None] - reference_age[None, :]) / h) ** 2)
    comparisons = (reference_gap[None, :] < query_gap[:, None]).astype(float)
    comparisons += .5 * (reference_gap[None, :] == query_gap[:, None])
    denominator = weights.sum(axis=1)
    assert np.all(denominator > 0)
    result = 100 * (weights * comparisons).sum(axis=1) / denominator
    effective_n = denominator ** 2 / np.sum(weights ** 2, axis=1)
    assert np.isfinite(result).all()
    return result, effective_n


def operator(threads, seed):
    import phate
    return phate.PHATE(n_components=2, knn=50, t=50, decay=40, gamma=1,
                       n_pca=100, n_landmark=2000, knn_dist='euclidean',
                       mds_dist='euclidean', mds='metric', mds_solver='smacof',
                       n_jobs=threads, random_state=seed, verbose=0)


def fit_once(X, age_gnpc, pc, threads, seed):
    started = time.monotonic()
    op = operator(threads, seed)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        xy = np.asarray(op.fit_transform(X))
    assert xy.shape == (len(X), 2) and np.isfinite(xy).all()
    with contextlib.redirect_stdout(io.StringIO()):
        raw, curve, projection, distance, curve_pt, report = pc.fit_principal_curve(
            xy, max_iter=1000)
    assert report['converged'], 'Principal curve did not converge'
    flip = correlation(age_gnpc, raw[:len(age_gnpc)]) < 0
    pt = 1 - raw if flip else raw
    if flip:
        curve, curve_pt = curve[::-1], (1 - curve_pt)[::-1]
    assert np.isfinite(pt).all() and pt.min() == 0 and pt.max() == 1
    audit = dict(seconds=time.monotonic() - started, seed=seed,
                 orientation_uses_GNPC_age_only=True,
                 principal_curve_converged=bool(report['converged']),
                 principal_curve_iterations=report['iterations'],
                 warnings=sorted(set(str(w.message) for w in caught)))
    return pt, xy, curve, curve_pt, op, audit

