"""Age-local percentile sensitivity across complete cumulative joint fits.

No study data are included. Inputs and outputs belong in an approved environment.
"""
import argparse
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def age_rank_diagnostics(G, metrics, results):
    """Compare age-local GNPC ranks as well as overall score ranks."""
    source = Path(__file__).resolve().parent
    gm = load_module(source / 'GAM_methods.py', 'local_gam_methods')
    fixed = load_module(source / 'fixed_GAM_methods.py', 'local_fixed_methods')
    from scipy.interpolate import PchipInterpolator
    ages = G.age.to_numpy()
    weights = np.exp(-.5 * ((ages[:, None] - ages[None, :]) / 5.) ** 2)
    # Each GNPC diagnostic query is excluded from its own reference CDF.
    np.fill_diagonal(weights, 0.)
    denominator = weights.sum(axis=1)
    out = []
    for k in range(135):
        pt = np.load(results / f'steps/step_{k:03d}.npz')['pseudotime'][:1232]
        curve, _ = fixed.fixed_gam(ages, pt, gm, total_edf=6)
        gap = pt - PchipInterpolator(curve.age, curve.fitted_pseudotime)(ages)
        mid = (gap[None, :] < gap[:, None]).astype(float)
        mid += .5 * (gap[None, :] == gap[:, None])
        percentiles = 100 * (weights * mid).sum(axis=1) / denominator
        if k == 0:
            reference_percentiles = percentiles.copy()
        change = np.abs(percentiles - reference_percentiles)
        out.append(dict(added_people=k,
                        GNPC_within_age_percentile_rho_vs_step0=float(spearmanr(reference_percentiles, percentiles).statistic),
                        GNPC_mean_abs_within_age_percentile_change_pp=float(change.mean()),
                        GNPC_p95_abs_within_age_percentile_change_pp=float(np.quantile(change, .95))))
    return metrics.merge(pd.DataFrame(out), on='added_people', validate='one_to_one')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    metadata = pd.read_csv(args.results / 'PRIVATE_GNPC_final.csv')
    metrics = pd.read_csv(args.results / 'incremental_GNPC_metrics.csv')
    assert len(metadata) == 1232
    assert metrics.added_people.tolist() == list(range(135))
    result = age_rank_diagnostics(metadata, metrics, args.results)
    args.output.mkdir(parents=True, exist_ok=False)
    result.to_csv(args.output / 'within_age_reference_sensitivity.csv', index=False)


if __name__ == '__main__':
    main()
