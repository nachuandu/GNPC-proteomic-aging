"""Cumulative participant-addition audit with fixed GNPC-selected features.

Each fit starts from scratch. Each participant contributes all eligible visits.
This is a joint-fit sensitivity analysis, not frozen-reference validation.
"""
import argparse
import contextlib
import hashlib
import importlib.metadata
import importlib.util
import io
import json
import pickle
import platform
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import phate
from scipy.interpolate import PchipInterpolator
from scipy.stats import rankdata, spearmanr
from threadpoolctl import threadpool_limits


MODEL_SEED = 20261012
ORDER_SEED = 20261018
LABELS = np.array(['Bottom third', 'Middle third', 'Top third'])
CUTS = np.array([100 / 3, 200 / 3])


def save_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


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


def derive_tables(G, I, pt, gm, fixed):
    G, I = G.copy(), I.copy()
    ng = len(G)
    G['pseudotime'] = pt[:ng]
    I['pseudotime'] = pt[ng:]
    curve, gam_report = fixed.fixed_gam(G.age.to_numpy(), pt[:ng], gm, total_edf=6)
    predict = PchipInterpolator(curve.age, curve.fitted_pseudotime, extrapolate=False)
    for frame in [G, I]:
        frame['age_expected_ptime'] = predict(frame.age)
        frame['gap'] = frame.pseudotime - frame.age_expected_ptime
        assert np.isfinite(frame[['pseudotime', 'age_expected_ptime', 'gap']]).all().all()
    I['percentile'], I['reference_effective_n'] = percentile(
        G.age.to_numpy(), G.gap.to_numpy(), I.age.to_numpy(), I.gap.to_numpy())
    ordered = I.sort_values(['participant', 'visit'], kind='stable')
    baseline = ordered.groupby('participant', sort=False).head(1).copy()
    last = ordered.groupby('participant', sort=False).tail(1).copy()
    assert baseline.visit.eq(1).all()
    membership = baseline[['participant']].copy()
    membership['baseline_tertile_group'] = LABELS[np.digitize(baseline.percentile, CUTS)]
    I = ordered.merge(membership, on='participant', validate='many_to_one')
    columns = ['participant', 'visit', 'age', 'pseudotime', 'gap', 'percentile']
    pairs = baseline[columns].merge(last[columns], on='participant', validate='one_to_one',
                                   suffixes=('_baseline', '_last'))
    pairs = pairs[pairs.visit_last.gt(pairs.visit_baseline)].copy()
    counts = np.zeros((3, 3), dtype=int)
    np.add.at(counts, (np.digitize(pairs.percentile_baseline, CUTS),
                      np.digitize(pairs.percentile_last, CUTS)), 1)
    denominator = counts.sum(axis=1)
    percentages = np.divide(100. * counts, denominator[:, None],
                            out=np.zeros((3, 3)), where=denominator[:, None] > 0)
    baseline_n = np.bincount(np.digitize(baseline.percentile, CUTS), minlength=3)
    report = dict(n_people=int(len(baseline)), n_visits=int(len(I)), n_pairs=int(len(pairs)),
                  baseline_tertile_counts=baseline_n.tolist(),
                  paired_baseline_tertile_counts=denominator.tolist(),
                  transition_counts=counts.tolist(), transition_row_percentages=percentages.tolist(),
                  rank_persistence_rho=correlation(pairs.percentile_baseline, pairs.percentile_last)
                  if len(pairs) >= 3 else None,
                  percentile_age_bandwidth_years=5., GNPC_GAM_total_EDF=6,
                  minimum_reference_effective_n=float(I.reference_effective_n.min()) if len(I) else None)
    return G, I, curve, pairs, report, gam_report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--last-step', type=int, default=134)
    parser.add_argument('--seed-controls', action='store_true')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    steps = args.output / 'steps'
    steps.mkdir(exist_ok=True)
    manifest = json.loads((args.input / 'input_manifest.json').read_text())
    for name, expected in manifest['input_sha256'].items():
        assert hashlib.sha256((args.input / name).read_bytes()).hexdigest() == expected, name
    Xg = np.load(args.input / 'GNPC_X.npy', allow_pickle=False)
    Xi = np.load(args.input / 'IADRC_reference_scaled.npy', allow_pickle=False)
    G = pd.read_csv(args.input / 'GNPC_reference.csv', dtype={'sample_id': str})
    I = pd.read_csv(args.input / 'IADRC_visit_metadata.csv', dtype={'sample_id': str})
    people = pd.read_csv(args.input / 'people.csv')
    chosen = people[people.original_240.eq('Yes') &
                    people.followup_status.eq('Later all healthy flags') &
                    people.all_record_status.eq('All recorded flags healthy')].copy()
    history = pd.read_csv(args.input / 'visits.csv')
    assert len(chosen) == 134 and chosen.participant.is_unique
    assert history.loc[history.participant.isin(chosen.participant),
                       'healthy_current'].astype(str).str.lower().eq('true').all()
    mask = I.participant.isin(chosen.participant).to_numpy()
    assert Xg.shape == (1232, 227) and Xi.shape == (651, 227)
    assert mask.sum() == 443 and I.loc[mask, 'participant'].nunique() == 134
    assert I.loc[mask, 'flag_interpretation'].eq('H').all()
    assert G.sample_id.is_unique and I.sample_id.is_unique
    assert not set(G.sample_id).intersection(I.sample_id)
    assert not I.loc[mask, ['participant', 'visit']].duplicated().any()
    assert np.isfinite(Xg).all() and np.isfinite(Xi).all()
    ordered_people = np.sort(chosen.participant.to_numpy())
    order = np.random.default_rng(ORDER_SEED).permutation(ordered_people)
    selected_rank = {p: k + 1 for k, p in enumerate(order)}
    order_table = pd.DataFrame({'addition_step': np.arange(1, 135), 'participant': order})
    order_table['eligible_visits'] = order_table.participant.map(I.groupby('participant').size())
    order_table.to_csv(args.output / 'PRIVATE_participant_addition_order.csv', index=False)
    pc = load_module(Path(__file__).resolve().parent / 'principal_curve_methods.py', 'curve_methods')
    gm = load_module(Path(__file__).resolve().parent / 'GAM_methods.py', 'gam_methods')
    fixed = load_module(Path(__file__).resolve().parent / 'fixed_GAM_methods.py', 'fixed_methods')
    protocol = dict(n_GNPC=1232, selected_IADRC_people=134, selected_IADRC_visits=443,
                    n_analytes=227, preprocessing='Original GNPC log2-imputed mean and sample SD; no new scaling or feature selection',
                    model_seed=MODEL_SEED, participant_order_seed=ORDER_SEED,
                    participant_order='One prespecified random permutation; not selected using scores, ages or outcomes',
                    input_row_order='Original GNPC rows followed by included IADRC rows in original metadata order',
                    fit_at_each_step='Fresh PHATE and principal-curve fit to GNPC plus all visits of cumulative participants',
                    PHATE=operator(args.threads, MODEL_SEED).get_params(),
                    principal_curve='Existing custom six-knot single curve; native geometry; max_iter=1000; tolerance=1e-6',
                    age_use='Age-informed fixed features and orientation using GNPC age; no age input to PHATE or curve',
                    no_independent_validation=True, no_aging_rate_inference=True,
                    repeated_visits_enter_joint_graph=True, last_step_equivalent_to_one_shot_with_identical_set_order_seed=True,
                    versions={n: importlib.metadata.version(n) for n in ['phate', 'numpy', 'scipy', 'scikit-learn', 'pandas']},
                    python=platform.python_version(), input_sha256=manifest['input_sha256'])
    save_json(args.output / 'protocol.json', protocol)
    metrics = []
    started = time.monotonic()
    with threadpool_limits(limits=args.threads):
        for k in range(args.last_step + 1):
            use = I.participant.map(selected_rank).le(k).to_numpy()
            indices = np.flatnonzero(use)
            target = steps / f'step_{k:03d}.npz'
            audit_path = steps / f'step_{k:03d}.json'
            if target.exists() and audit_path.exists():
                result = np.load(target, allow_pickle=False)
                np.testing.assert_array_equal(result['iadrc_indices'], indices)
                pt = result['pseudotime']
                audit = json.loads(audit_path.read_text())
            else:
                X = np.vstack([Xg, Xi[indices]])
                pt, xy, curve, curve_pt, op, audit = fit_once(
                    X, G.age.to_numpy(), pc, args.threads, MODEL_SEED)
                np.savez_compressed(target, pseudotime=pt, coordinates=xy, curve=curve,
                                    curve_ptime=curve_pt, iadrc_indices=indices)
                save_json(audit_path, audit)
                if k in [0, 134]:
                    with (args.output / f'PRIVATE_PHATE_step_{k:03d}.pkl').open('wb') as stream:
                        pickle.dump(op, stream)
            if k == 0:
                baseline_pt = pt[:1232].copy()
                baseline_rank = rankdata(baseline_pt) / 1232 * 100
            gnpc_pt = pt[:1232]
            rank_change = np.abs(rankdata(gnpc_pt) / 1232 * 100 - baseline_rank)
            row = dict(added_people=k, added_visits=len(indices), joint_rows=len(pt),
                       GNPC_vs_step0_rho=correlation(baseline_pt, gnpc_pt),
                       GNPC_vs_saved_historical_rho=correlation(G.pseudotime, gnpc_pt),
                       GNPC_age_rho=correlation(G.age, gnpc_pt),
                       GNPC_mean_abs_rank_change_pp=float(rank_change.mean()),
                       GNPC_p95_abs_rank_change_pp=float(np.quantile(rank_change, .95)),
                       GNPC_max_abs_rank_change_pp=float(rank_change.max()),
                       GNPC_same_tertile_fraction=float(np.mean(
                           np.digitize(baseline_rank, CUTS) == np.digitize(rankdata(gnpc_pt) / 1232 * 100, CUTS))),
                       fit_seconds=audit['seconds'])
            metrics.append(row)
            pd.DataFrame(metrics).to_csv(args.output / 'incremental_GNPC_metrics.csv', index=False)
            save_json(args.output / 'status.json', dict(status='RUNNING', completed_step=k,
                                                       requested_last_step=args.last_step,
                                                       elapsed_seconds=time.monotonic() - started))
            print(json.dumps(row), flush=True)
        if args.last_step == 134:
            Gf, If, age_curve, pairs, summary, gam_report = derive_tables(G, I.iloc[indices], pt, gm, fixed)
            assert (len(If), If.participant.nunique(), len(pairs)) == (443, 134, 133)
            Gf.to_csv(args.output / 'PRIVATE_GNPC_final.csv', index=False)
            If.to_csv(args.output / 'PRIVATE_IADRC_final.csv', index=False)
            pairs.to_csv(args.output / 'PRIVATE_pairs_final.csv', index=False)
            age_curve.to_csv(args.output / 'GNPC_final_age_curve.csv', index=False)
            save_json(args.output / 'final_GAM_report.json', gam_report)
            summary.update(final_GNPC_vs_new_GNPC_only_rho=metrics[-1]['GNPC_vs_step0_rho'],
                           final_GNPC_vs_saved_historical_rho=metrics[-1]['GNPC_vs_saved_historical_rho'],
                           minimum_incremental_GNPC_rho=min(r['GNPC_vs_step0_rho'] for r in metrics),
                           final_GNPC_mean_abs_rank_change_pp=metrics[-1]['GNPC_mean_abs_rank_change_pp'],
                           final_GNPC_p95_abs_rank_change_pp=metrics[-1]['GNPC_p95_abs_rank_change_pp'],
                           no_participants_outside_134_entered_any_refit=True,
                           original_inputs_unchanged=True, total_seconds=time.monotonic() - started)
            save_json(args.output / 'analysis_summary.json', summary)
            if args.seed_controls:
                controls = []
                for seed in [20261013, 20261014]:
                    for endpoint in [0, 134]:
                        X = Xg if endpoint == 0 else np.vstack([Xg, Xi[indices]])
                        control_pt, _, _, _, _, audit = fit_once(X, G.age.to_numpy(), pc, args.threads, seed)
                        reference = baseline_pt if endpoint == 0 else pt[:1232]
                        controls.append(dict(seed=seed, added_people=endpoint,
                                             GNPC_rho_vs_primary_same_endpoint=correlation(reference, control_pt[:1232]),
                                             fit_seconds=audit['seconds']))
                        np.savez_compressed(steps / f'seed_{seed}_endpoint_{endpoint:03d}.npz', pseudotime=control_pt)
                pd.DataFrame(controls).to_csv(args.output / 'algorithmic_seed_controls.csv', index=False)
        for name, expected in manifest['input_sha256'].items():
            assert hashlib.sha256((args.input / name).read_bytes()).hexdigest() == expected
    save_json(args.output / 'status.json', dict(status='COMPLETE', completed_step=args.last_step,
                                               elapsed_seconds=time.monotonic() - started))


if __name__ == '__main__':
    main()
