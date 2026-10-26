"""Independent GNPC + one IADRC visit refits; no cumulative inclusion.

Individual-level outputs remain private. Fit-specific 0-1 scores are never
pooled as if they came from a common embedding. Longitudinal summaries use
age-local GNPC-reference percentiles derived separately within each fit.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import importlib.metadata
import importlib.util
import json
import multiprocessing as mp
from pathlib import Path
import time

import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator
from scipy.stats import spearmanr
from threadpoolctl import threadpool_limits

MODEL_SEED = 20261012
CUTS = np.array([100 / 3, 200 / 3])
LABELS = np.array(['Bottom third', 'Middle third', 'Top third'])
STATE = {}
CODE = Path(__file__).resolve().parent


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_json(path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rho(x, y):
    value = float(spearmanr(x, y).statistic)
    assert np.isfinite(value)
    return value


def selected_input(folder):
    manifest = json.loads((folder / 'input_manifest.json').read_text())
    for name, expected in manifest['input_sha256'].items():
        assert sha(folder / name) == expected, name
    Xg = np.load(folder / 'GNPC_X.npy', allow_pickle=False)
    Xi = np.load(folder / 'IADRC_reference_scaled.npy', allow_pickle=False)
    G = pd.read_csv(folder / 'GNPC_reference.csv', dtype={'sample_id': str})
    I = pd.read_csv(folder / 'IADRC_visit_metadata.csv', dtype={'sample_id': str})
    people = pd.read_csv(folder / 'people.csv')
    chosen = people[people.original_240.eq('Yes') &
                    people.followup_status.eq('Later all healthy flags') &
                    people.all_record_status.eq('All recorded flags healthy')]
    history = pd.read_csv(folder / 'visits.csv')
    assert len(chosen) == 134 and chosen.participant.is_unique
    assert history.loc[history.participant.isin(chosen.participant),
                       'healthy_current'].astype(str).str.lower().eq('true').all()
    use = I.participant.isin(chosen.participant).to_numpy()
    assert Xg.shape == (1232, 227) and Xi.shape == (651, 227)
    assert use.sum() == 443 and I.loc[use, 'participant'].nunique() == 134
    assert I.loc[use, 'flag_interpretation'].eq('H').all()
    assert not I.loc[use, ['participant', 'visit']].duplicated().any()
    assert G.sample_id.is_unique and I.sample_id.is_unique
    assert not set(G.sample_id).intersection(I.sample_id)
    assert np.isfinite(Xg).all() and np.isfinite(Xi).all()
    indices = np.flatnonzero(use)
    return manifest, Xg, Xi, G, I, indices


def age_reference(age, score, gm, fixed, weights):
    curve, report = fixed.fixed_gam(age, score, gm, total_edf=6)
    prediction = PchipInterpolator(curve.age, curve.fitted_pseudotime,
                                  extrapolate=False)
    residual = score - prediction(age)
    mid = (residual[None, :] < residual[:, None]).astype(float)
    mid += .5 * (residual[None, :] == residual[:, None])
    percentiles = 100 * np.sum(weights * mid, axis=1) / weights.sum(axis=1)
    assert np.isfinite(percentiles).all()
    return curve, prediction, residual, percentiles, report


def initialize(input_path, helper_path, output_path, threads, fingerprint):
    folder, output = Path(input_path), Path(output_path)
    _, Xg, Xi, G, I, _ = selected_input(folder)
    helper = load_module(Path(helper_path), 'single_visit_original_helper')
    gm = load_module(CODE / 'GAM_methods.py', 'single_visit_gam')
    fixed = load_module(CODE / 'fixed_GAM_methods.py', 'single_visit_fixed')
    pc = load_module(CODE / 'principal_curve_methods.py', 'single_visit_pc')
    reference = np.load(output / 'GNPC_only.npz', allow_pickle=False)
    ages = G.age.to_numpy()
    weights = np.exp(-.5 * ((ages[:, None] - ages[None, :]) / 5.) ** 2)
    np.fill_diagonal(weights, 0.)
    STATE.update(Xg=Xg, Xi=Xi, G=G, I=I, helper=helper, gm=gm,
                 fixed=fixed, pc=pc, reference=reference, ages=ages,
                 weights=weights, output=output, threads=threads,
                 fingerprint=fingerprint)


def fit_case(ordinal, index):
    s = STATE
    target = s['output'] / 'cases' / f'case_{ordinal:03d}.npz'
    audit_path = target.with_suffix('.json')
    if target.exists() and audit_path.exists():
        saved = json.loads(audit_path.read_text())
        assert saved['iadrc_input_index'] == index
        assert saved['protocol_fingerprint'] == s['fingerprint']
        return saved
    X = np.vstack([s['Xg'], s['Xi'][index:index + 1]])
    assert X.shape == (1233, 227)
    started = time.monotonic()
    with threadpool_limits(limits=s['threads']):
        pt, _, _, _, _, audit = s['helper'].fit_once(
            X, s['ages'], s['pc'], s['threads'], MODEL_SEED)
        curve, prediction, residual, gp, _ = age_reference(
            s['ages'], pt[:1232], s['gm'], s['fixed'], s['weights'])
        visit = s['I'].iloc[index]
        query_residual = float(pt[-1] - prediction(float(visit.age)))
        query_pct, effective_n = s['helper'].percentile(
            s['ages'], residual, np.array([float(visit.age)]),
            np.array([query_residual]), h=5.)
    change = np.abs(gp - s['reference']['gnpc_age_percentile'])
    record = dict(case_index=ordinal, iadrc_input_index=index, joint_rows=1233,
                  included_IADRC_visits=1, included_IADRC_people=1,
                  protocol_fingerprint=s['fingerprint'],
                  GNPC_score_rho=rho(s['reference']['gnpc_score'], pt[:1232]),
                  GNPC_within_age_percentile_rho=rho(
                      s['reference']['gnpc_age_percentile'], gp),
                  GNPC_mean_absolute_percentile_shift=float(change.mean()),
                  GNPC_p95_absolute_percentile_shift=float(np.quantile(change, .95)),
                  GNPC_max_absolute_percentile_shift=float(change.max()),
                  query_fit_specific_score=float(pt[-1]),
                  query_percentile=float(query_pct[0]),
                  query_reference_effective_n=float(effective_n[0]),
                  query_fit_specific_residual=query_residual,
                  principal_curve_converged=audit['principal_curve_converged'],
                  principal_curve_iterations=audit['principal_curve_iterations'],
                  warnings=audit['warnings'], seconds=time.monotonic() - started)
    assert 0 <= record['query_percentile'] <= 100
    np.savez_compressed(target, gnpc_score=pt[:1232], gnpc_age_percentile=gp,
                        iadrc_input_index=np.array([index]),
                        query_score=pt[-1:], query_percentile=query_pct)
    write_json(audit_path, record)
    return record


def aggregate(I, indices, records, output):
    metrics = pd.DataFrame(records).sort_values('case_index')
    assert metrics.case_index.tolist() == list(range(443))
    assert metrics.iadrc_input_index.tolist() == indices.tolist()
    V = I.iloc[indices].copy().reset_index(drop=True)
    V['case_index'] = metrics.case_index.to_numpy()
    V['fit_specific_proteomic_aging_score'] = metrics.query_fit_specific_score.to_numpy()
    V['percentile'] = metrics.query_percentile.to_numpy()
    V['reference_effective_n'] = metrics.query_reference_effective_n.to_numpy()
    V = V.sort_values(['participant', 'visit'], kind='stable')
    baseline = V.groupby('participant', sort=False).head(1).copy()
    last = V.groupby('participant', sort=False).tail(1).copy()
    assert len(baseline) == 134 and baseline.visit.eq(1).all()
    baseline['baseline_tertile_group'] = LABELS[np.digitize(baseline.percentile, CUTS)]
    V = V.merge(baseline[['participant', 'baseline_tertile_group']],
                on='participant', validate='many_to_one')
    columns = ['participant', 'visit', 'age', 'percentile', 'case_index']
    pairs = baseline[columns].merge(last[columns], on='participant',
                                  validate='one_to_one', suffixes=('_baseline', '_last'))
    pairs = pairs[pairs.visit_last.gt(pairs.visit_baseline)].copy()
    assert len(pairs) == pairs.participant.nunique() == 133
    counts = np.zeros((3, 3), dtype=int)
    np.add.at(counts, (np.digitize(pairs.percentile_baseline, CUTS),
                      np.digitize(pairs.percentile_last, CUTS)), 1)
    assert counts.sum() == 133
    n = counts.sum(axis=1)
    percentages = np.divide(100 * counts, n[:, None], out=np.zeros((3, 3)),
                            where=n[:, None] > 0)
    metrics.to_csv(output / 'single_visit_metrics.csv', index=False)
    V.to_csv(output / 'PRIVATE_IADRC_visit_scores.csv', index=False)
    pairs.to_csv(output / 'PRIVATE_pairs.csv', index=False)
    baseline.to_csv(output / 'PRIVATE_baseline.csv', index=False)
    worst = metrics.sort_values(['GNPC_within_age_percentile_rho', 'case_index']).iloc[0]
    summary = dict(n_GNPC=1232, n_IADRC_people=134, n_IADRC_visits=443,
                   n_independent_query_refits=443, n_total_fits=444, n_pairs=133,
                   joint_rows_in_every_query_fit=1233,
                   baseline_tertile_counts=np.bincount(
                       np.digitize(baseline.percentile, CUTS), minlength=3).tolist(),
                   paired_baseline_tertile_counts=n.tolist(),
                   transition_counts=counts.tolist(),
                   transition_row_percentages=percentages.tolist(),
                   rank_persistence_rho=rho(pairs.percentile_baseline, pairs.percentile_last),
                   mean_absolute_baseline_last_percentile_change=float(np.mean(
                       np.abs(pairs.percentile_last - pairs.percentile_baseline))),
                   GNPC_score_rho_min=float(metrics.GNPC_score_rho.min()),
                   GNPC_score_rho_median=float(metrics.GNPC_score_rho.median()),
                   GNPC_score_rho_max=float(metrics.GNPC_score_rho.max()),
                   GNPC_age_percentile_rho_min=float(metrics.GNPC_within_age_percentile_rho.min()),
                   GNPC_age_percentile_rho_median=float(metrics.GNPC_within_age_percentile_rho.median()),
                   GNPC_age_percentile_rho_max=float(metrics.GNPC_within_age_percentile_rho.max()),
                   GNPC_mean_shift_median=float(metrics.GNPC_mean_absolute_percentile_shift.median()),
                   GNPC_mean_shift_max=float(metrics.GNPC_mean_absolute_percentile_shift.max()),
                   GNPC_p95_shift_median=float(metrics.GNPC_p95_absolute_percentile_shift.median()),
                   GNPC_p95_shift_max=float(metrics.GNPC_p95_absolute_percentile_shift.max()),
                   worst_concordance_case_index=int(worst.case_index),
                   all_principal_curves_converged=bool(metrics.principal_curve_converged.all()),
                   min_query_reference_effective_n=float(V.reference_effective_n.min()),
                   note='Percentiles are pooled; raw fit-specific scores are not pooled. Visits are clustered within people. Single-visit joint refitting is not frozen-reference prediction.')
    write_json(output / 'analysis_summary.json', summary)
    print('SUMMARY ' + json.dumps(summary), flush=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--input', type=Path, required=True)
    p.add_argument('--helper', type=Path, default=CODE / 'trajectory_backend.py')
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--workers', type=int, default=4)
    p.add_argument('--threads', type=int, default=2)
    p.add_argument('--limit', type=int, default=443)
    args = p.parse_args()
    assert 1 <= args.limit <= 443
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / 'cases').mkdir(exist_ok=True)
    manifest, Xg, Xi, G, I, indices = selected_input(args.input)
    protocol = dict(mode='Independent single-visit joint refits; no cumulative inclusion',
                    n_GNPC=1232, query_visits_per_fit=1, query_fit_rows=1233,
                    selected_IADRC_people=134, selected_IADRC_visits=443, n_analytes=227,
                    input_sha256=manifest['input_sha256'], script_sha256=sha(Path(__file__)),
                    original_helper_sha256=sha(args.helper), model_seed=MODEL_SEED,
                    helper_sha256={name:sha(CODE / name) for name in ['GAM_methods.py', 'fixed_GAM_methods.py', 'principal_curve_methods.py']},
                    PHATE=dict(knn=50, t=50, decay=40, gamma=1, n_pca=100,
                               n_landmark=2000, mds='metric', mds_solver='smacof',
                               knn_dist='euclidean', mds_dist='euclidean'),
                    principal_curve='Existing six-knot natural-spline curve; native geometry; tolerance 1e-6; max_iter 1000',
                    preprocessing='Fixed GNPC-selected 227 features and original GNPC standardization; no refitting of feature selection or scaling',
                    age_reference='GNPC-only age curve within each refit; unconstrained natural cubic spline, total EDF 6',
                    percentile='Gaussian age weights, bandwidth 5 years; weighted empirical mid-distribution',
                    GNPC_diagnostic='Each GNPC query excluded from its own percentile reference',
                    baseline_groups='GNPC-reference cutoffs 100/3 and 200/3; not equally sized IADRC groups',
                    no_frozen_reference_validation_claim=True, no_aging_rate_claim=True,
                    pool_fit_specific_scores=False, inferential_unit='Participant, not visit',
                    workers=args.workers, threads_per_worker=args.threads,
                    versions={n: importlib.metadata.version(n) for n in
                              ['numpy', 'scipy', 'pandas', 'phate', 'scikit-learn']})
    fingerprint = hashlib.sha256(json.dumps(protocol, sort_keys=True).encode()).hexdigest()
    protocol['fingerprint'] = fingerprint
    prior = args.output / 'protocol.json'
    if prior.exists():
        assert json.loads(prior.read_text())['fingerprint'] == fingerprint
    else:
        write_json(prior, protocol)
    helper = load_module(args.helper, 'single_visit_baseline_helper')
    pc = load_module(CODE / 'principal_curve_methods.py', 'single_visit_base_pc')
    gm = load_module(CODE / 'GAM_methods.py', 'single_visit_base_gam')
    fixed = load_module(CODE / 'fixed_GAM_methods.py', 'single_visit_base_fixed')
    started = time.monotonic()
    reference_path = args.output / 'GNPC_only.npz'
    if not reference_path.exists():
        with threadpool_limits(limits=args.threads):
            score, xy, curve, curve_pt, _, audit = helper.fit_once(
                Xg, G.age.to_numpy(), pc, args.threads, MODEL_SEED)
            ages = G.age.to_numpy()
            weights = np.exp(-.5 * ((ages[:, None] - ages[None, :]) / 5.) ** 2)
            np.fill_diagonal(weights, 0.)
            age_curve, _, _, gp, gam_report = age_reference(ages, score, gm, fixed, weights)
        np.savez_compressed(reference_path, gnpc_score=score, gnpc_age_percentile=gp,
                            coordinates=xy, principal_curve=curve, curve_score=curve_pt)
        age_curve.to_csv(args.output / 'GNPC_only_age_curve.csv', index=False)
        G = G.copy()
        G['proteomic_aging_score'] = score
        G['within_age_percentile_LOO'] = gp
        G.to_csv(args.output / 'PRIVATE_GNPC_reference.csv', index=False)
        write_json(args.output / 'GNPC_only_fit_audit.json', audit)
        write_json(args.output / 'GNPC_only_age_model.json', gam_report)
        print('GNPC_ONLY_READY ' + json.dumps(audit), flush=True)
    records, failures = [], []
    with ProcessPoolExecutor(max_workers=args.workers, mp_context=mp.get_context('spawn'),
            initializer=initialize,
            initargs=(str(args.input), str(args.helper), str(args.output), args.threads, fingerprint)) as executor:
        futures = {executor.submit(fit_case, k, int(index)): k
                   for k, index in enumerate(indices[:args.limit])}
        for future in as_completed(futures):
            ordinal = futures[future]
            try:
                record = future.result()
                records.append(record)
            except Exception as exc:
                failures.append(dict(case_index=ordinal, error=repr(exc)))
            if (len(records) + len(failures)) % 10 == 0 or len(records) + len(failures) == args.limit:
                write_json(args.output / 'status.json', dict(status='RUNNING',
                    completed_visits=len(records), failed_visits=failures,
                    requested_visits=args.limit, seconds=time.monotonic() - started))
                print(json.dumps(dict(completed_visits=len(records),
                    failed_visits=len(failures), seconds=time.monotonic() - started)), flush=True)
    write_json(args.output / 'failures.json', failures)
    assert not failures, failures
    pd.DataFrame(records).sort_values('case_index').to_csv(
        args.output / 'single_visit_metrics.csv', index=False)
    if args.limit == 443:
        aggregate(I, indices, records, args.output)
    for name, expected in manifest['input_sha256'].items():
        assert sha(args.input / name) == expected, name
    write_json(args.output / 'status.json', dict(
        status='COMPLETE' if args.limit == 443 else 'SMOKE_TEST_COMPLETE',
        completed_visits=len(records), requested_visits=args.limit,
        all_query_fits_have_exactly_one_visit=True, original_inputs_unchanged=True,
        seconds=time.monotonic() - started))


if __name__ == '__main__':
    main()
