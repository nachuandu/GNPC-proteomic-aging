"""Fit a monotonic age GAM and compute score minus expected score.

Spline count and penalty are REQUIRED: historical GNPC/ST16 experiments used
different settings. Supply the parameters corresponding to the saved analysis.
This is not the unconstrained EDF-6 age reference used for IADRC.
"""
import argparse
import json
import numpy as np
import pandas as pd
from analysis_core import output_folder


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scores', required=True)
    p.add_argument('--n-splines', required=True, type=int)
    p.add_argument('--lam', required=True, type=float)
    p.add_argument('--output', required=True)
    a = p.parse_args()
    from pygam import LinearGAM, s
    df = pd.read_csv(a.scores, dtype={'sample_id': str})
    if not df.sample_id.is_unique or df.sample_id.isna().any():
        raise ValueError('Unique sample_id required.')
    if not np.isfinite(df[['age', 'proteomic_aging_score']].to_numpy()).all():
        raise ValueError('Finite age and positively oriented score required.')
    if df.age.corr(df.proteomic_aging_score, method='spearman') <= 0:
        raise ValueError('Orient scores positively before deriving deviations.')
    if a.n_splines < 4 or a.lam < 0:
        raise ValueError('Invalid spline count or penalty.')
    model = LinearGAM(s(0, n_splines=a.n_splines, spline_order=3,
                        constraints='monotonic_inc'), lam=a.lam).fit(
                            df[['age']].to_numpy(), df.proteomic_aging_score.to_numpy())
    df['age_expected_score'] = model.predict(df[['age']].to_numpy())
    df['proteomic_deviation_score'] = df.proteomic_aging_score - df.age_expected_score
    out = output_folder(a.output)
    df.to_csv(out / 'PRIVATE_scores_and_deviations.csv', index=False)
    (out / 'age_model.json').write_text(json.dumps({'n_splines': a.n_splines,
        'lam': a.lam, 'spline_order': 3, 'monotonic': True,
        'deviation': 'score minus fitted score; no residual standardization'}, indent=2))


if __name__ == '__main__':
    main()
