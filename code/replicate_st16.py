"""Directional concordance of GNPC classes in an independently scored ST16 cohort.

Mapping is an explicit approved input, not inferred by dropping numeric suffixes.
All duplicate/many-to-one mappings are classified as ambiguous by this adapter.
"""
import argparse
import numpy as np
import pandas as pd
from analysis_core import output_folder


def concordance(gnpc, st16, mapping):
    for table in [gnpc, st16]:
        if not table.feature.is_unique:
            raise ValueError('Association tables require unique feature IDs.')
    required = ['gnpc_feature', 'st16_feature']
    mapping = mapping[required].drop_duplicates().dropna()
    unique = mapping[~mapping.gnpc_feature.duplicated(False) &
                     ~mapping.st16_feature.duplicated(False)]
    selected = gnpc[gnpc['class'].isin(['Class I', 'Class II'])].copy()
    joined = selected.merge(unique, left_on='feature', right_on='gnpc_feature',
                            how='left', validate='one_to_one')
    other = st16[['feature', 'delta_rho']].rename(columns={'feature': 'st16_feature',
                                                        'delta_rho': 'delta_rho_ST16'})
    joined = joined.merge(other, on='st16_feature', how='left', validate='many_to_one')
    matched = np.isfinite(joined.delta_rho) & np.isfinite(joined.delta_rho_ST16)
    joined['status'] = 'Unmatched or ambiguous'
    product = joined.delta_rho * joined.delta_rho_ST16
    joined.loc[matched & product.gt(0), 'status'] = 'Concordant direction'
    joined.loc[matched & product.lt(0), 'status'] = 'Opposite direction'
    joined.loc[matched & product.eq(0), 'status'] = 'Zero difference'
    rows = []
    for name, frame in joined.groupby('class'):
        counts = frame.status.value_counts()
        denominator = int(frame.status.isin(['Concordant direction', 'Opposite direction']).sum())
        n = int(counts.get('Concordant direction', 0))
        rows.append({'class': name, 'n_total': len(frame), 'n_concordant': n,
                     'n_opposite': int(counts.get('Opposite direction', 0)),
                     'n_unmatched_or_ambiguous': int(counts.get('Unmatched or ambiguous', 0)),
                     'n_zero_difference': int(counts.get('Zero difference', 0)),
                     'n_direction_evaluable': denominator,
                     'percent_concordant': 100*n/denominator if denominator else np.nan})
    return joined, pd.DataFrame(rows)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--gnpc', required=True, help='GNPC global associations')
    p.add_argument('--st16', required=True, help='ST16 global associations')
    p.add_argument('--mapping', required=True)
    p.add_argument('--output', required=True)
    a = p.parse_args()
    detail, summary = concordance(pd.read_csv(a.gnpc), pd.read_csv(a.st16), pd.read_csv(a.mapping))
    out = output_folder(a.output)
    detail.to_csv(out / 'PRIVATE_directional_concordance.csv', index=False)
    summary.to_csv(out / 'directional_concordance_counts.csv', index=False)


if __name__ == '__main__':
    main()
