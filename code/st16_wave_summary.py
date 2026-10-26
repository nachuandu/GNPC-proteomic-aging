"""Summarize ST16 assay R2 using frozen GNPC gene-label Wave assignments.

No new Wave selection is performed in ST16. Every matching assay is retained;
multiple ST16 probes for the same exact label remain separate observations.
"""
import argparse
import pandas as pd
from analysis_core import output_folder


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--associations', required=True)
    p.add_argument('--annotation', required=True)
    p.add_argument('--gnpc-waves', required=True)
    p.add_argument('--output', required=True)
    a = p.parse_args()
    ann = pd.read_csv(a.annotation, keep_default_na=False)[['feature', 'gene_label']]
    waves = pd.read_csv(a.gnpc_waves)[['gene_label', 'wave']]
    if not ann.feature.is_unique or not waves.gene_label.is_unique:
        raise ValueError('Ambiguous annotations or gene-label Wave membership.')
    table = pd.read_csv(a.associations).merge(ann, on='feature', validate='many_to_one')
    table = table.merge(waves[waves.wave.isin(['Wave 1','Wave 2','Wave 3','Wave 4'])],
                        on='gene_label', validate='many_to_one')
    out = output_folder(a.output)
    table.to_csv(out / 'PRIVATE_ST16_frozen_wave_R2.csv', index=False)
    table.groupby(['wave', 'age_bin']).R2_deviation.agg(['count','median','mean']).to_csv(
        out / 'ST16_wave_R2_summary.csv')


if __name__ == '__main__':
    main()
