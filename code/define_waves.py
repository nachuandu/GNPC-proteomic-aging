"""Reconstruct exact-label Wave membership, including non-Wave selected patterns."""
import argparse
import pandas as pd
from analysis_core import wave_membership, output_folder


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--associations', required=True)
    p.add_argument('--annotation', required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--r2', type=float, default=.2)
    p.add_argument('--q', type=float, default=.05)
    a = p.parse_args()
    members, counts = wave_membership(pd.read_csv(a.associations), pd.read_csv(a.annotation), a.r2, a.q)
    out = output_folder(a.output)
    members.to_csv(out / 'PRIVATE_gene_label_membership.csv', index=False)
    counts.to_csv(out / 'pattern_counts.csv', index=False)


if __name__ == '__main__':
    main()
