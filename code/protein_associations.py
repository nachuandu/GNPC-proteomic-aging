"""Global Spearman associations or within-age GAM R2 + Pearson/BH screening."""
import argparse
from analysis_core import (read_matrix, align_metadata, global_associations,
                           within_age_associations, output_folder)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode', choices=['global', 'within-age'])
    p.add_argument('--matrix', required=True)
    p.add_argument('--metadata', required=True)
    p.add_argument('--output', required=True)
    a = p.parse_args()
    X = read_matrix(a.matrix)
    columns = ['age', 'proteomic_aging_score'] if a.mode == 'global' else ['age', 'proteomic_deviation_score']
    meta = align_metadata(X, a.metadata, columns)
    result = (global_associations if a.mode == 'global' else within_age_associations)(X, meta)
    result.to_csv(output_folder(a.output) / 'PRIVATE_associations.csv', index=False)


if __name__ == '__main__':
    main()
