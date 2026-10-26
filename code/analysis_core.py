"""Shared numerical routines. No study records or fitted parameters are bundled."""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

BINS = ('50_60', '60_70', '70_80', '80_90')
WAVES = {'1111': 'Wave 1', '0111': 'Wave 2', '0011': 'Wave 3', '0001': 'Wave 4'}


def bh(p):
    """BH over the complete supplied family; undefined tests count as P=1."""
    p = np.asarray(p, dtype=float).copy()
    p[~np.isfinite(p)] = 1.
    if np.any((p < 0) | (p > 1)):
        raise ValueError('P values must lie in [0,1].')
    if not len(p):
        return p
    order = np.argsort(p, kind='stable')
    q = np.empty_like(p)
    q[order] = np.minimum(1, np.minimum.accumulate(
        (p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1])
    return q


def read_matrix(path):
    table = pd.read_csv(path, dtype={'sample_id': str})
    if table.columns[0] != 'sample_id':
        raise ValueError('First column must be sample_id; remaining columns are assays.')
    if table.sample_id.isna().any() or not table.sample_id.is_unique:
        raise ValueError('Missing or duplicate sample_id.')
    matrix = table.set_index('sample_id').astype(float)
    if not matrix.columns.is_unique or not np.isfinite(matrix.to_numpy()).all():
        raise ValueError('Matrix requires unique assays and finite processed values.')
    return matrix


def align_metadata(matrix, path, required):
    meta = pd.read_csv(path, dtype={'sample_id': str})
    if not meta.sample_id.is_unique or meta.sample_id.isna().any():
        raise ValueError('Metadata sample_id must be unique and nonmissing.')
    meta = meta.set_index('sample_id').loc[matrix.index]
    if not np.isfinite(meta[list(required)].to_numpy(float)).all():
        raise ValueError('Required metadata must be finite.')
    return meta


def correlation(x, y, method='pearson'):
    x, y = np.asarray(x, float), np.asarray(y, float)
    good = np.isfinite(x) & np.isfinite(y)
    x, y = x[good], y[good]
    if len(x) < 3 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return np.nan, 1., len(x)
    r, p = (spearmanr if method == 'spearman' else pearsonr)(x, y)
    return float(r), float(p), len(x)


def global_associations(matrix, meta):
    rows = []
    for feature in matrix:
        y = matrix[feature].to_numpy()
        ra, pa, n = correlation(meta.age, y, 'spearman')
        rs, ps, _ = correlation(meta.proteomic_aging_score, y, 'spearman')
        rows.append((feature, n, ra, pa, rs, ps))
    out = pd.DataFrame(rows, columns=['feature', 'n', 'rho_age', 'p_age', 'rho_score', 'p_score'])
    out['q_age'], out['q_score'] = bh(out.p_age), bh(out.p_score)
    out['delta_rho'] = out.rho_score - out.rho_age
    out['class'] = classify(out.rho_age, out.delta_rho)
    return out


def classify(rho_age, delta_rho):
    """Signed rho difference, then absolute cutoff; not difference of |rho|."""
    age, delta = np.asarray(rho_age), np.asarray(delta_rho)
    result = np.full(len(age), 'Unclassified', dtype=object)
    result[(np.abs(age) <= .2) & (np.abs(delta) > .1)] = 'Class I'
    result[(np.abs(age) >= .3) & (np.abs(delta) > .1)] = 'Class II'
    return result


def gam_r2(x, y):
    """Historical Gaussian GAM explained-deviance R2, NOT its significance test."""
    from pygam import LinearGAM, s
    x, y = np.asarray(x, float), np.asarray(y, float)
    good = np.isfinite(x) & np.isfinite(y)
    x, y = x[good], y[good]
    if len(x) < 5 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return np.nan
    model = LinearGAM(s(0, n_splines=5, spline_order=3), lam=.6).fit(x[:, None], y)
    return float(model.statistics_['pseudo_r2']['explained_deviance'])


def within_age_associations(matrix, meta):
    tables = []
    for label in BINS:
        lower, upper = map(int, label.split('_'))
        use = (meta.age >= lower) & (meta.age < upper)
        if use.sum() < 5:
            raise ValueError(f'Insufficient observations in {label}; do not silently drop a bin.')
        x = meta.loc[use, 'proteomic_deviation_score'].to_numpy()
        ages = meta.loc[use, 'age'].to_numpy()
        rows = []
        for feature in matrix:
            y = matrix.loc[use, feature].to_numpy()
            r, p, n = correlation(x, y)
            rows.append((feature, label, n, gam_r2(x, y), r, p, gam_r2(ages, y)))
        out = pd.DataFrame(rows, columns=['feature', 'age_bin', 'n', 'R2_deviation',
                                         'pearson_r', 'p_pearson', 'R2_age'])
        out['q_pearson'] = bh(out.p_pearson)
        out['selected'] = (out.R2_deviation >= .2) & (out.q_pearson < .05)
        tables.append(out)
    return pd.concat(tables, ignore_index=True)


def wave_membership(associations, annotations, r2=.2, q=.05):
    """Screen assays, OR within exact nonempty gene labels, retain all patterns."""
    ann = annotations[['feature', 'gene_label']].copy()
    ann['gene_label'] = ann.gene_label.fillna('').astype(str).str.strip()
    if not ann.feature.is_unique:
        raise ValueError('Exactly one annotation label per assay is required.')
    if associations.duplicated(['feature', 'age_bin']).any():
        raise ValueError('Duplicate assay/bin entries.')
    table = associations.copy()
    if set(table.age_bin) != set(BINS):
        raise ValueError('Expected all four age bins.')
    expected = pd.MultiIndex.from_product([table.feature.unique(), BINS])
    if len(table) != len(expected):
        raise ValueError('Every assay must have a row for every bin.')
    table = table.merge(ann, on='feature', how='left', validate='many_to_one')
    if table.gene_label.isna().any():
        raise ValueError('Unannotated assays; supply empty labels explicitly if necessary.')
    table['pass'] = (table.R2_deviation >= r2) & (table.q_pearson < q)
    grouped = table[table.gene_label.ne('')].groupby(['gene_label', 'age_bin'])['pass'].any()
    membership = grouped.unstack('age_bin').reindex(columns=BINS).astype(bool)
    membership['pattern'] = membership.apply(lambda row: ''.join('1' if v else '0' for v in row), axis=1)
    membership['wave'] = membership.pattern.map(WAVES).fillna('Other selected pattern')
    membership.loc[membership.pattern.eq('0000'), 'wave'] = 'Not selected'
    counts = membership.groupby(['pattern', 'wave']).size().rename('n_gene_labels').reset_index()
    return membership.reset_index(), counts


def output_folder(path):
    path = Path(path)
    path.mkdir(parents=True, exist_ok=False)
    return path
