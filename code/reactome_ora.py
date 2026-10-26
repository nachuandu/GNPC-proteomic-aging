"""Optional offline over-representation analysis with an explicit tested universe.

Rewritten utility, NOT a claim to reproduce the historical Enrichr output.
Supply a versioned Reactome GMT and the actual eligible background gene labels.
This script makes no network requests and sends no gene lists to services.
"""
import argparse
import hashlib
import json
import pandas as pd
from scipy.stats import hypergeom
from analysis_core import BINS, bh, output_folder


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--membership', required=True)
    p.add_argument('--background', required=True, help='CSV with gene_label column')
    p.add_argument('--gmt', required=True)
    p.add_argument('--release', required=True)
    p.add_argument('--output', required=True)
    a = p.parse_args()
    from pathlib import Path
    universe = set(pd.read_csv(a.background).gene_label.dropna().astype(str)) - {''}
    if not universe:
        raise ValueError('Empty background.')
    membership = pd.read_csv(a.membership).set_index('gene_label')
    pathways = {}
    for line in Path(a.gmt).read_text().splitlines():
        fields = line.split('\t')
        if len(fields) >= 3:
            if fields[0] in pathways:
                raise ValueError('Duplicate pathway ID in GMT.')
            genes = set(fields[2:]) & universe
            if genes:
                pathways[fields[0]] = genes
    rows = []
    for label in BINS:
        values = membership[label].astype(str).str.lower()
        if not values.isin(['true','false','1','0']).all():
            raise ValueError('Membership values must be boolean or 0/1.')
        selected = set(membership.index[values.isin(['true','1'])])
        if not selected <= universe:
            raise ValueError('Selected genes outside background; review universe first.')
        group = []
        for pathway, genes in pathways.items():
            overlap = len(genes & selected)
            pvalue = hypergeom.sf(overlap-1, len(universe), len(genes), len(selected))
            group.append({'age_bin':label, 'pathway':pathway, 'protein_count':overlap,
                          'pathway_background_n':len(genes), 'selected_n':len(selected),
                          'background_n':len(universe), 'p':float(pvalue)})
        frame = pd.DataFrame(group)
        if frame.empty:
            raise ValueError('No pathways overlap the background.')
        frame['q'] = bh(frame.p)
        rows.append(frame)
    out = output_folder(a.output)
    pd.concat(rows).to_csv(out / 'PRIVATE_Reactome_ORA.csv', index=False)
    (out / 'Reactome_provenance.json').write_text(json.dumps({'release':a.release,
        'gmt_sha256':hashlib.sha256(Path(a.gmt).read_bytes()).hexdigest(),
        'method':'hypergeometric upper-tail; BH over all testable pathways within each bin',
        'historical_Enrichr_equivalence_verified':False}, indent=2))


if __name__ == '__main__':
    main()
