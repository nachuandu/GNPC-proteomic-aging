"""Build the restricted input bundle for run_single_visit.py using explicit IDs.

Eligibility and clinical flags must be supplied from the approved clinical
review. No diagnosis codes, missing flags or approximate sample matches are
automatically converted into healthy status. Outputs must remain private.
"""
import argparse
import hashlib
import json
import shutil
import numpy as np
import pandas as pd
from analysis_core import read_matrix, align_metadata, output_folder


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['gnpc', 'gnpc-metadata', 'iadrc-rfu', 'features', 'scaling',
                 'visit-metadata', 'people', 'history', 'output']:
        p.add_argument('--' + name, required=True)
    a = p.parse_args()
    features = pd.read_csv(a.features).feature
    if len(features) != 227 or not features.is_unique:
        raise ValueError('Expected the frozen 227-feature list.')
    G = read_matrix(a.gnpc).loc[:, features]
    gm = align_metadata(G, a.gnpc_metadata, ['age'])
    I = read_matrix(a.iadrc_rfu).loc[:, features]
    im = align_metadata(I, a.visit_metadata, ['age', 'visit'])
    scaling = pd.read_csv(a.scaling).set_index('feature')
    if not scaling.index.is_unique:
        raise ValueError('Duplicate scaling parameters.')
    scaling = scaling.loc[features, ['mean_log2', 'sd_log2']]
    if not np.isfinite(scaling.to_numpy()).all() or (scaling.sd_log2 <= 0).any():
        raise ValueError('Invalid fixed GNPC scaling parameters.')
    if G.shape != (1232, 227) or I.shape != (651, 227) or (I.to_numpy() <= 0).any():
        raise ValueError('Expected 1232 GNPC and 651 eligible positive-RFU IADRC rows.')
    # Preserve the completed IADRC analysis: log2(RFU), without an added offset.
    Xi = (np.log2(I.to_numpy()) - scaling.mean_log2.to_numpy()) / scaling.sd_log2.to_numpy()
    out = output_folder(a.output)
    np.save(out / 'GNPC_X.npy', G.to_numpy())
    np.save(out / 'IADRC_reference_scaled.npy', Xi)
    gm.reset_index().to_csv(out / 'GNPC_reference.csv', index=False)
    im.reset_index().to_csv(out / 'IADRC_visit_metadata.csv', index=False)
    shutil.copyfile(a.people, out / 'people.csv')
    shutil.copyfile(a.history, out / 'visits.csv')
    names = ['GNPC_X.npy', 'IADRC_reference_scaled.npy', 'GNPC_reference.csv',
             'IADRC_visit_metadata.csv', 'people.csv', 'visits.csv']
    manifest = {name: hashlib.sha256((out / name).read_bytes()).hexdigest() for name in names}
    (out / 'input_manifest.json').write_text(json.dumps({'input_sha256': manifest}, indent=2))
    # Reuse the final runner's exact inclusion and count checks before declaring success.
    from run_single_visit import selected_input
    selected_input(out)
    print('Input checks passed: 134 people and 443 single-visit queries.')


if __name__ == '__main__':
    main()
