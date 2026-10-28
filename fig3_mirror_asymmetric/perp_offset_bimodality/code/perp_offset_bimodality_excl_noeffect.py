#!/usr/bin/env python
"""perp_offset_bimodality.py rerun without the no-effect CNVs.

A locus is dropped when neither its deletion nor its duplication effect reaches p < 0.05
(two-sided, from beta/SE). The filter is applied before anything is fitted, so the mirror
line, the offsets, the mixture fits and the bootstrap nulls all use the remaining loci
only. Everything else (inputs, mirror fit, EM, nulls, seeds) is perp_offset_bimodality.py.

Outputs go to outputs/excl_noeffect_p05/ with the same file names as the original run;
make the figure with
    python code/make_figure.py --variant excl_noeffect_p05
"""
from __future__ import annotations
import argparse, pickle
import numpy as np
import pandas as pd
from scipy.stats import norm

from perp_offset_bimodality import HERE, analyse

VARIANT = 'excl_noeffect_p05'
P_CUTOFF = 0.05


def drop_no_effect(dat, name, dropped):
    pD = 2 * norm.sf(np.abs(dat.Deletion_beta / dat.Deletion_SE))
    pU = 2 * norm.sf(np.abs(dat.Duplication_beta / dat.Duplication_SE))
    null = (pD >= P_CUTOFF) & (pU >= P_CUTOFF)
    cols = [c for c in ['CNV_ID', 'Dose_response_class', 'Gene_count', 'CNV_length_Mb',
                        'Deletion_beta', 'Duplication_beta'] if c in dat.columns]
    dropped.append(dat.loc[null, cols].assign(trait=name, p_DEL=pD[null], p_DUP=pU[null]))
    return dat.loc[~null].reset_index(drop=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--nboot', type=int, default=5000)
    ap.add_argument('--starts', type=int, default=25)
    a = ap.parse_args()
    root = HERE.parent
    dropped = []
    H = drop_no_effect(pd.read_csv(root / 'CNV_height_BMI_corrected_directional_models_v2' /
                                   'height' / 'inputs' / 'sacks_height_all_loci_exact.csv'),
                       'Height', dropped)
    M = drop_no_effect(pd.read_csv(root / 'CNV_height_BMI_corrected_directional_models_v2' /
                                   'bmi' / 'inputs' / 'sacks_bmi_complete_reciprocal_loci_exact.csv'),
                       'BMI', dropped)
    rows, per = [], []
    res = [analyse(H, 0.30, 'Height', a.nboot, a.starts, rows, per),
           analyse(M, 11 / 14, 'BMI', a.nboot, a.starts, rows, per),
           analyse(H, 0.30, 'Height', a.nboot, a.starts, rows, per, mirror_subset=True),
           analyse(M, 11 / 14, 'BMI', a.nboot, a.starts, rows, per, mirror_subset=True)]
    out = HERE / 'outputs' / VARIANT
    out.mkdir(parents=True, exist_ok=True)
    tab = pd.DataFrame(rows)
    tab.to_csv(out / 'bimodality_test.csv', index=False)
    pd.concat(per).to_csv(out / 'perp_offsets_per_locus.csv', index=False)
    pd.concat(dropped).to_csv(out / 'excluded_no_effect_loci.csv', index=False)

    def tag_of(r):
        return f'{r["name"]}_{"mirrorclass" if r["mirror_subset"] else "fulldata"}'
    np.savez_compressed(
        out / 'null_LRT_draws.npz',
        **{f'{tag_of(r)}_{t}_{k}': r['fits'][t][k]
           for r in res for t in ('free_var', 'equal_var') for k in ('nullA', 'nullB')})
    (out / '_res.pkl').write_bytes(pickle.dumps(res))
    pd.set_option('display.width', 250)
    print(pd.concat(dropped).to_string())
    print(tab.T.to_string())
