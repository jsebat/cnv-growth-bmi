#!/usr/bin/env python
"""Measurement-error-aware bimodality test of the height perpendicular offsets.

Same offsets and the same separately fitted mirror line as perp_offset_bimodality.py
(line calibrated on the mirror-class loci only), but each locus's own sampling SE is
part of both fits:
  one group   x_i ~ N(mu, tau^2 + s_i^2)
  two groups  x_i ~ pi N(m1, tau^2 + s_i^2) + (1 - pi) N(m2, tau^2 + s_i^2)
tau is the true scatter of loci about their line; s_i is the locus's offset SE. Noisy
loci count for less, so their scatter is not mistaken for a second group.

p is a parametric bootstrap: offsets are redrawn from the fitted one-group model, each
locus with its own SE, and refitted the same way.

Two versions, both written to outputs/me/<version>/ in the layout make_height_figure.py
reads (perp_offsets_per_locus.csv, bimodality_test.csv, _res.pkl):
  all_loci          the 33 loci, offsets from outputs/
  excl_noeffect_p05 the 24 loci left after dropping those where neither DEL nor DUP
                    reaches p < 0.05, offsets from outputs/excl_noeffect_p05/
Run from the folder above code/:  python code/me_bimodality.py
"""
from __future__ import annotations
import argparse, pickle
import numpy as np
import pandas as pd

from perp_offset_bimodality import HERE, SEED, me_ll1, me_ll2, me_posterior_low

VERSIONS = {'all_loci': HERE / 'outputs',
            'excl_noeffect_p05': HERE / 'outputs' / 'excl_noeffect_p05'}

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--nboot', type=int, default=5000)
    ap.add_argument('--starts', type=int, default=25)
    a = ap.parse_args()
    for version, src in VERSIONS.items():
        po = pd.read_csv(src / 'perp_offsets_per_locus.csv')
        po = po[po.trait == 'Height']
        po = po.tail(len(po) // 2).reset_index(drop=True)      # mirror-class calibration
        t = pd.read_csv(src / 'bimodality_test.csv')
        rho = float(t[(t.trait == 'Height') & (t.mirror_fit == 'mirror-class only')].rho_full_fit.iloc[0])
        x = po.perp_offset.to_numpy(float); s = po.perp_offset_SE.to_numpy(float); n = len(x)

        a1, mu, tau1 = me_ll1(x[None, :], s)
        b2, par = me_ll2(x[None, :], s, n_starts=40, seed=SEED)
        obs = float(2 * (b2 - a1)[0]); mu, tau1 = float(mu[0]), float(tau1[0])
        rg = np.random.default_rng(SEED + 7)
        null_x = rg.normal(mu, np.sqrt(tau1 ** 2 + s ** 2), size=(a.nboot, n))
        na, _, _ = me_ll1(null_x, s)
        nb, _ = me_ll2(null_x, s, n_starts=a.starts, seed=SEED + 2)
        null = 2 * (nb - na)
        p = (1 + int((null >= obs).sum())) / (a.nboot + 1)

        pi_, m1, m2, tau2 = par[0]
        if m1 > m2: pi_, m1, m2 = 1 - pi_, m2, m1
        post = me_posterior_low(x[None, :], s, par)[0]
        cls = po.Dose_response_class.str.lower().to_numpy()
        out = HERE / 'outputs' / 'me' / version
        out.mkdir(parents=True, exist_ok=True)
        po.assign(P_lower_group=post).to_csv(out / 'perp_offsets_per_locus.csv', index=False)
        row = dict(trait='Height', mirror_fit='mirror-class only', alternative='equal_var',
                   test='measurement-error-aware', n_loci=n, rho_full_fit=rho,
                   LRT=obs, p_bootstrap_nullB=p, one_group_mean=mu, one_group_tau=tau1,
                   minor_weight=float(min(pi_, 1 - pi_)), comp1_mean=m1, comp2_mean=m2,
                   tau=tau2, separation=m2 - m1,
                   classes_agree=int(((post > .5) == (cls == 'asymmetric')).sum()),
                   null_median_LRT=float(np.median(null)),
                   null_p95_LRT=float(np.quantile(null, .95)), n_boot=a.nboot,
                   n_starts=a.starts)
        pd.DataFrame([row]).to_csv(out / 'bimodality_test.csv', index=False)
        # the layout make_height_figure.py expects; 'me' carries what panel A needs
        res = [dict(name='Height', mirror_subset=True, x=x, se=s, rho=rho, mu=mu, sd=tau1,
                    cls=cls, me=dict(S=s, one=(mu, tau1), two=(pi_, m1, m2, tau2)),
                    fits={'equal_var': dict(obs=obs, nullB=null, nullA=null,
                                            par=(pi_, m1, tau2, m2, tau2), pB=p, pA=p)})]
        (out / '_res.pkl').write_bytes(pickle.dumps(res))
        print(f'{version:<18s} n={n} rho={rho:.3f}  LR={obs:.3f}  p={p:.4f}  '
              f'groups {m1:+.3f} / {m2:+.3f} (low weight {pi_:.2f}), tau {tau2:.3f} '
              f'(one group: tau {tau1:.3f})  classes agree {row["classes_agree"]}/{n}  '
              f'null median {row["null_median_LRT"]:.2f}', flush=True)
